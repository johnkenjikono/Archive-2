import argparse
import contextlib
import csv
import glob
import io
import subprocess
import sys
import shutil
import os
import zipfile
import json
import random
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from datetime import datetime

# For outgroup download
from steps.download_outgroup import download_ncbi_fasta, find_outgroup_species, find_typestrain

# Imports for the tree and ecosim steps
from steps.move_largest_numeric import move_largest_numeric_to_top
from steps.run_trees import make_trees_batch, write_error_log
from steps.reroot_tree import reroot_tree_by_first_fasta
from steps.rarefaction import create_rarefaction_fastas
from steps.run_ecosim import run_ecosim_batch, _resolve_ecosim_jar
from steps.parsing import summarize_ecotypes_in_folder

# Falls back to the bundled tools/ecosim.jar / tools/ dir (where bin/ lives) when env vars are unset.
ECOSIM_JAR = _resolve_ecosim_jar()
# tools/bin holds macOS arm64 binaries; setup.sh builds Linux ones into tools/linux/bin.
ECOSIM_DIR = os.environ.get("ECOSIM_DIR") or str(
    Path(__file__).resolve().parent / "tools" / ("linux" if sys.platform.startswith("linux") else "")
)

MAX_GFF_FILES = 201

def _cleanup(path, label=None):
    p = str(path)
    name = label or os.path.basename(p)
    try:
        if os.path.isdir(p):
            shutil.rmtree(p)
            print(f"Removed: {name}")
        elif os.path.isfile(p):
            os.remove(p)
            print(f"Removed: {name}")
    except Exception as e:
        print(f"Warning: could not remove {name}: {e}")

def setup_directories(base_dir: Path):
    input_dir = base_dir / "input"
    bakta_dir = base_dir / "intermediate_bakta"
    roary_dir = base_dir / "output_roary"
    
    for d in [input_dir, bakta_dir, roary_dir]:
        d.mkdir(parents=True, exist_ok=True)
        
    return input_dir, bakta_dir, roary_dir

def get_random_accessions(species_name, sample_size=200, random_seed=42):
    cmd = [
        "datasets", "summary", "genome", "taxon", species_name,
        "--assembly-level", "chromosome,complete",
        "--as-json-lines"
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError:
        print(f"Failed to fetch assembly list for {species_name}")
        return []

    accessions = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError:
            continue

        accession = record.get("accession") or record.get("assembly_accession")
        if accession:
            accessions.append(accession)

    # Sort the deduplicated accessions so the sampling population is deterministic
    # (independent of NCBI's response ordering).
    # GCA_x and GCF_x are the GenBank/RefSeq copies of one assembly; keep one (GCF sorts last, wins).
    # Duplicate identical genomes crash Panaroo (KeyError in collapse_families).
    accessions = sorted({a[4:]: a for a in sorted(set(accessions))}.values())
    if not accessions:
        return []

    if len(accessions) <= sample_size:
        print(
            f"Found {len(accessions)} distinct assemblies (fewer than or equal to {sample_size}); downloading all."
        )
        return accessions

    # Seed so the same species + seed always yields the same assemblies.
    selected = random.Random(random_seed).sample(accessions, sample_size)
    print(f"Selected {sample_size} random distinct assemblies out of {len(accessions)} total (seed={random_seed}).")
    return selected

def download_accessions_zip(accessions, zip_filename):
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as temp_file:
        for accession in accessions:
            temp_file.write(f"{accession}\n")
        accession_file = temp_file.name

    cmd = [
        "datasets", "download", "genome", "accession",
        "--inputfile", accession_file,
        "--include", "genome",
        "--filename", zip_filename
    ]

    try:
        subprocess.run(cmd, check=True)
        return True
    except subprocess.CalledProcessError:
        return False
    finally:
        if os.path.exists(accession_file):
            os.remove(accession_file)

def process_species_genome(species_name, sample_size, output_dir, random_seed=42):
    zip_filename = f"{species_name.replace(' ', '_')}.zip"
    
    print(f"\n============================================================")
    print(f"--- Downloading genomes for {species_name} ---")
    print(f"============================================================\n")

    selected_accessions = get_random_accessions(species_name, sample_size, random_seed)
    if not selected_accessions:
        print(f"No chromosome/complete assemblies found for {species_name}")
        sys.exit(1)

    if not download_accessions_zip(selected_accessions, zip_filename):
        print(f"Failed to download selected assemblies for {species_name}")
        sys.exit(1)

    extract_path = "temp_ext"
    with zipfile.ZipFile(zip_filename, 'r') as zip_ref:
        zip_ref.extractall(extract_path)

    data_path = os.path.join(extract_path, "ncbi_dataset", "data")
    
    if os.path.exists(data_path):
        os.makedirs(output_dir, exist_ok=True)
        for root, dirs, files in os.walk(data_path):
            for file in files:
                if file.endswith(".fna"):
                    old_file_path = os.path.join(root, file)
                    accession = os.path.basename(root)
                    new_file_path = os.path.join(output_dir, f"{accession}.fna")
                    shutil.move(old_file_path, new_file_path)
                    print(f"Moved: {accession}.fna to {output_dir}/")

    if os.path.exists(zip_filename):
        os.remove(zip_filename)
    if os.path.exists(extract_path):
        shutil.rmtree(extract_path)
        
    print(f"Download complete. Files are in {output_dir}/\n")

def find_fasta_files(input_dir: Path):
    extensions = ("*.fasta", "*.fa", "*.fna")
    files = []
    for ext in extensions:
        files.extend(input_dir.glob(ext))
    return files

def find_gff_files(bakta_dir: Path):
    """Recursively find all .gff3 files in the bakta directory."""
    return list(bakta_dir.glob("**/*.gff3"))

def run_bakta(fasta_file: Path, bakta_dir: Path, db_path: str, threads: int):
    sample_name = fasta_file.stem
    output_dir = bakta_dir / sample_name
    expected_gff = output_dir / f"{sample_name}.gff3"

    if expected_gff.exists():
        print(f"Skipping {fasta_file.name}: Bakta output already exists at {expected_gff}")
        return expected_gff
    
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        # --no-capture-output streams Bakta's log live instead of dumping it when the genome finishes.
        "conda", "run", "--no-capture-output", "-n", "bakta_env", "bakta",
        "--db", str(db_path),
        "--output", str(output_dir),
        "--prefix", sample_name,
        "--threads", str(threads),
        # Panaroo only reads CDS features: skip the slow ncRNA (Infernal) searches and plots.
        "--skip-ncrna", "--skip-ncrna-region", "--skip-plot",
        "--force",
        str(fasta_file)
    ]
    
    print(f"============================================================")
    print(f"Running Bakta on {fasta_file.name}...")
    print(f"Command: {' '.join(cmd)}")
    print(f"============================================================")
    
    try:
        subprocess.run(cmd, check=True)
        print(f"Successfully processed {fasta_file.name}\n")
        return expected_gff
    except subprocess.CalledProcessError as e:
        print(f"Error running Bakta on {fasta_file.name}: {e}\n", file=sys.stderr)
        return None

def run_panaroo(gff_files: list[Path], roary_dir: Path, threads: int, error_log=None):
    # ponytail: dir/arg still named roary_* so the downstream core_gene_alignment.aln
    # path keeps working unchanged.
    if not gff_files:
        error_msg = "No GFF3 files found to run Panaroo."
        print(error_msg, file=sys.stderr)
        if error_log:
            write_error_log(error_log, f"SKIP: {error_msg}")
        return

    run_output_dir = roary_dir / "results"

    if error_log is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        error_log = str(roary_dir / f"panaroo_errors_{timestamp}.log")

    # Initialize error log
    try:
        os.makedirs(os.path.dirname(error_log), exist_ok=True)
        open(error_log, "w").close()
    except Exception as e:
        print(f"Warning: Could not initialize error log at {error_log}: {e}")

    if run_output_dir.exists():
        print(f"Removing previous Panaroo results directory at {run_output_dir}...")
        shutil.rmtree(run_output_dir)
    run_output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        # Wrapper works around a panaroo <=1.8.0 crash; see steps/panaroo_patched.py.
        "conda", "run", "-n", "panaroo_env", "python",
        str(Path(__file__).parent / "steps" / "panaroo_patched.py"),
        "-i", *[str(f) for f in gff_files],
        "-o", str(run_output_dir),
        "--clean-mode", "strict",
        "-a", "core",            # build core-gene alignment (core_gene_alignment.aln)
        "-t", str(threads),
    ]

    print(f"============================================================")
    print(f"Running Panaroo on {len(gff_files)} files...")
    print(f"Command: panaroo -i *.gff3 -o {run_output_dir} --clean-mode strict -a core -t {threads}")
    print(f"============================================================")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        success_msg = f"Successfully ran Panaroo! Results available in: {run_output_dir}"
        print(f"\n✅ {success_msg}")
        write_error_log(error_log, success_msg)
    except subprocess.CalledProcessError as e:
        error_msg = f"Panaroo failed with exit code {e.returncode}"
        print(f"\n❌ {error_msg}", file=sys.stderr)
        # Full stderr to the log: the real traceback sits between tqdm bars and conda's command echo.
        write_error_log(error_log, f"{error_msg} | stderr:\n{e.stderr}")
        print(f"Error log: {error_log}")
        sys.exit(1)
    except Exception as e:
        error_msg = f"Exception running Panaroo: {str(e)}"
        print(f"\n❌ {error_msg}", file=sys.stderr)
        write_error_log(error_log, error_msg)
        print(f"Error log: {error_log}")
        sys.exit(1)

def run_tree_pipeline(input_fasta, start_step, outgroup_id, base_dir, output_tree_name=None, threads=8):
    if not os.path.exists(input_fasta):
        print(f" Error: Input fasta '{input_fasta}' does not exist.")
        sys.exit(1)

    species_name = os.path.basename(input_fasta).replace(".fasta", "").replace(".aln", "")
    print("\n==========================================")
    print(f" Starting Tree Pipeline for {species_name} at Step {start_step}")
    print("==========================================\n")

    # Define directories and files
    temp_tree_rdy = os.path.join(base_dir, f"pipeline_temp_{species_name}", "tree_rdy_fastas")
    temp_trees_final = os.path.join(base_dir, f"pipeline_temp_{species_name}", "trees_final")
    final_rerooted_dir = os.path.join(base_dir, "rerooted_trees")
    rarefaction_out_dir = os.path.join(base_dir, f"rarefaction_fastas_{species_name}")
    rarefaction_tree_dir = os.path.join(base_dir, f"rarefaction_trees_{species_name}")
    temp_rare_trees = os.path.join(base_dir, f"pipeline_temp_{species_name}", "rarefaction_trees_unrooted")
    ecosim_out_dir = os.path.join(base_dir, f"ecosim_output_{species_name}")

    os.makedirs(temp_tree_rdy, exist_ok=True)
    os.makedirs(temp_trees_final, exist_ok=True)
    os.makedirs(final_rerooted_dir, exist_ok=True)

    sorted_fasta = os.path.join(temp_tree_rdy, f"{species_name}_sorted.fasta")
    unrooted_tree = os.path.join(temp_trees_final, f"{species_name}_sorted_tree.nwk")
    final_tree_name = f"{output_tree_name}.nwk" if output_tree_name else f"{species_name}.nwk"
    rerooted_tree = os.path.join(final_rerooted_dir, final_tree_name)
    aligned_fasta = os.path.join(temp_tree_rdy, f"{species_name}_aligned.fasta")

    # --- Step 4: Move largest numeric outgroup to top ---
    if start_step <= 4:
        print(f"--- 4. Placing outgroup at the top ---")
        moved_ok = move_largest_numeric_to_top(input_fasta, sorted_fasta, target_id=outgroup_id)
        if not moved_ok:
            print("Tree pipeline requires a Panaroo core-gene alignment with equal-length sequences.")
            print("Please inspect the Panaroo output in output_roary/results/core_gene_alignment.aln.")
            sys.exit(1)
    else:
        print("--- Skipping Step 4: Using existing sorted FASTA ---")

    # --- Step 5: Run FastTree ---
    if start_step <= 5:
        if os.path.exists(unrooted_tree) and os.path.getsize(unrooted_tree) > 0:
            print(f"\n--- Skipping Step 5: Existing tree found for {species_name} ---")
        else:
            print(f"\n--- 5. Building unrooted tree using FastTree ---")
            make_trees_batch(final_folder=temp_tree_rdy, tree_folder=temp_trees_final, threads=threads)

            if not os.path.exists(unrooted_tree) or os.path.getsize(unrooted_tree) == 0:
                print("FastTree failed to produce an output tree. Exiting pipeline.")
                sys.exit(1)
    else:
        print("\n--- Skipping Step 5: Using existing unrooted tree ---")

    # --- Step 6: Reroot Tree ---
    if start_step <= 6:
        print(f"\n--- 6. Rerooting tree using outgroup ---")
        reroot_tree_by_first_fasta(
            sorted_fasta,
            unrooted_tree,
            rerooted_tree,
            outgroup_name=outgroup_id,
        )
        if not os.path.exists(rerooted_tree):
            print("Rerooting tree failed. Exiting pipeline.")
            sys.exit(1)
    else:
        print("\n--- Skipping Step 6: Using existing rerooted tree ---")

    # --- Step 7: Rarefaction Fasta Creation ---
    if start_step <= 7:
        print(f"\n--- 7. Creating Rarefaction FASTAs ---")
        create_rarefaction_fastas(input_fasta=sorted_fasta, output_folder=rarefaction_out_dir)
        print(f"\n--- 7b. Building one tree per rarefaction replicate ---")
        make_trees_batch(
            final_folder=rarefaction_out_dir,
            tree_folder=temp_rare_trees,
            threads=threads,
        )
        os.makedirs(rarefaction_tree_dir, exist_ok=True)
        rerooted = 0
        for fasta in sorted(glob.glob(os.path.join(rarefaction_out_dir, "*.fasta"))):
            stem = Path(fasta).stem
            unrooted = os.path.join(temp_rare_trees, f"{stem}_tree.nwk")
            if not (os.path.exists(unrooted) and os.path.getsize(unrooted) > 0):
                continue
            # Rerooting is chatty and runs 100x; only show the log when it fails.
            log = io.StringIO()
            with contextlib.redirect_stdout(log):
                ok = reroot_tree_by_first_fasta(
                    fasta,
                    unrooted,
                    os.path.join(rarefaction_tree_dir, f"{stem}_tree.nwk"),
                    outgroup_name=outgroup_id,
                )
            if ok:
                rerooted += 1
            else:
                print(f"  {stem}: reroot failed, falling back to the full tree")
                print(log.getvalue())
        print(f"Rerooted {rerooted} sub-alignment trees in {rarefaction_tree_dir}")
    else:
        print("\n--- Skipping Step 7: Using existing rarefaction FASTAs and trees ---")

    # --- Step 8: Run EcoSim ---
    if start_step <= 8:
        print(f"\n--- 8. Running EcoSim ---")
        if not (ECOSIM_JAR and os.path.exists(ECOSIM_JAR)):
            print("EcoSim jar not found. Set ECOSIM_JAR: export ECOSIM_JAR=~/ecosim/ecosim.jar")
            sys.exit(1)  # keep rarefaction FASTAs + tree so the run can resume with --start-step 8
        succeeded = run_ecosim_batch(
            fasta_dir=rarefaction_out_dir,
            tree_dir=rarefaction_tree_dir,
            full_tree_path=rerooted_tree,
            output_dir=ecosim_out_dir,
            ecosim_jar=ECOSIM_JAR,
            ecosim_dir=ECOSIM_DIR,
            memory_gb=6,
            threads=threads,
        )
        if not succeeded:
            print("EcoSim produced no results; keeping inputs for a rerun with --start-step 8.")
            sys.exit(1)
        _cleanup(rarefaction_out_dir, "rarefaction FASTAs")
        _cleanup(rarefaction_tree_dir, "rarefaction trees")
        _cleanup(rerooted_tree, "rerooted tree")
    else:
        print("\n--- Skipping Step 8 ---")

    # --- Step 9: Parse Results ---
    if start_step <= 9:
        print(f"\n--- 9. Parsing EcoSim Output ---")
        summary = summarize_ecotypes_in_folder(ecosim_out_dir)
        print("\n📊 Final Ecotype Counts per File:")
        if summary:
            for filename, count in sorted(summary.items()):
                print(f"{filename}: {count} ecotypes")
            csv_path = os.path.join(ecosim_out_dir, "ecotype_summary.csv")
            with open(csv_path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["file", "ecotype_count"])
                for filename, count in sorted(summary.items()):
                    writer.writerow([filename, count])
            print(f"Summary saved to {csv_path}")
        else:
            print("No ecotypes found or error in parsing.")

    print(f"\nPipeline for {species_name} completed successfully!")
    
    # Cleanup temporary directories
    if start_step <= 9:
        try:
            shutil.rmtree(os.path.join(base_dir, f"pipeline_temp_{species_name}"))
            print("Cleaned up temporary tree directory.")
        except Exception as e:
            print(f"Warning: Could not clean up temporary directory pipeline_temp_{species_name}: {e}")

def load_species_outgroups_from_csv(csv_path):
    """Load species and outgroup pairs from the summary CSV, starting at row 2."""
    rows = []

    with open(csv_path, newline="", encoding="utf-8-sig") as csv_file:
        reader = csv.reader(csv_file)
        next(reader, None)  # skip header row

        for row_number, row in enumerate(reader, start=2):
            if not row or all(not cell.strip() for cell in row):
                continue

            species_name = row[0].strip() if len(row) > 0 else ""
            outgroup_name = row[3].strip() if len(row) > 3 else ""

            if not species_name:
                continue

            rows.append((row_number, species_name, outgroup_name))

    return rows

def run_pipeline_for_species(species_name, outgroup_name, args):
    """Run the full pipeline for one species/outgroup pair."""
    # Create a subfolder for this run using the species name (spaces replaced with underscores)
    run_folder = Path(args.workdir).resolve() / species_name.replace(" ", "_")
    run_folder.mkdir(parents=True, exist_ok=True)
    print(f"Base repository directory: {run_folder}")

    base_dir = run_folder
    input_dir, bakta_dir, roary_dir = setup_directories(base_dir)

    if args.setup_only:
        print("Setup complete. Place your .fasta files into the input directory and run the script again without the --setup-only flag.")
        return

    core_alignment_path = None
    outgroup_fasta_path = None

    # --- Step 1: Download ---
    if args.start_step <= 1:
        if not shutil.which("datasets"):
            print("Error: NCBI 'datasets' CLI tool not found. Please install it (e.g. via conda) and ensure it's in your PATH.", file=sys.stderr)
            sys.exit(1)
        process_species_genome(species_name, args.sample_size, str(input_dir), args.seed)

    # --- Outgroup: explicit arg, or auto-detect ---
    if not outgroup_name and not args.no_auto_outgroup:
        print(f"\n--- Auto-detecting outgroup for {species_name} ---")
        og_acc, og_org = find_outgroup_species(species_name)
        if og_acc:
            outgroup_name = og_acc
            print(f"Auto-detected outgroup: {og_org} ({og_acc})")
        else:
            print("Could not auto-detect outgroup; continuing without one.")

    # --- Outgroup Download (only needed when it will be annotated in Step 2) ---
    if outgroup_name and args.start_step <= 2:
        print(f"\n--- Downloading outgroup genome: {outgroup_name} ---")
        outgroup_fasta_path = download_ncbi_fasta(outgroup_name, str(input_dir))
        if not outgroup_fasta_path:
            print(f"Error: Could not download outgroup genome for {outgroup_name}", file=sys.stderr)
            sys.exit(1)
        print(f"Outgroup genome downloaded to {outgroup_fasta_path}")

    # --- Step 2: Bakta ---
    if args.start_step <= 2:
        if not args.db:
            print("Error: --db (Bakta database path) is required if starting at Step 1 or 2.", file=sys.stderr)
            sys.exit(1)

        fasta_files = find_fasta_files(input_dir)
        if not fasta_files:
            print(f"No FASTA files found in {input_dir}. Please add some '.fasta', '.fa', or '.fna' files or run Step 1.")
            sys.exit(1)

        # Ensure deterministic ordering and cap the number of files Bakta will process
        fasta_files = sorted(fasta_files)
        if len(fasta_files) > MAX_GFF_FILES:
            print(f"Found {len(fasta_files)} FASTA files; limiting to first {MAX_GFF_FILES} for Bakta processing.")
            fasta_files = fasta_files[:MAX_GFF_FILES]

        # Bakta scales poorly past a few threads, so run several genomes at once.
        jobs = args.bakta_jobs or max(1, args.threads // 4)
        threads_per_job = max(1, args.threads // jobs)
        print(f"Found {len(fasta_files)} FASTA files to process ({jobs} parallel Bakta job(s), {threads_per_job} thread(s) each).")
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            results = list(pool.map(lambda f: run_bakta(f, bakta_dir, args.db, threads_per_job), fasta_files))
        failed = sum(r is None for r in results)
        if failed:
            print(f"Warning: Bakta failed on {failed}/{len(fasta_files)} genome(s).", file=sys.stderr)
        _cleanup(input_dir, "input genomes")

    # --- Step 3: Panaroo ---
    if args.start_step <= 3:
        gff_files = find_gff_files(bakta_dir)
        if not gff_files:
            print(f"No GFF3 files found in {bakta_dir}. Cannot run Panaroo.", file=sys.stderr)
            sys.exit(1)

        # Ensure deterministic ordering and cap the number of GFF files Panaroo will process
        gff_files = sorted(gff_files)
        if len(gff_files) > MAX_GFF_FILES:
            print(f"Found {len(gff_files)} GFF3 files; limiting to first {MAX_GFF_FILES} for Panaroo processing.")
            gff_files = gff_files[:MAX_GFF_FILES]

        # Setup error logging for Panaroo
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        error_log = str(roary_dir / f"panaroo_errors_{timestamp}.log")

        run_panaroo(gff_files, roary_dir, args.threads, error_log)
        core_alignment_path = roary_dir / "results" / "core_gene_alignment.aln"
        _cleanup(bakta_dir, "Bakta annotations")
    else:
        # If starting from step 4+, require an input fasta or use a default one if it exists
        if args.input_fasta:
            core_alignment_path = Path(args.input_fasta)
        else:
            default_path = roary_dir / "results" / "core_gene_alignment.aln"
            if default_path.exists():
                core_alignment_path = default_path
            else:
                core_alignment_path = None

    # --- Steps 4-9: Tree Pipeline ---
    if args.start_step <= 9:
        if not core_alignment_path or not core_alignment_path.exists():
            print(f"\nError: Could not find core gene alignment file. Please provide --input-fasta or run previous steps.")
            sys.exit(1)
        # Outgroup ID priority: explicit flag > outgroup genome > typestrain
        if args.outgroup_id:
            resolved_outgroup_id = args.outgroup_id
        elif outgroup_fasta_path or (outgroup_name and args.start_step > 2):
            resolved_outgroup_id = "outgroup"
        else:
            print(f"\n--- Auto-detecting typestrain for {species_name} ---")
            ts_acc, ts_org = find_typestrain(species_name)
            if ts_acc:
                print(f"No external outgroup; using typestrain {ts_org} ({ts_acc}) as tree root.")
            else:
                print("Could not identify typestrain from NCBI.")
            resolved_outgroup_id = ts_acc

        run_tree_pipeline(
            input_fasta=str(core_alignment_path),
            start_step=max(4, args.start_step),
            outgroup_id=resolved_outgroup_id,
            base_dir=str(base_dir),
            output_tree_name=species_name,
            threads=args.threads
        )
        if args.start_step <= 3:
            _cleanup(roary_dir, "Panaroo output")

def run_csv_batch(csv_path, args):
    """Run the pipeline for every species/outgroup row in the CSV file."""
    rows = load_species_outgroups_from_csv(csv_path)
    if not rows:
        print(f"No species rows found in CSV file: {csv_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(rows)} species row(s) from {csv_path}")

    succeeded = 0
    failed = 0

    for row_number, species_name, outgroup_name in rows:
        print("\n" + "=" * 80)
        print(f"Row {row_number}: {species_name}")
        print(f"Outgroup: {outgroup_name or '(none specified)'}")
        print("=" * 80)

        try:
            run_pipeline_for_species(species_name, outgroup_name, args)
            succeeded += 1
        except SystemExit as exc:
            failed += 1
            code = exc.code if isinstance(exc.code, int) else 1
            print(f"\nRow {row_number} failed with exit code {code}; continuing to next row.")
        except Exception as exc:
            failed += 1
            print(f"\nRow {row_number} failed: {exc}; continuing to next row.")

    print("\n" + "=" * 80)
    print(f"CSV batch complete. Succeeded: {succeeded}, Failed: {failed}")
    print("=" * 80)

def main():
    parser = argparse.ArgumentParser(
        description="A complete end-to-end pipeline from downloading genomes to phylogenetic trees and EcoSim."
    )

    parser.add_argument("--csv-file", help="CSV file with Species Name and Outgroup columns to process row by row from row 2 onward.")
    
    # Download Step arguments
    parser.add_argument("--species", help="Species name to download from NCBI (e.g., 'Treponema paraluiscuniculi'). Required if starting at Step 1.")
    parser.add_argument("--sample-size", type=int, default=200, help="Number of random distinct assemblies to download if starting at Step 1 (default: 200)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducible assembly sampling (default: 42)")
    
    # Bakta/Panaroo arguments
    parser.add_argument("-d", "--db", help="Path to the Bakta database. Required if starting at Step 2.")
    parser.add_argument("-w", "--workdir", default=".", help="Base directory for the pipeline (default: current directory)")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Total threads for Bakta, Panaroo, tree building and EcoSim (default: 8)")
    parser.add_argument("--bakta-jobs", type=int, default=None, help="Genomes to annotate in parallel (default: threads // 4). Lower it if RAM is tight.")
    
    # Control flow arguments
    parser.add_argument("--start-step", type=int, default=1, choices=range(1, 10), 
                        help="Step to start at (1=Download, 2=Bakta, 3=Panaroo, 4=Sort, 5=FastTree, 6=Reroot, 7=Rarefaction, 8=EcoSim, 9=Parse)")
    
    # Tree arguments
    parser.add_argument("--input-fasta", help="Optional: Path to core gene alignment FASTA. Required if starting from Step 4 or later.")
    parser.add_argument("--outgroup-id", default=None, help="Optional outgroup ID to force rooting (supports GCA/GCF and .fna/.fasta suffixes)")
    
    parser.add_argument("--setup-only", action="store_true", help="Only create the directory structure and exit")

    # Outgroup genome argument
    parser.add_argument("--outgroup", help="NCBI accession or species name for outgroup genome to download. Auto-detected from genus if omitted.")
    parser.add_argument("--no-auto-outgroup", action="store_true", help="Disable automatic outgroup detection; run without an outgroup unless --outgroup is set.")
    
    args = parser.parse_args()

    if args.csv_file:
        run_csv_batch(args.csv_file, args)
        return
    

    if not args.species:
        print("Error: --species is required unless --csv-file is provided.", file=sys.stderr)
        sys.exit(1)

    run_pipeline_for_species(args.species, args.outgroup, args)

if __name__ == "__main__":
    main()
