import os
import subprocess
import shutil
import tempfile
from Bio.SeqIO.FastaIO import SimpleFastaParser
from datetime import datetime

def validate_sequences(fasta_path):
    """
    Validate that all sequences are aligned (same length).
    Returns: (num_seqs, seq_length) or (None, None) if invalid
    """
    try:
        # Stream lengths only; don't hold the whole alignment in memory.
        with open(fasta_path) as handle:
            lengths = [len(seq) for _, seq in SimpleFastaParser(handle)]
        if not lengths:
            print(f"  ⚠️  No sequences found in {os.path.basename(fasta_path)}")
            return None, None

        # Check that all sequences have same length
        distinct = set(lengths)
        if len(distinct) > 1:
            print(f"  ⚠️  ERROR: Sequences are not aligned!")
            print(f"     Found {len(distinct)} different lengths: {sorted(distinct)}")
            return None, None

        return len(lengths), lengths[0]
    except Exception as e:
        print(f"  ⚠️  ERROR reading FASTA: {e}")
        return None, None

def write_error_log(error_file, error_message):
    """Write error message to a log file with timestamp."""
    with open(error_file, "a") as f:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"[{timestamp}] {error_message}\n")

def find_fasttree_executable():
    """Find VeryFastTree/FastTree on PATH or common install locations."""
    exe = (shutil.which("veryfasttree") or shutil.which("fasttree")
           or shutil.which("FastTree") or shutil.which("FastTreeMP"))
    if exe:
        return exe
    for path in ["/usr/local/bin/veryfasttree", "/usr/local/bin/fasttree",
                 "/opt/local/bin/veryfasttree"]:
        if os.path.exists(path):
            return path
    raise FileNotFoundError(
        "VeryFastTree not found. Install it and ensure it's on PATH.\n"
        "macOS: brew install veryfasttree\n"
        "Fallback: brew install fasttree"
    )

def make_trees_batch(final_folder="tree_rdy_fastas", tree_folder="trees_final", 
                     stable_mode=True, fasttree_exe=None, max_retries=1,
                     max_sequences_per_run=300, error_log=None, threads=None):
    """
    Build phylogenetic trees using FastTree.
    
    Args:
        final_folder: Input FASTA directory
        tree_folder: Output tree directory
        stable_mode: If True, use conservative flags for stability. If False, use -fastest.
        fasttree_exe: Path to FastTree executable (auto-detect if None)
        max_retries: Number of attempts per file (FastTree is deterministic, so retries rarely help)
        max_sequences_per_run: Max sequences per FASTA file (uses first X sequences, None = no limit)
        error_log: Path to error log file (auto-generated if None)
        threads: Threads for tree building (default: all CPUs)
    """
    os.makedirs(tree_folder, exist_ok=True)

    if error_log is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        error_log = os.path.join(tree_folder, f"tree_building_errors_{timestamp}.log")

    try:
        # Clear error log at start of run
        open(error_log, "w").close()
    except Exception as e:
        print(f"Warning: Could not initialize error log at {error_log}: {e}")

    try:
        if fasttree_exe is None:
            fasttree_exe = find_fasttree_executable()
        print(f"Using tree builder: {fasttree_exe}")
    except Exception as e:
        error_msg = f"FATAL: Could not find FastTree executable: {e}"
        print(f"  ❌ {error_msg}")
        write_error_log(error_log, error_msg)
        raise

    # Support values are never used downstream (EcoSim only needs topology + branch lengths),
    # so skip computing them.
    if stable_mode:
        # Conservative flags: better for stability with diverse sequences
        base_flags = ["-nt", "-gtr", "-gamma", "-nosupport"]
        print("Running in STABLE mode (recommended for large/diverse datasets)")
    else:
        # Faster but less stable flags
        base_flags = ["-nt", "-speediest", "-gtr", "-nosupport"]
        print("Running in FAST mode (less stable)")

    threads = threads or os.cpu_count() or 1
    # VeryFastTree defaults to 1 thread; OpenMP FastTree builds read OMP_NUM_THREADS.
    if "veryfasttree" in os.path.basename(fasttree_exe).lower():
        base_flags += ["-threads", str(threads)]
    run_env = {**os.environ, "OMP_NUM_THREADS": str(threads)}
    print(f"Using {threads} thread(s)")

    if max_sequences_per_run is not None:
        print(f"Sequence cap per FASTA: {max_sequences_per_run}")

    failed_files = []
    succeeded_files = []
    skipped_files = []
    truncated_files = []
    processed_sequences = 0

    for file in sorted(os.listdir(final_folder)):
        if not file.endswith(".fasta"):
            continue

        input_path = os.path.join(final_folder, file)
        output_path = os.path.join(tree_folder, file.replace(".fasta", "_tree.nwk"))

        # Skip only if output exists and is non-empty
        if os.path.exists(output_path):
            if os.path.getsize(output_path) > 0:
                print(f"✓ Skipping {file} (already built)")
                skipped_files.append(file)
                continue
            print(f"⚠️  Found empty tree file for {file}; rebuilding...")
            os.remove(output_path)

        # Validate input
        num_seqs, seq_len = validate_sequences(input_path)
        if num_seqs is None:
            error_msg = f"SKIP: Validation failed for {file}"
            write_error_log(error_log, error_msg)
            skipped_files.append(file)
            continue

        run_input_path = input_path
        run_num_seqs = num_seqs
        temp_subset_path = None

        if max_sequences_per_run is not None and num_seqs > max_sequences_per_run:
            # Build on the first X sequences when FASTA is larger than requested cap.
            run_num_seqs = max_sequences_per_run
            fd, temp_subset_path = tempfile.mkstemp(suffix=".fasta", prefix="tree_subset_")
            os.close(fd)

            with open(input_path) as src, open(temp_subset_path, "w") as dst:
                for i, (title, seq) in enumerate(SimpleFastaParser(src)):
                    if i >= max_sequences_per_run:
                        break
                    dst.write(f">{title}\n{seq}\n")
            run_input_path = temp_subset_path
            truncated_files.append((file, num_seqs, run_num_seqs))
            print(
                f"✂️  {file}: using first {run_num_seqs} of {num_seqs} sequences "
                f"(cap={max_sequences_per_run})"
            )

        print(f"Building tree for {file} ({run_num_seqs} sequences, length={seq_len})...")

        cmd = [fasttree_exe] + base_flags + [run_input_path]

        # Retry mechanism for crashes
        for attempt in range(1, max_retries + 1):
            try:
                with open(output_path, "w") as f:
                    result = subprocess.run(
                        cmd,
                        stdout=f,
                        stderr=subprocess.PIPE,
                        text=True,
                        env=run_env,
                    )

                if result.returncode == 0:
                    print(f"  ✅ Tree written to {output_path}")
                    succeeded_files.append(file)
                    processed_sequences += run_num_seqs
                    break
                else:
                    error_msg = (
                        f"FastTree failed for {file} "
                        f"(Attempt {attempt}/{max_retries}, Exit code: {result.returncode})"
                    )
                    if result.stderr:
                        error_msg += f" | stderr: {result.stderr[:300]}"
                    print(f"  ❌ {error_msg}")
                    write_error_log(error_log, error_msg)
                    if os.path.exists(output_path):
                        os.remove(output_path)

                    if attempt == max_retries:
                        failed_files.append(f"{file} (all retries exhausted)")

            except Exception as e:
                error_msg = f"Exception for {file} (Attempt {attempt}/{max_retries}): {str(e)}"
                print(f"  ⚠️  {error_msg}")
                write_error_log(error_log, error_msg)
                if os.path.exists(output_path):
                    os.remove(output_path)
                if attempt == max_retries:
                    failed_files.append(f"{file} ({str(e)})")

        if temp_subset_path and os.path.exists(temp_subset_path):
            os.remove(temp_subset_path)

    # Summary
    print("\n" + "="*50)
    print("Tree building complete!")
    print(f"  ✅ Succeeded: {len(succeeded_files)}")
    print(f"  🧬 Sequences processed: {processed_sequences}")
    print(f"  ⏭️  Skipped: {len(skipped_files)} (already built)")
    if truncated_files:
        print(f"  ✂️  Truncated by sequence cap: {len(truncated_files)}")
        for name, original_n, used_n in truncated_files:
            print(f"     - {name}: used first {used_n}/{original_n}")
    if failed_files:
        print(f"  ❌ Failed: {len(failed_files)}")
        for f in failed_files:
            print(f"     - {f}")
    print("="*50)

    # Write summary to error log
    if failed_files or skipped_files:
        summary = "\n=== FINAL SUMMARY ===\n"
        summary += f"Succeeded: {len(succeeded_files)}\n"
        summary += f"Failed: {len(failed_files)}\n"
        if failed_files:
            summary += "Failed files:\n"
            for f in failed_files:
                summary += f"  - {f}\n"
        if skipped_files:
            summary += f"Skipped: {len(skipped_files)}\n"
        write_error_log(error_log, summary)
        print(f"Error log written to: {error_log}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build phylogenetic trees with VeryFastTree/FastTree")
    parser.add_argument("--fast", action="store_true", help="Use faster but less stable mode")
    parser.add_argument(
        "--max-sequences",
        type=int,
        default=None,
        help="Maximum sequences to use per FASTA file (takes first X sequences)",
    )
    args = parser.parse_args()

    if args.max_sequences is not None and args.max_sequences <= 0:
        raise ValueError("--max-sequences must be greater than 0")

    make_trees_batch(stable_mode=not args.fast, max_sequences_per_run=args.max_sequences)
