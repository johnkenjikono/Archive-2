"""Optional alternative population/ecotype delimitation methods, run beside EcoSim.

Tree/alignment methods (after step 6, on the deduplicated core alignment and rerooted tree):
    TreeCluster (Balaban et al. 2019), mPTP (Kapli et al. 2017), fastbaps (Tonkin-Hill et al. 2019)
Genome methods (in step 2, on the assemblies, which are deleted after Bakta by default):
    PopPUNK (Lees et al. 2019), PopCOGenT (Arevalo et al. 2019)

Each writes <workdir>/delimitation/<method>_membership.csv in ecotype_membership.csv's columns, plus
the tool's raw output in <workdir>/delimitation/<method>/. TreeCluster writes one per threshold
(treecluster_t0.01_membership.csv, ...). find_memberships() lists these files together with GMYC's and
MED's; --niche then tests each method's clusters against host, niche and country, beside EcoSim's
ecotypes (method_niche_association.csv).
"""
import csv
import glob
import os
import shutil
import subprocess
from pathlib import Path

from steps.parsing import load_clone_map

OUTGROUP_STEM = "outgroup"  # download_outgroup.py saves the outgroup genome as outgroup.fna
# max_clade thresholds (substitutions/site) for --treecluster without values: no single cut suits every
# species, so a ladder from near-clonal lineages (0.001) to deep within-species splits (0.03).
TREECLUSTER_LADDER = (0.001, 0.003, 0.01, 0.03)

FASTBAPS_R = """
args <- commandArgs(trailingOnly = TRUE)
suppressMessages(library(fastbaps))
sparse <- optimise_prior(import_fasta_sparse_nt(args[1]), type = "optimise.symmetric")
part <- best_baps_partition(sparse, fast_baps(sparse))
write.csv(data.frame(genome = colnames(sparse$snp.matrix), cluster = part), args[2], row.names = FALSE, quote = FALSE)
"""


def write_membership(out_csv, method, cluster_of, clone_map=None, outgroup=None):
    """Write {genome: cluster} as ecotype_membership.csv rows, expanding collapsed clones and
    dropping the outgroup. Returns out_csv."""
    members = load_clone_map(clone_map)
    rows = [(str(c), taxon, rep) for rep, c in cluster_of.items() if rep != outgroup
            for taxon in members.get(rep, [rep])]
    size = {c: sum(r[0] == c for r in rows) for c, _, _ in rows}
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "gene_count", "trial", "ecotype", "ecotype_size", "taxon", "representative"])
        # gene_count left blank: load_membership int()s it when the file is not the EcoSim full-core XML
        writer.writerows([method, "", "", c, size[c], taxon, rep] for c, taxon, rep in rows)
    print(f"{method}: {len(size)} clusters over {len(rows)} genomes. Saved {out_csv}")
    return out_csv


def _run(cmd, **kwargs):
    if not shutil.which(cmd[0]):
        raise RuntimeError(f"{cmd[0]} not found on PATH (conda install -c bioconda -c conda-forge "
                           "treecluster poppunk r-fastbaps 'setuptools<81'; mptp: build from "
                           "github.com/Pas-Kapli/mptp)")
    run = subprocess.run(cmd, capture_output=True, text=True, **kwargs)
    if run.returncode:
        raise RuntimeError(f"{' '.join(map(str, cmd))} failed:\n{(run.stderr or run.stdout)[-1500:]}")
    return run


def treecluster_thresholds(values):
    """--treecluster's thresholds: the values given, or TREECLUSTER_LADDER when the flag has none."""
    return sorted(set(values)) if values else list(TREECLUSTER_LADDER)


def run_treecluster(tree, out_dir, threshold, outgroup=None, clone_map=None):
    """TreeCluster max_clade: clades whose max pairwise tip distance is <= threshold (subs/site).
    Written as method treecluster_t<threshold>, so several thresholds sit side by side."""
    raw = Path(out_dir) / "treecluster"
    raw.mkdir(parents=True, exist_ok=True)
    method = f"treecluster_t{threshold:g}"
    membership = _fresh(out_dir, method)
    out = raw / f"clusters_t{threshold:g}.tsv"
    _run(["TreeCluster.py", "-i", str(tree), "-o", str(out), "-t", str(threshold), "-m", "max_clade"])
    with open(out) as f:
        rows = list(csv.reader(f, delimiter="\t"))[1:]
    # -1 marks a singleton: give each its own cluster id
    cluster_of = {name: c if c != "-1" else f"singleton_{name}" for name, c in rows}
    return write_membership(membership, method, cluster_of, clone_map, outgroup)


def _family(method):
    """treecluster_t0.01 -> treecluster; med_<species> -> med; others are their own family."""
    return "treecluster" if method.startswith("treecluster_t") else "med" if method.startswith("med") else method


def find_memberships(base_dir, ecosim_out_dir, families=None):
    """{method: membership csv} for every delimitation besides EcoSim that a run wrote: delimitation/,
    GMYC's (in the EcoSim folder) and MED's (med_output_<species>/). families: only these methods
    ("treecluster", "fastbaps", "gmyc", ...), e.g. the ones this run asked for, so files an earlier run
    left in the folder are not tested again."""
    base_dir, ecosim_out_dir = Path(base_dir), Path(ecosim_out_dir)
    found = {p.name[:-len("_membership.csv")]: p for p in sorted((base_dir / "delimitation").glob("*_membership.csv"))}
    if (ecosim_out_dir / "gmyc_membership.csv").exists():
        found["gmyc"] = ecosim_out_dir / "gmyc_membership.csv"
    for med in sorted(base_dir.glob("med_output_*/ecotype_membership.csv")):
        found["med" if "med" not in found else f"med_{med.parent.name[len('med_output_'):]}"] = med
    return {m: p for m, p in found.items() if families is None or _family(m) in families}


def clear_treecluster(out_dir):
    """Remove the TreeCluster memberships of an earlier run, whose thresholds may differ from this one's."""
    for old in Path(out_dir).glob("treecluster_t*_membership.csv"):
        old.unlink()


def _fresh(out_dir, method):
    """The membership path for method, with any earlier run's file removed, so a failure leaves none behind."""
    path = Path(out_dir) / f"{method}_membership.csv"
    path.unlink(missing_ok=True)
    return path


def _ingroup(genomes):
    """The assemblies without the downloaded outgroup: the genome methods model one species, and another
    species' distances would set their thresholds (PopPUNK split 24 C. trachomatis into 10 clusters with
    C. muridarum in its fit, 3 without)."""
    return [g for g in genomes if Path(g).stem != OUTGROUP_STEM]


def parse_mptp(txt):
    """{tip: species number} from mPTP's output file ('Species 1:' then one tip per line)."""
    cluster_of, current = {}, None
    for line in Path(txt).read_text().splitlines():
        line = line.strip()
        if line.startswith("Species ") and line.endswith(":"):
            current = line[len("Species "):-1]
        elif line and current:
            cluster_of[line] = current
    return cluster_of


def run_mptp(tree, out_dir, outgroup=None, clone_map=None):
    """mPTP maximum-likelihood delimitation, rooted on and cropping the outgroup."""
    raw = Path(out_dir) / "mptp"
    raw.mkdir(parents=True, exist_ok=True)
    membership = _fresh(out_dir, "mptp")
    cmd = ["mptp", "--ml", "--multi", "--tree_file", str(tree), "--output_file", str(raw / "mptp")]
    if outgroup:
        cmd += ["--outgroup", outgroup, "--outgroup_crop"]
    _run(cmd)
    return write_membership(membership, "mptp", parse_mptp(raw / "mptp.txt"), clone_map, outgroup)


def _without(fasta, drop, out_path):
    """A copy of an alignment without the records named in drop. Returns how many records were left out."""
    keep, dropped = True, 0
    with open(fasta) as src, open(out_path, "w") as dst:
        for line in src:
            if line.startswith(">"):
                keep = line[1:].split()[0] not in drop if line[1:].split() else True
                dropped += not keep
            if keep:
                dst.write(line)
    return dropped


def run_fastbaps(fasta, out_dir, outgroup=None, clone_map=None):
    """fastbaps (R) best BAPS partition of the core alignment, on the ingroup only. Another species' sequence
    takes over BAPS's partition: with the C. muridarum outgroup in the alignment, 23 C. trachomatis genomes
    split 22 + 1; without it, into 4 populations that refine PopPUNK's 3 lineages."""
    raw = Path(out_dir) / "fastbaps"
    raw.mkdir(parents=True, exist_ok=True)
    membership = _fresh(out_dir, "fastbaps")
    if outgroup:
        ingroup = raw / "ingroup.fasta"
        _without(fasta, {outgroup}, ingroup)
        fasta = ingroup
    out = raw / "fastbaps_clusters.csv"
    _run(["Rscript", "-e", FASTBAPS_R, str(fasta), str(out)])
    with open(out, newline="") as f:
        cluster_of = {r["genome"]: r["cluster"] for r in csv.DictReader(f)}
    return write_membership(membership, "fastbaps", cluster_of, clone_map, outgroup)


def run_poppunk(genomes, out_dir, threads=8):
    """PopPUNK database + DBSCAN model fit on the assemblies; clusters named by genome file stem."""
    raw = Path(out_dir) / "poppunk"
    raw.mkdir(parents=True, exist_ok=True)
    membership = _fresh(out_dir, "poppunk")
    genomes = _ingroup(genomes)
    rfile = raw / "rfiles.txt"
    # PopPUNK rewrites '.' and other characters in names, so it gets g0, g1, ... and they are mapped back.
    stem_of = {f"g{i}": Path(g).stem for i, g in enumerate(genomes)}
    rfile.write_text("".join(f"g{i}\t{Path(g).resolve()}\n" for i, g in enumerate(genomes)))
    db = raw / "db"
    _run(["poppunk", "--create-db", "--r-files", str(rfile), "--output", str(db), "--threads", str(threads),
          "--overwrite"])
    _run(["poppunk", "--fit-model", "dbscan", "--ref-db", str(db), "--output", str(db), "--threads", str(threads),
          "--overwrite"])
    with open(db / "db_clusters.csv", newline="") as f:
        cluster_of = {stem_of[r["Taxon"]]: r["Cluster"] for r in csv.DictReader(f)}
    return write_membership(membership, "poppunk", cluster_of, outgroup=OUTGROUP_STEM)


def run_popcogent(genomes, out_dir, threads=8):
    """PopCOGenT gene-flow units (Main_cluster) from pairwise genome alignments.

    Not on conda: needs a PopCOGenT checkout and its conda env (see its README). Set
    POPCOGENT_DIR (the src/PopCOGenT folder), MUGSY_PATH, MUGSY_ENV (mugsyenv.sh), INFOMAP_PATH,
    and POPCOGENT_ENV if its conda env is not named PopCOGenT.
    """
    membership = _fresh(out_dir, "popcogent")
    env = {k: os.environ.get(k) for k in ("POPCOGENT_DIR", "MUGSY_PATH", "MUGSY_ENV", "INFOMAP_PATH")}
    missing = [k for k, v in env.items() if not v]
    if missing:
        raise RuntimeError(f"set {', '.join(missing)} (see the PopCOGenT README for installing it)")
    raw = (Path(out_dir) / "popcogent").resolve()
    genome_dir, align_dir = raw / "genomes", raw / "alignments"
    for d in (genome_dir, align_dir):
        d.mkdir(parents=True, exist_ok=True)
    # PopCOGenT globs one extension and names strains by file stem, so link every genome in as <stem>.fasta.
    for g in _ingroup(genomes):
        link = genome_dir / f"{Path(g).stem}.fasta"
        link.unlink(missing_ok=True)
        link.symlink_to(Path(g).resolve())
    src = Path(env["POPCOGENT_DIR"])
    conda_py = ["conda", "run", "--no-capture-output", "-n", os.environ.get("POPCOGENT_ENV", "PopCOGenT"), "python"]
    _run(conda_py + [str(src / "get_alignment_and_length_bias.py"), "--genome_dir", f"{genome_dir}/",
                     "--genome_ext", ".fasta", "--alignment_dir", str(align_dir), "--mugsy_path", env["MUGSY_PATH"],
                     "--mugsy_env", env["MUGSY_ENV"], "--base_name", "popcogent", "--final_output_dir", str(raw),
                     "--num_threads", str(threads)], cwd=raw)
    # cluster.py writes infomap_out/ into its working directory.
    _run(conda_py + [str(src / "cluster.py"), "--base_name", "popcogent", "--length_bias_file",
                     str(raw / "popcogent.length_bias.txt"), "--output_directory", str(raw),
                     "--infomap_path", env["INFOMAP_PATH"]], cwd=raw)
    tab = glob.glob(str(raw / "popcogent_*.cluster.tab.txt"))
    if not tab:
        raise RuntimeError(f"PopCOGenT wrote no cluster table in {raw}")
    with open(tab[0], newline="") as f:
        cluster_of = {r["Strain"]: r["Main_cluster"] for r in csv.DictReader(f, delimiter="\t")}
    return write_membership(membership, "popcogent", cluster_of, outgroup=OUTGROUP_STEM)
