"""Single-threshold GMYC species delimitation (R package splits) on the rerooted core-genome tree.

Writes gmyc_clusters.csv + gmyc_summary.csv (from steps/gmyc.R) and gmyc_membership.csv in
ecotype_membership.csv's columns, so load_membership() and the niche/ANI summaries can read it.

    python -m steps.gmyc tree.nwk out_dir [--outgroup outgroup] [--clone-map clone_groups.csv]
"""
import csv
import shutil
import subprocess
from pathlib import Path

from steps.parsing import load_clone_map

R_SCRIPT = Path(__file__).with_name("gmyc.R")
INSTALL_HINT = ('Rscript -e \'install.packages("ape"); '
                'install.packages("splits", repos=c("http://R-Forge.R-project.org", "https://cloud.r-project.org"))\'')


def check_r():
    """Raise RuntimeError if Rscript, ape or splits is missing."""
    rscript = shutil.which("Rscript")
    if not rscript:
        raise RuntimeError(f"Rscript not found on PATH. Install R, then: {INSTALL_HINT}")
    probe = subprocess.run([rscript, "-e", 'quit(status = !all(sapply(c("ape", "splits"), requireNamespace, '
                                           'quietly = TRUE)))'], capture_output=True)
    if probe.returncode:
        raise RuntimeError(f"R packages ape and/or splits are missing. Install with: {INSTALL_HINT}")
    return rscript


def run_gmyc(tree_path, out_dir, outgroup=None, clone_map=None):
    """Run GMYC on tree_path, dropping the outgroup tip. clone_map (dedup_clones.py's clone_groups.csv)
    gives collapsed genomes their representative's cluster. Returns the gmyc_membership.csv path."""
    rscript = check_r()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    run = subprocess.run([rscript, str(R_SCRIPT), str(tree_path), str(out_dir)] + ([outgroup] if outgroup else []),
                         capture_output=True, text=True)
    if run.returncode:
        raise RuntimeError(f"gmyc.R failed on {tree_path}:\n{run.stderr.strip()}")

    with open(out_dir / "gmyc_clusters.csv", newline="") as f:
        cluster_of = {r["genome"]: r["cluster"] for r in csv.DictReader(f)}
    members = load_clone_map(clone_map)
    rows = [(cluster, taxon, rep) for rep, cluster in cluster_of.items() for taxon in members.get(rep, [rep])]
    size = {c: sum(r[0] == c for r in rows) for c in cluster_of.values()}
    out = out_dir / "gmyc_membership.csv"
    with open(out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "gene_count", "trial", "ecotype", "ecotype_size", "taxon", "representative"])
        # gene_count left blank: load_membership int()s it when the file is not the EcoSim full-core XML
        writer.writerows(["gmyc", "", "", c, size[c], taxon, rep] for c, taxon, rep in rows)
    with open(out_dir / "gmyc_summary.csv", newline="") as f:
        s = next(csv.DictReader(f))
    print(f"GMYC: {s['entities']} entities (CI {s['entities_ci_low']}-{s['entities_ci_high']}), "
          f"{s['clusters']} clusters, LR test p = {float(s['lr_test_p']):.3g}. Saved {out}")
    return out


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tree")
    parser.add_argument("out_dir")
    parser.add_argument("--outgroup", default="outgroup", help="tip to drop before GMYC ('' to keep every tip)")
    parser.add_argument("--clone-map")
    args = parser.parse_args()
    run_gmyc(args.tree, args.out_dir, args.outgroup or None, args.clone_map)
