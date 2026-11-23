"""
Which accessory genes go with which ecotype?

Reads ecotype_membership.csv from an EcoSim output folder and Panaroo's gene_presence_absence.csv
(the pipeline keeps it by default as <run>/pangenome/gene_presence_absence.csv; --keep pangenome
is the item that holds it). Genome columns are matched to membership taxa by name, so collapsed
clones are covered and the outgroup (not in the membership file) is ignored.

Writes into the EcoSim output folder:
  gene_ecotype_association.csv  per accessory gene: permutation p (same test as niche_metadata) and BH q
  ecotype_unique_genes.csv      genes in nearly all of one ecotype and nearly none of the other genomes

  python -m steps.gene_content ecosim_output_dir pangenome/gene_presence_absence.csv
"""
import argparse
import csv
import random
from collections import Counter
from pathlib import Path

try:
    from steps.niche_metadata import load_membership, association
except ImportError:  # run as a script: python steps/gene_content.py
    from niche_metadata import load_membership, association

# Gene, Non-unique Gene name, Annotation in Panaroo's own CSV. Later columns are kept only if they name a
# membership taxon, so the _roary.csv's extra fixed columns (No. isolates, ...) are skipped too.
PANAROO_FIXED_COLUMNS = 3
NEARLY_ALL = 0.9  # "unique": present in >= this fraction of the ecotype and in <= 1 - this outside it


def load_presence(path, taxa):
    """Return ([(gene, annotation, {taxon present})], taxa found in the table)."""
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        cols = [i for i in range(PANAROO_FIXED_COLUMNS, len(header)) if header[i] in taxa]
        genes = [(r[0], r[2], {header[i] for i in cols if r[i]}) for r in reader]
    return genes, {header[i] for i in cols}


def bh_qvalues(pvalues):
    """Benjamini-Hochberg q-values, in input order."""
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i], reverse=True)
    q, running = [0.0] * len(pvalues), 1.0
    for rank, i in zip(range(len(pvalues), 0, -1), order):
        running = min(running, pvalues[i] * len(pvalues) / rank)
        q[i] = running
    return q


def _fmt(counts, sizes):
    return " ".join(f"{e}:{counts[e]}/{sizes[e]}" for e in sorted(sizes, key=lambda e: (len(e), e)))


def run_gene_content(ecosim_out_dir, presence_absence_path):
    """Run both analyses on one EcoSim output folder. Returns (association CSV, unique-genes CSV), or None."""
    out = Path(ecosim_out_dir)
    membership_csv = out / "ecotype_membership.csv"
    if not membership_csv.exists() or not Path(presence_absence_path).exists():
        print(f"Need {membership_csv} and {presence_absence_path}; skipping gene content analysis.")
        return None
    xml_used, ecotype_of = load_membership(membership_csv)
    genes, found = load_presence(presence_absence_path, set(ecotype_of))
    ecotype_of = {a: e for a, e in ecotype_of.items() if a in found}
    sizes = Counter(ecotype_of.values())
    print(f"Gene content on {xml_used}: {len(ecotype_of)} genomes, {len(sizes)} ecotypes, {len(genes)} genes")

    rows, p_cache = [], {}
    for gene, annotation, present in genes:
        present = present & ecotype_of.keys()
        if len(present) <= 1 or len(present) == len(ecotype_of):
            continue  # absent, in one genome only, or core: nothing to associate
        key = frozenset(present)
        if key not in p_cache:  # accessory genes often share a pattern; same data and seed give the same p
            p_cache[key] = association(ecotype_of, {a: "yes" if a in present else "no" for a in ecotype_of})[5]
        counts = Counter(ecotype_of[a] for a in present)
        rows.append([gene, annotation, len(present), _fmt(counts, sizes), p_cache[key]])
    q = bh_qvalues([r[4] for r in rows])
    rows = sorted((r + [round(x, 4)] for r, x in zip(rows, q)), key=lambda r: (r[5], r[4], r[0]))
    assoc_path = out / "gene_ecotype_association.csv"
    with open(assoc_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "annotation", "genomes_with_gene", "present_by_ecotype", "permutation_p", "q_bh"])
        writer.writerows(rows)
    print(f"Saved {assoc_path} ({len(rows)} accessory genes, {len(p_cache)} distinct patterns)")

    # Ecotypes with a single genome are skipped: "all of it" means nothing for one genome.
    unique_path = out / "ecotype_unique_genes.csv"
    n_unique = 0
    with open(unique_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ecotype", "gene", "annotation", "fraction_inside", "fraction_outside"])
        for gene, annotation, present in genes:
            counts = Counter(ecotype_of[a] for a in present if a in ecotype_of)
            for eco, size in sizes.items():
                if size < 2:
                    continue
                inside = counts[eco] / size
                outside = (sum(counts.values()) - counts[eco]) / (len(ecotype_of) - size) if len(ecotype_of) > size else 0
                if inside >= NEARLY_ALL and outside <= 1 - NEARLY_ALL:
                    writer.writerow([eco, gene, annotation, round(inside, 3), round(outside, 3)])
                    n_unique += 1
    print(f"Saved {unique_path} ({n_unique} ecotype-specific genes)")
    return assoc_path, unique_path


def _self_test():
    import tempfile
    rng = random.Random(1)
    taxa = [f"GCA_{i}" for i in range(30)]
    eco = {t: "1" if i < 15 else "2" for i, t in enumerate(taxa)}
    random_genes = {a: rng.random() < 0.5 for a in taxa}
    tied = lambda t: eco[t] == "1"  # noqa: E731
    table = {  # gene -> (annotation, present?)
        "tied": ("tied to 1", tied), "rand": ("random", lambda t: random_genes[t]),
        "core": ("core", lambda t: True), "single": ("one genome", lambda t: t == "GCA_3"), "none": ("none", lambda t: False),
    }
    with tempfile.TemporaryDirectory() as d:
        pa = Path(d) / "gene_presence_absence.csv"
        with open(pa, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["Gene", "Non-unique Gene name", "Annotation"] + taxa + ["outgroup"])  # Panaroo's layout
            for g, (ann, has) in table.items():
                writer.writerow([g, "", ann] + [g if has(t) else "" for t in taxa] + [g])
        with open(Path(d) / "ecotype_membership.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["file", "gene_count", "trial", "ecotype", "ecotype_size", "taxon", "representative"])
            for t in taxa:
                writer.writerow(["full_core_genome_results.xml", "full", "", eco[t], 15, t, t])
        assoc_path, unique_path = run_gene_content(d, pa)
        assoc = list(csv.DictReader(open(assoc_path)))
        assert [r["gene"] for r in assoc] == ["tied", "rand"], assoc  # core, single and none skipped; sorted by q
        assert assoc[0]["present_by_ecotype"] == "1:15/15 2:0/15" and float(assoc[0]["q_bh"]) < 0.01, assoc[0]
        assert float(assoc[1]["q_bh"]) > 0.05, assoc[1]
        assert [(r["ecotype"], r["gene"]) for r in csv.DictReader(open(unique_path))] == [("1", "tied")]
    assert [round(x, 3) for x in bh_qvalues([0.01, 0.04, 0.03, 0.2])] == [0.04, 0.053, 0.053, 0.2]
    print("self-test passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ecosim_out_dir", nargs="?")
    parser.add_argument("presence_absence", nargs="?")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        _self_test()
    elif args.ecosim_out_dir and args.presence_absence:
        run_gene_content(args.ecosim_out_dir, args.presence_absence)
    else:
        parser.error("give ecosim_out_dir and presence_absence, or --self-test")
