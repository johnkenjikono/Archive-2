import csv
import re
import statistics
import xml.etree.ElementTree as ET
from collections import Counter
from math import comb
from pathlib import Path

try:
    from steps.rarefaction import FULL_CORE_NAME
except ImportError:  # run as a script: python steps/parsing.py
    from rarefaction import FULL_CORE_NAME

_RAREFACTION_RE = re.compile(r"_g(\d+)_t(\d+)")
FULL_CORE_XML = f"{FULL_CORE_NAME}_results.xml"

def replicate_key(filename):
    """(gene_count, trial) from a rarefaction file name; ("full", "") for the full core genome."""
    if Path(filename).name.startswith(FULL_CORE_NAME):
        return "full", ""
    match = _RAREFACTION_RE.search(filename)
    return match.groups() if match else ("", "")

def count_ecotypes_in_file(file_path):
    try:
        # Count only demarcated ecotypes (matches post_processing.py). A bare
        # './/ecotype' also matches ecotype elements outside the demarcation
        # block and inflates the count.
        return len(ET.parse(file_path).getroot().findall('.//demarcation/ecotypes/ecotype'))
    except ET.ParseError as e:
        print(f"Error parsing {file_path}: {e}")
        return None

def summarize_ecotypes_in_folder(folder_path):
    summary = {}
    for xml_file in Path(folder_path).glob('*.xml'):
        count = count_ecotypes_in_file(xml_file)
        if count is not None:
            summary[xml_file.name] = count
    return summary

def load_clone_map(clone_map_path):
    """Read the clone_groups.csv from dedup_clones.py into {representative: [members]}."""
    groups = {}
    if clone_map_path and Path(clone_map_path).exists():
        with open(clone_map_path, newline="") as handle:
            for row in csv.DictReader(handle):
                groups.setdefault(row["representative"], []).append(row["member"])
    return groups

def parse_ecosim_xml(file_path):
    """
    Parse one EcoSim result XML into a dict with the outgroup, the newick tree
    EcoSim ran on, the hill-climb npop/omega/sigma, and the demarcated ecotypes
    as [{"number", "members": [...]}, ...]. Returns None if the XML is invalid.
    """
    try:
        root = ET.parse(file_path).getroot()
    except ET.ParseError as e:
        print(f"Error parsing {file_path}: {e}")
        return None

    outgroup = root.find("./phylogeny/outgroup")
    tree = root.find("./phylogeny/tree")
    hillclimb = root.find("./hillclimb/result")
    ecotypes = [
        {"number": int(eco.get("number", i + 1)),
         "members": [m.get("name") for m in eco.findall("member")]}
        for i, eco in enumerate(root.findall("./demarcation/ecotypes/ecotype"))
    ]
    return {
        "file": Path(file_path).name,
        "outgroup": outgroup.get("value") if outgroup is not None else None,
        "tree": tree.get("value") if tree is not None else None,
        "npop": hillclimb.get("npop") if hillclimb is not None else None,
        "omega": hillclimb.get("omega") if hillclimb is not None else None,
        "sigma": hillclimb.get("sigma") if hillclimb is not None else None,
        "ecotypes": ecotypes,
    }

def get_ecotype_membership(file_path, clone_map=None):
    """
    {taxon: ecotype_number} for one EcoSim XML. With a clone_map, genomes that
    were collapsed before EcoSim get their representative's ecotype.
    """
    result = parse_ecosim_xml(file_path)
    if result is None:
        return None
    clone_map = clone_map or {}
    return {taxon: eco["number"]
            for eco in result["ecotypes"]
            for member in eco["members"]
            for taxon in clone_map.get(member, [member])}

def write_membership_csv(folder_path, clone_map_path=None, csv_path=None):
    """
    Write taxon -> ecotype for every EcoSim XML in folder_path to
    ecotype_membership.csv (one row per replicate and genome, clones expanded).
    Returns the CSV path, or None if there was nothing to parse.
    """
    clone_map = load_clone_map(clone_map_path)
    rep_of = {m: rep for rep, members in clone_map.items() for m in members}

    rows = []
    for xml_file in sorted(Path(folder_path).glob("*.xml")):
        result = parse_ecosim_xml(xml_file)
        if result is None:
            continue
        gene_count, trial = replicate_key(result["file"])
        for eco in result["ecotypes"]:
            genomes = [t for m in eco["members"] for t in clone_map.get(m, [m])]
            for taxon in genomes:
                rows.append([result["file"], gene_count, trial, eco["number"],
                             len(genomes), taxon, rep_of.get(taxon, taxon)])

    if not rows:
        return None
    csv_path = csv_path or str(Path(folder_path) / "ecotype_membership.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "gene_count", "trial", "ecotype", "ecotype_size",
                         "taxon", "representative"])
        writer.writerows(rows)
    return csv_path

def _adjusted_rand_index(labels_a, labels_b):
    """Adjusted Rand index of two {taxon: cluster} maps, over the taxa they share."""
    shared = labels_a.keys() & labels_b.keys()
    n = len(shared)
    if n < 2:
        return None
    pairs = sum(comb(c, 2) for c in Counter((labels_a[t], labels_b[t]) for t in shared).values())
    rows = sum(comb(c, 2) for c in Counter(labels_a[t] for t in shared).values())
    cols = sum(comb(c, 2) for c in Counter(labels_b[t] for t in shared).values())
    expected = rows * cols / comb(n, 2)
    maximum = (rows + cols) / 2
    if maximum == expected:  # both partitions all-singletons or one cluster
        return 1.0
    return (pairs - expected) / (maximum - expected)

def compare_to_full_core(folder_path):
    """
    Compare every rarefaction replicate with the full core genome run
    (full_core_genome_results.xml). Writes, into folder_path:

    - rarefaction_comparison.csv: per replicate, its ecotype count, the full-core
      count, the adjusted Rand index between the two demarcations, and the share of
      full-core ecotypes it recovers with exactly the same members.
    - rarefaction_by_gene_count.csv: the same, summarised per gene count.
    - ecotype_support.csv: for each full-core ecotype, how often each gene count
      recovers it exactly (a confidence level for that ecotype).

    Members are compared as EcoSim saw them (clone representatives), since every
    replicate uses the same deduplicated sequences. Returns the written paths,
    or [] when there is no full-core result.
    """
    folder = Path(folder_path)
    full = parse_ecosim_xml(folder / FULL_CORE_XML) if (folder / FULL_CORE_XML).exists() else None
    if not full or not full["ecotypes"]:
        return []

    full_labels = {m: eco["number"] for eco in full["ecotypes"] for m in eco["members"]}
    full_sets = {eco["number"]: frozenset(eco["members"]) for eco in full["ecotypes"]}

    replicates = []  # (gene_count, trial, file, n_ecotypes, ari, recovered_set)
    for xml_file in sorted(folder.glob("*.xml")):
        gene_count, trial = replicate_key(xml_file.name)
        if gene_count in ("", "full"):
            continue
        result = parse_ecosim_xml(xml_file)
        if result is None:
            continue
        labels = {m: eco["number"] for eco in result["ecotypes"] for m in eco["members"]}
        sets = {frozenset(eco["members"]) for eco in result["ecotypes"]}
        recovered = {num for num, members in full_sets.items() if members in sets}
        replicates.append((int(gene_count), int(trial), xml_file.name, len(result["ecotypes"]),
                           _adjusted_rand_index(full_labels, labels), recovered))
    if not replicates:
        return []
    replicates.sort()
    n_full = len(full["ecotypes"])

    def fmt(x):
        return "" if x is None else f"{x:.4f}"

    per_replicate = folder / "rarefaction_comparison.csv"
    with open(per_replicate, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "gene_count", "trial", "ecotype_count", "full_core_ecotype_count",
                         "difference", "adjusted_rand_index", "full_core_ecotypes_recovered"])
        for gene_count, trial, name, count, ari, recovered in replicates:
            writer.writerow([name, gene_count, trial, count, n_full, count - n_full,
                             fmt(ari), fmt(len(recovered) / n_full)])

    gene_counts = sorted({r[0] for r in replicates})
    by_count = {g: [r for r in replicates if r[0] == g] for g in gene_counts}
    per_count = folder / "rarefaction_by_gene_count.csv"
    with open(per_count, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene_count", "replicates", "mean_ecotypes", "median_ecotypes", "min_ecotypes",
                         "max_ecotypes", "full_core_ecotype_count", "mean_adjusted_rand_index",
                         "mean_full_core_ecotypes_recovered"])
        for g, reps in by_count.items():
            counts = [r[3] for r in reps]
            aris = [r[4] for r in reps if r[4] is not None]
            writer.writerow([g, len(reps), fmt(statistics.mean(counts)), statistics.median(counts),
                             min(counts), max(counts), n_full,
                             fmt(statistics.mean(aris)) if aris else "",
                             fmt(statistics.mean(len(r[5]) / n_full for r in reps))])

    support = folder / "ecotype_support.csv"
    with open(support, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ecotype", "size"] + [f"support_g{g}" for g in gene_counts] + ["members"])
        for eco in full["ecotypes"]:
            num = eco["number"]
            writer.writerow([num, len(eco["members"])]
                            + [fmt(sum(num in r[5] for r in by_count[g]) / len(by_count[g])) for g in gene_counts]
                            + [";".join(eco["members"])])

    return [per_replicate, per_count, support]

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Ecotype counts and taxon -> ecotype membership from EcoSim XMLs")
    parser.add_argument("folder", nargs="?", default="ecosim_results", help="Folder of EcoSim XML files")
    parser.add_argument("--clone-map", default=None,
                        help="clone_groups.csv from dedup_clones.py (default: <folder>/clone_groups.csv if present)")
    args = parser.parse_args()

    summary = summarize_ecotypes_in_folder(args.folder)
    print("Ecotype Counts per File:")
    for filename, count in sorted(summary.items()):
        print(f"{filename}: {count} ecotypes")

    clone_map = args.clone_map or str(Path(args.folder) / "clone_groups.csv")
    path = write_membership_csv(args.folder, clone_map)
    if path:
        print(f"Membership saved to {path}")
    for path in compare_to_full_core(args.folder):
        print(f"Saved {path}")
