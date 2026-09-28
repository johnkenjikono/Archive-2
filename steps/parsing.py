import csv
import re
import xml.etree.ElementTree as ET
from pathlib import Path

_RAREFACTION_RE = re.compile(r"_g(\d+)_t(\d+)")

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
        match = _RAREFACTION_RE.search(result["file"])
        gene_count, trial = match.groups() if match else ("", "")
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
