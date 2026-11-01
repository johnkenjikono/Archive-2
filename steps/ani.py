"""
Do EcoSim ecotypes hold up under whole-genome ANI?

run_fastani runs FastANI all-vs-all on the downloaded genomes (step 2, before they are
deleted); summarize_ani then compares ANI within each ecotype to ANI between ecotypes.

Writes ecotype_ani.csv into the EcoSim output folder: one row per ecotype, an "all" row
pooling every within- and between-ecotype pair, and an "outgroup" row for reference.
FastANI does not report pairs below ~80% ANI, so those pairs are simply absent.

  python steps/ani.py fastani.tsv ecosim_output_dir
"""
import csv
import statistics
import subprocess
import tempfile
from collections import defaultdict
from itertools import combinations
from pathlib import Path

try:
    from steps.niche_metadata import load_membership
except ImportError:  # run as a script: python steps/ani.py
    from niche_metadata import load_membership


def run_fastani(genome_paths, out_tsv, threads=8):
    """All-vs-all FastANI on genome_paths into out_tsv. Returns out_tsv, or None on failure."""
    with tempfile.NamedTemporaryFile("w", suffix=".txt") as genome_list:
        genome_list.write("\n".join(str(p) for p in genome_paths) + "\n")
        genome_list.flush()
        cmd = ["fastANI", "--ql", genome_list.name, "--rl", genome_list.name, "-t", str(threads), "-o", str(out_tsv)]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except FileNotFoundError:
            print("fastANI not found; install it (conda install -c bioconda fastani) or drop --ani.")
            return None
        except subprocess.CalledProcessError as e:
            print(f"fastANI failed ({e.returncode}): {(e.stderr or '')[-500:]}")
            return None
    return out_tsv


def load_ani(ani_tsv):
    """{frozenset({a, b}): ANI} keyed by genome file stem, both directions averaged."""
    both = defaultdict(list)
    with open(ani_tsv) as f:
        for row in csv.reader(f, delimiter="\t"):
            a, b = Path(row[0]).stem, Path(row[1]).stem
            if a != b:
                both[frozenset((a, b))].append(float(row[2]))
    return {pair: statistics.mean(values) for pair, values in both.items()}


def _stats(values):
    return (len(values), f"{statistics.mean(values):.3f}", f"{min(values):.3f}") if values else (0, "", "")


def summarize_ani(ani_tsv, ecosim_out_dir, outgroup="outgroup"):
    """Within- vs between-ecotype ANI. Returns the CSV path, or None without membership data."""
    membership = Path(ecosim_out_dir) / "ecotype_membership.csv"
    if not membership.exists():
        print(f"No {membership.name}; skipping the ANI summary.")
        return None
    xml_used, ecotype_of = load_membership(membership)
    ani = load_ani(ani_tsv)

    members = defaultdict(list)
    for genome, eco in ecotype_of.items():
        members[eco].append(genome)

    def pairs(group_a, group_b):
        return [ani[p] for a in group_a for b in group_b if a != b and (p := frozenset((a, b))) in ani]

    out_csv = Path(ecosim_out_dir) / "ecotype_ani.csv"
    all_within, all_between = [], []
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ecotype", "genomes", "within_pairs", "within_mean_ani", "within_min_ani",
                         "between_pairs", "between_mean_ani", "nearest_ecotype", "nearest_mean_ani",
                         "mean_ani_to_outgroup", "ecotype_file"])
        for eco in sorted(members, key=lambda e: (len(e), e)):
            genomes = members[eco]
            others = [g for e, gs in members.items() if e != eco for g in gs]
            within = [ani[p] for p in map(frozenset, combinations(genomes, 2)) if p in ani]
            between = pairs(genomes, others)
            all_within += within
            all_between += between
            nearest = max(((e, statistics.mean(v)) for e in members if e != eco if (v := pairs(genomes, members[e]))),
                          key=lambda item: item[1], default=("", None))
            to_outgroup = pairs(genomes, [outgroup])
            writer.writerow([eco, len(genomes), *_stats(within), len(between),
                             f"{statistics.mean(between):.3f}" if between else "",
                             nearest[0], "" if nearest[1] is None else f"{nearest[1]:.3f}",
                             f"{statistics.mean(to_outgroup):.3f}" if to_outgroup else "", xml_used])
        # every between pair was counted from both ecotypes' side
        writer.writerow(["all", len(ecotype_of), *_stats(all_within), len(all_between) // 2,
                         f"{statistics.mean(all_between):.3f}" if all_between else "", "", "",
                         "", xml_used])
        to_outgroup = pairs(list(ecotype_of), [outgroup])
        writer.writerow(["outgroup", 1, 0, "", "", len(to_outgroup),
                         f"{statistics.mean(to_outgroup):.3f}" if to_outgroup else "", "", "", "", xml_used])
    if all_within and all_between:
        print(f"ANI within ecotypes {statistics.mean(all_within):.2f}% vs. between {statistics.mean(all_between):.2f}%")
    return out_csv


def _self_test():
    with tempfile.TemporaryDirectory() as d:
        tsv, folder = Path(d) / "ani.tsv", Path(d)
        rows = {("A1", "A2"): 99.9, ("B1", "B2"): 99.8, ("A1", "B1"): 98.0, ("A1", "B2"): 98.0,
                ("A2", "B1"): 98.0, ("A2", "B2"): 98.0}
        with open(tsv, "w") as f:
            for (a, b), v in rows.items():
                f.write(f"/x/{a}.fna\t/x/{b}.fna\t{v}\t10\t10\n/x/{b}.fna\t/x/{a}.fna\t{v}\t10\t10\n")
            f.write("/x/A1.fna\t/x/outgroup.fna\t85.0\t5\t10\n/x/A3.fna\t/x/A1.fna\t99.0\t10\t10\n")
        with open(folder / "ecotype_membership.csv", "w") as f:
            f.write("file,gene_count,trial,ecotype,ecotype_size,taxon,representative\n")
            for taxon, eco in (("A1", 1), ("A2", 1), ("A3", 1), ("B1", 2), ("B2", 2)):
                f.write(f"full_core_genome_results.xml,full,,{eco},2,{taxon},{taxon}\n")
        rows = list(csv.reader(open(summarize_ani(tsv, folder))))
        # A3 matches only A1 (other pairs below FastANI's cutoff): 2 within pairs, mean (99.9 + 99) / 2
        assert rows[1][:9] == ["1", "3", "2", "99.450", "99.000", "4", "98.000", "2", "98.000"], rows[1]
        assert rows[1][9] == "85.000", rows[1]
        assert rows[3][:7] == ["all", "5", "3", "99.567", "99.000", "4", "98.000"], rows[3]
    print("self-test passed")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ani_tsv", nargs="?")
    parser.add_argument("ecosim_out_dir", nargs="?")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        _self_test()
    elif args.ani_tsv and args.ecosim_out_dir:
        print(f"Saved {summarize_ani(args.ani_tsv, args.ecosim_out_dir)}")
    else:
        parser.error("give ani_tsv and ecosim_out_dir, or --self-test")
