"""
Minimum Entropy Decomposition (MED; Eren et al. 2015, ISME J 9:968-979) on the core-genome alignment,
as a delimitation method to compare with EcoSim.

Input: the deduplicated core alignment the tree is built from (step 4's
pipeline_temp_<species>/tree_rdy_fastas/<species>_sorted.fasta), one row per genome, outgroup included.
Output, in out_dir:
  ecotype_membership.csv  same columns as EcoSim's (file = "med"), so ani.py / niche_metadata.py read it
  med_summary.txt         node count and the parameters used

  python -m steps.med alignment.fasta out_dir [--outgroup ID] [--clone-map clone_groups.csv]
  python -m steps.med --self-test

Canonical MED (oligotyping's `decompose`, kept here):
  - Shannon entropy (log2) per alignment column over A, C, G, T and '-'; other characters (N, IUPAC)
    count toward the column size but are no class. Entropy < 1e-5 is 0.
  - A node is split on its top d columns (default 3) whose entropy exceeds m (default 0.3); every
    distinct combination of bases at those columns becomes a child node. Ties: lower column first.
  - Splitting repeats until no column in a node exceeds m. Nodes smaller than M are removed as outliers;
    outliers are relocated to the closest final node (decompose's optional -R).

Adapted for genomes (decompose assumes millions of reads, many identical):
  - Each genome is one sequence of abundance 1 (clones were collapsed upstream, as for EcoSim). Clone
    counts are not used as abundances, so outbreak oversampling does not drive the entropy.
  - M is a genome count: a node with fewer than M genomes is an outlier. decompose tests the abundance
    of the most abundant unique read instead, which is always 1 here. Default M keeps the canonical
    max(4, round(N / 5000)), i.e. 4 for any realistic genome set.
  - A node is finalized unless its split gives at least two children of M or more genomes (the analogue
    of decompose's "second most abundant unique read < M" rule, which would stop every split here).
  - m is static: decompose's normalized m is scaled by the most abundant read and would hit its 0.05
    floor for every node.
  - Dropped, as they measure read abundance or read error: the competing-unique-sequence-ratio and
    density stops, maximum-variation outlier removal (BLAST), and homopolymer-split merging.
  - Outliers are relocated by Hamming distance on variable columns to each final node's consensus
    (decompose uses BLAST to the node's most abundant read, which does not exist here).
  - Only variable columns are kept; invariant columns have zero entropy and never change the result.
"""
import argparse
import csv
import tempfile
from pathlib import Path

import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser

try:
    from steps.parsing import load_clone_map
except ImportError:  # run as a script: python steps/med.py
    from parsing import load_clone_map

CLASSES = np.frombuffer(b"ACGT-", dtype=np.uint8)


def _is_outgroup(name, outgroup):
    # matched the way step 4 (move_largest_numeric) finds it
    return outgroup is not None and (name == outgroup or outgroup in name)


def read_variable_columns(alignment_fasta, outgroup=None):
    """(genome names, alignment length, uint8 matrix of the variable columns), outgroup excluded."""
    def rows():
        with open(alignment_fasta) as handle:
            for title, seq in SimpleFastaParser(handle):
                name = title.split(None, 1)[0]
                if not _is_outgroup(name, outgroup):
                    yield name, np.frombuffer(seq.upper().encode(), dtype=np.uint8)

    # two passes so only the variable columns are ever held for all genomes
    names, first, variable = [], None, None
    for name, row in rows():
        names.append(name)
        if first is None:
            first, variable = row, np.zeros(len(row), dtype=bool)
        elif len(row) != len(first):
            raise ValueError(f"{name} is {len(row)} columns long, expected {len(first)}: not an alignment")
        else:
            variable |= row != first
    if not names:
        raise ValueError(f"no ingroup sequences in {alignment_fasta}")
    return names, len(first), np.array([row[variable] for _, row in rows()]).reshape(len(names), -1)


def column_entropy(block):
    n = len(block)
    h = np.zeros(block.shape[1])
    for c in CLASSES:
        p = (block == c).sum(axis=0) / n
        h -= np.where(p > 0, p * np.log2(np.where(p > 0, p, 1)), 0)
    h[h < 1e-5] = 0
    return h


def decompose(matrix, m=0.3, d=3, M=4):
    """Return (final nodes, outliers) as arrays of row indices; outliers not yet relocated."""
    if len(matrix) < M:
        raise ValueError(f"{len(matrix)} genomes is fewer than M={M}")
    final, outliers, todo = [], [], [np.arange(len(matrix))]
    while todo:
        idx = todo.pop()
        if len(idx) < M:
            outliers.extend(idx)
            continue
        h = column_entropy(matrix[idx])
        disc = [j for j in np.argsort(-h, kind="stable")[:d] if h[j] > m]
        if not disc:
            final.append(idx)
            continue
        _, label = np.unique(matrix[np.ix_(idx, disc)], axis=0, return_inverse=True)
        label = label.ravel()
        children = [idx[label == k] for k in range(label.max() + 1)]
        if sum(len(c) >= M for c in children) < 2:
            final.append(idx)
        else:
            todo += children
    return final, outliers


def relocate(matrix, final, outliers):
    """Add each outlier to the final node whose consensus is nearest (Hamming on variable columns)."""
    alphabet = np.unique(matrix)
    consensus = [alphabet[np.argmax([(matrix[idx] == a).sum(axis=0) for a in alphabet], axis=0)] for idx in final]
    nearest = [int(np.argmin([(matrix[o] != rep).sum() for rep in consensus])) for o in outliers]
    return [np.concatenate([idx, [o for o, k in zip(outliers, nearest) if k == i]]).astype(int)
            for i, idx in enumerate(final)]


def run_med(alignment_fasta, out_dir, outgroup=None, clone_map=None, m=0.3, d=3, M=None):
    """Write ecotype_membership.csv and med_summary.txt into out_dir; return the number of MED nodes."""
    out_dir = Path(out_dir)
    if any(out_dir.glob("*_results.xml")):
        raise ValueError(f"{out_dir} holds EcoSim results; give MED its own folder so its membership "
                         "does not overwrite EcoSim's")
    names, length, matrix = read_variable_columns(alignment_fasta, outgroup)
    M = max(4, round(len(names) / 5000)) if M is None else M
    final, outliers = decompose(matrix, m, d, M)
    nodes = sorted(relocate(matrix, final, outliers), key=lambda idx: (-len(idx), idx.min()))

    clones = load_clone_map(clone_map)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "ecotype_membership.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "gene_count", "trial", "ecotype", "ecotype_size", "taxon", "representative"])
        for number, idx in enumerate(nodes, 1):
            genomes = [(t, names[i]) for i in sorted(idx) for t in clones.get(names[i], [names[i]])]
            writer.writerows(["med", "", "", number, len(genomes), t, rep] for t, rep in genomes)

    summary = {"alignment": alignment_fasta, "outgroup_excluded": outgroup or "", "genomes": len(names),
               "alignment_length": length, "variable_columns": matrix.shape[1], "m": m, "d": d, "M": M,
               "med_nodes": len(nodes), "outliers_relocated": len(outliers),
               "node_sizes": " ".join(str(len(idx)) for idx in nodes)}
    (out_dir / "med_summary.txt").write_text("".join(f"{k}\t{v}\n" for k, v in summary.items()))
    print(f"MED: {len(nodes)} nodes from {len(names)} genomes ({matrix.shape[1]} variable of {length} columns, "
          f"m={m} d={d} M={M}, {len(outliers)} outliers relocated)")
    return len(nodes)


def _self_test():
    try:
        from steps.niche_metadata import load_membership
    except ImportError:
        from niche_metadata import load_membership
    rng = np.random.default_rng(0)
    base = rng.choice(list("ACGT"), 300)
    seqs = {"outgroup": "".join(rng.choice(list("ACGT"), 300))}
    for group, allele in zip("ABC", "ACG"):
        for i in range(6):
            s = base.copy()
            s[10] = allele                        # three-way column separating the groups
            s[20] = "T" if group == "A" else "A"  # A-specific SNP
            s[30] = "T" if group == "B" else "A"  # B-specific SNP
            s[100 + 10 * "ABC".index(group) + i] = "-" if i == 2 else "N" if i == 3 else "C"  # private noise
            seqs[f"{group}{i}"] = "".join(s)
    # X is an A genome with a private allele at the three-way column: alone at the root split, so an
    # outlier (1 < M) that must be relocated to A.
    s = np.array(list(seqs["A0"]))
    s[10], s[100] = "T", base[100]
    seqs["X"] = "".join(s)
    with tempfile.TemporaryDirectory() as tmp:
        fasta, clone_csv = Path(tmp) / "aln.fasta", Path(tmp) / "clone_groups.csv"
        fasta.write_text("".join(f">{k} desc\n{v}\n" for k, v in seqs.items()))
        clone_csv.write_text("representative,member\nB1,B1\nB1,B1dup\n")
        assert run_med(str(fasta), Path(tmp) / "med", outgroup="outgroup", clone_map=clone_csv) == 3
        chosen, ecotype_of = load_membership(Path(tmp) / "med" / "ecotype_membership.csv")
        groups = {}
        for taxon, eco in ecotype_of.items():
            groups.setdefault(eco, set()).add(taxon)
        assert chosen == "med" and "outgroup" not in ecotype_of
        assert groups == {"1": {f"A{i}" for i in range(6)} | {"X"},
                          "2": {f"B{i}" for i in range(6)} | {"B1dup"},
                          "3": {f"C{i}" for i in range(6)}}, groups
        summary = (Path(tmp) / "med" / "med_summary.txt").read_text()
        # each group of 6 has private SNPs (entropy 0.65 > m) that must not split it: children of 1 < M
        assert "outliers_relocated\t1\n" in summary and "node_sizes\t7 6 6\n" in summary, summary
    print("self-test passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("alignment_fasta", nargs="?")
    parser.add_argument("out_dir", nargs="?")
    parser.add_argument("--outgroup", help="ID of the outgroup sequence to exclude (step 4 matching)")
    parser.add_argument("--clone-map", help="clone_groups.csv, to list collapsed clones under their representative")
    parser.add_argument("-m", type=float, default=0.3, help="minimum entropy to split on (default 0.3)")
    parser.add_argument("-d", type=int, default=3, help="discriminant columns per split (default 3)")
    parser.add_argument("-M", type=int, help="minimum node size in genomes (default max(4, N/5000))")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        _self_test()
    elif args.alignment_fasta and args.out_dir:
        run_med(args.alignment_fasta, args.out_dir, args.outgroup, args.clone_map, args.m, args.d, args.M)
    else:
        parser.error("give alignment_fasta and out_dir, or --self-test")
