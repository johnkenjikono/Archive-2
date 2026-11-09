"""
Minimum Entropy Decomposition (Eren et al. 2015, ISME J) on a core alignment,
compared with EcoSim's full-core-genome ecotypes. Standalone: not part of the pipeline.

    python steps/med.py Species/pangenome/core_gene_alignment.aln \
        Species/ecosim_output_core_gene_alignment [--min-size 2] [--min-entropy 0.0965]

The pipeline deletes the core alignment at the end unless run with --keep core_alignment.
Writes med_nodes.csv (taxon, med_node, ecotype) into the EcoSim folder and prints the
node/ecotype counts and their adjusted Rand index.

MED splits a group of sequences on the alignment column with the highest Shannon
entropy and repeats on each part until no column is above --min-entropy.
Amplicon MED relies on read abundance; here each genome is one sequence, so an
allele carried by fewer than --min-size genomes in a node is not counted as
variation (otherwise every distinct genome ends up its own node).
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser

try:
    from steps.parsing import _adjusted_rand_index
    from steps.rarefaction import FULL_CORE_NAME
except ImportError:
    from parsing import _adjusted_rand_index
    from rarefaction import FULL_CORE_NAME

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
MISSING = 0  # gaps (Panaroo's padding for absent genes), N and other ambiguity codes


def load_variable_columns(path, keep=None):
    """Names and a (genomes x variable columns) uint8 array; non-ACGT become MISSING."""
    names, rows = [], []
    with open(path) as handle:
        for header, seq in SimpleFastaParser(handle):
            name = header.split()[0]
            if keep is None or name in keep:
                names.append(name)
                rows.append(np.frombuffer(seq.upper().encode(), dtype=np.uint8))
    # ponytail: whole alignment in memory (~genomes x length bytes, a few x that while filtering);
    # stream it in column chunks if 200+ genome alignments run out of RAM.
    aln = np.vstack(rows)
    aln[~np.isin(aln, BASES)] = MISSING
    variable = sum((aln == b).any(axis=0) for b in BASES) >= 2
    return names, aln[:, variable]


def base_counts(block):
    return np.stack([(block == b).sum(axis=0) for b in BASES])  # 4 x columns


def column_entropy(block, min_size):
    """Shannon entropy (bits) of each column, ignoring alleles seen in fewer than min_size rows."""
    counts = base_counts(block)
    counts[counts < min_size] = 0
    p = counts / np.maximum(counts.sum(axis=0), 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return -np.nansum(p * np.log2(p), axis=0)


def decompose(aln, min_size=2, min_entropy=0.0965):
    """MED nodes as a list of row-index arrays."""
    final, todo = [], [np.arange(len(aln))]
    while todo:
        idx = todo.pop()
        block = aln[idx]
        entropy = column_entropy(block, min_size)
        if entropy.size == 0 or entropy.max() <= min_entropy:
            final.append(idx)
            continue
        # ponytail: one discriminant column per split (MED's -d 1); more would split faster, not differently.
        alleles = block[:, entropy.argmax()]
        kept = [b for b in BASES if (alleles == b).sum() >= min_size]  # >= 2, since entropy > 0
        groups = [list(idx[alleles == b]) for b in kept]
        consensus = [BASES[base_counts(aln[g]).argmax(axis=0)] for g in groups]
        # Rows with a gap or a rare allele here join the child whose consensus they match best.
        for row in idx[~np.isin(alleles, kept)]:
            seq = aln[row]
            seen = seq != MISSING
            best = min(range(len(groups)), key=lambda k: np.count_nonzero(seen & (seq != consensus[k])))
            groups[best].append(row)
        todo.extend(np.array(g) for g in groups)
    return sorted(final, key=len, reverse=True)


def full_core_ecotypes(ecosim_dir):
    """{taxon: ecotype} for the full core genome from ecotype_membership.csv (clones included)."""
    with open(Path(ecosim_dir) / "ecotype_membership.csv", newline="") as handle:
        rows = [r for r in csv.DictReader(handle) if r["file"] == f"{FULL_CORE_NAME}_results.xml"]
    if not rows:
        raise SystemExit("No full core genome ecotypes in ecotype_membership.csv (run without --no-full-core).")
    return {r["taxon"]: r["ecotype"] for r in rows}


def compare(alignment, ecosim_dir, min_size=2, min_entropy=0.0965):
    ecotypes = full_core_ecotypes(ecosim_dir)
    names, aln = load_variable_columns(alignment, keep=ecotypes.keys())
    missing = ecotypes.keys() - set(names)
    if missing:
        print(f"warning: {len(missing)} genomes with an ecotype are not in the alignment, e.g. {sorted(missing)[:3]}")
    nodes = {names[i]: n for n, idx in enumerate(decompose(aln, min_size, min_entropy), 1) for i in idx}

    out = Path(ecosim_dir) / "med_nodes.csv"
    with open(out, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["taxon", "med_node", "ecotype"])
        writer.writerows(sorted((t, n, ecotypes[t]) for t, n in nodes.items()))

    ari = _adjusted_rand_index(nodes, ecotypes)
    print(f"{len(nodes)} genomes, {aln.shape[1]} variable columns")
    print(f"EcoSim ecotypes: {len(set(ecotypes[t] for t in nodes))}   MED nodes: {len(set(nodes.values()))}")
    print(f"adjusted Rand index: {ari:.3f}" if ari is not None else "adjusted Rand index: n/a (< 2 genomes)")
    print(f"wrote {out}")
    return nodes, ecotypes, ari


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("alignment", help="core_gene_alignment.aln (FASTA)")
    parser.add_argument("ecosim_dir", help="folder holding ecotype_membership.csv")
    parser.add_argument("--min-size", type=int, default=2,
                        help="alleles in fewer genomes than this are ignored, so nodes have at least this many (default 2)")
    parser.add_argument("--min-entropy", type=float, default=0.0965,
                        help="stop splitting when no column's entropy is above this, in bits (default 0.0965)")
    args = parser.parse_args()
    compare(args.alignment, args.ecosim_dir, args.min_size, args.min_entropy)
