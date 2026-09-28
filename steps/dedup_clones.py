import csv
import os
import sys

import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser

try:
    from steps.reroot_tree import _resolve_outgroup_name
except ImportError:  # run as a script: python steps/dedup_clones.py
    from reroot_tree import _resolve_outgroup_name

# Per-site divergence at or below which two genomes count as clones.
DEFAULT_CLONE_THRESHOLD = 1e-5
# Columns compared at a time; most non-clones exceed the mismatch limit within the first chunk.
_CHUNK = 1 << 16


def _mismatches_within(seq, reps, limit):
    """
    Mismatch counts between seq and each array in reps, or None for reps with more
    than limit mismatches. Stops comparing a rep once it passes the limit.
    """
    counts = [0] * len(reps)
    alive = list(range(len(reps)))
    for start in range(0, len(seq), _CHUNK):
        if not alive:
            break
        chunk = seq[start:start + _CHUNK]
        still = []
        for i in alive:
            counts[i] += int(np.count_nonzero(reps[i][start:start + _CHUNK] != chunk))
            if counts[i] <= limit:
                still.append(i)
        alive = still
    alive = set(alive)
    return [counts[i] if i in alive else None for i in range(len(reps))]


def collapse_identical_sequences(input_fasta, output_fasta, clone_map_path, outgroup_id=None,
                                 threshold=DEFAULT_CLONE_THRESHOLD):
    """
    Collapse clones (near-identical aligned sequences) to one representative each.

    Two sequences are clones when their per-site divergence (differing alignment
    columns / alignment length, case-insensitive; a gap against a base counts as a
    difference) is at most threshold. threshold=0 collapses only exact duplicates.
    Groups are built greedily in input order: each sequence joins the closest
    existing representative within the threshold, or becomes a new one. Members are
    therefore within threshold of their representative, and at most 2x threshold of
    each other.

    Clones add near-zero branches to the tree and extra taxa to EcoSim without adding
    information.

    The root genome (outgroup_id, matched the same way step 6 matches it, or the
    first record when there is no outgroup id) is always kept as its own group:
    EcoSim never assigns the outgroup to an ecotype, so genomes collapsed into it
    would vanish from the results, and collapsing it away would break rerooting.
    Genomes that are clones of the root are grouped among themselves instead.

    input_fasta and output_fasta may be the same path. clone_map_path gets one
    row per input genome (representative, member, group_size, differences,
    divergence) so parsing can give collapsed genomes their representative's ecotype.
    Returns (n_input, n_representatives), or None on failure.
    """
    with open(input_fasta) as handle:
        records = [(title.split(None, 1)[0], seq) for title, seq in SimpleFastaParser(handle)]
    if not records:
        print(f"❌ No sequences found in {input_fasta}")
        return None
    length = len(records[0][1])
    if any(len(seq) != length for _, seq in records):
        print(f"❌ Sequences in {input_fasta} are not all the same length; cannot compare them.")
        return None

    ids = [seq_id for seq_id, _ in records]
    root = (_resolve_outgroup_name(outgroup_id, ids) if outgroup_id else None) or ids[0]
    # divergence <= threshold  <=>  mismatches <= threshold * length
    limit = int(threshold * length + 1e-9)

    group_of_seq = {}  # upper-cased sequence -> (representative id, mismatches); exact repeats skip comparing
    rep_arrays = []    # encoded representative sequences (ingroup only)
    rep_ids = []
    members_of = {}    # representative id -> [(member id, mismatches)], in input order
    representatives = []
    root_array = None
    for seq_id, seq in records:
        upper = seq.upper()
        if seq_id == root:
            root_array = np.frombuffer(upper.encode(), dtype=np.uint8)
            members_of[seq_id] = [(seq_id, 0)]
            representatives.append((seq_id, seq))
            continue
        if upper in group_of_seq:
            rep, diffs = group_of_seq[upper]
            members_of[rep].append((seq_id, diffs))
            continue
        arr = np.frombuffer(upper.encode(), dtype=np.uint8)
        best = None
        if limit and rep_arrays:
            counts = _mismatches_within(arr, rep_arrays, limit)
            hits = [(c, i) for i, c in enumerate(counts) if c is not None]
            best = min(hits) if hits else None
        if best is None:
            group_of_seq[upper] = (seq_id, 0)
            rep_arrays.append(arr)
            rep_ids.append(seq_id)
            members_of[seq_id] = [(seq_id, 0)]
            representatives.append((seq_id, seq))
        else:
            rep = rep_ids[best[1]]
            group_of_seq[upper] = (rep, best[0])
            members_of[rep].append((seq_id, best[0]))

    # Write to a temp file first so the input can be replaced in place.
    tmp_out = f"{output_fasta}.tmp"
    with open(tmp_out, "w") as out_f:
        for seq_id, seq in representatives:
            out_f.write(f">{seq_id}\n{seq}\n")
    os.replace(tmp_out, output_fasta)

    os.makedirs(os.path.dirname(os.path.abspath(clone_map_path)), exist_ok=True)
    with open(clone_map_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["representative", "member", "group_size", "differences", "divergence"])
        for rep, _ in representatives:
            for member, diffs in members_of[rep]:
                writer.writerow([rep, member, len(members_of[rep]), diffs, f"{diffs / length:.2e}"])

    n_in, n_out = len(records), len(representatives)
    groups = [(rep, [m for m, _ in ms]) for rep, ms in members_of.items() if len(ms) > 1]
    rule = "identical" if limit == 0 else f"divergence <= {threshold:g} (<= {limit} of {length} sites)"
    print(f"✅ Collapsed {n_in} sequences into {n_out} representatives, clones = {rule} "
          f"({len(groups)} clone group(s), {n_in - n_out} sequence(s) removed; root {root} kept)")
    for rep, members in groups[:10]:
        shown = ", ".join(members[1:4]) + ("..." if len(members) > 4 else "")
        print(f"   {rep}: {len(members)} clones ({shown})")
    if len(groups) > 10:
        print(f"   ... and {len(groups) - 10} more group(s); see {clone_map_path}")
    if rep_arrays:
        near_root = [rep_ids[i] for i, c in enumerate(_mismatches_within(root_array, rep_arrays, limit))
                     if c is not None]
        for twin in near_root:
            print(f"ℹ️  {twin} is a clone of the root {root}; kept as an ingroup representative "
                  f"({len(members_of[twin])} genome(s)).")

    return n_in, n_out


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Collapse clones (near-identical aligned sequences)")
    parser.add_argument("input_fasta")
    parser.add_argument("output_fasta")
    parser.add_argument("clone_map", help="Output clone_groups.csv")
    parser.add_argument("outgroup_id", nargs="?", default=None)
    parser.add_argument("--threshold", type=float, default=DEFAULT_CLONE_THRESHOLD,
                        help=f"Max per-site divergence for clones (default: {DEFAULT_CLONE_THRESHOLD:g}; 0 = identical only)")
    args = parser.parse_args()

    if collapse_identical_sequences(args.input_fasta, args.output_fasta, args.clone_map,
                                    args.outgroup_id, args.threshold) is None:
        sys.exit(1)
