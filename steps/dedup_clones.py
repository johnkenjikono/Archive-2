import csv
import os
import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

try:
    from steps.reroot_tree import _resolve_outgroup_name
except ImportError:  # run as a script: python steps/dedup_clones.py
    from reroot_tree import _resolve_outgroup_name


def collapse_identical_sequences(input_fasta, output_fasta, clone_map_path, outgroup_id=None):
    """
    Collapse identical aligned sequences (clones) to one representative each.

    Identical core genomes add zero-length branches to the tree and extra taxa to
    EcoSim without adding information. Each clone group keeps its first-seen member.

    The root genome (outgroup_id, matched the same way step 6 matches it, or the
    first record when there is no outgroup id) is always kept as its own group:
    EcoSim never assigns the outgroup to an ecotype, so genomes collapsed into it
    would vanish from the results, and collapsing it away would break rerooting.
    Genomes identical to the root form a separate ingroup clone group instead.

    input_fasta and output_fasta may be the same path. clone_map_path gets one
    row per input genome (representative, member, group_size) so parsing can give
    collapsed genomes their representative's ecotype.
    Returns (n_input, n_representatives), or None on failure.
    """
    with open(input_fasta) as handle:
        records = [(title.split(None, 1)[0], seq) for title, seq in SimpleFastaParser(handle)]
    if not records:
        print(f"❌ No sequences found in {input_fasta}")
        return None

    ids = [seq_id for seq_id, _ in records]
    root = (_resolve_outgroup_name(outgroup_id, ids) if outgroup_id else None) or ids[0]

    rep_for_seq = {}   # sequence -> representative id (ingroup only)
    members_of = {}    # representative id -> [member ids], in input order
    representatives = []
    root_seq = None
    for seq_id, seq in records:
        if seq_id == root:
            root_seq = seq.upper()
            members_of[seq_id] = [seq_id]
            representatives.append((seq_id, seq))
            continue
        rep = rep_for_seq.setdefault(seq.upper(), seq_id)
        if rep == seq_id:
            members_of[seq_id] = [seq_id]
            representatives.append((seq_id, seq))
        else:
            members_of[rep].append(seq_id)

    # Write to a temp file first so the input can be replaced in place.
    tmp_out = f"{output_fasta}.tmp"
    with open(tmp_out, "w") as out_f:
        for seq_id, seq in representatives:
            out_f.write(f">{seq_id}\n{seq}\n")
    os.replace(tmp_out, output_fasta)

    os.makedirs(os.path.dirname(os.path.abspath(clone_map_path)), exist_ok=True)
    with open(clone_map_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["representative", "member", "group_size"])
        for rep, _ in representatives:
            for member in members_of[rep]:
                writer.writerow([rep, member, len(members_of[rep])])

    n_in, n_out = len(records), len(representatives)
    groups = [(rep, m) for rep, m in members_of.items() if len(m) > 1]
    print(f"✅ Collapsed {n_in} sequences into {n_out} unique haplotypes "
          f"({len(groups)} clone group(s), {n_in - n_out} sequence(s) removed; root {root} kept)")
    for rep, members in groups[:10]:
        shown = ", ".join(members[1:4]) + ("..." if len(members) > 4 else "")
        print(f"   {rep}: {len(members)} identical ({shown})")
    if len(groups) > 10:
        print(f"   ... and {len(groups) - 10} more group(s); see {clone_map_path}")
    twin = rep_for_seq.get(root_seq)
    if twin:
        print(f"ℹ️  {twin} is identical to the root {root}; kept as an ingroup representative "
              f"({len(members_of[twin])} genome(s)).")

    return n_in, n_out


if __name__ == "__main__":
    if len(sys.argv) not in (4, 5):
        print("Usage: python dedup_clones.py <input_fasta> <output_fasta> <clone_groups.csv> [outgroup_id]")
        sys.exit(1)

    outgroup = sys.argv[4] if len(sys.argv) == 5 else None
    if collapse_identical_sequences(sys.argv[1], sys.argv[2], sys.argv[3], outgroup) is None:
        sys.exit(1)
