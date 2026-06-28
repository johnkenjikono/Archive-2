import sys
import re
from Bio import SeqIO
from Bio import Phylo

def _normalize_outgroup_id(text):
    cleaned = text.strip().split()[0]
    cleaned = re.sub(r"\.(fasta|fa|fna)$", "", cleaned, flags=re.IGNORECASE)
    return cleaned


def _resolve_outgroup_name(requested_outgroup, leaf_names):
    if not requested_outgroup:
        return None

    normalized = _normalize_outgroup_id(requested_outgroup)
    leaf_set = set(leaf_names)

    if normalized in leaf_set:
        return normalized

    if normalized.startswith("GCF_"):
        swapped = "GCA_" + normalized[4:]
        if swapped in leaf_set:
            return swapped
    elif normalized.startswith("GCA_"):
        swapped = "GCF_" + normalized[4:]
        if swapped in leaf_set:
            return swapped

    accession_match = re.search(r"(GC[AF]_\d+\.\d+)", normalized)
    accession = accession_match.group(1) if accession_match else normalized
    for leaf in leaf_names:
        if accession in str(leaf):
            return leaf

    return None


def reroot_tree_by_first_fasta(fasta_file, tree_file, output_rooted_tree, outgroup_name=None):
    """
    Root a phylogenetic tree using the first sequence in a FASTA file as outgroup.
    """
    # 1. Read the outgroup name from the top of the FASTA file
    print(f"Reading FASTA file to find outgroup: {fasta_file}")

    if outgroup_name:
        print(f"Using requested outgroup: {outgroup_name}")
    else:
        try:
            # We only need the very first record
            first_record = next(SeqIO.parse(fasta_file, "fasta"))
        except StopIteration:
            print("❌ Fasta file is empty!")
            return False

        outgroup_name = first_record.id
        print(f"✅ Found outgroup at the top of the FASTA: {outgroup_name}")

# 2. Open the unrooted tree from FastTree
#     
    print(f"Reading unrooted tree: {tree_file}")
    try:
        tree = Phylo.read(tree_file, "newick")
    except Exception as e:
        print(f"❌ Error reading tree file: {e}")
        print(f"   Make sure FastTree completed successfully: {tree_file}")
        return False

    # Get all leaf names for debugging
    leaf_names = [clade.name for clade in tree.get_terminals()]

    resolved_outgroup = _resolve_outgroup_name(outgroup_name, leaf_names)
    if not resolved_outgroup:
        print(f"❌ Outgroup '{outgroup_name}' could not be matched to any tree leaf.")
        print(f"   Tree contains {len(leaf_names)} sequences.")
        return False

    if resolved_outgroup != outgroup_name:
        print(f"ℹ️  Matched outgroup '{outgroup_name}' to tree leaf '{resolved_outgroup}'")

    # 3. Root the tree with that outgroup
    print(f"Rooting the tree with outgroup '{resolved_outgroup}'...")
    try:
        # BioPython requires the exact leaf name to root
        tree.root_with_outgroup({"name": resolved_outgroup})
    except ValueError as e:
        print(f"❌ Error rooting tree: {e}")
        print(f"   Outgroup '{resolved_outgroup}' not found in tree!")
        print(f"   Tree contains {len(leaf_names)} sequences.")
        print(f"   Tip: Sequence IDs must match exactly (spaces, special chars count)")

        # Try to help find similar names
        close_matches = [name for name in leaf_names if outgroup_name.split()[0] in str(name)]
        if close_matches:
            print(f"   Did you mean one of these?")
            for match in close_matches[:5]:
                print(f"     - {match}")
        return False
    except Exception as e:
        print(f"❌ Unexpected error while rooting tree: {e}")
        return False

    # 4. Save the new rooted tree
    print(f"Saving rerooted tree to {output_rooted_tree}")
    try:
        Phylo.write(tree, output_rooted_tree, "newick")
        print("✅ Done! Tree successfully rerooted and ready for EcoSim!")
        return True
    except Exception as e:
        print(f"❌ Error writing tree file: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) not in (4, 5):
        print("Usage: python reroot_tree.py <sorted_fasta_file> <unrooted_nwk_file> <output_rooted_nwk> [outgroup_id]")
        print("Example: python reroot_tree.py Bordetella_pertussis_sorted.fasta trees_final/Bordetella_pertussis_tree.nwk rerooted_trees/Bordetella_pertussis_rerooted_tree.nwk GCF_000217655.1.fna")
        sys.exit(1)

    fasta_file = sys.argv[1]
    tree_file = sys.argv[2]
    out_tree = sys.argv[3]
    outgroup = sys.argv[4] if len(sys.argv) == 5 else None
    reroot_tree_by_first_fasta(fasta_file, tree_file, out_tree, outgroup_name=outgroup)