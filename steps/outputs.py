"""The outputs a pipeline run can keep, shared by pipeline.py (--keep/--discard) and the UI checkboxes."""
import argparse

# What a run leaves on disk; everything else is deleted once the next step has used it.
# --keep / --discard change the set (see resolve_keep).
KEEP_ITEMS = {
    "genomes":        (False, "downloaded genome .fna files (input/)"),
    "annotations":    (False, "Bakta annotations (intermediate_bakta/)"),
    "pangenome":      (True,  "Panaroo gene presence/absence table + summary statistics (pangenome/)"),
    "core_alignment": (False, "core gene alignment + gene coordinates (pangenome/)"),
    "panaroo":        (False, "the whole Panaroo output folder (output_roary/)"),
    "rerooted_tree":  (True,  "rooted full core-genome tree (rerooted_trees/<Species>.nwk)"),
    "rarefaction":    (False, "rarefaction sub-alignments and their trees"),
    "ecosim_xml":     (True,  "raw EcoSim XML per replicate (needed to rerun steps 9-10)"),
    "ecotypes":       (True,  "ecotype membership, full-core comparison and ecotype support tables"),
    "tree_plots":     (True,  "trees colored by ecotype (PNG + iTOL colors)"),
    "pie_charts":     (True,  "pie charts of genomes per ecotype"),
}
DEFAULT_KEEP = {name for name, (default, _) in KEEP_ITEMS.items() if default}

def _parse_keep_list(raw_values):
    names = set()
    for raw in raw_values or []:
        for name in (n.strip() for n in raw.split(",")):
            if not name:
                continue
            if name == "all":
                names |= set(KEEP_ITEMS)
            elif name in KEEP_ITEMS:
                names.add(name)
            else:
                raise argparse.ArgumentTypeError(
                    f"unknown item '{name}' (choose from: all, {', '.join(KEEP_ITEMS)})")
    return names

def resolve_keep(keep_values, discard_values):
    """Defaults minus --discard, plus --keep (so '--discard all --keep rerooted_tree' works)."""
    return (DEFAULT_KEEP - _parse_keep_list(discard_values)) | _parse_keep_list(keep_values)
