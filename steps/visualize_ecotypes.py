"""Draw EcoSim phylogenies with each ecotype colored (PNG + iTOL TREE_COLORS file).

Uses the tree embedded in each EcoSim XML by default: step 7b gives every
rarefaction replicate its own tree, and step 8 deletes the tree files after EcoSim.
"""
import os
import re
import glob
from io import StringIO

import matplotlib
matplotlib.use("Agg")  # headless (Modal containers have no display)
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.patches import Patch
from Bio import Phylo

try:
    from steps.parsing import parse_ecosim_xml, load_clone_map, FULL_CORE_XML
except ImportError:  # run as a script: python steps/visualize_ecotypes.py
    from parsing import parse_ecosim_xml, load_clone_map, FULL_CORE_XML

OUTGROUP_COLOR = "#000000"
SINGLETON_COLOR = "#8c8c8c"
UNASSIGNED_COLOR = "#d0d0d0"
_GENE_COUNT_RE = re.compile(r"_g(\d+)_t(\d+)")


def _ecotype_palette(n):
    """n distinct hex colors, skipping the grays reserved for singletons/outgroup."""
    # Saturated tab20 shades first (== tab10), then the light shades, so sister
    # ecotypes don't get two tints of the same hue until more than ~9 are needed.
    tab20 = [mcolors.to_hex(plt.get_cmap("tab20")(i)) for i in range(20)]
    base = tab20[0::2] + tab20[1::2]
    for name in ("tab20b", "tab20c"):
        cmap = plt.get_cmap(name)
        base.extend(mcolors.to_hex(cmap(i)) for i in range(cmap.N))
    base = [c for c in base if len({c[1:3], c[3:5], c[5:7]}) > 1]  # drop pure grays
    if n <= len(base):
        return base[:n]
    hsv = plt.get_cmap("hsv")
    return [mcolors.to_hex(hsv(i / n)) for i in range(n)]


def _ecotype_colors(result):
    """{ecotype number: color} for multi-member ecotypes; single-member ones use SINGLETON_COLOR."""
    multi = [eco for eco in result["ecotypes"] if len(eco["members"]) > 1]
    return dict(zip((eco["number"] for eco in multi), _ecotype_palette(len(multi))))


def plot_ecotype_pie(xml_path, out_png, clone_map_path=None):
    """
    Pie chart of how the genomes split into ecotypes, colored like the tree plot.

    Slices count genomes (collapsed clones expanded), so an ecotype of one unique
    sequence shared by ten genomes counts ten. Ecotypes of a single unique sequence
    are pooled into one gray slice, matching the tree plot's gray tips.
    """
    result = parse_ecosim_xml(xml_path)
    if result is None or not result["ecotypes"]:
        return False
    clone_map = load_clone_map(clone_map_path)
    palette = _ecotype_colors(result)

    def genomes(eco):
        return sum(len(clone_map.get(m, [m])) for m in eco["members"])

    slices = sorted(((genomes(eco), eco["number"]) for eco in result["ecotypes"] if eco["number"] in palette),
                    reverse=True)
    singles = [eco for eco in result["ecotypes"] if eco["number"] not in palette]
    sizes = [n for n, _ in slices]
    labels = [f"E{num} ({n})" for n, num in slices]
    colors = [palette[num] for _, num in slices]
    if singles:
        sizes.append(sum(genomes(eco) for eco in singles))
        labels.append(f"{len(singles)} single-member ecotypes ({sizes[-1]})")
        colors.append(SINGLETON_COLOR)
    total = sum(sizes)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.pie(sizes, colors=colors, startangle=90, counterclock=False,
           wedgeprops={"edgecolor": "white", "linewidth": 0.8},
           autopct=lambda pct: f"{pct:.0f}%" if pct >= 3 else "", pctdistance=0.8)
    ax.set_title(f"{result['file']}: {len(result['ecotypes'])} ecotypes, {total} genomes")
    ax.legend(labels, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=8,
              title="Ecotype (genomes)")
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_ecotype_tree(xml_path, out_png, tree_path=None, clone_map_path=None, itol_path=None):
    """
    Draw the phylogeny with each EcoSim ecotype colored.

    Ecotypes that form a clade are colored at their common ancestor, so the whole
    clade shares one color; non-monophyletic ecotypes color only their tips.
    Single-member ecotypes are gray and the outgroup is black. Leaf labels show
    the ecotype number and how many clones were collapsed into the leaf.
    """
    result = parse_ecosim_xml(xml_path)
    if result is None:
        return False
    if tree_path:
        tree = Phylo.read(tree_path, "newick")
    elif result["tree"]:
        tree = Phylo.read(StringIO(result["tree"]), "newick")
    else:
        raise ValueError(f"No tree in {result['file']} and no tree_path given.")

    clone_map = load_clone_map(clone_map_path)
    leaves = {leaf.name: leaf for leaf in tree.get_terminals()}
    ecotype_of = {m: eco["number"] for eco in result["ecotypes"] for m in eco["members"]}
    missing = [m for m in ecotype_of if m not in leaves]
    if missing:
        print(f"⚠️  {len(missing)} ecotype member(s) not in tree, e.g. {missing[:3]}")

    multi = [eco for eco in result["ecotypes"] if len(eco["members"]) > 1]
    palette = _ecotype_colors(result)

    itol_rows = []
    for eco in result["ecotypes"]:
        members = [leaves[m] for m in eco["members"] if m in leaves]
        if not members:
            continue
        color = palette.get(eco["number"], SINGLETON_COLOR)
        mrca = tree.common_ancestor(members) if len(members) > 1 else members[0]
        if {t.name for t in mrca.get_terminals()} <= set(eco["members"]):
            mrca.color = color  # Phylo.draw passes a clade's color down to its descendants
            if len(members) > 1:
                itol_rows.append(f"{members[0].name}|{members[-1].name}\tclade\t{color}\tnormal\t2")
        else:
            for leaf in members:
                leaf.color = color
        itol_rows.extend(f"{leaf.name}\trange\t{color}\tEcotype {eco['number']}" for leaf in members)

    outgroup = result["outgroup"]
    if outgroup in leaves:
        leaves[outgroup].color = OUTGROUP_COLOR

    def label_for(clade):
        if not clade.is_terminal():
            return ""
        name = clade.name
        extra = len(clone_map.get(name, [name])) - 1
        if name == outgroup:
            tag = "[outgroup]"
        elif name in ecotype_of:
            tag = f"[E{ecotype_of[name]}]"
        else:
            tag = "[unassigned]"
        return f"{name} {tag}" + (f" (+{extra} clones)" if extra > 0 else "")

    def color_for(name):
        if name == outgroup:
            return OUTGROUP_COLOR
        if name in ecotype_of:
            return palette.get(ecotype_of[name], SINGLETON_COLOR)
        return UNASSIGNED_COLOR

    label_colors = {label_for(leaf): color_for(name) for name, leaf in leaves.items()}

    fig, ax = plt.subplots(figsize=(11, max(4, 0.22 * len(leaves))))
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False,
               label_func=label_for, label_colors=label_colors)
    ax.set_title(f"{result['file']}: {len(result['ecotypes'])} ecotypes "
                 f"(npop={result['npop']}), {len(leaves)} leaves")
    ax.set_ylabel("")
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)

    handles = [Patch(color=palette[eco["number"]], label=f"E{eco['number']} ({len(eco['members'])})")
               for eco in multi[:40]]
    handles.append(Patch(color=SINGLETON_COLOR, label="single-member ecotype"))
    if outgroup in leaves:
        handles.append(Patch(color=OUTGROUP_COLOR, label="outgroup"))
    # Below the figure: tip labels run past the axes and would overlap a legend on the right.
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.0),
               ncol=min(5, len(handles)), frameon=False, fontsize=8, title="Ecotype (unique seqs)")

    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)

    if itol_path:
        with open(itol_path, "w") as handle:
            handle.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n" + "\n".join(itol_rows) + "\n")
    return True


def plot_folder(ecosim_dir, output_dir, clone_map_path=None, gene_count="max", trees=True, pies=True):
    """
    Plot EcoSim XMLs in ecosim_dir into output_dir: a tree (PNG + *_itol_colors.txt)
    and/or a pie chart (*_pie.png) per XML.

    gene_count="max" plots the full core genome result when there is one, otherwise
    the replicates with the most genes (the best-resolved demarcations); an int
    plots that gene count; None plots all.
    Returns the number of XMLs plotted.
    """
    xml_files = sorted(glob.glob(os.path.join(ecosim_dir, "*.xml")))
    full_core = os.path.join(ecosim_dir, FULL_CORE_XML)
    counts = {f: int(m.group(1)) for f in xml_files if (m := _GENE_COUNT_RE.search(os.path.basename(f)))}
    if gene_count == "max" and full_core in xml_files:
        xml_files = [full_core]
    elif gene_count == "max" and counts:
        gene_count = max(counts.values())
    if isinstance(gene_count, int):
        xml_files = [f for f in xml_files if counts.get(f) == gene_count]

    plotted = 0
    for xml_file in xml_files:
        stem = os.path.splitext(os.path.basename(xml_file))[0]
        try:
            ok = True
            if trees:
                ok = plot_ecotype_tree(xml_file, os.path.join(output_dir, f"{stem}.png"),
                                       clone_map_path=clone_map_path,
                                       itol_path=os.path.join(output_dir, f"{stem}_itol_colors.txt"))
            if pies:
                ok = plot_ecotype_pie(xml_file, os.path.join(output_dir, f"{stem}_pie.png"),
                                      clone_map_path=clone_map_path) and ok
            plotted += ok
        except Exception as e:
            print(f"❌ Failed to plot {os.path.basename(xml_file)}: {e}")
    return plotted


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Color EcoSim phylogenies by ecotype")
    parser.add_argument("input", help="EcoSim XML file, or an ecosim_output_* folder")
    parser.add_argument("--tree", default=None, help="Newick tree to draw instead of the one in the XML (single file only)")
    parser.add_argument("--clone-map", default=None,
                        help="clone_groups.csv (default: clone_groups.csv in the XML's folder if present)")
    parser.add_argument("--output", default=None, help="Output PNG or folder (default: <folder>/ecotype_plots)")
    parser.add_argument("--gene-count", default="max",
                        help="Folder mode: 'max' (default: full core genome if present, else most genes), 'all', or e.g. 20")
    parser.add_argument("--no-trees", action="store_true", help="Skip the colored tree plots")
    parser.add_argument("--no-pies", action="store_true", help="Skip the ecotype pie charts")
    args = parser.parse_args()

    folder = args.input if os.path.isdir(args.input) else os.path.dirname(os.path.abspath(args.input))
    clone_map = args.clone_map or os.path.join(folder, "clone_groups.csv")
    output = args.output or os.path.join(folder, "ecotype_plots")

    if os.path.isdir(args.input):
        gc = None if args.gene_count == "all" else ("max" if args.gene_count == "max" else int(args.gene_count))
        n = plot_folder(args.input, output, clone_map, gc, trees=not args.no_trees, pies=not args.no_pies)
        print(f"Saved {n} plot(s) to {output}")
    else:
        out_png = output if output.endswith(".png") else os.path.join(
            output, os.path.splitext(os.path.basename(args.input))[0] + ".png")
        if not args.no_trees:
            plot_ecotype_tree(args.input, out_png, args.tree, clone_map,
                              itol_path=os.path.splitext(out_png)[0] + "_itol_colors.txt")
            print(f"Saved {out_png}")
        if not args.no_pies:
            pie_png = os.path.splitext(out_png)[0] + "_pie.png"
            plot_ecotype_pie(args.input, pie_png, clone_map)
            print(f"Saved {pie_png}")
