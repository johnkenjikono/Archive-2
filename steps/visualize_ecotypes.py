"""Draw EcoSim phylogenies with each ecotype colored (PNG + iTOL TREE_COLORS file).

Uses the tree embedded in each EcoSim XML by default: step 7b gives every
rarefaction replicate its own tree, and step 8 deletes the tree files after EcoSim.

For the full core genome result it also draws two figures that add rarefaction
support, the share of random gene samples that recovered each ecotype exactly:
*_summary.png (one wedge per ecotype) and *_tracks.png (the whole tree with data tracks).
"""
import csv
import math
import os
import re
import glob
from io import StringIO

import matplotlib
matplotlib.use("Agg")  # headless (Modal containers have no display)
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Polygon, Rectangle
from Bio import Phylo

try:
    from steps.parsing import parse_ecosim_xml, load_clone_map, FULL_CORE_XML
    from steps.download_outgroup import OUTGROUP_INFO, outgroup_label as read_outgroup_label
except ImportError:  # run as a script: python steps/visualize_ecotypes.py
    from parsing import parse_ecosim_xml, load_clone_map, FULL_CORE_XML
    from download_outgroup import OUTGROUP_INFO, outgroup_label as read_outgroup_label

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


def plot_ecotype_tree(xml_path, out_png, tree_path=None, clone_map_path=None, itol_path=None,
                      outgroup_label=None):
    """
    Draw the phylogeny with each EcoSim ecotype colored.

    Ecotypes that form a clade are colored at their common ancestor, so the whole
    clade shares one color; non-monophyletic ecotypes color only their tips.
    Single-member ecotypes are gray and the outgroup is black. Leaf labels show
    the ecotype number and how many clones were collapsed into the leaf.
    outgroup_label (e.g. "Listeria innocua ATCC 33090 (GCF_028596125.1)") names the
    outgroup tip and legend entry; the pipeline renames that genome to plain "outgroup".
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
            if outgroup_label:
                return f"{outgroup_label} [outgroup]"
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
        handles.append(Patch(color=OUTGROUP_COLOR,
                             label=f"outgroup: {outgroup_label}" if outgroup_label else "outgroup"))
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


# --- Rarefaction support figures (full core genome result only) -------------------------------
# A support cell is the share of random gene samples whose EcoSim run produced an ecotype with
# exactly this ecotype's members; parsing.compare_to_full_core writes them to ecotype_support.csv.
SUPPORT_CSV = "ecotype_support.csv"
MIXED_COLOR = "#52514e"  # branches spanning several ecotypes
SUPPORT_CMAP = LinearSegmentedColormap.from_list(
    "support", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
# Tips end at different depths, so a faint dotted line runs from each tip to its label or strip.
_LEADER = dict(colors="#d3d2cc", linewidths=0.45, linestyles=(0, (1, 1.6)))


def load_rarefaction_support(ecosim_dir):
    """
    What the support figures need from an EcoSim output folder, or None without ecotype_support.csv:
    {"support": {ecotype: {gene_count: fraction}}, "gene_counts": [...], "replicates": {gene_count: n},
    "pool": genes the samples were drawn from, or None for runs made before rarefaction recorded it}.
    """
    path = os.path.join(ecosim_dir, SUPPORT_CSV)
    if not os.path.exists(path):
        return None
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return None
    gene_counts = sorted(int(k[len("support_g"):]) for k in rows[0] if k.startswith("support_g"))
    support = {int(r["ecotype"]): {g: float(r[f"support_g{g}"] or 0) for g in gene_counts} for r in rows}
    replicates = {}
    by_count = os.path.join(ecosim_dir, "rarefaction_by_gene_count.csv")
    if os.path.exists(by_count):
        with open(by_count, newline="") as f:
            replicates = {int(r["gene_count"]): int(r["replicates"]) for r in csv.DictReader(f)}
    pool = None
    manifest = os.path.join(ecosim_dir, "rarefaction_genes.csv")
    if os.path.exists(manifest):
        with open(manifest, newline="") as f:
            pool = next((int(r["sampled_from"]) for r in csv.DictReader(f)
                         if r.get("trial") and r.get("sampled_from")), None)
    return {"support": support, "gene_counts": gene_counts, "replicates": replicates, "pool": pool}


def _support_header(rarefaction):
    """The support columns' header, with the run's own sample count and gene pool where known."""
    reps = {rarefaction["replicates"].get(g) for g in rarefaction["gene_counts"]}
    # Equal counts: "5 random". Unequal (e.g. --each-gene) or unknown: per-column n in the ticks.
    n = f"{reps.pop()} random " if len(reps) == 1 and None not in reps else ""
    drawn = f", drawn from {rarefaction['pool']:,} core genes" if rarefaction["pool"] else ""
    return f"% of {n}gene samples that\nrecovered this exact ecotype\n(columns: genes per sample{drawn})"


def _support_ticks(rarefaction):
    reps = rarefaction["replicates"]
    uniform = len({reps.get(g) for g in rarefaction["gene_counts"]}) == 1
    return [str(g) if uniform or g not in reps else f"{g}\nn={reps[g]}" for g in rarefaction["gene_counts"]]


def _layout(root):
    """{id(clade): x} as distance from root and {id(clade): y} as tip index, internal nodes midway."""
    xs, ys = {id(root): 0.0}, {}
    for clade in root.find_clades(order="preorder"):
        for child in clade.clades:
            xs[id(child)] = xs[id(clade)] + (child.branch_length or 0.0)
    for i, tip in enumerate(root.get_terminals()):
        ys[id(tip)] = float(i)
    for clade in root.find_clades(order="postorder"):
        if clade.clades:
            child_ys = [ys[id(child)] for child in clade.clades]
            ys[id(clade)] = (min(child_ys) + max(child_ys)) / 2
    return xs, ys


def _runs(values):
    """[value, first index, last index] for each run of equal consecutive values."""
    out = []
    for i, v in enumerate(values):
        if out and out[-1][0] == v:
            out[-1][2] = i
        else:
            out.append([v, i, i])
    return out


def _spread(targets, gap, lo, hi):
    """Nudge sorted label positions at least `gap` apart, staying within [lo, hi]."""
    pos = list(targets)
    for i in range(1, len(pos)):
        pos[i] = max(pos[i], pos[i - 1] + gap)
    pos[-1] = min(pos[-1], hi)
    for i in range(len(pos) - 2, -1, -1):
        pos[i] = min(pos[i], pos[i + 1] - gap)
    pos[0] = max(pos[0], lo)
    return pos


def _scale_bar_length(xmax):
    """A round branch length, about a sixth of the tree's depth."""
    target = xmax / 6
    base = 10 ** math.floor(math.log10(target))
    return max(m * base for m in (1, 2, 5) if m * base <= target)


def _genomes_text(n):
    return f"{n} genome" + ("" if n == 1 else "s")


def _support_figure_data(xml_path, clone_map_path):
    """Everything both support figures draw from: the ingroup tree, ecotypes, colors and genome counts."""
    result = parse_ecosim_xml(xml_path)
    if result is None or not result["ecotypes"] or not result["tree"]:
        return None
    tree = Phylo.read(StringIO(result["tree"]), "newick")
    tips = [t for t in tree.get_terminals() if t.name != result["outgroup"]]
    root = tree.common_ancestor(tips) if len(tips) > 1 else tips[0]
    ecotype_of = {m: eco["number"] for eco in result["ecotypes"] for m in eco["members"]}
    clone_map = load_clone_map(clone_map_path)
    palette = _ecotype_colors(result)
    sets = {}  # id(clade) -> ecotypes of its tips (None = unassigned); per call, as ids are reused
    for clade in root.find_clades(order="postorder"):
        sets[id(clade)] = ({ecotype_of.get(clade.name)} if clade.is_terminal()
                           else set().union(*(sets[id(child)] for child in clade.clades)))

    def color_of(clade):
        found = sets[id(clade)]
        if len(found) != 1:
            return MIXED_COLOR
        eco = next(iter(found))
        return UNASSIGNED_COLOR if eco is None else palette.get(eco, SINGLETON_COLOR)

    genomes = {eco["number"]: sum(len(clone_map.get(m, [m])) for m in eco["members"]) for eco in result["ecotypes"]}
    unique = {eco["number"]: len(eco["members"]) for eco in result["ecotypes"]}
    return {"result": result, "root": root, "sets": sets, "color_of": color_of, "palette": palette,
            "clone_map": clone_map, "ecotype_of": ecotype_of, "genomes": genomes, "unique": unique}


def _support_titles(fig, data, outgroup_label, what):
    height = fig.get_figheight()
    title = (f"{len(data['result']['ecotypes'])} ecotypes, {sum(data['genomes'].values())} genomes "
             f"(full core genome)")
    fig.text(0.015, 1 - 0.12 / height, title, ha="left", va="top", fontsize=13, weight="bold")
    fig.text(0.015, 1 - 0.42 / height, f"{what}  Rooted with outgroup {outgroup_label or 'outgroup'}; "
             "outgroup branch not drawn.", ha="left", va="top", fontsize=9, color=MIXED_COLOR)


def plot_ecotype_summary(xml_path, out_png, rarefaction, clone_map_path=None, outgroup_label=None):
    """
    One wedge per ecotype (height scaled to its distinct genomes) beside its rarefaction support.
    A compact view of the full-core result for talks; an ecotype that is split across the tree
    shows as numbered parts, with its support row on the first.
    """
    data = _support_figure_data(xml_path, clone_map_path)
    if data is None:
        return False
    root, sets, color_of = data["root"], data["sets"], data["color_of"]
    xs, _ = _layout(root)
    units = []  # the largest clades whose tips all share one ecotype, in tree order

    def walk(clade):
        if len(sets[id(clade)]) == 1:
            units.append(clade)
        else:
            for child in clade.clades:
                walk(child)
    walk(root)
    unit_ids = {id(u) for u in units}
    inside = {id(d) for u in units for d in u.find_clades()}
    ys, heights, y = {}, {}, 0.0
    for unit in units:
        h = 0.5 + 0.06 * len(unit.get_terminals())
        ys[id(unit)], heights[id(unit)] = y + h / 2, h
        y += h + 0.35

    def set_y(clade):
        if id(clade) not in unit_ids:
            child_ys = [set_y(child) for child in clade.clades]
            ys[id(clade)] = (min(child_ys) + max(child_ys)) / 2
        return ys[id(clade)]
    set_y(root)
    ytot = y
    xmax = max(xs[id(t)] for t in root.get_terminals()) or 1.0

    top_in, bottom_in = 1.75, 0.55
    fig_h = top_in + bottom_in + 0.43 * len(units)
    fig = plt.figure(figsize=(14, fig_h))
    rect_y, rect_h = bottom_in / fig_h, (fig_h - top_in - bottom_in) / fig_h
    ax = fig.add_axes([0.03, rect_y, 0.44, rect_h])
    segs, cols = [], []
    for clade in root.find_clades(order="preorder"):
        if id(clade) in inside or not clade.clades:
            continue
        child_ys = [ys[id(child)] for child in clade.clades]
        segs.append([(xs[id(clade)], min(child_ys)), (xs[id(clade)], max(child_ys))]); cols.append(MIXED_COLOR)
        for child in clade.clades:
            segs.append([(xs[id(clade)], ys[id(child)]), (xs[id(child)], ys[id(child)])])
            cols.append(color_of(child) if id(child) in unit_ids else MIXED_COLOR)
    lines = LineCollection(segs, colors=cols, linewidths=1.6)
    lines.set_capstyle("round")
    ax.add_collection(lines)

    parts = {}
    for unit in units:
        parts.setdefault(next(iter(sets[id(unit)])), []).append(unit)
    rows = []
    for unit in units:
        eco = next(iter(sets[id(unit)]))
        x0, yc, h = xs[id(unit)], ys[id(unit)], heights[id(unit)]
        x1 = max(xs[id(t)] for t in unit.get_terminals())
        if unit.is_terminal():
            ax.plot([x0], [yc], "o", ms=6, color=color_of(unit), mec="white", mew=1.5)
        else:
            ax.add_patch(Polygon([(x0, yc), (x1, yc - 0.45 * h), (x1, yc + 0.45 * h)], closed=True,
                                 facecolor=color_of(unit), edgecolor="white", linewidth=1.5))
        k = parts[eco].index(unit)
        rows.append((yc, eco, f"  (part {k + 1} of {len(parts[eco])})" if len(parts[eco]) > 1 else "", k == 0))
    bar = _scale_bar_length(xmax)
    ax.plot([0, bar], [ytot + 0.25, ytot + 0.25], color="black", lw=1.5)
    ax.text(bar / 2, ytot + 0.45, f"{bar:g} subst./site", ha="center", va="top", fontsize=8, color=MIXED_COLOR)
    ax.set_xlim(-xmax * 0.02, xmax * 1.02)
    ax.set_ylim(ytot + 1.0, -0.4)
    ax.axis("off")

    labels = fig.add_axes([0.47, rect_y, 0.25, rect_h], sharey=ax)
    labels.axis("off")
    for yc, eco, suffix, _ in rows:
        if eco is None:
            labels.text(0.02, yc, "unassigned" + suffix, fontsize=9, va="center", color=MIXED_COLOR)
            continue
        uniq = data["unique"][eco]
        detail = _genomes_text(data["genomes"][eco]) + (f" · {uniq} unique" if uniq != data["genomes"][eco] else "")
        labels.text(0.02, yc, f"E{eco}", fontsize=11, weight="bold", va="center")
        labels.text(0.17, yc, detail + suffix, fontsize=9, va="center", color=MIXED_COLOR)

    gene_counts = rarefaction["gene_counts"]
    grid = fig.add_axes([0.73, rect_y, 0.24, rect_h], sharey=ax)
    for yc, eco, _, first in rows:
        if not first or eco not in rarefaction["support"]:
            continue
        for j, g in enumerate(gene_counts):
            v = rarefaction["support"][eco][g]
            grid.add_patch(Rectangle((j + 0.06, yc - 0.42), 0.88, 0.84, facecolor=SUPPORT_CMAP(v), edgecolor="white"))
            grid.text(j + 0.5, yc, f"{v:.0%}", ha="center", va="center", fontsize=8,
                      color="white" if v >= 0.55 else "black")
    grid.set_xlim(0, len(gene_counts))
    grid.set_xticks([j + 0.5 for j in range(len(gene_counts))])
    grid.set_xticklabels(_support_ticks(rarefaction))
    grid.xaxis.set_ticks_position("top")
    grid.tick_params(length=0, labelsize=9, colors=MIXED_COLOR)
    grid.set_xlabel(_support_header(rarefaction), fontsize=9, color=MIXED_COLOR, labelpad=8, linespacing=1.3)
    grid.xaxis.set_label_position("top")
    grid.set_yticks([])
    for spine in grid.spines.values():
        spine.set_visible(False)

    _support_titles(fig, data, outgroup_label, "One wedge per ecotype, its height scaled to distinct genomes.")
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return True


def plot_ecotype_tracks(xml_path, out_png, rarefaction, clone_map_path=None, outgroup_label=None):
    """
    The whole tree without tip names, beside tracks for ecotype, clones collapsed into each tip and
    rarefaction support. A paper-style figure; genome names are in ecotype_membership.csv.
    """
    data = _support_figure_data(xml_path, clone_map_path)
    if data is None:
        return False
    root, color_of, ecotype_of = data["root"], data["color_of"], data["ecotype_of"]
    tips = root.get_terminals()
    n_tips = len(tips)
    xs, ys = _layout(root)
    xmax = max(xs[id(t)] for t in tips) or 1.0
    extra = {t.name: len(data["clone_map"].get(t.name, [t.name])) - 1 for t in tips}

    top_in, bottom_in = 1.8, 0.55
    fig_h = max(6.0, top_in + bottom_in + 0.035 * n_tips)
    fig = plt.figure(figsize=(15, fig_h))
    grid_spec = fig.add_gridspec(1, 5, width_ratios=[6.5, 0.35, 1.5, 1.1, 2.4], wspace=0.04, left=0.02, right=0.98,
                                 top=1 - top_in / fig_h, bottom=bottom_in / fig_h)
    ax_tree, ax_strip, ax_label, ax_clone, ax_grid = (fig.add_subplot(grid_spec[0, i]) for i in range(5))

    segs, cols = [], []
    for clade in root.find_clades(order="preorder"):
        if not clade.clades:
            continue
        child_ys = [ys[id(child)] for child in clade.clades]
        segs.append([(xs[id(clade)], min(child_ys)), (xs[id(clade)], max(child_ys))]); cols.append(color_of(clade))
        for child in clade.clades:
            segs.append([(xs[id(clade)], ys[id(child)]), (xs[id(child)], ys[id(child)])]); cols.append(color_of(child))
    lines = LineCollection(segs, colors=cols, linewidths=0.9)
    lines.set_capstyle("round")
    ax_tree.add_collection(lines)
    ax_tree.add_collection(LineCollection([[(xs[id(t)], ys[id(t)]), (xmax * 1.01, ys[id(t)])] for t in tips], **_LEADER))
    bar = _scale_bar_length(xmax)
    ax_tree.plot([0, bar], [n_tips + 1.5, n_tips + 1.5], color="black", lw=1.5)
    ax_tree.text(bar / 2, n_tips + 2.5, f"{bar:g} subst./site", ha="center", va="top", fontsize=8, color=MIXED_COLOR)
    ax_tree.set_xlim(-xmax * 0.01, xmax * 1.01)

    blocks = _runs([ecotype_of.get(t.name) for t in tips])
    gene_counts = rarefaction["gene_counts"]
    for eco, i0, i1 in blocks:
        strip_color = UNASSIGNED_COLOR if eco is None else data["palette"].get(eco, SINGLETON_COLOR)
        ax_strip.add_patch(Rectangle((0, i0 - 0.5), 1, i1 - i0 + 1, facecolor=strip_color, edgecolor="white", linewidth=1.2))
        for j, g in enumerate(gene_counts):
            v = rarefaction["support"].get(eco, {}).get(g, 0.0)
            ax_grid.add_patch(Rectangle((j, i0 - 0.5), 1, i1 - i0 + 1, facecolor=SUPPORT_CMAP(v),
                                        edgecolor="white", linewidth=1.0))
    # Keep block labels about one text line apart, with a short connector back to the block.
    row_in = fig_h * (1 - (top_in + bottom_in) / fig_h) / (n_tips + 4)
    mids = [(i0 + i1) / 2 for _, i0, i1 in blocks]
    for (eco, _, _), ym, yl in zip(blocks, mids, _spread(mids, 0.135 / row_in, 0, n_tips + 1)):
        ax_label.plot([0, 0.05, 0.09], [ym, ym, yl], color=SINGLETON_COLOR, lw=0.7)
        if eco is None:
            ax_label.text(0.11, yl, "unassigned", fontsize=8.5, va="center", color=MIXED_COLOR)
            continue
        ax_label.text(0.11, yl, f"E{eco}", fontsize=9.5, weight="bold", va="center")
        ax_label.text(0.36, yl, _genomes_text(data["genomes"][eco]), fontsize=8.5, va="center", color=MIXED_COLOR)
    ax_clone.barh([ys[id(t)] for t in tips], [extra[t.name] for t in tips], height=0.8, color=SINGLETON_COLOR)
    ax_clone.set_xlim(0, max(extra.values()) + 1)

    for axis in (ax_tree, ax_strip, ax_label, ax_clone, ax_grid):
        axis.set_ylim(n_tips + 3, -1)
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_visible(False)
    for axis in (ax_tree, ax_strip, ax_label):
        axis.set_xticks([])
    ax_strip.set_xlim(0, 1)
    ax_label.set_xlim(0, 1)
    ax_clone.tick_params(axis="x", labelsize=8, colors=MIXED_COLOR, length=0)
    ax_clone.set_title("Clones collapsed\ninto each tip", fontsize=9, color=MIXED_COLOR)
    ax_grid.set_xlim(0, len(gene_counts))
    ax_grid.set_xticks([j + 0.5 for j in range(len(gene_counts))])
    ax_grid.set_xticklabels(_support_ticks(rarefaction), fontsize=9, color=MIXED_COLOR)
    ax_grid.tick_params(length=0)
    ax_grid.xaxis.set_ticks_position("top")
    ax_grid.set_title(_support_header(rarefaction), fontsize=9, color=MIXED_COLOR, pad=16, linespacing=1.3)
    scale = plt.cm.ScalarMappable(cmap=SUPPORT_CMAP, norm=plt.Normalize(0, 1))
    bar_axis = fig.colorbar(scale, ax=ax_grid, orientation="horizontal", fraction=0.03, pad=0.01,
                            format=lambda v, _: f"{v:.0%}")
    bar_axis.ax.tick_params(labelsize=8, colors=MIXED_COLOR, length=0)
    bar_axis.outline.set_visible(False)

    _support_titles(fig, data, outgroup_label, "Tree with data tracks; genome names are in ecotype_membership.csv.")
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return True


def plot_folder(ecosim_dir, output_dir, clone_map_path=None, gene_count="max", trees=True, pies=True,
                outgroup_label=None, summary=True, tracks=True):
    """
    Plot EcoSim XMLs in ecosim_dir into output_dir: a tree (PNG + *_itol_colors.txt)
    and/or a pie chart (*_pie.png) per XML. outgroup_label names the outgroup on the trees.
    For the full core genome result, when step 9 wrote ecotype_support.csv, summary and
    tracks add *_summary.png and *_tracks.png, which show each ecotype's rarefaction support.

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

    rarefaction = None
    if (summary or tracks) and full_core in xml_files:
        rarefaction = load_rarefaction_support(ecosim_dir)
        if rarefaction is None:
            print(f"No {SUPPORT_CSV} (no rarefaction replicates?); skipping the summary and tracks figures.")

    plotted = 0
    for xml_file in xml_files:
        stem = os.path.splitext(os.path.basename(xml_file))[0]
        try:
            ok = True
            if trees:
                ok = plot_ecotype_tree(xml_file, os.path.join(output_dir, f"{stem}.png"),
                                       clone_map_path=clone_map_path,
                                       itol_path=os.path.join(output_dir, f"{stem}_itol_colors.txt"),
                                       outgroup_label=outgroup_label)
            if pies:
                ok = plot_ecotype_pie(xml_file, os.path.join(output_dir, f"{stem}_pie.png"),
                                      clone_map_path=clone_map_path) and ok
            if rarefaction and xml_file == full_core:
                if summary:
                    ok = plot_ecotype_summary(xml_file, os.path.join(output_dir, f"{stem}_summary.png"), rarefaction,
                                              clone_map_path=clone_map_path, outgroup_label=outgroup_label) and ok
                if tracks:
                    ok = plot_ecotype_tracks(xml_file, os.path.join(output_dir, f"{stem}_tracks.png"), rarefaction,
                                             clone_map_path=clone_map_path, outgroup_label=outgroup_label) and ok
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
    parser.add_argument("--no-summary", action="store_true",
                        help="Skip the one-wedge-per-ecotype summary with rarefaction support (full core result)")
    parser.add_argument("--no-tracks", action="store_true",
                        help="Skip the tree-with-tracks figure with rarefaction support (full core result)")
    parser.add_argument("--outgroup-label", default=None,
                        help=f"Name to show for the outgroup (default: read from {OUTGROUP_INFO} in the species folder)")
    args = parser.parse_args()

    folder = args.input if os.path.isdir(args.input) else os.path.dirname(os.path.abspath(args.input))
    clone_map = args.clone_map or os.path.join(folder, "clone_groups.csv")
    output = args.output or os.path.join(folder, "ecotype_plots")
    # The species folder is the parent of ecosim_output_*.
    og_label = args.outgroup_label or read_outgroup_label(
        os.path.join(os.path.dirname(os.path.abspath(folder)), OUTGROUP_INFO))

    if os.path.isdir(args.input):
        gc = None if args.gene_count == "all" else ("max" if args.gene_count == "max" else int(args.gene_count))
        n = plot_folder(args.input, output, clone_map, gc, trees=not args.no_trees, pies=not args.no_pies,
                        outgroup_label=og_label, summary=not args.no_summary, tracks=not args.no_tracks)
        print(f"Saved {n} plot(s) to {output}")
    else:
        out_png = output if output.endswith(".png") else os.path.join(
            output, os.path.splitext(os.path.basename(args.input))[0] + ".png")
        if not args.no_trees:
            plot_ecotype_tree(args.input, out_png, args.tree, clone_map,
                              itol_path=os.path.splitext(out_png)[0] + "_itol_colors.txt",
                              outgroup_label=og_label)
            print(f"Saved {out_png}")
        if not args.no_pies:
            pie_png = os.path.splitext(out_png)[0] + "_pie.png"
            plot_ecotype_pie(args.input, pie_png, clone_map)
            print(f"Saved {pie_png}")
        support = load_rarefaction_support(folder) if os.path.basename(args.input) == FULL_CORE_XML else None
        for flag, suffix, plot in ((args.no_summary, "_summary.png", plot_ecotype_summary),
                                   (args.no_tracks, "_tracks.png", plot_ecotype_tracks)):
            if support and not flag:
                extra_png = os.path.splitext(out_png)[0] + suffix
                plot(args.input, extra_png, support, clone_map, outgroup_label=og_label)
                print(f"Saved {extra_png}")
