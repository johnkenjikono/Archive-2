"""
Interactive ecotype tree: one self-contained HTML file, no server or internet needed.

Opens on the full core genome EcoSim tree (or the replicate with the most genes) with
each ecotype collapsed to a wedge, next to how often rarefaction subsets recover it.
Clicking an ecotype zooms into it with every strain drawn and named; "All strains"
shows the whole tree expanded. Hovering or clicking a strain shows its ecotype, the
clones collapsed into it, and whatever metadata earlier steps saved next to the XML:
  genome_metadata.csv  host / niche / country / year      (--niche)
  genome_sources.csv   BioProject, submitter, PubMed ids (--papers)
  papers.csv           paper title, link and Claude summary (--papers)
  ecotype_tree_support.csv  SH-like support of each ecotype's clade
  ecotype_support.csv  share of rarefaction subsets recovering each ecotype

  python steps/interactive_tree.py ecosim_output_dir
"""
import csv
import html
import json
import re
from io import StringIO
from pathlib import Path

from Bio import Phylo

try:
    from steps.parsing import parse_ecosim_xml, load_clone_map, FULL_CORE_XML
    from steps.visualize_ecotypes import _ecotype_colors, SINGLETON_COLOR
except ImportError:  # run as a script: python steps/interactive_tree.py
    from parsing import parse_ecosim_xml, load_clone_map, FULL_CORE_XML
    from visualize_ecotypes import _ecotype_colors, SINGLETON_COLOR

BRANCH_COLOR = "#555555"


def _pick_xml(folder):
    full = folder / FULL_CORE_XML
    if full.exists():
        return full
    by_genes = [(int(m.group(1)), -int(m.group(2)), f) for f in folder.glob("*_results.xml")
                if (m := re.search(r"_g(\d+)_t(\d+)", f.name))]
    return max(by_genes)[2] if by_genes else None


def _rows(path, key):
    if not path.exists():
        return {}
    with open(path, newline="") as f:
        return {row[key]: row for row in csv.DictReader(f)}


def _run_parameters(folder):
    path = folder.parent / "run_parameters.json"
    return json.loads(path.read_text()) if path.exists() else {}


def write_interactive_tree(ecosim_out_dir, out_html=None):
    """Write the HTML viewer. Without an EcoSim XML it draws step 6's rooted tree with no ecotypes.
    Returns its path, or None when there is no tree at all."""
    folder = Path(ecosim_out_dir).resolve()
    xml_path = _pick_xml(folder)
    result = parse_ecosim_xml(xml_path) if xml_path else None
    if result and result["tree"]:
        tree = Phylo.read(StringIO(result["tree"]), "newick")
    else:
        # ponytail: one species per run folder, so rerooted_trees/ holds just its tree
        xml_path = next(iter(sorted((folder.parent / "rerooted_trees").glob("*.nwk"))), None)
        if not xml_path:
            print(f"No EcoSim XML or rooted tree for {folder}; skipping the interactive tree.")
            return None
        print(f"No EcoSim XML with a tree in {folder}; drawing {xml_path.name} without ecotypes.")
        tree = Phylo.read(xml_path, "newick")
        # reroot_tree.py roots on the outgroup, leaving it the one tip directly under the root.
        root_tips = [c.name for c in tree.root.clades if c.is_terminal()]
        result = {"file": xml_path.name, "outgroup": root_tips[0] if len(root_tips) == 1 else None, "ecotypes": []}
    out_html = Path(out_html or folder / "ecotype_plots" / f"{xml_path.stem}_interactive.html")

    clones = load_clone_map(folder / "clone_groups.csv")
    metadata = _rows(folder / "genome_metadata.csv", "accession")
    sources = _rows(folder / "genome_sources.csv", "accession")
    papers = _rows(folder / "papers.csv", "pmid")
    support = _rows(folder / "ecotype_tree_support.csv", "ecotype")
    # ecotype_support.csv may be left over from an earlier run whose ecotype numbers mean other strains:
    # trust a row only for the full core tree and only when its members are this ecotype's members.
    members_key = {str(eco["number"]): ";".join(eco["members"]) for eco in result["ecotypes"]}
    recovery = {eco: row for eco, row in _rows(folder / "ecotype_support.csv", "ecotype").items()
                if xml_path.name == FULL_CORE_XML and row.get("members") == members_key.get(eco)}
    params = _run_parameters(folder)
    ecotype_of = {m: eco["number"] for eco in result["ecotypes"] for m in eco["members"]}
    palette = _ecotype_colors(result)

    # The outgroup's long branch would squeeze the whole tree against one edge, and it has no ecotype.
    outgroup = result["outgroup"]
    if outgroup and any(t.name == outgroup for t in tree.get_terminals()):
        tree.prune(outgroup)
        tree.root.branch_length = 0  # the old root's other branch would otherwise indent the whole tree

    depths = tree.depths()
    has_lengths = any(depths.values())
    if not has_lengths:  # no branch lengths: fall back to node depth
        depths = tree.depths(unit_branch_lengths=True)
    tips = tree.get_terminals()
    tip_index = {tip: i for i, tip in enumerate(tips)}

    def color_of(name):
        return palette.get(ecotype_of.get(name), SINGLETON_COLOR)

    # A branch takes its ecotype's color when every tip below it is in that one (multi-member) ecotype.
    # uniform[clade] is the one ecotype every tip below shares; unassigned tips never merge with anything.
    branch_color, uniform = {}, {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            uniform[clade] = ecotype_of.get(clade.name, clade)
            branch_color[clade] = palette.get(ecotype_of.get(clade.name), BRANCH_COLOR)
        else:
            below = {branch_color[c] for c in clade.clades}
            branch_color[clade] = below.pop() if len(below) == 1 else BRANCH_COLOR
            kinds = {uniform[c] for c in clade.clades}
            uniform[clade] = kinds.pop() if len(kinds) == 1 else None

    def genomes_of(names):
        return sum(len(clones.get(name, [name])) for name in names)

    gene_counts = sorted(int(k[len("support_g"):]) for k in next(iter(recovery.values()), {})
                         if k.startswith("support_g"))

    # One overview row per collapsed clade: the highest node whose tips all share an ecotype.
    nodes, rows, index_of = [], [], {}
    crowns_of = {}
    for clade in tree.find_clades(order="preorder"):
        index_of[clade] = len(nodes)
        nodes.append({"p": -1, "x": round(depths[clade], 8), "t": tip_index.get(clade, -1),
                      "c": branch_color[clade], "r": -1, "h": False})
    for clade in tree.find_clades(order="preorder"):
        node = nodes[index_of[clade]]
        for child in clade.clades:
            nodes[index_of[child]]["p"] = index_of[clade]
        parent = nodes[node["p"]] if node["p"] >= 0 else None
        if parent and (parent["h"] or parent["r"] >= 0):
            node["h"] = True
        elif uniform[clade] is not None:
            node["r"] = len(rows)
            below = clade.get_terminals()
            eco = ecotype_of.get(below[0].name)
            crowns_of[eco] = crowns_of.get(eco, 0) + 1
            rows.append({"eco": eco, "color": color_of(below[0].name),
                         "first": tip_index[below[0]], "last": tip_index[below[-1]],
                         "name": below[0].name if eco is None else None,
                         "tip_x": round(max(depths[t] for t in below), 8)})
    members_of = {eco["number"]: eco["members"] for eco in result["ecotypes"]}
    for row in rows:
        eco = row["eco"]
        members = members_of.get(eco, [row["name"]])
        row["unique"], row["genomes"] = len(members), genomes_of(members)
        row["split"] = eco is not None and crowns_of[eco] > 1
        row["sh"] = (support.get(str(eco)) or {}).get("sh_support", "")
        saved = recovery.get(str(eco))
        row["support"] = [float(saved[f"support_g{g}"]) for g in gene_counts] if saved else None

    tip_info = []
    for tip in tips:
        name = tip.name
        genomes = clones.get(name, [name])
        eco = ecotype_of.get(name)
        tip_info.append({
            "name": name,
            "color": color_of(name),
            "ecotype": eco,
            "ecotype_size": len(members_of.get(eco, [])),
            "support": (support.get(str(eco)) or {}).get("sh_support", ""),
            "genomes": [{
                "accession": g,
                **{k: metadata.get(g, {}).get(k, "") for k in ("host", "niche", "country", "year")},
                **{k: sources.get(g, {}).get(k, "") for k in ("bioproject", "bioproject_title", "submitter")},
                "papers": [{k: papers[p][k] for k in ("pmid", "title", "journal", "year", "url", "summary")}
                           for p in sources.get(g, {}).get("pmids", "").split(";") if p in papers],
            } for g in genomes],
        })

    species = params.get("species") or folder.parent.name.replace("_", " ")
    source = "full core genome" if xml_path.name == FULL_CORE_XML else result["file"]
    n_genomes = genomes_of(t.name for t in tips)
    ecotypes = f"{len(result['ecotypes'])} ecotypes" if result["ecotypes"] else "no EcoSim ecotypes"
    title = f"{species} · {n_genomes} genomes · {ecotypes} ({source})"
    # EcoSim only knows the outgroup as "outgroup"; the run parameters hold its real name.
    root_name = params.get("outgroup")
    subtitle = f"Rooted with outgroup{f' {root_name}' if root_name else ''}; outgroup branch not drawn" if outgroup else ""
    data = {"title": title, "subtitle": subtitle, "gene_counts": gene_counts, "has_lengths": has_lengths,
            "nodes": nodes, "rows": rows, "tips": tip_info}
    page = (_TEMPLATE
            .replace("__TITLE__", html.escape(title))
            .replace("__SUBTITLE__", html.escape(subtitle))
            .replace("__DATA__", json.dumps(data, ensure_ascii=False).replace("</", "<\\/")))
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(page)
    return out_html


_TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>__TITLE__</title>
<style>
html,body{height:100%}
body{margin:0;font:13px system-ui,sans-serif;color:#222;background:#fcfcfb;display:flex;height:100vh}
#main{flex:1;min-width:0;display:flex;flex-direction:column}
header{padding:10px 14px 6px}
h1{font-size:15px;margin:0 0 3px}
#subtitle{margin:0 0 7px;color:#666;font-size:12px}
nav a{display:inline-block;padding:3px 10px;margin-right:4px;border:1px solid #ccc;border-radius:4px;color:#222;
  text-decoration:none;background:#fff}
nav a[aria-current=page]{background:#222;border-color:#222;color:#fff}
nav span{margin-left:8px;color:#666}
#tree{flex:1;min-height:0;overflow:auto}
svg{display:block}
svg text{font-size:11px;fill:#444}
.b{fill:none;stroke-width:1.3}
.band{fill:transparent;cursor:pointer}
.band:hover,.sel{fill:#00000012}
[data-z]{cursor:pointer}
[data-z]:hover .lab{text-decoration:underline}
[data-z]:focus,[data-i]:focus{outline:2px solid #0b5cad;outline-offset:-1px}
.lab{font-weight:700;fill:#222}
.axis{fill:#777}
.strain{font:11px ui-monospace,monospace;fill:#222;pointer-events:none}
.dim{opacity:.35}
#panel{width:340px;flex:none;border-left:1px solid #ddd;padding:12px 16px;overflow:auto;background:#fafafa}
#panel h2{font-size:15px;margin:0 0 8px}
.genome{border-top:1px solid #e3e3e3;padding:8px 0}
.paper{margin:6px 0 0 10px}
.muted{color:#777}
#panel a{color:#0b5cad}
</style></head><body>
<div id="main">
<header><h1>__TITLE__</h1><p id="subtitle">__SUBTITLE__</p>
<nav><a id="to-overview" href="#">Ecotypes</a><a id="to-all" href="#all">All strains</a><span id="where"></span></nav></header>
<div id="tree"></div>
</div>
<div id="panel"></div>
<script>
const DATA = __DATA__;
const {nodes: NODES, rows: ROWS, tips: TIPS, gene_counts: GENES} = DATA;
const box = document.getElementById("tree"), panel = document.getElementById("panel");
const esc = s => String(s ?? "").replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const kids = NODES.map(() => []);
NODES.forEach((n, i) => { if (n.p >= 0) kids[n.p].push(i); });
const MAX_X = Math.max(...NODES.map(n => n.x)) || 1;
const CLONES = TIPS.map(t => t.genomes.length - 1), MAX_CLONES = Math.max(0, ...CLONES);
const CELL = 38, HEAT_W = GENES.length ? GENES.length * CELL + 8 : 0;
const NAMES_W = 230, MIN_NAME_ROW = 11, DETAIL_TOP = 52, DETAIL_BOTTOM = 46;
let pinned = null, view = null, current = null, frame = null;

// One hue, light to dark: 0% is still a visible cell, so "never recovered" differs from "no data".
function heat(v) {
  const a = [214, 228, 251], b = [18, 48, 95];
  return `rgb(${a.map((c, i) => Math.round(c + (b[i] - c) * v)).join()})`;
}
const pct = v => `${Math.round(v * 100)}%`;
const genomesText = r => `${r.genomes} genome${r.genomes === 1 ? "" : "s"}` +
  (r.unique !== r.genomes ? ` · ${r.unique} unique` : "");
const rowLabel = r => r.eco === null ? r.name : `E${r.eco}`;
const target = r => r.eco === null ? "#all" : `#E${r.eco}`;

function heatHeader(x, y) {
  if (!GENES.length) return "";
  let h = `<text x="${x + HEAT_W / 2}" y="${y - 30}" text-anchor="middle">Recovered in rarefaction</text>` +
    `<text x="${x + HEAT_W / 2}" y="${y - 17}" text-anchor="middle">subsets (genes per subset)</text>`;
  GENES.forEach((g, c) => h += `<text class="axis" x="${x + 4 + (c + .5) * CELL}" y="${y - 3}" text-anchor="middle">${g}</text>`);
  return h;
}

function scaleBar(sx, x, y) {
  if (!DATA.has_lengths) return "";
  const raw = 80 / sx, pow = 10 ** Math.floor(Math.log10(raw));
  const v = [5, 2, 1].map(m => m * pow).find(m => m <= raw) ?? pow;
  return `<path d="M${x},${y}h${v * sx}" stroke="#222" stroke-width="1.3"/>` +
    `<text x="${x}" y="${y + 14}">${+v.toPrecision(1)} subst./site</text>`;
}

function overview() {
  const W = Math.max(box.clientWidth, 640), H = box.clientHeight;
  const top = 48, bottom = 40, GAP = 6, MIN_ROW = 21;
  const heightOf = (u, r) => Math.max(u * r.unique + GAP, MIN_ROW);
  let lo = 0, hi = 3;  // px per unique sequence, as large as fits (capped so few rows stay wedges)
  for (let k = 0; k < 24; k++) {
    const u = (lo + hi) / 2;
    if (ROWS.reduce((s, r) => s + heightOf(u, r), 0) <= H - top - bottom) lo = u; else hi = u;
  }
  const unit = lo;
  const heatX = W - HEAT_W - 14, labelX = heatX - 210, treeR = labelX - 16, treeL = 16;
  const sx = (treeR - treeL) / MAX_X, px = x => treeL + x * sx;
  let y = top, out = heatHeader(heatX, top);
  const rowMid = ROWS.map(r => { const m = y + heightOf(unit, r) / 2; y += heightOf(unit, r); return m; });
  const ny = NODES.map(() => 0);
  for (let i = NODES.length - 1; i >= 0; i--) {
    const n = NODES[i];
    if (n.h) continue;
    ny[i] = n.r >= 0 ? rowMid[n.r] : (ny[kids[i][0]] + ny[kids[i][kids[i].length - 1]]) / 2;
  }
  ROWS.forEach((r, ri) => {
    const m = rowMid[ri], rh = heightOf(unit, r);
    out += `<g data-z="${target(r)}" tabindex="0" role="link" aria-label="Zoom to ${esc(rowLabel(r))}: ${genomesText(r)}">` +
      `<rect class="band" x="0" y="${m - rh / 2}" width="${W}" height="${rh}"/>` +
      `<text class="lab" x="${labelX}" y="${m + 4}" style="font-size:13px">${esc(rowLabel(r))}</text>` +
      `<text x="${labelX + 38}" y="${m + 4}">${genomesText(r)}${r.split ? " · split clade" : ""}</text>`;
    (r.support || []).forEach((v, c) => {
      const cy = m - 8, cx = heatX + 4 + c * CELL;
      out += `<rect x="${cx + 1}" y="${cy}" width="${CELL - 2}" height="16" fill="${heat(v)}"><title>${esc(rowLabel(r))}: recovered in ${pct(v)} of ${GENES[c]}-gene subsets</title></rect>` +
        `<text x="${cx + CELL / 2}" y="${m + 4}" text-anchor="middle" style="fill:${v >= .6 ? "#fff" : "#223"};font-size:10px;pointer-events:none">${pct(v)}</text>`;
    });
    out += `</g>`;
  });
  NODES.forEach((n, i) => {
    if (n.h) return;
    const x = px(n.x);
    if (n.p >= 0) out += `<path class="b" pointer-events="none" d="M${px(NODES[n.p].x)},${ny[i]}H${x}" stroke="${n.r >= 0 ? n.c : "#444"}"/>`;
    if (n.r >= 0) {
      const r = ROWS[n.r], m = ny[i];
      if (r.unique === 1) out += `<circle pointer-events="none" cx="${x}" cy="${m}" r="3" fill="${r.color}"/>`;
      else {
        const half = Math.max(unit * r.unique, 5) / 2, end = Math.max(px(r.tip_x), x + 4);
        out += `<path pointer-events="none" d="M${x},${m}L${end},${m - half}V${m + half}Z" fill="${r.color}"/>`;
      }
    } else if (kids[i].length) {
      out += `<path class="b" pointer-events="none" d="M${x},${ny[kids[i][0]]}V${ny[kids[i][kids[i].length - 1]]}" stroke="#444"/>`;
    }
  });
  out += scaleBar(sx, treeL, y + 14);
  return `<svg width="${W}" height="${y + bottom}" role="img" aria-label="${esc(DATA.title)}">${out}</svg>`;
}

function detail(d) {
  // Fixed tracks need about 530px; below that the pane scrolls sideways rather than squeezing the tree away.
  const W = Math.max(box.clientWidth, 750 + d.names), top = DETAIL_TOP, bottom = DETAIL_BOTTOM;
  const span = d.y1 - d.y0, plotH = Math.max(box.clientHeight - top - bottom, 60);
  const row = Math.max(plotH / span, d.minRow), H = top + span * row + bottom;
  const heatX = W - HEAT_W - 14, cloneW = MAX_CLONES ? 92 : 0, cloneX = heatX - cloneW - 8;
  const labelX = cloneX - 150, stripX = labelX - 40, STRIP = 24;
  const treeL = 16, treeR = stripX - 8 - d.names;
  const sx = (treeR - treeL) / (d.x1 - d.x0), px = x => treeL + (x - d.x0) * sx, py = t => top + (t - d.y0) * row;
  const named = row >= MIN_NAME_ROW;
  const ny = NODES.map(() => 0);
  for (let i = NODES.length - 1; i >= 0; i--) {
    const n = NODES[i];
    ny[i] = n.t >= 0 ? n.t : (ny[kids[i][0]] + ny[kids[i][kids[i].length - 1]]) / 2;
  }
  const dim = t => d.focus && (t < d.focus[0] || t > d.focus[1]);
  let branches = "", names = "", bands = "", tracks = "";
  NODES.forEach((n, i) => {
    const x = px(n.x), y = py(ny[i]);
    if (n.p >= 0) branches += `<path class="b" d="M${px(NODES[n.p].x)},${y}H${x}" stroke="${n.c}"/>`;
    if (kids[i].length) branches += `<path class="b" d="M${x},${py(ny[kids[i][0]])}V${py(ny[kids[i][kids[i].length - 1]])}" stroke="${n.c}"/>`;
    if (n.t < 0) return;
    const t = TIPS[n.t], c = CLONES[n.t];
    if (py(n.t) < top - row || py(n.t) > H - bottom + row) return;
    branches += `<path d="M${x + 2},${y}H${stripX - 2}" stroke="#ececec" stroke-width=".7"/>`;
    bands += `<rect class="band${pinned === n.t ? " sel" : ""}" data-i="${n.t}" x="0" y="${y - row / 2}" width="${stripX - 2}" height="${row}"` +
      (named ? ` tabindex="0" role="button" aria-label="Show details of ${esc(t.name)}"` : "") + `/>`;
    if (named) names += `<text class="strain${dim(n.t) ? " dim" : ""}" x="${x + 5}" y="${y + 4}">${esc(t.name)}${c ? ` (+${c} clone${c > 1 ? "s" : ""})` : ""}</text>`;
    if (c) tracks += `<rect x="${cloneX}" y="${y - Math.min(row, 3) / 2}" width="${c / MAX_CLONES * (cloneW - 6)}" height="${Math.min(row, 3)}" fill="#888"><title>${esc(t.name)}: ${c} clone${c > 1 ? "s" : ""} collapsed</title></rect>`;
  });
  // Labels sit beside their strip, pushed apart where small ecotypes would overlap, with a leader line.
  let lastY = -Infinity;
  const labels = [];
  ROWS.forEach((r, ri) => {
    const a = py(r.first - .5), b = py(r.last + .5);
    if (b <= top || a >= H - bottom) return;
    const mid = (Math.max(a, top) + Math.min(b, H - bottom)) / 2;  // of the part in view
    const ly = Math.max(mid, lastY + 15);
    lastY = ly;
    labels.push({r, ri, a, b, mid, ly});
  });
  for (let k = labels.length - 1, floor = H - bottom; k >= 0; k--) {
    labels[k].ly = Math.min(labels[k].ly, floor);
    floor = labels[k].ly - 15;
  }
  labels.forEach(({r, a, b, mid, ly}) => {
    const faded = d.focus && (r.first < d.focus[0] || r.last > d.focus[1]);
    tracks += `<g data-z="${target(r)}" tabindex="0" role="link" aria-label="Zoom to ${esc(rowLabel(r))}: ${genomesText(r)}"${faded ? ' class="dim"' : ""}>` +
      `<rect x="${stripX}" y="${a + .5}" width="${STRIP}" height="${Math.max(b - a - 1, 1)}" fill="${r.color}"/>` +
      `<path d="M${stripX + STRIP + 1},${mid}L${labelX - 3},${ly - 4}" stroke="#999" stroke-width=".7"/>` +
      `<text class="lab" x="${labelX}" y="${ly}">${esc(rowLabel(r))}</text>` +
      `<text x="${labelX + 32}" y="${ly}">${r.genomes} genome${r.genomes === 1 ? "" : "s"}</text>`;
    (r.support || []).forEach((v, c) => {
      tracks += `<rect x="${heatX + 4 + c * CELL + 1}" y="${a + .5}" width="${CELL - 2}" height="${Math.max(b - a - 1, 1)}" fill="${heat(v)}"><title>${esc(rowLabel(r))}: recovered in ${pct(v)} of ${GENES[c]}-gene subsets</title></rect>`;
    });
    tracks += `</g>`;
  });
  let chrome = heatHeader(heatX, top - 4);
  if (MAX_CLONES) chrome += `<text x="${cloneX + cloneW / 2 - 3}" y="${top - 21}" text-anchor="middle">Clones collapsed</text>` +
    `<text x="${cloneX + cloneW / 2 - 3}" y="${top - 7}" text-anchor="middle">into each tip</text>` +
    `<text class="axis" x="${cloneX}" y="${H - bottom + 18}" text-anchor="middle">0</text>` +
    `<text class="axis" x="${cloneX + cloneW - 6}" y="${H - bottom + 18}" text-anchor="middle">${MAX_CLONES}</text>`;
  if (GENES.length) {
    const lx = heatX + 4, lw = GENES.length * CELL;
    for (let k = 0; k < lw; k += 2) chrome += `<rect x="${lx + k}" y="${H - bottom + 8}" width="2.5" height="8" fill="${heat(k / lw)}"/>`;
    [0, .5, 1].forEach(v => chrome += `<text class="axis" x="${lx + v * lw}" y="${H - bottom + 28}" text-anchor="middle">${pct(v)}</text>`);
  }
  return `<svg width="${W}" height="${H}" role="img" aria-label="${esc(DATA.title)}">` +
    `<defs><clipPath id="pane"><rect x="0" y="${top}" width="${stripX - 2}" height="${H - top - bottom}"/></clipPath></defs>` +
    `<clipPath id="rows"><rect x="0" y="${top}" width="${W}" height="${H - top - bottom}"/></clipPath>` +
    `<g clip-path="url(#pane)"><g pointer-events="none">${branches}</g>${bands}${names}</g>` +
    `<g clip-path="url(#rows)">${tracks}</g>${chrome}` +
    scaleBar(sx, treeL, H - bottom + 14) + `</svg>`;
}

function restoringFocus(html) {
  const f = document.activeElement, keep = f && (f.dataset?.z || f.dataset?.i);
  box.innerHTML = html;
  if (keep != null) box.querySelector(`[data-z="${f.dataset.z}"],[data-i="${f.dataset.i}"]`)?.focus();
}
const draw = () => restoringFocus(view === "overview" ? overview() : detail(current));

function domainOf(eco) {
  const rows = ROWS.filter(r => r.eco === eco);
  const first = Math.min(...rows.map(r => r.first)), last = Math.max(...rows.map(r => r.last));
  const crowns = NODES.filter(n => n.r >= 0 && ROWS[n.r].eco === eco);
  const left = Math.min(...crowns.map(n => NODES[n.p]?.x ?? n.x));
  const right = Math.max(...rows.map(r => r.tip_x));
  const pad = Math.max((right - left) * .04, MAX_X * .002);
  // At most 22px a row: a small ecotype shows its neighbours (faded) for context, not giant rows.
  const fit = Math.max(Math.ceil((box.clientHeight - DETAIL_TOP - DETAIL_BOTTOM) / 22), last - first + 3);
  const mid = (first + last) / 2;
  let y0 = mid - fit / 2, y1 = mid + fit / 2;
  if (y0 < -.5) { y1 += -.5 - y0; y0 = -.5; }
  if (y1 > TIPS.length - .5) { y0 = Math.max(-.5, y0 - (y1 - TIPS.length + .5)); y1 = TIPS.length - .5; }
  return {x0: left - pad, x1: right + pad, y0, y1, names: NAMES_W, minRow: 14, focus: [first, last]};
}
const allDomain = () => ({x0: 0, x1: MAX_X * 1.01, y0: -.5, y1: TIPS.length - .5, names: 0, minRow: 1.5, focus: null});

function animateTo(to) {
  clearTimeout(frame);
  const from = current || allDomain(), still = matchMedia("(prefers-reduced-motion: reduce)").matches;
  const start = performance.now(), T = 320;
  (function step() {  // a timer, not requestAnimationFrame: the last frame still lands in a throttled tab
    const k = still ? 1 : Math.min((performance.now() - start) / T, 1), e = k * k * (3 - 2 * k);
    current = {...to};
    if (k < 1) {
      for (const p of ["x0", "x1", "y0", "y1", "names", "minRow"]) current[p] = from[p] + (to[p] - from[p]) * e;
      frame = setTimeout(step, 16);
    }
    draw();
  })();
}

function route() {
  const hash = location.hash.slice(1), eco = /^E([0-9]+)$/.exec(hash), n = eco && +eco[1];
  const known = eco && ROWS.some(r => r.eco === n);
  // Without ecotypes the overview is one row per strain, so open on the full tree instead.
  const next = known ? n : hash === "all" || !ROWS.some(r => r.eco !== null) ? "all" : "overview";
  const was = view;
  view = next;
  document.getElementById("to-overview").toggleAttribute("aria-current", false);
  document.getElementById("to-all").toggleAttribute("aria-current", false);
  if (next === "overview") document.getElementById("to-overview").setAttribute("aria-current", "page");
  if (next === "all") document.getElementById("to-all").setAttribute("aria-current", "page");
  document.getElementById("where").textContent = known ? `Zoomed to ecotype E${n}` : "";
  clearTimeout(frame);
  if (next === "overview") { current = null; draw(); box.scrollTop = 0; return showHelp(); }
  if (was === "overview" || was === null) current = allDomain();
  if (next === "all") animateTo(allDomain()); else animateTo(domainOf(n));
  box.scrollTop = 0;
  if (known) showEcotype(n);
  else if (pinned === null) showHelp();
}

function showHelp() {
  if (pinned !== null) return show(pinned);
  panel.innerHTML = `<p class="muted">Click an ecotype to zoom in and see every strain.</p>
    <p class="muted">Hover or click a strain for its ecotype, clones, metadata and source papers. Click keeps the selection.</p>`;
}
function showEcotype(eco) {
  if (pinned !== null && TIPS[pinned].ecotype === eco) return show(pinned);
  pinned = null;
  const rows = ROWS.filter(r => r.eco === eco), r = rows[0];
  let h = `<h2>Ecotype ${eco}</h2><p>${genomesText(r)}` + (r.sh ? ` &middot; clade support ${esc(r.sh)}` : "") +
    (r.split ? ` &middot; not a single clade in this tree` : "") + `</p>`;
  if (r.support) h += `<p>Recovered exactly in: ` + r.support.map((v, c) => `${pct(v)} of ${GENES[c]}-gene subsets`).join(", ") + `.</p>`;
  panel.innerHTML = h + `<p class="muted">Hover or click a strain for its clones, metadata and source papers.</p>`;
}
function show(i) {
  const t = TIPS[i];
  let h = `<h2>${esc(t.name)}</h2>`;
  if (t.ecotype) h += `<p><b>Ecotype ${t.ecotype}</b> &middot; ${t.ecotype_size} unique sequence(s)` +
    (t.support ? ` &middot; clade support ${esc(t.support)}` : "") + `</p>`;
  h += `<p class="muted">${t.genomes.length} genome(s) in this tip</p>`;
  for (const g of t.genomes) {
    h += `<div class="genome"><b><a href="https://www.ncbi.nlm.nih.gov/datasets/genome/${esc(g.accession)}/" target="_blank">${esc(g.accession)}</a></b>`;
    const meta = ["host", "niche", "country", "year"].filter(k => g[k] && g[k] !== "unknown").map(k => `${k}: ${esc(g[k])}`);
    if (meta.length) h += `<br>${meta.join(" &middot; ")}`;
    if (g.bioproject) h += `<br>BioProject <a href="https://www.ncbi.nlm.nih.gov/bioproject/${esc(g.bioproject)}" target="_blank">${esc(g.bioproject)}</a> ${esc(g.bioproject_title)}`;
    if (g.submitter) h += `<br><span class="muted">Submitted by ${esc(g.submitter)}</span>`;
    for (const p of g.papers) h += `<div class="paper"><a href="${esc(p.url)}" target="_blank">${esc(p.title)}</a>
      <span class="muted">${esc(p.journal)} ${esc(p.year)}</span>${p.summary ? `<br>${esc(p.summary)}` : ""}</div>`;
    if (g.bioproject && !g.papers.length) h += `<div class="paper muted">No paper linked to this BioProject in NCBI.</div>`;
    h += `</div>`;
  }
  panel.innerHTML = h;
}

function activate(el) {
  if (el.dataset.z) {
    if (location.hash === el.dataset.z) route();
    else location.hash = el.dataset.z;
    return;
  }
  const i = +el.dataset.i;
  pinned = pinned === i ? null : i;
  draw();
  if (pinned !== null) show(pinned);
}
box.addEventListener("click", e => { const el = e.target.closest("[data-z],[data-i]"); if (el) activate(el); });
box.addEventListener("mouseover", e => {
  const el = e.target.closest("[data-i]");
  if (el && pinned === null) show(+el.dataset.i);
});
document.addEventListener("keydown", e => {
  const el = e.target.closest?.("[data-z],[data-i]");
  if (el && (e.key === "Enter" || e.key === " ")) { e.preventDefault(); activate(el); }
  if (e.key === "Escape" && location.hash) location.hash = "";
});
addEventListener("hashchange", route);
addEventListener("resize", () => {
  if (typeof view === "number") current = domainOf(view); else if (view === "all") current = allDomain();
  draw();
});
route();
</script></body></html>
"""


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ecosim_out_dir")
    parser.add_argument("--output", default=None, help="HTML path (default: ecotype_plots/<xml>_interactive.html)")
    args = parser.parse_args()
    path = write_interactive_tree(args.ecosim_out_dir, args.output)
    if path:
        print(f"Saved {path}")
