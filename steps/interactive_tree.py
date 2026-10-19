"""
Interactive ecotype tree: one self-contained HTML file, no server or internet needed.

Draws the full core genome EcoSim tree (or the replicate with the most genes) with
branches colored by ecotype. Hovering or clicking a tip shows its ecotype, the clones
collapsed into it, and whatever metadata earlier steps saved next to the XML:
  genome_metadata.csv  host / niche / country / year      (--niche)
  genome_sources.csv   BioProject, submitter, PubMed ids (--papers)
  papers.csv           paper title, link and Claude summary (--papers)
  ecotype_tree_support.csv  SH-like support of each ecotype's clade

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
    from steps.visualize_ecotypes import _ecotype_colors, OUTGROUP_COLOR, SINGLETON_COLOR
except ImportError:  # run as a script: python steps/interactive_tree.py
    from parsing import parse_ecosim_xml, load_clone_map, FULL_CORE_XML
    from visualize_ecotypes import _ecotype_colors, OUTGROUP_COLOR, SINGLETON_COLOR

ROW = 14       # px per tip
WIDTH = 760    # px for the deepest tip
BRANCH_COLOR = "#9a9a9a"


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


def write_interactive_tree(ecosim_out_dir, out_html=None):
    """Write the HTML viewer. Returns its path, or None when there is no EcoSim XML with a tree."""
    folder = Path(ecosim_out_dir)
    xml_path = _pick_xml(folder)
    result = parse_ecosim_xml(xml_path) if xml_path else None
    if not result or not result["tree"]:
        print(f"No EcoSim XML with a tree in {folder}; skipping the interactive tree.")
        return None
    out_html = Path(out_html or folder / "ecotype_plots" / f"{xml_path.stem}_interactive.html")

    tree = Phylo.read(StringIO(result["tree"]), "newick")
    clones = load_clone_map(folder / "clone_groups.csv")
    metadata = _rows(folder / "genome_metadata.csv", "accession")
    sources = _rows(folder / "genome_sources.csv", "accession")
    papers = _rows(folder / "papers.csv", "pmid")
    support = _rows(folder / "ecotype_tree_support.csv", "ecotype")
    ecotype_of = {m: eco["number"] for eco in result["ecotypes"] for m in eco["members"]}
    palette = _ecotype_colors(result)
    outgroup = result["outgroup"]

    def color_of(name):
        if name == outgroup:
            return OUTGROUP_COLOR
        return palette.get(ecotype_of.get(name), SINGLETON_COLOR)

    # Layout: x = distance from the root, y = tip order; an internal node sits midway between its children.
    depths = tree.depths()
    if not any(depths.values()):  # no branch lengths: fall back to node depth
        depths = tree.depths(unit_branch_lengths=True)
    scale = WIDTH / (max(depths.values()) or 1)
    y, tips = {}, tree.get_terminals()
    for i, tip in enumerate(tips):
        y[tip] = 10 + i * ROW
    for clade in tree.get_nonterminals(order="postorder"):
        y[clade] = (y[clade.clades[0]] + y[clade.clades[-1]]) / 2

    # A branch takes its ecotype's color when every tip below it is in that one (multi-member) ecotype.
    branch_color = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            branch_color[clade] = color_of(clade.name) if clade.name == outgroup or ecotype_of.get(clade.name) in palette else BRANCH_COLOR
        else:
            below = {branch_color[c] for c in clade.clades}
            branch_color[clade] = below.pop() if len(below) == 1 else BRANCH_COLOR

    lines, labels, tip_info = [], [], []
    for clade in tree.find_clades():
        x = 10 + depths[clade] * scale
        for child in clade.clades:
            cx = 10 + depths[child] * scale
            lines.append(f'<path d="M{x:.1f},{y[child]:.1f}H{cx:.1f}" stroke="{branch_color[child]}"/>')
        if clade.clades:
            lines.append(f'<path d="M{x:.1f},{y[clade.clades[0]]:.1f}V{y[clade.clades[-1]]:.1f}" '
                         f'stroke="{branch_color[clade]}"/>')
    for i, tip in enumerate(tips):
        name = tip.name
        genomes = clones.get(name, [name])
        eco = ecotype_of.get(name)
        label = name + (" [outgroup]" if name == outgroup else f" [E{eco}]" if eco else "") \
            + (f" (+{len(genomes) - 1} clones)" if len(genomes) > 1 else "")
        x = 10 + depths[tip] * scale
        labels.append(f'<text x="{x + 5:.1f}" y="{y[tip] + 4:.1f}" fill="{color_of(name)}" data-i="{i}">'
                      f'{html.escape(label)}</text>')
        tip_info.append({
            "name": name,
            "ecotype": eco,
            "ecotype_size": len(next((e["members"] for e in result["ecotypes"] if e["number"] == eco), [])),
            "support": (support.get(str(eco)) or {}).get("sh_support", ""),
            "outgroup": name == outgroup,
            "genomes": [{
                "accession": g,
                **{k: metadata.get(g, {}).get(k, "") for k in ("host", "niche", "country", "year")},
                **{k: sources.get(g, {}).get(k, "") for k in ("bioproject", "bioproject_title", "submitter")},
                "papers": [{k: papers[p][k] for k in ("pmid", "title", "journal", "year", "url", "summary")}
                           for p in sources.get(g, {}).get("pmids", "").split(";") if p in papers],
            } for g in genomes],
        })

    height = 20 + len(tips) * ROW
    width = WIDTH + 420
    title = f"{result['file']}: {len(result['ecotypes'])} ecotypes, {len(tips)} tips"
    page = (_TEMPLATE
            .replace("__TITLE__", html.escape(title))
            .replace("__SVG__", f'<svg width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
                                f'<g class="br">{"".join(lines)}</g><g class="tip">{"".join(labels)}</g></svg>')
            .replace("__DATA__", json.dumps(tip_info).replace("</", "<\\/")))
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(page)
    return out_html


_TEMPLATE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>__TITLE__</title>
<style>
body{margin:0;font:13px system-ui,sans-serif;color:#222;background:#fff;display:flex;height:100vh}
#tree{flex:1;overflow:auto;padding:8px}
#panel{width:380px;border-left:1px solid #ddd;padding:12px 16px;overflow:auto;background:#fafafa}
h1{font-size:15px;margin:0 0 8px}
.br path{fill:none;stroke-width:1.4}
.tip text{font:11px ui-monospace,monospace;cursor:pointer}
.tip text:hover,.tip text.sel{font-weight:bold;text-decoration:underline}
.genome{border-top:1px solid #e3e3e3;padding:8px 0}
.paper{margin:6px 0 0 10px}
.muted{color:#777}
a{color:#0b5cad}
</style></head><body>
<div id="tree"><h1>__TITLE__</h1>__SVG__</div>
<div id="panel"><p class="muted">Hover or click a tip to see its ecotype, clones, metadata and source papers.
Click keeps the selection.</p></div>
<script>
const TIPS = __DATA__;
const panel = document.getElementById("panel");
let pinned = null;
const esc = s => String(s ?? "").replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
function show(i) {
  const t = TIPS[i];
  let h = `<h1>${esc(t.name)}</h1>`;
  if (t.outgroup) h += `<p>Outgroup (roots the tree; not part of any ecotype).</p>`;
  else if (t.ecotype) h += `<p><b>Ecotype ${t.ecotype}</b> &middot; ${t.ecotype_size} unique sequence(s)` +
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
document.querySelectorAll(".tip text").forEach(el => {
  el.addEventListener("mouseenter", () => { if (pinned === null) show(el.dataset.i); });
  el.addEventListener("click", () => {
    document.querySelectorAll(".tip text.sel").forEach(s => s.classList.remove("sel"));
    pinned = pinned === el.dataset.i ? null : el.dataset.i;
    if (pinned !== null) { el.classList.add("sel"); show(pinned); }
  });
});
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
