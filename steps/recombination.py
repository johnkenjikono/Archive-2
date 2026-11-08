"""
Step 4d (optional, --mask-recombination): mask recombinant regions with Gubbins.

Gubbins (Croucher et al. 2015) walks the tree branch by branch and flags stretches whose SNP
density says they were imported by recombination rather than built up by mutation. EcoSim
models divergence by mutation alone, so masking those stretches leaves the clonal frame it
assumes.

Gubbins finds those stretches from where SNPs sit along the genome, and its manual warns
that Roary/Panaroo core-gene concatenations do not work as input: Panaroo concatenates
genes in file-name order, which scatters genome neighbours. So the core genes are first put
back in the order they sit in one reference genome (Bakta numbers locus tags along the
chromosome, and Panaroo's gene_presence_absence.csv lists them), Gubbins runs on that, and
its blocks are mapped back onto Panaroo's gene order, which steps 5-7 and the gene header use.

Only the ingroup goes to Gubbins. The outgroup is another species, and long branches are
where Gubbins is most prone to false positives (its manual's advice); it is left unmasked.

Masked sites become '-'. EcoSim never compares sequences: it bins on the tree's branch
lengths and takes only the alignment length and outgroup name from the FASTA. VeryFastTree
treats gaps as missing data, so a masked site counts as unknown, not as a difference.

Writes into out_dir:
  recombination_predictions.gff   Gubbins blocks (coordinates in reference gene order)
  recombination_per_branch.tsv    Gubbins per-branch statistics (r/m etc.; Gubbins names it .csv but writes tabs)
  recombination_masking.csv       per genome: blocks, masked sites, masked fraction

  python steps/recombination.py sorted.fasta core_alignment_header.embl gene_presence_absence.csv out_dir [outgroup_id]
"""
import bisect
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from Bio import Phylo
from Bio.SeqIO.FastaIO import SimpleFastaParser

try:
    from steps.rarefaction import load_gene_coordinates
    from steps.reroot_tree import _resolve_outgroup_name
except ImportError:  # run as a script: python steps/recombination.py
    from rarefaction import load_gene_coordinates
    from reroot_tree import _resolve_outgroup_name

# Gubbins pins Python 3.10 and numpy 1.23.0, so it lives in its own conda env
# (environment-gubbins.yml) rather than in the ecotype env.
GUBBINS_ENV = "gubbins"
PREFIX = "gubbins"
# Bakta locus tags (PREFIX_01234) count up along each contig; Panaroo's re-found genes
# ("12_refound_34") carry no position.
_LOCUS_RE = re.compile(r"_(\d+)(?:_pseudo)?$")
_TAXA_RE = re.compile(r'taxa="([^"]*)"')
_NOT_GUBBINS_BASE = re.compile(r"[^ACGTN-]")


def gubbins_command():
    """run_gubbins.py from the `gubbins` conda env when it exists, else from PATH; None if neither.
    The env wins, as in pipeline._probe_version, so the version recorded is the one that ran."""
    conda = shutil.which("conda")
    if conda:
        try:
            listing = subprocess.run([conda, "env", "list", "--json"], capture_output=True, text=True, timeout=120)
            envs = json.loads(listing.stdout).get("envs", [])
        except (OSError, subprocess.SubprocessError, ValueError):
            envs = []
        if any(Path(env).name == GUBBINS_ENV for env in envs):
            return [conda, "run", "--no-capture-output", "-n", GUBBINS_ENV, "run_gubbins.py"]
    exe = shutil.which("run_gubbins.py")
    return [exe] if exe else None


def find_presence_absence(alignment_path, base_dir):
    """Panaroo's gene_presence_absence.csv: next to the alignment, else the saved pangenome/ copy."""
    for candidate in (Path(alignment_path).resolve().parent / "gene_presence_absence.csv",
                      Path(base_dir) / "pangenome" / "gene_presence_absence.csv"):
        if candidate.exists():
            return candidate
    return None


def _locus_number(cell):
    """Position of a gene in one genome from its locus tag(s), or None (absent, re-found, unparsable)."""
    nums = [int(m.group(1)) for tag in cell.split(";")
            if tag and "refound" not in tag and (m := _LOCUS_RE.search(tag))]
    return min(nums) if nums else None


def genome_order(gene_names, presence_csv, genomes):
    """
    Order core genes as they sit in one reference genome.

    The reference is the genome in `genomes` placing the most of gene_names, then the one with
    the fewest multi-copy genes (placed by their first copy), then the first by name. Returns
    (reference, ordered gene_names, n_unplaced); genes it lacks keep their relative order at the end.
    """
    with open(presence_csv, newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        wanted = set(gene_names)
        rows = {row[0]: row for row in reader if row and row[0] in wanted}
    missing = [g for g in gene_names if g not in rows]
    if missing:
        raise ValueError(f"{len(missing)} core gene(s) (e.g. {missing[0]}) are not in {presence_csv}; "
                         "it must come from the same Panaroo run as the alignment.")
    columns = {name: i for i, name in enumerate(header) if name in set(genomes)}
    if not columns:
        raise ValueError(f"None of the alignment's genomes are columns of {presence_csv}.")
    cells = {genome: {g: rows[g][i] if i < len(rows[g]) else "" for g in gene_names} for genome, i in columns.items()}
    placed = {genome: {g: n for g, cell in genome_cells.items() if (n := _locus_number(cell)) is not None}
              for genome, genome_cells in cells.items()}
    multi_copy = {genome: sum(";" in cell for cell in genome_cells.values()) for genome, genome_cells in cells.items()}
    reference = min(placed, key=lambda genome: (-len(placed[genome]), multi_copy[genome], genome))
    keys = placed[reference]
    ordered = sorted(keys, key=keys.get) + [g for g in gene_names if g not in keys]
    return reference, ordered, len(gene_names) - len(keys)


def _segments(genes, order):
    """[(start in reordered alignment, start in original alignment, length)] for genes taken in order."""
    by_name = {name: (start, end) for name, start, end in genes}
    segments, position = [], 0
    for name in order:
        start, end = by_name[name]
        segments.append((position, start, end - start))
        position += end - start
    return segments


def _to_original(block_start, block_end, segments, segment_starts):
    """Original-alignment intervals [a, b) covered by reordered interval [block_start, block_end)."""
    out = []
    i = max(0, bisect.bisect_right(segment_starts, block_start) - 1)
    while i < len(segments) and segments[i][0] < block_end:
        new_start, orig_start, length = segments[i]
        lo, hi = max(block_start, new_start), min(block_end, new_start + length)
        if lo < hi:
            out.append((orig_start + lo - new_start, orig_start + hi - new_start))
        i += 1
    return out


def _merge(intervals):
    merged = []
    for a, b in sorted(intervals):
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return merged


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_gubbins(alignment, work_dir, threads, seed, command=None):
    """
    Run Gubbins on alignment inside work_dir (it writes to its working directory).

    Outputs from an earlier run on an identical alignment are reused, so a restarted run
    does not repeat it. Returns (gff, per-branch csv, final tree) paths; raises on failure.
    """
    work_dir = Path(work_dir)
    gff = work_dir / f"{PREFIX}.recombination_predictions.gff"
    stats = work_dir / f"{PREFIX}.per_branch_statistics.csv"
    tree = work_dir / f"{PREFIX}.final_tree.tre"
    stamp = work_dir / f"{PREFIX}.input.sha256"
    checksum = _sha256(alignment)
    if all(p.exists() for p in (gff, stats, tree)) and stamp.exists() and stamp.read_text().strip() == checksum:
        print(f"Reusing Gubbins results for this alignment in {work_dir}")
        return gff, stats, tree
    command = command or gubbins_command()
    if not command:
        raise RuntimeError("run_gubbins.py not found: create the gubbins env "
                           "(conda env create -f environment-gubbins.yml) or put Gubbins on PATH.")
    for old in work_dir.glob(f"{PREFIX}.*"):
        old.unlink()
    cmd = [*command, str(Path(alignment).resolve()),
           "--prefix", PREFIX,
           "--threads", str(threads),
           # Same tree builder as step 5; Gubbins fits the reconstruction model with IQ-TREE.
           "--tree-builder", "veryfasttree",
           # Gubbins drops genomes with more than 25% gaps by default, which would silently leave
           # them unmasked; keep every genome and check below that all made it into the tree.
           "--filter-percentage", "100",
           "--seed", str(seed)]
    print(f"Command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=work_dir)
    # Gubbins exits 0 on some errors (e.g. an invalid model), so check for its outputs.
    missing = [p.name for p in (gff, stats, tree) if not p.exists()]
    if missing:
        raise RuntimeError(f"Gubbins finished without writing {', '.join(missing)}; see its log above.")
    stamp.write_text(checksum + "\n")
    return gff, stats, tree


def read_blocks(gff):
    """[(start, end, taxa)] from a Gubbins GFF: 1-based inclusive coordinates."""
    blocks = []
    with open(gff) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            taxa = _TAXA_RE.search(fields[8])
            blocks.append((int(fields[3]), int(fields[4]), taxa.group(1).split() if taxa else []))
    return blocks


def mask_recombinant_regions(fasta_path, gene_header, presence_csv, out_dir, work_dir,
                             outgroup_id=None, threads=8, seed=42, command=None):
    """
    Mask Gubbins-detected recombination in fasta_path, rewriting it in place.

    Returns {"reference", "blocks", "masked_sites", "genes_emptied"} or None on failure.
    """
    out_dir, work_dir = Path(out_dir), Path(work_dir)
    with open(fasta_path) as handle:
        records = [(title.split(None, 1)[0], seq) for title, seq in SimpleFastaParser(handle)]
    ids = [seq_id for seq_id, _ in records]
    root = (_resolve_outgroup_name(outgroup_id, ids) if outgroup_id else None) or ids[0]
    ingroup = [(seq_id, seq) for seq_id, seq in records if seq_id != root]
    if len(ingroup) < 3:
        print(f"Only {len(ingroup)} ingroup genome(s); Gubbins needs 3 or more, so nothing was masked.")
        return {"reference": None, "blocks": 0, "masked_sites": 0, "genes_emptied": 0}

    length = len(records[0][1])
    genes = load_gene_coordinates(gene_header)
    if not genes or genes[-1][2] != length:
        print(f"❌ {gene_header} does not describe this {length} bp alignment.")
        return None
    # Panaroo labels each block with its gene file name (abrB.aln); the gene table uses abrB.
    names = [name[:-4] if name.endswith(".aln") else name for name, _, _ in genes]
    genes = [(name, start, end) for name, (_, start, end) in zip(names, genes)]
    try:
        reference, order, unplaced = genome_order(names, presence_csv, [seq_id for seq_id, _ in ingroup])
    except (OSError, ValueError) as exc:
        print(f"❌ Could not put the core genes in genome order: {exc}")
        return None
    print(f"Core genes put in genome order using {reference} ({len(names) - unplaced} of {len(names)} genes placed"
          + (f"; {unplaced} it lacks go at the end" if unplaced else "") + ")")
    segments = _segments(genes, order)
    segment_starts = [s[0] for s in segments]

    work_dir.mkdir(parents=True, exist_ok=True)
    gubbins_input = work_dir / "ingroup_genome_order.fasta"
    by_name = {name: (start, end) for name, start, end in genes}
    with open(gubbins_input, "w") as out_f:
        for seq_id, seq in ingroup:
            reordered = "".join(seq[by_name[g][0]:by_name[g][1]] for g in order)
            # Gubbins accepts only ACGTN and '-'; other IUPAC codes become N.
            out_f.write(f">{seq_id}\n{_NOT_GUBBINS_BASE.sub('N', reordered.upper())}\n")

    try:
        gff, stats, tree = run_gubbins(gubbins_input, work_dir, threads, seed, command)
    except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"❌ Gubbins failed: {exc}")
        return None
    leaves = {leaf.name for leaf in Phylo.read(str(tree), "newick").get_terminals()}
    dropped = sorted({seq_id for seq_id, _ in ingroup} - leaves)
    if dropped:
        print(f"❌ Gubbins left {len(dropped)} genome(s) out of its tree ({', '.join(dropped[:5])}"
              f"{'...' if len(dropped) > 5 else ''}), so they would go unmasked.")
        return None

    blocks = read_blocks(gff)
    intervals = {}  # genome -> original-alignment intervals
    block_count = {}
    for start, end, taxa in blocks:
        spans = _to_original(start - 1, end, segments, segment_starts)
        for taxon in taxa:
            intervals.setdefault(taxon, []).extend(spans)
            block_count[taxon] = block_count.get(taxon, 0) + 1
    unknown = sorted(set(intervals) - leaves)
    if unknown:
        print(f"❌ Gubbins reported blocks in genomes not in the alignment: {', '.join(unknown[:5])}")
        return None

    masked, any_column, emptied = {}, [], set()
    sequences = dict(records)
    for taxon, spans in intervals.items():
        merged = _merge(spans)
        seq = bytearray(sequences[taxon].encode())
        before = bytes(seq)
        for a, b in merged:
            seq[a:b] = b"-" * (b - a)
        # Genes now all gaps in this genome drop out of step 7's sampling pool.
        for name, start, end in genes:
            if not seq[start:end].strip(b"-") and before[start:end].strip(b"-"):
                emptied.add(name)
        sequences[taxon] = seq.decode()
        masked[taxon] = sum(b - a for a, b in merged)
        any_column.extend(merged)
    columns_masked = sum(b - a for a, b in _merge(any_column))

    tmp_out = f"{fasta_path}.tmp"
    with open(tmp_out, "w") as out_f:
        for seq_id, _ in records:
            out_f.write(f">{seq_id}\n{sequences[seq_id]}\n")
    os.replace(tmp_out, fasta_path)

    out_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(gff, out_dir / "recombination_predictions.gff")
    shutil.copy2(stats, out_dir / "recombination_per_branch.tsv")
    summary_csv = out_dir / "recombination_masking.csv"
    with open(summary_csv, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["genome", "blocks", "masked_sites", "masked_fraction"])
        for seq_id, _ in ingroup:
            writer.writerow([seq_id, block_count.get(seq_id, 0), masked.get(seq_id, 0),
                             f"{masked.get(seq_id, 0) / length:.4g}"])
        # The "all" row counts every Gubbins block, and alignment columns masked in at least one genome.
        writer.writerow(["all", len(blocks), columns_masked, f"{columns_masked / length:.4g}"])

    total = sum(masked.values())
    print(f"✅ Gubbins found {len(blocks)} recombination block(s); masked {total:,} sites across "
          f"{len(masked)} of {len(ingroup)} ingroup genome(s) ({total / (length * len(ingroup)):.2%} of ingroup sites; "
          f"{columns_masked / length:.1%} of columns in at least one genome). Root {root} left unmasked.")
    if emptied:
        print(f"   {len(emptied)} core gene(s) are now entirely masked in at least one genome; "
              "step 7 leaves them out of its sampling pool.")
    print(f"   Saved {summary_csv.name}, recombination_predictions.gff and recombination_per_branch.tsv in {out_dir}")
    return {"reference": reference, "blocks": len(blocks), "masked_sites": total, "genes_emptied": len(emptied)}


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Mask Gubbins-detected recombination in a Panaroo core alignment (in place)")
    parser.add_argument("fasta", help="Sorted core alignment (rewritten in place)")
    parser.add_argument("gene_header", help="core_alignment_header.embl from the same Panaroo run")
    parser.add_argument("presence_csv", help="gene_presence_absence.csv from the same Panaroo run")
    parser.add_argument("out_dir", help="Where the GFF and summary CSVs go")
    parser.add_argument("outgroup_id", nargs="?", default=None)
    parser.add_argument("--work-dir", default=None, help="Gubbins working folder (default: <out_dir>/gubbins_work)")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    if mask_recombinant_regions(args.fasta, args.gene_header, args.presence_csv, args.out_dir,
                                args.work_dir or Path(args.out_dir) / "gubbins_work",
                                args.outgroup_id, args.threads) is None:
        sys.exit(1)
