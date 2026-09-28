import csv
import os
import random
import re
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Any feature key: Panaroo writes "feature"; Roary also writes "misc_feature" for group_* genes.
_FEATURE_RE = re.compile(r"^FT   (\S+)\s+(\d+)\.\.(\d+)")
_LABEL_RE = re.compile(r"^FT\s+/label=(.+)$")
# The whole core alignment as one extra replicate, so every gene subset can be
# compared against the complete core genome.
FULL_CORE_NAME = "full_core_genome"


def load_gene_coordinates(header_embl):
    """
    Parse a Panaroo/Roary core_alignment_header.embl into [(label, start, end), ...].

    Coordinates are 0-based, end-exclusive slices into the concatenated
    core_gene_alignment.aln, one entry per core gene, in alignment order.
    """
    genes = []
    current = None
    with open(header_embl) as handle:
        for line in handle:
            feature = _FEATURE_RE.match(line)
            if feature:
                if current:
                    genes.append(tuple(current))
                current = [f"gene_{len(genes) + 1}", int(feature.group(2)) - 1, int(feature.group(3))]
                continue
            label = _LABEL_RE.match(line)
            if label and current:
                current[0] = label.group(1).strip()
    if current:
        genes.append(tuple(current))
    return genes


def find_gene_header(alignment_path):
    """
    Locate the gene-coordinate header Panaroo wrote next to a core alignment:
    core_gene_alignment.aln -> core_alignment_header.embl,
    core_gene_alignment_filtered.aln -> core_alignment_filtered_header.embl.
    """
    folder = os.path.dirname(os.path.abspath(alignment_path))
    name = "core_alignment_filtered_header.embl" if "filtered" in os.path.basename(alignment_path) \
        else "core_alignment_header.embl"
    candidate = os.path.join(folder, name)
    return candidate if os.path.exists(candidate) else None


def create_rarefaction_fastas(input_fasta,
                              output_folder="rarefaction_fastas",
                              gene_header=None,
                              gene_counts=None,
                              trials_per_count=20,
                              random_seed=42,
                              manifest_path=None,
                              skip_absent_genes=True,
                              each_gene=False,
                              include_full=True):
    """
    Build rarefaction replicates by sampling whole core genes, not fixed-width
    windows: each sim_species_g{N}_t{T}.fasta concatenates N distinct genes
    (drawn without replacement) using the boundaries in gene_header.

    Panaroo pads a core gene with '-' in genomes that lack it. With
    skip_absent_genes, genes that are all gaps in any genome (outgroup included)
    are not sampled, so a small replicate never holds blank genomes that would
    look identical to EcoSim or leave the root without data.

    With each_gene, the 1-gene level uses every core gene once (one replicate per
    gene) instead of trials_per_count random genes, so genes can be compared.
    With include_full, the whole alignment is also written as
    full_core_genome.fasta (every core gene, gaps included).

    The genes used by every replicate are written to manifest_path
    (default: <output_folder>/rarefaction_genes.csv).
    """
    if gene_counts is None:
        gene_counts = [1, 3, 7, 20, 100]

    if not gene_header or not os.path.exists(gene_header):
        raise FileNotFoundError(
            "Gene coordinate file not found. Rarefaction samples real core genes and needs the "
            "core_alignment_header.embl that Panaroo writes next to core_gene_alignment.aln."
        )

    rng = random.Random(random_seed)
    os.makedirs(output_folder, exist_ok=True)

    # Load and validate input sequences
    # Plain (id, str) pairs: slicing str is far cheaper than Bio.Seq in the inner loop.
    with open(input_fasta) as handle:
        sequences = [(title.split(None, 1)[0], seq) for title, seq in SimpleFastaParser(handle)]
    if not sequences:
        raise ValueError("Input FASTA file is empty or not found.")
    seq_len = len(sequences[0][1])
    if any(len(seq) != seq_len for _, seq in sequences):
        raise ValueError("Not all sequences are the same length. Please align them first.")

    genes = load_gene_coordinates(gene_header)
    if not genes:
        raise ValueError(f"No gene features found in {gene_header}")
    # The genes must tile the alignment exactly; anything else means a header from
    # another run or an entry format this parser doesn't understand.
    ends = [0] + [end for _, _, end in genes[:-1]]
    if any(start != prev for (_, start, _), prev in zip(genes, ends)) or genes[-1][2] != seq_len:
        raise ValueError(
            f"Genes in {os.path.basename(gene_header)} do not tile the {seq_len} bp alignment; "
            "the header must come from the same Panaroo run as the alignment."
        )
    print(f"Loaded {len(genes)} core genes from {os.path.basename(gene_header)}")

    if skip_absent_genes:
        present = [g for g in genes if all(seq[g[1]:g[2]].strip("-") for _, seq in sequences)]
        if len(present) < len(genes):
            print(f"Skipping {len(genes) - len(present)} gene(s) that are all gaps in at least one genome; "
                  f"sampling from {len(present)}.")
        genes = present

    manifest = []
    for gene_count in gene_counts:
        if gene_count > len(genes):
            print(f"⚠️  Skipping g={gene_count}: only {len(genes)} core genes available.")
            continue
        if each_gene and gene_count == 1:
            draws = [[i] for i in range(len(genes))]
            print(f"g=1: one replicate per core gene ({len(draws)} replicates).")
        else:
            draws = [sorted(rng.sample(range(len(genes)), gene_count)) for _ in range(trials_per_count)]
        for trial, chosen in enumerate(draws, start=1):
            name = f"sim_species_g{gene_count}_t{trial}"
            output_file = os.path.join(output_folder, f"{name}.fasta")

            with open(output_file, "w") as f_out:
                for seq_id, seq in sequences:
                    combined = ''.join(seq[genes[i][1]:genes[i][2]] for i in chosen)
                    f_out.write(f">{seq_id}\n{combined}\n")

            manifest.append([name, gene_count, trial,
                             sum(genes[i][2] - genes[i][1] for i in chosen),
                             ";".join(genes[i][0] for i in chosen)])

    if include_full:
        # Same sequences as the input, so step 7b can reuse the step-5 tree for it.
        with open(os.path.join(output_folder, f"{FULL_CORE_NAME}.fasta"), "w") as f_out:
            for seq_id, seq in sequences:
                f_out.write(f">{seq_id}\n{seq}\n")
        all_genes = load_gene_coordinates(gene_header)
        manifest.append([FULL_CORE_NAME, len(all_genes), "", seq_len, "all"])

    manifest_path = manifest_path or os.path.join(output_folder, "rarefaction_genes.csv")
    os.makedirs(os.path.dirname(os.path.abspath(manifest_path)), exist_ok=True)
    with open(manifest_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "gene_count", "trial", "alignment_bp", "genes"])
        writer.writerows(manifest)

    print(f"✅ Done. Generated {len(manifest)} FASTA files in '{output_folder}' (genes listed in {manifest_path}).")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create gene-level rarefaction FASTAs from a Panaroo core alignment")
    parser.add_argument("input_fasta", help="core_gene_alignment.aln (or a sorted/deduplicated copy)")
    parser.add_argument("--gene-header", default=None,
                        help="core_alignment_header.embl (default: next to the alignment)")
    parser.add_argument("--output-folder", default="rarefaction_fastas")
    parser.add_argument("--gene-counts", type=int, nargs="+", default=None)
    parser.add_argument("--trials", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--keep-absent-genes", action="store_true",
                        help="Also sample genes that are all gaps in some genome")
    parser.add_argument("--each-gene", action="store_true",
                        help="At the 1-gene level, use every core gene once instead of random trials")
    parser.add_argument("--no-full-core", action="store_true",
                        help="Don't also write the whole alignment as full_core_genome.fasta")
    args = parser.parse_args()

    create_rarefaction_fastas(
        input_fasta=args.input_fasta,
        output_folder=args.output_folder,
        gene_header=args.gene_header or find_gene_header(args.input_fasta),
        gene_counts=args.gene_counts,
        trials_per_count=args.trials,
        random_seed=args.seed,
        skip_absent_genes=not args.keep_absent_genes,
        each_gene=args.each_gene,
        include_full=not args.no_full_core,
    )
