import os
import random
from Bio.SeqIO.FastaIO import SimpleFastaParser

def create_rarefaction_fastas(input_fasta,
                              output_folder="rarefaction_fastas",
                              gene_length=1000, 
                              gene_counts=None, 
                              trials_per_count=20, 
                              random_seed=42):
    if gene_counts is None:
        gene_counts = [1, 3, 7, 20, 100]

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

    # Valid window start positions are 0 .. seq_len - gene_length (inclusive).
    max_start = seq_len - gene_length
    if max_start < 0:
        raise ValueError(
            f"gene_length ({gene_length}) exceeds alignment length ({seq_len}); "
            "cannot sample gene windows."
        )
    num_positions = max_start + 1  # +1 fixes the old off-by-one that dropped the last start
    if max(gene_counts) > num_positions:
        raise ValueError(
            f"Requested up to {max(gene_counts)} gene windows but only {num_positions} "
            f"distinct start positions exist (length {seq_len}, gene_length {gene_length})."
        )

    # Main loop
    for gene_count in gene_counts:
        for trial in range(1, trials_per_count + 1):
            starts = sorted(rng.sample(range(num_positions), gene_count))
            output_file = os.path.join(output_folder, f"sim_species_g{gene_count}_t{trial}.fasta")

            with open(output_file, "w") as f_out:
                for record_id, seq in sequences:
                    combined = ''.join(seq[start:start + gene_length] for start in starts)
                    f_out.write(f">{record_id}\n{combined}\n")

    print(f"✅ Done. Generated {len(gene_counts) * trials_per_count} FASTA files in '{output_folder}'.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python Rarefaction_fasta_creation.py <input.fasta> [output_folder]")
        sys.exit(1)
    create_rarefaction_fastas(
        input_fasta=sys.argv[1],
        output_folder=sys.argv[2] if len(sys.argv) > 2 else "rarefaction_fastas",
    )
