import os
import random
from Bio import SeqIO

def create_rarefaction_fastas(input_fasta,
                              output_folder="rarefaction_fastas",
                              gene_length=1000, 
                              gene_counts=None, 
                              trials_per_count=20, 
                              random_seed=42):
    if gene_counts is None:
        gene_counts = [1, 3, 7, 20, 100]

    random.seed(random_seed)
    os.makedirs(output_folder, exist_ok=True)

    # Load and validate input sequences
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    if not sequences:
        raise ValueError("Input FASTA file is empty or not found.")
    seq_len = len(sequences[0].seq)
    if any(len(record.seq) != seq_len for record in sequences):
        raise ValueError("Not all sequences are the same length. Please align them first.")

    # Main loop
    for gene_count in gene_counts:
        for trial in range(1, trials_per_count + 1):
            starts = sorted(random.sample(range(seq_len - gene_length), gene_count))
            output_file = os.path.join(output_folder, f"sim_species_g{gene_count}_t{trial}.fasta")

            with open(output_file, "w") as f_out:
                for record in sequences:
                    combined = ''.join(str(record.seq[start:start + gene_length]) for start in starts)
                    f_out.write(f">{record.id}\n{combined}\n")

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
