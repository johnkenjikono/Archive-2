import sys
from Bio import SeqIO

def move_outgroup_to_top(input_fasta, output_fasta, outgroup_name):
    """
    Reads a FASTA file, finds the sequence matching the outgroup name,
    and writes a new FASTA file with the outgroup as the very first sequence.
    """
    print(f"Reading {input_fasta}...")
    
    outgroup_record = None
    other_records = []
    
    # Read through the FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        # We check if the outgroup name is in the record ID or description
        # (Roary usually sets the ID to the base filename of the Bakta .gff)
        if outgroup_name.lower() in record.id.lower() or outgroup_name.lower() in record.description.lower():
            outgroup_record = record
        else:
            other_records.append(record)
            
    if outgroup_record is None:
        print(f"❌ Could not find any sequence matching '{outgroup_name}' in the file.")
        print("Please check the exact spelling of your outgroup in the FASTA file.")
        return
        
    # Put outgroup first
    final_records = [outgroup_record] + other_records
    
    print(f"✅ Found outgroup: {outgroup_record.id}")
    print(f"Writing {len(final_records)} sequences to {output_fasta}...")
    
    # Write the new file
    with open(output_fasta, "w") as out_f:
        SeqIO.write(final_records, out_f, "fasta")
        
    print("Done! The outgroup is now at the top of the file.")

if __name__ == "__main__":
    # Example usage:
    # python move_outgroup_to_top.py input.fasta output.fasta MyOutgroupName
    
    if len(sys.argv) != 4:
        print("Usage: python move_outgroup_to_top.py <input_fasta> <output_fasta> <outgroup_name>")
        print("Example: python move_outgroup_to_top.py core_gene_alignment.aln core_gene_sorted.fasta Strain_XYZ")
        sys.exit(1)
        
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    outgroup_name = sys.argv[3]
    
    move_outgroup_to_top(input_fasta, output_fasta, outgroup_name)
