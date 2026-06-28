import sys
import re
from Bio import SeqIO

def get_numeric_value(text):
    """
    Extracts all digits from a string and converts them to a single integer.
    Example: 'GCA_001687385.1_ASM168738v1_genomic' becomes 168738511687381
    """
    digits = re.sub(r'\D', '', text)
    if not digits:
        return -1
    return int(digits)

def _normalize_record_id(text):
    return text.strip().split()[0]


def move_largest_numeric_to_top(input_fasta, output_fasta, target_id=None):
    print(f"Reading {input_fasta}...")
    
    records = list(SeqIO.parse(input_fasta, "fasta"))
    if not records:
        print("❌ No sequences found in the file!")
        return False
    
    # CRITICAL: Validate alignment (all sequences must be same length for tree building)
    lengths = set(len(record.seq) for record in records)
    if len(lengths) > 1:
        print(f"⚠️  WARNING: Sequences have different lengths: {sorted(lengths)}")
        print(f"    This will cause FastTree to CRASH!")
        print(f"    Please inspect the Roary core-gene alignment before proceeding")
        return False
    
    seq_len = list(lengths)[0]
    print(f"✅ All sequences aligned ({len(records)} sequences, length={seq_len})")

    selected_record = None
    if target_id:
        normalized_target = _normalize_record_id(target_id)
        for record in records:
            record_id = _normalize_record_id(record.id)
            if record_id == normalized_target or normalized_target in record_id:
                selected_record = record
                break

        if selected_record:
            print(f"✅ Found target sequence to move to top: {selected_record.id}")

    if selected_record is None:
        # Fall back to the record with the largest numeric value in its ID.
        selected_record = max(records, key=lambda r: get_numeric_value(r.id))
        print(f"✅ Found sequence with largest numerical value: {selected_record.id}")

    # Put the selected record at the top by filtering it out from the others and prepending it.
    other_records = [r for r in records if r.id != selected_record.id]
    final_records = [selected_record] + other_records
    
    print(f"Writing {len(final_records)} sequences to {output_fasta}...")
    
    # Write the new file
    with open(output_fasta, "w") as out_f:
        SeqIO.write(final_records, out_f, "fasta")
        
    print("✅ Done! The selected sequence is now at the top of the file.")
    return True

if __name__ == "__main__":
    # Example usage:
    # python move_largest_numeric.py input.fasta output.fasta
    
    if len(sys.argv) != 3:
        print("Usage: python move_largest_numeric.py <input_fasta> <output_fasta>")
        print("Example: python move_largest_numeric.py Bordetella_pertussis.fasta Bordetella_pertussis_sorted.fasta")
        sys.exit(1)
        
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    
    move_largest_numeric_to_top(input_fasta, output_fasta)
