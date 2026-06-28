import os
import subprocess
import tempfile
from pathlib import Path

def download_ncbi_fasta(identifier, output_dir):
    """
    Download a genome FASTA from NCBI given an accession or species name.
    Returns the path to the downloaded FASTA file, or None if failed.
    """
    os.makedirs(output_dir, exist_ok=True)
    # If identifier looks like an accession, use it directly; else treat as species name
    if identifier.startswith("GCA_") or identifier.startswith("GCF_"):
        acc_list = [identifier]
    else:
        # Use datasets summary to get one accession for the species
        cmd = [
            "datasets", "summary", "genome", "taxon", identifier,
            "--assembly-level", "chromosome,complete",
            "--as-json-lines"
        ]
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            for line in result.stdout.splitlines():
                if 'accession' in line or 'assembly_accession' in line:
                    import json
                    record = json.loads(line)
                    acc = record.get("accession") or record.get("assembly_accession")
                    if acc:
                        acc_list = [acc]
                        break
            else:
                return None
        except Exception:
            return None
    # Write accessions to temp file
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as temp_file:
        for acc in acc_list:
            temp_file.write(f"{acc}\n")
        accession_file = temp_file.name
    zip_filename = os.path.join(output_dir, "outgroup_download.zip")
    cmd = [
        "datasets", "download", "genome", "accession",
        "--inputfile", accession_file,
        "--include", "genome",
        "--filename", zip_filename
    ]
    try:
        subprocess.run(cmd, check=True)
        import zipfile
        with zipfile.ZipFile(zip_filename, 'r') as zip_ref:
            zip_ref.extractall(output_dir)
        data_path = os.path.join(output_dir, "ncbi_dataset", "data")
        for root, dirs, files in os.walk(data_path):
            for file in files:
                if file.endswith(".fna"):
                    fasta_path = os.path.join(root, file)
                    out_path = os.path.join(output_dir, "outgroup.fna")
                    os.rename(fasta_path, out_path)
                    return out_path
        return None
    except Exception:
        return None
    finally:
        if os.path.exists(accession_file):
            os.remove(accession_file)
        if os.path.exists(zip_filename):
            os.remove(zip_filename)
