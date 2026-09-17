import json
import os
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path


def _datasets_summary(taxon, reference_only=False):
    """Run `datasets summary genome taxon` and yield parsed JSON records."""
    cmd = [
        "datasets", "summary", "genome", "taxon", taxon,
        "--assembly-level", "chromosome,complete",
        "--as-json-lines",
    ]
    if reference_only:
        # Server-side filter: avoids pulling metadata for every assembly in the taxon.
        cmd.append("--reference")
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError:
        return
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            yield json.loads(line)
        except json.JSONDecodeError:
            continue


def find_typestrain(species_name):
    """
    Return (accession, organism_name) for the typestrain of *species_name*.

    Uses the NCBI RefSeq "reference genome" designation as a proxy for the
    type strain — this is correct for the vast majority of species.
    Returns (None, None) if nothing is found.
    """
    for record in _datasets_summary(species_name, reference_only=True):
        org_name = record.get("organism", {}).get("organism_name", "")
        refseq_cat = record.get("assembly_info", {}).get("refseq_category", "")
        accession = record.get("accession")
        if accession and refseq_cat == "reference genome":
            return accession, org_name
    return None, None


def find_outgroup_species(species_name):
    """
    Return (accession, organism_name) for a suitable outgroup: a reference
    genome from the same genus but a different species.

    Prefers RefSeq reference genomes; falls back to any complete genome in
    the genus if no reference is available.
    Returns (None, None) if nothing is found.
    """
    genus = species_name.split()[0]
    species_lower = species_name.lower()

    def other_species(records):
        for record in records:
            org_name = record.get("organism", {}).get("organism_name", "")
            accession = record.get("accession")
            if not accession or not org_name:
                continue
            if not org_name.lower().startswith(genus.lower() + " "):
                continue
            if org_name.lower().startswith(species_lower):
                continue
            yield record, accession, org_name

    # Cheap query first: reference genomes only.
    for record, accession, org_name in other_species(_datasets_summary(genus, reference_only=True)):
        if record.get("assembly_info", {}).get("refseq_category") == "reference genome":
            return accession, org_name

    # Fallback: first complete genome of another species in the genus.
    for _, accession, org_name in other_species(_datasets_summary(genus)):
        return accession, org_name
    return None, None


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
        acc_list = None
        for record in _datasets_summary(identifier):
            acc = record.get("accession")
            if acc:
                acc_list = [acc]
                break
        if not acc_list:
            return None

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as tmp:
        for acc in acc_list:
            tmp.write(f"{acc}\n")
        accession_file = tmp.name

    zip_filename = os.path.join(output_dir, "outgroup_download.zip")
    cmd = [
        "datasets", "download", "genome", "accession",
        "--inputfile", accession_file,
        "--include", "genome",
        "--filename", zip_filename,
    ]
    try:
        subprocess.run(cmd, check=True)
        with zipfile.ZipFile(zip_filename, "r") as zip_ref:
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
        shutil.rmtree(os.path.join(output_dir, "ncbi_dataset"), ignore_errors=True)
