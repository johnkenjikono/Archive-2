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


# Lives in the species folder: input/ (and outgroup.fna with it) is deleted after Bakta by default,
# so this is the only record later steps such as plotting have of which genome the outgroup was.
OUTGROUP_INFO = "outgroup.json"


def _write_outgroup_info(data_path, accession, info_path):
    """Record which genome became outgroup.fna, from the report NCBI ships inside the download.

    That report uses camelCase keys (organismName) where `datasets summary` uses snake_case
    (organism_name), so accept either. Falls back to the accession alone if it can't be read.
    """
    info = {"accession": accession, "organism": "", "strain": ""}
    try:
        with open(os.path.join(data_path, "assembly_data_report.jsonl")) as f:
            record = json.loads(f.readline())
        org = record.get("organism") or {}
        names = org.get("infraspecificNames") or org.get("infraspecific_names") or {}
        info.update(accession=record.get("accession") or accession,
                    organism=org.get("organismName") or org.get("organism_name") or "",
                    strain=names.get("strain") or "")
    except (OSError, ValueError):
        pass
    with open(info_path, "w") as f:
        json.dump(info, f, indent=2)


def outgroup_label(info_path):
    """'Listeria innocua ATCC 33090 (GCF_028596125.1)' from an outgroup.json, or None if absent."""
    try:
        with open(info_path) as f:
            info = json.load(f)
    except (OSError, ValueError):
        return None
    organism, strain, accession = info.get("organism", ""), info.get("strain", ""), info.get("accession", "")
    name = organism if strain and strain in organism else " ".join(p for p in (organism, strain) if p)
    if name and accession:
        return f"{name} ({accession})"
    return name or accession or None


def download_ncbi_fasta(identifier, output_dir, info_path=None):
    """
    Download a genome FASTA from NCBI given an accession or species name.
    Returns the path to the downloaded FASTA file, or None if failed.
    With info_path, also records which genome it was (see OUTGROUP_INFO).
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
                    if info_path:
                        try:  # a missing label must never fail the download itself
                            _write_outgroup_info(data_path, acc_list[0], info_path)
                        except OSError as e:
                            print(f"Warning: could not record the outgroup's identity: {e}")
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
