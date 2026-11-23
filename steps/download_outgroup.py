import json
import os
import random
import shutil
import statistics
import subprocess
import tempfile
import zipfile
from itertools import combinations
from pathlib import Path

from steps.ani import load_ani, run_fastani


def _datasets_summary(taxon, reference_only=False, from_type=False):
    """Run `datasets summary genome taxon` and yield parsed JSON records."""
    cmd = [
        "datasets", "summary", "genome", "taxon", taxon,
        "--assembly-level", "chromosome,complete",
        "--as-json-lines",
    ]
    if reference_only:
        # Server-side filter: avoids pulling metadata for every assembly in the taxon.
        cmd.append("--reference")
    if from_type:
        cmd.append("--from-type")
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


def outgroup_candidates(species_name):
    """
    Return {accession: organism_name}: one type-strain genome for every other species in the genus.

    Uses NCBI's type-material tag; only if the genus has no tagged genome at all does it fall
    back to RefSeq reference genomes as the type-strain proxy. Subspecies count as their species.
    """
    words = species_name.lower().split()
    if len(words) < 2:
        raise RuntimeError(f"'{species_name}' is not a 'Genus species' name")
    genus, species = words[:2]
    records = list(_datasets_summary(genus, from_type=True))
    if not records:
        records = [r for r in _datasets_summary(genus, reference_only=True)
                   if r.get("assembly_info", {}).get("refseq_category") == "reference genome"]
    # One genome per species: the RefSeq reference copy if there is one, then GCF over GCA.
    records.sort(key=lambda r: (r.get("assembly_info", {}).get("refseq_category") != "reference genome",
                                not r.get("accession", "").startswith("GCF_")))
    picked = {}
    for record in records:
        org_name = record.get("organism", {}).get("organism_name", "")
        words = org_name.lower().split()
        accession = record.get("accession")
        if accession and len(words) >= 2 and words[0] == genus and words[1] != species:
            picked.setdefault(words[1], (accession, org_name))
    return dict(picked.values())


def download_genomes(accessions, output_dir):
    """Download genome FASTAs into output_dir as <accession>.fna. Returns {accession: Path}."""
    if not accessions:
        return {}
    os.makedirs(output_dir, exist_ok=True)
    accession_file = os.path.join(output_dir, "accessions.txt")
    zip_filename = os.path.join(output_dir, "genomes_download.zip")
    with open(accession_file, "w") as f:
        f.write("\n".join(accessions) + "\n")
    cmd = ["datasets", "download", "genome", "accession", "--inputfile", accession_file,
           "--include", "genome", "--filename", zip_filename]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        with zipfile.ZipFile(zip_filename) as zip_ref:
            zip_ref.extractall(output_dir)
        paths = {}
        for fna in Path(output_dir, "ncbi_dataset", "data").glob("*/*.fna"):
            paths[fna.parent.name] = fna.rename(Path(output_dir) / f"{fna.parent.name}.fna")
        return paths
    except (subprocess.CalledProcessError, OSError, zipfile.BadZipFile) as e:
        raise RuntimeError(f"could not download genomes from NCBI: {e}")
    finally:
        for path in (accession_file, zip_filename):
            if os.path.exists(path):
                os.remove(path)
        shutil.rmtree(os.path.join(output_dir, "ncbi_dataset"), ignore_errors=True)


# shortcut: the ingroup minimum ANI and the candidate ranking use a 20-genome ingroup sample, so the sample minimum can sit above the true one; measure all ingroup pairs if a too-close outgroup ever slips through
ANI_SAMPLE = 20


def first_passing(ranked, ani_to_ingroup, ingroup_min):
    """
    Return (accession, ANI values) for the first of *ranked* (closest first) whose ANI to every
    ingroup genome is measured and below *ingroup_min*, or (None, None).

    fastANI reports nothing for pairs below ~80% ANI. A missing value (None) therefore means the
    distance was not measured, so that candidate is rejected rather than assumed distant enough.
    """
    for accession in ranked:
        values = ani_to_ingroup(accession)
        if values and all(v is not None and v < ingroup_min for v in values):
            return accession, values
    return None, None


def _fastani(queries, refs, out_tsv, threads):
    if not run_fastani(queries, out_tsv, threads, refs=refs):
        raise RuntimeError("fastANI failed, so the outgroup cannot be measured")
    return load_ani(out_tsv)


def select_outgroup(species_name, ingroup_paths, work_dir, threads=8, seed=42):
    """
    Measured outgroup choice: of the type strains of the other species in the genus, the one with
    the highest fastANI to the ingroup whose ANI to every ingroup genome is still below the lowest
    ANI between two ingroup genomes.

    Returns (accession, organism_name, selection), where selection records the measurements
    (for outgroup.json). Raises RuntimeError when nothing qualifies; there is no unmeasured fallback.
    """
    candidates = outgroup_candidates(species_name)
    if not candidates:
        raise RuntimeError(f"no type-strain genome of another {species_name.split()[0]} species on NCBI")
    ingroup = sorted(Path(p) for p in ingroup_paths if Path(p).stem != "outgroup")
    if len(ingroup) < 2:
        raise RuntimeError(f"need at least 2 ingroup genomes to measure against, found {len(ingroup)}")
    sample = random.Random(seed).sample(ingroup, min(ANI_SAMPLE, len(ingroup)))
    rest = [g for g in ingroup if g not in sample]

    with tempfile.TemporaryDirectory(dir=work_dir) as tmp:
        tmp = Path(tmp)
        print(f"Measuring {len(candidates)} type-strain candidate(s) against {len(sample)} ingroup genome(s) with fastANI.")
        genomes = download_genomes(list(candidates), tmp)
        if not genomes:
            raise RuntimeError("none of the candidate genomes downloaded")
        within = _fastani(sample, sample, tmp / "ingroup.tsv", threads)
        pairs = [within.get(frozenset((a.stem, b.stem))) for a, b in combinations(sample, 2)]
        if None in pairs:
            raise RuntimeError("some ingroup genome pairs fall below fastANI's ~80% range, "
                               "so the ingroup's minimum ANI cannot be measured")
        ingroup_min = min(pairs)

        to_sample = _fastani(list(genomes.values()), sample, tmp / "candidates.tsv", threads)
        sampled = {acc: [to_sample.get(frozenset((acc, g.stem))) for g in sample] for acc in genomes}
        mean_ani = {acc: statistics.mean(m) for acc, v in sampled.items() if (m := [x for x in v if x is not None])}

        def ani_to_ingroup(acc):
            if not rest or first_passing([acc], sampled.get, ingroup_min)[0] is None:
                return sampled[acc]  # already decided by the sample
            ani = _fastani([genomes[acc]], rest, tmp / f"{acc}.tsv", threads)
            return sampled[acc] + [ani.get(frozenset((acc, g.stem))) for g in rest]

        chosen, values = first_passing(sorted(mean_ani, key=mean_ani.get, reverse=True), ani_to_ingroup, ingroup_min)

    selection = {
        "method": "fastANI: closest type strain of another species in the genus whose ANI to every "
                  "ingroup genome is below the ingroup's minimum pairwise ANI",
        "ingroup_genomes": len(ingroup),
        "ingroup_sample": len(sample),
        "ingroup_min_ani": round(ingroup_min, 3),
        "chosen_max_ani": round(max(values), 3) if chosen else None,
        # None: not downloaded, or below fastANI's ~80% range to every sampled ingroup genome
        "candidate_mean_ani": {acc: {"organism": org, "mean_ani_to_sample": round(mean_ani[acc], 3) if acc in mean_ani else None}
                               for acc, org in candidates.items()},
    }
    if not chosen:
        closest = ", ".join(f"{candidates[a]} {mean_ani[a]:.2f}" for a in sorted(mean_ani, key=mean_ani.get, reverse=True)[:3])
        raise RuntimeError(f"no candidate is measurably further from every ingroup genome than the ingroup's "
                           f"minimum pairwise ANI ({ingroup_min:.2f}); closest by mean ANI: {closest or 'none within fastANI range'}")
    return chosen, candidates[chosen], selection


# Lives in the species folder: input/ (and outgroup.fna with it) is deleted after Bakta by default,
# so this is the only record later steps such as plotting have of which genome the outgroup was.
OUTGROUP_INFO = "outgroup.json"


def _write_outgroup_info(data_path, accession, info_path, selection=None):
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
    if selection:
        info["selection"] = selection
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


def download_ncbi_fasta(identifier, output_dir, info_path=None, selection=None):
    """
    Download a genome FASTA from NCBI given an accession or species name.
    Returns the path to the downloaded FASTA file, or None if failed.
    With info_path, also records which genome it was, and how it was chosen (see OUTGROUP_INFO).
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
                            _write_outgroup_info(data_path, acc_list[0], info_path, selection)
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
