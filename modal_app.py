"""Run the ecotype pipeline on Modal.

One-time:  modal run modal_app.py::download_bakta_db            (add --light for the small db)
One run:   modal run modal_app.py --species "Treponema paraluiscuniculi" --extra "--sample-size 5"
Batch:     modal run modal_app.py --csv-file species.csv         (one container per species, in parallel)
Results:   modal volume get ecotype-results Treponema_paraluiscuniculi/ecosim_output_core_gene_alignment .
"""
import shlex
import subprocess
import sys
from pathlib import Path

import modal

HERE = Path(__file__).parent
CPUS = 16
MEMORY_MB = 64 * 1024  # each parallel Bakta job loads its own db indexes
TIMEOUT = 24 * 60 * 60  # Modal's maximum

image = (
    # Miniforge gives `conda run -n <env>`, which pipeline.py uses for Bakta and Panaroo.
    modal.Image.from_registry("condaforge/miniforge3:latest")
    .apt_install("git", "make", "gcc", "gfortran")
    .run_commands(
        "conda create -y -n bakta_env -c conda-forge -c bioconda bakta",
        "conda create -y -n panaroo_env -c conda-forge -c bioconda panaroo mafft",
        "conda install -y -n base -c conda-forge -c bioconda ncbi-datasets-cli veryfasttree openjdk",
        "conda clean -afy",
    )
    .pip_install("biopython==1.87", "pandas==2.2.3", "openpyxl==3.1.5")
    # tools/bin is macOS arm64; build the Linux helpers ecosim.jar calls (it looks in <cwd>/bin).
    .add_local_file(HERE / "tools" / "build_ecosim_linux.sh", "/opt/build_ecosim_linux.sh", copy=True)
    .run_commands("bash /opt/build_ecosim_linux.sh /opt/ecosim")
    .env({"ECOSIM_DIR": "/opt/ecosim"})
    # Only the code: local species result folders stay off the image.
    .add_local_file(HERE / "pipeline.py", "/repo/pipeline.py")
    .add_local_dir(HERE / "steps", "/repo/steps", ignore=["__pycache__"])
    .add_local_file(HERE / "tools" / "ecosim.jar", "/repo/tools/ecosim.jar")
)

app = modal.App("ecotype-pipeline", image=image)
bakta_db = modal.Volume.from_name("bakta-db", create_if_missing=True)
results = modal.Volume.from_name("ecotype-results", create_if_missing=True)


@app.function(volumes={"/bakta": bakta_db}, timeout=TIMEOUT)
def download_bakta_db(light: bool = False):
    """Download the Bakta database into the bakta-db volume (/bakta/db or /bakta/db-light)."""
    subprocess.run(
        ["conda", "run", "--no-capture-output", "-n", "bakta_env", "bakta_db", "download",
         "--output", "/bakta", "--type", "light" if light else "full"],
        check=True,
    )
    bakta_db.commit()


@app.function(
    volumes={"/bakta": bakta_db, "/results": results},
    cpu=CPUS, memory=MEMORY_MB, timeout=TIMEOUT,
)
def run_species(species: str, outgroup: str = "", db: str = "db", extra: str = ""):
    """Run pipeline.py for one species; results land in the ecotype-results volume."""
    cmd = [sys.executable, "/repo/pipeline.py", "--species", species,
           "--db", f"/bakta/{db}", "--threads", str(CPUS), "--workdir", "/results",
           *shlex.split(extra)]
    if outgroup:
        cmd += ["--outgroup", outgroup]
    try:
        subprocess.run(cmd, check=True)
    finally:
        results.commit()  # keep partial output on failure so --start-step can resume


@app.local_entrypoint()
def main(species: str = "", csv_file: str = "", outgroup: str = "", db: str = "db", extra: str = ""):
    """--db is db or db-light; --extra is passed to pipeline.py, e.g. "--sample-size 5 --start-step 4"."""
    if csv_file:
        from pipeline import load_species_outgroups_from_csv

        rows = load_species_outgroups_from_csv(csv_file)
        outcomes = run_species.starmap(
            [(sp, og, db, extra) for _, sp, og in rows], return_exceptions=True
        )
        for (row, sp, _), outcome in zip(rows, outcomes):
            print(f"Row {row} {sp}: {'failed: ' + str(outcome) if isinstance(outcome, Exception) else 'ok'}")
    elif species:
        run_species.remote(species, outgroup, db, extra)
    else:
        sys.exit("Pass --species or --csv-file")
