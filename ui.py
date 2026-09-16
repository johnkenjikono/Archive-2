from __future__ import annotations

import os
import queue
import shutil
import signal
import subprocess
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

FASTA_SUFFIXES = (".fna", ".fasta", ".fa")

Mode = Literal["species", "csv", "local"]


@dataclass
class FormValues:
    mode: Mode
    species: str = ""
    outgroup: str = ""
    sample_size: int = 200
    seed: int = 42
    csv_path: str = ""
    local_folder: str = ""
    workdir: str = "."
    db: str = ""
    threads: int = 12
    bakta_jobs: int | None = None
    start_step: int = 1
    input_fasta: str = ""
    outgroup_id: str = ""
    setup_only: bool = False


def build_command(values: FormValues, python_exe: str, pipeline_py: Path) -> list[str]:
    argv = [python_exe, str(pipeline_py)]
    if values.mode == "csv":
        argv.extend(["--csv-file", values.csv_path])
        argv.extend(["--sample-size", str(values.sample_size), "--seed", str(values.seed)])
    elif values.mode == "species":
        argv.extend(["--species", values.species])
        if values.outgroup.strip():
            argv.extend(["--outgroup", values.outgroup.strip()])
        argv.extend(["--sample-size", str(values.sample_size), "--seed", str(values.seed)])
    else:
        argv.extend(["--species", values.species])

    argv.extend(
        [
            "--workdir",
            values.workdir,
            "--threads",
            str(values.threads),
            "--start-step",
            str(values.start_step),
        ]
    )
    if values.db.strip():
        argv.extend(["--db", values.db.strip()])
    if values.bakta_jobs is not None:
        argv.extend(["--bakta-jobs", str(values.bakta_jobs)])
    if values.input_fasta.strip():
        argv.extend(["--input-fasta", values.input_fasta.strip()])
    if values.outgroup_id.strip():
        argv.extend(["--outgroup-id", values.outgroup_id.strip()])
    if values.setup_only:
        argv.append("--setup-only")
    return argv


def species_folder_name(species: str) -> str:
    return species.replace(" ", "_")


def local_input_dir(values: FormValues) -> Path:
    return Path(values.workdir) / species_folder_name(values.species) / "input"


def results_dir(values: FormValues) -> Path:
    if values.mode == "csv":
        return Path(values.workdir)
    return Path(values.workdir) / species_folder_name(values.species)


def find_ecotype_summaries(root: Path) -> list[Path]:
    direct = sorted(root.glob("ecosim_output_*/ecotype_summary.csv"))
    nested = sorted(root.glob("*/ecosim_output_*/ecotype_summary.csv"))
    seen: set[Path] = set()
    out: list[Path] = []
    for path in direct + nested:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            out.append(path)
    return out


def copy_local_fastas(src_folder: Path, dest_input: Path) -> list[Path]:
    if not src_folder.is_dir():
        raise FileNotFoundError(f"Genome folder does not exist: {src_folder}")
    dest_input.mkdir(parents=True, exist_ok=True)
    copied: list[Path] = []
    for path in sorted(src_folder.iterdir()):
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES:
            dest = dest_input / path.name
            shutil.copy2(path, dest)
            copied.append(dest)
    return copied


def default_alignment_path(values: FormValues) -> Path:
    return (
        Path(values.workdir)
        / species_folder_name(values.species)
        / "output_roary"
        / "results"
        / "core_gene_alignment.aln"
    )


def validate_form(values: FormValues) -> list[str]:
    errors: list[str] = []
    if values.mode == "species":
        if not values.species.strip():
            errors.append("Enter a species name.")
    elif values.mode == "csv":
        if not values.csv_path.strip():
            errors.append("Choose a CSV file.")
        elif not Path(values.csv_path).is_file():
            errors.append("CSV file does not exist.")
    elif values.mode == "local":
        if not values.species.strip():
            errors.append("Enter a species name.")
        if not values.local_folder.strip():
            errors.append("Choose a folder of genome files.")
        elif not Path(values.local_folder).is_dir():
            errors.append("Genome folder does not exist.")
        if values.start_step < 2:
            errors.append(
                "Local genome folder runs start at step 2 (Bakta). Choose 2 or later."
            )

    if values.start_step <= 2:
        if not values.db.strip():
            errors.append("Choose a Bakta database folder (required for steps 1–2).")
        elif not Path(values.db).is_dir():
            errors.append("Bakta database folder does not exist.")

    if values.input_fasta.strip():
        if not Path(values.input_fasta).is_file():
            errors.append("Core-alignment FASTA does not exist.")
    elif values.start_step >= 4 and values.mode != "csv":
        if not default_alignment_path(values).is_file():
            errors.append("Choose a core-alignment FASTA, or lower the start step.")
    return errors


class PipelineRunner:
    def __init__(self) -> None:
        self._proc: subprocess.Popen[str] | None = None
        self._thread: threading.Thread | None = None
        self.was_stopped = False

    @property
    def running(self) -> bool:
        return self._proc is not None and self._proc.poll() is None

    def start(self, argv: list[str], cwd: Path, events: queue.Queue) -> None:
        if self.running:
            return
        self.was_stopped = False
        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"
        self._proc = subprocess.Popen(
            argv,
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=True,
            env=env,
        )
        self._thread = threading.Thread(
            target=self._pump, args=(events,), daemon=True
        )
        self._thread.start()

    def stop(self) -> None:
        proc = self._proc
        if proc is None or proc.poll() is not None:
            return
        self.was_stopped = True
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        except ProcessLookupError:
            return
        deadline = time.time() + 5
        while time.time() < deadline:
            if proc.poll() is not None:
                return
            time.sleep(0.1)
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except ProcessLookupError:
            return

    def _pump(self, events: queue.Queue) -> None:
        proc = self._proc
        assert proc is not None
        try:
            assert proc.stdout is not None
            for line in proc.stdout:
                events.put(("line", line.rstrip("\n")))
        finally:
            code = proc.wait()
            events.put(("done", code))
