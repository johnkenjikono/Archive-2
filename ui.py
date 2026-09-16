from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

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


def species_folder_name(species: str) -> str:
    return species.replace(" ", "_")


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
