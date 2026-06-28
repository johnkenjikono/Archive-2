#!/usr/bin/env python3
import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def is_fasta_like(file_path: Path) -> bool:
    with file_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            return stripped.startswith(">")
    return False


def parse_clustal_alignment(input_aln: Path) -> list[tuple[str, str]]:
    sequences: dict[str, list[str]] = {}

    with input_aln.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue

            upper = stripped.upper()
            if upper.startswith("CLUSTAL") or upper.startswith("MUSCLE"):
                continue

            if stripped.startswith(("*", ":", ".")):
                continue

            if line.startswith((" ", "\t")):
                continue

            parts = stripped.split()
            if len(parts) < 2:
                continue

            seq_id = parts[0]
            seq_chunk = parts[1]

            sequences.setdefault(seq_id, []).append(seq_chunk)

    joined = [(seq_id, "".join(chunks)) for seq_id, chunks in sequences.items()]
    if not joined:
        raise ValueError("No sequence blocks found in CLUSTAL alignment.")

    lengths = {len(seq) for _, seq in joined}
    if len(lengths) != 1:
        raise ValueError("Parsed CLUSTAL sequences have inconsistent lengths.")

    return joined


def write_fasta(records: list[tuple[str, str]], output_fasta: Path) -> None:
    with output_fasta.open("w", encoding="utf-8") as handle:
        for seq_id, sequence in records:
            handle.write(f">{seq_id}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def aln_to_fasta(input_aln: Path, output_fasta: Path) -> None:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    if is_fasta_like(input_aln):
        shutil.copyfile(input_aln, output_fasta)
        return

    clustal_records = parse_clustal_alignment(input_aln)
    write_fasta(clustal_records, output_fasta)


def run_archiv_pipeline(archiv_pipeline_py: Path, input_fasta: Path, start_step: int) -> int:
    cmd = [
        sys.executable,
        str(archiv_pipeline_py),
        str(input_fasta),
        "--start-step",
        str(start_step),
    ]

    print("Running Archiv pipeline:")
    print(" ".join(cmd))

    result = subprocess.run(cmd, cwd=archiv_pipeline_py.parent)
    return result.returncode


def run_main_pipeline(
    main_pipeline_py: Path,
    db_path: Path,
    workdir: Path,
    threads: int,
) -> int:
    cmd = [
        sys.executable,
        str(main_pipeline_py),
        "--db",
        str(db_path),
        "--workdir",
        str(workdir),
        "--threads",
        str(threads),
    ]

    print("Running Bakta -> Roary pipeline:")
    print(" ".join(cmd))

    result = subprocess.run(cmd, cwd=main_pipeline_py.parent)
    return result.returncode


def main() -> int:
    repo_root = Path(__file__).resolve().parent
    default_main_pipeline = repo_root / "pipeline.py"
    default_workdir = repo_root
    default_input_aln = default_workdir / "output_roary" / "results" / "core_gene_alignment.aln"
    default_output_fasta = repo_root / "Archiv" / "core_gene_alignment.fasta"
    default_archiv_pipeline = repo_root / "Archiv" / "pipeline.py"

    parser = argparse.ArgumentParser(
        description=(
            "Run Bakta->Roary (pipeline.py), then convert core_gene_alignment.aln to FASTA, "
            "then run Archiv/pipeline.py."
        )
    )
    parser.add_argument(
        "--db",
        type=Path,
        required=True,
        help="Path to Bakta database for the first-stage Bakta->Roary pipeline",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads for Bakta/Roary stage (default: 1)",
    )
    parser.add_argument(
        "--workdir",
        type=Path,
        default=default_workdir,
        help=f"Working directory for first-stage pipeline (default: {default_workdir})",
    )
    parser.add_argument(
        "--main-pipeline",
        type=Path,
        default=default_main_pipeline,
        help=f"Path to top-level pipeline.py (default: {default_main_pipeline})",
    )
    parser.add_argument(
        "--skip-main-pipeline",
        action="store_true",
        help="Skip Bakta->Roary stage and use an existing alignment file",
    )
    parser.add_argument(
        "--input-aln",
        type=Path,
        default=default_input_aln,
        help=f"Path to Roary alignment file (default: {default_input_aln})",
    )
    parser.add_argument(
        "--output-fasta",
        type=Path,
        default=default_output_fasta,
        help=f"Path to write FASTA file (default: {default_output_fasta})",
    )
    parser.add_argument(
        "--archiv-pipeline",
        type=Path,
        default=default_archiv_pipeline,
        help=f"Path to Archiv pipeline.py (default: {default_archiv_pipeline})",
    )
    parser.add_argument(
        "--start-step",
        type=int,
        default=1,
        choices=range(1, 7),
        help="Archiv pipeline step to start from (1-6, default 1)",
    )
    args = parser.parse_args()

    main_pipeline = args.main_pipeline.resolve()
    workdir = args.workdir.resolve()
    db_path = args.db.resolve()
    input_aln = args.input_aln.resolve()
    output_fasta = args.output_fasta.resolve()
    archiv_pipeline = args.archiv_pipeline.resolve()

    if not args.skip_main_pipeline:
        if not main_pipeline.exists():
            print(f"Error: top-level pipeline not found: {main_pipeline}")
            return 1

        if not db_path.exists():
            print(f"Error: Bakta DB path not found: {db_path}")
            return 1

        first_stage_code = run_main_pipeline(
            main_pipeline_py=main_pipeline,
            db_path=db_path,
            workdir=workdir,
            threads=args.threads,
        )
        if first_stage_code != 0:
            print("Error: Bakta->Roary stage failed; stopping.")
            return first_stage_code

    if not input_aln.exists():
        print(f"Error: alignment file not found: {input_aln}")
        return 1

    if not archiv_pipeline.exists():
        print(f"Error: Archiv pipeline not found: {archiv_pipeline}")
        return 1

    print(f"Converting alignment to FASTA:\n  {input_aln}\n  -> {output_fasta}")
    try:
        aln_to_fasta(input_aln, output_fasta)
    except Exception as exc:
        print(f"Error converting alignment: {exc}")
        return 1

    return run_archiv_pipeline(archiv_pipeline, output_fasta, args.start_step)


if __name__ == "__main__":
    raise SystemExit(main())
