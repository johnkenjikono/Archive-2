"""
Fill blank outgroups in a pipeline batch CSV (column 1 = species, column 4 = outgroup).

Rows that already have an outgroup are never touched, so re-running is safe.
Each blank is chosen by fastANI (steps/download_outgroup.select_outgroup) against a small
sample of the species' NCBI genomes, downloaded to a temporary folder and then deleted.
How it was chosen goes in column 5 ("Outgroup Source"), which pipeline.py ignores.

  python fill_outgroups.py species.csv
"""
import argparse
import csv
import sys
import tempfile
from pathlib import Path

from pipeline import get_random_accessions
from steps.download_outgroup import ANI_SAMPLE, download_genomes, select_outgroup
from steps.dotenv import load_dotenv

SPECIES_COL, OUTGROUP_COL, SOURCE_COL = 0, 3, 4


def choose_outgroup(species_name):
    """Return (accession, source). accession is '' if no candidate could be measured as suitable."""
    with tempfile.TemporaryDirectory() as tmp:
        try:
            ingroup = download_genomes(get_random_accessions(species_name, ANI_SAMPLE), tmp)
            acc, org, sel = select_outgroup(species_name, list(ingroup.values()), tmp)
        except RuntimeError as exc:
            return "", f"not chosen: {exc}"
    return acc, (f"fastANI: {org}; max ANI to ingroup {sel['chosen_max_ani']} < ingroup minimum "
                 f"{sel['ingroup_min_ani']} ({sel['ingroup_sample']} NCBI genomes sampled)")


def fill_rows(rows, chooser=choose_outgroup):
    """Fill blank outgroups in *rows* (header first) in place. Returns the number of rows filled."""
    header = rows[0]
    header += [""] * (SOURCE_COL + 1 - len(header))
    header[SOURCE_COL] = header[SOURCE_COL] or "Outgroup Source"

    filled = 0
    for row in rows[1:]:
        if not row or not row[SPECIES_COL].strip():
            continue
        row += [""] * (SOURCE_COL + 1 - len(row))
        if row[OUTGROUP_COL].strip():
            continue
        species = row[SPECIES_COL].strip()
        acc, source = chooser(species)
        row[OUTGROUP_COL], row[SOURCE_COL] = acc, source
        print(f"{species}: {acc or '(none)'}  [{source}]")
        filled += bool(acc)
    return filled


def _self_test():
    rows = [["Species Name", "Genomes", "Notes", "Outgroup"],
            ["Bacillus subtilis", "", "", "Bacillus licheniformis"],
            ["Treponema paraluiscuniculi"],
            []]
    fake = lambda species: ("GCF_1", "fastANI: x")
    assert fill_rows(rows, chooser=fake) == 1
    assert rows[0][4] == "Outgroup Source"
    assert rows[1][3] == "Bacillus licheniformis" and rows[1][4] == ""  # existing outgroup untouched
    assert rows[2][3] == "GCF_1" and rows[2][4] == "fastANI: x"
    assert fill_rows(rows, chooser=fake) == 0  # idempotent
    print("self-test passed")


def main():
    parser = argparse.ArgumentParser(description="Fill blank outgroups (column 4) in a pipeline batch CSV.")
    parser.add_argument("csv_file", nargs="?", help="Batch CSV (column 1 = species, column 4 = outgroup)")
    parser.add_argument("-o", "--output", help="Output CSV (default: <input>_outgroups.csv)")
    parser.add_argument("--self-test", action="store_true", help="Run the offline self-test and exit")
    args = parser.parse_args()

    if args.self_test:
        _self_test()
        return
    if not args.csv_file:
        parser.error("csv_file is required")

    in_path = Path(args.csv_file)
    out_path = Path(args.output) if args.output else in_path.with_name(f"{in_path.stem}_outgroups.csv")

    with open(in_path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.reader(f))
    if len(rows) < 2:
        print(f"No species rows found in {in_path}", file=sys.stderr)
        sys.exit(1)

    filled = fill_rows(rows)

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows(rows)
    print(f"\nFilled {filled} outgroup(s) measured by fastANI. Wrote {out_path}")
    print("Review the 'Outgroup Source' column, then run: python pipeline.py --csv-file", out_path)


if __name__ == "__main__":
    load_dotenv()
    main()
