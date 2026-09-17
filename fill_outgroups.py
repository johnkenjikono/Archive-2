"""
Fill blank outgroups in a pipeline batch CSV (column 1 = species, column 4 = outgroup).

Rows that already have an outgroup are never touched, so re-running is safe.
Each filled row records how it was chosen in column 5 ("Outgroup Source"),
which pipeline.py ignores.

  python fill_outgroups.py species.csv            # NCBI heuristic only (deterministic)
  python fill_outgroups.py species.csv --llm      # Claude picks from NCBI candidates
"""
import argparse
import csv
import json
import sys
from pathlib import Path

from steps.download_outgroup import _datasets_summary, find_outgroup_species

SPECIES_COL, OUTGROUP_COL, SOURCE_COL = 0, 3, 4
LLM_MODEL = "claude-opus-5"


def outgroup_candidates(species_name):
    """Return {organism_name: accession} for other species in the genus, reference genomes first."""
    genus = species_name.split()[0].lower()
    candidates = {}
    for reference_only in (True, False):
        for record in _datasets_summary(genus, reference_only=reference_only):
            org = record.get("organism", {}).get("organism_name", "")
            acc = record.get("accession")
            if acc and org.lower().startswith(genus + " ") and not org.lower().startswith(species_name.lower()):
                candidates.setdefault(org, acc)
        if candidates:
            break  # ponytail: full-genus query only when no reference genomes exist; it is slow for big genera
    return candidates


def pick_with_claude(species_name, candidates):
    """Return (accession, source) chosen by Claude from *candidates*, or (None, reason) on failure."""
    import anthropic  # only needed with --llm

    names = sorted(candidates)
    response = anthropic.Anthropic().beta.messages.create(
        model=LLM_MODEL,
        max_tokens=4000,
        betas=["server-side-fallback-2026-07-01"],
        fallbacks="default",
        output_config={
            "effort": "medium",
            "format": {
                "type": "json_schema",
                "schema": {
                    "type": "object",
                    "properties": {
                        "organism_name": {"type": "string", "enum": names},
                        "reason": {"type": "string"},
                    },
                    "required": ["organism_name", "reason"],
                    "additionalProperties": False,
                },
            },
        },
        messages=[{"role": "user", "content": (
            f"Ingroup species: {species_name}\n"
            "Choose the best outgroup for rooting a core-genome phylogeny of this species: "
            "the closest relative that is clearly outside the ingroup, so core genes still align.\n"
            "Candidates (NCBI genomes in the same genus):\n" + "\n".join(names) +
            "\n\nGive a one-sentence reason."
        )}],
    )
    if response.stop_reason != "end_turn":
        return None, f"llm stopped: {response.stop_reason}"
    text = next(b.text for b in response.content if b.type == "text")
    pick = json.loads(text)
    org = pick["organism_name"]
    if org not in candidates:  # enum should guarantee this; never trust an unverified accession
        return None, f"llm picked unknown organism: {org}"
    return candidates[org], f"{response.model}: {org}: {pick['reason']}"


def choose_outgroup(species_name, use_llm):
    """Return (accession, source). accession is '' if nothing was found."""
    if use_llm:
        candidates = outgroup_candidates(species_name)
        if candidates:
            try:
                acc, source = pick_with_claude(species_name, candidates)
            except Exception as exc:  # API/network errors: fall back, but record it
                acc, source = None, f"llm error: {exc}"
            if acc:
                return acc, source
            failure = source
        else:
            failure = "no candidates in genus"
        acc, org = find_outgroup_species(species_name)
        return (acc or ""), f"heuristic ({failure}): {org}" if acc else f"not found ({failure})"

    acc, org = find_outgroup_species(species_name)
    return (acc, f"heuristic: {org}") if acc else ("", "not found")


def fill_rows(rows, use_llm, chooser=choose_outgroup):
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
        acc, source = chooser(species, use_llm)
        row[OUTGROUP_COL], row[SOURCE_COL] = acc, source
        print(f"{species}: {acc or '(none)'}  [{source}]")
        filled += bool(acc)
    return filled


def _self_test():
    rows = [["Species Name", "Genomes", "Notes", "Outgroup"],
            ["Bacillus subtilis", "", "", "Bacillus licheniformis"],
            ["Treponema paraluiscuniculi"],
            []]
    fake = lambda species, use_llm: ("GCF_1", f"{'llm' if use_llm else 'heuristic'}: x")
    assert fill_rows(rows, use_llm=True, chooser=fake) == 1
    assert rows[0][4] == "Outgroup Source"
    assert rows[1][3] == "Bacillus licheniformis" and rows[1][4] == ""  # existing outgroup untouched
    assert rows[2][3] == "GCF_1" and rows[2][4] == "llm: x"
    assert fill_rows(rows, use_llm=False, chooser=fake) == 0  # idempotent
    print("self-test passed")


def main():
    parser = argparse.ArgumentParser(description="Fill blank outgroups (column 4) in a pipeline batch CSV.")
    parser.add_argument("csv_file", nargs="?", help="Batch CSV (column 1 = species, column 4 = outgroup)")
    parser.add_argument("-o", "--output", help="Output CSV (default: <input>_outgroups.csv)")
    parser.add_argument("--llm", action="store_true",
                        help=f"Let Claude ({LLM_MODEL}) choose from NCBI candidates. Needs ANTHROPIC_API_KEY "
                             "and `pip install anthropic`. Off by default: the NCBI heuristic is used.")
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

    filled = fill_rows(rows, args.llm)

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f).writerows(rows)
    print(f"\nFilled {filled} outgroup(s) using {'Claude' if args.llm else 'NCBI heuristic'}. Wrote {out_path}")
    print("Review the 'Outgroup Source' column, then run: python pipeline.py --csv-file", out_path)


if __name__ == "__main__":
    main()
