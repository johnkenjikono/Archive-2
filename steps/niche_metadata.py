"""
Do ecotypes line up with where the genomes came from?

Reads ecotype_membership.csv from an EcoSim output folder, fetches each genome's
NCBI BioSample metadata, has Claude normalize the free-text host / isolation
source fields into a small set of labels, and tests ecotype x label association.

Writes into the EcoSim output folder:
  genome_metadata.csv           one row per genome: host, niche, country, year + the raw_* fields Claude saw
  ecotype_x_<variable>.csv      ecotype x category counts (host, niche, country)
  ecotype_niche_association.csv chi-squared, Cramer's V and permutation p-value per variable

Needs ANTHROPIC_API_KEY and `pip install anthropic`.

  python steps/niche_metadata.py ecosim_output_dir "Streptococcus pyogenes"
"""
import argparse
import csv
import json
import random
import re
import subprocess
import tempfile
from collections import Counter
from math import sqrt
from pathlib import Path

try:
    from steps.parsing import FULL_CORE_XML
except ImportError:  # run as a script: python steps/niche_metadata.py
    from parsing import FULL_CORE_XML

LLM_MODEL = "claude-sonnet-5-5"
# Fields Claude sees; country and year are parsed without it.
LLM_FIELDS = ("host", "host_disease", "isolation_source", "sample_type", "source_type",
              "env_broad_scale", "env_local_scale", "env_medium", "host_tissue_sampled", "host_body_site")
CHUNK = 200  # distinct field combinations per Claude request
PERMUTATIONS = 2000
UNKNOWN = "unknown"


def load_membership(membership_csv):
    """Return (xml file used, {accession: ecotype}): the full-core run, else the most-genes replicate."""
    with open(membership_csv, newline="") as f:
        rows = list(csv.DictReader(f))
    files = {r["file"]: r for r in rows}
    if FULL_CORE_XML in files:
        chosen = FULL_CORE_XML
    else:  # ponytail: no full-core run, so use the replicate with the most genes (lowest trial on ties)
        chosen = max(files, key=lambda name: (int(files[name]["gene_count"] or 0), -int(files[name]["trial"] or 0)))
    return chosen, {r["taxon"]: r["ecotype"] for r in rows if r["file"] == chosen}


def fetch_biosamples(accessions):
    """Return {accession: {attribute: value}} from `datasets summary genome accession`."""
    with tempfile.NamedTemporaryFile("w", suffix=".txt") as acc_file:
        acc_file.write("\n".join(accessions))
        acc_file.flush()
        out = subprocess.run(["datasets", "summary", "genome", "accession", "--inputfile", acc_file.name,
                              "--as-json-lines"], check=True, capture_output=True, text=True).stdout
    meta = {}
    for line in out.splitlines():
        record = json.loads(line)
        biosample = record.get("assembly_info", {}).get("biosample", {})
        meta[record["accession"]] = {a["name"]: a.get("value", "") for a in biosample.get("attributes", [])}
    return meta


def country_and_year(attrs):
    country = attrs.get("geo_loc_name", "").split(":")[0].strip()
    year = re.search(r"\b(1[89]|20)\d\d\b", attrs.get("collection_date", ""))
    missing = {"", "missing", "not collected", "not provided", "not applicable", "unknown"}
    return (UNKNOWN if country.lower() in missing else country), (year.group(0) if year else "")


def _llm_key(attrs):
    return tuple((k, attrs[k]) for k in LLM_FIELDS if attrs.get(k))


def label_with_claude(species, keys):
    """Return {key: (host, niche)} for each distinct tuple of (field, value) pairs."""
    import anthropic  # only needed for this step

    client = anthropic.Anthropic()
    labels = {}
    for start in range(0, len(keys), CHUNK):
        chunk = keys[start:start + CHUNK]
        vocabulary = sorted({niche for _, niche in labels.values()})
        lines = [f"{i}: " + ("; ".join(f"{k}={v}" for k, v in key) or "(no metadata)") for i, key in enumerate(chunk)]
        response = client.beta.messages.create(
            model=LLM_MODEL,
            max_tokens=16000,
            betas=["server-side-fallback-2026-07-01"],
            fallbacks="default",
            output_config={
                "effort": "low",
                "format": {"type": "json_schema", "schema": {
                    "type": "object",
                    "properties": {"labels": {"type": "array", "items": {
                        "type": "object",
                        "properties": {"id": {"type": "integer"}, "host": {"type": "string"},
                                       "niche": {"type": "string"}},
                        "required": ["id", "host", "niche"],
                        "additionalProperties": False,
                    }}},
                    "required": ["labels"],
                    "additionalProperties": False,
                }},
            },
            messages=[{"role": "user", "content": (
                f"These are NCBI BioSample metadata for {species} genomes, one distinct record per line.\n"
                "For each line give:\n"
                "- host: the host as a common lowercase name (human, cattle, chicken, pig, ...), "
                "'environment' if not from a host, or 'unknown'.\n"
                "- niche: where in or around the host the isolate came from, as a short lowercase label "
                "(e.g. blood, respiratory tract, skin/soft tissue, gastrointestinal, urogenital, food, soil, water). "
                "Use as few distinct labels as the data allow, about 10 at most, so they can be compared. "
                "Use 'unknown' when the fields are missing or say nothing about the source.\n"
                + (f"Reuse these niche labels where they fit: {', '.join(vocabulary)}\n" if vocabulary else "")
                + "Do not guess beyond what the fields say.\n\n" + "\n".join(lines)
            )}],
        )
        if response.stop_reason != "end_turn":
            raise RuntimeError(f"Claude stopped early: {response.stop_reason}")
        text = next(b.text for b in response.content if b.type == "text")
        for item in json.loads(text)["labels"]:
            if 0 <= item["id"] < len(chunk):
                labels[chunk[item["id"]]] = (item["host"].strip().lower() or UNKNOWN,
                                             item["niche"].strip().lower() or UNKNOWN)
    return labels


def chi2(pairs):
    """Chi-squared statistic of the contingency table built from (row, col) pairs."""
    n = len(pairs)
    cells = Counter(pairs)
    rows = Counter(r for r, _ in pairs)
    cols = Counter(c for _, c in pairs)
    return sum((cells[(r, c)] - rows[r] * cols[c] / n) ** 2 / (rows[r] * cols[c] / n)
               for r in rows for c in cols)


def association(ecotype_of, category_of, seed=42):
    """(n, ecotypes, categories, chi2, Cramer's V, permutation p) over genomes with a known category."""
    pairs = [(ecotype_of[a], category_of[a]) for a in ecotype_of if category_of.get(a, UNKNOWN) != UNKNOWN]
    n_rows = len({r for r, _ in pairs})
    n_cols = len({c for _, c in pairs})
    if n_rows < 2 or n_cols < 2:
        return len(pairs), n_rows, n_cols, None, None, None
    observed = chi2(pairs)
    # Permutation p-value: correct for the sparse tables small ecotypes give, where the chi2 approximation is not.
    rng = random.Random(seed)
    eco, cat = [r for r, _ in pairs], [c for _, c in pairs]
    hits = 0
    for _ in range(PERMUTATIONS):
        rng.shuffle(cat)
        hits += chi2(list(zip(eco, cat))) >= observed - 1e-9
    v = sqrt(observed / (len(pairs) * (min(n_rows, n_cols) - 1)))
    return len(pairs), n_rows, n_cols, round(observed, 3), round(v, 3), round((hits + 1) / (PERMUTATIONS + 1), 4)


def write_crosstab(path, ecotype_of, category_of):
    table = Counter((ecotype_of[a], category_of.get(a, UNKNOWN)) for a in ecotype_of)
    categories = sorted({c for _, c in table})
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ecotype"] + categories)
        for eco in sorted({e for e, _ in table}, key=lambda e: (len(e), e)):
            writer.writerow([eco] + [table[(eco, c)] for c in categories])


def run_niche_analysis(ecosim_out_dir, species, labeler=label_with_claude, fetcher=fetch_biosamples):
    """Run the whole step on one EcoSim output folder. Returns the association CSV path, or None."""
    out = Path(ecosim_out_dir)
    membership_csv = out / "ecotype_membership.csv"
    if not membership_csv.exists():
        print(f"No {membership_csv.name} in {out}; skipping niche analysis.")
        return None
    xml_used, ecotype_of = load_membership(membership_csv)
    ecotype_of = {a: e for a, e in ecotype_of.items() if a.startswith(("GCA_", "GCF_"))}
    print(f"Niche analysis on {xml_used}: {len(ecotype_of)} genomes")

    meta = fetcher(sorted(ecotype_of))
    keys = sorted({_llm_key(meta.get(a, {})) for a in ecotype_of})
    print(f"Asking {LLM_MODEL} to label {len(keys)} distinct metadata record(s)")
    labels = labeler(species, keys)

    columns = {"host": {}, "niche": {}, "country": {}}
    with open(out / "genome_metadata.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["accession", "ecotype", "host", "niche", "country", "year", *(f"raw_{k}" for k in LLM_FIELDS)])
        for acc in sorted(ecotype_of):
            attrs = meta.get(acc, {})
            host, niche = labels.get(_llm_key(attrs), (UNKNOWN, UNKNOWN))
            country, year = country_and_year(attrs)
            columns["host"][acc], columns["niche"][acc], columns["country"][acc] = host, niche, country
            writer.writerow([acc, ecotype_of[acc], host, niche, country, year, *(attrs.get(k, "") for k in LLM_FIELDS)])

    assoc_path = out / "ecotype_niche_association.csv"
    with open(assoc_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["variable", "ecotype_file", "genomes_with_data", "ecotypes", "categories",
                         "chi2", "cramers_v", "permutation_p"])
        for variable, category_of in columns.items():
            write_crosstab(out / f"ecotype_x_{variable}.csv", ecotype_of, category_of)
            stats = association(ecotype_of, category_of)
            writer.writerow([variable, xml_used, *stats])
            print(f"  {variable}: n={stats[0]}, V={stats[4]}, p={stats[5]}")
    print(f"Saved {assoc_path}")
    return assoc_path


def _self_test():
    # Perfect association: ecotype 1 all blood, ecotype 2 all skin.
    eco = {f"GCA_{i}": "1" if i < 10 else "2" for i in range(20)}
    niche = {a: "blood" if e == "1" else "skin" for a, e in eco.items()}
    n, r, c, stat, v, p = association(eco, niche)
    assert (n, r, c, stat, v) == (20, 2, 2, 20.0, 1.0) and p < 0.01, (n, r, c, stat, v, p)
    # Unknowns are dropped; one category left means no test.
    assert association(eco, {a: "blood" if e == "1" else UNKNOWN for a, e in eco.items()})[3] is None
    # No association: each ecotype half blood, half skin.
    _, _, _, stat, v, p = association(eco, {a: "blood" if i % 2 else "skin" for i, a in enumerate(eco)})
    assert stat == 0 and v == 0 and p == 1
    assert country_and_year({"geo_loc_name": "United Kingdom: Scotland", "collection_date": "2018"}) == ("United Kingdom", "2018")
    assert country_and_year({"geo_loc_name": "missing", "collection_date": "not collected"}) == (UNKNOWN, "")
    print("self-test passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ecosim_out_dir", nargs="?")
    parser.add_argument("species", nargs="?")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        _self_test()
    elif args.ecosim_out_dir and args.species:
        run_niche_analysis(args.ecosim_out_dir, args.species)
    else:
        parser.error("give ecosim_out_dir and species, or --self-test")
