#!/usr/bin/env python
"""Run the pipeline from an NCBI link: python run_from_link.py <url|taxid|accession|name> [pipeline.py args]"""
import json
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse


def parse_input(text):
    """Return ('accession'|'taxon', value) from an NCBI URL, bare taxid/accession, or species name."""
    text = text.strip()
    m = re.search(r"GC[AF]_\d+(?:\.\d+)?", text)
    if m:
        return "accession", m.group(0)
    if not re.match(r"https?://|(www\.)?ncbi\.nlm\.nih\.gov", text):
        return "taxon", text  # bare taxid or species name
    url = urlparse(text if "//" in text else "https://" + text)
    q = {k: v[0] for k, v in parse_qs(url.query).items()}
    m = re.search(r"/taxonomy/(\d+)", url.path)
    for value in (m and m.group(1), q.get("id"), q.get("taxon"), q.get("name"), q.get("term")):
        if value:
            return "taxon", unquote(value).strip('"')
    sys.exit(f"Could not find a taxid, accession or name in: {text}")


def ncbi_lookup(*args):
    """Run `datasets summary ...` and return the first JSON record."""
    cmd = ["datasets", "summary", *args, "--as-json-lines"]
    try:
        out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
        return json.loads(out.splitlines()[0])
    except (OSError, subprocess.CalledProcessError, IndexError, json.JSONDecodeError) as e:
        sys.exit(f"NCBI lookup failed ({' '.join(cmd[:5])}): {e}. Is `datasets` on PATH (conda activate ecotype)?")


def resolve_species(text, lookup=ncbi_lookup):
    kind, value = parse_input(text)
    if kind == "accession":
        value = str(lookup("genome", "accession", value)["organism"]["tax_id"])
    tax = lookup("taxonomy", "taxon", value)["taxonomy"]
    species = tax.get("classification", {}).get("species", {}).get("name")
    if not species:
        sys.exit(f"'{text}' resolves to {tax.get('current_scientific_name', {}).get('name', value)} "
                 f"(rank {tax.get('rank', '?').lower()}), which is not a species or below. Link a species page.")
    return species


def _self_test():
    sp = {"classification": {"species": {"name": "Streptococcus pyogenes"}}, "rank": "SPECIES"}
    genus = {"classification": {"genus": {"name": "Streptococcus"}}, "rank": "GENUS",
             "current_scientific_name": {"name": "Streptococcus"}}
    seen = []

    def fake(*args):
        seen.append(args)
        if args[0] == "genome":
            return {"organism": {"tax_id": 160490}}
        return {"taxonomy": genus if args[2] == "1301" else sp}

    base = "https://www.ncbi.nlm.nih.gov/"
    for text, want in [
        (base + "Taxonomy/Browser/wwwtax.cgi?id=1314", ("taxon", "1314")),
        (base + "Taxonomy/Browser/wwwtax.cgi?mode=Info&id=1314&lvl=3", ("taxon", "1314")),
        (base + "Taxonomy/Browser/wwwtax.cgi?name=Streptococcus+pyogenes", ("taxon", "Streptococcus pyogenes")),
        (base + "datasets/taxonomy/1314/", ("taxon", "1314")),
        (base + "datasets/genome/GCF_000006785.2/", ("accession", "GCF_000006785.2")),
        (base + "datasets/genome/?taxon=1314", ("taxon", "1314")),
        (base + "genome/?term=Streptococcus%20pyogenes", ("taxon", "Streptococcus pyogenes")),
        ("1314", ("taxon", "1314")),
        ("GCF_000006785.2", ("accession", "GCF_000006785.2")),
        ("Streptococcus pyogenes", ("taxon", "Streptococcus pyogenes")),
    ]:
        assert parse_input(text) == want, (text, parse_input(text))
    assert resolve_species("GCF_000006785.2", fake) == "Streptococcus pyogenes"
    assert seen[-1] == ("taxonomy", "taxon", "160490")
    try:
        resolve_species(base + "datasets/taxonomy/1301/", fake)
        raise AssertionError("genus should be rejected")
    except SystemExit as e:
        assert "not a species" in str(e)
    print("self-test passed")


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == "--self-test":
        _self_test() if len(sys.argv) > 1 else sys.exit(__doc__)
        sys.exit()
    species = resolve_species(sys.argv[1])
    print(f"Resolved species: {species}", flush=True)
    sys.exit(subprocess.call([sys.executable, str(Path(__file__).with_name("pipeline.py")),
                              "--species", species, *sys.argv[2:]]))
