"""
NCBI access for the isolate metadata catalogue, through the public Datasets v2 API (no key needed):

  harvest()       every assembly report for a taxon at the chosen assembly levels, paged and resumable,
                  kept as compact JSON lines (accession, organism, assembly facts, BioSample attributes)
  resolve_taxa()  NCBI Taxonomy records (scientific name, rank, lineage) for names or taxids, cached

Requests carry a generic User-Agent and nothing about who runs them.
"""
import datetime
import gzip
import http.client
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
USER_AGENT = "isolate-catalogue/0.1"
LEVELS = ("complete_genome", "chromosome")
PAGE_SIZE = 1000
TAXA_BATCH = 200  # names or taxids per taxonomy request


def fetch_json(url, payload=None, retries=6, timeout=300, sleep=time.sleep):
    """GET (or POST payload as JSON) and parse the reply, retrying throttling, server and network errors."""
    data = json.dumps(payload).encode() if payload is not None else None
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    if data is not None:
        headers["Content-Type"] = "application/json"
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.load(response)
        except urllib.error.HTTPError as e:
            if e.code not in (429, 500, 502, 503, 504) or attempt == retries - 1:
                raise
        except (urllib.error.URLError, OSError, http.client.HTTPException, json.JSONDecodeError):
            if attempt == retries - 1:
                raise
        sleep(min(2 ** attempt * 2, 120))
    raise RuntimeError("unreachable")


def report_url(taxon, levels=LEVELS, source="genbank", page_size=PAGE_SIZE, page_token=None, extra=()):
    params = [("filters.assembly_source", source), ("page_size", page_size)]
    params += [("filters.assembly_level", level) for level in levels]
    params += list(extra)
    if page_token:
        params.append(("page_token", page_token))
    return f"{API}/genome/taxon/{taxon}/dataset_report?{urllib.parse.urlencode(params)}"


def compact(report):
    """The fields the catalogue keeps from one Datasets assembly report."""
    organism = report.get("organism") or {}
    names = organism.get("infraspecific_names") or {}
    info = report.get("assembly_info") or {}
    biosample = info.get("biosample") or {}
    stats = report.get("assembly_stats") or {}
    checkm = report.get("checkm_info") or {}
    ani = report.get("average_nucleotide_identity") or {}
    atypical = info.get("atypical") or {}
    attributes = {}
    for attribute in biosample.get("attributes") or []:
        name, value = attribute.get("name"), attribute.get("value", "")
        if name and value and name not in attributes:
            attributes[name] = value
    return {
        "accession": report.get("accession", ""),
        "paired_accession": report.get("paired_accession", ""),
        "tax_id": organism.get("tax_id"),
        "organism_name": organism.get("organism_name", ""),
        "strain": names.get("strain") or biosample.get("strain", ""),
        "isolate": names.get("isolate") or biosample.get("isolate", ""),
        "serovar": biosample.get("serovar", ""),
        "assembly_level": info.get("assembly_level", ""),
        "assembly_status": info.get("assembly_status", ""),
        "release_date": info.get("release_date", ""),
        "submitter": info.get("submitter", ""),
        "bioproject": info.get("bioproject_accession", ""),
        "refseq_category": info.get("refseq_category", ""),
        "atypical": bool(atypical.get("is_atypical")),
        "atypical_warnings": "; ".join(atypical.get("warnings") or []),
        "type_material": bool(report.get("type_material")),
        "biosample": biosample.get("accession", ""),
        "biosample_package": biosample.get("package", ""),
        "biosample_title": (biosample.get("description") or {}).get("title", ""),
        "checkm_completeness": checkm.get("completeness"),
        "checkm_contamination": checkm.get("contamination"),
        "taxonomy_check": ani.get("taxonomy_check_status", ""),
        "total_length": stats.get("total_sequence_length"),
        "gc_percent": stats.get("gc_percent"),
        "attributes": attributes,
    }


def _write_json(path, data):
    tmp = Path(f"{path}.tmp{os.getpid()}")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=1, sort_keys=True)
    os.replace(tmp, path)


def _read_json(path, default):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return default


def _now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def harvest(out_dir, taxon=2, levels=LEVELS, source="genbank", page_size=PAGE_SIZE, fetch=fetch_json,
            log=print, restart=False):
    """Download every assembly report, one gzipped JSON-lines file per page, then the list of
    metagenome-assembled genomes (MAGs) among them. Resumes from harvest_state.json after an interruption;
    restart=True starts over (needed if NCBI expired the saved page token)."""
    out = Path(out_dir)
    pages = out / "pages"
    pages.mkdir(parents=True, exist_ok=True)
    state_path = out / "harvest_state.json"
    query = {"taxon": str(taxon), "levels": list(levels), "source": source}
    state = {} if restart else _read_json(state_path, {})
    if state.get("query") == query and state.get("done"):
        log(f"Harvest already complete: {state['records']:,} assemblies in {state['pages']} pages ({state['finished']})")
        return state
    if state.get("query") != query:
        for old in pages.glob("page_*.jsonl.gz"):
            old.unlink()
        state = {"query": query, "pages": 0, "records": 0, "token": None, "started": _now(), "api": API}
        _write_json(state_path, state)
    while True:
        data = fetch(report_url(taxon, levels, source, page_size, state["token"]))
        reports = data.get("reports") or []
        number = state["pages"] + 1
        page = pages / f"page_{number:05d}.jsonl.gz"
        tmp = Path(f"{page}.tmp")
        with gzip.open(tmp, "wt", encoding="utf-8") as f:
            for report in reports:
                f.write(json.dumps(compact(report), ensure_ascii=False, separators=(",", ":")) + "\n")
        os.replace(tmp, page)
        state.update(pages=number, records=state["records"] + len(reports), token=data.get("next_page_token"),
                     total_count=data.get("total_count", state.get("total_count")))
        _write_json(state_path, state)
        if number % 10 == 0 or not state["token"]:
            log(f"  page {number}: {state['records']:,} of {state.get('total_count') or '?':,} assemblies")
        if not state["token"] or not reports:
            break
    mags = fetch_accessions(taxon, levels, source, fetch=fetch,
                            extra=[("filters.is_metagenome_derived", "metagenome_derived_only")])
    (out / "mag_accessions.txt").write_text("".join(f"{a}\n" for a in sorted(mags)), encoding="utf-8")
    state.update(done=True, finished=_now(), mags=len(mags), token=None)
    _write_json(state_path, state)
    if state.get("total_count") and state["records"] != state["total_count"]:
        log(f"  note: NCBI reported {state['total_count']:,} assemblies but {state['records']:,} were returned")
    log(f"Harvested {state['records']:,} assemblies ({len(mags):,} metagenome-assembled) into {pages}")
    return state


def fetch_accessions(taxon, levels=LEVELS, source="genbank", extra=(), fetch=fetch_json):
    """Accessions only, for a filtered query (e.g. the MAGs)."""
    accessions, token = [], None
    while True:
        data = fetch(report_url(taxon, levels, source, PAGE_SIZE, token,
                                extra=[*extra, ("returned_content", "ASSM_ACC")]))
        accessions += [r["accession"] for r in data.get("reports") or [] if r.get("accession")]
        token = data.get("next_page_token")
        if not token or not data.get("reports"):
            return accessions


def iter_assemblies(out_dir):
    """Every harvested assembly (compact dicts), page by page."""
    for page in sorted((Path(out_dir) / "pages").glob("page_*.jsonl.gz")):
        with gzip.open(page, "rt", encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    yield json.loads(line)


def mag_accessions(out_dir):
    path = Path(out_dir) / "mag_accessions.txt"
    return set(path.read_text(encoding="utf-8").split()) if path.exists() else set()


def _slim(taxonomy):
    classification = {rank: {"id": v.get("id"), "name": v.get("name", "")}
                      for rank, v in (taxonomy.get("classification") or {}).items() if isinstance(v, dict)}
    return {"tax_id": taxonomy.get("tax_id"),
            "name": (taxonomy.get("current_scientific_name") or {}).get("name", ""),
            "rank": (taxonomy.get("rank") or "").lower(),
            "parents": taxonomy.get("parents") or [],
            "classification": classification}


def resolve_taxa(queries, cache_path, fetch=fetch_json, batch=TAXA_BATCH, log=print):
    """{query: record} for NCBI Taxonomy names (scientific or common, any case) or taxids. A record is
    {"tax_id", "name", "rank", "parents", "classification"}; None when NCBI does not recognise the query,
    and {"ambiguous": [records]} when one name matches several taxa. Answers are cached in cache_path."""
    cache_path = Path(cache_path)
    cache = _read_json(cache_path, {})
    todo = sorted({str(q) for q in queries} - set(cache))
    for start in range(0, len(todo), batch):
        chunk = todo[start:start + batch]
        data = fetch(f"{API}/taxonomy/dataset_report", payload={"taxons": chunk})
        found = {}
        for report in data.get("reports") or []:
            query = (report.get("query") or [None])[0]
            if query is not None and report.get("taxonomy"):
                found.setdefault(str(query), []).append(_slim(report["taxonomy"]))
        for query in chunk:
            records = {r["tax_id"]: r for r in found.get(query, [])}
            if not records:
                cache[query] = None
            elif len(records) == 1:
                cache[query] = next(iter(records.values()))
            else:
                cache[query] = {"ambiguous": list(records.values())}
        if (start // batch) % 20 == 19 or start + batch >= len(todo):
            _write_json(cache_path, cache)
            log(f"  taxonomy: {min(start + batch, len(todo)):,} of {len(todo):,} new queries resolved")
    return {str(q): cache.get(str(q)) for q in queries}
