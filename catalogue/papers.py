"""
Pilot: what the papers linked to the catalogue's genomes add to their BioSample metadata.

NCBI links many genomes' BioProjects to PubMed papers. `papers-pilot` reads a sample of those papers with Claude
Fable 5.2 and measures how many genomes gain something the catalogue lacks. A gain is one of:
  - a host or source where the catalogue says unknown;
  - a country or collection year;
  - facts BioSample fields rarely carry: disease, outbreak, isolation method or phenotype.

The pilot has four steps:

  sample   Genomes are taken in random order. Each one whose BioProject links a paper adds one of its papers,
           until there are `--sample` distinct papers. Because the sample is over genomes, its gain rate is
           a rate over the genomes that have a paper.
  texts    The PubMed Central full text where NCBI serves it (methods-like sections first, then tables, the
           other sections and figure captions); otherwise the PubMed title and abstract. Cached in
           <workdir>/papers/.
  read     One Batch API request per paper, holding the text, the genomes NCBI links to it, and the curation
           vocabulary. The JSON schema limits hosts and sources to the vocabulary. Every value must quote the
           paper, and a value whose quote is not in the text is dropped.
  compare  Each linked genome's answer (its own, else the paper's study-wide one) against its catalogue row.

Links to PubMed come from E-utilities (free) and are cached in <workdir>/papers/links.json. The requests go
out in waves, each checked beforehand against the budget at its worst case (every request using all of
max_tokens), so the budget holds even though a batch cannot be stopped part-way.
"""
import collections
import csv
import gzip
import json
import random
import re
import unicodedata
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

try:
    from steps.metadata_curation import BATCH_DISCOUNT, PRICES, usage_cost
    from steps.papers import _eutils, fetch_bioproject_papers, fetch_pubmed
except ImportError:  # pragma: no cover - run from another directory
    from metadata_curation import BATCH_DISCOUNT, PRICES, usage_cost
    from papers import _eutils, fetch_bioproject_papers, fetch_pubmed

MODEL, EFFORT = "claude-fable-5-2", "high"
MAX_TOKENS = 12_000   # per paper; also what the budget's worst case assumes
RETRY_MAX_TOKENS = 20_000  # for papers sent again with --retry-model (one ran out of room at 12,000)
MAX_CHARS = 100_000   # paper text per request (about 25k tokens)
MAX_GENOMES = 30      # linked genomes named per paper; the sampled ones are always among them
WAVE = 100            # papers per batch
FACTS = ("host", "source", "country", "collection_year", "disease", "isolation_method", "outbreak", "phenotype")
CORE = ("host", "source", "country", "collection_year")  # the catalogue's own fields
IDCONV = "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/"
USER_AGENT = "isolate-catalogue/0.1"
METHODS = re.compile(r"method|material|sampl|isolat|strain|bacteri|collect|study design|population", re.I)


# --- links and sample ---------------------------------------------------------------------------------
def _save(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, sort_keys=True))
    tmp.replace(path)


def ensure_links(projects, path, fetch=fetch_bioproject_papers):
    """{bioproject: [pmid, ...]}, looked up where not cached (a project NCBI does not know, or that links no
    paper, is cached as [])."""
    path = Path(path)
    links = json.loads(path.read_text()) if path.exists() else {}
    todo = sorted(p for p in projects if p and p not in links)
    for start in range(0, len(todo), 1000):
        chunk = todo[start:start + 1000]
        found = fetch(chunk)
        for p in chunk:
            links[p] = sorted(set(found.get(p, ("", []))[1]))
        _save(path, links)
    return links


def sample_papers(genome_project, links_path, n, seed=1, fetch=fetch_bioproject_papers):
    """Genomes in random order; each whose BioProject links a paper adds one of its papers, until n distinct
    papers. Returns ({pmid: [sampled accessions]}, links, genomes tried, genomes tried that have a paper)."""
    rng = random.Random(seed)
    order = sorted(genome_project)
    rng.shuffle(order)
    chosen, tried, with_paper, links = {}, 0, 0, {}
    for start in range(0, len(order), 500):
        block = order[start:start + 500]
        links = ensure_links({genome_project[g] for g in block}, links_path, fetch)
        for g in block:
            tried += 1
            pmids = links.get(genome_project[g]) or []
            if pmids:
                with_paper += 1
                chosen.setdefault(rng.choice(pmids), []).append(g)
                if len(chosen) >= n:
                    return chosen, links, tried, with_paper
    return chosen, links, tried, with_paper


# --- texts ---------------------------------------------------------------------------------------------
def _get_json(url):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response)


def pmc_ids(pmids):
    """{pmid: PMCID} for the papers in PubMed Central."""
    out, pmids = {}, sorted(pmids)
    for start in range(0, len(pmids), 150):
        query = urllib.parse.urlencode({"ids": ",".join(pmids[start:start + 150]), "format": "json",
                                        "tool": "isolate-catalogue"})
        for record in _get_json(f"{IDCONV}?{query}").get("records", []):
            if record.get("pmcid") and record.get("pmid"):
                out[str(record["pmid"])] = record["pmcid"]
    return out


def _text(element):
    return re.sub(r"\s+", " ", " ".join(element.itertext())).strip() if element is not None else ""


def jats_text(xml, methods_only=False):
    """Plain text of a PMC article (JATS XML): title, abstract, methods-like sections, tables, the other
    sections, figure captions. None without a body: outside its open-access subset PMC serves only the front
    matter. methods_only: only the methods-like sections and the tables, where a paper says where its isolates
    came from (None when it has neither)."""
    root = ET.fromstring(xml)
    article = root if root.tag == "article" else root.find(".//article")
    body = article.find("body") if article is not None else None
    if body is None or not _text(body):
        return None
    tables = []
    for wrap in article.iter("table-wrap"):
        rows = [" | ".join(_text(cell) for cell in row) for row in wrap.iter("tr")]
        tables.append("Table. " + _text(wrap.find("caption")) + "\n" + "\n".join(rows))
    figures = ["Figure. " + _text(fig.find("caption")) for fig in article.iter("fig") if fig.find("caption") is not None]
    for parent in list(article.iter()):  # tables and figures are listed on their own
        for child in list(parent):
            if child.tag in ("table-wrap", "fig"):
                parent.remove(child)
    sections = [(sec.findtext("title") or "", _text(sec)) for sec in body.findall("sec")] or [("", _text(body))]
    first = [text for title, text in sections if METHODS.search(title)]
    rest = [text for title, text in sections if not METHODS.search(title)]
    if methods_only:
        return "\n\n".join(p for p in (*first, *tables) if p.strip()) or None
    title = _text(article.find(".//article-meta//article-title"))
    abstract = " ".join(_text(a) for a in article.findall(".//article-meta/abstract"))
    parts = [f"Title: {title}", f"Abstract: {abstract}", *first, *tables, *rest, *figures]
    return "\n\n".join(p for p in parts if p.strip())


def paper_texts(pmids, folder, pmc_lookup=pmc_ids, fetch_article=None, fetch_abstracts=fetch_pubmed):
    """{pmid: {"kind": "full text" | "abstract" | "none", "text"}}, cached one file per paper in folder."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    fetch_article = fetch_article or (lambda pmcid: _eutils("efetch", db="pmc", id=pmcid.replace("PMC", ""),
                                                            retmode="xml"))
    out = {p: json.loads((folder / f"{p}.json").read_text()) for p in pmids if (folder / f"{p}.json").exists()}
    todo = sorted(set(pmids) - set(out))
    if not todo:
        return out
    pmcs = pmc_lookup(todo)
    for pmid in todo:
        if pmid in pmcs:
            try:
                text = jats_text(fetch_article(pmcs[pmid]))
            except ET.ParseError:
                text = None
            if text:
                out[pmid] = {"kind": "full text", "pmcid": pmcs[pmid], "text": text[:MAX_CHARS]}
    abstracts = fetch_abstracts([p for p in todo if p not in out])
    for pmid in todo:
        if pmid not in out:
            a = abstracts.get(pmid) or {}
            text = f"Title: {a.get('title', '')}\n\nAbstract: {a.get('abstract', '')}".strip()
            out[pmid] = {"kind": "abstract" if a.get("abstract") else "none", "text": text}
        (folder / f"{pmid}.json").write_text(json.dumps(out[pmid]))
    return out


# --- requests ------------------------------------------------------------------------------------------
FIELD_HELP = {
    "host": "- host: the host the isolate was taken from, as a key from the host list below.\n",
    "source": "- source: the material or body site sampled, as a key from the source list below.\n",
    "country": "- country: the country of collection, as an INSDC country name (USA, United Kingdom, China, ...).\n",
    "collection_year": "- collection_year: the year the isolate was collected, not published; a range such as "
                       "2010-2012 is fine.\n",
    "disease": "- disease: the disease or clinical condition of the host the isolate came from.\n",
    "isolation_method": "- isolation_method: how the isolate was obtained, such as an enrichment or a selective "
                        "medium.\n",
    "outbreak": "- outbreak: the outbreak, investigation or named strain collection the isolate belongs to.\n",
    "phenotype": "- phenotype: traits the paper reports for the isolate, such as antimicrobial resistance or "
                 "virulence.\n"}


def system_prompt(vocab, facts=FACTS):
    return (
        "You read a published paper about bacteria and report what it says about where the isolates behind the "
        "listed genomes came from. NCBI links these genomes to the paper through their BioProject. The paper may "
        "describe them one by one (by strain or isolate name, in the text or a table), all together (\"all "
        "isolates came from raw milk\"), or not at all.\n\n"
        "Report:\n"
        "- study: what the paper states for all of its isolates together. Leave a field empty (unknown for host "
        "and source) unless the paper states it for all of them.\n"
        "- genomes: one entry for each listed genome the paper identifies individually, matched by strain or "
        "isolate name or by accession, with what the paper states about that isolate. Leave out the others.\n"
        "- describes: whether the paper describes the listed genomes individually, collectively, both, or not at "
        "all.\n\n"
        "Fields:\n" + "".join(FIELD_HELP[f] for f in facts) +
        "Each field has a quote: the exact words from the paper that support it, copied verbatim (at most about "
        "20 words; quotes are checked against the text). Leave both empty when the paper does not say. Answer "
        "from the paper only: what the species usually does, or what a strain name suggests, is not evidence.\n\n"
        + vocab.prompt_listing())


def answer_schema(vocab, accessions, facts=FACTS):
    fact = {}
    for f in facts:
        enum = {"host": sorted(vocab.hosts), "source": sorted(vocab.terms)}.get(f)
        fact[f] = {"type": "string", "enum": enum} if enum else {"type": "string"}
        fact[f"{f}_quote"] = {"type": "string"}
    names = list(fact)
    item = {"type": "object", "additionalProperties": False, "required": ["accession", "name_in_paper", *names],
            "properties": {"accession": {"type": "string", "enum": sorted(accessions)},
                           "name_in_paper": {"type": "string"}, **fact}}
    return {"type": "object", "additionalProperties": False, "required": ["describes", "study", "genomes", "certain"],
            "properties": {"describes": {"type": "string",
                                         "enum": ["individually", "collectively", "both", "not_at_all"]},
                           "study": {"type": "object", "additionalProperties": False, "required": names,
                                     "properties": fact},
                           "genomes": {"type": "array", "items": item},
                           "certain": {"type": "boolean"}}}


def build_request(vocab, pmid, paper, genomes, more=0, facts=FACTS, max_tokens=MAX_TOKENS):
    """Messages API parameters for one paper; genomes: [{"accession", "organism", "strain", "isolate"}].
    facts: the facts to ask for (NARROW_FACTS for the narrow request)."""
    lines = [f"- {g['accession']}: {g['organism']}; strain {g['strain'] or '-'}; isolate {g['isolate'] or '-'}"
             for g in genomes]
    user = (f"Paper (PMID {pmid}, {paper['kind']}):\n{paper['text']}\n\n"
            "Genomes NCBI links to this paper (accession: organism; strain; isolate):\n" + "\n".join(lines)
            + (f"\n(and {more:,} more genomes from the same BioProject, not listed)" if more else ""))
    schema = answer_schema(vocab, [g["accession"] for g in genomes], facts)
    return {"model": MODEL, "max_tokens": max_tokens,
            "system": [{"type": "text", "text": system_prompt(vocab, facts),
                        "cache_control": {"type": "ephemeral", "ttl": "1h"}}],
            "messages": [{"role": "user", "content": user}],
            "output_config": {"effort": EFFORT, "format": {"type": "json_schema", "schema": schema}}}


# The narrow request, for papers both models refused (pathogen and toxin papers): only the catalogue's own
# facts, from the methods sections and tables, without the results and discussion about toxins and virulence.
NARROW_FACTS = CORE
NARROW_MAX_TOKENS = 6_000
NARROW_WAVE = 12  # small waves: each is checked against the budget at its worst case after the last one's spend


def narrow_text(paper, fetch_article=None):
    """{"kind", "text"} for the narrow request: a full text's methods-like sections and tables (fetched again
    from PMC), an abstract as it is; None when a full text has neither."""
    if paper["kind"] != "full text":
        return paper
    fetch_article = fetch_article or (lambda pmcid: _eutils("efetch", db="pmc", id=pmcid.replace("PMC", ""),
                                                            retmode="xml"))
    text = jats_text(fetch_article(paper["pmcid"]), methods_only=True)
    return {"kind": "methods and tables", "pmcid": paper["pmcid"], "text": text[:MAX_CHARS]} if text else None


def worst_case_usd(requests):
    """Upper bound for a wave at Batch prices: input at 3 characters a token, every request using max_tokens,
    each at its own model's prices."""
    usd = 0.0
    for r in requests:
        price_in, price_out = PRICES[r["params"]["model"]]
        usd += (len(json.dumps(r["params"])) / 3.0 * price_in + r["params"]["max_tokens"] * price_out) / 1e6
    return usd * BATCH_DISCOUNT


# --- answers and comparison ------------------------------------------------------------------------------
def _norm(s):
    s = unicodedata.normalize("NFKC", s or "").lower()
    s = re.sub(r"[‐-―−]", "-", s)
    s = re.sub(r"[‘’‚‛′`]", "'", s)
    s = re.sub(r"[“”„″]", '"', s)
    return re.sub(r"\s+", " ", s).strip()


COUNTRY_ALIASES = {"czech republic": "czechia", "united states": "usa", "united states of america": "usa",
                   "us": "usa", "uk": "united kingdom", "great britain": "united kingdom", "england": "united kingdom",
                   "scotland": "united kingdom", "wales": "united kingdom", "korea": "south korea",
                   "republic of korea": "south korea", "the netherlands": "netherlands", "holland": "netherlands",
                   "russian federation": "russia", "viet nam": "viet nam", "vietnam": "viet nam",
                   "turkiye": "turkey", "türkiye": "turkey", "ivory coast": "cote d'ivoire"}


def _country(name):
    name = _norm(name).strip(" .")
    return COUNTRY_ALIASES.get(name, name)


def _years(text):
    """(first, last) year named in a value ("2015", "2010-2012", "from 1930 to 1964"), or None."""
    years = [int(y) for y in re.findall(r"\b(1[89]\d\d|20\d\d)\b", text or "")]
    return (min(years), max(years)) if years else None


def checked(answer, text):
    """({fact: (value, quote)} for the facts whose value is stated and whose quote is in the paper text, how many
    stated facts failed the quote check)."""
    text, out, failed = _norm(text), {}, 0
    for f in FACTS:
        value, quote = (answer.get(f) or "").strip(), (answer.get(f"{f}_quote") or "").strip()
        if not value or value == "unknown":
            continue
        if quote and _norm(quote) in text:
            out[f] = (value, quote)
        else:
            failed += 1
    return out, failed


def parse_message(message):
    """The answer dict from a finished message, or None with the reason."""
    if getattr(message, "stop_reason", None) != "end_turn":
        return None, str(getattr(message, "stop_reason", "no message"))
    try:
        return json.loads(next(b.text for b in message.content if b.type == "text")), "ok"
    except (StopIteration, ValueError, TypeError):
        return None, "unparsable"


def compare(row, facts, vocab):
    """Per fact: none (the paper does not say), new (the catalogue lacks it), same, same niche (source), or
    differs; disease is "in BioSample" when the BioSample already names one. facts: {fact: value or (value,
    quote)}."""
    out = {}
    for f in FACTS:
        value = facts.get(f)
        value = value[0] if isinstance(value, tuple) else value
        if value is None:
            out[f] = "none"
            continue
        if f in ("host", "source"):
            table = vocab.hosts if f == "host" else vocab.terms
            label = table[value]["label"]
            ours = row["host" if f == "host" else "source_term"]
            if ours in ("", "unknown"):
                out[f] = "new"
            elif ours == label:
                out[f] = "same"
            elif f == "source" and row.get("niche") == vocab.niches.get(vocab.terms[value]["niche"]):
                out[f] = "same niche"
            else:
                out[f] = "differs"
        elif f == "country":
            ours = row.get("country", "")
            same = _country(ours) == _country(value)
            out[f] = "new" if ours in ("", "unknown") else "same" if same else "differs"
        elif f == "collection_year":
            ours, span = row.get("collection_year", ""), _years(value)
            if not ours:
                out[f] = "new"
            elif span and span[0] == span[1] == int(ours):
                out[f] = "same"
            elif span and span[0] <= int(ours) <= span[1]:
                out[f] = "within range"  # a study period that includes the catalogue's year
            else:
                out[f] = "differs"
        elif f == "disease":
            out[f] = "in BioSample" if row.get("raw_host_disease") else "new"
        else:
            out[f] = "new"
    return out


def facts_for(answer, accession, text):
    """The facts for one genome, {fact: (value, quote)}: its own entry where the paper identifies it, over the
    study-wide ones. Returns (facts, "genome" or "study", quotes not found)."""
    study, failed = checked(answer.get("study") or {}, text)
    own = next((g for g in answer.get("genomes") or [] if g.get("accession") == accession), None)
    if own is None:
        return study, "study", failed
    mine, failed_own = checked(own, text)
    return {**study, **mine}, "genome", failed + failed_own


def summarize(rows, weight_key="sampled"):
    """Shares over the rows of genomes: any gain, a gain in the catalogue's own fields, in the extra fields,
    per fact, and disagreements with the catalogue."""
    rows = [r for r in rows if r[weight_key]] if weight_key else rows
    n = len(rows)
    if not n:
        return {"genomes": 0}
    def share(pred):
        k = sum(1 for r in rows if pred(r))
        return {"genomes": k, "share": round(k / n, 3)}
    out = {"genomes": n,
           "any gain": share(lambda r: any(r[f] == "new" for f in FACTS)),
           "gain in host, source, country or year": share(lambda r: any(r[f] == "new" for f in CORE)),
           "gain in disease, method, outbreak or phenotype": share(lambda r: any(r[f] == "new" for f in FACTS[4:])),
           "differs from the catalogue (host, source, country or year)": share(
               lambda r: any(r[f] == "differs" for f in CORE)),
           "facts from the genome's own entry": share(lambda r: r["origin"] == "genome")}
    out["by fact"] = {f: dict(collections.Counter(r[f] for r in rows)) for f in FACTS}
    return out


def catalogue_rows(release_tsv, accessions):
    """The catalogue rows for these accessions, streamed from the release table."""
    want, rows = set(accessions), {}
    with gzip.open(release_tsv, "rt", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["accession"] in want:
                rows[row["accession"]] = row
    return rows
