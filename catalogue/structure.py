"""
Pilot: each species' published population structure, as a yardstick for its ecotypes.

For each species, `structure-pilot` does four things:

  papers   PubMed's most relevant papers on its population structure (reviews first), as their abstracts, plus
           the PubMed Central full text of up to FULL_TEXTS of them where NCBI serves it.
  read     One Batch API request per species. Claude Fable 5.2 lists the groups the papers define (subspecies,
           lineages, phylogroups, clades, clonal complexes, serotypes, ecotypes), how they nest, what defines
           them, the reference strains the papers name in them, and the ecology the papers report. Every group
           needs a verbatim quote from a listed paper; quotes not found in that paper are dropped.
  match    The named reference strains are matched to the catalogue's genomes of that species by strain or
           isolate name. That gives genomes with a known, published group: a yardstick for EcoSim, TreeCluster,
           fastbaps and PopPUNK.
  report   reports/structure_pilot_<time>.json, with *_groups.tsv and *_genomes.tsv.

Paper lists and texts are cached in <workdir>/papers/ (shared with the papers pilot), answers in
<workdir>/structure/answers/.
"""
import collections
import csv
import gzip
import json
import re
from pathlib import Path

try:
    from catalogue import papers
    from steps.metadata_curation import BATCH_DISCOUNT, PRICES
    from steps.papers import _eutils
except ImportError:  # pragma: no cover
    import papers
    from metadata_curation import BATCH_DISCOUNT, PRICES
    from papers import _eutils

MODEL, EFFORT, MAX_TOKENS = papers.MODEL, "high", 16_000
PAPERS, REVIEWS, FULL_TEXTS = 15, 5, 2
FULL_TEXT_CHARS = 40_000  # per full text in a request
KINDS = ["subspecies", "genomospecies or genomovar", "lineage", "phylogroup", "clade", "clonal complex",
         "sequence type", "serotype or serogroup", "ecotype", "other"]
TOPIC = ('("population structure"[tiab] OR "population genetics"[tiab] OR "population genomics"[tiab] '
         'OR phylogenomic*[tiab] OR "multilocus sequence typing"[tiab] OR "clonal complex*"[tiab] OR phylogroup*[tiab] '
         'OR ecotype*[tiab] OR genomovar*[tiab] OR lineages[tiab] OR "genomic diversity"[tiab])')


def search(species, reviews=REVIEWS, total=PAPERS, esearch=None, name=None):
    """PMIDs of the papers most relevant to a species' population structure: reviews with the species in the
    title first, then other papers with it in the title, then in the abstract. name: the name to search by
    instead, e.g. the genus of a species that is its genus's only one (papers say "Prochlorococcus"); a species
    with fewer than 5 such papers also falls back to its genus."""
    species = name or species
    esearch = esearch or (lambda term, n: json.loads(_eutils("esearch", db="pubmed", term=term, retmax=n,
                                                             sort="relevance", retmode="json"))
                          ["esearchresult"]["idlist"])
    found = []
    for term, n in ((f'"{species}"[ti] AND {TOPIC} AND review[pt]', reviews), (f'"{species}"[ti] AND {TOPIC}', total),
                    (f'"{species}"[tiab] AND {TOPIC}', total)):
        found += [p for p in esearch(term, n) if p not in found]
        if len(found) >= total:
            break
    if len(found) < 5:
        found += [p for p in esearch(f'"{species.split()[0]}"[ti] AND {TOPIC}', total) if p not in found]
    return found[:total]


def system_prompt():
    return (
        "You compile the published population structure of one bacterial species from the papers below. List the "
        "named groups the papers define within the species: subspecies, genomospecies or genomovars, lineages, "
        "phylogroups, clades, clonal complexes, sequence types, serotypes or serogroups, ecotypes. Nest them through "
        "parent (the name of the enclosing group, or empty for a top-level group). For each group give:\n"
        "- kind, and defined_by: how the papers define it (MLST, cgMLST, core-genome SNP phylogeny, serotyping, ...);\n"
        "- reference_strains: strains the papers name as members, as written (type strains, reference genomes);\n"
        "- ecology: hosts, sources, niches or traits the papers associate with the group, or empty;\n"
        "- evidence: one or more quotes copied verbatim from a listed paper (at most about 25 words each, checked "
        "against the text), each with that paper's PMID.\n"
        "Use only what the listed papers say, not your own knowledge. Keep to the main structure (at most about 25 "
        "groups): the top-level groups, and below them the subgroups the papers treat as major. If the papers define "
        "no structure, return no groups. summary: two or three sentences on how the species is structured. certain: "
        "false when the papers disagree or say little.")


def answer_schema(pmids):
    evidence = {"type": "object", "additionalProperties": False, "required": ["pmid", "quote"],
                "properties": {"pmid": {"type": "string", "enum": sorted(pmids)}, "quote": {"type": "string"}}}
    group = {"type": "object", "additionalProperties": False,
             "required": ["name", "kind", "parent", "defined_by", "reference_strains", "ecology", "evidence"],
             "properties": {"name": {"type": "string"}, "kind": {"type": "string", "enum": KINDS},
                            "parent": {"type": "string"}, "defined_by": {"type": "string"},
                            "reference_strains": {"type": "array", "items": {"type": "string"}},
                            "ecology": {"type": "string"},
                            "evidence": {"type": "array", "items": evidence}}}
    return {"type": "object", "additionalProperties": False, "required": ["summary", "groups", "certain"],
            "properties": {"summary": {"type": "string"}, "groups": {"type": "array", "items": group},
                           "certain": {"type": "boolean"}}}


def species_texts(pmids, texts, full_texts=FULL_TEXTS):
    """{pmid: {"kind", "text"}} for one request: abstracts, and the first full_texts full texts, shortened."""
    out, full = {}, 0
    for pmid in pmids:
        paper = texts.get(pmid)
        if not paper or paper["kind"] == "none":
            continue
        if paper["kind"] == "full text" and full < full_texts:
            out[pmid] = {"kind": "full text", "text": paper["text"][:FULL_TEXT_CHARS]}
            full += 1
        elif paper["kind"] == "full text":  # beyond the full-text allowance: its title and abstract only
            head = paper["text"].split("\n\n")[:2]
            out[pmid] = {"kind": "abstract", "text": "\n\n".join(head)}
        else:
            out[pmid] = {"kind": paper["kind"], "text": paper["text"]}
    return out


def build_request(species, chosen):
    """Messages API parameters for one species; chosen: species_texts' {pmid: {"kind", "text"}}."""
    user = f"Species: {species}\n\n" + "\n\n".join(
        f"Paper PMID {pmid} ({paper['kind']}):\n{paper['text']}" for pmid, paper in chosen.items())
    return {"model": MODEL, "max_tokens": MAX_TOKENS,
            "system": [{"type": "text", "text": system_prompt(), "cache_control": {"type": "ephemeral", "ttl": "1h"}}],
            "messages": [{"role": "user", "content": user}],
            "output_config": {"effort": EFFORT, "format": {"type": "json_schema", "schema": answer_schema(chosen)}}}


def worst_case_usd(requests):
    price_in, price_out = PRICES[MODEL]
    return sum((len(json.dumps(r["params"])) / 3.0 * price_in + r["params"]["max_tokens"] * price_out) / 1e6
               for r in requests) * BATCH_DISCOUNT


def verified(answer, chosen):
    """The answer's groups with only the evidence whose quote is in the cited paper; a group left with no
    evidence is marked unsupported. Returns (groups, quotes dropped)."""
    norm = {pmid: papers._norm(paper["text"]) for pmid, paper in chosen.items()}
    groups, dropped = [], 0
    for g in answer.get("groups") or []:
        kept = [e for e in g.get("evidence") or [] if papers._norm(e.get("quote")) and
                papers._norm(e["quote"]) in norm.get(e.get("pmid"), "")]
        dropped += len(g.get("evidence") or []) - len(kept)
        groups.append({**g, "evidence": kept, "supported": bool(kept)})
    return groups, dropped


def _strain_key(name):
    return re.sub(r"[^a-z0-9]", "", (name or "").lower())


_TAXON = re.compile(r"^(?:[A-Z][a-z]*\.?\s+[a-z]+\s+)(?:subsp\.\s+[a-z]+\s+)?(?:(?:str|strain|var)\.?\s+)?")


def strain_names(name):
    """The names one reference strain is written under: "Bacillus subtilis subsp. subtilis str. 168" is 168,
    "KT2440 (mt-2 = DSM 6125 = ATCC 47054)" is also mt-2, DSM 6125 and ATCC 47054, and a type strain's "T"
    mark is dropped ("Kh7 T" is Kh7)."""
    out = []
    for part in re.split(r"[()=;]", name or ""):
        part = _TAXON.sub("", part.strip())
        part = re.sub(r"\s+T$|\s+MAG$", "", part).strip()
        if part:
            out.append(part)
    return out


def catalogue_strains(release_tsv, species_names):
    """{species: {strain key: [accessions]}} from the release table's strain and isolate columns."""
    want, out = set(species_names), {s: collections.defaultdict(list) for s in species_names}
    with gzip.open(release_tsv, "rt", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["species"] in want:
                for name in [row["strain"], *re.split(r"[;,]", row["isolate"])]:
                    key = _strain_key(name)
                    if key:
                        out[row["species"]][key].append(row["accession"])
    return out


def match_strains(groups, strains):
    """For each group, the catalogue accessions of the reference strains it names: {group name: {strain:
    [accessions]}}. A strain matches a genome whose strain or isolate name is the same, ignoring case, spaces
    and punctuation; a name of fewer than 3 characters is not matched (too ambiguous)."""
    out = {}
    for g in groups:
        hits = {}
        for name in g.get("reference_strains") or []:
            accessions = {acc for alias in strain_names(name) for key in [_strain_key(alias)]
                          if len(key) >= 3 and key in strains for acc in strains[key]}
            if accessions:
                hits[name] = sorted(accessions)
        out[g["name"]] = hits
    return out
