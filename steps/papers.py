"""
Where did each genome come from, and what did the study that produced it do?

For every genome in ecotype_membership.csv: its NCBI BioProject and submitter, the PubMed
papers NCBI links to that BioProject, and a short Claude summary of what each paper did
(from its title and abstract). Many BioProjects have no linked paper (surveillance
uploads, sequencing centres); those genomes still get their BioProject link.

Writes into the EcoSim output folder:
  genome_sources.csv  accession, ecotype, bioproject, bioproject_title, submitter, pmids
  papers.csv          pmid, title, journal, year, doi, url, genomes, summary

Needs ANTHROPIC_API_KEY for the summaries (rows are still written without one).

  python steps/papers.py ecosim_output_dir "Treponema pallidum"
"""
import argparse
import csv
import json
import subprocess
import tempfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

try:
    from steps.niche_metadata import load_membership
except ImportError:  # run as a script: python steps/papers.py
    from niche_metadata import load_membership

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
LLM_MODEL = "claude-opus-5-5"
BATCH = 200          # ids per E-utilities request
SUMMARY_CHUNK = 25   # papers per Claude request


def _eutils(tool, **params):
    """POST one E-utilities request (POST, so long id lists fit). Sleeps to stay under NCBI's 3 requests/s."""
    data = urllib.parse.urlencode(params, doseq=True).encode()
    with urllib.request.urlopen(f"{EUTILS}/{tool}.fcgi", data=data, timeout=120) as response:
        body = response.read()
    time.sleep(0.4)
    return body


def fetch_assembly_projects(accessions):
    """{accession: (bioproject, submitter)} from `datasets summary genome accession`."""
    with tempfile.NamedTemporaryFile("w", suffix=".txt") as acc_file:
        acc_file.write("\n".join(accessions))
        acc_file.flush()
        out = subprocess.run(["datasets", "summary", "genome", "accession", "--inputfile", acc_file.name,
                              "--as-json-lines"], check=True, capture_output=True, text=True).stdout
    projects = {}
    for line in out.splitlines():
        record = json.loads(line)
        info = record.get("assembly_info", {})
        projects[record["accession"]] = (info.get("bioproject_accession", ""), info.get("submitter", ""))
    return projects


def fetch_bioproject_papers(bioprojects):
    """{bioproject accession: (title, [pmid, ...])} for the BioProjects NCBI knows."""
    found = {}
    bioprojects = sorted(bioprojects)
    for start in range(0, len(bioprojects), BATCH):
        chunk = bioprojects[start:start + BATCH]
        uids = json.loads(_eutils("esearch", db="bioproject", term=" OR ".join(chunk),
                                  retmax=len(chunk) * 2, retmode="json"))["esearchresult"]["idlist"]
        if not uids:
            continue
        summary = json.loads(_eutils("esummary", db="bioproject", id=",".join(uids), retmode="json"))["result"]
        acc_of = {uid: summary[uid]["project_acc"] for uid in summary.get("uids", [])}
        # One id= per uid keeps a separate linkset per BioProject.
        links = json.loads(_eutils("elink", dbfrom="bioproject", db="pubmed", id=uids, retmode="json"))
        pmids = {}
        for linkset in links.get("linksets", []):
            for db in linkset.get("linksetdbs", []):
                if db.get("linkname") == "bioproject_pubmed":
                    pmids[linkset["ids"][0]] = db["links"]
        for uid, acc in acc_of.items():
            if acc in chunk:
                found[acc] = (summary[uid].get("project_title", ""), pmids.get(uid, []))
    return found


def fetch_pubmed(pmids):
    """{pmid: {title, journal, year, doi, abstract}} from PubMed."""
    papers = {}
    pmids = sorted(pmids)
    for start in range(0, len(pmids), BATCH):
        root = ET.fromstring(_eutils("efetch", db="pubmed", id=",".join(pmids[start:start + BATCH]), retmode="xml"))
        for article in root.iter("PubmedArticle"):
            pmid = article.findtext(".//MedlineCitation/PMID")
            doi = next((a.text for a in article.iter("ArticleId") if a.get("IdType") == "doi"), "")
            papers[pmid] = {
                "title": "".join(article.find(".//ArticleTitle").itertext()).strip()
                         if article.find(".//ArticleTitle") is not None else "",
                "journal": article.findtext(".//Journal/Title") or "",
                "year": article.findtext(".//JournalIssue/PubDate/Year")
                        or (article.findtext(".//JournalIssue/PubDate/MedlineDate") or "")[:4],
                "doi": doi or "",
                "abstract": " ".join("".join(t.itertext()).strip() for t in article.iter("AbstractText")),
            }
    return papers


def summarize_with_claude(species, papers):
    """{pmid: summary} - what each study did, from its title and abstract only."""
    import anthropic  # only needed for this step

    print(f"Asking {LLM_MODEL} to summarise {len(papers)} paper(s)")
    client = anthropic.Anthropic()
    summaries = {}
    pmids = sorted(papers)
    for start in range(0, len(pmids), SUMMARY_CHUNK):
        chunk = pmids[start:start + SUMMARY_CHUNK]
        listing = "\n\n".join(f"PMID {p}\nTitle: {papers[p]['title']}\nAbstract: {papers[p]['abstract'] or '(none)'}"
                              for p in chunk)
        response = client.beta.messages.create(
            model=LLM_MODEL,
            max_tokens=16000,
            betas=["server-side-fallback-2026-07-01"],
            fallbacks="default",
            output_config={
                "effort": "low",
                "format": {"type": "json_schema", "schema": {
                    "type": "object",
                    "properties": {"summaries": {"type": "array", "items": {
                        "type": "object",
                        "properties": {"pmid": {"type": "string"}, "summary": {"type": "string"}},
                        "required": ["pmid", "summary"],
                        "additionalProperties": False,
                    }}},
                    "required": ["summaries"],
                    "additionalProperties": False,
                }},
            },
            messages=[{"role": "user", "content": (
                f"These papers are linked in NCBI to {species} genomes we analysed. For each, write 2-3 plain "
                "sentences on what the study did: which isolates it sequenced (source, place, time), what it "
                "set out to answer, and its main finding. Use only the title and abstract given; if the "
                "abstract is missing, say what the title implies and that no abstract was available.\n\n" + listing
            )}],
        )
        if response.stop_reason != "end_turn":
            raise RuntimeError(f"Claude stopped early: {response.stop_reason}")
        text = next(b.text for b in response.content if b.type == "text")
        for item in json.loads(text)["summaries"]:
            if item["pmid"] in chunk:
                summaries[item["pmid"]] = item["summary"].strip()
    return summaries


def run_paper_lookup(ecosim_out_dir, species, summarizer=summarize_with_claude):
    """Run the whole step on one EcoSim output folder. Returns the papers.csv path, or None."""
    out = Path(ecosim_out_dir)
    membership_csv = out / "ecotype_membership.csv"
    if not membership_csv.exists():
        print(f"No {membership_csv.name} in {out}; skipping the paper lookup.")
        return None
    _, ecotype_of = load_membership(membership_csv)
    accessions = sorted(a for a in ecotype_of if a.startswith(("GCA_", "GCF_")))

    projects = fetch_assembly_projects(accessions)
    bioprojects = fetch_bioproject_papers({p for p, _ in projects.values() if p})
    with open(out / "genome_sources.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["accession", "ecotype", "bioproject", "bioproject_title", "submitter", "pmids"])
        for acc in accessions:
            project, submitter = projects.get(acc, ("", ""))
            title, pmids = bioprojects.get(project, ("", []))
            writer.writerow([acc, ecotype_of[acc], project, title, submitter, ";".join(pmids)])

    genomes_of = defaultdict(list)
    for acc in accessions:
        for pmid in bioprojects.get(projects.get(acc, ("", ""))[0], ("", []))[1]:
            genomes_of[pmid].append(acc)
    papers = fetch_pubmed(genomes_of) if genomes_of else {}
    print(f"{len(accessions)} genomes, {len({p for p, _ in projects.values()})} BioProjects, {len(papers)} linked paper(s)")
    summaries = {}
    if papers:
        try:
            summaries = summarizer(species, papers)
        except Exception as exc:  # the citation table is still useful without summaries
            print(f"Warning: paper summaries failed: {exc}")

    papers_csv = out / "papers.csv"
    with open(papers_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["pmid", "title", "journal", "year", "doi", "url", "genomes", "summary"])
        for pmid, paper in sorted(papers.items()):
            writer.writerow([pmid, paper["title"], paper["journal"], paper["year"], paper["doi"],
                             f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/", len(genomes_of[pmid]),
                             summaries.get(pmid, "")])
    print(f"Saved {out / 'genome_sources.csv'} and {papers_csv}")
    return papers_csv


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ecosim_out_dir")
    parser.add_argument("species")
    parser.add_argument("--no-summaries", action="store_true", help="Skip the Claude summaries")
    args = parser.parse_args()
    run_paper_lookup(args.ecosim_out_dir, args.species,
                     summarizer=(lambda species, papers: {}) if args.no_summaries else summarize_with_claude)
