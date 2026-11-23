"""
Isolate metadata catalogue: curated host and isolation-source labels, each tied to the words it came from,
for every bacterial genome NCBI holds at complete-genome or chromosome level. See catalogue/README.md.

  python -m catalogue.build harvest                         # NCBI assembly reports + BioSample attributes (free)
  python -m catalogue.build prepare                         # taxonomy + ontology lookups, records, rules, cost (free)
  python -m catalogue.build run --budget 10 --sample 500    # pilot: the model tiers on 500 records (Batch API)
  python -m catalogue.build run --budget 150                # every record the rules and the cache leave open
  python -m catalogue.build verify --budget 30 --sample 500 # pilot: Claude Fable 5.2 cross-checks 500 records
  python -m catalogue.build verify --budget 200 --effort high  # Fable cross-checks every model label, judges disputes
  python -m catalogue.build rejudge --budget 100 --sample 300  # Opus 5.5 judges 300 of those disputes again
  python -m catalogue.build papers-pilot --budget 146 --sample 200  # Fable reads 200 linked papers: what they add
  python -m catalogue.build assemble                        # catalogue table, per-species summary, report
  python -m catalogue.build status                          # what each stage has done and what it has cost

Everything lives under --workdir (default ./catalogue_work). Only `run`, `verify`, `rejudge` and `papers-pilot`
call the Claude API; they read ANTHROPIC_API_KEY from the environment or from a .env file in the current
directory or the repository root.
"""
import argparse
import collections
import csv
import datetime
import gzip
import json
import math
import os
import random
import re
import subprocess
import sys
from pathlib import Path

from catalogue.batch import BatchTransport, BudgetExceeded, Ledger, estimate_usd
from catalogue.crosscheck import (EFFORTS, FABLE, JUDGE_EFFORT, PARTS, CrossCheck, compare_judges, judged_items,
                                  open_records, summarize, summarize_judges, third_party_sample, write_disagreements,
                                  write_judge_comparison)
from catalogue import papers
from catalogue.hosts import HostResolver
from catalogue.ontology import OntologyResolver, codes_to_resolve, resolve_codes
from catalogue.ncbi import harvest, iter_assemblies, mag_accessions, resolve_taxa
from steps.metadata_curation import (BATCH_DISCOUNT, HOST_FIELDS, PRICES, SOURCE_FIELDS, TIERS, AnswerCache,
                                     CuratedLabeler, Vocabulary, normalize, usage_cost)
from steps.niche_metadata import LLM_FIELDS, country_and_year


def foreign_modules():
    """steps/ has no __init__.py, so Python joins every steps/ folder on sys.path into one package. Run from
    inside another checkout of this repository (whose folder comes first on sys.path), that checkout's
    steps/ modules silently replace this one's. Returns the files of any that came from elsewhere."""
    repo = Path(__file__).resolve().parents[1]
    files = [sys.modules[name].__file__ for name in ("steps.metadata_curation", "steps.niche_metadata")]
    return [f for f in files if Path(f).resolve().parents[1] != repo]


CATALOGUE_VERSION = "0.1"
# BioSample fields the curation reads: the pipeline's fields plus non-standard ones that, across the
# harvested genomes, carry host or source information (habitat, ENA's names for the MIxS environment
# fields, the metagenome a MAG came from, ...). Names are lowercased with spaces and hyphens as "_".
EXTRA_FIELDS = ("host_common_name", "isolation_site", "source", "isolation_source_category", "host_body_habitat",
                "plant_body_site", "tissue", "food_source", "habitat", "environment", "metagenome_source",
                "coll_site_geo_feat")
CATALOGUE_FIELDS = (*LLM_FIELDS, *EXTRA_FIELDS)
SOURCE_PRIORITY = (*SOURCE_FIELDS, "isolation_site", "source", "isolation_source_category", "host_body_habitat",
                   "plant_body_site", "tissue", "food_source", "habitat", "environment", "metagenome_source",
                   "coll_site_geo_feat")
HOST_PRIORITY = (*HOST_FIELDS, "host_common_name")
ALIASES = {"local_environmental_context": "env_local_scale",
           "broad_scale_environmental_context": "env_broad_scale",
           "environmental_medium": "env_medium",
           "metagenomic_source": "metagenome_source",
           "isolation_source_host_associated": "isolation_source"}
RANKS = ("species", "genus", "family", "order", "class", "phylum")


def canonical(attributes):
    """BioSample attributes under canonical names (lowercase, '_' for spaces and hyphens, ENA aliases)."""
    out = {}
    for name, value in attributes.items():
        out.setdefault(re.sub(r"[\s\-]+", "_", name.strip().lower()), value)
    for alias, target in ALIASES.items():
        if alias in out and target not in out:
            out[target] = out[alias]
    return out


def record_key(attributes):
    """The distinct-record key the curation works on: (field, value) pairs in CATALOGUE_FIELDS order."""
    attrs = canonical(attributes)
    return tuple((f, attrs[f]) for f in CATALOGUE_FIELDS if attrs.get(f))


class Paths:
    def __init__(self, workdir):
        self.root = Path(workdir).resolve()
        self.ncbi = self.root / "ncbi"
        self.taxa = self.root / "taxa_cache.json"
        self.ontology = self.root / "ontology_cache.json"
        self.cache = self.root / "curation_cache.json"
        self.crosscheck = self.root / "crosscheck_cache.json"
        self.reports = self.root / "reports"
        self.release = self.root / "release"


def _now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def _require_harvest(paths):
    state_path = paths.ncbi / "harvest_state.json"
    if not state_path.exists() or not json.loads(state_path.read_text()).get("done"):
        sys.exit("No complete harvest yet: run `python -m catalogue.build harvest` first.")
    return json.loads(state_path.read_text())


def load_env():
    """Read ANTHROPIC_API_KEY (and anything else) from .env in the current directory, then the repo root."""
    from steps.dotenv import ENV_FILE, load_dotenv
    for path in (Path.cwd() / ".env", ENV_FILE):
        load_dotenv(path)


def _require_key():
    load_env()
    if not (os.environ.get("ANTHROPIC_API_KEY") or os.environ.get("ANTHROPIC_AUTH_TOKEN")):
        sys.exit("No ANTHROPIC_API_KEY. Put `ANTHROPIC_API_KEY=...` in a .env file in this directory (it is "
                 "gitignored) or export it, then rerun.")


def make_labeler(paths, vocab, resolver, client=None, transport=None, log=print):
    return CuratedLabeler(vocab=vocab, client=client, cache_path=paths.cache, fields=CATALOGUE_FIELDS, log=log,
                          host_resolver=resolver, transport=transport, source_fields=SOURCE_PRIORITY,
                          host_fields=HOST_PRIORITY, source_from_host=True,
                          term_resolver=load_ontology(paths, vocab))


def load_ontology(paths, vocab):
    """The ontology codes `prepare` placed (an empty resolver before it has run)."""
    return OntologyResolver(vocab, json.loads(paths.ontology.read_text()) if paths.ontology.exists() else {})


def resolve_all(paths, vocab, log=print, fetch=None):
    """Look up every host name, host taxid and organism taxid in NCBI Taxonomy (cached). Returns
    (HostResolver, taxa)."""
    host_values, host_taxids, organism_ids = set(), set(), set()
    for assembly in iter_assemblies(paths.ncbi):
        attrs = canonical(assembly["attributes"])
        for field in HOST_PRIORITY:
            if vocab.informative(attrs.get(field, "")):
                host_values.add(attrs[field])
        if str(attrs.get("host_taxid", "")).strip():
            host_taxids.add(attrs["host_taxid"])
        if assembly.get("tax_id"):
            organism_ids.add(str(assembly["tax_id"]))
    queries = HostResolver.queries(host_values, host_taxids) | organism_ids
    log(f"Taxonomy: {len(host_values):,} distinct host names, {len(host_taxids):,} host taxids, "
        f"{len(organism_ids):,} organism taxids ({len(queries):,} lookups, cached in {paths.taxa.name})")
    kwargs = {"fetch": fetch} if fetch else {}
    taxa = resolve_taxa(queries, paths.taxa, log=log, **kwargs)
    return HostResolver(vocab, taxa), taxa


def load_resolver(paths, vocab):
    taxa = json.loads(paths.taxa.read_text()) if paths.taxa.exists() else {}
    return HostResolver(vocab, taxa), taxa


def distinct_keys(paths):
    counts = collections.Counter(record_key(a["attributes"]) for a in iter_assemblies(paths.ncbi))
    return counts


def split_by_rules(labeler, keys):
    """(keys settled by the rules, {cache key: (key, normalized fields)} still open, open cache keys already
    answered)."""
    records, key_to_ck = open_records(labeler, keys)
    settled = [key for key in keys if key not in key_to_ck]
    pending = {ck: record for ck, record in records.items() if not labeler.cache.get(ck)}
    return settled, pending, set(records) - set(pending)


# --- cost estimate ----------------------------------------------------------------------------------
def estimate_costs(labeler, open_records, escalate, dispute, expected=True):
    """USD at list price for the model tiers on open_records ({rid: fields}): Haiku on all of them,
    two Sonnet runs on the escalated share, Opus on the disputed share of those. expected=False gives
    the conservative figure the budget checks use."""
    chars_per_token = 3.6 if expected else 3.0
    out_per_record = ({"haiku": 150, "sonnet": 250, "opus": 500} if expected
                      else {"haiku": 300, "sonnet": 450, "opus": 900})
    ids = list(open_records)

    def tier_cost(tier, rids, disputes=None):
        cfg = labeler.tiers[tier]
        price_in, price_out = PRICES[cfg["model"]]
        chars = 0
        for start in range(0, len(rids), cfg["chunk"]):
            part = {i: open_records[i] for i in rids[start:start + cfg["chunk"]]}
            chars += len(json.dumps(labeler.build_request(tier, part, disputes and {i: disputes for i in part})))
        return (chars / chars_per_token * price_in + len(rids) * out_per_record[tier] * price_out) / 1e6

    sample = random.Random(1).sample(ids, min(len(ids), 3000))  # Sonnet/Opus sizes from a sample, then scaled
    scale = len(ids) / max(len(sample), 1)
    fake_dispute = {"host": ({"host": "unknown", "host_field": "", "host_quote": ""},) * 2}
    haiku = tier_cost("haiku", ids)
    sonnet = 2 * tier_cost("sonnet", sample) * scale * escalate
    opus = tier_cost("opus", sample, fake_dispute) * scale * escalate * dispute
    return {"haiku": haiku, "sonnet x2": sonnet, "opus": opus, "total": haiku + sonnet + opus}


# --- commands --------------------------------------------------------------------------------------
def cmd_harvest(args, paths, log):
    harvest(paths.ncbi, taxon=args.taxon, levels=args.levels.split(","), log=log, restart=args.restart)


def cmd_prepare(args, paths, log):
    state = _require_harvest(paths)
    vocab = Vocabulary()
    resolver, taxa = resolve_all(paths, vocab, log=log)
    counts = distinct_keys(paths)
    codes = codes_to_resolve(make_labeler(paths, vocab, resolver, log=lambda *_: None), counts)
    resolve_codes(codes, paths.ontology, vocab, log=log)
    labeler = make_labeler(paths, vocab, resolver, log=lambda *_: None)  # with the codes just placed
    placed = sum(1 for c in codes if labeler.term_resolver.vocab_key(c))
    settled, pending, cached = split_by_rules(labeler, list(counts))
    genomes_settled = sum(counts[k] for k in settled)
    no_info = sum(n for k, n in counts.items() if not any(vocab.informative(v) for _, v in k))
    open_genomes = collections.Counter()
    unresolved = collections.Counter()
    for cache_key, (key, norm) in pending.items():
        open_genomes[cache_key] += counts[key]
    for key, n in counts.items():
        fields = {k: v for k, v in key if vocab.informative(v)}
        if labeler.rule_source(fields) is None:
            first = next((fields[f] for f in SOURCE_PRIORITY if vocab.informative(fields.get(f, ""))), "")
            unresolved[normalize(first)[:80]] += n
    open_records = {i: norm for i, (_, (key, norm)) in enumerate(pending.items())}
    scenarios = {
        "expected (30% escalate past Haiku, 30% of those disputed)": estimate_costs(labeler, open_records, 0.3, 0.3),
        "high (every record escalates, half disputed)": estimate_costs(labeler, open_records, 1.0, 0.5),
        "budget check ceiling (high, conservative token counts)": estimate_costs(labeler, open_records, 1.0, 0.5,
                                                                                 expected=False),
    }
    hosts_resolved = sum(1 for q, r in taxa.items() if r and "ambiguous" not in r)
    summary = {
        "prepared": _now(), "assemblies": sum(counts.values()), "mags": len(mag_accessions(paths.ncbi)),
        "distinct_records": len(counts), "records_settled_by_rules": len(settled),
        "genomes_settled_by_rules": genomes_settled, "genomes_without_any_source_or_host_field": no_info,
        "open_records": len(pending), "open_records_already_cached": len(cached),
        "genomes_in_open_records": sum(open_genomes.values()),
        "taxonomy_lookups": len(taxa), "taxonomy_resolved": hosts_resolved,
        "ontology_codes": len(codes), "ontology_codes_placed_in_vocabulary": placed,
        "cost_usd_list": {name: {k: round(v, 2) for k, v in s.items()} for name, s in scenarios.items()},
        "cost_usd_batch": {name: round(s["total"] * BATCH_DISCOUNT, 2) for name, s in scenarios.items()},
        "harvest": {k: state.get(k) for k in ("records", "finished", "query")},
    }
    paths.reports.mkdir(parents=True, exist_ok=True)
    (paths.reports / "prepare.json").write_text(json.dumps(summary, indent=1))
    lines = [f"# Catalogue preparation, {summary['prepared']}", "",
             f"- Assemblies: {summary['assemblies']:,} ({summary['mags']:,} metagenome-assembled)",
             f"- Distinct BioSample records (field combinations): {len(counts):,}",
             f"- Settled by the rules: {len(settled):,} records covering {genomes_settled:,} genomes "
             f"(of which {no_info:,} genomes have no source or host field at all, so they are unknown by rule)",
             f"- Open for the models: {len(pending):,} records covering {summary['genomes_in_open_records']:,} genomes"
             f" ({len(cached):,} more already answered in the cache)",
             f"- Ontology codes in the deciding source field: {len(codes):,}, of which {placed:,} lead to a vocabulary "
             f"term within {labeler.term_resolver.max_steps} steps up ENVO, UBERON, FoodOn or PO",
             "", "| Scenario | Haiku | Sonnet x2 | Opus | Total, list | Total, Batch API |", "|---|---|---|---|---|---|"]
    for name, s in scenarios.items():
        lines.append(f"| {name} | ${s['haiku']:.2f} | ${s['sonnet x2']:.2f} | ${s['opus']:.2f} | ${s['total']:.2f} "
                     f"| ${s['total'] * BATCH_DISCOUNT:.2f} |")
    lines += ["", "Most common source values the rules could not map (genomes):", ""]
    lines += [f"- {n:,}  {v or '(no source field; host only)'}" for v, n in unresolved.most_common(40)]
    (paths.reports / "prepare_report.md").write_text("\n".join(lines) + "\n")
    for line in lines[:9]:
        log(line)
    for name, s in scenarios.items():
        log(f"  {name}: ${s['total']:.2f} list, ${s['total'] * BATCH_DISCOUNT:.2f} batched")
    log(f"Report: {paths.reports / 'prepare_report.md'}")
    return summary


def cmd_run(args, paths, log):
    _require_harvest(paths)
    _require_key()
    if args.sync and (not args.sample or args.sample > 200):
        sys.exit("--sync (full price, no batching) is only for small tests: use it with --sample 200 or fewer.")
    vocab = Vocabulary()
    resolver, _ = load_resolver(paths, vocab)
    if not paths.taxa.exists():
        resolver, _ = resolve_all(paths, vocab, log=log)
    counts = distinct_keys(paths)
    keys = sorted(counts)
    probe = make_labeler(paths, vocab, resolver, log=lambda *_: None)
    _, pending, _ = split_by_rules(probe, keys)
    label = "full"
    if args.sample:
        chosen = random.Random(args.seed).sample(sorted(pending), min(args.sample, len(pending)))
        keys = sorted(pending[c][0] for c in chosen)
        label = f"pilot{len(chosen)}"
    log(f"{len(pending):,} open records; this run sends {len(keys) if args.sample else len(pending):,} "
        f"to the models (budget ${args.budget:.2f})")
    import anthropic
    client = anthropic.Anthropic()
    transport = None if args.sync else BatchTransport(paths.root, args.budget, client=client, log=log,
                                                      poll_seconds=args.poll, stage_prefix=f"{label} ")
    labeler = make_labeler(paths, vocab, resolver, client=client, transport=transport, log=log)
    started = _now()
    try:
        details = labeler.curate(keys)
    except BudgetExceeded as e:
        sys.exit(f"Stopped before spending more: {e}")
    sent = len(keys) if args.sample else len(pending)
    tiers = collections.Counter()
    open_keys = {v[0] for v in pending.values()}
    for key, d in details.items():
        if args.sample or key in open_keys:
            tiers[("source", d["source_tier"])] += 1
            tiers[("host", d["host_tier"])] += 1
    ledger = Ledger(paths.root / "spend.jsonl")
    report = {"run": label, "started": started, "finished": _now(), "records_sent": sent,
              "stats": dict(labeler.stats), "usd_this_run": round(labeler.cost(), 4),
              "usd_ledger_total": round(ledger.total, 4),
              "tiers": {f"{p} {t}": n for (p, t), n in sorted(tiers.items())}}
    if args.sample and sent:
        per_record = labeler.cost() / sent
        report["usd_per_record"] = round(per_record, 6)
        report["projected_usd_for_all_open_records"] = round(per_record * len(pending), 2)
    paths.reports.mkdir(parents=True, exist_ok=True)
    out = paths.reports / f"run_{label}_{started.replace(':', '')}.json"
    out.write_text(json.dumps(report, indent=1))
    log(json.dumps(report, indent=1))
    log(f"Report: {out}")
    return report


def cmd_verify(args, paths, log):
    _require_harvest(paths)
    _require_key()
    vocab = Vocabulary()
    resolver, _ = load_resolver(paths, vocab)
    counts = distinct_keys(paths)
    main = make_labeler(paths, vocab, resolver, log=lambda *_: None)
    records, _ = open_records(main, sorted(counts))
    total = len(records)
    label = f"all-{args.effort}"
    if args.sample:
        chosen = random.Random(args.seed).sample(sorted(records), min(args.sample, total))
        records = {ck: records[ck] for ck in chosen}
        label = f"pilot{len(records)}-{args.effort}"
    log(f"{total:,} records were left to the models; this run cross-checks {len(records):,} with {FABLE} "
        f"(budget ${args.budget:.2f} for everything spent in this work directory)")
    import anthropic
    client = anthropic.Anthropic()
    transport = BatchTransport(paths.root, args.budget, client=client, log=log, poll_seconds=args.poll,
                               stage_prefix=f"verify {label} ", state_name="crosscheck_batches.json")
    xc = CrossCheck(main, paths.crosscheck, effort=args.effort, client=client, transport=transport, log=log)
    started = _now()
    try:
        asked = xc.label(records)
        label_usd = xc.cost()
        results = xc.check(records)  # Fable's labels are cached now: this compares and judges
    except BudgetExceeded as e:
        sys.exit(f"Stopped before spending more: {e}")
    except anthropic.BadRequestError as e:
        sys.exit(f"The API rejected a {FABLE} request: {e.message}")
    judge_usd = xc.cost() - label_usd
    compared = [r for r in results.values() if r is not None]
    disputed = sum(1 for r in compared if any(part["status"].startswith(("judge", "disagrees"))
                                             for part in r["parts"].values()))
    others = {o: xc.effort_agreement(records, o) for o in EFFORTS if o != args.effort}
    report = {"run": label, "started": started, "finished": _now(), "model": FABLE, "effort": args.effort,
              "open records": total, **summarize(results),
              "records judged (model-decided disagreement)": disputed,
              "same labels as Fable at other efforts": {o: a for o, a in others.items() if a},
              "stats": dict(xc.stats), "usd_this_run": round(xc.cost(), 4),
              "usd_ledger_total": round(Ledger(paths.root / "spend.jsonl").total, 4)}
    if asked:
        report["label usd per record"] = round(label_usd / asked, 6)
        report["projected label usd, every open record"] = round(label_usd / asked * total, 2)
    if xc.stats["judged"] and compared:
        per_judged = judge_usd / xc.stats["judged"]
        report["judge usd per judged record"] = round(per_judged, 6)
        report["projected judge usd, every open record"] = round(per_judged * disputed / len(compared) * total, 2)
    paths.reports.mkdir(parents=True, exist_ok=True)
    stem = f"verify_{label}_{started.replace(':', '')}"
    (paths.reports / f"{stem}.json").write_text(json.dumps(report, indent=1))
    write_disagreements(results, paths.reports / f"{stem}_disagreements.tsv", vocab)
    for p in PARTS:
        outcomes = report["parts"][p]["all"]
        n = sum(outcomes.values())
        agree = sum(v for k, v in outcomes.items() if k in ("agrees", "same niche"))
        log(f"  {p}: Fable agrees on {agree:,} of {n:,} ({agree / max(n, 1):.0%}); "
            + ", ".join(f"{k} {v:,}" for k, v in outcomes.items()))
    log(f"  cost this run ${xc.cost():.2f}; spent in this work directory ${report['usd_ledger_total']:.2f}")
    log(f"Report: {paths.reports / stem}.json (disagreements in {stem}_disagreements.tsv)")
    return report


def cmd_rejudge(args, paths, log):
    """A second model judges a random sample of the disputes Fable's judge settled, from the same prompts, so the
    two judges differ only in the model: a check on Fable's judge favouring Fable's own answers. Nothing it
    decides changes the catalogue."""
    if args.model == FABLE or args.model not in PRICES:
        sys.exit(f"--model must be another priced model than {FABLE}: {', '.join(m for m in PRICES if m != FABLE)}")
    _require_harvest(paths)
    _require_key()
    vocab = Vocabulary()
    resolver, _ = load_resolver(paths, vocab)
    counts = distinct_keys(paths)
    main = make_labeler(paths, vocab, resolver, log=lambda *_: None)
    records, _ = open_records(main, sorted(counts))
    first = CrossCheck(main, paths.crosscheck, lookup=tuple(reversed(EFFORTS)), log=lambda *_: None)
    results = first.check(records, offline=True)
    items = judged_items(results)
    sample, eligible = third_party_sample(items, args.model, args.sample, args.seed)
    log(f"Fable's judge settled {len(items):,} records; {eligible:,} have no curated answer from {args.model}; "
        f"{len(sample):,} sampled for it to judge again (budget ${args.budget:.2f} for everything spent in this "
        "work directory)")
    import anthropic
    client = anthropic.Anthropic()
    label = f"{args.model}-{len(sample)}"
    transport = BatchTransport(paths.root, args.budget, client=client, log=log, poll_seconds=args.poll,
                               stage_prefix=f"rejudge {label} ", state_name="crosscheck_batches.json")
    second = CrossCheck(main, paths.crosscheck, client=client, transport=transport, log=log, judge_model=args.model)
    started = _now()
    try:
        decisions = second.judge(sample)
    except BudgetExceeded as e:
        sys.exit(f"Stopped before spending more: {e}")
    except anthropic.BadRequestError as e:
        sys.exit(f"The API rejected a {args.model} request: {e.message}")
    rows = compare_judges(results, sample, decisions)
    agreement = summarize_judges(rows)
    report = {"run": f"rejudge {label}", "started": started, "finished": _now(), "first judge": FABLE,
              "second judge": args.model, "effort": JUDGE_EFFORT, "seed": args.seed,
              "records Fable's judge settled": len(items), f"records with no curated answer from {args.model}": eligible,
              "records sampled": len(sample), "records the second judge decided": sum(ck in decisions for ck in sample),
              "second judge answered by": dict(collections.Counter(r["second"].model for r in rows)),
              "agreement": agreement, "stats": dict(second.stats), "usd_this_run": round(second.cost(), 4),
              "usd_ledger_total": round(Ledger(paths.root / "spend.jsonl").total, 4)}
    paths.reports.mkdir(parents=True, exist_ok=True)
    stem = f"rejudge_{label}_{started.replace(':', '')}"
    (paths.reports / f"{stem}.json").write_text(json.dumps(report, indent=1))
    write_judge_comparison(rows, paths.reports / f"{stem}.tsv", vocab)
    for name in ("all", *PARTS):
        a = agreement.get(name)
        if not a:
            continue
        took = (a["took Fable's, first judge"], a["took Fable's, second judge"])
        only = (a["first took Fable's, second kept curation"], a["first kept curation, second took Fable's"])
        log(f"  {name}: same decision on {a['same decision']:,} of {a['parts']:,} parts "
            f"({a['same decision share']:.0%}); took Fable's answer: Fable's judge {took[0]:,}, {args.model} "
            f"{took[1]:,}; only Fable's judge took it over the curation's {only[0]:,}, only {args.model} "
            f"{only[1]:,} (McNemar p = {a['McNemar exact p']:.3g})")
    log(f"  cost this run ${second.cost():.2f}; spent in this work directory ${report['usd_ledger_total']:.2f}")
    log(f"Report: {paths.reports / stem}.json (each part in {stem}.tsv)")
    return report


def cmd_papers_pilot(args, paths, log):
    """Read a sample of the papers NCBI links to the catalogue's genomes and measure what they add (see
    catalogue/papers.py). Needs `assemble` first: genomes are compared with the release table."""
    _require_harvest(paths)
    release = paths.release / "isolate_catalogue.tsv.gz"
    if not release.exists():
        sys.exit("No catalogue yet: run `python -m catalogue.build assemble` first.")
    vocab = Vocabulary()
    folder = paths.root / "papers"
    genome_project, project_genomes = {}, collections.defaultdict(list)
    for a in iter_assemblies(paths.ncbi):
        if a.get("bioproject"):
            genome_project[a["accession"]] = a["bioproject"]
            project_genomes[a["bioproject"]].append(a["accession"])
    log(f"Sampling genomes until {args.sample} distinct linked papers (budget ${args.budget:.2f} for everything "
        "spent in this work directory)")
    chosen, links, tried, with_paper = papers.sample_papers(genome_project, folder / "links.json", args.sample,
                                                            args.seed)
    log(f"  {tried:,} genomes tried, {with_paper:,} with a linked paper ({with_paper / tried:.0%}); "
        f"{len(chosen):,} distinct papers")
    texts = papers.paper_texts(sorted(chosen), folder / "texts")
    kinds = collections.Counter(t["kind"] for t in texts.values())
    log("  texts: " + ", ".join(f"{k} {v:,}" for k, v in kinds.most_common()))
    pmid_projects = collections.defaultdict(set)
    for project, pmids in links.items():
        for pmid in pmids:
            if pmid in chosen:
                pmid_projects[pmid].add(project)
    rng = random.Random(args.seed)
    listed, more = {}, {}
    for pmid, sampled in chosen.items():
        first = list(dict.fromkeys(sampled))
        others = sorted({g for p in pmid_projects[pmid] for g in project_genomes[p]} - set(first))
        rng.shuffle(others)
        listed[pmid] = (first + others)[:max(papers.MAX_GENOMES, len(first))]
        more[pmid] = len(first) + len(others) - len(listed[pmid])
    rows = papers.catalogue_rows(release, {g for gs in listed.values() for g in gs})
    requests, genome_lists = {}, {}
    for pmid in sorted(chosen):
        if texts[pmid]["kind"] == "none":
            continue
        genome_lists[pmid] = genomes = [{"accession": g, "organism": rows[g]["organism_name"],
                                         "strain": rows[g]["strain"], "isolate": rows[g]["isolate"]}
                                        for g in listed[pmid] if g in rows]
        requests[pmid] = {"custom_id": f"pmid{pmid}",
                          "params": papers.build_request(vocab, pmid, texts[pmid], genomes, more[pmid])}
    worst = papers.worst_case_usd(list(requests.values()))
    sampled = sum(len(set(s)) for s in chosen.values())
    log(f"  {len(requests):,} papers to read, {sampled:,} sampled genomes, "
        f"{sum(len(listed[p]) for p in requests):,} genomes listed; worst case ${worst:.2f} at Batch prices "
        "(every request using all of max_tokens)")
    if args.dry_run:
        return {"papers": len(requests), "texts": dict(kinds), "worst_case_usd": round(worst, 2)}
    # Each paper's answer is kept in papers/answers/, so a rerun (or a new report) pays nothing for it again. The
    # waves are always the same papers in the same order, so their batches resume instead of being sent again.
    answers_dir = folder / "answers"
    answers_dir.mkdir(parents=True, exist_ok=True)
    saved = {p: json.loads((answers_dir / f"{p}.json").read_text())
             for p in requests if (answers_dir / f"{p}.json").exists()}

    def keep(pmid, answer, reason, usd, how):
        saved[pmid] = {"answer": answer, "reason": reason, "usd": usd, "how": how}
        (answers_dir / f"{pmid}.json").write_text(json.dumps(saved[pmid]))

    smoke = next(iter(sorted(requests, key=lambda p: len(texts[p]["text"]))), None)
    rest = [p for p in sorted(requests) if p != smoke]
    ledger_before = Ledger(paths.root / "spend.jsonl").total  # resumed batches are not charged twice

    def retryable():
        return [p for p in sorted(requests) if args.retry_model and p in saved and saved[p]["answer"] is None
                and saved[p]["reason"] in ("refusal", "max_tokens") and "first_reason" not in saved[p]]

    def refused():  # still refused after any retry, and not yet asked the narrow way
        return [p for p in sorted(requests) if args.narrow_refused and p in saved and saved[p]["answer"] is None
                and saved[p]["reason"] == "refusal" and "narrow_reason" not in saved[p]]

    if any(p not in saved for p in requests) or retryable() or refused():
        _require_key()
        import anthropic
        client = anthropic.Anthropic()
        ledger = Ledger(paths.root / "spend.jsonl")
        # One paper synchronously first, at full price, so a request the API rejects stops before a batch.
        if smoke not in saved:
            worst = papers.worst_case_usd([requests[smoke]]) / BATCH_DISCOUNT
            if ledger.total + worst > args.budget:
                sys.exit(f"Stopped before spending more: one paper could cost ${worst:.2f}; ${ledger.total:.2f} of "
                         f"the ${args.budget:.2f} budget is already spent.")
            try:
                message = client.messages.create(**requests[smoke]["params"])
            except anthropic.BadRequestError as e:
                sys.exit(f"The API rejected the paper request: {e.message}")
            cost = usage_cost(papers.MODEL, message.usage)
            ledger.add(stage="papers pilot first paper", kind="sync", model=papers.MODEL, usd=round(cost, 9),
                       requests=1)
            keep(smoke, *papers.parse_message(message), cost, "sync")
            log(f"  first paper (PMID {smoke}) answered synchronously for ${cost:.3f}")
        transport = BatchTransport(paths.root, args.budget, client=client, log=log, poll_seconds=args.poll,
                                   stage_prefix="papers pilot ", state_name="papers_batches.json")
        recorder = CuratedLabeler(vocab=vocab, log=lambda *_: None)
        try:
            for start in range(0, len(rest), papers.WAVE):
                wave = rest[start:start + papers.WAVE]
                if all(p in saved for p in wave):
                    continue
                got = transport.run_requests(f"wave {start // papers.WAVE + 1}", [requests[p] for p in wave],
                                             recorder, estimate=papers.worst_case_usd([requests[p] for p in wave]))
                for pmid in wave:
                    kind, message = got.get(f"pmid{pmid}", ("missing", None))
                    ok = kind == "succeeded" and message is not None
                    answer, reason = papers.parse_message(message) if ok else (None, kind)
                    keep(pmid, answer, reason, usage_cost(papers.MODEL, message.usage, batch=True) if ok else 0.0,
                         "batch")
            # Papers Fable refused (its safety screening declines some pathogen and toxin papers) or cut off,
            # once more with another model and more room. Batches have no server-side refusal fallback.
            retry = retryable()
            if retry:
                again = [{"custom_id": f"pmid{p}", "params": {**requests[p]["params"], "model": args.retry_model,
                                                              "max_tokens": papers.RETRY_MAX_TOKENS}} for p in retry]
                got = transport.run_requests(f"retry {args.retry_model}", again, recorder,
                                             estimate=papers.worst_case_usd(again))
                for p in retry:
                    kind, message = got.get(f"pmid{p}", ("missing", None))
                    ok = kind == "succeeded" and message is not None
                    answer, reason = papers.parse_message(message) if ok else (None, kind)
                    cost = usage_cost(args.retry_model, message.usage, batch=True) if ok else 0.0
                    saved[p] = {**saved[p], "first_reason": saved[p]["reason"], "answer": answer, "reason": reason,
                                "usd": saved[p]["usd"] + cost, "how": f"retried with {args.retry_model}"}
                    (answers_dir / f"{p}.json").write_text(json.dumps(saved[p]))
            # Papers still refused, asked the narrow way: the catalogue's own facts from the methods and tables.
            narrow = refused()
            texts_narrow = {p: papers.narrow_text(texts[p]) for p in narrow}
            for p in narrow:
                if texts_narrow[p] is None:  # a full text with neither methods-like sections nor tables
                    saved[p] = {**saved[p], "narrow_reason": "no methods or tables"}
                    (answers_dir / f"{p}.json").write_text(json.dumps(saved[p]))
            todo = [p for p in narrow if texts_narrow[p] is not None]
            for start in range(0, len(todo), papers.NARROW_WAVE):
                wave = todo[start:start + papers.NARROW_WAVE]
                again = [{"custom_id": f"pmid{p}", "params": papers.build_request(
                    vocab, p, texts_narrow[p], genome_lists[p], more[p], facts=papers.NARROW_FACTS,
                    max_tokens=papers.NARROW_MAX_TOKENS)} for p in wave]
                got = transport.run_requests(f"narrow {start // papers.NARROW_WAVE + 1}", again, recorder,
                                             estimate=papers.worst_case_usd(again))
                for p in wave:
                    kind, message = got.get(f"pmid{p}", ("missing", None))
                    ok = kind == "succeeded" and message is not None
                    answer, reason = papers.parse_message(message) if ok else (None, kind)
                    cost = usage_cost(papers.MODEL, message.usage, batch=True) if ok else 0.0
                    saved[p] = {**saved[p], "narrow_reason": reason, "answer": answer, "reason": reason,
                                "usd": saved[p]["usd"] + cost, "how": "narrow request"}
                    (answers_dir / f"{p}.json").write_text(json.dumps(saved[p]))
        except BudgetExceeded as e:
            log(f"Stopped before spending more: {e}")
    spent_now = Ledger(paths.root / "spend.jsonl").total - ledger_before
    usd = {p: s["usd"] for p, s in saved.items()}
    out_rows, failures, describes, quote_failures = [], collections.Counter(), collections.Counter(), 0
    for pmid in sorted(requests):
        answer, reason = (saved[pmid]["answer"], saved[pmid]["reason"]) if pmid in saved else (None, "not sent")
        if answer is None:
            failures[reason] += 1
            continue
        describes[answer.get("describes", "")] += 1
        for g in listed[pmid]:
            if g not in rows:
                continue
            facts, origin, failed = papers.facts_for(answer, g, texts[pmid]["text"])
            quote_failures += failed
            outcome = papers.compare(rows[g], facts, vocab)
            out_rows.append({"pmid": pmid, "kind": texts[pmid]["kind"], "accession": g, "sampled": g in chosen[pmid],
                             "answered_by": saved[pmid]["how"], "origin": origin, **outcome,
                             **{f"{f}_value": facts[f][0] if f in facts else "" for f in papers.FACTS},
                             **{f"{f}_quote": facts[f][1] if f in facts else "" for f in papers.FACTS},
                             **{f"catalogue_{c}": rows[g][c] for c in ("host", "source_term", "country",
                                                                         "collection_year")}})
    answered = len(requests) - sum(failures.values())
    per_kind = collections.defaultdict(list)
    for pmid, cost in usd.items():
        per_kind[texts[pmid]["kind"]].append(cost)
    report = {"run": "papers pilot", "finished": _now(), "model": papers.MODEL, "effort": papers.EFFORT,
              "seed": args.seed, "genomes tried": tried, "genomes with a linked paper": with_paper,
              "share of genomes with a linked paper": round(with_paper / tried, 3),
              "papers": len(chosen), "texts": dict(kinds), "papers answered": answered,
              "papers without an answer": dict(failures), "how papers describe their genomes": dict(describes),
              "Fable's first answer": dict(collections.Counter(s.get("first_reason", s["reason"]) for s in saved.values())),
              "answered by": dict(collections.Counter(s["how"] for s in saved.values() if s["answer"] is not None)),
              "quotes not found in the paper (facts dropped)": quote_failures,
              "sampled genomes": papers.summarize(out_rows),
              "sampled genomes, papers answered after a retry": papers.summarize(
                  [r for r in out_rows if r["answered_by"].startswith("retried")]),
              "narrow request": dict(collections.Counter(s["narrow_reason"] for s in saved.values()
                                                         if "narrow_reason" in s)),
              "sampled genomes, papers answered by the narrow request": papers.summarize(
                  [r for r in out_rows if r["answered_by"] == "narrow request"]),
              "sampled genomes, full text": papers.summarize([r for r in out_rows if r["kind"] == "full text"]),
              "sampled genomes, abstract only": papers.summarize([r for r in out_rows if r["kind"] == "abstract"]),
              "all listed genomes": papers.summarize(out_rows, weight_key=None),
              "usd_this_run": round(spent_now, 4), "usd_these_papers": round(sum(usd.values()), 4),
              "usd per paper": {k: round(sum(v) / len(v), 4) for k, v in per_kind.items()},
              "usd_ledger_total": round(Ledger(paths.root / "spend.jsonl").total, 4)}
    paths.reports.mkdir(parents=True, exist_ok=True)
    stem = f"papers_pilot_{_now().replace(':', '')}"
    (paths.reports / f"{stem}.json").write_text(json.dumps(report, indent=1))
    if out_rows:
        with open(paths.reports / f"{stem}.tsv", "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(out_rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(out_rows)
    s = report["sampled genomes"]
    if s.get("genomes"):
        log(f"  sampled genomes compared: {s['genomes']:,}; any gain {s['any gain']['share']:.0%}; in host, source, "
            f"country or year {s['gain in host, source, country or year']['share']:.0%}; in disease, method, "
            f"outbreak or phenotype {s['gain in disease, method, outbreak or phenotype']['share']:.0%}")
    log(f"  papers answered {answered:,} of {len(requests):,}; cost of these papers ${report['usd_these_papers']:.2f} "
        f"(${report['usd_this_run']:.2f} this run); spent in this work directory ${report['usd_ledger_total']:.2f}")
    log(f"Report: {paths.reports / stem}.json (each genome in {stem}.tsv)")
    return report


def _lineage(taxa, tax_id):
    record = taxa.get(str(tax_id)) if tax_id else None
    if not record or "ambiguous" in record:
        return {}
    classification = record.get("classification") or {}
    return {rank: classification.get(rank) or {} for rank in RANKS}


def _entropy(counter):
    total = sum(counter.values())
    return -sum(n / total * math.log2(n / total) for n in counter.values() if n) if total else 0.0


def _clean(value):
    return re.sub(r"[\t\r\n]+", " ", "" if value is None else str(value)).strip()


COLUMNS = ["accession", "paired_accession", "biosample", "bioproject", "organism_name", "tax_id",
           "species", "species_taxid", "genus", "family", "order", "class", "phylum", "strain", "isolate", "serovar",
           "assembly_level", "release_date", "is_mag", "atypical", "type_material", "checkm_completeness",
           "checkm_contamination", "taxonomy_check",
           "host", "host_taxid", "host_ncbi_name", "host_ncbi_taxid", "host_evidence_field", "host_quote", "host_tier",
           "host_model", "host_confidence",
           "source_category", "mixs_extension", "niche", "source_term", "source_term_id", "body_site", "food",
           "environment", "source_evidence_field", "source_quote", "source_tier", "source_model", "source_confidence",
           "host_crosscheck", "source_crosscheck", "fable_host", "fable_source_term", "crosscheck_model",
           "country", "country_insdc", "collection_year",
           *(f"raw_{f}" for f in CATALOGUE_FIELDS), "raw_geo_loc_name", "raw_collection_date", "vocabulary_version"]


def cmd_assemble(args, paths, log):
    state = _require_harvest(paths)
    vocab = Vocabulary()
    resolver, taxa = load_resolver(paths, vocab)
    counts = distinct_keys(paths)
    labeler = make_labeler(paths, vocab, resolver, log=log)
    details = labeler.curate(sorted(counts), offline=True)
    # Fable's cross-check, from its cache only (the highest effort available per record), with the judge's
    # corrections applied.
    xc = CrossCheck(labeler, paths.crosscheck, lookup=tuple(reversed(EFFORTS)), log=lambda *_: None)
    open_recs, key_to_ck = open_records(labeler, sorted(counts))
    checked = xc.check(open_recs, offline=True) if paths.crosscheck.exists() else {}
    finals = {}

    def crosschecked(rkey):
        """(catalogue columns, cross-check columns) for one record."""
        if rkey not in finals:
            d, ck = details[rkey], key_to_ck.get(rkey)
            r = checked.get(ck) if ck else None
            if r is None:
                status = "rules only" if ck is None else "not checked"
                cross = {"host_crosscheck": status, "source_crosscheck": status}
            else:
                parts = r["parts"]
                if any(parts[p]["final"] is not parts[p]["ours"] for p in PARTS):
                    own = xc.curated(rkey, r["norm"])  # this key's own rule answers for the part not judged
                    final = {p: parts[p]["final"] if parts[p]["final"] is not parts[p]["ours"] else own[p]
                             for p in PARTS}
                    d = labeler._final(final["host"], final["source"])
                fable = {p: parts[p]["fable"][0] if parts[p]["fable"] else None for p in PARTS}
                cross = {"host_crosscheck": parts["host"]["status"], "source_crosscheck": parts["source"]["status"],
                         "fable_host": vocab.hosts[fable["host"]]["label"] if fable["host"] else "",
                         "fable_source_term": vocab.terms[fable["source"]]["label"] if fable["source"] else "",
                         "crosscheck_model": f"{r['fable_model']}, {r['fable_effort']} effort"}
            finals[rkey] = (d, cross)
        return finals[rkey]

    mags = mag_accessions(paths.ncbi)
    paths.release.mkdir(parents=True, exist_ok=True)
    table = paths.release / "isolate_catalogue.tsv.gz"
    species = {}
    totals = collections.Counter()
    tmp = Path(f"{table}.tmp")
    with gzip.open(tmp, "wt", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for a in iter_assemblies(paths.ncbi):
            attrs = canonical(a["attributes"])
            d, cross = crosschecked(record_key(a["attributes"]))
            lineage = _lineage(taxa, a.get("tax_id"))
            country, year = country_and_year(a["attributes"])
            host_taxid, host_name = resolver.exact(attrs, HOST_PRIORITY)
            row = {**{k: a.get(k) for k in ("accession", "paired_accession", "biosample", "bioproject",
                                            "organism_name", "tax_id", "strain", "isolate", "serovar",
                                            "assembly_level", "release_date", "atypical", "type_material",
                                            "checkm_completeness", "checkm_contamination", "taxonomy_check")},
                   **{rank: lineage.get(rank, {}).get("name", "") for rank in RANKS},
                   "species_taxid": lineage.get("species", {}).get("id", ""),
                   "is_mag": a["accession"] in mags, **d, **cross,
                   "host_ncbi_name": host_name, "host_ncbi_taxid": host_taxid,
                   "country": country, "country_insdc": vocab.country(country), "collection_year": year,
                   **{f"raw_{f}": attrs.get(f, "") for f in CATALOGUE_FIELDS},
                   "raw_geo_loc_name": attrs.get("geo_loc_name", ""),
                   "raw_collection_date": attrs.get("collection_date", "")}
            writer.writerow([_clean(row.get(c, "")) for c in COLUMNS])
            totals["genomes"] += 1
            totals[f"source_tier {d['source_tier']}"] += 1
            totals[f"host_tier {d['host_tier']}"] += 1
            totals[f"source_crosscheck {cross['source_crosscheck']}"] += 1
            totals[f"host_crosscheck {cross['host_crosscheck']}"] += 1
            totals["source known"] += d["source_category"] != "unknown"
            totals["host known"] += d["host"] != "unknown"
            sp_key = row["species_taxid"] or row["organism_name"]
            s = species.setdefault(sp_key, {"species": row["species"] or row["organism_name"],
                                            "species_taxid": row["species_taxid"], "genus": row["genus"],
                                            "family": row["family"], "phylum": row["phylum"],
                                            "genomes": 0, "isolates": 0, "complete": 0, "type_material": 0,
                                            "categories": collections.Counter(), "niches": collections.Counter(),
                                            "hosts": collections.Counter(), "countries": set(), "years": [],
                                            "bioprojects": set()})
            s["genomes"] += 1
            s["isolates"] += not row["is_mag"]
            s["complete"] += a.get("assembly_level") == "Complete Genome"
            s["type_material"] += bool(a.get("type_material"))
            if d["source_category"] != "unknown":
                s["categories"][d["source_category"]] += 1
                s["niches"][d["niche"]] += 1
            if d["host"] != "unknown":
                s["hosts"][d["host"]] += 1
            if row["country_insdc"]:
                s["countries"].add(row["country_insdc"])
            if year:
                s["years"].append(int(year))
            if a.get("bioproject"):
                s["bioprojects"].add(a["bioproject"])
    os.replace(tmp, table)

    def top(counter, n=5):
        return "; ".join(f"{k} {v}" for k, v in counter.most_common(n))

    summary_path = paths.release / "species_summary.tsv"
    with open(summary_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["species", "species_taxid", "genus", "family", "phylum", "genomes", "isolates",
                         "complete_genomes", "type_material", "genomes_with_source", "source_categories",
                         "source_entropy_bits", "top_source_categories", "niches", "top_niches", "host_groups",
                         "top_hosts", "countries", "first_year", "last_year", "bioprojects"])
        for s in sorted(species.values(), key=lambda s: (-s["genomes"], s["species"])):
            writer.writerow([_clean(v) for v in (
                s["species"], s["species_taxid"], s["genus"], s["family"], s["phylum"], s["genomes"], s["isolates"],
                s["complete"], s["type_material"], sum(s["categories"].values()), len(s["categories"]),
                f"{_entropy(s['categories']):.3f}", top(s["categories"]), len(s["niches"]), top(s["niches"]),
                len(s["hosts"]), top(s["hosts"]), len(s["countries"]), min(s["years"], default=""),
                max(s["years"], default=""), len(s["bioprojects"]))])
    ledger = Ledger(paths.root / "spend.jsonl")
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True,
                            cwd=Path(__file__).resolve().parent).stdout.strip()
    manifest = {"catalogue_version": CATALOGUE_VERSION, "assembled": _now(), "code_commit": commit or None,
                "ncbi_harvest": {k: state.get(k) for k in ("query", "started", "finished", "records", "mags", "api")},
                "vocabulary_version": vocab.version, "rules_revision": vocab.rules_revision,
                "models": {**{t: c["model"] for t, c in TIERS.items()}, **({"crosscheck": FABLE} if checked else {})},
                "genomes": totals["genomes"], "species": len(species), "usd_spent": round(ledger.total, 2),
                "counts": dict(sorted(totals.items()))}
    (paths.release / "manifest.json").write_text(json.dumps(manifest, indent=1))
    g = max(totals["genomes"], 1)
    log(f"Catalogue: {totals['genomes']:,} genomes, {len(species):,} species -> {table}")
    log(f"  source known for {totals['source known']:,} ({totals['source known'] / g:.0%}), host known for "
        f"{totals['host known']:,} ({totals['host known'] / g:.0%})")
    log("  source tiers: " + ", ".join(f"{k.split(' ', 1)[1]} {v:,}" for k, v in sorted(totals.items())
                                        if k.startswith("source_tier")))
    if checked:
        log("  source cross-check: " + ", ".join(f"{k.split(' ', 1)[1]} {v:,}" for k, v in sorted(totals.items())
                                              if k.startswith("source_crosscheck")))
    log(f"  per-species summary -> {summary_path}")
    return manifest


def cmd_status(args, paths, log):
    state = json.loads((paths.ncbi / "harvest_state.json").read_text()) if (paths.ncbi / "harvest_state.json").exists() else {}
    log(f"Harvest: {state.get('records', 0):,} assemblies, {'complete' if state.get('done') else 'not complete'}"
        + (f" ({state.get('finished')})" if state.get("done") else ""))
    if paths.taxa.exists():
        log(f"Taxonomy cache: {len(json.loads(paths.taxa.read_text())):,} lookups")
    if paths.ontology.exists():
        log(f"Ontology cache: {len(json.loads(paths.ontology.read_text())):,} ENVO/UBERON/FoodOn/PO terms")
    cache = AnswerCache(paths.cache)
    log(f"Model answers cached: {len(cache.entries):,} distinct records")
    if paths.crosscheck.exists():
        log(f"Cross-check answers and judge decisions cached: {len(AnswerCache(paths.crosscheck).entries):,}")
    batches = {}
    for name in ("batches.json", "crosscheck_batches.json"):
        if (paths.root / name).exists():
            batches.update(json.loads((paths.root / name).read_text()))
    for key, stage in batches.items():
        done = "finished " + stage["finished"] if stage.get("finished") else "waiting"
        log(f"  stage {stage['name']}: {stage['requests']} requests, batches {', '.join(stage['batches'])} ({done})")
    ledger = Ledger(paths.root / "spend.jsonl")
    log(f"Spent so far: ${ledger.total:.2f} ({len(ledger.entries)} ledger entries)")
    load_env()
    if os.environ.get("ANTHROPIC_API_KEY") or os.environ.get("ANTHROPIC_AUTH_TOKEN"):
        # Batches the API knows about but batches.json does not, e.g. a creation whose reply was lost.
        import anthropic
        tracked = {b for stage in batches.values() for b in stage["batches"]}
        try:
            recent = list(anthropic.Anthropic().messages.batches.list(limit=50))
        except anthropic.APIError as e:
            log(f"Could not list batches from the API ({type(e).__name__})")
            return
        untracked = [b for b in recent if b.id not in tracked]
        log(f"API: {len(recent)} recent batch(es), {len(untracked)} not tracked by this work directory")
        for b in untracked:
            log(f"  untracked {b.id} created {b.created_at:%Y-%m-%d %H:%M} {b.processing_status}, "
                f"{b.request_counts.succeeded} succeeded")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", default="catalogue_work", help="where everything is kept (default ./catalogue_work)")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("harvest", help="download NCBI assembly reports and BioSample attributes (free)")
    p.add_argument("--taxon", default="2", help="NCBI taxon (default 2, Bacteria)")
    p.add_argument("--levels", default="complete_genome,chromosome", help="assembly levels, comma-separated")
    p.add_argument("--restart", action="store_true", help="start the harvest over")
    sub.add_parser("prepare", help="taxonomy and ontology lookups, distinct records, rules tier and cost estimate "
                                   "(free)")
    p = sub.add_parser("run", help="send the records the rules leave open to the model tiers (costs money)")
    p.add_argument("--budget", type=float, required=True, help="hard cap in USD for everything spent in this workdir")
    p.add_argument("--sample", type=int, help="pilot: only this many randomly chosen open records")
    p.add_argument("--seed", type=int, default=1, help="random seed for --sample (default 1)")
    p.add_argument("--sync", action="store_true", help="small tests only: call the API directly at full price")
    p.add_argument("--poll", type=int, default=60, help="seconds between batch status checks (default 60)")
    p = sub.add_parser("verify", help=f"cross-check the model labels with {FABLE} and judge disagreements "
                                      "(costs money; run after `run`)")
    p.add_argument("--budget", type=float, required=True, help="hard cap in USD for everything spent in this workdir")
    p.add_argument("--effort", choices=EFFORTS, default="low",
                   help="effort for Fable's own labels (default low); the judge always runs at high")
    p.add_argument("--sample", type=int, help="pilot: only this many randomly chosen open records")
    p.add_argument("--seed", type=int, default=1, help="random seed for --sample (default 1)")
    p.add_argument("--poll", type=int, default=60, help="seconds between batch status checks (default 60)")
    p = sub.add_parser("rejudge", help="a second model judges a sample of the disputes Fable's judge settled, to "
                                       "check it for self-preference (costs money; run after `verify`)")
    p.add_argument("--budget", type=float, required=True, help="hard cap in USD for everything spent in this workdir")
    p.add_argument("--model", default="claude-opus-5-5", help="the second judge (default claude-opus-5-5)")
    p.add_argument("--sample", type=int, default=300, help="records to judge again (default 300)")
    p.add_argument("--seed", type=int, default=1, help="random seed for --sample (default 1)")
    p.add_argument("--poll", type=int, default=60, help="seconds between batch status checks (default 60)")
    p = sub.add_parser("papers-pilot", help="read a sample of the papers linked to the genomes with Claude Fable 5.2 "
                                            "and measure what they add (costs money; run after `assemble`)")
    p.add_argument("--budget", type=float, required=True, help="hard cap in USD for everything spent in this workdir")
    p.add_argument("--sample", type=int, default=200, help="distinct papers to read (default 200)")
    p.add_argument("--seed", type=int, default=1, help="random seed for the genome order (default 1)")
    p.add_argument("--poll", type=int, default=60, help="seconds between batch status checks (default 60)")
    p.add_argument("--narrow-refused", action="store_true",
                   help="ask the papers still refused once more, for host, source, country and year only, from their "
                        "methods sections and tables")
    p.add_argument("--retry-model", help="send the papers Fable refused or cut off once more with this model "
                                          "(e.g. claude-opus-5-5)")
    p.add_argument("--dry-run", action="store_true",
                   help="sample, fetch the texts and price the requests without calling the Claude API (free)")
    sub.add_parser("assemble", help="write the catalogue table, per-species summary and manifest")
    sub.add_parser("status", help="progress and spend")
    args = parser.parse_args(argv)
    foreign = foreign_modules()
    if foreign:
        parser.error(f"steps/ modules were loaded from another checkout of the repository ({', '.join(foreign)}). "
                     f"Run this from its own checkout's root: {Path(__file__).resolve().parents[1]}")
    paths = Paths(args.workdir)

    def log(message):
        print(message, flush=True)

    commands = {"harvest": cmd_harvest, "prepare": cmd_prepare, "run": cmd_run, "verify": cmd_verify,
                "rejudge": cmd_rejudge, "papers-pilot": cmd_papers_pilot, "assemble": cmd_assemble,
                "status": cmd_status}
    return commands[args.command](args, paths, log)


if __name__ == "__main__":
    main()
