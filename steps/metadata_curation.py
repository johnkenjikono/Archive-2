"""
Tiered, evidence-linked curation of isolate source metadata.

Adapted from the Microbiome Repo design (Olm lab): rules first, models only where they are
needed, answers limited to a fixed vocabulary, and every value tied to the text it came from.
For each distinct NCBI BioSample record (the fields niche_metadata.py reads):

  1. Rules, no model: missing values dropped, country and year parsed, host and source mapped
     straight to the vocabulary when a field equals a known synonym, or names a vocabulary
     ontology id whose label agrees with the id.
  2. Haiku maps each remaining record to the vocabulary and quotes the words it used.
  3. Records Haiku is unsure of, or whose quote is not in the record, go to two independent
     Sonnet runs; a host or source is kept when both runs agree and both quotes check out.
  4. Opus settles the Sonnet disagreements. A value it cannot back with a quote is withheld
     as unknown rather than guessed.

The vocabulary is steps/metadata_vocabulary.json: hosts with NCBI taxon ids, and source terms
with ENVO / UBERON / PO / FoodOn ids, each grouped into a fixed niche and source category so
labels mean the same thing for every species. Model answers are cached per distinct record
(keyed on the normalized record, the vocabulary version and the models), so a string such as
"blood" or "soil" is classified once across all species.

Writes the usual niche outputs (see niche_metadata.py) with extra columns in genome_metadata.csv
(source category, MIxS extension, ontology term and id, body site / food / environment, evidence
field, quote, tier, model and confidence for host and source, vocabulary version) and an extra
source_category row in ecotype_niche_association.csv.

  python steps/metadata_curation.py ecosim_output_dir "Listeria monocytogenes"   # needs ANTHROPIC_API_KEY
  python steps/metadata_curation.py --dry-run summary.jsonl [...]                # rules only, no API calls
  python steps/metadata_curation.py --verify-vocabulary                          # check ids against EBI OLS
"""
import argparse
import csv
import hashlib
import json
import os
import re
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

try:
    from steps.niche_metadata import LLM_FIELDS, association, run_niche_analysis, write_crosstab
except ImportError:  # run as a script: python steps/metadata_curation.py
    from niche_metadata import LLM_FIELDS, association, run_niche_analysis, write_crosstab

VOCAB_PATH = Path(__file__).with_name("metadata_vocabulary.json")
CACHE_NAME = "metadata_curation_cache.json"

# The one place to change the model behind each tier. fallbacks: server-side refusal fallback
# (Claude API; not offered for Haiku). chunk: distinct records per request.
TIERS = {
    "haiku": {"model": "claude-haiku-5-5", "effort": "low", "chunk": 60, "fallbacks": False},
    "sonnet": {"model": "claude-sonnet-5-5", "effort": "low", "chunk": 40, "fallbacks": True},
    "opus": {"model": "claude-opus-5-5", "effort": "medium", "chunk": 20, "fallbacks": True},
}
# USD per million input / output tokens (list prices, Nov 2026), for the cost line printed after a run.
PRICES = {"claude-haiku-5-5": (0.10, 0.50), "claude-sonnet-5-5": (2.00, 10.00), "claude-opus-5-5": (4.00, 20.00),
          "claude-fable-5-2": (8.00, 40.00)}
# Cache reads as a fraction of the input price (0.025x on Fable 5.2, 0.05x on Opus 5.5 and Sonnet 5.5, 0.1x
# elsewhere); cache writes cost 1.25x (5-minute entries) or 2x (1-hour entries). The Batch API halves every rate.
CACHE_READ = {"claude-haiku-5-5": 0.10, "claude-sonnet-5-5": 0.05, "claude-opus-5-5": 0.05, "claude-fable-5-2": 0.025}
BATCH_DISCOUNT = 0.5
MAX_TOKENS = 16000
FALLBACK_BETA = "server-side-fallback-2026-07-01"

# Where a source can come from, in priority order: the first informative field decides.
SOURCE_FIELDS = ("isolation_source", "host_body_site", "host_tissue_sampled", "env_medium",
                 "sample_type", "source_type", "env_local_scale", "env_broad_scale")
HOST_FIELDS = ("host",)
# Host-field values that name no particular host ("animal", "plant", "soil"): on a record the rules cannot read,
# they do not settle the host, because the model reading the record may find a more specific one.
VAGUE_HOSTS = ("other_animal", "other_plant", "no_host")
# Engine tiers, not calibrated probabilities.
CONFIDENCE = {"rule": "high", "haiku": "medium", "sonnet x2": "medium", "opus": "low",
              "withheld": "none", "no evidence": "none"}
CURATION_COLUMNS = ("source_category", "mixs_extension", "source_term", "source_term_id", "body_site",
                    "food", "environment", "source_evidence_field", "source_quote", "source_tier",
                    "source_model", "source_confidence", "host_taxid", "host_evidence_field", "host_quote",
                    "host_tier", "host_model", "host_confidence", "country_insdc", "vocabulary_version")
_ID_PREFIXES = {"envo": "ENVO", "uberon": "UBERON", "foodon": "FOODON", "po": "PO", "ncbitaxon": "NCBITaxon"}
_BRACKETED_ID = re.compile(r"^(.*?)[\[(]?\s*\b(envo|uberon|foodon|po|ncbitaxon)[:_](\d+)\s*[\])]?\s*$", re.I)


def normalize(value):
    """Lowercase, collapse whitespace, trim surrounding quotes and trailing punctuation."""
    return re.sub(r"\s+", " ", str(value or "")).strip().strip("\"'").rstrip(".;,").strip().lower()


def _squash(text):
    return re.sub(r"\s+", " ", str(text or "")).strip().lower()


class Vocabulary:
    def __init__(self, path=VOCAB_PATH):
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        self.data = data
        self.version = data["version"]  # the answer space; cached model answers are keyed on it
        self.rules_revision = data.get("rules_revision", 0)  # synonyms and skip lists only
        self.terms = {t["key"]: t for t in data["terms"]}
        self.hosts = {h["key"]: h for h in data["hosts"]}
        self.niches = data["niches"]
        self.categories = data["source_categories"]
        self.missing = {normalize(v) for v in data["missing_values"]}
        self.uninformative = {normalize(v) for v in data["uninformative_values"]}
        self.countries = {c.lower(): c for c in data["countries"] + data["historical_countries"]}
        self.term_by_id = {t["id"].upper(): k for k, t in self.terms.items() if t["id"]}
        self.term_synonyms, self.host_synonyms = {}, {}
        for key, term in self.terms.items():
            # Ontology labels are not synonyms: "fruit" is PO's plant part, but in a record it usually means food.
            for s in (term["label"], *term["synonyms"]):
                if normalize(s) and key != "unknown":
                    self.term_synonyms.setdefault(normalize(s), key)
        for key, host in self.hosts.items():
            for s in (host["ncbi_name"], *host["synonyms"]):
                if normalize(s):
                    self.host_synonyms.setdefault(normalize(s), key)
        self.host_by_taxid = {str(h["taxid"]): k for k, h in self.hosts.items() if h["taxid"]}

    def problems(self):
        """Internal consistency problems (no network): an empty list means the file is usable."""
        out = []
        for key, term in self.terms.items():
            if term["niche"] not in self.niches:
                out.append(f"term {key}: unknown niche {term['niche']}")
            if term["kind"] not in ("body_site", "plant_part", "food", "environment", "laboratory", "unknown"):
                out.append(f"term {key}: unknown kind {term['kind']}")
            if term["kind"] != "body_site" and term.get("category") not in self.categories:
                out.append(f"term {key}: unknown category {term.get('category')}")
            if term["id"] and not re.fullmatch(r"(ENVO|UBERON|FOODON|PO):\d+", term["id"]):
                out.append(f"term {key}: malformed id {term['id']}")
            if term["id"] and not term["ontology_label"]:
                out.append(f"term {key}: id without its verified ontology_label")
        for key, host in self.hosts.items():
            if host["kind"] not in ("human", "animal", "plant", "fungus", "none", "unknown"):
                out.append(f"host {key}: unknown kind {host['kind']}")
        for kind, entries, names in (("term", self.terms, lambda t: (t["label"], *t["synonyms"])),
                                     ("host", self.hosts, lambda h: (h["ncbi_name"], *h["synonyms"]))):
            seen = {}
            for key, entry in entries.items():
                for s in names(entry):
                    if normalize(s) and seen.setdefault(normalize(s), key) != key:
                        out.append(f"{kind} name {s!r} maps to both {seen[normalize(s)]} and {key}")
        ids = Counter(t["id"] for t in self.terms.values() if t["id"])
        out += [f"id {i} used by {n} terms" for i, n in ids.items() if n > 1]
        return out

    def informative(self, value):
        n = normalize(value)
        return bool(n) and n not in self.missing and n not in self.uninformative

    def match_term(self, value):
        """Vocabulary term for one field value, or None. A bracketed ontology id counts only when its
        label agrees: submitters attach wrong ids (e.g. "soil [ENVO:00002007]", which is sediment)."""
        if "|" in value:
            return None
        bracketed = _BRACKETED_ID.match(value.strip())
        if bracketed:
            label = normalize(bracketed.group(1))
            key = self.term_by_id.get(f"{_ID_PREFIXES[bracketed.group(2).lower()]}:{bracketed.group(3)}".upper())
            if key:
                term = self.terms[key]
                names = {normalize(term["ontology_label"]), *(normalize(s) for s in (term["label"], *term["synonyms"]))}
                return key if not label or label in names else None
            return self.term_synonyms.get(label) if label else None
        return self.term_synonyms.get(normalize(value))

    def match_host(self, value, loose=False):
        """Host key for a value naming a host. loose (for the host field only) also tries a name in
        parentheses and the binomial of a longer name ("capra hircus coreanae" -> goat)."""
        n = normalize(value)
        if n in self.host_synonyms or not loose:
            return self.host_synonyms.get(n)
        inner = re.search(r"\(([^)]*)\)", n)
        if inner and normalize(inner.group(1)) in self.host_synonyms:
            return self.host_synonyms[normalize(inner.group(1))]
        return self.host_synonyms.get(" ".join(n.split()[:2])) if len(n.split()) > 2 else None

    def category_of(self, term_key, host_key):
        """(source category key, MIxS extension) for a term; body sites take theirs from the host."""
        term, host = self.terms[term_key], self.hosts[host_key]
        if term["kind"] != "body_site":
            return term["category"], term.get("mixs") or ""
        if host["kind"] == "human":
            return "human_associated", term.get("human_mixs", "HumanAssociated")
        if host["kind"] == "animal":
            return "animal_associated", "HostAssociated"
        if host["kind"] == "plant":
            return "plant_associated", "PlantAssociated"
        return "host_associated", "HostAssociated"

    def country(self, country):
        return self.countries.get(str(country).lower(), "")

    def prompt_listing(self):
        hosts = "\n".join(f"- {k}: {h['label']}" + (f" ({h['ncbi_name']})" if h["ncbi_name"] else "")
                          for k, h in self.hosts.items())
        terms = "\n".join(f"- {k}: {t['label']} [{self.niches[t['niche']]}]" for k, t in self.terms.items())
        return f"Hosts:\n{hosts}\n\nSources (key: meaning [group]):\n{terms}"


class Answer:
    """One curated value: a vocabulary key plus the evidence and tier that produced it."""
    __slots__ = ("key", "field", "quote", "tier", "model")

    def __init__(self, key, field="", quote="", tier="no evidence", model=""):
        self.key, self.field, self.quote, self.tier, self.model = key, field, quote, tier, model

    def as_list(self):
        return [self.key, self.field, self.quote, self.tier, self.model]


def usage_cost(model, usage, batch=False):
    """USD for one response's usage (an SDK Usage object or a dict with the same fields)."""
    def get(obj, name):
        value = obj.get(name) if isinstance(obj, dict) else getattr(obj, name, None)
        return value or 0
    if usage is None or model not in PRICES:
        return 0.0
    price_in, price_out = PRICES[model]
    written = get(usage, "cache_creation_input_tokens")
    breakdown = get(usage, "cache_creation")
    long_lived = get(breakdown, "ephemeral_1h_input_tokens") if breakdown else 0
    usd = (get(usage, "input_tokens") * price_in
           + (written - long_lived) * price_in * 1.25 + long_lived * price_in * 2.0
           + get(usage, "cache_read_input_tokens") * price_in * CACHE_READ.get(model, 0.1)
           + get(usage, "output_tokens") * price_out) / 1e6
    return usd * (BATCH_DISCOUNT if batch else 1.0)


def quote_ok(fields, field, quote, value):
    """A non-unknown answer must quote words that appear in the field it names."""
    if value == "unknown":
        return True
    return bool(field in fields and _squash(quote) and _squash(quote) in _squash(fields[field]))


class AnswerCache:
    """Model answers per distinct record, shared across species. Saving merges with what is on disk,
    so parallel runs only ever lose entries (to be classified again), never corrupt the file."""

    def __init__(self, path):
        self.path = Path(path) if path else None
        self.entries = self._read()
        self.new = {}

    def _read(self):
        try:
            with open(self.path, encoding="utf-8") as f:
                return json.load(f).get("entries", {})
        except (TypeError, FileNotFoundError, json.JSONDecodeError):
            return {}

    def get(self, key):
        return self.entries.get(key)

    def put(self, key, value):
        self.entries[key] = self.new[key] = value

    def save(self):
        if not self.path or not self.new:
            return
        merged = {**self._read(), **self.new}
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(f".tmp{os.getpid()}")
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump({"entries": merged}, f, ensure_ascii=False, indent=0, sort_keys=True)
        os.replace(tmp, self.path)
        self.entries, self.new = merged, {}


class CuratedLabeler:
    """niche_metadata labeler: called with (species, keys) and returns {key: (host, niche)}; the full
    curated record for each key is kept in .details for the extra genome_metadata.csv columns.

    host_resolver (optional) maps host names and taxids the vocabulary does not list onto a vocabulary
    host through the NCBI lineage (catalogue/hosts.py); it needs .vocab_key_for_name(value) and
    .vocab_key_for_taxid(taxid), each returning a host key or None. transport (optional) runs the model
    tiers some other way, e.g. through the Message Batches API (catalogue/batch.py); it needs
    .call_many(labeler, jobs) returning one {rid: answer} dict per (tier, records, disputes) job.
    source_fields and host_fields set where the rules look, in priority order (the catalogue adds
    non-standard BioSample fields such as habitat or host_common_name)."""
    name = "tiered curation (rules, Haiku, 2x Sonnet, Opus)"

    def __init__(self, vocab=None, client=None, cache_path=None, tiers=None, fields=LLM_FIELDS, log=print,
                 host_resolver=None, transport=None, source_fields=SOURCE_FIELDS, host_fields=HOST_FIELDS,
                 source_from_host=False):
        self.vocab = vocab or Vocabulary()
        self.client = client
        self.tiers = tiers or TIERS
        self.fields = fields
        self.source_fields = source_fields
        self.host_fields = host_fields
        # A record that names a host but no source at all: label the source "host sample, site not stated"
        # (or "plant, part not stated") from the host field, instead of unknown.
        self.source_from_host = source_from_host
        self.cache = AnswerCache(cache_path)
        self.log = log
        self.host_resolver = host_resolver
        self.transport = transport
        self.details = {}
        self.usage = Counter()
        self.spend = Counter()  # USD per model, from each response's usage
        self.stats = Counter()
        self.uncacheable = set()

    def __call__(self, species, keys):
        self.curate(keys)
        return {k: (d["host"], d["niche"]) for k, d in self.details.items()}

    # --- rules -------------------------------------------------------------------------------
    def rule_source(self, fields):
        for f in self.source_fields:
            value = fields.get(f, "")
            if not self.vocab.informative(value):
                continue
            key = self.vocab.match_term(value)
            if key:
                return Answer(key, f, value.strip(), "rule")
            host = self.vocab.match_host(value)  # e.g. isolation_source "human": a host, site not stated
            if host and self.vocab.hosts[host]["kind"] in ("human", "animal"):
                return Answer("host_sample", f, value.strip(), "rule")
            if host and self.vocab.hosts[host]["kind"] == "plant":
                return Answer("plant_unspecified", f, value.strip(), "rule")
            return None  # the first informative field is not understood: a model reads it
        return Answer("unknown")

    def rule_host(self, fields, source):
        taxid = fields.get("host_taxid", "").strip()
        if taxid in self.vocab.host_by_taxid:
            return Answer(self.vocab.host_by_taxid[taxid], "host_taxid", taxid, "rule")
        if taxid and self.host_resolver:
            key = self.host_resolver.vocab_key_for_taxid(taxid)
            if key:
                return Answer(key, "host_taxid", taxid, "rule")
        for hf in self.host_fields:  # the first informative host field decides
            if self.vocab.informative(fields.get(hf, "")):
                key = self.vocab.match_host(fields[hf], loose=True)
                if key in VAGUE_HOSTS and source is None:  # host "plant", "animal" or "soil" on a sample the rules
                    return None  # cannot read: the model reading it may find the host ("isolated from cucumber")
                if not key and self.host_resolver:  # a taxon the vocabulary does not list: its NCBI lineage decides
                    key = self.host_resolver.vocab_key_for_name(fields[hf])
                return Answer(key, hf, fields[hf].strip(), "rule") if key else None
        if source is None:  # a model reads the sample's description anyway ("Rectum of pig"), so it also names
            return None  # the host, rather than a generic word in another field (source_type "animal")
        for f in self.source_fields:  # no host field: a field that is just a host name ("human", "bovine")
            key = self.vocab.match_host(fields.get(f, "")) if self.vocab.informative(fields.get(f, "")) else None
            if key:
                return Answer(key, f, fields[f].strip(), "rule")
        term = self.vocab.terms[source.key]
        if term["kind"] in ("food", "environment") and term.get("category") != "plant_associated":
            return Answer("no_host", source.field, source.quote, "rule")  # food, soil, water: no living host
        return Answer("unknown")  # a body site, a plant part or a lab strain: the host is not stated

    # --- model tiers -------------------------------------------------------------------------
    def _schema(self, record_fields):
        field_enum = [*record_fields, ""]
        item = {"type": "object", "additionalProperties": False,
                "required": ["id", "host", "host_field", "host_quote", "source", "source_field", "source_quote",
                             "certain"],
                "properties": {"id": {"type": "integer"},
                               "host": {"type": "string", "enum": list(self.vocab.hosts)},
                               "host_field": {"type": "string", "enum": field_enum},
                               "host_quote": {"type": "string"},
                               "source": {"type": "string", "enum": list(self.vocab.terms)},
                               "source_field": {"type": "string", "enum": field_enum},
                               "source_quote": {"type": "string"},
                               "certain": {"type": "boolean"}}}
        return {"type": "object", "additionalProperties": False, "required": ["labels"],
                "properties": {"labels": {"type": "array", "items": item}}}

    def _system(self):
        return (
            "You map NCBI BioSample metadata of bacterial isolates onto a fixed vocabulary. Each record is one "
            "distinct combination of BioSample fields, written as field = value.\n\n"
            "For each record return:\n"
            "- host: the host the isolate was taken from. Name a species only when a field names it; use a group "
            "entry (bird, other_mammal, other_animal, other_plant, ...) when the text names only a group; no_host "
            "when the sample was not taken from a living host (food, soil, water, a building); unknown when the "
            "fields do not say, including laboratory strains of unstated origin.\n"
            "- source: the one source term that best describes the material or body site sampled. Prefer the most "
            "specific term the text supports; use a 'not stated' term when the text gives only the broad kind; "
            "unknown when the fields say nothing about the source.\n"
            "- host_field, source_field: the field each answer rests on; host_quote, source_quote: the exact words "
            "from that field's value that support it, copied verbatim (they are checked against the record). Leave "
            "field and quote empty only for unknown.\n"
            "- certain: false when the record is ambiguous or contradictory, or the answer needs knowledge beyond "
            "the text (for example a local food name you are not sure of).\n\n"
            "Answer from the record text only. A disease name alone does not establish a body site, and where a "
            "species usually lives is not evidence about this isolate. Ignore values that describe how the isolate "
            "was grown, such as 'whole organism' or 'pure culture'.\n\n" + self.vocab.prompt_listing())

    def _record_text(self, rid, fields, disputes=None):
        lines = [f"Record {rid}:"] + [f"  {k} = {v}" for k, v in fields.items()]
        for part, (a, b) in (disputes or {}).items():
            def show(ans):
                return "no valid answer" if not ans else (
                    f"{ans[part]} (field {ans[part + '_field'] or '-'}, quote \"{ans[part + '_quote']}\")")
            lines.append(f"  disputed {part}: reviewer A says {show(a)}; reviewer B says {show(b)}")
        return "\n".join(lines)

    def _call(self, tier, records, disputes=None):
        """Ask one tier about records ({rid: fields}); returns {rid: answer dict}. A chunk that hits
        max_tokens is split; a refused or unparsable chunk returns no answers (escalated by the caller)."""
        cfg = self.tiers[tier]
        out = {}
        ids = list(records)
        for start in range(0, len(ids), cfg["chunk"]):
            out.update(self._call_chunk(tier, {i: records[i] for i in ids[start:start + cfg["chunk"]]}, disputes))
        return out

    def _call_many(self, jobs):
        """Run (tier, records, disputes) jobs, through the transport when one is set (a batch transport can
        send independent jobs, such as the two Sonnet runs, together). Returns one answer dict per job."""
        if self.transport is not None:
            return self.transport.call_many(self, jobs)
        return [self._call(tier, records, disputes) for tier, records, disputes in jobs]

    def build_request(self, tier, records, disputes=None, batch=False):
        """Messages API parameters for one chunk of records ({rid: fields}). Batch requests mark the
        system prompt for the 1-hour cache, since a batch can spread its requests over a long window."""
        cfg = self.tiers[tier]
        record_fields = sorted({f for fields in records.values() for f in fields})
        intro = ("Two independent reviewers disagreed on the parts marked 'disputed'. Decide those parts from the "
                 "record text alone, or answer unknown when the text does not clearly support either.\n\n"
                 if disputes else "")
        user = intro + "\n\n".join(self._record_text(i, f, (disputes or {}).get(i)) for i, f in records.items())
        cache = {"type": "ephemeral", "ttl": "1h"} if batch else {"type": "ephemeral"}
        return dict(model=cfg["model"], max_tokens=MAX_TOKENS,
                    system=[{"type": "text", "text": self._system(), "cache_control": cache}],
                    messages=[{"role": "user", "content": user}],
                    output_config={"effort": cfg["effort"],
                                   "format": {"type": "json_schema", "schema": self._schema(record_fields)}})

    def parse_response(self, tier, records, response):
        """(answers, status) for one chunk's response. status is "ok", "max_tokens" (split and retry),
        or why nothing usable came back (the records are then escalated by the caller)."""
        if response.stop_reason == "max_tokens":
            return None, "max_tokens"
        if response.stop_reason != "end_turn":
            return None, str(response.stop_reason)
        try:
            text = next(b.text for b in response.content if b.type == "text")
            items = json.loads(text)["labels"]
        except (StopIteration, ValueError, KeyError, TypeError):
            return None, "unparsable"
        answers = {item["id"]: item for item in items if isinstance(item, dict) and item.get("id") in records}
        model = getattr(response, "model", None)  # who answered: a refusal fallback can differ from the tier's model
        if isinstance(model, str) and model:
            for item in answers.values():
                item["_model"] = model
        return answers, "ok"

    def record_usage(self, model, usage, batch=False):
        self.usage[(model, "in")] += sum((getattr(usage, f, 0) or 0) for f in (
            "input_tokens", "cache_read_input_tokens", "cache_creation_input_tokens"))
        self.usage[(model, "out")] += getattr(usage, "output_tokens", 0) or 0
        self.spend[model] += usage_cost(model, usage, batch=batch)

    def _call_chunk(self, tier, records, disputes):
        cfg = self.tiers[tier]
        kwargs = self.build_request(tier, records, disputes)
        if self.client is None:
            import anthropic  # only needed when a model tier runs
            self.client = anthropic.Anthropic()
        if cfg.get("fallbacks"):
            response = self.client.beta.messages.create(betas=[FALLBACK_BETA], fallbacks="default", **kwargs)
        else:
            response = self.client.messages.create(**kwargs)
        # A fallback can bill a second model; response.model is the one that finished.
        self.record_usage(getattr(response, "model", None) or cfg["model"], getattr(response, "usage", None))
        answers, status = self.parse_response(tier, records, response)
        if status == "max_tokens" and len(records) > 1:
            ids = list(records)
            half = len(ids) // 2
            return {**self._call_chunk(tier, {i: records[i] for i in ids[:half]}, disputes),
                    **self._call_chunk(tier, {i: records[i] for i in ids[half:]}, disputes)}
        if status != "ok":
            self.log(f"  {cfg['model']} returned no usable answer ({status}) for {len(records)} record(s); "
                     "escalating them")
            return {}
        return answers

    def _valid(self, fields, ans, part):
        if not ans or ans.get(part) not in (self.vocab.hosts if part == "host" else self.vocab.terms):
            return False
        return quote_ok(fields, ans[f"{part}_field"], ans[f"{part}_quote"], ans[part])

    def _cache_key(self, fields):
        payload = [sorted(fields.items()), self.vocab.version, {t: c["model"] for t, c in self.tiers.items()}]
        return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()

    def run_models(self, pending):
        """pending: {rid: normalized fields}. Returns {rid: {"host": Answer, "source": Answer}}."""
        done = {}
        if not pending:
            return done
        haiku, = self._call_many([("haiku", pending, None)])
        to_sonnet = {}
        for rid, fields in pending.items():
            ans = haiku.get(rid)
            if ans and ans["certain"] and self._valid(fields, ans, "host") and self._valid(fields, ans, "source"):
                model = ans.get("_model") or self.tiers["haiku"]["model"]
                done[rid] = {p: Answer(ans[p], ans[f"{p}_field"], ans[f"{p}_quote"], "haiku", model)
                             for p in ("host", "source")}
            else:
                to_sonnet[rid] = fields
        self.stats["haiku accepted"] += len(done)
        if to_sonnet:
            # Two independent runs, the second with the records in the opposite order.
            first, second = self._call_many([("sonnet", to_sonnet, None),
                                             ("sonnet", dict(reversed(list(to_sonnet.items()))), None)])
            disputes = {}
            for rid, fields in to_sonnet.items():
                a, b = first.get(rid), second.get(rid)
                done[rid] = {}
                for p in ("host", "source"):
                    if self._valid(fields, a, p) and self._valid(fields, b, p) and a[p] == b[p]:
                        model = a.get("_model") or self.tiers["sonnet"]["model"]
                        done[rid][p] = Answer(a[p], a[f"{p}_field"], a[f"{p}_quote"], "sonnet x2", model)
                    else:
                        disputes.setdefault(rid, {})[p] = (a if self._valid(fields, a, p) else None,
                                                           b if self._valid(fields, b, p) else None)
            self.stats["sonnet records"] += len(to_sonnet)
            if disputes:
                opus, = self._call_many([("opus", {rid: to_sonnet[rid] for rid in disputes}, disputes)])
                model = self.tiers["opus"]["model"]
                for rid, parts in disputes.items():
                    ans = opus.get(rid)
                    if ans is None:  # the request failed rather than the evidence: try again next run
                        self.uncacheable.add(rid)
                    answered_by = (ans or {}).get("_model") or model
                    for p in parts:
                        if self._valid(to_sonnet[rid], ans, p) and ans[p] != "unknown":
                            done[rid][p] = Answer(ans[p], ans[f"{p}_field"], ans[f"{p}_quote"], "opus", answered_by)
                        else:
                            done[rid][p] = Answer("unknown", tier="withheld", model=answered_by)
                self.stats["opus records"] += len(disputes)
        return done

    # --- whole step --------------------------------------------------------------------------
    def curate(self, keys, offline=False):
        """Curate every key (a tuple of (field, value) pairs, as niche_metadata builds them). offline uses
        the rules and the cache only: records no model has seen yet come back unknown ("not curated") and
        nothing is written to the cache."""
        rules, pending, members = {}, {}, {}
        for key in keys:
            fields = {k: v for k, v in key if self.vocab.informative(v)}
            source = self.rule_source(fields)
            host = self.rule_host(fields, source)
            if source is not None and host is not None:
                rules[key] = (host, source)
                continue
            norm = {k: normalize(v) for k, v in fields.items()}
            cache_key = self._cache_key(norm)
            members.setdefault(cache_key, []).append((key, host, source))
            if cache_key not in pending and not self.cache.get(cache_key):
                pending[cache_key] = norm
        self.stats["rules"] += len(rules)
        self.stats["cache hits"] += len([k for k in members if k not in pending])
        self.log(f"Curating {len(keys)} distinct record(s): {len(rules)} by rules, "
                 f"{len(members) - len(pending)} from the cache, "
                 + (f"{len(pending)} not curated yet (offline)" if offline else f"{len(pending)} for the models"))
        ids = {cache_key: i for i, cache_key in enumerate(pending)}
        answers = {} if offline else self.run_models({ids[k]: fields for k, fields in pending.items()})
        fresh = {}
        for cache_key, rid in ids.items():
            ans = answers.get(rid, {})
            missing = Answer("unknown", tier="not curated" if offline else "withheld")
            fresh[cache_key] = {p: (ans[p] if p in ans else missing).as_list() for p in ("host", "source")}
            if rid not in self.uncacheable and not offline:
                self.cache.put(cache_key, fresh[cache_key])
        self.cache.save()
        for key, (host, source) in rules.items():
            self.details[key] = self._final(host, self._host_only_source(host, source))
        for cache_key, group in members.items():
            cached = fresh.get(cache_key) or self.cache.get(cache_key) or {}
            model_host, model_source = (Answer(*cached[p]) if p in cached else Answer("unknown", tier="withheld")
                                        for p in ("host", "source"))
            for key, host, source in group:
                # A rule answer beats a model answer; "no evidence" from the rules defers to the model.
                host = host if host is not None and host.tier == "rule" else model_host
                source = source if source is not None and source.tier == "rule" else model_source
                self.details[key] = self._final(host, source)
        self.log("  " + ", ".join(f"{k}: {v}" for k, v in sorted(self.stats.items())))
        if self.usage:
            self.log(f"  API usage ≈ ${self.cost():.2f} ({self.usage_text()})")
        return self.details

    def _host_only_source(self, host, source):
        """With source_from_host: a rule-settled record that names a living host but says nothing about the
        sample is host-associated with the site not stated; the host field is its evidence."""
        if (not self.source_from_host or source.key != "unknown" or source.tier != "no evidence"
                or host.key in ("unknown", "no_host")):
            return source
        term = "plant_unspecified" if self.vocab.hosts[host.key]["kind"] == "plant" else "host_sample"
        return Answer(term, host.field, host.quote, host.tier, host.model)

    def cost(self):
        return sum(self.spend.values())

    def usage_text(self):
        models = sorted({m for m, _ in self.usage})
        return "; ".join(f"{m}: {self.usage[(m, 'in')]:,} in / {self.usage[(m, 'out')]:,} out" for m in models)

    def _final(self, host, source):
        vocab = self.vocab
        term, host_entry = vocab.terms[source.key], vocab.hosts[host.key]
        category, mixs = vocab.category_of(source.key, host.key)
        term_text = f"{term['label']} [{term['id']}]" if term["id"] else term["label"]
        return {
            "host": host_entry["label"],
            "niche": vocab.niches[term["niche"]],
            "source_category": vocab.categories[category]["label"],
            "mixs_extension": mixs,
            "source_term": term["label"],
            "source_term_id": term["id"] or "",
            "body_site": term_text if term["kind"] in ("body_site", "plant_part") else "",
            "food": term_text if term["kind"] == "food" else "",
            "environment": term_text if term["kind"] in ("environment", "laboratory") else "",
            "source_evidence_field": source.field,
            "source_quote": source.quote,
            "source_tier": source.tier,
            "source_model": source.model,
            "source_confidence": CONFIDENCE.get(source.tier, "none") if source.key != "unknown" else "none",
            "host_taxid": host_entry["taxid"] or "",
            "host_evidence_field": host.field,
            "host_quote": host.quote,
            "host_tier": host.tier,
            "host_model": host.model,
            "host_confidence": CONFIDENCE.get(host.tier, "none") if host.key != "unknown" else "none",
            "vocabulary_version": vocab.version,
        }


def _key_from_row(row, fields=LLM_FIELDS):
    """The labeler key niche_metadata built for this genome, rebuilt from its raw_* columns."""
    return tuple((k, row[f"raw_{k}"]) for k in fields if row.get(f"raw_{k}"))


def add_curation_columns(table, labeler):
    """Append the curated columns to genome_metadata.csv in place."""
    with open(table, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        header, rows = reader.fieldnames, list(reader)
    blank = labeler._final(Answer("unknown"), Answer("unknown"))
    for row in rows:
        detail = labeler.details.get(_key_from_row(row, labeler.fields), blank)
        row.update({c: detail[c] for c in CURATION_COLUMNS if c in detail})
        row["country_insdc"] = labeler.vocab.country(row.get("country", ""))
    tmp = Path(f"{table}.tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[*header, *(c for c in CURATION_COLUMNS if c not in header)])
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, table)
    return rows


def test_source_category(out, rows):
    """Add a source_category row (and its crosstab) next to niche_metadata's host / niche / country tests."""
    assoc_path = out / "ecotype_niche_association.csv"
    ecotype_of = {r["accession"]: r["ecotype"] for r in rows if r.get("ecotype")}
    if not ecotype_of or not assoc_path.exists():
        return None
    category_of = {r["accession"]: r["source_category"] for r in rows}  # "unknown" is left out of the test
    write_crosstab(out / "ecotype_x_source_category.csv", ecotype_of, category_of)
    stats = association(ecotype_of, category_of)
    with open(assoc_path, newline="", encoding="utf-8") as f:
        existing = list(csv.reader(f))
    xml_used = existing[1][1] if len(existing) > 1 else ""
    with open(assoc_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerows([r for r in existing if not r or r[0] != "source_category"])
        writer.writerow(["source_category", xml_used, *stats])
    print(f"  source_category: n={stats[0]}, V={stats[4]}, p={stats[5]}")
    return stats


def default_cache_path(ecosim_out_dir):
    """<workdir>/metadata_curation_cache.json: the folder holding every species' run (the results
    volume on Modal), so classified strings are shared across species."""
    return Path(ecosim_out_dir).resolve().parent.parent / CACHE_NAME


def run_curated_niche_analysis(ecosim_out_dir, species, cache_path=None, client=None, vocab=None, fetcher=None):
    """niche_metadata.run_niche_analysis with the tiered labeler, then the extra columns and test."""
    out = Path(ecosim_out_dir)
    labeler = CuratedLabeler(vocab=vocab, client=client, cache_path=cache_path or default_cache_path(out))
    kwargs = {"fetcher": fetcher} if fetcher else {}
    result = run_niche_analysis(out, species, labeler=labeler, **kwargs)
    table = out / "genome_metadata.csv"
    if result is None or not table.exists():
        return result
    rows = add_curation_columns(table, labeler)
    test_source_category(out, rows)
    tiers = Counter(r["source_tier"] for r in rows)
    print(f"Curated metadata for {len(rows)} genomes (source tiers: "
          + ", ".join(f"{t} {n}" for t, n in tiers.most_common()) + f"); vocabulary {labeler.vocab.version}")
    return result


# --- offline helpers ---------------------------------------------------------------------------
def attrs_from_summary(record):
    biosample = record.get("assembly_info", {}).get("biosample", {})
    return {a["name"]: a.get("value", "") for a in biosample.get("attributes", [])}


def dry_run(paths, vocab=None, log=print):
    """Rules tier only, on `datasets summary genome ... --as-json-lines` files: how much the rules settle,
    and how many distinct records would go to the models. Makes no API calls."""
    vocab = vocab or Vocabulary()
    results = {}
    for path in paths:
        labeler = CuratedLabeler(vocab=vocab, log=lambda *_: None)
        genomes = [attrs_from_summary(json.loads(line)) for line in open(path, encoding="utf-8") if line.strip()]
        keys = [tuple((k, a[k]) for k in LLM_FIELDS if a.get(k)) for a in genomes]
        by_rules, labelled, pending, unresolved = 0, 0, set(), Counter()
        source_tiers = Counter()
        for key in keys:
            fields = {k: v for k, v in key if vocab.informative(v)}
            source = labeler.rule_source(fields)
            host = labeler.rule_host(fields, source)
            if source is not None and host is not None:
                by_rules += 1
                labelled += source.key != "unknown"  # the rest have no informative field at all
                source_tiers[vocab.terms[source.key]["niche"]] += 1
            else:
                pending.add(json.dumps(sorted({k: normalize(v) for k, v in fields.items()}.items())))
                if source is None:
                    first = next((fields[f] for f in SOURCE_FIELDS if vocab.informative(fields.get(f, ""))), "")
                    unresolved[f"source: {normalize(first)}"] += 1
                if host is None and fields.get("host"):  # with no host field the model reads host and source together
                    unresolved[f"host: {normalize(fields['host'])}"] += 1
        system_tokens = len(labeler._system()) / 4
        record_tokens = sum(len(labeler._record_text(0, dict(json.loads(p)))) for p in pending) / 4
        results[Path(path).name] = {"genomes": len(keys), "distinct_records": len(set(keys)),
                                    "resolved_by_rules": by_rules, "source_labelled_by_rules": labelled,
                                    "distinct_for_models": len(pending),
                                    "system_tokens": round(system_tokens), "record_tokens": round(record_tokens),
                                    "top_unresolved": unresolved.most_common(15),
                                    "rule_niches": source_tiers.most_common(10)}
        r = results[Path(path).name]
        log(f"{Path(path).name}: {r['genomes']} genomes, {r['distinct_records']} distinct records; rules settle "
            f"{by_rules} genomes ({by_rules / max(len(keys), 1):.0%}; {labelled} with a source, the rest have no "
            f"source information); {len(pending)} distinct records for the models")
    return results


def estimate_cost(n_records, system_tokens, record_tokens, escalate=0.3, dispute=0.3, out_per_record=120):
    """Rough USD for one species: Haiku on every record, two Sonnet runs on the escalated share, Opus on
    the disputed share of those. Output allows for quotes plus some thinking at low effort."""
    def run(model, records, chunk):
        if not records:
            return 0.0
        calls = -(-records // chunk)
        tokens_in = calls * system_tokens + record_tokens * records / max(n_records, 1)
        price_in, price_out = PRICES[model]
        return tokens_in / 1e6 * price_in + records * out_per_record / 1e6 * price_out
    haiku = run(TIERS["haiku"]["model"], n_records, TIERS["haiku"]["chunk"])
    sonnet = 2 * run(TIERS["sonnet"]["model"], round(n_records * escalate), TIERS["sonnet"]["chunk"])
    opus = run(TIERS["opus"]["model"], round(n_records * escalate * dispute), TIERS["opus"]["chunk"])
    return haiku, sonnet, opus


def ols_lookup(curie, timeout=60):
    """The EBI OLS4 term records for one OBO id (empty when the id does not exist)."""
    url = (f"https://www.ebi.ac.uk/ols4/api/ontologies/{curie.split(':')[0].lower()}/terms?"
           f"{urllib.parse.urlencode({'obo_id': curie})}")
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return json.load(response).get("_embedded", {}).get("terms", [])


def verify_vocabulary(vocab=None, lookup=ols_lookup):
    """Check every ontology id and host taxon id against EBI OLS4. Returns a list of problems."""
    vocab = vocab or Vocabulary()
    checks = [(t["id"], t["ontology_label"]) for t in vocab.terms.values() if t["id"]]
    checks += [(f"NCBITaxon:{h['taxid']}", h["ncbi_name"]) for h in vocab.hosts.values() if h["taxid"]]
    problems = []
    for curie, expected in checks:
        try:
            terms = lookup(curie)
        except Exception as exc:  # noqa: BLE001 - report and keep checking the rest
            problems.append(f"{curie}: lookup failed ({exc})")
            continue
        if not terms:
            problems.append(f"{curie}: not found in OLS")
        elif terms[0]["label"].lower() != expected.lower():
            problems.append(f"{curie}: OLS label {terms[0]['label']!r}, vocabulary says {expected!r}")
        elif terms[0].get("is_obsolete"):
            problems.append(f"{curie}: obsolete")
    return problems


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ecosim_out_dir", nargs="?")
    parser.add_argument("species", nargs="?")
    parser.add_argument("--cache", help=f"answer cache (default: {CACHE_NAME} two folders above ecosim_out_dir)")
    parser.add_argument("--dry-run", nargs="+", metavar="SUMMARY_JSONL",
                        help="rules only on `datasets summary genome ... --as-json-lines` files; no API calls")
    parser.add_argument("--verify-vocabulary", action="store_true", help="check every id against EBI OLS4")
    args = parser.parse_args()
    if args.verify_vocabulary:
        found = Vocabulary().problems() + verify_vocabulary()
        print("\n".join(found) or "vocabulary OK")
        sys.exit(1 if found else 0)
    elif args.dry_run:
        for name, r in dry_run(args.dry_run).items():
            haiku, sonnet, opus = estimate_cost(r["distinct_for_models"], r["system_tokens"], r["record_tokens"])
            print(f"  {name}: est. ${haiku + sonnet + opus:.2f} (Haiku ${haiku:.3f}, Sonnet x2 ${sonnet:.3f}, "
                  f"Opus ${opus:.3f}) at 30% escalation / 30% disputes")
            print("  most common unresolved:", r["top_unresolved"][:10])
    elif args.ecosim_out_dir and args.species:
        run_curated_niche_analysis(args.ecosim_out_dir, args.species, cache_path=args.cache)
    else:
        parser.error("give ecosim_out_dir and species, --dry-run files, or --verify-vocabulary")
