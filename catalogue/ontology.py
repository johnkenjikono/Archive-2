"""
Ontology codes in BioSample source fields, placed in the curation vocabulary through the ontology itself.

Many BioSamples give the sample only as a code, alone or with its label: env_medium "ENVO:00002262",
"clay soil [ENVO:00002262]", "uberon_0004983". The vocabulary lists one id per term, so the rules understood a
code only when it was a term's own id. For any other ENVO, UBERON, FoodOn or PO code (looked up on EBI OLS4):

  1. The code's own name counts as text: if its ontology label or a synonym is a vocabulary synonym, that term
     (ocean water, ENVO:00002149 -> seawater; fecal material, ENVO:00002003 -> feces).
  2. Otherwise the code climbs its ontology, nearest first, to the closest specific term, by id or by label.
     ENVO and FoodOn codes climb is_a links only (clay soil is a soil). UBERON and PO codes also climb part_of
     links (mucosa of vagina, UBERON:0004983, is part of the vagina). In ENVO, part_of joins a material to the
     place it lies in (estuarine mud is part of an estuary), not to what it is made of. The climb stays inside
     these four ontologies: NCBI Taxonomy's top node is called "root". Two kinds of term are never reached by
     climbing:
       - A broad term ("water, type not stated", "environmental sample, type not stated"). Estuarine mud is a
         kind of environmental material, but a model reading "estuarine mud" will say sediment.
       - The terms the vocabulary lists under not_placed_by_descent (ENVO's house includes cow sheds).

A code is left to the models, as before, in four cases:
  - no term is found within MAX_STEPS steps;
  - the nearest terms tie;
  - a label given with the code is neither the code's own label nor one of its synonyms (submitters attach wrong
    ids);
  - the codes in one field lead to different terms.

`prepare` looks the codes up once and caches the ontology fragment in <workdir>/ontology_cache.json; the rules
then read only that cache.
"""
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

try:
    from steps.metadata_curation import normalize
except ImportError:  # pragma: no cover - run from another directory
    from metadata_curation import normalize

MAX_STEPS = 3  # beyond three steps up, the nearest vocabulary term is usually too general to be useful
PREFIXES = {"envo": "ENVO", "uberon": "UBERON", "foodon": "FOODON", "po": "PO"}
# OLS4 links per ontology: "parents" are is_a parents; "hierarchicalParents" add part_of (and, in ENVO, links
# such as "determined by" that lead out of the ontology).
RELATION = {"uberon": "hierarchicalParents", "po": "hierarchicalParents", "envo": "parents", "foodon": "parents"}
OLS = "https://www.ebi.ac.uk/ols4/api/ontologies"
USER_AGENT = "isolate-catalogue/0.1"
_CODE = re.compile(r"^(?P<pre>.*?)[\[(]?\s*\b(?P<prefix>envo|uberon|foodon|po)[:_](?P<num>\d+)\s*[\])]?\s*"
                   r"(?:\((?P<post>[^()]*)\))?\s*$", re.I)
_URL = re.compile(r"^https?://purl\.obolibrary\.org/obo/(envo|uberon|foodon|po)_(\d+)$", re.I)
_FAILED = object()


def parse_codes(value):
    """[(CURIE, label)] for a value made only of ENVO, UBERON, FoodOn or PO codes, each with at most a label
    ("ENVO:00002262", "clay soil [ENVO:00002262]", "ENVO_00005750 (grassland soil)", several joined by | or ;).
    [] for anything else, including codes from other ontologies."""
    out = []
    for part in re.split(r"[|;]", value):
        part = part.strip()
        if not part:
            continue
        url = _URL.match(part)
        match = None if url else _CODE.match(part)
        if url:
            out.append((f"{PREFIXES[url.group(1).lower()]}:{url.group(2)}", ""))
        elif match:
            label = (match.group("pre") or "").strip(" _-:") or (match.group("post") or "").strip()
            out.append((f"{PREFIXES[match.group('prefix').lower()]}:{match.group('num')}", label))
        else:
            return []
    return out


def _get(url, timeout):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.load(response)


def relation(curie):
    """The OLS4 parent link a code climbs (see RELATION)."""
    return RELATION.get(curie.split(":")[0].lower(), "parents")


def ols_node(curie, timeout=60):
    """One term from EBI OLS4: {"label", "synonyms", "relation", "parents": [[CURIE, label], ...]} with the
    direct parents it climbs (relation(curie)) inside ENVO, UBERON, FoodOn and PO, or None when OLS does not
    know the id."""
    ontology = curie.split(":")[0].lower()
    try:
        found = _get(f"{OLS}/{ontology}/terms?{urllib.parse.urlencode({'obo_id': curie})}", timeout)
    except urllib.error.HTTPError as e:
        if e.code == 404:  # how OLS answers for an id it does not know (mistyped, too few digits, retired)
            return None
        raise
    terms = found.get("_embedded", {}).get("terms", [])
    term = next((t for t in terms if t.get("is_defining_ontology")), terms[0] if terms else None)
    if term is None:
        return None
    parents = []
    href = term.get("_links", {}).get(relation(curie), {}).get("href")
    if href:
        listed = _get(href, timeout).get("_embedded", {}).get("terms", [])
        parents = [[t["obo_id"], t.get("label") or ""] for t in listed
                   if (t.get("obo_id") or "").split(":")[0].lower() in PREFIXES]
    return {"label": term.get("label") or "", "synonyms": term.get("synonyms") or [], "relation": relation(curie),
            "parents": parents}


def current(graph, curie):
    """A cached term, or None when it is missing, unknown to OLS, or was fetched with other parent links."""
    node = graph.get(curie)
    return node if node and node.get("relation") == relation(curie) else None


def codes_to_resolve(labeler, keys):
    """The codes the source rule may look up: those in the first informative source field of each record (the
    field that decides the source), where the vocabulary alone does not settle it."""
    codes = set()
    for key in keys:
        fields = {k: v for k, v in key if labeler.vocab.informative(v)}
        value = next((fields[f] for f in labeler.source_fields if f in fields), None)
        if value and labeler.vocab.match_term(value) is None:
            codes.update(curie for curie, _ in parse_codes(value))
    return codes


def resolve_codes(codes, path, vocab, log=print, fetch=None, max_steps=MAX_STEPS, workers=6):
    """Look up each code and its ancestors, up to max_steps above it or a vocabulary term, caching every term in
    path (an unknown id is cached as null; a failed lookup is not cached, so the next run tries again). Returns
    the cached graph {CURIE: node or None}."""
    path = Path(path)
    graph = json.loads(path.read_text()) if path.exists() else {}
    fetch = fetch or ols_node

    def safe(curie):
        try:
            return fetch(curie)
        except Exception:  # noqa: BLE001 - network trouble: not cached, tried again next run
            return _FAILED

    # Each code climbs one step at a time until it names a term itself or a step reaches one (as
    # OntologyResolver.vocab_key searches), so only the terms that search needs are looked up. Every code itself
    # is looked up, for its label and synonyms.
    resolver = OntologyResolver(vocab, graph, max_steps)
    active, fetched, failed = {c: {c} for c in codes}, 0, 0
    for step in range(max_steps):  # the terms up to max_steps - 1 above a code name those up to max_steps above
        todo = sorted({n for nodes in active.values() for n in nodes
                       if not (n in graph and (graph[n] is None or current(graph, n)))})
        if todo:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                for curie, node in zip(todo, pool.map(safe, todo)):
                    if node is _FAILED:
                        failed += 1
                    else:
                        graph[curie] = node
                        fetched += 1
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(".tmp")
            tmp.write_text(json.dumps(graph, indent=0, sort_keys=True))
            tmp.replace(path)
        climbing = {}
        for code, nodes in active.items():
            if step == 0 and resolver.own_term(code):
                continue
            parents = {(p, label) for n in nodes if current(graph, n) for p, label in graph[n]["parents"]}
            if parents and not any(resolver.ancestor_term(p, label) for p, label in parents):
                climbing[code] = {p for p, _ in parents}
        active = climbing
    log(f"Ontology codes: {len(codes):,} to place, {fetched:,} terms looked up on EBI OLS4"
        + (f", {failed:,} lookups failed (tried again next run)" if failed else "")
        + f" (cached in {path.name})")
    return graph


class OntologyResolver:
    """graph: {CURIE: {"label", "synonyms", "parents"} or None}, from resolve_codes."""

    def __init__(self, vocab, graph, max_steps=MAX_STEPS):
        self.vocab, self.graph, self.max_steps = vocab, graph, max_steps
        # Terms a code may be but never climbs to: the broad "... not stated" ones, and those the vocabulary
        # means more narrowly than the ontology does.
        self.not_by_descent = {k for k, t in vocab.terms.items() if t["label"].endswith("not stated")}
        self.not_by_descent |= set(vocab.not_placed_by_descent)

    def own_term(self, curie):
        """The term a code names itself: its id is the term's, or its label or a synonym is a vocabulary synonym
        (all of them the same term)."""
        if curie.upper() in self.vocab.term_by_id:
            return self.vocab.term_by_id[curie.upper()]
        node = current(self.graph, curie) or {}
        names = (node.get("label"), *(node.get("synonyms") or ()))
        keys = {self.vocab.term_synonyms.get(normalize(n)) for n in names if n} - {None}
        return keys.pop() if len(keys) == 1 else None

    def ancestor_term(self, curie, label):
        """The specific term an ancestor of a code is, by id or label, or None."""
        key = self.vocab.term_by_id.get(curie.upper()) or self.vocab.term_synonyms.get(normalize(label or ""))
        return None if key in self.not_by_descent else key

    def vocab_key(self, curie):
        """The term curie names itself, else the nearest specific term above it; None when none is within
        max_steps or the nearest tie."""
        return self.own_term(curie) or self.climb(curie)

    def climb(self, curie):
        """The nearest specific term above curie, or None when none is within max_steps or the nearest tie."""
        frontier, seen = [curie], {curie}
        for _ in range(self.max_steps):
            parents = []
            for c in frontier:
                for p, label in (current(self.graph, c) or {}).get("parents", ()):
                    if p not in seen:
                        seen.add(p)
                        parents.append((p, label))
            hits = {self.ancestor_term(p, label) for p, label in parents} - {None}
            if hits:
                return hits.pop() if len(hits) == 1 else None
            frontier = [p for p, _ in parents]
        return None

    def key_for_value(self, value, others=()):
        """The vocabulary term for a field value made only of ontology codes, or None (see the module doc).
        others: the record's other informative source fields. A code naming a term itself counts like any
        synonym, but climbing is an inference: it is used only when those fields hold nothing but codes too.
        Where there is free text ("orchard" coded in env_medium, plant_body_site "leaf"), a model reads the
        whole record, as before."""
        codes = parse_codes(value)
        if not codes:
            return None
        climb = all(parse_codes(v) for v in others)
        keys = set()
        for curie, label in codes:
            node = current(self.graph, curie)
            if not node:
                return None
            names = {normalize(node["label"]), *(normalize(s) for s in node.get("synonyms") or ())}
            if label and normalize(label) not in names:
                return None
            keys.add(self.own_term(curie) or (self.climb(curie) if climb else None))
        return keys.pop() if len(keys) == 1 and None not in keys else None
