"""
Independent cross-check of the catalogue's model labels with Claude Fable 5.2.

Fable labels every record the curation sent to the models again, with the same prompt, vocabulary, answer
schema and quote check, without seeing the first answers. Each part (host, source) is then compared with the
curation's final answer:

  agrees         the same host, or the same source term
  same niche     (source only) a different term in the same niche: kept as curated, not judged
  disagrees      anything else, including a label against unknown

A disagreement on a part the models decided goes to a judge (Fable at high effort) that sees the record and both
answers, in an order that alternates between records, through the curation's own "two reviewers disagreed"
prompt. The judge's answer, or unknown when it cannot back one with a quote from the record, replaces the
curated one. A disagreement on a part the rules decided is reported but not judged: rules are deterministic,
so a wrong rule is fixed in the vocabulary instead. Fable's labels and the judge's decisions are cached in
<workdir>/crosscheck_cache.json, keyed on the record, the vocabulary version, the model and the effort.

Because the judge is also Fable, it may favour Fable's answers. To check for that, a second model can judge a sample
of the settled disputes (`build.py rejudge`). It gets the same prompts, with the answers in the same order. Its
decisions are cached apart and never change the catalogue. The sample leaves out records whose curated answer came
from the second model, so that model is a third party to both answers.
"""
import csv
import hashlib
import json
import math
import random
from collections import Counter

try:
    from steps.metadata_curation import Answer, AnswerCache, CuratedLabeler, normalize
except ImportError:  # pragma: no cover - run from another directory
    from metadata_curation import Answer, AnswerCache, CuratedLabeler, normalize

FABLE = "claude-fable-5-2"
EFFORTS = ("low", "medium", "high", "xhigh")
PARTS = ("host", "source")
# Records per request: fewer at higher effort, so thinking plus answers stay well inside max_tokens.
LABEL_CHUNK = {"low": 40, "medium": 40, "high": 25, "xhigh": 20}
JUDGE_EFFORT, JUDGE_CHUNK = "high", 15
AGREEING = ("agrees", "same niche")
JUDGED = ("judge kept curation", "judge took Fable's", "judge: unknown", "judge chose another")


def tiers(effort, judge_model=FABLE):
    return {"fable": {"model": FABLE, "effort": effort, "chunk": LABEL_CHUNK[effort], "fallbacks": True},
            "judge": {"model": judge_model, "effort": JUDGE_EFFORT, "chunk": JUDGE_CHUNK, "fallbacks": True}}


def open_records(main, keys):
    """The records the rules leave to the models: ({main cache key: (key, normalized fields)},
    {key: main cache key}). Several raw keys can share one normalized record."""
    records, key_to_ck = {}, {}
    for key in keys:
        fields = {k: v for k, v in key if main.vocab.informative(v)}
        source = main.rule_source(fields)
        if source is not None and main.rule_host(fields, source) is not None:
            continue
        norm = {k: normalize(v) for k, v in fields.items()}
        ck = main._cache_key(norm)
        records.setdefault(ck, (key, norm))
        key_to_ck[key] = ck
    return records, key_to_ck


class CrossCheck:
    """main is the catalogue's CuratedLabeler: its rules and its answer cache are the curation being checked.
    effort is the effort Fable labels at; lookup lists the efforts whose cached labels count, in order of
    preference (assemble prefers the highest effort available). judge_model is Fable for the cross-check itself,
    or another model judging the same disputes a second time (see judge)."""

    def __init__(self, main, cache_path, effort="low", lookup=None, client=None, transport=None, log=print,
                 judge_model=FABLE):
        self.main, self.vocab, self.effort, self.log = main, main.vocab, effort, log
        self.lookup = tuple(lookup or (effort,))
        self.judge_model = judge_model
        self.labeler = CuratedLabeler(vocab=self.vocab, client=client, tiers=tiers(effort, judge_model),
                                      fields=main.fields, log=log, transport=transport)
        self.cache = AnswerCache(cache_path)
        self.stats = Counter()

    def _key(self, kind, norm, extra, model=FABLE):
        payload = [kind, sorted(norm.items()), self.vocab.version, model, extra]
        return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()

    def cost(self):
        return self.labeler.cost()

    # --- Fable's own labels ------------------------------------------------------------------------
    def label(self, records):
        """Ask Fable about the records ({main cache key: (key, norm)}) it has not labelled at this effort.
        Returns how many were asked."""
        todo = {ck: norm for ck, (_, norm) in records.items()
                if not self.cache.get(self._key("label", norm, self.effort))}
        self.log(f"Fable 5.2 at {self.effort} effort: {len(records) - len(todo):,} of {len(records):,} record(s) "
                 f"already labelled, {len(todo):,} to ask")
        if not todo:
            return 0
        self.stats["labelled"] += len(todo)
        ids = {ck: i for i, ck in enumerate(todo)}
        answers, = self.labeler._call_many([("fable", {ids[ck]: norm for ck, norm in todo.items()}, None)])
        for ck, norm in todo.items():
            ans = answers.get(ids[ck])
            if ans is None:  # the request failed rather than the evidence: ask again next time
                self.stats["fable gave no answer"] += 1
                continue
            entry = {"model": ans.get("_model") or FABLE, "effort": self.effort, "certain": bool(ans.get("certain"))}
            for p in PARTS:  # an answer whose quote is not in the record does not count
                entry[p] = [ans[p], ans[f"{p}_field"], ans[f"{p}_quote"]] if self.labeler._valid(norm, ans, p) else None
            self.cache.put(self._key("label", norm, self.effort), entry)
        self.cache.save()
        return len(todo)

    def fable_answer(self, norm):
        for effort in self.lookup:
            entry = self.cache.get(self._key("label", norm, effort))
            if entry:
                return entry
        return None

    def effort_agreement(self, records, other):
        """Per part, how often Fable's labels at this effort match its labels at another effort."""
        counts = {p: Counter() for p in PARTS}
        for _, norm in records.values():
            mine = self.cache.get(self._key("label", norm, self.effort))
            theirs = self.cache.get(self._key("label", norm, other))
            if mine and theirs:
                for p in PARTS:
                    same = (mine[p] or [None])[0] == (theirs[p] or [None])[0]
                    counts[p]["same" if same else "different"] += 1
        return {p: dict(c) for p, c in counts.items() if c}

    # --- comparison --------------------------------------------------------------------------------
    def curated(self, key, norm):
        """The curation's final answer per part, as the catalogue reports it, or None before `run` answered."""
        cached = self.main.cache.get(self.main._cache_key(norm))
        if not cached:
            return None
        fields = {k: v for k, v in key if self.vocab.informative(v)}
        source = self.main.rule_source(fields)
        rule = {"source": source, "host": self.main.rule_host(fields, source)}
        # As in CuratedLabeler.curate: a rule answer beats a model answer.
        return {p: rule[p] if rule[p] is not None and rule[p].tier == "rule" else Answer(*cached[p]) for p in PARTS}

    def compare(self, part, ours, theirs):
        """agrees, same niche, disagrees, or no valid answer (Fable's answer failed the quote check)."""
        if theirs is None:
            return "no valid answer"
        if ours.key == theirs[0]:
            return "agrees"
        terms = self.vocab.terms
        if (part == "source" and "unknown" not in (ours.key, theirs[0])
                and terms[ours.key]["niche"] == terms[theirs[0]]["niche"]):
            return "same niche"
        return "disagrees"

    def check(self, records, offline=False):
        """Label records with Fable (unless offline), compare with the curation and judge the model-decided
        disagreements (offline: only decisions already cached). Returns {main cache key: result}, with None
        where the curation or Fable has not answered yet. result["parts"][part] has "ours" (Answer), "fable"
        ([key, field, quote], or None), "status" and "final" (the Answer the catalogue reports)."""
        if not offline:
            self.label(records)
        results, disputed = {}, {}
        for ck, (key, norm) in records.items():
            ours, fable = self.curated(key, norm), self.fable_answer(norm)
            if ours is None or fable is None:
                results[ck] = None
                continue
            result = results[ck] = {"norm": norm, "fable_model": fable["model"], "fable_effort": fable["effort"],
                                    "fable_certain": fable["certain"], "parts": {}}
            for p in PARTS:
                status = self.compare(p, ours[p], fable[p])
                if status == "disagrees" and ours[p].tier == "rule":
                    status = "rule disputed"
                elif status == "disagrees":
                    disputed.setdefault(ck, {})[p] = (ours[p], fable[p])
                result["parts"][p] = {"ours": ours[p], "fable": fable[p], "status": status, "final": ours[p]}
        decisions = self.judge({ck: (results[ck]["norm"], parts) for ck, parts in disputed.items()}, offline)
        for ck, parts in disputed.items():
            for p, (ours, theirs) in parts.items():
                part, decision = results[ck]["parts"][p], decisions.get(ck, {}).get(p)
                if decision is None:
                    part["status"] = "disagrees (not judged)"
                elif decision.key == ours.key:
                    part["status"] = "judge kept curation"
                else:
                    part["final"] = decision
                    part["status"] = ("judge: unknown" if decision.key == "unknown" else
                                      "judge took Fable's" if decision.key == theirs[0] else "judge chose another")
        return results

    def judge(self, items, offline=False):
        """items: {main cache key: (norm, {part: (our Answer, Fable's [key, field, quote])})}. Returns
        {main cache key: {part: Answer}} for the records decided, from the cache or (unless offline) by asking.
        The answers' order comes from Fable's judge key whatever the judge model, so a second judge sees each
        record exactly as Fable's judge did; only its cache key carries its own model."""
        decisions, todo = {}, {}
        for ck, (norm, parts) in items.items():
            shown = {p: [ours.key, ours.field, ours.quote, *theirs] for p, (ours, theirs) in parts.items()}
            order = self._key("judge", norm, [JUDGE_EFFORT, shown])
            jkey = order if self.judge_model == FABLE else self._key("judge", norm, [JUDGE_EFFORT, shown],
                                                                     self.judge_model)
            cached = self.cache.get(jkey)
            if cached:
                decisions[ck] = {p: Answer(*v) for p, v in cached.items()}
            elif not offline:
                todo[ck] = (norm, parts, jkey, order)
        if not todo:
            return decisions
        self.log(f"Judge ({self.judge_model}, {JUDGE_EFFORT} effort): {len(todo):,} record(s) with a disagreement")
        self.stats["judged"] += len(todo)
        ids = {ck: i for i, ck in enumerate(todo)}
        disputes = {}
        for ck, (norm, parts, jkey, order) in todo.items():
            shown = disputes[ids[ck]] = {}
            for p, (ours, theirs) in parts.items():
                a = {p: ours.key, f"{p}_field": ours.field, f"{p}_quote": ours.quote}
                b = {p: theirs[0], f"{p}_field": theirs[1], f"{p}_quote": theirs[2]}
                shown[p] = (a, b) if int(order[:8], 16) % 2 else (b, a)  # neither side is always reviewer A
        answers, = self.labeler._call_many([("judge", {ids[ck]: todo[ck][0] for ck in todo}, disputes)])
        tier = "fable judge" if self.judge_model == FABLE else "second judge"
        for ck, (norm, parts, jkey, order) in todo.items():
            ans = answers.get(ids[ck])
            if ans is None:  # the request failed: judge again next time
                self.stats["judge gave no answer"] += 1
                continue
            model = ans.get("_model") or self.judge_model
            decisions[ck] = {p: (Answer(ans[p], ans[f"{p}_field"], ans[f"{p}_quote"], tier, model)
                                 if self.labeler._valid(norm, ans, p) and ans[p] != "unknown"
                                 else Answer("unknown", tier="withheld by judge", model=model)) for p in parts}
            self.cache.put(jkey, {p: a.as_list() for p, a in decisions[ck].items()})
        self.cache.save()
        return decisions


def summarize(results):
    """Report counts: per part, the comparison outcome overall and by the curation's tier."""
    out = {"records": len(results), "awaiting the curation's or Fable's answer": sum(r is None for r in results.values()),
           "fable not certain": sum(1 for r in results.values() if r and not r["fable_certain"]), "parts": {}}
    for p in PARTS:
        by_tier = {}
        for r in results.values():
            if r is not None:
                part = r["parts"][p]
                by_tier.setdefault(part["ours"].tier, Counter())[part["status"]] += 1
        total = sum(by_tier.values(), Counter())
        out["parts"][p] = {"all": dict(total.most_common()),
                           "by curation tier": {t: dict(c.most_common()) for t, c in sorted(by_tier.items())}}
    return out


def _labeller(vocab):
    def label(key, part):
        return (vocab.hosts if part == "host" else vocab.terms)[key]["label"] if key else ""
    return label


def write_disagreements(results, path, vocab):
    """One row per part where Fable did not agree: the record, both answers, the outcome and the final answer."""
    label = _labeller(vocab)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["record", "part", "status", "curated", "curated_tier", "curated_quote", "fable",
                         "fable_quote", "final", "final_tier", "final_quote"])
        for r in results.values():
            if r is None:
                continue
            for p in PARTS:
                part = r["parts"][p]
                if part["status"] in AGREEING:
                    continue
                ours, fable, final = part["ours"], part["fable"] or [None, "", ""], part["final"]
                record = "; ".join(f"{k}={v}" for k, v in r["norm"].items())
                writer.writerow([record, p, part["status"], label(ours.key, p), ours.tier, ours.quote,
                                 label(fable[0], p), fable[2], label(final.key, p), final.tier, final.quote])


# --- a second judge ------------------------------------------------------------------------------------
def judged_items(results):
    """The records whose model-decided disagreements the judge settled, in the form judge() takes:
    {main cache key: (norm, {part: (our Answer, Fable's answer)})}, with the same parts check() sent."""
    items = {}
    for ck, r in results.items():
        if r is None:
            continue
        parts = {p: (part["ours"], part["fable"]) for p, part in r["parts"].items() if part["status"] in JUDGED}
        if parts:
            items[ck] = (r["norm"], parts)
    return items


def third_party_sample(items, model, n, seed=1):
    """A random sample of n judged records with no curated answer from model, so that model judges them as a
    third party. Returns (sample, how many records qualified)."""
    eligible = sorted(ck for ck, (_, parts) in items.items() if all(ours.model != model for ours, _ in parts.values()))
    chosen = random.Random(seed).sample(eligible, min(n, len(eligible)))
    return {ck: items[ck] for ck in chosen}, len(eligible)


def side(decision, ours, theirs):
    """Which answer a judge's decision took: the curation's, Fable's, neither (unknown) or another term."""
    if decision.key == ours.key:
        return "curation"
    if decision.key == theirs[0]:
        return "Fable's"
    return "unknown" if decision.key == "unknown" else "another term"


def compare_judges(results, items, second):
    """One row per disputed part both judges decided. results: check()'s, which holds the first (Fable) judge's
    decisions; items: the judged records asked again; second: the second judge's decisions, from judge()."""
    rows = []
    for ck, (norm, parts) in items.items():
        if ck not in second:
            continue
        for p, (ours, theirs) in parts.items():
            first, other = results[ck]["parts"][p]["final"], second[ck][p]
            rows.append({"norm": norm, "part": p, "ours": ours, "fable": theirs,
                         "fable_certain": results[ck]["fable_certain"], "first": first, "second": other,
                         "first_side": side(first, ours, theirs), "second_side": side(other, ours, theirs)})
    return rows


def mcnemar(b, c):
    """Exact two-sided McNemar p-value for b and c discordant pairs (binomial test of b of b + c at 1/2)."""
    n = b + c
    if n == 0:
        return 1.0
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(min(b, c) + 1)) / 2 ** n)


def summarize_judges(rows):
    """How often the second judge decided as the first: overall, per part, per curation tier and by whether
    Fable was certain. A first judge that favours Fable's answers would take them, where the second judge keeps
    the curation's, far more often than the reverse; the McNemar test compares those two counts."""
    def block(rs):
        same = sum(r["first"].key == r["second"].key for r in rs)
        sides = Counter((r["first_side"], r["second_side"]) for r in rs)
        b, c = sides[("Fable's", "curation")], sides[("curation", "Fable's")]
        return {"parts": len(rs), "same decision": same, "same decision share": round(same / len(rs), 3),
                "took Fable's, first judge": sum(r["first_side"] == "Fable's" for r in rs),
                "took Fable's, second judge": sum(r["second_side"] == "Fable's" for r in rs),
                "first took Fable's, second kept curation": b, "first kept curation, second took Fable's": c,
                "McNemar exact p": round(mcnemar(b, c), 4),
                "first judge / second judge": {f"{x} / {y}": n for (x, y), n in sides.most_common()}}

    tiers_seen = sorted({r["ours"].tier for r in rows})
    groups = {"all": rows, **{p: [r for r in rows if r["part"] == p] for p in PARTS},
              **{f"curation tier {t}": [r for r in rows if r["ours"].tier == t] for t in tiers_seen},
              "Fable certain": [r for r in rows if r["fable_certain"]],
              "Fable not certain": [r for r in rows if not r["fable_certain"]]}
    return {name: block(rs) for name, rs in groups.items() if rs}


def write_judge_comparison(rows, path, vocab):
    """One row per part: the record, the curation's and Fable's answers, and both judges' decisions."""
    label = _labeller(vocab)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["record", "part", "curated", "curated_tier", "curated_quote", "fable", "fable_quote",
                         "fable_certain", "first_judge", "first_side", "second_judge", "second_side",
                         "second_quote", "second_model"])
        for r in rows:
            p, ours, fable, second = r["part"], r["ours"], r["fable"], r["second"]
            writer.writerow(["; ".join(f"{k}={v}" for k, v in r["norm"].items()), p, label(ours.key, p), ours.tier,
                             ours.quote, label(fable[0], p), fable[2], r["fable_certain"], label(r["first"].key, p),
                             r["first_side"], label(second.key, p), r["second_side"], second.quote, second.model])
