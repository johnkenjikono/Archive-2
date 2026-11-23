"""Run from the repo root: python -m unittest catalogue.test_catalogue   (offline: NCBI and the API are faked)"""
import argparse
import csv
import gzip
import json
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path
from types import SimpleNamespace

from catalogue import build
from catalogue.batch import BatchTransport, BudgetExceeded, Ledger, estimate_usd
from catalogue.crosscheck import FABLE, CrossCheck, open_records, summarize, write_disagreements
from catalogue.hosts import HostResolver, candidates
from catalogue.ncbi import compact, harvest, iter_assemblies, mag_accessions, resolve_taxa
from steps.metadata_curation import CuratedLabeler, Vocabulary, usage_cost
from steps.test_metadata_curation import HAIKU, OPUS, SONNET, FakeClient, key, label


def quiet(*_):
    pass


class FakeBatches:
    """messages.batches stand-in: answers every request through a FakeClient when the batch is created.
    force(custom_id, params) may return "errored" / "expired" (no message) or a stop_reason to stamp on it."""

    def __init__(self, fake, force=None):
        self.fake, self.force = fake, force or (lambda custom_id, params: None)
        self.store, self.created = {}, 0

    def create(self, requests):
        self.created += 1
        batch_id = f"msgbatch_{self.created}"
        results = []
        for request in requests:
            forced = self.force(request["custom_id"], request["params"])
            if forced in ("errored", "expired", "canceled"):
                results.append(SimpleNamespace(custom_id=request["custom_id"],
                                               result=SimpleNamespace(type=forced, message=None)))
                continue
            message = self.fake._create(**request["params"])
            message.model = request["params"]["model"]
            if forced:
                message.stop_reason = forced
            results.append(SimpleNamespace(custom_id=request["custom_id"],
                                           result=SimpleNamespace(type="succeeded", message=message)))
        self.store[batch_id] = results
        return SimpleNamespace(id=batch_id, processing_status="in_progress")

    def retrieve(self, batch_id):
        types = [r.result.type for r in self.store[batch_id]]
        counts = SimpleNamespace(processing=0, **{t: types.count(t) for t in ("succeeded", "errored", "expired",
                                                                               "canceled")})
        return SimpleNamespace(id=batch_id, processing_status="ended", request_counts=counts)

    def results(self, batch_id):
        return iter(self.store[batch_id])


def batch_client(fake, force=None):
    return SimpleNamespace(messages=SimpleNamespace(create=fake._create, batches=FakeBatches(fake, force)),
                           beta=fake.beta)


def good(fields, rid, certain=True):
    if "chicken" in fields["isolation_source"]:
        return label(rid, "no_host", "isolation_source", "raw chicken meat", "chicken_meat", "isolation_source",
                     "raw chicken meat", certain)
    return label(rid, "human", "isolation_source", "human", "host_sample", "isolation_source", "human listeriosis",
                 certain)


class BatchTransportTest(unittest.TestCase):
    """The tiers through the Batch API: same answers as the synchronous path, plus the safety rails."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)
        self.keys = [key(isolation_source="pre-made raw chicken meat"), key(isolation_source="human listeriosis")]

    def curate(self, answer, client, budget=10.0, cache="cache.json"):
        transport = BatchTransport(self.dir, budget, client=client, sleep=quiet, log=quiet)
        labeler = CuratedLabeler(client=client, cache_path=self.dir / cache, transport=transport, log=quiet)
        return labeler.curate(self.keys), transport, labeler

    @staticmethod
    def disputed(model, rid, fields, call):
        a = good(fields, rid, certain=model != HAIKU)
        if model == SONNET and call == 2 and rid == 0:
            a.update(source="poultry_meat", source_quote="chicken meat")
        return a

    def test_same_answers_as_the_synchronous_path(self):
        fake = FakeClient(self.disputed)
        client = batch_client(fake)
        details, transport, labeler = self.curate(self.disputed, client)
        sync = CuratedLabeler(client=FakeClient(self.disputed), cache_path=self.dir / "sync.json", log=quiet)
        self.assertEqual(details, sync.curate(self.keys))
        self.assertEqual(details[self.keys[0]]["source_tier"], "opus")
        self.assertEqual(fake.models_called(), [HAIKU, SONNET, SONNET, OPUS])
        self.assertEqual(client.messages.batches.created, 3)  # Haiku, both Sonnet runs together, Opus
        params = client.messages.batches.store["msgbatch_1"][0].result.message
        self.assertNotIn("fallbacks", fake.calls[0])  # not allowed in batch requests
        self.assertEqual(fake.calls[0]["system"][0]["cache_control"], {"type": "ephemeral", "ttl": "1h"})
        self.assertIsNotNone(params)
        ledger = Ledger(self.dir / "spend.jsonl")
        expected = sum(usage_cost(m, SimpleNamespace(input_tokens=100, output_tokens=50), batch=True)
                       for m in (HAIKU, SONNET, SONNET, OPUS))
        self.assertAlmostEqual(ledger.total, expected)
        self.assertAlmostEqual(labeler.cost(), expected)
        self.assertEqual(len(json.loads((self.dir / "batches.json").read_text())), 3)

    def test_an_interrupted_run_resumes_its_batches_instead_of_paying_again(self):
        fake = FakeClient(self.disputed)
        client = batch_client(fake)
        first, _, _ = self.curate(self.disputed, client)
        total = Ledger(self.dir / "spend.jsonl").total
        again, _, _ = self.curate(self.disputed, client, cache="fresh_cache.json")  # same open records
        self.assertEqual(again, first)
        self.assertEqual(client.messages.batches.created, 3)
        self.assertAlmostEqual(Ledger(self.dir / "spend.jsonl").total, total)  # each batch charged once

    def test_budget_stops_before_anything_is_sent(self):
        client = batch_client(FakeClient(self.disputed))
        with self.assertRaises(BudgetExceeded):
            self.curate(self.disputed, client, budget=1e-6)
        self.assertEqual(client.messages.batches.created, 0)

    def test_refused_sonnet_chunk_is_retried_synchronously_with_fallbacks(self):
        answer = lambda m, rid, f, c: good(f, rid, certain=m != HAIKU)  # noqa: E731
        fake = FakeClient(answer)
        client = batch_client(fake, lambda cid, params: "refusal" if params["model"] == SONNET and cid == "j0c0" else None)
        details, _, _ = self.curate(answer, client)
        self.assertEqual(details[self.keys[0]]["source_tier"], "sonnet x2")
        retried = [c for c in fake.calls if c.get("fallbacks") == "default"]
        self.assertEqual([c["model"] for c in retried], [SONNET])

    def test_refused_haiku_chunk_escalates_without_a_retry(self):
        answer = lambda m, rid, f, c: good(f, rid)  # noqa: E731
        fake = FakeClient(answer)
        client = batch_client(fake, lambda cid, params: "refusal" if params["model"] == HAIKU else None)
        details, _, _ = self.curate(answer, client)
        self.assertEqual(fake.models_called(), [HAIKU, SONNET, SONNET])
        self.assertEqual(details[self.keys[1]]["host_tier"], "sonnet x2")

    def test_errored_requests_are_retried(self):
        answer = lambda m, rid, f, c: good(f, rid)  # noqa: E731
        fake = FakeClient(answer)
        client = batch_client(fake, lambda cid, params: "errored" if params["model"] == HAIKU else None)
        details, _, _ = self.curate(answer, client)
        self.assertEqual(details[self.keys[0]]["source_tier"], "haiku")
        self.assertEqual(fake.models_called(), [HAIKU])  # only the synchronous retry produced a message
        self.assertEqual(Ledger(self.dir / "spend.jsonl").entries[-1]["kind"], "sync")

    def test_dropped_connections_while_polling_are_waited_out(self):
        answer = lambda m, rid, f, c: good(f, rid)  # noqa: E731
        client = batch_client(FakeClient(answer))
        batches, failures = client.messages.batches, {"left": 2}
        real_retrieve = batches.retrieve

        def flaky_retrieve(batch_id):
            if failures["left"]:
                failures["left"] -= 1
                raise ConnectionError("dropped")
            return real_retrieve(batch_id)

        batches.retrieve = flaky_retrieve
        with unittest.mock.patch("catalogue.batch._transient_errors", lambda: (ConnectionError,)):
            details, _, _ = self.curate(answer, client)
        self.assertEqual(details[self.keys[0]]["source_tier"], "haiku")
        self.assertEqual(failures["left"], 0)

    def test_estimate_is_an_upper_bound(self):
        fake = FakeClient(lambda m, rid, f, c: good(f, rid))
        labeler = CuratedLabeler(client=fake, log=quiet)
        request = {"params": labeler.build_request("haiku", {0: {"isolation_source": "blood"}}, batch=True)}
        self.assertGreater(estimate_usd([request]), usage_cost(HAIKU, SimpleNamespace(input_tokens=2000,
                                                                                      output_tokens=100), batch=True))


class RecordingFake(FakeClient):
    """A FakeClient whose answer function can also look at the request being answered (self.request)."""

    def _create(self, **kwargs):
        self.request = kwargs
        return super()._create(**kwargs)


def curated_answer(model, rid, fields, call):
    """The curation's (Haiku's) certain answers for the four cross-check records."""
    src = fields.get("isolation_source", "")
    if "chicken" in src:
        return label(rid, "no_host", "isolation_source", "raw chicken meat", "chicken_meat", "isolation_source",
                     "raw chicken meat")
    if "cheese" in src:
        return label(rid, "no_host", "isolation_source", "cheese brine", "cheese", "isolation_source", "cheese")
    if "blood" in src:
        return label(rid, "human", "isolation_source", "patient", "host_sample", "isolation_source", "patient")
    return label(rid, "human", "host", "homo sapiens")  # "ward 7 bench": source unknown


class CrossCheckTest(unittest.TestCase):
    """Four records the curation sent to the models, checked by a fake Fable:
    chicken meat (agrees), cheese brine (another dairy term: same niche), blood of a patient (Fable says blood,
    the curation 'host sample': judged), and a human host by rule that Fable calls cattle (rule disputed)."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)
        self.vocab = Vocabulary()
        niche = self.vocab.terms["cheese"]["niche"]
        self.other_dairy = next(k for k, t in self.vocab.terms.items() if t["niche"] == niche and k != "cheese")
        self.keys = [key(isolation_source="pre-made raw chicken meat"), key(isolation_source="cheese brine from a dairy"),
                     key(isolation_source="blood of a patient on the ward"),
                     key(host="Homo sapiens", isolation_source="ward 7 bench")]
        self.main = CuratedLabeler(vocab=self.vocab, client=FakeClient(curated_answer), cache_path=self.dir / "main.json",
                                   log=quiet)
        self.main.curate(self.keys)
        self.records, self.key_to_ck = open_records(self.main, self.keys)

    def fable(self, effort="low", transport=None, client=None):
        fake = client or RecordingFake(None)
        if client is None:
            fake.answer = lambda model, rid, fields, call: self.fable_answer(fake.request, rid, fields)
        xc = CrossCheck(self.main, self.dir / "crosscheck.json", effort=effort, client=fake, transport=transport,
                        log=quiet)
        return xc, fake

    def fable_answer(self, request, rid, fields):
        src = fields.get("isolation_source", "")
        judging = "disputed" in request["messages"][0]["content"]
        if "chicken" in src:
            return curated_answer(None, rid, fields, 0)
        if "cheese" in src:
            return label(rid, "no_host", "isolation_source", "cheese brine", self.other_dairy, "isolation_source",
                         "cheese brine")
        if "blood" in src:
            return label(rid, "human", "isolation_source", "patient", "blood", "isolation_source",
                         "blood" if not judging else "blood of a patient")
        return label(rid, "cattle", "host", "homo sapiens")

    def statuses(self, results):
        by_source = {self.records[ck][1].get("isolation_source"): r for ck, r in results.items()}
        return {src: (r["parts"]["host"]["status"], r["parts"]["source"]["status"]) for src, r in by_source.items()}

    def test_outcomes_and_what_is_judged(self):
        xc, fake = self.fable()
        results = xc.check(self.records)
        self.assertEqual(self.statuses(results), {
            "pre-made raw chicken meat": ("agrees", "agrees"),
            "cheese brine from a dairy": ("agrees", "same niche"),
            "blood of a patient on the ward": ("agrees", "judge took Fable's"),
            "ward 7 bench": ("rule disputed", "agrees")})
        label_call, judge_call = fake.calls
        self.assertEqual((label_call["model"], label_call["output_config"]["effort"]), (FABLE, "low"))
        self.assertEqual((judge_call["model"], judge_call["output_config"]["effort"]), (FABLE, "high"))
        prompt = judge_call["messages"][0]["content"]
        self.assertIn("disputed source", prompt)  # only the model-decided disagreement is judged
        self.assertNotIn("disputed host", prompt)
        self.assertIn("blood of a patient on the ward", prompt)
        self.assertNotIn("ward 7 bench", prompt)
        self.assertEqual((judge_call["betas"], judge_call["fallbacks"]), (["server-side-fallback-2026-07-01"], "default"))
        blood = next(r for r in results.values() if "blood" in r["norm"]["isolation_source"])
        final = blood["parts"]["source"]["final"]
        self.assertEqual((final.key, final.tier, final.quote), ("blood", "fable judge", "blood of a patient"))
        counts = summarize(results)["parts"]["source"]
        self.assertEqual(counts["by curation tier"]["haiku"]["judge took Fable's"], 1)
        tsv = self.dir / "disagreements.tsv"
        write_disagreements(results, tsv, self.vocab)
        with open(tsv, encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(sorted((r["part"], r["status"]) for r in rows),
                         [("host", "rule disputed"), ("source", "judge took Fable's")])

    def test_answers_and_decisions_are_cached_per_effort(self):
        self.fable()[0].check(self.records)
        again, fake = self.fable(client=FakeClient(lambda *a: self.fail("cached answers must not be asked again")))
        self.assertEqual(self.statuses(again.check(self.records)), self.statuses(again.check(self.records, offline=True)))
        self.assertEqual(fake.calls, [])
        high, fake = self.fable("high")
        high.check(self.records)
        self.assertEqual([c["output_config"]["effort"] for c in fake.calls], ["high"])  # labels only: same dispute
        self.assertEqual(high.effort_agreement(self.records, "low")["source"], {"same": 4})

    def test_not_curated_records_wait(self):
        unseen = {**self.records, **open_records(self.main, [key(isolation_source="mystery pond scum")])[0]}
        results = self.fable()[0].check(unseen)
        self.assertEqual(sum(r is None for r in results.values()), 1)

    def test_through_the_batch_api(self):
        fake = RecordingFake(None)
        fake.answer = lambda model, rid, fields, call: self.fable_answer(fake.request, rid, fields)
        transport = BatchTransport(self.dir, 10.0, client=batch_client(fake), sleep=quiet, log=quiet,
                                   state_name="crosscheck_batches.json")
        xc, _ = self.fable(transport=transport, client=fake)
        self.assertEqual(self.statuses(xc.check(self.records))["blood of a patient on the ward"],
                         ("agrees", "judge took Fable's"))
        self.assertTrue((self.dir / "crosscheck_batches.json").exists())
        self.assertFalse((self.dir / "batches.json").exists())
        entries = Ledger(self.dir / "spend.jsonl").entries
        self.assertEqual({e["model"] for e in entries}, {FABLE})
        self.assertGreater(xc.cost(), 0)


class TransportStateTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)

    def test_saving_keeps_stages_another_process_saved(self):
        first = BatchTransport(self.dir, 10.0, client=None, log=quiet)
        second = BatchTransport(self.dir, 10.0, client=None, log=quiet)
        first.state["a"] = {"name": "first"}
        first._save_state()
        second.state["b"] = {"name": "second"}
        second._save_state()
        self.assertEqual(set(json.loads((self.dir / "batches.json").read_text())), {"a", "b"})

    def test_budget_counts_what_other_commands_spent(self):
        transport = BatchTransport(self.dir, 10.0, client=None, log=quiet)
        transport._check_budget(1.0, "fits")
        Ledger(self.dir / "spend.jsonl").add(stage="another command", kind="batch", model="m", usd=9.5)
        with self.assertRaises(BudgetExceeded):
            transport._check_budget(1.0, "no longer fits")


class ModelAttributionTest(unittest.TestCase):
    def test_the_answering_model_is_recorded(self):
        labeler = CuratedLabeler(log=quiet)
        response = SimpleNamespace(stop_reason="end_turn", model="claude-opus-4-8", content=[SimpleNamespace(
            type="text", text=json.dumps({"labels": [label(0, "human", "host", "homo sapiens")]}))])
        answers, status = labeler.parse_response("sonnet", {0: {}}, response)
        self.assertEqual((status, answers[0]["_model"]), ("ok", "claude-opus-4-8"))


class RulesRevisionTest(unittest.TestCase):
    def test_rule_changes_keep_paid_answers_valid(self):
        """Synonyms and skip lists (rules_revision) are not part of the answer cache key; the version is."""
        labeler = CuratedLabeler(log=quiet)
        fields = {"isolation_source": "pre-made raw chicken meat"}
        before = labeler._cache_key(fields)
        labeler.vocab.rules_revision += 1
        self.assertEqual(labeler._cache_key(fields), before)
        labeler.vocab.version += ".x"
        self.assertNotEqual(labeler._cache_key(fields), before)


class OfflineTest(unittest.TestCase):
    def test_offline_curation_neither_calls_models_nor_caches(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        cache = Path(tmp.name) / "cache.json"
        client = FakeClient(lambda *a: self.fail("offline must not call a model"))
        labeler = CuratedLabeler(client=client, cache_path=cache, log=quiet)
        keys = [key(isolation_source="pre-made raw chicken meat"), key(isolation_source="blood", host="Homo sapiens")]
        details = labeler.curate(keys, offline=True)
        self.assertEqual(details[keys[0]]["source_tier"], "not curated")
        self.assertEqual((details[keys[1]]["source_tier"], details[keys[1]]["host"]), ("rule", "human"))
        self.assertFalse(cache.exists() and json.loads(cache.read_text())["entries"])


class HostOnlyTest(unittest.TestCase):
    """A record with a host and no source field: host-associated, site not stated (catalogue option only)."""

    def test_source_from_host(self):
        keys = [key(host="Homo sapiens"), key(host="Oryza sativa"), key(host="Homo sapiens", isolation_source="blood")]
        details = CuratedLabeler(log=quiet, source_from_host=True).curate(keys, offline=True)
        human, rice, blood = (details[k] for k in keys)
        self.assertEqual((human["source_term"], human["source_category"], human["niche"]),
                         ("host sample, site not stated", "human-associated", "host sample, site not stated"))
        self.assertEqual((human["source_tier"], human["source_evidence_field"], human["source_quote"]),
                         ("rule", "host", "Homo sapiens"))
        self.assertEqual((rice["source_term"], rice["source_category"]), ("plant, part not stated", "plant-associated"))
        self.assertEqual(blood["source_term"], "blood")  # a stated source is never replaced
        # A culture descriptor is not a source (the pilot's Haiku read this one as "laboratory strain").
        pure = key(host="Homo sapiens", sample_type="Monoclonal strain")
        labeled = CuratedLabeler(log=quiet, source_from_host=True).curate([pure], offline=True)[pure]
        self.assertEqual((labeled["source_term"], labeled["source_tier"]), ("host sample, site not stated", "rule"))
        pipeline = CuratedLabeler(log=quiet).curate(keys[:1], offline=True)  # off by default
        self.assertEqual(pipeline[keys[0]]["source_term"], "unknown")


TAXA = {
    "ixodes scapularis": {"tax_id": 6945, "name": "Ixodes scapularis", "rank": "species",
                          "parents": [1, 131567, 2759, 33208, 6656, 6935]},
    "oryza sativa": {"tax_id": 4530, "name": "Oryza sativa", "rank": "species", "parents": [1, 131567, 2759, 33090]},
    "cancer": {"tax_id": 6754, "name": "Cancer", "rank": "genus", "parents": [1, 131567, 2759, 33208]},
    "escherichia coli": {"tax_id": 562, "name": "Escherichia coli", "rank": "species", "parents": [1, 131567, 2]},
    "honeybee": {"tax_id": 7460, "name": "Apis mellifera", "rank": "species", "parents": [1, 2759, 33208, 50557]},
    "9825": {"tax_id": 9825, "name": "Sus scrofa domesticus", "rank": "subspecies",
             "parents": [1, 2759, 33208, 40674, 9823]},
    "homo sapiens": {"tax_id": 9606, "name": "Homo sapiens", "rank": "species", "parents": [1, 2759, 33208, 40674]},
    "rat": {"ambiguous": [{"tax_id": 10116}, {"tax_id": 10114}]},
}


class HostResolverTest(unittest.TestCase):
    def setUp(self):
        self.resolver = HostResolver(Vocabulary(), TAXA)

    def test_candidates(self):
        self.assertIn("oryza sativa", candidates("Oryza sativa L."))
        self.assertEqual(candidates("pig; slaughterhouse")[:2], ["pig; slaughterhouse", "pig"])
        self.assertIn("cattle", candidates("Bos taurus (cattle)"))
        self.assertEqual(candidates("missing")[0], "missing")  # informative() filters these before lookup

    def test_lineage_places_any_taxon_in_the_vocabulary(self):
        self.assertEqual(self.resolver.vocab_key_for_name("Ixodes scapularis"), "tick")
        self.assertEqual(self.resolver.vocab_key_for_name("Oryza sativa L."), "rice")
        self.assertEqual(self.resolver.vocab_key_for_name("honeybee"), "insect")
        self.assertEqual(self.resolver.vocab_key_for_taxid("9825"), "pig")

    def test_words_that_are_genus_names_bacteria_and_ambiguity_are_refused(self):
        self.assertIsNone(self.resolver.vocab_key_for_name("cancer"))
        self.assertIsNone(self.resolver.vocab_key_for_name("Escherichia coli"))
        self.assertIsNone(self.resolver.vocab_key_for_name("rat"))

    def test_exact_host(self):
        self.assertEqual(self.resolver.exact({"host_taxid": "9825", "host": "Homo sapiens"}),
                         (9825, "Sus scrofa domesticus"))
        self.assertEqual(self.resolver.exact({"host_common_name": "honeybee"}, ("host", "host_common_name")),
                         (7460, "Apis mellifera"))
        self.assertEqual(self.resolver.exact({"host": "cancer"}), ("", ""))

    def test_rules_use_the_resolver(self):
        labeler = CuratedLabeler(log=quiet, host_resolver=self.resolver)
        fields = {"isolation_source": "blood", "host": "Ixodes scapularis"}
        host = labeler.rule_host(fields, labeler.rule_source(fields))
        self.assertEqual((host.key, host.tier), ("tick", "rule"))
        self.assertIsNone(CuratedLabeler(log=quiet).rule_host(fields, None))  # without it the models decide


class CheckoutTest(unittest.TestCase):
    def test_steps_modules_must_come_from_this_checkout(self):
        self.assertEqual(build.foreign_modules(), [])
        stray = SimpleNamespace(__file__="/elsewhere/repo/steps/niche_metadata.py")
        with unittest.mock.patch.dict(sys.modules, {"steps.niche_metadata": stray}):
            self.assertEqual(build.foreign_modules(), ["/elsewhere/repo/steps/niche_metadata.py"])


class FieldTest(unittest.TestCase):
    def test_canonical_names_and_ena_aliases(self):
        attrs = build.canonical({"local environmental context": "marine biome", "Isolation Site": "Soil",
                                 "Source": "BLOOD", "metagenomic source": "soil metagenome", "host": "Homo sapiens"})
        self.assertEqual((attrs["env_local_scale"], attrs["isolation_site"], attrs["source"],
                          attrs["metagenome_source"]), ("marine biome", "Soil", "BLOOD", "soil metagenome"))
        attrs = build.canonical({"env_local_scale": "kept", "local environmental context": "alias"})
        self.assertEqual(attrs["env_local_scale"], "kept")

    def test_record_key_order_and_content(self):
        k = build.record_key({"habitat": "coral reef", "isolation_source": "seawater", "strain": "X1", "host": ""})
        self.assertEqual(k, (("isolation_source", "seawater"), ("habitat", "coral reef")))


def fake_report(accession, tax_id, attributes, level="Complete Genome"):
    return {"accession": accession, "paired_accession": accession.replace("GCA", "GCF"),
            "organism": {"tax_id": tax_id, "organism_name": f"Organism {tax_id}", "infraspecific_names": {"strain": "S1"}},
            "assembly_info": {"assembly_level": level, "release_date": "2020-01-01", "bioproject_accession": "PRJNA1",
                              "biosample": {"accession": f"SAMN{tax_id}",
                                            "attributes": [{"name": k, "value": v} for k, v in attributes.items()]}},
            "checkm_info": {"completeness": 99.5, "contamination": 0.2},
            "average_nucleotide_identity": {"taxonomy_check_status": "OK"},
            "assembly_stats": {"total_sequence_length": "3000000", "gc_percent": 38.0}}


class NcbiTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)

    def test_compact_keeps_what_the_catalogue_needs(self):
        c = compact(fake_report("GCA_1.1", 1639, {"isolation_source": "blood", "host": "Homo sapiens"}))
        self.assertEqual((c["accession"], c["paired_accession"], c["tax_id"], c["strain"], c["biosample"]),
                         ("GCA_1.1", "GCF_1.1", 1639, "S1", "SAMN1639"))
        self.assertEqual(c["attributes"], {"isolation_source": "blood", "host": "Homo sapiens"})
        self.assertEqual((c["checkm_completeness"], c["taxonomy_check"]), (99.5, "OK"))

    def test_harvest_pages_and_resumes(self):
        pages = {None: {"reports": [fake_report("GCA_1.1", 1, {})], "next_page_token": "t2", "total_count": 2},
                 "t2": {"reports": [fake_report("GCA_2.1", 2, {})], "total_count": 2}}
        calls = []

        def fetch(url, payload=None, fail_on=None):
            calls.append(url)
            if "returned_content=ASSM_ACC" in url:
                return {"reports": [{"accession": "GCA_2.1"}]}
            token = url.split("page_token=")[1] if "page_token=" in url else None
            if token == fail_on:
                raise OSError("network down")
            return pages[token]

        with self.assertRaises(OSError):
            harvest(self.dir, fetch=lambda u, p=None: fetch(u, p, fail_on="t2"), log=quiet)
        state = harvest(self.dir, fetch=fetch, log=quiet)
        self.assertTrue(state["done"])
        self.assertEqual([a["accession"] for a in iter_assemblies(self.dir)], ["GCA_1.1", "GCA_2.1"])
        self.assertEqual(mag_accessions(self.dir), {"GCA_2.1"})
        self.assertEqual(sum("page_token" not in u and "ASSM_ACC" not in u for u in calls), 1)  # page 1 fetched once
        self.assertEqual(harvest(self.dir, fetch=lambda *a: self.fail("already done"), log=quiet)["records"], 2)

    def test_taxonomy_lookups_are_batched_cached_and_flag_ambiguity(self):
        posts = []

        def fetch(url, payload=None):
            posts.append(payload["taxons"])
            reports = []
            for q in payload["taxons"]:
                if q == "homo sapiens":
                    reports.append({"query": [q], "taxonomy": {"tax_id": 9606, "rank": "SPECIES", "parents": [1, 2759],
                                                               "current_scientific_name": {"name": "Homo sapiens"},
                                                               "classification": {"species": {"id": 9606,
                                                                                              "name": "Homo sapiens"}}}})
                elif q == "rat":
                    for t in (10116, 10114):
                        reports.append({"query": [q], "taxonomy": {"tax_id": t, "parents": [1],
                                                                   "current_scientific_name": {"name": str(t)}}})
                else:
                    reports.append({"query": [q], "errors": [{"reason": "not recognized"}]})
            return {"reports": reports}

        cache = self.dir / "taxa.json"
        out = resolve_taxa(["homo sapiens", "rat", "bogus"], cache, fetch=fetch, batch=2, log=quiet)
        self.assertEqual((out["homo sapiens"]["tax_id"], out["homo sapiens"]["rank"]), (9606, "species"))
        self.assertIn("ambiguous", out["rat"])
        self.assertIsNone(out["bogus"])
        self.assertEqual(len(posts), 2)
        resolve_taxa(["homo sapiens", "bogus"], cache, fetch=lambda *a, **k: self.fail("cached"), log=quiet)


class AssembleTest(unittest.TestCase):
    """The assemble command on a three-genome work directory."""

    def test_catalogue_and_species_summary(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        paths = build.Paths(tmp.name)
        (paths.ncbi / "pages").mkdir(parents=True)
        reports = [fake_report("GCA_1.1", 1639, {"isolation_source": "blood", "host": "Homo sapiens",
                                                 "geo_loc_name": "USA: CA", "collection_date": "2019-03"}),
                   fake_report("GCA_2.1", 1639, {"isolation_source": "pre-made raw chicken meat"}),
                   fake_report("GCA_3.1", 562, {"metagenomic source": "soil metagenome"}, level="Chromosome")]
        with gzip.open(paths.ncbi / "pages" / "page_00001.jsonl.gz", "wt") as f:
            for r in reports:
                f.write(json.dumps(compact(r)) + "\n")
        (paths.ncbi / "harvest_state.json").write_text(json.dumps({"done": True, "records": 3, "query": {}}))
        (paths.ncbi / "mag_accessions.txt").write_text("GCA_3.1\n")
        lineage = lambda sid, name: {"tax_id": sid, "name": name, "rank": "species", "parents": [1, 2],  # noqa: E731
                                     "classification": {"species": {"id": sid, "name": name},
                                                        "genus": {"id": 0, "name": name.split()[0]},
                                                        "phylum": {"id": 0, "name": "Bacillota"}}}
        paths.taxa.write_text(json.dumps({"1639": lineage(1639, "Listeria monocytogenes"),
                                          "562": lineage(562, "Escherichia coli"), **TAXA}))
        manifest = build.cmd_assemble(argparse.Namespace(), paths, quiet)
        self.assertEqual((manifest["genomes"], manifest["species"]), (3, 2))
        with gzip.open(paths.release / "isolate_catalogue.tsv.gz", "rt") as f:
            rows = {r["accession"]: r for r in csv.DictReader(f, delimiter="\t")}
        blood, chicken, mag = rows["GCA_1.1"], rows["GCA_2.1"], rows["GCA_3.1"]
        self.assertEqual((blood["host"], blood["source_category"], blood["source_tier"], blood["host_ncbi_taxid"]),
                         ("human", "human-associated", "rule", "9606"))
        self.assertEqual((blood["species"], blood["country_insdc"], blood["collection_year"]),
                         ("Listeria monocytogenes", "USA", "2019"))
        self.assertEqual(chicken["source_tier"], "not curated")  # nothing in the cache yet, and assemble never calls
        self.assertEqual((mag["is_mag"], mag["raw_metagenome_source"]), ("True", "soil metagenome"))
        with open(paths.release / "species_summary.tsv") as f:
            summary = {r["species"]: r for r in csv.DictReader(f, delimiter="\t")}
        listeria = summary["Listeria monocytogenes"]
        self.assertEqual((listeria["genomes"], listeria["genomes_with_source"], listeria["top_source_categories"]),
                         ("2", "1", "human-associated 1"))
        self.assertEqual((blood["host_crosscheck"], chicken["source_crosscheck"]), ("rules only", "not checked"))

    def test_judge_corrections_reach_the_catalogue(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        paths = build.Paths(tmp.name)
        (paths.ncbi / "pages").mkdir(parents=True)
        reports = [fake_report("GCA_1.1", 1639, {"isolation_source": "blood of a patient on the ward"}),
                   fake_report("GCA_2.1", 1639, {"isolation_source": "pre-made raw chicken meat"}),
                   fake_report("GCA_3.1", 1639, {"isolation_source": "blood", "host": "Homo sapiens"})]
        with gzip.open(paths.ncbi / "pages" / "page_00001.jsonl.gz", "wt") as f:
            for r in reports:
                f.write(json.dumps(compact(r)) + "\n")
        (paths.ncbi / "harvest_state.json").write_text(json.dumps({"done": True, "records": 3, "query": {}}))
        vocab = Vocabulary()
        keys = sorted(build.distinct_keys(paths))
        main = build.make_labeler(paths, vocab, build.load_resolver(paths, vocab)[0], client=FakeClient(curated_answer),
                                  log=quiet)
        main.curate(keys)
        fake = RecordingFake(None)

        def fable(model, rid, fields, call):
            if "blood" in fields["isolation_source"]:
                judging = "disputed" in fake.request["messages"][0]["content"]
                return label(rid, "human", "isolation_source", "patient", "blood", "isolation_source",
                             "blood of a patient" if judging else "blood")
            return curated_answer(model, rid, fields, call)

        fake.answer = fable
        CrossCheck(main, paths.crosscheck, client=fake, log=quiet).check(open_records(main, keys)[0])
        manifest = build.cmd_assemble(argparse.Namespace(), paths, quiet)
        with gzip.open(paths.release / "isolate_catalogue.tsv.gz", "rt") as f:
            rows = {r["accession"]: r for r in csv.DictReader(f, delimiter="\t")}
        judged, chicken, ruled = rows["GCA_1.1"], rows["GCA_2.1"], rows["GCA_3.1"]
        self.assertEqual((judged["source_term"], judged["source_tier"], judged["source_quote"]),
                         ("blood", "fable judge", "blood of a patient"))
        self.assertEqual((judged["source_category"], judged["host"], judged["host_tier"]),
                         ("human-associated", "human", "haiku"))
        self.assertEqual((judged["source_crosscheck"], judged["fable_source_term"], judged["crosscheck_model"]),
                         ("judge took Fable's", "blood", f"{FABLE}, low effort"))
        self.assertEqual((chicken["source_term"], chicken["source_crosscheck"], chicken["host_crosscheck"]),
                         ("chicken meat", "agrees", "agrees"))
        self.assertEqual(ruled["source_crosscheck"], "rules only")
        self.assertEqual(manifest["models"]["crosscheck"], FABLE)


if __name__ == "__main__":
    unittest.main()
