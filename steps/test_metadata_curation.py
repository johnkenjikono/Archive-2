"""Run from the repo root: python -m unittest steps.test_metadata_curation
(set ONLINE_TESTS=1 to also check every vocabulary id against EBI OLS)"""
import csv
import json
import os
import re
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from steps.metadata_curation import (CURATION_COLUMNS, TIERS, Answer, CuratedLabeler, Vocabulary, quote_ok,
                                     run_curated_niche_analysis, verify_vocabulary)
from steps.niche_metadata import LLM_FIELDS
from steps.parsing import FULL_CORE_XML

HAIKU, SONNET, OPUS = (TIERS[t]["model"] for t in ("haiku", "sonnet", "opus"))


class FakeClient:
    """Stands in for anthropic.Anthropic(). answer(model, rid, fields, call) returns one label dict, or None
    to leave the record out; stop(model, n_records) can force a stop_reason. Records every request."""

    def __init__(self, answer, stop=None):
        self.answer, self.stop = answer, stop or (lambda model, n: "end_turn")
        self.calls = []
        self.messages = SimpleNamespace(create=self._create)
        self.beta = SimpleNamespace(messages=SimpleNamespace(create=self._create))

    def _create(self, **kwargs):
        self.calls.append(kwargs)
        model = kwargs["model"]
        records = {}
        for block in re.split(r"\n\n(?=Record \d+:)", kwargs["messages"][0]["content"]):
            match = re.search(r"Record (\d+):\n((?:  .*\n?)*)", block)
            if match:
                records[int(match.group(1))] = dict(
                    re.match(r"  (\S+) = (.*)", line).groups() for line in match.group(2).splitlines()
                    if re.match(r"  \S+ = ", line))
        call = sum(1 for c in self.calls if c["model"] == model)
        labels = [a for rid, fields in records.items() if (a := self.answer(model, rid, fields, call))]
        text = json.dumps({"labels": labels})
        return SimpleNamespace(stop_reason=self.stop(model, len(records)),
                               content=[SimpleNamespace(type="text", text=text)],
                               usage=SimpleNamespace(input_tokens=100, output_tokens=50))

    def models_called(self):
        return [c["model"] for c in self.calls]


def label(rid, host="unknown", host_field="", host_quote="", source="unknown", source_field="", source_quote="",
          certain=True):
    return {"id": rid, "host": host, "host_field": host_field, "host_quote": host_quote, "source": source,
            "source_field": source_field, "source_quote": source_quote, "certain": certain}


def key(**fields):
    return tuple((k, fields[k]) for k in LLM_FIELDS if fields.get(k))


class VocabularyTest(unittest.TestCase):
    def setUp(self):
        self.vocab = Vocabulary()

    def test_file_is_consistent(self):
        self.assertEqual(self.vocab.problems(), [])

    def test_every_id_is_well_formed_and_carries_its_ontology_label(self):
        for term in self.vocab.terms.values():
            if term["id"]:
                self.assertRegex(term["id"], r"^(ENVO|UBERON|FOODON|PO):\d+$")
                self.assertTrue(term["ontology_label"], term["key"])

    def test_unknown_or_relabelled_ids_are_reported(self):
        labels = {t["id"]: t["ontology_label"] for t in self.vocab.terms.values() if t["id"]}
        labels.update({f"NCBITaxon:{h['taxid']}": h["ncbi_name"] for h in self.vocab.hosts.values() if h["taxid"]})
        labels["ENVO:00001998"] = "sediment"  # the id now names another term
        del labels["UBERON:0000178"]  # the id is not in the ontology
        problems = verify_vocabulary(self.vocab, lookup=lambda c: [{"label": labels[c]}] if c in labels else [])
        self.assertEqual(len(problems), 2, problems)
        self.assertTrue(problems[0].startswith("UBERON:0000178: not found"), problems)
        self.assertIn("ENVO:00001998: OLS label 'sediment'", problems[1])

    @unittest.skipUnless(os.environ.get("ONLINE_TESTS"), "set ONLINE_TESTS=1 to query EBI OLS")
    def test_every_id_exists_in_ols_with_its_label(self):
        self.assertEqual(verify_vocabulary(self.vocab), [])


class RulesTest(unittest.TestCase):
    def setUp(self):
        self.labeler = CuratedLabeler(log=lambda *_: None)

    def rules(self, **fields):
        fields = {k: v for k, v in fields.items() if self.labeler.vocab.informative(v)}
        source = self.labeler.rule_source(fields)
        return source, self.labeler.rule_host(fields, source)

    def test_synonyms_and_hosts(self):
        source, host = self.rules(isolation_source="Blood", host="Homo sapiens")
        self.assertEqual((source.key, source.field, source.tier), ("blood", "isolation_source", "rule"))
        self.assertEqual((host.key, host.field), ("human", "host"))
        self.assertEqual(self.rules(isolation_source="CSF")[0].key, "cerebrospinal_fluid")

    def test_missing_and_growth_values_are_skipped(self):
        source, host = self.rules(isolation_source="missing", sample_type="whole organism", env_medium="soil")
        self.assertEqual((source.key, source.field), ("soil", "env_medium"))
        self.assertEqual(host.key, "no_host")  # soil has no living host

    def test_ontology_id_must_agree_with_its_label(self):
        self.assertEqual(self.rules(env_medium="rhizosphere [ENVO:00005801]")[0].key, "rhizosphere")
        self.assertEqual(self.rules(env_medium="ENVO:01001969")[0].key, "lab_culture")  # a bare id is taken as given
        self.assertIsNone(self.rules(env_medium="soil [ENVO:00002007]")[0])  # 00002007 is sediment
        self.assertIsNone(self.rules(env_medium="cropland [ENVO:00000096]")[0])  # an id outside the vocabulary

    def test_unknown_first_field_goes_to_the_models(self):
        self.assertIsNone(self.rules(isolation_source="pre-made raw chicken meat", env_medium="soil")[0])

    def test_host_names(self):
        self.assertEqual(self.rules(isolation_source="milk", host="Capra hircus coreanae")[1].key, "goat")
        self.assertEqual(self.rules(isolation_source="milk", host="caprinae (caprine)")[1].key, "goat")
        source, host = self.rules(isolation_source="human")
        self.assertEqual((source.key, host.key, host.field), ("host_sample", "human", "isolation_source"))
        self.assertIsNone(self.rules(isolation_source="blood", host="Microtus obscurus")[1])

    def test_a_sample_the_models_read_leaves_vague_or_indirect_hosts_to_them(self):
        # Fable's cross-check found these rule hosts beaten by the record's own text: the model that reads
        # "Rectum of pig" names the pig, where source_type "animal" or host "plant" name nothing in particular.
        for fields in ({"isolation_source": "Rectum of pig", "source_type": "animal"},
                       {"isolation_source": "sweet potato leaf", "sample_type": "environmental"},
                       {"isolation_source": "isolated from cucumber", "host": "plant"},
                       {"isolation_source": "citrus rhizosphere soil", "host": "soil"}):
            source, host = self.rules(**fields)
            self.assertEqual((source, host), (None, None), fields)
        # A named host still decides, and a sample the rules can read keeps its rule host.
        self.assertEqual(self.rules(isolation_source="Rectum of pig", host="Sus scrofa")[1].key, "pig")
        self.assertEqual(self.rules(isolation_source="feces", source_type="animal")[1].key, "other_animal")
        self.assertEqual(self.rules(isolation_source="soil", host="soil")[1].key, "no_host")

    def test_no_host_field(self):
        self.assertEqual(self.rules(isolation_source="cheese")[1].key, "no_host")  # food
        for source in ("blood", "laboratory strain", "rhizosphere"):  # body site, lab strain, plant: host not stated
            host = self.rules(isolation_source=source)[1]
            self.assertEqual((host.key, host.tier), ("unknown", "no evidence"), source)

    def test_body_site_category_follows_the_host(self):
        vocab = self.labeler.vocab
        self.assertEqual(vocab.category_of("feces", "human"), ("human_associated", "HumanGut"))
        self.assertEqual(vocab.category_of("feces", "cattle"), ("animal_associated", "HostAssociated"))
        self.assertEqual(vocab.category_of("feces", "unknown"), ("host_associated", "HostAssociated"))
        self.assertEqual(vocab.category_of("cheese", "unknown"), ("food", "FoodHumanFoods"))


class QuoteTest(unittest.TestCase):
    def test_quote_must_be_in_the_named_field(self):
        fields = {"isolation_source": "Pre-made  raw chicken meat", "host": "missing"}
        self.assertTrue(quote_ok(fields, "isolation_source", "raw CHICKEN meat", "chicken_meat"))
        self.assertFalse(quote_ok(fields, "host", "raw chicken meat", "chicken_meat"))
        self.assertFalse(quote_ok(fields, "isolation_source", "raw pork", "pork"))
        self.assertFalse(quote_ok(fields, "isolation_source", "", "chicken_meat"))
        self.assertTrue(quote_ok(fields, "", "", "unknown"))


class TierTest(unittest.TestCase):
    """Two records the rules cannot settle: 0 = 'pre-made raw chicken meat', 1 = 'human listeriosis'."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.cache = Path(tmp.name) / "cache.json"
        self.keys = [key(isolation_source="pre-made raw chicken meat"), key(isolation_source="human listeriosis")]

    def curate(self, answer, stop=None, cache=None):
        client = FakeClient(answer, stop)
        labeler = CuratedLabeler(client=client, cache_path=cache or self.cache, log=lambda *_: None)
        return labeler.curate(self.keys), client

    @staticmethod
    def good(fields, rid, certain=True):
        if "chicken" in fields["isolation_source"]:
            return label(rid, "no_host", "isolation_source", "raw chicken meat", "chicken_meat", "isolation_source",
                         "raw chicken meat", certain)
        return label(rid, "human", "isolation_source", "human", "host_sample", "isolation_source", "human listeriosis",
                     certain)

    def test_certain_haiku_with_good_quotes_is_accepted(self):
        details, client = self.curate(lambda m, rid, f, c: self.good(f, rid))
        self.assertEqual(client.models_called(), [HAIKU])
        chicken = details[self.keys[0]]
        self.assertEqual((chicken["source_term"], chicken["source_tier"], chicken["source_confidence"]),
                         ("chicken meat", "haiku", "medium"))
        self.assertEqual(chicken["source_category"], "food")
        self.assertEqual(details[self.keys[1]]["host"], "human")
        self.assertNotIn("betas", client.calls[0])  # no server-side fallback for Haiku

    def test_bad_quote_escalates_and_two_agreeing_sonnets_decide(self):
        def answer(model, rid, fields, call):
            a = self.good(fields, rid)
            if model == HAIKU:
                a["source_quote"] = "made up words"
            return a
        details, client = self.curate(answer)
        self.assertEqual(client.models_called(), [HAIKU, SONNET, SONNET])
        self.assertEqual(details[self.keys[0]]["source_tier"], "sonnet x2")
        self.assertEqual(client.calls[1]["fallbacks"], "default")

    def test_sonnet_disagreement_goes_to_opus(self):
        def answer(model, rid, fields, call):
            a = self.good(fields, rid, certain=model != HAIKU)
            if model == SONNET and call == 2 and rid == 0:
                a.update(source="poultry_meat", source_quote="chicken meat")
            return a
        details, client = self.curate(answer)
        self.assertEqual(client.models_called(), [HAIKU, SONNET, SONNET, OPUS])
        chicken, human = details[self.keys[0]], details[self.keys[1]]
        self.assertEqual((chicken["source_term"], chicken["source_tier"], chicken["source_confidence"]),
                         ("chicken meat", "opus", "low"))
        self.assertEqual(chicken["host_tier"], "sonnet x2")  # the host was agreed, only the source disputed
        self.assertEqual(human["source_tier"], "sonnet x2")
        self.assertIn("disputed source", client.calls[3]["messages"][0]["content"])

    def test_opus_unknown_or_unquoted_is_withheld(self):
        for opus_answer in ("unknown", "bad quote"):
            def answer(model, rid, fields, call, opus_answer=opus_answer):
                a = self.good(fields, rid, certain=False)
                if model == SONNET and call == 2:
                    a.update(source="meat", source_quote="meat")
                if model == OPUS:
                    a.update(source="unknown", source_field="", source_quote="") if opus_answer == "unknown" \
                        else a.update(source_quote="not in the record")
                return a
            self.cache.unlink(missing_ok=True)
            details, _ = self.curate(answer)
            chicken = details[self.keys[0]]
            self.assertEqual((chicken["source_term"], chicken["source_tier"], chicken["niche"]),
                             ("unknown", "withheld", "unknown"), opus_answer)

    def test_answers_are_cached_across_runs(self):
        self.curate(lambda m, rid, f, c: self.good(f, rid))
        details, client = self.curate(lambda m, rid, f, c: self.fail("cached records must not be re-asked"))
        self.assertEqual(client.calls, [])
        self.assertEqual(details[self.keys[0]]["source_tier"], "haiku")
        self.assertIn("raw chicken meat", self.cache.read_text())

    def test_failed_requests_are_withheld_but_not_cached(self):
        stop = lambda model, n: "refusal" if model == OPUS else "end_turn"  # noqa: E731
        def answer(model, rid, fields, call):
            a = self.good(fields, rid, certain=False)
            if model == SONNET and call == 2:
                a.update(source="meat", source_quote="meat", host="unknown", host_field="", host_quote="")
            return a
        details, _ = self.curate(answer, stop)
        self.assertEqual(details[self.keys[0]]["source_tier"], "withheld")
        details, client = self.curate(lambda m, rid, f, c: self.good(f, rid))
        self.assertEqual(details[self.keys[0]]["source_tier"], "haiku")  # asked again, not served from the cache

    def test_max_tokens_splits_the_chunk(self):
        stop = lambda model, n: "max_tokens" if n > 1 else "end_turn"  # noqa: E731
        details, client = self.curate(lambda m, rid, f, c: self.good(f, rid), stop)
        self.assertEqual(client.models_called(), [HAIKU, HAIKU, HAIKU])
        self.assertEqual(details[self.keys[1]]["source_tier"], "haiku")

    def test_answers_are_limited_to_the_vocabulary(self):
        _, client = self.curate(lambda m, rid, f, c: self.good(f, rid))
        item = client.calls[0]["output_config"]["format"]["schema"]["properties"]["labels"]["items"]
        self.assertIn("chicken_meat", item["properties"]["source"]["enum"])
        self.assertIn("no_host", item["properties"]["host"]["enum"])
        self.assertEqual(item["properties"]["source_field"]["enum"], ["isolation_source", ""])


class CsvTest(unittest.TestCase):
    """The whole step on a fake run: niche_metadata's outputs keep their columns and gain the curated ones."""

    def test_genome_metadata_is_extended_not_changed(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        out = Path(tmp.name) / "work" / "Species" / "ecosim_output_x"
        out.mkdir(parents=True)
        meta = {"GCA_1.1": {"isolation_source": "blood", "host": "Homo sapiens", "geo_loc_name": "USA: CA",
                            "collection_date": "2019"},
                "GCA_2.1": {"isolation_source": "cheese", "geo_loc_name": "france", "collection_date": "2020-05"},
                "GCA_3.1": {"isolation_source": "blood", "host": "Homo sapiens"},
                "GCA_4.1": {"isolation_source": "soil", "geo_loc_name": "Atlantis"}}
        with open(out / "ecotype_membership.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["file", "gene_count", "trial", "ecotype", "taxon"])
            for acc, eco in zip(meta, ["1", "2", "1", "2"]):
                writer.writerow([FULL_CORE_XML, "", "", eco, acc])
        client = FakeClient(lambda *a: self.fail("everything here is settled by the rules"))
        run_curated_niche_analysis(out, "Species", client=client, fetcher=lambda accs: {a: meta[a] for a in accs})

        with open(out / "genome_metadata.csv", newline="") as f:
            reader = csv.DictReader(f)
            rows = {r["accession"]: r for r in reader}
            header = reader.fieldnames
        self.assertEqual(header[:6], ["accession", "ecotype", "host", "niche", "country", "year"])
        self.assertEqual(header[6:6 + len(LLM_FIELDS)], [f"raw_{k}" for k in LLM_FIELDS])
        self.assertEqual(header[6 + len(LLM_FIELDS):], list(CURATION_COLUMNS))
        blood, cheese, soil = rows["GCA_1.1"], rows["GCA_2.1"], rows["GCA_4.1"]
        self.assertEqual((blood["host"], blood["niche"], blood["source_category"], blood["body_site"]),
                         ("human", "blood", "human-associated", "blood [UBERON:0000178]"))
        self.assertEqual((blood["mixs_extension"], blood["host_taxid"], blood["source_confidence"]),
                         ("HumanAssociated", "9606", "high"))
        self.assertEqual((cheese["host"], cheese["niche"], cheese["food"]),
                         ("no host (food or environment)", "dairy", "cheese [FOODON:00001013]"))
        self.assertEqual((cheese["country"], cheese["country_insdc"]), ("france", "France"))
        self.assertEqual((soil["country_insdc"], soil["vocabulary_version"]), ("", Vocabulary().version))

        with open(out / "ecotype_niche_association.csv", newline="") as f:
            variables = [r["variable"] for r in csv.DictReader(f)]
        self.assertEqual(variables, ["host", "niche", "country", "source_category"])
        self.assertTrue((out / "ecotype_x_source_category.csv").exists())
        self.assertEqual(client.calls, [])


if __name__ == "__main__":
    unittest.main()
