"""Run from the repo root: python -m unittest steps.test_adaptml"""
import csv
import json
import random
import tempfile
import unittest
from pathlib import Path

from steps.adaptml import run_adaptml
from steps.niche_metadata import load_membership

N = 30  # tips per clade


def clade(name, rng):
    """N tips joined in random order by short branches."""
    parts = [f"{name}{i}:0.002" for i in range(N)]
    while len(parts) > 1:
        a, b = parts.pop(rng.randrange(len(parts))), parts.pop(rng.randrange(len(parts)))
        parts.append(f"({a},{b}):0.002")
    return parts[0].rsplit(":", 1)[0]


def run(folder, labels, seed=42):
    """Two clades on long stems plus an outgroup; labels = {tip: environment}. Returns (run json, genome rows)."""
    rng = random.Random(1)
    folder = Path(folder)
    (folder / "tree.nwk").write_text(f"(({clade('a', rng)}:0.05,{clade('b', rng)}:0.05):0.1,outgroup:0);\n")
    with open(folder / "genome_metadata.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["accession", "ecotype", "host", "niche"])
        writer.writerows([tip, "1", "human", env] for tip, env in labels.items())
    out = run_adaptml(folder / "tree.nwk", folder / "genome_metadata.csv", folder / "out", outgroup="outgroup",
                      seed=seed)
    with open(folder / "out" / "adaptml_genomes.csv", newline="") as f:
        rows = {r["taxon"]: r for r in csv.DictReader(f)}
    return json.loads((folder / "out" / "adaptml_run.json").read_text()), rows, out


def two_niches(noise, rng):
    """Clade a in blood, clade b on skin, each tip swapped with probability `noise`."""
    return {f"{c}{i}": (home if rng.random() >= noise else other)
            for c, home, other in (("a", "blood", "skin"), ("b", "skin", "blood")) for i in range(N)}


class AdaptMLTest(unittest.TestCase):
    def test_two_clades_in_two_environments(self):
        labels = two_niches(0.1, random.Random(2))
        del labels["a3"]            # no row in the CSV
        labels["b4"] = "unknown"    # niche_metadata.py's missing value
        with tempfile.TemporaryDirectory() as tmp:
            info, rows, membership = run(tmp, labels)
            self.assertEqual(info["habitats"], 2)
            self.assertEqual(info["populations"], 2)
            self.assertNotIn("outgroup", rows)
            habitat_a = {rows[f"a{i}"]["habitat"] for i in range(N)}
            habitat_b = {rows[f"b{i}"]["habitat"] for i in range(N)}
            self.assertEqual((len(habitat_a), len(habitat_b)), (1, 1))  # unlabelled tips follow their clade
            self.assertNotEqual(habitat_a, habitat_b)
            # The membership file reads like EcoSim's, and the populations are the two clades.
            _, population_of = load_membership(membership)
            self.assertEqual(len(population_of), 2 * N)
            self.assertEqual(len({population_of[f"a{i}"] for i in range(N)}), 1)
            self.assertNotEqual(population_of["a0"], population_of["b0"])
            with open(Path(tmp) / "out" / "adaptml_habitats.csv", newline="") as f:
                habitats = list(csv.DictReader(f))
            low, high = sorted(float(h["blood"]) for h in habitats)
            self.assertTrue(low < 0.2 and high > 0.8, (low, high))

    def test_shuffled_environments_rarely_give_clade_structure(self):
        # Populations are tested at the 0.95 level, so a chance imbalance between the two clades is
        # sometimes picked up: bound how often, rather than demand never.
        structured = 0
        for seed in range(20):
            labels = two_niches(0.1, random.Random(2))
            values = list(labels.values())
            random.Random(seed).shuffle(values)
            with tempfile.TemporaryDirectory() as tmp:
                _, rows, _ = run(tmp, dict(zip(labels, values)), seed=seed)
            # Structured: some habitat holds most of one clade and little of the other.
            structured += any(abs(sum(rows[f"a{i}"]["habitat"] == h for i in range(N))
                                  - sum(rows[f"b{i}"]["habitat"] == h for i in range(N))) > N / 2
                              for h in {r["habitat"] for r in rows.values()})
        self.assertLessEqual(structured, 3)

    def test_needs_two_environments(self):
        with tempfile.TemporaryDirectory() as tmp, self.assertRaises(ValueError):
            run(tmp, {f"{c}{i}": "blood" for c in "ab" for i in range(N)})


if __name__ == "__main__":
    unittest.main()
