"""Run from the repo root: python -m unittest steps.test_outgroup"""
import random
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from steps import download_outgroup
from steps.download_outgroup import select_outgroup

INGROUP = [Path(f"/x/g{i}.fna") for i in (1, 2, 3)]
WITHIN = {("g1", "g2"): 99.5, ("g1", "g3"): 99.0, ("g2", "g3"): 99.2}  # ingroup minimum 99.0
CANDIDATES = {"GCF_A": "Genus closest", "GCF_B": "Genus farther"}


def run(ani, sample_size=20):
    """select_outgroup with NCBI and fastANI faked; ani maps (candidate, genome) to ANI, absent = unreported."""
    table = {frozenset(k): v for k, v in {**WITHIN, **ani}.items()}

    def fake_fastani(queries, out_tsv, threads=8, refs=None):
        with open(out_tsv, "w") as f:
            for q in queries:
                for r in refs or queries:
                    if (v := table.get(frozenset((Path(q).stem, Path(r).stem)))) is not None:
                        f.write(f"{q}\t{r}\t{v}\t10\t10\n")
        return out_tsv

    with tempfile.TemporaryDirectory() as d, \
            mock.patch.object(download_outgroup, "outgroup_candidates", return_value=CANDIDATES), \
            mock.patch.object(download_outgroup, "download_genomes",
                              side_effect=lambda accs, out: {a: Path(out) / f"{a}.fna" for a in accs}), \
            mock.patch.object(download_outgroup, "run_fastani", side_effect=fake_fastani), \
            mock.patch.object(download_outgroup, "ANI_SAMPLE", sample_size):
        return select_outgroup("Genus species", INGROUP, d)


def both(a, b):
    """ANI of candidates A and B to every ingroup genome."""
    return {**{("GCF_A", f"g{i}"): a for i in (1, 2, 3)}, **{("GCF_B", f"g{i}"): b for i in (1, 2, 3)}}


class SelectOutgroupTest(unittest.TestCase):
    def test_closest_passing_wins(self):
        acc, org, sel = run(both(98.5, 90.0))
        self.assertEqual((acc, org), ("GCF_A", "Genus closest"))
        self.assertEqual((sel["ingroup_min_ani"], sel["chosen_max_ani"]), (99.0, 98.5))
        self.assertEqual(sel["candidate_mean_ani"]["GCF_B"]["mean_ani_to_sample"], 90.0)

    def test_closest_failing_is_skipped(self):
        ani = both(98.5, 90.0)
        ani[("GCF_A", "g1")] = 99.1  # closer to g1 than g1 and g3 are to each other
        self.assertEqual(run(ani)[0], "GCF_B")

    def test_unreported_ani_is_not_assumed_distant(self):
        ani = both(98.5, 90.0)
        del ani[("GCF_A", "g2")]  # fastANI silent: below ~80%, i.e. not measured
        self.assertEqual(run(ani)[0], "GCF_B")

    def test_none_passing_raises(self):
        with self.assertRaisesRegex(RuntimeError, "minimum pairwise ANI"):
            run(both(99.3, 99.6))
        with self.assertRaisesRegex(RuntimeError, "none within fastANI range"):
            run({})  # every candidate too divergent to measure

    def test_sample_pass_is_verified_against_every_ingroup_genome(self):
        left_out = (set(INGROUP) - set(random.Random(42).sample(sorted(INGROUP), 2))).pop().stem
        ani = both(98.5, 90.0)
        ani[("GCF_A", left_out)] = 99.9  # only the full check sees this
        self.assertEqual(run(ani, sample_size=2)[0], "GCF_B")


if __name__ == "__main__":
    unittest.main()
