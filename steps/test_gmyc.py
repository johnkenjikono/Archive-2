"""Run from the repo root: python -m unittest steps.test_gmyc"""
import csv
import tempfile
import unittest
from pathlib import Path

from steps.gmyc import check_r, run_gmyc
from steps.niche_metadata import load_membership

try:
    check_r()
    R_MISSING = None
except RuntimeError as e:
    R_MISSING = str(e)


def clade(name, n):
    """n tips with short, staggered coalescent-like branches below a long stem."""
    tips = [f"{name}{i}" for i in range(n)]
    newick, depth = tips[0], 0.0
    for i, tip in enumerate(tips[1:], 1):
        h = 0.001 * i
        newick = f"({newick}:{h - depth:.4f},{tip}:{h:.4f})"
        depth = h
    return newick, depth


def synthetic_tree():
    """Three clades of 6 genomes separated by deep branches, plus a tip named outgroup."""
    (a, ha), (b, hb), (c, hc) = clade("A", 6), clade("B", 6), clade("C", 6)
    ab = f"({a}:{0.3 - ha:.4f},{b}:{0.3 - hb:.4f})"
    return f"(({ab}:0.2,{c}:{0.5 - hc:.4f}):0.5,outgroup:1.0);"


@unittest.skipIf(R_MISSING, f"GMYC test needs R with ape and splits: {R_MISSING}")
class GmycTest(unittest.TestCase):
    def test_three_clades(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "tree.nwk").write_text(synthetic_tree())
            with open(d / "clone_groups.csv", "w") as f:
                f.write("representative,member\nA0,A0\nA0,A0_clone\n")
            out = run_gmyc(d / "tree.nwk", d / "gmyc", outgroup="outgroup", clone_map=d / "clone_groups.csv")

            _, cluster_of = load_membership(out)
            self.assertNotIn("outgroup", cluster_of)
            self.assertEqual(len(cluster_of), 19)  # 18 tips + the collapsed clone
            self.assertEqual(cluster_of["A0_clone"], cluster_of["A0"])
            groups = {letter: {cluster_of[f"{letter}{i}"] for i in range(6)} for letter in "ABC"}
            self.assertTrue(all(len(g) == 1 for g in groups.values()), groups)
            self.assertEqual(len(set.union(*groups.values())), 3, groups)
            with open(d / "gmyc" / "gmyc_summary.csv", newline="") as f:
                self.assertEqual(next(csv.DictReader(f))["entities"], "3")

    def test_missing_outgroup_is_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            (Path(d) / "tree.nwk").write_text(synthetic_tree())
            with self.assertRaisesRegex(RuntimeError, "not a tip"):
                run_gmyc(Path(d) / "tree.nwk", d, outgroup="nope")


if __name__ == "__main__":
    unittest.main()
