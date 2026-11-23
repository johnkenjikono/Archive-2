"""Run from the repo root: python -m unittest steps.test_delimitation
Tool-backed tests skip when TreeCluster.py, mptp or Rscript+fastbaps is not on PATH."""
import csv
import random
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from steps.delimitation import parse_mptp, run_fastbaps, run_mptp, run_treecluster, write_membership
from steps.niche_metadata import load_membership

# Three clades of 4 tips, far apart, plus a tip named outgroup.
CLADE = "({0}0:0.001,({0}1:0.001,({0}2:0.001,{0}3:0.002):0.001):0.001)"
TREE = f"((({CLADE.format('A')}:0.2,{CLADE.format('B')}:0.2):0.1,{CLADE.format('C')}:0.3):0.5,outgroup:1.0);"


def has_fastbaps():
    return shutil.which("Rscript") and not subprocess.run(
        ["Rscript", "-e", "quit(status = !requireNamespace('fastbaps', quietly = TRUE))"]).returncode


def same_partition(test, membership_csv):
    _, ecotype_of = load_membership(membership_csv)
    test.assertNotIn("outgroup", ecotype_of)
    test.assertEqual(len(ecotype_of), 12)
    for clade in "ABC":
        test.assertEqual(len({ecotype_of[f"{clade}{i}"] for i in range(4)}), 1, clade)
    test.assertEqual(len({ecotype_of[f"{c}0"] for c in "ABC"}), 3)


class WriteMembershipTest(unittest.TestCase):
    def test_expands_clones_and_drops_outgroup(self):
        with tempfile.TemporaryDirectory() as d:
            clones = Path(d) / "clone_groups.csv"
            clones.write_text("representative,member\nA0,A0\nA0,A9\n")
            out = write_membership(Path(d) / "m.csv", "x", {"A0": "1", "B0": "2", "outgroup": "3"}, clones, "outgroup")
            with open(out, newline="") as f:
                rows = {r["taxon"]: r for r in csv.DictReader(f)}
            self.assertEqual(set(rows), {"A0", "A9", "B0"})
            self.assertEqual((rows["A9"]["representative"], rows["A9"]["ecotype_size"]), ("A0", "2"))

    def test_parse_mptp(self):
        with tempfile.TemporaryDirectory() as d:
            txt = Path(d) / "mptp.txt"
            txt.write_text("Number of delimited species: 2\n\nSpecies 1:\nA0\nA1\n\nSpecies 2:\nB0\n")
            self.assertEqual(parse_mptp(txt), {"A0": "1", "A1": "1", "B0": "2"})


class ToolTest(unittest.TestCase):
    def setUp(self):
        self.dir = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.dir)
        (self.dir / "tree.nwk").write_text(TREE)

    @unittest.skipUnless(shutil.which("TreeCluster.py"), "TreeCluster.py not on PATH")
    def test_treecluster(self):
        same_partition(self, run_treecluster(self.dir / "tree.nwk", self.dir, 0.05, "outgroup"))

    @unittest.skipUnless(shutil.which("mptp"), "mptp not on PATH")
    def test_mptp(self):
        same_partition(self, run_mptp(self.dir / "tree.nwk", self.dir, "outgroup"))

    @unittest.skipUnless(has_fastbaps(), "Rscript with fastbaps not available")
    def test_fastbaps(self):
        rng = random.Random(1)
        seqs = {c: [rng.choice("ACGT") for _ in range(300)] for c in "ABC"}
        with open(self.dir / "aln.fasta", "w") as f:
            for c in "ABC":
                for i in range(4):
                    s = list(seqs[c])
                    s[rng.randrange(300)] = rng.choice("ACGT")  # a little within-clade variation
                    f.write(f">{c}{i}\n{''.join(s)}\n")
        same_partition(self, run_fastbaps(self.dir / "aln.fasta", self.dir))


if __name__ == "__main__":
    unittest.main()
