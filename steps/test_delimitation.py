"""Run from the repo root: python -m unittest steps.test_delimitation
Tool-backed tests skip when TreeCluster.py, mptp or Rscript+fastbaps is not on PATH."""
import csv
import random
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from unittest import mock

from steps import delimitation
from steps.delimitation import (TREECLUSTER_LADDER, find_memberships, parse_mptp, run_fastbaps, run_mptp,
                                run_poppunk, run_treecluster, treecluster_thresholds, write_membership)
from steps.niche_metadata import load_membership, run_niche_analysis
from steps.parsing import FULL_CORE_XML

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

    def test_treecluster_thresholds(self):
        self.assertEqual(treecluster_thresholds([]), list(TREECLUSTER_LADDER))  # --treecluster alone: the ladder
        self.assertEqual(treecluster_thresholds([0.02, 0.005, 0.02]), [0.005, 0.02])


def write_ecosim_membership(path, ecotype_of):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "gene_count", "trial", "ecotype", "ecotype_size", "taxon", "representative"])
        writer.writerows([FULL_CORE_XML, 900, 0, e, 0, a, a] for a, e in ecotype_of.items())


class PopPUNKTest(unittest.TestCase):
    def test_the_outgroup_stays_out_of_the_fit(self):
        """PopPUNK models one species: with C. muridarum in its fit it split 24 C. trachomatis into 10
        clusters, without it into 3. The outgroup must not reach --create-db at all."""
        seen = {}

        def fake_run(cmd, **kwargs):
            if "--create-db" in cmd:
                seen["rfile"] = Path(cmd[cmd.index("--r-files") + 1]).read_text()
            else:  # --fit-model: one cluster for every genome it was given
                db = Path(cmd[cmd.index("--output") + 1])
                names = [line.split("\t")[0] for line in seen["rfile"].splitlines()]
                (db / "db_clusters.csv").write_text("Taxon,Cluster\n" + "".join(f"{n},1\n" for n in names))

        with tempfile.TemporaryDirectory() as d, mock.patch.object(delimitation, "_run", fake_run):
            genomes = [Path(d) / f"{name}.fna" for name in ("GCA_1.1", "outgroup", "GCA_2.1")]
            (Path(d) / "delimitation" / "poppunk" / "db").mkdir(parents=True)
            out = run_poppunk(genomes, Path(d) / "delimitation", threads=1)
            self.assertNotIn("outgroup", seen["rfile"])
            self.assertEqual(len(seen["rfile"].splitlines()), 2)
            self.assertEqual(sorted(load_membership(out)[1]), ["GCA_1.1", "GCA_2.1"])


class FastbapsIngroupTest(unittest.TestCase):
    def test_the_outgroup_is_left_out_of_the_alignment(self):
        seen = {}

        def fake_run(cmd, **kwargs):  # Rscript -e <script> <alignment> <out csv>
            seen["fasta"] = Path(cmd[3]).read_text()
            Path(cmd[4]).write_text("genome,cluster\nGCA_1.1,1\nGCA_2.1,2\n")

        with tempfile.TemporaryDirectory() as d, mock.patch.object(delimitation, "_run", fake_run):
            aln = Path(d) / "aln.fasta"
            aln.write_text(">GCA_1.1\nACGT\n>outgroup\nTTTT\n>GCA_2.1 a description\nACGA\n")
            out = run_fastbaps(aln, d, outgroup="outgroup")
            self.assertEqual(seen["fasta"], ">GCA_1.1\nACGT\n>GCA_2.1 a description\nACGA\n")
            self.assertEqual(sorted(load_membership(out)[1]), ["GCA_1.1", "GCA_2.1"])
            run_fastbaps(aln, d)  # without an outgroup (a type strain root) every record stays
            self.assertIn(">outgroup", seen["fasta"])


class StaleOutputTest(unittest.TestCase):
    def test_a_failed_rerun_leaves_no_earlier_membership(self):
        def broken(cmd, **kwargs):
            raise RuntimeError("Rscript not found")

        with tempfile.TemporaryDirectory() as d, mock.patch.object(delimitation, "_run", broken):
            old = write_membership(Path(d) / "fastbaps_membership.csv", "fastbaps", {"GCA_1.1": "1"})
            with self.assertRaises(RuntimeError):
                run_fastbaps(Path(d) / "aln.fasta", d)
            self.assertFalse(old.exists())  # the niche tests would otherwise read the earlier run's clusters
            write_membership(Path(d) / "treecluster_t0.5_membership.csv", "treecluster_t0.5", {"GCA_1.1": "1"})
            delimitation.clear_treecluster(d)
            self.assertEqual(list(Path(d).glob("treecluster_*")), [])


class NicheTest(unittest.TestCase):
    """Every method's clusters reach the ecotype x niche tests, beside EcoSim's ecotypes."""

    def test_find_memberships_and_test_each_method(self):
        with tempfile.TemporaryDirectory() as d:
            run = Path(d)
            ecosim, delim, med = run / "ecosim_output_x", run / "delimitation", run / "med_output_x"
            for folder in (ecosim, delim, med):
                folder.mkdir()
            # Six genomes: EcoSim lumps them, TreeCluster splits blood (1) from soil (2).
            genomes = [f"GCA_{i}.1" for i in range(6)]
            write_ecosim_membership(ecosim / "ecotype_membership.csv", {a: "1" for a in genomes})
            write_membership(delim / "treecluster_t0.01_membership.csv", "treecluster_t0.01",
                             {a: "1" if i < 3 else "2" for i, a in enumerate(genomes)})
            write_membership(med / "ecotype_membership.csv", "med", {a: "1" for a in genomes})
            found = find_memberships(run, ecosim)
            self.assertEqual(sorted(found), ["med", "treecluster_t0.01"])
            # Only the methods a run asked for: another run's files in the folder are not tested again.
            self.assertEqual(sorted(find_memberships(run, ecosim, {"treecluster", "fastbaps"})), ["treecluster_t0.01"])

            def fetcher(accessions):
                return {a: {"isolation_source": "blood" if a < "GCA_3" else "soil"} for a in accessions}

            def labeler(species, keys):
                return {k: ("human", dict(k)["isolation_source"]) for k in keys}

            run_niche_analysis(ecosim, "X x", labeler=labeler, fetcher=fetcher, memberships=found)
            with open(ecosim / "method_niche_association.csv", newline="") as f:
                rows = {(r["method"], r["variable"]): r for r in csv.DictReader(f)}
            self.assertEqual({m for m, _ in rows}, {"ecosim", "treecluster_t0.01", "med"})
            split = rows[("treecluster_t0.01", "niche")]
            self.assertEqual((split["clusters"], split["cramers_v"]), ("2", "1.0"))
            self.assertEqual(rows[("ecosim", "niche")]["chi2"], "")  # one ecotype: nothing to test
            with open(ecosim / "method_niche_crosstab.csv", newline="") as f:
                table = [r for r in csv.DictReader(f) if r["method"] == "treecluster_t0.01" and r["variable"] == "niche"]
            self.assertEqual(sorted((r["cluster"], r["category"], r["genomes"]) for r in table),
                             [("1", "blood", "3"), ("2", "soil", "3")])
            with open(ecosim / "genome_metadata.csv", newline="") as f:
                self.assertEqual(len(list(csv.DictReader(f))), 6)


class ToolTest(unittest.TestCase):
    def setUp(self):
        self.dir = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.dir)
        (self.dir / "tree.nwk").write_text(TREE)

    @unittest.skipUnless(shutil.which("TreeCluster.py"), "TreeCluster.py not on PATH")
    def test_treecluster(self):
        same_partition(self, run_treecluster(self.dir / "tree.nwk", self.dir, 0.05, "outgroup"))
        # Thresholds sit side by side: at 0.0001 every tip is its own cluster.
        fine = run_treecluster(self.dir / "tree.nwk", self.dir, 0.0001, "outgroup")
        self.assertEqual(sorted(p.name for p in self.dir.glob("treecluster_*_membership.csv")),
                         ["treecluster_t0.0001_membership.csv", "treecluster_t0.05_membership.csv"])
        self.assertEqual(len(set(load_membership(fine)[1].values())), 12)

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
