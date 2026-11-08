"""Run from the repo root: python -m unittest steps.test_oversample"""
import csv
import os
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from pipeline import oversampled_size, run_tree_pipeline
from steps.dedup_clones import trim_representatives

CLONE_MAP_HEADER = ["representative", "member", "group_size", "differences", "divergence"]


def write_fasta(path, ids):
    Path(path).write_text("".join(f">{i}\nACGT\n" for i in ids))


def fasta_ids(path):
    return [line[1:] for line in Path(path).read_text().splitlines() if line.startswith(">")]


def write_clone_map(path, groups):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(CLONE_MAP_HEADER)
        for rep, members in groups.items():
            for member in members:
                writer.writerow([rep, member, len(members), 0, "0.00e+00"])


def clone_map_members(path):
    with open(path, newline="") as handle:
        return {row["member"]: row["representative"] for row in csv.DictReader(handle)}


class OversampledSizeTest(unittest.TestCase):
    def test_sizes(self):
        for sample_size, percent, want in [(200, 15, 230), (10, 15, 12), (200, 0, 200), (1, 15, 2)]:
            self.assertEqual(oversampled_size(sample_size, percent), want)


class TrimRepresentativesTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.fasta = Path(tmp.name) / "sorted.fasta"
        self.clone_map = Path(tmp.name) / "clone_groups.csv"

    def test_keeps_root_plus_sample_size_in_input_order(self):
        ids = ["outgroup", "a", "b", "c", "d", "e"]
        write_fasta(self.fasta, ids)
        write_clone_map(self.clone_map, {i: [i] for i in ids})

        trim_representatives(self.fasta, self.clone_map, 3, "outgroup")

        kept = fasta_ids(self.fasta)
        self.assertEqual(len(kept), 4)
        self.assertEqual(kept[0], "outgroup")
        self.assertEqual(kept, [i for i in ids if i in kept])

    def test_dropping_a_representative_drops_its_clones_from_the_map(self):
        groups = {"outgroup": ["outgroup"], "a": ["a", "a2", "a3"], "b": ["b", "b2"], "c": ["c"]}
        write_fasta(self.fasta, list(groups))
        write_clone_map(self.clone_map, groups)

        trim_representatives(self.fasta, self.clone_map, 2, "outgroup")

        members = clone_map_members(self.clone_map)
        kept = set(fasta_ids(self.fasta))
        self.assertEqual(set(members.values()), kept)
        want = {m for rep in kept for m in groups[rep]}
        self.assertEqual(set(members), want)

    def test_nothing_changes_at_or_under_the_limit(self):
        ids = ["outgroup", "a", "b"]
        write_fasta(self.fasta, ids)
        write_clone_map(self.clone_map, {i: [i] for i in ids})

        trim_representatives(self.fasta, self.clone_map, 2, "outgroup")

        self.assertEqual(fasta_ids(self.fasta), ids)
        self.assertEqual(len(clone_map_members(self.clone_map)), 3)

    def test_same_seed_same_genomes_and_seeds_differ(self):
        ids = ["outgroup"] + [f"g{i}" for i in range(30)]

        def kept(seed):
            write_fasta(self.fasta, ids)
            write_clone_map(self.clone_map, {i: [i] for i in ids})
            trim_representatives(self.fasta, self.clone_map, 10, "outgroup", seed=seed)
            return fasta_ids(self.fasta)

        self.assertEqual(kept(1), kept(1))
        self.assertNotEqual(kept(1), kept(2))

    def test_first_record_is_root_without_an_outgroup_id(self):
        ids = ["first", "a", "b", "c"]
        write_fasta(self.fasta, ids)
        write_clone_map(self.clone_map, {i: [i] for i in ids})

        trim_representatives(self.fasta, self.clone_map, 1, None)

        kept = fasta_ids(self.fasta)
        self.assertEqual(len(kept), 2)
        self.assertEqual(kept[0], "first")


class Step4Test(unittest.TestCase):
    """Real steps 4 and 4b: 5 genomes, 2 of them clones, so 4 distinct sequences remain."""

    def run_step_4(self, sample_size, dedup=True):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        base = Path(tmp.name)
        alignment = base / "core_gene_alignment.aln"
        alignment.write_text(">outgroup\nAAAAAAAA\n>a\nCCCCCCCC\n>a2\nCCCCCCCC\n>b\nGGGGGGGG\n>c\nTTTTTTTT\n")
        (base / "core_alignment_header.embl").write_text("")
        run_tree_pipeline(str(alignment), 4, "outgroup", str(base), threads=1, dedup=dedup, stop_step=4,
                          sample_size=sample_size)
        folder = base / "pipeline_temp_core_gene_alignment" / "tree_rdy_fastas"
        return fasta_ids(folder / "core_gene_alignment_sorted.fasta"), base

    def test_trims_after_collapsing_clones(self):
        kept, base = self.run_step_4(sample_size=2)
        self.assertEqual(len(kept), 3)
        self.assertEqual(kept[0], "outgroup")
        clone_map = base / "ecosim_output_core_gene_alignment" / "clone_groups.csv"
        self.assertEqual(set(clone_map_members(clone_map).values()), set(kept))

    def test_clones_do_not_count_towards_the_sample(self):
        kept, _ = self.run_step_4(sample_size=3)
        self.assertEqual(sorted(kept), ["a", "b", "c", "outgroup"])

    def test_no_sample_size_keeps_every_representative(self):
        kept, _ = self.run_step_4(sample_size=None)
        self.assertEqual(sorted(kept), ["a", "b", "c", "outgroup"])


class CommandLineTest(unittest.TestCase):
    """pipeline.py as a user runs it: which flags switch the trim and the extra download on."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.work = Path(tmp.name)

    def pipeline(self, *flags, env=None):
        repo = Path(__file__).resolve().parent.parent
        subprocess.run([sys.executable, str(repo / "pipeline.py"), "--species", "Test species", "--workdir", str(self.work),
                        *flags], check=True, capture_output=True, text=True, env=env)

    def run_from_step_4(self, *flags):
        results = self.work / "Test_species" / "output_roary" / "results"
        results.mkdir(parents=True)
        alignment = results / "core_gene_alignment.aln"
        alignment.write_text(">outgroup\nAAAA\n>a\nCCCC\n>b\nGGGG\n>c\nTTTT\n")
        (results / "core_alignment_header.embl").write_text("")
        self.pipeline("--start-step", "4", "--stop-step", "4", "--no-auto-outgroup", "--outgroup-id", "outgroup",
                      "--sample-size", "1", *flags)
        return fasta_ids(self.work / "Test_species" / "pipeline_temp_core_gene_alignment" / "tree_rdy_fastas"
                         / "core_gene_alignment_sorted.fasta")

    def test_resume_at_step_4_keeps_every_genome(self):
        # Step 1 did not run here, so nothing was oversampled and there is nothing to trim back.
        self.assertEqual(len(self.run_from_step_4()), 4)

    def test_oversample_zero_keeps_every_genome(self):
        self.assertEqual(len(self.run_from_step_4("--oversample", "0")), 4)

    def test_no_dedup_keeps_every_genome(self):
        self.assertEqual(len(self.run_from_step_4("--no-dedup")), 4)

    def test_input_fasta_is_not_trimmed(self):
        own = self.work / "own"
        own.mkdir()
        (own / "core_gene_alignment.aln").write_text(">outgroup\nAAAA\n>a\nCCCC\n>b\nGGGG\n>c\nTTTT\n")
        (own / "core_alignment_header.embl").write_text("")
        kept = self.run_from_step_4("--input-fasta", str(own / "core_gene_alignment.aln"))
        self.assertEqual(len(kept), 4)

    def downloaded_count(self, *flags):
        """Step 1 against a fake `datasets` CLI listing 300 assemblies; returns how many were requested."""
        bin_dir = self.work / "bin"
        bin_dir.mkdir()
        fake = bin_dir / "datasets"
        record = self.work / "requested.txt"
        fake.write_text(f"""#!/bin/sh
if [ "$1" = summary ]; then
  i=0; while [ $i -lt 300 ]; do echo '{{"accession":"GCA_'$i'.1"}}'; i=$((i+1)); done
else
  while [ $# -gt 0 ]; do [ "$1" = --inputfile ] && cp "$2" "{record}"; shift; done
  exit 1
fi
""")
        fake.chmod(fake.stat().st_mode | stat.S_IEXEC)
        env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}"}
        with self.assertRaises(subprocess.CalledProcessError):  # the fake download fails, after recording the list
            self.pipeline("--stop-step", "1", "--sample-size", "200", *flags, env=env)
        return len(record.read_text().split())

    def test_step_1_requests_the_oversampled_count(self):
        self.assertEqual(self.downloaded_count(), 230)

    def test_no_dedup_requests_only_the_sample_size(self):
        self.assertEqual(self.downloaded_count("--no-dedup"), 200)


if __name__ == "__main__":
    unittest.main()
