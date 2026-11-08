"""Run from the repo root: python -m unittest steps.test_recombination"""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pipeline import run_tree_pipeline
from steps.recombination import _segments, _to_original, genome_order, mask_recombinant_regions

# Panaroo order (file names) is A, B, C; in the genomes the genes run C, A, B.
GENES = [("geneA", 4), ("geneB", 6), ("geneC", 5)]
SEQS = {
    "outgroup": "TTTTTTTTTTTTTTT",
    "g1": "ACGTACGTACGTACG",
    "g2": "AAGTACGTACGTACG",
    "g3": "ACCTACGTACGTACG",
    "g4": "ACGRACGTACGTACC",  # R: an IUPAC code Gubbins would reject
}
PRESENCE = {
    "geneA": {"outgroup": "OUT_00001", "g1": "AAA_00020", "g2": "BBB_00020", "g3": "CCC_00020", "g4": "DDD_00020"},
    "geneB": {"outgroup": "OUT_00002", "g1": "AAA_00030", "g2": "BBB_00030", "g3": "CCC_00030", "g4": "DDD_00030"},
    # g2's copy was re-found by Panaroo (no position) and g3 has two, so g1 and g4 place every gene.
    "geneC": {"outgroup": "OUT_00003", "g1": "AAA_00010", "g2": "7_refound_3", "g3": "CCC_00005;CCC_00090",
              "g4": ""},
}

FAKE_GUBBINS = """
import json, shutil, sys
from pathlib import Path
config = json.loads(Path(sys.argv[1]).read_text())
args = sys.argv[2:]
alignment, prefix = args[0], args[args.index("--prefix") + 1]
with open(config["record"], "a") as log:
    log.write(json.dumps(args) + "\\n")
shutil.copy(alignment, config["record"] + ".input.fasta")
names = [line[1:].strip() for line in open(alignment) if line.startswith(">")]
with open(prefix + ".recombination_predictions.gff", "w") as gff:
    gff.write("##gff-version 3\\n##sequence-region SEQUENCE 1 15\\n")
    for start, end, taxa in config["blocks"]:
        gff.write(f'SEQUENCE\\tGUBBINS\\tCDS\\t{start}\\t{end}\\t0.000\\t.\\t0\\tnode="internal_1";'
                  f'neg_log_likelihood="1.0";taxa="  {" ".join(taxa)}";snp_count="3";\\n')
Path(prefix + ".per_branch_statistics.csv").write_text("Node\\tTotal SNPs\\n")
kept = [n for n in names if n not in config.get("drop", [])]
Path(prefix + ".final_tree.tre").write_text("(" + ",".join(n + ":0.1" for n in kept) + ");\\n")
"""


def write_panaroo_results(folder, seqs=SEQS):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "core_gene_alignment.aln").write_text("".join(f">{n}\n{s}\n" for n, s in seqs.items()))
    lines, start = ["ID   Genome standard; DNA; PRO; 1234 BP.", "XX", "FH   Key             Location/Qualifiers", "FH"], 1
    for name, length in GENES:
        lines += [f"FT   feature         {start}..{start + length - 1}",
                  f"FT                   /label={name}.aln", f"FT                   /locus_tag={name}.aln"]
        start += length
    (folder / "core_alignment_header.embl").write_text("\n".join(lines) + "\n")
    genomes = list(seqs)
    with open(folder / "gene_presence_absence.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Gene", "Non-unique Gene name", "Annotation", *genomes])
        for gene, cells in PRESENCE.items():
            writer.writerow([gene, "", "", *(cells.get(g, "") for g in genomes)])
    return folder / "core_gene_alignment.aln"


def fasta(path):
    out, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
        else:
            out[name] = line
    return out


class FakeGubbins:
    """A run_gubbins.py stand-in: writes Gubbins-style outputs for the given blocks into its cwd."""

    def __init__(self, folder, blocks, drop=()):
        self.record = folder / "gubbins_calls.txt"
        script = folder / "fake_run_gubbins.py"
        script.write_text(FAKE_GUBBINS)
        config = folder / "fake_gubbins.json"
        config.write_text(json.dumps({"record": str(self.record), "blocks": blocks, "drop": list(drop)}))
        self.command = [sys.executable, str(script), str(config)]

    def calls(self):
        return [json.loads(line) for line in self.record.read_text().splitlines()] if self.record.exists() else []

    def input(self):
        return fasta(f"{self.record}.input.fasta")


class GenomeOrderTest(unittest.TestCase):
    def test_orders_genes_by_the_best_placed_genome(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            write_panaroo_results(folder)
            reference, order, unplaced = genome_order([g for g, _ in GENES], folder / "gene_presence_absence.csv",
                                                      ["g1", "g2", "g3", "g4"])
        self.assertEqual(reference, "g1")  # g1 and g3 place all three, but g3 has a two-copy gene
        self.assertEqual(order, ["geneC", "geneA", "geneB"])
        self.assertEqual(unplaced, 0)

    def test_unknown_gene_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            write_panaroo_results(folder)
            with self.assertRaises(ValueError):
                genome_order(["geneA", "geneZ"], folder / "gene_presence_absence.csv", ["g1"])

    def test_blocks_map_back_across_gene_boundaries(self):
        genes = [("geneA", 0, 4), ("geneB", 4, 10), ("geneC", 10, 15)]
        segments = _segments(genes, ["geneC", "geneA", "geneB"])
        starts = [s[0] for s in segments]
        # Reordered [3, 11) is the end of C, all of A and the start of B.
        self.assertEqual(_to_original(3, 11, segments, starts), [(13, 15), (0, 4), (4, 6)])


class MaskTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.tmp = Path(tmp.name)
        self.results = self.tmp / "results"
        aln = write_panaroo_results(self.results)
        self.fasta = self.tmp / "sorted.fasta"
        self.fasta.write_text(aln.read_text())
        self.out = self.tmp / "ecosim_out"

    def mask(self, gubbins, seqs_path=None):
        return mask_recombinant_regions(seqs_path or self.fasta, self.results / "core_alignment_header.embl",
                                        self.results / "gene_presence_absence.csv", self.out, self.tmp / "work",
                                        "outgroup", threads=2, command=gubbins.command)

    def test_masks_blocks_in_original_gene_order(self):
        # Gubbins coordinates are in genome order C(1-5) A(6-9) B(10-15).
        gubbins = FakeGubbins(self.tmp, [[6, 9, ["g1", "g2"]], [4, 11, ["g3"]]])
        result = self.mask(gubbins)
        masked = fasta(self.fasta)
        self.assertEqual(masked["g1"], "----" + SEQS["g1"][4:])
        self.assertEqual(masked["g3"], "------" + SEQS["g3"][6:13] + "--")
        self.assertEqual(masked["g4"], SEQS["g4"])  # untouched, IUPAC code kept
        self.assertEqual(masked["outgroup"], SEQS["outgroup"])
        self.assertEqual(result["blocks"], 2)
        self.assertEqual(result["genes_emptied"], 1)  # geneA, all gaps in g1, g2 and g3

        sent = gubbins.input()
        self.assertNotIn("outgroup", sent)  # the outgroup never goes to Gubbins
        self.assertEqual(sent["g1"], SEQS["g1"][10:] + SEQS["g1"][:4] + SEQS["g1"][4:10])  # C, A, B
        self.assertEqual(sent["g4"][5 + 3], "N")  # the R, now behind geneC's 5 bases
        args = gubbins.calls()[0]
        self.assertEqual(args[args.index("--filter-percentage") + 1], "100")

        with open(self.out / "recombination_masking.csv", newline="") as handle:
            rows = {row["genome"]: row for row in csv.DictReader(handle)}
        self.assertEqual((rows["g1"]["blocks"], rows["g1"]["masked_sites"]), ("1", "4"))
        self.assertEqual(rows["g3"]["masked_sites"], "8")
        self.assertEqual(rows["g4"]["masked_sites"], "0")
        self.assertEqual((rows["all"]["blocks"], rows["all"]["masked_sites"]), ("2", "8"))
        self.assertTrue((self.out / "recombination_predictions.gff").exists())
        self.assertTrue((self.out / "recombination_per_branch.tsv").exists())

    def test_genome_left_out_of_gubbins_tree_fails(self):
        before = self.fasta.read_text()
        self.assertIsNone(self.mask(FakeGubbins(self.tmp, [[6, 9, ["g1"]]], drop=["g4"])))
        self.assertEqual(self.fasta.read_text(), before)

    def test_reuses_gubbins_results_for_the_same_alignment(self):
        gubbins = FakeGubbins(self.tmp, [[6, 9, ["g1"]]])
        original = self.fasta.read_text()
        self.mask(gubbins)
        self.fasta.write_text(original)
        self.mask(gubbins)
        self.assertEqual(len(gubbins.calls()), 1)

    def test_fewer_than_three_ingroup_genomes_skips_gubbins(self):
        small = self.tmp / "small.fasta"
        small.write_text("".join(f">{n}\n{SEQS[n]}\n" for n in ("outgroup", "g1", "g2")))
        gubbins = FakeGubbins(self.tmp, [[6, 9, ["g1"]]])
        self.assertEqual(self.mask(gubbins, small)["blocks"], 0)
        self.assertEqual(gubbins.calls(), [])


class PipelineStepTest(unittest.TestCase):
    """Step 4 of run_tree_pipeline with the fake Gubbins."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.base = Path(tmp.name)
        self.alignment = write_panaroo_results(self.base / "output_roary" / "results")
        self.gubbins = FakeGubbins(self.base, [[6, 9, ["g1"]]])
        self.temp = self.base / "pipeline_temp_core_gene_alignment"
        self.out = self.base / "ecosim_output_core_gene_alignment"

    def run_step_4(self, mask):
        with mock.patch("steps.recombination.gubbins_command", return_value=self.gubbins.command):
            run_tree_pipeline(str(self.alignment), 4, "outgroup", str(self.base), threads=1, stop_step=4,
                              mask_recombination=mask)
        return fasta(self.temp / "tree_rdy_fastas" / "core_gene_alignment_sorted.fasta")

    def cached_tree(self):
        tree = self.temp / "trees_final" / "core_gene_alignment_sorted_tree.nwk"
        tree.parent.mkdir(parents=True, exist_ok=True)
        tree.write_text("(outgroup:1,g1:1,g2:1,g3:1,g4:1);\n")
        return tree

    def test_step_4_masks_after_collapsing_clones(self):
        sorted_fasta = self.run_step_4(mask=True)
        self.assertEqual(list(sorted_fasta)[0], "outgroup")
        self.assertEqual(sorted_fasta["g1"][:4], "----")
        self.assertTrue((self.out / "recombination_masking.csv").exists())

    def test_switching_masking_off_drops_cached_trees_and_outputs(self):
        self.run_step_4(mask=True)
        tree = self.cached_tree()
        self.assertEqual(self.run_step_4(mask=False)["g1"], SEQS["g1"])
        self.assertFalse(tree.exists())
        self.assertFalse((self.out / "recombination_masking.csv").exists())

    def test_switching_masking_on_drops_cached_trees(self):
        self.run_step_4(mask=False)
        tree = self.cached_tree()
        self.run_step_4(mask=True)
        self.assertFalse(tree.exists())

    def test_same_setting_keeps_cached_trees(self):
        self.run_step_4(mask=True)
        tree = self.cached_tree()
        self.run_step_4(mask=True)
        self.assertTrue(tree.exists())
        self.assertEqual(len(self.gubbins.calls()), 1)  # Gubbins results reused as well

    def test_missing_gene_table_stops_before_step_4(self):
        (self.alignment.parent / "gene_presence_absence.csv").unlink()
        with self.assertRaises(SystemExit):
            self.run_step_4(mask=True)


class CommandLineTest(unittest.TestCase):
    def test_flag_without_gubbins_is_refused_up_front(self):
        repo = Path(__file__).resolve().parent.parent
        with tempfile.TemporaryDirectory() as tmp:
            # Only the interpreter's folder on PATH: no conda, no run_gubbins.py.
            env = {**os.environ, "PATH": str(Path(sys.executable).parent)}
            proc = subprocess.run([sys.executable, str(repo / "pipeline.py"), "--species", "Test species",
                                   "--workdir", tmp, "--start-step", "4", "--mask-recombination"],
                                  capture_output=True, text=True, env=env)
        self.assertEqual(proc.returncode, 2)
        self.assertIn("needs Gubbins", proc.stderr)


if __name__ == "__main__":
    unittest.main()
