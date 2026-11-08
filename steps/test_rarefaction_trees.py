"""Run from the repo root: python -m unittest steps.test_rarefaction_trees"""
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pipeline
from pipeline import run_tree_pipeline

SEQS = {"outgroup": "TTTTTTTTTT", "a": "ACGTACGTAC", "b": "AAGTACGTAC", "c": "ACCTACGTAC"}
TREE = "(outgroup:1,(a:0.1,(b:0.1,c:0.1):0.1):0.1);\n"


def fake_make_trees(final_folder, tree_folder, **_):
    """Stands in for VeryFastTree: one tree per FASTA, named as make_trees_batch names them."""
    Path(tree_folder).mkdir(parents=True, exist_ok=True)
    for fasta in Path(final_folder).glob("*.fasta"):
        (Path(tree_folder) / f"{fasta.stem}_tree.nwk").write_text(TREE)


def fake_reroot(fasta, unrooted, rooted, outgroup_name=None):
    shutil.copy(unrooted, rooted)
    return True


class RarefactionTreeModeTest(unittest.TestCase):
    """Steps 7-8 with tree building and EcoSim mocked out."""

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.base = Path(tmp.name)
        results = self.base / "output_roary" / "results"
        results.mkdir(parents=True)
        self.alignment = results / "core_gene_alignment.aln"
        self.alignment.write_text("".join(f">{n}\n{s}\n" for n, s in SEQS.items()))
        (results / "core_alignment_header.embl").write_text(
            "FT   feature         1..4\nFT                   /label=geneA.aln\n"
            "FT   feature         5..10\nFT                   /label=geneB.aln\n")
        # What steps 4-6 leave behind.
        temp = self.base / "pipeline_temp_core_gene_alignment"
        (temp / "tree_rdy_fastas").mkdir(parents=True)
        (temp / "tree_rdy_fastas" / "core_gene_alignment_sorted.fasta").write_text(self.alignment.read_text())
        (temp / "trees_final").mkdir()
        (temp / "trees_final" / "core_gene_alignment_sorted_tree.nwk").write_text(TREE)
        self.rooted = self.base / "rerooted_trees" / "core_gene_alignment.nwk"
        self.rooted.parent.mkdir()
        self.rooted.write_text(TREE)
        self.draw_trees = self.base / "rarefaction_trees_core_gene_alignment"
        self.out = self.base / "ecosim_output_core_gene_alignment"
        jar = self.base / "ecosim.jar"
        jar.write_text("")
        self.patches = [mock.patch.object(pipeline, "ECOSIM_JAR", str(jar)),
                        mock.patch.object(pipeline, "make_trees_batch", side_effect=fake_make_trees),
                        mock.patch.object(pipeline, "reroot_tree_by_first_fasta", side_effect=fake_reroot),
                        mock.patch.object(pipeline, "run_ecosim_batch", return_value=1)]

    def run_steps(self, start_step=7, **kwargs):
        with self.patches[0], self.patches[1] as make_trees, self.patches[2], self.patches[3] as ecosim:
            run_tree_pipeline(str(self.alignment), start_step, "outgroup", str(self.base), threads=1, stop_step=8,
                              gene_counts=[1], replicates=2, keep={"rarefaction", "rerooted_tree"}, **kwargs)
        return make_trees, ecosim

    def test_default_builds_one_tree_per_draw(self):
        make_trees, ecosim = self.run_steps()
        make_trees.assert_called_once()
        kwargs = ecosim.call_args.kwargs
        self.assertEqual(kwargs["tree_dir"], str(self.draw_trees))
        self.assertEqual(kwargs["full_tree_path"], str(self.rooted))
        self.assertEqual(sorted(p.name for p in self.draw_trees.glob("*.nwk")),
                         ["full_core_genome_tree.nwk", "sim_species_g1_t1_tree.nwk", "sim_species_g1_t2_tree.nwk"])
        self.assertEqual(len(list((self.base / "rarefaction_fastas_core_gene_alignment").glob("*.fasta"))), 3)

    def test_shared_tree_runs_every_draw_against_the_full_tree(self):
        make_trees, ecosim = self.run_steps(shared_tree=True)
        make_trees.assert_not_called()
        kwargs = ecosim.call_args.kwargs
        self.assertIsNone(kwargs["tree_dir"])
        self.assertEqual(kwargs["full_tree_path"], str(self.rooted))
        self.assertFalse(self.draw_trees.exists())

    def test_shared_tree_step_7_removes_per_draw_trees_from_an_earlier_run(self):
        self.run_steps()
        self.run_steps(shared_tree=True)
        self.assertFalse(self.draw_trees.exists())

    def test_shared_tree_resume_at_step_8_reruns_results_made_on_per_draw_trees(self):
        self.run_steps()
        self.out.mkdir(exist_ok=True)
        stale = self.out / "sim_species_g1_t1_results.xml"
        stale.write_text("<xml/>")
        _, ecosim = self.run_steps(start_step=8, shared_tree=True)
        self.assertIsNone(ecosim.call_args.kwargs["tree_dir"])
        self.assertFalse(stale.exists())
        self.assertFalse(self.draw_trees.exists())

    def test_default_resume_at_step_8_without_per_draw_trees_refuses(self):
        # Falling back to the full tree for every draw would make the rarefaction vacuous.
        self.run_steps(shared_tree=True)
        with self.assertRaises(SystemExit):
            self.run_steps(start_step=8)

    def test_default_resume_at_step_8_with_per_draw_trees_uses_them(self):
        self.run_steps()
        _, ecosim = self.run_steps(start_step=8)
        self.assertEqual(ecosim.call_args.kwargs["tree_dir"], str(self.draw_trees))


if __name__ == "__main__":
    unittest.main()
