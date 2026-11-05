"""Run from the repo root: python -m unittest steps.test_interactive_tree"""
import json
import re
import tempfile
import unittest
from pathlib import Path

from steps.interactive_tree import write_interactive_tree
from steps.parsing import FULL_CORE_XML

# outgroup is sister to everything; E1 = (a,b) and E2 = c are clades, E3 = d + e is split by c2 (E4).
TREE = "(outgroup:0.5,((a:0.1,b:0.2):0.1,(c:0.3,(d:0.1,(c2:0.1,e:0.2):0.1):0.1):0.1):0.1);"
ECOTYPES = {1: ["a", "b"], 2: ["c"], 3: ["d", "e"], 4: ["c2"]}


def write_run(folder, support=True):
    ecotypes = "".join(
        f'<ecotype number="{n}">' + "".join(f'<member name="{m}"/>' for m in members) + "</ecotype>"
        for n, members in ECOTYPES.items())
    (folder / FULL_CORE_XML).write_text(
        f'<root><phylogeny><outgroup value="outgroup"/><tree value="{TREE}"/></phylogeny>'
        f"<demarcation><ecotypes>{ecotypes}</ecotypes></demarcation></root>")
    (folder / "clone_groups.csv").write_text(
        "representative,member,group_size,differences,divergence\n"
        "a,a,3,0,0\na,a2,3,0,0\na,a3,3,0,0\nb,b,1,0,0\n")
    (folder / "genome_metadata.csv").write_text("accession,host,niche,country,year\na2,cow,soil,Peru,2001\n")
    if support:
        (folder / "ecotype_support.csv").write_text(
            "ecotype,size,support_g1,support_g20,members\n"
            "1,2,0.0000,0.2500,a;b\n2,1,0.5000,1.0000,c\n3,2,0.0000,0.0000,d;e\n4,1,0.0000,0.0000,c2\n")


def tree_data(folder):
    page = write_interactive_tree(folder).read_text()
    return json.loads(re.search(r"const DATA = (.*);\n", page).group(1))


class InteractiveTreeDataTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.folder = Path(tmp.name) / "Listeria_monocytogenes" / "ecosim_output_x"
        self.folder.mkdir(parents=True)

    def row(self, data, eco):
        return next(r for r in data["rows"] if r["eco"] == eco)

    def test_outgroup_branch_is_not_drawn(self):
        write_run(self.folder)
        data = tree_data(self.folder)
        self.assertEqual([t["name"] for t in data["tips"]], ["a", "b", "c", "d", "c2", "e"])
        self.assertEqual(len(data["nodes"]), 11)
        self.assertIn("outgroup branch not drawn", data["subtitle"])
        self.assertEqual(data["nodes"][0]["x"], 0)  # the tree starts at the left edge, not one pruned branch in

    def test_subtitle_names_the_outgroup_only_when_the_run_recorded_it(self):
        write_run(self.folder)
        self.assertNotIn("outgroup outgroup", tree_data(self.folder)["subtitle"])
        (self.folder.parent / "run_parameters.json").write_text('{"outgroup": "Listeria innocua"}')
        self.assertIn("Rooted with outgroup Listeria innocua;", tree_data(self.folder)["subtitle"])

    def test_overview_row_counts_unique_sequences_and_expanded_genomes(self):
        write_run(self.folder)
        row = self.row(tree_data(self.folder), 1)
        self.assertEqual((row["unique"], row["genomes"]), (2, 4))
        self.assertEqual((row["first"], row["last"]), (0, 1))

    def test_each_clade_collapses_to_one_crown_with_its_tips_hidden(self):
        write_run(self.folder)
        data = tree_data(self.folder)
        crown = next(i for i, n in enumerate(data["nodes"]) if n["r"] == data["rows"].index(self.row(data, 1)))
        below = [n for n in data["nodes"] if n["p"] == crown]
        self.assertEqual(len(below), 2)
        self.assertTrue(all(n["h"] for n in below))
        self.assertFalse(data["nodes"][crown]["h"])

    def test_non_monophyletic_ecotype_gets_a_row_per_clade(self):
        write_run(self.folder)
        rows = [r for r in tree_data(self.folder)["rows"] if r["eco"] == 3]
        self.assertEqual([(r["first"], r["last"]) for r in rows], [(3, 3), (5, 5)])
        self.assertTrue(all(r["split"] for r in rows))
        self.assertFalse(self.row(tree_data(self.folder), 1)["split"])

    def test_rarefaction_recovery_follows_gene_count_order(self):
        write_run(self.folder)
        data = tree_data(self.folder)
        self.assertEqual(data["gene_counts"], [1, 20])
        self.assertEqual(self.row(data, 2)["support"], [0.5, 1.0])

    def test_no_support_table_means_no_recovery_columns(self):
        write_run(self.folder, support=False)
        data = tree_data(self.folder)
        self.assertEqual(data["gene_counts"], [])
        self.assertIsNone(self.row(data, 1)["support"])

    def test_stale_support_rows_with_other_members_are_ignored(self):
        write_run(self.folder)
        (self.folder / "ecotype_support.csv").write_text(
            "ecotype,size,support_g1,support_g20,members\n1,2,0.9,0.9,x;y\n2,1,0.5,1.0,c\n")
        data = tree_data(self.folder)
        self.assertIsNone(self.row(data, 1)["support"])
        self.assertEqual(self.row(data, 2)["support"], [0.5, 1.0])

    def test_no_matching_support_rows_means_no_recovery_columns(self):
        write_run(self.folder)
        (self.folder / "ecotype_support.csv").write_text("ecotype,size,support_g1,members\n1,2,0.9,x;y\n")
        data = tree_data(self.folder)
        self.assertEqual(data["gene_counts"], [])
        self.assertIsNone(self.row(data, 1)["support"])

    def test_support_is_only_used_with_the_full_core_tree(self):
        write_run(self.folder)
        (self.folder / FULL_CORE_XML).rename(self.folder / "core_g20_t1_results.xml")
        data = tree_data(self.folder)
        self.assertEqual(data["gene_counts"], [])
        self.assertIsNone(self.row(data, 1)["support"])

    def test_strain_keeps_clones_and_environmental_metadata(self):
        write_run(self.folder)
        tip = tree_data(self.folder)["tips"][0]
        self.assertEqual(tip["ecotype"], 1)
        clone = next(g for g in tip["genomes"] if g["accession"] == "a2")
        self.assertEqual((clone["host"], clone["niche"], clone["country"], clone["year"]),
                         ("cow", "soil", "Peru", "2001"))

    def test_title_names_species_genomes_and_ecotypes(self):
        write_run(self.folder)
        title = tree_data(self.folder)["title"]
        self.assertEqual(title, "Listeria monocytogenes · 8 genomes · 4 ecotypes (full core genome)")


if __name__ == "__main__":
    unittest.main()
