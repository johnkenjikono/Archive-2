"""Run from the repo root: python -m unittest steps.test_med"""
import csv
import tempfile
import unittest
from pathlib import Path

from steps.med import compare

CLADES = {"A": "AAAAAAAAAA", "B": "CCCCCCCCCC", "C": "GGGGGGGGGG"}


class MedTest(unittest.TestCase):
    def test_recovers_clades_despite_singletons_and_gaps(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        aln, members = Path(tmp.name) / "core.aln", Path(tmp.name) / "ecotype_membership.csv"

        records, rows = [">outgroup\n" + "T" * 30 + "\n"], []
        for c, (clade, block) in enumerate(CLADES.items()):
            for i in range(5):
                seq = list(block + "ACGTACGTAC" * 2)
                seq[10 + 5 * c + i] = "N" if i == 4 else "T"  # a different column per genome
                if i == 0:
                    seq[0:5] = "-----"  # a missing gene
                records.append(f">{clade}{i}\n{''.join(seq)}\n")
                rows.append(["full_core_genome_results.xml", "full", "", clade, 5, f"{clade}{i}", f"{clade}{i}"])
        aln.write_text("".join(records))
        with open(members, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["file", "gene_count", "trial", "ecotype", "ecotype_size", "taxon", "representative"])
            writer.writerows(rows)

        nodes, ecotypes, ari = compare(aln, tmp.name)

        self.assertNotIn("outgroup", nodes)
        self.assertEqual(len(set(nodes.values())), 3)
        self.assertEqual(ari, 1.0)
        self.assertTrue((Path(tmp.name) / "med_nodes.csv").exists())


if __name__ == "__main__":
    unittest.main()
