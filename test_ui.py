import tempfile
import unittest
from pathlib import Path

from ui import FormValues, default_alignment_path, species_folder_name, validate_form


class SpeciesFolderTests(unittest.TestCase):
    def test_spaces_become_underscores(self):
        self.assertEqual(species_folder_name("Bacillus subtilis"), "Bacillus_subtilis")


class ValidateFormTests(unittest.TestCase):
    def test_species_mode_requires_name(self):
        errors = validate_form(FormValues(mode="species", start_step=3, db=""))
        self.assertIn("Enter a species name.", errors)

    def test_csv_mode_requires_existing_file(self):
        errors = validate_form(FormValues(mode="csv", csv_path="/no/such.csv", start_step=3))
        self.assertIn("CSV file does not exist.", errors)

    def test_local_mode_rejects_step_one(self):
        with tempfile.TemporaryDirectory() as tmp:
            errors = validate_form(
                FormValues(
                    mode="local",
                    species="Bacillus subtilis",
                    local_folder=tmp,
                    start_step=1,
                    db=tmp,
                )
            )
        self.assertIn(
            "Local genome folder runs start at step 2 (Bakta). Choose 2 or later.",
            errors,
        )

    def test_db_required_when_start_step_at_most_2(self):
        errors = validate_form(
            FormValues(mode="species", species="Bacillus subtilis", start_step=2, db="")
        )
        self.assertIn("Choose a Bakta database folder (required for steps 1–2).", errors)

    def test_db_must_exist_when_required(self):
        errors = validate_form(
            FormValues(
                mode="species",
                species="Bacillus subtilis",
                start_step=1,
                db="/no/such/db",
            )
        )
        self.assertIn("Bakta database folder does not exist.", errors)

    def test_alignment_required_for_step_4_without_default_file(self):
        errors = validate_form(
            FormValues(
                mode="species",
                species="Bacillus subtilis",
                start_step=4,
                workdir=".",
                input_fasta="",
            )
        )
        self.assertIn(
            "Choose a core-alignment FASTA, or lower the start step.",
            errors,
        )

    def test_alignment_ok_when_default_path_exists(self):
        with tempfile.TemporaryDirectory() as tmp:
            aln = (
                Path(tmp)
                / "Bacillus_subtilis"
                / "output_roary"
                / "results"
                / "core_gene_alignment.aln"
            )
            aln.parent.mkdir(parents=True)
            aln.write_text(">a\nACGT\n")
            values = FormValues(
                mode="species",
                species="Bacillus subtilis",
                start_step=4,
                workdir=tmp,
                input_fasta="",
            )
            self.assertEqual(default_alignment_path(values), aln)
            self.assertEqual(validate_form(values), [])

    def test_filled_alignment_must_exist(self):
        errors = validate_form(
            FormValues(
                mode="species",
                species="Bacillus subtilis",
                start_step=4,
                input_fasta="/no/such.aln",
            )
        )
        self.assertIn("Core-alignment FASTA does not exist.", errors)

    def test_csv_step_4_without_alignment_is_allowed(self):
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as handle:
            csv_path = handle.name
        try:
            errors = validate_form(
                FormValues(mode="csv", csv_path=csv_path, start_step=4, input_fasta="")
            )
        finally:
            Path(csv_path).unlink()
        self.assertEqual(errors, [])


if __name__ == "__main__":
    unittest.main()
