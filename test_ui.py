import tempfile
import unittest
from pathlib import Path

from ui import FormValues, build_command, default_alignment_path, species_folder_name, validate_form


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


class BuildCommandTests(unittest.TestCase):
    def test_species_command_includes_species_and_outgroup(self):
        argv = build_command(
            FormValues(
                mode="species",
                species="Bacillus subtilis",
                outgroup="Bacillus licheniformis",
                sample_size=200,
                seed=42,
                workdir="/tmp/work",
                db="/tmp/db",
                threads=8,
                start_step=1,
            ),
            python_exe="/usr/bin/python",
            pipeline_py=Path("/repo/pipeline.py"),
        )
        self.assertEqual(argv[0], "/usr/bin/python")
        self.assertEqual(argv[1], "/repo/pipeline.py")
        self.assertIn("--species", argv)
        self.assertEqual(argv[argv.index("--species") + 1], "Bacillus subtilis")
        self.assertIn("--outgroup", argv)
        self.assertNotIn("--csv-file", argv)
        self.assertIn("--sample-size", argv)
        self.assertIn("--seed", argv)

    def test_csv_command_omits_species(self):
        argv = build_command(
            FormValues(
                mode="csv",
                csv_path="/data/batch.csv",
                sample_size=50,
                seed=7,
                workdir=".",
                threads=12,
                start_step=1,
            ),
            python_exe="python",
            pipeline_py=Path("pipeline.py"),
        )
        self.assertIn("--csv-file", argv)
        self.assertEqual(argv[argv.index("--csv-file") + 1], "/data/batch.csv")
        self.assertNotIn("--species", argv)
        self.assertNotIn("--outgroup", argv)

    def test_local_command_skips_sample_flags(self):
        argv = build_command(
            FormValues(
                mode="local",
                species="Treponema paraluiscuniculi",
                local_folder="/genomes",
                start_step=2,
                workdir=".",
                threads=12,
            ),
            python_exe="python",
            pipeline_py=Path("pipeline.py"),
        )
        self.assertIn("--species", argv)
        self.assertGreaterEqual(int(argv[argv.index("--start-step") + 1]), 2)
        self.assertNotIn("--sample-size", argv)
        self.assertNotIn("--seed", argv)
        self.assertNotIn("--csv-file", argv)

    def test_setup_only_flag(self):
        argv = build_command(
            FormValues(mode="species", species="X", start_step=3, setup_only=True),
            python_exe="python",
            pipeline_py=Path("pipeline.py"),
        )
        self.assertIn("--setup-only", argv)

    def test_empty_optional_flags_omitted(self):
        argv = build_command(
            FormValues(
                mode="species",
                species="X",
                start_step=3,
                bakta_jobs=None,
                input_fasta="",
                outgroup_id="",
                db="",
            ),
            python_exe="python",
            pipeline_py=Path("pipeline.py"),
        )
        self.assertNotIn("--bakta-jobs", argv)
        self.assertNotIn("--input-fasta", argv)
        self.assertNotIn("--outgroup-id", argv)
        self.assertNotIn("--db", argv)


if __name__ == "__main__":
    unittest.main()
