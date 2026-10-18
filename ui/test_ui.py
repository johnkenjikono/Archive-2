import queue
import sys
import tempfile
import time
import unittest
from io import StringIO
from pathlib import Path
from unittest import mock

from ui import (
    DEFAULT_KEEP,
    PipelineRunner,
    FormValues,
    build_command,
    build_fill_command,
    filled_csv_path,
    copy_local_fastas,
    default_alignment_path,
    find_bakta_dbs,
    find_ecotype_summaries,
    local_input_dir,
    resolve_workdir,
    results_dir,
    species_folder_name,
    validate_form,
    format_elapsed,
    is_progress_tick,
    load_settings,
    parse_number,
    progress_total,
    save_settings,
    line_kind,
    parse_step_label,
    step_from_line,
    step_label,
)


class StepStatusTests(unittest.TestCase):
    def test_step_from_pipeline_headers(self):
        self.assertEqual(step_from_line("--- Downloading genomes for Bacillus subtilis ---"), 1)
        self.assertEqual(step_from_line("Running Bakta on GCF_1.fna..."), 2)
        self.assertEqual(step_from_line("Running Panaroo on 12 files..."), 3)
        self.assertEqual(step_from_line("--- 4b. Collapsing clones (per-site divergence <= 1e-05) ---"), 4)
        self.assertEqual(step_from_line("--- 10. Plotting ecotypes (trees: yes, pie charts: yes) ---"), 10)
        self.assertIsNone(step_from_line("--- Skipping Step 5: Using existing unrooted tree ---"))
        self.assertIsNone(step_from_line("Removed: EcoSim XML files"))

    def test_line_kind(self):
        self.assertEqual(line_kind("--- 8. Running EcoSim ---"), "step")
        self.assertEqual(line_kind("Error: Could not download outgroup genome"), "error")
        self.assertEqual(line_kind("Warning: Bakta failed on 1/10 genome(s)."), "error")
        self.assertEqual(line_kind("Warning: could not remove x"), "warn")
        self.assertIsNone(line_kind("Saved ecotype_summary.csv"))

    def test_step_label_round_trip(self):
        self.assertEqual(parse_step_label(step_label(7)), 7)

    def test_format_elapsed(self):
        self.assertEqual(format_elapsed(3725.9), "1:02:05")


class ProgressTests(unittest.TestCase):
    def test_totals_from_bakta_and_ecosim(self):
        self.assertEqual(progress_total("Found 12 FASTA files to process (3 parallel Bakta job(s), 4 thread(s) each)."), 12)
        self.assertEqual(progress_total("📁 Found 101 FASTA file(s) in /x/rarefaction"), 101)
        self.assertIsNone(progress_total("Found 300 FASTA files; limiting to first 201 for Bakta processing."))

    def test_ticks(self):
        self.assertTrue(is_progress_tick("Successfully processed GCF_1.fna"))
        self.assertTrue(is_progress_tick("Skipping GCF_1.fna: Bakta output already exists at x"))
        self.assertTrue(is_progress_tick("✅ EcoSim output saved: x.xml"))
        self.assertTrue(is_progress_tick("❌ EcoSim error for a.fasta: boom"))
        self.assertFalse(is_progress_tick("Running Bakta on GCF_1.fna..."))


class ParseNumberTests(unittest.TestCase):
    def test_blank_uses_default(self):
        self.assertEqual(parse_number("  ", "Threads", "12"), 12)
        self.assertEqual(parse_number("", "Clone threshold", "1e-5", float), 1e-5)

    def test_error_names_field_and_value(self):
        with self.assertRaises(ValueError) as ctx:
            parse_number("1o", "Threads", "12")
        self.assertEqual(str(ctx.exception), 'Threads must be a whole number (got "1o").')
        with self.assertRaises(ValueError) as ctx:
            parse_number("abc", "Clone threshold", "1e-5", float)
        self.assertEqual(str(ctx.exception), 'Clone threshold must be a number (got "abc").')


class SettingsTests(unittest.TestCase):
    def test_round_trip_and_bad_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "s.json"
            self.assertEqual(load_settings(path), {})
            save_settings(path, {"workdir": "/w", "db": "/db"})
            self.assertEqual(load_settings(path), {"workdir": "/w", "db": "/db"})
            path.write_text("not json")
            self.assertEqual(load_settings(path), {})
            path.write_text('["a list"]')
            self.assertEqual(load_settings(path), {})
            save_settings(Path(tmp) / "missing" / "s.json", {"db": "x"})  # must not raise


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

    def test_alignment_ok_when_saved_copy_exists(self):
        with tempfile.TemporaryDirectory() as tmp:
            aln = Path(tmp) / "Bacillus_subtilis" / "pangenome" / "core_gene_alignment.aln"
            aln.parent.mkdir(parents=True)
            aln.write_text(">a\nACGT\n")
            values = FormValues(mode="species", species="Bacillus subtilis", start_step=4, workdir=tmp)
            self.assertEqual(validate_form(values), [])

    def test_step_8_needs_no_alignment(self):
        values = FormValues(mode="species", species="Bacillus subtilis", start_step=8, workdir="/no/such")
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
        self.assertIn("--workdir", argv)
        self.assertEqual(argv[argv.index("--workdir") + 1], "/tmp/work")
        self.assertIn("--threads", argv)
        self.assertEqual(argv[argv.index("--threads") + 1], "8")
        self.assertIn("--start-step", argv)
        self.assertEqual(argv[argv.index("--start-step") + 1], "1")

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
        self.assertNotIn("--keep", argv)
        self.assertNotIn("--discard", argv)
        self.assertNotIn("--no-full-core", argv)
        self.assertNotIn("--each-gene", argv)
        self.assertNotIn("--replicates", argv)

    def test_keep_and_discard_are_differences_from_defaults(self):
        keep = (DEFAULT_KEEP - {"pie_charts", "tree_plots"}) | {"core_alignment"}
        argv = build_command(
            FormValues(mode="species", species="X", start_step=3, keep=frozenset(keep)),
            python_exe="python",
            pipeline_py=Path("pipeline.py"),
        )
        self.assertEqual(argv[argv.index("--discard") + 1], "pie_charts,tree_plots")
        self.assertEqual(argv[argv.index("--keep") + 1], "core_alignment")

    def test_rarefaction_options(self):
        argv = build_command(
            FormValues(mode="species", species="X", start_step=3,
                       full_core=False, each_gene=True, replicates=100),
            python_exe="python",
            pipeline_py=Path("pipeline.py"),
        )
        self.assertIn("--no-full-core", argv)
        self.assertIn("--each-gene", argv)
        self.assertEqual(argv[argv.index("--replicates") + 1], "100")

    def test_clone_threshold_only_when_changed(self):
        default = build_command(FormValues(mode="species", species="X"), "python", Path("pipeline.py"))
        self.assertNotIn("--clone-threshold", default)
        argv = build_command(FormValues(mode="species", species="X", clone_threshold=0.0), "python", Path("pipeline.py"))
        self.assertEqual(argv[argv.index("--clone-threshold") + 1], "0")


class FillOutgroupsCommandTests(unittest.TestCase):
    def test_llm_flag_only_when_checked(self):
        off = build_fill_command("s.csv", False, "python", Path("fill_outgroups.py"))
        on = build_fill_command("s.csv", True, "python", Path("fill_outgroups.py"))
        self.assertEqual(off, ["python", "fill_outgroups.py", "s.csv"])
        self.assertEqual(on, ["python", "fill_outgroups.py", "s.csv", "--llm"])

    def test_filled_csv_path_matches_script_default(self):
        self.assertEqual(filled_csv_path("/data/species.csv"), Path("/data/species_outgroups.csv"))


class LocalCopyTests(unittest.TestCase):
    def test_copies_fasta_extensions_non_recursive(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "src"
            nested = src / "nested"
            dest = Path(tmp) / "dest"
            src.mkdir()
            nested.mkdir()
            (src / "a.fna").write_text(">a\nA\n")
            (src / "b.fasta").write_text(">b\nA\n")
            (src / "c.fa").write_text(">c\nA\n")
            (src / "skip.txt").write_text("nope")
            (nested / "hidden.fna").write_text(">h\nA\n")
            copied = copy_local_fastas(src, dest)
            names = sorted(path.name for path in copied)
            self.assertEqual(names, ["a.fna", "b.fasta", "c.fa"])
            self.assertTrue((src / "a.fna").exists())
            self.assertFalse((dest / "hidden.fna").exists())

    def test_overwrite_same_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "src"
            dest = Path(tmp) / "dest"
            src.mkdir()
            dest.mkdir()
            (src / "a.fna").write_text("new")
            (dest / "a.fna").write_text("old")
            copy_local_fastas(src, dest)
            self.assertEqual((dest / "a.fna").read_text(), "new")

    def test_missing_source_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                copy_local_fastas(Path(tmp) / "missing", Path(tmp) / "dest")


class ResultPathTests(unittest.TestCase):
    def test_local_input_dir_uses_species_folder(self):
        values = FormValues(
            mode="local",
            species="Bacillus subtilis",
            workdir="/work",
        )
        self.assertEqual(
            local_input_dir(values),
            Path("/work/Bacillus_subtilis/input"),
        )

    def test_results_dir_csv_is_workdir(self):
        self.assertEqual(
            results_dir(FormValues(mode="csv", workdir="/out")),
            Path("/out"),
        )

    def test_results_dir_species_uses_folder_name(self):
        self.assertEqual(
            results_dir(FormValues(mode="species", species="Bacillus subtilis", workdir="/out")),
            Path("/out/Bacillus_subtilis"),
        )

    def test_find_ecotype_summaries_in_species_and_batch_layouts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            one = root / "ecosim_output_x" / "ecotype_summary.csv"
            two = root / "Bacillus_subtilis" / "ecosim_output_y" / "ecotype_summary.csv"
            one.parent.mkdir(parents=True)
            two.parent.mkdir(parents=True)
            one.write_text("file,ecotype_count\n")
            two.write_text("file,ecotype_count\n")
            found = {path.resolve() for path in find_ecotype_summaries(root)}
            self.assertEqual(found, {one.resolve(), two.resolve()})


class ResolveWorkdirTests(unittest.TestCase):
    def test_relative_path_resolves_against_cwd(self):
        resolved = resolve_workdir(".")
        self.assertTrue(Path(resolved).is_absolute())
        self.assertEqual(Path(resolved), Path(".").expanduser().resolve())

    def test_expands_user_home(self):
        resolved = resolve_workdir("~")
        self.assertEqual(Path(resolved), Path.home().resolve())


class PipelineRunnerTests(unittest.TestCase):
    def test_start_sets_python_unbuffered_env(self):
        proc = mock.MagicMock()
        proc.poll.return_value = None
        proc.stdout = StringIO("")
        proc.wait.return_value = 0
        events: queue.Queue = queue.Queue()
        with mock.patch("ui.subprocess.Popen", return_value=proc) as popen:
            runner = PipelineRunner()
            runner.start([sys.executable, "-c", "pass"], Path("."), events)
            env = popen.call_args.kwargs["env"]
        self.assertEqual(env["PYTHONUNBUFFERED"], "1")
        self.assertIn("PATH", env)
        deadline = time.time() + 2
        while time.time() < deadline:
            try:
                kind, _payload = events.get(timeout=0.1)
            except queue.Empty:
                continue
            if kind == "done":
                break

    def test_streams_output_and_exit_code(self):
        events: queue.Queue = queue.Queue()
        runner = PipelineRunner()
        runner.start(
            [sys.executable, "-c", "print('hello-ui'); raise SystemExit(0)"],
            cwd=Path("."),
            events=events,
        )
        lines = []
        code = None
        deadline = time.time() + 10
        while time.time() < deadline:
            try:
                kind, payload = events.get(timeout=0.2)
            except queue.Empty:
                continue
            if kind == "line":
                lines.append(payload)
            elif kind == "done":
                code = payload
                break
        self.assertFalse(runner.running)
        self.assertEqual(code, 0)
        self.assertTrue(any("hello-ui" in line for line in lines))
        self.assertTrue(runner._proc.stdout.closed)

    def test_stop_kills_long_process(self):
        events: queue.Queue = queue.Queue()
        runner = PipelineRunner()
        runner.start(
            [sys.executable, "-c", "import time; time.sleep(30)"],
            cwd=Path("."),
            events=events,
        )
        time.sleep(0.3)
        self.assertTrue(runner.running)
        runner.stop()
        deadline = time.time() + 10
        code = None
        while time.time() < deadline:
            try:
                kind, payload = events.get(timeout=0.2)
            except queue.Empty:
                continue
            if kind == "done":
                code = payload
                break
        self.assertFalse(runner.running)
        self.assertTrue(runner.was_stopped)
        self.assertIsNotNone(code)
        self.assertNotEqual(code, 0)

    def test_second_start_ignored_while_running(self):
        events: queue.Queue = queue.Queue()
        runner = PipelineRunner()
        runner.start(
            [sys.executable, "-c", "import time; time.sleep(5)"],
            cwd=Path("."),
            events=events,
        )
        time.sleep(0.2)
        first_proc = runner._proc
        runner.start(
            [sys.executable, "-c", "print('should-not-run')"],
            cwd=Path("."),
            events=events,
        )
        self.assertIs(runner._proc, first_proc)
        runner.stop()
        deadline = time.time() + 10
        while time.time() < deadline:
            try:
                kind, _payload = events.get(timeout=0.2)
            except queue.Empty:
                continue
            if kind == "done":
                break


class FindBaktaDbsTests(unittest.TestCase):
    def test_finds_db_and_db_light_only_with_version_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "bakta_db"
            for name in ("db", "db-light", "empty"):
                (root / name).mkdir(parents=True)
            (root / "db" / "version.json").write_text("{}")
            (root / "db-light" / "version.json").write_text("{}")
            found = find_bakta_dbs([root, root])
            self.assertEqual(found, [str((root / "db").resolve()), str((root / "db-light").resolve())])

    def test_missing_root_finds_nothing(self):
        self.assertEqual(find_bakta_dbs([Path("/no/such/bakta_db")]), [])


if __name__ == "__main__":
    unittest.main()
