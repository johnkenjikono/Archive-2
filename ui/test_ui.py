import queue
import sys
import tempfile
import time
import unittest
from io import StringIO
from pathlib import Path
from unittest import mock

from ui import (
    PipelineRunner,
    FormValues,
    build_command,
    build_fill_command,
    filled_csv_path,
    copy_local_fastas,
    default_alignment_path,
    find_ecotype_summaries,
    local_input_dir,
    resolve_workdir,
    results_dir,
    species_folder_name,
    validate_form,
)


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


if __name__ == "__main__":
    unittest.main()
