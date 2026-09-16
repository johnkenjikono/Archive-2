from __future__ import annotations

import os
import queue
import shutil
import signal
import subprocess
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

FASTA_SUFFIXES = (".fna", ".fasta", ".fa")

Mode = Literal["species", "csv", "local"]


@dataclass
class FormValues:
    mode: Mode
    species: str = ""
    outgroup: str = ""
    sample_size: int = 200
    seed: int = 42
    csv_path: str = ""
    local_folder: str = ""
    workdir: str = "."
    db: str = ""
    threads: int = 12
    bakta_jobs: int | None = None
    start_step: int = 1
    input_fasta: str = ""
    outgroup_id: str = ""
    setup_only: bool = False


def build_command(values: FormValues, python_exe: str, pipeline_py: Path) -> list[str]:
    argv = [python_exe, str(pipeline_py)]
    if values.mode == "csv":
        argv.extend(["--csv-file", values.csv_path])
        argv.extend(["--sample-size", str(values.sample_size), "--seed", str(values.seed)])
    elif values.mode == "species":
        argv.extend(["--species", values.species])
        if values.outgroup.strip():
            argv.extend(["--outgroup", values.outgroup.strip()])
        argv.extend(["--sample-size", str(values.sample_size), "--seed", str(values.seed)])
    else:
        argv.extend(["--species", values.species])

    argv.extend(
        [
            "--workdir",
            values.workdir,
            "--threads",
            str(values.threads),
            "--start-step",
            str(values.start_step),
        ]
    )
    if values.db.strip():
        argv.extend(["--db", values.db.strip()])
    if values.bakta_jobs is not None:
        argv.extend(["--bakta-jobs", str(values.bakta_jobs)])
    if values.input_fasta.strip():
        argv.extend(["--input-fasta", values.input_fasta.strip()])
    if values.outgroup_id.strip():
        argv.extend(["--outgroup-id", values.outgroup_id.strip()])
    if values.setup_only:
        argv.append("--setup-only")
    return argv


def species_folder_name(species: str) -> str:
    return species.replace(" ", "_")


def local_input_dir(values: FormValues) -> Path:
    return Path(values.workdir) / species_folder_name(values.species) / "input"


def results_dir(values: FormValues) -> Path:
    if values.mode == "csv":
        return Path(values.workdir)
    return Path(values.workdir) / species_folder_name(values.species)


def find_ecotype_summaries(root: Path) -> list[Path]:
    direct = sorted(root.glob("ecosim_output_*/ecotype_summary.csv"))
    nested = sorted(root.glob("*/ecosim_output_*/ecotype_summary.csv"))
    seen: set[Path] = set()
    out: list[Path] = []
    for path in direct + nested:
        resolved = path.resolve()
        if resolved not in seen:
            seen.add(resolved)
            out.append(path)
    return out


def copy_local_fastas(src_folder: Path, dest_input: Path) -> list[Path]:
    if not src_folder.is_dir():
        raise FileNotFoundError(f"Genome folder does not exist: {src_folder}")
    dest_input.mkdir(parents=True, exist_ok=True)
    copied: list[Path] = []
    for path in sorted(src_folder.iterdir()):
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES:
            dest = dest_input / path.name
            shutil.copy2(path, dest)
            copied.append(dest)
    return copied


def default_alignment_path(values: FormValues) -> Path:
    return (
        Path(values.workdir)
        / species_folder_name(values.species)
        / "output_roary"
        / "results"
        / "core_gene_alignment.aln"
    )


def validate_form(values: FormValues) -> list[str]:
    errors: list[str] = []
    if values.mode == "species":
        if not values.species.strip():
            errors.append("Enter a species name.")
    elif values.mode == "csv":
        if not values.csv_path.strip():
            errors.append("Choose a CSV file.")
        elif not Path(values.csv_path).is_file():
            errors.append("CSV file does not exist.")
    elif values.mode == "local":
        if not values.species.strip():
            errors.append("Enter a species name.")
        if not values.local_folder.strip():
            errors.append("Choose a folder of genome files.")
        elif not Path(values.local_folder).is_dir():
            errors.append("Genome folder does not exist.")
        if values.start_step < 2:
            errors.append(
                "Local genome folder runs start at step 2 (Bakta). Choose 2 or later."
            )

    if values.start_step <= 2:
        if not values.db.strip():
            errors.append("Choose a Bakta database folder (required for steps 1–2).")
        elif not Path(values.db).is_dir():
            errors.append("Bakta database folder does not exist.")

    if values.input_fasta.strip():
        if not Path(values.input_fasta).is_file():
            errors.append("Core-alignment FASTA does not exist.")
    elif values.start_step >= 4 and values.mode != "csv":
        if not default_alignment_path(values).is_file():
            errors.append("Choose a core-alignment FASTA, or lower the start step.")
    return errors


class PipelineRunner:
    def __init__(self) -> None:
        self._proc: subprocess.Popen[str] | None = None
        self._thread: threading.Thread | None = None
        self.was_stopped = False

    @property
    def running(self) -> bool:
        return self._proc is not None and self._proc.poll() is None

    def start(self, argv: list[str], cwd: Path, events: queue.Queue) -> None:
        if self.running:
            return
        self.was_stopped = False
        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"
        self._proc = subprocess.Popen(
            argv,
            cwd=str(cwd),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=True,
            env=env,
        )
        self._thread = threading.Thread(
            target=self._pump, args=(events,), daemon=True
        )
        self._thread.start()

    def stop(self) -> None:
        proc = self._proc
        if proc is None or proc.poll() is not None:
            return
        self.was_stopped = True
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        except ProcessLookupError:
            return
        deadline = time.time() + 5
        while time.time() < deadline:
            if proc.poll() is not None:
                return
            time.sleep(0.1)
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except ProcessLookupError:
            return

    def _pump(self, events: queue.Queue) -> None:
        proc = self._proc
        assert proc is not None
        try:
            assert proc.stdout is not None
            for line in proc.stdout:
                events.put(("line", line.rstrip("\n")))
        finally:
            code = proc.wait()
            events.put(("done", code))


import shlex
import sys

REPO_ROOT = Path(__file__).resolve().parent
PIPELINE_PY = REPO_ROOT / "pipeline.py"


def quote_command(argv: list[str]) -> str:
    return shlex.join(argv)


def main() -> None:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk

    class App(tk.Tk):
        def __init__(self) -> None:
            super().__init__()
            self.title("Ecotype discovery pipeline")
            self.minsize(720, 640)
            self.runner = PipelineRunner()
            self.events: queue.Queue = queue.Queue()
            self._form_widgets: list[tk.Widget] = []
            self._pending_values = FormValues(mode="species")
            self._build()
            self.protocol("WM_DELETE_WINDOW", self._on_close)
            self.after(100, self._drain_events)

        def _track(self, widget: tk.Widget) -> tk.Widget:
            self._form_widgets.append(widget)
            return widget

        def _row(self, parent, row: int, label: str, var: tk.Variable, browse: str | None = None):
            ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=4, pady=2)
            entry = self._track(ttk.Entry(parent, textvariable=var, width=48))
            entry.grid(row=row, column=1, sticky="ew", padx=4, pady=2)
            if browse == "file":
                btn = self._track(
                    ttk.Button(parent, text="Browse", command=lambda: var.set(filedialog.askopenfilename() or var.get()))
                )
                btn.grid(row=row, column=2, padx=4)
            elif browse == "dir":
                btn = self._track(
                    ttk.Button(parent, text="Browse", command=lambda: var.set(filedialog.askdirectory() or var.get()))
                )
                btn.grid(row=row, column=2, padx=4)
            parent.columnconfigure(1, weight=1)

        def _build(self) -> None:
            self.mode_var = tk.StringVar(value="species")
            self.species_var = tk.StringVar()
            self.outgroup_var = tk.StringVar()
            self.sample_size_var = tk.StringVar(value="200")
            self.seed_var = tk.StringVar(value="42")
            self.csv_var = tk.StringVar()
            self.local_folder_var = tk.StringVar()
            self.workdir_var = tk.StringVar(value=".")
            self.db_var = tk.StringVar()
            self.threads_var = tk.StringVar(value="12")
            self.bakta_jobs_var = tk.StringVar()
            self.start_step_var = tk.StringVar(value="1")
            self.input_fasta_var = tk.StringVar()
            self.outgroup_id_var = tk.StringVar()
            self.setup_only_var = tk.BooleanVar(value=False)

            mode = ttk.LabelFrame(self, text="Start mode")
            mode.pack(fill="x", padx=8, pady=6)
            for value, text in (
                ("species", "Species name"),
                ("csv", "Batch CSV"),
                ("local", "Local genome folder"),
            ):
                rb = self._track(
                    ttk.Radiobutton(
                        mode,
                        text=text,
                        value=value,
                        variable=self.mode_var,
                        command=self._show_mode,
                    )
                )
                rb.pack(side="left", padx=8, pady=4)

            self.species_frame = ttk.Frame(self)
            self._row(self.species_frame, 0, "Species", self.species_var)
            self._row(self.species_frame, 1, "Outgroup", self.outgroup_var)
            self._row(self.species_frame, 2, "Sample size", self.sample_size_var)
            self._row(self.species_frame, 3, "Random seed", self.seed_var)

            self.csv_frame = ttk.Frame(self)
            self._row(self.csv_frame, 0, "CSV file", self.csv_var, browse="file")
            self._row(self.csv_frame, 1, "Sample size", self.sample_size_var)
            self._row(self.csv_frame, 2, "Random seed", self.seed_var)

            self.local_frame = ttk.Frame(self)
            self._row(self.local_frame, 0, "Species", self.species_var)
            self._row(self.local_frame, 1, "Genome folder", self.local_folder_var, browse="dir")

            options = ttk.LabelFrame(self, text="Options")
            self.options_frame = options
            options.pack(fill="x", padx=8, pady=6)
            self._row(options, 0, "Working directory", self.workdir_var, browse="dir")
            self._row(options, 1, "Bakta database folder", self.db_var, browse="dir")
            self._row(options, 2, "Threads", self.threads_var)
            self._row(options, 3, "Bakta parallel jobs", self.bakta_jobs_var)
            ttk.Label(options, text="Start step").grid(row=4, column=0, sticky="w", padx=4, pady=2)
            step = self._track(
                ttk.Combobox(
                    options,
                    textvariable=self.start_step_var,
                    values=[str(n) for n in range(1, 10)],
                    width=8,
                    state="readonly",
                )
            )
            step.grid(row=4, column=1, sticky="w", padx=4, pady=2)
            self._row(options, 5, "Core-alignment FASTA", self.input_fasta_var, browse="file")
            self._row(options, 6, "Outgroup ID", self.outgroup_id_var)
            setup = self._track(
                ttk.Checkbutton(options, text="Setup only", variable=self.setup_only_var)
            )
            setup.grid(row=7, column=1, sticky="w", padx=4, pady=2)

            run = ttk.LabelFrame(self, text="Run")
            run.pack(fill="both", expand=True, padx=8, pady=6)
            btns = ttk.Frame(run)
            btns.pack(fill="x")
            self._track(ttk.Button(btns, text="Run pipeline", command=self._on_run)).pack(side="left", padx=4, pady=4)
            ttk.Button(btns, text="Stop", command=self._on_stop).pack(side="left", padx=4, pady=4)
            self.log = tk.Text(run, height=16, wrap="word", state="disabled")
            self.log.pack(fill="both", expand=True, padx=4, pady=4)
            self._show_mode()

        def _show_mode(self) -> None:
            self.species_frame.pack_forget()
            self.csv_frame.pack_forget()
            self.local_frame.pack_forget()
            mode = self.mode_var.get()
            if mode == "species":
                self.species_frame.pack(fill="x", padx=8, before=self.options_frame)
            elif mode == "csv":
                self.csv_frame.pack(fill="x", padx=8, before=self.options_frame)
            else:
                self.local_frame.pack(fill="x", padx=8, before=self.options_frame)
                if int(self.start_step_var.get() or "1") < 2:
                    self.start_step_var.set("2")

        def _collect(self) -> FormValues:
            jobs_raw = self.bakta_jobs_var.get().strip()
            bakta_jobs = int(jobs_raw) if jobs_raw else None
            return FormValues(
                mode=self.mode_var.get(),  # type: ignore[arg-type]
                species=self.species_var.get().strip(),
                outgroup=self.outgroup_var.get().strip(),
                sample_size=int(self.sample_size_var.get() or "200"),
                seed=int(self.seed_var.get() or "42"),
                csv_path=self.csv_var.get().strip(),
                local_folder=self.local_folder_var.get().strip(),
                workdir=self.workdir_var.get().strip() or ".",
                db=self.db_var.get().strip(),
                threads=int(self.threads_var.get() or "12"),
                bakta_jobs=bakta_jobs,
                start_step=int(self.start_step_var.get() or "1"),
                input_fasta=self.input_fasta_var.get().strip(),
                outgroup_id=self.outgroup_id_var.get().strip(),
                setup_only=bool(self.setup_only_var.get()),
            )

        def _set_form_enabled(self, enabled: bool) -> None:
            state = "normal" if enabled else "disabled"
            combo_state = "readonly" if enabled else "disabled"
            for widget in self._form_widgets:
                try:
                    if isinstance(widget, ttk.Combobox):
                        widget.configure(state=combo_state)
                    else:
                        widget.configure(state=state)
                except tk.TclError:
                    pass

        def _log(self, line: str) -> None:
            try:
                self.log.configure(state="normal")
                self.log.insert("end", line + "\n")
                self.log.see("end")
                self.log.configure(state="disabled")
            except tk.TclError:
                return

        def _on_run(self) -> None:
            if self.runner.running:
                return
            if not PIPELINE_PY.is_file():
                messagebox.showerror("Missing pipeline", f"Could not find {PIPELINE_PY}")
                return
            try:
                values = self._collect()
            except ValueError:
                messagebox.showerror("Cannot run", "Threads, sample size, seed, start step, and Bakta jobs must be integers.")
                return
            errors = validate_form(values)
            if errors:
                for error in errors:
                    self._log(error)
                messagebox.showerror("Cannot run", "\n".join(errors))
                return
            if values.mode == "local" and values.start_step <= 2:
                try:
                    copied = copy_local_fastas(
                        Path(values.local_folder), local_input_dir(values)
                    )
                except OSError as exc:
                    self._log(str(exc))
                    messagebox.showerror("Cannot copy genomes", str(exc))
                    return
                if not copied:
                    msg = "No .fna, .fasta, or .fa files in that folder."
                    self._log(msg)
                    messagebox.showerror("Cannot copy genomes", msg)
                    return
                self._log(f"Copied {len(copied)} genome file(s) to {local_input_dir(values)}")
            argv = build_command(values, sys.executable, PIPELINE_PY)
            self._pending_values = values
            self._log(quote_command(argv))
            self._set_form_enabled(False)
            self.runner.start(argv, REPO_ROOT, self.events)

        def _on_stop(self) -> None:
            if self.runner.running:
                self._log("Stopping the run…")
                self.runner.stop()

        def _drain_events(self) -> None:
            try:
                while True:
                    kind, payload = self.events.get_nowait()
                    if kind == "line":
                        self._log(str(payload))
                    elif kind == "done":
                        self._on_done(int(payload))
            except queue.Empty:
                pass
            except tk.TclError:
                pass
            self.after(100, self._drain_events)

        def _on_done(self, code: int) -> None:
            self._set_form_enabled(True)
            if self.runner.was_stopped:
                self._log("Run stopped. Partial output was left on disk.")
                return
            if code == 0:
                out = results_dir(self._pending_values)
                self._log(f"Pipeline finished. Results: {out}")
                for summary in find_ecotype_summaries(out):
                    self._log(f"Ecotype summary: {summary}")
                return
            self._log("Pipeline failed — see the messages above.")

        def _on_close(self) -> None:
            if self.runner.running:
                if not messagebox.askyesno(
                    "Stop and quit",
                    "A run is still going. Stop it and quit?",
                ):
                    return
                self.runner.stop()
            self.destroy()

    App().mainloop()


if __name__ == "__main__":
    main()
