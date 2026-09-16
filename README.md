# Ecotype Discovery Pipeline

Give it a bacterial species name. It downloads complete genomes from NCBI, builds a core-genome phylogeny, and runs **Ecotype Simulation (EcoSim)** to estimate how many ecotypes the species contains, and how stable that estimate is as more genes are sampled.

```
NCBI Datasets ─► Bakta ─► Panaroo ─► VeryFastTree ─► reroot ─► rarefaction ─► EcoSim ─► ecotype counts
  (download)    (annotate) (core genes)  (tree)     (outgroup)  (gene subsets)   (demarcate)   (CSV)
```

Run it from the desktop window (`bash run_ui.sh`) or from the terminal, one species at a time or a whole CSV of species.

---

## Contents

- [Quick start](#quick-start)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Graphical interface](#graphical-interface)
  - [One species](#one-species)
  - [Many species from a CSV](#many-species-from-a-csv)
- [How the pipeline works](#how-the-pipeline-works)
- [Outputs](#outputs)
- [Performance tuning](#performance-tuning)
- [Running steps on their own](#running-steps-on-their-own)
- [Repository layout](#repository-layout)
- [Troubleshooting](#troubleshooting)
- [Citations](#citations)

---

## Quick start

```bash
bash setup.sh                                   # conda envs + Python venv + tool checks
source venv/bin/activate                        # fish: source venv/bin/activate.fish

conda run -n bakta_env bakta_db download --output ~/bakta_db   # one-time database download

bash run_ui.sh                                  # desktop window; no CLI flags
```

In the window, choose **Species name**, enter the species, Browse to the Bakta `db` folder, set a small **Sample size** (try 5 first), and press **Run pipeline**.

The same run from the terminal:

```bash
python pipeline.py \
  --species "Treponema paraluiscuniculi" \
  --sample-size 5 \
  --db ~/bakta_db/db \
  --threads 8
```

The ecotype counts land in `./Treponema_paraluiscuniculi/ecosim_output_core_gene_alignment/ecotype_summary.csv`.

> A full run of 200 genomes takes hours, mostly in Bakta and EcoSim.

---

## Requirements

| Tool | Used in | Installed by `setup.sh`? |
|---|---|---|
| [Conda / Miniforge](https://github.com/conda-forge/miniforge) | Bakta and Panaroo environments | No, install it first |
| Python ≥ 3.10 | The pipeline and the desktop window (stdlib `tkinter`) | Creates `venv/` from `requirements.txt` |
| [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) (`datasets`) | Step 1, outgroup/type-strain lookup | Tries `conda install ncbi-datasets-cli` |
| [Bakta](https://github.com/oschwengers/bakta) + its database | Step 2 | Creates `bakta_env`. You download the database |
| [Panaroo](https://github.com/gtonkinhill/panaroo) + MAFFT | Step 3 | Creates `panaroo_env` |
| [VeryFastTree](https://github.com/citiususc/veryfasttree) (or FastTree) | Step 5 | Tries `brew` (macOS) or `apt` (Linux) |
| Java 8+ | Step 8 (EcoSim) | No, only checked |

**EcoSim is bundled.** `ecosim.jar` and its helper binaries in `bin/` ship with the repo. The binaries in `bin/` are **macOS arm64 (Apple Silicon)** builds. On any other platform, supply EcoSim binaries built for that platform (see [Configuration](#configuration)).

---

## Installation

### 1. Environments and dependencies

```bash
bash setup.sh
```

This script:
- creates the `bakta_env` and `panaroo_env` conda environments (skips any that already exist)
- creates `venv/` and installs `biopython`, `pandas` and `openpyxl`
- checks for `datasets`, `java` and `veryfasttree`, and installs what it can

### 2. Bakta database (one-time)

```bash
conda activate bakta_env
bakta_db download --output /path/to/bakta_db              # full database (recommended)
# bakta_db download --output /path/to/bakta_db --type light  # much smaller, less precise
conda deactivate
```

Pass the **`db` (or `db-light`) subfolder** to the pipeline, e.g. `--db /path/to/bakta_db/db`.

### 3. Activate the venv before each session

```bash
source venv/bin/activate        # bash / zsh
source venv/bin/activate.fish   # fish
```

`bash run_ui.sh` does this for you when `venv/` exists.

### Configuration

No configuration is needed if you use the bundled EcoSim. To use your own copy, set:

| Variable | Default | Purpose |
|---|---|---|
| `ECOSIM_JAR` | `./ecosim.jar` | Path to the EcoSim jar |
| `ECOSIM_DIR` | repo root | Working directory for EcoSim (must contain its `bin/` helpers) |

---

## Usage

### Graphical interface

The window collects the same options as `pipeline.py` and runs it for you. You do not need to assemble a command.

```bash
bash run_ui.sh
```

That script activates `venv/` if it exists, then starts `ui.py`. With the venv already active you can also run `python ui.py`.

**Start mode** (only the fields for the chosen mode are shown):

| Mode | What you provide | What it runs |
|---|---|---|
| **Species name** | Species, optional outgroup, sample size, random seed | NCBI download, then the full pipeline |
| **Batch CSV** | A CSV file (column 1 = species, column 4 = optional outgroup) | One row after another, same as `--csv-file` |
| **Local genome folder** | Species (names the results folder) and a folder of `.fna` / `.fasta` / `.fa` files | Copies those genomes into `input/` and starts at step 2 (Bakta). Sample size and seed are hidden. Start step cannot be 1 |

**Options** (always on the screen; Browse buttons pick files and folders):

| Window field | Same as |
|---|---|
| Working directory | `--workdir` (default: the folder you launched from, resolved to an absolute path) |
| Bakta database folder | `--db` (required for steps 1–2; point at the `db` or `db-light` subfolder) |
| Threads | `--threads` (default 12) |
| Bakta parallel jobs | `--bakta-jobs` (leave blank to use `threads // 4`) |
| Start step | `--start-step` (1–9) |
| Core-alignment FASTA | `--input-fasta` (needed when starting at step 4+ unless that file already exists in the species folder) |
| Outgroup ID | `--outgroup-id` (force a tree leaf) |
| Setup only | `--setup-only` (create folders and stop) |

**Run:** **Run pipeline** starts one job at a time. The log is the same output as the terminal, including the quoted command so you can replay it. **Stop** ends the job and leaves partial output on disk so you can raise **Start step** and continue. Closing the window while a run is going asks before it stops the job.

Results land under `<workdir>/<Species_name>/` (or `<workdir>/` for a CSV batch). If `ecotype_summary.csv` was written, the log prints its path.

**Local genome folder** still contacts NCBI to pick an outgroup when the run includes annotation (start step 2). There is no checkbox for `--no-auto-outgroup`; use the CLI if you need that. For a resume from an existing alignment, use **Start step** 4+ and **Core-alignment FASTA**.

### One species

```bash
python pipeline.py --species "Bacillus subtilis" --db /path/to/bakta_db/db --threads 12
```

### Many species from a CSV

```bash
python pipeline.py --csv-file species.csv --db /path/to/bakta_db/db
```

The first row is treated as a header. **Column 1** is the species and **column 4** is an optional outgroup (an accession or species name). Other columns are ignored. Leave the outgroup blank to auto-detect one.

```csv
Species Name,Genomes,Notes,Outgroup
Bacillus subtilis,,,Bacillus licheniformis
Treponema paraluiscuniculi,,,
```

If one species fails, the pipeline logs the failure and moves on to the next row. A success/failure tally prints at the end.

### Resuming

Every step can be skipped with `--start-step N`. For example, to rebuild the tree and everything after it from an existing core alignment:

```bash
python pipeline.py --species "Bacillus subtilis" --start-step 4 \
  --input-fasta path/to/core_gene_alignment.aln \
  --outgroup-id GCF_000009045.1 --no-auto-outgroup
```

When resuming, pass both `--outgroup-id` and `--no-auto-outgroup` to skip the NCBI lookups. Use `--outgroup-id outgroup` if an outgroup genome was annotated in the original run.

If EcoSim fails or can't be found, the rarefaction FASTAs and rooted tree are **kept**, so `--start-step 8` can pick up from there.

### All options

```
python pipeline.py [OPTIONS]

Input (one required)
  --species TEXT           Species to download, e.g. "Bacillus subtilis"
  --csv-file PATH          Batch file: column 1 = species, column 4 = outgroup

Download (step 1)
  --sample-size INT        Assemblies to randomly sample (default: 200)
  --seed INT               Sampling seed, for reproducible genome sets (default: 42)

Annotation / pan-genome (steps 2–3)
  -d, --db PATH            Bakta database (required when starting at step 1 or 2)
  --bakta-jobs INT         Genomes annotated in parallel (default: threads // 4)

Outgroup / rooting
  --outgroup TEXT          Accession or species name to download as the outgroup
  --no-auto-outgroup       Don't auto-detect an outgroup
  --outgroup-id TEXT       Root on this exact leaf ID (overrides everything else)

Execution
  -w, --workdir PATH       Where per-species folders are created (default: .)
  -t, --threads INT        Total threads for Bakta, Panaroo, tree building and EcoSim (default: 12)
  --start-step {1..9}      Step to start from (default: 1)
  --input-fasta PATH       Core alignment to use when starting at step 4+
  --setup-only             Create the folder structure and exit
```

---

## How the pipeline works

| # | Step | What happens |
|---|---|---|
| 1 | **Download** | Lists every *chromosome* or *complete* assembly for the species and draws a random sample of `--sample-size` with a fixed seed, so the same seed always gives the same genomes. Then downloads the `.fna` files. |
| – | **Outgroup** | Uses `--outgroup` if given. Otherwise picks a RefSeq reference genome from another species in the same genus, falling back to any complete genome in the genus. The outgroup is annotated alongside the other genomes and appears in the alignment as `outgroup`. |
| 2 | **Bakta** | Annotates every genome (up to 201), several at a time. Non-coding RNA searches and plots are skipped because Panaroo only uses coding genes. |
| 3 | **Panaroo** | Builds the pan-genome in `strict` clean mode and a core-gene alignment with MAFFT (`core_gene_alignment.aln`). |
| 4 | **Sort** | Checks that all sequences are the same length and moves the outgroup to the top of the alignment. |
| 5 | **Tree** | Builds a maximum-likelihood tree with VeryFastTree (`-nt -gtr -gamma -nosupport`) on all threads. |
| 6 | **Reroot** | Roots the tree on the outgroup with Biopython. GCA/GCF accession variants are matched automatically. |
| 7 | **Rarefaction** | Creates 100 sub-alignments by concatenating randomly placed 1,000 bp windows: 1, 3, 7, 20 and 100 windows, 20 replicates each. This shows how the ecotype estimate changes with the amount of sequence. |
| 8 | **EcoSim** | Runs `ecosim.jar` (demarcation mode, no GUI, 12 GB heap) on each sub-alignment against the full rooted tree. |
| 9 | **Parse** | Counts the demarcated ecotypes in each EcoSim XML file and writes `ecotype_summary.csv`. |

**How the root is chosen, in order of priority:** `--outgroup-id` → the downloaded outgroup genome → the species' RefSeq reference genome (a proxy for the type strain) → the sequence whose accession number is largest.

---

## Outputs

Each species gets its own folder, `<workdir>/<Species_name>/`. **Intermediates are deleted once the next step has used them**, to keep disk use down on large batches:

| Path | Lifetime |
|---|---|
| `input/*.fna` | Deleted after Bakta |
| `intermediate_bakta/<accession>/` | Deleted after Panaroo |
| `output_roary/results/core_gene_alignment.aln` | Deleted at the end of a full run. Copy it out if you want to keep it |
| `output_roary/panaroo_errors_<timestamp>.log` | Deleted with the folder above |
| `pipeline_temp_<alignment>/` (sorted FASTA, unrooted tree) | Deleted at the end |
| `rerooted_trees/<Species name>.nwk` | Deleted after EcoSim succeeds |
| `rarefaction_fastas_<alignment>/` | Deleted after EcoSim succeeds |
| **`ecosim_output_<alignment>/`** | **Kept: the final results** |

`<alignment>` is the alignment file name without its extension, which is `core_gene_alignment` in a normal run.

The final results folder contains:

```
ecosim_output_core_gene_alignment/
├── sim_species_g1_t1_results.xml     # one EcoSim result per rarefaction replicate
├── ...
├── sim_species_g100_t20_results.xml
└── ecotype_summary.csv               # file,ecotype_count
```

For a spreadsheet of ecotype membership (which strains belong to which ecotype), run:

```bash
python post_processing.py Species_name/ecosim_output_core_gene_alignment
# → .../parsed_results/*.xlsx
```

---

## Performance tuning

| Knob | Effect |
|---|---|
| `--threads` | Used by every step. VeryFastTree and EcoSim scale well with more threads. |
| `--bakta-jobs` | Genomes annotated at the same time. Each job gets `threads / jobs` threads. Lower it if you run out of RAM, since every Bakta job loads its own database indexes. |
| `--sample-size` | Bakta time grows linearly with genome count; Panaroo and tree building grow faster than that. |
| `memory_gb` in `pipeline.py` | EcoSim Java heap (default 12 GB). Raise it for very large alignments. |

Steps 2 and 5 skip genomes and trees whose output already exists, so an interrupted run can be restarted without redoing that work.

---

## Running steps on their own

Every module also works as a standalone script:

```bash
# Put the outgroup (or the largest accession number) first
python move_largest_numeric.py alignment.aln sorted.fasta

# Build trees for every *.fasta in ./tree_rdy_fastas → ./trees_final
python run_trees.py [--fast] [--max-sequences N]

# Root a tree
python reroot_tree.py sorted.fasta unrooted.nwk rooted.nwk [outgroup_id]

# Rarefaction sub-alignments
python Rarefaction_fasta_creation.py sorted.fasta [output_folder]

# EcoSim on a folder of FASTAs
python run_ecosim.py rarefaction_fastas/ --full-tree-path rooted.nwk \
  --output-dir ecosim_results [--memory-gb 24]

# Summaries
python parsing.py                        # counts ecotypes in ./ecosim_results
python post_processing.py ecosim_results # per-ecotype membership spreadsheets
```

---

## Repository layout

```
pipeline.py                    Orchestrator: CLI, steps 1–9, batch mode
ui.py                          Desktop window: same options, live log, Start/Stop
run_ui.sh                      Activate venv if present, then python ui.py
test_ui.py                     Unit tests for the window's validation and commands
download_outgroup.py           NCBI outgroup / type-strain lookup and download
move_largest_numeric.py        Step 4: alignment check and outgroup ordering
run_trees.py                   Step 5: VeryFastTree/FastTree wrapper
reroot_tree.py                 Step 6: outgroup rooting
Rarefaction_fasta_creation.py  Step 7: gene-window sub-alignments
run_ecosim.py                  Step 8: EcoSim batch runner
parsing.py                     Step 9: ecotype counts
post_processing.py             Optional: XML → Excel membership tables
ecosim.jar, bin/               Bundled EcoSim and its native helpers (macOS arm64)
setup.sh                       Environment setup
requirements.txt               Python dependencies
*.ipynb                        Original exploratory notebooks (Windows paths; reference only)
Treponema_paraluiscuniculi/    Example outputs from a small test run
```

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `bash run_ui.sh` / `python ui.py` fails with `No module named '_tkinter'` | The venv's Python was built without Tk. On macOS with Homebrew: `brew install python-tk`, then recreate `venv/` with `bash setup.sh` |
| `python: command not found` from `run_ui.sh` | Run `bash setup.sh` first, or activate `venv/` and use `python ui.py` |
| Local genomes copied but Bakta says no FASTA files | Check **Working directory** — genomes are copied into `<workdir>/<Species_name>/input/` |
| `NCBI 'datasets' CLI tool not found` | `conda install -c conda-forge ncbi-datasets-cli` |
| `--db ... is required` | Starting at step 1 or 2 needs `--db /path/to/bakta_db/db` |
| Bakta fails right away | Check that `--db` points to the `db`/`db-light` subfolder, and that `conda run -n bakta_env bakta --version` works |
| Bakta is killed or the machine swaps | Lower `--bakta-jobs` |
| `Sequences have different lengths` | The input isn't an alignment. Use Panaroo's `core_gene_alignment.aln` or align with MAFFT first |
| `Outgroup ... could not be matched to any tree leaf` | Pass the exact leaf name with `--outgroup-id` (leaf names are the genome file names, e.g. `GCF_000217655.1` or `outgroup`) |
| `gene_length (1000) exceeds alignment length` | The core alignment is too short for rarefaction, usually because too few genes are shared. Check the genome set and outgroup |
| `EcoSim jar not found` / `Java not found` | Install Java 8+. Set `ECOSIM_JAR` if you moved the jar |
| EcoSim fails on Linux or Intel Macs | The bundled `bin/` helpers are Apple Silicon only. Point `ECOSIM_DIR` at a directory with binaries for your platform |
| A batch row failed | Look for `Row N failed` in the output, fix the problem, then rerun that species with `--species` |

---

## Citations

If this pipeline contributes to published work, please cite the tools it runs:

- **Bakta**: Schwengers O. *et al.* (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. *Microbial Genomics*.
- **Panaroo**: Tonkin-Hill G. *et al.* (2020). Producing polished prokaryotic pangenomes with the Panaroo pipeline. *Genome Biology*.
- **VeryFastTree**: Piñeiro C., Abuín J.M., Pichel J.C. (2020). Very Fast Tree: speeding up the estimation of phylogenies for large alignments through parallelization and vectorization strategies. *Bioinformatics*.
- **FastTree 2**: Price M.N., Dehal P.S., Arkin A.P. (2010). FastTree 2: approximately maximum-likelihood trees for large alignments. *PLoS ONE*.
- **Ecotype Simulation**: Koeppel A. *et al.* (2008). Identifying the fundamental units of bacterial diversity: a paradigm shift to incorporate ecology into bacterial systematics. *PNAS*.
- **NCBI Datasets**: O'Leary N.A. *et al.* (2024). Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. *Scientific Data*.
