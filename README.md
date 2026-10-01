# Ecotype Discovery Pipeline

Give it a bacterial species name. It downloads complete genomes from NCBI, builds a core-genome phylogeny, and runs **Ecotype Simulation (EcoSim)** to estimate how many ecotypes the species contains, and how stable that estimate is as more genes are sampled.

```
NCBI Datasets ─► Bakta ─► Panaroo ─► VeryFastTree ─► reroot ─► rarefaction ─► EcoSim ─► ecotype counts
  (download)    (annotate) (core genes)  (tree)     (outgroup)  (gene subsets)   (demarcate)   (CSV)
```

Run it from the desktop window (`bash run_ui.sh`) or from the terminal, one species at a time or a whole CSV of species.

---

## Contents

- [The idea in plain language](#the-idea-in-plain-language)
- [Getting started, step by step](#getting-started-step-by-step)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
  - [Graphical interface](#graphical-interface)
  - [One species](#one-species)
  - [Many species from a CSV](#many-species-from-a-csv)
  - [On Modal (cloud Linux)](#on-modal-cloud-linux)
- [How the pipeline works](#how-the-pipeline-works)
- [Outputs](#outputs)
- [Performance tuning](#performance-tuning)
- [Running steps on their own](#running-steps-on-their-own)
- [Repository layout](#repository-layout)
- [Troubleshooting](#troubleshooting)
- [Citations](#citations)

---

## The idea in plain language

**The question.** A bacterial species like *Streptococcus pyogenes* is not one uniform population. Inside it there are groups of cells that live in slightly different ways (different hosts, tissues, or niches) and that compete mainly with their own group. These groups are called **ecotypes**. This pipeline asks: *how many ecotypes does this species contain, and which genomes belong to which?*

**How we answer it, in everyday terms.**

| Step | What it does | Analogy |
|---|---|---|
| Download | Fetches complete genomes of the species from NCBI, plus one close relative (the **outgroup**) | Collecting specimens, plus one cousin to use as a reference point |
| Annotate (Bakta) | Finds and labels every gene in each genome | Labeling every sentence in every book |
| Core genes (Panaroo) | Finds the genes that **every** genome has and lines them up | Keeping only the chapters that appear in every book |
| Tree (VeryFastTree) | Builds a family tree from differences in those shared genes | Drawing who is related to whom |
| Reroot | Places the root of the tree at the outgroup | Deciding which end of the tree is "oldest" |
| Rarefaction | Repeats the analysis using only 1, 3, 7, 20 or 100 random genes | Checking whether a small sample gives the same answer as the whole book |
| EcoSim | Cuts the tree into ecotypes using a model of how populations diversify and are purged by selection | The actual counting step |
| Summary and plots | Writes tables and colored trees | The figures for your report |

**Why rarefaction matters.** If the ecotype count stays the same whether you use 7 genes or 100, the estimate is stable. If it keeps climbing, you have not sampled enough of the genome to trust it yet.

**Words you will see.**

- **Core genome**: genes present in every genome analysed.
- **Outgroup**: a related but different species, used only to anchor the root of the tree. It is never counted as part of an ecotype.
- **Clones**: genomes with (nearly) identical core genes. The pipeline merges them so they do not inflate the tree, then gives each clone its representative's ecotype in the final tables.
- **Replicate**: one random draw of genes. Each replicate gets its own tree and its own EcoSim run.
- **Terminal / command line**: the text window where you type commands. Everything below that starts with `bash`, `python`, `conda` or `source` goes there.

---

## Getting started, step by step

This assumes macOS or Linux and no prior setup. Budget about 30 minutes for installing, plus the database download. Setup commands are run once. Only steps 4 and 5 repeat each session.

### 0. Check you have what you need

| Need | Why | How to check |
|---|---|---|
| A terminal | Setup and optional command-line use | Open **Terminal** (macOS) or your Linux terminal |
| **Conda or Miniforge** | Installs Bakta and Panaroo | `conda --version`. If not found, install [Miniforge](https://github.com/conda-forge/miniforge) first, then open a new terminal |
| Python 3.10 or newer | Runs the pipeline | `python3 --version` |
| Git | Gets the code | `git --version` |
| Internet | Downloads genomes from NCBI | n/a |
| Disk space | Full Bakta database is about 70 GB; `light` is much smaller (less precise). Each run needs more for genomes and intermediates | `df -h .` |
| RAM and cores | 8 cores and 16 GB is enough for a small test run | `nproc` (Linux) or `sysctl -n hw.ncpu` (macOS) |

### 1. Get the code

```bash
git clone https://github.com/johnkenjikono/Archive-2.git
cd Archive-2
```

All commands below are run from inside the new `Archive-2` folder.

### 2. Install everything (one time)

```bash
bash setup.sh
```

This takes several minutes. It creates three conda environments (`bakta_env`, `panaroo_env`, `ecotools`) and a Python environment (`venv/`), and on Linux it also builds EcoSim's helper programs. On Arch Linux it asks for your password to install Java and a compiler. Lines starting with `✓` are good. A `⚠` line tells you what is missing and usually how to fix it. It is safe to rerun `setup.sh`: it skips anything that already exists.

Check it worked:

```bash
source venv/bin/activate        # fish shell: source venv/bin/activate.fish
datasets --version              # NCBI download tool
java -version                   # needed by EcoSim
conda run -n bakta_env bakta --version
```

All three should print a version number, not an error.

### 3. Download the Bakta database (one time)

Bakta needs a reference database to name genes. This is the largest download.

```bash
conda run -n bakta_env bakta_db download --output ~/bakta_db               # full, ~70 GB, recommended
# conda run -n bakta_env bakta_db download --output ~/bakta_db --type light  # smaller, good for a first test
```

When it finishes, the database is in `~/bakta_db/db` (or `~/bakta_db/db-light`). **Remember this path.** You give the `db` or `db-light` folder to the pipeline, not `bakta_db` itself.

### 4. Activate the environment (every new terminal)

```bash
cd Archive-2
source venv/bin/activate        # fish shell: source venv/bin/activate.fish
```

Your prompt now starts with `(venv)`. If you skip this, you will see `NCBI 'datasets' CLI tool not found`. `bash run_ui.sh` does this for you.

### 5. Run a small test first

Do not start with the default of 200 genomes: that takes hours. Use **10 genomes** on a species with a close, well-known relative. This one takes about 25 minutes on an 8-core laptop with `db-light`.

**Option A: the window (no commands to remember)**

```bash
bash run_ui.sh
```

1. Set **Start mode** to **Species name**.
2. Species: `Treponema pallidum`. Outgroup: `Treponema paraluiscuniculi`.
3. **Bakta database folder**: click Browse and choose `~/bakta_db/db` (or `db-light`).
4. **Sample size**: `10`.
5. Click **Run pipeline**. The log scrolls as each step runs. **Stop** ends the run and keeps what was written.

**Option B: the terminal**

```bash
python pipeline.py \
  --species "Treponema pallidum" \
  --outgroup "Treponema paraluiscuniculi" \
  --sample-size 10 \
  --db ~/bakta_db/db \
  --threads 8
```

Set `--threads` to the number of cores you can spare. If your machine slows to a crawl or a step is "Killed", add `--bakta-jobs 1`.

**What you should see:** the log announces each step in order (download, Bakta, Panaroo, tree, rarefaction, EcoSim, parsing, plotting). Bakta and EcoSim take the longest. If the run is interrupted, rerun the same command: finished genomes and trees are skipped.

### 6. Read the results

Everything is in a folder named after the species, inside the folder you ran from:

```
Treponema_pallidum/ecosim_output_core_gene_alignment/
```

Start with these:

| File | What it tells you |
|---|---|
| `ecotype_plots/full_core_genome_results.png` | The tree with each ecotype in its own color. This is the main figure |
| `ecotype_plots/full_core_genome_results_pie.png` | How many genomes fall in each ecotype |
| `ecotype_summary.csv` | Ecotype count for every replicate, including the full core genome |
| `rarefaction_by_gene_count.csv` | How the count changes with 1, 3, 7, 20, 100 genes. Stable numbers mean a trustworthy estimate |
| `ecotype_membership.csv` | Which genome belongs to which ecotype |

Open the CSVs in Excel, Numbers or LibreOffice. The full list of files is under [Outputs](#outputs).

### 7. Run your own species

1. Pick a species with **several complete genomes** on NCBI and a close relative to use as the outgroup. With very few genomes and a distant outgroup, Panaroo can find 0 core genes and the run stops before the tree is built.
2. Run the same command or window settings with your species. Raise `--sample-size` gradually (10, then 50, then more). Time and memory grow with every extra genome.
3. For many species at once, use a CSV: see [Many species from a CSV](#many-species-from-a-csv).
4. To resume a run that stopped partway through, raise **Start step** or add `--start-step N`: see [Resuming](#resuming).

### If something goes wrong

Read the last red or error line of the log, then look it up in [Troubleshooting](#troubleshooting). The most common causes are: environment not activated (step 4), `--db` pointing at `bakta_db` instead of `bakta_db/db`, and a species with too few genomes.

---

## Requirements

| Tool | Used in | Installed by `setup.sh`? |
|---|---|---|
| [Conda / Miniforge](https://github.com/conda-forge/miniforge) | Bakta and Panaroo environments | No, install it first |
| Python ≥ 3.10 | The pipeline and the desktop window (stdlib `tkinter`) | Creates `venv/` from `requirements.txt` |
| [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) (`datasets`) | Step 1, outgroup/type-strain lookup | Yes: `ecotools` conda env, linked into `venv/bin` |
| [Bakta](https://github.com/oschwengers/bakta) + its database | Step 2 | Creates `bakta_env`. You download the database |
| [Panaroo](https://github.com/gtonkinhill/panaroo) + MAFFT | Step 3 | Creates `panaroo_env` |
| [VeryFastTree](https://github.com/citiususc/veryfasttree) (or FastTree) | Step 5 | Linux: `ecotools` conda env, linked into `venv/bin`. macOS: `brew` |
| Java 8+ | Step 8 (EcoSim) | Arch: `pacman` (with `gcc-fortran` for the EcoSim helpers). Elsewhere only checked |

**EcoSim is bundled.** `tools/ecosim.jar` and its helper binaries in `tools/bin/` ship with the repo. The binaries in `tools/bin/` are **macOS arm64 (Apple Silicon)** builds. On Linux, `setup.sh` builds them into `tools/linux/bin`, which the pipeline picks up automatically. On any other platform, supply EcoSim binaries built for that platform (see [Configuration](#configuration)).

---

## Installation

### 1. Environments and dependencies

```bash
bash setup.sh
```

This script:
- on Arch Linux, installs `gcc-fortran`, `jre-openjdk`, `git` and `base-devel` with `pacman` if any are missing (asks for sudo)
- creates the `bakta_env`, `panaroo_env` and `ecotools` (`datasets`, VeryFastTree) conda environments from conda-forge + bioconda only (skips any that already exist)
- creates `venv/`, installs `requirements.txt`, and links `datasets` and `veryfasttree` into `venv/bin`
- on Linux, builds EcoSim's helpers into `tools/linux/bin`

### 2. Bakta database (one-time)

```bash
conda run -n bakta_env bakta_db download --output bakta_db                # full database (recommended, ~70 GB unpacked)
# conda run -n bakta_env bakta_db download --output bakta_db --type light  # much smaller, less precise
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
| `ECOSIM_JAR` | `tools/ecosim.jar` | Path to the EcoSim jar |
| `ECOSIM_DIR` | `tools/` | Working directory for EcoSim (must contain its `bin/` helpers) |

---

## Usage

### Graphical interface

The window collects the same options as `pipeline.py` and runs it for you. You do not need to assemble a command.

```bash
bash run_ui.sh
```

That script activates `venv/` if it exists, then starts `ui/ui.py`. With the venv already active you can also run `python ui/ui.py`.

**Start mode** (only the fields for the chosen mode are shown):

| Mode | What you provide | What it runs |
|---|---|---|
| **Species name** | Species, optional outgroup, sample size, random seed | NCBI download, then the full pipeline |
| **Batch CSV** | A CSV file (column 1 = species, column 4 = optional outgroup) | One row after another, same as `--csv-file`. **Fill blank outgroups** runs `fill_outgroups.py` first (add `--llm` with **Use Claude (LLM)**) and switches the CSV field to the filled file for review |
| **Local genome folder** | Species (names the results folder) and a folder of `.fna` / `.fasta` / `.fa` files | Copies those genomes into `input/` and starts at step 2 (Bakta). Sample size and seed are hidden. Start step cannot be 1 |

**Options** (always on the screen; Browse buttons pick files and folders):

| Window field | Same as |
|---|---|
| Working directory | `--workdir` (default: the folder you launched from, resolved to an absolute path) |
| Bakta database folder | `--db` (required for steps 1–2; point at the `db` or `db-light` subfolder) |
| Threads | `--threads` (default 12) |
| Bakta parallel jobs | `--bakta-jobs` (leave blank to use `threads // 4`) |
| Start step | `--start-step` (1–10) |
| Core-alignment FASTA | `--input-fasta` (needed when starting at step 4+ unless that file already exists in the species folder) |
| Outgroup ID | `--outgroup-id` (force a tree leaf) |
| Setup only | `--setup-only` (create folders and stop) |
| **Analysis** box | Clone threshold (`--clone-threshold`), Full core genome comparison (`--no-full-core` when unchecked), every gene at the 1-gene level (`--each-gene`), replicates per gene count (`--replicates`) |
| **Outputs to keep** box | One checkbox per `--keep` item (see [Outputs](#outputs)). Unchecking a default adds `--discard`, checking an extra one adds `--keep` |

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

To choose outgroups once and review them before a batch run, fill the blanks ahead of time:

```bash
python fill_outgroups.py species.csv          # NCBI heuristic (deterministic)
python fill_outgroups.py species.csv --llm    # Claude picks from NCBI genomes in the genus
```

This writes `species_outgroups.csv`. Rows that already have an outgroup are left alone. A new **Outgroup Source** column records how each outgroup was chosen (heuristic, or the model name and its reason), and the pipeline ignores that column. `--llm` needs `pip install anthropic` and `ANTHROPIC_API_KEY`, and costs about $0.01–0.06 per species. If the Claude call fails, the script uses the heuristic for that row and records that in the source column.

### On Modal (cloud Linux)

`modal_app.py` runs the same `pipeline.py` in a Linux container on [Modal](https://modal.com). The image holds both conda envs, `datasets`, VeryFastTree, Java and Linux builds of the EcoSim helpers (built by `tools/build_ecosim_linux.sh`). The Bakta database and all results live in Modal Volumes (`bakta-db`, `ecotype-results`), so they survive between runs and `--start-step` resumes work.

```bash
pip install -r requirements.txt && modal setup                 # one-time: installs the client, logs in
modal run modal_app.py::download_bakta_db                      # one-time; add --light for db-light (then pass --db db-light)

modal run modal_app.py --species "Treponema pallidum" --outgroup "Treponema paraluiscuniculi" --extra "--sample-size 10"
modal run modal_app.py --csv-file species.csv                  # one container per row, all in parallel

modal volume get ecotype-results Treponema_pallidum/ecosim_output_core_gene_alignment .
```

`--extra` passes any other `pipeline.py` flags, e.g. `--extra "--start-step 8 --outgroup-id outgroup"`. Each container gets 16 CPUs and 64 GB RAM (`CPUS` / `MEMORY_MB` at the top of `modal_app.py`) and a 24-hour limit, which is Modal's maximum. If a run fails, whatever it wrote is kept in the volume.

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
  --sample-size INT        Assemblies to randomly sample (default: 200; more = longer runs, more memory)
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
  --start-step {1..10}     Step to start from (default: 1)
  --no-dedup               Keep every genome (skip clone collapsing in step 4)
  --clone-threshold FLOAT  Max per-site divergence for two genomes to count as clones (default: 1e-5; 0 = identical only)
  --input-fasta PATH       Core alignment to use when starting at steps 4-7 (core_alignment_header.embl must be next to it)
  --setup-only             Create the folder structure and exit

Rarefaction (steps 7-9)
  --gene-counts N [N ...]  Core genes per replicate (default: 1 3 7 20 100)
  --replicates INT         Random replicates per gene count (default: 20); each one is an EcoSim run
  --each-gene              At the 1-gene level, run every core gene once instead of random replicates
  --no-full-core           Skip the EcoSim run on the complete core genome (and the comparison to it)

Outputs (see "Outputs")
  --keep ITEM[,ITEM]       Also keep these (repeatable; "all" keeps everything)
  --discard ITEM[,ITEM]    Drop these defaults (repeatable; "all" drops everything). --keep wins
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
| 4b | **Collapse clones** | Genomes whose core alignments differ at no more than 10⁻⁵ of sites (`--clone-threshold`; 0 means identical only) are collapsed to one representative. A gap against a base counts as a difference. Groups are built greedily in alignment order: each genome joins the closest representative within the threshold, so members are within the threshold of their representative (and at most twice it of each other). `clone_groups.csv` records each genome's differences from its representative. The root genome is never collapsed, and genomes identical to it form their own group, because EcoSim never places the outgroup in an ecotype. Identical sequences give EcoSim zero-length branches and extra taxa without adding information. Every genome is listed in `clone_groups.csv`, and collapsed genomes get their representative's ecotype in the parsed results. `--no-dedup` turns this off. |
| 5 | **Tree** | Builds a maximum-likelihood tree with VeryFastTree (`-nt -gtr -gamma -nosupport`) on all threads. |
| 6 | **Reroot** | Roots the tree on the outgroup with Biopython. GCA/GCF accession variants are matched automatically. |
| 7 | **Rarefaction** | Creates 100 sub-alignments (change with `--gene-counts` / `--replicates`; `--each-gene` runs every core gene once at the 1-gene level) by concatenating whole core genes drawn at random without replacement: 1, 3, 7, 20 and 100 genes, 20 replicates each. Gene boundaries come from Panaroo's `core_alignment_header.embl`, which must sit next to the alignment (checked before step 4 starts). Genes that are all gaps in any genome, which is how Panaroo pads a gene a genome lacks, are not sampled. The genes in each replicate are listed in `rarefaction_genes.csv`. This shows how the ecotype estimate changes with the number of genes. The whole core alignment is also added as one more replicate, `full_core_genome` (every core gene), so each subset can be compared with it; its tree is the step-5 tree. `--no-full-core` turns this off. |
| 7b | **Rarefaction trees** | Builds and reroots one tree per sub-alignment. EcoSim's binning is driven by the tree, so reusing the full-alignment tree for every replicate made the rarefaction curve flat by construction. |
| 8 | **EcoSim** | Runs `ecosim.jar` (demarcation mode, no GUI, 6 GB heap) on each sub-alignment against its own rooted tree, falling back to the full rooted tree for any replicate whose tree could not be built or rerooted. |
| 9 | **Parse** | Counts the demarcated ecotypes in each EcoSim XML file and writes `ecotype_summary.csv`, plus `ecotype_membership.csv` listing which genome is in which ecotype for every replicate (clones expanded). It then compares every replicate with the full core genome: `rarefaction_comparison.csv` (per replicate), `rarefaction_by_gene_count.csv` (per gene count) and `ecotype_support.csv` (how often each full-core ecotype is recovered exactly at each gene count). |
| 10 | **Plot** | Draws the tree EcoSim used for the full core genome (or, with `--no-full-core`, each 100-gene replicate) with every ecotype in its own color (`ecotype_plots/*.png`), plus an iTOL `TREE_COLORS` file per tree for the [iTOL](https://itol.embl.de) viewer, and a pie chart of genomes per ecotype (`*_pie.png`) in the same colors. |

**How the root is chosen, in order of priority:** `--outgroup-id` → the downloaded outgroup genome → the species' RefSeq reference genome (a proxy for the type strain) → the sequence whose accession number is largest.

---

## Outputs

Each species gets its own folder, `<workdir>/<Species_name>/`. Intermediates are deleted once the next step has used them, to keep disk use down on large batches. **What is kept is set with `--keep` and `--discard`** (or the **Outputs to keep** checkboxes in the window):

| Item | Default | What it keeps |
|---|---|---|
| `genomes` | deleted after Bakta | `input/*.fna` |
| `annotations` | deleted after Panaroo | `intermediate_bakta/<accession>/` |
| `pangenome` | **kept** | `pangenome/gene_presence_absence.csv`, `.Rtab`, `gene_presence_absence_roary.csv`, `summary_statistics.txt` |
| `core_alignment` | deleted at the end | `pangenome/core_gene_alignment.aln` + `core_alignment_header.embl`. A later `--start-step 4` finds it there |
| `panaroo` | deleted at the end | the whole `output_roary/` folder |
| `rerooted_tree` | **kept** | `rerooted_trees/<Species name>.nwk`, the rooted full core-genome tree |
| `rarefaction` | deleted after EcoSim | `rarefaction_fastas_<alignment>/`, `rarefaction_trees_<alignment>/` |
| `ecosim_xml` | **kept** | the EcoSim XML per replicate. Needed to rerun steps 9–10 |
| `ecotypes` | **kept** | `ecotype_membership.csv` and the three full-core comparison tables |
| `tree_plots` | **kept** | `ecotype_plots/*.png` + iTOL color files |
| `pie_charts` | **kept** | `ecotype_plots/*_pie.png` |

`ecotype_summary.csv`, `clone_groups.csv` and `rarefaction_genes.csv` are always kept, and `pipeline_temp_<alignment>/` (sorted FASTA, unrooted tree) is always deleted at the end. Every run also writes `run_parameters.json` with its settings and the kept/discarded items.

```bash
python pipeline.py --species "Bacillus subtilis" --db bakta_db/db --keep core_alignment   # defaults + the alignment
python pipeline.py --species "Bacillus subtilis" --db bakta_db/db --discard all --keep rerooted_tree,ecotypes
python pipeline.py --species "Bacillus subtilis" --db bakta_db/db --keep all               # keep every intermediate
```

`<alignment>` is the alignment file name without its extension, which is `core_gene_alignment` in a normal run.

The final results folder contains:

```
ecosim_output_core_gene_alignment/
├── sim_species_g1_t1_results.xml     # one EcoSim result per rarefaction replicate
├── ...
├── sim_species_g100_t20_results.xml
├── full_core_genome_results.xml      # EcoSim on every core gene
├── ecotype_summary.csv               # file,ecotype_count
├── rarefaction_comparison.csv        # per replicate vs. full core: count difference, adjusted Rand index, ecotypes recovered
├── rarefaction_by_gene_count.csv     # the same, averaged per gene count
├── ecotype_support.csv               # per full-core ecotype: share of replicates that recover it exactly, per gene count
├── ecotype_membership.csv            # file,gene_count,trial,ecotype,ecotype_size,taxon,representative
├── clone_groups.csv                  # representative,member,group_size (one row per genome)
├── rarefaction_genes.csv             # which core genes went into each replicate
└── ecotype_plots/
    ├── full_core_genome_results.png               # tree colored by ecotype
    ├── full_core_genome_results_itol_colors.txt
    └── full_core_genome_results_pie.png           # genomes per ecotype
```

To plot other replicates, or to re-plot after the run:

```bash
python steps/visualize_ecotypes.py Species_name/ecosim_output_core_gene_alignment --gene-count 20
```

`ecotype_membership.csv` already lists which genome is in which ecotype, including collapsed clones. For a spreadsheet with one row per ecotype (it lists only the representatives EcoSim saw, not collapsed clones), run:

```bash
python steps/post_processing.py Species_name/ecosim_output_core_gene_alignment
# → .../parsed_results/*.xlsx
```

---

## Performance tuning

| Knob | Effect |
|---|---|
| `--threads` | Used by every step. VeryFastTree and EcoSim scale well with more threads. |
| `--bakta-jobs` | Genomes annotated at the same time. Each job gets `threads / jobs` threads. Lower it if you run out of RAM, since every Bakta job loads its own database indexes. |
| `--sample-size` | Bakta time grows linearly with genome count; Panaroo and tree building grow faster than that. |
| `memory_gb` in `pipeline.py` | EcoSim Java heap (default 6 GB). Raise it for very large alignments. |

Steps 2 and 5 skip genomes and trees whose output already exists, so an interrupted run can be restarted without redoing that work.

---

## Running steps on their own

Every module also works as a standalone script:

```bash
# Put the outgroup (or the largest accession number) first
python steps/move_largest_numeric.py alignment.aln sorted.fasta

# Build trees for every *.fasta in ./tree_rdy_fastas → ./trees_final
python steps/run_trees.py [--fast] [--max-sequences N]

# Root a tree
python steps/reroot_tree.py sorted.fasta unrooted.nwk rooted.nwk [outgroup_id]

# Rarefaction sub-alignments
python steps/rarefaction.py sorted.fasta [output_folder]

# EcoSim on a folder of FASTAs
python steps/run_ecosim.py rarefaction_fastas/ --full-tree-path rooted.nwk \
  --output-dir ecosim_results [--memory-gb 24]

# Summaries
python steps/parsing.py                        # counts ecotypes in ./ecosim_results
python steps/post_processing.py ecosim_results # per-ecotype membership spreadsheets
```

---

## Repository layout

```
pipeline.py        Orchestrator: CLI, steps 1–10, batch mode (start here)
run_ui.sh          Activate venv if present, then launch the desktop window
setup.sh           Environment setup
modal_app.py       Run the pipeline on Modal (see "On Modal")
requirements.txt   Python dependencies

steps/             One module per pipeline step, each also runnable on its own
  download_outgroup.py      NCBI outgroup / type-strain lookup and download
  move_largest_numeric.py   Step 4: alignment check and outgroup ordering
  run_trees.py              Step 5: VeryFastTree/FastTree wrapper
  reroot_tree.py            Step 6: outgroup rooting
  rarefaction.py            Step 7: gene-window sub-alignments
  run_ecosim.py             Step 8: EcoSim batch runner
  parsing.py                Step 9: ecotype counts
  post_processing.py        Optional: XML → Excel membership tables

ui/                Desktop window
  ui.py                     Same options as pipeline.py, live log, Start/Stop
  test_ui.py                Unit tests for the window's validation and commands

tools/             Bundled EcoSim
  ecosim.jar, bin/          EcoSim and its native helpers (macOS arm64)
  build_ecosim_linux.sh     Builds the Linux helpers (used by setup.sh and modal_app.py)
```

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `bash run_ui.sh` / `python ui/ui.py` fails with `No module named '_tkinter'` | The venv's Python was built without Tk. On macOS with Homebrew: `brew install python-tk`, then recreate `venv/` with `bash setup.sh` |
| `python: command not found` from `run_ui.sh` | Run `bash setup.sh` first, or activate `venv/` and use `python ui/ui.py` |
| Local genomes copied but Bakta says no FASTA files | Check **Working directory** — genomes are copied into `<workdir>/<Species_name>/input/` |
| `NCBI 'datasets' CLI tool not found` | Activate the venv (`source venv/bin/activate`); `setup.sh` links `datasets` there. Or rerun `bash setup.sh` |
| `--db ... is required` | Starting at step 1 or 2 needs `--db /path/to/bakta_db/db` |
| Bakta fails right away | Check that `--db` points to the `db`/`db-light` subfolder, and that `conda run -n bakta_env bakta --version` works |
| Bakta is killed or the machine swaps | Lower `--bakta-jobs` |
| `Sequences have different lengths` | The input isn't an alignment. Use Panaroo's `core_gene_alignment.aln` or align with MAFFT first |
| `Outgroup ... could not be matched to any tree leaf` | Pass the exact leaf name with `--outgroup-id` (leaf names are the genome file names, e.g. `GCF_000217655.1` or `outgroup`) |
| `gene_length (1000) exceeds alignment length` | The core alignment is too short for rarefaction, usually because too few genes are shared. Check the genome set and outgroup |
| `EcoSim jar not found` / `Java not found` | Install Java 8+. Set `ECOSIM_JAR` if you moved the jar |
| EcoSim fails on Linux or Intel Macs | The bundled `tools/bin/` helpers are Apple Silicon only. On Linux, `bash setup.sh` builds them into `tools/linux/bin` (via `tools/build_ecosim_linux.sh`). Elsewhere, point `ECOSIM_DIR` at a directory with binaries for your platform |
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
