# Ecotype Discovery Pipeline

End-to-end bioinformatics pipeline for bacterial ecotype discovery:

**NCBI download → Bakta annotation → Roary pan-genome → FastTree phylogeny → EcoSim demarcation**

---

## Prerequisites

| Requirement | Notes |
|---|---|
| [Conda / Miniforge](https://github.com/conda-forge/miniforge) | Required for Bakta and Roary environments |
| Python ≥ 3.10 | With `venv` support |
| [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) | For genome downloads (Step 1) |
| VeryFastTree or FastTree | For tree building (Step 5) |
| Java JRE/JDK | For EcoSim (Step 8) |
| Bakta database | Downloaded separately (see setup below) |
| EcoSim `.jar` | **Bundled** as `ecosim.jar` in repo root — set `ECOSIM_JAR` env var to point to it |

---

## Installation

```bash
bash setup.sh
```

This creates two conda environments (`bakta_env`, `roary_env`), a Python `venv`, and checks all tool dependencies.

**Activate the venv before running the pipeline:**
```bash
source venv/bin/activate        # bash/zsh
source venv/bin/activate.fish   # fish
```

**Download the Bakta database** (one-time, ~30 GB):
```bash
conda activate bakta_env
bakta_db download --output /path/to/bakta_db
conda deactivate
```

**Configure EcoSim** (skip if not running Step 8):

`ecosim.jar` is bundled at the repo root, so you can point directly to it:
```bash
export ECOSIM_JAR="$(pwd)/ecosim.jar"
export ECOSIM_DIR="$(pwd)"             # working dir for the jar
```

Or set a custom path if you have your own copy:
```bash
export ECOSIM_JAR=/path/to/ecosim.jar
export ECOSIM_DIR=/path/to/ecosim/
```

---

## Running the Pipeline

### Graphical interface

```bash
bash run_ui.sh
```

Or, with the venv already active: `python ui.py`.

Pick **Species name**, **Batch CSV**, or **Local genome folder**, fill in the same options as the CLI, and press **Run pipeline**. The log is the same output as the terminal. Results land under `<workdir>/<Species_name>/` (or `<workdir>/` for a CSV batch).

### Single species

```bash
python pipeline.py \
  --species "Bacillus subtilis" \
  --db /path/to/bakta_db \
  --threads 12
```

### Batch mode (CSV)

```bash
python pipeline.py --csv-file batch.csv --db /path/to/bakta_db
```

CSV format (header row skipped, uses columns 1 and 4):
```
Species Name,col2,col3,Outgroup
Bacillus subtilis,,,Bacillus licheniformis
Treponema paraluiscuniculi,,,
```

### Resume from a later step

```bash
# Resume from tree building (skips download, Bakta, Roary)
python pipeline.py --species "Bacillus subtilis" --start-step 4

# Resume from EcoSim using an existing alignment
python pipeline.py --species "Bacillus subtilis" --start-step 8 \
  --input-fasta /path/to/core_gene_alignment.aln
```

---

## Pipeline Steps

| Step | Name | What it does |
|---|---|---|
| 1 | Download | Samples up to 200 chromosome/complete assemblies from NCBI, downloads `.fna` files |
| 2 | Bakta | Annotates each genome, produces `.gff3` files (capped at 201 genomes) |
| 3 | Roary | Builds pan-genome, outputs `core_gene_alignment.aln` (capped at 201 GFF files) |
| 4 | Sort | Moves outgroup sequence to top of alignment FASTA |
| 5 | FastTree | Builds unrooted phylogenetic tree (`-gtr -gamma -boot 50`) |
| 6 | Reroot | Roots tree on outgroup using BioPython |
| 7 | Rarefaction | Creates subsampled FASTAs for gene-count sensitivity analysis |
| 8 | EcoSim | Runs Ecotype Simulation (Java) on each rarefaction FASTA |
| 9 | Parse | Counts ecotypes across EcoSim XML outputs |

---

## All CLI Flags

```
python pipeline.py [OPTIONS]

Input:
  --species TEXT          Species name for NCBI download (e.g. "Bacillus subtilis")
  --csv-file PATH         CSV file for batch processing (columns 1=species, 4=outgroup)

Download:
  --sample-size INT       Max assemblies to download (default: 200)

Annotation:
  -d, --db PATH           Path to Bakta database (required for steps 1-2)

Execution:
  -w, --workdir PATH      Base working directory (default: .)
  -t, --threads INT       Threads for Bakta and Roary (default: 12)
  --start-step INT        Resume from step 1-9 (default: 1)
  --setup-only            Create directories and exit

Tree:
  --input-fasta PATH      Core alignment FASTA (required if --start-step >= 4)
  --outgroup-id TEXT      Force a specific outgroup leaf ID for tree rooting
  --outgroup TEXT         NCBI accession or species name to download as outgroup
  --no-auto-outgroup      Disable automatic outgroup detection
```

---

## Output Structure

All output for a run is placed under `<workdir>/<Species_name>/`:

```
<Species_name>/
├── input/                          # Downloaded .fna files
├── intermediate_bakta/             # Per-genome Bakta annotation dirs
│   └── <accession>/
│       └── <accession>.gff3
├── output_roary/
│   └── results/
│       └── core_gene_alignment.aln # Core gene alignment (input to tree step)
├── rerooted_trees/
│   └── <Species_name>.nwk          # Final rooted tree
├── rarefaction_fastas_<species>/   # Subsampled FASTAs for EcoSim
└── ecosim_output_<species>/        # EcoSim XML results
```

Temporary directories (`pipeline_temp_<species>/`) are cleaned up automatically after a successful run.

---

## Known Limits & Gotchas

- **Roary hard cap:** 201 GFF files maximum. The pipeline truncates automatically if more are found.
- **Sequences must be aligned** before tree building (identical length per sequence). Roary's `core_gene_alignment.aln` satisfies this automatically. If providing your own `--input-fasta`, run MAFFT first: `mafft input.fasta > aligned.fasta`
- **Outgroup auto-detection:** If `--outgroup` is not provided, the pipeline queries NCBI for a reference genome from the same genus. If none is found, the typestrain is used as root instead.
- **EcoSim memory:** Defaults to 12 GB Java heap (`-Xmx12G`). For large datasets, set a larger value by editing `run_ecosim.py`.

---

## Running Individual Scripts

Each module can also be run standalone:

```bash
# Build trees from an existing aligned FASTA directory
python run_trees.py [--fast] [--max-sequences N]

# Run EcoSim batch
python run_ecosim.py <fasta_dir> --full-tree-path <tree.nwk> --output-dir <out>

# Reroot a tree manually
python reroot_tree.py sorted.fasta unrooted.nwk rerooted.nwk [outgroup_id]

# Create rarefaction FASTAs
python Rarefaction_fasta_creation.py input.fasta [output_folder]

# Parse EcoSim results
python parsing.py   # reads from ./ecosim_results/ by default
```
