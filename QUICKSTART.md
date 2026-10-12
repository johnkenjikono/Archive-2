# Quick Start Guide

The short version of [README.md](README.md). This branch targets macOS.

## 🚀 First Time Setup

You need [Miniforge](https://github.com/conda-forge/miniforge) (or any conda) and git. The conda environment provides Python and everything else.

```bash
git clone https://github.com/johnkenjikono/Archive-2.git
cd Archive-2

# One time (Apple Silicon needs the CONDA_SUBDIR prefix; see README)
CONDA_SUBDIR=osx-64 conda env create -f environment.yml

# Bakta database: light (~4 GB unpacked) for a first test, or drop --type light for the full ~70 GB one
conda run -n ecotype bakta_db download --output ~/bakta_db --type light
```

## 🏃 Running the Pipeline

Activate the environment in every new terminal:

```bash
conda activate ecotype
```

```bash
# Desktop UI (no CLI flags)
python ui/ui.py

# Small test run: 10 genomes, about 25 minutes on an 8-core laptop
python pipeline.py --species "Treponema pallidum" --outgroup "Treponema paraluiscuniculi" \
    --sample-size 10 --db ~/bakta_db/db-light --threads 8

# Many species from a CSV (column 1 = species, column 4 = optional outgroup)
python pipeline.py --csv-file species.csv --db ~/bakta_db/db-light

# Resume a stopped run from a later step
python pipeline.py --species "Treponema pallidum" --db ~/bakta_db/db-light --start-step 5

# All options
python pipeline.py --help
```

`--db` takes the `db` or `db-light` subfolder, not `bakta_db` itself.

## 🔧 Running steps on their own

```bash
python steps/move_largest_numeric.py input.fasta output.fasta
python steps/run_trees.py [--fast] [--max-sequences N]
python steps/reroot_tree.py sorted.fasta tree.nwk output.nwk
python steps/dedup_clones.py sorted.fasta sorted.fasta clone_groups.csv   # --threshold, default 1e-5
python steps/rarefaction.py core_gene_alignment.aln                        # needs core_alignment_header.embl alongside
python steps/parsing.py ecosim_output_X                                    # ecotype counts + membership CSV
python steps/visualize_ecotypes.py ecosim_output_X                         # trees colored by ecotype
```

## 📁 Output Files

Everything lands in `<Species_name>/`:

- `ecosim_output_core_gene_alignment/` - EcoSim results, `ecotype_summary.csv`, `rarefaction_by_gene_count.csv`, `ecotype_membership.csv`, and `ecotype_plots/`
- `pangenome/` - Panaroo gene presence/absence table and summary statistics
- `run_parameters.json` - settings, tool and database versions, and which outputs were kept

Choose what is kept with `--keep ITEM` / `--discard ITEM` (e.g. `--keep core_alignment`, `--discard all --keep rerooted_tree`).

## ❌ Troubleshooting

| Issue | Solution |
|-------|----------|
| `VeryFastTree/FastTree not found` or `datasets` not found | `conda activate ecotype` |
| EcoSim fails on Intel Mac | The bundled `tools/bin` helpers are Apple Silicon only |
| `No module named '_tkinter'` | The `ecotype` environment is not active. Run `conda activate ecotype` |
| Panaroo finds 0 core genes | Use a closer outgroup or more genomes (`--sample-size`) |
| A step is "Killed" | Out of memory: add `--bakta-jobs 1` or lower `--threads` |

More in the [README's Troubleshooting section](README.md#troubleshooting).
