# Quick Start Guide

The short version of [README.md](README.md). Works on Linux (Arch, Debian/Ubuntu, Fedora) and macOS.

## 🚀 First Time Setup

You need [Miniforge](https://github.com/conda-forge/miniforge) (or any conda), Python 3.11+, and git.

```bash
git clone https://github.com/johnkenjikono/Archive-2.git
cd Archive-2

# Conda envs, venv/, NCBI datasets, VeryFastTree, EcoSim helpers (Linux).
# On Linux it asks for sudo if Java or the compilers are missing.
bash setup.sh

# Bakta database: light (~4 GB unpacked) for a first test, or drop --type light for the full ~70 GB one
conda run -n bakta_env bakta_db download --output ~/bakta_db --type light
```

## 🏃 Running the Pipeline

Activate the venv in every new terminal (`run_ui.sh` does it for you):

```bash
source venv/bin/activate          # fish: source venv/bin/activate.fish
```

```bash
# Desktop UI (no CLI flags)
bash run_ui.sh

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
| `NCBI 'datasets' CLI tool not found` | `source venv/bin/activate` first |
| `setup.sh`: Python too old | Install Python 3.11+, or run `setup.sh` inside a conda env with Python 3.11+ |
| EcoSim fails on Linux | `bash setup.sh` builds the Linux helpers into `tools/linux/bin` (the bundled `tools/bin` is macOS arm64 only) |
| `No module named '_tkinter'` | Install Tk (`tk` on Arch, `python3-tk` on Debian/Ubuntu, `python3-tkinter` on Fedora, `brew install python-tk` on macOS), delete `venv/`, rerun `bash setup.sh` |
| Panaroo finds 0 core genes | Use a closer outgroup or more genomes (`--sample-size`) |
| A step is "Killed" | Out of memory: add `--bakta-jobs 1` or lower `--threads` |

More in the [README's Troubleshooting section](README.md#troubleshooting).
