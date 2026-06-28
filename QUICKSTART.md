# Quick Start Guide

## 🚀 First Time Setup

```bash
cd /Users/carlofedolfi/Downloads/move_largest_numeric

# Run setup script (optional - you already have VeryFastTree installed)
bash setup.sh

# Activate virtual environment
source venv/bin/activate.fish
```

## ✅ Before Running

1. **Ensure sequences are aligned** (all same length)
   - If not aligned, use: `mafft input.fasta > aligned.fasta`

2. **VeryFastTree is already installed** ✓
   - If you need to reinstall: `brew install veryfasttree`

3. **Set EcoSim jar path** (only if using EcoSim step)
   ```bash
   export ECOSIM_JAR=/path/to/ecosim.jar
   ```
   If needed, also set the working directory:
   ```bash
   export ECOSIM_DIR=/path/to/ecosim
   ```

## 🏃 Running the Pipeline

```bash
# Full pipeline (all steps)
python pipeline.py your_sequences.fasta

# Resume from step 3 (if previous steps completed)
python pipeline.py your_sequences.fasta --start-step 3

# Run individual tools
python move_largest_numeric.py input.fasta output.fasta
python run_trees.py                          # Build trees
python reroot_tree.py sorted.fasta tree.nwk output.nwk
```

## 📊 What to Expect

**For 100-1000 sequences:**
- Sorting & validation: seconds
- Tree building: minutes (depends on sequence length & alignment)
- Rooting: seconds
- EcoSim: varies (hours for large datasets)

**If FastTree crashes:**
1. Check sequence alignment: `python move_largest_numeric.py input.fasta sorted.fasta`
2. Look for unaligned sequences error message
3. Align sequences first or reduce dataset size

---

## 🔧 Configuration

### FastTree Stability vs Speed

**STABLE MODE (default, recommended):**
```python
make_trees_batch(stable_mode=True)  # Using -gtr -gamma flags
```

**FAST MODE (less stable):**
```python
make_trees_batch(stable_mode=False)  # Using -speediest flag
python run_trees.py --fast
```

### Memory Settings

Edit in `pipeline.py`:
```python
# For EcoSim (default 12GB)
run_ecosim_batch(..., memory_gb=24)

# For large datasets set timeout
make_trees_batch(fasttree_exe="fasttree")
```

---

## ❌ Troubleshooting

| Issue | Solution |
|-------|----------|
| `VeryFastTree/FastTree not found` | `brew install veryfasttree` (macOS) or `apt install veryfasttree` (Linux) |
| `Sequences not aligned` | Use `mafft your.fasta > aligned.fasta` |
| `FastTree crashes` | Ensure sequences have same length |
| `EcoSim jar not found` | Set `export ECOSIM_JAR=/path/to/ecosim.jar` and ensure Java is installed |
| `Tree rooting failed` | Check sequence IDs match exactly in FASTA and tree file |

---

## 📁 Output Files

After running pipeline, you'll get:
- `pipeline_temp_{species}/tree_rdy_fastas/` - Sorted sequences
- `pipeline_temp_{species}/trees_final/` - Unrooted trees
- `rerooted_trees/` - Rooted trees
- `rarefaction_fastas_{species}/` - Subsampled alignments
- `ecosim_output_{species}/` - EcoSim results (XML)

---

## 📖 More Details

See [IMPROVEMENTS.md](IMPROVEMENTS.md) for detailed explanation of improvements and advanced usage.
