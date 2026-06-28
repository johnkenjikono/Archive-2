# Pipeline Improvements Guide

## Critical Issues Fixed

### 1. **FastTree Crashes** ✅
**Problem:** FastTree was crashing with unaligned sequences or aggressive flags  
**Solutions Added:**
- **Sequence alignment validation** - Checked before tree building
- **Stable mode (default)** - Conservative flags (`-gtr -gamma`) for stability, with `-speediest` option for speed
- **Error capture** - Now shows VeryFastTree/FastTree error messages (stderr)
- **Retry logic** - Automatically retries failed trees up to 2 times
- **Timeout protection** - 10-minute timeout per file to prevent hangs
- **Progress reporting** - Better status messages and summary at end
- **VeryFastTree support** - Optimized for faster tree building on modern hardware

### 2. **Cross-Platform Path Issues** ✅
**Problem:** Hardcoded Windows paths (`C:\ecosim\`, `C:/ecosim/`) don't work on macOS/Linux  
**Solutions:**
- **Tree builder auto-detection** - Searches PATH for VeryFastTree, FastTree, or FastTreeMP
- **Environment variable support** - Set `ECOSIM_JAR` and `ECOSIM_DIR` for your system
- **Better error messages** - Tells you exactly how to set paths if executables not found

### 3. **Missing Data Validation** ✅
**Problem:** No checks for misaligned sequences before tree building  
**Solutions:**
- [move_largest_numeric.py](move_largest_numeric.py): Now validates all sequences same length
- [run_trees.py](run_trees.py): Validates input alignments before calling FastTree

---

## How to Use the Improved Pipeline

### Installation

```bash
cd /Users/carlofedolfi/Downloads/move_largest_numeric
source venv/bin/activate.fish  # Activate virtual environment

# VeryFastTree should already be installed (faster than FastTree!)
# If not: brew install veryfasttree
```

### Configuration (if using EcoSim)

For macOS/Linux, set environment variables:

```bash
export ECOSIM_JAR=/path/to/ecosim.jar
export ECOSIM_DIR=/path/to/ecosim/directory
```

### Running the Pipeline

**Basic usage (recommended - stable mode):**
```bash
python pipeline.py your_input.fasta
```

**For large diverse datasets (very stable but slower):**
- Edit [run_trees.py](run_trees.py) line with `make_trees_batch()` call
- Change to: `make_trees_batch(stable_mode=True)`

**For smaller, less diverse datasets (faster):**
```bash
python run_trees.py --fast
```

### Running Individual Steps

**Step 1: Sort sequences by numeric value**
```bash
python move_largest_numeric.py input.fasta output.fasta
```

**Step 2: Build trees**
```bash
python run_trees.py
```

**Step 3: Reroot trees**
```bash
python reroot_tree.py sorted.fasta tree.nwk output.nwk
```

---

## Performance Tips for Large Datasets

1. **Verify alignment first** - Before running pipeline, ensure sequences are aligned
   ```bash
   # Check with:
   python -c "from Bio import SeqIO; recs = list(SeqIO.parse('your.fasta', 'fasta')); lengths = set(len(r.seq) for r in recs); print(f'Unique lengths: {lengths}')"
   ```

2. **Start with stable mode** - If getting crashes, use default stable mode
   ```python
   make_trees_batch(stable_mode=True)  # Conservative flags
   ```

3. **Memory for EcoSim** - Set memory if running out:
   ```python
   run_ecosim_batch(memory_gb=24)  # Increase as needed
   ```

4. **Skip completed steps** - Resume from specific step:
   ```bash
   python pipeline.py input.fasta --start-step 3
   ```

---

## Troubleshooting

### FastTree still crashing?
1. Check that sequences are aligned (all same length)
2. Try stable mode: `make_trees_batch(stable_mode=True)`
3. Check FastTree is installed: `which fasttree`

### EcoSim path not found?
```bash
export ECOSIM_JAR=/path/to/ecosim.jar
export ECOSIM_DIR=/path/to/ecosim
```

### Sequence length mismatch error?
Use an alignment tool first (MAFFT, muscle, clustalw):
```bash
mafft your_sequences.fasta > aligned.fasta
```

### Tree rooting failed?
Make sure sequence IDs match exactly between FASTA and tree file (no extra spaces/characters).

---

## Key Changes Summary

| File | Changes |
|------|---------|
| [run_trees.py](run_trees.py) | Auto-detect FastTree, validate alignment, capture errors, retry logic, stable/fast modes |
| [pipeline.py](pipeline.py) | Environment variable support for EcoSim paths, better error messages |
| [run_ecosim.py](run_ecosim.py) | Environment variable support for paths |
| [move_largest_numeric.py](move_largest_numeric.py) | Alignment validation before processing |
| [reroot_tree.py](reroot_tree.py) | Better error messages, sequence name matching hints |
| [requirements.txt](requirements.txt) | Created with pinned versions |

---

## Next Steps

1. **Test with a small dataset first** to ensure everything works on your system
2. **Set FastTree in PATH** or verify it's installed: `brew install fasttree` (macOS)
3. **Configure EcoSim paths** if you plan to use that step
4. **Monitor the first few runs** to catch any platform-specific issues
