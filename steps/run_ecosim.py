import os
import shutil
import subprocess
from pathlib import Path


def _resolve_ecosim_jar(ecosim_jar=None):
    if ecosim_jar:
        return str(Path(ecosim_jar).expanduser().resolve())

    env_jar = os.environ.get("ECOSIM_JAR")
    if env_jar:
        return str(Path(env_jar).expanduser().resolve())

    default_jar = Path(__file__).resolve().parent.parent / "tools" / "ecosim.jar"
    if default_jar.exists():
        return str(default_jar.resolve())

    return None


def _iter_fasta_files(fasta_dir):
    base = Path(fasta_dir)
    if base.is_file():
        return [base]

    if not base.is_dir():
        raise FileNotFoundError(f"FASTA directory not found: {fasta_dir}")

    files = []
    for pattern in ("*.fasta", "*.fa", "*.fna"):
        files.extend(base.glob(pattern))
    return sorted(files)


def _resolve_tree_path(fasta_path, tree_dir=None, full_tree_path=None):
    if tree_dir:
        tree_root = Path(tree_dir)
        if not tree_root.is_dir():
            raise FileNotFoundError(f"Tree directory not found: {tree_dir}")

        stem = fasta_path.stem
        candidates = [
            tree_root / f"{stem}_tree.nwk",
            tree_root / f"{stem}.nwk",
            tree_root / f"{stem}_tree.newick",
            tree_root / f"{stem}.newick",
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate.resolve()

    if full_tree_path:
        tree_path = Path(full_tree_path)
        if not tree_path.exists():
            raise FileNotFoundError(f"Tree file not found: {full_tree_path}")
        return tree_path.resolve()

    raise FileNotFoundError(
        f"No tree found for {fasta_path.name}. Provide --tree-dir or --full-tree-path."
    )


def run_ecosim_batch(
    fasta_dir,
    tree_dir=None,
    output_dir="ecosim_results",
    ecosim_jar=None,
    ecosim_dir=None,
    memory_gb=12,
    full_tree_path=None,
    threads=None,
):
    """Run EcoSim.jar on a batch of FASTA files.

    Parameters
    ----------
    fasta_dir : str or Path
        Directory containing FASTA files, or a single FASTA file.
    tree_dir : str or Path, optional
        Directory containing per-FASTA tree files.
    output_dir : str or Path
        Directory where XML results are written.
    ecosim_jar : str, optional
        Path to ecosim.jar. Falls back to ECOSIM_JAR or a local ecosim.jar.
    ecosim_dir : str, optional
        Working directory for the EcoSim process.
    memory_gb : int
        Java heap size in gigabytes.
    full_tree_path : str or Path, optional
        Tree file to reuse for every FASTA if no per-file tree exists.
    threads : int, optional
        EcoSim worker threads (EcoSim defaults to all CPUs).

    Returns the number of FASTA files EcoSim succeeded on.
    """

    if not shutil.which("java"):
        raise FileNotFoundError(
            "Java not found. Install a JRE/JDK and make sure `java` is on PATH."
        )

    ecosim_jar = _resolve_ecosim_jar(ecosim_jar)
    if not ecosim_jar or not Path(ecosim_jar).exists():
        raise FileNotFoundError(
            "EcoSim jar not found. Set ECOSIM_JAR or pass ecosim_jar=/path/to/ecosim.jar"
        )

    output_path = Path(output_dir).expanduser().resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    fasta_files = _iter_fasta_files(fasta_dir)
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in {fasta_dir}")

    if tree_dir is None and full_tree_path is None:
        raise ValueError("Either tree_dir or full_tree_path must be provided.")

    if ecosim_dir is not None:
        ecosim_dir = str(Path(ecosim_dir).expanduser().resolve())

    succeeded = 0
    failed = 0

    print(f"📁 Found {len(fasta_files)} FASTA file(s) in {Path(fasta_dir)}")
    print(f"☕ Using EcoSim jar: {ecosim_jar}")
    if ecosim_dir:
        print(f"📂 Working directory: {ecosim_dir}")

    for fasta_file in fasta_files:
        fasta_path = fasta_file.resolve()
        try:
            tree_path = _resolve_tree_path(
                fasta_path,
                tree_dir=tree_dir,
                full_tree_path=full_tree_path,
            )
        except FileNotFoundError as e:
            print(f"⚠️  {e}")
            failed += 1
            continue

        xml_output = output_path / f"{fasta_path.stem}_results.xml"

        print(f"🌱 Running EcoSim for {fasta_path.name}...")
        cmd = [
            "java",
            f"-Xmx{memory_gb}G",
            "-jar",
            ecosim_jar,
            f"-s={fasta_path}",
            f"-p={tree_path}",
            f"-o={xml_output}",
            "-r",
            "-n",
        ]
        if threads:
            cmd.append(f"-t={threads}")

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                cwd=ecosim_dir,
            )

            if result.stdout.strip():
                print(result.stdout.strip())
            print(f"✅ EcoSim output saved: {xml_output}")
            succeeded += 1
        except subprocess.CalledProcessError as e:
            failed += 1
            stderr = (e.stderr or "").strip()
            stdout = (e.stdout or "").strip()
            if stdout:
                print(stdout)
            if stderr:
                print(f"❌ EcoSim error for {fasta_path.name}: {stderr[:500]}")
            else:
                print(f"❌ EcoSim error for {fasta_path.name}: {e}")

    print("\n" + "=" * 50)
    print("EcoSim analysis complete!")
    print(f"  ✅ Succeeded: {succeeded}")
    print(f"  ❌ Failed: {failed}")
    print(f"  📊 Output directory: {output_path}")
    print("=" * 50)
    return succeeded


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Run EcoSim.jar on a batch of FASTA files"
    )
    parser.add_argument("fasta_dir", help="Directory containing FASTA files")
    parser.add_argument(
        "--tree-dir",
        default=None,
        help="Directory containing per-FASTA tree files (optional)",
    )
    parser.add_argument(
        "--full-tree-path",
        default=None,
        help="Single tree file to use for every FASTA if no per-file tree exists",
    )
    parser.add_argument(
        "--output-dir",
        default="ecosim_results",
        help="Directory to write EcoSim XML results (default: ecosim_results)",
    )
    parser.add_argument(
        "--ecosim-jar",
        default=None,
        help="Path to ecosim.jar (defaults to ECOSIM_JAR or ./ecosim.jar)",
    )
    parser.add_argument(
        "--ecosim-dir",
        default=None,
        help="Working directory for EcoSim (optional)",
    )
    parser.add_argument(
        "--memory-gb",
        type=int,
        default=12,
        help="Java heap size in gigabytes (default: 12)",
    )

    args = parser.parse_args()

    run_ecosim_batch(
        fasta_dir=args.fasta_dir,
        tree_dir=args.tree_dir,
        output_dir=args.output_dir,
        ecosim_jar=args.ecosim_jar,
        ecosim_dir=args.ecosim_dir,
        memory_gb=args.memory_gb,
        full_tree_path=args.full_tree_path,
    )
