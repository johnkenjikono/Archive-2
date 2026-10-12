"""Load KEY=value lines from the repo's .env into os.environ. Variables already set in the shell win."""
import os
from pathlib import Path

ENV_FILE = Path(__file__).resolve().parent.parent / ".env"


def load_env(path=ENV_FILE):
    if not path.exists():
        return
    for line in path.read_text().splitlines():
        key, sep, value = line.partition("=")
        key = key.strip().removeprefix("export ").strip()
        value = value.strip().strip("'\"")
        # Skip blanks so an empty ANTHROPIC_API_KEY= doesn't shadow other credentials.
        if sep and key and value and not key.startswith("#"):
            os.environ.setdefault(key, value)


if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        f = Path(d) / ".env"
        f.write_text("# comment\nA_TEST=1\nexport B_TEST='two'\nC_TEST=\nPATH=nope\n")
        load_env(f)
    assert os.environ["A_TEST"] == "1" and os.environ["B_TEST"] == "two"
    assert "C_TEST" not in os.environ and os.environ["PATH"] != "nope"
    print("self-test passed")
