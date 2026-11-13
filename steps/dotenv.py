"""Load KEY=VALUE lines from the repo's .env into os.environ, if the file exists."""
import os
from pathlib import Path

ENV_FILE = Path(__file__).resolve().parent.parent / ".env"


def load_dotenv(path=ENV_FILE):
    """Set each variable from *path* unless it is already set in the shell."""
    # shortcut: no multi-line values or ${VAR} expansion; switch to python-dotenv if those are ever needed
    if not path.is_file():
        return
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.removeprefix("export ").split("=", 1)
        os.environ.setdefault(key.strip(), value.strip().strip("'\""))
