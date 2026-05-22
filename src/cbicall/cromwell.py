import os
from pathlib import Path
from typing import Optional


def find_cromwell_jar_on_path() -> Optional[Path]:
    """Return the first cromwell*.jar found in PATH, preferring versioned jars."""
    candidates = []
    for entry in os.environ.get("PATH", "").split(os.pathsep):
        if not entry:
            continue
        directory = Path(entry)
        if not directory.is_dir():
            continue
        candidates.extend(sorted(directory.glob("cromwell-*.jar")))
        candidates.extend(sorted(directory.glob("cromwell.jar")))
        candidates.extend(sorted(directory.glob("cromwell*.jar")))
    seen = set()
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if candidate.is_file():
            return candidate
    return None
