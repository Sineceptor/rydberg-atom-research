"""Verify delivered files. New runs and notebook autosaves are not baselines."""
from pathlib import Path
import hashlib


def verify(root=None):
    root = Path(root) if root else Path(__file__).resolve().parent
    manifest = root / "SHA256SUMS.txt"
    failures = []
    lines = manifest.read_text().splitlines()
    for line in lines:
        expected, name = line.split("  ", 1)
        path = root / name
        actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else "MISSING"
        if actual != expected:
            failures.append(name)
    if failures:
        raise AssertionError("Delivered files differ: " + ", ".join(failures))
    print(f"Verified {len(lines)} delivered files. New runs are separate.")


if __name__ == "__main__":
    verify()
