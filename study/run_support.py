"""Keep delivered evidence separate from new calculations.

Each managed run is built in its own staging folder. Only successful runs are
renamed to a visible final folder. Failed work remains clearly marked .failed.
No latest pointer, source CSV or delivered figure is replaced.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import platform
import uuid

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run_managed(worker, label, runs_root=None):
    runs_root = Path(runs_root) if runs_root is not None else ROOT / "runs"
    runs_root.mkdir(parents=True, exist_ok=True)
    ident = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:8]
    final = runs_root / (ident + "-" + label)
    staging = runs_root / ("." + final.name + ".incomplete")
    staging.mkdir()
    metadata = {"status": "running", "label": label, "started_utc": ident[:16],
                "python": platform.python_version(),
                "source_sha256": {p.name: digest(p) for p in sorted(ROOT.glob("*.py"))},
                "scope": "Recalculation only. Guides and delivered baseline are not rewritten."}
    (staging / "run_metadata.json").write_text(json.dumps(metadata, indent=2))
    try:
        result = worker(staging)
        metadata["status"] = "complete"
        metadata["completed_utc"] = datetime.now(timezone.utc).isoformat()
        (staging / "run_metadata.json").write_text(json.dumps(metadata, indent=2))
        lines = [f"{digest(p)}  {p.relative_to(staging).as_posix()}" for p in sorted(staging.rglob("*")) if p.is_file()]
        (staging / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n")
        os.rename(staging, final)
        print(f"Completed run saved to: {final}", flush=True)
        return result
    except BaseException as exc:
        metadata.update(status="failed", error_type=type(exc).__name__, error=str(exc))
        (staging / "run_metadata.json").write_text(json.dumps(metadata, indent=2))
        failed = runs_root / (final.name + ".failed")
        os.rename(staging, failed)
        print(f"Run failed; diagnostic files kept at: {failed}", flush=True)
        raise
