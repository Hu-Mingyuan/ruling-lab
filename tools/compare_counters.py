#!/usr/bin/env python3
"""Run the independent public counters, then compare all Laurent coefficients.

This validation utility is not imported by either counting algorithm.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parent.parent
RULING = ROOT / "downloads/direct-ruling-counter"
TROPICAL = ROOT / "downloads/tropical-tree-counter"


def source_hashes() -> dict[str, str]:
    return {
        path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for directory in (RULING, TROPICAL)
        for path in sorted(directory.rglob("*"))
        if path.suffix in {".py", ".js"} and not path.name.startswith("test_")
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="polygon JSON accepted by both counters")
    parser.add_argument("--output", type=Path, required=True, help="complete validation JSON")
    parser.add_argument("--jobs", type=int, default=1, help="parallel ruling genus searches")
    parser.add_argument("--seeds", type=int, nargs="+", default=[17471, 98189])
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    source = args.source.resolve()
    output = args.output.resolve()
    if not output.parent.is_dir():
        parser.error("the output directory must already exist")
    started = time.perf_counter()
    hashes = source_hashes()
    # Keep temporary files beside the requested report, where the caller has
    # chosen write access. Restricted environments may give subprocesses a
    # different view of the system temporary directory.
    with tempfile.TemporaryDirectory(prefix="ruling-tropical-check-", dir=output.parent) as temporary:
        directory = Path(temporary)
        ruling_path = directory / "ruling.json"
        completed = subprocess.run(
            [sys.executable, "-B", str(RULING / "ruling_polygon.py"), str(source),
             "--jobs", str(args.jobs)],
            check=True, capture_output=True, text=True,
        )
        ruling = json.loads(completed.stdout)
        ruling_path.write_text(completed.stdout, encoding="utf-8")
        print(f"Ruling enumeration complete: {len(ruling['polygons'])} polygons", flush=True)
        runs = []
        all_match = True
        for seed in args.seeds:
            tropical_path = directory / f"tropical-{seed}.json"
            completed = subprocess.run(
                [sys.executable, "-B", str(TROPICAL / "tropical_tree.py"), str(source),
                 "--seed", str(seed), "--compare-rulings", str(ruling_path),
                 "--output", str(tropical_path)],
                check=False,
            )
            if completed.returncode not in (0, 1):
                raise RuntimeError(f"tropical counter failed for seed {seed}")
            result = json.loads(tropical_path.read_text(encoding="utf-8"))
            all_match = all_match and completed.returncode == 0
            runs.append({"seed": seed, "exit_status": completed.returncode, "result": result})
        if source_hashes() != hashes:
            raise RuntimeError("counter source changed during the validation run")
        report = {
            "schema": "independent-counter-comparison-v1",
            "input": json.loads(source.read_text(encoding="utf-8")),
            "source_sha256": hashes,
            "ruling": ruling,
            "tropical_runs": runs,
            "all_match": all_match,
            "elapsed_seconds": time.perf_counter() - started,
        }
        output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"{len(ruling['polygons'])} polygons, {len(runs)} tropical seeds: "
              f"{'all match' if all_match else 'MISMATCH'}; report: {args.output}")
        return 0 if all_match else 1


if __name__ == "__main__":
    raise SystemExit(main())
