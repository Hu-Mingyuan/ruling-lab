#!/usr/bin/env python3
"""Run the bundled polygon through every admissible genus."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent


def run_counter(source: Path, *arguments: str) -> dict[str, object]:
    completed = subprocess.run(
        [
            sys.executable,
            str(HERE / "ruling_polygon.py"),
            str(source),
            *arguments,
        ],
        cwd=HERE,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        sys.stderr.write(completed.stderr or completed.stdout)
        raise RuntimeError("direct ruling subprocess failed")
    return json.loads(completed.stdout)


def main() -> int:
    payload = run_counter(HERE / "example_polygon.json")

    assert payload["schema"] == "direct-ruling-counts-v1"
    assert payload["uses_tropical_count"] is False
    assert len(payload["polygons"]) == 1

    polygon = payload["polygons"][0]
    shear = int(polygon["counting_shear"])
    assert polygon["vertical_direction_source"] == [-shear, 1]
    assert polygon["vertical_direction_count"] == [0, 1]
    assert polygon["sweep_covector_source"] == [1, shear]
    assert polygon["sweep_covector_count"] == [1, 0]
    for start, finish in zip(
        polygon["vertices"],
        (*polygon["vertices"][1:], polygon["vertices"][0]),
    ):
        edge = (finish[0] - start[0], finish[1] - start[1])
        vertical = polygon["vertical_direction_source"]
        assert edge[0] * vertical[1] - edge[1] * vertical[0] != 0
    interior = int(polygon["interior_lattice_points"])
    genera = polygon["genera"]
    assert [record["genus"] for record in genera] == list(range(interior + 1))
    for record in genera:
        assert record["crossings"] > 0
        assert record["search_nodes"] > 0
        if record["phase_cardinality"] == "FINITE":
            assert record["total"] == record["all_disk"] + record["annular"]

    square_payload = run_counter(
        HERE / "example_vertical_polygon.json",
        "--genus",
        "0",
    )
    square = square_payload["polygons"][0]
    assert square["vertical_direction_source"] == [-1, 1]
    assert square["sweep_selection"] == "automatic-generic"
    assert [record["genus"] for record in square["genera"]] == [0]

    print(
        "direct-ruling smoke test passed: "
        f"{len(genera)} genera plus an automatically sheared square"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
