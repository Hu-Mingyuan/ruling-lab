#!/usr/bin/env python3
"""Run the exact direct disk/annular ruling search."""

from __future__ import annotations

import json
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from multigeodesic_annular_input import exact_annular_payload
from verify_triangle_rulings import Geodesic, Point


@dataclass(frozen=True)
class DirectRulingCount:
    all_disk: int
    annular: int | str
    total: int | str
    phase_cardinality: str
    elapsed_seconds: float
    audit: dict[str, object]


def _node_executable(explicit: str | Path | None = None) -> Path:
    if explicit is not None:
        candidate = Path(explicit).expanduser().resolve()
        if candidate.is_file():
            return candidate
        raise FileNotFoundError(f"Node.js executable does not exist: {candidate}")
    on_path = shutil.which("node") or shutil.which("node.exe")
    if on_path:
        return Path(on_path).resolve()
    bundled = (
        Path.home()
        / ".cache"
        / "codex-runtimes"
        / "codex-primary-runtime"
        / "dependencies"
        / "node"
        / "bin"
        / "node.exe"
    )
    if bundled.is_file():
        return bundled
    raise FileNotFoundError(
        "Node.js was not found; install Node.js or pass --node /path/to/node"
    )


def count_multigeodesic_rulings(
    geodesics: Iterable[Geodesic],
    *,
    boundary_components: int,
    target_genus: int,
    basepoint: Point,
    tau: tuple[int, int] | None = None,
    sigma: tuple[int, int] | None = None,
    node_executable: str | Path | None = None,
    incremental: bool = True,
    name: str = "primitive multigeodesic arrangement",
) -> DirectRulingCount:
    """Exhaust the direct diagram, including exact no-cutoff annular phases."""

    payload = exact_annular_payload(
        geodesics,
        boundary_components=boundary_components,
        target_genus=target_genus,
        basepoint=basepoint,
        masks=None,
        tau=tau,
        sigma=sigma,
        name=name,
    )
    script_name = (
        "generic_incremental_ruling_search.js"
        if incremental
        else "generic_annular_exact_search.js"
    )
    script = Path(__file__).with_name("tmp") / script_name
    started = time.perf_counter()
    completed = subprocess.run(
        [str(_node_executable(node_executable)), str(script)],
        input=json.dumps(payload),
        text=True,
        capture_output=True,
        check=False,
        cwd=Path(__file__).parent,
    )
    elapsed = time.perf_counter() - started
    if completed.returncode != 0:
        raise RuntimeError(
            "direct ruling engine failed: "
            + (completed.stderr.strip() or f"exit code {completed.returncode}")
        )
    audit = json.loads(completed.stdout)
    audit["sweepCovector"] = payload["tau"]
    audit["transverseCovector"] = payload["sigma"]
    audit["verticalDirection"] = payload["verticalDirection"]
    audit["sweepSelection"] = payload["sweepSelection"]
    audit["typeBEpsilon"] = payload["typeBEpsilon"]
    phase = str(audit["phaseCardinality"])
    return DirectRulingCount(
        all_disk=int(audit["exactAllDiskRulingCount"]),
        annular=(
            int(audit["exactAnnularRulingCount"])
            if phase == "FINITE"
            else str(audit["exactAnnularRulingCount"])
        ),
        total=(
            int(audit["exactTotalRulingCount"])
            if phase == "FINITE"
            else str(audit["exactTotalRulingCount"])
        ),
        phase_cardinality=phase,
        elapsed_seconds=elapsed,
        audit=audit,
    )
