#!/usr/bin/env python3
"""Regression test for simultaneous crossings in the download engine."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "downloads" / "direct-ruling-counter"
sys.path.insert(0, str(ENGINE))

from fig3_direct_ruling_counts import build_polygon_ruling_diagram  # noqa: E402
from multigeodesic_annular_input import exact_annular_payload  # noqa: E402
from verify_multigeodesic_rulings import (  # noqa: E402
    make_parallel_geodesic_arrangement,
)


def main() -> int:
    label = "polygon-c848bca89e"
    vertices = ((0, 2), (2, 0), (4, 1), (2, 2))
    diagram = build_polygon_ruling_diagram(
        label,
        vertices,
        gallery_name="generic-sweep-regression",
    )
    arrangement = make_parallel_geodesic_arrangement(
        diagram.geodesics,
        name=label,
    )
    x_levels = [vertex.position[0] for vertex in arrangement.vertices]
    assert len(x_levels) != len(set(x_levels)), "fixture must have tied x-levels"

    payload = exact_annular_payload(
        diagram.geodesics,
        boundary_components=len(diagram.geodesics),
        target_genus=0,
        basepoint=diagram.alternating_chambers[0].basepoint,
        name=label,
    )
    tau = tuple(payload["tau"])
    assert tau != (1, 0)
    assert payload["sweepSelection"] == "automatic-generic"

    raw = payload["arrangement"]
    period = int(raw["L"])
    levels = [
        (tau[0] * int(vertex["x"]) + tau[1] * int(vertex["y"])) % period
        for vertex in raw["vertices"]
    ]
    assert len(levels) == len(set(levels))
    assert raw["crossingOrder"] == list(arrangement.crossing_order)

    ordered = sorted(levels)
    gaps = [
        finish - start
        for start, finish in zip(
            ordered,
            (*ordered[1:], ordered[0] + period),
        )
    ]
    epsilon = Fraction(str(payload["typeBEpsilon"]))
    assert 0 < epsilon <= Fraction(min(gaps), 4)

    try:
        exact_annular_payload(
            diagram.geodesics,
            boundary_components=len(diagram.geodesics),
            target_genus=0,
            basepoint=diagram.alternating_chambers[0].basepoint,
            tau=(1, 0),
            sigma=(0, 1),
            name=label,
        )
    except ValueError as error:
        assert "same sweep level" in str(error)
    else:
        raise AssertionError("an explicitly nongeneric sweep must be rejected")

    print(f"generic sweep regression passed: tau={tau}, epsilon={epsilon}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
