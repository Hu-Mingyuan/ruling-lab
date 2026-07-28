#!/usr/bin/env python3
"""Build and count a direct ruling diagram from polygon vertex coordinates."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from fractions import Fraction
from math import gcd
import json
from pathlib import Path
import random
from typing import Iterable, Sequence

from direct_ruling_count import count_multigeodesic_rulings
from multigeodesic_alternating_audit import TorusChamber, enumerate_torus_chambers
from verify_multigeodesic_rulings import make_parallel_geodesic_arrangement
from verify_triangle_rulings import Geodesic

PointI = tuple[int, int]
Direction = tuple[int, int]
ORIGIN_DENOMINATOR = 1009
COUNT_VERTICAL_DIRECTION: Direction = (0, 1)
COUNT_SWEEP_COVECTOR: Direction = (1, 0)


@dataclass(frozen=True)
class PolygonRulingDiagram:
    label: str
    vertices: tuple[PointI, ...]
    geodesics: tuple[Geodesic, ...]
    shear: int
    vertical_direction_source: Direction
    vertical_direction_count: Direction
    sweep_covector_source: Direction
    sweep_covector_count: Direction
    origin_residues: tuple[int, ...]
    origin_attempt: int
    crossings: int
    chambers: int
    alternating_chambers: tuple[TorusChamber, ...]


def _primitive_polygon_edges(vertices: Sequence[PointI]) -> tuple[Direction, ...]:
    directions: list[Direction] = []
    for start, finish in zip(vertices, vertices[1:] + vertices[:1]):
        dx = finish[0] - start[0]
        dy = finish[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        if length == 0:
            raise ValueError("consecutive polygon vertices coincide")
        directions.extend([(dx // length, dy // length)] * length)
    if tuple(map(sum, zip(*directions))) != (0, 0):
        raise AssertionError("expanded polygon edges are not balanced")
    return tuple(directions)


def _transverse_shear(directions: Sequence[Direction]) -> int:
    """Choose a count frame whose vertical direction misses every edge.

    For ``C_k(x,y)=(x+k*y,y)``, the count-frame vertical vector is
    ``(0,1)`` and its source-frame pullback is ``(-k,1)``.  Hence

    ``det((dx,dy),(-k,1)) = dx+k*dy``.
    """

    for magnitude in range(0, 1 + len(directions)):
        candidates = (0,) if magnitude == 0 else (magnitude, -magnitude)
        for shear in candidates:
            if all(dx + shear * dy != 0 for dx, dy in directions):
                return shear
    raise AssertionError("failed to find a transverse integral shear")


def build_polygon_ruling_diagram(
    label: str,
    vertices: Iterable[Sequence[int]],
    *,
    gallery_name: str = "lattice polygon",
    max_origin_attempts: int = 100,
) -> PolygonRulingDiagram:
    """Construct one deterministic generic diagram with an alternating face."""

    polygon = tuple((int(point[0]), int(point[1])) for point in vertices)
    primitive = _primitive_polygon_edges(polygon)
    shear = _transverse_shear(primitive)
    directions = tuple((dx + shear * dy, dy) for dx, dy in primitive)
    vertical_direction_source = (-shear, 1)
    sweep_covector_source = (1, shear)
    if any(
        dx * vertical_direction_source[1]
        - dy * vertical_direction_source[0] == 0
        for dx, dy in primitive
    ):
        raise AssertionError("source vertical direction is parallel to an edge")
    if any(
        dx * COUNT_VERTICAL_DIRECTION[1]
        - dy * COUNT_VERTICAL_DIRECTION[0] == 0
        for dx, dy in directions
    ):
        raise AssertionError("count vertical direction is parallel to an edge")
    if any(
        COUNT_SWEEP_COVECTOR[0] * dx
        + COUNT_SWEEP_COVECTOR[1] * dy == 0
        for dx, dy in directions
    ):
        raise AssertionError("count sweep covector vanishes on an edge")
    generator = random.Random(1000 + sum(ord(character) for character in label))

    last_error: Exception | None = None
    for attempt in range(max_origin_attempts):
        residues = tuple(
            generator.sample(range(ORIGIN_DENOMINATOR), len(directions))
        )
        lines = tuple(
            Geodesic(
                f"{label}{index}",
                direction,
                (Fraction(0), Fraction(residues[index], ORIGIN_DENOMINATOR)),
            )
            for index, direction in enumerate(directions)
        )
        try:
            arrangement = make_parallel_geodesic_arrangement(
                lines, name=f"{gallery_name}({label})"
            )
            chambers = enumerate_torus_chambers(lines, arrangement=arrangement)
        except (AssertionError, ValueError) as error:
            last_error = error
            continue
        alternating = tuple(chamber for chamber in chambers if chamber.alternating)
        if alternating:
            return PolygonRulingDiagram(
                label=label,
                vertices=polygon,
                geodesics=lines,
                shear=shear,
                vertical_direction_source=vertical_direction_source,
                vertical_direction_count=COUNT_VERTICAL_DIRECTION,
                sweep_covector_source=sweep_covector_source,
                sweep_covector_count=COUNT_SWEEP_COVECTOR,
                origin_residues=residues,
                origin_attempt=attempt,
                crossings=len(arrangement.vertices),
                chambers=len(chambers),
                alternating_chambers=alternating,
            )
    detail = f": {last_error}" if last_error is not None else ""
    raise RuntimeError(
        f"no generic alternating diagram for {gallery_name}({label}) after "
        f"{max_origin_attempts} attempts{detail}"
    )


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def count_case(
    record: dict[str, object],
    *,
    gallery_name: str = "lattice polygon",
    target_genus: int = 0,
    incremental: bool = True,
    include_certificates: bool = False,
) -> dict[str, object]:
    identifier = record.get("id", record.get("label"))
    if identifier is None:
        raise ValueError("each polygon record needs an 'id' or 'label'")
    label = str(identifier)
    diagram = build_polygon_ruling_diagram(  # type: ignore[arg-type]
        label,
        record["vertices"],
        gallery_name=gallery_name,
    )
    chamber = diagram.alternating_chambers[0]
    result = count_multigeodesic_rulings(
        diagram.geodesics,
        boundary_components=len(diagram.geodesics),
        target_genus=target_genus,
        basepoint=chamber.basepoint,
        incremental=incremental,
        name=f"{gallery_name}({label})",
    )
    audit = result.audit
    output = {
        "label": label,
        "gallery": gallery_name,
        "engine": "incremental" if incremental else "flat-reference",
        "target_genus": target_genus,
        "polygon_vertices_ccw": [list(point) for point in diagram.vertices],
        "lattice_perimeter": len(diagram.geodesics),
        "crossings": diagram.crossings,
        "shear": diagram.shear,
        "vertical_direction_source": list(
            diagram.vertical_direction_source
        ),
        "vertical_direction_count": list(
            diagram.vertical_direction_count
        ),
        "sweep_covector_source": list(diagram.sweep_covector_source),
        "sweep_covector_count": list(diagram.sweep_covector_count),
        "origin_denominator": ORIGIN_DENOMINATOR,
        "origin_residues": list(diagram.origin_residues),
        "origin_attempt": diagram.origin_attempt,
        "chambers": diagram.chambers,
        "alternating_chambers": len(diagram.alternating_chambers),
        "basepoint": [_fraction_text(value) for value in chamber.basepoint],
        "all_disk": result.all_disk,
        "annular": result.annular,
        "total": result.total,
        "phase_cardinality": result.phase_cardinality,
        "elapsed_seconds": result.elapsed_seconds,
        "search_nodes": audit.get("searchNodes"),
        "branch_attempts": audit.get("branchAttempts"),
        "completed_leaves": audit.get("completedLeaves"),
        "all_disk_leaves": audit.get("allDiskLeaves"),
        "annular_leaves": audit.get("annularLeaves"),
        "masks_requested": audit.get("masksRequested"),
        "masks_actually_traced": audit.get("masksActuallyTraced"),
        "exact_phase_solver": audit.get("diagnostics", {}).get("solver", {}),
        "all_disk_masks": audit.get("allDiskMasks", []),
    }
    if include_certificates:
        output["annular_results"] = audit.get("annularResults", [])
    return output


def _load_records(source: Path) -> list[dict[str, object]]:
    payload = json.loads(source.read_text(encoding="utf-8"))
    records = payload["polygons"]
    if not isinstance(records, list):
        raise ValueError("polygon source has no list-valued 'polygons' field")
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "labels",
        nargs="*",
        help="optional polygon ids; default: every record",
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--genus", type=int, default=0)
    parser.add_argument("--flat-reference", action="store_true")
    arguments = parser.parse_args()
    if arguments.jobs < 1:
        parser.error("--jobs must be positive")
    if arguments.genus < 0:
        parser.error("--genus must be nonnegative")
    records = _load_records(arguments.source)
    if arguments.labels:
        wanted = set(arguments.labels)
        records = [
            record
            for record in records
            if str(record.get("id", record.get("label"))) in wanted
        ]
        missing = wanted - {
            str(record.get("id", record.get("label"))) for record in records
        }
        if missing:
            parser.error("unknown labels: " + ", ".join(sorted(missing)))

    incremental = not arguments.flat_reference
    gallery_name = arguments.source.stem
    if arguments.jobs == 1:
        for record in records:
            print(
                json.dumps(
                    count_case(
                        record,
                        gallery_name=gallery_name,
                        target_genus=arguments.genus,
                        incremental=incremental,
                    ),
                    separators=(",", ":"),
                ),
                flush=True,
            )
        return 0

    with ThreadPoolExecutor(max_workers=arguments.jobs) as executor:
        futures = {
            executor.submit(
                count_case,
                record,
                gallery_name=gallery_name,
                target_genus=arguments.genus,
                incremental=incremental,
            ): record
            for record in records
        }
        for future in as_completed(futures):
            print(json.dumps(future.result(), separators=(",", ":")), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
