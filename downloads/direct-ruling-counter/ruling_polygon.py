#!/usr/bin/env python3
"""Count direct rulings of convex lattice polygons in every possible genus.

This is the reader-facing command-line entry point.  It validates the polygon,
constructs a deterministic primitive multigeodesic diagram, and calls the
direct disk/annular ruling search.  It does not import a tropical counter or a
table of previously known answers.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from math import gcd
from pathlib import Path
from typing import Iterable, Sequence

from fig3_direct_ruling_counts import count_case


Point = tuple[int, int]


def _signed_double_area(vertices: Sequence[Point]) -> int:
    return sum(
        start[0] * finish[1] - start[1] * finish[0]
        for start, finish in zip(vertices, (*vertices[1:], vertices[0]))
    )


def _cross(origin: Point, left: Point, right: Point) -> int:
    return (
        (left[0] - origin[0]) * (right[1] - origin[1])
        - (left[1] - origin[1]) * (right[0] - origin[0])
    )


def _canonical_translation_key(vertices: Sequence[Point]) -> str:
    candidates: list[str] = []
    for sequence in (tuple(vertices), tuple(reversed(vertices))):
        for shift in range(len(sequence)):
            rotated = sequence[shift:] + sequence[:shift]
            origin_x, origin_y = rotated[0]
            candidates.append(
                ";".join(
                    f"{x - origin_x},{y - origin_y}" for x, y in rotated
                )
            )
    return min(candidates)


def _stable_id(vertices: Sequence[Point]) -> str:
    digest = hashlib.sha256(
        _canonical_translation_key(vertices).encode("ascii")
    ).hexdigest()[:12]
    return f"polygon-{digest}"


def _validate_vertices(raw_vertices: object) -> tuple[Point, ...]:
    if not isinstance(raw_vertices, list):
        raise ValueError("'vertices' must be a JSON list")
    vertices: list[Point] = []
    for point in raw_vertices:
        if (
            not isinstance(point, list)
            or len(point) != 2
            or any(isinstance(value, bool) or not isinstance(value, int)
                   for value in point)
        ):
            raise ValueError("every vertex must be a pair of JSON integers")
        vertices.append((point[0], point[1]))
    if len(vertices) > 1 and vertices[0] == vertices[-1]:
        vertices.pop()
    if len(vertices) < 3:
        raise ValueError("a polygon needs at least three vertices")
    if len(set(vertices)) != len(vertices):
        raise ValueError("polygon vertices must be distinct")

    turns = [
        _cross(
            vertices[index],
            vertices[(index + 1) % len(vertices)],
            vertices[(index + 2) % len(vertices)],
        )
        for index in range(len(vertices))
    ]
    if any(turn == 0 for turn in turns):
        raise ValueError("consecutive polygon edges must not be collinear")
    if not (all(turn > 0 for turn in turns) or all(turn < 0 for turn in turns)):
        raise ValueError(
            "vertices must be listed cyclically around a strictly convex polygon"
        )
    if _signed_double_area(vertices) < 0:
        vertices.reverse()
    return tuple(vertices)


def _polygon_invariants(vertices: Sequence[Point]) -> tuple[int, int, int]:
    area2 = _signed_double_area(vertices)
    boundary = sum(
        gcd(abs(finish[0] - start[0]), abs(finish[1] - start[1]))
        for start, finish in zip(vertices, (*vertices[1:], vertices[0]))
    )
    numerator = area2 - boundary + 2
    if area2 <= 0 or numerator < 0 or numerator % 2:
        raise ValueError("the input failed the lattice-polygon Pick check")
    return area2, boundary, numerator // 2


def _load_records(source: Path) -> list[dict[str, object]]:
    raw = json.loads(source.read_text(encoding="utf-8"))
    if isinstance(raw, dict) and "polygons" in raw:
        raw_records = raw["polygons"]
    elif isinstance(raw, dict) and "vertices" in raw:
        raw_records = [raw]
    elif isinstance(raw, list):
        raw_records = [{"vertices": raw}]
    else:
        raise ValueError(
            "input must be {'vertices': [...]}, {'polygons': [...]}, "
            "or a bare vertex list"
        )
    if not isinstance(raw_records, list) or not raw_records:
        raise ValueError("'polygons' must be a nonempty list")

    records: list[dict[str, object]] = []
    ids: set[str] = set()
    for raw_record in raw_records:
        if not isinstance(raw_record, dict):
            raise ValueError("each polygon record must be a JSON object")
        vertices = _validate_vertices(raw_record.get("vertices"))
        area2, boundary, interior = _polygon_invariants(vertices)
        identifier = str(raw_record.get("id") or _stable_id(vertices))
        if identifier in ids:
            raise ValueError(f"duplicate polygon id: {identifier}")
        ids.add(identifier)
        records.append(
            {
                "id": identifier,
                "vertices": [list(point) for point in vertices],
                "double_area": area2,
                "boundary_lattice_points": boundary,
                "interior_lattice_points": interior,
            }
        )
    return records


def _count_one(
    record: dict[str, object],
    genus: int,
    include_certificates: bool,
) -> tuple[str, int, dict[str, object]]:
    result = count_case(
        record,
        gallery_name="direct polygon input",
        target_genus=genus,
        incremental=True,
        include_certificates=include_certificates,
    )
    return str(record["id"]), genus, result


def count_records(
    records: Iterable[dict[str, object]],
    *,
    selected_genus: int | None,
    jobs: int,
    include_certificates: bool,
) -> dict[str, object]:
    records = list(records)
    tasks = [
        (record, genus)
        for record in records
        for genus in (
            [selected_genus]
            if selected_genus is not None
            else range(int(record["interior_lattice_points"]) + 1)
        )
    ]
    completed: dict[tuple[str, int], dict[str, object]] = {}
    if jobs == 1:
        for record, genus in tasks:
            identifier, actual_genus, result = _count_one(
                record, genus, include_certificates
            )
            completed[(identifier, actual_genus)] = result
    else:
        with ThreadPoolExecutor(max_workers=jobs) as executor:
            futures = {
                executor.submit(
                    _count_one, record, genus, include_certificates
                ): (record, genus)
                for record, genus in tasks
            }
            for future in as_completed(futures):
                identifier, actual_genus, result = future.result()
                completed[(identifier, actual_genus)] = result

    output_records: list[dict[str, object]] = []
    for record in records:
        identifier = str(record["id"])
        genera = []
        genus_values = (
            [selected_genus]
            if selected_genus is not None
            else range(int(record["interior_lattice_points"]) + 1)
        )
        for genus in genus_values:
            result = completed[(identifier, genus)]
            genus_record: dict[str, object] = {
                "genus": genus,
                "all_disk": result["all_disk"],
                "annular": result["annular"],
                "total": result["total"],
                "phase_cardinality": result["phase_cardinality"],
                "elapsed_seconds": result["elapsed_seconds"],
                "crossings": result["crossings"],
                "search_nodes": result["search_nodes"],
            }
            if include_certificates:
                genus_record["all_disk_masks"] = result["all_disk_masks"]
                genus_record["annular_results"] = result["annular_results"]
            genera.append(genus_record)
        first_genus = int(genera[0]["genus"])
        diagram = completed[(identifier, first_genus)]
        output_records.append({
            **record,
            "counting_shear": diagram["shear"],
            "vertical_direction_source": diagram[
                "vertical_direction_source"
            ],
            "vertical_direction_count": diagram[
                "vertical_direction_count"
            ],
            "sweep_covector_source": diagram["sweep_covector_source"],
            "sweep_covector_count": diagram["sweep_covector_count"],
            "sweep_selection": diagram["sweep_selection"],
            "transverse_covector_count": diagram[
                "transverse_covector_count"
            ],
            "type_b_epsilon": diagram["type_b_epsilon"],
            "genera": genera,
        })

    return {
        "schema": "direct-ruling-counts-v1",
        "method": (
            "primitive multigeodesic arrangement; rollback smoothing DFS; "
            "exact disk and unbounded annular audits"
        ),
        "uses_tropical_count": False,
        "polygons": output_records,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument(
        "--genus",
        type=int,
        help=(
            "compute only this genus; omit to compute every genus from 0 "
            "through the number of interior lattice points"
        ),
    )
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--certificates", action="store_true")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if arguments.genus is not None and arguments.genus < 0:
        parser.error("--genus must be nonnegative")
    if arguments.jobs < 1:
        parser.error("--jobs must be positive")

    try:
        records = _load_records(arguments.source)
        if arguments.genus is not None:
            too_large = [
                str(record["id"])
                for record in records
                if arguments.genus > int(record["interior_lattice_points"])
            ]
            if too_large:
                parser.error(
                    "requested genus exceeds the arithmetic genus for: "
                    + ", ".join(too_large)
                )
        payload = count_records(
            records,
            selected_genus=arguments.genus,
            jobs=arguments.jobs,
            include_certificates=arguments.certificates,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as error:
        parser.error(str(error))

    text = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    if arguments.output is None:
        print(text, end="")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(text, encoding="utf-8")
        print(arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
