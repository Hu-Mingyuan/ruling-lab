#!/usr/bin/env python3
"""Build exact torus arrangements from primitive geodesics."""

from __future__ import annotations

from fractions import Fraction
from functools import cmp_to_key
from itertools import combinations
from math import gcd
from typing import Iterable

import ruling_smoothing as rs
from verify_triangle_rulings import (
    Geodesic,
    Point,
    Vector,
    _det,
    _intersection_parameters,
    _ray_compare,
)

COUNT_VERTICAL_DIRECTION: Vector = (0, 1)
COUNT_SWEEP_COVECTOR: Vector = (1, 0)


def _sweep_value(direction: Vector) -> int:
    return (
        COUNT_SWEEP_COVECTOR[0] * direction[0]
        + COUNT_SWEEP_COVECTOR[1] * direction[1]
    )


def _make_geodesic_arrangement(
    geodesics: Iterable[Geodesic],
    *,
    name: str = "primitive multigeodesic arrangement",
    allow_parallel: bool,
) -> rs.TorusArrangement:
    """Build an exact torus half-edge arrangement of ``geodesics``.

    Directions are the oriented boundary classes.  They must be primitive,
    balanced, and have nonzero first coordinate so the chosen vertical sweep
    is transverse to every strand.  When ``allow_parallel`` is true, parallel
    primitive geodesics are allowed provided their images in the torus are
    disjoint.  (Coincident parallel geodesics do not define a four-valent
    generic arrangement.)
    """

    lines = tuple(geodesics)
    if len(lines) < 2:
        raise ValueError("at least two geodesics are required")
    if len({line.name for line in lines}) != len(lines):
        raise ValueError("geodesic names must be distinct")
    for line in lines:
        if line.direction == (0, 0):
            raise ValueError("a geodesic direction must be nonzero")
        if gcd(abs(line.direction[0]), abs(line.direction[1])) != 1:
            raise ValueError("every geodesic direction must be primitive")
        if _det(line.direction, COUNT_VERTICAL_DIRECTION) == 0:
            raise ValueError(
                "a geodesic is parallel to the count-frame vertical "
                "direction v=(0,1)"
            )
        if _sweep_value(line.direction) == 0:
            raise AssertionError(
                "the sweep covector must be nonzero on every geodesic"
            )
    if tuple(sum(line.direction[j] for line in lines) for j in range(2)) != (0, 0):
        raise ValueError("the oriented boundary classes are not balanced")

    records: list[dict[str, object]] = []
    parameters: list[dict[int, Fraction]] = [dict() for _ in lines]
    positions_seen: set[Point] = set()
    expected_crossings = 0
    for left_index, right_index in combinations(range(len(lines)), 2):
        left = lines[left_index]
        right = lines[right_index]
        determinant = int(_det(left.direction, right.direction))
        if determinant == 0:
            if not allow_parallel:
                raise ValueError("parallel geodesics are not supported")
            transverse_offset = _det(
                left.direction,
                (
                    right.origin[0] - left.origin[0],
                    right.origin[1] - left.origin[1],
                ),
            )
            if transverse_offset.denominator == 1:
                raise ValueError(
                    "parallel geodesics have coincident torus images: "
                    f"{left.name!r} and {right.name!r}"
                )
            continue
        expected_crossings += abs(determinant)
        on_left = _intersection_parameters(left, right)
        on_right = _intersection_parameters(right, left)
        if set(on_left) != set(on_right):
            raise AssertionError("the two parameterizations disagree")
        for point in sorted(on_left):
            if point in positions_seen:
                raise ValueError(
                    "the chosen origins create a nongeneric triple crossing "
                    f"at {point}"
                )
            positions_seen.add(point)
            vertex = len(records)
            records.append(
                {
                    "point": point,
                    "families": (left_index, right_index),
                    "label": f"{left.name}{right.name}_{point}",
                }
            )
            parameters[left_index][vertex] = on_left[point]
            parameters[right_index][vertex] = on_right[point]

    if len(records) != expected_crossings:
        raise AssertionError("wrong number of pairwise crossings")

    name_to_family = {line.name: family for family, line in enumerate(lines)}
    port_of: dict[tuple[int, int, int], int] = {}
    vertices: list[rs.TorusVertex] = []
    for vertex, record in enumerate(records):
        first, second = record["families"]  # type: ignore[misc]
        rays: list[tuple[str, int, Vector]] = []
        for family in (first, second):
            direction = lines[family].direction
            rays.append((lines[family].name, 1, direction))
            rays.append(
                (lines[family].name, -1, (-direction[0], -direction[1]))
            )
        rays.sort(key=cmp_to_key(_ray_compare))
        ports = tuple(4 * vertex + index for index in range(4))
        for port, (line_name, sign, _direction) in zip(ports, rays):
            port_of[(vertex, name_to_family[line_name], sign)] = port
        future_ports = tuple(
            port_of[
                (
                    vertex,
                    family,
                    1 if _sweep_value(lines[family].direction) > 0 else -1,
                )
            ]
            for family in (first, second)
        )
        vertices.append(
            rs.TorusVertex(
                index=vertex,
                position=record["point"],  # type: ignore[arg-type]
                ports=ports,
                future_ports=future_ports,  # type: ignore[arg-type]
                label=str(record["label"]),
            )
        )

    port_count = 4 * len(vertices)
    twins = [-1] * port_count
    shifts: list[tuple[int, int] | None] = [None] * port_count
    strands = [""] * port_count
    signs = [0] * port_count
    for vertex, record in enumerate(records):
        for family in record["families"]:  # type: ignore[union-attr]
            signs[port_of[(vertex, family, 1)]] = 1
            signs[port_of[(vertex, family, -1)]] = -1

    for family, line in enumerate(lines):
        occurrences = sorted(
            (parameter, vertex)
            for vertex, parameter in parameters[family].items()
        )
        expected = sum(
            abs(int(_det(line.direction, other.direction)))
            for index, other in enumerate(lines)
            if index != family
        )
        if len(occurrences) != expected:
            raise AssertionError(f"{line.name} has the wrong crossing count")
        if not occurrences:
            raise ValueError(
                f"geodesic {line.name!r} has no crossings; vertex-free "
                "components are not represented by TorusArrangement"
            )
        for index, (parameter, origin_vertex) in enumerate(occurrences):
            next_parameter, target_vertex = occurrences[(index + 1) % len(occurrences)]
            delta = next_parameter - parameter
            if delta <= 0:
                delta += 1
            origin = records[origin_vertex]["point"]
            target = records[target_vertex]["point"]
            lifted_target = (
                origin[0] + delta * line.direction[0],  # type: ignore[index]
                origin[1] + delta * line.direction[1],  # type: ignore[index]
            )
            raw_shift = (
                lifted_target[0] - target[0],  # type: ignore[index]
                lifted_target[1] - target[1],  # type: ignore[index]
            )
            if any(value.denominator != 1 for value in raw_shift):
                raise AssertionError("a geodesic edge has fractional deck shift")
            shift = int(raw_shift[0]), int(raw_shift[1])
            outgoing = port_of[(origin_vertex, family, 1)]
            incoming = port_of[(target_vertex, family, -1)]
            if twins[outgoing] >= 0 or twins[incoming] >= 0:
                raise AssertionError("a port belongs to two geodesic edges")
            twins[outgoing] = incoming
            twins[incoming] = outgoing
            shifts[outgoing] = shift
            shifts[incoming] = (-shift[0], -shift[1])
            strands[outgoing] = strands[incoming] = line.name

    if any(twin < 0 for twin in twins) or any(shift is None for shift in shifts):
        raise AssertionError("not every port was joined")
    halfedges = tuple(
        rs.TorusHalfEdge(
            index=port,
            origin=port // 4,
            twin=twins[port],
            deck_shift=shifts[port],  # type: ignore[arg-type]
            coorientation_sign=signs[port],
            strand=strands[port],
        )
        for port in range(port_count)
    )
    crossing_order = tuple(
        sorted(
            range(len(vertices)),
            key=lambda vertex: (
                vertices[vertex].position[0],
                vertices[vertex].position[1],
                vertex,
            ),
        )
    )
    arrangement = rs.TorusArrangement(
        periods=(Fraction(1), Fraction(1)),
        vertices=tuple(vertices),
        halfedges=halfedges,
        crossing_order=crossing_order,
        name=name,
    )
    errors = arrangement.validation_errors()
    if errors:
        raise AssertionError("invalid arrangement: " + "; ".join(errors))
    return arrangement


def make_geodesic_arrangement(
    geodesics: Iterable[Geodesic],
    *,
    name: str = "primitive multigeodesic arrangement",
) -> rs.TorusArrangement:
    """Build the exact arrangement of pairwise nonparallel geodesics.

    This is the original strict API.  Use
    :func:`make_parallel_geodesic_arrangement` for several translated copies
    of the same primitive boundary direction.
    """

    return _make_geodesic_arrangement(
        geodesics,
        name=name,
        allow_parallel=False,
    )


def make_parallel_geodesic_arrangement(
    geodesics: Iterable[Geodesic],
    *,
    name: str = "primitive multigeodesic arrangement with parallel families",
) -> rs.TorusArrangement:
    """Build an exact generic arrangement allowing parallel families.

    Every direction is still required to be primitive and transverse to the
    count-frame vertical direction ``v=(0,1)`` (equivalently, the sweep
    covector ``tau=(1,0)`` is nonzero on it).  Parallel geodesics must be
    distinct translates on the torus; nonparallel pairs must have no triple
    (or higher) intersection.  The total oriented homology class must vanish.
    """

    return _make_geodesic_arrangement(
        geodesics,
        name=name,
        allow_parallel=True,
    )
