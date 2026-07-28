#!/usr/bin/env python3
"""Construct exact alternating chambers for primitive torus geodesics."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from math import floor
from typing import Iterable

from verify_triangle_rulings import Geodesic, Point, Vector
import verify_multigeodesic_rulings as vmr

Descriptor = tuple[int, str, int]


def determinant(left: Vector, right: Vector) -> int | Fraction:
    return left[0] * right[1] - left[1] * right[0]


def line_value(line: Geodesic, point: Point) -> Fraction:
    """Return ``det(direction, point-origin)`` independently of the target."""

    return Fraction(
        determinant(
            line.direction,
            (point[0] - line.origin[0], point[1] - line.origin[1]),
        )
    )


def level_intersection(
    left: Geodesic,
    left_level: Fraction,
    right: Geodesic,
    right_level: Fraction,
) -> Point:
    """Solve ``det(v_i,x-o_i)=level_i`` by Cramer's rule."""

    det = determinant(left.direction, right.direction)
    if det == 0:
        raise ValueError("parallel lifts do not have a unique intersection")
    left_constant = left_level + determinant(left.direction, left.origin)
    right_constant = right_level + determinant(right.direction, right.origin)
    return (
        Fraction(
            left_constant * right.direction[0]
            - left.direction[0] * right_constant,
            det,
        ),
        Fraction(
            left_constant * right.direction[1]
            - left.direction[1] * right_constant,
            det,
        ),
    )


def strip_bounds(
    lines: tuple[Geodesic, ...], point: Point
) -> tuple[tuple[int, int], ...]:
    result: list[tuple[int, int]] = []
    for line in lines:
        value = line_value(line, point)
        if value.denominator == 1:
            raise ValueError("point lies on the lifted arrangement")
        lower = floor(value)
        result.append((lower, lower + 1))
    return tuple(result)


def _clip_halfplane(
    polygon: tuple[Point, ...],
    line: Geodesic,
    level: int,
    lower_side: bool,
) -> tuple[Point, ...]:
    """Exact Sutherland--Hodgman clipping by one lifted-line half-plane."""

    def height(point: Point) -> Fraction:
        value = line_value(line, point) - level
        return value if lower_side else -value

    output: list[Point] = []
    for current, following in zip(polygon, polygon[1:] + polygon[:1]):
        current_height = height(current)
        following_height = height(following)
        current_inside = current_height >= 0
        following_inside = following_height >= 0
        if current_inside:
            output.append(current)
        if current_inside != following_inside:
            parameter = current_height / (current_height - following_height)
            output.append(
                (
                    current[0] + parameter * (following[0] - current[0]),
                    current[1] + parameter * (following[1] - current[1]),
                )
            )

    cleaned: list[Point] = []
    for point in output:
        if not cleaned or point != cleaned[-1]:
            cleaned.append(point)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return tuple(cleaned)


def _twice_area(polygon: tuple[Point, ...]) -> Fraction:
    return sum(
        left[0] * right[1] - left[1] * right[0]
        for left, right in zip(polygon, polygon[1:] + polygon[:1])
    )


@dataclass(frozen=True)
class IndependentCell:
    vertices: tuple[Point, ...]
    sides: tuple[Descriptor, ...]
    orientation_signs: tuple[int, ...]

    @property
    def alternating(self) -> bool:
        signs = self.orientation_signs
        return (
            len(signs) >= 4
            and len(signs) % 2 == 0
            and all(signs[index] == -signs[(index + 1) % len(signs)]
                    for index in range(len(signs)))
        )


def independent_lifted_cell(
    geodesics: Iterable[Geodesic], basepoint: Point
) -> IndependentCell:
    """Construct the lifted chamber by exact polygon clipping."""

    lines = tuple(geodesics)
    bounds = strip_bounds(lines, basepoint)
    descriptors = tuple(
        (family, side, bounds[family][0 if side == "lower" else 1])
        for family in range(len(lines))
        for side in ("lower", "upper")
    )

    intersections = tuple(
        level_intersection(
            lines[left[0]], Fraction(left[2]),
            lines[right[0]], Fraction(right[2]),
        )
        for left, right in combinations(descriptors, 2)
        if determinant(
            lines[left[0]].direction,
            lines[right[0]].direction,
        ) != 0
    )
    radius = 1 + max(
        abs(coordinate)
        for point in intersections
        for coordinate in point
    )
    polygon: tuple[Point, ...] = (
        (-radius, -radius),
        (radius, -radius),
        (radius, radius),
        (-radius, radius),
    )
    for family, line in enumerate(lines):
        lower, upper = bounds[family]
        polygon = _clip_halfplane(polygon, line, lower, True)
        polygon = _clip_halfplane(polygon, line, upper, False)
    if len(polygon) < 3 or _twice_area(polygon) <= 0:
        raise AssertionError("strip intersection is not a CCW bounded polygon")

    sides: list[Descriptor] = []
    for vertex, following in zip(polygon, polygon[1:] + polygon[:1]):
        active = []
        for descriptor in descriptors:
            family, _side, level = descriptor
            if (
                line_value(lines[family], vertex) == level
                and line_value(lines[family], following) == level
            ):
                active.append(descriptor)
        if len(active) != 1:
            raise AssertionError(
                f"polygon edge has {len(active)} supporting lifted lines"
            )
        sides.append(active[0])

    # Along a lower boundary the strip is to the left of the original
    # direction, so its CCW boundary orientation agrees with that direction.
    signs = tuple(1 if side == "lower" else -1 for _, side, _ in sides)
    return IndependentCell(polygon, tuple(sides), signs)


def _local_sector_basepoints(
    lines: tuple[Geodesic, ...], crossing: Point
) -> tuple[Point, ...]:
    active = tuple(
        family
        for family, line in enumerate(lines)
        if line_value(line, crossing).denominator == 1
    )
    if len(active) != 2:
        raise ValueError("crossing is not generic and double")
    first, second = active
    results: list[Point] = []
    zero_origin = (Fraction(0), Fraction(0))
    first_linear = Geodesic("first-linear", lines[first].direction, zero_origin)
    second_linear = Geodesic("second-linear", lines[second].direction, zero_origin)
    for first_sign, second_sign in product((-1, 1), repeat=2):
        direction = level_intersection(
            first_linear,
            Fraction(first_sign),
            second_linear,
            Fraction(second_sign),
        )
        epsilon_bounds = [Fraction(1, 4)]
        for family, line in enumerate(lines):
            if family in active:
                continue
            value = line_value(line, crossing)
            fractional = value - floor(value)
            distance = min(fractional, 1 - fractional)
            slope = Fraction(determinant(line.direction, direction))
            if slope:
                epsilon_bounds.append(distance / (2 * abs(slope)))
        epsilon = min(epsilon_bounds)
        point = (
            crossing[0] + epsilon * direction[0],
            crossing[1] + epsilon * direction[1],
        )
        if any(line_value(line, point).denominator == 1 for line in lines):
            raise AssertionError("local sector perturbation landed on a line")
        results.append(point)
    return tuple(results)


def _canonical_bounds_and_translation(
    lines: tuple[Geodesic, ...], point: Point
) -> tuple[tuple[int, ...], tuple[int, int]]:
    """Normalize strip bounds using the first unimodular direction pair.

    A unimodular pair is only a convenient coordinate choice; it need not be
    the first two input geodesics, and parallel copies may occur anywhere in
    the list.  Arrangements whose direction lattice has no unimodular pair
    are handled separately by :func:`_deck_equivalent_bounds`.
    """

    bounds = tuple(lower for lower, _upper in strip_bounds(lines, point))
    pair = next(
        (
            (left, right)
            for left, right in combinations(range(len(lines)), 2)
            if abs(
                int(determinant(
                    lines[left].direction,
                    lines[right].direction,
                ))
            ) == 1
        ),
        None,
    )
    if pair is None:
        raise ValueError("the directions contain no unimodular pair")
    first_index, second_index = pair
    first = lines[first_index].direction
    second = lines[second_index].direction
    det = int(determinant(first, second))
    first_bound = bounds[first_index]
    second_bound = bounds[second_index]
    x_shift_num = -first_bound * second[0] + first[0] * second_bound
    y_shift_num = first[1] * second_bound - first_bound * second[1]
    if x_shift_num % det or y_shift_num % det:
        raise AssertionError("normalizing deck translation is not integral")
    shift = x_shift_num // det, y_shift_num // det
    normalized = tuple(
        bound + int(determinant(line.direction, shift))
        for line, bound in zip(lines, bounds)
    )
    if (
        normalized[first_index] != 0
        or normalized[second_index] != 0
    ):
        raise AssertionError("failed to normalize the chosen strip levels")
    return normalized, shift


def _deck_equivalent_bounds(
    lines: tuple[Geodesic, ...],
    left: tuple[int, ...],
    right: tuple[int, ...],
) -> bool:
    """Return whether two strip-bound tuples differ by an integral deck shift."""

    pair = next(
        (
            (first, second)
            for first, second in combinations(range(len(lines)), 2)
            if determinant(
                lines[first].direction,
                lines[second].direction,
            ) != 0
        ),
        None,
    )
    if pair is None:
        raise ValueError("the geodesic directions do not span the plane")
    first_index, second_index = pair
    first = lines[first_index].direction
    second = lines[second_index].direction
    det = int(determinant(first, second))
    first_delta = right[first_index] - left[first_index]
    second_delta = right[second_index] - left[second_index]
    x_num = first_delta * second[0] - first[0] * second_delta
    y_num = first_delta * second[1] - first[1] * second_delta
    if x_num % det or y_num % det:
        return False
    shift = x_num // det, y_num // det
    return all(
        right_bound
        == left_bound + int(determinant(line.direction, shift))
        for line, left_bound, right_bound in zip(lines, left, right)
    )


def _fractional_point(point: Point) -> Point:
    return point[0] - floor(point[0]), point[1] - floor(point[1])


@dataclass(frozen=True)
class TorusChamber:
    signature: tuple[int, ...]
    basepoint: Point
    cell: IndependentCell

    @property
    def alternating(self) -> bool:
        return self.cell.alternating


def enumerate_torus_chambers(
    geodesics: Iterable[Geodesic],
    *,
    arrangement: object | None = None,
) -> tuple[TorusChamber, ...]:
    """Enumerate every torus chamber from crossing sectors, exactly.

    A lift is the convex intersection of strips, hence a disk.  Every face
    has a crossing vertex.  The four sectors at every crossing therefore hit
    every chamber.  Strip-bound tuples modulo deck translation distinguish
    the lifted disks.
    """

    lines = tuple(geodesics)
    if arrangement is None:
        parallel = any(
            determinant(left.direction, right.direction) == 0
            for left, right in combinations(lines, 2)
        )
        if parallel:
            builder = getattr(vmr, "make_parallel_geodesic_arrangement", None)
            if builder is None:
                raise ValueError(
                    "parallel chamber enumeration requires "
                    "make_parallel_geodesic_arrangement"
                )
            arrangement = builder(lines)
        else:
            arrangement = vmr.make_geodesic_arrangement(lines)
    chambers: dict[tuple[int, ...], TorusChamber] = {}
    nonunimodular: list[tuple[tuple[int, ...], TorusChamber]] = []
    for vertex in arrangement.vertices:
        for raw_point in _local_sector_basepoints(lines, vertex.position):
            basepoint = _fractional_point(raw_point)
            try:
                signature, _shift = _canonical_bounds_and_translation(
                    lines, raw_point
                )
            except ValueError as error:
                if "no unimodular pair" not in str(error):
                    raise
                signature = tuple(
                    lower for lower, _upper in strip_bounds(lines, raw_point)
                )
                if any(
                    _deck_equivalent_bounds(lines, prior, signature)
                    for prior, _chamber in nonunimodular
                ):
                    continue
                cell = independent_lifted_cell(lines, basepoint)
                chamber = TorusChamber(signature, basepoint, cell)
                nonunimodular.append((signature, chamber))
            else:
                if signature not in chambers:
                    cell = independent_lifted_cell(lines, basepoint)
                    chambers[signature] = TorusChamber(
                        signature, basepoint, cell
                    )

    crossing_count = sum(
        abs(int(determinant(left.direction, right.direction)))
        for left, right in combinations(lines, 2)
    )
    if len(arrangement.vertices) != crossing_count:
        raise AssertionError("crossing count is wrong")
    # A generic 4-valent cellular torus graph has E=2V and F=E-V=V.
    chamber_values = list(chambers.values()) + [
        chamber for _signature, chamber in nonunimodular
    ]
    if len(chamber_values) != crossing_count:
        raise AssertionError(
            f"crossing sectors found {len(chamber_values)} chambers, "
            f"but Euler requires {crossing_count}"
        )
    return tuple(sorted(chamber_values, key=lambda chamber: chamber.signature))
