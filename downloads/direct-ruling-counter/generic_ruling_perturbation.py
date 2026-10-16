#!/usr/bin/env python3
"""Separate sweep events by a certified exact piecewise-linear isotopy.

The original convex chamber decomposition is triangulated from its supplied
interior basepoints.  Crossing vertices move horizontally and all chamber
basepoints stay fixed.  Positive oriented areas on every triangle certify a
periodic orientation-preserving PL map of degree one, hence a homeomorphism
of the torus.  In particular no additional crossing is introduced.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import replace
from fractions import Fraction
from math import lcm
from typing import Iterable, TYPE_CHECKING

import ruling_smoothing as rs

if TYPE_CHECKING:
    from multigeodesic_alternating_audit import TorusChamber


def perturb_arrangement(
    arrangement: rs.TorusArrangement,
    chambers: Iterable[TorusChamber],
) -> tuple[rs.TorusArrangement, Fraction, int]:
    """Return an isotopic straight-edge arrangement with distinct x events.

    The half-edge graph, deck shifts, cyclic port order, and future ports are
    unchanged.  Callers must derive tangent vectors from the returned edge
    endpoints: after perturbation a strand is piecewise linear, not geodesic.
    The chambers must be the complete original convex chamber decomposition.
    The returned step is zero when the input already has distinct x values.
    """

    vertices = arrangement.vertices
    if len({vertex.position[0] for vertex in vertices}) == len(vertices):
        return arrangement, Fraction(0), 0
    errors = arrangement.validation_errors()
    if errors:
        raise ValueError("invalid input arrangement: " + "; ".join(errors))
    period_x, period_y = arrangement.periods
    if any(not (0 <= vertex.position[0] < period_x
                and 0 <= vertex.position[1] < period_y)
           for vertex in vertices):
        raise ValueError("crossing representatives must lie in the fundamental rectangle")
    position_to_index = {vertex.position: vertex.index for vertex in vertices}
    if len(position_to_index) != len(vertices):
        raise ValueError("crossing positions must be distinct")
    count = len(vertices)
    coefficients = tuple(index + 1 for index in range(count))
    bounds = [(period_x - max(vertex.position[0] for vertex in vertices)) / (4 * count)]
    distinct_x = sorted({vertex.position[0] for vertex in vertices})
    bounds.extend((right - left) / (4 * count)
                  for left, right in zip(distinct_x, distinct_x[1:]))

    def representative(point: rs.Point) -> rs.Point:
        return point[0] % period_x, point[1] % period_y

    def boundary_key(left: rs.Point, right: rs.Point) -> tuple[int, int, int, int]:
        try:
            origin = position_to_index[representative(left)]
            target = position_to_index[representative(right)]
        except KeyError as error:
            raise ValueError("a chamber corner is not an arrangement crossing") from error
        origin_position = vertices[origin].position
        target_position = vertices[target].position
        shift_x = ((right[0] - left[0])
                   - (target_position[0] - origin_position[0])) / period_x
        shift_y = ((right[1] - left[1])
                   - (target_position[1] - origin_position[1])) / period_y
        if shift_x.denominator != 1 or shift_y.denominator != 1:
            raise ValueError("a chamber edge has a fractional deck shift")
        return origin, target, int(shift_x), int(shift_y)

    expected_sides: Counter[tuple[int, int, int, int]] = Counter()
    for edge in arrangement.halfedges:
        target = arrangement.halfedges[edge.twin].origin
        expected_sides[(edge.origin, target, *edge.deck_shift)] += 1
        delta_x = (vertices[target].position[0] + edge.deck_shift[0] * period_x
                   - vertices[edge.origin].position[0])
        if delta_x == 0:
            raise ValueError("an original edge is tangent to the vertical sweep")
        derivative = coefficients[target] - coefficients[edge.origin]
        if derivative:
            bounds.append(abs(delta_x) / (2 * abs(derivative)))

    actual_sides: Counter[tuple[int, int, int, int]] = Counter()
    # (old oriented area, derivative) is sufficient because y never changes.
    fan_areas: list[tuple[Fraction, Fraction]] = []
    for chamber in chambers:
        center_x, center_y = chamber.basepoint
        polygon = chamber.cell.vertices
        if len(polygon) < 3:
            raise ValueError("a chamber must be a polygon")
        for left, right in zip(polygon, polygon[1:] + polygon[:1]):
            key = boundary_key(left, right)
            actual_sides[key] += 1
            origin, target = key[:2]
            area = ((left[0] - center_x) * (right[1] - center_y)
                    - (left[1] - center_y) * (right[0] - center_x))
            if area <= 0:
                raise ValueError("a chamber basepoint is not strictly inside its CCW fan")
            derivative = (coefficients[origin] * (right[1] - center_y)
                          - (left[1] - center_y) * coefficients[target])
            fan_areas.append((area, derivative))
            if derivative:
                bounds.append(area / (2 * abs(derivative)))
    if actual_sides != expected_sides:
        raise ValueError("chambers do not cover each directed arrangement edge exactly once")
    if not fan_areas:
        raise ValueError("the chamber decomposition is empty")

    bound = min(bounds)
    if bound <= 0:
        raise ValueError("no positive perturbation bound exists")
    scale = lcm(period_x.denominator, period_y.denominator,
                *(coordinate.denominator for vertex in vertices
                  for coordinate in vertex.position))
    inverse_bound = Fraction(1, 1) / (scale * bound)
    denominator_factor = inverse_bound.numerator // inverse_bound.denominator + 1
    step = Fraction(1, scale * denominator_factor)
    if not 0 < step < bound:
        raise AssertionError("the chosen perturbation does not satisfy its exact bound")
    moved = tuple(replace(vertex, position=(vertex.position[0]
                                           + coefficients[vertex.index] * step,
                                           vertex.position[1]))
                  for vertex in vertices)
    if len({vertex.position[0] for vertex in moved}) != count:
        raise AssertionError("perturbed sweep events are not distinct")
    if any(not 0 <= vertex.position[0] < period_x for vertex in moved):
        raise AssertionError("a perturbation crossed the fundamental rectangle boundary")
    if any(area + step * derivative <= 0 for area, derivative in fan_areas):
        raise AssertionError("a perturbed fan triangle changed orientation")
    future = {port for vertex in vertices for port in vertex.future_ports}
    for edge in arrangement.halfedges:
        target = arrangement.halfedges[edge.twin].origin
        delta_x = (moved[target].position[0] + edge.deck_shift[0] * period_x
                   - moved[edge.origin].position[0])
        if (delta_x > 0) != (edge.index in future) or delta_x == 0:
            raise AssertionError("a perturbation changed a future/past edge direction")
    result = replace(arrangement, vertices=moved,
                     crossing_order=tuple(sorted(range(count),
                                                 key=lambda index: moved[index].position[0])))
    if result.validation_errors():
        raise AssertionError("the perturbed arrangement failed structural validation")
    return result, step, count
