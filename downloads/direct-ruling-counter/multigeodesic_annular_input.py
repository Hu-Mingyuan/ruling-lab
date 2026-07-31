#!/usr/bin/env python3
"""Serialize a primitive multigeodesic diagram for the exact annular solver."""

from __future__ import annotations

from fractions import Fraction
from math import floor, gcd, lcm
from typing import Iterable, Sequence

import ruling_smoothing as rs
from verify_multigeodesic_rulings import make_parallel_geodesic_arrangement
from verify_triangle_rulings import Geodesic, Point


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def _crossing_levels(
    arrangement: object,
    tau: tuple[int, int],
) -> tuple[Fraction, ...]:
    """Return the crossing levels in the circle ``R/Z``."""

    return tuple(
        (tau[0] * vertex.position[0] + tau[1] * vertex.position[1]) % 1
        for vertex in arrangement.vertices
    )


def _sweep_agreements(
    lines: Sequence[Geodesic],
    tau: tuple[int, int],
) -> set[int]:
    """Compare ``tau`` with the x-sweep on every geodesic direction."""

    agreements: set[int] = set()
    for line in lines:
        dx, dy = line.direction
        if dx == 0:
            raise ValueError(
                f"the input count frame is tangent to geodesic {line.name!r}"
            )
        sweep = tau[0] * dx + tau[1] * dy
        if sweep == 0:
            return set()
        agreements.add(
            (1 if sweep > 0 else -1) * (1 if dx > 0 else -1)
        )
    return agreements


def _unimodular_companion(tau: tuple[int, int]) -> tuple[int, int]:
    """Choose ``sigma`` with ``det(tau, sigma)=1``."""

    a, b = tau
    old_r, remainder = abs(a), abs(b)
    old_s, coefficient_s = 1, 0
    old_t, coefficient_t = 0, 1
    while remainder:
        quotient = old_r // remainder
        old_r, remainder = remainder, old_r - quotient * remainder
        old_s, coefficient_s = coefficient_s, old_s - quotient * coefficient_s
        old_t, coefficient_t = coefficient_t, old_t - quotient * coefficient_t
    bezout_a = old_s if a >= 0 else -old_s
    bezout_b = old_t if b >= 0 else -old_t
    if a * bezout_a + b * bezout_b != 1:
        raise ValueError("tau must be a primitive covector")
    return -bezout_b, bezout_a


def _automatic_generic_sweep(
    lines: Sequence[Geodesic],
    arrangement: object,
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Find a primitive integral perturbation of the x-sweep.

    We retain ``(1,0)`` whenever it is already generic.  Otherwise the
    deterministic search tries ``(M,1)`` and ``(M,-1)`` in increasing order.
    Large ``M`` lies in the x-sweep chamber because every input direction has
    nonzero x-component.
    """

    search_span = max(4096, 16 * len(arrangement.vertices) ** 2)
    def candidates() -> Iterable[tuple[int, int]]:
        yield 1, 0
        for sign in (1, -1):
            chamber_threshold = max(
                (-Fraction(sign * line.direction[1], line.direction[0])
                 for line in lines),
                default=Fraction(0),
            )
            first_coefficient = max(1, floor(chamber_threshold) + 1)
            yield from (
                (coefficient, sign)
                for coefficient in range(
                    first_coefficient,
                    first_coefficient + search_span,
                )
            )

    for candidate in candidates():
        if _sweep_agreements(lines, candidate) != {1}:
            continue
        levels = _crossing_levels(arrangement, candidate)
        if len(levels) == len(set(levels)):
            return candidate, _unimodular_companion(candidate)
    raise ValueError(
        "no generic integral sweep was found in the x-sweep chamber; "
        "perturb the geodesic translates before counting"
    )


def exact_annular_payload(
    geodesics: Iterable[Geodesic],
    *,
    boundary_components: int,
    target_genus: int,
    basepoint: Point,
    masks: Sequence[int | str] | None = None,
    tau: tuple[int, int] | None = None,
    sigma: tuple[int, int] | None = None,
    name: str = "primitive multigeodesic arrangement",
) -> dict[str, object]:
    """Return JSON-ready data for ``tmp/generic_annular_exact_search.js``."""

    lines = tuple(geodesics)
    # Polygon edges of lattice length > 1 contribute several translated
    # primitive geodesics with the same (or opposite) direction.  The generic
    # parallel-family builder is therefore the right common entry point; it
    # also handles the pairwise-nonparallel special case.
    arrangement = make_parallel_geodesic_arrangement(lines, name=name)
    automatic_sweep = tau is None
    if automatic_sweep:
        if sigma is not None:
            raise ValueError("sigma cannot be supplied when tau is automatic")
        tau, sigma = _automatic_generic_sweep(lines, arrangement)
    else:
        if (
            not isinstance(tau, tuple)
            or len(tau) != 2
            or any(isinstance(entry, bool) or not isinstance(entry, int)
                   for entry in tau)
        ):
            raise ValueError("tau must be a pair of integers")
        if gcd(abs(tau[0]), abs(tau[1])) != 1:
            raise ValueError("tau must be a primitive covector")
        if sigma is None:
            sigma = _unimodular_companion(tau)
        elif (
            not isinstance(sigma, tuple)
            or len(sigma) != 2
            or any(isinstance(entry, bool) or not isinstance(entry, int)
                   for entry in sigma)
        ):
            raise ValueError("sigma must be a pair of integers")

    assert tau is not None and sigma is not None
    if abs(tau[0] * sigma[1] - tau[1] * sigma[0]) != 1:
        raise ValueError("tau and sigma must form a unimodular covector basis")
    vertical_direction = (-tau[1], tau[0])
    sweep_agreements = _sweep_agreements(lines, tau)
    if not sweep_agreements:
        raise ValueError(
            f"vertical direction {vertical_direction} is parallel to a geodesic"
        )
    if len(sweep_agreements) != 1:
        raise ValueError(
            "tau is incompatible with the arrangement's vertical direction; "
            "their future directions do not agree up to global reversal"
        )
    crossing_levels = _crossing_levels(arrangement, tau)
    if len(crossing_levels) != len(set(crossing_levels)):
        raise ValueError(
            "tau puts distinct crossings on the same sweep level; omit tau "
            "to select an automatic generic perturbation"
        )

    direction_by_name = {line.name: line.direction for line in lines}
    scale = 1
    for vertex in arrangement.vertices:
        scale = lcm(
            scale,
            vertex.position[0].denominator,
            vertex.position[1].denominator,
        )

    scaled_levels = sorted(int(level * scale) for level in crossing_levels)
    if scaled_levels:
        cyclic_gaps = [
            finish - start
            for start, finish in zip(
                scaled_levels,
                (*scaled_levels[1:], scaled_levels[0] + scale),
            )
        ]
        minimum_gap = min(cyclic_gaps)
    else:
        minimum_gap = scale
    type_b_epsilon = min(Fraction(1, 10), Fraction(minimum_gap, 4))
    vertices = []
    for vertex in arrangement.vertices:
        vertices.append(
            {
                "index": vertex.index,
                "label": vertex.label,
                "x": int(vertex.position[0] * scale),
                "y": int(vertex.position[1] * scale),
                "ports": list(vertex.ports),
                "futurePorts": list(vertex.future_ports),
                "nonswitchPairs": [list(pair) for pair in rs.nonswitch_pairs(vertex)],
                "switchPairs": [
                    list(pair) for pair in rs.switch_pairs(arrangement, vertex)
                ],
                "orientation": (
                    "H"
                    if rs.switch_orientation(arrangement, vertex) == "horizontal"
                    else "V"
                ),
            }
        )

    port_vector_x: list[int] = []
    port_vector_y: list[int] = []
    for edge in arrangement.halfedges:
        direction = direction_by_name[edge.strand]
        port_vector_x.append(edge.coorientation_sign * direction[0])
        port_vector_y.append(edge.coorientation_sign * direction[1])

    edges = []
    for edge in arrangement.halfedges:
        if edge.index > edge.twin:
            continue
        edges.append(
            {
                "a": edge.index,
                "b": edge.twin,
                "sx": edge.deck_shift[0],
                "sy": edge.deck_shift[1],
                "family": edge.strand,
            }
        )

    target_chi = 2 - 2 * target_genus - boundary_components
    scaled_basepoint = (
        Fraction(basepoint[0]) * scale,
        Fraction(basepoint[1]) * scale,
    )
    raw_arrangement = {
        "name": name,
        "L": scale,
        "vertices": vertices,
        "edges": edges,
        "portCount": len(arrangement.halfedges),
        "portVertex": [edge.origin for edge in arrangement.halfedges],
        "portVectorX": port_vector_x,
        "portVectorY": port_vector_y,
        "portSign": [edge.coorientation_sign for edge in arrangement.halfedges],
        "portTwin": [edge.twin for edge in arrangement.halfedges],
        "portShiftX": [edge.deck_shift[0] for edge in arrangement.halfedges],
        "portShiftY": [edge.deck_shift[1] for edge in arrangement.halfedges],
        # This is a stable mask-bit encoding, not the order of the perturbed
        # sweep.  Keeping it unchanged preserves certificate identifiers.
        "crossingOrder": list(arrangement.crossing_order),
        "boundaryComponents": boundary_components,
        "targetGenus": target_genus,
        "targetChi": target_chi,
        "requiredH": -target_chi,
    }
    payload: dict[str, object] = {
        "arrangement": raw_arrangement,
        "basePoint": [_fraction_text(value) for value in scaled_basepoint],
        "tau": list(tau),
        "sigma": list(sigma),
        "verticalDirection": list(vertical_direction),
        "sweepSelection": "automatic-generic" if automatic_sweep else "explicit",
        # Crossing levels are distinct integers modulo L after scaling.  This
        # probe is at most one quarter of their minimum cyclic separation.
        "typeBEpsilon": _fraction_text(type_b_epsilon),
        "balanceTarget": 0,
        "requireBalance": True,
        "collectAnswers": True,
    }
    if masks is not None:
        payload["masks"] = [
            value if isinstance(value, str) else f"0x{value:x}" for value in masks
        ]
    return payload
