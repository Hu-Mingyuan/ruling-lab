#!/usr/bin/env python3
"""Serialize a primitive multigeodesic diagram for the exact annular solver."""

from __future__ import annotations

from fractions import Fraction
from math import gcd, lcm
from typing import Iterable, Sequence

import ruling_smoothing as rs
from verify_multigeodesic_rulings import make_parallel_geodesic_arrangement
from verify_triangle_rulings import Geodesic, Point


def _fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def exact_annular_payload(
    geodesics: Iterable[Geodesic],
    *,
    boundary_components: int,
    target_genus: int,
    basepoint: Point,
    masks: Sequence[int | str] | None = None,
    tau: tuple[int, int] = (1, 0),
    sigma: tuple[int, int] = (0, 1),
    name: str = "primitive multigeodesic arrangement",
    arrangement: rs.TorusArrangement | None = None,
) -> dict[str, object]:
    """Return JSON-ready data for ``tmp/generic_annular_exact_search.js``."""

    lines = tuple(geodesics)
    if any(isinstance(entry, bool) or not isinstance(entry, int)
           for entry in (*tau, *sigma)):
        raise ValueError("tau and sigma must be integral covectors")
    if gcd(abs(tau[0]), abs(tau[1])) != 1:
        raise ValueError("tau must be a primitive covector")
    if abs(tau[0] * sigma[1] - tau[1] * sigma[0]) != 1:
        raise ValueError("tau and sigma must form a unimodular covector basis")
    vertical_direction = (-tau[1], tau[0])
    sweep_agreements: set[int] = set()
    for line in lines:
        sweep = tau[0] * line.direction[0] + tau[1] * line.direction[1]
        if sweep == 0:
            raise ValueError(
                f"vertical direction {vertical_direction} is parallel to "
                f"geodesic {line.name!r}"
            )
        # The arrangement uses tau=(1,0), whose kernel is the count-frame
        # vertical direction (0,1).  A different tau is safe only in the same
        # (or globally opposite) covector chamber, so that every local
        # horizontal/vertical switch type is unchanged.
        sweep_agreements.add(
            (1 if sweep > 0 else -1)
            * (1 if line.direction[0] > 0 else -1)
        )
    if len(sweep_agreements) != 1:
        raise ValueError(
            "tau is incompatible with the arrangement's vertical direction; "
            "their future directions do not agree up to global reversal"
        )
    same_sweep_orientation = next(iter(sweep_agreements)) > 0
    # Polygon edges of lattice length > 1 contribute several translated
    # primitive geodesics with the same (or opposite) direction.  The generic
    # parallel-family builder is therefore the right common entry point; it
    # also handles the pairwise-nonparallel special case.
    if arrangement is None:
        arrangement = make_parallel_geodesic_arrangement(lines, name=name)
    errors = arrangement.validation_errors()
    if errors:
        raise ValueError("invalid prepared arrangement: " + "; ".join(errors))
    if arrangement.periods != (1, 1):
        raise ValueError("the serialized arrangement must use the unit torus")
    if {edge.strand for edge in arrangement.halfedges} != {line.name for line in lines}:
        raise ValueError("prepared arrangement has different geodesic labels")
    scale = 1
    for vertex in arrangement.vertices:
        scale = lcm(
            scale,
            vertex.position[0].denominator,
            vertex.position[1].denominator,
        )

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
        # A genuine generic perturbation bends the original geodesics.
        # Tangents must therefore come from the actual lifted edge, not
        # from its original straight-line homology class.
        start = arrangement.vertices[edge.origin].position
        finish = arrangement.vertices[arrangement.halfedges[edge.twin].origin].position
        dx = int((finish[0] + edge.deck_shift[0] - start[0]) * scale)
        dy = int((finish[1] + edge.deck_shift[1] - start[1]) * scale)
        divisor = gcd(abs(dx), abs(dy))
        if divisor == 0:
            raise ValueError("prepared arrangement has a zero-length edge")
        direction = (dx // divisor, dy // divisor)
        future = edge.index in arrangement.vertices[edge.origin].future_ports
        if not direction[0] or (direction[0] > 0) != future:
            raise ValueError("prepared edge tangent disagrees with its future port")
        sweep = tau[0] * direction[0] + tau[1] * direction[1]
        expected_positive = ((direction[0] > 0) if same_sweep_orientation
                             else (direction[0] < 0))
        if not sweep or (sweep > 0) != expected_positive:
            raise ValueError("tau is incompatible with a perturbed edge tangent")
        port_vector_x.append(direction[0])
        port_vector_y.append(direction[1])

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
        # All vertices and edge endpoints are integral after scaling.  A
        # tau-step of 1/10 therefore stays strictly inside each incident edge.
        "typeBEpsilon": "1/10",
        "balanceTarget": 0,
        "requireBalance": True,
        "collectAnswers": True,
    }
    if masks is not None:
        payload["masks"] = [
            value if isinstance(value, str) else f"0x{value:x}" for value in masks
        ]
    return payload
