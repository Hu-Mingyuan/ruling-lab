#!/usr/bin/env python3
"""Minimal exact half-edge geometry used by the direct ruling counter."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction

DeckVector = tuple[int, int]
Point = tuple[Fraction, Fraction]


def _neg(value: DeckVector) -> DeckVector:
    return -value[0], -value[1]


@dataclass(frozen=True)
class TorusHalfEdge:
    """One directed edge of a finite torus graph.

    ``deck_shift`` is the translation of the target representative relative
    to the origin representative.  ``coorientation_sign`` is ``+1`` when the
    coorientation lies to the left of this directed edge and ``-1`` when it
    lies to the right.  Reversing a half-edge reverses both values.
    """

    index: int
    origin: int
    twin: int
    deck_shift: DeckVector
    coorientation_sign: int
    strand: str = ""


@dataclass(frozen=True)
class TorusVertex:
    """A four-valent crossing with ports in counterclockwise order."""

    index: int
    position: Point
    ports: tuple[int, int, int, int]
    future_ports: tuple[int, int]
    label: str = ""
    horizontal: int | None = None
    diagonal: int | None = None


@dataclass(frozen=True)
class TorusArrangement:
    """Finite cooriented four-valent arrangement on a rectangular torus."""

    periods: Point
    vertices: tuple[TorusVertex, ...]
    halfedges: tuple[TorusHalfEdge, ...]
    crossing_order: tuple[int, ...]
    name: str = ""

    def validation_errors(self) -> tuple[str, ...]:
        """Return structural errors rather than failing later while tracing."""

        errors: list[str] = []
        if self.periods[0] <= 0 or self.periods[1] <= 0:
            errors.append("the two torus periods must be positive")
        if tuple(vertex.index for vertex in self.vertices) != tuple(
            range(len(self.vertices))
        ):
            errors.append("vertex indices must be consecutive")
        if tuple(edge.index for edge in self.halfedges) != tuple(
            range(len(self.halfedges))
        ):
            errors.append("half-edge indices must be consecutive")
        if sorted(self.crossing_order) != list(range(len(self.vertices))):
            errors.append("crossing_order must contain every vertex once")

        port_occurrences: Counter[int] = Counter()
        for vertex in self.vertices:
            if len(set(vertex.ports)) != 4:
                errors.append(f"vertex {vertex.index} does not have four ports")
            for port in vertex.ports:
                if not 0 <= port < len(self.halfedges):
                    errors.append(f"vertex {vertex.index} has invalid port {port}")
                    continue
                port_occurrences[port] += 1
                if self.halfedges[port].origin != vertex.index:
                    errors.append(
                        f"port {port} has the wrong origin for vertex {vertex.index}"
                    )
            if (
                len(set(vertex.future_ports)) != 2
                or not set(vertex.future_ports) <= set(vertex.ports)
            ):
                errors.append(
                    f"vertex {vertex.index} must have exactly two future ports"
                )

        for edge in self.halfedges:
            if edge.coorientation_sign not in {-1, 1}:
                errors.append(
                    f"half-edge {edge.index} has a non-binary coorientation sign"
                )
            if not 0 <= edge.twin < len(self.halfedges):
                errors.append(f"half-edge {edge.index} has invalid twin")
                continue
            twin = self.halfedges[edge.twin]
            if twin.twin != edge.index:
                errors.append(f"half-edge {edge.index} has a non-reciprocal twin")
            if twin.deck_shift != _neg(edge.deck_shift):
                errors.append(f"half-edge {edge.index} has inconsistent deck shift")
            if twin.coorientation_sign != -edge.coorientation_sign:
                errors.append(
                    f"half-edge {edge.index} has inconsistent coorientation"
                )
            if port_occurrences[edge.index] != 1:
                errors.append(
                    f"half-edge {edge.index} occurs in {port_occurrences[edge.index]} ports"
                )

        if not errors:
            future_ports = {
                port
                for vertex in self.vertices
                for port in vertex.future_ports
            }
            for edge in self.halfedges:
                if (edge.index in future_ports) == (edge.twin in future_ports):
                    errors.append(
                        f"half-edge {edge.index} and its twin must point to "
                        "opposite sweep directions"
                    )
            for vertex in self.vertices:
                if any(
                    (left in set(vertex.future_ports))
                    == (right in set(vertex.future_ports))
                    for left, right in nonswitch_pairs(vertex)
                ):
                    errors.append(
                        f"an original strand at vertex {vertex.index} is not "
                        "transverse to the sweep"
                    )
                candidates = _switch_candidates(vertex)
                admissible = [
                    pairs
                    for pairs in candidates
                    if all(
                        self.halfedges[left].coorientation_sign
                        == -self.halfedges[right].coorientation_sign
                        for left, right in pairs
                    )
                ]
                if len(admissible) != 1:
                    errors.append(
                        f"vertex {vertex.index} has {len(admissible)} "
                        "coorientation-preserving switch resolutions"
                    )
        return tuple(dict.fromkeys(errors))


def _switch_candidates(
    vertex: TorusVertex,
) -> tuple[
    tuple[tuple[int, int], tuple[int, int]],
    tuple[tuple[int, int], tuple[int, int]],
]:
    p0, p1, p2, p3 = vertex.ports
    return (((p0, p1), (p2, p3)), ((p0, p3), (p1, p2)))


def nonswitch_pairs(
    vertex: TorusVertex,
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Pair opposite ports, following the two original strands."""

    p0, p1, p2, p3 = vertex.ports
    return (p0, p2), (p1, p3)


def switch_pairs(
    arrangement: TorusArrangement, vertex: TorusVertex
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Derive the unique switch resolution from the coorientations.

    If a circuit arrives through ``left`` and leaves through ``right``, its
    coorientation sign is preserved precisely when the signs on the two
    outward ports are opposite.  Exactly one of the two adjacent resolutions
    has this property for a valid cooriented crossing.
    """

    admissible = [
        pairs
        for pairs in _switch_candidates(vertex)
        if all(
            arrangement.halfedges[left].coorientation_sign
            == -arrangement.halfedges[right].coorientation_sign
            for left, right in pairs
        )
    ]
    if len(admissible) != 1:
        raise ValueError(
            f"vertex {vertex.index} does not have a unique switch resolution"
        )
    return admissible[0]


def switch_orientation(
    arrangement: TorusArrangement, vertex: TorusVertex
) -> str:
    """Return ``vertical`` for two tip arcs and ``horizontal`` otherwise."""

    future = set(vertex.future_ports)
    patterns = []
    for left, right in switch_pairs(arrangement, vertex):
        left_future = left in future
        right_future = right in future
        patterns.append(left_future == right_future)
    if patterns == [True, True]:
        return "vertical"
    if patterns == [False, False]:
        return "horizontal"
    return "invalid"
