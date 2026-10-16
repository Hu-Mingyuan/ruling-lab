#!/usr/bin/env python3
"""Independent geometric regressions for the exact sweep perturbation."""

from __future__ import annotations

from fractions import Fraction as F
from itertools import product
from math import gcd
import unittest

from fig3_direct_ruling_counts import build_polygon_ruling_diagram
from generic_ruling_perturbation import perturb_arrangement
from multigeodesic_alternating_audit import enumerate_torus_chambers
from multigeodesic_annular_input import exact_annular_payload
from verify_multigeodesic_rulings import make_parallel_geodesic_arrangement
from verify_triangle_rulings import Geodesic


def square_lines():
    return (
        Geodesic("a", (1, 1), (F(0), F(0))),
        Geodesic("b", (1, -1), (F(0), F(1, 4))),
        Geodesic("c", (-1, -1), (F(0), F(1, 2))),
        Geodesic("d", (-1, 1), (F(0), F(3, 4))),
    )


def edge_segments(arrangement):
    for edge in arrangement.halfedges:
        if edge.index > edge.twin:
            continue
        origin = arrangement.vertices[edge.origin].position
        target = arrangement.vertices[arrangement.halfedges[edge.twin].origin].position
        yield origin, (target[0] + edge.deck_shift[0] * arrangement.periods[0],
                       target[1] + edge.deck_shift[1] * arrangement.periods[1])


def determinant(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(point, a, b):
    return (determinant(a, b, point) == 0
            and min(a[0], b[0]) <= point[0] <= max(a[0], b[0])
            and min(a[1], b[1]) <= point[1] <= max(a[1], b[1]))


class ExactPerturbationTests(unittest.TestCase):
    def setUp(self):
        self.lines = square_lines()
        self.original = make_parallel_geodesic_arrangement(self.lines)
        self.chambers = enumerate_torus_chambers(self.lines, arrangement=self.original)

    def test_symmetric_square_ties_are_separated_without_changing_graph(self):
        original = self.original
        self.assertLess(len({v.position[0] for v in original.vertices}), len(original.vertices))
        moved, step, count = perturb_arrangement(original, self.chambers)
        self.assertGreater(step, 0)
        self.assertEqual(count, len(original.vertices))
        self.assertEqual(moved.halfedges, original.halfedges)
        self.assertEqual(len({v.position[0] for v in moved.vertices}), count)
        for before, after in zip(original.vertices, moved.vertices):
            self.assertEqual(before.position[1], after.position[1])
            self.assertEqual(before.ports, after.ports)
            self.assertEqual(before.future_ports, after.future_ports)
        self.assertEqual(moved.validation_errors(), ())
        # A second application is exactly a no-op.
        again, second_step, second_count = perturb_arrangement(moved, ())
        self.assertIs(again, moved)
        self.assertEqual((second_step, second_count), (F(0), 0))

    def test_deck_translates_gain_no_extra_intersections(self):
        moved, _, _ = perturb_arrangement(self.original, self.chambers)
        segments = tuple(edge_segments(moved))
        # This fixture's segments have displacement at most one in both
        # coordinates.  These translates exhaust every possible overlap.
        for first, (a, b) in enumerate(segments):
            for second in range(first, len(segments)):
                c0, d0 = segments[second]
                for sx, sy in product(range(-2, 3), repeat=2):
                    if first == second and (sx, sy) == (0, 0):
                        continue
                    c, d = (c0[0] + sx, c0[1] + sy), (d0[0] + sx, d0[1] + sy)
                    cross1, cross2 = determinant(a, b, c), determinant(a, b, d)
                    cross3, cross4 = determinant(c, d, a), determinant(c, d, b)
                    self.assertFalse(cross1 * cross2 < 0 and cross3 * cross4 < 0)
                    for point, left, right in ((a, c, d), (b, c, d), (c, a, b), (d, a, b)):
                        if on_segment(point, left, right):
                            self.assertIn(point, (left, right))
                    shared = set((a, b)) & set((c, d))
                    self.assertLess(len(shared), 2, "two different lifted edges overlap")

    def test_all_old_chamber_basepoints_remain_off_edges(self):
        moved, _, _ = perturb_arrangement(self.original, self.chambers)
        for chamber in self.chambers:
            for a, b in edge_segments(moved):
                for sx, sy in product(range(-2, 3), repeat=2):
                    point = chamber.basepoint[0] + sx, chamber.basepoint[1] + sy
                    self.assertFalse(on_segment(point, a, b))

    def test_incomplete_chamber_certificate_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "cover each directed"):
            perturb_arrangement(self.original, self.chambers[:-1])

    def test_existing_strict_sweep_returns_original_object(self):
        lines = (
            Geodesic("a", (1, 0), (F(0), F(0))),
            Geodesic("b", (1, 1), (F(0), F(1, 7))),
            Geodesic("c", (-2, -1), (F(0), F(3, 7))),
        )
        arrangement = make_parallel_geodesic_arrangement(lines)
        self.assertEqual(len({v.position[0] for v in arrangement.vertices}), len(arrangement.vertices))
        moved, step, count = perturb_arrangement(arrangement, ())
        self.assertIs(moved, arrangement)
        self.assertEqual((step, count), (F(0), 0))

    def test_public_builder_separates_structurally_tied_polygon(self):
        diagram = build_polygon_ruling_diagram(
            "forced-tie", [(0, 0), (2, -1), (4, 0), (2, 1)]
        )
        self.assertGreater(diagram.perturbation_step, 0)
        self.assertEqual(diagram.perturbed_vertices, diagram.crossings)
        self.assertEqual(len({vertex.position[0] for vertex in diagram.arrangement.vertices}),
                         diagram.crossings)
        original = make_parallel_geodesic_arrangement(diagram.geodesics)
        self.assertLess(len({vertex.position[0] for vertex in original.vertices}),
                        diagram.crossings)

    def test_serialized_tangents_follow_actual_perturbed_edges(self):
        diagram = build_polygon_ruling_diagram(
            "forced-tie", [(0, 0), (2, -1), (4, 0), (2, 1)]
        )
        payload = exact_annular_payload(
            diagram.geodesics, boundary_components=len(diagram.geodesics),
            target_genus=0, basepoint=diagram.alternating_chambers[0].basepoint,
            arrangement=diagram.arrangement,
        )
        raw = payload["arrangement"]
        original_directions = {line.name: line.direction for line in diagram.geodesics}
        has_bent_edge = False
        for edge in diagram.arrangement.halfedges:
            start = diagram.arrangement.vertices[edge.origin].position
            finish = diagram.arrangement.vertices[
                diagram.arrangement.halfedges[edge.twin].origin
            ].position
            dx, dy = finish[0] + edge.deck_shift[0] - start[0], finish[1] + edge.deck_shift[1] - start[1]
            vx, vy = raw["portVectorX"][edge.index], raw["portVectorY"][edge.index]
            self.assertEqual(vx * dy - vy * dx, 0)
            self.assertGreater(vx * dx + vy * dy, 0)
            self.assertEqual(gcd(abs(vx), abs(vy)), 1)
            old_x, old_y = original_directions[edge.strand]
            has_bent_edge |= old_x * dy - old_y * dx != 0
        self.assertTrue(has_bent_edge, "fixture must detect stale geodesic tangent serialization")

    def test_source_1009_parallelogram_selects_safe_origin_denominator(self):
        # Primitive basis vectors (1009,1),(1008,1) have determinant one.
        # The previous fixed denominator 1009 made opposite copies coincide.
        diagram = build_polygon_ruling_diagram(
            "source-1009", [(0, 0), (1009, 1), (2017, 2), (1008, 1)]
        )
        self.assertEqual(diagram.shear, 0)
        self.assertGreater(diagram.origin_denominator, 1009)
        self.assertTrue(all(gcd(abs(line.direction[0]), diagram.origin_denominator) == 1
                            for line in diagram.geodesics))
        self.assertEqual(diagram.crossings, 4)
        self.assertTrue(diagram.alternating_chambers)


if __name__ == "__main__":
    unittest.main()
