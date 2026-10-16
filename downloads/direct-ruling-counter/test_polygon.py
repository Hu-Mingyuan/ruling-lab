"""Regression tests for polygon input and count-only orchestration."""

import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import direct_ruling_count as direct
import ruling_polygon as polygon


class PolygonInputTests(unittest.TestCase):
    def test_star_with_consistent_local_turns_is_rejected(self):
        vertices = [[0, 3], [2, -3], [-3, 1], [3, 1], [-2, -3]]
        self.assertTrue(all(
            polygon._cross(vertices[i], vertices[(i + 1) % 5],
                           vertices[(i + 2) % 5]) < 0
            for i in range(5)
        ))
        with self.assertRaisesRegex(ValueError, "strictly convex"):
            polygon._validate_vertices(vertices)

    def test_clockwise_closed_polygon_is_normalized(self):
        vertices = polygon._validate_vertices(
            [[0, 0], [0, 2], [3, 2], [3, 0], [0, 0]]
        )
        self.assertEqual(polygon._polygon_invariants(vertices), (12, 10, 2))
        self.assertGreater(polygon._signed_double_area(vertices), 0)

    def test_collinear_and_concave_inputs_are_rejected(self):
        for vertices in (
            [[0, 0], [1, 0], [2, 0], [0, 2]],
            [[0, 0], [3, 0], [1, 1], [3, 3], [0, 3]],
        ):
            with self.subTest(vertices=vertices), self.assertRaises(ValueError):
                polygon._validate_vertices(vertices)


class CountOrchestrationTests(unittest.TestCase):
    def test_count_only_and_certificate_flags_reach_engine(self):
        engine_result = {
            "phaseCardinality": "FINITE",
            "exactAllDiskRulingCount": 3,
            "exactAnnularRulingCount": 0,
            "exactTotalRulingCount": 3,
        }
        arrangement = object()
        for include in (False, True):
            with self.subTest(certificates=include), \
                    patch.object(direct, "exact_annular_payload", return_value={}) as payload, \
                    patch.object(direct, "_node_executable", return_value=Path("node")), \
                    patch.object(direct.subprocess, "run", return_value=SimpleNamespace(
                        returncode=0, stdout=json.dumps(engine_result), stderr=""
                    )) as run:
                result = direct.count_multigeodesic_rulings(
                    [], boundary_components=3, target_genus=0, basepoint=(0, 0),
                    include_certificates=include, arrangement=arrangement,
                )
                sent = json.loads(run.call_args.kwargs["input"])
                self.assertEqual(sent["countOnly"], not include)
                self.assertEqual(sent["collectAnswers"], include)
                self.assertIs(payload.call_args.kwargs["arrangement"], arrangement)
                self.assertEqual(result.total, 3)

    def test_diagram_is_built_once_per_polygon_with_parallel_genera(self):
        records = [
            {"id": "a", "vertices": [[0, 0], [1, 0], [2, 3]],
             "interior_lattice_points": 1},
            {"id": "b", "vertices": [[0, 0], [1, 0], [2, 5]],
             "interior_lattice_points": 2},
        ]
        for jobs in (1, 3):
            diagrams = {record["id"]: object() for record in records}
            seen = []

            def count(record, *, target_genus, diagram, **kwargs):
                seen.append((record["id"], target_genus, diagram))
                return {
                    "all_disk": 1, "annular": 0, "total": 1,
                    "phase_cardinality": "FINITE", "elapsed_seconds": 0,
                    "crossings": 1, "search_nodes": 1, "shear": 0,
                    "vertical_direction_source": [0, 1],
                    "vertical_direction_count": [0, 1],
                    "sweep_covector_source": [1, 0],
                    "sweep_covector_count": [1, 0],
                    "origin_denominator": 1009, "origin_residues": [1, 2, 3],
                    "origin_attempt": 0, "basepoint": ["1/3", "1/4"],
                    "chambers": 1, "alternating_chambers": 1,
                }

            with self.subTest(jobs=jobs), \
                    patch.object(polygon, "build_polygon_ruling_diagram",
                                 side_effect=lambda label, *_args, **_kwargs: diagrams[label]) as build, \
                    patch.object(polygon, "count_case", side_effect=count):
                result = polygon.count_records(
                    records, selected_genus=None, jobs=jobs,
                    include_certificates=False,
                )
                self.assertEqual(build.call_count, 2)
                self.assertEqual(len(seen), 5)
                self.assertTrue(all(diagram is diagrams[label]
                                    for label, _genus, diagram in seen))
                self.assertEqual(result["algorithm_version"], direct.ALGORITHM_VERSION)
                self.assertEqual(result["polygons"][0]["basepoint"], ["1/3", "1/4"])


if __name__ == "__main__":
    unittest.main()
