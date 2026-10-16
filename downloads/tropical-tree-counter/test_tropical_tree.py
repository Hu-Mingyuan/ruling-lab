"""Regression tests for the independent tropical counter and its CLI."""
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import tropical_tree as tt

TRIANGLE = [[-1, -1], [1, 0], [0, 1]]
SIMPLEX = [[0, 0], [1, 0], [0, 1]]
SQUARE = [[0, 0], [1, 0], [1, 1], [0, 1]]
PACKAGE = Path(__file__).resolve().parent


def polynomial(result):
    return {term["t_power"]: term["coefficient"] for term in result["refined_terms"]}


def ruling_payload(vertices=TRIANGLE):
    return {"polygons": [{"id": "example", "vertices": deepcopy(vertices),
                         "genera": [{"genus": 0, "total": 3, "phase_cardinality": "FINITE"},
                                    {"genus": 1, "total": 1, "phase_cardinality": "FINITE"}]}]}


class GeometryTests(unittest.TestCase):
    def test_clockwise_and_closing_vertex(self):
        raw = list(reversed(TRIANGLE))
        raw.append(raw[0])
        self.assertEqual(tt._polygon_key(raw), tt._polygon_key(TRIANGLE))
        self.assertEqual(polynomial(tt.count_tropical(raw)), polynomial(tt.count_tropical(TRIANGLE)))

    def test_strict_convexity_rejects_stars_and_bad_coordinates(self):
        invalid = [
            [[0, 3], [2, -3], [-3, 1], [3, 1], [-2, -3]],
            [[0, 0], [2, 0], [1, 1], [2, 2], [0, 2]],
            [[0, 0], [1, 0], [2, 0], [0, 1]],
            [[0, 0], [True, 0], [0, 1]],
            [[0, 0], [1.0, 0], [0, 1]],
            [[0, 0], [1, 0]],
        ]
        for vertices in invalid:
            with self.subTest(vertices=vertices), self.assertRaises(ValueError):
                tt.validate_vertices(vertices)

    def test_input_forms_and_duplicate_ids(self):
        for payload in (TRIANGLE, {"vertices": TRIANGLE}, {"polygons": [{"vertices": TRIANGLE}]}):
            self.assertEqual(len(tt.parse_records(payload)), 1)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            tt.parse_records({"polygons": [{"id": "x", "vertices": SIMPLEX},
                                            {"id": "x", "vertices": TRIANGLE}]})

    def test_end_limit_checked_before_enumeration(self):
        with self.assertRaisesRegex(ValueError, "exceeding max_ends"):
            tt.count_tropical([[0, 0], [1000000000, 0], [0, 1]])
        with self.assertRaises(ValueError):
            tt.count_tropical(SQUARE, max_ends=3)
        self.assertEqual(tt.count_tropical(SQUARE, max_ends=4)["unrefined_count"], 1)


class TropicalTests(unittest.TestCase):
    def test_independent_known_counts_two_seeds(self):
        cases = [
            (SIMPLEX, {0: 1}), (SQUARE, {0: 1}),
            (TRIANGLE, {-2: 1, 0: 1, 2: 1}),
            ([[0, 0], [1, 0], [2, 5]], {-4: 1, -2: 1, 0: 1, 2: 1, 4: 1}),
            ([[0, 0], [3, 0], [0, 3]], {-2: 1, 0: 7, 2: 1}),
            ([[0, 0], [3, 0], [3, 2], [0, 2]], {-4: 1, -2: 8, 0: 30, 2: 8, 4: 1}),
        ]
        with patch.object(tt, "_ruling_terms", side_effect=AssertionError("not part of the count")):
            for vertices, expected in cases:
                for seed in (17471, 98189):
                    with self.subTest(vertices=vertices, seed=seed):
                        result = tt.count_tropical(vertices, seed=seed)
                        self.assertEqual(polynomial(result), expected)
                        self.assertEqual(result["unrefined_count"], sum(expected.values()))
                        self.assertEqual(sum(end["moment"] for end in result["ends"]), 0)
                        self.assertEqual(sum(end["root"] for end in result["ends"]), 1)
                        self.assertTrue(result["ends"][-1]["root"])

    def test_explicit_rational_moments(self):
        result = tt.count_tropical(TRIANGLE, moments=[Fraction(1, 2), Fraction(2, 3)])
        self.assertIsNone(result["seed"])
        self.assertEqual(result["moment_source"], "explicit")
        self.assertEqual([end["moment"] for end in result["ends"]], ["1/2", "2/3", "-7/6"])
        self.assertEqual(polynomial(result), {-2: 1, 0: 1, 2: 1})

    def test_invalid_moment_input(self):
        for moments in ([1], [1, 2, 3], [1.5, 2], [True, 2]):
            with self.subTest(moments=moments), self.assertRaises(ValueError):
                tt.count_tropical(TRIANGLE, moments=moments)

    def test_balanced_subset_is_rejected(self):
        with self.assertRaisesRegex(tt.NonGenericMomentsError, "balanced subset"):
            tt.count_tropical(SQUARE, moments=[0, 1, 0, -1])

    def test_coincident_root_end_is_rejected(self):
        rectangle = [[0, 0], [1, 0], [1, 2], [0, 2]]
        with self.assertRaisesRegex(tt.NonGenericMomentsError, "parallel ends"):
            tt.count_tropical(rectangle, moments=[1, 2, 3, 4, -5, -5])

    def test_parallel_aggregate_support_is_rejected(self):
        hexagon = [[0, 0], [0, 1], [-1, 2], [-2, 2], [-2, 1], [-1, 0]]
        with self.assertRaisesRegex(tt.NonGenericMomentsError, "parallel subset supports"):
            tt.count_tropical(hexagon, moments=[1, 3, 2, 10, 20, -36])

    def test_zero_length_edge_is_rejected(self):
        pentagon = [[0, 0], [1, 0], [2, 1], [1, 2], [0, 1]]
        with self.assertRaisesRegex(tt.NonGenericMomentsError, "zero-length"):
            tt.count_tropical(pentagon, moments=[2, 5, 1, 10])


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.result = tt.count_records(tt.parse_records({"vertices": TRIANGLE}))

    def test_complete_finite_comparison(self):
        ruling = ruling_payload()
        vertices = ruling["polygons"][0]["vertices"]
        ruling["polygons"][0]["vertices"] = list(reversed(vertices[1:] + vertices[:1]))
        comparison = tt.compare_rulings(self.result, ruling)
        self.assertTrue(comparison["all_match"])

    def test_missing_genus_is_not_zero(self):
        ruling = ruling_payload()
        ruling["polygons"][0]["genera"].pop()
        with self.assertRaisesRegex(ValueError, "every genus"):
            tt.compare_rulings(self.result, ruling)

    def test_reject_infinite_or_noninteger_totals_and_duplicate_genera(self):
        for field, value in (("phase_cardinality", "INFINITE"), ("total", 3.0),
                             ("total", True), ("total", -1), ("total", "3")):
            ruling = ruling_payload()
            ruling["polygons"][0]["genera"][0][field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                tt.compare_rulings(self.result, ruling)
        ruling = ruling_payload()
        ruling["polygons"][0]["genera"].append(deepcopy(ruling["polygons"][0]["genera"][0]))
        with self.assertRaisesRegex(ValueError, "distinct"):
            tt.compare_rulings(self.result, ruling)

    def test_comparison_requires_same_vertices(self):
        ruling = ruling_payload([[x + 1, y] for x, y in TRIANGLE])
        with self.assertRaisesRegex(ValueError, "same vertices"):
            tt.compare_rulings(self.result, ruling)

    def test_coefficient_mismatch_is_reported(self):
        ruling = ruling_payload()
        ruling["polygons"][0]["genera"][0]["total"] = 4
        comparison = tt.compare_rulings(self.result, ruling)
        self.assertFalse(comparison["all_match"])
        self.assertEqual(comparison["polygons"][0]["differences"],
                         [{"t_power": 0, "tropical": 1, "ruling": 2}])


class CommandLineTests(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, str(PACKAGE / "tropical_tree.py"), *map(str, args)],
                              capture_output=True, text=True, check=False)

    def test_examples_batch_and_output_file(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as directory:
            output = Path(directory) / "result.json"
            run = self.run_cli(PACKAGE / "examples.json", "--output", output)
            self.assertEqual(run.returncode, 0, run.stderr)
            result = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(result["schema"], "tropical-tree-counts-v1")
            self.assertFalse(result["uses_ruling_counts"])
            self.assertEqual([row["unrefined_count"] for row in result["polygons"]], [1, 5, 48])

    def test_comparison_exit_codes(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as directory:
            path = Path(directory) / "ruling.json"
            for condition, code in (("matching", 0), ("mismatch", 1), ("missing", 2)):
                ruling = ruling_payload()
                if condition == "mismatch":
                    ruling["polygons"][0]["genera"][0]["total"] = 4
                elif condition == "missing":
                    ruling["polygons"][0]["genera"].pop()
                path.write_text(json.dumps(ruling), encoding="utf-8")
                run = self.run_cli(PACKAGE / "example_polygon.json", "--compare-rulings", path)
                self.assertEqual(run.returncode, code, run.stderr)
                if condition != "missing":
                    self.assertEqual(json.loads(run.stdout)["comparison"]["all_match"], code == 0)
                else:
                    self.assertIn("every genus", run.stderr)

    def test_cli_max_ends_error(self):
        run = self.run_cli(PACKAGE / "examples.json", "--max-ends", 3)
        self.assertEqual(run.returncode, 2)
        self.assertIn("exceeding max_ends", run.stderr)


if __name__ == "__main__":
    unittest.main()
