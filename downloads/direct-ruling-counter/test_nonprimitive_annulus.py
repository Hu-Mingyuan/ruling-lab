"""Definition-level annular witnesses whose full period must be retained."""

import json
from pathlib import Path
import subprocess
import unittest

from direct_ruling_count import _node_executable
from fig3_direct_ruling_counts import build_polygon_ruling_diagram
from multigeodesic_annular_input import exact_annular_payload


class FullAnnularPeriodTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        diagram = build_polygon_ruling_diagram(
            "validation-S7", [[0, 5], [1, 2], [4, 0], [3, 3]]
        )
        cls.payload = exact_annular_payload(
            diagram.geodesics, boundary_components=4, target_genus=1,
            basepoint=diagram.alternating_chambers[0].basepoint,
            arrangement=diagram.arrangement,
        )
        # These are directly checked smoothing witnesses, not values inferred
        # from a tropical polynomial.  The first has period (2,-2); the second
        # has primitive period (1,-1).
        cls.payload["masks"] = ["0x4303080", "0x4330800"]
        cls.engine = Path(__file__).parent / "tmp/generic_annular_exact_search.js"

    def audit(self, tau, sigma, certificates):
        payload = {**self.payload, "tau": tau, "sigma": sigma,
                   "collectAnswers": certificates, "countOnly": not certificates}
        run = subprocess.run(
            [str(_node_executable()), str(self.engine)], input=json.dumps(payload),
            text=True, capture_output=True, check=True, timeout=30,
        )
        return json.loads(run.stdout)

    def test_full_period_and_primitive_witnesses_in_both_sweep_directions(self):
        for tau, sigma in (([1, 0], [0, 1]), ([4, 1], [1, 0]), ([-1, 0], [0, 1])):
            with self.subTest(tau=tau):
                result = self.audit(tau, sigma, True)
                self.assertEqual(result["phaseCardinality"], "FINITE")
                self.assertEqual(result["exactTotalRulingCount"], 2)
                witnesses = [entry["witnesses"][0] for entry in result["results"]]
                self.assertEqual([w["annularHomologies"] for w in witnesses],
                                 [[[2, -2]], [[1, -1]]])
                for witness in witnesses:
                    self.assertEqual(witness["baseBalance"], 0)
                    self.assertEqual(witness["ribbonGenus"], 1)
                    self.assertEqual(witness["ribbonBoundaryComponents"], 4)
                    self.assertEqual(witness["typeBCompanionOrders"],
                                     {"same": 2, "opposite": 0, "inconsistent": 0})
                self.assertEqual(witnesses[0]["annularPhaseShifts"], [[0, 0]])

    def test_count_only_preserves_nonprimitive_witness_count(self):
        result = self.audit([1, 0], [0, 1], False)
        self.assertEqual(result["exactTotalRulingCount"], 2)
        self.assertEqual(result["results"], [])


if __name__ == "__main__":
    unittest.main()
