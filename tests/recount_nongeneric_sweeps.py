#!/usr/bin/env python3
"""Recount atlas diagrams whose deterministic x-sweep has tied crossings."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "downloads" / "direct-ruling-counter"
sys.path.insert(0, str(ENGINE))

from direct_ruling_count import count_multigeodesic_rulings  # noqa: E402
from fig3_direct_ruling_counts import build_polygon_ruling_diagram  # noqa: E402
from verify_multigeodesic_rulings import (  # noqa: E402
    make_parallel_geodesic_arrangement,
)


def has_tied_times(arrangement: object) -> bool:
    times = [vertex.position[0] for vertex in arrangement.vertices]
    return len(times) != len(set(times))


def generic_covector(arrangement: object, lines: tuple[object, ...]) -> tuple[int, int]:
    """Choose (M,1) in the x-sweep chamber with distinct S1 event times."""

    for coefficient in range(1, 100_000):
        if any(
            (coefficient * line.direction[0] + line.direction[1])
            * line.direction[0] <= 0
            for line in lines
        ):
            continue
        times = [
            (coefficient * vertex.position[0] + vertex.position[1]) % 1
            for vertex in arrangement.vertices
        ]
        if len(times) == len(set(times)):
            return coefficient, 1
    raise RuntimeError("failed to find a generic covector in the x-sweep chamber")


def load_cases() -> list[dict[str, object]]:
    text = (ROOT / "data" / "cases.js").read_text(encoding="utf-8")
    prefix = "window.RULING_CASES = "
    payload = text.split(prefix, 1)[1].strip()
    return json.loads(payload[:-1])


def annular_mask_multiplicities(entries: list[dict[str, object]]) -> dict[str, int]:
    return {
        str(entry["mask"]): int(entry["accepted"])
        for entry in entries
        if int(entry["accepted"]) > 0
    }


def recorded_annular_multiplicities(
    entries: list[dict[str, object]],
) -> dict[str, int]:
    return {
        str(entry["mask"]): int(entry["multiplicity"])
        for entry in entries
        if entry["sector"] == "annular"
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("ids", nargs="*")
    parser.add_argument("--genus", type=int)
    arguments = parser.parse_args()
    selected = set(arguments.ids)
    comparisons = []
    for record in load_cases():
        label = str(record["id"])
        if selected and label not in selected:
            continue
        old_diagram = build_polygon_ruling_diagram(
            label,
            record["vertices"],
            gallery_name="sweep-recount",
        )
        old_arrangement = make_parallel_geodesic_arrangement(
            old_diagram.geodesics,
            name=label,
        )
        if not has_tied_times(old_arrangement):
            continue
        tau = generic_covector(old_arrangement, old_diagram.geodesics)
        sigma = (-1, 0)
        chamber = old_diagram.alternating_chambers[0]
        genus_results = []
        for old in record["genera"]:
            genus = int(old["genus"])
            if arguments.genus is not None and genus != arguments.genus:
                continue
            result = count_multigeodesic_rulings(
                old_diagram.geodesics,
                boundary_components=len(old_diagram.geodesics),
                target_genus=genus,
                basepoint=chamber.basepoint,
                incremental=True,
                name=f"sweep-generic-recount({label})",
                tau=tau,
                sigma=sigma,
            )
            counts = old["counts"]
            old_triplet = [counts["allDisk"], counts["annular"], counts["total"]]
            new_triplet = [result.all_disk, result.annular, result.total]
            rulings = old["rulings"]
            old_disk_masks = sorted(
                entry["mask"] for entry in rulings if entry["sector"] == "all-disk"
            )
            new_disk_masks = sorted(result.audit.get("allDiskMasks", []))
            old_annular_masks = recorded_annular_multiplicities(rulings)
            new_annular_masks = annular_mask_multiplicities(
                result.audit.get("annularResults", [])
            )
            genus_results.append(
                {
                    "genus": genus,
                    "old": old_triplet,
                    "new": new_triplet,
                    "totalEqual": old_triplet[2] == new_triplet[2],
                    "sectorCountsEqual": old_triplet == new_triplet,
                    "diskMasksEqual": old_disk_masks == new_disk_masks,
                    "annularMasksEqual": old_annular_masks == new_annular_masks,
                    "oldOnlyDiskMasks": sorted(
                        set(old_disk_masks) - set(new_disk_masks)
                    ),
                    "newOnlyDiskMasks": sorted(
                        set(new_disk_masks) - set(old_disk_masks)
                    ),
                    "oldOnlyAnnularMasks": {
                        mask: multiplicity
                        for mask, multiplicity in old_annular_masks.items()
                        if new_annular_masks.get(mask) != multiplicity
                    },
                    "newOnlyAnnularMasks": {
                        mask: multiplicity
                        for mask, multiplicity in new_annular_masks.items()
                        if old_annular_masks.get(mask) != multiplicity
                    },
                }
            )
        comparisons.append(
            {
                "id": label,
                "originAttempt": old_diagram.origin_attempt,
                "genericTau": list(tau),
                "genericSigma": list(sigma),
                "genera": genus_results,
            }
        )
        print(
            json.dumps(comparisons[-1], ensure_ascii=False),
            file=sys.stderr,
            flush=True,
        )

    output = {
        "diagramsRecounted": len(comparisons),
        "allTotalsEqual": all(
            genus["totalEqual"]
            for comparison in comparisons
            for genus in comparison["genera"]
        ),
        "allSectorCountsEqual": all(
            genus["sectorCountsEqual"]
            for comparison in comparisons
            for genus in comparison["genera"]
        ),
        "allCertificateSetsEqual": all(
            genus["diskMasksEqual"] and genus["annularMasksEqual"]
            for comparison in comparisons
            for genus in comparison["genera"]
        ),
        "comparisons": comparisons,
    }
    print(json.dumps(output, indent=2, ensure_ascii=False))
    return 0 if output["allTotalsEqual"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
