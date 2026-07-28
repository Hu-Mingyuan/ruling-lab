#!/usr/bin/env python3
"""Export the one-interior-point ruling atlas as browser-ready JavaScript."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
SITE_ROOT = HERE.parent
PROJECT_ROOT = SITE_ROOT.parent
DEFAULT_POLYGONS = PROJECT_ROOT / "one_interior_lattice_polygons.json"
DEFAULT_COUNTS = PROJECT_ROOT / "one_interior_ruling_counts.json"
DEFAULT_MANIFEST = (
    SITE_ROOT / "assets" / "rulings" / "one-interior" / "manifest.json"
)
DEFAULT_OUTPUT = SITE_ROOT / "data" / "cases.js"


def load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def arrangement_shear(vertices: object) -> int:
    if not isinstance(vertices, list) or len(vertices) < 3:
        raise ValueError("polygon vertices must be a list")
    polygon = [
        (int(point[0]), int(point[1]))
        for point in vertices
        if isinstance(point, list) and len(point) == 2
    ]
    if len(polygon) != len(vertices):
        raise ValueError("malformed polygon vertex")
    directions: list[tuple[int, int]] = []
    for start, end in zip(polygon, polygon[1:] + polygon[:1]):
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        length = math.gcd(abs(dx), abs(dy))
        if length == 0:
            raise ValueError("polygon has a repeated vertex")
        directions.extend([(dx // length, dy // length)] * length)
    for magnitude in range(0, 1 + len(directions)):
        candidates = (0,) if magnitude == 0 else (magnitude, -magnitude)
        for shear in candidates:
            if all(dx + shear * dy != 0 for dx, dy in directions):
                return shear
    raise AssertionError("failed to find a transverse integral shear")


def browser_drawings(
    identifier: str,
    genus_number: int,
    manifest_records: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    record = manifest_records.get(identifier)
    if record is None:
        return []
    raw_drawings = record.get("drawings", [])
    if not isinstance(raw_drawings, list):
        raise ValueError(f"invalid drawings for {identifier}")

    output: list[dict[str, object]] = []
    sector_indices: dict[str, int] = {}
    for raw in raw_drawings:
        if not isinstance(raw, dict) or int(raw.get("genus", -1)) != genus_number:
            continue
        sector = (
            "annular"
            if str(raw["sector"]).lower().startswith("annular")
            else "all-disk"
        )
        sector_indices[sector] = sector_indices.get(sector, 0) + 1
        if sector == "annular":
            identifier_text = f"A{sector_indices[sector]:02d}"
        else:
            identifier_text = f"D{sector_indices[sector]:02d}"
        relative_file = Path(str(raw["file"])).as_posix()
        disk_eyes = int(raw["disk_eyes"])
        annular_eyes = int(raw["annular_eyes"])
        output.append(
            {
                "identifier": identifier_text,
                "sector": sector,
                "src": (
                    "assets/rulings/one-interior/" + relative_file
                ),
                "mask": raw["mask"],
                "witness": raw.get("witness"),
                "diskEyes": disk_eyes,
                "annularEyes": annular_eyes,
                "switches": int(raw["switches"]),
                "sl2z": raw["sl2z"],
                "profile": f"D={disk_eyes} · A={annular_eyes}",
                "multiplicity": 1,
            }
        )
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--polygons", type=Path, default=DEFAULT_POLYGONS)
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()

    polygon_payload = load_json(arguments.polygons)
    count_payload = load_json(arguments.counts)
    source_records = {
        str(record["id"]): record
        for record in polygon_payload["polygons"]  # type: ignore[index]
    }
    count_records = {
        str(record["id"]): record
        for record in count_payload["polygons"]  # type: ignore[index]
    }
    if source_records.keys() != count_records.keys() or len(source_records) != 16:
        raise AssertionError("polygon and count files must contain the same 16 ids")

    manifest: dict[str, object] = {}
    if arguments.manifest.is_file():
        manifest = load_json(arguments.manifest)
    if manifest and manifest.get("uses_tropical_count") is not False:
        raise ValueError("drawing manifest must record uses_tropical_count=false")
    raw_manifest_polygons = manifest.get("polygons", [])
    if not isinstance(raw_manifest_polygons, list):
        raise ValueError("manifest 'polygons' must be a list")
    manifest_records: dict[str, dict[str, object]] = {}
    for record in raw_manifest_polygons:
        if not isinstance(record, dict):
            raise ValueError("malformed drawing manifest polygon")
        manifest_records[str(record["id"])] = record
    if manifest_records and manifest_records.keys() != source_records.keys():
        raise AssertionError("drawing manifest must contain the same 16 ids")

    cases: list[dict[str, object]] = []
    for identifier, source in source_records.items():
        counted = count_records[identifier]
        manifest_record = manifest_records.get(identifier, {})
        coordinate_frames = manifest_record.get("coordinate_frames", {})
        if coordinate_frames and not isinstance(coordinate_frames, dict):
            raise ValueError(f"invalid coordinate frames for {identifier}")
        if coordinate_frames:
            count_from_source = coordinate_frames["count_from_source"]
            display_from_count = coordinate_frames["display_from_count"]
            display_from_source = coordinate_frames["display_from_source"]
            vertical_direction = coordinate_frames["vertical_direction"]
            if (
                not isinstance(count_from_source, list)
                or count_from_source[0][0] != 1
                or count_from_source[1] != [0, 1]
            ):
                raise ValueError(f"invalid counting shear for {identifier}")
            counting_shear = int(count_from_source[0][1])
        else:
            counting_shear = arrangement_shear(source["vertices"])
            count_from_source = [[1, counting_shear], [0, 1]]
            display_from_count = [[1, 0], [0, 1]]
            display_from_source = count_from_source
            vertical_direction = [0, 1]
        genera = []
        for genus in counted["genera"]:  # type: ignore[index]
            genus_number = int(genus["genus"])
            drawings = browser_drawings(
                identifier, genus_number, manifest_records
            )
            display_counts = genus
            manifest_counts = manifest_record.get("counts")
            if isinstance(manifest_counts, dict):
                candidate = manifest_counts.get(f"genus_{genus_number}")
                if isinstance(candidate, dict):
                    if int(candidate["total"]) != int(genus["total"]):
                        raise AssertionError(
                            f"{identifier} genus {genus_number}: "
                            "display and direct totals disagree"
                        )
                    display_counts = candidate
            genera.append(
                {
                    "genus": genus_number,
                    "counts": {
                        "allDisk": display_counts["all_disk"],
                        "annular": display_counts["annular"],
                        "total": display_counts["total"],
                    },
                    "phaseCardinality": genus["phase_cardinality"],
                    "rulings": drawings,
                }
            )
        cases.append(
            {
                "id": identifier,
                "vertices": source["vertices"],
                "arrangementShear": counting_shear,
                "displaySl2z": display_from_count,
                "displayFromSource": display_from_source,
                "verticalDirection": vertical_direction,
                "displayName": manifest_record.get("display_name"),
                "lambdaModel": manifest_record.get("lambda_model"),
                "doubleArea": source["double_area"],
                "boundaryLatticePoints": source["boundary_lattice_points"],
                "interiorLatticePoints": len(source["interior_lattice_points"]),
                "genera": genera,
                "provenance": {
                    "method": count_payload["method"],
                    "usesTropicalCount": False,
                },
            }
        )

    cases.sort(
        key=lambda item: (
            len(item["vertices"]),
            item["boundaryLatticePoints"],
            item["vertices"],
        )
    )
    for item in cases:
        for genus in item["genera"]:
            counts = genus["counts"]
            if counts["total"] != counts["allDisk"] + counts["annular"]:
                raise AssertionError(
                    f"{item['id']} genus {genus['genus']}: sectors do not sum"
                )
            drawings = genus["rulings"]
            if drawings:
                if any(
                    drawing["sl2z"] != item["displaySl2z"]
                    for drawing in drawings
                ):
                    raise AssertionError(
                        f"{item['id']}: drawings do not share one display frame"
                    )
                multiplicity = sum(
                    int(drawing.get("multiplicity", 1)) for drawing in drawings
                )
                if multiplicity != int(counts["total"]):
                    raise AssertionError(
                        f"{item['id']} genus {genus['genus']}: "
                        "drawing multiplicities do not sum"
                    )

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        "/* Generated by tools/export_cases.py. */\n"
        "window.RULING_CASES = "
        + json.dumps(cases, indent=2, ensure_ascii=False)
        + ";\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "cases": len(cases),
                "genera": sum(len(case["genera"]) for case in cases),
                "drawings": sum(
                    len(genus["rulings"])
                    for case in cases
                    for genus in case["genera"]
                ),
                "usesTropicalCount": False,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
