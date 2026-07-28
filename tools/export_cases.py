#!/usr/bin/env python3
"""Export both interior-point ruling collections as browser-ready JavaScript."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


HERE = Path(__file__).resolve().parent
SITE_ROOT = HERE.parent
PROJECT_ROOT = SITE_ROOT.parent
DEFAULT_OUTPUT = SITE_ROOT / "data" / "cases.js"
EXPECTED_COLLECTION_MULTIPLICITY = {
    "interior-1": 112,
    "interior-2": 1489,
}


@dataclass(frozen=True)
class Collection:
    key: str
    interior_points: int
    expected_polygons: int
    asset_directory: str
    polygons: Path
    counts: Path
    manifest: Path


def default_collections() -> tuple[Collection, Collection]:
    return (
        Collection(
            key="interior-1",
            interior_points=1,
            expected_polygons=16,
            asset_directory="one-interior",
            polygons=PROJECT_ROOT / "one_interior_lattice_polygons.json",
            counts=PROJECT_ROOT / "one_interior_ruling_counts.json",
            manifest=(
                SITE_ROOT
                / "assets"
                / "rulings"
                / "one-interior"
                / "manifest.json"
            ),
        ),
        Collection(
            key="interior-2",
            interior_points=2,
            expected_polygons=45,
            asset_directory="two-interior",
            polygons=PROJECT_ROOT / "two_interior_lattice_polygons.json",
            counts=PROJECT_ROOT / "two_interior_symmetric_ruling_counts.json",
            manifest=(
                SITE_ROOT
                / "assets"
                / "rulings"
                / "two-interior"
                / "manifest.json"
            ),
        ),
    )


def load_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def records_by_id(
    payload: dict[str, object], path: Path
) -> dict[str, dict[str, object]]:
    raw_records = payload.get("polygons")
    if not isinstance(raw_records, list):
        raise ValueError(f"{path}: 'polygons' must be a list")
    records: dict[str, dict[str, object]] = {}
    for raw_record in raw_records:
        if not isinstance(raw_record, dict) or "id" not in raw_record:
            raise ValueError(f"{path}: malformed polygon record")
        identifier = str(raw_record["id"])
        if identifier in records:
            raise ValueError(f"{path}: duplicate id {identifier}")
        records[identifier] = raw_record
    return records


def load_manifest_records(
    collection: Collection,
) -> dict[str, dict[str, object]]:
    if not collection.manifest.is_file():
        raise FileNotFoundError(
            f"{collection.key}: drawing manifest is not ready: "
            f"{collection.manifest}"
        )
    payload = load_json(collection.manifest)
    if payload.get("uses_tropical_count") is not False:
        raise ValueError(
            f"{collection.manifest}: uses_tropical_count must be false"
        )
    return records_by_id(payload, collection.manifest)


def integer_vector(value: object, label: str) -> list[int]:
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(f"{label} must be a two-component vector")
    try:
        return [int(value[0]), int(value[1])]
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be an integer vector") from error


def integer_matrix(value: object, label: str) -> list[list[int]]:
    if (
        not isinstance(value, list)
        or len(value) != 2
        or not all(isinstance(row, list) and len(row) == 2 for row in value)
    ):
        raise ValueError(f"{label} must be a 2 by 2 matrix")
    try:
        return [
            [int(value[0][0]), int(value[0][1])],
            [int(value[1][0]), int(value[1][1])],
        ]
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be an integer matrix") from error


def determinant(matrix: list[list[int]]) -> int:
    return (
        matrix[0][0] * matrix[1][1]
        - matrix[0][1] * matrix[1][0]
    )


def matrix_product(
    left: list[list[int]], right: list[list[int]]
) -> list[list[int]]:
    return [
        [
            left[row][0] * right[0][column]
            + left[row][1] * right[1][column]
            for column in range(2)
        ]
        for row in range(2)
    ]


def coordinate_frame_fields(
    identifier: str, manifest_record: dict[str, object]
) -> dict[str, object]:
    raw_frames = manifest_record.get("coordinate_frames")
    if not isinstance(raw_frames, dict):
        raise ValueError(
            f"{identifier}: manifest coordinate_frames are required"
        )

    count_from_source = integer_matrix(
        raw_frames.get("count_from_source"),
        f"{identifier}: count_from_source",
    )
    display_from_count = integer_matrix(
        raw_frames.get("display_from_count"),
        f"{identifier}: display_from_count",
    )
    display_from_source = integer_matrix(
        raw_frames.get("display_from_source"),
        f"{identifier}: display_from_source",
    )
    vertical_direction_source = integer_vector(
        raw_frames.get("vertical_direction_source"),
        f"{identifier}: vertical_direction_source",
    )
    vertical_direction_count = integer_vector(
        raw_frames.get("vertical_direction_count"),
        f"{identifier}: vertical_direction_count",
    )
    vertical_direction = integer_vector(
        raw_frames.get("vertical_direction"),
        f"{identifier}: vertical_direction",
    )
    sweep_covector_source = integer_vector(
        raw_frames.get("sweep_covector_source"),
        f"{identifier}: sweep_covector_source",
    )
    sweep_covector_count = integer_vector(
        raw_frames.get("sweep_covector_count"),
        f"{identifier}: sweep_covector_count",
    )

    if (
        count_from_source[0][0] != 1
        or count_from_source[1] != [0, 1]
    ):
        raise ValueError(f"{identifier}: invalid counting shear")
    counting_shear = count_from_source[0][1]
    if determinant(display_from_count) != 1:
        raise ValueError(
            f"{identifier}: display_from_count must lie in SL(2,Z)"
        )
    if display_from_source != matrix_product(
        display_from_count, count_from_source
    ):
        raise ValueError(
            f"{identifier}: display_from_source is not the composed frame"
        )
    if vertical_direction_source != [-counting_shear, 1]:
        raise ValueError(
            f"{identifier}: invalid source vertical direction"
        )
    if vertical_direction_count != [0, 1]:
        raise ValueError(
            f"{identifier}: invalid count vertical direction"
        )
    expected_vertical_direction = [
        display_from_count[0][1],
        display_from_count[1][1],
    ]
    if vertical_direction != expected_vertical_direction:
        raise ValueError(
            f"{identifier}: vertical direction was not transported"
        )
    if sweep_covector_source != [1, counting_shear]:
        raise ValueError(f"{identifier}: invalid source sweep covector")
    if sweep_covector_count != [1, 0]:
        raise ValueError(f"{identifier}: invalid count sweep covector")

    return {
        "arrangementShear": counting_shear,
        "displaySl2z": display_from_count,
        "displayFromSource": display_from_source,
        "sourceVerticalDirection": vertical_direction_source,
        "countVerticalDirection": vertical_direction_count,
        "verticalDirection": vertical_direction,
        "sourceSweepCovector": sweep_covector_source,
        "countSweepCovector": sweep_covector_count,
    }


def browser_drawings(
    identifier: str,
    genus_number: int,
    manifest_record: dict[str, object],
    asset_directory: str,
) -> list[dict[str, object]]:
    raw_drawings = manifest_record.get("drawings")
    if not isinstance(raw_drawings, list):
        raise ValueError(f"{identifier}: manifest drawings must be a list")

    output: list[dict[str, object]] = []
    sector_indices: dict[str, int] = {}
    for raw in raw_drawings:
        if not isinstance(raw, dict):
            raise ValueError(f"{identifier}: malformed drawing record")
        if int(raw.get("genus", -1)) != genus_number:
            continue

        sector = (
            "annular"
            if str(raw.get("sector", "")).lower().startswith("annular")
            else "all-disk"
        )
        sector_indices[sector] = sector_indices.get(sector, 0) + 1
        prefix = "A" if sector == "annular" else "D"
        identifier_text = f"{prefix}{sector_indices[sector]:02d}"
        relative_file = Path(str(raw["file"])).as_posix()
        disk_eyes = int(raw["disk_eyes"])
        annular_eyes = int(raw["annular_eyes"])
        multiplicity = int(raw.get("multiplicity", 1))
        if multiplicity <= 0:
            raise ValueError(
                f"{identifier} genus {genus_number}: "
                "drawing multiplicity must be positive"
            )
        sl2z = integer_matrix(
            raw.get("sl2z"),
            f"{identifier} genus {genus_number}: drawing sl2z",
        )
        if determinant(sl2z) != 1:
            raise ValueError(
                f"{identifier} genus {genus_number}: "
                "drawing sl2z must lie in SL(2,Z)"
            )
        output.append(
            {
                "identifier": identifier_text,
                "sector": sector,
                "src": (
                    f"assets/rulings/{asset_directory}/{relative_file}"
                ),
                "mask": raw.get("mask"),
                "witness": raw.get("witness"),
                "diskEyes": disk_eyes,
                "annularEyes": annular_eyes,
                "switches": int(raw["switches"]),
                "sl2z": sl2z,
                "profile": f"D={disk_eyes} · A={annular_eyes}",
                "multiplicity": multiplicity,
            }
        )
    return output


def export_collection(collection: Collection) -> list[dict[str, object]]:
    polygon_payload = load_json(collection.polygons)
    count_payload = load_json(collection.counts)
    if count_payload.get("uses_tropical_count") is not False:
        raise ValueError(
            f"{collection.counts}: uses_tropical_count must be false"
        )

    source_records = records_by_id(polygon_payload, collection.polygons)
    count_records = records_by_id(count_payload, collection.counts)
    expected = collection.expected_polygons
    if (
        source_records.keys() != count_records.keys()
        or len(source_records) != expected
    ):
        raise AssertionError(
            f"{collection.key}: polygon and count files must contain "
            f"the same {expected} ids"
        )

    manifest_records = load_manifest_records(collection)
    if (
        manifest_records.keys() != source_records.keys()
        or len(manifest_records) != expected
    ):
        raise AssertionError(
            f"{collection.key}: drawing manifest must contain the same "
            f"{expected} ids"
        )

    cases: list[dict[str, object]] = []
    for identifier, source in source_records.items():
        counted = count_records[identifier]
        manifest_record = manifest_records[identifier]
        frame_fields = coordinate_frame_fields(identifier, manifest_record)

        stored_interior = source.get("interior_lattice_points")
        interior_count = (
            len(stored_interior)
            if isinstance(stored_interior, list)
            else int(stored_interior)
        )
        if interior_count != collection.interior_points:
            raise AssertionError(
                f"{identifier}: expected {collection.interior_points} "
                "interior lattice points"
            )

        raw_genera = counted.get("genera")
        if not isinstance(raw_genera, list):
            raise ValueError(f"{identifier}: genera must be a list")
        manifest_counts = manifest_record.get("counts")
        if not isinstance(manifest_counts, dict):
            raise ValueError(f"{identifier}: manifest counts are missing")

        genera: list[dict[str, object]] = []
        for raw_genus in raw_genera:
            if not isinstance(raw_genus, dict):
                raise ValueError(f"{identifier}: malformed genus record")
            genus_number = int(raw_genus["genus"])
            recorded = manifest_counts.get(f"genus_{genus_number}")
            if not isinstance(recorded, dict):
                raise ValueError(
                    f"{identifier}: manifest genus {genus_number} is missing"
                )

            direct_total = int(raw_genus["total"])
            display_counts = {
                "allDisk": int(recorded["all_disk"]),
                "annular": int(recorded["annular"]),
                "total": int(recorded["total"]),
            }
            if display_counts["total"] != direct_total:
                raise AssertionError(
                    f"{identifier} genus {genus_number}: "
                    "manifest and direct totals disagree"
                )
            if (
                display_counts["allDisk"] + display_counts["annular"]
                != display_counts["total"]
            ):
                raise AssertionError(
                    f"{identifier} genus {genus_number}: "
                    "manifest sectors do not sum"
                )

            drawings = browser_drawings(
                identifier,
                genus_number,
                manifest_record,
                collection.asset_directory,
            )
            if any(
                drawing["sl2z"] != frame_fields["displaySl2z"]
                for drawing in drawings
            ):
                raise AssertionError(
                    f"{identifier}: drawings do not share one display frame"
                )

            weighted_by_sector = {"all-disk": 0, "annular": 0}
            for drawing in drawings:
                weighted_by_sector[str(drawing["sector"])] += int(
                    drawing["multiplicity"]
                )
            if (
                weighted_by_sector["all-disk"]
                != display_counts["allDisk"]
                or weighted_by_sector["annular"]
                != display_counts["annular"]
            ):
                raise AssertionError(
                    f"{identifier} genus {genus_number}: "
                    "drawing multiplicities disagree with sector counts"
                )
            if (
                sum(weighted_by_sector.values())
                != display_counts["total"]
            ):
                raise AssertionError(
                    f"{identifier} genus {genus_number}: "
                    "drawing multiplicities do not sum"
                )

            genera.append(
                {
                    "genus": genus_number,
                    "counts": display_counts,
                    "phaseCardinality": raw_genus["phase_cardinality"],
                    "rulings": drawings,
                }
            )

        expected_genera = list(range(collection.interior_points + 1))
        if sorted(int(entry["genus"]) for entry in genera) != expected_genera:
            raise AssertionError(
                f"{identifier}: expected genera {expected_genera}"
            )

        cases.append(
            {
                "id": identifier,
                "collection": collection.key,
                "vertices": source["vertices"],
                **frame_fields,
                "fixedDisplayFrame": True,
                "displayName": manifest_record.get("display_name"),
                "lambdaModel": manifest_record.get("lambda_model"),
                "doubleArea": int(source["double_area"]),
                "boundaryLatticePoints": int(
                    source["boundary_lattice_points"]
                ),
                "interiorLatticePoints": interior_count,
                "genera": genera,
                "drawingsReady": True,
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
    collection_multiplicity = sum(
        int(drawing["multiplicity"])
        for item in cases
        for genus in item["genera"]
        for drawing in genus["rulings"]
    )
    expected_multiplicity = EXPECTED_COLLECTION_MULTIPLICITY[collection.key]
    if collection_multiplicity != expected_multiplicity:
        raise AssertionError(
            f"{collection.key}: expected weighted multiplicity "
            f"{expected_multiplicity}, got {collection_multiplicity}"
        )
    return cases


def configured_collections(arguments: argparse.Namespace) -> Sequence[Collection]:
    defaults = default_collections()
    return (
        Collection(
            **{
                **defaults[0].__dict__,
                "polygons": arguments.one_polygons,
                "counts": arguments.one_counts,
                "manifest": arguments.one_manifest,
            }
        ),
        Collection(
            **{
                **defaults[1].__dict__,
                "polygons": arguments.two_polygons,
                "counts": arguments.two_counts,
                "manifest": arguments.two_manifest,
            }
        ),
    )


def main() -> int:
    one, two = default_collections()
    parser = argparse.ArgumentParser()
    parser.add_argument("--one-polygons", type=Path, default=one.polygons)
    parser.add_argument("--one-counts", type=Path, default=one.counts)
    parser.add_argument("--one-manifest", type=Path, default=one.manifest)
    parser.add_argument("--two-polygons", type=Path, default=two.polygons)
    parser.add_argument("--two-counts", type=Path, default=two.counts)
    parser.add_argument("--two-manifest", type=Path, default=two.manifest)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()

    try:
        cases = [
            item
            for collection in configured_collections(arguments)
            for item in export_collection(collection)
        ]
    except FileNotFoundError as error:
        parser.error(str(error))
    if len(cases) != 61 or len({str(item["id"]) for item in cases}) != 61:
        raise AssertionError("combined atlas must contain 61 unique polygons")

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        "/* Generated by tools/export_cases.py. */\n"
        "window.RULING_CASES = "
        + json.dumps(cases, indent=2, ensure_ascii=False)
        + ";\n",
        encoding="utf-8",
    )
    summary = {
        "cases": len(cases),
        "collections": {
            key: sum(item["collection"] == key for item in cases)
            for key in ("interior-1", "interior-2")
        },
        "genera": sum(len(case["genera"]) for case in cases),
        "drawings": sum(
            len(genus["rulings"])
            for case in cases
            for genus in case["genera"]
        ),
        "ruling_multiplicity": sum(
            int(drawing["multiplicity"])
            for case in cases
            for genus in case["genera"]
            for drawing in genus["rulings"]
        ),
        "usesTropicalCount": False,
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
