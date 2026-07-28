#!/usr/bin/env python3
"""Refresh and reproducibly archive the reader-facing direct ruling engine."""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path
import re
import shutil
import zipfile


SITE_ROOT = Path(__file__).resolve().parent.parent
PROJECT_ROOT = SITE_ROOT.parent
DOWNLOADS = SITE_ROOT / "downloads"
PACKAGE_NAME = "direct-ruling-counter"
PACKAGE_ROOT = DOWNLOADS / PACKAGE_NAME
ARCHIVE = DOWNLOADS / f"{PACKAGE_NAME}.zip"

ENGINE_FILES = (
    "ruling_polygon.py",
    "fig3_direct_ruling_counts.py",
    "direct_ruling_count.py",
    "multigeodesic_annular_input.py",
    "multigeodesic_alternating_audit.py",
    "verify_multigeodesic_rulings.py",
    "verify_triangle_rulings.py",
    "ruling_smoothing.py",
    "tmp/generic_incremental_ruling_search.js",
    "tmp/generic_annular_exact_search.js",
    "tmp/p2_annular_exact_search.js",
    "tmp/p2o4_annular_search.js",
    "tmp/p2o4_smoothing_search.js",
    "tmp/annular_phase_exact.js",
)

SUPPORT_FILES = (
    "README.md",
    "example_polygon.json",
    "smoke_test.py",
)

EXCLUDED_BASENAMES = {
    "tropical_count.py",
    "blomme_recursion.py",
    "one_interior_ruling_counts.json",
    "fig2_direct_ruling_counts.json",
    "fig3_direct_ruling_counts.json",
    "fig10_direct_ruling_counts.json",
    "fig12_direct_ruling_counts.json",
}

JS_LOCAL_REQUIRE = re.compile(r"""require\(["']\./([^"']+\.js)["']\)""")

DIRECT_COPY_FILES = (
    "ruling_polygon.py",
    "tmp/generic_incremental_ruling_search.js",
    "tmp/generic_annular_exact_search.js",
    "tmp/p2_annular_exact_search.js",
    "tmp/p2o4_annular_search.js",
    "tmp/p2o4_smoothing_search.js",
    "tmp/annular_phase_exact.js",
)

PRUNED_MODULES = {
    "fig3_direct_ruling_counts.py": {
        "header": """#!/usr/bin/env python3
\"\"\"Build and count a direct ruling diagram from polygon vertex coordinates.\"\"\"

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from fractions import Fraction
from math import gcd
import json
from pathlib import Path
import random
from typing import Iterable, Sequence

from direct_ruling_count import count_multigeodesic_rulings
from multigeodesic_alternating_audit import TorusChamber, enumerate_torus_chambers
from verify_multigeodesic_rulings import make_parallel_geodesic_arrangement
from verify_triangle_rulings import Geodesic

PointI = tuple[int, int]
Direction = tuple[int, int]
ORIGIN_DENOMINATOR = 1009
COUNT_VERTICAL_DIRECTION: Direction = (0, 1)
COUNT_SWEEP_COVECTOR: Direction = (1, 0)
""",
        "definitions": (
            "PolygonRulingDiagram",
            "_primitive_polygon_edges",
            "_transverse_shear",
            "build_polygon_ruling_diagram",
            "_fraction_text",
            "count_case",
            "_load_records",
        ),
        "replacements": (
            ('gallery_name: str = "Fig. 3"', 'gallery_name: str = "lattice polygon"'),
        ),
        "footer": """
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "labels",
        nargs="*",
        help="optional polygon ids; default: every record",
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--genus", type=int, default=0)
    parser.add_argument("--flat-reference", action="store_true")
    arguments = parser.parse_args()
    if arguments.jobs < 1:
        parser.error("--jobs must be positive")
    if arguments.genus < 0:
        parser.error("--genus must be nonnegative")
    records = _load_records(arguments.source)
    if arguments.labels:
        wanted = set(arguments.labels)
        records = [
            record
            for record in records
            if str(record.get("id", record.get("label"))) in wanted
        ]
        missing = wanted - {
            str(record.get("id", record.get("label"))) for record in records
        }
        if missing:
            parser.error("unknown labels: " + ", ".join(sorted(missing)))

    incremental = not arguments.flat_reference
    gallery_name = arguments.source.stem
    if arguments.jobs == 1:
        for record in records:
            print(
                json.dumps(
                    count_case(
                        record,
                        gallery_name=gallery_name,
                        target_genus=arguments.genus,
                        incremental=incremental,
                    ),
                    separators=(",", ":"),
                ),
                flush=True,
            )
        return 0

    with ThreadPoolExecutor(max_workers=arguments.jobs) as executor:
        futures = {
            executor.submit(
                count_case,
                record,
                gallery_name=gallery_name,
                target_genus=arguments.genus,
                incremental=incremental,
            ): record
            for record in records
        }
        for future in as_completed(futures):
            print(json.dumps(future.result(), separators=(",", ":")), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
""",
    },
    "direct_ruling_count.py": {
        "header": """#!/usr/bin/env python3
\"\"\"Run the exact direct disk/annular ruling search.\"\"\"

from __future__ import annotations

import json
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from multigeodesic_annular_input import exact_annular_payload
from verify_triangle_rulings import Geodesic, Point
""",
        "definitions": (
            "DirectRulingCount",
            "_node_executable",
            "count_multigeodesic_rulings",
        ),
    },
    "multigeodesic_annular_input.py": {
        "header": """#!/usr/bin/env python3
\"\"\"Serialize a primitive multigeodesic diagram for the exact annular solver.\"\"\"

from __future__ import annotations

from fractions import Fraction
from math import gcd, lcm
from typing import Iterable, Sequence

import ruling_smoothing as rs
from verify_multigeodesic_rulings import make_parallel_geodesic_arrangement
from verify_triangle_rulings import Geodesic, Point
""",
        "definitions": (
            "_fraction_text",
            "exact_annular_payload",
        ),
        "replacements": (
            (
                "# also handles the pairwise-nonparallel special case used by Fig. 3(a).",
                "# also handles the pairwise-nonparallel special case.",
            ),
        ),
    },
    "multigeodesic_alternating_audit.py": {
        "header": """#!/usr/bin/env python3
\"\"\"Construct exact alternating chambers for primitive torus geodesics.\"\"\"

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from math import floor
from typing import Iterable

from verify_triangle_rulings import Geodesic, Point, Vector
import verify_multigeodesic_rulings as vmr

Descriptor = tuple[int, str, int]
""",
        "definitions": (
            "determinant",
            "line_value",
            "level_intersection",
            "strip_bounds",
            "_clip_halfplane",
            "_twice_area",
            "IndependentCell",
            "independent_lifted_cell",
            "_local_sector_basepoints",
            "_canonical_bounds_and_translation",
            "_deck_equivalent_bounds",
            "_fractional_point",
            "TorusChamber",
            "enumerate_torus_chambers",
        ),
    },
    "verify_multigeodesic_rulings.py": {
        "header": """#!/usr/bin/env python3
\"\"\"Build exact torus arrangements from primitive geodesics.\"\"\"

from __future__ import annotations

from fractions import Fraction
from functools import cmp_to_key
from itertools import combinations
from math import gcd
from typing import Iterable

import ruling_smoothing as rs
from verify_triangle_rulings import (
    Geodesic,
    Point,
    Vector,
    _det,
    _intersection_parameters,
    _ray_compare,
)

COUNT_VERTICAL_DIRECTION: Vector = (0, 1)
COUNT_SWEEP_COVECTOR: Vector = (1, 0)
""",
        "definitions": (
            "_sweep_value",
            "_make_geodesic_arrangement",
            "make_geodesic_arrangement",
            "make_parallel_geodesic_arrangement",
        ),
    },
    "verify_triangle_rulings.py": {
        "header": """#!/usr/bin/env python3
\"\"\"Exact primitive-geodesic types and intersection geometry.\"\"\"

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import floor, gcd

Point = tuple[Fraction, Fraction]
Vector = tuple[int, int]
""",
        "definitions": (
            "_det",
            "_mod_one",
            "_mod_point",
            "Geodesic",
            "_intersection_parameters",
            "_ray_compare",
        ),
    },
    "ruling_smoothing.py": {
        "header": """#!/usr/bin/env python3
\"\"\"Minimal exact half-edge geometry used by the direct ruling counter.\"\"\"

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction

DeckVector = tuple[int, int]
Point = tuple[Fraction, Fraction]
""",
        "definitions": (
            "_neg",
            "TorusHalfEdge",
            "TorusVertex",
            "TorusArrangement",
            "_switch_candidates",
            "nonswitch_pairs",
            "switch_pairs",
            "switch_orientation",
        ),
    },
}


def _definition_source(source: str, node: ast.AST) -> str:
    decorators = getattr(node, "decorator_list", ())
    start = min(
        [node.lineno, *(decorator.lineno for decorator in decorators)]  # type: ignore[attr-defined]
    )
    end = node.end_lineno  # type: ignore[attr-defined]
    return "\n".join(source.splitlines()[start - 1:end])


def _write_pruned_module(relative: str, specification: dict[str, object]) -> None:
    source_path = PROJECT_ROOT / relative
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=relative)
    requested = tuple(specification["definitions"])
    nodes = {
        node.name: node
        for node in tree.body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
    }
    missing = [name for name in requested if name not in nodes]
    if missing:
        raise RuntimeError(
            f"cannot export {relative}; definitions disappeared: "
            + ", ".join(missing)
        )
    sections = [
        str(specification["header"]).rstrip(),
        *(_definition_source(source, nodes[name]) for name in requested),
    ]
    footer = str(specification.get("footer", "")).strip()
    if footer:
        sections.append(footer)
    output = "\n\n\n".join(sections) + "\n"
    for old, new in specification.get("replacements", ()):
        output = output.replace(old, new)
    destination = PACKAGE_ROOT / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(output, encoding="utf-8")


def refresh_sources() -> None:
    """Export only the audited live-computation closure."""

    PACKAGE_ROOT.mkdir(parents=True, exist_ok=True)
    for relative in DIRECT_COPY_FILES:
        source = PROJECT_ROOT / relative
        destination = PACKAGE_ROOT / relative
        if not source.is_file():
            raise FileNotFoundError(f"missing engine source: {source}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        if relative == "tmp/p2_annular_exact_search.js":
            text = destination.read_text(encoding="utf-8")
            text = text.replace(
                "The old finite-K audit is useful as a regression oracle",
                "The old finite-K audit is useful as a regression cross-check",
            )
            destination.write_text(text, encoding="utf-8")
    for relative, specification in PRUNED_MODULES.items():
        _write_pruned_module(relative, specification)

    stale = PACKAGE_ROOT / "verify_local_quantum_integer.py"
    if stale.is_file():
        stale.unlink()


def verify_python_closure() -> None:
    project_modules = {
        path.stem
        for path in PROJECT_ROOT.glob("*.py")
        if path.is_file()
    }
    packaged_modules = {
        Path(relative).stem
        for relative in ENGINE_FILES
        if relative.endswith(".py")
    }
    missing: set[str] = set()
    for relative in ENGINE_FILES:
        if not relative.endswith(".py"):
            continue
        tree = ast.parse(
            (PACKAGE_ROOT / relative).read_text(encoding="utf-8"),
            filename=relative,
        )
        # Only top-level imports are required by the reader-facing execution
        # path.  Some copied geometry modules retain lazy regression adapters
        # for the development tree; those are never called by ruling_polygon.
        for node in tree.body:
            if isinstance(node, ast.Import):
                roots = {alias.name.split(".", 1)[0] for alias in node.names}
            elif isinstance(node, ast.ImportFrom) and node.module:
                roots = {node.module.split(".", 1)[0]}
            else:
                continue
            missing.update(
                root
                for root in roots
                if root in project_modules and root not in packaged_modules
            )
    if missing:
        raise RuntimeError(
            "Python dependency closure is incomplete: " + ", ".join(sorted(missing))
        )


def verify_javascript_closure() -> None:
    packaged = {
        Path(relative).name
        for relative in ENGINE_FILES
        if relative.endswith(".js")
    }
    missing: set[str] = set()
    for relative in ENGINE_FILES:
        if not relative.endswith(".js"):
            continue
        text = (PACKAGE_ROOT / relative).read_text(encoding="utf-8")
        missing.update(
            dependency
            for dependency in JS_LOCAL_REQUIRE.findall(text)
            if dependency not in packaged
        )
    if missing:
        raise RuntimeError(
            "JavaScript dependency closure is incomplete: "
            + ", ".join(sorted(missing))
        )


def package_members() -> list[Path]:
    relative_paths = sorted(
        {Path(relative) for relative in (*SUPPORT_FILES, *ENGINE_FILES)},
        key=lambda path: path.as_posix(),
    )
    for relative in relative_paths:
        source = PACKAGE_ROOT / relative
        if not source.is_file():
            raise FileNotFoundError(f"missing package member: {source}")
        if source.name in EXCLUDED_BASENAMES:
            raise RuntimeError(f"excluded file entered package: {source.name}")
    return relative_paths


def write_archive(members: list[Path]) -> str:
    DOWNLOADS.mkdir(parents=True, exist_ok=True)
    temporary = ARCHIVE.with_name(ARCHIVE.name + ".tmp")
    with zipfile.ZipFile(
        temporary,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=9,
    ) as archive:
        for relative in members:
            data = (PACKAGE_ROOT / relative).read_bytes()
            archive_name = (Path(PACKAGE_NAME) / relative).as_posix()
            info = zipfile.ZipInfo(archive_name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            info.create_system = 3
            archive.writestr(info, data, compress_type=zipfile.ZIP_DEFLATED)
    temporary.replace(ARCHIVE)
    return hashlib.sha256(ARCHIVE.read_bytes()).hexdigest()


def main() -> int:
    refresh_sources()
    verify_python_closure()
    verify_javascript_closure()
    members = package_members()
    digest = write_archive(members)
    print(f"{ARCHIVE} ({ARCHIVE.stat().st_size} bytes)")
    print(f"sha256={digest}")
    print(f"members={len(members)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
