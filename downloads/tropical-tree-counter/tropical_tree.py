#!/usr/bin/env python3
"""Count refined rational tropical trees from a convex lattice polygon.

The subset dynamic program uses labelled primitive ends with prescribed
moments. It does not use a ruling count or an assumed ruling correspondence.
Only the optional, subsequent comparison evaluates a supplied ruling polynomial.
Python 3.10 or later; standard library only.
"""
from __future__ import annotations

import argparse
from bisect import bisect_left
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from math import comb, gcd
from pathlib import Path
import random
import time
from typing import Iterable, Sequence

ALGORITHM_VERSION = "exact-tropical-tree-dp-v1"
DEFAULT_SEED = 17471
DEFAULT_MAX_ENDS = 16
Point = tuple[int, int]
Laurent = dict[int, int]


class NonGenericMomentsError(ValueError):
    """The prescribed moments lie on a detected degeneracy wall."""


def _integer(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _area2(vertices: Sequence[Point]) -> int:
    return sum(x * yy - xx * y for (x, y), (xx, yy)
               in zip(vertices, (*vertices[1:], vertices[0])))


def validate_vertices(raw: object) -> tuple[Point, ...]:
    """Validate a cyclic strictly convex integer polygon and orient it CCW."""
    if not isinstance(raw, (list, tuple)):
        raise ValueError("vertices must be a list of integer pairs")
    vertices = []
    for point in raw:
        if (not isinstance(point, (list, tuple)) or len(point) != 2
                or not all(_integer(value) for value in point)):
            raise ValueError("every vertex must be a pair of integers")
        vertices.append(tuple(point))
    if len(vertices) > 1 and vertices[0] == vertices[-1]:
        vertices.pop()
    if len(vertices) < 3 or len(set(vertices)) != len(vertices):
        raise ValueError("a polygon needs at least three distinct vertices")
    area2 = _area2(vertices)
    orientation = 1 if area2 > 0 else -1
    for index, (x, y) in enumerate(vertices):
        xx, yy = vertices[(index + 1) % len(vertices)]
        for other, (px, py) in enumerate(vertices):
            if other in (index, (index + 1) % len(vertices)):
                continue
            if orientation * ((xx - x) * (py - y) - (yy - y) * (px - x)) <= 0:
                raise ValueError("vertices must be cyclically ordered and strictly convex")
    if area2 < 0:
        vertices.reverse()
    return tuple(vertices)


def _invariants(vertices: Sequence[Point]) -> tuple[int, int, int]:
    area2 = _area2(vertices)
    boundary = sum(gcd(abs(xx - x), abs(yy - y))
                   for (x, y), (xx, yy) in zip(vertices, (*vertices[1:], vertices[0])))
    numerator = area2 - boundary + 2
    if area2 <= 0 or numerator < 0 or numerator % 2:
        raise ValueError("polygon failed the lattice Pick check")
    return area2, boundary, numerator // 2


def _stable_id(vertices: Sequence[Point]) -> str:
    # Match the direct ruling CLI's translation/cyclic/reversal-stable IDs.
    keys = []
    for sequence in (tuple(vertices), tuple(reversed(vertices))):
        for index in range(len(sequence)):
            rotated = sequence[index:] + sequence[:index]
            ox, oy = rotated[0]
            keys.append(";".join(f"{x - ox},{y - oy}" for x, y in rotated))
    return "polygon-" + hashlib.sha256(min(keys).encode("ascii")).hexdigest()[:12]


def _polygon_key(vertices: Sequence[Point]) -> tuple[Point, ...]:
    # Comparisons require the same actual vertices, not merely a translation.
    vertices = validate_vertices(vertices)
    return min(vertices[index:] + vertices[:index] for index in range(len(vertices)))


def parse_records(payload: object) -> list[dict[str, object]]:
    """Accept the same polygon/batch JSON forms as ruling_polygon.py."""
    if isinstance(payload, dict) and "polygons" in payload:
        raw_records = payload["polygons"]
    elif isinstance(payload, dict) and "vertices" in payload:
        raw_records = [payload]
    elif isinstance(payload, list):
        raw_records = [{"vertices": payload}]
    else:
        raise ValueError("input needs vertices, a polygons list, or a bare vertex list")
    if not isinstance(raw_records, list) or not raw_records:
        raise ValueError("polygons must be a nonempty list")
    records = []
    identifiers = set()
    for record in raw_records:
        if not isinstance(record, dict):
            raise ValueError("each polygon record must be an object")
        vertices = validate_vertices(record.get("vertices"))
        identifier = str(record.get("id") or _stable_id(vertices))
        if identifier in identifiers:
            raise ValueError(f"duplicate polygon id: {identifier}")
        identifiers.add(identifier)
        records.append({"id": identifier, "vertices": [list(point) for point in vertices]})
    return records


def _add(left: Laurent, right: Laurent) -> Laurent:
    answer = left.copy()
    for power, coefficient in right.items():
        answer[power] = answer.get(power, 0) + coefficient
    return {power: coefficient for power, coefficient in answer.items() if coefficient}


def _multiply(left: Laurent, right: Laurent) -> Laurent:
    answer: dict[int, int] = defaultdict(int)
    for lp, lc in left.items():
        for rp, rc in right.items():
            answer[lp + rp] += lc * rc
    return dict(answer)


def _quantum(multiplicity: int) -> Laurent:
    return {power: 1 for power in range(1 - multiplicity, multiplicity, 2)}


def _terms(polynomial: Laurent) -> list[dict[str, int]]:
    return [{"t_power": power, "coefficient": coefficient}
            for power, coefficient in sorted(polynomial.items()) if coefficient]


def _moment_text(value: Fraction) -> int | str:
    return value.numerator if value.denominator == 1 else str(value)


def _nongeneric(detail: str) -> None:
    raise NonGenericMomentsError(detail + "; choose another seed or explicit moments")


def count_tropical(
    vertices: Sequence[Point], *, seed: int = DEFAULT_SEED,
    moments: Sequence[int | Fraction] | None = None,
    max_ends: int = DEFAULT_MAX_ENDS,
) -> dict[str, object]:
    """Compute the exact refined count using prescribed labelled end moments.

    Explicit moments may specify b-1 non-root ends, or all b ends with sum zero.
    An end in direction u has moment det(u, x). The last end is the root.
    A detected non-generic configuration raises NonGenericMomentsError.
    """
    started = time.perf_counter()
    polygon = validate_vertices(vertices)
    area2, boundary, interior = _invariants(polygon)
    if not _integer(max_ends) or max_ends < 3:
        raise ValueError("max_ends must be an integer at least 3")
    if boundary > max_ends:
        raise ValueError(f"polygon has {boundary} ends, exceeding max_ends={max_ends}; "
                         "explicitly raise --max-ends to allow this exponential computation")
    if not _integer(seed):
        raise ValueError("seed must be an integer")
    ends = []
    for (x, y), (xx, yy) in zip(polygon, (*polygon[1:], polygon[0])):
        dx, dy = xx - x, yy - y
        weight = gcd(abs(dx), abs(dy))
        ends.extend([(dy // weight, -dx // weight)] * weight)
    n = boundary - 1
    if moments is None:
        mu = [Fraction(value) for value in
              random.Random(seed).sample(range(-10**12, 10**12), n)]
        mu.append(-sum(mu))
    else:
        if not isinstance(moments, (list, tuple)) or len(moments) not in (n, boundary):
            raise ValueError(f"explicit moments need {n} or {boundary} entries")
        if any(isinstance(value, bool) or not isinstance(value, (int, Fraction))
               for value in moments):
            raise ValueError("explicit moments must be integers or fractions.Fraction values")
        mu = [Fraction(value) for value in moments]
        if len(mu) == n:
            mu.append(-sum(mu))
        elif sum(mu) != 0:
            raise ValueError("the moments of all ends must sum to zero")
    seen = set()
    for direction, moment in zip(ends, mu):
        key = (direction, moment)
        if key in seen:
            _nongeneric("coincident parallel ends, including the root")
        seen.add(key)

    size = 1 << n
    directions = [(0, 0)] * size
    sums = [Fraction(0)] * size
    for subset in range(1, size):
        bit = subset & -subset
        index = bit.bit_length() - 1
        previous = subset ^ bit
        directions[subset] = (directions[previous][0] + ends[index][0],
                              directions[previous][1] + ends[index][1])
        sums[subset] = sums[previous] + mu[index]
        # Complements cover proper balanced subsets containing the root.
        if directions[subset] == (0, 0) and sums[subset] == 0:
            _nongeneric(f"proper balanced subset {subset:#x} has zero moment sum")

    states: list[tuple[list[Fraction], list[Laurent]] | None] = [None] * size
    stored_states = candidate_splits = nonparallel_splits = 0

    def query(subset: int, x: Fraction, y: Fraction) -> Laurent:
        if subset & (subset - 1) == 0:
            return {0: 1}
        state = states[subset]
        if state is None:
            return {}
        values, suffix = state
        position = directions[subset][0] * x + directions[subset][1] * y
        index = bisect_left(values, position)
        if index < len(values) and values[index] == position:
            _nongeneric(f"zero-length bounded edge at subset {subset:#x}")
        return suffix[index]

    for subset in sorted(range(1, size), key=int.bit_count):
        if subset & (subset - 1) == 0 or directions[subset] == (0, 0):
            continue
        nodes: dict[Fraction, Laurent] = {}
        first_bit = subset & -subset
        left = (subset - 1) & subset
        while left:
            # Exactly one order of every unordered binary split is used.
            if left & first_bit:
                right = subset ^ left
                ax, ay = directions[left]
                bx, by = directions[right]
                determinant = ax * by - ay * bx
                candidate_splits += 1
                if determinant == 0:
                    if ((ax, ay) != (0, 0) and (bx, by) != (0, 0)
                            and ax * sums[right] == bx * sums[left]
                            and ay * sums[right] == by * sums[left]):
                        _nongeneric(f"coincident parallel subset supports {left:#x}, {right:#x}")
                elif ((left & (left - 1) == 0 or states[left] is not None)
                      and (right & (right - 1) == 0 or states[right] is not None)):
                    nonparallel_splits += 1
                    x = (bx * sums[left] - ax * sums[right]) / determinant
                    y = (by * sums[left] - ay * sums[right]) / determinant
                    lp = query(left, x, y)
                    if lp:
                        rp = query(right, x, y)
                        if rp:
                            weight = _multiply(_multiply(lp, rp), _quantum(abs(determinant)))
                            position = directions[subset][0] * x + directions[subset][1] * y
                            nodes[position] = _add(nodes.get(position, {}), weight)
            left = (left - 1) & subset
        if nodes:
            values = sorted(nodes)
            suffix: list[Laurent] = [{} for _ in range(len(values) + 1)]
            for index in range(len(values) - 1, -1, -1):
                suffix[index] = _add(suffix[index + 1], nodes[values[index]])
            states[subset] = values, suffix
            stored_states += len(values)
    polynomial = states[-1][1][0] if states[-1] is not None else {}
    return {
        "id": _stable_id(polygon), "vertices": [list(point) for point in polygon],
        "double_area": area2, "boundary_lattice_points": boundary,
        "interior_lattice_points": interior,
        "seed": seed if moments is None else None,
        "moment_source": "seed" if moments is None else "explicit",
        "ends": [{"index": index, "direction": list(direction),
                  "moment": _moment_text(moment), "root": index == n}
                 for index, (direction, moment) in enumerate(zip(ends, mu))],
        "refined_terms": _terms(polynomial), "unrefined_count": sum(polynomial.values()),
        "elapsed_seconds": time.perf_counter() - started,
        "diagnostics": {"root_index": n, "nonroot_subsets": size - 1,
                        "stored_root_states": stored_states,
                        "candidate_splits": candidate_splits,
                        "nonparallel_splits": nonparallel_splits,
                        "detected_degeneracies": 0, "max_ends": max_ends},
    }


def count_records(records: Iterable[dict[str, object]], *, seed: int = DEFAULT_SEED,
                  max_ends: int = DEFAULT_MAX_ENDS) -> dict[str, object]:
    polygons = []
    for record in records:
        result = count_tropical(record["vertices"], seed=seed, max_ends=max_ends)
        result["id"] = str(record.get("id") or result["id"])
        polygons.append(result)
    return {"schema": "tropical-tree-counts-v1", "algorithm_version": ALGORITHM_VERSION,
            "variable": "t=q^(1/2)",
            "method": "exact labelled tropical tree subset dynamic program with prescribed end moments",
            "uses_ruling_counts": False, "polygons": polygons}


def _ruling_terms(genera: dict[int, int]) -> Laurent:
    # This algebraic transformation is used only AFTER the independent count.
    polynomial: Laurent = {}
    for genus, count in genera.items():
        for index in range(2 * genus + 1):
            power = 2 * genus - 2 * index
            coefficient = count * (-1) ** index * comb(2 * genus, index)
            polynomial[power] = polynomial.get(power, 0) + coefficient
    return {power: coefficient for power, coefficient in polynomial.items() if coefficient}


def compare_rulings(tropical_payload: dict[str, object], ruling_payload: object) -> dict[str, object]:
    """Compare already computed tropical terms with complete finite ruling data."""
    if not isinstance(ruling_payload, dict) or not isinstance(ruling_payload.get("polygons"), list):
        raise ValueError("ruling comparison input needs a polygons list")
    by_polygon = {}
    identifiers = set()
    for record in ruling_payload["polygons"]:
        if not isinstance(record, dict):
            raise ValueError("each ruling polygon must be an object")
        key = _polygon_key(record.get("vertices"))
        if key in by_polygon:
            raise ValueError("ruling comparison has duplicate polygon vertices")
        identifier = record.get("id")
        if identifier is not None:
            if not isinstance(identifier, str) or identifier in identifiers:
                raise ValueError("ruling polygon ids must be unique strings")
            identifiers.add(identifier)
        by_polygon[key] = record
    comparisons = []
    for polygon in tropical_payload["polygons"]:
        vertices = validate_vertices(polygon["vertices"])
        record = by_polygon.get(_polygon_key(vertices))
        if record is None:
            raise ValueError(f"no ruling polygon with the same vertices as {polygon['id']}")
        rows = record.get("genera")
        if not isinstance(rows, list):
            raise ValueError("ruling polygon needs a genera list")
        genera = {}
        for row in rows:
            if not isinstance(row, dict):
                raise ValueError("each ruling genus must be an object")
            genus, total = row.get("genus"), row.get("total")
            if not _integer(genus) or genus < 0 or genus in genera:
                raise ValueError("ruling genera must be distinct nonnegative integers")
            if row.get("phase_cardinality") != "FINITE" or not _integer(total) or total < 0:
                raise ValueError("each ruling genus needs FINITE phase_cardinality and a nonnegative integer total")
            genera[genus] = total
        interior = _invariants(vertices)[2]
        if set(genera) != set(range(interior + 1)):
            raise ValueError(f"ruling data for {polygon['id']} must contain every genus 0..{interior}, exactly once")
        expected = _ruling_terms(genera)
        actual = {term["t_power"]: term["coefficient"] for term in polygon["refined_terms"]}
        differences = [{"t_power": power, "tropical": actual.get(power, 0),
                        "ruling": expected.get(power, 0)}
                       for power in sorted(set(actual) | set(expected))
                       if actual.get(power, 0) != expected.get(power, 0)]
        comparisons.append({"id": polygon["id"], "ruling_id": record.get("id"),
                            "match": not differences, "tropical_terms": _terms(actual),
                            "ruling_terms": _terms(expected), "differences": differences})
    return {"formula": "N_Delta(q)=sum_g total[g]*(t-t^(-1))^(2g), t=q^(1/2)",
            "all_match": all(row["match"] for row in comparisons), "polygons": comparisons}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--max-ends", type=int, default=DEFAULT_MAX_ENDS)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare-rulings", type=Path)
    arguments = parser.parse_args(argv)
    try:
        records = parse_records(json.loads(arguments.source.read_text(encoding="utf-8")))
        payload = count_records(records, seed=arguments.seed, max_ends=arguments.max_ends)
        # Reading and comparing ruling totals occurs only after every tree count.
        if arguments.compare_rulings is not None:
            ruling_payload = json.loads(arguments.compare_rulings.read_text(encoding="utf-8"))
            payload["comparison"] = compare_rulings(payload, ruling_payload)
        text = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
        if arguments.output is None:
            print(text, end="")
        else:
            arguments.output.parent.mkdir(parents=True, exist_ok=True)
            arguments.output.write_text(text, encoding="utf-8")
    except (OSError, ValueError) as error:
        parser.error(str(error))
    return 0 if payload.get("comparison", {}).get("all_match", True) else 1


if __name__ == "__main__":
    raise SystemExit(main())
