#!/usr/bin/env python3
"""Exact primitive-geodesic types and intersection geometry."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import floor, gcd

Point = tuple[Fraction, Fraction]
Vector = tuple[int, int]


def _det(left: tuple[int | Fraction, int | Fraction],
         right: tuple[int | Fraction, int | Fraction]) -> Fraction:
    return Fraction(left[0]) * Fraction(right[1]) - Fraction(left[1]) * Fraction(right[0])


def _mod_one(value: Fraction) -> Fraction:
    return value - floor(value)


def _mod_point(point: Point) -> Point:
    return _mod_one(point[0]), _mod_one(point[1])


@dataclass(frozen=True)
class Geodesic:
    name: str
    direction: Vector
    origin: Point


def _intersection_parameters(left: Geodesic, right: Geodesic) -> dict[Point, Fraction]:
    """Return every intersection point, indexed by the parameter on ``left``."""

    determinant = int(_det(left.direction, right.direction))
    if determinant == 0:
        raise ValueError("parallel geodesics are not supported")
    if gcd(abs(right.direction[0]), abs(right.direction[1])) != 1:
        raise ValueError("the right geodesic must be primitive")
    delta = (
        right.origin[0] - left.origin[0],
        right.origin[1] - left.origin[1],
    )
    offset = _det(delta, right.direction)
    answer: dict[Point, Fraction] = {}
    for residue in range(abs(determinant)):
        parameter = _mod_one((offset + residue) / determinant)
        point = _mod_point(
            (
                left.origin[0] + parameter * left.direction[0],
                left.origin[1] + parameter * left.direction[1],
            )
        )
        if point in answer:
            raise AssertionError("intersection formula produced a duplicate")
        answer[point] = parameter
    return answer


def _ray_compare(left: tuple[str, int, Vector], right: tuple[str, int, Vector]) -> int:
    """Counterclockwise comparison of two nonzero integer rays."""

    lv = left[2]
    rv = right[2]
    left_upper = lv[1] > 0 or (lv[1] == 0 and lv[0] > 0)
    right_upper = rv[1] > 0 or (rv[1] == 0 and rv[0] > 0)
    if left_upper != right_upper:
        return -1 if left_upper else 1
    cross = lv[0] * rv[1] - lv[1] * rv[0]
    if cross:
        return -1 if cross > 0 else 1
    return 0
