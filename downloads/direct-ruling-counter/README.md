# Direct Ruling Counter

This download counts rulings directly from a convex lattice polygon.  It
constructs a primitive multigeodesic diagram, searches its switch smoothings
with rollback pruning, and checks both disk and annular eyes exactly.

The computation does **not** call a tropical counter and does not look up
previously known answers.

## Requirements

- Python 3.10 or newer
- Node.js 18 or newer, available as `node` on `PATH`

The Python code uses only the standard library; there is no `pip install`
step.

## Quick start

The bundled example is the triangle with vertices

```text
(-1,-1), (1,0), (0,1).
```

It has the unique interior lattice point `(0,0)`.  From this directory, run:

```console
python ruling_polygon.py example_polygon.json
```

With no genus option, the program computes every genus from `0` through the
number of interior lattice points.  To compute only one genus:

```console
python ruling_polygon.py example_polygon.json --genus 0
```

Useful options include:

```console
python ruling_polygon.py example_polygon.json --jobs 4
python ruling_polygon.py example_polygon.json --certificates
python ruling_polygon.py example_polygon.json --output result.json
python smoke_test.py
```

## Input format

A single polygon may be written as:

```json
{
  "vertices": [[-1, -1], [1, 0], [0, 1]]
}
```

The vertices must be integer points listed cyclically around a strictly
convex polygon.  Clockwise input is accepted and normalized.  The first
vertex may be repeated at the end, but it is not required.

For a batch, use:

```json
{
  "polygons": [
    {"id": "triangle-a", "vertices": [[-1, -1], [1, 0], [0, 1]]},
    {"id": "triangle-b", "vertices": [[-2, -1], [1, 0], [0, 1]]}
  ]
}
```

The output separates `all_disk`, `annular`, and `total` counts for each
genus.  `--certificates` also retains the smoothing masks and annular
certificates; those records can be substantially larger.

## Performance

The outer ruling search is exponential in the number of crossings, although
rollback connectivity and topology pruning remove many branches.  Small
polygons usually finish quickly; polygons with long boundary or many
crossings can require much more time and memory.  Use `--genus` when only one
genus is needed and `--jobs` when counting independent polygons or genera on
a multicore machine.

## Reproducibility and scope

The polygon determines a stable identifier when no `id` is supplied.  The
diagram construction uses deterministic rational translates.  Counts are
reported as JSON, including the shear, translates, basepoint, crossing count,
and search diagnostics used for that run.

No license file is included in this preview package.  Please contact the
repository owner about reuse until a license is selected.
