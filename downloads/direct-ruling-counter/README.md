# Direct Ruling Counter

This download counts rulings directly from a convex lattice polygon.  It
constructs a primitive multigeodesic diagram, searches its switch smoothings
with rollback pruning, and checks both disk and annular eyes exactly.

The computation does **not** call a tropical counter and does not look up
previously known answers.

Before counting, the engine checks whether the standard x-sweep has one
crossing on each critical slice.  If it does not, the engine deterministically
chooses a primitive integral covector in the same transverse chamber whose
crossing levels are all distinct.  Type-B local probes are then chosen below
one quarter of the minimum cyclic gap.  The effective sweep covector,
vertical direction, transverse covector, and probe size are included in the
JSON output.

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
number of interior lattice points.  Any single genus in that range may be
computed with `--genus`:

```console
python ruling_polygon.py example_polygon.json --genus 0
```

The second example has vertical edges in the input coordinates:

```console
python ruling_polygon.py example_vertical_polygon.json --genus 0
```

For that square the program automatically uses the source-frame vertical
direction `v=(-1,1)` rather than the invalid choice `v=(0,1)`.

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

The allowed range is

```text
0 <= genus <= number of interior lattice points.
```

The initial count frame is automatic.  If an edge has primitive direction
`(dx,dy)`, the program finds an integral shear `k` for which `dx+k*dy` is
nonzero on every edge; in the original coordinates its vertical direction is
`v=(-k,1)`.  If the resulting crossing slices are not generic, the exact
annular input layer then selects a primitive integral sweep covector in the
same chamber with distinct crossing levels.

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

Before constructing the ruling diagram, the program chooses an integral
vertical direction that is not parallel to any polygon edge.  Concretely, it
chooses a shear `C_k(x,y)=(x+k*y,y)`.  The vertical direction is then
`(-k,1)` in the input coordinates and `(0,1)` in the counting coordinates.
The dual sweep covectors are `(1,k)` and `(1,0)`, respectively.  These four
vectors are included in the JSON output, and the constructor rejects the
input arrangement if transversality fails.

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
reported as JSON, including the shear, translates, basepoint, effective
generic sweep, crossing count, and search diagnostics used for that run.

No license file is included in this preview package.  Please contact the
repository owner about reuse until a license is selected.
