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
Without `--certificates`, the engine counts without constructing or retaining
the corresponding witnesses. The diagram is constructed once per polygon
and reused for all requested genera.

Before constructing the ruling diagram, the program chooses an integral
vertical direction that is not parallel to any polygon edge.  Concretely, it
chooses a shear `C_k(x,y)=(x+k*y,y)`.  The vertical direction is then
`(-k,1)` in the input coordinates and `(0,1)` in the counting coordinates.
The dual sweep covectors are `(1,k)` and `(1,0)`, respectively.  These four
vectors are included in the JSON output, and the constructor rejects the
input arrangement if transversality fails.

The counter also separates crossing events that have the same sweep
coordinate. It applies a deterministic rational horizontal perturbation and
certifies, using positive chamber triangulations, that this is an isotopy
fixing the chamber basepoints. All crossing x coordinates are then distinct;
the cyclic port order and future/past directions are preserved. The engine
uses the resulting piecewise-linear edge tangents. This perturbation concerns
the counting diagram only; it does not change the atlas illustrations.

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
reported as JSON, including the vertical direction, shear, translates,
basepoint, crossing count, and search diagnostics used for that run.
The output also records `algorithm_version`, `perturbation_step`, and
`perturbed_vertices`. Translate denominators are chosen coprime to every
horizontal strand direction, so parallel copies remain distinct.

Python geometry uses exact rational arithmetic. The JavaScript search uses
checked integer coordinates, with BigInt for exact turn and area signs and
for the annular phase solver. Inputs exceeding its supported exact integer
range are rejected explicitly instead of silently overflowing. The search
is still exponential; a large polygon may be impractical even when its
coordinates are supported.

Run the regression tests from this directory with:

```console
python -m unittest discover -p "test_*.py" -v
node test_engine_regressions.js
python smoke_test.py
```

No license file is included in this preview package.  Please contact the
repository owner about reuse until a license is selected.
