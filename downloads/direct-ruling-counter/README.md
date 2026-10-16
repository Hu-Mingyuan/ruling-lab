# Direct Ruling Counter

This download counts connected rulings directly from a convex lattice polygon. It
constructs a primitive multigeodesic diagram, searches its switch smoothings
with rollback pruning, and checks both disk and annular eyes exactly.

The computation does **not** call a tropical counter and does not look up
previously known answers.

See the [definition and independence audit](https://github.com/Hu-Mingyuan/ruling-lab/blob/one-interior-preview/docs/ruling-definition-audit.md)
for the correspondence between the ruling conditions and the implementation.

Use this package when you want the numbers of rulings in each genus, or
individual smoothing certificates. For the independent refined tropical tree
calculation, use the separate
[Tropical Tree Counter](../tropical-tree-counter/README.md) package. Both
downloads and their source are available on the
[Ruling Atlas website](https://hu-mingyuan.github.io/ruling-lab/).

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

The example returns total counts **3 in genus 0 and 1 in genus 1**, so its
ruling polynomial is `R_Delta(z) = 3 + z^2`. The independent tropical answer
is `N_Delta(q) = q + 1 + q^(-1)`. Substituting
`z = q^(1/2) - q^(-1/2)` makes these polynomials equal. In particular,
`N_Delta(1)` corresponds to the genus-zero ruling count, not the sum of the
counts over all genera. Such finite comparisons test the refined ruling
conjecture; they are not an all-polygon proof.

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

In the output, `polygons[].genera` lists `genus`, `all_disk`, `annular`, and
`total`. Compare **total** when comparing different realizations: the
all-disk/annular split depends on the chosen diagram. `phase_cardinality`
records whether the exact annular phase calculation is finite; a non-finite
result is reported explicitly rather than replaced by a bounded search.

For an automated polynomial comparison, save the complete output without
`--genus` and follow the two-package example in the tropical counter's README.

## Files in this package

| File | Purpose |
| --- | --- |
| `ruling_polygon.py` | Main command-line interface: polygon input and counts by genus. |
| `example_polygon.json` | Small ready-to-run input. |
| `fig3_direct_ruling_counts.py`, `generic_ruling_perturbation.py` | Construct the counting diagram and separate coincident sweep events. |
| `direct_ruling_count.py`, `multigeodesic_annular_input.py` | Send the exact diagram to the search engine and return its results. |
| `ruling_smoothing.py`, `multigeodesic_alternating_audit.py`, `verify_*.py` | Supporting geometry and topology routines; these are not additional required commands. |
| `tmp/*.js` | Required disk/annular search and exact phase solver. Despite the folder name, keep these files. |
| `smoke_test.py`, `test_*.py`, `test_engine_regressions.js` | Installation checks and regression tests. |

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

Version `direct-ruling-search-v3` allows nonprimitive full periods of annular
boundaries. The complete deck quotient `Z^2 / Z*h` is used for phases,
same-lift checks, and translated-eye multiplicity. Version 2 incorrectly
discarded some such rulings: the stress polygon
`[(0,5),(1,2),(4,0),(3,3)]` had genus-one count 175 instead of 196. The repair
comes from explicit missing ruling certificates and the full-period deck
calculation, independently of tropical enumeration. The output's `count_scope`
field records the connected-surface convention.

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
node test_annular_periods.js
python smoke_test.py
```

No license file is included in this preview package.  Please contact the
repository owner about reuse until a license is selected.
