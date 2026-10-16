# Refined Tropical Tree Counter

This standalone program computes the refined tropical Laurent polynomial
`N_Delta(q)` of a convex lattice polygon from rational tropical trees with
fixed generic boundary moments. It uses exact arithmetic and a subset dynamic
program. It does not call the ruling counter or look up known answers.

The trees have abstract genus zero; their images in the plane may have
crossings. This is not a program for enumerating tropical curves in each
genus. Its refined polynomial can instead be compared with the
genus-indexed ruling counts through the conjectural identity below.

## Requirements and quick start

Use Python 3.10 or newer. Only the Python standard library is required;
there is no `pip install` step and no Node.js requirement for this package.
Extract the ZIP and, from inside `tropical-tree-counter`, run:

```console
python tropical_tree.py example_polygon.json --seed 17471 --output tropical.json
```

Use `python3` instead of `python` if that is your Python 3 command. Without
`--output`, the program prints JSON to standard output. Available options:

```console
python tropical_tree.py --help
python tropical_tree.py examples.json --seed 98189 --output batch.json
python tropical_tree.py example_polygon.json --compare-rulings ruling.json --output comparison.json
```

The last command requires a complete output file produced independently by
the direct ruling counter. See the two-package example below.

## Input

The input format is shared with the direct ruling counter:

```json
{
  "vertices": [[-1, -1], [1, 0], [0, 1]]
}
```

Vertices must be integer points listed cyclically around a strictly convex
polygon. Clockwise input is accepted and normalized. You may repeat the first
vertex at the end. A bare list of vertices is also accepted. For a batch:

```json
{
  "polygons": [
    {"id": "triangle", "vertices": [[-1, -1], [1, 0], [0, 1]]},
    {"id": "square", "vertices": [[0, 0], [1, 0], [1, 1], [0, 1]]}
  ]
}
```

Supply unique string IDs, or omit them to generate stable IDs. Use the same
polygon input for both counters when comparing their outputs.

## Output and variable convention

The top-level JSON contains `schema`, `algorithm_version`, `variable`,
`uses_ruling_counts`, and a `polygons` array. Each polygon records:

| Field | Meaning |
| --- | --- |
| `id`, `vertices` | Polygon identifier and normalized vertices |
| `double_area`, `boundary_lattice_points`, `interior_lattice_points` | Lattice-polygon invariants |
| `refined_terms` | Nonzero Laurent terms, each with integer `t_power` and `coefficient` |
| `unrefined_count` | Sum of the coefficients, equal to `N_Delta(1)` |
| `seed`, `moment_source` | The seed used to generate moments, or `null` and `explicit` for supplied moments |
| `ends` | Every labelled end's index, primitive direction, exact moment, and root flag |
| `elapsed_seconds`, `diagnostics` | Timing and enumeration diagnostics |

The stored variable is **`t = q^(1/2)`**. Thus a term with `t_power: 2`
means `q`, not `q^2`. The bundled triangle returns these terms (ordering
does not affect their meaning):

```json
[
  {"t_power": -2, "coefficient": 1},
  {"t_power": 0, "coefficient": 1},
  {"t_power": 2, "coefficient": 1}
]
```

This is `t^2 + 1 + t^(-2) = q + 1 + q^(-1)`, with
`unrefined_count: 3`. The output retains the actual end moments so the
enumeration can be reproduced, rather than relying only on a seed.

## What is counted

For each counterclockwise polygon edge `(dx, dy)`, let
`l = gcd(abs(dx), abs(dy))`. The boundary list contains `l` labelled copies
of the primitive outward normal `(dy/l, -dx/l)`. These directions sum to zero;
the number of ends is the boundary lattice-point count.

An end with direction `u_i` is constrained to a supporting line
`det(u_i, x) = m_i`, where the signed moments satisfy `sum_i m_i = 0`.
The CLI generates rationally exact moments from `--seed` and fixes one end
as the root. Parallel ends remain individually labelled by their fixed
boundary constraints. The count is **not divided by factorials** for repeated
directions.

At a trivalent vertex, with two outgoing integral direction vectors `u, v`,
the multiplicity is the quantum integer

```text
[abs(det(u,v))]_q = (t^abs(det(u,v)) - t^-abs(det(u,v))) / (t - t^-1),
t = q^(1/2).
```

The weight of a tree is the product of its vertex multiplicities, and
`N_Delta(q)` is the sum over the admissible trees. The program builds rooted
subtrees from subsets of the nonroot ends. It considers each unordered split
once, solves the moment-line intersections exactly, checks that bounded
edges have positive length, and accumulates the refined weights. It does
not enumerate rulings as an intermediate step.

This is a fixed-boundary enumerative problem. It should not be confused with
the usual tropical count through a collection of interior point constraints.

## Compare with direct ruling counts

Download the separate direct ruling package as well. Place the two extracted
directories beside each other. **From their common parent directory**, run:

```console
python direct-ruling-counter/ruling_polygon.py direct-ruling-counter/example_polygon.json --output ruling.json
python tropical-tree-counter/tropical_tree.py direct-ruling-counter/example_polygon.json --seed 17471 --compare-rulings ruling.json --output tropical.json
```

The first command additionally requires Node.js 18+ on `PATH`. Do not restrict
it with `--genus`: the comparison needs all genera from `0` through the number
of interior lattice points, and all total counts must be finite.

For `R_Delta(z) = sum_g r_(Delta,g) z^(2g)`, the conjectural comparison is

```text
N_Delta(q) = R_Delta(q^(1/2) - q^(-1/2))
          = sum_g r_(Delta,g) (q + q^(-1) - 2)^g.
```

The example has ruling totals `[3, 1]` in ascending genus, giving `3 + z^2`.
After substitution this is `q + 1 + q^(-1)`, exactly the tropical answer.
At `q = 1`, the identity gives the genus-zero ruling count, not the sum over
all genera. The disk/annular split of the ruling count depends on the chosen
diagram realization; the comparison uses each genus's total.

Tropical enumeration finishes before the comparison reads the ruling output.
Missing genera, nonfinite counts, and incompatible polygon data are rejected.
A comparison adds a report to the JSON; it does not replace the tropical
polynomial with the ruling answer. Agreement is a finite computational check
of the conjecture, not a proof for all polygons.

| Exit status | Meaning |
| --- | --- |
| `0` | Computation succeeded; a requested comparison also matched |
| `1` | Both calculations were available, but the comparison mismatched |
| `2` | Invalid input, unsupported size, or a degenerate moment choice |

## Genericity, performance, and reproducibility

The search is exponential in the number of boundary ends. The default limit
is 16 ends; an explicit option such as `--max-ends 18` raises it. This is a
guard against accidentally starting a very large computation, not a promise
that every 16-end polygon finishes quickly. Increasing the limit can require
substantial time and memory.

All line intersections and edge-length comparisons use exact rational
arithmetic. A degenerate moment choice, including an encountered zero-length
bounded edge, is reported as an error. The program does not silently discard
that degeneration and present the remaining sum as a complete count. Retry
with another seed, for example `--seed 98189`. Repeating a calculation with
different generic moments is a useful additional check.

The degeneracy checks are conservative: they can reject a local subtree
configuration that does not extend to a full tree. They are not a necessary
and sufficient description of the entire generic moment locus. A rejected
seed is a request to choose different moments, not a computed zero count.

For scripted use, import the independent core function from this directory:

```python
from tropical_tree import count_tropical

vertices = [[-1, -1], [1, 0], [0, 1]]
result = count_tropical(vertices, seed=17471, max_ends=16)
explicit = count_tropical(vertices, moments=[1, 3, -4])
print(result["refined_terms"])
```

`moments` accepts integers or `fractions.Fraction` values, in the end order
recorded by the program. Supply one moment per end with sum zero, or omit the
last (root) moment and let the program complete that sum. Explicit moments
set the output seed to `null`; they are still checked for degeneracy.

## Files and tests

- `tropical_tree.py`: standalone CLI, exact tropical enumeration, and optional
  comparison with ruling output.
- `example_polygon.json`: the small triangle used in this README.
- `examples.json`: a small batch of sample inputs.
- `test_tropical_tree.py`: computational and input-validation regression tests.
- `README.md`: this usage guide and mathematical conventions.

From inside this directory, run:

```console
python -m unittest -v test_tropical_tree.py
```

No license file is included in this preview package. Please contact the
repository owner about reuse until a license is selected.
