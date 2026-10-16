# Downloadable programs

[Project guide](../README.md) · [Browse the atlas](https://hu-mingyuan.github.io/ruling-lab/)

Choose the program for the quantity you want to compute. Both accept convex
lattice polygons in the same JSON format, use only standard libraries, and
run locally after extracting their ZIP. The atlas itself displays saved
diagrams and does not run these computations in the browser.

| Program | Purpose | Requirements | Download and documentation |
| --- | --- | --- | --- |
| Direct ruling counter | Count disk and annular rulings in every requested genus | Python 3.10+ and Node.js 18+ | [ZIP](https://hu-mingyuan.github.io/ruling-lab/downloads/direct-ruling-counter.zip) · [README and source](direct-ruling-counter/README.md) |
| Refined tropical tree counter | Compute the refined Laurent polynomial from fixed boundary moments; optionally compare it with independently computed ruling totals | Python 3.10+ | [ZIP](https://hu-mingyuan.github.io/ruling-lab/downloads/tropical-tree-counter.zip) · [README and source](tropical-tree-counter/README.md) |

Each ZIP extracts to a directory with the corresponding name. The two source
directories in this folder are the readable versions of the same code.

## Quick start

From inside the extracted `direct-ruling-counter` directory:

```console
python ruling_polygon.py example_polygon.json --output ruling.json
```

From inside the extracted `tropical-tree-counter` directory:

```console
python tropical_tree.py example_polygon.json --seed 17471 --output tropical.json
```

To compare both outputs, extract both packages beside each other and use the
[complete two-package example](tropical-tree-counter/README.md#compare-with-direct-ruling-counts).
The tropical enumeration does not use ruling counts; the optional comparison
is performed afterward. The ruling counter does not call the tropical code.

The package READMEs explain input formats, output fields, variable conventions,
options, tests, and computational limits. Both searches can become expensive
on large polygons. For the scope of the existing checks, see the
[validation notes](../docs/validation.md).
