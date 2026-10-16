# Ruling Atlas

[Browse the atlas](https://hu-mingyuan.github.io/ruling-lab/) to see ruling
diagrams for convex lattice polygons. This repository also supplies two
independent command-line counters:

| Component | What it does | Requirements | Instructions |
| --- | --- | --- | --- |
| Static atlas | Displays saved polygons, ruling counts, and SVG diagrams | Web browser | [Website](https://hu-mingyuan.github.io/ruling-lab/) |
| Direct ruling counter | Constructs a torus diagram and counts connected rulings, including disk and annular eyes, in every requested genus | Python 3.10+ and Node.js 18+ | [README and source](downloads/direct-ruling-counter/README.md) |
| Refined tropical tree counter | Counts rational tropical trees with fixed boundary moments and returns the refined Laurent polynomial; optionally compares it with independently computed rulings | Python 3.10+ | [README and source](downloads/tropical-tree-counter/README.md) |

Both counters use only standard libraries. No Python or npm packages need to
be installed. The website displays precomputed data; it does not run either
counter in the browser. The ruling counter does not import the tropical
counter, and the tropical enumeration does not use ruling counts.

See the [definition audit](docs/ruling-definition-audit.md) for the actual
counting path and the correction of a nonprimitive-annulus omission in version 2.

## Download and run

Download and extract the package you need:

- [Direct ruling counter ZIP](https://hu-mingyuan.github.io/ruling-lab/downloads/direct-ruling-counter.zip)
- [Refined tropical tree counter ZIP](https://hu-mingyuan.github.io/ruling-lab/downloads/tropical-tree-counter.zip)

From inside the extracted `direct-ruling-counter` directory:

```console
python ruling_polygon.py example_polygon.json --output ruling.json
```

From inside the extracted `tropical-tree-counter` directory:

```console
python tropical_tree.py example_polygon.json --seed 17471 --output tropical.json
```

Use `python3` in place of `python` if that is the Python 3 command on your
system. Run either entry point with `--help` for its options.

### Input and output

Both programs accept the same polygon JSON, for example:

```json
{
  "vertices": [[-1, -1], [1, 0], [0, 1]]
}
```

List integer vertices cyclically around a strictly convex polygon. Clockwise
input is accepted. A batch has the form
`{"polygons": [{"id": "example", "vertices": [[-1,-1],[1,0],[0,1]]}]}`.
The package READMEs give the complete format and options.

The ruling output gives `all_disk`, `annular`, and `total` for each genus,
together with the exact diagram parameters and search diagnostics. The
example has total counts `[3, 1]` in ascending genus, hence
`R_Delta(z) = 3 + z^2`. The disk/annular split depends on the chosen diagram
realization; compare the total counts when checking the refined identity.
The counting diagram receives a certified small perturbation to separate
vertical crossing coordinates. The saved atlas pictures retain their
existing realizations.

The tropical output gives `refined_terms` as integer powers of
`t = q^(1/2)` with integer coefficients. The example gives
`t^2 + 1 + t^(-2) = q + 1 + q^(-1)`. It also records
`unrefined_count = N_Delta(1)`, the boundary directions and moments, the
actual seed, and search diagnostics. "Tree" refers to the abstract
genus-zero tropical source; its plane image can have crossings. It does not
mean an enumeration of tropical curves in every genus.

### Compare the two independent calculations

Extract both ZIPs into the same parent directory, so that
`direct-ruling-counter/` and `tropical-tree-counter/` are siblings.
Run these commands **from that parent directory**, using the same input
file in both commands:

```console
python direct-ruling-counter/ruling_polygon.py direct-ruling-counter/example_polygon.json --output ruling.json
python tropical-tree-counter/tropical_tree.py direct-ruling-counter/example_polygon.json --seed 17471 --compare-rulings ruling.json --output tropical.json
```

The comparison checks the conjectural identity

```text
N_Delta(q) = R_Delta(q^(1/2) - q^(-1/2))
          = sum_g r_(Delta,g) (q + q^(-1) - 2)^g.
```

Thus `3 + z^2` corresponds to `q + 1 + q^(-1)`. In particular,
`N_Delta(1) = r_(Delta,0)`, not the sum of the ruling counts over all genera.
The tropical count is completed independently before the comparison reads
the ruling results. The latter must contain every genus and finite total
counts. Exit status is `0` for a successful calculation (and a matching
comparison, when requested), `1` for a mismatch, or `2` for invalid input
or a degenerate moment choice. A finite comparison tests the conjecture;
it is not a proof for all polygons.

Both searches are exponential. The tropical program defaults to at most
16 boundary ends; increasing `--max-ends` explicitly can require substantial
time and memory. A degenerate moment choice is reported as an error; use a
different `--seed` rather than treating an incomplete enumeration as a count.
For the ruling program, `--genus` restricts the requested genus and
`--jobs` runs independent counting tasks concurrently.

## Repository guide

- [All downloads](downloads/README.md): choose a program and find its source.
- [Atlas contents and geometric conventions](docs/atlas-conventions.md):
  polygon collections, drawing representatives, vertical directions, and orbit sizes.
- [Computational validation](docs/validation.md): checks performed and their limits.
- [Definition and independence audit](docs/ruling-definition-audit.md):
  how the ruling conditions are implemented and the exact counting scope.
- [Higher-genus stress inputs](examples/stress-polygons.json): four irregular
  polygons with five or six interior lattice points.

Repository layout:

- `downloads/direct-ruling-counter/`: standalone ruling program, example,
  tests, and usage guide.
- `downloads/tropical-tree-counter/`: standalone tropical program, example,
  tests, and mathematical conventions.
- `downloads/*.zip`: matching standalone release archives.
- `data/` and `assets/`: precomputed atlas data and illustrations.
- `examples/`: additional inputs shared by both counters.
- Root HTML, JavaScript, and CSS: the static browser interface.
- `tools/`: atlas export and packaging utilities; some require the separate
  research repository and are not needed to use either downloaded counter.
- `tests/validate.mjs`: checks the static atlas data, figures, and download links.

Serve the repository root with any static HTTP server:

```console
python -m http.server 8000
```

Then open `http://localhost:8000/`. From a full checkout of this repository,
validate the site with:

```console
node tests/validate.mjs
```

The two package READMEs describe their own computational tests. Rebuilding
the atlas with `tools/export_cases.py` additionally requires the source
polygon catalogs and direct-count audits from the parent research repository;
it is not part of the standalone quick start.

## Contributing and release checks

Keep the two counters independent. A change to one counter should include
its relevant regression checks and updated usage instructions. Before a
release, run the package tests, regenerate its matching ZIP, and run the
static-site validator above; check that the ZIP contains the current source
and can run after extraction. Updating a counter does not require exporting
or replacing the saved atlas drawings.

Rebuild both ZIPs directly from the checked-in source with:

```console
python tools/package_counters.py
```

This packaging command is standalone and reproducible. The older
`tools/package_engine.py` refreshes files from the separate parent research
tree; it is not the release command for these checked-in counters.

The site is published separately from the author's main homepage. Every HTML
page carries `noindex` metadata to discourage search-engine indexing; this
is not access control. No license file is included in these preview packages;
please contact the repository owner about reuse.
