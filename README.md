# Rational Ruling Atlas — private preview

This repository is an unpublished static atlas of direct ruling counts and
diagrams for lattice polygons. It is separate from the public homepage and has
not been deployed.

The home page is a collection index. Its first link opens the collection
“Polygons with one interior lattice point.”

This preview contains the 16 lattice-equivalence classes of convex lattice
polygons with exactly one interior lattice point. Every entry is identified
only by its vertex coordinates, centered at the unique interior point `(0,0)`.

For each polygon the site records:

- its lattice diagram and vertex coordinates;
- area, boundary lattice-point count, and interior lattice-point count;
- its genus-zero rational ruling count;
- all-disk and annular sectors for the displayed deterministic realization;
- one SVG diagram for every rational ruling.

The preview displays 96 rational rulings in total. Its counts and diagrams are
generated from direct ruling certificates; the atlas exporter does not import
or consult a tropical counter.

## Run the private site locally

Serve the repository root with any static HTTP server, for example:

```text
python -m http.server 8000
```

Then open `http://localhost:8000/`.

## Downloadable counter

`downloads/direct-ruling-counter.zip` contains the standalone command-line
counter used to produce the atlas. It requires Python 3.10+ and Node.js 18+,
with no third-party Python or npm packages.

```text
python ruling_polygon.py example_polygon.json
```

The package is intentionally distributed without a license while this
repository remains private. Choose and add a license before public release.

## Rebuild and validate

From the parent `Nodal_Curves` working directory:

```text
python ruling-lab/tools/export_cases.py
node ruling-lab/tests/validate.mjs
```

The one-interior-point source data and direct-count audit live in the parent
research repository. This repository stores the static browser data and
generated SVG assets.
