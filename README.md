# Ruling Atlas

This repository is a static atlas of direct ruling counts and diagrams for
lattice polygons. It is published separately from the author's main homepage.
Every HTML page carries `noindex` metadata to discourage search-engine indexing;
this is not access control.
The repository also contains a standalone general Newton-polygon counter.

The home page links to three collections:

- the 16 lattice-equivalence classes with one interior lattice point;
- the 45 lattice-equivalence classes with two interior lattice points;
- polygons with more interior points, currently containing `O(4) on P^2`
  and `O(2,4) on P^1 x P^1`.

Every entry is identified by its vertex coordinates. The two interior-point
collections list every ruling separately; the more-interior-points collection
instead preserves symmetry representatives and labels every orbit size. Its
landing view first lists the two polygons, and selecting either card reveals
only that polygon's details.

For each polygon the site records:

- its lattice diagram and vertex coordinates;
- area, boundary lattice-point count, and interior lattice-point count;
- its ruling counts in every possible genus;
- all-disk and annular sectors for the displayed deterministic realization;
- one SVG diagram for every ruling in the interior-point collections;
- the single `SL(2,Z)`-equivalent polygon used to draw every ruling of the
  selected lattice-equivalence class.

Selecting a polygon first shows its common drawing representative. The ruling
galleries then appear from highest genus to lowest genus as independent
expandable sections.

For a fixed polygon, all-disk and annular diagrams share one fundamental square
and one fixed Lambda. Annular eyes are normalized only in lifted coordinates
and mapped back before rendering. Switches are the black points in the ruling
diagrams. Each entry labels its chosen `SL(2,Z)` representative, the matrix
from the catalog coordinates, and the vertical direction. The selected
representative also draws that direction as an arrow beside the polygon.
The counter first chooses a source-frame vertical direction `(-k,1)` whose
determinant with every primitive polygon edge is nonzero; after the counting
shear it is `(0,1)`.  The arrow is the image of this same vector in the
display frame, not a separately chosen drawing direction.

The plane cubic entry is explicitly labeled `O(3) on P^2`. It uses
`Conv{(-1,-1), (2,-1), (-1,2)}` with vertical direction `(-3,2)`. Its exact
symmetric three-family realization has a horizontal, diagonal, and vertical
family; its genus-one ruling consists of 18 equal triangular eyes and all 27
crossings are switches. In this displayed realization the genus-zero split is
9 all-disk and 0 annular.

The `O(2,2) on P^1 x P^1` entry uses the standard centered square
`Conv{(-1,-1), (1,-1), (1,1), (-1,1)}`. Its exact four-by-four Lambda grid
has four evenly spaced horizontal and four evenly spaced vertical components.
Its displayed vertical direction is `(-1,2)`, transverse to both families.
The genus-one ruling has eight equal square eyes and all 16 crossings are
switches.

The `O(4) on P^2` entry uses
`Conv{(-2,-1), (2,-1), (-2,3)}` for every genus, with displayed vertical
direction `(-1,2)`. This is the direction used by the ruling diagram itself:
the covector `(2,1)` is positive on the displayed future tangents
`(1,0)`, `(0,1)`, and `(1,-1)`. Its page contains one genus-three diagram,
one genus-two symmetry representative, six genus-one symmetry
representatives, and thirteen genus-zero symmetry representatives. The single
genus-two orbit has size 16; an exhaustive exact annular search without
translation reduction or a node cutoff gives zero additional genus-two
rulings. The six genus-one orbit sizes sum to 104, with zero annular rulings.
The genus-zero orbit sizes sum to 304 (288 all-disk and 16 annular). The
figures use green Lambda, red and blue eyes, black switches, and darker fills
where eyes overlap.

The `O(2,4) on P^1 x P^1` entry uses the rectangle
`Conv{(0,0), (4,0), (4,2), (0,2)}` in every genus. Its source-to-count matrix
is `[[1,-1],[1,0]]`, and its source-frame vertical direction is `(1,1)`.
Every figure is drawn on the same exact 8-by-4 rectangular Lambda. A complete
direct annular-phase audit gives the counts 1, 16, 96, and 256 in genera 3,
2, 1, and 0, with zero annular rulings in this equally spaced realization.
Quotienting by the 32 line-preserving diagram symmetries leaves 1, 1, 7, and
11 displayed representatives, respectively.

The centrally symmetric six-vertex polygon
`Conv{(-1,0), (0,-1), (1,-1), (1,0), (0,1), (-1,1)}` uses exact separated
phases. Its three parallel pairs have gaps `1/3`, `1/2`, and `1/2`, so no two
components are visually coincident. Recounting this realization directly gives
five all-disk and one annular genus-zero ruling.

The five-vertex class
`Conv{(-1,-1), (1,-1), (1,0), (0,1), (-1,1)}` uses a reflection-symmetric
exact realization: the horizontal and vertical phases are both
`{0,1/3,2/3}`, and the diagonal phase is `1/6`.

The two-interior collection was recounted after choosing exact
symmetry-oriented phases for each of its 45 classes. In particular, the five
triangle examples use congruent red and blue triangular eyes in their
highest-genus rulings. The `O(2,3) on P^1 x P^1` entry keeps the catalog
rectangle `Conv{(0,0), (3,0), (3,2), (0,2)}` and draws every ruling on one
fixed 6-by-4 rectangular Lambda; its highest-genus ruling has six red and six
blue equal square eyes.

The two census collections display all 1,601 rulings separately: 112 for the
one-interior collection and 1,489 for the two-interior collection. The
other-polygons page adds 41 saved symmetry representatives for `O(4)` and
`O(2,4)`. Every page records the ruling polynomial
`R_Delta(z) = sum_g r_(Delta,g) z^(2g)`. The census counts and diagrams are
generated from direct ruling certificates; the atlas exporter does not import
or consult a tropical counter.

## Run the site locally

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
python ruling_polygon.py example_polygon.json --genus 1
```

The input may be any strictly convex lattice Newton polygon with cyclically
listed integer vertices. The program automatically chooses an integral
vertical direction transverse to every edge, applies a generic sweep
perturbation when crossing levels coincide, and accepts every genus
`0 <= g <= I(Delta)`. Omitting `--genus` computes the complete range.

For an edge with primitive direction `(dx,dy)`, the program chooses an
integral shear `k` such that `dx+k*dy` is nonzero for every edge. In the source
coordinates this gives the transverse vertical direction `v=(-k,1)`. If two
crossings still occur on one sweep level, the annular solver automatically
chooses a primitive generic covector in the same chamber. The JSON result
records the effective `v`, sweep covector, transverse covector, and local
probe size.

No license file is currently included. Contact the repository owner about
reuse until a license is selected.

## Rebuild and validate

From the repository root:

```text
python tools/export_cases.py
node tests/validate.mjs
```

The source polygon catalogs and direct-count audits live in the parent research
repository. This repository stores the static browser data and generated SVG
assets.
