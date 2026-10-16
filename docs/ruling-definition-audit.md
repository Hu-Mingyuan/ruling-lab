# What the direct ruling counter actually counts

The ruling counter enumerates smoothings of a torus diagram. It does not
invoke the tropical tree algorithm, read the atlas answer tables, or obtain
its coefficients by inverting the conjectural polynomial identity.
Computational independence and correctness of the mathematical implementation
are separate questions. This document records both the definition checks and
the scope of the implementation.

## Counting path and independence

The public input is a convex lattice polygon, with an optional requested
genus. The computation follows this path:

1. `ruling_polygon.py` and `fig3_direct_ruling_counts.py` construct rational
   geodesic translates, a transverse sweep, an alternating base chamber,
   and a certified small perturbation separating crossing coordinates.
2. `multigeodesic_annular_input.py` sends the resulting geometry and target
   Euler characteristic to `generic_incremental_ruling_search.js`.
3. The search makes a switch/no-switch choice at each crossing. Rollback
   connectivity tracks partial smoothing circuits and rejects impossible
   branches.
4. Completed circuits are tested as disk boundaries or matched in pairs as
   annular boundaries. Local coorientation, normality, balance, and surface
   topology are checked. Annular deck translations are solved using exact
   integer constraints, without a finite translation cutoff.
5. The engine adds one for each accepted disk configuration and the exact
   number of accepted annular configurations. The Python interface reports
   these two counts and their sum.

The runtime dependency closure consists of Python standard-library modules,
the Python files in this package, seven package JavaScript files, and the Node
built-ins `fs` and `node:crypto`. The `uses_tropical_count: false` output field
is descriptive metadata; it was not used as evidence for independence.

In a clean-extraction audit, the ruling ZIP was unpacked alone, without atlas
data or the tropical package. Python and Node read hooks restricted execution
to the interpreter, package source, input polygon, and standard input. Three
polygons, all their genera, and both certificate/count-only modes completed
successfully: 12 Node counting processes in total. Static imports, subprocess
arguments, and file-reading sites were also inspected.
The [machine-readable audit](independence-validation.json) records source and
archive hashes, imports, the guarded counting runs, and their outputs.

The tropical package's optional comparison reads ruling output only **after**
the tree enumeration has finished. Comparisons are not part of ruling
enumeration.

## Definition-to-code correspondence

Paths below are relative to
[`downloads/direct-ruling-counter/`](../downloads/direct-ruling-counter/).

| Mathematical condition | Implementation |
| --- | --- |
| A transverse vertical direction and at most one crossing on each vertical circle | `fig3_direct_ruling_counts.py`, `generic_ruling_perturbation.py`: rational shear and certified piecewise-linear perturbation |
| Equivariant path decomposition with allowed switch/non-switch resolutions | `ruling_smoothing.py`; `tmp/generic_incremental_ruling_search.js`: binary smoothing enumeration |
| Each path has constant upward/downward coorientation; paired paths have opposite coorientation | `traceCircuits` and `perfectMatchings` in `tmp/p2o4_annular_search.js` |
| Disk eyes have the prescribed two tips; annular eyes are monotone unbounded strips in the universal cover | `traceCircuits`, disk geometry checks, and `compileAnnularEyeResidue` in `tmp/annular_phase_exact.js` |
| Red coorientation points outward and blue coorientation inward | Circuit side/color data, `localAudit`, and the annular facing inequalities |
| Allowed local A/B configurations and type-B companion order | `localAudit`; exact `compileTypeBSameOrder` constraints on slices before and after each B switch |
| Blue multiplicity minus red multiplicity equals the chosen balance target at the basepoint | Exact translated-disk multiplicity and `compilePhaseBalance`; the polygon CLI uses target zero in its alternating chamber |
| Surface obtained by gluing the eyes at switches; orientability, boundary count, and genus | Parity connectivity in `localAudit` and the independent boundary-cycle calculation in `ribbonTopology` |

The internal H/V crossing names are not the paper's horizontal/vertical A
names: code `V` denotes a tip-producing switch. There are equally many disk
eyes and such switches, so `chi = disks - switches = -H`. For a connected
surface with `b` boundary components this gives `H = 2*g + b - 2`.

## Nonprimitive annular periods

An annular eye is a strip in the universal cover. Its torus projection need
not be embedded. In particular, its boundary's full deck period `h` need
not be primitive. Dividing `h` by its coordinate gcd can identify distinct
lifted boundaries and discard valid rulings.

For a primitive sweep covector `tau`, choose integral vectors `e,f` forming
a lattice basis with `tau(e)=1` and `tau(f)=0`. Orient the full boundary period
so that

```text
h = q*e + r*f,     q = tau(h) > 0.
```

Every class of `Z^2 / Z*h` has a unique representative

```text
a*e + z*f,         0 <= a < q,     z an integer.
```

Indeed, subtracting `floor(tau(D)/q)*h` from a deck vector `D` gives this
form; uniqueness follows by applying `tau` to the difference of two such
representatives. This parameterization retains the torsion when `h` is
nonprimitive. Matching boundaries must have the same **full** period,
not merely parallel primitive directions. Two translated boundaries are
the same lift precisely when their deck difference lies in `Z*h`.

Exact annular certificates record the phase basis and actual deck shifts.
Their stage is `connected-local-exact-phase`. The older `stripOk` field refers
to a stricter adjacent-strip check on the torus and is not the acceptance
criterion for these universal-cover eyes.

The version-2 implementation incorrectly required primitive annular periods.
The stress polygon `[(0,5),(1,2),(4,0),(3,3)]` exposed the problem: genus one
was reported as 175 instead of 196. Removing the restriction alone is not a
valid repair; deck phases, same-lift checks, and multiplicities must all use
the full-period quotient above.

An independent rational-geometry reconstruction of a missing smoothing
(`0x4303080` in the saved stress-test diagram) finds two disk eyes and one red
annular eye with full period `(2,-2)`. Its boundary graphs are strictly
separated, both B switches have the required companion order on both sides,
all A switches have the required colors and convexity, and the basepoint
balance is zero. The eye adjacency graph is connected, with two disks,
six switches, and four boundary components: `chi = 2 - 6 = -4`, hence genus
one. The half-period `(1,-1)` is explicitly not a symmetry of either boundary.
This is a geometric witness, not a
coefficient inserted from the tropical answer.

The [standalone witness checker](../examples/nonprimitive-witness/verify_nonprimitive_witness.py)
reads only [four explicit covering-space paths and the original geodesics](../examples/nonprimitive-witness/nonprimitive-witness-geometry.json).
It does not import either counter. It checks the complete edge decomposition,
disk and strip geometry, full period, all switch and nonswitch local conditions,
basepoint balance, and connected surface topology using rational arithmetic.
From the repository root:

```console
python examples/nonprimitive-witness/verify_nonprimitive_witness.py --output witness-check.json
```

The focused `test_nonprimitive_annulus.py` regression also passes this
nonprimitive smoothing and a primitive smoothing through the public flat
engine with three sweep covectors, including a reversed sweep. The separate
`test_annular_periods.js` compares translated-strip multiplicity with an
independent exact geometric enumeration.

## Scope and verification limits

Following the author's convention, the counter counts connected ruling
surfaces. `localAudit` checks connectedness and orientability before using
`chi = 2 - 2*g - b` to check the genus.

Finite regressions and tropical comparisons do not prove exhaustiveness for
all polygons. The program reports unsupported arithmetic and unresolved or
infinite phase systems explicitly. The searches remain exponential.
