# Computational validation

The ruling and tropical programs enumerate different mathematical objects.
The tropical program does not use ruling results when enumerating trees;
comparison is a separate step after the tropical polynomial is computed.
These are finite computational checks of the refined ruling conjecture,
not a proof for all polygons.

## Atlas checks (2026-10-16)

The released `exact-tropical-tree-dp-v1` program was run independently with
moment seeds `17471` and `98189` on all 63 atlas polygons: the 16 one-interior
classes, the 45 two-interior classes, and the `O(4)` triangle and `O(2,4)`
rectangle. All **126 Laurent-polynomial comparisons** matched the stored
ruling polynomials after `z = t - t^(-1)`, with `t = q^(1/2)`.
No detected moment degeneracy occurred in these runs.

After the full-period annular repair in `direct-ruling-search-v3`, all genera
of the 61 census polygons were freshly recounted again, giving **167 matching
genus counts**, with unchanged runtime sources throughout the run. See the
[counts and source hashes](atlas-validation.json).

The `O(2,4)` rectangle was also freshly recounted with version 3 as
`[256, 96, 16, 1]` in ascending genus; see its [run record](o24-validation.json).
Thus the atlas regression comprises **62 polygons and 171 genus counts**.
The additional fresh `O(4)` ruling run was stopped before completion;
it is not included in the fresh ruling totals. Its tropical comparison above uses the
previously stored atlas polynomial `[304, 104, 16, 1]`.

The ruling update also passed 40 actual enumerations on four perturbed
diagrams, comparing incremental versus flat enumeration and count-only
versus certificate mode. These checks included a nonzero annular sector.
The package regression tests cover input convexity, exact generic
perturbation, integer overflow, and certificate suppression. Tests of the
tropical package cover direct polynomial answers, moment degeneracies,
input rejection, complete-genus comparison, and command-line exit status.

The version-3 tests additionally exercise full deck classes with positive
and negative sweep, distinguish nonprimitive torsion classes, and compare
50 translated-strip multiplicities with independent exact geometric counts.
An actual missing ruling is checked by a
[standalone geometry verifier](../examples/nonprimitive-witness/verify_nonprimitive_witness.py).
The [definition audit](ruling-definition-audit.md) explains the error and repair.

## Higher-genus stress tests

The four [stress inputs](../examples/stress-polygons.json) have five or six
interior lattice points. All **26 genus counts** were freshly enumerated with
`direct-ruling-search-v3`; both tropical moment seeds `17471` and `98189`
matched every Laurent coefficient in all four cases. The complete
[results and source hashes](stress-validation.json) are included.

| Case | Vertices | Interior points | Ruling coefficients, genus 0 upward |
| --- | --- | ---: | --- |
| Q5 | `(0,1), (1,0), (5,3), (5,4)` | 5 | `35, 105, 112, 54, 12, 1` |
| P5 | `(0,0), (1,0), (2,1), (3,4), (2,5)` | 5 | `84, 181, 156, 65, 13, 1` |
| S7 | `(0,5), (1,2), (4,0), (3,3)` | 6 | `49, 196, 294, 210, 77, 14, 1` |
| T13 | `(0,0), (1,0), (4,13)` | 6 | `13, 91, 182, 156, 65, 13, 1` |

The initial version-2 run matched Q5, P5, and T13, but returned 175 for S7
in genus one. The version-3 repair recovers 21 explicit nonprimitive-annulus
rulings there. The old 112 disk certificates and 63 annular certificates are
preserved, giving `112 + 84 = 196`. The other S7 genera are unchanged.
This discrepancy was retained and investigated; the ruling coefficients were
not replaced with values inferred from tropical enumeration.

Reproduce all their genus counts and both tropical moment choices with:

```console
python tools/compare_counters.py examples/stress-polygons.json --jobs 2 --output stress-results.json
```

This utility invokes the two public programs separately, then compares the
complete Laurent polynomials. It records the input, both outputs, seeds,
runtime-source hashes, and timings. It is not imported by either counter.
The test can take several minutes.

## Reproduce a comparison

From the repository root, with Python 3.10+ and Node.js 18+ available:

```console
python downloads/direct-ruling-counter/ruling_polygon.py downloads/direct-ruling-counter/example_polygon.json --output ruling.json
python downloads/tropical-tree-counter/tropical_tree.py downloads/direct-ruling-counter/example_polygon.json --seed 17471 --compare-rulings ruling.json --output tropical.json
```

Repeat the second command with `--seed 98189` for an independent generic
moment choice. The comparison checks every Laurent coefficient and requires
complete finite ruling totals for all genera. It does not compare the
realization-dependent disk/annular split.
