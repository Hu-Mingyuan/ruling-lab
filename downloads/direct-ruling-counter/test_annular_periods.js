"use strict";

// Independent lattice and multiplicity checks for nonprimitive annular
// periods. The oracle below counts actual translated planar strips; it does
// not use the phase solver's floor/ceiling multiplicity formula.
const assert = require("node:assert/strict");
const phase = require("./tmp/annular_phase_exact.js");
const core = require("./tmp/p2o4_annular_search.js");

const periods = [[2, -2], [3, 0], [0, -2], [2, 4], [-2, -4], [1, 2]];
const sweeps = [[1, 0], [0, 1], [1, 1], [2, -1], [-1, -1]];
let latticeChecks = 0;
for (const [hx, hy] of periods) {
  for (const [tx, ty] of sweeps) {
    if (tx * hx + ty * hy === 0) continue;
    const basis = phase.phaseResidueData(hx, hy, tx, ty);
    for (let label = -8n; label <= 8n; label += 1n) {
      const shift = phase.phaseDeckShift(basis, label);
      assert.equal(phase.deckPhaseClass(basis, ...shift), label);
      for (const multiple of [-2n, 0n, 3n]) {
        assert.equal(
          phase.deckPhaseClass(basis,
            shift[0] + multiple * BigInt(hx),
            shift[1] + multiple * BigInt(hy)),
          label,
          "adding a full circuit period preserves the lifted-path class"
        );
      }
      latticeChecks += 1;
    }
    for (let dx = -2; dx <= 2; dx += 1) {
      for (let dy = -2; dy <= 2; dy += 1) {
        const encoded = phase.deckPhaseClass(basis, dx, dy);
        const representative = phase.phaseDeckShift(basis, encoded);
        const rx = BigInt(dx) - representative[0];
        const ry = BigInt(dy) - representative[1];
        assert.equal(BigInt(hx) * ry - BigInt(hy) * rx, 0n);
        assert.equal(hx !== 0 ? rx % BigInt(hx) : ry % BigInt(hy), 0n);
        assert.equal(core.deckInCircuitPeriod(hx, hy, dx, dy), encoded === 0n);
        latticeChecks += 1;
      }
    }
  }
}
assert.equal(core.deckInCircuitPeriod(2, -2, 1, -1), false,
  "a primitive parallel translation is not the actual period");

// L=10 is the torus lattice scale. This monotone polygonal line has least
// deck period (2,0), rather than (1,0): its heights at x=5 and x=15 differ.
// Translating by (1,0) gives a different lifted boundary, although the two
// paths can cross. Such intersections between different lifts are allowed.
const lower = [[0, 0], [5, 1], [10, 0], [15, 2], [20, 0]];
const upper = lower.map(([x, y]) => [x, y + 2]);
const segments = (points) => points.slice(0, -1).map(
  (point, index) => [...point, ...points[index + 1]]
);
const circuits = [
  { id: 0, kind: "essential", hx: 2, hy: 0, rawHx: 2, rawHy: 0,
    traceSide: 1, segments: segments(lower), visits: [] },
  { id: 1, kind: "essential", hx: 2, hy: 0, rawHx: 2, rawHy: 0,
    traceSide: -1, segments: segments(upper), visits: [] },
];

// All queried x coordinates are integers. Every segment has horizontal
// length 5, so this returns exactly five times its rational height.
function graphHeightTimesFive(points, x) {
  x = ((x % 20n) + 20n) % 20n;
  for (let index = 0; index + 1 < points.length; index += 1) {
    const [ax, ay] = points[index].map(BigInt);
    const [bx, by] = points[index + 1].map(BigInt);
    if (ax <= x && x < bx) return 5n * ay + (x - ax) * (by - ay);
  }
  throw new Error("point outside periodic graph");
}

function bruteMultiplicity([shiftX, shiftY]) {
  let count = 0n;
  // The eye's translates are indexed by (a,b), a=0,1 and b integral,
  // modulo its true period (2,0). Heights before translation are between
  // 0 and 4, so |b| > |shiftY|+2 cannot enclose the base point (3,3).
  const limit = (shiftY < 0n ? -shiftY : shiftY) + 2n;
  for (let a = 0n; a < 2n; a += 1n) {
    for (let b = -limit; b <= limit; b += 1n) {
      const low = graphHeightTimesFive(lower, 3n - a * 10n) + b * 50n;
      const high = graphHeightTimesFive(upper, 3n - (a + shiftX) * 10n)
        + (b + shiftY) * 50n;
      if (low < 15n && 15n < high) count += 1n;
    }
  }
  return count;
}

function boundariesFaceWithoutMeeting([shiftX, shiftY]) {
  // Translation by shiftX*10 preserves the breakpoint grid. The difference
  // is linear between these points, so testing the breakpoints is exact.
  return [0n, 5n, 10n, 15n].every((x) => (
    graphHeightTimesFive(upper, x - shiftX * 10n) + shiftY * 50n
      > graphHeightTimesFive(lower, x)
  ));
}

let multiplicityChecks = 0;
for (const tau of [[1, 0], [-1, 0], [1, 1], [-1, -1], [2, 1]]) {
  const sigma = [-tau[1], tau[0]];
  const basis = phase.phaseResidueData(2, 0, ...tau);
  for (let residue = 0n; residue < basis.q; residue += 1n) {
    const descriptor = phase.compileAnnularEyeResidue({
      arr: { L: 10, vertices: [] }, circuits, pair: [0, 1], color: 1,
      residue, tau, sigma, basePoint: [3, 3],
    });
    assert.equal(descriptor.ok, true);
    for (let z = -3n; z <= 3n; z += 1n) {
      const shift = phase.phaseDeckShift(basis, residue + basis.q * z);
      if (!boundariesFaceWithoutMeeting(shift)) continue;
      const expected = bruteMultiplicity(shift);
      const computed = descriptor.multiplicity.coefficient * z
        + descriptor.multiplicity.constant;
      assert.equal(computed, expected,
        `translate multiplicity for tau=${tau}, residue=${residue}, z=${z}`);
      multiplicityChecks += 1;
    }
  }
}
assert.equal(bruteMultiplicity([0n, 0n]), 1n);
assert.equal(bruteMultiplicity([0n, 1n]), 3n,
  "one transverse deck step adds one lift for each of the two sheets");
assert.ok(multiplicityChecks >= 40);
console.log(`annular period checks passed: ${latticeChecks} lattice checks, ${multiplicityChecks} exact translate counts`);
