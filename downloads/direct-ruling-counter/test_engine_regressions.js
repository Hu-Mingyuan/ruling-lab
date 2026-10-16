#!/usr/bin/env node
"use strict";

const assert = require("node:assert/strict");
const core = require("./tmp/p2o4_annular_search.js");
const exact = require("./tmp/p2_annular_exact_search.js");
const phase = require("./tmp/annular_phase_exact.js");
const flat = require("./tmp/generic_annular_exact_search.js");
const incremental = require("./tmp/generic_incremental_ruling_search.js");
const { exactTurnSign } = require("./tmp/exact_number_guards.js");

// Exact public-builder fixture: label audit-triangle, polygon
// [(0,0),(1,k),(2,2k+1)], residues [397,840,909] / 1009.
function triangle(k) {
  const L = 1009;
  const coordinates = [[566, 397], [994, 397], [138, 978]];
  const portTwin = [6, 11, 4, 9, 2, 10, 0, 8, 7, 3, 5, 1];
  const portVertex = Array.from({ length: 12 }, (_, p) => Math.floor(p / 4));
  const portVectorX = [1, 1, -1, -1, 1, 2, -1, -2, 2, 1, -2, -1];
  const portVectorY = [0, 1, 0, -1, 0, 1, 0, -1, 1, 1, -1, -1];
  const portShiftX = [0, 1, -1, 0, 1, 2, 0, 0, 0, 0, -2, -1];
  const portShiftY = [0, 0, 0, -1, 0, 0, 0, -1, 1, 1, 0, 0];
  const floors = coordinates.map(([x, y]) => Math.floor((y + k * x) / L));
  const vertices = coordinates.map(([x, y], index) => {
    const p = 4 * index;
    return { index, x, y: (y + k * x) % L, ports: [p, p + 1, p + 2, p + 3],
      futurePorts: [p, p + 1], nonswitchPairs: [[p, p + 2], [p + 1, p + 3]],
      switchPairs: index === 0 ? [[p, p + 3], [p + 1, p + 2]]
        : [[p, p + 1], [p + 2, p + 3]], orientation: index === 0 ? "H" : "V" };
  });
  const shiftedY = portShiftY.map((y, p) => y + k * portShiftX[p]
    + floors[portVertex[portTwin[p]]] - floors[portVertex[p]]);
  return {
    arrangement: { L, vertices, portCount: 12, portVertex, portTwin, portVectorX,
      portVectorY: portVectorY.map((y, p) => y + k * portVectorX[p]),
      portShiftX, portShiftY: shiftedY,
      portSign: [1, 1, -1, -1, 1, -1, -1, 1, -1, 1, 1, -1],
      edges: portTwin.flatMap((b, a) => a < b
        ? [{ a, b, sx: portShiftX[a], sy: shiftedY[a] }] : []),
      crossingOrder: [2, 0, 1], targetChi: -1, targetGenus: 0, boundaryComponents: 3 },
    basePoint: ["123/2", `${(579 + k * 246) % (4 * L)}/4`],
    tau: [1, 0], sigma: [0, 1], collectAnswers: true,
  };
}

for (const k of [0, 40000, 100000]) {
  const payload = triangle(k);
  const expected = flat.run(payload);
  for (const order of ["crossing", "reverse", "hfirst", "vfirst"]) {
    const got = incremental.run(payload, { order });
    assert.equal(got.exactTotalRulingCount, 1);
    assert.equal(got.exactTotalRulingCount, expected.exactTotalRulingCount);
    assert.equal(got.exactAnnularRulingCount, expected.exactAnnularRulingCount);
  }
  for (const engine of [flat, incremental]) {
    const countOnly = engine.run({ ...payload, collectAnswers: false });
    assert.equal(countOnly.exactTotalRulingCount, 1);
    assert.deepEqual(countOnly.allDiskMasks || countOnly.allDiskResults, []);
    assert.deepEqual(countOnly.annularResults || countOnly.results, []);
  }
}

// A valid explicit O(12) standard ruling disproves the former root prune.
// Check the exact same suffix bounds without enumerating 2^432 masks.
const large = core.buildArrangement(12, 55);
const order = large.vertices.map((v) => v.index);
const bounds = incremental.remainingSwitchCounts(large, order);
assert.equal(bounds.remainingH[0], 144);
assert.equal(bounds.remainingV[0], 288);
assert.ok(bounds.remainingH[0] >= -large.targetChi);
const choice = new Uint8Array(large.vertices.length).fill(1);
const trace = core.traceCircuits(large, choice);
assert.equal(trace.circuits.length, 288);
const options = { tau: [2, 1], sigma: [1, 0], basePoint: ["37/100", "73/100"] };
const witness = exact.auditAllDiskLocalExact(core, large, choice, trace, options);
assert.equal(witness.ok, true);
assert.equal(witness.ribbonGenus, 55);
assert.equal(witness.ribbonBoundaryComponents, 36);
assert.equal(trace.circuits.reduce((sum, c) => sum + BigInt(c.color)
  * phase.translatedDiskMultiplicityExact(large, c, options.basePoint), 0n), 0n);

const dsu = new incremental.RollbackVoltageDSU(40001);
assert.equal(dsu.unite(0, 40000, 0, 40000), true);
assert.deepEqual(dsu.relation(0, 40000), [0, 40000]);
assert.equal(exactTurnSign(1, 0, 1, 129), 1);
assert.equal(exactTurnSign(100000000, 99999999, 100000001, 100000000), 1);
const turns = core.buildArrangement(1, 0);
turns.portVectorY = Float64Array.from(turns.portVectorY, (y) => 129 * y);
const turnTrace = core.traceCircuits(turns, new Uint8Array(3).fill(1));
for (const turn of turnTrace.turnAt) assert.equal(Math.abs(turn), 1);

const unsafe = triangle(0);
unsafe.arrangement.L = Number.MAX_SAFE_INTEGER + 1;
assert.throws(() => incremental.run(unsafe), /safe integer/);
const enormous = triangle(0);
enormous.arrangement.L = Number.MAX_SAFE_INTEGER;
assert.throws(() => flat.run(enormous), /exact integer geometry range/);
console.log("Engine regressions passed: shears, orders, count-only, O(12) witness, integer ranges.");
