#!/usr/bin/env node
"use strict";

/*
 * Exact no-cutoff annular audit for a small, explicitly supplied torus
 * arrangement.  Unlike the P^2 outer enumerator, this driver makes no use of
 * translation symmetry or a degree-specific rollback search: it simply
 * enumerates the requested smoothing masks and delegates every unbounded
 * relative deck phase to annular_phase_exact.js.
 */

const fs = require("fs");
const {
  auditAllDiskLocalExact,
  makeExactAuditor,
} = require("./p2_annular_exact_search.js");
const { translatedDiskMultiplicityExact } = require("./annular_phase_exact.js");

function hydrateArrangement(raw) {
  const arr = { ...raw };
  const signed16 = ["portVertex", "portTwin"];
  const signed32 = ["portVectorX", "portVectorY", "portShiftX", "portShiftY"];
  for (const name of signed16) arr[name] = Int16Array.from(raw[name]);
  for (const name of signed32) arr[name] = Int32Array.from(raw[name]);
  arr.portSign = Int8Array.from(raw.portSign);
  arr.vertices = raw.vertices.map((vertex) => ({
    ...vertex,
    ports: vertex.ports.slice(),
    futurePorts: new Set(vertex.futurePorts),
    nonswitchPairs: vertex.nonswitchPairs.map((pair) => pair.slice()),
    switchPairs: vertex.switchPairs.map((pair) => pair.slice()),
  }));
  arr.edges = raw.edges.map((edge) => ({ ...edge }));
  arr.crossingOrder = (raw.crossingOrder || raw.vertices.map((vertex) => vertex.index)).slice();
  return arr;
}

function diagnosticsObject() {
  return {
    infinite: [],
    unknown: [],
    solver: {
      calls: 0,
      finite: 0,
      empty: 0,
      infinite: 0,
      unknown: 0,
      activeSets: 0,
      integerNodes: 0,
      maxActiveSets: 0,
      maxIntegerNodes: 0,
    },
  };
}

function parseMask(value) {
  if (typeof value === "number") return BigInt(value);
  if (typeof value === "string") return BigInt(value);
  throw new Error(`bad smoothing mask ${value}`);
}

function choiceFromOrderedMask(arr, mask) {
  const choice = new Uint8Array(arr.vertices.length);
  for (let bit = 0; bit < arr.crossingOrder.length; bit += 1) {
    if (((mask >> BigInt(bit)) & 1n) !== 0n) choice[arr.crossingOrder[bit]] = 1;
  }
  return choice;
}

function maskHex(mask, width) {
  return `0x${mask.toString(16).padStart(Math.ceil(width / 4), "0")}`;
}

function increment(object, key, amount = 1) {
  object[key] = (object[key] || 0) + amount;
}

function* masksWithHorizontalCount(arr, requiredH) {
  const horizontalBits = [];
  const verticalBits = [];
  for (let bit = 0; bit < arr.crossingOrder.length; bit += 1) {
    const vertex = arr.vertices[arr.crossingOrder[bit]];
    (vertex.orientation === "H" ? horizontalBits : verticalBits).push(bit);
  }
  if (requiredH < 0 || requiredH > horizontalBits.length) return;
  const bases = [];
  function choose(index, remaining, mask) {
    if (remaining === 0) {
      bases.push(mask);
      return;
    }
    if (horizontalBits.length - index < remaining) return;
    choose(index + 1, remaining - 1, mask | (1n << BigInt(horizontalBits[index])));
    choose(index + 1, remaining, mask);
  }
  choose(0, requiredH, 0n);
  const verticalLimit = 1n << BigInt(verticalBits.length);
  for (const base of bases) {
    for (let subset = 0n; subset < verticalLimit; subset += 1n) {
      let mask = base;
      for (let index = 0; index < verticalBits.length; index += 1) {
        if (((subset >> BigInt(index)) & 1n) !== 0n) {
          mask |= 1n << BigInt(verticalBits[index]);
        }
      }
      yield mask;
    }
  }
}

function auditAllDisk(search, arr, choice, options) {
  const trace = search.traceCircuits(arr, choice);
  if (trace === null) return { kind: "invalid", stage: "rich-trace" };
  if (trace.circuits.some((circuit) => circuit.kind === "essential")) {
    return { kind: "annular", trace };
  }
  const local = auditAllDiskLocalExact(search, arr, choice, trace, options);
  if (!local.ok) return { kind: "all-disk", ok: false, stage: local.stage };
  const balance = trace.circuits.reduce((sum, circuit) => (
    sum + BigInt(circuit.color)
      * translatedDiskMultiplicityExact(arr, circuit, options.basePoint)
  ), 0n);
  if (options.requireBalance !== false && balance !== BigInt(options.balanceTarget)) {
    return { kind: "all-disk", ok: false, stage: "base-unbalanced", balance };
  }
  return { kind: "all-disk", ok: true, trace, local, balance };
}

function run(payload) {
  const arr = hydrateArrangement(payload.arrangement);
  const diagnostics = diagnosticsObject();
  const search = require("./p2o4_annular_search.js");
  const audit = makeExactAuditor(
    search, diagnostics,
  );
  const options = {
    collectAnswers: payload.collectAnswers !== false,
    requireBalance: payload.requireBalance !== false,
    balanceTarget: payload.balanceTarget === undefined ? 0n : BigInt(payload.balanceTarget),
    basePoint: payload.basePoint,
    tau: payload.tau || [1, 0],
    sigma: payload.sigma || [0, 1],
    typeBEpsilon: payload.typeBEpsilon || "1/10",
    phaseSolverOptions: payload.phaseSolverOptions,
  };
  const size = arr.crossingOrder.length;
  if (size > 52) throw new Error("the direct mask driver is intended for at most 52 crossings");
  const requiredH = payload.requiredHorizontalSwitches === undefined
    ? -arr.targetChi : payload.requiredHorizontalSwitches;
  const fullMaskSpace = 2 ** size;
  const explicitMasks = payload.masks === undefined ? null : payload.masks.map(parseMask);
  const masks = explicitMasks === null
    ? masksWithHorizontalCount(arr, requiredH)
    : explicitMasks;
  const results = [];
  const rejected = Object.create(null);
  let eulerCompatibleMasks = 0;
  let pairableMasks = 0;
  let exactAllDiskRulingCount = 0;
  let exactAnnularRulingCount = 0;
  const allDiskResults = [];
  for (const mask of masks) {
    const choice = choiceFromOrderedMask(arr, mask);
    const horizontal = arr.vertices.reduce((count, vertex) => (
      count + (choice[vertex.index] && vertex.orientation === "H" ? 1 : 0)
    ), 0);
    if (requiredH !== null && horizontal !== requiredH) {
      increment(rejected, "wrong-horizontal-switch-count");
      continue;
    }
    eulerCompatibleMasks += 1;
    const diskAudit = auditAllDisk(search, arr, choice, options);
    if (diskAudit.kind === "all-disk") {
      if (!diskAudit.ok) {
        increment(rejected, `all-disk-${diskAudit.stage}`);
        continue;
      }
      exactAllDiskRulingCount += 1;
      allDiskResults.push({
        mask: maskHex(mask, size),
        diskEyes: diskAudit.local.disks,
        switches: diskAudit.local.switches,
        localTypes: diskAudit.local.localTypes,
        balance: diskAudit.balance,
        ribbonBoundaryComponents: diskAudit.local.ribbonBoundaryComponents,
        ribbonGenus: diskAudit.local.ribbonGenus,
      });
      continue;
    }
    if (diskAudit.kind === "invalid") {
      increment(rejected, diskAudit.stage);
      continue;
    }
    const result = audit(arr, choice, {
      ...options,
      precomputedTrace: diskAudit.trace,
    });
    if (!result.ok) {
      increment(rejected, result.stage);
      continue;
    }
    pairableMasks += 1;
    exactAnnularRulingCount += result.connected;
    for (const [stage, count] of Object.entries(result.rejected)) {
      increment(rejected, stage, count);
    }
    results.push({
      mask: maskHex(mask, size),
      diskEyes: result.diskEyes,
      essentialCircuits: result.essentialCircuits,
      annularEyes: result.annularEyes,
      matchingConfigurations: result.matchingConfigurations,
      phaseResidueConfigurations: result.colorConfigurations,
      accepted: result.connected,
      distribution: result.connectedDistribution,
      witnesses: result.witnesses,
    });
  }
  const phaseCardinality = diagnostics.infinite.length > 0
    ? "INFINITE" : diagnostics.unknown.length > 0 ? "UNKNOWN" : "FINITE";
  return {
    method: "direct-smoothing-plus-exact-integer-annular-phase",
    usesTropicalCount: false,
    finitePhaseCutoffUsed: false,
    masksRequested: explicitMasks === null ? fullMaskSpace : explicitMasks.length,
    masksActuallyTraced: eulerCompatibleMasks,
    eulerCompatibleMasks,
    pairableMasks,
    phaseCardinality,
    exactAllDiskRulingCount,
    exactAnnularRulingCount: phaseCardinality === "FINITE"
      ? exactAnnularRulingCount : phaseCardinality,
    exactTotalRulingCount: phaseCardinality === "FINITE"
      ? exactAllDiskRulingCount + exactAnnularRulingCount : phaseCardinality,
    allDiskResults,
    results,
    rejected,
    diagnostics,
  };
}

function main(argv) {
  const path = argv[0];
  const source = path ? fs.readFileSync(path, "utf8") : fs.readFileSync(0, "utf8");
  const payload = JSON.parse(source);
  process.stdout.write(`${JSON.stringify(run(payload), (_, entry) => (
    typeof entry === "bigint" ? entry.toString() : entry
  ), 2)}\n`);
}

if (require.main === module) {
  try { main(process.argv.slice(2)); }
  catch (error) {
    process.stderr.write(`${error.stack || error}\n`);
    process.exitCode = 1;
  }
}

module.exports = { hydrateArrangement, run };
