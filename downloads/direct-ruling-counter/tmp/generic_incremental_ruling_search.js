#!/usr/bin/env node
"use strict";

/*
 * Standalone prototype: rollback search for a JSON-supplied primitive
 * multigeodesic arrangement.
 *
 * This intentionally does not modify either ruling_smoothing.py or the
 * production annular search.  External geodesic edges are installed in a
 * Z^2-voltage DSU once.  A DFS then adds the two local smoothing arcs at one
 * crossing at a time.  A component is checked as soon as it closes:
 *
 *   contractible  => exactly two vertical tips;
 *   essential     => zero vertical tips and primitive homology.
 *
 * Active zero-lift inequalities reject a nonswitch self-intersection or a
 * same-lift self-switch as soon as its two ports become connected.  Closed
 * essential circuits are also balanced incrementally by primitive homology
 * and canonical coorientation side.  Only surviving leaves are passed to
 * the existing full disk/annular auditors.
 */

const fs = require("fs");
const {
  hydrateArrangement,
} = require("./generic_annular_exact_search.js");
const {
  auditAllDiskLocalExact,
  makeExactAuditor,
} = require("./p2_annular_exact_search.js");
const core = require("./p2o4_annular_search.js");
const {
  translatedDiskMultiplicityExact,
} = require("./annular_phase_exact.js");

function gcd(left, right) {
  left = Math.abs(left);
  right = Math.abs(right);
  while (right !== 0) [left, right] = [right, left % right];
  return left;
}

class RollbackVoltageDSU {
  constructor(size) {
    this.parent = Int16Array.from({ length: size }, (_, index) => index);
    this.componentSize = new Int16Array(size);
    this.componentSize.fill(1);
    this.potentialX = new Int16Array(size);
    this.potentialY = new Int16Array(size);
    this.tips = new Int8Array(size);
    // 0=open, 1=closed disk, 2=closed essential.
    this.closedKind = new Int8Array(size);
    this.closedHx = new Int16Array(size);
    this.closedHy = new Int16Array(size);
    this.closedSide = new Int8Array(size);
    this.zeroMate = new Int16Array(size);
    this.zeroKind = new Int8Array(size);
    this.zeroMate.fill(-1);
    this.stack = [];
    this.components = size;
    this.closedCycles = 0;
    this.diskCycles = 0;
    this.essentialCycles = 0;
    this.lastFailure = null;
    this.lastClosedRoot = -1;
  }

  mark() {
    return this.stack.length;
  }

  find(node) {
    let x = 0;
    let y = 0;
    let current = node;
    while (this.parent[current] !== current) {
      x += this.potentialX[current];
      y += this.potentialY[current];
      current = this.parent[current];
    }
    return [current, x, y];
  }

  relation(left, right) {
    const a = this.find(left);
    const b = this.find(right);
    if (a[0] !== b[0]) return null;
    return [b[1] - a[1], b[2] - a[2]];
  }

  scanZeroConstraints() {
    for (let port = 0; port < this.zeroMate.length; port += 1) {
      const mate = this.zeroMate[port];
      if (mate < port) continue;
      const relation = this.relation(port, mate);
      if (relation === null || relation[0] !== 0 || relation[1] !== 0) continue;
      this.lastFailure = this.zeroKind[port] === 2
        ? "same-lift-self-switch" : "nonswitch-self-intersection";
      return false;
    }
    return true;
  }

  activateZeroConstraint(left, right, kind) {
    if (this.zeroMate[left] !== -1 || this.zeroMate[right] !== -1) {
      throw new Error("a port acquired two zero-lift constraints");
    }
    this.zeroMate[left] = right;
    this.zeroMate[right] = left;
    this.zeroKind[left] = kind;
    this.zeroKind[right] = kind;
    return this.scanZeroConstraints();
  }

  deactivateZeroConstraint(left, right) {
    if (this.zeroMate[left] !== right || this.zeroMate[right] !== left) {
      throw new Error("zero-lift constraint rollback mismatch");
    }
    this.zeroMate[left] = -1;
    this.zeroMate[right] = -1;
    this.zeroKind[left] = 0;
    this.zeroKind[right] = 0;
  }

  unite(left, right, dx, dy, addedTips = 0) {
    const a = this.find(left);
    const b = this.find(right);
    let root = a[0];
    let child = b[0];
    let relativeX = a[1] + dx - b[1];
    let relativeY = a[2] + dy - b[2];
    this.lastClosedRoot = -1;

    if (root === child) {
      if (this.closedKind[root] !== 0) {
        this.lastFailure = "duplicate-cycle";
        return false;
      }
      const tipCount = this.tips[root] + addedTips;
      let kind;
      if (relativeX === 0 && relativeY === 0) {
        if (tipCount !== 2) {
          this.lastFailure = "disk-wrong-tip-count";
          return false;
        }
        kind = 1;
      } else {
        if (tipCount !== 0) {
          this.lastFailure = "essential-has-tip";
          return false;
        }
        if (gcd(relativeX, relativeY) !== 1) {
          this.lastFailure = "essential-nonprimitive";
          return false;
        }
        kind = 2;
      }
      this.stack.push({ kind: "close", root, oldTips: this.tips[root] });
      this.tips[root] = tipCount;
      this.closedKind[root] = kind;
      if (kind === 2) {
        let hx = relativeX;
        let hy = relativeY;
        let side = 1;
        if (hx < 0 || (hx === 0 && hy < 0)) {
          hx = -hx;
          hy = -hy;
          side = -1;
        }
        this.closedHx[root] = hx;
        this.closedHy[root] = hy;
        this.closedSide[root] = side;
        this.essentialCycles += 1;
      } else {
        this.diskCycles += 1;
      }
      this.closedCycles += 1;
      this.lastClosedRoot = root;
      return true;
    }

    if (this.closedKind[root] !== 0 || this.closedKind[child] !== 0) {
      this.lastFailure = "join-closed-cycle";
      return false;
    }
    const tipCount = this.tips[root] + this.tips[child] + addedTips;
    if (tipCount > 2) {
      this.lastFailure = "tip-overflow";
      return false;
    }
    if (this.componentSize[root] < this.componentSize[child]) {
      [root, child] = [child, root];
      relativeX = -relativeX;
      relativeY = -relativeY;
    }
    this.stack.push({
      kind: "union",
      root,
      child,
      oldSize: this.componentSize[root],
      oldTips: this.tips[root],
    });
    this.parent[child] = root;
    this.potentialX[child] = relativeX;
    this.potentialY[child] = relativeY;
    this.componentSize[root] += this.componentSize[child];
    this.tips[root] = tipCount;
    this.components -= 1;
    return this.scanZeroConstraints();
  }

  rollback(mark) {
    while (this.stack.length > mark) {
      const record = this.stack.pop();
      if (record.kind === "close") {
        const kind = this.closedKind[record.root];
        if (kind === 1) this.diskCycles -= 1;
        else this.essentialCycles -= 1;
        this.closedCycles -= 1;
        this.tips[record.root] = record.oldTips;
        this.closedKind[record.root] = 0;
        this.closedHx[record.root] = 0;
        this.closedHy[record.root] = 0;
        this.closedSide[record.root] = 0;
      } else {
        this.parent[record.child] = record.child;
        this.potentialX[record.child] = 0;
        this.potentialY[record.child] = 0;
        this.componentSize[record.root] = record.oldSize;
        this.tips[record.root] = record.oldTips;
        this.components += 1;
      }
    }
  }

  openTipCounts() {
    const counts = [0, 0, 0];
    for (let node = 0; node < this.parent.length; node += 1) {
      if (this.parent[node] === node && this.closedKind[node] === 0) {
        counts[this.tips[node]] += 1;
      }
    }
    return counts;
  }
}

function addExternalEdges(arr, dsu) {
  for (const edge of arr.edges) {
    if (!dsu.unite(edge.a, edge.b, edge.sx, edge.sy)) {
      throw new Error(`external edge failed: ${dsu.lastFailure}`);
    }
  }
  dsu.stack.length = 0;
}

function addLocal(arr, dsu, vertex, switched, closedRoots) {
  const pairs = switched ? vertex.switchPairs : vertex.nonswitchPairs;
  const tip = switched && vertex.orientation === "V" ? 1 : 0;
  for (const [left, right] of pairs) {
    // The adapter guarantees opposite coorientation signs in every local
    // pair.  Orienting from sign - to sign + fixes traceSide=+1.
    const incoming = arr.portSign[left] < 0 ? left : right;
    const outgoing = incoming === left ? right : left;
    const before = dsu.closedCycles;
    if (!dsu.unite(incoming, outgoing, 0, 0, tip)) return false;
    if (dsu.closedCycles !== before) closedRoots.push(dsu.lastClosedRoot);
  }
  return true;
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

function exactAllDiskAudit(arr, choice, trace, payload) {
  if (trace.circuits.some((circuit) => circuit.kind !== "disk")) {
    return { ok: false, stage: "not-all-disk" };
  }
  const structural = auditAllDiskLocalExact(core, arr, choice, trace, {
    tau: payload.tau,
    sigma: payload.sigma,
    typeBEpsilon: payload.typeBEpsilon,
    phaseSolverOptions: payload.phaseSolverOptions,
  });
  if (!structural.ok) return structural;
  let balance = 0n;
  for (const circuit of trace.circuits) {
    balance += BigInt(circuit.color)
      * translatedDiskMultiplicityExact(arr, circuit, payload.basePoint);
  }
  const target = payload.balanceTarget === undefined ? 0n : BigInt(payload.balanceTarget);
  if (balance !== target) return { ok: false, stage: "base-unbalanced", balance };
  return { ok: true, structural, balance };
}

function orderFromName(arr, name) {
  const crossing = arr.crossingOrder.slice();
  if (name === "crossing") return crossing;
  if (name === "reverse") return crossing.reverse();
  if (name === "hfirst") {
    return crossing.sort((left, right) => (
      (arr.vertices[left].orientation === "H" ? 0 : 1)
      - (arr.vertices[right].orientation === "H" ? 0 : 1)
    ));
  }
  if (name === "vfirst") {
    return crossing.sort((left, right) => (
      (arr.vertices[left].orientation === "V" ? 0 : 1)
      - (arr.vertices[right].orientation === "V" ? 0 : 1)
    ));
  }
  throw new Error(`unknown DFS order ${name}`);
}

function run(payload, options = {}) {
  const arr = hydrateArrangement(payload.arrangement);
  const order = orderFromName(arr, options.order || "crossing");
  const requiredH = payload.requiredHorizontalSwitches === undefined
    ? -arr.targetChi : Number(payload.requiredHorizontalSwitches);
  const remainingH = new Int8Array(order.length + 1);
  const remainingV = new Int8Array(order.length + 1);
  for (let index = order.length - 1; index >= 0; index -= 1) {
    const horizontal = arr.vertices[order[index]].orientation === "H";
    remainingH[index] = remainingH[index + 1] + (horizontal ? 1 : 0);
    remainingV[index] = remainingV[index + 1] + (horizontal ? 0 : 1);
  }

  const diagnostics = diagnosticsObject();
  const annularAudit = makeExactAuditor(core, diagnostics);
  const collectAnswers = payload.collectAnswers === true;
  const annularOptions = {
    collectAnswers,
    requireBalance: payload.requireBalance !== false,
    balanceTarget: payload.balanceTarget === undefined ? 0n : BigInt(payload.balanceTarget),
    basePoint: payload.basePoint,
    tau: payload.tau || [1, 0],
    sigma: payload.sigma || [0, 1],
    typeBEpsilon: payload.typeBEpsilon || "1/10",
    phaseSolverOptions: payload.phaseSolverOptions,
  };

  const dsu = new RollbackVoltageDSU(arr.portCount);
  addExternalEdges(arr, dsu);
  const choice = new Uint8Array(arr.vertices.length);
  const essentialBalance = new Map();
  let essentialImbalance = 0;
  const stats = {
    order: options.order || "crossing",
    searchNodes: 0,
    branchAttempts: 0,
    fullTraceCalls: 0,
    completedLeaves: 0,
    allDiskLeaves: 0,
    annularLeaves: 0,
    allDiskRulings: 0,
    annularRulings: 0,
    allDiskMasks: [],
    annularResults: [],
    earlyRejected: Object.create(null),
    leafRejected: Object.create(null),
    diagnostics,
    elapsedMs: 0,
  };
  const started = Date.now();
  const count = (table, key, amount = 1) => {
    table[key] = (table[key] || 0) + amount;
  };
  function addEssentialRoots(roots) {
    const signatures = [];
    for (const root of roots) {
      if (dsu.closedKind[root] !== 2) continue;
      const key = `${dsu.closedHx[root]},${dsu.closedHy[root]}`;
      const side = dsu.closedSide[root];
      const before = essentialBalance.get(key) || 0;
      const after = before + side;
      essentialImbalance += Math.abs(after) - Math.abs(before);
      if (after === 0) essentialBalance.delete(key);
      else essentialBalance.set(key, after);
      signatures.push([key, side]);
    }
    return signatures;
  }
  function removeEssentialRoots(signatures) {
    for (let index = signatures.length - 1; index >= 0; index -= 1) {
      const [key, side] = signatures[index];
      const before = essentialBalance.get(key) || 0;
      const after = before - side;
      essentialImbalance += Math.abs(after) - Math.abs(before);
      if (after === 0) essentialBalance.delete(key);
      else essentialBalance.set(key, after);
    }
  }

  function visit(index, horizontalSwitches, verticalSwitches, mask) {
    stats.searchNodes += 1;
    if (horizontalSwitches > requiredH
        || horizontalSwitches + remainingH[index] < requiredH) {
      count(stats.earlyRejected, "horizontal-count-bound");
      return;
    }
    if (dsu.diskCycles > verticalSwitches
        || dsu.diskCycles > verticalSwitches + remainingV[index]) {
      count(stats.earlyRejected, "disk-V-bound");
      return;
    }
    const openTips = dsu.openTipCounts();
    if (essentialImbalance > openTips[0]) {
      count(stats.earlyRejected, "essential-pairing-bound");
      return;
    }
    if (index === order.length) {
      stats.completedLeaves += 1;
      if (dsu.components !== dsu.closedCycles) {
        count(stats.leafRejected, "incomplete-cycle-cover");
        return;
      }
      if (dsu.diskCycles !== verticalSwitches) {
        count(stats.leafRejected, "disk-V-mismatch");
        return;
      }
      if (essentialImbalance !== 0) {
        count(stats.leafRejected, "essential-unpairable");
        return;
      }
      stats.fullTraceCalls += 1;
      const trace = core.traceCircuits(arr, choice);
      if (trace === null) {
        throw new Error(`incremental survivor ${mask.toString(16)} failed full trace`);
      }
      if (dsu.essentialCycles === 0) {
        stats.allDiskLeaves += 1;
        const result = exactAllDiskAudit(arr, choice, trace, payload);
        if (!result.ok) {
          count(stats.leafRejected, `all-disk-${result.stage}`);
          return;
        }
        stats.allDiskRulings += 1;
        stats.allDiskMasks.push(`0x${mask.toString(16).padStart(4, "0")}`);
      } else {
        stats.annularLeaves += 1;
        const result = annularAudit(arr, choice, {
          ...annularOptions,
          precomputedTrace: trace,
        });
        if (!result.ok) {
          count(stats.leafRejected, result.stage);
          return;
        }
        stats.annularRulings += result.connected;
        if (collectAnswers && result.connected > 0) {
          stats.annularResults.push({
            mask: `0x${mask.toString(16).padStart(Math.ceil(arr.crossingOrder.length / 4), "0")}`,
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
        for (const [stage, amount] of Object.entries(result.rejected)) {
          count(stats.leafRejected, stage, amount);
        }
      }
      return;
    }

    const vertexId = order[index];
    const vertex = arr.vertices[vertexId];
    const horizontal = vertex.orientation === "H";
    const mustSwitch = horizontal
      && horizontalSwitches + remainingH[index] === requiredH;
    const cannotSwitch = horizontal && horizontalSwitches === requiredH;
    for (const switched of [false, true]) {
      if ((switched && cannotSwitch) || (!switched && mustSwitch)) continue;
      stats.branchAttempts += 1;
      const mark = dsu.mark();
      choice[vertexId] = switched ? 1 : 0;
      const pairs = switched ? vertex.switchPairs : vertex.nonswitchPairs;
      const zeroLeft = pairs[0][0];
      const zeroRight = pairs[1][0];
      const zeroKind = switched ? 2 : 1;
      const closedRoots = [];
      let signatures = [];
      let okay = dsu.activateZeroConstraint(zeroLeft, zeroRight, zeroKind);
      if (!okay) count(stats.earlyRejected, dsu.lastFailure);
      if (okay) {
        okay = addLocal(arr, dsu, vertex, switched, closedRoots);
        if (!okay) count(stats.earlyRejected, dsu.lastFailure);
      }
      if (okay) {
        signatures = addEssentialRoots(closedRoots);
        const nextMask = switched
          ? mask | (1n << BigInt(arr.crossingOrder.indexOf(vertexId))) : mask;
        visit(
          index + 1,
          horizontalSwitches + (switched && horizontal ? 1 : 0),
          verticalSwitches + (switched && !horizontal ? 1 : 0),
          nextMask,
        );
      }
      removeEssentialRoots(signatures);
      dsu.rollback(mark);
      dsu.deactivateZeroConstraint(zeroLeft, zeroRight);
      choice[vertexId] = 0;
    }
  }

  visit(0, 0, 0, 0n);
  stats.elapsedMs = Date.now() - started;
  const phaseCardinality = diagnostics.infinite.length > 0
    ? "INFINITE" : diagnostics.unknown.length > 0 ? "UNKNOWN" : "FINITE";
  stats.method = "incremental-rollback-smoothing-plus-exact-integer-annular-phase";
  stats.usesTropicalCount = false;
  stats.finitePhaseCutoffUsed = false;
  stats.phaseCardinality = phaseCardinality;
  stats.exactAllDiskRulingCount = stats.allDiskRulings;
  stats.exactAnnularRulingCount = phaseCardinality === "FINITE"
    ? stats.annularRulings : phaseCardinality;
  stats.exactTotalRulingCount = phaseCardinality === "FINITE"
    ? stats.allDiskRulings + stats.annularRulings : phaseCardinality;
  return stats;
}

function main(argv) {
  let order = "crossing";
  let path = null;
  for (let index = 0; index < argv.length; index += 1) {
    if (argv[index] === "--order") order = argv[++index];
    else if (path === null) path = argv[index];
    else throw new Error(`unknown argument ${argv[index]}`);
  }
  const source = path ? fs.readFileSync(path, "utf8") : fs.readFileSync(0, "utf8");
  const result = run(JSON.parse(source), { order });
  process.stdout.write(`${JSON.stringify(result, (_, entry) => (
    typeof entry === "bigint" ? entry.toString() : entry
  ), 2)}\n`);
}

if (require.main === module) {
  try {
    main(process.argv.slice(2));
  } catch (error) {
    process.stderr.write(`${error.stack || error}\n`);
    process.exitCode = 1;
  }
}

module.exports = { RollbackVoltageDSU, run };
