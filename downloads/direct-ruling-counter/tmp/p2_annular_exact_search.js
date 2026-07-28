#!/usr/bin/env node
"use strict";

/*
 * Independent no-K integration of the annular phase compiler.
 *
 * This file intentionally does not modify p2o4_annular_search.js.  It loads
 * that searcher's well-tested mask/DSU machinery, substitutes only the leaf
 * audit in memory, and solves every annular phase family as an exact integer
 * polyhedron.  Degree 4 is the local regression target.  Larger degrees
 * should be run on a remote machine.
 */

const {
  IntegerPolyhedron,
  compileAnnularEyeResidue,
  compileTypeBSameOrder,
  installCompiledConstraints,
  phaseResidueData,
  solveCompiledPhaseSystem,
  translatedDiskMultiplicityExact,
  stringifyExactResult,
} = require("./annular_phase_exact.js");

function loadSearchInternals() {
  return require("./p2o4_annular_search.js");
}

function safeNumber(value, label) {
  if (typeof value === "number") {
    if (!Number.isSafeInteger(value)) throw new Error(`${label} is not a safe integer`);
    return value;
  }
  const answer = Number(value);
  if (!Number.isSafeInteger(answer) || BigInt(answer) !== value) {
    throw new Error(`${label} exceeds the safe integer range`);
  }
  return answer;
}

function jsonInteger(value) {
  value = BigInt(value);
  const number = Number(value);
  return Number.isSafeInteger(number) && BigInt(number) === value
    ? number : value.toString();
}

function increment(object, key, amount = 1) {
  object[key] = (object[key] || 0) + amount;
}

function isUnknownStage(stage) {
  return typeof stage === "string" && stage.includes("unknown");
}

function residueProducts(qValues, callback, prefix = [], index = 0) {
  if (index === qValues.length) {
    callback(prefix.slice());
    return;
  }
  for (let residue = 0; residue < qValues[index]; residue += 1) {
    prefix.push(residue);
    residueProducts(qValues, callback, prefix, index + 1);
    prefix.pop();
  }
}

function compileAllTypeB(
  arr, trace, descriptors, staticData, options = {},
) {
  const annularByEye = new Map();
  for (let index = 0; index < descriptors.length; index += 1) {
    annularByEye.set(
      staticData.circuitToEye[descriptors[index].left], descriptors[index],
    );
  }
  const compiled = [];
  for (const deferred of staticData.deferredTypeB) {
    const result = compileTypeBSameOrder({
      arr,
      circuits: trace.circuits,
      vertex: deferred.vertex,
      circuitIds: deferred.circuitIds,
      eyeIds: deferred.ids,
      eyeCircuits: staticData.eyeCircuits,
      annularByEye,
      dimension: descriptors.length,
      tau: options.tau,
      sigma: options.sigma,
      epsilon: options.typeBEpsilon,
    });
    if (!result.ok) {
      return { ok: false, stage: result.stage, vertex: deferred.vertex.index };
    }
    compiled.push(result);
  }
  return { ok: true, constraints: compiled, countB: compiled.length };
}

/*
 * Exact all-disk local audit shared by the generic flat-mask and incremental
 * drivers.  localAudit supplies eye incidence, convexity, connectivity,
 * orientability, Euler, and ribbon checks.  Its legacy type-B order probe is
 * specialized to tau=2x+y, so defer every B switch and compile the same-order
 * tests with the caller's exact tau/sigma/epsilon instead.  With no annular
 * phase variables, the resulting polyhedron has dimension zero and contains
 * the single empty vector exactly when all constant inequalities hold.
 */
function auditAllDiskLocalExact(search, arr, choice, trace, options = {}) {
  if (trace.circuits.some((circuit) => circuit.kind !== "disk")) {
    return { ok: false, stage: "not-all-disk" };
  }
  const structural = search.localAudit(
    arr,
    choice,
    trace,
    [],
    0n,
    {
      ...options,
      structuralOnly: true,
      requireBalance: false,
      requireNormalB: false,
      requireAdjacentStrips: false,
      stripAudit: { ok: false, stage: "all-disk" },
      selfSwitchMode: "quotient",
    },
    null,
  );
  if (!structural.ok) return structural;

  const typeB = compileAllTypeB(arr, trace, [], structural, options);
  if (!typeB.ok) {
    return {
      ...typeB,
      stage: "non-normal-B",
      exactTypeBStage: typeB.stage,
    };
  }
  if (typeB.countB !== structural.localTypes.B) {
    throw new Error("exact all-disk type-B compiler missed a B switch");
  }
  const normality = new IntegerPolyhedron(0);
  for (const compiled of typeB.constraints) {
    installCompiledConstraints(normality, compiled.constraints);
  }
  const solved = normality.solve(options.phaseSolverOptions || {});
  if (solved.status !== "FINITE" || solved.count !== 1n) {
    if (solved.status === "UNKNOWN") {
      return {
        ok: false,
        stage: "exact-type-B-unknown",
        reason: solved.reason,
      };
    }
    return {
      ok: false,
      stage: "non-normal-B",
      exactTypeBStage: "constant-same-order-failed",
    };
  }
  return {
    ...structural,
    ok: true,
    stage: "connected-local-exact-all-disk",
    typeBCompanionOrders: {
      same: structural.localTypes.B,
      opposite: 0,
      inconsistent: 0,
    },
    ribbonBoundaryComponents: structural.ribbon.boundaryComponents,
    ribbonGenus: structural.ribbon.genus,
    ribbonCycleLengths: structural.ribbon.cycleLengths.slice(),
  };
}

function exactResultFromStructure(
  arr, choice, trace, matching, colors, descriptors, structural, phaseVector,
  balanceTarget,
) {
  const selfSwitchDecks = structural.selfSwitchDecks.map((entry) => {
    if (entry[1] !== "annular-phase") return entry.slice();
    const vertex = arr.vertices[entry[0]];
    const circuitId = trace.circuitOf[vertex.switchPairs[0][0]];
    const eyeId = structural.circuitToEye[circuitId];
    const annularIndex = eyeId - structural.disks;
    if (annularIndex < 0 || annularIndex >= phaseVector.length) {
      throw new Error("annular self-switch lost its exact phase variable");
    }
    return [entry[0], "annular-phase", jsonInteger(phaseVector[annularIndex])];
  });
  return {
    ok: true,
    stage: "connected-local-exact-phase",
    eyes: structural.eyes,
    disks: structural.disks,
    annuli: matching.length,
    switches: structural.switches,
    chi: structural.chi,
    baseBalance: jsonInteger(balanceTarget),
    localTypes: { ...structural.localTypes },
    typeBCompanionOrders: {
      same: structural.localTypes.B,
      opposite: 0,
      inconsistent: 0,
    },
    matching: matching.map((pair) => pair.slice()),
    colors: colors.map((color) => (color > 0 ? "B" : "R")),
    stripOk: false,
    stripStage: "exact-annular-phase-polyhedron",
    regionOk: true,
    annularPhases: phaseVector.map(jsonInteger),
    annularPhaseShifts: descriptors.map((descriptor, index) => (
      descriptor.phaseStep.map((coordinate) => jsonInteger(
        coordinate * phaseVector[index],
      ))
    )),
    ribbonBoundaryComponents: structural.ribbon.boundaryComponents,
    ribbonGenus: structural.ribbon.genus,
    ribbonCycleLengths: structural.ribbon.cycleLengths.slice(),
    selfSwitches: structural.selfSwitches,
    selfSwitchDecks,
  };
}

function choiceMaskHex(choice) {
  let mask = 0n;
  for (let index = 0; index < choice.length; index += 1) {
    if (choice[index]) mask |= 1n << BigInt(index);
  }
  return `0x${mask.toString(16)}`;
}

function reindexDescriptor(base, variable, dimension) {
  if (!base.ok) return base;
  const constraints = base.localConstraints.map((constraint) => {
    const coefficients = new Array(dimension).fill(0n);
    coefficients[variable] = constraint.coefficient;
    return { ...constraint, coefficients };
  });
  for (const excluded of base.forbiddenZ) {
    const coefficients = new Array(dimension).fill(0n);
    coefficients[variable] = 1n;
    constraints.push({
      coefficients,
      operator: "!=",
      rhs: excluded,
      label: "annular-boundaries-meet",
    });
  }
  return { ...base, variable, dimension, constraints };
}

function makeExactAuditor(search, globalDiagnostics) {
  return function auditMaskExact(arr, choice, options) {
    const trace = options.precomputedTrace || search.traceCircuits(arr, choice);
    if (trace === null) return { ok: false, stage: "rich-trace" };
    const diskCircuits = trace.circuits.filter((circuit) => circuit.kind === "disk");
    const essential = trace.circuits.filter((circuit) => circuit.kind === "essential");
    if (essential.length === 0) return { ok: false, stage: "all-disk" };
    const matchings = search.perfectMatchings(trace.circuits, null);
    if (matchings.length === 0) return { ok: false, stage: "unpairable" };
    const diskBase = diskCircuits.reduce((sum, circuit) => (
      sum + BigInt(circuit.color)
        * translatedDiskMultiplicityExact(arr, circuit, options.basePoint)
    ), 0n);
    const totals = {
      ok: true,
      stage: "pairable-exact-phase",
      diskEyes: diskCircuits.length,
      essentialCircuits: essential.length,
      annularEyes: essential.length / 2,
      matchingConfigurations: matchings.length,
      colorConfigurations: 0,
      boundaryDisjoint: 0,
      stripValid: 0,
      embeddedValid: 0,
      locallyValid: 0,
      connected: 0,
      rejected: Object.create(null),
      connectedDistribution: Object.create(null),
      witnesses: [],
    };

    // Pair/color/residue geometry is repeatedly reused across different
    // perfect matchings of the same essential circuits.
    const descriptorCache = new Map();

    for (const matching of matchings) {
      const qValues = matching.map((pair) => {
        const circuit = trace.circuits[pair[0]];
        const tau = (options.tau || [2, 1]).map((entry) => BigInt(entry));
        return safeNumber(
          phaseResidueData(
            BigInt(circuit.hx), BigInt(circuit.hy), tau[0], tau[1],
          ).q,
          "annular phase residue count",
        );
      });
      if (qValues.some((q) => !Number.isSafeInteger(q) || q <= 0)) {
        increment(totals.rejected, "annular-zero-sweep");
        continue;
      }
      const colorLimit = 1n << BigInt(matching.length);
      for (let colorMask = 0n; colorMask < colorLimit; colorMask += 1n) {
        const colors = matching.map((_, index) => (
          ((colorMask >> BigInt(index)) & 1n) === 0n ? -1 : 1
        ));
        // This is genuinely phase-independent and is therefore evaluated
        // once per matching/color, before the finite residue loop.
        const structural = search.localAudit(
          arr,
          choice,
          trace,
          matching,
          colorMask,
          {
            ...options,
            structuralOnly: true,
            requireBalance: false,
            requireNormalB: false,
            requireAdjacentStrips: false,
            stripAudit: { ok: false, stage: "structural-only" },
            selfSwitchMode: "quotient",
          },
          null,
        );
        if (!structural.ok) {
          increment(totals.rejected, structural.stage);
          continue;
        }
        residueProducts(qValues, (residues) => {
          totals.colorConfigurations += 1;
          const descriptors = [];
          for (let eye = 0; eye < matching.length; eye += 1) {
            const cacheKey = `${matching[eye][0]},${matching[eye][1]}`
              + `:${colors[eye]}:${residues[eye]}`;
            let base = descriptorCache.get(cacheKey);
            if (base === undefined) {
              base = compileAnnularEyeResidue({
                arr,
                circuits: trace.circuits,
                pair: matching[eye],
                color: colors[eye],
                residue: residues[eye],
                variable: 0,
                dimension: 1,
                tau: options.tau,
                sigma: options.sigma,
                basePoint: options.basePoint,
                solverOptions: options.phaseSolverOptions,
              });
              descriptorCache.set(cacheKey, base);
            }
            const descriptor = reindexDescriptor(base, eye, matching.length);
            if (!descriptor.ok) {
              if (isUnknownStage(descriptor.stage)) {
                globalDiagnostics.unknown.push({
                  maskHex: choiceMaskHex(choice),
                  matching: matching.map((pair) => pair.slice()),
                  colors: colors.slice(),
                  residues: residues.slice(),
                  eye,
                  stage: descriptor.stage,
                  reason: descriptor.reason || descriptor.stage,
                });
                increment(totals.rejected, "exact-phase-unknown");
              } else increment(totals.rejected, descriptor.stage);
              return;
            }
            descriptors.push(descriptor);
          }
          totals.boundaryDisjoint += 1;
          totals.embeddedValid += 1;
          const typeB = compileAllTypeB(
            arr, trace, descriptors, structural, options,
          );
          if (!typeB.ok) {
            increment(totals.rejected, typeB.stage);
            return;
          }
          if (typeB.countB !== structural.localTypes.B) {
            throw new Error("exact type-B compiler missed a static B switch");
          }
          const solved = solveCompiledPhaseSystem({
            descriptors,
            typeB: typeB.constraints,
            diskContribution: diskBase,
            balanceTarget: options.balanceTarget === undefined
              ? 0n : BigInt(options.balanceTarget),
            requireBalance: options.requireBalance !== false,
            solverOptions: options.phaseSolverOptions,
          });
          globalDiagnostics.solver.calls += 1;
          globalDiagnostics.solver.activeSets += solved.activeSets || 0;
          globalDiagnostics.solver.integerNodes += solved.integerNodes || 0;
          globalDiagnostics.solver.maxActiveSets = Math.max(
            globalDiagnostics.solver.maxActiveSets, solved.activeSets || 0,
          );
          globalDiagnostics.solver.maxIntegerNodes = Math.max(
            globalDiagnostics.solver.maxIntegerNodes, solved.integerNodes || 0,
          );
          globalDiagnostics.solver[solved.status.toLowerCase()] += 1;
          if (solved.status === "UNKNOWN") {
            globalDiagnostics.unknown.push({
              maskHex: choiceMaskHex(choice),
              matching: matching.map((pair) => pair.slice()),
              colors: colors.slice(),
              residues: residues.slice(),
              reason: solved.reason,
            });
            increment(totals.rejected, "exact-phase-unknown");
            return;
          }
          if (solved.status === "INFINITE") {
            globalDiagnostics.infinite.push({
              maskHex: choiceMaskHex(choice),
              matching: matching.map((pair) => pair.slice()),
              colors: colors.slice(),
              residues: residues.slice(),
              point: solved.point.slice(),
              ray: solved.ray.slice(),
            });
            increment(totals.rejected, "infinite-phase-family");
            return;
          }
          if (solved.status === "EMPTY") {
            increment(totals.rejected, "exact-phase-empty");
            return;
          }
          if (!solved.collectedAll) {
            globalDiagnostics.unknown.push({ reason: "finite phase points not all collected" });
            increment(totals.rejected, "exact-phase-unknown");
            return;
          }
          if (BigInt(solved.points.length) !== solved.count) {
            throw new Error("finite exact phase count/point list mismatch");
          }
          for (const point of solved.points) {
            const phaseVector = point.map((z, eye) => (
              descriptors[eye].residue + descriptors[eye].q * z
            ));
            const finalResult = exactResultFromStructure(
              arr, choice, trace, matching, colors, descriptors, structural, phaseVector,
              options.balanceTarget === undefined ? 0n : BigInt(options.balanceTarget),
            );
            // The old finite-K audit is useful as a regression cross-check, but it
            // contains a floating-point type-B probe and is therefore never
            // part of the exact acceptance decision.
            if (options.verifyLegacyPhaseAudit === true) {
              const legacyPhaseVector = phaseVector.map((phase) => (
                safeNumber(phase, "legacy exact phase")
              ));
              const legacy = search.localAudit(
                arr,
                choice,
                trace,
                matching,
                colorMask,
                {
                  ...options,
                  // localAudit's legacy balance target is hard-coded to 0.
                  // The exact system above has already checked a nonzero d.
                  requireBalance: (options.balanceTarget === undefined
                    || BigInt(options.balanceTarget) === 0n)
                    && options.requireBalance !== false,
                  selfSwitchMode: "quotient",
                },
                legacyPhaseVector,
              );
              if (!legacy.ok) {
                throw new Error(`exact compiler/legacy audit mismatch: ${legacy.stage}`);
              }
            }
            totals.locallyValid += 1;
            totals.connected += 1;
            const local = finalResult.localTypes;
            const key = `D${finalResult.disks}-R${finalResult.annuli}`
              + `-A${local.A}-B${local.B}-V${local.V}-S${finalResult.selfSwitches}`;
            increment(totals.connectedDistribution, key);
            if (options.collectAnswers === true || totals.witnesses.length < 4) {
              totals.witnesses.push(finalResult);
            }
          }
        });
      }
    }
    return totals;
  };
}

function runExact(options = {}) {
  const search = loadSearchInternals();
  const diagnostics = {
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
  const auditMaskOverride = makeExactAuditor(search, diagnostics);
  const result = search.enumerate({
    degree: options.degree === undefined ? 4 : options.degree,
    genus: options.genus === undefined ? 0 : options.genus,
    order: options.order || "yrev",
    translationSymmetry: options.translationSymmetry !== false,
    selfSwitchMode: "quotient",
    collectAnswers: options.collectAnswers !== false,
    maxNodes: options.maxNodes,
    prefix: options.prefix || "",
    requireBalance: options.requireBalance !== false,
    essentialPairingPruning: options.essentialPairingPruning !== false,
    verifyIncremental: options.verifyIncremental === true,
    phaseSolverOptions: options.phaseSolverOptions,
    auditMaskOverride,
  });
  result.phaseMode = "exact-no-K";
  result.countScope = "annular-rulings-only";
  result.finitePhaseCutoffUsed = false;
  result.strictCompleteScope = "encoded-p2-annular-model";
  result.outerDefinitionEquivalence = "not-yet-formally-proved";
  result.definitionComplete = false;
  result.definitionCompletenessReason =
    "the general Newton-polygon adapter and the local/annular geometric equivalence lemmas are not yet formalized";
  result.balanceTarget = options.balanceTarget === undefined
    ? 0 : jsonInteger(options.balanceTarget);
  result.strictlyDecided = diagnostics.unknown.length === 0 && !result.aborted;
  result.strictComplete = diagnostics.infinite.length === 0
    && diagnostics.unknown.length === 0
    && !result.aborted;
  result.phaseSearchExhaustive = diagnostics.unknown.length === 0
    && !result.aborted;
  result.encodedModelComplete = result.strictComplete;
  result.infinitePhaseFamilies = diagnostics.infinite;
  result.unknownPhaseFamilies = diagnostics.unknown;
  result.exactPhaseSolver = diagnostics.solver;
  if (diagnostics.infinite.length !== 0) {
    result.phaseCardinality = "INFINITE";
    result.exactAnnularRulingCount = "INFINITE";
    result.exactRulingCount = "INFINITE";
  } else if (diagnostics.unknown.length !== 0 || result.aborted) {
    result.phaseCardinality = "UNKNOWN";
    result.exactAnnularRulingCount = "UNKNOWN";
    result.exactRulingCount = "UNKNOWN";
  } else {
    result.phaseCardinality = "FINITE";
    result.exactAnnularRulingCount = result.connectedConfigurations;
    result.exactRulingCount = result.connectedConfigurations;
  }
  return result;
}

function main(argv) {
  const options = {};
  for (let index = 0; index < argv.length; index += 1) {
    const argument = argv[index];
    if (argument === "--degree") options.degree = Number.parseInt(argv[++index], 10);
    else if (argument === "--genus") options.genus = Number.parseInt(argv[++index], 10);
    else if (argument === "--order") options.order = argv[++index];
    else if (argument === "--prefix") options.prefix = argv[++index];
    else if (argument === "--max-nodes") {
      options.maxNodes = Number.parseInt(argv[++index], 10);
    } else if (argument === "--balance-target") {
      options.balanceTarget = BigInt(argv[++index]);
    } else if (argument === "--no-translation-symmetry") {
      options.translationSymmetry = false;
    } else if (argument === "--no-answers") options.collectAnswers = false;
    else if (argument === "--help" || argument === "-h") {
      process.stdout.write(
        "Usage: node tmp/p2_annular_exact_search.js [--degree 4] [--genus 0] "
          + "[--order yrev] [--prefix 010] [--max-nodes N] "
          + "[--balance-target D] "
          + "[--no-translation-symmetry] [--no-answers]\n",
      );
      return;
    } else throw new Error(`unknown argument ${argument}`);
  }
  if ((options.degree || 4) > 4) {
    process.stderr.write(
      "warning: degree > 4 can be expensive; run the full job on the VPS\n",
    );
  }
  process.stdout.write(`${stringifyExactResult(runExact(options))}\n`);
}

if (require.main === module) {
  try {
    main(process.argv.slice(2));
  } catch (error) {
    process.stderr.write(`${error.stack || error}\n`);
    process.exitCode = 1;
  }
}

module.exports = { auditAllDiskLocalExact, runExact, makeExactAuditor };
