"use strict";

/*
 * Exact integer-polyhedron helper for annular deck phases.
 *
 * The intended input is the constraint system obtained after fixing an
 * annular matching, its colors, and the finite residue vector
 *
 *     k_i = r_i + q_i z_i,       z_i in Z.
 *
 * All arithmetic below is BigInt/rational arithmetic.  There is no floating
 * point tolerance and, in particular, no artificial box |z_i| <= K.
 *
 * This is deliberately a small-dimensional solver.  It finds recession
 * directions and vertices by exact active-set enumeration.  That is a good
 * fit for the phase problem (the dimension is the number of annular eyes),
 * but it is not meant to compete with a general-purpose MILP package in high
 * dimension.  Resource limits produce UNKNOWN, never a finite/infinite claim.
 */

class PolyhedralLimitError extends Error {
  constructor(message) {
    super(message);
    this.name = "PolyhedralLimitError";
  }
}

function absBigInt(value) {
  return value < 0n ? -value : value;
}

function gcdBigInt(left, right) {
  left = absBigInt(left);
  right = absBigInt(right);
  while (right !== 0n) [left, right] = [right, left % right];
  return left;
}

function lcmBigInt(left, right) {
  if (left === 0n || right === 0n) return 0n;
  return absBigInt((left / gcdBigInt(left, right)) * right);
}

function floorDiv(numerator, denominator) {
  if (denominator <= 0n) throw new Error("floorDiv needs a positive denominator");
  let quotient = numerator / denominator;
  if (numerator < 0n && numerator % denominator !== 0n) quotient -= 1n;
  return quotient;
}

function ceilDiv(numerator, denominator) {
  if (denominator <= 0n) throw new Error("ceilDiv needs a positive denominator");
  return -floorDiv(-numerator, denominator);
}

function makeFraction(numerator, denominator = 1n) {
  numerator = BigInt(numerator);
  denominator = BigInt(denominator);
  if (denominator === 0n) throw new Error("zero rational denominator");
  if (denominator < 0n) {
    numerator = -numerator;
    denominator = -denominator;
  }
  const divisor = gcdBigInt(numerator, denominator);
  return { n: numerator / divisor, d: denominator / divisor };
}

function asFraction(value) {
  if (typeof value === "bigint") return { n: value, d: 1n };
  if (typeof value === "number") {
    if (!Number.isSafeInteger(value)) {
      throw new Error(`only safe integer Numbers are exact, got ${value}`);
    }
    return { n: BigInt(value), d: 1n };
  }
  if (typeof value === "string") {
    const pieces = value.split("/");
    if (pieces.length === 1) return makeFraction(BigInt(pieces[0]));
    if (pieces.length === 2) return makeFraction(BigInt(pieces[0]), BigInt(pieces[1]));
    throw new Error(`bad rational string ${value}`);
  }
  if (value !== null && typeof value === "object"
      && value.n !== undefined && value.d !== undefined) {
    return makeFraction(value.n, value.d);
  }
  throw new Error(`cannot convert ${String(value)} to an exact rational`);
}

function fractionAdd(left, right) {
  return makeFraction(left.n * right.d + right.n * left.d, left.d * right.d);
}

function fractionSubtract(left, right) {
  return makeFraction(left.n * right.d - right.n * left.d, left.d * right.d);
}

function fractionMultiply(left, right) {
  return makeFraction(left.n * right.n, left.d * right.d);
}

function fractionDivide(left, right) {
  if (right.n === 0n) throw new Error("division by zero rational");
  return makeFraction(left.n * right.d, left.d * right.n);
}

function fractionNegate(value) {
  return { n: -value.n, d: value.d };
}

function fractionCompare(left, right) {
  const difference = left.n * right.d - right.n * left.d;
  return difference < 0n ? -1 : difference > 0n ? 1 : 0;
}

function fractionIsZero(value) {
  return value.n === 0n;
}

function fractionFloor(value) {
  return floorDiv(value.n, value.d);
}

function fractionCeil(value) {
  return ceilDiv(value.n, value.d);
}

function zeroFraction() {
  return { n: 0n, d: 1n };
}

function dotBigInt(coefficients, point) {
  let answer = 0n;
  for (let index = 0; index < coefficients.length; index += 1) {
    answer += coefficients[index] * point[index];
  }
  return answer;
}

function dotFraction(coefficients, point) {
  let answer = zeroFraction();
  for (let index = 0; index < coefficients.length; index += 1) {
    if (coefficients[index] === 0n) continue;
    answer = fractionAdd(
      answer,
      fractionMultiply(makeFraction(coefficients[index]), point[index]),
    );
  }
  return answer;
}

function integerVectorFromRational(vector) {
  let denominator = 1n;
  for (const entry of vector) denominator = lcmBigInt(denominator, entry.d);
  const integers = vector.map((entry) => entry.n * (denominator / entry.d));
  let divisor = 0n;
  for (const entry of integers) divisor = gcdBigInt(divisor, entry);
  if (divisor === 0n) return integers;
  const primitive = integers.map((entry) => entry / divisor);
  const first = primitive.find((entry) => entry !== 0n);
  return first !== undefined && first < 0n
    ? primitive.map((entry) => -entry)
    : primitive;
}

function cloneFractionMatrix(matrix) {
  return matrix.map((row) => row.map((entry) => ({ n: entry.n, d: entry.d })));
}

function rref(input, coefficientColumns) {
  const matrix = cloneFractionMatrix(input);
  const rows = matrix.length;
  const columns = rows === 0 ? coefficientColumns : matrix[0].length;
  const pivots = [];
  let pivotRow = 0;
  for (let column = 0; column < coefficientColumns && pivotRow < rows; column += 1) {
    let selected = pivotRow;
    while (selected < rows && fractionIsZero(matrix[selected][column])) selected += 1;
    if (selected === rows) continue;
    [matrix[pivotRow], matrix[selected]] = [matrix[selected], matrix[pivotRow]];
    const pivot = matrix[pivotRow][column];
    for (let j = column; j < columns; j += 1) {
      matrix[pivotRow][j] = fractionDivide(matrix[pivotRow][j], pivot);
    }
    for (let row = 0; row < rows; row += 1) {
      if (row === pivotRow || fractionIsZero(matrix[row][column])) continue;
      const factor = matrix[row][column];
      for (let j = column; j < columns; j += 1) {
        matrix[row][j] = fractionSubtract(
          matrix[row][j], fractionMultiply(factor, matrix[pivotRow][j]),
        );
      }
    }
    pivots.push(column);
    pivotRow += 1;
  }
  return { matrix, pivots, rank: pivots.length };
}

function rationalRows(integerRows) {
  return integerRows.map((row) => row.map((entry) => makeFraction(entry)));
}

function rankOfIntegerRows(rows, dimension) {
  if (rows.length === 0) return 0;
  return rref(rationalRows(rows), dimension).rank;
}

function nullspaceIntegerRows(rows, dimension) {
  if (dimension === 0) return [];
  const reduced = rref(rationalRows(rows), dimension);
  const pivotSet = new Set(reduced.pivots);
  const basis = [];
  for (let free = 0; free < dimension; free += 1) {
    if (pivotSet.has(free)) continue;
    const vector = Array.from({ length: dimension }, zeroFraction);
    vector[free] = makeFraction(1n);
    for (let row = 0; row < reduced.pivots.length; row += 1) {
      vector[reduced.pivots[row]] = fractionNegate(reduced.matrix[row][free]);
    }
    basis.push(vector);
  }
  return basis;
}

function solveUnique(integerRows, integerRightSides, dimension) {
  const augmented = integerRows.map((row, index) => (
    row.map((entry) => makeFraction(entry)).concat([makeFraction(integerRightSides[index])])
  ));
  const reduced = rref(augmented, dimension);
  for (const row of reduced.matrix) {
    let allZero = true;
    for (let column = 0; column < dimension; column += 1) {
      if (!fractionIsZero(row[column])) {
        allZero = false;
        break;
      }
    }
    if (allZero && !fractionIsZero(row[dimension])) return null;
  }
  if (reduced.rank !== dimension) return null;
  const answer = Array.from({ length: dimension }, zeroFraction);
  for (let row = 0; row < reduced.pivots.length; row += 1) {
    answer[reduced.pivots[row]] = reduced.matrix[row][dimension];
  }
  return answer;
}

function equalitiesConsistent(polyhedron) {
  if (polyhedron.equalities.length === 0) return true;
  const augmented = polyhedron.equalities.map((constraint) => (
    constraint.a.map((entry) => makeFraction(entry)).concat([makeFraction(constraint.b)])
  ));
  const reduced = rref(augmented, polyhedron.dimension);
  for (const row of reduced.matrix) {
    let allZero = true;
    for (let column = 0; column < polyhedron.dimension; column += 1) {
      if (!fractionIsZero(row[column])) {
        allZero = false;
        break;
      }
    }
    if (allZero && !fractionIsZero(row[polyhedron.dimension])) return false;
  }
  return true;
}

function pointSatisfiesContinuous(polyhedron, point) {
  for (const constraint of polyhedron.equalities) {
    if (fractionCompare(dotFraction(constraint.a, point), makeFraction(constraint.b)) !== 0) {
      return false;
    }
  }
  for (const constraint of polyhedron.inequalities) {
    if (fractionCompare(dotFraction(constraint.a, point), makeFraction(constraint.b)) > 0) {
      return false;
    }
  }
  return true;
}

function pointSatisfiesInteger(polyhedron, point, includeDisequalities = true) {
  for (const constraint of polyhedron.equalities) {
    if (dotBigInt(constraint.a, point) !== constraint.b) return false;
  }
  for (const constraint of polyhedron.inequalities) {
    if (dotBigInt(constraint.a, point) > constraint.b) return false;
  }
  if (includeDisequalities) {
    for (const constraint of polyhedron.disequalities) {
      if (dotBigInt(constraint.a, point) === constraint.b) return false;
    }
  }
  return true;
}

function raySatisfiesRecession(polyhedron, ray) {
  if (ray.every((entry) => entry === 0n)) return false;
  for (const constraint of polyhedron.equalities) {
    if (dotBigInt(constraint.a, ray) !== 0n) return false;
  }
  for (const constraint of polyhedron.inequalities) {
    if (dotBigInt(constraint.a, ray) > 0n) return false;
  }
  return true;
}

function forEachCombination(size, choose, state, callback) {
  if (choose < 0 || choose > size) return true;
  const current = [];
  function visit(start) {
    if (current.length === choose) {
      state.activeSets += 1;
      if (state.activeSets > state.options.maxActiveSets) {
        throw new PolyhedralLimitError(
          `active-set limit ${state.options.maxActiveSets} exceeded`,
        );
      }
      return callback(current.slice()) !== false;
    }
    const remaining = choose - current.length;
    for (let value = start; value <= size - remaining; value += 1) {
      current.push(value);
      if (!visit(value + 1)) return false;
      current.pop();
    }
    return true;
  }
  return visit(0);
}

function findRecessionDirection(polyhedron, state) {
  const n = polyhedron.dimension;
  if (n === 0) return null;
  const equalityRows = polyhedron.equalities.map((constraint) => constraint.a);
  const inequalityRows = polyhedron.inequalities.map((constraint) => constraint.a);
  const equalityRank = rankOfIntegerRows(equalityRows, n);
  const affineDimension = n - equalityRank;
  if (affineDimension === 0) return null;

  // First extract the lineality space E v=0, A v=0.  Such a direction can
  // be followed both forwards and backwards.
  const lineality = nullspaceIntegerRows(equalityRows.concat(inequalityRows), n);
  if (lineality.length !== 0) {
    const ray = integerVectorFromRational(lineality[0]);
    if (!raySatisfiesRecession(polyhedron, ray)
        || !raySatisfiesRecession(polyhedron, ray.map((entry) => -entry))) {
      throw new Error("internal lineality verification failed");
    }
    return { ray, kind: "lineality" };
  }

  // A nonzero pointed rational cone has an extreme ray.  In an affine
  // space of dimension d, an extreme ray contains d-1 independent active
  // facets.  Enumerating those active sets is exact and very effective for
  // the small annular-eye dimensions.
  const needed = affineDimension - 1;
  let answer = null;
  forEachCombination(inequalityRows.length, needed, state, (indices) => {
    const rows = equalityRows.concat(indices.map((index) => inequalityRows[index]));
    const nullspace = nullspaceIntegerRows(rows, n);
    if (nullspace.length !== 1) return true;
    const candidate = integerVectorFromRational(nullspace[0]);
    for (const sign of [1n, -1n]) {
      const ray = candidate.map((entry) => sign * entry);
      if (raySatisfiesRecession(polyhedron, ray)) {
        answer = { ray, kind: "pointed" };
        return false;
      }
    }
    return true;
  });
  return answer;
}

function enumerateVertices(polyhedron, state) {
  const n = polyhedron.dimension;
  const equalityRows = polyhedron.equalities.map((constraint) => constraint.a);
  const equalityRightSides = polyhedron.equalities.map((constraint) => constraint.b);
  const equalityRank = rankOfIntegerRows(equalityRows, n);
  const affineDimension = n - equalityRank;
  const vertices = [];
  const seen = new Set();
  const addVertex = (point) => {
    if (point === null || !pointSatisfiesContinuous(polyhedron, point)) return;
    const key = point.map((entry) => `${entry.n}/${entry.d}`).join(",");
    if (!seen.has(key)) {
      seen.add(key);
      vertices.push(point);
    }
  };
  if (affineDimension === 0) {
    addVertex(solveUnique(equalityRows, equalityRightSides, n));
    return vertices;
  }
  forEachCombination(polyhedron.inequalities.length, affineDimension, state, (indices) => {
    const rows = equalityRows.concat(indices.map((index) => polyhedron.inequalities[index].a));
    const rightSides = equalityRightSides.concat(
      indices.map((index) => polyhedron.inequalities[index].b),
    );
    addVertex(solveUnique(rows, rightSides, n));
    return true;
  });
  return vertices;
}

function integerBoundsFromVertices(vertices, dimension) {
  if (vertices.length === 0) return null;
  const lower = [];
  const upper = [];
  for (let variable = 0; variable < dimension; variable += 1) {
    let minimum = vertices[0][variable];
    let maximum = vertices[0][variable];
    for (let index = 1; index < vertices.length; index += 1) {
      if (fractionCompare(vertices[index][variable], minimum) < 0) {
        minimum = vertices[index][variable];
      }
      if (fractionCompare(vertices[index][variable], maximum) > 0) {
        maximum = vertices[index][variable];
      }
    }
    lower.push(fractionCeil(minimum));
    upper.push(fractionFloor(maximum));
    if (lower[variable] > upper[variable]) return null;
  }
  return { lower, upper };
}

function partialFeasible(polyhedron, assigned, bounds) {
  function intervalFor(coefficients) {
    let minimum = 0n;
    let maximum = 0n;
    for (let variable = 0; variable < coefficients.length; variable += 1) {
      const coefficient = coefficients[variable];
      if (assigned[variable] !== null) {
        minimum += coefficient * assigned[variable];
        maximum += coefficient * assigned[variable];
      } else if (coefficient >= 0n) {
        minimum += coefficient * bounds.lower[variable];
        maximum += coefficient * bounds.upper[variable];
      } else {
        minimum += coefficient * bounds.upper[variable];
        maximum += coefficient * bounds.lower[variable];
      }
    }
    return [minimum, maximum];
  }
  for (const constraint of polyhedron.inequalities) {
    const [minimum] = intervalFor(constraint.a);
    if (minimum > constraint.b) return false;
  }
  for (const constraint of polyhedron.equalities) {
    const [minimum, maximum] = intervalFor(constraint.a);
    if (constraint.b < minimum || constraint.b > maximum) return false;
    let fixed = 0n;
    let divisor = 0n;
    for (let variable = 0; variable < constraint.a.length; variable += 1) {
      if (assigned[variable] === null) divisor = gcdBigInt(divisor, constraint.a[variable]);
      else fixed += constraint.a[variable] * assigned[variable];
    }
    if (divisor === 0n) {
      if (fixed !== constraint.b) return false;
    } else if ((constraint.b - fixed) % divisor !== 0n) return false;
  }
  return true;
}

/*
 * Cheap exact presolve for the phase systems produced below.  Almost all
 * geometric constraints first give one-variable facing bounds; type-B then
 * couples at most two phase variables and balance contributes one equality.
 * Running a general active-set search before using those bounds was the main
 * cost of the no-K O(4) run.
 *
 * Every bound derived here is a necessary consequence of an original
 * inequality.  Therefore an inconsistency is an exact emptiness certificate,
 * and, when every variable is bounded, enumeration of the resulting integer
 * box is exhaustive.  Systems not settled by this presolve fall back to the
 * general recession/vertex solver below.
 */
function intervalPresolve(polyhedron) {
  const dimension = polyhedron.dimension;
  const lower = new Array(dimension).fill(null);
  const upper = new Array(dimension).fill(null);
  let impossible = false;
  let changed = false;

  function tightenLower(variable, value) {
    if (lower[variable] === null || value > lower[variable]) {
      lower[variable] = value;
      changed = true;
    }
    if (upper[variable] !== null && lower[variable] > upper[variable]) impossible = true;
  }

  function tightenUpper(variable, value) {
    if (upper[variable] === null || value < upper[variable]) {
      upper[variable] = value;
      changed = true;
    }
    if (lower[variable] !== null && lower[variable] > upper[variable]) impossible = true;
  }

  function propagateInequality(coefficients, rightSide) {
    let globalMinimum = 0n;
    let globalMinimumFinite = true;
    for (let variable = 0; variable < dimension; variable += 1) {
      const coefficient = coefficients[variable];
      if (coefficient === 0n) continue;
      const endpoint = coefficient > 0n ? lower[variable] : upper[variable];
      if (endpoint === null) {
        globalMinimumFinite = false;
        break;
      }
      globalMinimum += coefficient * endpoint;
    }
    if (globalMinimumFinite && globalMinimum > rightSide) {
      impossible = true;
      return;
    }

    for (let variable = 0; variable < dimension && !impossible; variable += 1) {
      const coefficient = coefficients[variable];
      if (coefficient === 0n) continue;
      let otherMinimum = 0n;
      let finite = true;
      for (let other = 0; other < dimension; other += 1) {
        if (other === variable) continue;
        const otherCoefficient = coefficients[other];
        if (otherCoefficient === 0n) continue;
        const endpoint = otherCoefficient > 0n ? lower[other] : upper[other];
        if (endpoint === null) {
          finite = false;
          break;
        }
        otherMinimum += otherCoefficient * endpoint;
      }
      if (!finite) continue;
      const residual = rightSide - otherMinimum;
      if (coefficient > 0n) {
        tightenUpper(variable, floorDiv(residual, coefficient));
      } else {
        tightenLower(variable, ceilDiv(-residual, -coefficient));
      }
    }
  }

  function equalityDivisibility(constraint) {
    let fixed = 0n;
    let divisor = 0n;
    for (let variable = 0; variable < dimension; variable += 1) {
      const coefficient = constraint.a[variable];
      if (lower[variable] !== null && upper[variable] !== null
          && lower[variable] === upper[variable]) {
        fixed += coefficient * lower[variable];
      } else {
        divisor = gcdBigInt(divisor, coefficient);
      }
    }
    if (divisor === 0n) {
      if (fixed !== constraint.b) impossible = true;
    } else if ((constraint.b - fixed) % divisor !== 0n) impossible = true;
  }

  // Propagation is monotone.  The cap only limits presolve work; stopping
  // early leaves sound (possibly loose) bounds and lets the exact fallback
  // handle anything not already decided.
  for (let round = 0; round < 1024 && !impossible; round += 1) {
    changed = false;
    for (const constraint of polyhedron.inequalities) {
      propagateInequality(constraint.a, constraint.b);
      if (impossible) break;
    }
    for (const constraint of polyhedron.equalities) {
      propagateInequality(constraint.a, constraint.b);
      if (impossible) break;
      propagateInequality(
        constraint.a.map((entry) => -entry), -constraint.b,
      );
      if (impossible) break;
      equalityDivisibility(constraint);
      if (impossible) break;
    }
    if (!changed) break;
  }
  return {
    impossible,
    bounded: !impossible && lower.every((entry) => entry !== null)
      && upper.every((entry) => entry !== null),
    bounds: { lower, upper },
  };
}

function enumerateIntegerPointsInBounds(
  polyhedron, state, bounds, stopAfter = null, includeDisequalities = true,
) {
  const order = Array.from({ length: polyhedron.dimension }, (_, index) => index)
    .sort((left, right) => {
      const leftWidth = bounds.upper[left] - bounds.lower[left];
      const rightWidth = bounds.upper[right] - bounds.lower[right];
      return leftWidth < rightWidth ? -1 : leftWidth > rightWidth ? 1 : left - right;
    });
  const assigned = new Array(polyhedron.dimension).fill(null);
  const points = [];
  let count = 0n;
  let stopped = false;
  function visit(depth) {
    if (stopped || !partialFeasible(polyhedron, assigned, bounds)) return;
    state.integerNodes += 1;
    if (state.integerNodes > state.options.maxIntegerNodes) {
      throw new PolyhedralLimitError(
        `integer-node limit ${state.options.maxIntegerNodes} exceeded`,
      );
    }
    if (depth === order.length) {
      const point = assigned.slice();
      if (!pointSatisfiesInteger(polyhedron, point, includeDisequalities)) return;
      count += 1n;
      if (points.length < state.options.maxCollectedPoints) points.push(point);
      if (stopAfter !== null && count >= stopAfter) stopped = true;
      return;
    }
    const variable = order[depth];
    for (let value = bounds.lower[variable]; value <= bounds.upper[variable]; value += 1n) {
      assigned[variable] = value;
      visit(depth + 1);
      if (stopped) break;
    }
    assigned[variable] = null;
  }
  visit(0);
  return { count, points };
}

function enumerateBoundedIntegerPoints(
  polyhedron, state, stopAfter = null, includeDisequalities = true,
) {
  const vertices = enumerateVertices(polyhedron, state);
  const bounds = integerBoundsFromVertices(vertices, polyhedron.dimension);
  if (bounds === null) return { count: 0n, points: [], vertices: vertices.length };
  const result = enumerateIntegerPointsInBounds(
    polyhedron, state, bounds, stopAfter, includeDisequalities,
  );
  return { ...result, vertices: vertices.length };
}

function clonePolyhedron(polyhedron) {
  return {
    dimension: polyhedron.dimension,
    impossible: polyhedron.impossible,
    inequalities: polyhedron.inequalities.map((constraint) => ({
      a: constraint.a.slice(), b: constraint.b, label: constraint.label,
    })),
    equalities: polyhedron.equalities.map((constraint) => ({
      a: constraint.a.slice(), b: constraint.b, label: constraint.label,
    })),
    disequalities: polyhedron.disequalities.map((constraint) => ({
      a: constraint.a.slice(), b: constraint.b, label: constraint.label,
    })),
  };
}

function appendIntegerInequality(polyhedron, coefficients, rightSide, label) {
  const next = clonePolyhedron(polyhedron);
  next.inequalities.push({ a: coefficients.slice(), b: rightSide, label });
  return next;
}

function findIntegerPointModuloRecession(polyhedron, state, depth = 0) {
  if (depth > polyhedron.dimension + 1) {
    throw new Error("recession fundamental-domain recursion did not decrease dimension");
  }
  if (!equalitiesConsistent(polyhedron)) return null;
  const recession = findRecessionDirection(polyhedron, state);
  if (recession === null) {
    // Fundamental-domain normalization is valid for the linear
    // equalities/inequalities.  A finite collection of disequality
    // hyperplanes is not translation invariant (and could delete the chosen
    // representative while leaving the rest of its ray intact), so find the
    // seed in the closed polyhedron and move it beyond all excluded values
    // only after returning to the original recession ray.
    const bounded = enumerateBoundedIntegerPoints(polyhedron, state, 1n, false);
    return bounded.count === 0n ? null : bounded.points[0];
  }
  const { ray, kind } = recession;
  if (kind === "lineality") {
    const coordinate = ray.findIndex((entry) => entry !== 0n);
    const modulus = absBigInt(ray[coordinate]);
    const lowerCoefficients = new Array(polyhedron.dimension).fill(0n);
    lowerCoefficients[coordinate] = -1n;
    let branch = appendIntegerInequality(
      polyhedron, lowerCoefficients, 0n, `lineality-fundamental-lower-${depth}`,
    );
    const upperCoefficients = new Array(polyhedron.dimension).fill(0n);
    upperCoefficients[coordinate] = 1n;
    branch = appendIntegerInequality(
      branch, upperCoefficients, modulus - 1n, `lineality-fundamental-upper-${depth}`,
    );
    return findIntegerPointModuloRecession(branch, state, depth + 1);
  }

  const blockers = [];
  for (let index = 0; index < polyhedron.inequalities.length; index += 1) {
    const constraint = polyhedron.inequalities[index];
    const drift = dotBigInt(constraint.a, ray);
    if (drift < 0n) blockers.push({ index, constraint, drift });
  }
  if (blockers.length === 0) {
    throw new Error("pointed recession ray has no strict blocking inequality");
  }
  // Every integer point can be translated backwards along ray until one of
  // these finitely many inequalities blocks the next step.  The resulting
  // point lies in one of the exact one-step slabs below.
  for (const blocker of blockers) {
    const lower = blocker.constraint.b + blocker.drift + 1n;
    const coefficients = blocker.constraint.a.map((entry) => -entry);
    const branch = appendIntegerInequality(
      polyhedron,
      coefficients,
      -lower,
      `pointed-fundamental-${depth}-${blocker.index}`,
    );
    const point = findIntegerPointModuloRecession(branch, state, depth + 1);
    if (point !== null) return point;
  }
  return null;
}

function movePoint(point, ray, amount) {
  return point.map((entry, index) => entry + amount * ray[index]);
}

function avoidFutureDisequalities(polyhedron, point, ray) {
  let lastForbidden = -1n;
  for (const constraint of polyhedron.disequalities) {
    const slope = dotBigInt(constraint.a, ray);
    const offset = constraint.b - dotBigInt(constraint.a, point);
    if (slope === 0n) {
      if (offset === 0n) throw new Error("seed violates a phase disequality");
      continue;
    }
    if (offset % slope !== 0n) continue;
    const parameter = offset / slope;
    if (parameter >= 0n && parameter > lastForbidden) lastForbidden = parameter;
  }
  const tailStart = lastForbidden + 1n;
  return { point: movePoint(point, ray, tailStart), tailStart };
}

function findSeedForInfiniteRay(polyhedron, ray, state) {
  // Hyperplanes transverse to the ray are met at most once and can be
  // skipped by moving far enough.  Hyperplanes parallel to the ray must be
  // avoided by the seed itself.  Turn each such disequality into its exact
  // two half-space branches before fundamentalizing the recession cone.
  const invariant = polyhedron.disequalities.filter(
    (constraint) => dotBigInt(constraint.a, ray) === 0n,
  );
  const base = clonePolyhedron(polyhedron);
  base.disequalities = [];
  function visit(index, current) {
    if (index === invariant.length) {
      return findIntegerPointModuloRecession(current, state);
    }
    const excluded = invariant[index];
    const below = appendIntegerInequality(
      current,
      excluded.a,
      excluded.b - 1n,
      `invariant-disequality-below-${index}`,
    );
    const belowPoint = visit(index + 1, below);
    if (belowPoint !== null) return belowPoint;
    const above = appendIntegerInequality(
      current,
      excluded.a.map((entry) => -entry),
      -excluded.b - 1n,
      `invariant-disequality-above-${index}`,
    );
    return visit(index + 1, above);
  }
  return visit(0, base);
}

function defaultOptions(options = {}) {
  return {
    maxActiveSets: options.maxActiveSets === undefined ? 2_000_000 : options.maxActiveSets,
    maxIntegerNodes: options.maxIntegerNodes === undefined ? 5_000_000 : options.maxIntegerNodes,
    maxCollectedPoints: options.maxCollectedPoints === undefined
      ? 100_000 : options.maxCollectedPoints,
  };
}

function solveNormalizedPolyhedron(polyhedron, options = {}) {
  const state = {
    options: defaultOptions(options),
    activeSets: 0,
    integerNodes: 0,
  };
  try {
    if (polyhedron.impossible || !equalitiesConsistent(polyhedron)) {
      return {
        status: "EMPTY",
        count: 0n,
        activeSets: state.activeSets,
        integerNodes: state.integerNodes,
      };
    }
    const presolved = intervalPresolve(polyhedron);
    if (presolved.impossible) {
      return {
        status: "EMPTY",
        count: 0n,
        presolve: "interval-infeasible",
        activeSets: state.activeSets,
        integerNodes: state.integerNodes,
      };
    }
    if (presolved.bounded) {
      const finite = enumerateIntegerPointsInBounds(
        polyhedron, state, presolved.bounds, null, true,
      );
      return {
        status: finite.count === 0n ? "EMPTY" : "FINITE",
        count: finite.count,
        points: finite.points,
        collectedAll: finite.count <= BigInt(state.options.maxCollectedPoints),
        vertices: null,
        presolve: "interval-bounded",
        activeSets: state.activeSets,
        integerNodes: state.integerNodes,
      };
    }
    const recession = findRecessionDirection(polyhedron, state);
    if (recession === null) {
      const finite = enumerateBoundedIntegerPoints(polyhedron, state, null);
      return {
        status: finite.count === 0n ? "EMPTY" : "FINITE",
        count: finite.count,
        points: finite.points,
        collectedAll: finite.count <= BigInt(state.options.maxCollectedPoints),
        vertices: finite.vertices,
        activeSets: state.activeSets,
        integerNodes: state.integerNodes,
      };
    }

    const seed = findSeedForInfiniteRay(polyhedron, recession.ray, state);
    if (seed === null) {
      return {
        status: "EMPTY",
        count: 0n,
        unboundedRelaxation: true,
        activeSets: state.activeSets,
        integerNodes: state.integerNodes,
      };
    }
    const tail = avoidFutureDisequalities(polyhedron, seed, recession.ray);
    if (!pointSatisfiesInteger(polyhedron, tail.point, true)
        || !raySatisfiesRecession(polyhedron, recession.ray)) {
      throw new Error("internal infinite-certificate verification failed");
    }
    return {
      status: "INFINITE",
      point: tail.point,
      ray: recession.ray,
      rayKind: recession.kind,
      tailStart: tail.tailStart,
      activeSets: state.activeSets,
      integerNodes: state.integerNodes,
    };
  } catch (error) {
    if (!(error instanceof PolyhedralLimitError)) throw error;
    return {
      status: "UNKNOWN",
      reason: error.message,
      activeSets: state.activeSets,
      integerNodes: state.integerNodes,
    };
  }
}

function normalizeLinearConstraint(dimension, coefficients, operator, rightSide, label = null) {
  if (!Array.isArray(coefficients) || coefficients.length !== dimension) {
    throw new Error(`constraint needs ${dimension} coefficients`);
  }
  const rationalCoefficients = coefficients.map(asFraction);
  const rationalRightSide = asFraction(rightSide);
  let denominator = rationalRightSide.d;
  for (const coefficient of rationalCoefficients) {
    denominator = lcmBigInt(denominator, coefficient.d);
  }
  let a = rationalCoefficients.map((coefficient) => (
    coefficient.n * (denominator / coefficient.d)
  ));
  let b = rationalRightSide.n * (denominator / rationalRightSide.d);
  let relation = operator;
  if (relation === "<") {
    relation = "<=";
    b -= 1n;
  } else if (relation === ">") {
    relation = ">=";
    b += 1n;
  }
  if (relation === ">=") {
    relation = "<=";
    a = a.map((entry) => -entry);
    b = -b;
  }
  if (!new Set(["<=", "=", "!="]).has(relation)) {
    throw new Error(`unknown linear relation ${operator}`);
  }

  let divisor = 0n;
  for (const entry of a) divisor = gcdBigInt(divisor, entry);
  if (divisor === 0n) return { a, b, relation, label, constant: true };
  if (relation === "<=") {
    b = floorDiv(b, divisor);
    a = a.map((entry) => entry / divisor);
  } else if (b % divisor === 0n) {
    b /= divisor;
    a = a.map((entry) => entry / divisor);
  } else if (relation === "=") {
    return { a, b, relation, label, impossible: true };
  } else {
    return { a, b, relation, label, tautology: true };
  }
  if (relation !== "<=") {
    const first = a.find((entry) => entry !== 0n);
    if (first !== undefined && first < 0n) {
      a = a.map((entry) => -entry);
      b = -b;
    }
  }
  return { a, b, relation, label };
}

class IntegerPolyhedron {
  constructor(dimension) {
    if (!Number.isInteger(dimension) || dimension < 0) {
      throw new Error(`bad integer-polyhedron dimension ${dimension}`);
    }
    this.dimension = dimension;
    this.inequalities = [];
    this.equalities = [];
    this.disequalities = [];
    this.impossible = false;
    this._inequalityIndex = new Map();
    this._equalityKeys = new Set();
    this._disequalityKeys = new Set();
  }

  add(coefficients, operator, rightSide, label = null) {
    const constraint = normalizeLinearConstraint(
      this.dimension, coefficients, operator, rightSide, label,
    );
    if (constraint.impossible) {
      this.impossible = true;
      return this;
    }
    if (constraint.tautology) return this;
    if (constraint.constant) {
      const holds = constraint.relation === "<="
        ? 0n <= constraint.b
        : constraint.relation === "="
          ? constraint.b === 0n
          : constraint.b !== 0n;
      if (!holds) this.impossible = true;
      return this;
    }
    const stored = { a: constraint.a, b: constraint.b, label: constraint.label };
    const coefficientKey = stored.a.join(",");
    if (constraint.relation === "<=") {
      // Parallel inequalities with the same normal contribute only their
      // tightest right side.  Facing sheets and the two type-B probes often
      // generate exact duplicates; retaining all of them makes active-set
      // enumeration needlessly combinatorial.
      const existing = this._inequalityIndex.get(coefficientKey);
      if (existing === undefined) {
        this._inequalityIndex.set(coefficientKey, this.inequalities.length);
        this.inequalities.push(stored);
      } else if (stored.b < this.inequalities[existing].b) {
        this.inequalities[existing] = stored;
      }
    } else if (constraint.relation === "=") {
      const key = `${coefficientKey}=${stored.b}`;
      if (!this._equalityKeys.has(key)) {
        this._equalityKeys.add(key);
        this.equalities.push(stored);
      }
    } else {
      const key = `${coefficientKey}!=${stored.b}`;
      if (!this._disequalityKeys.has(key)) {
        this._disequalityKeys.add(key);
        this.disequalities.push(stored);
      }
    }
    return this;
  }

  solve(options = {}) {
    return solveNormalizedPolyhedron(this, options);
  }

  verifyPoint(point) {
    if (!Array.isArray(point) || point.length !== this.dimension) return false;
    const integerPoint = point.map((entry) => BigInt(entry));
    return pointSatisfiesInteger(this, integerPoint, true);
  }

  verifyRay(ray) {
    if (!Array.isArray(ray) || ray.length !== this.dimension) return false;
    const integerRay = ray.map((entry) => BigInt(entry));
    return raySatisfiesRecession(this, integerRay);
  }
}

function phaseResidueData(hx, hy, tauX, tauY) {
  hx = BigInt(hx);
  hy = BigInt(hy);
  tauX = BigInt(tauX);
  tauY = BigInt(tauY);
  const tauH = tauX * hx + tauY * hy;
  if (tauH === 0n) throw new Error("sweep covector is tangent to annular homology");
  const q = absBigInt(tauH);

  // Extended Euclid, with det(h,p)=hx*p_y-hy*p_x=1.
  function extendedGcd(a, b) {
    const signA = a < 0n ? -1n : 1n;
    const signB = b < 0n ? -1n : 1n;
    let oldR = absBigInt(a);
    let r = absBigInt(b);
    let oldS = 1n;
    let s = 0n;
    let oldT = 0n;
    let t = 1n;
    while (r !== 0n) {
      const quotient = oldR / r;
      [oldR, r] = [r, oldR - quotient * r];
      [oldS, s] = [s, oldS - quotient * s];
      [oldT, t] = [t, oldT - quotient * t];
    }
    return [oldR, oldS * signA, oldT * signB];
  }
  const [divisor, xCoefficient, yCoefficient] = extendedGcd(hx, hy);
  if (divisor !== 1n) throw new Error("annular homology must be primitive");
  const px = -yCoefficient;
  const py = xCoefficient;
  const tauP = tauX * px + tauY * py;
  const signTauH = tauH < 0n ? -1n : 1n;
  const wx = q * px - signTauH * tauP * hx;
  const wy = q * py - signTauH * tauP * hy;
  if (tauX * wx + tauY * wy !== 0n || hx * wy - hy * wx !== q) {
    throw new Error("internal phase-period identity failed");
  }
  return {
    q,
    phaseStep: [px, py],
    fixedSlicePeriod: [wx, wy],
    tauH,
  };
}

function fractionScale(value, scalar) {
  return makeFraction(value.n * BigInt(scalar), value.d);
}

function exactLinearForm(covector, x, y) {
  return fractionAdd(
    fractionScale(x, covector[0]),
    fractionScale(y, covector[1]),
  );
}

function exactPointFromIntegers(x, y) {
  return { x: makeFraction(x), y: makeFraction(y) };
}

function exactPointAddVector(point, dx, dy) {
  return {
    x: fractionAdd(point.x, makeFraction(dx)),
    y: fractionAdd(point.y, makeFraction(dy)),
  };
}

function exactPointDifference(right, left) {
  return {
    x: fractionSubtract(right.x, left.x),
    y: fractionSubtract(right.y, left.y),
  };
}

function exactCross(vx, vy, vector) {
  return fractionSubtract(
    fractionScale(vector.y, vx),
    fractionScale(vector.x, vy),
  );
}

function exactSign(value) {
  return value.n < 0n ? -1 : value.n > 0n ? 1 : 0;
}

function exactPointEqual(left, right) {
  return fractionCompare(left.x, right.x) === 0
    && fractionCompare(left.y, right.y) === 0;
}

function exactPointOnSegment(point, left, right) {
  const segment = exactPointDifference(right, left);
  const offset = exactPointDifference(point, left);
  const cross = fractionSubtract(
    fractionMultiply(segment.x, offset.y),
    fractionMultiply(segment.y, offset.x),
  );
  if (exactSign(cross) !== 0) return false;
  const between = (value, first, second) => {
    const low = fractionCompare(first, second) <= 0 ? first : second;
    const high = fractionCompare(first, second) <= 0 ? second : first;
    return fractionCompare(low, value) <= 0 && fractionCompare(value, high) <= 0;
  };
  return between(point.x, left.x, right.x) && between(point.y, left.y, right.y);
}

function pointInPolygonExact(polygon, point) {
  let inside = false;
  for (let index = 0, previous = polygon.length - 1;
    index < polygon.length; previous = index, index += 1) {
    const left = exactPointFromIntegers(polygon[previous][0], polygon[previous][1]);
    const right = exactPointFromIntegers(polygon[index][0], polygon[index][1]);
    if (exactPointOnSegment(point, left, right)) {
      return { inside: false, onBoundary: true };
    }
    const leftAbove = fractionCompare(left.y, point.y) > 0;
    const rightAbove = fractionCompare(right.y, point.y) > 0;
    if (leftAbove === rightAbove) continue;
    const alpha = fractionDivide(
      fractionSubtract(point.y, left.y),
      fractionSubtract(right.y, left.y),
    );
    const hitX = fractionAdd(
      left.x,
      fractionMultiply(alpha, fractionSubtract(right.x, left.x)),
    );
    if (fractionCompare(point.x, hitX) < 0) inside = !inside;
  }
  return { inside, onBoundary: false };
}

/*
 * Count disk-eye translates containing the chosen base point without the
 * floating point ray cast used by the legacy audit.  A base point on a
 * boundary is nongeneric input, so it is rejected instead of being assigned
 * to one side by a tolerance.
 */
function translatedDiskMultiplicityExact(arr, circuit, basePoint) {
  if (circuit.kind !== "disk") throw new Error("disk multiplicity needs a disk circuit");
  if (!Array.isArray(circuit.polygon) || circuit.polygon.length < 3) {
    throw new Error("disk circuit has no polygon");
  }
  const base = normalizeExactPoint(
    basePoint,
    [makeFraction(37n, 100n), makeFraction(73n, 100n)],
  );
  const period = BigInt(arr.L);
  if (period <= 0n) throw new Error("arrangement period must be positive");
  const xs = circuit.polygon.map((point) => BigInt(point[0]));
  const ys = circuit.polygon.map((point) => BigInt(point[1]));
  const xmin = xs.reduce((left, right) => (left < right ? left : right));
  const xmax = xs.reduce((left, right) => (left > right ? left : right));
  const ymin = ys.reduce((left, right) => (left < right ? left : right));
  const ymax = ys.reduce((left, right) => (left > right ? left : right));
  const translationBounds = (minimum, maximum, coordinate) => ({
    lower: fractionCeil(fractionDivide(
      fractionSubtract(makeFraction(minimum), coordinate), makeFraction(period),
    )),
    upper: fractionFloor(fractionDivide(
      fractionSubtract(makeFraction(maximum), coordinate), makeFraction(period),
    )),
  });
  const xBounds = translationBounds(xmin, xmax, base[0]);
  const yBounds = translationBounds(ymin, ymax, base[1]);
  let count = 0n;
  for (let xShift = xBounds.lower; xShift <= xBounds.upper; xShift += 1n) {
    for (let yShift = yBounds.lower; yShift <= yBounds.upper; yShift += 1n) {
      const point = {
        x: fractionAdd(base[0], makeFraction(xShift * period)),
        y: fractionAdd(base[1], makeFraction(yShift * period)),
      };
      const location = pointInPolygonExact(circuit.polygon, point);
      if (location.onBoundary) {
        throw new Error("the balancing base point lies on a disk-eye boundary");
      }
      if (location.inside) count += 1n;
    }
  }
  return count;
}

function normalizeCovector(covector, fallback) {
  const value = covector === undefined ? fallback : covector;
  if (!Array.isArray(value) || value.length !== 2) {
    throw new Error("a planar covector must have two entries");
  }
  return value.map((entry) => BigInt(entry));
}

function normalizeExactPoint(point, fallback) {
  const value = point === undefined ? fallback : point;
  if (!Array.isArray(value) || value.length !== 2) {
    throw new Error("an exact planar point must have two entries");
  }
  return value.map(asFraction);
}

/*
 * Hit one specified lift of a tau-monotone essential circuit with the exact
 * slice tau(x)=level.  Longitudinal translates are solved by floor/ceil
 * formulas, so this routine has no search window.
 */
function liftedLineHitExact(arr, circuit, shiftX, shiftY, level, tau) {
  const period = BigInt(arr.L);
  shiftX = BigInt(shiftX);
  shiftY = BigInt(shiftY);
  const rawPeriodX = period * BigInt(circuit.rawHx);
  const rawPeriodY = period * BigInt(circuit.rawHy);
  const rawPeriodTau = tau[0] * rawPeriodX + tau[1] * rawPeriodY;
  if (rawPeriodTau === 0n) return null;
  const reverse = rawPeriodTau < 0n;
  const periodX = reverse ? -rawPeriodX : rawPeriodX;
  const periodY = reverse ? -rawPeriodY : rawPeriodY;
  const periodTau = reverse ? -rawPeriodTau : rawPeriodTau;
  const hits = [];
  for (const rawSegment of circuit.segments) {
    const segment = rawSegment.map((entry) => BigInt(entry));
    const [rawAx, rawAy, rawBx, rawBy] = segment;
    const firstX = (reverse ? rawBx : rawAx) + shiftX * period;
    const firstY = (reverse ? rawBy : rawAy) + shiftY * period;
    const secondX = (reverse ? rawAx : rawBx) + shiftX * period;
    const secondY = (reverse ? rawAy : rawBy) + shiftY * period;
    const ta = tau[0] * firstX + tau[1] * firstY;
    const tb = tau[0] * secondX + tau[1] * secondY;
    if (tb <= ta) return null;
    const lower = fractionDivide(
      fractionSubtract(level, makeFraction(tb)),
      makeFraction(periodTau),
    );
    const upper = fractionDivide(
      fractionSubtract(level, makeFraction(ta)),
      makeFraction(periodTau),
    );
    const firstLongitudinal = fractionFloor(lower) + 1n;
    // Use the half-open segment convention [first, second).  With the old
    // strict-open convention, a slice through a shared circuit vertex was
    // counted by neither adjacent segment and a genuinely tau-monotone
    // circuit was misreported as `annular-phase-nonmonotone`.  The lower
    // endpoint is now included and the upper endpoint remains excluded, so
    // every such vertex contributes exactly one hit.
    const lastLongitudinal = fractionFloor(upper);
    for (let longitudinal = firstLongitudinal;
      longitudinal <= lastLongitudinal; longitudinal += 1n) {
      const shiftedTa = ta + longitudinal * periodTau;
      const alpha = fractionDivide(
        fractionSubtract(level, makeFraction(shiftedTa)),
        makeFraction(tb - ta),
      );
      const baseX = firstX + longitudinal * periodX;
      const baseY = firstY + longitudinal * periodY;
      hits.push({
        x: fractionAdd(makeFraction(baseX), fractionScale(alpha, secondX - firstX)),
        y: fractionAdd(makeFraction(baseY), fractionScale(alpha, secondY - firstY)),
        // Preserve the original trace orientation for the coorientation test.
        vx: rawBx - rawAx,
        vy: rawBy - rawAy,
      });
    }
  }
  return hits.length === 1 ? hits[0] : null;
}

function finiteSliceHitsExact(circuit, coordinateShiftX, coordinateShiftY, level, tau) {
  coordinateShiftX = BigInt(coordinateShiftX);
  coordinateShiftY = BigInt(coordinateShiftY);
  const hits = [];
  for (const rawSegment of circuit.segments) {
    const [rawAx, rawAy, rawBx, rawBy] = rawSegment.map((entry) => BigInt(entry));
    const ax = rawAx + coordinateShiftX;
    const ay = rawAy + coordinateShiftY;
    const bx = rawBx + coordinateShiftX;
    const by = rawBy + coordinateShiftY;
    const ta = tau[0] * ax + tau[1] * ay;
    const tb = tau[0] * bx + tau[1] * by;
    if (ta === tb) continue;
    const fromStart = fractionDivide(
      fractionSubtract(level, makeFraction(ta)),
      makeFraction(tb - ta),
    );
    // Follow the trace orientation and assign a shared endpoint to the
    // segment that starts there: alpha in [0,1).  This is the finite-circuit
    // counterpart of liftedLineHitExact's half-open convention.
    if (fractionCompare(fromStart, makeFraction(0n)) < 0
        || fractionCompare(fromStart, makeFraction(1n)) >= 0) continue;
    const alpha = fromStart;
    hits.push({
      x: fractionAdd(makeFraction(ax), fractionScale(alpha, bx - ax)),
      y: fractionAdd(makeFraction(ay), fractionScale(alpha, by - ay)),
      vx: rawBx - rawAx,
      vy: rawBy - rawAy,
    });
  }
  return hits;
}

function exactDeckCoordinate(arr, visit) {
  const vertex = arr.vertices[visit[0]];
  const numeratorX = BigInt(visit[3]) - BigInt(vertex.x);
  const numeratorY = BigInt(visit[4]) - BigInt(vertex.y);
  const period = BigInt(arr.L);
  if (numeratorX % period !== 0n || numeratorY % period !== 0n) {
    throw new Error("circuit visit has a nonintegral deck coordinate");
  }
  return [numeratorX / period, numeratorY / period];
}

function forbiddenAnnularPhasesExact(arr, circuits, pair) {
  const [leftId, rightId] = pair;
  const left = circuits[leftId];
  const right = circuits[rightId];
  const rightByVertex = new Map();
  for (const visit of right.visits) {
    if (!rightByVertex.has(visit[0])) rightByVertex.set(visit[0], []);
    rightByVertex.get(visit[0]).push(visit);
  }
  const forbidden = new Set();
  for (const leftVisit of left.visits) {
    for (const rightVisit of rightByVertex.get(leftVisit[0]) || []) {
      const [leftX, leftY] = exactDeckCoordinate(arr, leftVisit);
      const [rightX, rightY] = exactDeckCoordinate(arr, rightVisit);
      const dx = leftX - rightX;
      const dy = leftY - rightY;
      const phase = BigInt(left.hx) * dy - BigInt(left.hy) * dx;
      forbidden.add(phase.toString());
    }
  }
  return [...forbidden].map((entry) => BigInt(entry));
}

function constraintInVariable(dimension, variable, coefficient, operator, rhs, label) {
  const coefficients = new Array(dimension).fill(0n);
  coefficients[variable] = coefficient;
  return { coefficients, operator, rhs, label };
}

function installCompiledConstraints(polyhedron, constraints) {
  for (const constraint of constraints) {
    polyhedron.add(
      constraint.coefficients,
      constraint.operator,
      constraint.rhs,
      constraint.label,
    );
  }
  return polyhedron;
}

function annularMultiplicityAtPhaseExact(
  arr, circuits, pair, phase, basePoint, tau, sigma, residueData,
) {
  const [leftId, rightId] = pair;
  const left = circuits[leftId];
  const right = circuits[rightId];
  const [px, py] = residueData.phaseStep;
  const q = residueData.q;
  const period = BigInt(arr.L);
  const level = exactLinearForm(tau, basePoint[0], basePoint[1]);
  const baseSigma = exactLinearForm(sigma, basePoint[0], basePoint[1]);
  let multiplicity = 0n;
  const orderSlopes = [];
  for (let residue = 0n; residue < q; residue += 1n) {
    const leftShiftX = residue * px;
    const leftShiftY = residue * py;
    const rightShiftX = leftShiftX + phase * px;
    const rightShiftY = leftShiftY + phase * py;
    const nextShiftX = (residue + q) * px;
    const nextShiftY = (residue + q) * py;
    const leftHit = liftedLineHitExact(
      arr, left, leftShiftX, leftShiftY, level, tau,
    );
    const rightHit = liftedLineHitExact(
      arr, right, rightShiftX, rightShiftY, level, tau,
    );
    const nextLeftHit = liftedLineHitExact(
      arr, left, nextShiftX, nextShiftY, level, tau,
    );
    if (leftHit === null || rightHit === null || nextLeftHit === null) {
      return { ok: false, stage: "annular-phase-nonmonotone" };
    }
    const leftSigma = exactLinearForm(sigma, leftHit.x, leftHit.y);
    const rightSigma = exactLinearForm(sigma, rightHit.x, rightHit.y);
    const nextLeftSigma = exactLinearForm(sigma, nextLeftHit.x, nextLeftHit.y);
    const slicePeriod = fractionSubtract(nextLeftSigma, leftSigma);
    if (exactSign(slicePeriod) === 0) {
      return { ok: false, stage: "annular-phase-bad-transverse-coordinate" };
    }
    const difference = fractionSubtract(rightSigma, leftSigma);
    if (exactSign(difference) === 0) {
      return { ok: false, stage: "annular-phase-boundaries-meet" };
    }
    orderSlopes.push(exactSign(difference) * exactSign(slicePeriod));
    const low = fractionCompare(leftSigma, rightSigma) < 0 ? leftSigma : rightSigma;
    const high = fractionCompare(leftSigma, rightSigma) < 0 ? rightSigma : leftSigma;
    const magnitude = exactSign(slicePeriod) > 0
      ? slicePeriod : fractionNegate(slicePeriod);
    const lower = fractionDivide(fractionSubtract(baseSigma, high), magnitude);
    const upper = fractionDivide(fractionSubtract(baseSigma, low), magnitude);
    const count = fractionCeil(upper) - fractionFloor(lower) - 1n;
    if (count > 0n) multiplicity += count;
  }
  if (orderSlopes.some((entry) => entry !== orderSlopes[0])) {
    return { ok: false, stage: "annular-phase-inconsistent-sheet-order" };
  }
  return {
    ok: true,
    multiplicity,
    zSlope: BigInt(orderSlopes[0]) * q,
  };
}

/*
 * Compile one annular eye and one residue class k=r+qz into exact linear
 * constraints.  The returned constraints are ready to add to an
 * IntegerPolyhedron of `dimension` phase variables.
 */
function compileAnnularEyeResidue(options) {
  const {
    arr,
    circuits,
    pair,
    color,
    variable = 0,
    dimension = 1,
  } = options;
  if (color !== 1 && color !== -1) throw new Error("annular color must be +1 or -1");
  if (!Number.isInteger(variable) || variable < 0 || variable >= dimension) {
    throw new Error("bad annular phase variable index");
  }
  const tau = normalizeCovector(options.tau, [2n, 1n]);
  const sigma = normalizeCovector(options.sigma, [1n, 0n]);
  const basePoint = normalizeExactPoint(
    options.basePoint,
    [makeFraction(37n, 100n), makeFraction(73n, 100n)],
  );
  const [leftId, rightId] = pair;
  const left = circuits[leftId];
  const right = circuits[rightId];
  if (left.kind !== "essential" || right.kind !== "essential"
      || left.hx !== right.hx || left.hy !== right.hy) {
    throw new Error("annular pair must have the same essential primitive homology");
  }
  const residueData = phaseResidueData(left.hx, left.hy, tau[0], tau[1]);
  const q = residueData.q;
  const residue = BigInt(options.residue);
  if (residue < 0n || residue >= q) {
    throw new Error(`phase residue ${residue} is outside [0,${q})`);
  }
  const [px, py] = residueData.phaseStep;
  const [wx, wy] = residueData.fixedSlicePeriod;
  const period = BigInt(arr.L);
  const level = exactLinearForm(tau, basePoint[0], basePoint[1]);
  const localConstraints = [];
  for (let sheet = 0n; sheet < q; sheet += 1n) {
    const leftShiftX = sheet * px;
    const leftShiftY = sheet * py;
    const rightShiftX = (sheet + residue) * px;
    const rightShiftY = (sheet + residue) * py;
    const leftHit = liftedLineHitExact(
      arr, left, leftShiftX, leftShiftY, level, tau,
    );
    const rightHit = liftedLineHitExact(
      arr, right, rightShiftX, rightShiftY, level, tau,
    );
    if (leftHit === null || rightHit === null) {
      return {
        ok: false,
        stage: "annular-phase-nonmonotone",
        pair: pair.slice(),
        residue,
        q,
      };
    }
    const difference = exactPointDifference(rightHit, leftHit);
    const differenceSlope = exactPointFromIntegers(period * wx, period * wy);
    const leftConstant = exactCross(leftHit.vx, leftHit.vy, difference);
    const leftCoefficient = exactCross(leftHit.vx, leftHit.vy, differenceSlope);
    const rightConstant = fractionNegate(
      exactCross(rightHit.vx, rightHit.vy, difference),
    );
    const rightCoefficient = fractionNegate(
      exactCross(rightHit.vx, rightHit.vy, differenceSlope),
    );
    const leftInside = BigInt(color * left.traceSide);
    const rightInside = BigInt(color * right.traceSide);
    localConstraints.push({
      coefficient: fractionScale(leftCoefficient, leftInside),
      operator: ">",
      rhs: fractionNegate(fractionScale(leftConstant, leftInside)),
      label: `annular-facing-left-sheet-${sheet}`,
    });
    localConstraints.push({
      coefficient: fractionScale(rightCoefficient, rightInside),
      operator: ">",
      rhs: fractionNegate(fractionScale(rightConstant, rightInside)),
      label: `annular-facing-right-sheet-${sheet}`,
    });
  }

  const forbiddenPhases = forbiddenAnnularPhasesExact(arr, circuits, pair);
  const forbiddenZ = forbiddenPhases
    .filter((phase) => (phase - residue) % q === 0n)
    .map((phase) => (phase - residue) / q);
  const local = new IntegerPolyhedron(1);
  for (const constraint of localConstraints) {
    local.add([constraint.coefficient], constraint.operator, constraint.rhs, constraint.label);
  }
  for (const excluded of forbiddenZ) local.add([1n], "!=", excluded, "boundary-meeting");
  const localSolverOptions = { ...(options.solverOptions || {}) };
  if (localSolverOptions.maxCollectedPoints !== undefined
      && localSolverOptions.maxCollectedPoints < 1) {
    localSolverOptions.maxCollectedPoints = 1;
  }
  const localResult = local.solve(localSolverOptions);
  if (localResult.status === "UNKNOWN") {
    return {
      ok: false,
      stage: "annular-phase-polyhedron-unknown",
      reason: localResult.reason,
      pair: pair.slice(),
      residue,
      q,
    };
  }
  if (localResult.status === "EMPTY") {
    return {
      ok: false,
      stage: "annular-phase-incompatible-sides",
      pair: pair.slice(),
      residue,
      q,
    };
  }
  if (localResult.status === "FINITE"
      && (!Array.isArray(localResult.points) || localResult.points.length === 0)) {
    return {
      ok: false,
      stage: "annular-phase-polyhedron-unknown",
      reason: "finite local phase cell did not retain a sample point",
      pair: pair.slice(),
      residue,
      q,
    };
  }
  const sampleZ = localResult.status === "INFINITE"
    ? localResult.point[0] : localResult.points[0][0];
  const samplePhase = residue + q * sampleZ;
  const sampleMultiplicity = annularMultiplicityAtPhaseExact(
    arr, circuits, pair, samplePhase, basePoint, tau, sigma, residueData,
  );
  if (!sampleMultiplicity.ok) {
    return {
      ok: false,
      stage: sampleMultiplicity.stage,
      pair: pair.slice(),
      residue,
      q,
    };
  }
  const multiplicityCoefficient = sampleMultiplicity.zSlope;
  const multiplicityConstant = sampleMultiplicity.multiplicity
    - multiplicityCoefficient * sampleZ;
  const constraints = localConstraints.map((constraint) => constraintInVariable(
    dimension,
    variable,
    constraint.coefficient,
    constraint.operator,
    constraint.rhs,
    constraint.label,
  ));
  for (const excluded of forbiddenZ) {
    constraints.push(constraintInVariable(
      dimension, variable, 1n, "!=", excluded, "annular-boundaries-meet",
    ));
  }
  return {
    ok: true,
    stage: "annular-phase-affine",
    pair: pair.slice(),
    left: leftId,
    right: rightId,
    color,
    variable,
    dimension,
    residue,
    q,
    phaseStep: residueData.phaseStep.slice(),
    fixedSlicePeriod: residueData.fixedSlicePeriod.slice(),
    constraints,
    localConstraints,
    forbiddenPhases,
    forbiddenZ,
    sampleZ,
    samplePhase,
    multiplicity: {
      coefficient: multiplicityCoefficient,
      constant: multiplicityConstant,
    },
    tau,
    sigma,
    basePoint,
  };
}

function matchingVisitExact(circuit, vertex, pair) {
  const key = pair.slice().sort((left, right) => left - right);
  const visits = circuit.visits.filter((visit) => {
    if (visit[0] !== vertex.index) return false;
    const local = [visit[1], visit[2]].sort((left, right) => left - right);
    return local[0] === key[0] && local[1] === key[1];
  });
  return visits.length === 1 ? visits[0] : null;
}

function lookupAnnularDescriptor(descriptors, eyeId) {
  if (descriptors instanceof Map) return descriptors.get(eyeId) || null;
  if (Array.isArray(descriptors)) return descriptors[eyeId] || null;
  if (descriptors !== null && typeof descriptors === "object") {
    return descriptors[eyeId] || null;
  }
  return null;
}

/*
 * Compile the normal same-order condition at one already-classified type-B
 * switch.  This function intentionally does NOT impose convex/outwardness:
 * the caller classifies B by the article's rule (exactly one incident corner
 * convex), so the inward type-B pictures remain allowed.
 */
function compileTypeBSameOrder(options) {
  const {
    arr,
    circuits,
    vertex,
    circuitIds,
    eyeIds,
    eyeCircuits,
    annularByEye,
    dimension,
  } = options;
  const tau = normalizeCovector(options.tau, [2n, 1n]);
  const sigma = normalizeCovector(options.sigma, [1n, 0n]);
  const epsilon = asFraction(
    options.epsilon === undefined ? makeFraction(1n, 10n) : options.epsilon,
  );
  if (circuitIds.length !== 2 || eyeIds.length !== 2) {
    throw new Error("a type-B compiler needs two incident paths");
  }
  const entries = [];
  for (let index = 0; index < 2; index += 1) {
    const circuitId = circuitIds[index];
    const circuit = circuits[circuitId];
    const visit = matchingVisitExact(circuit, vertex, vertex.switchPairs[index]);
    if (visit === null) return { ok: false, stage: "type-B-visit-inconsistent" };
    entries.push({ circuitId, circuit, visit, eyeId: eyeIds[index] });
  }
  const common = exactPointFromIntegers(vertex.x, vertex.y);
  const switchLevel = exactLinearForm(tau, common.x, common.y);
  const constraints = [];
  const diagnostics = [];
  for (const side of [-1n, 1n]) {
    const level = fractionAdd(switchLevel, fractionScale(epsilon, side));
    const paths = [];
    for (const entry of entries) {
      const incoming = entry.visit[1];
      const outgoing = entry.visit[2];
      const candidates = [
        {
          port: incoming,
          vx: BigInt(arr.portVectorX[incoming]),
          vy: BigInt(arr.portVectorY[incoming]),
        },
        {
          port: outgoing,
          vx: BigInt(arr.portVectorX[outgoing]),
          vy: BigInt(arr.portVectorY[outgoing]),
        },
      ].filter((ray) => {
        const sweep = tau[0] * ray.vx + tau[1] * ray.vy;
        return sweep !== 0n && (sweep < 0n ? -1n : 1n) === side;
      });
      if (candidates.length !== 1) {
        return { ok: false, stage: "type-B-ray-inconsistent" };
      }
      const ray = candidates[0];
      const raySweep = tau[0] * ray.vx + tau[1] * ray.vy;
      const alpha = fractionDivide(fractionScale(epsilon, side), makeFraction(raySweep));
      const near = {
        x: fractionAdd(common.x, fractionScale(alpha, ray.vx)),
        y: fractionAdd(common.y, fractionScale(alpha, ray.vy)),
      };
      const nearSigma = exactLinearForm(sigma, near.x, near.y);
      const companionCoefficients = new Array(dimension).fill(0n);
      let companionConstant;
      const boundaries = eyeCircuits[entry.eyeId];
      if (boundaries.length === 1) {
        const shiftX = BigInt(vertex.x) - BigInt(entry.visit[3]);
        const shiftY = BigInt(vertex.y) - BigInt(entry.visit[4]);
        const hits = finiteSliceHitsExact(
          entry.circuit, shiftX, shiftY, level, tau,
        );
        if (hits.length !== 2) {
          return { ok: false, stage: "type-B-disk-companion-inconsistent" };
        }
        const nearIndex = hits.findIndex((hit) => exactPointEqual(hit, near));
        if (nearIndex < 0) {
          return { ok: false, stage: "type-B-near-hit-inconsistent" };
        }
        const companion = hits[1 - nearIndex];
        companionConstant = exactLinearForm(sigma, companion.x, companion.y);
      } else {
        const descriptor = lookupAnnularDescriptor(annularByEye, entry.eyeId);
        if (descriptor === null || descriptor.ok !== true) {
          return { ok: false, stage: "type-B-missing-annular-descriptor" };
        }
        const period = BigInt(arr.L);
        const numeratorX = BigInt(vertex.x) - BigInt(entry.visit[3]);
        const numeratorY = BigInt(vertex.y) - BigInt(entry.visit[4]);
        if (numeratorX % period !== 0n || numeratorY % period !== 0n) {
          return { ok: false, stage: "type-B-nondeck-visit" };
        }
        const baseX = numeratorX / period;
        const baseY = numeratorY / period;
        let mateId;
        let orientation;
        if (entry.circuitId === descriptor.left) {
          mateId = descriptor.right;
          orientation = 1n;
        } else if (entry.circuitId === descriptor.right) {
          mateId = descriptor.left;
          orientation = -1n;
        } else return { ok: false, stage: "type-B-eye-boundary-inconsistent" };
        const phaseX = orientation * descriptor.residue * descriptor.phaseStep[0];
        const phaseY = orientation * descriptor.residue * descriptor.phaseStep[1];
        const hit = liftedLineHitExact(
          arr,
          circuits[mateId],
          baseX + phaseX,
          baseY + phaseY,
          level,
          tau,
        );
        if (hit === null) {
          return { ok: false, stage: "type-B-annular-companion-inconsistent" };
        }
        companionConstant = exactLinearForm(sigma, hit.x, hit.y);
        const coefficient = sigma[0] * period * orientation
          * descriptor.fixedSlicePeriod[0]
          + sigma[1] * period * orientation * descriptor.fixedSlicePeriod[1];
        companionCoefficients[descriptor.variable] = coefficient;
      }
      paths.push({ nearSigma, companionConstant, companionCoefficients });
    }
    const nearDifference = fractionSubtract(paths[1].nearSigma, paths[0].nearSigma);
    const nearSign = exactSign(nearDifference);
    if (nearSign === 0) return { ok: false, stage: "type-B-near-order-tie" };
    const farConstant = fractionSubtract(
      paths[1].companionConstant,
      paths[0].companionConstant,
    );
    const farCoefficients = paths[1].companionCoefficients.map(
      (entry, index) => fractionSubtract(
        asFraction(entry), asFraction(paths[0].companionCoefficients[index]),
      ),
    );
    const signedCoefficients = farCoefficients.map(
      (entry) => fractionScale(entry, BigInt(nearSign)),
    );
    const signedConstant = fractionScale(farConstant, BigInt(nearSign));
    constraints.push({
      coefficients: signedCoefficients,
      operator: ">",
      rhs: fractionNegate(signedConstant),
      label: `type-B-same-order-${side < 0n ? "before" : "after"}`,
    });
    diagnostics.push({ side, nearDifference, farCoefficients, farConstant });
  }
  return {
    ok: true,
    stage: "type-B-affine",
    constraints,
    diagnostics,
    inwardAllowed: true,
  };
}

function compilePhaseBalance(descriptors, diskContribution = 0n, target = 0n) {
  if (!Array.isArray(descriptors) || descriptors.length === 0) {
    throw new Error("phase balance needs at least one annular descriptor");
  }
  const dimension = descriptors[0].dimension;
  const coefficients = new Array(dimension).fill(0n);
  let constant = BigInt(diskContribution);
  for (const descriptor of descriptors) {
    if (descriptor.ok !== true || descriptor.dimension !== dimension) {
      throw new Error("inconsistent annular descriptor in phase balance");
    }
    const sign = BigInt(descriptor.color);
    coefficients[descriptor.variable] += sign * descriptor.multiplicity.coefficient;
    constant += sign * descriptor.multiplicity.constant;
  }
  return {
    coefficients,
    operator: "=",
    rhs: BigInt(target) - constant,
    label: "base-balance",
    constant,
    target: BigInt(target),
  };
}

function solveCompiledPhaseSystem(options) {
  const { descriptors, typeB = [] } = options;
  if (!Array.isArray(descriptors) || descriptors.length === 0) {
    throw new Error("compiled phase system needs annular descriptors");
  }
  for (const descriptor of descriptors) {
    if (descriptor.ok !== true) {
      const status = typeof descriptor.stage === "string"
          && descriptor.stage.includes("unknown") ? "UNKNOWN" : "EMPTY";
      return {
        status,
        stage: descriptor.stage,
        reason: descriptor.reason || descriptor.stage,
        count: 0n,
      };
    }
  }
  for (const compiled of typeB) {
    if (compiled.ok !== true) {
      const status = typeof compiled.stage === "string"
          && compiled.stage.includes("unknown") ? "UNKNOWN" : "EMPTY";
      return {
        status,
        stage: compiled.stage,
        reason: compiled.reason || compiled.stage,
        count: 0n,
      };
    }
  }
  const dimension = descriptors[0].dimension;
  const polyhedron = new IntegerPolyhedron(dimension);
  for (const descriptor of descriptors) {
    if (descriptor.dimension !== dimension) {
      throw new Error("inconsistent annular descriptor dimension");
    }
    installCompiledConstraints(polyhedron, descriptor.constraints);
  }
  for (const compiled of typeB) {
    installCompiledConstraints(polyhedron, compiled.constraints);
  }
  if (options.requireBalance !== false) {
    const balance = compilePhaseBalance(
      descriptors,
      options.diskContribution === undefined ? 0n : options.diskContribution,
      options.balanceTarget === undefined ? 0n : options.balanceTarget,
    );
    installCompiledConstraints(polyhedron, [balance]);
  }
  const result = polyhedron.solve(options.solverOptions || {});
  return { ...result, polyhedron };
}

function stringifyExactResult(value) {
  return JSON.stringify(value, (_, entry) => (
    typeof entry === "bigint" ? entry.toString() : entry
  ), 2);
}

module.exports = {
  IntegerPolyhedron,
  PolyhedralLimitError,
  makeFraction,
  normalizeLinearConstraint,
  phaseResidueData,
  liftedLineHitExact,
  forbiddenAnnularPhasesExact,
  translatedDiskMultiplicityExact,
  compileAnnularEyeResidue,
  compileTypeBSameOrder,
  compilePhaseBalance,
  solveCompiledPhaseSystem,
  installCompiledConstraints,
  solveNormalizedPolyhedron,
  stringifyExactResult,
  // Exported for certificate/unit-test verification, not as a stable API.
  _internal: {
    findRecessionDirection,
    pointSatisfiesInteger,
    raySatisfiesRecession,
  },
};
