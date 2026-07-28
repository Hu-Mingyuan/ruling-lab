#!/usr/bin/env node
"use strict";

/*
 * Annular-aware smoothing audit for the P^2,O(d) three-geodesic-family
 * arrangement.  This is deliberately separate from the all-disk search.
 * A partial circuit may close either as
 *
 *   - a contractible disk boundary with exactly two V tips, or
 *   - a primitive essential annular boundary with zero V tips.
 *
 * At a leaf, essential circuits are grouped by unoriented homology and
 * opposite coorientation side.  Every perfect matching and both colors of
 * every proposed annular eye are then checked against the local A/B/V/N
 * rules, embedded-annulus boundary disjointness, connectedness, signed
 * orientability, and the requested genus/Euler characteristic.  The script reports subsequent
 * topology/balance stages separately instead of promoting unresolved phase
 * choices to rulings.
 */

const { buildArrangement, buildOrder, maskHex } = require("./p2o4_smoothing_search.js");

function buildAnnularOrder(arr, name) {
  if (name !== "yrev") return buildOrder(arr, name);
  return arr.vertices.map((vertex) => vertex.index).sort((a, b) => (
    arr.vertices[b].y - arr.vertices[a].y
      || arr.vertices[a].x - arr.vertices[b].x
      || a - b
  ));
}

function translationInverses(arr) {
  const d = arr.d;
  const byLabel = new Map(arr.vertices.map((vertex) => (
    [`${vertex.type}:${vertex.i}:${vertex.j}`, vertex.index]
  )));
  const inverses = [];
  const mod = (value) => ((value % d) + d) % d;
  for (let a = 0; a < d; a += 1) {
    for (let b = 0; b < d; b += 1) {
      const inverse = new Int16Array(arr.vertices.length);
      for (const vertex of arr.vertices) {
        let i;
        let j;
        if (vertex.type === "AB") {
          i = mod(vertex.i + b);
          j = mod(vertex.j + a);
        } else if (vertex.type === "AC") {
          i = mod(vertex.i + b);
          j = mod(vertex.j + a + b);
        } else {
          i = mod(vertex.i + a);
          j = mod(vertex.j + a + b);
        }
        const target = byLabel.get(`${vertex.type}:${i}:${j}`);
        if (target === undefined) throw new Error("translation lost a crossing");
        inverse[target] = vertex.index;
      }
      inverses.push(inverse);
    }
  }
  return inverses;
}

// Return zero for a noncanonical mask.  For a canonical mask return the
// number of distinct translated masks, so auditing one representative and
// multiplying by this weight counts the whole mask orbit exactly, even when
// the mask has a translation stabilizer.
function canonicalTranslationWeight(choice, inverses) {
  const orbit = new Set();
  for (const inverse of inverses) {
    let translatedMask = 0n;
    let comparison = 0;
    for (let index = 0; index < choice.length; index += 1) {
      const translated = choice[inverse[index]];
      if (translated) translatedMask |= 1n << BigInt(index);
      if (comparison === 0 && choice[index] !== translated) {
        comparison = choice[index] < translated ? -1 : 1;
      }
    }
    if (comparison > 0) return 0;
    orbit.add(translatedMask.toString(16));
  }
  return orbit.size;
}

function horizontalResidueFeasibility(arr, order) {
  const d = arr.d;
  const stride = d * d;
  const width = (arr.requiredH + 1) * stride;
  const feasible = Array.from({ length: order.length + 1 }, () => new Uint8Array(width));
  feasible[order.length][0] = 1;
  for (let index = order.length - 1; index >= 0; index -= 1) {
    const vertex = arr.vertices[order[index]];
    feasible[index].set(feasible[index + 1]);
    if (vertex.orientation !== "H") continue;
    if (vertex.type !== "BC") throw new Error("centroid symmetry expects BC horizontal crossings");
    for (let count = 0; count < arr.requiredH; count += 1) {
      for (let ri = 0; ri < d; ri += 1) {
        for (let rj = 0; rj < d; rj += 1) {
          if (!feasible[index + 1][count * stride + ri * d + rj]) continue;
          const nextI = (ri + vertex.i) % d;
          const nextJ = (rj + vertex.j) % d;
          feasible[index][(count + 1) * stride + nextI * d + nextJ] = 1;
        }
      }
    }
  }
  return feasible;
}

function gcd(a, b) {
  a = Math.abs(a);
  b = Math.abs(b);
  while (b !== 0) [a, b] = [b, a % b];
  return a;
}

class AnnularRollbackDSU {
  constructor(n) {
    this.parent = new Int16Array(n);
    this.size = new Int16Array(n);
    this.potX = new Int16Array(n);
    this.potY = new Int16Array(n);
    this.tips = new Int8Array(n);
    // Signed area in the canonical circuit orientation whose outgoing ports
    // have coorientation sign +1.  A V tip fixes this sign before closure.
    this.areaSign = new Int8Array(n);
    // 0=open, 1=closed disk, 2=closed essential.
    this.closedKind = new Int8Array(n);
    this.diskColor = new Int8Array(n);
    this.diskInteriorSide = new Int8Array(n);
    this.diskTraceSide = new Int8Array(n);
    this.closedRawHx = new Int32Array(n);
    this.closedRawHy = new Int32Array(n);
    this.closedTraceSide = new Int8Array(n);
    this.forbiddenMask = Array(n).fill(0n);
    // Intrusive member lists let a union inspect only the smaller component.
    // They are used by the dynamic zero-lift inequalities below.  Each port
    // can be an endpoint of at most one active local inequality because a
    // crossing is either a switch or a nonswitch, never both.
    this.memberHead = new Int16Array(n);
    this.memberTail = new Int16Array(n);
    this.memberNext = new Int16Array(n);
    this.zeroConstraintMate = new Int16Array(n);
    this.zeroConstraintKind = new Int8Array(n);
    this.memberNext.fill(-1);
    this.zeroConstraintMate.fill(-1);
    for (let i = 0; i < n; i += 1) {
      this.parent[i] = i;
      this.size[i] = 1;
      this.memberHead[i] = i;
      this.memberTail[i] = i;
    }
    const capacity = 4 * n + 16;
    this.stackKind = new Int8Array(capacity);
    this.stackChild = new Int16Array(capacity);
    this.stackRoot = new Int16Array(capacity);
    this.stackSize = new Int16Array(capacity);
    this.stackTips = new Int8Array(capacity);
    this.stackAreaSign = new Int8Array(capacity);
    this.stackClosed = new Int8Array(capacity);
    this.stackMaskRoot = Array(capacity).fill(0n);
    this.stackMaskChild = Array(capacity).fill(0n);
    this.stackMemberTail = new Int16Array(capacity);
    this.top = 0;
    this.components = n;
    this.closedCycles = 0;
    this.diskCycles = 0;
    this.essentialCycles = 0;
    // Number of open path components carrying 0, 1, or 2 V tips.
    this.openTipCounts = new Int16Array(3);
    this.openTipCounts[0] = n;
    this.lastFailure = "";
  }

  find(a) {
    let x = 0;
    let y = 0;
    let node = a;
    while (this.parent[node] !== node) {
      x += this.potX[node];
      y += this.potY[node];
      node = this.parent[node];
    }
    return [node, x, y];
  }

  relation(a, b) {
    const fa = this.find(a);
    const fb = this.find(b);
    if (fa[0] !== fb[0]) return null;
    return [fb[1] - fa[1], fb[2] - fa[2]];
  }

  activateZeroConstraint(a, b, kind) {
    if (this.zeroConstraintMate[a] !== -1 || this.zeroConstraintMate[b] !== -1) {
      throw new Error("a port acquired two local zero-lift constraints");
    }
    this.zeroConstraintMate[a] = b;
    this.zeroConstraintMate[b] = a;
    this.zeroConstraintKind[a] = kind;
    this.zeroConstraintKind[b] = kind;
    const relation = this.relation(a, b);
    if (relation !== null && relation[0] === 0 && relation[1] === 0) {
      this.lastFailure = kind === 2
        ? "same-lift-self-switch"
        : "nonswitch-self-intersection";
      return false;
    }
    return true;
  }

  deactivateZeroConstraint(a, b) {
    if (this.zeroConstraintMate[a] !== b || this.zeroConstraintMate[b] !== a) {
      throw new Error("mismatched zero-lift constraint rollback");
    }
    this.zeroConstraintMate[a] = -1;
    this.zeroConstraintMate[b] = -1;
    this.zeroConstraintKind[a] = 0;
    this.zeroConstraintKind[b] = 0;
  }

  newlyViolatesZeroConstraint(childRoot) {
    for (let port = this.memberHead[childRoot]; port !== -1; port = this.memberNext[port]) {
      const mate = this.zeroConstraintMate[port];
      if (mate === -1) continue;
      const relation = this.relation(port, mate);
      if (relation === null || relation[0] !== 0 || relation[1] !== 0) continue;
      this.lastFailure = this.zeroConstraintKind[port] === 2
        ? "same-lift-self-switch"
        : "nonswitch-self-intersection";
      return true;
    }
    return false;
  }

  unite(a, b, dx, dy, addedTips = 0, requiredAreaSign = 0, closureTraceSide = 0) {
    const fa = this.find(a);
    const fb = this.find(b);
    let ra = fa[0];
    let rb = fb[0];
    let relX = fa[1] + dx - fb[1];
    let relY = fa[2] + dy - fb[2];
    if (ra === rb) {
      if (this.closedKind[ra] !== 0) {
        this.lastFailure = "duplicate-cycle";
        return false;
      }
      const newTips = this.tips[ra] + addedTips;
      if (requiredAreaSign !== 0 && this.areaSign[ra] !== 0
          && requiredAreaSign !== this.areaSign[ra]) {
        this.lastFailure = "conflicting-V-convexity";
        return false;
      }
      let kind;
      if (relX === 0 && relY === 0) {
        if (newTips !== 2) {
          this.lastFailure = "disk-wrong-tip-count";
          return false;
        }
        kind = 1;
      } else {
        if (newTips !== 0) {
          this.lastFailure = "essential-has-tip";
          return false;
        }
        if (gcd(relX, relY) !== 1) {
          this.lastFailure = "essential-nonprimitive";
          return false;
        }
        kind = 2;
      }
      const at = this.top++;
      this.stackKind[at] = 0;
      this.stackRoot[at] = ra;
      this.stackTips[at] = this.tips[ra];
      this.stackAreaSign[at] = this.areaSign[ra];
      this.stackClosed[at] = this.closedKind[ra];
      this.openTipCounts[this.tips[ra]] -= 1;
      this.tips[ra] = newTips;
      if (requiredAreaSign !== 0) this.areaSign[ra] = requiredAreaSign;
      this.closedKind[ra] = kind;
      this.closedRawHx[ra] = relX;
      this.closedRawHy[ra] = relY;
      this.closedTraceSide[ra] = closureTraceSide;
      this.lastClosedRoot = ra;
      this.closedCycles += 1;
      if (kind === 1) this.diskCycles += 1;
      else this.essentialCycles += 1;
      return true;
    }
    if (this.closedKind[ra] !== 0 || this.closedKind[rb] !== 0) {
      this.lastFailure = "join-closed-cycle";
      return false;
    }
    const newTips = this.tips[ra] + this.tips[rb] + addedTips;
    if (newTips > 2) {
      this.lastFailure = "tip-overflow";
      return false;
    }
    let newAreaSign = this.areaSign[ra] || this.areaSign[rb] || requiredAreaSign;
    for (const sign of [this.areaSign[ra], this.areaSign[rb], requiredAreaSign]) {
      if (sign !== 0 && sign !== newAreaSign) {
        this.lastFailure = "conflicting-V-convexity";
        return false;
      }
    }
    if ((this.forbiddenMask[ra] & this.forbiddenMask[rb]) !== 0n) {
      this.lastFailure = "one-eye-at-switch";
      return false;
    }
    if (this.size[ra] < this.size[rb]) {
      [ra, rb] = [rb, ra];
      relX = -relX;
      relY = -relY;
    }
    const at = this.top++;
    this.stackKind[at] = 1;
    this.stackChild[at] = rb;
    this.stackRoot[at] = ra;
    this.stackSize[at] = this.size[ra];
    this.stackTips[at] = this.tips[ra];
    this.stackAreaSign[at] = this.areaSign[ra];
    this.stackClosed[at] = this.closedKind[ra];
    this.stackMaskRoot[at] = this.forbiddenMask[ra];
    this.stackMemberTail[at] = this.memberTail[ra];
    this.openTipCounts[this.tips[ra]] -= 1;
    this.openTipCounts[this.tips[rb]] -= 1;
    this.openTipCounts[newTips] += 1;
    this.parent[rb] = ra;
    this.potX[rb] = relX;
    this.potY[rb] = relY;
    this.size[ra] += this.size[rb];
    this.tips[ra] = newTips;
    this.areaSign[ra] = newAreaSign;
    this.forbiddenMask[ra] |= this.forbiddenMask[rb];
    this.memberNext[this.memberTail[ra]] = this.memberHead[rb];
    this.memberTail[ra] = this.memberTail[rb];
    this.components -= 1;
    if (this.newlyViolatesZeroConstraint(rb)) return false;
    return true;
  }

  addInequality(a, b, bit) {
    const ra = this.find(a)[0];
    const rb = this.find(b)[0];
    if (ra === rb) {
      this.lastFailure = "one-eye-at-switch";
      return false;
    }
    const at = this.top++;
    this.stackKind[at] = 2;
    this.stackRoot[at] = ra;
    this.stackChild[at] = rb;
    this.stackMaskRoot[at] = this.forbiddenMask[ra];
    this.stackMaskChild[at] = this.forbiddenMask[rb];
    this.forbiddenMask[ra] |= bit;
    this.forbiddenMask[rb] |= bit;
    return true;
  }

  rollback(mark) {
    while (this.top > mark) {
      const at = --this.top;
      const root = this.stackRoot[at];
      if (this.stackKind[at] === 0) {
        const kind = this.closedKind[root];
        if (kind === 1) this.diskCycles -= 1;
        else if (kind === 2) this.essentialCycles -= 1;
        this.closedCycles -= 1;
        this.openTipCounts[this.stackTips[at]] += 1;
        this.tips[root] = this.stackTips[at];
        this.areaSign[root] = this.stackAreaSign[at];
        this.closedKind[root] = this.stackClosed[at];
        this.diskColor[root] = 0;
        this.diskInteriorSide[root] = 0;
        this.diskTraceSide[root] = 0;
        this.closedRawHx[root] = 0;
        this.closedRawHy[root] = 0;
        this.closedTraceSide[root] = 0;
      } else if (this.stackKind[at] === 1) {
        const child = this.stackChild[at];
        const mergedTips = this.tips[root];
        this.parent[child] = child;
        this.potX[child] = 0;
        this.potY[child] = 0;
        this.size[root] = this.stackSize[at];
        this.tips[root] = this.stackTips[at];
        this.areaSign[root] = this.stackAreaSign[at];
        this.closedKind[root] = this.stackClosed[at];
        this.forbiddenMask[root] = this.stackMaskRoot[at];
        const oldTail = this.stackMemberTail[at];
        this.memberNext[oldTail] = -1;
        this.memberTail[root] = oldTail;
        this.openTipCounts[mergedTips] -= 1;
        this.openTipCounts[this.stackTips[at]] += 1;
        this.openTipCounts[this.tips[child]] += 1;
        this.components += 1;
      } else {
        const child = this.stackChild[at];
        this.forbiddenMask[root] = this.stackMaskRoot[at];
        this.forbiddenMask[child] = this.stackMaskChild[at];
      }
    }
  }
}

function addExternalEdges(arr, dsu) {
  for (const edge of arr.edges) {
    if (!dsu.unite(edge.a, edge.b, edge.sx, edge.sy)) {
      throw new Error(`external edge failed: ${dsu.lastFailure}`);
    }
  }
  dsu.top = 0;
}

function addLocal(arr, dsu, vertexId, switched, closedRoots = null) {
  const vertex = arr.vertices[vertexId];
  const pairs = switched ? vertex.switchPairs : vertex.nonswitchPairs;
  const tip = switched && vertex.orientation === "V" ? 1 : 0;
  for (const [left, right] of pairs) {
    const incoming = arr.portSign[left] < 0 ? left : right;
    const outgoing = arr.portSign[left] > 0 ? left : right;
    let requiredAreaSign = 0;
    if (tip) {
      const turn = (-arr.portVectorX[incoming]) * arr.portVectorY[outgoing]
        - (-arr.portVectorY[incoming]) * arr.portVectorX[outgoing];
      requiredAreaSign = Math.sign(turn);
      if (requiredAreaSign === 0) throw new Error("vertical tip has zero turn");
    }
    const before = dsu.closedCycles;
    if (!dsu.unite(incoming, outgoing, 0, 0, tip, requiredAreaSign, 1)) return false;
    if (closedRoots !== null && dsu.closedCycles !== before) {
      closedRoots.push(dsu.lastClosedRoot);
    }
  }
  return true;
}

function closedCircuitGeometry(arr, dsu, choice, root) {
  const start = dsu.memberHead[root];
  if (start < 0) throw new Error("closed disk root has no port");
  let current = start;
  let liftX = 0;
  let liftY = 0;
  let area2 = 0;
  let traceSide = 0;
  for (let guard = 0; ; guard += 1) {
    if (guard > arr.portCount) throw new Error("closed circuit trace overflow");
    if (guard > 0 && current === start) break;
    const vertex = arr.vertices[arr.portVertex[current]];
    const pairs = choice[vertex.index] ? vertex.switchPairs : vertex.nonswitchPairs;
    let outgoing = -1;
    for (const [left, right] of pairs) {
      if (left === current) outgoing = right;
      else if (right === current) outgoing = left;
    }
    if (outgoing < 0) throw new Error("closed circuit trace lost its local pair");
    const side = arr.portSign[outgoing];
    if (traceSide === 0) traceSide = side;
    else if (traceSide !== side) throw new Error("closed circuit has inconsistent trace side");
    const ax = vertex.x + liftX * arr.L;
    const ay = vertex.y + liftY * arr.L;
    liftX += arr.portShiftX[outgoing];
    liftY += arr.portShiftY[outgoing];
    const next = arr.portTwin[outgoing];
    const target = arr.vertices[arr.portVertex[next]];
    const bx = target.x + liftX * arr.L;
    const by = target.y + liftY * arr.L;
    area2 += ax * by - ay * bx;
    current = next;
  }
  return { liftX, liftY, area2, traceSide };
}

function closedDiskGeometry(arr, dsu, choice, root) {
  const geometry = closedCircuitGeometry(arr, dsu, choice, root);
  const { liftX, liftY, area2, traceSide } = geometry;
  if (liftX !== 0 || liftY !== 0 || area2 === 0) {
    throw new Error("invalid closed disk geometry");
  }
  const interiorSide = area2 > 0 ? 1 : -1;
  const color = traceSide === interiorSide ? 1 : -1;
  return { color, interiorSide, traceSide };
}

function closedEssentialSignatures(arr, dsu, choice, closedRoots) {
  const signatures = [];
  for (const root of closedRoots) {
    if (dsu.closedKind[root] !== 2) continue;
    const rawHx = dsu.closedRawHx[root];
    const rawHy = dsu.closedRawHy[root];
    const traceSide = dsu.closedTraceSide[root];
    if (gcd(rawHx, rawHy) !== 1 || traceSide === 0) {
      throw new Error("closed essential circuit is not primitive");
    }
    const [hx, hy, canonicalSide] = canonicalHomology(
      rawHx, rawHy, traceSide,
    );
    signatures.push({ key: `${hx},${hy}`, canonicalSide });
  }
  return signatures;
}

function annotateClosedDisks(arr, dsu, choice, closedRoots, verifyGeometry = false) {
  let addedDisk = false;
  for (const root of closedRoots) {
    if (dsu.closedKind[root] !== 1) continue;
    // A disk has exactly two V tips.  Each tip fixes the signed area in the
    // canonical orientation (outgoing-port sign +1), and unite() has already
    // rejected disagreeing signs.  Thus areaSign is exactly the disk color;
    // choosing canonical traceSide=+1 gives interiorSide=areaSign.  No full
    // circuit trace is needed on the hot path.
    const canonicalAreaSign = dsu.areaSign[root];
    if (canonicalAreaSign === 0) {
      throw new Error("closed disk did not acquire a V-tip area sign");
    }
    if (verifyGeometry) {
      const geometry = closedDiskGeometry(arr, dsu, choice, root);
      if (geometry.color !== canonicalAreaSign) {
        throw new Error("V-tip area sign disagrees with closed disk geometry");
      }
    }
    dsu.diskColor[root] = canonicalAreaSign;
    dsu.diskInteriorSide[root] = canonicalAreaSign;
    dsu.diskTraceSide[root] = 1;
    addedDisk = true;
  }
  return addedDisk;
}

function horizontalConstraint(arr, vertex) {
  return vertex.switchPairs.map(([left, right]) => {
    const incoming = arr.portSign[left] < 0 ? left : right;
    const outgoing = arr.portSign[left] > 0 ? left : right;
    const turn = (-arr.portVectorX[incoming]) * arr.portVectorY[outgoing]
      - (-arr.portVectorY[incoming]) * arr.portVectorX[outgoing];
    return [left, Math.sign(turn)];
  });
}

// The active A/B pictures at a horizontal switch say that at least one of
// the two incident eyes is convex.  Once both partial components have signs
// forced by their V tips, this clause is final and can be rejected early.
function violatesHorizontalClause(dsu, constraints) {
  for (const [[left, leftWanted], [right, rightWanted]] of constraints) {
    const leftSign = dsu.areaSign[dsu.find(left)[0]];
    const rightSign = dsu.areaSign[dsu.find(right)[0]];
    if (leftSign !== 0 && rightSign !== 0
        && leftSign !== leftWanted && rightSign !== rightWanted) return true;
  }
  return false;
}

function diskTurnAtPair(arr, pair, traceSide) {
  const outgoing = arr.portSign[pair[0]] === traceSide ? pair[0] : pair[1];
  const incoming = outgoing === pair[0] ? pair[1] : pair[0];
  return (-arr.portVectorX[incoming]) * arr.portVectorY[outgoing]
    - (-arr.portVectorY[incoming]) * arr.portVectorX[outgoing];
}

// Exact early version of the disk-disk part of localAudit.  It is invoked
// only when a new disk closes, so every local datum used here is final.
function partialClosedDiskAudit(
  arr, dsu, choice, assigned, allowQuotientSelfSwitch,
  closedRoots = null, vertexStamp = null, stamp = 0,
) {
  let vertices = arr.vertices;
  if (closedRoots !== null) {
    const candidateIds = [];
    for (const root of closedRoots) {
      if (dsu.closedKind[root] !== 1) continue;
      for (let port = dsu.memberHead[root]; port !== -1; port = dsu.memberNext[port]) {
        const vertexId = arr.portVertex[port];
        if (vertexStamp[vertexId] === stamp) continue;
        vertexStamp[vertexId] = stamp;
        candidateIds.push(vertexId);
      }
    }
    vertices = candidateIds.map((id) => arr.vertices[id]);
  }
  for (const vertex of vertices) {
    if (!assigned[vertex.index]) continue;
    const switched = choice[vertex.index] !== 0;
    const pairs = switched ? vertex.switchPairs : vertex.nonswitchPairs;
    const roots = pairs.map((pair) => dsu.find(pair[0])[0]);
    if (dsu.closedKind[roots[0]] !== 1 || dsu.closedKind[roots[1]] !== 1) continue;
    const colors = roots.map((root) => dsu.diskColor[root]);
    if (!switched) {
      const overlap = [...occupiedSectors(arr, vertex, pairs[0][0], colors[0])]
        .filter((sector) => occupiedSectors(arr, vertex, pairs[1][0], colors[1]).has(sector));
      if (overlap.length !== 1) return "invalid-nonswitch";
      continue;
    }
    if (roots[0] === roots[1] && !allowQuotientSelfSwitch) return "one-eye-at-switch";
    const convex = pairs.map((pair, index) => (
      diskTurnAtPair(arr, pair, dsu.diskTraceSide[roots[index]])
        * dsu.diskInteriorSide[roots[index]] > 0
    ));
    if (vertex.orientation === "V") {
      if (colors[0] === colors[1] || !convex[0] || !convex[1]) return "invalid-V";
    } else if (colors[0] !== colors[1]) {
      if (!convex[0] || !convex[1]) return "invalid-A";
    } else if (convex[0] === convex[1]) {
      return "invalid-B";
    }
  }
  return null;
}

function forcedNonswitch(dsu, constraints) {
  for (const [left, right] of constraints) {
    const relation = dsu.relation(left, right);
    if (relation !== null && relation[0] === 0 && relation[1] === 0) return true;
  }
  return false;
}

function canonicalHomology(x, y, side) {
  if (x < 0 || (x === 0 && y < 0)) return [-x, -y, -side];
  return [x, y, side];
}

function modInteger(value, period) {
  value %= period;
  return value < 0 ? value + period : value;
}

function traceCircuits(arr, choice) {
  const seen = new Uint8Array(arr.portCount);
  const circuitOf = new Int16Array(arr.portCount);
  const turnAt = new Int8Array(arr.portCount);
  circuitOf.fill(-1);
  const circuits = [];
  for (let start = 0; start < arr.portCount; start += 1) {
    if (seen[start]) continue;
    const id = circuits.length;
    let current = start;
    let liftX = 0;
    let liftY = 0;
    let tips = 0;
    let area2 = 0;
    const polygon = [];
    const segments = [];
    const visits = [];
    // Record every lifted visit, not only nonswitch visits.  A later
    // nonswitch at an already visited lift is a genuine self-intersection
    // even when the earlier visit happened to be switched.
    const liftedVertices = new Set();
    for (let guard = 0; ; guard += 1) {
      if (guard > arr.portCount) return null;
      if (guard > 0 && current === start) break;
      if (seen[current]) return null;
      const v = arr.portVertex[current];
      const vertex = arr.vertices[v];
      const switched = choice[v] !== 0;
      const key = `${v}:${liftX}:${liftY}`;
      if (!switched && liftedVertices.has(key)) return null;
      liftedVertices.add(key);
      const pairs = switched ? vertex.switchPairs : vertex.nonswitchPairs;
      let outgoing = -1;
      for (const [left, right] of pairs) {
        if (left === current) outgoing = right;
        else if (right === current) outgoing = left;
      }
      if (outgoing < 0) return null;
      const ax = vertex.x + liftX * arr.L;
      const ay = vertex.y + liftY * arr.L;
      polygon.push([ax, ay]);
      visits.push([v, current, outgoing, ax, ay]);
      seen[current] = 1;
      seen[outgoing] = 1;
      circuitOf[current] = id;
      circuitOf[outgoing] = id;
      const turn = (-arr.portVectorX[current]) * arr.portVectorY[outgoing]
        - (-arr.portVectorY[current]) * arr.portVectorX[outgoing];
      turnAt[current] = turn;
      turnAt[outgoing] = turn;
      if (switched && vertex.orientation === "V") tips += 1;
      liftX += arr.portShiftX[outgoing];
      liftY += arr.portShiftY[outgoing];
      const next = arr.portTwin[outgoing];
      const target = arr.vertices[arr.portVertex[next]];
      const bx = target.x + liftX * arr.L;
      const by = target.y + liftY * arr.L;
      segments.push([ax, ay, bx, by]);
      area2 += ax * by - ay * bx;
      current = next;
    }
    const signs = new Set(visits.map((visit) => arr.portSign[visit[2]]));
    if (signs.size !== 1) return null;
    const traceSide = [...signs][0];
    if (liftX === 0 && liftY === 0) {
      if (tips !== 2 || area2 === 0) return null;
      const areaSign = area2 > 0 ? 1 : -1;
      const color = traceSide === areaSign ? 1 : -1;
      circuits.push({
        id, kind: "disk", color, area2, interiorSide: areaSign,
        polygon, segments, visits, hx: 0, hy: 0, rawHx: 0, rawHy: 0,
        traceSide,
      });
    } else {
      if (tips !== 0 || gcd(liftX, liftY) !== 1) return null;
      const [hx, hy, canonicalSide] = canonicalHomology(liftX, liftY, traceSide);
      circuits.push({
        id, kind: "essential", color: 0, area2: 0, interiorSide: 0,
        polygon, segments, visits, hx, hy, rawHx: liftX, rawHy: liftY,
        traceSide, canonicalSide,
      });
    }
  }
  return { circuits, circuitOf, turnAt };
}

/*
 * Exact annular-strip audit in a transverse sweep.
 *
 * Put tau=2x+y and sigma=-x+2y.  Every future A/B/C edge has strictly
 * positive tau increment.  Consequently an essential zero-tip smoothing
 * circuit is a monotone proper line after choosing its future orientation.
 * At a generic tau-slice its deck translates have a well-defined cyclic
 * sigma-order (period 5L).  For primitive homology h, a quotient circuit
 * contributes q=|2h_x+h_y| entries to that order; the entry phase is the
 * coset det(h,t) mod q of the translating deck vector t.
 *
 * An embedded annular eye of the chosen color has no relative-deck-phase
 * ambiguity: from either boundary, the color chooses an interior side, and
 * the first boundary met on that side must be its mate.  Requiring this on
 * one slice in every interval between crossing events simultaneously
 * (1) chooses a constant relative phase, (2) excludes a non-adjacent or
 * interleaved pairing, and (3) detects an intervening strand after any
 * nonswitch order exchange.
 */
function auditAnnularStrips(arr, trace, matching, colorMask) {
  if (matching.length === 0) return { ok: true, phases: [] };

  const essentialGroups = new Map();
  for (const circuit of trace.circuits) {
    if (circuit.kind !== "essential") continue;
    const key = `${circuit.hx},${circuit.hy}`;
    if (!essentialGroups.has(key)) essentialGroups.set(key, []);
    essentialGroups.get(key).push(circuit.id);
  }

  const mate = new Int16Array(trace.circuits.length);
  const pairIndex = new Int16Array(trace.circuits.length);
  mate.fill(-1);
  pairIndex.fill(-1);
  const sideFuture = new Int8Array(trace.circuits.length);
  for (let i = 0; i < matching.length; i += 1) {
    const color = ((colorMask >> BigInt(i)) & 1n) !== 0n ? 1 : -1;
    const [a, b] = matching[i];
    mate[a] = b;
    mate[b] = a;
    pairIndex[a] = i;
    pairIndex[b] = i;
    for (const id of [a, b]) {
      const circuit = trace.circuits[id];
      const rawQ = 2 * circuit.rawHx + circuit.rawHy;
      if (rawQ === 0) return { ok: false, stage: "annular-zero-sweep" };
      const futureSign = rawQ > 0 ? 1 : -1;
      sideFuture[id] = color * circuit.traceSide * futureSign;
    }
    if (sideFuture[a] !== -sideFuture[b]) {
      return { ok: false, stage: "annular-sides-not-facing" };
    }
  }

  // Crossing coordinates are integral in this model.  Midpoints of the
  // distinct tau residues avoid every event, including simultaneous ones.
  const residues = [...new Set(arr.vertices.map(
    (vertex) => modInteger(2 * vertex.x + vertex.y, arr.L),
  ))].sort((a, b) => a - b);
  const samples = [];
  for (let i = 0; i < residues.length; i += 1) {
    const left = residues[i];
    let right = residues[(i + 1) % residues.length];
    if (i + 1 === residues.length) right += arr.L;
    samples.push(modInteger((left + right) / 2, arr.L));
  }

  const phaseByPair = new Int16Array(matching.length);
  phaseByPair.fill(-1);

  for (const [key, ids] of essentialGroups.entries()) {
    const [canonicalHx, canonicalHy] = key.split(",").map(Number);
    const canonicalQ = 2 * canonicalHx + canonicalHy;
    if (canonicalQ === 0) return { ok: false, stage: "annular-zero-sweep" };
    const futureHx = canonicalQ > 0 ? canonicalHx : -canonicalHx;
    const futureHy = canonicalQ > 0 ? canonicalHy : -canonicalHy;
    const q = Math.abs(canonicalQ);

    // Normalize every fundamental lifted polygon to the common future
    // orientation.  Reversing a segment list does not change its deck phase;
    // it only makes tau strictly increasing for the slice calculation.
    const futureSegments = new Map();
    for (const id of ids) {
      const circuit = trace.circuits[id];
      const rawQ = 2 * circuit.rawHx + circuit.rawHy;
      let segments = circuit.segments;
      if (rawQ < 0) {
        segments = segments.slice().reverse().map(
          ([ax, ay, bx, by]) => [bx, by, ax, ay],
        );
      }
      if (segments.some(([ax, ay, bx, by]) => 2 * bx + by <= 2 * ax + ay)) {
        return { ok: false, stage: "annular-nonmonotone" };
      }
      futureSegments.set(id, segments);
    }

    for (const sample of samples) {
      const entries = [];
      for (const id of ids) {
        const phasesSeen = new Set();
        for (const [ax, ay, bx, by] of futureSegments.get(id)) {
          const tauA = 2 * ax + ay;
          const tauB = 2 * bx + by;
          // A translate (0,k) changes tau by kL and sigma by 2kL.
          // There is at most one such k because an arrangement edge has
          // tau-length < L; the while form also remains correct generically.
          const kMin = Math.floor((sample - tauB) / arr.L) + 1;
          const kMax = Math.ceil((sample - tauA) / arr.L) - 1;
          for (let k = kMin; k <= kMax; k += 1) {
            const translatedTauA = tauA + k * arr.L;
            const alpha = (sample - translatedTauA) / (tauB - tauA);
            if (!(alpha > 0 && alpha < 1)) continue;
            const sigmaA = -ax + 2 * ay;
            const sigmaB = -bx + 2 * by;
            const sigma = modInteger(
              sigmaA + alpha * (sigmaB - sigmaA) + 2 * k * arr.L,
              5 * arr.L,
            );
            const phase = modInteger(futureHx * k, q);
            if (phasesSeen.has(phase)) {
              return { ok: false, stage: "annular-duplicate-slice-phase" };
            }
            phasesSeen.add(phase);
            entries.push({ id, phase, sigma });
          }
        }
        if (phasesSeen.size !== q) {
          return { ok: false, stage: "annular-slice-count" };
        }
      }
      entries.sort((a, b) => a.sigma - b.sigma || a.id - b.id || a.phase - b.phase);
      for (let i = 1; i < entries.length; i += 1) {
        if (Math.abs(entries[i].sigma - entries[i - 1].sigma) < 1e-9) {
          return { ok: false, stage: "annular-slice-coincidence" };
        }
      }

      for (let index = 0; index < entries.length; index += 1) {
        const here = entries[index];
        const direction = sideFuture[here.id];
        const neighbor = entries[modInteger(index + direction, entries.length)];
        if (neighbor.id !== mate[here.id]) {
          return { ok: false, stage: "annular-nonadjacent" };
        }
        if (sideFuture[neighbor.id] !== -direction) {
          return { ok: false, stage: "annular-sides-not-facing" };
        }
        const pair = pairIndex[here.id];
        // Store phase as (positive-side circuit phase) minus (negative-side
        // circuit phase), independent of which endpoint is visited first.
        const delta = direction > 0
          ? modInteger(neighbor.phase - here.phase, q)
          : modInteger(here.phase - neighbor.phase, q);
        if (phaseByPair[pair] < 0) phaseByPair[pair] = delta;
        else if (phaseByPair[pair] !== delta) {
          return { ok: false, stage: "annular-phase-monodromy" };
        }
      }
    }
  }
  return { ok: true, phases: [...phaseByPair] };
}

const ANNULAR_STRIP_FAILURES = new Set([
  "annular-zero-sweep", "annular-sides-not-facing", "annular-nonmonotone",
  "annular-duplicate-slice-phase", "annular-slice-count",
  "annular-slice-coincidence", "annular-nonadjacent",
  "annular-phase-monodromy",
]);

const ANNULAR_REGION_FAILURES = new Set([
  "annular-boundaries-meet", "annular-complement-not-two",
  "annular-boundary-nonseparating", "annular-incompatible-sides",
  "annular-phase-nonmonotone", "annular-phase-boundaries-meet",
  "annular-phase-incompatible-sides", "annular-phase-range-too-small",
]);

function perfectMatchings(circuits, incompatible = null) {
  const groups = new Map();
  for (const circuit of circuits) {
    if (circuit.kind !== "essential") continue;
    const key = `${circuit.hx},${circuit.hy}`;
    if (!groups.has(key)) groups.set(key, { negative: [], positive: [] });
    groups.get(key)[circuit.canonicalSide > 0 ? "positive" : "negative"].push(circuit.id);
  }
  let all = [[]];
  for (const group of groups.values()) {
    if (group.negative.length !== group.positive.length) return [];
    const local = [];
    const used = new Uint8Array(group.positive.length);
    const pairs = [];
    function permute(index) {
      if (index === group.negative.length) {
        local.push(pairs.map((pair) => pair.slice()));
        return;
      }
      for (let j = 0; j < group.positive.length; j += 1) {
        if (used[j]) continue;
        const left = group.negative[index];
        const right = group.positive[j];
        if (incompatible !== null
            && incompatible.has(left < right ? `${left},${right}` : `${right},${left}`)) {
          continue;
        }
        used[j] = 1;
        pairs.push([left, right]);
        permute(index + 1);
        pairs.pop();
        used[j] = 0;
      }
    }
    permute(0);
    const combined = [];
    for (const prefix of all) for (const suffix of local) combined.push(prefix.concat(suffix));
    all = combined;
  }
  return all;
}

function annularPairIncompatibilities(arr, choice, trace) {
  const incompatible = new Set();
  for (const vertex of arr.vertices) {
    const pairs = choice[vertex.index] ? vertex.switchPairs : vertex.nonswitchPairs;
    const left = trace.circuitOf[pairs[0][0]];
    const right = trace.circuitOf[pairs[1][0]];
    if (left === right
        || trace.circuits[left].kind !== "essential"
        || trace.circuits[right].kind !== "essential") continue;
    incompatible.add(left < right ? `${left},${right}` : `${right},${left}`);
  }
  return incompatible;
}

class ParityDSU {
  constructor(n) {
    this.parent = new Int16Array(n);
    this.parity = new Int8Array(n);
    for (let i = 0; i < n; i += 1) this.parent[i] = i;
  }
  find(a) {
    if (this.parent[a] === a) return [a, 0];
    const [root, parity] = this.find(this.parent[a]);
    this.parity[a] ^= parity;
    this.parent[a] = root;
    return [root, this.parity[a]];
  }
  unite(a, b, wanted) {
    const fa = this.find(a);
    const fb = this.find(b);
    if (fa[0] === fb[0]) return (fa[1] ^ fb[1]) === wanted;
    this.parent[fb[0]] = fa[0];
    this.parity[fb[0]] = fa[1] ^ fb[1] ^ wanted;
    return true;
  }
}

function occupiedSectors(arr, vertex, port, color) {
  const index = vertex.ports.indexOf(port);
  const left = new Set([index, (index + 1) & 3]);
  const right = new Set([(index + 2) & 3, (index + 3) & 3]);
  const cooriented = arr.portSign[port] === 1 ? left : right;
  const opposite = arr.portSign[port] === 1 ? right : left;
  return color === 1 ? cooriented : opposite;
}

class RegionDSU {
  constructor(n) {
    this.parent = new Int16Array(n);
    this.size = new Int16Array(n);
    for (let i = 0; i < n; i += 1) {
      this.parent[i] = i;
      this.size[i] = 1;
    }
  }
  find(a) {
    let root = a;
    while (this.parent[root] !== root) root = this.parent[root];
    while (this.parent[a] !== a) {
      const next = this.parent[a];
      this.parent[a] = root;
      a = next;
    }
    return root;
  }
  unite(a, b) {
    a = this.find(a);
    b = this.find(b);
    if (a === b) return;
    if (this.size[a] < this.size[b]) [a, b] = [b, a];
    this.parent[b] = a;
    this.size[a] += this.size[b];
  }
}

function pointInPolygon(polygon, px, py) {
  let inside = false;
  for (let i = 0, j = polygon.length - 1; i < polygon.length; j = i, i += 1) {
    const [xi, yi] = polygon[i];
    const [xj, yj] = polygon[j];
    if ((yi > py) !== (yj > py)) {
      const hit = ((xj - xi) * (py - yi)) / (yj - yi) + xi;
      if (px < hit) inside = !inside;
    }
  }
  return inside;
}

function translatedDiskMultiplicity(arr, circuit, px = 0.37, py = 0.73) {
  const polygon = circuit.polygon;
  const xmin = Math.min(...polygon.map((point) => point[0]));
  const xmax = Math.max(...polygon.map((point) => point[0]));
  const ymin = Math.min(...polygon.map((point) => point[1]));
  const ymax = Math.max(...polygon.map((point) => point[1]));
  let count = 0;
  const kMin = Math.ceil((xmin - px) / arr.L);
  const kMax = Math.floor((xmax - px) / arr.L);
  const lMin = Math.ceil((ymin - py) / arr.L);
  const lMax = Math.floor((ymax - py) / arr.L);
  for (let k = kMin; k <= kMax; k += 1) {
    for (let l = lMin; l <= lMax; l += 1) {
      if (pointInPolygon(polygon, px + k * arr.L, py + l * arr.L)) count += 1;
    }
  }
  return count;
}

function bezout(a, b) {
  let oldR = Math.abs(a);
  let r = Math.abs(b);
  let oldS = 1;
  let s = 0;
  let oldT = 0;
  let t = 1;
  while (r !== 0) {
    const q = Math.floor(oldR / r);
    [oldR, r] = [r, oldR - q * r];
    [oldS, s] = [s, oldS - q * s];
    [oldT, t] = [t, oldT - q * t];
  }
  return [oldR, oldS * Math.sign(a || 1), oldT * Math.sign(b || 1)];
}

function phaseShift(hx, hy, phase) {
  const [divisor, xCoefficient, yCoefficient] = bezout(hx, hy);
  if (divisor !== 1) throw new Error("annular homology is not primitive");
  // hx * sy - hy * sx = phase.
  return [-yCoefficient * phase, xCoefficient * phase];
}

function visitLift(arr, visit) {
  const vertex = arr.vertices[visit[0]];
  const x = (visit[3] - vertex.x) / arr.L;
  const y = (visit[4] - vertex.y) / arr.L;
  if (!Number.isInteger(x) || !Number.isInteger(y)) {
    throw new Error("circuit visit does not have integral deck coordinates");
  }
  return [x, y];
}

function forbiddenAnnularPhases(arr, circuits, pair) {
  const [leftId, rightId] = pair;
  const left = circuits[leftId];
  const right = circuits[rightId];
  const visitsByVertex = new Map();
  for (const visit of right.visits) {
    if (!visitsByVertex.has(visit[0])) visitsByVertex.set(visit[0], []);
    visitsByVertex.get(visit[0]).push(visit);
  }
  const forbidden = new Set();
  for (const leftVisit of left.visits) {
    const rightVisits = visitsByVertex.get(leftVisit[0]) || [];
    const [leftX, leftY] = visitLift(arr, leftVisit);
    for (const rightVisit of rightVisits) {
      const [rightX, rightY] = visitLift(arr, rightVisit);
      const sx = leftX - rightX;
      const sy = leftY - rightY;
      forbidden.add(left.hx * sy - left.hy * sx);
    }
  }
  return forbidden;
}

function annularPhaseCandidates(circuit, windingRange) {
  // A transverse tau-slice meets q=|tau(h)| sheets of the quotient circuit.
  // Its finite sheet label is r mod q, while changing the relative lift by
  // one full transverse deck period changes det(h,s) by q.  Thus the full
  // integer phase is k=r+qz.  `windingRange` bounds z, not k.
  const q = Math.abs(2 * circuit.hx + circuit.hy);
  if (q === 0) return [];
  const phases = [];
  for (let z = -windingRange; z <= windingRange; z += 1) {
    for (let residue = 0; residue < q; residue += 1) {
      phases.push(residue + q * z);
    }
  }
  return phases;
}

// All slice levels used below are deliberately decimal rationals (the base
// point has two decimal places and the type-B probe uses one).  Converting
// them once to a fixed decimal denominator lets the longitudinal bounds and
// the annular multiplicity be computed with strict, exact inequalities.
const SLICE_SCALE = 1000000000n;

function bigGcd(a, b) {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b !== 0n) [a, b] = [b, a % b];
  return a;
}

function exactFraction(numerator, denominator = 1n) {
  if (denominator === 0n) throw new Error("zero exact denominator");
  if (denominator < 0n) {
    numerator = -numerator;
    denominator = -denominator;
  }
  const divisor = bigGcd(numerator, denominator);
  return { n: numerator / divisor, d: denominator / divisor };
}

function exactFromInteger(value) {
  if (!Number.isSafeInteger(value)) throw new Error("unsafe exact integer");
  return { n: BigInt(value), d: 1n };
}

function exactFromNumber(value) {
  if (!Number.isFinite(value)) throw new Error("nonfinite exact number");
  if (Number.isSafeInteger(value)) return exactFromInteger(value);
  return exactFraction(BigInt(Math.round(value * Number(SLICE_SCALE))), SLICE_SCALE);
}

function exactAdd(left, right) {
  return exactFraction(left.n * right.d + right.n * left.d, left.d * right.d);
}

function exactSubtract(left, right) {
  return exactFraction(left.n * right.d - right.n * left.d, left.d * right.d);
}

function exactMultiplyInteger(value, integer) {
  return exactFraction(value.n * BigInt(integer), value.d);
}

function exactDivide(left, right) {
  return exactFraction(left.n * right.d, left.d * right.n);
}

function exactCompare(left, right) {
  const difference = left.n * right.d - right.n * left.d;
  return difference < 0n ? -1 : difference > 0n ? 1 : 0;
}

function exactSign(value) {
  return value.n < 0n ? -1 : value.n > 0n ? 1 : 0;
}

function exactFloor(value) {
  let quotient = value.n / value.d;
  if (value.n < 0n && value.n % value.d !== 0n) quotient -= 1n;
  return quotient;
}

function exactCeil(value) {
  return -exactFloor({ n: -value.n, d: value.d });
}

function exactToNumber(value) {
  return Number(value.n) / Number(value.d);
}

function exactEqual(left, right) {
  return exactCompare(left, right) === 0;
}

function liftedLineHit(arr, circuit, shiftX, shiftY, level) {
  const rawPeriodX = arr.L * circuit.rawHx;
  const rawPeriodY = arr.L * circuit.rawHy;
  const rawPeriodTau = 2 * rawPeriodX + rawPeriodY;
  if (rawPeriodTau === 0) return null;
  const reverse = rawPeriodTau < 0;
  const periodX = reverse ? -rawPeriodX : rawPeriodX;
  const periodY = reverse ? -rawPeriodY : rawPeriodY;
  const periodTau = reverse ? -rawPeriodTau : rawPeriodTau;
  const exactLevel = exactFromNumber(level);
  const hits = [];
  for (const [rawAx, rawAy, rawBx, rawBy] of circuit.segments) {
    const firstX = (reverse ? rawBx : rawAx) + shiftX * arr.L;
    const firstY = (reverse ? rawBy : rawAy) + shiftY * arr.L;
    const secondX = (reverse ? rawAx : rawBx) + shiftX * arr.L;
    const secondY = (reverse ? rawAy : rawBy) + shiftY * arr.L;
    const ax = firstX;
    const ay = firstY;
    const bx = secondX;
    const by = secondY;
    const ta = 2 * ax + ay;
    const tb = 2 * bx + by;
    if (tb <= ta) return null;

    // ta+nT < level < tb+nT.  These strict bounds avoid both an arbitrary
    // search window and double-counting a vertex when a caller accidentally
    // chooses a nongeneric slice.
    const lower = exactDivide(
      exactSubtract(exactLevel, exactFromInteger(tb)),
      exactFromInteger(periodTau),
    );
    const upper = exactDivide(
      exactSubtract(exactLevel, exactFromInteger(ta)),
      exactFromInteger(periodTau),
    );
    const firstLongitudinal = exactFloor(lower) + 1n;
    const lastLongitudinal = exactCeil(upper) - 1n;
    for (let longitudinal = firstLongitudinal;
      longitudinal <= lastLongitudinal; longitudinal += 1n) {
      const shiftedTa = BigInt(ta) + longitudinal * BigInt(periodTau);
      const alpha = exactDivide(
        exactSubtract(exactLevel, exactFraction(shiftedTa)),
        exactFromInteger(tb - ta),
      );
      const baseX = BigInt(ax) + longitudinal * BigInt(periodX);
      const baseY = BigInt(ay) + longitudinal * BigInt(periodY);
      const xExact = exactAdd(
        exactFraction(baseX),
        exactMultiplyInteger(alpha, bx - ax),
      );
      const yExact = exactAdd(
        exactFraction(baseY),
        exactMultiplyInteger(alpha, by - ay),
      );
      hits.push({
        x: exactToNumber(xExact),
        y: exactToNumber(yExact),
        xExact,
        yExact,
        // Coorientation signs in traceCircuits use the original circuit
        // orientation, even when the monotone slice calculation reverses it.
        vx: rawBx - rawAx,
        vy: rawBy - rawAy,
      });
    }
  }
  return hits.length === 1 ? hits[0] : null;
}

function exactCrossSign(vx, vy, dx, dy) {
  return exactSign(exactSubtract(
    exactMultiplyInteger(dy, vx),
    exactMultiplyInteger(dx, vy),
  ));
}

function phaseAnnularRegion(arr, circuits, pair, color, phase, px = 0.37, py = 0.73) {
  const [leftId, rightId] = pair;
  const left = circuits[leftId];
  const right = circuits[rightId];
  const [shiftX, shiftY] = phaseShift(left.hx, left.hy, phase);
  const [stepX, stepY] = phaseShift(left.hx, left.hy, 1);
  const level = 2 * px + py;
  const q = Math.abs(2 * left.hx + left.hy);
  if (q === 0) return { ok: false, stage: "annular-phase-nonmonotone" };
  const exactPx = exactFromNumber(px);
  let multiplicity = 0n;

  // Every eye translate has label n=residue+q*a.  Increasing a by one
  // translates both slice hits by the same tau-zero deck period P, while
  // leaving their relative order and coorientation unchanged.  It therefore
  // suffices to inspect q residues and count all a in one exact integer
  // interval for each residue.
  for (let residue = 0; residue < q; residue += 1) {
    const translate = residue;
    const tx = translate * stepX;
    const ty = translate * stepY;
    const leftHit = liftedLineHit(arr, left, tx, ty, level);
    const rightHit = liftedLineHit(arr, right, tx + shiftX, ty + shiftY, level);
    const nextTx = (translate + q) * stepX;
    const nextTy = (translate + q) * stepY;
    const nextLeftHit = liftedLineHit(arr, left, nextTx, nextTy, level);
    const nextRightHit = liftedLineHit(
      arr, right, nextTx + shiftX, nextTy + shiftY, level,
    );
    if (leftHit === null || rightHit === null
        || nextLeftHit === null || nextRightHit === null) {
      return { ok: false, stage: "annular-phase-nonmonotone" };
    }
    const dx = exactSubtract(rightHit.xExact, leftHit.xExact);
    const dy = exactSubtract(rightHit.yExact, leftHit.yExact);
    if (exactSign(dx) === 0 && exactSign(dy) === 0) {
      return { ok: false, stage: "annular-phase-boundaries-meet" };
    }
    const leftToward = exactCrossSign(leftHit.vx, leftHit.vy, dx, dy);
    const rightToward = exactCrossSign(rightHit.vx, rightHit.vy,
      exactMultiplyInteger(dx, -1), exactMultiplyInteger(dy, -1));
    const leftInside = color * left.traceSide;
    const rightInside = color * right.traceSide;
    if (leftToward !== leftInside || rightToward !== rightInside) {
      return { ok: false, stage: "annular-phase-incompatible-sides" };
    }

    const leftPeriod = exactSubtract(nextLeftHit.xExact, leftHit.xExact);
    const rightPeriod = exactSubtract(nextRightHit.xExact, rightHit.xExact);
    if (!exactEqual(leftPeriod, rightPeriod) || exactSign(leftPeriod) === 0) {
      return { ok: false, stage: "annular-phase-nonmonotone" };
    }
    const lowX = exactCompare(leftHit.xExact, rightHit.xExact) < 0
      ? leftHit.xExact : rightHit.xExact;
    const highX = exactCompare(leftHit.xExact, rightHit.xExact) < 0
      ? rightHit.xExact : leftHit.xExact;
    const periodMagnitude = exactSign(leftPeriod) > 0
      ? leftPeriod : exactMultiplyInteger(leftPeriod, -1);

    // After replacing a by sign(P)*a, both boundaries move by +a|P|.
    // Thus (px-high)/|P| < a < (px-low)/|P|, with both inequalities strict.
    const lower = exactDivide(
      exactSubtract(exactPx, highX), periodMagnitude,
    );
    const upper = exactDivide(
      exactSubtract(exactPx, lowX), periodMagnitude,
    );
    const count = exactCeil(upper) - exactFloor(lower) - 1n;
    if (count > 0n) multiplicity += count;
  }
  if (multiplicity > BigInt(Number.MAX_SAFE_INTEGER)) {
    throw new Error("annular multiplicity exceeds the safe integer range");
  }
  return {
    ok: true,
    containsBase: multiplicity > 0n,
    multiplicity: Number(multiplicity),
    phase,
    shift: [shiftX, shiftY],
  };
}

/*
 * The four sectors at a crossing are numbered by the CCW interval from
 * port q to port q+1.  Across an arrangement edge, its two side sectors are
 * identified.  Locally, a selected smoothed boundary arc cuts the vertex
 * disk into the two cyclic intervals between its endpoint ports; when there
 * is no selected arc, all four sectors communicate.
 *
 * Applying this to the two proposed essential boundary circuits computes
 * the complement on the torus without choosing a deck phase.  A genuine
 * embedded annulus has exactly two components, and the color/coorientation
 * side of both boundaries must select the same component.
 */
function annularRegion(arr, circuits, pair, color, baseSector) {
  const selected = new Set(pair);
  const arcsAtVertex = Array.from({ length: arr.vertices.length }, () => []);
  for (const circuitId of pair) {
    for (const visit of circuits[circuitId].visits) {
      arcsAtVertex[visit[0]].push([visit[1], visit[2], circuitId]);
    }
  }
  if (arcsAtVertex.some((arcs) => arcs.length > 1)) {
    return { ok: false, stage: "annular-boundaries-meet" };
  }

  const regions = new RegionDSU(4 * arr.vertices.length);
  const sectorNode = (vertex, sector) => 4 * vertex + ((sector + 4) & 3);
  for (const edge of arr.edges) {
    const va = arr.portVertex[edge.a];
    const vb = arr.portVertex[edge.b];
    const ia = arr.vertices[va].ports.indexOf(edge.a);
    const ib = arr.vertices[vb].ports.indexOf(edge.b);
    regions.unite(sectorNode(va, ia), sectorNode(vb, ib - 1));
    regions.unite(sectorNode(va, ia - 1), sectorNode(vb, ib));
  }
  function cyclicInterval(start, end) {
    const answer = [];
    for (let q = start; q !== end; q = (q + 1) & 3) answer.push(q);
    return answer;
  }
  function joinGroup(vertex, group) {
    for (let i = 1; i < group.length; i += 1) {
      regions.unite(sectorNode(vertex, group[0]), sectorNode(vertex, group[i]));
    }
  }
  for (const vertex of arr.vertices) {
    const arcs = arcsAtVertex[vertex.index];
    if (arcs.length === 0) {
      joinGroup(vertex.index, [0, 1, 2, 3]);
      continue;
    }
    const ia = vertex.ports.indexOf(arcs[0][0]);
    const ib = vertex.ports.indexOf(arcs[0][1]);
    joinGroup(vertex.index, cyclicInterval(ia, ib));
    joinGroup(vertex.index, cyclicInterval(ib, ia));
  }
  const roots = new Set();
  for (let node = 0; node < 4 * arr.vertices.length; node += 1) roots.add(regions.find(node));
  if (roots.size !== 2) return { ok: false, stage: "annular-complement-not-two" };

  let insideRoot = -1;
  let outsideRoot = -1;
  for (const circuitId of pair) {
    const circuit = circuits[circuitId];
    const visit = circuit.visits[0];
    const vertex = arr.vertices[visit[0]];
    const outgoingIndex = vertex.ports.indexOf(visit[2]);
    const interiorSide = color * circuit.traceSide;
    const inside = sectorNode(vertex.index, interiorSide > 0 ? outgoingIndex : outgoingIndex - 1);
    const outside = sectorNode(vertex.index, interiorSide > 0 ? outgoingIndex - 1 : outgoingIndex);
    const thisInside = regions.find(inside);
    const thisOutside = regions.find(outside);
    if (thisInside === thisOutside) return { ok: false, stage: "annular-boundary-nonseparating" };
    if (insideRoot >= 0 && insideRoot !== thisInside) {
      return { ok: false, stage: "annular-incompatible-sides" };
    }
    if (outsideRoot >= 0 && outsideRoot !== thisOutside) {
      return { ok: false, stage: "annular-incompatible-sides" };
    }
    insideRoot = thisInside;
    outsideRoot = thisOutside;
  }
  return {
    ok: true,
    insideRoot,
    outsideRoot,
    containsBase: regions.find(baseSector) === insideRoot,
    selected,
  };
}

function arrangementBaseSector(arr, px = 0.37, py = 0.73) {
  const seen = new Uint8Array(arr.portCount);
  for (let seed = 0; seed < arr.portCount; seed += 1) {
    if (seen[seed]) continue;
    let current = seed;
    let liftX = 0;
    let liftY = 0;
    const polygon = [];
    const sectors = [];
    for (let guard = 0; ; guard += 1) {
      if (guard > arr.portCount) throw new Error("face trace overflow");
      if (guard > 0 && current === seed) break;
      if (seen[current]) throw new Error("face trace has a tail");
      seen[current] = 1;
      const vertexId = arr.portVertex[current];
      const vertex = arr.vertices[vertexId];
      const portIndex = vertex.ports.indexOf(current);
      polygon.push([vertex.x + liftX * arr.L, vertex.y + liftY * arr.L]);
      sectors.push(4 * vertexId + portIndex);
      liftX += arr.portShiftX[current];
      liftY += arr.portShiftY[current];
      const incoming = arr.portTwin[current];
      const target = arr.vertices[arr.portVertex[incoming]];
      const incomingIndex = target.ports.indexOf(incoming);
      current = target.ports[(incomingIndex + 3) & 3];
    }
    if (liftX !== 0 || liftY !== 0) throw new Error("essential arrangement face");
    const xmin = Math.min(...polygon.map((point) => point[0]));
    const xmax = Math.max(...polygon.map((point) => point[0]));
    const ymin = Math.min(...polygon.map((point) => point[1]));
    const ymax = Math.max(...polygon.map((point) => point[1]));
    const kMin = Math.ceil((xmin - px) / arr.L);
    const kMax = Math.floor((xmax - px) / arr.L);
    const lMin = Math.ceil((ymin - py) / arr.L);
    const lMax = Math.floor((ymax - py) / arr.L);
    for (let k = kMin; k <= kMax; k += 1) {
      for (let l = lMin; l <= lMax; l += 1) {
        if (pointInPolygon(polygon, px + k * arr.L, py + l * arr.L)) return sectors[0];
      }
    }
  }
  throw new Error("base point is not in an arrangement face");
}

function strictSliceHits(segments, level, shiftX = 0, shiftY = 0) {
  const hits = [];
  for (const [rawAx, rawAy, rawBx, rawBy] of segments) {
    const ax = rawAx + shiftX;
    const ay = rawAy + shiftY;
    const bx = rawBx + shiftX;
    const by = rawBy + shiftY;
    const ta = 2 * ax + ay;
    const tb = 2 * bx + by;
    if (!((ta < level && level < tb) || (tb < level && level < ta))) continue;
    const alpha = (level - ta) / (tb - ta);
    hits.push({ x: ax + alpha * (bx - ax), y: ay + alpha * (by - ay) });
  }
  return hits;
}

/* All deck translates of an essential circuit at one tau-slice, reduced to
 * the few representatives nearest referenceX.  Writing a deck vector as
 * (m,n), its tau shift is L(2m+n) and its x shift is Lm. */
function periodicSliceHits(arr, circuit, level, referenceX) {
  const hits = [];
  for (const [ax, ay, bx, by] of circuit.segments) {
    const ta = 2 * ax + ay;
    const tb = 2 * bx + by;
    const low = Math.min(ta, tb);
    const high = Math.max(ta, tb);
    const rMin = Math.floor((level - high) / arr.L) + 1;
    const rMax = Math.ceil((level - low) / arr.L) - 1;
    for (let r = rMin; r <= rMax; r += 1) {
      const target = level - r * arr.L;
      const alpha = (target - ta) / (tb - ta);
      if (!(alpha > 0 && alpha < 1)) continue;
      const baseX = ax + alpha * (bx - ax);
      const center = Math.round((referenceX - baseX) / arr.L);
      for (let m = center - 2; m <= center + 2; m += 1) {
        hits.push(baseX + m * arr.L);
      }
    }
  }
  hits.sort((a, b) => a - b);
  return hits.filter((value, index) => index === 0 || Math.abs(value - hits[index - 1]) > 1e-8);
}

function matchingVisit(circuit, vertex, pair) {
  const key = [...pair].sort((a, b) => a - b);
  const visits = circuit.visits.filter((visit) => {
    if (visit[0] !== vertex.index) return false;
    const local = [visit[1], visit[2]].sort((a, b) => a - b);
    return local[0] === key[0] && local[1] === key[1];
  });
  return visits.length === 1 ? visits[0] : null;
}

/*
 * Full type-B normal order, including disk/annulus mixtures.  On each of the
 * two nearby tau-slices, locate the incident boundary of both eyes and its
 * companion boundary.  For an annular eye the color-selected side chooses
 * the nearest translate of the paired essential circuit.  The active
 * type-B pictures in main.tex are exactly the cases where the two incident
 * paths and the two companion paths occur in the same transverse order.
 */
function typeBCompanionRelation(
  arr, circuits, vertex, circuitIds, eyeIds, eyeCircuits, interiorSide,
  annularPhaseByEye,
) {
  const entries = [];
  for (let index = 0; index < 2; index += 1) {
    const circuitId = circuitIds[index];
    const circuit = circuits[circuitId];
    const visit = matchingVisit(circuit, vertex, vertex.switchPairs[index]);
    if (visit === null) return "inconsistent";
    entries.push({ circuitId, circuit, visit, eyeId: eyeIds[index] });
  }
  const commonX = vertex.x;
  const commonY = vertex.y;
  const switchLevel = 2 * commonX + commonY;
  const relations = new Set();
  for (const side of [-1, 1]) {
    const level = switchLevel + side * 0.1;
    const paths = [];
    for (const entry of entries) {
      const incoming = entry.visit[1];
      const outgoing = entry.visit[2];
      const candidates = [
        {
          port: incoming,
          vx: arr.portVectorX[incoming],
          vy: arr.portVectorY[incoming],
          traceVx: -arr.portVectorX[incoming],
          traceVy: -arr.portVectorY[incoming],
        },
        {
          port: outgoing,
          vx: arr.portVectorX[outgoing],
          vy: arr.portVectorY[outgoing],
          traceVx: arr.portVectorX[outgoing],
          traceVy: arr.portVectorY[outgoing],
        },
      ].filter((ray) => Math.sign(2 * ray.vx + ray.vy) === side);
      if (candidates.length !== 1) return "inconsistent";
      const ray = candidates[0];
      const alpha = (level - switchLevel) / (2 * ray.vx + ray.vy);
      const nearX = commonX + alpha * ray.vx;
      let companionX;
      const boundaries = eyeCircuits[entry.eyeId];
      if (boundaries.length === 1) {
        const shiftX = commonX - entry.visit[3];
        const shiftY = commonY - entry.visit[4];
        const hits = strictSliceHits(entry.circuit.segments, level, shiftX, shiftY)
          .sort((left, right) => Math.abs(left.x - nearX) - Math.abs(right.x - nearX));
        if (hits.length !== 2 || Math.abs(hits[0].x - nearX) > 1e-7) return "inconsistent";
        companionX = hits[1].x;
      } else {
        const mateId = boundaries[0] === entry.circuitId ? boundaries[1] : boundaries[0];
        const phaseData = annularPhaseByEye === null
          ? null
          : annularPhaseByEye[entry.eyeId];
        if (phaseData !== null && phaseData !== undefined) {
          const baseShiftX = (commonX - entry.visit[3]) / arr.L;
          const baseShiftY = (commonY - entry.visit[4]) / arr.L;
          if (!Number.isInteger(baseShiftX) || !Number.isInteger(baseShiftY)) {
            return "inconsistent";
          }
          let mateShiftX;
          let mateShiftY;
          if (entry.circuitId === phaseData.left) {
            mateShiftX = baseShiftX + phaseData.shift[0];
            mateShiftY = baseShiftY + phaseData.shift[1];
          } else if (entry.circuitId === phaseData.right) {
            mateShiftX = baseShiftX - phaseData.shift[0];
            mateShiftY = baseShiftY - phaseData.shift[1];
          } else return "inconsistent";
          const hit = liftedLineHit(arr, circuits[mateId], mateShiftX, mateShiftY, level);
          if (hit === null) return "inconsistent";
          companionX = hit.x;
        } else {
          const dtTrace = 2 * ray.traceVx + ray.traceVy;
          const direction = -interiorSide[entry.circuitId] * Math.sign(dtTrace);
          const hits = periodicSliceHits(arr, circuits[mateId], level, nearX)
            .filter((x) => direction * (x - nearX) > 1e-8)
            .sort((left, right) => direction * (left - right));
          if (hits.length === 0) return "inconsistent";
          companionX = hits[0];
        }
      }
      paths.push({ nearX, companionX });
    }
    const nearDifference = paths[1].nearX - paths[0].nearX;
    const farDifference = paths[1].companionX - paths[0].companionX;
    if (Math.abs(nearDifference) < 1e-8 || Math.abs(farDifference) < 1e-8) {
      return "inconsistent";
    }
    relations.add(nearDifference * farDifference > 0 ? "same" : "opposite");
  }
  return relations.size === 1 ? [...relations][0] : "inconsistent";
}

function ribbonTopology(arr, choice, circuits, interiorSide, switchKind, chi) {
  const dartsByBoundary = circuits.map(() => []);
  const dartsAtVertex = Array.from({ length: arr.vertices.length }, () => []);
  let dartCount = 0;
  for (const circuit of circuits) {
    let switchVertices = circuit.visits
      .filter((visit) => choice[visit[0]] !== 0)
      .map((visit) => visit[0]);
    // Normalize every eye-boundary component so that the eye lies to its
    // left.  This is area orientation for disks and color-selected strip
    // orientation for essential annular boundaries.
    if (interiorSide[circuit.id] < 0) switchVertices = switchVertices.reverse();
    for (const vertex of switchVertices) {
      dartsByBoundary[circuit.id].push(dartCount);
      dartsAtVertex[vertex].push(dartCount);
      dartCount += 1;
    }
  }
  const untouched = dartsByBoundary.filter((darts) => darts.length === 0).length;
  const adjacency = Array.from({ length: 2 * dartCount }, () => []);
  function join(left, right) {
    adjacency[left].push(right);
    adjacency[right].push(left);
  }
  for (const darts of dartsByBoundary) {
    if (darts.length === 0) continue;
    for (let i = 0; i < darts.length; i += 1) {
      join(2 * darts[i] + 1, 2 * darts[(i + 1) % darts.length]);
    }
  }
  for (let vertex = 0; vertex < arr.vertices.length; vertex += 1) {
    if (!choice[vertex]) continue;
    const darts = dartsAtVertex[vertex];
    if (darts.length !== 2) return { ok: false, stage: "bad-switch-dart-count" };
    const [left, right] = darts;
    if (switchKind[vertex] === 2) {
      join(2 * left, 2 * right + 1);
      join(2 * left + 1, 2 * right);
    } else {
      join(2 * left, 2 * right);
      join(2 * left + 1, 2 * right + 1);
    }
  }
  if (adjacency.some((neighbors) => neighbors.length !== 2)) {
    return { ok: false, stage: "bad-ribbon-flag-degree" };
  }
  const seen = new Uint8Array(adjacency.length);
  const cycleLengths = [];
  for (let seed = 0; seed < adjacency.length; seed += 1) {
    if (seen[seed]) continue;
    let previous = -1;
    let current = seed;
    let length = 0;
    while (!seen[current]) {
      seen[current] = 1;
      length += 1;
      const next = adjacency[current][0] === previous
        ? adjacency[current][1]
        : adjacency[current][0];
      previous = current;
      current = next;
    }
    cycleLengths.push(length);
  }
  cycleLengths.sort((a, b) => a - b);
  const boundaryComponents = untouched + cycleLengths.length;
  const numerator = 2 - boundaryComponents - chi;
  const genus = numerator >= 0 && numerator % 2 === 0 ? numerator / 2 : null;
  if (boundaryComponents !== arr.boundaryComponents || genus !== arr.targetGenus) {
    return {
      ok: false,
      stage: "ribbon-topology",
      boundaryComponents,
      genus,
      cycleLengths,
      untouched,
    };
  }
  return { ok: true, boundaryComponents, genus, cycleLengths, untouched };
}

function localAudit(arr, choice, trace, matching, colorMask, options = {}, phaseVector = null) {
  const requireNormalB = options.requireNormalB !== false;
  const requireBalance = options.requireBalance !== false;
  // The eye incidence, A/B/V/N kind, connectivity, orientability and ribbon
  // topology depend on the matching and colors, but not on the relative deck
  // phases of annular eyes.  The exact no-cutoff phase solver uses this mode
  // once per matching/color before compiling the remaining linear phase
  // constraints.  Default callers still follow the original full audit.
  const structuralOnly = options.structuralOnly === true;
  const phaseMode = phaseVector !== null || structuralOnly;
  const allowQuotientSelfSwitch = options.selfSwitchMode === "quotient" || phaseMode;
  const { circuits, circuitOf, turnAt } = trace;
  const circuitToEye = new Int16Array(circuits.length);
  const interiorSide = new Int8Array(circuits.length);
  circuitToEye.fill(-1);
  const eyeColors = [];
  const eyeKinds = [];
  const eyeCircuits = [];
  const annularPhaseByEye = [];
  let eyes = 0;
  let disks = 0;
  for (const circuit of circuits) {
    if (circuit.kind !== "disk") continue;
    circuitToEye[circuit.id] = eyes;
    interiorSide[circuit.id] = circuit.interiorSide;
    eyeColors.push(circuit.color);
    eyeKinds.push("disk");
    eyeCircuits.push([circuit.id]);
    annularPhaseByEye.push(null);
    eyes += 1;
    disks += 1;
  }
  for (let i = 0; i < matching.length; i += 1) {
    const color = ((colorMask >> BigInt(i)) & 1n) !== 0n ? 1 : -1;
    const [left, right] = matching[i];
    circuitToEye[left] = eyes;
    circuitToEye[right] = eyes;
    interiorSide[left] = color * circuits[left].traceSide;
    interiorSide[right] = color * circuits[right].traceSide;
    eyeColors.push(color);
    eyeKinds.push("annular");
    eyeCircuits.push([left, right]);
    const phase = phaseMode ? (structuralOnly ? 0 : phaseVector[i]) : null;
    annularPhaseByEye.push(phaseMode ? {
      left,
      right,
      phase,
      shift: phaseShift(circuits[left].hx, circuits[left].hy, phase),
    } : null);
    eyes += 1;
  }

  // The two proposed boundary circuits of one embedded annulus must be
  // disjoint in the resolved arrangement, including at switches.
  if (!phaseMode) {
    for (const vertex of arr.vertices) {
      const pairs = choice[vertex.index] ? vertex.switchPairs : vertex.nonswitchPairs;
      const left = circuitToEye[circuitOf[pairs[0][0]]];
      const right = circuitToEye[circuitOf[pairs[1][0]]];
      if (left === right && eyeKinds[left] === "annular") {
        return { ok: false, stage: "annular-boundaries-meet" };
      }
    }
  }

  // The no-intervening-boundary sweep can either be retained as a diagnostic
  // (eyes may overlap/nest) or imposed as the stricter adjacent-strip model.
  const stripAudit = options.stripAudit
    || auditAnnularStrips(arr, trace, matching, colorMask);
  if (options.requireAdjacentStrips === true && !stripAudit.ok) {
    return { ...stripAudit, stripOk: false, regionOk: false };
  }

  let baseBalance = 0;
  if (!structuralOnly) {
    for (const circuit of circuits) {
      if (circuit.kind === "disk") {
        baseBalance += circuit.color * translatedDiskMultiplicity(arr, circuit);
      }
    }
  }
  // Structural-only and exact-phase audits do not use the legacy fixed
  // arrangement base sector.  Computing it anyway can create a spurious
  // failure when that legacy probe happens to lie on a generic input line.
  let baseSector = null;
  if (!structuralOnly && !phaseMode && matching.length > 0) {
    baseSector = arr.baseSector === undefined
      ? (arr.baseSector = arrangementBaseSector(arr))
      : arr.baseSector;
  }
  for (let i = 0; !structuralOnly && i < matching.length; i += 1) {
    const color = ((colorMask >> BigInt(i)) & 1n) !== 0n ? 1 : -1;
    const region = phaseMode
      ? (options.phaseRegions === undefined
        ? phaseAnnularRegion(arr, circuits, matching[i], color, phaseVector[i])
        : options.phaseRegions[i])
      : annularRegion(arr, circuits, matching[i], color, baseSector);
    if (!region.ok) {
      return { ...region, stripOk: stripAudit.ok, regionOk: false };
    }
    if (phaseMode) baseBalance += color * region.multiplicity;
    else if (region.containsBase) baseBalance += color;
  }
  if (requireBalance && baseBalance !== 0) {
    return {
      ok: false,
      stage: "base-unbalanced",
      baseBalance,
      stripOk: stripAudit.ok,
      regionOk: true,
    };
  }

  const topology = new ParityDSU(eyes);
  const connectivity = new ParityDSU(eyes);
  const localTypes = { A: 0, B: 0, V: 0, N: 0 };
  const typeBCompanionOrders = { same: 0, opposite: 0, inconsistent: 0 };
  const deferredTypeB = [];
  const switchKind = new Int8Array(arr.vertices.length);
  let selfSwitches = 0;
  const selfSwitchDecks = [];
  for (const vertex of arr.vertices) {
    if (!choice[vertex.index]) {
      localTypes.N += 1;
      const pairs = vertex.nonswitchPairs;
      const ids = pairs.map((pair) => circuitToEye[circuitOf[pair[0]]]);
      const occupied = [...occupiedSectors(arr, vertex, pairs[0][0], eyeColors[ids[0]])]
        .filter((sector) => occupiedSectors(arr, vertex, pairs[1][0], eyeColors[ids[1]]).has(sector));
      if (occupied.length !== 1) return { ok: false, stage: "invalid-nonswitch" };
      continue;
    }
    const pairs = vertex.switchPairs;
    const circuitIds = pairs.map((pair) => circuitOf[pair[0]]);
    const ids = circuitIds.map((circuit) => circuitToEye[circuit]);
    if (ids[0] === ids[1]) {
      if (!allowQuotientSelfSwitch) return { ok: false, stage: "one-eye-at-switch" };
      if (circuitIds[0] !== circuitIds[1]) {
        if (!phaseMode || eyeKinds[ids[0]] !== "annular") {
          return { ok: false, stage: "self-switch-two-boundaries" };
        }
        selfSwitches += 1;
        selfSwitchDecks.push([vertex.index, "annular-phase", annularPhaseByEye[ids[0]].phase]);
      } else {
        const selfCircuit = circuits[circuitIds[0]];
        const visits = pairs.map((pair) => matchingVisit(selfCircuit, vertex, pair));
        if (visits.some((visit) => visit === null)) {
          return { ok: false, stage: "self-switch-visit-inconsistent" };
        }
        const deckX = (visits[1][3] - visits[0][3]) / arr.L;
        const deckY = (visits[1][4] - visits[0][4]) / arr.L;
        if (!Number.isInteger(deckX) || !Number.isInteger(deckY)) {
          return { ok: false, stage: "self-switch-nondeck-offset" };
        }
        // A disk lift has no longitudinal period, so only the zero deck
        // offset is the same lift.  An essential boundary is periodic along
        // h; there its invariant lift label is det(h,deck), and any multiple
        // of h is the same lifted boundary.
        const deckPhase = selfCircuit.kind === "essential"
          ? selfCircuit.hx * deckY - selfCircuit.hy * deckX
          : null;
        const sameLift = selfCircuit.kind === "disk"
          ? deckX === 0 && deckY === 0
          : deckPhase === 0;
        if (sameLift) {
          return { ok: false, stage: "same-lift-self-switch" };
        }
        selfSwitches += 1;
        selfSwitchDecks.push([vertex.index, deckX, deckY, deckPhase]);
      }
    }
    const convex = pairs.map((pair, index) => (
      turnAt[pair[0]] * interiorSide[circuitIds[index]] > 0
    ));
    let kind;
    if (vertex.orientation === "V") {
      if (eyeColors[ids[0]] === eyeColors[ids[1]] || !convex[0] || !convex[1]) {
        return { ok: false, stage: "invalid-V" };
      }
      kind = "V";
    } else if (eyeColors[ids[0]] !== eyeColors[ids[1]]) {
      if (!convex[0] || !convex[1]) return { ok: false, stage: "invalid-A" };
      kind = "A";
    } else {
      if (convex[0] === convex[1]) return { ok: false, stage: "invalid-B" };
      kind = "B";
      if (structuralOnly) {
        deferredTypeB.push({
          vertex,
          circuitIds: circuitIds.slice(),
          ids: ids.slice(),
        });
      } else {
        const relation = typeBCompanionRelation(
          arr, circuits, vertex, circuitIds, ids, eyeCircuits, interiorSide,
          phaseMode ? annularPhaseByEye : null,
        );
        typeBCompanionOrders[relation] += 1;
        if (requireNormalB && relation !== "same") {
          return {
            ok: false,
            stage: "non-normal-B",
            relation,
            stripOk: stripAudit.ok,
            regionOk: true,
            typeBCompanionOrders,
          };
        }
      }
    }
    localTypes[kind] += 1;
    switchKind[vertex.index] = kind === "B" ? 2 : 1;
    connectivity.unite(ids[0], ids[1], 0);
    if (!topology.unite(ids[0], ids[1], kind === "B" ? 0 : 1)) {
      return { ok: false, stage: "nonorientable" };
    }
  }
  const roots = new Set();
  for (let eye = 0; eye < eyes; eye += 1) roots.add(connectivity.find(eye)[0]);
  if (roots.size !== 1) return { ok: false, stage: "disconnected" };
  const switches = localTypes.A + localTypes.B + localTypes.V;
  const chi = disks - switches;
  if (chi !== arr.targetChi) return { ok: false, stage: "wrong-euler" };
  const ribbon = ribbonTopology(arr, choice, circuits, interiorSide, switchKind, chi);
  if (!ribbon.ok) return ribbon;
  if (structuralOnly) {
    return {
      ok: true,
      stage: "phase-independent-structural",
      eyes,
      disks,
      annuli: matching.length,
      switches,
      chi,
      localTypes,
      circuitToEye,
      interiorSide,
      eyeColors,
      eyeKinds,
      eyeCircuits,
      deferredTypeB,
      ribbon,
      selfSwitches,
      selfSwitchDecks,
    };
  }
  return {
    ok: true,
    stage: "connected-local",
    eyes,
    disks,
    annuli: matching.length,
    switches,
    chi,
    baseBalance,
    localTypes,
    typeBCompanionOrders,
    matching: matching.map((pair) => pair.slice()),
    colors: matching.map((_, i) => (((colorMask >> BigInt(i)) & 1n) !== 0n ? "B" : "R")),
    stripOk: stripAudit.ok,
    stripStage: stripAudit.ok ? "adjacent-strip" : stripAudit.stage,
    regionOk: true,
    annularPhases: phaseMode ? phaseVector.slice() : (stripAudit.ok ? stripAudit.phases : null),
    annularPhaseShifts: phaseMode
      ? matching.map((_, index) => annularPhaseByEye[disks + index].shift.slice())
      : null,
    ribbonBoundaryComponents: ribbon.boundaryComponents,
    ribbonGenus: ribbon.genus,
    ribbonCycleLengths: ribbon.cycleLengths,
    selfSwitches,
    selfSwitchDecks,
  };
}

function auditMask(arr, choice, options = {}) {
  const trace = traceCircuits(arr, choice);
  if (trace === null) return { ok: false, stage: "rich-trace" };
  const diskBase = trace.circuits.reduce((sum, circuit) => (
    circuit.kind === "disk"
      ? sum + circuit.color * translatedDiskMultiplicity(arr, circuit)
      : sum
  ), 0);
  const essential = trace.circuits.filter((circuit) => circuit.kind === "essential");
  if (essential.length === 0) return { ok: false, stage: "all-disk" };
  const phaseRange = options.annularPhaseRange;
  const phaseMode = Number.isInteger(phaseRange) && phaseRange >= 0;
  const matchings = perfectMatchings(
    trace.circuits,
    phaseMode ? null : annularPairIncompatibilities(arr, choice, trace),
  );
  if (matchings.length === 0) return { ok: false, stage: "unpairable" };
  const totals = {
    ok: true,
    stage: "pairable",
    diskEyes: trace.circuits.length - essential.length,
    essentialCircuits: essential.length,
    annularEyes: essential.length / 2,
    matchingConfigurations: 0,
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
  function phaseProducts(lists, callback, prefix = [], index = 0) {
    if (index === lists.length) {
      callback(prefix.slice());
      return;
    }
    for (const phase of lists[index]) {
      prefix.push(phase);
      phaseProducts(lists, callback, prefix, index + 1);
      prefix.pop();
    }
  }
  function recordResult(result, stripAudit, phaseModeForResult) {
    if (result.stage !== "annular-boundaries-meet"
        && result.stage !== "annular-phase-boundaries-meet") {
      totals.boundaryDisjoint += 1;
    }
    if (!phaseModeForResult
        && result.stage !== "annular-boundaries-meet" && stripAudit.ok) {
      totals.stripValid += 1;
    }
    if (!ANNULAR_REGION_FAILURES.has(result.stage)
        && (!options.requireAdjacentStrips || stripAudit.ok)) {
      totals.embeddedValid += 1;
    }
    if (!result.ok) {
      totals.rejected[result.stage] = (totals.rejected[result.stage] || 0) + 1;
      if (result.stage === "disconnected") totals.locallyValid += 1;
      return;
    }
    totals.locallyValid += 1;
    totals.connected += 1;
    const t = result.localTypes;
    const key = `D${result.disks}-R${result.annuli}-A${t.A}-B${t.B}-V${t.V}`
      + `-S${result.selfSwitches}`;
    totals.connectedDistribution[key] = (totals.connectedDistribution[key] || 0) + 1;
    if (options.collectAnswers === true || totals.witnesses.length < 4) {
      totals.witnesses.push(result);
    }
  }
  for (const matching of matchings) {
    const phaseLists = phaseMode ? matching.map((pair) => {
      const forbidden = forbiddenAnnularPhases(arr, trace.circuits, pair);
      const allowed = [];
      for (const phase of annularPhaseCandidates(trace.circuits[pair[0]], phaseRange)) {
        if (!forbidden.has(phase)) allowed.push(phase);
      }
      return allowed;
    }) : [null];
    if (phaseMode && phaseLists.some((list) => list.length === 0)) continue;
    if (phaseMode) {
      // Region geometry depends only on one paired eye, its color and its
      // phase.  Precompute it once per candidate instead of recomputing it
      // for every Cartesian product of the other annular eyes.
      const entryLists = matching.map((pair, index) => {
        const entries = [];
        for (const phase of phaseLists[index]) {
          for (let bit = 0; bit < 2; bit += 1) {
            const color = bit === 0 ? -1 : 1;
            const region = phaseAnnularRegion(
              arr, trace.circuits, pair, color, phase,
            );
            entries.push({
              phase,
              bit,
              region,
              contribution: region.ok ? color * region.multiplicity : null,
            });
          }
        }
        return entries;
      });
      totals.matchingConfigurations += phaseLists.reduce(
        (product, list) => product * list.length, 1,
      );
      totals.colorConfigurations += entryLists.reduce(
        (product, list) => product * list.length, 1,
      );

      const validLists = entryLists.map((list) => list.filter((entry) => entry.region.ok));
      // Account exactly for combinations rejected by the first bad annular
      // eye, matching the order used by localAudit, without enumerating them.
      const suffixRaw = new Array(entryLists.length + 1).fill(1);
      for (let i = entryLists.length - 1; i >= 0; i -= 1) {
        suffixRaw[i] = suffixRaw[i + 1] * entryLists[i].length;
      }
      let validPrefix = 1;
      for (let i = 0; i < entryLists.length; i += 1) {
        const invalidByStage = Object.create(null);
        for (const entry of entryLists[i]) {
          if (entry.region.ok) continue;
          invalidByStage[entry.region.stage] = (invalidByStage[entry.region.stage] || 0) + 1;
        }
        for (const [stage, count] of Object.entries(invalidByStage)) {
          const amount = validPrefix * count * suffixRaw[i + 1];
          totals.rejected[stage] = (totals.rejected[stage] || 0) + amount;
          if (stage !== "annular-phase-boundaries-meet") {
            totals.boundaryDisjoint += amount;
          }
        }
        validPrefix *= validLists[i].length;
      }
      if (validLists.some((list) => list.length === 0)) continue;

      let productLists = validLists;
      let enumerateProducts = phaseProducts;
      if (options.requireBalance !== false && options.phaseBalancePruning !== false) {
        // Exact suffix dynamic programming in the integer base multiplicity.
        // suffixCounts[i][s] is the number of tuples from eyes i..end whose
        // total annular contribution is s.  It both counts all unbalanced
        // tuples for statistics and prunes them before the costly local/type-B
        // audit.
        const suffixCounts = new Array(validLists.length + 1);
        suffixCounts[validLists.length] = new Map([[0, 1]]);
        for (let i = validLists.length - 1; i >= 0; i -= 1) {
          const counts = new Map();
          for (const entry of validLists[i]) {
            for (const [tailSum, tailCount] of suffixCounts[i + 1]) {
              const sum = entry.contribution + tailSum;
              counts.set(sum, (counts.get(sum) || 0) + tailCount);
            }
          }
          suffixCounts[i] = counts;
        }
        const target = -diskBase;
        const validCount = validLists.reduce((product, list) => product * list.length, 1);
        const balancedCount = suffixCounts[0].get(target) || 0;
        const unbalancedCount = validCount - balancedCount;
        if (unbalancedCount > 0) {
          totals.rejected["base-unbalanced"] =
            (totals.rejected["base-unbalanced"] || 0) + unbalancedCount;
          totals.boundaryDisjoint += unbalancedCount;
          if (!options.requireAdjacentStrips) totals.embeddedValid += unbalancedCount;
        }
        enumerateProducts = (lists, callback) => {
          const prefix = [];
          function visit(index, sum) {
            if (index === lists.length) {
              callback(prefix.slice());
              return;
            }
            for (const entry of lists[index]) {
              const nextSum = sum + entry.contribution;
              if (!suffixCounts[index + 1].has(target - nextSum)) continue;
              prefix.push(entry);
              visit(index + 1, nextSum);
              prefix.pop();
            }
          }
          visit(0, 0);
        };
        if (balancedCount === 0) continue;
      }
      enumerateProducts(productLists, (entries) => {
        const phaseVector = entries.map((entry) => entry.phase);
        let colorMask = 0n;
        for (let i = 0; i < entries.length; i += 1) {
          if (entries[i].bit !== 0) colorMask |= 1n << BigInt(i);
        }
        const stripAudit = {
          ok: false,
          stage: "explicit-annular-phase",
          phases: phaseVector,
        };
        const result = localAudit(arr, choice, trace, matching, colorMask, {
          ...options,
          stripAudit,
          phaseRegions: entries.map((entry) => entry.region),
        }, phaseVector);
        recordResult(result, stripAudit, true);
      });
      continue;
    }

    totals.matchingConfigurations += 1;
    {
      const colorLimit = 1n << BigInt(matching.length);
      for (let colorMask = 0n; colorMask < colorLimit; colorMask += 1n) {
        totals.colorConfigurations += 1;
        const stripAudit = auditAnnularStrips(arr, trace, matching, colorMask);
        const result = localAudit(arr, choice, trace, matching, colorMask, {
          ...options,
          stripAudit,
        }, null);
        recordResult(result, stripAudit, false);
      }
    }
  }
  return totals;
}

function enumerate(options = {}) {
  const degree = options.degree || 4;
  const targetGenus = options.genus === undefined ? 0 : options.genus;
  const selfSwitchMode = options.selfSwitchMode || "strong";
  if (!new Set(["strong", "leaf", "quotient"]).has(selfSwitchMode)) {
    throw new Error(`unknown self-switch mode ${selfSwitchMode}`);
  }
  const annularPhaseRange = options.annularPhaseRange;
  if (annularPhaseRange !== undefined
      && (!Number.isInteger(annularPhaseRange) || annularPhaseRange < 0)) {
    throw new Error("annular phase range must be a nonnegative integer");
  }
  if (!Number.isInteger(targetGenus) || targetGenus < 0) {
    throw new Error(`genus must be a nonnegative integer, got ${targetGenus}`);
  }
  const arr = buildArrangement(degree);
  arr.targetGenus = targetGenus;
  arr.targetChi = 2 - 2 * targetGenus - arr.boundaryComponents;
  // Every disk eye has two V tips and every V switch supplies two tips,
  // whereas annular eyes have neither.  Thus D=V and chi=D-(H+V)=-H.
  arr.requiredH = -arr.targetChi;
  const horizontalCrossings = arr.vertices.filter((vertex) => vertex.orientation === "H").length;
  if (arr.requiredH > horizontalCrossings) {
    throw new Error(
      `genus ${targetGenus} requires ${arr.requiredH} horizontal switches, but only ${horizontalCrossings} exist`,
    );
  }
  const order = buildAnnularOrder(arr, options.order || "xy");
  const prefix = options.prefix || "";
  if (typeof prefix !== "string" || !/^[01]*$/.test(prefix) || prefix.length > order.length) {
    throw new Error(`prefix must be a 0/1 string of length at most ${order.length}`);
  }
  const useTranslationSymmetry = options.translationSymmetry === true;
  const useCentroidSymmetry = useTranslationSymmetry && gcd(arr.requiredH, arr.d) === 1;
  const inverses = useTranslationSymmetry && !useCentroidSymmetry
    ? translationInverses(arr)
    : null;
  const residueFeasible = useCentroidSymmetry
    ? horizontalResidueFeasibility(arr, order)
    : null;
  const remainingH = new Int16Array(order.length + 1);
  const remainingV = new Int16Array(order.length + 1);
  for (let i = order.length - 1; i >= 0; i -= 1) {
    const isH = arr.vertices[order[i]].orientation === "H";
    remainingH[i] = remainingH[i + 1] + (isH ? 1 : 0);
    remainingV[i] = remainingV[i + 1] + (isH ? 0 : 1);
  }
  const dsu = new AnnularRollbackDSU(arr.portCount);
  addExternalEdges(arr, dsu);
  const choice = new Uint8Array(arr.vertices.length);
  const assigned = new Uint8Array(arr.vertices.length);
  const horizontalConstraints = [];
  const diskAuditVertexStamp = new Int32Array(arr.vertices.length);
  let diskAuditEpoch = 0;
  const maxNodes = options.maxNodes === undefined ? Infinity : options.maxNodes;
  if (!(maxNodes === Infinity || (Number.isInteger(maxNodes) && maxNodes > 0))) {
    throw new Error(`maxNodes must be a positive integer, got ${maxNodes}`);
  }
  const stats = {
    degree,
    targetGenus,
    targetEulerCharacteristic: arr.targetChi,
    boundaryComponents: arr.boundaryComponents,
    requireNormalB: options.requireNormalB !== false,
    requireBalance: options.requireBalance !== false,
    requireAdjacentStrips: options.requireAdjacentStrips === true,
    selfSwitchMode,
    annularPhaseRange: annularPhaseRange === undefined ? null : annularPhaseRange,
    phaseBalancePruning: options.phaseBalancePruning !== false,
    essentialPairingPruning: options.essentialPairingPruning !== false,
    translationSymmetry: useTranslationSymmetry,
    translationSymmetryMode: useCentroidSymmetry ? "horizontal-centroid" : (
      useTranslationSymmetry ? "lex-leaf" : "none"
    ),
    crossings: arr.vertices.length,
    requiredHorizontalSwitches: arr.requiredH,
    prefix,
    maxNodes: Number.isFinite(maxNodes) ? maxNodes : null,
    aborted: false,
    searchNodes: 0,
    completedMasks: 0,
    symmetryCanonicalMasks: 0,
    symmetrySkippedMasks: 0,
    pairableMasks: 0,
    masksWithAnnuli: 0,
    matchingConfigurations: 0,
    colorConfigurations: 0,
    boundaryDisjointConfigurations: 0,
    stripValidConfigurations: 0,
    embeddedAnnulusConfigurations: 0,
    locallyValidConfigurations: 0,
    connectedConfigurations: 0,
    masksLocallyValid: 0,
    masksConnected: 0,
    earlyRejected: Object.create(null),
    leafRejected: Object.create(null),
    configurationRejected: Object.create(null),
    distribution: Object.create(null),
    connectedDistribution: Object.create(null),
    witnesses: [],
    answers: [],
    elapsedMs: 0,
  };
  const started = Date.now();
  function count(table, key) { table[key] = (table[key] || 0) + 1; }
  const essentialBalance = new Map();
  let essentialImbalance = 0;
  function addEssentialSignatures(closedRoots) {
    const signatures = closedEssentialSignatures(arr, dsu, choice, closedRoots);
    for (const signature of signatures) {
      const before = essentialBalance.get(signature.key) || 0;
      const after = before + signature.canonicalSide;
      essentialImbalance += Math.abs(after) - Math.abs(before);
      if (after === 0) essentialBalance.delete(signature.key);
      else essentialBalance.set(signature.key, after);
    }
    return signatures;
  }
  function removeEssentialSignatures(signatures) {
    for (let i = signatures.length - 1; i >= 0; i -= 1) {
      const signature = signatures[i];
      const before = essentialBalance.get(signature.key) || 0;
      const after = before - signature.canonicalSide;
      essentialImbalance += Math.abs(after) - Math.abs(before);
      if (after === 0) essentialBalance.delete(signature.key);
      else essentialBalance.set(signature.key, after);
    }
  }

  function search(index, hSwitches, vSwitches, hSumI, hSumJ) {
    if (stats.aborted) return;
    stats.searchNodes += 1;
    if (options.verifyIncremental === true) {
      const open = dsu.openTipCounts[0] + dsu.openTipCounts[1] + dsu.openTipCounts[2];
      if (open !== dsu.components - dsu.closedCycles
          || [...dsu.openTipCounts].some((amount) => amount < 0)) {
        throw new Error("incremental open-component accounting drifted");
      }
    }
    if (stats.searchNodes >= maxNodes) {
      stats.aborted = true;
      return;
    }
    if (hSwitches > arr.requiredH || hSwitches + remainingH[index] < arr.requiredH) {
      count(stats.earlyRejected, "horizontal-count-bound");
      return;
    }
    if (useCentroidSymmetry) {
      const needed = arr.requiredH - hSwitches;
      const d = arr.d;
      const wantedI = ((-hSumI % d) + d) % d;
      const wantedJ = ((-hSumJ % d) + d) % d;
      const stride = d * d;
      if (needed < 0
          || !residueFeasible[index][needed * stride + wantedI * d + wantedJ]) {
        count(stats.earlyRejected, "translation-centroid-bound");
        return;
      }
    }
    if (dsu.diskCycles > vSwitches || dsu.diskCycles > vSwitches + remainingV[index]) {
      count(stats.earlyRejected, "disk-V-bound");
      return;
    }
    const remaining = order.length - index;
    if (dsu.components - 2 * remaining > dsu.closedCycles + 2 * remaining) {
      count(stats.earlyRejected, "cycle-cover-bound");
      return;
    }
    if (index === order.length) {
      stats.completedMasks += 1;
      if (dsu.components !== dsu.closedCycles) {
        count(stats.leafRejected, "incomplete-cycle-cover");
        return;
      }
      if (dsu.diskCycles !== vSwitches) {
        count(stats.leafRejected, "disk-V-mismatch");
        return;
      }
      if (dsu.essentialCycles === 0) {
        count(stats.leafRejected, "all-disk");
        return;
      }
      if (options.verifyIncremental === true) {
        const traced = traceCircuits(arr, choice);
        if (traced === null) throw new Error("incremental leaf failed a complete trace");
        const tracedBalance = new Map();
        for (const circuit of traced.circuits) {
          if (circuit.kind !== "essential") continue;
          const key = `${circuit.hx},${circuit.hy}`;
          const value = (tracedBalance.get(key) || 0) + circuit.canonicalSide;
          if (value === 0) tracedBalance.delete(key);
          else tracedBalance.set(key, value);
        }
        const left = [...essentialBalance].sort();
        const right = [...tracedBalance].sort();
        if (JSON.stringify(left) !== JSON.stringify(right)) {
          throw new Error(`incremental essential signature drifted: ${JSON.stringify(left)} != ${JSON.stringify(right)}`);
        }
      }
      const orbitWeight = useCentroidSymmetry
        ? arr.d * arr.d
        : (useTranslationSymmetry ? canonicalTranslationWeight(choice, inverses) : 1);
      if (orbitWeight === 0) {
        stats.symmetrySkippedMasks += 1;
        return;
      }
      stats.symmetryCanonicalMasks += 1;
      const result = options.auditMaskOverride
        ? options.auditMaskOverride(arr, choice, options)
        : auditMask(arr, choice, options);
      if (!result.ok) {
        count(stats.leafRejected, result.stage);
        return;
      }
      stats.pairableMasks += orbitWeight;
      stats.masksWithAnnuli += orbitWeight;
      stats.matchingConfigurations += orbitWeight * result.matchingConfigurations;
      stats.colorConfigurations += orbitWeight * result.colorConfigurations;
      stats.boundaryDisjointConfigurations += orbitWeight * result.boundaryDisjoint;
      stats.stripValidConfigurations += orbitWeight * result.stripValid;
      stats.embeddedAnnulusConfigurations += orbitWeight * result.embeddedValid;
      stats.locallyValidConfigurations += orbitWeight * result.locallyValid;
      stats.connectedConfigurations += orbitWeight * result.connected;
      for (const [stage, amount] of Object.entries(result.rejected)) {
        stats.configurationRejected[stage] = (stats.configurationRejected[stage] || 0)
          + orbitWeight * amount;
      }
      for (const [family, amount] of Object.entries(result.connectedDistribution)) {
        stats.connectedDistribution[family] = (stats.connectedDistribution[family] || 0)
          + orbitWeight * amount;
      }
      if (result.locallyValid > 0) stats.masksLocallyValid += orbitWeight;
      if (result.connected > 0) stats.masksConnected += orbitWeight;
      const key = `D${result.diskEyes}-A${result.annularEyes}-E${result.essentialCircuits}`;
      stats.distribution[key] = (stats.distribution[key] || 0) + orbitWeight;
      if (result.connected > 0) {
        const common = {
          maskHex: maskHex(choice),
          diskEyes: result.diskEyes,
          annularEyes: result.annularEyes,
          translationOrbitSize: orbitWeight,
        };
        if (options.collectAnswers === true) {
          for (const answer of result.witnesses) stats.answers.push({ ...common, result: answer });
        }
        if (stats.witnesses.length < 32) {
          stats.witnesses.push({ ...common, first: result.witnesses[0] });
        }
      }
      return;
    }
    const v = order[index];
    const vertex = arr.vertices[v];
    const isH = vertex.orientation === "H";
    const mustSwitch = isH && hSwitches + remainingH[index] === arr.requiredH;
    const cannotSwitch = isH && hSwitches === arr.requiredH;
    const forcedChoice = index < prefix.length ? Number(prefix[index]) : -1;
    if (!mustSwitch && forcedChoice !== 1) {
      const mark = dsu.top;
      choice[v] = 0;
      assigned[v] = 1;
      const zeroA = vertex.nonswitchPairs[0][0];
      const zeroB = vertex.nonswitchPairs[1][0];
      const closedRoots = [];
      let essentialSignatures = [];
      if (!dsu.activateZeroConstraint(zeroA, zeroB, 1)) {
        count(stats.earlyRejected, dsu.lastFailure);
      }
      else if (!addLocal(arr, dsu, v, false, closedRoots)) {
        count(stats.earlyRejected, dsu.lastFailure);
      }
      else {
        essentialSignatures = addEssentialSignatures(closedRoots);
        if (options.essentialPairingPruning !== false
            && essentialImbalance > dsu.openTipCounts[0]) {
          count(stats.earlyRejected, "essential-pairing-bound");
        } else if (violatesHorizontalClause(dsu, horizontalConstraints)) {
          count(stats.earlyRejected, "horizontal-local-clause");
        } else {
          const addedDisk = annotateClosedDisks(
            arr, dsu, choice, closedRoots, options.verifyIncremental === true,
          );
          const partialFailure = addedDisk
            ? partialClosedDiskAudit(
              arr, dsu, choice, assigned, selfSwitchMode === "quotient",
              closedRoots, diskAuditVertexStamp, ++diskAuditEpoch,
            )
            : null;
          if (partialFailure !== null) count(stats.earlyRejected, partialFailure);
          else search(index + 1, hSwitches, vSwitches, hSumI, hSumJ);
        }
      }
      removeEssentialSignatures(essentialSignatures);
      dsu.rollback(mark);
      dsu.deactivateZeroConstraint(zeroA, zeroB);
      assigned[v] = 0;
    }
    if (!cannotSwitch && forcedChoice !== 0) {
      const mark = dsu.top;
      choice[v] = 1;
      assigned[v] = 1;
      if (isH) horizontalConstraints.push(horizontalConstraint(arr, vertex));
      const zeroA = vertex.switchPairs[0][0];
      const zeroB = vertex.switchPairs[1][0];
      const zeroActive = selfSwitchMode !== "quotient"
        || dsu.activateZeroConstraint(zeroA, zeroB, 2);
      const closedRoots = [];
      let essentialSignatures = [];
      if (!zeroActive) {
        count(stats.earlyRejected, dsu.lastFailure);
      }
      else if (!addLocal(arr, dsu, v, true, closedRoots)) {
        count(stats.earlyRejected, dsu.lastFailure);
      }
      else if (selfSwitchMode === "strong" && !dsu.addInequality(
        vertex.switchPairs[0][0], vertex.switchPairs[1][0], 1n << BigInt(v),
      )) count(stats.earlyRejected, dsu.lastFailure);
      else {
        essentialSignatures = addEssentialSignatures(closedRoots);
        if (options.essentialPairingPruning !== false
            && essentialImbalance > dsu.openTipCounts[0]) {
          count(stats.earlyRejected, "essential-pairing-bound");
        } else if (violatesHorizontalClause(dsu, horizontalConstraints)) {
          count(stats.earlyRejected, "horizontal-local-clause");
        } else {
          const addedDisk = annotateClosedDisks(
            arr, dsu, choice, closedRoots, options.verifyIncremental === true,
          );
          const partialFailure = addedDisk
            ? partialClosedDiskAudit(
              arr, dsu, choice, assigned, selfSwitchMode === "quotient",
              closedRoots, diskAuditVertexStamp, ++diskAuditEpoch,
            )
            : null;
          if (partialFailure !== null) count(stats.earlyRejected, partialFailure);
          else search(
            index + 1,
            hSwitches + (isH ? 1 : 0),
            vSwitches + (isH ? 0 : 1),
            hSumI + (isH ? vertex.i : 0),
            hSumJ + (isH ? vertex.j : 0),
          );
        }
      }
      removeEssentialSignatures(essentialSignatures);
      dsu.rollback(mark);
      if (selfSwitchMode === "quotient") dsu.deactivateZeroConstraint(zeroA, zeroB);
      if (isH) horizontalConstraints.pop();
      choice[v] = 0;
      assigned[v] = 0;
    }
  }
  search(0, 0, 0, 0, 0);
  stats.elapsedMs = Date.now() - started;
  return stats;
}

function main(argv) {
  let degree = 4;
  let genus = 0;
  let order = "xy";
  let maxNodes;
  let translationSymmetry = false;
  let prefix = "";
  let selfSwitchMode = "strong";
  let requireAdjacentStrips = false;
  let phaseBalancePruning = true;
  let annularPhaseRange;
  for (let i = 0; i < argv.length; i += 1) {
    if (argv[i] === "--degree") degree = Number.parseInt(argv[++i], 10);
    else if (argv[i] === "--genus") genus = Number.parseInt(argv[++i], 10);
    else if (argv[i] === "--order") order = argv[++i];
    else if (argv[i] === "--max-nodes") maxNodes = Number.parseInt(argv[++i], 10);
    else if (argv[i] === "--translation-symmetry") translationSymmetry = true;
    else if (argv[i] === "--prefix") prefix = argv[++i];
    else if (argv[i] === "--self-switch-mode") selfSwitchMode = argv[++i];
    else if (argv[i] === "--annular-phase-range") {
      annularPhaseRange = Number.parseInt(argv[++i], 10);
    }
    else if (argv[i] === "--adjacent-strips") requireAdjacentStrips = true;
    else if (argv[i] === "--no-phase-balance-pruning") phaseBalancePruning = false;
    else if (argv[i] === "--help" || argv[i] === "-h") {
      process.stdout.write("Usage: node tmp/p2o4_annular_search.js [--degree D] [--genus G] [--order xy|yx|yrev|group|native|type] [--prefix 010...] [--max-nodes N] [--translation-symmetry] [--self-switch-mode strong|leaf|quotient] [--annular-phase-range K] [--adjacent-strips]\n");
      return;
    } else throw new Error(`unknown argument ${argv[i]}`);
  }
  process.stdout.write(`${JSON.stringify(enumerate({
    degree, genus, order, prefix, maxNodes, translationSymmetry, selfSwitchMode,
    annularPhaseRange, requireAdjacentStrips, phaseBalancePruning,
  }), null, 2)}\n`);
}

if (require.main === module) {
  try { main(process.argv.slice(2)); }
  catch (error) {
    process.stderr.write(`${error.stack || error}\n`);
    process.exitCode = 1;
  }
}

module.exports = {
  buildArrangement,
  enumerate,
  traceCircuits,
  auditMask,
  auditAnnularStrips,
  phaseAnnularRegion,
  perfectMatchings,
  localAudit,
  translatedDiskMultiplicity,
};
