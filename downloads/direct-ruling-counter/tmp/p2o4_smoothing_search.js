#!/usr/bin/env node
"use strict";

/*
 * Span-free all-disk smoothing search for the degree-d plane triangle
 * Delta = conv((0,0),(d,0),(0,d)).  The three families A,B,C contain d
 * cooriented primitive geodesics.  We use the square torus of period 5d:
 *
 *   A_i : y = 5i,             future tangent ( 1, 0), sign +
 *   B_j : x = 5j+1,           future tangent ( 0, 1), sign -
 *   C_k : x+y = 5k+2,         future tangent ( 1,-1), sign -.
 *
 * The offsets 0,1,2 remove triple intersections without changing the
 * isotopy class.  AB and AC crossings have the tip (vertical) switch;
 * BC crossings have the through-going (horizontal) switch.  For an all-disk
 * ruling, D=V and chi=D-(V+H)=-H.  Since the original Legendrian has 3d
 * components, a connected ruling surface of genus g has
 *
 *     chi = 2 - 2g - 3d,       H = 3d + 2g - 2.
 *
 * The rollback DSU records full Z^2 deck voltage and the number of tip arcs
 * in every partial component.  It rejects noncontractible cycles, cycles
 * with other than two tips, tip overflow, and forced same-lift nonswitch
 * self-crossings before the expensive geometric/local inspection.
 *
 * At final inspection we check every condition used by the O(2,b) audit:
 * embedded contractible disk eyes, consistent red/blue color, all local
 * nonswitch sectors, convex/convex type A, convex/reflex type B in either
 * order (so inward B is allowed), connected orientable skeleton, the requested
 * target genus, base-point balance, and exact edge coverage.
 */

const crypto = require("node:crypto");

function cross(ax, ay, bx, by) {
  return ax * by - ay * bx;
}

function binomial(n, k) {
  if (k < 0 || k > n) return 0;
  k = Math.min(k, n - k);
  let answer = 1;
  for (let i = 1; i <= k; i += 1) answer = (answer * (n - k + i)) / i;
  return answer;
}

class WeightedRollbackDSU {
  constructor(n) {
    this.parent = new Int16Array(n);
    this.size = new Int16Array(n);
    this.potX = new Int16Array(n);
    this.potY = new Int16Array(n);
    this.tips = new Int8Array(n);
    this.closed = new Uint8Array(n);
    for (let i = 0; i < n; i += 1) {
      this.parent[i] = i;
      this.size[i] = 1;
    }
    this.components = n;
    this.closedCycles = 0;
    // A search node performs at most two unions per arrangement vertex.
    // Typed rollback storage avoids millions of short-lived object records.
    const capacity = 4 * n + 16;
    this.stackKind = new Int8Array(capacity);
    this.stackChild = new Int16Array(capacity);
    this.stackRoot = new Int16Array(capacity);
    this.stackSize = new Int16Array(capacity);
    this.stackTips = new Int8Array(capacity);
    this.stackClosed = new Uint8Array(capacity);
    this.top = 0;
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

  unite(a, b, dx, dy, addedTips = 0) {
    const fa = this.find(a);
    const fb = this.find(b);
    let ra = fa[0];
    let rb = fb[0];
    let relX = fa[1] + dx - fb[1];
    let relY = fa[2] + dy - fb[2];

    if (ra === rb) {
      if (relX !== 0 || relY !== 0) {
        this.lastFailure = "noncontractible-cycle";
        return false;
      }
      if (this.closed[ra]) {
        this.lastFailure = "duplicate-cycle";
        return false;
      }
      const newTips = this.tips[ra] + addedTips;
      if (newTips !== 2) {
        this.lastFailure = "closed-wrong-tip-count";
        return false;
      }
      const at = this.top++;
      this.stackKind[at] = 0;
      this.stackRoot[at] = ra;
      this.stackTips[at] = this.tips[ra];
      this.stackClosed[at] = this.closed[ra];
      this.tips[ra] = newTips;
      this.closed[ra] = 1;
      this.closedCycles += 1;
      return true;
    }

    if (this.closed[ra] || this.closed[rb]) {
      this.lastFailure = "join-closed-cycle";
      return false;
    }
    const newTips = this.tips[ra] + this.tips[rb] + addedTips;
    if (newTips > 2) {
      this.lastFailure = "tip-overflow";
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
    this.stackClosed[at] = this.closed[ra];
    this.parent[rb] = ra;
    this.potX[rb] = relX;
    this.potY[rb] = relY;
    this.size[ra] += this.size[rb];
    this.tips[ra] = newTips;
    this.components -= 1;
    return true;
  }

  rollback(mark) {
    while (this.top > mark) {
      const at = --this.top;
      const root = this.stackRoot[at];
      if (this.stackKind[at] === 0) {
        this.tips[root] = this.stackTips[at];
        this.closed[root] = this.stackClosed[at];
        this.closedCycles -= 1;
      } else {
        const child = this.stackChild[at];
        this.parent[child] = child;
        this.potX[child] = 0;
        this.potY[child] = 0;
        this.size[root] = this.stackSize[at];
        this.tips[root] = this.stackTips[at];
        this.closed[root] = this.stackClosed[at];
        this.components += 1;
      }
    }
  }
}

class SmallParityDSU {
  constructor(n) {
    this.parent = new Int16Array(n);
    this.parity = new Int8Array(n);
    for (let i = 0; i < n; i += 1) this.parent[i] = i;
  }

  find(a) {
    if (this.parent[a] === a) return [a, 0];
    const [root, p] = this.find(this.parent[a]);
    this.parity[a] ^= p;
    this.parent[a] = root;
    return [root, this.parity[a]];
  }

  unite(a, b, wantedParity) {
    const fa = this.find(a);
    const fb = this.find(b);
    if (fa[0] === fb[0]) return (fa[1] ^ fb[1]) === wantedParity;
    this.parent[fb[0]] = fa[0];
    this.parity[fb[0]] = fa[1] ^ fb[1] ^ wantedParity;
    return true;
  }
}

function mod(value, period) {
  value %= period;
  return value < 0 ? value + period : value;
}

function buildArrangement(d, targetGenus = 0) {
  if (!Number.isInteger(d) || d < 1) throw new Error("degree must be a positive integer");
  if (!Number.isInteger(targetGenus) || targetGenus < 0) {
    throw new Error("target genus must be a nonnegative integer");
  }
  const L = 5 * d;
  const vertices = [];
  const vertexByLabel = new Map();

  function addVertex(type, i, j, x, y) {
    const index = vertices.length;
    const vertex = {
      index,
      type,
      i,
      j,
      x: mod(x, L),
      y: mod(y, L),
      portNames: null,
      portOf: Object.create(null),
    };
    vertices.push(vertex);
    vertexByLabel.set(`${type}:${i}:${j}`, index);
  }

  for (let i = 0; i < d; i += 1) {
    for (let j = 0; j < d; j += 1) addVertex("AB", i, j, 5 * j + 1, 5 * i);
  }
  for (let i = 0; i < d; i += 1) {
    for (let k = 0; k < d; k += 1) addVertex("AC", i, k, 5 * k + 2 - 5 * i, 5 * i);
  }
  for (let j = 0; j < d; j += 1) {
    for (let k = 0; k < d; k += 1) addVertex("BC", j, k, 5 * j + 1, 5 * (k - j) + 1);
  }

  // CCW port orders.  Opposite entries are the original through-strands.
  const orders = {
    AB: ["A+", "B+", "A-", "B-"],
    AC: ["A+", "C-", "A-", "C+"],
    BC: ["B+", "C-", "B-", "C+"],
  };
  const portCount = 4 * vertices.length;
  const portVertex = new Int16Array(portCount);
  const portVectorX = new Int8Array(portCount);
  const portVectorY = new Int8Array(portCount);
  const portSign = new Int8Array(portCount);
  const portTwin = new Int16Array(portCount);
  const portShiftX = new Int8Array(portCount);
  const portShiftY = new Int8Array(portCount);
  portTwin.fill(-1);
  const vectors = {
    "A+": [1, 0, 1], "A-": [-1, 0, -1],
    "B+": [0, 1, -1], "B-": [0, -1, 1],
    "C+": [1, -1, -1], "C-": [-1, 1, 1],
  };
  for (const vertex of vertices) {
    vertex.portNames = orders[vertex.type];
    vertex.ports = [];
    for (let q = 0; q < 4; q += 1) {
      const port = 4 * vertex.index + q;
      const name = vertex.portNames[q];
      vertex.ports.push(port);
      vertex.portOf[name] = port;
      portVertex[port] = vertex.index;
      portVectorX[port] = vectors[name][0];
      portVectorY[port] = vectors[name][1];
      portSign[port] = vectors[name][2];
    }
    vertex.nonswitchPairs = [[vertex.ports[0], vertex.ports[2]], [vertex.ports[1], vertex.ports[3]]];
    const candidates = [
      [[vertex.ports[0], vertex.ports[1]], [vertex.ports[2], vertex.ports[3]]],
      [[vertex.ports[0], vertex.ports[3]], [vertex.ports[1], vertex.ports[2]]],
    ];
    const admissible = candidates.filter((pairs) => pairs.every(
      ([a, b]) => portSign[a] === -portSign[b],
    ));
    if (admissible.length !== 1) throw new Error(`nonunique switch at ${vertex.index}`);
    vertex.switchPairs = admissible[0];
    const future = new Set(vertex.portNames.map((name, q) => name.endsWith("+") ? vertex.ports[q] : -1));
    vertex.futurePorts = future;
    const pattern = vertex.switchPairs.map(([a, b]) => future.has(a) === future.has(b));
    vertex.orientation = pattern.every(Boolean) ? "V" : pattern.every((x) => !x) ? "H" : "invalid";
    if (vertex.orientation === "invalid") throw new Error(`invalid switch orientation at ${vertex.index}`);
  }

  const edges = [];
  function addEdge(origin, plusName, target, minusName, sx, sy, family) {
    const a = vertices[origin].portOf[plusName];
    const b = vertices[target].portOf[minusName];
    if (portTwin[a] !== -1 || portTwin[b] !== -1) throw new Error("port used by two external edges");
    portTwin[a] = b;
    portTwin[b] = a;
    portShiftX[a] = sx;
    portShiftY[a] = sy;
    portShiftX[b] = -sx;
    portShiftY[b] = -sy;
    edges.push({ a, b, sx, sy, family });
  }

  for (let i = 0; i < d; i += 1) {
    const line = [];
    for (let j = 0; j < d; j += 1) line.push(vertexByLabel.get(`AB:${i}:${j}`));
    for (let k = 0; k < d; k += 1) line.push(vertexByLabel.get(`AC:${i}:${k}`));
    line.sort((a, b) => vertices[a].x - vertices[b].x);
    for (let q = 0; q < line.length; q += 1) {
      const origin = line[q];
      const target = line[(q + 1) % line.length];
      addEdge(origin, "A+", target, "A-", vertices[target].x <= vertices[origin].x ? 1 : 0, 0, `A${i}`);
    }
  }
  for (let j = 0; j < d; j += 1) {
    const line = [];
    for (let i = 0; i < d; i += 1) line.push(vertexByLabel.get(`AB:${i}:${j}`));
    for (let k = 0; k < d; k += 1) line.push(vertexByLabel.get(`BC:${j}:${k}`));
    line.sort((a, b) => vertices[a].y - vertices[b].y);
    for (let q = 0; q < line.length; q += 1) {
      const origin = line[q];
      const target = line[(q + 1) % line.length];
      addEdge(origin, "B+", target, "B-", 0, vertices[target].y <= vertices[origin].y ? 1 : 0, `B${j}`);
    }
  }
  for (let k = 0; k < d; k += 1) {
    const line = [];
    for (let i = 0; i < d; i += 1) line.push(vertexByLabel.get(`AC:${i}:${k}`));
    for (let j = 0; j < d; j += 1) line.push(vertexByLabel.get(`BC:${j}:${k}`));
    line.sort((a, b) => vertices[a].x - vertices[b].x);
    for (let q = 0; q < line.length; q += 1) {
      const origin = line[q];
      const target = line[(q + 1) % line.length];
      const sx = vertices[target].x <= vertices[origin].x ? 1 : 0;
      const delta = vertices[target].x + sx * L - vertices[origin].x;
      const numerator = vertices[origin].y - delta - vertices[target].y;
      if (numerator % L !== 0) throw new Error(`nonintegral C shift at ${origin}`);
      addEdge(origin, "C+", target, "C-", sx, numerator / L, `C${k}`);
    }
  }
  if (edges.length !== 6 * d * d) throw new Error("wrong external edge count");
  if ([...portTwin].some((twin) => twin < 0)) throw new Error("unpaired port");

  return {
    d,
    L,
    vertices,
    edges,
    portCount,
    portVertex,
    portVectorX,
    portVectorY,
    portSign,
    portTwin,
    portShiftX,
    portShiftY,
    boundaryComponents: 3 * d,
    targetGenus,
    targetChi: 2 - 2 * targetGenus - 3 * d,
    requiredH: 3 * d + 2 * targetGenus - 2,
  };
}

function addExternalEdges(arr, dsu) {
  for (const edge of arr.edges) {
    if (!dsu.unite(edge.a, edge.b, edge.sx, edge.sy, 0)) {
      throw new Error(`external edge failed: ${dsu.lastFailure}`);
    }
  }
  dsu.top = 0;
}

function addLocalPairing(arr, dsu, vertexId, switched) {
  const vertex = arr.vertices[vertexId];
  const pairs = switched ? vertex.switchPairs : vertex.nonswitchPairs;
  const tip = switched && vertex.orientation === "V" ? 1 : 0;
  for (const [left, right] of pairs) {
    if (!dsu.unite(left, right, 0, 0, tip)) return false;
  }
  if (!switched) {
    const relation = dsu.relation(pairs[0][0], pairs[1][0]);
    if (relation !== null && relation[0] === 0 && relation[1] === 0) {
      dsu.lastFailure = "forced-nonswitch-self-intersection";
      return false;
    }
  }
  return true;
}

/*
 * A nonswitch contributes two straight arcs.  They are forbidden from ever
 * becoming part of the same eye at the same lift.  Checking only when the
 * nonswitch is inserted misses a later union that joins the two partial
 * paths; this was the source of almost every expensive leaf rejection in
 * the d=3 profile.  There are at most 48 constraints for d=4, so the direct
 * scan is cheap compared with completing those branches.
 */
function hasForcedNonswitchSelfIntersection(dsu, constraints) {
  for (let i = 0; i < constraints.length; i += 1) {
    const pair = constraints[i];
    const relation = dsu.relation(pair[0], pair[1]);
    if (relation !== null && relation[0] === 0 && relation[1] === 0) return true;
  }
  return false;
}

/*
 * At a switch the two smoothing branches are the boundaries of two distinct
 * eyes.  Once their partial paths enter the same DSU component they can never
 * split again, so this is a monotone rejection exactly analogous to the
 * nonswitch self-intersection test.  Unlike the nonswitch test, any deck
 * relation is forbidden here: the quotient eye orbit itself must be distinct.
 */
function hasForcedSwitchSelfAttachment(dsu, constraints) {
  for (let i = 0; i < constraints.length; i += 1) {
    const pair = constraints[i];
    if (dsu.relation(pair[0], pair[1]) !== null) return true;
  }
  return false;
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

function translatedMultiplicity(arr, eye) {
  const px = 0.37;
  const py = 0.73;
  let count = 0;
  const kMin = Math.ceil((eye.xmin - px) / arr.L);
  const kMax = Math.floor((eye.xmax - px) / arr.L);
  const lMin = Math.ceil((eye.ymin - py) / arr.L);
  const lMax = Math.floor((eye.ymax - py) / arr.L);
  for (let k = kMin; k <= kMax; k += 1) {
    for (let l = lMin; l <= lMax; l += 1) {
      if (pointInPolygon(eye.polygon, px + arr.L * k, py + arr.L * l)) count += 1;
    }
  }
  return count;
}

function occupiedSectors(arr, vertex, port, color) {
  const i = vertex.ports.indexOf(port);
  const left = new Set([i, (i + 1) & 3]);
  const right = new Set([(i + 2) & 3, (i + 3) & 3]);
  const cooriented = arr.portSign[port] === 1 ? left : right;
  const opposite = arr.portSign[port] === 1 ? right : left;
  return color === 1 ? cooriented : opposite;
}

function intersection(left, right) {
  return [...left].filter((value) => right.has(value));
}

function sliceIntersections(polygon, level) {
  const answer = [];
  for (let i = 0; i < polygon.length; i += 1) {
    const start = polygon[i];
    const end = polygon[(i + 1) % polygon.length];
    const a = 2 * start[0] + start[1];
    const b = 2 * end[0] + end[1];
    if ((a < level && level <= b) || (b < level && level <= a)) {
      const ratio = (level - a) / (b - a);
      answer.push([
        start[0] + ratio * (end[0] - start[0]),
        start[1] + ratio * (end[1] - start[1]),
      ]);
    }
  }
  return answer;
}

function typeBCompanionRelation(arr, eyes, vertex, ids) {
  const pairs = vertex.switchPairs;
  const entries = [];
  for (let index = 0; index < 2; index += 1) {
    const key = [...pairs[index]].sort((a, b) => a - b);
    const matches = eyes[ids[index]].visits.filter((visit) => {
      if (visit[0] !== vertex.index) return false;
      const local = [visit[1], visit[2]].sort((a, b) => a - b);
      return local[0] === key[0] && local[1] === key[1];
    });
    if (matches.length !== 1) return "inconsistent";
    entries.push(matches[0]);
  }
  const commonX = entries[0][3];
  const commonY = entries[0][4];
  const polygons = entries.map((entry, index) => {
    const shiftX = commonX - entry[3];
    const shiftY = commonY - entry[4];
    return eyes[ids[index]].polygon.map(([x, y]) => [x + shiftX, y + shiftY]);
  });
  const switchLevel = 2 * commonX + commonY;
  const relations = new Set();
  for (const side of [-1, 1]) {
    const hits = polygons.map((polygon) => (
      sliceIntersections(polygon, switchLevel + side * 0.1)
    ));
    if (hits.some((points) => points.length !== 2)) return "inconsistent";
    const split = hits.map((points) => [...points].sort((left, right) => {
      const leftDistance = (left[0] - commonX) ** 2 + (left[1] - commonY) ** 2;
      const rightDistance = (right[0] - commonX) ** 2 + (right[1] - commonY) ** 2;
      return leftDistance - rightDistance;
    }));
    const nearDifference = split[1][0][0] - split[0][0][0];
    const farDifference = split[1][1][0] - split[0][1][0];
    if (nearDifference === 0 || farDifference === 0) return "inconsistent";
    relations.add(nearDifference * farDifference > 0 ? "same" : "opposite");
  }
  return relations.size === 1 ? [...relations][0] : "inconsistent";
}

function inspectChoices(arr, choice, options = {}) {
  const requireBalance = options.requireBalance !== false;
  const requireNormalB = options.requireNormalB === true;
  const requireDistinctSwitchEyes = options.requireDistinctSwitchEyes === true;
  const seen = new Uint8Array(arr.portCount);
  const eyeOf = new Int16Array(arr.portCount);
  const convexAt = new Uint8Array(arr.portCount);
  eyeOf.fill(-1);
  const eyes = [];

  for (let start = 0; start < arr.portCount; start += 1) {
    if (seen[start]) continue;
    const eyeIndex = eyes.length;
    const startVertex = arr.portVertex[start];
    let liftX = 0;
    let liftY = 0;
    let current = start;
    let area2 = 0;
    let tips = 0;
    let xmin = arr.vertices[startVertex].x;
    let xmax = xmin;
    let ymin = arr.vertices[startVertex].y;
    let ymax = ymin;
    const polygon = [];
    const visits = [];
    const liftedVertices = new Set();

    for (let guard = 0; ; guard += 1) {
      if (guard > arr.portCount) return { ok: false, stage: "trace-branch" };
      if (guard > 0 && current === start) break;
      if (seen[current]) return { ok: false, stage: "trace-branch" };
      const v = arr.portVertex[current];
      const vertex = arr.vertices[v];
      const switched = choice[v] !== 0;
      const key = `${v}:${liftX}:${liftY}`;
      if (liftedVertices.has(key) && !switched) return { ok: false, stage: "self-intersection" };
      liftedVertices.add(key);
      const pairs = switched ? vertex.switchPairs : vertex.nonswitchPairs;
      let outgoing = -1;
      for (const [a, b] of pairs) {
        if (a === current) outgoing = b;
        else if (b === current) outgoing = a;
      }
      if (outgoing < 0) return { ok: false, stage: "missing-local-mate" };
      seen[current] = 1;
      seen[outgoing] = 1;
      eyeOf[current] = eyeIndex;
      eyeOf[outgoing] = eyeIndex;
      const ax = vertex.x + liftX * arr.L;
      const ay = vertex.y + liftY * arr.L;
      polygon.push([ax, ay]);
      if (switched && vertex.orientation === "V") tips += 1;
      visits.push([v, current, outgoing, ax, ay]);

      liftX += arr.portShiftX[outgoing];
      liftY += arr.portShiftY[outgoing];
      const next = arr.portTwin[outgoing];
      const target = arr.vertices[arr.portVertex[next]];
      const bx = target.x + liftX * arr.L;
      const by = target.y + liftY * arr.L;
      area2 += ax * by - ay * bx;
      xmin = Math.min(xmin, bx);
      xmax = Math.max(xmax, bx);
      ymin = Math.min(ymin, by);
      ymax = Math.max(ymax, by);
      current = next;
    }
    if (liftX !== 0 || liftY !== 0) return { ok: false, stage: "noncontractible" };
    if (area2 === 0) return { ok: false, stage: "zero-area" };
    if (tips !== 2) return { ok: false, stage: "not-two-V-endpoints" };
    const areaSign = area2 > 0 ? 1 : -1;
    let color = 0;
    for (const [v, incoming, outgoing] of visits) {
      const edgeColor = arr.portSign[outgoing] === areaSign ? 1 : -1;
      if (color !== 0 && color !== edgeColor) return { ok: false, stage: "color-inconsistent" };
      color = edgeColor;
      if (choice[v]) {
        const inX = -arr.portVectorX[incoming];
        const inY = -arr.portVectorY[incoming];
        const turn = cross(inX, inY, arr.portVectorX[outgoing], arr.portVectorY[outgoing]);
        const convex = (turn > 0) === (area2 > 0);
        convexAt[incoming] = convex ? 1 : 0;
        convexAt[outgoing] = convex ? 1 : 0;
      }
    }
    eyes.push({
      color,
      area2,
      tips,
      polygon,
      xmin,
      xmax,
      ymin,
      ymax,
      degree: 0,
      // Preserve the cyclic order of switch darts around this disk.  Euler
      // characteristic alone is not enough: different half-twisted band
      // gluings can change the actual boundary count and hence the genus.
      // Circuit seeds give arbitrary orientation.  Normalize every disk to
      // the same (positive/CCW) orientation before using its rotation system.
      switchDarts: (() => {
        const darts = visits.filter(([v]) => choice[v]).map(([v]) => v);
        if (area2 < 0) darts.reverse();
        return darts;
      })(),
      visits,
    });
  }

  // Every port belongs to one eye, and the two ends of every arrangement
  // edge must inherit the same eye/color.  This is the exact edge-cover audit.
  for (let p = 0; p < arr.portCount; p += 1) {
    if (!seen[p] || eyeOf[p] < 0) return { ok: false, stage: "uncovered-edge-port" };
    if (eyeOf[p] !== eyeOf[arr.portTwin[p]]) return { ok: false, stage: "edge-color-conflict" };
  }

  let switches = 0;
  for (const bit of choice) switches += bit;
  if (eyes.length - switches !== arr.targetChi) return { ok: false, stage: "euler" };

  const connectivity = new SmallParityDSU(eyes.length);
  const orientation = new SmallParityDSU(eyes.length);
  const localTypes = { A: 0, B: 0, V: 0, N: 0 };
  const typeBCompanionOrders = { same: 0, opposite: 0, inconsistent: 0 };
  const switchKind = new Int8Array(arr.vertices.length); // 1=A/V, 2=B
  let facingNonswitches = 0;
  let selfSwitches = 0;
  for (const vertex of arr.vertices) {
    if (!choice[vertex.index]) {
      localTypes.N += 1;
      const ids = vertex.nonswitchPairs.map((pair) => eyeOf[pair[0]]);
      if (ids.some((id) => id < 0)) return { ok: false, stage: "missing-eye" };
      const occupied = intersection(
        occupiedSectors(arr, vertex, vertex.nonswitchPairs[0][0], eyes[ids[0]].color),
        occupiedSectors(arr, vertex, vertex.nonswitchPairs[1][0], eyes[ids[1]].color),
      );
      if (occupied.length !== 1) return { ok: false, stage: "invalid-nonswitch-sector" };
      const past = new Set(vertex.ports.filter((p) => !vertex.futurePorts.has(p)));
      const pastSectors = [];
      for (let q = 0; q < 4; q += 1) {
        if (past.has(vertex.ports[q]) && past.has(vertex.ports[(q + 1) & 3])) pastSectors.push(q);
      }
      if (pastSectors.length !== 1) return { ok: false, stage: "invalid-sweep-sector" };
      if (eyes[ids[0]].color === eyes[ids[1]].color && occupied[0] === pastSectors[0]) {
        facingNonswitches += 1;
      }
      continue;
    }

    const ids = vertex.switchPairs.map((pair) => eyeOf[pair[0]]);
    const convex = vertex.switchPairs.map((pair) => convexAt[pair[0]] !== 0);
    if (ids.some((id) => id < 0)) return { ok: false, stage: "missing-eye" };
    if (ids[0] === ids[1]) selfSwitches += 1;
    let kind;
    if (vertex.orientation === "V") {
      if (eyes[ids[0]].color === eyes[ids[1]].color || !convex[0] || !convex[1]) {
        return { ok: false, stage: "invalid-V" };
      }
      kind = "V";
      localTypes.V += 1;
    } else if (eyes[ids[0]].color !== eyes[ids[1]].color) {
      if (!convex[0] || !convex[1]) return { ok: false, stage: "inward-A" };
      kind = "A";
      localTypes.A += 1;
    } else {
      // Both convex/reflex cyclic orders are allowed: in particular the
      // inward type-B picture is not discarded.
      if (convex[0] === convex[1]) return { ok: false, stage: "invalid-B-sectors" };
      kind = "B";
      localTypes.B += 1;
      const companionRelation = typeBCompanionRelation(arr, eyes, vertex, ids);
      typeBCompanionOrders[companionRelation] += 1;
      if (requireNormalB && companionRelation !== "same") {
        return {
          ok: false,
          stage: "non-normal-B",
          details: { crossing: vertex.index, companionRelation },
        };
      }
    }
    eyes[ids[0]].degree += 1;
    eyes[ids[1]].degree += 1;
    switchKind[vertex.index] = kind === "B" ? 2 : 1;
    connectivity.unite(ids[0], ids[1], 0);
    if (!orientation.unite(ids[0], ids[1], kind === "B" ? 0 : 1)) {
      return { ok: false, stage: "nonorientable" };
    }
  }

  const roots = new Set();
  for (let i = 0; i < eyes.length; i += 1) roots.add(connectivity.find(i)[0]);
  if (roots.size !== 1) return { ok: false, stage: "skeleton-disconnected" };

  /*
   * Exact ribbon boundary trace.  A dart is one switch corner in the cyclic
   * order around an eye, with side flags 0/1.  The disk boundary joins
   * (dart,1) to (next,0).  A/V half-twisted bands join equal labels; a type-B
   * convex/reflex attachment reverses the side label and joins opposite
   * labels.  The resulting 2-regular flag graph is the actual boundary of
   * the ruling surface.  Circuit seeds orient eyes arbitrarily, so the eye
   * rotations above were first normalized by signed area; omitting that
   * normalization gives a false genus obstruction.
   */
  const dartsByEye = eyes.map(() => []);
  const dartsAtVertex = Array.from({ length: arr.vertices.length }, () => []);
  let dartCount = 0;
  for (let eye = 0; eye < eyes.length; eye += 1) {
    for (const vertex of eyes[eye].switchDarts) {
      dartsByEye[eye].push(dartCount);
      dartsAtVertex[vertex].push(dartCount);
      dartCount += 1;
    }
  }
  const flagAdj = Array.from({ length: 2 * dartCount }, () => []);
  function joinFlag(left, right) {
    flagAdj[left].push(right);
    flagAdj[right].push(left);
  }
  for (const darts of dartsByEye) {
    if (darts.length === 0) return { ok: false, stage: "eye-without-switch" };
    for (let i = 0; i < darts.length; i += 1) {
      joinFlag(2 * darts[i] + 1, 2 * darts[(i + 1) % darts.length]);
    }
  }
  for (let v = 0; v < arr.vertices.length; v += 1) {
    if (!choice[v]) continue;
    const darts = dartsAtVertex[v];
    if (darts.length !== 2) return { ok: false, stage: "bad-switch-dart-count" };
    const [left, right] = darts;
    if (switchKind[v] === 2) {
      joinFlag(2 * left, 2 * right + 1);
      joinFlag(2 * left + 1, 2 * right);
    } else {
      joinFlag(2 * left, 2 * right);
      joinFlag(2 * left + 1, 2 * right + 1);
    }
  }
  if (flagAdj.some((neighbors) => neighbors.length !== 2)) {
    return { ok: false, stage: "bad-ribbon-flag-degree" };
  }
  const flagSeen = new Uint8Array(flagAdj.length);
  const ribbonCycleLengths = [];
  for (let seed = 0; seed < flagAdj.length; seed += 1) {
    if (flagSeen[seed]) continue;
    let previous = -1;
    let current = seed;
    let length = 0;
    while (!flagSeen[current]) {
      flagSeen[current] = 1;
      length += 1;
      const next = flagAdj[current][0] === previous
        ? flagAdj[current][1]
        : flagAdj[current][0];
      previous = current;
      current = next;
    }
    ribbonCycleLengths.push(length);
  }
  ribbonCycleLengths.sort((a, b) => a - b);
  const ribbonBoundaryComponents = ribbonCycleLengths.length;
  const genusNumerator = 2 - ribbonBoundaryComponents - (eyes.length - switches);
  const ribbonGenus = genusNumerator >= 0 && genusNumerator % 2 === 0
    ? genusNumerator / 2
    : null;
  if (ribbonBoundaryComponents !== arr.boundaryComponents || ribbonGenus !== arr.targetGenus) {
    return {
      ok: false,
      stage: "ribbon-topology",
      details: { ribbonBoundaryComponents, ribbonGenus, ribbonCycleLengths },
    };
  }

  if (requireDistinctSwitchEyes && selfSwitches !== 0) {
    return { ok: false, stage: "switch-self-attachment", details: { selfSwitches } };
  }

  let baseBalance = 0;
  for (const eye of eyes) baseBalance += eye.color * translatedMultiplicity(arr, eye);
  if (requireBalance && baseBalance !== 0) return { ok: false, stage: "base-unbalanced" };

  return {
    ok: true,
    stage: "accepted",
    eyes: eyes.length,
    switches,
    chi: eyes.length - switches,
    baseBalance,
    facingNonswitches,
    selfSwitches,
    localTypes,
    typeBCompanionOrders,
    eyeDegrees: eyes.map((eye) => eye.degree).sort((a, b) => a - b),
    eyeAreas2: eyes.map((eye) => Math.abs(eye.area2)).sort((a, b) => a - b),
    edgeOrbits: arr.edges.length,
    coveredEdgeOrbits: arr.edges.length,
    ribbonBoundaryComponents,
    ribbonGenus,
    ribbonCycleLengths,
  };
}

function buildOrder(arr, name) {
  const order = arr.vertices.map((vertex) => vertex.index);
  if (name === "group") {
    order.sort((a, b) => {
      const ah = arr.vertices[a].orientation === "H" ? 0 : 1;
      const bh = arr.vertices[b].orientation === "H" ? 0 : 1;
      return ah - bh || a - b;
    });
  } else if (name === "type") {
    order.sort((a, b) => arr.vertices[a].type.localeCompare(arr.vertices[b].type) || a - b);
  } else if (name === "xy") {
    order.sort((a, b) => arr.vertices[a].x - arr.vertices[b].x || arr.vertices[a].y - arr.vertices[b].y || a - b);
  } else if (name === "yx") {
    order.sort((a, b) => arr.vertices[a].y - arr.vertices[b].y || arr.vertices[a].x - arr.vertices[b].x || a - b);
  } else if (name !== "native") {
    throw new Error(`unknown order ${name}`);
  }
  return order;
}

function maskHex(choice) {
  let mask = 0n;
  for (let i = 0; i < choice.length; i += 1) if (choice[i]) mask |= 1n << BigInt(i);
  return `0x${mask.toString(16).padStart(Math.ceil(choice.length / 4), "0")}`;
}

function familyKey(result, areas = false) {
  const t = result.localTypes;
  let key = `D${result.eyes}-A${t.A}-B${t.B}-V${t.V}-N${t.N}`
    + `-F${result.facingNonswitches}-deg[${result.eyeDegrees.join(".")}]`;
  if (areas) key += `-area2[${result.eyeAreas2.join(".")}]`;
  return key;
}

function incrementFamily(map, key, example) {
  if (!map.has(key)) map.set(key, { family: key, count: 0, examples: [] });
  const value = map.get(key);
  value.count += 1;
  if (value.examples.length < 4) value.examples.push(example);
}

function enumerate(options = {}) {
  const d = options.degree || 4;
  const targetGenus = options.targetGenus ?? 0;
  const maxSearchNodes = options.maxSearchNodes ?? Number.POSITIVE_INFINITY;
  if (!(maxSearchNodes > 0)) throw new Error("max search nodes must be positive");
  const arr = buildArrangement(d, targetGenus);
  const orderName = options.order || "xy";
  const order = buildOrder(arr, orderName);
  const collectAnswers = options.collectAnswers === true;
  const requireBalance = options.requireBalance !== false;
  const requireNormalB = options.requireNormalB === true;
  const requireDistinctSwitchEyes = options.requireDistinctSwitchEyes === true;
  const remainingH = new Int16Array(order.length + 1);
  const remainingV = new Int16Array(order.length + 1);
  for (let i = order.length - 1; i >= 0; i -= 1) {
    const isH = arr.vertices[order[i]].orientation === "H";
    remainingH[i] = remainingH[i + 1] + (isH ? 1 : 0);
    remainingV[i] = remainingV[i + 1] + (isH ? 0 : 1);
  }
  const dsu = new WeightedRollbackDSU(arr.portCount);
  addExternalEdges(arr, dsu);
  const choice = new Uint8Array(arr.vertices.length);
  const nonswitchConstraints = [];
  const switchConstraints = [];
  const families = new Map();
  const areaFamilies = new Map();
  const answers = [];
  const stats = {
    polygon: `P2 O(${d})`,
    arrangement: `${d}+${d}+${d} geodesics`,
    crossings: arr.vertices.length,
    edgeOrbits: arr.edges.length,
    boundaryComponents: arr.boundaryComponents,
    targetGenus: arr.targetGenus,
    targetEulerCharacteristic: arr.targetChi,
    horizontalCrossings: d * d,
    verticalCrossings: 2 * d * d,
    requiredHorizontalSwitches: arr.requiredH,
    theoremCandidateMasks: binomial(d * d, arr.requiredH) * (2 ** (2 * d * d)),
    order: orderName,
    searchNodes: 0,
    maxSearchNodes: Number.isFinite(maxSearchNodes) ? maxSearchNodes : null,
    truncated: false,
    requireDistinctSwitchEyes,
    completedMasks: 0,
    inspectCalls: 0,
    earlyRejected: Object.create(null),
    rejected: Object.create(null),
    accepted: 0,
    elapsedMs: 0,
    families: [],
    areaFamilies: [],
  };
  const started = Date.now();
  function count(table, key) {
    table[key] = (table[key] || 0) + 1;
  }

  function search(index, hSwitches, vSwitches) {
    if (stats.truncated) return;
    if (stats.searchNodes >= maxSearchNodes) {
      stats.truncated = true;
      return;
    }
    stats.searchNodes += 1;
    if (hSwitches > arr.requiredH || hSwitches + remainingH[index] < arr.requiredH) {
      count(stats.earlyRejected, "horizontal-count-bound");
      return;
    }
    const remaining = order.length - index;
    const minCircuits = Math.max(dsu.closedCycles, dsu.components - 2 * remaining);
    const maxCircuits = dsu.components;
    const minV = Math.max(vSwitches, dsu.closedCycles);
    const maxV = vSwitches + remainingV[index];
    if (Math.max(minCircuits, minV) > Math.min(maxCircuits, maxV)) {
      count(stats.earlyRejected, "circuit-V-bound");
      return;
    }
    if (index === order.length) {
      stats.completedMasks += 1;
      if (dsu.components !== dsu.closedCycles) {
        count(stats.earlyRejected, "incomplete-cycle-cover");
        return;
      }
      if (dsu.components !== vSwitches) {
        count(stats.earlyRejected, "eye-V-mismatch");
        return;
      }
      stats.inspectCalls += 1;
      const result = inspectChoices(arr, choice, {
        requireBalance,
        requireNormalB,
        requireDistinctSwitchEyes,
      });
      if (!result.ok) {
        count(stats.rejected, result.stage);
        return;
      }
      const hex = maskHex(choice);
      result.maskHex = hex;
      stats.accepted += 1;
      incrementFamily(families, familyKey(result, false), hex);
      incrementFamily(areaFamilies, familyKey(result, true), hex);
      if (collectAnswers) answers.push(result);
      return;
    }

    const v = order[index];
    const isH = arr.vertices[v].orientation === "H";
    const mustSwitchH = isH && hSwitches + remainingH[index] === arr.requiredH;
    const cannotSwitchH = isH && hSwitches === arr.requiredH;
    if (!mustSwitchH) {
      const mark = dsu.top;
      choice[v] = 0;
      const vertex = arr.vertices[v];
      nonswitchConstraints.push([
        vertex.nonswitchPairs[0][0],
        vertex.nonswitchPairs[1][0],
      ]);
      if (!addLocalPairing(arr, dsu, v, false)) {
        count(stats.earlyRejected, dsu.lastFailure);
      } else if (hasForcedNonswitchSelfIntersection(dsu, nonswitchConstraints)) {
        count(stats.earlyRejected, "forced-nonswitch-self-intersection-late");
      } else if (requireDistinctSwitchEyes
        && hasForcedSwitchSelfAttachment(dsu, switchConstraints)) {
        count(stats.earlyRejected, "forced-switch-self-attachment");
      } else {
        search(index + 1, hSwitches, vSwitches);
      }
      nonswitchConstraints.pop();
      dsu.rollback(mark);
    }
    if (!cannotSwitchH) {
      const mark = dsu.top;
      choice[v] = 1;
      const vertex = arr.vertices[v];
      if (requireDistinctSwitchEyes) {
        switchConstraints.push([vertex.switchPairs[0][0], vertex.switchPairs[1][0]]);
      }
      if (addLocalPairing(arr, dsu, v, true)) {
        if (hasForcedNonswitchSelfIntersection(dsu, nonswitchConstraints)) {
          count(stats.earlyRejected, "forced-nonswitch-self-intersection-late");
        } else if (requireDistinctSwitchEyes
          && hasForcedSwitchSelfAttachment(dsu, switchConstraints)) {
          count(stats.earlyRejected, "forced-switch-self-attachment");
        } else {
          search(index + 1, hSwitches + (isH ? 1 : 0), vSwitches + (isH ? 0 : 1));
        }
      } else count(stats.earlyRejected, dsu.lastFailure);
      if (requireDistinctSwitchEyes) switchConstraints.pop();
      dsu.rollback(mark);
      choice[v] = 0;
    }
  }

  search(0, 0, 0);
  stats.elapsedMs = Date.now() - started;
  stats.families = [...families.values()].sort((a, b) => b.count - a.count || a.family.localeCompare(b.family));
  stats.areaFamilies = [...areaFamilies.values()].sort((a, b) => b.count - a.count || a.family.localeCompare(b.family));
  answers.sort((a, b) => a.maskHex.localeCompare(b.maskHex));
  const canonical = answers.map((answer) => answer.maskHex).join("\n");
  stats.maskSetSha256 = collectAnswers ? crypto.createHash("sha256").update(canonical).digest("hex") : null;
  return { stats, answers };
}

function main(argv) {
  let degree = 4;
  let targetGenus = 0;
  let order = "xy";
  let json = false;
  let list = false;
  let requireBalance = true;
  let requireNormalB = false;
  let requireDistinctSwitchEyes = false;
  let maxSearchNodes = Number.POSITIVE_INFINITY;
  for (let i = 0; i < argv.length; i += 1) {
    const arg = argv[i];
    if (arg === "--degree") degree = Number.parseInt(argv[++i], 10);
    else if (arg === "--genus") targetGenus = Number.parseInt(argv[++i], 10);
    else if (arg === "--order") order = argv[++i];
    else if (arg === "--json") json = true;
    else if (arg === "--list") list = true;
    else if (arg === "--no-balance") requireBalance = false;
    else if (arg === "--normal-b") requireNormalB = true;
    else if (arg === "--distinct-switch-eyes") requireDistinctSwitchEyes = true;
    else if (arg === "--max-nodes") maxSearchNodes = Number.parseInt(argv[++i], 10);
    else if (arg === "--help" || arg === "-h") {
      process.stdout.write("Usage: node tmp/p2o4_smoothing_search.js [--degree D] [--genus G] [--order native|group|type|xy|yx] [--json] [--list] [--no-balance] [--normal-b] [--distinct-switch-eyes] [--max-nodes N]\n");
      return;
    } else throw new Error(`unknown argument ${arg}`);
  }
  const result = enumerate({
    degree,
    targetGenus,
    order,
    requireBalance,
    requireNormalB,
    requireDistinctSwitchEyes,
    maxSearchNodes,
    collectAnswers: json || list,
  });
  if (json) process.stdout.write(`${JSON.stringify(result, null, 2)}\n`);
  else {
    process.stdout.write(`${JSON.stringify(result.stats, null, 2)}\n`);
    if (list) {
      for (const answer of result.answers) {
        const t = answer.localTypes;
        process.stdout.write(`${answer.maskHex} D=${answer.eyes} A=${t.A} B=${t.B} V=${t.V} N=${t.N} deg=[${answer.eyeDegrees.join(",")}]\n`);
      }
    }
  }
}

if (require.main === module) {
  try {
    main(process.argv.slice(2));
  } catch (error) {
    process.stderr.write(`${error.stack || error}\n`);
    process.exitCode = 1;
  }
}

module.exports = { buildArrangement, enumerate, inspectChoices, buildOrder, maskHex };
