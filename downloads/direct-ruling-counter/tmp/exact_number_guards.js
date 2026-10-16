"use strict";

// Geometry is serialized as JSON integers. Reject an unsupported magnitude
// before a Number or a typed array can silently change the exact diagram.
function safeInteger(value, label) {
  if (!Number.isSafeInteger(value)) {
    throw new RangeError(`${label} must be a safe integer`);
  }
  return value;
}

function exactTurnSign(ax, ay, bx, by) {
  const cross = BigInt(ax) * BigInt(by) - BigInt(ay) * BigInt(bx);
  return cross > 0n ? 1 : cross < 0n ? -1 : 0;
}

function validateArrangementNumbers(raw) {
  const size = safeInteger(raw.portCount, "portCount");
  if (size < 0 || size > 0x7fffffff || raw.vertices.length > 0x7fffffff) {
    throw new RangeError("arrangement indices exceed the signed 32-bit range");
  }
  if (safeInteger(raw.L, "arrangement period") <= 0) {
    throw new RangeError("arrangement period must be positive");
  }
  for (const name of ["targetChi", "targetGenus", "boundaryComponents"]) {
    safeInteger(raw[name], name);
  }
  let sx = 0n;
  let sy = 0n;
  for (const name of ["portVertex", "portTwin", "portVectorX", "portVectorY",
    "portShiftX", "portShiftY", "portSign"]) {
    if (raw[name].length !== size) throw new Error(`${name} has the wrong length`);
    for (const value of raw[name]) safeInteger(value, name);
  }
  const abs = (value) => value < 0n ? -value : value;
  for (let port = 0; port < size; port += 1) {
    if (raw.portTwin[port] < 0 || raw.portTwin[port] >= size
        || raw.portVertex[port] < 0 || raw.portVertex[port] >= raw.vertices.length) {
      throw new RangeError("half-edge index is outside the arrangement");
    }
    if (raw.portSign[port] !== -1 && raw.portSign[port] !== 1) {
      throw new RangeError("half-edge coorientation sign must be -1 or 1");
    }
    sx += abs(BigInt(raw.portShiftX[port]));
    sy += abs(BigInt(raw.portShiftY[port]));
  }
  let x = 0n;
  let y = 0n;
  for (const vertex of raw.vertices) {
    x = x > abs(BigInt(safeInteger(vertex.x, "vertex x")))
      ? x : abs(BigInt(vertex.x));
    y = y > abs(BigInt(safeInteger(vertex.y, "vertex y")))
      ? y : abs(BigInt(vertex.y));
    safeInteger(vertex.index, "vertex index");
  }
  for (const edge of raw.edges) {
    for (const key of ["a", "b", "sx", "sy"]) safeInteger(edge[key], `edge ${key}`);
    if (edge.a < 0 || edge.a >= size || edge.b < 0 || edge.b >= size
        || raw.portTwin[edge.a] !== edge.b) {
      throw new RangeError("edge endpoints disagree with the half-edge indices");
    }
    if (edge.sx !== raw.portShiftX[edge.a] || edge.sy !== raw.portShiftY[edge.a]) {
      throw new Error("edge and half-edge deck shifts disagree");
    }
  }
  const limit = BigInt(Number.MAX_SAFE_INTEGER);
  // Every traced path uses each port at most once. These bounds also cover
  // differences of lifted coordinates and the structural self-switch test.
  const liftedX = x + BigInt(raw.L) * sx;
  const liftedY = y + BigInt(raw.L) * sy;
  if (4n * (sx + sy) > limit || 2n * liftedX > limit
      || 2n * liftedY > limit || 8n * sx * sy > limit) {
    throw new RangeError("diagram exceeds the exact integer geometry range");
  }
}

function validatePayloadNumbers(payload) {
  for (const name of ["tau", "sigma"]) {
    if (payload[name] !== undefined) {
      for (const entry of payload[name]) safeInteger(entry, name);
    }
  }
  if (payload.requiredHorizontalSwitches !== undefined) {
    safeInteger(payload.requiredHorizontalSwitches, "requiredHorizontalSwitches");
  }
}

module.exports = { safeInteger, exactTurnSign, validateArrangementNumbers, validatePayloadNumbers };
