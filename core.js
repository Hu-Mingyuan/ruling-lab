(() => {
  "use strict";

  function parseVertices(text) {
    const trimmed = text.trim();
    if (!trimmed) {
      throw new Error("Enter at least three lattice vertices.");
    }

    let value;
    try {
      value = JSON.parse(
        trimmed
          .replace(/\(/g, "[")
          .replace(/\)/g, "]")
      );
    } catch (_error) {
      const numbers = trimmed.match(/-?\d+/g);
      if (!numbers || numbers.length % 2 !== 0) {
        throw new Error(
          "Could not read the vertices. Use pairs such as 0,0; 3,0; 3,1."
        );
      }
      value = [];
      for (let index = 0; index < numbers.length; index += 2) {
        value.push([Number(numbers[index]), Number(numbers[index + 1])]);
      }
    }

    if (!Array.isArray(value)) {
      throw new Error("The polygon must be a list of coordinate pairs.");
    }
    const vertices = value.map((point) => {
      if (
        !Array.isArray(point) ||
        point.length !== 2 ||
        !Number.isInteger(point[0]) ||
        !Number.isInteger(point[1])
      ) {
        throw new Error("Every vertex must be a pair of integers.");
      }
      return [point[0], point[1]];
    });

    if (
      vertices.length > 1 &&
      samePoint(vertices[0], vertices[vertices.length - 1])
    ) {
      vertices.pop();
    }
    return vertices;
  }

  function validatePolygon(vertices) {
    if (vertices.length < 3) {
      throw new Error("A polygon needs at least three distinct vertices.");
    }
    const unique = new Set(vertices.map((point) => point.join(",")));
    if (unique.size !== vertices.length) {
      throw new Error("The vertex list contains a repeated point.");
    }
    if (signedDoubleArea(vertices) === 0) {
      throw new Error("The listed points have zero area.");
    }
  }

  function signedDoubleArea(vertices) {
    return vertices.reduce((sum, point, index) => {
      const next = vertices[(index + 1) % vertices.length];
      return sum + point[0] * next[1] - point[1] * next[0];
    }, 0);
  }

  function polygonKey(vertices) {
    const candidates = [];
    for (const sequence of [vertices, [...vertices].reverse()]) {
      for (let shift = 0; shift < sequence.length; shift += 1) {
        const rotated = sequence
          .slice(shift)
          .concat(sequence.slice(0, shift));
        const [originX, originY] = rotated[0];
        candidates.push(
          rotated
            .map(([x, y]) => `${x - originX},${y - originY}`)
            .join(";")
        );
      }
    }
    return candidates.sort()[0];
  }

  function samePoint(left, right) {
    return left[0] === right[0] && left[1] === right[1];
  }

  window.RulingLabCore = Object.freeze({
    parseVertices,
    polygonKey,
    signedDoubleArea,
    validatePolygon,
  });
})();
