(() => {
  "use strict";

  function signedDoubleArea(vertices) {
    return vertices.reduce((sum, point, index) => {
      const next = vertices[(index + 1) % vertices.length];
      return sum + point[0] * next[1] - point[1] * next[0];
    }, 0);
  }

  function pointOnSegment(point, start, finish) {
    const cross =
      (finish[0] - start[0]) * (point[1] - start[1]) -
      (finish[1] - start[1]) * (point[0] - start[0]);
    if (cross !== 0) {
      return false;
    }
    return (
      point[0] >= Math.min(start[0], finish[0]) &&
      point[0] <= Math.max(start[0], finish[0]) &&
      point[1] >= Math.min(start[1], finish[1]) &&
      point[1] <= Math.max(start[1], finish[1])
    );
  }

  function pointLocation(vertices, point) {
    for (let index = 0; index < vertices.length; index += 1) {
      if (
        pointOnSegment(
          point,
          vertices[index],
          vertices[(index + 1) % vertices.length]
        )
      ) {
        return "boundary";
      }
    }
    let inside = false;
    for (let index = 0, previous = vertices.length - 1;
      index < vertices.length;
      previous = index, index += 1) {
      const currentPoint = vertices[index];
      const previousPoint = vertices[previous];
      const crossesRay =
        (currentPoint[1] > point[1]) !== (previousPoint[1] > point[1]);
      if (!crossesRay) {
        continue;
      }
      const intersectionX =
        ((previousPoint[0] - currentPoint[0]) *
          (point[1] - currentPoint[1])) /
          (previousPoint[1] - currentPoint[1]) +
        currentPoint[0];
      if (point[0] < intersectionX) {
        inside = !inside;
      }
    }
    return inside ? "interior" : "outside";
  }

  function formatVertices(vertices) {
    return `Conv{${vertices.map(([x, y]) => `(${x},${y})`).join(", ")}}`;
  }

  function formatArea(doubleArea) {
    return doubleArea % 2 === 0
      ? String(doubleArea / 2)
      : `${doubleArea}/2`;
  }

  window.RulingLabCore = Object.freeze({
    formatArea,
    formatVertices,
    pointLocation,
    signedDoubleArea,
  });
})();
