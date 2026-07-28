import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import vm from "node:vm";
import { fileURLToPath } from "node:url";


const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "..");

function runBrowserScript(filename, context) {
  const source = fs.readFileSync(path.join(root, filename), "utf8");
  vm.runInNewContext(source, context, { filename });
}

function enumerateInteriorPoints(core, vertices) {
  const xs = vertices.map(([x]) => x);
  const ys = vertices.map(([, y]) => y);
  const interior = [];
  for (let x = Math.min(...xs); x <= Math.max(...xs); x += 1) {
    for (let y = Math.min(...ys); y <= Math.max(...ys); y += 1) {
      if (core.pointLocation(vertices, [x, y]) === "interior") {
        interior.push([x, y]);
      }
    }
  }
  return interior;
}

const coreContext = { window: {} };
runBrowserScript("core.js", coreContext);
const core = coreContext.window.RulingLabCore;

assert.equal(core.signedDoubleArea([[-1, -1], [1, 0], [0, 1]]), 3);
assert.equal(core.formatArea(3), "3/2");
assert.equal(
  core.formatVertices([[-1, -1], [1, 0], [0, 1]]),
  "Conv{(-1,-1), (1,0), (0,1)}"
);

const caseContext = { window: {} };
runBrowserScript("data/cases.js", caseContext);
const cases = caseContext.window.RULING_CASES;
assert.equal(cases.length, 16, "the atlas contains all 16 classes");
assert.equal(new Set(cases.map((item) => item.id)).size, 16, "unique ids");

let drawingCount = 0;
let rulingMultiplicity = 0;
for (const item of cases) {
  assert.match(item.id, /^polygon-[0-9a-f]{10}$/);
  assert.ok(item.vertices.length >= 3);
  assert.ok(core.signedDoubleArea(item.vertices) > 0, `${item.id}: CCW`);
  assert.equal(item.interiorLatticePoints, 1, `${item.id}: stored I`);
  assert.deepEqual(
    enumerateInteriorPoints(core, item.vertices),
    [[0, 0]],
    `${item.id}: unique interior point`
  );
  assert.equal(
    item.doubleArea,
    item.boundaryLatticePoints + 2 * item.interiorLatticePoints - 2,
    `${item.id}: Pick's theorem`
  );
  assert.match(item.provenance.method, /rollback smoothing DFS/);
  assert.equal(item.provenance.usesTropicalCount, false);

  const genera = [...item.genera].sort((left, right) => left.genus - right.genus);
  assert.deepEqual(
    genera.map((entry) => entry.genus),
    [0, 1],
    `${item.id}: all possible genera`
  );
  assert.equal(genera[1].counts.total, 1, `${item.id}: standard ruling`);

  for (const entry of genera) {
    assert.equal(
      entry.counts.total,
      entry.counts.allDisk + entry.counts.annular,
      `${item.id}, genus ${entry.genus}: sector sum`
    );
    assert.ok(entry.rulings.length > 0, `${item.id}: drawings generated`);
    const multiplicity = entry.rulings.reduce(
      (sum, drawing) => sum + Number(drawing.multiplicity || 1),
      0
    );
    assert.equal(
      multiplicity,
      entry.counts.total,
      `${item.id}, genus ${entry.genus}: drawing multiplicities`
    );
    drawingCount += entry.rulings.length;
    rulingMultiplicity += multiplicity;

    for (const drawing of entry.rulings) {
      const svgPath = path.join(root, drawing.src);
      assert.ok(fs.existsSync(svgPath), `${drawing.src} exists`);
      const svg = fs.readFileSync(svgPath, "utf8");
      assert.match(svg, /<svg\b[^>]*xmlns="http:\/\/www\.w3\.org\/2000\/svg"/);
      assert.match(svg, /fill-opacity:(?:0?\.)30/);
      assert.doesNotMatch(svg, /\bDing\b|Fig(?:ure)?\.?\s*\d/i);
    }
  }
}

assert.equal(rulingMultiplicity, 112, "96 rational + 16 standard rulings");
assert.ok(drawingCount > 0 && drawingCount <= rulingMultiplicity);

const index = fs.readFileSync(path.join(root, "index.html"), "utf8");
assert.match(index, /name="robots" content="noindex, nofollow, noarchive"/);
assert.match(index, /src="core\.js"/);
assert.match(index, /src="data\/cases\.js"/);
assert.match(index, /src="app\.js"/);
assert.doesNotMatch(index, /https?:\/\/[^"]+\.js/);
assert.doesNotMatch(index, /\bDing\b|Fig(?:ure)?\.?\s*\d/i);

const packageRoot = path.join(root, "downloads", "direct-ruling-counter");
const archive = path.join(root, "downloads", "direct-ruling-counter.zip");
assert.ok(fs.existsSync(archive), "downloadable archive exists");
assert.ok(fs.statSync(archive).size > 0, "downloadable archive is nonempty");
for (const relative of [
  "README.md",
  "ruling_polygon.py",
  "example_polygon.json",
  "smoke_test.py",
]) {
  assert.ok(fs.existsSync(path.join(packageRoot, relative)), `${relative} exists`);
}
for (const entry of fs.readdirSync(packageRoot, { recursive: true })) {
  assert.doesNotMatch(String(entry), /tropical/i);
}

console.log(
  `ruling-lab validation passed: ${cases.length} polygons, ` +
  `${drawingCount} SVG representatives, multiplicity ${rulingMultiplicity}`
);
