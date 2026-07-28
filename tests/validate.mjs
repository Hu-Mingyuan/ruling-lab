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

const coreContext = { window: {} };
runBrowserScript("core.js", coreContext);
const core = coreContext.window.RulingLabCore;

assert.deepEqual(
  Array.from(
    core.parseVertices("[(0,0),(3,0),(3,1),(0,2)]"),
    (point) => Array.from(point)
  ),
  [[0, 0], [3, 0], [3, 1], [0, 2]]
);
assert.deepEqual(
  Array.from(
    core.parseVertices("0,0; 3,0; 3,1; 0,2; 0,0"),
    (point) => Array.from(point)
  ),
  [[0, 0], [3, 0], [3, 1], [0, 2]]
);
assert.equal(
  core.polygonKey([[0, 0], [3, 0], [3, 1], [0, 2]]),
  core.polygonKey([[8, -4], [8, -2], [11, -3], [11, -4]])
);
assert.equal(
  core.polygonKey([[0, 0], [4, 0], [1, 2]]),
  core.polygonKey([[5, 5], [6, 7], [9, 5]])
);
assert.throws(
  () => core.validatePolygon([[0, 0], [1, 0], [2, 0]]),
  /zero area/
);

const caseContext = { window: {} };
runBrowserScript("data/cases.js", caseContext);
const cases = caseContext.window.RULING_CASES;
assert.equal(cases.length, 2);

for (const item of cases) {
  assert.equal(
    item.counts.total,
    item.counts.allDisk + item.counts.annular,
    `${item.id}: ruling sectors`
  );
  const allDiskMultiplicity = item.rational
    .filter((drawing) => drawing.annularEyes === 0)
    .reduce((sum, drawing) => sum + drawing.multiplicity, 0);
  const annularMultiplicity = item.rational
    .filter((drawing) => drawing.annularEyes > 0)
    .reduce((sum, drawing) => sum + drawing.multiplicity, 0);
  assert.equal(allDiskMultiplicity, item.counts.allDisk);
  assert.equal(annularMultiplicity, item.counts.annular);

  for (const drawing of [item.standard, ...item.rational]) {
    const svgPath = path.join(root, drawing.src);
    assert.ok(fs.existsSync(svgPath), `${drawing.src} exists`);
    const svg = fs.readFileSync(svgPath, "utf8");
    assert.match(svg, /<svg xmlns="http:\/\/www\.w3\.org\/2000\/svg"/);
    assert.match(svg, /fill-opacity:\.30/);
  }
}

const index = fs.readFileSync(path.join(root, "index.html"), "utf8");
assert.match(index, /name="robots" content="noindex, nofollow, noarchive"/);
assert.match(index, /src="core\.js"/);
assert.match(index, /src="data\/cases\.js"/);
assert.match(index, /src="app\.js"/);
assert.doesNotMatch(index, /https?:\/\/[^"]+\.js/);

console.log("ruling-lab validation passed");
