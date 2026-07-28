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

function gcd(left, right) {
  let a = Math.abs(left);
  let b = Math.abs(right);
  while (b !== 0) {
    [a, b] = [b, a % b];
  }
  return a;
}

function expectedArrangementShear(vertices) {
  const directions = [];
  for (let index = 0; index < vertices.length; index += 1) {
    const start = vertices[index];
    const end = vertices[(index + 1) % vertices.length];
    const dx = end[0] - start[0];
    const dy = end[1] - start[1];
    const length = gcd(dx, dy);
    for (let copy = 0; copy < length; copy += 1) {
      directions.push([dx / length, dy / length]);
    }
  }
  for (let magnitude = 0; magnitude <= directions.length; magnitude += 1) {
    const candidates = magnitude === 0 ? [0] : [magnitude, -magnitude];
    for (const shear of candidates) {
      if (directions.every(([dx, dy]) => dx + shear * dy !== 0)) {
        return shear;
      }
    }
  }
  throw new Error("no transverse shear");
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
assert.equal(
  core.formatRulingPolynomial([
    { genus: 0, counts: { total: 8 } },
    { genus: 1, counts: { total: 1 } },
  ]),
  "z² + 8",
  "the ruling polynomial is written from highest genus to lowest"
);
const annularTransform = core.totalDrawingTransform(
  1,
  [[0, 1], [-1, 2]]
);
assert.deepEqual(
  JSON.parse(JSON.stringify(annularTransform)),
  [[0, 1], [-1, 1]],
  "the displayed transform is A times the arrangement shear"
);
assert.deepEqual(
  JSON.parse(JSON.stringify(core.transformVertices(
    [[-3, -2], [1, 0], [1, 2]],
    annularTransform
  ))),
  [[-2, 1], [0, -1], [2, 1]],
  "the displayed SL2-equivalent polygon uses the full drawing transform"
);

const caseContext = { window: {} };
runBrowserScript("data/cases.js", caseContext);
const cases = caseContext.window.RULING_CASES;
assert.equal(cases.length, 16, "the atlas contains all 16 classes");
assert.equal(new Set(cases.map((item) => item.id)).size, 16, "unique ids");

let drawingCount = 0;
let rulingMultiplicity = 0;
let rationalMultiplicity = 0;
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
  assert.equal(
    item.arrangementShear,
    expectedArrangementShear(item.vertices),
    `${item.id}: exported drawing shear`
  );

  const genera = [...item.genera].sort((left, right) => left.genus - right.genus);
  const drawingRepresentatives = new Set();
  assert.deepEqual(
    genera.map((entry) => entry.genus),
    [0, 1],
    `${item.id}: all possible genera`
  );
  assert.equal(genera[1].counts.total, 1, `${item.id}: genus-one ruling`);
  assert.equal(
    core.formatRulingPolynomial(genera),
    `z² + ${genera[0].counts.total}`,
    `${item.id}: ruling polynomial`
  );

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
    if (entry.genus === 0) {
      rationalMultiplicity += multiplicity;
    }

    for (const drawing of entry.rulings) {
      assert.doesNotMatch(
        drawing.identifier,
        /standard/i,
        `${item.id}: ordinary ruling identifier`
      );
      assert.doesNotMatch(
        drawing.sector,
        /standard/i,
        `${item.id}: ordinary ruling sector`
      );
      assert.equal(drawing.sl2z.length, 2, `${item.id}: SL2 matrix rows`);
      const [[a, b], [c, d]] = drawing.sl2z;
      assert.equal(a * d - b * c, 1, `${item.id}: determinant one`);
      assert.deepEqual(
        JSON.parse(JSON.stringify(drawing.sl2z)),
        [[1, 0], [0, 1]],
        `${item.id}: fixed global frame`
      );
      const drawingTransform = core.totalDrawingTransform(
        item.arrangementShear,
        drawing.sl2z
      );
      drawingRepresentatives.add(JSON.stringify(
        core.transformVertices(item.vertices, drawingTransform)
      ));
      const svgPath = path.join(root, drawing.src);
      assert.ok(fs.existsSync(svgPath), `${drawing.src} exists`);
      const svg = fs.readFileSync(svgPath, "utf8");
      assert.match(svg, /<svg\b[^>]*xmlns="http:\/\/www\.w3\.org\/2000\/svg"/);
      assert.match(svg, /fill-opacity:(?:0?\.)30/);
      assert.equal(
        (svg.match(/<circle\b/g) || []).length,
        drawing.switches,
        `${drawing.src}: black switch points`
      );
      assert.match(
        svg,
        /&quot;sl2z&quot;:\[\[1,0\],\[0,1\]\]/,
        `${drawing.src}: identity global frame`
      );
      if (entry.counts.annular > 0 && drawing.sector === "annular") {
        assert.match(svg, /fixed fundamental square/);
      }
      assert.doesNotMatch(svg, /\bDing\b|Fig(?:ure)?\.?\s*\d/i);
    }
  }
  assert.equal(
    drawingRepresentatives.size,
    1,
    `${item.id}: one representative for every ruling`
  );
}

assert.equal(rulingMultiplicity, 112, "96 genus-zero + 16 genus-one rulings");
assert.equal(rationalMultiplicity, 96, "the preview contains 96 genus-zero rulings");
assert.ok(drawingCount > 0 && drawingCount <= rulingMultiplicity);

const index = fs.readFileSync(path.join(root, "index.html"), "utf8");
assert.match(index, /name="robots" content="noindex, nofollow, noarchive"/);
assert.match(index, /href="one-interior\.html"/);
assert.match(index, /Polygons with one interior lattice point/);
assert.match(index, /112 rulings/);
assert.doesNotMatch(index, /https?:\/\/[^"]+\.js/);
assert.doesNotMatch(index, /\bDing\b|Fig(?:ure)?\.?\s*\d/i);
assert.match(index, /Rational Ruling Atlas/);
assert.doesNotMatch(index, /src="(?:core|app)\.js"|src="data\/cases\.js"/);

const oneInterior = fs.readFileSync(path.join(root, "one-interior.html"), "utf8");
assert.match(oneInterior, /src="core\.js"/);
assert.match(oneInterior, /src="data\/cases\.js"/);
assert.match(oneInterior, /src="app\.js"/);
assert.match(oneInterior, /Polygons with one interior lattice point/);
assert.doesNotMatch(oneInterior, /https?:\/\/[^"]+\.js/);
assert.doesNotMatch(oneInterior, /\bDing\b|Fig(?:ure)?\.?\s*\d/i);
assert.match(oneInterior, /<th scope="col">Genus<\/th>/);
assert.match(oneInterior, /genera&nbsp;0 and&nbsp;1/);
assert.match(oneInterior, /Polygon used for every ruling/);
assert.match(oneInterior, /fixed\s+Λ/);
assert.doesNotMatch(oneInterior, /diagram-representative/);
assert.doesNotMatch(oneInterior, /standard ruling/i);
assert.doesNotMatch(oneInterior, /id="standard-ruling"/);
assert.match(oneInterior, /id="genus-sections"/);
assert.match(oneInterior, /id="ruling-polynomial"/);
assert.match(oneInterior, /Switches are the black points/);

const app = fs.readFileSync(path.join(root, "app.js"), "utf8");
assert.doesNotMatch(app, /standard ruling/i);
assert.doesNotMatch(app, /standardEntry|standardRuling/);
assert.match(app, /document\.createElement\("details"\)/);
assert.match(app, /sortedGenera\.map\(\(entry\) => renderGenusPanel\(item, entry\)\)/);
assert.doesNotMatch(app, /genus-links|showGenus/);
assert.equal(
  (app.match(/right\.genus - left\.genus/g) || []).length,
  2,
  "atlas counts and ruling galleries are ordered from high genus to low"
);
assert.match(app, /totalDrawingTransform/);
assert.match(app, /fixedTransform/);
assert.match(app, /formatRulingPolynomial/);
assert.match(app, /Polygon used for every ruling/);

const readme = fs.readFileSync(path.join(root, "README.md"), "utf8");
assert.doesNotMatch(readme, /standard ruling/i);
assert.match(readme, /highest genus to lowest genus/);

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
