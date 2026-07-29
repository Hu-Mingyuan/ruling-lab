import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import vm from "node:vm";
import { fileURLToPath } from "node:url";


const here = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(here, "..");
const cjkText =
  /[\p{Script=Han}\p{Script=Hiragana}\p{Script=Katakana}\p{Script=Hangul}]/u;

for (const filename of [
  "index.html",
  "one-interior.html",
  "two-interior.html",
  "other-polygons.html",
  "app.js",
  "core.js",
  "styles.css",
  "data/cases.js",
]) {
  assert.doesNotMatch(
    fs.readFileSync(path.join(root, filename), "utf8"),
    cjkText,
    `${filename}: public-facing text must not contain CJK characters`
  );
}

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

function uniqueTorusSwitchCount(svg) {
  const frames = [...svg.matchAll(
    /<rect x="(-?\d+(?:\.\d+)?)" y="(-?\d+(?:\.\d+)?)" width="(\d+(?:\.\d+)?)" height="(\d+(?:\.\d+)?)" style="fill:none;stroke:var\(--(?:viz-series-3|lambda)/g
  )];
  assert.ok(frames.length > 0, "the SVG has a torus frame");
  const frame = frames.at(-1);
  const origin = [Number(frame[1]), Number(frame[2])];
  const period = [Number(frame[3]), Number(frame[4])];
  const points = new Set();
  for (const match of svg.matchAll(
    /<circle cx="(-?\d+(?:\.\d+)?)" cy="(-?\d+(?:\.\d+)?)"/g
  )) {
    const normalize = (value, axis) => {
      const offset = Number(value) - origin[axis];
      const reduced = ((offset % period[axis]) + period[axis]) % period[axis];
      return Math.round((origin[axis] + reduced) * 1e6) / 1e6;
    };
    points.add(`${normalize(match[1], 0)},${normalize(match[2], 1)}`);
  }
  return points.size;
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
assert.equal(cases.length, 61, "the atlas contains all 61 classes");
assert.equal(new Set(cases.map((item) => item.id)).size, 61, "unique ids");

const collectionExpectations = new Map([
  [
    "interior-1",
    {
      interiorPoints: 1,
      polygonCount: 16,
      genera: [0, 1],
      rulingMultiplicity: 112,
    },
  ],
  [
    "interior-2",
    {
      interiorPoints: 2,
      polygonCount: 45,
      genera: [0, 1, 2],
      rulingMultiplicity: 1489,
    },
  ],
]);
const collectionStats = new Map(
  [...collectionExpectations].map(([key]) => [
    key,
    { polygons: 0, drawings: 0, multiplicity: 0 },
  ])
);

let drawingCount = 0;
let rulingMultiplicity = 0;
const drawingPaths = new Set();
for (const item of cases) {
  const expectation = collectionExpectations.get(item.collection);
  assert.ok(expectation, `${item.id}: known collection`);
  const stats = collectionStats.get(item.collection);
  stats.polygons += 1;

  assert.match(item.id, /^polygon-[0-9a-f]{10}$/);
  assert.ok(item.vertices.length >= 3);
  assert.ok(core.signedDoubleArea(item.vertices) > 0, `${item.id}: CCW`);
  assert.equal(
    item.interiorLatticePoints,
    expectation.interiorPoints,
    `${item.id}: stored interior-point count`
  );
  const interiorPoints = enumerateInteriorPoints(core, item.vertices);
  assert.equal(
    interiorPoints.length,
    expectation.interiorPoints,
    `${item.id}: enumerated interior-point count`
  );
  if (item.collection === "interior-1") {
    assert.deepEqual(
      interiorPoints,
      [[0, 0]],
      `${item.id}: unique interior point`
    );
  }
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
  assert.equal(item.displaySl2z.length, 2, `${item.id}: display matrix rows`);
  const [[displayA, displayB], [displayC, displayD]] = item.displaySl2z;
  assert.equal(
    displayA * displayD - displayB * displayC,
    1,
    `${item.id}: display matrix determinant`
  );
  assert.deepEqual(
    JSON.parse(JSON.stringify(item.displayFromSource)),
    JSON.parse(JSON.stringify(core.totalDrawingTransform(
      item.arrangementShear,
      item.displaySl2z
    ))),
    `${item.id}: count and display frames compose`
  );
  const expectedVerticalDirection = [displayB, displayD];
  assert.deepEqual(
    JSON.parse(JSON.stringify(item.verticalDirection)),
    expectedVerticalDirection,
    `${item.id}: vertical direction is transported to the display frame`
  );
  assert.deepEqual(
    JSON.parse(JSON.stringify(item.sourceVerticalDirection)),
    [item.arrangementShear === 0 ? 0 : -item.arrangementShear, 1],
    `${item.id}: source-frame vertical direction`
  );
  assert.deepEqual(
    JSON.parse(JSON.stringify(item.countVerticalDirection)),
    [0, 1],
    `${item.id}: count-frame vertical direction`
  );
  assert.deepEqual(
    JSON.parse(JSON.stringify(item.sourceSweepCovector)),
    [1, item.arrangementShear],
    `${item.id}: source-frame sweep covector`
  );
  assert.deepEqual(
    JSON.parse(JSON.stringify(item.countSweepCovector)),
    [1, 0],
    `${item.id}: count-frame sweep covector`
  );
  const [
    [displaySourceA, displaySourceB],
    [displaySourceC, displaySourceD],
  ] = item.displayFromSource;
  assert.deepEqual(
    JSON.parse(JSON.stringify([
      displaySourceA * item.sourceVerticalDirection[0]
        + displaySourceB * item.sourceVerticalDirection[1],
      displaySourceC * item.sourceVerticalDirection[0]
        + displaySourceD * item.sourceVerticalDirection[1],
    ])),
    JSON.parse(JSON.stringify(item.verticalDirection)),
    `${item.id}: source vertical direction reaches the displayed direction`
  );
  const displayedVertices = core.transformVertices(
    item.vertices,
    item.displayFromSource
  );
  for (let index = 0; index < displayedVertices.length; index += 1) {
    const start = displayedVertices[index];
    const finish = displayedVertices[(index + 1) % displayedVertices.length];
    const edge = [finish[0] - start[0], finish[1] - start[1]];
    assert.notEqual(
      item.verticalDirection[0] * edge[1]
        - item.verticalDirection[1] * edge[0],
      0,
      `${item.id}: vertical direction is transverse to every Lambda family`
    );
  }

  const genera = [...item.genera].sort((left, right) => left.genus - right.genus);
  const drawingRepresentatives = new Set();
  assert.deepEqual(
    genera.map((entry) => entry.genus),
    expectation.genera,
    `${item.id}: all possible genera`
  );
  assert.equal(
    genera.at(-1).counts.total,
    1,
    `${item.id}: top-genus ruling`
  );
  if (item.collection === "interior-1") {
    assert.equal(
      core.formatRulingPolynomial(genera),
      `z² + ${genera[0].counts.total}`,
      `${item.id}: ruling polynomial`
    );
  } else {
    assert.equal(
      core.formatRulingPolynomial(genera),
      `z⁴ + ${genera[1].counts.total}z² + ${genera[0].counts.total}`,
      `${item.id}: ruling polynomial`
    );
  }

  for (const entry of genera) {
    assert.equal(
      entry.counts.total,
      entry.counts.allDisk + entry.counts.annular,
      `${item.id}, genus ${entry.genus}: sector sum`
    );
    assert.ok(entry.rulings.length > 0, `${item.id}: drawings generated`);
    const weightedSectors = { "all-disk": 0, annular: 0 };
    const multiplicity = entry.rulings.reduce((sum, drawing) => {
      const weight = Number(drawing.multiplicity || 1);
      assert.equal(
        weight,
        1,
        `${item.id}: every ruling has its own SVG`
      );
      assert.ok(
        Object.hasOwn(weightedSectors, drawing.sector),
        `${item.id}: known ruling sector`
      );
      weightedSectors[drawing.sector] += weight;
      return sum + weight;
    }, 0);
    assert.equal(
      multiplicity,
      entry.counts.total,
      `${item.id}, genus ${entry.genus}: drawing multiplicities`
    );
    assert.equal(
      weightedSectors["all-disk"],
      entry.counts.allDisk,
      `${item.id}, genus ${entry.genus}: weighted all-disk drawings`
    );
    assert.equal(
      weightedSectors.annular,
      entry.counts.annular,
      `${item.id}, genus ${entry.genus}: weighted annular drawings`
    );
    drawingCount += entry.rulings.length;
    stats.drawings += entry.rulings.length;
    rulingMultiplicity += multiplicity;
    stats.multiplicity += multiplicity;

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
        JSON.parse(JSON.stringify(item.displaySl2z)),
        `${item.id}: shared polygon display frame`
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
      assert.ok(
        !drawingPaths.has(drawing.src),
        `${drawing.src}: one card per SVG representative`
      );
      drawingPaths.add(drawing.src);
      const svg = fs.readFileSync(svgPath, "utf8");
      assert.doesNotMatch(svg, cjkText, `${drawing.src}: English-only SVG`);
      assert.match(svg, /<svg\b[^>]*xmlns="http:\/\/www\.w3\.org\/2000\/svg"/);
      assert.match(
        svg,
        /&quot;uses_tropical_count&quot;:false/,
        `${drawing.src}: direct-ruling provenance`
      );
      assert.match(
        svg,
        /fill-opacity:0?\.(?:28|30)/,
        `${drawing.src}: translucent eye fill`
      );
      assert.equal(
        uniqueTorusSwitchCount(svg),
        drawing.switches,
        `${drawing.src}: unique black switch points on the torus`
      );
      assert.ok(
        svg.includes(
          `&quot;sl2z&quot;:${JSON.stringify(drawing.sl2z)}`
        ),
        `${drawing.src}: embedded shared display frame`
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

for (const [collection, expectation] of collectionExpectations) {
  const stats = collectionStats.get(collection);
  assert.equal(
    stats.polygons,
    expectation.polygonCount,
    `${collection}: polygon count`
  );
  assert.equal(
    stats.multiplicity,
    expectation.rulingMultiplicity,
    `${collection}: weighted ruling count`
  );
  assert.equal(
    stats.drawings,
    stats.multiplicity,
    `${collection}: one SVG for every ruling`
  );
}

const planeCubic = cases.find((item) => item.id === "polygon-a2cd47c761");
assert.ok(planeCubic, "the O(3) case is present");
assert.equal(planeCubic.displayName, "O(3) on P^2");
assert.deepEqual(
  JSON.parse(JSON.stringify(planeCubic.displaySl2z)),
  [[2, -3], [-1, 2]]
);
assert.deepEqual(
  JSON.parse(JSON.stringify(planeCubic.displayFromSource)),
  [[2, -1], [-1, 1]]
);
assert.deepEqual(
  JSON.parse(JSON.stringify(planeCubic.verticalDirection)),
  [-3, 2]
);
assert.deepEqual(
  JSON.parse(JSON.stringify(core.transformVertices(
    planeCubic.vertices,
    planeCubic.displayFromSource
  ))),
  [[-1, -1], [2, -1], [-1, 2]],
  "O(3) uses the standard centered triangle with a vertical edge"
);
const planeCubicGenusZero = planeCubic.genera.find(
  (entry) => entry.genus === 0
);
const planeCubicGenusOne = planeCubic.genera.find(
  (entry) => entry.genus === 1
);
assert.deepEqual(
  JSON.parse(JSON.stringify(planeCubicGenusZero.counts)),
  { allDisk: 9, annular: 0, total: 9 },
  "the symmetric O(3) realization has nine all-disk rational rulings"
);
assert.equal(planeCubicGenusOne.rulings.length, 1);
assert.equal(planeCubicGenusOne.rulings[0].mask, "0x7ffffff");
assert.equal(planeCubicGenusOne.rulings[0].diskEyes, 18);
assert.equal(planeCubicGenusOne.rulings[0].switches, 27);
assert.equal(
  core.formatRulingPolynomial(planeCubic.genera),
  "z² + 9"
);
const planeCubicTopSvg = fs.readFileSync(
  path.join(root, planeCubicGenusOne.rulings[0].src),
  "utf8"
);
assert.match(planeCubicTopSvg, /18 equal triangular eyes/);

const o22 = cases.find((item) => item.id === "polygon-10c88eafad");
assert.ok(o22, "the O(2,2) case is present");
assert.equal(o22.displayName, "O(2,2) on P^1 x P^1");
assert.deepEqual(
  JSON.parse(JSON.stringify(o22.displaySl2z)),
  [[1, -1], [-1, 2]]
);
assert.deepEqual(
  JSON.parse(JSON.stringify(o22.displayFromSource)),
  [[1, 0], [-1, 1]]
);
assert.deepEqual(
  JSON.parse(JSON.stringify(o22.verticalDirection)),
  [-1, 2]
);
assert.deepEqual(
  JSON.parse(JSON.stringify(core.transformVertices(
    o22.vertices,
    o22.displayFromSource
  ))),
  [[-1, -1], [1, -1], [1, 1], [-1, 1]],
  "O(2,2) uses the standard centered square"
);
const o22GenusZero = o22.genera.find((entry) => entry.genus === 0);
const o22GenusOne = o22.genera.find((entry) => entry.genus === 1);
assert.deepEqual(
  JSON.parse(JSON.stringify(o22GenusZero.counts)),
  { allDisk: 8, annular: 0, total: 8 }
);
assert.equal(o22GenusOne.rulings[0].mask, "0xffff");
assert.equal(o22GenusOne.rulings[0].diskEyes, 8);
assert.equal(o22GenusOne.rulings[0].switches, 16);
const o22TopSvg = fs.readFileSync(
  path.join(root, o22GenusOne.rulings[0].src),
  "utf8"
);
assert.match(o22TopSvg, /eight equal square eyes/);

const fig2c = cases.find((item) => item.id === "polygon-dcf2998504");
assert.ok(fig2c, "the symmetric two-interior triangle is present");
assert.deepEqual(
  JSON.parse(JSON.stringify(fig2c.displayFromSource)),
  [[1, 0], [0, 1]]
);
assert.deepEqual(
  JSON.parse(JSON.stringify(fig2c.verticalDirection)),
  [0, 1]
);
const fig2cGenusZero = fig2c.genera.find((entry) => entry.genus === 0);
const fig2cGenusTwo = fig2c.genera.find((entry) => entry.genus === 2);
assert.deepEqual(
  JSON.parse(JSON.stringify(fig2cGenusZero.counts)),
  { allDisk: 8, annular: 8, total: 16 }
);
assert.equal(fig2cGenusTwo.rulings[0].mask, "0xffffff");
assert.equal(fig2cGenusTwo.rulings[0].diskEyes, 16);
assert.equal(fig2cGenusTwo.rulings[0].switches, 24);

const o23 = cases.find((item) => item.id === "polygon-8ee8c76d7d");
assert.ok(o23, "the O(2,3) rectangle is present");
assert.equal(o23.displayName, "O(2,3) on P^1 x P^1");
assert.deepEqual(
  JSON.parse(JSON.stringify(o23.displayFromSource)),
  [[1, 0], [0, 1]]
);
assert.deepEqual(
  JSON.parse(JSON.stringify(o23.verticalDirection)),
  [-1, 1]
);
assert.deepEqual(
  JSON.parse(JSON.stringify(core.transformVertices(
    o23.vertices,
    o23.displayFromSource
  ))),
  [[0, 0], [3, 0], [3, 2], [0, 2]],
  "O(2,3) keeps the catalog rectangle"
);
const o23GenusZero = o23.genera.find((entry) => entry.genus === 0);
const o23GenusOne = o23.genera.find((entry) => entry.genus === 1);
const o23GenusTwo = o23.genera.find((entry) => entry.genus === 2);
assert.deepEqual(
  JSON.parse(JSON.stringify(o23GenusZero.counts)),
  { allDisk: 48, annular: 0, total: 48 }
);
assert.deepEqual(
  JSON.parse(JSON.stringify(o23GenusOne.counts)),
  { allDisk: 12, annular: 0, total: 12 }
);
assert.equal(o23GenusTwo.rulings[0].mask, "0xffffff");
assert.equal(o23GenusTwo.rulings[0].diskEyes, 12);
assert.equal(o23GenusTwo.rulings[0].switches, 24);
const o23TopSvg = fs.readFileSync(
  path.join(root, o23GenusTwo.rulings[0].src),
  "utf8"
);
assert.match(o23TopSvg, /6-by-4 rectangular Lambda/);

const separatedHexagon = cases.find(
  (item) => item.id === "polygon-e95d7ae246"
);
assert.ok(separatedHexagon, "the separated hexagon is present");
assert.deepEqual(
  JSON.parse(JSON.stringify(separatedHexagon.displayFromSource)),
  [[1, 0], [0, 1]]
);
assert.deepEqual(
  JSON.parse(JSON.stringify(separatedHexagon.verticalDirection)),
  [1, 1]
);
const hexagonGenusZero = separatedHexagon.genera.find(
  (entry) => entry.genus === 0
);
const hexagonGenusOne = separatedHexagon.genera.find(
  (entry) => entry.genus === 1
);
assert.deepEqual(
  JSON.parse(JSON.stringify(hexagonGenusZero.counts)),
  { allDisk: 5, annular: 1, total: 6 },
  "the separated display has five disk and one annular rational ruling"
);
assert.equal(
  hexagonGenusZero.rulings.filter(
    (drawing) => drawing.sector === "annular"
  ).length,
  1
);
assert.equal(hexagonGenusOne.rulings[0].mask, "0x0bf7");
assert.equal(hexagonGenusOne.rulings[0].switches, 10);

const symmetricPentagon = cases.find(
  (item) => item.id === "polygon-6a605367ef"
);
assert.ok(symmetricPentagon, "the reflection-symmetric pentagon is present");
assert.deepEqual(
  JSON.parse(JSON.stringify(symmetricPentagon.displayFromSource)),
  [[1, 0], [0, 1]]
);
assert.deepEqual(
  JSON.parse(JSON.stringify(symmetricPentagon.verticalDirection)),
  [1, 1]
);
const pentagonGenusOne = symmetricPentagon.genera.find(
  (entry) => entry.genus === 1
);
assert.equal(pentagonGenusOne.rulings[0].mask, "0x7ce7");
assert.equal(pentagonGenusOne.rulings[0].diskEyes, 4);
assert.equal(pentagonGenusOne.rulings[0].switches, 11);

assert.equal(rulingMultiplicity, 1601, "112 one-point + 1489 two-point rulings");
assert.equal(drawingCount, 1601, "all 1601 rulings are displayed separately");
assert.equal(
  drawingPaths.size,
  drawingCount,
  "one file per displayed SVG representative"
);

const index = fs.readFileSync(path.join(root, "index.html"), "utf8");
assert.match(index, /name="robots" content="noindex, nofollow, noarchive"/);
assert.match(index, /href="one-interior\.html"/);
assert.match(index, /href="two-interior\.html"/);
assert.match(index, /href="other-polygons\.html"/);
assert.match(index, /Polygons with one interior lattice point/);
assert.match(index, /Polygons with two interior lattice points/);
assert.match(index, /Some other polygons/);
assert.match(index, /112 rulings/);
assert.match(index, /1,489 rulings/);
assert.match(index, /2 polygons/);
assert.match(index, /41 saved figures/);
assert.match(
  index,
  /general convex lattice polygon can be computed with\s+the <a href="downloads\/direct-ruling-counter\.zip">accompanying code<\/a>/
);
assert.doesNotMatch(index, /https?:\/\/[^"]+\.js/);
assert.doesNotMatch(index, /\bDing\b|Fig(?:ure)?\.?\s*\d/i);
assert.match(index, /Ruling Atlas/);
assert.doesNotMatch(index, /Rational Ruling Atlas/);
assert.doesNotMatch(index, /src="(?:core|app)\.js"|src="data\/cases\.js"/);

const oneInterior = fs.readFileSync(path.join(root, "one-interior.html"), "utf8");
const twoInterior = fs.readFileSync(path.join(root, "two-interior.html"), "utf8");
const otherPolygons = fs.readFileSync(
  path.join(root, "other-polygons.html"),
  "utf8"
);
for (const [filename, document, collection] of [
  ["one-interior.html", oneInterior, "interior-1"],
  ["two-interior.html", twoInterior, "interior-2"],
]) {
  assert.match(
    document,
    /name="robots" content="noindex, nofollow, noarchive"/,
    `${filename}: noindex`
  );
  assert.match(
    document,
    /<title>[^<]*Ruling Atlas<\/title>/,
    `${filename}: site title`
  );
  assert.doesNotMatch(document, /Rational Ruling Atlas/);
  assert.match(document, new RegExp(`data-collection="${collection}"`));
  assert.match(document, /src="core\.js"/);
  assert.match(document, /src="data\/cases\.js"/);
  assert.match(document, /src="app\.js"/);
  assert.doesNotMatch(document, /https?:\/\/[^"]+\.js/);
  assert.doesNotMatch(document, /\bDing\b|Fig(?:ure)?\.?\s*\d/i);
  assert.match(document, /<th scope="col">Genus<\/th>/);
  assert.match(document, /Polygon used for every ruling/);
  assert.match(document, /fixed\s+Λ/);
  assert.doesNotMatch(document, /diagram-representative/);
  assert.doesNotMatch(document, /standard ruling/i);
  assert.doesNotMatch(document, /id="standard-ruling"/);
  assert.match(document, /id="genus-sections"/);
  assert.match(document, /id="ruling-polynomial"/);
  assert.match(document, /Switches are the black points/);
  assert.match(document, /id="vertical-direction-diagram"/);
  assert.match(document, /class="vertical-direction-arrow"/);
  assert.match(document, /id="vertical-direction-label"/);
  assert.match(document, />vertical\s+v = \(0,1\)<\/span>/);
}
assert.match(oneInterior, /Polygons with one interior lattice point/);
assert.match(oneInterior, /genera&nbsp;0 and&nbsp;1/);
assert.match(twoInterior, /Polygons with two interior lattice points/);
assert.match(twoInterior, /genera&nbsp;0, 1, and&nbsp;2/);
assert.match(twoInterior, /45 polygons\s+·\s+1,489 rulings/);
assert.match(otherPolygons, /Other polygons/);
assert.match(otherPolygons, /O\(4\) on P\^2/);
assert.match(otherPolygons, /O\(2,4\) on P\^1 × P\^1/);
assert.match(otherPolygons, /symmetry orbit/i);
assert.match(otherPolygons, /Genus 3/);
assert.match(otherPolygons, /Genus 2/);
assert.match(otherPolygons, /Genus 1/);
assert.match(otherPolygons, /Genus 0/);
assert.match(otherPolygons, /16 rulings/);
assert.match(otherPolygons, /104 rulings/);
assert.match(otherPolygons, /304 rational rulings/);
assert.match(otherPolygons, /256 rational rulings/);
assert.match(otherPolygons, /Conv\{\(-2,-1\), \(2,-1\), \(-2,3\)\}/);
assert.match(otherPolygons, /vertical direction: \(-1,2\)/);
assert.match(otherPolygons, /z⁶ \+ 16z⁴ \+ 104z² \+ 304/);
assert.match(
  otherPolygons,
  /Conv\{\(0,0\), \(4,0\), \(4,2\), \(0,2\)\}/
);
assert.match(otherPolygons, /source-to-count matrix: \[\[1,-1\], \[1,0\]\]/);
assert.match(otherPolygons, /vertical direction: \(1,1\)/);
assert.match(otherPolygons, /z⁶ \+ 16z⁴ \+ 96z² \+ 256/);
assert.match(
  otherPolygons,
  /G2D01 · symmetry orbit ×16<\/span><span class="diagram-profile">D28 · 42 switches · A\/B\/V 14\/0\/28/
);

const o4FigureDirectory = path.join(
  root,
  "assets",
  "rulings",
  "other",
  "o4-p2"
);
const o4GenusTwoFigures = [
  ["O4_genus2_D01_x16.svg", 16],
];
const o4GenusOneFigures = [
  ["O4_genus1_D01_x16.svg", 16],
  ["O4_genus1_D02_x32.svg", 32],
  ["O4_genus1_D03_x16.svg", 16],
  ["O4_genus1_D04_x16.svg", 16],
  ["O4_genus1_D05_x16.svg", 16],
  ["O4_genus1_D06_x8.svg", 8],
];
const o4GenusZeroFigures = [
  ["O4_rational_D01_x16.svg", 16],
  ["O4_rational_D02_x32.svg", 32],
  ["O4_rational_D03_x32.svg", 32],
  ["O4_rational_D04_x32.svg", 32],
  ["O4_rational_D05_x16.svg", 16],
  ["O4_rational_D06_x32.svg", 32],
  ["O4_rational_D07_x16.svg", 16],
  ["O4_rational_D08_x32.svg", 32],
  ["O4_rational_D09_x32.svg", 32],
  ["O4_rational_D10_x16.svg", 16],
  ["O4_rational_D11_x16.svg", 16],
  ["O4_rational_D12_x16.svg", 16],
  ["O4_rational_R01_x16.svg", 16],
];
const o4Figures = [
  ["O4_standard_bipartite_ruling.svg", 1],
  ...o4GenusTwoFigures,
  ...o4GenusOneFigures,
  ...o4GenusZeroFigures,
];
assert.deepEqual(
  fs.readdirSync(o4FigureDirectory).sort(),
  ["manifest.json", ...o4Figures.map(([filename]) => filename)].sort(),
  "the O(4) page contains exactly 21 SVG symmetry representatives and its manifest"
);
for (const [filename, multiplicity] of o4Figures) {
  const figurePath = path.join(o4FigureDirectory, filename);
  const svg = fs.readFileSync(figurePath, "utf8");
  assert.doesNotMatch(svg, cjkText, `${filename}: English-only SVG`);
  assert.match(svg, /<svg\b[^>]*xmlns="http:\/\/www\.w3\.org\/2000\/svg"/);
  assert.ok(
    otherPolygons.includes(
      `assets/rulings/other/o4-p2/${filename}`
    ),
    `${filename}: linked from Other polygons`
  );
  assert.match(
    otherPolygons,
    new RegExp(`symmetry orbit ×${multiplicity}`)
  );
}
assert.equal(
  o4GenusTwoFigures.reduce(
    (sum, [, multiplicity]) => sum + multiplicity,
    0
  ),
  16,
  "the genus-two symmetry orbit represents all 16 rulings"
);
assert.equal(
  o4GenusOneFigures.reduce(
    (sum, [, multiplicity]) => sum + multiplicity,
    0
  ),
  104,
  "the six genus-one symmetry orbits represent all 104 rulings"
);
assert.equal(
  o4GenusZeroFigures.reduce(
    (sum, [, multiplicity]) => sum + multiplicity,
    0
  ),
  304,
  "the 13 rational symmetry orbits represent all 304 genus-zero rulings"
);
const o4Manifest = JSON.parse(
  fs.readFileSync(path.join(o4FigureDirectory, "manifest.json"), "utf8")
);
assert.equal(o4Manifest.uses_tropical_count, false);
assert.deepEqual(o4Manifest.vertical_direction, [-1, 2]);
assert.deepEqual(o4Manifest.sweep_covector, [2, 1]);
assert.equal(o4Manifest.counts.genus_2_all_disk, 16);
assert.equal(o4Manifest.counts.genus_2_annular, 0);
assert.equal(o4Manifest.counts.genus_1_all_disk, 104);
assert.equal(o4Manifest.counts.genus_1_annular, 0);
assert.equal(o4Manifest.counts.genus_0_all_disk, 288);
assert.equal(o4Manifest.counts.genus_0_annular, 16);
assert.equal(o4Manifest.genus_2_annular_audit.strict_complete, true);
assert.equal(o4Manifest.genus_2_annular_audit.phase_search_exhaustive, true);
assert.equal(o4Manifest.genus_2_annular_audit.finite_phase_cutoff_used, false);
assert.equal(o4Manifest.genus_2_annular_audit.translation_symmetry, false);
assert.equal(o4Manifest.genus_2_annular_audit.max_nodes, null);
assert.equal(o4Manifest.genus_1_annular_audit.strict_complete, true);
assert.equal(o4Manifest.genus_1_annular_audit.phase_search_exhaustive, true);
assert.equal(o4Manifest.genus_1_annular_audit.finite_phase_cutoff_used, false);
assert.equal(o4Manifest.drawings.length, 21);

const o24FigureDirectory = path.join(
  root,
  "assets",
  "rulings",
  "other",
  "o24-p1xp1"
);
const o24GenusThreeFigures = [
  ["O24_genus3_G3D01_x1.svg", 1],
];
const o24GenusTwoFigures = [
  ["O24_genus2_G2D01_x16.svg", 16],
];
const o24GenusOneFigures = [
  ["O24_genus1_G1D01_x16.svg", 16],
  ["O24_genus1_G1D02_x16.svg", 16],
  ["O24_genus1_G1D03_x8.svg", 8],
  ["O24_genus1_G1D04_x16.svg", 16],
  ["O24_genus1_G1D05_x16.svg", 16],
  ["O24_genus1_G1D06_x16.svg", 16],
  ["O24_genus1_G1D07_x8.svg", 8],
];
const o24GenusZeroFigures = [
  ["O24_rational_D01_x16.svg", 16],
  ["O24_rational_D02_x32.svg", 32],
  ["O24_rational_D03_x16.svg", 16],
  ["O24_rational_D04_x32.svg", 32],
  ["O24_rational_D05_x32.svg", 32],
  ["O24_rational_D06_x16.svg", 16],
  ["O24_rational_D07_x32.svg", 32],
  ["O24_rational_D08_x32.svg", 32],
  ["O24_rational_D09_x16.svg", 16],
  ["O24_rational_D10_x16.svg", 16],
  ["O24_rational_D11_x16.svg", 16],
];
const o24Figures = [
  ...o24GenusThreeFigures,
  ...o24GenusTwoFigures,
  ...o24GenusOneFigures,
  ...o24GenusZeroFigures,
];
assert.deepEqual(
  fs.readdirSync(o24FigureDirectory).sort(),
  ["manifest.json", ...o24Figures.map(([filename]) => filename)].sort(),
  "the O(2,4) page contains exactly 20 SVG symmetry representatives and its manifest"
);
const o24Manifest = JSON.parse(
  fs.readFileSync(path.join(o24FigureDirectory, "manifest.json"), "utf8")
);
const o24DrawingByFilename = new Map(
  o24Manifest.drawings.map((drawing) => [drawing.filename, drawing])
);
for (const [filename, multiplicity] of o24Figures) {
  const figurePath = path.join(o24FigureDirectory, filename);
  const svg = fs.readFileSync(figurePath, "utf8");
  const drawing = o24DrawingByFilename.get(filename);
  assert.ok(drawing, `${filename}: listed in the O(2,4) manifest`);
  assert.doesNotMatch(svg, cjkText, `${filename}: English-only SVG`);
  assert.match(svg, /<svg\b[^>]*xmlns="http:\/\/www\.w3\.org\/2000\/svg"/);
  assert.match(svg, /viewBox="0 0 332\.0 180\.0"/);
  assert.match(svg, /uses_tropical_count[^<]*false/);
  assert.equal(drawing.orbit_size, multiplicity);
  assert.equal(drawing.sector, "all-disk");
  assert.equal(
    uniqueTorusSwitchCount(svg),
    drawing.switches,
    `${filename}: unique black switch count`
  );
  assert.ok(
    otherPolygons.includes(
      `assets/rulings/other/o24-p1xp1/${filename}`
    ),
    `${filename}: linked from Other polygons`
  );
  assert.match(
    otherPolygons,
    new RegExp(`symmetry orbit ×${multiplicity}`)
  );
}
for (const [figures, expected, description] of [
  [o24GenusThreeFigures, 1, "genus-three"],
  [o24GenusTwoFigures, 16, "genus-two"],
  [o24GenusOneFigures, 96, "genus-one"],
  [o24GenusZeroFigures, 256, "genus-zero"],
]) {
  assert.equal(
    figures.reduce((sum, [, multiplicity]) => sum + multiplicity, 0),
    expected,
    `the O(2,4) ${description} symmetry orbits are exhaustive`
  );
}
assert.equal(o24Manifest.uses_tropical_count, false);
assert.equal(o24Manifest.uses_preknown_answers, false);
assert.deepEqual(o24Manifest.polygon, [[0, 0], [4, 0], [4, 2], [0, 2]]);
assert.deepEqual(o24Manifest.count_from_source, [[1, -1], [1, 0]]);
assert.deepEqual(o24Manifest.source_from_count, [[0, 1], [-1, 1]]);
assert.deepEqual(o24Manifest.vertical_direction_source, [1, 1]);
assert.deepEqual(o24Manifest.vertical_direction_count, [0, 1]);
assert.deepEqual(o24Manifest.sweep_covector_source, [1, -1]);
assert.deepEqual(o24Manifest.sweep_covector_count, [1, 0]);
assert.equal(o24Manifest.symmetry_group.order, 32);
assert.equal(o24Manifest.counts.genus_3_all_disk, 1);
assert.equal(o24Manifest.counts.genus_3_annular, 0);
assert.equal(o24Manifest.counts.genus_2_all_disk, 16);
assert.equal(o24Manifest.counts.genus_2_annular, 0);
assert.equal(o24Manifest.counts.genus_1_all_disk, 96);
assert.equal(o24Manifest.counts.genus_1_annular, 0);
assert.equal(o24Manifest.counts.genus_0_all_disk, 256);
assert.equal(o24Manifest.counts.genus_0_annular, 0);
assert.equal(o24Manifest.drawings.length, 20);
assert.equal(
  o24Manifest.drawings.filter((drawing) => drawing.genus === 3).length,
  1
);
assert.equal(
  o24Manifest.drawings.filter((drawing) => drawing.genus === 2).length,
  1
);
assert.equal(
  o24Manifest.drawings.filter((drawing) => drawing.genus === 1).length,
  7
);
assert.equal(
  o24Manifest.drawings.filter((drawing) => drawing.genus === 0).length,
  11
);
for (const genus of [0, 1, 2, 3]) {
  const audit = o24Manifest.exact_annular_audits[`genus_${genus}`];
  assert.equal(audit.strict_complete, true);
  assert.equal(audit.phase_search_exhaustive, true);
  assert.equal(audit.finite_phase_cutoff_used, false);
  assert.equal(audit.phase_cardinality, "FINITE");
  assert.equal(audit.exact_annular_ruling_count, 0);
  assert.ok(audit.search_nodes > 0);
  assert.ok(audit.exact_phase_solver_calls > 0);
  assert.equal(audit.empty_phase_systems, audit.exact_phase_solver_calls);
}

const app = fs.readFileSync(path.join(root, "app.js"), "utf8");
assert.doesNotMatch(app, /standard ruling/i);
assert.doesNotMatch(app, /standardEntry|standardRuling/);
assert.match(app, /dataset\.collection/);
assert.match(app, /item\.collection/);
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
assert.match(app, /SL₂\(ℤ\) representative/);
assert.match(app, /matrix from catalog coordinates/);
assert.match(app, /vertical direction/);
assert.match(app, /item\.displaySl2z/);
assert.match(app, /renderDirectionIndicator\(item\.verticalDirection\)/);
assert.match(app, /Math\.atan2\(x, y\)/);

const styles = fs.readFileSync(path.join(root, "styles.css"), "utf8");
assert.match(styles, /\.selected-representative/);
assert.match(styles, /\.vertical-direction-arrow::before/);
assert.match(styles, /rotate\(var\(--direction-angle\)\)/);

const readme = fs.readFileSync(path.join(root, "README.md"), "utf8");
assert.doesNotMatch(readme, /standard ruling/i);
assert.match(readme, /highest genus to lowest genus/);
assert.match(readme, /O\(3\) on P\^2/);
assert.match(readme, /18 equal triangular eyes/);
assert.match(readme, /O\(2,4\) on P\^1 x P\^1/);
assert.match(readme, /zero annular rulings/);

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
