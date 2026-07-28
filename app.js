(() => {
  "use strict";

  const {
    formatArea,
    formatMatrix,
    formatRulingPolynomial,
    formatVertices,
    pointLocation,
    totalDrawingTransform,
    transformVertices,
  } = window.RulingLabCore;
  const cases = Array.isArray(window.RULING_CASES)
    ? window.RULING_CASES
    : [];
  const caseById = new Map(cases.map((item) => [item.id, item]));

  const gallery = document.querySelector("#polygon-gallery");
  const report = document.querySelector("#report");
  const reportTitle = document.querySelector("#report-title");
  const selectedPolygon = document.querySelector("#selected-polygon");
  const selectedVertices = document.querySelector("#selected-vertices");
  const countTableBody = document.querySelector("#count-table-body");
  const rulingPolynomial = document.querySelector("#ruling-polynomial");
  const genusSections = document.querySelector("#genus-sections");
  const cardTemplate = document.querySelector("#polygon-card-template");
  const rulingTemplate = document.querySelector("#ruling-card-template");
  const backButton = document.querySelector("#back-to-atlas");
  const atlasStatus = document.querySelector("#atlas-status");

  let selectedId = null;
  const cardById = new Map();

  renderAtlas();
  restoreFromHash(false);

  window.addEventListener("hashchange", () => restoreFromHash(false));
  backButton.addEventListener("click", () => {
    selectedId = null;
    report.hidden = true;
    markSelectedCard();
    history.pushState(null, "", window.location.pathname + window.location.search);
    document.querySelector("#atlas-title").scrollIntoView({
      behavior: "smooth",
      block: "start",
    });
  });

  function renderAtlas() {
    const fragment = document.createDocumentFragment();
    for (const item of cases) {
      const card = cardTemplate.content.firstElementChild.cloneNode(true);
      const image = card.querySelector(".polygon-card-image");
      const coordinates = card.querySelector(".polygon-card-coordinates");
      const counts = card.querySelector(".polygon-card-counts");
      const displayTransform = totalDrawingTransform(
        Number(item.arrangementShear || 0),
        item.displaySl2z
      );
      const drawingVertices = transformVertices(
        item.vertices,
        displayTransform
      );
      const displayLabel = item.displayName
        ? `${item.displayName}: `
        : "";

      card.dataset.caseId = item.id;
      card.setAttribute(
        "aria-label",
        `Open all rulings for ${displayLabel}${formatVertices(drawingVertices)}`
      );
      image.innerHTML = polygonSvg(drawingVertices, { compact: true });
      coordinates.textContent = [
        item.displayName || null,
        `SL₂(ℤ) representative: ${formatVertices(drawingVertices)}`,
        `vertical direction: ${formatDirection(item.verticalDirection)}`,
      ].filter(Boolean).join("\n");
      counts.textContent = [...item.genera]
        .sort((left, right) => right.genus - left.genus)
        .map((entry) => `g=${entry.genus}: ${entry.counts.total}`)
        .join(" · ");
      card.addEventListener("click", () => selectCase(item.id, true));
      cardById.set(item.id, card);
      fragment.append(card);
    }
    gallery.replaceChildren(fragment);
    atlasStatus.textContent =
      `${cases.length} polygons · ` +
      `${cases.reduce(
        (sum, item) => sum + item.genera.reduce(
          (inner, entry) => inner + Number(entry.counts.total),
          0
        ),
        0
      )} rulings`;
  }

  function restoreFromHash(shouldScroll) {
    const id = decodeURIComponent(window.location.hash.replace(/^#/, ""));
    if (caseById.has(id)) {
      selectCase(id, false, shouldScroll);
    }
  }

  function selectCase(id, updateHash, shouldScroll = true) {
    const item = caseById.get(id);
    if (!item) {
      return;
    }
    selectedId = id;
    markSelectedCard();
    renderReport(item);
    if (updateHash) {
      history.pushState(null, "", `#${encodeURIComponent(id)}`);
    }
    if (shouldScroll) {
      report.scrollIntoView({ behavior: "smooth", block: "start" });
      reportTitle.focus({ preventScroll: true });
    }
  }

  function markSelectedCard() {
    for (const [id, card] of cardById) {
      card.classList.toggle("is-selected", id === selectedId);
      card.setAttribute("aria-pressed", String(id === selectedId));
    }
  }

  function renderReport(item) {
    const fixedTransform = totalDrawingTransform(
      Number(item.arrangementShear || 0),
      item.displaySl2z
    );
    const drawingVertices = transformVertices(item.vertices, fixedTransform);
    const coordinateText = formatVertices(drawingVertices);
    reportTitle.textContent = [
      item.displayName,
      `Polygon used for every ruling: ${coordinateText}`,
    ].filter(Boolean).join(" — ");
    selectedVertices.textContent = [
      `SL₂(ℤ) representative: ${coordinateText}`,
      `catalog representative: ${formatVertices(item.vertices)}`,
      `matrix from catalog coordinates: ${formatMatrix(fixedTransform)}`,
      `vertical direction: ${formatDirection(item.verticalDirection)}`,
    ].join("\n");
    selectedPolygon.innerHTML = polygonSvg(drawingVertices, { compact: false });
    document.querySelector("#polygon-area").textContent =
      formatArea(item.doubleArea);
    document.querySelector("#boundary-count").textContent =
      item.boundaryLatticePoints;
    document.querySelector("#interior-count").textContent =
      item.interiorLatticePoints;

    const sortedGenera = [...item.genera].sort(
      (left, right) => right.genus - left.genus
    );
    countTableBody.replaceChildren(
      ...sortedGenera.map((entry) => {
        const row = document.createElement("tr");
        const genusCell = document.createElement("th");
        genusCell.scope = "row";
        genusCell.textContent = entry.genus;
        row.append(genusCell);
        for (const value of [
          entry.counts.allDisk,
          entry.counts.annular,
          entry.counts.total,
        ]) {
          const cell = document.createElement("td");
          cell.textContent = value;
          row.append(cell);
        }
        return row;
      })
    );
    rulingPolynomial.textContent = formatRulingPolynomial(sortedGenera);

    genusSections.replaceChildren(
      ...sortedGenera.map((entry) => renderGenusPanel(item, entry))
    );
    report.hidden = false;
  }

  function renderGenusPanel(item, entry) {
    const panel = document.createElement("details");
    const summary = document.createElement("summary");
    const heading = document.createElement("span");
    const title = document.createElement("strong");
    const description = document.createElement("span");
    const summaryCounts = document.createElement("span");
    const rulingGallery = document.createElement("div");

    panel.className = "genus-panel";
    heading.className = "genus-heading";
    title.textContent = `Genus ${entry.genus}`;
    description.textContent = "rulings";
    summaryCounts.className = "genus-summary-counts";
    summaryCounts.textContent =
      `${entry.counts.total} total · ` +
      `${entry.counts.allDisk} disk · ${entry.counts.annular} annular`;
    heading.append(title, description);
    summary.append(heading, summaryCounts);
    rulingGallery.className = "ruling-gallery";

    const rulings = Array.isArray(entry.rulings) ? entry.rulings : [];
    if (rulings.length === 0) {
      const note = document.createElement("p");
      note.className = "empty-gallery";
      note.textContent = "Ruling diagrams are being generated for this entry.";
      rulingGallery.append(note);
    } else {
      rulingGallery.append(
        ...rulings.map((ruling, index) =>
          renderRulingCard(item, entry, ruling, index)
        )
      );
    }
    panel.append(summary, rulingGallery);
    return panel;
  }

  function renderRulingCard(item, genusEntry, ruling, index) {
    const figure = rulingTemplate.content.firstElementChild.cloneNode(true);
    const image = figure.querySelector("img");
    const label = figure.querySelector(".diagram-label");
    const profile = figure.querySelector(".diagram-profile");
    const multiplicity = Number(ruling.multiplicity || 1);

    const transform = totalDrawingTransform(
      Number(item.arrangementShear || 0),
      ruling.sl2z
    );
    const drawingVertices = transformVertices(item.vertices, transform);

    image.src = ruling.src;
    image.alt =
      `${formatVertices(drawingVertices)}, genus ${genusEntry.genus}, ` +
      `${ruling.sector || ruling.family || "ruling"} ` +
      `${index + 1} of ${genusEntry.rulings.length}`;
    const sectorLabel = String(ruling.sector || "")
      .toLowerCase()
      .startsWith("annular")
      ? "Annular"
      : "All-disk";
    label.textContent =
      `${sectorLabel} · ${ruling.identifier || `R${index + 1}`}` +
      (multiplicity > 1 ? ` ×${multiplicity}` : "");
    profile.textContent = [
      ruling.profile,
      Number.isInteger(ruling.switches)
        ? `${ruling.switches} switches`
        : null,
    ].filter(Boolean).join(" · ");
    return figure;
  }

  function polygonSvg(vertices, { compact }) {
    const width = compact ? 300 : 440;
    const height = compact ? 220 : 330;
    const padding = compact ? 27 : 38;
    const xs = vertices.map((point) => point[0]);
    const ys = vertices.map((point) => point[1]);
    const minX = Math.min(...xs);
    const maxX = Math.max(...xs);
    const minY = Math.min(...ys);
    const maxY = Math.max(...ys);
    const spanX = Math.max(1, maxX - minX);
    const spanY = Math.max(1, maxY - minY);
    const scale = Math.min(
      (width - 2 * padding) / spanX,
      (height - 2 * padding) / spanY
    );
    const offsetX = (width - spanX * scale) / 2;
    const offsetY = (height - spanY * scale) / 2;
    const project = ([x, y]) => [
      offsetX + (x - minX) * scale,
      height - offsetY - (y - minY) * scale,
    ];
    const polygonPoints = vertices
      .map((point) => project(point).join(","))
      .join(" ");

    const grid = [];
    for (let x = Math.ceil(minX); x <= Math.floor(maxX); x += 1) {
      const [projectedX] = project([x, minY]);
      grid.push(
        `<line x1="${projectedX}" y1="${padding / 2}" ` +
        `x2="${projectedX}" y2="${height - padding / 2}"/>`
      );
    }
    for (let y = Math.ceil(minY); y <= Math.floor(maxY); y += 1) {
      const [, projectedY] = project([minX, y]);
      grid.push(
        `<line x1="${padding / 2}" y1="${projectedY}" ` +
        `x2="${width - padding / 2}" y2="${projectedY}"/>`
      );
    }

    const dots = [];
    for (let x = Math.ceil(minX); x <= Math.floor(maxX); x += 1) {
      for (let y = Math.ceil(minY); y <= Math.floor(maxY); y += 1) {
        const location = pointLocation(vertices, [x, y]);
        if (location === "outside") {
          continue;
        }
        const [dotX, dotY] = project([x, y]);
        const isInterior = location === "interior";
        dots.push(
          `<circle cx="${dotX}" cy="${dotY}" ` +
          `r="${isInterior ? 5.2 : 3.6}" ` +
          `fill="${isInterior ? "#9f2f25" : "#17202a"}" ` +
          (isInterior
            ? 'stroke="#fffdf8" stroke-width="1.8"'
            : "") +
          "/>"
        );
      }
    }

    return (
      `<svg viewBox="0 0 ${width} ${height}" role="img" ` +
      `aria-label="Lattice polygon ${escapeAttribute(formatVertices(vertices))}">` +
      `<g stroke="#d9d1c2" stroke-width="1">${grid.join("")}</g>` +
      `<polygon points="${polygonPoints}" fill="rgba(159,47,37,.13)" ` +
      `stroke="#9f2f25" stroke-width="${compact ? 2.6 : 3.2}"/>` +
      `<g>${dots.join("")}</g>` +
      `</svg>`
    );
  }

  function escapeAttribute(value) {
    return value
      .replaceAll("&", "&amp;")
      .replaceAll('"', "&quot;")
      .replaceAll("<", "&lt;")
      .replaceAll(">", "&gt;");
  }

  function formatDirection(direction) {
    const [x, y] = Array.isArray(direction) ? direction : [0, 1];
    return `(${Number(x)},${Number(y)})`;
  }
})();
