(() => {
  "use strict";

  const { formatVertices, pointLocation } = window.RulingLabCore;
  const cards = [...document.querySelectorAll(".more-polygon-card")];
  const reports = new Map(
    [...document.querySelectorAll(".more-polygon-report")]
      .map((report) => [report.id, report])
  );
  const atlasTitle = document.querySelector("#more-atlas-title");

  for (const card of cards) {
    const vertices = JSON.parse(card.dataset.vertices);
    card.querySelector(".polygon-card-image").innerHTML =
      polygonSvg(vertices);
    card.addEventListener("click", () => {
      selectReport(card.dataset.reportId, true);
    });
  }

  for (const button of document.querySelectorAll(".more-back-to-atlas")) {
    button.addEventListener("click", () => {
      clearSelection(true, true);
    });
  }

  restoreFromHash(false);
  window.addEventListener("hashchange", () => restoreFromHash(false));
  window.addEventListener("popstate", () => restoreFromHash(false));

  function restoreFromHash(shouldScroll) {
    const id = decodeURIComponent(window.location.hash.replace(/^#/, ""));
    if (reports.has(id)) {
      selectReport(id, false, shouldScroll);
      return;
    }
    clearSelection(false, false);
  }

  function selectReport(id, updateHash, shouldScroll = true) {
    const selected = reports.get(id);
    if (!selected) {
      return;
    }
    for (const [reportId, report] of reports) {
      report.hidden = reportId !== id;
    }
    for (const card of cards) {
      const active = card.dataset.reportId === id;
      card.classList.toggle("is-selected", active);
      card.setAttribute("aria-pressed", String(active));
      card.setAttribute("aria-expanded", String(active));
    }
    if (updateHash) {
      history.pushState(null, "", `#${encodeURIComponent(id)}`);
    }
    if (shouldScroll) {
      selected.scrollIntoView({ behavior: "smooth", block: "start" });
      selected.querySelector("h2").focus({ preventScroll: true });
    }
  }

  function clearSelection(updateHistory, shouldScroll) {
    for (const report of reports.values()) {
      report.hidden = true;
    }
    for (const card of cards) {
      card.classList.remove("is-selected");
      card.setAttribute("aria-pressed", "false");
      card.setAttribute("aria-expanded", "false");
    }
    if (updateHistory) {
      history.pushState(
        null,
        "",
        window.location.pathname + window.location.search
      );
    }
    if (shouldScroll) {
      atlasTitle.scrollIntoView({ behavior: "smooth", block: "start" });
      atlasTitle.focus({ preventScroll: true });
    }
  }

  function polygonSvg(vertices) {
    const width = 300;
    const height = 220;
    const padding = 27;
    const xs = vertices.map(([x]) => x);
    const ys = vertices.map(([, y]) => y);
    const minimumX = Math.min(...xs);
    const maximumX = Math.max(...xs);
    const minimumY = Math.min(...ys);
    const maximumY = Math.max(...ys);
    const spanX = Math.max(1, maximumX - minimumX);
    const spanY = Math.max(1, maximumY - minimumY);
    const scale = Math.min(
      (width - 2 * padding) / spanX,
      (height - 2 * padding) / spanY
    );
    const offsetX = (width - spanX * scale) / 2;
    const offsetY = (height - spanY * scale) / 2;
    const project = ([x, y]) => [
      offsetX + (x - minimumX) * scale,
      height - offsetY - (y - minimumY) * scale,
    ];
    const polygonPoints = vertices
      .map((point) => project(point).join(","))
      .join(" ");

    const grid = [];
    for (let x = Math.ceil(minimumX); x <= Math.floor(maximumX); x += 1) {
      const [projectedX] = project([x, minimumY]);
      grid.push(
        `<line x1="${projectedX}" y1="${padding / 2}" ` +
        `x2="${projectedX}" y2="${height - padding / 2}"/>`
      );
    }
    for (let y = Math.ceil(minimumY); y <= Math.floor(maximumY); y += 1) {
      const [, projectedY] = project([minimumX, y]);
      grid.push(
        `<line x1="${padding / 2}" y1="${projectedY}" ` +
        `x2="${width - padding / 2}" y2="${projectedY}"/>`
      );
    }

    const dots = [];
    for (let x = Math.ceil(minimumX); x <= Math.floor(maximumX); x += 1) {
      for (let y = Math.ceil(minimumY); y <= Math.floor(maximumY); y += 1) {
        const location = pointLocation(vertices, [x, y]);
        if (location === "outside") {
          continue;
        }
        const [dotX, dotY] = project([x, y]);
        const interior = location === "interior";
        dots.push(
          `<circle cx="${dotX}" cy="${dotY}" ` +
          `r="${interior ? 5.2 : 3.6}" ` +
          `fill="${interior ? "#9f2f25" : "#17202a"}" ` +
          (interior
            ? 'stroke="#fffdf8" stroke-width="1.8"'
            : "") +
          "/>"
        );
      }
    }

    return (
      `<svg viewBox="0 0 ${width} ${height}" role="img" ` +
      `aria-label="Lattice polygon ${escapeAttribute(
        formatVertices(vertices)
      )}">` +
      `<g stroke="#d9d1c2" stroke-width="1">${grid.join("")}</g>` +
      `<polygon points="${polygonPoints}" fill="rgba(159,47,37,.13)" ` +
      'stroke="#9f2f25" stroke-width="2.6"/>' +
      `<g>${dots.join("")}</g>` +
      "</svg>"
    );
  }

  function escapeAttribute(value) {
    return value
      .replaceAll("&", "&amp;")
      .replaceAll('"', "&quot;")
      .replaceAll("<", "&lt;")
      .replaceAll(">", "&gt;");
  }
})();
