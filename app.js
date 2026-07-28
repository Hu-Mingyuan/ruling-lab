(() => {
  "use strict";

  const {
    parseVertices,
    polygonKey,
    validatePolygon,
  } = window.RulingLabCore;
  const cases = Array.isArray(window.RULING_CASES)
    ? window.RULING_CASES
    : [];
  const caseById = new Map(cases.map((item) => [item.id, item]));
  const caseByKey = new Map(
    cases.map((item) => [polygonKey(item.vertices), item])
  );

  const input = document.querySelector("#polygon-input");
  const analyzeButton = document.querySelector("#analyze-button");
  const message = document.querySelector("#input-message");
  const result = document.querySelector("#result");
  const exampleButtons = document.querySelectorAll("[data-case]");

  analyzeButton.addEventListener("click", analyze);
  input.addEventListener("keydown", (event) => {
    if ((event.ctrlKey || event.metaKey) && event.key === "Enter") {
      analyze();
    }
  });
  exampleButtons.forEach((button) => {
    button.addEventListener("click", () => {
      const item = caseById.get(button.dataset.case);
      if (!item) {
        return;
      }
      input.value = JSON.stringify(item.vertices);
      analyze();
    });
  });

  function analyze() {
    hideMessage();
    let vertices;
    try {
      vertices = parseVertices(input.value);
      validatePolygon(vertices);
    } catch (error) {
      result.hidden = true;
      showMessage(error.message);
      return;
    }

    const item = caseByKey.get(polygonKey(vertices));
    if (!item) {
      result.hidden = true;
      showMessage(
        "This first private version currently recognizes two verified cases. " +
          "The arbitrary-polygon ruling engine will be connected in the next stage."
      );
      return;
    }
    renderCase(item, vertices);
  }

  function renderCase(item, enteredVertices) {
    document.querySelector("#case-title").textContent = item.title;
    document.querySelector("#case-subtitle").textContent = item.description;
    document.querySelector("#canonical-vertices").textContent =
      JSON.stringify(item.vertices);
    document.querySelector("#all-disk-count").textContent =
      item.counts.allDisk;
    document.querySelector("#annular-count").textContent =
      item.counts.annular;
    document.querySelector("#total-count").textContent = item.counts.total;
    document.querySelector("#polygon-preview").innerHTML =
      polygonSvg(enteredVertices);

    const standardImage = document.querySelector("#standard-image");
    standardImage.src = item.standard.src;
    standardImage.alt = `${item.title}, standard ruling`;
    document.querySelector("#standard-profile").textContent =
      item.standard.profile;
    document.querySelector("#standard-caption").textContent =
      `${item.standard.diskEyes} disk eyes · ` +
      `${item.standard.switches} switches · genus ${item.standard.genus}`;

    document.querySelector("#symmetry-note").textContent = item.symmetry;
    const gallery = document.querySelector("#rational-gallery");
    gallery.replaceChildren(
      ...item.rational.map((drawing) => {
        const figure = document.createElement("figure");
        const image = document.createElement("img");
        const caption = document.createElement("figcaption");
        const badge = document.createElement("span");
        const profile = document.createElement("span");

        image.src = drawing.src;
        image.alt =
          `${item.title}, rational ruling representative ` +
          `${drawing.identifier}, multiplicity ${drawing.multiplicity}`;
        badge.className = "diagram-label";
        badge.textContent =
          `${drawing.identifier} ×${drawing.multiplicity}`;
        profile.textContent =
          `${drawing.profile} · ${drawing.switches} switches`;
        caption.append(badge, profile);
        figure.append(image, caption);
        return figure;
      })
    );

    result.hidden = false;
    result.scrollIntoView({ behavior: "smooth", block: "start" });
  }

  function polygonSvg(vertices) {
    const width = 320;
    const height = 230;
    const padding = 32;
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
    const points = vertices
      .map((point) => project(point).join(","))
      .join(" ");

    const grid = [];
    for (let x = Math.ceil(minX); x <= Math.floor(maxX); x += 1) {
      const [[px]] = [project([x, minY])];
      grid.push(
        `<line x1="${px}" y1="${padding / 2}" x2="${px}" ` +
          `y2="${height - padding / 2}"/>`
      );
    }
    for (let y = Math.ceil(minY); y <= Math.floor(maxY); y += 1) {
      const [, py] = project([minX, y]);
      grid.push(
        `<line x1="${padding / 2}" y1="${py}" ` +
          `x2="${width - padding / 2}" y2="${py}"/>`
      );
    }

    const dots = vertices
      .map((point) => {
        const [x, y] = project(point);
        return `<circle cx="${x}" cy="${y}" r="4.5"/>`;
      })
      .join("");

    return (
      `<svg viewBox="0 0 ${width} ${height}" role="img" ` +
      `aria-label="Newton polygon preview">` +
      `<g stroke="#d9d1c2" stroke-width="1">${grid.join("")}</g>` +
      `<polygon points="${points}" fill="rgba(159,47,37,.13)" ` +
      `stroke="#9f2f25" stroke-width="3"/>` +
      `<g fill="#17202a">${dots}</g>` +
      `</svg>`
    );
  }

  function showMessage(text) {
    message.textContent = text;
    message.hidden = false;
  }

  function hideMessage() {
    message.hidden = true;
    message.textContent = "";
  }

  analyze();
})();
