let expressionData = [];
let annotationData = [];
let svgOriginalText = "";
let currentGene = null;
let currentGwasResults = [];
let currentExpressionTable = [];

const tissues = [
  "Androecium", "Bracteole", "Embryo", "Endosperm",
  "Flower_1", "Flower_2", "Flower_3", "Flower_4", "Flower_5",
  "Flower_bud_1", "Flower_bud_2", "Flower_bud_3", "Flower_bud_4",
  "Germinating_Seedling", "Gynoecium", "Mature_leaf", "Nodule",
  "Pedicel", "Petal", "Pod_Shell", "Root", "Root_Hair", "Root_tip",
  "SAM", "Seed_10_dap", "Seed_20_dap", "Seed_30_dap", "Seed_5_dap",
  "Seed_Coat", "Sepal", "Shoot", "Young_leaf"
];

const tissueGroups = {
  "All tissues": tissues,
  "Root/Nodule": ["Root", "Root_Hair", "Root_tip", "Nodule"],
  "Leaf/Shoot/SAM": ["Mature_leaf", "Young_leaf", "Shoot", "SAM"],
  "Flower organs": [
    "Flower_1", "Flower_2", "Flower_3", "Flower_4", "Flower_5",
    "Flower_bud_1", "Flower_bud_2", "Flower_bud_3", "Flower_bud_4",
    "Androecium", "Gynoecium", "Petal", "Sepal", "Pedicel", "Bracteole"
  ],
  "Seed/Pod": [
    "Seed_5_dap", "Seed_10_dap", "Seed_20_dap", "Seed_30_dap",
    "Seed_Coat", "Embryo", "Endosperm", "Pod_Shell"
  ]
};

const palettes = {
  "Yellow-orange-red": ["#ffffcc", "#ffeda0", "#feb24c", "#f03b20", "#bd0026"],
  "White-yellow-red": ["#f7fbff", "#ffffb2", "#fd8d3c", "#bd0026"],
  "White-orange-red": ["#fff7ec", "#fdd49e", "#fc8d59", "#b30000"],
  "Blue-white-red": ["#2166ac", "#f7f7f7", "#b2182b"],
  "Purple-yellow": ["#2d004b", "#762a83", "#f7f7f7", "#fdb863", "#e66101"],
  "Viridis": ["#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"],
  "Magma": ["#000004", "#3b0f70", "#8c2981", "#de4968", "#fe9f6d", "#fcfdbf"],
  "Plasma": ["#0d0887", "#6a00a8", "#b12a90", "#e16462", "#fca636", "#f0f921"],
  "Cividis": ["#00204c", "#31446b", "#666970", "#958f78", "#c6ba7c", "#ffea46"],
  "Green-yellow-red": ["#006837", "#78c679", "#ffffbf", "#fdae61", "#a50026"],
  "Light-blue-dark-blue": ["#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"]
};

document.addEventListener("DOMContentLoaded", init);

async function init() {
  setupTabs();
  populatePaletteOptions();
  populateGroupOptions();

  expressionData = await loadCsv("data/FPKM_File_RK.csv");
  annotationData = await loadCsv("data/gene_annotation.csv");
  svgOriginalText = await fetch("assets/Chickpea_gene_expression_atlas_RK.svg").then(r => r.text());

  expressionData.forEach(row => {
    tissues.forEach(t => row[t] = Number(row[t]) || 0);
  });

  populateGeneOptions();

  currentGene = expressionData[0].Gene_ID;
  document.getElementById("geneInput").value = currentGene;

  setupEvents();
  renderAtlas();
  renderBarPlot();
}

function loadCsv(path) {
  return new Promise((resolve, reject) => {
    Papa.parse(path, {
      download: true,
      header: true,
      skipEmptyLines: true,
      complete: results => resolve(results.data),
      error: err => reject(err)
    });
  });
}

function setupTabs() {
  document.querySelectorAll(".tab-btn").forEach(btn => {
    btn.addEventListener("click", () => {
      document.querySelectorAll(".tab-btn").forEach(b => b.classList.remove("active"));
      document.querySelectorAll(".tab-content").forEach(tab => tab.classList.remove("active"));

      btn.classList.add("active");
      document.getElementById(btn.dataset.tab).classList.add("active");
    });
  });
}

function populatePaletteOptions() {
  const sel = document.getElementById("paletteSelect");
  Object.keys(palettes).forEach(name => {
    const opt = document.createElement("option");
    opt.value = name;
    opt.textContent = name;
    sel.appendChild(opt);
  });
}

function populateGroupOptions() {
  ["barGroup", "summaryGroup", "heatmapGroup"].forEach(id => {
    const sel = document.getElementById(id);
    Object.keys(tissueGroups).forEach(name => {
      const opt = document.createElement("option");
      opt.value = name;
      opt.textContent = name;
      sel.appendChild(opt);
    });
  });
}

function populateGeneOptions() {
  const geneIds = expressionData.map(row => row.Gene_ID).sort();

  const datalist = document.getElementById("geneList");
  const barGene = document.getElementById("barGene");

  geneIds.forEach(gene => {
    const opt1 = document.createElement("option");
    opt1.value = gene;
    datalist.appendChild(opt1);

    const opt2 = document.createElement("option");
    opt2.value = gene;
    opt2.textContent = gene;
    barGene.appendChild(opt2);
  });
}

function setupEvents() {
  document.getElementById("loadGeneBtn").addEventListener("click", () => {
    const gene = document.getElementById("geneInput").value.trim();
    const found = expressionData.find(row => row.Gene_ID === gene);

    if (!found) {
      document.getElementById("geneStatus").textContent = "Gene not found in expression matrix.";
      document.getElementById("geneStatus").style.color = "#b2182b";
      return;
    }

    currentGene = gene;
    document.getElementById("geneStatus").textContent = "Gene loaded.";
    document.getElementById("geneStatus").style.color = "#1b7837";
    renderAtlas();
  });

  ["maxFpkm", "paletteSelect", "strokeMode", "customStrokeColor", "strokeWidth"].forEach(id => {
    document.getElementById(id).addEventListener("input", renderAtlas);
  });

  document.getElementById("maxFpkm").addEventListener("input", e => {
    document.getElementById("maxFpkmLabel").textContent = e.target.value;
  });

  document.getElementById("strokeWidth").addEventListener("input", e => {
    document.getElementById("strokeWidthLabel").textContent = e.target.value;
  });

  document.getElementById("downloadSvgBtn").addEventListener("click", downloadCurrentSvg);
  document.getElementById("downloadPngBtn").addEventListener("click", () => downloadAtlasImage("png"));
  document.getElementById("downloadJpegBtn").addEventListener("click", () => downloadAtlasImage("jpeg"));

  document.getElementById("barGene").addEventListener("change", renderBarPlot);
  document.getElementById("barGroup").addEventListener("change", renderBarPlot);
  document.getElementById("barLog").addEventListener("change", renderBarPlot);

  document.getElementById("summaryBtn").addEventListener("click", renderSummaryPlot);
  document.getElementById("heatmapBtn").addEventListener("click", renderHeatmap);
  document.getElementById("searchGwasBtn").addEventListener("click", runGwasSearch);
  document.getElementById("downloadGwasCsvBtn").addEventListener("click", () => downloadCsv(currentGwasResults, "GWAS_candidate_genes.csv"));

  document.getElementById("tableBtn").addEventListener("click", renderExpressionTable);
  document.getElementById("downloadExprCsvBtn").addEventListener("click", () => downloadCsv(currentExpressionTable, "selected_gene_expression_table.csv"));
}

function getCurrentGeneRow() {
  return expressionData.find(row => row.Gene_ID === currentGene);
}

function hexToRgb(hex) {
  hex = hex.replace("#", "");
  if (hex.length === 3) {
    hex = hex.split("").map(x => x + x).join("");
  }
  const num = parseInt(hex, 16);
  return {
    r: (num >> 16) & 255,
    g: (num >> 8) & 255,
    b: num & 255
  };
}

function rgbToHex(r, g, b) {
  return "#" + [r, g, b].map(x => {
    const h = Math.round(x).toString(16);
    return h.length === 1 ? "0" + h : h;
  }).join("");
}

function interpolateColor(c1, c2, t) {
  const a = hexToRgb(c1);
  const b = hexToRgb(c2);
  return rgbToHex(
    a.r + (b.r - a.r) * t,
    a.g + (b.g - a.g) * t,
    a.b + (b.b - a.b) * t
  );
}

function expressionToColor(value, maxValue, palette) {
  value = Math.max(0, Math.min(Number(value) || 0, maxValue));
  const scaled = value / maxValue;

  const n = palette.length - 1;
  const pos = scaled * n;
  const i = Math.min(Math.floor(pos), n - 1);
  const t = pos - i;

  return interpolateColor(palette[i], palette[i + 1], t);
}

function getStrokeSettings() {
  const mode = document.getElementById("strokeMode").value;
  const width = document.getElementById("strokeWidth").value;

  if (mode === "No stroke") return { color: "none", width: 0 };
  if (mode === "Black stroke") return { color: "#000000", width };
  if (mode === "Gray stroke") return { color: "#333333", width };
  if (mode === "White stroke") return { color: "#FFFFFF", width };
  return { color: document.getElementById("customStrokeColor").value, width };
}

function generateColoredSvgText() {
  const geneRow = getCurrentGeneRow();
  if (!geneRow) return svgOriginalText;

  const parser = new DOMParser();
  const doc = parser.parseFromString(svgOriginalText, "image/svg+xml");
  const maxFpkm = Number(document.getElementById("maxFpkm").value);
  const palette = palettes[document.getElementById("paletteSelect").value];
  const stroke = getStrokeSettings();

  tissues.forEach(tissue => {
    const node = doc.getElementById(tissue);
    if (node) {
      const color = expressionToColor(geneRow[tissue], maxFpkm, palette);
      node.setAttribute("style", `fill:${color};stroke:${stroke.color};stroke-width:${stroke.width};`);
      node.setAttribute("fill", color);
      node.setAttribute("stroke", stroke.color);
      node.setAttribute("stroke-width", stroke.width);
    }
  });

  return new XMLSerializer().serializeToString(doc);
}

function renderAtlas() {
  const svgText = generateColoredSvgText();
  document.getElementById("svgContainer").innerHTML = svgText;
}

function downloadCurrentSvg() {
  const svgText = generateColoredSvgText();
  downloadText(svgText, `${currentGene}_expression_atlas.svg`, "image/svg+xml");
}

async function downloadAtlasImage(format) {
  const finalWidth = 4500;
  const finalHeight = 3250;
  const svgText = generateColoredSvgText();

  const canvas = document.createElement("canvas");
  canvas.width = finalWidth;
  canvas.height = finalHeight;
  const ctx = canvas.getContext("2d");

  ctx.fillStyle = "white";
  ctx.fillRect(0, 0, finalWidth, finalHeight);

  const atlasImg = await svgToImage(svgText);

  const topMargin = 70;
  const sideMargin = 100;
  const legendHeight = 260;
  const gap = 70;
  const bottomMargin = 70;

  const availableWidth = finalWidth - 2 * sideMargin;
  const availableHeight = finalHeight - topMargin - legendHeight - gap - bottomMargin;

  const scale = Math.min(availableWidth / atlasImg.width, availableHeight / atlasImg.height);
  const drawWidth = atlasImg.width * scale;
  const drawHeight = atlasImg.height * scale;

  const x = (finalWidth - drawWidth) / 2;
  const y = topMargin;

  ctx.drawImage(atlasImg, x, y, drawWidth, drawHeight);

  drawLegend(ctx, finalWidth, finalHeight, legendHeight);

  const mime = format === "jpeg" ? "image/jpeg" : "image/png";
  const ext = format === "jpeg" ? "jpeg" : "png";

  canvas.toBlob(blob => {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = `${currentGene}_expression_atlas.${ext}`;
    a.click();
    URL.revokeObjectURL(a.href);
  }, mime, 0.95);
}

function svgToImage(svgText) {
  return new Promise((resolve, reject) => {
    const blob = new Blob([svgText], { type: "image/svg+xml;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const img = new Image();

    img.onload = () => {
      URL.revokeObjectURL(url);
      resolve(img);
    };

    img.onerror = reject;
    img.src = url;
  });
}

function drawLegend(ctx, finalWidth, finalHeight, legendHeight) {
  const palette = palettes[document.getElementById("paletteSelect").value];
  const maxFpkm = Number(document.getElementById("maxFpkm").value);

  const legendWidth = 1400;
  const barHeight = 95;
  const x0 = (finalWidth - legendWidth) / 2;
  const y0 = finalHeight - legendHeight + 30;

  const grad = ctx.createLinearGradient(x0, y0, x0 + legendWidth, y0);
  palette.forEach((color, i) => {
    grad.addColorStop(i / (palette.length - 1), color);
  });

  ctx.fillStyle = grad;
  ctx.fillRect(x0, y0, legendWidth, barHeight);

  ctx.strokeStyle = "#333333";
  ctx.lineWidth = 2;
  ctx.strokeRect(x0, y0, legendWidth, barHeight);

  ctx.fillStyle = "#111111";
  ctx.font = "34px Arial";
  ctx.textAlign = "center";

  const ticks = [0, maxFpkm * 0.25, maxFpkm * 0.5, maxFpkm * 0.75, maxFpkm];
  ticks.forEach(tick => {
    const x = x0 + (tick / maxFpkm) * legendWidth;
    ctx.beginPath();
    ctx.moveTo(x, y0 + barHeight);
    ctx.lineTo(x, y0 + barHeight + 15);
    ctx.stroke();
    ctx.fillText(Math.round(tick), x, y0 + barHeight + 55);
  });

  ctx.font = "40px Arial";
  ctx.fillText("FPKM expression", finalWidth / 2, y0 + barHeight + 115);
}

function renderBarPlot() {
  const gene = document.getElementById("barGene").value || currentGene;
  const group = document.getElementById("barGroup").value || "All tissues";
  const useLog = document.getElementById("barLog").checked;

  const row = expressionData.find(r => r.Gene_ID === gene);
  if (!row) return;

  const selectedTissues = tissueGroups[group];

  const x = selectedTissues.map(t => useLog ? Math.log2(row[t] + 1) : row[t]);
  const y = selectedTissues;

  Plotly.newPlot("barPlot", [{
    x,
    y,
    type: "bar",
    orientation: "h",
    marker: { color: "#2c7fb8" }
  }], {
    title: `Expression profile of ${gene}`,
    xaxis: { title: useLog ? "log2(FPKM + 1)" : "FPKM" },
    yaxis: { automargin: true },
    margin: { l: 160, r: 30, t: 70, b: 60 }
  }, { responsive: true });
}

function parseGeneText(text) {
  return [...new Set(text.split(/[\s,;]+/).map(x => x.trim()).filter(Boolean))];
}

function renderSummaryPlot() {
  const genes = parseGeneText(document.getElementById("summaryGenes").value);
  const group = document.getElementById("summaryGroup").value || "All tissues";
  const useLog = document.getElementById("summaryLog").checked;
  const selectedTissues = tissueGroups[group];

  const rows = expressionData.filter(r => genes.includes(r.Gene_ID));
  if (rows.length === 0) return alert("No valid genes found.");

  const means = [];
  const medians = [];

  selectedTissues.forEach(tissue => {
    const values = rows.map(r => useLog ? Math.log2(r[tissue] + 1) : r[tissue]);
    means.push(mean(values));
    medians.push(median(values));
  });

  Plotly.newPlot("summaryPlot", [
    { x: means, y: selectedTissues, type: "bar", orientation: "h", name: "Mean" },
    { x: medians, y: selectedTissues, type: "bar", orientation: "h", name: "Median" }
  ], {
    title: "Mean and median expression across selected genes",
    barmode: "group",
    xaxis: { title: useLog ? "log2(FPKM + 1)" : "FPKM" },
    margin: { l: 160, r: 30, t: 70, b: 60 }
  }, { responsive: true });
}

function renderHeatmap() {
  const genes = parseGeneText(document.getElementById("heatmapGenes").value);
  const group = document.getElementById("heatmapGroup").value || "All tissues";
  const scaleType = document.getElementById("heatmapScale").value;
  const selectedTissues = tissueGroups[group];

  const rows = expressionData.filter(r => genes.includes(r.Gene_ID));
  if (rows.length === 0) return alert("No valid genes found.");

  const z = rows.map(row => {
    let values = selectedTissues.map(t => {
      if (scaleType === "Raw FPKM") return row[t];
      return Math.log2(row[t] + 1);
    });

    if (scaleType === "Row-scaled Z-score") {
      const m = mean(values);
      const sd = std(values);
      values = values.map(v => sd === 0 ? 0 : (v - m) / sd);
    }

    return values;
  });

  Plotly.newPlot("heatmapPlot", [{
    z,
    x: selectedTissues,
    y: rows.map(r => r.Gene_ID),
    type: "heatmap",
    colorscale: [
      [0, "#f7fbff"],
      [0.33, "#ffffb2"],
      [0.66, "#fd8d3c"],
      [1, "#bd0026"]
    ]
  }], {
    title: "Gene expression heatmap",
    margin: { l: 160, r: 30, t: 70, b: 150 },
    xaxis: { tickangle: -45 }
  }, { responsive: true });
}

function normalizeChr(x) {
  return String(x).replace(/chr/gi, "").replace(/^0+/, "");
}

function runGwasSearch() {
  const snpId = document.getElementById("snpId").value.trim();
  const chr = normalizeChr(document.getElementById("gwasChr").value);
  const pos = Number(document.getElementById("snpPos").value);
  const upstream = Number(document.getElementById("upstream").value);
  const downstream = Number(document.getElementById("downstream").value);

  const start = Math.max(1, pos - upstream);
  const end = pos + downstream;

  currentGwasResults = annotationData
    .filter(g => normalizeChr(g.Chr) === chr && Number(g.End) >= start && Number(g.Start) <= end)
    .map(g => {
      const midpoint = (Number(g.Start) + Number(g.End)) / 2;
      let direction = "SNP within/overlapping gene";
      if (Number(g.End) < pos) direction = "Upstream";
      if (Number(g.Start) > pos) direction = "Downstream";

      return {
        SNP_ID: snpId,
        Chr: g.Chr,
        SNP_position: pos,
        Region_start: start,
        Region_end: end,
        Gene_ID: g.Gene_ID,
        Gene_raw: g.Gene_raw,
        Start: g.Start,
        End: g.End,
        Strand: g.Strand,
        Gene_length_bp: g.Gene_length_bp,
        Distance_from_SNP_bp: Math.round(Math.abs(midpoint - pos)),
        Direction: direction,
        Annotation: g.Annotation,
        Dbxref: g.Dbxref,
        Present_in_expression_matrix: expressionData.some(r => r.Gene_ID === g.Gene_ID) ? "Yes" : "No"
      };
    })
    .sort((a, b) => a.Distance_from_SNP_bp - b.Distance_from_SNP_bp);

  document.getElementById("gwasSummary").innerHTML =
    `<strong>Search region:</strong> Chr${chr}: ${start.toLocaleString()} - ${end.toLocaleString()}<br>
     <strong>SNP position:</strong> ${pos.toLocaleString()}<br>
     <strong>Candidate genes found:</strong> ${currentGwasResults.length}`;

  renderTable("gwasTable", currentGwasResults);
}

function renderExpressionTable() {
  const genes = parseGeneText(document.getElementById("tableGenes").value);
  currentExpressionTable = expressionData.filter(r => genes.includes(r.Gene_ID));

  if (currentExpressionTable.length === 0) return alert("No valid genes found.");

  renderTable("exprTable", currentExpressionTable);
}

function renderTable(tableId, data) {
  const table = document.getElementById(tableId);
  table.innerHTML = "";

  if (!data || data.length === 0) {
    table.innerHTML = "<tr><td>No data found.</td></tr>";
    return;
  }

  const keys = Object.keys(data[0]);

  const thead = document.createElement("thead");
  const trh = document.createElement("tr");

  keys.forEach(k => {
    const th = document.createElement("th");
    th.textContent = k;
    trh.appendChild(th);
  });

  thead.appendChild(trh);
  table.appendChild(thead);

  const tbody = document.createElement("tbody");

  data.forEach(row => {
    const tr = document.createElement("tr");
    keys.forEach(k => {
      const td = document.createElement("td");
      td.textContent = row[k] ?? "";
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  });

  table.appendChild(tbody);
}

function mean(arr) {
  return arr.reduce((a, b) => a + b, 0) / arr.length;
}

function median(arr) {
  const s = [...arr].sort((a, b) => a - b);
  const mid = Math.floor(s.length / 2);
  return s.length % 2 ? s[mid] : (s[mid - 1] + s[mid]) / 2;
}

function std(arr) {
  const m = mean(arr);
  const variance = mean(arr.map(v => (v - m) ** 2));
  return Math.sqrt(variance);
}

function downloadText(text, filename, mime) {
  const blob = new Blob([text], { type: mime });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  a.click();
  URL.revokeObjectURL(a.href);
}

function downloadCsv(data, filename) {
  if (!data || data.length === 0) return alert("No data to download.");

  const csv = Papa.unparse(data);
  downloadText(csv, filename, "text/csv");
}