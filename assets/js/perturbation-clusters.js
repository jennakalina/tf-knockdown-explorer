/* ============================================================
   Page logic for perturbation-clusters.html
   - Scatter plot of TF knockdown lines (V1/V2 projection), every
     point labeled, colored by cluster.
   - Selecting a cluster (click a point or use the dropdown) loads
     the DAM/DAP data (lazily, once) and renders two heatmaps:
     log2FC by RNAi line, restricted to the union of DAMs/DAPs that
     are significant (adjusted p-value < 0.05) in at least one RNAi
     line belonging to a TF in that cluster.
   ============================================================ */

// Fixed cluster color palette, per the user's spec — do not cycle/reorder.
const CLUSTER_COLORS = {
  1: "#319e77",
  2: "#d95f01",
  3: "#7570b3",
  4: "#e7298a",
  5: "#66a61d",
  6: "#e6ab00",
};

// Control-genotype lines that show up in the cluster assignment but are
// excluded from heatmap construction since they're control backgrounds, not
// TF knockdowns (see README for the "weird cases" found while wiring this
// up — notably Attp40, which DOES have real DAM rows under the literal
// target name "Attp40" and would otherwise show up as an ordinary column).
const CONTROL_EXCLUDE = new Set(["w1118", "Attp2", "Attp40"]);

const DATA_FILES = {
  clusters: "data/clusters.data.js",
  dams: "data/dams.data.js",
  daps: "data/daps.data.js",
};

const state = {
  clusters: null, // { columns, rows }
  dams: null,
  daps: null,
  damsPromise: null,
  dapsPromise: null,
  selectedCluster: null,
  byCluster: null, // Map<cluster, [tf, ...]>
};

const CLUSTER_COL_IDX = { tf: 0, v1: 1, v2: 2, cluster: 3 };
const DA_COL_IDX = { rnai: 0, name: 1, log2fc: 2, pvalue: 3, padj: 4, target: 5, reg: 6 };
const SIG_PADJ_THRESHOLD = 0.05;

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
  }[c]));
}

function fmtNum(n, digits = 3) {
  if (n === null || n === undefined || Number.isNaN(n)) return "NA";
  if (n === 0) return "0";
  if (Math.abs(n) < 0.001) return n.toExponential(2);
  return n.toFixed(digits);
}

// Pick "nice" round tick values (1/2/5 × a power of ten) spanning [lo, hi].
function niceTicksRange(lo, hi, count = 5) {
  const rawStep = (hi - lo) / (count - 1);
  const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const norm = rawStep / mag;
  const niceNorm = norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10;
  const step = niceNorm * mag;
  const start = Math.ceil(lo / step) * step;
  const ticks = [];
  for (let v = start; v <= hi + step * 1e-6; v += step) {
    ticks.push(+v.toFixed(6));
  }
  return ticks;
}

function niceTicks(half, count = 5) {
  return niceTicksRange(-half, half, count);
}

// ---------- diverging color scale for the heatmaps ----------
// Blue for negative log2FC, white for 0, red for positive — the same
// blue/red hex as --status-down / --status-up respectively, which now
// matches the rest of the site's convention after the 2026-08-21 color
// flip (down/discordant = blue, up/concordant = red). This chart still
// encodes raw magnitude/polarity rather than an up/down significance
// call, but the color mapping itself is consistent with everywhere else.
const HEAT_NEG = [42, 120, 214]; // #2a78d6
const HEAT_POS = [227, 73, 72]; // #e34948
const HEAT_MID = [255, 255, 255];

function lerpRgb(a, b, t) {
  return a.map((v, i) => Math.round(v + (b[i] - v) * t));
}

function log2fcColor(value, maxAbs) {
  if (!maxAbs) return "rgb(255,255,255)";
  const t = Math.max(-1, Math.min(1, value / maxAbs));
  const rgb = t < 0 ? lerpRgb(HEAT_MID, HEAT_NEG, -t) : lerpRgb(HEAT_MID, HEAT_POS, t);
  return `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
}

function textColorForBg(rgbString) {
  const m = rgbString.match(/\d+/g).map(Number);
  const lum = (0.299 * m[0] + 0.587 * m[1] + 0.114 * m[2]) / 255;
  return lum > 0.6 ? "#0b0b0b" : "#ffffff";
}

// ---------- data loading (script-tag based, not fetch — see README) ----------

function loadDataFile(key) {
  window.__TF_DATA = window.__TF_DATA || {};
  if (window.__TF_DATA[key]) return Promise.resolve(window.__TF_DATA[key]);
  return new Promise((resolve, reject) => {
    const script = document.createElement("script");
    script.src = DATA_FILES[key];
    script.onload = () => {
      const data = window.__TF_DATA && window.__TF_DATA[key];
      if (!data) { reject(new Error(`Loaded ${DATA_FILES[key]} but no data found for "${key}".`)); return; }
      resolve(data);
    };
    script.onerror = () => reject(new Error(`Failed to load data file: ${DATA_FILES[key]}`));
    document.head.appendChild(script);
  });
}

// ---------- scatter plot ----------

function buildByCluster() {
  const map = new Map();
  state.clusters.rows.forEach((r) => {
    const cluster = r[CLUSTER_COL_IDX.cluster];
    if (!map.has(cluster)) map.set(cluster, []);
    map.get(cluster).push(r[CLUSTER_COL_IDX.tf]);
  });
  state.byCluster = map;
}

function populateClusterDropdown() {
  const select = document.getElementById("cluster-select");
  const clusterNums = [...state.byCluster.keys()].sort((a, b) => a - b);
  select.innerHTML =
    `<option value="">Select a cluster…</option>` +
    clusterNums.map((c) => `<option value="${c}">Cluster ${c}</option>`).join("");
}

function renderLegend() {
  const el = document.getElementById("cluster-legend");
  const clusterNums = [...state.byCluster.keys()].sort((a, b) => a - b);
  el.innerHTML = clusterNums
    .map(
      (c) =>
        `<span class="legend-item"><span class="legend-swatch" style="background:${CLUSTER_COLORS[c]};"></span> Cluster ${c}</span>`
    )
    .join("");
}

function renderScatter() {
  const wrap = document.getElementById("cluster-scatter-wrap");
  const rows = state.clusters.rows;

  const v1s = rows.map((r) => r[CLUSTER_COL_IDX.v1]);
  const v2s = rows.map((r) => r[CLUSTER_COL_IDX.v2]);
  const xHalf = Math.max(Math.abs(Math.min(...v1s)), Math.abs(Math.max(...v1s))) * 1.15 || 1;
  const yHalf = Math.max(Math.abs(Math.min(...v2s)), Math.abs(Math.max(...v2s))) * 1.15 || 1;

  const W = 640, H = 480;
  const marginLeft = 56, marginRight = 24, marginTop = 20, marginBottom = 48;
  const plotW = W - marginLeft - marginRight;
  const plotH = H - marginTop - marginBottom;

  const xScale = (x) => marginLeft + ((x + xHalf) / (2 * xHalf)) * plotW;
  const yScale = (y) => marginTop + (1 - (y + yHalf) / (2 * yHalf)) * plotH;

  const zeroX = xScale(0).toFixed(1);
  const zeroY = yScale(0).toFixed(1);

  const xTicks = niceTicks(xHalf);
  const yTicks = niceTicks(yHalf);

  const xTickSvg = xTicks
    .map((v) => {
      const tx = xScale(v).toFixed(1);
      const bottom = H - marginBottom;
      return `<line class="tick-line" x1="${tx}" y1="${bottom}" x2="${tx}" y2="${bottom + 6}" />
      <text class="tick-label" x="${tx}" y="${bottom + 18}" text-anchor="middle">${fmtNum(v, v === 0 ? 0 : 2)}</text>`;
    })
    .join("");

  const yTickSvg = yTicks
    .map((v) => {
      const ty = yScale(v).toFixed(1);
      return `<line class="tick-line" x1="${marginLeft - 6}" y1="${ty}" x2="${marginLeft}" y2="${ty}" />
      <text class="tick-label" x="${marginLeft - 10}" y="${(+ty + 3.5).toFixed(1)}" text-anchor="end">${fmtNum(v, v === 0 ? 0 : 2)}</text>`;
    })
    .join("");

  const pointsSvg = rows
    .map((r) => {
      const tf = r[CLUSTER_COL_IDX.tf];
      const cluster = r[CLUSTER_COL_IDX.cluster];
      const cx = xScale(r[CLUSTER_COL_IDX.v1]).toFixed(1);
      const cy = yScale(r[CLUSTER_COL_IDX.v2]).toFixed(1);
      const color = CLUSTER_COLORS[cluster] || "#898781";
      return `<g>
        <circle class="cluster-point" data-cluster="${cluster}" cx="${cx}" cy="${cy}" r="5" style="fill:${color};" tabindex="0" aria-label="${escapeHtml(tf)}, cluster ${cluster}"></circle>
        <text class="cluster-point-label" data-cluster="${cluster}" x="${(+cx + 6).toFixed(1)}" y="${(+cy + 3).toFixed(1)}">${escapeHtml(tf)}</text>
      </g>`;
    })
    .join("");

  wrap.innerHTML = `
    <div class="scatter-wrap">
      <svg class="scatter-svg" id="cluster-svg" viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg">
        <rect class="plot-frame" x="${marginLeft}" y="${marginTop}" width="${plotW}" height="${plotH}" />
        <line class="axis-line" x1="${marginLeft}" y1="${zeroY}" x2="${W - marginRight}" y2="${zeroY}" />
        <line class="axis-line" x1="${zeroX}" y1="${marginTop}" x2="${zeroX}" y2="${H - marginBottom}" />
        ${xTickSvg}
        ${yTickSvg}
        ${pointsSvg}
        <text class="axis-title" x="${W / 2}" y="${H - 10}" text-anchor="middle">Dimension 1</text>
        <text class="axis-title" x="${14}" y="${H / 2}" text-anchor="middle" transform="rotate(-90 14 ${H / 2})">Dimension 2</text>
      </svg>
    </div>
  `;

  document.getElementById("cluster-svg").addEventListener("click", (e) => {
    const circle = e.target.closest("circle[data-cluster]");
    if (!circle) return;
    selectCluster(+circle.dataset.cluster);
  });

  document.getElementById("cluster-svg").addEventListener("keydown", (e) => {
    if (e.key !== "Enter" && e.key !== " ") return;
    const circle = e.target.closest("circle[data-cluster]");
    if (!circle) return;
    e.preventDefault();
    selectCluster(+circle.dataset.cluster);
  });
}

function updateScatterHighlight() {
  const selected = state.selectedCluster;
  document.querySelectorAll("#cluster-svg .cluster-point").forEach((el) => {
    const isSelected = selected !== null && +el.dataset.cluster === selected;
    const isOther = selected !== null && !isSelected;
    el.classList.toggle("is-selected", isSelected);
    el.classList.toggle("is-dimmed", isOther);
  });
  document.querySelectorAll("#cluster-svg .cluster-point-label").forEach((el) => {
    const isOther = selected !== null && +el.dataset.cluster !== selected;
    el.classList.toggle("is-dimmed", isOther);
  });
}

// ---------- heatmaps ----------

function computeHeatmap(dataset, clusterTFs) {
  const idx = DA_COL_IDX;
  const tfSet = new Set(clusterTFs);
  const filtered = dataset.rows.filter((r) => tfSet.has(r[idx.target]));
  if (!filtered.length) return { rowNames: [], cols: [], valueMap: new Map(), maxAbs: 0 };

  const sigNames = new Set();
  filtered.forEach((r) => {
    if (r[idx.padj] < SIG_PADJ_THRESHOLD) sigNames.add(r[idx.name]);
  });

  const colMap = new Map(); // rnai -> target
  filtered.forEach((r) => colMap.set(r[idx.rnai], r[idx.target]));
  const cols = [...colMap.entries()]
    .map(([rnai, target]) => ({ rnai, target }))
    .sort((a, b) => a.target.localeCompare(b.target) || a.rnai.localeCompare(b.rnai));

  const rowNames = [...sigNames].sort((a, b) => a.localeCompare(b));

  const valueMap = new Map();
  let maxAbs = 0;
  filtered.forEach((r) => {
    if (sigNames.has(r[idx.name])) {
      const v = r[idx.log2fc];
      valueMap.set(r[idx.name] + " " + r[idx.rnai], v);
      if (Math.abs(v) > maxAbs) maxAbs = Math.abs(v);
    }
  });

  return { rowNames, cols, valueMap, maxAbs: maxAbs || 1 };
}

function renderHeatmap(containerId, heatmapData, datasetLabel) {
  const wrap = document.getElementById(containerId);
  if (!heatmapData || !heatmapData.rowNames.length) {
    wrap.innerHTML = `<div class="empty-state">No significant ${datasetLabel} (adjusted p &lt; 0.05) found for this cluster.</div>`;
    return;
  }
  const { rowNames, cols, valueMap, maxAbs } = heatmapData;

  const groups = [];
  cols.forEach((c) => {
    const last = groups[groups.length - 1];
    if (last && last.target === c.target) last.rnais.push(c.rnai);
    else groups.push({ target: c.target, rnais: [c.rnai] });
  });

  const headerRow1 = groups
    .map((g) => `<th class="heatmap-col-tf" colspan="${g.rnais.length}">${escapeHtml(g.target)}</th>`)
    .join("");
  const headerRow2 = cols.map((c) => `<th class="heatmap-col-rnai">${escapeHtml(c.rnai)}</th>`).join("");

  const bodyRows = rowNames
    .map((name) => {
      const cells = cols
        .map((c) => {
          const v = valueMap.get(name + " " + c.rnai);
          if (v === undefined) return `<td class="heatmap-cell"></td>`;
          const bg = log2fcColor(v, maxAbs);
          const fg = textColorForBg(bg);
          const title = `${escapeHtml(name)} × ${escapeHtml(c.rnai)} (${escapeHtml(c.target)}): log2FC ${fmtNum(v, 2)}`;
          return `<td class="heatmap-cell" style="background:${bg};color:${fg};" title="${title}">${fmtNum(v, 2)}</td>`;
        })
        .join("");
      return `<tr><td class="heatmap-row-label" title="${escapeHtml(name)}">${escapeHtml(name)}</td>${cells}</tr>`;
    })
    .join("");

  const negSwatch = log2fcColor(-maxAbs, maxAbs);
  const posSwatch = log2fcColor(maxAbs, maxAbs);

  wrap.innerHTML = `
    <div class="heatmap-meta">${rowNames.length.toLocaleString()} ${datasetLabel} significant (adjusted p &lt; 0.05) in at least one RNAi line in this cluster &times; ${cols.length} RNAi line${cols.length === 1 ? "" : "s"} across ${groups.length} TF${groups.length === 1 ? "" : "s"}</div>
    <div class="heatmap-scroll">
      <table class="heatmap-table">
        <thead>
          <tr><th class="heatmap-corner"></th>${headerRow1}</tr>
          <tr><th class="heatmap-corner heatmap-corner-2"></th>${headerRow2}</tr>
        </thead>
        <tbody>${bodyRows}</tbody>
      </table>
    </div>
    <div class="heatmap-legend">
      <span>${fmtNum(-maxAbs, 2)}</span>
      <span class="heatmap-legend-bar" style="background:linear-gradient(to right, ${negSwatch} 0%, #ffffff 50%, ${posSwatch} 100%);"></span>
      <span>${fmtNum(maxAbs, 2)}</span>
      <span style="margin-left: 6px;">log2FC</span>
    </div>
  `;
}

function renderHeatmapsForCluster(cluster) {
  const clusterTFs = (state.byCluster.get(cluster) || []).filter((tf) => !CONTROL_EXCLUDE.has(tf));

  document.getElementById("dam-heatmap-panel").style.display = "";
  document.getElementById("dap-heatmap-panel").style.display = "";
  document.getElementById("heatmap-hint").style.display = "none";

  const damHeatmap = computeHeatmap(state.dams, clusterTFs);
  const dapHeatmap = computeHeatmap(state.daps, clusterTFs);
  renderHeatmap("dam-heatmap-wrap", damHeatmap, "DAMs");
  renderHeatmap("dap-heatmap-wrap", dapHeatmap, "DAPs");
}

function selectCluster(cluster) {
  state.selectedCluster = cluster;
  document.getElementById("cluster-select").value = String(cluster);
  updateScatterHighlight();

  document.getElementById("dam-heatmap-wrap").innerHTML = `<div class="loading-state">Loading differential analysis data…</div>`;
  document.getElementById("dap-heatmap-wrap").innerHTML = `<div class="loading-state">Loading differential analysis data…</div>`;
  document.getElementById("dam-heatmap-panel").style.display = "";
  document.getElementById("dap-heatmap-panel").style.display = "";
  document.getElementById("heatmap-hint").style.display = "none";

  if (!state.damsPromise) {
    state.damsPromise = loadDataFile("dams").then((data) => { state.dams = data; return data; });
  }
  if (!state.dapsPromise) {
    state.dapsPromise = loadDataFile("daps").then((data) => { state.daps = data; return data; });
  }

  Promise.all([state.damsPromise, state.dapsPromise])
    .then(() => {
      // A newer selection may have started while these were in flight.
      if (state.selectedCluster === cluster) renderHeatmapsForCluster(cluster);
    })
    .catch((err) => {
      document.getElementById("dam-heatmap-wrap").innerHTML = `<div class="empty-state">Couldn't load differential analysis data: ${escapeHtml(err.message)}</div>`;
      document.getElementById("dap-heatmap-wrap").innerHTML = "";
      console.error(err);
    });
}

function clearSelection() {
  state.selectedCluster = null;
  updateScatterHighlight();
  document.getElementById("dam-heatmap-panel").style.display = "none";
  document.getElementById("dap-heatmap-panel").style.display = "none";
  document.getElementById("heatmap-hint").style.display = "";
}

// ---------- init ----------

async function init() {
  renderNav("perturbation-clusters");

  document.getElementById("cluster-select").addEventListener("change", (e) => {
    const val = e.target.value;
    if (val === "") clearSelection();
    else selectCluster(+val);
  });

  try {
    const data = await loadDataFile("clusters");
    state.clusters = data;
    buildByCluster();
    populateClusterDropdown();
    renderLegend();
    renderScatter();
  } catch (err) {
    document.getElementById("cluster-scatter-wrap").innerHTML =
      `<div class="empty-state">Couldn't load cluster data: ${escapeHtml(err.message)}</div>`;
    console.error(err);
  }
}

document.addEventListener("DOMContentLoaded", init);
