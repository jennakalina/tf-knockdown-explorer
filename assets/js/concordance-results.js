/* ============================================================
  Page logic for concordance-results.html
Two modes:
  - Single analyte search: filter conc_results by a protein or
metabolite substring, sorted by padj.
- PMI pair search: pick one exact protein + one exact metabolite,
show the matching row's stats plus a quadrant scatter plot
     built from the per-sample scaled metabolite/protein values.
   ============================================================ */

const DIR_LABEL = { c: "Concordant", d: "Discordant", ns: "Not sig." };
const DATA_FILES = {
  conc: "data/conc.data.js",
  scaled: "data/scaled.data.js",
};

const state = {
  mode: "single",
  conc: null, // { columns, rows }
  scaled: null, // { samples, targets, metCols, metValues, protCols, protValues }
  pairIndex: null, // Map "proteinmetabolite" -> row
  proteinSet: null,
  metaboliteSet: null,
  analyteQuery: "",
  analyteFiltered: [],
  analyteTable: null,
  showLabels: true,
};

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

// ---------- data loading (script-tag based, not fetch — see README for why) ----------
  
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

// ---------- mode toggle ----------
  
  function setMode(mode) {
    state.mode = mode;
    document.querySelectorAll("#mode-toggle button").forEach((btn) => {
      btn.classList.toggle("active", btn.dataset.mode === mode);
    });
    document.getElementById("single-panel").style.display = mode === "single" ? "" : "none";
    document.getElementById("pair-panel").style.display = mode === "pair" ? "" : "none";
    
    if (mode === "pair" && !state.scaled) {
      document.getElementById("pair-hint").textContent = "Loading sample data…";
      loadDataFile("scaled").then((data) => {
        state.scaled = data;
        document.getElementById("pair-hint").textContent = "Select an exact protein and metabolite to see their result.";
        maybeRenderPair();
      }).catch((err) => {
        document.getElementById("pair-hint").textContent = "Couldn't load sample data: " + err.message;
      });
    }
  }

// ---------- single analyte search ----------
  
  function rowToObj(r) {
    return { protein: r[0], metabolite: r[1], concordance: r[2], variance: r[3], pvalue: r[4], padj: r[5], direction: r[6] };
  }

function applyAnalyteFilter() {
  const q = state.analyteQuery.trim().toLowerCase();
  if (!q) { state.analyteFiltered = []; return; }
  state.analyteFiltered = state.conc.rows.filter(
    (r) => r[0].toLowerCase().includes(q) || r[1].toLowerCase().includes(q)
  );
}

function analyteColumns() {
  return [
    { key: "protein", label: "Protein", flex: 1.2, render: (r) => `<strong>${escapeHtml(r.protein)}</strong>` },
    { key: "metabolite", label: "Metabolite", flex: 1.4, render: (r) => escapeHtml(r.metabolite) },
    { key: "concordance", label: "Concordance", flex: 0.9, render: (r) => fmtNum(r.concordance) },
    { key: "variance", label: "Variance", flex: 0.9, render: (r) => fmtNum(r.variance, 4) },
    { key: "pvalue", label: "p-value", flex: 0.9, render: (r) => fmtNum(r.pvalue, 4) },
    { key: "padj", label: "p.adj", flex: 0.9, render: (r) => fmtNum(r.padj, 4) },
    { key: "direction", label: "Direction", flex: 0.9, render: (r) => `<span class="reg-pill reg-${r.direction}">${DIR_LABEL[r.direction] || r.direction}</span>` },
  ];
}

function renderAnalyteResults() {
  const container = document.getElementById("analyte-results-table");
  const countEl = document.getElementById("analyte-result-count");
  
  if (!state.analyteQuery.trim()) {
    container.innerHTML = `<div class="empty-state">Type a protein or metabolite name above to search.</div>`;
    countEl.textContent = "";
    state.analyteTable = null;
    return;
  }
  
  if (state.analyteFiltered.length === 0) {
    container.innerHTML = `<div class="empty-state">No rows match "${escapeHtml(state.analyteQuery)}".</div>`;
    countEl.textContent = "0 rows";
    state.analyteTable = null;
    return;
  }
  
  countEl.textContent = `${state.analyteFiltered.length.toLocaleString()} rows`;
  
  if (!state.analyteTable) {
    state.analyteTable = new VirtualTable({
      container,
      rowHeight: 40,
      columns: analyteColumns(),
      getRowCount: () => state.analyteFiltered.length,
      getRow: (i) => rowToObj(state.analyteFiltered[i]),
      rowClass: (row) => (row.direction === "c" ? "reg-c" : row.direction === "d" ? "reg-d" : ""),
    });
  }
  state.analyteTable.refresh();
}

// ---------- PMI pair search ----------
  
  function buildPairIndex() {
    const map = new Map();
    const proteinSet = new Set();
    const metaboliteSet = new Set();
    for (const r of state.conc.rows) {
      map.set(r[0] + "" + r[1], r);
      proteinSet.add(r[0]);
      metaboliteSet.add(r[1]);
    }
    state.pairIndex = map;
    state.proteinSet = proteinSet;
    state.metaboliteSet = metaboliteSet;
  }

function populateDatalists() {
  const proteinList = document.getElementById("protein-datalist");
  const metaboliteList = document.getElementById("metabolite-datalist");
  const proteins = [...state.proteinSet].sort();
  const metabolites = [...state.metaboliteSet].sort();
  proteinList.innerHTML = proteins.map((p) => `<option value="${escapeHtml(p)}"></option>`).join("");
  metaboliteList.innerHTML = metabolites.map((m) => `<option value="${escapeHtml(m)}"></option>`).join("");
}

function maybeRenderPair() {
  const proteinInput = document.getElementById("protein-input");
  const metaboliteInput = document.getElementById("metabolite-input");
  const hint = document.getElementById("pair-hint");
  const protein = proteinInput.value.trim();
  const metabolite = metaboliteInput.value.trim();
  
  const proteinValid = state.proteinSet && state.proteinSet.has(protein);
  const metaboliteValid = state.metaboliteSet && state.metaboliteSet.has(metabolite);
  
  proteinInput.classList.toggle("invalid", protein.length > 0 && !proteinValid);
  metaboliteInput.classList.toggle("invalid", metabolite.length > 0 && !metaboliteValid);
  
  if (!proteinValid || !metaboliteValid) {
    hint.textContent = "Select an exact protein and metabolite from the list.";
    hint.classList.remove("is-invalid");
    document.getElementById("pair-output").innerHTML =
      `<div class="empty-state">Choose a protein and a metabolite above.</div>`;
    return;
  }
  
  const row = state.pairIndex.get(protein + "" + metabolite);
  if (!row) {
    hint.textContent = `No concordance result found for ${protein} × ${metabolite}.`;
    hint.classList.add("is-invalid");
    document.getElementById("pair-output").innerHTML = `<div class="empty-state">No data for this pair.</div>`;
    return;
  }
  
  hint.textContent = "";
  hint.classList.remove("is-invalid");
  renderPairOutput(protein, metabolite, rowToObj(row));
}

function renderPairOutput(protein, metabolite, result) {
  const container = document.getElementById("pair-output");
  
  const statTiles = `
  <div class="stat-tile-row">
    <div class="stat-tile"><span class="stat-tile-label">Concordance</span><span class="stat-tile-value">${fmtNum(result.concordance)}</span></div>
      <div class="stat-tile"><span class="stat-tile-label">Variance</span><span class="stat-tile-value">${fmtNum(result.variance, 4)}</span></div>
        <div class="stat-tile"><span class="stat-tile-label">p-value</span><span class="stat-tile-value">${fmtNum(result.pvalue, 4)}</span></div>
          <div class="stat-tile"><span class="stat-tile-label">p.adj</span><span class="stat-tile-value">${fmtNum(result.padj, 4)}</span></div>
            <div class="stat-tile"><span class="stat-tile-label">Direction</span><span class="stat-tile-value"><span class="reg-pill reg-${result.direction}">${DIR_LABEL[result.direction] || result.direction}</span></span></div>
              </div>
              `;
            
            if (!state.scaled) {
              container.innerHTML = statTiles + `<div class="empty-state">Loading sample data for the quadrant plot…</div>`;
              return;
            }
            
            container.innerHTML = statTiles + `
            <div class="scatter-panel">
              <div class="scatter-panel-header">
              <span class="scatter-title">Concordance quadrant plot: ${escapeHtml(protein)} vs ${escapeHtml(metabolite)}</span>
              <label class="checkbox-field">
              <input type="checkbox" id="labels-toggle" ${state.showLabels ? "checked" : ""} />
              Show point labels
            </label>
              </div>
              <div class="scatter-wrap" id="scatter-wrap"></div>
              <div class="legend">
              <span class="legend-item"><span class="legend-swatch swatch-up"></span> Concordant</span>
              <span class="legend-item"><span class="legend-swatch swatch-down"></span> Discordant</span>
              </div>
              </div>
              `;
            
            renderScatter(protein, metabolite, result);
            
            document.getElementById("labels-toggle").addEventListener("change", (e) => {
              state.showLabels = e.target.checked;
              document.querySelectorAll(".scatter-svg .point-label").forEach((el) => {
                el.style.display = state.showLabels ? "" : "none";
              });
            });
}

// Pick "nice" round tick values (1/2/5 × a power of ten) spanning
// [-half, half], symmetric around 0 so 0 itself is always a tick.
function niceTicks(half, count = 5) {
  const rawStep = (half * 2) / (count - 1);
  const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const norm = rawStep / mag;
  const niceNorm = norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10;
  const step = niceNorm * mag;
  const maxK = Math.ceil(half / step);
  const ticks = [];
  for (let k = -maxK; k <= maxK; k++) {
    const v = +(k * step).toFixed(6);
    if (v >= -half - step * 1e-6 && v <= half + step * 1e-6) ticks.push(v);
  }
  return ticks;
}

function renderScatter(protein, metabolite, result) {
  const wrap = document.getElementById("scatter-wrap");
  const metIdx = state.scaled.metCols.indexOf(metabolite);
  const protIdx = state.scaled.protCols.indexOf(protein);
  
  if (metIdx === -1 || protIdx === -1) {
    wrap.innerHTML = `<div class="empty-state">Couldn't find per-sample values for this pair.</div>`;
    return;
  }
  
  const points = state.scaled.samples.map((sample, i) => {
    const x = state.scaled.metValues[i][metIdx];
    const y = state.scaled.protValues[i][protIdx];
    const target = state.scaled.targets[i];
    const quadrant = (x >= 0 && y >= 0) || (x <= 0 && y <= 0) ? "c" : "d";
    return { sample, target, x, y, quadrant };
  });
  
  const xs = points.map((p) => p.x);
  const ys = points.map((p) => p.y);
  const xHalf = Math.max(Math.abs(Math.min(...xs)), Math.abs(Math.max(...xs))) * 1.15 || 1;
  const yHalf = Math.max(Math.abs(Math.min(...ys)), Math.abs(Math.max(...ys))) * 1.15 || 1;
  
  // Layout
  const W = 640, H = 440;
  const marginLeft = 64, marginRight = 24, marginTop = 20, marginBottom = 56;
  const plotW = W - marginLeft - marginRight;
  const plotH = H - marginTop - marginBottom;
  
  const xScale = (x) => marginLeft + ((x + xHalf) / (2 * xHalf)) * plotW;
  const yScale = (y) => marginTop + (1 - (y + yHalf) / (2 * yHalf)) * plotH; // flip: up = smaller svg-y
  
  const circles = points.map((p) => {
    const cx = xScale(p.x).toFixed(1);
    const cy = yScale(p.y).toFixed(1);
    const cls = p.quadrant === "c" ? "pt-c" : "pt-d";
    const labelSvg = `<text class="point-label" x="${(+cx + 6).toFixed(1)}" y="${(+cy + 3).toFixed(1)}" style="${state.showLabels ? "" : "display:none;"}">${escapeHtml(p.target)}</text>`;
    return `<g>
      <circle class="${cls}" cx="${cx}" cy="${cy}" r="5"><title>${escapeHtml(p.target)} (${escapeHtml(p.sample)})\nMetabolite log2FC: ${fmtNum(p.x)}\nProtein log2FC: ${fmtNum(p.y)}</title></circle>
        ${labelSvg}
      </g>`;
  }).join("");
  
  const zeroX = xScale(0).toFixed(1);
  const zeroY = yScale(0).toFixed(1);
  
  const xTicks = niceTicks(xHalf);
  const yTicks = niceTicks(yHalf);
  
  const xTickSvg = xTicks.map((v) => {
    const tx = xScale(v).toFixed(1);
    const bottom = H - marginBottom;
    return `<line class="tick-line" x1="${tx}" y1="${bottom}" x2="${tx}" y2="${bottom + 6}" />
      <text class="tick-label" x="${tx}" y="${bottom + 18}" text-anchor="middle">${fmtNum(v, v === 0 ? 0 : 2)}</text>`;
  }).join("");
  
  const yTickSvg = yTicks.map((v) => {
    const ty = yScale(v).toFixed(1);
    return `<line class="tick-line" x1="${marginLeft - 6}" y1="${ty}" x2="${marginLeft}" y2="${ty}" />
      <text class="tick-label" x="${marginLeft - 10}" y="${(+ty + 3.5).toFixed(1)}" text-anchor="end">${fmtNum(v, v === 0 ? 0 : 2)}</text>`;
  }).join("");
  
  const svg = `
  <svg class="scatter-svg" viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg">
    <rect class="plot-frame" x="${marginLeft}" y="${marginTop}" width="${plotW}" height="${plotH}" />
    
    <line class="axis-line" x1="${marginLeft}" y1="${zeroY}" x2="${W - marginRight}" y2="${zeroY}" />
    <line class="axis-line" x1="${zeroX}" y1="${marginTop}" x2="${zeroX}" y2="${H - marginBottom}" />
    
    ${xTickSvg}
  ${yTickSvg}
  
  ${circles}
  
  <text class="annotation-text" x="${marginLeft + 10}" y="${marginTop + 16}">Concordance = ${fmtNum(result.concordance)}</text>
    <text class="annotation-text" x="${marginLeft + 10}" y="${marginTop + 32}">Adjusted p-value = ${fmtNum(result.padj, 4)}</text>
    
    <text class="axis-title" x="${W / 2}" y="${H - 14}" text-anchor="middle">Metabolite log2FC: ${escapeHtml(metabolite)}</text>
    <text class="axis-title" x="${16}" y="${H / 2}" text-anchor="middle" transform="rotate(-90 16 ${H / 2})">Protein log2FC: ${escapeHtml(protein)}</text>
    </svg>
    `;
  
  wrap.innerHTML = svg;
}

// ---------- init ----------
  
  async function init() {
    renderNav("concordance-results");
    
    document.querySelectorAll("#mode-toggle button").forEach((btn) => {
      btn.addEventListener("click", () => setMode(btn.dataset.mode));
    });
    
    document.getElementById("analyte-search-input").addEventListener("input", (e) => {
      state.analyteQuery = e.target.value;
      applyAnalyteFilter();
      renderAnalyteResults();
    });
    
    document.getElementById("protein-input").addEventListener("input", maybeRenderPair);
    document.getElementById("metabolite-input").addEventListener("input", maybeRenderPair);
    
    document.getElementById("analyte-results-table").innerHTML = `<div class="loading-state">Loading concordance results…</div>`;
    try {
      state.conc = await loadDataFile("conc");
      buildPairIndex();
      populateDatalists();
      document.getElementById("analyte-results-table").innerHTML =
        `<div class="empty-state">Type a protein or metabolite name above to search.</div>`;
    } catch (err) {
      document.getElementById("analyte-results-table").innerHTML =
        `<div class="empty-state">Couldn't load the data file: ${escapeHtml(err.message)}</div>`;
      console.error(err);
    }
  }

document.addEventListener("DOMContentLoaded", init);    const an = a === null || a === undefined || Number.isNaN(a) ? Infinity : a;
    const bn = b === null || b === undefined || Number.isNaN(b) ? Infinity : b;
    return an - bn;
  }
  return String(a).localeCompare(String(b));
}
 
function sortRows(rows, idxMap, typeMap, key, dir) {
  const idx = idxMap[key];
  const type = typeMap[key] || "string";
  const factor = dir === "desc" ? -1 : 1;
  return rows.slice().sort((a, b) => factor * compareValues(a[idx], b[idx], type));
}
 
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
 
// ---------- data loading (script-tag based, not fetch — see README for why) ----------
 
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
 
// ---------- mode toggle ----------
 
function setMode(mode) {
  state.mode = mode;
  document.querySelectorAll("#mode-toggle button").forEach((btn) => {
    btn.classList.toggle("active", btn.dataset.mode === mode);
  });
  document.getElementById("single-panel").style.display = mode === "single" ? "" : "none";
  document.getElementById("pair-panel").style.display = mode === "pair" ? "" : "none";
 
  if (mode === "pair" && !state.scaled) {
    document.getElementById("pair-hint").textContent = "Loading sample data…";
    loadDataFile("scaled").then((data) => {
      state.scaled = data;
      document.getElementById("pair-hint").textContent = "Select an exact protein and metabolite to see their result.";
      maybeRenderPair();
    }).catch((err) => {
      document.getElementById("pair-hint").textContent = "Couldn't load sample data: " + err.message;
    });
  }
}
 
// ---------- single analyte search ----------
 
function rowToObj(r) {
  return { protein: r[0], metabolite: r[1], concordance: r[2], variance: r[3], pvalue: r[4], padj: r[5], direction: r[6] };
}
 
function applyAnalyteFilter() {
  const q = state.analyteQuery.trim().toLowerCase();
  if (!q) { state.analyteFiltered = []; return; }
  const filtered = state.conc.rows.filter(
    (r) => r[0].toLowerCase().includes(q) || r[1].toLowerCase().includes(q)
  );
  state.analyteFiltered = sortRows(filtered, CONC_COL_IDX, CONC_COL_TYPE, state.analyteSort.key, state.analyteSort.dir);
}
 
function refreshAnalyte() {
  applyAnalyteFilter();
  renderAnalyteResults();
}
 
function handleAnalyteSort(key) {
  if (state.analyteSort.key === key) {
    state.analyteSort = { key, dir: state.analyteSort.dir === "asc" ? "desc" : "asc" };
  } else {
    state.analyteSort = { key, dir: "asc" };
  }
  if (state.analyteTable) state.analyteTable.setSort(state.analyteSort.key, state.analyteSort.dir);
  refreshAnalyte();
}
 
function analyteColumns() {
  return [
    { key: "protein", label: "Protein", flex: 1.2, render: (r) => `<strong>${escapeHtml(r.protein)}</strong>` },
    { key: "metabolite", label: "Metabolite", flex: 1.4, render: (r) => escapeHtml(r.metabolite) },
    { key: "concordance", label: "Concordance", flex: 0.9, render: (r) => fmtNum(r.concordance) },
    { key: "variance", label: "Variance", flex: 0.9, render: (r) => fmtNum(r.variance, 4) },
    { key: "pvalue", label: "p-value", flex: 0.9, render: (r) => fmtNum(r.pvalue, 4) },
    { key: "padj", label: "p.adj", flex: 0.9, render: (r) => fmtNum(r.padj, 4) },
    { key: "direction", label: "Direction", flex: 0.9, render: (r) => `<span class="reg-pill reg-${r.direction}">${DIR_LABEL[r.direction] || r.direction}</span>` },
  ];
}
 
function renderAnalyteResults() {
  const container = document.getElementById("analyte-results-table");
  const countEl = document.getElementById("analyte-result-count");
 
  if (!state.analyteQuery.trim()) {
    container.innerHTML = `<div class="empty-state">Type a protein or metabolite name above to search.</div>`;
    countEl.textContent = "";
    state.analyteTable = null;
    return;
  }
 
  if (state.analyteFiltered.length === 0) {
    container.innerHTML = `<div class="empty-state">No rows match "${escapeHtml(state.analyteQuery)}".</div>`;
    countEl.textContent = "0 rows";
    state.analyteTable = null;
    return;
  }
 
  countEl.textContent = `${state.analyteFiltered.length.toLocaleString()} rows`;
 
  if (!state.analyteTable) {
    state.analyteTable = new VirtualTable({
      container,
      rowHeight: 40,
      columns: analyteColumns(),
      getRowCount: () => state.analyteFiltered.length,
      getRow: (i) => rowToObj(state.analyteFiltered[i]),
      rowClass: (row) => (row.direction === "c" ? "reg-c" : row.direction === "d" ? "reg-d" : ""),
      onSort: handleAnalyteSort,
      sortState: state.analyteSort,
    });
  }
  state.analyteTable.refresh();
}
 
// ---------- PMI pair search ----------
 
function buildPairIndex() {
  const map = new Map();
  const proteinSet = new Set();
  const metaboliteSet = new Set();
  for (const r of state.conc.rows) {
    map.set(r[0] + "" + r[1], r);
    proteinSet.add(r[0]);
    metaboliteSet.add(r[1]);
  }
  state.pairIndex = map;
  state.proteinSet = proteinSet;
  state.metaboliteSet = metaboliteSet;
}
 
function populateDatalists() {
  const proteinList = document.getElementById("protein-datalist");
  const metaboliteList = document.getElementById("metabolite-datalist");
  const proteins = [...state.proteinSet].sort();
  const metabolites = [...state.metaboliteSet].sort();
  proteinList.innerHTML = proteins.map((p) => `<option value="${escapeHtml(p)}"></option>`).join("");
  metaboliteList.innerHTML = metabolites.map((m) => `<option value="${escapeHtml(m)}"></option>`).join("");
}
 
function maybeRenderPair() {
  const proteinInput = document.getElementById("protein-input");
  const metaboliteInput = document.getElementById("metabolite-input");
  const hint = document.getElementById("pair-hint");
  const protein = proteinInput.value.trim();
  const metabolite = metaboliteInput.value.trim();
 
  const proteinValid = state.proteinSet && state.proteinSet.has(protein);
  const metaboliteValid = state.metaboliteSet && state.metaboliteSet.has(metabolite);
 
  proteinInput.classList.toggle("invalid", protein.length > 0 && !proteinValid);
  metaboliteInput.classList.toggle("invalid", metabolite.length > 0 && !metaboliteValid);
 
  if (!proteinValid || !metaboliteValid) {
    hint.textContent = "Select an exact protein and metabolite from the list.";
    hint.classList.remove("is-invalid");
    document.getElementById("pair-output").innerHTML =
      `<div class="empty-state">Choose a protein and a metabolite above.</div>`;
    return;
  }
 
  const row = state.pairIndex.get(protein + "" + metabolite);
  if (!row) {
    hint.textContent = `No concordance result found for ${protein} × ${metabolite}.`;
    hint.classList.add("is-invalid");
    document.getElementById("pair-output").innerHTML = `<div class="empty-state">No data for this pair.</div>`;
    return;
  }
 
  hint.textContent = "";
  hint.classList.remove("is-invalid");
  renderPairOutput(protein, metabolite, rowToObj(row));
}
 
function renderPairOutput(protein, metabolite, result) {
  const container = document.getElementById("pair-output");
 
  const statTiles = `
    <div class="stat-tile-row">
      <div class="stat-tile"><span class="stat-tile-label">Concordance</span><span class="stat-tile-value">${fmtNum(result.concordance)}</span></div>
      <div class="stat-tile"><span class="stat-tile-label">Variance</span><span class="stat-tile-value">${fmtNum(result.variance, 4)}</span></div>
      <div class="stat-tile"><span class="stat-tile-label">p-value</span><span class="stat-tile-value">${fmtNum(result.pvalue, 4)}</span></div>
      <div class="stat-tile"><span class="stat-tile-label">p.adj</span><span class="stat-tile-value">${fmtNum(result.padj, 4)}</span></div>
      <div class="stat-tile"><span class="stat-tile-label">Direction</span><span class="stat-tile-value"><span class="reg-pill reg-${result.direction}">${DIR_LABEL[result.direction] || result.direction}</span></span></div>
    </div>
  `;
 
  if (!state.scaled) {
    container.innerHTML = statTiles + `<div class="empty-state">Loading sample data for the quadrant plot…</div>`;
    return;
  }
 
  container.innerHTML = statTiles + `
    <div class="scatter-panel">
      <div class="scatter-panel-header">
        <span class="scatter-title">Concordance quadrant plot: ${escapeHtml(protein)} vs ${escapeHtml(metabolite)}</span>
        <label class="checkbox-field">
          <input type="checkbox" id="labels-toggle" ${state.showLabels ? "checked" : ""} />
          Show point labels
        </label>
      </div>
      <div class="scatter-wrap" id="scatter-wrap"></div>
      <div class="legend">
        <span class="legend-item"><span class="legend-swatch swatch-up"></span> Concordant</span>
        <span class="legend-item"><span class="legend-swatch swatch-down"></span> Discordant</span>
      </div>
    </div>
  `;
 
  renderScatter(protein, metabolite, result);
 
  document.getElementById("labels-toggle").addEventListener("change", (e) => {
    state.showLabels = e.target.checked;
    document.querySelectorAll(".scatter-svg .point-label").forEach((el) => {
      el.style.display = state.showLabels ? "" : "none";
    });
  });
}
 
// Pick "nice" round tick values (1/2/5 × a power of ten) spanning
// [-half, half], symmetric around 0 so 0 itself is always a tick.
function niceTicks(half, count = 5) {
  const rawStep = (half * 2) / (count - 1);
  const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const norm = rawStep / mag;
  const niceNorm = norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10;
  const step = niceNorm * mag;
  const maxK = Math.ceil(half / step);
  const ticks = [];
  for (let k = -maxK; k <= maxK; k++) {
    const v = +(k * step).toFixed(6);
    if (v >= -half - step * 1e-6 && v <= half + step * 1e-6) ticks.push(v);
  }
  return ticks;
}
 
function renderScatter(protein, metabolite, result) {
  const wrap = document.getElementById("scatter-wrap");
  const metIdx = state.scaled.metCols.indexOf(metabolite);
  const protIdx = state.scaled.protCols.indexOf(protein);
 
  if (metIdx === -1 || protIdx === -1) {
    wrap.innerHTML = `<div class="empty-state">Couldn't find per-sample values for this pair.</div>`;
    return;
  }
 
  const points = state.scaled.samples.map((sample, i) => {
    const x = state.scaled.metValues[i][metIdx];
    const y = state.scaled.protValues[i][protIdx];
    const target = state.scaled.targets[i];
    const quadrant = (x >= 0 && y >= 0) || (x <= 0 && y <= 0) ? "c" : "d";
    return { sample, target, x, y, quadrant };
  });
 
  const xs = points.map((p) => p.x);
  const ys = points.map((p) => p.y);
  const xHalf = Math.max(Math.abs(Math.min(...xs)), Math.abs(Math.max(...xs))) * 1.15 || 1;
  const yHalf = Math.max(Math.abs(Math.min(...ys)), Math.abs(Math.max(...ys))) * 1.15 || 1;
 
  // Layout
  const W = 640, H = 440;
  const marginLeft = 64, marginRight = 24, marginTop = 20, marginBottom = 56;
  const plotW = W - marginLeft - marginRight;
  const plotH = H - marginTop - marginBottom;
 
  const xScale = (x) => marginLeft + ((x + xHalf) / (2 * xHalf)) * plotW;
  const yScale = (y) => marginTop + (1 - (y + yHalf) / (2 * yHalf)) * plotH; // flip: up = smaller svg-y
 
  const circles = points.map((p) => {
    const cx = xScale(p.x).toFixed(1);
    const cy = yScale(p.y).toFixed(1);
    const cls = p.quadrant === "c" ? "pt-c" : "pt-d";
    const labelSvg = `<text class="point-label" x="${(+cx + 6).toFixed(1)}" y="${(+cy + 3).toFixed(1)}" style="${state.showLabels ? "" : "display:none;"}">${escapeHtml(p.target)}</text>`;
    return `<g>
      <circle class="${cls}" cx="${cx}" cy="${cy}" r="5"><title>${escapeHtml(p.target)} (${escapeHtml(p.sample)})\nMetabolite log2FC: ${fmtNum(p.x)}\nProtein log2FC: ${fmtNum(p.y)}</title></circle>
      ${labelSvg}
    </g>`;
  }).join("");
 
  const zeroX = xScale(0).toFixed(1);
  const zeroY = yScale(0).toFixed(1);
 
  const xTicks = niceTicks(xHalf);
  const yTicks = niceTicks(yHalf);
 
  const xTickSvg = xTicks.map((v) => {
    const tx = xScale(v).toFixed(1);
    const bottom = H - marginBottom;
    return `<line class="tick-line" x1="${tx}" y1="${bottom}" x2="${tx}" y2="${bottom + 6}" />
      <text class="tick-label" x="${tx}" y="${bottom + 18}" text-anchor="middle">${fmtNum(v, v === 0 ? 0 : 2)}</text>`;
  }).join("");
 
  const yTickSvg = yTicks.map((v) => {
    const ty = yScale(v).toFixed(1);
    return `<line class="tick-line" x1="${marginLeft - 6}" y1="${ty}" x2="${marginLeft}" y2="${ty}" />
      <text class="tick-label" x="${marginLeft - 10}" y="${(+ty + 3.5).toFixed(1)}" text-anchor="end">${fmtNum(v, v === 0 ? 0 : 2)}</text>`;
  }).join("");
 
  const svg = `
    <svg class="scatter-svg" viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg">
      <rect class="plot-frame" x="${marginLeft}" y="${marginTop}" width="${plotW}" height="${plotH}" />
 
      <line class="axis-line" x1="${marginLeft}" y1="${zeroY}" x2="${W - marginRight}" y2="${zeroY}" />
      <line class="axis-line" x1="${zeroX}" y1="${marginTop}" x2="${zeroX}" y2="${H - marginBottom}" />
 
      ${xTickSvg}
      ${yTickSvg}
 
      ${circles}
 
      <text class="annotation-text" x="${marginLeft + 10}" y="${marginTop + 16}">Concordance = ${fmtNum(result.concordance)}</text>
      <text class="annotation-text" x="${marginLeft + 10}" y="${marginTop + 32}">Adjusted p-value = ${fmtNum(result.padj, 4)}</text>
 
      <text class="axis-title" x="${W / 2}" y="${H - 14}" text-anchor="middle">Metabolite log2FC: ${escapeHtml(metabolite)}</text>
      <text class="axis-title" x="${16}" y="${H / 2}" text-anchor="middle" transform="rotate(-90 16 ${H / 2})">Protein log2FC: ${escapeHtml(protein)}</text>
    </svg>
  `;
 
  wrap.innerHTML = svg;
}
 
// ---------- init ----------
 
async function init() {
  renderNav("concordance-results");
 
  document.querySelectorAll("#mode-toggle button").forEach((btn) => {
    btn.addEventListener("click", () => setMode(btn.dataset.mode));
  });
 
  document.getElementById("analyte-search-input").addEventListener("input", (e) => {
    state.analyteQuery = e.target.value;
    refreshAnalyte();
  });
 
  document.getElementById("protein-input").addEventListener("input", maybeRenderPair);
  document.getElementById("metabolite-input").addEventListener("input", maybeRenderPair);
 
  document.getElementById("analyte-results-table").innerHTML = `<div class="loading-state">Loading concordance results…</div>`;
  try {
    state.conc = await loadDataFile("conc");
    buildPairIndex();
    populateDatalists();
    document.getElementById("analyte-results-table").innerHTML =
      `<div class="empty-state">Type a protein or metabolite name above to search.</div>`;
  } catch (err) {
    document.getElementById("analyte-results-table").innerHTML =
      `<div class="empty-state">Couldn't load the data file: ${escapeHtml(err.message)}</div>`;
    console.error(err);
  }
}
 
document.addEventListener("DOMContentLoaded", init);
   if (n === null || n === undefined || Number.isNaN(n)) return "NA";
  if (n === 0) return "0";
  if (Math.abs(n) < 0.001) return n.toExponential(2);
  return n.toFixed(digits);
}

// ---------- data loading (script-tag based, not fetch — see README for why) ----------

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

// ---------- mode toggle ----------

function setMode(mode) {
  state.mode = mode;
  document.querySelectorAll("#mode-toggle button").forEach((btn) => {
    btn.classList.toggle("active", btn.dataset.mode === mode);
  });
  document.getElementById("single-panel").style.display = mode === "single" ? "" : "none";
  document.getElementById("pair-panel").style.display = mode === "pair" ? "" : "none";

  if (mode === "pair" && !state.scaled) {
    document.getElementById("pair-hint").textContent = "Loading sample data…";
    loadDataFile("scaled").then((data) => {
      state.scaled = data;
      document.getElementById("pair-hint").textContent = "Select an exact protein and metabolite to see their result.";
      maybeRenderPair();
    }).catch((err) => {
      document.getElementById("pair-hint").textContent = "Couldn't load sample data: " + err.message;
    });
  }
}

// ---------- single analyte search ----------

function rowToObj(r) {
  return { protein: r[0], metabolite: r[1], concordance: r[2], variance: r[3], pvalue: r[4], padj: r[5], direction: r[6] };
}

function applyAnalyteFilter() {
  const q = state.analyteQuery.trim().toLowerCase();
  if (!q) { state.analyteFiltered = []; return; }
  state.analyteFiltered = state.conc.rows.filter(
    (r) => r[0].toLowerCase().includes(q) || r[1].toLowerCase().includes(q)
  );
}

function analyteColumns() {
  return [
    { key: "protein", label: "Protein", flex: 1.2, render: (r) => `<strong>${escapeHtml(r.protein)}</strong>` },
    { key: "metabolite", label: "Metabolite", flex: 1.4, render: (r) => escapeHtml(r.metabolite) },
    { key: "concordance", label: "Concordance", flex: 0.9, render: (r) => fmtNum(r.concordance) },
    { key: "variance", label: "Variance", flex: 0.9, render: (r) => fmtNum(r.variance, 4) },
    { key: "pvalue", label: "p-value", flex: 0.9, render: (r) => fmtNum(r.pvalue, 4) },
    { key: "padj", label: "p.adj", flex: 0.9, render: (r) => fmtNum(r.padj, 4) },
    { key: "direction", label: "Direction", flex: 0.9, render: (r) => `<span class="reg-pill reg-${r.direction}">${DIR_LABEL[r.direction] || r.direction}</span>` },
  ];
}

function renderAnalyteResults() {
  const container = document.getElementById("analyte-results-table");
  const countEl = document.getElementById("analyte-result-count");

  if (!state.analyteQuery.trim()) {
    container.innerHTML = `<div class="empty-state">Type a protein or metabolite name above to search.</div>`;
    countEl.textContent = "";
    state.analyteTable = null;
    return;
  }

  if (state.analyteFiltered.length === 0) {
    container.innerHTML = `<div class="empty-state">No rows match "${escapeHtml(state.analyteQuery)}".</div>`;
    countEl.textContent = "0 rows";
    state.analyteTable = null;
    return;
  }

  countEl.textContent = `${state.analyteFiltered.length.toLocaleString()} rows`;

  if (!state.analyteTable) {
    state.analyteTable = new VirtualTable({
      container,
      rowHeight: 40,
      columns: analyteColumns(),
      getRowCount: () => state.analyteFiltered.length,
      getRow: (i) => rowToObj(state.analyteFiltered[i]),
      rowClass: (row) => (row.direction === "c" ? "reg-c" : row.direction === "d" ? "reg-d" : ""),
    });
  }
  state.analyteTable.refresh();
}

// ---------- PMI pair search ----------

function buildPairIndex() {
  const map = new Map();
  const proteinSet = new Set();
  const metaboliteSet = new Set();
  for (const r of state.conc.rows) {
    map.set(r[0] + " " + r[1], r);
    proteinSet.add(r[0]);
    metaboliteSet.add(r[1]);
  }
  state.pairIndex = map;
  state.proteinSet = proteinSet;
  state.metaboliteSet = metaboliteSet;
}

function populateDatalists() {
  const proteinList = document.getElementById("protein-datalist");
  const metaboliteList = document.getElementById("metabolite-datalist");
  const proteins = [...state.proteinSet].sort();
  const metabolites = [...state.metaboliteSet].sort();
  proteinList.innerHTML = proteins.map((p) => `<option value="${escapeHtml(p)}"></option>`).join("");
  metaboliteList.innerHTML = metabolites.map((m) => `<option value="${escapeHtml(m)}"></option>`).join("");
}

function maybeRenderPair() {
  const proteinInput = document.getElementById("protein-input");
  const metaboliteInput = document.getElementById("metabolite-input");
  const hint = document.getElementById("pair-hint");
  const protein = proteinInput.value.trim();
  const metabolite = metaboliteInput.value.trim();

  const proteinValid = state.proteinSet && state.proteinSet.has(protein);
  const metaboliteValid = state.metaboliteSet && state.metaboliteSet.has(metabolite);

  proteinInput.classList.toggle("invalid", protein.length > 0 && !proteinValid);
  metaboliteInput.classList.toggle("invalid", metabolite.length > 0 && !metaboliteValid);

  if (!proteinValid || !metaboliteValid) {
    hint.textContent = "Select an exact protein and metabolite from the list.";
    hint.classList.remove("is-invalid");
    document.getElementById("pair-output").innerHTML =
      `<div class="empty-state">Choose a protein and a metabolite above.</div>`;
    return;
  }

  const row = state.pairIndex.get(protein + " " + metabolite);
  if (!row) {
    hint.textContent = `No concordance result found for ${protein} × ${metabolite}.`;
    hint.classList.add("is-invalid");
    document.getElementById("pair-output").innerHTML = `<div class="empty-state">No data for this pair.</div>`;
    return;
  }

  hint.textContent = "";
  hint.classList.remove("is-invalid");
  renderPairOutput(protein, metabolite, rowToObj(row));
}

function renderPairOutput(protein, metabolite, result) {
  const container = document.getElementById("pair-output");

  const statTiles = `
    <div class="stat-tile-row">
      <div class="stat-tile"><span class="stat-tile-label">Concordance</span><span class="stat-tile-value">${fmtNum(result.concordance)}</span></div>
      <div class="stat-tile"><span class="stat-tile-label">Variance</span><span class="stat-tile-value">${fmtNum(result.variance, 4)}</span></div>
      <div class="stat-tile"><span class="stat-tile-label">p-value</span><span class="stat-tile-value">${fmtNum(result.pvalue, 4)}</span></div>
      <div class="stat-tile"><span class="stat-tile-label">p.adj</span><span class="stat-tile-value">${fmtNum(result.padj, 4)}</span></div>
      <div class="stat-tile"><span class="stat-tile-label">Direction</span><span class="stat-tile-value"><span class="reg-pill reg-${result.direction}">${DIR_LABEL[result.direction] || result.direction}</span></span></div>
    </div>
  `;

  if (!state.scaled) {
    container.innerHTML = statTiles + `<div class="empty-state">Loading sample data for the quadrant plot…</div>`;
    return;
  }

  container.innerHTML = statTiles + `
    <div class="scatter-panel">
      <div class="scatter-panel-header">
        <span class="scatter-title">Concordance quadrant plot: ${escapeHtml(protein)} vs ${escapeHtml(metabolite)}</span>
        <label class="checkbox-field">
          <input type="checkbox" id="labels-toggle" ${state.showLabels ? "checked" : ""} />
          Show point labels
        </label>
      </div>
      <div class="scatter-wrap" id="scatter-wrap"></div>
      <div class="legend">
        <span class="legend-item"><span class="legend-swatch swatch-up"></span> Concordant</span>
        <span class="legend-item"><span class="legend-swatch swatch-down"></span> Discordant</span>
      </div>
    </div>
  `;

  renderScatter(protein, metabolite, result);

  document.getElementById("labels-toggle").addEventListener("change", (e) => {
    state.showLabels = e.target.checked;
    document.querySelectorAll(".scatter-svg .point-label").forEach((el) => {
      el.style.display = state.showLabels ? "" : "none";
    });
  });
}

// Pick "nice" round tick values (1/2/5 × a power of ten) spanning
// [-half, half], symmetric around 0 so 0 itself is always a tick.
function niceTicks(half, count = 5) {
  const rawStep = (half * 2) / (count - 1);
  const mag = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const norm = rawStep / mag;
  const niceNorm = norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10;
  const step = niceNorm * mag;
  const maxK = Math.ceil(half / step);
  const ticks = [];
  for (let k = -maxK; k <= maxK; k++) {
    const v = +(k * step).toFixed(6);
    if (v >= -half - step * 1e-6 && v <= half + step * 1e-6) ticks.push(v);
  }
  return ticks;
}

function renderScatter(protein, metabolite, result) {
  const wrap = document.getElementById("scatter-wrap");
  const metIdx = state.scaled.metCols.indexOf(metabolite);
  const protIdx = state.scaled.protCols.indexOf(protein);

  if (metIdx === -1 || protIdx === -1) {
    wrap.innerHTML = `<div class="empty-state">Couldn't find per-sample values for this pair.</div>`;
    return;
  }

  const points = state.scaled.samples.map((sample, i) => {
    const x = state.scaled.metValues[i][metIdx];
    const y = state.scaled.protValues[i][protIdx];
    const target = state.scaled.targets[i];
    const quadrant = (x >= 0 && y >= 0) || (x <= 0 && y <= 0) ? "c" : "d";
    return { sample, target, x, y, quadrant };
  });

  const xs = points.map((p) => p.x);
  const ys = points.map((p) => p.y);
  const xHalf = Math.max(Math.abs(Math.min(...xs)), Math.abs(Math.max(...xs))) * 1.15 || 1;
  const yHalf = Math.max(Math.abs(Math.min(...ys)), Math.abs(Math.max(...ys))) * 1.15 || 1;

  // Layout
  const W = 640, H = 440;
  const marginLeft = 64, marginRight = 24, marginTop = 20, marginBottom = 56;
  const plotW = W - marginLeft - marginRight;
  const plotH = H - marginTop - marginBottom;

  const xScale = (x) => marginLeft + ((x + xHalf) / (2 * xHalf)) * plotW;
  const yScale = (y) => marginTop + (1 - (y + yHalf) / (2 * yHalf)) * plotH; // flip: up = smaller svg-y

  const circles = points.map((p) => {
    const cx = xScale(p.x).toFixed(1);
    const cy = yScale(p.y).toFixed(1);
    const cls = p.quadrant === "c" ? "pt-c" : "pt-d";
    const labelSvg = `<text class="point-label" x="${(+cx + 6).toFixed(1)}" y="${(+cy + 3).toFixed(1)}" style="${state.showLabels ? "" : "display:none;"}">${escapeHtml(p.target)}</text>`;
    return `<g>
      <circle class="${cls}" cx="${cx}" cy="${cy}" r="5"><title>${escapeHtml(p.target)} (${escapeHtml(p.sample)})\nMetabolite log2FC: ${fmtNum(p.x)}\nProtein log2FC: ${fmtNum(p.y)}</title></circle>
      ${labelSvg}
    </g>`;
  }).join("");

  const zeroX = xScale(0).toFixed(1);
  const zeroY = yScale(0).toFixed(1);

  const xTicks = niceTicks(xHalf);
  const yTicks = niceTicks(yHalf);

  const xTickSvg = xTicks.map((v) => {
    const tx = xScale(v).toFixed(1);
    const bottom = H - marginBottom;
    return `<line class="tick-line" x1="${tx}" y1="${bottom}" x2="${tx}" y2="${bottom + 6}" />
      <text class="tick-label" x="${tx}" y="${bottom + 18}" text-anchor="middle">${fmtNum(v, v === 0 ? 0 : 2)}</text>`;
  }).join("");

  const yTickSvg = yTicks.map((v) => {
    const ty = yScale(v).toFixed(1);
    return `<line class="tick-line" x1="${marginLeft - 6}" y1="${ty}" x2="${marginLeft}" y2="${ty}" />
      <text class="tick-label" x="${marginLeft - 10}" y="${(+ty + 3.5).toFixed(1)}" text-anchor="end">${fmtNum(v, v === 0 ? 0 : 2)}</text>`;
  }).join("");

  const svg = `
    <svg class="scatter-svg" viewBox="0 0 ${W} ${H}" xmlns="http://www.w3.org/2000/svg">
      <rect class="plot-frame" x="${marginLeft}" y="${marginTop}" width="${plotW}" height="${plotH}" />

      <line class="axis-line" x1="${marginLeft}" y1="${zeroY}" x2="${W - marginRight}" y2="${zeroY}" />
      <line class="axis-line" x1="${zeroX}" y1="${marginTop}" x2="${zeroX}" y2="${H - marginBottom}" />

      ${xTickSvg}
      ${yTickSvg}

      ${circles}

      <text class="annotation-text" x="${marginLeft + 10}" y="${marginTop + 16}">Concordance = ${fmtNum(result.concordance)}</text>
      <text class="annotation-text" x="${marginLeft + 10}" y="${marginTop + 32}">Adjusted p-value = ${fmtNum(result.padj, 4)}</text>

      <text class="axis-title" x="${W / 2}" y="${H - 14}" text-anchor="middle">Metabolite log2FC: ${escapeHtml(metabolite)}</text>
      <text class="axis-title" x="${16}" y="${H / 2}" text-anchor="middle" transform="rotate(-90 16 ${H / 2})">Protein log2FC: ${escapeHtml(protein)}</text>
    </svg>
  `;

  wrap.innerHTML = svg;
}

// ---------- init ----------

async function init() {
  renderNav("concordance-results");

  document.querySelectorAll("#mode-toggle button").forEach((btn) => {
    btn.addEventListener("click", () => setMode(btn.dataset.mode));
  });

  document.getElementById("analyte-search-input").addEventListener("input", (e) => {
    state.analyteQuery = e.target.value;
    applyAnalyteFilter();
    renderAnalyteResults();
  });

  document.getElementById("protein-input").addEventListener("input", maybeRenderPair);
  document.getElementById("metabolite-input").addEventListener("input", maybeRenderPair);

  document.getElementById("analyte-results-table").innerHTML = `<div class="loading-state">Loading concordance results…</div>`;
  try {
    state.conc = await loadDataFile("conc");
    buildPairIndex();
    populateDatalists();
    document.getElementById("analyte-results-table").innerHTML =
      `<div class="empty-state">Type a protein or metabolite name above to search.</div>`;
  } catch (err) {
    document.getElementById("analyte-results-table").innerHTML =
      `<div class="empty-state">Couldn't load the data file: ${escapeHtml(err.message)}</div>`;
    console.error(err);
  }
}

document.addEventListener("DOMContentLoaded", init);
