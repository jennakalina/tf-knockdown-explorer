/* ============================================================
   Page logic for differential-analysis.html
   Loads the DAMs (metabolites) and DAPs (proteins) result sets
   and wires up dataset toggle, search, "significant only" filter,
   and the virtualized results table.
   ============================================================ */

const REG_LABEL = { up: "Up", down: "Down", ns: "Not sig.", na: "NA" };

const DATASETS = {
  dams: { url: "data/dams.data.js", label: "Metabolites (DAMs)", nameCol: "Metabolite" },
  daps: { url: "data/daps.data.js", label: "Proteins (DAPs)", nameCol: "Protein" },
};

const state = {
  active: "dams",
  data: {}, // { dams: {columns, rows}, daps: {...} }
  filtered: [], // currently-displayed rows for the active dataset
  query: "",
  sigOnly: false,
  table: null,
};

function fmtNum(n, digits = 3) {
  if (n === null || n === undefined || Number.isNaN(n)) return "NA";
  if (n === 0) return "0";
  if (Math.abs(n) < 0.001) return n.toExponential(2);
  return n.toFixed(digits);
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;",
  }[c]));
}

// Data files are loaded as plain <script> tags (not fetch/XHR) so the site
// works when opened directly from disk (file://) as well as when hosted —
// fetch() of local files is blocked by CORS under file://, but a <script src>
// load is not subject to that restriction.
function loadDataset(key) {
  if (state.data[key]) return Promise.resolve(state.data[key]);

  const globalKey = key; // matches window.__TF_DATA.<key> set by data/*.data.js
  window.__TF_DATA = window.__TF_DATA || {};
  if (window.__TF_DATA[globalKey]) {
    state.data[key] = window.__TF_DATA[globalKey];
    return Promise.resolve(state.data[key]);
  }

  return new Promise((resolve, reject) => {
    const cfg = DATASETS[key];
    const script = document.createElement("script");
    script.src = cfg.url;
    script.onload = () => {
      const data = window.__TF_DATA && window.__TF_DATA[globalKey];
      if (!data) {
        reject(new Error(`Loaded ${cfg.url} but no data found for "${globalKey}".`));
        return;
      }
      state.data[key] = data;
      resolve(data);
    };
    script.onerror = () => reject(new Error(`Failed to load data file: ${cfg.url}`));
    document.head.appendChild(script);
  });
}

function applyFilters() {
  const ds = state.data[state.active];
  if (!ds) return;

  const tokens = state.query.trim().toLowerCase().split(/\s+/).filter(Boolean);
  const idx = { rnai: 0, name: 1, log2fc: 2, pvalue: 3, padj: 4, target: 5, reg: 6 };

  let rows = ds.rows;

  if (tokens.length) {
    rows = rows.filter((r) => {
      const haystack = (r[idx.target] + " " + r[idx.name] + " " + r[idx.rnai]).toLowerCase();
      return tokens.every((t) => haystack.includes(t));
    });
  }

  if (state.sigOnly) {
    rows = rows.filter((r) => r[idx.reg] === "up" || r[idx.reg] === "down");
  }

  state.filtered = rows;
}

function rowToObj(r) {
  return { rnai: r[0], name: r[1], log2fc: r[2], pvalue: r[3], padj: r[4], target: r[5], reg: r[6] };
}

function buildColumns() {
  const nameLabel = DATASETS[state.active].nameCol;
  return [
    { key: "target", label: "Target (TF)", flex: 1.1, render: (r) => `<strong>${escapeHtml(r.target)}</strong>` },
    { key: "name", label: nameLabel, flex: 1.6, render: (r) => escapeHtml(r.name) },
    { key: "rnai", label: "RNAi line", flex: 1.3, render: (r) => escapeHtml(r.rnai) },
    { key: "log2fc", label: "log2FC", flex: 0.8, render: (r) => fmtNum(r.log2fc) },
    { key: "pvalue", label: "p-value", flex: 0.9, render: (r) => fmtNum(r.pvalue, 4) },
    { key: "padj", label: "p.adj", flex: 0.9, render: (r) => fmtNum(r.padj, 4) },
    {
      key: "reg",
      label: "Regulation",
      flex: 0.9,
      render: (r) => `<span class="reg-pill reg-${r.reg}">${REG_LABEL[r.reg] || r.reg}</span>`,
    },
  ];
}

function updateResultCount() {
  const el = document.getElementById("result-count");
  const total = state.data[state.active] ? state.data[state.active].rows.length : 0;
  el.textContent = `${state.filtered.length.toLocaleString()} of ${total.toLocaleString()} rows`;
}

function renderTable() {
  const container = document.getElementById("results-table");
  if (state.filtered.length === 0) {
    container.innerHTML = `<div class="empty-state">No rows match your search.</div>`;
    state.table = null;
    return;
  }
  if (!state.table || state.table._active !== state.active) {
    state.table = new VirtualTable({
      container,
      rowHeight: 40,
      columns: buildColumns(),
      getRowCount: () => state.filtered.length,
      getRow: (i) => rowToObj(state.filtered[i]),
      rowClass: (row) => (row.reg === "up" ? "reg-up" : row.reg === "down" ? "reg-down" : ""),
    });
    state.table._active = state.active;
  }
  state.table.refresh();
}

function refresh() {
  applyFilters();
  updateResultCount();
  renderTable();
}

function showLoadError(err) {
  const container = document.getElementById("results-table");
  container.innerHTML = `<div class="empty-state">Couldn't load the data file: ${escapeHtml(err.message)}</div>`;
  console.error(err);
}

function setActiveDataset(key) {
  state.active = key;
  document.querySelectorAll(".dataset-toggle button").forEach((btn) => {
    btn.classList.toggle("active", btn.dataset.key === key);
  });
  const container = document.getElementById("results-table");
  container.innerHTML = `<div class="loading-state">Loading ${DATASETS[key].label}…</div>`;
  loadDataset(key).then(refresh).catch(showLoadError);
}

async function init() {
  renderNav("differential-analysis");

  document.querySelectorAll(".dataset-toggle button").forEach((btn) => {
    btn.addEventListener("click", () => setActiveDataset(btn.dataset.key));
  });

  const searchInput = document.getElementById("search-input");
  searchInput.addEventListener("input", (e) => {
    state.query = e.target.value;
    refresh();
  });

  const sigCheckbox = document.getElementById("sig-only");
  sigCheckbox.addEventListener("change", (e) => {
    state.sigOnly = e.target.checked;
    refresh();
  });

  document.getElementById("results-table").innerHTML =
    `<div class="loading-state">Loading ${DATASETS[state.active].label}…</div>`;
  try {
    await loadDataset(state.active);
    refresh();
  } catch (err) {
    showLoadError(err);
  }
}

document.addEventListener("DOMContentLoaded", init);
