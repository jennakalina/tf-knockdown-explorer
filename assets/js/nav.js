/* ============================================================
   Single source of truth for the site's tabs.
   Add a new dataset/tab by adding one entry here — it will
   automatically show up in the top nav AND on the home page.
   ============================================================ */

const SITE_TABS = [
  {
    id: "differential-analysis",
    label: "Browse Differential Analysis Results",
    href: "differential-analysis.html",
    description:
      "Search differentially abundant metabolite and protein results by knockdown target, metabolite/protein name, or RNAi line.",
  },
  {
    id: "concordance-results",
    label: "Browse Concordance Results",
    href: "concordance-results.html",
    description:
      "Search concordance results by single analyte or by PMI (protein-metabolite interaction) pair.",
  },
  {
    id: "perturbation-clusters",
    label: "Perturbation Clusters",
    href: "perturbation-clusters.html",
    description:
      "Explore perturbation clusters of TF knockdown lines, including per-cluster DAM and DAP log2FC heatmaps.",
  },
  // Add future tabs here, e.g.:
  // {
  //   id: "another-dataset",
  //   label: "Browse Another Dataset",
  //   href: "another-dataset.html",
  //   description: "One-line description shown on the home page.",
  // },
];

function renderNav(activeId) {
  const nav = document.getElementById("site-nav");
  if (!nav) return;

  const links = SITE_TABS.map(
    (tab) =>
      `<a href="${tab.href}" class="nav-link${
        tab.id === activeId ? " active" : ""
      }">${tab.label}</a>`
  ).join("");

  nav.innerHTML = `
    <a href="index.html" class="nav-brand">TF Knockdown Data Explorer</a>
    <div class="nav-links">${links}</div>
  `;
}
