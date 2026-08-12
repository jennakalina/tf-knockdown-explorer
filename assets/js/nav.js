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
      "Search differential metabolite (DAMs) and protein (DAPs) abundance results across TF knockdown RNAi lines, by target gene, metabolite, or protein.",
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
