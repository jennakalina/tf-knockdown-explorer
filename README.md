# TF Knockdown Data Explorer

A static website for browsing TF knockdown datasets. No build step, no
server-side code — just HTML/CSS/JS and JSON data files, so it can be hosted
directly on GitHub Pages.

## What's here

```
index.html                        Home page — links to every tab
differential-analysis.html        "Browse Differential Analysis Results" tab
assets/css/style.css              Shared styles (light + dark mode)
assets/js/nav.js                  Tab registry — add new tabs here
assets/js/virtual-table.js        Reusable virtualized table (handles large row counts smoothly)
assets/js/differential-analysis.js  Page logic for the differential analysis tab
data/dams.data.js                 Differential metabolite results (17,150 rows)
data/daps.data.js                 Differential protein results (95,680 rows)
```

**Note on the data files:** they're `.js` files (not `.json`) loaded via a plain
`<script src="...">` tag rather than `fetch()`. That's deliberate — `fetch()`
of a local file is blocked by the browser's CORS policy when you open the
HTML directly from disk (`file://...`), which is what caused search to show
nothing when testing locally before deploying. A `<script src>` load isn't
subject to that restriction, so the site now works identically whether you
open `differential-analysis.html` straight from disk or serve it over
GitHub Pages. If you ever want to test with a local server anyway (optional,
not required): `python3 -m http.server` from this folder, then visit
`http://localhost:8000`.

## Deploying to GitHub Pages

1. Create a new repository on GitHub (or reuse an existing one), e.g. `tf-knockdown-explorer`.
2. Copy all the files in this folder into the repository (keep the folder structure as-is).
3. Commit and push:
   ```
   git init
   git add .
   git commit -m "Initial TF knockdown data explorer"
   git branch -M main
   git remote add origin https://github.com/<your-username>/<your-repo>.git
   git push -u origin main
   ```
4. On GitHub, go to the repo's **Settings → Pages**.
5. Under "Build and deployment", set **Source** to "Deploy from a branch", pick the `main` branch and `/ (root)` folder, then **Save**.
6. GitHub will give you a URL like `https://<your-username>.github.io/<your-repo>/` — the site is live there within a minute or two.

No build tools, npm install, or CI step is needed — everything is static and works as-is once pushed.

## Adding another tab later

1. Add a new HTML page at the root (copy `differential-analysis.html` as a starting point).
2. Add one entry to the `SITE_TABS` array in `assets/js/nav.js` — the new tab link will automatically appear in the top nav on every page and as a card on the home page.
3. Drop any new data file(s) into `data/`.

## About the data

Both `dams.data.js` and `daps.data.js` are compact conversions of the
original `DAMs_all.csv` and `daps_for_interface.csv` files, one row per
RNAi line × metabolite (or protein) combination, sorted by adjusted p-value
(most significant first). Each file just assigns a plain JS object to
`window.__TF_DATA.dams` / `window.__TF_DATA.daps`. Each row has:

| field | meaning |
|---|---|
| `rnai` | RNAi line ID |
| `name` | metabolite or protein name |
| `log2fc` | log2 fold change |
| `pvalue` | raw p-value |
| `padj` | adjusted p-value |
| `target` | the TF gene knocked down |
| `reg` | `up`, `down`, or `ns` (not significant) |

The page loads whichever dataset is selected (DAMs or DAPs) as a `<script>`
tag the first time it's needed, then searches/filters entirely in the
browser — no backend required. If you regenerate these files from updated
CSVs, keep the same `window.__TF_DATA.<key> = {...}` wrapper and column
order, or update `assets/js/differential-analysis.js` to match.
