/* ============================================================
   Minimal virtualized table.
   Renders only the rows currently in (or near) the viewport,
   so browsing tens of thousands of rows stays smooth without
   pagination — the scrollbar reflects the true full row count.
   ============================================================ */
 
class VirtualTable {
  /**
   * @param {Object} opts
   * @param {HTMLElement} opts.container - scrollable viewport element
   * @param {number} [opts.rowHeight]
   * @param {Array<{key:string, label:string, flex?:number, render?:(row)=>string, sortable?:boolean}>} opts.columns
   * @param {() => number} opts.getRowCount
   * @param {(index:number) => any} opts.getRow
   * @param {(row:any) => string} [opts.rowClass]
   * @param {(key:string) => void} [opts.onSort] - called with a column key when its (sortable) header is clicked
   * @param {{key:string, dir:'asc'|'desc'}} [opts.sortState] - initial sort indicator to display
   */
  constructor({ container, rowHeight = 38, columns, getRowCount, getRow, rowClass, onSort, sortState }) {
    this.container = container;
    this.rowHeight = rowHeight;
    this.columns = columns;
    this.getRowCount = getRowCount;
    this.getRow = getRow;
    this.rowClass = rowClass || (() => "");
    this.onSort = onSort || null;
    this.sortState = sortState || null;
    this._build();
  }
 
  _build() {
    this.container.innerHTML = "";
    this.container.classList.add("vtable-viewport");
 
    const table = document.createElement("div");
    table.className = "vtable";
 
    this.header = document.createElement("div");
    this.header.className = "vtable-header";
    this._renderHeader();
    if (this.onSort) {
      this.header.addEventListener("click", (e) => {
        const cell = e.target.closest(".vtable-th");
        if (!cell) return;
        this.onSort(cell.dataset.key);
      });
    }
    table.appendChild(this.header);
 
    const body = document.createElement("div");
    body.className = "vtable-body";
    table.appendChild(body);
 
    this.container.appendChild(table);
    this.body = body;
 
    this.spacer = document.createElement("div");
    this.spacer.className = "vtable-spacer";
    body.appendChild(this.spacer);
 
    this.rowsLayer = document.createElement("div");
    this.rowsLayer.className = "vtable-rows";
    this.rowsLayer.style.position = "absolute";
    this.rowsLayer.style.top = "0";
    this.rowsLayer.style.left = "0";
    this.rowsLayer.style.right = "0";
    body.appendChild(this.rowsLayer);
 
    this._onScroll = () => this._render();
    this._onResize = () => this._render();
    this.container.addEventListener("scroll", this._onScroll);
    window.addEventListener("resize", this._onResize);
  }
 
  _renderHeader() {
    this.header.innerHTML = this.columns
      .map((c) => {
        const sortable = this.onSort && c.sortable !== false;
        const isActive = this.sortState && this.sortState.key === c.key;
        const arrow = isActive ? (this.sortState.dir === "asc" ? " ▲" : " ▼") : "";
        const cls = "vtable-cell" + (sortable ? " vtable-th" : "") + (isActive ? " vtable-th-active" : "");
        return `<div class="${cls}" style="flex:${c.flex || 1}"${sortable ? ` data-key="${c.key}"` : ""}>${c.label}<span class="sort-arrow">${arrow}</span></div>`;
      })
      .join("");
  }
 
  // Update the header's sort indicator (e.g. after the page re-sorts its
  // own data and calls refresh()) without rebuilding the whole table.
  setSort(key, dir) {
    this.sortState = key ? { key, dir } : null;
    this._renderHeader();
  }
 
  refresh() {
    const total = this.getRowCount();
    this.spacer.style.height = `${total * this.rowHeight}px`;
    // Reset scroll to top when the underlying data set changes shape
    // (e.g. a new search) so we don't render past the end.
    if (this.container.scrollTop > total * this.rowHeight) {
      this.container.scrollTop = 0;
    }
    this._render();
  }
 
  _render() {
    const total = this.getRowCount();
    const viewportHeight = this.container.clientHeight || 1;
    const scrollTop = this.container.scrollTop;
    const buffer = 10;
 
    const startIdx = Math.max(0, Math.floor(scrollTop / this.rowHeight) - buffer);
    const visibleCount = Math.ceil(viewportHeight / this.rowHeight) + buffer * 2;
    const endIdx = Math.min(total, startIdx + visibleCount);
 
    this.rowsLayer.style.transform = `translateY(${startIdx * this.rowHeight}px)`;
 
    let html = "";
    for (let i = startIdx; i < endIdx; i++) {
      const row = this.getRow(i);
      if (!row) continue;
      const cls = this.rowClass(row);
      html += `<div class="vtable-row ${cls}" style="height:${this.rowHeight}px">`;
      for (const c of this.columns) {
        const content = c.render ? c.render(row) : row[c.key];
        html += `<div class="vtable-cell" style="flex:${c.flex || 1}">${content}</div>`;
      }
      html += `</div>`;
    }
    this.rowsLayer.innerHTML = html;
  }
 
  destroy() {
    this.container.removeEventListener("scroll", this._onScroll);
    window.removeEventListener("resize", this._onResize);
  }
}
     header.innerHTML = this.columns
      .map((c) => `<div class="vtable-cell" style="flex:${c.flex || 1}">${c.label}</div>`)
      .join("");
    table.appendChild(header);

    const body = document.createElement("div");
    body.className = "vtable-body";
    table.appendChild(body);

    this.container.appendChild(table);
    this.body = body;

    this.spacer = document.createElement("div");
    this.spacer.className = "vtable-spacer";
    body.appendChild(this.spacer);

    this.rowsLayer = document.createElement("div");
    this.rowsLayer.className = "vtable-rows";
    this.rowsLayer.style.position = "absolute";
    this.rowsLayer.style.top = "0";
    this.rowsLayer.style.left = "0";
    this.rowsLayer.style.right = "0";
    body.appendChild(this.rowsLayer);

    this._onScroll = () => this._render();
    this._onResize = () => this._render();
    this.container.addEventListener("scroll", this._onScroll);
    window.addEventListener("resize", this._onResize);
  }

  refresh() {
    const total = this.getRowCount();
    this.spacer.style.height = `${total * this.rowHeight}px`;
    // Reset scroll to top when the underlying data set changes shape
    // (e.g. a new search) so we don't render past the end.
    if (this.container.scrollTop > total * this.rowHeight) {
      this.container.scrollTop = 0;
    }
    this._render();
  }

  _render() {
    const total = this.getRowCount();
    const viewportHeight = this.container.clientHeight || 1;
    const scrollTop = this.container.scrollTop;
    const buffer = 10;

    const startIdx = Math.max(0, Math.floor(scrollTop / this.rowHeight) - buffer);
    const visibleCount = Math.ceil(viewportHeight / this.rowHeight) + buffer * 2;
    const endIdx = Math.min(total, startIdx + visibleCount);

    this.rowsLayer.style.transform = `translateY(${startIdx * this.rowHeight}px)`;

    let html = "";
    for (let i = startIdx; i < endIdx; i++) {
      const row = this.getRow(i);
      if (!row) continue;
      const cls = this.rowClass(row);
      html += `<div class="vtable-row ${cls}" style="height:${this.rowHeight}px">`;
      for (const c of this.columns) {
        const content = c.render ? c.render(row) : row[c.key];
        html += `<div class="vtable-cell" style="flex:${c.flex || 1}">${content}</div>`;
      }
      html += `</div>`;
    }
    this.rowsLayer.innerHTML = html;
  }

  destroy() {
    this.container.removeEventListener("scroll", this._onScroll);
    window.removeEventListener("resize", this._onResize);
  }
}
