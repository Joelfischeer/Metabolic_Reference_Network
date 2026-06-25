import json
import re
import networkx as nx


# ---------------------------------------------------------------------------
# Organ colour palette
# ---------------------------------------------------------------------------

ORGAN_COLORS = {
    "Brain": "#7c3aed",
    "Heart": "#dc2626",
    "Liver": "#92400e",
    "Kidney": "#ea580c",
    "Muscle": "#1d4ed8",
    "Pancreas": "#0f766e",
    "WAT": "#b45309",
    "Thyroid": "#15803d",
    "Adrenal Glands": "#9f1239",
    "Lung": "#0284c7",
    "Spleen": "#6b21a8",
    "Bone Marrow": "#374151",
    "Small Intestine": "#c2410c",
    "Colon": "#78350f",
}
DEFAULT_NODE_COLOR = "#2563eb"


def _parse_edge_text(text: str) -> dict:
    """Extract structured fields from free-text edge description."""
    result = {"connection_type": "", "key_players": [], "notes": "", "sources": []}
    if not text:
        return result

    m = re.search(r"Type:\s*(.+?)(?:\n|$)", text, re.IGNORECASE)
    if m:
        result["connection_type"] = m.group(1).strip()

    m = re.search(r"Key Players:\s*(.+?)(?:\n|$)", text, re.IGNORECASE)
    if m:
        raw = m.group(1).strip()
        result["key_players"] = [p.strip() for p in raw.split(",") if p.strip()]

    m = re.search(r"Notes?:\s*(.*?)(?:Sources?:|$)", text, re.IGNORECASE | re.DOTALL)
    if m:
        result["notes"] = m.group(1).strip()

    m = re.search(r"Sources?:\s*(.*?)$", text, re.IGNORECASE | re.DOTALL)
    if m:
        raw_src = m.group(1).strip()
        result["sources"] = [
            s.strip().lstrip("- ").strip()
            for s in re.split(r"\n-|\n", raw_src)
            if s.strip().lstrip("- ").strip()
        ]

    return result


def export_network_to_cytoscape_dashboard(
    graph: nx.Graph,
    filename: str = "network_dashboard.html",
    directed: bool | None = None,
    include_legend: bool = True,
    title: str = "Metabolic Reference Network",
):
    """
    Export a NetworkX graph to an interactive Cytoscape.js dashboard.

    Features
    --------
    - Click edge → sidebar with key-player chips, notes, sources, PubMed papers
    - Click node → sidebar with organ description + connected edges summary
    - Search bar — type any molecule/keyword → matching edges highlight, results listed
    - Layout switcher (force / circle / grid)
    - Organ-specific node colours
    - Edge thickness reflects evidence strength
    """

    if directed is None:
        directed = graph.is_directed()

    elements = []

    # --- Nodes ---
    for node, attrs in graph.nodes(data=True):
        color = attrs.get("color") or ORGAN_COLORS.get(str(node), DEFAULT_NODE_COLOR)
        elements.append({
            "data": {
                "id": str(node),
                "label": str(node),
                "description": attrs.get("description", ""),
                "color": color,
            }
        })

    # --- Edges ---
    for u, v, attrs in graph.edges(data=True):
        raw_desc = attrs.get("description", "")
        parsed = _parse_edge_text(raw_desc)

        # Enrich from merged_data if it exists on the edge
        merged = attrs.get("merged_data", {})
        key_players = (
            merged.get("key_players_raw")
            or parsed["key_players"]
        )
        all_kp_categories = merged.get("key_players_merged", {})
        all_kp_counts     = merged.get("key_players_counts", {})
        pubmed = merged.get("pubmed", {})

        connection_type = (
            merged.get("connection_type")
            or parsed["connection_type"]
        )
        notes = merged.get("notes") or parsed["notes"]
        sources = merged.get("sources") or parsed["sources"]
        ai_description = merged.get("ai_description", "")

        edge_color = attrs.get("color", "#64748b")

        elements.append({
            "data": {
                "id": f"{u}__{v}",
                "source": str(u),
                "target": str(v),
                "color": edge_color,
                "description": raw_desc,
                "connection_type": connection_type,
                "key_players": key_players,
                "key_players_metabolites": all_kp_categories.get("metabolites", []),
                "key_players_hormones": all_kp_categories.get("hormones", []),
                "key_players_proteins": all_kp_categories.get("proteins", []),
                "key_players_counts_metabolites": all_kp_counts.get("metabolites", {}),
                "key_players_counts_hormones": all_kp_counts.get("hormones", {}),
                "key_players_counts_proteins": all_kp_counts.get("proteins", {}),
                "notes": notes,
                "sources": sources,
                "ai_description": ai_description,
                "pubmed_n": pubmed.get("n_papers", 0),
                "pubmed_papers": pubmed.get("papers", []),
                "pubmed_query": pubmed.get("query", ""),
            }
        })

    elements_json = json.dumps(elements, ensure_ascii=False)
    arrow_shape = "triangle" if directed else "none"

    legend_html = ""
    if include_legend:
        legend_html = """
<div id="legend">
  <div class="legend-title">Edge Legend</div>
  <div class="legend-item"><span class="legend-dot" style="background:#22c55e"></span>In both networks</div>
  <div class="legend-item"><span class="legend-dot" style="background:#ef4444"></span>Missing in input</div>
  <div class="legend-item"><span class="legend-dot" style="background:#f97316"></span>Only in reference</div>
</div>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<script src="https://unpkg.com/cytoscape@3.29.2/dist/cytoscape.min.js"></script>
<style>
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}

  body {{
    font-family: 'Inter', system-ui, -apple-system, sans-serif;
    background: #0f172a;
    color: #e2e8f0;
    height: 100vh;
    display: flex;
    flex-direction: column;
    overflow: hidden;
  }}

  /* ---- Top bar ---- */
  #topbar {{
    display: flex;
    align-items: center;
    gap: 12px;
    padding: 10px 20px;
    background: #1e293b;
    border-bottom: 1px solid #334155;
    flex-shrink: 0;
    z-index: 10;
  }}

  #topbar h1 {{
    font-size: 16px;
    font-weight: 700;
    color: #f1f5f9;
    white-space: nowrap;
  }}

  #search-wrap {{
    position: relative;
    flex: 1;
    max-width: 400px;
  }}

  #search-icon {{
    position: absolute;
    left: 10px;
    top: 50%;
    transform: translateY(-50%);
    color: #94a3b8;
    pointer-events: none;
    font-size: 14px;
  }}

  #search-input {{
    width: 100%;
    padding: 7px 12px 7px 32px;
    background: #0f172a;
    border: 1px solid #334155;
    border-radius: 8px;
    color: #e2e8f0;
    font-size: 14px;
    outline: none;
    transition: border-color 0.2s;
  }}
  #search-input:focus {{ border-color: #6366f1; }}
  #search-input::placeholder {{ color: #64748b; }}

  #search-count {{
    font-size: 12px;
    color: #64748b;
    white-space: nowrap;
  }}

  .ctrl-btn {{
    padding: 6px 12px;
    background: #334155;
    border: none;
    border-radius: 6px;
    color: #cbd5e1;
    font-size: 13px;
    cursor: pointer;
    transition: background 0.15s, color 0.15s;
    white-space: nowrap;
  }}
  .ctrl-btn:hover {{ background: #475569; color: #f1f5f9; }}
  .ctrl-btn.active {{ background: #6366f1; color: white; }}

  /* ---- Main area ---- */
  #main {{
    display: flex;
    flex: 1;
    overflow: hidden;
  }}

  /* ---- Graph canvas ---- */
  #cy-wrap {{
    flex: 1;
    position: relative;
    overflow: hidden;
  }}

  #cy {{
    width: 100%;
    height: 100%;
  }}

  /* ---- Sidebar ---- */
  #sidebar {{
    width: 340px;
    min-width: 260px;
    max-width: 480px;
    background: #1e293b;
    border-left: 1px solid #334155;
    display: flex;
    flex-direction: column;
    overflow: hidden;
    flex-shrink: 0;
  }}

  #sidebar-header {{
    padding: 14px 16px 10px;
    border-bottom: 1px solid #334155;
    flex-shrink: 0;
  }}

  #sidebar-title {{
    font-size: 13px;
    font-weight: 600;
    color: #94a3b8;
    text-transform: uppercase;
    letter-spacing: 0.06em;
    margin-bottom: 4px;
  }}

  #sidebar-name {{
    font-size: 17px;
    font-weight: 700;
    color: #f1f5f9;
    word-break: break-word;
  }}

  #sidebar-body {{
    flex: 1;
    overflow-y: auto;
    padding: 14px 16px;
  }}

  #sidebar-body::-webkit-scrollbar {{ width: 5px; }}
  #sidebar-body::-webkit-scrollbar-track {{ background: transparent; }}
  #sidebar-body::-webkit-scrollbar-thumb {{ background: #334155; border-radius: 4px; }}

  .section-label {{
    font-size: 11px;
    font-weight: 600;
    color: #64748b;
    text-transform: uppercase;
    letter-spacing: 0.07em;
    margin: 14px 0 6px;
  }}
  .section-label:first-child {{ margin-top: 0; }}

  .connection-type-badge {{
    display: inline-block;
    padding: 3px 10px;
    border-radius: 999px;
    font-size: 12px;
    font-weight: 600;
    background: #312e81;
    color: #a5b4fc;
    margin-bottom: 2px;
  }}

  .chips-wrap {{
    display: flex;
    flex-wrap: wrap;
    gap: 5px;
  }}

  .chip {{
    display: inline-flex;
    align-items: center;
    padding: 3px 9px;
    border-radius: 999px;
    font-size: 12px;
    font-weight: 500;
    cursor: default;
    transition: opacity 0.15s;
  }}
  .chip:hover {{ opacity: 0.85; }}

  .chip-hormone   {{ background: #451a03; color: #fdba74; border: 1px solid #92400e; }}
  .chip-metabolite {{ background: #042f2e; color: #5eead4; border: 1px solid #0f766e; }}
  .chip-protein   {{ background: #1e1b4b; color: #a5b4fc; border: 1px solid #3730a3; }}
  .chip-general   {{ background: #1e293b; color: #94a3b8; border: 1px solid #334155; }}

  .notes-text {{
    font-size: 13px;
    color: #cbd5e1;
    line-height: 1.6;
    white-space: pre-wrap;
  }}

  .source-item {{
    font-size: 12px;
    color: #6366f1;
    margin-bottom: 3px;
    word-break: break-all;
    display: block;
    text-decoration: none;
  }}
  .source-item:hover {{ color: #818cf8; }}

  .paper-card {{
    background: #0f172a;
    border: 1px solid #334155;
    border-radius: 8px;
    padding: 10px;
    margin-bottom: 8px;
  }}
  .paper-title {{
    font-size: 13px;
    font-weight: 600;
    color: #e2e8f0;
    margin-bottom: 4px;
    line-height: 1.4;
  }}
  .paper-meta {{
    font-size: 11px;
    color: #64748b;
    margin-bottom: 4px;
  }}
  .paper-abstract {{
    font-size: 12px;
    color: #94a3b8;
    line-height: 1.5;
    display: -webkit-box;
    -webkit-line-clamp: 3;
    -webkit-box-orient: vertical;
    overflow: hidden;
  }}
  .paper-link {{
    font-size: 11px;
    color: #6366f1;
    text-decoration: none;
    display: inline-block;
    margin-top: 4px;
  }}
  .paper-link:hover {{ color: #818cf8; }}

  .empty-state {{
    text-align: center;
    padding: 40px 20px;
    color: #475569;
    font-size: 14px;
    line-height: 1.6;
  }}
  .empty-state svg {{ margin: 0 auto 12px; display: block; opacity: 0.4; }}

  /* ---- Search results panel ---- */
  #search-results {{
    display: none;
    position: absolute;
    top: 52px;
    left: 0;
    right: 0;
    background: #1e293b;
    border-bottom: 1px solid #334155;
    max-height: 200px;
    overflow-y: auto;
    z-index: 20;
    padding: 8px 0;
  }}

  .search-result-item {{
    padding: 8px 20px;
    cursor: pointer;
    font-size: 13px;
    color: #cbd5e1;
    display: flex;
    align-items: center;
    gap: 8px;
    transition: background 0.1s;
  }}
  .search-result-item:hover {{ background: #334155; }}
  .search-result-item .match-type {{ font-size: 11px; color: #64748b; }}

  /* ---- Legend ---- */
  #legend {{
    position: absolute;
    bottom: 16px;
    left: 16px;
    background: #1e293b;
    border: 1px solid #334155;
    border-radius: 10px;
    padding: 10px 14px;
    font-size: 12px;
    z-index: 5;
  }}
  .legend-title {{ font-weight: 600; color: #94a3b8; margin-bottom: 6px; font-size: 11px; text-transform: uppercase; letter-spacing: 0.05em; }}
  .legend-item {{ display: flex; align-items: center; gap: 6px; margin-bottom: 3px; color: #cbd5e1; }}
  .legend-dot {{ width: 10px; height: 10px; border-radius: 50%; flex-shrink: 0; }}

  /* ---- Hover tooltip ---- */
  #tooltip {{
    position: fixed;
    display: none;
    background: #1e293b;
    border: 1px solid #334155;
    padding: 6px 10px;
    border-radius: 6px;
    font-size: 12px;
    color: #e2e8f0;
    pointer-events: none;
    z-index: 100;
    max-width: 220px;
    box-shadow: 0 4px 16px rgba(0,0,0,0.5);
  }}

  /* ---- Highlighted search match ---- */
  .search-matched {{ border: 2px solid #fbbf24 !important; }}
</style>
</head>
<body>

<!-- Top bar -->
<div id="topbar">
  <h1>🧬 {title}</h1>

  <div id="search-wrap">
    <span id="search-icon">🔍</span>
    <input id="search-input" type="text"
           placeholder="Search molecule, hormone, organ…"
           autocomplete="off">
    <div id="search-results"></div>
  </div>

  <span id="search-count"></span>

  <button class="ctrl-btn" onclick="runLayout('cose')" title="Force-directed layout">⚡ Force</button>
  <button class="ctrl-btn" onclick="runLayout('circle')" title="Circle layout">⭕ Circle</button>
  <button class="ctrl-btn" onclick="runLayout('grid')" title="Grid layout">⊞ Grid</button>
  <button class="ctrl-btn" onclick="cy.fit()" title="Fit all">⤢ Fit</button>
</div>

<!-- Main area -->
<div id="main">
  <div id="cy-wrap">
    <div id="cy"></div>
    {legend_html}
  </div>

  <div id="sidebar">
    <div id="sidebar-header">
      <div id="sidebar-title">Select a node or edge</div>
      <div id="sidebar-name">—</div>
    </div>
    <div id="sidebar-body">
      <div class="empty-state">
        <svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
          <circle cx="12" cy="12" r="10"/>
          <path d="M12 8v4M12 16h.01"/>
        </svg>
        Click any <strong>node</strong> or <strong>edge</strong> to see details.<br><br>
        Use the search bar above to find connections by molecule name.
      </div>
    </div>
  </div>
</div>

<div id="tooltip"></div>

<script>
const elements = {elements_json};

// ── Build lookup maps ─────────────────────────────────────────────────────
const edgeMap = {{}};   // id → element data
const nodeMap = {{}};

elements.forEach(el => {{
  if (el.data.source) edgeMap[el.data.id] = el.data;
  else nodeMap[el.data.id] = el.data;
}});

// ── Cytoscape init ────────────────────────────────────────────────────────
const cy = cytoscape({{
  container: document.getElementById('cy'),
  elements,
  style: [
    {{
      selector: 'node',
      style: {{
        'background-color': 'data(color)',
        'label': 'data(label)',
        'color': '#f1f5f9',
        'text-valign': 'center',
        'text-halign': 'center',
        'font-size': '11px',
        'font-weight': '600',
        'width': 52,
        'height': 52,
        'border-width': 2,
        'border-color': 'rgba(255,255,255,0.15)',
        'text-outline-color': 'rgba(0,0,0,0.6)',
        'text-outline-width': 2,
        'opacity': 1,
        'transition-property': 'opacity, border-color, width, height',
        'transition-duration': '0.15s',
      }}
    }},
    {{
      selector: 'node:selected, node.highlighted',
      style: {{
        'border-width': 3,
        'border-color': '#fbbf24',
        'width': 60,
        'height': 60,
      }}
    }},
    {{
      selector: 'edge',
      style: {{
        'width': 2.5,
        'line-color': 'data(color)',
        'target-arrow-color': 'data(color)',
        'target-arrow-shape': '{arrow_shape}',
        'curve-style': 'bezier',
        'opacity': 1,
        'transition-property': 'opacity, width, line-color',
        'transition-duration': '0.15s',
      }}
    }},
    {{
      selector: 'edge:selected, edge.highlighted',
      style: {{
        'width': 5,
        'line-color': '#fbbf24',
        'target-arrow-color': '#fbbf24',
        'opacity': 1,
      }}
    }},
    {{
      selector: '.faded',
      style: {{ 'opacity': 0.08 }}
    }},
    {{
      selector: '.search-hit',
      style: {{
        'line-color': '#fbbf24',
        'target-arrow-color': '#fbbf24',
        'width': 4,
        'opacity': 1,
      }}
    }},
    {{
      selector: 'node.search-hit',
      style: {{
        'border-color': '#fbbf24',
        'border-width': 3,
        'opacity': 1,
      }}
    }},
  ],
  layout: {{
    name: 'cose',
    animate: true,
    padding: 80,
    nodeRepulsion: 15000,
    idealEdgeLength: 180,
    edgeElasticity: 200,
    gravity: 0.25,
    randomize: false,
  }}
}});

// ── Helpers ───────────────────────────────────────────────────────────────

function runLayout(name) {{
  cy.layout({{ name, animate: true, padding: 60 }}).run();
}}

function chipsHtml(items, cls, counts) {{
  if (!items || items.length === 0) return '<span style="color:#475569;font-size:12px">None identified</span>';
  return '<div class="chips-wrap">' +
    items.map(t => {{
      const n = counts && counts[t];
      const badge = n ? `<span style="margin-left:4px;background:rgba(0,0,0,0.25);border-radius:8px;padding:0 5px;font-size:10px;font-weight:700">${{n}}×</span>` : '';
      return `<span class="chip ${{cls}}" title="${{n ? n + ' mentions' : ''}}">${{escHtml(t)}}${{badge}}</span>`;
    }}).join('') +
    '</div>';
}}

function escHtml(s) {{
  return String(s)
    .replace(/&/g, '&amp;').replace(/</g, '&lt;')
    .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}}

function section(label, content) {{
  return `<div class="section-label">${{label}}</div>${{content}}`;
}}

// ── Sidebar rendering ─────────────────────────────────────────────────────

function showEdgeSidebar(data) {{
  document.getElementById('sidebar-title').textContent = 'Connection';
  document.getElementById('sidebar-name').textContent =
    data.source + ' ↔ ' + data.target;

  let html = '';

  if (data.connection_type) {{
    html += section('Type',
      `<span class="connection-type-badge">${{escHtml(data.connection_type)}}</span>`);
  }}

  // Key players — hormones, metabolites, proteins
  const hormones    = data.key_players_hormones    || [];
  const metabolites = data.key_players_metabolites || [];
  const proteins    = data.key_players_proteins    || [];
  const rawKP       = data.key_players             || [];
  const cntH = data.key_players_counts_hormones    || {{}};
  const cntM = data.key_players_counts_metabolites || {{}};
  const cntP = data.key_players_counts_proteins    || {{}};

  // If we have categorised players, show them; otherwise fall back to raw list
  if (hormones.length || metabolites.length || proteins.length) {{
    html += section('Key Players',
      '<div style="display:flex;flex-direction:column;gap:8px">' +
      (hormones.length    ? '<div><span style="font-size:11px;color:#64748b;margin-bottom:3px;display:block">Hormones</span>'              + chipsHtml(hormones,    'chip-hormone',    cntH) + '</div>' : '') +
      (metabolites.length ? '<div><span style="font-size:11px;color:#64748b;margin-bottom:3px;display:block">Metabolites</span>'           + chipsHtml(metabolites, 'chip-metabolite', cntM) + '</div>' : '') +
      (proteins.length    ? '<div><span style="font-size:11px;color:#64748b;margin-bottom:3px;display:block">Proteins / Transporters</span>' + chipsHtml(proteins,    'chip-protein',    cntP) + '</div>' : '') +
      '</div>'
    );
  }} else if (rawKP.length) {{
    html += section('Key Players', chipsHtml(rawKP, 'chip-general', {{}}));
  }}

  if (data.ai_description) {{
    // Turn [PMID 12345678] into clickable PubMed links
    const linkedSummary = escHtml(data.ai_description).replace(
      /\[PMID\s+(\d+)\]/gi,
      (_, pmid) => `<a href="https://pubmed.ncbi.nlm.nih.gov/${{pmid}}/" target="_blank" rel="noopener"
         style="color:#818cf8;text-decoration:underline dotted">[PMID ${{pmid}}]</a>`
    );
    html += section('Literature Summary',
      `<div style="line-height:1.7;border-left:3px solid #6366f1;padding-left:10px;color:#cbd5e1;font-size:13px">${{linkedSummary}}</div>`);
  }}

  if (data.notes) {{
    html += section('Notes', `<p class="notes-text">${{escHtml(data.notes)}}</p>`);
  }}

  if (data.sources && data.sources.length) {{
    html += section('Sources',
      data.sources.filter(s => s.trim()).map(s => {{
        const href = s.match(/https?:[/][/][^ ]+/) ? s.match(/https?:[/][/][^ ]+/)[0] : '#';
        return `<a class="source-item" href="${{href}}" target="_blank" rel="noopener">${{escHtml(s)}}</a>`;
      }}).join('')
    );
  }}

  // PubMed section
  if (data.pubmed_n > 0) {{
    let papersHtml = `<div style="font-size:12px;color:#64748b;margin-bottom:8px">`+
      `Found ${{data.pubmed_n}} papers in last 5 years</div>`;
    if (data.pubmed_papers && data.pubmed_papers.length) {{
      papersHtml += data.pubmed_papers.slice(0, 4).map(p => {{
        const link = p.doi ? `https://doi.org/${{p.doi}}` : `https://pubmed.ncbi.nlm.nih.gov/${{p.pmid}}/`;
        return `<div class="paper-card">
          <div class="paper-title">${{escHtml(p.title)}}</div>
          <div class="paper-meta">PMID: ${{p.pmid}} · ${{p.year}}</div>
          ${{p.abstract ? `<div class="paper-abstract">${{escHtml(p.abstract)}}</div>` : ''}}
          <a class="paper-link" href="${{link}}" target="_blank" rel="noopener">→ View on PubMed</a>
        </div>`;
      }}).join('');
    }}
    html += section('PubMed Literature', papersHtml);
  }}

  document.getElementById('sidebar-body').innerHTML = html || '<p style="color:#475569;font-size:13px">No additional information available.</p>';
}}

function showNodeSidebar(data) {{
  document.getElementById('sidebar-title').textContent = 'Organ';
  document.getElementById('sidebar-name').textContent = data.label;

  let html = '';
  if (data.description) {{
    html += section('Description',
      `<p class="notes-text">${{escHtml(data.description)}}</p>`);
  }}

  // List connected edges
  const node = cy.getElementById(data.id);
  const edges = node.connectedEdges();
  if (edges.length) {{
    html += section(`Connections (${{edges.length}})`,
      '<div style="display:flex;flex-direction:column;gap:4px">' +
      edges.map(e => {{
        const other = e.source().id() === data.id ? e.target().id() : e.source().id();
        const ct = e.data('connection_type') || '';
        return `<div style="display:flex;justify-content:space-between;align-items:center;padding:5px 8px;background:#0f172a;border-radius:6px;cursor:pointer"
          onclick="selectEdge('${{e.id()}}')" >
          <span style="font-size:13px;color:#cbd5e1">${{escHtml(other)}}</span>
          ${{ct ? `<span style="font-size:11px;color:#64748b">${{escHtml(ct)}}</span>` : ''}}
        </div>`;
      }}).join('') +
      '</div>'
    );
  }}

  document.getElementById('sidebar-body').innerHTML = html ||
    '<p style="color:#475569;font-size:13px">No description available.</p>';
}}

function selectEdge(edgeId) {{
  cy.elements().removeClass('faded highlighted');
  const edge = cy.getElementById(edgeId);
  if (edge.length) {{
    cy.elements().addClass('faded');
    edge.removeClass('faded').addClass('highlighted');
    edge.connectedNodes().removeClass('faded').addClass('highlighted');
    showEdgeSidebar(edge.data());
  }}
}}

// ── Click events ──────────────────────────────────────────────────────────

cy.on('tap', 'edge', function(evt) {{
  const el = evt.target;
  cy.elements().removeClass('faded highlighted search-hit');
  cy.elements().addClass('faded');
  el.removeClass('faded').addClass('highlighted');
  el.connectedNodes().removeClass('faded').addClass('highlighted');
  showEdgeSidebar(el.data());
}});

cy.on('tap', 'node', function(evt) {{
  const el = evt.target;
  cy.elements().removeClass('faded highlighted search-hit');
  cy.elements().addClass('faded');
  el.removeClass('faded').addClass('highlighted');
  el.connectedEdges().removeClass('faded').addClass('highlighted');
  el.connectedNodes().removeClass('faded');
  showNodeSidebar(el.data());
}});

cy.on('tap', function(evt) {{
  if (evt.target === cy) {{
    cy.elements().removeClass('faded highlighted search-hit');
    document.getElementById('sidebar-title').textContent = 'Select a node or edge';
    document.getElementById('sidebar-name').textContent = '—';
    document.getElementById('sidebar-body').innerHTML = `<div class="empty-state">
      <svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
        <circle cx="12" cy="12" r="10"/><path d="M12 8v4M12 16h.01"/>
      </svg>
      Click any <strong>node</strong> or <strong>edge</strong> to see details.
    </div>`;
  }}
}});

// ── Tooltip ───────────────────────────────────────────────────────────────

const tooltip = document.getElementById('tooltip');

cy.on('mouseover', 'node', function(evt) {{
  const d = evt.target.data();
  tooltip.textContent = d.label;
  tooltip.style.display = 'block';
}});

cy.on('mouseover', 'edge', function(evt) {{
  const d = evt.target.data();
  const kp = (d.key_players || []).slice(0, 4).join(', ');
  tooltip.innerHTML = `<strong>${{escHtml(d.source)}} ↔ ${{escHtml(d.target)}}</strong>${{kp ? '<br><span style="color:#94a3b8">' + escHtml(kp) + '</span>' : ''}}`;
  tooltip.style.display = 'block';
}});

cy.on('mousemove', function(evt) {{
  tooltip.style.left = (evt.originalEvent.clientX + 12) + 'px';
  tooltip.style.top = (evt.originalEvent.clientY + 12) + 'px';
}});

cy.on('mouseout', 'node, edge', function() {{
  tooltip.style.display = 'none';
}});

// ── Search ────────────────────────────────────────────────────────────────

const searchInput = document.getElementById('search-input');
const searchResults = document.getElementById('search-results');
const searchCount = document.getElementById('search-count');

function buildSearchIndex() {{
  const index = [];
  cy.edges().forEach(e => {{
    const d = e.data();
    const text = [
      d.description || '',
      (d.key_players || []).join(' '),
      (d.key_players_hormones || []).join(' '),
      (d.key_players_metabolites || []).join(' '),
      (d.key_players_proteins || []).join(' '),
      d.notes || '',
      d.connection_type || '',
    ].join(' ').toLowerCase();
    index.push({{ el: e, text, label: d.source + ' ↔ ' + d.target, type: 'edge' }});
  }});
  cy.nodes().forEach(n => {{
    const d = n.data();
    const text = (d.description || '').toLowerCase() + ' ' + d.label.toLowerCase();
    index.push({{ el: n, text, label: d.label, type: 'node' }});
  }});
  return index;
}}

const searchIndex = buildSearchIndex();

searchInput.addEventListener('input', function() {{
  const q = this.value.trim().toLowerCase();

  if (!q) {{
    cy.elements().removeClass('faded search-hit');
    searchResults.style.display = 'none';
    searchCount.textContent = '';
    return;
  }}

  // Find matches
  const hits = searchIndex.filter(entry => entry.text.includes(q));
  const hitEls = cy.collection(hits.map(h => h.el));
  const edgeHits = hits.filter(h => h.type === 'edge');
  const nodeHits = hits.filter(h => h.type === 'node');

  // Visual feedback
  cy.elements().removeClass('search-hit').addClass('faded');
  hitEls.removeClass('faded').addClass('search-hit');
  // Also un-fade nodes connected to hit edges
  edgeHits.forEach(h => h.el.connectedNodes().removeClass('faded'));

  searchCount.textContent = `${{edgeHits.length}} edge${{edgeHits.length !== 1 ? 's' : ''}}`;

  // Dropdown results list
  if (hits.length === 0) {{
    searchResults.innerHTML = '<div class="search-result-item" style="color:#64748b">No results found</div>';
  }} else {{
    searchResults.innerHTML = hits.slice(0, 20).map(h => {{
      const icon = h.type === 'edge' ? '🔗' : '⬤';
      return `<div class="search-result-item" data-id="${{h.el.id()}}" data-type="${{h.type}}">
        <span>${{icon}}</span>
        <span>${{escHtml(h.label)}}</span>
        <span class="match-type">${{h.type}}</span>
      </div>`;
    }}).join('');
    searchResults.querySelectorAll('.search-result-item[data-id]').forEach(item => {{
      item.addEventListener('click', () => {{
        const id = item.dataset.id;
        const type = item.dataset.type;
        searchResults.style.display = 'none';
        if (type === 'edge') {{
          selectEdge(id);
        }} else {{
          cy.getElementById(id).trigger('tap');
        }}
        cy.center(cy.getElementById(id));
      }});
    }});
  }}
  searchResults.style.display = 'block';
}});

// Hide dropdown on outside click
document.addEventListener('click', e => {{
  if (!document.getElementById('search-wrap').contains(e.target)) {{
    searchResults.style.display = 'none';
  }}
}});
searchInput.addEventListener('focus', () => {{
  if (searchInput.value.trim()) searchResults.style.display = 'block';
}});

</script>
</body>
</html>"""

    with open(filename, "w", encoding="utf-8") as f:
        f.write(html)

    graph_type = "Directed" if directed else "Undirected"
    print(f"[✔] Dashboard saved to {filename} ({graph_type})")
