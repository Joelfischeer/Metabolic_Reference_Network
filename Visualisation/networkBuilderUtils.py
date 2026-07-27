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
    start_layout: str = "cose",
    info_panel_tabs: list | None = None,
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
                "id":              str(node),
                "label":           str(node),
                "description":     attrs.get("description", ""),
                "llm_description": attrs.get("llm_description", ""),
                "llm_papers":      attrs.get("llm_papers", []),
                "color":           color,
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

    # ── Info panel tabs ──────────────────────────────────────────────────────
    _default_tabs = [
        {
            "id": "search", "label": "Literature Search",
            "content": """
        <div class="info-h2">Overview</div>
        <p class="info-p">
          For each organ pair (and each individual organ), the pipeline queries the
          <strong>PubMed / NCBI E-utilities API</strong> to retrieve recent research
          papers on their metabolic interaction. All searches are cached locally so
          the pipeline can be resumed without repeating completed queries.
        </p>
        <div class="info-h2">Search Parameters</div>
        <div class="info-stat-grid">
          <div class="info-stat"><div class="info-stat-val">10 yrs</div><div class="info-stat-lbl">Look-back window (curated edges)</div></div>
          <div class="info-stat"><div class="info-stat-val">5 yrs</div><div class="info-stat-lbl">Look-back window (general network)</div></div>
          <div class="info-stat"><div class="info-stat-val">200</div><div class="info-stat-lbl">Max papers fetched per pair</div></div>
          <div class="info-stat" id="stat-edges-with-papers"><div class="info-stat-val" id="stat-val-edges">—</div><div class="info-stat-lbl">Edges with ≥1 paper (this network)</div></div>
        </div>
        <div class="info-h2">Query Construction</div>
        <p class="info-p">Each organ has a set of <strong>aliases</strong> and an <strong>MeSH term</strong>. For a pair A–B:</p>
        <div class="info-code">(MeSH_A OR aliases_A) AND (MeSH_B OR aliases_B)
AND (metabolism OR metabolic OR substrate
    OR glucose OR "fatty acid" OR insulin
    OR energy OR oxidation OR hormone
    OR gluconeogenesis)</div>
        <p class="info-p">If fewer than 5 papers are returned, a <strong>cascade fallback</strong> broadens the query.</p>
        <div class="info-h2">Key Player Extraction</div>
        <p class="info-p">All abstracts are scanned against three curated vocabulary lists:</p>
        <ul style="margin:0 0 10px 16px;padding:0;color:#cbd5e1">
          <li style="margin-bottom:4px"><strong style="color:#fdba74">Hormones</strong> — insulin, glucagon, cortisol, leptin, GLP-1, IGF-1, FGF21, irisin, …</li>
          <li style="margin-bottom:4px"><strong style="color:#5eead4">Metabolites</strong> — glucose, lactate, fatty acids, bile acids, amino acids, ATP, …</li>
          <li style="margin-bottom:4px"><strong style="color:#a5b4fc">Proteins</strong> — GLUT1–5, AMPK, mTOR, PPARs, PGC-1α, CPT1, LPL, …</li>
        </ul>
        <p class="info-p">Key players are <strong>ranked by mention count</strong> and shown as chips. Use the <em>Key Players</em> threshold to filter by count.</p>
""",
        },
        {
            "id": "llm", "label": "LLM Descriptions",
            "content": """
        <div class="info-h2">Model</div>
        <p class="info-p">Descriptions are generated by a <strong>locally running Ollama model</strong> (<code>north-mini-code-1.0</code>) via <code>http://localhost:11434</code>.</p>
        <div class="info-h2">Edge (Axis) Descriptions</div>
        <p class="info-p">For each organ–organ pair the top <strong>5 papers</strong> are sent to the model:</p>
        <div class="info-code">Write a concise scientific summary (3–4 sentences)
explaining the metabolic and hormonal basis of the
[Organ A]–[Organ B] axis based on the papers below.
Cite papers inline: [1] or [2,3].
Write flowing prose — no bullet points, no headings,
no closing notes or meta-commentary.</div>
        <div class="info-h2">Organ Descriptions</div>
        <div class="info-code">Write exactly 5 sentences describing the metabolic
function of the [Organ]. Focus on substrates consumed
or produced and key hormones that regulate its metabolism.</div>
        <div class="info-h2">Caching &amp; Resumability</div>
        <p class="info-p">Every description is saved after generation so the pipeline can be interrupted and resumed without losing progress.</p>
""",
        },
    ]
    tabs = info_panel_tabs if info_panel_tabs is not None else _default_tabs
    info_tabs_html = "\n      ".join(
        f'<button class="info-tab{"  active" if i == 0 else ""}" '
        f'onclick="showInfoSection(\'{t["id"]}\')">{t["label"]}</button>'
        for i, t in enumerate(tabs)
    )
    info_sections_html = "\n      ".join(
        f'<div id="info-section-{t["id"]}" class="info-section{"  active" if i == 0 else ""}">'
        f'{t["content"]}</div>'
        for i, t in enumerate(tabs)
    )

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

  /* ---- Info panel (left) ---- */
  #info-panel {{
    width: 0;
    min-width: 0;
    background: #1e293b;
    border-right: 1px solid #334155;
    display: flex;
    flex-direction: column;
    overflow: hidden;
    flex-shrink: 0;
    transition: width 0.25s ease, min-width 0.25s ease;
  }}
  #info-panel.open {{
    width: 360px;
    min-width: 260px;
  }}

  #info-panel-header {{
    padding: 14px 16px 10px;
    border-bottom: 1px solid #334155;
    flex-shrink: 0;
    display: flex;
    align-items: flex-start;
    justify-content: space-between;
    gap: 8px;
  }}

  #info-panel-title {{
    font-size: 17px;
    font-weight: 700;
    color: #f1f5f9;
  }}

  #info-panel-close {{
    background: none;
    border: none;
    color: #64748b;
    font-size: 18px;
    cursor: pointer;
    line-height: 1;
    padding: 2px 4px;
    flex-shrink: 0;
  }}
  #info-panel-close:hover {{ color: #f1f5f9; }}

  #info-tabs {{
    display: flex;
    gap: 0;
    padding: 10px 16px 0;
    border-bottom: 1px solid #334155;
    flex-shrink: 0;
  }}

  .info-tab {{
    background: none;
    border: none;
    border-bottom: 2px solid transparent;
    color: #64748b;
    font-size: 12px;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.06em;
    padding: 6px 10px 8px;
    cursor: pointer;
    transition: color 0.15s, border-color 0.15s;
    white-space: nowrap;
  }}
  .info-tab:hover {{ color: #cbd5e1; }}
  .info-tab.active {{ color: #a5b4fc; border-bottom-color: #6366f1; }}

  #info-body {{
    flex: 1;
    overflow-y: auto;
    padding: 16px;
    font-size: 13px;
    color: #cbd5e1;
    line-height: 1.65;
  }}
  #info-body::-webkit-scrollbar {{ width: 5px; }}
  #info-body::-webkit-scrollbar-track {{ background: transparent; }}
  #info-body::-webkit-scrollbar-thumb {{ background: #334155; border-radius: 4px; }}

  .info-section {{ display: none; }}
  .info-section.active {{ display: block; }}

  .info-h2 {{
    font-size: 12px;
    font-weight: 700;
    color: #94a3b8;
    text-transform: uppercase;
    letter-spacing: 0.07em;
    margin: 18px 0 6px;
  }}
  .info-h2:first-child {{ margin-top: 0; }}

  .info-p {{
    margin: 0 0 10px;
    color: #cbd5e1;
  }}

  .info-stat-grid {{
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 8px;
    margin-bottom: 12px;
  }}

  .info-stat {{
    background: #0f172a;
    border: 1px solid #334155;
    border-radius: 8px;
    padding: 10px 12px;
  }}
  .info-stat-val {{
    font-size: 22px;
    font-weight: 700;
    color: #a5b4fc;
    line-height: 1.1;
  }}
  .info-stat-lbl {{
    font-size: 11px;
    color: #64748b;
    margin-top: 3px;
  }}

  .info-code {{
    background: #0f172a;
    border: 1px solid #334155;
    border-radius: 6px;
    padding: 8px 10px;
    font-family: monospace;
    font-size: 11px;
    color: #7dd3fc;
    margin: 6px 0 10px;
    white-space: pre-wrap;
    word-break: break-all;
  }}

  /* ---- info button in topbar ---- */
  #info-btn {{
    display: flex;
    align-items: center;
    justify-content: center;
    width: 30px;
    height: 30px;
    border-radius: 50%;
    background: #1e293b;
    border: 1px solid #334155;
    color: #94a3b8;
    font-size: 14px;
    font-weight: 700;
    cursor: pointer;
    flex-shrink: 0;
    transition: border-color 0.2s, color 0.2s, background 0.2s;
  }}
  #info-btn:hover {{ border-color: #6366f1; color: #a5b4fc; }}
  #info-btn.active {{ background: #1e1b4b; border-color: #6366f1; color: #a5b4fc; }}

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

  /* ---- Cytoscape cursor passthrough ---- */
  #cy {{ cursor: default; }}

  /* ---- Highlighted search match ---- */
  .search-matched {{ border: 2px solid #fbbf24 !important; }}

  /* ---- Key Players filter ---- */
  #kp-filter-wrap {{
    position: relative;
  }}

  #kp-filter-btn {{
    display: flex;
    align-items: center;
    gap: 6px;
    padding: 6px 12px;
    background: #1e293b;
    border: 1px solid #334155;
    border-radius: 8px;
    color: #cbd5e1;
    font-size: 13px;
    cursor: pointer;
    white-space: nowrap;
    transition: border-color 0.2s, color 0.2s;
  }}
  #kp-filter-btn:hover {{ border-color: #6366f1; color: #e2e8f0; }}
  #kp-filter-btn.active {{ border-color: #6366f1; color: #a5b4fc; background: #1e1b4b; }}

  #kp-filter-dropdown {{
    display: none;
    position: absolute;
    top: calc(100% + 6px);
    right: 0;
    background: #1e293b;
    border: 1px solid #334155;
    border-radius: 10px;
    padding: 14px 16px;
    z-index: 200;
    min-width: 220px;
    box-shadow: 0 8px 24px rgba(0,0,0,0.5);
  }}
  #kp-filter-dropdown.open {{ display: block; }}

  #kp-filter-dropdown label {{
    display: block;
    font-size: 12px;
    color: #94a3b8;
    margin-bottom: 8px;
    font-weight: 500;
    text-transform: uppercase;
    letter-spacing: 0.05em;
  }}

  #kp-filter-dropdown .kp-input-row {{
    display: flex;
    align-items: center;
    gap: 8px;
  }}

  #kp-threshold-input {{
    flex: 1;
    padding: 6px 10px;
    background: #0f172a;
    border: 1px solid #334155;
    border-radius: 6px;
    color: #e2e8f0;
    font-size: 14px;
    width: 80px;
  }}
  #kp-threshold-input:focus {{ outline: none; border-color: #6366f1; }}

  #kp-filter-dropdown .kp-hint {{
    margin-top: 8px;
    font-size: 11px;
    color: #475569;
    line-height: 1.4;
  }}

  #kp-threshold-badge {{
    display: none;
    background: #6366f1;
    color: #fff;
    border-radius: 10px;
    padding: 1px 6px;
    font-size: 11px;
    font-weight: 700;
    margin-left: 2px;
  }}
</style>
</head>
<body>

<!-- Top bar -->
<div id="topbar">
  <button id="info-btn" onclick="toggleInfoPanel()" title="How it works">i</button>

  <h1>{title}</h1>

  <div id="search-wrap">
    <span id="search-icon">🔍</span>
    <input id="search-input" type="text"
           placeholder="Search molecule, hormone, organ…"
           autocomplete="off">
    <div id="search-results"></div>
  </div>

  <span id="search-count"></span>

  <div id="kp-filter-wrap">
    <button id="kp-filter-btn" onclick="toggleKpDropdown(event)" title="Filter key players by minimum mention count">
      Key Players <span id="kp-threshold-badge"></span> ▾
    </button>
    <div id="kp-filter-dropdown">
      <label>Min. mentions per key player</label>
      <div class="kp-input-row">
        <input id="kp-threshold-input" type="number" min="0" step="1" value="0"
               placeholder="0 = show all">
        <button class="ctrl-btn" onclick="resetKpThreshold()" title="Reset to 0">✕ Reset</button>
      </div>
      <div class="kp-hint">Key players below threshold are hidden. Search only highlights edges where the molecule meets this count.</div>
    </div>
  </div>

  <button class="ctrl-btn" onclick="runLayout('circle')" title="Circle layout">⭕ Circle</button>
  <button class="ctrl-btn" onclick="runLayout('grid')" title="Grid layout">⊞ Grid</button>
  <button class="ctrl-btn" onclick="cy.fit()" title="Fit all">⤢ Fit</button>
</div>

<!-- Main area -->
<div id="main">

  <!-- Info panel (left) -->
  <div id="info-panel">
    <div id="info-panel-header">
      <div id="info-panel-title">How It Works</div>
      <button id="info-panel-close" onclick="toggleInfoPanel()" title="Close">✕</button>
    </div>
    <div id="info-tabs">
      {info_tabs_html}
    </div>
    <div id="info-body">
      {info_sections_html}
    </div><!-- /info-body -->
  </div><!-- /info-panel -->

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
      selector: 'node.hovered',
      style: {{
        'border-width': 3,
        'border-color': '#94a3b8',
        'width': 58,
        'height': 58,
      }}
    }},
    {{
      selector: 'edge.hovered',
      style: {{
        'width': 5,
        'line-color': '#94a3b8',
        'target-arrow-color': '#94a3b8',
        'opacity': 1,
      }}
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
    name: '{start_layout}',
    animate: {str(start_layout != 'circle').lower()},
    padding: 80,
    nodeRepulsion: 15000,
    idealEdgeLength: 180,
    edgeElasticity: 200,
    gravity: 0.25,
    randomize: false,
  }}
}});

// ── Helpers ───────────────────────────────────────────────────────────────

let kpThreshold = 0;
let currentEdgeData = null;  // last edge displayed in sidebar

function runLayout(name) {{
  cy.layout({{ name, animate: true, padding: 60 }}).run();
}}

function chipsHtml(items, cls, counts) {{
  if (!items || items.length === 0) return '<span style="color:#475569;font-size:12px">None identified</span>';
  const visible = kpThreshold > 0
    ? items.filter(t => (counts && counts[t] || 0) >= kpThreshold)
    : items;
  if (visible.length === 0) return `<span style="color:#475569;font-size:12px">None with ≥${{kpThreshold}} mentions</span>`;
  return '<div class="chips-wrap">' +
    visible.map(t => {{
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
  currentEdgeData = data;
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
    const papers = data.pubmed_papers || [];

    // Replace [1], [2,3], [1-3] etc. with superscript links and collect cited indices
    const citedIndices = new Set();
    const linkedSummary = escHtml(data.ai_description).replace(
      /\[(\d[\d,\s-]*)\]/g,
      (match, inner) => {{
        // Parse individual numbers out of "1,2" or "1-3" or "1, 2"
        const nums = [];
        inner.split(/[,\s]+/).forEach(part => {{
          const dash = part.match(/^(\d+)-(\d+)$/);
          if (dash) {{
            for (let i = parseInt(dash[1]); i <= parseInt(dash[2]); i++) nums.push(i);
          }} else {{
            const n = parseInt(part);
            if (!isNaN(n)) nums.push(n);
          }}
        }});
        nums.forEach(n => citedIndices.add(n));
        const links = nums.map(n => {{
          const p = papers[n - 1];
          const href = p ? (p.doi ? `https://doi.org/${{p.doi}}` : `https://pubmed.ncbi.nlm.nih.gov/${{p.pmid}}/`) : '#';
          return `<a href="${{href}}" target="_blank" rel="noopener"
            style="color:#818cf8;text-decoration:none;font-size:10px;vertical-align:super;font-weight:700">${{n}}</a>`;
        }}).join(',');
        return `[${{links}}]`;
      }}
    );

    // Build bibliography only for cited papers
    const bibEntries = [];
    citedIndices.forEach(n => {{
      const p = papers[n - 1];
      if (!p) return;
      const href = p.doi ? `https://doi.org/${{p.doi}}` : `https://pubmed.ncbi.nlm.nih.gov/${{p.pmid}}/`;
      bibEntries.push({{ n, p, href }});
    }});
    bibEntries.sort((a, b) => a.n - b.n);

    const bibHtml = bibEntries.length ? bibEntries.map(entry =>
      `<div style="display:flex;gap:6px;margin-bottom:6px;font-size:11px;color:#94a3b8">
        <span style="color:#818cf8;font-weight:700;flex-shrink:0">${{entry.n}}.</span>
        <span><a href="${{entry.href}}" target="_blank" rel="noopener"
          style="color:#cbd5e1;text-decoration:underline dotted">${{escHtml(entry.p.title)}}</a>
          <span style="color:#475569"> · ${{entry.p.year}} · PMID ${{entry.p.pmid}}</span></span>
      </div>`
    ).join('') : '';

    const summaryBlock =
      `<div style="line-height:1.7;color:#cbd5e1;font-size:13px;margin-bottom:${{bibHtml ? '12px' : '0'}}">${{linkedSummary}}</div>` +
      (bibHtml ? `<div style="border-top:1px solid #334155;padding-top:8px">${{bibHtml}}</div>` : '');

    html += section('Literature Summary',
      `<div style="border-left:3px solid #6366f1;padding-left:10px">${{summaryBlock}}</div>`);
  }}


  // Show total paper count as a compact note (no cards)
  if (data.pubmed_n > 0 && !data.ai_description) {{
    html += section('PubMed Literature',
      `<div style="font-size:12px;color:#64748b">Found ${{data.pubmed_n}} papers — run LLM step to generate summary.</div>`);
  }}

  document.getElementById('sidebar-body').innerHTML = html || '<p style="color:#475569;font-size:13px">No additional information available.</p>';
}}

function showNodeSidebar(data) {{
  currentEdgeData = null;
  document.getElementById('sidebar-title').textContent = 'Organ';
  document.getElementById('sidebar-name').textContent = data.label;

  let html = '';

  if (data.llm_description) {{
    const papers = data.llm_papers || [];
    const citedIndices = new Set();
    const linkedSummary = escHtml(data.llm_description).replace(
      /\[(\d[\d,\s-]*)\]/g,
      (match, inner) => {{
        const nums = [];
        inner.split(/[,\s]+/).forEach(part => {{
          const dash = part.match(/^(\d+)-(\d+)$/);
          if (dash) {{
            for (let i = parseInt(dash[1]); i <= parseInt(dash[2]); i++) nums.push(i);
          }} else {{
            const n = parseInt(part);
            if (!isNaN(n)) nums.push(n);
          }}
        }});
        nums.forEach(n => citedIndices.add(n));
        const links = nums.map(n => {{
          const p = papers[n - 1];
          const href = p ? (p.doi ? `https://doi.org/${{p.doi}}` : `https://pubmed.ncbi.nlm.nih.gov/${{p.pmid}}/`) : '#';
          return `<a href="${{href}}" target="_blank" rel="noopener"
            style="color:#818cf8;text-decoration:none;font-size:10px;vertical-align:super;font-weight:700">${{n}}</a>`;
        }}).join(',');
        return `[${{links}}]`;
      }}
    );
    const bibEntries = [];
    citedIndices.forEach(n => {{
      const p = papers[n - 1];
      if (!p) return;
      const href = p.doi ? `https://doi.org/${{p.doi}}` : `https://pubmed.ncbi.nlm.nih.gov/${{p.pmid}}/`;
      bibEntries.push({{ n, p, href }});
    }});
    bibEntries.sort((a, b) => a.n - b.n);
    const bibHtml = bibEntries.map(entry =>
      `<div style="display:flex;gap:6px;margin-bottom:6px;font-size:11px;color:#94a3b8">
        <span style="color:#818cf8;font-weight:700;flex-shrink:0">${{entry.n}}.</span>
        <span><a href="${{entry.href}}" target="_blank" rel="noopener"
          style="color:#cbd5e1;text-decoration:underline dotted">${{escHtml(entry.p.title)}}</a>
          <span style="color:#475569"> · ${{entry.p.year}} · PMID ${{entry.p.pmid}}</span></span>
      </div>`
    ).join('');
    const summaryBlock =
      `<div style="line-height:1.7;color:#cbd5e1;font-size:13px;margin-bottom:${{bibHtml ? '12px' : '0'}}">${{linkedSummary}}</div>` +
      (bibHtml ? `<div style="border-top:1px solid #334155;padding-top:8px">${{bibHtml}}</div>` : '');
    html += section('Metabolic Summary',
      `<div style="border-left:3px solid #0ea5e9;padding-left:10px">${{summaryBlock}}</div>`);
  }} else if (data.description) {{
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
  const edge = cy.getElementById(edgeId);
  if (!edge.length) return;
  if (!searchInput.value.trim()) {{
    cy.elements().removeClass('faded highlighted');
    cy.elements().addClass('faded');
    edge.removeClass('faded').addClass('highlighted');
    edge.connectedNodes().removeClass('faded').addClass('highlighted');
  }}
  showEdgeSidebar(edge.data());
}}

// ── Click events ──────────────────────────────────────────────────────────

cy.on('tap', 'edge', function(evt) {{
  const el = evt.target;
  if (!searchInput.value.trim()) {{
    cy.elements().removeClass('faded highlighted');
    cy.elements().addClass('faded');
    el.removeClass('faded').addClass('highlighted');
    el.connectedNodes().removeClass('faded').addClass('highlighted');
  }}
  showEdgeSidebar(el.data());
}});

cy.on('tap', 'node', function(evt) {{
  const el = evt.target;
  if (!searchInput.value.trim()) {{
    cy.elements().removeClass('faded highlighted');
    cy.elements().addClass('faded');
    el.removeClass('faded').addClass('highlighted');
    el.connectedEdges().removeClass('faded').addClass('highlighted');
    el.connectedNodes().removeClass('faded');
  }}
  showNodeSidebar(el.data());
}});

cy.on('tap', function(evt) {{
  if (evt.target === cy) {{
    if (!searchInput.value.trim()) {{
      cy.elements().removeClass('faded highlighted');
    }}
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

// ── Hover effects ─────────────────────────────────────────────────────────

cy.on('mouseover', 'node, edge', function(evt) {{
  evt.target.addClass('hovered');
  document.getElementById('cy').style.cursor = 'pointer';
}});

cy.on('mouseout', 'node, edge', function(evt) {{
  evt.target.removeClass('hovered');
  document.getElementById('cy').style.cursor = 'default';
}});

// ── Search ────────────────────────────────────────────────────────────────

const searchInput = document.getElementById('search-input');
const searchResults = document.getElementById('search-results');
const searchCount = document.getElementById('search-count');

function buildSearchIndex() {{
  const index = [];
  cy.edges().forEach(e => {{
    const d = e.data();
    index.push({{ el: e, data: d, label: d.source + ' ↔ ' + d.target, type: 'edge' }});
  }});
  cy.nodes().forEach(n => {{
    const d = n.data();
    index.push({{ el: n, data: d, label: d.label, type: 'node' }});
  }});
  return index;
}}

function getSearchText(entry) {{
  if (entry.type === 'node') return '';   // key-player search only matches edges
  const d = entry.data;
  return [
    ...(d.key_players_hormones    || []),
    ...(d.key_players_metabolites || []),
    ...(d.key_players_proteins    || []),
    ...(d.key_players             || []),
  ].join(' ').toLowerCase();
}}

// Return the mention count for a query term on an edge (checks all categories).
// Handles partial matches: "bile" finds max count among "bile acid", "bile salt", etc.
function getTermCount(d, q) {{
  const allCounts = Object.assign(
    {{}},
    d.key_players_counts_hormones    || {{}},
    d.key_players_counts_metabolites || {{}},
    d.key_players_counts_proteins    || {{}}
  );
  let best = 0;
  for (const [term, n] of Object.entries(allCounts)) {{
    if (term.toLowerCase().includes(q)) best = Math.max(best, n);
  }}
  return best;
}}

const searchIndex = buildSearchIndex();

function runSearch() {{
  const q = searchInput.value.trim().toLowerCase();

  if (!q) {{
    cy.nodes().removeClass('search-hit faded');
    cy.edges().removeClass('faded search-hit');
    searchResults.style.display = 'none';
    searchCount.textContent = '';
    return;
  }}

  // Only edges match (key-player search).
  // If a threshold is set, the searched term must appear ≥ kpThreshold times on the edge.
  const edgeHits = searchIndex.filter(entry => {{
    if (entry.type !== 'edge') return false;
    if (!getSearchText(entry).includes(q)) return false;
    if (kpThreshold > 0 && getTermCount(entry.data, q) < kpThreshold) return false;
    return true;
  }});
  const hitEls = cy.collection(edgeHits.map(h => h.el));

  // Fade non-matching edges; nodes always stay visible
  cy.nodes().removeClass('search-hit faded');
  cy.edges().removeClass('search-hit').addClass('faded');
  hitEls.removeClass('faded').addClass('search-hit');

  searchCount.textContent = `${{edgeHits.length}} edge${{edgeHits.length !== 1 ? 's' : ''}}`;

  // Dropdown results list — show count for the searched term
  if (edgeHits.length === 0) {{
    searchResults.innerHTML = '<div class="search-result-item" style="color:#64748b">No matching edges</div>';
  }} else {{
    searchResults.innerHTML = edgeHits.slice(0, 20).map(h => {{
      const n = getTermCount(h.data, q);
      const badge = n ? `<span style="margin-left:auto;background:#1e3a5f;color:#7dd3fc;border-radius:8px;padding:1px 7px;font-size:11px;font-weight:700">${{n}}×</span>` : '';
      return `
      <div class="search-result-item" data-id="${{h.el.id()}}">
        <span>🔗</span>
        <span style="flex:1">${{escHtml(h.label)}}</span>
        ${{badge}}
      </div>`;
    }}).join('');
    searchResults.querySelectorAll('.search-result-item[data-id]').forEach(item => {{
      item.addEventListener('click', () => {{
        searchResults.style.display = 'none';
        showEdgeSidebar(cy.getElementById(item.dataset.id).data());
        cy.center(cy.getElementById(item.dataset.id));
      }});
    }});
  }}
  searchResults.style.display = 'block';
}}

searchInput.addEventListener('input', runSearch);

// Hide search dropdown on outside click
document.addEventListener('click', e => {{
  if (!document.getElementById('search-wrap').contains(e.target)) {{
    searchResults.style.display = 'none';
  }}
  if (!document.getElementById('kp-filter-wrap').contains(e.target)) {{
    document.getElementById('kp-filter-dropdown').classList.remove('open');
  }}
}});
searchInput.addEventListener('focus', () => {{
  if (searchInput.value.trim()) searchResults.style.display = 'block';
}});

// ── Info panel ───────────────────────────────────────────────────────────

function toggleInfoPanel() {{
  const panel = document.getElementById('info-panel');
  const btn   = document.getElementById('info-btn');
  panel.classList.toggle('open');
  btn.classList.toggle('active');
}}

function showInfoSection(name) {{
  document.querySelectorAll('.info-section').forEach(s => s.classList.remove('active'));
  document.querySelectorAll('.info-tab').forEach(t => t.classList.remove('active'));
  document.getElementById('info-section-' + name).classList.add('active');
  event.currentTarget.classList.add('active');
}}

// Fill dynamic stat: edges with ≥1 paper
(function() {{
  const count = cy.edges().filter(e => (e.data('pubmed_n') || 0) > 0).length;
  const total = cy.edges().length;
  document.getElementById('stat-val-edges').textContent = count + ' / ' + total;
}})();

// ── Key Players filter ────────────────────────────────────────────────────

function toggleKpDropdown(e) {{
  e.stopPropagation();
  document.getElementById('kp-filter-dropdown').classList.toggle('open');
}}

function applyKpThreshold(val) {{
  kpThreshold = Math.max(0, parseInt(val) || 0);
  const btn   = document.getElementById('kp-filter-btn');
  const badge = document.getElementById('kp-threshold-badge');
  if (kpThreshold > 0) {{
    badge.textContent = kpThreshold;
    badge.style.display = 'inline';
    btn.classList.add('active');
  }} else {{
    badge.style.display = 'none';
    btn.classList.remove('active');
  }}
  // Re-render sidebar if an edge is open
  if (currentEdgeData) showEdgeSidebar(currentEdgeData);
  // Re-run search with new threshold
  runSearch();
}}

function resetKpThreshold() {{
  document.getElementById('kp-threshold-input').value = '0';
  applyKpThreshold(0);
}}

document.getElementById('kp-threshold-input').addEventListener('input', function() {{
  applyKpThreshold(this.value);
}});

</script>
</body>
</html>"""

    with open(filename, "w", encoding="utf-8") as f:
        f.write(html)

    graph_type = "Directed" if directed else "Undirected"
    print(f"[✔] Dashboard saved to {filename} ({graph_type})")
