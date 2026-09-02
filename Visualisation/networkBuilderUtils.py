import base64
import gzip
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


def _html_attr_escape(text: str) -> str:
    """
    Minimal escaping for embedding text in a double-quoted HTML attribute.
    Not named using the `html` module to avoid colliding with the local
    variable `html` (the assembled page string) inside the exporter below.
    """
    return (text.replace("&", "&amp;")
                .replace('"', "&quot;")
                .replace("<", "&lt;")
                .replace(">", "&gt;"))


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
    filename: str | None = "network_dashboard.html",
    directed: bool | None = None,
    include_legend: bool = True,
    title: str = "Metabolic Reference Network",
    start_layout: str = "cose",
    info_panel_tabs: list | None = None,
    extra_topbar_buttons: list | None = None,
    total_possible_edges: int | None = None,
    comparison_toggle_label: str | None = None,
    comparison_lit_label: str | None = None,
    comparison_ref_label: str = "Reference Network",
    threshold_control: dict | None = None,
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

    extra_topbar_buttons (optional): list of {"icon", "label", "url"} dicts.
    Each renders as a top-bar button next to the info ("i") button; clicking
    it opens `url` in a full-screen overlay layer (iframe) above everything,
    including the graph — the graph itself is never resized or pushed.
    Omit (default) for the standard topbar with no behavior change.

    total_possible_edges (optional): denominator for the "Edges"/"Robust
    edges" stat in the info panel (id="stat-val-edges"). When given, the
    stat reads "<edges actually in the graph> / <total_possible_edges>" —
    e.g. robust edges out of every organ pair considered. Omit (default) to
    keep the original behavior: "<edges with >=1 paper> / <edges in the
    graph>", which is what a graph containing ALL candidate edges (not just
    ones that already passed a threshold) wants. Edges flagged
    only_reference=True on their merged_data are excluded from this stat's
    numerator regardless (they're not "in" the base network — see below).

    comparison_toggle_label (optional): renders a topbar toggle button with
    this label. Clicking it flips every edge between its normal color/type
    (edge's own `color` + `connection_type`) and a comparison view computed
    live from two per-edge merged_data flags: `bootstrap_mean` (float — this
    edge's literature-robustness score) and `is_ref_edge` (bool — whether
    this organ pair is a predefined reference-network edge, with
    `ref_n_papers` for its paper count). An edge with bootstrap_mean below
    the current threshold (see threshold_control) is hidden unless
    comparison mode is on AND is_ref_edge is True, in which case it's shown
    in amber as "only in the reference network". No separate page or iframe
    is involved — this recolors the same graph in place. Requires
    threshold_control to be meaningful; omit (default) to skip the toggle.
    comparison_lit_label/comparison_ref_label name the two networks being
    compared in the generated labels/descriptions (e.g. "Healthy
    Literature-Based Network" / "Reference Metabolic Network").

    threshold_control (optional): {"default": float, "elbow": float,
    "label": str} — renders a topbar slider that live-adjusts the
    literature-robustness threshold edges are filtered against (compared to
    each edge's merged_data `bootstrap_mean`), instead of baking a single
    threshold in at build time. "default" is the slider's initial value
    (typically the elbow-resolved MIN_BOOTSTRAP_MEAN used when the data was
    built); "elbow" is marked as a reset target; "label" captions the
    control (e.g. "Otsuka–Ochiai threshold"). The slider's min/max are
    derived at runtime from the actual bootstrap_mean values present. Edges
    with no bootstrap_mean set are always treated as robust (unaffected by
    the slider) — this only filters edges from callers that opt in. The
    "Robust edges" stat (id="stat-val-edges") tracks the live count.

    Returns the assembled HTML as a string. Pass filename=None to get the
    string back without writing a file (e.g. to embed it inline elsewhere
    via <iframe srcdoc>) — otherwise it's written to `filename` as before.
    """

    if directed is None:
        directed = graph.is_directed()

    elements = []
    edge_details: dict = {}

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

        # Support both the reference-network structure (key_players_merged / key_players_counts)
        # and the bootstrap-network structure (key_players_hormones / metabolites / proteins
        # and key_players_counts_hormones / metabolites / proteins stored directly).
        all_kp_categories = merged.get("key_players_merged", {})
        all_kp_counts     = merged.get("key_players_counts", {})
        kp_hormones    = merged.get("key_players_hormones",    all_kp_categories.get("hormones", []))
        kp_metabolites = merged.get("key_players_metabolites", all_kp_categories.get("metabolites", []))
        kp_proteins    = merged.get("key_players_proteins",    all_kp_categories.get("proteins", []))
        kp_cnt_hormones    = merged.get("key_players_counts_hormones",    all_kp_counts.get("hormones", {}))
        kp_cnt_metabolites = merged.get("key_players_counts_metabolites", all_kp_counts.get("metabolites", {}))
        kp_cnt_proteins    = merged.get("key_players_counts_proteins",    all_kp_counts.get("proteins", {}))
        kp_bootstrap       = merged.get("key_players_bootstrap", False)

        # Paper / pubmed data — reference network nests under "pubmed" key;
        # bootstrap networks store papers / n_papers_found / pubmed_query directly.
        pubmed = merged.get("pubmed", {})
        pubmed_papers = pubmed.get("papers", []) or merged.get("papers", [])
        pubmed_n      = pubmed.get("n_papers", 0) or merged.get("n_papers_found", 0)
        pubmed_query  = pubmed.get("query", "")   or merged.get("pubmed_query", "")

        connection_type = (
            merged.get("connection_type")
            or parsed["connection_type"]
        )
        connection_type_secondary = merged.get("connection_type_secondary", "")
        notes = merged.get("notes") or parsed["notes"]
        sources = merged.get("sources") or parsed["sources"]
        ai_description = merged.get("ai_description", "")

        edge_color = attrs.get("color", "#64748b")

        # Adaptive-threshold / comparison-toggle fields (see threshold_control
        # and comparison_toggle_label docstrings above) — bootstrap_mean is
        # None (edge always treated as robust, unaffected by the slider) for
        # callers that don't set it.
        bootstrap_mean = merged.get("bootstrap_mean", None)
        is_ref_edge    = bool(merged.get("is_ref_edge", False))
        ref_n_papers   = merged.get("ref_n_papers", 0)

        edge_id = f"{u}__{v}"
        elements.append({
            "data": {
                "id": edge_id,
                "source": str(u),
                "target": str(v),
                "color": edge_color,
                "connection_type":           connection_type,
                "connection_type_secondary": connection_type_secondary,
                "bootstrapMean": bootstrap_mean,
                "isRefEdge":     is_ref_edge,
                "refNPapers":    ref_n_papers,
                # Key players stay eager (not in edge_details below) because
                # the search bar indexes every edge's key players at page
                # load, before any edge has been clicked.
                "key_players": key_players,
                "key_players_metabolites":       kp_metabolites,
                "key_players_hormones":          kp_hormones,
                "key_players_proteins":          kp_proteins,
                "key_players_counts_metabolites": kp_cnt_metabolites,
                "key_players_counts_hormones":    kp_cnt_hormones,
                "key_players_counts_proteins":    kp_cnt_proteins,
                "key_players_bootstrap":          kp_bootstrap,
                "pubmed_n": pubmed_n,
            }
        })

        # Heavy, click-only fields (raw description, notes/sources, full
        # paper list with abstracts, the LLM summary and its PubMed query)
        # are kept OUT of the main elements payload — that's what Cytoscape
        # parses eagerly for every edge at page load. Instead they're
        # gzip-compressed into a single side blob (see edge_details_gz_b64
        # below) and lazily decompressed + attached only when the user
        # actually opens an edge's sidebar (see getEdgeDetails() /
        # showEdgeSidebar() in the generated JS). For networks with large
        # per-edge paper lists (bootstrap networks scanning every organ
        # pair) this is the dominant driver of file size.
        edge_details[edge_id] = {
            "description": raw_desc,
            "notes": notes,
            "sources": sources,
            "ai_description": ai_description,
            "pubmed_papers": pubmed_papers,
            "pubmed_query": pubmed_query,
        }

    elements_json = json.dumps(elements, ensure_ascii=False)
    edge_details_json = json.dumps(edge_details, ensure_ascii=False)
    edge_details_gz_b64 = base64.b64encode(
        gzip.compress(edge_details_json.encode("utf-8"), compresslevel=9)
    ).decode("ascii")
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

    total_possible_edges_js = (
        "null" if total_possible_edges is None else str(int(total_possible_edges))
    )

    # Optional extra top-bar buttons, each opening a full-screen overlay
    # layer (iframe) — used e.g. to surface a sibling HTML file (bootstrap
    # overview, comparison network) without resizing or hiding the graph.
    #
    # Each button is either:
    #   {"url": "..."}   an external/relative link — iframe navigates via src.
    #   {"html": "..."}  inline content, however large — stored in a
    #                     <script type="application/json"> element (immune to
    #                     the ~2MB hard URL-length limit that data: URIs hit
    #                     when assigned to iframe.src) and applied via
    #                     iframe.srcdoc, which has no such limit.
    extra_topbar_buttons = extra_topbar_buttons or []
    _layer_content_blocks = []
    _topbar_button_tags = []
    for _i, b in enumerate(extra_topbar_buttons):
        label_js = json.dumps(b["label"])
        if "html" in b:
            content_id = f"layer-content-{_i}"
            # </script within the payload would otherwise prematurely close
            # this tag — escaping the slash is valid JSON (`\/` decodes to
            # `/`) and neutralizes it for the outer HTML tokenizer.
            escaped = re.sub(r"</(script)", r"<\\/\1", json.dumps(b["html"]),
                             flags=re.IGNORECASE)
            _layer_content_blocks.append(
                f'<script type="application/json" id="{content_id}">{escaped}</script>'
            )
            onclick = f"openLayerOverlay(null, {label_js}, {json.dumps(content_id)})"
        else:
            onclick = f"openLayerOverlay({json.dumps(b['url'])}, {label_js}, null)"
        # onclick is single-quoted: json.dumps() always produces double-quoted
        # JS string literals, so the attribute delimiter must be the other kind.
        _topbar_button_tags.append(
            f'<button class="ctrl-btn" title="{_html_attr_escape(b["label"])}" '
            f"onclick='{onclick}'>{b.get('icon', '')} {b['label']}</button>"
        )
    topbar_extra_buttons_html = "\n  ".join(_topbar_button_tags)
    layer_content_html = "\n".join(_layer_content_blocks)

    # Adaptive-threshold slider + comparison toggle: both act on the SAME
    # graph in place (no separate page/iframe). The slider live-filters
    # edges by merged_data.bootstrap_mean vs. a threshold; the toggle
    # switches between normal styling and a "vs reference" view computed
    # from is_ref_edge/ref_n_papers. See threshold_control and
    # comparison_toggle_label in the docstring above.
    threshold_html = ""
    comparison_toggle_html = ""
    # Fallback when neither feature is requested: a one-time stat fill
    # matching the pre-adaptive-threshold behavior.
    topbar_dynamic_js = f"""
document.getElementById('stat-val-edges').textContent = (function() {{
  const totalPossible = {total_possible_edges_js};
  const edges = cy.edges();
  if (totalPossible !== null) return edges.length + ' / ' + totalPossible;
  const withPapers = edges.filter(e => (e.data('pubmed_n') || 0) > 0).length;
  return withPapers + ' / ' + edges.length;
}})();
"""
    if threshold_control or comparison_toggle_label:
        default_t  = threshold_control.get("default") if threshold_control else None
        elbow_t    = threshold_control.get("elbow", default_t) if threshold_control else None
        unit_label = threshold_control.get("label", "Threshold") if threshold_control else ""

        if threshold_control:
            threshold_html = f"""
<div class="ctrl-group" style="display:flex;align-items:center;gap:6px;margin-right:8px;padding-right:8px;border-right:1px solid #334155">
  <label for="threshold-slider" style="font-size:0.72rem;color:#94a3b8;white-space:nowrap">{_html_attr_escape(unit_label)}</label>
  <input type="range" id="threshold-slider" step="0.0001" style="width:110px;vertical-align:middle">
  <span id="threshold-value" style="font-size:0.72rem;color:#e2e8f0;min-width:56px;display:inline-block;font-variant-numeric:tabular-nums"></span>
  <button type="button" id="threshold-reset" class="ctrl-btn" title="Reset to elbow suggestion" style="padding:2px 7px;font-size:0.68rem">↺ elbow</button>
</div>"""

        if comparison_toggle_label:
            comparison_toggle_html = (
                f'<button id="comparison-toggle-btn" class="ctrl-btn" '
                f'title="{_html_attr_escape(comparison_toggle_label)}" '
                f'onclick="toggleComparisonMode()">🔀 {comparison_toggle_label}</button>'
            )

        default_t_js = "null" if default_t is None else repr(float(default_t))
        elbow_t_js   = "null" if elbow_t   is None else repr(float(elbow_t))

        topbar_dynamic_js = f"""
const LIT_LABEL = {json.dumps(comparison_lit_label or "Literature-Based Network")};
const REF_LABEL = {json.dumps(comparison_ref_label)};
const COLOR_SHARED = '#10b981', COLOR_ONLY_LIT = '#38bdf8', COLOR_ONLY_REF = '#f59e0b';
let comparisonMode = false;
let currentThreshold = {default_t_js};

// Capture each edge's own styling once, before any threshold/comparison
// filtering runs, so normal mode can always restore it.
cy.edges().forEach(e => e.data('baseColor', e.data('color')));

function edgeIsRobust(d) {{
  return (d.bootstrapMean === null || d.bootstrapMean === undefined)
    ? true : (d.bootstrapMean >= currentThreshold);
}}

function applyThresholdAndComparison() {{
  cy.edges().forEach(e => {{
    const d = e.data();
    const robust = edgeIsRobust(d);
    const isRef  = !!d.isRefEdge;
    let visible, color, connType, desc = null;

    if (!comparisonMode) {{
      visible  = robust;
      color    = d.baseColor;
      connType = d.connection_type;
    }} else if (robust && isRef) {{
      visible = true; color = COLOR_SHARED; connType = 'SHARED — in both networks';
      desc = `Shared with the ${{REF_LABEL}} (reference papers: ${{d.refNPapers}}).`;
    }} else if (robust && !isRef) {{
      visible = true; color = COLOR_ONLY_LIT; connType = `ONLY in ${{LIT_LABEL}}`;
      desc = `Robust in the ${{LIT_LABEL}} bootstrap, but not a predefined edge in the ${{REF_LABEL}}.`;
    }} else if (!robust && isRef) {{
      visible = true; color = COLOR_ONLY_REF; connType = `ONLY in ${{REF_LABEL}}`;
      desc = `Predefined reference edge, below the ${{LIT_LABEL}} threshold (${{d.refNPapers}} reference papers).`;
    }} else {{
      visible = false; color = d.baseColor; connType = d.connection_type;
    }}

    e.style('display', visible ? 'element' : 'none');
    e.data('color', color);
    e.data('connection_type', connType);
    if (desc !== null) e.data('ai_description', desc);
  }});
  updateEdgeStat();
  if (typeof currentEdgeData !== 'undefined' && currentEdgeData) {{
    const refreshed = cy.getElementById(currentEdgeData.id);
    if (refreshed && refreshed.length) showEdgeSidebar(refreshed.data());
  }}
}}

function updateEdgeStat() {{
  const el = document.getElementById('stat-val-edges');
  if (!el) return;
  const totalPossible = {total_possible_edges_js};
  const robustEdges = cy.edges().filter(e => edgeIsRobust(e.data()));
  if (totalPossible !== null) {{
    el.textContent = robustEdges.length + ' / ' + totalPossible;
  }} else {{
    const withPapers = robustEdges.filter(e => (e.data('pubmed_n') || 0) > 0).length;
    el.textContent = withPapers + ' / ' + robustEdges.length;
  }}
}}

function toggleComparisonMode() {{
  comparisonMode = !comparisonMode;
  const btn = document.getElementById('comparison-toggle-btn');
  if (btn) btn.classList.toggle('active', comparisonMode);
  applyThresholdAndComparison();
}}

(function() {{
  const slider = document.getElementById('threshold-slider');
  if (slider) {{
    const means = cy.edges().map(e => e.data('bootstrapMean'))
      .filter(m => m !== null && m !== undefined);
    const minMean = means.length ? Math.min(...means) : 0;
    const maxMean = means.length ? Math.max(...means) : 1;
    const valueEl = document.getElementById('threshold-value');
    slider.min = minMean; slider.max = maxMean; slider.value = currentThreshold;
    valueEl.textContent = currentThreshold.toFixed(5);
    slider.addEventListener('input', () => {{
      currentThreshold = parseFloat(slider.value);
      valueEl.textContent = currentThreshold.toFixed(5);
      applyThresholdAndComparison();
    }});
    const resetBtn = document.getElementById('threshold-reset');
    if (resetBtn) resetBtn.addEventListener('click', () => {{
      const elbow = {elbow_t_js};
      if (elbow === null) return;
      currentThreshold = elbow; slider.value = elbow;
      valueEl.textContent = elbow.toFixed(5);
      applyThresholdAndComparison();
    }});
  }}
}})();

applyThresholdAndComparison();
"""

    layer_overlay_css = ""
    layer_overlay_markup = ""
    layer_overlay_js = ""
    if extra_topbar_buttons:
        layer_overlay_css = """
  .layer-overlay {
    display: none;
    position: fixed;
    inset: 0;
    z-index: 1000;
    background: #0f172a;
    flex-direction: column;
  }
  .layer-overlay.active { display: flex; }
  .layer-overlay-header {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 10px 20px;
    background: #1e293b;
    border-bottom: 1px solid #334155;
    flex-shrink: 0;
  }
  .layer-overlay-header h2 {
    font-size: 15px;
    font-weight: 700;
    color: #f1f5f9;
  }
  .layer-overlay-close {
    background: none;
    border: none;
    color: #94a3b8;
    font-size: 20px;
    cursor: pointer;
    line-height: 1;
    padding: 4px 8px;
  }
  .layer-overlay-close:hover { color: #f1f5f9; }
  .layer-overlay-iframe {
    flex: 1;
    border: none;
    width: 100%;
  }
"""
        layer_overlay_markup = """
<div id="layer-overlay" class="layer-overlay">
  <div class="layer-overlay-header">
    <h2 id="layer-overlay-title"></h2>
    <button class="layer-overlay-close" onclick="closeLayerOverlay()" title="Close">✕</button>
  </div>
  <iframe id="layer-overlay-iframe" class="layer-overlay-iframe" src="about:blank" title="Overlay content"></iframe>
</div>"""
        layer_overlay_js = """
function openLayerOverlay(url, title, contentId) {
  document.getElementById('layer-overlay-title').textContent = title;
  const frame = document.getElementById('layer-overlay-iframe');
  if (contentId) {
    // Inline content, however large: read the JSON-encoded payload and
    // assign via srcdoc, which (unlike iframe.src on a data: URI) has no
    // practical URL-length ceiling.
    const dataEl = document.getElementById(contentId);
    let payload = '<p style="font-family:sans-serif;color:#f87171;padding:20px">Content not found.</p>';
    if (dataEl) {
      try { payload = JSON.parse(dataEl.textContent); }
      catch (e) { payload = '<p style="font-family:sans-serif;color:#f87171;padding:20px">Failed to load content: ' + e.message + '</p>'; }
    }
    frame.removeAttribute('src');
    frame.srcdoc = payload;
  } else {
    frame.removeAttribute('srcdoc');
    frame.src = url;
  }
  document.getElementById('layer-overlay').classList.add('active');
}
function closeLayerOverlay() {
  document.getElementById('layer-overlay').classList.remove('active');
  const frame = document.getElementById('layer-overlay-iframe');
  frame.removeAttribute('srcdoc');
  frame.src = 'about:blank';
}
"""

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
    width: 760px;
    min-width: 500px;
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
    flex-wrap: wrap;
    gap: 2px 4px;
    padding: 10px 16px 0;
    border-bottom: 1px solid #334155;
    flex-shrink: 0;
  }}

  .info-tab {{
    background: none;
    border: none;
    border-bottom: 2px solid transparent;
    color: #64748b;
    font-size: 10.5px;
    font-weight: 600;
    text-transform: uppercase;
    letter-spacing: 0.06em;
    padding: 6px 8px 8px;
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
    transition: width 0.25s ease, min-width 0.25s ease;
  }}
  #sidebar.closed {{
    width: 0;
    min-width: 0;
    border-left: none;
  }}

  #sidebar-header {{
    padding: 14px 16px 10px;
    border-bottom: 1px solid #334155;
    flex-shrink: 0;
    display: flex;
    flex-direction: column;
    gap: 2px;
  }}
  #sidebar-header-row {{
    display: flex;
    align-items: flex-start;
    justify-content: space-between;
    gap: 8px;
  }}

  #sidebar-close {{
    background: none;
    border: none;
    color: #64748b;
    font-size: 18px;
    cursor: pointer;
    line-height: 1;
    padding: 2px 4px;
    flex-shrink: 0;
  }}
  #sidebar-close:hover {{ color: #f1f5f9; }}

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

  #sidebar-tabs {{
    display: flex;
    border-bottom: 1px solid #334155;
    flex-shrink: 0;
    padding: 0 8px;
    gap: 2px;
    background: #1e293b;
  }}

  .sidebar-tab {{
    padding: 8px 12px;
    font-size: 12px;
    font-weight: 600;
    color: #64748b;
    cursor: pointer;
    border-bottom: 2px solid transparent;
    white-space: nowrap;
    transition: color 0.15s, border-color 0.15s;
    background: none;
    border-top: none;
    border-left: none;
    border-right: none;
  }}

  .sidebar-tab:hover {{ color: #cbd5e1; }}
  .sidebar-tab.active {{ color: #818cf8; border-bottom-color: #818cf8; }}

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
{layer_overlay_css}
</style>
</head>
<body>

<!-- Top bar -->
<div id="topbar">
  <button id="info-btn" onclick="toggleInfoPanel()" title="How it works">i</button>
  {threshold_html}
  {comparison_toggle_html}
  {topbar_extra_buttons_html}
  {layer_content_html}

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
      <div id="sidebar-header-row">
        <div id="sidebar-title">Select a node or edge</div>
        <button id="sidebar-close" onclick="closeSidebar()" title="Close">✕</button>
      </div>
      <div id="sidebar-name">—</div>
    </div>
    <div id="sidebar-tabs"></div>
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
{layer_overlay_markup}
<script>
const elements = {elements_json};

// Heavy per-edge fields (full paper list, notes/sources, LLM summary) live
// here gzip-compressed instead of in `elements` above, and are decompressed
// once, lazily, on first edge click — see getEdgeDetails() below.
const EDGE_DETAILS_GZ_B64 = "{edge_details_gz_b64}";
let _edgeDetailsPromise = null;
async function getEdgeDetails() {{
  if (_edgeDetailsPromise) return _edgeDetailsPromise;
  _edgeDetailsPromise = (async () => {{
    if (typeof DecompressionStream === 'undefined') return {{}};
    try {{
      const binary = atob(EDGE_DETAILS_GZ_B64);
      const bytes = new Uint8Array(binary.length);
      for (let i = 0; i < binary.length; i++) bytes[i] = binary.charCodeAt(i);
      const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
      const text = await new Response(stream).text();
      return JSON.parse(text);
    }} catch (err) {{
      console.warn('Edge detail decompression failed:', err);
      return {{}};
    }}
  }})();
  return _edgeDetailsPromise;
}}

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

function chipsHtml(items, cls, counts, unit) {{
  unit = unit || '×';
  const isBootstrap = unit === '%';
  if (!items || items.length === 0) return '<span style="color:#475569;font-size:12px">None identified</span>';
  const visible = kpThreshold > 0
    ? items.filter(t => (counts && counts[t] || 0) >= kpThreshold)
    : items;
  if (visible.length === 0) return `<span style="color:#475569;font-size:12px">None with ≥${{kpThreshold}} ${{isBootstrap ? '% bootstrap freq' : 'mentions'}}</span>`;
  return '<div class="chips-wrap">' +
    visible.map(t => {{
      const n = counts && counts[t];
      const label = isBootstrap ? (n != null ? n + '%' : '') : (n ? n + '×' : '');
      const title = isBootstrap
        ? (n != null ? `appeared in ${{n}}% of bootstrap iterations` : '')
        : (n ? n + ' mentions' : '');
      const badge = label ? `<span style="margin-left:4px;background:rgba(0,0,0,0.25);border-radius:8px;padding:0 5px;font-size:10px;font-weight:700">${{label}}</span>` : '';
      return `<span class="chip ${{cls}}" title="${{title}}">${{escHtml(t)}}${{badge}}</span>`;
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

// ── Sidebar tab helpers ───────────────────────────────────────────────────

let _sidebarTabs = {{}};   // {{ id: htmlString }}
let _activeSidebarTab = '';

function _renderSidebarTabs(tabs, activeId) {{
  _sidebarTabs = tabs;
  _activeSidebarTab = activeId;
  const tabBar = document.getElementById('sidebar-tabs');
  tabBar.innerHTML = Object.keys(tabs).map(id =>
    `<button class="sidebar-tab${{id === activeId ? ' active' : ''}}"
       onclick="switchSidebarTab('${{id}}')">${{id}}</button>`
  ).join('');
  document.getElementById('sidebar-body').innerHTML = tabs[activeId] ||
    '<p style="color:#475569;font-size:13px">No information available.</p>';
}}

function switchSidebarTab(id) {{
  _activeSidebarTab = id;
  document.querySelectorAll('.sidebar-tab').forEach(b =>
    b.classList.toggle('active', b.textContent === id));
  document.getElementById('sidebar-body').innerHTML = _sidebarTabs[id] ||
    '<p style="color:#475569;font-size:13px">No information available.</p>';
}}

// ── Linked summary helper (citations → superscript links + bibliography) ──

function _linkedSummaryHtml(text, papers, accentColor) {{
  if (!text) return '';

  function parseNums(inner) {{
    const nums = [];
    inner.split(/[,\\s]+/).forEach(part => {{
      const dash = part.match(/^(\\d+)-(\\d+)$/);
      if (dash) {{
        for (let i = parseInt(dash[1]); i <= parseInt(dash[2]); i++) nums.push(i);
      }} else {{
        const n = parseInt(part);
        if (!isNaN(n)) nums.push(n);
      }}
    }});
    return nums;
  }}

  const citePattern = /\\[(\\d[\\d,\\s-]*)\\]/g;

  // First pass: collect every cited raw paper index, in first-appearance
  // order, so the citation numbers shown in the text and bibliography
  // always run 1..n with no gaps — regardless of which original positions
  // (e.g. 1, 4, 5 out of a larger given set) the model actually cited.
  const rawOrder = [];
  const seenRaw = new Set();
  let m;
  while ((m = citePattern.exec(text)) !== null) {{
    parseNums(m[1]).forEach(n => {{
      if (!seenRaw.has(n)) {{ seenRaw.add(n); rawOrder.push(n); }}
    }});
  }}
  const displayNum = {{}};
  rawOrder.forEach((n, i) => {{ displayNum[n] = i + 1; }});

  citePattern.lastIndex = 0;
  const linked = escHtml(text).replace(citePattern, (match, inner) => {{
    const links = parseNums(inner).map(n => {{
      const p = papers[n - 1];
      const href = p ? (p.doi ? `https://doi.org/${{p.doi}}` : `https://pubmed.ncbi.nlm.nih.gov/${{p.pmid}}/`) : '#';
      const shown = displayNum[n] || n;
      return `<a href="${{href}}" target="_blank" rel="noopener"
        style="color:#818cf8;text-decoration:none;font-weight:700">${{shown}}</a>`;
    }}).join(',');
    return `[${{links}}]`;
  }});

  const bibEntries = [];
  rawOrder.forEach(n => {{
    const p = papers[n - 1];
    if (!p) return;
    const href = p.doi ? `https://doi.org/${{p.doi}}` : `https://pubmed.ncbi.nlm.nih.gov/${{p.pmid}}/`;
    bibEntries.push({{ n: displayNum[n], p, href }});
  }});
  bibEntries.sort((a, b) => a.n - b.n);
  const bibHtml = bibEntries.length ? bibEntries.map(e =>
    `<div style="display:flex;gap:6px;margin-bottom:6px;font-size:11px;color:#94a3b8">
      <span style="color:#818cf8;font-weight:700;flex-shrink:0">${{e.n}}.</span>
      <span><a href="${{e.href}}" target="_blank" rel="noopener"
        style="color:#cbd5e1;text-decoration:underline dotted">${{escHtml(e.p.title)}}</a>
        <span style="color:#475569"> · ${{e.p.year}} · PMID ${{e.p.pmid}}</span></span>
    </div>`
  ).join('') : '';
  const summaryBlock =
    `<div style="line-height:1.7;color:#cbd5e1;font-size:13px;margin-bottom:${{bibHtml ? '12px' : '0'}}">${{linked}}</div>` +
    (bibHtml ? `<div style="border-top:1px solid #334155;padding-top:8px">${{bibHtml}}</div>` : '');
  return `<div style="border-left:3px solid ${{accentColor}};padding-left:10px">${{summaryBlock}}</div>`;
}}

// ── Sidebar rendering ─────────────────────────────────────────────────────

function renderEdgeSidebarContent(data) {{
  document.getElementById('sidebar-title').textContent = 'Connection';
  document.getElementById('sidebar-name').textContent =
    data.source + ' ↔ ' + data.target;

  // ── Tab: Summary ──────────────────────────────────────────────────────
  let summaryHtml = '';

  if (data.connection_type) {{
    const secAttr = data.connection_type_secondary
      ? ` title="Alternative: ${{escHtml(data.connection_type_secondary)}}"`
      : '';
    summaryHtml += section('Type',
      `<span class="connection-type-badge"${{secAttr}} style="cursor:${{data.connection_type_secondary?'help':'default'}}">${{escHtml(data.connection_type)}}</span>` +
      (data.connection_type_secondary
        ? `<div style="margin-top:4px;font-size:11px;color:#64748b">Alternative: ${{escHtml(data.connection_type_secondary)}}</div>`
        : ''));
  }}

  if (data.ai_description) {{
    summaryHtml += section('Literature Summary',
      _linkedSummaryHtml(data.ai_description, data.pubmed_papers || [], '#6366f1'));
  }} else if (!data._detailsLoaded && data.pubmed_n > 0) {{
    summaryHtml += section('PubMed Literature',
      `<div style="font-size:12px;color:#64748b">Loading details…</div>`);
  }} else if (data.pubmed_n > 0) {{
    summaryHtml += section('PubMed Literature',
      `<div style="font-size:12px;color:#64748b">Found ${{data.pubmed_n}} papers — run LLM step to generate summary.</div>`);
  }}

  // ── Tab: Key Players ──────────────────────────────────────────────────
  const hormones    = data.key_players_hormones    || [];
  const metabolites = data.key_players_metabolites || [];
  const proteins    = data.key_players_proteins    || [];
  const rawKP       = data.key_players             || [];
  const cntH = data.key_players_counts_hormones    || {{}};
  const cntM = data.key_players_counts_metabolites || {{}};
  const cntP = data.key_players_counts_proteins    || {{}};
  const kpUnit = data.key_players_bootstrap ? '%' : '×';
  const kpSectionTitle = data.key_players_bootstrap ? 'Key Players (bootstrap frequency)' : 'Key Players';

  let kpHtml = '';
  if (hormones.length || metabolites.length || proteins.length) {{
    const kpCountRow =
      '<div style="display:flex;gap:12px;margin-bottom:8px">' +
      (hormones.length    ? `<span style="font-size:11px;color:#64748b"><span style="color:#fdba74;font-weight:700">${{hormones.length}}</span> hormone${{hormones.length!==1?'s':''}}</span>` : '') +
      (metabolites.length ? `<span style="font-size:11px;color:#64748b"><span style="color:#5eead4;font-weight:700">${{metabolites.length}}</span> metabolite${{metabolites.length!==1?'s':''}}</span>` : '') +
      (proteins.length    ? `<span style="font-size:11px;color:#64748b"><span style="color:#a5b4fc;font-weight:700">${{proteins.length}}</span> protein${{proteins.length!==1?'s':''}}</span>` : '') +
      '</div>';
    kpHtml += section(kpSectionTitle,
      kpCountRow +
      '<div style="display:flex;flex-direction:column;gap:8px">' +
      (hormones.length    ? '<div><span style="font-size:11px;color:#64748b;margin-bottom:3px;display:block">Hormones</span>'               + chipsHtml(hormones,    'chip-hormone',    cntH, kpUnit) + '</div>' : '') +
      (metabolites.length ? '<div><span style="font-size:11px;color:#64748b;margin-bottom:3px;display:block">Metabolites</span>'            + chipsHtml(metabolites, 'chip-metabolite', cntM, kpUnit) + '</div>' : '') +
      (proteins.length    ? '<div><span style="font-size:11px;color:#64748b;margin-bottom:3px;display:block">Proteins / Transporters</span>' + chipsHtml(proteins,    'chip-protein',    cntP, kpUnit) + '</div>' : '') +
      '</div>'
    );
  }} else if (rawKP.length) {{
    kpHtml += section('Key Players',
      `<div style="font-size:11px;color:#64748b;margin-bottom:6px"><span style="color:#e2e8f0;font-weight:700">${{rawKP.length}}</span> key player${{rawKP.length!==1?'s':''}}</div>` +
      chipsHtml(rawKP, 'chip-general', {{}}, kpUnit));
  }} else {{
    kpHtml = '<p style="color:#475569;font-size:13px">No key players found.</p>';
  }}

  const tabs = {{}};
  tabs['Summary']     = summaryHtml || '<p style="color:#475569;font-size:13px">No summary available.</p>';
  tabs['Key Players'] = kpHtml;
  _renderSidebarTabs(tabs, 'Summary');
}}

function showEdgeSidebar(data) {{
  currentEdgeData = data;
  document.getElementById('sidebar').classList.remove('closed');
  renderEdgeSidebarContent(data);

  // Heavy fields (papers, LLM summary, notes/sources) aren't in `data` yet —
  // see edge_details_gz_b64 / getEdgeDetails() above. Fetch (or reuse the
  // already-decompressed cache) and re-render in place once attached.
  if (!data._detailsLoaded) {{
    getEdgeDetails().then(details => {{
      const heavy = details[data.id];
      const e = cy.getElementById(data.id);
      if (!e.length) return;
      if (heavy) Object.keys(heavy).forEach(k => e.data(k, heavy[k]));
      e.data('_detailsLoaded', true);
      if (currentEdgeData && currentEdgeData.id === data.id) {{
        currentEdgeData = e.data();
        renderEdgeSidebarContent(currentEdgeData);
      }}
    }});
  }}
}}

function showNodeSidebar(data) {{
  currentEdgeData = null;
  document.getElementById('sidebar').classList.remove('closed');
  document.getElementById('sidebar-title').textContent = 'Organ';
  document.getElementById('sidebar-name').textContent = data.label;

  // ── Tab: Summary (LLM description) ───────────────────────────────────
  let summaryHtml = '';
  if (data.llm_description) {{
    summaryHtml += section('Metabolic Summary',
      _linkedSummaryHtml(data.llm_description, data.llm_papers || [], '#0ea5e9'));
  }}

  // ── Tab: Description (curated text) ──────────────────────────────────
  let descHtml = '';
  if (data.description) {{
    descHtml += section('Description',
      `<p class="notes-text">${{escHtml(data.description)}}</p>`);
  }}

  // ── Tab: Connections list ─────────────────────────────────────────────
  let connHtml = '';
  const node = cy.getElementById(data.id);
  const edges = node.connectedEdges();
  if (edges.length) {{
    connHtml += section(`Connections (${{edges.length}})`,
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

  const tabs = {{}};
  if (summaryHtml) tabs['Summary']     = summaryHtml;
  if (descHtml)    tabs['Description'] = descHtml;
  if (connHtml)    tabs['Connections'] = connHtml;

  if (Object.keys(tabs).length === 0) {{
    document.getElementById('sidebar-tabs').innerHTML = '';
    document.getElementById('sidebar-body').innerHTML =
      '<p style="color:#475569;font-size:13px">No description available.</p>';
  }} else {{
    _renderSidebarTabs(tabs, Object.keys(tabs)[0]);
  }}
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
    document.getElementById('sidebar-tabs').innerHTML = '';
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
{layer_overlay_js}

function closeSidebar() {{
  document.getElementById('sidebar').classList.add('closed');
}}

function showInfoSection(name) {{
  document.querySelectorAll('.info-section').forEach(s => s.classList.remove('active'));
  document.querySelectorAll('.info-tab').forEach(t => t.classList.remove('active'));
  document.getElementById('info-section-' + name).classList.add('active');
  event.currentTarget.classList.add('active');
}}

// Fill dynamic "Robust edges" stat and wire up the threshold slider /
// comparison toggle, if either was requested (see threshold_control and
// comparison_toggle_label in the docstring above). When neither is used,
// this falls back to a plain one-time stat fill matching the old behavior.
{topbar_dynamic_js}

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

    if filename:
        with open(filename, "w", encoding="utf-8") as f:
            f.write(html)
        graph_type = "Directed" if directed else "Undirected"
        print(f"[✔] Dashboard saved to {filename} ({graph_type})")

    return html
