"""
Build the interactive multi-layer literature reference network HTML.

Each layer toggle shows/hides edges of that type independently.
Edge colors are fixed (determined by the layer with most papers for that edge)
and do not change when layers are toggled.

Usage (standalone):
    uv run python Literature_Reference_Network/lit_ref_viz.py
    uv run python Literature_Reference_Network/lit_ref_viz.py --min-papers 10
"""

import json
import sys
import argparse
from pathlib import Path

HERE = Path(__file__).parent.parent
sys.path.insert(0, str(HERE))

DEFAULT_SEARCH_RESULTS = HERE / "metabolic_data" / "lit_ref_results.json"
DEFAULT_OUTPUT_HTML    = HERE / "metabolic_data" / "literature_reference_network.html"
DEFAULT_MIN_PAPERS     = 5

NODE_COLORS: dict[str, str] = {
    "Liver":          "#16a34a",
    "Pancreas":       "#dc2626",
    "WAT":            "#ca8a04",
    "Muscle":         "#2563eb",
    "Brain":          "#7c3aed",
    "Heart":          "#db2777",
    "Kidney":         "#0891b2",
    "Adrenal Glands": "#c2410c",
    "Thyroid":        "#059669",
    "Small Intestine":"#d97706",
    "Colon":          "#7c2d12",
    "Lung":           "#1d4ed8",
    "Bone Marrow":    "#6b21a8",
    "Spleen":         "#065f46",
}

LAYER_KP_LABEL: dict[str, str] = {
    "metabolic":  "Key Metabolites",
    "hormonal":   "Key Hormones",
    "immune":     "Key Immune Molecules",
    "neural":     "Key Neural Molecules",
    "mechanical": "Key Mechanical Signals",
    "undefined":  "Key Molecules",
}

LAYER_CHIP_STYLE: dict[str, str] = {
    "metabolic":  "background:#042f2e;color:#5eead4;border:1px solid #0f766e",
    "hormonal":   "background:#431407;color:#fdba74;border:1px solid #92400e",
    "immune":     "background:#450a0a;color:#fca5a5;border:1px solid #991b1b",
    "neural":     "background:#2e1065;color:#d8b4fe;border:1px solid #6d28d9",
    "mechanical": "background:#082f49;color:#7dd3fc;border:1px solid #0369a1",
    "undefined":  "background:#1e293b;color:#cbd5e1;border:1px solid #475569",
}


def build_lit_ref_viz(
    search_results: dict,
    node_metadata: dict,
    output_html: "str | Path" = DEFAULT_OUTPUT_HTML,
    min_papers: int = DEFAULT_MIN_PAPERS,
    title: str = "Literature Reference Network",
) -> None:
    from Literature_Reference_Network.lit_ref_search import (
        LAYER_COLORS, LAYER_LABELS, ALL_LAYERS,
        extract_layer_key_players,
    )

    output_html = Path(output_html)
    output_html.parent.mkdir(parents=True, exist_ok=True)

    organs = sorted(node_metadata.keys())

    # ── Detect cache format ───────────────────────────────────────────────────
    old_format = any(k.count("|") == 2 for k in search_results)

    # ── Build elements ────────────────────────────────────────────────────────
    elements = []
    for organ in organs:
        desc = node_metadata.get(organ, "")
        elements.append({"data": {
            "id":          organ,
            "label":       organ,
            "description": desc if isinstance(desc, str) else str(desc),
            "color":       NODE_COLORS.get(organ, "#64748b"),
        }})

    if old_format:
        # ── Legacy per-layer cache format ─────────────────────────────────────
        pair_layers: dict[tuple, dict] = {}
        for key, data in search_results.items():
            if key.count("|") != 2:
                continue
            o1, o2, layer = data["organ1"], data["organ2"], data["layer"]
            pair = (min(o1, o2), max(o1, o2))
            pair_layers.setdefault(pair, {})[layer] = {
                "n_papers": data["n_papers"],
                "papers":   data["papers"][:5],
            }

        for (o1, o2), layers in pair_layers.items():
            qualifying = {
                lname: ldata for lname, ldata in layers.items()
                if ldata["n_papers"] >= min_papers
            }
            if not qualifying:
                continue
            layer_fields: dict = {}
            for lname, ldata in qualifying.items():
                kp, kp_cnt = extract_layer_key_players(ldata["papers"], lname)
                layer_fields[f"lyr_{lname}_n"]         = ldata["n_papers"]
                layer_fields[f"lyr_{lname}_papers"]    = ldata["papers"]
                layer_fields[f"lyr_{lname}_kp"]        = kp
                layer_fields[f"lyr_{lname}_kp_counts"] = kp_cnt
            elements.append({"data": {
                "id": f"{o1}--{o2}", "source": o1, "target": o2,
                "color": "#60a5fa", "layers": list(qualifying.keys()),
                **layer_fields,
            }})

    else:
        # ── New co-occurrence cache format ────────────────────────────────────
        for key, data in search_results.items():
            if key.count("|") != 1:
                continue
            o1, o2 = data["organ1"], data["organ2"]
            if data.get("n_papers_cooccur", 0) < min_papers:
                continue

            layer_papers_raw = data.get("layer_papers", {})
            qualifying_layers = [
                lyr for lyr in ALL_LAYERS
                if len(layer_papers_raw.get(lyr, [])) >= min_papers
            ]
            if not qualifying_layers:
                continue

            layer_fields = {}
            for lname in qualifying_layers:
                papers      = layer_papers_raw.get(lname, [])[:5]
                kp, kp_cnt  = extract_layer_key_players(papers, lname)
                layer_fields[f"lyr_{lname}_n"]         = data["layer_counts"].get(lname, 0)
                layer_fields[f"lyr_{lname}_papers"]    = papers
                layer_fields[f"lyr_{lname}_kp"]        = kp
                layer_fields[f"lyr_{lname}_kp_counts"] = kp_cnt

            elements.append({"data": {
                "id":     f"{o1}--{o2}",
                "source": o1,
                "target": o2,
                "color":  "#60a5fa",
                "layers": qualifying_layers,
                "n_cooccur": data.get("n_papers_cooccur", 0),
                **layer_fields,
            }})

    elements_json   = json.dumps(elements, ensure_ascii=False)
    layer_meta_json = json.dumps({
        name: {"color": LAYER_COLORS[name], "label": LAYER_LABELS[name]}
        for name in ALL_LAYERS
    })
    layer_names_json  = json.dumps(ALL_LAYERS)
    layer_chip_json   = json.dumps(LAYER_CHIP_STYLE)
    layer_kp_lbl_json = json.dumps(LAYER_KP_LABEL)

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
    background: #0f172a; color: #e2e8f0;
    height: 100vh; display: flex; flex-direction: column; overflow: hidden;
  }}

  /* ── Top bar ── */
  #topbar {{
    display: flex; align-items: center; gap: 10px;
    padding: 10px 20px; background: #1e293b;
    border-bottom: 1px solid #334155; flex-shrink: 0; z-index: 10;
    flex-wrap: wrap;
  }}
  #topbar h1 {{ font-size: 15px; font-weight: 700; color: #f1f5f9; white-space: nowrap; }}

  #info-btn {{
    width: 30px; height: 30px; border-radius: 50%;
    background: #1e293b; border: 1px solid #334155;
    color: #94a3b8; font-size: 14px; font-weight: 700;
    cursor: pointer; flex-shrink: 0;
    display: flex; align-items: center; justify-content: center;
    transition: border-color 0.2s, color 0.2s;
  }}
  #info-btn:hover {{ border-color: #6366f1; color: #a5b4fc; }}
  #info-btn.active {{ background: #1e1b4b; border-color: #6366f1; color: #a5b4fc; }}

  .sep {{ width: 1px; height: 22px; background: #334155; flex-shrink: 0; }}

  /* ── Layer toggles ── */
  #layer-toggles {{ display: flex; gap: 6px; align-items: center; flex-wrap: wrap; }}
  .layer-toggle {{
    display: flex; align-items: center; gap: 5px;
    padding: 5px 12px; border-radius: 999px;
    font-size: 12px; font-weight: 600;
    cursor: pointer; border: 2px solid;
    transition: opacity 0.2s, background 0.2s;
    user-select: none; white-space: nowrap;
  }}
  .layer-toggle .dot {{ width: 8px; height: 8px; border-radius: 50%; flex-shrink: 0; }}
  .layer-toggle.off {{
    opacity: 0.35; border-color: #334155 !important;
    background: #1e293b !important; color: #64748b !important;
  }}
  .layer-toggle.off .dot {{ background: #334155 !important; }}

  /* ── Search ── */
  #search-wrap {{ position: relative; flex: 1; max-width: 300px; min-width: 140px; }}
  #search-icon {{ position: absolute; left: 9px; top: 50%; transform: translateY(-50%); color: #94a3b8; pointer-events: none; font-size: 13px; }}
  #search-input {{
    width: 100%; padding: 6px 10px 6px 28px;
    background: #0f172a; border: 1px solid #334155;
    border-radius: 8px; color: #e2e8f0; font-size: 13px; outline: none;
  }}
  #search-input:focus {{ border-color: #6366f1; }}
  #search-input::placeholder {{ color: #64748b; }}
  #search-results {{
    display: none; position: absolute; top: calc(100% + 4px); left: 0; right: 0;
    background: #1e293b; border: 1px solid #334155; border-radius: 8px;
    z-index: 300; max-height: 240px; overflow-y: auto;
    box-shadow: 0 8px 24px rgba(0,0,0,0.5);
  }}
  .sr-item {{
    padding: 8px 14px; cursor: pointer; font-size: 13px; color: #cbd5e1;
    display: flex; align-items: center; gap: 8px; transition: background 0.1s;
  }}
  .sr-item:hover {{ background: #334155; }}
  #search-count {{ font-size: 12px; color: #64748b; white-space: nowrap; }}

  /* ── Layout btns ── */
  .ctrl-btn {{
    padding: 5px 10px; background: #0f172a; border: 1px solid #334155;
    border-radius: 6px; color: #94a3b8; font-size: 12px; cursor: pointer;
    white-space: nowrap; transition: background 0.15s;
  }}
  .ctrl-btn:hover {{ background: #475569; color: #f1f5f9; }}

  /* ── Main layout ── */
  #main {{ flex: 1; display: flex; overflow: hidden; }}

  /* ── Info panel (left) ── */
  #info-panel {{
    width: 0; min-width: 0; background: #1e293b;
    border-right: 1px solid #334155;
    display: flex; flex-direction: column; overflow: hidden; flex-shrink: 0;
    transition: width 0.25s ease, min-width 0.25s ease;
  }}
  #info-panel.open {{ width: 360px; min-width: 260px; }}
  #info-panel-header {{
    padding: 14px 16px 10px; border-bottom: 1px solid #334155; flex-shrink: 0;
    display: flex; align-items: flex-start; justify-content: space-between; gap: 8px;
  }}
  #info-panel-title {{ font-size: 17px; font-weight: 700; color: #f1f5f9; }}
  #info-panel-close {{
    background: none; border: none; color: #64748b;
    font-size: 18px; cursor: pointer; line-height: 1; padding: 2px 4px; flex-shrink: 0;
  }}
  #info-panel-close:hover {{ color: #f1f5f9; }}
  #info-tabs {{
    display: flex; padding: 10px 16px 0;
    border-bottom: 1px solid #334155; flex-shrink: 0; flex-wrap: wrap;
  }}
  .info-tab {{
    background: none; border: none; border-bottom: 2px solid transparent;
    color: #64748b; font-size: 12px; font-weight: 600;
    text-transform: uppercase; letter-spacing: 0.06em;
    padding: 6px 10px 8px; cursor: pointer; white-space: nowrap;
    transition: color 0.15s, border-color 0.15s;
  }}
  .info-tab:hover {{ color: #cbd5e1; }}
  .info-tab.active {{ color: #a5b4fc; border-bottom-color: #6366f1; }}
  #info-body {{
    flex: 1; overflow-y: auto; padding: 16px;
    font-size: 13px; color: #cbd5e1; line-height: 1.65;
  }}
  #info-body::-webkit-scrollbar {{ width: 5px; }}
  #info-body::-webkit-scrollbar-thumb {{ background: #334155; border-radius: 4px; }}
  .info-section {{ display: none; }}
  .info-section.active {{ display: block; }}
  .info-h2 {{
    font-size: 12px; font-weight: 700; color: #94a3b8;
    text-transform: uppercase; letter-spacing: 0.07em; margin: 18px 0 6px;
  }}
  .info-h2:first-child {{ margin-top: 0; }}
  .info-p {{ margin: 0 0 10px; }}
  .info-stat-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 8px; margin-bottom: 12px; }}
  .info-stat {{
    background: #0f172a; border: 1px solid #334155; border-radius: 8px; padding: 10px 12px;
  }}
  .info-stat-val {{ font-size: 22px; font-weight: 700; color: #a5b4fc; line-height: 1.1; }}
  .info-stat-lbl {{ font-size: 11px; color: #64748b; margin-top: 3px; }}
  .info-code {{
    background: #0f172a; border: 1px solid #334155; border-radius: 6px;
    padding: 8px 10px; font-family: monospace; font-size: 11px; color: #7dd3fc;
    margin: 6px 0 10px; white-space: pre-wrap; word-break: break-all;
  }}
  .layer-row {{
    display: flex; align-items: center; gap: 8px;
    padding: 7px 0; border-bottom: 1px solid #1e293b;
  }}
  .layer-dot {{ width: 10px; height: 10px; border-radius: 50%; flex-shrink: 0; }}
  .layer-desc {{ font-size: 12px; color: #94a3b8; line-height: 1.45; }}

  /* ── Graph ── */
  #cy-wrap {{ flex: 1; position: relative; overflow: hidden; }}
  #cy {{ width: 100%; height: 100%; cursor: default; }}

  /* ── Sidebar ── */
  #sidebar {{
    width: 340px; min-width: 260px; background: #1e293b;
    border-left: 1px solid #334155;
    display: flex; flex-direction: column; overflow: hidden; flex-shrink: 0;
  }}
  #sidebar-header {{ padding: 14px 16px 10px; border-bottom: 1px solid #334155; flex-shrink: 0; }}
  #sidebar-title {{
    font-size: 12px; font-weight: 600; color: #94a3b8;
    text-transform: uppercase; letter-spacing: 0.06em; margin-bottom: 4px;
  }}
  #sidebar-name {{ font-size: 17px; font-weight: 700; color: #f1f5f9; word-break: break-word; }}
  #layer-tabs {{
    display: flex; border-bottom: 1px solid #334155; flex-shrink: 0;
    padding: 8px 16px 0; flex-wrap: wrap;
  }}
  .layer-tab {{
    background: none; border: none;
    border-bottom: 2px solid transparent;
    font-size: 11px; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em;
    padding: 4px 8px 6px; cursor: pointer; color: #64748b;
    transition: color 0.15s, border-color 0.15s; white-space: nowrap;
  }}
  .layer-tab:hover {{ color: #cbd5e1; }}
  .layer-tab.active {{ color: #f1f5f9; border-bottom-color: var(--ltab-color); }}
  #sidebar-body {{
    flex: 1; overflow-y: auto; padding: 14px 16px;
    font-size: 13px; color: #cbd5e1; line-height: 1.6;
  }}
  #sidebar-body::-webkit-scrollbar {{ width: 5px; }}
  #sidebar-body::-webkit-scrollbar-thumb {{ background: #334155; border-radius: 4px; }}
  .sec-label {{
    font-size: 11px; font-weight: 600; color: #64748b;
    text-transform: uppercase; letter-spacing: 0.07em; margin: 14px 0 6px;
  }}
  .sec-label:first-child {{ margin-top: 0; }}
  .layer-stat-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 6px; margin-bottom: 10px; }}
  .layer-stat {{
    background: #0f172a; border-radius: 8px; padding: 9px 11px; border-left: 3px solid;
  }}
  .layer-stat-n {{ font-size: 20px; font-weight: 700; line-height: 1; }}
  .layer-stat-lbl {{ font-size: 11px; color: #64748b; margin-top: 2px; }}
  .layer-stat.absent {{ opacity: 0.3; border-left-color: #334155 !important; }}
  .layer-stat.absent .layer-stat-n {{ color: #475569; }}
  .chips-wrap {{ display: flex; flex-wrap: wrap; gap: 5px; margin-bottom: 10px; }}
  .chip {{
    display: inline-flex; align-items: center; padding: 3px 9px;
    border-radius: 999px; font-size: 12px; font-weight: 500;
  }}
  .paper-item {{
    background: #0f172a; border: 1px solid #1e293b;
    border-radius: 7px; padding: 8px 10px; margin-bottom: 6px;
  }}
  .paper-title {{ font-size: 12px; color: #e2e8f0; line-height: 1.4; margin-bottom: 3px; }}
  .paper-meta {{ font-size: 11px; color: #64748b; }}
  .paper-link {{ color: #818cf8; text-decoration: none; font-size: 11px; }}
  .paper-link:hover {{ color: #a5b4fc; text-decoration: underline; }}
  .empty-state {{
    text-align: center; color: #475569; padding: 40px 20px; line-height: 1.7;
  }}
  .empty-state svg {{ margin: 0 auto 12px; display: block; opacity: 0.4; }}
  .faded {{ opacity: 0.06 !important; }}
  .highlighted {{ opacity: 1 !important; }}
  edge.search-hit {{ line-color: #fbbf24 !important; width: 5 !important; }}
</style>
</head>
<body>

<div id="topbar">
  <button id="info-btn" onclick="toggleInfoPanel()" title="How it works">i</button>
  <h1>{title}</h1>
  <div class="sep"></div>
  <div id="layer-toggles"></div>
  <div class="sep"></div>
  <div id="search-wrap">
    <span id="search-icon">🔍</span>
    <input id="search-input" type="text" placeholder="Search papers…" autocomplete="off">
    <div id="search-results"></div>
  </div>
  <span id="search-count"></span>
  <div class="sep"></div>
  <label style="display:flex;align-items:center;gap:6px;font-size:12px;color:#94a3b8;white-space:nowrap">
    Min papers
    <input id="min-papers-input" type="number" min="0" step="10" value="0"
      style="width:70px;padding:4px 7px;background:#0f172a;border:1px solid #334155;
             border-radius:6px;color:#e2e8f0;font-size:13px;outline:none"
      oninput="applyLayerVisibility();if(searchInput.value.trim())runSearch()">
  </label>
  <div class="sep"></div>
  <button class="ctrl-btn" onclick="runLayout('circle')">⭕ Circle</button>
  <button class="ctrl-btn" onclick="runLayout('grid')">⊞ Grid</button>
  <button class="ctrl-btn" onclick="cy.fit()">⤢ Fit</button>
</div>

<div id="main">

  <!-- Info panel (left) -->
  <div id="info-panel">
    <div id="info-panel-header">
      <div id="info-panel-title">How It Works</div>
      <button id="info-panel-close" onclick="toggleInfoPanel()">✕</button>
    </div>
    <div id="info-tabs">
      <button class="info-tab active" onclick="showInfoSection('search', event)">Literature Search</button>
      <button class="info-tab" onclick="showInfoSection('layers', event)">Layer Types</button>
    </div>
    <div id="info-body">

      <div id="info-section-search" class="info-section active">
        <div class="info-h2">Overview</div>
        <p class="info-p">
          For every organ-organ pair (91 pairs from 14 organs) the pipeline runs
          <strong>one broad PubMed query</strong>, then applies two post-processing steps:
          co-occurrence filtering and layer classification.
        </p>
        <div class="info-h2">Search Parameters</div>
        <div class="info-stat-grid">
          <div class="info-stat"><div class="info-stat-val">91</div><div class="info-stat-lbl">Organ pairs searched</div></div>
          <div class="info-stat"><div class="info-stat-val">91</div><div class="info-stat-lbl">Total queries (1 per pair)</div></div>
          <div class="info-stat"><div class="info-stat-val">10 yrs</div><div class="info-stat-lbl">PubMed look-back window</div></div>
          <div class="info-stat"><div class="info-stat-val">200</div><div class="info-stat-lbl">Max papers per pair</div></div>
        </div>
        <div class="info-h2">Step 1 — PubMed Query</div>
        <p class="info-p">Each query retrieves papers mentioning both organs plus a broad interaction filter:</p>
        <div class="info-code">(MeSH_A OR aliases_A)
AND (MeSH_B OR aliases_B)
AND (axis OR crosstalk OR interaction
    OR regulation OR signaling OR coupling)</div>
        <div class="info-h2">Step 2 — Co-occurrence Filter</div>
        <p class="info-p">
          A paper is kept only when the two organs co-occur in <strong>at least one</strong> of:
        </p>
        <ul style="margin:0 0 10px 16px;padding:0;color:#cbd5e1;line-height:1.9">
          <li><strong>Same sentence</strong> — both organ names appear in the same sentence of the abstract or title</li>
          <li><strong>Hyphenated compound</strong> — e.g. "liver-gut axis", "gut-liver crosstalk"</li>
          <li><strong>Merged compound word</strong> — e.g. "hepatorenal", "cardiorenal", "enterohepatic"</li>
        </ul>
        <div class="info-h2">Step 3 — Layer Classification</div>
        <p class="info-p">
          Each surviving paper is scanned for layer-specific mechanism phrases.
          A paper is assigned to every layer whose phrases it contains.
          Papers matching none of the five layer sets are placed in
          <strong>Undefined</strong>.
        </p>
        <div class="info-h2">Key Player Extraction</div>
        <p class="info-p">
          For each edge × layer, paper abstracts are scanned against a curated
          vocabulary. Terms are ranked by mention count and shown as chips.
        </p>
      </div>

      <div id="info-section-layers" class="info-section">
        <div class="info-h2">Communication Layers &amp; Paper Selection Criteria</div>
        <p class="info-p">Each layer button switches the visible edge set.
        A paper is counted for a layer only when its title or abstract contains
        <strong>at least one</strong> of the required phrases below, confirming
        it describes organ–organ communication of that specific type.</p>

        <div class="layer-row">
          <div class="layer-dot" style="background:#c084fc"></div>
          <div style="flex:1">
            <strong style="color:#d8b4fe">Neural</strong>
            <div class="layer-desc">Autonomic (sympathetic/parasympathetic), vagal, neuroendocrine, and neuronal signalling between organs.</div>
            <div class="info-code" style="margin-top:5px">autonomic innervation · sympathetic innervation
parasympathetic innervation · vagal regulation
vagal efferent · vagal afferent · neuroendocrine axis
neural regulation · autonomic regulation
nerve-mediated · neural crosstalk · neuronal regulation</div>
          </div>
        </div>

        <div class="layer-row">
          <div class="layer-dot" style="background:#fb923c"></div>
          <div style="flex:1">
            <strong style="color:#fdba74">Hormonal</strong>
            <div class="layer-desc">Circulating endocrine axes — peptide hormones, steroid hormones, hepatokines, adipokines, myokines, gut hormones, and other blood-borne signals.</div>
            <div class="info-code" style="margin-top:5px">endocrine axis · hormonal axis · circulating hormone
endocrine crosstalk · hormonal crosstalk
hormonal regulation · endocrine regulation
blood-borne signal · hormonal communication
endocrine communication · hormonal signaling · endocrine signaling</div>
          </div>
        </div>

        <div class="layer-row">
          <div class="layer-dot" style="background:#f87171"></div>
          <div style="flex:1">
            <strong style="color:#fca5a5">Immune</strong>
            <div class="layer-desc">Cytokine and chemokine signalling, immune-cell trafficking, inflammatory mediators, and immunomodulatory crosstalk.</div>
            <div class="info-code" style="margin-top:5px">immune crosstalk · inflammatory crosstalk
cytokine-mediated · immune-mediated
inflammatory mediator · immune regulation
cytokine signaling · immunomodulatory
immune communication · cytokine crosstalk · immune axis</div>
          </div>
        </div>

        <div class="layer-row">
          <div class="layer-dot" style="background:#34d399"></div>
          <div style="flex:1">
            <strong style="color:#5eead4">Metabolic</strong>
            <div class="layer-desc">Substrate exchange and inter-organ metabolic crosstalk — glucose, lipids, amino acids, ketone bodies, bile acids and other metabolic signals.</div>
            <div class="info-code" style="margin-top:5px">metabolic crosstalk · inter-organ · organ crosstalk
organ communication · substrate exchange
metabolic communication · metabolic axis
metabolic interplay · metabolic interaction · metabolic relay</div>
          </div>
        </div>

        <div class="layer-row">
          <div class="layer-dot" style="background:#60a5fa"></div>
          <div style="flex:1">
            <strong style="color:#93c5fd">Mechanical</strong>
            <div class="layer-desc">Hemodynamic coupling, blood-flow-mediated signalling, pressure and stretch sensing, baroreceptor pathways, and mechanotransduction.</div>
            <div class="info-code" style="margin-top:5px">hemodynamic coupling · blood flow-mediated · flow-mediated
pressure sensing · stretch sensing · pressure-mediated
baroreceptor · mechanotransduction · mechanosensing
vascular coupling · mechanical coupling · hemodynamic regulation</div>
          </div>
        </div>
        <div class="layer-row" style="border-bottom:none">
          <div class="layer-dot" style="background:#94a3b8"></div>
          <div style="flex:1">
            <strong style="color:#cbd5e1">Undefined</strong>
            <div class="layer-desc">Co-occurring papers that do not match any of the five mechanism phrase sets. May describe novel, mixed, or indirect inter-organ relationships not yet captured by the current classification.</div>
          </div>
        </div>
      </div>

    </div>
  </div>

  <div id="cy-wrap"><div id="cy"></div></div>

  <div id="sidebar">
    <div id="sidebar-header">
      <div id="sidebar-title">Select a node or edge</div>
      <div id="sidebar-name">—</div>
    </div>
    <div id="layer-tabs" style="display:none"></div>
    <div id="sidebar-body">
      <div class="empty-state">
        <svg width="44" height="44" viewBox="0 0 24 24" fill="none"
             stroke="currentColor" stroke-width="1.5">
          <circle cx="12" cy="12" r="10"/><path d="M12 8v4M12 16h.01"/>
        </svg>
        Click any <strong>node</strong> or <strong>edge</strong> to see details.<br><br>
        Toggle communication layers using the coloured buttons above.
      </div>
    </div>
  </div>
</div>

<script>
const elements    = {elements_json};
const LAYER_META  = {layer_meta_json};
const LAYER_NAMES = {layer_names_json};
const CHIP_STYLE  = {layer_chip_json};
const KP_LABEL    = {layer_kp_lbl_json};

// Only one layer active at a time (exclusive select)
let activeLayer = LAYER_NAMES[0];
let currentEdgeData = null;
let currentLayerTab = null;

// ── Cytoscape ──────────────────────────────────────────────────────────────────
const cy = cytoscape({{
  container: document.getElementById('cy'),
  elements,
  style: [
    {{
      selector: 'node',
      style: {{
        'background-color': 'data(color)', 'label': 'data(label)',
        'color': '#f1f5f9', 'text-valign': 'center', 'text-halign': 'center',
        'font-size': '11px', 'font-weight': '600', 'width': 52, 'height': 52,
        'border-width': 2, 'border-color': 'rgba(255,255,255,0.15)',
        'text-outline-color': 'rgba(0,0,0,0.6)', 'text-outline-width': 2,
        'transition-property': 'opacity, border-color, width, height',
        'transition-duration': '0.15s',
      }}
    }},
    {{
      selector: 'node.highlighted',
      style: {{ 'border-width': 3, 'border-color': '#fbbf24', 'width': 60, 'height': 60 }}
    }},
    {{
      selector: 'node.hovered',
      style: {{ 'border-width': 3, 'border-color': '#94a3b8', 'width': 58, 'height': 58 }}
    }},
    {{
      selector: 'edge',
      style: {{
        'width': 2.5, 'line-color': 'data(color)',
        'target-arrow-color': 'data(color)', 'target-arrow-shape': 'none',
        'curve-style': 'bezier', 'opacity': 1,
        'transition-property': 'opacity, width', 'transition-duration': '0.15s',
      }}
    }},
    {{
      selector: 'edge.highlighted',
      style: {{ 'width': 5, 'opacity': 1 }}
    }},
    {{
      selector: 'edge.hovered',
      style: {{ 'width': 5, 'opacity': 1 }}
    }},
    {{
      selector: 'edge.faded',
      style: {{ 'opacity': 0.05 }}
    }},
    {{
      selector: 'edge.hidden',
      style: {{ 'display': 'none' }}
    }},
    {{
      selector: 'edge.search-hit',
      style: {{ 'width': 5, 'opacity': 1 }}
    }},
  ],
  layout: {{ name: 'circle', animate: false, padding: 60 }},
}});

// ── Layer toggle buttons (exclusive single-select) ────────────────────────────
LAYER_NAMES.forEach(name => {{
  const meta = LAYER_META[name];
  const btn  = document.createElement('button');
  btn.className = 'layer-toggle';
  btn.id = 'ltoggle-' + name;
  btn.style.background    = meta.color + '22';
  btn.style.borderColor   = meta.color;
  btn.style.color         = meta.color;
  btn.innerHTML = `<span class="dot" style="background:${{meta.color}}"></span>${{meta.label}}`;
  btn.title     = `Show ${{meta.label}} connections`;
  btn.onclick   = () => toggleLayer(name);
  if (name !== activeLayer) btn.classList.add('off');
  document.getElementById('layer-toggles').appendChild(btn);
}});

function toggleLayer(name) {{
  activeLayer = name;
  LAYER_NAMES.forEach(l => {{
    const btn = document.getElementById('ltoggle-' + l);
    if (l === activeLayer) btn.classList.remove('off');
    else btn.classList.add('off');
  }});
  applyLayerVisibility();
  if (searchInput.value.trim()) runSearch();
}}

function applyLayerVisibility() {{
  const minPapers = parseInt(document.getElementById('min-papers-input').value) || 0;
  cy.edges().forEach(edge => {{
    const layers = edge.data('layers') || [];
    const n = edge.data('lyr_' + activeLayer + '_n') || 0;
    const visible = layers.includes(activeLayer) && n >= minPapers;
    if (visible) {{
      edge.removeClass('hidden');
    }} else {{
      edge.addClass('hidden');
      edge.removeClass('faded highlighted search-hit');
    }}
  }});
  cy.nodes().removeClass('hidden');
}}

applyLayerVisibility();

// ── Helpers ────────────────────────────────────────────────────────────────────
function runLayout(name) {{ cy.layout({{ name, animate: true, padding: 60 }}).run(); }}

function escHtml(s) {{
  return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}}

function chipsHtml(items, layer) {{
  if (!items || items.length === 0)
    return '<span style="color:#475569;font-size:12px">None identified</span>';
  const style = CHIP_STYLE[layer] || '';
  return '<div class="chips-wrap">' +
    items.map(t => `<span class="chip" style="${{style}}">${{escHtml(t)}}</span>`).join('') +
    '</div>';
}}

function papersHtml(papers) {{
  if (!papers || papers.length === 0)
    return '<p style="color:#475569;font-size:12px">No papers stored for this layer.</p>';
  return papers.map(p => {{
    const href = p.doi
      ? `https://doi.org/${{p.doi}}`
      : `https://pubmed.ncbi.nlm.nih.gov/${{p.pmid}}/`;
    return `<div class="paper-item">
      <div class="paper-title">${{escHtml(p.title || 'No title')}}</div>
      <div class="paper-meta">${{p.year || ''}} · <a class="paper-link" href="${{href}}"
        target="_blank" rel="noopener">PMID ${{p.pmid || '?'}}</a></div>
    </div>`;
  }}).join('');
}}

// ── Sidebar ────────────────────────────────────────────────────────────────────
function showEdgeSidebar(data) {{
  currentEdgeData = data;
  document.getElementById('sidebar-title').textContent = 'Connection';
  document.getElementById('sidebar-name').textContent  = data.source + ' ↔ ' + data.target;

  const layers = data.layers || [];
  if (!currentLayerTab || !layers.includes(currentLayerTab)) {{
    currentLayerTab = layers.includes(activeLayer) ? activeLayer : (layers[0] || null);
  }}

  // Build tabs
  const tabsEl = document.getElementById('layer-tabs');
  tabsEl.style.display = 'flex';
  tabsEl.innerHTML = layers.map(l => {{
    const meta = LAYER_META[l];
    return `<button class="layer-tab ${{l === currentLayerTab ? 'active' : ''}}"
      style="--ltab-color:${{meta.color}}"
      onclick="switchLayerTab('${{l}}')">${{meta.label}}</button>`;
  }}).join('');

  renderEdgeLayerContent(data);
}}

function switchLayerTab(layer) {{
  currentLayerTab = layer;
  showEdgeSidebar(currentEdgeData);
}}

function renderEdgeLayerContent(data) {{
  const layers = data.layers || [];
  const layer  = currentLayerTab;
  const meta   = layer ? LAYER_META[layer] : null;
  let html = '';

  // Paper count grid for all 5 layers
  html += `<div class="sec-label">Papers found per layer</div><div class="layer-stat-grid">`;
  LAYER_NAMES.forEach(l => {{
    const n   = data['lyr_' + l + '_n'] || 0;
    const m   = LAYER_META[l];
    const has = layers.includes(l);
    html += `<div class="layer-stat ${{has ? '' : 'absent'}}" style="border-left-color:${{m.color}}">
      <div class="layer-stat-n" style="color:${{has ? m.color : '#475569'}}">${{n}}</div>
      <div class="layer-stat-lbl">${{m.label}}</div>
    </div>`;
  }});
  html += '</div>';

  // Key players for selected tab
  if (layer && layers.includes(layer)) {{
    const kp     = data['lyr_' + layer + '_kp'] || [];
    const papers = data['lyr_' + layer + '_papers'] || [];
    html += `<div class="sec-label" style="color:${{meta.color}}">${{KP_LABEL[layer] || 'Key Players'}}</div>`;
    html += chipsHtml(kp, layer);
    html += `<div class="sec-label" style="color:${{meta.color}}">Top Papers</div>`;
    html += papersHtml(papers);
  }}

  document.getElementById('sidebar-body').innerHTML = html;
}}

function showNodeSidebar(data) {{
  currentEdgeData = null;
  currentLayerTab = null;
  document.getElementById('sidebar-title').textContent = 'Organ';
  document.getElementById('sidebar-name').textContent  = data.label;
  document.getElementById('layer-tabs').style.display  = 'none';

  const edges = cy.getElementById(data.id).connectedEdges().filter(e => (e.data('layers')||[]).includes(activeLayer));
  const layerCounts = {{}};
  LAYER_NAMES.forEach(l => layerCounts[l] = 0);
  edges.forEach(e => (e.data('layers') || []).forEach(l => layerCounts[l]++));

  let html = '';
  if (data.description) {{
    html += `<div class="sec-label">Description</div>
    <p style="color:#cbd5e1;line-height:1.65;margin-bottom:12px">${{escHtml(data.description)}}</p>`;
  }}
  html += `<div class="sec-label">Connections by layer (visible)</div><div class="layer-stat-grid">`;
  LAYER_NAMES.forEach(l => {{
    const n = layerCounts[l];
    const m = LAYER_META[l];
    html += `<div class="layer-stat ${{n === 0 ? 'absent' : ''}}" style="border-left-color:${{m.color}}">
      <div class="layer-stat-n" style="color:${{n > 0 ? m.color : '#475569'}}">${{n}}</div>
      <div class="layer-stat-lbl">${{m.label}}</div>
    </div>`;
  }});
  html += '</div>';
  document.getElementById('sidebar-body').innerHTML = html;
}}

// ── Click events ───────────────────────────────────────────────────────────────
cy.on('tap', 'edge', function(evt) {{
  const el = evt.target;
  if (el.hasClass('hidden')) return;
  if (!searchInput.value.trim()) {{
    cy.edges().removeClass('faded highlighted');
    cy.nodes().removeClass('highlighted');
    cy.edges().not('.hidden').addClass('faded');
    el.removeClass('faded').addClass('highlighted');
    el.connectedNodes().addClass('highlighted');
  }}
  showEdgeSidebar(el.data());
}});

cy.on('tap', 'node', function(evt) {{
  const el = evt.target;
  if (!searchInput.value.trim()) {{
    cy.edges().removeClass('faded highlighted');
    cy.nodes().removeClass('highlighted');
    cy.edges().not('.hidden').addClass('faded');
    el.connectedEdges().not('.hidden').removeClass('faded').addClass('highlighted');
    el.addClass('highlighted');
    el.connectedNodes().addClass('highlighted');
  }}
  showNodeSidebar(el.data());
}});

cy.on('tap', function(evt) {{
  if (evt.target === cy) {{
    if (!searchInput.value.trim()) {{
      cy.edges().removeClass('faded highlighted');
      cy.nodes().removeClass('highlighted');
    }}
    currentEdgeData = null;
    document.getElementById('sidebar-title').textContent = 'Select a node or edge';
    document.getElementById('sidebar-name').textContent  = '—';
    document.getElementById('layer-tabs').style.display  = 'none';
    document.getElementById('sidebar-body').innerHTML = `<div class="empty-state">
      <svg width="44" height="44" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
        <circle cx="12" cy="12" r="10"/><path d="M12 8v4M12 16h.01"/></svg>
      Click any <strong>node</strong> or <strong>edge</strong> to see details.
    </div>`;
  }}
}});

// ── Hover effects ──────────────────────────────────────────────────────────────
cy.on('mouseover', 'node', function(evt) {{
  evt.target.addClass('hovered');
  document.getElementById('cy').style.cursor = 'pointer';
}});
cy.on('mouseover', 'edge', function(evt) {{
  if (!evt.target.hasClass('hidden')) {{
    evt.target.addClass('hovered');
    document.getElementById('cy').style.cursor = 'pointer';
  }}
}});
cy.on('mouseout', 'node, edge', function(evt) {{
  evt.target.removeClass('hovered');
  document.getElementById('cy').style.cursor = 'default';
}});

// ── Info panel ─────────────────────────────────────────────────────────────────
function toggleInfoPanel() {{
  document.getElementById('info-panel').classList.toggle('open');
  document.getElementById('info-btn').classList.toggle('active');
}}

function showInfoSection(name, evt) {{
  document.querySelectorAll('.info-section').forEach(s => s.classList.remove('active'));
  document.querySelectorAll('.info-tab').forEach(t => t.classList.remove('active'));
  document.getElementById('info-section-' + name).classList.add('active');
  if (evt) evt.currentTarget.classList.add('active');
}}

// ── Search ─────────────────────────────────────────────────────────────────────
const searchInput   = document.getElementById('search-input');
const searchResults = document.getElementById('search-results');
const searchCount   = document.getElementById('search-count');

const searchIndex = [];
cy.edges().forEach(e => {{
  const d = e.data();
  const kpText = [];
  LAYER_NAMES.forEach(l => {{
    (d['lyr_' + l + '_kp'] || []).forEach(k => kpText.push(k.toLowerCase()));
  }});
  searchIndex.push({{
    el: e, data: d,
    kpText: kpText.join(' '),
    label: d.source + ' ↔ ' + d.target,
  }});
}});

// Return the mention count of a query in the active layer's kp_counts for an edge.
// Handles partial matches (e.g. "bile" finds "bile acid", "bile salt").
function getKpTermCount(d, q) {{
  const counts = d['lyr_' + activeLayer + '_kp_counts'] || {{}};
  let best = 0;
  for (const [term, n] of Object.entries(counts)) {{
    if (term.toLowerCase().includes(q)) best = Math.max(best, n);
  }}
  return best;
}}

function runSearch() {{
  const q = searchInput.value.trim().toLowerCase();
  if (!q) {{
    cy.edges().not('.hidden').removeClass('faded search-hit');
    cy.nodes().removeClass('faded');
    searchResults.style.display = 'none';
    searchCount.textContent = '';
    return;
  }}
  const minPapers = parseInt(document.getElementById('min-papers-input').value) || 0;
  const hits = searchIndex.filter(entry => {{
    const d = entry.data;
    const layers = d.layers || [];
    if (!layers.includes(activeLayer)) return false;
    const layerN = d['lyr_' + activeLayer + '_n'] || 0;
    if (layerN < minPapers) return false;
    if (!entry.kpText.includes(q)) return false;
    // Threshold: the searched term must appear ≥ minPapers times as a key player
    if (minPapers > 0 && getKpTermCount(d, q) < minPapers) return false;
    return true;
  }});
  cy.nodes().removeClass('faded');
  cy.edges().not('.hidden').removeClass('search-hit').addClass('faded');
  hits.forEach(h => h.el.removeClass('faded').addClass('search-hit'));
  searchCount.textContent = `${{hits.length}} edge${{hits.length !== 1 ? 's' : ''}}`;
  if (hits.length === 0) {{
    searchResults.innerHTML = '<div class="sr-item" style="color:#64748b">No results</div>';
  }} else {{
    searchResults.innerHTML = hits.slice(0, 20).map(h => {{
      const n = getKpTermCount(h.data, q);
      const badge = n ? `<span style="margin-left:auto;background:#1e3a5f;color:#7dd3fc;border-radius:8px;padding:1px 7px;font-size:11px;font-weight:700">${{n}}×</span>` : '';
      return `<div class="sr-item" data-id="${{h.el.id()}}" style="display:flex;align-items:center;gap:8px">🔗 <span style="flex:1">${{escHtml(h.label)}}</span>${{badge}}</div>`;
    }}).join('');
    searchResults.querySelectorAll('.sr-item[data-id]').forEach(item => {{
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
searchInput.addEventListener('focus', () => {{ if (searchInput.value.trim()) searchResults.style.display = 'block'; }});
document.addEventListener('click', e => {{
  if (!document.getElementById('search-wrap').contains(e.target)) searchResults.style.display = 'none';
  if (!document.getElementById('info-panel').contains(e.target) &&
      !document.getElementById('info-btn').contains(e.target)) {{
    // leave panel open — user closes with ✕ or info-btn
  }}
}});
</script>
</body>
</html>"""

    output_html.write_text(html, encoding="utf-8")
    n_edges = sum(1 for el in elements if "source" in el["data"])
    n_nodes = sum(1 for el in elements if "source" not in el["data"])
    print(f"[✔] Saved to {output_html}  ({n_nodes} organs, {n_edges} edges)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-papers", type=int, default=DEFAULT_MIN_PAPERS)
    args = parser.parse_args()

    if not DEFAULT_SEARCH_RESULTS.exists():
        print(f"[!] No search results at {DEFAULT_SEARCH_RESULTS}")
        print("    Run:  uv run python run_lit_ref.py")
        return

    import json as _j
    with open(DEFAULT_SEARCH_RESULTS, encoding="utf-8") as f:
        results = _j.load(f)

    from Data_Loader.load_data import load_node_metadata_from_csv
    node_meta = load_node_metadata_from_csv(str(HERE / "metabolic_data" / "organ_data.csv"))

    build_lit_ref_viz(results, node_meta, DEFAULT_OUTPUT_HTML, args.min_papers)


if __name__ == "__main__":
    main()
