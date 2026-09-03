"""
run_comparison.py
=================
Compares the condition's literature-based network (bootstrapped from PubMed
per organ) against the Reference Network Only Metabolic, restricted to the
organ pairs marked in that condition's cohort connections CSV:
  healthy -> healthy/healthy_cohort_connections.csv
  obese   -> obese/obese_cohort_connections.csv

Each edge falls into one of three categories:
  SHARED      — present in both networks
  ONLY_B      — robust in the literature-based network but not a predefined reference edge
  ONLY_REF    — predefined reference edge with papers, but below the literature-based network threshold

Edge colours in the visualization:
  Teal   (#10b981) — SHARED
  Blue   (#38bdf8) — ONLY_B
  Amber  (#f59e0b) — ONLY_REF

Run from the Metabolic_Reference_Network/ directory:
    uv run -m Edge_cosine_met_reference_network.run_comparison --condition healthy
    uv run -m Edge_cosine_met_reference_network.run_comparison --condition obese

    --viz-only   rebuild HTML from existing comparison JSON without reloading
"""

import sys
import csv
import json
import argparse
import importlib.util
from pathlib import Path

HERE = Path(__file__).parent
ROOT = HERE.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Windows console/redirected-output encoding defaults to cp1252, which can't
# encode the arrows/ellipses used in progress prints below — force UTF-8 so
# the run doesn't crash mid-way through (e.g. when stdout is piped to a file).
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8", errors="replace")

from threshold_utils import Elbow, resolve_min_bootstrap_mean

REF_JSON   = ROOT / "reference_network_only_metabolic" / "metabolic_literature_results.json"
REF_EDGES_CSV = ROOT / "metabolic_data" / "connection_data.csv"

CONDITION_CONFIGS = {
    "healthy": HERE / "healthy" / "config_healthy.py",
    "obese":   HERE / "obese"   / "config_obese.py",
}

# Cohort connections CSV that defines which organ pairs count as reference
# edges for each condition.
COHORT_CONNECTIONS_CSV = {
    "healthy": HERE / "healthy" / "healthy_cohort_connections.csv",
    "obese":   HERE / "obese"   / "obese_cohort_connections.csv",
}

COLOR_SHARED   = "#10b981"   # teal
COLOR_ONLY_B   = "#38bdf8"   # sky-blue
COLOR_ONLY_REF = "#f59e0b"   # amber


# ── Config loader ─────────────────────────────────────────────────────────────

def _load_config(condition: str):
    path = CONDITION_CONFIGS[condition]
    spec = importlib.util.spec_from_file_location(f"_cfg_{condition}", path)
    mod  = importlib.util.module_from_spec(spec)
    mod.Elbow = Elbow  # lets the config file write `MIN_BOOTSTRAP_MEAN = Elbow`
    spec.loader.exec_module(mod)
    return mod


# ── Edge set helpers ──────────────────────────────────────────────────────────

def _canonical(o1: str, o2: str) -> tuple:
    return (min(o1, o2), max(o1, o2))


def load_cohort_edge_filter(csv_path: Path) -> set[tuple[str, str]]:
    """
    Parse a 0/1 adjacency matrix CSV.  Returns the set of canonical organ
    pairs where the cell value is 1.  Both upper and lower triangle are
    accepted; diagonal and empty cells are ignored.
    """
    pairs: set[tuple[str, str]] = set()
    with open(csv_path, encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        col_organs = header[1:]
        for row in reader:
            if not row:
                continue
            row_organ = row[0].strip()
            for col_idx, val in enumerate(row[1:], start=0):
                val = val.strip()
                if val == "1" and col_idx < len(col_organs):
                    col_organ = col_organs[col_idx].strip()
                    if row_organ and col_organ and row_organ != col_organ:
                        pairs.add((min(row_organ, col_organ),
                                   max(row_organ, col_organ)))
    return pairs


def load_cohort_organs(csv_path: Path) -> set[str]:
    """Return the organ names from the header row of a cohort connections CSV."""
    with open(csv_path, encoding="utf-8") as f:
        header = next(csv.reader(f))
    return {h.strip() for h in header[1:] if h.strip()}


def load_ref_edges(ref_json_path: Path, allowed_pairs: set[tuple[str, str]],
                   min_papers: int = 1) -> dict:
    """
    Return a dict of canonical edge -> ref data for all predefined edges
    that are present in allowed_pairs (the condition's cohort connections
    CSV) and have at least min_papers in the reference metabolic network.
    """
    if not ref_json_path.exists():
        print(f"[!] Reference results not found: {ref_json_path}")
        return {}
    with open(ref_json_path, encoding="utf-8") as f:
        raw = json.load(f)

    edges = {}
    for key, data in raw.items():
        o1, o2 = data["organ1"], data["organ2"]
        pair = _canonical(o1, o2)
        if pair not in allowed_pairs:
            continue
        if data.get("n_papers_found", 0) >= min_papers:
            edges[pair] = data
    return edges


def load_network_edges(condition: str, cfg, allowed_organs: set[str]) -> dict:
    """
    Return a dict of canonical edge -> bootstrap data for all literature-based
    network edges between two cohort organs that meet the MIN_BOOTSTRAP_MEAN
    threshold.
    """
    boot_json = HERE / condition / f"bootstrap_results_{condition}.json"
    if not boot_json.exists():
        print(f"[!] Bootstrap results not found: {boot_json}")
        print(f"    Run:  uv run -m Edge_cosine_met_reference_network.run_network "
              f"--condition {condition}")
        return {}
    with open(boot_json, encoding="utf-8") as f:
        raw = json.load(f)

    # Resolve MIN_BOOTSTRAP_MEAN = Elbow using this condition's cohort-restricted
    # bootstrap means, and mutate cfg in place so callers see the real number too.
    cohort_means = [data.get("mean", 0) for data in raw.values()
                    if data.get("organ1") in allowed_organs
                    and data.get("organ2") in allowed_organs]
    cfg.MIN_BOOTSTRAP_MEAN = resolve_min_bootstrap_mean(cfg.MIN_BOOTSTRAP_MEAN, cohort_means)

    edges = {}
    for key, data in raw.items():
        o1, o2 = data["organ1"], data["organ2"]
        if o1 not in allowed_organs or o2 not in allowed_organs:
            continue
        if data.get("mean", 0) >= cfg.MIN_BOOTSTRAP_MEAN:
            edges[_canonical(o1, o2)] = data
    return edges


def load_network_search(condition: str) -> dict:
    """Load per-organ search results for the literature-based network.  Keyed by organ name."""
    search_json = HERE / condition / f"search_results_{condition}.json"
    if not search_json.exists():
        return {}
    with open(search_json, encoding="utf-8") as f:
        return json.load(f)


# ── Comparison logic ──────────────────────────────────────────────────────────

def compare(ref_edges: dict, b_edges: dict) -> dict:
    """
    Returns a dict of canonical edge -> category ('shared'|'only_b'|'only_ref').
    """
    all_edges = set(ref_edges) | set(b_edges)
    categories = {}
    for edge in all_edges:
        in_ref = edge in ref_edges
        in_b   = edge in b_edges
        if in_ref and in_b:
            categories[edge] = "shared"
        elif in_b:
            categories[edge] = "only_b"
        else:
            categories[edge] = "only_ref"
    return categories


# ── Visualization ─────────────────────────────────────────────────────────────

def _kp_from_bootstrap(kp_ranked: dict, cat: str):
    """Extract term list and integer-percentage count dict from bootstrap key-player data."""
    entries = kp_ranked.get(cat, [])
    terms   = [e["term"] for e in entries]
    counts  = {e["term"]: int(round(e["freq"] * 100)) for e in entries}
    return terms, counts


def compute_comparison_data(condition: str, cfg, min_ref_papers: int = 1) -> dict | None:
    """
    Loads and classifies edges for the condition vs reference comparison.
    Returns None if the cohort CSV or reference/bootstrap data isn't ready,
    else a dict with categories/ref_edges/b_edges/b_search/allowed_organs/
    cohort_csv_name — everything run_network.py needs to recolor its own
    robust-network graph in place (comparison toggle) instead of rendering a
    separate comparison page.
    """
    cohort_csv = COHORT_CONNECTIONS_CSV[condition]
    if not cohort_csv.exists():
        return None
    allowed_pairs  = load_cohort_edge_filter(cohort_csv)
    allowed_organs = load_cohort_organs(cohort_csv)

    ref_edges = load_ref_edges(REF_JSON, allowed_pairs, min_papers=min_ref_papers)
    b_edges   = load_network_edges(condition, cfg, allowed_organs)
    b_search  = load_network_search(condition)
    if not ref_edges and not b_edges:
        return None

    categories = compare(ref_edges, b_edges)
    return {
        "categories":      categories,
        "ref_edges":       ref_edges,
        "b_edges":         b_edges,
        "b_search":        b_search,
        "allowed_organs":  allowed_organs,
        "cohort_csv_name": cohort_csv.name,
    }


def build_comparison_html_string(condition: str, cfg, min_ref_papers: int = 1) -> str | None:
    """
    Full comparison-vs-reference pipeline (load, classify, render) returning
    the dashboard as a standalone HTML string. Kept for CLI/standalone use;
    the robust network itself no longer calls this — it consumes
    compute_comparison_data() directly and recolors its own graph in place
    via the comparison toggle. Returns None if data isn't ready.
    """
    data = compute_comparison_data(condition, cfg, min_ref_papers=min_ref_papers)
    if data is None:
        return None
    return _build_comparison_html_string(data["categories"], data["ref_edges"], data["b_edges"],
                                         data["b_search"], cfg, data["cohort_csv_name"],
                                         data["allowed_organs"])


def _build_comparison_html_string(categories: dict, ref_edges: dict, b_edges: dict,
                                  b_search: dict, cfg, cohort_csv_name: str,
                                  allowed_organs: set[str]) -> str | None:
    try:
        import networkx as nx
        from Visualisation.networkBuilderUtils import export_network_to_cytoscape_dashboard

        G = nx.Graph()
        for organ in allowed_organs:
            G.add_node(organ)

        for (o1, o2), category in categories.items():
            color = {
                "shared":   COLOR_SHARED,
                "only_b":   COLOR_ONLY_B,
                "only_ref": COLOR_ONLY_REF,
            }[category]

            label_map = {
                "shared":   f"SHARED — in both networks",
                "only_b":   f"ONLY in {cfg.VIZ_LABEL} Literature-Based Network",
                "only_ref": "ONLY in Reference Metabolic Network",
            }

            # Pull papers and stats from whichever network has them
            if category in ("shared", "only_b"):
                b_data  = b_edges.get((o1, o2), {})
                papers  = b_data.get("papers", [])
                extra_stats = (
                    f"Otsuka–Ochiai coefficient: {b_data.get('mean', 0):.5f} ± "
                    f"{b_data.get('std', 0):.5f}  |  "
                    f"Papers (co-occur): {b_data.get('n_cooccur_total', 0)} / "
                    f"{b_data.get('n_found', 0)} found"
                )
                ref_data  = ref_edges.get((o1, o2), {})
                ref_stats = (f"  |  Reference n_papers: {ref_data.get('n_papers_found', '—')}"
                             if category == "shared" else "")
                description = f"{label_map[category]}\n{extra_stats}{ref_stats}"

                kp_ranked = b_data.get("key_players_bootstrap", {})
                h_terms, h_counts = _kp_from_bootstrap(kp_ranked, "hormones")
                m_terms, m_counts = _kp_from_bootstrap(kp_ranked, "metabolites")
                p_terms, p_counts = _kp_from_bootstrap(kp_ranked, "proteins")
                has_kp = bool(kp_ranked)
            else:  # only_ref
                ref_data    = ref_edges.get((o1, o2), {})
                papers      = ref_data.get("papers", [])
                description = (
                    f"{label_map[category]}\n"
                    f"Reference n_papers: {ref_data.get('n_papers_found', 0)}"
                )
                h_terms, h_counts = [], {}
                m_terms, m_counts = [], {}
                p_terms, p_counts = [], {}
                has_kp = False

            G.add_edge(o1, o2)
            G.edges[o1, o2]['color']       = color
            G.edges[o1, o2]['description'] = description
            G.edges[o1, o2]['merged_data'] = {
                "n_papers_found":   (b_edges.get((o1, o2), {}).get("n_found", 0)
                                     or ref_edges.get((o1, o2), {}).get("n_papers_found", 0)),
                "n_papers_cooccur": b_edges.get((o1, o2), {}).get("n_cooccur_total", 0),
                "bootstrap_mean":   b_edges.get((o1, o2), {}).get("mean", 0),
                "papers":           papers,
                "key_players_hormones":           h_terms,
                "key_players_metabolites":        m_terms,
                "key_players_proteins":           p_terms,
                "key_players_counts_hormones":    h_counts,
                "key_players_counts_metabolites": m_counts,
                "key_players_counts_proteins":    p_counts,
                "key_players_bootstrap":          has_kp,
                "connection_type":  label_map[category],
                "pubmed_query":     (b_search.get(o1, {}).get("query", "")
                                     or ref_edges.get((o1, o2), {}).get("pubmed_query", "")),
            }

        return export_network_to_cytoscape_dashboard(
            graph=G,
            filename=None,
            include_legend=False,
            title=f"Network Comparison — {cfg.VIZ_LABEL} vs Reference Metabolic",
            info_panel_tabs=_info_tabs(categories, cfg, cohort_csv_name),
        )
    except Exception as e:
        import traceback
        print(f"[!] Comparison visualization failed: {e}")
        traceback.print_exc()
        return None


def _organ_breakdown_tab(categories: dict, cfg) -> dict:
    """Per-organ tab: a horizontal stacked bar per organ showing how many of
    its connections are shared, only in the literature-based network, or
    only in the reference network."""
    from collections import defaultdict

    lit_label = f"{cfg.VIZ_LABEL} Literature-Based Network"
    counts: dict[str, dict[str, int]] = defaultdict(
        lambda: {"shared": 0, "only_b": 0, "only_ref": 0})
    for (o1, o2), c in categories.items():
        counts[o1][c] += 1
        counts[o2][c] += 1

    organs    = sorted(counts.keys())
    max_total = max((sum(v.values()) for v in counts.values()), default=0) or 1
    bar_w     = 200  # px, full bar width at max_total

    rows = []
    for organ in organs:
        c = counts[organ]
        total = sum(c.values())
        scale = bar_w * (total / max_total) if total else 0
        seg_shared = (c["shared"]   / total) * scale if total else 0
        seg_b      = (c["only_b"]   / total) * scale if total else 0
        seg_ref    = (c["only_ref"] / total) * scale if total else 0
        rows.append(f"""
        <div style="display:flex;align-items:center;gap:8px;margin-bottom:7px">
          <div style="width:108px;font-size:0.74rem;color:#cbd5e1;text-align:right;flex-shrink:0">{organ}</div>
          <div style="display:flex;height:14px;background:#0f172a;border-radius:3px;overflow:hidden;width:{bar_w}px;flex-shrink:0">
            <div style="width:{seg_shared:.1f}px;background:{COLOR_SHARED}" title="Shared: {c['shared']}"></div>
            <div style="width:{seg_b:.1f}px;background:{COLOR_ONLY_B}" title="Only in {lit_label}: {c['only_b']}"></div>
            <div style="width:{seg_ref:.1f}px;background:{COLOR_ONLY_REF}" title="Only in reference: {c['only_ref']}"></div>
          </div>
          <div style="font-size:0.72rem;color:#64748b;flex-shrink:0;white-space:nowrap">
            {c['shared']} / {c['only_b']} / {c['only_ref']}
          </div>
        </div>""")

    return {
        "id": "organ_breakdown",
        "label": "Per-Organ Breakdown",
        "content": f"""
    <div class="info-h2">Connections per Organ</div>
    <p class="info-p">
      For each organ, how many of its connections are <strong>shared</strong>
      between the <strong>{lit_label}</strong> and the reference network, found
      <strong>only in the {lit_label}</strong>, or <strong>only in the
      reference network</strong>. The trailing numbers indicate Shared / Only-Literature / Only-Reference.
    </p>
    <div style="display:flex;align-items:center;gap:14px;margin-bottom:14px;font-size:0.74rem;color:#94a3b8;flex-wrap:wrap">
      <span><span style="display:inline-block;width:12px;height:10px;background:{COLOR_SHARED};border-radius:2px;vertical-align:middle;margin-right:4px"></span>Shared</span>
      <span><span style="display:inline-block;width:12px;height:10px;background:{COLOR_ONLY_B};border-radius:2px;vertical-align:middle;margin-right:4px"></span>Only {lit_label}</span>
      <span><span style="display:inline-block;width:12px;height:10px;background:{COLOR_ONLY_REF};border-radius:2px;vertical-align:middle;margin-right:4px"></span>Only Reference</span>
    </div>
    {''.join(rows) if rows else '<p style="color:#64748b;font-size:0.82rem">No organs with connections.</p>'}
""",
    }


def _info_tabs(categories: dict, cfg, cohort_csv_name: str) -> list[dict]:
    shared   = [(o1, o2) for (o1, o2), c in categories.items() if c == "shared"]
    only_b   = [(o1, o2) for (o1, o2), c in categories.items() if c == "only_b"]
    only_ref = [(o1, o2) for (o1, o2), c in categories.items() if c == "only_ref"]
    lit_label = f"{cfg.VIZ_LABEL} Literature-Based Network"

    def _pair_list(pairs: list, color: str) -> str:
        if not pairs:
            return '<p style="color:#64748b;font-size:0.82rem">None</p>'
        items = "".join(
            f'<li style="margin-bottom:3px">{o1} ↔ {o2}</li>'
            for o1, o2 in sorted(pairs)
        )
        return f'<ul style="margin:0 0 10px 16px;padding:0;color:{color};font-size:0.82rem">{items}</ul>'

    comparison_tab = {
        "id": "comparison",
        "label": "Comparison",
        "content": f"""
    <div class="info-h2">What This Shows</div>
    <p class="info-p">
      This visualization compares the <strong>{lit_label}</strong>
      against the <strong>Reference Metabolic Network</strong>, which is based on
      the partial correlation network of healthy subjects (n= 241) with BMI above 24.4 kg/m^2.
      In this case {len(shared) + len(only_ref)} predefined organ-pair edges.
    </p>
    <p class="info-p" style="color:#94a3b8">
      The counts and lists below are a at the elbow-suggested Otsuka-Ochiai threshold.
    </p>
    <div class="info-h2">Edge Legend</div>
    <div style="display:flex;flex-direction:column;gap:8px;margin-bottom:16px">
      <div style="display:flex;align-items:center;gap:10px">
        <div style="width:28px;height:6px;border-radius:3px;background:{COLOR_SHARED}"></div>
        <span style="color:#e2e8f0;font-size:0.83rem">
          <strong>Shared</strong> — in both networks ({len(shared)} edges)
        </span>
      </div>
      <div style="display:flex;align-items:center;gap:10px">
        <div style="width:28px;height:6px;border-radius:3px;background:{COLOR_ONLY_B}"></div>
        <span style="color:#e2e8f0;font-size:0.83rem">
          <strong>Only in {lit_label}</strong> — robust in bootstrap,
          not a predefined reference edge ({len(only_b)} edges)
        </span>
      </div>
      <div style="display:flex;align-items:center;gap:10px">
        <div style="width:28px;height:6px;border-radius:3px;background:{COLOR_ONLY_REF}"></div>
        <span style="color:#e2e8f0;font-size:0.83rem">
          <strong>Only in Reference</strong> — predefined edge with papers,
          below {cfg.VIZ_LABEL} Otsuka–Ochiai threshold ({len(only_ref)} edges)
        </span>
      </div>
    </div>
    <div class="info-h2">{lit_label} Threshold</div>
    <p class="info-p">
      An edge appears in the {lit_label} when its mean Otsuka–Ochiai coefficient
      ≥ <strong>{cfg.MIN_BOOTSTRAP_MEAN}</strong>
      (across {cfg.N_BOOTSTRAP} × {int(cfg.SAMPLE_FRACTION*100)}% random samples).
    </p>
    <div class="info-h2">Shared Edges ({len(shared)})</div>
    {_pair_list(shared, COLOR_SHARED)}
    <div class="info-h2">Only in {lit_label} ({len(only_b)})</div>
    {_pair_list(only_b, COLOR_ONLY_B)}
    <div class="info-h2">Only in Reference Network ({len(only_ref)})</div>
    {_pair_list(only_ref, COLOR_ONLY_REF)}
""",
    }

    return [comparison_tab, _organ_breakdown_tab(categories, cfg)]


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--condition", required=True,
                        choices=list(CONDITION_CONFIGS),
                        help="Which condition's literature-based network to compare against the reference.")
    parser.add_argument("--min-ref-papers", type=int, default=1,
                        help="Min papers in reference network for an edge to count (default: 1).")
    args = parser.parse_args()

    cfg = _load_config(args.condition)

    cohort_csv = COHORT_CONNECTIONS_CSV[cfg.CONDITION_NAME]
    if not cohort_csv.exists():
        print(f"[!] Cohort connections file not found: {cohort_csv}")
        sys.exit(1)
    allowed_pairs  = load_cohort_edge_filter(cohort_csv)
    allowed_organs = load_cohort_organs(cohort_csv)

    print(f"\n[i] Condition : {cfg.VIZ_LABEL}")
    print(f"[i] Reference cohort filter: {cohort_csv.name} → "
          f"{len(allowed_organs)} organs, {len(allowed_pairs)} pairs")
    print(f"[i] Reference min papers: {args.min_ref_papers}\n")

    # Load both networks (load_network_edges resolves MIN_BOOTSTRAP_MEAN =
    # Elbow in-place on cfg, so the threshold print below shows the real value)
    ref_edges = load_ref_edges(REF_JSON, allowed_pairs, min_papers=args.min_ref_papers)
    b_edges   = load_network_edges(cfg.CONDITION_NAME, cfg, allowed_organs)
    b_search  = load_network_search(cfg.CONDITION_NAME)

    print(f"[i] {cfg.VIZ_LABEL} Literature-Based Network threshold : mean ≥ {cfg.MIN_BOOTSTRAP_MEAN}")

    if not ref_edges and not b_edges:
        print("[!] Both networks are empty — nothing to compare.")
        return

    print(f"[i] Reference edges  : {len(ref_edges)}")
    print(f"[i] {cfg.VIZ_LABEL} Literature-Based Network edges  : {len(b_edges)}")

    # Classify
    categories = compare(ref_edges, b_edges)
    shared   = sum(1 for c in categories.values() if c == "shared")
    only_b   = sum(1 for c in categories.values() if c == "only_b")
    only_ref = sum(1 for c in categories.values() if c == "only_ref")

    print(f"\n[i] Shared              : {shared}")
    print(f"[i] Only in Literature-Based Network : {only_b}")
    print(f"[i] Only in Reference   : {only_ref}")

    # Print details
    if shared:
        print("\n  Shared edges:")
        for (o1, o2), c in sorted(categories.items()):
            if c == "shared":
                print(f"    {o1} <-> {o2}")
    if only_b:
        print("\n  Only in Literature-Based Network:")
        for (o1, o2), c in sorted(categories.items()):
            if c == "only_b":
                print(f"    {o1} <-> {o2}")
    if only_ref:
        print("\n  Only in Reference:")
        for (o1, o2), c in sorted(categories.items()):
            if c == "only_ref":
                print(f"    {o1} <-> {o2}")

    print(f"\n[i] This comparison view is embedded in "
          f"robust_network_{cfg.CONDITION_NAME}.html (no separate file is "
          f"written). Rebuild it with:\n"
          f"    uv run -m Edge_cosine_met_reference_network.run_network "
          f"--condition {cfg.CONDITION_NAME} --viz-only")


if __name__ == "__main__":
    main()
