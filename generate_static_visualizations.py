"""
Static, high-resolution PNG exports of the organ cross-talk chord (Sankey)
diagram and a matching papers-per-organ bar chart, for each of the 5
dashboard cases.

Reuses the exact same organ_papers/links computation as each dashboard's
in-browser "Literature Statistics" tab (Edge_cosine_met_reference_network,
Edge_cosine_general_reference_network, reference_network_only_metabolic),
and the same ORGAN_COLORS palette, so the static exports visually match the
interactive versions.

Run from the project root:
    python generate_static_visualizations.py
Outputs to visualizations/<case>/sankey.png and papers_per_organ.png
"""

import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Patch

from Visualisation.networkBuilderUtils import ORGAN_COLORS, DEFAULT_NODE_COLOR

HERE    = Path(__file__).resolve().parent
OUT_DIR = HERE / "visualizations"
DPI     = 600
# The combined multi-row/column overview grids are physically much larger
# (up to ~22x30 in) than the standalone per-case plots, so a naive DPI bump
# balloons file size fast: 900 DPI here produced 300-470 megapixel files
# that PIL (and very likely PowerPoint, Word, and most browsers) refuse to
# open by default. 400 DPI keeps the largest grid to ~95-120 megapixels --
# comfortably under common image-library safety limits (~180 megapixels)
# while still ~1.3x sharper than these grids' original 300 DPI.
GRID_DPI = 400

sys.path.insert(0, str(HERE / "reference_network_only_metabolic"))
import config as _ref_cfg  # noqa: E402 -- needed for CONNECTION_TYPES labels


# ── Data loaders ──────────────────────────────────────────────────────────

def _load_edge_filter(csv_path: Path) -> list[tuple[str, str]]:
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
                if val.strip() == "1" and col_idx < len(col_organs):
                    col_organ = col_organs[col_idx].strip()
                    if row_organ and col_organ and row_organ != col_organ:
                        pairs.add((min(row_organ, col_organ), max(row_organ, col_organ)))
    return sorted(pairs)


def _load_cosine_case(folder: str, condition: str):
    """met/general pipelines: per-organ search + bootstrap co-occurrence."""
    base       = HERE / folder / condition
    search_fp  = base / f"search_results_{condition}.json"
    boot_fp    = base / f"bootstrap_results_{condition}.json"
    if not search_fp.exists() or not boot_fp.exists():
        missing = search_fp.name if not search_fp.exists() else boot_fp.name
        raise FileNotFoundError(f"{missing} not found in {base} — pipeline hasn't finished running yet")

    search = json.loads(search_fp.read_text(encoding="utf-8"))
    boot   = json.loads(boot_fp.read_text(encoding="utf-8"))

    organs = sorted(search.keys())
    organ_papers = {o: search.get(o, {}).get("n_found", 0) for o in organs}

    links = []
    for b in boot.values():
        o1, o2 = b.get("organ1"), b.get("organ2")
        n = b.get("n_cooccur_total", 0)
        if o1 in organ_papers and o2 in organ_papers and n > 0:
            if organ_papers[o1] < organ_papers[o2]:
                o1, o2 = o2, o1
            links.append({"source": o1, "target": o2, "value": n})
    return organ_papers, links


def _load_reference_case(results_filename: str = "metabolic_literature_results.json",
                          csv_filename: str = "all_organ_connections.csv"):
    """reference_network_only_metabolic: per-pair search, organ totals are
    the union of PMIDs across all pairs touching that organ.
    healthy/obese pass their own condition-filtered results file + cohort
    CSV (see run_metabolic_lit_search.py --condition healthy|obese) instead
    of the "all" results filtered down to that scope's pairs."""
    base = HERE / "reference_network_only_metabolic"
    results_fp = base / results_filename
    csv_fp     = base / csv_filename
    if not results_fp.exists() or not csv_fp.exists():
        raise FileNotFoundError(f"required files not found in {base}")

    results = json.loads(results_fp.read_text(encoding="utf-8"))
    pairs   = _load_edge_filter(csv_fp)

    results_by_pair = {}
    for v in results.values():
        o1, o2 = v.get("organ1", ""), v.get("organ2", "")
        if o1 and o2:
            results_by_pair[(min(o1, o2), max(o1, o2))] = v

    organs = sorted({o for pair in pairs for o in pair})
    organ_pmids = {o: set() for o in organs}
    for pair in pairs:
        data = results_by_pair.get(pair, {})
        if data.get("n_papers_found", 0) <= 0:
            continue
        pmids = {p.get("pmid") for p in data.get("papers", []) if p.get("pmid")}
        organ_pmids[pair[0]].update(pmids)
        organ_pmids[pair[1]].update(pmids)
    organ_papers = {o: len(organ_pmids.get(o, set())) for o in organs}

    links = []
    for pair in pairs:
        o1, o2 = pair
        data = results_by_pair.get(pair, {})
        n = data.get("n_papers_found", 0)
        if n <= 0:
            continue
        if organ_papers.get(o1, 0) < organ_papers.get(o2, 0):
            o1, o2 = o2, o1
        links.append({"source": o1, "target": o2, "value": n})
    return organ_papers, links


_COSINE_SUBTITLE = "Total Papers Found per Organ (direct per-organ search, same-sentence organ + crosstalk-keyword co-mentions)"
_REFERENCE_SUBTITLE = "Total Papers Found per Organ (union of papers across all its predefined organ-organ axes)"

CASES = {
    "met_healthy":         (lambda: _load_cosine_case("Edge_cosine_met_reference_network", "healthy"),
                             "Metabolic Network — Healthy", _COSINE_SUBTITLE),
    "met_obese":           (lambda: _load_cosine_case("Edge_cosine_met_reference_network", "obese"),
                             "Metabolic Network — Obese", _COSINE_SUBTITLE),
    "general_healthy":     (lambda: _load_cosine_case("Edge_cosine_general_reference_network", "healthy"),
                             "Metabolic + Hormonal Network — Healthy", _COSINE_SUBTITLE),
    "general_obese":       (lambda: _load_cosine_case("Edge_cosine_general_reference_network", "obese"),
                             "Metabolic + Hormonal Network — Obese", _COSINE_SUBTITLE),
    "reference_metabolic": (_load_reference_case,
                             "Metabolic Reference Network — All Connections", _REFERENCE_SUBTITLE),
    # healthy/obese use their own condition-filtered results file + cohort
    # CSV -- not the "all" results filtered down to that scope's pairs.
    "reference_metabolic_healthy": (
        lambda: _load_reference_case("metabolic_literature_results_healthy.json",
                                      "healthy_cohort_connections.csv"),
        "Metabolic Reference Network — Healthy Connections", _REFERENCE_SUBTITLE),
    "reference_metabolic_obese": (
        lambda: _load_reference_case("metabolic_literature_results_obese.json",
                                      "obese_cohort_connections.csv"),
        "Metabolic Reference Network — Obese Connections", _REFERENCE_SUBTITLE),
}

# Row order + label for the combined reference-network overview grid
# (draw_reference_overview_grid) -- all/healthy/obese only.
_REFERENCE_OVERVIEW_ROWS = {
    "reference_metabolic":         "All Connections",
    "reference_metabolic_healthy": "Healthy Connections",
    "reference_metabolic_obese":   "Obese Connections",
}


# ── Connection-type data (reference_network_only_metabolic only) ───────────
# LLM-classified connection type per organ pair (metabolic_connection_types.json,
# generated by run_metabolic_lit_search.py / Literature_Search/llm_connection_type.py).
# A pair may have 0-3 types (majority vote across 3 LLM passes); a pair with
# zero types is bucketed under UNCLASSIFIED_KEY so bar totals stay honest
# about every pair in scope, not just the classified ones.

UNCLASSIFIED_KEY   = "__unclassified__"
UNCLASSIFIED_LABEL = "No type assigned"

CONNECTION_TYPE_LABELS = {k: v["label"] for k, v in _ref_cfg.CONNECTION_TYPES.items()}
CONNECTION_TYPE_LABELS[UNCLASSIFIED_KEY] = UNCLASSIFIED_LABEL

# Which edge-filter CSV + LLM classification file defines each of the three
# scopes to plot. healthy/obese use their own condition-filtered
# classification (run_metabolic_lit_search.py --condition healthy|obese) --
# not the "all" classification filtered down to that scope's pairs.
CONNECTION_TYPE_SCOPES = {
    "all":     ("all_organ_connections.csv",      "All Connections",
                "metabolic_connection_types.json"),
    "healthy": ("healthy_cohort_connections.csv",  "Healthy Connections",
                "metabolic_connection_types_healthy.json"),
    "obese":   ("obese_cohort_connections.csv",    "Obese Connections",
                "metabolic_connection_types_obese.json"),
}


def _load_connection_types(
    filename: str = "metabolic_connection_types.json",
) -> dict[tuple[str, str], list[str]]:
    """{(organ1, organ2) sorted: [type_key, ...]} for every classified pair
    in the given classification file (defaults to the "all" case)."""
    types_fp = HERE / "reference_network_only_metabolic" / filename
    if not types_fp.exists():
        raise FileNotFoundError(f"{types_fp.name} not found — run the matching LLM "
                                 f"connection-type step first (uv run python "
                                 f"run_metabolic_lit_search.py [--condition healthy|obese])")
    raw = json.loads(types_fp.read_text(encoding="utf-8"))
    pair_types = {}
    for key, entry in raw.items():
        o1, o2 = key.split("|", 1)
        pair_types[(min(o1, o2), max(o1, o2))] = entry.get("types", [])
    return pair_types


def _connection_type_counts(pair_types: dict, pairs: list[tuple[str, str]]):
    """Tally type assignments for one scope's pairs.

    Returns (type_count, organ_type_count):
      type_count        -- {type_key: n} across all `pairs`
      organ_type_count  -- {organ: {type_key: n}}, each pair contributing to
                            both of its organs (a connection touches both)
    A pair with no assigned types counts once toward UNCLASSIFIED_KEY.
    """
    type_count = {k: 0 for k in CONNECTION_TYPE_LABELS}
    organ_type_count: dict[str, dict[str, int]] = {}

    for o1, o2 in pairs:
        types = pair_types.get((o1, o2)) or [UNCLASSIFIED_KEY]
        for t in types:
            type_count[t] = type_count.get(t, 0) + 1
            for o in (o1, o2):
                od = organ_type_count.setdefault(o, {})
                od[t] = od.get(t, 0) + 1
    return type_count, organ_type_count


def _type_colors(type_keys_in_order: list[str]) -> dict[str, tuple]:
    """One distinct color per type key (tab20 palette) + fixed gray for
    UNCLASSIFIED_KEY, so 'unclassified' reads as neutral in every plot."""
    palette = plt.cm.tab20.colors
    colors = {k: palette[i % len(palette)] for i, k in enumerate(type_keys_in_order)}
    colors[UNCLASSIFIED_KEY] = (0.72, 0.72, 0.76, 1.0)
    return colors


# ── Connection-type bar plot (one bar per type, this scope's total count) ──
# Core drawing logic lives in _draw_* helpers that take an existing Axes, so
# both the standalone single-scope PNGs and the combined 3-row grid (see
# draw_connection_type_grid) render identical bars from one code path.

def _draw_connection_type_bars(ax, type_count: dict, order: list[str], colors: dict,
                                label_fontsize: int = 10):
    values = [type_count.get(k, 0) for k in order]
    names  = [CONNECTION_TYPE_LABELS.get(k, k) for k in order]
    bar_colors = [colors[k] for k in order]

    bars = ax.barh(names, values, color=bar_colors, edgecolor="white", linewidth=0.6)
    ax.invert_yaxis()
    ax.set_xlabel("Times assigned", fontsize=label_fontsize + 1, color="#1e293b")
    ax.tick_params(colors="#1e293b", labelsize=label_fontsize)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#cbd5e1")
    ax.grid(axis="x", color="#e2e8f0", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)

    max_val = max(values) if values else 1
    for bar, v in zip(bars, values):
        ax.text(bar.get_width() + max_val * 0.01, bar.get_y() + bar.get_height() / 2,
                f"{v:,}", va="center", fontsize=label_fontsize, color="#1e293b")


def draw_connection_type_barplot(type_count: dict, order: list[str], colors: dict,
                                  title: str, subtitle: str, out_path: Path):
    fig, ax = plt.subplots(figsize=(11, max(5, 0.5 * len(order))), dpi=DPI)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    _draw_connection_type_bars(ax, type_count, order, colors)
    ax.set_title(f"{title}\n{subtitle}", fontsize=15, pad=14, color="#1e293b")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  [ok] {out_path.relative_to(HERE)}")


# ── Per-organ stacked connection-type bar plot ──────────────────────────────

def _draw_organ_stacked_bars(ax, organ_type_count: dict, order: list[str], colors: dict,
                              show_legend: bool = False, label_fontsize: int = 10):
    organ_totals = {o: sum(d.values()) for o, d in organ_type_count.items()}
    organs = sorted(organ_totals, key=lambda o: -organ_totals[o])

    left = np.zeros(len(organs))
    for k in order:
        vals = np.array([organ_type_count.get(o, {}).get(k, 0) for o in organs], dtype=float)
        if vals.sum() == 0:
            continue
        ax.barh(organs, vals, left=left, color=colors[k], edgecolor="white",
                linewidth=0.4, label=CONNECTION_TYPE_LABELS.get(k, k))
        left += vals

    ax.invert_yaxis()
    ax.set_xlabel("Connections assigned", fontsize=label_fontsize + 1, color="#1e293b")
    ax.tick_params(colors="#1e293b", labelsize=label_fontsize)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#cbd5e1")
    ax.grid(axis="x", color="#e2e8f0", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)

    max_total = max(organ_totals.values(), default=1) or 1
    for i, o in enumerate(organs):
        ax.text(organ_totals[o] + max_total * 0.01, i, f"{organ_totals[o]:,}",
                va="center", fontsize=label_fontsize, color="#1e293b")

    if show_legend:
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9, frameon=False)


def draw_organ_connection_type_stacked_barplot(organ_type_count: dict, order: list[str],
                                                colors: dict, title: str, subtitle: str,
                                                out_path: Path):
    if not organ_type_count:
        print(f"  [!] no organ connection-type data — skipping {out_path.name}")
        return

    n_organs = len(organ_type_count)
    fig, ax = plt.subplots(figsize=(13, max(5, 0.5 * n_organs)), dpi=DPI)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    _draw_organ_stacked_bars(ax, organ_type_count, order, colors, show_legend=True)
    ax.set_title(f"{title}\n{subtitle}", fontsize=15, pad=14, color="#1e293b")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  [ok] {out_path.relative_to(HERE)}")


# ── Combined 3-row grid: All / Healthy / Obese, each with both plot types ──

def draw_connection_type_grid(scope_rows: list[tuple[str, dict, dict]], order: list[str],
                               colors: dict, out_path: Path):
    """scope_rows: ordered list of (row_label, type_count, organ_type_count),
    one row per scope (top to bottom). Left column = frequency barplot,
    right column = per-organ stacked barplot, sharing one legend."""
    n_rows = len(scope_rows)
    max_organs = max((len(oc) for _, _, oc in scope_rows), default=1)
    row_h = max(4.2, 0.4 * max_organs)

    fig, axes = plt.subplots(n_rows, 2, figsize=(19, row_h * n_rows), dpi=GRID_DPI,
                              gridspec_kw={"width_ratios": [1, 1.3]})
    fig.patch.set_facecolor("white")
    if n_rows == 1:
        axes = axes.reshape(1, 2)

    for row_idx, (row_label, type_count, organ_type_count) in enumerate(scope_rows):
        ax_freq, ax_organ = axes[row_idx]
        ax_freq.set_facecolor("white")
        ax_organ.set_facecolor("white")

        _draw_connection_type_bars(ax_freq, type_count, order, colors, label_fontsize=9)
        ax_freq.set_title(f"{row_label}, Type frequency", fontsize=13, pad=8, color="#1e293b")

        _draw_organ_stacked_bars(ax_organ, organ_type_count, order, colors,
                                  show_legend=False, label_fontsize=9)
        ax_organ.set_title(f"{row_label}, Per organ", fontsize=13, pad=8, color="#1e293b")

    legend_handles = [Patch(facecolor=colors[k], edgecolor="white", label=CONNECTION_TYPE_LABELS.get(k, k))
                       for k in order]
    fig.legend(handles=legend_handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
               fontsize=10, frameon=False, title="Connection type", title_fontsize=11)

    fig.suptitle("Connection Types by Condition",
                 fontsize=18, y=0.995, color="#1e293b")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.subplots_adjust(hspace=0.35)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=GRID_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  [ok] {out_path.relative_to(HERE)}")


# ── Chord / Sankey diagram ───────────────────────────────────────────────

def _hex_to_rgb(h: str):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


def _blend(c1: str, c2: str, t: float = 0.5):
    a, b = _hex_to_rgb(c1), _hex_to_rgb(c2)
    return tuple(a[i] + (b[i] - a[i]) * t for i in range(3))


def _draw_chord_diagram_ax(ax, organ_papers: dict, links: list) -> bool:
    """Core chord/sankey drawing onto an existing Axes (no title, no figure
    creation/save). Returns False (and draws nothing) if there are no
    cross-mention links -- caller decides whether to skip entirely
    (standalone export) or show a placeholder (grid)."""
    organs_with_links = sorted(
        {o for l in links for o in (l["source"], l["target"])},
        key=lambda o: -organ_papers.get(o, 0),
    )
    if not organs_with_links or not links:
        return False

    totals = {o: 0 for o in organs_with_links}
    for l in links:
        totals[l["source"]] += l["value"]
        totals[l["target"]] += l["value"]
    grand = sum(totals.values()) or 1

    GAP_DEG = 1.4 if len(organs_with_links) > 1 else 0
    usable_deg = 360 - GAP_DEG * len(organs_with_links)
    R = 1.0
    ARC_W = 16   # linewidth in points

    def color_for(o):
        return ORGAN_COLORS.get(o, DEFAULT_NODE_COLOR)

    # deg follows the canvas convention used by the HTML dashboard (-90 =
    # top, increasing clockwise). The y-component is negated so the shape
    # renders correctly in matplotlib's y-up coordinate system directly,
    # without needing an axis-inversion hack that would desync label
    # rotation (a screen-space angle, unaffected by axis inversion) from
    # element position (a data-space angle, which IS affected by it).
    def point(deg, radius):
        rad = np.radians(deg)
        return radius * np.cos(rad), -radius * np.sin(rad)

    angle = -90.0
    organ_range = {}
    for o in organs_with_links:
        span = (totals[o] / grand) * usable_deg
        organ_range[o] = (angle, angle + span)
        angle += span + GAP_DEG

    # Ribbons (drawn first, under the arcs)
    cursor = {o: organ_range[o][0] for o in organs_with_links}
    for l in sorted(links, key=lambda l: l["source"] + l["target"]):
        span = (l["value"] / grand) * usable_deg
        a0 = cursor[l["source"]]; cursor[l["source"]] += span
        b0 = cursor[l["target"]]; cursor[l["target"]] += span
        a1, b1 = a0 + span, b0 + span

        x0, y0   = point(a0, R)
        x0b, y0b = point(a1, R)
        x1, y1   = point(b0, R)
        x1b, y1b = point(b1, R)

        verts = [(x0, y0), (0, 0), (x1, y1), (x1b, y1b), (0, 0), (x0b, y0b), (x0, y0)]
        codes = [MplPath.MOVETO, MplPath.CURVE3, MplPath.CURVE3,
                 MplPath.LINETO, MplPath.CURVE3, MplPath.CURVE3, MplPath.CLOSEPOLY]
        patch = PathPatch(MplPath(verts, codes),
                          facecolor=_blend(color_for(l["source"]), color_for(l["target"])),
                          edgecolor="none", alpha=0.55, zorder=1)
        ax.add_patch(patch)

    # Arcs
    for o in organs_with_links:
        start, end = organ_range[o]
        if end <= start:
            continue
        theta = np.radians(np.linspace(start, end, 60))
        ax.plot(R * np.cos(theta), -R * np.sin(theta),
                color=color_for(o), linewidth=ARC_W, solid_capstyle="butt", zorder=2)

    # Labels: organs with a small paper share get thin arc slices that can
    # sit only a couple of degrees apart, so placing each label at its own
    # slice's exact midpoint (as the HTML/JS canvas version does) makes
    # neighbouring labels collide once rendered at print resolution. Keep
    # each label's natural angle as a starting point, then push apart any
    # that are closer than MIN_GAP_DEG so text never overlaps; a short
    # leader line reconnects a nudged label back to its actual arc.
    FONT_SIZE   = 11
    LABEL_R     = R + 0.22
    MIN_GAP_DEG = 5.5

    drawn = [o for o in organs_with_links if organ_range[o][1] > organ_range[o][0]]
    raw_mids = [(organ_range[o][0] + organ_range[o][1]) / 2 for o in drawn]

    label_mids = list(raw_mids)
    for _ in range(200):
        moved = False
        for i in range(1, len(label_mids)):
            gap = label_mids[i] - label_mids[i - 1]
            if gap < MIN_GAP_DEG:
                push = (MIN_GAP_DEG - gap) / 2
                label_mids[i]     += push
                label_mids[i - 1] -= push
                moved = True
        if not moved:
            break

    for o, raw_mid, mid in zip(drawn, raw_mids, label_mids):
        if abs(mid - raw_mid) > 0.3:
            ax_x, ax_y = point(raw_mid, R + 0.015)
            lead_x, lead_y = point(mid, LABEL_R - 0.03)
            ax.plot([ax_x, lead_x], [ax_y, lead_y],
                    color="#94a3b8", linewidth=0.8, zorder=2)

        lx, ly = point(mid, LABEL_R)
        flip = 90 < (mid % 360) < 270
        # Matplotlib's `rotation` is a screen-space angle, counterclockwise
        # from horizontal — the negative of our canvas-style `deg`, since
        # point() negates y to go from canvas to plot coordinates. This
        # keeps the label pointing straight outward from its own node
        # (horizontal at the 3/9 o'clock positions, vertical at 12/6
        # o'clock, and the node's own angle everywhere in between) instead
        # of drifting out of sync once the axis is no longer flipped.
        rot = -mid + (180 if flip else 0)
        ax.text(lx, ly, o, rotation=rot, rotation_mode="anchor",
                ha="right" if flip else "left", va="center",
                fontsize=FONT_SIZE, color="#1e293b", zorder=3)

    ax.set_xlim(-1.75, 1.75)
    ax.set_ylim(-1.75, 1.75)
    ax.set_aspect("equal")
    ax.axis("off")
    return True


def draw_chord_diagram(organ_papers: dict, links: list, title: str, out_path: Path):
    fig, ax = plt.subplots(figsize=(11, 11), dpi=DPI)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    if not _draw_chord_diagram_ax(ax, organ_papers, links):
        print(f"  [!] no cross-mention links — skipping chord diagram for {out_path.name}")
        plt.close(fig)
        return

    ax.set_title(f"{title}\nOrgan Cross-Talk — Sankey Diagram", fontsize=16, pad=18, color="#1e293b")
    fig.text(0.5, 0.02, "Ribbon thickness = number of same-sentence co-occurring papers for that organ pair.",
              ha="center", fontsize=10, color="#64748b")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  [ok] {out_path.relative_to(HERE)}")


# ── Bar plot ──────────────────────────────────────────────────────────────

def _draw_barplot_ax(ax, organ_papers: dict, label_fontsize: int = 10) -> bool:
    """Core papers-per-organ bar drawing onto an existing Axes (no title,
    no figure creation/save). Returns False if there's no data."""
    if not organ_papers:
        return False

    organs = sorted(organ_papers, key=lambda o: -organ_papers[o])
    values = [organ_papers[o] for o in organs]
    colors = [ORGAN_COLORS.get(o, DEFAULT_NODE_COLOR) for o in organs]

    bars = ax.barh(organs, values, color=colors, edgecolor="white", linewidth=0.6)
    ax.invert_yaxis()
    ax.set_xlabel("Papers found", fontsize=label_fontsize + 2, color="#1e293b")
    ax.tick_params(colors="#1e293b", labelsize=label_fontsize + 1)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#cbd5e1")
    ax.grid(axis="x", color="#e2e8f0", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)

    max_val = max(values) if values else 1
    for bar, v in zip(bars, values):
        ax.text(bar.get_width() + max_val * 0.01, bar.get_y() + bar.get_height() / 2,
                f"{v:,}", va="center", fontsize=label_fontsize, color="#1e293b")
    return True


def draw_barplot(organ_papers: dict, title: str, subtitle: str, out_path: Path):
    if not organ_papers:
        print(f"  [!] no organ paper counts — skipping bar plot for {out_path.name}")
        return

    fig, ax = plt.subplots(figsize=(11, max(5, 0.5 * len(organ_papers))), dpi=DPI)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    _draw_barplot_ax(ax, organ_papers, label_fontsize=10)
    ax.set_title(f"{title}\n{subtitle}", fontsize=15, pad=14, color="#1e293b")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  [ok] {out_path.relative_to(HERE)}")


# ── Combined sankey + papers-per-organ grid (reference_network_only_metabolic) ──

def draw_reference_overview_grid(scope_rows: list[tuple[str, dict, list]], out_path: Path):
    """scope_rows: ordered list of (row_label, organ_papers, links), one row
    per condition (top to bottom). Left column = sankey diagram, right
    column = papers-per-organ bar chart, reusing the same _draw_*_ax
    helpers as the standalone exports so every panel renders identically to
    its standalone counterpart."""
    n_rows = len(scope_rows)
    max_organs = max((len(op) for _, op, _ in scope_rows), default=1)
    row_h = max(9, 0.45 * max_organs)

    fig, axes = plt.subplots(n_rows, 2, figsize=(22, row_h * n_rows), dpi=GRID_DPI,
                              gridspec_kw={"width_ratios": [1.6, 1], "hspace": 0.4,
                                           "top": 0.96, "bottom": 0.02})
    fig.patch.set_facecolor("white")
    if n_rows == 1:
        axes = axes.reshape(1, 2)

    for row_idx, (row_label, organ_papers, links) in enumerate(scope_rows):
        ax_sankey, ax_bar = axes[row_idx]
        ax_sankey.set_facecolor("white")
        ax_bar.set_facecolor("white")

        if not _draw_chord_diagram_ax(ax_sankey, organ_papers, links):
            ax_sankey.text(0.5, 0.5, "No cross-mention links", ha="center", va="center",
                            transform=ax_sankey.transAxes, color="#94a3b8", fontsize=12)
            ax_sankey.axis("off")
        else:
            # Small extra headroom above the top organ label (which sits
            # close to the y=1.75 data limit) so the title never collides
            # with it, without wasting so much vertical range that the
            # circle itself shrinks (aspect="equal" fits the larger of
            # xlim/ylim span into the available box).
            ax_sankey.set_ylim(-1.75, 2.05)
        ax_sankey.set_title(row_label, fontsize=14, pad=22, color="#1e293b")

        if not _draw_barplot_ax(ax_bar, organ_papers, label_fontsize=9):
            ax_bar.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax_bar.transAxes, color="#94a3b8", fontsize=12)
            ax_bar.axis("off")
        ax_bar.set_title(f"{row_label}, Papers per organ", fontsize=13, pad=18, color="#1e293b")

    # subplots_adjust (via gridspec_kw top/bottom/hspace above) instead of
    # tight_layout -- tight_layout's margin math is unreliable here because
    # the sankey axes use aspect="equal" (matplotlib warns about this).
    fig.suptitle("Papers per Organ by Condition", fontsize=18, y=0.985, color="#1e293b")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=GRID_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  [ok] {out_path.relative_to(HERE)}")


# ── Main ─────────────────────────────────────────────────────────────────

def main():
    reference_overview_rows: dict[str, tuple] = {}  # case_name -> (organ_papers, links)

    for case_name, (loader, title, subtitle) in CASES.items():
        print(f"[i] {case_name}")
        try:
            organ_papers, links = loader()
        except FileNotFoundError as exc:
            print(f"  [!] skipped: {exc}")
            continue

        case_dir = OUT_DIR / case_name
        draw_chord_diagram(organ_papers, links, title, case_dir / "sankey.png")
        draw_barplot(organ_papers, title, subtitle, case_dir / "papers_per_organ_barplot.png")

        if case_name in _REFERENCE_OVERVIEW_ROWS:
            reference_overview_rows[case_name] = (organ_papers, links)

    if reference_overview_rows:
        print("[i] reference_metabolic_overview_grid")
        grid_rows = [
            (_REFERENCE_OVERVIEW_ROWS[name], *reference_overview_rows[name])
            for name in ("reference_metabolic", "reference_metabolic_healthy", "reference_metabolic_obese")
            if name in reference_overview_rows
        ]
        draw_reference_overview_grid(grid_rows, OUT_DIR / "reference_metabolic_overview_grid.png")

    # Connection-type plots: how often each LLM-assigned type occurs, and
    # how those types distribute per organ, for all/healthy/obese scopes.
    # healthy/obese each use their own condition-filtered classification
    # file (see run_metabolic_lit_search.py --condition), not the "all"
    # classification filtered down to that scope's pairs.
    try:
        all_pair_types = _load_connection_types()
    except FileNotFoundError as exc:
        print(f"[!] connection-type plots skipped: {exc}")
        return

    ref_dir = HERE / "reference_network_only_metabolic"
    all_pairs = _load_edge_filter(ref_dir / "all_organ_connections.csv")
    all_type_count, _ = _connection_type_counts(all_pair_types, all_pairs)
    # Fixed category order (by "all"-scope frequency) reused across all three
    # scopes' plots so the same type sits at the same row/color in each.
    order = sorted((k for k in CONNECTION_TYPE_LABELS if k != UNCLASSIFIED_KEY),
                   key=lambda k: -all_type_count.get(k, 0))
    order.append(UNCLASSIFIED_KEY)
    colors = _type_colors(order[:-1])

    grid_rows = []  # (row_label, type_count, organ_type_count), top-to-bottom

    for scope_name, (csv_name, scope_label, types_filename) in CONNECTION_TYPE_SCOPES.items():
        print(f"[i] connection_types_{scope_name}")
        csv_fp = ref_dir / csv_name
        if not csv_fp.exists():
            print(f"  [!] skipped: {csv_fp.name} not found")
            continue
        try:
            scope_pair_types = _load_connection_types(types_filename)
        except FileNotFoundError as exc:
            print(f"  [!] skipped: {exc}")
            continue

        pairs = _load_edge_filter(csv_fp)
        type_count, organ_type_count = _connection_type_counts(scope_pair_types, pairs)
        plot_title = f"Metabolic Reference Network — {scope_label}"
        case_dir = OUT_DIR / f"connection_types_{scope_name}"

        draw_connection_type_barplot(
            type_count, order, colors, plot_title,
            f"How often each connection type was assigned ({len(pairs)} organ pairs, "
            f"LLM-classified, up to 3 types per pair)",
            case_dir / "connection_type_frequency_barplot.png")
        draw_organ_connection_type_stacked_barplot(
            organ_type_count, order, colors, plot_title,
            "Connection types touching each organ (stacked; a pair with multiple "
            "types counts toward each)",
            case_dir / "connection_type_per_organ_stacked_barplot.png")

        grid_rows.append((scope_label, type_count, organ_type_count))

    print("[i] connection_types_grid")
    draw_connection_type_grid(grid_rows, order, colors,
                               OUT_DIR / "connection_types_grid.png")


if __name__ == "__main__":
    main()
