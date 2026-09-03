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
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch

from Visualisation.networkBuilderUtils import ORGAN_COLORS, DEFAULT_NODE_COLOR

HERE    = Path(__file__).resolve().parent
OUT_DIR = HERE / "visualizations"
DPI     = 300


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


def _load_reference_case():
    """reference_network_only_metabolic: per-pair search, organ totals are
    the union of PMIDs across all pairs touching that organ."""
    base = HERE / "reference_network_only_metabolic"
    results_fp = base / "metabolic_literature_results.json"
    csv_fp     = base / "healthy_cohort_connections.csv"
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
                             "Metabolic Reference Network", _REFERENCE_SUBTITLE),
}


# ── Chord / Sankey diagram ───────────────────────────────────────────────

def _hex_to_rgb(h: str):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


def _blend(c1: str, c2: str, t: float = 0.5):
    a, b = _hex_to_rgb(c1), _hex_to_rgb(c2)
    return tuple(a[i] + (b[i] - a[i]) * t for i in range(3))


def draw_chord_diagram(organ_papers: dict, links: list, title: str, out_path: Path):
    organs_with_links = sorted(
        {o for l in links for o in (l["source"], l["target"])},
        key=lambda o: -organ_papers.get(o, 0),
    )
    if not organs_with_links or not links:
        print(f"  [!] no cross-mention links — skipping chord diagram for {out_path.name}")
        return

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

    fig, ax = plt.subplots(figsize=(11, 11), dpi=DPI)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

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
    ax.set_title(f"{title}\nOrgan Cross-Talk — Sankey Diagram", fontsize=16, pad=18, color="#1e293b")
    fig.text(0.5, 0.02, "Ribbon thickness = number of same-sentence co-occurring papers for that organ pair.",
              ha="center", fontsize=10, color="#64748b")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  [ok] {out_path.relative_to(HERE)}")


# ── Bar plot ──────────────────────────────────────────────────────────────

def draw_barplot(organ_papers: dict, title: str, subtitle: str, out_path: Path):
    if not organ_papers:
        print(f"  [!] no organ paper counts — skipping bar plot for {out_path.name}")
        return

    organs = sorted(organ_papers, key=lambda o: -organ_papers[o])
    values = [organ_papers[o] for o in organs]
    colors = [ORGAN_COLORS.get(o, DEFAULT_NODE_COLOR) for o in organs]

    fig, ax = plt.subplots(figsize=(11, max(5, 0.5 * len(organs))), dpi=DPI)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    bars = ax.barh(organs, values, color=colors, edgecolor="white", linewidth=0.6)
    ax.invert_yaxis()
    ax.set_xlabel("Papers found", fontsize=12, color="#1e293b")
    ax.set_title(f"{title}\n{subtitle}", fontsize=15, pad=14, color="#1e293b")
    ax.tick_params(colors="#1e293b", labelsize=11)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#cbd5e1")
    ax.grid(axis="x", color="#e2e8f0", linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)

    max_val = max(values) if values else 1
    for bar, v in zip(bars, values):
        ax.text(bar.get_width() + max_val * 0.01, bar.get_y() + bar.get_height() / 2,
                f"{v:,}", va="center", fontsize=10, color="#1e293b")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  [ok] {out_path.relative_to(HERE)}")


# ── Main ─────────────────────────────────────────────────────────────────

def main():
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


if __name__ == "__main__":
    main()
