"""
generate_key_player_pieplots.py
================================
For each of the 5 dashboard cases, generates a figure with three pie charts
side by side (Hormones | Metabolites | Proteins/Transporters).

Each pie shows the fraction of total mentions each key player contributes,
aggregated across ALL organ-organ connections in that case:

  • Cosine cases (met_healthy/obese, general_healthy/obese):
      reads bootstrap_results_*.json — each pair entry has
      key_players_bootstrap.{hormones,metabolites,proteins}, each item a
      {term, mean, freq} dict; `mean` (bootstrap-estimated co-occurrence
      count) is summed across all pairs per term.

  • reference_metabolic:
      reads metabolic_literature_results.json — each pair entry has
      key_players.{hormones,metabolites,proteins}_counts dicts
      (term → raw count); raw counts are summed across all pairs per term.

Top 10 key players per category get individually named and coloured slices;
everything below rank 10 is merged into a single grey "Other" slice.

Run from the project root (Metabolic_Reference_Network/ folder):
    python generate_key_player_pieplots.py
Output: visualizations/<case>/key_players_pie.png
"""

import json
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE    = Path(__file__).resolve().parent
OUT_DIR = HERE / "visualizations"
DPI     = 600
# The combined 3x3 overview grid is physically much larger (22x27 in) than
# the standalone per-case pie figures, so a naive DPI bump balloons file
# size fast: 900 DPI here produced a ~470 megapixel file that PIL (and very
# likely PowerPoint, Word, and most browsers) refuse to open by default.
# 400 DPI keeps it to ~95 megapixels -- comfortably under common
# image-library safety limits (~180 megapixels).
GRID_DPI = 400
TOP_N   = 15

_TAB20      = matplotlib.colormaps["tab20"]
OTHER_COLOR = "#64748b"

CATEGORY_LABELS = {
    "hormones":    "Hormones",
    "metabolites": "Metabolites",
    "proteins":    "Proteins / Transporters",
}


# ── Data loaders ──────────────────────────────────────────────────────────────

def _aggregate_bootstrap(bootstrap_path: Path) -> dict[str, dict[str, float]]:
    """
    Cosine-pipeline bootstrap_results JSON → {category: {term: summed_mean}}.
    Each pair entry has key_players_bootstrap with lists of {term, mean, freq};
    we sum `mean` across all pairs per term.
    """
    data = json.loads(bootstrap_path.read_text(encoding="utf-8"))
    totals: dict[str, dict[str, float]] = {
        "hormones":    defaultdict(float),
        "metabolites": defaultdict(float),
        "proteins":    defaultdict(float),
    }
    for entry in data.values():
        kp = entry.get("key_players_bootstrap", {})
        for cat in ("hormones", "metabolites", "proteins"):
            for item in kp.get(cat, []):
                term = (item.get("term") or "").strip()
                mean = float(item.get("mean") or 0.0)
                if term:
                    totals[cat][term] += mean
    return {cat: dict(v) for cat, v in totals.items()}


def _aggregate_reference(results_path: Path) -> dict[str, dict[str, float]]:
    """
    Reference-network metabolic_literature_results JSON →
    {category: {term: summed_raw_count}}.
    Each pair entry has key_players with *_counts dicts (term → raw count).
    Aggregates every pair present in the file -- for the healthy/obese
    scopes, pass the condition-specific results file (already restricted to
    that cohort's pairs and that condition's filtered papers) rather than
    filtering the "all" file down to a pair subset.
    """
    data = json.loads(results_path.read_text(encoding="utf-8"))
    totals: dict[str, dict[str, float]] = {
        "hormones":    defaultdict(float),
        "metabolites": defaultdict(float),
        "proteins":    defaultdict(float),
    }
    cat_keys = (
        ("hormones",    "hormones_counts"),
        ("metabolites", "metabolites_counts"),
        ("proteins",    "proteins_counts"),
    )
    for entry in data.values():
        kp = entry.get("key_players", {})
        for cat, key in cat_keys:
            for term, count in (kp.get(key) or {}).items():
                term = term.strip()
                if term:
                    totals[cat][term] += float(count)
    return {cat: dict(v) for cat, v in totals.items()}


# ── Case registry ─────────────────────────────────────────────────────────────

CASES = {
    "met_healthy": {
        "title":  "Metabolic Network — Healthy",
        "loader": lambda: _aggregate_bootstrap(
            HERE / "Edge_cosine_met_reference_network" / "healthy"
                 / "bootstrap_results_healthy.json"
        ),
    },
    "met_obese": {
        "title":  "Metabolic Network — Obese",
        "loader": lambda: _aggregate_bootstrap(
            HERE / "Edge_cosine_met_reference_network" / "obese"
                 / "bootstrap_results_obese.json"
        ),
    },
    "general_healthy": {
        "title":  "Metabolic + Hormonal Network — Healthy",
        "loader": lambda: _aggregate_bootstrap(
            HERE / "Edge_cosine_general_reference_network" / "healthy"
                 / "bootstrap_results_healthy.json"
        ),
    },
    "general_obese": {
        "title":  "Metabolic + Hormonal Network — Obese",
        "loader": lambda: _aggregate_bootstrap(
            HERE / "Edge_cosine_general_reference_network" / "obese"
                 / "bootstrap_results_obese.json"
        ),
    },
    "reference_metabolic": {
        "title":  "Metabolic Reference Network — All Connections",
        "loader": lambda: _aggregate_reference(
            HERE / "reference_network_only_metabolic"
                 / "metabolic_literature_results.json"
        ),
    },
    # healthy/obese use their own condition-filtered results file (see
    # run_metabolic_lit_search.py --condition healthy|obese) -- papers
    # already restricted to that condition's keyword filter, key players
    # recomputed from that smaller set -- not the "all" results filtered
    # down to that scope's pairs.
    "reference_metabolic_healthy": {
        "title":  "Metabolic Reference Network — Healthy Connections",
        "loader": lambda: _aggregate_reference(
            HERE / "reference_network_only_metabolic"
                 / "metabolic_literature_results_healthy.json"
        ),
    },
    "reference_metabolic_obese": {
        "title":  "Metabolic Reference Network — Obese Connections",
        "loader": lambda: _aggregate_reference(
            HERE / "reference_network_only_metabolic"
                 / "metabolic_literature_results_obese.json"
        ),
    },
}


# ── Plot helpers ──────────────────────────────────────────────────────────────

def _top_n_plus_other(
    counts: dict[str, float], n: int = TOP_N
) -> tuple[list[str], list[float]]:
    """Sort by descending count, name the top n, collapse the rest to 'Other'."""
    sorted_items = sorted(counts.items(), key=lambda x: -x[1])
    top   = sorted_items[:n]
    other = sorted_items[n:]
    labels = [t for t, _ in top]
    values = [v for _, v in top]
    if other:
        labels.append("Other")
        values.append(sum(v for _, v in other))
    return labels, values


def _draw_pie_ax(ax: plt.Axes, counts: dict[str, float], category: str) -> None:
    ax.set_facecolor("white")
    title = CATEGORY_LABELS[category]

    if not counts or sum(counts.values()) == 0:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, color="#94a3b8", fontsize=13)
        ax.set_title(title, fontsize=14, pad=12, color="#1e293b", fontweight="bold")
        ax.axis("off")
        return

    labels, values = _top_n_plus_other(counts)
    total = sum(values)

    has_other = labels[-1] == "Other"
    n_named   = len(labels) - (1 if has_other else 0)
    colors    = [_TAB20(i / max(n_named - 1, 1)) for i in range(n_named)]
    if has_other:
        colors.append(OTHER_COLOR)

    wedges, _ = ax.pie(
        values,
        colors=colors,
        startangle=90,
        wedgeprops={"linewidth": 0.8, "edgecolor": "white"},
        counterclock=False,
    )

    legend_labels = [
        f"{lbl}  ({v / total * 100:.1f}%)"
        for lbl, v in zip(labels, values)
    ]
    ax.legend(
        wedges,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.04),
        ncol=2,
        fontsize=9,
        frameon=True,
        framealpha=0.95,
        edgecolor="#e2e8f0",
        handlelength=1.2,
        handletextpad=0.5,
        columnspacing=1.0,
    )
    ax.set_title(title, fontsize=14, pad=12, color="#1e293b", fontweight="bold")


def draw_key_player_grid(
    scope_rows: list[tuple[str, dict[str, dict[str, float]]]], out_path: Path
) -> None:
    """scope_rows: ordered list of (row_label, kp_data), one row per
    condition (top to bottom). Columns are Hormones | Metabolites |
    Proteins/Transporters, reusing _draw_pie_ax so every pie renders
    identically to its standalone counterpart."""
    n_rows = len(scope_rows)
    fig, axes = plt.subplots(n_rows, 3, figsize=(22, 9 * n_rows), dpi=GRID_DPI)
    fig.patch.set_facecolor("white")
    if n_rows == 1:
        axes = axes.reshape(1, 3)

    for row_idx, (row_label, kp_data) in enumerate(scope_rows):
        for col_idx, cat in enumerate(("hormones", "metabolites", "proteins")):
            _draw_pie_ax(axes[row_idx, col_idx], kp_data.get(cat, {}), cat)

        # Row label to the left of the first subplot in this row.
        axes[row_idx, 0].text(
            -0.18, 0.5, row_label, transform=axes[row_idx, 0].transAxes,
            rotation=90, ha="center", va="center", fontsize=15, fontweight="bold",
            color="#1e293b",
        )

    fig.suptitle("Key Player Distribution by Condition", fontsize=20, y=1.0, color="#1e293b")
    fig.tight_layout(rect=[0.015, 0, 1, 0.98])

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=GRID_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  [ok] {out_path.relative_to(HERE)}")


def draw_key_player_figure(
    kp_data: dict[str, dict[str, float]], title: str, out_path: Path
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(22, 9))
    fig.patch.set_facecolor("white")

    for ax, cat in zip(axes, ("hormones", "metabolites", "proteins")):
        _draw_pie_ax(ax, kp_data.get(cat, {}), cat)

    fig.suptitle(
        f"{title}\nKey Player Distribution Across All Organ–Organ Connections",
        fontsize=15,
        y=1.01,
        color="#1e293b",
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  [ok] {out_path.relative_to(HERE)}")


# ── Main ──────────────────────────────────────────────────────────────────────

# case_name -> row label, for the combined reference-network grid figure
# (all/healthy/obese only -- the cosine-pipeline cases aren't part of it).
_REFERENCE_GRID_ROWS = {
    "reference_metabolic":         "All Connections",
    "reference_metabolic_healthy": "Healthy Connections",
    "reference_metabolic_obese":   "Obese Connections",
}


def main() -> None:
    grid_data: dict[str, dict[str, dict[str, float]]] = {}

    for case_name, cfg in CASES.items():
        print(f"[i] {case_name}")
        try:
            kp_data = cfg["loader"]()
        except FileNotFoundError as exc:
            print(f"  [!] skipped: {exc}")
            continue

        out_path = OUT_DIR / case_name / "key_players_pie.png"
        draw_key_player_figure(kp_data, cfg["title"], out_path)

        if case_name in _REFERENCE_GRID_ROWS:
            grid_data[case_name] = kp_data

    print("\n[ok] All key player pie plots done.")

    if grid_data:
        print("[i] reference_metabolic_key_players_grid")
        grid_rows = [
            (_REFERENCE_GRID_ROWS[name], grid_data[name])
            for name in ("reference_metabolic", "reference_metabolic_healthy", "reference_metabolic_obese")
            if name in grid_data
        ]
        draw_key_player_grid(grid_rows, OUT_DIR / "reference_metabolic_key_players_grid.png")


if __name__ == "__main__":
    main()
