"""
generate_muscle_condition_comparison.py
========================================
Illustrative hub-and-spoke figure comparing Muscle's organ-organ connections
under the Healthy vs. Obese condition filters. Muscle sits at the center of
each panel with its connected organs arranged around it; edge thickness and
label encode the number of supporting papers. Organ colors match
Visualisation.networkBuilderUtils.ORGAN_COLORS so nodes are consistent with
every other figure/dashboard in this project.

Run from the Metabolic_Reference_Network/ directory:
    uv run python generate_muscle_condition_comparison.py
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from Visualisation.networkBuilderUtils import ORGAN_COLORS, DEFAULT_NODE_COLOR

HERE = Path(__file__).parent
OUT_DIR = HERE / "visualizations"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RESULTS_HEALTHY = HERE / "reference_network_only_metabolic" / "metabolic_literature_results_healthy.json"
RESULTS_OBESE   = HERE / "reference_network_only_metabolic" / "metabolic_literature_results_obese.json"

DPI = 400
CENTER_ORGAN = "Muscle"


def _muscle_edges(results_path: Path) -> list[tuple[str, int]]:
    """Return [(other_organ, n_papers)] for every predefined Muscle edge in
    this condition's cohort, including ones with 0 papers -- the cohort edge
    exists (it's part of the partial-correlation network) even where no
    condition-specific literature evidence was found for it. Sorted by
    descending paper count."""
    data = json.load(open(results_path, encoding="utf-8"))
    edges = []
    for edge in data.values():
        o1, o2 = edge.get("organ1", ""), edge.get("organ2", "")
        if CENTER_ORGAN not in (o1, o2):
            continue
        n = edge.get("n_papers_found", 0)
        other = o2 if o1 == CENTER_ORGAN else o1
        edges.append((other, n))
    edges.sort(key=lambda x: -x[1])
    return edges


def _draw_hub_ax(ax, edges: list[tuple[str, int]], title: str, subtitle: str):
    center_color = ORGAN_COLORS.get(CENTER_ORGAN, DEFAULT_NODE_COLOR)
    n = len(edges)
    radius = 1.0

    nonzero = [c for _, c in edges if c > 0]
    max_papers = max(nonzero, default=1)
    # sqrt scaling keeps small-count edges visible next to much larger ones
    lw = lambda c: 1.5 + 9.0 * np.sqrt(c / max_papers)
    node_r = lambda c: 0.16 + 0.14 * np.sqrt(c / max_papers)

    angles = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, n, endpoint=False)
    positions = {org: (radius * np.cos(a), radius * np.sin(a)) for (org, _), a in zip(edges, angles)}

    # edges (drawn first, under the nodes)
    for org, count in edges:
        x, y = positions[org]
        if count == 0:
            # predefined cohort pair with no supporting literature found --
            # draw as a muted dashed line, not a colored/weighted real edge
            ax.plot([0, x], [0, y], color="#9ca3af", linewidth=1.4,
                    alpha=0.6, linestyle=(0, (4, 3)), solid_capstyle="round",
                    zorder=1)
            mx, my = x * 0.56, y * 0.56
            ax.text(mx, my, "0", fontsize=9.5, fontweight="bold", style="italic",
                    ha="center", va="center", color="#6b7280", zorder=3,
                    bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                              edgecolor="none", alpha=0.85))
            continue
        color = ORGAN_COLORS.get(org, DEFAULT_NODE_COLOR)
        ax.plot([0, x], [0, y], color=color, linewidth=lw(count),
                alpha=0.55, solid_capstyle="round", zorder=1)
        mx, my = x * 0.56, y * 0.56
        ax.text(mx, my, str(count), fontsize=10, fontweight="bold",
                ha="center", va="center", color="#1f2937", zorder=3,
                bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                          edgecolor="none", alpha=0.85))

    # spoke organ nodes
    for org, count in edges:
        x, y = positions[org]
        color = ORGAN_COLORS.get(org, DEFAULT_NODE_COLOR)
        r = node_r(count) if count > 0 else 0.14
        face_alpha = 1.0 if count > 0 else 0.35
        ax.add_patch(plt.Circle((x, y), r, facecolor=color, edgecolor="white",
                                 linewidth=1.8, alpha=face_alpha, zorder=4))
        ang = np.arctan2(y, x)
        label_r = radius + r + 0.28
        lx, ly = label_r * np.cos(ang), label_r * np.sin(ang)
        # align text away from the node instead of centering on top of it
        if np.cos(ang) > 0.35:
            ha = "left"
        elif np.cos(ang) < -0.35:
            ha = "right"
        else:
            ha = "center"
        va = "bottom" if np.sin(ang) > 0.35 else ("top" if np.sin(ang) < -0.35 else "center")
        label_color = color if count > 0 else "#9ca3af"
        label_style = "normal" if count > 0 else "italic"
        ax.text(lx, ly, org, fontsize=10.5, fontweight="bold", ha=ha,
                va=va, color=label_color, style=label_style, zorder=5)

    # center node
    ax.add_patch(plt.Circle((0, 0), 0.30, facecolor=center_color,
                             edgecolor="white", linewidth=2.2, zorder=6))
    ax.text(0, 0, CENTER_ORGAN, fontsize=11.5, fontweight="bold", ha="center",
            va="center", color="white", zorder=7)

    ax.set_title(title, fontsize=14, fontweight="bold", pad=14)
    ax.text(0, -2.05, subtitle, fontsize=9.5, ha="center", va="center",
            color="#6b7280", style="italic")

    lim = 2.35
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    healthy_edges = _muscle_edges(RESULTS_HEALTHY)
    obese_edges = _muscle_edges(RESULTS_OBESE)

    def subtitle(edges):
        n_evidence = sum(1 for _, c in edges if c > 0)
        total = sum(c for _, c in edges)
        return f"{n_evidence}/{len(edges)} cohort pairs with evidence, {total} papers"

    fig, axes = plt.subplots(1, 2, figsize=(11, 6.2), dpi=DPI)
    _draw_hub_ax(axes[0], healthy_edges, "Healthy", subtitle(healthy_edges))
    _draw_hub_ax(axes[1], obese_edges, "Obese", subtitle(obese_edges))

    fig.suptitle("Muscle organ-organ connections: Healthy vs. Obese",
                  fontsize=16, fontweight="bold", y=0.99)
    fig.text(0.5, 0.02,
              "Edge thickness and label = number of supporting papers "
              "(same-sentence + multi-organ filtered, condition-keyword subselected)",
              fontsize=8.5, ha="center", color="#6b7280")

    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    out_path = OUT_DIR / "muscle_healthy_vs_obese.png"
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] {out_path}")


if __name__ == "__main__":
    main()
