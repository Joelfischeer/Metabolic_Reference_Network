"""
generate_overview_statistics.py
================================
Produces three output files:

  visualizations/overview_statistics.png
      4-panel overview for the 14-organ Reference Metabolic Network only:
        A. Edge counts — No Condition / Healthy / Obese vs. cohort
        B. Paper-count distributions per condition (box + jitter)
        C. Literature vs cohort edge overlap (Healthy and Obese)
        D. OO coefficient distributions per condition (violin)

  visualizations/cosine_network_statistics.png
      2-panel overview for the 11-organ OO-threshold cosine networks:
        A. Edge counts per case (total / edges defined / cohort)
        B. OO bootstrap-mean violin plots with elbow thresholds

  visualizations/overview_table.png
      Formatted summary table with two sections:
        Top    — 14-organ Reference Metabolic Network
        Bottom — 11-organ Cosine Networks

Run from the project root (Metabolic_Reference_Network/ folder):
    uv run python generate_overview_statistics.py
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

HERE    = Path(__file__).resolve().parent
VIZ_DIR = HERE / "visualizations"
VIZ_DIR.mkdir(exist_ok=True)

sys.path.insert(0, str(HERE))
from threshold_utils import kneedle_elbow

# ── colour palette ─────────────────────────────────────────────────────────────
C_NOCOND  = "#a78bfa"   # violet  — No Condition
C_HEALTHY = "#22d3ee"   # cyan    — Healthy
C_OBESE   = "#fbbf24"   # amber   — Obese
C_COHORT  = "#64748b"   # slate   — Cohort / reference
C_TOTAL   = "#94a3b8"   # light slate — total pairs searched

CASE_COLORS = {
    "met_healthy": "#38bdf8",
    "met_obese":   "#fb923c",
}

# ── data loaders ───────────────────────────────────────────────────────────────

def load_bootstrap(path: Path) -> tuple[dict, float, list[float]]:
    data   = json.loads(path.read_text(encoding="utf-8"))
    means  = sorted(v["mean"] for v in data.values())
    thresh = kneedle_elbow(means)
    return data, thresh, means


def load_lit(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def count_cohort_pairs(csv_path: Path) -> int:
    if not csv_path.exists():
        return 0
    import csv as _csv
    count = 0
    with open(csv_path, encoding="utf-8") as f:
        reader = _csv.reader(f)
        header = next(reader)
        cols   = header[1:]
        for row in reader:
            if not row:
                continue
            for j, val in enumerate(row[1:]):
                if val.strip() == "1" and j < len(cols):
                    ro, co = row[0].strip(), cols[j].strip()
                    if ro and co and ro != co and ro < co:
                        count += 1
    return count


# ── load all data ──────────────────────────────────────────────────────────────
print("[i] Loading cosine bootstrap data …")
MH_data, MH_thresh, MH_means = load_bootstrap(
    HERE / "Edge_cosine_met_reference_network/healthy/bootstrap_results_healthy.json")
MO_data, MO_thresh, MO_means = load_bootstrap(
    HERE / "Edge_cosine_met_reference_network/obese/bootstrap_results_obese.json")

print("[i] Loading reference-metabolic literature results …")
REF_ALL = load_lit(HERE / "reference_network_only_metabolic/metabolic_literature_results.json")

print("[i] Computing condition-filtered edge sets …")
from run_metabolic_lit_search import (
    CONDITION_KEYWORDS, load_edge_filter, EDGE_FILTER_CSV,
    HEALTHY_FILTER_CSV, OBESE_FILTER_CSV,
    _compile_patterns, _filter_by_condition,
)
cond_pats   = {c: _compile_patterns(CONDITION_KEYWORDS[c]) for c in ("healthy", "obese")}
all_pairs   = set(load_edge_filter(EDGE_FILTER_CSV))

lit_h_pairs: set = set()
lit_o_pairs: set = set()
for pair in all_pairs:
    key   = f"{pair[0]}|{pair[1]}"
    rev   = f"{pair[1]}|{pair[0]}"
    entry = REF_ALL.get(key) or REF_ALL.get(rev) or {}
    papers = entry.get("papers", [])
    if _filter_by_condition(papers, cond_pats["healthy"]):
        lit_h_pairs.add(pair)
    if _filter_by_condition(papers, cond_pats["obese"]):
        lit_o_pairs.add(pair)

COHORT_H = set(load_edge_filter(HEALTHY_FILTER_CSV)) if HEALTHY_FILTER_CSV.exists() else set()
COHORT_O = set(load_edge_filter(OBESE_FILTER_CSV))   if OBESE_FILTER_CSV.exists() else set()
cos_cohort_h = count_cohort_pairs(
    HERE / "Edge_cosine_met_reference_network/healthy/healthy_cohort_connections.csv")
cos_cohort_o = count_cohort_pairs(
    HERE / "Edge_cosine_met_reference_network/obese/obese_cohort_connections.csv")


print("[i] Computing paper-count distributions …")
ppc_nocond  = [v.get("n_papers_found", 0) for v in REF_ALL.values()]

def _cond_papers(condition: str) -> list[int]:
    pats = cond_pats[condition]
    out  = []
    for pair in all_pairs:
        key  = f"{pair[0]}|{pair[1]}"; rev = f"{pair[1]}|{pair[0]}"
        entry = REF_ALL.get(key) or REF_ALL.get(rev) or {}
        out.append(len(_filter_by_condition(entry.get("papers", []), pats)))
    return out

ppc_healthy = _cond_papers("healthy")
ppc_obese   = _cond_papers("obese")

print("[i] Computing OO coefficients for 14-organ network …")
def _oo_for_condition(pair_papers: dict) -> list[float]:
    """pair_papers: {pair -> [paper, ...]}  →  list of OO values."""
    organ_pmids: dict = {}
    for (o1, o2), papers in pair_papers.items():
        pmids = {p.get("pmid") for p in papers if p.get("pmid")}
        organ_pmids.setdefault(o1, set()).update(pmids)
        organ_pmids.setdefault(o2, set()).update(pmids)
    organ_n = {o: len(s) for o, s in organ_pmids.items()}
    oo_vals = []
    for (o1, o2), papers in pair_papers.items():
        n_co  = len([p for p in papers if p.get("pmid")])
        n1, n2 = organ_n.get(o1, 0), organ_n.get(o2, 0)
        denom = (n1 * n2) ** 0.5 if n1 > 0 and n2 > 0 else 0
        oo_vals.append((n_co / denom) if denom > 0 else 0.0)
    return sorted(oo_vals)

def _pair_papers_for_cond(condition: str | None) -> dict:
    pats = cond_pats[condition] if condition else None
    out  = {}
    for pair in all_pairs:
        key  = f"{pair[0]}|{pair[1]}"; rev = f"{pair[1]}|{pair[0]}"
        entry = REF_ALL.get(key) or REF_ALL.get(rev) or {}
        papers = entry.get("papers", [])
        if pats is not None:
            papers = _filter_by_condition(papers, pats)
        out[pair] = papers
    return out

def _oo_threshold_pairs(condition: str) -> tuple[set, float]:
    """OO-elbow-selected pairs for the 14-organ network under a given condition."""
    pair_pp = _pair_papers_for_cond(condition)
    organ_pmids: dict = {}
    for (o1, o2), papers in pair_pp.items():
        pmids = {p.get("pmid") for p in papers if p.get("pmid")}
        organ_pmids.setdefault(o1, set()).update(pmids)
        organ_pmids.setdefault(o2, set()).update(pmids)
    organ_n = {o: len(s) for o, s in organ_pmids.items()}
    oo_map: dict = {}
    for (o1, o2), papers in pair_pp.items():
        n_co = len([p for p in papers if p.get("pmid")])
        n1, n2 = organ_n.get(o1, 0), organ_n.get(o2, 0)
        denom = (n1 * n2) ** 0.5 if n1 > 0 and n2 > 0 else 0
        oo_map[(o1, o2)] = (n_co / denom) if denom > 0 else 0.0
    nonzero = sorted(v for v in oo_map.values() if v > 0)
    thresh  = kneedle_elbow(nonzero) if len(nonzero) >= 3 else 0.0
    selected = {pair for pair, oo in oo_map.items() if oo >= thresh and oo > 0}
    return selected, thresh

oo_nocond  = _oo_for_condition(_pair_papers_for_cond(None))
oo_healthy = _oo_for_condition(_pair_papers_for_cond("healthy"))
oo_obese   = _oo_for_condition(_pair_papers_for_cond("obese"))

oo_h_pairs, oo_h_thresh = _oo_threshold_pairs("healthy")
oo_o_pairs, oo_o_thresh = _oo_threshold_pairs("obese")

print("[i] Computing literature vs cohort overlap …")
cohort_organs = {o for p in COHORT_H | COHORT_O for o in p}
# Use OO-threshold-selected pairs (14-organ) restricted to cohort organ set for comparison
oo_h_11   = {p for p in oo_h_pairs if p[0] in cohort_organs and p[1] in cohort_organs}
oo_o_11   = {p for p in oo_o_pairs if p[0] in cohort_organs and p[1] in cohort_organs}
shared_h  = oo_h_11 & COHORT_H;  only_lit_h = oo_h_11 - COHORT_H;  only_coh_h = COHORT_H - oo_h_11
shared_o  = oo_o_11 & COHORT_O;  only_lit_o = oo_o_11 - COHORT_O;  only_coh_o = COHORT_O - oo_o_11

ref_nocond_n = sum(1 for v in REF_ALL.values() if v.get("n_papers_found", 0) > 0)
MH_edges = sum(1 for v in MH_data.values() if v["mean"] >= MH_thresh and v["mean"] > 0)
MO_edges = sum(1 for v in MO_data.values() if v["mean"] >= MO_thresh and v["mean"] > 0)

print("[i] Data loaded.\n")


# ══════════════════════════════════════════════════════════════════════════════
#  Figure 1 — 14-organ Reference Metabolic Network overview
# ══════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(20, 16), facecolor="white")
gs  = fig.add_gridspec(2, 2, hspace=0.44, wspace=0.35,
                        left=0.07, right=0.97, top=0.93, bottom=0.07)
ax_edge = fig.add_subplot(gs[0, 0])
ax_ppc  = fig.add_subplot(gs[0, 1])
ax_cmp  = fig.add_subplot(gs[1, 0])
ax_oo   = fig.add_subplot(gs[1, 1])

# ── Panel A: Edge counts (14-organ, 3 conditions) ─────────────────────────────
cond_labels    = ["No Condition", "Healthy\n(OO ≥ elbow)", "Obese\n(OO ≥ elbow)"]
total_searched = [91, 91, 91]
edges_lit      = [ref_nocond_n, len(oo_h_pairs), len(oo_o_pairs)]
cohort_edges   = [91, len(COHORT_H), len(COHORT_O)]
cond_cols      = [C_NOCOND, C_HEALTHY, C_OBESE]

x = np.arange(3)
w = 0.26
bar_tot = ax_edge.bar(x - w, total_searched, width=w, label="Total pairs searched",
                      color=C_TOTAL, alpha=0.75, edgecolor="white", linewidth=0.5)
bar_lit = ax_edge.bar(x,     edges_lit,      width=w,
                      label="Edges defined (≥1 paper / OO ≥ elbow)",
                      color=cond_cols, alpha=0.85, edgecolor="white", linewidth=0.5)
bar_coh = ax_edge.bar(x + w, cohort_edges,   width=w, label="Cohort-defined edges",
                      color=C_COHORT, alpha=0.75, edgecolor="white", linewidth=0.5)

for bar_grp in (bar_tot, bar_lit, bar_coh):
    for bar in bar_grp:
        h = bar.get_height()
        ax_edge.text(bar.get_x() + bar.get_width() / 2, h + 0.5,
                     str(int(h)), ha="center", va="bottom", fontsize=8, color="#334155")

ax_edge.set_xticks(x)
ax_edge.set_xticklabels(cond_labels, fontsize=10)
ax_edge.set_ylabel("Number of edges / pairs", fontsize=10)
ax_edge.set_title("A  Edge counts — Reference Metabolic Network (14 organs, 91 pairs)",
                  fontsize=11, fontweight="bold", loc="left", color="#1e293b", pad=8)
ax_edge.legend(fontsize=8.5, framealpha=0.9)
ax_edge.set_ylim(0, 110)
ax_edge.spines[["top", "right"]].set_visible(False)

# ── Panel B: Paper-count distributions per condition ──────────────────────────
bp_data   = [ppc_nocond, ppc_healthy, ppc_obese]
bp_labels = ["No Condition\n(91 pairs)", "Healthy\n(91 pairs)", "Obese\n(91 pairs)"]
bp_cols   = [C_NOCOND, C_HEALTHY, C_OBESE]

bp = ax_ppc.boxplot(bp_data, patch_artist=True, widths=0.45,
                     medianprops=dict(color="#1e293b", linewidth=2),
                     whiskerprops=dict(color="#475569"),
                     capprops=dict(color="#475569"),
                     flierprops=dict(marker="o", markersize=3,
                                     markerfacecolor="#94a3b8", markeredgewidth=0))
for patch, col in zip(bp["boxes"], bp_cols):
    patch.set_facecolor(col); patch.set_alpha(0.6)

rng = np.random.default_rng(42)
for i, (dat, col) in enumerate(zip(bp_data, bp_cols), start=1):
    jitter = rng.uniform(-0.18, 0.18, len(dat))
    ax_ppc.scatter(np.full(len(dat), i) + jitter, dat,
                   color=col, alpha=0.35, s=12, edgecolors="none", zorder=2)

for i, dat in enumerate(bp_data, start=1):
    med = int(np.median(dat))
    pct = 100 * sum(1 for v in dat if v > 0) / len(dat)
    ax_ppc.text(i, max(dat) * 1.04, f"med={med}\n{pct:.0f}% >0",
                ha="center", va="bottom", fontsize=7.5, color="#334155")

ax_ppc.set_xticklabels(bp_labels, fontsize=9)
ax_ppc.set_ylabel("Papers passing filter (per pair)", fontsize=10)
ax_ppc.set_title("B  Paper-count distribution per condition",
                  fontsize=11, fontweight="bold", loc="left", color="#1e293b", pad=8)
ax_ppc.spines[["top", "right"]].set_visible(False)

# ── Panel C: Literature vs cohort edge overlap ─────────────────────────────────
x2 = np.arange(2)
w2 = 0.45
bsh = ax_cmp.bar(x2, [len(shared_h),   len(shared_o)],  width=w2,
                 label="Shared (OO-threshold & cohort)", color="#22c55e", alpha=0.85)
bol = ax_cmp.bar(x2, [len(only_lit_h), len(only_lit_o)], width=w2,
                 bottom=[len(shared_h), len(shared_o)],
                 label="OO-threshold-only", color=C_HEALTHY, alpha=0.85)
boc = ax_cmp.bar(x2, [len(only_coh_h), len(only_coh_o)], width=w2,
                 bottom=[len(shared_h)+len(only_lit_h), len(shared_o)+len(only_lit_o)],
                 label="Cohort-only", color=C_OBESE, alpha=0.85)

for bars, vals, bottoms in (
        (bsh, [len(shared_h),   len(shared_o)],   [0, 0]),
        (bol, [len(only_lit_h), len(only_lit_o)],  [len(shared_h), len(shared_o)]),
        (boc, [len(only_coh_h), len(only_coh_o)],
              [len(shared_h)+len(only_lit_h), len(shared_o)+len(only_lit_o)])):
    for bar, val, bot in zip(bars, vals, bottoms):
        if val > 0:
            ax_cmp.text(bar.get_x() + bar.get_width()/2, bot + val/2,
                        str(val), ha="center", va="center",
                        fontsize=10, fontweight="bold", color="white")

for i, (ln, cn) in enumerate([(len(oo_h_11), len(COHORT_H)),
                               (len(oo_o_11), len(COHORT_O))]):
    tot = len(shared_h if i==0 else shared_o) + \
          len(only_lit_h if i==0 else only_lit_o) + \
          len(only_coh_h if i==0 else only_coh_o)
    ax_cmp.text(i, tot + 0.8, f"OO: {ln}  Cohort: {cn}",
                ha="center", fontsize=8, color="#334155")

ax_cmp.set_xticks(x2)
ax_cmp.set_xticklabels(["Healthy\n(within 11-organ cohort set)",
                         "Obese\n(within 11-organ cohort set)"], fontsize=9.5)
ax_cmp.set_ylabel("Number of edges", fontsize=10)
ax_cmp.set_title("C  OO-threshold vs cohort edge overlap",
                  fontsize=11, fontweight="bold", loc="left", color="#1e293b", pad=8)
ax_cmp.legend(fontsize=8.5, framealpha=0.9)
ax_cmp.spines[["top", "right"]].set_visible(False)

# ── Panel D: OO coefficient distribution per condition (violin) ────────────────
violin_data = [oo_nocond, oo_healthy, oo_obese]
viol_cols   = [C_NOCOND, C_HEALTHY, C_OBESE]
viol_labels = ["No Condition\n(91 pairs)", "Healthy\n(91 pairs)", "Obese\n(91 pairs)"]
viol_thresh = [None, oo_h_thresh, oo_o_thresh]  # elbow for healthy/obese

parts = ax_oo.violinplot(violin_data, positions=[1, 2, 3], widths=0.55,
                          showmedians=True, showextrema=False)
for pc, col in zip(parts["bodies"], viol_cols):
    pc.set_facecolor(col); pc.set_alpha(0.5)
parts["cmedians"].set_color("#1e293b"); parts["cmedians"].set_linewidth(1.5)

for i, (dat, col, thresh) in enumerate(zip(violin_data, viol_cols, viol_thresh), start=1):
    median = float(np.median(dat))
    ax_oo.text(i + 0.33, median, f"{median:.4f}", fontsize=7.5, va="center", color=col)
    n_pos = sum(1 for v in dat if v > 0)
    ax_oo.text(i, max(dat) * 1.02 if dat else 0.01,
               f"{n_pos}/{len(dat)} pairs >0",
               ha="center", va="bottom", fontsize=7.5, color="#334155")
    if thresh is not None:
        ax_oo.hlines(thresh, i - 0.35, i + 0.35,
                     colors=col, linewidths=2, linestyles="--", zorder=3)
        ax_oo.text(i - 0.38, thresh, f"elbow\n{thresh:.4f}",
                   fontsize=7, va="center", ha="right", color=col)

ax_oo.set_xticks([1, 2, 3])
ax_oo.set_xticklabels(viol_labels, fontsize=9)
ax_oo.set_ylabel("OO coefficient", fontsize=10)
ax_oo.set_title("D  OO coefficient distribution per condition (dashed = elbow threshold)",
                fontsize=11, fontweight="bold", loc="left", color="#1e293b", pad=8)
ax_oo.spines[["top", "right"]].set_visible(False)

fig.suptitle("Reference Metabolic Network (14 organs) — Overview Statistics",
             fontsize=14, fontweight="bold", color="#1e293b", y=0.97)

out1 = VIZ_DIR / "overview_statistics.png"
fig.savefig(out1, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"[ok] {out1.relative_to(HERE)}")


# ══════════════════════════════════════════════════════════════════════════════
#  Figure 2 — 11-organ cosine network statistics (Met-Healthy / Met-Obese)
# ══════════════════════════════════════════════════════════════════════════════
fig2, (ax2_edge, ax2_oo) = plt.subplots(1, 2, figsize=(14, 8), facecolor="white")
fig2.subplots_adjust(left=0.08, right=0.97, top=0.88, bottom=0.12, wspace=0.38)

# ── Panel A: Edge counts ───────────────────────────────────────────────────────
case_labels = ["Met-Healthy", "Met-Obese"]
total_cos   = [55, 55]
edges_cos   = [MH_edges, MO_edges]
cohort_cos  = [cos_cohort_h, cos_cohort_o]
cos_cols    = [CASE_COLORS["met_healthy"], CASE_COLORS["met_obese"]]

xc = np.arange(2)
wc = 0.26
bc_tot = ax2_edge.bar(xc - wc, total_cos,  width=wc, label="Total pairs searched",
                      color=C_TOTAL, alpha=0.75, edgecolor="white", linewidth=0.5)
bc_lit = ax2_edge.bar(xc,       edges_cos, width=wc, label="Edges defined (OO ≥ elbow)",
                      color=cos_cols, alpha=0.85, edgecolor="white", linewidth=0.5)
bc_coh = ax2_edge.bar(xc + wc,  cohort_cos, width=wc, label="Cohort-defined edges",
                      color=C_COHORT, alpha=0.75, edgecolor="white", linewidth=0.5)

for bg in (bc_tot, bc_lit, bc_coh):
    for bar in bg:
        h = bar.get_height()
        ax2_edge.text(bar.get_x() + bar.get_width()/2, h + 0.5, str(int(h)),
                      ha="center", va="bottom", fontsize=9, color="#334155")

ax2_edge.set_xticks(xc)
ax2_edge.set_xticklabels(case_labels, fontsize=11)
ax2_edge.set_ylabel("Number of edges / pairs", fontsize=10)
ax2_edge.set_title("A  Edge counts — Metabolic cosine network (11 organs, 55 pairs)",
                   fontsize=11, fontweight="bold", loc="left", color="#1e293b", pad=8)
ax2_edge.legend(fontsize=9, framealpha=0.9)
ax2_edge.set_ylim(0, 68)
ax2_edge.spines[["top", "right"]].set_visible(False)

# ── Panel B: OO violin distributions ──────────────────────────────────────────
vd_cos     = [MH_means, MO_means]
vl_cos     = ["Met-Healthy\n(n=55)", "Met-Obese\n(n=55)"]
thresh_cos = [MH_thresh, MO_thresh]

vp = ax2_oo.violinplot(vd_cos, positions=[1, 2], widths=0.55,
                        showmedians=True, showextrema=False)
for pc, col in zip(vp["bodies"], cos_cols):
    pc.set_facecolor(col); pc.set_alpha(0.5)
vp["cmedians"].set_color("#1e293b"); vp["cmedians"].set_linewidth(1.5)

for i, (thresh, col, n_edges, dat) in enumerate(
        zip(thresh_cos, cos_cols, [MH_edges, MO_edges], vd_cos), start=1):
    ax2_oo.hlines(thresh, i - 0.35, i + 0.35,
                  colors=col, linewidths=2, linestyles="--")
    ax2_oo.text(i + 0.38, thresh, f"elbow {thresh:.4f}", fontsize=8, va="center", color=col)
    ax2_oo.text(i, max(dat) * 1.02,
                f"{n_edges}/55 edges",
                ha="center", va="bottom", fontsize=8, color="#334155")

ax2_oo.set_xticks([1, 2])
ax2_oo.set_xticklabels(vl_cos, fontsize=10)
ax2_oo.set_ylabel("Bootstrap mean OO coefficient", fontsize=10)
ax2_oo.set_title("B  OO coefficient distributions (dashed = elbow threshold)",
                 fontsize=11, fontweight="bold", loc="left", color="#1e293b", pad=8)
ax2_oo.spines[["top", "right"]].set_visible(False)

fig2.suptitle("Metabolic Cosine Network (11 organs) — Overview Statistics",
              fontsize=14, fontweight="bold", color="#1e293b", y=0.97)

out2 = VIZ_DIR / "cosine_network_statistics.png"
fig2.savefig(out2, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig2)
print(f"[ok] {out2.relative_to(HERE)}")


# ══════════════════════════════════════════════════════════════════════════════
#  Figure 3 — summary table (two sections)
# ══════════════════════════════════════════════════════════════════════════════
table_rows_ref = [
    ["No Condition", "14 organs", "91", str(ref_nocond_n),      "91 (all pairs)", "≥1 paper (same-sentence filter)"],
    ["Healthy",      "14 organs", "91", str(len(oo_h_pairs)),   str(len(COHORT_H)), f"OO ≥ {oo_h_thresh:.4f} (elbow)"],
    ["Obese",        "14 organs", "91", str(len(oo_o_pairs)),   str(len(COHORT_O)), f"OO ≥ {oo_o_thresh:.4f} (elbow)"],
]
table_rows_cos = [
    ["Met-Healthy", "11 organs", "55", str(MH_edges), str(cos_cohort_h), f"OO ≥ {MH_thresh:.4f}"],
    ["Met-Obese",   "11 organs", "55", str(MO_edges), str(cos_cohort_o), f"OO ≥ {MO_thresh:.4f}"],
]

col_headers = ["Case / Condition", "Organ set", "Pairs\nsearched",
               "Edges\n(literature)", "Edges\n(cohort)", "Edge criterion"]
col_widths  = [0.17, 0.10, 0.09, 0.12, 0.10, 0.34]

all_rows = table_rows_ref + [None] + table_rows_cos  # None = separator

fig3, ax = plt.subplots(figsize=(16, 6.5), facecolor="white")
ax.axis("off")
xs = np.cumsum([0] + col_widths[:-1])

hdr_color_ref = "#4338ca"  # indigo for reference network header
hdr_color_cos = "#0f766e"  # teal for cosine header

def _draw_section_header(ax, label, y_top, color):
    ax.add_patch(plt.Rectangle((0, y_top), 1.0, 0.06,
                                transform=ax.transAxes, color=color, clip_on=False))
    ax.text(0.5, y_top + 0.03, label, transform=ax.transAxes,
            ha="center", va="center", fontsize=10, fontweight="bold", color="white")

def _draw_table_header(ax, y_top, bg="#1e293b"):
    for x, w, h in zip(xs, col_widths, col_headers):
        ax.add_patch(plt.Rectangle((x, y_top), w, 0.055,
                                    transform=ax.transAxes, color=bg, clip_on=False))
        ax.text(x + w/2, y_top + 0.027, h, transform=ax.transAxes,
                ha="center", va="center", fontsize=8.5, fontweight="bold", color="white")

def _draw_rows(ax, rows, y_start, row_h):
    row_colors = ["#f8fafc", "#f1f5f9"]
    ri = 0
    for row in rows:
        y_top = y_start - ri * row_h
        bg = row_colors[ri % 2]
        for x, w, cell in zip(xs, col_widths, row):
            ax.add_patch(plt.Rectangle((x, y_top - row_h), w, row_h,
                                        transform=ax.transAxes, facecolor=bg, clip_on=False,
                                        linewidth=0.5, edgecolor="#e2e8f0"))
            ax.text(x + w/2, y_top - row_h/2, cell, transform=ax.transAxes,
                    ha="center", va="center", fontsize=8.5, color="#1e293b")
        ri += 1

# Layout: reference section at top, cosine section below with a gap
row_h = 0.09
sec1_hdr_y = 0.97
sec1_tbl_y = 0.91   # top of column headers
sec1_rows_y = 0.855 # top of first data row
# reference: 3 rows
sec2_hdr_y  = sec1_rows_y - 3 * row_h - 0.05
sec2_tbl_y  = sec2_hdr_y - 0.06
sec2_rows_y = sec2_tbl_y - 0.055
# cosine: 4 rows

_draw_section_header(ax, "Reference Metabolic Network (14 organs, 91 pairs)", sec1_hdr_y, hdr_color_ref)
_draw_table_header(ax, sec1_tbl_y)
_draw_rows(ax, table_rows_ref, sec1_rows_y, row_h)

_draw_section_header(ax, "OO-Threshold Cosine Networks (11 organs, 55 pairs)", sec2_hdr_y, hdr_color_cos)
_draw_table_header(ax, sec2_tbl_y, bg="#134e4a")
_draw_rows(ax, table_rows_cos, sec2_rows_y, row_h)

# footnote
note_y = sec2_rows_y - 4 * row_h - 0.04
notes = [
    f"Healthy OO-vs-cohort (within 11-organ set): Shared {len(shared_h)}  |  OO-only {len(only_lit_h)}  |  Cohort-only {len(only_coh_h)}  (OO elbow: {oo_h_thresh:.4f})",
    f"Obese OO-vs-cohort   (within 11-organ set): Shared {len(shared_o)}  |  OO-only {len(only_lit_o)}  |  Cohort-only {len(only_coh_o)}  (OO elbow: {oo_o_thresh:.4f})",
    f"Paper counts (No Condition, 91 pairs): median {int(np.median(ppc_nocond))}  |  mean {np.mean(ppc_nocond):.1f}  |  max {max(ppc_nocond)}  |  pairs with ≥1 paper {ref_nocond_n}/91",
]
for li, line in enumerate(notes):
    ax.text(0.5, note_y - li * 0.045, line, transform=ax.transAxes,
            ha="center", va="top", fontsize=8, color="#334155", style="italic")

ax.set_title("Metabolic Reference Network — Summary Statistics Table",
             fontsize=13, fontweight="bold", color="#1e293b", pad=14)

out3 = VIZ_DIR / "overview_table.png"
fig3.savefig(out3, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig3)
print(f"[ok] {out3.relative_to(HERE)}")

print("\n[ok] All statistics figures saved.")
