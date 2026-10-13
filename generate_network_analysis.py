"""
generate_network_analysis.py
============================
Compares four networks across two conditions:

  Healthy:
    - Pre-defined cohort network  (healthy_cohort_connections.csv)
    - OO-threshold reference       (reference metabolic, healthy papers, elbow cutoff)

  Obese:
    - Pre-defined cohort network  (obese_cohort_connections.csv)
    - OO-threshold reference       (reference metabolic, obese papers, elbow cutoff)

Metrics computed per network (Barabási Network Science definitions):
  N        : number of nodes (organs with at least one edge)
  E        : number of edges
  Density  : E / (N*(N-1)/2)
  <k>      : mean degree  = 2E/N
  <k²>     : mean squared degree
  κ        : degree heterogeneity = <k²>/<k>   (Barabási Ch.4)

Output:
  visualizations/network_analysis.png   – comparison figure
"""

import csv
import json
import sys
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import numpy as np
from scipy import stats as scipy_stats
from scipy.optimize import curve_fit
import pandas as pd

HERE    = Path(__file__).resolve().parent
OUT_DIR = HERE / "visualizations_all_papers"
DPI     = 300
sys.path.insert(0, str(HERE))


# ── Graph helpers ─────────────────────────────────────────────────────────────

def weighted_edges_from_pc_csv(
    pc_path: Path, sig_path: Path
) -> tuple[list[tuple[str, str, float]], list[str]]:
    """Load partial-correlation edges with |r| weights for significant pairs (sig==1).
    Returns (edges_weighted, organs).
    """
    pc  = pd.read_csv(pc_path,  index_col=0)
    sig = pd.read_csv(sig_path, index_col=0)
    organs = list(pc.columns)
    edges  = []
    for i, o1 in enumerate(organs):
        for j, o2 in enumerate(organs):
            if j > i and sig.loc[o1, o2] == 1:
                edges.append((o1, o2, abs(float(pc.loc[o1, o2]))))
    return sorted(edges), organs


def edges_from_csv(csv_path: Path) -> tuple[list[tuple[str, str]], list[str]]:
    """Read upper-triangular adjacency CSV.
    Returns (edges, organs) where organs is the full node list from the header."""
    with open(csv_path, encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        col_organs = [h.strip() for h in header[1:]]
        edges = []
        for row in reader:
            if not row:
                continue
            row_organ = row[0].strip()
            for col_idx, val in enumerate(row[1:]):
                if val.strip() == "1" and col_idx < len(col_organs):
                    col_organ = col_organs[col_idx]
                    if row_organ and col_organ and row_organ != col_organ:
                        edges.append((min(row_organ, col_organ), max(row_organ, col_organ)))
    return sorted(set(edges)), col_organs


def network_metrics(edges: list[tuple[str, str]],
                    node_universe: list[str] | None = None) -> dict:
    """Compute all network metrics for an undirected, unweighted graph.
    node_universe fixes the node set (including isolated nodes); if None,
    only nodes appearing in edges are counted."""
    nodes = sorted(node_universe) if node_universe else sorted({o for pair in edges for o in pair})
    N = len(nodes)
    E = len(edges)

    degree: dict[str, int] = defaultdict(int)
    for o1, o2 in edges:
        degree[o1] += 1
        degree[o2] += 1

    k_vals = [degree[n] for n in nodes]
    k_mean  = np.mean(k_vals) if k_vals else 0.0
    k2_mean = np.mean([k ** 2 for k in k_vals]) if k_vals else 0.0
    kappa   = k2_mean / k_mean if k_mean > 0 else 0.0
    max_E   = N * (N - 1) / 2
    density = E / max_E if max_E > 0 else 0.0

    return {
        "nodes":   nodes,
        "edges":   edges,
        "degree":  degree,
        "N":       N,
        "E":       E,
        "density": density,
        "k_mean":  k_mean,
        "k2_mean": k2_mean,
        "kappa":   kappa,
    }


def compute_disorder(metrics: dict) -> dict:
    """
    Fit an exponential decay to the ranked node-strength distribution.

    With unweighted edges (w = 1 for all edges), strength equals degree.
    Nodes with strength = 0 (isolated in the network) are excluded.

    Ranked sequence: s_(1) >= s_(2) >= ... >= s_(M)
    Fit:             s_(r) = A * exp(-beta * r)
    Linearised OLS:  ln(s_(r)) = ln(A) - beta * r

    Returns a dict with:
        beta      : disorder exponent (positive = decaying distribution)
        amp       : amplitude A
        r2        : coefficient of determination R² of the log-linear fit
        ranks     : rank array [1..M] used in the fit
        strengths : sorted strengths [s_(1)..s_(M)] used in the fit
    """
    degree = metrics["degree"]
    nodes  = metrics["nodes"]

    # Strength = degree (unweighted: all present edges contribute w = 1)
    strengths_pos = sorted(
        [degree.get(n, 0) for n in nodes if degree.get(n, 0) > 0],
        reverse=True,
    )

    if len(strengths_pos) < 2:
        return {"beta": 0.0, "amp": 0.0, "r2": 0.0, "ranks": [], "strengths": []}

    ranks  = np.arange(1, len(strengths_pos) + 1, dtype=float)
    log_s  = np.log(np.array(strengths_pos, dtype=float))

    # OLS on ln(s) = ln(A) - beta*r  →  coeffs[0] = -beta, coeffs[1] = ln(A)
    coeffs = np.polyfit(ranks, log_s, 1)
    beta   = float(-coeffs[0])
    amp    = float(np.exp(coeffs[1]))

    # R²
    y_pred  = np.polyval(coeffs, ranks)
    ss_res  = float(np.sum((log_s - y_pred) ** 2))
    ss_tot  = float(np.sum((log_s - np.mean(log_s)) ** 2))
    r2      = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return {
        "beta":      beta,
        "amp":       amp,
        "r2":        r2,
        "ranks":     ranks.tolist(),
        "strengths": strengths_pos,
    }


# ── Full condition network ────────────────────────────────────────────────────

def full_condition_network(condition: str, all_organs: list[str],
                           restrict_to: list[str] | None = None) -> tuple[list[tuple[str, str]], dict]:
    """
    Build a full condition network (no OO threshold) from metabolic literature results.

    restrict_to: if given, only pairs between these organs are used, and the
    per-organ paper totals behind the OO score are computed from those pairs
    only (so the OO values do not depend on organs outside the set).

    Returns (edges, oo_map) where:
      - edges: sorted list of (o1, o2) pairs with ≥1 relevant paper
      - oo_map: dict with (o1,o2) and (o2,o1) both mapped to the OO score
    """
    from run_metabolic_lit_search import (
        CONDITION_KEYWORDS, _compile_patterns, _filter_by_condition,
    )

    ref_dir      = HERE / "reference_network_only_metabolic"
    results_path = ref_dir / "metabolic_literature_results.json"
    data = json.loads(results_path.read_text(encoding="utf-8"))

    if condition != "all":
        pats = _compile_patterns(CONDITION_KEYWORDS[condition])

    pair_papers: dict[tuple, list] = {}
    for key, entry in data.items():
        parts = key.split("|", 1)
        if len(parts) != 2:
            continue
        o1, o2 = min(parts), max(parts)
        if restrict_to is not None and (o1 not in restrict_to or o2 not in restrict_to):
            continue
        papers = entry.get("papers", [])
        if condition != "all":
            papers = _filter_by_condition(papers, pats)
        pair_papers[(o1, o2)] = papers

    organ_pmids: dict[str, set] = {}
    for (o1, o2), papers in pair_papers.items():
        pmids = {p.get("pmid") for p in papers if p.get("pmid")}
        organ_pmids.setdefault(o1, set()).update(pmids)
        organ_pmids.setdefault(o2, set()).update(pmids)
    organ_n = {o: len(s) for o, s in organ_pmids.items()}

    oo_map: dict[tuple, float] = {}
    edges: list[tuple[str, str]] = []
    for (o1, o2), papers in pair_papers.items():
        n_co = len([p for p in papers if p.get("pmid")])
        if n_co < 1:
            continue
        n1, n2 = organ_n.get(o1, 0), organ_n.get(o2, 0)
        denom   = (n1 * n2) ** 0.5 if n1 > 0 and n2 > 0 else 0
        score   = (n_co / denom) if denom > 0 else 0.0
        oo_map[(o1, o2)] = score
        oo_map[(o2, o1)] = score
        edges.append((o1, o2))

    return sorted(edges), oo_map


def compute_disorder_weighted(nodes: list[str], strength: dict) -> dict:
    """
    Fit an exponential decay to the ranked node-strength distribution.

    Like compute_disorder() but takes a strength dict {organ: float}
    instead of using metrics["degree"].

    Ranked sequence: s_(1) >= s_(2) >= ... >= s_(M)
    Fit:             s_(r) = A * exp(-beta * r)
    Linearised OLS:  ln(s_(r)) = ln(A) - beta * r

    Returns a dict with:
        beta      : disorder exponent (positive = decaying distribution)
        amp       : amplitude A
        r2        : coefficient of determination R² of the log-linear fit
        ranks     : rank array [1..M] used in the fit
        strengths : sorted strengths [s_(1)..s_(M)] used in the fit
    """
    strengths_pos = sorted(
        [strength.get(n, 0.0) for n in nodes if strength.get(n, 0.0) > 0],
        reverse=True,
    )

    if len(strengths_pos) < 2:
        return {"beta": 0.0, "amp": 0.0, "r2": 0.0, "ranks": [], "strengths": []}

    ranks  = np.arange(1, len(strengths_pos) + 1, dtype=float)
    log_s  = np.log(np.array(strengths_pos, dtype=float))

    coeffs = np.polyfit(ranks, log_s, 1)
    beta   = float(-coeffs[0])
    amp    = float(np.exp(coeffs[1]))

    y_pred  = np.polyval(coeffs, ranks)
    ss_res  = float(np.sum((log_s - y_pred) ** 2))
    ss_tot  = float(np.sum((log_s - np.mean(log_s)) ** 2))
    r2      = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return {
        "beta":      beta,
        "amp":       amp,
        "r2":        r2,
        "ranks":     ranks.tolist(),
        "strengths": strengths_pos,
    }


# ── New metric functions (provided) ──────────────────────────────────────────

def calculate_network_density(
    weight_matrix: np.ndarray,
    directed: bool = False,
    include_weights: bool = True
) -> float:
    weight_matrix = np.array(weight_matrix, dtype=float)

    if not include_weights:
        num_nodes = int(weight_matrix.shape[0])
        num_edges = np.count_nonzero(weight_matrix) - np.count_nonzero(np.diag(weight_matrix))
        if directed:
            max_edges = num_nodes * (num_nodes - 1)
        else:
            max_edges = num_nodes * (num_nodes - 1) / 2
        return num_edges / max_edges if max_edges > 0 else 0.0

    max_weight = np.max(weight_matrix)
    if max_weight == 0:
        return 0.0

    norm_matrix = weight_matrix / max_weight
    num_nodes = int(norm_matrix.shape[0])

    node_strength = np.sum(np.abs(norm_matrix))
    if not directed:
        node_strength /= 2

    max_edges = num_nodes * (num_nodes - 1) if directed else num_nodes * (num_nodes - 1) / 2
    density = node_strength / max_edges

    return density * 100


def calculate_disorder(network: pd.DataFrame):
    edges_count = network.sum(axis=1).to_dict()

    def exp_decay(x, a, b, c):
        return a * np.exp(-b * x) + c

    x = np.arange(len(network.columns))
    y = list(edges_count.values())

    try:
        p0 = (float(max(y)), 1.0, float(min(y)))
        popt, _ = curve_fit(exp_decay, x, y, p0=p0, maxfev=5000)
    except (RuntimeError, ValueError):
        popt = (1.0, 0.0, float(np.mean(y)))

    x_fit = np.linspace(min(x), max(x), 200)
    y_fit = exp_decay(x_fit, *popt)

    return edges_count, x, y, x_fit, y_fit


def _compute_disorder_edges(edges_weighted: list, nodes: list) -> dict:
    """
    Wrapper: build a weight-matrix from edges, sort by strength, call
    calculate_disorder(), and return the disorder dict used by plotting functions.

    edges_weighted: list of (o1, o2) for binary or (o1, o2, w) for weighted.
    """
    n   = len(nodes)
    idx = {o: i for i, o in enumerate(nodes)}
    W   = np.zeros((n, n), dtype=float)
    for e in edges_weighted:
        w = e[2] if len(e) == 3 else 1.0
        i, j = idx.get(e[0], -1), idx.get(e[1], -1)
        if i >= 0 and j >= 0:
            W[i, j] = w
            W[j, i] = w

    df    = pd.DataFrame(W, index=nodes, columns=nodes)
    order = df.sum(axis=1).sort_values(ascending=False).index
    df    = df.loc[order, order]

    def _exp_decay(x, a, b, c):
        return a * np.exp(-b * x) + c

    _ec, x, y, x_fit, y_fit = calculate_disorder(df)

    try:
        p0   = (float(max(y)), 1.0, float(min(y)))
        popt, _ = curve_fit(_exp_decay, x, y, p0=p0, maxfev=5000)
    except (RuntimeError, ValueError):
        popt = (1.0, 0.0, float(np.mean(y)))

    y_pred = _exp_decay(x, *popt)
    ss_res = float(np.sum((np.array(y, dtype=float) - y_pred) ** 2))
    ss_tot = float(np.sum((np.array(y, dtype=float) - np.mean(y)) ** 2))
    r2     = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return {
        "beta":      float(popt[1]),
        "amp":       float(popt[0]),
        "offset":    float(popt[2]),
        "r2":        float(r2),
        "ranks":     x.tolist(),
        "strengths": list(y),
        "x_fit":     x_fit.tolist(),
        "y_fit":     y_fit.tolist(),
    }


def _weighted_density(edges_weighted: list, nodes: list) -> float:
    """Build weight matrix and return weighted density via calculate_network_density()."""
    n   = len(nodes)
    idx = {o: i for i, o in enumerate(nodes)}
    W   = np.zeros((n, n), dtype=float)
    for e in edges_weighted:
        w = e[2] if len(e) == 3 else 1.0
        i, j = idx.get(e[0], -1), idx.get(e[1], -1)
        if i >= 0 and j >= 0:
            W[i, j] = w
            W[j, i] = w
    return calculate_network_density(W, directed=False, include_weights=True)


# ── OO-threshold network ──────────────────────────────────────────────────────

def oo_threshold_edges_all() -> tuple[list[tuple[str, str]], float]:
    """
    OO threshold on all 91 pairs with no condition keyword filtering —
    uses every paper found per pair regardless of health condition.
    """
    from threshold_utils import kneedle_elbow

    ref_dir      = HERE / "reference_network_only_metabolic"
    results_path = ref_dir / "metabolic_literature_results.json"
    data = json.loads(results_path.read_text(encoding="utf-8"))

    pair_papers: dict[tuple, list] = {}
    for key, entry in data.items():
        parts = key.split("|", 1)
        if len(parts) != 2:
            continue
        o1, o2 = min(parts), max(parts)
        pair_papers[(o1, o2)] = entry.get("papers", [])

    organ_pmids: dict[str, set] = {}
    for (o1, o2), papers in pair_papers.items():
        pmids = {p.get("pmid") for p in papers if p.get("pmid")}
        organ_pmids.setdefault(o1, set()).update(pmids)
        organ_pmids.setdefault(o2, set()).update(pmids)
    organ_n = {o: len(s) for o, s in organ_pmids.items()}

    oo_map: dict[tuple, float] = {}
    for (o1, o2), papers in pair_papers.items():
        n_co  = len([p for p in papers if p.get("pmid")])
        n1, n2 = organ_n.get(o1, 0), organ_n.get(o2, 0)
        denom  = (n1 * n2) ** 0.5 if n1 > 0 and n2 > 0 else 0
        oo_map[(o1, o2)] = (n_co / denom) if denom > 0 else 0.0

    nonzero = sorted(v for v in oo_map.values() if v > 0)
    thresh  = kneedle_elbow(nonzero) if len(nonzero) >= 3 else 0.0

    selected = [(o1, o2) for (o1, o2), oo in oo_map.items() if oo >= thresh and oo > 0]
    return sorted(selected), thresh


def oo_threshold_edges(condition: str) -> tuple[list[tuple[str, str]], float]:
    """
    Compute OO coefficient per pair for the 14-organ reference network filtered
    by condition keywords, apply the kneedle elbow threshold, and return
    (selected_edges, threshold).
    """
    from threshold_utils import kneedle_elbow
    from run_metabolic_lit_search import (
        CONDITION_KEYWORDS, _compile_patterns, _filter_by_condition,
    )

    ref_dir      = HERE / "reference_network_only_metabolic"
    results_path = ref_dir / "metabolic_literature_results.json"
    data = json.loads(results_path.read_text(encoding="utf-8"))
    pats = _compile_patterns(CONDITION_KEYWORDS[condition])

    pair_papers: dict[tuple, list] = {}
    for key, entry in data.items():
        parts = key.split("|", 1)
        if len(parts) != 2:
            continue
        o1, o2 = min(parts), max(parts)
        pair_papers[(o1, o2)] = _filter_by_condition(entry.get("papers", []), pats)

    organ_pmids: dict[str, set] = {}
    for (o1, o2), papers in pair_papers.items():
        pmids = {p.get("pmid") for p in papers if p.get("pmid")}
        organ_pmids.setdefault(o1, set()).update(pmids)
        organ_pmids.setdefault(o2, set()).update(pmids)
    organ_n = {o: len(s) for o, s in organ_pmids.items()}

    oo_map: dict[tuple, float] = {}
    for (o1, o2), papers in pair_papers.items():
        n_co  = len([p for p in papers if p.get("pmid")])
        n1, n2 = organ_n.get(o1, 0), organ_n.get(o2, 0)
        denom  = (n1 * n2) ** 0.5 if n1 > 0 and n2 > 0 else 0
        oo_map[(o1, o2)] = (n_co / denom) if denom > 0 else 0.0

    nonzero = sorted(v for v in oo_map.values() if v > 0)
    thresh  = kneedle_elbow(nonzero) if len(nonzero) >= 3 else 0.0

    selected = [(o1, o2) for (o1, o2), oo in oo_map.items() if oo >= thresh and oo > 0]
    return sorted(selected), thresh


# ── Plotting ──────────────────────────────────────────────────────────────────

CONDITION_COLORS = {
    "all_full":       "#64748b",
    "healthy_predef": "#38bdf8",
    "healthy_full":   "#0284c7",
    "obese_predef":   "#fb923c",
    "obese_full":     "#c2410c",
}

NETWORK_LABELS = {
    "all_full":       "No Condition\n(Full network)",
    "healthy_predef": "Healthy\n(pre-defined)",
    "healthy_full":   "Healthy\n(Full network)",
    "obese_predef":   "Obese\n(pre-defined)",
    "obese_full":     "Obese\n(Full network)",
}


def _metric_bar(ax, values: list[float], labels: list[str], colors: list,
                title: str, ylabel: str, annotate_fmt: str = "{:.2f}"):
    bars = ax.bar(labels, values, color=colors, edgecolor="white", linewidth=0.8, width=0.55)
    ax.set_title(title, fontsize=11, pad=8, color="#1e293b", fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=9, color="#475569")
    ax.tick_params(labelsize=8, colors="#1e293b")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#cbd5e1")
    ax.grid(axis="y", color="#e2e8f0", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    max_v = max(values) if values else 1
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max_v * 0.02,
                annotate_fmt.format(v), ha="center", va="bottom", fontsize=8, color="#1e293b")


def _degree_dist(ax, metrics: dict, color: str, label: str):
    degree = metrics["degree"]
    k_vals = sorted(degree.values())
    k_unique, counts = np.unique(k_vals, return_counts=True)
    freq = counts / counts.sum()
    ax.bar(k_unique, freq, color=color, edgecolor="white", linewidth=0.5,
           alpha=0.85, label=label, width=0.4)
    ax.set_xlabel("Degree k", fontsize=9, color="#475569")
    ax.set_ylabel("P(k)", fontsize=9, color="#475569")
    ax.tick_params(labelsize=8, colors="#1e293b")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#cbd5e1")
    ax.grid(axis="y", color="#e2e8f0", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def _ranked_strength_ax(ax, metrics: dict, dis: dict, color: str):
    """Plot ranked strength distribution with exponential fit on a log-linear scale."""
    ranks     = dis["ranks"]
    strengths = dis["strengths"]

    if not ranks:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, color="#94a3b8", fontsize=10)
        ax.axis("off")
        return

    # Scatter — log scale: filter non-positive values
    pos_mask = [s > 0 for s in strengths]
    r_pos = [r for r, m in zip(ranks, pos_mask) if m]
    s_pos = [s for s, m in zip(strengths, pos_mask) if m]
    if r_pos:
        ax.semilogy(r_pos, s_pos, "o", color=color, markersize=5, alpha=0.9,
                    markeredgecolor="white", markeredgewidth=0.5, zorder=3)

    # Fit line using pre-computed x_fit / y_fit
    x_fit = dis.get("x_fit", [])
    y_fit = dis.get("y_fit", [])
    if x_fit and y_fit:
        xf = np.array(x_fit); yf = np.array(y_fit)
        pos = yf > 0
        if pos.any():
            ax.semilogy(xf[pos], yf[pos], "--", color=color, linewidth=1.8, alpha=0.75, zorder=2)

    ax.set_xlabel("Rank r", fontsize=9, color="#475569")
    ax.set_ylabel("Strength s(r)  [log]", fontsize=9, color="#475569")
    ax.tick_params(labelsize=8, colors="#1e293b")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#cbd5e1")
    ax.grid(axis="y", color="#e2e8f0", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def _strength_bar_ax(ax, metrics: dict, dis: dict, color: str):
    """Bar plot of organ strengths sorted high-to-low with exponential fit overlaid."""
    strength_dict = metrics.get("strength", metrics["degree"])
    nodes  = metrics["nodes"]

    # Sort organs by strength descending (include zero-strength nodes)
    sorted_nodes  = sorted(nodes, key=lambda n: strength_dict.get(n, 0), reverse=True)
    strengths_all = [strength_dict.get(n, 0) for n in sorted_nodes]
    x             = np.arange(len(sorted_nodes))

    # Shorten labels to fit: replace spaces with newline
    short_labels = [n.replace(" ", "\n") for n in sorted_nodes]

    bars = ax.bar(x, strengths_all, color=color, edgecolor="white",
                  linewidth=0.6, alpha=0.82, zorder=2, width=0.65)

    # Exponential fit line using pre-computed x_fit / y_fit from calculate_disorder
    x_fit = dis.get("x_fit", [])
    y_fit = dis.get("y_fit", [])
    if x_fit and y_fit:
        ax.plot(x_fit, y_fit, "--", color="#1e293b", linewidth=1.8,
                alpha=0.75, zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels(short_labels, fontsize=6.5, rotation=45, ha="right",
                       color="#1e293b", multialignment="center")
    ylabel = "OO-score strength" if "strength" in metrics else "Strength (degree)"
    ax.set_ylabel(ylabel, fontsize=9, color="#475569")
    ax.tick_params(axis="y", labelsize=8, colors="#1e293b")
    ax.tick_params(axis="x", length=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#cbd5e1")
    ax.grid(axis="y", color="#e2e8f0", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.set_xlim(-0.6, len(sorted_nodes) - 0.4)


_SHORT = {
    "all_full":       "All-full",
    "healthy_predef": "H-pre",
    "healthy_full":   "H-full",
    "obese_predef":   "O-pre",
    "obese_full":     "O-full",
}


def _sig_label(p: float) -> str:
    if np.isnan(p):   return "—"
    if p < 0.001:     return "***"
    if p < 0.01:      return "**"
    if p < 0.05:      return "*"
    return "ns"


def compute_pairwise_stats(all_metrics: dict, keys: list[str]) -> list[dict]:
    """
    For every pair of networks:
      - Welch's t-test on full degree sequences  (density comparison)
      - Mann-Whitney U on positive-strength sequences  (disorder comparison)
    """
    results = []
    for i, k1 in enumerate(keys):
        for k2 in keys[i + 1:]:
            m1, m2 = all_metrics[k1], all_metrics[k2]

            # Degree sequences (all nodes in their respective universes)
            deg1 = [m1["degree"].get(n, 0) for n in m1["nodes"]]
            deg2 = [m2["degree"].get(n, 0) for n in m2["nodes"]]
            t_stat, p_t = scipy_stats.ttest_ind(deg1, deg2, equal_var=False)

            # Strength sequences for MWU — exclude zero-strength nodes
            s1 = m1.get("strength_seq", [v for v in [m1["degree"].get(n, 0) for n in m1["nodes"]] if v > 0])
            s2 = m2.get("strength_seq", [v for v in [m2["degree"].get(n, 0) for n in m2["nodes"]] if v > 0])
            if s1 and s2:
                u_stat, p_u = scipy_stats.mannwhitneyu(s1, s2, alternative="two-sided")
            else:
                u_stat, p_u = float("nan"), float("nan")

            results.append({
                "label":         f"{_SHORT[k1]} vs {_SHORT[k2]}",
                "k1": k1, "k2": k2,
                "delta_density": m1["density"] - m2["density"],
                "t_stat":        float(t_stat),
                "p_t":           float(p_t),
                "delta_beta":    m1["disorder"]["beta"] - m2["disorder"]["beta"],
                "u_stat":        float(u_stat),
                "p_u":           float(p_u),
            })
    return results


def _stats_table(ax, results: list[dict]):
    """Draw the pairwise statistical comparison table."""
    ax.axis("off")

    col_headers = [
        "Comparison",
        "Δ Density", "t-stat", "p (t-test)", "sig",
        "Δ β", "U-stat", "p (MWU)", "sig",
    ]
    col_widths = [0.18, 0.09, 0.09, 0.10, 0.06, 0.09, 0.09, 0.10, 0.06]

    y_header = 0.96
    x0 = 0.02
    x = x0
    for header, w in zip(col_headers, col_widths):
        ax.text(x + w / 2, y_header, header, ha="center", va="center",
                fontsize=8.5, fontweight="bold", color="#1e293b",
                transform=ax.transAxes)
        x += w

    ax.plot([0.01, 0.99], [y_header - 0.04, y_header - 0.04],
            color="#94a3b8", linewidth=0.8, transform=ax.transAxes)

    # Colour-code significance
    sig_colors = {"***": "#dc2626", "**": "#ea580c", "*": "#ca8a04", "ns": "#64748b", "—": "#94a3b8"}

    n_rows = len(results)
    row_h  = (y_header - 0.09) / n_rows

    for row_i, r in enumerate(results):
        y = y_header - 0.07 - row_i * row_h - row_h / 2
        bg = "#f8fafc" if row_i % 2 == 0 else "#f1f5f9"
        rect = mpatches.FancyBboxPatch(
            (0.01, y - row_h / 2 + 0.004), 0.98, row_h - 0.008,
            boxstyle="square,pad=0", linewidth=0,
            facecolor=bg, alpha=0.6, transform=ax.transAxes,
        )
        ax.add_patch(rect)

        sig_t = _sig_label(r["p_t"])
        sig_u = _sig_label(r["p_u"])

        cells = [
            r["label"],
            f"{r['delta_density']*100:+.1f}%",
            f"{r['t_stat']:.2f}",
            f"{r['p_t']:.3f}" if r["p_t"] >= 0.001 else "<0.001",
            sig_t,
            f"{r['delta_beta']*100:+.2f}%",
            f"{int(r['u_stat'])}" if not np.isnan(r["u_stat"]) else "—",
            f"{r['p_u']:.3f}" if not np.isnan(r["p_u"]) and r["p_u"] >= 0.001 else ("<0.001" if not np.isnan(r["p_u"]) else "—"),
            sig_u,
        ]

        x = x0
        for ci, (cell, w) in enumerate(zip(cells, col_widths)):
            color = "#1e293b"
            fw    = "normal"
            if ci == 4:   # sig_t
                color = sig_colors.get(cell, "#1e293b")
                fw    = "bold" if cell not in ("ns", "—") else "normal"
            elif ci == 8:  # sig_u
                color = sig_colors.get(cell, "#1e293b")
                fw    = "bold" if cell not in ("ns", "—") else "normal"
            ax.text(x + w / 2, y, cell, ha="center", va="center",
                    fontsize=8, color=color, fontweight=fw,
                    transform=ax.transAxes)
            x += w

    # Vertical separator between density block and disorder block
    sep_x = x0 + sum(col_widths[:5])
    ax.plot([sep_x, sep_x], [0.02, y_header - 0.03],
            color="#cbd5e1", linewidth=0.8, transform=ax.transAxes)

    # Group headers
    density_mid = x0 + sum(col_widths[:1]) + sum(col_widths[1:5]) / 2
    disorder_mid = x0 + sum(col_widths[:5]) + sum(col_widths[5:]) / 2
    ax.text(density_mid, y_header + 0.03, "Density  (Student t-test on degree sequences)",
            ha="center", va="bottom", fontsize=8.5, color="#475569", style="italic",
            transform=ax.transAxes)
    ax.text(disorder_mid, y_header + 0.03, "Disorder  (Mann-Whitney U on strength sequences)",
            ha="center", va="bottom", fontsize=8.5, color="#475569", style="italic",
            transform=ax.transAxes)

    ax.text(0.99, 0.01, "* p<0.05   ** p<0.01   *** p<0.001   ns not significant",
            ha="right", va="bottom", fontsize=7.5, color="#94a3b8",
            transform=ax.transAxes)


def _summary_table(ax, rows: list[dict], keys: list[str]):
    """Draw a styled summary table onto an existing Axes."""
    ax.axis("off")
    col_headers = ["Network", "N", "E", "Density (%)", "⟨k⟩", "⟨k²⟩", "κ", "β (%) (disorder)"]
    col_widths   = [0.22, 0.06, 0.06, 0.09, 0.09, 0.09, 0.09, 0.16]

    # Header row
    y_header = 0.94
    x = 0.02
    for header, w in zip(col_headers, col_widths):
        ax.text(x + w / 2, y_header, header, ha="center", va="center",
                fontsize=8.5, fontweight="bold", color="#1e293b",
                transform=ax.transAxes)
        x += w

    ax.plot([0.01, 0.99], [y_header - 0.05, y_header - 0.05],
            color="#94a3b8", linewidth=0.8, transform=ax.transAxes)

    row_h = 0.17
    for row_i, key in enumerate(keys):
        m = rows[key]
        y = y_header - 0.1 - row_i * row_h
        bg_color = "#f8fafc" if key == "all_full" else ("#f0f9ff" if "healthy" in key else "#fff7ed")
        rect = mpatches.FancyBboxPatch(
            (0.01, y - row_h / 2 + 0.005), 0.98, row_h - 0.01,
            boxstyle="square,pad=0", linewidth=0,
            facecolor=bg_color, alpha=0.6, transform=ax.transAxes,
        )
        ax.add_patch(rect)

        cells = [
            NETWORK_LABELS[key].replace("\n", " "),
            str(m["N"]),
            str(m["E"]),
            f"{m['density']*100:.1f}",
            f"{m['k_mean']:.2f}",
            f"{m['k2_mean']:.2f}",
            f"{m['kappa']:.2f}",
            f"{m['disorder']['beta']*100:.2f}",
        ]
        x = 0.02
        for cell, w in zip(cells, col_widths):
            ax.text(x + w / 2, y, cell, ha="center", va="center",
                    fontsize=8, color="#1e293b", transform=ax.transAxes)
            x += w

    ax.plot([0.01, 0.99], [y_header - 0.055, y_header - 0.055],
            color="#94a3b8", linewidth=0.5, transform=ax.transAxes)


def draw_analysis_figure(all_metrics: dict, pairwise: list[dict], out_path: Path):
    keys        = ["all_full", "healthy_predef", "healthy_full", "obese_predef", "obese_full"]
    colors      = [CONDITION_COLORS[k] for k in keys]
    label_short = ["All-full", "H-pre", "H-full", "O-pre", "O-full"]

    # badge background per network
    badge_colors = {
        "all_full":       "#f1f5f9",
        "healthy_predef": "#e0f2fe",
        "healthy_full":   "#e0f2fe",
        "obese_predef":   "#fff7ed",
        "obese_full":     "#fff7ed",
    }

    fig = plt.figure(figsize=(24, 28), dpi=DPI)
    fig.patch.set_facecolor("white")

    # Layout: metrics table | metric bars | strength bar plots | stats table
    gs = gridspec.GridSpec(4, 5, figure=fig, height_ratios=[0.9, 1.1, 1.4, 1.3],
                           hspace=0.65, wspace=0.42)

    # ── Summary table ────────────────────────────────────────────────────────
    ax_table = fig.add_subplot(gs[0, :])
    _summary_table(ax_table, all_metrics, keys)
    ax_table.set_title("Network Metrics Summary", fontsize=13, pad=10,
                        color="#1e293b", fontweight="bold", loc="left",
                        transform=ax_table.transAxes)

    # ── Metric bar plots (5 metrics, each with 5 bars) ────────────────────────
    metric_specs = [
        ("N (nodes)",    [all_metrics[k]["N"]                         for k in keys], "{:.0f}"),
        ("E (edges)",    [all_metrics[k]["E"]                         for k in keys], "{:.0f}"),
        ("Density (%)",  [all_metrics[k]["density"] * 100             for k in keys], "{:.1f}"),
        ("κ (degree heterogeneity)", [all_metrics[k]["kappa"]         for k in keys], "{:.2f}"),
        ("β (%) (disorder)", [all_metrics[k]["disorder"]["beta"] * 100 for k in keys], "{:.2f}"),
    ]
    for col_i, (title, vals, fmt) in enumerate(metric_specs):
        ax = fig.add_subplot(gs[1, col_i])
        ax.set_facecolor("white")
        _metric_bar(ax, vals, label_short, colors, title, "", annotate_fmt=fmt)

    # ── Ranked strength bar plots (organs ordered high → low, fit overlaid) ────
    bar_axes = []
    for col_i, key in enumerate(keys):
        ax = fig.add_subplot(gs[2, col_i])
        ax.set_facecolor("white")
        dis = all_metrics[key]["disorder"]
        _strength_bar_ax(ax, all_metrics[key], dis, CONDITION_COLORS[key])
        ax.set_title(NETWORK_LABELS[key].replace("\n", " — ") + "\nStrength per Organ",
                     fontsize=9.5, color="#1e293b", pad=6)
        ax.text(0.97, 0.97,
                f"β = {dis['beta']*100:.2f}%",
                ha="right", va="top", fontsize=8.5, color="#1e293b",
                transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.3", fc=badge_colors[key], ec="none"))
        bar_axes.append(ax)

    # Sync y-axis across all strength bar plots
    y_max_bar = max(ax.get_ylim()[1] for ax in bar_axes)
    for ax in bar_axes:
        ax.set_ylim(0, y_max_bar)
        if ax is not bar_axes[0]:
            ax.set_ylabel("")

    # ── Pairwise statistical tests ─────────────────────────────────────────────
    ax_stats = fig.add_subplot(gs[3, :])
    ax_stats.set_facecolor("white")
    _stats_table(ax_stats, pairwise)
    ax_stats.set_title("Pairwise Network Comparisons",
                        fontsize=13, pad=10, color="#1e293b", fontweight="bold",
                        loc="left", transform=ax_stats.transAxes)

    fig.suptitle("Metabolic Reference Network — Topology Comparison\n"
                 "No Condition · Pre-defined Cohort Networks · Full Literature Networks",
                 fontsize=15, y=0.995, color="#1e293b", fontweight="bold")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] {out_path.relative_to(HERE)}")


# ── Restricted comparison figure ─────────────────────────────────────────────

def draw_restricted_comparison(all_metrics, h_oo_map, o_oo_map,
                                h_organs, o_organs,
                                h_predef_edges, o_predef_edges,
                                out_path):
    """
    Compare full OO networks restricted to predefined node sets vs predefined networks.
    """
    BETA_PREDEF = {"healthy": 0.13, "obese": 0.17}

    comparisons = {}
    for cond, predef_organs, predef_edges, oo_map in [
        ("healthy", h_organs, h_predef_edges, h_oo_map),
        ("obese",   o_organs, o_predef_edges, o_oo_map),
    ]:
        node_set = set(predef_organs)

        # Restrict full OO edges to predefined nodes
        restricted_edges = [(o1, o2) for o1, o2 in [
            (min(a, b), max(a, b)) for (a, b) in oo_map.keys() if a != b
            ] if o1 in node_set and o2 in node_set]
        # deduplicate
        restricted_edges = sorted(set(restricted_edges))

        # Compute OO-score strength for restricted network
        oo_strength = {n: 0.0 for n in predef_organs}
        for (o1, o2) in restricted_edges:
            score = oo_map.get((o1, o2), oo_map.get((o2, o1), 0.0))
            oo_strength[o1] += score
            oo_strength[o2] += score

        # Build restricted full metrics (unweighted topology)
        restr_weighted = [(o1, o2, oo_map.get((o1, o2), oo_map.get((o2, o1), 0.0)))
                          for o1, o2 in restricted_edges]
        restr_metrics = network_metrics(restricted_edges, predef_organs)
        restr_metrics["strength"] = oo_strength
        restr_metrics["strength_seq"] = sorted([v for v in oo_strength.values() if v > 0], reverse=True)
        restr_metrics["disorder"] = _compute_disorder_edges(restr_weighted, predef_organs)

        # Predefined metrics (recompute with correct node universe)
        predef_metrics = network_metrics(predef_edges, predef_organs)
        predef_metrics["strength"] = {n: float(predef_metrics["degree"].get(n, 0)) for n in predef_organs}
        predef_metrics["strength_seq"] = sorted([v for v in predef_metrics["strength"].values() if v > 0], reverse=True)
        predef_metrics["disorder"] = _compute_disorder_edges(list(predef_edges), predef_organs)

        # Statistical tests
        deg_restr  = [restr_metrics["degree"].get(n, 0) for n in predef_organs]
        deg_predef = [predef_metrics["degree"].get(n, 0) for n in predef_organs]
        t_stat, p_t = scipy_stats.ttest_ind(deg_restr, deg_predef, equal_var=False)

        s_restr  = restr_metrics["strength_seq"]
        s_predef = predef_metrics["strength_seq"]
        if s_restr and s_predef:
            u_stat, p_u = scipy_stats.mannwhitneyu(s_restr, s_predef, alternative="two-sided")
        else:
            u_stat, p_u = float("nan"), float("nan")

        comparisons[cond] = {
            "restr":   restr_metrics,
            "predef":  predef_metrics,
            "organs":  predef_organs,
            "t_stat":  float(t_stat), "p_t": float(p_t),
            "u_stat":  float(u_stat), "p_u": float(p_u),
            "n_restr_edges": len(restricted_edges),
            "n_predef_edges": len(predef_edges),
        }

    COND_COLOR_RESTR = {"healthy": "#0284c7", "obese": "#c2410c"}
    COND_COLOR_PRE   = {"healthy": "#38bdf8", "obese": "#fb923c"}
    DARK = "#1e293b"; MID = "#475569"; LIGHT = "#94a3b8"

    fig = plt.figure(figsize=(20, 14), dpi=DPI)
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.6, 1.0], hspace=0.55, wspace=0.35)

    for col_i, cond in enumerate(("healthy", "obese")):
        c = comparisons[cond]
        ax = fig.add_subplot(gs[0, col_i])
        ax.set_facecolor("white")
        organs = c["organs"]
        # Sort by predefined strength descending
        organs_sorted = sorted(organs, key=lambda n: c["predef"]["strength"].get(n, 0), reverse=True)
        restr_s  = [c["restr"]["strength"].get(n, 0)  for n in organs_sorted]
        predef_s = [c["predef"]["strength"].get(n, 0) for n in organs_sorted]
        labels   = [n.replace(" ", "\n") for n in organs_sorted]
        x        = np.arange(len(labels))
        w = 0.38
        ax.bar(x - w/2, predef_s, w, color=COND_COLOR_PRE[cond],   alpha=0.85,
               label=f"Predefined  (β={c['predef']['disorder']['beta']*100:.2f}%)", zorder=2)
        ax.bar(x + w/2, restr_s,  w, color=COND_COLOR_RESTR[cond], alpha=0.85,
               label=f"Full OO (restricted)  (β={c['restr']['disorder']['beta']*100:.2f}%)", zorder=2)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=7, rotation=45, ha="right")
        ax.set_ylabel("Strength", fontsize=9, color=DARK)
        ax.tick_params(colors=MID)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"): ax.spines[sp].set_color(LIGHT)
        ax.grid(axis="y", color="#e2e8f0", linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)
        cond_label = cond.capitalize()
        ax.set_title(
            f"{cond_label}: Full OO (restricted) vs Pre-defined\n"
            f"Nodes: {len(organs)}  ·  Pre-defined edges: {c['n_predef_edges']}  ·  "
            f"Full OO restricted edges: {c['n_restr_edges']}",
            fontsize=10, color=DARK, pad=6, fontweight="bold",
        )
        ax.legend(fontsize=8, framealpha=0.8)

    # Stats table (bottom row, spanning both columns)
    ax_stats = fig.add_subplot(gs[1, :])
    ax_stats.axis("off")

    col_headers = [
        "Comparison", "Nodes",
        "Predef edges", "Full-OO edges (restricted)",
        "Predef density", "Full-OO density", "Δ Density", "t-stat", "p (t-test)", "sig",
        "Predef β", "Full-OO β", "Δ β", "U-stat", "p (MWU)", "sig",
    ]
    col_w = [0.11, 0.04, 0.07, 0.11, 0.07, 0.07, 0.06, 0.05, 0.07, 0.04,
             0.05, 0.05, 0.05, 0.05, 0.07, 0.04]

    y0 = 0.92; x0 = 0.01
    for h, w in zip(col_headers, col_w):
        ax_stats.text(x0 + w/2, y0, h, ha="center", va="center",
                      fontsize=7.5, fontweight="bold", color=DARK, transform=ax_stats.transAxes)
        x0 += w

    ax_stats.plot([0.01, 0.99], [y0 - 0.06, y0 - 0.06], color="#94a3b8", lw=0.8,
                  transform=ax_stats.transAxes)

    sig_colors = {"***": "#dc2626", "**": "#ea580c", "*": "#ca8a04", "ns": "#64748b", "—": "#94a3b8"}

    for row_i, cond in enumerate(("healthy", "obese")):
        c = comparisons[cond]
        rm = c["restr"]; pm = c["predef"]
        y = y0 - 0.18 - row_i * 0.22
        bg = "#f0f9ff" if cond == "healthy" else "#fff7ed"
        ax_stats.add_patch(mpatches.FancyBboxPatch(
            (0.01, y - 0.09), 0.98, 0.18, boxstyle="square,pad=0",
            linewidth=0, facecolor=bg, alpha=0.6, transform=ax_stats.transAxes))

        delta_d = rm["density"] - pm["density"]
        delta_b = rm["disorder"]["beta"] - pm["disorder"]["beta"]
        sig_t = _sig_label(c["p_t"]); sig_u = _sig_label(c["p_u"])

        cells = [
            f"{cond.capitalize()} full-OO vs predef",
            str(len(c["organs"])),
            str(c["n_predef_edges"]), str(c["n_restr_edges"]),
            f"{pm['density']*100:.1f}%", f"{rm['density']*100:.1f}%", f"{delta_d*100:+.1f}%",
            f"{c['t_stat']:.2f}",
            f"{c['p_t']:.3f}" if c["p_t"] >= 0.001 else "<0.001",
            sig_t,
            f"{pm['disorder']['beta']*100:.2f}%", f"{rm['disorder']['beta']*100:.2f}%", f"{delta_b*100:+.2f}%",
            f"{int(c['u_stat'])}" if not np.isnan(c["u_stat"]) else "—",
            f"{c['p_u']:.3f}" if not np.isnan(c["p_u"]) and c["p_u"] >= 0.001 else ("<0.001" if not np.isnan(c["p_u"]) else "—"),
            sig_u,
        ]
        x0 = 0.01
        for ci, (cell, w) in enumerate(zip(cells, col_w)):
            color = DARK; fw = "normal"
            if ci == 9:  color = sig_colors.get(cell, DARK); fw = "bold" if cell not in ("ns","—") else "normal"
            if ci == 15: color = sig_colors.get(cell, DARK); fw = "bold" if cell not in ("ns","—") else "normal"
            ax_stats.text(x0 + w/2, y, cell, ha="center", va="center",
                          fontsize=7.5, color=color, fontweight=fw, transform=ax_stats.transAxes)
            x0 += w

    ax_stats.set_title("Restricted Comparison: Full OO Network (restricted to predefined nodes) vs Pre-defined Network",
                        fontsize=11, pad=8, color=DARK, fontweight="bold", loc="left",
                        transform=ax_stats.transAxes)
    ax_stats.text(0.99, -0.02, "* p<0.05   ** p<0.01   *** p<0.001   ns not significant",
                  ha="right", va="bottom", fontsize=7.5, color="#94a3b8",
                  transform=ax_stats.transAxes)

    fig.suptitle(
        "Metabolic Reference Network — Restricted Node Comparison\n"
        "Full OO Network (restricted to predefined nodes) vs Pre-defined Cohort Networks",
        fontsize=13, y=0.995, color=DARK, fontweight="bold",
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] {out_path.relative_to(HERE)}")


def _compute_lit_vs_pc(lit_data: dict, pc_data: dict, global_paper_max: float):
    """
    Density and disorder for the literature networks (OO-weighted and
    normalised paper count) and the partial-correlation (|r|) reference
    network of one condition, all restricted to the organs of the PC network
    so that every scheme covers the same organ pairs.

    lit_data: {"oo_edges": [(o1,o2,w)], "paper_edges": [(o1,o2,count)], ...}
    pc_data:  {"pc_edges": [(o1,o2,|r|)], "nodes": [organs]}

    Density = sum of edge weights / number of possible pairs.
    Disorder = β of a·exp(-β·x)+c fitted on the per-organ strengths sorted by
    that scheme's own strength (descending).
    """
    organs = list(pc_data["nodes"])
    oset   = set(organs)
    pairs  = [(a, b) for i, a in enumerate(organs) for b in organs[i + 1:]]
    max_edges = len(pairs)

    def restrict(edges, scale=1.0):
        m = {}
        for e in edges:
            if len(e) == 3 and e[0] in oset and e[1] in oset:
                m[frozenset((e[0], e[1]))] = e[2] * scale
        return m

    schemes = {
        "lit_oo":   restrict(lit_data["oo_edges"]),
        "lit_norm": restrict(lit_data["paper_edges"], 1.0 / global_paper_max),
        "pc":       restrict(pc_data["pc_edges"]),
    }

    out = {}
    for key, wmap in schemes.items():
        vec = [wmap.get(frozenset(p), 0.0) for p in pairs]
        strength = {o: 0.0 for o in organs}
        for k, w in wmap.items():
            a, b = tuple(k)
            strength[a] += w
            strength[b] += w
        order = sorted(organs, key=lambda o: strength[o], reverse=True)
        vals  = [strength[o] for o in order]
        beta, _, _, r2, xf, yf = _fit_to_sequence(vals)
        out[key] = {
            "n_edges": len(wmap), "density": sum(vec) / max_edges, "vec": vec,
            "strength": strength, "order": order, "vals": vals,
            "beta": beta, "r2": r2, "xf": xf, "yf": yf,
        }
    return out, max_edges


def draw_pc_strength_figure(pc_metrics: dict, out_path: Path) -> None:
    """
    Per-organ strength distribution (|r| weights, p < 0.05 edges) of the
    healthy and obese partial-correlation reference networks, organs sorted by
    their own strength, with the exponential disorder fit.

    pc_metrics: {"healthy": metrics["pc"], "obese": metrics["pc"]} where each
    value comes from _compute_lit_vs_pc().
    """
    COND_LABELS = {"healthy": "Healthy", "obese": "Obese"}
    C_PC = {"healthy": "#059669", "obese": "#7c3aed"}
    DARK = "#1e293b"; MID = "#475569"; LIGHT = "#94a3b8"

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.2), dpi=DPI,
                             gridspec_kw={"wspace": 0.22})
    fig.patch.set_facecolor("white")

    for ax, cond in zip(axes, ("healthy", "obese")):
        m = pc_metrics[cond]
        clr = C_PC[cond]
        ax.set_facecolor("white")
        for sp in ("top", "right"):   ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"): ax.spines[sp].set_color(LIGHT)
        ax.grid(axis="y", color="#e2e8f0", linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)

        x = np.arange(len(m["order"]))
        ax.bar(x, m["vals"], color=clr, alpha=0.85, zorder=2, width=0.70)
        ax.plot(m["xf"], m["yf"], "--", color=DARK, linewidth=1.8, alpha=0.75, zorder=3)
        ax.set_xticks(x)
        ax.set_xticklabels([n.replace(" ", "\n") for n in m["order"]], fontsize=7,
                           rotation=45, ha="right", color=DARK, multialignment="center")
        ax.tick_params(axis="y", labelsize=8, colors=DARK)
        ax.tick_params(axis="x", length=0)
        ax.set_xlim(-0.6, len(m["order"]) - 0.4)
        ax.set_ylabel("Strength (sum of |r| per organ)", fontsize=9, color=MID)
        ax.set_title(
            f"{COND_LABELS[cond]} PC reference  ·  {m['n_edges']} edges  ·  "
            f"density {m['density']*100:.1f}%\nβ = {m['beta']*100:.2f}%  (R²={m['r2']:.3f})",
            fontsize=10, color=clr, pad=6, fontweight="bold",
        )

    fig.suptitle("Strength Distribution and Disorder of the Partial-Correlation Reference Networks",
                 fontsize=12.5, y=1.02, color=DARK, fontweight="bold")
    fig.text(0.99, -0.08,
             "Edges: partial correlations with p < 0.05 (uncorrected), weight = |r|  ·  "
             "organs sorted by own strength (descending)  ·  dashed line: a·exp(−β·x)+c fit",
             ha="right", va="bottom", fontsize=7, color=LIGHT)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] {out_path.relative_to(HERE)}")


def draw_lit_strength_figure(lit_metrics: dict, out_path: Path) -> None:
    """
    Per-organ strength distribution of the literature networks (restricted to the
    PC organs) of the healthy and obese condition, with the exponential disorder
    fit a·exp(-β·x)+c -- the literature counterpart of draw_pc_strength_figure().

    lit_metrics: {cond: metrics} where metrics is the first element returned by
    _compute_lit_vs_pc(); rows = condition, columns = OO-weighted / normalised
    paper count.
    """
    COND_LABELS = {"healthy": "Healthy", "obese": "Obese"}
    C_LIT = {"healthy": "#0284c7", "obese": "#c2410c"}
    SCHEMES = (("lit_oo", "OO-weighted", "Strength (sum of OO per organ)"),
               ("lit_norm", "normalised paper count", "Strength (sum of normalised paper counts per organ)"))
    DARK = "#1e293b"; MID = "#475569"; LIGHT = "#94a3b8"

    fig, axes = plt.subplots(2, 2, figsize=(15, 10), dpi=DPI,
                             gridspec_kw={"wspace": 0.22, "hspace": 0.55})
    fig.patch.set_facecolor("white")

    for ri, cond in enumerate(("healthy", "obese")):
        for ci, (key, scheme_label, ylabel) in enumerate(SCHEMES):
            ax = axes[ri, ci]
            m = lit_metrics[cond][key]
            clr = C_LIT[cond]
            ax.set_facecolor("white")
            for sp in ("top", "right"):   ax.spines[sp].set_visible(False)
            for sp in ("left", "bottom"): ax.spines[sp].set_color(LIGHT)
            ax.grid(axis="y", color="#e2e8f0", linewidth=0.7, zorder=0)
            ax.set_axisbelow(True)

            x = np.arange(len(m["order"]))
            ax.bar(x, m["vals"], color=clr, alpha=0.85 if ci == 0 else 0.55, zorder=2, width=0.70)
            ax.plot(m["xf"], m["yf"], "--", color=DARK, linewidth=1.8, alpha=0.75, zorder=3)
            ax.set_xticks(x)
            ax.set_xticklabels([n.replace(" ", "\n") for n in m["order"]], fontsize=7,
                               rotation=45, ha="right", color=DARK, multialignment="center")
            ax.tick_params(axis="y", labelsize=8, colors=DARK)
            ax.tick_params(axis="x", length=0)
            ax.set_xlim(-0.6, len(m["order"]) - 0.4)
            ax.set_ylabel(ylabel, fontsize=9, color=MID)
            ax.set_title(
                f"{COND_LABELS[cond]} literature ({scheme_label})  ·  {m['n_edges']} edges  ·  "
                f"density {m['density']*100:.1f}%\nβ = {m['beta']*100:.2f}%  (R²={m['r2']:.3f})",
                fontsize=10, color=clr, pad=6, fontweight="bold",
            )

    fig.suptitle("Strength Distribution and Disorder of the Literature Networks "
                 "(restricted to the organs of the partial-correlation networks)",
                 fontsize=12.5, y=0.98, color=DARK, fontweight="bold")
    fig.text(0.99, 0.02,
             "Literature network rebuilt from the PC organs only (Colon excluded; OO and normalisation recomputed)  ·  "
             "organs sorted by own strength (descending)  ·  dashed line: a·exp(−β·x)+c fit",
             ha="right", va="bottom", fontsize=7, color=LIGHT)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] {out_path.relative_to(HERE)}")


def compute_matrix_correlations(lit_restricted: dict, pc_full_abs: dict, global_paper_max: float,
                                drop_zero: bool = False):
    """
    Direct matrix comparison of the literature networks (OO-weighted and
    normalised paper count, restricted to the PC organs) with the full
    partial-correlation matrix (|r|, all pairs) of the same condition.

    Only the upper triangle (without diagonal) is used -- the networks are
    undirected, so the lower triangle holds no extra information.
    Returns (rows, vectors): rows = list of result dicts, vectors =
    {(cond, scheme): (pair_labels, x_literature, y_pc)}.
    """
    rows, vectors = [], {}
    for cond in ("healthy", "obese"):
        pc_df = pc_full_abs[cond]
        organs = list(pc_df.columns)
        pairs = [(a, b) for i, a in enumerate(organs) for b in organs[i + 1:]]
        y = np.array([pc_df.loc[a, b] for a, b in pairs], dtype=float)
        oset = set(organs)
        for scheme, key, scale in (("OO-weighted", "oo_edges", 1.0),
                                   ("normalised paper count", "paper_edges", 1.0 / global_paper_max)):
            wmap = {frozenset((e[0], e[1])): e[2] * scale
                    for e in lit_restricted[cond][key] if e[0] in oset and e[1] in oset}
            x = np.array([wmap.get(frozenset(p), 0.0) for p in pairs], dtype=float)
            n_all = len(pairs)
            if drop_zero:
                # keep only pairs that have an edge (weight > 0) in BOTH networks
                keep = (x > 0) & (y > 0)
                pairs_k = [p for p, k in zip(pairs, keep) if k]
                x_k, y_k = x[keep], y[keep]
            else:
                pairs_k, x_k, y_k = pairs, x, y
            r, p_r = scipy_stats.pearsonr(x_k, y_k)
            rho, p_rho = scipy_stats.spearmanr(x_k, y_k)
            rows.append({
                "variant": "zero-weight pairs removed" if drop_zero else "all pairs",
                "condition": cond, "literature_weighting": scheme, "n_organs": len(organs),
                "n_pairs_total": n_all, "n_pairs": len(pairs_k),
                "n_literature_edges": int((x > 0).sum()),
                "pearson_r": float(r), "pearson_p": float(p_r),
                "spearman_rho": float(rho), "spearman_p": float(p_rho),
            })
            vectors[(cond, scheme)] = (pairs_k, x_k, y_k)
    return rows, vectors


def draw_matrix_correlation_figure(rows: list, vectors: dict, out_path: Path,
                                   drop_zero: bool = False) -> None:
    """2x2 scatter: literature weight (x) vs |partial correlation| (y) per organ pair."""
    C_LIT = {"healthy": "#0284c7", "obese": "#c2410c"}
    LABEL = {"healthy": "Healthy", "obese": "Obese"}
    DARK = "#1e293b"; MID = "#475569"; LIGHT = "#94a3b8"
    schemes = ("OO-weighted", "normalised paper count")
    res = {(r["condition"], r["literature_weighting"]): r for r in rows}

    def _p(p):
        return "p < 0.001" if p < 0.001 else f"p = {p:.3f}"

    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=DPI, gridspec_kw={"wspace": 0.22, "hspace": 0.32})
    fig.patch.set_facecolor("white")
    for ri, cond in enumerate(("healthy", "obese")):
        for ci, scheme in enumerate(schemes):
            ax = axes[ri, ci]
            pairs, x, y = vectors[(cond, scheme)]
            r = res[(cond, scheme)]
            clr = C_LIT[cond]
            for sp in ("top", "right"):   ax.spines[sp].set_visible(False)
            for sp in ("left", "bottom"): ax.spines[sp].set_color(LIGHT)
            ax.grid(color="#e2e8f0", linewidth=0.7, zorder=0)
            ax.set_axisbelow(True)
            ax.scatter(x, y, s=34, color=clr, alpha=0.75, edgecolor="white", linewidth=0.5, zorder=3)
            if x.std() > 0:
                k, b = np.polyfit(x, y, 1)
                xs = np.linspace(x.min(), x.max(), 50)
                ax.plot(xs, k * xs + b, "--", color=DARK, linewidth=1.5, alpha=0.8, zorder=2)
            # label the pairs that stand out on either axis
            idx = set(np.argsort(-y)[:3]) | set(np.argsort(-x)[:3])
            for i in idx:
                ax.annotate(f"{pairs[i][0]}–{pairs[i][1]}", (x[i], y[i]), xytext=(4, 4),
                            textcoords="offset points", fontsize=6.5, color=MID)
            ax.set_xlabel(f"Literature weight ({scheme})", fontsize=9, color=MID)
            ax.set_ylabel("|partial correlation| (PET)", fontsize=9, color=MID)
            ax.tick_params(labelsize=8, colors=DARK)
            ax.set_title(f"{LABEL[cond]} — {scheme}  (n = {r['n_pairs']} pairs)\n"
                         f"Pearson r = {r['pearson_r']:+.2f} ({_p(r['pearson_p'])})  ·  "
                         f"Spearman ρ = {r['spearman_rho']:+.2f} ({_p(r['spearman_p'])})",
                         fontsize=10, color=clr, pad=6, fontweight="bold")
    n = rows[0]["n_pairs_total"]
    suffix = " (Pairs Without a Literature Edge Removed)" if drop_zero else ""
    fig.suptitle("Literature-Based vs Partial-Correlation (PET) Networks: Pair-wise Correlation of the Network Matrices"
                 + suffix, fontsize=12.5, y=0.97, color=DARK, fontweight="bold")
    zero_note = ("Pairs with weight 0 in either network are removed (here: pairs without a literature edge)"
                 if drop_zero else "Pairs with no literature edge have weight 0")
    fig.text(0.99, 0.045,
             f"Upper triangle of the {rows[0]['n_organs']}×{rows[0]['n_organs']} matrices ({n} organ pairs, undirected)  ·  "
             "PC: full partial-correlation matrix, absolute values, all pairs  ·  literature rebuilt from the PC organs only "
             f"(Colon excluded)\n{zero_note}  ·  the pairs are not independent "
             "(each organ appears in 12 pairs), so p-values are approximate",
             ha="right", va="bottom", fontsize=7, color=LIGHT)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] {out_path.relative_to(HERE)}")


def draw_lit_vs_pc_table(cond: str, metrics: dict, max_edges: int,
                         global_paper_max: float, out_path: Path) -> None:
    """
    Table comparing the literature-based network (OO-weighted and normalised
    paper count) with the partial-correlation reference of one condition, in
    density and disorder, with statistics.

    Section 1: edges, density, disorder β and R² per network.
    Section 2: literature − PC differences with
      · paired t-test on the per-pair weight vectors (density),
      · Mann-Whitney U on the positive per-organ strengths (disorder).
    """
    label = {"healthy": "Healthy", "obese": "Obese"}[cond]
    DARK = "#1e293b"; MID = "#475569"; LIGHT = "#94a3b8"
    col_lit = {"healthy": "#0284c7", "obese": "#c2410c"}[cond]
    col_pc  = {"healthy": "#059669", "obese": "#7c3aed"}[cond]
    row_bg  = {"lit_oo": "#f0f9ff" if cond == "healthy" else "#fff7ed",
               "lit_norm": "#f8fafc", "pc": "#f0fdf4" if cond == "healthy" else "#faf5ff"}
    sig_col = {"***": "#dc2626", "**": "#ea580c", "*": "#ca8a04", "ns": MID, "—": LIGHT}

    names = {"lit_oo": f"Literature {label} (OO-weighted)",
             "lit_norm": f"Literature {label} (norm. paper count)",
             "pc": f"{label} PC reference (|r|)"}
    colors = {"lit_oo": col_lit, "lit_norm": MID, "pc": col_pc}

    def _p_str(p):
        if np.isnan(p): return "—"
        return f"{p:.3f}" if p >= 0.001 else "<0.001"

    def _delta_clr(v):
        return "#059669" if v > 0 else ("#dc2626" if v < 0 else MID)

    comps = []
    for lit_key in ("lit_oo", "lit_norm"):
        a, b = metrics[lit_key], metrics["pc"]
        t, p_t = scipy_stats.ttest_rel(a["vec"], b["vec"])
        pos_a = [v for v in a["vals"] if v > 0]
        pos_b = [v for v in b["vals"] if v > 0]
        if pos_a and pos_b:
            u, p_u = scipy_stats.mannwhitneyu(pos_a, pos_b, alternative="two-sided")
        else:
            u, p_u = float("nan"), float("nan")
        comps.append({
            "label": f"{names[lit_key].split(' (')[0]} ({'OO' if lit_key == 'lit_oo' else 'norm. paper'}) vs PC",
            "dd": a["density"] - b["density"], "t": float(t), "p_t": float(p_t),
            "db": a["beta"] - b["beta"], "u": float(u), "p_u": float(p_u),
        })

    fig, axes = plt.subplots(2, 1, figsize=(13, 5.6), dpi=DPI,
                             gridspec_kw={"height_ratios": [1, 0.8], "hspace": 0.55})
    fig.patch.set_facecolor("white")

    # ── Section 1 ────────────────────────────────────────────────────────────
    ax1 = axes[0]; ax1.axis("off")
    ax1.set_title(f"Section 1 — {label}: density and disorder per network "
                  f"({len(metrics['pc']['order'])} shared organs, {max_edges} possible pairs)",
                  fontsize=9, color=DARK, fontweight="bold", loc="left", pad=4)
    hdrs1  = ["Network", "Edges", "Density (%)", "Disorder β (%)", "R² (fit)"]
    col_x1 = [0.00, 0.40, 0.52, 0.68, 0.85]
    col_w1 = [0.40, 0.12, 0.16, 0.17, 0.15]
    y0 = 0.84
    for h, cx, cw in zip(hdrs1, col_x1, col_w1):
        ax1.text(cx + cw/2, y0, h, ha="center", va="center", fontsize=8.5,
                 fontweight="bold", color=DARK, transform=ax1.transAxes)
    ax1.plot([0.01, 0.99], [y0 - 0.13, y0 - 0.13], color=DARK, lw=0.8, transform=ax1.transAxes)
    for ri, key in enumerate(("lit_oo", "lit_norm", "pc")):
        m = metrics[key]
        y = y0 - 0.28 - ri * 0.28
        ax1.add_patch(mpatches.FancyBboxPatch(
            (0.01, y - 0.12), 0.98, 0.24, boxstyle="square,pad=0",
            linewidth=0, facecolor=row_bg[key], alpha=0.7, transform=ax1.transAxes,
            clip_on=False))
        cells = [names[key], str(m["n_edges"]), f"{m['density']*100:.1f}%",
                 f"{m['beta']*100:.2f}%", f"{m['r2']:.3f}"]
        clrs  = [colors[key], DARK, colors[key], colors[key], DARK]
        for cell, cx, cw, clr in zip(cells, col_x1, col_w1, clrs):
            ax1.text(cx + cw/2, y, cell, ha="center", va="center", fontsize=9, color=clr,
                     fontweight="bold" if clr not in (DARK, MID) else "normal",
                     transform=ax1.transAxes)

    # ── Section 2 ────────────────────────────────────────────────────────────
    ax2 = axes[1]; ax2.axis("off")
    ax2.set_title("Section 2 — Literature minus PC reference  "
                  "(density: paired t-test on per-pair weights · disorder: Mann-Whitney U on organ strengths)",
                  fontsize=9, color=DARK, fontweight="bold", loc="left", pad=4)
    hdrs2  = ["Comparison", "Δ density\n(%)", "t", "p (t-test)", "sig",
              "Δβ\n(%)", "U (MWU)", "p (MWU)", "sig"]
    col_x2 = [0.00, 0.31, 0.41, 0.49, 0.60, 0.67, 0.77, 0.86, 0.94]
    col_w2 = [0.31, 0.10, 0.08, 0.11, 0.07, 0.10, 0.09, 0.08, 0.06]
    y0b = 0.84
    for h, cx, cw in zip(hdrs2, col_x2, col_w2):
        ax2.text(cx + cw/2, y0b, h, ha="center", va="center", fontsize=8.5,
                 fontweight="bold", color=DARK, transform=ax2.transAxes,
                 multialignment="center")
    ax2.plot([0.01, 0.99], [y0b - 0.16, y0b - 0.16], color=DARK, lw=0.8, transform=ax2.transAxes)
    for ri, c in enumerate(comps):
        y = y0b - 0.34 - ri * 0.30
        ax2.add_patch(mpatches.FancyBboxPatch(
            (0.01, y - 0.13), 0.98, 0.26, boxstyle="square,pad=0",
            linewidth=0, facecolor=["#f1f5f9", "#e9f0f7"][ri], alpha=0.6, transform=ax2.transAxes))
        s_t, s_u = _sig_label(c["p_t"]), _sig_label(c["p_u"])
        cells = [c["label"], f"{c['dd']*100:+.1f}%", f"{c['t']:+.2f}", _p_str(c["p_t"]), s_t,
                 f"{c['db']*100:+.2f}%", "—" if np.isnan(c["u"]) else str(int(c["u"])),
                 _p_str(c["p_u"]), s_u]
        clrs = [DARK, _delta_clr(c["dd"]), DARK, DARK, sig_col.get(s_t, DARK),
                _delta_clr(c["db"]), DARK, DARK, sig_col.get(s_u, DARK)]
        for cell, cx, cw, clr in zip(cells, col_x2, col_w2, clrs):
            ax2.text(cx + cw/2, y, cell, ha="center", va="center", fontsize=8.5, color=clr,
                     fontweight="bold" if clr not in (DARK, MID) else "normal",
                     transform=ax2.transAxes)

    fig.text(0.99, 0.0,
             "Literature network rebuilt from the PC organs only (Colon excluded; OO and normalisation recomputed)  ·  "
             "PC: partial correlations, p < 0.05 uncorrected, weight |r|\n"
             f"norm. paper count = count ÷ max ({global_paper_max:.0f})  ·  "
             "weights are on different scales (OO, normalised count, |r|)  ·  "
             "Δ = literature − PC  ·  * p<0.05  ** p<0.01  *** p<0.001",
             ha="right", va="bottom", fontsize=6.8, color=LIGHT)
    fig.suptitle(f"{label}: Literature-Based Network vs Partial-Correlation Reference — Density and Disorder",
                 fontsize=12, y=1.02, color=DARK, fontweight="bold")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] {out_path.relative_to(HERE)}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ref_dir = HERE / "reference_network_only_metabolic"

    print("[i] Loading pre-defined networks from CSV …")
    _, all_organs             = edges_from_csv(ref_dir / "all_organ_connections.csv")
    h_predef_edges, h_organs  = edges_from_csv(ref_dir / "healthy_cohort_connections.csv")
    o_predef_edges, o_organs  = edges_from_csv(ref_dir / "obese_cohort_connections.csv")

    print("[i] Computing full condition networks …")
    all_full_edges, all_oo_map = full_condition_network("all",     all_organs)
    h_full_edges,   h_oo_map   = full_condition_network("healthy", all_organs)
    o_full_edges,   o_oo_map   = full_condition_network("obese",   all_organs)
    print(f"    All     full: {len(all_full_edges)} edges")
    print(f"    Healthy full: {len(h_full_edges)} edges")
    print(f"    Obese   full: {len(o_full_edges)} edges")

    print("[i] Computing unweighted metrics …")
    all_metrics = {
        "all_full":       network_metrics(all_full_edges,  all_organs),
        "healthy_full":   network_metrics(h_full_edges,    all_organs),
        "obese_full":     network_metrics(o_full_edges,    all_organs),
        "healthy_predef": network_metrics(h_predef_edges,  h_organs),
        "obese_predef":   network_metrics(o_predef_edges,  o_organs),
    }

    # Add OO-score strength for full networks
    for key, oo_map in [("all_full", all_oo_map), ("healthy_full", h_oo_map), ("obese_full", o_oo_map)]:
        m = all_metrics[key]
        oo_strength = {n: 0.0 for n in m["nodes"]}
        for (o1, o2), score in oo_map.items():
            if o1 in oo_strength: oo_strength[o1] += score
            if o2 in oo_strength: oo_strength[o2] += score
        # Each pair added twice (o1,o2) and (o2,o1), so halve
        oo_strength = {n: v / 2 for n, v in oo_strength.items()}
        m["strength"] = oo_strength
        m["strength_seq"] = sorted([v for v in oo_strength.values() if v > 0], reverse=True)

    # Predefined: strength = degree (unweighted)
    for key in ("healthy_predef", "obese_predef"):
        m = all_metrics[key]
        m["strength"] = {n: float(m["degree"].get(n, 0)) for n in m["nodes"]}
        m["strength_seq"] = sorted([v for v in m["strength"].values() if v > 0], reverse=True)

    print("[i] Computing weighted density and disorder …")
    # Build weighted edge lists for each network
    _all_full_w = [(o1, o2, all_oo_map[(o1, o2)]) for o1, o2 in all_full_edges]
    _h_full_w   = [(o1, o2, h_oo_map[(o1, o2)])   for o1, o2 in h_full_edges]
    _o_full_w   = [(o1, o2, o_oo_map[(o1, o2)])   for o1, o2 in o_full_edges]
    _h_pre_w    = list(h_predef_edges)
    _o_pre_w    = list(o_predef_edges)

    _full_edge_map = {
        "all_full":       (_all_full_w, all_organs),
        "healthy_full":   (_h_full_w,   all_organs),
        "obese_full":     (_o_full_w,   all_organs),
        "healthy_predef": (_h_pre_w,    h_organs),
        "obese_predef":   (_o_pre_w,    o_organs),
    }

    for key, (ew, nodes) in _full_edge_map.items():
        m = all_metrics[key]
        m["disorder"] = _compute_disorder_edges(ew, nodes)

    for key, m in all_metrics.items():
        dis = m["disorder"]
        print(f"  {NETWORK_LABELS[key].replace(chr(10), ' '):38s}  "
              f"N={m['N']:2d}  E={m['E']:2d}  "
              f"density={m['density']*100:.1f}%  "
              f"<k>={m['k_mean']:.2f}  <k²>={m['k2_mean']:.2f}  κ={m['kappa']:.2f}  "
              f"β={dis['beta']*100:.2f}%  R²={dis['r2']:.3f}")

    keys = ["all_full", "healthy_predef", "healthy_full", "obese_predef", "obese_full"]
    print("[i] Computing pairwise statistical tests …")
    pairwise = compute_pairwise_stats(all_metrics, keys)
    for r in pairwise:
        print(f"  {r['label']:22s}  "
              f"t={r['t_stat']:+6.2f} p_t={r['p_t']:.3f} ({_sig_label(r['p_t'])})  "
              f"U={r['u_stat']:.0f} p_u={r['p_u']:.3f} ({_sig_label(r['p_u'])})")

    print("[i] Drawing main figure …")
    draw_analysis_figure(all_metrics, pairwise, OUT_DIR / "network_analysis.png")

    print("[i] Computing paper-count networks …")
    all_pc_edges, all_pc_map = full_condition_network_papercounts("all",     all_organs)
    h_pc_edges,   h_pc_map   = full_condition_network_papercounts("healthy", all_organs)
    o_pc_edges,   o_pc_map   = full_condition_network_papercounts("obese",   all_organs)

    pc_conditions = {
        "all": {
            "paper_edges": [(o1, o2, all_pc_map[(o1, o2)]) for o1, o2 in all_pc_edges],
            "oo_edges":    _all_full_w,
            "nodes":       all_organs,
        },
        "healthy": {
            "paper_edges": [(o1, o2, h_pc_map[(o1, o2)]) for o1, o2 in h_pc_edges],
            "oo_edges":    _h_full_w,
            "nodes":       all_organs,
        },
        "obese": {
            "paper_edges": [(o1, o2, o_pc_map[(o1, o2)]) for o1, o2 in o_pc_edges],
            "oo_edges":    _o_full_w,
            "nodes":       all_organs,
        },
    }
    print("[i] Drawing paper-count vs OO-coefficient comparison figure …")
    draw_papercount_vs_oo_figure(
        pc_conditions,
        OUT_DIR / "network_analysis_papercount_vs_oo.png",
    )

    print("[i] Drawing disorder comparison table …")
    draw_disorder_table(
        pc_conditions,
        OUT_DIR / "network_analysis_disorder_table.png",
    )

    print("[i] Drawing density comparison table …")
    draw_density_table(
        pc_conditions,
        OUT_DIR / "network_analysis_density_table.png",
    )

    print("[i] Drawing disorder-types comparison table …")
    draw_disorder_types_table(
        pc_conditions,
        OUT_DIR / "network_analysis_disorder_types_table.png",
    )

    print("[i] Drawing disorder-types 3×3 figure …")
    draw_disorder_types_figure(
        pc_conditions,
        OUT_DIR / "network_analysis_disorder_types_figure.png",
    )

    print("[i] Drawing restricted comparison figure …")
    draw_restricted_comparison(
        all_metrics, h_oo_map, o_oo_map,
        h_organs, o_organs,
        h_predef_edges, o_predef_edges,
        OUT_DIR / "network_analysis_restricted.png",
    )

    print("[i] Loading PC reference networks …")
    h_pc_edges, h_pc_organs = weighted_edges_from_pc_csv(
        ref_dir / "healthy_cohort_PC_connections.csv",
        ref_dir / "healthy_cohort_PC_significance.csv",
    )
    o_pc_edges, o_pc_organs = weighted_edges_from_pc_csv(
        ref_dir / "obese_cohort_PC_connections.csv",
        ref_dir / "obese_cohort_PC_significance.csv",
    )
    print(f"    Healthy PC: {len(h_pc_edges)} edges, {len(h_pc_organs)} organs")
    print(f"    Obese   PC: {len(o_pc_edges)} edges, {len(o_pc_organs)} organs")

    pc_ref_conditions = {
        "healthy": {"pc_edges": h_pc_edges, "nodes": h_pc_organs},
        "obese":   {"pc_edges": o_pc_edges, "nodes": o_pc_organs},
    }

    print("[i] Computing density and disorder: literature vs PC reference …")
    # The literature networks are rebuilt from the PC network's organs only:
    # pairs with other organs are dropped BEFORE OO (organ totals) and the
    # paper-count normalisation (global max) are computed, so neither depends
    # on organs outside the PC network.
    lit_restricted = {}
    for cond in ("healthy", "obese"):
        organs_c = list(pc_ref_conditions[cond]["nodes"])
        oo_edges_c, oo_map_c = full_condition_network(cond, all_organs, restrict_to=organs_c)
        paper_edges_c, paper_map_c = full_condition_network_papercounts(cond, all_organs, restrict_to=organs_c)
        lit_restricted[cond] = {
            "oo_edges":    [(o1, o2, oo_map_c[(o1, o2)]) for o1, o2 in oo_edges_c],
            "paper_edges": [(o1, o2, paper_map_c[(o1, o2)]) for o1, o2 in paper_edges_c],
            "nodes":       organs_c,
        }
        print(f"    {cond}: {len(organs_c)} PC organs -> {len(oo_edges_c)} literature edges among them")
    # shared normalisation constant for both conditions: max pair count in the
    # restricted networks (all-papers network restricted to the same organs)
    _, all_map_r = full_condition_network_papercounts("all", all_organs,
                                                      restrict_to=list(pc_ref_conditions["healthy"]["nodes"]))
    pc_global_max = max(
        [v for v in all_map_r.values()]
        + [e[2] for c in ("healthy", "obese") for e in lit_restricted[c]["paper_edges"]],
        default=1.0,
    )
    lit_vs_pc = {}
    for cond in ("healthy", "obese"):
        lit_vs_pc[cond] = _compute_lit_vs_pc(lit_restricted[cond], pc_ref_conditions[cond], pc_global_max)

    print("[i] Drawing PC reference strength distribution figure …")
    draw_pc_strength_figure(
        {c: lit_vs_pc[c][0]["pc"] for c in ("healthy", "obese")},
        OUT_DIR / "network_analysis_pc_strength_distribution.png",
    )
    print("[i] Correlating literature and PC network matrices (upper triangle) …")
    pc_full_abs = {c: pd.read_csv(ref_dir / f"{c}_cohort_PC_connections.csv", index_col=0).abs()
                   for c in ("healthy", "obese")}
    stat_dir = HERE / "statistics" / "network_analysis"
    stat_dir.mkdir(parents=True, exist_ok=True)
    all_corr_rows, pair_tbl = [], []
    for drop_zero, fig_name in ((False, "network_analysis_lit_vs_pc_correlation.png"),
                                (True, "network_analysis_lit_vs_pc_correlation_nonzero.png")):
        corr_rows, corr_vectors = compute_matrix_correlations(lit_restricted, pc_full_abs, pc_global_max,
                                                              drop_zero=drop_zero)
        all_corr_rows += corr_rows
        variant = corr_rows[0]["variant"]
        for (cond, scheme), (pairs, x, y) in corr_vectors.items():
            for (a, b), xv, yv in zip(pairs, x, y):
                pair_tbl.append({"variant": variant, "condition": cond, "literature_weighting": scheme,
                                 "organ_1": a, "organ_2": b, "literature_weight": xv,
                                 "abs_partial_correlation": yv})
        draw_matrix_correlation_figure(corr_rows, corr_vectors, OUT_DIR / fig_name, drop_zero=drop_zero)
    corr_df = pd.DataFrame(all_corr_rows)
    corr_df.to_csv(stat_dir / "lit_vs_pc_matrix_correlation.csv", index=False, encoding="utf-8-sig")
    corr_df.to_excel(stat_dir / "lit_vs_pc_matrix_correlation.xlsx", index=False)
    pd.DataFrame(pair_tbl).to_csv(stat_dir / "lit_vs_pc_matrix_correlation_pairs.csv",
                                  index=False, encoding="utf-8-sig")
    print(corr_df[["variant", "condition", "literature_weighting", "n_pairs", "pearson_r", "pearson_p",
                   "spearman_rho", "spearman_p"]].to_string(index=False))

    print("[i] Drawing literature strength distribution figure …")
    draw_lit_strength_figure(
        {c: lit_vs_pc[c][0] for c in ("healthy", "obese")},
        OUT_DIR / "network_analysis_lit_strength_distribution.png",
    )
    for cond in ("healthy", "obese"):
        print(f"[i] Drawing {cond} literature vs PC table …")
        metrics, max_edges = lit_vs_pc[cond]
        draw_lit_vs_pc_table(cond, metrics, max_edges, pc_global_max,
                             OUT_DIR / f"network_analysis_{cond}_lit_vs_pc_table.png")


def full_condition_network_papercounts(condition: str, all_organs: list[str],
                                       restrict_to: list[str] | None = None) -> tuple[list, dict]:
    """
    Like full_condition_network() but edge weight = raw paper count (PMIDs) per pair.
    restrict_to: if given, only pairs between these organs are used.
    Returns (edges, paper_map) where paper_map[(o1,o2)] = paper_map[(o2,o1)] = count.
    """
    from run_metabolic_lit_search import (
        CONDITION_KEYWORDS, _compile_patterns, _filter_by_condition,
    )

    ref_dir      = HERE / "reference_network_only_metabolic"
    results_path = ref_dir / "metabolic_literature_results.json"
    data = json.loads(results_path.read_text(encoding="utf-8"))

    if condition != "all":
        pats = _compile_patterns(CONDITION_KEYWORDS[condition])

    paper_map: dict[tuple, int] = {}
    edges: list[tuple[str, str]] = []
    for key, entry in data.items():
        parts = key.split("|", 1)
        if len(parts) != 2:
            continue
        o1, o2 = min(parts), max(parts)
        if restrict_to is not None and (o1 not in restrict_to or o2 not in restrict_to):
            continue
        papers = entry.get("papers", [])
        if condition != "all":
            papers = _filter_by_condition(papers, pats)
        n_co = len([p for p in papers if p.get("pmid")])
        if n_co < 1:
            continue
        paper_map[(o1, o2)] = n_co
        paper_map[(o2, o1)] = n_co
        edges.append((o1, o2))

    return sorted(edges), paper_map


def draw_density_table(conditions_data: dict, out_path: Path) -> None:
    """
    Table comparing OO-weighted and paper-count-weighted (normalised) density
    across the three conditions, plus pairwise condition comparisons.
    """
    COND_LABELS = {"all": "No Condition", "healthy": "Healthy", "obese": "Obese"}
    DARK = "#1e293b"; MID = "#475569"; LIGHT = "#94a3b8"
    C_UNW   = "#64748b"
    C_OO    = {"all": "#64748b", "healthy": "#0284c7", "obese": "#c2410c"}
    C_PAPER = {"all": "#7c3aed", "healthy": "#059669", "obese": "#ea580c"}
    ROW_BG  = {"all": "#f8fafc", "healthy": "#f0f9ff", "obese": "#fff7ed"}
    PAIR_BG = ["#f1f5f9", "#e9f0f7", "#f5f0fb"]

    def _node_strengths(weighted_edges, nodes):
        s = {n: 0.0 for n in nodes}
        for e in weighted_edges:
            w = e[2] if len(e) == 3 else 1.0
            if e[0] in s: s[e[0]] += w
            if e[1] in s: s[e[1]] += w
        return s

    global_paper_max = max(
        (e[2] for cond in ("all", "healthy", "obese")
         for e in conditions_data[cond]["paper_edges"] if len(e) == 3),
        default=1.0,
    )

    rows = {}
    for cond in ("all", "healthy", "obese"):
        d     = conditions_data[cond]
        nodes = d["nodes"]
        N         = len(nodes)
        max_edges = N * (N - 1) / 2

        n_edges          = len(d["paper_edges"])

        oo_weight_sum    = sum(e[2] for e in d["oo_edges"] if len(e) == 3)
        density_oo_w     = oo_weight_sum / max_edges if max_edges > 0 else 0.0

        paper_weight_sum = sum(e[2] / global_paper_max for e in d["paper_edges"] if len(e) == 3)
        density_paper_w  = paper_weight_sum / max_edges if max_edges > 0 else 0.0

        rows[cond] = {
            "n_edges":   n_edges,
            "d_oo_w":    density_oo_w,
            "d_paper_w": density_paper_w,
        }

    # ── Per-edge weight vectors for paired t-tests ──────────────────────────
    ref_nodes = conditions_data["all"]["nodes"]
    all_pairs = [(a, b) for i, a in enumerate(ref_nodes) for b in ref_nodes[i + 1:]]

    def _oo_vec(cond):
        wmap = {(min(e[0], e[1]), max(e[0], e[1])): e[2]
                for e in conditions_data[cond]["oo_edges"] if len(e) == 3}
        return [wmap.get((min(a, b), max(a, b)), 0.0) for a, b in all_pairs]

    def _paper_vec(cond):
        wmap = {(min(e[0], e[1]), max(e[0], e[1])): e[2] / global_paper_max
                for e in conditions_data[cond]["paper_edges"] if len(e) == 3}
        return [wmap.get((min(a, b), max(a, b)), 0.0) for a, b in all_pairs]

    oo_vecs    = {c: _oo_vec(c)    for c in ("all", "healthy", "obese")}
    paper_vecs = {c: _paper_vec(c) for c in ("all", "healthy", "obese")}

    pairs = [
        ("all",     "healthy", "No Condition vs Healthy"),
        ("all",     "obese",   "No Condition vs Obese"),
        ("healthy", "obese",   "Healthy vs Obese"),
    ]
    pair_rows = []
    for c1, c2, label in pairs:
        t_oo,    p_oo    = scipy_stats.ttest_rel(oo_vecs[c1],    oo_vecs[c2])
        t_paper, p_paper = scipy_stats.ttest_rel(paper_vecs[c1], paper_vecs[c2])
        pair_rows.append({
            "label":      label,
            "dd_oo_w":    rows[c1]["d_oo_w"]    - rows[c2]["d_oo_w"],
            "t_oo":       float(t_oo),
            "p_oo":       float(p_oo),
            "dd_paper_w": rows[c1]["d_paper_w"] - rows[c2]["d_paper_w"],
            "t_paper":    float(t_paper),
            "p_paper":    float(p_paper),
        })

    # ── Figure ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(13, 5.5), dpi=DPI,
                             gridspec_kw={"height_ratios": [1, 0.85], "hspace": 0.55})
    fig.patch.set_facecolor("white")

    def _delta_clr(v):
        return "#059669" if v > 0 else ("#dc2626" if v < 0 else MID)

    def _sig_clr(p):
        return "#059669" if p < 0.05 else MID

    # ── Section 1: per-condition density ────────────────────────────────────
    ax1 = axes[0]; ax1.axis("off")
    ax1.set_title("Section 1 — Density per condition across two weighting schemes "
                  "(paired t-test: OO vs paper-count weighted, within each condition)",
                  fontsize=9, color=DARK, fontweight="bold", loc="left", pad=4)

    within = {}
    for c in ("all", "healthy", "obese"):
        t_w, p_w = scipy_stats.ttest_rel(oo_vecs[c], paper_vecs[c])
        within[c] = (float(t_w), float(p_w))

    hdrs1  = ["Condition", "Edges", "OO weighted (%)", "Paper weighted (%)",
              "t\n(OO vs Paper)", "p-value\n(OO vs Paper)", "sig"]
    col_x1 = [0.00, 0.17, 0.26, 0.41, 0.58, 0.71, 0.85]
    col_w1 = [0.17, 0.09, 0.15, 0.17, 0.13, 0.14, 0.10]

    y0 = 0.84
    for h, cx, cw in zip(hdrs1, col_x1, col_w1):
        ax1.text(cx + cw/2, y0, h, ha="center", va="center",
                 fontsize=8.5, fontweight="bold", color=DARK,
                 transform=ax1.transAxes, multialignment="center")
    ax1.plot([0.01, 0.99], [y0 - 0.13, y0 - 0.13], color=DARK, lw=0.8, transform=ax1.transAxes)

    for ri, cond in enumerate(("all", "healthy", "obese")):
        r  = rows[cond]
        y  = y0 - 0.28 - ri * 0.28
        ax1.add_patch(mpatches.FancyBboxPatch(
            (0.01, y - 0.12), 0.98, 0.24, boxstyle="square,pad=0",
            linewidth=0, facecolor=ROW_BG[cond], alpha=0.7, transform=ax1.transAxes))

        cells1  = [COND_LABELS[cond], str(r["n_edges"]),
                   f"{r['d_oo_w']*100:.1f}%",
                   f"{r['d_paper_w']*100:.1f}%",
                   f"{within[cond][0]:+.2f}",
                   f"{within[cond][1]:.3g}",
                   _sig_label(within[cond][1])]
        colors1 = [DARK, DARK, C_OO[cond], C_PAPER[cond], DARK, DARK,
                   _sig_clr(within[cond][1])]
        for cell, cx, cw, clr in zip(cells1, col_x1, col_w1, colors1):
            ax1.text(cx + cw/2, y, cell, ha="center", va="center",
                     fontsize=9, color=clr,
                     fontweight="bold" if clr not in (DARK, MID) else "normal",
                     transform=ax1.transAxes)

    # ── Section 2: pairwise condition differences + paired t-tests ──────────
    ax2 = axes[1]; ax2.axis("off")
    ax2.set_title(
        "Section 2 — Pairwise condition differences in density (paired t-test on per-edge weights)",
        fontsize=9, color=DARK, fontweight="bold", loc="left", pad=4)

    hdrs2  = ["Comparison",
              "Δ OO\n(%)",  "t\n(OO)",  "p\n(OO)",  "sig\n(OO)",
              "Δ Paper\n(%)", "t\n(Paper)", "p\n(Paper)", "sig\n(Paper)"]
    col_x2 = [0.00, 0.20, 0.30, 0.39, 0.49, 0.55, 0.65, 0.74, 0.84]
    col_w2 = [0.20, 0.10, 0.09, 0.10, 0.06, 0.10, 0.09, 0.10, 0.06]

    y0b = 0.84
    for h, cx, cw in zip(hdrs2, col_x2, col_w2):
        ax2.text(cx + cw/2, y0b, h, ha="center", va="center",
                 fontsize=8, fontweight="bold", color=DARK,
                 transform=ax2.transAxes, multialignment="center")
    ax2.plot([0.01, 0.99], [y0b - 0.15, y0b - 0.15], color=DARK, lw=0.8, transform=ax2.transAxes)

    for ri, pr in enumerate(pair_rows):
        y  = y0b - 0.30 - ri * 0.28
        ax2.add_patch(mpatches.FancyBboxPatch(
            (0.01, y - 0.12), 0.98, 0.24, boxstyle="square,pad=0",
            linewidth=0, facecolor=PAIR_BG[ri], alpha=0.6, transform=ax2.transAxes))

        cells2  = [pr["label"],
                   f"{pr['dd_oo_w']*100:+.1f}%",
                   f"{pr['t_oo']:+.2f}",
                   f"{pr['p_oo']:.3f}",
                   _sig_label(pr["p_oo"]),
                   f"{pr['dd_paper_w']*100:+.1f}%",
                   f"{pr['t_paper']:+.2f}",
                   f"{pr['p_paper']:.3f}",
                   _sig_label(pr["p_paper"])]
        colors2 = [DARK,
                   _delta_clr(pr["dd_oo_w"]),
                   DARK,
                   DARK,
                   _sig_clr(pr["p_oo"]),
                   _delta_clr(pr["dd_paper_w"]),
                   DARK,
                   DARK,
                   _sig_clr(pr["p_paper"])]
        fweights2 = ["normal", "bold", "normal", "normal", "bold",
                     "bold", "normal", "normal", "bold"]
        for cell, cx, cw, clr, fw in zip(cells2, col_x2, col_w2, colors2, fweights2):
            ax2.text(cx + cw/2, y, cell, ha="center", va="center",
                     fontsize=8.5, color=clr, fontweight=fw,
                     transform=ax2.transAxes)

    fig.text(0.99, 0.00,
             f"OO weighted: mean OO score across all possible pairs  ·  "
             f"Paper weighted: edge weight / global max ({global_paper_max:.0f} papers), "
             f"then averaged across all possible pairs  ·  "
             f"t-test: paired on per-edge weight vectors ({len(all_pairs)} pairs)  ·  "
             f"green = positive / significant, red = negative / ns",
             ha="right", va="bottom", fontsize=7, color=LIGHT)

    fig.suptitle(
        "Network Density Comparison: OO-Weighted vs Normalised Paper-Count-Weighted",
        fontsize=12, y=1.01, color=DARK, fontweight="bold",
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] {out_path.relative_to(HERE)}")


def draw_disorder_types_figure(conditions_data: dict, out_path: Path) -> None:
    """
    3×3 figure: rows = weighting scheme (unweighted / OO / paper-count),
                 cols = condition (No Condition / Healthy / Obese).
    Each panel: per-organ strength bars sorted by own strength descending,
    exponential fit overlay, β and R² in title.
    """
    COND_LABELS = {"all": "No Condition", "healthy": "Healthy", "obese": "Obese"}
    C_UNW   = {"all": "#64748b", "healthy": "#64748b", "obese": "#64748b"}
    C_OO    = {"all": "#64748b", "healthy": "#0284c7", "obese": "#c2410c"}
    C_PAPER = {"all": "#7c3aed", "healthy": "#059669", "obese": "#ea580c"}
    DARK = "#1e293b"; MID = "#475569"; LIGHT = "#94a3b8"
    SCHEME_LABELS = ["OO-coefficient", "Normalised paper count\n(÷ global max)"]

    cond_order   = ["all", "healthy", "obese"]
    scheme_order = ["oo", "norm"]

    def _node_strengths(weighted_edges, nodes):
        s = {n: 0.0 for n in nodes}
        for e in weighted_edges:
            w = e[2] if len(e) == 3 else 1.0
            if e[0] in s: s[e[0]] += w
            if e[1] in s: s[e[1]] += w
        return s

    global_paper_max = max(
        (e[2] for cond in cond_order for e in conditions_data[cond]["paper_edges"] if len(e) == 3),
        default=1.0,
    )

    # Pre-compute all panels
    panels = {}
    for cond in cond_order:
        d     = conditions_data[cond]
        nodes = d["nodes"]

        oo_str   = _node_strengths(d["oo_edges"], nodes)
        norm_str = _node_strengths(
            [(e[0], e[1], e[2] / global_paper_max) for e in d["paper_edges"] if len(e) == 3],
            nodes,
        )

        for scheme, s_dict in [("oo", oo_str), ("norm", norm_str)]:
            organs = sorted(nodes, key=lambda n: s_dict[n], reverse=True)
            vals   = [s_dict[n] for n in organs]
            beta, amp, off, r2, xf, yf = _fit_to_sequence(vals)
            panels[(cond, scheme)] = {
                "organs": organs, "vals": vals,
                "beta": beta, "r2": r2,
                "xf": xf, "yf": yf,
            }

    # ── Figure ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(
        2, 3, figsize=(22, 11), dpi=DPI,
        gridspec_kw={"hspace": 0.65, "wspace": 0.30},
    )
    fig.patch.set_facecolor("white")

    color_map = {"oo": C_OO, "norm": C_UNW}
    y_labels  = {
        "oo":   "OO-score strength",
        "norm": "Normalised paper-count strength\n(sum of count/global_max per organ)",
    }

    for row_i, scheme in enumerate(scheme_order):
        for col_i, cond in enumerate(cond_order):
            ax  = axes[row_i][col_i]
            p   = panels[(cond, scheme)]
            clr = color_map[scheme][cond]

            ax.set_facecolor("white")
            for sp in ("top", "right"):   ax.spines[sp].set_visible(False)
            for sp in ("left", "bottom"): ax.spines[sp].set_color(LIGHT)
            ax.grid(axis="y", color="#e2e8f0", linewidth=0.7, zorder=0)
            ax.set_axisbelow(True)

            x = np.arange(len(p["organs"]))
            ax.bar(x, p["vals"], color=clr, alpha=0.85, zorder=2, width=0.70)
            ax.plot(p["xf"], p["yf"], "--", color="#1e293b",
                    linewidth=1.8, alpha=0.75, zorder=3)

            labels = [n.replace(" ", "\n") for n in p["organs"]]
            ax.set_xticks(x)
            ax.set_xticklabels(labels, fontsize=6.0, rotation=45, ha="right",
                               color=DARK, multialignment="center")
            ax.tick_params(axis="y", labelsize=8, colors=DARK)
            ax.tick_params(axis="x", length=0)
            ax.set_xlim(-0.6, len(p["organs"]) - 0.4)

            if col_i == 0:
                ax.set_ylabel(y_labels[scheme], fontsize=9, color=MID)

            cond_prefix = f"{COND_LABELS[cond]}\n" if row_i == 0 else ""
            ax.set_title(
                cond_prefix +
                f"{SCHEME_LABELS[row_i]}  ·  β = {p['beta']*100:.2f}%  (R²={p['r2']:.3f})",
                fontsize=9.5, color=clr, pad=5, fontweight="bold",
            )

        axes[row_i][0].annotate(
            SCHEME_LABELS[row_i],
            xy=(0, 0.5), xycoords="axes fraction",
            xytext=(-0.28, 0.5), textcoords="axes fraction",
            fontsize=9, color=MID, ha="center", va="center",
            rotation=90, annotation_clip=False,
        )

    fig.suptitle(
        "Disorder Across Edge-Weighting Schemes and Conditions",
        fontsize=13, y=0.995, color=DARK, fontweight="bold",
    )

    fig.text(0.99, 0.005,
             "Each panel sorted by its own strength (descending); the fit is on this order  ·  "
             f"Normalised paper count: each edge ÷ global max ({global_paper_max:.0f} papers)  ·  "
             "Dashed line: a·exp(−β·x)+c fit",
             ha="right", va="bottom", fontsize=7, color=LIGHT)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] {out_path.relative_to(HERE)}")


def draw_disorder_types_table(conditions_data: dict, out_path: Path) -> None:
    """
    Table comparing disorder (β) computed from three edge-weighting schemes
    (unweighted / OO-weighted / paper-count-weighted) across the three conditions.
    β is computed as in draw_disorder_types_figure: _fit_to_sequence on the
    per-node strength sequence of each scheme, sorted by that scheme's own
    strength descending.
    """
    COND_LABELS = {"all": "No Condition", "healthy": "Healthy", "obese": "Obese"}
    DARK = "#1e293b"; MID = "#475569"; LIGHT = "#94a3b8"
    C_UNW   = "#64748b"
    C_OO    = {"all": "#64748b", "healthy": "#0284c7", "obese": "#c2410c"}
    C_PAPER = {"all": "#7c3aed", "healthy": "#059669", "obese": "#ea580c"}
    ROW_BG  = {"all": "#f8fafc", "healthy": "#f0f9ff", "obese": "#fff7ed"}
    PAIR_BG = ["#f1f5f9", "#e9f0f7", "#f5f0fb"]

    def _node_strengths(weighted_edges, nodes):
        s = {n: 0.0 for n in nodes}
        for e in weighted_edges:
            w = e[2] if len(e) == 3 else 1.0
            if e[0] in s: s[e[0]] += w
            if e[1] in s: s[e[1]] += w
        return s

    global_paper_max = max(
        (e[2] for cond in ("all", "healthy", "obese")
         for e in conditions_data[cond]["paper_edges"] if len(e) == 3),
        default=1.0,
    )

    rows = {}
    for cond in ("all", "healthy", "obese"):
        d     = conditions_data[cond]
        nodes = d["nodes"]

        oo_str   = _node_strengths(d["oo_edges"], nodes)
        norm_str = _node_strengths(
            [(e[0], e[1], e[2] / global_paper_max) for e in d["paper_edges"] if len(e) == 3],
            nodes,
        )

        organs_oo   = sorted(nodes, key=lambda n: oo_str[n],   reverse=True)
        organs_norm = sorted(nodes, key=lambda n: norm_str[n], reverse=True)

        oo_beta,   _, _, oo_r2,   _, _ = _fit_to_sequence([oo_str[n]   for n in organs_oo])
        norm_beta, _, _, norm_r2, _, _ = _fit_to_sequence([norm_str[n] for n in organs_norm])

        rows[cond] = {
            "n_edges":   len(d["paper_edges"]),
            "oo_beta":   oo_beta,   "oo_r2":   oo_r2,
            "norm_beta": norm_beta, "norm_r2": norm_r2,
        }

    pairs = [
        ("all",     "healthy", "No Condition vs Healthy"),
        ("all",     "obese",   "No Condition vs Obese"),
        ("healthy", "obese",   "Healthy vs Obese"),
    ]
    pair_rows = []
    for c1, c2, label in pairs:
        pair_rows.append({
            "label":   label,
            "dd_oo":   rows[c1]["oo_beta"]   - rows[c2]["oo_beta"],
            "dd_norm": rows[c1]["norm_beta"]  - rows[c2]["norm_beta"],
        })

    # ── Figure ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(11, 5.8), dpi=DPI,
                             gridspec_kw={"height_ratios": [1, 0.85], "hspace": 0.55})
    fig.patch.set_facecolor("white")

    def _delta_clr(v):
        return "#059669" if v > 0 else ("#dc2626" if v < 0 else MID)

    # ── Section 1: per-condition disorder ────────────────────────────────────
    ax1 = axes[0]; ax1.axis("off")
    ax1.set_title("Section 1 — Disorder (β) per condition across two edge-weighting schemes",
                  fontsize=9, color=DARK, fontweight="bold", loc="left", pad=4)

    hdrs1  = ["Condition", "Edges",
              "OO weighted β (%)", "OO weighted R²",
              "Norm. paper\nβ (%)", "Norm. paper\nR²"]
    col_x1 = [0.00, 0.17, 0.28, 0.49, 0.64, 0.83]
    col_w1 = [0.17, 0.11, 0.21, 0.15, 0.19, 0.17]

    y0 = 0.84
    for h, cx, cw in zip(hdrs1, col_x1, col_w1):
        ax1.text(cx + cw/2, y0, h, ha="center", va="center",
                 fontsize=8.5, fontweight="bold", color=DARK,
                 transform=ax1.transAxes, multialignment="center")
    ax1.plot([0.01, 0.99], [y0 - 0.13, y0 - 0.13], color=DARK, lw=0.8, transform=ax1.transAxes)

    for ri, cond in enumerate(("all", "healthy", "obese")):
        r  = rows[cond]
        y  = y0 - 0.28 - ri * 0.28
        ax1.add_patch(mpatches.FancyBboxPatch(
            (0.01, y - 0.12), 0.98, 0.24, boxstyle="square,pad=0",
            linewidth=0, facecolor=ROW_BG[cond], alpha=0.7, transform=ax1.transAxes))

        cells1  = [COND_LABELS[cond], str(r["n_edges"]),
                   f"{r['oo_beta']*100:.2f}%",   f"{r['oo_r2']:.3f}",
                   f"{r['norm_beta']*100:.2f}%", f"{r['norm_r2']:.3f}"]
        colors1 = [DARK, DARK,
                   C_OO[cond], DARK,
                   C_UNW, DARK]
        for cell, cx, cw, clr in zip(cells1, col_x1, col_w1, colors1):
            ax1.text(cx + cw/2, y, cell, ha="center", va="center",
                     fontsize=9, color=clr,
                     fontweight="bold" if clr not in (DARK, MID, C_UNW) else "normal",
                     transform=ax1.transAxes)

    # ── Section 2: pairwise condition β differences ───────────────────────────
    ax2 = axes[1]; ax2.axis("off")
    ax2.set_title("Section 2 — Pairwise condition differences in disorder β (first − second)",
                  fontsize=9, color=DARK, fontweight="bold", loc="left", pad=4)

    hdrs2  = ["Comparison", "Δ OO weighted β (%)", "Δ Norm. paper β (%)"]
    col_x2 = [0.00, 0.42, 0.70]
    col_w2 = [0.42, 0.28, 0.30]

    y0b = 0.84
    for h, cx, cw in zip(hdrs2, col_x2, col_w2):
        ax2.text(cx + cw/2, y0b, h, ha="center", va="center",
                 fontsize=8.5, fontweight="bold", color=DARK,
                 transform=ax2.transAxes, multialignment="center")
    ax2.plot([0.01, 0.99], [y0b - 0.13, y0b - 0.13], color=DARK, lw=0.8, transform=ax2.transAxes)

    for ri, pr in enumerate(pair_rows):
        y  = y0b - 0.28 - ri * 0.28
        ax2.add_patch(mpatches.FancyBboxPatch(
            (0.01, y - 0.12), 0.98, 0.24, boxstyle="square,pad=0",
            linewidth=0, facecolor=PAIR_BG[ri], alpha=0.6, transform=ax2.transAxes))

        cells2  = [pr["label"],
                   f"{pr['dd_oo']*100:+.2f}%",
                   f"{pr['dd_norm']*100:+.2f}%"]
        colors2 = [DARK,
                   _delta_clr(pr["dd_oo"]),
                   _delta_clr(pr["dd_norm"])]
        for cell, cx, cw, clr in zip(cells2, col_x2, col_w2, colors2):
            ax2.text(cx + cw/2, y, cell, ha="center", va="center",
                     fontsize=9, color=clr,
                     fontweight="bold" if clr not in (DARK, MID) else "normal",
                     transform=ax2.transAxes)

    fig.text(0.99, 0.00,
             f"OO weighted: sum of OO scores per organ  ·  "
             f"Norm. paper: each edge ÷ global max ({global_paper_max:.0f}), then summed per organ  ·  "
             "Each scheme sorted by its own strength (descending) before fitting  ·  "
             "green = first condition higher β, red = lower",
             ha="right", va="bottom", fontsize=7, color=LIGHT)

    fig.suptitle(
        "Network Disorder (β) Comparison: OO-Weighted vs Normalised Paper Count",
        fontsize=12, y=1.01, color=DARK, fontweight="bold",
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] {out_path.relative_to(HERE)}")


def draw_disorder_table(conditions_data: dict, out_path: Path) -> None:
    """
    Table with two sections:
    Section 1 — disorder β and R² per condition for normalised paper-count vs OO networks,
                 plus MWU test comparing the two weighting schemes.
    Section 2 — pairwise condition comparisons (MWU) within each weighting scheme.
    """
    COND_LABELS = {"all": "No Condition", "healthy": "Healthy", "obese": "Obese"}
    DARK = "#1e293b"; MID = "#475569"; LIGHT = "#94a3b8"
    C_PAPER = {"all": "#7c3aed", "healthy": "#059669", "obese": "#ea580c"}
    C_OO    = {"all": "#64748b", "healthy": "#0284c7", "obese": "#c2410c"}
    ROW_BG  = {"all": "#f8fafc", "healthy": "#f0f9ff", "obese": "#fff7ed"}
    PAIR_BG = ["#f1f5f9", "#e9f0f7", "#f5f0fb"]

    def _node_strengths(weighted_edges, nodes):
        s = {n: 0.0 for n in nodes}
        for e in weighted_edges:
            w = e[2] if len(e) == 3 else 1.0
            if e[0] in s: s[e[0]] += w
            if e[1] in s: s[e[1]] += w
        return s

    global_paper_max = max(
        (e[2] for cond in ("all", "healthy", "obese")
         for e in conditions_data[cond]["paper_edges"] if len(e) == 3),
        default=1.0,
    )

    # ── Compute per-condition metrics ────────────────────────────────────────
    rows = {}
    for cond in ("all", "healthy", "obese"):
        d     = conditions_data[cond]
        nodes = d["nodes"]

        oo_str    = _node_strengths(d["oo_edges"],    nodes)
        paper_str = _node_strengths(
            [(e[0], e[1], e[2] / global_paper_max) for e in d["paper_edges"] if len(e) == 3],
            nodes,
        )

        organs = sorted(nodes, key=lambda n: oo_str[n], reverse=True)
        oo_beta,    _, _, oo_r2,    _, _ = _fit_to_sequence([oo_str[n]    for n in organs])
        organs_paper = sorted(nodes, key=lambda n: paper_str[n], reverse=True)
        paper_beta, _, _, paper_r2, _, _ = _fit_to_sequence([paper_str[n] for n in organs_paper])

        paper_pos = sorted([paper_str[n] for n in nodes if paper_str[n] > 0], reverse=True)
        oo_pos    = sorted([oo_str[n]    for n in nodes if oo_str[n]    > 0], reverse=True)
        if paper_pos and oo_pos:
            u_raw_oo, p_raw_oo = scipy_stats.mannwhitneyu(paper_pos, oo_pos, alternative="two-sided")
        else:
            u_raw_oo, p_raw_oo = float("nan"), float("nan")

        rows[cond] = {
            "n_edges":    len(d["paper_edges"]),
            "oo_beta":    oo_beta,    "oo_r2":    oo_r2,
            "paper_beta": paper_beta, "paper_r2": paper_r2,
            "delta_beta": oo_beta - paper_beta,
            "u_raw_oo": float(u_raw_oo), "p_raw_oo": float(p_raw_oo),
            "paper_pos": paper_pos,
            "oo_pos":    oo_pos,
        }

    # ── Pairwise condition comparisons ───────────────────────────────────────
    pairs = [
        ("all", "healthy", "No Condition vs Healthy"),
        ("all", "obese",   "No Condition vs Obese"),
        ("healthy", "obese", "Healthy vs Obese"),
    ]
    pair_rows = []
    for c1, c2, label in pairs:
        p1 = rows[c1]["paper_pos"]; p2 = rows[c2]["paper_pos"]
        o1 = rows[c1]["oo_pos"];    o2 = rows[c2]["oo_pos"]
        if p1 and p2:
            u_p, pv_p = scipy_stats.mannwhitneyu(p1, p2, alternative="two-sided")
        else:
            u_p, pv_p = float("nan"), float("nan")
        if o1 and o2:
            u_o, pv_o = scipy_stats.mannwhitneyu(o1, o2, alternative="two-sided")
        else:
            u_o, pv_o = float("nan"), float("nan")

        db_paper = rows[c1]["paper_beta"] - rows[c2]["paper_beta"]
        db_oo    = rows[c1]["oo_beta"]    - rows[c2]["oo_beta"]
        pair_rows.append({
            "label": label,
            "db_paper": db_paper, "u_p": float(u_p), "p_p": float(pv_p),
            "db_oo":    db_oo,    "u_o": float(u_o), "p_o": float(pv_o),
        })

    # ── Figure ───────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(13, 5.8), dpi=DPI,
                             gridspec_kw={"height_ratios": [1, 1], "hspace": 0.55})
    fig.patch.set_facecolor("white")

    sig_col = {"***": "#dc2626", "**": "#ea580c", "*": "#ca8a04", "ns": MID, "—": LIGHT}

    def _p_str(p):
        if np.isnan(p): return "—"
        return f"{p:.3f}" if p >= 0.001 else "<0.001"

    def _u_str(u):
        return "—" if np.isnan(u) else str(int(u))

    # ── Section 1: per-condition β summary ──────────────────────────────────
    ax1 = axes[0]; ax1.axis("off")
    ax1.set_title("Section 1 — Disorder (β) per condition: Normalised paper count vs OO weighting",
                  fontsize=9, color=DARK, fontweight="bold", loc="left", pad=4)

    hdrs1  = ["Condition", "Edges", "Norm. paper\nβ (%)", "Norm. paper\nR²",
              "OO β (%)", "OO R²", "Δβ (OO − Norm.\npaper, %)", "U (MWU)", "p (MWU)", "Sig"]
    col_x1 = [0.00, 0.11, 0.20, 0.30, 0.40, 0.50, 0.60, 0.72, 0.82, 0.91]
    col_w1 = [0.11, 0.09, 0.10, 0.10, 0.10, 0.10, 0.12, 0.10, 0.09, 0.09]

    y0 = 0.85
    for h, cx, cw in zip(hdrs1, col_x1, col_w1):
        ax1.text(cx + cw/2, y0, h, ha="center", va="center",
                 fontsize=8.5, fontweight="bold", color=DARK,
                 transform=ax1.transAxes, multialignment="center")
    ax1.plot([0.01, 0.99], [y0 - 0.12, y0 - 0.12], color=DARK, lw=0.8, transform=ax1.transAxes)

    for ri, cond in enumerate(("all", "healthy", "obese")):
        r  = rows[cond]
        y  = y0 - 0.26 - ri * 0.26
        ax1.add_patch(mpatches.FancyBboxPatch(
            (0.01, y - 0.11), 0.98, 0.22, boxstyle="square,pad=0",
            linewidth=0, facecolor=ROW_BG[cond], alpha=0.7, transform=ax1.transAxes))

        cells1  = [COND_LABELS[cond], str(r["n_edges"]),
                   f"{r['paper_beta']*100:.2f}%", f"{r['paper_r2']:.3f}",
                   f"{r['oo_beta']*100:.2f}%",    f"{r['oo_r2']:.3f}",
                   f"{r['delta_beta']*100:+.2f}%",
                   _u_str(r["u_raw_oo"]), _p_str(r["p_raw_oo"]),
                   _sig_label(r["p_raw_oo"])]
        colors1 = [DARK, DARK,
                   C_PAPER[cond], DARK,
                   C_OO[cond],    DARK,
                   "#059669" if r["delta_beta"] > 0 else "#dc2626",
                   DARK, DARK, sig_col.get(_sig_label(r["p_raw_oo"]), DARK)]
        for cell, cx, cw, clr in zip(cells1, col_x1, col_w1, colors1):
            ax1.text(cx + cw/2, y, cell, ha="center", va="center",
                     fontsize=8.5, color=clr,
                     fontweight="bold" if clr not in (DARK, MID) else "normal",
                     transform=ax1.transAxes)

    # ── Section 2: pairwise condition comparisons ────────────────────────────
    ax2 = axes[1]; ax2.axis("off")
    ax2.set_title("Section 2 — Pairwise condition comparisons (MWU on sorted positive strength sequences)",
                  fontsize=9, color=DARK, fontweight="bold", loc="left", pad=4)

    hdrs2  = ["Comparison",
              "Δβ Norm.\npaper (%)", "U Norm.\npaper", "p Norm.\npaper", "Sig Norm.\npaper",
              "Δβ OO (%)",  "U OO",  "p OO",  "Sig OO"]
    col_x2 = [0.00, 0.22, 0.34, 0.44, 0.54, 0.63, 0.75, 0.85, 0.91]
    col_w2 = [0.22, 0.12, 0.10, 0.10, 0.09, 0.12, 0.10, 0.09, 0.09]

    y0b = 0.85
    for h, cx, cw in zip(hdrs2, col_x2, col_w2):
        ax2.text(cx + cw/2, y0b, h, ha="center", va="center",
                 fontsize=8.5, fontweight="bold", color=DARK,
                 transform=ax2.transAxes, multialignment="center")
    ax2.plot([0.01, 0.99], [y0b - 0.12, y0b - 0.12], color=DARK, lw=0.8, transform=ax2.transAxes)

    for ri, pr in enumerate(pair_rows):
        y  = y0b - 0.26 - ri * 0.26
        ax2.add_patch(mpatches.FancyBboxPatch(
            (0.01, y - 0.11), 0.98, 0.22, boxstyle="square,pad=0",
            linewidth=0, facecolor=PAIR_BG[ri], alpha=0.6, transform=ax2.transAxes))

        sig_p = _sig_label(pr["p_p"]); sig_o = _sig_label(pr["p_o"])
        cells2  = [pr["label"],
                   f"{pr['db_paper']*100:+.2f}%", _u_str(pr["u_p"]), _p_str(pr["p_p"]), sig_p,
                   f"{pr['db_oo']*100:+.2f}%",    _u_str(pr["u_o"]), _p_str(pr["p_o"]), sig_o]
        colors2 = [DARK,
                   "#059669" if pr["db_paper"] > 0 else "#dc2626", DARK, DARK,
                   sig_col.get(sig_p, DARK),
                   "#059669" if pr["db_oo"] > 0 else "#dc2626",    DARK, DARK,
                   sig_col.get(sig_o, DARK)]
        for cell, cx, cw, clr in zip(cells2, col_x2, col_w2, colors2):
            ax2.text(cx + cw/2, y, cell, ha="center", va="center",
                     fontsize=8.5, color=clr,
                     fontweight="bold" if clr not in (DARK, MID) else "normal",
                     transform=ax2.transAxes)

    fig.text(0.99, 0.00,
             "β from a·exp(−β·x)+c fit on each scheme's own strength-ordered sequence  ·  "
             f"paper count ÷ global max ({global_paper_max:.0f} papers)  ·  "
             "Δβ = first condition − second  ·  * p<0.05  ** p<0.01  *** p<0.001  ns not significant",
             ha="right", va="bottom", fontsize=7, color=LIGHT)

    fig.suptitle(
        "Disorder (β) Comparison: Normalised Paper-Count vs OO-Coefficient Networks",
        fontsize=12, y=1.01, color=DARK, fontweight="bold",
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] {out_path.relative_to(HERE)}")


def _fit_to_sequence(values: list) -> tuple:
    """
    Fit a*exp(-b*x)+c to an ordered sequence (0-indexed x).
    Returns (beta, amp, offset, r2, x_fit, y_fit).
    """
    def _exp(x, a, b, c):
        return a * np.exp(-b * x) + c

    x = np.arange(len(values), dtype=float)
    y = np.array(values, dtype=float)
    try:
        p0 = (float(max(y)), 1.0, float(min(y)))
        popt, _ = curve_fit(_exp, x, y, p0=p0, maxfev=5000)
    except (RuntimeError, ValueError):
        popt = (1.0, 0.0, float(np.mean(y)))

    x_fit = np.linspace(0, len(values) - 1, 200)
    y_fit = _exp(x_fit, *popt)

    y_pred = _exp(x, *popt)
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2     = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    return float(popt[1]), float(popt[0]), float(popt[2]), r2, x_fit, y_fit


def draw_papercount_vs_oo_figure(
    conditions_data: dict,
    out_path: Path,
) -> None:
    """
    Compare normalised paper-count networks vs OO-coefficient networks across conditions.
    Row 0: OO-score strength (actual values, organs sorted by OO strength desc).
    Row 1: Normalised paper-count strength (count/global_max, same organ order as row 0).
    Row 2: Stats table.
    """
    COND_LABELS = {"all": "No Condition", "healthy": "Healthy", "obese": "Obese"}
    C_PAPER = {"all": "#7c3aed", "healthy": "#059669", "obese": "#ea580c"}
    C_OO    = {"all": "#64748b", "healthy": "#0284c7", "obese": "#c2410c"}
    DARK = "#1e293b"; MID = "#475569"; LIGHT = "#94a3b8"

    cond_order = ["all", "healthy", "obese"]

    # ── Compute metrics ──────────────────────────────────────────────────────
    def _node_strengths(weighted_edges, nodes):
        s = {n: 0.0 for n in nodes}
        for e in weighted_edges:
            w = e[2] if len(e) == 3 else 1.0
            if e[0] in s: s[e[0]] += w
            if e[1] in s: s[e[1]] += w
        return s

    # Global max paper-count edge weight across all conditions (for comparable weighted density)
    global_paper_max = max(
        (e[2] for cond in cond_order for e in conditions_data[cond]["paper_edges"] if len(e) == 3),
        default=1.0,
    )

    results = {}
    for cond in cond_order:
        d     = conditions_data[cond]
        nodes = d["nodes"]

        oo_str   = _node_strengths(d["oo_edges"], nodes)
        norm_str = _node_strengths(
            [(e[0], e[1], e[2] / global_paper_max) for e in d["paper_edges"] if len(e) == 3],
            nodes,
        )

        # Sort organs by OO strength descending — shared order for both rows
        organs = sorted(nodes, key=lambda n: oo_str[n], reverse=True)

        oo_vals   = [oo_str[n]  for n in organs]
        norm_vals = [norm_str[n] for n in organs]

        # Fit exponential on the displayed (OO-order) sequence for each
        oo_beta,   oo_amp,   oo_off,   oo_r2,   oo_xf,   oo_yf   = _fit_to_sequence(oo_vals)
        norm_beta, norm_amp, norm_off, norm_r2, norm_xf, norm_yf = _fit_to_sequence(norm_vals)

        # Statistical tests on per-node strength sequences (all nodes)
        oo_seq   = [oo_str[n]  for n in nodes]
        norm_seq = [norm_str[n] for n in nodes]
        t_stat, p_t = scipy_stats.ttest_ind(norm_seq, oo_seq, equal_var=False)

        norm_pos = sorted([v for v in norm_seq if v > 0], reverse=True)
        oo_pos   = sorted([v for v in oo_seq   if v > 0], reverse=True)
        if norm_pos and oo_pos:
            u_stat, p_u = scipy_stats.mannwhitneyu(norm_pos, oo_pos, alternative="two-sided")
        else:
            u_stat, p_u = float("nan"), float("nan")

        N         = len(nodes)
        max_edges = N * (N - 1) / 2

        oo_weight_sum    = sum(e[2] for e in d["oo_edges"] if len(e) == 3)
        density_oo_w     = oo_weight_sum / max_edges if max_edges > 0 else 0.0

        paper_weight_sum = sum(e[2] / global_paper_max for e in d["paper_edges"] if len(e) == 3)
        density_paper_w  = paper_weight_sum / max_edges if max_edges > 0 else 0.0

        results[cond] = {
            "organs":          organs,
            "oo_vals":         oo_vals,   "norm_vals":  norm_vals,
            "oo_xf":           oo_xf,     "oo_yf":      oo_yf,
            "norm_xf":         norm_xf,   "norm_yf":    norm_yf,
            "oo_beta":         oo_beta,   "oo_r2":      oo_r2,
            "norm_beta":       norm_beta, "norm_r2":    norm_r2,
            "t_stat":          float(t_stat), "p_t":    float(p_t),
            "u_stat":          float(u_stat), "p_u":    float(p_u),
            "density_oo_w":    density_oo_w,
            "density_paper_w": density_paper_w,
        }

    # ── Figure ───────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(22, 18), dpi=DPI)
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(
        3, 3, figure=fig,
        height_ratios=[1.4, 1.4, 0.8],
        hspace=0.60, wspace=0.30,
    )

    def _style_ax(ax):
        ax.set_facecolor("white")
        for sp in ("top", "right"):  ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"): ax.spines[sp].set_color(LIGHT)
        ax.grid(axis="y", color="#e2e8f0", linewidth=0.7, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(axis="y", labelsize=8, colors=DARK)
        ax.tick_params(axis="x", length=0)

    for col_i, cond in enumerate(cond_order):
        r      = results[cond]
        organs = r["organs"]
        x      = np.arange(len(organs))
        labels = [n.replace(" ", "\n") for n in organs]

        # Row 0 — OO-score strength (actual values)
        ax0 = fig.add_subplot(gs[0, col_i])
        _style_ax(ax0)
        ax0.bar(x, r["oo_vals"], color=C_OO[cond], alpha=0.85, zorder=2, width=0.65)
        ax0.plot(r["oo_xf"], r["oo_yf"], "--", color="#1e293b", linewidth=1.8, alpha=0.75, zorder=3)
        ax0.set_xticks(x)
        ax0.set_xticklabels(labels, fontsize=6.5, rotation=45, ha="right",
                            color=DARK, multialignment="center")
        ax0.set_ylabel("OO-score strength", fontsize=9, color=MID)
        ax0.set_xlim(-0.6, len(organs) - 0.4)
        ax0.set_title(
            f"{COND_LABELS[cond]}  —  OO-coefficient strength\n"
            f"{len(conditions_data[cond]['oo_edges'])} edges  ·  "
            f"OO density {r['density_oo_w']*100:.1f}%  ·  "
            f"β = {r['oo_beta']*100:.2f}%  (R²={r['oo_r2']:.3f})",
            fontsize=9.5, color=C_OO[cond], pad=6, fontweight="bold",
        )

        # Row 1 — Normalised paper-count strength (same organ order as row 0)
        ax1 = fig.add_subplot(gs[1, col_i])
        _style_ax(ax1)
        ax1.bar(x, r["norm_vals"], color=C_PAPER[cond], alpha=0.85, zorder=2, width=0.65)
        ax1.plot(r["norm_xf"], r["norm_yf"], "--", color="#1e293b", linewidth=1.8, alpha=0.75, zorder=3)
        ax1.set_xticks(x)
        ax1.set_xticklabels(labels, fontsize=6.5, rotation=45, ha="right",
                            color=DARK, multialignment="center")
        ax1.set_ylabel("Norm. paper-count strength\n(sum of count/global_max per organ)", fontsize=9, color=MID)
        ax1.set_xlim(-0.6, len(organs) - 0.4)
        ax1.set_title(
            f"{COND_LABELS[cond]}  —  Normalised paper-count strength\n"
            f"β = {r['norm_beta']*100:.2f}%  (R²={r['norm_r2']:.3f})",
            fontsize=9.5, color=C_PAPER[cond], pad=6, fontweight="bold",
        )

    # ── Stats table ───────────────────────────────────────────────────────────
    ax_s = fig.add_subplot(gs[2, :])
    ax_s.axis("off")

    col_hdrs = [
        "Condition", "Edges",
        "Density\nOO weighted", "Density\npaper weighted",
        "OO β (%)", "Norm. paper β (%)",
        "U-stat", "p (MWU)", "sig (disorder)",
        "t-stat", "p (t-test)", "sig (strength)",
    ]
    col_w = [0.12, 0.06, 0.11, 0.12, 0.09, 0.11, 0.07, 0.08, 0.10, 0.07, 0.08, 0.10]

    y0 = 0.92
    x0 = 0.01
    for h, cw in zip(col_hdrs, col_w):
        ax_s.text(x0 + cw / 2, y0, h, ha="center", va="center",
                  fontsize=7.5, fontweight="bold", color=DARK,
                  transform=ax_s.transAxes)
        x0 += cw
    ax_s.plot([0.01, 0.99], [y0 - 0.07, y0 - 0.07],
              color=LIGHT, lw=0.8, transform=ax_s.transAxes)

    sig_colors = {"***": "#dc2626", "**": "#ea580c", "*": "#ca8a04",
                  "ns": "#64748b", "—": LIGHT}
    row_bgs = {"all": "#f8fafc", "healthy": "#f0f9ff", "obese": "#fff7ed"}

    for row_i, cond in enumerate(cond_order):
        r   = results[cond]
        y   = y0 - 0.20 - row_i * 0.24
        ax_s.add_patch(mpatches.FancyBboxPatch(
            (0.01, y - 0.10), 0.98, 0.20, boxstyle="square,pad=0",
            linewidth=0, facecolor=row_bgs[cond], alpha=0.6,
            transform=ax_s.transAxes))

        sig_u = _sig_label(r["p_u"])
        sig_t = _sig_label(r["p_t"])

        cells = [
            COND_LABELS[cond],
            str(len(conditions_data[cond]["paper_edges"])),
            f"{r['density_oo_w']*100:.1f}%",
            f"{r['density_paper_w']*100:.1f}%",
            f"{r['oo_beta']*100:.2f}%",
            f"{r['norm_beta']*100:.2f}%",
            f"{int(r['u_stat'])}" if not np.isnan(r["u_stat"]) else "—",
            f"{r['p_u']:.3f}"     if not np.isnan(r["p_u"]) and r["p_u"] >= 0.001 else (
                "<0.001" if not np.isnan(r["p_u"]) else "—"),
            sig_u,
            f"{r['t_stat']:.2f}",
            f"{r['p_t']:.3f}"     if r["p_t"] >= 0.001 else "<0.001",
            sig_t,
        ]
        x0 = 0.01
        for ci, (cell, cw) in enumerate(zip(cells, col_w)):
            clr = DARK; fw = "normal"
            if ci == 8:  clr = sig_colors.get(cell, DARK); fw = "bold" if cell not in ("ns","—") else "normal"
            if ci == 11: clr = sig_colors.get(cell, DARK); fw = "bold" if cell not in ("ns","—") else "normal"
            ax_s.text(x0 + cw / 2, y, cell, ha="center", va="center",
                      fontsize=7.5, color=clr, fontweight=fw,
                      transform=ax_s.transAxes)
            x0 += cw

    ax_s.text(0.99, -0.06,
              "Organs sorted by OO-score strength (descending) — same order in both rows  ·  "
              "β fitted on displayed organ order  ·  "
              f"Norm. paper weighted density: each edge / global max ({global_paper_max:.0f} papers)  ·  "
              "OO weighted density: mean OO score across all possible pairs  ·  "
              "t-test on per-node strength sequences  ·  MWU on sorted positive strengths  ·  "
              "* p<0.05  ** p<0.01  *** p<0.001  ns not significant",
              ha="right", va="bottom", fontsize=7, color=LIGHT,
              transform=ax_s.transAxes)

    fig.suptitle(
        "Normalised Paper-Count Networks vs OO-Coefficient Networks — Strength & Disorder Comparison",
        fontsize=13, y=0.995, color=DARK, fontweight="bold",
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] {out_path.relative_to(HERE)}")


if __name__ == "__main__":
    main()
