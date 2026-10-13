"""
Key player analyses for the literature-based reference network (all-papers scope:
all 91 organ pairs, papers filtered by condition keywords, no OO threshold --
the same scope as visualizations_all_papers/).

1. Condition specificity
   Paper-level enrichment of every key player, healthy vs obese
   (papers matching BOTH condition keyword sets are excluded from this
   comparison), and each condition vs the rest of the network. Fisher's exact
   test, Benjamini-Hochberg FDR within each category, log2 odds ratio with a
   0.5 continuity correction.
   -> statistics/key_players/key_player_condition_enrichment.xlsx / .csv
   -> visualizations_all_papers/key_player_analysis/key_player_enrichment_healthy_vs_obese.png

2. Mediators per organ and per axis
   Organ x player and axis x player matrices (unique papers per organ / pair),
   top-3 players per axis.
   -> statistics/key_players/key_player_organ_axis_<condition>.xlsx
   -> visualizations_all_papers/key_player_analysis/organ_<share|lift>_<condition>.png
   -> visualizations_all_papers/key_player_analysis/axis_share_<condition>.png

Counting unit: unique papers (a paper counts once per player, however many
organ pairs it matches), unlike the pie plots, which sum per-pair counts.

Run: uv run python generate_key_player_analysis.py
"""
from __future__ import annotations

import json
import math
import time
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from Literature_Search.pubmed_search import extract_key_players
from run_metabolic_lit_search import CONDITION_KEYWORDS, _compile_patterns, _filter_by_condition

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "reference_network_only_metabolic" / "metabolic_literature_results.json"
STAT_DIR = HERE / "statistics" / "key_players"
FIG_DIR = HERE / "visualizations_all_papers" / "key_player_analysis"
AXIS_DIR = FIG_DIR / "axis_share"          # axis heatmaps (with and without insulin/glucose)

CATS = ("hormones", "metabolites", "proteins")
CAT_LABEL = {"hormones": "Hormones", "metabolites": "Metabolites", "proteins": "Proteins / Transporters"}
COND_LABEL = {"all": "No condition (all papers)", "healthy": "Healthy", "obese": "Obese"}
C_HEALTHY, C_OBESE = "#2563a8", "#c2570c"
MIN_PAPERS = 5          # min papers (healthy + obese) for a term to be tested
TOP_BARS = 15
TOP_COLS = 12
TOP_AXES = 20
MIN_MENTIONS_LIFT = 30  # min mentions in a category for an organ's enrichment to be shown
DPI = 300

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


# ── data ─────────────────────────────────────────────────────────────────────

def load():
    data = json.loads(RESULTS.read_text(encoding="utf-8"))
    papers, pair_pmids = {}, {}
    for key, e in data.items():
        o1, o2 = e.get("organ1"), e.get("organ2")
        if not o1 or not o2:
            continue
        pair = (min(o1, o2), max(o1, o2))
        ids = set()
        for p in e.get("papers", []):
            pm = p.get("pmid")
            if pm:
                papers[pm] = p
                ids.add(pm)
        pair_pmids[pair] = ids
    return papers, pair_pmids


def paper_players(papers: dict) -> dict[str, dict[str, set]]:
    out = {}
    for i, (pm, p) in enumerate(papers.items()):
        kp = extract_key_players([p])
        out[pm] = {c: set(kp.get(f"{c}_counts", {})) for c in CATS}
        if (i + 1) % 5000 == 0:
            print(f"   players extracted for {i + 1}/{len(papers)} papers", flush=True)
    return out


# ── statistics ───────────────────────────────────────────────────────────────

def _logfact(n: int) -> np.ndarray:
    return np.concatenate([[0.0], np.cumsum(np.log(np.arange(1, n + 1)))])


def fisher_two_sided(a: int, n1: int, c: int, n2: int, lf: np.ndarray) -> float:
    """2x2 [[a, n1-a], [c, n2-c]]; two-sided exact p (sum of tables no likelier than observed)."""
    k, n = a + c, n1 + n2
    lo, hi = max(0, k - n2), min(k, n1)
    xs = np.arange(lo, hi + 1)
    lp = (lf[n1] - lf[xs] - lf[n1 - xs]) + (lf[n2] - lf[k - xs] - lf[n2 - k + xs]) - (lf[n] - lf[k] - lf[n - k])
    obs = lp[a - lo]
    return float(min(1.0, np.exp(lp[lp <= obs + 1e-9]).sum()))


def bh(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    run = 1.0
    for rank, idx in zip(range(n, 0, -1), order[::-1]):
        run = min(run, p[idx] * n / rank)
        q[idx] = run
    return q


def enrichment(pp: dict, set_a: set, set_b: set, label_a: str, label_b: str, min_papers=MIN_PAPERS):
    """Terms enriched in set_a vs set_b, per category. Positive log2OR = more in set_a."""
    na, nb = len(set_a), len(set_b)
    lf = _logfact(na + nb)
    rows = []
    for cat in CATS:
        ca, cb = Counter(), Counter()
        for pm in set_a:
            ca.update(pp[pm][cat])
        for pm in set_b:
            cb.update(pp[pm][cat])
        terms = [t for t in set(ca) | set(cb) if ca[t] + cb[t] >= min_papers]
        res = []
        for t in terms:
            a, c = ca[t], cb[t]
            p = fisher_two_sided(a, na, c, nb, lf)
            lor = math.log2(((a + .5) / (na - a + .5)) / ((c + .5) / (nb - c + .5)))
            res.append((t, a, c, p, lor))
        if not res:
            continue
        q = bh(np.array([r[3] for r in res]))
        for (t, a, c, p, lor), qq in zip(res, q):
            rows.append({"category": cat, "term": t,
                         f"n_{label_a}": a, f"N_{label_a}": na, f"pct_{label_a}": round(100 * a / na, 2),
                         f"n_{label_b}": c, f"N_{label_b}": nb, f"pct_{label_b}": round(100 * c / nb, 2),
                         "log2_odds_ratio": round(lor, 3), "p_value": p, "q_value_BH": qq,
                         "significant_q<0.05": bool(qq < 0.05),
                         "enriched_in": label_a if lor > 0 else label_b})
    return pd.DataFrame(rows)


# ── figure 1: diverging enrichment bars ──────────────────────────────────────

def plot_enrichment(df: pd.DataFrame, out: Path, n_h: int, n_o: int):
    fig, axes = plt.subplots(1, 3, figsize=(27, 11))
    for ax, cat in zip(axes, CATS):
        d = df[(df.category == cat) & (df["significant_q<0.05"])]
        top_h = d[d.log2_odds_ratio > 0].nlargest(TOP_BARS, "log2_odds_ratio")
        top_o = d[d.log2_odds_ratio < 0].nsmallest(TOP_BARS, "log2_odds_ratio")
        sel = pd.concat([top_h.sort_values("log2_odds_ratio"), top_o.sort_values("log2_odds_ratio", ascending=False)])
        sel = pd.concat([top_o.sort_values("log2_odds_ratio", ascending=False),
                         top_h.sort_values("log2_odds_ratio")]).reset_index(drop=True)
        ax.set_title(f"{CAT_LABEL[cat]}\n{len(top_h)} healthy-enriched · {len(top_o)} obese-enriched "
                     f"(of {int((df.category == cat).sum())} tested)", fontsize=12, fontweight="bold")
        if sel.empty:
            ax.text(0.5, 0.5, "no significant terms (q < 0.05)", ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            continue
        y = np.arange(len(sel))
        cols = [C_HEALTHY if v > 0 else C_OBESE for v in sel.log2_odds_ratio]
        ax.barh(y, sel.log2_odds_ratio, color=cols, height=0.75)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r.term}   ({r.pct_healthy:.1f}% vs {r.pct_obese:.1f}%)"
                            for r in sel.itertuples()], fontsize=9)
        ax.axvline(0, color="#334155", lw=0.8)
        lim = max(abs(sel.log2_odds_ratio).max() * 1.1, 1)
        ax.set_xlim(-lim, lim)
        ax.set_xlabel("log2 odds ratio  (← obese-enriched | healthy-enriched →)")
    fig.suptitle(f"Condition-specific key players: healthy vs obese literature "
                 f"(papers exclusive to each condition: healthy n={n_h:,}, obese n={n_o:,})\n"
                 f"Fisher's exact test, Benjamini–Hochberg q < 0.05; labels = % of papers mentioning the term (healthy vs obese)",
                 fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  [ok]", out.relative_to(HERE))


# ── figure 2: organ / axis heatmaps ──────────────────────────────────────────

def count_terms(pmids, pp, cat) -> Counter:
    c = Counter()
    for pm in pmids:
        c.update(pp[pm][cat])
    return c


def build_matrices(cond: str, pair_pmids: dict, members: set | None, pp: dict):
    """Return organ and axis count matrices for one condition."""
    pairs = {p: (ids if members is None else ids & members) for p, ids in pair_pmids.items()}
    pairs = {p: ids for p, ids in pairs.items() if ids}
    organ_ids = defaultdict(set)
    for (o1, o2), ids in pairs.items():
        organ_ids[o1] |= ids
        organ_ids[o2] |= ids
    net_ids = set().union(*pairs.values())
    out = {"pairs": pairs, "organ_ids": dict(organ_ids), "net_ids": net_ids}
    out["net"] = {c: count_terms(net_ids, pp, c) for c in CATS}
    out["organ"] = {o: {c: count_terms(ids, pp, c) for c in CATS} for o, ids in organ_ids.items()}
    out["axis"] = {p: {c: count_terms(ids, pp, c) for c in CATS} for p, ids in pairs.items()}
    return out


def _heat(ax, M, rows, cols, cmap, vmin, vmax, title, fmt, cbar_label):
    im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=55, ha="right", fontsize=9)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=9)
    ax.set_title(title, fontsize=12, fontweight="bold")
    mid = (vmin + vmax) / 2
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if np.isnan(v):
                continue
            ax.text(j, i, fmt(v), ha="center", va="center", fontsize=8,
                    color="white" if abs(v - mid) > (vmax - vmin) * 0.28 else "#1e293b")
    cb = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label(cbar_label, fontsize=8)
    for s in ax.spines.values():
        s.set_visible(False)


def heat_data(mats: dict, ref: dict, exclude: dict | None = None) -> dict:
    """Heatmap matrices for one condition, laid out on the rows and columns of the
    reference condition (obese): its organ order, its top axes and its top players
    per category, so that the figures of all conditions can be compared directly.
    Rows/columns that do not occur in this condition stay blank (n = 0).

    exclude: {category: {terms}} to leave out (e.g. insulin, glucose). Excluded
    terms are dropped from the columns (the next-ranked player moves in) and from
    the denominator of the axis shares, so the remaining players add up to 100%."""
    exclude = exclude or {}
    net = mats["net"]
    organs = sorted(ref["organ_ids"], key=lambda o: -len(ref["organ_ids"][o]))
    axes_sorted = sorted(ref["pairs"], key=lambda p: -len(ref["pairs"][p]))[:TOP_AXES]
    cols = {cat: [t for t, _ in ref["net"][cat].most_common(TOP_COLS + len(exclude.get(cat, ())))
                  if t not in exclude.get(cat, ())][:TOP_COLS] for cat in CATS}
    out = {
        "organs": organs, "axes": axes_sorted, "cols": cols,
        "org_labels": [f"{o} ({len(mats['organ_ids'].get(o, ())):,})" for o in organs],
        "ax_labels": [f"{a}–{b} ({len(mats['pairs'].get((a, b), ())):,})" for a, b in axes_sorted],
        "organ_share": {}, "organ_lift": {}, "axis_share": {},
    }
    for cat in CATS:
        tot_net = sum(net[cat].values())
        Ms = np.full((len(organs), TOP_COLS), np.nan)
        Ml = np.full((len(organs), TOP_COLS), np.nan)
        for i, o in enumerate(organs):
            c = mats["organ"].get(o, {}).get(cat, Counter())
            tot = sum(c.values())
            for j, t in enumerate(cols[cat]):
                share = 100 * c.get(t, 0) / tot if tot else np.nan
                Ms[i, j] = share
                base = 100 * net[cat].get(t, 0) / tot_net if tot_net else 0
                # too few mentions or player absent from the network -> no meaningful ratio (blank)
                Ml[i, j] = np.log2(max(share, 0.1) / base) if tot >= MIN_MENTIONS_LIFT and base else np.nan
        Ma = np.full((len(axes_sorted), TOP_COLS), np.nan)
        for i, p in enumerate(axes_sorted):
            c = mats["axis"].get(p, {}).get(cat, Counter())
            tot = sum(v for t, v in c.items() if t not in exclude.get(cat, ()))
            for j, t in enumerate(cols[cat]):
                Ma[i, j] = 100 * c.get(t, 0) / tot if tot else np.nan
        out["organ_share"][cat], out["organ_lift"][cat], out["axis_share"][cat] = Ms, Ml, Ma
    return out


def _fmt_share(v):
    return f"{v:.0f}" if v >= 10 else f"{v:.1f}"


def draw_heatmaps(cond: str, data: dict, vmax: dict, kinds=("organ_share", "organ_lift", "axis_share"),
                  suffix: str = "", variant_note: str = "", panel_suffix: dict | None = None):
    """vmax[(kind, cat)]: common colour limit across conditions for the share heatmaps.
    The axis-share figures go to their own subfolder (AXIS_DIR); suffix/variant_note/panel_suffix
    label a variant (e.g. without insulin and glucose)."""
    panel_suffix = panel_suffix or {}
    ref_note = f"columns: top {TOP_COLS} players of the obese network (identical in all conditions)"
    for kind in kinds:
        is_axis = kind == "axis_share"
        rows = data["ax_labels"] if is_axis else data["org_labels"]
        fig, axes = plt.subplots(1, 3, figsize=(24, 9 if is_axis else 7.6))
        for ax, cat in zip(axes, CATS):
            M, cols = data[kind][cat], data["cols"][cat]
            if kind == "organ_lift":
                _heat(ax, M, rows, cols, "RdBu_r", -3, 3, CAT_LABEL[cat] + panel_suffix.get(cat, ""),
                      lambda v: f"{v:+.1f}", "log2 (organ share / network share)")
            else:
                _heat(ax, M, rows, cols, "YlGnBu", 0, vmax[(kind, cat)], CAT_LABEL[cat] + panel_suffix.get(cat, ""),
                      _fmt_share, "% of axis' mentions" if is_axis else "% of organ's mentions")
        if kind == "organ_share":
            title = (f"Share of each organ's key player mentions — {COND_LABEL[cond]}\n"
                     f"Rows: organs in the order of the obese figure (n unique papers in this condition); "
                     f"{ref_note}; cell text = % share")
        elif kind == "organ_lift":
            title = (f"Enrichment of each player in an organ relative to the whole network — {COND_LABEL[cond]}\n"
                     f"Rows: organs in the order of the obese figure (n unique papers); {ref_note}; "
                     f"cell text and colour = log2 enrichment (blank: < {MIN_MENTIONS_LIFT} mentions in the category "
                     f"or player absent; absent players floored at 0.1% share)")
        else:
            title = (f"Share of key player mentions per organ–organ axis — {COND_LABEL[cond]}\n"
                     f"Rows: the {TOP_AXES} axes with most papers in the obese network (n unique papers in this "
                     f"condition; blank row = no papers); {ref_note}")
            if variant_note:
                title += f"\n{variant_note}"
        fig.suptitle(title, fontsize=12.5, y=1.0)
        fig.tight_layout(rect=[0, 0, 1, 0.93 if not is_axis else 0.94])
        out_dir = AXIS_DIR if is_axis else FIG_DIR
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / f"{kind}{suffix}_{cond}.png"
        for attempt in range(5):        # a PNG that was just viewed can be locked for a moment on Windows
            try:
                fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
                break
            except OSError:
                if attempt == 4:
                    raise
                time.sleep(2)
        plt.close(fig)
        print("  [ok]", out.relative_to(HERE))


def write_tables(cond: str, mats: dict):
    net = mats["net"]
    sheets = {}
    # organ x player (counts and shares), per category
    for cat in CATS:
        terms = [t for t, _ in net[cat].most_common(40)]
        organs = sorted(mats["organ_ids"], key=lambda o: -len(mats["organ_ids"][o]))
        cnt = pd.DataFrame({t: [mats["organ"][o][cat].get(t, 0) for o in organs] for t in terms}, index=organs)
        tot = pd.Series({o: sum(mats["organ"][o][cat].values()) for o in organs})
        cnt.insert(0, "n_papers_organ", [len(mats["organ_ids"][o]) for o in organs])
        cnt.insert(1, "total_mentions", tot.values)
        cnt.index.name = "organ"
        sheets[f"organ_{cat}_counts"] = cnt.reset_index()
        shr = (cnt[terms].div(tot.values, axis=0) * 100).round(2)
        shr.insert(0, "n_papers_organ", cnt["n_papers_organ"])
        shr.index.name = "organ"
        sheets[f"organ_{cat}_share_pct"] = shr.reset_index()
    # top-3 per axis
    rows = []
    for p in sorted(mats["pairs"], key=lambda p: -len(mats["pairs"][p])):
        r = {"organ_1": p[0], "organ_2": p[1], "n_papers": len(mats["pairs"][p])}
        for cat in CATS:
            top = mats["axis"][p][cat].most_common(3)
            r[f"top3_{cat}"] = "; ".join(f"{t} ({n})" for t, n in top) if top else ""
            r[f"n_distinct_{cat}"] = len(mats["axis"][p][cat])
        rows.append(r)
    sheets["axis_top3_players"] = pd.DataFrame(rows)
    with pd.ExcelWriter(STAT_DIR / f"key_player_organ_axis_{cond}.xlsx") as xw:
        for name, df in sheets.items():
            df.to_excel(xw, sheet_name=name[:31], index=False)
            ws = xw.sheets[name[:31]]
            ws.freeze_panes = "B2"
            for i, col in enumerate(df.columns, start=1):
                w = max(len(str(col)), *(len(str(v)) for v in df[col].head(100))) + 2
                ws.column_dimensions[ws.cell(1, i).column_letter].width = min(w, 70)
    sheets["axis_top3_players"].to_csv(STAT_DIR / f"key_player_axis_top3_{cond}.csv", index=False, encoding="utf-8-sig")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    STAT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    print("[i] loading papers")
    papers, pair_pmids = load()
    print(f"    {len(papers):,} unique papers, {len(pair_pmids)} pairs")
    print("[i] extracting key players per paper")
    pp = paper_players(papers)

    members = {}
    for cond in ("healthy", "obese"):
        pats = _compile_patterns(CONDITION_KEYWORDS[cond])
        members[cond] = {pm for pm, p in papers.items() if _filter_by_condition([p], pats)}
    allp = set(papers)
    both = members["healthy"] & members["obese"]
    h_only, o_only = members["healthy"] - both, members["obese"] - both
    print(f"    healthy {len(members['healthy']):,} | obese {len(members['obese']):,} | both {len(both):,} "
          f"| healthy-only {len(h_only):,} | obese-only {len(o_only):,}")

    # 1. condition specificity
    print("[i] condition enrichment")
    hv = enrichment(pp, h_only, o_only, "healthy", "obese")
    hv.insert(0, "comparison", "healthy_only vs obese_only")
    h_rest = enrichment(pp, members["healthy"], allp - members["healthy"], "healthy", "rest")
    h_rest.insert(0, "comparison", "healthy vs rest of network")
    o_rest = enrichment(pp, members["obese"], allp - members["obese"], "obese", "rest")
    o_rest.insert(0, "comparison", "obese vs rest of network")
    summary = []
    for name, df in (("healthy_vs_obese", hv), ("healthy_vs_rest", h_rest), ("obese_vs_rest", o_rest)):
        for cat in CATS:
            d = df[df.category == cat]
            summary.append({"comparison": name, "category": cat, "terms_tested": len(d),
                            "significant_q<0.05": int(d["significant_q<0.05"].sum())})
    with pd.ExcelWriter(STAT_DIR / "key_player_condition_enrichment.xlsx") as xw:
        pd.DataFrame(summary).to_excel(xw, sheet_name="summary", index=False)
        for name, df in (("healthy_vs_obese", hv), ("healthy_vs_rest", h_rest), ("obese_vs_rest", o_rest)):
            df.sort_values(["category", "q_value_BH"]).to_excel(xw, sheet_name=name, index=False)
            xw.sheets[name].freeze_panes = "A2"
    hv.sort_values(["category", "q_value_BH"]).to_csv(STAT_DIR / "key_player_enrichment_healthy_vs_obese.csv",
                                                       index=False, encoding="utf-8-sig")
    print(pd.DataFrame(summary).to_string(index=False))
    for cat in CATS:
        d = hv[(hv.category == cat) & hv["significant_q<0.05"]]
        for side in ("healthy", "obese"):
            s = d[d.enriched_in == side].sort_values("q_value_BH").head(6)
            print(f"  {cat} {side}-enriched:", [(r.term, round(r.log2_odds_ratio, 2)) for r in s.itertuples()])
    plot_enrichment(hv, FIG_DIR / "key_player_enrichment_healthy_vs_obese.png", len(h_only), len(o_only))

    # 2. organ / axis: build all conditions first, then draw them on the layout of the obese condition
    all_mats = {}
    for cond, mem in (("all", None), ("healthy", members["healthy"]), ("obese", members["obese"])):
        print(f"[i] organ/axis matrices: {cond}")
        all_mats[cond] = build_matrices(cond, pair_pmids, mem, pp)
        write_tables(cond, all_mats[cond])
    ref = all_mats["obese"]
    heat = {cond: heat_data(m, ref) for cond, m in all_mats.items()}
    vmax = {}
    for kind in ("organ_share", "axis_share"):
        for cat in CATS:
            vmax[(kind, cat)] = float(np.nanmax([np.nanmax(heat[c][kind][cat]) for c in heat]))
    for cond in ("all", "healthy", "obese"):
        draw_heatmaps(cond, heat[cond], vmax)

    # axis heatmaps without the two players that dominate almost every axis
    # (insulin among hormones, glucose among metabolites): the next-ranked obese players
    # move in as columns and shares are renormalised without the removed players
    excl = {"hormones": {"Insulin"}, "metabolites": {"glucose"}}
    heat_x = {cond: heat_data(m, ref, exclude=excl) for cond, m in all_mats.items()}
    vmax_x = {("axis_share", cat): float(np.nanmax([np.nanmax(heat_x[c]["axis_share"][cat]) for c in heat_x]))
              for cat in CATS}
    for cond in ("all", "healthy", "obese"):
        draw_heatmaps(cond, heat_x[cond], vmax_x, kinds=("axis_share",), suffix="_no_insulin_glucose",
                      variant_note="Without insulin (hormones) and glucose (metabolites); shares renormalised to the "
                                   "remaining players, next-ranked obese players added as columns",
                      panel_suffix={"hormones": " (without insulin)", "metabolites": " (without glucose)"})


if __name__ == "__main__":
    main()
