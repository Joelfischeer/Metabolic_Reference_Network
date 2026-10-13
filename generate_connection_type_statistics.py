"""
Connection-type statistics for the metabolic reference network.

Builds detailed statistic tables for every scope that is also plotted by
generate_static_visualizations.py (same pair selection, same classification
file), so the tables match the PNGs one-to-one:

  all            all 91-pair edge set, all papers
  healthy        OO-elbow-threshold pairs, healthy-filtered papers
  obese          OO-elbow-threshold pairs, obese-filtered papers
  healthy_all    pairs with >=1 healthy paper (no threshold)
  obese_all      pairs with >=1 obese paper (no threshold)

Output: statistics/connection_types/
  connection_type_stats_<scope>.xlsx   (one workbook per scope, several sheets)
  connection_type_stats_<scope>_*.csv  (same tables as CSV)
  connection_type_comparison.xlsx/.csv (type counts side by side, all scopes)

Run: uv run python generate_connection_type_statistics.py
"""
from __future__ import annotations

import json
import math
from collections import Counter
from itertools import combinations
from pathlib import Path

import pandas as pd

import generate_static_visualizations as gsv
from generate_static_visualizations import (
    CONNECTION_TYPE_LABELS, UNCLASSIFIED_KEY, UNCLASSIFIED_LABEL,
    _condition_pairs_all, _load_connection_types, _load_edge_filter,
    _oo_pairs_for_condition,
)
from run_metabolic_lit_search import (
    CONDITION_KEYWORDS, _compile_patterns, _filter_by_condition,
)

HERE = Path(__file__).resolve().parent
REF_DIR = HERE / "reference_network_only_metabolic"
OUT = HERE / "statistics" / "connection_types"

# scope -> (label, condition or None, pair selection: "edges" | "threshold" | "any")
SCOPES = {
    "all":         ("All papers (no condition)",       None,      "edges"),
    "healthy":     ("Healthy (OO threshold)",          "healthy", "threshold"),
    "obese":       ("Obese (OO threshold)",            "obese",   "threshold"),
    "healthy_all": ("Healthy (all papers, no thr.)",   "healthy", "any"),
    "obese_all":   ("Obese (all papers, no thr.)",     "obese",   "any"),
    # PET partial-correlation cohort edges, classified with the condition's own
    # LLM file (run_metabolic_lit_search.py --condition healthy|obese)
    "healthy_cohort": ("Healthy (PET cohort edges)",   "healthy", "cohort"),
    "obese_cohort":   ("Obese (PET cohort edges)",     "obese",   "cohort"),
}


def _load_pair_papers(condition: str | None) -> dict[tuple, list[dict]]:
    data = json.loads((REF_DIR / "metabolic_literature_results.json").read_text(encoding="utf-8"))
    pats = _compile_patterns(CONDITION_KEYWORDS[condition]) if condition else None
    out = {}
    for key, entry in data.items():
        parts = key.split("|", 1)
        if len(parts) != 2:
            continue
        papers = entry.get("papers", [])
        if pats is not None:
            papers = _filter_by_condition(papers, pats)
        out[(min(parts), max(parts))] = papers
    return out


def _oo_values(pair_papers: dict) -> dict[tuple, float]:
    """Overlap coefficient-style OO: n_co / sqrt(n_organ1 * n_organ2)."""
    organ_pmids: dict[str, set] = {}
    for (o1, o2), papers in pair_papers.items():
        pm = {p.get("pmid") for p in papers if p.get("pmid")}
        organ_pmids.setdefault(o1, set()).update(pm)
        organ_pmids.setdefault(o2, set()).update(pm)
    oo = {}
    for (o1, o2), papers in pair_papers.items():
        n_co = len([p for p in papers if p.get("pmid")])
        d = math.sqrt(len(organ_pmids.get(o1, ())) * len(organ_pmids.get(o2, ())))
        oo[(o1, o2)] = n_co / d if d > 0 else 0.0
    return oo


def _shannon(counts: list[int]) -> float:
    tot = sum(counts)
    if tot == 0:
        return 0.0
    return abs(-sum(c / tot * math.log2(c / tot) for c in counts if c))


def _label(k: str) -> str:
    return CONNECTION_TYPE_LABELS.get(k, k)


def build_scope(scope: str, label: str, condition: str | None, mode: str):
    if mode == "cohort":
        pair_types = _load_connection_types(f"metabolic_connection_types_{condition}.json")
        pair_papers = _load_pair_papers(condition)
    else:
        pair_types = _load_connection_types()
        pair_papers = _load_pair_papers(condition)
    oo = _oo_values(pair_papers)

    if mode == "edges":
        pairs = _load_edge_filter(REF_DIR / "all_organ_connections.csv")
    elif mode == "cohort":
        pairs = _load_edge_filter(REF_DIR / f"{condition}_cohort_connections.csv")
    elif mode == "threshold":
        pairs = _oo_pairs_for_condition(condition)
    else:
        pairs = _condition_pairs_all(condition)
    pairs = sorted(set(pairs))
    n_pairs = len(pairs)

    # ── pair-level table ─────────────────────────────────────────────────────
    rows = []
    for p in pairs:
        types = pair_types.get(p, [])
        papers = pair_papers.get(p, [])
        pmids = {x.get("pmid") for x in papers if x.get("pmid")}
        years = [int(x["year"]) for x in papers if str(x.get("year", "")).isdigit()]
        rows.append({
            "organ_1": p[0], "organ_2": p[1],
            "n_types": len(types),
            "types": "; ".join(_label(t) for t in types) if types else UNCLASSIFIED_LABEL,
            "n_papers": len(pmids),
            "OO": round(oo.get(p, 0.0), 4),
            "median_year": float(pd.Series(years).median()) if years else None,
            "latest_year": max(years) if years else None,
        })
    pair_df = pd.DataFrame(rows).sort_values(["n_papers", "OO"], ascending=False).reset_index(drop=True)

    # ── type summary ─────────────────────────────────────────────────────────
    keys = list(CONNECTION_TYPE_LABELS)
    n_classified = sum(1 for p in pairs if pair_types.get(p))
    total_assign = sum(len(pair_types.get(p, [])) for p in pairs)
    t_rows = []
    for k in keys:
        if k == UNCLASSIFIED_KEY:
            tp = [p for p in pairs if not pair_types.get(p)]
        else:
            tp = [p for p in pairs if k in pair_types.get(p, [])]
        if not tp:
            t_rows.append({"type": _label(k), "n_pairs": 0})
            continue
        counts = [len({x.get("pmid") for x in pair_papers.get(p, []) if x.get("pmid")}) for p in tp]
        pmids = set()
        for p in tp:
            pmids |= {x.get("pmid") for x in pair_papers.get(p, []) if x.get("pmid")}
        organs = Counter(o for p in tp for o in p)
        solo = sum(1 for p in tp if len(pair_types.get(p, [])) == 1)
        top = sorted(tp, key=lambda p: -len({x.get("pmid") for x in pair_papers.get(p, [])}))[:3]
        s = pd.Series(counts, dtype=float)
        t_rows.append({
            "type": _label(k),
            "n_pairs": len(tp),
            "pct_of_pairs": round(100 * len(tp) / n_pairs, 1) if n_pairs else 0,
            "pct_of_assignments": round(100 * len(tp) / total_assign, 1)
                                  if total_assign and k != UNCLASSIFIED_KEY else None,
            "n_pairs_exclusive": solo,
            "n_pairs_shared_with_other_type": len(tp) - solo if k != UNCLASSIFIED_KEY else None,
            "n_unique_papers": len(pmids),
            "papers_per_pair_mean": round(s.mean(), 2),
            "papers_per_pair_median": float(s.median()),
            "papers_per_pair_min": int(s.min()),
            "papers_per_pair_max": int(s.max()),
            "papers_per_pair_sd": round(s.std(), 2) if len(s) > 1 else None,
            "OO_mean": round(sum(oo.get(p, 0) for p in tp) / len(tp), 4),
            "n_organs_involved": len(organs),
            "most_frequent_organs": ", ".join(f"{o} ({c})" for o, c in organs.most_common(3)),
            "top_pairs_by_papers": "; ".join(f"{a}-{b}" for a, b in top),
        })
    type_df = pd.DataFrame(t_rows).sort_values("n_pairs", ascending=False).reset_index(drop=True)

    # ── organ x type matrix ──────────────────────────────────────────────────
    organs = sorted({o for p in pairs for o in p})
    cols = [k for k in keys if k != UNCLASSIFIED_KEY] + [UNCLASSIFIED_KEY]
    m = pd.DataFrame(0, index=organs, columns=[_label(k) for k in cols])
    for p in pairs:
        for t in (pair_types.get(p) or [UNCLASSIFIED_KEY]):
            for o in p:
                m.loc[o, _label(t)] += 1
    m = m.loc[:, (m != 0).any(axis=0)]
    n_conn = Counter(o for p in pairs for o in p)
    m.insert(0, "n_connections", [n_conn[o] for o in m.index])
    m["n_type_assignments"] = m.drop(columns="n_connections").sum(axis=1)
    tcols = [c for c in m.columns if c not in ("n_connections", "n_type_assignments")]
    m["n_distinct_types"] = (m[tcols] > 0).sum(axis=1)
    m["shannon_diversity_bits"] = [round(_shannon(list(m.loc[o, tcols])), 3) for o in m.index]
    m["dominant_type"] = m[tcols].idxmax(axis=1)
    m["dominant_type_share_pct"] = (100 * m[tcols].max(axis=1) / m["n_type_assignments"]).round(1)
    m = m.sort_values("n_connections", ascending=False)
    m.index.name = "organ"
    m = m.reset_index()

    # ── type co-occurrence matrix (pairs carrying both types) ────────────────
    used = [k for k in cols if k != UNCLASSIFIED_KEY and any(k in pair_types.get(p, []) for p in pairs)]
    co = pd.DataFrame(0, index=[_label(k) for k in used], columns=[_label(k) for k in used])
    for p in pairs:
        ts = [t for t in pair_types.get(p, []) if t in used]
        for t in ts:
            co.loc[_label(t), _label(t)] += 1
        for a, b in combinations(ts, 2):
            co.loc[_label(a), _label(b)] += 1
            co.loc[_label(b), _label(a)] += 1
    co.index.name = "type (diagonal = pairs with type)"
    co = co.reset_index()

    # ── overview ─────────────────────────────────────────────────────────────
    ntypes = Counter(len(pair_types.get(p, [])) for p in pairs)
    all_pm = set()
    for p in pairs:
        all_pm |= {x.get("pmid") for x in pair_papers.get(p, []) if x.get("pmid")}
    pc = pair_df["n_papers"]
    ov = [
        ("Scope", label),
        ("Condition filter", condition or "none"),
        ("Pair selection", {"edges": "edge set (all_organ_connections.csv)",
                            "threshold": "OO kneedle elbow threshold",
                            "any": ">=1 condition paper",
                            "cohort": "PET cohort edges (own LLM classification)"}[mode]),
        ("Organ pairs in scope", n_pairs),
        ("Organs in scope", len(organs)),
        ("Pairs with >=1 type", n_classified),
        ("Pairs without type", n_pairs - n_classified),
        ("% pairs classified", round(100 * n_classified / n_pairs, 1) if n_pairs else 0),
        ("Total type assignments", total_assign),
        ("Mean types per pair (all pairs)", round(total_assign / n_pairs, 2) if n_pairs else 0),
        ("Mean types per classified pair", round(total_assign / n_classified, 2) if n_classified else 0),
        ("Pairs with 0 / 1 / 2 / 3 types", " / ".join(str(ntypes.get(i, 0)) for i in range(4))),
        ("Distinct types used", len(used)),
        ("Unique papers in scope", len(all_pm)),
        ("Papers per pair: mean", round(pc.mean(), 2)),
        ("Papers per pair: median", float(pc.median())),
        ("Papers per pair: max", int(pc.max())),
        ("Pairs with 0 papers", int((pc == 0).sum())),
        ("Shannon diversity of types (bits)",
         round(_shannon([sum(1 for p in pairs if t in pair_types.get(p, [])) for t in used]), 3)),
    ]
    ov_df = pd.DataFrame(ov, columns=["metric", "value"])
    return {"overview": ov_df, "type_summary": type_df, "organ_x_type": m,
            "type_cooccurrence": co, "pair_details": pair_df}, type_df


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    comp = {}
    for scope, (label, cond, mode) in SCOPES.items():
        print(f"[i] {scope}")
        try:
            tables, type_df = build_scope(scope, label, cond, mode)
        except FileNotFoundError as exc:
            print(f"  [!] skipped: {exc}")
            continue
        with pd.ExcelWriter(OUT / f"connection_type_stats_{scope}.xlsx") as xw:
            for name, df in tables.items():
                df.to_excel(xw, sheet_name=name[:31], index=False)
                ws = xw.sheets[name[:31]]
                ws.freeze_panes = "A2"
                for i, col in enumerate(df.columns, start=1):
                    width = max(len(str(col)), *(len(str(v)) for v in df[col].head(200))) + 2
                    ws.column_dimensions[ws.cell(1, i).column_letter].width = min(width, 60)
                df.to_csv(OUT / f"connection_type_stats_{scope}_{name}.csv", index=False, encoding="utf-8-sig")
        comp[label] = type_df.set_index("type")[["n_pairs", "pct_of_pairs", "n_unique_papers"]]
        print(tables["overview"].to_string(index=False, header=False))
        print()
        print(type_df[["type", "n_pairs", "pct_of_pairs", "n_pairs_exclusive", "n_unique_papers",
                       "papers_per_pair_median"]].to_string(index=False))
        print()

    if comp:
        wide = pd.concat(comp, axis=1)
        wide = wide.fillna(0).sort_values((list(comp)[0], "n_pairs"), ascending=False)
        wide.columns = [f"{a} | {b}" for a, b in wide.columns]
        wide.index.name = "type"
        wide.to_csv(OUT / "connection_type_comparison.csv", encoding="utf-8-sig")
        wide.to_excel(OUT / "connection_type_comparison.xlsx")
        print("[ok] comparison:", OUT / "connection_type_comparison.xlsx")


if __name__ == "__main__":
    main()
