"""
Compute partial correlation networks for healthy and obese cohorts from FDG-PET data.

Patients are split at BMI 25 (<=25 = healthy, >25 = obese).
Missing organ values are imputed with within-group column means before fitting.

Outputs:
  *_cohort_connections.csv        — binary adjacency matrix (p < 0.05, uncorrected)
                                    replaces the literature-based binary reference CSVs
  *_cohort_PC_connections.csv     — full partial correlation matrix (r-values)
  *_cohort_PC_significance.csv    — binary significance mask (p < 0.05, uncorrected)
  *_cohort_PC_pvalues.csv         — raw p-values

Canonical organ names match all_organ_connections.csv / the HTML visualisation.
"""
import numpy as np
import pandas as pd
from scipy.stats import t as t_dist
from pathlib import Path

HERE = Path(__file__).parent
DATA_FILE = HERE / "reference_network_only_metabolic" / "Data FDG healthy 20250910 anonym.csv"
OUT_DIR   = HERE / "reference_network_only_metabolic"

ORGAN_COLS = [
    "AdrGl max cB",
    "Brain cB",
    "intestine cB",
    "heart cB",
    "Kidney cB",
    "Liver cB",
    "Lungs max cB",
    "Pancreas cB",
    "Gluteus cB",
    "Spleen cB",
    "Fat sub cB",
    "BM max cB",
    "Thyroid cB",
]

# Canonical names matching all_organ_connections.csv and the HTML visualisation
ORGAN_NAMES = [
    "Adrenal Glands",
    "Brain",
    "Small Intestine",
    "Heart",
    "Kidney",
    "Liver",
    "Lung",
    "Pancreas",
    "Muscle",
    "Spleen",
    "WAT",
    "Bone Marrow",
    "Thyroid",
]

# Full 14-organ list used in the HTML (Colon is absent from FDG data → always 0)
HTML_ORGANS = [
    "Adrenal Glands", "Bone Marrow", "Brain", "Colon", "Heart",
    "Kidney", "Liver", "Lung", "Muscle", "Pancreas",
    "Small Intestine", "Spleen", "Thyroid", "WAT",
]


def partial_correlation_network(
    data: pd.DataFrame,
    calc_significance: str | float = "significance",
    alpha: float = 0.05,
) -> dict:
    cov = np.cov(data.values, rowvar=False)
    precision = np.linalg.pinv(cov)
    d = np.sqrt(np.diag(precision))
    partial_corr = -precision / np.outer(d, d)
    np.fill_diagonal(partial_corr, 1.0)
    partial_corr_df = pd.DataFrame(partial_corr, index=data.columns, columns=data.columns)

    results = {"partial_corr": partial_corr_df}

    n, p = data.shape
    dof = n - p

    if isinstance(calc_significance, str):
        # Clamp to avoid divide-by-zero at |r|=1 (diagonal handled by fill_diagonal)
        r2_clamp = np.clip(partial_corr ** 2, 0, 1 - 1e-10)
        t_stats  = partial_corr * np.sqrt(dof / (1 - r2_clamp))
        p_values = 2 * (1 - t_dist.cdf(np.abs(t_stats), df=dof))
        np.fill_diagonal(p_values, 1.0)
        pval_df = pd.DataFrame(p_values, index=data.columns, columns=data.columns)
        results["p_values"] = pval_df

        if calc_significance == "significance":
            sig_mask = (p_values < alpha).astype(int)
        elif calc_significance == "bonferroni":
            m = p * (p - 1) / 2
            sig_mask = (p_values < alpha / m).astype(int)
        else:
            raise ValueError("calc_significance must be 'significance' or 'bonferroni'")

        np.fill_diagonal(sig_mask, 0)
        results["edges"] = pd.DataFrame(sig_mask, index=data.columns, columns=data.columns)

    elif isinstance(calc_significance, (float, int)):
        edge_mask = np.abs(partial_corr) >= calc_significance
        np.fill_diagonal(edge_mask, 0)
        results["edges"] = pd.DataFrame(edge_mask.astype(int), index=data.columns, columns=data.columns)

    return results


def _to_binary_adjacency_csv(sig_df: pd.DataFrame) -> str:
    """Render significance DataFrame as a 0/1 full adjacency CSV string."""
    organs = list(sig_df.columns)
    lines  = ["," + ",".join(organs)]
    for row in organs:
        cells = [row]
        for col in organs:
            if row == col:
                cells.append("0")
            else:
                cells.append(str(int(sig_df.loc[row, col])))
        lines.append(",".join(cells))
    return "\n".join(lines)


def _to_html_adjacency_csv(sig_df: pd.DataFrame, pc_df: pd.DataFrame) -> str:
    """
    Expand a 13-organ PC significance matrix to the 14-organ HTML format.
    Significant edges carry the actual r-value; non-significant and Colon are 0.
    """
    pc_organs = list(sig_df.columns)
    lines = ["," + ",".join(HTML_ORGANS)]
    for row in HTML_ORGANS:
        cells = [row]
        for col in HTML_ORGANS:
            if row == col:
                cells.append("0")
            elif row in pc_organs and col in pc_organs and sig_df.loc[row, col] == 1:
                cells.append(str(float(pc_df.loc[row, col])))
            else:
                cells.append("0")
        lines.append(",".join(cells))
    return "\n".join(lines)


def main():
    df = pd.read_csv(DATA_FILE)
    print(f"[i] Loaded {len(df)} rows from {DATA_FILE.name}")

    df = df[df["BMI"].notna()].copy()
    print(f"[i] {len(df)} rows have BMI values")

    col_map = dict(zip(ORGAN_COLS, ORGAN_NAMES))

    healthy_raw = df[df["BMI"] <= 25][ORGAN_COLS].rename(columns=col_map).copy()
    obese_raw   = df[df["BMI"] >  25][ORGAN_COLS].rename(columns=col_map).copy()
    print(f"[i] Healthy (BMI ≤ 25): {len(healthy_raw)} patients")
    print(f"[i] Obese   (BMI > 25): {len(obese_raw)} patients")

    healthy = healthy_raw.fillna(healthy_raw.mean())
    obese   = obese_raw.fillna(obese_raw.mean())

    html_csvs = {}

    for group_name, data in [("healthy", healthy), ("obese", obese)]:
        print(f"\n[i] Computing partial correlation network — {group_name} (n={len(data)})")
        res = partial_correlation_network(data, calc_significance="significance", alpha=0.05)

        pc  = res["partial_corr"]
        sig = res["edges"]
        pv  = res["p_values"]
        n_edges = int(sig.values[np.triu_indices_from(sig.values, k=1)].sum())
        print(f"    Significant edges (p < 0.05): {n_edges}")

        # Weighted r-value matrix
        pc_path = OUT_DIR / f"{group_name}_cohort_PC_connections.csv"
        pc.to_csv(pc_path)
        print(f"[ok] Partial correlation matrix → {pc_path.relative_to(HERE)}")

        # Binary significance mask (p < 0.05)
        sig_path = OUT_DIR / f"{group_name}_cohort_PC_significance.csv"
        sig.to_csv(sig_path)
        print(f"[ok] Significance mask          → {sig_path.relative_to(HERE)}")

        # p-values
        pval_path = OUT_DIR / f"{group_name}_cohort_PC_pvalues.csv"
        pv.to_csv(pval_path)
        print(f"[ok] p-values                   → {pval_path.relative_to(HERE)}")

        # Binary adjacency CSV replacing the literature reference
        binary_csv = _to_binary_adjacency_csv(sig)
        conn_path  = OUT_DIR / f"{group_name}_cohort_connections.csv"
        conn_path.write_text(binary_csv, encoding="utf-8")
        print(f"[ok] Binary adjacency (reference) → {conn_path.relative_to(HERE)}")

        html_csvs[group_name] = _to_html_adjacency_csv(sig, pc)

        # Print significant edges
        organs = list(sig.columns)
        print(f"    Significant edges (p < 0.05):")
        for i, a in enumerate(organs):
            for j, b in enumerate(organs):
                if j > i and sig.loc[a, b] == 1:
                    print(f"      {a} — {b}: r={pc.loc[a,b]:.4f}, p={pv.loc[a,b]:.4f}")

    # Patch HTML with new PC-based cohort references
    html_path = OUT_DIR / "metabolic_literature_network.html"
    if html_path.exists():
        _patch_html(html_path, html_csvs["healthy"], html_csvs["obese"])
        print(f"\n[ok] Patched HTML cohort references → {html_path.relative_to(HERE)}")
    else:
        print(f"\n[!] HTML not found at {html_path} — skipping patch")


def _patch_html(html_path: Path, h_csv: str, o_csv: str) -> None:
    """Replace _HC and _OC template-literal content in the generated HTML."""
    text = html_path.read_text(encoding="utf-8")

    import re
    h_escaped = h_csv.replace("`", "\\`")
    o_escaped = o_csv.replace("`", "\\`")

    # Replace  var _HC = `...`;
    text = re.sub(
        r"var _HC = `[^`]*`",
        f"var _HC = `{h_escaped}`",
        text,
        flags=re.DOTALL,
    )
    # Replace  var _OC = `...`;
    text = re.sub(
        r"var _OC = `[^`]*`",
        f"var _OC = `{o_escaped}`",
        text,
        flags=re.DOTALL,
    )

    html_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
