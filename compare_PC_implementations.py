"""
Compare three partial correlation implementations on the FDG-PET cohort data.

1. Current (covariance pinv)  — precision matrix via pseudo-inverse of covariance
2. Correlation pinv           — precision matrix via pseudo-inverse of correlation matrix
3. Regression-based           — classical OLS residual approach (Lawley 1940)

For each pair (i,j), the regression approach computes:
  e_i = residuals of regressing X_i on all other variables except X_j
  e_j = residuals of regressing X_j on all other variables except X_i
  partial_corr(i,j) = pearson(e_i, e_j)
This is the closed-form definition of partial correlation.
"""

import numpy as np
import pandas as pd
from scipy.stats import t as t_dist, pearsonr
from pathlib import Path

HERE = Path(__file__).parent
DATA_FILE = HERE / "reference_network_only_metabolic" / "Data FDG healthy 20250910 anonym.csv"

ORGAN_COLS = [
    "AdrGl max cB", "Brain cB", "intestine cB", "heart cB", "Kidney cB",
    "Liver cB", "Lungs max cB", "Pancreas cB", "Gluteus cB", "Spleen cB",
    "Fat sub cB", "BM max cB", "Thyroid cB",
]
ORGAN_NAMES = [
    "Adrenal Glands", "Brain", "Small intestine", "Heart", "Kidney",
    "Liver", "Lungs", "Pancreas", "Muscle", "Spleen",
    "WAT", "Bone Marrow", "Thyroid",
]


# ── Implementation 1: current (covariance pinv) ──────────────────────────────
def pc_covariance_pinv(X: np.ndarray) -> np.ndarray:
    cov = np.cov(X, rowvar=False)
    precision = np.linalg.pinv(cov)
    d = np.sqrt(np.diag(precision))
    pc = -precision / np.outer(d, d)
    np.fill_diagonal(pc, 1.0)
    return pc


# ── Implementation 2: correlation pinv ───────────────────────────────────────
def pc_correlation_pinv(X: np.ndarray) -> np.ndarray:
    # Standardise to unit variance before inverting → precision = R^{-1}
    corr = np.corrcoef(X, rowvar=False)
    precision = np.linalg.pinv(corr)
    d = np.sqrt(np.diag(precision))
    pc = -precision / np.outer(d, d)
    np.fill_diagonal(pc, 1.0)
    return pc


# ── Implementation 3: regression-based (OLS residuals) ───────────────────────
def pc_regression(X: np.ndarray) -> np.ndarray:
    n, p = X.shape
    pc = np.eye(p)
    for i in range(p):
        for j in range(i + 1, p):
            # Conditioning set: all variables except i and j
            others = [k for k in range(p) if k != i and k != j]
            Z = X[:, others]  # (n, p-2)

            # Add intercept
            Z1 = np.column_stack([np.ones(n), Z])

            # Residuals of X_i ~ Z (OLS)
            beta_i = np.linalg.lstsq(Z1, X[:, i], rcond=None)[0]
            e_i = X[:, i] - Z1 @ beta_i

            # Residuals of X_j ~ Z (OLS)
            beta_j = np.linalg.lstsq(Z1, X[:, j], rcond=None)[0]
            e_j = X[:, j] - Z1 @ beta_j

            r, _ = pearsonr(e_i, e_j)
            pc[i, j] = r
            pc[j, i] = r

    return pc


# ── p-values (Bonferroni) ─────────────────────────────────────────────────────
def bonferroni_edges(pc: np.ndarray, n: int, alpha: float = 0.05) -> np.ndarray:
    p = pc.shape[0]
    dof = n - p
    # Clamp to avoid sqrt of negative (numerical noise at |r|~1)
    r2 = np.clip(pc ** 2, 0, 1 - 1e-10)
    t_stats = pc * np.sqrt(dof / (1 - r2))
    pvals = 2 * (1 - t_dist.cdf(np.abs(t_stats), df=dof))
    np.fill_diagonal(pvals, 1.0)
    m = p * (p - 1) / 2
    sig = (pvals < alpha / m).astype(int)
    np.fill_diagonal(sig, 0)
    return sig, pvals


# ── Comparison helpers ────────────────────────────────────────────────────────
def upper_tri(mat: np.ndarray) -> np.ndarray:
    idx = np.triu_indices_from(mat, k=1)
    return mat[idx]


def compare(name_a: str, mat_a: np.ndarray,
            name_b: str, mat_b: np.ndarray) -> None:
    diff = mat_a - mat_b
    ut   = upper_tri(diff)
    print(f"  {name_a} vs {name_b}:")
    print(f"    max |diff| = {np.max(np.abs(ut)):.6f}")
    print(f"    mean |diff| = {np.mean(np.abs(ut)):.6f}")
    print(f"    RMSE = {np.sqrt(np.mean(ut**2)):.6f}")


def sig_edges(sig: np.ndarray, names: list[str]) -> list[str]:
    out = []
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if j > i and sig[i, j]:
                out.append(f"{a} — {b}")
    return out


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    df = pd.read_csv(DATA_FILE)
    df = df[df["BMI"].notna()].copy()
    col_map = dict(zip(ORGAN_COLS, ORGAN_NAMES))

    METHODS = [
        ("Covariance pinv (current)", pc_covariance_pinv),
        ("Correlation pinv",          pc_correlation_pinv),
        ("Regression-based (OLS)",    pc_regression),
    ]

    for group_label, bmi_cond in [("HEALTHY (BMI ≤ 25)", df["BMI"] <= 25),
                                   ("OBESE   (BMI > 25)", df["BMI"] >  25)]:
        raw = df[bmi_cond][ORGAN_COLS].rename(columns=col_map)
        data = raw.fillna(raw.mean())
        X = data.values
        n, p = X.shape
        names = list(data.columns)

        print(f"\n{'='*60}")
        print(f"  {group_label}  (n={n}, p={p})")
        print(f"{'='*60}")

        results = {}
        for mname, mfunc in METHODS:
            pc = mfunc(X)
            sig, pvals = bonferroni_edges(pc, n)
            edges = sig_edges(sig, names)
            results[mname] = {"pc": pc, "sig": sig, "edges": edges}
            print(f"\n  [{mname}]")
            print(f"    Significant edges (Bonferroni): {len(edges)}")
            for e in edges:
                i, j = [(k, names.index(part)) for k, part in enumerate(names)
                         if part == e.split(" — ")[0] or part == e.split(" — ")[1]]
                # simpler:
                a, b = e.split(" — ")
                ia, ib = names.index(a), names.index(b)
                print(f"      {e}: r={pc[ia, ib]:.4f}, p={pvals[ia, ib]:.5f}")

        print(f"\n  -- Element-wise differences (upper triangle, {p*(p-1)//2} pairs) --")
        mnames = [m[0] for m in METHODS]
        mats   = [results[m]["pc"] for m in mnames]
        compare(mnames[0], mats[0], mnames[1], mats[1])
        compare(mnames[0], mats[0], mnames[2], mats[2])
        compare(mnames[1], mats[1], mnames[2], mats[2])

        print(f"\n  -- Edge agreement (Bonferroni) --")
        sets = {m: set(results[m]["edges"]) for m in mnames}
        for i, ma in enumerate(mnames):
            for mb in mnames[i+1:]:
                union = sets[ma] | sets[mb]
                inter = sets[ma] & sets[mb]
                only_a = sets[ma] - sets[mb]
                only_b = sets[mb] - sets[ma]
                print(f"  {ma} vs {mb}:")
                print(f"    shared={len(inter)}  only-first={len(only_a)}  only-second={len(only_b)}")
                if only_a: print(f"      only in first:  {', '.join(only_a)}")
                if only_b: print(f"      only in second: {', '.join(only_b)}")


if __name__ == "__main__":
    main()
