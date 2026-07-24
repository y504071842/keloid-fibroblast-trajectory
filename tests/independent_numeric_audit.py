from __future__ import annotations

import json
import math
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "results"
OUT = Path(os.getenv("KFT_NUMERIC_AUDIT_OUTPUT", AUDIT / "Independent_Numeric_Audit"))
OUT.mkdir(parents=True, exist_ok=True)


def cliff_delta(x: np.ndarray, y: np.ndarray) -> float:
    # x: keloid/recurrence, y: control/non-recurrence
    gt = sum(float(a > b) for a in x for b in y)
    lt = sum(float(a < b) for a in x for b in y)
    return (gt - lt) / (len(x) * len(y))


def auc_pairwise(y: np.ndarray, score: np.ndarray) -> float:
    pos = score[y == 1]
    neg = score[y == 0]
    wins = sum(float(a > b) + 0.5 * float(a == b) for a in pos for b in neg)
    return wins / (len(pos) * len(neg))


def average_precision(y: np.ndarray, score: np.ndarray) -> float:
    order = np.argsort(-score, kind="mergesort")
    yy = y[order]
    tp = np.cumsum(yy)
    precision = tp / np.arange(1, len(yy) + 1)
    return float(precision[yy == 1].sum() / yy.sum())


def rankdata_average(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and values[order[j]] == values[order[i]]:
            j += 1
        ranks[order[i:j]] = (i + 1 + j) / 2
        i = j
    return ranks


def exact_mann_whitney_p(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    values = np.concatenate([x, y])
    ranks = rankdata_average(values)
    if not np.allclose(ranks, np.round(ranks)):
        raise ValueError("Exact audit implementation requires no ties")
    m, n = len(x), len(y)
    rank_sum = int(round(ranks[:m].sum()))
    u = rank_sum - m * (m + 1) // 2
    # Dynamic program for the exact rank-sum distribution of choosing m ranks from 1..m+n.
    max_sum = m * (2 * (m + n) - m + 1) // 2
    dp = [[0] * (max_sum + 1) for _ in range(m + 1)]
    dp[0][0] = 1
    for rank in range(1, m + n + 1):
        for k in range(min(rank, m), 0, -1):
            for s in range(max_sum, rank - 1, -1):
                dp[k][s] += dp[k - 1][s - rank]
    total = math.comb(m + n, m)
    lower = sum(dp[m][: rank_sum + 1]) / total
    upper = sum(dp[m][rank_sum:]) / total
    return float(u), min(1.0, 2 * min(lower, upper))


def asymptotic_mann_whitney_p_with_continuity(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Reproduce two-sided R wilcox.test(..., exact=FALSE) for independent samples."""
    values = np.concatenate([x, y])
    ranks = rankdata_average(values)
    m, n = len(x), len(y)
    u = float(ranks[:m].sum() - m * (m + 1) / 2)
    _, tie_counts = np.unique(values, return_counts=True)
    total = m + n
    tie_term = float(np.sum(tie_counts**3 - tie_counts))
    variance = m * n / 12 * ((total + 1) - tie_term / (total * (total - 1)))
    if variance <= 0:
        return u, float("nan")
    centered = u - m * n / 2
    continuity = 0.5 * np.sign(centered)
    z = (centered - continuity) / math.sqrt(variance)
    p = math.erfc(abs(z) / math.sqrt(2))
    return u, float(p)


def validation_composition() -> dict:
    path = AUDIT / "GSE191067_standard_validation" / "standard_validation_cell_metadata.csv"
    df = pd.read_csv(path)
    df["standard_cluster"] = df["standard_cluster"].astype(int)
    ambiguous = {0, 2, 6}

    rows = []
    for analysis, sub in [
        ("all_cells", df),
        ("confident_clusters_only", df.loc[~df["standard_cluster"].isin(ambiguous)].copy()),
    ]:
        counts = (
            sub.groupby(["sample", "condition", "standard_program"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        totals = counts.groupby(["sample", "condition"], observed=True)["n_cells"].transform("sum")
        counts["proportion"] = counts["n_cells"] / totals
        counts.insert(0, "analysis", analysis)
        rows.append(counts)
    comp = pd.concat(rows, ignore_index=True)
    comp = comp.rename(columns={"standard_program": "state"})
    comp.to_csv(OUT / "validation_state_composition_sensitivity.csv", index=False)

    stats = []
    for (analysis, state), g in comp.groupby(["analysis", "state"], observed=True):
        piv = g[["sample", "condition", "proportion"]].drop_duplicates()
        k = piv.loc[piv.condition.eq("Keloid"), "proportion"].to_numpy()
        n = piv.loc[piv.condition.eq("Normal"), "proportion"].to_numpy()
        if len(k) and len(n):
            _, p = exact_mann_whitney_p(k, n)
            stats.append(
                {
                    "analysis": analysis,
                    "state": state,
                    "n_keloid_donors": len(k),
                    "n_normal_donors": len(n),
                    "median_keloid": float(np.median(k)),
                    "median_normal": float(np.median(n)),
                    "cliffs_delta": cliff_delta(k, n),
                    "exact_mann_whitney_p": p,
                }
            )
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(OUT / "validation_state_composition_sensitivity_stats.csv", index=False)
    return {
        "n_all": len(df),
        "n_confident": int((~df["standard_cluster"].isin(ambiguous)).sum()),
        "stats": stats_df.to_dict(orient="records"),
    }


def bulk_audit() -> dict:
    path = AUDIT / "bulk_recurrence" / "bulk_primary_scores.csv"
    df = pd.read_csv(path)
    reported_stats = pd.read_csv(AUDIT / "bulk_recurrence" / "bulk_ROC_AP_permutation_statistics.csv")
    reported_pairs = pd.read_csv(AUDIT / "bulk_recurrence" / "bulk_paired_bootstrap_AUC_comparisons.csv")
    y = df["recurrence"].to_numpy(dtype=int)
    rows = []
    for col in ["clean_D1_loss", "clean_D2", "five_hub"]:
        s = df[col].to_numpy(float)
        rec = s[y == 1]
        non = s[y == 0]
        _, mw_p = asymptotic_mann_whitney_p_with_continuity(rec, non)
        reported = reported_stats.loc[reported_stats["score"].eq(col)].iloc[0]
        auc = float(auc_pairwise(y, s))
        ap = float(average_precision(y, s))
        rows.append(
            {
                "score": col,
                "n": len(y),
                "events": int(y.sum()),
                "median_recurrence": float(np.median(rec)),
                "median_nonrecurrence": float(np.median(non)),
                "wilcoxon_asymptotic_continuity_p": mw_p,
                "cliffs_delta": cliff_delta(rec, non),
                "roc_auc": auc,
                "average_precision": ap,
                "matches_primary_wilcoxon": bool(np.isclose(mw_p, reported["wilcoxon_p"], atol=1e-12)),
                "matches_primary_auc": bool(np.isclose(auc, reported["auc"], atol=1e-12)),
                "matches_primary_average_precision": bool(np.isclose(ap, reported["average_precision"], atol=1e-12)),
            }
        )
    stats = pd.DataFrame(rows)
    stats.to_csv(OUT / "bulk_independent_point_estimates.csv", index=False)
    pair = reported_pairs.loc[reported_pairs["comparison"].eq("clean_D1loss_vs_fivehub")].iloc[0]
    observed_delta = auc_pairwise(y, df["clean_D1_loss"].to_numpy(float)) - auc_pairwise(
        y, df["five_hub"].to_numpy(float)
    )
    checks = stats[[
        "matches_primary_wilcoxon", "matches_primary_auc", "matches_primary_average_precision"
    ]].to_numpy(dtype=bool)
    return {
        "point_estimates": stats.to_dict(orient="records"),
        "paired_stratified_bootstrap": {
            "design": "5000 paired, recurrence-stratified bootstrap replicates; seed 20260711",
            "comparison": "clean_D1_loss minus five_hub",
            "independently_recomputed_observed_auc_difference": float(observed_delta),
            "reported_observed_auc_difference": float(pair["observed_auc_difference"]),
            "reported_95pct_ci": [float(pair["bootstrap_ci_low"]), float(pair["bootstrap_ci_high"])],
            "observed_difference_matches": bool(np.isclose(observed_delta, pair["observed_auc_difference"], atol=1e-12)),
        },
        "status": "PASS" if checks.all() and np.isclose(observed_delta, pair["observed_auc_difference"], atol=1e-12) else "FAIL",
    }


def main() -> None:
    result = {
        "validation_composition": validation_composition(),
        "bulk": bulk_audit(),
    }
    (OUT / "numeric_audit_summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
