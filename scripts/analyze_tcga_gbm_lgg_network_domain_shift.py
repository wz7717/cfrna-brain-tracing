from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = ROOT / "results" / "tcga_gbm_lgg_sample_mri_label_tracing_20260605"
NETWORK_FILE = INPUT_DIR / "tcga_gbm_lgg_sample_network_tracing.csv"
SUMMARY_FILE = INPUT_DIR / "tcga_gbm_lgg_sample_mri_label_tracing_summary.csv"
OUTDIR = ROOT / "results" / "tcga_gbm_vs_lgg_network_domain_shift_20260609"
RNG = np.random.default_rng(20260609)
BOOTSTRAP_N = 5000
PERMUTATION_N = 10000
LOW_CONF_THRESHOLDS = (0.001, 0.002, 0.005)
GROUP_ORDER = ("TCGA-GBM", "TCGA-LGG")


def configure_plots() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Microsoft YaHei", "SimHei", "Arial", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "figure.dpi": 140,
            "savefig.dpi": 220,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def bh_fdr(p_values: pd.Series) -> pd.Series:
    values = p_values.to_numpy(dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.clip(adjusted, 0, 1)
    return pd.Series(output, index=p_values.index)


def cliffs_delta(x: np.ndarray, y: np.ndarray) -> float:
    y_sorted = np.sort(y)
    less = np.searchsorted(y_sorted, x, side="left").sum()
    greater = (len(y_sorted) - np.searchsorted(y_sorted, x, side="right")).sum()
    return float((less - greater) / (len(x) * len(y)))


def bootstrap_median_difference(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    differences = np.empty(BOOTSTRAP_N, dtype=float)
    for index in range(BOOTSTRAP_N):
        bx = RNG.choice(x, size=len(x), replace=True)
        by = RNG.choice(y, size=len(y), replace=True)
        differences[index] = np.median(bx) - np.median(by)
    estimate = float(np.median(x) - np.median(y))
    low, high = np.quantile(differences, [0.025, 0.975])
    return estimate, float(low), float(high)


def odds_ratio_ci(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    cells = np.array([a, b, c, d], dtype=float)
    if np.any(cells == 0):
        cells += 0.5
    a2, b2, c2, d2 = cells
    odds_ratio = (a2 * d2) / (b2 * c2)
    se = np.sqrt(np.sum(1.0 / cells))
    low = np.exp(np.log(odds_ratio) - 1.96 * se)
    high = np.exp(np.log(odds_ratio) + 1.96 * se)
    return float(odds_ratio), float(low), float(high)


def cramers_v(table: np.ndarray) -> float:
    chi2 = stats.chi2_contingency(table, correction=False)[0]
    n = table.sum()
    return float(np.sqrt((chi2 / n) / min(table.shape[0] - 1, table.shape[1] - 1)))


def js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    p = p / p.sum()
    q = q / q.sum()
    midpoint = 0.5 * (p + q)
    return float(
        0.5 * stats.entropy(p, midpoint, base=2)
        + 0.5 * stats.entropy(q, midpoint, base=2)
    )


def permutation_distribution_pvalue(
    top1_networks: np.ndarray, group_labels: np.ndarray, networks: list[str]
) -> float:
    observed = pd.crosstab(group_labels, top1_networks).reindex(
        index=GROUP_ORDER, columns=networks, fill_value=0
    )
    observed_chi2 = stats.chi2_contingency(observed.to_numpy(), correction=False)[0]
    exceed = 0
    for _ in range(PERMUTATION_N):
        permuted = RNG.permutation(group_labels)
        table = pd.crosstab(permuted, top1_networks).reindex(
            index=GROUP_ORDER, columns=networks, fill_value=0
        )
        chi2 = stats.chi2_contingency(table.to_numpy(), correction=False)[0]
        exceed += chi2 >= observed_chi2
    return float((exceed + 1) / (PERMUTATION_N + 1))


def make_level_table(long_df: pd.DataFrame, level: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    id_column = "sample_id" if level == "sample" else "patient_barcode"
    if level == "patient":
        grouped = (
            long_df.groupby(["patient_barcode", "project_id", "network_id"], as_index=False)
            .agg(score=("score", "mean"), probability=("confidence", "mean"))
        )
    else:
        grouped = long_df.rename(columns={"confidence": "probability"}).copy()

    probability_sums = grouped.groupby(id_column)["probability"].transform("sum")
    grouped["probability"] = grouped["probability"] / probability_sums
    grouped["rank_probability"] = grouped.groupby(id_column)["probability"].rank(
        method="first", ascending=False
    )
    grouped["rank_score"] = grouped.groupby(id_column)["score"].rank(
        method="first", ascending=False
    )

    rows: list[dict[str, object]] = []
    for entity_id, frame in grouped.groupby(id_column, sort=False):
        frame = frame.sort_values("rank_probability")
        probabilities = frame["probability"].to_numpy(dtype=float)
        top = frame.iloc[0]
        second = frame.iloc[1]
        score_ranked = frame.sort_values("rank_score")
        entropy = stats.entropy(probabilities)
        normalized_entropy = entropy / np.log(len(probabilities))
        rows.append(
            {
                id_column: entity_id,
                "project_id": top["project_id"],
                "top1_network": top["network_id"],
                "top2_network": second["network_id"],
                "top1_probability": float(top["probability"]),
                "top2_probability": float(second["probability"]),
                "probability_margin": float(top["probability"] - second["probability"]),
                "raw_score_margin": float(
                    score_ranked.iloc[0]["score"] - score_ranked.iloc[1]["score"]
                ),
                "normalized_entropy": float(normalized_entropy),
                "effective_network_count": float(np.exp(entropy)),
                "low_conf_margin_lt_0p001": bool(
                    top["probability"] - second["probability"] < 0.001
                ),
                "low_conf_margin_lt_0p002": bool(
                    top["probability"] - second["probability"] < 0.002
                ),
                "low_conf_margin_lt_0p005": bool(
                    top["probability"] - second["probability"] < 0.005
                ),
            }
        )
    metrics = pd.DataFrame(rows)
    return grouped, metrics


def continuous_tests(metrics: pd.DataFrame, level: str) -> pd.DataFrame:
    rows = []
    for metric in (
        "normalized_entropy",
        "effective_network_count",
        "probability_margin",
        "raw_score_margin",
    ):
        gbm = metrics.loc[metrics["project_id"] == "TCGA-GBM", metric].to_numpy()
        lgg = metrics.loc[metrics["project_id"] == "TCGA-LGG", metric].to_numpy()
        mann = stats.mannwhitneyu(gbm, lgg, alternative="two-sided")
        difference, ci_low, ci_high = bootstrap_median_difference(gbm, lgg)
        rows.append(
            {
                "level": level,
                "metric": metric,
                "gbm_n": len(gbm),
                "lgg_n": len(lgg),
                "gbm_median": float(np.median(gbm)),
                "gbm_q1": float(np.quantile(gbm, 0.25)),
                "gbm_q3": float(np.quantile(gbm, 0.75)),
                "lgg_median": float(np.median(lgg)),
                "lgg_q1": float(np.quantile(lgg, 0.25)),
                "lgg_q3": float(np.quantile(lgg, 0.75)),
                "median_difference_gbm_minus_lgg": difference,
                "bootstrap_ci_low": ci_low,
                "bootstrap_ci_high": ci_high,
                "mann_whitney_u": float(mann.statistic),
                "p_value": float(mann.pvalue),
                "cliffs_delta_gbm_vs_lgg": cliffs_delta(gbm, lgg),
            }
        )
    result = pd.DataFrame(rows)
    result["p_fdr"] = bh_fdr(result["p_value"])
    return result


def low_confidence_tests(metrics: pd.DataFrame, level: str) -> pd.DataFrame:
    rows = []
    pooled_p10 = float(metrics["probability_margin"].quantile(0.10))
    threshold_specs = [(f"{value:g}", value) for value in LOW_CONF_THRESHOLDS]
    threshold_specs.append(("pooled_P10", pooled_p10))
    for label, threshold in threshold_specs:
        if label == "pooled_P10":
            low_flag = metrics["probability_margin"] < threshold
        else:
            column = f"low_conf_margin_lt_{str(threshold).replace('.', 'p')}"
            low_flag = metrics[column]
        gbm = low_flag.loc[metrics["project_id"] == "TCGA-GBM"]
        lgg = low_flag.loc[metrics["project_id"] == "TCGA-LGG"]
        a, b = int(gbm.sum()), int((~gbm).sum())
        c, d = int(lgg.sum()), int((~lgg).sum())
        fisher = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
        odds_ratio, ci_low, ci_high = odds_ratio_ci(a, b, c, d)
        rows.append(
            {
                "level": level,
                "threshold_label": label,
                "threshold": threshold,
                "gbm_low_n": a,
                "gbm_total": len(gbm),
                "gbm_low_fraction": a / len(gbm),
                "lgg_low_n": c,
                "lgg_total": len(lgg),
                "lgg_low_fraction": c / len(lgg),
                "odds_ratio_gbm_vs_lgg": odds_ratio,
                "odds_ratio_ci_low": ci_low,
                "odds_ratio_ci_high": ci_high,
                "fisher_p_value": float(fisher.pvalue),
            }
        )
    result = pd.DataFrame(rows)
    result["fisher_p_fdr"] = bh_fdr(result["fisher_p_value"])
    return result


def network_distribution_tests(
    metrics: pd.DataFrame, networks: list[str], level: str
) -> tuple[pd.DataFrame, dict[str, float]]:
    table = pd.crosstab(metrics["project_id"], metrics["top1_network"]).reindex(
        index=GROUP_ORDER, columns=networks, fill_value=0
    )
    active_networks = table.columns[table.sum(axis=0) > 0].tolist()
    active_table = table[active_networks]
    chi2, p_value, dof, _ = stats.chi2_contingency(
        active_table.to_numpy(), correction=False
    )
    gbm_counts = active_table.loc["TCGA-GBM"].to_numpy()
    lgg_counts = active_table.loc["TCGA-LGG"].to_numpy()
    overall = {
        "level": level,
        "chi_square": float(chi2),
        "degrees_of_freedom": int(dof),
        "chi_square_p_value": float(p_value),
        "permutation_p_value": permutation_distribution_pvalue(
            metrics["top1_network"].to_numpy(),
            metrics["project_id"].to_numpy(),
            active_networks,
        ),
        "cramers_v": cramers_v(active_table.to_numpy()),
        "jensen_shannon_divergence_bits": js_divergence(gbm_counts, lgg_counts),
        "active_top1_network_count": len(active_networks),
    }

    rows = []
    gbm_total = int(table.loc["TCGA-GBM"].sum())
    lgg_total = int(table.loc["TCGA-LGG"].sum())
    for network in networks:
        a = int(table.loc["TCGA-GBM", network])
        c = int(table.loc["TCGA-LGG", network])
        b, d = gbm_total - a, lgg_total - c
        fisher = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
        if a == 0 and c == 0:
            odds_ratio, ci_low, ci_high = np.nan, np.nan, np.nan
        else:
            odds_ratio, ci_low, ci_high = odds_ratio_ci(a, b, c, d)
        rows.append(
            {
                "level": level,
                "network_id": network,
                "gbm_top1_n": a,
                "gbm_total": gbm_total,
                "gbm_top1_fraction": a / gbm_total,
                "lgg_top1_n": c,
                "lgg_total": lgg_total,
                "lgg_top1_fraction": c / lgg_total,
                "fraction_difference_gbm_minus_lgg": a / gbm_total - c / lgg_total,
                "odds_ratio_gbm_vs_lgg": odds_ratio,
                "odds_ratio_ci_low": ci_low,
                "odds_ratio_ci_high": ci_high,
                "fisher_p_value": float(fisher.pvalue),
            }
        )
    result = pd.DataFrame(rows)
    result["fisher_p_fdr"] = bh_fdr(result["fisher_p_value"])
    return result, overall


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(OUTDIR / f"{stem}.png", bbox_inches="tight")
    fig.savefig(OUTDIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_top1_distribution(distribution: pd.DataFrame, level: str) -> None:
    frame = distribution.sort_values("gbm_top1_fraction", ascending=True)
    y = np.arange(len(frame))
    fig, ax = plt.subplots(figsize=(10, 6.5))
    width = 0.38
    ax.barh(y - width / 2, frame["gbm_top1_fraction"], height=width, label="GBM", color="#B43A3A")
    ax.barh(y + width / 2, frame["lgg_top1_fraction"], height=width, label="LGG", color="#287A78")
    ax.set_yticks(y, frame["network_id"])
    ax.set_xlabel("Top1 proportion")
    ax.set_title(f"GBM vs LGG Network Top1 distribution ({level} level)")
    ax.legend(frameon=False)
    ax.xaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    save_figure(fig, f"{level}_top1_network_distribution")


def plot_metric_box(metrics: pd.DataFrame, metric: str, ylabel: str, stem: str) -> None:
    gbm = metrics.loc[metrics["project_id"] == "TCGA-GBM", metric].to_numpy()
    lgg = metrics.loc[metrics["project_id"] == "TCGA-LGG", metric].to_numpy()
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    boxes = ax.boxplot(
        [gbm, lgg],
        tick_labels=["GBM", "LGG"],
        patch_artist=True,
        widths=0.55,
        showfliers=False,
    )
    for patch, color in zip(boxes["boxes"], ["#B43A3A", "#287A78"]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    jitter_gbm = RNG.normal(1, 0.045, size=len(gbm))
    jitter_lgg = RNG.normal(2, 0.045, size=len(lgg))
    ax.scatter(jitter_gbm, gbm, s=8, alpha=0.20, color="#6E2020", linewidths=0)
    ax.scatter(jitter_lgg, lgg, s=8, alpha=0.20, color="#164E4C", linewidths=0)
    ax.set_ylabel(ylabel)
    ax.set_title(f"Patient-level {ylabel}")
    save_figure(fig, stem)


def plot_low_confidence(low_conf: pd.DataFrame) -> None:
    frame = low_conf.loc[low_conf["level"] == "patient"].copy()
    x = np.arange(len(frame))
    width = 0.36
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.bar(
        x - width / 2,
        frame["gbm_low_fraction"],
        width=width,
        label="GBM",
        color="#B43A3A",
    )
    ax.bar(
        x + width / 2,
        frame["lgg_low_fraction"],
        width=width,
        label="LGG",
        color="#287A78",
    )
    labels = [
        "pooled P10" if label == "pooled_P10" else f"margin < {value:g}"
        for label, value in zip(frame["threshold_label"], frame["threshold"])
    ]
    ax.set_xticks(x, labels)
    ax.set_ylabel("Low-confidence proportion")
    ax.yaxis.set_major_formatter(lambda value, _: f"{value:.0%}")
    ax.set_title("Low-confidence threshold sensitivity")
    ax.legend(frameon=False)
    save_figure(fig, "patient_low_confidence_threshold_sensitivity")


def plot_mean_probability_heatmap(patient_long: pd.DataFrame, networks: list[str]) -> None:
    matrix = (
        patient_long.groupby(["project_id", "network_id"])["probability"]
        .mean()
        .unstack()
        .reindex(index=GROUP_ORDER, columns=networks)
    )
    fig, ax = plt.subplots(figsize=(11, 3.3))
    image = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="viridis")
    ax.set_yticks(np.arange(len(matrix.index)), matrix.index)
    ax.set_xticks(np.arange(len(networks)), networks, rotation=35, ha="right")
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix.iloc[row, column]
            ax.text(
                column,
                row,
                f"{value:.3f}",
                ha="center",
                va="center",
                fontsize=7,
                color="white" if value > matrix.to_numpy().mean() else "black",
            )
    fig.colorbar(image, ax=ax, label="Mean Network probability")
    ax.set_title("Patient-level mean Network probability")
    save_figure(fig, "patient_mean_network_probability_heatmap")


def plot_patient_sample_consistency(
    patient_tests: pd.DataFrame, sample_tests: pd.DataFrame
) -> None:
    patient = patient_tests.set_index("metric")["cliffs_delta_gbm_vs_lgg"]
    sample = sample_tests.set_index("metric")["cliffs_delta_gbm_vs_lgg"]
    metrics = list(patient.index)
    x = np.arange(len(metrics))
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(patient.loc[metrics], sample.loc[metrics], s=70, color="#355C9A")
    for metric, px, sx in zip(metrics, patient.loc[metrics], sample.loc[metrics]):
        ax.annotate(metric, (px, sx), xytext=(5, 5), textcoords="offset points", fontsize=8)
    low = min(patient.min(), sample.min()) - 0.03
    high = max(patient.max(), sample.max()) + 0.03
    ax.plot([low, high], [low, high], linestyle="--", color="#777777")
    ax.axhline(0, color="#BBBBBB", linewidth=0.8)
    ax.axvline(0, color="#BBBBBB", linewidth=0.8)
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    ax.set_xlabel("Patient-level Cliff's delta")
    ax.set_ylabel("Sample-level Cliff's delta")
    ax.set_title("Patient- and sample-level effect direction consistency")
    save_figure(fig, "patient_sample_effect_consistency")


def fmt_pct(value: float) -> str:
    return f"{100 * value:.1f}%"


def write_report(
    patient_metrics: pd.DataFrame,
    continuous: pd.DataFrame,
    low_conf: pd.DataFrame,
    distribution: pd.DataFrame,
    overall: dict[str, float],
) -> None:
    patient_tests = continuous.loc[continuous["level"] == "patient"].set_index("metric")
    patient_low = low_conf.loc[low_conf["level"] == "patient"].set_index("threshold")
    significant = distribution.loc[
        (distribution["level"] == "patient") & (distribution["fisher_p_fdr"] < 0.05)
    ].sort_values("fisher_p_fdr")
    top_counts = (
        patient_metrics.groupby(["project_id", "top1_network"])
        .size()
        .rename("n")
        .reset_index()
    )
    dominant = {}
    for group in GROUP_ORDER:
        row = top_counts.loc[top_counts["project_id"] == group].sort_values("n", ascending=False).iloc[0]
        dominant[group] = (row["top1_network"], int(row["n"]))

    entropy = patient_tests.loc["normalized_entropy"]
    margin = patient_tests.loc["probability_margin"]
    threshold = patient_low.loc[0.002]
    lines = [
        "# TCGA-GBM 与 TCGA-LGG Network 预测域偏移分析",
        "",
        "## 分析目的",
        "",
        "比较 GBM 与 LGG 的 Network 预测分布、预测熵、Top1/Top2 margin 和低置信度比例。"
        "本分析未使用 MRI 肿瘤位置真值，只描述疾病类型相关的表达域偏移，不评价脑区定位准确率。",
        "",
        "## 数据",
        "",
        f"- 患者级主分析：GBM {sum(patient_metrics.project_id == 'TCGA-GBM')} 人，"
        f"LGG {sum(patient_metrics.project_id == 'TCGA-LGG')} 人。",
        "- 样本级敏感性分析：GBM 285 个样本，LGG 516 个样本。",
        "- 每个样本均包含 10 个 Network 的 score 和归一化 confidence。",
        "",
        "## 主要结果",
        "",
        f"1. 患者级 Top1 Network 总体分布差异：置换检验 p={overall['permutation_p_value']:.4g}，"
        f"Cramér's V={overall['cramers_v']:.3f}，"
        f"Jensen-Shannon divergence={overall['jensen_shannon_divergence_bits']:.4f} bits。",
        f"2. GBM 最常见 Top1 为 {dominant['TCGA-GBM'][0]}（n={dominant['TCGA-GBM'][1]}）；"
        f"LGG 最常见 Top1 为 {dominant['TCGA-LGG'][0]}（n={dominant['TCGA-LGG'][1]}）。",
        f"3. 归一化熵中位数：GBM {entropy['gbm_median']:.4f}，LGG {entropy['lgg_median']:.4f}；"
        f"Cliff's delta={entropy['cliffs_delta_gbm_vs_lgg']:.3f}，FDR={entropy['p_fdr']:.4g}。",
        f"4. probability margin 中位数：GBM {margin['gbm_median']:.5f}，"
        f"LGG {margin['lgg_median']:.5f}；Cliff's delta={margin['cliffs_delta_gbm_vs_lgg']:.3f}，"
        f"FDR={margin['p_fdr']:.4g}。",
        f"5. margin<0.002 的低置信度比例：GBM {fmt_pct(threshold['gbm_low_fraction'])}，"
        f"LGG {fmt_pct(threshold['lgg_low_fraction'])}；"
        f"OR={threshold['odds_ratio_gbm_vs_lgg']:.3f} "
        f"(95% CI {threshold['odds_ratio_ci_low']:.3f}–{threshold['odds_ratio_ci_high']:.3f})，"
        f"FDR={threshold['fisher_p_fdr']:.4g}。",
        "",
        "## 逐 Network 差异",
        "",
    ]
    if significant.empty:
        lines.append("患者级逐 Network Fisher 检验经 BH-FDR 校正后未发现显著差异。")
    else:
        lines.extend(
            [
                "| Network | GBM比例 | LGG比例 | OR | FDR |",
                "| --- | ---: | ---: | ---: | ---: |",
            ]
        )
        for _, row in significant.iterrows():
            lines.append(
                f"| {row['network_id']} | {fmt_pct(row['gbm_top1_fraction'])} | "
                f"{fmt_pct(row['lgg_top1_fraction'])} | {row['odds_ratio_gbm_vs_lgg']:.3f} | "
                f"{row['fisher_p_fdr']:.4g} |"
            )

    lines.extend(
        [
            "",
            "## 解释边界",
            "",
            "- 差异表示 GBM 与 LGG bulk RNA-seq 投射到猕猴参考 Network 时的分布、集中度或置信度不同。",
            "- 可能来源包括肿瘤等级、肿瘤纯度、免疫浸润、坏死、增殖状态和跨物种表达域偏移。",
            "- 结果不能解释为 GBM/LGG 起源于某个 Network，也不能说明某组定位更准确。",
            "- 定位准确率必须等待同患者 MRI 分割和 atlas 配准真值。",
            "",
            "## 输出",
            "",
            "- `patient_level_metrics.csv`",
            "- `sample_level_metrics.csv`",
            "- `network_distribution_comparison.csv`",
            "- `group_statistical_tests.csv`",
            "- `low_confidence_threshold_sensitivity.csv`",
            "- `analysis_summary.json`",
            "- PNG/PDF 配套图表",
        ]
    )
    (OUTDIR / "tcga_gbm_vs_lgg_domain_shift_report_cn.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    configure_plots()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    network = pd.read_csv(NETWORK_FILE)
    summary = pd.read_csv(SUMMARY_FILE)[["sample_id", "patient_barcode", "project_id"]]
    network = network.drop(columns=["patient_barcode"], errors="ignore").merge(
        summary, on="sample_id", how="left", validate="many_to_one"
    )
    if network["project_id"].isna().any():
        raise ValueError("Some network rows lack project labels.")
    candidate_counts = network.groupby("sample_id")["network_id"].nunique()
    if candidate_counts.nunique() != 1 or candidate_counts.iloc[0] != 10:
        raise ValueError("Expected exactly 10 Network candidates per sample.")

    networks = (
        network.groupby("network_id")["confidence"].mean().sort_values(ascending=False).index.tolist()
    )
    sample_long, sample_metrics = make_level_table(network, "sample")
    patient_long, patient_metrics = make_level_table(network, "patient")

    sample_metrics.to_csv(OUTDIR / "sample_level_metrics.csv", index=False)
    patient_metrics.to_csv(OUTDIR / "patient_level_metrics.csv", index=False)
    patient_long.to_csv(OUTDIR / "patient_level_network_probabilities.csv", index=False)

    continuous_sample = continuous_tests(sample_metrics, "sample")
    continuous_patient = continuous_tests(patient_metrics, "patient")
    continuous = pd.concat([continuous_patient, continuous_sample], ignore_index=True)
    continuous.to_csv(OUTDIR / "group_statistical_tests.csv", index=False)

    low_sample = low_confidence_tests(sample_metrics, "sample")
    low_patient = low_confidence_tests(patient_metrics, "patient")
    low_conf = pd.concat([low_patient, low_sample], ignore_index=True)
    low_conf.to_csv(OUTDIR / "low_confidence_threshold_sensitivity.csv", index=False)

    distribution_patient, overall_patient = network_distribution_tests(
        patient_metrics, networks, "patient"
    )
    distribution_sample, overall_sample = network_distribution_tests(
        sample_metrics, networks, "sample"
    )
    distribution = pd.concat([distribution_patient, distribution_sample], ignore_index=True)
    distribution.to_csv(OUTDIR / "network_distribution_comparison.csv", index=False)

    plot_top1_distribution(distribution_patient, "patient")
    plot_top1_distribution(distribution_sample, "sample")
    plot_metric_box(
        patient_metrics,
        "normalized_entropy",
        "Normalized Shannon entropy",
        "patient_normalized_entropy",
    )
    plot_metric_box(
        patient_metrics,
        "probability_margin",
        "Top1-Top2 probability margin",
        "patient_probability_margin",
    )
    plot_metric_box(
        patient_metrics,
        "raw_score_margin",
        "Top1-Top2 raw score margin",
        "patient_raw_score_margin",
    )
    plot_low_confidence(low_conf)
    plot_mean_probability_heatmap(patient_long, networks)
    plot_patient_sample_consistency(continuous_patient, continuous_sample)

    summary_json = {
        "analysis_scope": (
            "Disease-domain shift only. No MRI truth labels were used and no localization "
            "accuracy is reported."
        ),
        "input_network_file": str(NETWORK_FILE),
        "n_networks": len(networks),
        "network_order": networks,
        "sample_counts": sample_metrics["project_id"].value_counts().to_dict(),
        "patient_counts": patient_metrics["project_id"].value_counts().to_dict(),
        "duplicate_patient_samples": (
            summary.groupby("patient_barcode").size().loc[lambda x: x > 1].to_dict()
        ),
        "patient_distribution_test": overall_patient,
        "sample_distribution_test": overall_sample,
        "low_confidence_primary_threshold": 0.002,
        "bootstrap_iterations": BOOTSTRAP_N,
        "permutation_iterations": PERMUTATION_N,
    }
    (OUTDIR / "analysis_summary.json").write_text(
        json.dumps(summary_json, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    write_report(
        patient_metrics,
        continuous,
        low_conf,
        distribution,
        overall_patient,
    )
    print(json.dumps(summary_json, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
