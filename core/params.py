from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class TraceParams:
    """cfRNA 脑区溯源分析的核心参数集。

    Attributes:
        use_value: 表达值变换方式。
            - 'log1p'（默认）：对数变换，适合大多数场景。
            - 'tpm'：原始 TPM，不做变换。
            - 'zscore'：先 log1p 再做 z-score 标准化。
        bootstrap_n: Bootstrap 重采样次数，用于估算置信区间与稳定性。
            设为 0 则跳过 bootstrap，分析速度更快。合理范围：50–200。
        bootstrap_gene_frac: 每次 bootstrap 随机采样的基因比例，范围 (0, 1]。
            默认 0.7，即每轮使用 70% 的 overlap 基因。
        l2: NNLS Simplex 优化中的 L2 正则化强度，范围 [0, ∞)。
            默认 1e-4，增大可抑制稀疏解，减小可提高拟合精度。
        topk: 返回结果的 top-K 脑区数量，默认 10。
            仅在 return_all=False 时生效。
        random_seed: Bootstrap 随机种子，保证结果可复现。
        atlas_id: 使用的参考图谱版本 ID，对应 atlas_versions 表中的 atlas_id。
            默认 1（legacy 图谱）。
        ensemble_alpha: Ensemble 模式下 corr 信号的混合权重，范围 [0, 1]。
            alpha=1 → 纯相关性；alpha=0 → 纯 NNLS；alpha=0.5（默认）为等权混合。
            其他信号（rank / marker / detect / support）的权重由 extra params 单独控制。
        return_all: 是否返回所有脑区的得分（不截断到 topk）。
            用于 benchmark 和可视化场景。
    """

    use_value: str = "log1p"
    bootstrap_n: int = 100
    bootstrap_gene_frac: float = 0.7
    l2: float = 1e-4
    topk: int = 10
    random_seed: int = 13
    atlas_id: int = 1
    ensemble_alpha: float = 0.5
    return_all: bool = False
