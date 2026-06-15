from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List


@dataclass(frozen=True)
class MethodSpec:
    code: str
    label: str
    short_label: str
    description: str
    best_for: str


METHOD_SPECS: Dict[str, MethodSpec] = {
    "ensemble": MethodSpec(
        code="ensemble",
        label="Multi-signal ensemble",
        short_label="Ensemble",
        description="Integrates correlation, NNLS/simplex fractions, rank concordance, signature/marker evidence, detection and support signals.",
        best_for="Default publish-grade analysis and benchmark reporting when reference/signature coverage is sufficient.",
    ),
    "correlation": MethodSpec(
        code="correlation",
        label="Correlation scoring",
        short_label="Correlation",
        description="Ranks regions by expression-profile similarity between the sample and each reference region.",
        best_for="Fast screening, sanity checks, and cases where contribution fractions are not required.",
    ),
    "nnls_simplex": MethodSpec(
        code="nnls_simplex",
        label="NNLS simplex deconvolution",
        short_label="NNLS simplex",
        description="Fits non-negative region fractions constrained to a simplex-like mixture and reports reconstruction error.",
        best_for="Mixture-style interpretation, fraction estimates, and bootstrap confidence intervals.",
    ),
    "marker_gene": MethodSpec(
        code="marker_gene",
        label="Marker-gene scoring",
        short_label="Marker genes",
        description="Uses region/signature marker genes as targeted evidence for source-region support.",
        best_for="Exploratory interpretation and signature hit inspection; not the preferred standalone benchmark method.",
    ),
}


METHOD_ALIASES = {
    "multi_signal": "ensemble",
    "multisignal": "ensemble",
    "integrated": "ensemble",
    "integrated_tracing": "ensemble",
    "ensemble_v2": "ensemble",
    "corr": "correlation",
    "pearson": "correlation",
    "spearman": "correlation",
    "nnls": "nnls_simplex",
    "simplex": "nnls_simplex",
    "nmf": "nnls_simplex",
    "least_squares": "nnls_simplex",
    "ls": "nnls_simplex",
    "marker": "marker_gene",
    "markers": "marker_gene",
    "marker_genes": "marker_gene",
}


def canonical_method(method: str | None) -> str:
    raw = (method or "ensemble").strip().lower().replace("-", "_").replace(" ", "_")
    return METHOD_ALIASES.get(raw, raw)


def method_label(method: str | None, short: bool = False) -> str:
    code = canonical_method(method)
    spec = METHOD_SPECS.get(code)
    if spec is None:
        return code
    return spec.short_label if short else spec.label


def method_choices(include_marker: bool = False) -> List[str]:
    codes = ["ensemble", "correlation", "nnls_simplex"]
    if include_marker:
        codes.append("marker_gene")
    return codes


def method_help_markdown(codes: Iterable[str] | None = None) -> str:
    rows = []
    for code in (codes or method_choices(include_marker=True)):
        spec = METHOD_SPECS[code]
        rows.append(f"- `{spec.code}`: **{spec.label}**. {spec.description} Recommended use: {spec.best_for}")
    return "\n".join(rows)
