from __future__ import annotations

from typing import Dict, List, Tuple
import numpy as np
from .models import trace_nnls_simplex


def bootstrap_nnls(A: np.ndarray, b: np.ndarray, regions: List[str], n: int, gene_frac: float, l2: float, seed: int):
    rng = np.random.default_rng(seed)
    G, R = A.shape
    m = max(10, int(G * gene_frac))
    W = np.zeros((n, R), dtype=float)
    top1 = np.zeros(n, dtype=int)
    for i in range(n):
        idx = rng.choice(G, size=m, replace=False)
        w_i, _ = trace_nnls_simplex(A[idx, :], b[idx], l2=l2)
        W[i, :] = w_i
        top1[i] = int(np.argmax(w_i))
    ci = {r: (float(np.quantile(W[:, j], 0.025)), float(np.quantile(W[:, j], 0.975))) for j, r in enumerate(regions)}
    stability = {r: float(np.mean(top1 == j)) for j, r in enumerate(regions)}
    return ci, stability
