from __future__ import annotations

from typing import Tuple
import numpy as np
from scipy.optimize import minimize


def apply_value_transform(x: np.ndarray, use_value: str) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if use_value in {'tpm', 'vsd'}:
        return x
    if use_value == 'log1p':
        return np.log1p(np.clip(x, 0, None))
    if use_value == 'zscore':
        y = np.log1p(np.clip(x, 0, None))
        if y.ndim == 2:
            return (y - y.mean(axis=0, keepdims=True)) / (y.std(axis=0, keepdims=True) + 1e-8)
        return (y - y.mean()) / (y.std() + 1e-8)
    return x


def zscore_safe(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    return (v - v.mean()) / (v.std() + 1e-8)


def softmax_confidence(scores: np.ndarray) -> np.ndarray:
    s = np.asarray(scores, dtype=float)
    s = s - np.max(s)
    p = np.exp(s)
    return p / (p.sum() + 1e-12)


def trace_corr(A: np.ndarray, b: np.ndarray) -> np.ndarray:
    b0 = b - b.mean()
    A0 = A - A.mean(axis=0, keepdims=True)
    num = (A0 * b0[:, None]).sum(axis=0)
    den = np.sqrt((A0 ** 2).sum(axis=0) * (b0 ** 2).sum() + 1e-12)
    corr = num / den
    return np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)


def trace_nnls_simplex(A: np.ndarray, b: np.ndarray, l2: float) -> Tuple[np.ndarray, float]:
    G, R = A.shape
    if R == 1:
        w = np.array([1.0], dtype=float)
        rmse = float(np.sqrt(np.mean((A @ w - b) ** 2)))
        return w, rmse

    w0 = np.full(R, 1.0 / R, dtype=float)

    def obj(w: np.ndarray) -> float:
        r = A @ w - b
        return 0.5 * float(r @ r) + 0.5 * float(l2) * float(w @ w)

    def grad(w: np.ndarray) -> np.ndarray:
        return A.T @ (A @ w - b) + float(l2) * w

    cons = ({'type': 'eq', 'fun': lambda w: np.sum(w) - 1.0},)
    bounds = [(0.0, 1.0) for _ in range(R)]
    res = minimize(obj, w0, method='SLSQP', jac=grad, bounds=bounds, constraints=cons, options={'maxiter': 300, 'ftol': 1e-10, 'disp': False})
    w = np.clip(res.x if res.success else w0, 0.0, None)
    s = w.sum()
    w = w / s if s > 0 else w0
    rmse = float(np.sqrt(np.mean((A @ w - b) ** 2)))
    return w.astype(float), rmse
