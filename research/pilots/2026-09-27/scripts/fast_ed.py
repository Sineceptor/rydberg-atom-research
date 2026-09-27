"""Vectorized ED for arbitrary 1D positions (pilot; bit i = atom i)."""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

def positions(N, eps):
    gaps = [(1 - eps) if i % 2 == 0 else (1 + eps) for i in range(N - 1)]
    return np.concatenate([[0.0], np.cumsum(gaps)])

def build(x, g, delta, omega=1.0):
    N = len(x); D = 2**N
    idx = np.arange(D, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(N)) & 1).astype(np.float64)
    diag = -delta * bits.sum(1)
    for i in range(N):
        for j in range(i + 1, N):
            diag += g / abs(x[j] - x[i])**6 * bits[:, i] * bits[:, j]
    rows = np.repeat(idx, N); cols = (idx[:, None] ^ (1 << np.arange(N))).ravel()
    H = csr_matrix((np.full(rows.size, omega / 2), (rows, cols)), shape=(D, D))
    H = H + csr_matrix((diag, (idx, idx)), shape=(D, D))
    return H, bits

def ground(x, g, delta, k=2):
    H, bits = build(x, g, delta)
    v0 = np.random.default_rng(7).normal(size=H.shape[0])
    E, V = eigsh(H, k=k, which="SA", tol=1e-10, v0=v0)
    o = np.argsort(E); E, V = E[o], V[:, o]
    w = np.abs(V[:, 0])**2
    return E, w @ bits
