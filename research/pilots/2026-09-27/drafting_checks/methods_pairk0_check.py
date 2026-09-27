"""Spot check of pilot pairk0.py precision: float32 diagonal + inversion-even start vector vs float64 + random start."""
import sys, json, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
sys.path.insert(0, "/home/user/rydberg-atom-research/research/pilots/2026-09-27/scripts"); sys.dont_write_bytecode = True
import pairk0
def eint64(Np, d, eps, g):
    s = 1.0 - eps; Na = 2 * Np; L = 2.0 * Np
    x = np.empty(Na); x[0::2] = 2 * np.arange(Np); x[1::2] = 2 * np.arange(Np) + s
    dx = np.abs(x[:, None] - x[None, :]); dx = np.minimum(dx, L - dx)
    V = np.zeros((Na, Na)); off = ~np.eye(Na, dtype=bool); V[off] = g / dx[off] ** 6
    D = d ** Np; c = np.arange(D, dtype=np.int64); pw = d ** (Np - 1)
    rep = c.copy(); rot = c.copy()
    for t in range(1, Np):
        rot = rot // d + (rot % d) * pw; np.minimum(rep, rot, out=rep)
    reps = np.nonzero(rep == c)[0]
    dig = np.empty((len(reps), Np), dtype=np.int8); tmp = reps.copy()
    for j in range(Np): dig[:, j] = tmp % d; tmp //= d
    occ = np.empty((len(reps), Na)); occ[:, 0::2] = ((dig == 1) | (dig == 3)); occ[:, 1::2] = ((dig == 2) | (dig == 3))
    return 0.5 * np.einsum('ij,ij->i', occ @ V, occ)
out = []
for Np, d, eps, g, delta in [(8, 4, 0.3, 11.39, 1.3257), (10, 3, 0.3, 11.39, 1.3257), (8, 4, 0.0, 11.39, 1.0785)]:
    B = pairk0.build(Np, d, eps, g)
    e64 = eint64(Np, d, eps, g)
    gs32, ge32, _ = pairk0.gaps(B, delta)
    B64 = dict(B); B64["Eint"] = e64
    # float64 diagonal, random start vector, parity via the same projection trick
    H = B64["Hoff"] + sp.diags(B64["Eint"] - delta * B64["nsum"])
    E, V = sla.eigsh(H, k=8, which="SA", tol=1e-12, v0=np.random.default_rng(1).normal(size=B["nrep"]))
    o = np.argsort(E); E = E[o]; V = V[:, o]
    M = np.diag(E) + 1e-7 * (V.T @ V[B["pidx"], :]); w, U = np.linalg.eigh(M); W = V @ U
    En = np.einsum('ij,ij->j', W, H @ W); par = np.einsum('ij,ij->j', W, W[B["pidx"], :]); o = np.argsort(En); En, par = En[o], par[o]
    odd = [e for e, p in zip(En, par) if p < -0.5]; even = [e for e, p in zip(En, par) if p > 0.5]
    gs64, ge64 = odd[0] - En[0], even[1] - En[0]
    out.append(dict(Np=Np, d=d, eps=eps, delta=delta, nrep=int(B["nrep"]), max_abs_Eint_err=float(np.max(abs(B["Eint"] - e64))),
                    gap_sigma_pilot=gs32, gap_sigma_f64=gs64, R_pilot=ge32 / gs32, R_f64=ge64 / gs64,
                    d_gap_sigma=gs32 - gs64, d_R=ge32 / gs32 - ge64 / gs64))
print(json.dumps(out, indent=1))
