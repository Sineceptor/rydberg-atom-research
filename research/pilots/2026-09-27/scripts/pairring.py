"""Periodic ring of strong pairs, pair-site ED with translation/inversion sector penalties.
States per pair: 0=00, 1=L(10), 2=R(01), [3=11 if nloc=4]. Units Omega=a=1."""
import numpy as np, sys, time
from scipy.sparse import csr_matrix, diags
from scipy.sparse.linalg import eigsh, LinearOperator

def build(Np, eps, g, delta, nloc=3, omega=1.0):
    s = 1.0 - eps; C = 2.0 * Np
    occ = {0: (0, 0), 1: (1, 0), 2: (0, 1), 3: (1, 1)}
    dim = nloc ** Np
    idx = np.arange(dim, dtype=np.int64)
    digs = np.empty((dim, Np), dtype=np.int8)
    t = idx.copy()
    for j in range(Np):
        digs[:, j] = t % nloc; t //= nloc
    nL = np.isin(digs, [1, 3]).astype(float); nR = np.isin(digs, [2, 3]).astype(float)
    pos = lambda j, a: 2.0 * j + (s if a else 0.0)
    def V(x, y):
        d = abs(x - y) % C; d = min(d, C - d); return g / d ** 6
    E = -delta * (nL.sum(1) + nR.sum(1))
    if nloc == 4:
        E += (g / s ** 6) * (nL * nR).sum(1)
    for i in range(Np):
        for j in range(i + 1, Np):
            for a, na in ((0, nL[:, i]), (1, nR[:, i])):
                for b, nb in ((0, nL[:, j]), (1, nR[:, j])):
                    E += V(pos(i, a), pos(j, b)) * na * nb
    rows, cols = [], []
    pw = nloc ** np.arange(Np)
    trans = [(0, 1), (0, 2)] + ([(1, 3), (2, 3)] if nloc == 4 else [])
    for j in range(Np):
        dj = digs[:, j]
        for a, b in trans:
            m = dj == a
            r = idx[m]; c = r + (b - a) * pw[j]
            rows += [r, c]; cols += [c, r]
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    H = csr_matrix((np.full(len(rows), 0.5 * omega), (rows, cols)), shape=(dim, dim)) + diags(E)
    # translation: new digit j = old digit j-1
    dT = np.roll(digs, 1, axis=1)
    permT = (dT.astype(np.int64) * pw).sum(1)
    swap = np.array([0, 2, 1, 3][:nloc], dtype=np.int8)
    dI = swap[digs[:, (-np.arange(Np)) % Np]]
    permI = (dI.astype(np.int64) * pw).sum(1)
    return H, permT, permI

def sector_levels(H, permT, permI, parity, k=2, lam=50.0):
    """Lowest k levels in the k=0 sector with inversion parity (+1/-1)."""
    dim = H.shape[0]
    invT = np.empty_like(permT); invT[permT] = np.arange(dim)
    def mv(v):
        v = v.ravel()
        out = H @ v
        out += lam * (v - 0.5 * (v[permT] + v[invT]))
        out += lam * 0.5 * (v - parity * v[permI])
        return out
    op = LinearOperator((dim, dim), matvec=mv, dtype=float)
    v0 = np.random.default_rng(1).normal(size=dim)
    # symmetrize v0 into the sector
    vals = eigsh(op, k=k, which='SA', tol=1e-10, v0=v0)[0]
    return np.sort(vals)

def levels(Np, eps, g, delta, nloc=3):
    H, pT, pI = build(Np, eps, g, delta, nloc)
    ev = sector_levels(H, pT, pI, +1, k=2)
    od = sector_levels(H, pT, pI, -1, k=1)
    E0, Ee, Es = ev[0], ev[1], od[0]
    return E0, Es - E0, Ee - E0
