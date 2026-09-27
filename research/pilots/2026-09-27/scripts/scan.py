import numpy as np, sys, time, json
from scipy.sparse import csr_matrix, identity
from scipy.sparse.linalg import eigsh
from scipy.optimize import brentq
sys.path.insert(0, '.')
from pairring import build

def perm_matrix(p):
    n = len(p); return csr_matrix((np.ones(n), (np.arange(n), p)), shape=(n, n))

def levels(Np, eps, g, delta, nloc=3, lam=40.0):
    H, pT, pI = build(Np, eps, g, delta, nloc)
    dim = H.shape[0]
    T = perm_matrix(pT); I = perm_matrix(pI); Id = identity(dim, format='csr')
    Pk = lam * (Id - 0.5 * (T + T.T))
    out = []
    for par, k in ((+1, 2), (-1, 1)):
        M = (H + Pk + lam * 0.5 * (Id - par * I)).tocsr()
        v0 = np.random.default_rng(1).normal(size=dim)
        out.append(np.sort(eigsh(M, k=k, which='SA', tol=1e-11, v0=v0)[0]))
    E0, Ee = out[0]; Es = out[1][0]
    return E0, Es - E0, Ee - E0

def point(eps, g, lo, hi, sizes=(8, 10), nloc=3):
    a, b = sizes
    f = lambda d: a * levels(a, eps, g, d, nloc)[1] - b * levels(b, eps, g, d, nloc)[1]
    fl, fh = f(lo), f(hi)
    res = {'eps': eps, 'g': g, 'lo': lo, 'hi': hi, 'f_lo': fl, 'f_hi': fh}
    if fl * fh > 0:
        res['status'] = 'no crossing'; return res
    dc = brentq(f, lo, hi, xtol=2e-3)
    res['delta_c'] = dc
    for N in (6, 8, 10):
        E0, gs, ge = levels(N, eps, g, dc, nloc)
        res[f'N{N}'] = {'Ngap_sigma': N * gs, 'Ngap_eps': N * ge, 'ratio': ge / gs, 'E0perpair': E0 / N}
    return res

if __name__ == '__main__':
    eps, g, lo, hi = map(float, sys.argv[1:5])
    t = time.time()
    r = point(eps, g, lo, hi)
    r['sec'] = time.time() - t
    print(json.dumps(r))
