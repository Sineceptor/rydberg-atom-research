"""Constrained ED on a periodic dimerized ring: no two excitations within 2 sites.
Diagnostics for Z4 -> PDW (pair density wave) splitting."""
import numpy as np, sys, time
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

def basis(N, excl=2):
    # enumerate bitstrings with cyclic gap > excl between ones
    # build linear strings with states tracked: config int, last-one position
    confs = np.array([0], dtype=np.int64); last = np.array([-100], dtype=np.int64); first = np.array([-1], dtype=np.int64)
    for i in range(N):
        # option 0
        c0, l0, f0 = confs, last, first
        ok = (i - last) > excl
        c1 = confs[ok] | (np.int64(1) << i); l1 = np.full(ok.sum(), i); f1 = np.where(first[ok] < 0, i, first[ok])
        confs = np.concatenate([c0, c1]); last = np.concatenate([l0, l1]); first = np.concatenate([f0, f1])
    # cyclic check: gap between last and first + N > excl
    ok = (first < 0) | ((first + N - last) > excl) | (first == last) & (N > excl)
    ok &= ~((first >= 0) & (first != last) & ((first + N - last) <= excl))
    return np.sort(confs[ok])

def positions(N, eps):
    gaps = np.where(np.arange(N) % 2 == 0, 1 - eps, 1 + eps)
    x = np.concatenate([[0.0], np.cumsum(gaps)[:-1]])
    return x, gaps.sum()

def Vmat(N, eps, g):
    x, L = positions(N, eps)
    d = x[:, None] - x[None, :]
    V = np.zeros((N, N))
    for m in range(-3, 4):
        dd = np.abs(d + m * L)
        with np.errstate(divide='ignore'):
            V += np.where(dd > 1e-9, g / dd**6, 0.0)
    np.fill_diagonal(V, 0)
    return V

def bits(confs, N):
    return ((confs[:, None] >> np.arange(N)) & 1).astype(np.int8)

def hamiltonian(confs, N, eps, g, delta, omega=1.0):
    V = Vmat(N, eps, g)
    B = bits(confs, N)
    diag = -delta * B.sum(1).astype(float)
    for i in range(N):
        for j in range(i + 1, N):
            if V[i, j] > 1e-8:
                diag += V[i, j] * (B[:, i] & B[:, j])
    dim = len(confs)
    rows = [np.arange(dim)]; cols = [np.arange(dim)]; vals = [diag]
    for i in range(N):
        tgt = confs ^ (np.int64(1) << i)
        idx = np.searchsorted(confs, tgt)
        idx[idx >= dim] = 0
        ok = confs[idx] == tgt
        rows.append(np.nonzero(ok)[0]); cols.append(idx[ok]); vals.append(np.full(ok.sum(), omega / 2))
    H = csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(dim, dim))
    return H, B

def perm_index(confs, N, perm):
    # perm[i] = image of site i
    B = (confs[:, None] >> np.arange(N)) & 1
    img = (B << np.array(perm, dtype=np.int64)).sum(1)
    return np.searchsorted(confs, img)

def run(N, eps, g, delta, k=6, excl=2):
    confs = basis(N, excl)
    H, B = hamiltonian(confs, N, eps, g, delta)
    E, vec = eigsh(H, k=k, which='SA', tol=1e-10)
    o = np.argsort(E); E = E[o]; vec = vec[:, o]
    T2 = perm_index(confs, N, [(i + 2) % N for i in range(N)])
    Is = perm_index(confs, N, [(1 - i) % N for i in range(N)])
    np_ = N // 2
    sgn = np.array([(-1) ** (i // 2) for i in range(N)])  # pair sublattice sign
    lr = np.array([1 if i % 2 == 0 else -1 for i in range(N)])
    Bf = B.astype(float)
    u = Bf @ sgn / np_          # pair density wave
    tau = Bf @ lr / np_         # L/R imbalance
    v = Bf @ (sgn * lr) / np_
    out = []
    for n in range(k):
        p = vec[:, n] ** 2
        t2 = vec[:, n] @ vec[T2, n]; ip = vec[:, n] @ vec[Is, n]
        out.append((E[n] - E[0], t2, ip, p @ u**2, p @ v**2, p @ tau**2, p @ Bf.sum(1) / N))
    return len(confs), E[0] / N, out

if __name__ == '__main__':
    N = int(sys.argv[1]); g = float(sys.argv[2]) ; delta = float(sys.argv[3])
    excl = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    for eps in [float(e) for e in sys.argv[4].split(',')]:
        t = time.time()
        dim, e0, out = run(N, eps, g, delta, excl=excl)
        print(f"N={N} dim={dim} g={g:.1f} D={delta} eps={eps:.2f} e0={e0:.5f} ({time.time()-t:.1f}s)")
        for (dE, t2, ip, u2, v2, t2_, rho) in out:
            print(f"   dE={dE:8.5f} T2={t2:+.2f} Is={ip:+.2f} u2={u2:.4f} v2={v2:.4f} tau2={t2_:.4f} rho={rho:.4f}")
