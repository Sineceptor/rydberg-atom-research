"""k=0, inversion-resolved ED of the dimerized Rydberg ring in the pair basis.
Pair states: 0=00, 1=10 (L), 2=01 (R), 3=11 (only if d=4)."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, sys, time

def build(Np, d, eps, g):
    s = 1.0 - eps
    Na = 2*Np; Lring = 2.0*Np
    x = np.empty(Na); x[0::2] = 2*np.arange(Np); x[1::2] = 2*np.arange(Np) + s
    dx = np.abs(x[:, None] - x[None, :]); dx = np.minimum(dx, Lring - dx)
    V = np.zeros((Na, Na)); off = ~np.eye(Na, dtype=bool); V[off] = g / dx[off]**6
    D = d**Np
    c = np.arange(D, dtype=np.int64)
    pw = d**(Np-1)
    rep = c.copy(); rot = c.copy()
    for t in range(1, Np):
        rot = rot // d + (rot % d) * pw
        np.minimum(rep, rot, out=rep)
    reps = np.nonzero(rep == c)[0]
    idx = np.full(D, -1, dtype=np.int64); idx[reps] = np.arange(len(reps))
    ridx = idx[rep]  # index of orbit for every config
    del rep, rot
    # periods
    per = np.zeros(len(reps), dtype=np.int64); rot = reps.copy(); done = np.zeros(len(reps), bool)
    for t in range(1, Np+1):
        rot = rot // d + (rot % d) * pw
        hit = (rot == reps) & ~done
        per[hit] = t; done |= hit
    # digits of reps
    dig = np.empty((len(reps), Np), dtype=np.int8); tmp = reps.copy()
    for j in range(Np):
        dig[:, j] = tmp % d; tmp //= d
    occ = np.empty((len(reps), Na), dtype=np.float32)
    occ[:, 0::2] = ((dig == 1) | (dig == 3)); occ[:, 1::2] = ((dig == 2) | (dig == 3))
    nsum = occ.sum(1).astype(np.float64)
    Eint = 0.5 * np.einsum('ij,ij->i', occ @ V.astype(np.float32), occ).astype(np.float64)
    # off-diagonal
    rows, cols, vals = [], [], []
    for j in range(Np):
        pwj = d**j
        sj = dig[:, j].astype(np.int64)
        for atom in (0, 1):
            if d == 4:
                if atom == 0: new = np.array([1, 0, 3, 2])[sj]
                else: new = np.array([2, 3, 0, 1])[sj]
                ok = np.ones(len(sj), bool)
            else:
                if atom == 0: tab = np.array([1, 0, -1]); 
                else: tab = np.array([2, -1, 0])
                new = tab[sj]; ok = new >= 0
            cnew = reps[ok] + (new[ok] - sj[ok]) * pwj
            tgt = ridx[cnew]; src = np.nonzero(ok)[0]
            rows.append(tgt); cols.append(src)
            vals.append(0.5 * np.sqrt(per[src] / per[tgt]))
    rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
    Hoff = sp.csr_matrix((vals, (rows, cols)), shape=(len(reps), len(reps)))
    # inversion: pair j -> pair (-j mod Np), swap L<->R
    sw = np.array([0, 2, 1, 3])[:d]
    pdig = sw[dig[:, (-np.arange(Np)) % Np]]
    pc = (pdig.astype(np.int64) * (d**np.arange(Np))).sum(1)
    pidx = ridx[pc]
    return dict(Hoff=Hoff, nsum=nsum, Eint=Eint, pidx=pidx, nrep=len(reps), Np=Np, d=d)

def spectrum(B, delta, omega=1.0, k=8):
    H = omega * B['Hoff'] + sp.diags(B['Eint'] - delta * B['nsum'])
    E, Vv = sla.eigsh(H, k=k, which='SA', tol=1e-10, v0=np.ones(B['nrep']))
    o = np.argsort(E); E = E[o]; Vv = Vv[:, o]
    PV = Vv[B['pidx'], :]
    Pm = Vv.T @ PV
    M = np.diag(E) + 1e-7 * Pm
    w, U = np.linalg.eigh(M)
    W = Vv @ U
    En = np.einsum('ij,ij->j', W, H @ W); par = np.einsum('ij,ij->j', W, W[B['pidx'], :])
    o = np.argsort(En)
    return En[o], par[o]

def gaps(B, delta, k=8):
    E, par = spectrum(B, delta, k=k)
    E0 = E[0]
    ev = [e for e, p in zip(E, par) if p > 0.5]; od = [e for e, p in zip(E, par) if p < -0.5]
    gs = od[0] - E0 if od else np.nan
    ge = ev[1] - E0 if len(ev) > 1 else np.nan
    return gs, ge, par[0]

if __name__ == '__main__':
    Np, d, eps, g = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    t = time.time(); B = build(Np, d, eps, g); print('build', B['nrep'], time.time() - t)
    for delta in [float(v) for v in sys.argv[5:]]:
        t = time.time(); gs, ge, p0 = gaps(B, delta)
        print(f'delta={delta:.3f} gap_sigma={gs:.5f} gap_eps={ge:.5f} ratio={ge/gs:.3f} P0={p0:+.2f} t={time.time()-t:.1f}')
