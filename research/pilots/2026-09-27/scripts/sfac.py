import sys, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
from pairk0 import build
def reps_digits(Np, d):
    D = d**Np; c = np.arange(D, dtype=np.int64); pw = d**(Np-1)
    rep = c.copy(); rot = c.copy()
    for t in range(1, Np):
        rot = rot // d + (rot % d) * pw; np.minimum(rep, rot, out=rep)
    reps = np.nonzero(rep == c)[0]
    dig = np.empty((len(reps), Np), dtype=np.int8); tmp = reps.copy()
    for j in range(Np): dig[:, j] = tmp % d; tmp //= d
    return dig
Np, d, eps, g = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
B = build(Np, d, eps, g); dig = reps_digits(Np, d)
Nocc = ((dig == 1) | (dig == 2)).astype(float) + 2*(dig == 3)
tau = (dig == 1).astype(float) - (dig == 2)
FN = np.abs(np.fft.fft(Nocc, axis=1))**2 / Np; FT = np.abs(np.fft.fft(tau, axis=1))**2 / Np
qs = 2*np.pi*np.arange(Np)/Np
for delta in [float(v) for v in sys.argv[5:]]:
    H = B['Hoff'] + sp.diags(B['Eint'] - delta * B['nsum'])
    E, V = sla.eigsh(H, k=1, which='SA', tol=1e-10, v0=np.ones(B['nrep']))
    p = V[:, 0]**2
    SN = p @ FN; ST = p @ FT; dens = (p @ Nocc.sum(1)) / (2*Np)
    i3 = Np // 3 if Np % 3 == 0 else None
    print(f'Np={Np} eps={eps} g={g} Delta={delta:5.2f} n={dens:.3f}  SN(pi)/Np={SN[Np//2]/Np:.3f}  '
          + (f'SN(2pi/3)/Np={SN[i3]/Np:.3f}  ' if i3 else '') + f'Stau(0)/Np={ST[0]/Np:.3f} Stau(pi)/Np={ST[Np//2]/Np:.3f}', flush=True)
