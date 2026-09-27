import sys, numpy as np
sys.path.insert(0, "/home/user/rydberg-atom-research")
from rydberg_chain.effective import tfim_couplings
from scipy.sparse import kron, identity, csr_matrix
from scipy.sparse.linalg import eigsh
g, delta = 11.39, 2.5
sz = csr_matrix(np.diag([1., -1.])); sx = csr_matrix(np.array([[0,1.],[1,0]]))
def op(o, k, M):
    out = identity(1, format='csr')
    for j in range(M):
        out = kron(out, o if j == k else identity(2, format='csr'), format='csr')
    return out
for M in [6, 8, 10]:
    for eps in [0.0, 0.3, 0.5, 0.6, 0.7]:
        c = tfim_couplings(eps, g, delta, mmax=M)
        H = csr_matrix((2**M, 2**M))
        Z = [op(sz, k, M) for k in range(M)]
        for mi, m in enumerate(c['m']):
            if m >= M: break
            J, B = c['J_m'][mi], c['B_m'][mi]
            for k in range(M - m):
                H = H + J * Z[k] @ Z[k+m] + B * (Z[k] - Z[k+m])
        for k in range(M):
            H = H - c['h'] * op(sx, k, M)
        E = np.sort(eigsh(H, k=6, which='SA')[0])
        print(f"pairs M={M} (N={2*M}) eps={eps}: |J1|/h={abs(c['J1'])/c['h']:.3f} TFIM gaps", np.round(E[1:]-E[0], 4))
