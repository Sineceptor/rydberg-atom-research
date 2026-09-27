import sys, time
import numpy as np
sys.path.insert(0, "/home/user/rydberg-atom-research")
from rydberg_chain.geometry import chain_positions, coupling_matrix
from rydberg_chain.ed import build_hamiltonian, lowest_states, diagonal_energies
from rydberg_chain.effective import transverse_field

g, delta = 11.39, 2.5
def wall_configs(N):
    # MIS configs: (10)^k (01)^(N/2-k), k=0..N/2 ; bit i = atom i
    out = []
    for k in range(N//2 + 1):
        occ = [1,0]*k + [0,1]*(N//2 - k)
        out.append(sum(b << i for i, b in enumerate(occ)))
    return out

def eff_model(N, eps, omega=1.0):
    """Second-order SW (symmetrized denominators) projected on the one-wall MIS manifold."""
    V = coupling_matrix(chain_positions(N, eps), g)
    Ecl = lambda s: -delta*bin(s).count('1') + sum(V[i,j] for i in range(N) for j in range(i+1,N) if (s>>i)&1 and (s>>j)&1)
    P = wall_configs(N); idx = {s:a for a,s in enumerate(P)}
    E0 = np.array([Ecl(s) for s in P]); M = len(P)
    Heff = np.diag(E0).astype(float)
    # intermediate single flips
    inter = {}
    for a, s in enumerate(P):
        for i in range(N):
            m = s ^ (1 << i)
            if m in idx: continue
            inter.setdefault(m, []).append(a)
    for m, lst in inter.items():
        Em = Ecl(m)
        for a in lst:
            for b in lst:
                Heff[a,b] += (omega/2)**2 * 0.5*(1/(E0[a]-Em) + 1/(E0[b]-Em))
    return Heff, P

for N in [12, 16, 20]:
    for eps in [0.0, 0.3, 0.5, 0.6]:
        t0 = time.time()
        V = coupling_matrix(chain_positions(N, eps), g)
        H = build_hamiltonian(V, 1.0, delta)
        k = 8
        E, vec, res = lowest_states(H, k=k, method="lobpcg" if N >= 16 else "auto")
        P = wall_configs(N)
        w = np.array([np.sum(vec[P, n]**2) for n in range(k)])
        pk = vec[P, 0]**2
        Heff, _ = eff_model(N, eps)
        Ee = np.linalg.eigvalsh(Heff)
        h = transverse_field(eps, g, delta)
        M = N//2 + 1
        print(f"N={N} eps={eps} time={time.time()-t0:.1f}s h={h:.4f}")
        print("  ED  gaps  :", np.round(E[1:k]-E[0], 4))
        print("  eff gaps  :", np.round(Ee[1:k]-Ee[0], 4))
        print("  ED E0, eff E0:", round(E[0],5), round(Ee[0],5))
        print("  weight on MIS wall manifold per state:", np.round(w, 3))
        print("  GS wall-position prob (k=0..N/2):", np.round(pk/pk.sum(), 3))
        print("  Heff diag (shifted):", np.round(np.diag(Heff)-Heff[1,1], 4))
        print("  Heff hops:", np.round(np.diag(Heff, 1), 4))
        sys.stdout.flush()
