"""Scratch numbers for paper/sections/02_methods.md (27 Sept 2026)."""
import sys, json, math
sys.path.insert(0, "/home/user/rydberg-atom-research"); sys.dont_write_bytecode = True
import numpy as np
from rydberg_chain import geometry, ed, observables as ob, effective
out = {}
# 1. tail bound for pair-site cutoff mmax=8 at N=20, eps=0.3, g=11.39
n, eps, g = 20, 0.3, 11.39
x = geometry.chain_positions(n, eps); V = geometry.coupling_matrix(x, g)
pair = np.arange(n) // 2
mask = np.triu(abs(pair[:, None] - pair[None, :]) > 8, 1)
out["tail_norm_N20_mmax8"] = float(V[mask].sum())
# 2. k=0 orbit counts (necklaces) for pair rings
def neck(L, d):
    return sum(math.gcd(t, L) and d ** math.gcd(t, L) for t in range(L)) // L
out["necklace"] = {f"Np={L},d={d}": neck(L, d) for L, d in [(10, 4), (12, 3), (16, 3), (18, 3), (20, 3)]}
out["full_dims"] = {f"N={2*L}": 2 ** (2 * L) for L in (10, 12, 16, 18, 20)}
# 3. CSR memory 12 (N+1) 2^N bytes, float64 data + int32 indices (int64 if nnz >= 2^31)
mem = {}
for N in (16, 18, 20, 22, 24, 26):
    nnz = (N + 1) * 2 ** N
    ib = 4 if nnz < 2 ** 31 else 8
    mem[N] = dict(nnz=nnz, matrix_GB=nnz * (8 + ib) / 1e9, vector_GB=8 * 2 ** N / 1e9)
out["csr_memory"] = mem
# 4. couplings quoted in text
out["V_s,V_w eps0.3 g11.39"] = geometry.strong_weak_couplings(0.3, 11.39)
out["V(2a) g=1291.5"] = 1291.5 / 64
out["h0 eps0.6 g11.39 D3.149"] = effective.transverse_field(0.6, 11.39, 3.149)
c = effective.tfim_couplings(0.6, 11.39, 3.149, mmax=2); out["J1,J2 eps0.6"] = (c["J_m"][0], c["J_m"][1])
# 5. usage snippet from the manuscript
x = geometry.chain_positions(12, eps=0.3)
V = geometry.coupling_matrix(x, g=11.39)
H = ed.build_hamiltonian(V, omega=1.0, delta=2.5)
E, v, res = ed.lowest_states(H, k=2)
out["snippet"] = dict(E=E.tolist(), res=res.tolist(), occ0=float(ob.site_occupations(v[:, 0], 12)[0]),
                      S6=ob.entanglement_entropy(v[:, 0], 12, 6))
print(json.dumps(out, indent=1, default=float))
