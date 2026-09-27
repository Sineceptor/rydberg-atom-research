"""Scratch check for the Methods draft: DMRG vs ED energy differences (same cases as tests/test_dmrg.py)."""
import json, logging, sys, time
import numpy as np
sys.path.insert(0, "/home/user/rydberg-atom-research")
sys.dont_write_bytecode = True
logging.disable(logging.WARNING)
from rydberg_chain import ed, geometry, observables as ob, dmrg
G, D = 11.39, 2.5
out = {"atom": [], "pair": [], "saved": None}
for n, eps in [(10, 0.0), (11, 0.3), (12, 0.3)]:
    x = geometry.chain_positions(n, eps)
    E, v, r = ed.lowest_states(ed.build_hamiltonian(geometry.coupling_matrix(x, G), 1.0, D), k=2)
    e, psi, obs = dmrg.atom_chain_ground_state(x, G, D, chi=64)
    out["atom"].append(dict(n=n, eps=eps, dE=e - E[0], occ=float(np.max(abs(obs["occupations"] - ob.site_occupations(v[:, 0], n)))), trunc=obs["final_trunc_err"]))
for P, eps, tight in [(6, 0.0, False), (6, 0.3, True), (7, 0.6, False)]:
    n = 2 * P
    x = geometry.chain_positions(n, eps)
    E, v, r = ed.lowest_states(ed.build_hamiltonian(geometry.coupling_matrix(x, G), 1.0, D), k=2)
    kw = {"max_S_err": 1e-10} if tight else {}
    e, psi, obs = dmrg.pair_chain_ground_state(eps, G, D, P, chi=64, **kw)
    out["pair"].append(dict(n=n, eps=eps, tight=tight, dE=e - E[0], occ=float(np.max(abs(obs["occupations"] - ob.site_occupations(v[:, 0], n)))), sweeps=obs["n_sweeps"]))
e, psi, obs = dmrg.atom_chain_ground_state(geometry.chain_positions(12), 8.0, 3.0, chi=64, variance=True)
E, v, r = ed.lowest_states(ed.chain_hamiltonian(12, 0.0, 8.0, 3.0), k=2)
out["saved"] = dict(dmrg=e, ed=float(E[0]), saved=-18.399305357538882, dmrg_minus_saved=e + 18.399305357538882, ed_minus_saved=float(E[0]) + 18.399305357538882, variance=obs["variance"], gap=float(E[1]-E[0]))
print(json.dumps(out, indent=1, default=float))
