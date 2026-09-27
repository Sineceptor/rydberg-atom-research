"""Spot check: Z2 critical point from (8,10)-pair ring crossing, full pair basis d=4 vs pilot d=3 (no |11>), eps=0.3, g=11.39."""
import sys, json, time
from scipy.optimize import brentq
sys.path.insert(0, "/home/user/rydberg-atom-research/research/pilots/2026-09-27/scripts"); sys.dont_write_bytecode = True
import pairk0
eps, g = 0.3, 11.39
res = {}
for d in (4, 3):
    t = time.time()
    B = {Np: pairk0.build(Np, d, eps, g) for Np in (8, 10)}
    cache = {}
    def gap(Np, D):
        k = (Np, round(D, 10))
        if k not in cache: cache[k] = pairk0.gaps(B[Np], D, k=8)
        return cache[k]
    f = lambda D: 10 * gap(10, D)[0] - 8 * gap(8, D)[0]
    Dc = brentq(f, 1.25, 1.40, xtol=1e-6)
    g10 = gap(10, Dc)
    res[f"d={d}"] = dict(Dc=Dc, R10=g10[1] / g10[0], Lgap10=10 * g10[0], dim10=int(B[10]["nrep"]), seconds=time.time() - t)
    print(d, res[f"d={d}"], flush=True)
res["Dc_d4_minus_d3"] = res["d=4"]["Dc"] - res["d=3"]["Dc"]
print(json.dumps(res, indent=1))
