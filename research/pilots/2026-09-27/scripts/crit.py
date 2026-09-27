import sys, numpy as np, time, json
from scipy.optimize import brentq
from pairk0 import build, gaps
eps, g, d = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
sizes = [int(v) for v in sys.argv[4].split(',')]
lo, hi = float(sys.argv[5]), float(sys.argv[6])
Bs = {Np: build(Np, d, eps, g) for Np in sizes}
cache = {}
def G(Np, D):
    key = (Np, round(D, 10))
    if key not in cache: cache[key] = gaps(Bs[Np], D, k=8)
    return cache[key]
out = {'eps': eps, 'g': g, 'd': d, 'pairs': []}
for a, b in zip(sizes[:-1], sizes[1:]):
    f = lambda D: b * G(b, D)[0] - a * G(a, D)[0]
    try:
        Dc = brentq(f, lo, hi, xtol=1e-4)
    except ValueError as e:
        print('no crossing', a, b, f(lo), f(hi)); continue
    ga, gb = G(a, Dc), G(b, Dc)
    rec = dict(sizes=(a, b), Dc=Dc, R_a=ga[1]/ga[0], R_b=gb[1]/gb[0], Lgap_b=b*gb[0])
    out['pairs'].append(rec); print(json.dumps(rec), flush=True)
json.dump(out, open(f'crit_eps{eps}_g{g}_d{d}.json', 'w'))
