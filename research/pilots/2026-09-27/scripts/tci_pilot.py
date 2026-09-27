import sys, time, json, numpy as np
from fast_ed import positions, ground
g = 11.39; out = []
deltas = np.round(np.arange(0.0, 3.01, 0.1), 3)
for eps in [0.0, 0.3, 0.45, 0.6, 0.7]:
    for N in [12, 14, 16]:
        t = time.time(); rows = []
        for d in deltas:
            E, n = ground(positions(N, eps), g, float(d))
            rows.append((float(d), float(n.mean()), float(E[1] - E[0])))
        dens = np.array([r[1] for r in rows]); gap = np.array([r[2] for r in rows])
        chi = np.gradient(dens, deltas)
        k = int(np.argmax(chi))
        out.append({"eps": eps, "N": N, "rows": rows})
        print(f"eps={eps} N={N}: max dn/dDelta={chi[k]:.3f} at Delta={deltas[k]:.2f}; min gap={gap.min():.4f} at Delta={deltas[int(np.argmin(gap))]:.2f}; t={time.time()-t:.0f}s", flush=True)
json.dump(out, open("tci_pilot.json", "w"))
