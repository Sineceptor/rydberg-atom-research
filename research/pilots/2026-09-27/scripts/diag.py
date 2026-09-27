import warnings, logging; logging.disable(logging.WARNING); warnings.filterwarnings("ignore")
import numpy as np, json
from dimer_model import ed_ground
N = 16
bits = (np.arange(2**N)[:, None] >> np.arange(N)) & 1
for eps in [0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.85]:
    E, v = ed_ground(N, eps, 11.39, 2.5)
    w = np.abs(v[:, 0])**2
    n = w @ bits
    top = np.argsort(w)[::-1][:3]
    print(f"eps={eps}: E0/N={E[0]/N:.5f} gap={E[1]-E[0]:.3e} dens={n.mean():.3f} centre nL-nR={n[6]-n[7]:+.3f},{n[8]-n[9]:+.3f} "
          f"top={[(''.join(map(str,bits[k])), round(float(w[k]),3)) for k in top]}")
