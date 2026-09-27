import numpy as np
from pbc_ed import run
import sys
sys.path.insert(0,'/home/user/rydberg-atom-research')
from rydberg_chain.effective import critical_eps
D=2.5
print("TFIM leading-order eps_c(Delta=2.5) =", critical_eps(11.39, D))
for eps in [0.48,0.51,0.54,0.57]:
    line=f"eps={eps:.2f} "
    for M in (8,10,12):
        _,e,gM,rho,m2,U=run(M,eps,[D])[0]
        line+=f"| M={M} U={U:.3f} m2={m2:.3f} gapM={gM:.3f} "
    print(line,flush=True)
