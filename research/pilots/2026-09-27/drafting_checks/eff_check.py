import numpy as np
from scipy.optimize import brentq
g=11.39
V=lambda r: g/r**6
ED={0.3:1.3257545685911245,0.45:1.8052382018518263,0.6:3.149354252011604,0.7:5.610715112526968}
for eps,dc in ED.items():
    s=1-eps; w=1+eps
    J1=0.25*(V(w)+V(2+s)-2*V(2))
    Vs=V(s)
    h0=lambda D: 0.25*(1/D+1/(Vs-D))
    IL=2*V(2); IR=V(w)+V(2+s)
    hd=lambda D: 0.125*(1/(D-IL)+1/(D-IR))+0.25/(Vs-D)
    D0=brentq(lambda D: h0(D)-J1, 1e-3, Vs/2)
    lo=max(IL,IR)+1e-6
    try: Dd=brentq(lambda D: hd(D)-J1, lo, 50)
    except Exception as e: Dd=np.nan
    print(f"eps={eps} J1={J1:.4f} Vw/4={V(w)/4:.4f} red={1-J1/(V(w)/4):.3f} IL={IL:.3f} IR={IR:.3f} Vs={Vs:.1f} ED={dc:.4f} LO={D0:.4f} ({(D0-dc)/dc*100:+.1f}%) dressed={Dd:.4f} ({(Dd-dc)/dc*100:+.1f}%)")
# J_m values, m=1,2 at several eps for g=11.39 and g=1291.5
for gg in [11.39,1291.5]:
    for eps in [0,0.3,0.5,0.7,0.9]:
        s=1-eps
        Jm=[0.25*(gg/(2*m-s)**6+gg/(2*m+s)**6-2*gg/(2*m)**6) for m in (1,2)]
        print(gg,eps,["%.4g"%j for j in Jm], "V(2a)=%.4g"%(gg/64))
