import numpy as np
from scipy.optimize import brentq
g=1291.5
V=lambda r: g/r**6
def J2(eps):
    s=1-eps; return 0.25*(V(4-s)+V(4+s)-2*V(4))
def hd(D,eps):
    s=1-eps; IL=2*V(4); IR=V(4-s)+V(4+s); Vs=V(s)
    return 0.125*(1/(D-IL)+1/(D-IR))+0.25/(Vs-D)
def h0(D,eps):
    s=1-eps; Vs=V(s); return 0.25*(1/D+1/(Vs-D))
for eps in [0.2,0.3,0.4,0.5,0.7]:
    print(eps, "J2=%.4f hd(3)=%.4f h0(3)=%.4f J2/hd=%.3f"%(J2(eps),hd(3,eps),h0(3,eps),J2(eps)/hd(3,eps)))
print("Delta=3: dressed eps_tau =", brentq(lambda e: J2(e)-hd(3,e), 0.05, 0.9), " leading eps_tau =", brentq(lambda e: J2(e)-h0(3,e), 0.05, 0.9))
for D in [1.5,2,3,5,6.5]:
    print("eps=0.5 D=",D," J2/hd=%.3f J2/h0=%.3f"%(J2(0.5)/hd(D,0.5), J2(0.5)/h0(D,0.5)))
s=0.5; IR=V(4-s)+V(4+s)
print("eps=0.5: dressed Delta_tau =", brentq(lambda D: J2(0.5)-hd(D,0.5), IR+0.01, 20), "leading:", brentq(lambda D: J2(0.5)-h0(D,0.5), 0.5, 20), "IL=",2*V(4),"IR=",IR)
