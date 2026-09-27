import numpy as np
g=11.39
V=lambda r: g/r**6
print("eps   s     w    V_w    J1     c_bdry   A_occ   h(D=2.5) h(D=1.5) Ren(D=2.5)")
for eps in [0.2,0.3,0.4,0.5,0.6,0.7,0.8]:
    s=1-eps; w=1+eps
    J=(V(w)+V(2+s)-2*V(2))/4
    c=(V(w)-V(2+s))/4
    A=(2*V(2)+V(2+s)+V(w))/4
    Vs=V(s)
    h=lambda D: 0.25*(1/D+1/(Vs-D))
    print(f"{eps:.2f} {s:.2f} {w:.2f} {V(w):6.3f} {J:6.3f} {c:7.3f} {A:7.3f} {h(2.5):7.3f} {h(1.5):7.3f}  {1-1/(2*2.5**2):.3f}")
# hardware: min spacing 4um, R_b from C6 (n=70 ~ 862690 MHz um^6) and Omega/2pi
C6=862690.
for Om in [1.0,1.5,2.5]:
    Rb=(C6/Om)**(1/6)
    for rb_a in [1.2,1.35,1.5]:
        a=Rb/rb_a; epsmax=1-4.0/a
        print(f"Omega/2pi={Om} MHz Rb={Rb:.2f}um Rb/a={rb_a} a={a:.2f}um eps_max={epsmax:.2f}")
