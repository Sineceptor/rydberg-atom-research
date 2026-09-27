import numpy as np, scipy.sparse as sp
from pairk0 import build, spectrum
import itertools
Np, d, eps, g, delta = 3, 4, 0.0, 11.39, 1.0
B = build(Np, d, eps, g)
H = B['Hoff']; print('herm', abs(H - H.T).max())
# full basis brute force
s = 1-eps; Na=2*Np; x=np.empty(Na); x[0::2]=2*np.arange(Np); x[1::2]=2*np.arange(Np)+s
dx=np.abs(x[:,None]-x[None,:]); dx=np.minimum(dx,2*Np-dx); V=np.zeros((Na,Na)); off=~np.eye(Na,dtype=bool); V[off]=g/dx[off]**6
confs=[]
for c in itertools.product([0,1],repeat=Na):
    c=np.array(c)
    if d==3 and any(c[2*j]&c[2*j+1] for j in range(Np)): continue
    confs.append(c)
confs=np.array(confs); ix={tuple(c):i for i,c in enumerate(confs)}
n=len(confs); Hf=np.zeros((n,n))
for i,c in enumerate(confs):
    Hf[i,i]=-delta*c.sum()+0.5*c@V@c
    for a in range(Na):
        c2=c.copy(); c2[a]^=1
        if tuple(c2) in ix: Hf[ix[tuple(c2)],i]+=0.5
Ef,Uf=np.linalg.eigh(Hf)
# translation by one pair: atom a -> a+2
T=np.zeros((n,n))
for i,c in enumerate(confs): T[ix[tuple(np.roll(c,2))],i]=1
kvals=[np.real(u@T@u) for u in Uf.T[:12]]
print(np.round(Ef[:12],5)); print(np.round(kvals,3))
Hd = (B['Hoff'] + sp.diags(B['Eint'] - delta*B['nsum'])).toarray()
print(np.round(np.linalg.eigvalsh(Hd)[:8],5))
