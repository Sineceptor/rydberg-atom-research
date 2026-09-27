import sys, numpy as np, time
sys.path.insert(0,'/home/user/rydberg-atom-research')
from rydberg_chain.ed import build_hamiltonian
from scipy.sparse.linalg import eigsh
import ced
N=16; g=1291.5; D=3.0
for eps in [0.0,0.3,0.5,0.7]:
    V=ced.Vmat(N,eps,g)
    H=build_hamiltonian(V, omega=1.0, delta=D)
    E,vec=eigsh(H.tocsr() if hasattr(H,'tocsr') else H,k=6,which='SA',tol=1e-11)
    o=np.argsort(E);E=E[o];vec=vec[:,o]
    confs=np.arange(2**N,dtype=np.int64)
    # bit convention check: assume bit i of index = atom i (may be reversed; symmetric ops handle both for T2 up to direction)
    T2=ced.perm_index(confs,N,[(i+2)%N for i in range(N)])
    Is=ced.perm_index(confs,N,[(1-i)%N for i in range(N)])
    dimc, e0c, outc = ced.run(N,eps,g,D)
    print(f"eps={eps} full E0/N={E[0]/N:.6f} constrained E0/N={e0c:.6f} diff/N={E[0]/N-e0c:.2e}")
    for n in range(6):
        print(f"   full dE={E[n]-E[0]:.5f} T2={vec[:,n]@vec[T2,n]:+.2f} Is={vec[:,n]@vec[Is,n]:+.2f} | constrained dE={outc[n][0]:.5f} T2={outc[n][1]:+.2f} Is={outc[n][2]:+.2f}")
