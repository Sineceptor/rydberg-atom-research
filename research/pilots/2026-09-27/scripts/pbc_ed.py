"""PBC ED of dimerized chain restricted to <=1 excitation per strong pair (3 states/pair).
Digits: 0 empty, 1 L excited, 2 R excited. Minimal-image distances on ring of length 2M."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla, sys, time
def build(M, eps, g, omega=1.0):
    s=1-eps; L=2.0*M
    dim=3**M
    idx=np.arange(dim)
    digs=np.empty((dim,M),dtype=np.int8)
    t=idx.copy()
    for k in range(M):
        digs[:,k]=t%3; t//=3
    occ=(digs>0).astype(np.int8)
    tau=np.where(digs==1,1,np.where(digs==2,-1,0)).astype(np.int8)
    pos=lambda k,d: 2.0*k+(0.0 if d==1 else s)
    Eint=np.zeros(dim)
    for k in range(M):
        for l in range(k+1,M):
            T=np.zeros((3,3))
            for a in (1,2):
                for b in (1,2):
                    r=abs(pos(k,a)-pos(l,b)); r=min(r,L-r)
                    T[a,b]=g/r**6
            Eint+=T[digs[:,k],digs[:,l]]
    rows=[];cols=[]
    for k in range(M):
        e=(digs[:,k]==0)
        src=idx[e]
        for d in (1,2):
            dst=src+d*3**k
            rows+= [src,dst]; cols+=[dst,src]
    rows=np.concatenate(rows); cols=np.concatenate(cols)
    Hoff=sp.csr_matrix((np.full(len(rows),omega/2),(rows,cols)),shape=(dim,dim))
    return Hoff, Eint, occ.sum(1), tau.sum(1)
def run(M, eps, deltas, g=11.39, k=2):
    Hoff,Eint,N,T=build(M,eps,g)
    out=[]
    v0=None
    for D in deltas:
        H=Hoff+sp.diags(Eint-D*N)
        w,v=sla.eigsh(H,k=k,which='SA',v0=v0,tol=1e-10)
        o=np.argsort(w); w=w[o]; v=v[:,o]
        p=np.abs(v[:,0])**2; v0=v[:,0]
        m=T/M
        m2=p@m**2; m4=p@m**4
        rho=p@N/(2*M)
        out.append((D,w[0]/M,(w[1]-w[0])*M,rho,m2,1-m4/(3*m2**2)))
    return out
if __name__=="__main__":
    eps=float(sys.argv[1]); Ms=[int(x) for x in sys.argv[2].split(',')]
    deltas=np.arange(float(sys.argv[3]),float(sys.argv[4])+1e-9,float(sys.argv[5]))
    res={}
    for M in Ms:
        t0=time.time(); res[M]=run(M,eps,deltas); print(f"# M={M} time {time.time()-t0:.1f}s",flush=True)
    print(f"# eps={eps}  columns: Delta | per M: rho, m2, Binder U, gap*M")
    for i,D in enumerate(deltas):
        line=f"{D:5.2f} "
        for M in Ms:
            _,e,gM,rho,m2,U=res[M][i]
            line+=f"| {rho:.3f} {m2:.3f} {U:6.3f} {gM:6.3f} "
        print(line)
