"""S14 item 7: gate 5's generic-pair denominators D_I = ||grad I|| * ||F|| confirmed
nonzero on the three pinned pairs at EVERY rung of both ladders.
Fields: the S 4.3 stream, generated faithfully (PCG64(20260901); agreement ladder ascending
then Control-B ladder ascending; per rung 20 scalar then 20 C^2 fields; level-major,
intertwiner, multiplet, component order; re/im per coefficient; unit L2 norm). Pairs
(0,3), (1,4), (2,5) of the scalar twenty as (psi, psi_dot). The cubic in grad E and F is
evaluated PREFLIGHT-GRADE on 500 quasi-uniform S3 sample points, adequate for a
nonzero-ness check and labeled as such; linear parts are exact in coefficients.
"""
import numpy as np, json, sys
from math import comb
rng_group = np.random.default_rng(0)

TAU=(1+np.sqrt(5))/2
def qmul(a,b):
    w1,x1,y1,z1=a; w2,x2,y2,z2=b
    return np.array([w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
                     w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2])
import itertools
def icosians():
    Q=[]
    for i in range(4):
        for s in (1.0,-1.0):
            q=np.zeros(4); q[i]=s; Q.append(q)
    for signs in range(16):
        Q.append(np.array([0.5 if (signs>>k)&1==0 else -0.5 for k in range(4)]))
    base=np.array([0.0,1.0,1/TAU,TAU])/2
    evens=[p for p in itertools.permutations(range(4))
           if sum(1 for a in range(4) for b in range(a+1,4) if p[a]>p[b])%2==0]
    for p in evens:
        for signs in range(8):
            sb=np.array([1.0,(-1)**(signs&1),(-1)**((signs>>1)&1),(-1)**((signs>>2)&1)])
            v=base*sb; Q.append(np.array([v[p.index(k)] for k in range(4)]))
    return np.array(Q)
G120=icosians()

def pi_l(l,q):
    a=complex(q[0],q[3]); b=complex(q[2],q[1])
    M=np.zeros((l+1,l+1),dtype=complex)
    for k in range(l+1):
        poly=np.zeros(l+1,dtype=complex)
        for i in range(l-k+1):
            for j in range(k+1):
                poly[i+j]+=comb(l-k,i)*a**(l-k-i)*b**i*comb(k,j)*(-np.conj(b))**(k-j)*np.conj(a)**j
        M[:,k]=poly
    return M

def chiV(l,q):
    c=np.clip(q[0],-1,1); th=np.arccos(c)
    if abs(np.sin(th))<1e-12: return (l+1)*(np.sign(c) if c!=0 else 1.0)**l
    return float(np.sin((l+1)*th)/np.sin(th))

def mult_R0(n): return round(float(np.mean([chiV(n,q) for q in G120])))
INV=[(n,mult_R0(n)) for n in range(0,61) if mult_R0(n)>0]
print(f"  R0 invariant levels <= 60: {INV}")

# intertwiner fixed vectors of pi_n(gamma)^T per level (orthonormalized if mult 2)
def intertwiners(n,m):
    P=np.zeros((n+1,n+1),dtype=complex)
    for q in G120: P+=pi_l(n,q).T
    P/=120
    U,S,Vh=np.linalg.svd(P)
    V=U[:, :m]
    # sign/phase fix
    for c in range(m):
        nz=np.argmax(np.abs(V[:,c])>1e-10)
        V[:,c]=V[:,c]*np.conj(V[nz,c])/abs(V[nz,c])
    return V
IV={n:intertwiners(n,m) for n,m in INV}
print("  intertwiners built")

# 500 quasi-uniform S3 sample points (preflight-grade cubic)
g=np.random.default_rng(777)
X=g.standard_normal((500,4)); X/=np.linalg.norm(X,axis=1,keepdims=True)
# per-level section matrices at nodes: W_{n,i}[node, j] = (iv_i^T pi_n(x))_j
W={}
for n,m in INV:
    mats=np.array([pi_l(n,x) for x in X])          # 500 x (n+1) x (n+1)
    for i in range(m):
        W[(n,i)]=np.einsum('r,xrj->xj', np.conj(IV[n][:,i]), mats)
print("  node section matrices built")

def basis_dims(N):
    return [(n,i,n+1) for n,m in INV if n<=N for i in range(m)]
def dimC(N): return sum(b[2] for b in basis_dims(N))

def spinT(n):
    s=n/2.0
    mvec=np.arange(s,-s-1,-1)
    J3=np.diag(mvec)
    ap=np.sqrt((s-mvec[1:])*(s+mvec[1:]+1))
    Jp=np.zeros((n+1,n+1)); Jp[np.arange(n),np.arange(1,n+1)]=ap
    Jm=Jp.T
    return [0.5*(Jp-Jm)+0j, 0.5j*(Jp+Jm), 1j*J3]

LADDERS=[24,32,40,48,36,44,52,60]
rng=np.random.Generator(np.random.PCG64(20260901))
results={}
for N in LADDERS:
    bd=basis_dims(N); dC=dimC(N)
    scal=[]
    for f in range(20):
        c=np.empty(dC,dtype=complex); pos=0
        for (n,i,w) in bd:
            for j in range(w):
                re=rng.standard_normal(); im=rng.standard_normal()
                c[pos]=re+1j*im; pos+=1
        c/=np.linalg.norm(c); scal.append(c)
    for f in range(20):                                   # consume the C^2 draws faithfully
        for (n,i,w) in bd:
            for j in range(w):
                for comp in range(2):
                    rng.standard_normal(); rng.standard_normal()
    lam=np.concatenate([[n*(n+2)]*w for (n,i,w) in bd]).astype(float)
    Ts=[]
    for a in range(3):
        blocks=[spinT(n)[a] for (n,i,w) in bd]
        from scipy.linalg import block_diag
        Ts.append(block_diag(*blocks))
    def field_at_nodes(c):
        v=np.zeros(500,dtype=complex); pos=0
        for (n,i,w) in bd:
            v+=W[(n,i)]@c[pos:pos+w]; pos+=w
        return v
    def cubic_proj(c):
        fv=field_at_nodes(c); gval=(np.abs(fv)**2)*fv
        out=np.empty(dimC(N),dtype=complex); pos=0
        for (n,i,w) in bd:
            out[pos:pos+w]=np.conj(W[(n,i)]).T@gval/500; pos+=w
        return out
    pairs=[(0,3),(1,4),(2,5)]
    row={}
    for (ia,ib) in pairs:
        psi, psid = scal[ia], scal[ib]
        cub=cubic_proj(psi)
        F=np.concatenate([psid, -lam*psi - cub])
        grads={
         "E": np.concatenate([lam*psi + cub, psid]),
         "Q": np.concatenate([-1j*psid, 1j*psi]),
        }
        for a in range(3):
            grads[f"M{a+1}"]=np.concatenate([-Ts[a].conj().T@psid, Ts[a]@psi])
        nF=np.linalg.norm(np.concatenate([F.real,F.imag]))
        for k,gv in grads.items():
            nG=np.linalg.norm(np.concatenate([gv.real,gv.imag]))
            row.setdefault(k,[]).append(float(nG*nF))
    results[N]={k:[f"{v:.4e}" for v in vals] for k,vals in row.items()}
    mn=min(min(float(v) for v in vals) for vals in row.values())
    print(f"  N={N:2d} dimC={dC:4d}  min D_I over 5 invariants x 3 pairs = {mn:.4e}")
allmin=min(min(float(x) for vals in results[N].values() for x in vals) for N in LADDERS)
json.dump(results,open("raw/item7_denominators.json","w"),indent=1)
print(f"  global min D_I = {allmin:.4e}  (threshold: nonzero; machine floor ~1e-13)")
ok = allmin > 1e-6
print(f"  ITEM 7 VERDICT: {'GREEN' if ok else 'RED'}")
sys.exit(0 if ok else 1)
