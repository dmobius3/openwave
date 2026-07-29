"""Exact (tolerance-free) GF(p) rank of the averaging projector on E_m (x) V_tau."""
from part2 import rank_gf, mat_mul, mat_eye, mat_add, mat_scal, mat_trace
from part1 import p
from core import sym_power, isotypic_split
import core

def explicit_irrep(G,a,dim_target,char_target):
    """Explicit matrices for the irreducible with the given exact character, from V_a."""
    pieces,rho=isotypic_split(G,a)
    for pc in pieces:
        if pc["dim"]!=dim_target: continue
        invm=pow(pc["mult"],p-2,p)
        th=[mat_trace(pc["Rall"][G.crep[ci]],p)*invm%p for ci in range(G.nc)]
        if pc["mult"]==1 and all(th[ci]==char_target[ci] for ci in range(G.nc)):
            return pc["Rall"]
    return None

def kron(A,B,pp):
    ra,ca=len(A),len(A[0]); rb,cb=len(B),len(B[0])
    return [[A[i//rb][j//cb]*B[i%rb][j%cb]%pp for j in range(ca*cb)] for i in range(ra*rb)]

def projector_rank_exact(G,m_,taumats,dtau,dimcap=300):
    """Full unreduced projector on E_m (x) V_tau; exact GF(p) rank. None if too big."""
    dimE=2*(m_*m_-1); dim=dimE*dtau
    if dim>dimcap: return {"dim":dim,"rank":None,"reason":"above exact-rank size cap"}
    N=G.N
    S1=[sym_power(g,m_,p) for g in G.elems]
    S2=[sym_power(g,m_-2,p) for g in G.elems]
    r1,r2=m_-1,m_+1; b1=(m_+1)*r1
    P=[[0]*dim for _ in range(dim)]
    for i in range(N):
        Eg=[[0]*dimE for _ in range(dimE)]
        K1=kron(S1[i],mat_eye(r1),p)
        for x in range(b1):
            for y in range(b1): Eg[x][y]=K1[x][y]
        K2=kron(S2[i],mat_eye(r2),p)
        for x in range(dimE-b1):
            for y in range(dimE-b1): Eg[b1+x][b1+y]=K2[x][y]
        T=kron(Eg,taumats[i],p)
        for x in range(dim):
            Px=P[x]; Tx=T[x]
            for y in range(dim): Px[y]=(Px[y]+Tx[y])%p
    invN=pow(N,p-2,p)
    P=[[v*invN%p for v in row] for row in P]
    rk=rank_gf(P,dim,p)
    trc=mat_trace(P,p)
    P2=mat_mul(P,P,p)
    idem=(P2==P)
    return {"dim":dim,"rank":rk,"trace_equals_rank":(trc%p==rk%p),
            "idempotent_exactly":idem,"reason":None}
