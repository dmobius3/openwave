"""Part 2: exact linear algebra over GF(p), polynomial root finding, cyclotomic ring Z[zeta_M]."""
from fractions import Fraction
from math import gcd
import random

# ---------------- linear algebra over GF(p) (exact) ----------------
def rref(rows, ncol, p):
    """In-place row reduce; returns (rows, pivot_cols)."""
    R=[r[:] for r in rows]; piv=[]; r=0
    for c in range(ncol):
        pr=None
        for i in range(r,len(R)):
            if R[i][c]%p: pr=i; break
        if pr is None: continue
        R[r],R[pr]=R[pr],R[r]
        inv=pow(R[r][c],p-2,p)
        R[r]=[(v*inv)%p for v in R[r]]
        for i in range(len(R)):
            if i!=r and R[i][c]%p:
                f=R[i][c]
                R[i]=[(a-f*b)%p for a,b in zip(R[i],R[r])]
        piv.append(c); r+=1
        if r==len(R): break
    return R[:r],piv

def nullspace(rows,ncol,p):
    """Basis of {x : rows.x = 0}, as list of length-ncol vectors."""
    R,piv=rref(rows,ncol,p)
    free=[c for c in range(ncol) if c not in piv]
    basis=[]
    for fc in free:
        v=[0]*ncol; v[fc]=1
        for i,pc in enumerate(piv):
            v[pc]=(-R[i][fc])%p
        basis.append(v)
    return basis

def rank_gf(rows,ncol,p):
    return len(rref(rows,ncol,p)[0])

def mat_mul(A,B,p):
    n=len(A); k=len(B); m=len(B[0])
    Bt=list(zip(*B))
    return [[sum(A[i][t]*Bt[j][t] for t in range(k))%p for j in range(m)] for i in range(n)]

def mat_eye(n): return [[1 if i==j else 0 for j in range(n)] for i in range(n)]
def mat_add(A,B,p): return [[(a+b)%p for a,b in zip(ra,rb)] for ra,rb in zip(A,B)]
def mat_scal(A,c,p): return [[(c*a)%p for a in ra] for ra in A]
def mat_trace(A,p): return sum(A[i][i] for i in range(len(A)))%p

# ---------------- polynomials over GF(p) (coeff list, low->high) ----------------
def pnorm(f,p):
    f=[c%p for c in f]
    while len(f)>1 and f[-1]==0: f.pop()
    return f
def padd(f,g,p):
    n=max(len(f),len(g)); return pnorm([ (f[i] if i<len(f) else 0)+(g[i] if i<len(g) else 0) for i in range(n)],p)
def psub(f,g,p):
    n=max(len(f),len(g)); return pnorm([ (f[i] if i<len(f) else 0)-(g[i] if i<len(g) else 0) for i in range(n)],p)
def pmul(f,g,p):
    if f==[0] or g==[0]: return [0]
    out=[0]*(len(f)+len(g)-1)
    for i,a in enumerate(f):
        if a:
            for j,b in enumerate(g): out[i+j]=(out[i+j]+a*b)%p
    return pnorm(out,p)
def pdivmod(f,g,p):
    f=pnorm(f,p); g=pnorm(g,p)
    if len(f)<len(g): return [0],f
    inv=pow(g[-1],p-2,p); q=[0]*(len(f)-len(g)+1); r=f[:]
    for i in range(len(f)-len(g),-1,-1):
        c=r[i+len(g)-1]*inv%p
        q[i]=c
        if c:
            for j in range(len(g)): r[i+j]=(r[i+j]-c*g[j])%p
    return pnorm(q,p),pnorm(r,p)
def pgcd(f,g,p):
    f=pnorm(f,p); g=pnorm(g,p)
    while g!=[0]:
        f,g=g,pdivmod(f,g,p)[1]
    if f!=[0]:
        inv=pow(f[-1],p-2,p); f=[c*inv%p for c in f]
    return f
def ppowmod(base,e,mod,p):
    r=[1]; b=pdivmod(base,mod,p)[1]
    while e:
        if e&1: r=pdivmod(pmul(r,b,p),mod,p)[1]
        e>>=1; b=pdivmod(pmul(b,b,p),mod,p)[1]
    return r

def poly_roots_split(f,p,rng):
    """All roots of f, assuming f is squarefree and splits completely over GF(p)."""
    f=pnorm(f,p)
    # strip zero roots
    roots=[]
    while len(f)>1 and f[0]==0:
        roots.append(0); f=pnorm(f[1:],p)
    def rec(h):
        h=pnorm(h,p)
        if len(h)<=1: return []
        if len(h)==2:            # c1 x + c0
            return [(-h[0]*pow(h[1],p-2,p))%p]
        while True:
            r=rng.randrange(p)
            t=ppowmod([r,1],(p-1)//2,h,p)
            g=pgcd(psub(t,[1],p),h,p)
            if 0<len(g)-1<len(h)-1:
                q=pdivmod(h,g,p)[0]
                return rec(g)+rec(q)
    return roots+rec(f)

def minpoly_of_matrix(M,p):
    """Minimal polynomial of square matrix M over GF(p)."""
    n=len(M); vecs=[]; P=mat_eye(n)
    powers=[P]
    for k in range(1,n+1):
        P=mat_mul(P,M,p); powers.append(P)
    # find least k with I,M,...,M^k dependent
    rows=[]
    for k in range(n+1):
        rows.append([powers[k][i][j] for i in range(n) for j in range(n)])
        R,piv=rref(rows,n*n,p)
        if len(R)<len(rows):
            ns=nullspace([[rows[t][c] for t in range(len(rows))] for c in range(n*n)],len(rows),p)
            assert ns
            return pnorm(ns[0],p)
    raise RuntimeError("minpoly failed")

# ---------------- exact cyclotomic ring Z[zeta_M] ----------------
_CYC={}
def cyclotomic_poly(n):
    """Phi_n as integer coeff list, low->high (exact, computed by division)."""
    if n in _CYC: return _CYC[n]
    num=[-1]+[0]*(n-1)+[1]
    for d in range(1,n):
        if n%d==0:
            num=_intpolydiv(num,cyclotomic_poly(d))
    _CYC[n]=num; return num
def _intpolydiv(num,den):
    num=num[:]; dn=len(den)-1; q=[0]*max(0,len(num)-dn)
    for i in range(len(num)-1,dn-1,-1):
        c=num[i]
        if c:
            q[i-dn]=c
            for j in range(dn+1): num[i-dn+j]-=c*den[j]
    assert all(c==0 for c in num[:dn]), "non-exact division"
    return q

class Cyc:
    """Integer/rational combination of zeta_M^k, stored as coeff list length M (x^M-1 model)."""
    __slots__=("M","c")
    def __init__(self,M,c=None):
        self.M=M; self.c=[0]*M if c is None else list(c)
    @staticmethod
    def mono(M,k,coef=1):
        z=Cyc(M); z.c[k%M]+=coef; return z
    def __add__(s,o): return Cyc(s.M,[a+b for a,b in zip(s.c,o.c)])
    def __sub__(s,o): return Cyc(s.M,[a-b for a,b in zip(s.c,o.c)])
    def __mul__(s,o):
        M=s.M; out=[0]*M
        for i,a in enumerate(s.c):
            if a:
                for j,b in enumerate(o.c):
                    if b: out[(i+j)%M]+=a*b
        return Cyc(M,out)
    def scal(s,n): return Cyc(s.M,[n*a for a in s.c])
    def conj(s):
        M=s.M; out=[0]*M
        for i,a in enumerate(s.c): out[(-i)%M]+=a
        return Cyc(M,out)
    def is_zero(s):
        return s.to_rational()==0 if s.rational() else all(v==0 for v in s.reduced())
    def reduced(s):
        r=[Fraction(v) for v in s.c]
        while len(r)>1 and r[-1]==0: r.pop()
        ph=[Fraction(v) for v in cyclotomic_poly(s.M)]
        dn=len(ph)-1
        for i in range(len(r)-1,dn-1,-1):
            c=r[i]
            if c:
                for j in range(dn+1): r[i-dn+j]-=c*ph[j]
        r=r[:dn] if len(r)>=dn else r
        while len(r)>1 and r[-1]==0: r.pop()
        return r
    def rational(s):
        return len(s.reduced())<=1
    def to_rational(s):
        rd=s.reduced()
        assert len(rd)<=1, "not rational: %r"%(rd,)
        return rd[0] if rd else Fraction(0)

if __name__=="__main__":
    p=10708457761
    rng=random.Random(7)
    # self-tests
    assert cyclotomic_poly(1)==[-1,1]
    assert cyclotomic_poly(4)==[1,0,1]
    assert cyclotomic_poly(6)==[1,-1,1]
    assert cyclotomic_poly(12)==[1,0,-1,0,1]
    assert cyclotomic_poly(8)==[1,0,0,0,1]
    z=Cyc.mono(5,1)
    s=Cyc(5)
    for k in range(5): s=s+Cyc.mono(5,k)
    assert s.to_rational()==0, s.reduced()
    t=(Cyc.mono(3,1)+Cyc.mono(3,2)).to_rational(); assert t==-1, t
    M=[[2,1,0],[0,2,0],[0,0,5]]
    mp_=minpoly_of_matrix(M,p)
    print("minpoly test:",mp_)
    print("roots:",sorted(poly_roots_split(pgcd(mp_,psub(ppowmod([0,1],p,mp_,p),[0,1],p),p),p,rng)))
    print("part2 self-tests OK")
