"""Part 1: exact field setup + exact group construction over GF(p). No tolerances."""
from math import gcd
from functools import reduce

def lcm(a,b): return a*b//gcd(a,b)

def is_prime(n):
    if n < 2: return False
    for q in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % q == 0: return n == q
    d, r = n-1, 0
    while d % 2 == 0: d//=2; r+=1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x = pow(a,d,n)
        if x in (1,n-1): continue
        for _ in range(r-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

# L must be divisible by every root-of-unity order we need:
#   C_n  n=1..20  -> n            BD_n n=2..12 -> lcm(2n,4) <= 24
#   2T 12, 2O 24, 2I 60,  plus 4 (i), 8 (sqrt2), 5 (sqrt5)
import os
LMAX = int(os.environ.get('AUDIT_LMAX','24'))
L = reduce(lcm, range(1,LMAX+1), 1)
L = lcm(L, 60)
k = 1
while True:
    p = k*L + 1
    if is_prime(p): break
    k += 1
print("L =", L, " p =", p, " k =", k, " (p-1)/L =", (p-1)//L)

def primitive_root(p):
    fac = set(); n = p-1; d = 2
    while d*d <= n:
        while n % d == 0: fac.add(d); n//=d
        d += 1
    if n > 1: fac.add(n)
    for g in range(2, 1000):
        if all(pow(g,(p-1)//f,p) != 1 for f in fac): return g
    raise RuntimeError

g0 = primitive_root(p)
zL = pow(g0, (p-1)//L, p)          # primitive L-th root of unity in GF(p)
assert pow(zL, L, p) == 1 and all(pow(zL, L//q, p) != 1 for q in (2,3,5,7,11,13,17,19,23))
def zeta(M):                        # primitive M-th root of unity
    assert L % M == 0
    return pow(zL, L//M, p)
I = zeta(4)                         # sqrt(-1)
assert (I*I) % p == p-1
z8 = zeta(8); SQRT2 = (z8 + pow(z8,7,p)) % p
assert SQRT2*SQRT2 % p == 2
z5 = zeta(5); SQRT5 = (z5 + pow(z5,4,p) - pow(z5,2,p) - pow(z5,3,p)) % p
assert SQRT5*SQRT5 % p == 5
inv2 = pow(2, p-2, p)
PHI = (1 + SQRT5) * inv2 % p        # golden ratio
assert (PHI*PHI - PHI - 1) % p == 0
PHIinv = (PHI - 1) % p
assert PHI*PHIinv % p == 1
print("field ok: i,sqrt2,sqrt5,phi all exact in GF(p)")

# ---- 2x2 matrices over GF(p) as flat 4-tuples (a,b,c,d) = [[a,b],[c,d]] ----
def mmul(X,Y):
    a,b,c,d = X; e,f,g,h = Y
    return ((a*e+b*g)%p,(a*f+b*h)%p,(c*e+d*g)%p,(c*f+d*h)%p)
EYE = (1,0,0,1)
def quat(a,b,c,d):
    """a+bi+cj+dk -> [[a+bi, c+di],[-c+di, a-bi]]  (spec's map)."""
    return ((a+b*I)%p,(c+d*I)%p,(-c+d*I)%p,(a-b*I)%p)
def det(X):
    a,b,c,d = X; return (a*d-b*c)%p
def tr(X):
    a,b,c,d = X; return (a+d)%p

def close_exact(gens, cap=100000):
    """Exact closure under right multiplication: equality is dict equality, no tolerance."""
    seen = {EYE:0}; elems=[EYE]; frontier=[EYE]
    while frontier:
        nf=[]
        for M in frontier:
            for G in gens:
                P = mmul(M,G)
                if P not in seen:
                    seen[P]=len(elems); elems.append(P); nf.append(P)
        frontier=nf
        if len(elems)>cap: raise RuntimeError("cap")
    return elems, seen

def build_groups(nc_max=20, nbd_max=12):
    out=[]
    for n in range(1,nc_max+1):
        zn = zeta(n)
        out.append(("C_%d"%n, [(zn,0,0,pow(zn,n-1,p) if n>1 else 1)], n))
    for n in range(2,nbd_max+1):
        z2n = zeta(2*n)
        a = (z2n,0,0,pow(z2n,2*n-1,p))
        b = (0,1,p-1,0)
        out.append(("BD_%d"%n, [a,b], 4*n))
    g2T=[quat(inv2,inv2,inv2,inv2), quat(0,1,0,0)]
    out.append(("2T", g2T, 24))
    isq2 = pow(SQRT2,p-2,p)
    out.append(("2O", g2T+[quat(isq2,isq2,0,0)], 48))
    out.append(("2I", [quat(inv2,inv2,inv2,inv2),
                       quat(PHI*inv2%p, PHIinv*inv2%p, inv2, 0)], 120))
    return out

if __name__ == "__main__":
    bad=[]
    for name,gens,exp in build_groups():
        for G in gens:
            assert det(G)==1, (name,"det!=1")
        els,_ = close_exact(gens)
        ok = (len(els)==exp)
        if not ok: bad.append((name,len(els),exp))
        print("%-6s built=%-4d expected=%-4d %s"%(name,len(els),exp,"OK" if ok else "MISMATCH"))
    print("MISMATCHES:", bad)
