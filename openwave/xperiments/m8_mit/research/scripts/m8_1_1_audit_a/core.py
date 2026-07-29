"""Core of the independent audit: exact GF(p) group/rep theory + exact Z[zeta_M] characters."""
from fractions import Fraction
from math import gcd
from functools import reduce
import random, sys
from part2 import (rref, nullspace, rank_gf, mat_mul, mat_eye, mat_add, mat_scal,
                   mat_trace, pnorm, padd, psub, pmul, pdivmod, pgcd, ppowmod,
                   poly_roots_split, minpoly_of_matrix, cyclotomic_poly, Cyc)
from part1 import (p, L, zeta, I, SQRT2, SQRT5, PHI, PHIinv, inv2,
                   mmul, EYE, quat, det, tr, build_groups)

RNG = random.Random(20260728)

def lcm(a,b): return a*b//gcd(a,b)

def mat_inv(A,pp):
    n=len(A)
    aug=[A[i][:]+[1 if j==i else 0 for j in range(n)] for i in range(n)]
    R,piv=rref(aug,2*n,pp)
    assert piv[:n]==list(range(n)), "singular"
    return [[R[i][n+j] for j in range(n)] for i in range(n)]

def close_with_words(gens):
    """Exact closure; returns elems, index dict, and word[i]=(parent,gen) with elems[i]=elems[parent]*gens[gen]."""
    idx={EYE:0}; elems=[EYE]; word=[None]; frontier=[0]
    while frontier:
        nf=[]
        for ii in frontier:
            M=elems[ii]
            for gi,G in enumerate(gens):
                P=mmul(M,G)
                if P not in idx:
                    idx[P]=len(elems); word.append((ii,gi)); elems.append(P); nf.append(idx[P])
        frontier=nf
    return elems, idx, word

def sym_power(g,a,pp):
    """Sym^a(g), monomial basis x^(a-k) y^k, column k = image of that basis vector."""
    al,be,ga,de=g
    from math import comb
    S=[[0]*(a+1) for _ in range(a+1)]
    for k in range(a+1):
        for i in range(a-k+1):
            cp=comb(a-k,i)*pow(al,a-k-i,pp)%pp*pow(ga,i,pp)%pp
            for j in range(k+1):
                cq=comb(k,j)*pow(be,k-j,pp)%pp*pow(de,j,pp)%pp
                S[i+j][k]=(S[i+j][k]+cp*cq)%pp
    return S

class Grp: pass

def build(name,gens,expected):
    G=Grp(); G.name=name
    elems,idx,word=close_with_words(gens)
    G.elems=elems; G.idx=idx; G.word=word; G.gens=gens
    N=len(elems); G.N=N; G.expected=expected; G.order_ok=(N==expected)
    # multiplication table (exact dict lookup, no tolerance)
    mul=[[idx[mmul(elems[i],elems[j])] for j in range(N)] for i in range(N)]
    G.mul=mul
    inv=[0]*N
    for i in range(N):
        r=[j for j in range(N) if mul[i][j]==0]; assert len(r)==1; inv[i]=r[0]
    G.inv=inv
    # element orders + power map
    order=[0]*N; powmap=[]
    for i in range(N):
        c=0; row=[0]
        for t in range(1,4*N+2):
            c=mul[c][i]; row.append(c)
            if c==0 and order[i]==0: order[i]=t
        order[i]=order[i] or 1
        powmap.append(row)
    G.order_el=order; G.powmap=powmap
    G.exponent=reduce(lcm,order,1)
    G.has_minus_I = ((p-1,0,0,p-1) in idx)
    # conjugacy classes
    seen=[False]*N; classes=[]
    for i in range(N):
        if seen[i]: continue
        orb=sorted({mul[mul[g][i]][inv[g]] for g in range(N)})
        for x in orb: seen[x]=True
        classes.append(orb)
    classes.sort(key=lambda c:(0 if 0 in c else 1, order[c[0]], len(c), c[0]))
    G.classes=classes; G.nc=len(classes)
    clsof=[0]*N
    for ci,c in enumerate(classes):
        for x in c: clsof[x]=ci
    G.clsof=clsof
    G.csize=[len(c) for c in classes]
    G.crep=[c[0] for c in classes]
    G.corder=[order[c[0]] for c in classes]
    G.class_sq=[clsof[mul[c[0]][c[0]]] for c in classes]
    G.class_inv=[clsof[inv[c[0]]] for c in classes]
    # eigenvalue exponent s_i : lambda = zeta_M^{s_i}, M = exponent
    M=G.exponent; G.M=M
    zM=zeta(M) if L%M==0 else None
    assert zM is not None
    G.zM=zM
    s=[]
    for ci in range(G.nc):
        t=tr(elems[G.crep[ci]]); n=G.corder[ci]; zn=zeta(n)
        found=None
        for m in range(n):
            lam=pow(zn,m,p)
            if (lam*lam - t*lam + 1)%p==0: found=m; break
        assert found is not None, (name,ci)
        s.append(found*(M//n)%M)
    G.s=s
    return G

def chi_sym_exact(G,a):
    """chi_a on each class as exact Cyc in Z[zeta_M]."""
    M=G.M; out=[]
    for ci in range(G.nc):
        z=Cyc(M)
        for k in range(a+1):
            z.c[((a-2*k)*G.s[ci])%M]+=1
        out.append(z)
    return out

def ip_exact(G,f,gfun,conj=True):
    """(1/|G|) sum_i |C_i| f_i * (conj) g_i  -> Fraction (must be exact)."""
    M=G.M; acc=Cyc(M)
    for ci in range(G.nc):
        t=f[ci]*(gfun[ci].conj() if conj else gfun[ci])
        acc=acc+t.scal(G.csize[ci])
    return acc.to_rational()/G.N

# -------- irreducible characters: isotypic split of V_a by class sums (exact GF(p)) --------
def restrict(G,rho_gens_all,B):
    """Given invariant subspace basis B (d x k), return the restricted rep of all elements."""
    d=len(B); k=len(B[0])
    Bt=[[B[i][j] for i in range(d)] for j in range(k)]
    R,piv=rref(Bt,d,p)
    assert len(R)==k
    Bsel=[[B[r][j] for j in range(k)] for r in piv]
    Binv=mat_inv(Bsel,p)
    out={}
    for gi,rg in enumerate(rho_gens_all):
        Y=mat_mul(rg,B,p)
        Ysel=[[Y[r][j] for j in range(k)] for r in piv]
        Rg=mat_mul(Binv,Ysel,p)
        assert mat_mul(rg,B,p)==mat_mul(B,Rg,p), "not invariant"
        out[gi]=Rg
    return out

def commutant_dim(gen_mats,k):
    """dim {X in M_k : X*A = A*X for all generators A}, exact over GF(p)."""
    rows=[]
    for A in gen_mats:
        for i in range(k):
            for j in range(k):
                v=[0]*(k*k)
                for t in range(k):
                    v[i*k+t]=(v[i*k+t]+A[t][j])%p     # (X A)_{ij}
                    v[t*k+j]=(v[t*k+j]-A[i][t])%p     # -(A X)_{ij}
                rows.append(v)
    return k*k - rank_gf(rows,k*k,p)

def isotypic_split(G,a):
    """Decompose V_a into isotypic components using class sums. Exact. Returns list of (k, R_all_elems)."""
    d=a+1
    rgens=[sym_power(g,a,p) for g in G.gens]
    rho=[None]*G.N; rho[0]=mat_eye(d)
    for i in range(1,G.N):
        par,gi=G.word[i]; rho[i]=mat_mul(rho[par],rgens[gi],p)
    G._rho_cache=getattr(G,"_rho_cache",{}); G._rho_cache[a]=rho
    subs=[mat_eye(d)]                       # list of d x k bases
    for ci in range(G.nc):
        Z=[[0]*d for _ in range(d)]
        for x in G.classes[ci]: Z=mat_add(Z,rho[x],p)
        new=[]
        for B in subs:
            k=len(B[0])
            if k==1: new.append(B); continue
            Bt=[[B[i][j] for i in range(d)] for j in range(k)]
            R,piv=rref(Bt,d,p); Bsel=[[B[r][j] for j in range(k)] for r in piv]
            Binv=mat_inv(Bsel,p)
            Y=mat_mul(Z,B,p); Ysel=[[Y[r][j] for j in range(k)] for r in piv]
            Mm=mat_mul(Binv,Ysel,p)
            assert mat_mul(Z,B,p)==mat_mul(B,Mm,p)
            f=minpoly_of_matrix(Mm,p)
            if len(f)<=2: new.append(B); continue     # scalar on this piece
            roots=poly_roots_split(f,p,RNG)
            assert len(roots)==len(f)-1, ("minpoly of class sum did not split",G.name,a,ci,f)
            for lam in roots:
                Ml=mat_add(Mm,mat_scal(mat_eye(k),(-lam)%p,p),p)
                ns=nullspace(Ml,k,p)
                if not ns: continue
                Bn=[[sum(B[i][t]*ns[c][t] for t in range(k))%p for c in range(len(ns))] for i in range(d)]
                new.append(Bn)
        subs=new
    out=[]
    for B in subs:
        k=len(B[0])
        rres=restrict(G,[rho[G.idx[g]] for g in G.gens],B)
        gen_mats=[rres[gi] for gi in range(len(G.gens))]
        m2=commutant_dim(gen_mats,k)
        m=int(round(m2**0.5)); assert m*m==m2,(G.name,a,k,m2)
        # full restricted rep by word propagation
        Rall=[None]*G.N; Rall[0]=mat_eye(k)
        for i in range(1,G.N):
            par,gi=G.word[i]; Rall[i]=mat_mul(Rall[par],gen_mats[gi],p)
        out.append({"k":k,"mult":m,"dim":k//m,"Rall":Rall})
        assert k%m==0
    return out,rho

def lift_character(G,theta_gfp,d):
    """theta_gfp[class] in GF(p) is an irreducible character; lift to exact Cyc in Z[zeta_M]."""
    M=G.M; zM=G.zM; invM=pow(M,p-2,p); out=[]
    for ci in range(G.nc):
        g=G.crep[ci]
        vals=[theta_gfp[G.clsof[G.powmap[g][t%G.order_el[g]]]] for t in range(M)]
        z=Cyc(M); tot=0
        for l in range(M):
            acc=0
            for t in range(M):
                acc=(acc+vals[t]*pow(zM,(-l*t)%M,p))%p
            n_l=acc*invM%p
            assert n_l<=d, ("lift failed",G.name,ci,l,n_l)
            z.c[l]=n_l; tot+=n_l
        assert tot==d, ("multiplicities do not sum to dim",G.name,ci,tot,d)
        out.append(z)
    return out
