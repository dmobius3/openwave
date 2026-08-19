"""
Pre-implementation validation: construct all 9 irreps of 2I.
V0-V5 = Sym^0 through Sym^5, V7 = sigma(Sym^1), V8 = sigma(Sym^2),
V6 extracted from Sym^6 by projection.
"""
from fractions import Fraction
import json
from math import comb

# ---- Q(phi) ----
class QG:
    __slots__ = ('a','b')
    def __init__(s, a=0, b=0):
        s.a = a if isinstance(a, Fraction) else Fraction(a)
        s.b = b if isinstance(b, Fraction) else Fraction(b)
    def __eq__(s,o):
        if isinstance(o,(int,Fraction)): return s.a==o and s.b==0
        return s.a==o.a and s.b==o.b
    def __hash__(s): return hash((s.a,s.b))
    def __add__(s,o):
        if isinstance(o,(int,Fraction)): return QG(s.a+o,s.b)
        return QG(s.a+o.a,s.b+o.b)
    def __radd__(s,o): return s+o
    def __sub__(s,o):
        if isinstance(o,(int,Fraction)): return QG(s.a-o,s.b)
        return QG(s.a-o.a,s.b-o.b)
    def __rsub__(s,o):
        if isinstance(o,(int,Fraction)): return QG(o-s.a,-s.b)
        return QG(o.a-s.a,o.b-s.b)
    def __neg__(s): return QG(-s.a,-s.b)
    def __mul__(s,o):
        if isinstance(o,(int,Fraction)): return QG(s.a*o,s.b*o)
        return QG(s.a*o.a+s.b*o.b, s.a*o.b+s.b*o.a+s.b*o.b)
    def __rmul__(s,o): return s*o
    def norm(s): return s.a*s.a+s.a*s.b-s.b*s.b
    def galois(s): return QG(s.a+s.b,-s.b)
    def __truediv__(s,o):
        if isinstance(o,(int,Fraction)): return QG(s.a/o,s.b/o)
        n=o.norm(); c=o.galois(); num=s*c; return QG(num.a/n,num.b/n)
    def is_zero(s): return s.a==0 and s.b==0
    def __repr__(s):
        if s.b == 0: return str(s.a)
        return f"({s.a}+{s.b}phi)"
    def to_float(s): return float(s.a)+float(s.b)*(1+5**0.5)/2

# ---- Q(phi,i) ----
class GC:
    __slots__ = ('re','im')
    def __init__(s, re=None, im=None):
        s.re = re if isinstance(re,QG) else QG(re if re is not None else 0)
        s.im = im if isinstance(im,QG) else QG(im if im is not None else 0)
    def __add__(s,o):
        if not isinstance(o,GC): o=GC(o if isinstance(o,QG) else QG(o))
        return GC(s.re+o.re,s.im+o.im)
    def __radd__(s,o): return s+o
    def __sub__(s,o):
        if not isinstance(o,GC): o=GC(o if isinstance(o,QG) else QG(o))
        return GC(s.re-o.re,s.im-o.im)
    def __rsub__(s,o):
        if not isinstance(o,GC): o=GC(o if isinstance(o,QG) else QG(o))
        return GC(o.re-s.re,o.im-s.im)
    def __neg__(s): return GC(-s.re,-s.im)
    def __mul__(s,o):
        if isinstance(o,(int,Fraction)): return GC(s.re*o,s.im*o)
        if isinstance(o,QG): return GC(s.re*o,s.im*o)
        return GC(s.re*o.re-s.im*o.im, s.re*o.im+s.im*o.re)
    def __rmul__(s,o): return s*o
    def conj(s): return GC(s.re,-s.im)
    def mod_sq(s): return s.re*s.re+s.im*s.im
    def __truediv__(s,o):
        if isinstance(o,(int,Fraction)): return GC(s.re/o,s.im/o)
        if isinstance(o,QG): return GC(s.re/o,s.im/o)
        d=o.mod_sq(); return GC((s.re*o.re+s.im*o.im)/d,(s.im*o.re-s.re*o.im)/d)
    def is_zero(s): return s.re.is_zero() and s.im.is_zero()
    def __eq__(s,o):
        if isinstance(o,(int,Fraction)): return s.re==o and s.im.is_zero()
        if isinstance(o,QG): return s.re==o and s.im.is_zero()
        return s.re==o.re and s.im==o.im
    def galois(s): return GC(s.re.galois(),s.im.galois())
    def __repr__(s): return f"({s.re}+{s.im}i)"

GC0=GC(); GC1=GC(QG(1)); GCI=GC(QG(0),QG(1))

def mz(m,n): return [[GC() for _ in range(n)] for _ in range(m)]
def mmul(A,B):
    m,p,n=len(A),len(A[0]),len(B[0]); C=mz(m,n)
    for i in range(m):
        for j in range(n):
            s=GC()
            for k in range(p): s=s+A[i][k]*B[k][j]
            C[i][j]=s
    return C
def mid(n):
    M=mz(n,n)
    for i in range(n): M[i][i]=GC1
    return M
def mtr(A): return sum((A[i][i] for i in range(len(A))),GC())
def mct(A):
    m,n=len(A),len(A[0]); return [[A[j][i].conj() for j in range(m)] for i in range(n)]
def madd(A,B): return [[A[i][j]+B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def msc(c,A): return [[c*A[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def meq(A,B): return all(A[i][j]==B[i][j] for i in range(len(A)) for j in range(len(A[0])))

# ---- Quaternion ----
class Q4:
    __slots__=('w','x','y','z')
    def __init__(s,w,x,y,z):
        s.w=w if isinstance(w,QG) else QG(w); s.x=x if isinstance(x,QG) else QG(x)
        s.y=y if isinstance(y,QG) else QG(y); s.z=z if isinstance(z,QG) else QG(z)
    def __eq__(s,o): return s.w==o.w and s.x==o.x and s.y==o.y and s.z==o.z
    def __hash__(s): return hash((s.w,s.x,s.y,s.z))
    def __mul__(s,o):
        return Q4(s.w*o.w-s.x*o.x-s.y*o.y-s.z*o.z,
                  s.w*o.x+s.x*o.w+s.y*o.z-s.z*o.y,
                  s.w*o.y-s.x*o.z+s.y*o.w+s.z*o.x,
                  s.w*o.z+s.x*o.y-s.y*o.x+s.z*o.w)
    def __neg__(s): return Q4(-s.w,-s.x,-s.y,-s.z)
    def sk(s):
        k=[]
        for c in [s.w,s.x,s.y,s.z]: k.append(int(c.a*2)); k.append(int(c.b*2))
        return tuple(k)

def pg(strs):
    comps=[]
    for s in strs:
        s=s.strip(); inner=s[1:-3]; parts=inner.replace(' ','').replace('*phi','p')
        toks=[]; cur=''
        for ch in parts:
            if ch in '+-' and cur: toks.append(cur); cur=ch
            else: cur+=ch
        if cur: toks.append(cur)
        a=0;b=0
        for tok in toks:
            if 'p' in tok:
                c2=tok.replace('p','')
                b=1 if c2 in ('','+') else (-1 if c2=='-' else int(c2))
            else: a=int(tok)
        comps.append(QG(Fraction(a,2),Fraction(b,2)))
    return Q4(*comps)

def eg(gp):
    g1=pg(gp['generators'][0]); g2=pg(gp['generators'][1])
    e=Q4(QG(1),QG(0),QG(0),QG(0)); elems={e,g1,g2}
    while True:
        new=set()
        for a in list(elems):
            for b in [g1,g2,-g1,-g2]:
                for c in [a*b,b*a]:
                    if c not in elems and c not in new: new.add(c)
        if not new: break
        elems.update(new)
    return sorted(elems,key=lambda q:q.sk())

def q2su2(q):
    return [[GC(q.w,q.x),GC(q.y,q.z)],[GC(-q.y,q.z),GC(q.w,-q.x)]]

def sym_pow(M,n):
    d=n+1; R=mz(d,d); a,b,g,dl=M[0][0],M[0][1],M[1][0],M[1][1]
    for k in range(d):
        for ja in range(n-k+1):
            for jb in range(k+1):
                j=ja+jb
                if j>n: continue
                coeff=GC(QG(comb(n-k,ja)*comb(k,jb)))
                val=coeff
                for _ in range(n-k-ja): val=val*a
                for _ in range(ja): val=val*g
                for _ in range(k-jb): val=val*b
                for _ in range(jb): val=val*dl
                R[j][k]=R[j][k]+val
    return R

def mgal(M): return [[M[i][j].galois() for j in range(len(M[0]))] for i in range(len(M))]

def det_gc(M):
    n=len(M); A=[[M[i][j] for j in range(n)] for i in range(n)]
    d=GC1
    for col in range(n):
        piv=None
        for row in range(col,n):
            if not A[row][col].is_zero(): piv=row; break
        if piv is None: return GC0
        if piv!=col: A[col],A[piv]=A[piv],A[col]; d=d*GC(QG(-1))
        d=d*A[col][col]; sc=A[col][col]
        for j in range(col,n): A[col][j]=A[col][j]/sc
        for row in range(col+1,n):
            f=A[row][col]
            if not f.is_zero():
                for j in range(col,n): A[row][j]=A[row][j]-f*A[col][j]
    return d


def main():
    with open('m8_5a_packet.json') as f: gp=json.load(f)
    with open('m8_8_construction_packet.json') as f: cp=json.load(f)

    se=eg(gp); assert len(se)==120
    s_id=cp['abstract_generators']['s']; t_id=cp['abstract_generators']['t']; e_id=119
    e2r={q:i for i,q in enumerate(se)}
    mt=[[e2r[se[i]*se[j]] for j in range(120)] for i in range(120)]

    print("Building SU(2) representations...")
    su2=[q2su2(q) for q in se]
    assert meq(su2[e_id],mid(2))
    for i,j in [(s_id,t_id),(t_id,s_id),(s_id,s_id)]:
        k=mt[i][j]; assert meq(mmul(su2[i],su2[j]),su2[k])
    print("SU(2) homomorphism: ✓")

    # Build Sym^0 through Sym^6
    print("\nBuilding symmetric powers...")
    sym_reps={}
    for n in range(7):
        sym_reps[n]=[sym_pow(su2[g],n) for g in range(120)]
        # Verify homomorphism
        for i,j in [(s_id,t_id),(0,1)]:
            k=mt[i][j]
            assert meq(mmul(sym_reps[n][i],sym_reps[n][j]),sym_reps[n][k]), f"Sym{n} hom fail"
        print(f"  Sym^{n} (dim {n+1}): ✓")

    # Build Galois conjugates of Sym^1 and Sym^2
    gal_reps={}
    for n in [1,2]:
        gal_reps[n]=[mgal(sym_reps[n][g]) for g in range(120)]
        for i,j in [(s_id,t_id),(0,1)]:
            k=mt[i][j]
            assert meq(mmul(gal_reps[n][i],gal_reps[n][j]),gal_reps[n][k]), f"GalSym{n} hom fail"
        print(f"  σ(Sym^{n}) (dim {n+1}): ✓")

    # Compute characters
    print("\nComputing characters...")
    chars = {}
    names = ['V0','V1','V2','V3','V4','V5','V7','V8']
    rep_sources = {
        'V0': sym_reps[0], 'V1': sym_reps[1], 'V2': sym_reps[2],
        'V3': sym_reps[3], 'V4': sym_reps[4], 'V5': sym_reps[5],
        'V7': gal_reps[1], 'V8': gal_reps[2]
    }
    for name, reps in rep_sources.items():
        ch = []
        for g in range(120):
            tr = mtr(reps[g])
            assert tr.im.is_zero(), f"{name}: complex character at {g}"
            ch.append(tr.re)
        chars[name] = ch

    # Compute V6 character: V6 = Sym^6 - V8
    ch_sym6 = []
    for g in range(120):
        tr = mtr(sym_reps[6][g])
        assert tr.im.is_zero()
        ch_sym6.append(tr.re)

    ch_v6 = [ch_sym6[g] - chars['V8'][g] for g in range(120)]
    chars['V6_expected'] = ch_v6

    # Verify V6 character has norm 1 (irreducible)
    norm_v6 = sum(ch_v6[g]*ch_v6[g] for g in range(120)) / 120
    print(f"  V6 character norm: {norm_v6} (should be 1)")
    assert norm_v6 == QG(1)

    # Verify V6 is orthogonal to all others
    for name in ['V0','V1','V2','V3','V4','V5','V7','V8']:
        ip = sum(ch_v6[g]*chars[name][g] for g in range(120)) / 120
        assert ip.is_zero(), f"<V6,{name}> = {ip}"
    print("  V6 character orthogonal to all others: ✓")

    # Construct V6 by projection from Sym^6
    print("\nExtracting V6 from Sym^6 by projection...")
    # P6 = (4/120) * sum_g chi_V6(g) * rho_Sym6(g)
    P6 = mz(7,7)
    for g in range(120):
        cv = ch_v6[g]  # QG value
        if cv.is_zero(): continue
        P6 = madd(P6, msc(GC(cv/120*4), sym_reps[6][g]))

    # Verify P6 is idempotent: P6^2 = P6
    P6sq = mmul(P6, P6)
    assert meq(P6, P6sq), "P6 is not idempotent!"
    print("  P6 idempotent: ✓")

    # Find rank of P6 (should be 4)
    # Find 4 independent columns
    rank_p6 = 0
    basis_cols = []
    for j in range(7):
        col = [P6[i][j] for i in range(7)]
        if all(c.is_zero() for c in col): continue
        basis_cols.append(j)
        rank_p6 += 1
        if rank_p6 == 4: break

    # Actually, we need to check linear independence. Let me find pivot columns.
    # Gaussian elimination on P6 to find rank and pivot columns.
    A = [[P6[i][j] for j in range(7)] for i in range(7)]
    pivots = []
    r = 0
    for col in range(7):
        piv = None
        for row in range(r, 7):
            if not A[row][col].is_zero():
                piv = row; break
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        sc = A[r][col]
        for j in range(7): A[r][j] = A[r][j] / sc
        for row in range(7):
            if row == r: continue
            f = A[row][col]
            if not f.is_zero():
                for j in range(7): A[row][j] = A[row][j] - f * A[r][j]
        pivots.append(col)
        r += 1

    print(f"  P6 rank: {r}, pivot cols: {pivots}")
    assert r == 4, f"Expected rank 4, got {r}"

    # Extract basis: columns of P6 at pivot positions
    C = [[P6[i][pivots[j]] for j in range(4)] for i in range(7)]  # 7x4 matrix

    # Compute C^dag C (4x4)
    Cdag = mct(C)  # 4x7
    CdC = mmul(Cdag, C)  # 4x4

    # Invert CdC
    # Build augmented matrix [CdC | I4]
    aug = [[CdC[i][j] for j in range(4)] + [GC1 if i==j else GC0 for j in range(4)] for i in range(4)]
    for col in range(4):
        piv = None
        for row in range(col, 4):
            if not aug[row][col].is_zero():
                piv = row; break
        assert piv is not None
        aug[col], aug[piv] = aug[piv], aug[col]
        sc = aug[col][col]
        for j in range(8): aug[col][j] = aug[col][j] / sc
        for row in range(4):
            if row == col: continue
            f = aug[row][col]
            if not f.is_zero():
                for j in range(8): aug[row][j] = aug[row][j] - f * aug[col][j]
    CdC_inv = [[aug[i][4+j] for j in range(4)] for i in range(4)]

    # Left pseudoinverse: C+ = (C^dag C)^{-1} C^dag  (4x7)
    Cplus = mmul(CdC_inv, Cdag)

    # Build V6 representations: for each g, R(g) = C+ * rho_Sym6(g) * C  (4x4)
    v6_reps = []
    for g in range(120):
        rC = mmul(sym_reps[6][g], C)  # 7x4
        Rg = mmul(Cplus, rC)  # 4x4
        v6_reps.append(Rg)

    # Verify homomorphism
    for i,j in [(s_id,t_id),(t_id,s_id),(0,1),(s_id,s_id)]:
        k = mt[i][j]
        assert meq(mmul(v6_reps[i],v6_reps[j]), v6_reps[k]), f"V6 hom fail at ({i},{j})"
    print("  V6 homomorphism: ✓")

    # Verify V6 character
    chars['V6'] = []
    for g in range(120):
        tr = mtr(v6_reps[g])
        assert tr.im.is_zero(), f"V6: complex char at {g}"
        chars['V6'].append(tr.re)
        assert tr.re == ch_v6[g], f"V6 char mismatch at {g}"
    print("  V6 character matches expected: ✓")

    # Add V6 to rep sources
    rep_sources['V6'] = v6_reps
    names = ['V0','V1','V2','V3','V4','V5','V6','V7','V8']

    # Full character orthogonality
    print("\nVerifying full character orthogonality (9x9)...")
    for i, ni in enumerate(names):
        for j, nj in enumerate(names):
            ip = sum(chars[ni][g]*chars[nj][g] for g in range(120)) / 120
            if i == j:
                assert ip == QG(1), f"<{ni},{nj}> = {ip}"
            else:
                assert ip.is_zero(), f"<{ni},{nj}> = {ip}"
    print("Character orthogonality: ✓")

    # Verify unitarity: invariant positive-definite Hermitian form exists
    print("\nVerifying unitarity (invariant Hermitian form)...")
    for name in names:
        reps = rep_sources[name]
        dim = len(reps[0])
        H = mz(dim,dim)
        for g in range(120):
            H = madd(H, mmul(mct(reps[g]), reps[g]))
        H = msc(GC(QG(Fraction(1,120))), H)
        Hdag = mct(H)
        assert meq(H, Hdag), f"{name}: H not Hermitian"
        for i in range(dim):
            assert H[i][i].im.is_zero(), f"{name}: H[{i},{i}] has imag part"
            assert H[i][i].re.to_float() > 0, f"{name}: H[{i},{i}] not positive"
        for gen_id, gen_name in [(s_id, 's'), (t_id, 't')]:
            lhs = mmul(mct(reps[gen_id]), mmul(H, reps[gen_id]))
            assert meq(lhs, H), f"{name}: H not invariant under {gen_name}"
        dH = det_gc(H)
        assert dH.im.is_zero(), f"{name}: det(H) complex"
        assert dH.re.to_float() > 0, f"{name}: det(H) not positive"
        print(f"  {name}: Hermitian ✓, invariant ✓, det>0 ✓")

    # Row signatures
    print("\nRow signatures:")
    st_id = mt[s_id][t_id]
    sigs = {}
    for name in names:
        dim_v = chars[name][e_id]
        cs = chars[name][s_id]
        ct = chars[name][t_id]
        cst = chars[name][st_id]
        sig = (dim_v, cs, ct, cst)
        sig_key = tuple((c.a,c.b) for c in sig)
        assert sig_key not in sigs, f"Duplicate sig: {name} matches {sigs[sig_key]}"
        sigs[sig_key] = name
        print(f"  {name}: dim={dim_v}, χ(s)={cs}, χ(t)={ct}, χ(st)={cst}")

    print("All 9 row signatures distinct: ✓")

    # Galois pairs
    print("\nGalois structure:")
    for name in names:
        gal_ch = [chars[name][g].galois() for g in range(120)]
        # Find which irrep matches
        for name2 in names:
            if all(gal_ch[g] == chars[name2][g] for g in range(120)):
                if name == name2:
                    print(f"  {name}: self-conjugate")
                else:
                    print(f"  {name} <-> {name2}")
                break

    # Dimensions
    dims = sorted(int(chars[n][e_id].a) for n in names)
    assert dims == [1,2,2,3,3,4,4,5,6]
    print(f"\nDimensions: {dims} ✓")
    print(f"Sum of squares: {sum(d*d for d in dims)} = 120 ✓")

    print("\n=== ALL REPRESENTATION CHECKS PASSED ===")

if __name__ == '__main__':
    main()
