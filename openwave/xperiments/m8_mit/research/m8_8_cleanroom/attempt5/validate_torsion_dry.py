"""
Pre-implementation validation: evaluate boundary maps per irrep,
verify dd=0 and acyclicity of twisted complexes.
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
    def to_float(s): return float(s.a)+float(s.b)*(1+5**0.5)/2

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
def mgal(M): return [[M[i][j].galois() for j in range(len(M[0]))] for i in range(len(M))]

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


def build_all_irreps(se, mt, s_id, t_id, e_id):
    """Build all 9 irreps. Returns dict name -> list of 120 matrices."""
    su2 = [q2su2(q) for q in se]

    sym = {}
    for n in range(7):
        sym[n] = [sym_pow(su2[g], n) for g in range(120)]

    reps = {}
    reps['V0'] = sym[0]
    reps['V1'] = sym[1]
    reps['V2'] = sym[2]
    reps['V3'] = sym[3]
    reps['V4'] = sym[4]
    reps['V5'] = sym[5]
    reps['V7'] = [mgal(sym[1][g]) for g in range(120)]
    reps['V8'] = [mgal(sym[2][g]) for g in range(120)]

    # V6 by projection from Sym^6
    # character of V6 = character of Sym6 - character of V8
    ch_sym6 = []
    ch_v8 = []
    for g in range(120):
        tr6 = mtr(sym[6][g])
        tr8 = mtr(reps['V8'][g])
        ch_sym6.append(tr6.re)
        ch_v8.append(tr8.re)
    ch_v6 = [ch_sym6[g] - ch_v8[g] for g in range(120)]

    P6 = mz(7,7)
    for g in range(120):
        cv = ch_v6[g]
        if cv.is_zero(): continue
        P6 = madd(P6, msc(GC(cv / 120 * 4), sym[6][g]))

    # Find pivot columns
    A = [[P6[i][j] for j in range(7)] for i in range(7)]
    pivots = []; r = 0
    for col in range(7):
        piv = None
        for row in range(r, 7):
            if not A[row][col].is_zero(): piv = row; break
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        sc = A[r][col]
        for j in range(7): A[r][j] = A[r][j] / sc
        for row in range(7):
            if row == r: continue
            f = A[row][col]
            if not f.is_zero():
                for j in range(7): A[row][j] = A[row][j] - f * A[r][j]
        pivots.append(col); r += 1
    assert r == 4

    C = [[P6[i][pivots[j]] for j in range(4)] for i in range(7)]
    Cdag = mct(C)
    CdC = mmul(Cdag, C)
    aug = [[CdC[i][j] for j in range(4)] + [GC1 if i==j else GC0 for j in range(4)] for i in range(4)]
    for col in range(4):
        piv = None
        for row in range(col, 4):
            if not aug[row][col].is_zero(): piv = row; break
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
    Cplus = mmul(CdC_inv, Cdag)

    v6_reps = []
    for g in range(120):
        rC = mmul(sym[6][g], C)
        v6_reps.append(mmul(Cplus, rC))
    reps['V6'] = v6_reps

    return reps


def rank_gc(M):
    """Compute rank of a matrix over GC by Gaussian elimination."""
    m = len(M); n = len(M[0])
    A = [[M[i][j] for j in range(n)] for i in range(m)]
    r = 0
    for col in range(n):
        piv = None
        for row in range(r, m):
            if not A[row][col].is_zero(): piv = row; break
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        sc = A[r][col]
        for j in range(n): A[r][j] = A[r][j] / sc
        for row in range(m):
            if row == r: continue
            f = A[row][col]
            if not f.is_zero():
                for j in range(n): A[row][j] = A[row][j] - f * A[r][j]
        r += 1
    return r


def main():
    with open('m8_5a_packet.json') as f: gp = json.load(f)
    with open('m8_8_construction_packet.json') as f: cp = json.load(f)

    se = eg(gp); assert len(se) == 120
    s_id = cp['abstract_generators']['s']
    t_id = cp['abstract_generators']['t']
    e_id = 119
    mt = [[se.index(se[i]*se[j]) for j in range(120)] for i in range(120)]

    # Faster lookup
    e2r = {q: i for i, q in enumerate(se)}
    mt = [[e2r[se[i]*se[j]] for j in range(120)] for i in range(120)]

    print("Building all 9 irreps...")
    reps = build_all_irreps(se, mt, s_id, t_id, e_id)

    # Parse boundary maps
    d1_raw = cp['boundary_maps']['d1']  # 2x1 over Z[2I]
    d2_raw = cp['boundary_maps']['d2']  # 2x2 over Z[2I]
    d3_raw = cp['boundary_maps']['d3']  # 1x2 over Z[2I]

    names = ['V0','V1','V2','V3','V4','V5','V6','V7','V8']

    def eval_gr_entry(entry, rho, d):
        """Evaluate a group-ring element sum_i c_i * g_i as sum_i c_i * rho(g_i)."""
        result = mz(d, d)
        for coeff, eid in entry:
            result = madd(result, msc(GC(QG(coeff)), rho[eid]))
        return result

    def eval_bmap(bmap_raw, rho, d):
        """Evaluate a boundary map over Z[2I] into a block matrix over C.
        bmap_raw has shape (rows_gr, cols_gr) where each entry is a list of (coeff, eid) pairs.
        Result is a (rows_gr*d, cols_gr*d) matrix."""
        rows_gr = len(bmap_raw)
        cols_gr = len(bmap_raw[0])
        M = mz(rows_gr * d, cols_gr * d)
        for bi in range(rows_gr):
            for bj in range(cols_gr):
                block = eval_gr_entry(bmap_raw[bi][bj], rho, d)
                for i in range(d):
                    for j in range(d):
                        M[bi*d + i][bj*d + j] = block[i][j]
        return M

    print("\nEvaluating twisted complexes...")
    for name in names:
        rho = reps[name]
        d = len(rho[0])  # dimension

        M1 = eval_bmap(d1_raw, rho, d)  # 2d x d
        M2 = eval_bmap(d2_raw, rho, d)  # 2d x 2d
        M3 = eval_bmap(d3_raw, rho, d)  # d x 2d

        # Check dimensions
        assert len(M1) == 2*d and len(M1[0]) == d, f"{name}: M1 dims wrong"
        assert len(M2) == 2*d and len(M2[0]) == 2*d, f"{name}: M2 dims wrong"
        assert len(M3) == d and len(M3[0]) == 2*d, f"{name}: M3 dims wrong"

        # Check dd=0: M3*M2 = 0 (d x 2d) and M2*M1 = 0 (2d x d)
        # Actually: boundary direction is right action, c -> c.d_k
        # C_3 -> C_2 -> C_1 -> C_0
        # d3: C_3 -> C_2 maps 1 cell to 2 cells, so M3 is 1*d rows, 2*d cols => d x 2d
        # d2: C_2 -> C_1 maps 2 cells to 2 cells, so M2 is 2*d rows, 2*d cols
        # d1: C_1 -> C_0 maps 2 cells to 1 cell, so M1 is 2*d rows, d cols
        # Chain condition: d_{k-1} . d_k = 0 as right action
        # c in C_3, c.d3 in C_2, (c.d3).d2 in C_1
        # So we need d3.d2 = 0 and d2.d1 = 0 (matrix products in right action)
        # With row vectors: d3 is d×2d, d2 is 2d×2d, so d3·d2 is d×2d... wait
        # Actually: row vector c (1×d for C_3), c·M3 (1×2d for C_2), (c·M3)·M2 (1×2d... wrong)
        # Hmm: C_3 has rank 1, C_2 has rank 2, C_1 has rank 2, C_0 has rank 1
        # After tensoring with V_rho (dim d):
        # C_3 ⊗ V: d-dimensional, C_2 ⊗ V: 2d-dimensional
        # d3: C_3⊗V -> C_2⊗V, so c (row of length d) maps to c·M3 (row of length 2d)
        # d2: C_2⊗V -> C_1⊗V, so c (row of length 2d) maps to c·M2 (row of length 2d)
        # d1: C_1⊗V -> C_0⊗V, so c (row of length 2d) maps to c·M1 (row of length d)
        #
        # Chain condition: d2·d1 should be zero as a map C_2⊗V -> C_0⊗V
        # In matrix form: (c·M2)·M1 = c·(M2·M1) = 0 for all c
        # So M2·M1 = 0 (2d×d matrix)
        # And d3·d2: (c·M3)·M2 = c·(M3·M2) = 0
        # So M3·M2 = 0 (d×2d matrix)

        prod32 = mmul(M3, M2)
        all_zero_32 = all(prod32[i][j].is_zero()
                         for i in range(d) for j in range(2*d))
        assert all_zero_32, f"{name}: M3·M2 ≠ 0!"

        prod21 = mmul(M2, M1)
        all_zero_21 = all(prod21[i][j].is_zero()
                         for i in range(2*d) for j in range(d))
        assert all_zero_21, f"{name}: M2·M1 ≠ 0!"

        # Compute ranks
        r3 = rank_gc(M3)
        r2 = rank_gc(M2)
        r1 = rank_gc(M1)

        if name == 'V0':
            # Trivial rep: all group elements map to [1].
            # Augmented complex: same as augmented homology.
            # M1 = [[s-1],[t-1]] = [[0],[0]] => rank 0
            # M2 should have rank 1 (since augmented d2 has det=-1)
            # M3 should have rank 0
            # NOT acyclic. This is correct.
            print(f"  {name} (dim {d}): ∂∂=0 ✓, ranks ({r3},{r2},{r1}), NON-ACYCLIC (trivial rep) ✓")
            assert r3 == 0 or r2 == 0 or r1 == 0, "V0 should be non-acyclic"
        else:
            # Nontrivial irreps should be acyclic
            # For acyclicity: rank(M3) + rank(M2) = 2d (exactness at C_2)
            #                 rank(M2) + rank(M1) = 2d (exactness at C_1)
            #                 rank(M3) = d (M3 is injective on d-dim space)
            #                 rank(M1) = d (M1 is surjective onto d-dim space)
            acyclic = (r3 == d) and (r1 == d) and (r2 == d)
            assert acyclic, f"{name}: NOT ACYCLIC! ranks=({r3},{r2},{r1}), expected ({d},{d},{d})"
            print(f"  {name} (dim {d}): ∂∂=0 ✓, ranks ({r3},{r2},{r1}) = ({d},{d},{d}) ACYCLIC ✓")

    print("\n=== ALL TWISTED-COMPLEX CHECKS PASSED ===")


if __name__ == '__main__':
    main()
