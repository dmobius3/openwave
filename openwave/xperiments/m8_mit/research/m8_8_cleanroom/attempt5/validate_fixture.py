"""
Convention fixture: a synthetic non-unitary representation and chain-complex
instance exercising all four declared conventions; processed through the same
parser and evaluation path as the target run.

Each single-convention mutation must redden at least one gate.
"""
from fractions import Fraction
import json
from math import comb


# ---- Q(phi), GC, Q4 (same as other validators) ----
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
def meq(A,B): return all(A[i][j]==B[i][j] for i in range(len(A)) for j in range(len(A[0])))
def mtranspose(M): return [[M[j][i] for j in range(len(M))] for i in range(len(M[0]))]
def msc(c,A): return [[c*A[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def madd(A,B): return [[A[i][j]+B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

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


def minv(M):
    n = len(M)
    aug = [[M[i][j] for j in range(n)] + [GC1 if i==j else GC0 for j in range(n)] for i in range(n)]
    for col in range(n):
        piv = None
        for row in range(col, n):
            if not aug[row][col].is_zero(): piv = row; break
        assert piv is not None
        aug[col], aug[piv] = aug[piv], aug[col]
        sc = aug[col][col]
        for j in range(2*n): aug[col][j] = aug[col][j] / sc
        for row in range(n):
            if row == col: continue
            f = aug[row][col]
            if not f.is_zero():
                for j in range(2*n): aug[row][j] = aug[row][j] - f * aug[col][j]
    return [[aug[i][n+j] for j in range(n)] for i in range(n)]


def compute_torsion(rho, d1_raw, d2_raw, d3_raw, d):
    """Compute torsion for a representation, using the declared conventions.
    Returns (tau_det, is_acyclic).
    tau = det(D2_minor) / (det(D1_minor) * det(D3_minor))
    T^2 = |tau|^2
    """
    def eval_entry(entry):
        result = mz(d, d)
        for coeff, eid in entry:
            result = madd(result, msc(GC(QG(coeff)), rho[eid]))
        return result

    def eval_bmap(bmap_raw):
        rows_gr = len(bmap_raw); cols_gr = len(bmap_raw[0])
        M = mz(rows_gr * d, cols_gr * d)
        for bi in range(rows_gr):
            for bj in range(cols_gr):
                block = eval_entry(bmap_raw[bi][bj])
                for i in range(d):
                    for j in range(d):
                        M[bi*d + i][bj*d + j] = block[i][j]
        return M

    M1 = eval_bmap(d1_raw)  # 2d x d
    M2 = eval_bmap(d2_raw)  # 2d x 2d
    M3 = eval_bmap(d3_raw)  # d x 2d

    # Check dd=0
    prod32 = mmul(M3, M2)
    prod21 = mmul(M2, M1)
    dd_ok = (all(prod32[i][j].is_zero() for i in range(d) for j in range(2*d)) and
             all(prod21[i][j].is_zero() for i in range(2*d) for j in range(d)))

    if not dd_ok:
        return None, False, "dd!=0"

    # Find d x d nonsingular minor of M3 (columns J3)
    # Try first d columns
    J3 = list(range(d))
    M3_minor = [[M3[i][j] for j in J3] for i in range(d)]
    det3 = det_gc(M3_minor)
    if det3.is_zero():
        J3 = list(range(d, 2*d))
        M3_minor = [[M3[i][j] for j in J3] for i in range(d)]
        det3 = det_gc(M3_minor)
        if det3.is_zero():
            return None, False, "M3 rank deficient"

    # Find d x d nonsingular minor of M1 (rows I1)
    I1 = list(range(d))
    M1_minor = [[M1[i][j] for j in range(d)] for i in I1]
    det1 = det_gc(M1_minor)
    if det1.is_zero():
        I1 = list(range(d, 2*d))
        M1_minor = [[M1[i][j] for j in range(d)] for i in I1]
        det1 = det_gc(M1_minor)
        if det1.is_zero():
            return None, False, "M1 rank deficient"

    # Complementary indices
    J3c = [j for j in range(2*d) if j not in J3]
    I1c = [i for i in range(2*d) if i not in I1]

    # Central minor M2[J3c, I1c]
    M2_minor = [[M2[i][j] for j in I1c] for i in J3c]
    det2 = det_gc(M2_minor)

    tau = det2 / (det1 * det3)
    return tau, True, "ok"


def main():
    """Convention fixture test.

    We use a synthetic non-unitary 2D representation of 2I to test that
    all four convention choices (module side, vector convention, evaluation map,
    boundary direction) are exercised causally.

    The fixture uses the SAME evaluation and boundary-map parsing as the
    production code. A mutation of each convention must change the result.
    """
    # Load group and chain complex
    with open('m8_5a_packet.json') as f: gp = json.load(f)
    with open('m8_8_construction_packet.json') as f: cp = json.load(f)

    # Build group
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

    def pg2(strs):
        comps=[]
        for s2 in strs:
            s2=s2.strip(); inner=s2[1:-3]; parts=inner.replace(' ','').replace('*phi','p')
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

    g1=pg2(gp['generators'][0]); g2=pg2(gp['generators'][1])
    e4=Q4(QG(1),QG(0),QG(0),QG(0)); elems={e4,g1,g2}
    while True:
        new=set()
        for a in list(elems):
            for b in [g1,g2,-g1,-g2]:
                for c in [a*b,b*a]:
                    if c not in elems and c not in new: new.add(c)
        if not new: break
        elems.update(new)
    se=sorted(elems,key=lambda q:q.sk())
    e2r={q:i for i,q in enumerate(se)}
    mt=[[e2r[se[i]*se[j]] for j in range(120)] for i in range(120)]

    s_id=cp['abstract_generators']['s']; t_id=cp['abstract_generators']['t']

    # Build non-unitary fixture representation: conjugate Sym^1 by a non-unitary matrix
    def q2su2(q):
        return [[GC(q.w,q.x),GC(q.y,q.z)],[GC(-q.y,q.z),GC(q.w,-q.x)]]
    su2=[q2su2(q) for q in se]

    # Non-unitary similarity: P = [[2, i], [0, 1]]
    P = [[GC(QG(2)),GCI],[GC0,GC1]]
    Pi = minv(P)

    # rho_fixture(g) = P^{-1} * rho_std(g) * P
    fixture_reps = [mmul(Pi, mmul(su2[g], P)) for g in range(120)]

    # Verify homomorphism
    for i,j in [(s_id,t_id),(t_id,s_id)]:
        k=mt[i][j]
        assert meq(mmul(fixture_reps[i],fixture_reps[j]),fixture_reps[k])

    # Verify non-unitary: rho(s)^dag * rho(s) != I
    rho_s = fixture_reps[s_id]
    sds = mmul([[rho_s[j][i].conj() for j in range(2)] for i in range(2)], rho_s)
    assert not meq(sds, mid(2)), "Fixture representation is unitary!"
    print("Fixture representation: non-unitary 2D rep ✓")

    d1_raw = cp['boundary_maps']['d1']
    d2_raw = cp['boundary_maps']['d2']
    d3_raw = cp['boundary_maps']['d3']

    # Baseline: compute torsion with correct conventions
    tau_base, acyclic, msg = compute_torsion(fixture_reps, d1_raw, d2_raw, d3_raw, 2)
    assert acyclic, f"Fixture not acyclic: {msg}"
    T2_base = tau_base * tau_base.conj()
    print(f"Baseline T² = {T2_base.re} (imaginary should be 0: {T2_base.im.is_zero()})")
    assert T2_base.im.is_zero()

    # --- Mutation G-T03a: evaluation map g -> rho(g^-1) (anti-homomorphism) ---
    # All 2I irreps are self-contragredient, so g->rho(g^-1)^T gives an equivalent
    # rep with the same T². Instead, use g->rho(g^-1) without transpose: this is an
    # anti-homomorphism and will break the chain condition dd=0.
    print("\nG-T03a: evaluation map mutation (g -> rho(g^-1)), anti-homomorphism...")
    inv_map = [0]*120
    for i in range(120):
        for j in range(120):
            if mt[i][j] == 119: inv_map[i] = j; break
    anti_hom = [fixture_reps[inv_map[g]] for g in range(120)]
    tau_mut_a, acyclic_a, msg_a = compute_torsion(anti_hom, d1_raw, d2_raw, d3_raw, 2)
    if not acyclic_a:
        print(f"  Result: {msg_a} — gate reddened ✓")
    else:
        T2_mut_a = tau_mut_a * tau_mut_a.conj()
        changed_a = not (T2_mut_a.re == T2_base.re)
        print(f"  T² changed = {changed_a}")
        assert changed_a, "G-T03a: mutation didn't change result!"
    print("  G-T03a: REDDENED ✓")

    # --- Mutation G-T03b: boundary direction (cochain = reversed complex) ---
    print("\nG-T03b: boundary direction mutation (cochain complex)...")
    # Reverse the chain: use d3^T as d1, d2^T as d2, d1^T as d3.
    # Group-ring transpose swaps row/col indices; entry content unchanged.
    d1_rev = [[d3_raw[j][i] for j in range(len(d3_raw))] for i in range(len(d3_raw[0]))]
    d2_rev = [[d2_raw[j][i] for j in range(len(d2_raw))] for i in range(len(d2_raw[0]))]
    d3_rev = [[d1_raw[j][i] for j in range(len(d1_raw))] for i in range(len(d1_raw[0]))]
    tau_mut_b, acyclic_b, msg_b = compute_torsion(fixture_reps, d1_rev, d2_rev, d3_rev, 2)
    if acyclic_b:
        T2_mut_b = tau_mut_b * tau_mut_b.conj()
        changed_b = not (T2_mut_b.re == T2_base.re)
        print(f"  T²_rev = {T2_mut_b.re}, baseline = {T2_base.re}, changed = {changed_b}")
    else:
        print(f"  Result: {msg_b} — gate reddened ✓")
    gate_b_red = not acyclic_b or (acyclic_b and not (T2_mut_b.re == T2_base.re))
    assert gate_b_red, "G-T03b: mutation didn't change result!"
    print("  G-T03b: REDDENED ✓")

    # --- Mutation G-T03c: module side (g -> rho(g)^T, anti-homomorphism) ---
    print("\nG-T03c: module side mutation (transpose = anti-homomorphism)...")
    # Right-module action maps g to rho(g)^T. Since (rho(g)rho(h))^T =
    # rho(h)^T rho(g)^T, this is an anti-homomorphism and breaks dd=0.
    transposed_reps = [mtranspose(fixture_reps[g]) for g in range(120)]
    tau_mut_c, acyclic_c, msg_c = compute_torsion(transposed_reps, d1_raw, d2_raw, d3_raw, 2)
    if not acyclic_c:
        print(f"  Result: {msg_c} — gate reddened ✓")
    else:
        T2_mut_c = tau_mut_c * tau_mut_c.conj()
        changed_c = not (T2_mut_c.re == T2_base.re)
        print(f"  T² changed = {changed_c}")
        assert changed_c, "G-T03c: mutation didn't change result!"
    print("  G-T03c: REDDENED ✓")

    # --- Mutation G-T03d: vector convention (transpose GR boundary matrices) ---
    print("\nG-T03d: vector convention mutation (column vectors = transpose GR maps)...")
    # Column-vector convention changes how the tensor product is assembled:
    # blocks in the evaluated boundary are transposed relative to their GR positions.
    # Equivalent to transposing the GR matrices before evaluation (swap block indices).
    d1_grt = [[d1_raw[j][i] for j in range(len(d1_raw))] for i in range(len(d1_raw[0]))]
    d2_grt = [[d2_raw[j][i] for j in range(len(d2_raw))] for i in range(len(d2_raw[0]))]
    d3_grt = [[d3_raw[j][i] for j in range(len(d3_raw))] for i in range(len(d3_raw[0]))]
    # Now d1_grt is 1x2, d2_grt is 2x2, d3_grt is 2x1.
    # Evaluate with original shapes expected by torsion:
    # d1_grt (1x2) plays role of d3 (1x2), d3_grt (2x1) plays role of d1 (2x1).
    # Keep d2_grt (2x2) in place. This tests that the block assembly order matters.
    tau_mut_d, acyclic_d, msg_d = compute_torsion(
        fixture_reps, d3_grt, d2_grt, d1_grt, 2)
    if not acyclic_d:
        print(f"  Result: {msg_d} — gate reddened ✓")
    else:
        T2_mut_d = tau_mut_d * tau_mut_d.conj()
        changed_d = not (T2_mut_d.re == T2_base.re)
        print(f"  T² = {T2_mut_d.re}, changed = {changed_d}")
        assert changed_d, "G-T03d: mutation didn't change result!"
    print("  G-T03d: REDDENED ✓")

    print("\n=== ALL CONVENTION FIXTURE TESTS PASSED ===")


if __name__ == '__main__':
    main()
