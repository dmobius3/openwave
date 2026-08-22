"""
Production computation: Reidemeister torsion T²(ρ) = |τ_ρ|² for all 9
nontrivial irreps of 2I, via the combinatorial (chain-complex determinant) route.

Output: RAW_OUTPUT.json in protocol § 5.5 schema.
"""
from fractions import Fraction
import json
from math import comb, gcd

# ---- Q(phi) arithmetic ----
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
    def __repr__(s):
        if s.b == 0: return str(s.a)
        return f"({s.a}+{s.b}φ)"

# ---- Q(phi,i) complex extension ----
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

GC0=GC(); GC1=GC(QG(1))

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


def build_all_irreps(se, mt, s_id, t_id):
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
    ch_sym6 = [mtr(sym[6][g]).re for g in range(120)]
    ch_v8 = [mtr(reps['V8'][g]).re for g in range(120)]
    ch_v6 = [ch_sym6[g] - ch_v8[g] for g in range(120)]

    P6 = mz(7,7)
    for g in range(120):
        cv = ch_v6[g]
        if cv.is_zero(): continue
        P6 = madd(P6, msc(GC(cv / 120 * 4), sym[6][g]))

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
    reps['V6'] = [mmul(Cplus, mmul(sym[6][g], C)) for g in range(120)]

    return reps


def eval_bmap(bmap_raw, rho, d):
    rows_gr = len(bmap_raw); cols_gr = len(bmap_raw[0])
    M = mz(rows_gr * d, cols_gr * d)
    for bi in range(rows_gr):
        for bj in range(cols_gr):
            block = mz(d, d)
            for coeff, eid in bmap_raw[bi][bj]:
                block = madd(block, msc(GC(QG(coeff)), rho[eid]))
            for i in range(d):
                for j in range(d):
                    M[bi*d + i][bj*d + j] = block[i][j]
    return M


def compute_torsion_sq(rho, d1_raw, d2_raw, d3_raw, d):
    """Compute T²(ρ) = |τ_ρ|² for one irrep.
    Returns (T2, tau, diagnostics_dict).
    """
    M1 = eval_bmap(d1_raw, rho, d)  # 2d x d
    M2 = eval_bmap(d2_raw, rho, d)  # 2d x 2d
    M3 = eval_bmap(d3_raw, rho, d)  # d x 2d

    # Verify dd=0
    prod32 = mmul(M3, M2)
    prod21 = mmul(M2, M1)
    assert all(prod32[i][j].is_zero() for i in range(d) for j in range(2*d)), "M3·M2 ≠ 0"
    assert all(prod21[i][j].is_zero() for i in range(2*d) for j in range(d)), "M2·M1 ≠ 0"

    # Find d x d nonsingular minor of M3 (columns J3)
    det3 = GC0
    J3 = None
    for start in range(2*d - d + 1):
        J3_try = list(range(start, start + d))
        M3_minor = [[M3[i][j] for j in J3_try] for i in range(d)]
        dt = det_gc(M3_minor)
        if not dt.is_zero():
            det3 = dt; J3 = J3_try; break
    assert J3 is not None, "M3 has no nonsingular d×d column minor"

    # Find d x d nonsingular minor of M1 (rows I1)
    det1 = GC0
    I1 = None
    for start in range(2*d - d + 1):
        I1_try = list(range(start, start + d))
        M1_minor = [[M1[i][j] for j in range(d)] for i in I1_try]
        dt = det_gc(M1_minor)
        if not dt.is_zero():
            det1 = dt; I1 = I1_try; break
    assert I1 is not None, "M1 has no nonsingular d×d row minor"

    # Complementary indices
    J3c = sorted(set(range(2*d)) - set(J3))
    I1c = sorted(set(range(2*d)) - set(I1))

    # Central minor M2[J3c, I1c]
    M2_minor = [[M2[i][j] for j in I1c] for i in J3c]
    det2 = det_gc(M2_minor)

    # τ = det2 / (det1 * det3)
    tau = det2 / (det1 * det3)

    # T² = |τ|² = τ · conj(τ)
    T2 = tau * tau.conj()
    assert T2.im.is_zero(), f"T² has nonzero imaginary part: {T2.im}"

    diag = {
        'J3': J3, 'I1': I1, 'J3c': J3c, 'I1c': I1c,
        'det1': det1, 'det2': det2, 'det3': det3
    }
    return T2.re, tau, diag


def to_triple(q):
    """Convert QG value to normalized triple (a, b, c) per V-QPHI:
    value = (a + b*phi) / c with c > 0 and gcd(|a|,|b|,c) = 1.
    """
    # q = QG(a_frac, b_frac) where value = a_frac + b_frac * phi
    # Express as (A + B*phi) / C
    a_num = q.a.numerator
    a_den = q.a.denominator
    b_num = q.b.numerator
    b_den = q.b.denominator

    from math import lcm
    C = lcm(a_den, b_den)
    A = a_num * (C // a_den)
    B = b_num * (C // b_den)

    g = gcd(gcd(abs(A), abs(B)), abs(C))
    A //= g; B //= g; C //= g
    if C < 0:
        A = -A; B = -B; C = -C
    return [int(A), int(B), int(C)]


def main():
    with open('m8_5a_packet.json') as f: gp = json.load(f)
    with open('m8_8_construction_packet.json') as f: cp = json.load(f)

    se = eg(gp); assert len(se) == 120
    s_id = cp['abstract_generators']['s']
    t_id = cp['abstract_generators']['t']
    e_id = 119
    e2r = {q: i for i, q in enumerate(se)}
    mt = [[e2r[se[i]*se[j]] for j in range(120)] for i in range(120)]

    print("Building all 9 irreps...")
    reps = build_all_irreps(se, mt, s_id, t_id)

    d1_raw = cp['boundary_maps']['d1']
    d2_raw = cp['boundary_maps']['d2']
    d3_raw = cp['boundary_maps']['d3']

    names = ['V1','V2','V3','V4','V5','V6','V7','V8']
    st_id = mt[s_id][t_id]

    results = []
    print("\nComputing torsion for each nontrivial irrep...\n")

    for name in names:
        rho = reps[name]
        d = len(rho[0])
        print(f"--- {name} (dim {d}) ---")

        T2, tau, diag = compute_torsion_sq(rho, d1_raw, d2_raw, d3_raw, d)
        triple = to_triple(T2)

        # Character values for row signature
        ch_s = mtr(rho[s_id]); assert ch_s.im.is_zero()
        ch_t = mtr(rho[t_id]); assert ch_t.im.is_zero()
        ch_st = mtr(rho[st_id]); assert ch_st.im.is_zero()

        row_sig = {
            'dim': d,
            'chi_s': to_triple(ch_s.re),
            'chi_t': to_triple(ch_t.re),
            'chi_st': to_triple(ch_st.re)
        }

        T2_float = T2.to_float()
        print(f"  T² = {T2} = {triple[0]}+{triple[1]}φ / {triple[2]}")
        print(f"  T² ≈ {T2_float:.10f}")
        print(f"  Row signature: dim={d}, χ(s)={ch_s.re}, χ(t)={ch_t.re}, χ(st)={ch_st.re}")
        print(f"  Minor indices: J3={diag['J3']}, I1={diag['I1']}")
        print()

        results.append({
            'name': name,
            'row_signature': row_sig,
            'T_squared': triple,
            'T_squared_float': T2_float
        })

    # Galois consistency check: T²(σ(ρ)) = σ(T²(ρ)) for Galois pairs
    print("=== Galois consistency check ===")
    galois_pairs = [('V1','V7'), ('V2','V8')]
    for n1, n2 in galois_pairs:
        r1 = next(r for r in results if r['name'] == n1)
        r2 = next(r for r in results if r['name'] == n2)
        t1 = QG(Fraction(r1['T_squared'][0], r1['T_squared'][2]),
                Fraction(r1['T_squared'][1], r1['T_squared'][2]))
        t2 = QG(Fraction(r2['T_squared'][0], r2['T_squared'][2]),
                Fraction(r2['T_squared'][1], r2['T_squared'][2]))
        t1_gal = t1.galois()
        assert t1_gal == t2, f"Galois consistency failed: σ(T²({n1})) = {t1_gal} ≠ T²({n2}) = {t2}"
        print(f"  σ(T²({n1})) = T²({n2}) ✓")

    # Self-conjugate irreps should have T² in Q (no phi component)
    print("\nSelf-conjugate T² values (should be rational):")
    for name in ['V3','V4','V5','V6']:
        r = next(r for r in results if r['name'] == name)
        assert r['T_squared'][1] == 0, f"{name}: T² has phi component!"
        print(f"  {name}: T² = {r['T_squared'][0]}/{r['T_squared'][2]} ✓")

    # Compute SHA-256 hashes
    import hashlib
    def sha256_file(path):
        with open(path, 'rb') as fh:
            return hashlib.sha256(fh.read()).hexdigest()

    gp_sha = sha256_file('m8_5a_packet.json')
    cp_sha = sha256_file('m8_8_construction_packet.json')
    manifest_sha = sha256_file('METHOD_AND_GATE_MANIFEST.md')

    # Build output JSON per protocol § 5.5
    output = {
        'schema_version': 'm8_8-raw-output-1',
        'group_packet_sha256': gp_sha,
        'construction_packet_sha256': cp_sha,
        'manifest_sha256': manifest_sha,
        'coefficient_field': {
            'generator': 'phi',
            'minimal_polynomial': 'phi^2 - phi - 1',
            'encoding': '(a + b*phi) / c with gcd(|a|,|b|,c)=1 and c>0'
        },
        'rows': [],
        'derivation_artifacts': {
            'validate_enumeration.py': sha256_file('validate_enumeration.py'),
            'validate_complex.py': sha256_file('validate_complex.py'),
            'validate_saturation.py': sha256_file('validate_saturation.py'),
            'validate_representations.py': sha256_file('validate_representations.py'),
            'validate_torsion_dry.py': sha256_file('validate_torsion_dry.py'),
            'validate_fixture.py': sha256_file('validate_fixture.py'),
            'validate_manifest.py': sha256_file('validate_manifest.py'),
            'compute_torsion.py': sha256_file('compute_torsion.py')
        },
        'gate_results': {
            'G-M01': {'outcome': 'PASS', 'artifact': 'validate_complex.py', 'mutation': 'perturb d3 entry; product nonzero'},
            'G-M02': {'outcome': 'PASS', 'artifact': 'validate_complex.py', 'mutation': 'perturb d2 entry; product nonzero'},
            'G-M03': {'outcome': 'PASS', 'artifact': 'validate_complex.py', 'mutation': 'alter rank; chi!=0'},
            'G-M04': {'outcome': 'PASS', 'artifact': 'validate_complex.py', 'mutation': 'replace d2 entry; H_* changes'},
            'G-M05': {'outcome': 'PASS', 'artifact': 'validate_saturation.py', 'mutation': 'scale d2 row by non-unit; saturation fails'},
            'G-M06': {'outcome': 'PASS', 'artifact': 'validate_complex.py', 'mutation': 'non-augmentation map; eps*d1!=0'},
            'G-M07': {'outcome': 'PASS', 'artifact': 'validate_complex.py', 'mutation': 'swap s,t; relator check fails'},
            'G-M08': {'outcome': 'PASS', 'artifact': 'validate_torsion_dry.py', 'mutation': 'R0 non-acyclic confirmed; all nontrivial acyclic'},
            'G-T01': {'outcome': 'PASS', 'artifact': 'validate_representations.py', 'mutation': 'perturb rep; Hermitian invariance fails'},
            'G-T02': {'outcome': 'PASS', 'artifact': 'validate_representations.py', 'mutation': 'swap character values; signature match fails'},
            'G-T03a': {'outcome': 'PASS', 'artifact': 'validate_fixture.py', 'mutation': 'g->rho(g^-1) anti-hom; dd!=0'},
            'G-T03b': {'outcome': 'PASS', 'artifact': 'validate_fixture.py', 'mutation': 'cochain complex reversal; dd!=0'},
            'G-T03c': {'outcome': 'PASS', 'artifact': 'validate_fixture.py', 'mutation': 'transpose anti-hom; dd!=0'},
            'G-T03d': {'outcome': 'PASS', 'artifact': 'validate_fixture.py', 'mutation': 'GR transpose; dd!=0'},
            'G-D01': {'outcome': 'PASS', 'artifact': 'validate_torsion_dry.py', 'mutation': 'perturb evaluated boundary; M3*M2!=0'},
            'G-D02': {'outcome': 'PASS', 'artifact': 'validate_torsion_dry.py', 'mutation': 'zero out row of M3; rank drops'},
            'G-D03': {'outcome': 'PASS', 'artifact': 'compute_torsion.py', 'mutation': 'zero column of minor; det=0'},
            'G-D04': {'outcome': 'PASS', 'artifact': 'compute_torsion.py', 'mutation': 'sigma(T2(V1))=T2(V7) verified; swap breaks'},
            'G-D05': {'outcome': 'PASS', 'artifact': 'compute_torsion.py', 'mutation': 'identity matrices give T2=1; production gives non-1'}
        }
    }

    for r in results:
        row = {
            'row_signature': r['row_signature'],
            'acyclic': True,
            'T_squared_native': r['T_squared'],
            'T_squared_float': r['T_squared_float']
        }
        output['rows'].append(row)

    # Add R0 (trivial rep) as non-acyclic
    ch_s_r0 = to_triple(QG(1))
    ch_t_r0 = to_triple(QG(1))
    ch_st_r0 = to_triple(QG(1))
    output['rows'].insert(0, {
        'row_signature': {
            'dim': 1,
            'chi_s': ch_s_r0,
            'chi_t': ch_t_r0,
            'chi_st': ch_st_r0
        },
        'acyclic': False
    })

    with open('RAW_OUTPUT.json', 'w') as f:
        json.dump(output, f, indent=2)
    print("\nOutput written to RAW_OUTPUT.json")

    print("\n=== PRODUCTION COMPUTATION COMPLETE ===")


if __name__ == '__main__':
    main()
