"""
Pre-implementation validation: saturation certificates for integral homology.
Proves H_*(C_*) = (Z, 0, 0, Z) integrally by finding unimodular minors.
"""
from fractions import Fraction
import json


class QGold:
    __slots__ = ('a', 'b')
    def __init__(self, a, b=0):
        self.a = Fraction(a); self.b = Fraction(b)
    def __eq__(self, other):
        if isinstance(other, (int, Fraction)): return self.a == other and self.b == 0
        return self.a == other.a and self.b == other.b
    def __hash__(self): return hash((self.a, self.b))
    def __add__(self, other):
        if isinstance(other, (int, Fraction)): return QGold(self.a + other, self.b)
        return QGold(self.a + other.a, self.b + other.b)
    def __radd__(self, other): return self.__add__(other)
    def __sub__(self, other):
        if isinstance(other, (int, Fraction)): return QGold(self.a - other, self.b)
        return QGold(self.a - other.a, self.b - other.b)
    def __neg__(self): return QGold(-self.a, -self.b)
    def __mul__(self, other):
        if isinstance(other, (int, Fraction)): return QGold(self.a * other, self.b * other)
        return QGold(self.a*other.a + self.b*other.b, self.a*other.b + self.b*other.a + self.b*other.b)
    def __rmul__(self, other): return self.__mul__(other)

class Quat:
    __slots__ = ('w','x','y','z')
    def __init__(self, w, x, y, z):
        self.w = w if isinstance(w, QGold) else QGold(w)
        self.x = x if isinstance(x, QGold) else QGold(x)
        self.y = y if isinstance(y, QGold) else QGold(y)
        self.z = z if isinstance(z, QGold) else QGold(z)
    def __eq__(self, other):
        return self.w==other.w and self.x==other.x and self.y==other.y and self.z==other.z
    def __hash__(self): return hash((self.w, self.x, self.y, self.z))
    def __mul__(self, other):
        w = self.w*other.w - self.x*other.x - self.y*other.y - self.z*other.z
        x = self.w*other.x + self.x*other.w + self.y*other.z - self.z*other.y
        y = self.w*other.y - self.x*other.z + self.y*other.w + self.z*other.x
        z = self.w*other.z + self.x*other.y - self.y*other.x + self.z*other.w
        return Quat(w, x, y, z)
    def __neg__(self): return Quat(-self.w, -self.x, -self.y, -self.z)
    def sort_key(self):
        key = []
        for c in [self.w, self.x, self.y, self.z]:
            A=c.a*2; B=c.b*2
            key.append(int(A)); key.append(int(B))
        return tuple(key)

def parse_gen(strs):
    comps = []
    for s in strs:
        s = s.strip(); inner = s[1:-3]
        parts = inner.replace(' ','').replace('*phi','p')
        tokens=[]; cur=''
        for ch in parts:
            if ch in '+-' and cur: tokens.append(cur); cur=ch
            else: cur+=ch
        if cur: tokens.append(cur)
        a=0; b=0
        for tok in tokens:
            if 'p' in tok:
                c=tok.replace('p','')
                b = 1 if c in ('','+') else (-1 if c=='-' else int(c))
            else: a=int(tok)
        comps.append(QGold(Fraction(a,2), Fraction(b,2)))
    return Quat(*comps)

def enum_group(gp):
    g1=parse_gen(gp['generators'][0]); g2=parse_gen(gp['generators'][1])
    e=Quat(QGold(1),QGold(0),QGold(0),QGold(0))
    elems={e,g1,g2}
    while True:
        new=set()
        for a in list(elems):
            for b in [g1,g2,-g1,-g2]:
                for c in [a*b, b*a]:
                    if c not in elems and c not in new: new.add(c)
        if not new: break
        elems.update(new)
    return sorted(elems, key=lambda q: q.sort_key())

def build_mult(se):
    n=len(se); e2r={q:i for i,q in enumerate(se)}
    return [[e2r[se[i]*se[j]] for j in range(n)] for i in range(n)]

def expand_gr_mat(mat, mult_table, n_group=120):
    rows_gr=len(mat); cols_gr=len(mat[0])
    inv=[0]*n_group
    for i in range(n_group):
        for j in range(n_group):
            if mult_table[i][j]==119:
                inv[i]=j; break
    Z_mat=[[0]*(cols_gr*n_group) for _ in range(rows_gr*n_group)]
    for bi in range(rows_gr):
        for bj in range(cols_gr):
            for eid,coeff in mat[bi][bj].items():
                for a in range(n_group):
                    b=mult_table[a][eid]
                    Z_mat[bi*n_group+a][bj*n_group+b]+=coeff
    return Z_mat

def parse_bmap(data):
    rows=[]
    for rd in data:
        row=[]
        for ed in rd:
            gr={}
            for c,eid in ed:
                gr[eid]=gr.get(eid,0)+c
                if gr[eid]==0: del gr[eid]
            row.append(gr)
        rows.append(row)
    return rows

def gauss_pivots(mat):
    m=len(mat)
    if m==0: return 0,[],[]
    n=len(mat[0])
    M=[[Fraction(mat[i][j]) for j in range(n)] for i in range(m)]
    rp=list(range(m)); pr=[]; pc=[]; r=0
    for col in range(n):
        piv=None
        for row in range(r,m):
            if M[row][col]!=0: piv=row; break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        rp[r],rp[piv]=rp[piv],rp[r]
        s=M[r][col]
        for j in range(n): M[r][j]/=s
        for row in range(m):
            if row==r: continue
            f=M[row][col]
            if f!=0:
                for j in range(n): M[row][j]-=f*M[r][j]
        pr.append(rp[r]); pc.append(col); r+=1
    return r,pr,pc

def det_frac(mat):
    n=len(mat)
    M=[[Fraction(mat[i][j]) for j in range(n)] for i in range(n)]
    d=Fraction(1)
    for col in range(n):
        piv=None
        for row in range(col,n):
            if M[row][col]!=0: piv=row; break
        if piv is None: return Fraction(0)
        if piv!=col: M[col],M[piv]=M[piv],M[col]; d=-d
        d*=M[col][col]; s=M[col][col]
        for j in range(col,n): M[col][j]/=s
        for row in range(col+1,n):
            f=M[row][col]
            if f!=0:
                for j in range(col,n): M[row][j]-=f*M[col][j]
    return d


def main():
    with open('m8_5a_packet.json') as f: gp=json.load(f)
    with open('m8_8_construction_packet.json') as f: cp=json.load(f)

    print("Setting up group...")
    se=enum_group(gp); mt=build_mult(se)
    d1=parse_bmap(cp['boundary_maps']['d1'])
    d2=parse_bmap(cp['boundary_maps']['d2'])
    d3=parse_bmap(cp['boundary_maps']['d3'])

    print("Expanding over Z...")
    d1Z=expand_gr_mat(d1,mt); d2Z=expand_gr_mat(d2,mt); d3Z=expand_gr_mat(d3,mt)
    print(f"d1_Z: {len(d1Z)}x{len(d1Z[0])}, d2_Z: {len(d2Z)}x{len(d2Z[0])}, d3_Z: {len(d3Z)}x{len(d3Z[0])}")

    # d3_Z: 120 x 240
    print("\n--- im(d3) saturation ---")
    r3,pr3,pc3=gauss_pivots(d3Z)
    print(f"rank(d3)={r3}")
    assert r3==119
    sub3=[[d3Z[pr3[i]][pc3[j]] for j in range(119)] for i in range(119)]
    det3=det_frac(sub3)
    print(f"|det|={abs(det3)}")
    assert abs(det3)==1
    print("SATURATED ✓")

    # d2_Z: 240 x 240
    print("\n--- im(d2) saturation ---")
    r2,pr2,pc2=gauss_pivots(d2Z)
    print(f"rank(d2)={r2}")
    assert r2==121
    sub2=[[d2Z[pr2[i]][pc2[j]] for j in range(121)] for i in range(121)]
    det2=det_frac(sub2)
    print(f"|det|={abs(det2)}")
    if abs(det2)==1:
        print("SATURATED ✓")
    else:
        print(f"WARNING: det={det2}")

    # d1_Z: 240 x 120
    print("\n--- im(d1) saturation ---")
    r1,pr1,pc1=gauss_pivots(d1Z)
    print(f"rank(d1)={r1}")
    assert r1==119
    sub1=[[d1Z[pr1[i]][pc1[j]] for j in range(119)] for i in range(119)]
    det1=det_frac(sub1)
    print(f"|det|={abs(det1)}")
    if abs(det1)==1:
        print("SATURATED ✓")
    else:
        print(f"WARNING: det={det1}")

    print("\n=== SATURATION SUMMARY ===")
    all_ok = abs(det1)==1 and abs(det2)==1 and abs(det3)==1
    for name,det in [("im(d3)",det3),("im(d2)",det2),("im(d1)",det1)]:
        status = "SATURATED" if abs(det)==1 else f"INDEX={abs(det)}"
        print(f"  {name}: {status}")

    if all_ok:
        print("\nAll boundary images are saturated (direct summands).")
        print("Combined with correct rational ranks (119, 121, 119):")
        print("  H_0 = Z^120/im(d1) free, rank 1 => Z")
        print("  H_1 = ker(d1)/im(d2) free, rank 0 => 0")
        print("  H_2 = ker(d2)/im(d3) free, rank 0 => 0")
        print("  H_3 = ker(d3) free, rank 1 => Z")
        print("\n  H_*(C_*) = (Z, 0, 0, Z): ESTABLISHED ✓")
    else:
        print("\nSome boundary images are NOT saturated. Need further analysis.")


if __name__=='__main__':
    main()
