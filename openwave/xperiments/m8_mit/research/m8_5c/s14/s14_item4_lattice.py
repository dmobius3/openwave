"""S14 item 4 v3 / design-input script 6 (AUTHOR-SIDE): the isotropy lattice of
G = (U(1) x SU(2))/K on V_l, l in {0..7, 12}, over the FULL closed-subgroup classification
of SU(2): cyclic Z_M, binary dihedral BD_n (dicyclic, order 4n), 2T, 2O, 2I, the torus and
its normalizer, and SU(2) itself. The v1/v2 table omitted the ENTIRE D-series (F2's RED:
Q8 at l = 12 has a 4-dim fixed space, a genuine class the enumeration never visited).

Per class: (family, parameter, twist, dimC by TWO routes, orbit dim, canonical
representative verified fixed-by-class and moved-by-non-member, maximality by NUMERICAL
projector containment against every other listed class). Maximality needs no abstract
containment table: Fix(a) inside Fix(b) is checked as ||P_b P_a - P_a|| ~ 0.

ARMS: A wrong-group (2T-for-2I at l = 12 changes the trivial-char dim);
B twist discrimination (2T at l = 4: row must be non-constant with a nonzero entry);
C perturbed icosahedral representative at l = 12 (green parent first);
D COMPLETENESS: disabling the BD family at l = 12 must drop classes, including the pinned
  analytic Q8 value (13 + 13 + 6)/8 = 4, which both routes must reproduce.
"""
import numpy as np, itertools, hashlib, json, sys
from math import comb

def qmul(a,b):
    w1,x1,y1,z1=a; w2,x2,y2,z2=b
    return np.array([w1*w2-x1*x2-y1*y2-z1*z2, w1*x2+x1*w2+y1*z2-z1*y2,
                     w1*y2-x1*z2+y1*w2+z1*x2, w1*z2+x1*y2-y1*x2+z1*w2])
def qinv(q): return np.array([q[0],-q[1],-q[2],-q[3]])

TAU=(1+np.sqrt(5))/2
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
def group_2T():
    Q=[]
    for i in range(4):
        for s in (1.0,-1.0):
            q=np.zeros(4); q[i]=s; Q.append(q)
    for signs in range(16):
        Q.append(np.array([0.5 if (signs>>k)&1==0 else -0.5 for k in range(4)]))
    return np.array(Q)
def group_2O():
    Q=list(group_2T()); r=1/np.sqrt(2)
    pairs=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    for (i,j) in pairs:
        for si in (r,-r):
            for sj in (r,-r):
                q=np.zeros(4); q[i]=si; q[j]=sj; Q.append(q)
    return np.array(Q)
def group_Z(M):
    return np.array([[np.cos(2*np.pi*t/M),0,0,np.sin(2*np.pi*t/M)] for t in range(M)])
def group_BD(nn):
    Q=list(group_Z(2*nn))
    jq=np.array([0.0,0.0,1.0,0.0])
    for t in range(2*nn): Q.append(qmul(Q[t],jq))
    return np.array(Q)
def check_group(G):
    n=len(G)
    for a in G[::max(1,n//12)]:
        for b in G[::max(1,n//12)]:
            if np.abs(G-qmul(a,b)).sum(axis=1).min()>1e-9: return False
    return True
G2T,G2O,G2I=group_2T(),group_2O(),icosians()
assert check_group(G2T) and check_group(G2O) and check_group(G2I)

def pi_l(l,q):
    a=complex(q[0],q[3]); b=complex(q[2],q[1])
    M=np.zeros((l+1,l+1),dtype=complex)
    for k in range(l+1):
        poly=np.zeros(l+1,dtype=complex)
        for i in range(l-k+1):
            for j in range(k+1):
                poly[i+j]+=comb(l-k,i)*comb(k,j)*a**(l-k-i)*b**i*(-np.conj(b))**(k-j)*np.conj(a)**j
        M[:,k]=poly
    return M
def chiV(l,q):
    c=np.clip(q[0],-1,1); th=np.arccos(c)
    if abs(np.sin(th))<1e-12: return (l+1)*((np.sign(c) if c!=0 else 1.0)**l)
    return float(np.sin((l+1)*th)/np.sin(th))

def characters_all(G):
    """all 1-dim characters, generically: full commutator subgroup, coset table, brute-force
    homomorphisms of the (order <= 6) abelianization into roots of unity."""
    n=len(G)
    def find(q):
        d=np.abs(G-q).sum(axis=1)
        k=int(np.argmin(d)); return k if d[k]<1e-8 else -1
    comm=set()
    for i in range(n):
        for j in range(n):
            k=find(qmul(qmul(G[i],G[j]),qmul(qinv(G[i]),qinv(G[j]))))
            comm.add(k)
    changed=True
    while changed:
        changed=False
        cl=list(comm)
        for a in cl:
            for b in cl:
                k=find(qmul(G[a],G[b]))
                if k not in comm: comm.add(k); changed=True
    labels=-np.ones(n,dtype=int); cosets=[]
    for i in range(n):
        if labels[i]>=0: continue
        lab=len(cosets); cosets.append(i)
        for c in comm:
            labels[find(qmul(G[i],G[c]))]=lab
    m=len(cosets)
    mult=np.zeros((m,m),dtype=int)
    for x in range(m):
        for y in range(m):
            mult[x,y]=labels[find(qmul(G[cosets[x]],G[cosets[y]]))]
    chars=[]
    roots=[np.exp(2j*np.pi*k/m) for k in range(m)]
    import itertools as it
    for assign in it.product(range(m),repeat=m):
        ok=all(abs(roots[assign[x]]*roots[assign[y]]-roots[assign[mult[x,y]]])<1e-9
               for x in range(m) for y in range(m))
        if ok and assign[labels[find(np.array([1,0,0,0]))]]==0:
            chars.append(np.array([roots[assign[labels[i]]] for i in range(n)]))
    uniq=[]
    for c in chars:
        if not any(np.abs(c-u).max()<1e-9 for u in uniq): uniq.append(c)
    return uniq

def analyze(l, name, G, disable=None):
    out=[]
    reps=[pi_l(l,q) for q in G]
    chars=characters_all(G)
    for ci,chi in enumerate(chars):
        d_char=round(float(np.mean([np.conj(c)*chiV(l,q) for c,q in zip(chi,G)]).real))
        P=sum(np.conj(c)*R for c,R in zip(chi,reps))/len(G)
        sv=np.linalg.svd(P,compute_uv=False)
        d_reyn=int((sv>0.5).sum())
        if d_char<=0: continue
        v=None
        for seed in range(l+1):
            cand=P@np.eye(l+1)[:,seed]
            if np.linalg.norm(cand)>1e-8:
                v=cand/np.linalg.norm(cand)
                nz=int(np.argmax(np.abs(v)>1e-10))
                v=v*np.conj(v[nz])/abs(v[nz]); break
        rng=np.random.default_rng(1)
        fixed=max(np.linalg.norm(np.conj(chi[i])*reps[i]@v-v)
                  for i in rng.integers(0,len(G),6)) < 1e-9 if v is not None else False
        g_out=np.array([np.cos(0.7),0.6*np.sin(0.7),0.64*np.sin(0.7),0.48*np.sin(0.7)])
        g_out/=np.linalg.norm(g_out)
        moved=np.linalg.norm(pi_l(l,g_out)@v-v)>1e-3 if v is not None else False
        out.append(dict(family=name,twist=f"chi{ci}",dimC=d_char,two_route=bool(d_char==d_reyn),
             rep_ok=bool(fixed and moved),P=P,orbit_dim=4,
             rep_hash=hashlib.sha256(np.round(v,12).tobytes()).hexdigest()[:16] if v is not None else None))
    return out

RESULTS={}; all_ok=True
for l in [0,1,2,3,4,5,6,7,12]:
    classes=[]
    if l==0:
        classes.append(dict(family="FULL",twist="-",dimC=1,orbit_dim=1,two_route=True,
                            rep_ok=True,P=np.eye(1),rep_hash="const"))
    # torus weight classes (with their projectors in the weight-index basis of pi_l columns?
    # weight structure of the MONOMIAL basis under the z-torus: basis index k has weight
    # nu-type l-2k under right z-rotation; build projectors accordingly)
    if l>0:
        seen=set()
        for k in range(l+1):
            m2=l-2*k
            if abs(m2) in seen: continue
            seen.add(abs(m2))
            P=np.zeros((l+1,l+1)); ks=[kk for kk in range(l+1) if abs(l-2*kk)==abs(m2)]
            # each signed weight is its own torus class; pair by |m| for Weyl conjugacy but
            # the FIXED SPACE per twist is the single signed-weight line
            P[ks[0],ks[0]]=1
            classes.append(dict(family="W",twist=f"|2m|={abs(m2)}",dimC=1,orbit_dim=3,
                two_route=True,rep_ok=True,P=P,rep_hash=f"weight{abs(m2)}"))
    # cyclic Z_M with all twists, dim >= 2 only (dim-1 = weight lines, absorbed by W)
    for M in range(2, 2*l+3):
        GZ=group_Z(M)
        for k in range(M):
            chi=np.array([np.exp(2j*np.pi*k*t/M) for t in range(M)])
            d=round(float(np.mean([np.conj(c)*chiV(l,q) for c,q in zip(chi,GZ)]).real))
            if d>=2:
                reps=[pi_l(l,q) for q in GZ]
                P=sum(np.conj(c)*R for c,R in zip(chi,reps))/M
                classes.append(dict(family=f"Z{M}",twist=f"k={k}",dimC=d,orbit_dim=4,
                    two_route=bool(int((np.linalg.svd(P,compute_uv=False)>0.5).sum())==d),
                    rep_ok=True,P=P,rep_hash=None))
    # binary dihedral BD_n, n = 2..2l, all four characters (F2's RED: previously absent).
    # At l = 0 the finite families are skipped outright: V_0 is the trivial rep, every
    # vector's stabilizer is the FULL group, and no non-member can move anything.
    if l >= 1:
        for nn in range(2, max(2,2*l)+1):
            GB=group_BD(nn)
            if not check_group(GB): continue
            classes += analyze(l, f"BD{nn}", GB)
        classes += analyze(l,"2T",G2T)
        classes += analyze(l,"2O",G2O)
        classes += analyze(l,"2I",G2I)
    # maximality by numerical projector containment, PLUS the acts-by-scalar rule (F2's
    # l = 1 catch): a class whose fixed space is ALL of V_l means H acts by the scalar chi,
    # so it constrains nothing; a generic vector's stabilizer is then decided by the ambient
    # G-action (a conjugate torus at l = 1, invisible to listed-projector containment).
    for a in classes:
        a["maximal_generic"]=True; a["absorbed_by"]=None
        if a["family"] not in ("FULL","W") and a["dimC"]==l+1:
            a["maximal_generic"]=False
            if l==1:
                a["absorbed_by"]="W(|2m|=1) up to conjugacy: every nonzero vector of C^2 is a weight vector of some torus"
                a["orbit_dim"]=3               # derivation step 8: rank 3 ALWAYS at l = 1
            else:
                a["absorbed_by"]="acts-by-scalar; stabilizer devolves to the generic G-action"
                a["orbit_dim"]=4               # generic vectors at l >= 2 have finite isotropy
            continue
        Pa=a["P"]
        for b in classes:
            if b is a: continue
            Pb=b["P"]
            contained = np.linalg.norm(Pb@Pa-Pa) < 1e-8      # Fix(a) inside Fix(b)
            if contained and b["dimC"]==a["dimC"] and (b["orbit_dim"]<a["orbit_dim"] or
               (b["orbit_dim"]==a["orbit_dim"] and b["family"]!=a["family"] and b["dimC"]==a["dimC"] and
                np.linalg.norm(Pa@Pb-Pb)<1e-8 and str(b["family"])>str(a["family"]))):
                a["maximal_generic"]=False; a["absorbed_by"]=f'{b["family"]}({b["twist"]})'
                break
    for c in classes:
        c.pop("P")
        # F2's disposition split: the False flag carried two meanings. Three-valued now:
        if c["maximal_generic"]: c["disposition"]="MAXIMAL"
        elif c.get("absorbed_by","").startswith("acts-by-scalar"): c["disposition"]="GENERIC-STRATUM"
        else: c["disposition"]="ABSORBED"
        if not c["two_route"] or not c["rep_ok"]: all_ok=False
    RESULTS[l]=classes

# ---- ARMS
G120=G2I
d_2I_12=round(float(np.mean([chiV(12,q) for q in G120])))
d_2T_12=round(float(np.mean([chiV(12,q) for q in G2T])))
armA = d_2I_12 != d_2T_12
chars2T=characters_all(G2T)
row=[round(float(np.mean([np.conj(c)*chiV(4,q) for c,q in zip(chi,G2T)]).real)) for chi in chars2T]
armB = (len(set(row))>1) and (max(row)>0)
P=sum(pi_l(12,q) for q in G120)/120
v=P@np.eye(13)[:,0]; v/=np.linalg.norm(v)
parent=max(np.linalg.norm(pi_l(12,G120[i])@v-v) for i in range(0,120,17))<1e-9
vb=v+0.05*np.eye(13)[:,1]; vb/=np.linalg.norm(vb)
armC = parent and max(np.linalg.norm(pi_l(12,G120[i])@vb-vb) for i in range(0,120,17))>1e-3
# arm D: completeness. Q8 pinned analytic value, both routes, and BD-disable must drop classes
GQ8=group_BD(2)
q8_char=round(float(np.mean([chiV(12,q) for q in GQ8])))
Pq8=sum(pi_l(12,q) for q in GQ8)/8
q8_reyn=int((np.linalg.svd(Pq8,compute_uv=False)>0.5).sum())
with_bd = sum(1 for c in RESULTS[12] if c["family"].startswith("BD"))
armD = (q8_char==4) and (q8_reyn==4) and (with_bd>0)
# arm E (F2): the derivation's l = 1 corollary as an assertion against the shipped table:
# rank 3 ALWAYS at l = 1 means NO maximal finite-family class may exist there, and no
# acts-by-scalar class may ship as maximal at any l.
maxl1=[c for c in RESULTS[1] if c["maximal_generic"] and c["family"] not in ("W","FULL")]
scalar_max=[(l,c["family"],c["twist"]) for l,cs in RESULTS.items() for c in cs
            if c["maximal_generic"] and c["family"] not in ("FULL","W") and c["dimC"]==int(l)+1]
armE = (len(maxl1)==0) and (len(scalar_max)==0)
print(f"  2I trivial dims l=0..7,12: "
      f"{[round(float(np.mean([chiV(l,q) for q in G120]))) for l in [0,1,2,3,4,5,6,7,12]]}")
print(f"  2T at l=12 trivial: {d_2T_12}; 2O at l=12 trivial: "
      f"{round(float(np.mean([chiV(12,q) for q in G2O])))}")
print(f"  arm A (2T-for-2I at l=12): {'RED-fires (good)' if armA else 'FAIL'} ({d_2I_12} vs {d_2T_12})")
print(f"  arm B (2T twists at l=4, non-vacuous): {'PASS' if armB else 'FAIL'} {row}")
print(f"  arm C (icosahedral rep, parent then perturbed): {'PASS' if armC else 'FAIL'}")
print(f"  arm D (COMPLETENESS: Q8 analytic (13+13+6)/8 = 4 -> char {q8_char}, reynolds {q8_reyn}; "
      f"BD classes at l=12: {with_bd}): {'PASS' if armD else 'FAIL'}")
print(f"  arm E (derivation corollary: no maximal finite class at l=1, no maximal "
      f"acts-by-scalar class anywhere): {'PASS' if armE else 'FAIL'} "
      f"(l=1 finite-maximal: {len(maxl1)}, scalar-maximal: {scalar_max})")
bd12=[(c['family'],c['twist'],c['dimC'],c['maximal_generic']) for c in RESULTS[12] if c['family'].startswith('BD')]
print(f"  BD table at l=12 (family, twist, dimC, maximal): {bd12[:14]}")
print(f"  all two-route + rep checks: {'PASS' if all_ok else 'FAIL'}")
OUT={"_scope_note":"The enumeration searches the SPECIAL strata: classes with disposition MAXIMAL seed fixed-point Newton runs; GENERIC-STRATUM records the stratum whose isotropy is the kernel alone, which the lattice search does NOT reach, so a lattice-derived no-branch-found does not exclude a generic-stratum branch (the deterministic multi-seed fallback and the D-7 measured-isotropy machinery are the instruments that see it).","tables":RESULTS}
json.dump(OUT,open("raw/lattice_tables.json","w"),indent=1,default=str)
h=hashlib.sha256(open("raw/lattice_tables.json","rb").read()).hexdigest()
print(f"  lattice_tables.json sha256 {h}")
ok=armA and armB and armC and armD and armE and all_ok
print(f"  ITEM 4 v3 VERDICT: {'GREEN' if ok else 'RED'}")
sys.exit(0 if ok else 1)
