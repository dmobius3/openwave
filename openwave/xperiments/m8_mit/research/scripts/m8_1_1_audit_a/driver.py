"""Driver: full independent recomputation of SPEC SHEET A + widened search."""
from fractions import Fraction
import json, sys, time
from part2 import (rref, nullspace, rank_gf, mat_mul, mat_eye, mat_add, mat_scal,
                   mat_trace, Cyc)
from part1 import p, zeta, mmul, EYE, quat, det, tr, build_groups
import core
from core import (build, chi_sym_exact, ip_exact, isotypic_split, lift_character,
                  sym_power, restrict, commutant_dim, mat_inv)

AMAX = 32          # symmetric powers V_0..V_32   (agent A used 14)
MMAX = 30          # E_m for m = 2..30            (agent A used 12)

def cyc_eq(x,y): return x.reduced()==y.reduced()

import cmath
def complex_val(z,M):
    v=sum(c*cmath.exp(2j*cmath.pi*l/M) for l,c in enumerate(z.c))
    return [round(v.real,15),round(v.imag,15)]

def analyse(name,gens,expected,verbose=True):
    t0=time.time()
    G=build(name,gens,expected)
    M=G.M
    # ---- exact symmetric-power characters
    chi=[chi_sym_exact(G,a) for a in range(AMAX+1)]
    # ---- irreducibles by isotypic splitting of V_a (exact, GF(p)), a ascending
    irr=[]            # list of dicts: char (list of Cyc), dim
    iso_mult={}       # (a, sigma_index) -> multiplicity from the module decomposition
    a=0
    while sum(x["dim"]**2 for x in irr)<G.N and a<=60:
        pieces,rho=isotypic_split(G,a)
        for pc in pieces:
            m=pc["mult"]; d=pc["dim"]; invm=pow(m,p-2,p)
            th=[mat_trace(pc["Rall"][G.crep[ci]],p)*invm%p for ci in range(G.nc)]
            ch=lift_character(G,th,d)
            k=None
            for si,x in enumerate(irr):
                if all(cyc_eq(ch[i],x["char"][i]) for i in range(G.nc)): k=si; break
            if k is None:
                irr.append({"char":ch,"dim":d}); k=len(irr)-1
            iso_mult[(a,k)]=m
        a+=1
    G.A0=a-1
    # ---- canonical order: trivial first, then (dim, exact-char key)
    triv=[i for i,x in enumerate(irr) if all(cyc_eq(x["char"][j],Cyc.mono(M,0)) for j in range(G.nc))]
    assert len(triv)==1,(name,triv)
    def key(x):
        return (x["dim"],)+tuple(tuple(map(str,x["char"][i].reduced())) for i in range(G.nc))
    rest=sorted([x for i,x in enumerate(irr) if i!=triv[0]],key=key)
    old=irr; irr=[old[triv[0]]]+rest
    remap={}
    for ni,x in enumerate(irr):
        for oi,y in enumerate(old):
            if x is y: remap[oi]=ni
    iso_mult={(aa,remap[ss]):v for (aa,ss),v in iso_mult.items()}
    nc=G.nc; assert len(irr)==nc,(name,len(irr),nc)
    dims=[x["dim"] for x in irr]
    sumsq=sum(d*d for d in dims)
    # ---- exact orthogonality
    orth_ok=True
    for i in range(nc):
        for j in range(nc):
            v=ip_exact(G,irr[i]["char"],irr[j]["char"])
            if v!=(1 if i==j else 0): orth_ok=False
    col_ok=True
    for i in range(nc):
        for j in range(nc):
            acc=Cyc(M)
            for sgm in range(nc): acc=acc+irr[sgm]["char"][i]*irr[sgm]["char"][j].conj()
            v=acc.to_rational()
            tgt=Fraction(G.N,G.csize[i]) if i==j else Fraction(0)
            if v!=tgt: col_ok=False
    # ---- branching, exact, a = 0..AMAX
    Mult=[[int(ip_exact(G,chi[aa],irr[s]["char"])) for s in range(nc)] for aa in range(AMAX+1)]
    iso_agree=all(Mult[aa][s]==v for (aa,s),v in iso_mult.items())
    dimchk=all(sum(Mult[aa][s]*dims[s] for s in range(nc))==aa+1 for aa in range(AMAX+1))
    # ---- adjacency, exact
    chi1=chi[1]
    Adj=[[int(ip_exact(G,[chi1[i]*irr[s]["char"][i] for i in range(nc)],irr[s2]["char"]))
          for s2 in range(nc)] for s in range(nc)]
    sym=all(Adj[i][j]==Adj[j][i] for i in range(nc) for j in range(nc))
    rowok=all(sum(Adj[s][s2]*dims[s2] for s2 in range(nc))==2*dims[s] for s in range(nc))
    # ---- distances / diameter
    def bfs(src):
        d=[-1]*nc; d[src]=0; fr=[src]
        while fr:
            nx=[]
            for u in fr:
                for v in range(nc):
                    if Adj[u][v] and d[v]<0: d[v]=d[u]+1; nx.append(v)
            fr=nx
        return d
    dist=bfs(0); diam=max(max(bfs(s)) for s in range(nc))
    # ---- T1 / T2
    T1=[]
    for s in range(nc):
        occ=[(aa,Mult[aa][s]) for aa in range(AMAX+1) if Mult[aa][s]]
        T1.append({"sigma":s,"dim":dims[s],"d":dist[s],
                   "least_a":occ[0][0] if occ else None,
                   "occurrences":[{"a":aa,"mult":mm} for aa,mm in occ]})
    viol=[{"sigma":s,"a":aa} for s in range(nc) for aa in range(AMAX+1)
          if Mult[aa][s] and (aa-dist[s])%2]
    # ---- mu tables
    def chiE(m_):
        return [chi[m_][i].scal(m_-1)+chi[m_-2][i].scal(m_+1) for i in range(nc)]
    taus=[]
    for s in range(nc):
        taus.append(("sigma%d"%s,irr[s]["char"],dims[s]))
    # T4: 2-dim irreps with det = 1
    T4=[]
    for s in range(nc):
        if dims[s]!=2: continue
        ch=irr[s]["char"]
        lam2=[(ch[i]*ch[i]-ch[G.class_sq[i]]).scal(Fraction(1,2)) for i in range(nc)]
        if not all(l.reduced()==[Fraction(1)] for l in lam2): continue
        Sch=[(ch[i]*ch[i]+ch[G.class_sq[i]]).scal(Fraction(1,2)) for i in range(nc)]
        cons=[]
        for s2 in range(nc):
            mm=ip_exact(G,Sch,irr[s2]["char"])
            assert mm.denominator==1
            if mm: cons.append({"sigma":s2,"dim":dims[s2],"d":dist[s2],"mult":int(mm)})
        assert sum(c["dim"]*c["mult"] for c in cons)==3,(name,s,cons)
        taus.append(("S(rho%d)"%s,Sch,3))
        T4.append({"rho_sigma":s,"rho_d":dist[s],"constituents":cons})
    mu={}
    for lab,ch,dt in taus:
        rows=[]
        for m_ in range(2,MMAX+1):
            cE=chiE(m_)
            va=ip_exact(G,cE,ch,conj=False)      # convention A: no conjugate (= spec's invariants)
            vb=ip_exact(G,cE,ch,conj=True)       # convention B
            assert va.denominator==1 and vb.denominator==1
            # third route: branching formula  (m-1)*mult(tau* in V_m) + (m+1)*mult(tau* in V_{m-2})
            f1=ip_exact(G,chi[m_],ch,conj=False); f2=ip_exact(G,chi[m_-2],ch,conj=False)
            vc=(m_-1)*f1+(m_+1)*f2
            rows.append({"m":m_,"convA":int(va),"convB":int(vb),"branch_formula":int(vc),
                         "agree":(va==vb==vc)})
        mu[lab]={"dim_tau":dt,"rows":rows}
    def firstq(lab):
        for r in mu[lab]["rows"]:
            if r["convA"]!=0: return r["m"]
        return None
    def leastk(ch):
        for k in range(AMAX+1):
            if ip_exact(G,chi[k],ch)!=0: return k
        return None
    for rec in T4:
        lab="S(rho%d)"%rec["rho_sigma"]; qq=firstq(lab)
        rec["q"]=qq; rec["q_squared"]=None if qq is None else qq*qq
        chS=dict((l,c) for l,c,_ in taus)[lab]
        kk=leastk(chS); rec["T7_least_k"]=kk
        rec["T7_k_k_plus_2"]=None if kk is None else kk*(kk+2)
    qt=firstq("sigma0"); kt=leastk(irr[0]["char"])
    T5={"tau":"trivial","q":qt,"q_squared":None if qt is None else qt*qt,
        "T7_least_k":kt,"T7_k_k_plus_2":None if kt is None else kt*(kt+2)}
    T6=[{"sigma":s,"dim":dims[s],"d":dist[s],"least_sym_level":T1[s]["least_a"],
         "e":firstq("sigma%d"%s)} for s in range(nc)]
    chartab_exact=[{"sigma":s,"dim":dims[s],"d":dist[s],
        "values_as_Z_zeta_M_coeffs":[[str(x) for x in irr[s]["char"][i].reduced()] for i in range(nc)],
        "values_numeric":[complex_val(irr[s]["char"][i],M) for i in range(nc)]} for s in range(nc)]
    res={"name":name,"character_table_exact":chartab_exact,"zeta_M_note":
         "row i = coefficients c_0..c_{deg Phi_M -1} of the value in the basis 1,z,...,z^(deg-1), z = exp(2*pi*i/M)","order_expected":expected,"order_built":G.N,"order_ok":G.order_ok,
         "exponent":M,"minus_I":G.has_minus_I,"n_classes":nc,"dims":dims,
         "sum_dim_sq":sumsq,"sum_dim_sq_ok":sumsq==G.N,
         "class_sizes":G.csize,"class_elem_orders":G.corder,
         "exact_row_orthogonality_ok":orth_ok,"exact_col_orthogonality_ok":col_ok,
         "A0_levels_needed_for_all_irreps":G.A0,
         "isotypic_vs_innerproduct_branching_agree":iso_agree,
         "branching_dimension_check":dimchk,
         "adjacency":Adj,"adjacency_symmetric":sym,"adjacency_row_weight_ok":rowok,
         "distance_vector":dist,"diameter":diam,
         "T1":T1,"T2":{"minus_I":G.has_minus_I,"n_violations":len(viol),
                       "parity_holds":not viol,"violations":viol[:12]},
         "T3_mu":mu,"T4":T4,"T5":T5,"T6":T6,
         "branching":Mult,"elapsed":round(time.time()-t0,2)}
    if verbose:
        print("%-6s N=%-3d nc=%-2d dims=%-30s diam=%d A0=%d  %.1fs"%(
            name,G.N,nc,str(dims),diam,G.A0,res["elapsed"]))
    return G,res,irr,chi,Mult,dist,mu,taus

if __name__=="__main__":
    G,res,*_=analyse("2T",build_groups()[-3][1],24)
    print(json.dumps({k:res[k] for k in ("dims","distance_vector","diameter","adjacency",
        "exact_row_orthogonality_ok","exact_col_orthogonality_ok","sum_dim_sq_ok",
        "isotypic_vs_innerproduct_branching_agree","branching_dimension_check",
        "adjacency_symmetric","adjacency_row_weight_ok","T4","T5")},indent=1))
    print("T6:",res["T6"])
    print("T1 least_a:",[t["least_a"] for t in res["T1"]])
