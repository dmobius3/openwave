import json, time
from fractions import Fraction
from driver import analyse, AMAX, MMAX
from part1 import build_groups, p
from proj import explicit_irrep, projector_rank_exact
from part2 import Cyc

GL = build_groups(nc_max=20, nbd_max=12)
results={}; keep={}
t0=time.time()
for name,gens,order in GL:
    G,res,irr,chi,Mult,dist,mu,taus = analyse(name,gens,order)
    results[name]=res; keep[name]=(G,irr,Mult,dist)
print("total %.1fs"%(time.time()-t0))
json.dump(results,open("_sweep_raw.json","w"))

# ---------- exact projector-rank checks (tolerance-free replacement for SVD) ----------
checks=[]
for gname,m_,dimt in [("BD_3",4,1),("BD_3",5,2),("2T",4,1),("2T",5,3),("2T",6,3),
                      ("2I",6,1),("2I",4,3),("2O",4,2),("2O",5,3),("BD_6",4,2),
                      ("C_5",4,1),("C_7",6,1),("BD_12",4,2),("2I",5,4)]:
    G,irr,Mult,dist=keep[gname]
    for s,x in enumerate(irr):
        if x["dim"]!=dimt: continue
        a=next((aa for aa in range(AMAX+1) if Mult[aa][s]==1),None)
        if a is None: continue
        ct=[None]*G.nc
        # GF(p) character values of this irrep, recovered from the exact Cyc values
        zM=G.zM
        for ci in range(G.nc):
            ct[ci]=sum(c*pow(zM,l,p) for l,c in enumerate(x["char"].__getitem__(ci).c))%p
        T=explicit_irrep(G,a,dimt,ct)
        if T is None: continue
        r=projector_rank_exact(G,m_,T,dimt)
        r.update({"group":gname,"m":m_,"sigma":s,"dim_tau":dimt,
                  "character_sum_convA":[row["convA"] for row in results[gname]["T3_mu"]["sigma%d"%s]["rows"] if row["m"]==m_][0]})
        r["agree"]= (r["rank"]==r["character_sum_convA"]) if r["rank"] is not None else None
        checks.append(r); break
for c in checks:
    print("EXACTRANK %-5s m=%-2d sigma%-2d dimtau=%d  dim=%-4s rank=%-5s charsum=%-3s agree=%s idem=%s"%(
        c["group"],c["m"],c["sigma"],c["dim_tau"],c["dim"],c["rank"],c["character_sum_convA"],
        c["agree"],c.get("idempotent_exactly")))
json.dump(checks,open("_projchecks.json","w"))
