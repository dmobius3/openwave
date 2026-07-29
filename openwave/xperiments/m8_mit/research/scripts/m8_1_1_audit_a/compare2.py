"""Rigorous adjacency comparison: search for an irrep relabelling consistent with every invariant."""
import json, itertools
MINE=json.load(open("_sweep_raw.json")); AG=json.load(open("../solverA/m8_1_1_coexact.json"))["groups"]
out={}
for name in AG:
    m=MINE[name]; a=AG[name]; n=len(m["dims"])
    Amult=a["branching_multiplicities_V_a"]
    sig=lambda dims,dist,br,s: (dims[s],dist[s],tuple(br[aa][s] for aa in range(15)))
    sM=[sig(m["dims"],m["distance_vector"],m["branching"],s) for s in range(n)]
    sA=[sig([c["dim"] for c in a["irreducible_characters"]],a["distance_vector"],Amult,s) for s in range(n)]
    if sorted(sM)!=sorted(sA):
        out[name]={"status":"INVARIANT MISMATCH"}; print(name,"INVARIANT MISMATCH"); continue
    cand=[[t for t in range(n) if sA[t]==sM[s]] for s in range(n)]
    MA=m["adjacency"]; AA=a["adjacency_A"]
    assign=[None]*n; used=[False]*n; found=[None]
    def bt(s):
        if found[0]: return
        if s==n: found[0]=assign[:]; return
        for t in cand[s]:
            if used[t]: continue
            ok=True
            for s2 in range(s):
                if MA[s][s2]!=AA[t][assign[s2]] or MA[s2][s]!=AA[assign[s2]][t]: ok=False; break
            if ok and MA[s][s]!=AA[t][t]: ok=False
            if ok:
                assign[s]=t; used[t]=True; bt(s+1); used[t]=False; assign[s]=None
    bt(0)
    if found[0] is None:
        out[name]={"status":"NO CONSISTENT RELABELLING -> adjacency graphs differ"}
        print(name,"ADJACENCY GRAPHS DIFFER")
    else:
        perm=found[0]
        chk=all(MA[i][j]==AA[perm[i]][perm[j]] for i in range(n) for j in range(n))
        # verify every other quantity under the SAME permutation
        errs=[]
        for s in range(n):
            t=perm[s]
            if m["distance_vector"][s]!=a["distance_vector"][t]: errs.append(("dist",s))
            if m["T1"][s]["least_a"]!=a["T1"][t]["least_a"]: errs.append(("least_a",s))
            if m["T6"][s]["e"]!=a["T6"][t]["e"]: errs.append(("e",s))
            mr={r["m"]:r["convA"] for r in m["T3_mu"]["sigma%d"%s]["rows"] if r["m"]<=12}
            ar={r["m"]:r["char_tau"] for r in a["T3_mu_tables"]["sigma%d"%t]["rows"]}
            sr={r["m"]:r["svd_rank"] for r in a["T3_mu_tables"]["sigma%d"%t]["rows"]}
            if mr!=ar: errs.append(("mu_char",s))
            if mr!=sr: errs.append(("mu_svd",s))
            for aa in range(15):
                if m["branching"][aa][s]!=Amult[aa][t]: errs.append(("branch",aa,s))
        out[name]={"status":"MATCH" if (chk and not errs) else "MISMATCH",
                   "adjacency_isomorphic":chk,"other_errors":errs[:10]}
        print("%-6s adjacency relabelling found, A identical=%s, other errors=%d"%(name,chk,len(errs)))
json.dump(out,open("_compare2.json","w"))
print("\nALL MATCH:", all(v.get("status")=="MATCH" for v in out.values()))
