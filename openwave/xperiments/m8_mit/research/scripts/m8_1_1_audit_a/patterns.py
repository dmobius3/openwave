import json, collections
M=json.load(open("_sweep_raw.json"))
AMAX=32; MMAX=30
fail=collections.defaultdict(list); cov=collections.Counter()
q2vals=collections.Counter(); q2detail=[]
maxe=0; maxleast=0; maxq=0
for name,g in M.items():
    n=len(g["dims"]); dist=g["distance_vector"]
    cov["groups"]+=1; cov["irreps"]+=n
    if not g["order_ok"]: fail["order"].append(name)
    if not g["sum_dim_sq_ok"]: fail["sum_dim_sq"].append(name)
    if not g["exact_row_orthogonality_ok"]: fail["row_orth"].append(name)
    if not g["exact_col_orthogonality_ok"]: fail["col_orth"].append(name)
    if not g["isotypic_vs_innerproduct_branching_agree"]: fail["iso_vs_ip"].append(name)
    if not g["branching_dimension_check"]: fail["dimcheck"].append(name)
    if not g["adjacency_symmetric"]: fail["A_sym"].append(name)
    if not g["adjacency_row_weight_ok"]: fail["A_rowweight"].append(name)
    if g["diameter"]!=max(dist): fail["diam_eq_max_dist_from_trivial"].append(name)
    for s in range(n):
        la=g["T1"][s]["least_a"]; e=g["T6"][s]["e"]; d=dist[s]
        maxleast=max(maxleast,la if la is not None else -1); maxe=max(maxe,e if e is not None else -1)
        if la!=d: fail["P1_least_a_eq_d"].append((name,s,la,d))
        exp_e = d if d>=2 else d+2
        if e!=exp_e: fail["P2_e_formula"].append((name,s,d,e,exp_e))
        cov["mu_rows"]+=len(g["T3_mu"]["sigma%d"%s]["rows"])
        for r in g["T3_mu"]["sigma%d"%s]["rows"]:
            if not r["agree"]: fail["P5_mu_three_methods"].append((name,s,r))
    par_ok = (g["T2"]["n_violations"]==0)
    if par_ok != g["minus_I"]: fail["P4_parity_iff_minusI"].append((name,g["minus_I"],par_ok))
    if g["T5"]["q"]!=2 or g["T5"]["q_squared"]!=4: fail["P6_trivial_q_is_2"].append((name,g["T5"]))
    for r in g["T4"]:
        cov["two_dim_det1_reps"]+=1
        q2vals[r["q_squared"]]+=1; maxq=max(maxq,r["q"])
        q2detail.append({"group":name,"rho_d":r["rho_d"],
                         "constituent_dims_d":[(c["dim"],c["d"],c["mult"]) for c in r["constituents"]],
                         "q":r["q"],"q2":r["q_squared"],
                         "T7_least_k":r["T7_least_k"],"T7_k_k_plus_2":r["T7_k_k_plus_2"]})
        if r["q"]!=r["T7_least_k"]: fail["P9_q_eq_leastk_for_T4"].append((name,r["q"],r["T7_least_k"]))
        # q should equal min over constituents of d, when that min >= 2, else +2
        md=min(c["d"] for c in r["constituents"])
        exp=md if md>=2 else md+2
        if r["q"]!=exp: fail["P8_q_from_constituent_distance"].append((name,r["q"],exp,md))
        for lab in ["S(rho%d)"%r["rho_sigma"]]:
            for row in g["T3_mu"][lab]["rows"]:
                if not row["agree"]: fail["P5_mu_three_methods"].append((name,lab,row))
                cov["mu_rows"]+=1
print("COVERAGE:",dict(cov))
print("max least_a over all groups:",maxleast," max e:",maxe," max q:",maxq,
      " (agent A ranges: a<=14, m<=12)")
print()
print("q^2 VALUE HISTOGRAM over every 2-dim det-1 irrep of every group:",dict(q2vals))
print("distinct q^2 values:",sorted(q2vals))
print()
for k,v in fail.items():
    print("PATTERN BREAK %-34s count=%d  examples=%s"%(k,len(v),v[:4]))
if not fail: print("NO PATTERN BROKE anywhere in the widened search space.")
print()
print("2-dim det-1 reps with q^2 != 4:")
for r in q2detail:
    if r["q2"]!=4: print("   ",r)
print()
print("count of 2-dim det-1 reps per group:")
print({n:len(M[n]["T4"]) for n in M if M[n]["T4"]})
json.dump({"coverage":dict(cov),"q2_hist":{str(k):v for k,v in q2vals.items()},
           "q2_detail":q2detail,"pattern_breaks":{k:v[:20] for k,v in fail.items()},
           "max_least_a":maxleast,"max_e":maxe,"max_q":maxq},open("_patterns.json","w"))
