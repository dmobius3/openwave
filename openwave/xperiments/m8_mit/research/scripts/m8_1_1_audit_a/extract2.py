import json
J=json.load(open("../solverA/m8_1_1_coexact.json"))
for n,b in J["groups"].items():
    dims=[c["dim"] for c in b["irreducible_characters"]]
    print("%-5s N=%-3d nc=%-2d dims=%-24s dist=%-26s diam=%d -I=%s"%(
      n,b["order_verified"],len(dims),dims,b["distance_vector"],b["graph_diameter"],b["minus_I_in_group"]))
    print("     T1 least_a =", [t["least_a"] for t in b["T1"]])
    print("     T6 e(sigma)=", [t["e"] for t in b["T6"]])
    print("     T2 holds=%s nviol=%d applicable=%s"%(b["T2"]["parity_rule_holds"],b["T2"]["n_violations"],b["T2"]["rule_applicable"]))
    for r in b["T4"]:
        print("     T4 rho=sigma%d cons=%s q=%s q2=%s leastk=%s k(k+2)=%s"%(
          r["rho_sigma"],[(c["sigma"],c["dim"],c["d"],c["mult"]) for c in r["constituents"]],
          r["q"],r["q_squared"],r["T7_least_k"],r["T7_k_k_plus_2"]))
    t5=b["T5"]; print("     T5 trivial q=%s q2=%s leastk=%s k(k+2)=%s"%(t5["q"],t5["q_squared"],t5["T7_least_k"],t5["T7_k_k_plus_2"]))
    # check all mu agree
    bad=[]
    for lab,tab in b["T3_mu_tables"].items():
        for row in tab["rows"]:
            if not row["agree"]: bad.append((lab,row["m"],row["char_tau"],row["char_tau_conj"],row["svd_rank"]))
    print("     mu disagreements:", bad if bad else "none")
