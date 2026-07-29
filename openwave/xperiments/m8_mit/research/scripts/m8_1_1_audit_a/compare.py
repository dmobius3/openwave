import json
MINE=json.load(open("_sweep_raw.json"))
A=json.load(open("../solverA/m8_1_1_coexact.json"))
AG=A["groups"]
diffs=[]; notes=[]
def D(msg): diffs.append(msg); print("DIFF:",msg)

for name in AG:
    m=MINE[name]; a=AG[name]
    if m["order_built"]!=a["order_verified"]: D("%s order %s vs %s"%(name,m["order_built"],a["order_verified"]))
    if m["exponent"]!=a["exponent"]: D("%s exponent %s vs %s"%(name,m["exponent"],a["exponent"]))
    if m["minus_I"]!=a["minus_I_in_group"]: D("%s -I %s vs %s"%(name,m["minus_I"],a["minus_I_in_group"]))
    ad=[c["dim"] for c in a["irreducible_characters"]]
    if sorted(m["dims"])!=sorted(ad): D("%s dim multiset %s vs %s"%(name,sorted(m["dims"]),sorted(ad)))
    if m["diameter"]!=a["graph_diameter"]: D("%s diameter %s vs %s"%(name,m["diameter"],a["graph_diameter"]))
    # match sigma indices by branching signature over a=0..14
    Amult=a["branching_multiplicities_V_a"]      # [a][sigma]
    sigA={s:tuple(Amult[aa][s] for aa in range(15)) for s in range(len(ad))}
    sigM={s:tuple(m["branching"][aa][s] for aa in range(15)) for s in range(len(m["dims"]))}
    if len(set(sigA.values()))!=len(sigA): notes.append("%s: agentA branching signatures not unique"%name)
    inv={v:k for k,v in sigA.items()}
    perm={}
    for s,v in sigM.items():
        if v in inv: perm[s]=inv[v]
        else: D("%s sigma%d branching column %s has no counterpart in agent A"%(name,s,v))
    if len(perm)!=len(sigM): continue
    # distances, T1, T6
    for s,t in perm.items():
        if m["distance_vector"][s]!=a["distance_vector"][t]:
            D("%s sigma%d(mine)~%d(A) distance %s vs %s"%(name,s,t,m["distance_vector"][s],a["distance_vector"][t]))
        if m["T1"][s]["least_a"]!=a["T1"][t]["least_a"]:
            D("%s sigma%d least_a %s vs %s"%(name,s,m["T1"][s]["least_a"],a["T1"][t]["least_a"]))
        if m["T6"][s]["e"]!=a["T6"][t]["e"]:
            D("%s sigma%d e %s vs %s"%(name,s,m["T6"][s]["e"],a["T6"][t]["e"]))
    # adjacency under the permutation
    for s1 in perm:
        for s2 in perm:
            if m["adjacency"][s1][s2]!=a["adjacency_A"][perm[s1]][perm[s2]]:
                D("%s A[%d][%d] %s vs %s"%(name,s1,s2,m["adjacency"][s1][s2],a["adjacency_A"][perm[s1]][perm[s2]]))
    # branching over the full common range a=0..14
    for aa in range(15):
        for s in perm:
            if m["branching"][aa][s]!=Amult[aa][perm[s]]:
                D("%s branching a=%d sigma%d %s vs %s"%(name,aa,s,m["branching"][aa][s],Amult[aa][perm[s]]))
    # T2
    if bool(m["T2"]["parity_holds"])!=bool(a["T2"]["parity_rule_holds"]) and a["T2"]["rule_applicable"]:
        D("%s T2 parity %s vs %s"%(name,m["T2"]["parity_holds"],a["T2"]["parity_rule_holds"]))
    # T4
    mt4=sorted([(m["distance_vector"][r["rho_sigma"]],
                 tuple(sorted((c["dim"],c["d"],c["mult"]) for c in r["constituents"])),
                 r["q"],r["q_squared"],r["T7_least_k"],r["T7_k_k_plus_2"]) for r in m["T4"]])
    at4=sorted([(a["distance_vector"][r["rho_sigma"]],
                 tuple(sorted((c["dim"],c["d"],c["mult"]) for c in r["constituents"])),
                 r["q"],r["q_squared"],r["T7_least_k"],r["T7_k_k_plus_2"]) for r in a["T4"]])
    if mt4!=at4: D("%s T4 %s vs %s"%(name,mt4,at4))
    # T5
    for k in ("q","q_squared","T7_least_k","T7_k_k_plus_2"):
        if m["T5"][k]!=a["T5"][k]: D("%s T5 %s %s vs %s"%(name,k,m["T5"][k],a["T5"][k]))
    # mu tables m=2..12 for every irreducible tau
    for s,t in perm.items():
        mr={r["m"]:r["convA"] for r in m["T3_mu"]["sigma%d"%s]["rows"] if r["m"]<=12}
        ar={r["m"]:r["char_tau"] for r in a["T3_mu_tables"]["sigma%d"%t]["rows"]}
        if mr!=ar: D("%s mu sigma%d mine=%s A=%s"%(name,s,mr,ar))
        arb={r["m"]:r["svd_rank"] for r in a["T3_mu_tables"]["sigma%d"%t]["rows"]}
        if mr!=arb: D("%s mu-vs-svdrank sigma%d %s vs %s"%(name,s,mr,arb))
    # T8 for 2I
    if name=="2I":
        for blk in a["T8_branching_2I"]:
            aa=blk["a"]
            got=sorted((c["dim"],c["d"],c["mult"]) for c in blk["constituents"])
            mine=sorted((m["dims"][s],m["distance_vector"][s],m["branching"][aa][s])
                        for s in range(len(m["dims"])) if m["branching"][aa][s])
            if got!=mine: D("2I T8 a=%d %s vs %s"%(aa,mine,got))
print()
print("TOTAL DIFFS vs agent A over the 18 spec groups:",len(diffs))
print("NOTES:",notes)
json.dump({"diffs":diffs,"notes":notes},open("_compare.json","w"))
