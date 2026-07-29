import json
MINE=json.load(open("_sweep_raw.json")); A=json.load(open("../solverA/m8_1_1_coexact.json"))
AG=A["groups"]; CMP=json.load(open("_compare2.json")); PAT=json.load(open("_patterns.json"))
V=[]
def add(claim,theirs,mine,verdict,note=""):
    V.append({"claim":claim,"agent_A":theirs,"audit":mine,"verdict":verdict,"note":note})

add("group orders: C_n=n (n=1..10), BD_n=4n (n=2..6), 2T=24, 2O=48, 2I=120",
    {n:AG[n]["order_verified"] for n in AG}, {n:MINE[n]["order_built"] for n in AG},
    "CONFIRMED","recomputed with exact GF(p) closure, zero tolerance; all 18 identical")
add("2I closure from the spec generators reaches 120 (no fallback list needed) and equals the explicit 120-quaternion list",
    "120, fallback_used=False, closure_equals_explicit_spec_list=True",
    "120 over GF(p); independent float rerun also 120 and set-identical to the explicit list",
    "CONFIRMED")
add("group exponents", {n:AG[n]["exponent"] for n in AG}, {n:MINE[n]["exponent"] for n in AG},"CONFIRMED")
add("-I in Gamma", {n:AG[n]["minus_I_in_group"] for n in AG}, {n:MINE[n]["minus_I"] for n in AG},"CONFIRMED")
add("irreducible dimensions (sum of squares = |Gamma|)",
    {n:[c["dim"] for c in AG[n]["irreducible_characters"]] for n in AG},
    {n:MINE[n]["dims"] for n in AG},"CONFIRMED",
    "my character table is built by a different algorithm and verified by EXACT row and column orthogonality")
add("adjacency matrix A (McKay graph), symmetric, sum_sigma' A[s][s'] dim s' = 2 dim s",
    "A as reported per group","identical up to the sigma relabelling found by search; "
    "A_symmetric and row-weight hold for all 34 groups","CONFIRMED",
    "agent A and I order complex-conjugate irrep pairs differently; a consistent relabelling reproduces A exactly for all 18 groups")
add("distance vector d(sigma) and graph diameter",
    {n:(AG[n]["distance_vector"],AG[n]["graph_diameter"]) for n in AG},
    {n:(MINE[n]["distance_vector"],MINE[n]["diameter"]) for n in AG},"CONFIRMED","matched under the same relabelling")
add("T1: least a with <chi_sigma, chi_a> != 0 equals d(sigma) for every sigma of every group",
    "holds for all 18 groups","holds for all 34 groups / 344 irreducibles, a<=32","CONFIRMED",
    "widened search found no counterexample")
add("T2: with -I in Gamma every occurrence of sigma in V_a has a = d(sigma) mod 2; without -I it fails",
    "holds for all 8 groups with -I; violations reported for the 10 without",
    "holds for all 21 groups with -I; fails for all 13 without, a<=32","CONFIRMED")
add("T3: the three routes to mu_tau(m) agree",
    "char_tau == char_tau_conj == svd_rank in every row",
    "convA == convB == branching formula in all 11165 rows over 34 groups, m<=30; "
    "and exact GF(p) projector rank agrees in 14 spot checks (dims 30..210), projectors exactly idempotent",
    "CONFIRMED","see the mutation notes: convA==convB is forced by chi_(E_m) being real, so that half of the agreement is not evidence")
add("T4: q and q^2 for every 2-dim determinant-one irreducible",
    "q^2 = 4 for 12 of the 13 such reps; 2I has one with q = 6, q^2 = 36",
    "q^2 = 4 for 40 of 41 such reps over 34 groups; the single exception is the same 2I rep, q = 6, q^2 = 36",
    "CONFIRMED","distinct q^2 values found: [4, 36]")
add("T5: q = 2, q^2 = 4 for the trivial twist, every Gamma",
    "q=2, q^2=4 for all 18","q=2, q^2=4 for all 34","CONFIRMED")
add("T6: e(sigma) table",{n:[t["e"] for t in AG[n]["T6"]] for n in AG},
    {n:[t["e"] for t in MINE[n]["T6"]] for n in AG},"CONFIRMED",
    "e(tau) = min{m>=2 : tau* in V_m or tau* in V_(m-2)} verified with zero violations over all 34 groups")
add("T7: least k and k(k+2)","T4 twists: k=2 -> 8, and 2I's exception k=6 -> 48; trivial: k=0 -> 0",
    "identical; over 34 groups k(k+2) = q^2 + 2q for every T4 twist (0 violations) but NOT for the trivial twist",
    "CONFIRMED","q^2 is never equal to k(k+2) in any of the 41+34 cases")
add("T8: branching of V_a|2I into irreducibles for a = 0..8",
    "as reported","identical (dim, d, multiplicity) for every a","CONFIRMED")
add("T9: rank of the full unreduced projector equals the character sum",
    "SVD rank == character sum in all 7 T9 cases and all 198 T9b cases, tol 1e-8",
    "exact GF(p) rank == character sum in 14 independent cases; every projector exactly idempotent over GF(p)",
    "CONFIRMED","tolerance is not load-bearing: smallest kept singular value ~1.0, largest dropped ~2.9e-15")
add("audit metric: min pairwise element distance over all groups = 0.43701602444882115",
    0.43701602444882115,0.43701602444882115,"CONFIRMED","independent float recomputation")
json.dump(V,open("_verdicts.json","w"),indent=1)
for v in V: print("%-9s %s"%(v["verdict"],v["claim"][:95]))
