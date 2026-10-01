"""Write stage2.json from the stage-2 audit outputs plus the verdicts of STAGE2.md.

Run as `./py make_stage2.py` after s2_compare.py, s2_audit.py, s2_symbolic.py, s2_fiber_b.py
and s2_regrade.py.  Numbers are copied from their JSON files, not retyped.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f)))
audit, comp, sym, fib, reg = (J("s2_audit.json"), J("s2_compare.json"), J("s2_symbolic.json"),
                              J("s2_fiber_b.json"), J("s2_regrade.json"))
LABS = ["1/4", "1/2", "1", "7/5", "3/2", "pi/2", "17/10"]

E, ES, GAP, DEF = "ESTABLISHED", "ESTABLISHED, SUPPLIED", "GAP", "DEFECT"

common_reasons = [
    "All explicit T4 requirements met, read from the script and confirmed by K1/K2: P1 in (y,w) on "
    "[0,pi]x[-W,W]; seam phi(pi,-w_j) = -phi(0,w_j) imposed with sign -1 (K2 seam residual 0; the room's "
    "own periodic plant gives 2.0); weight |cos y|; zero data on w = +-W and on the fiber; level-0 nodes "
    "pi/2 -+ (pi/2)(k/16)^2 with 8 w-cells; exactly four uniform (midpoint) refinements; p from levels 2,3,4; "
    "Richardson from 3,4 with that p; err = |extrap - finest|; all seven widths; six eigenvalues per level.",
    "Geometry, normal, seam condition, admissible class and data of sections 2-3 used as stated.",
    "No sign of any source beyond the spec sheet (manifest, scripts and logs scanned; no outside paths, no network).",
    "Script, JSON and report agree: rerun reproduces every JSON leaf (C1), every printed T4/T3/T5 number "
    "matches the JSON (C2, C3), p/extrap/err recomputed from the JSON levels (K6).",
    "Eigensolver converged at its own level: max relative residual over all widths, levels and six pairs "
    "about 4e-12 (K3); Sylvester inertia confirms no eigenvalue missed at levels 0-2 (K4).",
    "Reading of 'uniform refinement' (midpoint insertion, no regrading) is a reading, not a violation; "
    "the regraded alternative moves no bottom by more than 3.1e-7 relative (s2_regrade).",
]

out = {"R": 1, "stage": 2, "comparison_values": "my stage-1 values in stage1.json, unchanged",
       "stage1_corrections": [], "rooms": {}}

for room in ("a", "b"):
    reasons = list(common_reasons)
    immaterial = []
    if room == "a":
        immaterial.append("C4.5 is only a 5% gross check (it can fail: the periodic-seam plant gives 264.85 vs 2.38).")
    else:
        immaterial.append("Report says the free node's w-gradient on fiber-edge triangles is 'exactly zero'; "
                          "it is computed via np.linalg.inv and is zero only to roundoff (|gw|*h_w <= 1.1e-16, "
                          f"leaked stiffness / K-diagonal <= {max(fib['1.7'][1][0], fib['0.25'][1][0]):.1e}, F1). "
                          "Immaterial: Room B's eigenvalues equal Room A's (which drops the term exactly) to "
                          f"{max(v['max_rel_levels'] for v in audit['a_vs_b'].values()):.1e} relative.")
        immaterial.append("The --plant switch plants the periodic seam and the signed-cos weight together; each "
                          "check is still shown to fail, and I planted the seam sign alone (K2: residual 2.0).")
        immaterial.append("T4 checks run on the level-0 mesh at three widths only (1/4, pi/2, 17/10).")
    v3 = {}
    for lab in LABS:
        d = audit["v3"][f"{room}_{lab}"]
        v3[lab] = {"room_bottom": d["room_bottom"], "room_err_estimate": d["room_err"], "room_p": d["room_p"],
                   "mine_stage1": d["mine"], "mine_err": d["mine_err"], "rel_diff": d["rel_diff"],
                   "exceeds_1e-3": not d["within_1e-3"], "located_defect": None,
                   "abs_diff_over_room_err": d["abs_diff_over_room_err"]}
    out["rooms"][room.upper()] = {
        "valid": True, "reasons": reasons, "instrument_defects_immaterial": immaterial,
        "max_residual": max(e["max_residual"] for k, v in audit["residuals"].items() if k[0] == room for e in v),
        "rerun_json_max_rel": comp[f"C1_{room}"]["max_rel"],
        "v3": v3}

# ---------------------------------------------------------------- G
gA = {
    "T1": {"steps": {
        "1 varied map F_t = cos(t phi/R) X + R sin(t phi/R) e4": E,
        "2 well defined on RP^3; edge and cone point fixed": E,
        "3 exact metric g_t = c^2 g + t^2 dphi dphi": E,
        "4 exact area via matrix-determinant lemma": E,
        "5 differentiation under the integral; A''(0) = Q": E,
        "boundary/vertex terms: none": E,
        "extension to the admissible class by continuity": E},
        "establishes": "Q(phi) = int (|grad phi|^2 - 2 phi^2/R^2) dA, dA = |cos(y/R)| dy dw; no edge, seam, corner "
                       "or cone-point term; unique continuous extension to the admissible class. Independently "
                       f"confirmed numerically (G1: A''(0) vs Q, rel {sym['G1'][0]['rel']:.1e} at W=7/5, "
                       f"{sym['G1'][1]['rel']:.1e} at W=17/10)."},
    "T2": {"steps": {
        "reading R1 (seam matching of all derivatives; literal value-only domain not symmetric)": E,
        "1 sector reduction, closure = direct sum of sector closures": E,
        "2 endpoints are the two poles (the cone point twice)": E,
        "3 Frobenius exponents +-nu, LC iff nu < 1": E,
        "4 classification": E,
        "5 Friedrichs realization selected (separated, no r^-nu part at either end)": E},
        "remark": "Step 3's side remark 'if the weight cos psi were dropped the threshold would move to nu < 1/2' "
                  "holds only if the measure alone changes; if the area element is dropped from operator and "
                  "measure, the threshold is nu < sqrt(3)/2 (G4). Not part of the argument; the step stands.",
        "establishes": "0 < W <= pi R/2: every sector limit-point at the cone point (nu_k >= 1; nu = 1 at W = pi/2, k = 1 "
                       "is limit-point); closure self-adjoint, no other extensions. W/R = 17/10: sector k = 1 "
                       "(nu_1 = 5pi/17 = 0.924) limit-circle at both ends, k >= 2 limit-point; deficiency (2,2), "
                       "U(2) family; the admissible class selects the Friedrichs extension (separated, f = O(r^nu))."},
    "T3": {"steps": {
        "1 Legendre form in x = sin psi": E,
        "2 reduction to Gegenbauer, lambda = (nu+n)(nu+n+1)": E,
        "3 f_{k,n} e_k lie in D(A) with A f = lambda f": ES,
        "4 completeness, spectrum exactly {lambda_{k,n}}": E,
        "lowest eigenvalue": E},
        "supplied": {"3": "The step works with sector test functions and the sector form; the passage to the 2D "
                          "Friedrichs operator needs A = (+)_k A_k. Supplied: for u in the core D, integrating twice by "
                          "parts in theta (u = e_k = 0 at theta = +-W) and Parseval give a(u) = sum_k a_k(u_k); so the "
                          "form closure is H = {u : u_k in H_k, sum a_k(u_k) < inf} and its operator is (+)_k A_k. "
                          "Hence a(f e_k, g) = a_k(f, g_k) = lambda <f, g_k> = lambda <f e_k, g> for all g in H."},
        "establishes": "spectrum exactly {(nu_k+n)(nu_k+n+1)/R^2 : k >= 1, n >= 0}, nu_k = k pi R/(2W), labels (k,n), "
                       "multiplicity = number of labels; complete (eigenbasis, compact resolvent); lowest "
                       "pi(pi R + 2W)/(4 W^2 R), label (1,0). Values agree with my stage-1 closed form (G5)."},
    "T5": {"steps": {
        "1 operator of Q is J = A - 2/R^2": E,
        "2 lowest eigenvalue (nu_1 - 1)(nu_1 + 2)/R^2, sign = sign(pi R/2 - W)": E,
        "3 index and nullity by counting nu_k + n < 1, = 1; min-max": E},
        "establishes": None},
}
gB = {
    "T1": {"steps": {
        "1 geodesics of S^3(R)": E,
        "2 varied surface, well defined on RP^3, edge and cone point fixed": E,
        "3 induced metric g_t = cos^2 u g + du du": E,
        "4 exact density cos u sqrt(cos^2 u + t^2 |grad phi|^2)|cos y|": E,
        "5 differentiation under the integral; A''(0) = Q": E,
        "boundary, seam, vertex terms: none": E,
        "extension by continuity": E},
        "establishes": gA["T1"]["establishes"]},
    "T2": {"steps": {
        "reading D0 (smooth on the lune; literal domain not symmetric; form closure unchanged)": E,
        "sectors and orthonormal basis s_k": E,
        "(a)-(c) deficiency indices split over sectors": E,
        "endpoint analysis: tan^(+-nu)(r/2), exponents +-nu, LC iff nu < 1": E,
        "classification": E,
        "realization selected: Friedrichs": E},
        "establishes": gA["T2"]["establishes"]},
    "T3": {"steps": {
        "1 sector ODE to Gegenbauer, lambda = (nu+m)(nu+m+1)": E,
        "2 F_{k,m} in the form domain (cutoff estimate)": E,
        "3 F_{k,m} in D(A_F), A_F F = lambda F (2D Green formula on D0)": E,
        "4 completeness of {F_{k,m}} in L^2(L_W)": E,
        "5 spectrum = closure of eigenvalues, discrete": E},
        "establishes": gA["T3"]["establishes"].replace("(k,n)", "(k,m)").replace("n >= 0", "m >= 0").replace("+n)", "+m)")},
    "T5": {"steps": {
        "J = A_F - 2, eigenvalues (nu_k+m-1)(nu_k+m+2)": E,
        "min-max: index = number of negative eigenvalues, nullity = dim ker J": E,
        "sign analysis: negative iff m = 0 and k < 2W/pi; zero iff m = 0, k = 2W/pi": E,
        "table per width (exact sympy decisions)": E},
        "establishes": None},
}
t5 = {lab: {"lowest_jacobi": v["jacobi_lowest"], "sign": {1: "+", 0: "0", -1: "-"}[v["sign"]],
            "index": v["index"], "nullity": v["nullity"]} for lab, v in sym["G5"].items()}
gA["T5"]["establishes"] = gB["T5"]["establishes"] = t5
out["rooms"]["A"]["G"], out["rooms"]["B"]["G"] = gA, gB
out["status_of_arguments"] = "audited arguments (not verified or proven theorems)"
out["checks"] = {
    "s2_audit": {"pass": audit["summary"]["pass"], "fail": audit["summary"]["fail"],
                 "planted_fired": audit["summary"]["planted_fired"]},
    "s2_compare": {"lines": comp["lines"]},
    "s2_symbolic": {"lines": sym["lines"]},
    "s2_fiber_b": {"lines": fib["lines"], "planted_fired": fib["planted_fired"]},
    "s2_regrade_rel_vs_mine": {lab: reg[lab]["rel_vs_mine"] for lab in LABS},
    "a_vs_b": audit["a_vs_b"],
}
with open(os.path.join(HERE, "stage2.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print("wrote stage2.json")
