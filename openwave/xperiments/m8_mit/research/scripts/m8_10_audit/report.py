"""Compare the exact route with the explicit-function route, compute item 8, write audit_results.json."""
import json

import sympy as sp

import lib

C = lib.Checker("report")
R = lib.ROOM
alg = json.load(open(R / "algebra_primary.json"))
fp = json.load(open(R / "functional_primary.json"))
fs = json.load(open(R / "functional_secondary.json"))
i01 = json.load(open(R / "item01.json"))
cgd = json.load(open(R / "check_cg_D.json"))
i5 = json.load(open(R / "item5.json"))
i5b = json.load(open(R / "item5b.json"))   # marker: item5b.json

S = lambda s: sp.nsimplify(sp.sympify(s)) if s is not None else None
rays = ["R1", "R2", "R3", "R4", "R5"]


def fl(x):
    return float(sp.N(S(x), 30))


# ---------------- exact vs functional (both identifications)
maxdev = {}
for tag, fj in (("primary", fp), ("secondary", fs)):
    dev = 0.0
    dev_mut = 0.0
    for d in ("3", "4"):
        for r in rays:
            a = alg[d]["rays"][r]
            f = fj[d][r]
            pairs = [(fl(a["item7_lambda4_over_g2"]), f["along"]),
                     (fl(a["item7_lambda4_over_g2"]), f["lam4_identity"]),
                     (fl(a["item6_perp_norm2"]), f["perp_norm2"]),
                     (fl(a["lambda2_over_g"]), f["multiple"] if a["item3_is_multiple"] else fl(a["lambda2_over_g"]))]
            for n, x in a["item4_normPi_xi2_over_g2_by_n"].items():
                pairs.append((fl(x), f["xi_norms_by_n"][str(n)] if str(n) in f["xi_norms_by_n"] else f["xi_norms_by_n"][n]))
            for J, x in a["normPiN2_by_J"].items():
                pairs.append((fl(x), f["normPiN2_by_J"][str(J)]))
            if a["item3_is_multiple"]:
                pairs.append((fl(a["item3_multiple"]), f["multiple"]))
            for x, y in pairs:
                dev = max(dev, abs(x - y) / max(1.0, abs(x)))
            # mutant: compare against the other sector's value of lambda_4
            other = alg["4" if d == "3" else "3"]["rays"][r]
            dev_mut = max(dev_mut, abs(fl(other["item7_lambda4_over_g2"]) - f["along"]))
    maxdev[tag] = dev
    C.check(f"exact route == explicit-function route ({tag} identification), max rel dev {dev:.1e}",
            lambda z: z < 1e-9, dev, dev_mut, "compare lambda_4 against the other sector's exact value")

# ---------------- exact route, secondary identification: all exact scalars identical
als = json.load(open(R / "algebra_secondary.json"))
keys = ["item3_multiple", "item7_lambda4_over_g2", "item6_perp_norm2", "item6_along", "lambda2_over_g"]
ndiff = 0
for d in ("3", "4"):
    for r in rays:
        for k in keys:
            x, y = alg[d]["rays"][r][k], als[d]["rays"][r][k]
            if (x is None) != (y is None) or (x is not None and sp.simplify(S(x) - S(y)) != 0):
                ndiff += 1
        for n, x in alg[d]["rays"][r]["item4_normPi_xi2_over_g2_by_n"].items():
            if sp.simplify(S(x) - S(als[d]["rays"][r]["item4_normPi_xi2_over_g2_by_n"][n])) != 0:
                ndiff += 1
    for K, x in alg[d]["item2_kappa"].items():
        if sp.simplify(S(x) - S(als[d]["item2_kappa"][K])) != 0:
            ndiff += 1
mut = sum(sp.simplify(S(alg["3"]["rays"][r]["item7_lambda4_over_g2"]) -
                      S(als["4"]["rays"][r]["item7_lambda4_over_g2"])) != 0 for r in rays)
C.check("exact values identical under the secondary quaternion identification (items 2,3,4,6,7)",
        lambda z: z == 0, ndiff, mut, "compare sector-3 lambda_4 against sector-4 lambda_4")

# ---------------- item 8
item8 = {}
ratios = {}
mratios = {}
for r in rays:
    l3 = S(alg["3"]["rays"][r]["item7_lambda4_over_g2"])
    l4 = S(alg["4"]["rays"][r]["item7_lambda4_over_g2"])
    rat = sp.nsimplify(sp.radsimp(l3 / l4))
    ratios[r] = rat
    m3 = alg["3"]["rays"][r]["item3_multiple"]
    m4 = alg["4"]["rays"][r]["item3_multiple"]
    mr = None
    if m3 is not None and m4 is not None:
        mr = sp.nsimplify(sp.radsimp((S(m3) - 1) / (S(m4) - 1)))
    mratios[r] = mr
    item8[r] = {"lambda4_ratio_3_over_4": str(rat), "lambda4_ratio_float": float(sp.N(rat)),
                "multiple_minus_1_ratio": str(mr) if mr is not None else None,
                "multiple_minus_1_ratio_float": float(sp.N(mr)) if mr is not None else None}
    print(f"item 8 {r}: lambda4(3)/lambda4(4) = {rat} = {float(sp.N(rat)):.15g} ; "
          f"(mult-1) ratio = {mr}")
const_l4 = len({sp.nsimplify(x) for x in ratios.values()}) == 1
const_m = len({x for x in mratios.values() if x is not None}) == 1
item8["lambda4_single_c"] = const_l4
item8["multiple_minus_1_single_c"] = const_m
print("one c for lambda_4 at every ray:", const_l4, "; one c for (multiple - 1):", const_m)

# ---------------- assemble
out = {"item0": i01["item0"], "item1": i01["item1"], "item5": {"part_a": i5, "part_b": i5b}}
for d in ("3", "4"):
    sec = alg[d]
    out.setdefault("item2", {})[d] = {"kappa_K_exact": sec["item2_kappa"],
                                     "kappa_K_float": sec["item2_kappa_float"],
                                     "K_with_nonzero_right_invariant": sec["B_nonzero_K"]}
    for r in rays:
        a = sec["rays"][r]
        out.setdefault("item3", {}).setdefault(d, {})[r] = {
            "is_multiple": a["item3_is_multiple"], "multiple": a["item3_multiple"],
            "multiple_float": a["item3_multiple_float"], "normalised_v": a["v_normalised"],
            "w": a["w"]}
        out.setdefault("item4", {}).setdefault(d, {})[r] = {
            "normPi_n_xi_sq_over_g2": a["item4_normPi_xi2_over_g2_by_n"],
            "float": a["item4_float"], "normPi_N_sq_by_J": a["normPiN2_by_J"]}
        out.setdefault("item6", {}).setdefault(d, {})[r] = {
            "along_Phi_over_g": a["item6_along"], "along_float": a["item6_along_float"],
            "orthogonal_fibre_vector_over_g": a["item6_perp_fibre"],
            "orthogonal_norm_sq_over_g2": a["item6_perp_norm2"],
            "orthogonal_norm_sq_float": a["item6_perp_norm2_float"]}
        out.setdefault("item7", {}).setdefault(d, {})[r] = {
            "lambda4_over_g2": a["item7_lambda4_over_g2"], "float": a["item7_float"],
            "lambda2_over_g": a["lambda2_over_g"]}
        out.setdefault("item9", {}).setdefault(d, {})[r] = {
            "exact_coefficient_residual_max": a["item9_exact_coefficient_residual_max"],
            "pointwise_residual_casimir_primary": fp[d][r]["item9_residual_casimir_max"],
            "pointwise_residual_fd_primary": fp[d][r]["item9_residual_fd_max"],
            "pointwise_residual_casimir_secondary": fs[d][r]["item9_residual_casimir_max"],
            "pointwise_residual_fd_secondary": fs[d][r]["item9_residual_fd_max"],
            "scale_max_abs_N": fp[d][r]["item9_scale_maxN"]}
out["item8"] = item8
out["crosscheck_max_rel_dev"] = maxdev
allrec = (cgd["records"] + i01["checks"] + alg["checks"] + fp["checks"] + fs["checks"]
          + i5["checks"] + i5b["checks"] + C.records)
out["checks_summary"] = {"n": len(allrec), "n_pass": sum(r["status"] == "PASS" for r in allrec),
                         "failed": [r["label"] for r in allrec if r["status"] != "PASS"]}
out["algebra_flags"] = alg.get("flags", [])
json.dump(out, open(R / "audit_results.json", "w"), indent=1, default=str)
json.dump({"records": C.records}, open(R / "report_checks.json", "w"), indent=1)
print("checks:", out["checks_summary"]["n_pass"], "/", out["checks_summary"]["n"], "PASS;",
      "failed:", out["checks_summary"]["failed"])
