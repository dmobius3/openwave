"""Numeric identification of every engine number at two precisions, cross-checked against the exact route and the
float64 quadrature; writes results.json.
Rule: at 60 digits, best rational with denominator <= 1e20 accepted if |x - p/q| < 1e-50;
      at 110 digits, same bound, accepted if |x - p/q| < 1e-100; both must give the same p/q.
      |x| below the tolerance counts as numerically zero (exact zeros are established separately in exact_parts.py)."""
import json
from fractions import Fraction as Fr
import mpmath
from common import HERE, identify_rational

E = {p: json.loads((HERE / f"out_engine_{p}.json").read_text()) for p in (60, 110)}
X = json.loads((HERE / "out_exact.json").read_text())
Qd = json.loads((HERE / "out_quadcheck.json").read_text())
I01 = json.loads((HERE / "out_item01.json").read_text())
RL = {p: json.loads((HERE / f"out_relations_{p}.json").read_text()) for p in (60, 110)}
mp = mpmath.mp
BOUND = 10 ** 20
fails = []


def ident(key, getter):
    res = []
    for p in (60, 110):
        mp.dps = p
        x = mp.mpf(getter(E[p]))
        tol = mp.mpf(10) ** (-(p - 10))
        if abs(x) < tol:
            res.append(Fr(0))
            continue
        fr, err = identify_rational(x, mp, BOUND, tol)
        res.append(fr)
    if res[0] is None or res[0] != res[1]:
        fails.append((key, res))
        return None
    return res[0]


results = {"item0": {}, "item1": {}, "item2": {}, "item3": {}, "item4": {}, "item5": {}, "item6": {}, "item7": {},
           "item8": {}, "item9": {}, "checks": {}}
results["item0"] = {"norm2_q1": I01["norm2_q1"], "norm2_q2": I01["norm2_q2"], "order": I01["order"],
                    "derived_subgroup_order": I01["derived_order"], "equals_derived_subgroup": I01["perfect"],
                    "cg_33_3m3_60": "1/sqrt(924)", "cg_33_3m3_60_squared": I01["cg_33_3m3_60_squared"],
                    "cg_33_3m3_60_float": I01["cg_33_3m3_60_float"]}
results["item1"] = {"3dim": I01["mult_3dim"], "4dim": I01["mult_4dim"], "trivial_for_reference": I01["mult_trivial"],
                    "classes": I01["classes"], "character_3dim": I01["character_3dim"],
                    "character_4dim": I01["character_4dim"], "declaration": "DERIVED"}
for d in ("3", "4"):
    sx = X[f"sector_{d}"]
    r0n = ident(f"r0sq_{d}", lambda e: mp.mpf(e["item2"][d]["r0"]) ** 2)
    if str(r0n) != sx["r0_squared"]:
        fails.append(("r0sq mismatch", d))
    r6n2 = ident(f"r6sq_{d}", lambda e: mp.mpf(e["item2"][d]["r6"]) ** 2)
    results["item2"][f"{d}dim"] = {"r0": sx["r0_sign"] + "sqrt(" + sx["r0_squared"] + ")",
                                   "r6": sx["r6_sign"] + "sqrt(" + sx["r6_squared"] + ")",
                                   "r_L_for_L=1..5": "0 (T_L = 0 exactly)",
                                   "numeric_r6_squared_identified": str(r6n2),
                                   "T_L_norm2": {str(L): sx[f"T{L}_norm2"] for L in range(7)},
                                   "F_J": {str(J): sx[f"F_{J}"] for J in range(3, 10)}}
    if str(r6n2) != sx["r6_squared"]:
        fails.append(("r6sq mismatch", d))
results["item2"]["M_K_in_B_basis"] = RL[110]["W"]
assert RL[60]["W"] == RL[110]["W"]
results["item2"]["M_K_rank"] = E[110]["item2"]["M_K_rank"]

rays = ["R1", "R2", "R3", "R4", "R5"]
maxq = 0
for d in ("3", "4"):
    for R in rays:
        key = f"{d}:{R}"
        ex = X["assembled"][key]
        kap = ident(key + " kappa", lambda e: e["rays"][key]["kappa_re"])
        par = ident(key + " par", lambda e: mp.mpf(e["rays"][key]["parallel_resid"]) ** 2)
        if str(kap) != ex["kappa"]:
            fails.append((key, "kappa exact vs numeric", kap, ex["kappa"]))
        results["item3"][key] = {"parallel": par == 0, "multiple": ex["kappa"], "multiple_minus_1": ex["kappa_minus_1"],
                                 "parallel_resid_norm_110digits": E[110]["rays"][key]["parallel_resid"][:12],
                                 "exact_parallel_resid2": ex["parallel_resid2"]}
        xin = {}
        for n in ex["xi_norm2_over_g2"]:
            v = ident(key + f" xi{n}", lambda e: e["rays"][key]["xi_level_norm2_over_g2"][n])
            if str(v) != ex["xi_norm2_over_g2"][n]:
                fails.append((key, n, v, ex["xi_norm2_over_g2"][n]))
            xin[n] = str(v)
        results["item4"][key] = {"sector_levels": ex["sector_levels_n"],
                                 "norm2_over_g2_by_level": {n: xin[n] for n in xin if int(n) in ex["sector_levels_n"]},
                                 "level_6": "0 (xi orthogonal to the block by construction)",
                                 "non_sector_levels_n<=18": "0 (no sections there)"}
        results["item5"][key] = [int(n) for n in xin if int(n) in ex["sector_levels_n"] and xin[n] == "0"]
        al = ident(key + " along", lambda e: e["rays"][key]["DN_along_re"])
        al_im = ident(key + " along_im", lambda e: e["rays"][key]["DN_along_im"])
        # orthogonal norm^2: denominators near 1e20 -> use 110 and 160 digits with bound 1e40 (tol 1e-100, 1e-150)
        o_res = []
        for p, tl in ((110, 100), (160, 150)):
            mp.dps = p
            Ep = json.loads((HERE / f"out_engine_{p}.json").read_text())
            x = mp.mpf(Ep["rays"][key]["DN_orth_norm2"])
            if abs(x) < mp.mpf(10) ** (-tl):
                o_res.append(Fr(0))
            else:
                o_res.append(identify_rational(x, mp, 10 ** 40, mp.mpf(10) ** (-tl))[0])
        if o_res[0] is None or o_res[0] != o_res[1]:
            fails.append((key, "orth2", o_res))
        orth2 = o_res[0]
        blk = max(mp.mpf(E[p]["rays"][key]["DN_blockform_resid"]) for p in (110,))
        yp = E[110]["rays"][key]["DN_fibre_yperp"]
        mp.dps = 110
        # direction of the orthogonal fibre vector: compare with sin t v_2 - cos t v_{-3} (only meaningful when nonzero)
        yz = [mp.mpc(a, b) for a, b in yp]
        ynorm = mp.sqrt(sum(abs(z) ** 2 for z in yz))
        direction = "none (zero)"
        if ynorm > mp.mpf(10) ** (-90):
            ref = [mp.mpf(0)] * 7
            ref[5], ref[0] = 2 * mp.sqrt(3) / 5, -mp.sqrt(13) / 5          # sin t at m=2, -cos t at m=-3
            ov = abs(sum(r * z for r, z in zip(ref, yz))) / ynorm
            direction = f"along sin t v_2 - cos t v_-3: |overlap| = {mpmath.nstr(ov, 30)} (1 = parallel)"
        results["item6"][key] = {"along_Phi_coefficient_over_g": str(al), "along_imag": str(al_im),
                                 "orthogonal_norm2_over_g2": str(orth2),
                                 "orthogonal_direction": direction,
                                 "blockform_resid_110": mpmath.nstr(blk, 3)}
        results["item7"][key] = {"lambda4_over_g2": ex["lambda4_over_g2"], "sign": "negative",
                                 "numeric_identified": str(al),
                                 "a5_block_equation_consistent": orth2 == 0}
        if str(al) != ex["lambda4_over_g2"]:
            fails.append((key, "lambda4", al, ex["lambda4_over_g2"]))
        # quadrature agreement (float64)
        q = Qd[key]
        e110 = E[110]["rays"][key]
        dq = [abs(q["kappa"] - float(e110["kappa_re"])), abs(q["DN_along"] - float(e110["DN_along_re"])),
              abs(q["DN_orth"] - float(e110["DN_orth_norm"]))]
        dq += [abs(q["xi_level_norm2"][n] - float(e110["xi_level_norm2_over_g2"][n])) for n in q["xi_level_norm2"]]
        maxq = max(maxq, max(dq))
        E160 = json.loads((HERE / "out_engine_160.json").read_text())
        fmt = lambda s: mpmath.nstr(mpmath.mpf(s), 3)
        results["item9"][key] = {"order3_residual_60": fmt(E[60]["rays"][key]["order3_casimir_resid"]),
                                 "order3_residual_110": fmt(e110["order3_casimir_resid"]),
                                 "order3_residual_160": fmt(E160["rays"][key]["order3_casimir_resid"]),
                                 "control_residual_with_wrong_factor_3_110": fmt(e110["order3_casimir_resid_control_factor3"])}
        results["item3"][key]["parallel_resid_norm_110digits"] = fmt(e110["parallel_resid"])
        results["item3"][key]["parallel_resid_norm_160digits"] = fmt(E160["rays"][key]["parallel_resid"])
        if orth2 is not None and orth2 != 0:
            import sympy as sp
            o39 = sp.sqrt(sp.Rational(orth2.numerator, orth2.denominator) * 39)
            results["item6"][key]["orthogonal_norm_over_absg"] = f"({o39})/sqrt(39)" if o39.is_Rational else "not of the form q/sqrt(39)"
        results["item6"][key]["orth_norm_160digits"] = fmt(E160["rays"][key]["DN_orth_norm"])
results["item8"] = X["item8"]
results["checks"] = {"identification_failures": [str(f) for f in fails], "max_quadrature_vs_engine_abs_dev": maxq,
                     "exact_vs_engine110_maxdev": max(float(X["assembled"][k]["max_dev_vs_engine110"]) for k in X["assembled"]),
                     "engine_checks_60": {k: E[60][k] for k in E[60] if k.startswith("check")},
                     "engine_checks_110": {k: E[110][k] for k in E[110] if k.startswith("check")}}
print(json.dumps(results, indent=1))
print("IDENTIFICATION FAILURES:", fails)
print("max |quadrature - engine| (float64):", maxq)
(HERE / "results.json").write_text(json.dumps(results, indent=1))
