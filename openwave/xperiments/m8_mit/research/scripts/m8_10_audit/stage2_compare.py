"""Stage 2, task 1: set the other worker's values (solver_work/results.json, read only) beside mine,
exactly, item by item.  Also: independent exact check of its M_K -> B_K' table (item 2 reading) and of the
rank-4 claim, and the sign of the R5 orthogonal direction that it did not record.
Every comparison is a Checker line with a mutation."""
import json

import numpy as np
import sympy as sp

import lib

R = lib.ROOM
C = lib.Checker("stage2_compare")
T = json.load(open(R / "solver_work" / "results.json"))
Tx = json.load(open(R / "solver_work" / "out_exact.json"))
mine = json.load(open(R / "audit_results.json"))
alg = json.load(open(R / "algebra_primary.json"))
rays = ["R1", "R2", "R3", "R4", "R5"]
S = lambda s: sp.nsimplify(sp.sympify(s))
out = {}


def eqx(a, b):
    return sp.simplify(S(a) - S(b)) == 0


# ---------------- item 0
t0, m0 = T["item0"], mine["item0"]
rows0 = {"order": (t0["order"], m0["order"]), "derived": (t0["derived_subgroup_order"], m0["derived_subgroup_order"]),
         "perfect": (t0["equals_derived_subgroup"], m0["perfect"]),
         "norms": ((t0["norm2_q1"], t0["norm2_q2"]), ("1", "1")),
         "cg_squared": (S(t0["cg_33_3m3_60_squared"]), S(m0["cg_33_3m3_60"]) ** 2),
         "cg_sign": (t0["cg_33_3m3_60_float"] > 0, m0["cg_33_3m3_60_float"] > 0)}
ok0 = (rows0["order"][0] == rows0["order"][1] and rows0["derived"][0] == rows0["derived"][1]
       and rows0["perfect"][0] == rows0["perfect"][1] and sp.simplify(rows0["cg_squared"][0] - rows0["cg_squared"][1]) == 0
       and rows0["cg_sign"][0] == rows0["cg_sign"][1]
       and all(S(x.replace("+ 0*sqrt5", "").strip("()").split(" + ")[0]) == 1 for x in rows0["norms"][0]))
C.check("item 0: order, derived subgroup, norms, CG^2 and sign identical", lambda z: z, ok0,
        ok0 and (t0["order"] == 60), "pretend its order were 60")
out["item0"] = {k: [str(a), str(b)] for k, (a, b) in rows0.items()}

# ---------------- item 1
d1 = {d: [(T["item1"][f"{d}dim"][str(n)], mine["item1"]["homdims"][str(d)][str(n)]) for n in range(19)] for d in (3, 4)}
ok1 = all(a == b for d in d1 for a, b in d1[d])
C.check("item 1: all 38 Hom dimensions identical", lambda z: all(a == b for d in z for a, b in z[d]), d1,
        {3: d1[3], 4: [(a + (1 if i == 8 else 0), b) for i, (a, b) in enumerate(d1[4])]},
        "add 1 to its 4-dim n=8 entry")
out["item1"] = {"identical": ok1}

# ---------------- item 2
rows2 = {}
ok2 = True
for d in (3, 4):
    for K in ("0", "6"):
        theirs = T["item2"][f"{d}dim"][f"r{K}"]
        th = S(theirs.replace("sqrt(", "*sqrt(").lstrip("+").replace("-*", "-").lstrip("*"))
        mv = S(mine["item2"][str(d)]["kappa_K_exact"][K])
        rows2[f"{d}:K={K}"] = [str(th), str(mv)]
        ok2 &= sp.simplify(th - mv) == 0
C.check("item 2: coefficients of M_0, M_6 identical in both sectors (theirs sign*sqrt form)", lambda z: z, ok2,
        ok2 and rows2["3:K=6"][0] == rows2["4:K=6"][1], "compare its 3-dim r_6 with my 4-dim kappa_6")
out["item2"] = {"coefficients": rows2}

# its W table: M_K = sum_K' W_KK' B_K', B_K'(v) = [[v x v]_K' x Theta v]_3, exact check on Gaussian-integer v
W = T["item2"]["M_K_in_B_basis"]


def parseW(s):
    sign = -1 if s.startswith("-") else (0 if s.startswith("0") else 1)
    f = sp.Rational(s[s.index("(") + 1:s.index(")")])
    return sign * sp.sqrt(f)


rng = np.random.default_rng(11)
vx = [sp.Integer(int(a)) + sp.I * int(b) for a, b in rng.integers(-3, 4, (7, 2))]
tv = lib.theta(vx)
Bs = {Kp: lib.couple(lib.couple(vx, vx, Kp), tv, 3) for Kp in (0, 2, 4, 6)}
Ms = {K: lib.couple(lib.couple(vx, tv, K), vx, 3) for K in range(7)}
maxres = 0
for K in range(7):
    coeffs = [parseW(s) for s in W[str(K)]]
    rec = [sum(c * Bs[Kp][i] for c, Kp in zip(coeffs, (0, 2, 4, 6))) for i in range(7)]
    res = [sp.expand(a - b) for a, b in zip(Ms[K], rec)]
    maxres = max(maxres, max(float(abs(sp.N(x, 40))) for x in res))
    if K == 6:
        mut_res = max(float(abs(sp.N(sp.expand(a - b - Bs[6][i] / 10), 40))) for i, (a, b) in enumerate(zip(Ms[K], rec)))
C.check(f"its W table reproduces every M_K from B_0,B_2,B_4,B_6 with my CG, exact v (max |res| {maxres:.1e})",
        lambda z: z < 1e-30, maxres, mut_res, "add B_6/10 to its M_6 row")
# rank of span{M_K}: numeric on 8 random complex vectors
Mrows = []
for K in range(7):
    row = []
    for s in range(8):
        v = [complex(a, b) for a, b in np.random.default_rng(100 + s).normal(size=(7, 2))]
        vs_ = [sp.Float(z.real, 30) + sp.I * sp.Float(z.imag, 30) for z in v]
        mk = lib.couple(lib.couple(vs_, lib.theta(vs_), K), vs_, 3)
        row += [complex(sp.N(x, 20)) for x in mk]
    Mrows.append(row)
Mm = np.array(Mrows)
svals = np.linalg.svd(np.vstack([Mm.real, Mm.imag]).reshape(7, -1) if False else np.hstack([Mm.real, Mm.imag]),
                      compute_uv=False)
rank = int(np.sum(svals > 1e-8 * svals[0]))
C.check(f"rank of span{{M_0..M_6}} = 4 (singular values {', '.join(f'{x:.3g}' for x in svals)})",
        lambda r: r == 4, rank, rank + 3, "claim full rank 7")
out["item2"]["W_table_verified_max_res"] = maxres
out["item2"]["rank_M_K"] = rank
out["item2"]["singular_values"] = [float(x) for x in svals]

# ---------------- items 3, 4, 5, 6, 7
rows = {}
for d in ("3", "4"):
    for r in rays:
        key = f"{d}:{r}"
        a = alg[d]["rays"][r]
        t3, t4, t6, t7 = T["item3"][key], T["item4"][key], T["item6"][key], T["item7"][key]
        rec = {}
        rec["item3"] = (t3["parallel"] == a["item3_is_multiple"]) and eqx(t3["multiple"], a["item3_multiple"])
        lv = t4["norm2_over_g2_by_level"]
        rec["item4"] = all(eqx(x, a["item4_normPi_xi2_over_g2_by_n"][n]) for n, x in lv.items())
        nonsector = [n for n in a["item4_normPi_xi2_over_g2_by_n"] if int(n) not in t4["sector_levels"]]
        rec["item4_nonsector_zero_mine"] = all(S(a["item4_normPi_xi2_over_g2_by_n"][n]) == 0 for n in nonsector)
        rec["item4_levels_compared"] = sorted(int(n) for n in lv)
        my_zeros = sorted(int(n) for n, x in a["item4_normPi_xi2_over_g2_by_n"].items()
                          if int(n) in t4["sector_levels"] and int(n) != 6 and S(x) == 0)
        rec["item5"] = sorted(T["item5"][key]) == my_zeros
        rec["item5_zeros"] = my_zeros
        rec["item6_along"] = eqx(t6["along_Phi_coefficient_over_g"], a["item6_along"])
        rec["item6_orth2"] = eqx(t6["orthogonal_norm2_over_g2"], a["item6_perp_norm2"])
        rec["item7"] = eqx(t7["lambda4_over_g2"], a["item7_lambda4_over_g2"])
        rec["item7_consistent_flag"] = t7["a5_block_equation_consistent"] == (S(a["item6_perp_norm2"]) == 0)
        rows[key] = rec
for it in ("item3", "item4", "item4_nonsector_zero_mine", "item5", "item6_along", "item6_orth2", "item7",
           "item7_consistent_flag"):
    vals = {k: rows[k][it] for k in rows}
    C.check(f"{it}: identical in all 10 (sector, ray) cases", lambda z: all(z.values()), vals,
            {**vals, "3:R1": False}, "flip the 3-dim R1 comparison to False")
# exact-value mutation: the comparisons are not vacuous (other sector's value must differ)
mutcount = sum(eqx(T["item7"][f"3:{r}"]["lambda4_over_g2"], alg["4"]["rays"][r]["item7_lambda4_over_g2"]) for r in rays)
C.check("item 7 comparison discriminates (its 3-dim values never equal my 4-dim values)", lambda z: z == 0,
        mutcount, 5, "pretend all five matched")
out["items3to7"] = {k: {kk: (vv if not isinstance(vv, bool) else vv) for kk, vv in v.items()} for k, v in rows.items()}

# its R5 orthogonal norm given as (p)/sqrt(39): exact check against my exact norm^2
for d in ("3", "4"):
    s = T["item6"][f"{d}:R5"]["orthogonal_norm_over_absg"]
    p = S(s[1:s.index(")")])
    ok = sp.simplify(p ** 2 / 39 - S(alg[d]["rays"]["R5"]["item6_perp_norm2"])) == 0
    C.check(f"its R5 {d}-dim orthogonal norm {s} squared == my exact norm^2", lambda z: z, ok,
            sp.simplify(p ** 2 / 40 - S(alg[d]["rays"]["R5"]["item6_perp_norm2"])) == 0, "use sqrt(40) instead of sqrt(39)")

# sign of the R5 orthogonal fibre vector along e = sin t v_2 - cos t v_-3 (it did not record it)
st, ct = 2 * sp.sqrt(3) / 5, sp.sqrt(13) / 5
edir = [0] * 7
edir[5], edir[0] = st, -ct
signs = {}
for d in ("3", "4"):
    up = [S(x) for x in alg[d]["rays"]["R5"]["item6_perp_fibre"]]
    coef = sp.nsimplify(sp.radsimp(sum(a * b for a, b in zip(edir, up))))
    par = sp.simplify(sum((b - coef * a) ** 2 for a, b in zip(edir, up)))
    signs[d] = {"coefficient_along_e": str(coef), "float": float(coef), "residual_perp_to_e": str(par)}
    C.check(f"R5 {d}-dim: my orthogonal fibre vector is exactly parallel to sin t v_2 - cos t v_-3 (coef {float(coef):.6e})",
            lambda z: z == 0, par, sp.simplify(sum((b - coef * a - a / 7) ** 2 for a, b in zip(edir, up))),
            "subtract a wrong multiple")
out["R5_orthogonal_sign"] = signs

# ---------------- item 8
ok8 = all(eqx(T["item8"][r]["lambda4_ratio_3_over_4"], mine["item8"][r]["lambda4_ratio_3_over_4"])
          and eqx(T["item8"][r]["kappa_minus_1_ratio_3_over_4"], mine["item8"][r]["multiple_minus_1_ratio"]) for r in rays)
C.check("item 8: both ratio columns identical at all five rays", lambda z: z, ok8,
        ok8 and eqx(T["item8"]["R1"]["lambda4_ratio_3_over_4"], mine["item8"]["R2"]["lambda4_ratio_3_over_4"]),
        "compare its R1 ratio with my R2 ratio")
out["item8"] = {"identical": ok8}

# ---------------- item 9: its residual numbers (method differs; recorded only)
out["item9_theirs"] = T["item9"]
out["their_identification_failures"] = T["checks"]["identification_failures"]
out["checks"] = C.records
out["all_ok"] = C.all_ok()
json.dump(out, open(R / "stage2_compare.json", "w"), indent=1, default=str)
print(json.dumps(signs, indent=1))
print("ALL OK" if C.all_ok() else "SOME CHECKS FAILED")
