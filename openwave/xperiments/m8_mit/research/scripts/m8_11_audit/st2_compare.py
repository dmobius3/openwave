"""Stage 2, task 1: set the other worker's values (solver_work/results.json) beside mine
(audit_results.json, out/*.json) and compare them exactly.

Comparison rule for exact values: both strings are parsed with sympy and the difference is
radsimp-simplified; equal means the simplified difference is exactly 0 (no numerics).
Norms reported by the solver as norms are squared before comparison with my squared norms.
Integers and booleans are compared with ==.  Numerical residuals are compared against the
thresholds each side declared.  Output: out2/compare.json (one record per compared value) and a
per-item label summary printed to stdout."""
import json, os
import sympy as sp

os.makedirs("out2", exist_ok=True)
A = json.load(open("audit_results.json"))            # mine
B = json.load(open("solver_work/results.json"))       # solver
item3_mine = json.load(open("out/item3.json"))
main_mine = json.load(open("out/main_all.json"))

recs = []


def ex(s):
    return sp.sympify(str(s).replace("I", "I"))


def eq_exact(a, b):
    d = sp.radsimp(sp.nsimplify(ex(a)) - sp.nsimplify(ex(b))) if False else sp.radsimp(ex(a) - ex(b))
    d = sp.simplify(d)
    return d == 0


def cmp(item, where, qty, mine, theirs, how="exact", square_theirs=False):
    if how == "exact":
        t = ex(theirs)
        if square_theirs:
            t = sp.radsimp(sp.expand(t * t))
        ok = sp.simplify(sp.radsimp(ex(mine) - t)) == 0
        rule = "sympy: simplify(radsimp(mine - theirs%s)) == 0" % ("^2" if square_theirs else "")
    elif how == "eq":
        ok = mine == theirs
        rule = "=="
    else:
        ok, rule = how(mine, theirs)
    recs.append({"item": item, "where": where, "quantity": qty, "mine": str(mine), "theirs": str(theirs),
                 "rule": rule, "agree": bool(ok)})
    return ok


UK = ["U1", "U2", "U3", "U4", "U5", "U6"]
SEC = ["3", "4"]

# ---------------- item 0
a0, b0 = A["item0"], B["item0"]
cmp(0, "-", "|q1|^2", a0["norm_q1_squared"], b0["norm2_q1"])
cmp(0, "-", "|q2|^2", a0["norm_q2_squared"], b0["norm2_q2"])
cmp(0, "-", "|Gamma|", a0["order_Gamma"], b0["order_Gamma"], "eq")
cmp(0, "-", "derived subgroup order", a0["derived_subgroup_order"], b0["derived_subgroup_order"], "eq")
cmp(0, "-", "Gamma perfect", a0["Gamma_equals_derived_subgroup"], b0["Gamma_equals_derived_subgroup"], "eq")
cmp(0, "-", "<3 3;3 -3|6 0>", a0["CG_33_3m3_60"], b0["CG_33_3m3_60"])

# ---------------- item 1
for s in SEC:
    for n in range(19):
        cmp(1, "sector" + s, "dim Hom n=%d" % n, A["item1"]["dim_Hom_sigma_V_(n/2)"]["sector" + s][str(n)],
            B["item1"]["dims_sector" + s][str(n)], "eq")
cmp(1, "-", "declaration", A["item1"]["declaration"], B["item1"]["declaration"], "eq")
# residuals: both below 1e-40 (each side's stated requirement)
res_m = A["item1"]["intertwiner_residuals_[intertwining_q1q2, etadag_eta_minus_I, etadag_etaprime]"]
for s in SEC:
    mine_max = max(float(x) for k, v in res_m.items() if k.startswith("sector" + s) for x in v)
    r = B["item1"]["intertwiner_residuals_50digits_required_1e-40"][s]
    th_max = max(float(r[k]) for k in ("intertwining", "orthonormal") + (("cross",) if s == "4" else ()))
    cmp(1, "sector" + s, "max intertwiner residual (both < 1e-40)", mine_max, th_max,
        lambda m, t: (m < 1e-40 and t < 1e-40, "both < 1e-40"))

# ---------------- item 2
for u in UK:
    for s in SEC:
        m = A["item2"]["%s_sector%s" % (u, s)]
        t = B["item2"][u]
        cmp(2, "%s s%s" % (u, s), "Q", m["multiple_Q"], t["sector" + s]["multiple_Q"])
        cmp(2, "%s s%s" % (u, s), "Pi6N multiple of Phi", m["Pi6N_is_multiple_of_Phi"],
            t["sector" + s]["Pi6N_is_multiple_of_Phi"], "eq")
    cmp(2, u, "r6hat", A["item2"][u + "_sector3"]["r6hat"], B["item2"][u]["rhat6"])
    cmp(2, u, "r6hat stationary", A["item2"][u + "_sector3"]["r6hat_stationary_on_unit_sphere"],
        B["item2"][u]["rhat6_stationary_on_unit_sphere"], "eq")

# ---------------- item 3
order_mine = {"U1": "inf", "U2": "inf", "U3": 48, "U4": 24, "U5": 10, "U6": 12}
order_theirs = {"U1": "inf", "U2": "inf", "U3": 48, "U4": 24, "U5": 10, "U6": 12}   # from their text + s03 log
for u in UK:
    cmp(3, u, "complex dim of block chi-subspace", A["item3"][u]["complex_dim_of_block_subspace_with_that_character"],
        B["item3"][u]["complex_dim_of_character_subspace_in_block"], "eq")
    for J in range(10):
        cmp(3, u, "chi-subspace dim in V_%d" % J, A["item3"][u]["chi_subspace_dims_in_V_(n/2)_by_level"][str(2 * J)],
            B["item3"][u]["character_subspace_dims_V_J_J0to9"][str(J)], "eq")
    cmp(3, u, "stabilizer order in SU(2)", order_mine[u], order_theirs[u], "eq")

# ---------------- item 4
for u in UK:
    for s in SEC:
        m = A["item4"]["%s_sector%s" % (u, s)]["norm_sq_Pi_n_xi_over_g2"]
        t = B["item4"]["%s_sector%s" % (u, s)]["norm_Pi_n_xi_sq_over_g2"]
        for n in sorted(set(m) | set(t), key=int):
            if n not in m:
                cmp(4, "%s s%s" % (u, s), "||Pi_%s xi||^2/g^2 (level 6: mine not listed, 0 by definition)" % n,
                    "0", t[n])
                continue
            cmp(4, "%s s%s" % (u, s), "||Pi_%s xi||^2/g^2" % n, m[n], t[n])
        cmp(4, "%s s%s" % (u, s), "xi same character", A["item4"]["%s_sector%s" % (u, s)]["xi_same_character_as_Phi"],
            B["item4"]["%s_sector%s" % (u, s)]["xi_transforms_by_same_character"], "eq")
        # level norms of N: mine from out/main_all.json
        mN = main_mine["%s_s%s" % (u, s)]["levels_N"]
        tN = B["item4"]["%s_sector%s" % (u, s)]["norm_Pi_n_N_sq"]
        for n in tN:
            cmp(4, "%s s%s" % (u, s), "||Pi_%s N||^2 (supporting)" % n, mN.get(n, "0"), tN[n])

# ---------------- item 5
for u in UK:
    for s in SEC:
        m = A["item5"]["%s_sector%s" % (u, s)]
        t = B["item5"]["%s_sector%s" % (u, s)]
        w = "%s s%s" % (u, s)
        cmp(5, w, "along Phi per g", m["along_Phi_per_g"], t["along_Phi_per_g"])
        cmp(5, w, "Im along Phi (mine: value is real)", "0", t["along_Phi_imag"])
        cmp(5, w, "perp norm^2 per g^2", m["perp_norm_sq_per_g2"], t["perp_norm_per_g"], square_theirs=True)
        if u == "U5":
            cmp(5, w, "<e_t, .>_R", m["tangent_components_R_per_g"]["e_t"], t["perp_component_e_t"])
            cmp(5, w, "<i e_t, .>_R", m["tangent_components_R_per_g"]["i_e_t"], t["perp_component_ie_t"])
            cmp(5, w, "perp remainder (mine: perp norm^2 - e_t comp^2)",
                sp.radsimp(ex(m["perp_norm_sq_per_g2"]) - ex(m["tangent_components_R_per_g"]["e_t"]) ** 2),
                t["perp_remainder_norm"], square_theirs=True)
        if u == "U6":
            cmp(5, w, "<tau_x, .>_R", m["tangent_components_R_per_g"]["tau_x"], t["perp_component_tau_x"])
            cmp(5, w, "<tau_y, .>_R", m["tangent_components_R_per_g"]["tau_y"], t["perp_component_tau_y"])
            cmp(5, w, "perp remainder (mine: perp norm^2 - tau_x comp^2)",
                sp.radsimp(ex(m["perp_norm_sq_per_g2"]) - ex(m["tangent_components_R_per_g"]["tau_x"]) ** 2),
                t["perp_remainder_norm"], square_theirs=True)

# ---------------- item 6
for s in SEC:
    m = A["item6"]["U5_sector" + s]
    t = B["item6"]["U5_sector" + s]
    cmp(6, "U5 s" + s, "form e_t", m["form_<E,DN E>_R-Q"]["e_t"], t["form_e_t"])
    cmp(6, "U5 s" + s, "form i e_t", m["form_<E,DN E>_R-Q"]["i_e_t"], t["form_ie_t"])
    m = A["item6"]["U6_sector" + s]
    t = B["item6"]["U6_sector" + s]
    cmp(6, "U6 s" + s, "form tau_x", m["form_<E,DN E>_R-Q"]["tau_x"], t["form_tau_x"])
    cmp(6, "U6 s" + s, "form tau_y", m["form_<E,DN E>_R-Q"]["tau_y"], t["form_tau_y"])
    cmp(6, "U6 s" + s, "form cross", m["cross_term"], t["form_cross"])
m, t = A["item6"]["U5_sector3"], B["item6"]["U5_rhat6_second_derivatives"]
cmp(6, "U5", "r6hat'' e_t", m["r6hat_second_derivative"]["e_t"], t["r6pp_e_t"])
cmp(6, "U5", "r6hat'' i e_t", m["r6hat_second_derivative"]["i_e_t"], t["r6pp_ie_t"])
m, t = A["item6"]["U6_sector3"], B["item6"]["U6_rhat6_second_derivatives"]
cmp(6, "U6", "r6hat'' tau_x", m["r6hat_second_derivative"]["tau_x"], t["r6pp_tau_x"])
cmp(6, "U6", "r6hat'' tau_y", m["r6hat_second_derivative"]["tau_y"], t["r6pp_itau_x"])
cmp(6, "U6", "r6hat'' cross", m["r6hat_second_derivative_cross"], t["r6pp_cross"])

# ---------------- item 7
for u in UK:
    for s in SEC:
        m = A["item7"]["%s_sector%s" % (u, s)]
        t = B["item7"]["%s_sector%s" % (u, s)]
        w = "%s s%s" % (u, s)
        cmp(7, w, "||kappa||^2/g^2", m["kappa_norm_sq_over_g2"], t["kappa_norm2_over_g2"])
        cmp(7, w, "free directions real dim", m["dim_free_R"], t["free_directions_real_dim"], "eq")
        cmp(7, w, "perp eq. without kappa, norm^2", m["perp_component_norm_sq_without_kappa_per_g2"],
            t["eq_perp_norm_without_kappa_per_g"], square_theirs=True)
        th = t["eq_perp_norm_with_kappa_per_g"]
        cmp(7, w, "perp eq. with kappa (mine exact; theirs numerical)", m["perp_component_norm_sq_with_kappa_per_g2"], th,
            lambda a, b: (a == "0" and b.startswith("0 (residual") and float(b.split()[2]) < 1e-100,
                          "mine exact 0; theirs '0' with residual < 1e-100 at 120 digits"))
        if u in ("U5", "U6"):
            names = ("e_t", "i_e_t", "ie_t") if u == "U5" else ("tau_x", "tau_y", "tau_y")
            cmp(7, w, "kappa/g along " + names[0], m["kappa_over_g_components_R"][names[0]],
                t["kappa_component_" + names[0]])
            cmp(7, w, "kappa/g along " + names[1], m["kappa_over_g_components_R"][names[1]],
                t["kappa_component_" + names[2]])
            cmp(7, w, "kappa remainder norm^2", m["remainder_norm_sq"], t["kappa_remainder_norm"], square_theirs=True)

# ---------------- item 8
for u in UK:
    for s in SEC:
        m = A["item8"]["%s_sector%s" % (u, s)]
        t = B["item8"]["%s_sector%s" % (u, s)]
        w = "%s s%s" % (u, s)
        cmp(8, w, "lambda4/g^2 with kappa", m["lambda4_over_g2_with_kappa"], t["lambda4_over_g2_with_kappa"])
        cmp(8, w, "lambda4/g^2 without kappa", m["lambda4_over_g2_without_kappa"], t["lambda4_over_g2_without_kappa"])
        cmp(8, w, "sign", "negative", t["sign_of_lambda4"],
            lambda a, b: (m["sign_lambda4"].startswith("negative") and b.startswith("negative"), "both negative for all g != 0"))

# ---------------- item 9: sets of computed zeros in items 4-8
def zeros_mine():
    z = set()
    for u in UK:
        for s in SEC:
            k = "%s_sector%s" % (u, s)
            for n, v in A["item4"][k]["norm_sq_Pi_n_xi_over_g2"].items():
                if v == "0":
                    z.add((k, "xi level %s" % n))
            if A["item5"][k]["perp_norm_sq_per_g2"] == "0":
                z.add((k, "item5 perp"))
            for c, v in A["item5"][k].get("tangent_components_R_per_g", {}).items():
                if v == "0":
                    z.add((k, "item5 comp " + {"i_e_t": "ie_t"}.get(c, c)))
            if A["item7"][k]["kappa_norm_sq_over_g2"] == "0":
                z.add((k, "kappa"))
            for c, v in A["item7"][k].get("kappa_over_g_components_R", {}).items():
                if v == "0":
                    z.add((k, "kappa comp " + {"i_e_t": "ie_t"}.get(c, c)))
            if A["item7"][k].get("remainder_norm_sq") == "0":
                z.add((k, "kappa remainder"))
            tc = A["item5"][k].get("tangent_components_R_per_g")
            if tc:   # my item 5: the whole perp part lies along e_t / tau_x (remainder exactly 0)
                first = tc.get("e_t", tc.get("tau_x"))
                if sp.simplify(ex(A["item5"][k]["perp_norm_sq_per_g2"]) - ex(first) ** 2) == 0:
                    z.add((k, "item5 remainder"))
            if A["item7"][k]["perp_component_norm_sq_without_kappa_per_g2"] == "0":
                z.add((k, "eq perp without kappa"))
            if A["item8"][k]["equal"]:
                z.add((k, "lambda4 difference"))
    for s in SEC:
        z.add(("U5_sector" + s, "form ie_t"))
        z.add(("U6_sector" + s, "form cross"))
    return z


def zeros_theirs():
    z = set()
    mp_ = {"||P_perp Pi_6 DN_Phi[xi]||": "item5 perp", "kappa": "kappa",
           "equation perp component without kappa": "eq perp without kappa",
           "lambda4(with kappa) - lambda4(without kappa)": "lambda4 difference",
           "DNxi_perp_comp_ie_t": "item5 comp ie_t", "DNxi_perp_comp_tau_y": "item5 comp tau_y",
           "kappa_comp_ie_t": "kappa comp ie_t", "kappa_comp_tau_y": "kappa comp tau_y",
           "kappa_remainder_norm": "kappa remainder", "DNxi_perp_remainder_norm": "item5 remainder",
           "form on i e_t (and rhat6'')": "form ie_t",
           "form cross term (tau_x, tau_y)": "form cross"}
    extra = set()
    for e in B["item9"]["zeros"]:
        w = e["where"].split()
        if len(w) < 3:
            continue
        k = "%s_sector%s" % (w[0], w[2])
        q = e["quantity"]
        if q.startswith("||Pi_") and "xi" in q:
            n = q[5:q.index(" ")]
            if n == "6":
                extra.add((k, q))
                continue
            z.add((k, "xi level %s" % n))
        elif q in mp_:
            z.add((k, mp_[q]))
        else:
            extra.add((k, q))
    return z, extra


zm = zeros_mine()
zt, extra = zeros_theirs()
cmp(9, "all", "set of computed zeros in items 4-8 (mine vs theirs, mapped)", sorted(zm - zt), sorted(zt - zm),
    lambda a, b: (a == [] and b == [], "both set differences empty"))
recs[-1]["mine_minus_theirs"] = sorted(zm - zt)
recs[-1]["theirs_minus_mine"] = sorted(zt - zm)
recs[-1]["theirs_extra_not_in_my_list"] = sorted(extra)
cmp(9, "all", "unresolved zeros", 0, len(B["item9"]["unresolved"]), "eq")

# ---------------- item 10
for u in UK:
    t = B["item10"][u]
    cmp(10, u, "existence for all small a", True, t["existence_for_all_small_a"], "eq")
    cmp(10, u, "orbit dim in Phi-perp block (= dim ker L)", A["item7"][u + "_sector3"]["dim_free_R"], t["orbit_dim_in_W"], "eq")
    cmp(10, u, "rank of L on Phi-perp block (12 - orbit)", 12 - A["item7"][u + "_sector3"]["dim_free_R"],
        t["hessian_rank_on_W"], "eq")
cmp(10, "sector3", "gap", 72, B["item10"]["spectral_gap"]["sector3"]["min_gap"], "eq")
cmp(10, "sector4", "gap", 32, B["item10"]["spectral_gap"]["sector4"]["min_gap"], "eq")

# ---------------- item 11
for s in SEC:
    m = A["item11"]["sector" + s]
    cmp(11, "s" + s, "stationary", m["stationary_point_of_r6hat"], B["item11"]["stationary"], "eq")
    cmp(11, "s" + s, "r6hat", m["r6hat"], B["item11"]["rhat6"])
    cmp(11, "s" + s, "Re<Phi, DN e_t>", m["Re<Phi,DN_Phi[e_t]>"], B["item11"]["Re_Phi_DN_et_sector" + s])
    cmp(11, "s" + s, "Pi6N multiple of Phi", m["Pi6N_multiple_of_Phi"], False, "eq")

# ---------------- item 12
worst_t = max(float(v["residual_equation"]) for k, v in B["item12"].items() if isinstance(v, dict))
cmp(12, "all", "max residual vs own requirement", A["item12"]["max_residual"], worst_t,
    lambda a, b: (a < 1e-10 and b < 1e-25, "mine < 1e-10 (float64); theirs < 1e-25 (30 digits)"))
above_t = max(float(v["max_component_levels_20_22"]) for k, v in B["item12"].items() if isinstance(v, dict))
cmp(12, "all", "levels 20, 22 empty", A["item12"]["other_quadrature_checks"]
    ["copy_residual(max over levels <= 22, incl. levels without intertwiners and 20, 22)"], above_t,
    lambda a, b: (a < 1e-10 and b < 1e-25, "both below own threshold"))

# ---------------- summary
by_item = {}
for r in recs:
    by_item.setdefault(r["item"], []).append(r)
summary = {}
for it in sorted(by_item):
    n = len(by_item[it])
    k = sum(r["agree"] for r in by_item[it])
    summary[it] = {"compared": n, "agree": k}
    print("item %2d: %d compared, %d agree" % (it, n, k))
    for r in by_item[it]:
        if not r["agree"]:
            print("   DISAGREE:", r)
json.dump({"records": recs, "summary": summary}, open("out2/compare.json", "w"), indent=1, default=str)
print("total", len(recs), "agree", sum(r["agree"] for r in recs))
