"""Assemble audit_results.json (and out/tables.md) from the computed outputs in out/."""
import json
import sympy as sp

M = json.load(open("out/main_all.json"))
G = json.load(open("out/item01.json"))
I3 = json.load(open("out/item3.json"))
I12 = json.load(open("out/item12.json"))
V = json.load(open("out/verify.json"))
RES = json.load(open("out/item1_residuals.json"))


def sq(s):
    """exact sqrt of an exact nonnegative value given as string"""
    e = sp.sympify(s.replace("*I", "*I"))
    return str(sp.nsimplify(sp.sqrt(e))) if e != 0 else "0"


def fl(s):
    return float(sp.N(sp.sympify(s), 20))


dims = {d: {int(n): int(v) for n, v in G["dims"][d].items()} for d in G["dims"]}
sector_levels = {d: [n for n in range(0, 19) if dims[d][n] > 0] for d in dims}

out = {"item0": {"norm_q1_squared": G["norm_q1_sq"], "norm_q2_squared": G["norm_q2_sq"], "order_Gamma": G["order"],
                 "derived_subgroup_order": G["derived_order"], "Gamma_equals_derived_subgroup": G["perfect"],
                 "CG_33_3m3_60": G["cg_33_3m3_60"], "CG_33_3m3_60_simplified": "1/sqrt(924)"},
       "item1": {"declaration": "DERIVED", "dim_Hom_sigma_V_(n/2)": {f"sector{d}": {str(n): dims[d][n] for n in range(19)} for d in dims},
                 "intertwiner_residuals_[intertwining_q1q2, etadag_eta_minus_I, etadag_etaprime]": RES,
                 "precision": "mpmath 50 significant digits; required residual < 1e-40"},
       "item2": {}, "item3": {}, "item4": {}, "item5": {}, "item6": {}, "item7": {}, "item8": {}, "item11": {}, "item12": {}}
names = ["U1", "U2", "U3", "U4", "U5", "U6"]
for n in names:
    for d in ("3", "4"):
        r = M[f"{n}_s{d}"]
        key = f"{n}_sector{d}"
        out["item2"][key] = {"Pi6N_is_multiple_of_Phi": r["Pi6N_multiple_of_Phi"], "multiple_Q": r["Q_fibre_ratio"],
                             "r6hat": r["r6hat"], "r6hat_stationary_on_unit_sphere": r["r6hat_stationary"]}
        lv = {}
        for lev in sector_levels[d]:
            if lev == 6:
                continue
            lv[str(lev)] = r["xi_level_norm_sq_over_g2"].get(str(lev), "0")
        out["item4"][key] = {"norm_sq_Pi_n_xi_over_g2": lv, "xi_same_character_as_Phi": True}
        e5 = {"along_Phi_per_g": r["DNxi_along_Phi"], "perp_norm_sq_per_g2": r["DNxi_perp_norm_sq"],
              "perp_norm_per_g": sq(r["DNxi_perp_norm_sq"])}
        if "DNxi_tangent_components" in r:
            e5["tangent_components_R_per_g"] = r["DNxi_tangent_components"]
        out["item5"][key] = e5
        if n in ("U5", "U6"):
            out["item6"][key] = {"form_<E,DN E>_R-Q": r["qform"], "cross_term": r["qform_cross"],
                                 "r6hat_second_derivative": r["r6hat_dd"], "r6hat_second_derivative_cross": r["r6hat_dd_cross"]}
        e7 = {"kappa_norm_sq_over_g2": r["kappa_norm_sq_over_g2"],
              "free_directions": ("span of P_T(X.Phi), X in su(2): infinitesimal rotations about the x and y axes (real dim 2); "
                                  "the z-rotation acts on Phi by a phase" if n in ("U1", "U2") else
                                  "span of P_T(X.Phi), X in su(2): the three infinitesimal rotations (real dim 3)"),
              "dim_free_R": r["dim_ker_L_on_T"],
              "perp_component_norm_sq_without_kappa_per_g2": r["perp_eq_norm_sq_without_kappa"],
              "perp_component_norm_sq_with_kappa_per_g2": r["perp_eq_norm_sq_with_kappa"]}
        if "kappa_tangent_components" in r:
            e7["kappa_over_g_components_R"] = r["kappa_tangent_components"]
            e7["remainder_norm_sq"] = r["kappa_remainder_norm_sq"]
        out["item7"][key] = e7
        out["item8"][key] = {"lambda4_over_g2_with_kappa": r["lambda4_over_g2_with_kappa"],
                             "lambda4_over_g2_without_kappa": r["lambda4_over_g2_without_kappa"],
                             "equal": r["lambda4_over_g2_with_kappa"] == r["lambda4_over_g2_without_kappa"],
                             "lambda2_over_g": r["Q_fibre_ratio"], "sign_lambda4": "negative (for every real g != 0)",
                             "decimal": fl(r["lambda4_over_g2_with_kappa"])}
for d in ("3", "4"):
    r = M[f"U5q_s{d}"]
    out["item11"][f"sector{d}"] = {"stationary_point_of_r6hat": r["r6hat_stationary"], "r6hat": r["r6hat"],
                                   "Re<Phi,DN_Phi[e_t]>": V[f"item11_Re<Phi,DN e_t>_U5q_s{d}"],
                                   "Pi6N_multiple_of_Phi": r["Pi6N_multiple_of_Phi"]}
stab = {
    "U1": ("maximal torus T = {exp(-i th Jz)} (continuous, dim 1)", "exp(-i th Jz) acts on Phi by e^{3 i th}"),
    "U2": ("N(T) = T u T.ry(pi) (normalizer of the torus)", "T acts trivially; elements of T.ry(pi) act by -1"),
    "U3": ("binary octahedral group, order 48 (image: rotations of the octahedron with vertices +-e_z, (+-1,+-1,0)/sqrt2)",
           "+-1: +1 on the binary tetrahedral subgroup (e.g. r_(sqrt2,0,1)(2pi/3), rz(pi)), -1 on the rest (rz(pi/2), rx(pi))"),
    "U4": ("binary dihedral group of order 24 (image D6, 6-fold axis z)", "rz(pi/3) -> -1, ry(pi) -> +1, rx(pi) -> -1"),
    "U5": ("cyclic group of order 10, {+-rz(2 pi k/5)} (image C5 about z)", "rz(2 pi k/5) -> e^{4 pi i k/5}"),
    "U6": ("binary dihedral group of order 12 (image D3, 3-fold axis z)", "rz(2pi/3) -> +1, rx(pi) -> -1 (sign character of D3)"),
}
for n in names:
    out["item3"][n] = {"stabilizer": stab[n][0], "character": stab[n][1],
                       "complex_dim_of_block_subspace_with_that_character": I3[n]["block_chi_subspace_dim"],
                       "chi_subspace_dims_in_V_(n/2)_by_level": I3[n]["chi_subspace_dim_by_level"]}
out["item12"] = {"max_residual": I12["worst"]["item12_residual_max"], "required": 1e-10, "arithmetic": "float64",
                 "other_quadrature_checks": I12["worst"]}
out["verification"] = V
json.dump(out, open("audit_results.json", "w"), indent=1)
print("wrote audit_results.json")

# tables for the report
L = []
for d in ("3", "4"):
    L.append(f"\n#### Sector {d}\n")
    L.append("| U | Pi6N = Q Phi | Q = lambda2/g | r6hat | r6hat stationary | lambda4/g^2 (with = without kappa) | decimal |")
    L.append("|---|---|---|---|---|---|---|")
    for n in names:
        r = M[f"{n}_s{d}"]
        L.append(f"| {n} | {r['Pi6N_multiple_of_Phi']} | {r['Q_fibre_ratio']} | {r['r6hat']} | {r['r6hat_stationary']} | "
                 f"{r['lambda4_over_g2_with_kappa']} | {fl(r['lambda4_over_g2_with_kappa']):.6e} |")
L.append("\n### item 4 table\n")
for d in ("3", "4"):
    levs = [lv for lv in sector_levels[d] if lv != 6]
    L.append(f"\n#### Sector {d}: ||Pi_n xi||^2/g^2\n")
    L.append("| U | " + " | ".join(f"n={lv}" for lv in levs) + " |")
    L.append("|---|" + "---|" * len(levs))
    for n in names:
        r = M[f"{n}_s{d}"]["xi_level_norm_sq_over_g2"]
        L.append(f"| {n} | " + " | ".join(r.get(str(lv), "**0**") for lv in levs) + " |")
L.append("\n### item 5 table\n")
for d in ("3", "4"):
    L.append(f"\n#### Sector {d}\n")
    L.append("| U | <Phi, Pi6 DN xi>/g | perp norm^2 / g^2 | perp norm / g | tangent components /g |")
    L.append("|---|---|---|---|---|")
    for n in names:
        r = M[f"{n}_s{d}"]
        L.append(f"| {n} | {r['DNxi_along_Phi']} | {r['DNxi_perp_norm_sq']} | {sq(r['DNxi_perp_norm_sq'])} | {r.get('DNxi_tangent_components', '')} |")
L.append("\n### item 7 table\n")
for d in ("3", "4"):
    L.append(f"\n#### Sector {d}\n")
    L.append("| U | ||kappa||^2/g^2 | kappa/g components | remainder | free dirs (real dim) | perp eq. norm^2/g^2 without kappa | with kappa |")
    L.append("|---|---|---|---|---|---|---|")
    for n in names:
        r = M[f"{n}_s{d}"]
        L.append(f"| {n} | {r['kappa_norm_sq_over_g2']} | {r.get('kappa_tangent_components', '-')} | {r.get('kappa_remainder_norm_sq', '-')} | "
                 f"{r['dim_ker_L_on_T']} | {r['perp_eq_norm_sq_without_kappa']} | {r['perp_eq_norm_sq_with_kappa']} |")
open("out/tables.md", "w").write("\n".join(L))
print("\n".join(L))
