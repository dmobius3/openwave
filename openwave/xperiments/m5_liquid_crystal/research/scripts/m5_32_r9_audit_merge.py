"""Merge the M5.32 R9 audit checkpoints into ../data/m5_32_r9_audit.json.

Reads only ../checkpoints/m5_32_r9/*.json written by m5_32_r9_audit.py.
One entry per claim id (A1..A4, B1..B5) with verdict / own_number /
producer_number / rel_dev / note, plus summary and new_findings.
"""
from __future__ import annotations

import json
import os
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
CKPT = os.path.join(HERE, "..", "checkpoints", "m5_32_r9")


def load(tag):
    p = os.path.join(CKPT, f"{tag}.json")
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def rel(a, b):
    if a is None or b is None:
        return None
    d = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / d


def main():
    t0 = time.time()
    tube = load("tube")
    ring = load("ring")
    kin = load("kin")
    topo = load("topo")
    share = load("share")
    mut = load("mutate")
    dv = load("divergence")
    rctrl = load("relax_control")
    relax = load("relax") or load("relax_partial")
    probe = load("probe")
    ex = tube["exponents"]
    vd = tube["verdicts"]
    cuts = tube["rho_cuts"]
    ic = cuts.index(3.0)
    out = {}

    # ------------------------------ A1 ------------------------------
    h_z = vd["Q_I1sq"]["cyl_z"]["C4_h_exponent_vs_rho_cut"][ic]
    h_sh = vd["Q_I1sq"]["cyl_shift"]["C4_h_exponent_vs_rho_cut"][ic]
    L_z = vd["Q_I1sq"]["cyl_z"]["C4_L_exponent_vs_rho_cut"]
    ctrl = {fn: vd["Q_I1sq"][fn]["C4_h_exponent_vs_rho_cut"][ic]
            for fn in ("cyl_x", "cyl_off", "sph", "rand")}
    out["A1"] = {
        "verdict": "QUALIFIED",
        "own_number": h_z,
        "producer_number": 0.251,
        "rel_dev": rel(h_z, 0.251),
        "note": ("h-stability CONFIRMED and better than claimed: correctly "
                 "centered cut gives Q_I1sq C4 h-exponent +0.019 (uncut "
                 "-1.611), value stable to 1.3 pct across h; all four "
                 "volume-matched controls stay near -1.6, so the cure is "
                 "string-specific. IR-convergence is SUPPORTED by successive "
                 "increments that halve (0.0575 then 0.0288) but the stated "
                 "diagnostic is not self-consistent: the L-exponent drifts "
                 "0.115 to 0.400 as rho_cut runs 1.5 to 9."),
        "detail": {
            "IR_increment_test": {
                "definition": ("relative increment of C4 across the L ladder "
                               "48 -> 72 -> 96; a halving increment is a "
                               "converging sequence, a flat one is a power "
                               "law"),
                "cyl_z_rho3": [0.0575, 0.0288],
                "cyl_z_rho9": [0.2081, 0.0883],
                "uncut": [0.2671, 0.2055],
                "sphere_control_rho3": [0.4234, 0.2999],
                "reading": ("after the axis cut the increments halve at every "
                            "rho_cut, so the excised C4 converges; the uncut "
                            "and the volume-matched sphere do not")},
            "C4_absolute_L_ladder_cyl_z": {
                f"rho_{c:g}": [
                    tube["boxes"][k]["terms"]["Q_I1sq"][f"cyl_z_{c:g}"]["C4"]
                    for k in ("n32_L48", "n48_L72", "n64_L96")]
                for c in cuts},
            "C4_h_exponent_rho3_true_center": h_z,
            "C4_h_exponent_rho3_producer_mask": h_sh,
            "C4_h_exponent_uncut": vd["Q_I1sq"]["cyl_z"][
                "C4_h_exponent_vs_rho_cut"][0],
            "C4_h_exponent_volume_matched_controls_rho3": ctrl,
            "C4_L_exponent_vs_rho_cut": L_z,
            "rho_cuts": cuts,
            "Q_Fpair_C4_h_exponent_rho3_true_center":
                vd["Q_Fpair"]["cyl_z"]["C4_h_exponent_vs_rho_cut"][ic],
            "Q_Fpair_producer_number": 0.257,
            "C4_value_h1.5_vs_h0.75_rho3": [
                tube["boxes"]["n32_L48"]["terms"]["Q_I1sq"]["cyl_z_3"]["C4"],
                tube["boxes"]["n64_L48"]["terms"]["Q_I1sq"]["cyl_z_3"]["C4"]],
            "C4_value_h1.5_vs_h0.75_uncut": [
                tube["boxes"]["n32_L48"]["terms"]["Q_I1sq"]["cyl_z_0"]["C4"],
                tube["boxes"]["n64_L48"]["terms"]["Q_I1sq"]["cyl_z_0"]["C4"]],
            "producer_mask_defect": (
                "the producer's mask is built on a cell-centered grid "
                "x = (i - n//2) h while the density lives on the certified "
                "offset grid x = (i - (n-1)/2) h, so their cylinder is "
                "centered at (h/2, h/2), an axis offset of h/sqrt(2) that "
                "MOVES with the lattice. Reproduced to 6 significant figures "
                "at every rho_cut by the cyl_shift family."),
            "mutation_string_moved": (mut or {}).get("M2b_string_moved", {}),
            "point_core_vs_line_core": dv}}

    # ------------------------------ A2 ------------------------------
    a2 = vd["Q_I4sq"]["cyl_z"]["C4_h_exponent_vs_rho_cut"][ic]
    out["A2"] = {
        "verdict": "CONFIRMED",
        "own_number": a2,
        "producer_number": 5.653,
        "rel_dev": rel(a2, 5.653),
        "note": ("Q_I4sq stays undecidable: C4 h-exponent +4.27 at rho >= 3 "
                 "with the correct center (+5.62 reproducing the producer's "
                 "own mask, matching 5.653), and its L-exponent even changes "
                 "sign, 1.034 to -1.346, as rho_cut grows."),
        "detail": {
            "C4_h_exponent_vs_rho_cut_true_center":
                vd["Q_I4sq"]["cyl_z"]["C4_h_exponent_vs_rho_cut"],
            "C4_h_exponent_vs_rho_cut_producer_mask":
                vd["Q_I4sq"]["cyl_shift"]["C4_h_exponent_vs_rho_cut"],
            "C4_L_exponent_vs_rho_cut":
                vd["Q_I4sq"]["cyl_z"]["C4_L_exponent_vs_rho_cut"],
            "rho_cuts": cuts}}

    # ------------------------------ A3 ------------------------------
    hL = [ex["h"]["I1"][f"cyl_z_{c:g}"]["C2"][0] for c in cuts]
    LL = [ex["L"]["I1"][f"cyl_z_{c:g}"]["C2"] for c in cuts]
    LS = [ex["L"]["I1"][f"sph_{c:g}"]["C2"] for c in cuts]
    b = tube["boxes"]
    fullc2 = b["n32_L48"]["terms"]["I1"]["cyl_z_0"]["C2"]
    sh_z = (fullc2 - b["n32_L48"]["terms"]["I1"]["cyl_z_3"]["C2"]) / fullc2
    sh_s = (fullc2 - b["n32_L48"]["terms"]["I1"]["sph_3"]["C2"]) / fullc2
    out["A3"] = {
        "verdict": "CONFIRMED",
        "own_number": max(abs(v) for v in hL if v is not None),
        "producer_number": 0.131,
        "rel_dev": rel(max(abs(v) for v in hL if v is not None), 0.131),
        "note": ("the certified inertia is string-independent: kin's C2 "
                 "h-exponent stays within 0.062 of zero at every rho_cut, and "
                 "the L-exponent drift 1.063 to 1.443 is reproduced by the "
                 "equal-volume SPHERICAL control (1.063 to 1.461), so it is a "
                 "volume effect, not a string effect."),
        "detail": {
            "C2_h_exponent_vs_rho_cut": hL,
            "C2_L_exponent_vs_rho_cut_axis": LL,
            "C2_L_exponent_vs_rho_cut_sphere_control": LS,
            "rho_cuts": cuts,
            "kin_h1.5_vs_h0.75_at_L48": [
                -4.0 * b["n32_L48"]["terms"]["I1"]["cyl_z_0"]["C2"],
                -4.0 * b["n64_L48"]["terms"]["I1"]["cyl_z_0"]["C2"]],
            "kin_tube_share_rho_lt_3": sh_z,
            "kin_equal_volume_sphere_share_rho_lt_3": sh_s,
            "axis_excess_over_sphere_control": sh_z - sh_s}}

    # ------------------------------ A4 ------------------------------
    m3 = (mut or {}).get("M3_C6a_derivative_free", {})
    out["A4"] = {
        "verdict": "CONFIRMED",
        "own_number": m3.get("rel_dev_real_vs_closed"),
        "producer_number": 3.000,
        "rel_dev": rel(ex["L"]["Q_C6a"]["cyl_z_0"]["C4"], 3.000),
        "note": ("true but uninformative: C4(Q_C6a) = sum (tr(a0 eta a0 "
                 "eta))^2 h^3 contains NO spatial derivative (identity holds "
                 "to 1.3e-15), so L^3 and h^0 are automatic; the axis cut and "
                 "the equal-volume spherical cut give IDENTICAL exponents, so "
                 "the quantity carries no frame information at all."),
        "detail": {
            "C4_L_exponent_axis": [ex["L"]["Q_C6a"][f"cyl_z_{c:g}"]["C4"]
                                   for c in cuts],
            "C4_L_exponent_sphere": [ex["L"]["Q_C6a"][f"sph_{c:g}"]["C4"]
                                     for c in cuts],
            "C4_h_exponent_axis": [ex["h"]["Q_C6a"][f"cyl_z_{c:g}"]["C4"][0]
                                   for c in cuts],
            "derivative_free_identity": m3}}

    # ------------------------------ B1 ------------------------------
    b1 = ring["B1_frame_free_form"]
    out["B1"] = {
        "verdict": "CONFIRMED",
        "own_number": b1["delta_0"]["rel_to_field_scale_32"],
        "producer_number": 2.082e-17,
        "rel_dev": rel(b1["delta_0"]["rel_to_field_scale_32"], 2.082e-17),
        "note": ("identity holds at 2.0817e-17 relative. Derived: the delta > "
                 "0 residual of the SAME form is exactly delta (sup over "
                 "directions of delta |u2_i u2_j|), so 0.2998 is delta = 0.3 "
                 "under-sampled, not a separate number; the BEST frame-free "
                 "fit uses lam_perp = delta/2 and leaves exactly delta/2."),
        "detail": {
            "rows": b1,
            "residual_over_delta_at_delta_0.3":
                b1["delta_0.3"]["dev_over_delta"],
            "best_fit_residual_over_delta":
                b1["delta_0.3"]["best_dev_over_delta"],
            "units_caveat": ("the producer compares a RELATIVE number "
                             "(2.082e-17) with an ABSOLUTE one (0.2998); the "
                             "delta = 0.3 deviation relative to the field "
                             "scale 32 is 9.375e-3"),
            "derivation": ring["B1_derivation"]}}

    # ------------------------------ B2 ------------------------------
    r = ring["B2_rings"]
    out["B2"] = {
        "verdict": "CONFIRMED",
        "own_number": r["delta_0.3_z_12"]["spread_over_delta"],
        "producer_number": 0.5000,
        "rel_dev": rel(r["delta_0.3_z_12"]["spread_over_delta"], 0.5),
        "note": ("spread / delta = 0.5000000000000002 at rho = 1e-9, at z = "
                 "+12 and -12 and 0.5, for every delta > 0; at delta = 0 it "
                 "falls exactly as rho (ratio 1e-7 over 7 decades). The z "
                 "axis is the ONLY singular locus: x-axis rings shrink "
                 "linearly. QUALIFIER: at z = 0 the law is (1-delta)/2, which "
                 "does NOT vanish at delta = 0 (it is 0.5) and DOES vanish at "
                 "delta = 1; that ring measures the origin point core. A "
                 "GLOBAL box scan confirms the locus exhaustively: only the "
                 "rho < 1.5 bin keeps its neighbour jump when h halves (ratio "
                 "1.000), every other bin halves (0.679 to 0.511)."),
        "detail": {
            "spread_over_delta_z12": {
                k: v["spread_over_delta"] for k, v in r.items()
                if k.endswith("_z_12")},
            "delta_0_ratio_over_7_decades": r["delta_0_z_12"]["ratio_1e2_to_1e9"],
            "z_equals_0_rows": {k: v for k, v in r.items() if k.endswith("_z_0")},
            "other_axis": ring["B2_other_axis"],
            "singular_locus_scan": (probe or {}).get("singular_locus_scan"),
            "global_jump_scan": ring.get("B2_global_jump_scan"),
            "derivation": ring["B2_derivation"]}}

    # ------------------------------ B3 ------------------------------
    dr = kin["kin_over_delta2"]["drift_max_over_min"]
    out["B3"] = {
        "verdict": "QUALIFIED",
        "own_number": dr,
        "producer_number": 1.0,
        "rel_dev": rel(dr, 1.0),
        "note": ("kin = 0.0 EXACTLY at delta = 0 (a0 is identically zero), "
                 "confirmed. The delta^2 law is ASYMPTOTIC only: the local "
                 "exponent runs 1.9998 (delta 1e-4 to 3e-4) down to 1.946 "
                 "(0.1 to 0.3) and up to 2.650 (0.3 to 1), and kin/delta^2 "
                 "drifts by a factor 2.19 across the ladder."),
        "detail": {
            "kin_at_delta_0": kin["kin_at_delta_0"],
            "a0_at_delta_0_max_abs": kin["a0_at_delta_0_max_abs"],
            "rows": kin["rows"], "local_exponents": kin["local_exponents"],
            "kin_over_delta2": kin["kin_over_delta2"],
            "mutation": kin["mutation"]}}

    # ------------------------------ B4 ------------------------------
    cex = topo["smooth_orbit_counterexample"]
    kin_relaxed = ([rw["kin_frozen_a0_after"] for rw in relax["rows"]]
                   if relax and relax.get("rows") else None)
    b4 = {
        "verdict": "REFUTED",
        "own_number": (kin_relaxed[0] if kin_relaxed else
                       cex["rows"][0]["kin_radial_boost"]),
        "producer_number": 0.0,
        "rel_dev": 1.0,
        "note": ("broken twice. (i) At delta = 0 the radial boost K_1 has a "
                 "SMOOTH orbit and a0 = 33 sym(e0 n-hat) is nonzero, kin = "
                 "-6.094e6 and h-stable. (ii) FIRE resolves the line into a "
                 "finite core of physical radius ~3.6 to 4.0 (h-independent "
                 "to 10 pct) and the clock SURVIVES, kin = 351.17 at h = 1.5 "
                 "and 351.14 at h = 0.75, an h-convergence the rigid ansatz "
                 "never had (426.51 vs 445.22). A string-free hedgehog and a "
                 "nonzero clock therefore COEXIST. What survives is the "
                 "narrower statement: inside the eigenvalue-pinned family no "
                 "PERIODIC clock has a smooth orbit, which I verified "
                 "exhaustively over the whole algebra. The defect is physical "
                 "and computable, not an artifact: it is the class -1 of "
                 "pi_1(OPM) = Q8."),
        "detail": {
            "attack_a_assignment_enumeration": {
                "n_rotational_smooth_and_nonzero": topo[
                    "assignment_enumeration"]["n_rotational_smooth_and_nonzero"],
                "n_boost_smooth_and_nonzero": topo[
                    "assignment_enumeration"]["n_boost_smooth_and_nonzero"],
                "reading": topo["assignment_enumeration"]["reading"]},
            "attack_a_exhaustive_algebra": topo["attack_a_exhaustive_algebra"],
            "certified_catalog_cross_check": topo[
                "certified_catalog_cross_check"],
            "attack_b_rp1_lift": topo["line_field_lift"],
            "attack_b_poincare_hopf": topo["poincare_hopf"],
            "topology_of_the_order_parameter": {
                "discrete_stabilizer": topo["discrete_stabilizer"],
                "pi1": topo["pi1_preimage"],
                "line_class": topo["ansatz_line_class"]},
            "smooth_orbit_counterexample": cex,
            "stabilizer_tangents": topo["stabilizer_tangents"]}}
    if rctrl:
        gz = rctrl["dist_to_z_axis"]["gap_after"]
        gx = rctrl["dist_to_x_axis"]["gap_after"]
        b4["detail"]["attack_c_controls"] = {
            "dt_independence": {
                "dt0_0.01_E_u_end": (relax["rows"][0]["E_end"][0]
                                     if relax and relax.get("rows") else None),
                "dt0_0.02_E_u_end": rctrl["E_end"][0],
                "dt0_0.01_V4_end": (relax["rows"][0]["E_end"][1]
                                    if relax and relax.get("rows") else None),
                "dt0_0.02_V4_end": rctrl["E_end"][1],
                "reading": ("the two FIRE step sizes land on the same "
                            "endpoint, so the melt is not a stepper "
                            "artifact")},
            "axis_localization": {
                "rho": rctrl["dist_to_z_axis"]["rho"],
                "gap_vs_distance_to_z_axis": gz,
                "gap_vs_distance_to_x_axis": gx,
                "innermost_z": gz[0], "innermost_x": gx[0],
                "ratio": gx[0] / gz[0],
                "reading": ("measured on the SAME relaxed field: the biaxial "
                            "gap collapses to 0.166 next to the z axis but "
                            "stays flat at 0.294, the far-field value, next "
                            "to the non-singular x axis. The relaxation "
                            "resolves THE LINE, it does not soften the field "
                            "generically.")}}
    if relax and relax.get("rows"):
        b4["detail"]["attack_c_relaxation"] = relax
        rows = relax["rows"]

        def width_at(rw, frac=0.9):
            """interpolated rho where the biaxial gap reaches `frac` of its
            far-field value: the resolved core radius."""
            far = rw["gap_far_after"]
            mid, gap = rw["rho_mid"], rw["gap_after"]
            for i in range(1, len(mid)):
                if gap[i] >= frac * far:
                    g0, g1 = gap[i - 1], gap[i]
                    if g1 == g0:
                        return float(mid[i])
                    t = (frac * far - g0) / (g1 - g0)
                    return float(mid[i - 1] + t * (mid[i] - mid[i - 1]))
            return None
        b4["detail"]["attack_c_core_width"] = {
            "definition": "rho where the biaxial gap reaches 90 pct of far field",
            "h": [rw["h"] for rw in rows],
            "core_radius": [width_at(rw) for rw in rows],
            "core_radius_over_h": [
                (width_at(rw) / rw["h"] if width_at(rw) else None)
                for rw in rows],
            "gap_on_axis_over_far": [
                rw["gap_after"][0] / rw["gap_far_after"] for rw in rows],
            "verdict_rule": ("a core radius that is the same PHYSICAL length "
                             "at both h is a resolved physical core; one that "
                             "tracks h is a lattice artifact"),
            "core_radius_ratio_h0.75_over_h1.5": (
                width_at(rows[1]) / width_at(rows[0])
                if len(rows) > 1 and width_at(rows[0]) else None)}
        # THE smoothness test: a discontinuity keeps its near-axis ring
        # spread when h halves; a resolved core lets it shrink.
        b4["detail"]["attack_c_smoothness"] = {
            "h": [rw["h"] for rw in rows],
            "ring_spread_before_relax": [rw["ring_spread_before"] for rw in rows],
            "ring_spread_after_relax": [rw["ring_spread_after"] for rw in rows],
            "ratio_before": (rows[1]["ring_spread_before"]
                             / rows[0]["ring_spread_before"]
                             if len(rows) > 1 else None),
            "ratio_after": (rows[1]["ring_spread_after"]
                            / rows[0]["ring_spread_after"]
                            if len(rows) > 1 else None),
            "reading": ("BEFORE relaxation the near-axis ring spread does not "
                        "shrink when h halves (ratio 1.031): a genuine "
                        "discontinuity. AFTER relaxation it shrinks (ratio "
                        "0.702): the line has been resolved into a smooth "
                        "finite core. This is the test that decides attack "
                        "(c)."),
            "clock_h_convergence": {
                "kin_before_relax": [rw["kin_frozen_a0_before"] for rw in rows],
                "kin_after_relax": [rw["kin_frozen_a0_after"] for rw in rows],
                "reading": ("the rigid ansatz's clock inertia is h-dependent "
                            "(426.51 at h = 1.5, 445.22 at h = 0.75); after "
                            "the core resolves it becomes h-convergent "
                            "(351.17 and 351.14, 0.01 pct apart). Resolving "
                            "the defect does not kill the clock, it makes it "
                            "well defined.")}}
        b4["detail"]["attack_c_summary"] = {
            "gap_before_on_axis": [rw["gap_before"][0] for rw in rows],
            "gap_after_on_axis": [rw["gap_after"][0] for rw in rows],
            "gap_far_after": [rw["gap_far_after"] for rw in rows],
            "core_width_after": [rw["core_width_after"] for rw in rows],
            "h": [rw["h"] for rw in rows],
            "ring_spread_before": [rw["ring_spread_before"] for rw in rows],
            "ring_spread_after": [rw["ring_spread_after"] for rw in rows],
            "E_drop": [rw["E_drop"] for rw in rows],
            "V4_start": [rw["E_start"][1] for rw in rows],
            "V4_end": [rw["E_end"][1] for rw in rows],
            "leading_eigenvector_tilt_on_axis_after": [
                rw["tilt_after"][0] for rw in rows],
            "kin_frozen_a0_before": [rw["kin_frozen_a0_before"] for rw in rows],
            "kin_frozen_a0_after": [rw["kin_frozen_a0_after"] for rw in rows],
            "converged": [rw["stop"] for rw in rows],
            "mechanism": (
                "the relaxation does NOT escape into the third dimension (the "
                "leading spatial eigenvector stays aligned with n-hat to "
                "1e-3); it MELTS the biaxiality, driving the (delta, 0) "
                "eigenvalue pair toward degeneracy on the axis. V4 rises from "
                "1e-23 to a finite value, which is exactly the field leaving "
                "the eigenvalue-pinned manifold that B4's theorem assumes.")}
    out["B4"] = b4

    # ------------------------------ B5 ------------------------------
    s = share["rows"]
    exc = {f"delta_{d:g}": s[f"n32_delta{d:g}"]["Q_I1sq"]["C4_axis_excess_3"]
           for d in (0.003, 0.01, 0.03, 0.1, 0.3, 1.0)}
    out["B5"] = {
        "verdict": "REFUTED",
        "own_number": exc["delta_0.003"],
        "producer_number": 0.580,
        "rel_dev": rel(exc["delta_0.003"], 0.580),
        "note": ("REFUTED for delta <= 0.1. At delta = 0.003 the z-axis tube "
                 "(rho < 3) holds 0.5796 of Q_I1sq's C4 but the EQUAL-VOLUME "
                 "x-axis tube, which contains no string, holds 0.5868: the "
                 "null control captures MORE, so the axis excess is -0.0072. "
                 "The 0.58 share is the ORIGIN point core, not the line. The "
                 "claim holds only at delta >= 0.3 (excess +0.503) and delta "
                 "= 1 (+0.937)."),
        "detail": {
            "axis_excess_over_equal_volume_x_tube_rho3": exc,
            "z_tube_share_rho3": {
                f"delta_{d:g}": s[f"n32_delta{d:g}"]["Q_I1sq"][
                    "C4_cyl_z_share_lt_3"]
                for d in (0.003, 0.01, 0.03, 0.1, 0.3, 1.0)},
            "x_tube_share_rho3": {
                f"delta_{d:g}": s[f"n32_delta{d:g}"]["Q_I1sq"][
                    "C4_cyl_x_share_lt_3"]
                for d in (0.003, 0.01, 0.03, 0.1, 0.3, 1.0)},
            "sphere_share_rho3": {
                f"delta_{d:g}": s[f"n32_delta{d:g}"]["Q_I1sq"][
                    "C4_sph_share_lt_3"]
                for d in (0.003, 0.01, 0.03, 0.1, 0.3, 1.0)},
            "h_dependence_note": ("the producer's tube is rho < 1.01 h, a "
                                  "LATTICE region; at a fixed PHYSICAL tube "
                                  "the share still moves strongly with h "
                                  "(0.5796 at h = 1.5 to 0.8000 at h = 0.75 "
                                  "at delta = 0.003), and so does the null "
                                  "control, confirming an origin effect"),
            "rows": s}}

    verd = [out[k]["verdict"] for k in out]
    out["summary"] = (
        "Arm a's excision is a real and correctly targeted diagnostic: the "
        "h-divergence of the C5 omega^4 coefficient is cured ONLY by a cut on "
        "the z axis, and four volume-matched controls (x-axis tube, offset "
        "tube, spherical core, random cells) leave it at -1.6, while moving "
        "the hedgehog to x = 12 moves the cure with it. Two defects: the "
        "producer's mask is off-center by (h/2, h/2) so their stated cut is "
        "not the cut they ran and its center moves with the lattice, and the "
        "IR (L) exponent is not rho_cut independent, so 'IR-convergent' is "
        "not established. Arm b's theorem survives its two TOPOLOGICAL "
        "attacks and is in fact stronger than stated there: the z-axis line "
        "is the class -1 of pi_1(SO(1,3)+/D_2) = Q8, a topologically "
        "protected biaxial disclination rather than a coordinate artifact, "
        "while pi_2 = 0 means the hedgehog carries no point charge at all. "
        "But the headline B4 fails on the third attack, which is the one that "
        "decides the rung: relaxation resolves the line into a finite "
        "PHYSICAL core (radius 3.98 at h = 1.5, 3.58 at h = 0.75) by melting "
        "the biaxiality, the near-axis discontinuity disappears (ring spread "
        "ratio 1.031 before, 0.702 after), and the clock SURVIVES with an "
        "h-convergent inertia (351.17 and 351.14, against 426.51 and 445.22 "
        "rigid). A string-free hedgehog and a nonzero clock coexist, so the "
        "run's conclusion should be 'the defect is physical and must be "
        "resolved', not 'the ansatz is broken'. The exclusion survives only "
        "inside the rigid eigenvalue-pinned family, where I closed it "
        "exhaustively over the whole algebra. B5 is refuted below delta = "
        "0.3, where the quartic weight sits at the origin core, not the line.")
    out["n_confirmed"] = sum(v == "CONFIRMED" for v in verd)
    out["n_qualified"] = sum(v == "QUALIFIED" for v in verd)
    out["n_refuted"] = sum(v == "REFUTED" for v in verd)
    out["mutation_gates"] = {
        k: v.get("fires") for k, v in (mut or {}).items() if k.startswith("M")}
    out["all_mutations_fire"] = (mut or {}).get("all_mutations_fire")
    out["new_findings"] = [
        "PRODUCER MASK OFF-CENTER: arm a's excision mask is built on a "
        "cell-centered grid while the density lives on the certified offset "
        "grid, so the cylinder is centered at (h/2, h/2), an axis offset of "
        "h/sqrt(2) that MOVES with the lattice. Reproduced to 6 significant "
        "figures at every rho_cut. At rho_cut = 1.5 it removes 1 of the 4 "
        "innermost string columns. The correctly centered cut gives a BETTER "
        "h-exponent (+0.019 vs +0.251), so the conclusion strengthens, but "
        "the h-ladder compares two differently centered regions.",
        "THE LINE IS TOPOLOGICALLY PROTECTED, NOT A COORDINATE ARTIFACT: the "
        "stabilizer of d4 in SO(1,3)+ is the Klein four-group (order 4, "
        "verified by enumeration), so the order-parameter manifold has "
        "pi_1 = Q8 and the ansatz's 2 pi loop lifts to -1 in SU(2) (endpoint "
        "trace -2.0), a nontrivial class. The defect therefore cannot be "
        "combed away and must be resolved with a melted core; pi_2 = 0 "
        "additionally means the hedgehog carries NO topological point charge.",
        "B4 IS ALSO FALSE ON ITS OWN TERMS: at delta = 0 the radial boost "
        "K_1 has a "
        "SMOOTH orbit (only d00, d01, d11 move and the 0i row is d01 n-hat) "
        "with a0 = 33 sym(e0 n-hat) nonzero; measured kin = -6.094e6 at "
        "h = 1.5 and -6.255e6 at h = 0.75, h-stable to 2.6 percent. The "
        "exclusion is real only for PERIODIC clocks.",
        "B5 REFUTED BELOW delta = 0.3: the equal-volume x-axis tube, which "
        "contains no string, captures MORE of Q_I1sq's C4 than the z-axis "
        "tube at delta = 0.003 (0.5868 vs 0.5796, excess -0.0072). At small "
        "delta the omega^4 weight is at the ORIGIN point core; the spherical "
        "core control captures 0.7906.",
        "A4 IS VACUOUS AS EVIDENCE: C4(Q_C6a) equals sum (tr(a0 eta a0 eta))^2 "
        "h^3 exactly (1.3e-15), containing no spatial derivative, so its L^3 "
        "and h^0 exponents are automatic and it cannot report on the frame. "
        "The axis cut and the equal-volume spherical cut give identical "
        "exponents to three decimals.",
        "A1's IR LEG IS NOT CONVERGED: the C4 L-exponent after excision is "
        "not rho_cut independent, running 0.115, 0.123, 0.190, 0.249, 0.400 "
        "at rho_cut 1.5, 3, 4.5, 6, 9. The h leg IS clean and string-specific.",
        "B2's delta/2 LAW IS A z != 0 STATEMENT: on the equatorial ring the "
        "law is (1-delta)/2 instead, which does not vanish at delta = 0 (it "
        "is 0.5) and does vanish at delta = 1; that ring measures the origin "
        "point core rather than the axis.",
        "B3's delta^2 LAW IS ASYMPTOTIC, NOT EXACT: kin/delta^2 drifts by a "
        "factor 2.19 across the ladder and the local exponent runs 1.9998 at "
        "delta 1e-4 down to 1.946 at delta 0.1 to 0.3.",
        "EVERY TERM IS EXACTLY QUARTIC IN omega: the degree-4 fit residual is "
        "at most 2.8e-14 relative over every box, term and mask, so no "
        "degree-6 contamination enters any coefficient.",
        "ATTACK (a) CLOSED EXHAUSTIVELY, NOT BY SAMPLING: solving for every "
        "X in the 6-dimensional so(1,3) whose flow keeps the field frame-free "
        "gives a null space of dimension 2 at delta = 0, spanned by the "
        "perpendicular-pair rotation (tangent identically 0) and the radial "
        "boost (tangent 33, finite-orbit form violation 1.7e-13, not "
        "periodic). No smooth periodic clock exists anywhere in the algebra.",
        "TWO DIVERGENCES, ONLY ONE OF THEM THE STRING: E_static diverges as "
        "h^-1.03 at delta = 0, where the field is provably string-free, and "
        "h^-1.10 at delta = 0.3, so the omega^0 divergence is the hedgehog "
        "POINT core at the origin (|grad M| ~ 1/r gives a quartic density "
        "~ r^-4 and a radial integral ~ 1/h), not the line. The omega^4 "
        "divergence is the opposite: it vanishes at delta = 0 and is cured by "
        "a z-axis tube but not by an equal-volume spherical cut through the "
        "origin. Arm a addresses the LINE only; the point core is untouched "
        "by the excision and stays unresolved.",
        "TRAP IN THE CERTIFIED STACK: gen_catalog normalizes a0 by "
        "max(norm, 1e-300), so at delta = 0 its clock_local generator, whose "
        "RAW tangent norm is 6.5e-15, becomes a unit-norm noise field and "
        "reports a phantom kin of 2.25. Any delta -> 0 study routed through "
        "gen_catalog will see a clock that is not there. The raw norm "
        "confirms B4's rotational leg with the stack's own clock generator."]
    if relax and relax.get("rows"):
        out["new_findings"].insert(1, (
            "THE DEFECT IS PHYSICAL AND RESOLVABLE, NOT A BROKEN ANSATZ "
            "(attack c, the rung's decider). FIRE on the certified stack "
            "melts the biaxiality on the line rather than escaping into the "
            "third dimension: the on-axis gap drops 0.300 to 0.166 (h = 1.5) "
            "and to 0.077 (h = 0.75) while the far field holds at 0.297, V4 "
            "rises from 1.7e-23 to 0.26, and the leading eigenvector stays "
            "aligned with n-hat to 1.4e-3. The core radius is a PHYSICAL "
            "length, 3.98 at h = 1.5 and 3.58 at h = 0.75, not a lattice "
            "scale. The smoothness test is decisive: the near-axis ring "
            "spread does NOT shrink with h before relaxation (ratio 1.031) "
            "and DOES after it (0.702). Controls: the melt is absent about "
            "the non-singular x axis in the same field (gap 0.294 vs 0.166) "
            "and the endpoint is the same at dt0 = 0.01 and 0.02. Crucially "
            "the clock survives and becomes h-convergent: kin = 351.17 and "
            "351.14 at the two spacings, against 426.51 and 445.22 for the "
            "rigid ansatz."))
    out["caveats"] = [
        "SCOPE LIMIT ON ATTACK (c): the relaxations were run at g = 8.0, the "
        "certified stack's default LC_G and the value m5_21_8_b_lattice.relax "
        "hardcodes, NOT at the R9 working vacuum g = 32. This is legitimate "
        "for the defect question because E_u and kin are EXACTLY "
        "g-independent on this family (M_00 = -sg is a spatial constant and "
        "M_0i = 0, so the jets carry a zero 0-row and g cancels): E_u = "
        "62.8517443315 and kin = 426.5070121484 at both g = 8 and g = 32, to "
        "10 digits. But V4, the penalty the core must pay to melt, is NOT "
        "g-independent: at a displacement of 1e-3 off the pinned manifold V4 "
        "is 1.044e-3 at g = 8 and 4.157 at g = 32, a factor 3981 near the "
        "(g ratio)^6 = 4096 expected from the trace targets, and |grad|max "
        "goes 7.37 to 576.4. The melt trades u-energy against V4, so at "
        "g = 32 the equilibrium core will be TIGHTER and costlier than the "
        "radius 3.98 measured here. The existence of a resolved finite core "
        "and the survival of the clock are established at g = 8; their "
        "quantitative values at g = 32 are NOT measured and the next rung "
        "should measure them.",
        "Both relaxations stopped on max_iter, not on a force tolerance, so "
        "the core radii 3.98 and 3.58 are a snapshot of a still-descending "
        "flow and the 10 pct gap between them is partly the different "
        "iteration budgets (3000 at h = 1.5, 2000 at h = 0.75). The "
        "DIRECTION is unambiguous and the smoothness and clock-convergence "
        "tests do not depend on full convergence.",
        "The h-exponents in arm a are two-point slopes from the box pairs "
        "(32, 64) at L = 48 and (48, 96) at L = 72; both pairs are reported "
        "and they agree to about 0.03, so no single pair is carrying a "
        "verdict.",
        "The L-exponents come from three boxes only, which is why the "
        "increment test rather than the fitted exponent is used to judge IR "
        "convergence in A1."]
    # wall time of the whole audit, from the checkpoint file timestamps
    ts = [os.path.getmtime(os.path.join(CKPT, f))
          for f in os.listdir(CKPT) if f.endswith(".json")]
    out["runtime_s"] = round((max(ts) - min(ts)) if len(ts) > 1 else 0.0, 1)
    out["merge_s"] = round(time.time() - t0, 2)
    with open(os.path.join(DATA, "m5_32_r9_audit.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(json.dumps({k: {kk: out[k][kk] for kk in
                          ("verdict", "own_number", "producer_number",
                           "rel_dev")}
                      for k in ("A1", "A2", "A3", "A4",
                                "B1", "B2", "B3", "B4", "B5")}, indent=1,
                     default=float))
    print("confirmed", out["n_confirmed"], "qualified", out["n_qualified"],
          "refuted", out["n_refuted"])
    return out


if __name__ == "__main__":
    main()
