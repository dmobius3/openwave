"""M5.32 R7 audit: assemble ../data/m5_32_r7_audit.json from the checkpoints
written by m5_32_r7_audit.py. Verdicts are set here, each against the ONE
number that decides it. Run: python3 m5_32_r7_audit_merge.py
"""
from __future__ import annotations
import json, os, time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
CKPT = os.path.join(HERE, "..", "checkpoints", "m5_32_r7")


def ck(tag):
    p = os.path.join(CKPT, f"{tag}.json")
    return json.load(open(p)) if os.path.exists(p) else None


def V(verdict, own, prod, note, rel_dev=None):
    if rel_dev is None and isinstance(own, (int, float)) and isinstance(prod, (int, float)):
        d = max(abs(own), abs(prod), 1e-300)
        rel_dev = abs(own - prod) / d
    return {"verdict": verdict, "own_number": own, "producer_number": prod,
            "rel_dev": rel_dev, "note": note}


def build():
    t0 = time.time()
    a1 = ck("a1_covariance"); a2 = ck("a2_static"); an = ck("a2_anchors")
    a3 = ck("a3_clock_omega2"); a45 = ck("a4_a5_channels"); a6 = ck("a6_derrick")
    b48 = ck("b_box_n32_L48"); b72 = ck("b_box_n48_L72")
    b7 = ck("b7_boxes"); b7s = ck("b7_supplement"); c2l = ck("c2_ladder")
    R = {}

    # ------------------------------- A1 -------------------------------
    R["A1"] = V("CONFIRMED", a1["worst_drift"], 4.648239198203908e-09,
                ("Exactly covariant. Worst relative drift 1.5e-6 at rapidity 3; "
                 "at rapidity 1 I measure 3.9e-14, five decades below the "
                 "producer's 4.6e-9, so their figure is their own floor. The "
                 "drift is scale free (2.0e-9 to 3.7e-9 across four decades of "
                 "jet scale), a cancellation floor growing as cosh^4(rapidity), "
                 "not a covariance failure. All four mutants drift by order 1."))
    R["A1"]["mutant_drifts"] = a1["mutants"]
    R["A1"]["own_drift_at_rapidity_1"] = a1["transforms"]["mixed_rap1.0"]
    R["A1"]["rel_drift_vs_jet_scale"] = [r["rel_drift"] for r in
                                         a1["scale_growth"]["rows"]]

    # ------------------------------- A2 -------------------------------
    R["A2"] = V("QUALIFIED", an["max_E_KT_dressed"], 0.0,
                ("The identity is exact: 3.0e-17 worst relative K_T on 20 random "
                 "block-diagonal fields. The INFERENCE fails: 30 of the 44 stored "
                 "M5.32 R3 pair anchors carry a time row (up to 2.88) and E_KT up "
                 "to 3.9e4, and K_T shifts the pair interaction by -1.6e3 to "
                 "-4.7e3 per unit c2, growing with separation. The Coulomb "
                 "anchors are not c2 independent."),
                rel_dev=None)
    R["A2"]["own_number_identity_worst_rel"] = a2["random_block_diagonal"]["worst_rel"]
    R["A2"]["anchors_with_time_row"] = f"{an['n_with_time_row']}/{an['n_total']}"
    R["A2"]["K_T_pair_interaction_per_c2"] = an["K_T_pair_interaction_per_c2"]
    R["A2"]["undressed_anchor_max_abs_E_KT"] = an["max_E_KT_undressed"]
    R["A2"]["mutation_u_to_e1_E_KT"] = a2["mutation_u_to_e1"]["E_KT_mutant"]
    R["A2"]["degeneracy"] = {
        "t_star": a2["degeneracy_path"]["t_star"],
        "K_T_at_t": {str(r["t_over_tstar"]): r.get("K_T", "UNDEFINED")
                     for r in a2["degeneracy_path"]["rows"]},
        "note": ("K_T grows from 2.30 to 5382 as t -> t*, and past t* the "
                 "spectrum of M eta is COMPLEX so u, h, I1_h and K_T are all "
                 "undefined, not merely ill conditioned. Anchor min eigen gap "
                 "is 32, so no anchor is near it.")}

    # ------------------------------- A3 -------------------------------
    R["A3"] = V("CONFIRMED", a3["worst_rel_kappa_to_E_KT"], 6e-15,
                ("Nine omega samples, degree 4 fit: cubic and quartic "
                 "coefficients <= 4.1e-13, so no quartic contamination hides in "
                 "the producer's 3 point fit. The zero is structural, not a "
                 "cancellation: tr(h a0 h a0) and tr(eta a0 eta a0) are both "
                 "0.18 and the u frame a0 time row is 5e-18. Mutant u -> e1 "
                 "gives kappa 1.69e4."))
    R["A3"]["worst_abs_quartic_coeff"] = a3["worst_abs_quartic"]
    R["A3"]["kappa_mutant"] = a3["rows"][0]["kappa_mutant_u_to_e1"]

    # ------------------------------- A4 -------------------------------
    R["A4"] = V("CONFIRMED", a45["channels"]["boost_1"]["K0_per_c2"], 2178.0,
                ("K0 = 2(g + M_jj)^2 exactly (relative deviation 0.0 on boost 1, "
                 "2, 3) and jet independent to 4.5e-13 over three decades of jet "
                 "scale. My own antisymmetric probe flips the sign to -1922. The "
                 "exclusion is legitimate: d_0 M of a symmetric field is "
                 "symmetric, so 12 of the producer's 31 channels are not tangent "
                 "directions at all."))
    R["A4"]["antisym_probe_K0"] = a45["channels"]["antisym_probe"]["K0_per_c2"]
    R["A4"]["jet_independence_worst_abs"] = max(
        c["jet_independence_worst_abs"] for c in a45["channels"].values())

    # ------------------------------- A5 -------------------------------
    lstars = {nm: v["lambda_star_by_c2"] for nm, v in a45["lambda_star"].items()}
    R["A5"] = V("CONFIRMED", 0.5, 0.49999999999999645,
                ("lambda* = 1/2 on every channel for c2 from 0 to 1000. The "
                 "reason is structural: min eig of the omega^2 jet form is "
                 "EXACTLY affine in lambda with its zero at 1/2, while c2 K0 is "
                 "a jet independent constant, so any lambda < 1/2 fails at jet "
                 "scale sqrt(c2 K0 / |min eig|) = 0.11 for c2 = 0.1 on boost 1. "
                 "c2 can never replace lambda."))
    R["A5"]["lambda_star_by_channel_and_c2"] = lstars
    R["A5"]["missing_channel"] = a45["missing_channel"]

    # ------------------------------- A6 -------------------------------
    R["A6"] = V("QUALIFIED", a6["resolution"]["dilation_exponent_h0.75_R6"], 1.087300070117899,
                ("The Derrick exponent is +1 in the continuum, and the producer's "
                 "1.087 is a LATTICE artifact: halving h from 1.5 to 0.75 pulls "
                 "the same R = 6 dilation exponent to 1.021, and at R = 12 it is "
                 "1.034. But the R ladder has no exponent: its pair slopes run "
                 "monotonically 2.01 -> 0.73 over R = 3 to 48, so the quoted "
                 "1.13 to 1.31 is one window of a monotone trend, and above "
                 "R ~ 36 the localizer penalty falls below linear."))
    R["A6"]["dilation_exponent_h1.5_R6"] = a6["dilation_R6"]["loglog_exponent"]
    R["A6"]["dilation_exponent_h1.5_R12"] = a6["dilation_R12"]["loglog_exponent"]
    R["A6"]["R_ladder_pair_exponents"] = a6["R_ladder"]["pair_exponents"]
    R["A6"]["amp_exponent"] = a6["amp_exponent_R12"]
    R["A6"]["all_E_KT_positive"] = a6["positivity"]["all_positive"]
    return R, (b48, b72, b7, b7s, c2l), t0


def build_B(R, b48, b72, b7, b7s, c2l):
    r48, r72 = b48["rows"], b72["rows"]
    lin, dr = b7["linearity"], b7["drift"]

    def row(k):
        return r48[k], r72[k]

    # ------------------------------- B1 -------------------------------
    a, b = row("lam_0.75_c2_0")
    R["B1"] = V("CONFIRMED", a["omega_star"], 0.2001285645398461,
                ("With c2 = 0 the fixed J minimizer sits at the wall in both "
                 "boxes on MY 16 point R grid too (R* = 24 and 36, argmin index "
                 "15 of 15). omega* reproduces the producer to 1.1e-3 relative "
                 "with an independent stencil order, amp grid and minimizer."))
    R["B1"]["R_star_L48"] = a["R_star_grid"]; R["B1"]["R_star_L72"] = b["R_star_grid"]
    R["B1"]["amp_star_L48"] = a["amp_star"]

    # ------------------------------- B2 -------------------------------
    a3, b3 = row("lam_0.75_c2_0.03")
    a1_, b1_ = row("lam_0.75_c2_0.1")
    R["B2"] = V("QUALIFIED", a3["R_star_grid"], 14.412811041888046,
                ("An interior R* does appear at lambda = 0.75 for c2 in {0.03, "
                 "0.1} in BOTH boxes, reproduced independently. But it is the "
                 "dressing switching OFF, not a clock localizing: amp* falls "
                 "monotonically 0.02510 -> 0.02330 -> 0.02064 -> 0.01407 -> 0 "
                 "across c2 = 0, 0.01, 0.03, 0.1, 0.3, and the dressing is worth "
                 "-23.3%, -17.4%, -8.6%, -0.96%, 0.00% of E_J. The whole "
                 "interior window is the interval over which the dressing dies."))
    R["B2"]["amp_star_vs_c2_L48"] = {
        k.split("c2_")[1]: r48[k]["amp_star"] for k in
        ["lam_0.75_c2_0", "lam_0.75_c2_0.01", "lam_0.75_c2_0.03",
         "lam_0.75_c2_0.1", "lam_0.75_c2_0.3"]}
    R["B2"]["dressing_gain_frac_vs_c2_L48"] = {
        k.split("c2_")[1]: r48[k]["dressing_gain_frac"] for k in
        ["lam_0.75_c2_0", "lam_0.75_c2_0.01", "lam_0.75_c2_0.03",
         "lam_0.75_c2_0.1", "lam_0.75_c2_0.3"]}
    R["B2"]["R_star_box_drift_at_c2_0.03"] = b3["R_star_grid"] / a3["R_star_grid"] - 1.0
    R["B2"]["E_J_of_R_flatness_c2_0.03_L48"] = {
        "min_minus_wall": min(a3["E_J_of_R"]) - a3["E_J_of_R"][-1],
        "E_J": min(a3["E_J_of_R"]),
        "depth_below_wall_frac": (min(a3["E_J_of_R"]) - a3["E_J_of_R"][-1])
                                 / min(a3["E_J_of_R"]),
        "note": ("the L48 c2 = 0.03 curve is flat: E_J(R) - E_J(undressed) is "
                 "-7.20, -7.37, -7.38, -7.34, -7.26, -7.19, -7.08, -7.05 over "
                 "R = 10.5 to 24, so the interior minimum is 0.43% below the "
                 "wall and R* is only fixed to within a factor ~2. My finer grid "
                 "puts it at 13.5 where the producer's parabola gave 14.41")}
    R["B2"]["R_star_L72_at_c2_0.03"] = b3["R_star_grid"]
    R["B2"]["c2_0.1_R_star_both_boxes"] = [a1_["R_star_grid"], b1_["R_star_grid"]]

    # ------------------------------- B3 -------------------------------
    R["B3"] = V("CONFIRMED", r48["lam_0.75_c2_0.3"]["amp_star"], 0.0,
                ("amp* = 0 exactly at c2 >= 0.3 for lambda = 0.75 and at "
                 "c2 >= 0.1 for lambda = 1, in both boxes, with E_J falling back "
                 "onto the undressed value 86.29801769900872 to machine "
                 "precision. Localization is a window, not a trend."))
    R["B3"]["lam1_c2_0.1_amp_star"] = r48["lam_1_c2_0.1"]["amp_star"]
    R["B3"]["E_J_matches_undressed"] = abs(
        r48["lam_0.75_c2_0.3"]["E_J"] - r48["lam_0.75_c2_0.3"]["E_J_undressed"])

    # ------------------------------- B4 -------------------------------
    g7 = {"c2_interior_both_boxes_lam0.75": [], "factor": None}
    for c2 in ("0", "0.01", "0.03", "0.1", "0.3"):
        k = f"lam_0.75_c2_{c2}"
        if r48[k]["interior"] and r72[k]["interior"]:
            g7["c2_interior_both_boxes_lam0.75"].append(float(c2))
    if len(g7["c2_interior_both_boxes_lam0.75"]) > 1:
        g7["factor"] = (max(g7["c2_interior_both_boxes_lam0.75"])
                        / min(g7["c2_interior_both_boxes_lam0.75"]))
    if c2l:
        for tag, box in (("L48", "n32"), ("L72", "n48")):
            ints = sorted(float(k.split("c2_")[1]) for k, v
                          in c2l[box]["rows"].items()
                          if k.startswith("lam_0.75_") and v["interior"])
            g7[f"dense_ladder_interior_{tag}"] = [min(ints), max(ints)] if ints else []
        i48 = set(k for k, v in c2l["n32"]["rows"].items()
                  if k.startswith("lam_0.75_") and v["interior"])
        i72 = set(k for k, v in c2l["n48"]["rows"].items()
                  if k.startswith("lam_0.75_") and v["interior"])
        both = sorted(float(k.split("c2_")[1]) for k in (i48 & i72))
        g7["dense_ladder_interior_both"] = [min(both), max(both)] if both else []
        g7["dense_ladder_factor"] = (max(both) / min(both)) if both else None
    R["B4"] = V("QUALIFIED", g7.get("dense_ladder_factor", g7["factor"]),
                3.3333333333333335,
                ("G7_met = false is CONFIRMED but HALF THE REASON IS WRONG. The "
                 "producer's factor 3.33 is an artifact of a c2 ladder whose own "
                 "spacing is 3.33; on a 33 point ladder the interior-in-both "
                 "range at lambda = 0.75 is c2 in [0.0281, 0.1369], a factor "
                 "4.87, so the RANGE half of G7 is MET, not failed. G7 fails "
                 "solely on the drift gate: the smallest |omega* box drift| over "
                 "every c2 and lambda is 0.301, three times the 0.10 gate."))
    R["B4"]["range_detail"] = g7
    R["B4"]["min_abs_omega_drift"] = b7s["omega_drift_closest_to_zero"]
    R["B4"]["G7_met"] = False

    # ------------------------------- B5 -------------------------------
    R["B5"] = V("CONFIRMED", r48["lam_0.75_c2_-0.1"]["E_J"], -26.07022297523679,
                ("c2 = -0.1 reddens exactly as claimed: R* pinned at the wall in "
                 "both boxes, amp* grows 0.02510 -> 0.03858 (L48) and 0.02467 -> "
                 "0.04302 (L72), and E_J goes to -26.27 (L48) and -107.63 (L72), "
                 "the negative energy DEEPENING with box size, which is the "
                 "unbounded signature. My value matches the producer to 0.7%."))
    R["B5"]["E_J_L72"] = r72["lam_0.75_c2_-0.1"]["E_J"]
    R["B5"]["amp_star_L48"] = r48["lam_0.75_c2_-0.1"]["amp_star"]

    # ------------------------------- B6 -------------------------------
    R["B6"] = V("CONFIRMED", b48["undressed"]["lam_0.75"]["E_KT"],
                1.6253630521568044e-15,
                ("E_KT on the undressed hedgehog is 1.59e-15 (L48) and 1.84e-15 "
                 "(L72) against an internal trace scale of order 1.2e3, so "
                 "E_stat is c2 independent to 1e-18 relative; the vacuum has "
                 "zero jets so K_T is identically 0 there. Confirmed by the "
                 "block-diagonal identity of A2, which makes both exact rather "
                 "than numerical."))
    R["B6"]["E_KT_undressed_L72"] = b72["undressed"]["lam_0.75"]["E_KT"]

    # ------------------------------- B7 -------------------------------
    R["B7"] = V("QUALIFIED", lin["affine_fit"][0], 9.658318180039402,
                ("Every link of the mechanism is CONFIRMED and one word is not. "
                 "MEASURED: |a0| is exactly constant in r (radial exponent "
                 "4.8e-15, value 0.3*sqrt(2) everywhere), the hedgehog gradient "
                 "falls as r^-0.983, the kinetic density as r^-1.973, and the "
                 "shell integral is a flat plateau 16.25 +- 0.45 per unit r from "
                 "r = 5 to r = L/2. Only 6.0% of kin comes from r < 6 and 14.5% "
                 "from r < 12, so it is TAIL dominated, not core dominated. Four "
                 "boxes give kin = 9.66997 L - 37.78 with 3.1e-4 worst relative "
                 "residual, 11x better than the best power law. NOT an artifact: "
                 "fwd, bwd and sym stencils give identical kin to 1e-13, and the "
                 "kinetic density in the L = 48 box is BIT IDENTICAL to the "
                 "L = 120 box on every shared interior cell. REFUTED WORD: "
                 "'INDEPENDENT of c2 and of lambda' - the drift ranges over "
                 "-0.301 to -0.352, a 14.4% spread. The conclusion is untouched: "
                 "every value is 3x the gate."))
    R["B7"]["four_box_kin"] = {"L": lin["L"], "kin": lin["kin"]}
    R["B7"]["affine_fit_slope_intercept"] = lin["affine_fit"]
    R["B7"]["affine_worst_rel_residual"] = max(lin["affine_residual_rel"])
    R["B7"]["powerlaw_worst_rel_residual"] = max(lin["powerlaw_residual_rel"])
    R["B7"]["loglog_exponent_pairs"] = lin["loglog_exponent_pairs"]
    R["B7"]["producer_exponent_1.07_is_a_two_box_artifact"] = {
        "pairs": lin["loglog_exponent_pairs"],
        "note": ("1.0705, 1.0499, 1.0389 for L = 48/72, 72/96, 96/120: the "
                 "effective exponent decreases monotonically to 1 because the "
                 "affine intercept is negative (-37.78). The asymptotic exponent "
                 "is exactly 1, not 1.07.")}
    R["B7"]["omega_star_times_L"] = dr["omega_star_times_L"]
    R["B7"]["drift_vs_pure_1_over_L"] = {
        "measured": dr["drift_pairs"], "pure_1_over_L": dr["pure_1_over_L_prediction"]}
    R["B7"]["radial_exponents"] = {
        "a0_norm": b7s["a0_norm_radial_exponent"],
        "hedgehog_gradient": b7s["grad_exponent"],
        "kinetic_density": b7s["kin_density_exponent"],
        "fit_window_r": b7s["fit_window_r"]}
    R["B7"]["core_vs_tail"] = b7["shells"]["fraction_inside"]
    R["B7"]["plateau_kin_per_unit_r"] = [
        b7["shells"]["plateau_kin_per_unit_r_mean"],
        b7["shells"]["plateau_kin_per_unit_r_std"]]
    R["B7"]["boundary_artifact_test"] = {
        "max_abs_density_diff_interior_cells_L48_vs_L120":
            b7["cell_overlap"]["max_abs_diff_interior_cells"],
        "rel_dev_sum_interior": b7["cell_overlap"]["rel_dev_sum_interior"],
        "boundary_shell_only_diff": (b7["cell_overlap"]["boundary_shell_sum_big"]
                                     - b7["cell_overlap"]["boundary_shell_sum_small"])}
    R["B7"]["stencil_independence"] = {
        k: v["slope_per_unit_L"] for k, v in b7["stencils"].items()}
    R["B7"]["delta_ladder"] = {
        "slope_per_unit_L": {str(r["delta"]): r["slope_per_unit_L"]
                             for r in b7["delta_ladder"]},
        "scaling_exponent_in_delta": b7["delta_scaling_exponent"],
        "note": ("NEW: at delta = 0 the clock generator a0 is IDENTICALLY zero "
                 "and kin = 0 exactly; the extensive slope scales as "
                 "delta^1.9987. The whole IR-divergent clock inertia is carried "
                 "by the vacuum's delta entry.")}
    R["B7"]["omega_drift_by_channel"] = b7s["omega_drift_by_channel"]
    R["B7"]["omega_drift_spread_over_undressed"] = b7s["omega_drift_spread_over_undressed"]
    return R


NEW_FINDINGS = [
 {"id": "N1", "claim_touched": "B7",
  "finding": ("The kinetic inertia is EXACTLY affine in L, not a power law with "
              "exponent 1.07. Four boxes (L = 48, 72, 96, 120 at fixed h = 1.5) "
              "give kin = 9.66997 L - 37.782 with worst relative residual "
              "3.1e-4, while the best power law leaves 3.4e-3, an 11x worse fit. "
              "The effective log-log exponent falls monotonically 1.0705, "
              "1.0499, 1.0389 toward 1. The producer's 1.07 is a two-box "
              "artifact of the negative core intercept."),
  "numbers": {"slope": 9.669972189944541, "intercept": -37.78221100095749,
              "affine_residual": 3.061e-4, "powerlaw_residual": 3.363e-3}},
 {"id": "N2", "claim_touched": "B7",
  "finding": ("The extensive growth is NOT a boundary or stencil artifact, by a "
              "test the producer did not run: the kinetic density in the L = 48 "
              "box is BIT IDENTICAL (max absolute difference 0.0) to the same "
              "cells of the L = 120 box on every interior cell, the two grids "
              "sharing h = 1.5 and alignment. Only the one-cell outer shell "
              "differs (26.80 vs 28.93). Forward, backward and symmetric "
              "stencils give the same kin to 1e-13 and the same slope 9.6653."),
  "numbers": {"max_interior_cell_diff": 0.0, "slope_fwd": 9.665306445847623,
              "slope_bwd": 9.665306445847602, "slope_sym": 9.665306445847614}},
 {"id": "N3", "claim_touched": "B7",
  "finding": ("The IR divergence is carried entirely by the vacuum parameter "
              "delta. At delta = 0 the clock generator a0 = [G1, diag(g,1,delta,0)] "
              "is IDENTICALLY zero and kin = 0 exactly, and the extensive slope "
              "scales as delta^1.9987 over delta = 0.05 to 0.5. Setting delta "
              "toward 0 is therefore the only lever in the current vacuum that "
              "touches the divergence, and it destroys the clock with it."),
  "numbers": {"delta_exponent": 1.9986923769961822,
              "slope_at_delta_0.3": 9.665306445847614, "kin_at_delta_0": 0.0}},
 {"id": "N4", "claim_touched": "B7",
  "finding": ("Direct radial fits of the three links in the mechanism, none "
              "reported by the producer: |a0| is EXACTLY r independent (fitted "
              "radial exponent 4.8e-15, value 0.42426406588 = delta*sqrt(2) at "
              "every shell), the hedgehog gradient falls as r^-0.983 and the "
              "kinetic density as r^-1.973. The shell integral is a flat "
              "plateau of 16.25 +- 0.45 per unit r from r = 5 to r = L/2, and "
              "only 6.0% of kin sits inside r = 6, so the term is tail "
              "dominated. The core-dominated alternative reading is refuted."),
  "numbers": {"a0_exponent": 4.76e-15, "grad_exponent": -0.9832440392382437,
              "density_exponent": -1.9729389166917477,
              "fraction_r_lt_6": 0.060215005234976665,
              "fraction_r_lt_12": 0.14498640836111107}},
 {"id": "N5", "claim_touched": "B7",
  "finding": ("omega* x L converges: 11.254, 10.937, 10.781, 10.688 over the "
              "four boxes, with the increments halving, so omega* -> C/L with C "
              "about 10.6 at J = 200. The measured box drift tracks the pure "
              "1/L law to within 2 points in 100 (-0.352 vs -0.333, -0.261 vs "
              "-0.250, -0.207 vs -0.200) and converges onto it."),
  "numbers": {"omega_times_L": [11.2542, 10.9372, 10.7812, 10.6883]}},
 {"id": "N6", "claim_touched": "A2",
  "finding": ("30 of the 44 stored M5.32 R3 pair anchors in the campaign's own "
              "data folder carry a nonzero time row (up to 2.88, 8.9% of the "
              "field scale) and E_KT up to 3.94e4. K_T shifts the pair "
              "interaction energy by -1.6e3 (d = 10) to -4.7e3 (d = 24) per unit "
              "c2, monotonically in separation, against an R3 dressing-part Eint "
              "of order 3e3 to 5e3. So the Coulomb-force anchors are NOT c2 "
              "independent, even though the block-diagonal identity is exact."),
  "numbers": {"n_with_time_row": 30, "n_total": 44, "max_E_KT": 39384.876,
              "dEint_per_c2_d10": -1623.21, "dEint_per_c2_d24": -4729.65}},
 {"id": "N7", "claim_touched": "A2",
  "finding": ("Beyond the h_cov degeneracy locus t* = (g+1)/2 = 16.5 the "
              "spectrum of M eta goes COMPLEX, so u, h_cov, I1_h and K_T are all "
              "UNDEFINED, not merely ill conditioned; approaching it K_T grows "
              "2.30 -> 5382 as the eigen gap falls 32 -> 0.47. This is a defect "
              "of the whole h_cov family (the R2 lambda candidate included), not "
              "of K_T alone. No stored anchor is near it (minimum gap 32)."),
  "numbers": {"t_star": 16.5, "K_T_at_0.9999_tstar": 5381.6,
              "anchor_min_eigen_gap": 31.9972}},
 {"id": "N8", "claim_touched": "A5",
  "finding": ("The producer's 31 channel list contains no channel built from the "
              "clock flow the model actually uses. B8.a0_unit on the m5_21_8 "
              "hedgehog is symmetric AND block diagonal, so its K0 is EXACTLY 0: "
              "c2 gives literally zero omega^2 help on the only clock channel "
              "the lattice runs realize. 'c2 never hurts' is true there only "
              "because it also never acts."),
  "numbers": {"K0_hedgehog_clock": 0.0, "K0_boost_1": 2178.0}},
 {"id": "N9", "claim_touched": "A6",
  "finding": ("The Derrick exponent 1.087 is a lattice artifact, not the "
              "continuum value. Halving h from 1.5 to 0.75 pulls the same R = 6 "
              "dilation exponent from 1.0856 to 1.0206, and at R = 12 (a "
              "better-resolved dressing) it is 1.0336. The exact continuum "
              "exponent is 1, which the audit reaches from below as h shrinks. "
              "Separately the R ladder has no exponent at all: its pair slopes "
              "run monotonically 2.01, 1.68, 1.45, 1.29, 1.19, 1.12, 1.01, 0.73 "
              "over R = 3 to 48, so the localizer penalty is superlinear at "
              "small R and SUBLINEAR above R about 36."),
  "numbers": {"exp_h1.5_R6": 1.0856250850160385, "exp_h0.75_R6": 1.0206046354974883,
              "exp_h1.5_R12": 1.0336409429834257, "R_ladder_last_pair": 0.7321625648331863}},
 {"id": "N10", "claim_touched": "B2",
  "finding": ("The c2 = 0.03 interior minimum is not a resolved feature. On my "
              "16 point R grid the L48 curve E_J(R) - E_J(undressed) reads "
              "-7.20, -7.37, -7.38, -7.34, -7.26, -7.19, -7.08, -7.05 over "
              "R = 10.5 to 24: the interior minimum is 0.34 below the wall, "
              "0.43% of E_J, over an R interval spanning a factor 2.3. My grid "
              "puts R* at 13.5 where the producer's parabola gave 14.41, and "
              "the L72 box puts it at 12.0, an 11% box drift in R* itself."),
  "numbers": {"depth_below_wall": -0.3363, "depth_frac_of_E_J": -0.00426,
              "R_star_L48_mine": 13.5, "R_star_L48_producer": 14.4128,
              "R_star_L72_mine": 12.0}},
 {"id": "N11", "claim_touched": "B2",
  "finding": ("The interior window coincides exactly with the death of the "
              "dressing. Across c2 = 0, 0.01, 0.03, 0.1, 0.3 at lambda = 0.75, "
              "amp* falls monotonically 0.02510, 0.02330, 0.02064, 0.01407, 0 "
              "and the dressing is worth -23.3%, -17.4%, -8.6%, -0.96%, 0.00% "
              "of E_J. At c2 = 0.1, the deepest 'interior' point, the localized "
              "clock buys 0.96% of the total energy at 56% of the undressed "
              "amplitude. A localized minimum whose amplitude is collapsing is "
              "the undressed solution in disguise."),
  "numbers": {"amp_star_c2_0": 0.02510, "amp_star_c2_0.1": 0.01407,
              "gain_frac_c2_0": -0.2331, "gain_frac_c2_0.1": -0.0096}},
 {"id": "N12", "claim_touched": "A1",
  "finding": ("K_T's covariance floor is not 4.6e-9: at rapidity 1 I measure "
              "3.9e-14 with my own generator basis and transformation, five "
              "decades better than the producer's figure, so 4.6e-9 measures "
              "their implementation and not the term. The true floor grows as "
              "cosh^4(rapidity) from the tr(hAhA) cancellation, reaching 1.5e-6 "
              "at rapidity 3, and the RELATIVE drift is flat (2.0e-9 to 3.7e-9) "
              "across four decades of jet scale, which is the signature of a "
              "cancellation floor rather than a covariance failure."),
  "numbers": {"own_rap1": 3.93e-14, "producer_rap1": 4.648e-9, "own_rap3": 1.463e-6}},
 {"id": "N14", "claim_touched": "B4",
  "finding": ("The G7 coefficient-range failure is a c2-grid artifact. The "
              "producer's ladder {0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3} is "
              "spaced by exactly 3.33, so its measured range factor 3.33 is one "
              "grid step and cannot resolve the >= 4 bar. On a 33 point "
              "logarithmic ladder from 0.005 to 0.5 the interior-in-both-boxes "
              "window at lambda = 0.75 is c2 in [0.0281, 0.1369], a factor 4.87, "
              "so the range criterion is actually MET. G7 still fails, but "
              "solely on the omega* drift gate."),
  "numbers": {"dense_range_low": 0.028117, "dense_range_high": 0.136921,
              "dense_factor": 4.869687377743002, "producer_factor": 3.3333333,
              "n_ladder_points": 33}},
 {"id": "N13", "claim_touched": "A4",
  "finding": ("12 of the producer's 31 LP channels (clock_t, clock_tr, the six "
              "x_*_probe entries and the four aud_clock_random_*) are built from "
              "ANTISYMMETRIC a0 and are therefore not tangent directions of a "
              "symmetric matrix field at all. Excluding them is correct, but it "
              "means 39% of the channel list is inert, and the exclusion is what "
              "makes 'c2 never hurts' true, since every one of them carries "
              "K0 < 0 (my own probe: -1922)."),
  "numbers": {"n_antisymmetric": 12, "n_channels": 31, "own_probe_K0": -1922.0}},
]


def main():
    R, (b48, b72, b7, b7s, c2l), t0 = build()
    R = build_B(R, b48, b72, b7, b7s, c2l)
    n_c = sum(1 for v in R.values() if v["verdict"] == "CONFIRMED")
    n_q = sum(1 for v in R.values() if v["verdict"] == "QUALIFIED")
    n_r = sum(1 for v in R.values() if v["verdict"] == "REFUTED")
    out = {
        "task": "M5.32 R7 INDEPENDENT ADVERSARIAL AUDIT (class C4, the term K_T)",
        "auditor_scripts": ["m5_32_r7_audit.py", "m5_32_r7_audit_merge.py"],
        "producer_scripts_NOT_read": ["m5_32_r7_a_kt_form.py",
                                      "m5_32_r7_b_kt_lattice.py"],
        "oracles_used": ["m5_21_3_a_4d.py (certified stencil, e_parts, V4)",
                         "m5_21_8_b_lattice.py (hedgehog, a0_unit)",
                         "own K_T, own h_cov eigen-solve, own u-frame boost, own "
                         "finite differences, own energy assembly with the branch "
                         "order reversed, own amp grid, own two-stage minimizer"],
        "summary": (
            "The term K_T is exactly what the producer says it is, and the rung's "
            "headline B7 survives every attack I could mount. K_T is Lorentz "
            "covariant to the float cancellation floor, vanishes identically on "
            "block-diagonal fields, has a structurally (not numerically) zero "
            "omega^2 coefficient on the R4 dressing family, carries the constant "
            "K0 = 2(g + M_jj)^2 on boost channels, and cannot move the R2 form "
            "level no-go: lambda* stays 1/2 for c2 up to 1000. Three claims "
            "narrow. A2's identity is exact but its INFERENCE is false: 30 of the "
            "44 stored R3 pair anchors carry a time row and E_KT up to 3.9e4, so "
            "the Coulomb anchors do move with c2. A6's exponent 1.087 is a "
            "lattice artifact (1.021 at half the spacing) and the R ladder has no "
            "single exponent (2.01 down to 0.73). B2's interior R* is real but is "
            "the dressing switching off: over the whole interior window amp* "
            "collapses 0.0251 to 0.0141 and the dressing's worth falls from 23.3% "
            "to 0.96% of E_J, and at c2 = 0.03 the minimum is 0.43% deep on a "
            "plateau spanning a factor 2.3 in R. B7 is the strongest result of "
            "the rung: four boxes give kin = 9.66997 L - 37.78 with 3.1e-4 "
            "residual (11x better than any power law, so the asymptotic exponent "
            "is exactly 1, not 1.07); |a0| is exactly r independent, the gradient "
            "falls as r^-0.983 and the density as r^-1.973; the shell integral is "
            "a flat 16.25 per unit r out to the wall with only 6% inside r = 6; "
            "and the density is BIT IDENTICAL between the L48 and L120 boxes on "
            "every shared interior cell, which kills the boundary-artifact "
            "reading outright. The one word to strike is 'INDEPENDENT': the drift "
            "ranges -0.301 to -0.352 with c2 and lambda, a 14.4% spread, though "
            "every value is 3x the 0.10 gate so the conclusion is untouched. New: "
            "at delta = 0 the clock generator vanishes identically and kin is "
            "exactly 0, with the extensive slope scaling as delta^1.9987, so the "
            "whole IR divergence rides on the vacuum's delta entry. On B4 the "
            "verdict G7_met = false survives but half its stated reason does "
            "not: the producer's coefficient-range factor 3.33 is an artifact "
            "of a c2 ladder spaced by 3.33, and a 33 point ladder measures "
            "4.87, so the range criterion is actually MET and G7 fails only on "
            "the drift gate."),
        "n_confirmed": n_c, "n_qualified": n_q, "n_refuted": n_r,
        "claims": R,
        "new_findings": NEW_FINDINGS,
        "runtime_s": round(time.time() - t0, 1),
    }
    p = os.path.join(DATA, "m5_32_r7_audit.json")
    with open(p, "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(f"wrote {p}: {n_c} confirmed, {n_q} qualified, {n_r} refuted")
    for k, v in R.items():
        print(f"  {k:<4} {v['verdict']:<10} own={v['own_number']}")


if __name__ == "__main__":
    main()
