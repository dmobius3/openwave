"""P1A Qualification Runner — P1A.4a ladder qualification pass.

P1A.0: P0 freeze verification (rerun).
P1A.1-P1A.3: PASS, unchanged, not rerun.
P1A.4a: Cloud admissibility gate at fixed k=110, sector extraction,
        contamination qualification on qualified refinement family,
        holdouts with correct labels, discrimination test.
P1A.5: uninspected, deliberately.

Permanent status record:
  - Original per-sector contamination rulings: INVALID
  - First pooled replacement rulings: INVALID
  - Cause: calibration set never qualified as a common refinement family;
    K_STENCIL_RULE changed the discretization scheme across the ladder;
    extraction failures (||P_spec||=0) entered downstream calculations.
"""

import os
import sys
import time
import numpy as np

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, BASE_DIR)

from p0.group import (build_icosians, build_character_table,
                      LABELS, DIMS, MCKAY_DIST)
from p0.representations import build_all_representations
from p0.cloud import fibonacci_seeds_s3, build_orbit_cloud
from p0.bundle_operator import build_L_bundle
from p0.frozen_tolerances import GRID_PARAMS

from p1a.freeze_p0 import verify_p0_frozen
from p1a.mass_matrix import build_Mh_base, build_Mh_rho
from p1a.diagnostics import self_adjointness_residual
from p1a.subspace import extract_and_score
from p1a.ladder import analytic_riesz_params
from p1a.cloud_qualify import (K_FIXED, RBF_M, RBF_P,
                                CALIBRATION_CANDIDATES, HOLDOUT_FIXED,
                                FINER_HOLDOUT_SEARCH, FINER_HOLDOUT_C,
                                GATE, N_STABILITY_REF,
                                cloud_gate, is_extraction_ok,
                                compute_fill_distance)
from p1a.contamination import (numerical_floor, analytic_scale,
                                pooled_convergence_fit, envelope_value,
                                evaluate_entry)

N_PROD = GRID_PARAMS["n_seeds"]


def _fmt(v, fmt):
    if isinstance(v, (int, float, np.floating, np.integer)):
        return f"{v:{fmt}}"
    return str(v)


def _d_rho_map():
    return {label: MCKAY_DIST[i] for i, label in enumerate(LABELS)}


# ---------------------------------------------------------------------------
#  Sector extraction at one resolution — k=K_FIXED throughout
# ---------------------------------------------------------------------------

def _run_sector_extraction(X, oid, gid, W, n_seeds, elems, reps,
                           print_fn=print):
    h = compute_fill_distance(X)
    entries = {}
    for label in LABELS:
        idx = LABELS.index(label)
        d_fiber = DIMS[idx]
        d_rho = MCKAY_DIST[idx]
        expected_lambda = d_rho * (d_rho + 2)
        expected_dim = d_rho + 1
        try:
            L, _ = build_L_bundle(X, oid, gid, elems, reps[label],
                                   k=K_FIXED)
            Mh_diag = build_Mh_rho(W, d_fiber)
            eps_sa = self_adjointness_residual(L, Mh_diag)
            L_norm = float(np.linalg.norm(L, 2))
            centre, radius = analytic_riesz_params(label)
            ex, Q = extract_and_score(
                L, Mh_diag, expected_lambda, expected_dim, label=label,
                riesz_centre=centre, riesz_radius=radius)
            im_max = 0.0
            if "imaginary" in ex:
                im_max = ex["imaginary"]["max_im"]
            ok = is_extraction_ok(ex)
            entry = {
                "eps_SA": eps_sa,
                "cluster_position": ex.get("cluster_position"),
                "cluster_spread": ex.get("cluster_spread"),
                "im_max": im_max,
                "riesz_theta": ex.get("riesz", {}).get(
                    "theta_max_schur_vs_riesz"),
                "P_spec_norm": ex.get("riesz", {}).get("P_spec_norm"),
                "P_spec_norm_mh": ex.get("riesz", {}).get("P_spec_norm_mh"),
                "invariance_residual_rel": ex.get("projector", {}).get(
                    "invariance_residual_rel"),
                "L_norm": L_norm,
                "fill_distance": h,
                "pass": ex.get("pass", False),
                "extraction_ok": ok,
                "extraction_fail": not ok,
            }
            entries[label] = entry
            tag = "OK" if ok else "EXTRACTION_FAIL"
            print_fn(f"    {label}: eps_SA={eps_sa:.4e}, "
                     f"im={im_max:.2e}, ||L||={L_norm:.1f}, {tag}")
        except Exception as e:
            print_fn(f"    {label}: ERROR - {e}")
            entries[label] = {"error": str(e), "extraction_fail": True}
    return entries, h


# ---------------------------------------------------------------------------
#  Holdout runner (gate + extraction + evaluation)
# ---------------------------------------------------------------------------

def _run_holdout(n_seeds, elems, reps, drho_map, fit, print_fn=print):
    seeds = fibonacci_seeds_s3(n_seeds)
    X, oid, gid = build_orbit_cloud(seeds, elems)
    gp, gd = cloud_gate(X, oid, gid, elems, reps["R0"],
                         k=K_FIXED, m=RBF_M, p=RBF_P,
                         print_fn=lambda s: print_fn(f"  {s}"))
    print_fn(f"  Gate: {'PASS' if gp else 'FAIL'}")

    evals = {}
    h_r = gd["h"]
    if not gp:
        for label in LABELS:
            evals[(n_seeds, label)] = {
                "regime": "GATE_FAIL", "eligibility": "GATE_FAIL"}
        return gd, None, evals

    W = build_Mh_base(X, oid, n_seeds)
    entries, h_r = _run_sector_extraction(X, oid, gid, W, n_seeds,
                                          elems, reps, print_fn)
    for label, entry in entries.items():
        if entry.get("extraction_fail"):
            evals[(n_seeds, label)] = {
                "regime": "EXTRACTION_FAIL",
                "im_max": entry.get("im_max"), "h": h_r,
                "eligibility": "EXTRACTION_FAIL"}
            continue
        im_max = entry.get("im_max")
        L_norm = entry.get("L_norm")
        P_spec_norm = entry.get("P_spec_norm")
        if im_max is None or L_norm is None:
            evals[(n_seeds, label)] = {
                "regime": "ERROR", "eligibility": "NO_LABEL"}
            continue
        d_rho = drho_map[label]
        s_rho = analytic_scale(d_rho)
        ev = evaluate_entry(im_max, L_norm, P_spec_norm, s_rho, h_r, fit)
        evals[(n_seeds, label)] = ev
        print_fn(f"    {label}: im={im_max:.2e}, "
                 f"{ev['regime']}, {ev['eligibility']}")
    return gd, h_r, evals


# ---------------------------------------------------------------------------
#  Discrimination test
# ---------------------------------------------------------------------------

def _run_discrimination_test(cal_evals, fit, drho_map):
    if fit.get("status") != "ok" or fit.get("p", 0) <= 0:
        print("  No valid pooled fit; cannot run.")
        return {"status": "no_fit", "pass": False}
    test_key = None
    for key, ev in cal_evals.items():
        if ev.get("regime") == "ABOVE_FLOOR" and \
                ev.get("eligibility") == "QUALIFIED":
            test_key = key
            break
    if test_key is None:
        print("  No ABOVE_FLOOR QUALIFIED point available.")
        return {"status": "no_test_point", "pass": False}
    n_s, label = test_key
    ev = cal_evals[test_key]
    d_rho = drho_map[label]
    s_rho = analytic_scale(d_rho)
    E = envelope_value(s_rho, fit["A"], fit["p"], fit["r_max"], ev["h"])
    injected_im = 10.0 * E
    injected_ev = evaluate_entry(injected_im, ev["L_norm_2"],
                                  ev["P_spec_norm_used"], s_rho,
                                  ev["h"], fit)
    rejected = injected_ev.get("eligibility") == "NO_LABEL"
    print(f"  Source: {label}@n={n_s}")
    print(f"  Envelope E = {E:.2e}")
    print(f"  Injected im = {injected_im:.2e} (10x E)")
    print(f"  Result: {'REJECTED — PASS' if rejected else 'ACCEPTED — FAIL'}")
    return {
        "status": "ran",
        "source": f"{label}@n={n_s}",
        "original_im": ev["im_max"],
        "envelope": E,
        "injected_im": injected_im,
        "injected_I_over_E": injected_im / E if E > 0 else float('inf'),
        "result_eligibility": injected_ev.get("eligibility"),
        "rejected": rejected,
        "pass": rejected,
    }


# ---------------------------------------------------------------------------
#  Note generator
# ---------------------------------------------------------------------------

def produce_note(results, elapsed):
    L = []

    L.append("# M8.4 P1A Qualification Note — P1A.4a Ladder Qualification")
    L.append("")
    L.append(f"**Produced:** {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    L.append(f"**Python:** {sys.version.split()[0]}")
    L.append(f"**NumPy:** {np.__version__}")
    L.append(f"**SciPy:** {__import__('scipy').__version__}")
    L.append(f"**Total elapsed:** {elapsed:.1f}s")
    L.append("")

    # ── P1A.0 ──
    p0 = results["p1a_0"]
    L.append("## P1A.0: P0 Freeze Verification")
    L.append("")
    L.append(f"**{'PASS' if p0['pass'] else 'FAIL'}**")
    L.append("")
    for fname, info in sorted(p0.get("details", {}).items()):
        if isinstance(info, dict) and "match" in info:
            L.append(f"- `{fname}`: {'OK' if info['match'] else 'MISMATCH'}")
    L.append("")

    # ── P1A.1–P1A.3 ──
    L.append("## P1A.1–P1A.3: Unchanged")
    L.append("")
    L.append("**PASS**, unchanged, not rerun.")
    L.append("The invariant-subspace estimator remains qualified from the prior pass.")
    L.append("")

    # ── P1A.4a ──
    L.append("## P1A.4a: Resolution Family Qualification")
    L.append("")

    # History
    L.append("### History: both prior contamination rulings withdrawn")
    L.append("")
    L.append("**Original per-sector contamination rulings: INVALID.**")
    L.append("The per-sector power-law envelope `E_ρ(n) = C_ρ n^(−α_ρ)` was nine "
             "independent calibrations. Roundoff-dominated sectors had noise fitted "
             "(exponents −2.7 to +10.7, prefactors 5e-15 to 4e+10); genuinely "
             "contaminated sectors got plausible exponents and wide envelopes. "
             "All nine rulings were withdrawn before any target execution.")
    L.append("")
    L.append("**First pooled replacement rulings: INVALID.**")
    L.append("The pooled fit (p=19.02, exp(r_max)≈17,300) ran on unqualified "
             "ladder data. The calibration set was never qualified as a common "
             "refinement family: `K_STENCIL_RULE` changed the discretization "
             "scheme across most of the sequence, and extraction failures "
             "(||P_spec||=0) entered downstream calculations.")
    L.append("")
    L.append("**Cause:** `K_STENCIL_RULE = min(110, max(20, int(n*120*0.015)))` "
             "varied k with n. Below n≈62 the discretization SCHEME changed along "
             "the ladder, not merely its resolution. At n=30, k=54 gave "
             "||L||=15,629 vs k=110 giving ||L||=1,415. The ladder varied TWO "
             "parameters; only above n≈62 did k saturate at 110.")
    L.append("")
    L.append("> The original P1A.4 sequence was not a qualified ONE-PARAMETER "
             "refinement family: seed count did not monotonically control "
             "geometric resolution, and the stencil-size rule changed the "
             "discrete scheme across most of the sequence.")
    L.append("")

    # Gate parameters
    L.append("### Cloud Admissibility Gate")
    L.append("")
    L.append(f"Fixed stencil: k={K_FIXED}, m={RBF_M}, p={RBF_P}")
    L.append(f"Reference density: {N_STABILITY_REF} Fibonacci points on S³")
    L.append("")
    L.append("Gate thresholds (frozen before running):")
    L.append(f"- mesh ratio h/q ≤ {GATE['mesh_ratio_max']}")
    L.append(f"- stencil conditioning max ≤ {GATE['stencil_cond_max']:.0e}")
    L.append(f"- constant-function residual ||L·1||∞ ≤ {GATE['const_residual_max']:.0e}")
    L.append(f"- polynomial reproduction residual ≤ {GATE['reprod_max']}")
    L.append("")

    gate_results = results.get("gate_results", {})
    admissible_ns = results.get("admissible_ns", [])

    L.append("| n | N | h | q | h/q | cond_max | reprod | const_res "
             "| h²||L|| | PASS |")
    L.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for n in CALIBRATION_CANDIDATES:
        d = gate_results.get(n, {})
        L.append(
            f"| {n} | {d.get('N', '')} "
            f"| {d.get('h', 0):.4f} "
            f"| {d.get('q', 0):.6f} "
            f"| {d.get('mesh_ratio', 0):.1f} "
            f"| {d.get('stencil_cond_max', 0):.2e} "
            f"| {d.get('reprod_residual', 0):.2e} "
            f"| {d.get('const_residual', 0):.2e} "
            f"| {d.get('h2_L', 0):.2f} "
            f"| {'**PASS**' if d.get('pass') else 'FAIL'} |")
    L.append("")
    L.append(f"**Admissible:** {len(admissible_ns)}/{len(CALIBRATION_CANDIDATES)} "
             f"clouds passed the gate.")
    L.append(f"Admissible n: {admissible_ns}")
    L.append("")

    # Fill distance table
    fill_dists = results.get("fill_distances", {})
    if fill_dists:
        sorted_by_h = sorted(fill_dists.items(), key=lambda x: x[1])
        L.append("### Admissible family — fill distance (sorted by h)")
        L.append("")
        L.append("| n | h (geodesic) |")
        L.append("| --- | --- |")
        for n, h in sorted_by_h:
            L.append(f"| {n} | {h:.4f} |")
        L.append("")
        h_min = sorted_by_h[0][1]
        h_max = sorted_by_h[-1][1]
        L.append(f"h range: [{h_min:.4f}, {h_max:.4f}]")
        L.append("")

    # Sector extraction data
    conv = results.get("convergence", {})
    if conv and admissible_ns:
        sorted_ns = sorted(admissible_ns,
                           key=lambda n: fill_dists.get(n, 0))

        for tbl_title, tbl_key, tbl_fmt in [
            ("eps_SA (self-adjointness residual)", "eps_SA", ".4e"),
            ("Cluster position error |λ − λ_analytic|",
             "cluster_position", None),
            ("Cluster spread", "cluster_spread", ".4e"),
            ("im_max (imaginary contamination)", "im_max", ".2e"),
            ("P_spec norm", "P_spec_norm", ".1f"),
            ("||L||₂", "L_norm", ".1f"),
        ]:
            L.append(f"### {tbl_title}")
            L.append("")
            hdr = "| Sector |" + "".join(
                f" {n} |" for n in sorted_ns)
            L.append(hdr)
            L.append("| --- |" + " --- |" * len(sorted_ns))
            for label in LABELS:
                vals = []
                for n in sorted_ns:
                    raw = conv.get((n, label), {}).get(tbl_key)
                    if raw is None:
                        vals.append("N/A")
                    elif tbl_key == "cluster_position":
                        idx_l = LABELS.index(label)
                        lam_a = MCKAY_DIST[idx_l] * (MCKAY_DIST[idx_l] + 2)
                        vals.append(f"{abs(raw - lam_a):.2e}")
                    else:
                        vals.append(f"{raw:{tbl_fmt}}")
                L.append(f"| {label} | " + " | ".join(vals) + " |")
            L.append("")

        L.append("### Extraction status")
        L.append("")
        L.append("| Sector |" + "".join(
            f" {n} |" for n in sorted_ns))
        L.append("| --- |" + " --- |" * len(sorted_ns))
        for label in LABELS:
            vals = []
            for n in sorted_ns:
                e = conv.get((n, label), {})
                if e.get("extraction_fail"):
                    vals.append("**FAIL**")
                elif e.get("extraction_ok"):
                    vals.append("OK")
                else:
                    vals.append("?")
            L.append(f"| {label} | " + " | ".join(vals) + " |")
        L.append("")

    # ── Contamination qualification ──
    contam = results.get("contamination", {})
    cal_evals = contam.get("calibration_evals", {})
    fit = contam.get("fit", {})

    L.append("## Contamination Qualification")
    L.append("")

    # Floor classification
    L.append("### Floor classification (calibration)")
    L.append("")
    L.append("```")
    L.append("F = eps_mach * ||L||₂ * max(||P_spec||₂, 1)")
    L.append("EXTRACTION_FAIL if ||P_spec|| < 0.5")
    L.append("```")
    L.append("")
    L.append("| Sector | n | h | ||L||₂ | P_spec | Floor | im_max | Regime |")
    L.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for n in (sorted(admissible_ns, key=lambda x: fill_dists.get(x, 0))
              if admissible_ns else []):
        for label in LABELS:
            ev = cal_evals.get((n, label), {})
            if not ev or ev.get("regime") == "ERROR":
                continue
            floor_s = (f"{ev['floor']:.2e}" if ev.get("floor") is not None
                       else "N/A")
            L.append(
                f"| {label} | {n} "
                f"| {ev.get('h', 0):.4f} "
                f"| {_fmt(ev.get('L_norm_2', 'N/A'), '.1f')} "
                f"| {_fmt(ev.get('P_spec_norm_used', 'N/A'), '.1f')} "
                f"| {floor_s} "
                f"| {ev.get('im_max', 0):.2e} "
                f"| {ev.get('regime', 'N/A')} |")
    L.append("")

    n_floor = sum(1 for v in cal_evals.values()
                  if v.get("regime") == "FLOOR_LIMITED")
    n_above = sum(1 for v in cal_evals.values()
                  if v.get("regime") == "ABOVE_FLOOR")
    n_exfail = sum(1 for v in cal_evals.values()
                   if v.get("regime") == "EXTRACTION_FAIL")
    L.append(f"**Summary:** {n_floor} FLOOR_LIMITED, {n_above} ABOVE_FLOOR, "
             f"{n_exfail} EXTRACTION_FAIL")
    L.append("")

    # Pooled fit
    L.append("### Pooled convergence fit")
    L.append("")
    L.append("```")
    L.append("s_ρ = 1 + d_ρ(d_ρ + 2)       (analytic scale)")
    L.append("y = im_max / s_ρ              (normalized contamination)")
    L.append("y(h) ~ A h^p                  (one A, one p)")
    L.append("E_ρ(h) = s_ρ A h^p exp(r_max) (envelope)")
    L.append("```")
    L.append("")
    if fit.get("status") == "ok":
        L.append(f"**A** = {fit['A']:.4e}, **p** = {fit['p']:.2f}, "
                 f"**r_max** = {fit['r_max']:.4f}, "
                 f"**exp(r_max)** = {np.exp(fit['r_max']):.2f}, "
                 f"**n_points** = {fit['n_points']}")
        if fit["p"] <= 0:
            L.append("")
            L.append("**WARNING: p ≤ 0 — pooled convergence model failed.**")
    else:
        L.append(f"**Fit status:** {fit.get('status', 'unknown')}")
    L.append("")

    # Above-floor residuals
    cal_above = [(k, v) for k, v in sorted(cal_evals.items())
                 if v.get("regime") == "ABOVE_FLOOR"]
    if cal_above and fit.get("status") == "ok":
        L.append("### Calibration: above-floor points")
        L.append("")
        L.append("| Sector | n | h | y = I/s | log-residual "
                 "| within envelope |")
        L.append("| --- | --- | --- | --- | --- | --- |")
        tags = fit.get("tags", [])
        residuals = fit.get("residuals", [])
        tag_resid = dict(zip(tags, residuals)) if tags else {}
        for (n, label), ev in cal_above:
            s = ev.get("s_rho", 1)
            im = ev.get("im_max", 0)
            y = im / s if s > 0 else 0
            tag = f"{label}@{n}"
            r = tag_resid.get(tag, "N/A")
            within = ev.get("convergence", "N/A")
            L.append(
                f"| {label} | {n} "
                f"| {ev.get('h', 0):.4f} "
                f"| {y:.4e} "
                f"| {_fmt(r, '.4f') if isinstance(r, float) else r} "
                f"| {within} |")
        L.append("")

    # Holdout
    holdout_evals = contam.get("holdout_evals", {})
    holdout_fd = contam.get("holdout_fill_distances", {})
    finer_holdout = contam.get("finer_holdout")

    all_holdout_ns = list(HOLDOUT_FIXED)
    if finer_holdout:
        all_holdout_ns.append(finer_holdout["n"])

    h_min_cal = min(fill_dists.values()) if fill_dists else float('inf')

    L.append("### Holdout evaluation")
    L.append("")
    L.append("Holdout resolutions were predeclared before running. "
             "The pooled fit (A, p, r_max) was frozen on calibration data; "
             "no holdout data influenced the fit.")
    L.append("")

    for n_h in all_holdout_ns:
        h_h = holdout_fd.get(n_h)
        if h_h is not None and h_h >= h_min_cal:
            h_type = "interpolation (inside calibration h range)"
        elif h_h is not None:
            pct = (h_min_cal - h_h) / h_min_cal * 100
            if finer_holdout and n_h == finer_holdout["n"]:
                h_type = f"genuinely finer holdout ({pct:.1f}% finer)"
            elif pct < 5:
                h_type = f"shallow extrapolation ({pct:.1f}% finer)"
            else:
                h_type = f"extrapolation ({pct:.1f}% finer)"
        else:
            h_type = "gate failed"

        L.append(f"**n = {n_h}** — {h_type}")
        if h_h is not None:
            L.append(f"  h = {h_h:.4f}")
        L.append("")
        L.append("| Sector | im_max | Floor | Regime | Envelope "
                 "| I/E | Eligibility |")
        L.append("| --- | --- | --- | --- | --- | --- | --- |")
        for label in LABELS:
            ev = holdout_evals.get((n_h, label), {})
            if not ev:
                continue
            floor_s = (_fmt(ev.get("floor"), ".2e")
                       if ev.get("floor") is not None else "N/A")
            env_s = (_fmt(ev.get("envelope"), ".2e")
                     if ev.get("envelope") is not None else "N/A")
            ie_s = (_fmt(ev.get("I_over_E"), ".4f")
                    if ev.get("I_over_E") is not None else "N/A")
            L.append(
                f"| {label} "
                f"| {ev.get('im_max', 0):.2e} "
                f"| {floor_s} "
                f"| {ev.get('regime', 'N/A')} "
                f"| {env_s} "
                f"| {ie_s} "
                f"| {ev.get('eligibility', 'N/A')} |")
        L.append("")

    # Discrimination test
    discrim = contam.get("discrimination", {})
    L.append("### Discrimination test")
    L.append("")
    L.append("A deliberately contaminated record, exceeding the intended numerical "
             "regime, MUST be rejected. This proves the gate can go red without "
             "opening P1A.5.")
    L.append("")
    if discrim.get("status") == "ran":
        L.append(f"**Source:** {discrim['source']}")
        L.append(f"**Original im_max:** {discrim['original_im']:.2e}")
        L.append(f"**Envelope E:** {discrim['envelope']:.2e}")
        L.append(f"**Injected im_max:** {discrim['injected_im']:.2e} "
                 f"(10× envelope)")
        L.append(f"**I/E:** {discrim['injected_I_over_E']:.2f}")
        L.append(f"**Result:** {discrim['result_eligibility']}")
        L.append(f"**{'PASS — envelope rejected the injection' if discrim['pass'] else 'FAIL — envelope accepted the injection'}**")
    else:
        L.append(f"**Could not run:** {discrim.get('status', 'unknown')}")
    L.append("")

    # Production summary
    L.append("### Production summary (n = 60)")
    L.append("")
    discrim_pass = discrim.get("pass", False)
    prod_in_family = N_PROD in (admissible_ns or [])

    if not prod_in_family:
        L.append(f"**n={N_PROD} did not pass the cloud admissibility gate.** "
                 "No production ruling can be made.")
        L.append("")
    else:
        L.append("| Sector | Regime | Convergence | Holdout | Eligibility |")
        L.append("| --- | --- | --- | --- | --- |")
        sector_rulings = {}
        for label in LABELS:
            prod_ev = cal_evals.get((N_PROD, label), {})
            regime = prod_ev.get("regime", "N/A")
            conv_status = prod_ev.get("convergence", "N/A")

            holdout_ok = True
            holdout_detail = "OK"
            for n_h in all_holdout_ns:
                hev = holdout_evals.get((n_h, label), {})
                elig_h = hev.get("eligibility", "")
                if elig_h in ("NO_LABEL", "EXTRACTION_FAIL", "GATE_FAIL"):
                    holdout_ok = False
                    holdout_detail = f"{elig_h}@{n_h}"
                    break
                ie = hev.get("I_over_E")
                if isinstance(ie, float):
                    holdout_detail = f"I/E={ie:.4f}@{n_h}"

            if regime == "EXTRACTION_FAIL":
                elig = "EXTRACTION_FAIL"
            elif not discrim_pass:
                elig = "NO_LABEL"
            elif regime == "FLOOR_LIMITED":
                elig = "QUALIFIED" if holdout_ok else "NO_LABEL"
            elif fit.get("status") != "ok" or fit.get("p", 0) <= 0:
                elig = "NO_LABEL"
            elif not holdout_ok:
                elig = "NO_LABEL"
            else:
                elig = prod_ev.get("eligibility", "NO_LABEL")

            sector_rulings[label] = elig
            L.append(f"| {label} | {regime} | {conv_status} "
                     f"| {holdout_detail} | **{elig}** |")
        L.append("")

    # P1A.5
    L.append("## P1A.5: Manufactured Calibration")
    L.append("")
    L.append("**Uninspected, deliberately.**")
    L.append("Planted effect sizes must not become the scale against which "
             "\"acceptable\" contamination is unconsciously set.")
    L.append("")

    # Status table
    L.append("## Status Record")
    L.append("")
    L.append("| Item | Status |")
    L.append("| --- | --- |")
    L.append("| P1A.0 to P1A.3 | **PASS**, unchanged, estimator still qualified |")
    L.append("| P1A.4 raw measurements | retained |")
    L.append("| Original per-sector contamination rulings | **INVALID** |")
    L.append("| First pooled replacement rulings | **INVALID** |")
    L.append("| Cause | calibration set never qualified as common refinement "
             "family; K_STENCIL_RULE changed discretization scheme; "
             "extraction failures entered downstream |")
    L.append("| P1A.5 | uninspected, deliberately |")
    L.append("")

    # Final verdicts
    L.append("## Final Verdicts")
    L.append("")

    if not prod_in_family:
        L.append("**Cannot issue verdicts:** n=60 failed the cloud gate.")
    elif not discrim_pass:
        L.append("**Cannot issue verdicts:** discrimination test failed. "
                 "The pooled envelope cannot reject injected contamination, "
                 "so it proves nothing.")
    else:
        qualified = [lb for lb in LABELS
                     if sector_rulings.get(lb) == "QUALIFIED"]
        not_qualified = [lb for lb in LABELS
                         if sector_rulings.get(lb) != "QUALIFIED"]
        L.append("**Sector eligibility (P1A.4a):**")
        L.append("")
        for label in LABELS:
            idx_l = LABELS.index(label)
            d_rho = MCKAY_DIST[idx_l]
            elig = sector_rulings.get(label, "UNKNOWN")
            L.append(f"- **{label}** (d_ρ={d_rho}): {elig}")
        L.append("")
        if qualified:
            L.append(f"**Qualified sectors:** {', '.join(qualified)}")
        if not_qualified:
            L.append(f"**Not qualified / no label:** "
                     f"{', '.join(not_qualified)}")
        L.append("")
        if len(qualified) == len(LABELS):
            L.append("**P1A.4a: ALL NINE SECTORS QUALIFIED**")
        elif len(qualified) >= 7:
            L.append(f"**P1A.4a: {len(qualified)} OF {len(LABELS)} "
                     "SECTORS QUALIFIED**")
        else:
            L.append("**P1A.4a: QUALIFICATION INCOMPLETE**")

    return "\n".join(L)


# ---------------------------------------------------------------------------
#  main
# ---------------------------------------------------------------------------

def main():
    t0 = time.time()
    results = {}
    drho_map = _d_rho_map()

    # ===== P1A.0 =====
    print("P1A.0: Verifying P0 freeze...")
    p0_pass, p0_details = verify_p0_frozen()
    results["p1a_0"] = {"pass": p0_pass, "details": p0_details}
    print(f"  P1A.0: {'PASS' if p0_pass else 'FAIL'}")
    if not p0_pass:
        print("  STOP: P0 manifest mismatch.")
        elapsed = time.time() - t0
        note = produce_note(results, elapsed)
        _write_note(note)
        return 1

    # ===== P1A.1–P1A.3: unchanged =====
    print("\nP1A.1-P1A.3: PASS, unchanged, not rerun.")

    # ===== Build infrastructure =====
    print("\nBuilding group and representations...")
    elems = build_icosians()
    chi = build_character_table(elems)
    reps, bases = build_all_representations(elems, chi)

    # ===== Phase 1: Cloud Admissibility Gate =====
    print(f"\n{'='*60}")
    print(f"P1A.4a Phase 1: Cloud Admissibility Gate (k={K_FIXED})")
    print(f"  Candidates: {CALIBRATION_CANDIDATES}")
    print(f"{'='*60}")

    gate_results = {}
    admissible_ns = []
    clouds = {}
    for n_seeds in CALIBRATION_CANDIDATES:
        print(f"\n  --- n={n_seeds} ({n_seeds*120} nodes) ---")
        seeds = fibonacci_seeds_s3(n_seeds)
        X, oid, gid = build_orbit_cloud(seeds, elems)
        passed, diag = cloud_gate(X, oid, gid, elems, reps["R0"],
                                   k=K_FIXED, m=RBF_M, p=RBF_P,
                                   print_fn=lambda s: print(f"  {s}"))
        gate_results[n_seeds] = diag
        if passed:
            admissible_ns.append(n_seeds)
            W = build_Mh_base(X, oid, n_seeds)
            clouds[n_seeds] = (X, oid, gid, W)
            print("  GATE: PASS")
        else:
            print(f"  GATE: FAIL "
                  f"({', '.join(diag.get('fail_reasons', []))})")

    print(f"\nAdmissible: {len(admissible_ns)}/{len(CALIBRATION_CANDIDATES)}")
    print(f"  n = {admissible_ns}")

    if not admissible_ns:
        print("STOP: no admissible calibration clouds.")
        results["gate_results"] = gate_results
        results["admissible_ns"] = []
        results["convergence"] = {}
        results["fill_distances"] = {}
        results["contamination"] = {}
        elapsed = time.time() - t0
        note = produce_note(results, elapsed)
        _write_note(note)
        return 1

    # ===== Phase 2: Sector Extraction =====
    print(f"\n{'='*60}")
    print(f"P1A.4a Phase 2: Sector Extraction (k={K_FIXED})")
    print(f"{'='*60}")

    convergence = {}
    fill_distances = {}
    for n_seeds in admissible_ns:
        print(f"\n  --- n={n_seeds} seeds ---")
        X, oid, gid, W = clouds[n_seeds]
        entries, h = _run_sector_extraction(X, oid, gid, W, n_seeds,
                                            elems, reps)
        fill_distances[n_seeds] = h
        for label, entry in entries.items():
            convergence[(n_seeds, label)] = entry

    clouds.clear()

    results["gate_results"] = gate_results
    results["admissible_ns"] = admissible_ns
    results["convergence"] = convergence
    results["fill_distances"] = fill_distances

    # ===== Contamination Qualification =====
    print(f"\n{'='*60}")
    print("Contamination Qualification")
    print(f"{'='*60}")

    cal_evals = {}
    above_floor_points = []
    for n_seeds in admissible_ns:
        h = fill_distances[n_seeds]
        for label in LABELS:
            key = (n_seeds, label)
            entry = convergence.get(key, {})

            if entry.get("extraction_fail"):
                cal_evals[key] = {
                    "regime": "EXTRACTION_FAIL",
                    "im_max": entry.get("im_max"), "h": h,
                    "L_norm_2": entry.get("L_norm"),
                    "P_spec_norm_used": entry.get("P_spec_norm"),
                    "eligibility": "EXTRACTION_FAIL",
                }
                continue

            im_max = entry.get("im_max")
            L_norm = entry.get("L_norm")
            P_spec_norm = entry.get("P_spec_norm")

            if im_max is None or L_norm is None:
                cal_evals[key] = {"regime": "ERROR"}
                continue

            d_rho = drho_map[label]
            s_rho = analytic_scale(d_rho)
            floor = numerical_floor(L_norm, P_spec_norm)

            if floor is None:
                cal_evals[key] = {
                    "regime": "EXTRACTION_FAIL",
                    "floor": None, "im_max": im_max, "h": h,
                    "L_norm_2": L_norm,
                    "P_spec_norm_used": P_spec_norm,
                    "s_rho": s_rho,
                    "eligibility": "EXTRACTION_FAIL",
                }
                continue

            is_floor = im_max <= floor
            ev = {
                "regime": "FLOOR_LIMITED" if is_floor else "ABOVE_FLOOR",
                "floor": floor, "im_max": im_max,
                "L_norm_2": L_norm,
                "P_spec_norm_used": max(
                    float(P_spec_norm) if P_spec_norm else 1.0, 1.0),
                "s_rho": s_rho, "h": h,
            }
            cal_evals[key] = ev

            if not is_floor and im_max > 0:
                y = im_max / s_rho
                above_floor_points.append((h, y, f"{label}@{n_seeds}"))

    n_fl = sum(1 for v in cal_evals.values()
               if v.get("regime") == "FLOOR_LIMITED")
    n_ab = sum(1 for v in cal_evals.values()
               if v.get("regime") == "ABOVE_FLOOR")
    n_ef = sum(1 for v in cal_evals.values()
               if v.get("regime") == "EXTRACTION_FAIL")
    print(f"  {n_fl} FLOOR_LIMITED, {n_ab} ABOVE_FLOOR, "
          f"{n_ef} EXTRACTION_FAIL")

    fit = pooled_convergence_fit(above_floor_points)
    print(f"  Pooled fit: status={fit.get('status')}, "
          f"A={fit.get('A', 'N/A')}, p={fit.get('p', 'N/A')}, "
          f"r_max={fit.get('r_max', 'N/A')}, "
          f"n={fit.get('n_points', 0)}")

    if fit.get("status") == "ok":
        for key, ev in cal_evals.items():
            if ev.get("regime") == "ABOVE_FLOOR":
                n_s, label = key
                d_rho = drho_map[label]
                s_rho = analytic_scale(d_rho)
                E = envelope_value(s_rho, fit["A"], fit["p"],
                                   fit["r_max"], ev["h"])
                within = ev["im_max"] <= E
                ev["envelope"] = E
                ev["I_over_E"] = (ev["im_max"] / E if E > 0
                                  else float('inf'))
                ev["convergence"] = "PASS" if within else "FAIL"
                ev["eligibility"] = "QUALIFIED" if within else "NO_LABEL"
            elif ev.get("regime") == "FLOOR_LIMITED":
                ev["convergence"] = "N/A_FLOOR"
                ev["eligibility"] = "QUALIFIED"

    contam_results = {"calibration_evals": cal_evals, "fit": fit}

    # ===== Holdouts =====
    print(f"\n{'='*60}")
    print("Holdout Evaluation")
    print(f"{'='*60}")

    holdout_evals = {}
    holdout_fill_dists = {}
    holdout_gate = {}

    h_min_cal = min(fill_distances.values())

    for n_holdout in HOLDOUT_FIXED:
        print(f"\n  --- holdout n={n_holdout} ---")
        gd, h_r, evals = _run_holdout(n_holdout, elems, reps,
                                       drho_map, fit)
        holdout_gate[n_holdout] = gd
        if h_r is not None:
            holdout_fill_dists[n_holdout] = h_r
            if h_r >= h_min_cal:
                print(f"  h={h_r:.4f} — interpolation")
            else:
                pct = (h_min_cal - h_r) / h_min_cal * 100
                print(f"  h={h_r:.4f} — {pct:.1f}% finer")
        holdout_evals.update(evals)

    # Finer holdout
    print(f"\n  --- Genuinely finer holdout search ---")
    h_target = FINER_HOLDOUT_C * h_min_cal
    print(f"  h_min_cal={h_min_cal:.4f}, target h ≤ {h_target:.4f}")
    print(f"  Search: {FINER_HOLDOUT_SEARCH}")

    finer_holdout = None
    for n_s in FINER_HOLDOUT_SEARCH:
        print(f"\n  n={n_s}:")
        seeds = fibonacci_seeds_s3(n_s)
        X, oid, gid = build_orbit_cloud(seeds, elems)
        gp, gd = cloud_gate(X, oid, gid, elems, reps["R0"],
                             k=K_FIXED, m=RBF_M, p=RBF_P,
                             print_fn=lambda s: print(f"    {s}"))
        h_c = gd["h"]
        print(f"    gate={'PASS' if gp else 'FAIL'}, h={h_c:.4f}")

        if gp and h_c <= h_target:
            finer_holdout = {"n": n_s, "h": h_c, "gate_diag": gd}
            W = build_Mh_base(X, oid, n_s)
            entries, h_r = _run_sector_extraction(
                X, oid, gid, W, n_s, elems, reps)
            holdout_fill_dists[n_s] = h_r
            pct = (h_min_cal - h_r) / h_min_cal * 100
            print(f"  FOUND finer holdout: n={n_s}, "
                  f"h={h_r:.4f} ({pct:.1f}% finer)")
            for label, entry in entries.items():
                if entry.get("extraction_fail"):
                    holdout_evals[(n_s, label)] = {
                        "regime": "EXTRACTION_FAIL",
                        "im_max": entry.get("im_max"), "h": h_r,
                        "eligibility": "EXTRACTION_FAIL"}
                    continue
                im = entry.get("im_max")
                ln = entry.get("L_norm")
                ps = entry.get("P_spec_norm")
                if im is None or ln is None:
                    holdout_evals[(n_s, label)] = {
                        "regime": "ERROR", "eligibility": "NO_LABEL"}
                    continue
                d_rho = drho_map[label]
                s_rho = analytic_scale(d_rho)
                ev = evaluate_entry(im, ln, ps, s_rho, h_r, fit)
                holdout_evals[(n_s, label)] = ev
                print(f"    {label}: im={im:.2e}, "
                      f"{ev['regime']}, {ev['eligibility']}")
            break

    if finer_holdout is None:
        print("  No genuinely finer holdout found.")

    contam_results["holdout_evals"] = holdout_evals
    contam_results["holdout_fill_distances"] = holdout_fill_dists
    contam_results["holdout_gate"] = holdout_gate
    contam_results["finer_holdout"] = finer_holdout

    # ===== Discrimination Test =====
    print(f"\n{'='*60}")
    print("Discrimination Test")
    print(f"{'='*60}")

    discrim = _run_discrimination_test(cal_evals, fit, drho_map)
    contam_results["discrimination"] = discrim
    results["contamination"] = contam_results

    # ===== P1A.5 =====
    print("\nP1A.5: uninspected, deliberately.")
    results["p1a_5"] = {"status": "uninspected"}

    # ===== Note =====
    elapsed = time.time() - t0
    note = produce_note(results, elapsed)
    _write_note(note)
    print(f"\nTotal elapsed: {elapsed:.1f}s")
    return 0


def _write_note(note):
    note_path = os.path.join(BASE_DIR, "P1A_QUALIFICATION_NOTE.md")
    with open(note_path, 'w') as f:
        f.write(note)
    print(f"\nQualification note written to: {note_path}")


if __name__ == "__main__":
    sys.exit(main())
