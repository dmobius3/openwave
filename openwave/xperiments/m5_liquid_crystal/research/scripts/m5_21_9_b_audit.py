"""M5.21.9 independent adversarial audit (part b: the run claims).

Audits CL1-CL9 of the M5.21.9 run (fixed-J isorotation electron +
Larmor acceptance) with the auditor's OWN analysis. Imports ONLY the
certified instruments (m5_21_3_a_4d.py, audited M5.21.3;
m5_21_8_b_lattice.py, audited M5.21.8). Nothing is imported from the
task scripts m5_21_9_c/d/e/f; their JSONs and npz endpoints are read
as data. All fits (unwrap + slope), all 4x4 background arithmetic,
the leapfrog reproduction, and the omega = 0 discriminator runs are
re-implemented here from the equations in the run note.

Sections:
    CL1  delta = -0.3 statics recompute (own points + flip locator)
    CL2  relax collapse + V4 flat/fires probes (own perturbations)
    CL3  fixed-J internal consistency + npz kin recompute
    CL4  Legendre arithmetic recompute (+ alt midpoint convention)
    CL5  gate drift recompute + leapfrog failure-mode reproduction
    CL6  floor + perturbative ladder refits (own unwrap + lsq)
    CL7  saturation map verification from the saved rows
    CL8  differential refits (full + windows) + background F blocks
         + static cross-coupling + t = 0 torque split
    CL9  kin convention recompute (probe variant vs conjugation)
    X1   omega = 0 discriminator: reversed field on the static state

Run:  python3 m5_21_9_b_audit.py            (all sections, ~10 min)
      python3 m5_21_9_b_audit.py fast       (skips CL5/X1 dynamics)
Out:  ../data/m5_21_9_audit.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_s3 = importlib.util.spec_from_file_location(
    "ins4", os.path.join(HERE, "m5_21_3_a_4d.py"))
INS4 = importlib.util.module_from_spec(_s3)
_s3.loader.exec_module(INS4)
_s8 = importlib.util.spec_from_file_location(
    "lat8", os.path.join(HERE, "m5_21_8_b_lattice.py"))
LAT8 = importlib.util.module_from_spec(_s8)
_s8.loader.exec_module(LAT8)

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
RES = {"meta": {"date": "2026-07-20", "auditor": "independent",
                "instruments": ["m5_21_3_a_4d.py (certified)",
                                "m5_21_8_b_lattice.py (certified)"],
                "task_scripts_imported": "none"}}


def jload(name):
    with open(os.path.join(DATA, name)) as f:
        return json.load(f)


def say(*a):
    print(*a, flush=True)


# ================= own fit machinery =================
def phi_series(rows):
    """recompute phi from the saved J components (not the saved
    phi_xy field), own unwrap."""
    t = np.array([r["t"] for r in rows], dtype=float)
    jx = np.array([r["J"]["x"] for r in rows], dtype=float)
    jy = np.array([r["J"]["y"] for r in rows], dtype=float)
    ok = np.isfinite(jx) & np.isfinite(jy)
    phi = np.unwrap(np.arctan2(jy[ok], jx[ok]))
    return t[ok], phi, np.sqrt(jx[ok] ** 2 + jy[ok] ** 2)


def slope(t, y):
    """own least squares linear fit: slope + rms residual."""
    A = np.stack([t, np.ones_like(t)], axis=1)
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ c
    return float(c[0]), float(np.sqrt(np.mean(res ** 2)))


def fit_run(name, windows=((None, None), (None, 5.0), (5.0, None))):
    r = jload(name)
    t, phi, jxy = phi_series(r["rows"])
    out = {"eps": r["eps"], "B_meas": r["B_meas"],
           "om_start": r.get("om_start"),
           "jxy_first": float(jxy[0]), "jxy_last": float(jxy[-1]),
           "E_first": r["rows"][0]["E"], "E_last": r["rows"][-1]["E"]}
    for (a, b) in windows:
        m = np.ones_like(t, dtype=bool)
        if a is not None:
            m &= t >= a
        if b is not None:
            m &= t <= b
        sl, rms = slope(t[m], phi[m])
        if a is None and b is None:
            key = "Omega_full"
            out["fit_rms_rad"] = rms
        elif a is None:
            key = f"Omega_t<={b:g}"
        elif b is None:
            key = f"Omega_t>={a:g}"
        else:
            key = f"Omega_{a:g}_{b:g}"
        out[key] = sl
    return out


W3 = ((None, None), (0.0, 3.5), (3.0, 6.5), (6.5, 10.0))
W3K = ("Omega_0_3.5", "Omega_3_6.5", "Omega_6.5_10")


# ================= CL1 =================
def cl1():
    say("== CL1: delta = -0.3 statics ==")
    mc = jload("m5_21_9_lat_mcurve_dl-0.3_g8_n32.json")
    rows = mc["rows"]
    cfg = INS4.base_cfg(s=-1.0, g=8.0, n=32, L=48.0, delta=-0.3)
    own = {}
    for m in (0.0, 0.15, -0.15, 0.175, -0.175):
        M = LAT8.dressed(cfg, m)
        eu, ev = INS4.e_parts(M, cfg)
        k = INS4.kin_of(M, LAT8.a0_unit(cfg, m), cfg)
        own[f"{m:+.3f}"] = {"E_u": float(eu), "E_v": float(ev),
                            "kin": float(k)}
        say(f"  own m {m:+.3f} E_u {eu:.6f} E_v {ev:.2e} "
            f"kin {k:+.4f}")
    # compare against the saved rows
    bym = {round(r["m"], 3): r for r in rows}
    dmax = 0.0
    for m in (0.0, 0.15, -0.15, 0.175, -0.175):
        r = bym[round(m, 3)]
        o = own[f"{m:+.3f}"]
        dmax = max(dmax,
                   abs(o["E_u"] - r["E_u"]) / max(abs(r["E_u"]), 1),
                   abs(o["kin"] - r["kin"]) / max(abs(r["kin"]), 1))
    # evenness from the saved curve
    ev_gap = max(abs(rows[i]["E_u"] - rows[len(rows) - 1 - i]["E_u"])
                 / max(abs(rows[i]["E_u"]), 1e-12)
                 for i in range(len(rows) // 2))
    evmax = max(r["E_v"] for r in rows)
    # own parabolic vertex around the argmin
    ms = np.array([r["m"] for r in rows])
    Es = np.array([r["E_u"] + r["E_v"] for r in rows])
    i = int(np.argmin(Es))
    a, b, c = Es[i - 1], Es[i], Es[i + 1]
    mfit = ms[i] - 0.5 * (ms[1] - ms[0]) * (c - a) / (c - 2 * b + a)
    # kin sign-flip locator: own recompute between the grid points
    flip = {}
    for m in (0.155, 0.16, 0.165):
        M = LAT8.dressed(cfg, m)
        flip[f"{m:g}"] = float(
            INS4.kin_of(M, LAT8.a0_unit(cfg, m), cfg))
    say("  kin near flip:", flip)
    RES["CL1"] = {
        "verdict": "NUANCE",
        "claimed": {"E_u_mstar": -6.54, "kin_mstar": 39.9,
                    "flip": "past |m| ~ 0.16", "E_v": "~1e-22",
                    "even": True, "min_at": "|m| ~ 0.15"},
        "measured": {
            "own_points": own,
            "own_vs_json_relgap_max": float(dmax),
            "E_u_at_0.15": own["+0.150"]["E_u"],
            "kin_at_0.15": own["+0.150"]["kin"],
            "even_relgap_max": float(ev_gap),
            "E_v_family_max": float(evmax),
            "parab_vertex_own": float(abs(mfit)),
            "kin_near_flip": flip},
        "notes": [
            "own recompute matches the saved curve to machine "
            "precision; minima, sub-vacuum E_u, kin(+39.9), "
            "evenness (~1e-14 relative), E_v machine zero: all "
            "hold",
            "THE NUANCE: kin(0.155) = -7.4 already NEGATIVE, so "
            "the sign flip sits at |m| ~ 0.154 (linear interp), "
            "NOT 'past ~0.16' as claimed. The flip is only "
            "~0.003-0.015 in m above the minimum (parab vertex "
            "0.139 grid, analytic 0.1508): the 'kin positive at "
            "the minimum' margin is thin, and at the ANALYTIC "
            "minimum 0.1508 the interpolated kin is ~+32, still "
            "positive but a factor 23 below the m = 0 value"]}


# ================= CL2 =================
def cl2():
    say("== CL2: relax collapse + V4 probes ==")
    rn = jload("m5_21_9_lat_relax_dl-0.3_m-0.150_none0.json")
    rt = jload("m5_21_9_lat_relax_dl-0.3_m-0.150_tanh3.json")
    cfg = INS4.base_cfg(s=-1.0, g=8.0, delta=-0.3)
    M = LAT8.dressed(cfg, -0.15)
    eu0, ev0 = INS4.e_parts(M, cfg)
    c = 10  # interior off-center cell
    # probe A: local boost dressing at one cell (spectrum-preserving)
    r = 0.3
    K = np.zeros((4, 4)); K[0, 1] = K[1, 0] = 1.0
    K2 = np.zeros((4, 4)); K2[0, 0] = K2[1, 1] = 1.0
    Q = np.eye(4) + np.sinh(r) * K + (np.cosh(r) - 1.0) * K2
    MA = M.copy()
    MA[c, 16, 16] = Q @ MA[c, 16, 16] @ Q.T
    euA, evA = INS4.e_parts(MA, cfg)
    # probe B: additive off-spectrum bump at the same cell
    MB = M.copy()
    MB[c, 16, 16] = MB[c, 16, 16] + 0.3 * np.diag([0, 1.0, 0, 0])
    euB, evB = INS4.e_parts(MB, cfg)
    say(f"  E_v base {ev0:.3e} boost-dressed {evA:.3e} "
        f"off-spectrum {evB:.3e}")
    RES["CL2"] = {
        "verdict": "CONFIRMED",
        "claimed": {"free_dies": "non-finite, E_u ~ -6.3e266",
                    "tanh_dies_identically": True,
                    "boost_channel": "V4-flat on the boosted family",
                    "V4_fires_off_spectrum": True},
        "measured": {
            "none_stop": rn["stop"], "none_E_u_end": rn["E_u_end"],
            "tanh_stop": rt["stop"], "tanh_E_u_end": rt["E_u_end"],
            "identical_relgap": float(
                abs(rn["E_u_end"] - rt["E_u_end"])
                / abs(rn["E_u_end"])),
            "E_v_base": float(ev0),
            "E_v_boost_dressed_cell": float(evA),
            "E_v_off_spectrum_cell": float(evB),
            "E_u_shift_boost_cell": float(euA - eu0)},
        "notes": [
            "both relaxes stop 'non-finite' with E_u ~ -6.33e266; "
            "the two endpoints agree to 8e-9 relative: identical "
            "at the physics level, not bit-identical",
            "the JSON field finite=true refers to the M array (the "
            "ENERGY overflowed first); the note's 'dies non-finite' "
            "is the FIRE stop reason: wording holds",
            "local boost dressing of one cell leaves E_v at machine "
            "zero while an off-spectrum bump at the same cell fires "
            "V4 by ~17 orders: the boost channel is V4-flat and the "
            "collapse attribution is supported"]}


# ================= CL3 =================
def cl3():
    say("== CL3: fixed-J internal consistency ==")
    recs = {om: jload(f"m5_21_9_fixedj_om{om:g}.json")
            for om in (0.2, 0.5, 1.0)}
    rows = {}
    for om, r in recs.items():
        st, fn = r["start"], r["final"]
        J_chk = 2.0 * st["kin0"] * st["om_target"]
        omf_chk = st["J"] / (2.0 * fn["kin_final"])
        spec_d = max(abs(a - b) for a, b in
                     zip(st["spec_core0"], fn["spec_core"]))
        rows[f"om{om:g}"] = {
            "J": st["J"],
            "J_eq_2kin0om_gap": abs(st["J"] - J_chk),
            "omf_json": fn["omega_star_final"],
            "omf_eq_J_over_2kin_gap": abs(
                fn["omega_star_final"] - omf_chk),
            "kin_final": fn["kin_final"],
            "rel_move": fn["rel_move_from_start"],
            "spec_absdrift_max": float(spec_d),
            "iters": r["hist"][-1]["iter"]}
        say(f"  om {om}: omf {fn['omega_star_final']:.5f} "
            f"kin {fn['kin_final']:.5f} "
            f"rel_move {fn['rel_move_from_start']:.3e}")
    # own kin recompute on the om 0.5 endpoint npz
    Z = np.load(os.path.join(DATA, "m5_21_9_fixedj_om0.5_end.npz"))
    Mnp = Z["M"].astype(np.float64)
    _, cfg = INS4.load_p1("p1_s-1")
    a0 = INS4.gen_catalog(cfg, Mnp)["clock_local"]
    kin_own = float(INS4.kin_of(Mnp, a0, cfg))
    say(f"  own kin(om0.5 npz) {kin_own:.6f} vs json "
        f"{recs[0.5]['final']['kin_final']:.6f}")
    kins = [rows[f"om{o:g}"]["kin_final"] for o in (0.2, 0.5, 1.0)]
    oms = [rows[f"om{o:g}"]["omf_json"] for o in (0.2, 0.5, 1.0)]
    RES["CL3"] = {
        "verdict": "NUANCE",
        "claimed": {"rel_move_max": 1.3e-4,
                    "spec": [0.026, 0.283, 0.993, 8.000],
                    "kin_range": [0.297, 0.304],
                    "om_final": {"0.2": 0.1998, "0.5": 0.4961,
                                 "1.0": 0.978}},
        "measured": {
            "rungs": rows,
            "kin_own_npz_om0.5": kin_own,
            "kin_own_vs_json_relgap": abs(
                kin_own - recs[0.5]["final"]["kin_final"])
            / recs[0.5]["final"]["kin_final"],
            "kin_monotone_up": bool(kins[0] < kins[1] < kins[2]),
            "om_monotone_down_vs_target": bool(
                oms[0] < 0.2 and oms[1] < 0.5 and oms[2] < 1.0)},
        "notes": [
            "J = 2 kin0 om_target holds exactly; omega*_final = "
            "J/(2 kin_final) holds exactly; rel_move max 1.25e-4 "
            "<= 1.3e-4; core spectrum drift <= 6e-4 absolute",
            "own kin on the om0.5 endpoint npz reproduces kin_final "
            "to float32-quantization level",
            "THE NUANCE: the claimed omega*(0.2) = 0.1998 is the "
            "STALE first-pass (300-iteration) figure; the on-disk "
            "1200-iteration rerun gives 0.19923. The 0.5 and 1.0 "
            "figures are current. The direction (omega* falls, kin "
            "rises with J) is confirmed"]}


# ================= CL4 =================
def cl4():
    say("== CL4: Legendre recompute ==")
    recs = [jload(f"m5_21_9_fixedj_om{om:g}.json")
            for om in (0.2, 0.5, 1.0)]
    Js, Es, oms, evs = [], [], [], []
    for r in recs:
        J = r["start"]["J"]
        fn = r["final"]
        E = fn["E_u"] + fn["E_v"] + J * J / (4.0 * fn["kin_final"])
        Js.append(J); Es.append(E)
        oms.append(fn["omega_star_final"]); evs.append(fn["E_v"])
    ivals = []
    for i in (0, 1):
        dEdJ = (Es[i + 1] - Es[i]) / (Js[i + 1] - Js[i])
        om_mid = 0.5 * (oms[i] + oms[i + 1])
        # alternative midpoint: omega interpolated at the J midpoint
        Jm = 0.5 * (Js[i] + Js[i + 1])
        om_alt = np.interp(Jm, Js, oms)
        ivals.append({"dEdJ": float(dEdJ),
                      "om_mid": float(om_mid),
                      "ratio": float(dEdJ / om_mid),
                      "ratio_altmid": float(dEdJ / om_alt)})
        say(f"  interval {i}: dEdJ {dEdJ:.4f} om_mid {om_mid:.4f} "
            f"ratio {dEdJ / om_mid:.4f}")
    iters = [r["hist"][-1]["iter"] for r in recs]
    ev_tail = [r["hist"][-1]["E_v"] - r["hist"][-2]["E_v"]
               for r in recs]
    RES["CL4"] = {
        "verdict": "CONFIRMED",
        "claimed": {"ratio_lower": 0.996, "ratio_upper": 0.989,
                    "iters": 1200},
        "measured": {"E_tot": Es, "J": Js, "omega_star": oms,
                     "intervals": ivals, "iters": iters,
                     "E_v_last_window_drift": ev_tail},
        "notes": [
            "arithmetic reproduces 0.996 / 0.989 exactly from the "
            "JSONs with the stated E_tot definition, applied "
            "consistently across rungs",
            "alt midpoint convention (omega at the J midpoint) "
            "moves the ratios by < 0.1%: convention-robust",
            "fragility: all rungs stop at max_iter with E_v still "
            "declining ~0.015 per 300 iterations; the ~1% closure "
            "leans on the equal-depth cancellation of that common "
            "drift (E_v spread across rungs 0.0017, i.e. ~1-3% of "
            "the energy differences). Disclosed in the note"]}


# ================= CL5 =================
def leap_own(M, Mt, cfg, freew, gamw, dt, mask=True):
    """own leapfrog from the note's equations (M5.21.6 form)."""
    h3 = cfg["h"] ** 3
    F = -INS4.grad(M, cfg) / h3
    Mt = Mt + 0.5 * dt * F
    if mask:
        Mt = Mt * freew
    Mt = Mt / (1.0 + 0.5 * dt * gamw)
    M = INS4.sym4(M + dt * Mt)
    F = -INS4.grad(M, cfg) / h3
    Mt = Mt + 0.5 * dt * F
    if mask:
        Mt = Mt * freew
    Mt = Mt / (1.0 + 0.5 * dt * gamw)
    return M, Mt


def gamma_ramp(cfg, g_max=0.5):
    X, Y, Z = INS4.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    r0, r1 = 0.70 * cfg["L"] / 2, 0.95 * cfg["L"] / 2
    t = np.clip((r - r0) / (r1 - r0), 0.0, 1.0)
    return (g_max * t * t)[..., None, None]


def etot_own(M, Mt, cfg):
    eu, ev = INS4.e_parts(M, cfg)
    return float(eu + ev) + 0.5 * cfg["h"] ** 3 * float(
        np.sum(Mt * Mt))


def setup_own(eps, rev=False, env=True, om=None):
    """own setup mirroring the note's spec, on the current om0.2
    endpoint (the re-relaxed base)."""
    Z = np.load(os.path.join(DATA, "m5_21_9_fixedj_om0.2_end.npz"))
    M = Z["M"].astype(np.float64)
    _, cfg = INS4.load_p1("p1_s-1")
    rec = jload("m5_21_9_fixedj_om0.2.json")
    om = rec["final"]["omega_star_final"] if om is None else om
    a0 = INS4.gen_catalog(cfg, M)["clock_local"]
    MB = my_background(cfg, eps, rev)
    if env and eps != 0.0:
        MB = MB * INS4.envelope(cfg)[..., None, None]
    Mtot = INS4.sym4(M + MB)
    free = ~INS4.pin_shell(cfg["n"], cfg["h"])
    freew = free[..., None, None].astype(float)
    Mt = om * a0 * freew
    return Mtot, Mt, cfg, freew


def my_background(cfg, eps, rev=False):
    """own build: fwd = eps (x T_y - y T_x); rev = eps (x T_y +
    y T_x). T_i = E_0i + E_i0."""
    X, Y, Z = INS4.coords(cfg["n"], cfg["h"])
    Tx = np.zeros((4, 4)); Tx[0, 1] = Tx[1, 0] = 1.0
    Ty = np.zeros((4, 4)); Ty[0, 2] = Ty[2, 0] = 1.0
    s = 1.0 if rev else -1.0
    return eps * (X[..., None, None] * Ty
                  + s * Y[..., None, None] * Tx)


def cl5(run_dynamics=True):
    say("== CL5: gates ==")
    g = jload("m5_21_9_larmor_gates.json")
    e0, e1 = g["E0"] + g["KE0"], g["E1"] + g["KE1"]
    drift = abs(e1 - e0) / abs(e0)
    j0 = np.array([g["j0"][k] for k in "xyz"])
    j1 = np.array([g["j1"][k] for k in "xyz"])
    cosang = float(j0 @ j1 / (np.linalg.norm(j0)
                              * np.linalg.norm(j1)))
    ang = float(np.arccos(np.clip(cosang, -1, 1)))
    say(f"  drift {drift:.3e}  dir angle {ang:.3e} rad")
    repro = {}
    if run_dynamics is not False:
        done = {} if run_dynamics is True else dict(run_dynamics)
        # failure-mode reproduction on the CURRENT base, gamma = 0.
        # modes: mask = the certified form; free = no pinning at
        # all (a DIFFERENT conservative system: also conserves);
        # repin = M held on the shell but velocity NOT masked (the
        # plausible try-2 bug form: force-kicks pump KE forever)
        for (label, dt, steps, mode) in (
                ("masked_dt0.005_s400", 0.005, 400, "mask"),
                ("masked_dt0.02_s200", 0.02, 200, "mask"),
                ("unmasked_dt0.005_s400", 0.005, 400, "free"),
                ("unmasked_dt0.02_s200", 0.02, 200, "free"),
                ("repin_dt0.005_s400", 0.005, 400, "repin")):
            if label in done:
                repro[label] = done[label]
                continue
            M, Mt, cfg, freew = setup_own(0.0)
            pinw = 1.0 - freew
            Mpin = M * pinw
            gam0 = np.zeros_like(freew)
            E0 = etot_own(M, Mt, cfg)
            ok = True
            for k in range(steps):
                M, Mt = leap_own(M, Mt, cfg, freew, gam0, dt,
                                 mask=(mode == "mask"))
                if mode == "repin":
                    M = M * freew + Mpin
                if not np.isfinite(np.sum(Mt)):
                    ok = False
                    break
            d = abs(etot_own(M, Mt, cfg) - E0) / abs(E0) if ok \
                else float("inf")
            repro[label] = {"drift": d, "finite": ok}
            say(f"  repro {label}: drift {d:.3e} finite {ok}")
    RES["CL5"] = {
        "verdict": "NUANCE",
        "claimed": {"drift": 2.2e-8, "steps": 400, "dt": 0.005,
                    "dir_hold_pct": 1.0,
                    "try1": "dt=0.02 instability",
                    "try2": "unmasked-velocity bug"},
        "measured": {"drift_recomputed": float(drift),
                     "dir_angle_rad": ang,
                     "dir_change_frac": ang,
                     "j_norm_change_frac": float(
                         np.linalg.norm(j1) / np.linalg.norm(j0)
                         - 1.0),
                     "failure_mode_reproduction": repro},
        "notes": [
            "the gate numbers hold: drift recomputes to 2.23e-8 "
            "from the saved JSON, the auditor's own leapfrog "
            "reproduces 2.29e-8 on the production base, and the "
            "J direction moved 4.2e-4 rad (0.04%) < 1%",
            "the try-2 attribution is CONFIRMED QUANTITATIVELY: "
            "the repin form (shell held, velocity unmasked = the "
            "KE pump) reproduces 12.48% drift over 400 steps at "
            "dt = 0.005 vs the logged 12.5%",
            "THE NUANCE (try-1 attribution): the auditor's "
            "MASKED integrator at dt = 0.02 conserves to 1.6e-6 "
            "over 200 steps: dt = 0.02 is NOT unstable for this "
            "system, so the try-1 'dt instability' attribution "
            "is unsupported; try 1 was most plausibly the same "
            "masking bug at a coarser dt. The production choice "
            "dt = 0.005 is harmless but rests on a misdiagnosis",
            "a fully UNMASKED integrator also conserves (it is a "
            "different conservative system, free boundary): "
            "energy drift alone cannot distinguish it; only the "
            "mixed repin form pumps KE",
            "gates ran on the OLD base state (E0 6.7826 matches "
            "eps0_oldbase, not the re-relaxed 6.7409): the "
            "certificate transfers by construction, not by "
            "direct test, to the production base",
            "failed tries 1-2 are archived in the checkpoint "
            "only (49.8% / 12.5% not in data/)"]}


# ================= CL6 =================
def cl6():
    say("== CL6: floors + perturbative ladder ==")
    floors = {}
    for name in ("eps0_oldbase", "eps0v2"):
        floors[name] = fit_run(f"m5_21_9_larmor_{name}.json",
                               windows=W3)
        say(f"  {name}: Omega {floors[name]['Omega_full']:.4e} "
            f"rms {floors[name]['fit_rms_rad']:.2e} windows "
            + str([f"{floors[name][k]:.2e}" for k in W3K]))
    ladder = {}
    for name in ("eps0.002", "eps0.005", "eps0.01"):
        ladder[name] = fit_run(f"m5_21_9_larmor_{name}.json",
                               windows=W3)
        base = "eps0_oldbase" if \
            abs(ladder[name]["om_start"]
                - floors["eps0_oldbase"]["om_start"]) < 1e-12 \
            else "eps0v2"
        ladder[name]["floor_used"] = base
        ladder[name]["dev_from_floor"] = (
            ladder[name]["Omega_full"] - floors[base]["Omega_full"])
        say(f"  {name}: Omega {ladder[name]['Omega_full']:.4e} "
            f"dev {ladder[name]['dev_from_floor']:+.2e} "
            f"(floor {base})")
    devmax = max(abs(v["dev_from_floor"]) for v in ladder.values())
    curv = floors["eps0v2"]["Omega_6.5_10"] \
        / floors["eps0v2"]["Omega_0_3.5"]
    RES["CL6"] = {
        "verdict": "NUANCE",
        "claimed": {"Omega0_old": 9.06e-4, "Omega0_new": 9.23e-4,
                    "ladder_within": 1.1e-5, "B_max": 1.4e-4,
                    "floor_shape": "clean linear slope"},
        "measured": {"floors": floors, "ladder": ladder,
                     "dev_from_floor_max": float(devmax),
                     "floor_late_over_early_slope": float(curv)},
        "notes": [
            "own unwrap + lsq refit reproduces both floors "
            "(9.063e-4 / 9.232e-4); no LINEAR (Larmor-sign) "
            "response is resolvable anywhere in the ladder: that "
            "part of the conclusion stands",
            "PARTIAL REFUTATION (the +/-1.1e-5 figure): compared "
            "to its OWN base's floor, eps0.01 deviates -2.54e-5, "
            "2.3x beyond the claimed bar; the deviations follow "
            "a clean eps^2 law (-8.8e-7 / -6.2e-6 / -2.54e-5 at "
            "eps 0.002 / 0.005 / 0.01, ratio ~4 per eps "
            "doubling) continuing to the env0.03 symmetric "
            "slowdown (-1.8e-4): the QUADRATIC slowdown is "
            "already resolvable inside the 'perturbative' "
            "ladder. The claimed null holds only under a "
            "base-mismatched floor comparison",
            "THE NUANCE: the floor is NOT a clean linear drift. "
            "Windowed slopes run 3.7e-4 -> 9.3e-4 -> 1.37e-3 "
            "across t = 0..10 (late/early ~3.7x): an ACCELERATING "
            "self-drift. The 1.3e-3 rad fit residual is ~14% of "
            "the total phase advance, i.e. systematic curvature, "
            "not noise. The ladder null survives because every "
            "rung shares the same curvature shape",
            "eps0.002 ran on the OLD base (om_start 0.199810), "
            "eps0.005/0.01 on the re-relaxed base: the ladder "
            "straddles two bases and only stays clean because "
            "each rung is compared to its own floor",
            "'perturbative' refers to B_meas only: at eps 0.01 the "
            "background configuration energy is 10117 vs the "
            "state's 6.7 (1500x); the STATE sits in an enormous "
            "but weakly-curved background"]}


# ================= CL7 =================
def cl7():
    say("== CL7: saturation ==")
    out = {}
    for name in ("eps0.03", "eps0.1", "env0.03", "env0.1"):
        r = jload(f"m5_21_9_larmor_{name}.json")
        Es = [x["E"] for x in r["rows"]]
        fin = [e for e in Es if isinstance(e, float)
               and np.isfinite(e)]
        nnf = len(Es) - len(fin)
        out[name] = {"E0": Es[0], "E_last_finite": fin[-1],
                     "n_nonfinite_rows": nnf}
        say(f"  {name}: E0 {Es[0]:.4g} nonfinite rows {nnf}")
    RES["CL7"] = {
        "verdict": "CONFIRMED",
        "claimed": {"eps0.03_E0": 8.1e5, "eps0.1_E0": 9.1e7,
                    "both_nan": True, "env0.03": "53.4 -> 47.5",
                    "env0.1_E0": 5770, "env0.1_nan": True},
        "measured": out,
        "notes": [
            "raw eps 0.03 / 0.1: E0 8.12e5 / 9.12e7 with 11 / 19 "
            "non-finite rows (NaN mid-run): claim holds",
            "env0.03 runs 53.40 -> 47.46 with zero non-finite "
            "rows; env0.1 E0 5770 with 5 non-finite rows: the "
            "enveloped window claim holds"]}


# ================= CL8 =================
def cl8():
    say("== CL8: the differential ==")
    f_fwd = fit_run("m5_21_9_larmor_env0.03.json")
    f_neg = fit_run("m5_21_9_larmor_env-0.03.json")
    f_rev = fit_run("m5_21_9_larmor_env0.03rev.json")
    f_flo = fit_run("m5_21_9_larmor_eps0v2.json")
    B = f_fwd["B_meas"]
    Op, Om, O0 = (f_fwd["Omega_full"], f_rev["Omega_full"],
                  f_flo["Omega_full"])
    sym = 0.5 * (Op + Om) - O0
    asym = 0.5 * (Op - Om)
    gam = asym / B
    say(f"  Omega +B {Op:.4e}  -B {Om:.4e}  0 {O0:.4e}")
    say(f"  sym {sym:+.3e} asym {asym:+.3e} gamma {gam:.4f}")
    # windowed refits: does the antisymmetric part live early?
    asym_early = 0.5 * (f_fwd["Omega_t<=5"] - f_rev["Omega_t<=5"])
    asym_late = 0.5 * (f_fwd["Omega_t>=5"] - f_rev["Omega_t>=5"])
    say(f"  asym early(t<=5) {asym_early:+.3e} "
        f"late(t>=5) {asym_late:+.3e}")
    # row-level replica agreement (the minus-eps control)
    ra = jload("m5_21_9_larmor_env0.03.json")["rows"]
    rb = jload("m5_21_9_larmor_env-0.03.json")["rows"]
    rep_dphi = max(abs(x["phi_xy"] - y["phi_xy"])
                   for x, y in zip(ra, rb))
    rep_dE = max(abs(x["E"] - y["E"]) for x, y in zip(ra, rb))
    say(f"  replica rows: max|dphi| {rep_dphi:.2e} "
        f"max|dE| {rep_dE:.2e}")
    f_fwd3 = fit_run("m5_21_9_larmor_env0.03.json", windows=W3)
    f_rev3 = fit_run("m5_21_9_larmor_env0.03rev.json", windows=W3)
    asym_w = {k: 0.5 * (f_fwd3[k] - f_rev3[k]) for k in W3K}
    say("  asym windows:", {k: f"{v:+.2e}"
                            for k, v in asym_w.items()})
    # total-energy (E + KE) trajectories of the pair
    etot = {}
    for lab, nm in (("fwd", "env0.03"), ("rev", "env0.03rev")):
        rr = jload(f"m5_21_9_larmor_{nm}.json")["rows"]
        etot[lab] = {"first": rr[0]["E"] + rr[0]["KE"],
                     "last": rr[-1]["E"] + rr[-1]["KE"],
                     "KE_last": rr[-1]["KE"]}
    say("  E+KE fwd", etot["fwd"], " rev", etot["rev"])
    # background F blocks with own 4x4 arithmetic
    eps = 0.03
    Tx = np.zeros((4, 4)); Tx[0, 1] = Tx[1, 0] = 1.0
    Ty = np.zeros((4, 4)); Ty[0, 2] = Ty[2, 0] = 1.0
    Axf, Ayf = eps * Ty, -eps * Tx      # fwd: d_x MB, d_y MB
    Axr, Ayr = eps * Ty, eps * Tx       # rev
    Ff = Axf @ ETA @ Ayf - Ayf @ ETA @ Axf
    Fr = Axr @ ETA @ Ayr - Ayr @ ETA @ Axr
    flip_exact = bool(np.allclose(Fr, -Ff, atol=0))
    Bnorm = float(np.sqrt(np.sum(Ff * Ff)))
    zmask = np.zeros((4, 4), dtype=bool)
    zmask[1, 2] = zmask[2, 1] = True
    zblock = bool(np.all(np.abs(Ff[~zmask]) < 1e-15)
                  and abs(Ff[1, 2]) > 0)
    say(f"  F_rev = -F_fwd exact: {flip_exact}; "
        f"|F| {Bnorm:.6e} (= sqrt2 eps^2 {np.sqrt(2)*eps**2:.6e})")
    # static cross-coupling: E0 of the two setups, own recompute
    Mf, _, cfg, _ = setup_own(0.03, rev=False)
    Mr, _, _, _ = setup_own(0.03, rev=True)
    Ef = float(INS4.e_total(Mf, cfg))
    Er = float(INS4.e_total(Mr, cfg))
    # field self-energy on the vacuum (symmetry check)
    Mv = np.broadcast_to(INS4.vac4(cfg),
                         Mf.shape).copy()
    env = INS4.envelope(cfg)[..., None, None]
    Evf = float(INS4.e_total(
        INS4.sym4(Mv + my_background(cfg, 0.03) * env), cfg))
    Evr = float(INS4.e_total(
        INS4.sym4(Mv + my_background(cfg, 0.03, rev=True) * env),
        cfg))
    say(f"  E0 state+field: fwd {Ef:.4f} rev {Er:.4f} "
        f"(split {Ef - Er:+.4f})")
    say(f"  field-on-vacuum self-energy: fwd {Evf:.4f} "
        f"rev {Evr:.4f}")
    # t = 0 torque split (force is velocity-independent)
    tq = {}
    for lab, rv in (("fwd", False), ("rev", True), ("none", None)):
        if rv is None:
            M, Mt, cfg2, freew = setup_own(0.0)
        else:
            M, Mt, cfg2, freew = setup_own(0.03, rev=rv)
        F0 = -INS4.grad(M, cfg2) / cfg2["h"] ** 3 * freew
        fl = my_rot_flows(M)
        J0 = {k: float(np.sum(Mt * A)) for k, A in fl.items()}
        tau = {k: float(np.sum(F0 * A)) for k, A in fl.items()}
        jsq = J0["x"] ** 2 + J0["y"] ** 2
        tq[lab] = {"phidot_t0": (J0["x"] * tau["y"]
                                 - J0["y"] * tau["x"]) / jsq,
                   "tau": tau, "Jxy": float(np.sqrt(jsq))}
        say(f"  t=0 torque {lab}: phidot {tq[lab]['phidot_t0']:+.4e}")
    RES["CL8"] = {
        "verdict": "NUANCE",
        "claimed": {"Omega_pB": 8.249e-4, "Omega_mB": 6.525e-4,
                    "replica_repro": 1e-7, "sym": -1.7e-4,
                    "asym": 0.86e-4, "gamma": 0.068,
                    "grade": "CANDIDATE (self-declared)"},
        "measured": {
            "Omega_pB": Op, "Omega_negeps": f_neg["Omega_full"],
            "replica_gap": abs(Op - f_neg["Omega_full"]),
            "replica_rowlevel_max_dphi": float(rep_dphi),
            "replica_rowlevel_max_dE": float(rep_dE),
            "Omega_mB": Om, "Omega_0": O0,
            "sym": sym, "asym": asym, "gamma_candidate": gam,
            "asym_early_t<=5": float(asym_early),
            "asym_late_t>=5": float(asym_late),
            "asym_windows": {k: float(v)
                             for k, v in asym_w.items()},
            "pair_Etot": etot,
            "F_rev_is_minus_F_fwd": flip_exact,
            "F_norm_eq_sqrt2_eps2": Bnorm,
            "F_z_block_only": zblock,
            "E0_fwd": Ef, "E0_rev": Er,
            "E0_texture_split": Ef - Er,
            "field_self_energy_fwd_rev": [Evf, Evr],
            "E_end_fwd": f_fwd["E_last"], "E_end_rev": f_rev["E_last"],
            "jxy_last_fwd": f_fwd["jxy_last"],
            "jxy_last_rev": f_rev["jxy_last"],
            "t0_torque": tq},
        "notes": [
            "all three Omegas, the replica gap, the decomposition "
            "and gamma reproduce with the auditor's own fits",
            "the reversed background DOES flip F_xy exactly (own "
            "4x4 arithmetic) and B_meas = sqrt2 eps^2 checks",
            "the minus-eps replica is BIT-IDENTICAL row by row "
            "(dphi = 0, dJ = 0, dE one ulp): stronger than the "
            "claimed 1e-7, but it is an exact eps -> -eps "
            "symmetry of the implementation, so it certifies "
            "determinism and evenness in eps, NOT F-mediation "
            "specifically (any even-in-eps coupling reproduces "
            "identically)",
            "REFUTING EVIDENCE 1 (timing): the antisymmetric "
            "split is ABSENT early: asym(t 0-3.5) ~ +0.5e-5, "
            "BELOW the floor resolvability bar 1.1e-5, then grows "
            "to ~+2.9e-4 late (t 6.5-10), exactly where the two "
            "trajectories diverge dynamically. A genuine Larmor "
            "term Omega_L = gamma B would be present from t = 0 "
            "at constant rate. The measured split has the time "
            "signature of the divergence, not of a field coupling",
            "REFUTING EVIDENCE 2 (texture): the reversed "
            "background is a DIFFERENT tensor texture (x T_y + "
            "y T_x, shear-type), not a symmetry image of the "
            "forward one (x T_y - y T_x, rotation-type): it "
            "couples to the state's static anisotropy already at "
            "t = 0 (setup energies split by ~1.09 while the field "
            "self-energy on the vacuum is sign-even), so "
            "Omega(-B) is not the same system with B flipped",
            "REFUTING EVIDENCE 3 (energetics): E + KE is nearly "
            "conserved AND nearly equal in both runs (~53.3 vs "
            "~52.2 at t = 10); the note's 'E_end 47.5 vs 10.9' is "
            "the POTENTIAL split only: at the last snapshot the "
            "reversed run holds ~41 units of KE (violent "
            "ringing). The pair is in different dynamical "
            "regimes late-time; the late-window fit rides that",
            "the omega = 0 discriminator (X1) came out NULL: "
            "neither texture pumps any J at omega = 0 (jxy stays "
            "< 2e-19 vs the 0.017 the omega-runs carry), so the "
            "simplest static-texture-torque alternative is "
            "EXCLUDED and that refutation route did NOT fire; "
            "the surviving contamination channel is the "
            "late-time relaxation asymmetry itself, seeded by "
            "the 1.09-unit setup-energy split",
            "the claim's own CANDIDATE grade and contamination "
            "flag are honest, but the windowed refit shows the "
            "resolvable antisymmetric signal is consistent with "
            "ZERO at matched trajectories (early-window asym "
            "+0.5e-5 to +1.0e-5, at or below the floor bar): "
            "gamma = 0.068 should be read as an artifact-scale "
            "upper bound until a matched-relaxation ladder "
            "resolves a constant-rate component"]}


def my_rot_flows(M):
    """own unit-Frobenius global-rotation flows."""
    out = {}
    for nm, (a, b) in (("x", (2, 3)), ("y", (3, 1)), ("z", (1, 2))):
        G = np.zeros((4, 4))
        G[a, b], G[b, a] = -1.0, 1.0
        A = G @ M - M @ G.swapaxes(-1, -2)
        out[nm] = A / max(np.sqrt(np.sum(A * A)), 1e-300)
    return out


# ================= CL9 =================
def cl9():
    say("== CL9: kin conventions ==")
    M, cfg = INS4.load_p1("p1_s-1")
    a0 = INS4.gen_catalog(cfg, M)["clock_local"]
    kin_probe = float(INS4.kin_of(M, a0, cfg))
    # own conjugation-orbit tangent (the M5.21.3-audit-adopted one)
    w = INS4.envelope(cfg)[..., None, None]
    lam, V = np.linalg.eigh(M[..., 1:4, 1:4])
    vh = V[..., :, 2]
    W = np.zeros(vh.shape[:-1] + (4, 4))
    n1, n2, n3 = vh[..., 0], vh[..., 1], vh[..., 2]
    W[..., 1, 2], W[..., 1, 3] = -n3, n2
    W[..., 2, 1], W[..., 2, 3] = n3, -n1
    W[..., 3, 1], W[..., 3, 2] = -n2, n1
    ac = w * (W @ M + M @ W.swapaxes(-1, -2))
    ac = ac / np.sqrt(np.sum(ac * ac))
    kin_conj = float(INS4.kin_of(M, ac, cfg))
    say(f"  kin probe-variant {kin_probe:.6f}  "
        f"conjugation {kin_conj:.6f}")
    RES["CL9"] = {
        "verdict": "NUANCE",
        "claimed": {"kin0": 0.297,
                    "framing": "0.119 = different convention era, "
                               "not a discrepancy"},
        "measured": {"kin_probe_variant_own": kin_probe,
                     "kin_conjugation_own": kin_conj,
                     "m5_21_3_p2_row": 0.296955480141499,
                     "m5_21_3_confirm_conj_s-1": 0.1205488175625523},
        "notes": [
            "mechanically CONFIRMED: gen_catalog + kin_of give "
            "0.29696 on the s = -1 endpoint, then and now (single "
            "commit, no code drift); own conjugation-tangent "
            "recompute gives 0.1205 matching the confirm row",
            "THE FRAMING UNDERSTATES: 0.119 was not merely an "
            "older era. The M5.21.3 audit adopted a STRUCTURAL "
            "CATCH: a0 = GM - MG^T is the antisymmetric probe "
            "variant, the physical conjugation-orbit tangent is "
            "GM + MG^T, and 'the quotable inertia is the "
            "conjugation value (0.119, not 0.295)'. This run "
            "builds J on the probe variant, i.e. it reverts the "
            "prior audit's adopted convention without a stated "
            "physics justification. Internal consistency (CL3/CL4) "
            "is unaffected; any ABSOLUTE J or gamma calibration "
            "(hbar/2 bridge) inherits a factor ~2.5 ambiguity"]}


# ================= X1: omega = 0 discriminator =================
def x1(steps=800, dt=0.005, snap=40):
    say("== X1: omega = 0 discriminator ==")
    out = {}
    for lab, eps, rev in (("none", 0.0, False),
                          ("fwd", 0.03, False),
                          ("rev", 0.03, True)):
        M, Mt, cfg, freew = setup_own(eps, rev=rev, om=0.0)
        gam = gamma_ramp(cfg)
        rows = []
        t0 = time.time()
        for k in range(steps + 1):
            if k % snap == 0:
                fl = my_rot_flows(M)
                j = {kk: float(np.sum(Mt * A))
                     for kk, A in fl.items()}
                rows.append({"t": k * dt, "J": j,
                             "E": etot_own(M, Mt, cfg)})
            if k < steps:
                M, Mt = leap_own(M, Mt, cfg, freew, gam, dt)
        t = np.array([r["t"] for r in rows])
        jx = np.array([r["J"]["x"] for r in rows])
        jy = np.array([r["J"]["y"] for r in rows])
        jxy = np.sqrt(jx ** 2 + jy ** 2)
        # azimuth once J is established (skip the birth transient)
        est = jxy > 0.1 * max(jxy.max(), 1e-300)
        sl = float("nan")
        if est.sum() >= 4:
            phi = np.unwrap(np.arctan2(jy[est], jx[est]))
            sl, _ = slope(t[est], phi)
        out[lab] = {"jxy_max": float(jxy.max()),
                    "jxy_end": float(jxy[-1]),
                    "jz_end": rows[-1]["J"]["z"],
                    "phidot_est": sl,
                    "E0": rows[0]["E"], "E_end": rows[-1]["E"],
                    "wall_s": time.time() - t0}
        say(f"  {lab}: jxy_end {out[lab]['jxy_end']:.3e} "
            f"phidot {sl:+.3e} E {rows[0]['E']:.2f} -> "
            f"{rows[-1]['E']:.2f} [{out[lab]['wall_s']:.0f}s]")
    # scale references
    ref_jxy = 0.017     # the J carried by the omega-run state
    split = out["fwd"]["phidot_est"] - out["rev"]["phidot_est"]
    RES["X1"] = {
        "purpose": "reversed field on the omega = 0 state: a "
                   "fwd/rev azimuth split here needs no J, so it "
                   "is a non-Larmor channel",
        "measured": out,
        "fwd_rev_phidot_split": float(split)
        if np.isfinite(split) else None,
        "jxy_scale_of_omega_runs": ref_jxy,
        "claimed_asym_component": 0.86e-4,
        "finding": X1_FINDING}


X1_FINDING = (
    "NULL: all three omega = 0 runs (no field, fwd, rev) keep "
    "jxy at machine zero (< 2e-19), so no J is texture-pumped "
    "and no azimuth split exists without the isorotation; the "
    "phidot_est values on those zero vectors are noise. The "
    "static-texture alternative to the Larmor reading is "
    "excluded; the antisymmetric split requires the J-carrying "
    "state and lives in the late-time trajectory divergence "
    "(see CL8 windowed evidence)")


def main():
    np.seterr(over="ignore", invalid="ignore")
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    old = None
    if mode == "reuse":
        # recompute every cheap section fresh; carry the long
        # dynamics blocks (CL5 reproduction, X1) over from the
        # existing JSON of the full run
        with open(os.path.join(DATA, "m5_21_9_audit.json")) as f:
            old = json.load(f)
    cl1()
    cl2()
    cl3()
    cl4()
    if mode == "all":
        cl5(run_dynamics=True)
    elif old is not None:
        cl5(run_dynamics=old["CL5"]["measured"]
            ["failure_mode_reproduction"])
    else:
        cl5(run_dynamics=False)
    cl6()
    cl7()
    cl8()
    cl9()
    if mode == "all":
        x1()
    elif old is not None and "X1" in old:
        RES["X1"] = old["X1"]
        RES["X1"]["finding"] = X1_FINDING
    with open(os.path.join(DATA, "m5_21_9_audit.json"), "w") as f:
        json.dump(RES, f, indent=1)
    say("saved ../data/m5_21_9_audit.json")


if __name__ == "__main__":
    main()
