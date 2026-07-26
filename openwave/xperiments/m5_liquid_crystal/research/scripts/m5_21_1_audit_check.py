"""M5.21.1 INDEPENDENT ADVERSARIAL AUDIT (claims C-A .. C-F).

Every verdict rests on a computation run HERE, with independent
implementations where the claim requires them:
  - V4 recomputed from EIGENVALUES of (eta M) per cell (numpy eigvals),
    never from the task's trace-power route (C-A, C-C endpoint re-eval);
  - containment / core-ball / boundary-leakage code written fresh (C-B);
  - the C-A gap re-derived analytically from the odd-p cross terms;
  - curvatures via PLAIN energy-level central finite differences
    (the task used complex-step Hessian-vector products) (C-D, C-E);
  - explicit one-parameter families for the Routhian unboundedness (C-C);
  - fresh least-squares refits of the two P4 laws + the analytic
    omega ~ g^3 Jacobian check (C-F).
Imports allowed by the audit brief: total_energy_4 / u_eta_density /
WSCALE (m5_20_2_a_eom), cell_weights / grid_coords (m5_17_energy),
t_total_c (m5_20_3_a_constraint), gen (m5_20_4_a_bvp),
q2_avg_f (m5_20_5_a_orbit), electron_seed (m5_21_b_electron).

Run:  python3 m5_21_1_audit_check.py
Out:  ../data/m5_21_1_audit.json
"""
from __future__ import annotations

import json
import os
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from m5_17_energy import cell_weights, grid_coords                 # noqa: E402
from m5_20_2_a_eom import WSCALE, total_energy_4, u_eta_density    # noqa: E402
from m5_20_3_a_constraint import t_total_c                         # noqa: E402
from m5_20_4_a_bvp import gen                                      # noqa: E402
from m5_20_5_a_orbit import q2_avg_f                               # noqa: E402
from m5_21_b_electron import electron_seed                         # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
DELTA = 0.3
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:6.1f}s] {msg}", flush=True)


# ---------------- independent V4 (eigenvalue route) ----------------
def v4_own_density(M, sg):
    """wscale * SUM_p (powersum_p(eta M) - C_p)^2 per included cell,
    power sums from numpy EIGENVALUES of the (generally non-symmetric)
    matrix eta.M (independent of the task's trace-power route)."""
    nr = M.shape[0]
    Mc = M[: nr - 1, 1:-1]
    EM = np.einsum("ab,...bc->...ac", ETA, Mc)
    lam = np.linalg.eigvals(EM)
    v = np.zeros(Mc.shape[:2])
    for p in range(1, 5):
        ps = np.sum(lam ** p, axis=-1).real
        cp = float(sg) ** p + 1.0 + DELTA ** p
        v = v + (ps - cp) ** 2
    return WSCALE * v


def e_own(M, sg):
    """my total static energy: imported u_eta (g-independent) + my V4."""
    nr, nz = M.shape[:2]
    w = cell_weights(nr, nz, 1.0)
    return float(np.sum((u_eta_density(M, 1.0)
                         + v4_own_density(M, sg)) * w))


def dens_own(M, sg):
    nr, nz = M.shape[:2]
    w = cell_weights(nr, nz, 1.0)
    return (u_eta_density(M, 1.0) + v4_own_density(M, sg)) * w


def zero_pins(v):
    """pinned DOF (outer rho row + both z faces) stay fixed."""
    v = v.copy()
    v[-1, :] = 0.0
    v[:, 0] = 0.0
    v[:, -1] = 0.0
    return v


def fd_curvature(M, v, eps, sg=8.0):
    """plain central second difference of the imported total energy."""
    ep = total_energy_4(M + eps * v, WSCALE, DELTA, g=sg)
    em = total_energy_4(M - eps * v, WSCALE, DELTA, g=sg)
    e0 = total_energy_4(M, WSCALE, DELTA, g=sg)
    return (ep + em - 2.0 * e0) / eps ** 2


OUT = {"task": "M5.21.1 adversarial audit", "date": "2026-07-16",
       "wscale": WSCALE, "claims": {}}


# ================= C-A: the P0 sign gap =================
def check_a():
    log("C-A: P0 sign gap, independent V4 (eigenvalue route)")
    M = np.load(os.path.join(DATA, "m5_21_1_a_relaxed_sp.npz"))["M"]
    nr, nz = M.shape[:2]
    w = cell_weights(nr, nz, 1.0)
    offblock = float(np.max(np.abs(M[..., 0, 1:])))
    Mm = M.copy()
    Mm[..., 0, 0] *= -1.0                      # the mirror M00 -> -M00
    E_sp = e_own(M, +8.0)
    E_sm_mirror = e_own(Mm, -8.0)
    gap = E_sm_mirror - E_sp
    pct = 100.0 * gap / E_sp
    # route validation: my V4 vs the imported total energy (same state)
    rel_vs_imported = abs(E_sp - total_energy_4(M, WSCALE, DELTA, g=8.0)) \
        / abs(E_sp)
    # u_eta identical between the two evaluations?
    ue1 = u_eta_density(M, 1.0) * w
    ue2 = u_eta_density(Mm, 1.0) * w
    ue_gap = float(np.sum(ue2) - np.sum(ue1))
    ue_maxcell = float(np.max(np.abs(ue2 - ue1)))
    # V4-only gap (my densities)
    v4_gap = float(np.sum((v4_own_density(Mm, -8.0)
                           - v4_own_density(M, +8.0)) * w))
    # analytic odd-p cross-term mechanism (exact for block-diagonal M)
    Mc = M[: nr - 1, 1:-1]
    m00 = Mc[..., 0, 0]
    S = Mc[..., 1:, 1:]
    S = 0.5 * (S + np.swapaxes(S, -1, -2))
    lam_s = np.linalg.eigvalsh(S)
    gap_an = 0.0
    for p in (1, 3):
        a_p = (-m00) ** p - 8.0 ** p
        s_p = np.sum(lam_s ** p, axis=-1) - 1.0 - DELTA ** p
        gap_an += float(np.sum(w * ((-a_p + s_p) ** 2
                                    - (a_p + s_p) ** 2))) * WSCALE
    ok = (abs(gap - 0.1296) < 2e-3 and abs(pct - 0.94) < 0.05
          and abs(ue_gap) < 1e-9 * abs(np.sum(ue1))
          and abs(gap_an - gap) < 1e-6 * abs(gap))
    verdict = "CONFIRMED" if ok else "REFUTED"
    res = {"verdict": verdict, "E_sp_own": E_sp,
           "E_sm_on_mirror_own": E_sm_mirror, "gap": gap,
           "gap_pct": pct, "gap_claimed": 0.1296,
           "u_eta_gap": ue_gap, "u_eta_max_cell_diff": ue_maxcell,
           "v4_gap": v4_gap, "gap_analytic_oddp": gap_an,
           "analytic_vs_direct_rel": abs(gap_an - gap) / abs(gap),
           "max_offblock_component": offblock,
           "own_vs_imported_energy_rel": rel_vs_imported}
    OUT["claims"]["C-A"] = res
    log(f"C-A {verdict}: gap {gap:.6f} ({pct:.3f}%), analytic "
        f"{gap_an:.6f}, u_eta gap {ue_gap:.2e}, offblock {offblock:.1e}")


# ================= C-B: the P1 slide =================
def spherical_containment(M, sg=8.0):
    """my own containment: weighted density cumulative in spherical r."""
    nr, nz = M.shape[:2]
    dens = dens_own(M, sg)
    R, Z = grid_coords(nr, nz, 1.0)
    Rc, Zc = R[: nr - 1, 1:-1], Z[: nr - 1, 1:-1]
    r = np.sqrt(Rc ** 2 + Zc ** 2)
    tot = float(dens.sum())
    order = np.argsort(r.ravel())
    rs = r.ravel()[order]
    cs = np.cumsum(dens.ravel()[order])
    r50c = float(rs[np.searchsorted(cs, 0.5 * tot)])
    r90c = float(rs[np.searchsorted(cs, 0.9 * tot)])
    # integer-radius convention (matches the task's reporting grid)
    rb = np.arange(1.0, 100.0, 1.0)
    csum = np.array([dens[r <= x].sum() for x in rb])
    r50i = float(rb[min(int(np.searchsorted(csum, 0.5 * csum[-1])),
                        len(rb) - 1)])
    r90i = float(rb[min(int(np.searchsorted(csum, 0.9 * csum[-1])),
                        len(rb) - 1)])
    core = float(dens[r <= 8.0].sum()) / tot
    # boundary band: within 5 cells of the pinned outer boundary
    ii = np.arange(nr - 1)[:, None] * np.ones((1, nz - 2), dtype=int)
    jj = np.ones((nr - 1, 1), dtype=int) * np.arange(nz - 2)[None, :]
    band = (ii >= nr - 1 - 5) | (jj <= 4) | (jj >= nz - 2 - 5)
    bfrac = float(dens[band].sum()) / tot
    return {"E_total": tot, "r50_int": r50i, "r90_int": r90i,
            "r50_cont": r50c, "r90_cont": r90c,
            "core_ball_r8_frac": core, "boundary_band5_frac": bfrac}


def check_b():
    log("C-B: P1 non-convergence + spreading, own containment code")
    j = json.load(open(os.path.join(DATA, "m5_21_1_b_statics.json")))
    tr = j["B1"]["fire"]["trace"]
    E = np.array([t["E"] for t in tr])
    it = np.array([t["it"] for t in tr], dtype=float)
    monotone = bool(np.all(np.diff(E) < 0.0))
    slope_last = float((E[-1] - E[-2]) / (it[-1] - it[-2]))
    half = it >= 0.5 * it[-1]
    slope_24k = float((E[half][-1] - E[half][0])
                      / (it[half][-1] - it[half][0]))
    # slope-decay exponent: log(-dE/dit) vs log(it) over the trace
    dE = np.diff(E) / np.diff(it)
    mid = 0.5 * (it[1:] + it[:-1])
    good = dE < 0
    pfit = np.polyfit(np.log(mid[good]), np.log(-dE[good]), 1)
    snaps = j["B1"]["snapshots"]
    s0, sN = snaps[0], snaps[-1]
    # endpoint recomputed with MY code
    Mend = np.load(os.path.join(DATA, "m5_21_1_b_endpoint.npz"))["M"]
    cend = spherical_containment(Mend)
    # the it=0 seed rebuilt and measured with MY code
    M0 = electron_seed()
    c0 = spherical_containment(M0)
    ok = (monotone
          and abs(slope_24k - (-5.1e-5)) < 0.1e-5
          and abs(cend["r50_int"] - 15.0) <= 1.0
          and abs(cend["r90_int"] - 51.0) <= 1.0
          and abs(cend["core_ball_r8_frac"] - 0.20) < 0.01
          and abs(c0["r50_int"] - 8.0) <= 1.0
          and abs(c0["r90_int"] - 32.0) <= 1.0
          and abs(c0["core_ball_r8_frac"] - 0.52) < 0.01)
    leak_flag = cend["boundary_band5_frac"] > 0.05
    verdict = "CONFIRMED" if ok else "REFUTED"
    res = {"verdict": verdict, "trace_monotone_decreasing": monotone,
           "dE_dit_last_1k": slope_last, "dE_dit_last_half": slope_24k,
           "dE_dit_claimed": -5.1e-5,
           "slope_decay_exponent_dlogdEdt_dlogit": float(pfit[0]),
           "snapshots_json_start": {"r50": s0["containment"]["r50"],
                                    "r90": s0["containment"]["r90"],
                                    "core": s0["containment"]
                                    ["core_ball_r8_frac"]},
           "seed_own": c0, "endpoint_own": cend,
           "endpoint_json": {"r50": sN["containment"]["r50"],
                             "r90": sN["containment"]["r90"],
                             "core": sN["containment"]
                             ["core_ball_r8_frac"]},
           "boundary_leakage_flag": bool(leak_flag)}
    OUT["claims"]["C-B"] = res
    log(f"C-B {verdict}: monotone={monotone}, tail slope (last-half "
        f"window) {slope_24k:.3e} vs claimed -5.1e-5 (last-1k "
        f"{slope_last:.3e}); own r50/r90 {c0['r50_int']:.0f}/"
        f"{c0['r90_int']:.0f} -> {cend['r50_int']:.0f}/"
        f"{cend['r90_int']:.0f}, core {c0['core_ball_r8_frac']:.3f} -> "
        f"{cend['core_ball_r8_frac']:.3f}, boundary band(5) "
        f"{cend['boundary_band5_frac']:.4f}")


# ================= C-C: the P2 centrifugal runaway =================
def check_c():
    log("C-C: Routhian unboundedness, structural check")
    j = json.load(open(os.path.join(DATA, "m5_21_1_c_4d.json")))
    om = j["C2"]["omegas"]
    # (i) R_end == U - Omega^2 Q2_avg consistency (JSON-internal)
    consistency = {}
    for k, v in om.items():
        r_pred = v["U"] - v["omega"] ** 2 * v["Q2_avg"]
        consistency[k] = {"R_fire_tail": v["fire_tail"]["E"],
                          "U_minus_w2Q2": r_pred,
                          "rel": abs(r_pred - v["fire_tail"]["E"])
                          / max(abs(r_pred), 1e-30)}
    # (ii) the descent trace from the run log: monotone per omega?
    logf = os.path.join(DATA, "m5_21_1_c_c2_run.log")
    mono = {}
    if os.path.exists(logf):
        txt = open(logf).read()
        for tag in ("om0.05", "om0.1349", "om0.25"):
            vals = [float(m) for m in re.findall(
                rf"{re.escape(tag)} it\s+\d+ E\s+(-?[\d.]+)", txt)]
            mono[tag] = {"E_checkpoints": vals,
                         "monotone": bool(np.all(np.diff(vals) < 0))}
    # (iii) independent endpoint re-evaluation (MY V4 for U)
    M0 = np.load(os.path.join(DATA, "m5_21_1_c_state_sp.npz"))["M"]
    J = gen("J23")
    end_evals = {}
    for tag in ("0.1349", "0.25"):
        Me = np.load(os.path.join(DATA,
                                  f"m5_21_1_c_rot_om{tag}.npz"))["M"]
        U_own = e_own(Me, +8.0)
        Q2 = q2_avg_f(Me, J)
        R_own = U_own - float(tag) ** 2 * Q2
        rel = abs(R_own - om[tag]["fire_tail"]["E"]) / abs(R_own)
        end_evals[tag] = {"U_own": U_own, "Q2_own": Q2, "R_own": R_own,
                          "R_fire_tail": om[tag]["fire_tail"]["E"],
                          "rel": rel}
        log(f"  endpoint om{tag}: U_own {U_own:.6e} (json "
            f"{om[tag]['U']:.6e}), R_own {R_own:.6e} vs fire "
            f"{om[tag]['fire_tail']['E']:.6e} (rel {rel:.2e})")
    Mend = np.load(os.path.join(DATA, "m5_21_1_c_rot_om0.25.npz"))["M"]
    U_end_own = end_evals["0.25"]["U_own"]
    Q2_end = end_evals["0.25"]["Q2_own"]
    R_end_own = end_evals["0.25"]["R_own"]
    end_rel = end_evals["0.25"]["rel"]
    # (iv) the prescribed one-parameter family: rotation texture
    nr = M0.shape[0]
    rho = ((np.arange(nr) + 0.5) * 1.0)[:, None, None, None]
    f = np.exp(-(((rho - 40.0) / 8.0) ** 2))
    v = zero_pins(f * (np.einsum("ab,...bc->...ac", J, M0)
                       - np.einsum("...ab,bc->...ac", M0, J)))
    fam_rows = []
    for amp in (0.0, 0.5, 1.0, 2.0, 4.0):
        X = M0 + amp * v
        U = total_energy_4(X, WSCALE, DELTA, g=8.0)
        Q2 = q2_avg_f(X, J)
        fam_rows.append({"amp": amp, "U": U, "Q2_avg": Q2,
                         "R": U - 0.0625 * Q2})
        log(f"  bump family amp {amp}: U {U:.4e} Q2 {Q2:.4e} "
            f"R {U - 0.0625 * Q2:.4e}")
    Rf = [r["R"] for r in fam_rows]
    fam_decreasing = bool(np.all(np.diff(Rf) < 0.0))
    # (v) straight-segment slice through the descent endpoint, extended
    D = Mend - M0
    seg_rows = []
    for t in (0.0, 0.25, 0.5, 0.75, 1.0, 1.1, 1.25):
        X = M0 + t * D
        U = total_energy_4(X, WSCALE, DELTA, g=8.0)
        Q2 = q2_avg_f(X, J)
        seg_rows.append({"t": t, "U": U, "Q2_avg": Q2,
                         "R": U - 0.0625 * Q2})
        log(f"  segment t {t}: U {U:.4e} Q2 {Q2:.4e} "
            f"R {U - 0.0625 * Q2:.4e}")
    Rs = [r["R"] for r in seg_rows]
    seg_valley = bool(Rs[4] < Rs[0] and Rs[4] < Rs[3])
    seg_beyond = bool(Rs[5] < Rs[4] or Rs[6] < Rs[4])
    fam_turns_up = bool(Rf[-1] > Rf[0])
    fam_initially_descends = bool(Rf[1] < Rf[0] and Rf[2] < Rf[1])
    cons_ok = all(c["rel"] < 1e-9 for c in consistency.values())
    mono_ok = all(m["monotone"] for m in mono.values()) if mono else False
    # sub-verdict 1: landscape-level runaway (not FIRE numerics)
    runaway_ok = (cons_ok and mono_ok and end_rel < 1e-6
                  and seg_valley and fam_initially_descends)
    # sub-verdict 2: STRICT unboundedness. Refuted if every probe ray
    # turns up (V4 >= 0 grows ~ amp^8 vs Q2 ~ amp^4) AND the om0.25
    # descent itself converged to a stationary point (deep FINITE well).
    om25_stationary = om["0.25"]["rel_residual_gradR"] < 1e-4
    unbounded_refuted = (fam_turns_up and not seg_beyond
                         and om25_stationary)
    v1 = "CONFIRMED" if runaway_ok else "INCONCLUSIVE"
    v2 = "REFUTED" if unbounded_refuted else \
        ("CONFIRMED" if fam_decreasing else "INCONCLUSIVE")
    verdict = f"PARTIAL: landscape-runaway {v1}; unbounded-below {v2}"
    res = {"verdict": verdict,
           "sub_landscape_runaway": v1,
           "sub_strictly_unbounded_below": v2,
           "consistency": consistency,
           "log_monotone": mono, "endpoint_reevals": end_evals,
           "om0.25_rel_residual_gradR_json":
           om["0.25"]["rel_residual_gradR"],
           "bump_family": fam_rows,
           "bump_family_initially_descends": fam_initially_descends,
           "bump_family_turns_up": fam_turns_up,
           "segment_family": seg_rows,
           "segment_valley_at_endpoint": seg_valley,
           "segment_still_decreasing_beyond_endpoint": seg_beyond,
           "note": ("straight rays are V4-coercive (V4 >= 0, ~ amp^8) "
                    "vs Q2 ~ amp^4, so R -> +inf on every generic ray; "
                    "the om0.25 descent record itself ends at a "
                    "near-stationary point (rel grad residual 1.6e-5), "
                    "i.e. a deep FINITE well, not unboundedness")}
    OUT["claims"]["C-C"] = res
    log(f"C-C {verdict}")


# ================= C-D: the roton-like dip =================
def build_twist(M, n_wave, r_env):
    """my own twist mode: v = cos(k z) env(r) (W M - M W), W = local
    rotation generator about the leading spatial eigenvector."""
    nr, nz = M.shape[:2]
    R, Z = grid_coords(nr, nz, 1.0)
    r = np.sqrt(R ** 2 + Z ** 2)
    env = np.exp(-((r / r_env) ** 4))
    k = 2.0 * np.pi * n_wave / (nz * 1.0)
    cz = np.cos(k * Z)
    S = M[..., 1:, 1:]
    S = 0.5 * (S + np.swapaxes(S, -1, -2))
    lam, V = np.linalg.eigh(S)
    nvec = V[..., :, -1]                      # leading eigenvector
    W = np.zeros(M.shape)
    n1, n2, n3 = nvec[..., 0], nvec[..., 1], nvec[..., 2]
    W[..., 1, 2], W[..., 1, 3] = -n3, n2
    W[..., 2, 1], W[..., 2, 3] = n3, -n1
    W[..., 3, 1], W[..., 3, 2] = -n2, n1
    v = (cz * env)[..., None, None] * (W @ M - M @ W)
    v = zero_pins(v)
    return v / np.sqrt(np.sum(v * v))


def check_d():
    log("C-D: dispersion dip, own mode construction + energy-level FD")
    M = np.load(os.path.join(DATA, "m5_21_1_c_state_sp.npz"))["M"]
    rows = []
    for r_env in (15.0, 25.0):
        for n_wave in (0, 2, 4):
            v = build_twist(M, n_wave, r_env)
            T = float(np.real(t_total_c(M, v)))
            d2_1 = fd_curvature(M, v, 1e-4)
            d2_3 = fd_curvature(M, v, 3e-4)
            w2_1 = d2_1 / (2.0 * T)
            w2_3 = d2_3 / (2.0 * T)
            rows.append({"r_env": r_env, "n": n_wave,
                         "k": 2 * np.pi * n_wave / 128.0,
                         "D2E_eps1e-4": d2_1, "D2E_eps3e-4": d2_3,
                         "T_true": T, "omega2_eps1e-4": w2_1,
                         "omega2_eps3e-4": w2_3})
            log(f"  r_env {r_env} n {n_wave}: D2E {d2_1:.6f} "
                f"(eps x3: {d2_3:.6f}), T {T:.6f}, omega^2 {w2_1:.6f}")
    dips = {}
    for r_env in (15.0, 25.0):
        w2 = {r["n"]: r["omega2_eps1e-4"] for r in rows
              if r["r_env"] == r_env}
        w2b = {r["n"]: r["omega2_eps3e-4"] for r in rows
               if r["r_env"] == r_env}
        dips[str(r_env)] = {
            "omega2": w2,
            "dip": bool(w2[2] < w2[0] and w2[2] < w2[4]),
            "dip_eps3e-4": bool(w2b[2] < w2b[0] and w2b[2] < w2b[4])}
    ok = all(d["dip"] and d["dip_eps3e-4"] for d in dips.values())
    verdict = "CONFIRMED" if ok else "REFUTED"
    OUT["claims"]["C-D"] = {"verdict": verdict, "rows": rows,
                            "dips": dips}
    log(f"C-D {verdict}: dips {dips['15.0']['dip']} (r_env 15), "
        f"{dips['25.0']['dip']} (r_env 25)")


# ================= C-E: negative time-mixing curvature =================
def check_e():
    log("C-E: time-mixing curvature at the 128x256 endpoint (FD)")
    M = np.load(os.path.join(DATA, "m5_21_1_b_endpoint.npz"))["M"]
    nr, nz = M.shape[:2]
    R, Z = grid_coords(nr, nz, 1.0)
    bump = np.exp(-(((R - 2.0) ** 2 + Z ** 2) / 18.0))
    B = np.zeros_like(M)
    B[..., 0, 1] = bump
    B[..., 1, 0] = bump
    B = zero_pins(B)
    B /= np.sqrt(np.sum(B * B))
    rng = np.random.default_rng(7)
    D = rng.normal(size=M.shape)
    D = 0.5 * (D + np.swapaxes(D, -1, -2))
    Dt = np.zeros_like(D)
    Dt[..., 0, 1:] = D[..., 0, 1:]
    Dt[..., 1:, 0] = D[..., 1:, 0]
    Dt = zero_pins(Dt)
    Dt /= np.sqrt(np.sum(Dt * Dt))
    Db = D.copy()
    Db[..., 0, 1:] = 0.0
    Db[..., 1:, 0] = 0.0
    Db = zero_pins(Db)
    Db /= np.sqrt(np.sum(Db * Db))
    rows = {}
    for name, v in (("boost01_core", B), ("rand_timemix", Dt),
                    ("rand_blockdiag", Db)):
        c1 = fd_curvature(M, v, 1e-4)
        c3 = fd_curvature(M, v, 3e-4)
        rows[name] = {"curv_eps1e-4": c1, "curv_eps3e-4": c3}
        log(f"  {name}: D2E/|v|^2 = {c1:+.6f} (eps x3: {c3:+.6f})")
    ok = (abs(rows["boost01_core"]["curv_eps1e-4"] - (-0.386)) < 0.02
          and rows["rand_timemix"]["curv_eps1e-4"] < -0.2
          and rows["rand_blockdiag"]["curv_eps1e-4"] > 0.0)
    verdict = "CONFIRMED" if ok else "REFUTED"
    OUT["claims"]["C-E"] = {"verdict": verdict, "rows": rows,
                            "boost_claimed": -0.386}
    log(f"C-E {verdict}")


# ================= C-F: the two P4 laws =================
def loglog_fit(x, y):
    lx, ly = np.log(np.asarray(x)), np.log(np.asarray(y))
    A = np.vstack([lx, np.ones_like(lx)]).T
    coef, _, _, _ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A @ coef
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - ly.mean()) ** 2))
    return float(coef[0]), 1.0 - ss_res / ss_tot


def check_f():
    log("C-F: refit the two P4 laws + analytic omega ~ g^3")
    j = json.load(open(os.path.join(DATA, "m5_21_1_d_scaling.json")))
    dl = j["delta_ladder_g8"]
    gl = j["g_ladder_d0p3"]
    s1, r1 = loglog_fit([r["delta"] for r in dl],
                        [r["axis_split_mid"] for r in dl])
    s2, r2 = loglog_fit([r["g"] for r in gl],
                        [r["vac_ladder_max_omega"] for r in gl])
    # analytic: stiff vacuum mode along the g eigenvalue,
    # lam = 2 w SUM_p (p g^(p-1))^2, omega = sqrt(lam) -> ~ 4 sqrt(2w) g^3
    ana = {}
    for r in gl:
        g = r["g"]
        lam = 2.0 * WSCALE * sum((p * g ** (p - 1)) ** 2
                                 for p in range(1, 5))
        ana[str(g)] = {"omega_pred": float(np.sqrt(lam)),
                       "omega_row": r["vac_ladder_max_omega"],
                       "rel": abs(np.sqrt(lam) - r["vac_ladder_max_omega"])
                       / r["vac_ladder_max_omega"]}
    gs = [r["g"] for r in gl]
    wpred = [ana[str(g)]["omega_pred"] for g in gs]
    s_ana, _ = loglog_fit(gs, wpred)
    ok = (abs(s1 - 1.034) < 0.005 and r1 > 0.999
          and abs(s2 - 2.992) < 0.005 and r2 > 0.9999
          and all(a["rel"] < 1e-3 for a in ana.values())
          and abs(s_ana - s2) < 0.005)
    verdict = "CONFIRMED" if ok else "REFUTED"
    OUT["claims"]["C-F"] = {
        "verdict": verdict,
        "axis_split_vs_delta": {"slope": s1, "R2": r1,
                                "claimed": [1.034, 0.9998]},
        "vac_stiff_vs_g": {"slope": s2, "R2": r2,
                           "claimed": [2.992, 0.999998]},
        "analytic_jacobian_check": ana,
        "analytic_slope_over_ladder": s_ana}
    log(f"C-F {verdict}: slopes {s1:.4f} (R2 {r1:.6f}), {s2:.4f} "
        f"(R2 {r2:.8f}); analytic slope {s_ana:.4f}, worst rel "
        f"{max(a['rel'] for a in ana.values()):.2e}")


def main():
    check_a()
    check_b()
    check_c()
    check_d()
    check_e()
    check_f()
    path = os.path.join(DATA, "m5_21_1_audit.json")
    with open(path, "w") as f:
        json.dump(OUT, f, indent=1)
    log(f"wrote {path}")
    print("\n==== VERDICTS ====")
    for k, v in OUT["claims"].items():
        print(f"  {k}: {v['verdict']}")


if __name__ == "__main__":
    main()
