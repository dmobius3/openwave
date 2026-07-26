"""M5.21.1 phases P2 + P3: the 4D extension by MINIMIZATION (rotation-
sector orbit class, boost quarantined) + the Klein-Gordon twist limit.

HIS CHAIN (Duda 2026-07-15): "minimizing energy in 3D ... getting
static solution, then extending to 4D - where energy minimization
should lead to angular momentum + gravitational mass"; "for
perpendicular low energy twists (as quantum phase) you should get
Klein-Gordon-like equation".

OPERATIONALIZATION (no raw time integration anywhere: the M5.20.3
ill-posed-IVP lesson; dynamics enters as minimization, orbit classes,
and spectra only):
    P2  RIGID ROTATION-SECTOR ORBIT CLASS. M(t) = e^{Omega t G} M0
        e^{Omega t G^T}, G a (possibly boost-conjugated) so(1,3)
        generator. True kinetic on the orbit (the audited quartic L):
        T = Omega^2 Q2_avg(M, G), Q2_avg = <T_true(M, D_G M)>_phi
        (the M5.20.5 corrected phi-average, nphi = 5 band-limit
        EXACT, re-gated here on the hedgehog). Rotating-frame energy
        (Routhian) R(M; Omega) = U(M) - Omega^2 Q2_avg(M); its
        extremals gradU = Omega^2 gradQ2 are the rigid rotating
        equilibria. Reads: the a2x ALIGNMENT diagnostic FIRST
        (cos(g_kin, g_stat): whether ANY Omega can balance), then
        FIRE descent on R at fixed Omega; J = 2 Omega Q2_avg,
        mass split U vs T = Omega^2 Q2_avg.
        Boost quarantine: primary reads use PURE spatial rotations
        (chi = 0); chi-conjugated families are reported as labeled
        boost-adjacent context (the M5.20.4/5 free-period objects).
    P3  KG TWIST LIMIT, spectrally. Perpendicular low-energy twist
        modes v_k = cos(k z) env(r) [W(n), M], W(n) = local rotation
        about the director (the m5_21_c clock generator). Harmonic
        frequency under the TRUE kinetic metric:
            omega_true^2(k) = D2E(v_k) / (2 T_true(M, v_k)),
        (D2E = <v, Hess_E v> by complex-step HV on the audited
        gradient), plus the canonical-metric read for continuity.
        KG-like verdict: omega^2 = m^2 + c^2 k^2 fit quality.

GRID. P2 solves on the 64 x 128 axisym grid (the M5.20.3/4/5 series
precedent; the q2 library's audited module constants ARE this grid),
with the hedgehog re-relaxed per sign on it (c0). Cross-grid
consistency: single-eval a2x reads on the 128 x 256 P1 endpoint.

GATES (pre-registered)
    GC0  hedgehog re-relax per sign on 64 x 128: descent + q intact
    GC1  instrument re-gates ON THE HEDGEHOG background:
         grad_q2_avg == complex step (<= 1e-8); nphi 5 == 16
         (<= 1e-10); J12-commutant q2_avg == q2_of (<= 1e-12)
    GC2  KG mode sanity: T_true(M, v_k) > 0 (the twist is in the
         POSITIVE kinetic sector: quarantine holds for these modes)

Run:  python m5_21_1_c_4d.py c0|c1|c2 [omegas]|c3|all
Out:  ../data/m5_21_1_c_4d.json, ../data/m5_21_1_c_state_s[pm].npz,
      ../plots/m5_21_1_c_film_basic.png / _film_thermal.png /
      _a2x.png / _kg.png
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ARGV = sys.argv[1:]
sys.argv = sys.argv[:1]

from m5_16_axisym import pin_mask                                   # noqa: E402
from m5_17_energy import cell_weights, grid_coords                  # noqa: E402
from m5_20_2_a_eom import (G_T, H, WSCALE, grad_static_4,           # noqa: E402
                           total_energy_4)
from m5_20_3_a_constraint import _cs_dir, e_static_c, t_total_c     # noqa: E402
from m5_20_4_a_bvp import conj_gen, d_g, gen, grad_q2, q2_of        # noqa: E402
from m5_20_5_a_orbit import gens_phi, q2_avg_f, q2_avg_n            # noqa: E402
from m5_21_a_snap import eig_fields, spectral_amplitude             # noqa: E402
from m5_21_c_clockrun import local_gen                              # noqa: E402
from m5_film import film_strip                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
DELTA = 0.3
NPHI = 5
NR2, NZ2 = 64, 128                 # the M5.20-series P2 grid
W42 = cell_weights(NR2, NZ2, H)[..., None, None]
RHO42 = ((np.arange(NR2 - 1) + 0.5) * H)[:, None, None, None]
PIN2 = pin_mask(NR2, NZ2)
FREE42 = (~PIN2)[..., None, None].astype(float)
SIGNS = {"sp": +1.0, "sm": -1.0}
FAMS = {"J23^K2": ("J23", "K2", 0.5), "J23^K3": ("J23", "K3", 0.3),
        "J12^K1": ("J12", "K1", 0.3)}
OMEGAS = (0.05, 0.1349, 0.25)      # rung 0.1349 = the gap-ladder line


# ================= 64 x 128 stack =================
def seed64(delta, sg, r_c=4.0):
    R, Z = grid_coords(NR2, NZ2, H)
    r = np.sqrt(R ** 2 + Z ** 2)
    rs = np.where(r < 1e-12, 1e-12, r)
    n = np.stack([R / rs, np.zeros_like(R), Z / rs], axis=-1)
    m = np.broadcast_to(np.array([0.0, 1.0, 0.0]), n.shape)
    S = (n[..., :, None] * n[..., None, :]
         + delta * m[..., :, None] * m[..., None, :])
    a = (1.0 + delta) / 3.0
    w = (1.0 - np.exp(-((r / r_c) ** 2)))[..., None, None]
    M = np.zeros((NR2, NZ2, 4, 4))
    M[..., 1:4, 1:4] = w * S + (1.0 - w) * (a * np.eye(3))
    M[..., 0, 0] = -sg
    return M


def e_stat_s(M, sg):
    return float(e_static_c(M, WSCALE, DELTA, g=sg).real)


def g_stat_s(M, sg):
    return grad_static_4(M, WSCALE, DELTA, g=sg, w4=W42, rho4=RHO42)


def q_meridional64(M, r_m=10.0, npts=121):
    R, Z = grid_coords(NR2, NZ2, H)
    lam, V = eig_fields(M)
    rr = np.where(np.sqrt(R ** 2 + Z ** 2) < 1e-12, 1e-12,
                  np.sqrt(R ** 2 + Z ** 2))
    n = V[..., :, 0]
    ref = np.stack([R / rr, np.zeros_like(R), Z / rr], axis=-1)
    s = np.sign(np.sum(n * ref, axis=-1, keepdims=True))
    s[s == 0.0] = 1.0
    n = n * s
    th = np.linspace(1e-3, np.pi - 1e-3, npts)
    i = np.clip(np.round(r_m * np.sin(th) / H - 0.5).astype(int),
                0, NR2 - 1)
    j = np.clip(np.round(r_m * np.cos(th) / H + NZ2 / 2
                         - 0.5).astype(int), 0, NZ2 - 1)
    alpha = np.arctan2(n[i, j, 0], n[i, j, 2])
    dd = np.diff(alpha)
    dd = (dd + np.pi / 2) % np.pi - np.pi / 2
    return float(np.sum(dd) / np.pi)


def core64(M, r_avg=1.5):
    R, Z = grid_coords(NR2, NZ2, H)
    sel = np.sqrt(R ** 2 + Z ** 2) < r_avg
    lam, _ = eig_fields(M)
    lam_c = lam[sel].mean(axis=0)
    return {"core_lams": lam_c.tolist(),
            "core_spread": float(lam_c.max() - lam_c.min()),
            "core_mean": float(lam_c.mean())}


def fire64(M0, e_fn, g_fn, max_iter=6000, dt0=0.02, dt_max=0.2,
           log_every=1000, snaps=(), tag=""):
    w = np.ones((NR2, NZ2))
    w[: NR2 - 1, 1:-1] = cell_weights(NR2, NZ2, H)
    precond = (1.0 / w)[..., None, None]
    M = M0.copy()
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    hist, states = [], [{"it": 0, "t": 0.0, "M": M0.copy()}]
    F = -g_fn(M) * precond * FREE42
    f0 = float(np.max(np.abs(F)))
    t0 = time.time()
    for it in range(1, max_iter + 1):
        P = float(np.sum(F * v))
        if P > 0.0:
            n_up += 1
            vn = np.sqrt(np.sum(v * v))
            fn = np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * (F / max(fn, 1e-300)) * vn
            if n_up > 5:
                dt = min(dt * 1.1, dt_max)
                alpha *= 0.99
        else:
            v[:] = 0.0
            dt *= 0.5
            alpha = 0.1
            n_up = 0
        v += dt * F
        M += dt * v
        F = -g_fn(M) * precond * FREE42
        if it % log_every == 0 or it == max_iter:
            E = e_fn(M)
            hist.append({"it": it, "E": E,
                         "fmax": float(np.max(np.abs(F))), "dt": dt})
            print(f"  {tag} it {it:6d} E {E:12.6f} fmax "
                  f"{hist[-1]['fmax']:.3e} [{time.time() - t0:.0f}s]",
                  flush=True)
            if not np.isfinite(E):
                hist[-1]["stop"] = "non-finite (runaway)"
                break
        if it in snaps:
            states.append({"it": it, "t": float(it), "M": M.copy()})
    return M, states, {"fmax_seed": f0, "trace": hist}


def grad_q2_avg_g(M, G, w4, rho4, nphi=NPHI):
    """phi-averaged kinetic gradient with EXPLICIT grid weights (the
    m5_20_5 grad_q2_avg with the 64 x 128 module constants unpinned)."""
    out = None
    for Gp in gens_phi(G, nphi):
        g = grad_q2(M, Gp, w4=w4, rho4=rho4)
        out = g if out is None else out + g
    return out / nphi


# ================= c0: re-relax per sign + instrument gates ==========
def phase_c0(iters=6000):
    out = {}
    for tag, s in SIGNS.items():
        sg = s * G_T
        M0 = seed64(DELTA, sg)
        print(f"[C0 {tag}] 64x128 hedgehog relax ({iters} iters)",
              flush=True)
        M, _, rel = fire64(M0, lambda MM: e_stat_s(MM, sg),
                           lambda MM: g_stat_s(MM, sg),
                           max_iter=iters, tag=tag)
        np.savez_compressed(
            os.path.join(DATA, f"m5_21_1_c_state_{tag}.npz"), M=M)
        out[tag] = {"E": e_stat_s(M, sg), "q": q_meridional64(M),
                    "core": core64(M), "fire_tail": rel["trace"][-1],
                    "descent_ok": bool(rel["trace"][-1]["E"]
                                       < rel["trace"][0]["E"])}
        print(f"[C0 {tag}] E {out[tag]['E']:.5f} q {out[tag]['q']:.4f} "
              f"core {np.round(out[tag]['core']['core_lams'], 4)}",
              flush=True)
    # GC1 instrument gates on the sign-+ hedgehog background
    M = np.load(os.path.join(DATA, "m5_21_1_c_state_sp.npz"))["M"]
    rng = np.random.default_rng(3)
    G = conj_gen("J23", "K2", 0.5)
    worst = 0.0
    for _ in range(2):
        D = rng.normal(size=M.shape)
        D = 0.5 * (D + np.swapaxes(D, -1, -2))
        Dm = np.zeros_like(M)
        Dm[: NR2 - 1, 1:-1] = D[: NR2 - 1, 1:-1]
        lhs = float(np.sum(grad_q2_avg_g(M, G, W42, RHO42) * Dm))
        rhs = _cs_dir(lambda MM: q2_avg_n(MM, G, NPHI), M, Dm)
        worst = max(worst, abs(lhs - rhs) / max(abs(rhs), 1e-12))
    q16 = q2_avg_f(M, G, 16)
    gates = {
        "gc1_grad_cs": worst,
        "gc1_nphi_5_16": abs(q2_avg_f(M, G, 5) - q16) / abs(q16),
        "gc1_j12_rel": abs(q2_avg_f(M, gen("J12"), 5)
                           - q2_of(M, gen("J12")))
        / max(abs(q2_of(M, gen("J12"))), 1e-12)}
    gates["ok"] = bool(gates["gc1_grad_cs"] < 1e-8
                       and gates["gc1_nphi_5_16"] < 1e-10
                       and gates["gc1_j12_rel"] < 1e-12)
    out["GC1"] = gates
    print(f"[GC1] {'PASS' if gates['ok'] else 'FAIL'} cs "
          f"{gates['gc1_grad_cs']:.1e} nphi {gates['gc1_nphi_5_16']:.1e}"
          f" j12 {gates['gc1_j12_rel']:.1e}", flush=True)
    return out


# ================= c1: the a2x alignment read =================
def _sector_norms(G):
    t = float(np.sqrt(np.sum(G[..., 0, :] ** 2)
                      + np.sum(G[..., 1:, 0] ** 2)))
    s = float(np.sqrt(np.sum(G[..., 1:, 1:] ** 2)))
    return {"time_row": t, "spatial_block": s}


def a2x_row(M, G, sg, w4, rho4, free4):
    gk = grad_q2_avg_g(M, G, w4, rho4) * free4
    gs = grad_static_4(M, WSCALE, DELTA, g=sg, w4=w4, rho4=rho4) * free4
    nk = float(np.sqrt(np.sum(gk ** 2)))
    ns = float(np.sqrt(np.sum(gs ** 2)))
    dot = float(np.sum(gk * gs))
    cos = dot / max(nk * ns, 1e-300)
    w2_opt = dot / max(nk ** 2, 1e-300)
    return {"cos_align": cos, "norm_gkin": nk, "norm_gstat": ns,
            "omega2_optimal": w2_opt,
            "best_rel_residual_any_omega":
                float(np.sqrt(max(1.0 - cos ** 2, 0.0))),
            "sectors_gkin": _sector_norms(gk),
            "sectors_gstat": _sector_norms(gs)}


def phase_c1():
    out = {"states": {}}
    gens = {"J12_pure": gen("J12"), "J23_pure": gen("J23")}
    for name, (rot, boost, chi) in FAMS.items():
        gens[f"{name}_chi{chi}"] = conj_gen(rot, boost, chi)
    for tag, s in SIGNS.items():
        sg = s * G_T
        M = np.load(os.path.join(DATA,
                                 f"m5_21_1_c_state_{tag}.npz"))["M"]
        U = e_stat_s(M, sg)
        rows = {}
        for gname, G in gens.items():
            r = a2x_row(M, G, sg, W42, RHO42, FREE42)
            r["Q2_avg"] = q2_avg_f(M, G)
            r["root_omega"] = (float(np.sqrt(-U / r["Q2_avg"]))
                               if U * r["Q2_avg"] < 0 else None)
            rows[gname] = r
            print(f"[C1 {tag} {gname}] cos {r['cos_align']:+.4f} "
                  f"Q2 {r['Q2_avg']:+.4e} |gk| {r['norm_gkin']:.3e} "
                  f"|gs| {r['norm_gstat']:.3e} w2opt "
                  f"{r['omega2_optimal']:+.3e}", flush=True)
        # chi scan per family (boost-adjacent context, labeled)
        scans = {}
        for name, (rot, boost, _) in FAMS.items():
            chis = np.linspace(0.0, 2.4, 13)
            q2s = [q2_avg_f(M, conj_gen(rot, boost, float(c)))
                   for c in chis]
            cross = None
            for i in range(len(chis) - 1):
                if q2s[i] * q2s[i + 1] < 0:
                    cross = float(chis[i] - q2s[i]
                                  * (chis[i + 1] - chis[i])
                                  / (q2s[i + 1] - q2s[i]))
                    break
            scans[name] = {"chi": chis.tolist(), "Q2_avg": q2s,
                           "chi_crossing": cross}
            print(f"[C1 {tag} scan {name}] Q2(0) {q2s[0]:+.3e} "
                  f"crossing chi = {cross}", flush=True)
        out["states"][tag] = {"U": U, "rows": rows, "chi_scans": scans}
    # cross-grid consistency: the 128 x 256 P1 endpoint (sign +)
    fn = os.path.join(DATA, "m5_21_1_b_endpoint.npz")
    if os.path.exists(fn):
        Mb = np.load(fn)["M"]
        nrb, nzb = Mb.shape[:2]
        w4b = cell_weights(nrb, nzb, H)[..., None, None]
        rho4b = ((np.arange(nrb - 1) + 0.5) * H)[:, None, None, None]
        free4b = (~pin_mask(nrb, nzb))[..., None, None].astype(float)
        big = {}
        for gname in ("J23_pure", "J23^K2_chi0.5"):
            G = gens[gname]
            r = a2x_row(Mb, G, G_T, w4b, rho4b, free4b)
            r["Q2_avg"] = float(sum(
                t_total_c(Mb, d_g(Mb, Gp)).real
                for Gp in gens_phi(G, NPHI)) / NPHI)
            big[gname] = r
            print(f"[C1 big {gname}] cos {r['cos_align']:+.4f} "
                  f"Q2 {r['Q2_avg']:+.4e}", flush=True)
        out["big_grid_crosscheck"] = big
    return out


# ================= c2: rotating extremal (Routhian descent) ==========
def phase_c2(omegas=OMEGAS, iters=2500, film_omega=0.1349):
    out = {"omegas": {}, "note":
           "R(M; Omega) = U - Omega^2 Q2_avg on the PURE rotation "
           "J23 class (boost quarantined); FIRE descent from the c0 "
           "sign-+ state"}
    G = gen("J23")
    M_start = np.load(os.path.join(DATA, "m5_21_1_c_state_sp.npz"))["M"]
    U0 = e_stat_s(M_start, G_T)
    for om in omegas:
        def e_fn(MM):
            return e_stat_s(MM, G_T) - om ** 2 * q2_avg_f(MM, G)

        def g_fn(MM):
            return (g_stat_s(MM, G_T)
                    - om ** 2 * grad_q2_avg_g(MM, G, W42, RHO42))
        snaps = ((0, 300, 800, 1500, iters)
                 if abs(om - film_omega) < 1e-9 else ())
        print(f"[C2] Omega = {om} descent ({iters} iters)", flush=True)
        M, states, rel = fire64(M_start, e_fn, g_fn, max_iter=iters,
                                log_every=500, snaps=snaps,
                                tag=f"om{om}")
        U = e_stat_s(M, G_T)
        Q2 = q2_avg_f(M, G)
        gs = g_stat_s(M, G_T) * FREE42
        gr = (gs - om ** 2 * grad_q2_avg_g(M, G, W42, RHO42) * FREE42)
        rel_resid = (float(np.sqrt(np.sum(gr ** 2)))
                     / max(float(np.sqrt(np.sum(gs ** 2))), 1e-300))
        row = {"omega": om, "U_start": U0, "U": U, "Q2_avg": Q2,
               "T_rot": om ** 2 * Q2, "E_tot": U + om ** 2 * Q2,
               "J": 2.0 * om * Q2, "q": q_meridional64(M),
               "core": core64(M), "rel_residual_gradR": rel_resid,
               "fire_tail": rel["trace"][-1],
               "runaway_nonfinite": any("stop" in r
                                        for r in rel["trace"]),
               # centrifugal runaway: R unbounded below at this Omega
               # (U blowing up while R descends; the first run measured
               # exactly this at Omega >= 0.1349, flag added post-hoc)
               "centrifugal_runaway": bool(U > 10.0 * abs(U0) + 1.0)}
        out["omegas"][str(om)] = row
        print(f"[C2 om {om}] U {U:.5f} T {row['T_rot']:.5f} "
              f"J {row['J']:.4f} q {row['q']:.4f} relresid "
              f"{rel_resid:.3e} centrifugal_runaway="
              f"{row['centrifugal_runaway']}", flush=True)
        np.savez_compressed(
            os.path.join(DATA, f"m5_21_1_c_rot_om{om}.npz"), M=M)
        if snaps:
            for st in states:
                st["vac_spec"] = (G_T, 1.0, DELTA, 0.0)
            film_strip(states,
                       os.path.join(PLOTS, "m5_21_1_c_film_basic.png"),
                       template="basic", delta=DELTA, h=H, g=G_T,
                       wscale=WSCALE,
                       suptitle=f"M5.21.1-P2 Routhian descent at Omega "
                                f"= {om} (rotating-frame relax frames)")
            film_strip(states,
                       os.path.join(PLOTS,
                                    "m5_21_1_c_film_thermal.png"),
                       template="thermal", delta=DELTA, h=H, g=G_T,
                       suptitle=f"M5.21.1-P2 Routhian descent at Omega "
                                f"= {om} (thermal panel set)", step=6)
    # sign spot-check at the film omega
    om = film_omega
    Msm = np.load(os.path.join(DATA, "m5_21_1_c_state_sm.npz"))["M"]

    def e_fn_m(MM):
        return e_stat_s(MM, -G_T) - om ** 2 * q2_avg_f(MM, G)

    def g_fn_m(MM):
        return (g_stat_s(MM, -G_T)
                - om ** 2 * grad_q2_avg_g(MM, G, W42, RHO42))
    print(f"[C2 sm] sign-(-) spot-check at Omega = {om}", flush=True)
    Mm, _, relm = fire64(Msm, e_fn_m, g_fn_m, max_iter=iters,
                         log_every=500, tag="sm")
    Um, Q2m = e_stat_s(Mm, -G_T), q2_avg_f(Mm, G)
    out["sign_minus_spot"] = {
        "omega": om, "U": Um, "Q2_avg": Q2m, "J": 2.0 * om * Q2m,
        "q": q_meridional64(Mm), "fire_tail": relm["trace"][-1]}
    print(f"[C2 sm] U {Um:.5f} J {out['sign_minus_spot']['J']:.4f} "
          f"q {out['sign_minus_spot']['q']:.4f}", flush=True)
    return out


# ================= c3: KG twist dispersion =================
def hv_cs2(M0, V, sg, w4, rho4, free4, eps=1e-30):
    G = grad_static_4(M0.astype(complex) + 1j * eps * V, WSCALE, DELTA,
                      g=sg, w4=w4, rho4=rho4)
    return np.imag(G) / eps * free4


def phase_c3(ks_n=(0, 1, 2, 3, 4, 6, 8), r_env=20.0):
    out = {"grids": {}}
    specs = [("64x128", "m5_21_1_c_state_sp.npz", NR2, NZ2, W42, RHO42,
              FREE42)]
    fnb = os.path.join(DATA, "m5_21_1_b_endpoint.npz")
    if os.path.exists(fnb):
        nrb, nzb = 128, 256
        specs.append(("128x256", "m5_21_1_b_endpoint.npz", nrb, nzb,
                      cell_weights(nrb, nzb, H)[..., None, None],
                      ((np.arange(nrb - 1) + 0.5) * H)[:, None, None,
                                                       None],
                      (~pin_mask(nrb, nzb))[..., None,
                                            None].astype(float)))
    for gtag, fn, nr, nz, w4, rho4, free4 in specs:
        M = np.load(os.path.join(DATA, fn))["M"]
        R, Z = grid_coords(nr, nz, H)
        r = np.sqrt(R ** 2 + Z ** 2)
        env = np.exp(-((r / r_env) ** 4))
        _, V = eig_fields(M)
        W = local_gen(V[..., :, 0])
        base = W @ M - M @ W
        w2 = cell_weights(nr, nz, H)
        rows = []
        for n in ks_n:
            k = 2.0 * np.pi * n / (nz * H)
            prof = (env * np.cos(k * Z))[..., None, None]
            v = prof * base
            v = v * free4
            T2 = float(t_total_c(M, v).real)
            hv = hv_cs2(M, v, G_T, w4, rho4, free4)
            D2E = float(np.sum(v * hv))
            vc = v[: nr - 1, 1:-1]
            can = float(np.sum(w2 * np.sum(vc * vc, axis=(-2, -1))))
            rows.append({"n": n, "k": k, "T_true": T2, "D2E": D2E,
                         "w2_true": D2E / (2.0 * T2) if T2 > 0
                         else None,
                         "w2_canonical": D2E / (2.0 * 0.5 * can)})
            print(f"[C3 {gtag} n={n}] k {k:.4f} T_true {T2:.4e} "
                  f"D2E {D2E:.4e} w2_true "
                  f"{rows[-1]['w2_true']}", flush=True)
        ktab = np.array([r["k"] for r in rows])
        fits = {}
        for key in ("w2_true", "w2_canonical"):
            w2s = np.array([r[key] if r[key] is not None else np.nan
                            for r in rows])
            ok = np.isfinite(w2s)
            if ok.sum() >= 3:
                A = np.vstack([np.ones(ok.sum()), ktab[ok] ** 2]).T
                coef, res, *_ = np.linalg.lstsq(A, w2s[ok], rcond=None)
                pred = A @ coef
                ss = 1.0 - (np.sum((w2s[ok] - pred) ** 2)
                            / max(np.sum((w2s[ok]
                                          - w2s[ok].mean()) ** 2),
                                  1e-300))
                fits[key] = {"m2": float(coef[0]),
                             "c2": float(coef[1]), "R2": float(ss),
                             "m": float(np.sqrt(max(coef[0], 0.0))),
                             "c": float(np.sqrt(max(coef[1], 0.0)))}
        gc2 = all(r["T_true"] > 0 for r in rows)
        out["grids"][gtag] = {"rows": rows, "fits": fits,
                              "GC2_positive_kinetic": bool(gc2)}
        for key, f in fits.items():
            print(f"[C3 {gtag} fit {key}] m {f['m']:.4f} c {f['c']:.4f}"
                  f" R2 {f['R2']:.5f}", flush=True)
    # dispersion figure
    fig, axes = plt.subplots(1, len(out["grids"]),
                             figsize=(6.0 * len(out["grids"]), 4.2),
                             squeeze=False)
    for ax, (gtag, gg) in zip(axes[0], out["grids"].items()):
        ktab = np.array([r["k"] for r in gg["rows"]])
        for key, mk in (("w2_true", "o"), ("w2_canonical", "s")):
            w2s = np.array([r[key] if r[key] is not None else np.nan
                            for r in gg["rows"]])
            ax.plot(ktab ** 2, w2s, mk, label=key)
            if key in gg["fits"]:
                f = gg["fits"][key]
                kk = np.linspace(0, ktab.max() ** 2, 50)
                ax.plot(kk, f["m2"] + f["c2"] * kk, "--", lw=0.8)
        ax.set_xlabel("k^2")
        ax.set_ylabel("omega^2")
        ax.set_title(f"KG twist dispersion ({gtag})", fontsize=9)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_21_1_c_kg.png"), dpi=130)
    plt.close(fig)
    print("wrote plots/m5_21_1_c_kg.png")
    return out


def a2x_plot(res):
    if "C1" not in res:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    st = res["C1"]["states"]
    names = list(next(iter(st.values()))["rows"].keys())
    x = np.arange(len(names))
    for off, tag in ((-0.15, "sp"), (0.15, "sm")):
        cs = [st[tag]["rows"][n]["cos_align"] for n in names]
        axes[0].bar(x + off, cs, width=0.3, label=f"sign {tag}")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(names, rotation=30, fontsize=7)
    axes[0].axhline(0, color="k", lw=0.6)
    axes[0].set_title("a2x alignment cos(g_kin, g_stat) per generator",
                      fontsize=9)
    axes[0].legend(fontsize=7)
    for tag, ls in (("sp", "-"), ("sm", "--")):
        for name, sc in st[tag]["chi_scans"].items():
            axes[1].plot(sc["chi"], sc["Q2_avg"], ls, lw=1,
                         label=f"{tag} {name}")
    axes[1].axhline(0, color="k", lw=0.6)
    axes[1].set_xlabel("chi")
    axes[1].set_ylabel("Q2_avg")
    axes[1].set_title("Q2_avg(chi) per family (boost-adjacent context)",
                      fontsize=9)
    axes[1].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_21_1_c_a2x.png"), dpi=130)
    plt.close(fig)
    print("wrote plots/m5_21_1_c_a2x.png")


def main(which="all", arg1=None):
    os.makedirs(DATA, exist_ok=True)
    fn = os.path.join(DATA, "m5_21_1_c_4d.json")
    res = {"task": "M5.21.1", "phase": "P2+P3", "delta": DELTA,
           "grid_p2": [NR2, NZ2, H], "nphi": NPHI, "wscale": WSCALE}
    if os.path.exists(fn):
        with open(fn) as f:
            res.update(json.load(f))
    if which in ("c0", "all"):
        res["C0"] = phase_c0(int(arg1) if which == "c0" and arg1
                             else 6000)
    if which in ("c1", "all"):
        res["C1"] = phase_c1()
        a2x_plot(res)
    if which in ("c2", "all"):
        oms = (tuple(float(x) for x in arg1.split(","))
               if which == "c2" and arg1 else OMEGAS)
        res["C2"] = phase_c2(omegas=oms)
    if which in ("c3", "all"):
        res["C3"] = phase_c3()
    with open(fn, "w") as f:
        json.dump(res, f, indent=1, default=float)
    print("wrote data/m5_21_1_c_4d.json")
    return res


if __name__ == "__main__":
    main(ARGV[0] if ARGV else "all", ARGV[1] if len(ARGV) > 1 else None)
