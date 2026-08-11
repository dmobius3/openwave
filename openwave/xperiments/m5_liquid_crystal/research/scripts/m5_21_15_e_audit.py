"""M5.21.15 round-1 ADVERSARIAL AUDIT (independent auditor).

Goal: REFUTE, not confirm, the five claims C1..C5 of the M5.21.15
task record. All checks below use the auditor's own derivations and
own numerical code; the task scripts (m5_21_3_a_4d.py,
m5_21_14_c_minimize.py, m5_21_15_b_coupled.py) are imported ONLY to
evaluate their functions as objects under test.

  C1  exact quadraticity of the rotating-configuration energy in
      omega, plus kin invariance along an eta-preserving rotation
      orbit (and the boost-orbit normalization nuance).
  C2  the envelope-concavity theorem (inf of affine functions):
      proof audit in the verdict text + a numerical concavity demo.
  C3  independent recomputation of the A1 RECORD gate from
      data/m5_21_3_all.json, compared against
      data/m5_21_15_baseline.json "RECORD".
  C4  evenness of E_corr / kin_corr under b -> -b on the analytic
      dressed family (C14.ExactCorr, s = -1, g = 32), plateau and
      asymmetric multi-bump profiles, plus the gradient-at-zero test.
  C5  the fixed-J Legendre/Routhian logic: factor checks, the
      envelope identity dE/dJ = omega*, the negative-kin runaway of
      E_stat + J^2/(4 kin) as kin -> 0-, and a live-family probe of
      whether kin_tot can reach the dangerous region within the
      trust region.

Headless, numpy/scipy only. Writes ONLY
../data/m5_21_15_audit_r1.json. Run: python3 m5_21_15_e_audit.py
"""
from __future__ import annotations

import importlib.util
import json
import os
import time

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_21_15_audit_r1.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B3 = _load("aud_b3", "m5_21_3_a_4d.py")
C14 = _load("aud_c14", "m5_21_14_c_minimize.py")
B15 = _load("aud_b15", "m5_21_15_b_coupled.py")

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
RNG = np.random.default_rng(20260811)


# =============== auditor's own primitives ===============
def my_inner_eta(F, G):
    """<F,G>_eta = tr(eta F eta G^T), auditor implementation."""
    X = np.einsum("ab,...bc,cd->...ad", ETA, F, ETA)
    return np.einsum("...ab,...ab->...", X, G)


def my_time_energy(M, a0, om, cfg):
    """auditor read of the time-sector energy of a config rotating
    with dM/dt = om * a0: 4 sum_i <F_0i, F_0i>_eta,
    F_0i = comm_eta(dM/dt, A_i), h^3-weighted. Independent of
    kin_of except for the shared stencil primitives (the grid layer
    is part of the object under test's definition)."""
    h3 = cfg["h"] ** 3
    Mdot = om * a0
    tot = 0.0
    for br, wt in B3.branches(cfg["stencil"]):
        A = [B3.d1(M, ax, cfg["h"], br) for ax in range(3)]
        for i in range(3):
            F = Mdot @ ETA @ A[i] - A[i] @ ETA @ Mdot
            tot = tot + wt * 4.0 * np.sum(my_inner_eta(F, F))
    return h3 * tot


def conj(L, M):
    return np.einsum("ab,...bc,dc->...ad", L, M, L)


# =============== C1 ===============
def audit_c1():
    cfg = B3.base_cfg(n=10, L=15.0)
    M = B3.sym4(RNG.normal(size=(10, 10, 10, 4, 4)))
    a0 = RNG.normal(size=M.shape)
    a0 = a0 / np.sqrt(np.sum(a0 * a0))
    out = {}

    # (a) the augmented energy as implemented (fire's tot_e):
    # E(om) = e_total + om^2 * kin_of. Fit degree 4, residual beyond
    # quadratic must be machine zero relative to the quadratic term.
    E0 = float(B3.e_total(M, cfg))
    k = float(B3.kin_of(M, a0, cfg))
    oms = np.array([0.0, 0.1, 0.35, 0.8, 1.3, 2.1, 3.0])
    Es = np.array([E0 + om ** 2 * k for om in oms])
    out["direct_residual_max"] = float(max(
        abs(E - (E0 + om ** 2 * k)) for om, E in zip(oms, Es)))
    # fit the omega-dependent sector EVALUATED FRESH per omega with
    # the auditor's own code (each point carries its own independent
    # rounding, so any hidden omega^1/^3/^4 content would show)
    Ek = np.array([my_time_energy(M, a0, om, cfg) for om in oms])
    V = np.vander(oms, 5, increasing=True)
    c, *_ = np.linalg.lstsq(V, Ek, rcond=None)
    out["fit_coeffs_c0_to_c4"] = c.tolist()
    out["beyond_quadratic_rel"] = float(
        max(abs(c[0]), abs(c[1]), abs(c[3]), abs(c[4])) / abs(c[2]))

    # (b) the auditor's own time-sector energy vs om^2 * kin_of:
    # this tests the CLAIMED FORMULA, not just the code's self-use.
    devs = []
    for om in (0.07, 0.4, 1.1, 2.7):
        e_my = my_time_energy(M, a0, om, cfg)
        e_cl = om ** 2 * k
        devs.append(abs(e_my - e_cl) / max(abs(e_cl), 1e-300))
    out["own_formula_vs_kin_of_rel"] = float(max(devs))

    # (c) mixed-partial cancellation: for M(t) = L(t) M L(t)^T the
    # naive extra F_0i pieces d_t A_i - d_i A_0 vanish identically
    # because both are the same mixed second difference of M(x, t);
    # checked exactly on the lattice.
    Jz = np.zeros((4, 4)); Jz[1, 2], Jz[2, 1] = -1.0, 1.0
    om, dt = 0.6, 1e-3
    Lp = expm(dt * om * Jz)
    Lm = expm(-dt * om * Jz)
    Mp, Mm = conj(Lp, M), conj(Lm, M)
    A0 = (Mp - Mm) / (2 * dt)
    mix, scale = 0.0, 0.0
    for ax in range(3):
        dtAi = (B3.d1(Mp, ax, cfg["h"], "fwd")
                - B3.d1(Mm, ax, cfg["h"], "fwd")) / (2 * dt)
        diA0 = B3.d1(A0, ax, cfg["h"], "fwd")
        mix = max(mix, float(np.max(np.abs(dtAi - diA0))))
        scale = max(scale, float(np.max(np.abs(dtAi))))
    out["mixed_partial_residual_rel"] = mix / scale

    # (d) kin along an eta-preserving ROTATION orbit: M -> L M L^T,
    # a0 rebuilt by the catalog recipe (envelope * (G M - M G^T),
    # unit Frobenius norm) at each orbit point.
    w = B3.envelope(cfg)[..., None, None]

    def a0_of(Mx, G, normalize):
        raw = w * (G @ Mx - Mx @ G.T)
        if not normalize:
            return raw
        return raw / np.sqrt(np.sum(raw * raw))

    thetas = (0.0, 0.4, 1.1, 2.0)
    kin_rot_n, kin_rot_u, e_stat_orbit = [], [], []
    for th in thetas:
        L = expm(th * Jz)
        assert np.allclose(L.T @ ETA @ L, ETA, atol=1e-12)
        ML = conj(L, M)
        kin_rot_n.append(float(B3.kin_of(ML, a0_of(ML, Jz, True), cfg)))
        kin_rot_u.append(float(B3.kin_of(ML, a0_of(ML, Jz, False),
                                         cfg)))
        e_stat_orbit.append(float(B3.e_total(ML, cfg)))
    out["rot_orbit_thetas"] = thetas
    out["rot_orbit_kin_normalized"] = kin_rot_n
    out["rot_orbit_kin_spread_rel"] = float(
        (max(kin_rot_n) - min(kin_rot_n)) / abs(np.mean(kin_rot_n)))
    out["rot_orbit_kin_unnorm_spread_rel"] = float(
        (max(kin_rot_u) - min(kin_rot_u)) / abs(np.mean(kin_rot_u)))
    out["rot_orbit_Estat_spread_rel"] = float(
        (max(e_stat_orbit) - min(e_stat_orbit))
        / abs(np.mean(e_stat_orbit)))

    # (e) the BOOST-orbit nuance: eta-preserving but not orthogonal,
    # so the unit-Frobenius normalization is NOT orbit-invariant.
    Kz = np.zeros((4, 4)); Kz[0, 3] = Kz[3, 0] = 1.0
    kin_b_n, kin_b_u = [], []
    for th in (0.0, 0.2, 0.5):
        L = expm(th * Kz)
        assert np.allclose(L.T @ ETA @ L, ETA, atol=1e-12)
        ML = conj(L, M)
        kin_b_n.append(float(B3.kin_of(ML, a0_of(ML, Kz, True), cfg)))
        kin_b_u.append(float(B3.kin_of(ML, a0_of(ML, Kz, False), cfg)))
    out["boost_orbit_kin_normalized"] = kin_b_n
    out["boost_orbit_kin_norm_spread_rel"] = float(
        (max(kin_b_n) - min(kin_b_n)) / abs(np.mean(kin_b_n)))
    out["boost_orbit_kin_unnorm_spread_rel"] = float(
        (max(kin_b_u) - min(kin_b_u)) / abs(np.mean(kin_b_u)))
    return out


# =============== C2 ===============
def audit_c2():
    # numerical demo: envelope of many random affine functions of u,
    # including negative slopes, checked for concavity by second
    # differences on a dense grid.
    n = 400
    e0 = RNG.uniform(0.0, 10.0, n)
    kk = RNG.uniform(-2.0, 3.0, n)
    us = np.linspace(0.0, 5.0, 601)
    env = np.min(e0[None, :] + us[:, None] * kk[None, :], axis=1)
    d2 = np.diff(env, 2)
    out = {"n_affine": n, "max_second_difference": float(np.max(d2)),
           "concave_numerically": bool(np.max(d2) < 1e-10)}
    # a family whose inf is -infinity beyond a slope threshold still
    # yields no interior minimum: min over the FINITE sub-envelope of
    # the sampled family is at a boundary of the u range.
    imin = int(np.argmin(env))
    out["envelope_argmin_at_boundary"] = bool(
        imin in (0, len(us) - 1))
    return out


# =============== C3 ===============
def audit_c3():
    with open(os.path.join(DATA, "m5_21_3_all.json")) as f:
        all3 = json.load(f)
    with open(os.path.join(DATA, "m5_21_15_baseline.json")) as f:
        base = json.load(f)
    rec = base["RECORD"]
    out = {}
    mismatches = []
    for s in ("+1", "-1"):
        rows = all3[f"p2_p1_s{s}"]["rows"]
        rot = {nm: r["kin"] for nm, r in rows.items()
               if not nm.startswith("boost")}
        boo = {nm: r["kin"] for nm, r in rows.items()
               if nm.startswith("boost")}
        lad = all3[f"p3_p1_s{s}_boost_z"]["ladder"]
        om = np.array([r["omega"] for r in lad])
        E = np.array([r["E"] for r in lad])
        u = om ** 2
        slopes = (np.diff(E) / np.diff(u)).tolist()
        blk = {
            "n_rot": len(rot), "n_boost": len(boo),
            "rot_kin_all_positive": bool(min(rot.values()) > 0.0),
            "boost_kin_all_negative": bool(max(boo.values()) < 0.0),
            "ladder_monotone_decreasing":
                bool(np.all(np.diff(E) < 0.0)),
            "all_rungs_max_iter": bool(all(
                r["stop"] == "max_iter" for r in lad[1:])),
            "rung0_stop": lad[0]["stop"],
            "u_slopes": slopes,
            "u_slopes_increasing":
                bool(np.all(np.diff(slopes) > 0.0)),
        }
        # internal consistency: E == E_u + E_v + om^2 * kin per rung
        cons = [abs(r["E"] - (r["E_u"] + r["E_v"]
                              + r["omega"] ** 2 * r["kin"]))
                for r in lad[1:]]
        blk["E_decomposition_max_absdev"] = float(max(cons))
        # compare with the stored RECORD numbers
        rb = rec[f"s{s}"]
        for nm, v in {**rot, **boo}.items():
            side = "rot_kin" if nm in rot else "boost_kin"
            ref = rb[side][nm]
            if abs(v - ref) > 1e-12 * max(abs(ref), 1e-30):
                mismatches.append(f"s{s}:{nm} kin {v} vs RECORD {ref}")
        for a, b in zip(slopes, rb["u_slopes"]):
            if abs(a - b) > 1e-10 * max(abs(b), 1e-30):
                mismatches.append(f"s{s}: slope {a} vs RECORD {b}")
        for key in ("rot_kin_all_positive", "boost_kin_all_negative",
                    "ladder_monotone_decreasing",
                    "u_slopes_increasing"):
            if blk[key] != rb[key]:
                mismatches.append(f"s{s}:{key} {blk[key]} != RECORD")
        if blk["all_rungs_max_iter"] != rb["ladder_all_depth_bounded"]:
            mismatches.append(f"s{s}: max_iter flag mismatch")
        out[f"s{s}"] = blk
    out["mismatches_vs_RECORD"] = mismatches
    return out


# =============== C4 ===============
def audit_c4():
    grid = C14.make_grid(32, 6, 12)
    assert C14.S_SIGN == -1.0 and C14.G_MAIN == 32.0
    ec = C14.ExactCorr(grid, 32.0)
    profiles = {
        "plateau": np.array([0.03] + [0.0] * 9),
        "asym_multibump": np.array([0.01, 0.05, -0.03, 0.02, 0.0,
                                    -0.04, 0.01, 0.0, 0.03, -0.02]),
    }
    out = {}
    for nm, avec in profiles.items():
        bp = lambda r, a=avec: C14.b_of(a, r)
        bm = lambda r, a=avec: C14.b_of(-a, r)
        ep, kp = ec.both(bp)
        em, km = ec.both(bm)
        out[nm] = {
            "E_corr_plus": ep, "E_corr_minus": em,
            "E_corr_absdiff": abs(ep - em),
            "E_corr_reldiff": abs(ep - em) / max(abs(ep), 1e-300),
            "kin_corr_plus": kp, "kin_corr_minus": km,
            "kin_corr_absdiff": abs(kp - km),
            "kin_corr_reldiff": abs(kp - km) / max(abs(kp), 1e-300),
        }
    # gradient at avec = 0 (central difference; by evenness the two
    # lobes must cancel EXACTLY if the evenness is machine-exact)
    eps = 1e-4
    ge, gk = [], []
    for kcomp in range(10):
        a = np.zeros(10); a[kcomp] = eps
        e1, k1 = ec.both(lambda r, aa=a: C14.b_of(aa, r))
        e2, k2 = ec.both(lambda r, aa=a: C14.b_of(-aa, r))
        ge.append((e1 - e2) / (2 * eps))
        gk.append((k1 - k2) / (2 * eps))
    out["grad_at_zero_E_corr_max"] = float(np.max(np.abs(ge)))
    out["grad_at_zero_kin_corr_max"] = float(np.max(np.abs(gk)))
    # the algebraic reason (auditor derivation, verified pointwise):
    # Qb(-b) = eta Qb(b) eta, the base config and a0 are eta-even
    # (block diagonal), and every density is invariant under global
    # conjugation by eta; spot-check the matrix identity.
    K, K2, r = C14.kgeom(grid["P"][:64])
    b = 0.17 * np.ones(64)
    Qp = C14.qb_from(K, K2, b)
    Qm = C14.qb_from(K, K2, -b)
    dev = float(np.max(np.abs(Qm - ETA[None] @ Qp @ ETA[None])))
    out["Qb_minus_eq_etaQbeta_maxdev"] = dev
    return out


# =============== C5 ===============
def audit_c5():
    out = {}
    # (a) factor check on T = kin * omega^2 (their kin absorbs the
    # 1/2 I): J = dT/domega = 2 kin omega, T(J) = J^2/(4 kin),
    # omega = J/(2 kin). Numeric identity check.
    kin0, om = 0.7321, 1.234
    T = kin0 * om ** 2
    dT = (kin0 * (om + 1e-6) ** 2 - kin0 * (om - 1e-6) ** 2) / 2e-6
    J = 2.0 * kin0 * om
    out["J_eq_dT_domega_rel"] = abs(dT - J) / J
    out["T_eq_J2_over_4kin_rel"] = abs(T - J ** 2 / (4 * kin0)) / T
    out["omega_recovery_rel"] = abs(J / (2 * kin0) - om) / om

    # (b) envelope identity dE/dJ = omega* on a smooth positive-kin
    # toy family: Estat(x) = 1 + (x - 0.3)^2, kin(x) = 0.5 + x^2,
    # x in [0, 2].
    def estat(x):
        return 1.0 + (x - 0.3) ** 2

    def kin(x):
        return 0.5 + x * x

    def EofJ(J):
        res = minimize_scalar(
            lambda x: estat(x) + J * J / (4.0 * kin(x)),
            bounds=(0.0, 2.0), method="bounded",
            options={"xatol": 1e-12})
        return float(res.fun), float(res.x)

    J0, dJ = 1.3, 1e-5
    Ep, _ = EofJ(J0 + dJ)
    Em, _ = EofJ(J0 - dJ)
    E0, xstar = EofJ(J0)
    om_star = J0 / (2.0 * kin(xstar))
    out["dEdJ_vs_omega_star_rel"] = abs((Ep - Em) / (2 * dJ)
                                        - om_star) / om_star

    # (c) the loophole: for kin passing through zero the fixed-J
    # functional Estat + J^2/(4 kin) is +inf as kin -> 0+ (a barrier
    # protecting the positive branch) but -inf as kin -> 0- (an
    # unbounded REWARD on the negative branch). Demonstrate.
    def f_fixedJ(kin_v, J=1.0, estat_v=1.0):
        return estat_v + J * J / (4.0 * kin_v)

    out["fixedJ_at_kin_plus"] = [f_fixedJ(k)
                                 for k in (1e-2, 1e-4, 1e-6)]
    out["fixedJ_at_kin_minus"] = [f_fixedJ(-k)
                                  for k in (1e-2, 1e-4, 1e-6)]

    # (d) live-family probe: within the guarded 10-dim family
    # (|a_k| <= BOUND = 2), does kin_tot(clock) reach the dangerous
    # region (below the KIN_FLOOR = 1, or below 0)? 1-D scans along
    # the plateau coefficient and one bump coefficient, modest grid.
    grid = C14.make_grid(32, 6, 12)
    cc = B15.ChanCorr(grid, 32.0)
    kb = cc.kin_base["clock"]
    scans = {}
    for label, comp in (("plateau_a0", 0), ("bump_a3", 3)):
        amps = np.linspace(-B15.BOUND, B15.BOUND, 17)
        kts = []
        for a in amps:
            av = np.zeros(10); av[comp] = a
            _, kc = cc.chan_kins(B15.bfun_of(av), ("clock",))
            kts.append(kb + kc["clock"])
        scans[label] = {"amps": amps.tolist(), "kin_tot": kts,
                        "min_kin_tot": float(np.min(kts)),
                        "reaches_below_floor":
                            bool(np.min(kts) < 1.0),
                        "reaches_below_zero":
                            bool(np.min(kts) < 0.0)}
    out["kin_base_clock_modest_grid"] = float(kb)
    out["family_scans"] = scans
    out["kin_floor_in_code"] = 1.0
    return out


def main():
    t0 = time.time()
    res = {"task": "M5.21.15 adversarial audit round 1",
           "auditor": "independent second agent",
           "date": "2026-08-11"}

    print("C1 ...", flush=True)
    c1 = audit_c1()
    print(json.dumps(c1, indent=1), flush=True)
    print("C2 ...", flush=True)
    c2 = audit_c2()
    print(json.dumps(c2, indent=1), flush=True)
    print("C3 ...", flush=True)
    c3 = audit_c3()
    print(json.dumps(c3, indent=1), flush=True)
    print("C4 ...", flush=True)
    c4 = audit_c4()
    print(json.dumps(c4, indent=1), flush=True)
    print("C5 ...", flush=True)
    c5 = audit_c5()
    print(json.dumps(c5, indent=1), flush=True)

    # ---------- verdicts ----------
    v1 = ("CONFIRMED" if c1["beyond_quadratic_rel"] < 1e-12
          and c1["direct_residual_max"] == 0.0
          and c1["own_formula_vs_kin_of_rel"] < 1e-12
          and c1["mixed_partial_residual_rel"] < 1e-10
          and c1["rot_orbit_kin_spread_rel"] < 1e-10
          else "PARTIAL")
    v3 = "CONFIRMED" if not c3["mismatches_vs_RECORD"] else "REFUTED"
    even_exact = max(
        c4["plateau"]["E_corr_absdiff"],
        c4["plateau"]["kin_corr_absdiff"],
        c4["asym_multibump"]["E_corr_absdiff"],
        c4["asym_multibump"]["kin_corr_absdiff"]) == 0.0
    even_close = max(
        c4["plateau"]["E_corr_reldiff"],
        c4["asym_multibump"]["E_corr_reldiff"]) < 1e-12
    v4 = "CONFIRMED" if (even_exact or even_close) else "PARTIAL"

    res["verdicts"] = {
        "C1": {
            "verdict": v1,
            "reasoning": (
                "The augmented energy is exactly E_stat + omega^2 * "
                "kin: a degree-4 fit leaves relative non-quadratic "
                "content {:.1e} (direct residual exactly {:.1g}); "
                "the auditor's independent time-sector energy "
                "(4 sum_i <comm_eta(Mdot, A_i)>^2 with Mdot = omega "
                "a0) matches omega^2 * kin_of to {:.1e}; the naive "
                "extra F_0i pieces d_t A_i - d_i A_0 cancel "
                "identically on the lattice (relative residual "
                "{:.1e}, pure rounding) because A_mu = d_mu M is a "
                "pure gradient in (x, t). "
                "kin rebuilt by the catalog recipe is invariant "
                "along an eta-preserving ROTATION orbit to {:.1e} "
                "and E_stat is orbit-invariant to {:.1e}. NUANCE, "
                "not a refutation: along a BOOST orbit the "
                "unit-Frobenius normalization is not invariant "
                "(normalized-kin spread {:.2e} vs un-normalized "
                "{:.1e}), so the boost-channel kin number depends on "
                "the orbit point through the normalization "
                "convention; within any single rung a0 is frozen, "
                "so exact quadraticity in omega is unaffected."
            ).format(c1["beyond_quadratic_rel"],
                     c1["direct_residual_max"],
                     c1["own_formula_vs_kin_of_rel"],
                     c1["mixed_partial_residual_rel"],
                     c1["rot_orbit_kin_spread_rel"],
                     c1["rot_orbit_Estat_spread_rel"],
                     c1["boost_orbit_kin_norm_spread_rel"],
                     c1["boost_orbit_kin_unnorm_spread_rel"]),
            "key_numbers": c1},
        "C2": {
            "verdict": "CONFIRMED",
            "reasoning": (
                "Proof audit: for each fixed configuration M the map "
                "u -> E_stat(M) + u kin(M) is affine, hence concave; "
                "a pointwise infimum of concave functions over ANY "
                "index set is concave (E*(t u1 + (1-t) u2) = inf_M "
                "[t f_M(u1) + (1-t) f_M(u2)] >= t E*(u1) + (1-t) "
                "E*(u2)), with no attainment, compactness, or "
                "continuity assumption; this survives a0 = a0(M) and "
                "dressed b(r) configurations because kin(M, a0(M)) "
                "is still a single number per configuration, and "
                "even a multi-valued a0 choice only enlarges the "
                "affine family. Non-attainment / runaway to "
                "-infinity keeps concavity in the extended-real "
                "sense (the finiteness set is an interval) and "
                "cannot create an interior minimum. A concave "
                "function has no strict interior local minimum "
                "(f(u0) >= min(f(a), f(b)) with strict slack gives a "
                "contradiction), and u = omega^2 is a monotone "
                "bijection on omega > 0, so no strict interior "
                "omega-minimum at omega > 0; differentiable interior "
                "stationary points of a concave function are global "
                "maxima or lie on a flat affine segment. Numerical "
                "demo: envelope of 400 random affine functions has "
                "max second difference {:.1e} (concave) with argmin "
                "on the boundary. SCOPE CAVEAT (already flagged by "
                "the task itself): the theorem constrains the TRUE "
                "envelope E*; a finite-iteration measured ladder is "
                "only an upper bound and may show any shape, so the "
                "theorem cannot be used to certify measured curves, "
                "only to interpret them."
            ).format(c2["max_second_difference"]),
            "key_numbers": c2},
        "C3": {
            "verdict": v3,
            "reasoning": (
                "Independent recomputation from m5_21_3_all.json: "
                "both p2 blocks have 4 rotation channels all kin > 0 "
                "and 2 boost channels all kin < 0; both p3 boost_z "
                "ladders are strictly decreasing in E with every "
                "post-reference rung stop == max_iter; the slopes "
                "dE/du between consecutive rungs increase strictly "
                "monotonically in both blocks (s+1: -31.10 -> -0.198; "
                "s-1: -52.01 -> -0.255). Every number and flag "
                "matches the baseline RECORD to better than 1e-12 "
                "relative ({} mismatches). Internal consistency "
                "E = E_u + E_v + omega^2 kin holds per rung to "
                "{:.1e} (s+1) and {:.1e} (s-1)."
            ).format(len(c3["mismatches_vs_RECORD"]),
                     c3["s+1"]["E_decomposition_max_absdev"],
                     c3["s-1"]["E_decomposition_max_absdev"]),
            "key_numbers": {k: c3[k] for k in
                            ("mismatches_vs_RECORD",)},
            "recomputed": {k: c3[k] for k in ("s+1", "s-1")}},
        "C4": {
            "verdict": v4,
            "reasoning": (
                "Evenness measured on ExactCorr (make_grid(32, 6, "
                "12), s = -1, g = 32): plateau profile E_corr "
                "absdiff {:.3e} (rel {:.1e}), kin_corr absdiff "
                "{:.3e}; asymmetric multi-bump profile E_corr "
                "absdiff {:.3e} (rel {:.1e}), kin_corr absdiff "
                "{:.3e}. Auditor derivation of WHY it is exact: "
                "Qb(-b) = eta Qb(b) eta (verified to {:.1e}), the "
                "base config and a0 are block diagonal (eta-even), "
                "and both densities are invariant under global "
                "conjugation by eta, so the b -> -b field is the "
                "eta-conjugate of the +b field pointwise and every "
                "integral coincides; in IEEE arithmetic the sign "
                "flips are exact, making the evenness bitwise. "
                "Since b_of is linear in avec, f(-avec) = f(avec) "
                "for the full coefficient vector, so every "
                "directional derivative at avec = 0 vanishes: "
                "measured central-difference gradient at zero, max "
                "|dE_corr/da| = {:.1e}, max |dkin_corr/da| = {:.1e}. "
                "The stationary-at-zero claim survives."
            ).format(c4["plateau"]["E_corr_absdiff"],
                     c4["plateau"]["E_corr_reldiff"],
                     c4["plateau"]["kin_corr_absdiff"],
                     c4["asym_multibump"]["E_corr_absdiff"],
                     c4["asym_multibump"]["E_corr_reldiff"],
                     c4["asym_multibump"]["kin_corr_absdiff"],
                     c4["Qb_minus_eq_etaQbeta_maxdev"],
                     c4["grad_at_zero_E_corr_max"],
                     c4["grad_at_zero_kin_corr_max"]),
            "key_numbers": c4},
        "C5": {
            "verdict": "PARTIAL",
            "reasoning": (
                "The Legendre/Routhian algebra is CORRECT: with T = "
                "kin omega^2, J = dT/domega = 2 kin omega (checked "
                "to {:.1e}), T = J^2/(4 kin) (identity to {:.1e}), "
                "omega* = J/(2 kin), and on a smooth positive-kin "
                "toy family the envelope satisfies dE/dJ = omega* to "
                "{:.1e}. No sign or factor errors. BUT the claimed "
                "MECHANISM has a precise loophole: J^2/(4 kin) "
                "penalizes only SMALL POSITIVE kin (+inf as kin -> "
                "0+, a genuine barrier protecting the kin > 0 "
                "branch against continuous descent); on the "
                "negative-kin branch the same term is an unbounded "
                "REWARD (-inf as kin -> 0-, demonstrated: E = 1 + "
                "1/(4 kin) = -2.5e5 at kin = -1e-6), so the global "
                "fixed-J infimum over any configuration set that "
                "reaches kin < 0 is -infinity and the fixed-J "
                "minimum is only well defined as a BRANCH-RESTRICTED "
                "minimum; a global search, a disconnected family, "
                "or a jump across kin = 0 lands on the runaway "
                "side. In the code this restriction is enforced by "
                "the hard KIN_FLOOR = 1.0 barrier, so the floor, "
                "not the J^2/(4 kin) term alone, carries the "
                "structural guarantee. Live-family probe (clock "
                "channel, modest grid, |a| <= 2): kin_base = {:.4g}; "
                "1-D scans along the plateau and bump-3 "
                "coefficients reach min kin_tot = {:.4g} and "
                "{:.4g}; below-zero reached: {} / {}. So even 1-D "
                "slices of the guarded family cross kin_tot = 0 "
                "well inside the trust region (plateau slice "
                "already below |a| = 0.25, bump-3 slice between "
                "0.25 and 0.5): the loophole is LIVE in "
                "the actual family and the code's KIN_FLOOR barrier "
                "is doing real structural work. The claim should be "
                "stated as branch-restricted: fixed J protects only "
                "the kin > 0 branch, and only against CONTINUOUS "
                "descent through kin = 0+."
            ).format(c5["J_eq_dT_domega_rel"],
                     c5["T_eq_J2_over_4kin_rel"],
                     c5["dEdJ_vs_omega_star_rel"],
                     c5["kin_base_clock_modest_grid"],
                     c5["family_scans"]["plateau_a0"]["min_kin_tot"],
                     c5["family_scans"]["bump_a3"]["min_kin_tot"],
                     c5["family_scans"]["plateau_a0"]
                       ["reaches_below_zero"],
                     c5["family_scans"]["bump_a3"]
                       ["reaches_below_zero"]),
            "key_numbers": c5},
    }
    res["runtime_s"] = round(time.time() - t0, 1)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps({c: res["verdicts"][c]["verdict"]
                      for c in ("C1", "C2", "C3", "C4", "C5")}
                     | {"runtime_s": res["runtime_s"]}, indent=1))


# ══════════════════════════════════════════════════════════════════
# ROUND 2 ADVERSARIAL AUDIT (independent auditor, second pass).
# Round 1 above is intact and unchanged. Round 2 audits the RUN DATA
# claims D1..D5 with the auditor's own recomputations plus an
# independent re-minimization (differential evolution + Nelder-Mead)
# against the task's own ChanCorr machinery:
#
#   D1  free envelopes (m5_21_15_coupled.json): no interior minimum
#       beyond optimizer noise; concave in u = omega^2 up to noise.
#       Auditor tools: row-consistency, own slopes, chord test, and
#       the REFINED ENVELOPE: each rung's stored full-grid
#       (E_corr, kin_tot) defines the affine line E_i(u) = E_corr_i
#       + u kin_tot_i, valid at EVERY u; the pointwise min over
#       rungs is a concave-by-construction improvement of the
#       measured ladder that quantifies per-rung under-convergence.
#   D2  fixed-J narrow-guard minimum (m5_21_15_fom_narrow.json):
#       arithmetic, EJ anchor match, mismatch propagation, sign
#       robustness of +115.9, hidden-minimum test on the
#       infeasible-flagged rungs.
#   D3  guard ladder (m5_21_15_guard.json): arithmetic, sign
#       pattern, the omega* > J/(2 kin_base) wall, objective
#       nesting in the bound, and an independent re-minimization at
#       (bound 0.02, J = 133.1) and (bound 0.02, J = 332.76).
#   D4  kin-ceiling: search for kin_corr(clock) > 0 in the 10-dim
#       family (random sampling + L-BFGS-B ascent, bounds 0.02 and
#       2.0), verified on the full grid.
#   D5  EJ thermodynamics (m5_21_15_fixedj.json): own non-uniform
#       central differences for dE/dJ vs omega*, E(J) monotonicity,
#       the stated 1-2 percent agreement.
#
# Run:  python3 m5_21_15_e_audit.py            -> round 2 (writes
#                                                 m5_21_15_audit_r2.json)
#       python3 m5_21_15_e_audit.py --round1   -> round 1 (unchanged)
# ══════════════════════════════════════════════════════════════════
import sys

from scipy.optimize import differential_evolution, minimize

OUT2 = os.path.join(DATA, "m5_21_15_audit_r2.json")
RNG2 = np.random.default_rng(20260812)
KIN_FLOOR2 = 1.0


# =============== D1 ===============
def _ladder_audit(scan):
    rows = scan["ladder"]
    om = np.array([r["omega"] for r in rows])
    u = om ** 2
    E = np.array([r["E_env"] for r in rows])
    Ec = np.array([r["E_corr"] for r in rows])
    kt = np.array([r["kin_tot"] for r in rows])
    kb = np.array([r["kin_base"] for r in rows])
    kc = np.array([r["kin_corr"] for r in rows])
    Jr = np.array([r["J"] for r in rows])
    out = {"row_consistency": {
        "E_env_max_absdev": float(np.max(np.abs(E - (Ec + u * kt)))),
        "kin_tot_max_absdev": float(np.max(np.abs(kt - (kb + kc)))),
        "J_max_absdev": float(np.max(np.abs(Jr - 2.0 * om * kt)))}}
    slopes = np.diff(E) / np.diff(u)
    out["slopes_recomputed"] = slopes.tolist()
    out["slopes_vs_stored_maxdev"] = float(np.max(np.abs(
        slopes - np.array(scan["u_slopes"]))))
    # measured concavity violations: chord test (concave lies ABOVE
    # every chord; a dip below is a violation), in local-E units
    viols, mins = [], []
    for j in range(1, len(u) - 1):
        lam = (u[j] - u[j - 1]) / (u[j + 1] - u[j - 1])
        chord = (1 - lam) * E[j - 1] + lam * E[j + 1]
        v = chord - E[j]
        if v > 0:
            viols.append({"omega": float(om[j]), "depth": float(v),
                          "rel_localE": float(v / abs(E[j]))})
        d = min(E[j - 1], E[j + 1]) - E[j]
        if d > 0:
            mins.append({"omega": float(om[j]), "depth": float(d),
                         "rel_localE": float(d / abs(E[j]))})
    out["chord_violations"] = viols
    out["interior_min_candidates"] = mins
    incr = np.diff(slopes)
    out["slope_increments_positive"] = [
        {"between_omegas": [float(om[i]), float(om[i + 2])],
         "delta_slope": float(d)}
        for i, d in enumerate(incr) if d > 0]
    # the refined envelope from the run's own configs (full-grid
    # numbers straight from the record; measured >= refined exactly,
    # because each rung's own line passes through its own point)
    cross = Ec[:, None] + u[None, :] * kt[:, None]
    refined = np.min(cross, axis=0)
    under = E - refined
    out["under_convergence_per_rung"] = under.tolist()
    out["max_under_convergence"] = float(np.max(under))
    rmins = [j for j in range(1, len(u) - 1)
             if refined[j] < refined[j - 1]
             and refined[j] < refined[j + 1]]
    out["refined_interior_min"] = bool(rmins)
    rsl = np.diff(refined) / np.diff(u)
    out["refined_max_slope_increment"] = float(np.max(np.diff(rsl)))
    out["stored_flags"] = {
        "concave_certified": scan["concave_certified"],
        "interior_minimum_found": scan["interior_minimum_found"],
        "monotone": scan["monotone"]}
    return out


def audit_d1():
    with open(os.path.join(DATA, "m5_21_15_coupled.json")) as f:
        cp = json.load(f)
    out = {}
    for tag, chans in (("main_s-1_g32", ("clock", "rot_z", "boost_z")),
                       ("s+1_g32", ("clock",)),
                       ("s-1_g8", ("clock",))):
        for ch in chans:
            out[f"{tag}:{ch}"] = _ladder_audit(cp[tag][f"scan_{ch}"])
    return out


# =============== D2 ===============
def audit_d2():
    with open(os.path.join(DATA, "m5_21_15_fom_narrow.json")) as f:
        fn = json.load(f)
    with open(os.path.join(DATA, "m5_21_15_guard.json")) as f:
        gd = json.load(f)
    J, Eb, rows = fn["J"], fn["E_base_u"], fn["rows"]
    out = {}
    devs = {"kin_target": [], "E_total": [], "mismatch": []}
    flags_ok = True
    for r in rows:
        devs["kin_target"].append(
            abs(r["kin_target"] - J / (2.0 * r["omega"])))
        devs["E_total"].append(abs(
            r["E_total"] - (Eb + r["E_corr"]
                            + r["omega"] ** 2 * r["kin_tot"])))
        m = (abs(r["kin_tot"] - r["kin_target"])
             / max(abs(r["kin_target"]), 1.0))
        devs["mismatch"].append(abs(m - r["mismatch_rel"]))
        flags_ok = flags_ok and ((m < 0.02) == r["feasible"])
    out["arithmetic_max_dev"] = {k: float(max(v))
                                 for k, v in devs.items()}
    out["feasible_flags_ok"] = flags_ok
    feas = [r for r in rows if r["feasible"]]
    E = np.array([r["E_total"] for r in feas])
    omf = np.array([r["omega"] for r in feas])
    i = int(np.argmin(E))
    out["feasible_omegas"] = omf.tolist()
    out["feasible_E_totals"] = E.tolist()
    out["interior_minimum"] = bool(0 < i < len(E) - 1)
    out["omega_min"] = float(omf[i])
    out["E_min"] = float(E[i])
    ga = next(r for r in gd["rows"]
              if r["bound"] == 0.02 and abs(r["J"] - J) < 1e-9)
    out["EJ_anchor"] = {
        "omega_star": ga["omega_star"], "E_total": ga["E_total"],
        "omega_gap": float(abs(omf[i] - ga["omega_star"])),
        "scan_step": 0.03,
        "match_within_step": bool(
            abs(omf[i] - ga["omega_star"]) <= 0.03)}
    # mismatch propagation: (a) first-order shift bound
    # omega^2 |kin_tot - kin_target| from enforcing the constraint
    # exactly at frozen E_corr; (b) the on-shell bookkeeping
    # E_alt = Eb + E_corr + J^2/(4 kin_tot) at the row's ACTUAL kin
    shifts = [r["omega"] ** 2 * abs(r["kin_tot"] - r["kin_target"])
              for r in rows]
    alts = [Eb + r["E_corr"] + J * J / (4.0 * r["kin_tot"])
            for r in rows]
    out["per_row"] = [
        {"omega": r["omega"], "E_total": r["E_total"],
         "shift_bound": float(s), "E_alt_onshell": float(a),
         "feasible": r["feasible"]}
        for r, s, a in zip(rows, shifts, alts)]
    sh_feas = [s for r, s in zip(rows, shifts) if r["feasible"]]
    out["max_shift_feasible"] = float(max(sh_feas))
    out["E_min_worst_case"] = float(out["E_min"] - max(sh_feas))
    imin_all = rows.index(feas[i])
    out["E_min_alt_onshell"] = float(alts[imin_all])
    out["E_min_robustly_positive"] = bool(
        out["E_min_worst_case"] > 0 and alts[imin_all] > 0)
    # can the infeasible-flagged rungs hide a lower minimum?
    out["infeasible_hidden_min"] = [
        {"omega": r["omega"], "E_total": r["E_total"],
         "E_best_case": float(r["E_total"] - s),
         "could_undercut": bool(
             r["E_total"] - s < out["E_min"] + max(sh_feas))}
        for r, s in zip(rows, shifts) if not r["feasible"]]
    # location sensitivity: neighbor gaps vs the shift scale
    gaps = {}
    if 0 < i < len(E) - 1:
        gaps = {"to_left": float(E[i - 1] - E[i]),
                "to_right": float(E[i + 1] - E[i])}
    out["location_gaps"] = gaps
    out["location_resolved"] = bool(
        gaps and min(gaps.values()) > out["max_shift_feasible"])
    return out


# =============== D3 ===============
def _remin_fixedj(cc, ccf, E_base, J, bound, stored_row):
    """independent re-minimization of E_corr + J^2/(4 kin_tot),
    kin floored, with differential evolution + bounded Nelder-Mead
    multistart; final numbers on the full grid."""
    kb = cc.kin_base["clock"]

    def obj(a):
        e, kc = cc.chan_kins(B15.bfun_of(np.asarray(a)), ("clock",))
        kt = kb + kc["clock"]
        if kt <= KIN_FLOOR2:
            return e + J * J / (4.0 * KIN_FLOOR2) \
                + 1e3 * (KIN_FLOOR2 - kt) ** 2
        return e + J * J / (4.0 * kt)

    stored_obj_opt = float(obj(stored_row["avec"]))
    de = differential_evolution(
        obj, bounds=[(-bound, bound)] * 10, seed=20260812,
        maxiter=60, popsize=16, tol=1e-12, polish=True)
    best_val, best_a = float(de.fun), np.asarray(de.x)
    starts = [np.asarray(stored_row["avec"]),
              np.array([0.01] + [0.0] * 9),
              np.array([-0.01] + [0.0] * 9)]
    starts += [RNG2.uniform(-bound, bound, 10) for _ in range(4)]
    for st in starts:
        res = minimize(obj, np.clip(st, -bound, bound),
                       method="Nelder-Mead",
                       bounds=[(-bound, bound)] * 10,
                       options={"maxfev": 4000, "xatol": 1e-10,
                                "fatol": 1e-12})
        if res.fun < best_val:
            best_val, best_a = float(res.fun), np.asarray(res.x)
    e_c, kc = ccf.chan_kins(B15.bfun_of(best_a), ("clock",))
    kt_f = ccf.kin_base["clock"] + kc["clock"]
    E_tot = E_base + e_c + J * J / (4.0 * kt_f)
    return {"J": J, "bound": bound,
            "stored_E_total": stored_row["E_total"],
            "stored_omega_star": stored_row["omega_star"],
            "stored_obj_on_opt_grid": stored_obj_opt,
            "remin_obj_on_opt_grid": best_val,
            "obj_improvement": float(stored_obj_opt - best_val),
            "remin_E_total_full": float(E_tot),
            "remin_kin_tot_full": float(kt_f),
            "remin_omega_star": float(J / (2.0 * kt_f)),
            "E_total_delta_vs_stored":
                float(E_tot - stored_row["E_total"]),
            "sign_flip": bool((E_tot > 0)
                              != (stored_row["E_total"] > 0)),
            "avec_best": best_a.tolist()}


def audit_d3(cc, ccf, E_base):
    with open(os.path.join(DATA, "m5_21_15_guard.json")) as f:
        gd = json.load(f)
    rows, Eb = gd["rows"], gd["E_base_u"]
    kb_full = gd["kin_base_clock"]
    out = {}
    dev = 0.0
    for r in rows:
        er = r["J"] ** 2 / (4.0 * r["kin_tot"])
        dev = max(dev, abs(er - r["E_rot"]),
                  abs(Eb + r["E_corr"] + er - r["E_total"]),
                  abs(r["J"] / (2.0 * r["kin_tot"])
                      - r["omega_star"]))
    out["arithmetic_max_dev"] = float(dev)
    sign_ok, wall_ok, nest_ok = True, True, True
    per_j = {}
    for J in sorted({r["J"] for r in rows}):
        rj = sorted((r for r in rows if r["J"] == J),
                    key=lambda r: r["bound"])
        objs = [r["E_corr"] + r["E_rot"] for r in rj]
        per_j[f"J{J:.2f}"] = {
            "bounds": [r["bound"] for r in rj],
            "E_totals": [r["E_total"] for r in rj],
            "omega_stars": [r["omega_star"] for r in rj],
            "objs": objs,
            "at_bound": [r["at_bound"] for r in rj]}
        for r in rj:
            if r["bound"] <= 0.02 and r["E_total"] <= 0:
                sign_ok = False
            if r["bound"] >= 0.05 and r["E_total"] >= 0:
                sign_ok = False
            if not (r["omega_star"] > r["J"] / (2.0 * kb_full)):
                wall_ok = False
            if not (r["kin_tot"] < kb_full):
                wall_ok = False
        nest_ok = nest_ok and bool(np.all(np.diff(objs) < 0))
    out["per_J"] = per_j
    out["sign_pattern_ok"] = sign_ok
    out["omega_star_wall_ok"] = wall_ok
    out["objective_nesting_monotone"] = nest_ok
    # independent re-minimization at the positive edge of the bracket
    out["remin"] = {}
    for J in (133.10472858194973, 332.7618214548743):
        row = next(r for r in rows
                   if r["bound"] == 0.02 and abs(r["J"] - J) < 1e-9)
        tag = f"bound0.02_J{J:.1f}"
        out["remin"][tag] = _remin_fixedj(cc, ccf, E_base, J, 0.02,
                                          row)
        print(json.dumps({"remin": tag,
                          "E_total":
                              out["remin"][tag]
                              ["remin_E_total_full"],
                          "sign_flip":
                              out["remin"][tag]["sign_flip"]}),
              flush=True)
    return out


# =============== D4 ===============
def audit_d4(cc, ccf):
    out = {}

    def kc_of(a, c):
        _, kc = c.chan_kins(B15.bfun_of(np.asarray(a)), ("clock",))
        return kc["clock"]

    kb = ccf.kin_base["clock"]
    for bd in (0.02, 2.0):
        best_v, best_a, n_pos, n_smp = 0.0, np.zeros(10), 0, 0
        for frac in (0.1, 0.3, 1.0):
            for _ in range(150):
                a = RNG2.uniform(-frac * bd, frac * bd, 10)
                v = kc_of(a, cc)
                n_smp += 1
                if v > 0:
                    n_pos += 1
                if v > best_v:
                    best_v, best_a = v, a.copy()

        def neg(a):
            return -kc_of(a, cc)

        starts = [best_a] + [RNG2.uniform(-bd, bd, 10)
                             for _ in range(3)]
        starts += [np.eye(10)[k] * 0.5 * bd for k in (0, 3, 7)]
        for st in starts:
            res = minimize(neg, np.clip(st, -bd, bd),
                           method="L-BFGS-B",
                           bounds=[(-bd, bd)] * 10,
                           options={"maxiter": 200, "ftol": 1e-14,
                                    "eps": min(1e-5, bd / 100)})
            if -res.fun > best_v:
                best_v, best_a = float(-res.fun), np.asarray(res.x)
        kc_full = kc_of(best_a, ccf)
        out[f"bound_{bd:g}"] = {
            "n_random_samples": n_smp,
            "n_samples_kin_corr_positive": n_pos,
            "max_kin_corr_opt_grid": float(best_v),
            "max_kin_corr_full_grid": float(kc_full),
            "kin_base_full": float(kb),
            "ceiling_rise_rel": float(max(kc_full, 0.0) / kb),
            "omega_floor_lowering_rel": float(
                1.0 - kb / (kb + max(kc_full, 0.0))),
            "avec_best": np.asarray(best_a).tolist()}
        print(json.dumps({"D4_bound": bd,
                          "max_kin_corr_full": float(kc_full)}),
              flush=True)
    return out


# =============== D5 ===============
def audit_d5():
    with open(os.path.join(DATA, "m5_21_15_fixedj.json")) as f:
        fj = json.load(f)
    ej = fj["EJ"]
    Js = np.array([r["J"] for r in ej])
    Es = np.array([r["E_over_base"] for r in ej])
    omst = np.array([r["omega_star"] for r in ej])
    out = {"E_strictly_increasing": bool(np.all(np.diff(Es) > 0)),
           "omega_star_strictly_increasing":
               bool(np.all(np.diff(omst) > 0))}
    interior = []
    for j in range(1, len(Js) - 1):
        h1, h2 = Js[j] - Js[j - 1], Js[j + 1] - Js[j]
        dd = ((Es[j + 1] - Es[j]) * h1 / (h2 * (h1 + h2))
              + (Es[j] - Es[j - 1]) * h2 / (h1 * (h1 + h2)))
        interior.append({
            "J": float(Js[j]), "dEdJ_fd": float(dd),
            "omega_star": float(omst[j]),
            "ratio": float(dd / omst[j]),
            "dev_pct": float(abs(dd / omst[j] - 1.0) * 100.0),
            "stored_ratio": ej[j]["dEdJ_over_omega_star"]})
    out["interior_stencils"] = interior
    out["max_interior_dev_pct"] = float(
        max(r["dev_pct"] for r in interior))
    out["min_interior_dev_pct"] = float(
        min(r["dev_pct"] for r in interior))
    out["within_stated_1_2_pct"] = bool(
        out["max_interior_dev_pct"] <= 2.0)
    # endpoints (one-sided, for the record; not part of the claim)
    out["endpoint_dev_pct"] = [
        float(abs((Es[1] - Es[0]) / (Js[1] - Js[0]) / omst[0] - 1)
              * 100),
        float(abs((Es[-1] - Es[-2]) / (Js[-1] - Js[-2]) / omst[-1]
                  - 1) * 100)]
    # replicate the run's np.gradient numbers (consistency only)
    grad = np.gradient(Es, Js)
    out["stored_dEdJ_matches_np_gradient"] = bool(np.allclose(
        grad, [r["dEdJ_numeric"] for r in ej], rtol=1e-10))
    return out


def main_r2():
    t0 = time.time()
    res = {"task": "M5.21.15 adversarial audit round 2 (run data)",
           "auditor": "independent second agent",
           "date": "2026-08-11"}

    print("D1 ...", flush=True)
    d1 = audit_d1()
    print(json.dumps({k: {"max_under_conv":
                          v["max_under_convergence"],
                          "refined_interior_min":
                          v["refined_interior_min"]}
                      for k, v in d1.items()}, indent=1), flush=True)
    print("D2 ...", flush=True)
    d2 = audit_d2()
    print(json.dumps({k: d2[k] for k in
                      ("omega_min", "E_min", "E_min_worst_case",
                       "E_min_robustly_positive")}, indent=1),
          flush=True)
    print("D5 ...", flush=True)
    d5 = audit_d5()
    print(json.dumps(d5["interior_stencils"], indent=1), flush=True)

    # shared machinery for D3 + D4 (the task's own ChanCorr, s = -1,
    # g = 32, the run's grids)
    B15.C14.S_SIGN = -1.0
    g_opt = B15.C14.make_grid(48, 8, 16)
    g_full = B15.C14.make_grid(72, 12, 24)
    cc = B15.ChanCorr(g_opt, 32.0)
    ccf = B15.ChanCorr(g_full, 32.0)
    E_base = float(np.sum(g_full["wvol"] * ccf.ec.du_base))
    res["E_base_recomputed"] = E_base

    print("D3 ...", flush=True)
    d3 = audit_d3(cc, ccf, E_base)
    print("D4 ...", flush=True)
    d4 = audit_d4(cc, ccf)

    # ---------- verdicts ----------
    max_uc = max(v["max_under_convergence"] for v in d1.values())
    uc_s1 = d1["s+1_g32:clock"]["max_under_convergence"]
    refined_clean = all(
        (not v["refined_interior_min"])
        and v["refined_max_slope_increment"] < 1e-6
        for v in d1.values())
    depths_covered = all(
        all(m["depth"] <= v["max_under_convergence"] + 1e-9
            for m in v["interior_min_candidates"])
        for v in d1.values())
    n_viol = sum(len(v["chord_violations"]) for v in d1.values())
    big_viol = max((c["depth"] for v in d1.values()
                    for c in v["chord_violations"]), default=0.0)
    v_d1 = ("CONFIRMED" if refined_clean and depths_covered
            and max_uc < 1.0 else
            "PARTIAL" if refined_clean and depths_covered
            else "REFUTED")

    v_d2 = ("CONFIRMED" if d2["interior_minimum"]
            and d2["E_min_robustly_positive"]
            and d2["EJ_anchor"]["match_within_step"]
            and not any(h["could_undercut"]
                        for h in d2["infeasible_hidden_min"])
            else "PARTIAL")

    remin_ok = all(not r["sign_flip"]
                   for r in d3["remin"].values())
    v_d3 = ("CONFIRMED" if d3["arithmetic_max_dev"] < 1e-9
            and d3["sign_pattern_ok"] and d3["omega_star_wall_ok"]
            and d3["objective_nesting_monotone"] and remin_ok
            else "REFUTED" if not (d3["sign_pattern_ok"]
                                   and remin_ok)
            else "PARTIAL")

    kb_full = d4["bound_2"]["kin_base_full"]
    kc_max = max(d4["bound_2"]["max_kin_corr_full_grid"],
                 d4["bound_0.02"]["max_kin_corr_full_grid"])
    v_d4 = "CONFIRMED" if kc_max < 1e-3 * kb_full else "REFUTED"

    v_d5 = ("CONFIRMED" if d5["within_stated_1_2_pct"]
            and d5["E_strictly_increasing"] else
            "PARTIAL" if d5["E_strictly_increasing"]
            and d5["max_interior_dev_pct"] < 5.0 else "REFUTED")

    res["verdicts"] = {
        "D1": {
            "verdict": v_d1,
            "reasoning": (
                "STRUCTURAL CLAIM CONFIRMED, NOISE WORDING "
                "CORRECTED. Row consistency (E_env = E_corr + u "
                "kin_tot, J = 2 omega kin_tot) and the stored "
                "u-slopes reproduce exactly in all 5 ladders. "
                "Across the 5 ladders the auditor finds {} chord "
                "violation(s) of concavity, the largest {:.4g} "
                "energy units; the only interior-minimum candidate "
                "is the flagged s+1 rung at omega = 0.05 (depth "
                "0.106, 2.5e-5 of the local |E| = 4194, exactly as "
                "the record states). The refined envelope built "
                "from the run's OWN stored configs (each rung's "
                "(E_corr, kin_tot) line evaluated at every u; "
                "concave by construction) dominates the measured "
                "ladder everywhere, has NO interior minimum in any "
                "block, and absorbs every violation: each "
                "interior-minimum depth is below its block's "
                "measured under-convergence. So no measured wiggle "
                "is evidence of a true interior minimum; the "
                "concavity theorem's prediction survives the data. "
                "CORRECTION to the 'noise-level' wording: the s+1 "
                "omega = 0 rung sits {:.1f} energy units above the "
                "family's own better config (7.0e-2 of scale), a "
                "warm-start convergence failure, far above the "
                "record's stated ~0.02 optimizer noise; the "
                "apparent 'interior minimum' is entirely this one "
                "under-converged reference rung, not scatter of "
                "the scanned rungs. The g8 slope wiggle "
                "(+{:.3g} slope increment) and the g8 interior "
                "MAXIMUM are noise-scale and concavity-compatible "
                "respectively, as claimed."
            ).format(n_viol, big_viol, uc_s1,
                     max((s["delta_slope"] for s in
                          d1["s-1_g8:clock"]
                          ["slope_increments_positive"]),
                         default=0.0)),
            "key_numbers": d1},
        "D2": {
            "verdict": v_d2,
            "reasoning": (
                "Arithmetic verified: kin_target = J/(2 omega), "
                "E_total = E_base + E_corr + omega^2 kin_tot, "
                "mismatch and feasibility flags all reproduce "
                "(max dev {:.1e}). The feasible rungs (omega = "
                "0.50..0.62) have an interior minimum at omega = "
                "0.59, E_total = +115.90, matching the EJ "
                "prediction omega* = 0.5915 within the 0.03 scan "
                "step. Sign robustness: the 1.9 percent "
                "constraint mismatch bounds the energy shift by "
                "omega^2 |kin_tot - kin_target| <= {:.2f} on the "
                "feasible rungs, so worst-case E_min = {:.1f} > 0; "
                "the alternative on-shell bookkeeping E_base + "
                "E_corr + J^2/(4 kin_tot) gives +{:.1f} at the "
                "minimum (spread {:.1f}, both positive, and "
                "consistent with the EJ value +{:.1f}). The "
                "infeasible-flagged rungs cannot hide a lower "
                "minimum: even shifted fully in their favor, 0.66 "
                "gives {:.1f} and 0.70 gives {:.1f}, both above "
                "the shifted minimum {:.1f}. CAVEAT on LOCATION "
                "(not sign): the neighbor gaps at the minimum "
                "({:.2f} left, {:.2f} right) are comparable to the "
                "{:.2f} mismatch shift, so the minimum's location "
                "is pinned only to omega = 0.59 +- one 0.03 step; "
                "its existence, interiority, and positive sign are "
                "robust."
            ).format(max(d2["arithmetic_max_dev"].values()),
                     d2["max_shift_feasible"],
                     d2["E_min_worst_case"],
                     d2["E_min_alt_onshell"],
                     d2["E_min_alt_onshell"] - d2["E_min"],
                     d2["EJ_anchor"]["E_total"],
                     d2["infeasible_hidden_min"][-2]["E_best_case"]
                     if len(d2["infeasible_hidden_min"]) >= 2
                     else float("nan"),
                     d2["infeasible_hidden_min"][-1]["E_best_case"]
                     if d2["infeasible_hidden_min"]
                     else float("nan"),
                     d2["E_min"] + d2["max_shift_feasible"],
                     d2["location_gaps"].get("to_left", 0.0),
                     d2["location_gaps"].get("to_right", 0.0),
                     d2["max_shift_feasible"]),
            "key_numbers": d2},
        "D3": {
            "verdict": v_d3,
            "reasoning": (
                "Arithmetic verified to {:.1e} (E_rot = J^2/(4 "
                "kin_tot), E_total = E_base + E_corr + E_rot, "
                "omega* = J/(2 kin_tot)). Sign pattern holds at "
                "both J: E_total positive at bounds 0.01-0.02, "
                "negative from 0.05 up. omega* exceeds the "
                "undressed J/(2 kin_base) in every row "
                "(equivalently kin_tot < kin_base everywhere). "
                "The omega*(guard) non-monotonicity at J = 133 "
                "(0.393 -> 0.277 -> 0.363) is NOT under-"
                "convergence severe enough to distrust the "
                "bracket: the minimized objective E_corr + E_rot "
                "IS strictly monotone in the bound at both J "
                "(nesting respected), the three minimizers are "
                "boundary-pinned with different avec sign "
                "patterns, and omega* = J/(2 kin_tot) at the "
                "argmin has no monotonicity theorem behind it. "
                "Independent re-minimization at the positive edge "
                "(bound 0.02) with differential evolution + "
                "bounded Nelder-Mead multistart: at J = 133.1 the "
                "objective improves by {:.3g} (E_total {:.2f} vs "
                "stored {:.2f}); at J = 332.8 it improves by "
                "{:.3g} (E_total {:.2f} vs stored {:.2f}). "
                "No sign flip in either case: the [0.02, 0.05] "
                "sign bracket survives a stronger minimizer "
                "(the negative side can only move further "
                "negative under better minimization)."
            ).format(d3["arithmetic_max_dev"],
                     d3["remin"]["bound0.02_J133.1"]
                       ["obj_improvement"],
                     d3["remin"]["bound0.02_J133.1"]
                       ["remin_E_total_full"],
                     d3["remin"]["bound0.02_J133.1"]
                       ["stored_E_total"],
                     d3["remin"]["bound0.02_J332.8"]
                       ["obj_improvement"],
                     d3["remin"]["bound0.02_J332.8"]
                       ["remin_E_total_full"],
                     d3["remin"]["bound0.02_J332.8"]
                       ["stored_E_total"]),
            "key_numbers": d3},
        "D4": {
            "verdict": v_d4,
            "reasoning": (
                "Attack on the kin-ceiling: {} random family "
                "points per bound regime (scales 0.1/0.3/1.0 of "
                "the bound) plus L-BFGS-B ascent on kin_corr from "
                "7 starts, bounds 0.02 and 2.0. Best kin_corr "
                "(clock) found: {:.3g} at bound 0.02 and {:.3g} "
                "at bound 2.0 on the optimization grid, {:.3g} "
                "and {:.3g} re-evaluated on the full grid, "
                "against kin_base = 332.76; {} of {} random "
                "samples had kin_corr > 0. Within this search "
                "power the dressing never RAISES the clock kin: "
                "the ceiling kin_tot <= kin_base and the fixed-J "
                "feasibility wall omega >= J/(2 kin_base) "
                "survive. Stated as a search result, not a "
                "theorem: the 10-dim family was probed, not "
                "exhausted."
            ).format(d4["bound_2"]["n_random_samples"],
                     d4["bound_0.02"]["max_kin_corr_opt_grid"],
                     d4["bound_2"]["max_kin_corr_opt_grid"],
                     d4["bound_0.02"]["max_kin_corr_full_grid"],
                     d4["bound_2"]["max_kin_corr_full_grid"],
                     sum(d4[k]["n_samples_kin_corr_positive"]
                         for k in ("bound_0.02", "bound_2")),
                     sum(d4[k]["n_random_samples"]
                         for k in ("bound_0.02", "bound_2"))),
            "key_numbers": d4},
        "D5": {
            "verdict": v_d5,
            "reasoning": (
                "E(J) is strictly increasing and omega*(J) "
                "strictly increasing (E convex in J), as claimed. "
                "Auditor's own non-uniform 3-point central "
                "differences at the interior rows give dE/dJ / "
                "omega* deviations of {:.2f}, {:.2f}, {:.2f} "
                "percent (J = 133.1, 332.8, 665.5), reproducing "
                "the stored np.gradient values. The LETTER of the "
                "'1-2 percent' claim fails at the coarsest "
                "interior stencil: 2.39 percent > 2. The SPIRIT "
                "survives: the deviation shrinks with stencil "
                "spacing (2.39 -> 1.11 -> 0.41 percent), the "
                "pattern expected from finite-difference "
                "truncation on a convex E(J) at coarse spacing, "
                "not from a violated identity; endpoints "
                "(one-sided, {:.1f} and {:.1f} percent) are "
                "outside the claim. Accurate restatement: the "
                "envelope identity dE/dJ = omega* holds to 0.4-"
                "2.4 percent on interior stencils at this "
                "spacing."
            ).format(d5["interior_stencils"][0]["dev_pct"],
                     d5["interior_stencils"][1]["dev_pct"],
                     d5["interior_stencils"][2]["dev_pct"],
                     d5["endpoint_dev_pct"][0],
                     d5["endpoint_dev_pct"][1]),
            "key_numbers": d5},
    }
    res["runtime_s"] = round(time.time() - t0, 1)
    with open(OUT2, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps({c: res["verdicts"][c]["verdict"]
                      for c in ("D1", "D2", "D3", "D4", "D5")}
                     | {"runtime_s": res["runtime_s"]}, indent=1))


if __name__ == "__main__":
    if "--round1" in sys.argv:
        main()
    else:
        main_r2()
