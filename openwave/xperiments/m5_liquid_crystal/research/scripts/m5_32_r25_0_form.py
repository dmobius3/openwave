"""M5.32 R25-0: the form level of the strand rung, on the author's 2026-09-22 paper
post (v2.8), the 2026-09-22 19:24 UTC reply (v2.9) and our plan-time checks.

N = M eta, eta = diag(-1, 1, 1, 1), the certified biaxial vacuum diag(8, 1, delta, 0)
(the N-spectrum (-8, 1, delta, 0)), g = 8, delta = 0.3 unless said.
    E = 4 sum_{i<j} <F_ij, F_ij>_eta + V4,  F_ij = [d_i M, d_j M]_eta = A_i eta A_j - A_j eta A_i,
    <F, G>_eta = tr(eta F eta G^T),  V4 = w sum_{p<=4} (tr N^p - C_p)^2.
The paper's normalization (u = 4 sum abs([d_i M, d_j M])^2, W sum_p (tr M^p - C_p)^2) is this one
on the spatial block, so the paper's W is our w.

Checks (each carries a `fails_if` line in the JSON):

a   the block sector, exact in sympy: with only the transverse pair moving, S = s_0 I +
    b(rho) (cos(m phi) sigma_z + sin(m phi) sigma_x), the curvature density is
    u = 32 m^2 (b b' / rho)^2 (the 4x4 eta commutator and the 2x2 one agree: the winding costs
    nothing, only the core does), the trace differences are tr N^p - C_p = 0, 2 (f - f_0),
    6 s_0 (f - f_0), (f - f_0)(12 s_0^2 + 2 (f + f_0)) with f = b^2 (so E_u scales as delta^4
    exactly while V4 carries K(s_0) = 4 + 36 s_0^2 + 144 s_0^4 at leading order), and the
    Bogomolny bound: the 1D integrand minus 2 sqrt(8 w K) (f_0 - f) f' is a perfect square, the
    cross term integrates to sqrt(8 w K) f_0^2, so T >= pi sqrt(32 w K) b_0^4 m, attained by
    f = f_0 (1 - exp(-kappa rho^2)), kappa = sqrt(w K / 32) at m = 1 (the reply's result).
a2  the bound numerically, our own midpoint integrator with the EXACT block potential and the
    axis pinned at f(0) = 0: the 1D minimizer sits between the leading-K bound and the exact-
    vacuum-K bound at every (delta, w) and equals the bound to 1e-4 at delta 0.03; the trap:
    a minimizer with the axis value free lands below the bound by lifting f(0) (a singular
    core); the stack's own energy of the BPS profile on the slab against the 1D value at three
    spacings.
b   the exterior density: the pure-winding field (b = b_0 everywhere) has zero continuum
    density (a one-coordinate texture has F = 0); on the lattice it leaves a residue that falls
    with h and with rho; the BPS field's density beyond three core radii is exponentially small;
    the m = 2 BPS line costs twice the m = 1 line (the additivity in the index).
c   (mode stack) the E . B witness in the AUTHOR's definition (eb_outgoing.py, read, not run):
    C_mu nu = [d_mu M, d_nu M] with the clock jet d_t M = a_0 (the stack's generator catalog:
    `rot_z` the rigid rotation about z, `clock_local` about the local director), F_mu nu = the
    pair component of C in the local eigenframe (the rotation about the director), E_i = F_0i,
    B = (F_23, F_31, F_12). Exact: at fixed eigenvalues every jet is [Omega, M] and F is built
    from two one-forms, so E . B = 0 to round-off on random data (the reply's statement); an
    eigenvalue gradient in one jet switches it on. THE TRAP, documented: the full-trace four-form
    sum <[A_mu, A_nu], [A_rho, A_sigma]>_eta over the cyclic pairings vanishes for ANY symmetric
    jets by the Jacobi identity (tr([A,B][C,D]) = tr(A [B,[C,D]])), boosts included: the E . B
    question cannot be put to the certified action's own curvatures F = [dM, dM]_eta at all, it
    lives on the frame's projected component. On the lattice, the fixed-eigenvalue hedgehog's
    ratio must fall with h. Then the stored melted-core electrons (R21 `S1` W1 x 25 n 32 on the
    biaxial exterior; R22-1 `rad_pin` n 48 on the uniaxial one, halo cells only, where the pair
    frame exists) with both clocks, reported against the reply's 0.05 to 0.24 and its dipole
    parity of B_r, and the R24 boost ripple laid on the R21 field.
d   the twisting-vacuum pencil, exact in sympy at a point: about M(z) = R(tau z) D R(tau z)^T
    (the pair plane rotating about the director as z advances, d_z M = tau C, C = [G, D]) the
    quadratic part of the quartic Lagrangian for Phi cos(k . x - omega t) is
    4 tau^2 (omega^2 - k_x^2 - k_y^2) <[Phi, C]_eta, [Phi, C]_eta>_eta sin^2(...): no k_z, no
    tau beyond tau^2, omega = k_perp on every mode not commuting with C (the reply's L_2 up to
    normalization); the sign of the bracket for a spatial (rotation-type) and a time-space
    (boost-type) Phi is reported (the R24 witness lives here); a lattice second variation
    confirms the density on a plane wave to 1e-2; V4's own Hessian is a separate additive piece.
e   the spin gate, exact: on E_J = E_stat + J^2 / (4 C) with omega = dE/dJ the condition
    2 J omega / E = 1 gives J_*^2 = 4 C E_stat / 3, omega_* = sqrt(E_stat / (3 C)),
    E_rot / E = 1/4; the plan's sqrt(C E_stat) was the E_rot << E limit.
f   the frame-index bookkeeping of the R22-2 pair textures (the plan-time finding as a check):
    the pulled-back frame e_1 = normalize(d n / d Re zeta) is orthogonal to n, its singular set
    follows n = -z (the lower half-axis and a region below the pair, not the segment), and the
    line-field windings on axial loops (read only where abs(n_z) > 0.95) give the segment and
    the far half-axes their indices for the +- and the like pair.

Modes: run (a, a2, b, d, e, f; about 3 min) | stack (c; about 2 min, reads the stored R22-1 and
R21 end fields, SKIPPED if absent). Output data/m5_32_r25_0_form.json.
"""

import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r25_0_form.json")
G_ = 8.0
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
R0 = _load("m5_32_r20_0_class", "m5_32_r20_0_class.py")
S1 = _load("m5_32_r25_1_strand", "m5_32_r25_1_strand.py")
W1 = B3.W1
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ================= a: the block sector in sympy =================
def check_a():
    import sympy as sp

    s0, b0, rho, phi, w, K, kap = sp.symbols("s0 b0 rho phi w K kappa", positive=True)
    m = sp.symbols("m", positive=True, integer=True)
    sg = sp.symbols("sg", real=True)
    b = sp.Function("b")(rho)
    sz, sx = sp.Matrix([[1, 0], [0, -1]]), sp.Matrix([[0, 1], [1, 0]])
    S = s0 * sp.eye(2) + b * (sp.cos(m * phi) * sz + sp.sin(m * phi) * sx)
    d_rho = lambda X: X.diff(rho)  # noqa: E731
    d_phi = lambda X: X.diff(phi)  # noqa: E731
    dx = lambda X: sp.cos(phi) * d_rho(X) - sp.sin(phi) / rho * d_phi(X)  # noqa: E731
    dy = lambda X: sp.sin(phi) * d_rho(X) + sp.cos(phi) / rho * d_phi(X)  # noqa: E731
    Ax, Ay = dx(S), dy(S)
    F2 = Ax * Ay - Ay * Ax
    u2 = 4 * (F2 * F2.T).trace()
    target = 32 * m**2 * (b * b.diff(rho) / rho) ** 2
    ok_u2 = sp.simplify(sp.trigsimp(sp.expand(u2 - target))) == 0
    # the 4x4 embedding with the eta commutator and the eta inner product
    eta = sp.diag(-1, 1, 1, 1)
    M4 = sp.zeros(4, 4)
    M4[0, 0] = -sg
    M4[3, 3] = 1
    M4[1:3, 1:3] = S
    A4x, A4y = dx(M4), dy(M4)
    F4 = A4x * eta * A4y - A4y * eta * A4x
    u4 = 4 * (eta * F4 * eta * F4.T).trace()
    ok_u4 = sp.simplify(sp.trigsimp(sp.expand(u4 - target))) == 0
    # the trace differences on the block
    f, f0 = sp.symbols("f f0", positive=True)
    N4 = M4 * eta
    diffs = []
    for p in range(1, 5):
        trp = (N4**p).trace()
        cp = sg**p + 1 + (s0 + b0) ** p + (s0 - b0) ** p
        d = sp.expand(sp.trigsimp(sp.expand(trp - cp)))
        diffs.append(sp.simplify(d.subs({b: sp.sqrt(f), b0: sp.sqrt(f0)})))
    expected = [0, 2 * (f - f0), 6 * s0 * (f - f0), (f - f0) * (12 * s0**2 + 2 * (f + f0))]
    ok_tr = all(sp.simplify(a - e) == 0 for a, e in zip(diffs, expected))
    # the Bogomolny bound in 1D: the integrand minus the cross term is a perfect square
    fr = sp.Function("f")(rho)
    integrand = 8 * fr.diff(rho) ** 2 / rho + w * K * rho * (f0 - fr) ** 2
    cross = 2 * sp.sqrt(8 * w * K) * (f0 - fr) * fr.diff(rho)
    square = (sp.sqrt(8 / rho) * fr.diff(rho) - sp.sqrt(w * K * rho) * (f0 - fr)) ** 2
    ok_sq = sp.simplify(sp.expand(integrand - cross - square)) == 0
    ff = sp.symbols("ff", positive=True)
    ok_int = (
        sp.simplify(
            sp.integrate(2 * sp.sqrt(8 * w * K) * (f0 - ff), (ff, 0, f0))
            - sp.sqrt(8 * w * K) * f0**2
        )
        == 0
    )
    bound = 2 * sp.pi * sp.sqrt(8 * w * K) * f0**2
    ok_bound = (
        sp.simplify(bound - sp.pi * sp.sqrt(32 * w * K) * b0**4) == 0
        if False
        else sp.simplify(bound.subs(f0, b0**2) - sp.pi * sp.sqrt(32 * w * K) * b0**4) == 0
    )
    # the BPS profile solves the first-order equation with kappa = sqrt(w K / 32)
    fb = f0 * (1 - sp.exp(-kap * rho**2))
    eq = sp.sqrt(8 / rho) * fb.diff(rho) - sp.sqrt(w * K * rho) * (f0 - fb)
    ok_bps = sp.simplify(eq.subs(kap, sp.sqrt(w * K / 32))) == 0
    # the m scaling of the bound: u has m^2, the square root gives m
    out = {
        "u_2x2_is_32m2_bbprime_over_rho_sq": bool(ok_u2),
        "u_4x4_eta_equals_2x2": bool(ok_u4),
        "trace_differences": [str(d) for d in diffs],
        "trace_differences_ok": bool(ok_tr),
        "K_leading": "4 + 36 s0^2 + 144 s0^4 (the quartic (12 s0^2 + 2 (f + f0))^2 at f = 0)",
        "perfect_square_ok": bool(ok_sq),
        "cross_term_integral_ok": bool(ok_int),
        "bound_pi_sqrt32wK_b0^4_ok": bool(ok_bound),
        "bps_profile_solves_first_order_ok": bool(ok_bps),
        "m_scaling": "T_m = m T_1 (u carries m^2, the bound its square root)",
    }
    out["PASS"] = all(out[k] for k in out if k.endswith("_ok") or k.startswith("u_"))
    out["fails_if"] = (
        "any symbolic identity fails: the 2x2 or 4x4 density is not 32 m^2 (b b'/rho)^2, a trace difference "
        "departs from the stated polynomial, the integrand minus the cross term is not a perfect square, "
        "the cross term does not integrate to sqrt(8 w K) f0^2, or the BPS profile fails the first-order equation"
    )
    return out


# ================= a2: the bound numerically, our own integrator =================
def check_a2():
    from scipy.optimize import minimize

    out = {"rows": []}
    ok = True
    for delta in (0.3, 0.1, 0.03):
        for w1s in (6.25, 25.0, 100.0):
            w = W1 * w1s
            s0 = b0 = delta / 2.0
            f0 = b0**2
            lead = S1.T_bps(delta, w)
            exv = float(
                np.pi * np.sqrt(32 * w * (4 + 36 * s0**2 + (12 * s0**2 + 4 * f0) ** 2)) * b0**4
            )
            t1 = S1.T_1d_min(delta, w)
            # the trap: the axis value free
            rc = np.linspace(1e-3, 40.0, 801)
            f_init = f0 * (1 - np.exp(-S1.kappa_bps(delta, w) * rc**2))
            r = minimize(
                lambda x: S1.T_1d_of(np.concatenate([x, [f0]]), rc, f0, s0, w),
                f_init[:-1],
                method="L-BFGS-B",
                options={"maxiter": 3000, "ftol": 1e-14},
            )
            row = {
                "delta": delta,
                "w1s": w1s,
                "bound_lead_K": lead,
                "bound_exact_vacuum_K": exv,
                "T_1d_axis_pinned": t1,
                "T_1d_axis_free": float(r.fun),
                "axis_value_over_f0_when_free": float(r.x[0] / f0),
                "T_over_delta4": t1 / delta**4,
            }
            row["ok"] = bool(lead * (1 - 1e-5) <= t1 <= exv * (1 + 1e-5)) and (
                delta > 0.03 or abs(t1 / lead - 1) < 1e-4
            )
            ok = ok and row["ok"]
            out["rows"].append(row)
    # the stack's energy of the BPS profile on the slab
    stack = []
    for delta in (0.3, 0.1):
        w = W1 * 25.0
        t1 = S1.T_1d_min(delta, w)
        for n, L in ((32, 48.0), (48, 48.0), (96, 48.0)):
            h = L / n
            cfg = S1.slab_cfg(n, L, delta)
            M = S1.bps_field(n, S1.NZ, h, delta, w)
            eu, ev = S1.energy_parts(M, cfg, w)
            T = (eu + ev) / (S1.NZ * h)
            stack.append({"delta": delta, "h": h, "T_stack": T, "T_1d": t1, "ratio": T / t1})
            if h <= 1.0:
                ok = ok and abs(T / t1 - 1) < 0.015
    out["stack_bps_profile"] = stack
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "the axis-pinned 1D minimizer leaves [bound(leading K), bound(exact vacuum K)] at any (delta, w), or departs "
        "from the bound by over 1e-4 at delta 0.03, or the stack's energy of the BPS profile departs from the 1D "
        "value by over 1.5 percent at h <= 1"
    )
    return out


# ================= b: the exterior density on the stack =================
def check_b():
    delta, w = 0.3, W1 * 25.0
    out = {}
    res = {}
    for n, h in ((48, 1.0), (96, 0.5)):
        cfg = S1.slab_cfg(n, n * h, delta)
        Mw = S1.bps_field(n, S1.NZ, h, delta, w, kappa=1e6)
        ew = S1.density_u(Mw, cfg)[:, :, 1] / h**3  # per volume
        _, _, rho, _ = S1.slab_coords(n, h)
        band = {
            f"[{lo},{hi})": float(np.max(ew[(rho >= lo) & (rho < hi)]))
            for lo, hi in ((1, 2), (2, 4), (4, 6), (6, 9), (9, 14))
        }
        res[f"h{h:g}"] = band
    out["pure_winding_residue_per_volume_by_rho_band"] = res
    r1, r5 = res["h1"]["[4,6)"], res["h0.5"]["[4,6)"]
    out["residue_h_ratio_band_4_6"] = r1 / max(r5, 1e-300)
    out["residue_order_in_h"] = float(np.log2(r1 / max(r5, 1e-300)))
    out["residue_rho_falloff_h1"] = float(
        np.log(res["h1"]["[2,4)"] / res["h1"]["[9,14)"]) / np.log(11.5 / 3.0)
    )
    n, h = 48, 1.0
    cfg = S1.slab_cfg(n, 48.0, delta)
    Mb = S1.bps_field(n, S1.NZ, h, delta, w)
    eb = S1.density_u(Mb, cfg)[:, :, 1]
    _, _, rho, _ = S1.slab_coords(n, h)
    rc = float(np.sqrt(np.log(2.0) / S1.kappa_bps(delta, w)))
    out["bps_core_radius"] = rc
    out["bps_density_beyond_3rc_over_max"] = float(np.max(eb[rho > 3 * rc]) / np.max(eb))
    eu1, ev1 = S1.energy_parts(Mb, cfg, w)
    M2 = S1.bps_field(n, S1.NZ, h, delta, w, m=2)
    eu2, ev2 = S1.energy_parts(M2, cfg, w)
    out["additivity_T_m2_over_T_m1"] = (eu2 + ev2) / (eu1 + ev1)
    out["PASS"] = bool(
        out["residue_h_ratio_band_4_6"] > 4
        and out["bps_density_beyond_3rc_over_max"] < 1e-3
        and abs(out["additivity_T_m2_over_T_m1"] - 2) < 0.05
    )
    out["fails_if"] = (
        "the pure-winding residue in the band rho 4 to 6 does not fall by over 4 between h 1 and h 0.5 (it would be a "
        "density, not a residue), or the BPS field's density beyond three core radii exceeds 1e-3 of its maximum, "
        "or the m = 2 line does not cost twice the m = 1 line within 5 percent"
    )
    return out


# ================= c: the four-form witness =================
def _jets(M, h):
    return [(np.roll(M, -1, ax) - np.roll(M, 1, ax)) / (2 * h) for ax in range(3)]


def _trace_four_form(a0, A):
    """the FULL-TRACE four-form of the jet commutators (the trap): sum over the cyclic pairings of
    <[A_mu, A_nu]_eta, [A_rho, A_sigma]_eta>_eta; it vanishes for ANY symmetric jets by the Jacobi identity.
    """
    ce, ie = B3.comm_eta, B3.inner_eta
    F0 = [ce(a0, A[i]) for i in range(3)]
    F12, F23, F31 = ce(A[0], A[1]), ce(A[1], A[2]), ce(A[2], A[0])
    W = ie(F0[0], F23) + ie(F0[1], F31) + ie(F0[2], F12)
    N = sum(np.sum(F * F, axis=(-2, -1)) for F in F0 + [F12, F23, F31])
    return W, N


def _eb_fields(M, a0, h):
    """the author's E and B (eb_outgoing.py, read not run): C_mu nu = [d_mu M, d_nu M], F_mu nu = the pair
    component of C in the local eigenframe, e_low^T C e_mid (the rotation about the director), with the
    clock jet as d_t M; E_i = F_0i, B = (F_23, F_31, F_12). Here with the stack's eta commutator (equal to
    the plain one on block-diagonal fields, different once M_0i is switched on)."""
    A = [a0] + _jets(M, h)
    lam, V = np.linalg.eigh(M[..., 1:, 1:])
    e1 = np.zeros(M.shape[:-2] + (4,))
    e2 = np.zeros(M.shape[:-2] + (4,))
    e1[..., 1:] = V[..., :, 0]
    e2[..., 1:] = V[..., :, 1]

    def comp(X, Y):
        return np.einsum("...a,...ab,...b->...", e1, B3.comm_eta(X, Y), e2)

    E = np.stack([comp(A[0], A[i]) for i in (1, 2, 3)], -1)
    B = np.stack([comp(A[2], A[3]), comp(A[3], A[1]), comp(A[1], A[2])], -1)
    return E, B, lam


def _eb_reads(M, cfg, a0, mask, shells=(3.0, 4.5, 6.0, 9.0, 12.0)):
    h, n = cfg["h"], cfg["n"]
    E, B, lam = _eb_fields(M, a0, h)
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rhat = np.stack([X, Y, Z], -1) / np.maximum(r, 1e-300)[..., None]
    eb = np.sum(E * B, -1)
    en, bn = np.linalg.norm(E, axis=-1), np.linalg.norm(B, axis=-1)
    Br = np.sum(B * rhat, -1)
    out = {
        "sum_abs_EB_over_sum_EnBn": float(
            np.sum(np.abs(eb[mask])) / max(np.sum((en * bn)[mask]), 1e-300)
        ),
        "sum_EB_over_half_sum_F2": float(
            np.sum(eb[mask]) / max(0.5 * np.sum((en**2 + bn**2)[mask]), 1e-300)
        ),
        "sum_abs_EB_over_half_sum_F2": float(
            np.sum(np.abs(eb[mask])) / max(0.5 * np.sum((en**2 + bn**2)[mask]), 1e-300)
        ),
        "max_abs_EB_over_max_EnBn": float(
            np.max(np.abs(eb[mask])) / max(np.max((en * bn)[mask]), 1e-300)
        ),
        "pair_gap_min_in_mask": float(np.min((lam[..., 1] - lam[..., 0])[mask])),
        "shells": {},
    }
    for R in shells:
        sh = (np.abs(r - R) < 0.75 * h) & mask
        if not sh.any():
            continue
        north, south = sh & (Z > 0), sh & (Z < 0)
        bn_, bs_ = float(np.mean(Br[north])), float(np.mean(Br[south]))
        out["shells"][f"{R:g}"] = {
            "EB_over_EnBn": float(np.sum(eb[sh]) / max(np.sum((en * bn)[sh]), 1e-300)),
            "EB_mean_north": float(np.mean(eb[north])),
            "EB_mean_south": float(np.mean(eb[south])),
            "Br_mean_north": bn_,
            "Br_mean_south": bs_,
            "Br_even_over_odd": abs(bn_ + bs_) / max(abs(bn_ - bs_), 1e-300),
            "net_flux_Br_over_sum_absBr": float(
                np.sum(Br[sh]) / max(np.sum(np.abs(Br[sh])), 1e-300)
            ),
        }
    return out


def _witness(M, cfg, mask, clocks=("rot_z", "clock_local")):
    cat = B3.gen_catalog(cfg, M)
    out = {}
    for nm in clocks:
        a0 = cat[nm]
        rec = _eb_reads(M, cfg, a0, mask)
        W, N = _trace_four_form(a0, _jets(M, cfg["h"]))
        rec["trace_four_form_max_ratio"] = float(
            np.max(np.abs(W[mask])) / max(np.max(N[mask]), 1e-300)
        )
        rec["a0_null"] = bool(np.all(a0 == 0))
        out[nm] = rec
    return out


def _mask_of(cfg, margin=3, r_min=0.0):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    return (r < 0.5 * cfg["L"] - margin * h) & (r > r_min)


def check_c():
    R21 = _load("m5_32_r21_1_runs", "m5_32_r21_1_runs.py")
    R20 = R21.R20
    out = {}
    rng = np.random.default_rng(7)

    def so3(v):
        O = np.zeros((4, 4))
        O[1:, 1:] = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return O

    # (0) the trap: the full-trace four-form vanishes for ANY symmetric jets (boost entries included)
    worst_trace = 0.0
    for _ in range(300):
        M = B3.sym4(rng.standard_normal((4, 4)))
        J = [B3.sym4(rng.standard_normal((4, 4))) for _k in range(4)]
        W, N = _trace_four_form(J[0][None, None, None], [Jk[None, None, None] for Jk in J[1:]])
        worst_trace = max(worst_trace, float(abs(W).max() / max(N.max(), 1e-300)))
    out["trace_four_form_any_symmetric_jets_max_ratio"] = worst_trace
    out["trace_four_form_reading"] = (
        "tr([A,B][C,D]) = tr(A [B,[C,D]]) and the cyclic sum is the Jacobi identity, so the full-trace four-form "
        "of commutator curvatures is identically zero for any jets: an unfalsifiable witness, not the E . B"
    )
    # (i) the author's projected witness at fixed eigenvalues: jets [Omega, M], E . B must vanish
    worst_fixed, worst_melt = 0.0, 0.0
    for _ in range(300):
        Q, _r = np.linalg.qr(rng.standard_normal((3, 3)))
        M = np.zeros((4, 4))
        M[0, 0] = 8.0
        M[1:, 1:] = Q @ np.diag(np.sort(rng.uniform(0.1, 2.0, 3))) @ Q.T
        J = [
            so3(rng.standard_normal(3)) @ M - M @ so3(rng.standard_normal(3)).T * 0
            for _k in range(4)
        ]
        J = [Jk - Jk.T if False else B3.sym4(Jk) for Jk in J]
        Om = [so3(rng.standard_normal(3)) for _k in range(4)]
        J = [O @ M - M @ O for O in Om]
        E, B, _ = (
            _eb_fields(M[None, None, None], J[0][None, None, None], 1.0)
            if False
            else (None, None, None)
        )
        # by hand: the projected component with these exact jets (no lattice)
        lam, V = np.linalg.eigh(M[1:, 1:])
        e1 = np.concatenate([[0.0], V[:, 0]])
        e2 = np.concatenate([[0.0], V[:, 1]])
        comp = lambda X, Y: float(e1 @ B3.comm_eta(X, Y) @ e2)  # noqa: E731
        Ev = np.array([comp(J[0], J[i]) for i in (1, 2, 3)])
        Bv = np.array([comp(J[2], J[3]), comp(J[3], J[1]), comp(J[1], J[2])])
        worst_fixed = max(
            worst_fixed, abs(Ev @ Bv) / max(np.linalg.norm(Ev) * np.linalg.norm(Bv), 1e-300)
        )
        # with an eigenvalue gradient added to one spatial jet (a melt): E . B is allowed to be nonzero
        Jm = list(J)
        Jm[1] = Jm[1] + Q4 if False else Jm[1] + np.diag([0.0, 0.3, -0.2, 0.1])
        Ev = np.array([comp(Jm[0], Jm[i]) for i in (1, 2, 3)])
        Bv = np.array([comp(Jm[2], Jm[3]), comp(Jm[3], Jm[1]), comp(Jm[1], Jm[2])])
        worst_melt = max(
            worst_melt, abs(Ev @ Bv) / max(np.linalg.norm(Ev) * np.linalg.norm(Bv), 1e-300)
        )
    out["projected_EB_fixed_eigenvalues_max_ratio"] = worst_fixed
    out["projected_EB_with_an_eigenvalue_gradient_max_ratio"] = worst_melt
    # (ii) the lattice fixed-eigenvalue hedgehog at two spacings, the biaxial exterior (the pair frame defined)
    lat = {}
    for n, L in ((32, 48.0), (48, 48.0)):
        cfg = R21.cfg_of(n, L, G_, 0.3)
        M = R20.seed_axes(cfg, (1.0, 0.3, 0.0), r_c=1e-6)
        lat[f"h{cfg['h']:g}"] = _witness(M, cfg, _mask_of(cfg, r_min=3.0))
    out["lattice_fixed_eigenvalue_hedgehog_biaxial"] = lat
    orders = {}
    for ck in ("rot_z", "clock_local"):
        a = lat["h1.5"][ck]["sum_abs_EB_over_sum_EnBn"]
        b = lat["h1"][ck]["sum_abs_EB_over_sum_EnBn"]
        orders[ck] = float(np.log(a / max(b, 1e-300)) / np.log(1.5)) if a > 0 and b > 0 else None
    out["lattice_ratio_order_in_h"] = orders
    # (iii) the stored melted-core electrons
    stored = {}
    f21 = os.path.join(DATA, "m5_32_r21_1", "S1_Bseed_g8_d0.3_w25_n32.npz")
    f22 = os.path.join(DATA, "m5_32_r22_1", "rad_pin_d0.3_w25_n48_L48.npz")
    M21 = None
    if os.path.exists(f21):
        M = np.load(f21)["M"]
        cfg = R21.cfg_of(32, 48.0, G_, 0.3)
        stored["r21_S1_Bseed_w25_n32_biaxial"] = _witness(M, cfg, _mask_of(cfg))
        M21 = (M, cfg)
    else:
        stored["r21_S1_Bseed_w25_n32_biaxial"] = "SKIPPED: the local field is absent"
    if os.path.exists(f22):
        M = np.load(f22)["M"]
        cfg = R21.cfg_of(48, 48.0, G_, 0.3)
        lam = np.linalg.eigvalsh(M[..., 1:, 1:])
        halo = (
            lam[..., 1] - lam[..., 0]
        ) > 1e-2  # the pair frame is defined only where the pair is split
        stored["r22_1_rad_pin_n48_uniaxial_halo_only"] = _witness(M, cfg, _mask_of(cfg) & halo)
        stored["r22_1_rad_pin_n48_uniaxial_halo_only"]["halo_cells"] = int(
            (_mask_of(cfg) & halo).sum()
        )
    else:
        stored["r22_1_rad_pin_n48_uniaxial_halo_only"] = "SKIPPED: the local field is absent"
    out["stored_melted_core_electrons"] = stored
    out["reply_range_for_comparison"] = (
        "0.05 to 0.24 of abs(F)^2 with a melted core and the clock (evidence)"
    )
    # (iv) the R24 boost ripple on the biaxial R21 field
    if M21 is not None:
        from scipy.linalg import expm

        M, cfg = M21
        n, h = cfg["n"], cfg["h"]
        X, Y, Z = B3.coords(n, h)
        r = np.sqrt(X * X + Y * Y + Z * Z)
        mask = _mask_of(cfg)
        K = np.zeros((4, 4))
        K[0, 2] = K[2, 0] = 1.0  # a boost along y
        th = 0.05 * np.sin(np.pi / (2 * h) * X) * np.exp(-((r / 12.0) ** 4)) * mask
        Lam = np.array([expm(t * K) for t in th.ravel()]).reshape(th.shape + (4, 4))
        Mr = Lam @ M @ np.swapaxes(Lam, -1, -2)
        out["boost_ripple_on_r21_S1"] = _witness(Mr, cfg, mask)
    else:
        out["boost_ripple_on_r21_S1"] = "SKIPPED"
    ok = (
        worst_trace < 1e-13
        and worst_fixed < 1e-10
        and worst_melt > 1e-3
        and all(v is not None and v > 1.0 for v in orders.values())
    )
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "the full-trace four-form exceeds 1e-13 on random symmetric jets (it is the Jacobi identity), the projected "
        "E . B on fixed-eigenvalue jets exceeds 1e-10, an eigenvalue gradient does NOT switch it on (under 1e-3), "
        "or the lattice fixed-eigenvalue hedgehog's ratio does not fall with h (order under 1) for either clock; "
        "the stored-field and boost-ripple readings are reported against the reply's range, not gated"
    )
    return out


# ================= d: the twisting-vacuum pencil =================
def check_d():
    import sympy as sp

    om, kx, ky, kz, tau, eps, dl = sp.symbols("omega k_x k_y k_z tau epsilon delta", real=True)
    eta = sp.diag(-1, 1, 1, 1)
    D = sp.diag(8, 1, dl, 0)
    Gr = sp.zeros(4, 4)
    Gr[2, 3], Gr[3, 2] = (
        -1,
        1,
    )  # the rotation of the pair plane (indices 2, 3) about the director (index 1)
    C = Gr * D - D * Gr
    rng = np.random.default_rng(3)

    def rand_sym(kind):
        X = sp.zeros(4, 4)
        vals = [sp.Rational(int(v), 7) for v in rng.integers(-6, 7, 16)]
        k = 0
        for i in range(4):
            for j in range(i, 4):
                v = vals[k]
                k += 1
                if kind == "spatial" and (i == 0 or j == 0):
                    v = 0
                if kind == "timespace" and not (i == 0 and j > 0):
                    v = 0
                X[i, j] = X[j, i] = v
        return X

    def comm(A, B):
        return A * eta * B - B * eta * A

    def inner(F, Gm):
        return (eta * F * eta * Gm.T).trace()

    out = {"modes": {}}
    ok = True
    for kind in ("spatial", "timespace", "general"):
        Phi = rand_sym(kind)
        # the jets at a point where sin(k.x - omega t) = 1: d_t Phi = -eps omega Phi ... signs are squared away
        At = -eps * (-om) * Phi
        Ax, Ay, Az = eps * kx * Phi, eps * ky * Phi, tau * C + eps * kz * Phi
        kin = 4 * sum(inner(comm(At, Ai), comm(At, Ai)) for Ai in (Ax, Ay, Az))
        pot = 4 * sum(
            inner(comm(Ai, Aj), comm(Ai, Aj)) for Ai, Aj in ((Ax, Ay), (Ax, Az), (Ay, Az))
        )
        L = sp.expand(kin - pot)
        L2 = L.coeff(eps, 2)
        L0 = L.subs(eps, 0)
        bracket = inner(comm(Phi, C), comm(Phi, C))
        pred = 4 * tau**2 * (om**2 - kx**2 - ky**2) * bracket
        diff = sp.simplify(L2 - pred)
        rec = {
            "L0_background": str(sp.simplify(L0)),
            "L2_minus_prediction": str(diff),
            "bracket_<[Phi,C],[Phi,C]>_eta": str(sp.simplify(bracket)),
            "kz_free": bool(sp.simplify(L2).free_symbols.isdisjoint({kz})),
            "linear_in_eps_absent": bool(sp.simplify(L.coeff(eps, 1)) == 0),
        }
        rec["ok"] = bool(diff == 0 and rec["kz_free"] and rec["linear_in_eps_absent"])
        ok = ok and rec["ok"]
        out["modes"][kind] = rec
    # the sign of the bracket per mode class at delta 0.3 on random draws
    signs = {}
    for kind in ("spatial", "timespace"):
        vals = []
        for _ in range(20):
            Phi = rand_sym(kind)
            vals.append(float(inner(comm(Phi, C), comm(Phi, C)).subs(dl, sp.Rational(3, 10))))
        signs[kind] = {"min": min(vals), "max": max(vals)}
    out["bracket_sign_by_mode_class_delta_0.3"] = signs
    # the lattice second variation of the quartic term on a plane wave about the twisted vacuum
    n, h, tau_v, eps_v = 16, 1.0, 0.05, 1e-3
    cfg = B3.base_cfg(s=-1.0, n=n, L=float(n), delta=0.3)
    X, Y, Z = B3.coords(n, h)
    th = tau_v * Z
    Rz = np.zeros(th.shape + (4, 4))
    Rz[..., 0, 0] = Rz[..., 1, 1] = 1.0
    Rz[..., 2, 2] = Rz[..., 3, 3] = np.cos(th)
    Rz[..., 2, 3], Rz[..., 3, 2] = -np.sin(th), np.sin(th)
    Dn = np.diag([8.0, 1.0, 0.3, 0.0])
    Mbg = Rz @ Dn @ np.swapaxes(Rz, -1, -2)
    Phi0 = np.array(
        [[0.0, 0.0, 0.0, 0.0], [0.0, 0.2, 0.3, -0.1], [0.0, 0.3, -0.4, 0.5], [0.0, -0.1, 0.5, 0.1]]
    )
    kv = 2 * np.pi / n * np.array([1.0, 1.0, 1.0])
    wave = np.cos(kv[0] * X + kv[1] * Y + kv[2] * Z)
    inner_mask = np.zeros((n, n, n), dtype=bool)
    inner_mask[2:-2, 2:-2, 2:-2] = True

    def u_of(M):
        A = _jets(M, h)
        F = [B3.comm_eta(A[i], A[j]) for i, j in ((0, 1), (0, 2), (1, 2))]
        return 4 * sum(B3.inner_eta(Fk, Fk) for Fk in F)

    Mp = Mbg + eps_v * Phi0 * wave[..., None, None]
    Mm = Mbg - eps_v * Phi0 * wave[..., None, None]
    u2_lat = (u_of(Mp) + u_of(Mm) - 2 * u_of(Mbg)) / (2 * eps_v**2)
    # the prediction per cell: 4 tau^2 k_perp^2 <[Phi0, C],[Phi0, C]>_eta sin^2 (the omega 0 static case)
    Cn = Gr_np = np.zeros((4, 4))
    Gr_np[2, 3], Gr_np[3, 2] = -1.0, 1.0
    Cn = Gr_np @ Dn - Dn @ Gr_np
    # the bracket is invariant under the z rotation of both Phi0 and C only if Phi0 rotates too; use the local frame
    Phi_loc = np.swapaxes(Rz, -1, -2) @ Phi0 @ Rz  # Phi0 pulled back to the vacuum frame at each z
    Cc = B3.comm_eta(Phi_loc, np.broadcast_to(Cn, Phi_loc.shape))
    br = B3.inner_eta(Cc, Cc)
    sin2 = np.sin(kv[0] * X + kv[1] * Y + kv[2] * Z) ** 2
    # the central-difference jets carry the lattice wavenumber sin(k h) / h, not k
    ke = np.sin(kv * h) / h
    u2_pred = 4 * tau_v**2 * (ke[0] ** 2 + ke[1] ** 2) * br * sin2
    num, den = float(np.sum(u2_lat[inner_mask])), float(np.sum(u2_pred[inner_mask]))
    out["lattice_second_variation"] = {
        "sum_lattice": num,
        "sum_prediction_with_lattice_k": den,
        "ratio": num / den,
        "continuum_k_factor": float((ke[0] ** 2 + ke[1] ** 2) / (kv[0] ** 2 + kv[1] ** 2)),
    }
    ok_lat = abs(num / den - 1) < 1e-2
    out["PASS"] = bool(ok and ok_lat)
    out["fails_if"] = (
        "the symbolic quadratic Lagrangian departs from 4 tau^2 (omega^2 - k_x^2 - k_y^2) <[Phi, C], [Phi, C]>_eta "
        "for any mode class, a k_z or a linear term appears, or the lattice second variation of the quartic term on a "
        "plane wave departs from the prediction by over 1 percent"
    )
    out["note"] = (
        "V4's second variation is a separate additive piece (the R23-0 / R24-0 Hessians); the reply's L_2 = "
        "(1/2)(omega^2 - k_perp^2) abs([Phi, C])^2 is this pencil up to normalization and the eta metric"
    )
    return out


# ================= e: the spin gate =================
def check_e():
    import sympy as sp

    J, C, Es = sp.symbols("J C E_stat", positive=True)
    E = Es + J**2 / (4 * C)
    om = sp.diff(E, J)
    sol = sp.solve(sp.Eq(2 * J * om / E, 1), J)
    Js = [s for s in sol if s.is_positive is not False]
    Jstar = sp.simplify(Js[0])
    out = {
        "J_star": str(Jstar),
        "J_star_sq_over_C_E_stat": str(sp.simplify(Jstar**2 / (C * Es))),
        "omega_star": str(sp.simplify(om.subs(J, Jstar))),
        "E_rot_over_E": str(sp.simplify((J**2 / (4 * C) / E).subs(J, Jstar))),
        "plan_formula_limit": "sqrt(C E_stat) is the E_rot << E limit of the same condition",
    }
    out["PASS"] = bool(
        sp.simplify(Jstar**2 - sp.Rational(4, 3) * C * Es) == 0
        and sp.simplify(om.subs(J, Jstar) - sp.sqrt(Es / (3 * C))) == 0
        and sp.simplify((J**2 / (4 * C) / E).subs(J, Jstar) - sp.Rational(1, 4)) == 0
    )
    out["fails_if"] = (
        "J_*^2 is not 4 C E_stat / 3, omega_* is not sqrt(E_stat / (3 C)), or E_rot / E is not 1/4"
    )
    return out


# ================= f: the frame-index bookkeeping of the pair textures =================
def _director_pair(X, Y, Z, d, like):
    def w_of(zc):
        zk = Z - zc
        rk = np.sqrt(X * X + Y * Y + zk * zk)
        return (X + 1j * Y) / (rk + zk + 1e-300)

    w1, w2 = w_of(+d / 2), w_of(-d / 2)
    w = w1 * w2 if like else w1 * np.conj(w2)
    den = 1.0 + np.abs(w) ** 2
    return np.stack([2 * w.real, 2 * w.imag, 1.0 - np.abs(w) ** 2], -1) / den[..., None]


def _frame_of(n):
    a = n[..., 0] / (1.0 + n[..., 2])
    b = n[..., 1] / (1.0 + n[..., 2])
    D = 1.0 + a * a + b * b
    e = np.stack([2 * (1 - a * a + b * b), -4 * a * b, -4 * a], -1) / (D * D)[..., None]
    nrm = np.linalg.norm(e, axis=-1)
    return e / np.maximum(nrm, 1e-300)[..., None], nrm


def _winding_line(e_loop):
    th = np.arctan2(e_loop[:, 1], e_loop[:, 0])
    d = np.diff(np.concatenate([th, th[:1]]))
    d = (d + np.pi / 2) % np.pi - np.pi / 2
    return float(np.sum(d) / np.pi)  # half-turns


def check_f():
    out = {}
    d = 12.0
    rng = np.random.default_rng(1)
    P = rng.uniform(-20, 20, size=(20000, 3))
    ok = True
    for like in (False, True):
        tag = "like" if like else "opposite"
        n = _director_pair(P[:, 0], P[:, 1], P[:, 2], d, like)
        e, nrm = _frame_of(n)
        orth = float(np.max(np.abs(np.sum(e * n, -1))))
        small = nrm < 1e-2
        rho = np.hypot(P[:, 0], P[:, 1])
        rec = {
            "max_abs_e_dot_n": orth,
            "singular_points_of_20000": int(small.sum()),
            "singular_points_max_rho": float(rho[small].max()) if small.any() else 0.0,
            "singular_points_z_range": (
                [float(P[small, 2].min()), float(P[small, 2].max())] if small.any() else None
            ),
            "loops": {},
        }
        for zc, name in (
            (0.0, "segment_mid"),
            (3.0, "segment_upper"),
            (-3.0, "segment_lower"),
            (9.0, "above_charge_1"),
            (-9.0, "below_charge_2"),
            (20.0, "far_above"),
            (-20.0, "far_below"),
        ):
            ph = np.linspace(0, 2 * np.pi, 721, endpoint=False)
            Xl, Yl, Zl = np.cos(ph), np.sin(ph), np.full_like(ph, zc)
            nl = _director_pair(Xl, Yl, Zl, d, like)
            el, nr = _frame_of(nl)
            axial = bool(np.all(np.abs(nl[:, 2]) > 0.95))
            rec["loops"][name] = {
                "z": zc,
                "n_z_mean": float(nl[:, 2].mean()),
                "axial_read": axial,
                "half_turns": _winding_line(el) if axial else None,
                "min_frame_norm": float(nr.min()),
            }
        ok = ok and orth < 1e-9
        out[tag] = rec
    out["reading"] = (
        "the frame's singular set follows n = -z: for both product textures it is the lower half-axis and a region "
        "below the pair, not the segment; a pair seed with the strand on the segment is a designed object behind the "
        "puncture gate (each core sphere index 2, the far sphere 0 for a +- pair, 4 for a like pair)"
    )
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "the pulled-back frame is not orthogonal to n (over 1e-9); the windings are reported, not gated"
    )
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    res = {}
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            res = json.load(f)
    if mode == "run":
        for nm, fn in (
            ("a", check_a),
            ("a2", check_a2),
            ("b", check_b),
            ("d", check_d),
            ("e", check_e),
            ("f", check_f),
        ):
            log(f"check {nm}")
            res[nm] = fn()
            print(nm, json.dumps(res[nm], indent=1)[:2500], flush=True)
    elif mode == "stack":
        log("check c")
        res["c"] = check_c()
        print("c", json.dumps(res["c"], indent=1), flush=True)
    else:
        raise SystemExit(f"unknown mode {mode}")
    res["PASS"] = all(res[k].get("PASS", True) for k in res if isinstance(res[k], dict))
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1)
    print("R25-0 PASS:", res["PASS"])


if __name__ == "__main__":
    main()
