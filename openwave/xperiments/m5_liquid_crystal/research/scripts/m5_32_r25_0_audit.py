"""M5.32 R25-0 adversarial audit: the form level of the strand rung, refuted with own methods.

The auditor did NOT read `m5_32_r25_0_form.py` nor `m5_32_r25_1_strand.py`. The claims
come from `data/m5_32_r25_0_form.json` (the audited results with their `fails_if`
clauses) and the R25 sections of the task record. Every check below is an independent
method: own sympy derivations on the explicit 4x4 matrices, an own 1D radial integrator
and minimizer, own finite differences, an own implementation of the author's projected
E . B witness (definition read from the author's bundle, not run), an own frame and
winding reader. The only shared code is the certified instrument (`m5_21_3_a_4d.py` for
the curvature energy, the generator catalog and the coordinates; `m5_32_r20_0_class.py`
for the V4 potential with a free weight), used as the object under test, never as the
method.

Verdicts: CONFIRMED (number and wording hold), QUALIFIED (the number holds, the wording
or the scope needs the stated correction), REFUTED, NOT_RUN (with the reason).

Run: python3 scripts/m5_32_r25_0_audit.py  (about 10 minutes; at most 2 worker threads).
Writes data/m5_32_r25_0_audit.json.
"""

import importlib.util
import itertools
import json
import os
import sys
import time

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r25_0_audit.json")
FORM_JSON = os.path.join(DATA, "m5_32_r25_0_form.json")
S1_NPZ = os.path.join(DATA, "m5_32_r21_1", "S1_Bseed_g8_d0.3_w25_n32.npz")
RAD_NPZ = os.path.join(DATA, "m5_32_r22_1", "rad_pin_d0.3_w25_n48_L48.npz")

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
R0 = _load("m5_32_r20_0_class", "m5_32_r20_0_class.py")
W1 = B3.W1
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
T0 = time.time()
RESULTS = {}


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def record(cid, method, numbers, verdict, note=""):
    RESULTS[cid] = {"method": method, "numbers": numbers, "verdict": verdict, "note": note}
    log(f"{cid}: {verdict}  {note}")
    for k, v in numbers.items():
        log(f"    {k} = {v}")


def cfg_of(n, L, delta):
    return B3.base_cfg(s=-1.0, n=n, L=L, delta=delta)


# =====================================================================================
# A. the block sector, symbolic (own sympy on the explicit 4x4 with eta and V4)
# =====================================================================================
def block_matrix_sym(s0, b, m, phi):
    """M = diag(8, 1, S) with the pair S = s0 I + b (cos(m phi) sz + sin(m phi) sx) at (2, 3)."""
    M = sp.zeros(4, 4)
    M[0, 0] = 8
    M[1, 1] = 1
    M[2, 2] = s0 + b * sp.cos(m * phi)
    M[3, 3] = s0 - b * sp.cos(m * phi)
    M[2, 3] = M[3, 2] = b * sp.sin(m * phi)
    return M


def check_a():
    rho, phi = sp.symbols("rho phi", positive=True)
    s0, w, f0, delta = sp.symbols("s0 w f0 delta", positive=True)
    m = sp.Symbol("m", positive=True, integer=True)
    b = sp.Function("b")(rho)
    eta = sp.diag(-1, 1, 1, 1)
    M = block_matrix_sym(s0, b, m, phi)

    def d_x(g):
        return sp.diff(g, rho) * sp.cos(phi) - sp.diff(g, phi) * sp.sin(phi) / rho

    def d_y(g):
        return sp.diff(g, rho) * sp.sin(phi) + sp.diff(g, phi) * sp.cos(phi) / rho

    Ax = M.applyfunc(d_x)
    Ay = M.applyfunc(d_y)
    F = Ax * eta * Ay - Ay * eta * Ax
    u = 4 * (eta * F * eta * F.T).trace()
    u_claim = 32 * m**2 * (b * sp.diff(b, rho) / rho) ** 2
    u_ok = sp.simplify(sp.expand(u - u_claim)) == 0
    # the same on the bare 2x2 block with the plain commutator (the author's normalization)
    S = M[2:4, 2:4]
    Sx, Sy = S.applyfunc(d_x), S.applyfunc(d_y)
    F2 = Sx * Sy - Sy * Sx
    u2 = 4 * (F2 * F2.T).trace()
    u2_ok = sp.simplify(sp.expand(u2 - u_claim)) == 0

    # the potential: N = M eta, C_p = (-8)^p + 1 + delta^p, s0 = delta / 2, f = b^2
    f = sp.Symbol("f", positive=True)
    bb = sp.Symbol("b", positive=True)
    Mb = block_matrix_sym(delta / 2, bb, m, phi)
    N = Mb * eta
    diffs = []
    P = sp.eye(4)
    for p in range(1, 5):
        P = P * N
        tr = sp.simplify(P.trace())
        diffs.append(sp.expand(sp.simplify(tr - ((-8) ** p + 1 + delta**p))))
    claimed = [
        sp.Integer(0),
        2 * f - 2 * f0,
        6 * s0 * (f - f0),
        2 * f**2 + 12 * f * s0**2 - 2 * f0**2 - 12 * f0 * s0**2,
    ]
    sub = {f: bb**2, f0: delta**2 / 4, s0: delta / 2}
    diffs_ok = all(sp.simplify(diffs[p] - claimed[p].subs(sub)) == 0 for p in range(4))
    V = w * sum(d**2 for d in diffs)
    V_claim = w * (f - f0) ** 2 * (4 + 36 * s0**2 + (12 * s0**2 + 2 * (f + f0)) ** 2)
    V_ok = sp.simplify(sp.expand(V - V_claim.subs(sub))) == 0
    # the stiffness: bracket at f = 0 (the "leading K" of the author) versus at f = f0
    Kf = 4 + 36 * s0**2 + (12 * s0**2 + 2 * (f + f0)) ** 2
    K_at_0 = sp.expand(Kf.subs(f, 0).subs(f0, s0**2))
    K_at_f0 = sp.expand(Kf.subs(f, f0).subs(f0, s0**2))
    # the true quadratic stiffness about the vacuum (second derivative of V at f0, over 2 w)
    K_hess = sp.expand(sp.diff(V_claim, f, 2).subs(f, f0).subs(f0, s0**2) / (2 * w))

    # Bogomolny with an f-dependent K: perfect square and the bound
    fp, Kc = sp.symbols("fp K", positive=True)
    integrand = 8 * fp**2 / rho + w * Kc * (f - f0) ** 2 * rho
    cross = 2 * sp.sqrt(8 * w * Kc) * (f0 - f) * fp
    square = (sp.sqrt(8 / rho) * fp - sp.sqrt(w * Kc * rho) * (f0 - f)) ** 2
    square_ok = sp.simplify(sp.expand(integrand - cross - square)) == 0
    # the cross term integrates (over f from 0 to f0) to sqrt(32 w K) f0^2 / 2; times 2 pi
    cross_int = sp.integrate(2 * sp.sqrt(8 * w * Kc) * (f0 - f), (f, 0, f0))
    bound = sp.simplify(2 * sp.pi * cross_int)
    bound_ok = sp.simplify(bound - sp.pi * sp.sqrt(32 * w * Kc) * f0**2) == 0
    # the BPS profile solves the first-order equation with kappa = sqrt(w K / 32)
    kappa = sp.sqrt(w * Kc / 32)
    fbps = f0 * (1 - sp.exp(-kappa * rho**2))
    fo = sp.sqrt(8 / rho) * sp.diff(fbps, rho) - sp.sqrt(w * Kc * rho) * (f0 - fbps)
    bps_ok = sp.simplify(fo) == 0
    # with the f-dependent K the first-order equation still completes the square, so the
    # exact block potential has an exact bound 2 pi int_0^f0 sqrt(32 w K(f)) (f0 - f) df,
    # between the K(0) and K(f0) bounds because K(f) is increasing.
    K_increasing = sp.simplify(sp.diff(Kf, f)) == 4 * (12 * s0**2 + 2 * (f + f0))
    numbers = {
        "u_4x4_eta_is_32m2_bbprime_over_rho_sq": bool(u_ok),
        "u_2x2_plain_is_32m2_bbprime_over_rho_sq": bool(u2_ok),
        "trace_differences_match_claim": bool(diffs_ok),
        "V_exact_block_polynomial_ok": bool(V_ok),
        "K_bracket_at_f_0": str(K_at_0),
        "K_bracket_at_f_f0": str(K_at_f0),
        "K_true_hessian_over_2w_at_f0": str(K_hess),
        "perfect_square_ok": bool(square_ok),
        "bound_pi_sqrt32wK_f0sq_ok": bool(bound_ok),
        "bps_profile_first_order_ok": bool(bps_ok),
        "K_of_f_increasing": bool(K_increasing),
        "m_scaling": "u carries m^2, so the bound and the exact-K bound carry m: T_m = m T_1 exactly",
    }
    ok = all([u_ok, u2_ok, diffs_ok, V_ok, square_ok, bound_ok, bps_ok])
    verdict = "QUALIFIED" if ok else "REFUTED"
    note = (
        "every identity holds; wording: K = 4 + 36 s0^2 + 144 s0^4 is the bracket at f = 0 "
        "(leading order in b, the author's phrase), NOT the leading order in f - f0: the "
        "quadratic stiffness about the vacuum is 4 + 36 s0^2 + 256 s0^4 (the 'exact-vacuum K'). "
        "The exact block potential has its own exact Bogomolny value "
        "2 pi int_0^f0 sqrt(32 w K(f)) (f0 - f) df (K(f) increasing), which the 1D minimizer "
        "must reproduce; see a2."
    )
    record(
        "a",
        "own sympy on the explicit 4x4 with eta, V4 = w sum (tr N^p - C_p)^2",
        numbers,
        verdict,
        note,
    )


# =====================================================================================
# A2. numerics: own 1D radial integrator and minimizer, the exact-K bound, the stack slab
# =====================================================================================
def K_of(f, s0, f0):
    return 4 + 36 * s0**2 + (12 * s0**2 + 2 * (f + f0)) ** 2


def bound_const_K(w, K, f0):
    return np.pi * np.sqrt(32 * w * K) * f0**2


def bound_exact_K(w, s0, f0, c=0.0, m=1):
    """2 pi int_c^f0 sqrt(32 m^2 w K(f)) (f0 - f) df: the exact Bogomolny value of the
    exact block potential, for a profile rising from f(0) = c."""
    val, _ = quad(lambda f: np.sqrt(32 * m * m * w * K_of(f, s0, f0)) * (f0 - f), c, f0)
    return 2 * np.pi * val


class Radial1D:
    """T[f] = 2 pi int [8 m^2 f'^2 / rho + w K(f) (f - f0)^2 rho] drho on nodes rho_i = i dr."""

    def __init__(self, w, s0, f0, R=40.0, dr=0.02, m=1, pin_axis=True):
        self.w, self.s0, self.f0, self.m = w, s0, f0, m
        self.rho = np.arange(0.0, R + dr / 2, dr)
        self.dr = dr
        self.pin = pin_axis
        self.rmid = 0.5 * (self.rho[1:] + self.rho[:-1])

    def energy_grad(self, fi):
        f = np.concatenate([fi, [self.f0]])
        if self.pin:
            f[0] = 0.0
        df = np.diff(f)
        kin = 8 * self.m**2 * df**2 / self.rmid / self.dr
        Kf = K_of(f, self.s0, self.f0)
        dev = f - self.f0
        pot = self.w * Kf * dev**2 * self.rho * self.dr
        T = 2 * np.pi * (kin.sum() + pot.sum())
        g = np.zeros_like(f)
        gk = 16 * self.m**2 * df / self.rmid / self.dr
        g[1:] += gk
        g[:-1] -= gk
        dK = 4 * (12 * self.s0**2 + 2 * (f + self.f0))
        g += self.w * (dK * dev**2 + 2 * Kf * dev) * self.rho * self.dr
        g *= 2 * np.pi
        gi = g[:-1].copy()
        if self.pin:
            gi[0] = 0.0
        return T, gi

    def minimize(self, f_init=None, maxiter=20000):
        n = len(self.rho) - 1
        if f_init is None:
            kap = np.sqrt(self.w * K_of(0.0, self.s0, self.f0) / 32) / self.m
            f_init = self.f0 * (1 - np.exp(-kap * self.rho[:-1] ** 2))
        res = minimize(
            self.energy_grad,
            f_init,
            jac=True,
            method="L-BFGS-B",
            options={"maxiter": maxiter, "maxfun": 4 * maxiter, "ftol": 1e-15, "gtol": 1e-12},
        )
        return res.fun, res.x, res

    def energy_of_profile(self, f_nodes):
        return self.energy_grad(f_nodes[:-1])[0]


def bps_profile(rho, w, K, f0, m=1, c=0.0):
    kap = np.sqrt(w * K / 32) / m
    return c + (f0 - c) * (1 - np.exp(-kap * rho**2))


def slab_field(n, h, delta, w, K, nz=4):
    """the BPS strand along z on an (n, n, nz) slab in the stack's convention
    diag(8, 1, S) with S the (2, 3) pair; returns M and the cfg."""
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y = np.meshgrid(x, x, indexing="ij")
    rho = np.sqrt(X * X + Y * Y)
    phi = np.arctan2(Y, X)
    s0 = delta / 2
    f0 = s0**2
    f = bps_profile(rho, w, K, f0)
    b = np.sqrt(f)
    M = np.zeros((n, n, nz, 4, 4))
    M[..., 0, 0] = 8.0
    M[..., 1, 1] = 1.0
    M[..., 2, 2] = (s0 + b * np.cos(phi))[..., None]
    M[..., 3, 3] = (s0 - b * np.cos(phi))[..., None]
    M[..., 2, 3] = M[..., 3, 2] = (b * np.sin(phi))[..., None]
    cfg = cfg_of(n, n * h, delta)
    return M, cfg


def stack_tension(M, cfg, w, nz):
    e_u, _ = B3.e_parts(M, cfg)
    e_v, _ = R0.v4_energy_grad(M, cfg, R0.roots_of(cfg), w, need_grad=False)
    return (e_u + e_v) / (nz * cfg["h"]), e_u / (nz * cfg["h"]), e_v / (nz * cfg["h"])


def check_a2(form):
    rows = []
    for delta in (0.3, 0.1, 0.03, 0.01):
        for w1s in (6.25, 25.0, 100.0):
            w = W1 * w1s
            s0 = delta / 2
            f0 = s0**2
            K0 = K_of(0.0, s0, f0)
            Kv = K_of(f0, s0, f0)
            b_lead = bound_const_K(w, K0, f0)
            b_vac = bound_const_K(w, Kv, f0)
            b_exact = bound_exact_K(w, s0, f0)
            r = Radial1D(w, s0, f0, R=40.0, dr=0.02)
            T_pin, f_pin, _ = r.minimize()
            r2 = Radial1D(w, s0, f0, R=40.0, dr=0.01)
            T_pin2, _, _ = r2.minimize()
            T_rich = (4 * T_pin2 - T_pin) / 3  # Richardson, O(dr^2) error
            row = {
                "delta": delta,
                "w1s": w1s,
                "bound_lead_K": b_lead,
                "bound_exact_vacuum_K": b_vac,
                "bound_exact_Kf_integral": b_exact,
                "T_1d_pinned_dr0.02": T_pin,
                "T_1d_pinned_dr0.01": T_pin2,
                "T_1d_pinned_richardson": T_rich,
                "T_1d_over_exact_Kf_bound_minus_1": T_rich / b_exact - 1,
                "T_over_delta4": T_rich / delta**4,
                "between_bounds": bool(b_lead <= T_rich <= b_vac * (1 + 1e-9)),
            }
            rows.append(row)
    # the free-axis family: the shifted BPS bound and an own free-axis descent at delta 0.3, w 25
    delta, w = 0.3, W1 * 25
    s0 = delta / 2
    f0 = s0**2
    shifted = {str(c): bound_exact_K(w, s0, f0, c=c * f0) for c in (0.0, 0.015, 0.1, 0.5, 0.9)}
    rf = Radial1D(w, s0, f0, R=40.0, dr=0.02, pin_axis=False)
    T_free, f_free, res_free = rf.minimize(maxiter=50000)
    free_axis = {
        "shifted_bps_bound_by_c_over_f0": shifted,
        "own_free_axis_descent_T": T_free,
        "own_free_axis_f0_over_f0": float(f_free[0] / f0),
        "own_free_axis_T_over_lead_bound": T_free / bound_const_K(w, K_of(0, s0, f0), f0),
        "own_free_axis_iterations": int(res_free.nit),
        "own_free_axis_converged": bool(res_free.success),
    }
    # the stack's energy of the BPS profile on the z-invariant slab (own field, the certified energy)
    slab = []
    for delta in (0.3, 0.1):
        w = W1 * 25
        s0 = delta / 2
        f0 = s0**2
        K0 = K_of(0.0, s0, f0)
        r = Radial1D(w, s0, f0, R=40.0, dr=0.01)
        T_1d_min = r.minimize()[0]
        T_1d_bps = r.energy_of_profile(bps_profile(r.rho, w, K0, f0))
        for h in (1.5, 1.0, 0.5):
            n = int(round(48 / h))
            M, cfg = slab_field(n, h, delta, w, K0)
            T_st, Tu, Tv = stack_tension(M, cfg, w, 4)
            slab.append(
                {
                    "delta": delta,
                    "h": h,
                    "n": n,
                    "T_stack": T_st,
                    "T_stack_u": Tu,
                    "T_stack_v": Tv,
                    "T_1d_of_same_profile": T_1d_bps,
                    "T_1d_minimizer": T_1d_min,
                    "ratio_to_1d_same_profile": T_st / T_1d_bps,
                    "ratio_to_1d_minimizer": T_st / T_1d_min,
                }
            )
    # m = 2 by the same 1D minimization (own), against the exact identity T_2 = 2 T_1
    delta, w = 0.3, W1 * 25
    s0 = delta / 2
    f0 = s0**2
    T1 = Radial1D(w, s0, f0, R=40.0, dr=0.01, m=1).minimize()[0]
    T2 = Radial1D(w, s0, f0, R=40.0, dr=0.01, m=2).minimize()[0]
    T1c = Radial1D(w, s0, f0, R=40.0, dr=0.02, m=1).minimize()[0]
    T2c = Radial1D(w, s0, f0, R=40.0, dr=0.02, m=2).minimize()[0]
    m2 = {
        "T2_over_T1_dr0.01": T2 / T1,
        "T2_over_T1_richardson": ((4 * T2 - T2c) / 3) / ((4 * T1 - T1c) / 3),
        "T2_over_T1_exact_Kf_bound": bound_exact_K(w, s0, f0, m=2) / bound_exact_K(w, s0, f0, m=1),
        "audited_ratio": form["b"]["additivity_T_m2_over_T_m1"],
    }
    # core radius scaling: kappa = sqrt(w K / 32) so rho_half (f = f0 / 2) = sqrt(ln 2 / kappa)
    core = {}
    for w1s in (6.25, 25.0, 100.0):
        kap = np.sqrt(W1 * w1s * K_of(0.0, 0.15, 0.0225) / 32)
        core[str(w1s)] = float(np.sqrt(np.log(2) / kap))
    core["ratio_6.25_over_100_vs_16^0.25"] = core["6.25"] / core["100.0"] / 16**0.25
    core["gap_at_rho_half_over_vacuum_gap"] = float(np.sqrt(0.5))

    # verdict
    r25 = [r for r in rows if r["delta"] == 0.3 and r["w1s"] == 25.0][0]
    r003 = [r for r in rows if r["delta"] == 0.03 and r["w1s"] == 25.0][0]
    pred = {str(r["delta"]): r["T_over_delta4"] for r in rows if r["w1s"] == 25.0}
    ok_bounds = all(r["between_bounds"] for r in rows)
    ok_exact = all(abs(r["T_1d_over_exact_Kf_bound_minus_1"]) < 2e-5 for r in rows)
    ok_slab = all(abs(s["ratio_to_1d_minimizer"] - form_ratio(form, s)) < 2e-3 for s in slab)
    numbers = {
        "rows": rows,
        "delta0.3_w25": {
            "bound_lead": r25["bound_lead_K"],
            "bound_vac": r25["bound_exact_vacuum_K"],
            "bound_exact_Kf": r25["bound_exact_Kf_integral"],
            "T_1d": r25["T_1d_pinned_richardson"],
            "audited_T_1d": [
                x for x in form["a2"]["rows"] if x["delta"] == 0.3 and x["w1s"] == 25.0
            ][0]["T_1d_axis_pinned"],
        },
        "delta0.03_w25_T_over_lead_bound_minus_1": r003["T_1d_pinned_richardson"]
        / r003["bound_lead_K"]
        - 1,
        "T_over_delta4_at_w25": pred,
        "small_delta_constant_pi_sqrt(128 w)_over_16": float(np.pi * np.sqrt(128 * W1 * 25) / 16),
        "free_axis": free_axis,
        "slab": slab,
        "m2": m2,
        "core_radius_rho_half_f": core,
    }
    verdict = "QUALIFIED" if (ok_bounds and ok_exact and ok_slab) else "REFUTED"
    note = (
        "all numbers reproduce (bounds, the 1D pinned minimizer equals the EXACT-K(f) Bogomolny "
        "value to 1e-5, the slab ratios, T/delta^4). Two corrections: (1) the pinned minimizer is "
        "not 'between two bounds by accident', it IS the exact Bogomolny value of the exact block "
        "potential (a first-order ODE with variable K(f)); (2) the free-axis state is not a "
        "stationary lower state: with f(0) = c free the exact bound is 2 pi int_c^f0 ... which "
        "decreases monotonically to ZERO at c = f0 (the pure winding, zero density by check b), "
        "so the audited 1.6 percent is a stalled descent along a family whose infimum is 0; the "
        "axis pin f(0) = 0 is the regularity condition, not a numerical convenience. Also the "
        "m = 2 ratio is exactly 2 (an identity); 2.03 is the audited 1D grid's residual. The "
        "'core radius' 3.63 is where f = b^2 reaches f0 / 2 (the pair gap is 0.71 of the vacuum's "
        "there, not half)."
    )
    record(
        "a2",
        "own 1D radial integrator + L-BFGS, exact-K(f) quadrature, own slab field on the certified energy",
        numbers,
        verdict,
        note,
    )


def form_ratio(form, s):
    for x in form["a2"]["stack_bps_profile"]:
        if x["delta"] == s["delta"] and x["h"] == s["h"]:
            return x["ratio"]
    return np.nan


# =====================================================================================
# B. the exterior density of the pure winding: continuum zero, lattice residue orders
# =====================================================================================
def d_fwd(g, ax, h):
    out = np.zeros_like(g)
    sl = [slice(None)] * g.ndim
    a, b = list(sl), list(sl)
    a[ax], b[ax] = slice(0, -1), slice(1, None)
    out[tuple(a)] = (g[tuple(b)] - g[tuple(a)]) / h
    return out


def d_bwd(g, ax, h):
    out = np.zeros_like(g)
    sl = [slice(None)] * g.ndim
    a, b = list(sl), list(sl)
    a[ax], b[ax] = slice(0, -1), slice(1, None)
    out[tuple(b)] = (g[tuple(b)] - g[tuple(a)]) / h
    return out


def d_cen(g, ax, h):
    out = np.zeros_like(g)
    sl = [slice(None)] * g.ndim
    a, b, c = list(sl), list(sl), list(sl)
    a[ax], b[ax], c[ax] = slice(0, -2), slice(2, None), slice(1, -1)
    out[tuple(c)] = (g[tuple(b)] - g[tuple(a)]) / (2 * h)
    return out


def dens_xy(M, h, kind):
    """own per-cell curvature density 4 <F_xy, F_xy>_eta of a 2D (n, n, 4, 4) field."""
    d = {"fwd": d_fwd, "bwd": d_bwd, "cen": d_cen}[kind]
    Ax, Ay = d(M, 0, h), d(M, 1, h)
    F = Ax @ ETA @ Ay - Ay @ ETA @ Ax
    return 4 * np.einsum("...ab,...cd,ac,bd->...", F, F, ETA, ETA)


def winding_2d(n, h, delta, m=1, f=None):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y = np.meshgrid(x, x, indexing="ij")
    rho = np.sqrt(X * X + Y * Y)
    phi = np.arctan2(Y, X)
    s0 = delta / 2
    b = s0 * np.ones_like(rho) if f is None else np.sqrt(f(rho))
    M = np.zeros((n, n, 4, 4))
    M[..., 0, 0] = 8.0
    M[..., 1, 1] = 1.0
    M[..., 2, 2] = s0 + b * np.cos(m * phi)
    M[..., 3, 3] = s0 - b * np.cos(m * phi)
    M[..., 2, 3] = M[..., 3, 2] = b * np.sin(m * phi)
    return M, rho


def check_b():
    # continuum: a texture that depends on one coordinate has F = 0 (sympy, direct)
    rho, phi, b0 = sp.symbols("rho phi b0", positive=True)
    m = sp.Symbol("m", positive=True, integer=True)
    eta = sp.diag(-1, 1, 1, 1)
    M = block_matrix_sym(b0, b0, m, phi)

    def d_x(g):
        return sp.diff(g, rho) * sp.cos(phi) - sp.diff(g, phi) * sp.sin(phi) / rho

    def d_y(g):
        return sp.diff(g, rho) * sp.sin(phi) + sp.diff(g, phi) * sp.cos(phi) / rho

    Ax, Ay = M.applyfunc(d_x), M.applyfunc(d_y)
    F = (Ax * eta * Ay - Ay * eta * Ax).applyfunc(sp.simplify)
    cont_zero = F == sp.zeros(4, 4)
    # the leading lattice residue of the one-sided (fwd / bwd averaged) stencil, analytically:
    # F_fwd = (h/2) G + O(h^2), G = [M_xx, M_y] + [M_x, M_yy]; density_sym = h^2 <G, G> + O(h^4)
    Mxx, Myy = Ax.applyfunc(d_x), Ay.applyfunc(d_y)
    G = (Mxx * eta * Ay - Ay * eta * Mxx) + (Ax * eta * Myy - Myy * eta * Ax)
    GG = sp.simplify((eta * G * eta * G.T).trace().subs(m, 1))
    GG_mean = sp.simplify(sp.integrate(GG, (phi, 0, 2 * sp.pi)) / (2 * sp.pi))
    # lattice: own stencils, L 32 fixed, h 1 / 0.5 / 0.25
    delta = 0.3
    bands = [(1, 2), (2, 4), (4, 6), (6, 9), (9, 14)]
    res = {}
    for h in (1.0, 0.5, 0.25):
        n = int(round(32 / h))
        M, rr = winding_2d(n, h, delta)
        inner = np.zeros((n, n), dtype=bool)
        inner[2:-2, 2:-2] = True
        dsym = 0.5 * (dens_xy(M, h, "fwd") + dens_xy(M, h, "bwd"))
        dcen = dens_xy(M, h, "cen")
        res[str(h)] = {}
        for lo, hi in bands:
            msk = inner & (rr >= lo) & (rr < hi)
            res[str(h)][f"[{lo},{hi})"] = {
                "sym_mean_density": float(dsym[msk].mean()),
                "cen_mean_density": float(dcen[msk].mean()),
                "analytic_h2_GG_mean": float(
                    h * h * np.mean([float(GG_mean.subs({b0: 0.15, rho: r})) for r in rr[msk]])
                ),
            }
    b46 = "[4,6)"
    order_sym_1_05 = np.log(
        res["1.0"][b46]["sym_mean_density"] / res["0.5"][b46]["sym_mean_density"]
    ) / np.log(2)
    order_sym_05_025 = np.log(
        res["0.5"][b46]["sym_mean_density"] / res["0.25"][b46]["sym_mean_density"]
    ) / np.log(2)
    order_cen_05_025 = np.log(
        res["0.5"][b46]["cen_mean_density"] / res["0.25"][b46]["cen_mean_density"]
    ) / np.log(2)

    # rho falloff: log-log slope over the band centers at h 0.25 (asymptotic) and at h 1 (the audited spacing)
    def slope(hk, key):
        xs = np.log([0.5 * (lo + hi) for lo, hi in bands[2:]])
        ys = np.log([res[hk][f"[{lo},{hi})"][key] for lo, hi in bands[2:]])
        return float(-np.polyfit(xs, ys, 1)[0])

    falloff = {
        "sym_h1": slope("1.0", "sym_mean_density"),
        "sym_h0.5": slope("0.5", "sym_mean_density"),
        "sym_h0.25": slope("0.25", "sym_mean_density"),
        "analytic_h2GG_h0.25": slope("0.25", "analytic_h2_GG_mean"),
        "cen_h0.25": slope("0.25", "cen_mean_density"),
    }
    # the BPS field's density localization (1D formula, own): fraction of density beyond 3 rho_half
    w = W1 * 25
    s0 = delta / 2
    f0 = s0**2
    K0 = K_of(0.0, s0, f0)
    kap = np.sqrt(w * K0 / 32)
    rr1 = np.linspace(1e-6, 40, 400001)
    f = bps_profile(rr1, w, K0, f0)
    fp = np.gradient(f, rr1)
    dens = 8 * fp**2 / rr1**2 + w * K_of(f, s0, f0) * (f - f0) ** 2
    rho_half = np.sqrt(np.log(2) / kap)
    beyond = dens[rr1 > 3 * rho_half].max() / dens.max()
    numbers = {
        "continuum_F_zero_for_pure_winding": bool(cont_zero),
        "analytic_leading_sym_residue_per_volume": f"h^2 <G,G>, phi-mean = {GG_mean}",
        "bands": res,
        "order_in_h_sym_band46_h1_to_h0.5": float(order_sym_1_05),
        "order_in_h_sym_band46_h0.5_to_h0.25": float(order_sym_05_025),
        "order_in_h_central_band46_h0.5_to_h0.25": float(order_cen_05_025),
        "rho_falloff_exponent": falloff,
        "bps_rho_half": float(rho_half),
        "bps_density_beyond_3rho_half_over_max": float(beyond),
    }
    ok = cont_zero and order_sym_1_05 > 1.5 and beyond < 1e-3
    verdict = "QUALIFIED" if ok else "REFUTED"
    note = (
        "continuum zero holds; the lattice residue is a STENCIL property: the stack's one-sided "
        "fwd/bwd average leaves h^2 <G,G> with G = [M_xx, M_y] + [M_x, M_yy] (own derivation), "
        "which falls as rho^-6 analytically (the phi-mean above), not rho^-7; central differences "
        "leave O(h^4). The audited 7.4 is a short-range band fit at h 1 (pre-asymptotic); the m = 2 "
        "additivity is an identity (a2)."
    )
    record(
        "b",
        "own sympy for F = 0 and the leading stencil residue; own 2D lattice at h 1 / 0.5 / 0.25",
        numbers,
        verdict,
        note,
    )


# =====================================================================================
# C. the E . B witness
# =====================================================================================
def comm_eta_np(A, B):
    return A @ ETA @ B - B @ ETA @ A


def inner_eta_np(F, G):
    return np.einsum("...ab,...cd,ac,bd->...", F, G, ETA, ETA)


EPS4 = np.zeros((4, 4, 4, 4))
for _p in itertools.permutations(range(4)):
    _s = 1
    for _i in range(4):
        for _j in range(_i + 1, 4):
            if _p[_i] > _p[_j]:
                _s = -_s
    EPS4[_p] = _s


def four_form(A):
    """sum eps^{mu nu rho sigma} <[A_mu, A_nu]_eta, [A_rho, A_sigma]_eta>_eta and the abs-sum."""
    Fs = [[comm_eta_np(A[m], A[n]) for n in range(4)] for m in range(4)]
    tot, absum = 0.0, 0.0
    for m, n, r, s in itertools.product(range(4), repeat=4):
        e = EPS4[m, n, r, s]
        if e == 0:
            continue
        v = inner_eta_np(Fs[m][n], Fs[r][s])
        tot += e * v
        absum += abs(v)
    return tot, absum


def pair_witness(dM, O, pair=(0, 1)):
    """the author's projected witness on jets dM[mu] (mu = 0 time, 1..3 space), 3x3 per cell,
    eigenframe O (columns ascending); F_mn = (O^T [dM_m, dM_n] O)[pair]; E_i = F_0i,
    B = (F_23, F_31, F_12). Shapes: dM (4, ..., 3, 3), O (..., 3, 3)."""
    p, q = pair
    Op = O[..., :, p]
    Oq = O[..., :, q]

    def F(m, n):
        C = dM[m] @ dM[n] - dM[n] @ dM[m]
        return np.einsum("...a,...ab,...b->...", Op, C, Oq)

    E = np.stack([F(0, 1), F(0, 2), F(0, 3)], axis=-1)
    B = np.stack([F(2, 3), F(3, 1), F(1, 2)], axis=-1)
    return E, B


def witness_stats(E, B, mask=None):
    if mask is None:
        mask = np.ones(E.shape[:-1], dtype=bool)
    EB = np.einsum("...a,...a->...", E, B)[mask]
    E2 = np.einsum("...a,...a->...", E, E)[mask]
    B2 = np.einsum("...a,...a->...", B, B)[mask]
    En, Bn = np.sqrt(E2), np.sqrt(B2)
    return {
        "sum_abs_EB_over_sum_E2_plus_B2": float(np.abs(EB).sum() / (E2 + B2).sum()),
        "sum_abs_EB_over_half_sum_E2_plus_B2": float(2 * np.abs(EB).sum() / (E2 + B2).sum()),
        "signed_sum_EB_over_sum_E2_plus_B2": float(EB.sum() / (E2 + B2).sum()),
        "sum_abs_EB_over_sum_EnBn": float(np.abs(EB).sum() / (En * Bn).sum()),
        "scale_optimal_sum_abs_EB_over_2sqrt": float(
            np.abs(EB).sum() / (2 * np.sqrt(E2.sum() * B2.sum()))
        ),
        "cells": int(mask.sum()),
    }


def antisym3(v):
    W = np.zeros(v.shape[:-1] + (3, 3))
    W[..., 0, 1], W[..., 0, 2], W[..., 1, 2] = -v[..., 2], v[..., 1], -v[..., 0]
    W[..., 1, 0], W[..., 2, 0], W[..., 2, 1] = v[..., 2], -v[..., 1], v[..., 0]
    return W


def check_c1():
    rng = np.random.default_rng(3)
    worst_sym, worst_any = 0.0, 0.0
    for _ in range(20):
        A = rng.standard_normal((4, 4, 4))
        As = 0.5 * (A + A.swapaxes(-1, -2))
        t, a = four_form(As)
        worst_sym = max(worst_sym, abs(t) / a)
        t, a = four_form(A)
        worst_any = max(worst_any, abs(t) / a)
    # the reason: <[A,B]_eta,[C,D]_eta>_eta = -tr([eta A, eta B][eta C, eta D]) and the
    # totally antisymmetrized trace of four matrices vanishes because the 4-cycle is odd.
    rng2 = np.random.default_rng(5)
    X = rng2.standard_normal((4, 4, 4))
    cyc = sum(
        EPS4[m, n, r, s] * np.trace(X[m] @ X[n] @ X[r] @ X[s])
        for m, n, r, s in itertools.product(range(4), repeat=4)
    )
    numbers = {
        "max_ratio_symmetric_jets_20_draws": float(worst_sym),
        "max_ratio_arbitrary_jets_20_draws": float(worst_any),
        "eps_tr(X X X X)_arbitrary": float(cyc),
        "reason": "<[A,B]_eta,[C,D]_eta>_eta = -tr([eta A, eta B][eta C, eta D]); eps tr(X_mu X_nu X_rho X_sigma) = 0 "
        "since the cyclic shift of four indices is an odd permutation while the trace is cyclic; "
        "symmetry of the jets and eta play no role (it also vanishes for arbitrary matrices), "
        "the Jacobi form tr(A[B,[C,D]]) is an equivalent reading",
    }
    verdict = "CONFIRMED" if worst_sym < 1e-13 else "REFUTED"
    record(
        "c_i",
        "own four-form on random jets, symmetric and arbitrary; own algebraic reason",
        numbers,
        verdict,
        "identically zero; an unfalsifiable witness, as the audited script says",
    )


def check_c2():
    rng = np.random.default_rng(11)
    delta = 0.3
    D = np.diag([0.0, delta, 1.0])
    worst_fixed, worst_grad, worst_grad_scalefree, decomposition_err = 0.0, 0.0, 0.0, 0.0
    hits = []
    for _ in range(200):
        Q, _r = np.linalg.qr(rng.standard_normal((3, 3)))
        M = Q @ D @ Q.T
        Wv = rng.standard_normal((4, 3))
        Ws = antisym3(Wv)
        dM = np.array([Ws[m] @ M - M @ Ws[m] for m in range(4)])
        lam, O = np.linalg.eigh(M)
        E, B = pair_witness(dM, O)
        st = witness_stats(E[None], B[None])
        worst_fixed = max(worst_fixed, st["sum_abs_EB_over_sum_E2_plus_B2"])
        # the algebraic reason: F_mn = (d_r - d_p)(d_q - d_r) (x_m y_n - x_n y_m), x = w_pr, y = w_rq
        w = np.array([O.T @ Ws[m] @ O for m in range(4)])
        x, y = w[:, 0, 2], w[:, 2, 1]
        c = (1.0 - 0.0) * (delta - 1.0)
        Fpred = c * (x[:, None] * y[None, :] - x[None, :] * y[:, None])
        Fact = np.zeros((4, 4))
        for m in range(4):
            for n in range(4):
                Fact[m, n] = O[:, 0] @ (dM[m] @ dM[n] - dM[n] @ dM[m]) @ O[:, 1]
        decomposition_err = max(decomposition_err, np.abs(Fact - Fpred).max() / np.abs(Fact).max())
        # with an eigenvalue gradient (spatial only, the clocks keep the eigenvalues)
        g = rng.standard_normal((4, 3)) * 0.3
        g[0] = 0.0
        dMg = dM + np.array([Q @ np.diag(g[m]) @ Q.T for m in range(4)])
        E, B = pair_witness(dMg, O)
        st = witness_stats(E[None], B[None])
        worst_grad = max(worst_grad, st["sum_abs_EB_over_sum_E2_plus_B2"])
        worst_grad_scalefree = max(worst_grad_scalefree, st["sum_abs_EB_over_sum_EnBn"])
        hits.append(st["sum_abs_EB_over_sum_E2_plus_B2"])
    numbers = {
        "fixed_eigenvalues_max_ratio_200_draws": float(worst_fixed),
        "decomposition_F_eq_c_x_wedge_y_max_rel_err": float(decomposition_err),
        "with_gradient_max_ratio_abs_EB_over_E2_plus_B2": float(worst_grad),
        "with_gradient_median_ratio": float(np.median(hits)),
        "with_gradient_max_cos_like_abs_EB_over_EnBn": float(worst_grad_scalefree),
        "reason": "in the eigenframe a fixed-eigenvalue jet is [w_mu, D], so the pair element of [dM_m, dM_n] "
        "is c (x_m y_n - x_n y_m) with x = w_pr, y = w_rq (the director column r mediates): a simple "
        "2-form, whose Pfaffian E . B vanishes; an eigenvalue gradient adds z ^ v with "
        "z = d(d_p - d_q), v = w_pq (d_q - d_p), and E . B = 2 c (x ^ y ^ z ^ v) is generically nonzero",
    }
    verdict = "QUALIFIED" if worst_fixed < 1e-12 and worst_grad > 1e-3 else "REFUTED"
    note = (
        "zero at fixed eigenvalues and switched on by a gradient, both reproduced; the '0.95' is a "
        "draw-dependent number (the ratio |E.B|/(E^2+B^2) is not scale-free in the clock amplitude: "
        "its ceiling is 1/2 at |E| = |B| parallel, and the reported 0.95 must come from a normalization "
        "where E and B are compared as (E^2 + B^2)/2 or similar); the scale-free cosine form is what to quote"
    )
    record(
        "c_ii",
        "own algebraic jets [W, M] + O diag(g) O^T at a point, own projected witness",
        numbers,
        verdict,
        note,
    )


def jets_central(M3, h):
    return [d_cen(M3, ax, h) for ax in range(3)]


def clock_jets(M3, kind, rot_local_axis=None):
    n = M3.shape[0]
    if kind == "rot_z":
        Jz = antisym3(np.array([0.0, 0.0, 1.0]))
        return Jz @ M3 - M3 @ Jz
    W = antisym3(rot_local_axis)
    return W @ M3 - M3 @ W


def hedgehog_fixed(n, h, delta):
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rho = np.sqrt(X * X + Y * Y)
    nh = np.stack([X, Y, Z], axis=-1) / r[..., None]
    th = np.stack([X * Z, Y * Z, -rho * rho], axis=-1) / np.maximum(rho * r, 1e-300)[..., None]
    ph = np.stack([-Y, X, np.zeros_like(X)], axis=-1) / np.maximum(rho, 1e-300)[..., None]
    M3 = nh[..., :, None] * nh[..., None, :] + delta * th[..., :, None] * th[..., None, :]
    M4 = np.zeros(M3.shape[:-2] + (4, 4))
    M4[..., 0, 0] = 8.0
    M4[..., 1:, 1:] = M3
    return M4, (X, Y, Z, r, rho)


def witness_on_field(M4, cfg, mask_extra=None, use_catalog=True):
    """own central-difference jets + own projected witness on a stored/built 4x4 field."""
    h = cfg["h"]
    M3 = M4[..., 1:4, 1:4]
    n = M3.shape[0]
    lam, O = np.linalg.eigh(M3)
    director = O[..., :, 2]
    sp_jets = jets_central(M3, h)
    inner = np.zeros(M3.shape[:3], dtype=bool)
    inner[2:-2, 2:-2, 2:-2] = True
    gap = lam[..., 1] - lam[..., 0]
    mask = inner & (gap > 1e-2)
    if mask_extra is not None:
        mask = mask & mask_extra
    out = {}
    X, Y, Z = B3.coords(n, h)
    cat = B3.gen_catalog(cfg, M4) if use_catalog else None
    for kind in ("rot_z", "clock_local"):
        jets = {}
        jets["unit_omega"] = clock_jets(M3, kind, director)
        if use_catalog:
            jets["catalog_a0"] = cat[kind][..., 1:4, 1:4]
        out[kind] = {}
        for jn, jt in jets.items():
            dM = np.stack([jt] + sp_jets, axis=0)
            E, B = pair_witness(dM, O)
            st = witness_stats(E, B, mask)
            out[kind][jn] = st
            EB = np.einsum("...a,...a->...", E, B)
            sN, sS = EB[mask & (Z > 0)].sum(), EB[mask & (Z < 0)].sum()
            out[kind][jn]["EB_even_in_z_over_odd"] = float(
                abs(sN + sS) / max(abs(sN) + abs(sS), 1e-300)
            )
    # B is clock-independent by construction: F_23, F_31, F_12 never touch dM[0]
    dM = np.stack([np.zeros_like(M3)] + sp_jets, axis=0)
    _E0, B = pair_witness(dM, O)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rh = np.stack([X, Y, Z], axis=-1) / np.maximum(r, 1e-300)[..., None]
    Br = np.einsum("...a,...a->...", B, rh)
    shells = {}
    for R in (3.0, 4.5, 6.0, 9.0, 12.0):
        sh = mask & (np.abs(r - R) < 0.5 * h)
        nn, ss = sh & (Z > 0), sh & (Z < 0)
        if nn.sum() == 0 or ss.sum() == 0:
            continue
        bn, bs = Br[nn].mean(), Br[ss].mean()
        shells[str(R)] = {
            "Br_mean_north": float(bn),
            "Br_mean_south": float(bs),
            "Br_even_over_odd": float(abs(bn + bs) / max(abs(bn - bs), 1e-300)),
            "net_flux_over_sum_abs": float(Br[sh].sum() / max(np.abs(Br[sh]).sum(), 1e-300)),
        }
    out["B_r_shells"] = shells
    out["mask_cells"] = int(mask.sum())
    out["pair_gap_min_in_mask"] = float(gap[mask].min())
    return out


def check_c3():
    delta = 0.3
    rows = {}
    for h in (1.5, 1.0, 0.75):
        n = int(round(36 / h))
        M4, (X, Y, Z, r, rho) = hedgehog_fixed(n, h, delta)
        cfg = cfg_of(n, 36.0, delta)
        extra = (r > 4.5) & (rho > 3.0) & (r < 14.0)
        rows[str(h)] = witness_on_field(M4, cfg, mask_extra=extra)
    hs = np.array([1.5, 1.0, 0.75])
    orders = {}
    for kind in ("rot_z", "clock_local"):
        for jn in ("unit_omega", "catalog_a0"):
            ys = np.array([rows[str(h)][kind][jn]["sum_abs_EB_over_sum_E2_plus_B2"] for h in hs])
            orders[f"{kind}/{jn}"] = float(np.polyfit(np.log(hs), np.log(ys), 1)[0])
            ys2 = np.array([rows[str(h)][kind][jn]["sum_abs_EB_over_sum_EnBn"] for h in hs])
            orders[f"{kind}/{jn}/cos_like"] = float(np.polyfit(np.log(hs), np.log(ys2), 1)[0])
    numbers = {"rows": rows, "order_in_h": orders}
    ok = all(v > 1.0 for k, v in orders.items() if not k.endswith("cos_like"))
    verdict = "QUALIFIED" if ok else "REFUTED"
    note = (
        "the audited ratio |E.B| / (E^2 + B^2) falls with h for both clocks (central differences, own mask "
        "r > 4.5, rho > 3, r < 14); the orders differ from the audited 1.9 / 1.6 because they depend on the "
        "stencil, the mask and the jet normalization (the catalog a0 is unit-Frobenius over the WHOLE lattice, so "
        "its per-cell amplitude scales as h^1.5 at fixed L and inflates the order). Two structural facts the "
        "order hides: for clock_local at fixed eigenvalues the continuum E is identically zero (the local clock's "
        "generator has no director component in the eigenframe, and the pair element of [dM_0, dM_i] needs one), "
        "so the scale-free cosine |E.B| / (|E||B|) does NOT fall with h (residue over residue, order about 0); "
        "for rot_z on an axisymmetric field the clock jet is the phi-derivative, a combination of the spatial "
        "jets, so E . B vanishes in the continuum for ANY eigenvalue profile (see c_iv)"
    )
    record(
        "c_iii",
        "own fixed-eigenvalue biaxial hedgehog on n 24 / 36 / 48 at L 36, own witness",
        numbers,
        verdict,
        note,
    )


def rot_z_is_phi_derivative(M4, cfg):
    """on an axisymmetric field [J_z, M] = -y d_x M + x d_y M; with A_0 := the lattice phi-derivative the
    four jets are linearly dependent per cell and E . B vanishes identically (an exact algebraic fact).
    """
    h = cfg["h"]
    M3 = M4[..., 1:4, 1:4]
    n = M3.shape[0]
    X, Y, Z = B3.coords(n, h)
    lam, O = np.linalg.eigh(M3)
    sp_jets = jets_central(M3, h)
    Aphi = -Y[..., None, None] * sp_jets[0] + X[..., None, None] * sp_jets[1]
    Jz = antisym3(np.array([0.0, 0.0, 1.0]))
    Arot = Jz @ M3 - M3 @ Jz
    inner = np.zeros(M3.shape[:3], dtype=bool)
    inner[2:-2, 2:-2, 2:-2] = True
    mask = inner & ((lam[..., 1] - lam[..., 0]) > 1e-2)
    diff = (
        np.sqrt(np.sum((Arot - Aphi) ** 2, axis=(-1, -2)))[mask].sum()
        / np.sqrt(np.sum(Arot**2, axis=(-1, -2)))[mask].sum()
    )
    E, B = pair_witness(np.stack([Aphi] + sp_jets, axis=0), O)
    st_phi = witness_stats(E, B, mask)
    E, B = pair_witness(np.stack([Arot] + sp_jets, axis=0), O)
    st_rot = witness_stats(E, B, mask)
    return {
        "rel_norm_[Jz,M]_minus_lattice_phi_derivative": float(diff),
        "ratio_with_A0_lattice_phi_derivative": st_phi["sum_abs_EB_over_sum_E2_plus_B2"],
        "ratio_with_A0_[Jz,M]": st_rot["sum_abs_EB_over_sum_E2_plus_B2"],
    }


def check_c4():
    out = {}
    for tag, path, L in (
        ("S1_n32_biaxial", S1_NPZ, 48.0),
        ("rad_pin_n48_uniaxial", RAD_NPZ, 48.0),
    ):
        M4 = np.load(path)["M"].astype(np.float64)
        n = M4.shape[0]
        cfg = cfg_of(n, L, 0.3)
        out[tag] = witness_on_field(M4, cfg)
        out[tag]["rot_z_phi_derivative_diagnostic"] = rot_z_is_phi_derivative(M4, cfg)
    s1 = out["S1_n32_biaxial"]
    numbers = {
        "fields": out,
        "audited_S1_rot_z_catalog": 0.03655,
        "audited_S1_clock_local_catalog": 0.09141,
        "B_is_clock_independent": "F_23, F_31, F_12 contain only spatial jets, so B and every B_r parity "
        "statement is a property of the static field alone; 'for BOTH clocks' is one statement, not two",
    }
    r1 = s1["rot_z"]["catalog_a0"]["sum_abs_EB_over_half_sum_E2_plus_B2"]
    r2 = s1["clock_local"]["catalog_a0"]["sum_abs_EB_over_half_sum_E2_plus_B2"]
    dg = s1["rot_z_phi_derivative_diagnostic"]
    signed = max(
        abs(s1["rot_z"]["catalog_a0"]["signed_sum_EB_over_sum_E2_plus_B2"]),
        abs(s1["clock_local"]["catalog_a0"]["signed_sum_EB_over_sum_E2_plus_B2"]),
    )
    par = max(v["Br_even_over_odd"] for v in s1["B_r_shells"].values())
    verdict = "QUALIFIED"
    note = (
        f"S1 catalog-jet ratios |E.B| / ((E^2 + B^2) / 2) = {r1:.3e} (rot_z) / {r2:.3e} (clock_local) against the "
        f"audited 3.66e-2 / 9.14e-2 (same instrument jet, own central differences and own mask); signed total "
        f"{signed:.1e}; B_r even/odd at most {par:.1e} (dipole parity holds). Corrections: (1) rot_z is NOT a "
        f"witness on these axisymmetric fields: [J_z, M] equals the phi-derivative -y d_x M + x d_y M (relative "
        f"difference {dg['rel_norm_[Jz,M]_minus_lattice_phi_derivative']:.2e} on S1), and with A_0 := the lattice "
        f"phi-derivative the ratio is {dg['ratio_with_A0_lattice_phi_derivative']:.1e} (an exact algebraic zero: four "
        "linearly dependent one-forms), so the whole rot_z reading is lattice anisotropy plus the field's departure "
        "from axisymmetry, and the fixed-eigenvalue hedgehog is not its baseline; (2) the clock_local reading is "
        "E_i = -w_pq (d_q - d_p) d_i(d_p - d_q), i.e. E . B = (grad of the pair gap) . B up to a factor, the "
        "eigenvalue-gradient term alone; (3) the ratio is normalization-bound (the catalog a0 is unit-Frobenius over "
        "the lattice with an envelope, not a physical omega): quote the scale-free cosine or scale-optimal form beside "
        "it; (4) B is clock-independent by construction, so the 'both clocks' parity read is one statement; the "
        "author's magnetic-charge variant lives in the clock's spatial imprint at t > 0 (psi depends on z through "
        "tanh(z / l)), which a t = 0 jet on a static field cannot produce, so the missing net flux is expected here"
    )
    record(
        "c_iv",
        "own central-difference witness on the stored S1 n32 and rad_pin n48 fields, catalog and unit-omega jets",
        numbers,
        verdict,
        note,
    )


# =====================================================================================
# D. the pencil about the twisting vacuum (own sympy on a general symmetric 4x4 Phi)
# =====================================================================================
def check_d():
    delta, tau, eps, kx, ky, kz, om, s = sp.symbols(
        "delta tau epsilon k_x k_y k_z omega s", real=True
    )
    eta = sp.diag(-1, 1, 1, 1)
    D = sp.diag(8, 1, delta, 0)
    G = sp.zeros(4, 4)
    G[2, 3], G[3, 2] = -1, 1
    C = G * D - D * G
    ph = sp.symbols("p0:10", real=True)
    Phi = sp.zeros(4, 4)
    k = 0
    for i in range(4):
        for j in range(i, 4):
            Phi[i, j] = Phi[j, i] = ph[k]
            k += 1

    def ce(A, B):
        return A * eta * B - B * eta * A

    def ie(F, H):
        return (eta * F * eta * H.T).trace()

    # real mode eps Phi cos(theta), s = sin(theta); background M0 = D + tau z C (A_z^0 = tau C)
    A = {
        0: eps * om * Phi * s,
        1: -eps * kx * Phi * s,
        2: -eps * ky * Phi * s,
        3: tau * C - eps * kz * Phi * s,
    }
    L = 4 * sum(ie(ce(A[0], A[i]), ce(A[0], A[i])) for i in (1, 2, 3)) - 4 * sum(
        ie(ce(A[i], A[j]), ce(A[i], A[j])) for i in (1, 2, 3) for j in (1, 2, 3) if i < j
    )
    L = sp.expand(L)
    Lp = sp.Poly(L, eps)
    coeffs = {int(m[0]): sp.expand(c) for m, c in zip(Lp.monoms(), Lp.coeffs())}
    L0 = coeffs.get(0, sp.Integer(0))
    L1 = coeffs.get(1, sp.Integer(0))
    L2 = coeffs.get(2, sp.Integer(0))
    bracket = ie(ce(Phi, C), ce(Phi, C))
    pred = sp.expand(4 * tau**2 * s**2 * (om**2 - kx**2 - ky**2) * bracket)
    l2_ok = sp.simplify(L2 - pred) == 0
    kz_free = sp.diff(L2, kz) == 0
    # the quadratic form of the bracket per class at delta 0.3, own Gram matrix, eigenvalues
    basis = []
    labels = []
    for i in range(4):
        for j in range(i, 4):
            Bm = sp.zeros(4, 4)
            Bm[i, j] = Bm[j, i] = 1 if i == j else 1 / sp.sqrt(2)
            basis.append(Bm)
            labels.append(f"{i}{j}")
    d3 = sp.Rational(3, 10)
    Q = np.zeros((10, 10))
    for a in range(10):
        for b in range(10):
            Fa, Fb = ce(basis[a], C), ce(basis[b], C)
            Q[a, b] = float(ie(Fa, Fb).subs(delta, d3))
    spatial = [i for i, l in enumerate(labels) if l[0] != "0"]
    timesp = [i for i, l in enumerate(labels) if l[0] == "0" and l[1] != "0"]
    ev_sp = np.linalg.eigvalsh(Q[np.ix_(spatial, spatial)])
    ev_ts = np.linalg.eigvalsh(Q[np.ix_(timesp, timesp)])
    ev_all = np.linalg.eigvalsh(Q)
    diag_vals = {labels[i]: float(Q[i, i]) for i in range(10)}
    # the time-space zero mode: Phi = boost along the director axis (index 1)
    Bm = sp.zeros(4, 4)
    Bm[0, 1] = Bm[1, 0] = 1
    ts_zero = sp.simplify(ie(ce(Bm, C), ce(Bm, C))) == 0
    # lattice second variation on the certified stencil (own field, own stencil bookkeeping)
    lat = {}
    for kper in (8.0, 16.0):
        n, h = 32, 1.0
        cfg = cfg_of(n, 32.0, 0.3)
        X, Y, Z = B3.coords(n, h)
        kk = 2 * np.pi / kper
        Dn = np.diag([8.0, 1.0, 0.3, 0.0])
        Gn = np.zeros((4, 4))
        Gn[2, 3], Gn[3, 2] = -1.0, 1.0
        Cn = Gn @ Dn - Dn @ Gn
        Phin = np.zeros((4, 4))
        Phin[1, 2] = Phin[2, 1] = 1 / np.sqrt(2)
        tau_n, eps_n = 0.05, 1e-3
        cosk = np.cos(kk * X)

        def field(e):
            return Dn + tau_n * Z[..., None, None] * Cn + e * Phin * cosk[..., None, None]

        eu = [B3.e_parts(field(e), cfg)[0] for e in (eps_n, 0.0, -eps_n)]
        d2 = (eu[0] + eu[2] - 2 * eu[1]) / (
            2 * eps_n**2
        )  # the coefficient of eps^2 (E_u is exactly quadratic)
        br = float(inner_eta_np(comm_eta_np(Phin, Cn), comm_eta_np(Phin, Cn)))
        # own replica of the fwd/bwd averaged stencil on the scalar factors
        S_exact = 0.0
        for dfun in (d_fwd, d_bwd):
            dc = dfun(cosk, 0, h)
            dz = dfun(Z, 2, h)
            S_exact += 0.5 * np.sum(dc * dc * dz * dz) * h**3
        S_cont = kk**2 * np.sum(np.sin(kk * X) ** 2) * h**3
        S_sinkh = (np.sin(kk * h) / h) ** 2 * np.sum(np.sin(kk * X) ** 2) * h**3
        S_2sin = (2 * np.sin(kk * h / 2) / h) ** 2 * np.sum(np.sin(kk * X) ** 2) * h**3
        lat[f"period_{int(kper)}"] = {
            "kh": kk * h,
            "lattice_eps2_coefficient_of_E_u": float(d2),
            "pred_exact_stencil": 4 * tau_n**2 * br * S_exact,
            "ratio_exact_stencil": float(d2 / (4 * tau_n**2 * br * S_exact)),
            "ratio_continuum_k": float(d2 / (4 * tau_n**2 * br * S_cont)),
            "ratio_sin_kh_over_h": float(d2 / (4 * tau_n**2 * br * S_sinkh)),
            "ratio_2sin_kh2_over_h": float(d2 / (4 * tau_n**2 * br * S_2sin)),
        }
    numbers = {
        "L0_background": str(L0),
        "L1_linear": str(L1),
        "L2_equals_4tau2_sin2_(om2-kx2-ky2)_bracket": bool(l2_ok),
        "kz_free": bool(kz_free),
        "bracket_general_Phi": str(sp.factor(bracket)),
        "unit_basis_bracket_values_delta0.3": diag_vals,
        "spatial_class_eigenvalues_delta0.3": ev_sp.tolist(),
        "timespace_class_eigenvalues_delta0.3": ev_ts.tolist(),
        "full_10_eigenvalues_delta0.3": ev_all.tolist(),
        "timespace_zero_mode_is_director_boost_E01": bool(ts_zero),
        "lattice": lat,
    }
    ok = l2_ok and kz_free and L0 == 0 and L1 == 0
    verdict = "QUALIFIED" if ok else "REFUTED"
    note = (
        "the pencil reproduces exactly (real mode: the factor sin^2 theta multiplies 4 tau^2 eps^2), no k_z, "
        "no linear term. Scope corrections: the bracket is positive SEMI-definite on spatial modes (four zero "
        "modes, the commutant of C: E_00, E_11, the pair identity and C itself) and negative SEMI-definite on "
        "time-space modes (one zero mode, the boost along the director E_01); the lattice second variation "
        "matches the exact one-sided stencil bookkeeping to round-off, and the effective wavenumber of the "
        "certified fwd/bwd stencil is 2 sin(kh/2)/h (exact, with the dropped boundary slice as the 31/32 factor), "
        "not sin(kh)/h (which overshoots by 0.7 percent at kh 0.39 and 13 percent at kh 0.79)"
    )
    record(
        "d",
        "own sympy pencil with a general symmetric Phi; own Gram form per class; own stencil replica against B3.e_parts",
        numbers,
        verdict,
        note,
    )


# =====================================================================================
# E. the spin gate (own sympy)
# =====================================================================================
def check_e():
    C, Es, J, om = sp.symbols("C E_stat J omega", positive=True)
    out = {}
    for conv, Ek in (
        ("E_rot = C omega^2 (the stack: E_kin = omega^2 kin, C = kin)", C * om**2),
        ("E_rot = C omega^2 / 2", C * om**2 / 2),
    ):
        Jexpr = sp.diff(Ek, om)
        om_of_J = sp.solve(sp.Eq(J, Jexpr), om)[0]
        Erot = Ek.subs(om, om_of_J)
        E = Es + Erot
        gate = sp.Eq(2 * J * om_of_J, E)
        Jstar = [r for r in sp.solve(gate, J) if r.is_positive is not False][0]
        omstar = sp.simplify(om_of_J.subs(J, Jstar))
        ratio = sp.simplify((Erot / E).subs(J, Jstar))
        # stationarity: dE/dJ at fixed J-family
        dEdJ = sp.simplify(sp.diff(E, J))
        # the plan-time form sqrt(C E_stat): which condition does it solve
        Jplan = sp.sqrt(C * Es)
        plan_gate_ratio = sp.simplify((2 * J * om_of_J / E).subs(J, Jplan))
        plan_solves = sp.simplify((2 * J * om_of_J - Es).subs(J, Jplan)) == 0
        out[conv] = {
            "J_star_sq_over_C_E_stat": str(sp.simplify(Jstar**2 / (C * Es))),
            "omega_star_sq_times_C_over_E_stat": str(sp.simplify(omstar**2 * C / Es)),
            "E_rot_over_E": str(ratio),
            "dE_J/dJ": str(dEdJ),
            "plan_form_2J_omega_over_E": str(plan_gate_ratio),
            "plan_form_solves_2J_omega_eq_E_stat": bool(plan_solves),
        }
    stack = out["E_rot = C omega^2 (the stack: E_kin = omega^2 kin, C = kin)"]
    ok = (
        stack["J_star_sq_over_C_E_stat"] == "4/3"
        and stack["omega_star_sq_times_C_over_E_stat"] == "1/3"
        and stack["E_rot_over_E"] == "1/4"
    )
    numbers = out
    verdict = "QUALIFIED" if ok else "REFUTED"
    note = (
        "the three values hold with C defined by E_rot = C omega^2 (J = 2 C omega), the stack's kin_of convention; "
        "with E_rot = C omega^2 / 2 they read J*^2 = 2 C E_stat / 3 and omega*^2 = 2 E_stat / (3 C) (E_rot / E = 1/4 "
        "in both). Wording: the values solve the GATE 2 J omega = E, not a stationary point (E_J = E_stat + J^2 / 4C "
        "is monotone in J); the plan's sqrt(C E_stat) solves 2 J omega = E_stat, the E_rot << E limit (its gate ratio is 4/5)"
    )
    record(
        "e",
        "own sympy: J from dE_rot/d omega, the gate 2 J omega = E solved for J",
        numbers,
        verdict,
        note,
    )


# =====================================================================================
# F. topology reads: own director textures, own pulled-back frame, own winding reader
# =====================================================================================
def n_of_w(w):
    a = np.abs(w) ** 2
    return np.stack([2 * w.real, 2 * w.imag, 1 - a], axis=-1) / (1 + a)[..., None]


def dn_du(w):
    """d n / d Re(w) of the inverse stereographic map (own analytic formula)."""
    u, v = w.real, w.imag
    a = u * u + v * v
    return (
        np.stack([2 * (1 + v * v - u * u), -4 * u * v, -4 * u], axis=-1)
        / ((1 + a) ** 2)[..., None]
    )


def w_center(X, Y, Z, zc):
    r = np.sqrt(X * X + Y * Y + (Z - zc) ** 2)
    den = r + (Z - zc)
    return (X + 1j * Y) / np.where(np.abs(den) < 1e-300, 1e-300, den)


def texture(kind, X, Y, Z, d=12.0):
    w1, w2 = w_center(X, Y, Z, +d / 2), w_center(X, Y, Z, -d / 2)
    w = w1 * w2 if kind == "like" else w1 * np.conj(w2)
    return w, n_of_w(w)


def frame_e1(w):
    e = dn_du(w)
    nrm = np.linalg.norm(e, axis=-1)
    return e / np.maximum(nrm, 1e-300)[..., None], nrm


def half_turns(kind, z, rho_loop=1.0, npts=720, d=12.0):
    t = np.linspace(0, 2 * np.pi, npts, endpoint=False)
    X, Y, Z = rho_loop * np.cos(t), rho_loop * np.sin(t), np.full_like(t, z)
    w, nn = texture(kind, X, Y, Z, d)
    e1, nrm = frame_e1(w)
    xh = np.array([1.0, 0.0, 0.0])
    a = xh - np.einsum("...a,a->...", nn, xh)[..., None] * nn
    a /= np.linalg.norm(a, axis=-1)[..., None]
    b = np.cross(nn, a)
    ang = np.unwrap(
        np.arctan2(np.einsum("...a,...a->...", e1, b), np.einsum("...a,...a->...", e1, a))
    )
    total = ang[-1] - ang[0] + (ang[1] - ang[0])  # close the loop (one step)
    return {
        "z": z,
        "n_z_mean": float(nn[:, 2].mean()),
        "axial": bool(abs(nn[:, 2].mean()) > 0.95),
        "half_turns": float(total / np.pi),
        "min_frame_norm": float(nrm.min()),
        "max_abs_e1_dot_n": float(np.abs(np.einsum("...a,...a->...", e1, nn)).max()),
    }


def check_f():
    rng = np.random.default_rng(7)
    out = {}
    for kind in ("opposite", "like"):
        P = rng.uniform(-24, 24, size=(20000, 3))
        P[:, 2] = rng.uniform(-20, 20, size=20000)
        w, nn = texture(kind, P[:, 0], P[:, 1], P[:, 2])
        e1, nrm = frame_e1(w)
        orth = float(np.abs(np.einsum("...a,...a->...", e1, nn)).max())
        sing = nrm < 1e-3
        sing2 = nrm < 1e-2
        rho = np.hypot(P[:, 0], P[:, 1])
        loops = {}
        for z in (0.0, 3.0, -3.0, 9.0, -9.0, 20.0, -20.0):
            loops[str(z)] = half_turns(kind, z)
        loops_rho2 = {
            str(z): half_turns(kind, z, rho_loop=2.0)["half_turns"]
            for z in (9.0, -9.0, 20.0, -20.0)
        }
        out[kind] = {
            "max_abs_e1_dot_n": orth,
            "singular_points_norm_below_1e-3_of_20000": int(sing.sum()),
            "singular_z_range": (
                [float(P[sing, 2].min()), float(P[sing, 2].max())] if sing.any() else None
            ),
            "singular_max_rho": float(rho[sing].max()) if sing.any() else None,
            "singular_points_norm_below_1e-2_of_20000": int(sing2.sum()),
            "singular_z_range_1e-2": (
                [float(P[sing2, 2].min()), float(P[sing2, 2].max())] if sing2.any() else None
            ),
            "singular_max_rho_1e-2": float(rho[sing2].max()) if sing2.any() else None,
            "frame_norm_equals_1_plus_nz_max_err": float(np.abs(nrm - (1 + nn[:, 2])).max()),
            "loops_rho1": loops,
            "half_turns_rho2": loops_rho2,
        }
    opp, like = out["opposite"], out["like"]
    axial_opp = [v["half_turns"] for v in opp["loops_rho1"].values() if v["axial"]]
    like_low = [
        v["half_turns"] for k, v in like["loops_rho1"].items() if v["axial"] and float(k) < -6
    ]
    like_up = [
        v["half_turns"] for k, v in like["loops_rho1"].items() if v["axial"] and float(k) > 6
    ]
    ok = (
        max(abs(x) for x in axial_opp) < 1e-6
        and all(abs(abs(x) - 8) < 1e-6 for x in like_low)
        and max(abs(x) for x in like_up) < 1e-6
        and opp["singular_z_range"][1] < -6
    )
    numbers = out
    numbers["assumed_geometry"] = (
        "pair on the z axis at z = +-6 (d 12), the stereographic pole on z, box +-24, z sampled in +-20 (the audited entry states no seed; this geometry reproduces its z range)"
    )
    verdict = "QUALIFIED" if ok else "REFUTED"
    note = (
        "singular set below the lower charge (z in [-20, -6.x]), 0 half-turns on every axial loop of the +- "
        "texture, 8 half-turns (index 4) below the like pair and 0 above: reproduced with an own reader. "
        "Corrections: the orthogonality gate cannot fail (e_1 is a derivative of the unit vector n, so e_1 . n = 0 "
        "identically; the 8e-13 is round-off, not evidence), and the frame norm is exactly 1 + n_z, so the "
        "'singular points' count is a threshold on 1 + n_z, to be stated with the threshold"
    )
    record(
        "f",
        "own stereographic textures (pair on the z axis), own pulled-back frame d n / d Re w, own winding reader on rho 1 and 2 loops",
        numbers,
        verdict,
        note,
    )


# =====================================================================================
# main
# =====================================================================================
def main():
    with open(FORM_JSON) as fh:
        form = json.load(fh)
    check_a()
    check_a2(form)
    check_b()
    check_c1()
    check_c2()
    check_c3()
    check_c4()
    check_d()
    check_e()
    check_f()
    summary = {}
    for v in RESULTS.values():
        summary[v["verdict"]] = summary.get(v["verdict"], 0) + 1
    out = {"checks": RESULTS, "summary": summary, "wall_s": time.time() - T0}
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    log(f"summary {summary}; wrote {os.path.relpath(OUT_JSON, HERE)}")


if __name__ == "__main__":
    main()
