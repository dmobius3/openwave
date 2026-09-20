"""M5.32 R24-0: the form level of the 2026-09-20 14:12 UTC reply on the coordination thread.

N = M eta, the vacuum spectrum (-g, 1, delta, delta), g = 8, delta = 0.3.
V4 = w sum_p (tr N^p - C_p)^2, L the unique linear invariant of R23-0, V = V4 - c L.

Checks (each one carries a `fails_if` line in the JSON):

a  the degenerate pair has no quadratic stiffness under V4: on the split
   delta +/- eps, V4 / w = 4 eps^4 (1 + 9 delta^2 + 36 delta^4) + O(eps^6); the
   reply's 3x3 closed form is this leading order, not an exact identity.
b  the master curve of the collapse read: with rho = sqrt(beta) r^2 / 2 the
   decaying solution of eps'' = eps / r^2 + beta r^2 eps obeys
   eps r^(2 nu - 1/2) = const * rho^nu K_nu(rho), nu = sqrt(5) / 4, a plateau
   2^(nu - 1) Gamma(nu) at small rho; sqrt(rho) K_nu(rho) is not the curve.
   The reply's cutoff table, threshold estimate, physical point and sector depth.
c  the Riesz cluster projector on N: P_i = v_i (eta v_i)^T / (v_i^T eta v_i),
   P_23 = I - P_g - P_1, B = P_23 N P_23 - (1/2) tr(P_23 N) P_23. B = 0 on the
   vacuum orbit, tr B^2 = 2 eps^2 on the split orbit, smooth through the pair's
   coalescence, discontinuous at the lambda_1 = lambda_3 crossing.
d  the witness (the mechanism that rejected substrate-framework's P240 L1 axis):
   a static boost ripple N(x) = Lam(theta(x)) N_split Lam^-1 keeps every tr N^p
   (so V4 and L cost nothing) and has one derivative direction (so
   F = [dM, dM]_eta = 0); on it tr(dB dB) = -2 abs(B e_a)^2 theta'^2 for a
   boost along the unit vector a of the doublet plane, so kappa_0 Tr(dB dB)
   with the covariant trace is unbounded below on any split background. The
   same ripple under Frobenius and under h = eta + 2 u u (u the timelike
   eigenvector). A finite-amplitude chain confirms the closed form.
   d2: the ripple laid on a stored R23 end field (c 1e-3, n 32), against the
   certified action's own response to the same ripple (the baseline: L_cert is
   itself indefinite off the block-diagonal sector).
e  the halo equation's elastic part re-derived: the second variation in eps of
   the static quartic on M_sp = delta I + (1 - delta) rr^T + eps(r) T, T the
   split tensor in the (theta, phi) frame, against the R22-1 audit's
   8 (1 - delta)^2 / r^2 [2 eps_r^2 + 2 eps^2 / r^2]; the linear term in eps.
   e2: the same second variation on the stack's own curvature density, three
   spacings. e3: the jet Hessian for a general split field and the lowest
   angular eigenvalue (spin-2 harmonics), which fixes A in
   eps'' = A eps / r^2 + beta r^2 eps.

   e4: the doublet against the massless modes at second order (the director
   modes and the time-space modes): a total divergence and an identical zero.

   e5: the doublet against the massive rr and tt modes, integrated on the
   l = 2 sections, and the adiabatic estimate of A_eff(r) for the E-type half.

Modes: run (default) | mix (e4 and e5, sympy, about 15 minutes) | stack (d2 and e2, the checks that import the production stack). Output: data/m5_32_r24_0_form.json. Runtime: minutes.
"""

import importlib.util
import json
import os
import sys

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.special import gamma, kv

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r24_0_form.json")
G_, D_ = 8.0, 0.3
W1 = 0.000724023879  # the certified V4 weight; the main rows run W1 x 25
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
NU = np.sqrt(5.0) / 4.0
K_EL = 8 * (1 - D_) ** 2
M2 = (G_ + D_) * (1 - D_)
RNG = np.random.default_rng(24)


def beta_of(c):
    return M2 * c / (2 * K_EL)


# ================= a =================
def check_a():
    e, d, g = sp.symbols("eps delta g", positive=True)
    out = {}
    for nm, pmax, lead in (("3x3_p3", 3, [1.0]), ("3x3_p4", 4, [1.0])):
        V = sp.expand(
            sum(((d + e) ** p + (d - e) ** p - 2 * d**p) ** 2 for p in range(1, pmax + 1))
        )
        claim = sp.expand(4 * e**4 * (1 + 9 * d**2 + 36 * d**4))
        out[nm] = {"V_over_w": str(sp.factor(V)), "minus_reply_form": str(sp.factor(V - claim))}
    V4 = sp.expand(
        sum(
            ((-g) ** p + 1 + (d + e) ** p + (d - e) ** p - ((-g) ** p + 1 + 2 * d**p)) ** 2
            for p in range(1, 5)
        )
    )
    P = sp.Poly(V4, e)
    out["4x4_p4"] = {
        "V4_over_w": str(sp.factor(V4)),
        "lowest_power_of_eps": min(m[0] for m in P.monoms()),
    }
    lead = sp.simplify(P.coeff_monomial(e**4) - 4 * (1 + 9 * d**2 + 36 * d**4))
    out["PASS"] = bool(out["4x4_p4"]["lowest_power_of_eps"] == 4 and lead == 0)
    out["fails_if"] = (
        "V4 on the split had an eps^2 term, or its eps^4 coefficient were not 4 (1 + 9 d^2 + 36 d^4)"
    )
    out["note"] = (
        "the reply's 'exactly' holds at leading order only: p 4 adds 4 eps^6 (12 delta^2 + eps^2)"
    )
    return out


# ================= b =================
def halo_solution(beta, r_in=1.0, r_out=None):
    r_out = r_out or 6.5 * beta**-0.25
    f = lambda r: np.sqrt(r) * kv(NU, np.sqrt(beta) * r * r / 2)  # noqa: E731
    h = 1e-5 * r_out
    y0 = [f(r_out), (f(r_out + h) - f(r_out - h)) / (2 * h)]
    return solve_ivp(
        lambda r, y: [y[1], y[0] / r**2 + beta * r**2 * y[0]],
        [r_out, r_in],
        y0,
        rtol=1e-11,
        atol=1e-300,
        dense_output=True,
    )


def master(rho):
    """rho^nu K_nu(rho), normalized to 1 at rho -> 0."""
    return rho**NU * kv(NU, rho) / (2 ** (NU - 1) * gamma(NU))


def check_b():
    out = {"rows": []}
    worst = {"rho^nu K_nu": 0.0, "sqrt(rho) K_nu": 0.0}
    for c in (3e-4, 1e-3, 3e-3):
        beta = beta_of(c)
        sol = halo_solution(beta)
        rr = np.linspace(1.5, 3.0 * beta**-0.25, 40)
        y = sol.sol(rr)[0] * rr ** (2 * NU - 0.5)
        rho = np.sqrt(beta) * rr**2 / 2
        for nm, cur in (
            ("rho^nu K_nu", master(rho)),
            ("sqrt(rho) K_nu", np.sqrt(rho) * kv(NU, rho)),
        ):
            q = y / cur
            worst[nm] = max(worst[nm], float(q.max() / q.min() - 1))
    out["spread_of_solution_over_curve"] = worst
    out["spread_definition"] = (
        "max / min - 1 of solution over curve, r from 1.5 to 3 R_c, worst of c 3e-4, 1e-3, 3e-3"
    )
    # the same spread inside the read window of R24-1 (the audit's definition: r 6 to 18 at c 1e-3, factor 1.138)
    win = {}
    for c in (3e-4, 1e-3, 3e-3):
        beta = beta_of(c)
        sol = halo_solution(beta, r_out=max(20.0, 3.0 * beta**-0.25))
        rr = np.linspace(6.0, 18.0, 49)
        y = sol.sol(rr)[0] * rr ** (2 * NU - 0.5)
        rho = np.sqrt(beta) * rr**2 / 2
        q = y / (np.sqrt(rho) * kv(NU, rho))
        win[str(c)] = float(q.max() / q.min())
    out["sqrt_rho_K_nu_drift_factor_r6_to_r18"] = win
    rk = {3e-4: 8.92, 1e-3: 7.16, 3e-3: 5.12}
    for c in (3e-4, 1e-3, 3e-3, 1e-2, 5.1e-6):
        Rc = beta_of(c) ** -0.25
        out["rows"].append(
            {"c": c, "R_c": float(Rc), "R_K_over_R_c": rk[c] / Rc if c in rk else None}
        )
    out["threshold_c_at_ratio"] = {
        str(q): float(2 * K_EL / M2 * (4.44 / q) ** -4) for q in (1.09, 1.10, 1.18)
    }
    out["sector_depth_501_c2"] = {
        "c_3e-2": 501 * 9e-4,
        "c_5.1e-6": 501 * 5.1e-6**2,
        "barrier": 27.9,
    }
    out["PASS"] = bool(worst["rho^nu K_nu"] < 1e-6 and worst["sqrt(rho) K_nu"] > 0.05)
    out["fails_if"] = "the ODE solution times r^0.618 were not rho^nu K_nu(rho) up to one constant"
    return out


# ================= c =================
def lorentz(scale=0.6):
    A = RNG.normal(size=(4, 4)) * scale
    return expm(A - ETA @ A.T @ ETA)


def spectral(N):
    lam, v = np.linalg.eig(N)
    if np.abs(lam.imag).max() > 1e-9:
        raise ValueError("complex spectrum: the second signature sector")
    lam, v = lam.real, v.real
    nrm = np.array([v[:, i] @ ETA @ v[:, i] for i in range(4)])
    it = int(np.argmin(nrm))
    rest = [i for i in range(4) if i != it]
    i1 = rest[int(np.argmax(lam[rest]))]
    proj = lambda i: np.outer(v[:, i], ETA @ v[:, i]) / (v[:, i] @ ETA @ v[:, i])  # noqa: E731
    return lam, v, it, i1, proj


def riesz_B(N):
    lam, v, it, i1, proj = spectral(N)
    P = np.eye(4) - proj(it) - proj(i1)
    return P @ N @ P - 0.5 * np.trace(P @ N) * P


def riesz_B_field(M):
    """B on a field M[..., 4, 4] (symmetric, real spectrum of M eta assumed)."""
    N = M @ ETA
    lam, v = np.linalg.eig(N)
    lam, v = lam.real, v.real
    ev = np.einsum("...ai,ab,...bi->...i", v, ETA, v)
    it = np.argmin(ev, axis=-1)
    lam_sp = np.where(np.arange(4) == it[..., None], -np.inf, lam)
    i1 = np.argmax(lam_sp, axis=-1)
    P = np.broadcast_to(np.eye(4), N.shape).copy()
    for idx in (it, i1):
        vi = np.take_along_axis(v, idx[..., None, None], axis=-1)[..., 0]
        nrm = np.einsum("...a,ab,...b->...", vi, ETA, vi)
        P -= np.einsum("...a,...b->...ab", vi, vi @ ETA) / nrm[..., None, None]
    PN = P @ N @ P
    return PN - 0.5 * np.trace(P @ N, axis1=-2, axis2=-1)[..., None, None] * P


def naive_B(N):
    lam, v, it, i1, proj = spectral(N)
    u = v[:, it] / np.sqrt(abs(v[:, it] @ ETA @ v[:, it]))
    n = v[:, i1] / np.sqrt(abs(v[:, i1] @ ETA @ v[:, i1]))
    P = np.eye(4) - np.outer(u, u) - np.outer(n, n)
    return P @ N @ P - 0.5 * np.trace(P @ N) * P


def N_split(eps, Lm):
    return Lm @ np.diag([-G_, 1.0, D_ + eps, D_ - eps]) @ np.linalg.inv(Lm)


def check_c():
    out = {}
    wr = wn = ws = wt = 0.0
    for _ in range(200):
        Lm = lorentz()
        N = N_split(0.0, Lm)
        ws = max(ws, float(np.abs(N @ ETA - (N @ ETA).T).max()))
        wr = max(wr, float(np.abs(riesz_B(N)).max()))
        wn = max(wn, float(np.abs(naive_B(N)).max()))
        e = RNG.uniform(0.01, 0.2)
        B = riesz_B(N_split(e, Lm))
        wt = max(wt, float(abs(np.trace(B @ B) - 2 * e * e)))
    out["vacuum_orbit_200"] = {
        "M_asymmetry": ws,
        "max_abs_B_riesz": wr,
        "max_abs_B_naive_euclid_outer": wn,
    }
    out["split_orbit_trB2_minus_2eps2"] = wt
    Lm = lorentz()
    M0 = N_split(0.0, Lm) @ ETA
    X = RNG.normal(size=(4, 4))
    X = X + X.T
    ts = np.linspace(-1e-3, 1e-3, 21)
    Bs = np.array([riesz_B((M0 + t * X) @ ETA).ravel() for t in ts])
    fit = np.array([np.polyval(np.polyfit(ts, Bs[:, k], 3), ts) for k in range(16)]).T
    out["coalescence_cubic_fit_residual"] = float(np.abs(Bs - fit).max())
    out["coalescence_B_scale"] = float(np.abs(Bs).max())

    def P2(N):
        lam, v, it, i1, proj = spectral(N)
        rest = [i for i in range(4) if i not in (it, i1)]
        return proj(rest[int(np.argmax(lam[rest]))])

    Ps = np.array([P2((M0 + t * X) @ ETA).ravel() for t in ts])
    out["single_projector_jump_across_coalescence"] = float(np.abs(Ps[11] - Ps[9]).max())
    jump = riesz_B(np.diag([-G_, D_ + 0.002, D_ + 0.01, D_ - 0.01])) - riesz_B(
        np.diag([-G_, D_ + 0.02, D_ + 0.01, D_ - 0.01])
    )
    out["jump_of_B_near_lambda1_crossing"] = float(np.abs(jump).max())
    out["PASS"] = bool(
        wr < 1e-8
        and wt < 1e-8
        and out["coalescence_cubic_fit_residual"] < 1e-3 * out["coalescence_B_scale"] + 1e-5
        and out["single_projector_jump_across_coalescence"] > 0.1
        and out["jump_of_B_near_lambda1_crossing"] > 1e-3
    )
    out["fails_if"] = (
        "B were nonzero on a boosted vacuum, tr B^2 differed from 2 eps^2, or B jumped at the coalescence"
    )
    return out


# ================= d =================
def boost_gen(a):
    K = np.zeros((4, 4))
    K[0, 1:] = a
    K[1:, 0] = a
    return K


def contractions(dB, N):
    lam, v, it, i1, proj = spectral(N)
    u = v[:, it] / np.sqrt(abs(v[:, it] @ ETA @ v[:, it]))
    h = ETA + 2 * np.outer(ETA @ u, ETA @ u)
    S = dB @ ETA
    return {
        "covariant_trace": float(np.trace(dB @ dB)),
        "frobenius": float(np.trace(dB @ dB.T)),
        "h_metric": float(np.einsum("ab,bc,cd,da->", S, h, S, h)),
        "h_min_eig": float(np.linalg.eigvalsh(h).min()),
    }


def check_d():
    out = {"rows": []}
    worst = 0.0
    for _ in range(40):
        eps = RNG.uniform(0.02, 0.2)
        a = RNG.normal(size=3)
        a /= np.linalg.norm(a)
        th0 = RNG.uniform(-1.0, 1.0)
        R = lorentz(0.0)  # identity; the rotation of the frame is carried by a
        Ns = np.diag([-G_, 1.0, D_ + eps, D_ - eps])
        Bs = np.diag([0.0, 0.0, eps, -eps])
        K = boost_gen(a)
        f = lambda th: expm(th * K) @ Ns @ expm(-th * K)  # noqa: E731
        d = 1e-5
        dB = (riesz_B(f(th0 + d)) - riesz_B(f(th0 - d))) / (2 * d)
        con = contractions(dB, f(th0))
        closed = -2 * float(np.sum((Bs[1:, 1:] @ a) ** 2))
        inv = max(
            abs(
                np.trace(np.linalg.matrix_power(f(th0 + 0.3), p))
                - np.trace(np.linalg.matrix_power(Ns, p))
            )
            for p in range(1, 5)
        )
        worst = max(worst, abs(con["covariant_trace"] - closed), abs(con["h_metric"] + closed))
        out["rows"].append(
            dict(con, eps=eps, a=a.tolist(), theta0=th0, closed_form=closed, trNp_drift=float(inv))
        )
        del R
    out["worst_abs_deviation_from_closed_form"] = float(worst)
    out["all_covariant_negative"] = bool(
        all(r["covariant_trace"] < 0 for r in out["rows"] if abs(r["closed_form"]) > 1e-6)
    )
    out["all_h_positive"] = bool(
        all(r["h_metric"] > 0 and r["h_min_eig"] > 0 for r in out["rows"])
    )
    # a finite-amplitude chain: theta(x) = A sin(k x) on a periodic 1D lattice, central differences
    eps, h = 0.13, 1.0
    Ns = np.diag([-G_, 1.0, D_ + eps, D_ - eps])
    K = boost_gen(np.array([0.0, 1.0, 0.0]))
    chain = []
    for A, kk in ((0.05, 0.2), (0.3, 0.2), (0.3, 0.8), (1.0, 1.5)):
        n = int(round(2 * np.pi / kk / h)) * 4
        x = np.arange(n) * h
        k_lat = 2 * np.pi * 4 / (n * h)
        th = A * np.sin(k_lat * x)
        B = np.array([riesz_B(expm(t * K) @ Ns @ expm(-t * K)) for t in th])
        dB = (np.roll(B, -1, 0) - np.roll(B, 1, 0)) / (2 * h)
        T = float(np.einsum("xab,xba->", dB, dB) * h)
        dth = np.roll(th, -1) - np.roll(th, 1)
        # exact on the lattice: tr (B(t2) - B(t1))^2 = -2 eps^2 sinh^2(t2 - t1), at least as negative as theta'^2
        chain.append(
            {
                "A": A,
                "k": k_lat,
                "sum_tr_dBdB": T,
                "closed_form_lattice": float(-2 * eps**2 * np.sum(np.sinh(dth) ** 2) / (4 * h)),
                "continuum_form": float(-2 * eps**2 * np.sum((dth / (2 * h)) ** 2) * h),
            }
        )
    out["finite_amplitude_chain"] = chain
    out["chain_worst_relative"] = float(
        max(abs(q["sum_tr_dBdB"] / q["closed_form_lattice"] - 1) for q in chain)
    )
    out["PASS"] = bool(
        worst < 1e-6
        and out["all_covariant_negative"]
        and out["all_h_positive"]
        and out["chain_worst_relative"] < 1e-6
    )
    out["fails_if"] = (
        "tr(dB dB) on the ripple were non-negative, or differed from -2 abs(B a)^2 theta'^2"
    )
    out["reading"] = (
        "with the covariant trace the static energy kappa_0 sum_i tr(d_i B d_i B) is unbounded below on any split "
        "background (theta' is free, V4, L and F are blind to the ripple); the covariant positive contraction "
        "h = eta + 2 u u gives +2 abs(B a)^2 theta'^2 instead; Frobenius is positive but equals that value only at "
        "theta = 0 (the audit: 2 eps^2 theta'^2 cosh(4 theta) along a doublet axis, not Lorentz invariant). The closed "
        "form is the rest-frame statement: on a boosted split background tr(dB dB) = theta'^2 tr([K, B]^2) is "
        "indefinite (the audit), and the unboundedness stands"
    )
    return out


def check_d2():
    """the ripple on a stored R23 end field: the certified action's own response beside the new term's."""

    def _load(name, fname):
        spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        return mod

    CS = _load("m5_32_r23_1_cscan", "m5_32_r23_1_cscan.py")
    c = 1e-3
    tag = CS.job_tag(
        {"seed": "rad", "delta": 0.3, "w1s": 25.0, "c": c, "n": 32, "L": 48.0, "src": None}
    )
    f = os.path.join(CS.OUT_NPZ, tag + ".npz")
    if not os.path.exists(f):
        return {"status": "SKIPPED", "why": "the local R23 end field is absent", "tag": tag}
    M = np.load(f)["M"]
    cfg = CS.R21.cfg_of(32, 48.0, G_, 0.3)
    p = CS.R21.params_of(G_, 0.3)
    pot = ("v4", CS.R0.roots_of(cfg, degenerate=True), CS.W1 * 25.0)
    h = cfg["h"]
    X, Y, Z = CS.B3.coords(32, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    mask = r < 0.5 * 48.0 - 3 * h  # the ripple stays off the pinned shell
    E0 = CS.energy_grad(M, cfg, p, pot, c, need_grad=False)[0]
    parts0 = CS.R20.energy_parts(M, cfg, p, pot)

    def T_of(Mf):
        B = riesz_B_field(Mf)
        T = 0.0
        for ax in range(3):
            dB = (np.roll(B, -1, ax) - np.roll(B, 1, ax)) / (2 * h)
            T += np.einsum("...ab,...ba->...", dB, dB)[2:-2, 2:-2, 2:-2].sum()
        return float(T * h**3)

    T0 = T_of(M)
    out = {"tag": tag, "E0": float(E0), "T0_static_field": T0, "rows": []}
    K = boost_gen(np.array([0.0, 1.0, 0.0]))
    for A, kk in ((0.02, np.pi / (4 * h)), (0.02, np.pi / (2 * h)), (0.05, np.pi / (2 * h))):
        th = A * np.sin(kk * X) * np.exp(-((r / 12.0) ** 4)) * mask
        Lam = np.array([expm(t * K) for t in th.ravel()]).reshape(th.shape + (4, 4))
        Mr = Lam @ M @ np.swapaxes(Lam, -1, -2)
        E1 = CS.energy_grad(Mr, cfg, p, pot, c, need_grad=False)[0]
        parts1 = CS.R20.energy_parts(Mr, cfg, p, pot)
        dT = T_of(Mr) - T0
        row = {
            "A": A,
            "k": float(kk),
            "dE_certified_plus_cL": float(E1 - E0),
            "dE_curv": float(parts1["E_curv"] - parts0["E_curv"]),
            "dV4": float(parts1["V"] - parts0["V"]),
            "dT_per_kappa0": dT,
        }
        row["kappa0_star"] = (
            None
            if dT >= 0 or row["dE_certified_plus_cL"] <= 0
            else float(-row["dE_certified_plus_cL"] / dT)
        )
        out["rows"].append(row)
    out["reading"] = (
        "dT < 0 is the witness on a relaxed field; dE_curv is the certified action's own answer to the same ripple "
        "(the baseline: it is indefinite off the block-diagonal sector); kappa0_star, where defined, is the kappa_0 "
        "above which the ripple lowers the total"
    )
    return out


# ================= e =================
def check_e():
    r, th, ph, d, e0 = sp.symbols("r theta phi delta eps0", positive=True)
    eps = sp.Function("eps")(r)
    er = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])
    et = sp.Matrix([sp.cos(th) * sp.cos(ph), sp.cos(th) * sp.sin(ph), -sp.sin(th)])
    ep = sp.Matrix([-sp.sin(ph), sp.cos(ph), 0])
    T = et * et.T - ep * ep.T
    M = d * sp.eye(3) + (1 - d) * er * er.T + eps * T
    J = sp.Matrix(
        [[er[i], et[i] / r, ep[i] / (r * sp.sin(th))] for i in range(3)]
    )  # d q_a / d x_i
    dq = [M.diff(q) for q in (r, th, ph)]
    dM = [sum((J[i, a] * dq[a] for a in range(3)), sp.zeros(3, 3)) for i in range(3)]
    dens = 0
    for i in range(3):
        for j in range(i + 1, 3):
            F = dM[i] * dM[j] - dM[j] * dM[i]
            dens += (F * F.T).trace()
    dens = dens.subs(ph, 0)
    e_, e1 = sp.symbols("e e1")
    dens = dens.subs(sp.Derivative(eps, r), e1).subs(eps, e_)
    base = sp.simplify(dens.subs({e_: 0, e1: 0}))
    kappa = sp.simplify(
        8 * (1 - d) ** 4 / r**4 / base
    )  # the stack's normalization (R22-1 identity)
    dens = sp.expand(dens * kappa)
    lin_e = sp.simplify(dens.diff(e_).subs({e_: 0, e1: 0}))
    lin_e1 = sp.simplify(dens.diff(e1).subs({e_: 0, e1: 0}))
    q_e1e1 = sp.simplify(dens.diff(e1, 2).subs({e_: 0, e1: 0}) / 2)
    q_ee = sp.simplify(dens.diff(e_, 2).subs({e_: 0, e1: 0}) / 2)
    q_ee1 = sp.simplify(dens.diff(e_).diff(e1).subs({e_: 0, e1: 0}))
    avg = lambda f: sp.simplify(sp.integrate(f * sp.sin(th), (th, 0, sp.pi)) / 2)  # noqa: E731
    out = {
        "normalization_kappa": str(kappa),
        "linear_in_eps": str(lin_e),
        "linear_in_eps_r": str(lin_e1),
        "coef_eps_r^2": str(q_e1e1),
        "coef_eps^2": str(q_ee),
        "coef_eps_eps_r": str(q_ee1),
        "audit_form": "8 (1 - delta)^2 / r^2 [2 eps_r^2 + 2 eps^2 / r^2]",
    }
    try:
        out["sphere_average_linear_in_eps"] = str(avg(lin_e))
        out["sphere_average_coef_eps_r^2"] = str(avg(q_e1e1))
    except Exception as ex:  # noqa: BLE001
        out["sphere_average_error"] = repr(ex)
    # the eps^2 coefficient away from the frame's poles: its value at the equator and its pole behavior
    out["coef_eps^2_at_equator"] = str(sp.simplify(q_ee.subs(th, sp.pi / 2)))
    out["coef_eps^2_times_sin^2_at_pole"] = str(
        sp.limit(sp.simplify(q_ee * sp.sin(th) ** 2), th, 0)
    )
    closed = 16 * (1 - d) ** 2 * (2 * sp.cot(th) ** 2 - 1) / r**4
    out["coef_eps^2_closed_form"] = "16 (1 - delta)^2 (2 cot^2 theta - 1) / r^4"
    out["coef_eps^2_closed_form_residual"] = str(sp.simplify(sp.trigsimp(q_ee - closed)))
    out["coef_eps^2_closed_form_numeric_residual"] = float(
        max(abs(sp.N((q_ee - closed).subs({th: t, r: 1, d: 0.3}))) for t in (0.2, 0.7, 1.1, 1.5))
    )
    target = 16 * (1 - d) ** 2 / r**2
    out["coef_eps_r^2_equals_audit"] = bool(sp.simplify(q_e1e1 - target) == 0)
    out["coef_eps^2_equator_equals_audit"] = bool(
        sp.simplify(q_ee.subs(th, sp.pi / 2) - target / r**2) == 0
    )
    out["PASS"] = (
        True  # a derivation: the comparison lines above are the result, whichever way they fall
    )
    out["fails_if"] = "n/a: reported, not gated; the collapse read (R24-1) carries the consequence"
    return out


def check_e2():
    """the second variation of check e against the stack's own curvature density on a lattice ansatz."""

    def _load(name, fname):
        spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        return mod

    CS = sys.modules.get("m5_32_r23_1_cscan") or _load("m5_32_r23_1_cscan", "m5_32_r23_1_cscan.py")
    out = {"eps0": 0.02, "lattices": []}
    for n in (32, 48, 64):
        out["lattices"].append(_e2_lattice(CS, n, 48.0, out["eps0"]))
    dev = [
        [abs(b["stack_second_variation"] / b["this_derivation"] - 1) for b in q["bands"]]
        for q in out["lattices"]
    ]
    out["relative_deviation_by_lattice_and_band"] = dev
    out["PASS"] = bool(
        dev[-1][0] < 0.05
        and all(dev[-1][k] < dev[0][k] for k in range(3))
        and all(
            np.sign(b["stack_second_variation"]) == np.sign(b["this_derivation"])
            for q in out["lattices"]
            for b in q["bands"]
        )
    )
    out["fails_if"] = (
        "on the finest lattice the equatorial band differed from 16 (1 - delta)^2 (2 cot^2 theta - 1) eps^2 / r^4 by "
        "over 5 percent, a band's deviation did not shrink from h 1.5 to h 0.75, or a sign disagreed"
    )
    return out


def _e2_lattice(CS, n, L, e0):
    cfg = CS.R21.cfg_of(n, L, G_, D_)
    pot = (
        "v4",
        CS.R0.roots_of(cfg, degenerate=True),
        0.0,
    )  # the potential off: the curvature density alone
    h = cfg["h"]
    X, Y, Z = CS.B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z) + 1e-30
    ct = Z / r
    st = np.sqrt(np.maximum(1 - ct * ct, 1e-30))
    cp, sn = X / (r * st), Y / (r * st)
    er = np.stack([st * cp, st * sn, ct], -1)
    et = np.stack([ct * cp, ct * sn, -st], -1)
    ep = np.stack([-sn, cp, np.zeros_like(cp)], -1)

    def field(q):
        M = np.zeros(r.shape + (4, 4))
        M[..., 0, 0] = G_
        sp_ = D_ * np.eye(3) + (1 - D_) * np.einsum("...a,...b->...ab", er, er)
        sp_ = sp_ + q * (
            np.einsum("...a,...b->...ab", et, et) - np.einsum("...a,...b->...ab", ep, ep)
        )
        M[..., 1:, 1:] = sp_
        return M

    d0, dp, dm = (CS.R20.density(field(q), cfg, pot) for q in (0.0, e0, -e0))
    second = (dp + dm - 2 * d0) / 2
    first = (dp - dm) / 2
    out = {"n": n, "L": L, "h": h, "bands": []}
    for lo, hi in ((0.0, 0.3), (0.3, 0.6), (0.6, 0.8)):
        k = (r > 9.0) & (r < 18.0) & (np.abs(ct) >= lo) & (np.abs(ct) < hi)
        pred = h**3 * e0**2 * 16 * (1 - D_) ** 2 * (2 * ct[k] ** 2 / st[k] ** 2 - 1) / r[k] ** 4
        audit = h**3 * e0**2 * 16 * (1 - D_) ** 2 / r[k] ** 4
        out["bands"].append(
            {
                "abs_cos_theta": [lo, hi],
                "cells": int(k.sum()),
                "stack_second_variation": float(second[k].sum()),
                "this_derivation": float(pred.sum()),
                "r22_audit_form": float(audit.sum()),
                "stack_first_variation_over_second": float(first[k].sum() / abs(second[k].sum())),
                "base_density_over_identity": float(
                    d0[k].sum() / (h**3 * 8 * (1 - D_) ** 4 / r[k] ** 4).sum()
                ),
            }
        )
    return out


def check_e3():
    """the full second variation for a general split field psi = a T1 + b T2 (T1, T2 the two split tensors of the
    (theta, phi) frame), as a quadratic form in the jets, and the lowest angular eigenvalue of its angular part.
    """
    from scipy.linalg import eigh

    r, th, ph, d = sp.symbols("r theta phi delta", positive=True)
    er = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])
    et = sp.Matrix([sp.cos(th) * sp.cos(ph), sp.cos(th) * sp.sin(ph), -sp.sin(th)])
    ep = sp.Matrix([-sp.sin(ph), sp.cos(ph), 0])
    T1, T2 = et * et.T - ep * ep.T, et * ep.T + ep * et.T
    M0 = d * sp.eye(3) + (1 - d) * er * er.T
    J = sp.Matrix([[er[i], et[i] / r, ep[i] / (r * sp.sin(th))] for i in range(3)])
    jets = sp.symbols("a a_r a_t a_p b b_r b_t b_p")
    a, a_r, a_t, a_p, b, b_r, b_t, b_p = jets
    dq = [
        M0.diff(q)
        + (a_r, a_t, a_p)[k] * T1
        + a * T1.diff(q)
        + (b_r, b_t, b_p)[k] * T2
        + b * T2.diff(q)
        for k, q in enumerate((r, th, ph))
    ]
    dM = [sum((J[i, k] * dq[k] for k in range(3)), sp.zeros(3, 3)) for i in range(3)]
    dens = 0
    for i in range(3):
        for j in range(i + 1, 3):
            F = dM[i] * dM[j] - dM[j] * dM[i]
            dens += 4 * (F * F.T).trace()  # the stack's normalization (check e: kappa = 4)
    dens = sp.expand(dens.subs(ph, 0))
    zero = {y: 0 for y in jets}
    lin = [sp.simplify(dens.diff(x).subs(zero)) for x in jets]
    H = {}
    for i in range(8):
        for j in range(i, 8):
            v = sp.simplify(
                sp.trigsimp(dens.diff(jets[i]).diff(jets[j]).subs(zero) / 2 * r**4 / (1 - d) ** 2)
            )
            if v != 0:
                H[f"{jets[i]} {jets[j]}"] = v
    claim = {
        "a a": 32 / sp.sin(th) ** 2 - 48,
        "b b": 32 / sp.sin(th) ** 2 - 48,
        "a_r a_r": 16 * r**2,
        "b_r b_r": 16 * r**2,
        "a_t a_t": 8,
        "b_t b_t": 8,
        "a_p a_p": 8 / sp.sin(th) ** 2,
        "b_p b_p": 8 / sp.sin(th) ** 2,
        "a a_t": 48 / sp.tan(th),
        "b b_t": 48 / sp.tan(th),
        "a b_p": 16 * sp.cos(th) / sp.sin(th) ** 2,
        "a_p b": -16 * sp.cos(th) / sp.sin(th) ** 2,
        "a_t b_p": 24 / sp.sin(th),
        "a_p b_t": -24 / sp.sin(th),
    }
    resid = 0.0
    for k in set(H) | set(claim):
        f = H.get(k, 0) - claim.get(k, 0)
        resid = max(
            resid,
            (
                max(abs(float(sp.N(f.subs({th: t, r: 1.7})))) for t in (0.3, 0.9, 1.4, 2.5))
                if f != 0
                else 0.0
            ),
        )
    out = {
        "linear_terms_all_zero": bool(all(x == 0 for x in lin)),
        "hessian_entries_times_r4_over_(1-delta)^2": {k: str(v) for k, v in claim.items()},
        "hessian_residual_against_the_listed_form": resid,
    }

    def lowest(m, N=1500):
        t = np.linspace(0, np.pi, N + 2)[1:-1]
        h = t[1] - t[0]
        s_, c_ = np.sin(t), np.cos(t)
        W = 8 * m * m / s_**2 + 32 / s_**2 - 48 + 32 * m * c_ / s_**2
        C = 96 * c_ / s_ + 48 * m / s_
        A = np.zeros((N, N))
        sm = np.sin(0.5 * (t[1:] + t[:-1]))
        for i in range(N - 1):
            k = 8 * sm[i] / h
            A[i, i] += k
            A[i + 1, i + 1] += k
            A[i, i + 1] -= k
            A[i + 1, i] -= k
        A[0, 0] += 8 * np.sin(0.5 * h) / h
        A[-1, -1] += 8 * np.sin(np.pi - 0.5 * h) / h
        A += np.diag(W * s_ * h)
        for i in range(N):
            w = s_[i] * C[i] / 4
            if i + 1 < N:
                A[i, i + 1] += w
                A[i + 1, i] += w
            if i >= 1:
                A[i, i - 1] -= w
                A[i - 1, i] -= w
        return eigh(A, np.diag(s_ * h), eigvals_only=True, subset_by_index=[0, 1]).tolist()

    out["lowest_angular_eigenvalues_by_m"] = {str(m): lowest(m) for m in range(-3, 4)}
    lam = min(v[0] for v in out["lowest_angular_eigenvalues_by_m"].values())
    out["lowest_angular_eigenvalue"] = lam
    out["closed_form"] = "8 (l (l + 1) - 4) + 32 on spin-2 harmonics, l >= 2: 48, 96, 160"
    out["A_of_the_halo_equation"] = lam / 16
    out["nu"] = float(np.sqrt(1 + 4 * lam / 16) / 4)
    out["tail_exponent"] = float(0.5 - 2 * out["nu"])
    out["PASS"] = bool(out["linear_terms_all_zero"] and resid < 1e-9 and abs(lam - 48) < 0.05)
    out["fails_if"] = (
        "a linear term survived, the jet Hessian differed from the listed form, or the lowest mode were not 48"
    )
    out["reading"] = (
        "the radial part is the R22-1 audit's (16 eps_r^2); the angular part is 48 eps^2 / r^2 on the lowest spin-2 "
        "mode, not 16: eps'' = 3 eps / r^2 + beta r^2 eps, the same beta and cutoff, nu = sqrt(13) / 4, tail r^(-1.303)"
    )
    return out


def check_e4():
    """the split doublet against the modes the potential leaves massless, at second order (added at the R24 go:
    check e3 restricted the perturbation to the doublet, and a mixing with a massless mode could lower A).

    Spatial block: the director modes u, v (rhat -> rhat + u that + v phat). Full 4x4: the time-space modes
    w_r, w_t, w_p and M_00 (s), with the production contraction 4 F_ab F_cd eta^ac eta^bd.
    """
    r, th, ph, d = sp.symbols("r theta phi delta", positive=True)
    er = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])
    et = sp.Matrix([sp.cos(th) * sp.cos(ph), sp.cos(th) * sp.sin(ph), -sp.sin(th)])
    ep = sp.Matrix([-sp.sin(ph), sp.cos(ph), 0])
    T1, T2 = et * et.T - ep * ep.T, et * ep.T + ep * et.T
    U1, U2 = (1 - d) * (er * et.T + et * er.T), (1 - d) * (er * ep.T + ep * er.T)
    M0 = d * sp.eye(3) + (1 - d) * er * er.T
    J = sp.Matrix([[er[i], et[i] / r, ep[i] / (r * sp.sin(th))] for i in range(3)])

    def density(basis, jets, M_bg, eta, point=None):
        n = M_bg.shape[0]
        dq = []
        for k, q in enumerate((r, th, ph)):
            m = M_bg.diff(q)
            for f, T in enumerate(basis):
                m += jets[4 * f + 1 + k] * T + jets[4 * f] * T.diff(q)
            dq.append(m)
        dM = [
            sum((J[i, k] * dq[k] for k in range(3)), sp.zeros(n, n)).subs(ph, 0) for i in range(3)
        ]
        if point is not None:
            dM = [m.subs(point) for m in dM]
        dens = 0
        for i in range(3):
            for j in range(i + 1, 3):
                F = dM[i] * eta * dM[j] - dM[j] * eta * dM[i]
                dens += 4 * (F * eta * F.T * eta).trace()
        return sp.expand(dens)

    # spatial block, symbolic: the mixed part against a total divergence on the sphere
    jets = sp.symbols("a a_r a_t a_p b b_r b_t b_p u u_r u_t u_p v v_r v_t v_p")
    a, a_r, a_t, a_p, b, b_r, b_t, b_p, u, u_r, u_t, u_p, v, v_r, v_t, v_p = jets
    dens = density([T1, T2, U1, U2], jets, M0, sp.eye(3))
    split, direc = jets[:8], jets[8:]
    mixed = sum(
        dens.diff(x).diff(y).subs({z: 0 for z in jets}) * x * y for x in split for y in direc
    )
    div = (
        sp.cot(th) * (a * u + b * v)
        + a_t * u
        + a * u_t
        + b_t * v
        + b * v_t
        + (-a_p * v - a * v_p + b_p * u + b * u_p) / sp.sin(th)
    )
    resid = sp.simplify(sp.trigsimp(sp.expand(mixed - 16 * (1 - d) ** 3 / r**4 * div)))
    mixed_nonzero = sp.simplify(mixed) != 0
    # the full 4x4 at rational points: the blocks between the doublet or the director and the time-space modes
    g = sp.Integer(8)
    E = lambda i, j: sp.Matrix(
        4, 4, lambda p_, q_: 1 if (p_, q_) in ((i, j), (j, i)) else 0
    )  # noqa: E731

    def emb(T):
        m = sp.zeros(4, 4)
        m[1:, 1:] = T
        return m

    def tvec(e):
        m = sp.zeros(4, 4)
        for i in range(3):
            m[0, i + 1] = e[i]
            m[i + 1, 0] = e[i]
        return m

    M4 = emb(M0)
    M4[0, 0] = g
    basis4 = [emb(T1), emb(T2), emb(U1), emb(U2), tvec(er), tvec(et), tvec(ep), E(0, 0)]
    names = ["a", "b", "u", "v", "w_r", "w_t", "w_p", "s"]
    jets4 = sp.symbols(" ".join(f"{n_} {n_}_r {n_}_t {n_}_p" for n_ in names))
    eta4 = sp.diag(-1, 1, 1, 1)
    worst_w, worst_s = sp.Integer(0), sp.Integer(0)
    for pt in (
        {r: sp.Rational(17, 10), th: sp.Rational(9, 10), d: sp.Rational(3, 10)},
        {r: sp.Rational(23, 10), th: sp.Rational(21, 10), d: sp.Rational(3, 10)},
    ):
        d4 = density(basis4, jets4, M4, eta4, point=pt)
        zero = {z: 0 for z in jets4}
        for x in jets4[:16]:
            for y in jets4[16:28]:
                worst_w = max(worst_w, abs(sp.N(d4.diff(x).diff(y).subs(zero))))
            for y in jets4[28:]:
                worst_s = max(worst_s, abs(sp.N(d4.diff(x).diff(y).subs(zero))))
    out = {
        "spatial_mixed_block_is_nonzero_pointwise": bool(mixed_nonzero),
        "spatial_mixed_block_minus_16(1-delta)^3/r^4_div(W)": str(resid),
        "W": "W_theta = a u + b v, W_phi = b u - a v: the contraction of the split tensor with the director vector",
        "time_space_blocks_max_abs": float(worst_w),
        "M00_blocks_max_abs": float(worst_s),
        "reading": (
            "the doublet-director block is a total divergence on the sphere, so it integrates to zero at every r; the "
            "blocks with the time-space modes vanish identically (M -> T M T, T = diag(-1, 1, 1, 1), keeps the density, "
            "the background and the doublet, and flips w); what is left mixes only with modes V4 makes massive"
        ),
    }
    out["PASS"] = bool(resid == 0 and worst_w < 1e-12)
    out["fails_if"] = (
        "the doublet-director block differed from the listed divergence, or a block with the time-space modes survived"
    )
    return out


def check_e5():
    """the split doublet against the massive diagonal modes (rr and tt) at second order, and the adiabatic estimate of
    what they do to A (raised by the R24 audit, re-derived here). The l = 2, m = 0 sections: E-type a = alpha sin^2,
    B-type b = alpha sin^2; the scalars P(r) P_2 on rr and Q(r) P_2 on tt; M_00 is slaved through the V4 Hessian.
    """
    r, th, ph, d = sp.symbols("r theta phi delta", positive=True)
    er = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])
    et = sp.Matrix([sp.cos(th) * sp.cos(ph), sp.cos(th) * sp.sin(ph), -sp.sin(th)])
    ep = sp.Matrix([-sp.sin(ph), sp.cos(ph), 0])
    basis = [et * et.T - ep * ep.T, et * ep.T + ep * et.T, er * er.T, et * et.T + ep * ep.T]
    M0 = d * sp.eye(3) + (1 - d) * er * er.T
    J = sp.Matrix([[er[i], et[i] / r, ep[i] / (r * sp.sin(th))] for i in range(3)])
    names = ["a", "b", "p", "q"]
    jets = sp.symbols(" ".join(f"{n_} {n_}_r {n_}_t {n_}_p" for n_ in names))
    dq = []
    for k, x in enumerate((r, th, ph)):
        m = M0.diff(x)
        for f, T in enumerate(basis):
            m += jets[4 * f + 1 + k] * T + jets[4 * f] * T.diff(x)
        dq.append(m)
    dM = [sum((J[i, k] * dq[k] for k in range(3)), sp.zeros(3, 3)).subs(ph, 0) for i in range(3)]
    dens = 0
    for i in range(3):
        for j in range(i + 1, 3):
            F = dM[i] * dM[j] - dM[j] * dM[i]
            dens += 4 * (F * F.T).trace()
    dens = sp.expand(dens)
    zero = {z: 0 for z in jets}
    mixed = sum(dens.diff(x).diff(y).subs(zero) * x * y for x in jets[:8] for y in jets[8:])
    al, P, Q = [sp.Function(n_)(r) for n_ in ("alpha", "P", "Q")]
    P2 = (3 * sp.cos(th) ** 2 - 1) / 2
    integ = {}
    for label, (fa, fb) in (
        ("E", (al * sp.sin(th) ** 2, 0)),
        ("B", (0, al * sp.sin(th) ** 2)),
    ):
        sub = {}
        for n_, f in zip(names, (fa, fb, P * P2, Q * P2)):
            js = sp.symbols(f"{n_} {n_}_r {n_}_t {n_}_p")
            f = sp.sympify(f)
            sub.update(
                {js[0]: f, js[1]: sp.diff(f, r), js[2]: sp.diff(f, th), js[3]: sp.Integer(0)}
            )
        integ[label] = sp.simplify(
            sp.integrate(
                sp.expand(mixed.subs(sub, simultaneous=True) * sp.sin(th)), (th, 0, sp.pi)
            )
            * 2
            * sp.pi
        )
    claim_E = -sp.Rational(512, 5) * sp.pi * (1 - d) ** 2 * (P - Q) * al / r**4
    resid_E = sp.simplify(integ["E"] - claim_E)
    # eliminating (s00, P, Q) against the Hessian of V4 - c L on the diagonal modes:
    # dA = -192 (1 - delta)^2 v^T H^-1 v / r^4, v = (0, 1, -1)
    rows = []
    for c in (0.0, 1e-3, 3e-3, 8e-3):
        lam = np.array([-G_, 1.0, D_])
        dl = np.array([-1.0, 1.0, 2.0])
        H = np.zeros((3, 3))
        for pw in range(1, 5):
            gp = pw * lam ** (pw - 1) * dl
            H += 2 * W1 * 25.0 * np.outer(gp, gp)
        qp = lambda x: (x - 1) * (x - D_) + (x + G_) * (x - D_) + (x + G_) * (x - 1)  # noqa: E731
        H += -c * np.diag([qp(-G_), qp(1.0), 2 * qp(D_)])
        v = np.array([0.0, 1.0, -1.0])
        s_ = float(v @ np.linalg.solve(H, v))
        rows.append(
            {
                "c": c,
                "softest_massive_eigenvalue": float(np.linalg.eigvalsh(H)[0]),
                "v_Hinv_v": s_,
                "A_eff_E_type": {
                    str(rr): float(3 - 192 * (1 - D_) ** 2 * s_ / rr**4) for rr in (6, 9, 12, 18)
                },
            }
        )
    out = {
        "integrated_mixed_E_type": str(integ["E"]),
        "integrated_mixed_E_type_minus_claim": str(resid_E),
        "integrated_mixed_B_type": str(integ["B"]),
        "adiabatic_rows": rows,
        "reading": (
            "the E-type half of the l = 2 doublet couples to the massive (rr - tt) mode, the B-type half does not; "
            "to leading adiabatic order A_eff = 3 - 192 (1 - delta)^2 v H^-1 v / r^4 for the E-type half: about 2.99 at "
            "r 18, 2.93 at r 12, 2.77 at r 9, and the estimate stops being perturbative near r 6; the B-type half "
            "stays at 3. A = 3 is the far-field value, the inner third of the read window sits below it"
        ),
    }
    out["PASS"] = bool(resid_E == 0 and integ["B"] == 0)
    out["fails_if"] = (
        "the integrated E-type block differed from the listed form, or the B-type block did not vanish"
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
            ("b", check_b),
            ("c", check_c),
            ("d", check_d),
            ("e", check_e),
            ("e3", check_e3),
        ):
            res[nm] = fn()
            print(nm, json.dumps(res[nm], indent=1)[:3000], flush=True)
    elif mode == "mix":
        only = sys.argv[2] if len(sys.argv) > 2 else None
        for nm, fn in (("e4", check_e4), ("e5", check_e5)):
            if only and nm != only:
                continue
            res[nm] = fn()
            print(nm, json.dumps(res[nm], indent=1), flush=True)
    elif mode == "stack":
        for nm, fn in (("d2", check_d2), ("e2", check_e2)):
            res[nm] = fn()
            print(nm, json.dumps(res[nm], indent=1), flush=True)
    res["PASS"] = all(res[k].get("PASS", True) for k in res if isinstance(res[k], dict))
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1)
    print("R24-0 PASS:", res["PASS"])


if __name__ == "__main__":
    main()
