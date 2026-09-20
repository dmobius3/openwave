"""M5.32 R23-0: the form level of the second-order stiffness on the degenerate pair.

EQUATIONS FIRST
---------------
N = M eta, vacuum spectrum q = (-g, 1, delta, delta), V4 = w sum_p (tr N^p - C_p)^2.
(a) A potential V(t_1..t_4), t_p = tr N^p, has dV / d lambda_i = sum_p V_p p
    lambda_i^(p-1). The spectrum stays critical iff the cubic sum_p p V_p x^(p-1)
    vanishes at the three distinct eigenvalues, so the first derivatives have ONE
    free direction, a_p with sum_p p a_p x^(p-1) = (x + g)(x - 1)(x - delta).
(b) L = sum_p a_p (t_p - C_p) on the split (-g, 1, delta + eps, delta - eps) is
    exactly (delta + g)(delta - 1) eps^2 + eps^4 / 2, so V = V4 - c L (c >= 0)
    costs c (g + delta)(1 - delta) eps^2 - c eps^4 / 2 + V4's eps^4.
(c) The Hessian of V4 - c L on the ten symmetric entries at the vacuum: 5 zero
    modes (the Lorentz orbit: 6 generators minus the (2,3) rotation), the split
    doublet {M_23, (M_22 - M_33) / sqrt 2} at (g + delta)(1 - delta) c, three
    massive modes; c_crit where the lowest massive mode reaches zero.
(d) The vacuum as the global minimum of V4 - c L: on real spectra (many random
    starts) and on full symmetric 4x4 matrices (complex spectra allowed).
(e) The reply's M^2 / V_0 = 3.1e-3 turned into c for each candidate V_0 and for
    both conventions of the split amplitude.
(f) The halo equation from the R22-1 audit's second variation (check Q5),
        e_2 = (K / r^2) [2 eps_r^2 + 2 eps^2 / r^2],  K = 8 (1 - delta)^2,
        eps'' = eps / r^2 + beta r^2 eps,  beta = m2 / 2K,  m2 = (g + delta)(1 - delta) c:
    the scaling r -> r beta^(1/4) removes beta, so every radius of the decaying
    solution goes as c^(-1/4) exactly; the table gives what the R23-1 reader
    (eps r^0.618 falling to half its r 4.5 value) returns on that solution, the
    angular factor of the potential term left as f in [0.5, 2].
    Each check names the input that would make it fail.

AUDIT NOTES (m5_32_r23_0_audit.py, 25 own checks, 23 PASS; the two FAILs are wording)
(a) Criticality is the matrix condition Q(N) = 0 on the minimal polynomial. It equals
    the eigenvalue condition only where N is diagonalizable, which the vacuum orbit is:
    a symmetric M with the vacuum spectrum and a Jordan block has |grad L| =
    (g + delta)(1 - delta) 2e. The one-dimensional statement holds on the orbit
    (rank 3 of the 10 x 4 Jacobian on 200 boosted vacua).
(c) The orbit tangent space has rank 5; the (2,3) rotation generates a tangent of norm
    exactly 0, so a stabilizer removes a zero mode, it does not add one. c_crit by a
    root solve: 0.086281 (delta 0.3), 6.9021e-3 (row A), 1.04333e-6 (row B); the
    destabilizing mode is mostly M_11.
(d) The reading of the V < 0 states holds (a complex pair whose plane has signature
    (1,1); none within Frobenius distance 9.7 of M0, the first at 10.8; none on 200000
    block-diagonal matrices), with three corrections. The barrier is 27.90 = 1541 w per
    unit volume, not the straight-line 2e5: a block-diagonal path crosses -g with a
    spatial eigenvalue at -6.6 and lands on a twin vacuum diag(-delta, -g, 1, delta)
    with V = 0 for every c, which is a second global minimum INSIDE the instrument's
    sector and a saddle of the full problem (unstable along M_0i at -5.81 c). The pair
    is "next to delta" only at small c (0.2992 +/- 0.072i at c 3e-5, 0.2855 +/- 0.670i
    at c 3e-3). Frobenius distance to the orbit is not Lorentz invariant, so the
    invariant statement is the collision barrier. A descent from the seed never climbs
    27.9 per unit volume, so the rows stay in the vacuum's sector.
(f) Closed form eps = sqrt(r) K_nu(sqrt(beta) r^2 / 2), nu = sqrt(5) / 4; the reader
    radii are 18.512, 13.982, 10.977, 8.630, 7.154 (this script's values sit up to
    0.0024 higher from its reader grid). K and the angular factor were not audited.
(e) The re-solved V_0 is 0.00995457 at M_00 = 8.000348 (this script's grid gives
    0.00995476).

Regenerate: python m5_32_r23_0_form.py   (about 2 min)
"""

import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, "..", "data", "m5_32_r23_0_form.json")
G_, D_ = 8.0, 0.3
W1 = 0.000724023879
ALPHA = 1.0 / 137.035999
D_PHYS = 1.0 - (ALPHA / (64.0 * np.pi)) ** 0.25
C_SCAN = (3e-5, 1e-4, 3e-4, 1e-3, 3e-3)
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def a_of(g, d):
    x = sp.symbols("x")
    q = sp.Poly(sp.expand((x + g) * (x - 1) * (x - d)), x)
    return [q.coeff_monomial(x ** (p - 1)) / p for p in range(1, 5)]


def check_a():
    g, d = sp.symbols("g delta", positive=True)
    rows = [[p * xi ** (p - 1) for p in range(1, 5)] for xi in (-g, 1, d)]
    A = sp.Matrix(rows)
    ns = A.nullspace()
    a = sp.Matrix(a_of(g, d))
    member = sp.simplify(A * a) == sp.zeros(3, 1)
    prop = sp.simplify(ns[0] * a[3] - a * ns[0][3]) == sp.zeros(4, 1)
    # negative control: a_2 moved by 1 breaks criticality at all three eigenvalues
    b = a + sp.Matrix([0, 1, 0, 0])
    broken = [sp.simplify(v) for v in (A * b)]
    out = {
        "nullspace_dim": len(ns),
        "a_p": [str(sp.factor(v)) for v in a],
        "member_is_critical": bool(member),
        "nullspace_proportional_to_a": bool(prop),
        "negative_control_residuals": [str(v) for v in broken],
        "PASS": bool(len(ns) == 1 and member and prop and all(v != 0 for v in broken)),
        "fails_if": "a second independent null vector existed, or (x+g)(x-1)(x-delta) were not critical",
    }
    print("a:", out, flush=True)
    return out


def check_b():
    g, d, e = sp.symbols("g delta epsilon", positive=True)
    a = a_of(g, d)
    spec0, spec = (-g, 1, d, d), (-g, 1, d + e, d - e)
    L = sum(a[p - 1] * (sum(s**p for s in spec) - sum(s**p for s in spec0)) for p in range(1, 5))
    want = (d + g) * (d - 1) * e**2 + e**4 / 2
    ok = sp.simplify(sp.expand(L) - sp.expand(want)) == 0
    # L is first order along a generic (non-split) spectral direction only through squares: dL/dlambda_i = 0
    lam = sp.symbols("l0:4")
    Lg = sum(a[p - 1] * sum(s**p for s in lam) for p in range(1, 5))
    grad0 = [sp.simplify(sp.diff(Lg, v).subs(dict(zip(lam, spec0)))) for v in lam]
    out = {
        "L_on_split": str(sp.factor(sp.expand(L))),
        "PASS": bool(ok and all(v == 0 for v in grad0)),
        "coefficient_g8_d0.3": float((D_ + G_) * (1 - D_)),
        "coefficient_g8_dphys": float((D_PHYS + G_) * (1 - D_PHYS)),
        "fails_if": "the eps^2 coefficient differed from (g + delta)(1 - delta) or a first-order term survived",
    }
    print("b:", out)
    return out


IDX = [(i, j) for i in range(4) for j in range(i, 4)]


def v_of_M(M, w, c, d, a):
    spec0 = np.array([-G_, 1.0, d, d])
    N = M @ ETA
    P, t = np.eye(4), []
    for _ in range(4):
        P = P @ N
        t.append(np.trace(P))
    r = [t[k] - np.sum(spec0 ** (k + 1)) for k in range(4)]
    return w * sum(v * v for v in r) - c * sum(a[k] * r[k] for k in range(4))


def sym_of(v):
    M = np.zeros((4, 4), dtype=np.result_type(v, float))
    for k, (i, j) in enumerate(IDX):
        M[i, j] = M[j, i] = v[k] / (1.0 if i == j else np.sqrt(2.0))
    return M


def hessian_parts(d):
    """exact (sympy) Hessian on the ten orthonormal symmetric coordinates, H = w A - c B."""
    v = sp.symbols("v0:10")
    M = sp.zeros(4, 4)
    for k, (i, j) in enumerate(IDX):
        M[i, j] = M[j, i] = v[k] if i == j else v[k] / sp.sqrt(2)
    dq = sp.Rational(d).limit_denominator(10**12)
    N = (sp.diag(8, 1, dq, dq) + M) * sp.diag(-1, 1, 1, 1)
    spec0 = (-8, 1, dq, dq)
    a = a_of(sp.Integer(8), dq)
    P, r = sp.eye(4), []
    for p in range(1, 5):
        P = P * N
        r.append(sp.expand(P.trace()) - sum(s**p for s in spec0))
    zero = {s: 0 for s in v}
    A, B = sp.zeros(10, 10), sp.zeros(10, 10)
    for k in range(4):
        dr = [sp.diff(r[k], s).subs(zero) for s in v]
        for i_ in range(10):
            for j_ in range(10):
                A[i_, j_] += 2 * dr[i_] * dr[j_]  # r_k = 0 on the vacuum
                B[i_, j_] += a[k] * sp.diff(r[k], v[i_], v[j_]).subs(zero)
    return np.array(A.evalf(30), dtype=float), np.array(B.evalf(30), dtype=float)


def check_c():
    out = {"rows": []}
    ok = True
    k_ratio = (0.7 / (1 - D_PHYS)) ** 4
    parts = {D_: hessian_parts(D_), D_PHYS: hessian_parts(D_PHYS)}
    for d, w, c_used, name in (
        (D_, W1 * 25, C_SCAN, "delta 0.3, W1 x 25"),
        (D_PHYS, W1 * 25, (3.094e-05,), "delta phys, W1 x 25 (row A)"),
        (D_PHYS, W1 * 25 / k_ratio, (), "delta phys, scaled w (row B)"),
    ):
        A, B = parts[d]
        coef = (G_ + d) * (1 - d)
        for c in (0.0,) + tuple(c_used):
            ev = np.linalg.eigvalsh(w * A - c * B)
            scale = abs(ev).max()
            zero = int(np.sum(np.abs(ev) < 1e-10 * scale))
            dbl = [float(x) for x in ev if c > 0 and abs(x - coef * c) < 1e-6 * coef * c]
            massive = [
                float(x)
                for x in ev
                if abs(x) >= 1e-10 * scale and not (c > 0 and abs(x - coef * c) < 1e-6 * coef * c)
            ]
            out["rows"].append(
                {
                    "case": name,
                    "c": c,
                    "zero_modes": zero,
                    "doublet": dbl,
                    "doublet_closed_form": coef * c,
                    "massive": massive,
                    "min_massive_over_w": min(massive) / w,
                }
            )
            ok &= (zero == 5 and len(dbl) == 2 and min(massive) > 0) if c > 0 else zero == 7
        # c_crit: the first c where an eigenvalue falls below minus 1e-6 of the softest massive mode at c = 0
        soft = min(x for x in np.linalg.eigvalsh(w * A) if abs(x) > 1e-10 * abs(w * A).max())
        cs = np.geomspace(1e-3 * soft, 1e4 * soft, 8000)
        neg = [c_ for c_ in cs if np.linalg.eigvalsh(w * A - c_ * B)[0] < -1e-6 * soft]
        out[f"softest_massive_at_c0[{name}]"] = float(soft)
        out[f"c_crit[{name}]"] = float(neg[0]) if neg else None
        out[f"c_crit_over_c_max_used[{name}]"] = (
            float(neg[0] / max(c_used)) if (neg and c_used) else None
        )
    A, B = parts[D_PHYS]
    wB = W1 * 25 / k_ratio
    out["row_B_at_the_planned_c_3.094e-05"] = {
        "lowest_eigenvalue": float(np.linalg.eigvalsh(wB * A - 3.094e-05 * B)[0]),
        "reading": "negative: on the scaled w the vacuum is unstable at the planned c, so row B runs at c = 0 only "
        "(deviation 1); the largest stable c there cuts the halo at about the half-box",
    }
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "the zero-mode count differed from 5 (7 at c = 0), the doublet left (g + delta)(1 - delta) c, or a massive mode went negative inside the scan"
    )
    print("c:", {k: v for k, v in out.items() if k != "rows"})
    for r in out["rows"]:
        print(
            "   ",
            r["case"],
            r["c"],
            r["zero_modes"],
            [f"{x:.4g}" for x in r["doublet"]],
            [f"{x:.4g}" for x in r["massive"]],
        )
    return out


def check_d():
    rng = np.random.default_rng(11)
    out = {}
    ok = True
    for d, w, c_list, name in (
        (D_, W1 * 25, (0.0,) + C_SCAN, "delta 0.3, W1 x 25"),
        (D_PHYS, W1 * 25, (3.094e-05,), "row A"),
        (D_PHYS, W1 * 25 / ((0.7 / (1 - D_PHYS)) ** 4), (0.0,), "row B"),
    ):
        spec0 = np.array([-G_, 1.0, d, d])
        a = [float(v) for v in a_of(G_, d)]

        def v_spec(l, c_):
            r = [np.sum(l ** (k + 1)) - np.sum(spec0 ** (k + 1)) for k in range(4)]
            return w * sum(x * x for x in r) - c_ * sum(a[k] * r[k] for k in range(4))

        for c_ in c_list:
            best, best_l, below = np.inf, None, 0
            for _ in range(600):
                l0 = rng.uniform(-12, 12, size=4)
                res = minimize(v_spec, l0, args=(c_,), method="BFGS", options={"gtol": 1e-12})
                if res.fun < -1e-10 * w:
                    below += 1
                if res.fun < best:
                    best, best_l = float(res.fun), sorted(float(x) for x in res.x)
            # full symmetric 4x4 matrices: far random starts (any signature sector) and starts near the vacuum;
            # every end state is classified by the eta-norm of its eigenvectors (N is eta-self-adjoint, so a
            # complex pair spans a (1,1) plane and takes the timelike direction out of the real spectrum)
            sectors = {
                "vacuum sector (real spectrum, timelike eigenvalue near -g)": [],
                "other": [],
            }
            M0 = np.diag([G_, 1.0, d, d])
            for k_ in range(160):
                v0 = rng.uniform(-10, 10, size=10) if k_ < 80 else None
                start = (
                    sym_of(v0)
                    if k_ < 80
                    else M0 + sym_of(rng.normal(size=10) * (0.3 if k_ < 120 else 1.5))
                )
                res = minimize(
                    lambda v: v_of_M(sym_of(v), w, c_, d, a),
                    np.array([start[i, j] * (1.0 if i == j else np.sqrt(2.0)) for i, j in IDX]),
                    method="BFGS",
                    options={"gtol": 1e-12},
                )
                lam, vec = np.linalg.eig(sym_of(res.x) @ ETA)
                real = np.abs(lam.imag) < 1e-7
                tl = [
                    float(lam[i].real)
                    for i in range(4)
                    if real[i] and (vec[:, i].real @ ETA @ vec[:, i].real) < 0
                ]
                key = (
                    "vacuum sector (real spectrum, timelike eigenvalue near -g)"
                    if (real.all() and len(tl) == 1 and abs(tl[0] + G_) < 1.0)
                    else "other"
                )
                sectors[key].append(float(res.fun))
            closed = -c_ * c_ * sum(x * x for x in a) / (4.0 * w)
            out[f"{name} c={c_:g}"] = {
                "lowest_on_real_spectra": best,
                "spectrum_there": best_l,
                "real_spectra_minima_below_zero": below,
                "matrix_minima_in_vacuum_sector": len(
                    sectors["vacuum sector (real spectrum, timelike eigenvalue near -g)"]
                ),
                "lowest_in_vacuum_sector": min(
                    sectors["vacuum sector (real spectrum, timelike eigenvalue near -g)"],
                    default=None,
                ),
                "matrix_minima_in_other_sectors": len(sectors["other"]),
                "lowest_in_other_sectors": min(sectors["other"], default=None),
                "unconstrained_minimum_closed_form_-c^2_sum_a^2_/_4w": closed,
            }
            lv = out[f"{name} c={c_:g}"]["lowest_in_vacuum_sector"]
            ok &= below == 0 and (lv is None or lv > -1e-9 * w)
    out["reading"] = (
        "V4 - c L is lower than the vacuum by c^2 sum a_p^2 / 4w on matrices whose spectrum carries a complex "
        "pair next to delta (the timelike direction inside the pair): a different signature sector, which "
        "at c = 0 is exactly degenerate with the vacuum (the same characteristic polynomial) and is reached "
        "only by moving the timelike eigenvalue from -g to delta. The block-diagonal sector of the "
        "instrument has a real spectrum, where the vacuum is the global minimum"
    )
    out["PASS"] = bool(ok)
    out["fails_if"] = (
        "a minimum with V < 0 on real spectra, or inside the vacuum's signature sector on full matrices"
    )
    print("d:", json.dumps(out, indent=1))
    return out


def check_e():
    w = W1 * 25
    out = {}
    a_iso = (1 + 2 * D_) / 3
    spec0 = np.array([-G_, 1.0, D_, D_])

    def v4_spec(l):
        return w * sum((np.sum(np.array(l) ** p) - np.sum(spec0**p)) ** 2 for p in range(1, 5))

    cands = {
        "isotropic block at fixed trace, M_00 = g": v4_spec([-G_, a_iso, a_iso, a_iso]),
        "isotropic block at fixed trace, M_00 solved": min(
            v4_spec([u, a_iso, a_iso, a_iso]) for u in np.linspace(-8.2, -7.8, 40001)
        ),
        "zero spatial block, M_00 = g": v4_spec([-G_, 0.0, 0.0, 0.0]),
    }
    coef = (G_ + D_) * (1 - D_)
    for k, V0 in cands.items():
        M2 = 3.14e-3 * V0
        out[k] = {
            "V_0": float(V0),
            "c if V = (M^2 / 2) eps_full^2, eps_full = 2 eps": float(2 * M2 / coef),
            "c if V = M^2 eps_full^2": float(4 * M2 / coef),
            "c if V = (M^2 / 2) eps^2 (our eps)": float(M2 / (2 * coef)),
        }
    out["note"] = (
        "author-gated: which state defines V_0, and the convention of the split amplitude"
    )
    print("e:", json.dumps(out, indent=1))
    return out


def halo_profile(beta, r_max_scaled=7.0):
    """the decaying solution of eps'' = eps / r^2 + beta r^2 eps, integrated inward from the Gaussian tail."""
    Rs = beta ** (-0.25)
    r1 = r_max_scaled * Rs
    y1 = [1e-30, -np.sqrt(beta) * r1 * 1e-30]
    sol = solve_ivp(
        lambda r, y: [y[1], y[0] / r**2 + beta * r**2 * y[0]],
        (r1, 0.5),
        y1,
        rtol=1e-10,
        atol=1e-300,
        dense_output=True,
    )
    return sol.sol


def reader_radius(sol, r_ref=4.5, p_tail=0.6180339887498949):
    f = lambda r: sol(r)[0] * r**p_tail  # noqa: E731
    y0 = f(r_ref)
    rr = np.linspace(r_ref, 60.0, 20000)
    yy = np.array([f(r) for r in rr])
    k = np.where(yy < 0.5 * y0)[0]
    return float(rr[k[0]]) if len(k) else None


def check_f():
    K = 8 * (1 - D_) ** 2
    out = {"K": K, "rows": []}
    for c in C_SCAN:
        m2 = (G_ + D_) * (1 - D_) * c
        row = {"c": c, "R_scale (2K/m2)^(1/4)": float((2 * K / m2) ** 0.25)}
        for f_ang in (0.5, 1.0, 2.0):
            beta = f_ang * m2 / (2 * K)
            row[f"reader_radius_f{f_ang:g}"] = reader_radius(halo_profile(beta))
        out["rows"].append(row)

    # the exact exponent: a radius defined on the scaled solution (the maximum of eps r^0.618 has none; use the
    # point where the log slope of eps reaches -2)
    def slope_radius(beta):
        sol = halo_profile(beta)
        g_ = lambda r: r * sol(r)[1] / sol(r)[0] + 2.0  # noqa: E731
        return brentq(g_, 0.6, 6.5 * beta**-0.25)

    bs = [1e-5, 1e-4, 1e-3]
    rs = [slope_radius(b) for b in bs]
    expo = np.polyfit(np.log(bs), np.log(rs), 1)[0]
    rd = [r["reader_radius_f1"] for r in out["rows"]]
    out["intrinsic_radius_exponent_in_beta"] = float(expo)
    out["reader_log_slope_vs_c_f1_all_five"] = float(np.polyfit(np.log(C_SCAN), np.log(rd), 1)[0])
    out["reader_log_slope_vs_c_f1_lowest_three"] = float(
        np.polyfit(np.log(C_SCAN[:3]), np.log(rd[:3]), 1)[0]
    )
    out["PASS"] = bool(abs(expo + 0.25) < 1e-3)
    out["fails_if"] = "the intrinsic radius of the decaying solution did not scale as beta^(-1/4)"
    out["note"] = (
        "the reader's fixed r_ref 4.5 bends its slope away from -1/4 at large c, where the cutoff reaches "
        "the reference radius; the HALO_QUARTER window [-0.3, -0.2] is read against this table"
    )
    print("f:", json.dumps(out, indent=1))
    return out


def main():
    res = {
        "a": check_a(),
        "b": check_b(),
        "c": check_c(),
        "d": check_d(),
        "e": check_e(),
        "f": check_f(),
    }
    res["PASS"] = all(res[k].get("PASS", True) for k in res)
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1)
    print("R23-0 PASS:", res["PASS"])


if __name__ == "__main__":
    main()
