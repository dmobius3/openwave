"""M5.32 R22-0 adversarial audit: an independent attempt to refute the claims of
R22-0 (a, b, c), R22-3 and the instrument parts of R22-1 / R22-2.

RULES OF THIS FILE
------------------
Every quantity is re-implemented here (own one-sided and central stencils, own
eta-commutator, own sympy derivations); the production modules are imported
read-only and used only as the thing being compared against. A check prints
    C<id>.<n> PASS|FAIL <what> <numbers>
PASS means the claim survived the attack, FAIL means a concrete counter-number
was found. A PASS that still changes how a claim should be worded carries a
"note". All checks land in ../data/m5_32_r22_0_audit.json.

METHOD PER CLAIM
----------------
C1  sympy, general n on S^2 (sines and cosines of the two angles as algebraic symbols,
    polynomial remainder modulo s^2 + c^2 = 1), tangent jets a_i = alpha_i e_theta + beta_i e_phi:
    4 sum_{i<j} <F_ij, F_ij>_eta - 8 (1 - delta)^4 |E_top|^2 has remainder 0.
    With a radial part a_i -> a_i + gamma_i n the quartic gains the positive term
    32 (1 - delta)^4 sum_{i<j} |gamma_j t_i - gamma_i t_j|^2 (derived and verified).
    The author's dens() contraction is evaluated on the same jets (factor 4).
C2  the author's faber_identity.py is run unchanged in a subprocess; then the same
    field (own analytic derivatives through the complex w) is read at n 40, 80,
    160, 320 under nine norms, three margins and two references, and the table is
    searched for the posted pair (6.5 percent, 1.7 percent).
C3  a different smooth director field (spherical angles a(x), b(x), not the
    stereographic plane waves), own stencils, analytic |E_top|^2 = sin a |grad a x
    grad b|. Order of the fwd / bwd averaged density, the ratio to the central
    stencil, and the analytic leading truncation term: with D^{+/-} P = P_i +/-
    (h/2) P_ii + ..., F^{+/-} = F^c +/- (h/2) F_1 + (h^2/4) [P_ii, P_jj], so the
    average keeps h^2 ((1/4) |F_1|^2 + (1/2) F^c . [P_ii, P_jj]) which the central
    stencil does not have. The prediction is tested pointwise. The ratio is then
    read on three fields (own, the audited one, the hedgehog).
C4  the hedgehog shell sums with own stencil code.
C5  sympy at a point (R = 1 by covariance; the covariance itself is checked with
    random rotations in floating point): closed form, projector form, Maurer-Cartan,
    the uniaxial limit and the biaxial weights.
C6  own reader of the R21 field (own eigen-decomposition, own stencils, stack and
    central): shell ratio, axis share, <e r^4>, plus the radial trend of <e r^4>
    and the local value of G_large^2, which decide how "0.49" may be worded.
C7  sympy: the split series; the exact Hessian of V4 and Vspec on the vacuum as
    2 J^T J (no finite differences), its rank for symbolic (g, delta), the null
    directions, the nonzero eigenvalues at g 8, delta 0.3.
C8  all 48 signed permutations applied to a generic smooth 4x4 field; production
    B3.e_parts with the "sym" and the "cen" stencil; count of exact invariances.
C9  own retyping of the omega_E route; the control; the branch translation tested
    on a random block-diagonal field through production B3.e_parts (s = -1 against
    s = +1), curvature and V4 separately, with a negative control (M_0i not zero),
    and V4 on both branches for the two R21 fields themselves; the two R21 rows; then the
    sensitivity to the envelope radius, to the boost axis and to the lattice
    spacing (the unit Frobenius norm of the tangent is a cell sum, not an integral).
C10 energy_grad_n against production B3.e_parts / B3.grad through the chain rule,
    and a central finite-difference test, own field, box and seed.
C11 sympy: delta E_top = curl alpha with alpha_j = n . (delta n x d_j n), alpha . E = 0,
    E . grad n = 0; hence the wall term 2 oint alpha . (E x nu) and the natural
    condition E x nu = 0. G_D of the cube from a mixed eigenfunction series (sine
    modes across, closed sinh form along the pair axis), and from an own DST solve
    on a finer grid; the wall position of a cell-centered box from the effective cell
    counts of each field component.
X   findings outside the listed claims (docstring statements checked against the code).
C12 sympy: the flux-excluding sphere in a uniform field (dipole, outside and inside
    energies), the general multipole l, and the full multipole sum for a point
    charge at distance d against the dipole-only estimate.

Regenerate: python m5_32_r22_0_audit.py [--only C3,C9]   (about 6 min, under 3 GB)
"""

import importlib.util
import itertools
import json
import os
import subprocess
import sys
import time

import numpy as np
import sympy as sp

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
THEORY = os.path.join(HERE, "..", "..", "theory")
OUT_JSON = os.path.join(DATA, "m5_32_r22_0_audit.json")
AUTHOR_SCRIPT = os.path.join(THEORY, "r22_author_scripts", "files", "faber_identity.py")
EXT_008 = os.path.join(
    THEORY, "r22_mikulski_008", "reports", "008-i1-squared-clock", "results", "L_ladder"
)
R21_DIR = os.path.join(DATA, "m5_32_r21_1")
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
G8, DELTA = 8.0, 0.3
T0 = time.time()
CHECKS = []


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


R20 = _load("m5_32_r20_1_axes", "m5_32_r20_1_axes.py")
B3, R0 = R20.B3, R20.R0
PAIR22 = _load("m5_32_r22_2_pair", "m5_32_r22_2_pair.py")


def _js(v):
    if isinstance(v, dict):
        return {str(k): _js(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_js(x) for x in v]
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (sp.Basic,)):
        return str(v)
    return v


def record(cid, claim, method, value, expected, ok, note=""):
    verdict = "PASS" if ok else "FAIL"
    CHECKS.append(
        {
            "id": cid,
            "claim": claim,
            "method": method,
            "value": _js(value),
            "expected": _js(expected),
            "verdict": verdict,
            "note": note,
        }
    )
    tail = f" | note: {note}" if note else ""
    print(
        f"{cid} {verdict} {claim}: value {_js(value)} expected {_js(expected)}{tail}  "
        f"[{time.time() - T0:.0f}s]",
        flush=True,
    )


# ================= own lattice layer =================
def own_diff(f, ax, h, kind):
    """one-sided (fwd / bwd, zero where the neighbor is missing) or central (one-sided ends)."""
    out = np.zeros_like(f)
    lo = [slice(None)] * f.ndim
    hi = [slice(None)] * f.ndim
    if kind in ("fwd", "bwd"):
        lo[ax], hi[ax] = slice(0, -1), slice(1, None)
        tgt = lo if kind == "fwd" else hi
        out[tuple(tgt)] = (f[tuple(hi)] - f[tuple(lo)]) / h
        return out
    mid = [slice(None)] * f.ndim
    lo[ax], hi[ax], mid[ax] = slice(0, -2), slice(2, None), slice(1, -1)
    out[tuple(mid)] = (f[tuple(hi)] - f[tuple(lo)]) / (2 * h)
    for end, a, b in ((0, 1, 0), (-1, -1, -2)):
        e, p, q = [slice(None)] * f.ndim, [slice(None)] * f.ndim, [slice(None)] * f.ndim
        e[ax], p[ax], q[ax] = end, a, b
        out[tuple(e)] = (f[tuple(p)] - f[tuple(q)]) / h
    return out


def own_quartic(A):
    """4 sum_{i<j} tr(eta F eta F^T), F = A_i eta A_j - A_j eta A_i (eta = 1 for 3x3 blocks)."""
    eta = ETA if A[0].shape[-1] == 4 else np.eye(3)
    e = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = A[i] @ eta @ A[j] - A[j] @ eta @ A[i]
            e = e + 4.0 * np.einsum("...ab,...ab->...", eta @ F @ eta, F)
    return e


def own_density(M, h, stencil):
    kinds = ("fwd", "bwd") if stencil == "sym" else (stencil,)
    return sum(own_quartic([own_diff(M, ax, h, k) for ax in range(3)]) for k in kinds) / len(kinds)


def outer(u, v):
    return u[..., :, None] * v[..., None, :]


def grid(n_, L):
    h = L / n_
    x = (np.arange(n_) - (n_ - 1) / 2.0) * h
    return np.meshgrid(x, x, x, indexing="ij"), h


def rel_l2(d, r):
    return float(np.sqrt(np.sum(d * d) / np.sum(r * r)))


# ================= C1 =================
def c1():
    # sin / cos of the two angles as algebraic symbols, reduced modulo s^2 + c^2 = 1 (a unique
    # remainder: the two relations have coprime leading monomials, so they are a Groebner basis);
    # no trigonometric simplifier involved
    st, ct, sf, cf, dl = sp.symbols("s_theta c_theta s_phi c_phi delta", real=True)
    al, be, ga = (
        sp.symbols("alpha0:3", real=True),
        sp.symbols("beta0:3", real=True),
        sp.symbols("gamma0:3", real=True),
    )
    n = sp.Matrix([st * cf, st * sf, ct])
    e1 = sp.Matrix([ct * cf, ct * sf, -st])
    e2 = sp.Matrix([-sf, cf, 0])
    eta = sp.diag(-1, 1, 1, 1)

    def rem(expr):
        return sp.reduced(
            sp.expand(expr), [st**2 + ct**2 - 1, sf**2 + cf**2 - 1], st, sf, ct, cf
        )[1]

    frame_ok = all(
        rem(v) == 0
        for v in (n.dot(n) - 1, e1.dot(e1) - 1, e2.dot(e2) - 1, n.dot(e1), n.dot(e2), e1.dot(e2))
    )

    def total(jets):
        dM = []
        for a in jets:
            D = sp.zeros(4, 4)
            D[1:, 1:] = (1 - dl) * (a * n.T + n * a.T)
            dM.append(D)
        tot = 0
        for i in range(3):
            for j in range(i + 1, 3):
                F = dM[i] * eta * dM[j] - dM[j] * eta * dM[i]
                tot += 4 * (eta * F * eta * F.T).trace()
        return tot, dM

    def etop2(jets):
        return sum(n.dot(jets[i].cross(jets[j])) ** 2 for i, j in ((1, 2), (2, 0), (0, 1)))

    t = [al[i] * e1 + be[i] * e2 for i in range(3)]
    tot, dM = total(t)
    res = rem(tot - 8 * (1 - dl) ** 4 * etop2(t))
    record(
        "C1.1",
        "stack quartic = 8 (1-delta)^4 |E_top|^2 for tangent jets, symbolic",
        "sympy, general n on S^2 (polynomial remainder modulo s^2 + c^2 = 1), orthonormal frame "
        "verified",
        str(res),
        "0",
        res == 0 and frame_ok,
    )
    # the author's contraction 2 (T1 - T2) with D = eta dM on the same jets
    D = [eta * X for X in dM]
    Q = D[0] * D[0] + D[1] * D[1] + D[2] * D[2]
    T2 = sum((D[i] * D[j] * D[i] * D[j]).trace() for i in range(3) for j in range(3))
    res_a = rem(2 * ((Q * Q).trace() - T2) - 4 * (1 - dl) ** 4 * etop2(t))
    record(
        "C1.2",
        "the author's dens() contraction = 4 (1-delta)^4 |E_top|^2, so the stack sits at c4 = 2",
        "sympy on the posted formula 2 (T1 - T2)",
        str(res_a),
        "0",
        res_a == 0,
    )
    # radial component present
    a = [t[i] + ga[i] * n for i in range(3)]
    tot_r, _ = total(a)
    extra = sum(
        ((ga[j] * t[i] - ga[i] * t[j]).dot(ga[j] * t[i] - ga[i] * t[j]))
        for i in range(3)
        for j in range(i + 1, 3)
    )
    res_r = rem(tot_r - 8 * (1 - dl) ** 4 * (etop2(a) + 4 * extra))
    res_e = rem(etop2(a) - etop2(t))
    res_plain = rem(tot_r - 8 * (1 - dl) ** 4 * etop2(a))
    record(
        "C1.3",
        "with n . d_i n = gamma_i the quartic is 8 (1-delta)^4 (|E_top|^2 + 4 sum_{i<j} |gamma_j "
        "t_i - gamma_i t_j|^2)",
        "sympy; E_top itself does not see gamma; the plain identity must NOT hold there (negative "
        "control)",
        [str(res_r), str(res_e), "plain identity residual nonzero: " + str(res_plain != 0)],
        ["0", "0", "True"],
        res_r == 0 and res_e == 0 and res_plain != 0,
        note="the identity fails upward (positive excess) off the constraint; a one-sided lattice "
        "difference has "
        "gamma_i = -/+ (h/2) |d_i n|^2, one part of the C3 excess",
    )
    # floating point against the production algebra
    rng = np.random.default_rng(20260919)
    worst, fac = 0.0, []
    for _ in range(500):
        nn = rng.normal(size=3)
        nn /= np.linalg.norm(nn)
        d_ = rng.uniform(0.05, 0.95)
        jets = [v - nn * (nn @ v) for v in rng.normal(size=(3, 3))]
        A = []
        for v in jets:
            X = np.zeros((4, 4))
            X[1:, 1:] = (1 - d_) * (np.outer(v, nn) + np.outer(nn, v))
            A.append(X)
        own = float(own_quartic(A))
        prod = sum(
            4.0 * float(B3.inner_eta(B3.comm_eta(A[i], A[j]), B3.comm_eta(A[i], A[j])))
            for i in range(3)
            for j in range(i + 1, 3)
        )
        e2 = sum((nn @ np.cross(jets[i], jets[j])) ** 2 for i, j in ((1, 2), (2, 0), (0, 1)))
        worst = max(
            worst, abs(own - prod) / abs(prod), abs(prod - 8 * (1 - d_) ** 4 * e2) / abs(prod)
        )
        fac.append(prod / ((1 - d_) ** 4 * e2))
    record(
        "C1.4",
        "production comm_eta / inner_eta give the factor 8 on random tangent jets (own seed)",
        "own commutator against B3, 500 jets",
        {"max_rel": worst, "factor_min": min(fac), "factor_max": max(fac)},
        {"max_rel": "< 1e-10", "factor": 8},
        worst < 1e-10,
    )
    slope = 64 * np.pi * (1 - DELTA) ** 4
    record(
        "C1.5",
        "pair slope 64 pi (1-delta)^4 at delta 0.3",
        "arithmetic",
        slope,
        48.27,
        abs(slope - 48.27) < 0.01,
    )


# ================= C2 =================
def author_field(X, Y, Z, kscale=1.0):
    """The posted field with own analytic derivatives:
    rho_ij = 4 Im(conj(d_i w) d_j w) / (1 + |w|^2)^2."""
    rng = np.random.default_rng(0)
    md = [(rng.normal() + 1j * rng.normal(), rng.normal(size=3) * 0.6 * kscale) for _ in range(6)]
    w = sum(c * np.exp(1j * (k[0] * X + k[1] * Y + k[2] * Z)) for c, k in md)
    dw = [
        sum(1j * k[i] * c * np.exp(1j * (k[0] * X + k[1] * Y + k[2] * Z)) for c, k in md)
        for i in range(3)
    ]
    s = 1.0 + np.abs(w) ** 2
    n = np.stack([2 * w.real, 2 * w.imag, 2.0 - s], axis=-1) / s[..., None]
    rho = lambda i, j: 4.0 * np.imag(np.conj(dw[i]) * dw[j]) / s**2  # noqa: E731
    E = np.stack([rho(1, 2), rho(2, 0), rho(0, 1)], axis=-1)
    gmax = max(float((2 * np.abs(d) / s).max()) for d in dw)
    g4 = (sum(4 * np.abs(d) ** 2 for d in dw) / s**2) ** 2
    return n, E, gmax, g4


def c2_arrays(n_, L=10.0, slab=32):
    slab = slab if n_ <= 160 else 8
    """lhs (central quartic, the author's normalization), rhs from lattice E_top, rhs analytic;
    interior cells (one layer dropped, where the posted script's periodic roll wraps). Built in
    x slabs to bound memory."""
    h = L / n_
    x = (np.arange(n_) - n_ / 2 + 0.5) * h
    k4 = 4 * (1 - DELTA) ** 4
    lhs = np.zeros((n_ - 2, n_, n_))
    rl = np.zeros_like(lhs)
    ra = np.zeros_like(lhs)
    gmax = 0.0
    for i0 in range(1, n_ - 1, slab):
        i1 = min(i0 + slab, n_ - 1)
        X, Y, Z = np.meshgrid(x[i0 - 1 : i1 + 1], x, x, indexing="ij")
        n, E, gm, _ = author_field(X, Y, Z)
        gmax = max(gmax, gm)
        P = outer(n, n) * (1 - DELTA)
        A = [own_diff(P, ax, h, "cen") for ax in range(3)]
        dn = [own_diff(n, ax, h, "cen") for ax in range(3)]
        El = np.stack(
            [
                np.einsum("...a,...a->...", n, np.cross(dn[i], dn[j]))
                for i, j in ((1, 2), (2, 0), (0, 1))
            ],
            -1,
        )
        lhs[i0 - 1 : i1 - 1] = (own_quartic(A) / 2.0)[1:-1]
        rl[i0 - 1 : i1 - 1] = (k4 * np.sum(El * El, -1))[1:-1]
        ra[i0 - 1 : i1 - 1] = (k4 * np.sum(E * E, -1))[1:-1]
    return lhs, rl, ra, gmax * h


def c2():
    env = dict(os.environ, OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    out = subprocess.run(
        [sys.executable, AUTHOR_SCRIPT], capture_output=True, text=True, env=env, timeout=300
    ).stdout
    line = [ln for ln in out.splitlines() if "max|lhs-rhs|" in ln]
    printed = float(line[0].split("=")[1]) if line else None
    record(
        "C2.1",
        "the posted script, run unchanged, prints 7.36e-01 for item (1)",
        "subprocess, stdout parsed",
        printed,
        0.736,
        printed is not None and abs(printed - 0.736) < 5e-3,
    )
    table, gh = {}, {}
    for n_ in (40, 80, 160, 320):
        lhs, rl, ra, gh[n_] = c2_arrays(n_)
        for mg in (2, 4, 8):
            k = mg - 1  # one layer is already dropped along x
            I = (slice(k, -k) if k else slice(None), slice(mg, -mg), slice(mg, -mg))
            for ref_name, ref in (("lattice_rhs", rl), ("analytic_rhs", ra)):
                d = lhs[I] - ref[I]  # memory-lean: two work arrays, reused in place
                ar = np.abs(ref[I])
                row = {"sum_rel": float(abs(d.sum() / ref[I].sum()))}
                l2n, rms = float(np.vdot(d, d)), float(np.sqrt(np.vdot(d, d) / d.size))
                np.abs(d, out=d)
                amax = float(ar.max())
                row.update(
                    {
                        "max_over_max": float(d.max() / amax),
                        "L2": float(np.sqrt(l2n / np.vdot(ar, ar))),
                        "L1": float(d.sum() / ar.sum()),
                        "rms_over_max": rms / amax,
                        "mean_over_max": float(d.mean() / amax),
                        "weighted_rel": float(np.vdot(d, ar) / np.vdot(ar, ar)),
                    }
                )
                np.maximum(ar, 1e-300, out=ar)
                np.divide(d, ar, out=d)
                row.update(
                    {
                        "mean_pointwise_rel": float(d.mean()),
                        "median_pointwise_rel": float(np.median(d)),
                    }
                )
                table[(n_, mg, ref_name)] = row
                del d, ar
        del lhs, rl, ra
    own_max = table[(40, 2, "lattice_rhs")]["max_over_max"]
    record(
        "C2.2",
        "own re-implementation reproduces the printed max norm at n 40",
        "own central stencil and E_top",
        own_max,
        0.7359,
        abs(own_max - 0.7359) < 2e-3,
    )
    record(
        "C2.3",
        "the posted field is under-resolved at n 40: max |grad n| h about 1.3",
        "analytic |d_i n| = 2 |d_i w| / (1+|w|^2)",
        {"n40": gh[40], "n80": gh[80], "n160": gh[160], "n320": gh[320]},
        {"n40": 1.29},
        abs(gh[40] - 1.29) < 0.02,
    )
    hits_4080, near = [], []
    for (n_, mg, ref_name), row in table.items():
        if (2 * n_, mg, ref_name) not in table:
            continue
        nxt = table[(2 * n_, mg, ref_name)]
        for key in row:
            a, b = row[key], nxt[key]
            if abs(a - 0.065) < 0.008 and abs(b - 0.017) < 0.004:
                (hits_4080 if n_ == 40 else near).append(
                    {
                        "n": [n_, 2 * n_],
                        "margin": mg,
                        "ref": ref_name,
                        "norm": key,
                        "values": [a, b],
                    }
                )
    record(
        "C2.4",
        "no norm, mask or reference gives (6.5 percent, 1.7 percent) at n 40 and n 80",
        "9 norms x 3 margins x 2 references",
        {
            "hits_at_n40_n80": hits_4080,
            "n40_n80_rows_margin2_lattice": [
                table[(40, 2, "lattice_rhs")],
                table[(80, 2, "lattice_rhs")],
            ],
        },
        "no hit",
        len(hits_4080) == 0,
        note="the pair does come out of the same field on finer grids: "
        + "; ".join(
            f"n {q['n'][0]} to {q['n'][1]} {q['norm']} vs {q['ref']}: {q['values'][0]:.3f} to "
            f"{q['values'][1]:.3f}"
            for q in near
            if q["margin"] == 2
        ),
    )
    return {"table": {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in table.items()}, "near_hits": near}


# ================= C3 / C4 =================
K_OWN = np.array(
    [[0.31, -0.22, 0.17], [-0.12, 0.28, 0.35], [0.26, 0.19, -0.30], [0.20, -0.33, 0.11]]
)


def own_field(X, Y, Z):
    """n = (sin a cos b, sin a sin b, cos a), a and b sums of sines: analytic n, d_i n, d_ii n,
    E_top."""
    k = K_OWN
    p = [k[m, 0] * X + k[m, 1] * Y + k[m, 2] * Z for m in range(4)]
    a = 1.2 + 0.5 * np.sin(p[0] + 0.4) + 0.3 * np.cos(p[1] - 0.2)
    b = p[2] + 0.7 * np.sin(p[3] + 1.0)
    da = [
        0.5 * np.cos(p[0] + 0.4) * k[0, i] - 0.3 * np.sin(p[1] - 0.2) * k[1, i] for i in range(3)
    ]
    db = [k[2, i] + 0.7 * np.cos(p[3] + 1.0) * k[3, i] for i in range(3)]
    daa = [
        -0.5 * np.sin(p[0] + 0.4) * k[0, i] ** 2 - 0.3 * np.cos(p[1] - 0.2) * k[1, i] ** 2
        for i in range(3)
    ]
    dbb = [-0.7 * np.sin(p[3] + 1.0) * k[3, i] ** 2 for i in range(3)]
    sa, ca, sb, cb = np.sin(a), np.cos(a), np.sin(b), np.cos(b)
    z = np.zeros_like(sa)
    n = np.stack([sa * cb, sa * sb, ca], -1)
    na, nb = np.stack([ca * cb, ca * sb, -sa], -1), np.stack([-sa * sb, sa * cb, z], -1)
    nab, nbb = np.stack([-ca * sb, ca * cb, z], -1), np.stack([-sa * cb, -sa * sb, z], -1)
    ni = [na * da[i][..., None] + nb * db[i][..., None] for i in range(3)]
    nii = [
        -n * (da[i] ** 2)[..., None]
        + 2 * nab * (da[i] * db[i])[..., None]
        + nbb * (db[i] ** 2)[..., None]
        + na * daa[i][..., None]
        + nb * dbb[i][..., None]
        for i in range(3)
    ]
    E = sa[..., None] * np.cross(np.stack(da, -1), np.stack(db, -1))
    return n, ni, nii, E


def c3():
    rows = []
    for n_ in (32, 64, 128):
        (X, Y, Z), h = grid(n_, 10.0)
        n, ni, nii, E = own_field(X, Y, Z)
        ref = 8 * (1 - DELTA) ** 4 * np.sum(E * E, -1)
        M3 = DELTA * np.eye(3) + (1 - DELTA) * outer(n, n)
        es, ec = own_density(M3, h, "sym"), own_density(M3, h, "cen")
        k = max(2, int(round(0.08 * n_)))
        I = (slice(k, -k),) * 3
        row = {
            "n": n_,
            "h": h,
            "max_grad_n_h": float(max(np.linalg.norm(q, axis=-1).max() for q in ni) * h),
            "stack_L2": rel_l2((es - ref)[I], ref[I]),
            "central_L2": rel_l2((ec - ref)[I], ref[I]),
            "stack_sum_rel": float((es - ref)[I].sum() / ref[I].sum()),
            "central_sum_rel": float((ec - ref)[I].sum() / ref[I].sum()),
        }
        if n_ == 32:  # own stencil against the production density
            cfg = B3.base_cfg(s=-1.0, g=G8, n=n_, L=10.0, delta=DELTA)
            ep = (
                R20.density(
                    B3.embed34(M3, cfg), cfg, ("v4", R0.roots_of(cfg, degenerate=True), 0.0)
                )
                / h**3
            )
            row["own_vs_production_density_max_rel"] = float(
                np.abs(ep - es).max() / np.abs(es).max()
            )
        if n_ == 128:  # the analytic leading truncation term (central 64^3 block)
            S = (slice(32, 96),) * 3
            diff = (es - ec)[S]
            nS, niS, niiS = n[S], [q[S] for q in ni], [q[S] for q in nii]
            e2S = np.sum(E * E, -1)[S]
            del es, ec, M3, n, nii, E
            P1 = [(1 - DELTA) * (outer(niS[i], nS) + outer(nS, niS[i])) for i in range(3)]
            P2 = [
                (1 - DELTA) * (outer(niiS[i], nS) + 2 * outer(niS[i], niS[i]) + outer(nS, niiS[i]))
                for i in range(3)
            ]
            pred, pos, gam = 0.0, 0.0, 0.0
            for i in range(3):
                for j in range(i + 1, 3):
                    F0 = P1[i] @ P1[j] - P1[j] @ P1[i]
                    F1 = P2[i] @ P1[j] - P1[j] @ P2[i] + P1[i] @ P2[j] - P2[j] @ P1[i]
                    F2 = P2[i] @ P2[j] - P2[j] @ P2[i]
                    q1 = np.einsum("...ab,...ab->...", F1, F1)
                    pos = pos + h * h * q1
                    pred = pred + h * h * (q1 + 2.0 * np.einsum("...ab,...ab->...", F0, F2))
                    gi, gj = np.sum(niS[i] ** 2, -1), np.sum(niS[j] ** 2, -1)
                    v = gj[..., None] * niS[i] - gi[..., None] * niS[j]
                    gam = gam + 8 * (1 - DELTA) ** 4 * h * h * np.sum(v * v, -1)
            row["leading_term_rel_L2_misfit"] = rel_l2(diff - pred, diff)
            row["positive_definite_share"] = float(pos.sum() / diff.sum())
            row["radial_component_share"] = float(gam.sum() / diff.sum())
            row["mean_E2_over_mean_gradn4"] = float(
                e2S.mean() / (sum(np.sum(q * q, -1) for q in niS) ** 2).mean()
            )
        rows.append(row)
    hs = np.log([r["h"] for r in rows])
    order = float(np.polyfit(hs, np.log([r["stack_L2"] for r in rows]), 1)[0])
    order_c = float(np.polyfit(hs, np.log([r["central_L2"] for r in rows]), 1)[0])
    record(
        "C3.1",
        "own stencil code equals the production fwd / bwd density",
        "R20.density against own_density, n 32",
        rows[0]["own_vs_production_density_max_rel"],
        "< 1e-12",
        rows[0]["own_vs_production_density_max_rel"] < 1e-12,
    )
    record(
        "C3.2",
        "the fwd / bwd averaged density converges at order 2 on a different smooth field",
        "own angle field, L2 fit over n 32, 64, 128",
        {"stack_order": order, "central_order": order_c, "rows": rows},
        "1.5 to 2.5",
        1.5 <= order <= 2.5,
    )
    # the audited field, own code
    (X, Y, Z), h = grid(128, 10.0)
    n, E, gmax, g4 = author_field(X, Y, Z, kscale=0.125)
    ref = 8 * (1 - DELTA) ** 4 * np.sum(E * E, -1)
    M3 = DELTA * np.eye(3) + (1 - DELTA) * outer(n, n)
    I = (slice(10, -10),) * 3
    a_s, a_c = rel_l2((own_density(M3, h, "sym") - ref)[I], ref[I]), rel_l2(
        (own_density(M3, h, "cen") - ref)[I], ref[I]
    )
    aud = {
        "max_grad_n_h": gmax * h,
        "stack_L2": a_s,
        "central_L2": a_c,
        "ratio": a_s / a_c,
        "mean_E2_over_mean_gradn4": float(np.sum(E * E, -1)[I].mean() / g4[I].mean()),
    }
    record(
        "C3.3",
        "on the audited field: L2 16 percent (stack) against 0.12 percent (central) at max |grad "
        "n| h 0.035",
        "own stencils, own analytic reference through w",
        aud,
        {"stack_L2": 0.165, "central_L2": 0.00116},
        abs(a_s - 0.165) < 0.005 and abs(a_c - 0.00116) < 1e-4,
    )
    hed = c4_rows()
    ratios = {
        "own_field_n128": rows[-1]["stack_L2"] / rows[-1]["central_L2"],
        "audited_field_n128": aud["ratio"],
        "hedgehog_h1.5_abs_shell_sum": abs(hed[0]["stack_rel"] / hed[0]["central_rel"]),
        "hedgehog_h0.75_abs_shell_sum": abs(hed[-1]["stack_rel"] / hed[-1]["central_rel"]),
    }
    generic = ratios["own_field_n128"] > 50 and ratios["hedgehog_h1.5_abs_shell_sum"] > 50
    record(
        "C3.4",
        "the error constant of the stack density is roughly 100 to 150 times the central one "
        "(read as generic)",
        "the same ratio on three fields",
        ratios,
        "100 to 150 on every field",
        generic,
        note="the ratio is a property of the audited test field (mean |E_top|^2 / mean |grad n|^4 "
        "= "
        f"{aud['mean_E2_over_mean_gradn4']:.3f} there, {rows[-1]['mean_E2_over_mean_gradn4']:.3f} "
        f"on the own "
        "field); on the hedgehog the central stencil is the WORSE one",
    )
    record(
        "C3.5",
        "analytic leading term of (stack - central): h^2 sum_{i<j} 4 ((1/4)|F_1|^2 + (1/2) F_0 . "
        "[P_ii, P_jj])",
        "pointwise against the measured difference, own field, n 128",
        {
            k: rows[-1][k]
            for k in (
                "leading_term_rel_L2_misfit",
                "positive_definite_share",
                "radial_component_share",
            )
        },
        "misfit < 0.02",
        rows[-1]["leading_term_rel_L2_misfit"] < 0.02,
        note="the O(h) terms +/- (h/2) F_0 . F_1 cancel in the fwd / bwd energy average; the "
        "square (h^2/4)|F_1|^2 "
        "does not, is positive, and does not vanish where E_top does",
    )
    return hed


_C4 = {}


def c4_rows():
    if _C4:
        return _C4["rows"]
    rows = []
    for n_ in (32, 48, 64):
        (X, Y, Z), h = grid(n_, 48.0)
        r = np.sqrt(X * X + Y * Y + Z * Z)
        n = np.stack([X, Y, Z], -1) / r[..., None]
        M3 = DELTA * np.eye(3) + (1 - DELTA) * outer(n, n)
        m = (r > 4.0) & (r < 24.0 - 2 * h)
        an = float((8 * (1 - DELTA) ** 4 / r[m] ** 4).sum())
        rows.append(
            {
                "n": n_,
                "h": h,
                "stack_rel": float(own_density(M3, h, "sym")[m].sum() / an - 1),
                "central_rel": float(own_density(M3, h, "cen")[m].sum() / an - 1),
            }
        )
    _C4["rows"] = rows
    return rows


def c4():
    rows = c4_rows()
    got = [r["stack_rel"] for r in rows]
    exp = [0.050, 0.018, 0.0085]
    record(
        "C4.1",
        "hedgehog L 48, r_c 4: stencil residual +5.0 / +1.8 / +0.85 percent at h 1.5 / 1.0 / 0.75",
        "own stencil, same-cell analytic sum 8 (1-delta)^4 / r^4",
        rows,
        exp,
        all(abs(a - b) < 0.0007 for a, b in zip(got, exp)),
        note="central stencil on the same cells: "
        + ", ".join(f"{r['central_rel']:+.2%}" for r in rows),
    )


# ================= C5 =================
def hat_s(v):
    return sp.Matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])


def c5():
    l = sp.symbols("l1:4", real=True)
    W = [sp.Matrix(sp.symbols(f"w{i}_1:4", real=True)) for i in range(3)]
    Lam = sp.diag(*l)
    eta = sp.diag(-1, 1, 1, 1)
    dM = []
    for Wi in W:
        D = sp.zeros(4, 4)
        D[1:, 1:] = hat_s(Wi) * Lam - Lam * hat_s(Wi)
        dM.append(D)
    lhs = 0
    for i in range(3):
        for j in range(i + 1, 3):
            F = dM[i] * eta * dM[j] - dM[j] * eta * dM[i]
            lhs += 4 * (eta * F * eta * F.T).trace()
    Gs = [sp.prod([l[a] - l[b] for b in range(3) if b != a]) for a in range(3)]
    cur2 = [
        sum(W[i].cross(W[j])[a] ** 2 for i in range(3) for j in range(i + 1, 3)) for a in range(3)
    ]
    res = sp.expand(lhs - 8 * sum(Gs[a] ** 2 * cur2[a] for a in range(3)))
    record(
        "C5.1",
        "frame quartic = 8 sum_a G_a^2 sum_{i<j} rho_ij[e_a]^2, G_a = prod_{b != a} (l_a - l_b)",
        "sympy at a point, d_i M = [w_i, Lambda]",
        str(res),
        "0",
        res == 0,
    )
    ok_p, ok_mc = True, True
    for a in range(3):
        Ea = sp.zeros(3, 3)
        Ea[a, a] = 1
        ea = sp.Matrix([1 if b == a else 0 for b in range(3)])
        dP = [hat_s(Wi) * Ea - Ea * hat_s(Wi) for Wi in W]
        for i in range(3):
            for j in range(i + 1, 3):
                F = dP[i] * dP[j] - dP[j] * dP[i]
                ok_p &= sp.expand((F * F.T).trace() / 2 - W[i].cross(W[j])[a] ** 2) == 0
                mc = ea.dot((hat_s(W[i]) * ea).cross(hat_s(W[j]) * ea))
                ok_mc &= sp.expand(mc - W[i].cross(W[j])[a]) == 0
    record(
        "C5.2",
        "rho_ij[e_a]^2 = (1/2) tr(F F^T), F = [d_i P_a, d_j P_a], and rho_ij[e_a] = (W_i x W_j)^a",
        "sympy, all a and i < j",
        [bool(ok_p), bool(ok_mc)],
        [True, True],
        bool(ok_p and ok_mc),
    )
    d = sp.symbols("delta", real=True)
    uni = [sp.factor(G.subs({l[0]: 1, l[1]: d, l[2]: d})) for G in Gs]
    bia = [float(G.subs({l[0]: 1, l[1]: sp.Rational(3, 10), l[2]: 0}) ** 2) for G in Gs]
    record(
        "C5.3",
        "uniaxial limit G = ((1-delta)^2, 0, 0); biaxial G^2 = (0.49, 0.0441, 0.09) at (1, 0.3, "
        "0)",
        "substitution",
        {"uniaxial_G": [str(u) for u in uni], "biaxial_G2": bia},
        {"uniaxial_G": ["(delta - 1)**2", "0", "0"], "biaxial_G2": [0.49, 0.0441, 0.09]},
        uni[1] == 0
        and uni[2] == 0
        and sp.expand(uni[0] - (1 - d) ** 2) == 0
        and np.allclose(bia, [0.49, 0.0441, 0.09], atol=1e-12),
    )
    # covariance under a random rotation, floating point, production algebra
    rng = np.random.default_rng(919)
    worst = 0.0
    for _ in range(300):
        q, r = np.linalg.qr(rng.normal(size=(3, 3)))
        q = q * np.sign(np.diag(r))
        if np.linalg.det(q) < 0:
            q[:, 0] *= -1
        lam = rng.uniform(-1, 2, size=3)
        Wn = rng.normal(size=(3, 3))
        hat = lambda v: np.array(
            [[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0.0]]
        )  # noqa: E731
        A = []
        for Wi in Wn:
            X = np.zeros((4, 4))
            X[1:, 1:] = q @ (hat(Wi) @ np.diag(lam) - np.diag(lam) @ hat(Wi)) @ q.T
            A.append(X)
        prod = sum(
            4.0 * float(B3.inner_eta(B3.comm_eta(A[i], A[j]), B3.comm_eta(A[i], A[j])))
            for i in range(3)
            for j in range(i + 1, 3)
        )
        Gn = np.array([np.prod([lam[a] - lam[b] for b in range(3) if b != a]) for a in range(3)])
        c2_ = sum(np.cross(Wn[i], Wn[j]) ** 2 for i in range(3) for j in range(i + 1, 3))
        worst = max(worst, abs(prod - 8 * np.sum(Gn**2 * c2_)) / abs(prod))
    record(
        "C5.4",
        "the closed form holds for R not 1 (rotation covariance), production algebra",
        "300 random frames",
        worst,
        "< 1e-10",
        worst < 1e-10,
    )


# ================= C6 =================
def c6():
    out = {}
    for tag in ("S1_Bfree_g8_d0.3_w25_n32", "S1_Bseed_g8_d0.3_w25_n32"):
        M = np.load(os.path.join(R21_DIR, tag + ".npz"))["M"]
        (X, Y, Z), h = grid(M.shape[0], 48.0)
        r = np.sqrt(X * X + Y * Y + Z * Z)
        lam, vec = np.linalg.eigh(M[..., 1:, 1:])
        Gl = np.stack(
            [
                np.prod([lam[..., a] - lam[..., b] for b in range(3) if b != a], axis=0)
                for a in range(3)
            ],
            -1,
        )
        gap = np.minimum(lam[..., 1] - lam[..., 0], lam[..., 2] - lam[..., 1])
        rows = {}
        for st in ("sym", "cen"):
            e = own_density(M, h, st)
            # 8 G^2 rho^2 with rho^2 = (1/2) tr F F^T: own_quartic on P gives
            # 4 sum tr F F^T = 8 sum rho^2
            per = [
                Gl[..., a] ** 2 * own_density(outer(vec[..., :, a], vec[..., :, a]), h, st)
                for a in range(3)
            ]
            pred = per[0] + per[1] + per[2]
            for lo, hi in ((8, 10), (10, 12), (12, 14), (14, 16)):
                ok = (r >= lo) & (r < hi) & (np.abs(Z) < 0.8 * r) & (gap > 0.05)
                rows[f"{st} r{lo}-{hi}"] = {
                    "ratio": float(e[ok].sum() / pred[ok].sum()),
                    "share_large_axis": float(per[2][ok].sum() / pred[ok].sum()),
                    "mean_e_r4": float(np.mean(e[ok] * r[ok] ** 4)),
                    "mean_G_large_sq": float(np.mean(Gl[..., 2][ok] ** 2)),
                    "lam_mean": [float(lam[..., a][ok].mean()) for a in range(3)],
                }
        out[tag] = {
            "rows": rows,
            "max_abs_M0i": float(np.abs(M[..., 0, 1:]).max()),
            "M00_range": [float(M[..., 0, 0].min()), float(M[..., 0, 0].max())],
        }
    s = out["S1_Bfree_g8_d0.3_w25_n32"]["rows"]["sym r12-14"]
    record(
        "C6.1",
        "R21 S1_Bfree w25 n32, shell r 12 to 14 off axis: ratio about 1.10, large-axis share "
        "about 0.97, <e r^4> about 4.2",
        "own reader (own eigh, own fwd / bwd stencil)",
        s,
        {"ratio": 1.10, "share_large_axis": 0.97, "mean_e_r4": 4.2},
        abs(s["ratio"] - 1.10) < 0.02
        and abs(s["share_large_axis"] - 0.97) < 0.01
        and abs(s["mean_e_r4"] - 4.2) < 0.1,
    )
    rr = out["S1_Bfree_g8_d0.3_w25_n32"]["rows"]
    trend = {k: rr[k]["mean_e_r4"] for k in rr}
    cot2 = (np.arctanh(0.8) - 0.8) / 0.8
    ideal = {
        "8*0.49": 3.92,
        "plus_line_defect_share_weight_0.09": 3.92 + 8 * 0.09 * cot2,
        "plus_line_defect_share_weight_0.0441": 3.92 + 8 * 0.0441 * cot2,
        "uniaxial_8*(1-delta)^4": 8 * 0.7**4,
    }
    flat = max(trend[f"sym r{a}-{b}"] for a, b in ((8, 10), (10, 12), (12, 14), (14, 16))) / min(
        trend[f"sym r{a}-{b}"] for a, b in ((8, 10), (10, 12), (12, 14), (14, 16))
    )
    record(
        "C6.2",
        "the tail weight reads 0.49 (the biaxial G_large^2), not (1-delta)^4 = 0.2401",
        "<e r^4> on four shells, stack and central, against 8 x 0.49 = 3.92 and 8 x 0.2401 = 1.92",
        {"mean_e_r4": trend, "ideal": ideal, "max_over_min_across_shells": flat},
        "between 3.3 and 4.6, never near 1.92",
        all(3.0 < v < 4.7 for v in trend.values()),
        note="fair against 0.2401 (off by 2.2), but <e r^4> is not a plateau: it climbs from 3.4 "
        "(r 8 to 10) to 4.5 "
        "(r 14 to 16), the local G_large^2 runs 0.36 to 0.50, and 10 percent of the shell "
        "curvature is not frame "
        "rotation; '0.49 within about 15 percent' is what the field supports",
    )
    return out


# ================= C7 =================
def c7():
    g, d, e = sp.symbols("g delta epsilon", real=True)
    q = [-g, 1, d, d]
    cp = [sum(qi**p for qi in q) for p in range(1, 5)]
    lam = [-g, 1, d + e, d - e]
    V = sp.expand(sum((sum(li**p for li in lam) - cp[p - 1]) ** 2 for p in range(1, 5)))
    c2_, c4_ = sp.simplify(V.coeff(e, 2)), sp.simplify(V.coeff(e, 4))
    odd = (
        sp.simplify(V.coeff(e, 1)) == 0
        and sp.simplify(V.coeff(e, 3)) == 0
        and sp.simplify(V.coeff(e, 0)) == 0
    )
    ok = c2_ == 0 and sp.expand(c4_ - (4 + 36 * d**2 + 144 * d**4)) == 0 and odd
    record(
        "C7.1",
        "V4 along the split: eps^2 coefficient 0, eps^4 coefficient w (4 + 36 delta^2 + 144 "
        "delta^4)",
        "sympy, exact polynomial in eps",
        {"eps2": str(c2_), "eps4": str(sp.expand(c4_))},
        {"eps2": "0", "eps4": "144*delta**4 + 36*delta**2 + 4"},
        bool(ok),
    )
    # exact Hessians on the vacuum: every (tr N^p - C_p) and P(N) vanish there, so H = 2 J^T J
    basis = []
    for i in range(4):
        B = sp.zeros(4, 4)
        B[i, i] = 1
        basis.append(B)
    names = ["M00", "M11", "M22", "M33"]
    for i in range(4):
        for j in range(i + 1, 4):
            B = sp.zeros(4, 4)
            B[i, j] = B[j, i] = 1 / sp.sqrt(2)
            basis.append(B)
            names.append(f"M{i}{j}")
    eta = sp.diag(-1, 1, 1, 1)
    M0 = sp.diag(g, 1, d, d)  # N = M eta = diag(-g, 1, delta, delta): the code branch s = -1
    N0 = M0 * eta
    resid = [sp.simplify((N0**p).trace() - cp[p - 1]) for p in range(1, 5)]
    J = sp.Matrix(4, 10, lambda p, k: sp.simplify((p + 1) * (N0**p * basis[k] * eta).trace()))
    H4 = 2 * J.T * J
    rank4 = J.rank()
    # V_spec = tr[P(N)^2]
    fac = [N0 - qi * sp.eye(4) for qi in q]
    Pvac = sp.simplify(fac[0] * fac[1] * fac[2] * fac[3])
    dP = []
    for B in basis:
        X = B * eta
        t = sp.zeros(4, 4)
        for k in range(4):
            left = sp.eye(4)
            for i in range(k):
                left = left * fac[i]
            right = sp.eye(4)
            for i in range(k + 1, 4):
                right = right * fac[i]
            t += left * X * right
        dP.append(sp.simplify(t))
    Hs = sp.Matrix(10, 10, lambda a, b: sp.simplify(2 * (dP[a] * dP[b]).trace()))
    rank_s = Hs.rank()
    num = {g: 8, d: sp.Rational(3, 10)}
    ev4 = np.sort(np.linalg.eigvalsh(np.array(H4.subs(num), dtype=float)))
    evs = np.sort(np.linalg.eigvalsh(np.array(Hs.subs(num), dtype=float)))
    null4 = [names[k] for k in range(10) if all(sp.simplify(J[p, k]) == 0 for p in range(4))]
    split_dir = sp.simplify(J * sp.Matrix([0, 0, 1, -1, 0, 0, 0, 0, 0, 0]))
    nulls = [names[k] for k in range(10) if all(sp.simplify(Hs[a, k]) == 0 for a in range(10))]
    nz4, nzs = ev4[-3:], evs[-2:]
    closed_s = [2 * (9 * 0.49) ** 2, 2 * (9 * 8.3**2) ** 2]
    ok4 = (
        rank4 == 3
        and all(r == 0 for r in resid)
        and np.abs(ev4[:7]).max() < 1e-6
        and np.allclose(nz4, [3.0873, 35.0217, 8462877.68], rtol=2e-4)
    )
    record(
        "C7.2",
        "V4 Hessian on the vacuum: 7 zero modes of 10, nonzero about (3.09, 35.0, 8462878) at g "
        "8, delta 0.3, w 1",
        "sympy: H = 2 J^T J, J_pk = p tr(N^(p-1) B_k eta); rank for symbolic (g, delta)",
        {"rank_symbolic": rank4, "nonzero": nz4, "largest_zero": float(np.abs(ev4[:7]).max())},
        {"rank": 3, "nonzero": [3.09, 35.0, 8462878]},
        bool(ok4),
        note="zero modes: the six off-diagonal directions "
        + ", ".join(null4)
        + " (three boosts, two tilts of "
        "the eigenvalue-1 axis, and M23, the off-diagonal partner of the split) plus the diagonal "
        "split "
        f"(M22 - M33) / sqrt 2 (J . split = {list(split_dir)})",
    )
    oks = (
        rank_s == 2
        and Pvac == sp.zeros(4, 4)
        and np.abs(evs[:8]).max() < 1e-6
        and np.allclose(nzs, closed_s, rtol=1e-9)
        and np.allclose(nzs, [38.9, 768825], rtol=2e-3)
    )
    record(
        "C7.3",
        "Vspec Hessian: 8 zero modes, nonzero about 38.9 and 768825",
        "sympy: H_ab = 2 tr(dP_a dP_b); closed form 2 P'(q)^2 at the two simple roots",
        {"rank_symbolic": rank_s, "nonzero": nzs, "closed_form_2Pprime2": closed_s},
        {"rank": 2, "nonzero": [38.9, 768825]},
        bool(oks),
        note="zero modes: "
        + ", ".join(nulls)
        + " are null columns (M22, M33 included: the double root "
        "makes P'(delta) = 0, so both pair diagonals are flat)",
    )
    # production finite-difference Hessian against the exact one (own step, own basis order)
    cfg = B3.base_cfg(s=-1.0, g=G8, n=3, L=3.0, delta=DELTA)
    qn = R0.roots_of(cfg, degenerate=True)
    Mv = np.diag([G8, 1.0, DELTA, DELTA])
    nb = [np.array(B.tolist(), dtype=float) for B in basis]

    def gvac(Mx):
        return R0.v4_energy_grad(np.broadcast_to(Mx, (3, 3, 3, 4, 4)).copy(), cfg, qn, 1.0)[1][
            1, 1, 1
        ]

    Hfd = np.zeros((10, 10))
    for b_, Eb in enumerate(nb):
        dG = (gvac(Mv + 2e-5 * Eb) - gvac(Mv - 2e-5 * Eb)) / 4e-5
        for a_, Ea in enumerate(nb):
            Hfd[a_, b_] = np.sum(dG * Ea)
    dev = float(np.abs(Hfd - np.array(H4.subs(num), dtype=float)).max() / ev4[-1])
    record(
        "C7.4",
        "the production V4 gradient differentiates to the exact Hessian",
        "own FD step 2e-5 against sympy H",
        dev,
        "< 1e-8 of the top eigenvalue",
        dev < 1e-8,
    )


# ================= C8 =================
def cubic_group():
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            R = np.zeros((3, 3))
            for i in range(3):
                R[i, perm[i]] = signs[i]
            yield perm, signs, R


def transform_field(M, perm, signs, R):
    """M'(x) = R4 M(R^T x) R4^T on the cell-centered grid: (R^T x)_perm[i] = signs[i] x_i."""
    T = np.transpose(M, tuple(perm) + (3, 4))
    for i in range(3):
        if signs[i] < 0:
            T = np.flip(T, axis=i)
    R4 = np.eye(4)
    R4[1:, 1:] = R
    return np.ascontiguousarray(np.einsum("ab,...bc,dc->...ad", R4, T, R4))


def c8():
    n_, L = 12, 9.0
    (X, Y, Z), h = grid(n_, L)
    rng = np.random.default_rng(4242)
    par = {
        (a, b): (
            rng.normal(size=(2, 3)) * 0.5,
            rng.uniform(0, 6.28, size=2),
            rng.normal(size=2) * 0.4,
        )
        for a in range(4)
        for b in range(a, 4)
    }

    def field_at(Xq, Yq, Zq):
        Mq = np.zeros(Xq.shape + (4, 4))
        for (a, b), (k, p_, amp) in par.items():
            f = sum(
                amp[m] * np.sin(k[m, 0] * Xq + k[m, 1] * Yq + k[m, 2] * Zq + p_[m])
                for m in range(2)
            )
            Mq[..., a, b] = Mq[..., b, a] = f
        return Mq + np.diag([G8, 1.0, DELTA, 0.0])

    M = field_at(X, Y, Z)
    cfg = B3.base_cfg(s=-1.0, g=G8, n=n_, L=L, delta=DELTA)
    E0 = {st: float(B3.e_parts(M, cfg, st)[0]) for st in ("sym", "cen")}
    own0 = {st: float(own_density(M, h, st).sum() * h**3) for st in ("sym", "cen")}
    inv = {"sym": [], "cen": []}
    dev = {"sym": [], "cen": []}
    worst_tf = 0.0
    for perm, signs, R in cubic_group():
        Mt = transform_field(M, perm, signs, R)
        # the array transform against the definition M'(x) = R4 M(R^T x) R4^T,
        # evaluated analytically
        Q = [sum(R[i, j] * (X, Y, Z)[i] for i in range(3)) for j in range(3)]
        R4 = np.eye(4)
        R4[1:, 1:] = R
        worst_tf = max(
            worst_tf,
            float(np.abs(Mt - np.einsum("ab,...bc,dc->...ad", R4, field_at(*Q), R4)).max()),
        )
        for st in ("sym", "cen"):
            rel = abs(float(B3.e_parts(Mt, cfg, st)[0]) - E0[st]) / abs(E0[st])
            dev[st].append(rel)
            if rel < 1e-12:
                inv[st].append((perm, signs))
    all_equal = all(len(set(s)) == 1 for _, s in inv["sym"])
    ok = len(inv["sym"]) == 12 and len(inv["cen"]) == 48 and all_equal and worst_tf < 1e-12
    broken = [v for v in dev["sym"] if v >= 1e-12]
    record(
        "C8.1",
        "the fwd / bwd averaged functional is invariant under exactly 12 of the 48 cubic images, "
        "the central one under 48",
        "generic asymmetric smooth 4x4 field (own seed), production B3.e_parts on every image",
        {
            "sym_invariant": len(inv["sym"]),
            "cen_invariant": len(inv["cen"]),
            "the_12_have_equal_signs": bool(all_equal),
            "array_transform_vs_definition_max_abs": worst_tf,
            "sym_broken_rel_min": min(broken),
            "sym_broken_rel_max": max(broken),
            "own_vs_production_E": {st: abs(own0[st] - E0[st]) / abs(E0[st]) for st in E0},
        },
        {"sym_invariant": 12, "cen_invariant": 48},
        bool(ok),
        note="per plaquette (i, j) the fwd-fwd and bwd-bwd branches use the corner triangles at "
        "(0,0) and (1,1); a "
        "reflection of one of the two axes maps them to the corners (1,0) and (0,1), which the "
        "average does not "
        "contain; equal signs on all three axes (identity or full inversion) swap fwd and bwd as "
        "a whole",
    )


# ================= C9 =================
def omega_own(M, h, delta=DELTA, renv=10.0, gen=(0, 1)):
    n_ = M.shape[0]
    x = (np.arange(n_) - (n_ - 1) / 2.0) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z)
    env = np.exp(-((r / renv) ** 4)) if np.isfinite(renv) else np.ones_like(r)
    W = np.zeros((4, 4))
    W[gen[0], gen[1]] = W[gen[1], gen[0]] = 1.0
    a = env[..., None, None] * (W @ M + M @ W.T)
    nrm = float(np.sqrt(np.sum(a * a)))
    a0 = a / nrm
    xm = ETA @ M
    I4 = np.eye(4)
    q = (xm @ (xm - I4) @ (xm - delta * I4)) / (G8 * (G8 - 1) * (G8 - delta))
    Gm = ETA - 2.0 * q @ ETA
    ip = lambda F: np.einsum("...ab,...ac,...bd,...cd->...", F, Gm, Gm, F)  # noqa: E731
    i1, k1 = 0.0, 0.0
    for kind in ("fwd", "bwd"):
        A = [own_diff(M, ax, h, kind) for ax in range(3)]
        for i in range(3):
            k1 = k1 + 2.0 * ip(a0 @ ETA @ A[i] - A[i] @ ETA @ a0)
            for j in range(i + 1, 3):
                i1 = i1 + 2.0 * ip(A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i])
    C1, C2 = h**3 * float(np.sum(i1 * k1)), h**3 * float(np.sum(k1 * k1))
    return {
        "C1": C1,
        "C2": C2,
        "omega_E": float(np.sqrt(C1 / C2)) if C1 > 0 else None,
        "tangent_norm": nrm,
        "I1G_static": float(h**3 * np.sum(i1)),
    }


def c9():
    pol = np.load(os.path.join(EXT_008, "M_G_polished_N32.npz"))["M"]
    seed3 = np.load(
        os.path.join(EXT_008, "seeds", "m5_21_2b_end_A_T2_sym_e0_n32_d0.3_pinned.npz")
    )["M"].astype(float)
    with open(os.path.join(EXT_008, "chain_N32.json")) as f:
        prof = json.load(f)["profile"]
    N = 32
    mask = np.zeros((N, N, N), dtype=bool)
    wc = int(np.ceil(1.6 / 1.5))
    for ax in range(3):
        sl = [slice(None)] * 3
        sl[ax] = slice(0, wc)
        mask[tuple(sl)] = True
        sl[ax] = slice(N - wc, N)
        mask[tuple(sl)] = True
    seed4 = np.zeros((N, N, N, 4, 4))
    seed4[..., 1:, 1:] = seed3
    seed4[..., 0, 0] = -G8
    Mc = np.where(mask[..., None, None], seed4, 0.5 * (pol + np.swapaxes(pol, -1, -2)))
    ctl = omega_own(Mc, 1.5)
    ok = (
        abs(ctl["C1"] / prof["C1"] - 1) < 1e-9
        and abs(ctl["C2"] / prof["C2"] - 1) < 1e-9
        and abs(ctl["omega_E"] - prof["omega_pred_E"]) < 1e-9
    )
    record(
        "C9.1",
        "control: C1 9.0754e-06, C2 8.5648e-05, omega_E 0.3255186 on the committed N 32 field",
        "own retyping of the route, numbers read from chain_N32.json",
        {k: ctl[k] for k in ("C1", "C2", "omega_E")},
        {"C1": prof["C1"], "C2": prof["C2"], "omega_E": prof["omega_pred_E"]},
        ok,
    )
    # the branch translation on a random block-diagonal field, production energy on both branches
    n_, L = 10, 12.0
    (X, Y, Z), h = grid(n_, L)
    rng = np.random.default_rng(31337)
    Mb = np.zeros((n_, n_, n_, 4, 4))
    for a, b in [(0, 0)] + [(a, b) for a in range(1, 4) for b in range(a, 4)]:
        k = rng.normal(size=3) * 0.6
        Mb[..., a, b] = Mb[..., b, a] = 0.3 * np.sin(
            k[0] * X + k[1] * Y + k[2] * Z + rng.uniform(0, 6.28)
        )
    Mb += np.diag([G8, 1.0, DELTA, 0.0])  # code branch s = -1: M_00 = +g
    cm, cp_ = B3.base_cfg(s=-1.0, g=G8, n=n_, L=L, delta=DELTA), B3.base_cfg(
        s=+1.0, g=G8, n=n_, L=L, delta=DELTA
    )

    def both(Mx):
        Mf = Mx.copy()
        Mf[..., 0, 0] *= -1.0
        (eu_m, ev_m), (eu_p, ev_p) = B3.e_parts(Mx, cm), B3.e_parts(Mf, cp_)
        return abs(eu_m - eu_p) / abs(eu_m), abs(ev_m - ev_p) / abs(ev_m)

    curv_var, v4_var = both(Mb)  # M_00 varies, spatial block off target
    Mconst = Mb.copy()
    Mconst[..., 0, 0] = G8
    curv_const, v4_const = both(Mconst)  # M_00 = g exactly
    Mb2 = Mb.copy()
    bump = 0.2 * np.sin(0.5 * X + 0.3 * Y)
    Mb2[..., 0, 1] += bump
    Mb2[..., 1, 0] += bump
    neg = both(Mb2)[0]
    record(
        "C9.2",
        "block-diagonal static field: the curvature is unchanged under M_00 -> -M_00 (branch s = "
        "-1 to s = +1)",
        "random block-diagonal field with varying M_00, production B3.e_parts on both branches; "
        "negative control M_01 != 0",
        {
            "curvature_rel_diff_M00_varying": curv_var,
            "curvature_rel_diff_M00_constant": curv_const,
            "negative_control_with_M01": neg,
        },
        {"curvature": "< 1e-12", "negative_control": "> 1e-6"},
        curv_var < 1e-12 and curv_const < 1e-12 and neg > 1e-6,
    )
    r21_v4 = {}
    for tag in ("S1_Bseed_g8_d0.3_w25_n32", "S1_Bfree_g8_d0.3_w25_n32"):
        Mr = np.load(os.path.join(R21_DIR, tag + ".npz"))["M"]
        c_m, c_p = B3.base_cfg(s=-1.0, g=G8, n=32, L=48.0, delta=DELTA), B3.base_cfg(
            s=+1.0, g=G8, n=32, L=48.0, delta=DELTA
        )
        Mr_f = Mr.copy()
        Mr_f[..., 0, 0] *= -1.0
        v_m, v_p = B3.e_parts(Mr, c_m)[1], B3.e_parts(Mr_f, c_p)[1]
        r21_v4[tag] = {
            "V4_code_branch": float(v_m),
            "V4_translated": float(v_p),
            "rel_diff": float(abs(v_m - v_p) / abs(v_m)),
            "M00_minus_g_range": [
                float(Mr[..., 0, 0].min() - G8),
                float(Mr[..., 0, 0].max() - G8),
            ],
        }
    record(
        "C9.2b",
        "V4 is unchanged too, also where M_00 varies ('M_00 may vary', 'V4's odd powers flip "
        "together with their targets')",
        "the same random field through production V4 on both branches; then the two R21 fields",
        {
            "V4_rel_diff_M00_varying": v4_var,
            "V4_rel_diff_M00_constant": v4_const,
            "R21_fields": r21_v4,
        },
        "< 1e-12",
        v4_var < 1e-12,
        note="tr N^p - C_p = (u^p - sg^p) + (S_p - c_p) with u = N_00 and S_p the spatial trace; "
        "under the flip the first "
        "bracket changes sign for odd p, so the cross term 2 (u^p - sg^p)(S_p - c_p), p = 1 and "
        "3, changes sign: V4 "
        "is invariant only where M_00 = g exactly or the spatial block sits on target. On the R21 "
        "fields M_00 - g is "
        "only of order 1e-4, but u^3 and u^4 amplify it (4 g^3 = 2048) and V4 moves by 16 to 18 "
        "percent under the "
        "translation: the translated field is NOT a stationary point of the other branch's static "
        "energy. C1, C2 "
        "and omega_E do not contain V4, so the omega_E numbers themselves stand",
    )
    rows, sens = {}, {}
    for tag, claimed in (("S1_Bseed_g8_d0.3_w25_n32", 4.00), ("S1_Bfree_g8_d0.3_w25_n32", 4.59)):
        M = np.load(os.path.join(R21_DIR, tag + ".npz"))["M"].copy()
        off = float(np.abs(M[..., 0, 1:]).max())
        cfg = B3.base_cfg(s=-1.0, g=G8, n=32, L=48.0, delta=DELTA)
        e_eta = float(B3.e_parts(M, cfg)[0])
        M[..., 0, 0] *= -1.0
        base = omega_own(M, 1.5)
        rows[tag] = dict(
            base,
            max_abs_M0i=off,
            claimed=claimed,
            E_curv_eta=e_eta,
            I1G_over_E_curv_eta=base["I1G_static"] / e_eta,
        )
        sens[tag] = {
            "renv": {
                str(re): omega_own(M, 1.5, renv=re)["omega_E"]
                for re in (6.0, 8.0, 12.0, 16.0, 24.0, np.inf)
            },
            "boost_y": omega_own(M, 1.5, gen=(0, 2))["omega_E"],
            "boost_z": omega_own(M, 1.5, gen=(0, 3))["omega_E"],
        }
        sens[tag]["renv"]["10.0"] = base["omega_E"]
    okr = all(
        abs(rows[t]["omega_E"] - rows[t]["claimed"]) < 0.01 and rows[t]["max_abs_M0i"] < 1e-12
        for t in rows
    )
    record(
        "C9.3",
        "omega_E of the two R21 W1 x 25 rows: 4.00 (seed pin) and 4.59 (free)",
        "own route on the translated fields",
        {
            t: {
                k: rows[t][k]
                for k in ("omega_E", "C1", "C2", "max_abs_M0i", "I1G_over_E_curv_eta")
            }
            for t in rows
        },
        {"S1_Bseed": 4.00, "S1_Bfree": 4.59},
        okr,
        note="the read is formally the report's: M_0i = 0 exactly, and the G-metric I_1 equals "
        "the eta-metric curvature "
        "to 1e-4 on these fields; what the number means is limited by C9.2b, C9.4 and C9.5",
    )
    ctl_s = {
        str(re): omega_own(Mc, 1.5, renv=re)["omega_E"]
        for re in (6.0, 8.0, 12.0, 16.0, 24.0, np.inf)
    }
    ctl_s["10.0"] = ctl["omega_E"]
    t = "S1_Bseed_g8_d0.3_w25_n32"
    spread = sens[t]["renv"]["16.0"] / sens[t]["renv"]["8.0"]
    ratio = {k: sens[t]["renv"][k] / ctl_s[k] for k in ctl_s}
    record(
        "C9.4",
        "omega_E is a meaningful number on these fields (weak dependence on the envelope radius "
        "10)",
        "envelope radius 6 to infinity, boost axis x / y / z, the control swept the same way",
        {
            "fields": sens,
            "control": ctl_s,
            "omega_E(16)/omega_E(8)": spread,
            "ratio_to_control_by_renv": ratio,
        },
        "omega_E(16) / omega_E(8) within 1.25",
        spread < 1.25,
        note="omega_E runs 2.3 to 22 (seed pin) and 0.22 to 3.5 (the control itself) over renv 6 "
        "to infinity, about "
        "linear in renv; only the RATIO to the control is stable (6 to 13); boost-z differs from "
        "boost-x by "
        "15 percent (seed pin) because the seed's line defect is on z",
    )
    # the unit Frobenius norm is a cell sum: the same continuum field at half the spacing
    om = {}
    for n2 in (32, 64):
        cfg = R20.cfg_of(n2, 48.0, G8)
        Ms = R20.seed_axes(cfg, R20.OBJECTS["S1"]).copy()
        Ms[..., 0, 0] *= -1.0
        om[n2] = omega_own(Ms, cfg["h"])
    rat = om[64]["omega_E"] / om[32]["omega_E"]
    record(
        "C9.5",
        "omega_E does not depend on the lattice spacing of the field it is read on",
        "the analytic S1 seed sampled at n 32 and n 64 in the same L 48 box",
        {
            "omega_E_n32": om[32]["omega_E"],
            "omega_E_n64": om[64]["omega_E"],
            "ratio": rat,
            "2^1.5": 2**1.5,
            "tangent_norm_ratio": om[64]["tangent_norm"] / om[32]["tangent_norm"],
        },
        "ratio within 1 +/- 0.2",
        abs(rat - 1) < 0.2,
        note="a0 is normalized by the cell sum of a^2 (no h^3), so omega_E scales as h^(-3/2) at "
        "fixed continuum "
        "field (here 3.41 = 2.83 from the norm times 1.21 from the lattice); rows read at h other "
        "than 1.5 (the "
        "planned R22-1 n 48, L 48 rows, h 1.0: factor 1.84) are not comparable with the report's "
        "numbers unless "
        "rescaled by (h / 1.5)^(3/2)",
    )
    return {"rows": rows, "sensitivity": sens, "control_sensitivity": ctl_s}


# ================= C10 =================
def c10():
    n_, L, dl = 11, 13.2, 0.37
    (X, Y, Z), h = grid(n_, L)
    rng = np.random.default_rng(77)
    m = np.stack(
        [
            np.sin(0.4 * X + 0.2 * Z) + 0.3,
            np.cos(0.3 * Y - 0.25 * X),
            0.8 + 0.5 * np.sin(0.35 * Z + 0.1 * Y),
        ],
        -1,
    )
    m = m + 0.15 * rng.normal(size=m.shape)
    m *= rng.uniform(0.7, 1.4, size=m.shape[:3])[
        ..., None
    ]  # |m| not 1: the chain through m / |m| is exercised
    E, g = PAIR22.energy_grad_n(m, h, dl)
    nm = np.linalg.norm(m, axis=-1, keepdims=True)
    n = m / nm
    cfg = B3.base_cfg(s=-1.0, g=G8, n=n_, L=L, delta=dl)
    M = B3.embed34(dl * np.eye(3) + (1 - dl) * outer(n, n), cfg)
    E_prod = float(B3.e_parts(M, cfg)[0])
    E_own = float(own_density(M, h, "sym").sum() * h**3)
    Gc = B3.grad(M, cfg) - R0.v4_energy_grad(M, cfg, R0.roots_of(cfg), B3.W1)[1]
    gn = 2.0 * (1 - dl) * np.einsum("...ab,...b->...a", Gc[..., 1:, 1:], n)
    g_prod = (gn - n * np.sum(gn * n, -1, keepdims=True)) / nm
    rel_g = float(np.abs(g - g_prod).max() / np.abs(g_prod).max())
    worst_fd = 0.0
    for _ in range(6):
        v = rng.normal(size=m.shape)
        eps = 1e-5
        fd = (
            PAIR22.energy_grad_n(m + eps * v, h, dl, need_grad=False)[0]
            - PAIR22.energy_grad_n(m - eps * v, h, dl, need_grad=False)[0]
        ) / (2 * eps)
        worst_fd = max(worst_fd, abs(fd - float(np.sum(g * v))) / abs(fd))
    cells = [(0, 0, 0), (n_ - 1, 3, 5), (4, n_ - 1, 0), (5, 5, 5), (1, 0, n_ - 2)]
    worst_cell = 0.0
    for c in cells:
        for a in range(3):
            dm = np.zeros_like(m)
            dm[c + (a,)] = 1e-5
            fd = (
                PAIR22.energy_grad_n(m + dm, h, dl, need_grad=False)[0]
                - PAIR22.energy_grad_n(m - dm, h, dl, need_grad=False)[0]
            ) / 2e-5
            worst_cell = max(worst_cell, abs(fd - g[c + (a,)]) / np.abs(g).max())
    val = {
        "E_rel_vs_production": abs(E - E_prod) / abs(E_prod),
        "E_rel_vs_own": abs(E - E_own) / abs(E_own),
        "grad_rel_vs_production_chain": rel_g,
        "fd_directional_worst_rel": worst_fd,
        "fd_single_cell_worst": worst_cell,
    }
    record(
        "C10.1",
        "energy_grad_n equals the production energy and gradient (chain rule to n) and passes a "
        "finite-difference test",
        "own field (|m| not 1, noise), n 11, L 13.2, delta 0.37, seed 77; faces and corners "
        "included in the cell test",
        val,
        {"E": "< 1e-12", "grad": "< 1e-10", "fd": "< 1e-6"},
        val["E_rel_vs_production"] < 1e-12
        and val["E_rel_vs_own"] < 1e-12
        and rel_g < 1e-10
        and worst_fd < 1e-6
        and worst_cell < 1e-6,
    )


# ================= C11 =================
def green_series(a, x1, x2, mmax=241):
    """Dirichlet Green function of the cube [-a, a]^3 (G -> 1 / |x - x'|), sine modes in y, z,
    closed form in x.
    Both points on the axis y = z = 0 (only odd modes survive)."""
    Lb = 2.0 * a
    xa, xb = sorted((x1 + a, x2 + a))
    m = np.arange(1, mmax + 1, 2)
    gam = np.pi * np.sqrt(m[:, None] ** 2 + m[None, :] ** 2) / Lb
    A, B, C = gam * xa, gam * (Lb - xb), gam * Lb
    term = (
        np.exp(-gam * (xb - xa))
        * (1 - np.exp(-2 * A))
        * (1 - np.exp(-2 * B))
        / (2 * gam * (1 - np.exp(-2 * C)))
    )
    return float(4 * np.pi * (2.0 / Lb) ** 2 * term.sum())


def green_dst(a, n_, x1, x2):
    """own DST-II solve of the regular part on a cell-centered grid, ghost-cell walls at +/- a."""
    from scipy.fft import dstn, idstn

    h = 2 * a / n_
    c = (np.arange(n_) - (n_ - 1) / 2.0) * h
    rhs = np.zeros((n_, n_, n_))
    U, V = np.meshgrid(c, c, indexing="ij")
    for ax in range(3):
        for wall, idx in ((-a, 0), (a, n_ - 1)):
            p = [None] * 3
            o = [k for k in range(3) if k != ax]
            p[o[0]], p[o[1]], p[ax] = U, V, np.full_like(U, wall)
            gw = -1.0 / np.sqrt((p[0] - x1) ** 2 + p[1] ** 2 + p[2] ** 2)
            sl = [slice(None)] * 3
            sl[ax] = idx
            rhs[tuple(sl)] -= 2.0 * gw / h**2
    lam = -(2 - 2 * np.cos(np.pi * np.arange(1, n_ + 1) / n_)) / h**2
    psi = idstn(
        dstn(rhs, type=2) / (lam[:, None, None] + lam[None, :, None] + lam[None, None, :]), type=2
    )
    # trilinear read at (x2, 0, 0)
    val = psi
    for xq in (x2, 0.0, 0.0):
        t = (xq - c[0]) / h
        i0 = int(np.floor(t))
        f = t - i0
        val = val[i0] * (1 - f) + val[i0 + 1] * f
    return 1.0 / abs(x1 - x2) + float(val)


def c11():
    x, y, z, ep = sp.symbols("x y z epsilon", real=True)
    a, b, u, v = [sp.Function(nm)(x, y, z) for nm in ("a", "b", "u", "v")]
    X = (x, y, z)

    def n_of(aa, bb):
        return sp.Matrix([sp.sin(aa) * sp.cos(bb), sp.sin(aa) * sp.sin(bb), sp.cos(aa)])

    def etop(nn):
        d = [nn.diff(c) for c in X]
        return sp.Matrix([nn.dot(d[i].cross(d[j])) for i, j in ((1, 2), (2, 0), (0, 1))])

    n0 = n_of(a, b)
    dE = etop(n_of(a + ep * u, b + ep * v)).diff(ep).subs(ep, 0)
    dn = n_of(a + ep * u, b + ep * v).diff(ep).subs(ep, 0)
    alpha = sp.Matrix([n0.dot(dn.cross(n0.diff(c))) for c in X])
    curl = sp.Matrix(
        [
            alpha[2].diff(y) - alpha[1].diff(z),
            alpha[0].diff(z) - alpha[2].diff(x),
            alpha[1].diff(x) - alpha[0].diff(y),
        ]
    )
    E = etop(n0)
    r1 = sp.simplify(dE - curl)
    r2 = sp.simplify(alpha.dot(E))
    r3 = sp.simplify(sum((E[k] * n0.diff(X[k]) for k in range(3)), sp.zeros(3, 1)))
    ok = r1 == sp.zeros(3, 1) and r2 == 0 and r3 == sp.zeros(3, 1)
    record(
        "C11.1",
        "natural wall condition of int |E_top|^2 is E_top x nu = 0 (E_top normal to the wall)",
        "sympy with arbitrary a(x), b(x) and variation (u, v): delta E = curl alpha, alpha . E = "
        "0, E . grad n = 0",
        [str(list(r1)), str(r2), str(list(r3))],
        ["[0, 0, 0]", "0", "[0, 0, 0]"],
        bool(ok),
        note="delta E = 2 int alpha . curl E + 2 oint alpha . (E x nu); alpha spans the plane "
        "normal to E, so the wall "
        "term gives E x nu = 0 (equivalently d_nu n = 0); curl variations carry no flux, so the "
        "fixed flux adds "
        "no multiplier and no Neumann condition on the potential. The BULK equation is only (curl "
        "E) x E = 0 "
        "(force-free), weaker than curl E = 0: 8 pi k G_D is the lower bound over all "
        "divergence-free fields and "
        "is attained only if the grounded-box Coulomb field is realized by a director texture",
    )
    g24, g32 = green_series(24.0, 4.0, -4.0), green_series(32.0, 4.0, -4.0)
    conv = abs(green_series(24.0, 4.0, -4.0, mmax=161) - g24)
    record(
        "C11.2",
        "G_D of the cube, source (4,0,0), field point (-4,0,0): 0.08935 (half-width 24), 0.09802 "
        "(half-width 32)",
        "mixed eigenfunction series (odd sine modes to 241, closed sinh form along x)",
        {"half24": g24, "half32": g32, "truncation_change_161_to_241": conv},
        {"half24": 0.08935, "half32": 0.09802},
        abs(g24 - 0.08935) < 1e-4 and abs(g32 - 0.09802) < 1e-4 and conv < 1e-10,
    )
    prod = {
        "n48_L48": PAIR22.green_box(8.0, 48, 48.0)[0],
        "n32_L48": PAIR22.green_box(8.0, 32, 48.0)[0],
        "n64_L64": PAIR22.green_box(8.0, 64, 64.0)[0],
    }
    own_fine = {
        "half24_n96": green_dst(24.0, 96, 4.0, -4.0),
        "half32_n128": green_dst(32.0, 128, 4.0, -4.0),
    }
    dev = max(
        abs(prod["n48_L48"] - g24),
        abs(prod["n32_L48"] - g24),
        abs(prod["n64_L64"] - g32),
        abs(own_fine["half24_n96"] - g24),
        abs(own_fine["half32_n128"] - g32),
    )
    record(
        "C11.3",
        "green_box (the audited DST solve) and an own finer DST solve agree with the series",
        "three routes",
        {
            "green_box": prod,
            "own_dst_fine": own_fine,
            "series": {"half24": g24, "half32": g32},
            "max_abs_dev": dev,
        },
        "< 1e-4",
        dev < 1e-4,
    )
    # which wall does the lattice functional see? count the layers that carry each field component
    n_, L = 16, 16.0
    (X, Y, Z), h = grid(n_, L)
    al_, be_ = 0.02, 0.02

    def uniform_rho(u, w):  # cos a = 0.2 - al u, b = be w: rho_uw = al be exactly
        ca = 0.2 - al_ * u
        sa = np.sqrt(1 - ca * ca)
        return np.stack([sa * np.cos(be_ * w), sa * np.sin(be_ * w), ca], -1)

    layers = {}
    for name, nn in (
        ("E_z_from_xy_plaquettes", uniform_rho(X, Y)),
        ("E_x_from_yz_plaquettes", uniform_rho(Y, Z)),
    ):
        e = own_density(outer(nn, nn), h, "sym")
        layers[name] = float(e.sum() / e[n_ // 2, n_ // 2, n_ // 2])
    shift = {
        "half_L/2": g24,
        "half_L/2-h/2_at_h1": green_series(23.5, 4.0, -4.0),
        "half_L/2-h/2_at_h1.5": green_series(23.25, 4.0, -4.0),
        "free_space": 0.125,
    }
    rel = (shift["half_L/2-h/2_at_h1.5"] - g24) / g24
    ok = (
        abs(layers["E_z_from_xy_plaquettes"] - n_ * (n_ - 1) ** 2) < 2
        and abs(layers["E_x_from_yz_plaquettes"] - n_ * (n_ - 1) ** 2) < 2
    )
    record(
        "C11.4",
        "the wall of the cell-centered free box sits at L/2 (the choice made in psi_dirichlet)",
        "effective cell counts of a uniform-rho texture (n 16): a component is integrated over "
        "the full width L along "
        "its own axis and over L - h across it; then G_D at L/2 against L/2 - h/2",
        {
            "effective_cells": layers,
            "n*(n-1)^2": n_ * (n_ - 1) ** 2,
            "n^3": n_**3,
            "G_D": shift,
            "rel_change_h1.5": rel,
        },
        "each component: n (n-1)^2 cells",
        bool(ok),
        note="the wall-normal component of E_top is weighted out to L/2, the wall-tangential ones "
        "only to the outer "
        "cell centers; since the tangential field vanishes linearly at a conductor, this is an "
        "O(h^2) quadrature "
        "of the wall-at-L/2 problem (Hadamard: a half cell of normal-field energy equals a wall "
        "displacement "
        "of h/2 to first order). Putting the wall at L/2 - h/2 instead would lower G_D by "
        f"{-rel:.1%} at h 1.5: that is the size of the wall-position question, not a defect",
    )


# ================= C12 =================
def c12():
    r, c, E0, rc = sp.symbols("r c E_0 r_c", positive=True)
    p = -E0 * rc**3 / 2
    # flux-excluding sphere: phi = -E0 r c + p c / r^2 ; radial field at r_c must vanish
    phi = -E0 * r * c + p * c / r**2
    Er = -phi.diff(r)
    bc = sp.simplify(Er.subs(r, rc))
    # dipole field energy outside: |E_dip|^2 = p^2 (3 c^2 + 1) / r^6
    out_E = sp.integrate(
        sp.integrate(p**2 * (3 * c**2 + 1) / r**6 * 2 * sp.pi * r**2, (c, -1, 1)), (r, rc, sp.oo)
    )
    cross = sp.integrate(
        E0 * p * (3 * c**2 - 1) / r**3, (c, -1, 1)
    )  # E_0 . E_dip, angular integral
    in_E = sp.Rational(4, 3) * sp.pi * rc**3 * E0**2
    net = sp.simplify(out_E - in_E)
    ok = (
        bc == 0
        and sp.simplify(out_E - 2 * sp.pi / 3 * rc**3 * E0**2) == 0
        and cross == 0
        and sp.simplify(net + 2 * sp.pi / 3 * rc**3 * E0**2) == 0
    )
    record(
        "C12.1",
        "pinned ball of fixed uniform flux in a uniform field: outside +(2 pi / 3), inside -(4 pi "
        "/ 3), net "
        "-(2 pi / 3) k r_c^3 E_0^2 per ball",
        "sympy: induced dipole p = -E_0 r_c^3 / 2 from zero normal field, "
        "energy integrals, cross term",
        {
            "normal_field_at_r_c": str(bc),
            "outside": str(sp.simplify(out_E)),
            "cross_angular": str(cross),
            "net": str(net),
        },
        {"net": "-2*pi*E_0**2*r_c**3/3"},
        bool(ok),
    )
    # general multipole l, numerically:
    # net_l = -(l / (l + 1)) (4 pi / (2 l + 1)) a_l^2 r_c^(2 l + 1)
    from numpy.polynomial import legendre as npl
    from scipy.integrate import quad

    cg, wg = npl.leggauss(64)
    worst = 0.0
    for l_ in (1, 2, 3, 4):
        rcn = 2.5
        bl = l_ * rcn ** (2 * l_ + 1) / (l_ + 1)
        P = npl.Legendre.basis(l_)
        dP = P.deriv()

        def ang(rr, a_in, b_out):
            pr = (a_in * l_ * rr ** (l_ - 1) - b_out * (l_ + 1) * rr ** (-l_ - 2)) * P(cg)
            pt = (a_in * rr ** (l_ - 1) + b_out * rr ** (-l_ - 2)) * (-np.sqrt(1 - cg**2)) * dP(cg)
            return 2 * np.pi * np.sum(wg * (pr**2 + pt**2)) * rr**2

        e_out_tot = quad(lambda rr: ang(rr, 1.0, bl) - ang(rr, 1.0, 0.0), rcn, 20.0, limit=200)[0]
        e_tail = (
            (l_ + 1) * bl**2 * 20.0 ** (-2 * l_ - 1) * 4 * np.pi / (2 * l_ + 1)
        )  # induced part beyond r 20
        e_in = quad(lambda rr: ang(rr, 1.0, 0.0), 0.0, rcn)[0]
        closed = -(l_ / (l_ + 1)) * 4 * np.pi / (2 * l_ + 1) * rcn ** (2 * l_ + 1)
        worst = max(worst, abs((e_out_tot + e_tail - e_in) - closed) / abs(closed))
    k = 8 * (1 - DELTA) ** 4
    rows = {}
    for d_ in (6.0, 8.0, 12.0):
        dip = -(4 * np.pi / 3) * k * 2.5**3 / d_**4
        full = (
            -2
            * k
            * sum(
                (l_ / (l_ + 1))
                * 4
                * np.pi
                / (2 * l_ + 1)
                * 2.5 ** (2 * l_ + 1)
                / d_ ** (2 * l_ + 2)
                for l_ in range(1, 40)
            )
        )
        rows[str(d_)] = {
            "dipole_only": dip,
            "all_multipoles": full,
            "ratio": full / dip,
            "share_of_U_free": abs(full) / (8 * np.pi * k / d_),
        }
    okd = (
        abs(rows["6.0"]["dipole_only"] + 0.097) < 1e-3
        and abs(rows["8.0"]["dipole_only"] + 0.031) < 1e-3
        and abs(rows["12.0"]["dipole_only"] + 0.006) < 1e-3
    )
    record(
        "C12.2",
        "per pair -(4 pi / 3) k r_c^3 / d^4: -0.097 (d 6), -0.031 (d 8), -0.006 (d 12) at r_c 2.5",
        "arithmetic, then the full multipole sum of a point charge at distance d (net_l checked "
        "by quadrature, l 1 to 4)",
        {"rows": rows, "multipole_formula_quadrature_worst_rel": worst},
        {"6": -0.097, "8": -0.031, "12": -0.006},
        okd and worst < 1e-6,
        note="the cross term E_0 . E_dip is conditionally convergent and vanishes shell by shell "
        "(angular integral "
        "first), the order the docstring uses. The uniform-field (dipole) estimate is the l = 1 "
        "term; the field of the other core is not uniform over "
        "the ball, and the l >= 2 terms raise the artifact by 16 percent at d 6 (-0.113), 9 "
        "percent at d 8, "
        "4 percent at d 12; same sign, still charge-even",
    )


# ================= X: outside the listed claims =================
def cx():
    cfg = B3.base_cfg(s=-1.0, g=G8, n=2, L=3.0, delta=DELTA)
    m00 = float(B3.embed34(np.zeros((2, 2, 2, 3, 3)), cfg)[0, 0, 0, 0, 0])
    record(
        "X.1",
        "R22-0 (a) docstring: 'M_00 = -g constant (code branch s = -1, the R20 / R21 stack)'",
        "B3.embed34 on s = -1",
        m00,
        -G8,
        m00 == -G8,
        note="the code branch s = -1 has M_00 = +g (N_00 = -g), as the R22-3 docstring says; "
        "immaterial for the spatial "
        "quartic, but the two docstrings contradict each other",
    )


# ================= main =================
def main():
    only = None
    for i, a in enumerate(sys.argv):
        if a == "--only" and i + 1 < len(sys.argv):
            only = set(sys.argv[i + 1].split(","))
    extra = {}
    for name, fn in (
        ("C1", c1),
        ("C2", c2),
        ("C3", c3),
        ("C4", c4),
        ("C5", c5),
        ("C6", c6),
        ("C7", c7),
        ("C8", c8),
        ("C9", c9),
        ("C10", c10),
        ("C11", c11),
        ("C12", c12),
        ("X", cx),
    ):
        if only and name not in only:
            continue
        t = time.time()
        r = fn()
        if isinstance(r, (dict, list)):
            extra[name] = _js(r)
        print(f"--- {name} done in {time.time() - t:.1f}s", flush=True)
    n_pass = sum(c["verdict"] == "PASS" for c in CHECKS)
    n_fail = len(CHECKS) - n_pass
    res = {
        "task": "M5.32 R22-0 adversarial audit",
        "totals": {
            "checks": len(CHECKS),
            "PASS": n_pass,
            "FAIL": n_fail,
            "PASS_with_wording_note": sum(
                1 for c in CHECKS if c["verdict"] == "PASS" and c["note"]
            ),
        },
        "checks": CHECKS,
        "tables": extra,
        "runtime_s": time.time() - T0,
    }
    if only:
        print(f"partial run ({sorted(only)}): JSON not written")
    else:
        with open(OUT_JSON, "w") as f:
            json.dump(res, f, indent=1)
    print(
        f"TOTAL {len(CHECKS)} checks: PASS {n_pass}, FAIL {n_fail}  ({time.time() - T0:.0f}s) -> "
        f"{os.path.basename(OUT_JSON)}"
    )


if __name__ == "__main__":
    main()
