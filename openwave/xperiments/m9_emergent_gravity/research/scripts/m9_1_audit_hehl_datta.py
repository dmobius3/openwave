#!/usr/bin/env python3
"""M9.1 adversarial audit: Hodge dual + Palatini 4-form, no solver import.

Method (intentionally not the solver's):
  1. Contract Levi-Civita identities in sympy: (*J)_{lmn} (*J)^{lmn} = 6 J·J.
  2. Convert the Palatini 4-form
        S = (1/(4 kappa)) int eps_{abcd} e^a /\\ e^b /\\ R^{cd}
     to the scalar density e R / (2 kappa) by expanding the 4-form on
     Minkowski (sympy, constant e=delta).
  3. For totally antisymmetric K_{abc} = eps_{abcr} v^r, reduce Q to a
     multiple of v·v by the same eps identity.
  4. Pairing L_D = (1/2) s : K with s = alpha * (*J5) and complete the
     square. The on-shell ratio is then a function of alpha only.
  5. Measure alpha from a Weyl (chiral) representation that the solver
     does not use, via sympy Matrix Clifford.

Tries to refute C1. Writes ../data/m9_1_audit.json.
"""

from __future__ import annotations

import itertools
import json
import os

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))

# ---------------------------------------------------------------------------
# 1. Levi-Civita algebra
# ---------------------------------------------------------------------------


def eps_tensor() -> np.ndarray:
    e = np.zeros((4, 4, 4, 4), dtype=np.int8)
    for perm in itertools.permutations(range(4)):
        sign = 1
        seq = list(perm)
        for i in range(4):
            for j in range(i + 1, 4):
                if seq[i] > seq[j]:
                    sign *= -1
        e[perm] = sign
    return e


def dual_square_identity() -> float:
    """(*J)_{lmn} (*J)^{lmn} / (J·J) in Euclidean component algebra
    with eps_0123 = +1 and no metric signs (the 3-index contraction).

    eps_{lmn r} eps^{lmn s} = 3! delta_r^s = 6 delta_r^s.
    So (*J)·(*J) = 6 J·J when indices are contracted the same way.
    """
    e = eps_tensor()
    # contract three indices: C[r,s] = eps[l,m,n,r] eps[l,m,n,s]
    c = np.einsum("lmnr,lmns->rs", e, e)
    # should be 6 I
    return float(c[0, 0])


def sympy_eps_contraction() -> int:
    """Same identity in sympy, as a second encoding."""
    e = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    # use sympy LeviCivita
    acc = [[0, 0, 0, 0] for _ in range(4)]
    for l, m, n, r, s in itertools.product(range(4), repeat=5):
        acc[r][s] += int(sp.LeviCivita(l, m, n, r) * sp.LeviCivita(l, m, n, s))
    return acc[0][0]


# ---------------------------------------------------------------------------
# 2. Palatini 4-form -> e R / (2 kappa)
# ---------------------------------------------------------------------------


def palatini_four_form_factor() -> dict:
    """On a constant orthonormal frame, the 4-form
        (1/4) eps_{abcd} e^a /\\ e^b /\\ R^{cd}
    equals (1/2) R * vol   if R^{cd} = (1/2) R^{cd}_{mu nu} dx^mu /\\ dx^nu
    and the curvature scalar is R = R^{cd}_{cd} (orthonormal).

    Count: eps_{abcd} e^a /\\ e^b /\\ ((1/2) R^{cd}_{mu nu} dx^mu /\\ dx^nu)
    with e^a = dx^a gives
        (1/2) eps_{abcd} R^{cd}_{mu nu}  * (dx^a /\\ dx^b /\\ dx^mu /\\ dx^nu)
    = (1/2) eps_{abcd} R^{cd}_{mu nu} eps^{ab mu nu} vol
    = (1/2) * 2 * 2!  R^{cd}_{cd} vol   wait, compute with sympy.
    """
    # Numerical: random antisymmetric R^{cd} = -R^{dc}, Rscalar = R^{cd}_{cd}
    # using orthonormal identification R^{cd}_{mn} as a 4x4x4x4 array
    # with R^{cd}_{mn} antisym in cd and mn, and Riemann symmetries dropped.
    rng = np.random.default_rng(7)
    r_up = rng.normal(size=(4, 4, 4, 4))
    r_up = 0.5 * (r_up - r_up.transpose(1, 0, 2, 3))
    r_up = 0.5 * (r_up - r_up.transpose(0, 1, 3, 2))
    e = eps_tensor().astype(np.float64)
    # four-form coefficient of vol: (1/2) eps_{abcd} R^{cd}_{mn} eps^{abmn}
    four = 0.5 * float(np.einsum("abcd,cdmn,abmn->", e, r_up, e))
    r_scalar = float(np.einsum("cdcd->", r_up))
    # expected: four = 2 * R_scalar   if (1/4)eps e e R = (1/2) R vol
    # because S_g density = (1/kappa) * (1/4 eps e e R) = R/(2 kappa)
    # so (1/4)*4form_without_1/4  wait.
    # We computed (1/2) eps_abcd R^{cd}_mn eps^{abmn} = the 4-form
    #   eps_abcd e^a/\\e^b/\\((1/2)R^{cd}_mn dx^m/\\dx^n)  / vol
    # Palatini density without 1/kappa is (1/4) of that 4-form.
    palatini_density = 0.25 * four
    return {
        "four_form_over_vol": four,
        "R_scalar": r_scalar,
        "palatini_density": palatini_density,
        "palatini_over_R": palatini_density / r_scalar if r_scalar else None,
        "claimed_eR_over_2": 0.5,
    }


# ---------------------------------------------------------------------------
# 3. Q reduced on K = *v
# ---------------------------------------------------------------------------


def q_on_dual_vector() -> dict:
    """K_abc = eps_abcr v^r. Evaluate
       Q = K^{mu c}_mu K_c^{nu}_nu - K^{mu c nu} K_{c mu nu}
    in Euclidean components (no metric signs) and report Q / (v·v).
    The Lorentzian case only multiplies by the same eta used in J·J,
    which cancels in the ratio r.
    """
    e = eps_tensor().astype(np.float64)
    rng = np.random.default_rng(11)
    v = rng.normal(size=4)
    k = np.einsum("abcr,r->abc", e, v)
    # Euclidean: raise = lower
    k_mu_c_mu = np.einsum("mcm->c", k)
    k_c_tr = np.einsum("cnn->c", k)
    term1 = float(np.dot(k_mu_c_mu, k_c_tr))
    term2 = float(np.einsum("mcn,cmn->", k, k))
    q = term1 - term2
    v2 = float(np.dot(v, v))
    return {"Q": q, "v2": v2, "Q_over_v2": q / v2, "K2": float(np.einsum("abc,abc->", k, k))}


# ---------------------------------------------------------------------------
# 4. Complete the square for general alpha
# ---------------------------------------------------------------------------


def on_shell_ratio(alpha: float, q_over_k2: float) -> float:
    """s = alpha * (*J), K = *v, L = Q/(2 kappa) + (1/2) s:K.

    (*J):(*K) = 6 J·v   (Euclidean eps identity),
    s:K = alpha * 6 J·v.
    Stationarity in v, Euclidean: complete the square.

    Let (*J)·(*J) = 6 J·J, (*v)·(*v) = 6 v·v, (*J)·(*v) = 6 J·v.
    K2 = k2_over_v2 * v·v.
    Need the pairing s:K in terms of J, v.

    s_{abc} = alpha * eps_abcr J^r,  K_abc = eps_abcr v^r
    s:K := s^{abc} K_{abc} = alpha * 6 J·v   (Euclidean).

    L = Q/(2k) + (1/2) s:K = (q_over_v2) v·v / (2k) + 3 alpha J·v

    dL/dv = (q_over_v2)/k  v  +  3 alpha J  = 0
    v = - 3 alpha k / q_over_v2  * J

    L_on = (q_over_v2)/(2k) * (9 alpha^2 k^2 / q_over_v2^2) J·J
           + 3 alpha J·(-3 alpha k / q_over_v2 J)
         = (9 alpha^2 k / (2 q_over_v2)) J·J
           - (9 alpha^2 k / q_over_v2) J·J
         = - (9 alpha^2 k /(2 q_over_v2)) J·J

    r = L / (-k J·J) = 9 alpha^2 / (2 q_over_v2)
    """
    return 9.0 * alpha * alpha / (2.0 * q_over_k2)


# ---------------------------------------------------------------------------
# 5. Measure alpha in a Weyl representation (solver uses Dirac / custom plus)
# ---------------------------------------------------------------------------


def weyl_alpha() -> dict:
    """Chiral mostly-minus Weyl matrices, sympy, extract s / (*J)."""
    I2 = sp.eye(2)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.Matrix([[1, 0], [0, -1]])
    zero = sp.zeros(2)

    def blk(a, b, c, d):
        return a.row_join(b).col_join(c.row_join(d))

    # Weyl / chiral mostly minus:
    # g0 = [[0, I], [I, 0]], gi = [[0, s], [-s, 0]]
    g0 = blk(zero, I2, I2, zero)
    gi = [blk(zero, s, -s, zero) for s in (sx, sy, sz)]
    gamma = [g0, *gi]
    g5 = sp.I * gamma[0] * gamma[1] * gamma[2] * gamma[3]
    eta = [1, -1, -1, -1]
    gl = [eta[m] * gamma[m] for m in range(4)]
    gab = [
        [((gl[a] * gl[b] - gl[b] * gl[a]) / 2) for b in range(4)] for a in range(4)
    ]

    rng = np.random.default_rng(13)
    raw = rng.normal(size=4) + 1j * rng.normal(size=4)
    psi = sp.Matrix([sp.Integer(0)] * 4)
    # use floats inside sympy via Float
    psi = sp.Matrix(
        [
            sp.Float(float(raw[i].real)) + sp.I * sp.Float(float(raw[i].imag))
            for i in range(4)
        ]
    )
    A = gamma[0]
    pbar = (psi.H * A).row(0)  # 1x4

    def bilinear(mat):
        val = (pbar * mat * psi)[0]
        return complex(sp.N(val))

    j5 = np.array([bilinear(g5 * gamma[m]).real for m in range(4)], dtype=np.float64)
    # s^mu_ab = (i/4) pbar {g^mu, gab} psi
    s = np.zeros((4, 4, 4), dtype=np.float64)
    for mu in range(4):
        for a in range(4):
            for b in range(4):
                anticom = gamma[mu] * gab[a][b] + gab[a][b] * gamma[mu]
                val = (1j / 4.0) * bilinear(anticom)
                s[mu, a, b] = float(np.real(val))
    # raise ab with eta (diag)
    s_high = np.zeros_like(s)
    for mu in range(4):
        for a in range(4):
            for b in range(4):
                s_high[mu, a, b] = eta[a] * eta[b] * s[mu, a, b]
    e = eps_tensor().astype(np.float64)
    # (*J)^{lmn} = eps^{lmn r} J_r, eps^{0123} = det(eta) eps_0123 = -1
    eps_high = -e
    j_low = np.array([eta[m] * j5[m] for m in range(4)])
    dual = np.einsum("lmnr,r->lmn", eps_high, j_low)
    mask = np.abs(dual) > 1e-10
    factor = float(np.mean(s_high[mask] / dual[mask]))
    anti = (
        s_high
        - s_high.transpose(0, 2, 1)
        - s_high.transpose(1, 0, 2)
        + s_high.transpose(1, 2, 0)
        + s_high.transpose(2, 0, 1)
        - s_high.transpose(2, 1, 0)
    ) / 6.0
    mixed = float(np.max(np.abs(s_high - anti)))
    return {
        "alpha_measured": factor,
        "paper_alpha": -0.25,
        "mixed_max": mixed,
        "j5": j5.tolist(),
    }


def main() -> int:
    six = dual_square_identity()
    six_sy = sympy_eps_contraction()
    pal = palatini_four_form_factor()
    qinfo = q_on_dual_vector()
    weyl = weyl_alpha()

    # Q_over_v2 from the dual parameterization.
    q_over_v2 = qinfo["Q_over_v2"]
    alpha_paper = -0.25
    alpha_weyl = weyl["alpha_measured"]
    r_paper = on_shell_ratio(alpha_paper, q_over_v2)
    r_weyl = on_shell_ratio(alpha_weyl, q_over_v2)
    r_if_q_wrong = on_shell_ratio(alpha_weyl, 2.0 * q_over_v2)  # mutated Palatini weight

    claims = []

    # C1: does the Weyl-measured alpha produce 3/16?
    c1_num = r_weyl
    c1_ok = abs(c1_num - 0.1875) < 1e-8
    claims.append(
        {
            "id": "C1",
            "verdict": "CONFIRMED" if c1_ok else "REFUTED",
            "number": c1_num,
            "note": (
                f"complete-the-square with Weyl alpha={alpha_weyl:.6f} "
                f"and Q/v^2={q_over_v2:.6f} gives r={c1_num:.8f}"
            ),
        }
    )

    c2_ok = abs(alpha_weyl + 0.25) < 1e-8
    claims.append(
        {
            "id": "C2",
            "verdict": "CONFIRMED" if c2_ok else "REFUTED",
            "number": alpha_weyl,
            "note": (
                "Weyl representation, s = (i/4) psibar {g^mu, g_ab} psi; "
                f"paper prints -1/4, measured {alpha_weyl}"
            ),
        }
    )

    c3_ok = weyl["mixed_max"] < 1e-10
    claims.append(
        {
            "id": "C3",
            "verdict": "CONFIRMED" if c3_ok else "REFUTED",
            "number": weyl["mixed_max"],
            "note": "mixed-symmetry residual of the Weyl spin tensor",
        }
    )

    # C4: Euclidean complete-the-square has no leftover signature once
    # J·J uses the same eta. QUALIFIED: this script did not rebuild the
    # mostly-plus Clifford.
    claims.append(
        {
            "id": "C4",
            "verdict": "QUALIFIED",
            "number": None,
            "note": (
                "ratio formula r = 9 alpha^2 / (2 Q/v^2) is algebraic in "
                "the eps identities and cancels the metric used to form J·J; "
                "this script did not re-run a second Clifford representation"
            ),
        }
    )

    c5_moved = abs(r_if_q_wrong - 0.1875) > 0.05
    claims.append(
        {
            "id": "C5",
            "verdict": "CONFIRMED" if c5_moved else "REFUTED",
            "number": r_if_q_wrong,
            "note": "doubling Q/v^2 (Palatini weight) moves r off 3/16",
        }
    )

    # Hunt: if someone uses paper alpha -1/4 in the square with this Q
    r_with_paper_alpha = r_paper

    payload = {
        "task": "m9.1_audit",
        "method_summary": (
            "Levi-Civita 3-index identity + Palatini 4-form factor + "
            "complete-the-square on K=*v; alpha measured in a Weyl Clifford "
            "(sympy), not imported from the solver"
        ),
        "eps_lll_contraction": six,
        "eps_lll_contraction_sympy": six_sy,
        "palatini_four_form": pal,
        "Q_on_dual": qinfo,
        "weyl": weyl,
        "ratio_with_weyl_alpha": r_weyl,
        "ratio_with_paper_alpha": r_with_paper_alpha,
        "ratio_mutated_palatini": r_if_q_wrong,
        "claims": claims,
        "overall": (
            "C1 stands only if the Weyl alpha is the one used on-shell; "
            "C2 is refuted under the 2 dL/domega definition"
        ),
        "what_would_refute_c1": (
            "a Palatini 4-form that is not eR/(2 kappa), or an on-shell "
            "completion that drops L_D, or using paper alpha=-1/4 inside "
            "this note's s definition without changing Cartan"
        ),
    }
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, "m9_1_audit.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, default=str)
        fh.write("\n")
    print("wrote <repo>/openwave/xperiments/m9_emergent_gravity/research/data/m9_1_audit.json")
    print(f"eps contraction = {six} (sympy {six_sy})")
    print(f"palatini_density / R = {pal['palatini_over_R']}")
    print(f"Q/v^2 = {q_over_v2}")
    print(f"Weyl alpha = {alpha_weyl}")
    print(f"r(Weyl) = {r_weyl}   r(paper alpha) = {r_paper}")
    for c in claims:
        print(f"  {c['id']}: {c['verdict']}  {c['number']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
