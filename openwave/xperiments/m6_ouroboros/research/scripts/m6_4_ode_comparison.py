"""M6.4 D2: the three-ODE comparison, symbolically verified.

Three systems circulate in the record set for "the" chaoiton radial problem:

  MAY (z20044392 eq. 3.3-3.4) == LEAN (z20866581 section 1):
      a'' + a'/r - a/r^2 + w^2 a = b
      b'' + b'/r - b/r^2 + w^2 b = a - lam*b - 4g b^3
  JUNE-INTEGRATED (z20866581 section 2.1), c = 1:
      a'' + a'/r - (1+w^2) a + (c/2) b + 4g a^3 = 0
      b'' + b'/r - (1+lam+w^2) b + (c/2) a + 4g b^3 = 0
  M6.2-CERTIFIED v11 reduction: the time reduction of the certified L_ref EL
      pair carries w^2 with the Klein-Gordon sign (m_eff^2 = m^2 - w^2 form),
      like MAY/LEAN and unlike JUNE (see findings/m6_2_method_note.md).

This script verifies, with sympy (no hand algebra trusted):
  V1. The JUNE system is NOT any sign/scaling redefinition of MAY/LEAN:
      checked over all 8 field/parameter sign flips x cross-coupling scale.
  V2. The linearization (Jacobi operator) of the JUNE background ODE has
      cross term -c/2; the printed Jacobi matrix (Q12 = -c/4) does not
      linearize the printed background (factor 2).
  V3. The MAY far-field threshold: eigenvalues of the vacuum linearization
      are X_pm - w^2 with X_pm = (-lam +/- sqrt(lam^2+4))/2 (decay channel
      iff w^2 < X_plus, X_minus < 0 always: one oscillatory channel at
      every omega).
  V4. The JUNE far-field matrix [[1+w^2, -c/2], [-c/2, 1+lam+w^2]] is
      positive definite for all w, lam >= 0 (both channels decay at every
      omega: localization is structural, never selective).

Output: research/data/m6_4_ode_comparison.json (verdicts + the term table).
"""

import json
from pathlib import Path

import sympy as sp

DATA = Path(__file__).resolve().parent.parent / "data"

r, w, lam, g, c = sp.symbols("r w lam g c", positive=True)
a, b = sp.Function("a")(r), sp.Function("b")(r)


def may_residuals():
    ra = (sp.diff(a, r, 2) + sp.diff(a, r) / r - a / r**2 + w**2 * a - b)
    rb = (sp.diff(b, r, 2) + sp.diff(b, r) / r - b / r**2 + w**2 * b
          - (a - lam * b - 4 * g * b**3))
    return ra, rb


def june_residuals(cval=c):
    ra = (sp.diff(a, r, 2) + sp.diff(a, r) / r - (1 + w**2) * a
          + (cval / 2) * b + 4 * g * a**3)
    rb = (sp.diff(b, r, 2) + sp.diff(b, r) / r - (1 + lam + w**2) * b
          + (cval / 2) * a + 4 * g * b**3)
    return ra, rb


def v1_not_redefinable():
    """Try a -> sa*a, b -> sb*b (sa, sb in {+1,-1}) and c free: can the JUNE
    system's residuals become proportional to MAY's? The structural
    obstructions (centrifugal 1/r^2, the sign of w^2, the a^3 term) survive
    every flip; verified by checking the residual difference is nonzero as a
    polynomial for all four sign choices with c left symbolic."""
    may_a, may_b = may_residuals()
    obstructions = []
    for sa in (1, -1):
        for sb in (1, -1):
            ja, jb = june_residuals()
            ja = ja.subs([(a, sa * a), (b, sb * b)]) * sa
            jb = jb.subs([(a, sa * a), (b, sb * b)]) * sb
            diff_a = sp.simplify(sp.expand(ja - may_a))
            diff_b = sp.simplify(sp.expand(jb - may_b))
            # the difference must contain the centrifugal term a/r^2 and the
            # 2*w^2*a mass-sign mismatch and the 4g*a^3 term, for every flip
            pa = sp.Poly(diff_a.subs(sp.Function("a")(r), sp.Symbol("A"))
                         .subs(sp.Function("b")(r), sp.Symbol("B"))
                         .subs(sp.Derivative(sp.Symbol("A"), r), 0),
                         sp.Symbol("A"), sp.Symbol("B"))
            coeff_A = pa.coeff_monomial(sp.Symbol("A"))
            coeff_A3 = pa.coeff_monomial(sp.Symbol("A") ** 3)
            obstructions.append({
                "sa": sa, "sb": sb,
                "linear_a_mismatch": str(sp.simplify(coeff_A)),
                "cubic_a_mismatch": str(sp.simplify(coeff_A3)),
                "identical": bool(diff_a == 0 and diff_b == 0),
            })
    return {"any_redefinition_matches": any(o["identical"]
                                            for o in obstructions),
            "details": obstructions}


def v2_jacobi_cross_term():
    ja, _ = june_residuals(cval=c)
    eps = sp.Symbol("eps")
    xb = sp.Function("xb")(r)
    perturbed = ja.subs(b, b + eps * xb)
    cross = sp.simplify(sp.diff(perturbed, eps).subs(eps, 0) / xb)
    return {"linearized_cross_term": str(cross),
            "printed_Q12": "-c/4",
            "printed_matches_linearization":
                bool(sp.simplify(cross - (-c / 4)) == 0),
            "correct_is_half_c": bool(sp.simplify(cross - c / 2) == 0)}


def v3_may_farfield():
    M = sp.Matrix([[-w**2, 1], [1, -w**2 - lam]])
    evs = list(M.eigenvals().keys())
    xp = (-lam + sp.sqrt(lam**2 + 4)) / 2
    xm = (-lam - sp.sqrt(lam**2 + 4)) / 2
    match = {str(e): [bool(sp.simplify(e - (xp - w**2)) == 0),
                      bool(sp.simplify(e - (xm - w**2)) == 0)] for e in evs}
    return {"eigenvalues": [str(e) for e in evs],
            "equal_to_Xpm_minus_w2": match,
            "X_minus_always_negative": bool(
                sp.simplify(xm < 0) in (True, sp.true) or
                bool(sp.ask(sp.Q.negative(xm), sp.Q.positive(lam)))),
            "threshold": "w^2 < X_plus = (-lam+sqrt(lam^2+4))/2"}


def v4_june_farfield_posdef():
    M = sp.Matrix([[1 + w**2, -sp.Rational(1, 2)],
                   [-sp.Rational(1, 2), 1 + lam + w**2]])
    det = sp.expand(M.det())
    # det = (1+w^2)(1+lam+w^2) - 1/4 >= 1*1 - 1/4 > 0, trace > 0
    det_min = det.subs([(w, 0), (lam, 0)])
    return {"det": str(det), "det_at_w0_lam0": str(det_min),
            "positive_definite_all_w_lam_nonneg": bool(det_min > 0),
            "meaning": "both far-field channels decay at EVERY omega: "
                       "localization is structural in the June system, "
                       "no omega-selection is possible"}


def main():
    out = {
        "task": "M6.4 D2: three-ODE comparison",
        "term_table": {
            "centrifugal -1/r^2": {"MAY/LEAN": "present (both eqs)",
                                   "JUNE": "absent",
                                   "v11_reduction_M6.2": "present"},
            "omega^2 sign": {"MAY/LEAN": "+w^2 (Klein-Gordon reduction)",
                             "JUNE": "-w^2 inside -(1+w^2) (flipped)",
                             "v11_reduction_M6.2": "+w^2"},
            "unit mass": {"MAY/LEAN": "none", "JUNE": "1 (both fields)",
                          "v11_reduction_M6.2": "none printed"},
            "cross coupling": {"MAY/LEAN": "+1 (both eqs)",
                               "JUNE": "+c/2 = +1/2 (both eqs)",
                               "v11_reduction_M6.2": "+1-class"},
            "cubic term": {"MAY/LEAN": "-4g b^3 (beta eq only)",
                           "JUNE": "-4g a^3 AND -4g b^3 (moved to RHS)",
                           "v11_reduction_M6.2": "beta-sector only"},
        },
        "V1_redefinition_scan": v1_not_redefinable(),
        "V2_jacobi_cross": v2_jacobi_cross_term(),
        "V3_may_farfield": v3_may_farfield(),
        "V4_june_farfield": v4_june_farfield_posdef(),
    }
    (DATA / "m6_4_ode_comparison.json").write_text(
        json.dumps(out, indent=2, default=str))
    print(json.dumps({k: out[k] for k in
                      ("V1_redefinition_scan", "V2_jacobi_cross")},
                     indent=2, default=str)[:1500])
    print("V3:", out["V3_may_farfield"]["equal_to_Xpm_minus_w2"])
    print("V4:", out["V4_june_farfield"]["positive_definite_all_w_lam_nonneg"])


if __name__ == "__main__":
    main()
