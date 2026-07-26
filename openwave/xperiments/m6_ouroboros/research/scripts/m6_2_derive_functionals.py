"""M6.2 step D1-D3: derive H, Q, and the constrained EL system symbolically.

SYMBOLIC ONLY: this script outputs functionals and coefficient tables, never
an H/Q number. It runs BEFORE the pre-registration lock; the numeric gate
(m6_2_hq_decision.py) runs after.

Spec of record (M6.1 convention sheet): L_ref = -1/4 F^2 - 1/4 G^2 + J.A
- g (J.J)^2, mostly-plus signature (-,+,+,+), code units (mu^2 = 1), dynamics
(2.1) box A = J, (2.2) box J = A - 4g(J.J)J with box = d_t^2 - lap.

Reduction of record (production benchmark, archive/sandbox_v8):
azimuthal fields A = phihat a(rho) T_A(t), J = phihat b(rho) T_J(t),
temporal gauge, z-independent, 2D-cylindrical measure rho drho.

Steps:
  D0  Verify the curvilinear identities in Cartesian sympy:
      div(phihat f(rho)) = 0  and  |curl(phihat f(rho))|^2 = (f' + f/rho)^2.
  D1  T^00 of L_ref on the reduction; period-average; for BOTH phase
      conventions: MATCHED (T_A = T_J = cos wt) and PRINTED (cos/sin).
  D2  Noether charge of the U(1) phase rotation of the complexified fields
      (the only symmetry that can yield the "J-charge"; the REAL theory has
      no exact internal U(1): checked by transforming L under joint SO(2)
      rotation of (A, J)).  Coefficient DERIVED from the 1/4 kinetic
      normalization.
  D3  Constrained stationarity delta(H - lam*Q) = 0 for the period-averaged
      functionals; term-by-term comparison against the production ODE pair
      AND against the time-reduction of the pinned dynamics (2.1)/(2.2).

Outputs: data/m6_2_derivation.json + full printed derivation. Runtime ~1 min.
"""

import json
from pathlib import Path

import sympy as sp

OUT = Path(__file__).resolve().parent.parent / "data" / "m6_2_derivation.json"
report = {}

# ------------------------------------------------ D0: curvilinear identities
x, y, z = sp.symbols("x y z", real=True)
rho_c = sp.sqrt(x**2 + y**2)
f = sp.Function("f", positive=True)
phihat = sp.Matrix([-y / rho_c, x / rho_c, 0])
V = phihat * f(rho_c)
div_V = sum(sp.diff(V[i], c) for i, c in enumerate((x, y, z)))
curl_V = sp.Matrix([
    sp.diff(V[2], y) - sp.diff(V[1], z),
    sp.diff(V[0], z) - sp.diff(V[2], x),
    sp.diff(V[1], x) - sp.diff(V[0], y),
])
curl2 = sp.simplify(curl_V.dot(curl_V))
rho = sp.symbols("rho", positive=True)
expected_curl2 = (sp.Derivative(f(rho), rho).doit() + f(rho) / rho) ** 2
d0_div = sp.simplify(div_V) == 0
d0_curl = sp.simplify(curl2.subs(rho_c, rho) - expected_curl2) == 0
report["D0_identities"] = {"div_phihat_f_zero": bool(d0_div),
                           "curl2_equals_fprime_plus_f_over_rho_sq": bool(d0_curl)}

# ---------------------------------------- D1: T^00 on the reduction, averaged
# T^00 = 1/2 Adot^2 + 1/2 |curl A|^2 + 1/2 Jdot^2 + 1/2 |curl J|^2
#        - J.A + g (J.J)^2          [mostly-plus, temporal gauge; standard]
t_s, w = sp.symbols("t omega", positive=True)
a = sp.Function("alpha")(rho)
b = sp.Function("beta")(rho)
ap, bp = sp.Derivative(a, rho).doit(), sp.Derivative(b, rho).doit()

def T00(TA, TJ):
    TAd, TJd = sp.diff(TA, t_s), sp.diff(TJ, t_s)
    return (sp.Rational(1, 2) * a**2 * TAd**2
            + sp.Rational(1, 2) * (ap + a / rho) ** 2 * TA**2
            + sp.Rational(1, 2) * b**2 * TJd**2
            + sp.Rational(1, 2) * (bp + b / rho) ** 2 * TJ**2
            - a * b * TA * TJ
            + sp.Symbol("g", positive=True) * b**4 * TJ**4)

def period_avg(expr):
    T = 2 * sp.pi / w
    return sp.simplify(sp.integrate(expr, (t_s, 0, T)) / T)

g_s = sp.Symbol("g", positive=True)
conv = {}
for name, TA, TJ in [("matched", sp.cos(w * t_s), sp.cos(w * t_s)),
                     ("printed", sp.cos(w * t_s), sp.sin(w * t_s))]:
    avg = sp.expand(period_avg(T00(TA, TJ)))
    conv[name] = avg
    report[f"D1_T00_avg_{name}"] = str(avg)
# note: in <H> = Integral(<T00>) rho drho, the cross term 2 a a'/rho inside
# (a' + a/rho)^2 integrates to a boundary term (zero for localized fields):
cross_int = sp.integrate(2 * sp.Function("alpha")(rho)
                         * sp.Derivative(sp.Function("alpha")(rho), rho), rho)
report["D1_gradient_cross_term"] = (
    "int 2 a a'/rho * rho drho = a^2 | boundary = 0 for localized fields; "
    "so int (a'+a/rho)^2 rho drho = int (a'^2 + a^2/rho^2) rho drho")

# ------------------------------- D2: is there an exact internal U(1)? Noether Q
# (i) The REAL theory: joint SO(2) rotation (A,J) -> (A cos th - J sin th,
#     A sin th + J cos th). Kinetic terms invariant; check J.A and (J.J)^2:
th = sp.symbols("theta", real=True)
Av, Jv = sp.symbols("Acomp Jcomp", real=True)  # stand-ins for field values
JA_rot = sp.expand((Av * sp.sin(th) + Jv * sp.cos(th))
                   * (Av * sp.cos(th) - Jv * sp.sin(th)))
report["D2_real_theory_SO2"] = {
    "JA_invariant": bool(sp.simplify(JA_rot - Jv * Av) == 0),
    "conclusion": ("the real two-field theory has NO exact internal U(1): "
                   "J.A is not invariant under joint (A,J) rotation; "
                   "time translation gives H, spatial symmetries give L, P; "
                   "no Noether charge equals the production Q_J"),
}
# (ii) The complexified prescription (the only route to a 'J-charge'):
#      promote both fields to complex, physical field = Re[.], with
#      L_c = -1/4|F|^2 - 1/4|G|^2 + 1/2(J*.A + J.A*) - g |J.J|^2-class quartic,
#      U(1): (A,J) -> e^{i th}(A,J). Noether charge from the 1/4-normalized
#      kinetic term, computed on J = phihat b e^{i w t}, A = phihat a e^{i w t}:
#      pi_J = dL/d(Jdot*) = 1/2 Jdot  (from -1/4 |G|^2 time part 1/2|Jdot|^2)
#      Q = -i * int [pi_J . J* - pi_J* . J + (same for A)] d3x
bt, at_ = sp.symbols("b a", positive=True)
Jc, Ac = bt * sp.exp(sp.I * w * t_s), at_ * sp.exp(sp.I * w * t_s)
piJ, piA = sp.Rational(1, 2) * sp.diff(Jc, t_s), sp.Rational(1, 2) * sp.diff(Ac, t_s)
qdens = sp.simplify(-sp.I * (piJ * sp.conjugate(Jc) - sp.conjugate(piJ) * Jc
                             + piA * sp.conjugate(Ac) - sp.conjugate(piA) * Ac))
report["D2_noether_Q_complexified"] = {
    "Q_density": str(qdens),
    "reading": ("Q = w * int (a^2 + b^2) rho drho  [joint U(1)]; the J-only "
                "rotation is NOT a symmetry even of the complexified theory "
                "(J*.A breaks it), so a J-only charge w*int b^2 exists only "
                "if the coupling term is dropped: DERIVATION STATUS = "
                "convention, not theorem"),
    "production_Q": "int b^2 rho drho (coefficient 1, no w, J-only)",
}

# ---------------- D3: constrained EL of <H> - lam*Q vs the production ODE pair
lam = sp.Symbol("lambda", positive=True)
results_d3 = {}
for name in ("matched", "printed"):
    avg = conv[name]
    # build <H> integrand with the boundary-reduced gradient form
    integrand = sp.expand(avg.rewrite(sp.cos).simplify())
    # replace (a'+a/rho)^2 -> a'^2 + a^2/rho^2 (boundary identity, D1)
    integrand = sp.expand(integrand
                          .subs((ap + a / rho) ** 2, ap**2 + a**2 / rho**2)
                          .subs((bp + b / rho) ** 2, bp**2 + b**2 / rho**2))
    # Q candidates: joint (derived) and J-only (production-shaped)
    for qname, qdens_r in [("Qjoint_w(a2+b2)", w * (a**2 + b**2)),
                           ("QJonly_w_b2", w * b**2)]:
        F_int = (integrand - lam * qdens_r) * rho  # functional density * rho
        el_a = sp.expand(sp.simplify(
            sp.diff(F_int, a) - sp.diff(sp.diff(F_int, ap), rho)))
        el_b = sp.expand(sp.simplify(
            sp.diff(F_int, b) - sp.diff(sp.diff(F_int, bp), rho)))
        results_d3[f"{name}::{qname}"] = {"EL_alpha": str(sp.collect(el_a, rho)),
                                          "EL_beta": str(sp.collect(el_b, rho))}
report["D3_constrained_EL"] = results_d3

# time-reduction of the pinned dynamics (2.1)/(2.2) on the reduction,
# phase-matched e^{iwt} (the only phase assignment satisfying (2.1)):
#   (2.1): (-w^2 - lap_vec) a = b     -> lap_vec a + w^2 a = -b
#   (2.2): (-w^2 - lap_vec) b = a - 4g' b^3 (time-avg cubic: 3/4 factor if
#          real cos^3 projection; exact if complex |J|^2 J convention)
report["D3_time_reduction_of_dynamics"] = {
    "eq_2_1": "lap_vec alpha + w^2 alpha = -beta (phase-matched; sign of beta "
              "absorbable by beta -> -beta jointly in both equations)",
    "eq_2_2_real_cos_projection": "lap_vec beta + w^2 beta = -alpha + 3g beta^3 "
                                  "(cos^3 -> (3/4) cos fundamental-harmonic "
                                  "projection of the 4g beta^3 cubic)",
    "eq_2_2_complex_convention": "lap_vec beta + w^2 beta = -alpha + 4g beta^3 "
                                 "(quartic read as (J*.J)^2-class, no 3/4)",
    "production_ODE": "lap_vec alpha + w^2 alpha = +beta ; "
                      "lap_vec beta + w^2 beta = +alpha - lam*beta - 4g beta^3",
    "notes": ("production matches the complex-convention time reduction up to "
              "the joint beta-sign flip EXCEPT for the extra -lam*beta term, "
              "which does not arise from (2.1)/(2.2) at all: it is the "
              "constrained-minimization multiplier injected into the dynamics; "
              "a lam-term appearing ONLY in the beta equation is inconsistent "
              "with delta(H - lam Q)/delta alpha unless Q is J-only"),
}

# -------- D4: exhaustive convention enumeration, dynamics -> production ODE
# Time-reduce (2.1)/(2.2) on the reduction, phase-matched, over the sign/
# projection conventions, and compare term-by-term to the production pair
#   P1: lap a + w^2 a = +b
#   P2: lap b + w^2 b = +a - lam*b - 4g b^3
# Reduction: box = d_t^2 - lap  ->  (-w^2 - lap) f  for e^{iwt} fields.
#   (2.1) box A = J        ->  lap a + w^2 a = -b                   (always)
#   (2.2) box J = A - 4g (J.J) J,  (J.J) = s_JJ * b^2 (signature reading),
#         cubic projection factor p in {1 (complex |.|^2 reading),
#                                       3/4 (real cos^3 projection)}:
#         -w^2 b - lap b = a - 4 g p s_JJ b^3
#      ->  lap b + w^2 b = -a + 4 g p s_JJ b^3
# Optional joint flip b -> -b (relabeling freedom, applied to BOTH equations).
d4 = {}
for s_JJ in (1, -1):
    for p_name, p_val in (("complex", 1), ("real_proj", sp.Rational(3, 4))):
        for flip in (1, -1):
            # after b -> flip*b: eq1 RHS: -flip*b ; eq2: lap b + w^2 b =
            #   -flip*a + 4 g p s_JJ flip^3 b^3 -> divide flip:
            #   = -a*? ... do it directly:
            # eq1: lap a + w^2 a = -(flip) b
            # eq2 (divide both sides by flip): lap b + w^2 b
            #      = -(1/flip) a + 4 g p s_JJ flip^2 b^3
            c_b_in_eq1 = -flip
            c_a_in_eq2 = -sp.Rational(1, flip)
            c_cubic = 4 * p_val * s_JJ
            match = (c_b_in_eq1 == 1 and c_a_in_eq2 == 1 and c_cubic == -4)
            d4[f"sJJ={s_JJ}, proj={p_name}, flip={flip}"] = {
                "eq1_beta_coeff": str(c_b_in_eq1),
                "eq2_alpha_coeff": str(c_a_in_eq2),
                "eq2_cubic_coeff": f"{c_cubic}*g",
                "production_needs": "+1, +1, -4*g (plus the -lam*b term)",
                "matches_production_sans_lambda": bool(match),
            }
report["D4_convention_enumeration"] = d4
report["D4_conclusion"] = (
    "eq1 needs flip = -1 (joint b -> -b relabeling), giving eq2 alpha "
    "coefficient +1 automatically; the cubic then matches production's -4g "
    "ONLY for s_JJ = -1 with the complex projection, i.e. the mostly-MINUS "
    "reading of (J.J), the OPPOSITE of the signature M6.1 certified from the "
    "printed interaction signs; under the certified mostly-plus reading the "
    "production quartic has the wrong sign relative to the pinned dynamics; "
    "and the -lam*b term arises from NO reading of (2.1)/(2.2)")

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
print(f"\nwritten: {OUT}")
