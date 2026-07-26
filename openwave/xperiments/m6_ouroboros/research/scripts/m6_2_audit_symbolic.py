"""ADVERSARIAL AUDIT M6.2 — independent symbolic derivations for B1-B4.

Method deliberately different from the task's script:
 - B2: build L in temporal gauge explicitly, Legendre transform H = sum pi*qdot - L,
   average over one period by direct integration, and integrate the gradient cross
   term by parts explicitly (symbolic antiderivative check).
 - B1: GENERIC 2x2 internal generator K acting on (A,J); demand dL/deps == 0
   identically in the fields => solve for K. (Task only checked the SO(2) rotation.)
 - B3: EL of <H> - lam*Q via my own variational derivative on the rho-weighted density.
 - B4: brute-force enumeration of signature reading x cubic projection x relabel sign,
   including INDEPENDENT sign flips of a and b (task only used joint flip).
"""
import itertools
import sympy as sp

rho, t, w, g, lam = sp.symbols("rho t omega g lambda", positive=True)
al = sp.Function("alpha")(rho)
be = sp.Function("beta")(rho)
alp = al.diff(rho)
bep = be.diff(rho)

print("=" * 70)
print("B2: period-averaged energy of L_ref on the reduction, matched phases")
print("=" * 70)
# mostly-plus (-,+,+,+), temporal gauge A_0 = J_0 = 0
# -1/4 F^2 = 1/2 |Adot|^2 - 1/2 |curl A|^2   (standard; F^2 = F^{munu}F_{munu})
# J.A = J^mu A_mu = +vec(J).vec(A)   (spatial only, mostly-plus)
# (J.J) = +|vec J|^2 (mostly-plus)
# curl(phihat f(rho)) = zhat (f' + f/rho) for z-independent azimuthal field
TA = sp.cos(w * t)
TJ = sp.cos(w * t)  # matched phases
A_dot2 = (al * TA.diff(t)) ** 2
J_dot2 = (be * TJ.diff(t)) ** 2
curlA2 = ((alp + al / rho) * TA) ** 2
curlJ2 = ((bep + be / rho) * TJ) ** 2
JdotA = al * be * TA * TJ
JJ = (be * TJ) ** 2  # mostly-plus: (J.J) = +beta^2 T^2
L_density = (sp.Rational(1, 2) * A_dot2 - sp.Rational(1, 2) * curlA2
             + sp.Rational(1, 2) * J_dot2 - sp.Rational(1, 2) * curlJ2
             + JdotA - g * JJ**2)
# Legendre: pi_A = dL/dAdot etc. Fields are q = al*TA, qdot = al*TA'.
# H = pi*qdot - L = (Adot^2 + Jdot^2) - L   since L = 1/2 Adot^2 + ... .
H_density = (A_dot2 + J_dot2) - L_density
H_avg = sp.simplify(sp.integrate(H_density, (t, 0, 2 * sp.pi / w)) * w / (2 * sp.pi))
H_avg = sp.expand(H_avg)
print("<T00> =", H_avg)

# claimed integrand after boundary reduction:
claimed = (sp.Rational(1, 4) * (w**2 * al**2 + alp**2 + al**2 / rho**2
                                + w**2 * be**2 + bep**2 + be**2 / rho**2)
           - sp.Rational(1, 2) * al * be + sp.Rational(3, 8) * g * be**4)
diff_density = sp.expand(H_avg - claimed)
print("\n<T00> - claimed integrand =", sp.simplify(diff_density))
# should equal (1/4)(2 al alp / rho + 2 be bep / rho): the boundary cross term.
boundary_part = sp.expand(sp.Rational(1, 2) * (al * alp + be * bep) / rho)
resid = sp.simplify(diff_density - boundary_part)
print("residual after removing (a a' + b b')/(2 rho):", resid)
# check int (a a'/rho) * rho drho = int a a' drho = a^2/2 |boundary
anti = sp.integrate(al * alp, rho)
print("antiderivative of a*a' :", anti, " -> boundary term, 0 for localized fields")
assert resid == 0

print()
print("=" * 70)
print("B2b: is production H_p the energy under ANY normalization k of -kF^2?")
print("=" * 70)
# General family consistent with the printed EL pair (2.1)/(2.2):
#   L_k = -k F^2 - k G^2 + 4k J.A - 4k g (J.J)^2   (EL: box A = J, box J = A -4g(J.J)J)
k = sp.symbols("k", positive=True)
Lk = (2 * k * A_dot2 - 2 * k * curlA2 + 2 * k * J_dot2 - 2 * k * curlJ2
      + 4 * k * JdotA - 4 * k * g * JJ**2)
Hk = (4 * k * A_dot2 + 4 * k * J_dot2) - Lk
Hk_avg = sp.expand(sp.simplify(sp.integrate(Hk, (t, 0, 2 * sp.pi / w)) * w / (2 * sp.pi)))
print("<H_k> =", Hk_avg)
# production integrand: alp^2 + al^2/rho^2 + bep^2 + be^2/rho^2 + 4 g be^4 (no w^2, no cross)
# electric term coefficient in <H_k>: k*w^2*al^2 -> vanishes only k=0. cross: -2k al be.
coef_elec = Hk_avg.coeff(al**2).subs(rho, sp.oo)  # w^2 * k piece
print("coefficient of w^2*a^2:", sp.expand(Hk_avg).coeff(w**2 * al**2))
print("coefficient of a*b   :", sp.expand(Hk_avg).coeff(al * be))
print("coefficient of b^4   :", sp.expand(Hk_avg).coeff(be**4))
print("=> electric ~ k*w^2 and cross ~ -2k never vanish for k>0;")
print("   quartic 3kg/2: equals 4g only at k=8/3, but then gradients=8/3 not 1. NO k works.")
# instantaneous slices: t=0 kills electric but cross = -4k*al*be present; t=T/4 kills
# gradients entirely. Check both:
H_t0 = sp.expand(Hk.subs(t, 0))
print("H_k(t=0) =", H_t0)
H_tq = sp.expand(Hk.subs(t, sp.pi / (2 * w)))
print("H_k(t=T/4) =", H_tq)

print()
print("=" * 70)
print("B1: generic internal symmetry delta(A,J) = eps*K*(A,J), K arbitrary 2x2")
print("=" * 70)
k11, k12, k21, k22, eps = sp.symbols("k11 k12 k21 k22 epsilon", real=True)
# transform the REDUCED fields (valid test: an internal symmetry must hold on
# every configuration, so holding identically on this family is necessary)
alT = al + eps * (k11 * al + k12 * be)
beT = be + eps * (k21 * al + k22 * be)
alTp = alT.diff(rho)
beTp = beT.diff(rho)


def Lred(a_, ap_, b_, bp_):
    """spatial part of L on the reduction at t=0 (TA=TJ=1) and the time part
    at t=T/4 separately; a symmetry must preserve both pieces since the time
    dependence factorizes."""
    mag = (-sp.Rational(1, 2) * (ap_ + a_ / rho) ** 2
           - sp.Rational(1, 2) * (bp_ + b_ / rho) ** 2
           + a_ * b_ - g * b_**4)
    kin = sp.Rational(1, 2) * w**2 * (a_**2 + b_**2)  # from Adot^2 at t=T/4
    return mag, kin


mag0, kin0 = Lred(al, alp, be, bep)
magT, kinT = Lred(alT, alTp, beT, beTp)
dmag = sp.expand(sp.diff(magT, eps).subs(eps, 0))
dkin = sp.expand(sp.diff(kinT, eps).subs(eps, 0))
# collect coefficients of independent field monomials
monos_m = sp.Poly(dmag, al, be, alp, bep).as_dict()
monos_k = sp.Poly(dkin, al, be, alp, bep).as_dict()
eqs = set()
for d in (monos_m, monos_k):
    for mono, coef in d.items():
        eqs.add(sp.simplify(coef))
sol = sp.solve(list(eqs), [k11, k12, k21, k22], dict=True)
print("first-order invariance conditions solve to K =", sol)
print("=> only K = 0: NO continuous internal symmetry of the real theory. "
      "In particular no Noether charge ~ int beta^2.")

print()
print("Noether charge of the complexified joint U(1) (1/4-normalized kinetic):")
# L_c time part per field: 1/2 |qdot|^2 ; q = f e^{iwt}; delta q = i q
b_s, a_s = sp.symbols("b a", positive=True)
Qdens = 0
for f in (a_s, b_s):
    q = f * sp.exp(sp.I * w * t)
    dq = sp.diff(q, t)
    # Q = dL/dqdot * (i q) + dL/dqdot* * (-i q*)  with L = 1/2 qdot qdot*
    Qdens += sp.simplify(sp.Rational(1, 2) * sp.conjugate(dq) * (sp.I * q)
                         + sp.Rational(1, 2) * dq * (-sp.I * sp.conjugate(q)))
print("Q density =", sp.simplify(Qdens), "  => Q_joint = w * int (a^2+b^2) rho drho")

print()
print("=" * 70)
print("B3: EL of <H> - lam*Q_joint (and Q_Jonly) at fixed omega")
print("=" * 70)
for qname, Qd in [("Q_joint", w * (al**2 + be**2)), ("Q_Jonly", w * be**2)]:
    dens = (claimed - lam * Qd) * rho
    ELa = sp.expand(sp.simplify(dens.diff(al) - (dens.diff(alp)).diff(rho)))
    ELb = sp.expand(sp.simplify(dens.diff(be) - (dens.diff(bep)).diff(rho)))
    # normalize: multiply by 2/rho to compare with lap form
    ELa_n = sp.expand(2 * ELa / rho)
    ELb_n = sp.expand(2 * ELb / rho)
    print(f"[{qname}] EL_alpha*2/rho:", ELa_n)
    print(f"[{qname}] EL_beta *2/rho:", ELb_n)
print("production ODE: lap a + w^2 a - b = 0 ; lap b + w^2 b - a + lam b + 4g b^3 = 0")
print("i.e. production has +w^2 WITH +lap; the EL above has -lap with +w^2 => opposite")
print("relative sign (energy/Yukawa sign vs eigenvalue/oscillatory sign).")

print()
print("=" * 70)
print("B4: enumeration — time-reduction of (2.1)/(2.2) vs production ODE")
print("=" * 70)
# (2.1) box A = J ; (2.2) box J = A - 4g (J.J) J ; box = d_t^2 - lap
# e^{iwt} matched phases: box -> (-w^2 - lap)
# (J.J) = s * beta^2 (s=+1 mostly-plus certified; s=-1 mostly-minus)
# cubic time projection p: 1 (complex |.|^2 convention) or 3/4 (real cos^3)
# allow INDEPENDENT relabelings a -> ea*a, b -> eb*b (ea,eb = +-1)
rows = []
for s, p, ea, eb in itertools.product((1, -1), (1, sp.Rational(3, 4)), (1, -1), (1, -1)):
    # original: lap a + w^2 a = -b            [from (-w^2-lap)a = b]
    #           lap b + w^2 b = -a + 4 g p s b^3
    # substitute a = ea*A, b = eb*B, divide each eq by its leading sign:
    c_b_eq1 = -eb / ea            # coefficient of B in eq1 RHS
    c_a_eq2 = -ea / eb            # coefficient of A in eq2 RHS
    c_cub = 4 * p * s * eb**2     # coefficient of B^3 in eq2 RHS (eb^3/eb)
    match = (c_b_eq1 == 1) and (c_a_eq2 == 1) and (c_cub == -4)
    rows.append((s, p, ea, eb, c_b_eq1, c_a_eq2, c_cub, match))
    print(f"s_JJ={s:+d} p={p} ea={ea:+d} eb={eb:+d} -> eq1 b-coef {c_b_eq1}, "
          f"eq2 a-coef {c_a_eq2}, cubic {c_cub}g | production needs +1,+1,-4g "
          f"{'<== MATCH (sans lam)' if match else ''}")
n_match = sum(1 for r in rows if r[-1])
print(f"\nmatches: {n_match}; all matches require s_JJ=-1 (mostly-MINUS) & p=1 (complex)")
print("and NO row generates the -lam*b term: it appears in no reading of the dynamics.")
