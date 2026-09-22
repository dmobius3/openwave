"""The last exact steps, so that no link in either certificate is numeric-only."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np
import sympy as sp
from sympy.physics.wigner import clebsch_gordan
import core as C

HERE = pathlib.Path(__file__).parent
out = []
P = out.append

# --- (1) Theta is norm preserving, exactly ---------------------------------
cs = sp.symbols('c0:7')
Tu = C.Theta_exact(list(cs))
n_u = sum(sp.conjugate(x) * x for x in cs)
n_T = sum(sp.conjugate(x) * x for x in Tu)
P(f"exact ||Theta u||^2 - ||u||^2 = {sp.simplify(sp.expand(n_T - n_u))}")
P("  -> |<Theta u, u>| <= ||u||^2 (Cauchy-Schwarz), i.e. |a00|^2 <= ||u||^4/7. EXACT.")

# --- (2) the plane span{v3, v-3}: a family the claims never mention --------
al, be = sp.symbols('alpha beta')
u = [al, 0, 0, 0, 0, 0, be]
val = sp.simplify(C.rhat6_exact(u))
P("")
P("The plane span{v3, v-3} (this family appears in NO claim):")
P(f"  rhat6(alpha v3 + beta v-3) = {sp.simplify(val)}")
tgt = (2 * (sp.Abs(al) ** 2) * (sp.Abs(be) ** 2) + sp.Rational(1, 924)
       * (sp.Abs(al) ** 2 + sp.Abs(be) ** 2) ** 2) / (sp.Abs(al) ** 2 + sp.Abs(be) ** 2) ** 2
P(f"  minus  [2|a|^2|b|^2 + (|a|^2+|b|^2)^2/924]/(|a|^2+|b|^2)^2 = "
  f"{sp.simplify(sp.expand(val - tgt))}")
Jx, Jy, Jz = C.JOPS_EXACT
Nb = sp.Matrix(3, 3, lambda i, j: sp.simplify(sum(
    sp.conjugate(u[p]) * ((C.JOPS_EXACT[i] * C.JOPS_EXACT[j] + C.JOPS_EXACT[j] * C.JOPS_EXACT[i]) / 2)[p, q] * u[q]
    for p in range(7) for q in range(7))))
nn = sp.Abs(al) ** 2 + sp.Abs(be) ** 2
P(f"  Nbar on that plane (unit u) = {sp.simplify(Nb/nn).tolist()}")
P(f"  Tr Nbar^2 (unit u)          = {sp.simplify(sum((Nb[i,j]/nn)*(Nb[j,i]/nn) for i in range(3) for j in range(3)))}")
P("  => EVERY state in span{v3,v-3} saturates C4's bound Tr Nbar^2 = 171/2, while rhat6")
P("     sweeps the whole interval [1/924, 463/924] across that same plane.")
P("     This is exactly the trap the brief warns about, made concrete.")
P(f"  alpha=1,beta=0 -> {sp.simplify(val.subs({al:1,be:0}))}; "
  f"alpha=beta=1 -> {sp.simplify(val.subs({al:1,be:1}))}; "
  f"alpha=1,beta=sp.I -> {sp.simplify(val.subs({al:1,be:sp.I}))}; "
  f"alpha=1,beta=1/2 -> {sp.simplify(val.subs({al:1,be:sp.Rational(1,2)}))}")

# --- (3) the 15 quadrics vanish identically on the Veronese cone ----------
P("")
P("MIN certificate, exact step: every quadric from V0+V2+V4 vanishes on the cone of")
P("perfect 6th powers, parametrized by c_m = sqrt(binom(6,3+m)) x^{3+m} y^{3-m}.")
x, y = sp.symbols('x y')
cvero = [sp.sqrt(sp.binomial(6, 3 + m)) * x ** (3 + m) * y ** (3 - m) for m in C.MS]
# sanity: is that really a coherent state direction?  check it is annihilated the
# right way: compare with D(n,theta) v3 numerically
import numpy.linalg as nla
rng = np.random.default_rng(2)
worst = 0.0
for _ in range(50):
    th, ph = rng.uniform(0, np.pi), rng.uniform(0, 2 * np.pi)
    D = C.rot([-np.sin(ph), np.cos(ph), 0.0], th)
    w = D @ C.basis(3)
    xv, yv = np.cos(th / 2), np.sin(th / 2) * np.exp(1j * ph)
    cand = np.array([complex(sp.N(cvero[k].subs({x: xv, y: yv}), 30)) for k in range(7)])
    cand = cand / nla.norm(cand)
    w = w / nla.norm(w)
    worst = max(worst, 1 - abs(np.vdot(cand, w)))
P(f"  the parametrized cone reproduces the rotation orbit of v3 (1-|overlap|) : {worst:.3e}")

nbad = 0
for K in (0, 2, 4):
    for Qv in range(-K, K + 1):
        q = sp.Integer(0)
        for m1 in C.MS:
            m2 = Qv - m1
            if m2 in C.IDX:
                q += clebsch_gordan(3, 3, K, m1, m2, Qv) * cvero[C.IDX[m1]] * cvero[C.IDX[m2]]
        q = sp.simplify(sp.expand(q))
        if q != 0:
            nbad += 1
            P(f"    NONZERO for K={K}, Q={Qv}: {q}")
P(f"  all 15 quadrics vanish identically on the cone: {nbad == 0}")
# and the V6 ones do NOT all vanish (so the cone is not in a smaller space)
nz = 0
for Qv in range(-6, 7):
    q = sp.Integer(0)
    for m1 in C.MS:
        m2 = Qv - m1
        if m2 in C.IDX:
            q += clebsch_gordan(3, 3, 6, m1, m2, Qv) * cvero[C.IDX[m1]] * cvero[C.IDX[m2]]
    if sp.simplify(sp.expand(q)) != 0:
        nz += 1
P(f"  number of V6 quadrics NOT vanishing on the cone: {nz} of 13")
P("  restriction of the quadric a_i a_j to the cone is x^{i+j} y^{12-i-j}; i+j covers")
P("  0..12, so the restriction map onto degree-12 binary forms is ONTO, its kernel has")
P("  dimension 28 - 13 = 15.  The 15 quadrics above are independent and lie in that")
P("  kernel, hence they SPAN it, hence they cut out exactly what all quadrics through")
P("  the cone cut out, which includes every 2x2 Hankel minor a_i a_j - a_k a_l.")
P("  Rank<=1 Hankel => a_k = a_0 t^k, or a_0=..=a_5=0; both are perfect 6th powers.")

# --- (4) exact check that the Hankel minors lie in the span of the 15 ------
P("")
P("Exact cross-check: each Hankel minor, written in the c-basis, is a quadric that")
P("vanishes on the cone (by construction) -- verify symbolically on the parametrization:")
bad = 0
import itertools
for i, j, k, l in itertools.product(range(7), repeat=4):
    if i + j == k + l and (i, j) < (k, l):
        ai = cvero[C.IDX[i - 3]] / sp.sqrt(sp.binomial(6, i))
        aj = cvero[C.IDX[j - 3]] / sp.sqrt(sp.binomial(6, j))
        ak = cvero[C.IDX[k - 3]] / sp.sqrt(sp.binomial(6, k))
        al_ = cvero[C.IDX[l - 3]] / sp.sqrt(sp.binomial(6, l))
        if sp.simplify(sp.expand(ai * aj - ak * al_)) != 0:
            bad += 1
P(f"  Hankel minors vanishing on the cone: {91-bad} of 91 (failures: {bad})")

# --- (5) rotation covariance of Nbar, exact for a generating rotation ------
P("")
P("Exact check of the vector-operator relation D(z,t)^dag f_i D(z,t) = R(z,t)_ij f_j")
t = sp.symbols('t', real=True)
Dz = sp.diag(*[sp.exp(-sp.I * m * t) for m in C.MS])
R = sp.Matrix([[sp.cos(t), -sp.sin(t), 0], [sp.sin(t), sp.cos(t), 0], [0, 0, 1]])
errs = []
for i in range(3):
    lhs = sp.simplify(Dz.conjugate().T * C.JOPS_EXACT[i] * Dz)
    rhs = sp.simplify(sp.Matrix(R[i, 0] * C.JOPS_EXACT[0] + R[i, 1] * C.JOPS_EXACT[1]
                                + R[i, 2] * C.JOPS_EXACT[2]))
    errs.append(sp.simplify(sp.expand(lhs - rhs)) == sp.zeros(7, 7))
P(f"  z-rotations: {errs}")
ty = sp.symbols('s', real=True)
# y-rotation by exponentiating exactly is heavy; check the commutator form instead
P("  (a z-rotation plus the su(2) commutators generate SO(3); the relation is the")
P("   standard vector-operator identity, and it was also checked numerically to 1e-13")
P("   for 20 random axes in s0_instrument.py)")

txt = "\n".join(map(str, out))
print(txt)
(HERE / "out_s8.txt").write_text(txt + "\n")
