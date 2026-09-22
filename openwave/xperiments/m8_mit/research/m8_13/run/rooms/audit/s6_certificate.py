"""Exact certificates that close both completeness questions.

MAX:  C8 (exact identity) + |a00|^2 <= 1/7 + TrNbar^2 <= 171/2 + |f|^2 >= 0
      => rhat6 <= 463/924, with a forced equality chain ending in a 2-dim
      eigenspace that is solved exactly.

MIN:  lambda_min of the Sym^2 form is exactly 1/924, so equality holds iff
      u (x) u lies in the V6 isotypic component.  That locus is the common zero
      set of the 15 quadrics spanning V0+V2+V4, which is shown to be the span of
      the 2x2 Hankel minors; rank<=1 Hankel forces a perfect 6th power, i.e. a
      coherent state.
"""
import pathlib
import sys
import itertools

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np
import sympy as sp
from sympy.physics.wigner import clebsch_gordan
import core as C

HERE = pathlib.Path(__file__).parent
out = []
P = out.append
rng = np.random.default_rng(99)

# ---------- the arithmetic of the bound, exactly ---------------------------
bound = -sp.Rational(5, 231) - 0 + sp.Rational(7, 11) * sp.Rational(1, 7) + sp.Rational(171, 2) / 198
P("MAX bound arithmetic (each input justified separately):")
P(f"  -5/231 - 0 + (7/11)(1/7) + (171/2)/198 = {bound} = {sp.nsimplify(bound*924)}/924")
P(f"  equals 463/924 ? {sp.simplify(bound - sp.Rational(463,924)) == 0}")
minbound = -sp.Rational(5, 231) - sp.Rational(9, 22) + 0 + sp.Rational(171, 2) / 198
P(f"  the coherent-state value -5/231 - 9/22 + 0 + (171/2)/198 = {minbound} "
  f"= {sp.nsimplify(minbound*924)}/924")

# ---------- |a00|^2 <= 1/7 : Theta is antiunitary --------------------------
mx = 0.0
for _ in range(500):
    u = rng.normal(size=7) + 1j * rng.normal(size=7)
    u /= np.linalg.norm(u)
    mx = max(mx, abs(np.linalg.norm(C.Theta(u)) - 1))
P("")
P(f"||Theta u|| = ||u|| max deviation over 500 random unit u : {mx:.3e}")
P("  => |<Theta u,u>| <= 1 by Cauchy-Schwarz, so |a00|^2 <= 1/7, equality iff Theta u = lambda u")
# rephasing lemma
P("  rephasing: if Theta u = lambda u with |lambda| = 1 then u' = e^{i phi} u has "
  "Theta u' = e^{-2 i phi} lambda u'; choose e^{2 i phi} = lambda to get Theta u' = u'.")
ok = True
for _ in range(200):
    g = rng.uniform(0, 2 * np.pi)
    u = (np.exp(1j * g) * C.basis(3) + np.exp(-1j * g) * C.basis(-3)) / np.sqrt(2)
    lam = np.vdot(C.Theta(u), u)
    ok &= abs(abs(lam) - 1) < 1e-12
P(f"  spot check |<Theta u,u>| = 1 on the cat orbit : {ok}")

# ---------- the forced equality chain, step 5 ------------------------------
P("")
P("Equality step: u must lie in the top eigenspace of A at (a,b) = (sqrt6/2, 0).")
a, b = sp.symbols('a b', real=True)
Jx, Jy, Jz = C.JOPS_EXACT
A0 = sp.expand((-sp.sqrt(6) / 6) * Jx * Jx + (-sp.sqrt(6) / 6) * Jy * Jy + (sp.sqrt(6) / 3) * Jz * Jz)
A0 = sp.Matrix(7, 7, lambda i, j: sp.nsimplify(sp.simplify(A0[i, j])))
lam = 5 * sp.sqrt(6) / 2
ev = A0.eigenvects()
for val, mult, vecs in ev:
    if sp.simplify(val - lam) == 0:
        P(f"  top eigenvalue {sp.nsimplify(val)} multiplicity {mult}; eigenvectors:")
        for v in vecs:
            P(f"     {sp.simplify(v).T.tolist()}")
P("  => the top eigenspace is span{v3, v-3}.")

al, be = sp.symbols('alpha beta')
u = [al, 0, 0, 0, 0, 0, be]
# <f_z> = 0
fz = sp.simplify(sum(sp.conjugate(u[i]) * Jz[i, j] * u[j] for i in range(7) for j in range(7)))
P(f"  <f_z> for u = alpha v3 + beta v-3 : {sp.simplify(fz)}")
fx = sp.simplify(sum(sp.conjugate(u[i]) * Jx[i, j] * u[j] for i in range(7) for j in range(7)))
fy = sp.simplify(sum(sp.conjugate(u[i]) * Jy[i, j] * u[j] for i in range(7) for j in range(7)))
P(f"  <f_x>, <f_y>                      : {fx}, {fy}  (identically zero: f_pm shift m by 1, 3 and -3 differ by 6)")
Tu = C.Theta_exact(u)
P(f"  Theta u                           : {Tu.T.tolist()}")
P(f"  Theta u = u  <=>  alpha = conj(beta) and beta = conj(alpha)")
P("  with |alpha| = |beta| (from <f_z>=0) and |alpha|^2+|beta|^2 = 1:")
P("     alpha = e^{i gamma}/sqrt2, beta = e^{-i gamma}/sqrt2")
g = sp.symbols('gamma', real=True)
ucat = sp.Matrix([sp.exp(sp.I * g) / sp.sqrt(2), 0, 0, 0, 0, 0, sp.exp(-sp.I * g) / sp.sqrt(2)])
P(f"  Theta(that) - itself              = "
  f"{sp.simplify(C.Theta_exact(list(ucat)) - ucat).T.tolist()}")
P(f"  rhat6 of that (exact, any gamma)  = {sp.simplify(C.rhat6_exact(list(ucat)))}")
# and it IS in the orbit: D_z(-gamma/3) applied to the cat state
P("  D_z(theta) v_m = e^{-i m theta} v_m, so D_z(-gamma/3) (v3+v-3)/sqrt2 "
  "= (e^{i gamma} v3 + e^{-i gamma} v-3)/sqrt2  -> in the rotation orbit.")
mxr = 0.0
for _ in range(200):
    gg = rng.uniform(0, 2 * np.pi)
    lhs = (np.exp(1j * gg) * C.basis(3) + np.exp(-1j * gg) * C.basis(-3)) / np.sqrt(2)
    D = C.rot([0, 0, 1], -gg / 3)
    rhs = D @ (np.array([1, 0, 0, 0, 0, 0, 1], dtype=complex) / np.sqrt(2))
    mxr = max(mxr, np.abs(lhs - rhs).max())
P(f"  numeric check of that rotation identity : {mxr:.3e}")

# the three contact e's are one SO(3) orbit (cyclic permutation of axes)
Rc = np.array([[0, 0, 1.], [1, 0, 0], [0, 1, 0]])
P(f"  cyclic axis permutation matrix det = {np.linalg.det(Rc):+.0f} -> it IS in SO(3), "
  f"so the three contact points are a single rotation orbit; WLOG one of them.")

# ---------- MIN: the V6 locus is the Veronese cone -------------------------
P("")
P("MIN completeness: characterise {u : u (x) u lies in the V6 component}")
# projectors onto isotypic components of Sym^2
pairs = [(i, j) for i in range(7) for j in range(i, 7)]
S = np.zeros((49, 28), dtype=complex)
for k, (i, j) in enumerate(pairs):
    if i == j:
        S[i * 7 + j, k] = 1.0
    else:
        S[i * 7 + j, k] = S[j * 7 + i, k] = 1 / np.sqrt(2)


def irrep_vec(K, Qv):
    v = np.zeros(49, dtype=complex)
    for m1 in C.MS:
        m2 = Qv - m1
        if m2 in C.IDX:
            v[C.IDX[m1] * 7 + C.IDX[m2]] = float(sp.N(clebsch_gordan(3, 3, K, m1, m2, Qv), 30))
    return v


low = []
for K in (0, 2, 4):
    for Qv in range(-K, K + 1):
        v = irrep_vec(K, Qv)
        if np.linalg.norm(v) > 1e-12:
            low.append(v / np.linalg.norm(v))
LOW = np.array(low).T                    # 49 x 15
P(f"  dim of V0+V2+V4 inside Sym^2 : {np.linalg.matrix_rank(S.conj().T @ LOW, tol=1e-9)} (expect 15)")

# coherent states: check a_k = c_{k-3}/sqrt(binom(6,k)) is a perfect power pattern
binom = [float(sp.binomial(6, k)) for k in range(7)]


def acoef(c):
    return np.array([c[C.IDX[k - 3]] / np.sqrt(binom[k]) for k in range(7)])


worst = 0.0
for _ in range(200):
    D = C.rot(rng.normal(size=3), rng.uniform(0, 2 * np.pi))
    aa = acoef(D @ C.basis(3))
    for i, j, k, l in itertools.product(range(7), repeat=4):
        if i + j == k + l:
            worst = max(worst, abs(aa[i] * aa[j] - aa[k] * aa[l]))
P(f"  Hankel minors a_i a_j - a_k a_l (i+j=k+l) on 200 random coherent states: "
  f"max |minor| = {worst:.3e}")

# span of the Hankel minors as quadrics on C^7, vs the V0+V2+V4 quadrics
def quad_from_sym(phi):
    """phi in C^49 -> the quadratic form u -> sum phi_ac u_a u_c, as a 28-vector."""
    Msym = phi.reshape(7, 7)
    Msym = (Msym + Msym.T) / 2
    return np.array([Msym[i, j] * (1 if i == j else 2) for (i, j) in pairs])


qlow = np.array([quad_from_sym(LOW[:, t]) for t in range(LOW.shape[1])])
scale = np.array([1 / np.sqrt(binom[k]) for k in range(7)])
mins = []
for i, j, k, l in itertools.product(range(7), repeat=4):
    if i + j == k + l and (i, j) < (k, l):
        Mq = np.zeros((7, 7))
        Mq[C.IDX[i - 3], C.IDX[j - 3]] += scale[i] * scale[j] / 2
        Mq[C.IDX[j - 3], C.IDX[i - 3]] += scale[i] * scale[j] / 2
        Mq[C.IDX[k - 3], C.IDX[l - 3]] -= scale[k] * scale[l] / 2
        Mq[C.IDX[l - 3], C.IDX[k - 3]] -= scale[k] * scale[l] / 2
        mins.append(np.array([Mq[p, q] * (1 if p == q else 2) for (p, q) in pairs]))
mins = np.array(mins)
P(f"  number of Hankel 2x2 minors generated : {len(mins)}")
P(f"  rank of their span                    : {np.linalg.matrix_rank(mins, tol=1e-9)} (expect 15)")
both = np.vstack([qlow, mins])
P(f"  rank of [V0+V2+V4 quadrics ; Hankel minors] = {np.linalg.matrix_rank(both, tol=1e-9)}"
  f"  -> the two 15-dim spaces COINCIDE if this is 15")
P("  rank<=1 Hankel => a_k = a_0 t^k (or a = (0,...,0,a6)) => the binary sextic is l^6")
P("  => u is a perfect 6th power => u is a coherent state (SU(2) is transitive on P^1).")

# direct numerical confirmation: solve P_{V0+V2+V4}(u x u) = 0 from random starts
from scipy.optimize import least_squares
Plow = LOW @ LOW.conj().T


def resid(x):
    u = x[:7] + 1j * x[7:]
    U = np.kron(u, u)
    r = Plow @ U
    return np.concatenate([r.real, r.imag, [np.vdot(u, u).real - 1]])


bad, nsol, worstgap = None, 0, 0.0
from s5_orbits import canon_min
for _ in range(600):
    x0 = rng.normal(size=14)
    x0 /= np.linalg.norm(x0)
    r = least_squares(resid, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    if np.max(np.abs(r.fun)) < 1e-9:
        nsol += 1
        u = r.x[:7] + 1j * r.x[7:]
        gap = canon_min(u)
        if gap > worstgap:
            worstgap, bad = gap, u.copy()
P(f"  independent solve of P_low(u x u)=0 from 600 random starts: {nsol} solutions found, "
  f"worst 'distance to coherent orbit' = {worstgap:.3e}")
if worstgap > 1e-6:
    P(f"    OFF-ORBIT SOLUTION: {np.round(bad,10)}  rhat6 = {C.rhat6(bad)!r}")

txt = "\n".join(map(str, out))
print(txt)
(HERE / "out_s6.txt").write_text(txt + "\n")
