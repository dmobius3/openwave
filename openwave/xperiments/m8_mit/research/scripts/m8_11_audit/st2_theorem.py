"""Stage 2, task 5: checks of THEOREM.md's steps on MY stage-1 objects (out/main_all.pkl, exact MQ).

Parts
  S   the S lemma: derive the section-level action of S = R_z(pi/3) o (conjugation) in the P-trick
      representation, verify it pointwise (numerically) on sections of every level, verify that it
      preserves the sector, and check every consequence exactly on my computed objects.
  B   the bridge: d^2/ds^2 int|cos s Phi + sin s E|^4 at 0 = 4(<E, DN E>_R - Q), exactly, from a
      quartic interpolation of int|c Phi + s E|^4 (computed with N only, not DN).
  L4  the lambda_4 lemma: <Phi, DN_Phi v> = 2X + conj X for block v, X = <N(Phi), v>; X = 0 at kappa;
      sign ingredients (mu_n > 0, Q != 1, sum_{n != 6} ||Pi_n N||^2 = int|Phi|^6 - Q^2 > 0).
  E   E1/E2/T.b: chart actions, orbit direction at U5, no continuous symmetry of K_H at U6,
      Hess_T Q invertible, xi in X_H exactly (finite generators with exact matrices).
Every check prints PASS/FAIL and has a recorded mutant (the same comparison on a wrong object must fail)."""
import json, pickle, random
from fractions import Fraction as F
import numpy as np
from mq import MQ, ZERO, ONE, mq
from su2 import theta, D_exact, q_to_su2
from engine import block, N as Nsec, DN as DNsec, inner, compress, block_fibre, P as PROJ
from numD import D_of

st = pickle.load(open("out/main_all.pkl", "rb"))
GRP = pickle.load(open("out/group.pkl", "rb"))
LOG = []


def check(label, ok, mutant_ok, mutant_desc):
    """ok: result on the real object (must be True); mutant_ok: same test on a mutant (must be False)."""
    if mutant_ok is None:
        line = ("PASS " if ok else "FAIL ") + label + "   [no separate mutant: %s]" % mutant_desc
        caught = None
    else:
        caught = not mutant_ok
        line = ("PASS " if ok else "FAIL ") + label + "   [mutant: %s -> %s]" % (mutant_desc, "caught" if caught else "NOT caught")
    print(line)
    LOG.append({"label": label, "ok": bool(ok), "mutant": mutant_desc, "mutant_caught": caught})
    return ok


I = MQ.I
s3 = MQ.sqrt(3)
half = F(1, 2)


def cis_mpi3(m):
    """exp(i m pi/3) exactly."""
    t = {0: (ONE, ZERO), 1: (mq(half), s3 * half), 2: (mq(-half), s3 * half), 3: (mq(-1), ZERO),
         4: (mq(-half), -(s3 * half)), 5: (mq(half), -(s3 * half))}[m % 6]
    return t[0] + I * t[1]


# rotation r = R_z(pi/3) = diag(e^{-i pi/6}, e^{i pi/6}) in SU(2), exact
r_mat = [[s3 * half - I * half, ZERO], [ZERO, s3 * half + I * half]]


def conjD_r(J, x):
    """conj(D^J(r)) x, computed from my D_exact (not from the diagonal formula)."""
    D = D_exact(J, r_mat)
    return [sum((D[i][k].conj() * x[k] for k in range(2 * J + 1)), ZERO) for i in range(2 * J + 1)]


def S_term(J, x, Ycols):
    """S on a P-trick term x^T D^J(g) Y (Y given by its 7 columns): derived action
       x -> conj(D^J(r)) Theta x,   column b -> (-1)^b Theta(column -b)."""
    x2 = conjD_r(J, theta(x, J))
    Y2 = [None] * 7
    for b in range(-3, 4):
        col = theta(Ycols[-b + 3], J)
        Y2[b + 3] = [c * ((-1) ** (b % 2)) for c in col]
    return x2, Y2


def S_func(f):
    return {J: [S_term(J, x, Y) for x, Y in terms] for J, terms in f.items()}


def coeff(f):
    """canonical coefficient tensor per level: C[J][(m, k, a)] = sum x_m Y_a[k]."""
    out = {}
    for J, terms in f.items():
        C = {}
        for x, Y in terms:
            for m in range(2 * J + 1):
                if x[m].is_zero():
                    continue
                for a in range(7):
                    for k in range(2 * J + 1):
                        if Y[a][k].is_zero():
                            continue
                        key = (m, k, a)
                        C[key] = C.get(key, ZERO) + x[m] * Y[a][k]
        out[J] = {k: v for k, v in C.items() if not v.is_zero()}
    return out


def same(f1, f2, scale=ONE):
    c1, c2 = coeff(f1), coeff(f2)
    for J in set(c1) | set(c2):
        a, b = c1.get(J, {}), c2.get(J, {})
        for k in set(a) | set(b):
            if not (a.get(k, ZERO) - scale * b.get(k, ZERO)).is_zero():
                return False
    return True


def fscale(f, s):
    return {J: [([c * s for c in x], Y) for x, Y in t] for J, t in f.items()}


# ------------------------------------------------------------------------------------------ S lemma
print("=== S lemma: section-level derivation and consequences")
# (S0) the ingredient conj(P) = R P R exactly (R_{m,m'} = (-1)^m delta_{m,-m'}), so (Theta P) R = P
for d in (3, 4):
    P = PROJ[d]
    ok = all((P[i][k].conj() - P[6 - i][6 - k] * ((-1) ** ((i + k) % 2))).is_zero() for i in range(7) for k in range(7))
    Pm = [[P[i][k] * (1 if (i, k) != (0, 0) else 2) for k in range(7)] for i in range(7)]
    okm = all((Pm[i][k].conj() - Pm[6 - i][6 - k] * ((-1) ** ((i + k) % 2))).is_zero() for i in range(7) for k in range(7))
    check("sector %d: conj(P) = R P R exactly (P is a real-coefficient class-sum polynomial)" % d, ok, okm,
          "P with one entry doubled")
    Pcols = [[P[k][a] for k in range(7)] for a in range(7)]
    _, Y2 = S_term(3, [ONE] + [ZERO] * 6, Pcols)
    ok = all((Y2[a][k] - Pcols[a][k]).is_zero() for a in range(7) for k in range(7))
    check("sector %d: S maps the block right factor P to P exactly" % d, ok, None, "(implied by the line above)")

# (S1) pointwise, numerical: the derived action equals (S psi)(g) = conj(psi(r^-1 g)) R on every level,
#      preserves |.|, commutes with N pointwise, and keeps the sector: (S f)(g h) = (S f)(g) D^3(h), h = q1, q2
rng = np.random.default_rng(7)


def rand_su2():
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([[w - 1j * z, -y - 1j * x], [y - 1j * x, w + 1j * z]])


def evalf(f, g):
    v = np.zeros(7, complex)
    for J, terms in f.items():
        D = D_of(J, g)
        for x, Y in terms:
            xn = np.array([c.to_complex() for c in x])
            Yn = np.array([[c.to_complex() for c in col] for col in Y]).T      # (2J+1) x 7
            v += xn @ D @ Yn
    return v


Rm = np.zeros((7, 7))
for m in range(-3, 4):
    Rm[m + 3, -m + 3] = (-1) ** (m % 2)
rn = np.array([[c.to_complex() for c in row] for row in r_mat])
rinv = rn.conj().T
qn = {k: np.array([[c.to_complex() for c in row] for row in q_to_su2(GRP[k])]) for k in ("q1", "q2")}
worst = {"formula": 0.0, "norm": 0.0, "N": 0.0, "sector": 0.0, "mut": 0.0, "mut2": 0.0}
for key in (("U6", 3), ("U6", 4), ("U5", 3)):
    S_ = st[key]
    f = {3: [(S_["u"], [[PROJ[key[1]][k][a] for k in range(7)] for a in range(7)])]}
    for J, t in S_["xi"].items():
        f.setdefault(J, []).extend(t)
    Sf = S_func(f)
    for _ in range(4):
        g = rand_su2()
        lhs = evalf(Sf, g)
        rhs = np.conj(evalf(f, rinv @ g)) @ Rm
        worst["formula"] = max(worst["formula"], np.abs(lhs - rhs).max())
        worst["mut"] = max(worst["mut"], np.abs(lhs - np.conj(evalf(f, rn @ g)) @ Rm).max())   # mutant: r instead of r^-1
        worst["mut2"] = max(worst["mut2"], np.abs(lhs - np.conj(evalf(f, rinv @ g))).max())   # mutant: R dropped
        worst["norm"] = max(worst["norm"], abs(np.linalg.norm(lhs) - np.linalg.norm(evalf(f, rinv @ g))))
        Nf_S = np.linalg.norm(lhs) ** 2 * lhs
        f0 = evalf(f, rinv @ g)
        worst["N"] = max(worst["N"], np.abs(Nf_S - np.conj(np.linalg.norm(f0) ** 2 * f0) @ Rm).max())
        for h in qn.values():
            worst["sector"] = max(worst["sector"], np.abs(evalf(Sf, g @ h) - evalf(Sf, g) @ D_of(3, h)).max())
print("   pointwise deviations:", {k: float(v) for k, v in worst.items()})
check("derived fibre action of S equals conj(psi(r^-1 g)) R pointwise on all levels (dev < 1e-12)",
      worst["formula"] < 1e-12, worst["mut"] < 1e-12 or worst["mut2"] < 1e-12, "r used instead of r^-1 (degenerate at U6, whose stabilizer contains r^2); R dropped")
check("S preserves pointwise |psi| and commutes with N pointwise (dev < 1e-12)",
      worst["norm"] < 1e-12 and worst["N"] < 1e-12, None, "(covered by the formula mutant)")
check("S keeps the sector: (S f)(g h) = (S f)(g) D^3(h) for h = q1, q2 (dev < 1e-12)",
      worst["sector"] < 1e-12, None, "(see the conj(P) mutant)")

# (S2) chart actions, exactly, at generic complex z = 1/3 + 2i/5 and at z0
def u_of(z):
    v = [ZERO] * 7
    v[6] = ONE
    v[0] = ONE
    v[3] = z
    return v


def prop(a, b):
    """exact: a = c b for some scalar c; returns c or None."""
    idx = [i for i in range(7) if not b[i].is_zero()]
    c = a[idx[0]] / b[idx[0]]
    return c if all((a[i] - c * b[i]).is_zero() for i in range(7)) else None


z = mq(F(1, 3)) + I * F(2, 5)
zb = z.conj()
th = theta(u_of(z), 3)
c1 = prop(th, u_of(-zb))
c2 = prop(conjD_r(3, u_of(z)), u_of(-z))
c3 = prop(conjD_r(3, theta(u_of(z), 3)), u_of(zb))
check("Theta on the chart: Theta u(z) = %s * u(-conj z)" % c1, c1 is not None, prop(th, u_of(zb)) is not None,
      "compare with u(conj z)")
check("R_z(pi/3) on the chart: u(z) -> %s * u(-z)" % c2, c2 is not None, prop(conjD_r(3, u_of(z)), u_of(z)) is not None,
      "compare with u(z)")
check("S = R_z(pi/3) o Theta on the chart: u(z) -> %s * u(conj z) (no phase)" % c3, c3 is not None and c3 == ONE,
      prop(conjD_r(3, theta(u_of(z), 3)), u_of(z)) is not None, "compare with u(z)")

# (S3) consequences on my objects, exactly (U6, both sectors)
for d in (3, 4):
    S_ = st[("U6", d)]
    u = S_["u"]
    Phi = block(u, d)
    SPhi = S_func(Phi)
    check("U6 s%d: S Phi = Phi exactly" % d, same(SPhi, Phi), same(SPhi, fscale(Phi, I)), "compare with i Phi")
    xi = S_["xi"]
    Sxi = S_func(xi)
    check("U6 s%d: S xi = xi exactly at every level %s" % (d, sorted(2 * J for J in xi)), same(Sxi, xi),
          same(Sxi, fscale(xi, -1)), "compare with -xi")
    # block fibres: S acts on block fibres as w -> conj(D^3(r)) Theta w (right factor P fixed)
    Sb = lambda w: conjD_r(3, theta(w, 3))
    Gf = S_["Gf"]
    check("U6 s%d: forcing fibre Pi_6 DN_Phi[xi] is S-even" % d, all((a - b).is_zero() for a, b in zip(Sb(Gf), Gf)),
          all((a + b).is_zero() for a, b in zip(Sb(Gf), Gf)), "compare with minus")
    tx, ty = S_["tang"]["tau_x"], S_["tang"]["tau_y"]
    check("U6 s%d: S tau_x = tau_x and S tau_y = -tau_y exactly" % d,
          all((a - b).is_zero() for a, b in zip(Sb(tx), tx)) and all((a + b).is_zero() for a, b in zip(Sb(ty), ty)),
          all((a - b).is_zero() for a, b in zip(Sb(ty), ty)), "claim S tau_y = +tau_y")
    # L commutes with S (real 14x14 form): S_R L = L S_R
    def real(v):
        return [a.re() for a in v] + [a.im() for a in v]

    def cplx(c):
        return [c[k] + I * c[k + 7] for k in range(7)]
    SR = [real(Sb(cplx([ONE if i == k else ZERO for i in range(14)]))) for k in range(14)]   # columns
    Lm = S_["Lm"]
    SRm = [[SR[k][i] for k in range(14)] for i in range(14)]
    mul = lambda A, B: [[sum((A[i][t] * B[t][j] for t in range(14)), ZERO) for j in range(14)] for i in range(14)]
    LS, SL = mul(Lm, SRm), mul(SRm, Lm)
    ok = all((LS[i][j] - SL[i][j]).is_zero() for i in range(14) for j in range(14))
    Lmut = [row[:] for row in Lm]
    Lmut[0][10] = Lmut[0][10] + 1
    Lmut[10][0] = Lmut[10][0] + 1
    LSm, SLm = mul(Lmut, SRm), mul(SRm, Lmut)
    okm = all((LSm[i][j] - SLm[i][j]).is_zero() for i in range(14) for j in range(14))
    check("U6 s%d: L commutes with S exactly (so B(tau_x, tau_y) = 0)" % d, ok, okm, "L with a symmetric real/imag coupling added")
    kap = S_["kappa"]
    check("U6 s%d: kappa is S-even and its tau_y component is exactly 0" % d,
          all((a - b).is_zero() for a, b in zip(Sb(kap), kap)) and (sum((a.conj() * b for a, b in zip(ty, kap)), ZERO).re()).is_zero(),
          all((a - b).is_zero() for a, b in zip(Sb([c * I for c in kap]), [c * I for c in kap])), "i kappa")
    # the tau_y diagonal entry of L_T is nonzero: needed for 'hence v_y = 0'
    qyy = (sum((a.conj() * b for a, b in zip(ty, cplx([sum((Lm[i][k] * real(ty)[k] for k in range(14)), ZERO) for i in range(14)]))), ZERO).re()) * F(d, 7)
    check("U6 s%d: L_T(tau_y, tau_y) = %s != 0 (needed for v_y = 0)" % (d, qyy.pretty()), not qyy.is_zero(), None, "(value)")

# ------------------------------------------------------------------------------------------ bridge
print("=== bridge: Q'' = 4(<E, DN E>_R - Q) exactly")


def int4(v, d):
    f = block(v, d)
    return inner(f, compress(Nsec(f, levels={3})))


def Qpp(u, e, d):
    """d^2/ds^2 int|cos s Phi_u + sin s Phi_e|^4 at s = 0 from the quartic F(c, s) = sum a_k c^(4-k) s^k."""
    pts = [(1, 0), (0, 1), (1, 1), (1, -1), (1, 2)]
    vals = [int4([c * a + s * b for a, b in zip(u, e)], d) for c, s in pts]
    # solve the 5x5 rational system for a_0..a_4 (exact)
    A = [[F(c) ** (4 - k) * F(s) ** k for k in range(5)] for c, s in pts]
    # Gaussian elimination over MQ right-hand sides
    M = [[mq(x) for x in row] + [v] for row, v in zip(A, vals)]
    for col in range(5):
        piv = next(i for i in range(col, 5) if not M[i][col].is_zero())
        M[col], M[piv] = M[piv], M[col]
        inv = M[col][col].inv()
        M[col] = [x * inv for x in M[col]]
        for i in range(5):
            if i != col and not M[i][col].is_zero():
                fct = M[i][col]
                M[i] = [x - fct * y for x, y in zip(M[i], M[col])]
    a = [M[k][5] for k in range(5)]
    return a[2] * 2 - a[0] * 4


def qform(S_, e1, e2, d):
    Phi = block(S_["u"], d)
    DF = block_fibre(compress(DNsec(Phi, block(e2, d), levels={3})), d)
    return (sum((a.conj() * b for a, b in zip(e1, DF)), ZERO).re() - S_["Q"] * sum((a.conj() * b for a, b in zip(e1, e2)), ZERO).re()) * F(d, 7)


bridge = {}
for name, tangs in (("U5", ("e_t", "i_e_t")), ("U6", ("tau_x", "tau_y"))):
    for d in (3, 4):
        S_ = st[(name, d)]
        u = S_["u"]
        for tn in tangs:
            e = S_["tang"][tn]
            q2 = Qpp(u, e, d)
            fm = qform(S_, e, e, d)
            bridge["%s s%d %s" % (name, d, tn)] = {"Qpp": q2.pretty(), "4*form": (fm * 4).pretty()}
            check("%s s%d %s: Q'' = %s = 4 x form" % (name, d, tn, q2.pretty()), (q2 - fm * 4).is_zero(),
                  (q2 - (fm + F(1, 7)) * 4).is_zero() or (q2 - (fm + S_["Q"]) * 4).is_zero(), "4 x (form + 1/7), and 4 x (form + Q)")
        # cross term by polarization: 4 B(E1, E2) = [Q''(E+) - Q''(E-)]/2
        e1, e2 = S_["tang"][tangs[0]], S_["tang"][tangs[1]]
        h = MQ.sqrt(F(1, 2))
        ep = [(a + b) * h for a, b in zip(e1, e2)]
        em = [(a - b) * h for a, b in zip(e1, e2)]
        cross = (Qpp(u, ep, d) - Qpp(u, em, d)) * F(1, 8)
        fc = qform(S_, e1, e2, d)
        bridge["%s s%d cross" % (name, d)] = {"Qpp_polarized/4": cross.pretty(), "form": fc.pretty()}
        check("%s s%d cross (%s, %s): polarized Q''/4 = %s = form" % (name, d, tangs[0], tangs[1], cross.pretty()),
              (cross - fc).is_zero(), (cross - fc - 1).is_zero(), "form + 1")

# ------------------------------------------------------------------------------------------ lambda_4 lemma
print("=== lambda_4 lemma")
random.seed(3)
main = json.load(open("out/main_all.json"))
item1 = json.load(open("out/item01.json"))
lam = {}
for name in ("U1", "U2", "U3", "U4", "U5", "U6"):
    for d in (3, 4):
        S_ = st[(name, d)]
        u = S_["u"]
        Phi = block(u, d)
        Nf = S_["Nf"]
        # random block v, complex-orthogonal to Phi
        v = [mq(F(random.randint(-5, 5), 7)) + I * F(random.randint(-5, 5), 3) for _ in range(7)]
        # v deliberately NOT orthogonal to Phi: for v perp Phi, X = Q<Phi, v> = 0 by criticality
        fv = block(v, d)
        lhs = inner(Phi, compress(DNsec(Phi, fv, levels={3})))
        X = inner(Nf, fv)
        ok = (lhs - (X * 2 + X.conj())).is_zero()
        okm = (lhs - (X * 3)).is_zero()
        check("%s s%d: <Phi, DN_Phi v> = 2X + conj X for a random block v (X = %s)" % (name, d, X.pretty()[:30]),
              ok, okm, "3X")
        kf = block(S_["kappa"], d)
        Xk = inner(Nf, kf)
        check("%s s%d: X = <N(Phi), kappa> = 0 exactly" % (name, d), Xk.is_zero(), X.is_zero(), "X at the random v")
        levels = {int(n): mq(F(val)) if "sqrt" not in val else None for n, val in main["%s_s%d" % (name, d)]["levels_N"].items()}
        Q = S_["Q"]
        tot_off = sum((vv for n, vv in levels.items() if n != 6), ZERO)
        lam["%s_s%d" % (name, d)] = {"Q": Q.pretty(), "sum_{n!=6}||Pi_n N||^2": tot_off.pretty(),
                                     "int|Phi|^6 - Q^2": (sum(levels.values(), ZERO) - Q * Q).pretty()}
        check("%s s%d: sum_{n!=6} ||Pi_n N||^2 = int|Phi|^6 - Q^2 = %s > 0, Q = %s != 1" % (name, d, tot_off.pretty(), Q.pretty()),
              (tot_off - (sum(levels.values(), ZERO) - Q * Q)).is_zero() and tot_off.sign() > 0 and not (Q - 1).is_zero(),
              (tot_off - sum(levels.values(), ZERO)).is_zero(), "forget -Q^2")
A1 = json.load(open("audit_results.json"))["item1"]["dim_Hom_sigma_V_(n/2)"]
for s, near, gap in (("3", 10, 72), ("4", 8, 32)):
    lv = [int(n) for n, v in A1["sector" + s].items() if v > 0]
    below = [n for n in lv if n < 6]
    nearest = min(n for n in lv if n > 6)
    check("sector %s (T.a, lambda_4 sign): levels %s; none below 6; nearest other level %d, gap %d" % (s, lv, nearest, nearest * (nearest + 2) - 48),
          below == [] and nearest == near and nearest * (nearest + 2) - 48 == gap, None, "(read from my item 1)")

# ------------------------------------------------------------------------------------------ E1/E2/T.b
print("=== E1/E2/T.b")
# U5: the R_z orbit direction projected perpendicular to Phi is along i e_t
Jz = [[mq(m) if i == k else ZERO for k in range(7)] for i, m in enumerate(range(-3, 4))]
S_ = st[("U5", 3)]
u = S_["uhat"]
X = [sum((Jz[i][k].conj() * u[k] for k in range(7)), ZERO) * I for i in range(7)]
nu = sum((a.conj() * a for a in u), ZERO)
Xp = [a - sum((b.conj() * c for b, c in zip(u, X)), ZERO) / nu * b for a, b in zip(X, u)]
et = [c * MQ.sqrt(F(3, 7)) for c in S_["tang"]["e_t"]]      # unit fibre (tang is scaled by sqrt(7/d))
c = prop(Xp, [I * a for a in et])
cm = prop(Xp, et)
check("U5: P_perp(i conj(J_z) u) = %s * (i e_t) exactly, real multiple (orbit direction)" % (c.pretty() if c else None),
      c is not None and c.is_real(), cm is not None and cm.is_real(), "real multiple of e_t")
# U6: no element of su(2) maps K_H = span{v0, v3 + v-3} into itself
Jp = [[ZERO] * 7 for _ in range(7)]
for m in range(-3, 3):
    Jp[m + 4][m + 3] = MQ.sqrt((3 - m) * (3 + m + 1))
Jm = [[Jp[k][i] for k in range(7)] for i in range(7)]
Jx = [[(a + b) * half for a, b in zip(r1, r2)] for r1, r2 in zip(Jp, Jm)]
Jy = [[(a - b) * half * (-I) for a, b in zip(r1, r2)] for r1, r2 in zip(Jp, Jm)]
KH = [[ZERO] * 3 + [ONE] + [ZERO] * 3, [ONE] + [ZERO] * 5 + [ONE]]
def outside(v):
    # component outside K_H: coordinates m = +-1, +-2, and the antisymmetric part of (v3, v-3)
    return [v[1], v[2], v[4], v[5], v[6] - v[0]]
rows = []
for Jn in (Jx, Jy, Jz):
    col = []
    for w in KH:
        img = [sum((Jn[i][k].conj() * w[k] for k in range(7)), ZERO) * I for i in range(7)]
        col += outside(img)
    rows.append(col)
# real rank of the 3 columns (as real vectors): exact via 3x3 Gram determinant
realv = [[x for c in col for x in (c.re(), c.im())] for col in rows]
Gm = [[sum((a * b for a, b in zip(realv[i], realv[j])), ZERO) for j in range(3)] for i in range(3)]
det = (Gm[0][0] * (Gm[1][1] * Gm[2][2] - Gm[1][2] * Gm[2][1]) - Gm[0][1] * (Gm[1][0] * Gm[2][2] - Gm[1][2] * Gm[2][0])
       + Gm[0][2] * (Gm[1][0] * Gm[2][1] - Gm[1][1] * Gm[2][0]))
check("U6: no nonzero X in su(2) preserves K_H (Gram determinant of the leakage map = %s != 0)" % det.pretty(),
      not det.is_zero(), None, "(value)")
# Hess_T Q invertible on the slice (4 x item-6 forms), from the bridge values
for d in (3, 4):
    fe = qform(st[("U5", d)], st[("U5", d)]["tang"]["e_t"], st[("U5", d)]["tang"]["e_t"], d)
    S6 = st[("U6", d)]
    a = qform(S6, S6["tang"]["tau_x"], S6["tang"]["tau_x"], d)
    b = qform(S6, S6["tang"]["tau_y"], S6["tang"]["tau_y"], d)
    cc = qform(S6, S6["tang"]["tau_x"], S6["tang"]["tau_y"], d)
    check("s%d: Hess_T Q invertible on the slices: U5 %s; U6 det = %s" % (d, (fe * 4).pretty(), ((a * b - cc * cc) * 16).pretty()),
          not fe.is_zero() and not (a * b - cc * cc).is_zero(), None, "(values)")
# T.b: xi in X_H exactly for the finite generators with exact matrices (U3, U4, U6) and the torus weights (U1, U2)
def act_term_ok(J, terms, g, chi):
    D = D_exact(J, g)
    for x, Y in terms:
        pass
    C = coeff({J: terms})[J]
    # left action on the x index: C'[(m,k,a)] = sum_m' conj(D[m][m']) C[(m',k,a)]
    Cn = {}
    for (mm, k, a), val in C.items():
        for m in range(2 * J + 1):
            dd = D[m][mm].conj()
            if dd.is_zero():
                continue
            key = (m, k, a)
            Cn[key] = Cn.get(key, ZERO) + dd * val
    keys = set(C) | set(Cn)
    return all((Cn.get(k, ZERO) - chi * C.get(k, ZERO)).is_zero() for k in keys)
h2 = MQ.sqrt(F(1, 2))
gens = {
    "U3": [([[h2 - I * h2, ZERO], [ZERO, h2 + I * h2]], mq(-1)),                                   # Rz(pi/2)
           ([[h2, -I * h2 * (1 - I) * h2], [-I * h2 * (1 + I) * h2, h2]], mq(-1))],             # R((1,1,0)/sqrt2, pi/2)
    "U4": [([[s3 * half - I * half, ZERO], [ZERO, s3 * half + I * half]], mq(-1)),               # Rz(pi/3)
           ([[ZERO, -I], [-I, ZERO]], mq(-1))],                                                   # Rx(pi)
    "U6": [([[mq(half) - I * s3 * half, ZERO], [ZERO, mq(half) + I * s3 * half]], ONE),           # Rz(2pi/3)
           ([[ZERO, -I], [-I, ZERO]], mq(-1))],
    "U2": [([[ZERO, -I], [-I, ZERO]], mq(-1))],
}
for name, gl in gens.items():
    for d in (3, 4):
        xi = st[(name, d)]["xi"]
        ok = all(act_term_ok(J, t, g, chi) for g, chi in gl for J, t in xi.items())
        okm = all(act_term_ok(J, t, g, -chi) for g, chi in gl for J, t in xi.items())
        check("%s s%d: xi transforms by chi under the exact finite generators, all levels" % (name, d), ok, okm, "-chi")
for name, wt in (("U1", 3), ("U2", 0)):
    for d in (3, 4):
        C = coeff(st[(name, d)]["xi"])
        ok = all(m - J == wt for J, CJ in C.items() for (m, k, a) in CJ)
        okm = all(m - J == wt + 1 for J, CJ in C.items() for (m, k, a) in CJ)
        check("%s s%d: xi lives in torus weight %d (character of the torus) at all levels" % (name, d, wt), ok, okm, "weight + 1")

json.dump({"log": LOG, "bridge": bridge, "lambda4": lam}, open("out2/theorem_checks.json", "w"), indent=1, default=str)
print("total", len(LOG), "PASS", sum(l["ok"] for l in LOG), "mutants run", sum(l["mutant_caught"] is not None for l in LOG), "caught", sum(l["mutant_caught"] is True for l in LOG))
