"""s2b_04_step4_reconcile.py -- 3d, 3e, Step 4, and the reconcile with my stage-1 answer.

All objects rebuilt from the room's worklist conventions.
"""

from pathlib import Path
import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))
IDX = {m: i for i, m in enumerate(ORDER)}
S6 = sp.sqrt(6)
K = 15 / S6


def spin3():
    Jz = sp.diag(*[sp.Integer(m) for m in ORDER])
    Jp = sp.zeros(7, 7)
    Jm = sp.zeros(7, 7)
    for m in ORDER:
        if m + 1 <= 3:
            Jp[IDX[m + 1], IDX[m]] = sp.sqrt(12 - m * (m + 1))
        if m - 1 >= -3:
            Jm[IDX[m - 1], IDX[m]] = sp.sqrt(12 - m * (m - 1))
    return (Jp + Jm) / 2, (Jp - Jm) / (2 * sp.I), Jz


FX, FY, FZ = spin3()
F = [FX, FY, FZ]
Fn = [np.array(sp.matrix2numpy(X, dtype=complex)) for X in F]


def cg_str(m1, m2):
    M = m1 + m2
    if abs(M) > 6 or abs(m1) > 3 or abs(m2) > 3:
        return sp.Integer(0)
    f = sp.factorial
    return sp.sqrt(sp.Rational(f(6) * f(6) * f(6 + M) * f(6 - M),
                               f(12) * f(3 + m1) * f(3 - m1) * f(3 + m2) * f(3 - m2)))


def rho6_sq_vec(cv):
    """cv: dict m -> value"""
    th = {m: (-1) ** (3 - m) * sp.conjugate(cv[-m]) for m in range(-3, 4)}
    tot = 0
    for Q in range(-6, 7):
        s = sp.expand(sum(cg_str(m1, Q - m1) * cv[m1] * th[Q - m1]
                          for m1 in range(-3, 4) if abs(Q - m1) <= 3))
        tot += sp.expand(sp.re(s) ** 2 + sp.im(s) ** 2)
    return sp.simplify(sp.expand(tot))


def main():
    print("=========== 3e: the top eigenspace of (3 (n.f)^2 - 12)/sqrt6 ===========")
    Op = sp.simplify((3 * FZ * FZ - 12 * sp.eye(7)) / S6)
    print("  n = z:  operator is diagonal with entries (3m^2-12)/sqrt6 :")
    print("   ", [sp.radsimp(sp.simplify(Op[i, i])) for i in range(7)], " for m =", ORDER)
    ev = sp.simplify(sp.Matrix(Op).eigenvals())
    print("  eigenvalues with multiplicity:", {sp.radsimp(k): v for k, v in ev.items()})
    print("  top eigenvalue = 15/sqrt6 :", sp.simplify(sp.radsimp(Op[IDX[3], IDX[3]] - K)) == 0)
    print("  its eigenspace is exactly span{|3,3>, |3,-3>} (m = +-3 are the only m with 3m^2-12 = 15):",
          [m for m in ORDER if sp.simplify(3 * m ** 2 - 12 - 15) == 0])
    print("  so the top eigenspace is 2-dimensional. (The note asserts this with no computation.)")

    # A_E for E = (3 n n^T - I)/sqrt6 equals (3 (n.f)^2 - 12)/sqrt6
    rng = np.random.default_rng(11)
    worst = 0.0
    for _ in range(40):
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        E = (3 * np.outer(n, n) - np.eye(3)) / np.sqrt(6)
        AE = sum(E[i, k] * (Fn[i] @ Fn[k] + Fn[k] @ Fn[i]) / 2 for i in range(3) for k in range(3))
        nf = sum(n[i] * Fn[i] for i in range(3))
        worst = max(worst, np.max(np.abs(AE - (3 * nf @ nf - 12 * np.eye(7)) / np.sqrt(6))))
        # also: is E unit and traceless?
    print("  A_E = (3(n.f)^2 - 12)/sqrt6 for E = (3nn^T - I)/sqrt6 : max error %.2e over 40 n" % worst)
    Ez = (3 * np.outer([0, 0, 1.], [0, 0, 1.]) - np.eye(3)) / np.sqrt(6)
    print("  that E is traceless (%.1e) and unit (%.15f)" % (abs(np.trace(Ez)), np.linalg.norm(Ez)))
    print("  its diagonal is (2,-1,-1)/sqrt6 :", np.round(np.diag(Ez) * np.sqrt(6), 12))

    print("\n=========== 3d: from ||Q|| = 15/sqrt6 to the top eigenspace ===========")
    hexv = np.zeros(7, dtype=complex)
    hexv[IDX[3]] = hexv[IDX[-3]] = 1 / np.sqrt(2)
    Nb = np.array([[np.real(np.vdot(hexv, Fn[i] @ Fn[k] @ hexv)) for k in range(3)] for i in range(3)])
    Q = Nb - 4 * np.eye(3)
    print("  hexagon: Nbar = "); print(np.round(Nb, 12))
    print("  ||Q|| = %.15f    15/sqrt6 = %.15f" % (np.linalg.norm(Q), 15 / np.sqrt(6)))
    Est = Q / np.linalg.norm(Q)
    print("  E* = Q/||Q|| has eigenvalues %s ; times sqrt6 = %s"
          % (np.round(np.linalg.eigvalsh(Est), 12), np.round(np.linalg.eigvalsh(Est) * np.sqrt(6), 9)))
    AEs = sum(Est[i, k] * (Fn[i] @ Fn[k] + Fn[k] @ Fn[i]) / 2 for i in range(3) for k in range(3))
    w, V = np.linalg.eigh(AEs)
    print("  lambda_max(A_{E*}) = %.15f ; <u|A_{E*}|u> = %.15f"
          % (w.max(), np.real(np.vdot(hexv, AEs @ hexv))))
    proj = sum(abs(np.vdot(V[:, k], hexv)) ** 2 for k in range(7) if abs(w[k] - w.max()) < 1e-10)
    print("  weight of u in the top eigenspace of A_{E*} : %.15f (1 = fully inside)" % proj)
    print("  SUPPLIED argument: tr(rho A_{E*}) = ||Q|| = 15/sqrt6 and lambda_max(A_E) <= 15/sqrt6")
    print("  for every unit traceless E, so tr(rho A_{E*}) = lambda_max(A_{E*}); a pure state whose")
    print("  expectation equals the top eigenvalue lies in the top eigenspace, and E* must itself")
    print("  be an equality point, i.e. E* = (3nn^T - I)/sqrt6.  The note states the conclusion only.")

    print("\n=========== Step 4 ===========")
    al, be = sp.symbols("alpha beta")
    ra, rb, p1, p2 = sp.symbols("ra rb phi1 phi2", real=True, positive=True)
    cv = {m: sp.Integer(0) for m in range(-3, 4)}
    cv[3] = ra * sp.exp(sp.I * p1)
    cv[-3] = rb * sp.exp(sp.I * p2)
    uvec = sp.Matrix(7, 1, [cv[m] for m in ORDER])
    fexp = [sp.simplify(sp.expand((uvec.H * X * uvec)[0, 0])) for X in F]
    print("  <f> for u = alpha|3,3> + beta|3,-3> :", [sp.simplify(x) for x in fexp])
    print("    -> <f> = 3(|alpha|^2 - |beta|^2) z, as the note says :",
          sp.simplify(fexp[0]) == 0 and sp.simplify(fexp[1]) == 0
          and sp.simplify(fexp[2] - 3 * (ra ** 2 - rb ** 2)) == 0)
    # |a00|^2 when |alpha| = |beta|
    cv2 = dict(cv)
    cv2[3] = sp.exp(sp.I * p1) / sp.sqrt(2)
    cv2[-3] = sp.exp(sp.I * p2) / sp.sqrt(2)
    th = {m: (-1) ** (3 - m) * sp.conjugate(cv2[-m]) for m in range(-3, 4)}
    s = sp.simplify(sp.expand(sum(sp.conjugate(th[m]) * cv2[m] for m in range(-3, 4))))
    print("  |alpha|=|beta|=1/sqrt2 : <Theta u, u> =", s, " so |a00|^2 =",
          sp.simplify(sp.Abs(s) ** 2 / 7), "= 1/7 :", sp.simplify(sp.Abs(s) ** 2 / 7 - sp.Rational(1, 7)) == 0)
    print("    (the note says '|a00|^2 = 1/7 follows' with no computation)")
    # r6 on the family
    r6fam = sp.simplify(rho6_sq_vec(cv2))
    print("  ||rho6||^2 on the whole |alpha|=|beta| family (any phases) :", sp.nsimplify(r6fam),
          " = 463/924 :", sp.simplify(r6fam - sp.Rational(463, 924)) == 0)
    # unequal moduli
    print("  unequal moduli: r6 at |alpha|^2 = t")
    t = sp.Symbol("t", positive=True)
    cv3 = {m: sp.Integer(0) for m in range(-3, 4)}
    cv3[3] = sp.sqrt(t)
    cv3[-3] = sp.sqrt(1 - t)
    r6t = sp.simplify(sp.expand(rho6_sq_vec(cv3)))
    print("    r6(t) =", sp.nsimplify(sp.simplify(r6t)))
    print("    at t=1/2 :", sp.nsimplify(sp.simplify(r6t.subs(t, sp.Rational(1, 2)))),
          " at t=3/4 :", sp.nsimplify(sp.simplify(r6t.subs(t, sp.Rational(3, 4)))),
          " at t=1 :", sp.nsimplify(sp.simplify(r6t.subs(t, 1))))
    # relative phase is a rotation about n
    chi = sp.Symbol("chi", real=True)
    print("  rotation by chi about z sends |3,m> -> e^{-i m chi}|3,m>, so the relative phase")
    print("    phi1 - phi2 shifts by -6 chi; chi in [0, pi/3) covers the full circle :",
          sp.simplify(sp.exp(-sp.I * 6 * (2 * sp.pi / 6)) - 1) == 0)
    # the arithmetic
    print("  -5/231 + 1/11 + 171/396 =", sp.Rational(-5, 231) + sp.Rational(1, 11) + sp.Rational(171, 396),
          " = 463/924 :", sp.Rational(-5, 231) + sp.Rational(1, 11) + sp.Rational(171, 396) == sp.Rational(463, 924))
    print("  48 + 225/6 =", 48 + sp.Rational(225, 6), " = 171/2 :", 48 + sp.Rational(225, 6) == sp.Rational(171, 2))
    print("  (15/sqrt6)^2 =", sp.nsimplify(sp.simplify(K ** 2)))
    print("  signs: coefficient of |f|^2 is -1/22 < 0 (so the bound uses |f|^2 >= 0),")
    print("         of |a00|^2 is +7/11 > 0, of TrN^2 is +1/198 > 0 (so both upper bounds are used).")

    print("\n=========== RECONCILE with my stage-1 answer ===========")
    # my stage-1 constants n_J, recomputed here by an INDEPENDENT route:
    # r6 = sum_J n_J p_J with p_J(v_m) = |<3m;3m|J,2m>|^2, a triangular system on v3,v2,v1,v0
    Js = [0, 2, 4, 6]
    rowsP, rhs = [], []
    for m in (3, 2, 1, 0):
        row = []
        for J in Js:
            c = CG(3, m, 3, m, J, 2 * m).doit() if abs(2 * m) <= J else 0
            row.append(sp.nsimplify(sp.simplify(sp.Abs(c) ** 2)))
        rowsP.append(row)
        cv = {k: (sp.Integer(1) if k == m else sp.Integer(0)) for k in range(-3, 4)}
        rhs.append(rho6_sq_vec(cv))
    P = sp.Matrix(rowsP)
    print("  p_J(v_m) matrix (rows m=3,2,1,0; cols J=0,2,4,6):")
    sp.pprint(P)
    print("  row sums (must be 1):", [sp.simplify(sum(P.row(i))) for i in range(4)])
    nJ = P.solve(sp.Matrix(rhs))
    print("  solved n_J =", [sp.nsimplify(x) for x in nJ])
    print("  my stage-1 values (13/7, 65/84, 13/154, 1/924) :",
          [sp.simplify(nJ[i] - v) == 0 for i, v in enumerate(
              [sp.Rational(13, 7), sp.Rational(65, 84), sp.Rational(13, 154), sp.Rational(1, 924)])])

    # the link between my p_2 and the author's TrNbar^2, on the real form
    print("\n  on the real form {Theta u = u}: is TrNbar^2 = 48 + 126 p_2 and |f|^2 = 0 ?")
    tsym = sp.symbols("t0:7", real=True)
    cr = {0: sp.I * tsym[0]}
    for m, (pp, qq) in zip([1, 2, 3], [(tsym[1], tsym[2]), (tsym[3], tsym[4]), (tsym[5], tsym[6])]):
        cr[m] = pp + sp.I * qq
        cr[-m] = (-1) ** (3 + m) * sp.conjugate(cr[m])
    ur = sp.Matrix(7, 1, [cr[m] for m in ORDER])
    nrm = sp.expand(sum(sp.Abs(cr[m]) ** 2 for m in range(-3, 4)))
    fex = [sp.simplify(sp.expand((ur.H * X * ur)[0, 0])) for X in F]
    print("   |f|^2 on the real form :", sp.simplify(sum(sp.re(x) ** 2 + sp.im(x) ** 2 for x in fex)))
    Nb = sp.Matrix(3, 3, lambda i, k: sp.simplify(sp.re(sp.expand((ur.H * F[i] * F[k] * ur)[0, 0]))))
    trn2 = sp.simplify(sp.expand((Nb * Nb).trace()))
    p2 = 0
    for Qq in range(-2, 3):
        e = sp.expand(sum(sp.nsimplify(sp.simplify(CG(3, m1, 3, Qq - m1, 2, Qq).doit())) * cr[m1] * cr[Qq - m1]
                          for m1 in range(-3, 4) if abs(Qq - m1) <= 3))
        p2 += sp.expand(sp.re(e) ** 2 + sp.im(e) ** 2)
    p2 = sp.simplify(sp.expand(p2))
    print("   TrNbar^2 - (48*||u||^4 + 126*p_2) on the real form :",
          sp.simplify(sp.expand(trn2 - 48 * nrm ** 2 - 126 * p2)))
    print("   -> so the author's 'maximise TrNbar^2' and my stage-1 'maximise the quadrupole")
    print("      weight p_2' are the SAME optimisation, with ||Q||^2 = 126 p_2.")
    print("   check at the hexagon: 126 * 25/84 =", sp.Rational(126 * 25, 84), "= 75/2 = (15/sqrt6)^2 :",
          sp.Rational(126 * 25, 84) == sp.Rational(75, 2))

    print("\n  the author's equality set vs mine:")
    print("   author: u = (e^{i p1}|3,3>_n + e^{i p2}|3,-3>_n)/sqrt2, over n in S^2 and the phases")
    print("   mine  : {e^{i phi} D^3(g) (v3+v-3)/sqrt2}, i.e. Majorana constellation = regular hexagon")
    # show they coincide: rotate the hexagon to an arbitrary axis and compare with the author's family
    def Dm(nax, ang):
        Aop = nax[0] * Fn[0] + nax[1] * Fn[1] + nax[2] * Fn[2]
        w, V = np.linalg.eigh(Aop)
        return V @ np.diag(np.exp(-1j * ang * w)) @ V.conj().T
    rng2 = np.random.default_rng(77)
    err = 0.0
    for _ in range(200):
        nax = rng2.normal(size=3)
        nax /= np.linalg.norm(nax)
        ang = rng2.uniform(0, 2 * np.pi)
        ph = rng2.uniform(0, 2 * np.pi)
        w = np.exp(1j * ph) * (Dm(nax, ang) @ hexv)
        Nb2 = np.array([[np.real(np.vdot(w, Fn[i] @ Fn[k] @ w)) for k in range(3)] for i in range(3)])
        # r6 via the author's Step-1 identity, all pieces computed here
        f2 = sum(abs(np.vdot(w, Fn[i] @ w)) ** 2 for i in range(3)).real
        thv = np.array([(-1.0) ** (3 - m) * np.conjugate(w[IDX[-m]]) for m in ORDER])
        a00 = abs(np.vdot(thv, w)) ** 2 / 7
        r6 = -5 / 231 - f2 / 22 + (7 / 11) * a00 + np.trace(Nb2 @ Nb2) / 198
        err = max(err, abs(r6 - 463 / 924))
    print("   200 random rotations+phases of the hexagon: max |r6 - 463/924| = %.2e" % err)
    print("   every member of the author's family is a rotated, rephased hexagon and conversely,")
    print("   since the relative phase is absorbed by a rotation about n (checked above).")


if __name__ == "__main__":
    main()
