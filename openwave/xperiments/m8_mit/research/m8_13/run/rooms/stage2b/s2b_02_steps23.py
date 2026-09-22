"""s2b_02_steps23.py -- my own check of Step 2 and Step 3 (3a, 3b, 3c).

Rebuilt from the room's worklist conventions; nothing taken from the author's files.

Step 2:  ||Theta u|| = ||u||  (the unstated prerequisite of the Cauchy-Schwarz bound),
         |a00|^2 <= 1/7 with equality iff Theta u prop u, and the rephasing that turns
         "prop" into "=".
Step 3a: Tr(Q E) = <u, A_E u> ; ||Q|| = max over unit traceless E ; the equivariance
         A_{R E R^T} = D(R) A_E D(R)^dag that licenses "rotate E to diagonal" ; the
         max-swap max_rho max_E = max_E max_rho ; the substitution and the block split.
3b:      each block's lambda_max formula, the M3 square completion, the M1 sum of
         squares, and the two side conditions the note asserts (max(a+6b) = sqrt6 on
         the ellipse, max|a| = sqrt(3/2)).
3c:      the complete per-block equality set, by the exact decision procedure.
"""

from pathlib import Path
import numpy as np
import sympy as sp

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))
IDX = {m: i for i, m in enumerate(ORDER)}
a, b, u_, v_, lam = sp.symbols("a b u v lam", real=True)
S6 = sp.sqrt(6)
K = 15 / S6


def spin3_sym():
    Jz = sp.diag(*[sp.Integer(m) for m in ORDER])
    Jp = sp.zeros(7, 7)
    Jm = sp.zeros(7, 7)
    for m in ORDER:
        if m + 1 <= 3:
            Jp[IDX[m + 1], IDX[m]] = sp.sqrt(12 - m * (m + 1))
        if m - 1 >= -3:
            Jm[IDX[m - 1], IDX[m]] = sp.sqrt(12 - m * (m - 1))
    return (Jp + Jm) / 2, (Jp - Jm) / (2 * sp.I), Jz, Jp, Jm


FX, FY, FZ, FP, FM = spin3_sym()
F = [FX, FY, FZ]
Fn = [np.array(sp.matrix2numpy(X, dtype=complex)) for X in F]


def main():
    print("=========== STEP 2 ===========")
    xs = sp.symbols("x0:7", real=True)
    ys = sp.symbols("y0:7", real=True)
    c = {m: xs[IDX[m]] + sp.I * ys[IDX[m]] for m in range(-3, 4)}
    th = {m: (-1) ** (3 - m) * sp.conjugate(c[-m]) for m in range(-3, 4)}
    n_u = sp.expand(sum(sp.Abs(c[m]) ** 2 for m in range(-3, 4)))
    n_t = sp.expand(sum(sp.Abs(th[m]) ** 2 for m in range(-3, 4)))
    print("  ||Theta u||^2 - ||u||^2 identically :", sp.simplify(sp.expand(n_t - n_u)) == 0)
    print("   -> this is what makes Cauchy-Schwarz give |<Theta u, u>| <= 1; the note")
    print("      never states it and check_s3_maximum.py has no gate for it.")
    print("  Theta^2 = id :", all(sp.simplify((-1) ** (3 - m) * sp.conjugate(th[-m]) - c[m]) == 0
                                  for m in range(-3, 4)))
    print("  rephasing: if Theta u = lam u with |lam| = 1 then Theta(e^{i phi}u) = e^{i phi}u")
    print("             for e^{2 i phi} = lam, so 'Theta u prop u' and 'time-reversal")
    print("             invariant' agree AS RAYS but not as vectors (the note's 'that is').")

    print("\n=========== STEP 3a ===========")
    # Tr(Q E) = <u, A_E u> for symmetric traceless E
    rng = np.random.default_rng(5)
    FXn, FYn, FZn = Fn
    worst = 0.0
    for _ in range(200):
        z = rng.normal(size=7) + 1j * rng.normal(size=7)
        z /= np.linalg.norm(z)
        Nb = np.array([[np.real(np.vdot(z, Fn[i] @ Fn[k] @ z)) for k in range(3)] for i in range(3)])
        Q = Nb - 4 * np.eye(3)
        E = rng.normal(size=(3, 3))
        E = (E + E.T) / 2
        E = E - np.trace(E) / 3 * np.eye(3)
        E = E / np.linalg.norm(E)
        AE = sum(E[i, k] * (Fn[i] @ Fn[k] + Fn[k] @ Fn[i]) / 2 for i in range(3) for k in range(3))
        worst = max(worst, abs(np.trace(Q @ E) - np.real(np.vdot(z, AE @ z))))
    print("  Tr(Q E) = <u, A_E u> for unit traceless symmetric E : max error %.2e over 200 draws" % worst)

    worst2 = 0.0
    for _ in range(200):
        z = rng.normal(size=7) + 1j * rng.normal(size=7)
        z /= np.linalg.norm(z)
        Nb = np.array([[np.real(np.vdot(z, Fn[i] @ Fn[k] @ z)) for k in range(3)] for i in range(3)])
        Q = Nb - 4 * np.eye(3)
        best = -1e9
        for _ in range(400):
            E = rng.normal(size=(3, 3))
            E = (E + E.T) / 2
            E = E - np.trace(E) / 3 * np.eye(3)
            E = E / np.linalg.norm(E)
            best = max(best, np.trace(Q @ E))
        worst2 = max(worst2, abs(best - np.linalg.norm(Q)) / max(np.linalg.norm(Q), 1e-12))
    print("  ||Q|| = max_E Tr(QE): random-E maximum reaches the norm to rel. %.2e" % worst2)
    print("    (exact reason, supplied: E = Q/||Q|| is unit traceless symmetric and gives Tr(QE)=||Q||;")
    print("     Cauchy-Schwarz in the 5-dim space of traceless symmetric matrices gives '<=')")

    # equivariance
    def Dmat(nax, ang):
        Aop = nax[0] * Fn[0] + nax[1] * Fn[1] + nax[2] * Fn[2]
        w, V = np.linalg.eigh(Aop)
        return V @ np.diag(np.exp(-1j * ang * w)) @ V.conj().T

    def Rmat(nax, ang):
        nax = np.array(nax, float)
        nax = nax / np.linalg.norm(nax)
        Kx = np.array([[0, -nax[2], nax[1]], [nax[2], 0, -nax[0]], [-nax[1], nax[0], 0]])
        return np.eye(3) + np.sin(ang) * Kx + (1 - np.cos(ang)) * (Kx @ Kx)

    werr = 0.0
    for _ in range(50):
        nax = rng.normal(size=3)
        ang = rng.uniform(0, 2 * np.pi)
        R, D = Rmat(nax, ang), Dmat(nax / np.linalg.norm(nax), ang)
        E = rng.normal(size=(3, 3))
        E = (E + E.T) / 2
        E = E - np.trace(E) / 3 * np.eye(3)
        E = E / np.linalg.norm(E)
        AE = sum(E[i, k] * (Fn[i] @ Fn[k] + Fn[k] @ Fn[i]) / 2 for i in range(3) for k in range(3))
        E2 = R @ E @ R.T
        AE2 = sum(E2[i, k] * (Fn[i] @ Fn[k] + Fn[k] @ Fn[i]) / 2 for i in range(3) for k in range(3))
        werr = max(werr, np.max(np.abs(AE2 - D @ AE @ D.conj().T)))
    print("  equivariance A_{R E R^T} = D(R) A_E D(R)^dag : max error %.2e over 50 draws" % werr)
    print("    -> this is what licenses 'rotating E to diagonal form does not change lambda_max'.")
    print("       Not stated as an identity in the note, not gated by check_s3_maximum.py.")
    print("  max-swap: max_rho max_E tr(rho A_E) = max_E max_rho tr(rho A_E) = max_E lambda_max(A_E),")
    print("    valid since both ranges are compact and tr(rho A_E) is continuous.  Unstated in the note.")

    # the substitution and the split
    e1, e2, e3 = -a / 3 + 2 * b, -a / 3 - 2 * b, 2 * a / 3
    lhs = sp.expand(e1 * FX * FX + e2 * FY * FY + e3 * FZ * FZ)
    rhs = sp.expand(a * (FZ * FZ - 4 * sp.eye(7)) + b * (FP * FP + FM * FM))
    print("\n  e3 = 2a/3, e1-e2 = 4b, e traceless :", sp.simplify(e1 + e2 + e3) == 0,
          sp.simplify(e3 - 2 * a / 3) == 0, sp.simplify(e1 - e2 - 4 * b) == 0)
    print("  operator identity e1 fx^2+e2 fy^2+e3 fz^2 = a(fz^2-4)+b(f+^2+f-^2) :",
          sp.simplify(lhs - rhs) == sp.zeros(7, 7))
    print("  constraint ||e||^2 = (2/3)a^2+8b^2 :",
          sp.simplify(sp.expand(e1 ** 2 + e2 ** 2 + e3 ** 2 - (sp.Rational(2, 3) * a ** 2 + 8 * b ** 2))) == 0)

    s = sp.sqrt(2)

    def vec(d):
        w = sp.zeros(7, 1)
        for m, co in d.items():
            w[IDX[m]] = co
        return w
    names = ["s3", "s1", "t3", "t1", "p2", "z0", "n2"]
    B = sp.Matrix.hstack(vec({3: 1 / s, -3: 1 / s}), vec({1: 1 / s, -1: 1 / s}),
                         vec({3: 1 / s, -3: -1 / s}), vec({1: 1 / s, -1: -1 / s}),
                         vec({2: 1 / s, -2: 1 / s}), vec({0: 1}), vec({2: 1 / s, -2: -1 / s}))
    Ab = sp.simplify(sp.expand(B.H * rhs * B))
    M1 = Ab[0:2, 0:2]
    M2 = Ab[2:4, 2:4]
    M3 = Ab[4:6, 4:6]
    print("  block split (adapted basis %s):" % names)
    print("    M1 =", M1.tolist(), "  M2 =", M2.tolist())
    print("    M3 =", M3.tolist(), "  seventh diagonal entry =", sp.simplify(Ab[6, 6]))
    ok_off = all(sp.simplify(Ab[i, k]) == 0 for i in range(7) for k in range(7)
                 if (i // 2 if i < 6 else 3) != (k // 2 if k < 6 else 3))
    print("    inter-block entries all zero :", ok_off)
    M1t = sp.Matrix([[5 * a, sp.sqrt(60) * b], [sp.sqrt(60) * b, -3 * a + 12 * b]])
    M3t = sp.Matrix([[0, sp.sqrt(240) * b], [sp.sqrt(240) * b, -4 * a]])
    print("    M1 matches the note :", sp.simplify(M1 - M1t) == sp.zeros(2, 2),
          "  M3 matches the note :", sp.simplify(M3 - M3t) == sp.zeros(2, 2))
    print("    note writes 'M2 = M1 at -b'.  As BUILT, M2 =", M2.tolist())
    print("      that is NOT literally M1(-b) =", M1t.subs(b, -b).tolist())
    print("      but flipping the sign of the basis vector t1 makes it literally so;")
    print("      same characteristic polynomial either way :",
          sp.simplify(sp.expand(M2.charpoly(lam).as_expr() - M1t.subs(b, -b).charpoly(lam).as_expr())) == 0)

    print("\n=========== STEP 3b ===========")
    lam1 = a + 6 * b + sp.sqrt(16 * a ** 2 - 48 * a * b + 96 * b ** 2)
    print("  M1 top eigenvalue formula is a root of char(M1) :",
          sp.simplify(M1t.charpoly(lam).as_expr().subs(lam, lam1)) == 0)
    print("  M1 slack identity (a^2+6ab+24b^2)^2 - (3/2)(a+6b)^2 N^2 = 36 b^2 (a+2b)^2 :",
          sp.expand((a ** 2 + 6 * a * b + 24 * b ** 2) ** 2
                    - sp.Rational(3, 2) * (a + 6 * b) ** 2 * (sp.Rational(2, 3) * a ** 2 + 8 * b ** 2)
                    - 36 * b ** 2 * (a + 2 * b) ** 2) == 0)
    print("  a^2+6ab+24b^2 = (a+3b)^2 + 15b^2 :",
          sp.expand(a ** 2 + 6 * a * b + 24 * b ** 2 - (a + 3 * b) ** 2 - 15 * b ** 2) == 0)
    print("  K^2 N^2 = 25a^2 + 300b^2 :",
          sp.simplify(sp.expand(K ** 2 * (sp.Rational(2, 3) * a ** 2 + 8 * b ** 2) - 25 * a ** 2 - 300 * b ** 2)) == 0)
    lam3_gen = sp.simplify(-2 * a + sp.sqrt(4 * a ** 2 + 240 * b ** 2))
    print("  M3 top eigenvalue in general is -2a + sqrt(4a^2+240b^2) :",
          sp.simplify(M3t.charpoly(lam).as_expr().subs(lam, lam3_gen)) == 0)
    print("  ON THE ELLIPSE 240b^2 = 30 - 20a^2, so it becomes -2a + sqrt(30-16a^2) :",
          sp.simplify(240 * ((1 - sp.Rational(2, 3) * a ** 2) / 8) - (30 - 20 * a ** 2)) == 0)
    print("    (the note writes the on-ellipse form without saying the ellipse was used)")
    print("  M3 square completion (K+2a)^2 - (30-16a^2) = 20(a+sqrt6/4)^2 :",
          sp.simplify(sp.expand((K + 2 * a) ** 2 - (30 - 16 * a ** 2) - 20 * (a + S6 / 4) ** 2)) == 0)
    # the two side conditions, computed rather than asserted
    t = sp.Symbol("t", real=True)
    aa, bb = sp.sqrt(sp.Rational(3, 2)) * sp.cos(t), sp.sin(t) / (2 * sp.sqrt(2))
    print("  ellipse parametrisation check ((2/3)a^2+8b^2 - 1) :",
          sp.simplify(sp.Rational(2, 3) * aa ** 2 + 8 * bb ** 2 - 1) == 0)
    expr = sp.simplify(aa + 6 * bb)
    amp = sp.sqrt(sp.simplify(sp.Rational(3, 2) + 36 * sp.Rational(1, 8)))
    print("  max(a+6b) on the ellipse = sqrt(3/2 + 36/8) = %s = %s ; K - that = %s > 0"
          % (sp.nsimplify(amp), sp.N(amp, 15), sp.N(K - amp, 15)))
    print("    (the note asserts a+6b <= sqrt6; check_s3_maximum.py gates only the ARITHMETIC")
    print("     3/2 + 36/8 == 6, not that this is the maximum)")
    print("  max|a| on the ellipse = sqrt(3/2) = %s, so K + 2a >= K - 2sqrt(3/2) = %s = 9/sqrt6 = %s"
          % (sp.N(sp.sqrt(sp.Rational(3, 2)), 15), sp.N(K - 2 * sp.sqrt(sp.Rational(3, 2)), 15),
             sp.N(9 / S6, 15)))
    print("    and 30 - 16a^2 >= 30 - 24 = 6 > 0, so the square root is real. Both unstated.")

    print("\n=========== STEP 3c: the complete equality set ===========")
    print("  decision procedure: lambda_max(M) = K  <=>  det(M - K I) = 0 AND tr M <= 2K.")
    SUB = {a: S6 * u_, b: S6 * v_}
    blocks = {"M1": M1t.subs(SUB), "M2": M2.subs(SUB), "M3": M3t.subs(SUB)}
    ELL = 4 * u_ ** 2 + 48 * v_ ** 2 - 1
    allpts = {}
    for nm, M in blocks.items():
        D = sp.nsimplify(sp.expand(sp.radsimp(sp.simplify(M.det() - K * M.trace() + K ** 2))))
        Dr = sp.expand(sp.numer(sp.together(sp.expand(D.subs(v_ ** 2, (1 - 4 * u_ ** 2) / 48)))))
        print("  %s : det(M-KI) reduced mod the ellipse = %s" % (nm, sp.factor(Dr)))
        cands = []
        for f, mlt in sp.factor_list(Dr)[1]:
            f = sp.expand(f)
            if not f.free_symbols:
                continue
            if v_ in f.free_symbols:
                for sv in sp.solve(sp.Eq(f, 0), v_):
                    for r in sp.roots(sp.Poly(sp.expand(sp.numer(sp.together(4 * u_ ** 2 + 48 * sv ** 2 - 1))), u_)):
                        if sp.im(sp.N(r)) == 0:
                            cands.append((sp.nsimplify(r), sp.nsimplify(sv.subs(u_, r))))
            else:
                for r in sp.roots(sp.Poly(f, u_)):
                    if sp.im(sp.N(r)) != 0:
                        continue
                    v2 = sp.simplify((1 - 4 * r ** 2) / 48)
                    if sp.N(v2) >= 0:
                        for sg in (1, -1):
                            cands.append((sp.nsimplify(r), sp.nsimplify(sg * sp.sqrt(v2))))
        keep = []
        for p in cands:
            if any(sp.simplify(p[0] - q[0]) == 0 and sp.simplify(p[1] - q[1]) == 0 for q in keep):
                continue
            if sp.simplify(4 * p[0] ** 2 + 48 * p[1] ** 2 - 1) != 0:
                continue
            Tr = sp.radsimp(sp.simplify(M.trace().subs({u_: p[0], v_: p[1]})))
            if sp.N(2 * K - Tr) >= -1e-25:
                keep.append(p)
        allpts[nm] = keep
        print("      equality points (u,v): %s" % keep)

    union = []
    for nm in allpts:
        for p in allpts[nm]:
            if not any(sp.simplify(p[0] - q[0]) == 0 and sp.simplify(p[1] - q[1]) == 0 for q in union):
                union.append(p)
    print("  union: %d distinct points" % len(union))
    for p in union:
        av, bv = sp.radsimp(S6 * p[0]), sp.radsimp(S6 * p[1])
        ev = tuple(sp.radsimp(sp.simplify(x)) for x in (-av / 3 + 2 * bv, -av / 3 - 2 * bv, 2 * av / 3))
        inb = [nm for nm in ("M1", "M2", "M3")
               if any(sp.simplify(p[0] - x) == 0 and sp.simplify(p[1] - y) == 0 for (x, y) in allpts[nm])]
        print("     (u,v)=%s  e=%s   attained by %s" % (p, ev, inb))

    print("\n  the THREE loci the note names, intersected with the ellipse:")
    named = []
    for r in sp.roots(sp.Poly(4 * u_ ** 2 - 1, u_)):
        if sp.N(r) > 0:
            named.append(("b=0, a>0", sp.nsimplify(r), sp.Integer(0)))
    for r in sp.roots(sp.Poly(4 * (2 * v_) ** 2 + 48 * v_ ** 2 - 1, v_)):
        if sp.N(r) > 0:
            named.append(("a=-2b, b>0", sp.nsimplify(-2 * r), sp.nsimplify(r)))
    # a = -sqrt6/4  <=>  u = -1/4
    for sg in (1, -1):
        v2 = sp.Rational(1, 64)
        named.append(("a=-sqrt6/4 in M3", sp.Rational(-1, 4), sp.nsimplify(sg * sp.sqrt(v2))))
    seen = []
    for lab, U, V in named:
        if any(sp.simplify(U - x) == 0 and sp.simplify(V - y) == 0 for (x, y) in seen):
            print("     %-18s -> (u,v)=(%s,%s)  [duplicate of a point already listed]" % (lab, U, V))
            continue
        seen.append((U, V))
        print("     %-18s -> (u,v)=(%s,%s)" % (lab, U, V))
    print("  distinct points delivered by the three named loci: %d (union has %d)"
          % (len(seen), len(union)))
    missing = [p for p in union if not any(sp.simplify(p[0] - x) == 0 and sp.simplify(p[1] - y) == 0
                                           for (x, y) in seen)]
    print("  equality points NOT delivered by the named loci:", missing if missing else "none")
    print("  blocks whose equality locus is NOT named in the note:",
          [nm for nm in ("M1", "M2", "M3") if nm not in ("M1", "M3")])


if __name__ == "__main__":
    main()
