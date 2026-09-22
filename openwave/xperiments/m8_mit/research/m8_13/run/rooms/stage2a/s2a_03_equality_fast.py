"""s2a_03_equality_fast.py -- the equality set of step 3c, exhaustively.

(s2a_02_equality_set.py and s2a_02b_equality_set.py are kept but were KILLED: sympy's
`solve` on the radical / two-variable systems ran past the budget.  Here the
elimination is done explicitly and only univariate factorisation is asked of sympy,
which is instant.)

Setup, all rebuilt in s2a_01_blocks.py and re-derived here:
    a = sqrt(6) u ,  b = sqrt(6) v      ellipse (2/3)a^2+8b^2=1  <->  4u^2+48v^2=1
    L0 = 15/sqrt(6) = (5/2) sqrt(6)

For a 2x2 real symmetric M the statement  lambda_max(M) = L0  is EXACTLY equivalent to
        det(M - L0 I) = 0    AND    tr M <= 2 L0
(the first says L0 is an eigenvalue, the second says it is the larger one).  So the
equality set of a block is the intersection of a CONIC with the ELLIPSE, a
zero-dimensional system: reduce the conic modulo the ellipse (which eliminates v^2),
factor, and read off every root.  Nothing is found by inspection.
"""

from pathlib import Path
import numpy as np
import sympy as sp

HERE = Path(__file__).parent
u, v, L = sp.symbols("u v L", real=True)
S6 = sp.sqrt(6)
L0 = 15 / S6
ORDER = list(range(3, -4, -1))
IDX = {m: i for i, m in enumerate(ORDER)}

M1 = sp.Matrix([[5 * S6 * u, sp.sqrt(60) * S6 * v],
                [sp.sqrt(60) * S6 * v, -3 * S6 * u + 12 * S6 * v]])
M2 = sp.Matrix([[5 * S6 * u, sp.sqrt(60) * S6 * v],
                [sp.sqrt(60) * S6 * v, -3 * S6 * u - 12 * S6 * v]])
M3 = sp.Matrix([[0, sp.sqrt(240) * S6 * v], [sp.sqrt(240) * S6 * v, -4 * S6 * u]])
BLOCKS = {"M1": M1, "M2": M2, "M3": M3}


def e_of(uu, vv):
    a_, b_ = S6 * uu, S6 * vv
    return tuple(sp.radsimp(sp.simplify(x)) for x in (-a_ / 3 + 2 * b_, -a_ / 3 - 2 * b_, 2 * a_ / 3))


def lam_max_num(M, uu, vv):
    Mn = np.array(sp.matrix2numpy(M.subs({u: uu, v: vv}), dtype=float))
    return float(np.linalg.eigvalsh(Mn).max())


def main():
    print("L0 = 15/sqrt(6) = (5/2)sqrt(6) =", sp.N(L0, 25))
    print("The seventh eigenvalue of the operator is identically 0.")
    print("  L0 - 0 =", sp.N(L0, 20), "> 0, so the zero eigenvalue is never an equality point.\n")

    allpts = {}
    for name, M in BLOCKS.items():
        print("=== %s ===" % name)
        D = sp.expand(sp.radsimp(sp.simplify(M.det() - L0 * M.trace() + L0**2)))  # det(M - L0 I)
        D = sp.nsimplify(sp.expand(D))
        print("  det(M - L0 I) =", D)
        # reduce modulo the ellipse: 48 v^2 = 1 - 4 u^2
        Dr = sp.expand(D.subs(v**2, (1 - 4 * u**2) / 48))
        Dr = sp.nsimplify(sp.together(Dr))
        num = sp.expand(sp.numer(sp.together(Dr)))
        print("  reduced mod the ellipse (v^2 -> (1-4u^2)/48), numerator:")
        print("     ", sp.factor(num))
        facs = sp.factor_list(num)
        cands = []
        for f, mult in facs[1]:
            f = sp.expand(f)
            if f.free_symbols == set():
                continue
            print("     factor:", f, " (multiplicity %d)" % mult)
            # intersect factor = 0 with the ellipse
            if v in f.free_symbols:
                sol_v = sp.solve(sp.Eq(f, 0), v)
                for sv in sol_v:
                    poly = sp.Poly(sp.expand(sp.together(4 * u**2 + 48 * sv**2 - 1) * 1), u)
                    poly = sp.Poly(sp.expand(sp.numer(sp.together(4 * u**2 + 48 * sv**2 - 1))), u)
                    for r, m in sp.roots(poly).items():
                        if sp.im(sp.N(r)) != 0:
                            continue
                        cands.append((sp.nsimplify(r), sp.nsimplify(sv.subs(u, r))))
            else:
                for r, m in sp.roots(sp.Poly(f, u)).items():
                    if sp.im(sp.N(r)) != 0:
                        continue
                    vv2 = sp.simplify((1 - 4 * r**2) / 48)
                    for sgn in (1, -1):
                        if sp.N(vv2) < 0:
                            continue
                        cands.append((sp.nsimplify(r), sp.nsimplify(sgn * sp.sqrt(vv2))))
        # dedupe
        uniq = []
        for p in cands:
            if not any(sp.simplify(p[0] - q[0]) == 0 and sp.simplify(p[1] - q[1]) == 0 for q in uniq):
                uniq.append(p)
        keep = []
        for (uu, vv) in uniq:
            onell = sp.simplify(4 * uu**2 + 48 * vv**2 - 1) == 0
            Tr = sp.radsimp(sp.simplify(M.trace().subs({u: uu, v: vv})))
            larger = sp.N(2 * L0 - Tr) >= -1e-25
            lm = lam_max_num(M, float(uu), float(vv))
            print("     candidate (u,v)=(%s,%s): on ellipse %s, tr M = %s, tr <= 2L0 : %s, lam_max = %.15f"
                  % (uu, vv, onell, Tr, larger, lm))
            if onell and larger:
                keep.append((uu, vv))
        allpts[name] = keep
        print("  EQUALITY POINTS of %s : %s\n" % (name, keep))

    union = []
    for name in allpts:
        for p in allpts[name]:
            if not any(sp.simplify(p[0] - q[0]) == 0 and sp.simplify(p[1] - q[1]) == 0 for q in union):
                union.append(p)
    union.sort(key=lambda t: (-float(t[0]), -float(t[1])))

    perms = {"(2,-1,-1)/sqrt6": (2, -1, -1), "(-1,2,-1)/sqrt6": (-1, 2, -1), "(-1,-1,2)/sqrt6": (-1, -1, 2)}
    print("=== UNION over the three blocks: %d distinct points ===" % len(union))
    for (uu, vv) in union:
        ev = e_of(uu, vv)
        lab = [k for k, t in perms.items()
               if all(sp.simplify(ev[i] - sp.Integer(t[i]) / S6) == 0 for i in range(3))]
        print("  (u,v)=(%s,%s)  (a,b)=(%s,%s)  e=%s -> %s"
              % (uu, vv, sp.radsimp(S6 * uu), sp.radsimp(S6 * vv), ev, lab[0] if lab else "NOT a permutation"))

    print("\n=== which block attains L0 where ===")
    print("  %-14s %-24s %-5s %-5s %-5s" % ("(u,v)", "e", "M1", "M2", "M3"))
    for (uu, vv) in union:
        row = ["yes" if any(sp.simplify(uu - x) == 0 and sp.simplify(vv - y) == 0 for (x, y) in allpts[n])
               else "no" for n in ("M1", "M2", "M3")]
        print("  %-14s %-24s %-5s %-5s %-5s" % ("(%s, %s)" % (uu, vv), str(e_of(uu, vv)), *row))

    print("\n=== the two loci NAMED in the text, intersected with the ellipse ===")
    named = []
    # b = 0  <=> v = 0 ; a > 0 <=> u > 0
    for r in sp.roots(sp.Poly(4 * u**2 - 1, u)):
        if sp.N(r) > 0:
            named.append(("b = 0 with a > 0", sp.nsimplify(r), sp.Integer(0)))
    # a = -2b <=> u = -2v ; b > 0 <=> v > 0
    for r in sp.roots(sp.Poly(4 * (2 * v)**2 + 48 * v**2 - 1, v)):
        if sp.N(r) > 0:
            named.append(("a = -2b with b > 0", sp.nsimplify(-2 * r), sp.nsimplify(r)))
    for lab, uu, vv in named:
        print("  %-20s -> (u,v)=(%s,%s)  e=%s" % (lab, uu, vv, e_of(uu, vv)))
    print("  points delivered by the two named loci: %d" % len(named))
    missing = [p for p in union if not any(sp.simplify(p[0] - x) == 0 and sp.simplify(p[1] - y) == 0
                                           for (_, x, y) in named)]
    for p in missing:
        print("  NOT delivered by the named loci: (u,v)=(%s,%s), e=%s" % (p[0], p[1], e_of(*p)))

    print("\n=== lambda_max at the DISCARDED signs ===")
    for label, uu, vv in [("b=0, a<0", sp.Rational(-1, 2), sp.Integer(0)),
                          ("a=-2b, b<0 i.e. a=2b,b<0 reflected", sp.Rational(1, 4), sp.Rational(-1, 8))]:
        for name, M in BLOCKS.items():
            print("  %-36s %s: lam_max = %.15f   (L0 = %.15f)"
                  % (label, name, lam_max_num(M, float(uu), float(vv)), float(L0)))

    # ---------- independent dense sweep of the whole ellipse ----------
    print("\n=== independent dense numerical sweep of the whole ellipse ===")
    Jz = np.diag([float(m) for m in ORDER])
    Jp = np.zeros((7, 7))
    Jm = np.zeros((7, 7))
    for m in ORDER:
        if m + 1 <= 3:
            Jp[IDX[m + 1], IDX[m]] = np.sqrt(12 - m * (m + 1))
        if m - 1 >= -3:
            Jm[IDX[m - 1], IDX[m]] = np.sqrt(12 - m * (m - 1))
    FZ2 = Jz @ Jz
    PP = Jp @ Jp + Jm @ Jm
    L0f = 15 / np.sqrt(6)
    t = np.linspace(0, 2 * np.pi, 4000001)
    A = np.sqrt(1.5) * np.cos(t)
    B = np.sin(t) / (2 * np.sqrt(2))
    lam = np.empty_like(t)
    step = 200000
    for i in range(0, len(t), step):
        sl = slice(i, min(i + step, len(t)))
        Ops = (A[sl][:, None, None] * (FZ2 - 4 * np.eye(7))[None]
               + B[sl][:, None, None] * PP[None])
        lam[sl] = np.linalg.eigvalsh(Ops)[:, -1]
    print("  4,000,001 points on the ellipse; max lambda_max = %.15f   L0 = %.15f"
          % (lam.max(), L0f))
    print("  overshoot max(lam) - L0 = %.3e" % (lam.max() - L0f))
    near = np.where(L0f - lam < 1e-9)[0]
    # cluster the near-equality indices
    clusters = []
    for i in near:
        if clusters and i - clusters[-1][-1] <= 3:
            clusters[-1].append(i)
        else:
            clusters.append([i])
    if clusters and len(clusters) > 1 and clusters[0][0] == 0 and clusters[-1][-1] == len(t) - 1:
        clusters[0] = clusters[-1] + clusters[0]
        clusters.pop()
    print("  clusters of points within 1e-9 of L0: %d" % len(clusters))
    for cl in clusters:
        j = cl[len(cl) // 2]
        print("     t=%.9f  (a,b)=(%.12f, %.12f)  (u,v)=(%.9f, %.9f)  gap=%.2e"
              % (t[j], A[j], B[j], A[j] / np.sqrt(6), B[j] / np.sqrt(6), L0f - lam[j]))
    print("  exact points predicted, in (u,v): (0.5, 0), (-0.25, 0.125), (-0.25, -0.125)")


if __name__ == "__main__":
    main()
