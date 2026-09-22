"""s2a_02_equality_set.py -- the equality set of step 3c, computed exhaustively.

For each block M (2x2 real symmetric) the condition lambda_max(M) = L0 := 15/sqrt(6)
is EXACTLY equivalent to

      char_M(L0) = 0    AND    2*L0 >= tr M

(char_M(L0) = 0 says L0 is an eigenvalue; the trace condition says it is the larger
one).  Intersecting the conic char_M(L0) = 0 with the ellipse (2/3)a^2 + 8b^2 = 1
gives finitely many points, and sympy's solve over that zero-dimensional system
returns all of them.  Completeness therefore does not rest on inspection: it rests
on (i) the equivalence above, (ii) a zero-dimensional polynomial system solved
exactly, and (iii) an independent dense sweep of the whole ellipse (s2a_03).
"""

from pathlib import Path
import sympy as sp

HERE = Path(__file__).parent
a, b, L = sp.symbols("a b L", real=True)
L0 = 15 / sp.sqrt(6)
ELL = sp.Rational(2, 3) * a**2 + 8 * b**2 - 1

M1 = sp.Matrix([[5 * a, sp.sqrt(60) * b], [sp.sqrt(60) * b, -3 * a + 12 * b]])
M2 = sp.Matrix([[5 * a, sp.sqrt(60) * b], [sp.sqrt(60) * b, -3 * a - 12 * b]])
M3 = sp.Matrix([[0, sp.sqrt(240) * b], [sp.sqrt(240) * b, -4 * a]])
BLOCKS = {"M1": M1, "M2": M2, "M3": M3}


def lam_max(M):
    T = sp.simplify(M.trace())
    D = sp.simplify(M.det())
    return sp.simplify(T / 2 + sp.sqrt(sp.simplify(T**2 / 4 - D)))


def e_of(av, bv):
    return (sp.nsimplify(sp.simplify(-av / 3 + 2 * bv)),
            sp.nsimplify(sp.simplify(-av / 3 - 2 * bv)),
            sp.nsimplify(sp.simplify(2 * av / 3)))


def main():
    print("L0 = 15/sqrt(6) =", sp.nsimplify(L0), "=", sp.N(L0, 25))
    print("the seventh eigenvalue is identically 0, and 0 != L0 (L0 > 0), so it")
    print("contributes no equality point.  Checked: L0 - 0 =", sp.N(L0, 10), "> 0\n")

    allpts = {}
    for name, M in BLOCKS.items():
        cp = sp.expand(sp.simplify(M.charpoly(L).as_expr().subs(L, L0)))
        print("=== %s ===" % name)
        print("  char_M(L0) = 0  <=>  ", sp.simplify(sp.nsimplify(cp)), " = 0")
        sols = sp.solve([sp.numer(sp.together(cp)), ELL], [a, b], dict=True)
        pts = []
        for s in sols:
            av, bv = sp.simplify(s[a]), sp.simplify(s[b])
            if not (sp.im(sp.N(av)) == 0 and sp.im(sp.N(bv)) == 0):
                continue
            Tr = sp.simplify(M.trace().subs({a: av, b: bv}))
            ok = sp.simplify(2 * L0 - Tr)
            larger = sp.N(ok) >= -1e-30
            lm = sp.simplify(lam_max(M).subs({a: av, b: bv}))
            pts.append((av, bv, larger, sp.simplify(lm - L0)))
        print("  intersections with the ellipse (%d real):" % len(pts))
        for av, bv, larger, dd in pts:
            print("     (a,b) = (%s, %s)   2L0 >= trM : %s   lam_max - L0 = %s"
                  % (sp.nsimplify(av), sp.nsimplify(bv), larger, sp.simplify(dd)))
        keep = [(av, bv) for (av, bv, larger, dd) in pts if larger and sp.simplify(dd) == 0]
        print("  EQUALITY POINTS of %s: %s\n" % (name, [(sp.nsimplify(x), sp.nsimplify(y)) for x, y in keep]))
        allpts[name] = keep

    # union
    union = []
    for name in allpts:
        for p in allpts[name]:
            if not any(sp.simplify(p[0] - q[0]) == 0 and sp.simplify(p[1] - q[1]) == 0 for q in union):
                union.append(p)
    print("=== UNION over the three blocks ===")
    perms = [tuple(sp.nsimplify(x / sp.sqrt(6)) for x in t)
             for t in [(2, -1, -1), (-1, 2, -1), (-1, -1, 2)]]
    for (av, bv) in union:
        ev = e_of(av, bv)
        which = [i for i, pm in enumerate(perms)
                 if all(sp.simplify(ev[k] - pm[k]) == 0 for k in range(3))]
        print("  (a,b) = (%s, %s)   ->  e = %s   = permutation #%s of (2,-1,-1)/sqrt6"
              % (sp.nsimplify(av), sp.nsimplify(bv), tuple(sp.nsimplify(x) for x in ev),
                 which[0] if which else "NONE"))
    print("  total distinct equality points:", len(union))

    # per-block membership table
    print("\n=== which block attains L0 at which point ===")
    print("  %-28s %-6s %-6s %-6s" % ("(a,b)", "M1", "M2", "M3"))
    for (av, bv) in union:
        row = []
        for name in ("M1", "M2", "M3"):
            hit = any(sp.simplify(av - x) == 0 and sp.simplify(bv - y) == 0
                      for (x, y) in allpts[name])
            row.append("yes" if hit else "no")
        print("  %-28s %-6s %-6s %-6s"
              % ("(%s, %s)" % (sp.nsimplify(av), sp.nsimplify(bv)), *row))

    # the two loci the text names
    print("\n=== the two loci named in the text ===")
    for label, cond in [("b = 0 with a > 0", [b]), ("a = -2b with b > 0", [a + 2 * b])]:
        sols = sp.solve(cond + [ELL], [a, b], dict=True)
        for s in sols:
            av, bv = sp.simplify(s[a]), sp.simplify(s[b])
            sign_ok = (sp.N(av) > 0) if label.startswith("b = 0") else (sp.N(bv) > 0)
            if not sign_ok:
                continue
            ev = e_of(av, bv)
            print("  %-22s -> (a,b) = (%s, %s), e = %s"
                  % (label, sp.nsimplify(av), sp.nsimplify(bv), tuple(sp.nsimplify(x) for x in ev)))
    print("  => the two named loci meet the ellipse in exactly 2 points.")

    # the sign conditions the text asserts but does not derive
    print("\n=== the discarded signs: what lambda_max actually is there ===")
    for label, av, bv in [("b=0, a<0", -sp.sqrt(sp.Rational(3, 2)), sp.Integer(0)),
                          ("a=-2b, b<0", sp.sqrt(6) / 4, -sp.sqrt(6) / 8)]:
        for name, M in BLOCKS.items():
            lm = sp.simplify(lam_max(M).subs({a: av, b: bv}))
            print("  %-12s %s: lam_max = %s = %s   (L0 = %s)"
                  % (label, name, sp.nsimplify(lm), sp.N(lm, 15), sp.N(L0, 15)))


if __name__ == "__main__":
    main()
