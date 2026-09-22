"""s09_maximise.py -- maximise rhat_6 over the three-axes domain.

From s07 (verified exactly at three special configurations and to 1e-16 against
direct state evaluation in s08):

   rhat_6 = 20 * N / (77 * Den^2)   on the Theta-invariant (real-form) states,
   Den = 3 q1 - 2 p + 5
   N   = 15 q1^2 + 7 q2 - 20 p q1 + 9 p^2 + 71 q1 - 108 p + 30
   with a = n1.n2, b = n1.n3, c = n2.n3,
        q1 = a^2+b^2+c^2, q2 = a^2b^2+b^2c^2+c^2a^2, p = abc,
   over the ELLIPTOPE  E = { det Gram = 1 + 2p - q1 >= 0 }.

Stage 1: numerics (corroboration only).
Stage 2: exact critical-point analysis, interior and boundary.
"""

from pathlib import Path
import itertools
import numpy as np
import sympy as sp

HERE = Path(__file__).parent
a, b, c = sp.symbols("a b c", real=True)

q1 = a**2 + b**2 + c**2
q2 = a**2 * b**2 + b**2 * c**2 + c**2 * a**2
p = a * b * c
Den = 3 * q1 - 2 * p + 5
Nn = 15 * q1**2 + 7 * q2 - 20 * p * q1 + 9 * p**2 + 71 * q1 - 108 * p + 30
R = sp.Rational(20, 77) * Nn / Den**2
detG = 1 + 2 * p - q1


def main():
    Rf = sp.lambdify((a, b, c), R, "numpy")
    Df = sp.lambdify((a, b, c), detG, "numpy")

    # ---------- stage 1: numerics over the elliptope ----------
    print("STAGE 1 (corroboration only): grid + random search over the elliptope")
    rng = np.random.default_rng(2)
    best = (-1, None)
    # random unit vectors -> guaranteed feasible Gram points
    for _ in range(400000):
        V = rng.normal(size=(3, 3))
        V /= np.linalg.norm(V, axis=1)[:, None]
        g = (V[0] @ V[1], V[0] @ V[2], V[1] @ V[2])
        v = float(Rf(*g))
        if v > best[0]:
            best = (v, g)
    print("  random axes search max = %.15f  at Gram %s" % (best[0], np.round(best[1], 8)))
    print("  463/924 = %.15f" % (463 / 924))

    # local polish by projected coordinate search on axis angles
    def val_from_axes(x):
        V = x.reshape(3, 3)
        V = V / np.linalg.norm(V, axis=1)[:, None]
        return -float(Rf(V[0] @ V[1], V[0] @ V[2], V[1] @ V[2]))

    from scipy.optimize import minimize
    tops = []
    for t in range(300):
        x0 = rng.normal(size=9)
        r = minimize(val_from_axes, x0, method="Nelder-Mead",
                     options={"maxiter": 20000, "xatol": 1e-12, "fatol": 1e-14})
        V = r.x.reshape(3, 3)
        V = V / np.linalg.norm(V, axis=1)[:, None]
        g = (V[0] @ V[1], V[0] @ V[2], V[1] @ V[2])
        tops.append((-r.fun, g, float(Df(*g))))
    tops.sort(key=lambda t: -t[0])
    print("\n  top local maxima (value, Gram, detG):")
    for v, g, d in tops[:8]:
        print("    %.15f   %s   detG=%.2e" % (v, np.round(np.abs(np.sort(np.abs(g))), 10), d))
    print("  distinct plateaus:", sorted({round(t[0], 10) for t in tops}, reverse=True)[:6])

    # what are |a|,|b|,|c| at the best points?
    print("\n  |Gram entries| at the best points (sorted):")
    for v, g, d in tops[:5]:
        print("    %.12f  %s" % (v, np.round(np.sort(np.abs(g)), 10)))

    # ---------- stage 2: exact critical points ----------
    print("\nSTAGE 2: exact critical-point analysis")
    print("  R = 20 N / (77 Den^2), N and Den as above")

    # interior critical points: grad R = 0
    gr = [sp.simplify(sp.together(sp.diff(R, v))) for v in (a, b, c)]
    nums = [sp.numer(sp.together(g)) for g in gr]
    print("  interior: solving the three numerators of dR/d(a,b,c) = 0 ...")
    gb = sp.groebner(nums, a, b, c, order="lex")
    print("  Groebner basis length:", len(gb.exprs))
    sols = sp.solve(nums, [a, b, c], dict=True)
    print("  number of solutions found by solve():", len(sols))
    vals = []
    for s in sols:
        if all(v.is_real for v in s.values()) or True:
            try:
                av, bv, cv = [sp.nsimplify(s.get(x, x)) for x in (a, b, c)]
            except Exception:
                continue
            if any(x.free_symbols for x in (av, bv, cv)):
                continue
            if not all(sp.im(sp.N(x)) == 0 for x in (av, bv, cv)):
                continue
            d = sp.simplify(detG.subs({a: av, b: bv, c: cv}))
            if sp.N(d) < -1e-12:
                continue
            r = sp.simplify(R.subs({a: av, b: bv, c: cv}))
            vals.append((sp.N(r, 20), r, (av, bv, cv), sp.N(d, 10)))
    vals.sort(key=lambda t: -t[0])
    print("  feasible interior critical points (value, exact, point, detG):")
    for v in vals:
        print("   ", v)

    with open(HERE / "s09_out.txt", "w") as fh:
        fh.write("numeric max %.15f\n" % best[0])
        for v in vals:
            fh.write(str(v) + "\n")


if __name__ == "__main__":
    main()
