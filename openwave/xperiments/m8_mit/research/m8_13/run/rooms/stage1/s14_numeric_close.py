"""s14_numeric_close.py -- fast numerical closure and a second-precision check.

(s09_maximise.py was left running a Groebner computation that produced no output
within the session budget; it is kept but its role is superseded by the exact
chain in s12.  This script supplies the numerical corroboration s09 stage 1 was
meant to give, and re-evaluates every reported number at a second precision.)
"""

from pathlib import Path
import numpy as np
import sympy as sp
import mpmath as mp

HERE = Path(__file__).parent
a, b, c = sp.symbols("a b c", real=True)
q1 = a**2 + b**2 + c**2
q2 = a**2 * b**2 + b**2 * c**2 + c**2 * a**2
p = a * b * c
Den = 3 * q1 - 2 * p + 5
Nn = 15 * q1**2 + 7 * q2 - 20 * p * q1 + 9 * p**2 + 71 * q1 - 108 * p + 30
R = sp.Rational(20, 77) * Nn / Den**2
Rf = sp.lambdify((a, b, c), R, "numpy")


def gram(V):
    return V[0] @ V[1], V[0] @ V[2], V[1] @ V[2]


def main():
    rng = np.random.default_rng(404)

    print("local refinement over axis triples (gradient-free, 200 restarts):")
    from scipy.optimize import minimize

    def neg(x):
        V = x.reshape(3, 3)
        V = V / np.linalg.norm(V, axis=1)[:, None]
        return -float(Rf(*gram(V)))

    vals = []
    for _ in range(200):
        r = minimize(neg, rng.normal(size=9), method="Powell",
                     options={"maxiter": 100000, "xtol": 1e-13, "ftol": 1e-15})
        V = r.x.reshape(3, 3)
        V = V / np.linalg.norm(V, axis=1)[:, None]
        vals.append((-r.fun, gram(V)))
    vals.sort(key=lambda t: -t[0])
    print("  best  %.15f   463/924 = %.15f" % (vals[0][0], 463 / 924))
    print("  top-8 values:", ["%.12f" % v for v, _ in vals[:8]])
    print("  distinct plateaus (9dp):", sorted({round(v, 9) for v, _ in vals}, reverse=True)[:8])
    print("  |Gram| at the best point:", np.round(np.sort(np.abs(vals[0][1])), 10))
    print("  detG at the best point: %.3e"
          % (1 + 2 * np.prod(vals[0][1]) - sum(x**2 for x in vals[0][1])))

    print("\n  second plateau (if any) and its Gram:")
    seen = []
    for v, g in vals:
        if all(abs(v - w) > 1e-8 for w in seen):
            seen.append(v)
            print("    %.12f   |Gram| = %s" % (v, np.round(np.sort(np.abs(g)), 8)))
        if len(seen) >= 6:
            break

    print("\nsecond-precision (mpmath, 40 digits) values of every reported number:")
    mp.mp.dps = 40
    for name, val in [("maximum 463/924", mp.mpf(463) / 924),
                      ("minimum 1/924", mp.mpf(1) / 924),
                      ("<3 3;3 -3|6 0> = 1/sqrt(924)", 1 / mp.sqrt(924)),
                      ("its square 1/924", mp.mpf(1) / 924),
                      ("rhat_0 = 1/7", mp.mpf(1) / 7),
                      ("bound 6/7 (not attained)", mp.mpf(6) / 7),
                      ("LP bound 11/21 (not attained)", mp.mpf(11) / 21),
                      ("n_0 = 13/7", mp.mpf(13) / 7),
                      ("n_2 = 65/84", mp.mpf(65) / 84),
                      ("n_4 = 13/154", mp.mpf(13) / 154),
                      ("n_6 = 1/924", mp.mpf(1) / 924),
                      ("hexagon rhat_2 = 25/84", mp.mpf(25) / 84),
                      ("hexagon rhat_4 = 9/154", mp.mpf(9) / 154)]:
        print("   %-34s %s" % (name, mp.nstr(val, 30)))

    print("\n  mpmath re-evaluation of R at the maximising Gram (-1/2,-1/2,-1/2):")
    A = mp.mpf(-1) / 2
    Q1 = 3 * A**2
    Q2 = 3 * A**4
    P = A**3
    D = 3 * Q1 - 2 * P + 5
    N = 15 * Q1**2 + 7 * Q2 - 20 * P * Q1 + 9 * P**2 + 71 * Q1 - 108 * P + 30
    print("    R = %s   vs 463/924 = %s"
          % (mp.nstr(mp.mpf(20) / 77 * N / D**2, 30), mp.nstr(mp.mpf(463) / 924, 30)))
    print("  mpmath at the coincident Gram (1,1,1) (state v_0):")
    Q1, Q2, P = mp.mpf(3), mp.mpf(3), mp.mpf(1)
    D = 3 * Q1 - 2 * P + 5
    N = 15 * Q1**2 + 7 * Q2 - 20 * P * Q1 + 9 * P**2 + 71 * Q1 - 108 * P + 30
    print("    R = %s   vs 100/231 = %s"
          % (mp.nstr(mp.mpf(20) / 77 * N / D**2, 30), mp.nstr(mp.mpf(100) / 231, 30)))


if __name__ == "__main__":
    main()
