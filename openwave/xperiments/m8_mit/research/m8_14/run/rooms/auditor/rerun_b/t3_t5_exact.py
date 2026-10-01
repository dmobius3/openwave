"""Exact spectrum (T3) and index / nullity (T5) from the closed form, R = 1.

lambda_{k,m} = (nu_k + m)(nu_k + m + 1),  nu_k = k*pi/(2W),  k >= 1, m >= 0.
Jacobi operator J = -Delta - 2, so its eigenvalues are lambda_{k,m} - 2
= (nu_k + m - 1)(nu_k + m + 2); negative iff nu_k + m < 1, zero iff = 1.
Decisions of sign are made by sympy on exact expressions in pi.
Numerical values are evaluated at 30 and at 50 digits with mpmath and must agree.
Pass --plant to use nu_k = k*pi/W (a planted defect); the PASS check against
the hemisphere value lambda_1(W = pi/2) = 2 must then print FAIL.
"""
import json
import os
import sys
import sympy as sp
import mpmath as mp

PLANT = "--plant" in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
WIDTHS = {"1/4": sp.Rational(1, 4), "1/2": sp.Rational(1, 2), "1": sp.Integer(1),
          "7/5": sp.Rational(7, 5), "3/2": sp.Rational(3, 2), "pi/2": sp.pi/2,
          "17/10": sp.Rational(17, 10)}


def nu(k, W):
    return k*sp.pi/W if PLANT else k*sp.pi/(2*W)


def lowest(W, n=6, K=12, M=12):
    labs = [(k, m) for k in range(1, K + 1) for m in range(M)]
    vals = {(k, m): sp.expand((nu(k, W) + m)*(nu(k, W) + m + 1)) for (k, m) in labs}
    order = sorted(labs, key=lambda km: float(vals[km].evalf(30)))
    return [(km, vals[km]) for km in order[:n]]


res = {}
for lab, W in WIDTHS.items():
    low = lowest(W)
    entry = {"eigs": []}
    for (k, m), v in low:
        mp.mp.dps = 30
        a = mp.mpf(str(sp.N(v, 40)))
        mp.mp.dps = 50
        b = mp.mpf(str(sp.N(v, 60)))
        assert abs(a - b) < mp.mpf(10)**-28, (lab, k, m)
        entry["eigs"].append({"k": k, "m": m, "exact": str(v), "value": float(b)})
    # T5: index and nullity of J = -Delta - 2 on the admissible class.
    # nu_k + m < 1 forces m = 0; k ranges over k < 2W/pi (finite), k <= 2W/pi + 1 suffices.
    kmax = int(sp.floor(2*W/sp.pi)) + 2
    neg = [(k, m) for k in range(1, kmax + 1) for m in range(2)
           if sp.simplify(nu(k, W) + m - 1).is_negative]
    zero = [(k, m) for k in range(1, kmax + 1) for m in range(2)
            if sp.simplify(nu(k, W) + m - 1) == 0]
    j1 = sp.simplify((nu(1, W) - 1)*(nu(1, W) + 2))
    sign = "zero" if j1 == 0 else ("negative" if j1.is_negative else "positive")
    entry.update({"J_lowest_exact": str(j1), "J_lowest": float(sp.N(j1, 30)),
                  "J_lowest_sign": sign, "index": len(neg), "nullity": len(zero),
                  "negative_labels": neg, "null_labels": zero})
    res[lab] = entry
    print(f"W={lab:6s} nu_1={float(sp.N(nu(1, W))):.12f}  lowest six: " +
          ", ".join(f"({e['k']},{e['m']}) {e['value']:.12f}" for e in entry["eigs"]))
    print(f"          lambda_1(J) = {j1} = {entry['J_lowest']:.12f} ({sign}); index {len(neg)}, nullity {len(zero)}")

ok = sp.simplify(lowest(sp.pi/2)[0][1] - 2) == 0
print(f"{'PASS' if ok else 'FAIL'}  closed form gives hemisphere value lambda_1(pi/2) = 2")
if not PLANT:
    with open(os.path.join(HERE, "exact.json"), "w") as fh:
        json.dump(res, fh, indent=1)
