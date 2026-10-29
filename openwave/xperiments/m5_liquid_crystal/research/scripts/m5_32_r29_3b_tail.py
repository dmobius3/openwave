"""M5.32 R29-3b: the author's Delta_dir definition and tail normalization (reply of 2026-10-28, 16:00 UTC), checked on
this stack's own F, and the covariant partial inverse at E0.

EQUATIONS FIRST
---------------
Uniaxial hedgehog M = diag(m00, S), S = lam_perp I + Delta n n^T, n = x/r, Delta_dir = lam_long - lam_perp, static.
F_ij = [d_i M, d_j M]_eta = d_i M eta d_j M - d_j M eta d_i M  (the time row is zero, so eta acts as the identity on S).
The author's radial check: I1 = sum_{i<j} <F_ij, F_ij>_eta = 2 Delta^4 / r^4 and C_tail = R E(r > R) = 32 pi Delta^4 for
E = 4 I1 (the static energy density of -4 I1 with the record's sign), and d ln C_tail = 4 (d lam_long - d lam_perp)/Delta.
The covariant partial inverse E# = eta^-1 E^T (E eta^-1 E^T)^+ at E0 = diag(sqrt g, 1, sqrt delta, 0) against diag(1/sqrt g, 1, 1/sqrt delta, 0).
Run: python3 m5_32_r29_3b_tail.py -> ../data/m5_32_r29_3b_tail.json
"""
import json, os, sys
import numpy as np, sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
x, y, z = sp.symbols("x y z", real=True)
lp, Dl = sp.symbols("lambda_perp Delta", real=True)
r = sp.sqrt(x**2 + y**2 + z**2)
n = sp.Matrix([x, y, z]) / r
S = lp * sp.eye(3) + Dl * n * n.T
M = sp.zeros(4, 4); M[0, 0] = sp.Symbol("m00"); M[1:, 1:] = S
ETA = sp.diag(-1, 1, 1, 1)
A = [M.diff(v) for v in (x, y, z)]
def comm_eta(P, Q): return P * ETA * Q - Q * ETA * P
def inner_eta(F, G): return (ETA * F * ETA * G.T).trace()
I1 = 0
for i in range(3):
    for j in range(i + 1, 3):
        F = comm_eta(A[i], A[j]); I1 += inner_eta(F, F)      # eta^i eta^j = +1 for spatial pairs
I1_r = sp.simplify(I1.subs({x: 0, y: 0}).subs(z, sp.Symbol("r", positive=True)))
rr = sp.Symbol("r", positive=True)
out = {"I1_on_the_axis_(sympy)": str(I1_r), "author_claim": "2 Delta^4 / r^4", "ratio_to_2Delta^4/r^4": str(sp.simplify(I1_r / (2 * Dl**4 / rr**4)))}
# generic point check numerically
f = sp.lambdify((x, y, z, lp, Dl), I1, "numpy")
pts = np.random.default_rng(3).normal(size=(5, 3)) * 3
vals = [float(f(*p, 0.3, 0.7)) * float(np.sum(p**2))**2 / (2 * 0.7**4) for p in pts]
out["generic_points_ratio_I1_r^4/(2Delta^4)"] = vals
# the tail: E = 4 I1 (the record's -4 I1 as a positive static energy density), C_tail = R int_R^inf 4 pi r^2 4 I1 dr
Rs = sp.Symbol("R", positive=True)
C_tail = sp.simplify(Rs * sp.integrate(4 * sp.pi * rr**2 * 4 * (2 * Dl**4 / rr**4), (rr, Rs, sp.oo)))
out["C_tail_=_R_E(r>R)"] = str(C_tail); out["author_claim_C_tail"] = "32 pi Delta^4"
out["d_ln_C_tail"] = "4 (d lam_long - d lam_perp)/Delta_dir: on the uniaxial R23 vacuum (1, 0.3, 0.3), Delta_dir = 0.7: d ln C_tail / d ln lam_long = 4/0.7 = %.4f, d ln C_tail / d ln lam_perp = -4 x 0.3/0.7 = %.4f" % (4 / 0.7, -4 * 0.3 / 0.7)
# E#
g, d = 32.0, 0.3
E0 = np.diag([np.sqrt(g), 1.0, np.sqrt(d), 0.0]); eta = np.diag([-1.0, 1, 1, 1])
Es = np.linalg.inv(eta) @ E0.T @ np.linalg.pinv(E0 @ np.linalg.inv(eta) @ E0.T)
out["E_sharp_at_E0"] = np.round(Es, 6).tolist()
out["diag_pseudo_inverse"] = np.round(np.diag([1 / np.sqrt(g), 1.0, 1 / np.sqrt(d), 0.0]), 6).tolist()
out["E_sharp_equals_diag_pinv_max_abs_diff"] = float(np.abs(Es - np.diag([1 / np.sqrt(g), 1.0, 1 / np.sqrt(d), 0.0])).max())
out["note_on_E_sharp"] = "eta^-1 E^T (E eta^-1 E^T)^+ : with E0 diagonal, E eta^-1 E^T = diag(-g, 1, delta, 0) and its pseudo-inverse diag(-1/g, 1, 1/delta, 0); eta^-1 E0^T times that gives diag(-1)(sqrt g)(-1/g) = 1/sqrt g on the time leg"
json.dump(out, open(os.path.join(HERE, "..", "data", "m5_32_r29_3b_tail.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
