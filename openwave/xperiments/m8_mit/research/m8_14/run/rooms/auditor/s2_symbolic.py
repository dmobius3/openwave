"""Stage 2: independent checks of the claims used in grading T1, T2, T3 and T5 (R = 1).

Run as `./py s2_symbolic.py`.  Writes s2_symbolic.json.  Each PASS line has a planted twin
(marked [planted]) that must print FAIL.

  G1  T1 conclusion, computed independently of both rooms: the area of the varied surface
      F_t = cos(t phi) X + sin(t phi) e4, built in R^4 from the section-2 parametrization
      (Gram determinant of finite-difference tangents, Gauss quadrature on both halves), has
      A''(0) = Q(phi) = int (|grad phi|^2 - 2 phi^2) dA for an admissible test phi.
      Planted: compare with int |grad phi|^2 (no curvature term).
  G2  Room A T3 step 2: ((1-x^2) f')' for f = (1-x^2)^(nu/2) g equals the stated expression.
      Planted: drop the -nu g term.
  G3  sector ODE, independently: -((1-x^2) f')' + nu^2 f/(1-x^2) = lam f for f = (1-x^2)^(nu/2) C_m^(nu+1/2),
      lam = (nu+m)(nu+m+1), m = 0..4, at exact rational nu (and nu = 5*pi/17 numerically).
      Planted: lam = (nu+m)(nu+m+2).
  G4  Frobenius / limit-circle thresholds: indicial roots of the weighted sector operator are +-nu;
      with the weight cos psi removed from operator and measure, the roots are (1 +- sqrt(1+4 nu^2))/2
      and the limit-circle threshold is nu < sqrt(3)/2 (bears on Room A's side remark "nu < 1/2").
  G5  T3/T5 values: the lowest six (k,m) of (nu_k+m)(nu_k+m+1) agree with my stage-1 closed form
      (results_C.json) and with both rooms' JSON; sign, index and nullity of lam - 2, decided exactly.
      Planted: nu_k = k pi / W.
  G6  Room B section 0 remark (i): at W = pi/2 the two edges w = +-W are the same set in RP^3.
      Planted: W = 1.4 (they must then differ).
"""
import json
import math
import os

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
out, lines = {}, []


def report(name, ok, planted=False):
    s = f"{'[planted] ' if planted else ''}{'PASS' if ok else 'FAIL'}  {name}"
    print(s, flush=True)
    lines.append(s)
    return ok


# ---------------------------------------------------------------- G1
def bump(s, a=1.2):
    s = np.asarray(s, float)
    v = np.zeros_like(s)
    m = np.abs(s) < a
    v[m] = np.exp(-1.0 / (1 - (s[m] / a) ** 2))
    return v


def make_phi(W):
    k = math.pi / (2 * W)
    def u(psi, th):                      # smooth on the lune, zero near both poles, Dirichlet at th = +-W
        return bump(psi) * (1 + 0.3 * np.sin(psi) + 0.2 * np.cos(2 * psi)) * np.sin(k * (th + W)) \
            * (1 + 0.25 * np.sin(2 * k * (th + W)))
    def phi(y, w):
        return np.where(y <= math.pi / 2, u(y, w), -u(y - math.pi, -w))
    return phi


def Xmap(y, w):
    """section-2 parametrization on each half: X(y, w) for y <= pi/2, X(y, -w) for y > pi/2."""
    ws = np.where(y <= math.pi / 2, w, -w)
    return np.stack([np.cos(y) * np.cos(ws), np.cos(y) * np.sin(ws), np.sin(y), 0 * y])


def area(t, phi, Y, Wg, wq, h=1e-6):
    def F(y, w):
        p = phi(y, w)
        return np.cos(t * p) * Xmap(y, w) + np.sin(t * p) * np.array([0, 0, 0, 1.0])[:, None, None]
    Fy = (F(Y + h, Wg) - F(Y - h, Wg)) / (2 * h)
    Fw = (F(Y, Wg + h) - F(Y, Wg - h)) / (2 * h)
    E, G, Fm = (Fy * Fy).sum(0), (Fw * Fw).sum(0), (Fy * Fw).sum(0)
    return float(np.sum(wq * np.sqrt(np.maximum(E * G - Fm**2, 0))))


def g1(W, planted=False):
    phi = make_phi(W)
    gx, gw = np.polynomial.legendre.leggauss(160)
    ys, wys = [], []
    for a, b in ((0, math.pi / 2), (math.pi / 2, math.pi)):
        ys.append(0.5 * (b - a) * gx + 0.5 * (a + b)); wys.append(0.5 * (b - a) * gw)
    yq, wy = np.concatenate(ys), np.concatenate(wys)
    wq_, ww = W * gx, W * gw
    Y, Wg = np.meshgrid(yq, wq_, indexing="ij")
    wq = wy[:, None] * ww[None, :]
    A0 = area(0.0, phi, Y, Wg, wq)
    d2 = []
    for t in (0.02, 0.01):
        d2.append((area(t, phi, Y, Wg, wq) - 2 * A0 + area(-t, phi, Y, Wg, wq)) / t**2)
    A2 = d2[1] + (d2[1] - d2[0]) / 3                 # Richardson, O(t^2) error
    h = 1e-6
    py = (phi(Y + h, Wg) - phi(Y - h, Wg)) / (2 * h)
    pw = (phi(Y, Wg + h) - phi(Y, Wg - h)) / (2 * h)
    c = np.abs(np.cos(Y))
    P = phi(Y, Wg)
    grad2 = py**2 + np.where(c > 1e-300, pw**2 / np.maximum(c, 1e-300) ** 2, 0)
    Q = float(np.sum(wq * (grad2 - (0 if planted else 2) * P**2) * c))
    seam = float(np.max(np.abs(phi(np.full(5, math.pi), -np.linspace(-W, W, 5)) + phi(np.zeros(5), np.linspace(-W, W, 5)))))
    rel = abs(A2 - Q) / abs(Q)
    report(f"G1 W={W}: area {A0:.12f} (4W = {4*W}), A''(0) = {A2:.9f}, Q = {Q:.9f}, rel {rel:.1e}, "
           f"test phi seam residual {seam:.1e}", rel < 1e-5 and abs(A0 - 4 * W) < 1e-9, planted)
    return {"W": W, "area0": A0, "A2": A2, "Q": Q, "rel": rel}


out["G1"] = [g1(1.4), g1(1.7)]
out["G1_planted_fired"] = not g1(1.4, planted=True)["rel"] < 1e-5

# ---------------------------------------------------------------- G2, G3
x, nu = sp.symbols("x nu", positive=True)
g = sp.Function("g")(x)
f = (1 - x**2) ** (nu / 2) * g


def g2(planted=False):
    lhs = sp.diff((1 - x**2) * sp.diff(f, x), x)
    rhs = (1 - x**2) ** (nu / 2) * ((1 - x**2) * g.diff(x, 2) - 2 * (nu + 1) * x * g.diff(x)
                                    - (0 if planted else nu) * g + nu**2 * x**2 * g / (1 - x**2))
    ok = sp.simplify((lhs - rhs) / (1 - x**2) ** (nu / 2)) == 0
    return report("G2 Room A T3 step 2 identity for ((1-x^2) f')'", ok, planted)


out["G2"] = g2()
out["G2_planted_fired"] = not g2(True)


def g3(planted=False):
    ok = True
    for nv in (sp.Rational(3, 7), sp.Rational(1), sp.Rational(25, 11)):
        for m in range(5):
            F = (1 - x**2) ** (nv / 2) * sp.gegenbauer(m, nv + sp.Rational(1, 2), x)
            lam = (nv + m) * (nv + m + (2 if planted else 1))
            r = -sp.diff((1 - x**2) * sp.diff(F, x), x) + nv**2 * F / (1 - x**2) - lam * F
            ok &= sp.simplify(r / (1 - x**2) ** (nv / 2)) == 0
    nn = 5 * math.pi / 17
    for m in range(4):
        F = (1 - x**2) ** (sp.Float(nn, 30) / 2) * sp.gegenbauer(m, sp.Float(nn, 30) + sp.Rational(1, 2), x)
        lam = (nn + m) * (nn + m + (2 if planted else 1))
        r = -sp.diff((1 - x**2) * sp.diff(F, x), x) + nn**2 * F / (1 - x**2) - lam * F
        for xv in (-0.7, 0.1, 0.55):
            ok &= abs(float(r.subs(x, xv))) < 1e-20 * 1e6
    return report("G3 sector ODE solved by (1-x^2)^(nu/2) C_m^(nu+1/2) with lam=(nu+m)(nu+m+1)", ok, planted)


out["G3"] = g3()
out["G3_planted_fired"] = not g3(True)

# ---------------------------------------------------------------- G4
r, s = sp.symbols("r s")
ind_w = sp.expand((-(1 / r) * sp.diff(r * sp.diff(r**s, r), r) + nu**2 * r**s / r**2) * r**(2 - s))
roots_w = sorted(sp.solve(sp.Eq(ind_w, 0), s), key=str)
ind_u = sp.expand((-sp.diff(r**s, r, 2) + nu**2 * r**s / r**2) * r**(2 - s))
roots_u = sp.solve(sp.Eq(ind_u, 0), s)
sm = min(roots_u, key=lambda e: float(e.subs(nu, 1)))
thr = sp.solve(sp.Eq(2 * sm + 1, 0), nu)            # r^{2 s_-} integrable in dr iff 2 s_- > -1
report(f"G4 weighted indicial roots {roots_w}; unweighted roots {roots_u}; unweighted LC threshold nu < {thr}",
       set(roots_w) == {nu, -nu} and [sp.nsimplify(t) for t in thr] == [sp.sqrt(3) / 2])
# planted twin: the 3D-radial operator -(1/r^2)(r^2 f')' + nu^2 f/r^2 in place of the weighted
# sector operator; its roots are not +-nu, so the line must FAIL
ind_p = sp.expand((-(1 / r**2) * sp.diff(r**2 * sp.diff(r**s, r), r) + nu**2 * r**s / r**2) * r**(2 - s))
roots_p = sp.solve(sp.Eq(ind_p, 0), s)
out["G4_planted_fired"] = not report(f"G4 planted operator: indicial roots {roots_p} vs +-nu", set(roots_p) == {nu, -nu}, True)
out["G4"] = {"weighted_roots": [str(e) for e in roots_w], "unweighted_roots": [str(e) for e in roots_u],
             "unweighted_LC_threshold": [str(t) for t in thr],
             "measure_only_change_threshold": "nu < 1/2 (int r^{-2nu} dr < inf), roots kept at +-nu"}

# ---------------------------------------------------------------- G5
WID = {"1/4": sp.Rational(1, 4), "1/2": sp.Rational(1, 2), "1": sp.Integer(1), "7/5": sp.Rational(7, 5),
       "3/2": sp.Rational(3, 2), "pi/2": sp.pi / 2, "17/10": sp.Rational(17, 10)}
mineC = json.load(open(os.path.join(HERE, "results_C.json")))
ra = json.load(open(os.path.join(HERE, "room_a", "results.json")))
rb = json.load(open(os.path.join(HERE, "room_b", "results.json")))


def find_mine(lab):
    """my stage-1 closed-form lowest six for width lab (structure of results_C.json probed)."""
    d = mineC
    for key in ("widths", "closed_form"):
        if isinstance(d, dict) and key in d:
            d = d[key]
    e = d[lab]
    for key in ("lowest", "lowest6", "values", "eigs"):
        if isinstance(e, dict) and key in e:
            return [float(v) for v in e[key]]
    raise KeyError(lab)


def g5(planted=False):
    ok, tab = True, {}
    for lab, W in WID.items():
        nuk = lambda k: k * sp.pi / (W if planted else 2 * W)
        labs = [(k, m) for k in range(1, 14) for m in range(14)]
        val = {km: (nuk(km[0]) + km[1]) * (nuk(km[0]) + km[1] + 1) for km in labs}
        order = sorted(labs, key=lambda km: float(sp.N(val[km], 30)))
        six = [float(sp.N(val[km], 30)) for km in order[:6]]
        mine = find_mine(lab)
        a6 = ra["t3_closed_form"][lab]["lowest6"]
        b6 = [e["value"] for e in rb["t3_t5_exact"][lab]["eigs"]]
        dev = max(max(abs(p - q) / q for p, q in zip(six, mine)),
                  max(abs(p - q) / q for p, q in zip(six, a6)),
                  max(abs(p - q) / q for p, q in zip(six, b6)))
        j = sp.simplify(val[(1, 0)] - 2)
        sign = 0 if j == 0 else (1 if j.is_positive else -1)
        neg = sum(1 for km in labs if sp.simplify(val[km] - 2).is_negative)
        nul = sum(1 for km in labs if sp.simplify(val[km] - 2) == 0)
        ta = ra["t5"][lab]; tb = rb["t3_t5_exact"][lab]
        agree = (ta["negative_count"], ta["nullity"]) == (neg, nul) == (tb["index"], tb["nullity"])
        ok &= dev < 1e-12 and agree
        tab[lab] = {"lowest6": six, "labels": order[:6], "jacobi_lowest": float(sp.N(j, 20)), "sign": sign,
                    "index": neg, "nullity": nul, "max_rel_dev_vs_mine_and_rooms": dev, "rooms_T5_agree": agree}
    report("G5 lowest six vs my stage-1 closed form and both rooms; T5 sign/index/nullity vs both rooms", ok, planted)
    return ok, tab


ok5, out["G5"] = g5()
okp5, _ = g5(True)
out["G5_planted_fired"] = not okp5
for lab, v in out["G5"].items():
    print(f"   W={lab:6s} J_low={v['jacobi_lowest']:+.12f} sign {v['sign']:+d} index {v['index']} nullity {v['nullity']}")

# ---------------------------------------------------------------- G6
def g6(W, expect_same):
    yv = np.linspace(-math.pi / 2, math.pi / 2, 401)
    E1 = np.stack([np.cos(yv) * math.cos(W), np.cos(yv) * math.sin(W), np.sin(yv)], 1)
    E2 = np.stack([np.cos(yv) * math.cos(-W), np.cos(yv) * math.sin(-W), np.sin(yv)], 1)
    def dist(P, S):      # max over P of distance to the set S up to +-
        d = np.minimum(np.linalg.norm(P[:, None] - S[None], axis=2), np.linalg.norm(P[:, None] + S[None], axis=2))
        return float(d.min(1).max())
    d = max(dist(E1, E2), dist(E2, E1))
    return report(f"G6 W={W:.6f}: Hausdorff distance in RP^3 between the edges w=+W and w=-W: {d:.1e}",
                  (d < 1e-2) == expect_same and (d < 1e-2), not expect_same) , d


ok6, d6 = g6(math.pi / 2, True)
okp6, d6p = g6(1.4, False)
out["G6"] = {"pi/2": d6, "1.4": d6p}
out["G6_planted_fired"] = not okp6

out["lines"] = lines
with open(os.path.join(HERE, "s2_symbolic.json"), "w") as fh:
    json.dump(out, fh, indent=1, default=str)
