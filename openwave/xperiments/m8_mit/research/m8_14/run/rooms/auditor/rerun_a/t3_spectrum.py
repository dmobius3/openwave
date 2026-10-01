"""T2/T3/T5: symbolic checks of the separated problem and the closed-form spectrum (R = 1).

Lune picture (see RETURN.md): latitude psi in (-pi/2, pi/2), x = sin(psi), transverse sector k:
  e_k(w) = sin(nu_k (w + W)),  nu_k = k*pi/(2W)
  f_{k,n} = (1 - x^2)^(nu/2) C_n^(nu+1/2)(x),  lambda_{k,n} = (nu_k + n)(nu_k + n + 1).
Checks (each can fail; see --plant):
  C3.1  -((1-x^2) f')' + nu^2 f/(1-x^2) - lambda f == 0 symbolically in nu, n = 0..5.
  C3.2  the full separated function u = f(sin psi) sin(nu(w+W)) satisfies
        -Delta u = lambda u for ds^2 = dpsi^2 + cos^2(psi) dw^2, symbolic in nu, n = 0..3.
  C3.3  the y-picture function (phi = u on y<=pi/2, phi = -u(y-pi, -w) on y>pi/2) satisfies the
        seam condition phi(pi,-w) = -phi(0,w) with all y-derivatives up to order 3, k = 1..4, n = 0..2.
Plants: lam (use (nu+n)(nu+n+2)), seam (drop the minus sign in the lower branch).
Writes t3_out.json.
"""
import sys, os, json
import sympy as sp
import mpmath as mp

plant = sys.argv[sys.argv.index("--plant") + 1] if "--plant" in sys.argv else None
here = os.path.dirname(os.path.abspath(__file__))

x, psi, w, y = sp.symbols("x psi w y", real=True)
nu, W = sp.symbols("nu W", positive=True)


def lam(nuv, n):
    if plant == "lam":
        return (nuv + n) * (nuv + n + 2)
    return (nuv + n) * (nuv + n + 1)


ok_all = True


def report(name, ok):
    global ok_all
    ok_all &= ok
    print(f"{name}: {'PASS' if ok else 'FAIL'}")


ok = True
for n in range(6):
    g = sp.gegenbauer(n, nu + sp.Rational(1, 2), x)
    f = (1 - x**2) ** (nu / 2) * g
    res = -sp.diff((1 - x**2) * sp.diff(f, x), x) + nu**2 * f / (1 - x**2) - lam(nu, n) * f
    res = sp.simplify(sp.expand(sp.simplify(res / (1 - x**2) ** (nu / 2))))
    ok &= res == 0
report("C3.1 Gegenbauer eigenfunctions of the sector ODE, n=0..5, symbolic nu", ok)

ok = True
for n in range(4):
    g = sp.gegenbauer(n, nu + sp.Rational(1, 2), sp.sin(psi))
    u = sp.cos(psi) ** nu * g * sp.sin(nu * (w + W))
    lap = sp.diff(sp.cos(psi) * sp.diff(u, psi), psi) / sp.cos(psi) + sp.diff(u, w, 2) / sp.cos(psi) ** 2
    res = sp.simplify((-lap - lam(nu, n) * u) / (sp.cos(psi) ** nu * sp.sin(nu * (w + W))))
    ok &= sp.simplify(sp.expand_trig(res)) == 0
report("C3.2 separated eigenfunctions solve -Delta u = lambda u in (psi, w), n=0..3", ok)

ok = True
Wn = sp.Rational(7, 5)
for k in range(1, 5):
    nuk = k * sp.pi / (2 * Wn)
    for n in range(3):
        def u(p, ww):
            return sp.cos(p) ** nuk * sp.gegenbauer(n, nuk + sp.Rational(1, 2), sp.sin(p)) * sp.sin(nuk * (ww + Wn))
        upper = u(y, w)
        lower = (1 if plant == "seam" else -1) * u(y - sp.pi, -w)
        for d in range(4):
            lhs = sp.diff(lower, y, d).subs(y, sp.pi).subs(w, -w)
            rhs = -sp.diff(upper, y, d).subs(y, 0)
            for wv in (sp.Rational(-13, 10), sp.Rational(1, 3), sp.Rational(9, 10)):
                val = sp.N((lhs - rhs).subs(w, wv), 30)
                ok &= abs(val) < 1e-25
report("C3.3 seam condition with derivatives (k=1..4, n=0..2, d=0..3, W=7/5)", ok)

print("T3_ALL", "PASS" if ok_all else "FAIL")

# ---- closed-form numbers at the widths ----
widths = {"1/4": sp.Rational(1, 4), "1/2": sp.Rational(1, 2), "1": sp.Integer(1), "7/5": sp.Rational(7, 5),
          "3/2": sp.Rational(3, 2), "pi/2": sp.pi / 2, "17/10": sp.Rational(17, 10)}
out = {}
for lab, Wv in widths.items():
    levels = []
    for k in range(1, 40):
        for n in range(40):
            nuk = k * sp.pi / (2 * Wv)
            levels.append((float(lam(nuk, n)), k, n, sp.nsimplify(lam(nuk, n))))
    levels.sort(key=lambda t: t[0])
    six = levels[:6]
    nu1 = sp.pi / (2 * Wv)
    lc = [k for k in range(1, 40) if (k * sp.pi / (2 * Wv) < 1) == True]
    jac = [sp.N(lam(k * sp.pi / (2 * Wv), n) - 2, 20) for (_, k, n, _) in six]
    neg = sum(1 for (v, k, n, s) in levels if sp.N(s - 2, 50) < 0)
    nul = sum(1 for (v, k, n, s) in levels if sp.simplify(s - 2) == 0)
    out[lab] = {
        "nu1": float(sp.N(nu1, 20)),
        "nu1_exact": str(nu1),
        "limit_circle_sectors": lc,
        "lowest6": [float(sp.N(s, 20)) for (_, k, n, s) in six],
        "lowest6_labels": [[k, n] for (_, k, n, _) in six],
        "lowest6_exact": [str(sp.factor(s)) for (_, k, n, s) in six],
        "bottom_exact": str(sp.factor(six[0][3])),
        "jacobi_lowest": float(jac[0]),
        "jacobi_lowest6": [float(j) for j in jac],
        "jacobi_negative_count": neg,
        "jacobi_nullity": nul,
    }
    print(f"W={lab:6s} nu1={out[lab]['nu1']:.15f} LC sectors={lc} bottom={out[lab]['lowest6'][0]:.15f} "
          f"labels={out[lab]['lowest6_labels']} jac_low={out[lab]['jacobi_lowest']:+.15f} neg={neg} null={nul}")
    print("        lowest6 =", ", ".join(f"{v:.12f}" for v in out[lab]["lowest6"]))

if plant is None:
    with open(os.path.join(here, "t3_out.json"), "w") as fh:
        json.dump(out, fh, indent=1)
