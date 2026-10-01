"""Symbolic checks for T1, T2, T3 (R = 1 throughout, but T1 keeps R symbolic).

Every PASS line below compares a computed symbolic expression with zero.
Pass the argument --plant to plant a defect in each check (wrong curvature
constant, wrong eigenvalue, wrong exponent, wrong weight); the planted
checks must then print FAIL.
"""
import os
import sys
import sympy as sp

PLANT = "--plant" in sys.argv
out = []


def check(name, expr):
    ok = sp.simplify(expr) == 0
    line = f"{'PASS' if ok else 'FAIL'}  {name}"
    print(line)
    out.append(line)
    return ok


# ---------------------------------------------------------------- T1
# Varied map F = cos(u) X + R sin(u) e4, u = t*phi/R, X on S^3(R) with metric
# dy^2 + cos^2(y/R) dw^2.  Compute the induced metric from the 4-vector map.
y, w, t, R = sp.symbols("y w t R", positive=True)
phi = sp.Function("phi")(y, w)
X = sp.Matrix([R*sp.cos(y/R)*sp.cos(w/R), R*sp.cos(y/R)*sp.sin(w/R), R*sp.sin(y/R), 0])
e4 = sp.Matrix([0, 0, 0, 1])
u = t*phi/R
F = sp.cos(u)*X + R*sp.sin(u)*e4
Fy, Fw = F.diff(y), F.diff(w)
g = sp.Matrix([[Fy.dot(Fy), Fy.dot(Fw)], [Fw.dot(Fy), Fw.dot(Fw)]])
detg = sp.simplify(g.det())
# hand formula: det = cos^2(u) cos^2(y/R) [cos^2(u) + t^2 |grad phi|^2]
c = sp.cos(y/R)
grad2 = phi.diff(y)**2 + phi.diff(w)**2/c**2
hand = (1 if PLANT else sp.cos(u)**2)*c**2*(sp.cos(u)**2 + t**2*grad2)   # planted: drop cos^2 u
check("T1 induced-metric determinant equals hand formula", sp.expand(sp.simplify(detg - hand)))
# second t-derivative at 0 of density / |cos y| = cos u sqrt(cos^2 u + t^2 G)
G, f = sp.symbols("G f", real=True)
dens = sp.cos(t*f/R)*sp.sqrt(sp.cos(t*f/R)**2 + t**2*G)
d2 = sp.diff(dens, t, 2).subs(t, 0)
kappa = 3 if PLANT else 2          # planted defect: Ric(nu,nu) = 3/R^2
check("T1 d2/dt2 density at t=0 equals |grad phi|^2 - 2 phi^2/R^2", sp.simplify(d2 - (G - kappa*f**2/R**2)))
d1 = sp.diff(dens, t).subs(t, 0)
check("T1 first derivative of density at t=0 vanishes (totally geodesic)", d1 + (G if PLANT else 0))   # planted: spurious G

# ---------------------------------------------------------------- T2
# Sector ODE near the vertex, r = pi/2 - y: -(1/sin r)(sin r f')' + nu^2/sin^2 r f.
# Indicial equation at r = 0 from leading terms -(1/r)(r f')' + nu^2/r^2 f.
r, nu = sp.symbols("r nu", positive=True)
s = sp.symbols("s", real=True)
ind = sp.simplify((-(1/r)*sp.diff(r*sp.diff(r**s, r), r) + nu**2/r**2*r**s) / r**(s - 2))
roots = sp.solve(sp.Eq(ind, 0), s)
expected = {nu, -nu} if not PLANT else {nu, -nu + 1}
check("T2 indicial exponents are +nu, -nu", 0 if set(roots) == expected else 1)
# Exact zero-energy solutions of the sector ODE: tan^{+-nu}(r/2).
for sg in (1, -1):
    fz = sp.tan(r/2)**(sg*nu)
    resz = -sp.diff(sp.sin(r)*sp.diff(fz, r), r)/sp.sin(r) + (nu**2 + (1 if PLANT else 0))/sp.sin(r)**2*fz
    resz = sp.simplify(sp.expand_trig(sp.simplify(resz/fz)))
    check(f"T2 tan^({'+' if sg > 0 else '-'}nu)(r/2) solves the sector ODE at lambda=0", resz)
# L^2(r dr) norm of r^{-nu} near 0: integral of r^{1-2nu}; finite iff nu < 1
nu_test = 1 if PLANT else sp.Rational(1, 2)          # planted: nu = 1 is not integrable
I = sp.integrate(r**(1 - 2*nu_test), (r, 0, 1))      # nu = 1/2: finite
J = sp.integrate(r**(1 - 2*(sp.Rational(9, 10) if PLANT else 1)),
                 (r, sp.Rational(1, 10**6), 1))      # nu = 1: log growth
check("T2 r^-nu square integrable for nu=1/2", 0 if I.is_finite else 1)
check("T2 r^-1 norm on [eps,1] is log(1/eps) (diverges, nu=1 limit-point)", sp.simplify(J - sp.log(10**6)))

# ---------------------------------------------------------------- T3
# Colatitude theta: eigenfunction f = sin^nu(theta) C_m^{(nu+1/2)}(cos theta),
# eigenvalue (nu+m)(nu+m+1).  Check the ODE for m = 0..5 with symbolic nu.
th = sp.symbols("theta", positive=True)
for m in range(6):
    fm = sp.sin(th)**nu*sp.gegenbauer(m, nu + sp.Rational(1, 2), sp.cos(th))
    lam = (nu + m)*(nu + m + 1) + (1 if PLANT else 0)
    res = -sp.diff(sp.sin(th)*sp.diff(fm, th), th)/sp.sin(th) + nu**2/sp.sin(th)**2*fm - lam*fm
    res = sp.simplify(sp.expand_trig(res/sp.sin(th)**nu))
    # evaluate at several points with nu exact-rational to avoid simplification gaps
    vals = [sp.nsimplify(res.subs({nu: sp.Rational(a, 7), th: sp.Rational(b, 5)}).evalf(40)) for a in (3, 7, 12) for b in (1, 4, 11)]
    check(f"T3 sector ODE solved by sin^nu C_{m}^(nu+1/2), lambda=(nu+{m})(nu+{m+1})",
          0 if all(abs(v) < 1e-30 for v in vals) and sp.simplify(res) == 0 else 1)
# orthogonality in L^2(sin theta d theta), nu = 3/7, m != m'
nv = sp.Rational(3, 7)
x = sp.symbols("x")
wexp = nv + (sp.Rational(1, 2) if PLANT else 0)   # planted: wrong weight exponent
for (m1, m2) in [(0, 2), (1, 3), (2, 4)]:
    integ = sp.integrate((1 - x**2)**wexp*sp.gegenbauer(m1, nv + sp.Rational(1, 2), x)
                         * sp.gegenbauer(m2, nv + sp.Rational(1, 2), x), (x, -1, 1))
    check(f"T3 orthogonality m={m1},{m2} in weight (1-x^2)^nu", sp.simplify(integ))

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "out_symbolic_planted.txt" if PLANT else "out_symbolic.txt"), "w") as fh:
    fh.write("\n".join(out) + "\n")
