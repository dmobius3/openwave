"""T1: symbolic check of the second variation of area for M(W) in RP^3(R).

Varied map: F = cos(t*phi/R) X + R sin(t*phi/R) e4  (normal geodesic of S^3(R)).
Checks (each can fail; see --plant):
  C1.1  Gram matrix of F equals cos^2(t phi/R) g + t^2 dphi (x) dphi, exactly, both branches.
  C1.2  t^2 coefficient of sqrt(det g_t)/sqrt(det g) equals (|grad phi|^2 - 2 phi^2/R^2)/2.
  C1.3  seam: F(pi R, -w) = -F(0, w) when phi(pi R, -w) = -phi(0, w).
Plants: straight (F = X + t phi e4), seamsign (phi(pi R,-w) = +phi(0,w)).
"""
import sys
import sympy as sp

plant = sys.argv[sys.argv.index("--plant") + 1] if "--plant" in sys.argv else None

y, w, t, a = sp.symbols("y w t a", real=True)
R = sp.symbols("R", positive=True)
phi = sp.Function("phi")(y, w)
e4 = sp.Matrix([0, 0, 0, 1])


def X(yy, ww):
    return R * sp.Matrix([sp.cos(yy / R) * sp.cos(ww / R), sp.cos(yy / R) * sp.sin(ww / R), sp.sin(yy / R), 0])


def varied(Xv, f):
    if plant == "straight":
        return Xv + t * f * e4
    return sp.cos(t * f / R) * Xv + R * sp.sin(t * f / R) * e4


ok_all = True


def report(name, ok):
    global ok_all
    ok_all &= ok
    print(f"{name}: {'PASS' if ok else 'FAIL'}")


for branch, sgn in (("y<=piR/2", 1), ("y>piR/2", -1)):
    Xb = X(y, sgn * w)
    F = varied(Xb, phi)
    Fy, Fw = F.diff(y), F.diff(w)
    Gt = sp.Matrix([[Fy.dot(Fy), Fy.dot(Fw)], [Fw.dot(Fy), Fw.dot(Fw)]])
    Xy, Xw = Xb.diff(y), Xb.diff(w)
    g = sp.Matrix([[Xy.dot(Xy), Xy.dot(Xw)], [Xw.dot(Xy), Xw.dot(Xw)]])
    g_expected = sp.Matrix([[1, 0], [0, sp.cos(y / R) ** 2]])
    dphi = sp.Matrix([phi.diff(y), phi.diff(w)])
    c = sp.cos(t * phi / R)
    target = c**2 * g + t**2 * dphi * dphi.T
    diff = sp.simplify(sp.trigsimp(Gt - target))
    report(f"C1.1 [{branch}] g = ds^2 of sec.2 and g_t = cos^2 g + t^2 dphi dphi",
           diff == sp.zeros(2, 2) and sp.simplify(g - g_expected) == sp.zeros(2, 2))
    # area-density ratio sqrt(det g_t / det g), expanded to t^2
    ratio = sp.sqrt(sp.expand(Gt.det()) / g_expected.det())
    ser = sp.series(ratio, t, 0, 3).removeO()
    c2 = sp.simplify(ser.coeff(t, 2))
    grad2 = phi.diff(y) ** 2 + phi.diff(w) ** 2 / sp.cos(y / R) ** 2
    expected = (grad2 - 2 * phi**2 / R**2) / 2
    report(f"C1.2 [{branch}] t^2 coeff of dA_t/dA = (|grad phi|^2 - 2 phi^2/R^2)/2",
           sp.simplify(c2 - expected) == 0 and sp.simplify(ser.coeff(t, 1)) == 0
           and sp.simplify(ser.coeff(t, 0) - 1) == 0)

# seam: value a at (0,w) on branch 1, value -a at (pi R, -w) on branch 2
seam_val = a if plant == "seamsign" else -a
F1 = varied(X(0, w), a)
F2 = varied(X(sp.pi * R, -(-w)), seam_val)  # branch 2 parametrizes (y, w') by X(y, -w'); w' = -w
report("C1.3 seam: F(piR,-w) = -F(0,w) (same point of RP^3)",
       sp.simplify(F1 + F2) == sp.zeros(4, 1))

print("T1_ALL", "PASS" if ok_all else "FAIL")
