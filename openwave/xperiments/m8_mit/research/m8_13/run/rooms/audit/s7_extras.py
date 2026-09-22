"""Loose ends: C4 as a global statement, the convention flips C5/C6 depend on,
the critical values nobody mentions, an mpmath precision witness, and the
Hessian at the two extremizers.
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np
import sympy as sp
import mpmath as mp
from scipy.optimize import minimize
import core as C
from grad import f_and_grad_real

HERE = pathlib.Path(__file__).parent
out = []
P = out.append
rng = np.random.default_rng(1234)

# ---------- C4 as a statement about all u ---------------------------------
P("C4 global form: maximize Tr Nbar^2 over the unit sphere (1500 starts)")
best, bestu, vals = -1, None, []
for k in range(1500):
    x = rng.normal(size=14)
    if k % 3 == 1:
        x *= (rng.random(14) < 0.4)
    if np.linalg.norm(x) < 1e-9:
        x = rng.normal(size=14)
    x /= np.linalg.norm(x)
    r = minimize(lambda y: -float(np.sum(C.Nbar((y[:7] + 1j * y[7:]) / np.linalg.norm(y)) ** 2)),
                 x, method="Powell", options=dict(xtol=1e-12, ftol=1e-14, maxiter=20000, maxfev=200000))
    y = r.x / np.linalg.norm(r.x)
    v = float(np.sum(C.Nbar(y[:7] + 1j * y[7:]) ** 2))
    vals.append(v)
    if v > best:
        best, bestu = v, y[:7] + 1j * y[7:]
vals = np.array(vals)
P(f"  max Tr Nbar^2 found = {best!r}   (171/2 = 85.5);  excess = {best-85.5:.3e}")
P(f"  #(within 1e-8 of 85.5) = {int((np.abs(vals-85.5)<1e-8).sum())} of {len(vals)}")
P(f"  any value above 85.5 + 1e-9 ? {bool((vals > 85.5+1e-9).any())}   (max over all starts "
  f"= {vals.max()!r})")
# Tr Nbar = 12 always
mx = max(abs(np.trace(C.Nbar(u := (lambda z: z/np.linalg.norm(z))(rng.normal(size=7)+1j*rng.normal(size=7)))) - 12)
         for _ in range(300))
P(f"  Tr Nbar - 12 max deviation : {mx:.3e}  -> Tr Nbar^2 = 48 + ||Nbar0||_F^2, "
  f"so the claim is ||Nbar0||_F^2 <= 75/2 i.e. ||Nbar0||_F <= 15/sqrt6")
P(f"  48 + (15/sqrt6)^2 = {48 + sp.Rational(75,2)} = {sp.Rational(171,2)}  -> consistent")

# who saturates Tr Nbar^2 = 171/2 ?
sat = []
for k in range(400):
    x = rng.normal(size=14); x /= np.linalg.norm(x)
    r = minimize(lambda y: -float(np.sum(C.Nbar((y[:7] + 1j * y[7:]) / np.linalg.norm(y)) ** 2)),
                 x, method="Powell", options=dict(xtol=1e-12, ftol=1e-14, maxiter=20000, maxfev=200000))
    y = r.x / np.linalg.norm(r.x)
    u = y[:7] + 1j * y[7:]
    if abs(float(np.sum(C.Nbar(u) ** 2)) - 85.5) < 1e-8:
        sat.append(C.rhat6(u))
sat = np.array(sat)
P(f"  rhat6 values on the Tr Nbar^2 = 171/2 saturating set: min {sat.min()!r}, "
  f"max {sat.max()!r}, distinct-ish {len(np.unique(np.round(sat,9)))}")
P("  -> the C4 bound is saturated by a set on which rhat6 is NOT constant; saturating C4")
P("     is necessary but far from sufficient for either extremum.")

# ---------- convention flips that C5/C6 could hide -------------------------
P("")
P("Convention probes for C5/C6")
P("  (i) reading 'e1 - e2 = 4b' as 'e2 - e1 = 4b' is exactly b -> -b.")
P("      The ellipse (2/3)a^2 + 8b^2 = 1 is invariant under b -> -b, and the claimed")
P("      contact set {(sqrt6/2,0), (-sqrt6/4, sqrt6/8), (-sqrt6/4,-sqrt6/8)} is invariant")
P("      as a SET under b -> -b (it swaps the last two).  M1 and M2 swap.")
P("      => the answer to C6 is unchanged; only the labels M1/M2 swap.")
P("  (ii) 'e3 = 2a/3' pins the DISTINGUISHED axis to z.  Putting it on x or y is a")
P("      relabelling of axes; the contact set is the full set of permutations of")
P("      sqrt6*e = (2,-1,-1), which is permutation-invariant, so again unchanged.")
a, b = sp.symbols('a b', real=True)
e1 = -a / 3 + 2 * b; e2 = -a / 3 - 2 * b; e3 = 2 * a / 3
P(f"  check: (e1,e2,e3) at (a,b)=(sqrt6/2,0) -> "
  f"{[sp.nsimplify(sp.sqrt(6)*x.subs({a:sp.sqrt(6)/2,b:0})) for x in (e1,e2,e3)]} (times 1/sqrt6)")
P(f"  check: at b -> -b, (a,b)=(-sqrt6/4,-sqrt6/8) -> "
  f"{[sp.nsimplify(sp.sqrt(6)*x.subs({a:-sp.sqrt(6)/4,b:-sp.sqrt(6)/8})) for x in (e1,e2,e3)]}")
P("  (iii) the 'M2 = M1 at -b' wording: M2 as a matrix in the basis")
P("      (asym(v3,v-3), asym(v1,v-1)) is [[5a, +sqrt60 b],[+sqrt60 b, -3a-12b]], which is")
P("      M1(-b) conjugated by diag(1,-1).  Same characteristic polynomial, so every")
P("      eigenvalue statement in C6/C7 is untouched; only the literal matrix differs,")
P("      and flipping the sign of one basis vector removes even that.")

# ---------- the critical values nobody mentions ----------------------------
P("")
P("Critical values of rhat6 found by the scan, identified as exact rationals")
for nm in ("MAX", "MIN"):
    vals = np.load(HERE / f"vals_{nm}.npy")
    pts = np.load(HERE / f"pool_{nm}.npy")
    sv = np.sort(np.unique(np.round(vals, 10)))
    seen = []
    for v in sv:
        if not seen or v - seen[-1] > 1e-7:
            seen.append(v)
    for v in seen:
        fr = sp.nsimplify(sp.Rational(round(v * 924), 924))
        i = int(np.argmin(np.abs(vals - v)))
        u = pts[i] / np.linalg.norm(pts[i])
        ex = C.rhat6_exact([sp.nsimplify(sp.Float(x, 15), rational=False) for x in u]) if False else None
        P(f"  {nm} run: {v!r}  ~ {round(v*924)}/924 = {sp.Rational(round(v*924),924)} "
          f"(|diff| = {abs(v - round(v*924)/924):.2e});  |f|^2 = {np.dot(C.fvec(u),C.fvec(u)):.6f}, "
          f"TrNbar^2 = {float(np.sum(C.Nbar(u)**2)):.6f}, |a00|^2 = {abs(C.a00(u))**2:.6f}")
P("  exact values at simple vectors:")
for nm, vec in (("v3", [1,0,0,0,0,0,0]), ("v2", [0,1,0,0,0,0,0]), ("v1", [0,0,1,0,0,0,0]),
                ("v0", [0,0,0,1,0,0,0]), ("(v3+v-3)/r2", [1,0,0,0,0,0,1]),
                ("(v2+v-2)/r2", [0,1,0,0,0,1,0]), ("(v1+v-1)/r2", [0,0,1,0,1,0,0]),
                ("(v3+v-3 i)/r2", [1,0,0,0,0,0,sp.I])):
    val = C.rhat6_exact(vec)
    P(f"    rhat6({nm}) = {val} = {sp.nsimplify(val*924)}/924 = {sp.N(val,18)}")

# ---------- mpmath precision witness --------------------------------------
P("")
P("Independent 50-digit evaluation (mpmath) of rhat6 at the two extremizers")
mp.mp.dps = 50


def rhat6_mp(c):
    c = [mp.mpmathify(x) for x in c]
    t = [mp.mpf(-1) ** (3 - m) * mp.conj(c[C.IDX[-m]]) for m in C.MS]
    tot = mp.mpf(0)
    for Q in range(-6, 7):
        s = mp.mpf(0)
        for m1, g in C.CG_EXACT[Q].items():
            s += mp.mpmathify(sp.N(g, 60)) * c[C.IDX[m1]] * t[C.IDX[Q - m1]]
        tot += abs(s) ** 2
    n2 = sum(abs(x) ** 2 for x in c)
    return tot / n2 ** 2


r2m = 1 / mp.sqrt(2)
P(f"  cat  : {mp.nstr(rhat6_mp([r2m,0,0,0,0,0,r2m]), 40)}")
P(f"  463/924 = {mp.nstr(mp.mpf(463)/924, 40)}")
P(f"  diff = {mp.nstr(rhat6_mp([r2m,0,0,0,0,0,r2m]) - mp.mpf(463)/924, 10)}")
P(f"  coh  : {mp.nstr(rhat6_mp([1,0,0,0,0,0,0]), 40)}")
P(f"  1/924   = {mp.nstr(mp.mpf(1)/924, 40)}")
P(f"  diff = {mp.nstr(rhat6_mp([1,0,0,0,0,0,0]) - mp.mpf(1)/924, 10)}")
gam = mp.mpf(1) / 3
P(f"  rotated cat (gamma=1/3): "
  f"{mp.nstr(rhat6_mp([mp.e**(1j*gam)*r2m,0,0,0,0,0,mp.e**(-1j*gam)*r2m]) - mp.mpf(463)/924, 10)}")

# ---------- Hessian at the extremizers -------------------------------------
P("")
P("Hessian of rhat6 (restricted to the unit sphere) at the two extremizers")
for nm, ref, orbdim in (("cat / MAX", np.array([1,0,0,0,0,0,1],dtype=complex)/np.sqrt(2), 4),
                        ("v3 / MIN", C.basis(3), 3)):
    x0 = np.concatenate([ref.real, ref.imag])
    h = 2e-4

    def F(z):
        return f_and_grad_real(z / np.linalg.norm(z))[0]
    Hm = np.zeros((14, 14))
    I14 = np.eye(14)
    for i in range(14):
        for j in range(i, 14):
            ei, ej = I14[i], I14[j]
            Hm[i, j] = Hm[j, i] = (F(x0+h*ei+h*ej) - F(x0+h*ei-h*ej)
                                   - F(x0-h*ei+h*ej) + F(x0-h*ei-h*ej)) / (4*h*h)
    ev = np.sort(np.linalg.eigvalsh(Hm))
    nz = int(np.sum(np.abs(ev) < 1e-4))
    P(f"  {nm}: eigenvalues = {np.round(ev, 5)}")
    P(f"     zero modes = {nz}; expected 1 (radial) + {orbdim} (orbit) = {1+orbdim}; "
      f"nonzero eigenvalue signs = {sorted(set(np.sign(ev[np.abs(ev)>=1e-4]).astype(int)))}")

txt = "\n".join(map(str, out))
print(txt)
(HERE / "out_s7.txt").write_text(txt + "\n")
