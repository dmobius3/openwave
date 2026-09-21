"""S0 search aid, second part: do the literature's identities already settle the global extrema of r6 = ||rho_6||^2?

Checks, each armed:
  (E) r6 = C + x|f|^2 + y|a00|^2 + z TrN^2, Kawaguchi-Ueda's invariants (review Eq. 167, Table 7), with
      (C, x, y, z) = (-5/231, -1/22, 7/11, 1/198) as translated from the pair channels;
  (F) the bounds the argument needs, |a00|^2 <= 1/7 and TrN^2 <= 171/2 (stated in the review, not proved there),
      probed by multi-start local maximization, and a direct multi-start maximization and minimization of r6.
Not an S computation: no Hessian and no index is computed or reported.
"""
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import minimize
from translate_rho6 import rho_sq, vec, CG, ms, idx, j

fails = []


def gate(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        fails.append(name)


jz = np.diag([float(m) for m in ms]).astype(complex)
jp = np.zeros((7, 7), complex)
for m in ms[1:]:
    jp[idx[m + 1], idx[m]] = np.sqrt(j * (j + 1) - m * (m + 1))
jx, jy = (jp + jp.conj().T) / 2, (jp - jp.conj().T) / 2j
F3 = (jx, jy, jz)


def inv(u):
    u = u / np.linalg.norm(u)
    f2 = sum(abs(np.vdot(u, A @ u)) ** 2 for A in F3)
    a00 = abs(sum(CG[0, m, -m] * u[idx[m]] * u[idx[-m]] for m in ms)) ** 2
    N = np.array([[np.vdot(u, (A @ B + B @ A) @ u).real / 2 for B in F3] for A in F3])
    return f2, a00, float(np.sum(N * N)), float(np.trace(N))


rng = np.random.default_rng(7)
rand = lambda: rng.normal(size=7) + 1j * rng.normal(size=7)
coef = (Fr(-5, 231), Fr(-1, 22), Fr(7, 11), Fr(1, 198))
model = lambda u, c=coef: float(c[0]) + sum(float(ci) * vi for ci, vi in zip(c[1:], inv(u)[:3]))
S = [rand() for _ in range(60)]
print('(E) the Kawaguchi-Ueda form of r6')
gate('  TrN = f(f+1) = 12 on 60 random states (the review\'s constraint)', max(abs(inv(u)[3] - 12) for u in S) < 1e-12)
res = max(abs(model(u) - rho_sq(u)) for u in S)
gate(f'  r6 = -5/231 - |f|^2/22 + (7/11)|a00|^2 + TrN^2/198 on 60 random states (max residual {res:.1e})', res < 1e-12)
gate('  arm: z = 1/200 in place of 1/198 fails', max(abs(model(u, coef[:3] + (Fr(1, 200),)) - rho_sq(u)) for u in S) > 1e-4)
for name, u in [('coherent', vec(p3=1)), ('hexagon', vec(p3=1, m3=1)), ('octahedron', vec(p2=1, m2=1)), ('zonal', vec(p0=1))]:
    f2, a00, n2, _ = inv(u)
    print(f'  {name:10s} |f|^2 = {f2:.6f}  |a00|^2 = {a00:.6f} (1/7 = {1/7:.6f})  TrN^2 = {n2:.6f}  r6 = {rho_sq(u):.6f}')


def maximize(fun, starts=300):
    best = None
    for _ in range(starts):
        x0 = rng.normal(size=14)
        r = minimize(lambda x: -fun(x[:7] + 1j * x[7:]), x0, method='BFGS', options={'gtol': 1e-11})
        if best is None or -r.fun > best[0]:
            best = (-r.fun, (r.x[:7] + 1j * r.x[7:]) / np.linalg.norm(r.x))
    return best


print('(F) the bounds, and r6 itself, by multi-start local optimization (300 starts each)')
b_a00 = maximize(lambda u: inv(u)[1])
gate(f'  max |a00|^2 found = {b_a00[0]:.12f}, not above 1/7', b_a00[0] <= 1 / 7 + 1e-10)
b_n2 = maximize(lambda u: inv(u)[2])
gate(f'  max TrN^2 found = {b_n2[0]:.12f}, not above 171/2', b_n2[0] <= 85.5 + 1e-9)
gate('  arm: the probe reaches the bound (it can find a maximum), within 1e-8', abs(b_n2[0] - 85.5) < 1e-8 and abs(b_a00[0] - 1 / 7) < 1e-10)
b_max = maximize(rho_sq)
f2, a00, n2, _ = inv(b_max[1])
gate(f'  max r6 found = {b_max[0]:.12f} = 463/924 ({463/924:.12f}); maximizer has |f|^2={f2:.1e}, |a00|^2={a00:.6f}, TrN^2={n2:.6f}',
     abs(b_max[0] - 463 / 924) < 1e-10 and f2 < 1e-8 and abs(a00 - 1 / 7) < 1e-8 and abs(n2 - 85.5) < 1e-7)
b_min = maximize(lambda u: -rho_sq(u))
f2, a00, n2, _ = inv(b_min[1])
gate(f'  min r6 found = {-b_min[0]:.12f} = 1/924; minimizer has |f|^2 = {f2:.8f} (coherent: 9)', abs(-b_min[0] - 1 / 924) < 1e-10 and abs(f2 - 9) < 1e-6)
print(f'{len(fails)} FAILED: {fails}' if fails else 'ALL PASS')
