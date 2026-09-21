"""S0 lemma check: at spin 3 the hexagon (NOON) orbit is the unique global maximum of r6 = ||rho_6||^2.

Proof chain, each link checked exactly (sympy) or armed:
  (L1) r6 = -5/231 - |f|^2/22 + (7/11)|a00|^2 + TrN^2/198 exactly: the four functions span the 4-dim space of invariant
       quartics (independence determinant), and the four weight states fix the coefficients (exact rationals).
  (L2) TrN^2 = 48 + ||Q||^2 and max_rho ||Q|| = max over unit traceless e of lambda_max(e1 fx^2 + e2 fy^2 + e3 fz^2).
       Writing that operator as a(fz^2 - 4) + b(f+^2 + f-^2), with (2/3)a^2 + 8b^2 = 1, it splits into three 2x2 blocks
       and a zero; each block's top eigenvalue is at most 15/sqrt(6), by an exact square-completion (M3) and a sum of
       squares 36 b^2 (a + 2b)^2 (M1, M2). Equality only at e a permutation of (2,-1,-1)/sqrt(6).
  (L3) equality: support on m = +-3 along an axis; |f| = 0 then forces |alpha| = |beta|: the hexagon orbit.
Not an S computation: no Hessian, no index.
"""
import numpy as np
import sympy as sp
from translate_rho6 import rho_sq, vec, ms, idx, j
from check_extrema import inv, F3, jz, jp

fails = []


def gate(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        fails.append(name)


print('(L1) the invariant form, exactly')
R = sp.Rational
# weight states |3,m>: |f|^2 = m^2; a00 = 0 unless m = 0 (1/7); N = diag((12-m^2)/2, (12-m^2)/2, m^2)
rows, r6w = [], {3: R(1, 924), 2: R(36, 924), 1: R(225, 924), 0: R(400, 924)}
for m in (3, 2, 1, 0):
    trn2 = 2 * (R(12 - m * m, 2)) ** 2 + m ** 4
    rows.append([1, m * m, R(1, 7) if m == 0 else 0, trn2])
A = sp.Matrix(rows)
gate(f'  the four invariants are independent at the weight states (det = {A.det()})', A.det() != 0)
coef = A.solve(sp.Matrix([r6w[m] for m in (3, 2, 1, 0)]))
gate(f'  coefficients forced by the weight states: {list(coef)}', list(coef) == [R(-5, 231), R(-1, 22), R(7, 11), R(1, 198)])
for m, row in zip((3, 2, 1, 0), rows):                      # the weight-state inputs come from the aid, not from the paper alone
    f2, a00, n2, _ = inv(vec(**{('p' if m >= 0 else 'm') + str(abs(m)): 1}))
    gate(f'  weight state m={m}: aid gives |f|^2, |a00|^2, TrN^2 = {f2:.6f}, {a00:.6f}, {n2:.6f}; r6 = {rho_sq(vec(**{"p" + str(m): 1})) * 924:.6f}/924',
         abs(f2 - float(row[1])) < 1e-12 and abs(a00 - float(row[2])) < 1e-12 and abs(n2 - float(row[3])) < 1e-12
         and abs(rho_sq(vec(**{'p' + str(m): 1})) - float(r6w[m])) < 1e-12)

print('(L2) the eigenvalue lemma, exactly')
a, b, e1, e2, e3 = sp.symbols('a b e1 e2 e3', real=True)
K = 15 / sp.sqrt(6)
# the operator identity on V3, checked as 7x7 matrices with exact entries
Jz = sp.diag(*[m for m in ms])
Jp = sp.zeros(7, 7)
for m in ms[1:]:
    Jp[idx[m + 1], idx[m]] = sp.sqrt(j * (j + 1) - m * (m + 1))
Jm = Jp.T
Jx, Jy = (Jp + Jm) / 2, (Jp - Jm) / (2 * sp.I)
lhs = e1 * Jx ** 2 + e2 * Jy ** 2 + e3 * Jz ** 2
rhs = (sp.Rational(3, 2) * e3) * (Jz ** 2 - 4 * sp.eye(7)) + ((e1 - e2) / 4) * (Jp ** 2 + Jm ** 2)
gate('  e1 fx^2 + e2 fy^2 + e3 fz^2 = a(fz^2 - 4) + b(f+^2 + f-^2) on the traceless plane, a = 3e3/2, b = (e1-e2)/4',
     sp.simplify((lhs - rhs).subs(e1, -e2 - e3)) == sp.zeros(7, 7))
gate('  constraint: e unit traceless  <=>  (2/3)a^2 + 8b^2 = 1',
     sp.expand((e1 ** 2 + e2 ** 2 + e3 ** 2).subs(e1, -e2 - e3) - (sp.Rational(2, 3) * (sp.Rational(3, 2) * e3) ** 2
                                                                    + 8 * ((-e2 - e3 - e2) / 4) ** 2)) == 0)
Bm = rhs.subs({e3: 2 * a / 3}).subs({e1 - e2: 4 * b})
Bm = sp.Rational(3, 2) * (2 * a / 3) * (Jz ** 2 - 4 * sp.eye(7)) + b * (Jp ** 2 + Jm ** 2)
M1 = sp.Matrix([[5 * a, sp.sqrt(60) * b], [sp.sqrt(60) * b, -3 * a + 12 * b]])
M2 = sp.Matrix([[5 * a, sp.sqrt(60) * b], [sp.sqrt(60) * b, -3 * a - 12 * b]])
M3 = sp.Matrix([[0, sp.sqrt(240) * b], [sp.sqrt(240) * b, -4 * a]])
lam = sp.Symbol('lam')
cp = sp.factor(sp.expand((Bm - lam * sp.eye(7)).det()))
blocks = sp.factor(sp.expand(M1.charpoly(lam).as_expr() * M2.charpoly(lam).as_expr() * M3.charpoly(lam).as_expr() * (-lam)))
gate('  spectrum of the 7x7 operator = spectra of M1, M2, M3 and a zero (characteristic polynomials agree)', sp.expand(cp - blocks) == 0)
lam1 = a + 6 * b + sp.sqrt(16 * a ** 2 - 48 * a * b + 96 * b ** 2)
gate('  top eigenvalue of M1 is a + 6b + sqrt(16a^2 - 48ab + 96b^2)', sp.simplify(M1.charpoly(lam).as_expr().subs(lam, lam1)) == 0)
# M3: with 240 b^2 = 30 - 20 a^2 the bound is 20 (a + sqrt(6)/4)^2 >= 0
gate('  M3: (K + 2a)^2 - (30 - 16a^2) = 20 (a + sqrt(6)/4)^2', sp.simplify(sp.expand((K + 2 * a) ** 2 - (30 - 16 * a ** 2) - 20 * (a + sp.sqrt(6) / 4) ** 2)) == 0)
# M1: after homogenizing with N^2 = (2/3)a^2 + 8b^2, the squared inequality has slack 36 b^2 (a + 2b)^2
slack = sp.expand((a ** 2 + 6 * a * b + 24 * b ** 2) ** 2 - sp.Rational(3, 2) * (a + 6 * b) ** 2 * (sp.Rational(2, 3) * a ** 2 + 8 * b ** 2))
gate('  M1: (a^2 + 6ab + 24b^2)^2 - (3/2)(a + 6b)^2 N^2 = 36 b^2 (a + 2b)^2', sp.expand(slack - 36 * b ** 2 * (a + 2 * b) ** 2) == 0)
gate('  arm: the same identity with 24b^2 -> 23b^2 fails',
     sp.expand((a ** 2 + 6 * a * b + 23 * b ** 2) ** 2 - sp.Rational(3, 2) * (a + 6 * b) ** 2 * (sp.Rational(2, 3) * a ** 2 + 8 * b ** 2) - 36 * b ** 2 * (a + 2 * b) ** 2) != 0)
# the homogenized step itself: K^2 N^2 = 25a^2 + 300b^2
gate('  K^2 N^2 = 25a^2 + 300b^2', sp.expand(K ** 2 * (sp.Rational(2, 3) * a ** 2 + 8 * b ** 2) - 25 * a ** 2 - 300 * b ** 2) == 0)
# the side conditions the case analysis uses
gate('  M2 has the spectrum of M1 at b -> -b (equal characteristic polynomials), and the constraint is even in b',
     sp.expand(M2.charpoly(lam).as_expr() - M1.subs(b, -b).charpoly(lam).as_expr()) == 0)
gate('  arm: M2 and M1 themselves have different spectra (the b -> -b step is needed)',
     sp.expand(M2.charpoly(lam).as_expr() - M1.charpoly(lam).as_expr()) != 0)
gate('  K - a - 6b > 0 on the ellipse: max(a + 6b) = sqrt((3/2) + 36/8) = sqrt(6) < 15/sqrt(6)', sp.Rational(3, 2) + sp.Rational(36, 8) == 6 and sp.sqrt(6) < K)
gate('  K + 2a > 0 on the ellipse: K - 2 sqrt(3/2) = 9/sqrt(6)', sp.simplify(K - 2 * sp.sqrt(sp.Rational(3, 2)) - 9 / sp.sqrt(6)) == 0)
gate('  a^2 + 6ab + 24b^2 = (a + 3b)^2 + 15b^2 >= 0', sp.expand(a ** 2 + 6 * a * b + 24 * b ** 2 - (a + 3 * b) ** 2 - 15 * b ** 2) == 0)
# numeric sweep of the whole circle, an independent control on the exact chain
Jx_, Jy_, Jz_ = F3
Kf = 15 / np.sqrt(6)


def top_eig(t, extra=None):
    e = np.sqrt(2 / 3) * np.array([np.cos(t), np.cos(t - 2 * np.pi / 3), np.cos(t + 2 * np.pi / 3)])
    Bn = e[0] * Jx_ @ Jx_ + e[1] * Jy_ @ Jy_ + e[2] * Jz_ @ Jz_
    return np.linalg.eigvalsh(Bn if extra is None else Bn + extra).max()


th = np.linspace(0, 2 * np.pi, 20001)
top = np.array([top_eig(t) for t in th])
axes = np.array([0, 2 * np.pi / 3, 4 * np.pi / 3])
dist = np.min(np.abs(((th[:, None] - axes[None, :]) + np.pi) % (2 * np.pi) - np.pi), axis=1)
gate(f'  sweep of 20001 points: max lambda_max = {top.max():.12f} <= 15/sqrt(6) = {Kf:.12f}', top.max() <= Kf + 1e-12)
gate('  equality at the three axis directions exactly', all(abs(top_eig(t) - Kf) < 1e-12 for t in axes))
gate(f'  strictly below away from them: max over points > 0.01 rad from an axis = {top[dist > 0.01].max():.9f}', top[dist > 0.01].max() < Kf - 1e-5)
pert = 0.02 * (Jx_ @ Jx_)
gate('  arm: the same sweep on a perturbed operator (+0.02 fx^2) exceeds 15/sqrt(6)', max(top_eig(t, pert) for t in th[::50]) > Kf + 1e-3)

print('(L3) the equality case, on the computed objects')
hexa = vec(p3=np.sqrt(0.5), m3=np.sqrt(0.5))
coh = vec(p3=1)
f2h, a00h, n2h, _ = inv(hexa)
gate(f'  hexagon attains all three extremes: |f|^2 = {f2h:.1e}, |a00|^2 = {a00h:.12f}, TrN^2 = {n2h:.12f}', f2h < 1e-14 and abs(a00h - 1 / 7) < 1e-14 and abs(n2h - 85.5) < 1e-12)
fam = [vec(p3=np.cos(s), m3=np.sin(s)) for s in np.linspace(0.05, np.pi / 2 - 0.05, 9)]
gate('  every u on span{v3, v-3} has TrN^2 = 171/2 (the equality set of the TrN^2 bound)', all(abs(inv(u)[2] - 85.5) < 1e-12 for u in fam))
gate('  on that span, r6 < 463/924 unless |alpha| = |beta|', all(rho_sq(u) < 463 / 924 - 1e-6 for u in fam if abs(abs(u[idx[3]]) - abs(u[idx[-3]])) > 1e-3))
gate(f'  the value is 463/924: r6(hexagon) = {rho_sq(hexa) * 924:.12f}/924', abs(rho_sq(hexa) - 463 / 924) < 1e-14)

print('(L4) arithmetic the note quotes')
gate(f'  -5/231 + 1/11 + 171/396 = {R(-5, 231) + R(1, 11) + R(171, 396)} = 463/924', R(-5, 231) + R(1, 11) + R(171, 396) == R(463, 924))
gate(f'  (15/sqrt(6))^2 = {sp.nsimplify(K ** 2)} = 225/6, and 48 + 225/6 = {48 + R(225, 6)} = 171/2', sp.simplify(K ** 2 - R(225, 6)) == 0 and 48 + R(225, 6) == R(171, 2))
gate('  the M3 equality point is a = -√6/4, where 240 b^2 = 30 - 20 a^2 gives b^2 = 3/32', sp.simplify(30 - 20 * (sp.sqrt(6) / 4) ** 2 - 240 * R(3, 32)) == 0)
print(f'{len(fails)} FAILED: {fails}' if fails else 'ALL PASS')
