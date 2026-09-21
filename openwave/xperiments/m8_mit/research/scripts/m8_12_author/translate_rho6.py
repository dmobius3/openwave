"""S0 search aid: place M8.11's reduced quartic r6 = ||rho_6(u)||^2 in the literature's coordinates.

Not an S computation: no Hessian, no index. It checks three translations the search needs, each armed:
  (A) my rho_K implementation reproduces the paper's recorded values (Section 5.8), and fails with K=5 or u(x)u;
  (B) ||rho_6||^2 = sum_F c_F ||[u(x)u]_F||^2 (the spin-3 BEC pair channels), c_F rational;
  (C) ||rho_6||^2 = 13 * int Q_u(n) Q_u(-n) dmu(n) (antipodal Husimi overlap); the plain second moment must fail.
Then maps r6 into Diener-Ho's (alpha, beta, gamma) for both signs of g.
"""
import itertools
from fractions import Fraction as Fr
import numpy as np
from sympy.physics.wigner import clebsch_gordan as sCG

j = 3
ms = list(range(j, -j - 1, -1))            # basis order v_3 .. v_-3
idx = {m: i for i, m in enumerate(ms)}
CG = {}
for K in range(0, 2 * j + 1):
    for m1, m2 in itertools.product(ms, ms):
        if abs(m1 + m2) <= K:
            CG[K, m1, m2] = float(sCG(j, j, K, m1, m2, m1 + m2))
fails = []


def gate(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        fails.append(name)


def theta(u):                               # time reversal, (Theta u)_m = (-1)^(j-m) conj(u_-m)
    return np.array([(-1) ** (j - m) * np.conj(u[idx[-m]]) for m in ms])


def couple(a, b, K):
    return np.array([sum(CG[K, m1, Q - m1] * a[idx[m1]] * b[idx[Q - m1]]
                         for m1 in ms if abs(Q - m1) <= j and (K, m1, Q - m1) in CG) for Q in range(-K, K + 1)])


def rho_sq(u, K=6, tr=True):
    u = u / np.linalg.norm(u)
    v = couple(u, theta(u) if tr else u, K)
    return float(np.vdot(v, v).real)


def pair_sq(u, F):
    u = u / np.linalg.norm(u)
    v = couple(u, u, F)
    return float(np.vdot(v, v).real)


def vec(**c):
    u = np.zeros(7, complex)
    for k, x in c.items():
        u[idx[int(k.replace('m', '-').replace('p', ''))]] = x
    return u


s2 = np.sqrt(0.5)
t_pyr = np.arcsin(np.sqrt(12 / 25))
REC = [('coherent v3, 1/924', vec(p3=1), 1 / 924), ('zonal v0, 400/924', vec(p0=1), 400 / 924),
       ('octahedron, 288/924', vec(p2=s2, m2=s2), 288 / 924), ('hexagon, 463/924', vec(p3=s2, m3=s2), 463 / 924),
       ('pyramid s=12/25, 9/35', vec(p2=np.cos(t_pyr), m3=np.sin(t_pyr)), 9 / 35),
       ('prism z=sqrt(23/10), 200/903', vec(p3=1, p0=np.sqrt(2.3), m3=1), 200 / 903),
       ('octahedron in D3 chart z=i*sqrt(5/2), 24/77', vec(p3=1, p0=1j * np.sqrt(2.5), m3=1), 24 / 77)]
print('(A) the recorded values, Section 5.8')
for name, u, val in REC:
    gate(f'  {name}: ||rho_6||^2 = {val:.10f}', abs(rho_sq(u) - val) < 1e-12)
gate('  arm: K=5 misses the table', not all(abs(rho_sq(u, K=5) - v) < 1e-9 for _, u, v in REC))
gate('  arm: u(x)u in place of u(x)Theta(u) misses the table', not all(abs(rho_sq(u, tr=False) - v) < 1e-9 for _, u, v in REC))
for s in (0.1, 0.37, 0.8):
    t = np.arcsin(np.sqrt(s))
    gate(f'  pyramid line formula at s={s}', abs(rho_sq(vec(p2=np.cos(t), m3=np.sin(t))) - (-125 / 132 * s * s + 10 / 11 * s + 3 / 77)) < 1e-12)

print('(B) pair-channel decomposition')
rng = np.random.default_rng(20260919)
rand = lambda: rng.normal(size=7) + 1j * rng.normal(size=7)
train = [rand() for _ in range(40)]
A = np.array([[pair_sq(u, F) for F in (0, 2, 4, 6)] for u in train])
b = np.array([rho_sq(u) for u in train])
c = np.linalg.lstsq(A, b, rcond=None)[0]
cF = [Fr(x).limit_denominator(100000) for x in c]
print('  c_F (F = 0, 2, 4, 6):', [str(x) for x in cF])
test = [rand() for _ in range(40)]
res = max(abs(sum(float(cF[i]) * pair_sq(u, F) for i, F in enumerate((0, 2, 4, 6))) - rho_sq(u)) for u in test)
gate(f'  rational c_F reproduce ||rho_6||^2 on 40 fresh states (max residual {res:.1e})', res < 1e-12)
bad = list(cF); bad[1] += Fr(1, 1000)
res_bad = max(abs(sum(float(bad[i]) * pair_sq(u, F) for i, F in enumerate((0, 2, 4, 6))) - rho_sq(u)) for u in test)
gate('  arm: a perturbed c_2 fails', res_bad > 1e-6)

print('(C) Husimi reading')
xg, wg = np.polynomial.legendre.leggauss(40)            # cos(theta) nodes
ph = np.linspace(0, 2 * np.pi, 41)[:-1]
from math import comb


def coh(ct, p):                                         # spin-3 coherent state at (theta, phi)
    c2, s2_ = np.sqrt((1 + ct) / 2), np.sqrt((1 - ct) / 2)
    return np.array([np.sqrt(comb(2 * j, j - m)) * c2 ** (j + m) * s2_ ** (j - m) * np.exp(1j * (j - m) * p) for m in ms])


def husimi_moments(u):
    u = u / np.linalg.norm(u)
    anti = plain = 0.0
    for x, w in zip(xg, wg):
        for p in ph:
            q = abs(np.vdot(coh(x, p), u)) ** 2
            qa = abs(np.vdot(coh(-x, p + np.pi), u)) ** 2
            anti += w * q * qa
            plain += w * q * q
    norm = 2 * len(ph)                                  # int dmu = 1: sum w = 2, times 40 phi points
    return 13 * anti / norm, 13 * plain / norm


h_coh = husimi_moments(vec(p3=1))
gate(f'  control: 13*int Q^2 = 1 at the coherent state ({h_coh[1]:.12f}, Lieb bound attained)', abs(h_coh[1] - 1) < 1e-10)
errs = [abs(husimi_moments(u)[0] - rho_sq(u)) for u in test[:8]]
gate(f'  ||rho_6||^2 = 13*int Q(n)Q(-n) on 8 fresh states (max error {max(errs):.1e})', max(errs) < 1e-10)
gate('  arm: the plain second moment 13*int Q^2 fails off time reversal', min(abs(husimi_moments(u)[1] - rho_sq(u)) for u in test[:8]) > 1e-3)
for name, u, val in REC[:4]:
    a, p = husimi_moments(u)
    print(f'  {name}: 13*int Q(n)Q(-n) = {a:.10f}, 13*int Q^2 = {p:.10f}')

print('(D) Diener-Ho coordinates: E ~ alpha|Theta|^2 + beta sum|B_M|^2 + gamma <S>^2 + C, a_F - a_6 = sign(g) w6 (c_F - c_6)')
for sg in (+1, -1):
    d0, d2, d4 = (sg * (cF[i] - cF[3]) for i in range(3))   # a_F - a_6 up to the positive factor w6
    alpha = (d0 - Fr(21, 11) * d4) / 7
    beta = (d2 - Fr(18, 11) * d4) / 7
    gamma = -d4 / 11
    print(f'  g {"> 0" if sg > 0 else "< 0"}: (a0-a6, a2-a6, a4-a6) ~ ({d0}, {d2}, {d4}); alpha:beta:gamma ~ {alpha} : {beta} : {gamma}'
          f'; beta/gamma = {float(beta / gamma):.4f} (Cr-line: -0.892)')
print(f'{len(fails)} FAILED: {fails}' if fails else 'ALL PASS')
