"""M8.10 fresh derivation, step 3: the level norms at high precision (M810_DPS, default 60), identified as
rationals. Rational identification, not a symbolic derivation: the fractions are candidate exact values.

Same algebra as m810_levels, written as sparse Clebsch-Gordan loops in mpmath. Imports only
m810_core (for the exact Racah formula and the permutation parity) and mpmath.
"""
import itertools
import os
from fractions import Fraction
from functools import lru_cache

import mpmath as mp
from m810_core import cg_exact, _parity

mp.mp.dps = int(os.environ.get('M810_DPS', '60'))
PHI = (1 + mp.sqrt(5)) / 2
MS = [3, 2, 1, 0, -1, -2, -3]
IX = {m: i for i, m in enumerate(MS)}


@lru_cache(maxsize=None)
def cg(j1, m1, j2, m2, J, M):
    s, pre = cg_exact(j1, m1, j2, m2, J, M)
    return mp.mpf(0) if s == 0 else mp.mpf(s.numerator) / s.denominator * mp.sqrt(mp.mpf(pre.numerator) / pre.denominator)


def spin3():
    jp = mp.zeros(7, 7)
    for c in range(1, 7):
        jp[c - 1, c] = mp.sqrt(12 - MS[c] * (MS[c] + 1))
    jm = jp.T
    return (jp + jm) / 2, (jp - jm) * mp.mpc(0, -0.5), mp.diag([mp.mpf(m) for m in MS])


JX, JY, JZ = spin3()


def D3(q):
    a, v = q[0], q[1:]
    s = mp.sqrt(sum(x * x for x in v))
    if s < mp.mpf(10) ** -50:
        return mp.eye(7)
    th = 2 * mp.atan2(s, a)
    return mp.expm(mp.mpc(0, -1) * th * (v[0] / s * JX + v[1] / s * JY + v[2] / s * JZ))


def group():
    half = mp.mpf(1) / 2
    els = []
    for i in range(4):
        for sg in (1, -1):
            q = [mp.mpf(0)] * 4
            q[i] = mp.mpf(sg)
            els.append(q)
    els += [list(s) for s in itertools.product((half, -half), repeat=4)]
    base = (mp.mpf(0), half, half / PHI, half * PHI)
    for p in itertools.permutations(range(4)):
        if _parity(p) == 0:
            for sb, sc, sd in itertools.product((1, -1), repeat=3):
                v = (base[0], sb * base[1], sc * base[2], sd * base[3])
                els.append([v[p[k]] for k in range(4)])
    return els


CLS = [(mp.mpf(1), 3, 4), (mp.mpf(-1), -1, 0), (mp.mpf(-0.5), 0, 1), ((PHI - 1) / 2, 1 - PHI, -1), (-PHI / 2, PHI, -1)]


def char(sector, q):
    c = 2 * q[0] ** 2 - 1
    hit = [k for k in CLS if abs(k[0] - c) < mp.mpf(10) ** -(mp.mp.dps // 2)]   # classes are >= 0.19 apart
    assert len(hit) == 1
    return hit[0][1] if sector == '3p' else hit[0][2]


G = group()
assert len(G) == 120
DG = [D3(q) for q in G]
PROJ = {s: sum((char(s, q) * Dq for q, Dq in zip(G, DG)), mp.zeros(7, 7)) * (d / mp.mpf(120)) for s, d in (('3p', 3), ('4', 4))}


def level_norms(u, P):
    ut = [(-1) ** (3 + m) * mp.conj(u[IX[-m]]) for m in MS]
    Mm = [[(-1) ** (3 + a) * P[IX[-a], IX[b]] for b in MS] for a in MS]
    out, imag = {}, mp.mpf(0)
    for J in range(10):
        Xs, ys = [], []
        for J12 in range(7):
            if not abs(J12 - 3) <= J <= J12 + 3:
                continue
            y = [mp.mpc(0)] * (2 * J + 1)
            X = [[mp.mpc(0)] * (2 * J + 1) for _ in range(7)]
            for iM, Mv in enumerate(range(J, -J - 1, -1)):
                for m3 in MS:
                    M12 = Mv - m3
                    if abs(M12) > J12:
                        continue
                    c2 = cg(J12, M12, 3, m3, J, Mv)
                    if c2 == 0:
                        continue
                    for m1 in MS:
                        m2 = M12 - m1
                        if abs(m2) > 3:
                            continue
                        B = cg(3, m1, 3, m2, J12, M12) * c2
                        if B == 0:
                            continue
                        y[iM] += B * ut[IX[m1]] * u[IX[m2]] * u[IX[m3]]
                        mb = B * Mm[IX[m1]][IX[m2]]
                        for k in range(7):
                            X[k][iM] += mb * P[k, IX[m3]]
            Xs.append(X)
            ys.append(y)
        s = mp.mpc(0)
        for a in range(len(Xs)):
            for b in range(len(Xs)):
                tr = mp.fsum(mp.conj(Xs[b][k][i]) * Xs[a][k][i] for k in range(7) for i in range(2 * J + 1))
                s += tr * mp.fsum(mp.conj(ys[b][i]) * ys[a][i] for i in range(2 * J + 1))
        out[2 * J] = mp.re(s) / (2 * J + 1)
        imag = max(imag, abs(mp.im(s)))
    return out, imag


def as_fraction(x, maxden=10 ** 16):
    f = Fraction(mp.nstr(x, mp.mp.dps - 2)).limit_denominator(maxden)
    return f, abs(x - mp.mpf(f.numerator) / f.denominator)


RAYS = {'coherent v3': {3: 1}, 'zonal v0': {0: 1}, 'octahedron': {2: 1, -2: 1}, 'hexagon': {3: 1, -3: 1}}
TABLE = {'coherent v3': 1, 'zonal v0': 400, 'octahedron': 288, 'hexagon': 463}
DEN = {'3p': 1287, '4': 2288}


def den(n):
    return n * (n + 2) - 48


if __name__ == '__main__':
    worst = mp.mpf(0)
    for sector, d in (('3p', 3), ('4', 4)):
        P = PROJ[sector]
        print(f'\n=== sector {sector} ===')
        for name, spec in RAYS.items():
            u = [mp.mpc(0)] * 7
            for m, c in spec.items():
                u[IX[m]] = mp.mpc(c)
            scale = mp.sqrt(mp.mpf(7) / d / sum(abs(x) ** 2 for x in u))
            u = [x * scale for x in u]
            lv, imag = level_norms(u, P)
            exact = {}
            for n, x in lv.items():
                f, err = as_fraction(x)
                worst = max(worst, err)
                if f != 0:
                    exact[n] = f
            Q = 1 + Fraction(TABLE[name], DEN[sector])
            lam4 = -3 * sum(f / den(n) for n, f in exact.items() if n != 6)
            xi2 = sum(f / den(n) ** 2 for n, f in exact.items() if n != 6)
            print(f'  {name}: |P6 N|^2 = Q_d^2 exactly: {exact.get(6) == Q * Q}   imag {mp.nstr(imag, 3)}')
            print('    ||P_n N||^2: ' + ', '.join(f'n={n}: {f}' for n, f in exact.items() if n != 6))
            print(f'    lambda4/g^2 = {lam4}  = {float(lam4):.12e}')
            print(f'    ||xi||^2/g^2 = {xi2}  (||xi|| = {float(xi2) ** 0.5:.12e})')
            print('    ||P_n xi||^2/g^2: ' + ', '.join(f'n={n}: {f / den(n) ** 2}' for n, f in exact.items() if n != 6))
    print(f'\n  worst rational-identification residual: {mp.nstr(worst, 3)}  (needs < 1e-{mp.mp.dps - 15}, dps = {mp.mp.dps})')
