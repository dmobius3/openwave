"""M8.10 fresh derivation, step 5: why the frozen zeros are zero.

A ray [u] with stabilizer H (rotations of the fibre V3) has D^3(h)u = chi(h)u. Pi_n N(Phi) carries the
same character, so its fibre factor lies in the chi-isotypic part of V_J (n = 2J); if that part is zero,
the level vanishes identically and the zero is forced by H. Time reversal is antiunitary and not in H:
at a T-invariant ray, a level-16 zero that H does not force is the paper's Proposition 3.3 (M8.1.2's A5).
Diagnostic: whether ||Pi_n N||^2 = kappa(sector, J) ||[rho_6(u) (x) u]_J||^2 with kappa independent of the ray.
Writes out/zeros.json for the results generator.
"""
import json
import os
from math import cos, pi, sin

os.environ.setdefault('M810_DPS', '30')
import numpy as np
import m810_exact as E
from m810_core import cg, coupling, mvals, qmul, wigner_D

LEVELS = {'3p': [10, 14, 16, 18], '4': [8, 12, 14, 16, 18]}
T_INV = {'zonal v0', 'octahedron', 'hexagon'}
MS = mvals(3)
IX = {m: i for i, m in enumerate(MS)}


def rot(axis, angle):
    n = np.asarray(axis, float)
    return np.concatenate([[cos(angle / 2)], sin(angle / 2) * n / np.linalg.norm(n)])


POOL = [rot((0, 0, 1), 2 * pi / k) for k in (2, 3, 4, 5, 6, 20)]
POOL += [rot((cos(t), sin(t), 0), a) for t in np.arange(12) * pi / 12 for a in (pi, pi / 2)]


def closure(gens):
    e = np.array([1.0, 0.0, 0.0, 0.0])
    els, keys, frontier = [e], {tuple(e)}, [e]
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = qmul(x, g)
                k = tuple(np.round(y, 9) + 0.0)
                if k not in keys:
                    keys.add(k)
                    els.append(y)
                    new.append(y)
        frontier = new
        assert len(els) <= 2000
    return els


def stabilizer(u, pool=POOL):
    u = u / np.linalg.norm(u)
    H = closure([g for g in pool if abs(abs(np.vdot(u, wigner_D(3, g) @ u)) - 1) < 1e-10])
    chi = [np.vdot(u, wigner_D(3, h) @ u) for h in H]
    assert all(np.allclose(wigner_D(3, h) @ u, c * u) for h, c in zip(H, chi))
    return H, chi


def iso_dim(J, H, chi):
    v = sum(np.conj(c) * np.trace(wigner_D(J, h)) for h, c in zip(H, chi)) / len(H)
    assert abs(v - round(v.real)) < 1e-8, v
    return round(v.real)


def rho6(u):
    Pm, out = np.outer(u, u.conj()), np.zeros(13, complex)
    for a, N in enumerate(mvals(6)):
        for n in MS:
            if abs(N - n) <= 3:
                out[a] += cg(3, n, 3, N - n, 6, N) * (-1) ** (N - n) * Pm[IX[n], IX[n - N]]
    return out


def s(u, J):
    return np.linalg.norm(coupling(6, 3, J) @ np.kron(rho6(u), u)) ** 2


bad, record = [], {}
for sector, d in (('3p', 3), ('4', 4)):
    print(f'\n=== sector {sector} ===')
    kappa, factor_ok = {}, True
    for name, spec in E.RAYS.items():
        um = [E.mp.mpc(0)] * 7
        for m, c in spec.items():
            um[E.IX[m]] = E.mp.mpc(c)
        sc = E.mp.sqrt(E.mp.mpf(7) / d / sum(abs(x) ** 2 for x in um))
        lv, _ = E.level_norms([x * sc for x in um], E.PROJ[sector])
        lv = {n: float(lv[n]) for n in LEVELS[sector]}
        u = np.array([complex(x * sc) for x in um])
        zero = {n for n in LEVELS[sector] if abs(lv[n]) < 1e-20}
        H, chi = stabilizer(u)
        forced = {n for n in LEVELS[sector] if iso_dim(n // 2, H, chi) == 0}
        tz = {16} if name in T_INV else set()
        ok = zero == (forced | tz)
        bad += [] if ok else [f'{sector} {name}']
        record[f'{sector}|{name}'] = {'zero': sorted(zero), 'forced_by_H': sorted(forced), 'time_reversal_only': sorted(tz - forced),
                                      'su2_order_of_H': len(H)}
        print(f"  {name:12s} |H| = {len(H):3d} in SU(2)  zero levels {sorted(zero)}  forced by H {sorted(forced)}  "
              f"T only {sorted(tz - forced)}  [{'explained' if ok else 'UNEXPLAINED'}]")
        for n in LEVELS[sector]:
            sJ = s(u, n // 2)
            if sJ > 1e-14:
                kappa.setdefault(n, []).append(lv[n] / sJ)
            elif lv[n] > 1e-14:
                factor_ok = False
    spread = max(max(v) / min(v) - 1 for v in kappa.values())
    print(f'  diagnostic: ||Pi_n N||^2 / ||[rho_6(u) (x) u]_J||^2 is ray-independent: {factor_ok and spread < 1e-9} '
          f'(max spread {spread:.1e}); kappa by level: ' + ', '.join(f'{n}: {np.mean(v):.6e}' for n, v in kappa.items()))

u = np.zeros(7, complex)
u[1] = u[5] = 1
H, chi = stabilizer(u, pool=[g for g in POOL if abs(g[3] - sin(pi / 4)) < 1e-12])       # z-rotations by pi/2 only
arm = {n for n in LEVELS['4'] if iso_dim(n // 2, H, chi) == 0} | {16}
print(f'\n  arm: octahedron with H cut to C4 predicts zeros {sorted(arm)} in sector 4, which misses the computed '
      f'{record["4|octahedron"]["zero"]}: {"fails as required" if arm != set(record["4|octahedron"]["zero"]) else "DOES NOT FAIL"}')
bad += [] if arm != set(record['4|octahedron']['zero']) else ['C4 arm']
os.makedirs('out', exist_ok=True)
json.dump(record, open('out/zeros.json', 'w'), indent=1)
print('\n  ALL ZEROS EXPLAINED' if not bad else f'\n  FAILURES: {bad}')
