"""M8.10 fresh derivation, step 2b: parent checks in the paper's own conventions.

Paper (Zenodo tex): M_K(P)_N = sum_{n+n'=N} <3 n; 3 n' | K N> (-1)^{n'} P_{n,-n'} (line 227);
rho_K(u) = [u (x) Theta u]_K = M_K(u u^dagger), (Theta u)_{n'} = (-1)^{n'} conj(u_{-n'}) (240-257);
R_K = M_K(P); M_K(u) = [rho_K(u) (x) u]_3 (652, 776); w_K = (49/d^2) ||R_K(P)||^2/(2K+1) (834);
C3: the block projection of the cubic, in fibre form, is (d/7)||u||^2 u - (7-d)/sqrt(91) M_6(u) (981).
"""
from math import comb, sqrt

import numpy as np
from m810_core import binary_icosahedral, conj_matrix, coupling, isotypic_projector, mvals, cg, DIM

G = binary_icosahedral()
W = conj_matrix(3)
ms = mvals(3)
idx = {m: i for i, m in enumerate(ms)}
rng = np.random.default_rng(1810)
bad = []


def check(name, ok, detail=''):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        bad.append(name)


def multipole(Pm, K):
    out = np.zeros(2 * K + 1, complex)
    for a, N in enumerate(mvals(K)):
        for n in ms:
            n2 = N - n
            if abs(n2) <= 3:
                out[a] += cg(3, n, 3, n2, K, N) * (-1) ** n2 * Pm[idx[n], idx[-n2]]
    return out


def theta(u):
    return np.array([(-1) ** m * np.conj(u[idx[-m]]) for m in ms])


def bracket(x, jx, y, jy, J):
    return coupling(jx, jy, J) @ np.kron(x, y)


def rho(u, K):
    return multipole(np.outer(u, u.conj()), K)


def M(u, K):
    return bracket(rho(u, K), K, u, 3, 3)


def rand_u():
    u = rng.normal(size=7) + 1j * rng.normal(size=7)
    return u / np.linalg.norm(u)


us = [rand_u() for _ in range(4)]
check('rho_K(u) = [u (x) Theta u]_K, K = 0..6', all(np.allclose(rho(u, K), bracket(u, 3, theta(u), 3, K)) for u in us for K in range(7)))
check('M_0(u) = -||u||^2 u/sqrt(7)', all(np.allclose(M(u, 0), -u / sqrt(7)) for u in us))
check('<u, M_6(u)> = -(sqrt(91)/13) ||rho_6||^2',
      all(np.isclose(np.vdot(u, M(u, 6)), -sqrt(91) / 13 * np.linalg.norm(rho(u, 6)) ** 2) for u in us))
e = lambda m: np.eye(7)[idx[m]].astype(complex)
check('weight states: ||rho_6(v_m)||^2 = C(6,3+m)^2/C(12,6)',
      all(np.isclose(np.linalg.norm(rho(e(m), 6)) ** 2, comb(6, 3 + m) ** 2 / comb(12, 6)) for m in ms))
rays = {'coherent v3': e(3), 'zonal v0': e(0), 'octahedron': (e(2) + e(-2)) / sqrt(2), 'hexagon': (e(3) + e(-3)) / sqrt(2)}
vals = [924 * np.linalg.norm(rho(u, 6)) ** 2 for u in rays.values()]
check('paper table: 924 ||rho_6||^2 = 1, 400, 288, 463', np.allclose(vals, [1, 400, 288, 463]), str(np.round(vals, 9)))

C2 = {J: coupling(3, 3, J) for J in range(7)}
C3to3 = [coupling(J12, 3, 3) @ np.kron(C2[J12], np.eye(7)) for J12 in range(7)]


def level6_fibre(u, P):
    """Fibre vector w with (level-6 part of N(Phi'_u))(g) = P D(g)^dagger w, and the residual of T_b = w_b P."""
    Mm, ut = W.T @ P, W.T @ u.conj()
    L = np.einsum('ab,kc->kabc', Mm, P).reshape(7, 343)
    U3 = np.kron(np.kron(ut, u), u)
    T = sum(np.einsum('ka,b->kab', L @ C.T, C @ U3) for C in C3to3)
    wv = np.array([np.trace(T[:, :, b] @ P) for b in range(7)]) / np.trace(P).real
    return wv, max(np.abs(T[:, :, b] - wv[b] * P).max() for b in range(7))


for s, target in (('3p', 28 / 39), ('4', 21 / 52)):
    P, d = isotypic_projector(3, s, G), DIM[s]
    w = [49 / d ** 2 * np.linalg.norm(multipole(P, K)) ** 2 / (2 * K + 1) for K in range(7)]
    check(f'sector {s}: R_K(P) = 0 for K = 1..5, w_0 = 7, w_6 = {target:.9f}',
          np.allclose(w[1:6], 0) and np.isclose(w[0], 7) and np.isclose(w[6], target), str(np.round(w, 9)))
    errs, errs_bad, res = [], [], []
    for u in us:
        wv, r = level6_fibre(u, P)
        good = d / 7 * u - (7 - d) / sqrt(91) * M(u, 6)
        flipped = d / 7 * u + (7 - d) / sqrt(91) * M(u, 6)
        errs.append(np.abs(wv - good).max()); errs_bad.append(np.abs(wv - flipped).max()); res.append(r)
    check(f'C3, sector {s}: level-6 part of N = (d/7)||u||^2 u - (7-d)/sqrt(91) M_6(u), random u',
          max(errs) < 1e-12, f'(max err {max(errs):.1e}, equivariance residual {max(res):.1e})')
    check(f'  planted, sector {s}: the flipped M_6 sign fails', min(errs_bad) > 1e-3, f'(min err {min(errs_bad):.1e})')

print('\n  ALL PASS' if not bad else f'\n  FAILURES: {bad}')
