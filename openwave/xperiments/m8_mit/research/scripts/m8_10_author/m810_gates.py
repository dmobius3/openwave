"""M8.10 fresh derivation, step 4: the gates that can fail, each with a mutation arm.

(a) the a^3 range equation, with A applied through the SU(2) Casimir rather than the eigenvalue table;
(b) lambda_4 by a second route: the block part of DN_Phi[xi] is built directly as products of
    coefficient functions and paired with Phi, and its block component orthogonal to Phi must vanish.
A coefficient function is g -> X R^j(g) y with R(g) = D(g)^dagger; its conjugate is
g -> conj(X) W_j R^j(g) W_j^T conj(y). The cubic itself is rebuilt here by the same product
routine, so the level norms are recomputed along a second code path. g = 1 throughout.
"""
from math import sqrt

import numpy as np
from m810_core import binary_icosahedral, conj_matrix, coupling, isotypic_projector, spin_matrices, DIM

G = binary_icosahedral()
W = {j: conj_matrix(j) for j in range(10)}
_cp = {}


def cp(j1, j2, J):
    if (j1, j2, J) not in _cp:
        _cp[(j1, j2, J)] = coupling(j1, j2, J)
    return _cp[(j1, j2, J)]


def channels(a, b, c, Jout):
    return [cp(J12, c, Jout) @ np.kron(cp(a, b, J12), np.eye(2 * c + 1))
            for J12 in range(abs(a - b), a + b + 1) if abs(J12 - c) <= Jout <= J12 + c]


def conj_term(t):
    j, X, y = t
    return j, X.conj() @ W[j], W[j].T @ y.conj()


def product_level(tA, tB, tC, Jout):
    """Level-2Jout part of g -> [sum_k conj(A_k) B_k] C, as terms (Jout, X, y)."""
    a, XA, yA = conj_term(tA)
    b, XB, yB = tB
    c, XC, yC = tC
    Lam = np.einsum('ij,kl->kijl', XA.T @ XB, XC).reshape(7, -1)
    Y = np.kron(np.kron(yA, yB), yC)
    return [(Jout, Lam @ Ch.T, Ch @ Y) for Ch in channels(a, b, c, Jout)]


def norm2(terms):
    J = terms[0][0]
    return sum(np.trace(Xs.conj().T @ Xt) * np.vdot(ys, yt) for _, Xt, yt in terms for _, Xs, ys in terms).real / (2 * J + 1)


def fibre(terms, P, d):
    T = sum(np.einsum('ka,b->kab', X, y) for _, X, y in terms)
    w = np.array([np.trace(T[:, :, b] @ P) for b in range(7)]) / d
    return w, max(np.abs(T[:, :, b] - w[b] * P).max() for b in range(7))


def den(n):
    return n * (n + 2) - 48


def run(u, P, d, drop_conj_term=False, den_fn=den):
    Phi = (3, P.astype(complex), u)
    cubic = {J: product_level(Phi, Phi, Phi, J) for J in range(10)}
    lv = {2 * J: (norm2(t) if t else 0.0) for J, t in cubic.items()}
    xi = [(J, X, -y / den_fn(2 * J)) for J, ts in cubic.items() if J != 3 for (_, X, y) in ts]
    cas = 0.0
    for J, ts in cubic.items():
        if J == 3:
            continue
        C = sum(S @ S for S in spin_matrices(J))
        for (_, _, y), (_, _, yt) in zip(ts, [t for t in xi if t[0] == J]):
            cas = max(cas, np.abs((4 * C - 48 * np.eye(2 * J + 1)) @ yt + y).max())
    DN = []
    for t in xi:
        DN += product_level(Phi, t, Phi, 3) + product_level(Phi, Phi, t, 3)
        if not drop_conj_term:
            DN += product_level(t, Phi, Phi, 3)
    w, res = fibre(DN, P, d)
    route1 = -3 * sum(v / den(n) for n, v in lv.items() if n != 6)
    route2 = d / 7 * np.vdot(u, w)
    tang = np.linalg.norm(w - np.vdot(u, w) / np.vdot(u, u) * u)
    return lv, route1, route2, tang, res, cas


e = lambda m: np.eye(7)[3 - m].astype(complex)
RAYS = {'coherent v3': e(3), 'zonal v0': e(0), 'octahedron': (e(2) + e(-2)) / sqrt(2), 'hexagon': (e(3) + e(-3)) / sqrt(2)}
t = np.arcsin(sqrt(12 / 25))
PYRAMID = np.cos(t) * e(2) + np.sin(t) * e(-3)
bad = []
for s in ('3p', '4'):
    P, d = isotypic_projector(3, s, G), DIM[s]
    print(f'\n=== sector {s} ===')
    for name, v in RAYS.items():
        u = v * sqrt(7 / d)
        lv, r1, r2, tang, res, cas = run(u, P, d)
        ok = abs(r2 - r1) < 1e-12 and abs(r2.imag) < 1e-13 and tang < 1e-12 and res < 1e-12 and cas < 1e-12
        bad += [] if ok else [f'{s} {name}']
        print(f"  {name:12s} route1 {r1:.12e}  route2 {r2.real:.12e} (imag {abs(r2.imag):.0e})  tangential {tang:.0e}  "
              f"equivariance {res:.0e}  Casimir range {cas:.0e}  [{'PASS' if ok else 'FAIL'}]")
    u = RAYS['hexagon'] * sqrt(7 / d)
    r1, r2 = run(u, P, d)[1], run(u, P, d, drop_conj_term=True)[2]
    print(f'  arm: drop (xi^dagger Phi)Phi at the hexagon -> route2/route1 = {(r2 / r1).real:.6f} (must be 2/3)')
    cas_bad = run(u, P, d, den_fn=lambda n: den(n + 2) if n == 18 else den(n))[5]
    print(f'  arm: level-18 denominator shifted -> Casimir range residual {cas_bad:.2e} (must be large)')
    tang_pyr = run(PYRAMID * sqrt(7 / d), P, d)[3]
    print(f'  arm: pentagonal pyramid (Lemma 5.3(b) ray) -> tangential {tang_pyr:.3e} (the gate can fail off the four rays)')
    bad += [] if abs((r2 / r1).real - 2 / 3) < 1e-9 and cas_bad > 1e-3 else [f'{s} arms']

print('\n  ALL GATES PASS, ALL ARMS FAIL AS REQUIRED' if not bad else f'\n  FAILURES: {bad}')
