"""M8.10 fresh derivation, step 2: level decomposition of N(Phi) at the four rays, both sectors.

Phi' = P D(g)^dagger u is the V3-valued isometric image of the section (P = eta eta^dagger), with
||u||^2 = 7/d so that int |Phi'|^2 = 1. With R(g) = D(g)^dagger, u~ = W^T conj(u) and
L(x (x) y (x) z) = [(Wx)^T P y] P z, the cubic is N(Phi')(g) = L R(g)^{x3} (u~ (x) u (x) u).
Its level-2J part is L R^{x3} E_J U, E_J the J-isotypic projector of V3^{x3}; norms follow
from Schur orthogonality. Finite algebra (Clebsch-Gordan and 2I averages), no quadrature; float64 here,
60 and 100 digits in m810_exact.
"""
import numpy as np
from m810_core import binary_icosahedral, conj_matrix, coupling, isotypic_projector, DIM

G = binary_icosahedral()
W = conj_matrix(3)
C2 = {J: coupling(3, 3, J) for J in range(7)}
C3 = {}
for J12 in range(7):
    K = np.kron(C2[J12], np.eye(7))
    for J in range(abs(J12 - 3), J12 + 4):
        C3[(J12, J)] = coupling(J12, 3, J) @ K
RAYS = {'coherent v3': {0: 1.0}, 'zonal v0': {3: 1.0},
        'octahedron': {1: 1.0, 5: 1.0}, 'hexagon': {0: 1.0, 6: 1.0}}   # index i = 3 - m
TABLE = {'coherent v3': 1, 'zonal v0': 400, 'octahedron': 288, 'hexagon': 463}   # paper: 924 ||rho_6||^2
DEN = {'3p': 1287, '4': 2288}            # w_6(d)/924: Q_d = 1 + TABLE/DEN
LEVELS = {'3p': {6, 10, 14, 16, 18}, '4': {6, 8, 12, 14, 16, 18}}
T_INVARIANT = {'zonal v0', 'octahedron', 'hexagon'}
SEPT = {('4', 'hexagon'): (7.142118e-3, -5.852e-3), ('3p', 'hexagon'): (3.602677e-3, -4.186e-3)}


def ray_vector(spec, d):
    u = np.zeros(7, complex)
    for i, c in spec.items():
        u[i] = c
    return u * np.sqrt(7 / d) / np.linalg.norm(u)


def level_norms(u, P, use_W=True):
    Wm = W if use_W else np.eye(7)
    M, ut = Wm.T @ P, Wm.T @ u.conj()
    m2, U2 = M.reshape(-1), np.kron(ut, u)
    norm2 = (m2 @ C2[0].T @ C2[0] @ U2).real
    quartic = sum(np.linalg.norm(m2 @ C2[J].T) ** 2 * np.linalg.norm(C2[J] @ U2) ** 2 / (2 * J + 1) for J in range(7))
    L = np.einsum('ab,kc->kabc', M, P).reshape(7, 343)
    U3 = np.kron(U2, u)
    lv = {}
    for J in range(10):
        ks = [k for k in C3 if k[1] == J]
        X, y = [L @ C3[k].T for k in ks], [C3[k] @ U3 for k in ks]
        s = sum(np.trace(X[b].conj().T @ X[a]) * np.vdot(y[b], y[a]) for a in range(len(ks)) for b in range(len(ks)))
        lv[2 * J] = (s / (2 * J + 1)).real
    return norm2, quartic, lv


def den(n):
    return n * (n + 2) - 48


bad = []
def check(name, ok):
    if not ok:
        bad.append(name)
    return 'ok' if ok else 'FAIL'


for s in ('3p', '4'):
    P, d = isotypic_projector(3, s, G), DIM[s]
    print(f'\n=== sector {s} (d = {d}) ===')
    for name, spec in RAYS.items():
        u = ray_vector(spec, d)
        n2, Q, lv = level_norms(u, P)
        Qp = 1 + TABLE[name] / DEN[s]
        support = {n for n, v in lv.items() if v > 1e-12}
        lam4 = -3 * sum(v / den(n) for n, v in lv.items() if n != 6)
        xi = np.sqrt(sum(v / den(n) ** 2 for n, v in lv.items() if n != 6))
        print(f'  {name:12s} int|Phi|^2={n2:.12f} [{check(name + " norm", abs(n2 - 1) < 1e-12)}]  '
              f'Q_d={Q:.12f} vs paper {Qp:.12f} [{check(name + " Q_d", abs(Q - Qp) < 1e-12)}]  '
              f'|P6 N|^2-Q^2={lv[6] - Q * Q:+.1e} [{check(name + " critical", abs(lv[6] - Q * Q) < 1e-11)}]')
        print(f'    support {sorted(support)} [{check(name + " support", support <= LEVELS[s])}]  '
              f'level16={lv[16]:.3e} [{check(name + " level16", (lv[16] < 1e-12) == (name in T_INVARIANT))}]  '
              f'lambda4={lam4:.10f}  ||xi||={xi:.9e}')
        print('    ' + '  '.join(f'n={n}:{v:.10f}' for n, v in lv.items() if n != 6 and v > 1e-14))
        if (s, name) in SEPT:
            x0, l0 = SEPT[(s, name)]
            print(f'    September run: ||xi||={x0:.6e} (rel {xi / x0 - 1:+.1e})  lambda4 at a=0.1: {l0:.4e} (rel {lam4 / l0 - 1:+.1e})')
    u = ray_vector(RAYS['hexagon'], d)
    Qbad = level_norms(u, P, use_W=False)[1]
    print(f'  planted: dropping W from conj(D) gives Q_d={Qbad:.6f}, '
          f'{"mismatch as required" if abs(Qbad - (1 + 463 / DEN[s])) > 1e-6 else "NO MISMATCH (control failed)"}')
    if abs(Qbad - (1 + 463 / DEN[s])) <= 1e-6:
        bad.append('planted W control')

print('\n  ALL CHECKS PASS' if not bad else f'\n  FAILURES: {bad}')
