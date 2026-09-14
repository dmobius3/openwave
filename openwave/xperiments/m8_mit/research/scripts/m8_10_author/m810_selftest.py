"""Positive controls and planted failures for m810_core. Every PASS line has a failing twin."""
import numpy as np
from m810_core import (binary_icosahedral, qmul, su2, wigner_D, conj_matrix, coupling,
                       isotypic_projector, DIM)

rng = np.random.default_rng(810)
bad = []


def check(name, ok, detail=''):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        bad.append(name)


def rand_q():
    q = rng.normal(size=4)
    return q / np.linalg.norm(q)


G = binary_icosahedral()
keys = {tuple(np.round(q, 9)) for q in G}
closed = all(tuple(np.round(qmul(p, q), 9)) in keys for p in G for q in G)
check('2I: 120 distinct unit quaternions, closed under product',
      len(G) == 120 and len(keys) == 120 and np.allclose(np.linalg.norm(G, axis=1), 1) and closed)
G_odd = G.copy(); G_odd[-1] = G_odd[-1][[0, 2, 1, 3]]            # planted: one odd permutation
keys_odd = {tuple(np.round(q, 9)) for q in G_odd}
check('  planted: one odd-permutation element breaks closure',
      not all(tuple(np.round(qmul(p, q), 9)) in keys_odd for p in G_odd for q in G_odd))

p, q = rand_q(), rand_q()
check('U(pq) = U(p)U(q)', np.allclose(su2(qmul(p, q)), su2(p) @ su2(q)))
for j in (3, 9):
    Dp, Dq, Dpq = wigner_D(j, p), wigner_D(j, q), wigner_D(j, qmul(p, q))
    check(f'D^{j}: homomorphism and unitary', np.allclose(Dpq, Dp @ Dq) and np.allclose(Dp @ Dp.conj().T, np.eye(2 * j + 1)))
    W = conj_matrix(j)
    check(f'D^{j}: conj(D) = W D W^T', np.allclose(Dp.conj(), W @ Dp @ W.T))
check('  planted: D^3 of the reversed product fails the homomorphism',
      not np.allclose(wigner_D(3, qmul(q, p)), wigner_D(3, p) @ wigner_D(3, q)))

Cs = {J: coupling(3, 3, J) for J in range(7)}
check('CG (3x3): rows orthonormal, complete', all(np.allclose(C @ C.T, np.eye(2 * J + 1)) for J, C in Cs.items())
      and np.allclose(sum(C.T @ C for C in Cs.values()), np.eye(49)))
D3 = wigner_D(3, p)
check('CG intertwines: C_J (D3 x D3) = D^J C_J', all(np.allclose(C @ np.kron(D3, D3), wigner_D(J, p) @ C) for J, C in Cs.items()))
Cbad = Cs[4].copy(); Cbad[2, 10] *= -1
check('  planted: one flipped CG sign breaks the intertwining', not np.allclose(Cbad @ np.kron(D3, D3), wigner_D(4, p) @ Cbad))

P3p, P4 = isotypic_projector(3, '3p', G), isotypic_projector(3, '4', G)
for name, P in (('3p', P3p), ('4', P4)):
    ok = np.allclose(P @ P, P) and np.allclose(P, P.conj().T) and round(np.trace(P).real) == DIM[name]
    ok = ok and all(np.allclose(P @ wigner_D(3, h), wigner_D(3, h) @ P) for h in G)
    check(f'P_{name} on V3: Hermitian idempotent of rank {DIM[name]}, commutes with 2I', ok)
check('P_3p + P_4 = I on V3', np.allclose(P3p + P4, np.eye(7)))
check('  planted: the other 3-dim irrep (3, not 3p) has rank 0 on V3', abs(np.trace(isotypic_projector(3, '3', G))) < 1e-9)

expect = {'3p': {3: 1, 5: 1, 7: 1, 8: 1, 9: 1}, '4': {3: 1, 4: 1, 6: 1, 7: 1, 8: 1, 9: 2}}
for s in ('3p', '4'):
    got = {j: round(np.trace(isotypic_projector(j, s, G)).real / DIM[s]) for j in range(10)}
    got = {j: m for j, m in got.items() if m}
    check(f'sector {s}: multiplicities by level', got == expect[s], f'levels {[2 * j for j in got]} mult {list(got.values())}')

print('\n  ALL PASS' if not bad else f'\n  FAILURES: {bad}')
