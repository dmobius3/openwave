"""M8.11 fresh derivation, step 2: the trigonal prism (g = 1), with the S lemma first.

The locus (Surviving Ray Sec. 5.8): H = D3 = <R_z(2 pi/3), R_x(pi)> acts by -1 on v3 + v_-3 and v0, so the
chi-isotypic space is span{v3 + v_-3, v0}; chart u = v3 + z v0 + v_-3, prism at z = +sqrt(23/10). No continuous
symmetry-orbit direction lies inside this locus (R_z(theta) preserves it only when e^{6 i theta} = 1), so the
slice is the chart itself, with horizontal unit tangents tau_x (real) and tau_y = i tau_x, oriented by
increasing x and y. They are the S-even and S-odd directions.

The S lemma in three steps (F2's order), each gated:
  1. fibre side: S = R_z(pi/3) o Theta fixes the prism ray and acts on the chart as z -> conj(z); Theta alone,
     and the octahedron, are controls;
  2. bundle level: fibre conjugation C(f)(g) = W^T conj(f(g)) maps sections of E_sigma to sections iff
     W^T conj(P) W = P (route A); an explicit real orthonormal model of sigma, real orthogonal on all 120
     elements, is route B; and C and the left rotation each commute with the cubic, level by level;
  3. S xi = xi, checked level by level, and the forcing's S-odd component, computed and checked for zero.
Then: the in-locus second variation by two routes (the machinery, and (w6/4) x the paper's chart Hessian over the
Fubini-Study metric), the forcing, the tilt, the a^5 tangential residual with and without it, lambda_4 by two
routes and with the tilt, and the prism's xi level norms; the hexagon (z = 0, pinned) on the same locus is a
control. Float64; the author's exact route is step 3. Runs from a clean directory.
"""
import json
import os
import sys
from math import cos, pi, sin, sqrt

import numpy as np
from m810_core import mvals, wigner_D
from m811_ops import DIM, DN, G, PROJ, W, W6, block, conj_term, cubic_levels, den, fibre, norm2, pair, product_level, xi_terms

E = lambda m: np.eye(7)[3 - m].astype(complex)
MS = mvals(3)
X0 = sqrt(23 / 10)
FS = 200 / 1849                              # |d u_hat/dx|^2 = |horizontal d u_hat/dy|^2 at the prism, 2/(2 + x^2)^2
H_CHART = (184000 / 874577, 1600 / 20339)    # the paper's chart formula: d2/dx2, d2/dy2 at the prism (feasibility check)
RZ = np.array([cos(pi / 6), 0.0, 0.0, sin(pi / 6)])     # R_z(pi/3)
DRZ = {J: wigner_D(J, RZ) for J in range(10)}
bad, rec, count = [], {}, [0]


def check(name, ok, detail=''):
    count[0] += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ''))
    if not ok:
        bad.append(name)


def uhat(z):
    v = E(3) + z * E(0) + E(-3)
    return v / np.linalg.norm(v)


def paper_h(x, y):                            # the paper's ||rho_6||^2 on the D3 chart
    r2 = x * x + y * y
    return (100 * r2 * r2 - 20 * x * x + 148 * y * y + 463) / (231 * (r2 + 2) ** 2)


def theta(u):                                 # the paper's time reversal: (Theta u)_m = (-1)^m conj(u_{-m})
    return np.array([(-1) ** m * np.conj(u[MS.index(-m)]) for m in MS])


def S_fibre(u):
    return DRZ[3] @ theta(u)


def overlap(a, b):
    return abs(np.vdot(a, b)) / (np.linalg.norm(a) * np.linalg.norm(b))


def C_term(t):                                # fibre conjugation on a term: g -> W3^T conj(X R^J(g) y)
    j, Xc, yc = conj_term(t)
    return j, W[3].T @ Xc, yc


def L_term(t):                                # left rotation by R_z(pi/3): (L f)(g) = f(R^-1 g) acts on y by D^J(R)
    j, X, y = t
    return j, X, DRZ[j] @ y


def tensor(terms, J):                         # the coefficient tensor that determines a level-J function uniquely
    M = np.zeros((7, 2 * J + 1, 2 * J + 1), complex)
    for _, X, y in terms:
        M += X[:, :, None] * y[None, None, :]
    return M


def quartic(u, P, d):
    B = block(u, P)
    wN, res = fibre(product_level(B, B, B, 3), P, d)
    return pair(u, wN, d).real, wN, res


def tangential(u, w):
    return w - np.vdot(u, w) / np.vdot(u, u) * u


def hess_form(u, e, P, d, conj_on=True):
    w, _ = fibre(DN(block(u, P), [block(e, P)], conj_term_on=conj_on), P, d)
    return pair(e, w, d).real - quartic(u, P, d)[0]


def chart_tangent(x):                         # d u_hat/dx at real z = x: horizontal, real
    v = E(3) + x * E(0) + E(-3)
    n2 = 2 + x * x
    return (E(0) - x / n2 * v) / sqrt(n2)


def forcing(u, P, d):
    from m811_ops import xi_terms_supported
    Phi = block(u, P)
    cub = cubic_levels(Phi)
    xi = xi_terms_supported(cub, {3: '3p', 4: '4'}[d])     # the support gate runs in the main loop (F1)
    wF, _ = fibre(DN(Phi, xi), P, d)
    return Phi, cub, xi, wF


rng = np.random.default_rng(8111)
print('=== S lemma, step 1: the fibre side (acts on the ray in V3, sector-independent) ===')
u_pr = uhat(X0)
ov = overlap(u_pr, S_fibre(u_pr))
check('S = R_z(pi/3) o Theta fixes the prism ray', abs(ov - 1) < 1e-13, f'|<u, S u>| = {ov:.15f}')
zg = 0.7 + 0.4j
ok_chart = abs(overlap(uhat(np.conj(zg)), S_fibre(uhat(zg))) - 1) < 1e-13 and overlap(uhat(zg), S_fibre(uhat(zg))) < 0.99
check('S acts on the chart as z -> conj(z) (generic point z = 0.7 + 0.4i)', ok_chart)
check('arm: Theta alone does not fix the prism ray (it sends z -> -conj(z))', overlap(u_pr, theta(u_pr)) < 0.99,
      f'{overlap(u_pr, theta(u_pr)):.6f}')
u_oct = uhat(1j * sqrt(5 / 2))
check('arm: S does not fix the octahedron point z = i sqrt(5/2)', overlap(u_oct, S_fibre(u_oct)) < 0.99,
      f'{overlap(u_oct, S_fibre(u_oct)):.6f}')

for s in ('3p', '4'):
    P, d, w6 = PROJ[s], DIM[s], W6[s]
    print(f'\n=== sector {s} (d = {d}, w6 = {w6:.12f}) ===')
    r = rec[s] = {}

    # S lemma, step 2: bundle level.
    rA = np.abs(W[3].T @ P.conj() @ W[3] - P).max()
    check('S lemma 2, route A: W^T conj(P) W = P, so fibre conjugation preserves the sector', rA < 1e-13, f'{rA:.1e}')
    V = np.diag(np.exp(1j * rng.uniform(0, 2 * pi, 7)))
    Pt = V @ P @ V.conj().T
    check('arm: a phase-twisted projector fails route A', np.abs(W[3].T @ Pt.conj() @ W[3] - Pt).max() > 1e-3)
    Cv = lambda x: W[3].T @ x.conj()          # C on fibre vectors; C^2 = 1 since W^2 = 1 at integer spin
    basis = []
    for c in P.T:                             # C-fixed candidates c + Cc and i(c - Cc), Gram-Schmidt with real coefficients
        for b in (c + Cv(c), 1j * (c - Cv(c))):
            for q in basis:
                b = b - np.vdot(q, b) * q
            if np.linalg.norm(b) > 1e-9 and len(basis) < d:
                basis.append(b / np.linalg.norm(b))
    B = np.array(basis).T
    sig = [B.conj().T @ wigner_D(3, h) @ B for h in G]
    im = max(np.abs(x.imag).max() for x in sig)
    orth = max(np.abs(x.real @ x.real.T - np.eye(d)).max() for x in sig)
    fixed = max(np.abs(Cv(b) - b).max() for b in basis)
    check('S lemma 2, route B: a C-real orthonormal basis of the sector makes sigma real orthogonal on all 120 elements',
          len(basis) == d and fixed < 1e-13 and im < 1e-12 and orth < 1e-12, f'imag {im:.1e}, orthogonality {orth:.1e}')
    Bc = B @ np.diag([1j] + [1] * (d - 1))
    im_c = max(np.abs((Bc.conj().T @ wigner_D(3, h) @ Bc).imag).max() for h in G)
    check('arm: a basis that is not C-real gives a non-real sigma', im_c > 1e-3, f'imag {im_c:.3f}')
    ur = rng.normal(size=7) + 1j * rng.normal(size=7)
    ur *= sqrt(7 / d) / np.linalg.norm(ur)
    Nu = cubic_levels(block(ur, P))
    NCu = cubic_levels(block(Cv(ur), P))
    NLu = cubic_levels(block(DRZ[3] @ ur, P))
    cC = max(np.abs(tensor([C_term(t) for t in Nu[J]], J) - tensor(NCu[J], J)).max() for J in Nu if Nu[J])
    cL = max(np.abs(tensor([L_term(t) for t in Nu[J]], J) - tensor(NLu[J], J)).max() for J in Nu if Nu[J])
    check('S lemma 2: C commutes with the cubic, level by level (random block state)', cC < 1e-12, f'{cC:.1e}')
    check('S lemma 2: the left rotation commutes with the cubic, level by level', cL < 1e-12, f'{cL:.1e}')

    # Parents at the prism: the paper's chart formula, criticality in the whole block, Q.
    errs = [abs(quartic(sqrt(7 / d) * uhat(x + 1j * y), P, d)[0] - (1 + w6 * paper_h(x, y))) for x, y in ((X0, 0), (0.4, 0.9), (1.3, -0.2))]
    check('parent: Q on the chart = 1 + w6 * (paper chart formula), three points', max(errs) < 1e-12, f'max err {max(errs):.1e}')
    u0 = sqrt(7 / d) * u_pr
    Q0, wN0, res0 = quartic(u0, P, d)
    crit = np.linalg.norm(tangential(u0, wN0))
    check('parent: the prism is critical in the whole block', crit < 1e-12, f'{crit:.1e}')
    check('parent: Q at the prism = 1 + w6 * 200/903', abs(Q0 - (1 + w6 * 200 / 903)) < 1e-12, f'Q = {Q0:.15f}')

    # The section-level S, phase fixed so that S Phi = Phi.
    Phi, cub, xi, wF = forcing(u0, P, d)
    from m811_ops import LEVELS, forbidden_norms
    fb = forbidden_norms(cub, s)
    check("support (F1): the cubic vanishes at every level outside the sector's support", max(fb.values()) < 1e-24,
          ', '.join(f'{n}: {v:.0e}' for n, v in fb.items()))
    cub = {J: ts for J, ts in cub.items() if 2 * J in LEVELS[s]}
    SPhi = DRZ[3] @ Cv(u0)
    ph = np.vdot(u0, SPhi) / np.vdot(u0, u0)
    check('S fixes the prism at the section level: S Phi = (phase) Phi', abs(abs(ph) - 1) < 1e-13 and np.abs(SPhi - ph * u0).max() < 1e-13)
    S_term = lambda t: (lambda j, X, y: (j, X, np.conj(ph) * y))(*L_term(C_term(t)))
    sxi = max(np.abs(tensor([S_term(t) for t in xi if t[0] == J], J) - tensor([t for t in xi if t[0] == J], J)).max()
              for J in {t[0] for t in xi})
    check('S lemma 3: S xi = xi, level by level (xi is S-even, by uniqueness)', sxi < 1e-12, f'{sxi:.1e}')

    # The slice: tau_x (S-even), tau_y = i tau_x (S-odd).
    tx0 = chart_tangent(X0)
    check('Fubini-Study metric at the prism: |d u_hat/dx|^2 = 200/1849', abs(np.vdot(tx0, tx0).real - FS) < 1e-15)
    tau_x = sqrt(7 / d) * tx0 / np.linalg.norm(tx0)
    tau_y = 1j * tau_x
    Stx = np.conj(ph) * (DRZ[3] @ Cv(tau_x))
    Sty = np.conj(ph) * (DRZ[3] @ Cv(tau_y))
    check('slice: tau_x is S-even and tau_y is S-odd', np.abs(Stx - tau_x).max() < 1e-13 and np.abs(Sty + tau_y).max() < 1e-13)

    # In-locus second variation, two routes.
    Hxx, Hyy = hess_form(u0, tau_x, P, d), hess_form(u0, tau_y, P, d)
    Hxy = pair(tau_x, fibre(DN(Phi, [block(tau_y, P)]), P, d)[0], d).real
    Hp = [w6 / 4 * h / FS for h in H_CHART]
    check('two routes to the in-locus second variation: machinery = (w6/4) x paper chart Hessian / FS metric',
          abs(Hxx - Hp[0]) < 1e-12 and abs(Hyy - Hp[1]) < 1e-12, f'({Hxx:.12f}, {Hyy:.12f}) vs ({Hp[0]:.12f}, {Hp[1]:.12f})')
    check('no cross term between the S-even and S-odd directions', abs(Hxy) < 1e-13, f'{Hxy:.1e}')
    check('arm: without the Fubini-Study conversion the routes disagree', abs(Hxx - w6 / 4 * H_CHART[0]) > 1e-3)
    check('arm: dropping the conjugate term of DN breaks the agreement', abs(hess_form(u0, tau_x, P, d, conj_on=False) - Hp[0]) > 1e-3)

    # The forcing, its decomposition, and the S-odd component.
    tanF = tangential(u0, wF)
    zx = pair(tau_x, tanF, d)                 # real part along tau_x, imaginary part along tau_y
    outside = np.linalg.norm(tanF - zx * tau_x)
    check('forcing lies in the fixed kernel span{v3 + v_-3, v0}', outside < 1e-12, f'{outside:.1e}')
    check('S lemma 3: the forcing has no S-odd (tau_y) component', abs(zx.imag) < 1e-14, f'{zx.imag:.1e}')
    zg2 = X0 + 0.3j                           # a chart point S does not fix: the same component need not vanish
    ug = sqrt(7 / d) * uhat(zg2)
    wFg = forcing(ug, P, d)[3]
    tg = sqrt(7 / d) * ((E(0) - (zg2 / (2 + abs(zg2) ** 2)) * (E(3) + zg2 * E(0) + E(-3))))
    tg -= np.vdot(ug, tg) / np.vdot(ug, ug) * ug
    tg *= sqrt(7 / d) / np.linalg.norm(tg)
    comp = pair(tg, tangential(ug, wFg), d)
    check('arm: at the non-S-fixed chart point z = sqrt(23/10) + 0.3i the corresponding component is nonzero', abs(comp.imag) > 1e-6,
          f'{comp.imag:.3e}')
    uh = sqrt(7 / d) * uhat(0)                # the hexagon: pinned, on the same locus
    wFh = forcing(uh, P, d)[3]
    check('control: at the hexagon (z = 0, pinned) on the same locus the forcing has no tangential part',
          np.linalg.norm(tangential(uh, wFh)) < 1e-12, f'{np.linalg.norm(tangential(uh, wFh)):.1e}')

    # The tilt: solve the 2x2 system on the slice.
    Hm = np.array([[Hxx, Hxy], [Hxy, Hyy]])
    vx, vy = np.linalg.solve(Hm, -np.array([zx.real, zx.imag]))
    v = vx * tau_x + vy * tau_y
    wT, _ = fibre(DN(Phi, xi + [block(v, P)]), P, d)
    resid_with = np.linalg.norm(tangential(u0, wT) - Q0 * v)
    check('a^5 tangential residual vanishes with the tilt', resid_with < 1e-12, f'{resid_with:.1e}')
    nonzero = abs(zx.real) > 1e-10
    print(f"  the prism's forcing along tau_x is {'NONZERO' if nonzero else 'ZERO to working precision'}: {zx.real:+.15e}")
    if nonzero:
        check('arm: without the tilt the tangential residual is nonzero', np.linalg.norm(tanF) > 1e-6, f'{np.linalg.norm(tanF):.3e}')

    # lambda_4: two routes, and with the tilt.
    lv = {2 * J: (norm2(ts) if ts else 0.0) for J, ts in cub.items()}
    r1 = -3 * sum(val / den(n) for n, val in lv.items() if n != 6)
    r2, r2v = pair(u0, wF, d), pair(u0, wT, d)
    check('lambda_4, two routes agree', abs(r2 - r1) < 1e-12 and abs(r2.imag) < 1e-13, f'{r1:.15e}')
    check('lambda_4 lemma: including the tilt does not change it', abs(r2v - r2) < 1e-13, f'difference {abs(r2v - r2):.1e}')
    uc = sqrt(7 / d) * uhat(1.0)
    tc = chart_tangent(1.0)
    tc = sqrt(7 / d) * tc / np.linalg.norm(tc)
    Lc = pair(uc, fibre(DN(block(uc, P), [block(tc, P)]), P, d)[0], d).real
    check('arm: at the non-critical chart point z = 1 the lemma pairing is nonzero', abs(Lc) > 1e-3, f'{Lc:.6e}')

    per_level = {}
    for J in sorted(cub):
        xiJ = [t for t in xi if t[0] == J]
        if J != 3 and xiJ:
            per_level[2 * J] = pair(tau_x, tangential(u0, fibre(DN(Phi, xiJ), P, d)[0]), d).real
    support = {n: val for n, val in lv.items() if n != 6 and val > 1e-14}
    print('  forcing along tau_x by level of xi: ' + ', '.join(f'{n}: {val:+.6e}' for n, val in per_level.items()))
    print('  level norms ||P_n N||^2: ' + ', '.join(f'{n}: {val:.12e}' for n, val in support.items()))
    print(f'  Q = {Q0:.15f}   second variation (tau_x, tau_y) = ({Hxx:.15e}, {Hyy:.15e})')
    print(f'  forcing (tau_x, tau_y) = ({zx.real:+.15e}, {zx.imag:+.1e})   tilt (v_x, v_y) = ({vx:+.15e}, {vy:+.1e})   ||v||^2 = {vx * vx + vy * vy:.15e}')
    print(f'  lambda_4/g^2 = {r1:.15e}')
    r.update(Q=Q0, hess=[Hxx, Hyy, Hxy], forcing=[zx.real, zx.imag], tilt=[vx, vy], tilt_norm2=vx * vx + vy * vy, lambda4=r1,
             forcing_by_level=per_level, level_norms_N=support, xi_level_norms={n: val / den(n) ** 2 for n, val in support.items()})

os.makedirs('out', exist_ok=True)
json.dump(rec, open('out/prism.json', 'w'), indent=1)
print(f'\n  {count[0]} checks; ' + ('ALL PASS, ALL ARMS FIRE' if not bad else f'FAILURES: {bad}'))
sys.exit(1 if bad else 0)
