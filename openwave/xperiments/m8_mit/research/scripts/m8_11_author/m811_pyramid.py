"""M8.11 fresh derivation, step 1: the normalization bridge, then the pentagonal pyramid (g = 1).

F1's order: fix the operator normalization before any new number. For a block section Phi with int |Phi|^2 = 1
and a unit block section e orthogonal to Phi, the unit-speed great circle cos s Phi + sin s e carries
Q(s) = int |Phi(s)|^4. In the real Hilbert structure DQ[e] = 4 Re<N(Phi), e>, and at a critical point
d^2Q/ds^2 = 4 (Re<e, DN_Phi[e]> - Q). So the tilt operator Pi_T DN_Phi - Q is c Hess_T Q with c = 1/4
predicted; F2 read 1/2, the Wirtinger convention. c is checked by three routes through the shared algebra
(m811_ops): first order, second order, and at a pinned ray. They share code and are not independent
implementations; independence is for the author's exact route (step 3) and the blind agents.

Then at the pyramid (Surviving Ray Sec. 5.8: cos t v2 + sin t v_-3, sin^2 t = 12/25): Q against the paper's
restriction; criticality; the forcing Pi_T DN_Phi[xi] against M8.10's audited residual (a parent: #547 graded
it as D2); the tilt; the a^5 tangential residual with and without it; lambda_4 with and without it;
consistency with #547's blind R5 diagnostics (items 4 and 7, asked at go under the as-run worklist and
recorded apart from the verdicts, so reproduced here as ungraded values, not as parents); and the lambda_4
lemma's gate armed at the non-critical control point sin^2 t = 1/4. Float64; the author's exact route is
step 3. Every check exits nonzero on failure, every mutation arm must fire, and the script runs from a clean
directory.
"""
import json
import os
import sys
from math import asin, cos, sin, sqrt

import numpy as np
from m811_ops import DIM, DN, PROJ, W6, block, cubic_levels, den, fibre, norm2, pair, product_level, xi_terms

E = lambda m: np.eye(7)[3 - m].astype(complex)
T0 = asin(sqrt(12 / 25))                  # the pyramid
TC = asin(sqrt(1 / 4))                    # the non-critical control point on the same curve
AUDIT = {'3p': -7188839 * sqrt(91) / 81061695000,          # M8.10 method note 5.4: coefficient of the a^5
         '4': 19565553 * sqrt(273) / 384292480000}         # orthogonal part along the unit fibre vector e = sin t v2 - cos t v_-3
DIAG547 = {'3p': (-267786421 / 58544557500, {10: 77 / 5475600, 14: 281211 / 189112352000, 16: 553 / 1591200000,
                                              18: 30233 / 56458242360}),
           '4': (-4797453339 / 2497901120000, {8: 3591 / 346112000, 12: 63 / 83200000, 14: 361557 / 931014656000,
                                                16: 34839 / 147097600000, 18: 2371131 / 8029616691200})}
bad, rec, count = [], {}, [0]


def check(name, ok, detail=''):
    count[0] += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ''))
    if not ok:
        bad.append(name)


def paper_g(t):                           # the paper's restriction of ||rho_6||^2 to the pyramid's line
    s2 = sin(t) ** 2
    return -125 / 132 * s2 ** 2 + 10 / 11 * s2 + 3 / 77


def curve(t, d):
    return sqrt(7 / d) * (cos(t) * E(2) + sin(t) * E(-3))


def tangent(t, d):                        # unit section, d/dt of the curve, oriented by increasing t
    return sqrt(7 / d) * (-sin(t) * E(2) + cos(t) * E(-3))


def quartic(u, P, d):
    B = block(u, P)
    wN, res = fibre(product_level(B, B, B, 3), P, d)
    return pair(u, wN, d).real, wN, res


def tangential(u, w):
    return w - np.vdot(u, w) / np.vdot(u, u) * u


def hess_form(u, e, P, d, conj_on=True):
    """Re<e, DN_Phi[e]> - Q, for unit block sections u and e given by fibre vectors."""
    w, _ = fibre(DN(block(u, P), [block(e, P)], conj_term_on=conj_on), P, d)
    return pair(e, w, d).real - quartic(u, P, d)[0]


def d2_along(u, e, P, d, h=1e-2):
    f = lambda s: quartic(np.cos(s) * u + np.sin(s) * e, P, d)[0]
    return (-f(2 * h) + 16 * f(h) - 30 * f(0) + 16 * f(-h) - f(-2 * h)) / (12 * h * h)


rng = np.random.default_rng(8110)
for s in ('3p', '4'):
    P, d, w6 = PROJ[s], DIM[s], W6[s]
    print(f'\n=== sector {s} (d = {d}, w6 = {w6:.12f}) ===')
    r = rec[s] = {}

    # Parents: the paper's restriction and criticality.
    errs = [abs(quartic(curve(t, d), P, d)[0] - (1 + w6 * paper_g(t))) for t in (T0, 0.3, 1.1)]
    check('parent: Q along the line = 1 + w6 * (paper restriction), three points', max(errs) < 1e-12, f'max err {max(errs):.1e}')
    u0, et = curve(T0, d), tangent(T0, d)
    Q0, wN0, res0 = quartic(u0, P, d)
    crit = np.linalg.norm(tangential(u0, wN0))
    check('parent: the pyramid is critical, Pi_T N(Phi) = 0', crit < 1e-12, f'{crit:.1e}; block identification {res0:.1e}')
    check('parent: Q at the pyramid = 1 + w6 * 9/35', abs(Q0 - (1 + w6 * 9 / 35)) < 1e-12, f'Q = {Q0:.15f}')

    # Bridge, first order at the control point: Re<e, N(Phi)> = c1 * dQ/dt, dQ/dt = w6 * 115 sqrt(3)/528 from the paper.
    uc, ec = curve(TC, d), tangent(TC, d)
    Qc, wNc, _ = quartic(uc, P, d)
    c1 = pair(ec, wNc, d).real / (w6 * 115 * sqrt(3) / 528)
    check('bridge route 1, first order: Re<e, N(Phi)> = (1/4) dQ/dt at the control point', abs(c1 - 0.25) < 1e-12, f'c1 = {c1:.15f}')

    # Bridge, second order at the pyramid along e_t: Hess form = c2 * d2Q/dt2, d2Q/dt2 = w6 * (-104/55) from the paper.
    d2_num, d2_pap = d2_along(u0, et, P, d), w6 * (-104 / 55)
    check('two routes to d2Q/dt2 at the pyramid: finite difference = w6 * (-104/55)', abs(d2_num - d2_pap) < 1e-7,
          f'{d2_num:.12f} vs {d2_pap:.12f}')
    Ht = hess_form(u0, et, P, d)
    c2 = Ht / d2_pap
    check('bridge route 2, second order: Re<e, DN e> - Q = (1/4) d2Q/dt2 along e_t', abs(c2 - 0.25) < 1e-12, f'c2 = {c2:.15f}')
    Hperp, d2perp = hess_form(u0, 1j * et, P, d), d2_along(u0, 1j * et, P, d)
    check('slice: the direction i e_t inside the locus is null (rotation orbit)', abs(Hperp) < 1e-12 and abs(d2perp) < 1e-7,
          f'Hess form {Hperp:.1e}, d2Q {d2perp:.1e}')

    # Bridge route 3: at a pinned ray (the hexagon), random tangents.
    uh = sqrt(7 / d) * (E(3) + E(-3)) / sqrt(2)
    ch = []
    for _ in range(3):
        e = rng.normal(size=7) + 1j * rng.normal(size=7)
        e -= np.vdot(uh, e) / np.vdot(uh, uh) * uh
        e *= sqrt(7 / d) / np.linalg.norm(e)
        ch.append(hess_form(uh, e, P, d) / d2_along(uh, e, P, d))
    check('bridge route 3, at the hexagon, three random tangents: c = 1/4', max(abs(x - 0.25) for x in ch) < 1e-7,
          ', '.join(f'{x:.10f}' for x in ch))

    # The forcing Pi_T DN_Phi[xi] and its decomposition in the fixed kernel K_H = span{v2, v_-3}.
    Phi = block(u0, P)
    cub = cubic_levels(Phi)
    from m811_ops import LEVELS, forbidden_norms, xi_terms_supported
    fb = forbidden_norms(cub, s)
    check("support (F1): the cubic vanishes at every level outside the sector's support", max(fb.values()) < 1e-24,
          ', '.join(f'{n}: {v:.0e}' for n, v in fb.items()))
    xi = xi_terms_supported(cub, s)
    cub = {J: ts for J, ts in cub.items() if 2 * J in LEVELS[s]}
    wF, resF = fibre(DN(Phi, xi), P, d)
    tanF = tangential(u0, wF)
    z = pair(et, tanF, d)                 # coefficient along the unit SECTION e_t
    outside = np.linalg.norm(tanF - z * et)
    check('forcing lies in the fixed kernel: no component off span{v2, v_-3}', outside < 1e-12, f'{outside:.1e}')
    check('forcing has no component along the orbit direction i e_t', abs(z.imag) < 1e-13, f'{z.imag:.1e}')
    # M8.10's audit gives its coefficient along the UNIT FIBRE vector e, and |e_t|_fibre = sqrt(7/d). First pass,
    # recorded: comparing z with the audit directly FAILED in both sectors with ratio exactly sqrt(d/7) (0.654654 in
    # 3', 0.755929 in 4), the section-versus-fibre normalization F1's bridge gate exists to catch.
    e_fib = -et / np.linalg.norm(et)
    cf = np.vdot(e_fib, tanF)
    check("parent (#547 D2): the forcing reproduces M8.10's audited residual, in the audit's fibre basis (e = -e_t)",
          abs(cf.real - AUDIT[s]) < 1e-14 and abs(cf.imag) < 1e-14, f'{cf.real:+.15e} vs {AUDIT[s]:+.15e}')

    # The tilt: (Pi_T DN_Phi - Q) v = -Pi_T DN_Phi[xi] on the one-dimensional slice.
    beta = -z.real / Ht
    v = beta * et
    wT, _ = fibre(DN(Phi, xi + [block(v, P)]), P, d)
    resid_with = np.linalg.norm(tangential(u0, wT) - Q0 * v)
    resid_without = np.linalg.norm(tanF)
    check('a^5 tangential residual vanishes with the tilt', resid_with < 1e-12, f'{resid_with:.1e}')
    check("arm: without the tilt the residual is M8.10's nonzero value (fibre norm)",
          resid_without > 1e-5 and abs(resid_without - abs(AUDIT[s])) < 1e-14, f'{resid_without:.15e}')

    # lambda_4 by M8.10's two routes, and with the tilt included.
    lv = {2 * J: (norm2(ts) if ts else 0.0) for J, ts in cub.items()}
    r1 = -3 * sum(val / den(n) for n, val in lv.items() if n != 6)
    r2 = pair(u0, wF, d)
    r2v = pair(u0, wT, d)
    check('lambda_4, two routes agree (level sum and direct pairing)', abs(r2 - r1) < 1e-12 and abs(r2.imag) < 1e-13, f'{r1:.15e}')
    check('lambda_4 lemma: including the tilt does not change it', abs(r2v - r2) < 1e-13, f'difference {abs(r2v - r2):.1e}')

    # Consistency with #547's blind R5 diagnostics (asked at go, ungraded): reproduced, not used as parents.
    lam547, xi547 = DIAG547[s]
    mine = {n: val / den(n) ** 2 for n, val in lv.items() if n != 6 and val > 1e-14}
    xi_err = max(abs(mine[n] / xi547[n] - 1) for n in xi547) if set(mine) == set(xi547) else float('inf')
    check("consistency with #547's blind R5 diagnostic (item 4, ungraded): xi level norms, same support", xi_err < 1e-11,
          f'max rel err {xi_err:.1e}')
    check("consistency with #547's blind R5 diagnostic (item 7, ungraded): lambda_4", abs(r1 / lam547 - 1) < 1e-12,
          f'{lam547:.15e}')

    # The lemma's gate, armed: at the non-critical control point the tilt-type pairing is nonzero.
    Lc = pair(uc, fibre(DN(block(uc, P), [block(ec, P)]), P, d)[0], d).real
    check('arm: at sin^2 t = 1/4 (not critical) Re<Phi, DN_Phi[e_t]> is nonzero', abs(Lc) > 1e-3, f'{Lc:.6e}')

    # Arms on the bridge.
    miss = abs(Ht / (-104 / 55) - 0.25)
    check('arm: comparing without the w6 conversion fails', miss > 1e-3, f'off by {miss:.3e}')
    c_noconj = hess_form(u0, et, P, d, conj_on=False) / d2_pap
    check('arm: dropping the conjugate term of DN moves c off 1/4', abs(c_noconj - 0.25) > 1e-3, f'c = {c_noconj:.6f}')

    # Where the forcing comes from: its component along e_t, split by the level of xi (a decomposition, not a gate:
    # DN is linear in xi, so the parts reassemble identically and a reassembly check could not fail).
    per_level = {}
    for J in sorted(cub):
        xiJ = [t for t in xi if t[0] == J]
        if J != 3 and xiJ:
            per_level[2 * J] = pair(et, tangential(u0, fibre(DN(Phi, xiJ), P, d)[0]), d).real
    print('  forcing along e_t by level of xi: ' + ', '.join(f'{n}: {val:+.6e}' for n, val in per_level.items()))

    support = {n: val for n, val in lv.items() if n != 6 and val > 1e-14}
    print('  level norms ||P_n N||^2: ' + ', '.join(f'{n}: {val:.12e}' for n, val in support.items()))
    print(f'  Q = {Q0:.15f}   Hess form along e_t = {Ht:.15e}   forcing along e_t = {z.real:+.15e}')
    print(f'  tilt coefficient along e_t (per g, a^3 term, unit section) = {beta:+.15e}   ||v||^2 = {beta ** 2:.15e}')
    print(f'  lambda_4/g^2 = {r1:.15e}')
    r.update(Q=Q0, bridge_c=[c1, c2] + ch, hess_along_et=Ht, forcing_along_et_section=z.real, forcing_fibre_basis=cf.real,
             tilt_beta=beta, tilt_norm2=beta ** 2, lambda4=r1, forcing_by_level=per_level, level_norms_N=support,
             xi_level_norms={n: val / den(n) ** 2 for n, val in support.items()})

os.makedirs('out', exist_ok=True)         # clean-directory runs (F1 found the write crashed without a pre-existing out/)
json.dump(rec, open('out/pyramid.json', 'w'), indent=1)
print(f'\n  {count[0]} checks; ' + ('ALL PASS, ALL ARMS FIRE' if not bad else f'FAILURES: {bad}'))
sys.exit(1 if bad else 0)
