"""The negative control, computed rather than asserted: the point the worklist hands the rooms must really be
non-critical, and the criticality test must really fail there.

It reuses the step-3 machinery by running s1_H.py, so the control is measured with the same N, Hessian and
projections that produced the census, not with a second implementation of them.
Arms: at a critical point both quantities vanish, which is the green half of the pair.
Usage: python3 s1_control.py
"""
import runpy
import sys

import sympy as sp

ns = runpy.run_path('s1_H.py')
X, ms, idx, Nexpr, HessN = ns['X'], ns['ms'], ns['idx'], ns['Nexpr'], ns['HessN']
realvec, unit, Jx, Jy, Jz = ns['realvec'], ns['unit'], ns['Jx'], ns['Jy'], ns['Jz']
print('\n=== s1_control.py ===', flush=True)
fails = []


def gate(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)
    if not ok:
        fails.append(name)


def at(u):
    """(924 r6, the tangential gradient of r6 on the sphere, and per generator: |M d|, exactly, with grad N . d)."""
    xs = realvec(u)
    sub = {X[i]: xs[i] for i in range(14)}
    Nval = sp.simplify(Nexpr.subs(sub))
    gradN = sp.Matrix([sp.simplify(sp.diff(Nexpr, X[i]).subs(sub)) for i in range(14)])
    g = sp.Matrix([sp.simplify(c) for c in gradN - 4 * Nval * xs])  # degree 0, so this is already tangential
    M = sp.Matrix(14, 14, lambda a, b: sp.simplify(HessN[a, b].subs(sub) - (4 * Nval if a == b else 0)))
    per = []
    for name, d in (('i u', realvec(sp.I * u)), ('-i Jx u', realvec(-sp.I * (Jx * u))),
                    ('-i Jy u', realvec(-sp.I * (Jy * u))), ('-i Jz u', realvec(-sp.I * (Jz * u)))):
        sq = sp.simplify(((M * d).T * (M * d))[0, 0])
        per.append((name, sp.radsimp(sp.sqrt(sq)), sq, sp.simplify((gradN.T * d)[0, 0])))
    return sp.simplify(924 * Nval), sp.radsimp(sp.simplify(g.norm())), per


CONTROL = {3: sp.Integer(1), 1: sp.Integer(1), -2: sp.Integer(-1)}          # (v3 + v1 - v-2)/sqrt(3)
EX1 = {3: sp.Integer(2), 1: sp.Integer(1), -2: sp.Integer(-3)}              # worklist item 0
EX2 = {3: sp.Integer(1), 0: sp.I, -1: sp.Integer(2)}                        # worklist item 0

u_c = unit(CONTROL)
val, grad, per = at(u_c)
print(f'  control u = (v3 + v1 - v-2)/sqrt(3): 924*r6 = {val} = {float(val):.6f}', flush=True)
print(f'    tangential gradient norm = {grad} = {sp.N(grad, 20)}', flush=True)
for name, r, sq, dot in per:
    print(f'    |M d| at {name:8s} = {r} = {sp.N(r, 20)}   (grad N . d = {dot})', flush=True)
tot = sp.simplify(sum(sq for _, _, sq, _ in per))
print(f'    the four squares sum to {tot}', flush=True)
gate('  the control point is not critical: the tangential gradient is nonzero', grad != 0)
gate('  every orbit residual is nonzero, so the criticality test fails there', all(r != 0 for _, r, _, _ in per))
gate('  the residual at i u equals the gradient norm, as the phase identity gives',
     sp.simplify(per[0][1] - grad) == 0)
gate('  grad N is orthogonal to every orbit direction, so M and the second variation agree there',
     all(dot == 0 for _, _, _, dot in per))
gate('  the four squares sum to 17368/29403', sp.simplify(tot - sp.Rational(17368, 29403)) == 0)
gate('  the Jx and Jy squares sum to 11120/29403',
     sp.simplify(per[1][2] + per[2][2] - sp.Rational(11120, 29403)) == 0)
# the exact expressions N1 freezes, compared here rather than by string, since sympy prints an equal surd differently
FROZEN = {'i u': 2 * sp.sqrt(1002) / 297, '-i Jx u': 2 * sp.sqrt(312 * sp.sqrt(15) + 4170) / 297,
          '-i Jy u': 2 * sp.sqrt(4170 - 312 * sp.sqrt(15)) / 297, '-i Jz u': 2 * sp.sqrt(3684) / 297}
for name, r, _, _ in per:
    gate(f'  the residual at {name:8s} is N1\'s frozen {sp.sstr(FROZEN[name])}', sp.simplify(r - FROZEN[name]) == 0)
gate('  arm: a neighbouring surd is not equal to the frozen one',
     sp.simplify(2 * sp.sqrt(3685) / 297 - FROZEN['-i Jz u']) != 0)
gate('  its value is not a census value',
     val not in (1, 36, 225, 400, 288, 463, sp.Rational(1188, 5), sp.Rational(8800, 43)))

val0, grad0, per0 = at(unit({3: sp.Integer(1)}))
print(f'  arm, at the coherent state: 924*r6 = {val0}, gradient = {grad0}, residuals = {[str(r) for _, r, _, _ in per0]}', flush=True)
gate('  arm: at a critical point the gradient vanishes', grad0 == 0)
gate('  arm: at a critical point every orbit residual vanishes exactly', all(r == 0 for _, r, _, _ in per0))

for name, d in (('2v3 + v1 - 3v-2', EX1), ('v3 + i v0 + 2v-1', EX2)):
    v, g_, _ = at(unit(d))
    print(f'  worklist item 0 example {name}: 924*r6 = {v} = {float(v):.6f}, gradient = {float(g_):.3e}', flush=True)
    gate(f'    {name}: not critical either', g_ != 0)

print('ALL PASS' if not fails else f'{len(fails)} FAILED: {fails}')
sys.exit(1 if fails else 0)
