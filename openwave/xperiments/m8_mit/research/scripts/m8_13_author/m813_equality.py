"""The equality set of the nematic bound at spin 3, exactly: where does lambda_max(e1 fx^2 + e2 fy^2 + e3 fz^2) reach
15/sqrt(6) over unit traceless e, and on which states is TrNbar^2 = 171/2 attained.

Author-side check for the M8.13 pre-registration, disclosed there. It re-derives the equality set of S0_S3_MAXIMUM.md
step 3 by solving each block's equality condition exactly on the constraint ellipse, rather than reading it off a
numerical sweep, and checks the top eigenspace at each equality point. It does not grade the note's own tracing.
Usage: python3 m813_equality.py
"""
import sys

import sympy as sp

fails = []


def gate(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)
    if not ok:
        fails.append(name)


a, b, lam = sp.symbols('a b lambda', real=True)
R6, top = sp.sqrt(6), 15 / sp.sqrt(6)
ellipse = sp.Rational(2, 3) * a**2 + 8 * b**2 - 1

# spin-3 operators in the basis m = 3..-3, and the operator a(fz^2 - 4) + b(f+^2 + f-^2) of step 3
ms = list(range(3, -4, -1))
Jz = sp.diag(*[sp.Integer(m) for m in ms])
Jp = sp.zeros(7, 7)
for i, m in enumerate(ms[1:], start=1):
    Jp[i - 1, i] = sp.sqrt(12 - m * (m + 1))
Jm = Jp.T
Op = a * (Jz**2 - 4 * sp.eye(7)) + b * (Jp**2 + Jm**2)

print('(1) the three blocks, and where each reaches 15/sqrt(6)')
M1 = sp.Matrix([[5 * a, sp.sqrt(60) * b], [sp.sqrt(60) * b, -3 * a + 12 * b]])
M2 = M1.subs(b, -b)
M3 = sp.Matrix([[0, sp.sqrt(240) * b], [sp.sqrt(240) * b, -4 * a]])
cp_full = sp.factor(sp.expand((lam * sp.eye(7) - Op).det()))
cp_blocks = sp.factor(sp.expand(lam * (lam * sp.eye(2) - M1).det() * (lam * sp.eye(2) - M2).det() * (lam * sp.eye(2) - M3).det()))
gate('the 7x7 operator splits into a zero eigenvalue and the three blocks, exactly', sp.expand(cp_full - cp_blocks) == 0)

eq_points = {}
for name, M in (('M1', M1), ('M2', M2), ('M3', M3)):
    # Completeness by elimination, not by a solver: 15/sqrt(6) is an eigenvalue iff det(15/sqrt(6) I - M) = 0, a conic.
    # With the ellipse, Bezout bounds the common points by 4. The resultant in b is a polynomial in a of degree <= 4
    # whose real roots are found exactly by factoring over Q(sqrt 6); each root is then completed by solving for b.
    det_top = sp.expand((top * sp.eye(2) - M).det())
    res = sp.Poly(sp.resultant(det_top, ellipse, b), a)
    fac = sp.factor_list(res.as_expr(), extension=sp.sqrt(6))
    roots_a = []
    for f, _ in fac[1]:
        pf = sp.Poly(f, a)
        if pf.degree() == 1:
            roots_a.append(sp.simplify(-pf.all_coeffs()[1] / pf.all_coeffs()[0]))
        else:
            roots_a += [r for r in sp.Poly(f, a).all_roots() if r.is_real]
    cand = []
    for ra in roots_a:
        for rb in sp.Poly(sp.expand(ellipse.subs(a, ra)), b).all_roots():
            if rb.is_real and sp.simplify(det_top.subs({a: ra, b: rb})) == 0:
                cand.append((sp.radsimp(sp.simplify(ra)), sp.radsimp(sp.simplify(rb))))  # exact roots, never a numerical guess
    cand = sorted(set(cand), key=lambda p: (float(p[0]), float(p[1])))
    gate(f'{name}: the resultant is not zero, so the conic and the ellipse share no component', not res.is_zero)
    gate(f'{name}: two conics with no common component meet in at most 4 points (Bezout); the resultant has degree {res.degree()}',
         not res.is_zero and res.degree() <= 4)
    gate(f'{name}: {len(cand)} common points found, within the Bezout bound', len(cand) <= 4)
    # keep a point only where 15/sqrt(6) is the larger eigenvalue of the block, not the smaller
    pts, undecided = [], []
    for (av, bv) in cand:
        tr = sp.simplify((M.subs({a: av, b: bv})).trace())
        other = sp.simplify(tr - top)          # the block's other eigenvalue
        margin = sp.radsimp(sp.simplify(top - other))
        if margin.is_nonnegative is True:      # decided symbolically, never by a float
            pts.append((av, bv))
        elif margin.is_nonnegative is None:
            undecided.append((av, bv))
    gate(f'{name}: every candidate\'s larger-eigenvalue test is decided symbolically', not undecided)
    eq_points[name] = pts
    print(f'    {name}: reaches 15/sqrt(6) exactly at {pts}')
union = sorted({p for v in eq_points.values() for p in v}, key=lambda p: (float(p[0]), float(p[1])))
print(f'    union over the blocks: {union}')
want = sorted({(R6 / 2, sp.Integer(0)), (-R6 / 4, R6 / 8), (-R6 / 4, -R6 / 8)}, key=lambda p: (float(p[0]), float(p[1])))
gate('the equality set on the ellipse is exactly the three points (sqrt6/2, 0), (-sqrt6/4, +-sqrt6/8)',
     [tuple(sp.simplify(x - y) for x, y in zip(p, q)) for p, q in zip(union, want)] == [(0, 0)] * 3 and len(union) == 3)
gate('M1 reaches it at (sqrt6/2, 0) and (-sqrt6/4, sqrt6/8), the cases b = 0, a > 0 and a = -2b, b > 0',
     set(eq_points['M1']) == {(R6 / 2, 0), (-R6 / 4, R6 / 8)})
gate('M2 mirrors M1 under b -> -b', set(eq_points['M2']) == {(R6 / 2, 0), (-R6 / 4, -R6 / 8)})
gate('M3 reaches it only at a = -sqrt6/4, on both points of the ellipse there', set(eq_points['M3']) == {(-R6 / 4, R6 / 8), (-R6 / 4, -R6 / 8)})
gate('every equality point is reached by exactly two blocks, so no single block is load-bearing for the set',
     all(sum(p in v for v in eq_points.values()) == 2 for p in union))

print('(2) back to e, and the top eigenspace at each point')
for (av, bv) in union:
    e3 = sp.Rational(2, 3) * av
    e1 = sp.simplify((-e3 + 4 * bv) / 2)
    e2 = sp.simplify((-e3 - 4 * bv) / 2)
    e = [sp.simplify(x * R6) for x in (e1, e2, e3)]
    Ov = Op.subs({a: av, b: bv})
    evs = Ov.eigenvals()
    mult = evs.get(top, 0) + sum(k for v, k in evs.items() if sp.simplify(v - top) == 0 and v != top)
    null = (top * sp.eye(7) - Ov).nullspace()
    print(f'    (a, b) = ({av}, {bv}):  sqrt6 * e = {e};  15/sqrt6 has multiplicity {len(null)}')
    gate(f'      a permutation of (2, -1, -1), and a two-dimensional top eigenspace',
         sorted(e) == [-1, -1, 2] and len(null) == 2)
    if bv == 0:
        support = sorted({ms[i] for vec in null for i in range(7) if vec[i] != 0})
        gate('      at the z-axis point the top eigenspace is span{v3, v-3}', support == [-3, 3])

print('(3) arms')
M3p = M3 + sp.Rational(1, 50) * sp.Matrix([[1, 0], [0, 0]])
det_p = sp.expand((top * sp.eye(2) - M3p).det())
sols_p = sp.solve([det_p, ellipse], [a, b], dict=True)
gate('arm: a perturbed M3 still has equality points, and they move off a = -sqrt6/4',
     len(sols_p) > 0 and all(sp.simplify(s[a] + R6 / 4) != 0 for s in sols_p))
res0 = sp.Poly(sp.resultant(sp.expand(ellipse * (a + 1)), ellipse, b), a)
gate('arm: a conic sharing a component with the ellipse gives a zero resultant, which the Bezout guard refuses', res0.is_zero)
gate('arm: the comparison would catch a missing point', union[:2] != want)
print('ALL PASS' if not fails else f'{len(fails)} FAILED: {fails}')
sys.exit(1 if fails else 0)
