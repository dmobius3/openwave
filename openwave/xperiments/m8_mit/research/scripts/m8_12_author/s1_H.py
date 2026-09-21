"""S1 derivation, step 3 of 3: H, the exact transverse signature of r6 at the ten critical orbits of L0 and L1.

Object. For a unit critical u, H_u(e, e) = d^2/ds^2 r6(u cos s + e sin s) at s = 0, a self-adjoint form on
T_u = {e : Re<u, e> = 0} (13 real dimensions). Because r6 is homogeneous of degree 0 and u is critical, the full
gradient vanishes there, and on the sphere the form is the Euclidean Hessian of r6 at u. With ||u|| = 1 and
N(x) = ||rho_6||^2 that gives, on vectors orthogonal to u,

    H_u = Hess N(u) - 4 N(u) I .

The orbit directions O_u = span{i u, -i J_x u, -i J_y u, -i J_z u} are null; N_u is their orthogonal complement in
T_u, of real dimension 10 at the weight states and 9 elsewhere.

Frozen per orbit: the signature (n-, n0, n+) of H_u on N_u by exact inertia (primary); the exact characteristic
polynomial of the transverse operator (secondary); the eigenvalues with multiplicities.
Reading: the reduced Morse index is n- for g > 0 and n+ for g < 0.

Controls: M8.11's in-locus values at the pyramid and the prism, and the two global ends. Arms: dropping the -4N I
term must break the controls; at a non-critical point the orbit directions are not annihilated; the exact and
numerical signatures must agree.
"""
import itertools
import numpy as np
import sympy as sp
from sympy.physics.wigner import clebsch_gordan as CG

fails = []


def gate(name, ok):
    print(('PASS ' if ok else 'FAIL ') + str(name), flush=True)
    if not ok:
        fails.append(name)


ms = list(range(3, -4, -1))
idx = {m: i for i, m in enumerate(ms)}
X = sp.symbols('x0:14', real=True)
U = [X[i] + sp.I * X[i + 7] for i in range(7)]
CG6 = {(m1, m2): sp.nsimplify(CG(3, 3, 6, m1, m2, m1 + m2)) for m1 in ms for m2 in ms}
TU = [(-1) ** (3 - m) * sp.conjugate(U[idx[-m]]) for m in ms]
TU = [sp.expand(t.rewrite(sp.re)) if False else t for t in TU]
Nexpr = 0
for Q in range(-6, 7):
    c = sum(CG6[m1, Q - m1] * U[idx[m1]] * TU[idx[Q - m1]] for m1 in ms if abs(Q - m1) <= 3)
    Nexpr += sp.expand(sp.expand(c) * sp.conjugate(sp.expand(c)))
Nexpr = sp.expand(sp.simplify(sp.expand(Nexpr).rewrite(sp.re).expand()))
Nexpr = sp.expand(sp.re(Nexpr)) if Nexpr.has(sp.I) else sp.expand(Nexpr)
print('N built:', sp.count_ops(Nexpr), 'ops', flush=True)
HessN = sp.hessian(Nexpr, X)
print('Hessian built', flush=True)

Jz = sp.diag(*[sp.Integer(m) for m in ms])
Jp = sp.zeros(7, 7)
for m in ms[1:]:
    Jp[idx[m + 1], idx[m]] = sp.sqrt(3 * 4 - m * (m + 1))
Jx, Jy = (Jp + Jp.T) / 2, (Jp - Jp.T) / (2 * sp.I)


def realvec(u):
    return sp.Matrix([sp.re(c) for c in u] + [sp.im(c) for c in u])


def cplx(vec14):
    return sp.Matrix([vec14[i] + sp.I * vec14[i + 7] for i in range(7)])


def inertia(M):
    """Exact inertia of a symmetric matrix by symmetric elimination, with 2x2 blocks when a pivot vanishes."""
    M = sp.Matrix(M)
    n = M.rows
    neg = pos = zero = 0
    rows = list(range(n))
    while rows:
        d = [i for i in rows if sp.simplify(M[i, i]) != 0]
        if d:
            i = d[0]
            p = sp.simplify(M[i, i])
            pos, neg = (pos + 1, neg) if p > 0 else (pos, neg + 1)
            for a in rows:
                if a != i:
                    f = sp.simplify(M[a, i] / p)
                    for b in rows:
                        M[a, b] = sp.simplify(M[a, b] - f * M[i, b])
            rows.remove(i)
        else:
            off = [(i, k) for i in rows for k in rows if i < k and sp.simplify(M[i, k]) != 0]
            if not off:
                zero += len(rows)
                break
            i, k = off[0]
            pos += 1
            neg += 1                      # a 2x2 block with zero diagonal has determinant < 0
            for a in rows:
                if a not in (i, k):
                    f1 = sp.simplify(M[a, k] / M[i, k])
                    f2 = sp.simplify(M[a, i] / M[i, k])
                    for b in rows:
                        M[a, b] = sp.simplify(M[a, b] - f1 * M[i, b] - f2 * M[k, b])
            rows.remove(i)
            rows.remove(k)
    return neg, zero, pos


s5 = sp.sqrt(5)
ORB = {
    'coherent v3': {3: sp.Integer(1)},
    'v2': {2: sp.Integer(1)},
    'v1': {1: sp.Integer(1)},
    'zonal v0': {0: sp.Integer(1)},
    'octahedron': {2: sp.Integer(1), -2: sp.Integer(1)},
    'hexagon': {3: sp.Integer(1), -3: sp.Integer(1)},
    'pyramid': {3: sp.sqrt(12), -2: sp.sqrt(13)},
    'prism': {3: sp.Integer(1), 0: sp.sqrt(sp.Rational(23, 10)), -3: sp.Integer(1)},
    'C3 ray': {2: sp.Integer(1), -1: sp.Integer(2)},
    'D2 ray': {2: sp.Integer(1), 0: sp.I * sp.sqrt(sp.Rational(6, 5)), -2: sp.Integer(1)},
}


def unit(d):
    u = sp.Matrix([d.get(m, 0) for m in ms])
    return sp.simplify(u / sp.sqrt(sum(sp.Abs(c) ** 2 for c in u)))


def ortho(vecs):
    """Exact Gram-Schmidt that drops dependent vectors, which the orbit set is at the weight states."""
    out = []
    for w in vecs:
        w = sp.Matrix([sp.simplify(c) for c in w])
        for q in out:
            w = w - (q.T * w)[0, 0] * q
        w = sp.Matrix([sp.simplify(c) for c in w])
        n = sp.simplify(w.norm())
        if n != 0:
            out.append(sp.Matrix([sp.simplify(c / n) for c in w]))
    return out


def transverse(u):
    """H_u and an exact orthonormal basis B of N_u (columns, real 14-vectors)."""
    xs = realvec(u)
    sub = {X[i]: xs[i] for i in range(14)}
    Nval = sp.simplify(Nexpr.subs(sub))
    H = sp.Matrix(14, 14, lambda a, b: sp.simplify(HessN[a, b].subs(sub) - (4 * Nval if a == b else 0)))
    orbit = [realvec(sp.I * u)] + [realvec(-sp.I * (A * u)) for A in (Jx, Jy, Jz)]
    QB = ortho([xs] + orbit)
    P = sp.eye(14) - sum((q * q.T for q in QB), sp.zeros(14, 14))
    B = sp.Matrix.hstack(*ortho([P * sp.eye(14)[:, i] for i in range(14)]))
    return H, B, Nval, orbit


print('(H) the ten orbits')
rows = []
for name, d in ORB.items():
    u = unit(d)
    H, B, Nval, orbit = transverse(u)
    gate(f'  {name}: the orbit directions are null for H', all(sp.simplify((H * sp.Matrix(o)).norm()) == 0 for o in orbit))
    M = sp.simplify(B.T * H * B)
    k = M.rows
    neg, zero, pos = inertia(M)
    lam = sp.symbols('lam')
    cp = sp.factor(sp.simplify(M.charpoly(lam).as_expr()))
    ev = np.linalg.eigvalsh(np.array(M.evalf(30), dtype=float))
    nneg, nzer, npos = int(np.sum(ev < -1e-9)), int(np.sum(abs(ev) < 1e-9)), int(np.sum(ev > 1e-9))
    gate(f'  {name}: dim N_u = {k}; exact signature (n-, n0, n+) = ({neg}, {zero}, {pos}); numerical ({nneg}, {nzer}, {npos})',
         (neg, zero, pos) == (nneg, nzer, npos) and neg + zero + pos == k)
    print(f'    924*r6 = {sp.nsimplify(924 * Nval)}   charpoly = {cp}')
    print(f'    eigenvalues (30 digits): {np.round(ev, 6).tolist()}')
    if zero > 0:
        ker = M.nullspace()
        amb = [sp.Matrix([sp.simplify(c) for c in B * v_]) for v_ in ker]
        # the C4 line through this point is the flipped copy span{v-3, v1}, whose restriction has its vertex at this
        # endpoint; the first run of this gate pointed at v3, which is the same line at v-1, not at v1
        cand = [realvec(sp.Matrix([1 if mm == -3 else 0 for mm in ms])), realvec(sp.Matrix([sp.I if mm == -3 else 0 for mm in ms]))]
        proj = sp.Matrix.hstack(*[sp.Matrix([sp.simplify(c) for c in B * (B.T * cc)]) for cc in cand])
        both = sp.Matrix.hstack(*(amb + [proj[:, i] for i in range(proj.cols)]))
        gate(f'  {name}: the {zero}-dimensional kernel is the C4 line direction span(v-3, i v-3) (rank {both.rank()} = {len(amb)})',
             both.rank() == len(amb) == zero)
    rows.append((name, sp.nsimplify(924 * Nval), k, (neg, zero, pos), cp))

print('(P) parents and controls')
# M8.11 in-locus values: the pyramid along the C5 line, and the prism chart tangents
def second_along(u, w):
    """H(e, e) for e the unit tangent obtained by projecting the ambient direction w onto N_u, that is, orthogonal to
    the radial direction AND to the orbit directions. The first run projected out only the radial direction, which
    matches M8.11's tau_x but not its tau_y = i tau_x, and the prism's tau_y control caught it."""
    H, B, Nval, _ = transverse(u)
    e = sp.Matrix([sp.simplify(c) for c in B * (B.T * realvec(w))])
    n = sp.simplify(e.norm())
    e = sp.Matrix([sp.simplify(c / n) for c in e])
    return sp.simplify((e.T * H * e)[0, 0])


u_pyr = unit(ORB['pyramid'])
w_t = sp.Matrix([sp.Rational(1, 2) / sp.sqrt(sp.Rational(12, 25)) if m == 3 else (-sp.Rational(1, 2) / sp.sqrt(sp.Rational(13, 25)) if m == -2 else 0) for m in ms])
gate('  pyramid: H along the unit C5-line tangent is -104/55 (M8.11 T4)', sp.simplify(second_along(u_pyr, w_t) + sp.Rational(104, 55)) == 0)
u_pri = unit(ORB['prism'])
gate('  prism: H along the chart tangent v0 is 920/473 (M8.11 T4)', sp.simplify(second_along(u_pri, sp.Matrix([1 if m == 0 else 0 for m in ms])) - sp.Rational(920, 473)) == 0)
gate('  prism: H along the chart tangent i v0 is 8/11 (M8.11 T4)', sp.simplify(second_along(u_pri, sp.Matrix([sp.I if m == 0 else 0 for m in ms])) - sp.Rational(8, 11)) == 0)
w6 = {'3p': sp.Rational(28, 39), '4': sp.Rational(21, 52)}
gate('  the sector bridge gives M8.11 L_T at the pyramid: -56/165 and -21/110',
     sp.simplify(w6['3p'] / 4 * sp.Rational(-104, 55) + sp.Rational(56, 165)) == 0 and sp.simplify(w6['4'] / 4 * sp.Rational(-104, 55) + sp.Rational(21, 110)) == 0)
ends = {r[0]: r[3] for r in rows}
gate('  the two ends: n- = 0 at the coherent orbit and n+ = 0 at the hexagon', ends['coherent v3'][0] == 0 and ends['hexagon'][2] == 0)
gate('  the zonal is a saddle, as step 2 predicted (n- > 0 and n+ > 0)', ends['zonal v0'][0] > 0 and ends['zonal v0'][2] > 0)

print('(A) arms')
u_hex = unit(ORB['hexagon'])
xs = realvec(u_hex)
Hbad = sp.simplify(HessN.subs({X[i]: xs[i] for i in range(14)}))          # the -4 N I term dropped
_, Bh, _, _ = transverse(u_hex)
Mbad = sp.simplify(Bh.T * Hbad * Bh)
gate('  arm: dropping the -4 N I term changes the hexagon signature', inertia(Mbad) != ends['hexagon'])
u_nc = unit({2: sp.Integer(1), -1: sp.Integer(1)})                        # a non-critical point on the C3 line
xs = realvec(u_nc)
Hnc = sp.simplify(HessN.subs({X[i]: xs[i] for i in range(14)}) - 4 * sp.simplify(Nexpr.subs({X[i]: xs[i] for i in range(14)})) * sp.eye(14))
gate('  arm: at a non-critical point the orbit directions are not annihilated',
     any(sp.simplify((Hnc * sp.Matrix(realvec(-sp.I * (A * u_nc)))).norm()) != 0 for A in (Jx, Jy, Jz)))

print('(S) the census')
print(f'  {"orbit":14s} {"924*r6":>10s}  dim N_u  (n-, n0, n+)   index for g > 0   index for g < 0')
for name, val, k, sig, _ in rows:
    print(f'  {name:14s} {str(val):>10s}  {k:7d}  {str(sig):>12s}   {sig[0]:>13d}   {sig[2]:>14d}')
par = [r for r in rows if r[0] in ('coherent v3', 'v2', 'v1', 'zonal v0')]
print('  parity diagnostic (not a gate): indices at the weight states for g > 0:', [(r[0], r[3][0]) for r in par])
print(f'{len(fails)} FAILED: {fails}' if fails else 'ALL PASS')
