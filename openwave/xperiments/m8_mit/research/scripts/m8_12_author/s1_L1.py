"""S1 derivation, step 2 of 3: L1, the exact critical set of r6 on each of the seven fixed lines of L0, then the
union modulo rotation (author-side; no Hessian of the full problem; the in-line second derivative is used only for
the Poincare-Hopf index and is not a frozen claim).

Per line: the exact restriction of r6 (sympy), its complete critical set, exact coordinates and values.
  Cyclic lines u = sqrt(s) v_a + sqrt(1-s) e^{i phi} v_b: r6 is a quadratic q(s), independent of phi; critical points
  are the endpoints (weight states) and the vertex if it lies in (0, 1), a circle that is one orbit.
  Dihedral lines, chart u = b1 + z b2 (z = x + iy) plus z = infinity: all solutions of grad f = 0, exactly; the
  Poincare-Hopf indices must sum to chi(S^2) = 2 (indices from the chart Hessian's determinant sign; a degenerate
  point would get its winding number).
Controls: the paper's D3 chart formula and pyramid-line quadratic; the known values. Arms: the D3 formula fails with
u (x) u in place of u (x) Theta u; dropping a point breaks the index sum; a non-critical line point fails the
full-space gradient check.
Union: points with different exact values are different orbits; ties are resolved by a rotation search (same
orbit) or a rotation-invariant separator (different orbits): the invariant triple (|f|^2, |a00|^2, TrN^2) or the
isotropy orders.
"""
import itertools
import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import minimize
from sympy.physics.wigner import clebsch_gordan as CG

fails = []


def gate(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)
    if not ok:
        fails.append(name)


ms = list(range(3, -4, -1))
x, y, s_, t_ = sp.symbols('x y s t', real=True)
CG6 = {(m1, m2): sp.nsimplify(CG(3, 3, 6, m1, m2, m1 + m2)) for m1 in ms for m2 in ms}


def theta(u):
    return {m: (-1) ** (3 - m) * sp.conjugate(u.get(-m, 0)) for m in ms}


def rho6_sq(u, tr=True):
    """||[u (x) Theta u]_6||^2 for u given as {m: coefficient}, exact; tr=False gives u (x) u (the arm)."""
    w = theta(u) if tr else u
    tot = 0
    for Q in range(-6, 7):
        c = sum(CG6[m1, Q - m1] * u.get(m1, 0) * w.get(Q - m1, 0) for m1 in ms if abs(Q - m1) <= 3)
        tot += sp.expand(c * sp.conjugate(c))
    return sp.simplify(sp.expand(tot))


def norm_sq(u):
    return sp.expand(sum(u[m] * sp.conjugate(u[m]) for m in u))


# numeric side, for the full-space checks and the union
J3 = 3
idx = {m: i for i, m in enumerate(ms)}
Jz = np.diag([float(m) for m in ms]).astype(complex)
Jp = np.zeros((7, 7), complex)
for m in ms[1:]:
    Jp[idx[m + 1], idx[m]] = np.sqrt(J3 * (J3 + 1) - m * (m + 1))
Jx, Jy = (Jp + Jp.conj().T) / 2, (Jp - Jp.conj().T) / 2j
CGN = {(m1, m2): float(CG6[m1, m2]) for m1 in ms for m2 in ms}


def r6n(u):
    u = u / np.linalg.norm(u)
    tu = np.array([(-1) ** (3 - m) * np.conj(u[idx[-m]]) for m in ms])
    return float(sum(abs(sum(CGN[m1, Q - m1] * u[idx[m1]] * tu[idx[Q - m1]] for m1 in ms if abs(Q - m1) <= 3)) ** 2 for Q in range(-6, 7)))


def grad_norm(u, h=1e-6):
    u = u / np.linalg.norm(u)
    X = np.concatenate([u.real, u.imag])
    g = np.array([(r6n((X + h * e)[:7] + 1j * (X + h * e)[7:]) - r6n((X - h * e)[:7] + 1j * (X - h * e)[7:])) / (2 * h) for e in np.eye(14)])
    for d in (X, np.concatenate([-u.imag, u.real])):
        g -= g @ d / (d @ d) * d
    return float(np.linalg.norm(g))


def invariants(u):
    u = u / np.linalg.norm(u)
    f2 = sum(abs(np.vdot(u, A @ u)) ** 2 for A in (Jx, Jy, Jz))
    a00 = abs(sum((-1) ** (3 - m) / np.sqrt(7) * u[idx[m]] * u[idx[-m]] for m in ms)) ** 2
    N = np.array([[np.vdot(u, (A @ B + B @ A) @ u).real / 2 for B in (Jx, Jy, Jz)] for A in (Jx, Jy, Jz)])
    return np.array([f2, a00, float(np.sum(N * N))])


Dw = lambda w: expm(-1j * (w[0] * Jx + w[1] * Jy + w[2] * Jz))
rng = np.random.default_rng(7)


def same_orbit(u1, u2, starts=64):
    u1, u2 = u1 / np.linalg.norm(u1), u2 / np.linalg.norm(u2)
    f = lambda w: 1 - abs(np.vdot(u2, Dw(w) @ u1)) ** 2
    best = np.inf
    for _ in range(starts):
        r = minimize(f, rng.normal(size=3) * 2, method='BFGS', options={'gtol': 1e-14})
        best = min(best, r.fun)
        if best < 1e-13:
            break
    return best < 1e-12


def rot_orders(u):
    """Orders k in 2..6 of rotations fixing the ray, by a frame-independent axis search."""
    u = u / np.linalg.norm(u)
    grid = [(np.arccos(1 - 2 * (i + 0.5) / 1500), np.pi * (1 + 5 ** 0.5) * i) for i in range(1500)]
    ax = lambda th, ph: np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    out = []
    for k in range(2, 7):
        d = lambda g: 1 - abs(np.vdot(u, expm(-1j * 2 * np.pi / k * np.tensordot(ax(*g), np.array([Jx, Jy, Jz]), 1)) @ u))
        best = sorted((d(g), g) for g in grid)[:8]
        if any(minimize(d, np.array(g), method='Nelder-Mead', options={'xatol': 1e-12, 'fatol': 1e-15, 'maxiter': 3000}).fun < 1e-10 for _, g in best):
            out.append(k)
    return out


found = []          # (line, label, numeric unit vector, exact value of r6)


def add(line, label, u_exact, val):
    un = np.array([complex(sp.N(u_exact.get(m, 0), 30)) for m in ms])
    found.append((line, label, un / np.linalg.norm(un), sp.nsimplify(val)))


print('(C) the five cyclic lines: r6 = q(s), s = |coefficient of the first vector|^2')
CYC = {'C3 {v2,v-1}': (2, -1), 'C4 {v3,v-1}': (3, -1), 'C4 {v2,v-2}': (2, -2), 'C5 {v3,v-2}': (3, -2), 'C6 {v3,v-3}': (3, -3)}
ph = sp.symbols('phi', real=True)
r_ = sp.symbols('r', positive=True)


def cyc_q(a, b):
    """q(s) exactly: u = v_a + r e^{i phi} v_b with r > 0, so no square roots; then r^2 = (1 - s)/s, s = |c_a|^2."""
    u = {a: sp.Integer(1), b: r_ * sp.exp(sp.I * ph)}
    F = sp.simplify(rho6_sq(u) / norm_sq(u) ** 2)
    return F, sp.factor(sp.simplify(F.subs(r_, sp.sqrt((1 - s_) / s_))))


for name, (a, b) in CYC.items():
    F, q = cyc_q(a, b)
    gate(f'  {name}: r6 does not depend on the relative phase', sp.simplify(sp.diff(F, ph)) == 0)
    poly = sp.Poly(sp.expand(sp.cancel(q)), s_)
    gate(f'  {name}: r6 along the line is a polynomial of degree <= 2 in s', poly.degree() <= 2)
    A2, A1, A0 = [poly.coeff_monomial(s_ ** k) for k in (2, 1, 0)]
    print(f'    q(s) = {A2}*s**2 + {A1}*s + {A0}')
    add(name, f'v{a} (s = 1)', {a: 1}, A2 + A1 + A0)
    add(name, f'v{b} (s = 0)', {b: 1}, A0)
    if A2 != 0:
        sv = -A1 / (2 * A2)
        inside = 0 < sv < 1
        print(f'    vertex s* = {sv}  ({"interior: one critical circle, one orbit" if inside else "not interior: no interior critical point"})')
        if inside:
            add(name, f'interior circle, s* = {sv}', {a: sp.sqrt(sv), b: sp.sqrt(1 - sv)}, A2 * sv ** 2 + A1 * sv + A0)
    else:
        gate(f'  {name}: q is not constant along the line', A1 != 0)
q_c5 = sp.expand(sp.cancel(cyc_q(3, -2)[1]))
gate('  control: on the C5 line q(s) is the paper pyramid-line quadratic -125/132 s^2 + 10/11 s + 3/77 (up to the flip)',
     sp.simplify(q_c5 - (-sp.Rational(125, 132) * s_ ** 2 + sp.Rational(10, 11) * s_ + sp.Rational(3, 77))) == 0)

print('(D) the two dihedral lines, chart u = b1 + z b2, z = x + i y, and z = infinity')
z = x + sp.I * y
DIH = {'D3 {v3+v-3,v0}': ({3: 1, -3: 1}, {0: 1}), 'D2 {v2+v-2,v0}': ({2: 1, -2: 1}, {0: 1})}
for name, (b1, b2) in DIH.items():
    u = {m: b1.get(m, 0) + z * b2.get(m, 0) for m in ms if b1.get(m, 0) != 0 or b2.get(m, 0) != 0}
    f = sp.factor(sp.simplify(rho6_sq(u) / norm_sq(u) ** 2))
    print(f'    f(x, y) = {f}')
    if name.startswith('D3'):
        paper = (100 * (x ** 2 + y ** 2) ** 2 - 20 * x ** 2 + 148 * y ** 2 + 463) / (231 * (x ** 2 + y ** 2 + 2) ** 2)
        gate('  control: the D3 chart formula is the paper Section 5.8 formula', sp.simplify(f - paper) == 0)
        fbad = sp.simplify(rho6_sq(u, tr=False) / norm_sq(u) ** 2)
        gate('  arm: with u (x) u in place of u (x) Theta u the formula differs', sp.simplify(fbad - paper) != 0)
    num = sp.numer(sp.together(f))
    fx, fy = sp.numer(sp.together(sp.diff(f, x))), sp.numer(sp.together(sp.diff(f, y)))
    sols = sp.solve([fx, fy], [x, y], dict=True)
    real = []
    for so in sols:
        xv, yv = sp.nsimplify(sp.simplify(so[x])), sp.nsimplify(sp.simplify(so[y]))
        if xv.is_real and yv.is_real:
            real.append((xv, yv))
    real = sorted(set(real), key=lambda p: (float(p[0]), float(p[1])))
    # second exact route to completeness: f is even in x and in y, so fx = x*A(X, Y) and fy = y*B(X, Y) with X = x^2,
    # Y = y^2; the four cases are solved by elimination (resultants) and real-root counts, independently of solve()
    X_, Y_ = sp.symbols('X Y', nonnegative=True)
    gate('  f is even in x and in y separately', sp.simplify(f.subs(x, -x) - f) == 0 and sp.simplify(f.subs(y, -y) - f) == 0)
    A = sp.expand(sp.cancel(fx / x)).subs({x: sp.sqrt(X_), y: sp.sqrt(Y_)})
    B = sp.expand(sp.cancel(fy / y)).subs({x: sp.sqrt(X_), y: sp.sqrt(Y_)})
    A, B = sp.expand(A), sp.expand(B)
    gate('  fx/x and fy/y are polynomials in X = x^2, Y = y^2', A.free_symbols <= {X_, Y_} and B.free_symbols <= {X_, Y_} and A.is_polynomial(X_, Y_) and B.is_polynomial(X_, Y_))
    route2 = {(sp.Integer(0), sp.Integer(0))}
    for Xr in sp.Poly(A.subs(Y_, 0), X_).real_roots():                          # y = 0, A(X, 0) = 0, X > 0
        if Xr > 0:
            route2 |= {(sp.sqrt(Xr), sp.Integer(0)), (-sp.sqrt(Xr), sp.Integer(0))}
    for Yr in sp.Poly(B.subs(X_, 0), Y_).real_roots():                          # x = 0, B(0, Y) = 0, Y > 0
        if Yr > 0:
            route2 |= {(sp.Integer(0), sp.sqrt(Yr)), (sp.Integer(0), -sp.sqrt(Yr))}
    res = sp.resultant(sp.Poly(A, X_, Y_).as_expr(), sp.Poly(B, X_, Y_).as_expr(), Y_)   # off-axis: A = B = 0, X, Y > 0
    offaxis = []
    if res != 0:
        for Xr in sp.Poly(res, X_).real_roots():
            if Xr > 0:
                g = sp.gcd(sp.Poly(A.subs(X_, Xr), Y_), sp.Poly(B.subs(X_, Xr), Y_))
                for Yr in (sp.Poly(g, Y_).real_roots() if g.degree() > 0 else []):
                    if Yr > 0:
                        offaxis.append((Xr, Yr))
    gate(f'  off-axis case: the resultant is not identically zero and gives {len(offaxis)} off-axis critical points', res != 0)
    for Xr, Yr in offaxis:
        for sx, sy in itertools.product((1, -1), repeat=2):
            route2.add((sx * sp.sqrt(Xr), sy * sp.sqrt(Yr)))
    same = {(sp.nsimplify(a_), sp.nsimplify(b_)) for a_, b_ in route2} == {(sp.nsimplify(a_), sp.nsimplify(b_)) for a_, b_ in real}
    gate(f'  the elimination route finds the same {len(route2)} chart critical points as solve()', same)
    H = sp.hessian(f, (x, y))
    idx_sum = 0
    pts = []
    for xv, yv in real:
        det = sp.nsimplify(sp.simplify(H.subs({x: xv, y: yv}).det()))
        tr_ = sp.simplify(H.subs({x: xv, y: yv}).trace())
        ind = 1 if det > 0 else (-1 if det < 0 else None)
        kind = 'maximum' if det > 0 and tr_ < 0 else ('minimum' if det > 0 else ('saddle' if det < 0 else 'degenerate'))
        val = sp.nsimplify(sp.simplify(f.subs({x: xv, y: yv})))
        print(f'    z = {xv} + i*({yv}):  924*r6 = {sp.nsimplify(924 * val)}  in-line {kind}')
        gate(f'    nondegenerate in the line (index {ind})', ind is not None)
        idx_sum += ind or 0
        pts.append(ind)
        add(name, f'z = {xv} + i({yv}), in-line {kind}', {m: b1.get(m, 0) + (xv + sp.I * yv) * b2.get(m, 0) for m in ms}, val)
    # the point at infinity: chart w = 1/z, u ~ w b1 + b2
    w = x + sp.I * y
    uinf = {m: w * b1.get(m, 0) + b2.get(m, 0) for m in ms if b1.get(m, 0) != 0 or b2.get(m, 0) != 0}
    finf = sp.simplify(rho6_sq(uinf) / norm_sq(uinf) ** 2)
    ginf = [sp.simplify(sp.diff(finf, v_).subs({x: 0, y: 0})) for v_ in (x, y)]
    Hinf = sp.hessian(finf, (x, y)).subs({x: 0, y: 0})
    dinf, tinf = sp.simplify(Hinf.det()), sp.simplify(Hinf.trace())
    kinf = 'maximum' if dinf > 0 and tinf < 0 else ('minimum' if dinf > 0 else ('saddle' if dinf < 0 else 'degenerate'))
    gate(f'    z = infinity (v0) is critical in the line (gradient {ginf}), in-line {kinf}', all(g == 0 for g in ginf) and dinf != 0)
    idx_sum += 1 if dinf > 0 else -1
    add(name, f'z = infinity (v0), in-line {kinf}', {0: 1}, finf.subs({x: 0, y: 0}))
    gate(f'  {name}: Poincare-Hopf indices sum to {idx_sum} = chi(S^2) = 2 over {len(real) + 1} critical points', idx_sum == 2)
    gate(f'  arm: dropping one point breaks the sum ({idx_sum - pts[0]} != 2)', idx_sum - pts[0] != 2)

print('(F) every point found is critical on the whole sphere (numerical control; Palais guarantees it)')
worst = max(grad_norm(un) for _, _, un, _ in found)
gate(f'  max tangential gradient over {len(found)} points = {worst:.1e}', worst < 1e-7)
nc = np.zeros(7, complex); nc[idx[2]] = np.sqrt(0.5); nc[idx[-1]] = np.sqrt(0.5)
gate(f'  arm: a non-critical point on the C3 line (s = 1/2) fails it ({grad_norm(nc):.1e})', grad_norm(nc) > 1e-3)
gate('  every exact value agrees with the numerical value at its point', max(abs(float(v_) - r6n(un)) for _, _, un, v_ in found) < 1e-12)

print('(U) the union modulo rotation')
orbits = []          # list of [representative vector, value, [(line, label)]]
for line, label, un, val in found:
    placed = False
    for orb in orbits:
        if orb[1] == val and same_orbit(un, orb[0]):
            orb[2].append((line, label))
            placed = True
            break
    if not placed:
        orbits.append([un, val, [(line, label)]])
# ties in value between different orbits must be separated by a rotation invariant
for (i, o1), (j, o2) in itertools.combinations(enumerate(orbits), 2):
    if o1[1] == o2[1]:
        sep_inv = np.linalg.norm(invariants(o1[0]) - invariants(o2[0])) > 1e-6
        sep_iso = rot_orders(o1[0]) != rot_orders(o2[0])
        gate(f'  tie at 924*r6 = {sp.nsimplify(924 * o1[1])}: separated by the invariant triple ({sep_inv}) or the isotropy orders ({sep_iso})', sep_inv or sep_iso)
orbits.sort(key=lambda o: float(o[1]))
print(f'  {len(orbits)} distinct critical orbits on the seven lines and the six points:')
for un, val, where in orbits:
    iso = rot_orders(un)
    print(f'    924*r6 = {str(sp.nsimplify(924 * val)):>10s}  ({float(924 * val):9.4f})  rotation orders {iso}  on: ' + '; '.join(f'{l} [{lab}]' for l, lab in where))
gate('  arm: the union logic merges two rotated copies of the hexagon', same_orbit(np.eye(7)[0] + np.eye(7)[6], Dw([0.3, -1.1, 0.4]) @ (np.eye(7)[0] + np.eye(7)[6])))
print(f'{len(fails)} FAILED: {fails}' if fails else 'ALL PASS')
