"""Item 5: exact reasons for the vanishing levels.
For J != 3, Pi_J N(Phi_{sigma,u}) has coefficients M^{6,J}(u)_M (R^{6,J}_sigma)_{Ka} (see RETURN.md), and F_J = ||R^{6,J}_sigma||^2 != 0
at every sector level (exact_parts.py). So a sector level vanishes iff M^{6,J}(v) = 0 for the ray.
(a) exact symbolic evaluation of every component of M^{6,J}(v) at the zero levels (sympy, exact radicals);
(b) R2: CG parity rule <6 0; 3 0|J 0> = 0 for J even, and rho_6(v_0) has only an M = 0 component;
(c) R3: the ray is an eigenline of a 48-element subgroup H of SU(2) (a 45-degree z-rotated binary octahedral group, verified here by
    closure and by the eigenline test on all elements); M^{6,J}(R3) must lie in the chi-isotypic part of V_J, and its multiplicity,
    computed by characters, is 0 at J = 4, 5, 8;
(d) R4 at J = 8: explicit two-term CG sum, exact."""
import math
import numpy as np
import sympy as sp
from common import cg_exact, quat_to_su2, Dj_generic


def cgs(j1, m1, j2, m2, J, M):
    S, A = cg_exact(j1, m1, j2, m2, J, M)
    return sp.Integer(0) if S == 0 else sp.Rational(S) * sp.sqrt(sp.Rational(A))


def theta(u):
    return [(-1) ** ((i - 3) % 2) * sp.conjugate(u[6 - i]) for i in range(7)]


def couple(x, j1, y, j2, J):
    res = []
    for M in range(-J, J + 1):
        t = sp.Integer(0)
        for m1 in range(-j1, j1 + 1):
            m2 = M - m1
            if abs(m2) <= j2:
                t += cgs(j1, m1, j2, m2, J, M) * x[m1 + j1] * y[m2 + j2]
        res.append(sp.radsimp(sp.expand(t)))
    return res


e = lambda m: [sp.Integer(1) if i - 3 == m else sp.Integer(0) for i in range(7)]
R2 = e(0)
R3 = [x + y for x, y in zip(e(2), e(-2))]
R4 = [x + y for x, y in zip(e(3), e(-3))]
print("(a) exact components of M^{6,J}(v) at the vanishing levels:")
for name, v, Js in (("R2", R2, (4, 6, 8)), ("R3", R3, (4, 5, 8)), ("R4", R4, (8,))):
    rho6 = couple(v, 3, theta(v), 3, 6)
    for J in Js:
        comps = couple(rho6, 6, v, 3, J)
        print(f"   {name} J={J} (n={2 * J}): all {2 * J + 1} components exactly zero: {all(c == 0 for c in comps)}")
        assert all(c == 0 for c in comps)
# a check that the exact evaluation CAN return nonzero: the neighbouring levels
for name, v, J in (("R2", R2, 5), ("R3", R3, 7), ("R4", R4, 9)):
    comps = couple(couple(v, 3, theta(v), 3, 6), 6, v, 3, J)
    print(f"   control {name} J={J}: nonzero components exist: {any(c != 0 for c in comps)}")

print("(b) R2: rho_6(v_0) components:", [c for c in couple(R2, 3, theta(R2), 3, 6)])
print("    <6 0; 3 0 | J 0> for J = 3..9:", [str(cgs(6, 0, 3, 0, J, 0)) for J in range(3, 10)])

print("(d) R4, J = 8: rho_6(R4) =", couple(R4, 3, theta(R4), 3, 6))
rho = couple(R4, 3, theta(R4), 3, 6)
c1, c2 = cgs(6, 6, 3, -3, 8, 3), cgs(6, 0, 3, 3, 8, 3)
print(f"    M^(6,8)(R4)_3 = <6 6;3 -3|8 3> rho_6,6 + <6 0;3 3|8 3> rho_6,0 = ({c1})({rho[12]}) + ({c2})({rho[6]}) = "
      f"{sp.radsimp(sp.expand(c1 * rho[12] + c2 * rho[6]))}")
c1m, c2m = cgs(6, -6, 3, 3, 8, -3), cgs(6, 0, 3, -3, 8, -3)
print(f"    M^(6,8)(R4)_-3 = ({c1m})({rho[0]}) + ({c2m})({rho[6]}) = {sp.radsimp(sp.expand(c1m * rho[0] + c2m * rho[6]))}")

# (c) stabilizer argument for R3
s2 = 1 / math.sqrt(2)
units = []
basis = [np.array(v, float) for v in np.eye(4)]
for b in basis:
    units += [b, -b]
for sg in [(a, b, c, d) for a in (1, -1) for b in (1, -1) for c in (1, -1) for d in (1, -1)]:
    units.append(0.5 * np.array(sg, float))
for i in range(4):
    for k in range(i + 1, 4):
        for si in (1, -1):
            for sk in (1, -1):
                q = np.zeros(4)
                q[i], q[k] = si * s2, sk * s2
                units.append(q)


def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return np.array([a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2, a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                     a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2, a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2])


c = np.array([math.cos(math.pi / 8), math.sin(math.pi / 8), 0, 0])     # rotation by pi/4 about the quantization axis (quaternion i)
cinv = c * np.array([1, -1, -1, -1])
H = [qmul(qmul(c, u), cinv) for u in units]
# closure check (can fail)
key = lambda q: tuple(np.round(q, 9))
Hs = {key(q) for q in H}
closed = all(key(qmul(a, b)) in Hs for a in H for b in H)
print(f"(c) |H| = {len(Hs)}, closed under multiplication: {closed}")
assert closed and len(Hs) == 48
fact = lambda k: float(math.factorial(k))
v3 = np.array([complex(x) for x in R3]) / math.sqrt(2)
chis = []
worst = 0
for q in H:
    U = quat_to_su2(list(q), 1.0, 1j)
    D = np.array(Dj_generic(6, U, fact, math.sqrt, 0j))
    w = D @ v3
    lam = np.vdot(v3, w)
    worst = max(worst, np.linalg.norm(w - lam * v3))
    chis.append(lam)
print(f"    R3 is an eigenline of every element of H: max |D v - chi v| = {worst:.2e}")
assert worst < 1e-10
mult = {}
for J in range(3, 10):
    s = 0
    for q, ch in zip(H, chis):
        U = quat_to_su2(list(q), 1.0, 1j)
        D = np.array(Dj_generic(2 * J, U, fact, math.sqrt, 0j))
        s += np.conj(ch) * np.trace(D)
    mult[J] = s / len(H)
print("    multiplicity of chi in V_J, J=3..9:", {J: f"{mult[J].real:.10f}{mult[J].imag:+.1e}i" for J in mult})
for J in mult:
    assert abs(mult[J] - round(mult[J].real)) < 1e-9
print("    rounded:", {J: int(round(mult[J].real)) for J in mult})


# (e) stabilizers of R4 and R5 (used for items 3 and 6): generate numerically, test the eigenline, count chi-multiplicities
def stab_mult(gens, v, Js, label):
    key = lambda q: tuple(np.round(q, 9))
    elems = {key(np.array([1.0, 0, 0, 0])): np.array([1.0, 0, 0, 0])}
    frontier = list(elems.values())
    while frontier:
        new = []
        for x in frontier:
            for g in gens:
                y = qmul(x, g)
                if key(y) not in elems:
                    elems[key(y)] = y
                    new.append(y)
        frontier = new
    Hl = list(elems.values())
    vv = np.array(v, complex)
    vv = vv / np.linalg.norm(vv)
    ch, worst = [], 0
    for q in Hl:
        D = np.array(Dj_generic(6, quat_to_su2(list(q), 1.0, 1j), fact, math.sqrt, 0j))
        w = D @ vv
        lam = np.vdot(vv, w)
        worst = max(worst, np.linalg.norm(w - lam * vv))
        ch.append(lam)
    ms = {}
    for J in Js:
        s = sum(np.conj(c_) * np.trace(np.array(Dj_generic(2 * J, quat_to_su2(list(q), 1.0, 1j), fact, math.sqrt, 0j)))
                for q, c_ in zip(Hl, ch)) / len(Hl)
        assert abs(s - round(s.real)) < 1e-9
        ms[J] = int(round(s.real))
    print(f"(e) {label}: |H| = {len(Hl)}, eigenline max dev {worst:.1e}, chi-multiplicity in V_J: {ms}")
    assert worst < 1e-10
    return ms


zq = lambda a: np.array([math.cos(a), math.sin(a), 0, 0])          # acts on v_m by e^{2 i m a}
jq = np.array([0, 0, 1.0, 0])                                        # a half-turn about a horizontal axis
stab_mult([zq(math.pi / 6), jq], [float(x) for x in R4], (3, 8), "R4, <e^{i pi/6}, j> (binary dihedral)")
# R5 in index order i = m + 3: sin t at m = -3 (i = 0), cos t at m = 2 (i = 5)
stab_mult([zq(math.pi / 5)], [2 * math.sqrt(3) / 5, 0, 0, 0, 0, math.sqrt(13) / 5, 0], (3,), "R5, <e^{i pi/5}> (cyclic)")
