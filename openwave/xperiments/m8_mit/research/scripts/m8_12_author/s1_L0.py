"""S1 derivation, step 1 of 3: L0, the fixed-locus classification (author-side; no Hessian, no r6 critical data).

Claim to be shown: up to rotation, the character-fixed loci in P(V_3) of projective dimension <= 1 are exactly six
points (v3, v2, v1, v0, the octahedron, the hexagon) and seven lines (C3 on {v2,v-1}; C4 on {v3,v-1}; C4 on {v2,v-2};
C5 on {v3,v-2}; C6 on {v3,v-3}; D3 on {v3+v-3, v0}; D2 on {v2+v-2, v0}); every other character space of a finite
rotation group has projective dimension >= 2.

Routes, which must agree group by group:
  R1  characters: the 1-dim characters of each group, found by closure-consistency on its generators, and their
      multiplicities from chi_V3(theta) = 1 + 2 sum_{m=1..3} cos(m theta);
  R2  explicit joint eigenspaces of the generators' spin-3 matrices.
Classification: every point or line from R2 is carried to a canonical representative by an explicitly found rotation
(residual < 1e-12); the canonical classes are shown distinct by rotation invariants (ranges of the multipole norms).
Arms: D2 removed must lose its line class; T alone must give no line; a spin-7/2 character in R1 must break agreement;
the two C4 lines must not map to each other.
Groups: C_n and D_n for n up to 12 (n >= 7 sampled; the general argument is stated in the note), T, O, I; the
continuous subgroups give weight states only, by the same argument.
"""
import itertools
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize
from scipy.spatial.transform import Rotation

j = 3
ms = list(range(j, -j - 1, -1))
idx = {m: i for i, m in enumerate(ms)}
Jz = np.diag([float(m) for m in ms]).astype(complex)
Jp = np.zeros((7, 7), complex)
for m in ms[1:]:
    Jp[idx[m + 1], idx[m]] = np.sqrt(j * (j + 1) - m * (m + 1))
Jx, Jy = (Jp + Jp.conj().T) / 2, (Jp - Jp.conj().T) / 2j
J = (Jx, Jy, Jz)
fails = []


def gate(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok:
        fails.append(name)


def rotmat(n, th):
    n = np.asarray(n, float) / np.linalg.norm(n)
    K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return expm(th * K)


def axis_angle(R):
    """Axis and angle from scipy's rotation vector, which stays accurate at angle pi (an arccos of the trace does
    not: near -1 it loses half the digits, and the first run of this script failed on exactly that)."""
    w = Rotation.from_matrix(R).as_rotvec()
    th = float(np.linalg.norm(w))
    if th < 1e-12:
        return np.array([0, 0, 1.0]), 0.0
    return w / th, th


def D3(R):
    n, th = axis_angle(R)
    return expm(-1j * th * (n[0] * Jx + n[1] * Jy + n[2] * Jz))


def key(R):
    return (np.round(R, 8) + 0.0).tobytes()


def closure(gens, cap=200):
    els = {key(np.eye(3)): np.eye(3)}
    frontier = [np.eye(3)]
    while frontier:
        new = []
        for a in frontier:
            for g in gens:
                b = g @ a
                kb = key(b)
                if kb not in els:
                    els[kb] = b
                    new.append(b)
                    if len(els) > cap:
                        raise RuntimeError('closure exceeded cap')
        frontier = new
    return list(els.values())


def order(g):
    M = np.eye(3)
    for k in range(1, 61):
        M = g @ M
        if np.allclose(M, np.eye(3), atol=1e-9):
            return k
    raise RuntimeError('order > 60')


phi = (1 + 5 ** 0.5) / 2
GROUPS = {}
for n in range(1, 13):
    GROUPS[f'C{n}'] = [rotmat([0, 0, 1], 2 * np.pi / n)]
for n in range(2, 13):
    GROUPS[f'D{n}'] = [rotmat([0, 0, 1], 2 * np.pi / n), rotmat([1, 0, 0], np.pi)]
GROUPS['T'] = [rotmat([0, 0, 1], np.pi), rotmat([1, 0, 0], np.pi), rotmat([1, 1, 1], 2 * np.pi / 3)]
GROUPS['O'] = [rotmat([0, 0, 1], np.pi / 2), rotmat([1, 1, 1], 2 * np.pi / 3)]
GROUPS['I'] = [rotmat([0, 1, phi], 2 * np.pi / 5), rotmat([1, 1, 1], 2 * np.pi / 3)]
ORDER = {**{f'C{n}': n for n in range(1, 13)}, **{f'D{n}': 2 * n for n in range(2, 13)}, 'T': 12, 'O': 24, 'I': 60}


def chars_by_closure(gens, els):
    """R1: all 1-dim characters, as tuples of values on the generators, by closure-consistency. A character's value
    on a generator of order k is a k-th root of unity, so only those candidates are tried."""
    N = len(els)
    cands = [[np.exp(2j * np.pi * l / order(g)) for l in range(order(g))] for g in gens]
    found = []
    for vals in itertools.product(*cands):
        table = {key(np.eye(3)): (np.eye(3), 1.0 + 0j)}
        frontier = [(np.eye(3), 1.0 + 0j)]
        ok = True
        while frontier and ok:
            nxt = []
            for a, va in frontier:
                for g, vg in zip(gens, vals):
                    b, vb = g @ a, vg * va
                    kb = key(b)
                    if kb in table:
                        if abs(table[kb][1] - vb) > 1e-9:
                            ok = False
                            break
                    else:
                        table[kb] = (b, vb)
                        nxt.append((b, vb))
                if not ok:
                    break
            frontier = nxt
        if ok and len(table) == N:
            found.append((tuple(np.array(vals)), list(table.values())))
    return found


def chi_V(theta, spin=3.0):
    if spin == 3.0:
        return 1 + 2 * sum(np.cos(m * theta) for m in (1, 2, 3))
    ms_ = np.arange(-spin, spin + 1)
    return float(np.real(np.sum(np.exp(-1j * ms_ * theta))))


def route1(name, spin=3.0):
    gens = GROUPS[name]
    els = closure(gens)
    out = []
    for vals, table in chars_by_closure(gens, els):
        mult = sum(np.conj(v) * chi_V(axis_angle(e)[1], spin) for e, v in table) / len(table)
        out.append((vals, mult))
    return els, out


def route2(name):
    gens = GROUPS[name]
    Ds = [D3(g) for g in gens]
    spaces = []
    eigs = [np.unique(np.round(np.linalg.eigvals(Dg), 9)) for Dg in Ds]
    for lam in itertools.product(*eigs):
        A = np.vstack([Dg - l * np.eye(7) for Dg, l in zip(Ds, lam)])
        U, s, Vh = np.linalg.svd(A)
        null = Vh.conj().T[:, s < 1e-8] if len(s) == 7 else None
        k = int(np.sum(s < 1e-8))
        if k > 0:
            spaces.append((tuple(np.round(np.array(lam), 9)), k, Vh.conj().T[:, 7 - k:]))
    return spaces


gate('control: D3 of a half-turn about z has eigenvalues +-1 only (the angle-pi case)', np.allclose(sorted(np.round(np.linalg.eigvals(D3(rotmat([0, 0, 1], np.pi))).real, 9)), [-1] * 4 + [1] * 3))
print('(R) route agreement, group by group')
ALL = []
for name in GROUPS:
    els, r1 = route1(name)
    r2 = route2(name)
    gate_ok = len(els) == ORDER[name]
    m1 = [(np.array(v), int(round(m.real))) for v, m in r1 if round(m.real) > 0]
    frac = [m for v, m in r1 if abs(m - round(m.real)) > 1e-9]
    m2 = [(np.array(lam), k) for lam, k, _ in r2]
    match = [[i for i, (v1, k1) in enumerate(m1) if np.allclose(v1, lam, atol=1e-8)] for lam, k in m2]
    agree = (len(m1) == len(m2) and all(len(h) == 1 for h in match) and len({h[0] for h in match}) == len(m2)
             and all(m1[h[0]][1] == k for h, (lam, k) in zip(match, m2)))
    dims = sorted(k for _, k, _ in r2)
    gate(f'  {name:4s} order {len(els):2d}: {len(r1)} one-dim characters; multiplicities integral; routes agree; space dims {dims}',
         gate_ok and not frac and agree and sum(dims) <= 7)
    ALL += [(name, lam, k, B) for lam, k, B in r2]
_, r1bad = route1('D2', spin=2.0)
gate('arm: a spin-2 character table (a genuine SO(3) representation) in R1 does not reproduce D2 in V3', sorted(int(round(m.real)) for _, m in r1bad if round(m.real) > 0) != sorted(k for _, k, _ in route2('D2')))

v = lambda **c: sum(c[k] * np.eye(7)[idx[int(k.replace('m', '-').replace('p', ''))]] for k in c)
s2 = 2 ** -0.5
CAN_P = {'v3': v(p3=1), 'v2': v(p2=1), 'v1': v(p1=1), 'v0': v(p0=1), 'octahedron': v(p2=s2, m2=s2), 'hexagon': v(p3=s2, m3=s2)}
CAN_L = {'C3 {v2,v-1}': np.stack([v(p2=1), v(m1=1)], 1), 'C4 {v3,v-1}': np.stack([v(p3=1), v(m1=1)], 1),
         'C4 {v2,v-2}': np.stack([v(p2=1), v(m2=1)], 1), 'C5 {v3,v-2}': np.stack([v(p3=1), v(m2=1)], 1),
         'C6 {v3,v-3}': np.stack([v(p3=1), v(m3=1)], 1), 'D3 {v3+v-3,v0}': np.stack([v(p3=s2, m3=s2), v(p0=1)], 1),
         'D2 {v2+v-2,v0}': np.stack([v(p2=s2, m2=s2), v(p0=1)], 1)}
rng = np.random.default_rng(20260919)
Dw = lambda w: expm(-1j * (w[0] * Jx + w[1] * Jy + w[2] * Jz))


def carry_w(B, C, starts=48):
    k = B.shape[1]
    f = lambda w: k - np.linalg.norm(C.conj().T @ Dw(w) @ B) ** 2
    best, bw = np.inf, None
    for _ in range(starts):
        r = minimize(f, rng.normal(size=3) * 2, method='BFGS', options={'gtol': 1e-13})
        if r.fun < best:
            best, bw = r.fun, r.x
        if best < 1e-13:
            break
    return best, bw


def carry(B, C, starts=48):
    return carry_w(B, C, starts)[0]


print('(C) every point and line is carried to a canonical class; the classes are pairwise distinct by (N), so the class is unique')
classes = {'point': {}, 'line': {}}
for name, lam, k, B in ALL:
    if k not in (1, 2):
        continue
    kind, CAN = ('point', CAN_P) if k == 1 else ('line', CAN_L)
    hit, rec = [], ''
    for c, C in CAN.items():
        res, w = carry_w(B, C.reshape(7, -1))
        if res < 1e-12:
            recheck = B.shape[1] - np.linalg.norm(C.reshape(7, -1).conj().T @ Dw(w) @ B) ** 2
            hit, rec = [c], f'  rotation vector {np.round(w, 9).tolist()}, residual {recheck:.1e}'
            break
    gate(f'  {name:4s} {kind} -> {hit}{rec}', len(hit) == 1)
    if len(hit) == 1:
        classes[kind].setdefault(hit[0], []).append(name)
gate(f'the point classes found are exactly the six canonical points ({sorted(classes["point"])})', sorted(classes['point']) == sorted(CAN_P))
gate(f'the line classes found are exactly the seven canonical lines ({len(classes["line"])})', sorted(classes['line']) == sorted(CAN_L))
big = sorted({(name, k) for name, lam, k, B in ALL if k >= 3})
print('  character spaces of dimension >= 3 (projective dimension >= 2):', big)
gate('every other character space has projective dimension >= 2 (dimensions are only 1, 2, or >= 3)', all(k in (1, 2) or k >= 3 for _, _, k, _ in ALL))

print("(X) the exact distinctness argument, checked: each line's fixing group within O(2)_z, and the C4 pair")
print("    (a fixing rotation must fix a weight ray [v_a], a != 0, on a cyclic line, so it lies in SO(2)_z; on a dihedral line")
print("     it must fix the zonal ray [v0], so it lies in O(2)_z; the steps below are 6 degrees and 3 degrees)")
def scalar_on(B, R3):
    # a scalar of modulus 1: a zero matrix is also 'scalar' when D3(R3) maps the line to an orthogonal one, and the
    # first run of this block counted exactly those (60 spurious half-turn axes on three cyclic lines)
    M = B.conj().T @ D3(R3) @ B
    c = M[0, 0]
    return abs(abs(c) - 1) < 1e-10 and np.allclose(M, c * np.eye(M.shape[0]), atol=1e-10), c
gate("  arm: a half-turn that maps span(v2,v-1) to an orthogonal line is not counted as scalar",
     not scalar_on(np.stack([v(p2=1), v(m1=1)], 1), rotmat([1, 0, 0], np.pi))[0])
FIX = {'C3 {v2,v-1}': (3, 0), 'C4 {v3,v-1}': (4, 0), 'C4 {v2,v-2}': (4, 0), 'C5 {v3,v-2}': (5, 0), 'C6 {v3,v-3}': (6, 0),
       'D3 {v3+v-3,v0}': (3, 3), 'D2 {v2+v-2,v0}': (2, 2)}
for c, (n_rot, n_flip) in FIX.items():
    B = CAN_L[c]
    rots = sum(scalar_on(B, rotmat([0, 0, 1], 2 * np.pi * k / 60))[0] for k in range(1, 61))
    flips = sum(scalar_on(B, rotmat([np.cos(np.pi * k / 60), np.sin(np.pi * k / 60), 0], np.pi))[0] for k in range(60))
    gate(f"  {c}: {rots} rotations about z act as scalars (expected {n_rot}); {flips} horizontal half-turn axes do (expected {n_flip})",
         rots == n_rot and flips == n_flip)
B_bad = np.stack([v(p3=s2, m3=s2), v(p1=1)], 1)                      # a planted wrong line: v1 in place of v0
bad = (sum(scalar_on(B_bad, rotmat([0, 0, 1], 2 * np.pi * k / 60))[0] for k in range(1, 61)),
       sum(scalar_on(B_bad, rotmat([np.cos(np.pi * k / 60), np.sin(np.pi * k / 60), 0], np.pi))[0] for k in range(60)))
gate(f"  arm: the planted line span(v3+v-3, v1) does not pass as D3 (counts {bad}, not (3, 3))", bad != (3, 3))
ch_a = scalar_on(CAN_L['C4 {v3,v-1}'], rotmat([0, 0, 1], np.pi / 2))[1]
ch_b = scalar_on(CAN_L['C4 {v2,v-2}'], rotmat([0, 0, 1], np.pi / 2))[1]
gate(f"  the C4 generator acts by {np.round(ch_a, 9)} on span(v3,v-1) and by {np.round(ch_b, 9)} on span(v2,v-2): the classes (i, -i) and (-1) differ",
     abs(abs(ch_a.imag) - 1) < 1e-9 and abs(ch_a.real) < 1e-9 and abs(ch_b + 1) < 1e-9)
def f2(u):
    u = u / np.linalg.norm(u)
    return sum(abs(np.vdot(u, A @ u)) ** 2 for A in (Jx, Jy, Jz))
ok_a = all(abs(f2(np.sqrt(s_) * v(p3=1) + np.sqrt(1 - s_) * np.exp(0.7j) * v(m1=1)) - (4 * s_ - 1) ** 2) < 1e-12 for s_ in np.linspace(0, 1, 11))
ok_b = all(abs(f2(np.sqrt(s_) * v(p2=1) + np.sqrt(1 - s_) * np.exp(0.7j) * v(m2=1)) - 4 * (2 * s_ - 1) ** 2) < 1e-12 for s_ in np.linspace(0, 1, 11))
gate("  |f|^2 = (4s-1)^2 on span(v3,v-1), range [0,9]; |f|^2 = 4(2s-1)^2 on span(v2,v-2), range [0,4]", ok_a and ok_b)
print('(N) corroboration: ranges of rotation invariants on a grid')
from sympy.physics.wigner import clebsch_gordan as _CG
CGT = {(K, m1, m2): float(_CG(3, 3, K, m1, m2, m1 + m2)) for K in range(1, 7) for m1 in ms for m2 in ms if abs(m1 + m2) <= K}


def multipole_norms(u):
    u = u / np.linalg.norm(u)
    tu = np.array([(-1) ** (j - m) * np.conj(u[idx[-m]]) for m in ms])
    out = []
    for K in range(1, 7):
        tot = 0.0
        for Q in range(-K, K + 1):
            s_ = sum(CGT[K, m1, Q - m1] * u[idx[m1]] * tu[idx[Q - m1]] for m1 in ms if (K, m1, Q - m1) in CGT)
            tot += abs(s_) ** 2
        out.append(tot)
    return np.array(out)


pts = {c: multipole_norms(u) for c, u in CAN_P.items()}
gate('the six points have pairwise distinct multipole-norm vectors', min(np.linalg.norm(pts[a] - pts[b]) for a, b in itertools.combinations(pts, 2)) > 1e-3)
grid = [(np.cos(t), np.sin(t) * np.exp(1j * p)) for t in np.linspace(0, np.pi / 2, 25) for p in np.linspace(0, 2 * np.pi, 13)]
lranges = {}
for c, B in CAN_L.items():
    vals = np.array([multipole_norms(B @ np.array(z)) for z in grid])
    lranges[c] = np.concatenate([vals.min(0), vals.max(0)])
dmin = min(np.linalg.norm(lranges[a] - lranges[b]) for a, b in itertools.combinations(lranges, 2))
print('  invariant ranges on each line, ||rho_K||^2 for K = 1..6 (min | max):')
for c in CAN_L:
    print(f'    {c:16s} min {np.round(lranges[c][:6], 4)}  max {np.round(lranges[c][6:], 4)}')
pair = min(itertools.combinations(lranges, 2), key=lambda ab: np.linalg.norm(lranges[ab[0]] - lranges[ab[1]]))
gate(f'the closest pair {pair} is separated without K = 6 as well (K <= 5 only: {np.linalg.norm(lranges[pair[0]][[0,1,2,3,4,6,7,8,9,10]] - lranges[pair[1]][[0,1,2,3,4,6,7,8,9,10]]):.4f})',
     np.linalg.norm(lranges[pair[0]][[0, 1, 2, 3, 4, 6, 7, 8, 9, 10]] - lranges[pair[1]][[0, 1, 2, 3, 4, 6, 7, 8, 9, 10]]) > 1e-2)
gate(f'the seven lines have pairwise distinct invariant-range vectors (min separation {dmin:.3f})', dmin > 1e-2)
gate('arm: the two C4 lines are not carried to each other', carry(CAN_L['C4 {v3,v-1}'], CAN_L['C4 {v2,v-2}']) > 1e-3)

print('(A) arms on the enumeration')
def line_classes(names):
    out = set()
    for name, lam, k, B in ALL:
        if name in names and k == 2:
            for c, C in CAN_L.items():
                if carry(B, C, starts=16) < 1e-12:
                    out.add(c)
    return out
without_D2 = line_classes([nm for nm in GROUPS if nm != 'D2'])
gate(f'arm: with D2 removed its line class is lost ({len(without_D2)} classes)', 'D2 {v2+v-2,v0}' not in without_D2 and len(without_D2) == 6)
gate('arm: T alone gives no line', line_classes(['T']) == set())
print(f'{len(fails)} FAILED: {fails}' if fails else 'ALL PASS')
