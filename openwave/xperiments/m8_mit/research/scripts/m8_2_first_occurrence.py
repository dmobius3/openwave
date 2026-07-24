"""M8.2 V1c: McKay first-occurrence certification tables for M4 and M7.

Computes, from first principles (the 2I character table built by the McKay
recursion + SU(2) characters), the level at which each of the 8 nontrivial 2I
irreps first occurs in each family's twisted harmonic tower, for the three flat
connections sigma in {trivial, standard Q, Galois Q'}.

The certification target M8.5 must reproduce. Every number is script-derived;
the M4 table is cross-checked against mass-spectrum.md section 4 (the author's
own per-connection j_first table), which validates the whole computation.

Run: python3 m8_2_first_occurrence.py
"""
import math
import numpy as np

# ---- 2I irreps: R0..R8, dims (E8 marks), McKay distances -------------------
labels = ['R0', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8']
dims = [1, 2, 2, 3, 3, 4, 4, 5, 6]
dist = [0, 1, 7, 2, 6, 6, 3, 4, 5]
idx = {l: i for i, l in enumerate(labels)}

# ---- McKay graph = affine E8: tensoring by R1 (defining V_1) is adjacency ---
edges = [('R0', 'R1'), ('R1', 'R3'), ('R3', 'R6'), ('R6', 'R7'),
         ('R7', 'R8'), ('R8', 'R5'), ('R5', 'R2'), ('R8', 'R4')]
A = np.zeros((9, 9))
for a, b in edges:
    A[idx[a], idx[b]] = A[idx[b], idx[a]] = 1

# ---- 9 conjugacy classes of 2I: (SU(2) eigen-angle theta, class size) -------
classes = [(0.0, 1), (math.pi, 1), (math.pi / 2, 30), (math.pi / 3, 20),
           (2 * math.pi / 3, 20), (math.pi / 5, 12), (2 * math.pi / 5, 12),
           (3 * math.pi / 5, 12), (4 * math.pi / 5, 12)]
thetas = np.array([c[0] for c in classes])
sizes = np.array([c[1] for c in classes], float)


def chiV(n, th):
    """SU(2) character of V_n (dim n+1, spin n/2) at eigen-angle th."""
    s = math.sin(th)
    if abs(s) < 1e-12:                      # theta = 0 or pi
        return (n + 1) * math.cos(n * th)   # (n+1) or (n+1)(-1)^n
    return math.sin((n + 1) * th) / s


# ---- multiplicities m_n of V_n|_2I via V_{n+1} = V_1 (x) V_n - V_{n-1} -------
N = 24
m = [np.zeros(9) for _ in range(N)]
m[0][idx['R0']] = 1.0
m[1][idx['R1']] = 1.0
for n in range(1, N - 1):
    m[n + 1] = A.dot(m[n]) - m[n - 1]

# ---- solve for the irreducible character table chi[irrep][class] ------------
Mmat = np.array([m[n] for n in range(9)])           # rows n=0..8, cols irrep
chi = np.zeros((9, 9))
for c in range(9):
    rhs = np.array([chiV(n, thetas[c]) for n in range(9)])
    chi[:, c] = np.linalg.solve(Mmat, rhs)

# ---- VERIFY 1: orthonormality of the character table -----------------------
G = np.array([[np.sum(sizes * chi[i] * chi[j]) / 120.0 for j in range(9)]
              for i in range(9)])
ortho_err = np.max(np.abs(G - np.eye(9)))

# ---- VERIFY 2: Sym^2 Q = R3 (= V_2|_2I), and Sym^2 Q' = R4 (Galois partner) -
sym2Q_is_R3 = np.allclose(m[2], np.eye(9)[idx['R3']], atol=1e-6)
# R3 and R4 are the Galois pair Sym^2 Q, Sym^2 Q': equal on the rational
# (non-golden) classes, distinct (sqrt5-conjugate) on the four golden classes.
golden = [5, 6, 7, 8]                       # theta in {1,2,3,4}*pi/5
non_golden = [0, 1, 2, 3, 4]
r4_is_galois_of_r3 = (
    all(abs(chi[idx['R3']][c] - chi[idx['R4']][c]) < 1e-6 for c in non_golden)
    and all(abs(chi[idx['R3']][c] - chi[idx['R4']][c]) > 0.5 for c in golden))

tauchar = {'trivial': chi[idx['R0']], 'standard': chi[idx['R3']],
           'galois': chi[idx['R4']]}        # coefficient = Sym^2(sigma)


def mult(rho, tau_c, n):
    """multiplicity of irrep rho in V_n (x) tau."""
    v = np.array([chiV(n, thetas[c]) for c in range(9)])
    return np.sum(sizes * v * tau_c * chi[idx[rho]]) / 120.0


def first_level(rho, tau_c):
    for n in range(N - 1):
        if mult(rho, tau_c, n) > 0.5:
            return n
    return None


# ---- M4 (0-form) first-occurrence table: level n, j_first = n/2 -------------
order = ['R1', 'R3', 'R6', 'R7', 'R8', 'R4', 'R5', 'R2', 'R0']
m4 = {}
for rho in order:
    m4[rho] = {conn: first_level(rho, tauchar[conn]) for conn in tauchar}

# author's mass-spectrum.md section 4 table (j_first), for cross-check ---------
ms4_jfirst = {
    'R0': (0, 1, 3), 'R1': (0.5, 0.5, 2.5), 'R2': (3.5, 2.5, 1.5),
    'R3': (1, 0, 2), 'R4': (3, 2, 0), 'R5': (3, 2, 1),
    'R6': (1.5, 0.5, 1.5), 'R7': (2, 1, 1), 'R8': (2.5, 1.5, 0.5)}
m4_matches = all(
    (m4[r]['trivial'] / 2, m4[r]['standard'] / 2, m4[r]['galois'] / 2)
    == ms4_jfirst[r] for r in labels)

# ---- M7 (coexact 1-form) first-occurrence, all three connections ------------
# rho occurs in E_m (x) tau_sigma  <=>  E_m contains a constituent of
# rho (x) tau_sigma* = rho (x) tau_sigma (all 2I irreps self-dual, real chars).
# A constituent at McKay distance d enters the coexact tower E_m at level
# d (d>=2), 2 (d=0), 3 (d=1, bipartite since -I in 2I).  eigenvalue m^2/R^2.
def coexact_level(d):
    return 2 if d == 0 else (3 if d == 1 else d)


def m7_first_level(rho, tau_c):
    best = None
    for ci in range(9):                                   # constituents of rho (x) tau
        mlt = np.sum(sizes * chi[idx[rho]] * tau_c * chi[ci]) / 120.0
        if mlt > 0.5:
            lvl = coexact_level(dist[ci])
            best = lvl if best is None else min(best, lvl)
    return best


m7 = {r: {conn: m7_first_level(r, tauchar[conn]) for conn in tauchar}
      for r in order}
# consistency: the trivial column must equal the direct entry rule on rho itself
m7_trivial_ok = all(m7[r]['trivial'] == coexact_level(dist[idx[r]]) for r in labels)

# =========================== REPORT =========================================
print("=" * 70)
print("VERIFICATION")
print("=" * 70)
print(f"  character-table orthonormality error : {ortho_err:.2e}   "
      f"({'PASS' if ortho_err < 1e-6 else 'FAIL'})")
print(f"  Sym^2 Q  == R3 (= V_2|2I)            : {sym2Q_is_R3}")
print(f"  Sym^2 Q' == R4 (Galois partner of R3): {r4_is_galois_of_r3}")
print(f"  M4 table == mass-spectrum.md sec.4   : {m4_matches}   "
      f"({'PASS' if m4_matches else 'FAIL'})")
print(f"  M7 trivial col == direct coexact rule: {m7_trivial_ok}   "
      f"({'PASS' if m7_trivial_ok else 'FAIL'})")

print("\n" + "=" * 70)
print("M4 (0-form)  first-occurrence LEVEL n   [eigenvalue n(n+2)/R^2]")
print("            j_first = n/2 ;  cross-checked vs mass-spectrum sec.4")
print("=" * 70)
print(f"{'irrep':6}{'d':>3}{'trivial':>10}{'standard':>10}{'galois':>10}")
for r in order:
    t, s, g = m4[r]['trivial'], m4[r]['standard'], m4[r]['galois']
    print(f"{r:6}{dist[idx[r]]:>3}{t:>10}{s:>10}{g:>10}")

print("\n" + "=" * 70)
print("M7 (coexact 1-form)  first-occurrence LEVEL m   [eigenvalue m^2/R^2]")
print("  min coexact entry-level over constituents of rho (x) tau_sigma;")
print("  differs from M4 at R0 (d0) and R1 (d1) via the coexact shift")
print("=" * 70)
print(f"{'irrep':6}{'d':>3}{'trivial':>10}{'standard':>10}{'galois':>10}")
for r in order:
    t, s, g = m7[r]['trivial'], m7[r]['standard'], m7[r]['galois']
    print(f"{r:6}{dist[idx[r]]:>3}{t:>10}{s:>10}{g:>10}")
