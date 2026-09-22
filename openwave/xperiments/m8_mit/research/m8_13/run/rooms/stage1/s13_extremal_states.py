"""s13_extremal_states.py -- exact values at the two extremal orbits, the
p_6 = 1 <=> coherent equality analysis, and the stabiliser computation that
rules out rank-two maximisers.
"""

from pathlib import Path
from math import factorial
import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))
IDX = {m: i for i, m in enumerate(ORDER)}


def cgx(j1, m1, j2, m2, Jt, Q):
    if m1 + m2 != Q or abs(m1) > j1 or abs(m2) > j2 or abs(Q) > Jt:
        return sp.Integer(0)
    return sp.nsimplify(sp.simplify(CG(j1, m1, j2, m2, Jt, Q).doit()))


def TK(K):
    out = []
    for Q in range(-K, K + 1):
        M = sp.zeros(7, 7)
        for i, mm in enumerate(ORDER):
            for k, m1 in enumerate(ORDER):
                M[i, k] = cgx(3, m1, 3, -mm, K, Q) * (-1) ** (3 + mm)
        out.append(M)
    return out


def rvec(u, K):
    return [sp.simplify(sp.expand((u.H * M * u)[0, 0])) for M in TK(K)]


def main():
    print("=== exact multipole spectrum of the two extremal orbits ===")
    hexs = sp.zeros(7, 1)
    hexs[IDX[3]] = sp.Rational(1, 1) / sp.sqrt(2)
    hexs[IDX[-3]] = sp.Rational(1, 1) / sp.sqrt(2)
    coh = sp.zeros(7, 1)
    coh[IDX[3]] = 1

    for name, u in (("hexagon (v3+v-3)/sqrt2", hexs), ("coherent v3", coh)):
        print("\n  %s" % name)
        tot = 0
        for K in range(7):
            v = rvec(u, K)
            w = sp.simplify(sum(sp.Abs(x) ** 2 for x in v))
            tot += w
            print("     rhat_%d = %s = %s" % (K, w, sp.N(w, 20)))
        print("     sum = %s" % sp.simplify(tot))
        print("     rho_6 components =", [sp.simplify(x) for x in rvec(u, 6)])

    print("\n  Theta-invariance of the hexagon state:")
    th = sp.Matrix(7, 1, [(-1) ** (3 - m) * sp.conjugate(hexs[IDX[-m]]) for m in ORDER])
    print("     Theta h - h =", sp.simplify(th - hexs).T)
    print("  Theta of the coherent state:")
    thc = sp.Matrix(7, 1, [(-1) ** (3 - m) * sp.conjugate(coh[IDX[-m]]) for m in ORDER])
    print("     Theta v3 =", thc.T, "  (= v_-3, NOT a phase times v3)")
    print("     |<v3, Theta v3>| =", sp.simplify(sp.Abs((coh.H * thc)[0, 0])))

    # ---- p_6 = 1 <=> coherent, equality analysis of the Cauchy-Schwarz proof ----
    print("\n=== p_6 = 1 <=> coherent ===")

    def p6_of_b(bvals):
        """p_6 from the coefficient vector b_k = c_{k-3}, k = 0..6"""
        r = {}
        for n in range(13):
            s = 0
            for k in range(max(0, n - 6), min(6, n) + 1):
                s += bvals[k] * bvals[n - k] / sp.sqrt(sp.factorial(k) * sp.factorial(6 - k)
                                                       * sp.factorial(n - k) * sp.factorial(6 - n + k))
            r[n] = sp.simplify(s)
        num = sum(sp.Abs(r[n]) ** 2 * sp.factorial(n) * sp.factorial(12 - n) for n in range(13))
        den = 924 * (sum(sp.Abs(x) ** 2 for x in bvals)) ** 2
        return sp.simplify(num / den)

    q = sp.symbols("q", positive=True)
    bgeom = [q ** k * sp.sqrt(sp.binomial(6, k)) for k in range(7)]
    print("  b_k = q^k sqrt(C(6,k))  (i.e. P = (q z1 + z2)^6 up to scale):")
    print("     p_6 =", sp.simplify(p6_of_b(bgeom)))
    print("  b = e_6 (state v_3):      p_6 =", p6_of_b([0, 0, 0, 0, 0, 0, 1]))
    print("  b = e_0 (state v_-3):     p_6 =", p6_of_b([1, 0, 0, 0, 0, 0, 0]))
    print("  b = e_3 (state v_0):      p_6 =", p6_of_b([0, 0, 0, 1, 0, 0, 0]))
    print("  b = (1,0,0,0,0,0,1)/.. (hexagon):  p_6 =", p6_of_b([1, 0, 0, 0, 0, 0, 1]))
    print("  b = e_5 :                 p_6 =", p6_of_b([0, 0, 0, 0, 0, 1, 0]))

    # ---- stabiliser of rho_6(hexagon) ----
    print("\n=== stabiliser of rho_6(hexagon) in SO(3) ===")
    print("  rho_6(h) has nonzero components only at Q = -6, 0, +6:")
    print("     Q=+6 : 1/2,  Q=0 : 1/sqrt(924),  Q=-6 : 1/2")
    print("  - rotation by theta about z multiplies the Q component by exp(-i Q theta);")
    print("    invariance needs exp(-6 i theta) = 1, i.e. theta in (pi/3) Z  ->  C_6.")
    print("  - the pi rotation about x sends Q -> -Q; the vector is Q -> -Q symmetric  ->  D_6.")
    print("  - any closed subgroup H containing D_6 has an element of order 6, so H is")
    print("    C_6k, D_6k, SO(2)_z, O(2)_z or SO(3) (T,O,I have no order-6 element).")
    print("    C_12 fails (rotation by pi/6 sends the Q=+-6 components to their negatives),")
    print("    SO(2)_z/O(2)_z force all Q != 0 components to vanish, SO(3) forces rho_6 = 0.")
    print("    Hence Stab(rho_6(h)) = D_6, which is also the projective stabiliser of h.")

    # numerical confirmation: scan rotations, find those fixing rho_6(h)
    Jz = np.diag([float(m) for m in ORDER])
    Jp = np.zeros((7, 7), dtype=complex)
    Jm = np.zeros((7, 7), dtype=complex)
    for k, m in enumerate(ORDER):
        if m + 1 <= 3:
            Jp[ORDER.index(m + 1), k] = np.sqrt(12 - m * (m + 1))
        if m - 1 >= -3:
            Jm[ORDER.index(m - 1), k] = np.sqrt(12 - m * (m - 1))
    Jx, Jy = (Jp + Jm) / 2, (Jp - Jm) / 2j

    def D(n, t):
        A = n[0] * Jx + n[1] * Jy + n[2] * Jz
        w, V = np.linalg.eigh(A)
        return V @ np.diag(np.exp(-1j * t * w)) @ V.conj().T

    T6n = [np.array(sp.matrix2numpy(M, dtype=complex)) for M in TK(6)]
    hn = np.array([complex(x) for x in hexs]).ravel()

    def r6v(u):
        return np.array([np.einsum("i,ij,j->", u.conj(), M, u) for M in T6n])

    ref = r6v(hn)
    rng = np.random.default_rng(9)
    found = 0
    for _ in range(300000):
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        t = rng.uniform(0, 2 * np.pi)
        if np.linalg.norm(r6v(D(n, t) @ hn) - ref) < 1e-6:
            found += 1
    print("\n  random rotations (3e5) fixing rho_6(h) to 1e-6: %d" % found)
    print("  (D_6 is a measure-zero subset, so 0 is the expected count;")
    print("   the point of the scan is that no CONTINUOUS family shows up.)")
    # explicit check of the 12 elements of D_6
    ok = []
    for k in range(6):
        ok.append(np.linalg.norm(r6v(D([0, 0, 1], k * np.pi / 3) @ hn) - ref))
    for k in range(6):
        g = D([0, 0, 1], k * np.pi / 3) @ D([1, 0, 0], np.pi)
        ok.append(np.linalg.norm(r6v(g @ hn) - ref))
    print("  the 12 elements of D_6: max deviation %.3e" % max(ok))
    print("  rotation by pi/6 about z: deviation %.6f"
          % np.linalg.norm(r6v(D([0, 0, 1], np.pi / 6) @ hn) - ref))


if __name__ == "__main__":
    main()
