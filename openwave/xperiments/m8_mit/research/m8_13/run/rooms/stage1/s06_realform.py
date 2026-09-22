"""s06_realform.py -- reduction of the MAXIMUM to the real form, and the search
for a linear relation among p_2, p_4, p_6 there.

Step 1 (proved in the return, checked numerically here):
  Theta is an antiunitary involution (Theta^2 = 1), so V3 = V_R (+) i V_R with
  V_R = {u : Theta u = u} a real 7-dim SO(3)-invariant subspace (the real spin-3
  representation).  Writing u = x + i y with x,y in V_R,
      u (x) Theta u  has symmetric part  S = x x^T + y y^T   (real, PSD, tr S = 1)
  and rhat_6 = || Pi_6 S ||^2 because K = 6 is EVEN and so lives in Sym^2.
  S ranges exactly over {PSD, trace 1, rank <= 2}.  || Pi_6 S ||^2 is a CONVEX
  quadratic form in S, so its max over that set is its max over the convex hull
  {PSD, trace 1}, attained at an extreme point = rank one = y = 0 = Theta u = u.
  Hence the maximum of rhat_6 is attained on the real form, where rhat_6 = p_6.

Step 2 (this script): on V_R, are p_2, p_4, p_6 affinely dependent?  If the space
of SO(3)-invariant quartics on V_R is 2-dimensional then they are, and positivity
of the p_J alone bounds p_6.
"""

from pathlib import Path
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


def blocks_num():
    B = {}
    for Jt in [0, 2, 4, 6]:
        vs = []
        for Q in range(-Jt, Jt + 1):
            S = np.zeros((7, 7), dtype=complex)
            for m1 in range(-3, 4):
                m2 = Q - m1
                if abs(m2) <= 3:
                    S[IDX[m1], IDX[m2]] = complex(sp.N(cgx(3, m1, 3, m2, Jt, Q), 30))
            vs.append(S)
        B[Jt] = vs
    return B


def TK_num(K):
    out = []
    for Q in range(-K, K + 1):
        M = np.zeros((7, 7), dtype=complex)
        for i, mm in enumerate(ORDER):
            for k, m1 in enumerate(ORDER):
                M[i, k] = complex(sp.N(cgx(3, m1, 3, -mm, K, Q), 30)) * (-1.0) ** (3 + mm)
        out.append(M)
    return out


BL = blocks_num()
T = {K: TK_num(K) for K in range(7)}


def theta(u):
    out = np.zeros(7, dtype=complex)
    for i, m in enumerate(ORDER):
        out[i] = (-1.0) ** (3 - m) * np.conjugate(u[IDX[-m]])
    return out


def pJ(u):
    S = np.outer(u, u)
    return {J: float(sum(abs(np.sum(np.conjugate(v) * S)) ** 2 for v in BL[J]).real) for J in BL}


def rhatK(u, K):
    return float(sum(abs(np.einsum("i,ij,j->", u.conj(), M, u)) ** 2 for M in T[K]).real)


def random_real_state(rng):
    """uniform-ish random unit vector with Theta u = u (7 real parameters)."""
    t = rng.normal(size=7)
    c = {}
    c[0] = 1j * t[0]
    for m, (a, b) in zip([1, 2, 3], [(t[1], t[2]), (t[3], t[4]), (t[5], t[6])]):
        c[m] = a + 1j * b
        c[-m] = (-1) ** (3 + m) * np.conjugate(c[m])
    u = np.array([c[m] for m in ORDER], dtype=complex)
    return u / np.linalg.norm(u)


def main():
    rng = np.random.default_rng(11)

    # --- check the real-form parametrisation really gives Theta u = u ---
    devs = []
    for _ in range(50):
        u = random_real_state(rng)
        devs.append(np.max(np.abs(theta(u) - u)))
    print("real-form parametrisation: max |Theta u - u| over 50 samples = %.3e" % max(devs))

    # --- on V_R:  rhat_6 == p_6 ?  and p_0 == 1/7 ?  and odd rhat_K == 0 ? ---
    print("\non V_R:")
    for _ in range(4):
        u = random_real_state(rng)
        p = pJ(u)
        print("   p0=%.12f p2=%.12f p4=%.12f p6=%.12f | rhat6=%.12f | odd rhat: %.2e %.2e %.2e"
              % (p[0], p[2], p[4], p[6], rhatK(u, 6), rhatK(u, 1), rhatK(u, 3), rhatK(u, 5)))

    # --- is (1, p2, p4, p6) of rank 2 on V_R?  (=> one affine relation) ---
    rows = []
    for _ in range(60):
        u = random_real_state(rng)
        p = pJ(u)
        rows.append([1.0, p[2], p[4], p[6]])
    Mx = np.array(rows)
    sv = np.linalg.svd(Mx, compute_uv=False)
    print("\nrank test for span{1, p2, p4, p6} on V_R:")
    print("   singular values:", np.array2string(sv, precision=6))
    rank = int(np.sum(sv > 1e-9 * sv[0]))
    print("   numerical rank =", rank, " => number of independent affine relations =", 4 - rank)

    # null space -> the relation(s)
    U, s, Vt = np.linalg.svd(Mx)
    for k in range(rank, 4):
        v = Vt[k]
        v = v / v[np.argmax(np.abs(v))]
        print("   relation %d :  %s . (1,p2,p4,p6) = 0   ->" % (k, np.round(v, 10)),
              "  p6 = %s" % np.round(-v[:3] / v[3], 10) if abs(v[3]) > 1e-9 else "")

    # --- same test on the FULL complex sphere (should have rank 4: no relation) ---
    rows = []
    for _ in range(80):
        u = rng.normal(size=7) + 1j * rng.normal(size=7)
        u /= np.linalg.norm(u)
        p = pJ(u)
        rows.append([1.0, p[0], p[2], p[4], p[6]])
    sv2 = np.linalg.svd(np.array(rows), compute_uv=False)
    print("\nsame test on the full complex sphere, span{1,p0,p2,p4,p6}:")
    print("   singular values:", np.array2string(sv2, precision=6))
    print("   numerical rank =", int(np.sum(sv2 > 1e-9 * sv2[0])), "(4 expected: only p0+p2+p4+p6=1)")

    # --- the hexagon state, exactly and numerically ---
    hexs = np.zeros(7, dtype=complex)
    hexs[IDX[3]] = 1 / np.sqrt(2)
    hexs[IDX[-3]] = 1 / np.sqrt(2)
    print("\nhexagon state (v_3 + v_-3)/sqrt2 :")
    print("   Theta u - u :", np.max(np.abs(theta(hexs) - hexs)))
    p = pJ(hexs)
    print("   p =", {k: round(v, 12) for k, v in p.items()})
    print("   rhat_K =", {K: round(rhatK(hexs, K), 12) for K in range(7)})
    print("   rhat_6 = %.15f    463/924 = %.15f" % (rhatK(hexs, 6), 463 / 924))

    # --- numerical maximiser from s03: is it Theta-invariant up to phase? ---
    try:
        umax = np.load(HERE / "s03_umax.npy")
        # Theta u = e^{i a} u  <=>  |<u, Theta u>| = 1
        ov = abs(np.vdot(umax, theta(umax)))
        print("\ns03 numerical maximiser: |<u, Theta u>| = %.15f (1 => Theta-invariant up to phase)" % ov)
        print("   its p =", {k: round(v, 12) for k, v in pJ(umax).items()})
        umin = np.load(HERE / "s03_umin.npy")
        print("   s03 numerical minimiser: |<u,Theta u>| = %.12f" % abs(np.vdot(umin, theta(umin))))
        print("   its p =", {k: round(v, 12) for k, v in pJ(umin).items()})
    except FileNotFoundError:
        pass


if __name__ == "__main__":
    main()
