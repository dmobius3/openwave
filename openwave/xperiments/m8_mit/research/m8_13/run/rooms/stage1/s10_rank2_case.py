"""s10_rank2_case.py -- the case that the convexity argument leaves open.

The convexity argument shows max rhat_6 is ATTAINED at a rank-one S = x x^T,
i.e. at a Theta-invariant state.  It does NOT by itself show that every
maximiser is of that form.  Writing S = sum_i w_i x_i x_i^T (w_i > 0, sum w_i = 1,
x_i orthonormal real, at most two terms), and f(S) = ||Pi_6 S||^2,

    f(S) = || sum_i w_i Pi_6(x_i x_i^T) ||^2 <= ( sum_i w_i ||Pi_6(x_i x_i^T)|| )^2
                                             <= max_i f(x_i x_i^T)  <= M,

with equality throughout iff every x_i is a rank-one maximiser AND all the
vectors Pi_6(x_i x_i^T) are EQUAL.  So rank-two maximisers exist iff there are
two ORTHOGONAL real unit maximisers x, y with Pi_6(x x^T) = Pi_6(y y^T),
equivalently x x^T - y y^T in W_2 (+) W_4.

This script tests that numerically:
  (1) are all 60 stored numerical maximisers Theta-invariant up to a phase?
  (2) over the hexagon orbit, is there a pair x,y with x.y = 0 and
      ||Pi_6(xx^T) - Pi_6(yy^T)|| = 0 ?
"""

from pathlib import Path
import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG
from scipy.optimize import minimize

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))
IDX = {m: i for i, m in enumerate(ORDER)}


def T6_num():
    out = []
    for Q in range(-6, 7):
        M = np.zeros((7, 7), dtype=complex)
        for i, mm in enumerate(ORDER):
            for k, m1 in enumerate(ORDER):
                if m1 - mm == Q:
                    M[i, k] = float(sp.N(CG(3, m1, 3, -mm, 6, Q).doit(), 25)) * (-1.0) ** (3 + mm)
        out.append(M)
    return out


T6 = T6_num()


def theta(u):
    return np.array([(-1.0) ** (3 - m) * np.conjugate(u[IDX[-m]]) for m in ORDER])


def rh6(u):
    return float(sum(abs(np.einsum("i,ij,j->", u.conj(), M, u)) ** 2 for M in T6).real)


def rho6vec(u):
    """the 13 components rho_6(u)_Q = Pi_6 image, as a complex vector"""
    return np.array([np.einsum("i,ij,j->", u.conj(), M, u) for M in T6])


def main():
    # ---- (1) all stored maximisers ----
    top = np.load(HERE / "s03_top.npy")
    ovs = [abs(np.vdot(u, theta(u))) / np.vdot(u, u).real for u in top]
    vals = [rh6(u / np.linalg.norm(u)) for u in top]
    print("(1) %d stored near-maximisers" % len(top))
    print("    rhat_6 range: %.15f .. %.15f" % (min(vals), max(vals)))
    print("    |<u,Theta u>|/||u||^2 : min %.15f  max %.15f" % (min(ovs), max(ovs)))
    print("    -> all Theta-invariant up to phase" if min(ovs) > 1 - 1e-9 else "    -> SOME ARE NOT")

    bot = np.load(HERE / "s03_bot.npy")
    ovsb = [abs(np.vdot(u, theta(u))) / np.vdot(u, u).real for u in bot]
    print("    (minimisers for contrast: |<u,Theta u>| in %.4f .. %.4f)" % (min(ovsb), max(ovsb)))

    # ---- (2) two orthogonal real maximisers with equal Pi_6 image? ----
    # real form: 7 real parameters
    def real_state(t):
        cc = {0: 1j * t[0]}
        for m, (p, q) in zip([1, 2, 3], [(t[1], t[2]), (t[3], t[4]), (t[5], t[6])]):
            cc[m] = p + 1j * q
            cc[-m] = (-1) ** (3 + m) * np.conjugate(cc[m])
        u = np.array([cc[m] for m in ORDER], dtype=complex)
        return u / np.linalg.norm(u)

    M = 463 / 924

    def obj(z):
        x = real_state(z[:7])
        y = real_state(z[7:])
        pen = 0.0
        pen += (rh6(x) - M) ** 2 + (rh6(y) - M) ** 2
        pen += abs(np.vdot(x, y)) ** 2          # orthogonality in the real form
        pen += np.sum(np.abs(rho6vec(x) - rho6vec(y)) ** 2)
        return pen

    print("\n(2) searching for orthogonal real maximisers x,y with equal Pi_6 image")
    rng = np.random.default_rng(31)
    best = (1e9, None)
    for _ in range(300):
        z0 = rng.normal(size=14)
        r = minimize(obj, z0, method="Nelder-Mead",
                     options={"maxiter": 40000, "maxfev": 40000, "xatol": 1e-12, "fatol": 1e-16})
        if r.fun < best[0]:
            best = (r.fun, r.x)
    print("    best residual over 300 starts: %.3e" % best[0])
    x = real_state(best[1][:7])
    y = real_state(best[1][7:])
    print("    rhat6(x)=%.12f  rhat6(y)=%.12f  |<x,y>|=%.3e  ||drho||=%.3e"
          % (rh6(x), rh6(y), abs(np.vdot(x, y)), np.linalg.norm(rho6vec(x) - rho6vec(y))))

    # ---- (2b) drop orthogonality: are two DISTINCT real maximisers with the same
    #          Pi_6 image possible at all (the phases/rotations aside)? ----
    def obj2(z):
        x = real_state(z[:7])
        y = real_state(z[7:])
        pen = (rh6(x) - M) ** 2 + (rh6(y) - M) ** 2
        pen += np.sum(np.abs(rho6vec(x) - rho6vec(y)) ** 2)
        pen += 1.0 / (1e-6 + min(abs(np.vdot(x, y)) ** 2 - 1, 0) ** 2 * 0 + 1)  # placeholder
        return pen

    print("\n(2b) same Pi_6 image without the orthogonality demand:")
    print("     rho_6(x) determines rhat_6 = ||rho_6||^2, so equality of images is the")
    print("     only obstruction; if x,y real maximisers have rho_6(x)=rho_6(y) then")
    print("     any convex combination of xx^T,yy^T is also a maximiser.")
    # direct test: is the map x -> rho_6(x) injective on the hexagon orbit (mod sign)?
    hexs = np.zeros(7, dtype=complex)
    hexs[IDX[3]] = 1 / np.sqrt(2)
    hexs[IDX[-3]] = 1 / np.sqrt(2)
    print("     rho_6(hexagon) =", np.round(rho6vec(hexs), 10))
    print("     ||rho_6||^2 =", np.sum(np.abs(rho6vec(hexs)) ** 2))
    # rotate the hexagon and see whether rho_6 ever returns to the same vector
    # for a DIFFERENT state: that is exactly the stabiliser D_6 of the constellation.


if __name__ == "__main__":
    main()
