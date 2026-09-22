"""s03_numeric.py -- numerical corroboration only.

What this CAN do: exhibit states with particular values of rhat_6, and show that
many independent local searches stall at the same two values.
What this CANNOT do: establish a maximum, a minimum, or completeness of an
extremiser set.  Those are settled in s04/s05/s06 by exact algebra.

Method: rhat_6(u) = sum_Q |<u|T^6_Q|u>|^2 on the unit sphere of C^7.
Gradient wrt conj(u) is H(u) u with H(u) = sum_Q (conj(rho_Q) T_Q + rho_Q T_Q^H),
Hermitian, and rhat_6 = (1/2) <u|H(u)|u>.  Iterating u <- normalise((H(u)+sI)u)
ascends; using (-H(u)+sI) descends.
"""

from pathlib import Path
import numpy as np
from numpy.linalg import norm, eigh
import json
import sympy as sp
from sympy.physics.quantum.cg import CG

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))
rng = np.random.default_rng(20260922)


def build_T(Jt):
    T = []
    for Q in range(-Jt, Jt + 1):
        M = np.zeros((7, 7), dtype=complex)
        for i, mm in enumerate(ORDER):
            for k, m1 in enumerate(ORDER):
                if m1 - mm == Q:
                    M[i, k] = complex(sp.N(sp.simplify(CG(3, m1, 3, -mm, Jt, Q).doit()), 30)) * (-1.0) ** (3 + mm)
        T.append(M)
    return np.array(T)


T6 = build_T(6)
TALL = {Jt: build_T(Jt) for Jt in range(7)}


def rho(u, Ts):
    return np.einsum("i,qij,j->q", u.conj(), Ts, u)


def f6(u):
    u = u / norm(u)
    return float(np.sum(np.abs(rho(u, T6)) ** 2).real)


def weights(u):
    u = u / norm(u)
    return {Jt: float(np.sum(np.abs(rho(u, TALL[Jt])) ** 2).real) for Jt in range(7)}


def H(u):
    r = rho(u, T6)
    return np.einsum("q,qij->ij", r.conj(), T6) + np.einsum("q,qij->ij", r, T6.conj().transpose(0, 2, 1))


def iterate(u, sign, iters=4000, shift=6.0):
    u = u / norm(u)
    for _ in range(iters):
        Hu = sign * H(u) + shift * np.eye(7)
        w = Hu @ u
        nw = norm(w)
        if nw < 1e-14:
            break
        u = w / nw
    return u


def polish(u, sign):
    """A few Newton-free refinement rounds with a tiny shift, then report."""
    for sh in (6.0, 2.0, 0.5, 0.05):
        u = iterate(u, sign, iters=6000, shift=sh)
    return u


def rot_matrix(axis, theta):
    """D^3(n,theta) = exp(-i theta n.J) on the 7-dim space, built from J_z,J_pm."""
    Jz = np.diag([float(m) for m in ORDER])
    Jp = np.zeros((7, 7), dtype=complex)
    Jm = np.zeros((7, 7), dtype=complex)
    for k, m in enumerate(ORDER):
        if m + 1 <= 3:
            i = ORDER.index(m + 1)
            Jp[i, k] = np.sqrt(12 - m * (m + 1))
        if m - 1 >= -3:
            i = ORDER.index(m - 1)
            Jm[i, k] = np.sqrt(12 - m * (m - 1))
    Jx = (Jp + Jm) / 2
    Jy = (Jp - Jm) / (2j)
    n = np.array(axis, dtype=float)
    n = n / norm(n)
    Aop = n[0] * Jx + n[1] * Jy + n[2] * Jz
    w, V = eigh(Aop)
    return V @ np.diag(np.exp(-1j * theta * w)) @ V.conj().T


def main():
    out = {}

    # ---- invariance spot checks (numerical, 1e-14 level) ----
    u = rng.normal(size=7) + 1j * rng.normal(size=7)
    u /= norm(u)
    base = f6(u)
    devs = []
    for _ in range(20):
        R = rot_matrix(rng.normal(size=3), rng.uniform(0, 2 * np.pi))
        devs.append(abs(f6(R @ u) - base))
    ph = abs(f6(np.exp(1.234j) * u) - base)
    print("invariance: max |rhat6(Ru)-rhat6(u)| over 20 random rotations = %.3e" % max(devs))
    print("            |rhat6(e^{i phi}u)-rhat6(u)| = %.3e" % ph)
    print("            D^3(n,2pi) - I  max entry = %.3e"
          % np.max(np.abs(rot_matrix([0.3, -0.7, 0.2], 2 * np.pi) - np.eye(7))))
    out["invariance_rot"] = max(devs)

    # ---- global search ----
    best, worst = [], []
    for trial in range(400):
        u0 = rng.normal(size=7) + 1j * rng.normal(size=7)
        umax = polish(u0.copy(), +1)
        umin = polish(u0.copy(), -1)
        best.append((f6(umax), umax))
        worst.append((f6(umin), umin))
    best.sort(key=lambda t: -t[0])
    worst.sort(key=lambda t: t[0])

    vmax = best[0][0]
    vmin = worst[0][0]
    print("\n400 random starts:")
    print("  max found  %.15f" % vmax)
    print("  top-10     ", ["%.12f" % b[0] for b in best[:10]])
    print("  distinct high plateaus:", sorted({round(b[0], 9) for b in best[:120]}, reverse=True)[:8])
    print("  min found  %.15f" % vmin)
    print("  bot-10     ", ["%.12f" % b[0] for b in worst[:10]])
    print("  distinct low plateaus:", sorted({round(b[0], 9) for b in worst[:120]})[:8])

    print("\n  6/7 = %.15f   1/924 = %.15f" % (6 / 7, 1 / 924))
    print("  candidate max as fraction guesses:")
    for q in [sp.nsimplify(vmax, rational=False, tolerance=1e-9),
              sp.Rational(str(round(vmax, 10))).limit_denominator(100000)]:
        print("   ", q, float(q))
    print("  candidate min limit_denominator:",
          sp.Rational(str(round(vmin, 12))).limit_denominator(1000000))

    umax = best[0][1]
    umin = worst[0][1]
    print("\n  multipole weights at the maximiser:", {k: round(v, 10) for k, v in weights(umax).items()})
    print("  multipole weights at the minimiser:", {k: round(v, 10) for k, v in weights(umin).items()})
    print("\n  maximiser components (rotated to no particular frame):")
    for k, m in enumerate(ORDER):
        print("    m=%+d  %.10f %+.10fi   |.|=%.10f" % (m, umax[k].real, umax[k].imag, abs(umax[k])))
    print("  minimiser components:")
    for k, m in enumerate(ORDER):
        print("    m=%+d  %.10f %+.10fi   |.|=%.10f" % (m, umin[k].real, umin[k].imag, abs(umin[k])))

    out["max_numeric"] = vmax
    out["min_numeric"] = vmin
    np.save(HERE / "s03_umax.npy", umax)
    np.save(HERE / "s03_umin.npy", umin)
    np.save(HERE / "s03_top.npy", np.array([b[1] for b in best[:60]]))
    np.save(HERE / "s03_bot.npy", np.array([b[1] for b in worst[:60]]))
    with open(HERE / "s03_out.json", "w") as fh:
        json.dump(out, fh, indent=2)


if __name__ == "__main__":
    main()
