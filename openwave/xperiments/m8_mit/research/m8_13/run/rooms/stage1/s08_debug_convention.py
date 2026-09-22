"""s08_debug_convention.py -- isolate the failure of CHECK 4 in s07.

s07 CHECK 4 compared rhat_6 computed from a state built by state_from_axes()
against the axes formula, and they disagreed while CHECKs 1,2,3,5 all passed.
That points at state_from_axes (a convention slip), not at the formula.  Here I

  (a) test the formula WITHOUT any axes at all:  for a random Theta-invariant
      state u, is  rhat_6 = (13/49) <|P_u|^4>_{S^3} / <|P_u|^2>_{S^3}^2 ?
      This is the actual content of the derivation.
  (b) test  |q_n(z)|^2 = (1 - (m.n)^2)/4  for the quadratic q_n, fixing the
      spinor<->Bloch convention empirically instead of guessing it.
  (c) rebuild state_from_axes with the corrected convention and redo CHECK 4.
"""

from pathlib import Path
from math import factorial
import numpy as np
import sympy as sp
from sympy.physics.quantum.cg import CG

HERE = Path(__file__).parent
ORDER = list(range(3, -4, -1))
IDX = {m: i for i, m in enumerate(ORDER)}
rng = np.random.default_rng(101)


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


def rh6(u):
    return float(sum(abs(np.einsum("i,ij,j->", u.conj(), M, u)) ** 2 for M in T6).real)


def random_real_state():
    t = rng.normal(size=7)
    c = {0: 1j * t[0]}
    for m, (a, b) in zip([1, 2, 3], [(t[1], t[2]), (t[3], t[4]), (t[5], t[6])]):
        c[m] = a + 1j * b
        c[-m] = (-1) ** (3 + m) * np.conjugate(c[m])
    u = np.array([c[m] for m in ORDER], dtype=complex)
    return u / np.linalg.norm(u)


def P_of_u(u, z):
    """P_u(z) = sum_m c_m z1^(3+m) z2^(3-m) / sqrt((3+m)!(3-m)!)"""
    tot = 0j
    for i, m in enumerate(ORDER):
        tot += u[i] * z[0] ** (3 + m) * z[1] ** (3 - m) / np.sqrt(factorial(3 + m) * factorial(3 - m))
    return tot


def rand_spinors(n):
    z = rng.normal(size=(n, 2)) + 1j * rng.normal(size=(n, 2))
    z /= np.linalg.norm(z, axis=1)[:, None]
    return z


def bloch(z):
    return np.stack([2 * np.real(np.conj(z[:, 0]) * z[:, 1]),
                     2 * np.imag(np.conj(z[:, 0]) * z[:, 1]),
                     np.abs(z[:, 0]) ** 2 - np.abs(z[:, 1]) ** 2], axis=1)


def main():
    # ---------- (a) formula without axes ----------
    print("(a) rhat_6 =? (13/49) <|P|^4>/<|P|^2>^2  over S^3, Theta-invariant states")
    Z = rand_spinors(2000000)
    for _ in range(5):
        u = random_real_state()
        P = np.array([P_of_u(u, z) for z in Z[:200000]])
        a2 = np.mean(np.abs(P) ** 2)
        a4 = np.mean(np.abs(P) ** 4)
        print("    direct %.10f    MC formula %.10f" % (rh6(u), 13 / 49 * a4 / a2 ** 2))

    # ---------- (b) the quadratic q_n ----------
    print("\n(b) which quadratic q_n has |q_n(z)|^2 = (1-(m.n)^2)/4 ?")
    Zs = rand_spinors(6)
    ms = bloch(Zs)
    n = rng.normal(size=3)
    n /= np.linalg.norm(n)
    w = n[0] + 1j * n[1]
    cands = {
        "s07 used  (w/2)z1^2 - nz z1z2 - (wbar/2) z2^2":
            lambda z: (w / 2) * z[0] ** 2 - n[2] * z[0] * z[1] - (np.conj(w) / 2) * z[1] ** 2,
        "conjugate w swapped":
            lambda z: (np.conj(w) / 2) * z[0] ** 2 - n[2] * z[0] * z[1] - (w / 2) * z[1] ** 2,
        "with +nz":
            lambda z: (w / 2) * z[0] ** 2 + n[2] * z[0] * z[1] - (np.conj(w) / 2) * z[1] ** 2,
    }
    for name, f in cands.items():
        err = max(abs(abs(f(z)) ** 2 - (1 - (m @ n) ** 2) / 4) for z, m in zip(Zs, ms))
        print("    %-50s max err %.3e" % (name, err))

    # ---------- (c) corrected state_from_axes ----------
    print("\n(c) rebuild state_from_axes, redo CHECK 4")

    def state_from_axes(ns, conj_first):
        poly = np.array([1.0 + 0j])
        for n in ns:
            w = n[0] + 1j * n[1]
            ww = np.conj(w) if conj_first else w
            # coefficients ordered [z2^2, z1 z2, z1^2]
            q = np.array([-np.conj(ww) / 2, -n[2], ww / 2])
            poly = np.convolve(poly, q)
        c = np.array([poly[j] * np.sqrt(factorial(j) * factorial(6 - j)) for j in range(7)])
        c = c[::-1]
        return c / np.linalg.norm(c)

    c12, c13, c23 = sp.symbols("c12 c13 c23", real=True)
    txt = open(HERE / "s07_R.txt").read().splitlines()
    Fav = sp.sympify(txt[0].split("= ", 1)[1])
    F2av = sp.sympify(txt[1].split("= ", 1)[1])
    R = sp.Rational(13, 49) * F2av / Fav ** 2
    Rf = sp.lambdify((c12, c13, c23), R, "numpy")

    for conj_first in (False, True):
        print("   conj_first =", conj_first)
        for _ in range(4):
            ns = [v / np.linalg.norm(v) for v in rng.normal(size=(3, 3))]
            u = state_from_axes(ns, conj_first)
            th = np.array([(-1.0) ** (3 - m) * np.conjugate(u[IDX[-m]]) for m in ORDER])
            a = rh6(u)
            b = float(Rf(ns[0] @ ns[1], ns[0] @ ns[2], ns[1] @ ns[2]))
            print("      |<u,Theta u>| %.6f  direct %.12f  formula %.12f  diff %.2e"
                  % (abs(np.vdot(u, th)), a, b, abs(a - b)))


if __name__ == "__main__":
    main()
