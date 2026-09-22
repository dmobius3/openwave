"""Analytic gradient of rhat6, plus a finite-difference verification."""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import numpy as np
import core as C

# rho_Q = sum_{m1} T[Q][m1] * c_{m1} * conj(c_{m1-Q}),  T real
T = np.zeros((13, 7))          # T[Q+6, idx(m1)]
for Q in range(-6, 7):
    for m1, g in C.CG_NUM[Q].items():
        m2 = Q - m1
        eps = (-1.0) ** (3 - m2)
        # (Theta u)_{m2} = eps * conj(c_{-m2}),  -m2 = m1 - Q
        T[Q + 6, C.IDX[m1]] = g.real * eps


def rho_T(c):
    out = np.zeros(13, dtype=complex)
    for Q in range(-6, 7):
        s = 0j
        for m1 in C.MS:
            if (m1 - Q) in C.IDX and T[Q + 6, C.IDX[m1]] != 0.0:
                s += T[Q + 6, C.IDX[m1]] * c[C.IDX[m1]] * np.conj(c[C.IDX[m1 - Q]])
        out[Q + 6] = s
    return out


def N_and_grad(c):
    """N = sum_Q |rho_Q|^2 ; returns N and dN/dconj(c) (Wirtinger)."""
    r = rho_T(c)
    N = float(np.sum(np.abs(r) ** 2))
    g = np.zeros(7, dtype=complex)
    for Q in range(-6, 7):
        rq = r[Q + 6]
        for m1 in C.MS:
            m2i = m1 - Q
            if m2i in C.IDX and T[Q + 6, C.IDX[m1]] != 0.0:
                t = T[Q + 6, C.IDX[m1]]
                # d rho_Q / d conj(c_{m1-Q}) = t * c_{m1}
                g[C.IDX[m2i]] += t * c[C.IDX[m1]] * np.conj(rq)
                # d conj(rho_Q) / d conj(c_{m1}) = t * c_{m1-Q}
                g[C.IDX[m1]] += rq * t * c[C.IDX[m2i]]
    return N, g


def f_and_grad_real(x):
    """x in R^14 -> (rhat6, gradient in R^14)."""
    c = x[:7] + 1j * x[7:]
    n2 = float(np.vdot(c, c).real)
    N, gN = N_and_grad(c)
    val = N / n2 ** 2
    gW = gN / n2 ** 2 - 2.0 * N * c / n2 ** 3      # d/dconj(c)
    gr = 2.0 * gW                                  # real grad wrt (Re,Im)
    return val, np.concatenate([gr.real, gr.imag])


if __name__ == "__main__":
    rng = np.random.default_rng(7)
    worst_v, worst_g = 0.0, 0.0
    for _ in range(50):
        x = rng.normal(size=14)
        c = x[:7] + 1j * x[7:]
        worst_v = max(worst_v, abs(rho_T(c) - C.rho6(c)).max())
        v, g = f_and_grad_real(x)
        worst_v = max(worst_v, abs(v - C.rhat6(c)))
        h = 1e-6
        gfd = np.array([(f_and_grad_real(x + h * np.eye(14)[k])[0]
                         - f_and_grad_real(x - h * np.eye(14)[k])[0]) / (2 * h) for k in range(14)])
        worst_g = max(worst_g, np.abs(g - gfd).max())
    print(f"rho/value agreement with core : {worst_v:.3e}")
    print(f"gradient vs finite difference : {worst_g:.3e}")
