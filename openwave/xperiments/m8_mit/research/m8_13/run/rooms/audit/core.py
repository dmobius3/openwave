"""Audit core: everything built from conventions.md, nothing taken from claims.md.

Basis order used everywhere:  index k = 0..6  <->  m = 3-k, i.e. (v3,v2,v1,v0,v-1,v-2,v-3).
"""
import numpy as np
import sympy as sp
from sympy.physics.wigner import clebsch_gordan

J = 3
DIM = 2 * J + 1
MS = [J - k for k in range(DIM)]          # [3,2,1,0,-1,-2,-3]
IDX = {m: k for k, m in enumerate(MS)}


# ---------------------------------------------------------------- operators
def ops_exact():
    """J_z, J_+, J_-, J_x, J_y as exact sympy matrices, from conventions 1.2."""
    Jz = sp.zeros(DIM, DIM)
    Jp = sp.zeros(DIM, DIM)
    Jm = sp.zeros(DIM, DIM)
    for m in MS:
        Jz[IDX[m], IDX[m]] = m
        if m + 1 in IDX:
            Jp[IDX[m + 1], IDX[m]] = sp.sqrt(12 - m * (m + 1))
        if m - 1 in IDX:
            Jm[IDX[m - 1], IDX[m]] = sp.sqrt(12 - m * (m - 1))
    Jx = (Jp + Jm) / 2
    Jy = (Jp - Jm) / (2 * sp.I)
    return Jz, Jp, Jm, Jx, Jy


_Jz, _Jp, _Jm, _Jx, _Jy = ops_exact()
JZ = np.array(_Jz.evalf(30), dtype=complex)
JX = np.array(_Jx.evalf(30), dtype=complex)
JY = np.array(_Jy.evalf(30), dtype=complex)
JOPS = [JX, JY, JZ]
JOPS_EXACT = [_Jx, _Jy, _Jz]


# ---------------------------------------------------------------- Theta
def Theta(c):
    """(Theta u)_m = (-1)^(3-m) conj(c_{-m}).  c is indexed by k with m = 3-k."""
    out = np.zeros(DIM, dtype=complex)
    for m in MS:
        out[IDX[m]] = ((-1) ** (J - m)) * np.conj(c[IDX[-m]])
    return out


def Theta_exact(c):
    out = [0] * DIM
    for m in MS:
        out[IDX[m]] = sp.Integer((-1) ** (J - m)) * sp.conjugate(c[IDX[-m]])
    return sp.Matrix(out)


# ---------------------------------------------------------------- CG / rho6
def cg_table_exact():
    """CG[Q][m1] = <3 m1; 3 Q-m1 | 6 Q>, exact."""
    tab = {}
    for Q in range(-6, 7):
        row = {}
        for m1 in MS:
            m2 = Q - m1
            if m2 in IDX:
                row[m1] = clebsch_gordan(sp.Integer(3), sp.Integer(3), sp.Integer(6),
                                         sp.Integer(m1), sp.Integer(m2), sp.Integer(Q))
        tab[Q] = row
    return tab


CG_EXACT = cg_table_exact()
CG_NUM = {Q: {m1: complex(sp.N(v, 30)) for m1, v in row.items()} for Q, row in CG_EXACT.items()}


def rho6(c):
    """rho6(u)_Q, returned as a length-13 array indexed by Q = -6..6 -> idx Q+6."""
    t = Theta(c)
    out = np.zeros(13, dtype=complex)
    for Q in range(-6, 7):
        s = 0.0 + 0j
        for m1, coef in CG_NUM[Q].items():
            s += coef * c[IDX[m1]] * t[IDX[Q - m1]]
        out[Q + 6] = s
    return out


def rhat6(c):
    n2 = np.vdot(c, c).real
    return float(np.sum(np.abs(rho6(c)) ** 2) / n2 ** 2)


def rho6_exact(c):
    t = Theta_exact(c)
    out = []
    for Q in range(-6, 7):
        s = sp.Integer(0)
        for m1, coef in CG_EXACT[Q].items():
            s += coef * c[IDX[m1]] * t[IDX[Q - m1]]
        out.append(sp.simplify(s))
    return out


def rhat6_exact(c):
    c = sp.Matrix(c)
    n2 = sp.simplify(sum(sp.Abs(x) ** 2 for x in c))
    r = rho6_exact(list(c))
    num = sp.simplify(sum(sp.Abs(x) ** 2 for x in r))
    return sp.simplify(num / n2 ** 2)


# ---------------------------------------------------------------- invariants
def fvec(c):
    """f_i = <u, f_i u> (real for Hermitian f_i)."""
    return np.array([np.vdot(c, Op @ c).real for Op in JOPS])


def Nbar(c):
    """Nbar_ij = Re <u, f_i f_j u>."""
    return np.array([[np.vdot(c, JOPS[i] @ (JOPS[j] @ c)).real for j in range(3)]
                     for i in range(3)])


def a00(c):
    """a00 = <Theta u, u>/sqrt(7)."""
    return np.vdot(Theta(c), c) / np.sqrt(7.0)


def rot(nvec, theta):
    """D^3(n, theta) = exp(-i theta n.J), numeric."""
    import scipy.linalg as sla  # noqa
    n = np.asarray(nvec, dtype=float)
    n = n / np.linalg.norm(n)
    Aop = n[0] * JX + n[1] * JY + n[2] * JZ
    w, Vv = np.linalg.eigh(Aop)
    return Vv @ np.diag(np.exp(-1j * theta * w)) @ Vv.conj().T


# ---------------------------------------------------------------- basis vecs
def basis(m):
    v = np.zeros(DIM, dtype=complex)
    v[IDX[m]] = 1.0
    return v
