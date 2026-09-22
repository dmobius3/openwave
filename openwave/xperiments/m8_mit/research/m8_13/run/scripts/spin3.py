"""Spin-3 setup: operators, time reversal, the top multipole, the nematic tensor.

Everything exact (sympy). Written for the maintainer's own adjudication of the
uniqueness question, independently of the argument under audit and of the room.
"""

import sympy as sp
from sympy.physics.quantum.cg import CG

J = 3
MS = list(range(J, -J - 1, -1))  # 3, 2, 1, 0, -1, -2, -3
IDX = {m: i for i, m in enumerate(MS)}
N = len(MS)  # 7


def jz():
    return sp.diag(*[sp.Integer(m) for m in MS])


def jplus():
    M = sp.zeros(N, N)
    for m in MS:
        if m + 1 in IDX:
            M[IDX[m + 1], IDX[m]] = sp.sqrt(J * (J + 1) - m * (m + 1))
    return M


def jminus():
    return jplus().T


JZ = jz()
JP = jplus()
JM = jminus()
JX = (JP + JM) / 2
JY = (JP - JM) / (2 * sp.I)
FVEC = (JX, JY, JZ)


def theta(c, phase=lambda m: (-1) ** (J - m)):
    """(Theta u)_m = (-1)^(J-m) conj(c_{-m}).  `phase` is the arm hook."""
    return [sp.expand(phase(m) * sp.conjugate(c[IDX[-m]])) for m in MS]


def cg_table(sign=None):
    """<3 m1; 3 m2 | 6 Q>, Condon-Shortley.  `sign` is the arm hook: a callable
    (m1, m2, Q) -> +1/-1 multiplying the coefficient."""
    tab = {}
    for Q in range(-6, 7):
        for m1 in MS:
            m2 = Q - m1
            if m2 not in IDX:
                continue
            v = CG(J, m1, J, m2, 6, Q).doit()
            v = sp.nsimplify(sp.simplify(v))
            if sign is not None:
                v = v * sign(m1, m2, Q)
            tab[(m1, m2, Q)] = v
    return tab


CG6 = cg_table()


def rho6(c, tab=None, phase=None):
    tab = CG6 if tab is None else tab
    tu = theta(c) if phase is None else theta(c, phase)
    out = {}
    for Q in range(-6, 7):
        s = sp.Integer(0)
        for m1 in MS:
            m2 = Q - m1
            if m2 not in IDX:
                continue
            s += tab[(m1, m2, Q)] * c[IDX[m1]] * tu[IDX[m2]]
        out[Q] = sp.expand(s)
    return out


def norm2(c):
    return sp.expand(sum(sp.conjugate(x) * x for x in c))


def r6(c, tab=None, phase=None):
    r = rho6(c, tab=tab, phase=phase)
    num = sum(sp.expand(sp.conjugate(v) * v) for v in r.values())
    return sp.simplify(sp.expand(num) / sp.expand(norm2(c) ** 2))


def braket(c, Mat):
    """<u, Mat u> with u the column of coefficients c."""
    u = sp.Matrix(N, 1, list(c))
    return sp.expand((u.H * Mat * u)[0, 0])


def magnetization(c):
    """f = <u, f u> (unit u assumed by the caller); returns |f|^2 * ||u||^4 free
    form: the caller divides.  Here: the vector itself."""
    return [braket(c, A) for A in FVEC]


def a00(c):
    """<Theta u, u> / sqrt(7)."""
    tu = theta(c)
    return sp.expand(sum(sp.conjugate(tu[i]) * c[i] for i in range(N))) / sp.sqrt(7)


def nbar(c):
    """N_ij = Re <u, f_i f_j u>."""
    M = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            M[i, j] = sp.re(sp.expand(braket(c, FVEC[i] * FVEC[j])))
    return M


def nematic_operator(a, b):
    """a (f_z^2 - 4) + b (f_+^2 + f_-^2)."""
    return sp.expand(a * (JZ**2 - 4 * sp.eye(N)) + b * (JP**2 + JM**2))


def diag_operator(e1, e2, e3):
    """e1 f_x^2 + e2 f_y^2 + e3 f_z^2."""
    return sp.expand(e1 * JX**2 + e2 * JY**2 + e3 * JZ**2)


def blocks(a, b, s60=None, s240=None):
    s60 = sp.sqrt(60) * b if s60 is None else s60
    s240 = sp.sqrt(240) * b if s240 is None else s240
    M1 = sp.Matrix([[5 * a, s60], [s60, -3 * a + 12 * b]])
    M2 = sp.Matrix([[5 * a, -s60], [-s60, -3 * a - 12 * b]])
    M3 = sp.Matrix([[0, s240], [s240, -4 * a]])
    return {"M1": M1, "M2": M2, "M3": M3}


def evec(**kw):
    """Basis-coefficient vector from keyword m-labels, e.g. evec(p3=1, m3=1)."""
    c = [sp.Integer(0)] * N
    for k, v in kw.items():
        m = int(k[1:]) * (1 if k[0] == "p" else -1)
        c[IDX[m]] = sp.nsimplify(v)
    return c
