"""C8 as an exact polynomial identity.

Trick: treat c_a and d_a (standing for conj(c_a)) as 14 INDEPENDENT symbols.
Every object in C8 is the natural polarization of a bidegree-(2,2) form, so if
the difference vanishes as a polynomial in (c,d) the identity holds for every
u in C^7 (and the polarized statement is strictly stronger).
"""
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import sympy as sp
import core as C

HERE = pathlib.Path(__file__).parent
out = []
P = out.append

c = sp.symbols('c0:7')
d = sp.symbols('d0:7')

Jx, Jy, Jz = C.JOPS_EXACT
Jops = [Jx, Jy, Jz]


def sesq(M):
    """<u, M u>  ->  sum_a d_a M_ab c_b."""
    return sp.expand(sum(d[a] * M[a, b] * c[b] for a in range(7) for b in range(7)))


# ||u||^2
nrm2 = sp.expand(sum(d[a] * c[a] for a in range(7)))

# f_i = <u, f_i u>
fv = [sesq(M) for M in Jops]
f2 = sp.expand(sum(x ** 2 for x in fv))

# Nbar_ij = Re<u, f_i f_j u> = <u, {f_i,f_j}/2 u>
Nb = [[sesq((Jops[i] * Jops[j] + Jops[j] * Jops[i]) / 2) for j in range(3)] for i in range(3)]
TrN2 = sp.expand(sum(Nb[i][j] * Nb[j][i] for i in range(3) for j in range(3)))

# a00 = <Theta u, u>/sqrt 7 ;  (Theta u)_m = (-1)^(3-m) conj(c_{-m}) -> eps_m d_{-m}
a00 = sp.expand(sum(sp.Integer((-1) ** (3 - m)) * c[C.IDX[-m]] * c[C.IDX[m]] for m in C.MS)) / sp.sqrt(7)
a00b = sp.expand(sum(sp.Integer((-1) ** (3 - m)) * d[C.IDX[-m]] * d[C.IDX[m]] for m in C.MS)) / sp.sqrt(7)
absa2 = sp.expand(a00 * a00b)

# rho6 and ||rho6||^2
rho, rhob = {}, {}
for Q in range(-6, 7):
    s = sp.Integer(0)
    sb = sp.Integer(0)
    for m1, g in C.CG_EXACT[Q].items():
        m2 = Q - m1
        eps = sp.Integer((-1) ** (3 - m2))
        s += g * eps * c[C.IDX[m1]] * d[C.IDX[-m2]]
        sb += g * eps * d[C.IDX[m1]] * c[C.IDX[-m2]]
    rho[Q], rhob[Q] = sp.expand(s), sp.expand(sb)
R6 = sp.expand(sum(rho[Q] * rhob[Q] for Q in range(-6, 7)))

lhs = R6
rhs = sp.expand(-sp.Rational(5, 231) * nrm2 ** 2 - f2 / 22 + sp.Rational(7, 11) * absa2
                + TrN2 / sp.Integer(198))
diff = sp.expand(sp.simplify(sp.expand(lhs - rhs)))
P("C8 exact polarized test")
P(f"  lhs monomial count  = {len(lhs.as_ordered_terms())}")
P(f"  rhs monomial count  = {len(rhs.as_ordered_terms())}")
P(f"  lhs - rhs           = {diff}")
P(f"  is zero             = {diff == 0}")

# an independent numerical spot check at awkward points, incl. exact zeros,
# repeated entries, complex phases, near-degenerate
import numpy as np
from grad import f_and_grad_real
rng = np.random.default_rng(5)
tests = {
    "v3": [1, 0, 0, 0, 0, 0, 0],
    "v0": [0, 0, 0, 1, 0, 0, 0],
    "(v3+v-3)/r2": [1, 0, 0, 0, 0, 0, 1],
    "all ones": [1, 1, 1, 1, 1, 1, 1],
    "repeated+zeros": [2, 0, 2, 0, 2, 0, 2],
    "phases": [1, 1j, -1, -1j, 1, 1j, -1],
    "one tiny": [1, 1e-13, 0, 0, 0, 0, 1],
    "near-coherent": [1, 1e-7, 1e-7, 0, 0, 0, 0],
    "cplx generic": list(rng.normal(size=7) + 1j * rng.normal(size=7)),
    "v2 only": [0, 1, 0, 0, 0, 0, 0],
    "v1+v-1": [0, 0, 1, 0, 1, 0, 0],
    "v2+v-2": [0, 1, 0, 0, 0, 1, 0],
}
P("")
P("C8 numeric spot checks (lhs = rhat6*||u||^4, rhs = the C8 expression)")
P(f"{'point':<16}{'lhs':>24}{'rhs':>24}{'diff':>12}")
for nm, vec in tests.items():
    u = np.array(vec, dtype=complex)
    n2 = float(np.vdot(u, u).real)
    L = C.rhat6(u) * n2 ** 2
    fq = float(np.dot(C.fvec(u), C.fvec(u)))
    Nq = C.Nbar(u)
    R = (-5 / 231) * n2 ** 2 - fq / 22 + (7 / 11) * abs(C.a00(u)) ** 2 + float(np.sum(Nq * Nq)) / 198
    P(f"{nm:<16}{L:>24.15f}{R:>24.15f}{L-R:>12.2e}")

txt = "\n".join(map(str, out))
print(txt)
(HERE / "out_s2.txt").write_text(txt + "\n")
