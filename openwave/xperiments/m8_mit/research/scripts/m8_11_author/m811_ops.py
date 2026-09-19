"""M8.11 fresh derivation: section algebra on the flat bundle E_sigma over S^3/2I, levels 0 to 18.

A section is a list of terms (j, X, y) meaning g -> X R^j(g) y, with R(g) = D^j(g)^dagger, X a 7 x (2j+1)
matrix and y in V_j. The level-6 block in sector sigma is u -> (3, P, u), P the sigma-isotypic projector
on V_3, with int |Phi_u|^2 = (d/7) ||u||^2. The product, fibre and norm routines are M8.10's
(m810_gates.py), reused unchanged in substance; m810_core.py is a byte-identical copy of the file pinned
in #546 (sha256 9339d871...). Imports numpy and that file only.
"""
import numpy as np
from m810_core import binary_icosahedral, conj_matrix, coupling, isotypic_projector, DIM

G = binary_icosahedral()
W = {j: conj_matrix(j) for j in range(10)}
W6 = {'3p': 28 / 39, '4': 21 / 52}          # Q_sigma = 1 + w6(sigma) ||rho_6||^2 (paper Thm 5.1, M8.1.2 D7)
PROJ = {s: isotypic_projector(3, s, G) for s in ('3p', '4')}
_cp = {}


def cp(j1, j2, J):
    key = (j1, j2, J)
    if key not in _cp:
        _cp[key] = coupling(j1, j2, J)
    return _cp[key]


def channels(a, b, c, Jout):
    return [cp(J12, c, Jout) @ np.kron(cp(a, b, J12), np.eye(2 * c + 1))
            for J12 in range(abs(a - b), a + b + 1) if abs(J12 - c) <= Jout <= J12 + c]


def conj_term(t):
    j, X, y = t
    return j, X.conj() @ W[j], W[j].T @ y.conj()


def product_level(tA, tB, tC, Jout):
    """Level-2Jout part of g -> [sum_k conj(A_k) B_k] C, as terms (Jout, X, y)."""
    a, XA, yA = conj_term(tA)
    b, XB, yB = tB
    c, XC, yC = tC
    Lam = np.einsum('ij,kl->kijl', XA.T @ XB, XC).reshape(7, -1)
    Y = np.kron(np.kron(yA, yB), yC)
    return [(Jout, Lam @ Ch.T, Ch @ Y) for Ch in channels(a, b, c, Jout)]


def norm2(terms):
    """int |sum of terms|^2 for terms at one level (Schur orthogonality)."""
    J = terms[0][0]
    return sum(np.trace(Xs.conj().T @ Xt) * np.vdot(ys, yt) for _, Xt, yt in terms for _, Xs, ys in terms).real / (2 * J + 1)


def fibre(terms, P, d):
    """Fibre vector w of a level-6 section (it equals (3, P, w)), and the residual of that identification."""
    T = sum(np.einsum('ka,b->kab', X, y) for _, X, y in terms)
    w = np.array([np.trace(T[:, :, b] @ P) for b in range(7)]) / d
    return w, max(np.abs(T[:, :, b] - w[b] * P).max() for b in range(7))


def den(n):
    return n * (n + 2) - 48


def block(u, P):
    return (3, P.astype(complex), np.asarray(u, complex))


def pair(a, b, d):
    """Hermitian pairing int Phi_a^dagger Phi_b of two block sections given by their fibre vectors."""
    return d / 7 * np.vdot(a, b)


def cubic_levels(Phi):
    return {J: product_level(Phi, Phi, Phi, J) for J in range(10)}


def xi_terms(cubic, g=1.0):
    """The range correction xi = -g A^{-1} Pi_perp N(Phi), A = -Delta - 48, level by level."""
    return [(J, X, -g * y / den(2 * J)) for J, ts in cubic.items() if J != 3 for (_, X, y) in ts]


LEVELS = {'3p': {6, 10, 14, 16, 18}, '4': {6, 8, 12, 14, 16, 18}}


def forbidden_norms(cubic, sector):
    """Norms of the cubic's parts at levels outside the sector's support; they vanish exactly (F1's support gate)."""
    return {2 * J: (norm2(ts) if ts else 0.0) for J, ts in cubic.items() if 2 * J not in LEVELS[sector]}


def xi_terms_supported(cubic, sector, g=1.0):
    """xi restricted to the sector's own levels above 6; used only after forbidden_norms has been gated."""
    return [(J, X, -g * y / den(2 * J)) for J, ts in cubic.items() if J != 3 and 2 * J in LEVELS[sector] for (_, X, y) in ts]


def DN(Phi, ts, conj_term_on=True):
    """Level-6 terms of DN_Phi[t] = (Phi^dagger t) Phi + |Phi|^2 t + (t^dagger Phi) Phi, summed over the terms of t."""
    out = []
    for t in ts:
        out += product_level(Phi, t, Phi, 3) + product_level(Phi, Phi, t, 3)
        if conj_term_on:
            out += product_level(t, Phi, Phi, 3)
    return out
