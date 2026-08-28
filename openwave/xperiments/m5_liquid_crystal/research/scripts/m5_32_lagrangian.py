"""M5.32 R0: the covariant term REGISTRY of the 4x4 M5 Lagrangian.

Every candidate term has ONE definition here, with a sympy implementation
(the author's notebook conventions) and a numpy implementation (the
certified lattice stencil of m5_21_3_a_4d.py); the selftests prove the two
agree and that I1 reproduces the certified energies.

EQUATIONS FIRST
---------------
Field: M(x) a real symmetric 4x4 matrix per point. Internal metric
eta = diag(-1, 1, 1, 1) (index 0 = time). Jets A_mu = d_mu M, mu = 0..3
(mu = 0 the time derivative, A_0 = omega * a0 for a clock direction a0).

Curvature (the certified eta-commutator, m5_21_3_a_4d.py comm_eta):
    F_{mu nu} = A_mu eta A_nu - A_nu eta A_mu          (antisymmetric 4x4)
    F[mu, nu, a, b]: mu, nu DERIVATIVE indices (0..3), a, b INTERNAL
    indices (0..3), both antisymmetric pairs; the stored array holds the
    raw matrix entries F_{mu nu}[a, b].

Transformation law (M -> Lambda M Lambda^T with Lambda^T eta Lambda = eta,
x -> Lambda x): the matrix entries of M and F carry CONTRAVARIANT internal
indices (F -> Lambda F Lambda^T) while d_mu is COVARIANT (transforms with
Lambda^{-T} = eta Lambda eta). Therefore the invariant contraction rule is:
    derivative-derivative pair : eta^{mu nu}
    internal-internal pair     : eta_{ab}
    derivative-internal pair   : delta (Kronecker), NO eta
(equivalently: lower the internal indices with eta first, then raise
every index with eta; the "eta on all four raw indices" reading is NOT
covariant and is registered only as a control, I3_mixed_eta).

The inner product <F, G>_eta = tr(eta F eta G^T) = sum_ab eta_a eta_b F_ab G_ab.

Registered terms (densities; scalar contractions in the rule above):
    I1      = sum_{mu<nu} eta^mu eta^nu <F_mu nu, F_mu nu>_eta
            = (1/2) F_{mu nu a b} F^{mu nu a b}       (the certified action)
    I1_frob = sum_{mu<nu} eta^mu eta^nu tr(F_mu nu F_mu nu^T)
              (Frobenius internal metric; the M5.21.16 variant-A flip; control)
    I2      = F_{mu nu a b} F^{a b mu nu}
    I3      = F_{mu nu a b} F^{mu a nu b}
    I3_mixed_eta = I3 with eta on the mixed pairs too (non-covariant control)
    R_nu^a  = F_{mu nu}^{a mu}  (the derivative-internal trace; the four
              slot pairings R_reading_1..4 are enumerated in the selftest)
    I4      = R_{nu a} R^{nu a} = sum eta_nu eta_a R[nu,a]^2
    I5      = R_{nu a} R^{a nu} = sum R[nu,a] R[a,nu]
    I6      = R^2,  R = R_nu^nu = sum_nu R[nu,nu]
    V4      = w * sum_{p=1..4} (tr((M eta)^p) - C_p)^2,
              C_p = (s g)^p + 1 + delta^p        (the certified potential)

Vacuum (branch s): M_vac = diag(-s g, 1, delta, 0), (M eta) spectrum
(s g, 1, delta, 0). The author's notebook frame d = diag(g, 1, delta, 0)
with xi = diag(-1,1,1,1) is the s = -1 branch (M_vac = d at s = -1).

Legendre read (the Hamiltonian of a term): with A_0 = omega a0 every
registered term is exactly quadratic in omega, I(omega) = A + B omega +
C omega^2, and its Hamiltonian density is H_I = omega dI/domega - I =
C omega^2 - A. The certified functional is the Lagrangian
    L_cert = -4 * I1 - w * V4/w   ->   E_cert = 4 (U + omega^2 T) + V4
i.e. E_u = 4 h^3 sum_{i<j} <F_ij, F_ij>_eta,  kin = 4 h^3 sum_i <F_0i,F_0i>_eta
(m5_21_3_a_4d.py e_parts / kin_of). The FACTOR-4 BRIDGE: the notebook
H = sum_{mu<nu} <F,F>_eta carries no 4; the lattice energies carry 4 h^3.

Stencil: the certified sym stencil, 1/2(fwd + bwd) with the density
evaluated per branch and averaged (NOT the averaged derivative), exact
adjoints d1_adj, h = L/n, no pin in the energy.

SELFTESTS (python3 m5_32_lagrangian.py --selftest [--mutant eta_time_row]):
each a PASS line that can fail; the mutant flips the internal time-row
sign of the CONTRACTION metric (eta_int -> diag(+1,1,1,1)) in both
implementations and must redden at least one line.

Out: ../data/m5_32_r0_module_selftest.json
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import os
import subprocess
import sys
import time

import numpy as np
import sympy as sp
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PY = sys.executable


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")     # the certified stack

# ---------------- conventions (single source) ----------------
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])       # bracket metric (never mutated)
ETA_INT = np.diag([-1.0, 1.0, 1.0, 1.0])   # contraction metric (mutable)
XI = sp.diag(-1, 1, 1, 1)                  # notebook xi (bracket)
XI_INT = sp.diag(-1, 1, 1, 1)              # notebook contraction metric
W1 = B3.W1                                 # certified potential weight
MUTANT = None
# the stored certified field used by the selftests (the M5.21.3 seed npz
# is no longer on disk; this M5.21.11 end state has its E_u on record)
STORED3_NPZ = "m5_21_11_end_t11lad_A_n32_d0.3.npz"
STORED3_ROW = "m5_21_11_row_t11lad_A_n32_d0.3.json"


def load_stored3():
    return np.load(os.path.join(DATA, STORED3_NPZ))["M"].astype(np.float64)


def set_mutant(name):
    """eta_time_row: flip the internal time-row sign of the CONTRACTION
    metric in both implementations (the bracket keeps eta)."""
    global ETA_INT, XI_INT, MUTANT
    if name == "eta_time_row":
        ETA_INT = np.diag([1.0, 1.0, 1.0, 1.0])
        XI_INT = sp.diag(1, 1, 1, 1)
        MUTANT = name
    elif name:
        raise ValueError(name)


# ================= numpy layer =================
def F_of_A(A):
    """A: (4, ..., 4, 4) jets -> F: (..., 4, 4, 4, 4) = F[mu, nu, a, b]
    with F_{mu nu} = A_mu eta A_nu - A_nu eta A_mu (B3.comm_eta)."""
    AE = A @ ETA                                  # A_mu eta
    P = np.einsum("m...ab,n...bc->...mnac", AE, A, optimize=True)
    return P - P.swapaxes(-4, -3)


def _K_from_pattern(p1, p2, factor, int_metric="eta", mixed_metric="delta",
                    deriv_metric="eta"):
    """The 256x256 contraction matrix K such that
    density = sum F[i1..i4] K[i1..i4, j1..j4] F[j1..j4].
    p1, p2: 4-letter slot patterns (slots 0,1 derivative; 2,3 internal);
    a repeated letter is contracted with the metric selected by the slot
    types of its two occurrences."""
    slots = [(l, 1, k) for k, l in enumerate(p1)] + \
            [(l, 2, k) for k, l in enumerate(p2)]
    letters = {}
    for l, which, k in slots:
        letters.setdefault(l, []).append((which, k))
    metric = {}
    for l, occ in letters.items():
        assert len(occ) == 2, f"letter {l} must appear exactly twice"
        kinds = tuple(sorted("d" if k < 2 else "i" for _, k in occ))
        rule = {("d", "d"): deriv_metric, ("i", "i"): int_metric,
                ("d", "i"): mixed_metric}[kinds]
        if rule != "eta":
            metric[l] = np.ones(4)
        elif kinds == ("i", "i"):
            metric[l] = np.diag(ETA_INT)      # internal pair (mutable)
        else:
            metric[l] = np.diag(ETA)          # derivative or mixed pair
    K = np.zeros((4,) * 8)
    for vals in itertools.product(range(4), repeat=len(letters)):
        env = dict(zip(letters.keys(), vals))
        w = factor
        for l, v in env.items():
            w *= metric[l][v]
        i = tuple(env[l] for l in p1)
        j = tuple(env[l] for l in p2)
        K[i + j] += w
    return K.reshape(256, 256)


def density_from_K(F, K):
    f = F.reshape(F.shape[:-4] + (256,))
    return np.einsum("...i,ij,...j->...", f, K, f, optimize=True)


def dW_from_K(F, K):
    """W[mu nu a b] = d density / d F[mu nu a b] = (K + K^T) f."""
    f = F.reshape(F.shape[:-4] + (256,))
    W = f @ (K + K.T)
    return W.reshape(F.shape)


def dA_from_W(W, A):
    """chain W (d density/dF) back to dA_mu (all ordered (mu,nu) pairs):
    dF_{mu nu} = dA_mu eta A_nu + A_mu eta dA_nu - dA_nu eta A_mu
                 - A_nu eta dA_mu."""
    EA_T = (ETA @ A).swapaxes(-1, -2)            # (eta A_nu)^T
    AE_T = (A @ ETA).swapaxes(-1, -2)            # (A_nu eta)^T
    G = np.zeros_like(A)
    for mu in range(4):
        for nu in range(4):
            Wmn = W[..., mu, nu, :, :]
            G[mu] += Wmn @ EA_T[nu] - AE_T[nu] @ Wmn
            G[nu] += AE_T[mu] @ Wmn - Wmn @ EA_T[mu]
    return G


def R_readings_np(F):
    """the four derivative-internal slot pairings of the trace, plus the
    two within-type eta-traces (expected to vanish)."""
    out = {}
    out["R_reading_1 (slot0,slot3)"] = np.einsum("...mnam->...na", F)
    out["R_reading_2 (slot0,slot2)"] = np.einsum("...mnmb->...nb", F)
    out["R_reading_3 (slot1,slot3)"] = np.einsum("...mnan->...ma", F)
    out["R_reading_4 (slot1,slot2)"] = np.einsum("...mnnb->...mb", F)
    out["trace_deriv_pair eta^{mu nu} F_{mu nu ab}"] = \
        np.einsum("m,...mmab->...ab", np.diag(ETA), F)
    out["trace_int_pair eta_{ab} F_{mu nu ab}"] = \
        np.einsum("a,...mnaa->...mn", np.diag(ETA_INT), F)
    return out


def v4_traces_np(M):
    Me = M @ ETA
    P = Me
    t = []
    for p in range(4):
        if p:
            P = P @ Me
        t.append(np.einsum("...kk->...", P))
    return t


def c4_of(p):
    sg = p["s"] * p["g"]
    return tuple(sg ** k + 1.0 + p["delta"] ** k for k in range(1, 5))


def v4_density_np(A, M, p):
    """copied from B3.e_parts (lines 'Me = M @ ETA' .. 'vd = sum(...)')."""
    t = v4_traces_np(M)
    cp = c4_of(p)
    return p["w"] * sum((t[k] - cp[k]) ** 2 for k in range(4))


def v4_grad_np(M, p):
    """copied from B3.grad (the V4 part), W1 -> p['w']."""
    Me = M @ ETA
    pows = [np.broadcast_to(np.eye(4), M.shape).copy()]
    for k in range(1, 4):
        pows.append(pows[-1] @ Me)
    t = [np.einsum("...kk->...", P @ Me) for P in pows]
    cp = c4_of(p)
    GV = np.zeros_like(M)
    for k in range(1, 5):
        coef = 2.0 * p["w"] * (t[k - 1] - cp[k - 1]) * k
        X = ETA @ pows[k - 1]
        GV += coef[..., None, None] * X.swapaxes(-1, -2)
    return B3.sym4(GV)


# ================= sympy layer (notebook conventions) =================
G_SYM = [[sp.Symbol(f"G{m}_{j}", real=True) for j in (1, 2, 3)]
         for m in range(4)]
T_SYM = [[sp.Symbol(f"T{m}_{j}", real=True) for j in (1, 2, 3)]
         for m in range(4)]
g_s, delta_s, w_s = sp.symbols("g delta w", positive=True)
s_s = sp.Symbol("s", real=True)


def gamma_mu(m, boost_style="real"):
    """the notebook Gamma_mu (m5_21_16_a_symbolic.gamma_mu, verbatim)."""
    t1, t2, t3 = T_SYM[m]
    if boost_style == "imag":
        row0 = [sp.I * t1, sp.I * t2, sp.I * t3]
        col0 = [-sp.I * t1, -sp.I * t2, -sp.I * t3]
    else:
        row0 = [t1, t2, t3]
        col0 = [t1, t2, t3]
    r1, r2, r3 = G_SYM[m]
    Gm = sp.zeros(4, 4)
    for j in range(3):
        Gm[0, 1 + j] = row0[j]
        Gm[1 + j, 0] = col0[j]
    Gm[1, 2], Gm[1, 3] = -r3, r2
    Gm[2, 1], Gm[2, 3] = r3, -r1
    Gm[3, 1], Gm[3, 2] = -r2, r1
    return Gm


def coms(A, B):
    return A * XI * B - B * XI * A


def notebook_jets(dvec=(g_s, 1, delta_s, 0), boost_style="real"):
    """M_mu = coms(Gamma_mu, d): the notebook jets at the point M = d."""
    d = sp.diag(*dvec)
    return [coms(gamma_mu(m, boost_style), d) for m in range(4)], d


def F_of_jets_sym(Mmu):
    return {(mu, nu): coms(Mmu[mu], Mmu[nu])
            for mu in range(4) for nu in range(4)}


def _e(i):
    """internal-index metric (mutable under the mutant)."""
    return XI_INT[i, i]


def _d(i):
    """derivative-index metric (never mutated)."""
    return XI[i, i]


def I1_sym(F, M, p):
    tot = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Fm = F[(mu, nu)]
            tot += _d(mu) * _d(nu) * sp.trace(XI_INT * Fm * XI_INT * Fm.T)
    return tot


def I1_frob_sym(F, M, p):
    tot = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Fm = F[(mu, nu)]
            tot += _d(mu) * _d(nu) * sp.trace(Fm * Fm.T)
    return tot


def I2_sym(F, M, p):
    return sum(F[(m, n)][a, b] * F[(a, b)][m, n]
               for m in range(4) for n in range(4)
               for a in range(4) for b in range(4))


def I3_sym(F, M, p):
    return sum(_d(m) * _e(b) * F[(m, n)][a, b] * F[(m, a)][n, b]
               for m in range(4) for n in range(4)
               for a in range(4) for b in range(4))


def I3_mixed_eta_sym(F, M, p):
    return sum(_d(m) * _d(n) * _d(a) * _e(b)
               * F[(m, n)][a, b] * F[(m, a)][n, b]
               for m in range(4) for n in range(4)
               for a in range(4) for b in range(4))


def R_sym(F):
    return sp.Matrix(4, 4, lambda n, a: sum(F[(m, n)][a, m]
                                            for m in range(4)))


def I4_sym(F, M, p):
    R = R_sym(F)
    return sum(_d(n) * _e(a) * R[n, a] ** 2
               for n in range(4) for a in range(4))


def I5_sym(F, M, p):
    R = R_sym(F)
    return sum(R[n, a] * R[a, n] for n in range(4) for a in range(4))


def I6_sym(F, M, p):
    R = R_sym(F)
    return sp.trace(R) ** 2


def V4_sym(F, M, p):
    Me = M * XI
    P = sp.eye(4)
    tot = 0
    for k in range(1, 5):
        P = P * Me
        Ck = (p["s"] * p["g"]) ** k + 1 + p["delta"] ** k
        tot += (sp.trace(P) - Ck) ** 2
    return p["w"] * tot


# ================= the registry =================
class Term:
    def __init__(self, name, definition, sympy_fn, K=None, np_fn=None,
                 grad_fn=None, kind="curvature"):
        self.name = name
        self.definition = definition
        self.hash = hashlib.sha256(definition.encode()).hexdigest()[:12]
        self.sympy_fn = sympy_fn
        self._K = K
        self._np_fn = np_fn
        self._grad_fn = grad_fn
        self.kind = kind

    # --- numpy: density per cell from jets A (4,...,4,4) and M (...,4,4)
    def density(self, A, M, p):
        if self.kind == "curvature":
            return density_from_K(F_of_A(A), self._K())
        return self._np_fn(A, M, p)

    def dA(self, A, M, p):
        """d density / dA_mu (curvature terms), shape like A."""
        F = F_of_A(A)
        return dA_from_W(dW_from_K(F, self._K()), A)

    def dM_local(self, A, M, p):
        """d density / dM at fixed jets (potential terms only)."""
        if self.kind == "curvature":
            return np.zeros_like(M)
        return self._grad_fn(M, p)

    def sympy(self, F, M, p):
        return self.sympy_fn(F, M, p)


def _pat(p1, p2, factor, **kw):
    return lambda: _K_from_pattern(p1, p2, factor, **kw)


REGISTRY = {}


def _reg(t):
    REGISTRY[t.name] = t
    return t


_reg(Term("I1",
          "I1 = sum_{mu<nu} eta^mu eta^nu <F_mu nu, F_mu nu>_eta, "
          "<F,G>_eta = tr(eta F eta G^T), F_mu nu = A_mu eta A_nu - "
          "A_nu eta A_mu; = (1/2) F_{mu nu a b} F^{mu nu a b}",
          I1_sym, K=_pat("mnab", "mnab", 0.5)))
_reg(Term("I1_frob",
          "I1_frob = sum_{mu<nu} eta^mu eta^nu tr(F_mu nu F_mu nu^T) "
          "(Frobenius internal metric; M5.21.16 variant A; control)",
          I1_frob_sym, K=_pat("mnab", "mnab", 0.5, int_metric="delta")))
_reg(Term("I2",
          "I2 = F_{mu nu a b} F^{a b mu nu} (derivative pair of one F "
          "contracted with the internal pair of the other; mixed pairs "
          "use delta)",
          I2_sym, K=_pat("mnab", "abmn", 1.0)))
_reg(Term("I3",
          "I3 = F_{mu nu a b} F^{mu a nu b} (eta on the mu-mu and b-b "
          "pairs, delta on the mixed nu-a pairs)",
          I3_sym, K=_pat("mnab", "manb", 1.0)))
_reg(Term("I3_mixed_eta",
          "I3_mixed_eta = sum eta_mu eta_nu eta_a eta_b F[mu nu a b] "
          "F[mu a nu b] (eta on the mixed pairs; NON-covariant control)",
          I3_mixed_eta_sym, K=_pat("mnab", "manb", 1.0,
                                   mixed_metric="eta")))
_reg(Term("I4",
          "I4 = R_{nu a} R^{nu a} = sum eta_nu eta_a R[nu,a]^2, "
          "R[nu,a] = sum_mu F[mu nu a mu]",
          I4_sym, K=_pat("mnam", "pnap", 1.0)))
_reg(Term("I5",
          "I5 = R_{nu a} R^{a nu} = sum_{nu a} R[nu,a] R[a,nu], "
          "R[nu,a] = sum_mu F[mu nu a mu]",
          I5_sym, K=_pat("mnam", "panp", 1.0)))
_reg(Term("I6",
          "I6 = R^2, R = sum_nu R[nu,nu] = sum_{mu nu} F[mu nu nu mu]",
          I6_sym, K=_pat("mnnm", "pqqp", 1.0)))
_reg(Term("V4",
          "V4 = w sum_{p=1..4} (tr((M eta)^p) - C_p)^2, C_p = (s g)^p + 1 "
          "+ delta^p, M_vac = diag(-s g, 1, delta, 0)",
          V4_sym, np_fn=v4_density_np, grad_fn=v4_grad_np,
          kind="potential"))

# the certified Lagrangian: L = sum_k c_k I_k, E = sum_k c_k (C w^2 - A)
CERTIFIED_COEFFS = {"I1": -4.0, "V4": -1.0}


def default_params(**kw):
    p = {"s": 1.0, "g": B3.G_T, "delta": B3.DELTA0, "w": W1}
    p.update(kw)
    return p


# ================= lattice energies =================
def lattice_jets(M, cfg, a0=None, omega=0.0):
    """per stencil branch: (A (4,n,n,n,4,4), weight); A[0] = omega a0."""
    out = []
    for br, wt in B3.branches(cfg["stencil"]):
        A = np.zeros((4,) + M.shape, dtype=M.dtype)
        for ax in range(3):
            A[1 + ax] = B3.d1(M, ax, cfg["h"], br)
        if a0 is not None and omega != 0.0:
            A[0] = omega * a0
        out.append((A, wt))
    return out


def term_lagrangian(term, M, cfg, p, a0=None, omega=0.0):
    """h^3 sum over cells of the term density (the Lagrangian read)."""
    tot = 0.0
    for A, wt in lattice_jets(M, cfg, a0, omega):
        tot = tot + wt * np.sum(term.density(A, M, p))
    return cfg["h"] ** 3 * tot


def term_hamiltonian(term, M, cfg, p, a0=None, omega=0.0):
    """H = C omega^2 - A for I(omega) = A + B omega + C omega^2 (exact:
    every registered term is quadratic in omega)."""
    if a0 is None or omega == 0.0:
        return -term_lagrangian(term, M, cfg, p)
    lp = term_lagrangian(term, M, cfg, p, a0, omega)
    lm = term_lagrangian(term, M, cfg, p, a0, -omega)
    l0 = term_lagrangian(term, M, cfg, p)
    return 0.5 * (lp + lm) - 2.0 * l0


def term_energy(term, M, cfg, p, coef, a0=None, omega=0.0):
    return coef * term_hamiltonian(term, M, cfg, p, a0, omega)


def omega_decompose(term, M, cfg, p, a0):
    """(A, B, C) of the Lagrangian read I(omega) = A + B omega + C omega^2."""
    l0 = term_lagrangian(term, M, cfg, p)
    lp = term_lagrangian(term, M, cfg, p, a0, 1.0)
    lm = term_lagrangian(term, M, cfg, p, a0, -1.0)
    return l0, 0.5 * (lp - lm), 0.5 * (lp + lm) - l0


def term_grad_lagrangian(term, M, cfg, p, a0=None, omega=0.0):
    """d/dM of term_lagrangian (a0 FROZEN, the velocity-field read),
    symmetrized; exact adjoints of the certified stencil."""
    G = np.zeros_like(M)
    h3 = cfg["h"] ** 3
    for br, wt in B3.branches(cfg["stencil"]):
        A = np.zeros((4,) + M.shape, dtype=M.dtype)
        for ax in range(3):
            A[1 + ax] = B3.d1(M, ax, cfg["h"], br)
        if a0 is not None and omega != 0.0:
            A[0] = omega * a0
        if term.kind == "curvature":
            dA = term.dA(A, M, p)
            for ax in range(3):
                G += wt * B3.d1_adj(dA[1 + ax], ax, cfg["h"], br)
        else:
            G += wt * term.dM_local(A, M, p)
    return h3 * B3.sym4(G)


def term_grad_hamiltonian(term, M, cfg, p, a0=None, omega=0.0):
    if a0 is None or omega == 0.0:
        return -term_grad_lagrangian(term, M, cfg, p)
    gp = term_grad_lagrangian(term, M, cfg, p, a0, omega)
    gm = term_grad_lagrangian(term, M, cfg, p, a0, -omega)
    g0 = term_grad_lagrangian(term, M, cfg, p)
    return 0.5 * (gp + gm) - 2.0 * g0


def certified_energy(M, cfg, p, a0=None, omega=0.0):
    """E = 4 (U + omega^2 T) + V4 through the registry."""
    return sum(term_energy(REGISTRY[k], M, cfg, p, c, a0, omega)
               for k, c in CERTIFIED_COEFFS.items())


# ================= selftests =================
def _rel(a, b):
    return float(abs(a - b) / max(abs(b), 1e-300))


def _lorentz(rng, kind, scale=0.3):
    Gm = np.zeros((4, 4))
    if kind == "boost":
        v = scale * rng.normal(size=3)
        Gm[0, 1:] = v
        Gm[1:, 0] = v
    else:
        v = scale * rng.normal(size=3)
        Gm[1, 2], Gm[2, 1] = -v[2], v[2]
        Gm[1, 3], Gm[3, 1] = v[1], -v[1]
        Gm[2, 3], Gm[3, 2] = -v[0], v[0]
    L = expm(Gm)
    assert np.allclose(L.T @ ETA @ L, ETA, atol=1e-12)
    return L


def _random_jets(rng, npts, p, amp=0.5):
    cfg = B3.base_cfg(s=p["s"], g=p["g"], delta=p["delta"])
    M = B3.vac4(cfg)[None] + amp * B3.sym4(rng.normal(size=(npts, 4, 4)))
    A = B3.sym4(rng.normal(size=(4, npts, 4, 4)))
    return A, M


def st_certified(res, lines):
    """(a) I1 / V4 vs the certified stack on stored + regenerated fields."""
    t0 = time.time()
    # stored 3x3 end state (M5.21.11 ladder, n = 32, L = 48, sym stencil)
    # embedded at s = +1, g = 8: I1 must equal the stored E_u (the 3D
    # plain-commutator read; B3 gates g3) and the live certified e_parts
    cfg = B3.base_cfg(s=1.0)
    p = default_params(s=1.0)
    M4 = B3.embed34(load_stored3(), cfg)
    eu, ev = B3.e_parts(M4, cfg)
    e1 = term_energy(REGISTRY["I1"], M4, cfg, p, -4.0)
    e4 = term_energy(REGISTRY["V4"], M4, cfg, p, -1.0)
    with open(os.path.join(DATA, STORED3_ROW)) as f:
        e_u_stored = json.load(f)["E_u"]
    r = {"stored_field": STORED3_NPZ,
         "E_u_registry": float(e1), "E_u_certified_live": float(eu),
         "rel_u": _rel(e1, eu),
         "E_V_registry": float(e4), "E_V_certified_live": float(ev),
         "rel_V": _rel(e4, ev),
         "E_u_stored": e_u_stored,
         "rel_u_vs_stored": _rel(e1, e_u_stored)}
    # the M5.21.16 IDENT field regenerated (rng 21160, first draws)
    from scipy.ndimage import gaussian_filter
    rng = np.random.default_rng(21160)
    cfg16 = B3.base_cfg(n=16, L=24.0)
    M3 = np.stack([[gaussian_filter(rng.normal(size=(16,) * 3), 2.0)
                    for _ in range(3)] for _ in range(3)], axis=-1)
    M3 = M3.reshape(16, 16, 16, 3, 3)
    M3 = 0.5 * (M3 + M3.swapaxes(-1, -2))
    Mi = B3.embed34(M3, cfg16)
    with open(os.path.join(DATA, "m5_21_16_field.json")) as f:
        f16 = json.load(f)
    u_st = f16["IDENT"]["u_eta"]
    ui = term_energy(REGISTRY["I1"], Mi, cfg16, default_params(), -4.0)
    uf = term_energy(REGISTRY["I1_frob"], Mi, cfg16, default_params(),
                     -4.0)
    r["ident_u_registry"] = float(ui)
    r["ident_u_stored"] = u_st
    r["rel_ident_eta"] = _rel(ui, u_st)
    r["rel_ident_frob"] = _rel(uf, f16["IDENT"]["u_flip"])
    # CHAN: kin per channel on the analytic family (g=32, s=-1, n=32)
    B16 = _load("m5_21_16_b_field", "m5_21_16_b_field.py")
    cfgc = B3.base_cfg(s=-1.0, g=32.0, n=32, L=48.0)
    pc = default_params(s=-1.0, g=32.0)
    Mc = B16.lattice_family_M(cfgc, 32.0)
    a0s = B3.gen_catalog(cfgc, Mc)
    chan = {}
    worst = 0.0
    for nm in ("boost_z", "boost_x", "clock_local"):
        k_cert = float(B3.kin_of(Mc, a0s[nm], cfgc))
        e_static = term_energy(REGISTRY["I1"], Mc, cfgc, pc, -4.0)
        k_reg = term_energy(REGISTRY["I1"], Mc, cfgc, pc, -4.0,
                            a0s[nm], 1.0) - e_static
        k_frob = term_energy(REGISTRY["I1_frob"], Mc, cfgc, pc, -4.0,
                             a0s[nm], 1.0) - term_energy(
            REGISTRY["I1_frob"], Mc, cfgc, pc, -4.0)
        st = f16["CHAN"]["rows"][nm]
        chan[nm] = {"kin_registry": float(k_reg), "kin_certified": k_cert,
                    "kin_stored": st["kin_eta"],
                    "kin_frob_registry": float(k_frob),
                    "kin_frob_stored": st["kin_flip"],
                    "rel_live": _rel(k_reg, k_cert),
                    "rel_stored": _rel(k_reg, st["kin_eta"]),
                    "rel_frob_stored": _rel(k_frob, st["kin_flip"])}
        worst = max(worst, chan[nm]["rel_live"], chan[nm]["rel_stored"],
                    chan[nm]["rel_frob_stored"])
    r["chan"] = chan
    r["chan_worst_rel"] = worst
    # vacuum null, both branches
    vac = {}
    for s in (1.0, -1.0):
        c2 = B3.base_cfg(s=s, n=8, L=12.0)
        Mv = np.zeros((8, 8, 8, 4, 4)) + B3.vac4(c2)
        vac[f"s{int(s):+d}"] = float(abs(certified_energy(
            Mv, c2, default_params(s=s))))
    r["vacuum_abs"] = vac
    r["runtime_s"] = round(time.time() - t0, 1)
    res["certified"] = r
    lines.append(("I1 E_u vs certified e_parts (stored field, n=32)", r["rel_u"],
                  1e-12, r["rel_u"] <= 1e-12))
    lines.append(("V4 vs certified e_parts (stored field)", r["rel_V"], 1e-12,
                  r["rel_V"] <= 1e-12))
    lines.append(("I1 vs stored E_u (m5_21_11_row_t11lad_A_n32_d0.3)",
                  r["rel_u_vs_stored"], 1e-12,
                  r["rel_u_vs_stored"] <= 1e-12))
    lines.append(("I1 vs stored M5.21.16 IDENT u_eta (regen)",
                  r["rel_ident_eta"], 1e-12, r["rel_ident_eta"] <= 1e-12))
    lines.append(("I1 kin vs certified/stored CHAN (eta+frob, 3 ch.)",
                  worst, 1e-12, worst <= 1e-12))
    lines.append(("vacuum energy == 0 (both s)", max(vac.values()),
                  1e-16, max(vac.values()) <= 1e-16))


def st_sympy_vs_numpy(res, lines):
    """(b) every term: exact sympy at rational jets vs numpy, and the
    notebook N1 coefficient through the registry's I1."""
    t0 = time.time()
    rng = np.random.default_rng(3201)
    p = default_params(s=-1.0, g=32.0)
    out = {}
    worst = 0.0
    for trial in range(2):
        # random rational symmetric jets and point value
        def rsym():
            X = sp.zeros(4, 4)
            for i in range(4):
                for j in range(i, 4):
                    X[i, j] = X[j, i] = sp.Rational(
                        int(rng.integers(-9, 10)), int(rng.integers(1, 7)))
            return X
        Mmu = [rsym() for _ in range(4)]
        Mpt = sp.diag(-p["s"] * p["g"], 1, p["delta"], 0)
        Mpt = sp.Matrix(Mpt) + rsym() / 3
        F = F_of_jets_sym(Mmu)
        A_np = np.array([[[float(Mmu[m][i, j]) for j in range(4)]
                          for i in range(4)] for m in range(4)])[:, None]
        M_np = np.array([[float(Mpt[i, j]) for j in range(4)]
                         for i in range(4)])[None]
        ps = {"s": sp.Integer(int(p["s"])), "g": sp.Integer(32),
              "delta": sp.Rational(3, 10), "w": sp.Integer(1)}
        pn = dict(p, w=1.0)
        for nm, T in REGISTRY.items():
            vs = float(sp.nsimplify(T.sympy(F, Mpt, ps)))
            vn = float(T.density(A_np, M_np, pn)[0])
            rel = _rel(vn, vs)
            out.setdefault(nm, []).append(
                {"sympy": vs, "numpy": vn, "rel": rel})
            worst = max(worst, rel)
    res["sympy_vs_numpy"] = out
    lines.append(("sympy vs numpy, every term, 2 rational jets",
                  worst, 1e-12, worst <= 1e-12))
    # N1: the omega^2 leading coefficient of the I1 Hamiltonian read
    Mmu, d = notebook_jets()
    F = F_of_jets_sym(Mmu)
    H = sp.expand(-(sp.expand(I1_sym(F, d, None))))   # Lagrangian read
    om = G_SYM[0][0]
    # Hamiltonian = C w^2 - A
    Apart = H.coeff(om, 0)
    Cpart = H.coeff(om, 2)
    Bpart = H.coeff(om, 1)
    Hh = sp.expand(Cpart * om ** 2 - Apart)
    c2 = sp.expand(Hh.coeff(om, 2).subs(g_s, 1 / delta_s))
    lead = sp.expand(sp.series(c2, delta_s, 0, 1).removeO())
    SIX = sum(T_SYM[m][j] ** 2 for m in (1, 2, 3) for j in (1, 2))
    n1 = bool(sp.simplify(lead - (-2 * SIX)) == 0)
    Hf = sp.expand(-(I1_frob_sym(F, d, None)))
    Hhf = sp.expand(Hf.coeff(om, 2) * om ** 2 - Hf.coeff(om, 0))
    c2f = sp.expand(Hhf.coeff(om, 2).subs(g_s, 1 / delta_s))
    leadf = sp.expand(sp.series(c2f, delta_s, 0, 1).removeO())
    n4 = bool(sp.simplify(leadf - (2 * SIX)) == 0)
    with open(os.path.join(DATA, "m5_21_16_symbolic.json")) as f:
        s16 = json.load(f)
    res["notebook_N1"] = {"lead": str(lead), "stored": s16["N1_baseline_lead"],
                          "matches_minus2six": n1,
                          "G0_1_linear_term_via_other_Gamma0_entries": str(Bpart) != "0",
                          "frob_lead": str(leadf), "matches_plus2six": n4,
                          "runtime_s": round(time.time() - t0, 1)}
    lines.append(("notebook N1: I1 omega^2 lead == -2*SIX (sympy)",
                  0.0 if n1 else 1.0, 0.5, n1))
    lines.append(("notebook N4: I1_frob omega^2 lead == +2*SIX (sympy)",
                  0.0 if n4 else 1.0, 0.5, n4))


def st_gradient(res, lines):
    """(c) complex-step gradient of every numpy term (Lagrangian read,
    static and with a frozen clock a0, omega = 0.7)."""
    t0 = time.time()
    rng = np.random.default_rng(3202)
    cfg = B3.base_cfg(n=10, L=15.0, s=-1.0, g=32.0)
    p = default_params(s=-1.0, g=32.0)
    M = B3.vac4(cfg)[None, None, None] + 0.5 * B3.sym4(
        rng.normal(size=(10, 10, 10, 4, 4)))
    a0 = B3.sym4(rng.normal(size=M.shape))
    a0 /= np.sqrt(np.sum(a0 * a0))
    out = {}
    worst = 0.0
    for nm, T in REGISTRY.items():
        errs = []
        for (aa, om) in ((None, 0.0), (a0, 0.7)):
            G = term_grad_hamiltonian(T, M, cfg, p, aa, om)
            for _ in range(2):
                V = B3.sym4(rng.normal(size=M.shape))
                de_an = float(np.sum(G * V))
                t = 1e-30
                de = term_hamiltonian(T, M + 1j * t * V, cfg, p, aa,
                                      om).imag / t
                errs.append(_rel(de_an, de))
        out[nm] = float(max(errs))
        worst = max(worst, out[nm])
    res["gradient_complex_step"] = out
    res["gradient_runtime_s"] = round(time.time() - t0, 1)
    lines.append(("complex-step gradient, every term (static + clock)",
                  worst, 5e-9, worst <= 5e-9))


def st_invariance(res, lines):
    """(d) density-level invariance under random SO(1,3) boosts and SO(3)
    rotations acting on M by conjugation AND on the derivative index
    (A'_mu = (Lambda^{-T})_mu^nu Lambda A_nu Lambda^T); plus the lattice
    internal-conjugation test of I1 (the certified g1_so13 pattern)."""
    t0 = time.time()
    rng = np.random.default_rng(3203)
    p = default_params(s=-1.0, g=32.0)
    A, M = _random_jets(rng, 64, p)
    out = {}
    controls = {"I1_frob", "I3_mixed_eta"}
    worst_cov = 0.0
    min_ctrl = np.inf
    for kind in ("boost", "rotation"):
        for k in range(3):
            L = _lorentz(rng, kind)
            Linv_T = np.linalg.inv(L).T
            Mp = np.einsum("ab,nbc,dc->nad", L, M, L)
            Ap = np.einsum("mn,nxab->mxab", Linv_T,
                           np.einsum("ab,nxbc,dc->nxad", L, A, L))
            for nm, T in REGISTRY.items():
                d0 = T.density(A, M, p)
                d1 = T.density(Ap, Mp, p)
                drift = float(np.max(np.abs(d1 - d0))
                              / max(np.max(np.abs(d0)), 1e-300))
                key = f"{kind}_{k}"
                out.setdefault(nm, {})[key] = drift
                if nm in controls:
                    if kind == "boost":
                        min_ctrl = min(min_ctrl, drift)
                else:
                    worst_cov = max(worst_cov, drift)
    # lattice-level internal conjugation of I1 (derivative indices only
    # contracted among themselves for I1, so conjugation alone suffices)
    cfgB = B3.base_cfg(n=14, L=21.0)
    M3 = load_stored3()
    sub = M3[::2, ::2, ::2][:14, :14, :14]
    MB = B3.embed34(sub, cfgB)
    pB = default_params()
    E0 = term_energy(REGISTRY["I1"], MB, cfgB, pB, -4.0)
    Gm = np.zeros((4, 4))
    Gm[0, 1] = Gm[1, 0] = 0.11
    Gm[2, 3], Gm[3, 2] = -0.23, 0.23
    L = expm(Gm)
    ML = np.einsum("ab,...bc,dc->...ad", L, MB, L)
    lat = _rel(term_energy(REGISTRY["I1"], ML, cfgB, pB, -4.0), E0)
    with open(os.path.join(DATA, "m5_21_3_gates.json")) as f:
        g13 = json.load(f)["g1_so13"]
    res["invariance"] = {"per_term": out, "worst_covariant": worst_cov,
                         "min_control_boost_drift": float(min_ctrl),
                         "I1_lattice_conjugation_rel": lat,
                         "I1_lattice_conjugation_stored_g1_so13": g13,
                         "runtime_s": round(time.time() - t0, 1)}
    lines.append(("invariance: covariant terms, SO(1,3)+SO(3), density",
                  worst_cov, 1e-10, worst_cov <= 1e-10))
    lines.append(("invariance control: I1_frob / I3_mixed_eta boost drift",
                  float(min_ctrl), 1e-3, min_ctrl > 1e-3))
    lines.append(("I1 lattice internal-conjugation (g1_so13 pattern)",
                  lat, 1e-9, lat <= 1e-9))


def st_reduction(res, lines):
    """(f) the 3x3 reduction: static uniform-time-row field (block
    diagonal, vacuum time entry, no time dependence): per term vanishes /
    = c * I1 / differs; plus the plain-commutator 3D read of I1."""
    t0 = time.time()
    rng = np.random.default_rng(3204)
    from scipy.ndimage import gaussian_filter
    cfg = B3.base_cfg(n=16, L=24.0, s=-1.0, g=32.0)
    p = default_params(s=-1.0, g=32.0)
    vals = {nm: [] for nm in REGISTRY}
    e3s = []
    for trial in range(3):
        M3 = np.stack([[gaussian_filter(rng.normal(size=(16,) * 3), 2.0)
                        for _ in range(3)] for _ in range(3)], axis=-1)
        M3 = M3.reshape(16, 16, 16, 3, 3)
        M3 = 0.5 * (M3 + M3.swapaxes(-1, -2))
        M3 = M3 + np.diag([1.0, p["delta"], 0.0])
        M4 = B3.embed34(M3, cfg)
        for nm, T in REGISTRY.items():
            vals[nm].append(float(term_lagrangian(T, M4, cfg, p)))
        # plain 3D commutator read (B3.gates g3 pattern)
        h = cfg["h"]
        e3 = 0.0
        for br, wt in B3.branches("sym"):
            Ad = [B3.d1(M3, ax, h, br) for ax in range(3)]
            for i in range(3):
                for j in range(i + 1, 3):
                    C = Ad[i] @ Ad[j] - Ad[j] @ Ad[i]
                    e3 += wt * np.sum(np.einsum("...kl,...kl->...", C, C))
        e3s.append(float(e3 * h ** 3))
    col = {}
    i1 = np.array(vals["I1"])
    for nm in REGISTRY:
        v = np.array(vals[nm])
        if np.max(np.abs(v)) <= 1e-12 * np.max(np.abs(i1)):
            col[nm] = {"status": "vanishes identically", "ratio_to_I1": 0.0}
        else:
            ratios = v / i1
            if np.max(np.abs(ratios - ratios[0])) <= 1e-10 * abs(ratios[0]):
                col[nm] = {"status": f"= {ratios[0]:.12g} * I1",
                           "ratio_to_I1": float(ratios[0])}
            else:
                col[nm] = {"status": "differs from I1",
                           "ratios_to_I1": ratios.tolist()}
        col[nm]["values"] = v.tolist()
    g3 = max(_rel(a, b) for a, b in zip(vals["I1"], e3s))
    res["reduction_3x3"] = {"column": col, "I1_vs_plain_3D_rel": g3,
                            "runtime_s": round(time.time() - t0, 1)}
    lines.append(("3x3 reduction: I1 == plain-commutator 3D read",
                  g3, 1e-12, g3 <= 1e-12))
    lines.append(("3x3 reduction: I1_frob == I1 (identity)",
                  abs(col["I1_frob"].get("ratio_to_I1", 0.0) - 1.0),
                  1e-12, abs(col["I1_frob"].get("ratio_to_I1", 0.0) - 1.0)
                  <= 1e-12))


def st_R_pairings(res, lines):
    """the R_ac pairings: which vanish, which are +-R."""
    rng = np.random.default_rng(3205)
    p = default_params(s=-1.0, g=32.0)
    A, M = _random_jets(rng, 32, p)
    F = F_of_A(A)
    rd = R_readings_np(F)
    R1 = rd["R_reading_1 (slot0,slot3)"]
    scale = float(np.max(np.abs(R1)))
    out = {}
    for k, v in rd.items():
        mx = float(np.max(np.abs(v)))
        if mx <= 1e-13 * scale:
            out[k] = "vanishes identically"
        elif np.max(np.abs(v - R1)) <= 1e-13 * scale:
            out[k] = "= +R_reading_1"
        elif np.max(np.abs(v + R1)) <= 1e-13 * scale:
            out[k] = "= -R_reading_1"
        else:
            out[k] = "independent"
    sym = float(np.max(np.abs(R1 - R1.swapaxes(-1, -2))) / scale)
    out["R_reading_1 symmetric?"] = f"no (max asym {sym:.3g})"
    res["R_pairings"] = out
    ok = (out["trace_deriv_pair eta^{mu nu} F_{mu nu ab}"]
          == "vanishes identically"
          and out["trace_int_pair eta_{ab} F_{mu nu ab}"]
          == "vanishes identically"
          and all(out[k] in ("= +R_reading_1", "= -R_reading_1")
                  for k in list(rd)[1:4]))
    lines.append(("R pairings: within-type traces vanish, 4 mixed = +-R",
                  0.0 if ok else 1.0, 0.5, ok))


def st_mutation(res, lines):
    """(e) run the selftests under the eta_time_row mutant in a subprocess
    and count reddened lines."""
    t0 = time.time()
    cmd = [PY, os.path.abspath(__file__), "--selftest", "--mutant",
           "eta_time_row", "--no-mutation-stage", "--no-write"]
    pr = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
    reds = [ln for ln in pr.stdout.splitlines() if ln.startswith("FAIL")]
    res["mutation"] = {"mutant": "eta_time_row", "n_fail": len(reds),
                       "failed_lines": reds, "returncode": pr.returncode,
                       "runtime_s": round(time.time() - t0, 1)}
    lines.append(("mutation: eta_time_row mutant reddens >= 1 line",
                  float(len(reds)), 1.0, len(reds) >= 1))


def selftest(mutation_stage=True, write=True):
    t0 = time.time()
    res = {"mutant": MUTANT, "python": PY,
           "terms": {nm: {"definition": T.definition, "hash": T.hash,
                          "kind": T.kind} for nm, T in REGISTRY.items()},
           "certified_coeffs": CERTIFIED_COEFFS}
    lines = []
    for st in (st_certified, st_sympy_vs_numpy, st_gradient,
               st_invariance, st_reduction, st_R_pairings):
        st(res, lines)
    if mutation_stage:
        st_mutation(res, lines)
    res["selftest_lines"] = [{"name": n, "value": v, "threshold": t,
                              "pass": bool(ok)} for n, v, t, ok in lines]
    res["all_pass"] = all(ok for *_, ok in lines)
    res["runtime_s"] = round(time.time() - t0, 1)
    for n, v, t, ok in lines:
        print(f"{'PASS' if ok else 'FAIL'}  {n:58s} {v:.3e}  (thr {t:g})")
    print(f"ALL_PASS={res['all_pass']}  mutant={MUTANT}  "
          f"runtime {res['runtime_s']} s")
    if write:
        os.makedirs(DATA, exist_ok=True)
        with open(os.path.join(DATA, "m5_32_r0_module_selftest.json"),
                  "w") as f:
            json.dump(res, f, indent=1)
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--mutant", default=None)
    ap.add_argument("--no-mutation-stage", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    set_mutant(a.mutant)
    if a.list or not a.selftest:
        for nm, T in REGISTRY.items():
            print(f"{nm:14s} {T.hash}  {T.definition}")
    if a.selftest:
        r = selftest(mutation_stage=not a.no_mutation_stage,
                     write=not a.no_write)
        sys.exit(0 if r["all_pass"] else 1)
