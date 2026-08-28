"""M5.32 R1.a: the EXTENDED covariant term catalog (companion registry).

Companion of m5_32_lagrangian.py (imported, never modified): the same
registry shape (name, definition string, sympy callable, numpy callable,
hash), the same conventions (raw entries of M and F are CONTRAVARIANT
internal tensors, d_mu is covariant, eta = diag(-1, 1, 1, 1)), plus the
terms R1 adds on top of I1..I6:

EQUATIONS FIRST
---------------
Covariantized components (the R0 contraction rule, restated as one rule):
    F_{mu nu a b} := eta_{a a'} eta_{b b'} F[mu nu a' b']   (all indices low)
    eps_{i j k l} := the Levi-Civita SYMBOL (all indices low)
    every contracted pair is raised with eta^{..}
which reproduces the R0 rule (dd: eta, ii: eta, di: delta) and extends it
to epsilon slots (eps-derivative pair: eta, eps-internal pair: delta).
Under M -> L M L^T, A_mu -> (L^-T)_mu^nu L A_nu L^T with L in SO(1,3)
(det L = +1) every such contraction is invariant; with det L = -1 (a
parity or time reflection) the one-epsilon contractions flip sign
(pseudoscalars) while the eta-only ones do not.

(1) The parity-odd quadratic invariants: every full contraction of
    F (x) F with exactly one epsilon: the 4 epsilon slots take 4 of the 8
    F slots (C(8,4) = 70 choices; the slot order only fixes the sign) and
    the 4 remaining F slots are paired among themselves (3 pairings):
    210 patterns, reduced by the rank on generic jets to the basis
    E1..Er registered below (r measured, the external claim is FOUR).

(2) The C2 insertions (field-dependent internal metrics):
    u(M): the timelike unit eigenvector of N = M eta, N u = lambda u,
          u^a CONTRAVARIANT (N' = L N L^-1 so u' = L u), normalized
          u^T eta u = -1 (its sign is irrelevant: only u u^T enters).
          Why M eta and not eta M: M^{ab} is a (2,0) tensor, (M eta)^a_b
          is the (1,1) endomorphism whose eigenvectors are vectors; eta M
          has the covariant eigenvectors eta u (the same objects, lowered).
    h^{ab}(M) = eta^{ab} + 2 u^a u^b         (contravariant, the flip metric)
    h_{ab}    = eta_{ac} h^{cd} eta_{db} = eta_{ab} + 2 (eta u)_a (eta u)_b
    I1_h  = sum_{mu<nu} eta^mu eta^nu tr(h_cov F h_cov F^T)
          = sum_{mu<nu} eta^mu eta^nu h_{ac} h_{bd} F^{ab} F^{cd}
          (in the vacuum eigenframe u = e0, h_cov = 1: I1_h = I1_frob)
    J1    = sum_{mu<nu} eta^mu eta^nu tr(F eta M eta F eta M eta)
            [tr((F M)^2)-type; contravariant F, M joined by eta_ab]
    J2    = sum_{mu<nu} eta^mu eta^nu tr(F eta F eta M eta M eta)
            [tr(F^2 M^2)-type]
    Pgrad = sum_mu eta^{mu mu} q(d_mu P_t, d_mu P_t),  P_t = u u^T eta
            ((1,1) tensor; with u^T eta u = -1, P_t^2 = -P_t, so the
            projector proper is -P_t: the sign drops out of q(dP, dP)),
            q(X, Y) = X^a_b Y^c_d eta_{ac} eta^{bd} = sum_ab eta_a eta_b X[a,b] Y[a,b]
            d_mu u = sum_{k != 0} sigma_k (u_k^T eta A_mu eta u_0)
                     / (lambda_0 - lambda_k) u_k
            (first-order eigenvector perturbation of N = M eta, u_k the
            full eigenbasis with u_k^T eta u_j = sigma_j delta_kj,
            sigma = +-1; the u_0 component vanishes by the normalization)
            The search knob -kappa is the coefficient, not the term.

SELFTESTS (python3 m5_32_terms_ext.py --selftest [--mutant eta_time_row]):
covariance drift <= 1e-10 per new term (SO(1,3) boosts + SO(3) rotations),
the parity SIGN FLIP of every epsilon term (and the non-flip of the even
ones), sympy vs numpy on rational jets at a rationally boosted vacuum
point, I1_h == I1_frob in the vacuum eigenframe, the d_mu u perturbation
formula vs a central difference, the epsilon rank, and the mutant
reddening (the mutant flips the internal time-row sign of the contraction
metric of every new term).

Out: ../data/m5_32_r1_ext_selftest.json
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


L0 = _load("m5_32_lagrangian", "m5_32_lagrangian.py")   # the base registry
B3 = L0.B3

ETA = L0.ETA                                   # bracket metric (fixed)
ETA_X = np.diag([-1.0, 1.0, 1.0, 1.0])         # ext contraction metric (mutable)
XI = L0.XI
XI_X = sp.diag(-1, 1, 1, 1)
MUTANT = None


def set_mutant(name):
    global ETA_X, XI_X, MUTANT
    if name == "eta_time_row":
        ETA_X = np.diag([1.0, 1.0, 1.0, 1.0])
        XI_X = sp.diag(1, 1, 1, 1)
        MUTANT = name
    elif name:
        raise ValueError(name)


# ================= epsilon contraction matrices =================
def _eps4():
    E = np.zeros((4,) * 4)
    for perm in itertools.permutations(range(4)):
        sgn = 1
        for i in range(4):
            for j in range(i + 1, 4):
                if perm[i] > perm[j]:
                    sgn = -sgn
        E[perm] = sgn
    return E


EPS4 = _eps4()


def K_eps(p1, p2, pe, factor=1.0):
    """256x256 contraction matrix for F[p1] F[p2] eps[pe]: every letter
    appears exactly twice; slots 0,1 of p1/p2 are derivative, 2,3 internal;
    F-F pairs follow the R0 rule (dd eta, ii eta_int, di delta); eps-F
    pairs: eps-derivative eta, eps-internal delta."""
    occ = {}
    for k, l in enumerate(p1):
        occ.setdefault(l, []).append("d" if k < 2 else "i")
    for k, l in enumerate(p2):
        occ.setdefault(l, []).append("d" if k < 2 else "i")
    for l in pe:
        occ.setdefault(l, []).append("e")
    metric = {}
    for l, kinds in occ.items():
        assert len(kinds) == 2, (l, kinds)
        ks = tuple(sorted(kinds))
        assert ks != ("e", "e"), "eps-eps contraction vanishes"
        if ks == ("d", "d") or ks == ("d", "e"):
            metric[l] = np.diag(ETA)
        elif ks == ("i", "i"):
            metric[l] = np.diag(ETA_X)
        else:                                  # di, ei
            metric[l] = np.ones(4)
    letters = list(occ)
    K = np.zeros((4,) * 8)
    for vals in itertools.product(range(4), repeat=len(letters)):
        env = dict(zip(letters, vals))
        e = EPS4[tuple(env[l] for l in pe)]
        if e == 0.0:
            continue
        w = factor * e
        for l, v in env.items():
            w *= metric[l][v]
        i = tuple(env[l] for l in p1)
        j = tuple(env[l] for l in p2)
        K[i + j] += w
    return K.reshape(256, 256)


def eps_patterns():
    """the 210 one-epsilon patterns: (p1, p2, pe, label)."""
    slots = [("F1", "d", 0), ("F1", "d", 1), ("F1", "i", 2), ("F1", "i", 3),
             ("F2", "d", 0), ("F2", "d", 1), ("F2", "i", 2), ("F2", "i", 3)]
    out = []
    for S in itertools.combinations(range(8), 4):
        rest = [k for k in range(8) if k not in S]
        pairings = [((rest[0], rest[1]), (rest[2], rest[3])),
                    ((rest[0], rest[2]), (rest[1], rest[3])),
                    ((rest[0], rest[3]), (rest[1], rest[2]))]
        for pr in pairings:
            let = [None] * 8
            for k, s in enumerate(S):
                let[s] = "pqrs"[k]
            for k, (a, b) in enumerate(pr):
                let[a] = let[b] = "xy"[k]
            p1 = "".join(let[:4])
            p2 = "".join(let[4:])
            lab = ("eps[" + ",".join(f"{slots[s][0]}.{slots[s][1]}{slots[s][2]}"
                                     for s in S) + "] "
                   + " ".join(f"({slots[a][0]}.{slots[a][1]}{slots[a][2]}-"
                              f"{slots[b][0]}.{slots[b][1]}{slots[b][2]})"
                              for a, b in pr))
            out.append((p1, p2, "pqrs", lab))
    return out


# ================= sympy helpers =================
def sym_from_K(F, K):
    """f^T K f on a symbolic F dict (the (mu, nu) 4x4 blocks)."""
    f = []
    for mu in range(4):
        for nu in range(4):
            Fm = F[(mu, nu)]
            for a in range(4):
                for b in range(4):
                    f.append(Fm[a, b])
    Ks = sp.Matrix(K)
    tot = 0
    rows, cols = np.nonzero(K)
    for i, j in zip(rows, cols):
        tot += sp.nsimplify(K[i, j]) * f[i] * f[j]
    return sp.expand(tot)


def timelike_eig_sym(M):
    """(u0, lambda0, [u_k], [lambda_k], [sigma_k]) of N = M xi, exact."""
    N = M * XI
    us, lams = [], []
    for lam, mult, vecs in N.eigenvects():
        assert mult == 1
        v = vecs[0]
        n2 = (v.T * XI * v)[0, 0]
        v = v / sp.sqrt(sp.Abs(n2))
        us.append(sp.simplify(v))
        lams.append(lam)
    sig = [sp.sign((u.T * XI * u)[0, 0]) for u in us]
    k0 = [k for k in range(4) if sig[k] == -1]
    assert len(k0) == 1, "exactly one timelike eigenvector expected"
    k0 = k0[0]
    return us[k0], lams[k0], us, lams, sig


def du_sym(M, A):
    """d_mu u for jets A (list of 4 matrices), the perturbation formula."""
    u0, l0, us, lams, sig = timelike_eig_sym(M)
    out = []
    for mu in range(4):
        d = sp.zeros(4, 1)
        for k in range(4):
            if lams[k] == l0:
                continue
            c = sig[k] * (us[k].T * XI * A[mu] * XI * u0)[0, 0] / (l0 - lams[k])
            d += c * us[k]
        out.append(sp.simplify(d))
    return u0, out


def I1_h_sym(F, M, p):
    u0 = timelike_eig_sym(M)[0]
    hu = XI_X * u0
    h = XI_X + 2 * hu * hu.T
    tot = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Fm = F[(mu, nu)]
            tot += XI[mu, mu] * XI[nu, nu] * sp.trace(h * Fm * h * Fm.T)
    return sp.expand(tot)


def J1_sym(F, M, p):
    tot = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Fm = F[(mu, nu)]
            tot += XI[mu, mu] * XI[nu, nu] * sp.trace(
                Fm * XI_X * M * XI_X * Fm * XI_X * M * XI_X)
    return sp.expand(tot)


def J2_sym(F, M, p):
    tot = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Fm = F[(mu, nu)]
            tot += XI[mu, mu] * XI[nu, nu] * sp.trace(
                Fm * XI_X * Fm * XI_X * M * XI_X * M * XI_X)
    return sp.expand(tot)


def Pgrad_sym(F, M, p):
    A = F["jets"]
    u0, du = du_sym(M, A)
    tot = 0
    for mu in range(4):
        dP = du[mu] * u0.T * XI + u0 * du[mu].T * XI
        q = sum(XI_X[a, a] * XI_X[b, b] * dP[a, b] ** 2
                for a in range(4) for b in range(4))
        tot += XI[mu, mu] * q
    return sp.expand(tot)


# ================= numpy helpers =================
def timelike_eig_np(M):
    """batched: u0 (...,4), and the full basis U (...,4,4) columns, lams
    (...,4), sig (...,4), the timelike one at column index k0 (...)."""
    N = M @ ETA
    lam, V = np.linalg.eig(N)
    if np.max(np.abs(lam.imag)) > 1e-9 * max(np.max(np.abs(lam.real)), 1.0):
        raise ValueError("complex spectrum of M eta")
    lam = lam.real
    V = V.real
    n2 = np.einsum("...ak,a,...ak->...k", V, np.diag(ETA), V)
    V = V / np.sqrt(np.abs(n2))[..., None, :]
    sig = np.sign(n2)
    k0 = np.argmin(n2, axis=-1)
    u0 = np.take_along_axis(V, k0[..., None, None], axis=-1)[..., 0]
    return u0, V, lam, sig, k0


def du_np(A, M):
    """d_mu u (4, ..., 4) by the perturbation formula."""
    u0, V, lam, sig, k0 = timelike_eig_np(M)
    l0 = np.take_along_axis(lam, k0[..., None], axis=-1)[..., 0]
    eu0 = u0 @ ETA                                        # (eta u0)^T
    out = np.zeros(A.shape[:-1])
    for mu in range(4):
        # c_k = sig_k (u_k^T eta A eta u0) / (l0 - l_k)
        num = np.einsum("...ak,...ab,...b->...k", V, ETA @ A[mu], eu0)
        den = l0[..., None] - lam
        mask = np.abs(den) > 1e-12 * np.maximum(np.abs(l0)[..., None], 1.0)
        c = np.where(mask, sig * num / np.where(mask, den, 1.0), 0.0)
        out[mu] = np.einsum("...ak,...k->...a", V, c)
    return u0, out


def h_cov_np(M):
    u0 = timelike_eig_np(M)[0]
    hu = u0 @ ETA_X
    return ETA_X + 2.0 * hu[..., :, None] * hu[..., None, :]


def _pair_sum(F, fn):
    tot = 0.0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            tot = tot + ETA[mu, mu] * ETA[nu, nu] * fn(F[..., mu, nu, :, :])
    return tot


def I1_h_np(A, M, p):
    F = L0.F_of_A(A)
    h = h_cov_np(M)
    return _pair_sum(F, lambda Fm: np.einsum(
        "...ab,...bc,...cd,...ad->...", h, Fm, h, Fm))


def J1_np(A, M, p):
    F = L0.F_of_A(A)
    ME = ETA_X @ M @ ETA_X
    return _pair_sum(F, lambda Fm: np.einsum(
        "...ab,...bc,...cd,...da->...", Fm, ME, Fm, ME))


def J2_np(A, M, p):
    F = L0.F_of_A(A)
    ME = ETA_X @ M @ ETA_X
    return _pair_sum(F, lambda Fm: np.einsum(
        "...ab,...bc,...cd,...da->...", Fm, ETA_X @ Fm, ME, M @ ETA_X))


def Pgrad_np(A, M, p):
    u0, du = du_np(A, M)
    eu0 = u0 @ ETA
    tot = 0.0
    w = np.diag(ETA_X)
    for mu in range(4):
        edu = du[mu] @ ETA
        dP = du[mu][..., :, None] * eu0[..., None, :] \
            + u0[..., :, None] * edu[..., None, :]
        q = np.einsum("a,b,...ab,...ab->...", w, w, dP, dP)
        tot = tot + ETA[mu, mu] * q
    return tot


# ================= the extended registry =================
class TermX:
    """registry shape of m5_32_lagrangian.Term (duck-typed for its
    term_lagrangian / omega_decompose): density(A, M, p), sympy(F, M, p)."""

    def __init__(self, name, definition, sympy_fn, K=None, np_fn=None,
                 parity="even", kind="curvature_ext"):
        self.name = name
        self.definition = definition
        self.hash = hashlib.sha256(definition.encode()).hexdigest()[:12]
        self.sympy_fn = sympy_fn
        self._K = K
        self._np_fn = np_fn
        self.parity = parity
        self.kind = kind

    def density(self, A, M, p):
        if self._K is not None:
            return L0.density_from_K(L0.F_of_A(A), self._K())
        return self._np_fn(A, M, p)

    def sympy(self, F, M, p):
        if self._K is not None:
            return sym_from_K(F, self._K())
        return self.sympy_fn(F, M, p)


REGISTRY_EXT = {}


def _reg(t):
    REGISTRY_EXT[t.name] = t
    return t


def build_eps_basis(n_jets=60, seed=3210, tol=1e-9):
    """evaluate all 210 patterns on random jets, pick a greedy basis,
    reduce every pattern on it; returns the report dict and registers
    E1..Er."""
    rng = np.random.default_rng(seed)
    p = L0.default_params(s=-1.0, g=32.0)
    A, M = L0._random_jets(rng, n_jets, p)
    F = L0.F_of_A(A)
    pats = eps_patterns()
    vals = []
    Ks = []
    for p1, p2, pe, lab in pats:
        K = K_eps(p1, p2, pe)
        Ks.append(K)
        vals.append(L0.density_from_K(F, K))
    Vm = np.array(vals).T                                   # (jets, 210)
    scale = np.max(np.abs(Vm))
    nz = [k for k in range(len(pats)) if np.max(np.abs(Vm[:, k])) > tol * scale]
    sv = np.linalg.svd(Vm / scale, compute_uv=False)
    rank = int(np.sum(sv > tol))
    basis = []
    for k in nz:
        cand = basis + [k]
        r = np.linalg.matrix_rank(Vm[:, cand] / scale, tol=tol)
        if r == len(cand):
            basis.append(k)
        if len(basis) == rank:
            break
    Bm = Vm[:, basis]
    red = []
    hist = {}
    for k, (p1, p2, pe, lab) in enumerate(pats):
        if k not in nz:
            red.append({"pattern": lab, "status": "vanishes"})
            hist["0"] = hist.get("0", 0) + 1
            continue
        c, *_ = np.linalg.lstsq(Bm, Vm[:, k], rcond=None)
        res = float(np.max(np.abs(Bm @ c - Vm[:, k])) / scale)
        cr = [float(round(x, 8)) for x in c]
        red.append({"pattern": lab, "coeff_on_basis": cr, "residual": res})
        key = " ".join(f"{x:+g}E{i+1}" for i, x in enumerate(cr) if abs(x) > 1e-7)
        hist[key] = hist.get(key, 0) + 1
    for i, k in enumerate(basis):
        p1, p2, pe, lab = pats[k]
        _reg(TermX(f"E{i+1}", f"E{i+1} = F[{p1}] F[{p2}] eps[{pe}] "
                   f"({lab}; slots 0,1 derivative, 2,3 internal; one "
                   "Levi-Civita symbol; the R0 rule extended: eps-d eta, "
                   "eps-i delta)", None, K=(lambda K=Ks[k]: K),
                   parity="odd"))
    # rank together with the even six
    even = np.array([L0.REGISTRY[n].density(A, M, p)
                     for n in ("I1", "I2", "I3", "I4", "I5", "I6")]).T
    both = np.concatenate([even / np.max(np.abs(even)), Bm / scale], axis=1)
    sv2 = np.linalg.svd(both, compute_uv=False)
    return {"n_patterns": len(pats), "n_nonvanishing": len(nz),
            "rank_eps": rank, "singular_values": (sv[:rank + 2]).tolist(),
            "basis_patterns": [pats[k][3] for k in basis],
            "basis_slots": [(pats[k][0], pats[k][1], pats[k][2]) for k in basis],
            "rank_even6_plus_eps": int(np.sum(sv2 > tol)),
            "reduction_histogram": hist,
            "reductions": red}


EPS_REPORT = build_eps_basis()

_reg(TermX("I1_h",
           "I1_h = sum_{mu<nu} eta^mu eta^nu tr(h_cov F h_cov F^T), "
           "h_cov = eta + 2 (eta u)(eta u)^T, u = timelike unit eigenvector "
           "of M eta (contravariant, u^T eta u = -1); = I1_frob in the "
           "vacuum eigenframe", I1_h_sym, np_fn=I1_h_np))
_reg(TermX("J1",
           "J1 = sum_{mu<nu} eta^mu eta^nu tr(F eta M eta F eta M eta) "
           "(contravariant F, M joined by eta; tr((FM)^2)-type)",
           J1_sym, np_fn=J1_np))
_reg(TermX("J2",
           "J2 = sum_{mu<nu} eta^mu eta^nu tr(F eta F eta M eta M eta) "
           "(tr(F^2 M^2)-type)", J2_sym, np_fn=J2_np))
_reg(TermX("Pgrad",
           "Pgrad = sum_mu eta^{mu mu} q(d_mu P_t, d_mu P_t), P_t = u u^T eta, "
           "q(X,Y) = sum_ab eta_a eta_b X[a,b] Y[a,b], d_mu u by first-order "
           "eigenvector perturbation of M eta along A_mu (the search knob "
           "-kappa is the coefficient)", Pgrad_sym, np_fn=Pgrad_np,
           kind="projector_ext"))


def all_terms():
    """the full extended catalog: base six + controls + ext."""
    d = dict(L0.REGISTRY)
    d.update(REGISTRY_EXT)
    return d


# ================= selftests =================
def _rel(a, b):
    return float(abs(a - b) / max(abs(b), 1e-300))


def _transform(L, A, M):
    Linv_T = np.linalg.inv(L).T
    Mp = np.einsum("ab,...bc,dc->...ad", L, M, L)
    Ap = np.einsum("mn,n...ab->m...ab", Linv_T,
                   np.einsum("ab,n...bc,dc->n...ad", L, A, L))
    return Ap, Mp


def st_covariance(res, lines):
    rng = np.random.default_rng(3211)
    p = L0.default_params(s=-1.0, g=32.0)
    A, M = L0._random_jets(rng, 48, p)
    out = {}
    worst = 0.0
    for kind in ("boost", "rotation"):
        for k in range(3):
            Lm = L0._lorentz(rng, kind)
            Ap, Mp = _transform(Lm, A, M)
            for nm, T in REGISTRY_EXT.items():
                d0 = T.density(A, M, p)
                d1 = T.density(Ap, Mp, p)
                drift = float(np.max(np.abs(d1 - d0))
                              / max(np.max(np.abs(d0)), 1e-300))
                out.setdefault(nm, {})[f"{kind}_{k}"] = drift
                worst = max(worst, drift)
    res["covariance"] = {"per_term": out, "worst": worst}
    lines.append(("covariance: every ext term, SO(1,3)+SO(3), density",
                  worst, 1e-10, worst <= 1e-10))


def st_parity(res, lines):
    rng = np.random.default_rng(3212)
    p = L0.default_params(s=-1.0, g=32.0)
    A, M = L0._random_jets(rng, 48, p)
    out = {}
    worst_odd = 0.0
    worst_even = 0.0
    for lab, P in (("x_reflection", np.diag([1.0, -1.0, 1.0, 1.0])),
                   ("time_reversal", np.diag([-1.0, 1.0, 1.0, 1.0]))):
        Ap, Mp = _transform(P, A, M)
        for nm, T in all_terms().items():
            if T.kind == "potential":
                continue
            d0 = T.density(A, M, p)
            d1 = T.density(Ap, Mp, p)
            sc = max(np.max(np.abs(d0)), 1e-300)
            flip = float(np.max(np.abs(d1 + d0)) / sc)
            keep = float(np.max(np.abs(d1 - d0)) / sc)
            par = getattr(T, "parity", "even")
            out.setdefault(nm, {})[lab] = {"flip_residual": flip,
                                           "keep_residual": keep,
                                           "parity": par}
            if par == "odd":
                worst_odd = max(worst_odd, flip)
            else:
                worst_even = max(worst_even, keep)
    res["parity"] = {"per_term": out, "worst_odd_flip": worst_odd,
                     "worst_even_keep": worst_even}
    lines.append(("parity: every eps term flips sign (x-refl, T)",
                  worst_odd, 1e-10, worst_odd <= 1e-10))
    lines.append(("parity: every eta-only term keeps sign (x-refl, T)",
                  worst_even, 1e-10, worst_even <= 1e-10))


def rational_lorentz():
    """a rational SO(1,3) element: boost (cosh 5/4, sinh 3/4) in (0,1)
    times a rotation (3/5, 4/5) in (2,3)."""
    Bm = sp.eye(4)
    Bm[0, 0] = Bm[1, 1] = sp.Rational(5, 4)
    Bm[0, 1] = Bm[1, 0] = sp.Rational(3, 4)
    Rm = sp.eye(4)
    Rm[2, 2] = Rm[3, 3] = sp.Rational(3, 5)
    Rm[2, 3], Rm[3, 2] = sp.Rational(-4, 5), sp.Rational(4, 5)
    Lr = Bm * Rm
    assert Lr.T * XI * Lr == XI
    return Lr


def st_sympy_vs_numpy(res, lines):
    t0 = time.time()
    rng = np.random.default_rng(3213)
    Lr = rational_lorentz()
    out = {}
    worst = 0.0
    for trial in range(2):
        def rsym():
            X = sp.zeros(4, 4)
            for i in range(4):
                for j in range(i, 4):
                    X[i, j] = X[j, i] = sp.Rational(
                        int(rng.integers(-9, 10)), int(rng.integers(1, 7)))
            return X
        Mmu = [rsym() for _ in range(4)]
        D = sp.diag(32, 1, sp.Rational(3, 10), 0)       # s = -1 vacuum
        Mpt = Lr * D * Lr.T
        F = L0.F_of_jets_sym(Mmu)
        F["jets"] = Mmu
        A_np = np.array([[[float(Mmu[m][i, j]) for j in range(4)]
                          for i in range(4)] for m in range(4)])[:, None]
        M_np = np.array([[float(Mpt[i, j]) for j in range(4)]
                         for i in range(4)])[None]
        p = L0.default_params(s=-1.0, g=32.0, w=1.0)
        for nm, T in REGISTRY_EXT.items():
            vs = float(sp.nsimplify(T.sympy(F, Mpt, p)))
            vn = float(T.density(A_np, M_np, p)[0])
            rel = _rel(vn, vs)
            out.setdefault(nm, []).append({"sympy": vs, "numpy": vn, "rel": rel})
            worst = max(worst, rel)
    res["sympy_vs_numpy"] = {"per_term": out, "worst": worst,
                             "runtime_s": round(time.time() - t0, 1)}
    lines.append(("sympy vs numpy, every ext term, 2 rational jets",
                  worst, 1e-12, worst <= 1e-12))


def st_vacuum_frame(res, lines):
    """I1_h == I1_frob at M = vacuum (u = e0); and I1_h(boosted) ==
    I1_frob(unboosted) (the covariantized flip)."""
    rng = np.random.default_rng(3214)
    p = L0.default_params(s=-1.0, g=32.0)
    cfg = B3.base_cfg(s=-1.0, g=32.0)
    n = 40
    M = np.broadcast_to(B3.vac4(cfg), (n, 4, 4)).copy()
    A = B3.sym4(rng.normal(size=(4, n, 4, 4)))
    dh = REGISTRY_EXT["I1_h"].density(A, M, p)
    df = L0.REGISTRY["I1_frob"].density(A, M, p)
    de = L0.REGISTRY["I1"].density(A, M, p)
    r1 = float(np.max(np.abs(dh - df)) / np.max(np.abs(df)))
    Lm = L0._lorentz(rng, "boost", 0.7)
    Ap, Mp = _transform(Lm, A, M)
    dhb = REGISTRY_EXT["I1_h"].density(Ap, Mp, p)
    dfb = L0.REGISTRY["I1_frob"].density(Ap, Mp, p)
    r2 = float(np.max(np.abs(dhb - df)) / np.max(np.abs(df)))
    r3 = float(np.max(np.abs(dfb - df)) / np.max(np.abs(df)))
    res["vacuum_frame"] = {"I1_h_vs_I1_frob_at_vac": r1,
                           "I1_h_boosted_vs_I1_frob_unboosted": r2,
                           "I1_frob_boosted_vs_unboosted_(control, breaks)": r3,
                           "I1_h_vs_I1_at_vac_rel_diff": float(
                               np.max(np.abs(dh - de)) / np.max(np.abs(de)))}
    lines.append(("I1_h == I1_frob in the vacuum eigenframe", r1, 1e-12,
                  r1 <= 1e-12))
    lines.append(("I1_h(boosted) == I1_frob(vacuum frame) (covariantized flip)",
                  r2, 1e-10, r2 <= 1e-10))
    lines.append(("control: I1_frob itself breaks under the boost", r3,
                  1e-3, r3 > 1e-3))


def st_du(res, lines):
    rng = np.random.default_rng(3215)
    p = L0.default_params(s=-1.0, g=32.0)
    A, M = L0._random_jets(rng, 16, p)
    u0, du = du_np(A, M)
    worst = 0.0
    t = 1e-6
    for mu in range(4):
        up = timelike_eig_np(M + t * A[mu])[0]
        um = timelike_eig_np(M - t * A[mu])[0]
        # fix the sign freedom of eig by aligning to u0
        up = up * np.sign(np.einsum("na,na->n", up, u0))[:, None]
        um = um * np.sign(np.einsum("na,na->n", um, u0))[:, None]
        fd = (up - um) / (2 * t)
        worst = max(worst, float(np.max(np.abs(fd - du[mu]))
                                 / max(np.max(np.abs(du[mu])), 1e-300)))
    # normalization: u^T eta u = -1 and u^T eta du = 0
    n0 = float(np.max(np.abs(np.einsum("na,ab,nb->n", u0, ETA, u0) + 1.0)))
    n1 = float(np.max(np.abs(np.einsum("na,ab,mnb->mn", u0, ETA, du))))
    res["du"] = {"fd_vs_formula_rel": worst, "u_eta_u_plus_1": n0,
                 "u_eta_du": n1}
    lines.append(("d_mu u: perturbation formula vs central difference",
                  worst, 1e-6, worst <= 1e-6))
    lines.append(("u normalization: u eta u = -1, u eta du = 0",
                  max(n0, n1), 1e-10, max(n0, n1) <= 1e-10))


def st_eps_rank(res, lines):
    r = EPS_REPORT
    res["eps"] = {k: v for k, v in r.items() if k != "reductions"}
    res["eps"]["max_reduction_residual"] = max(
        x.get("residual", 0.0) for x in r["reductions"])
    ok = r["rank_even6_plus_eps"] == 6 + r["rank_eps"]
    lines.append(("eps: basis independent of the even six (rank 6 + r)",
                  0.0 if ok else 1.0, 0.5, ok))
    lines.append(("eps: every pattern reduces on the basis",
                  res["eps"]["max_reduction_residual"], 1e-10,
                  res["eps"]["max_reduction_residual"] <= 1e-10))


def st_mutation(res, lines):
    cmd = [PY, os.path.abspath(__file__), "--selftest", "--mutant",
           "eta_time_row", "--no-mutation-stage", "--no-write"]
    pr = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    reds = [ln for ln in pr.stdout.splitlines() if ln.startswith("FAIL")]
    res["mutation"] = {"mutant": "eta_time_row", "n_fail": len(reds),
                       "failed_lines": reds, "returncode": pr.returncode}
    lines.append(("mutation: eta_time_row mutant reddens >= 1 line",
                  float(len(reds)), 1.0, len(reds) >= 1))


def selftest(mutation_stage=True, write=True):
    t0 = time.time()
    res = {"mutant": MUTANT, "python": PY,
           "terms": {nm: {"definition": T.definition, "hash": T.hash,
                          "kind": T.kind, "parity": T.parity}
                     for nm, T in REGISTRY_EXT.items()}}
    lines = []
    for st in (st_covariance, st_parity, st_sympy_vs_numpy, st_vacuum_frame,
               st_du, st_eps_rank):
        st(res, lines)
    if mutation_stage:
        st_mutation(res, lines)
    res["selftest_lines"] = [{"name": n, "value": v, "threshold": t,
                              "pass": bool(ok)} for n, v, t, ok in lines]
    res["all_pass"] = all(ok for *_, ok in lines)
    res["runtime_s"] = round(time.time() - t0, 1)
    for n, v, t, ok in lines:
        print(f"{'PASS' if ok else 'FAIL'}  {n:62s} {v:.3e}  (thr {t:g})")
    print(f"ALL_PASS={res['all_pass']}  mutant={MUTANT}  "
          f"runtime {res['runtime_s']} s")
    if write:
        with open(os.path.join(DATA, "m5_32_r1_ext_selftest.json"), "w") as f:
            json.dump(res, f, indent=1)
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--mutant", default=None)
    ap.add_argument("--no-mutation-stage", action="store_true")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args()
    if a.mutant:
        set_mutant(a.mutant)
        EPS_REPORT = build_eps_basis()      # rebuild under the mutant metric
    for nm, T in REGISTRY_EXT.items():
        print(f"{nm:8s} {T.hash}  {T.definition[:90]}")
    print(f"eps rank = {EPS_REPORT['rank_eps']} "
          f"(even6+eps rank {EPS_REPORT['rank_even6_plus_eps']})")
    if a.selftest:
        r = selftest(mutation_stage=not a.no_mutation_stage,
                     write=not a.no_write)
        sys.exit(0 if r["all_pass"] else 1)
