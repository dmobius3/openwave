"""M5.32 R20-0 adversarial audit: an independent re-derivation of the six
claim groups of the R20-0 form-level script (m5_32_r20_0_class.py and its
data/m5_32_r20_0_class.json), written from the definitions with its own
methods (exact sympy Hessians by bivariate polynomial coefficients, exact
sympy gradients by univariate polynomial coefficients, eigenvalue power
sums instead of matrix powers, Horner evaluation of the expanded spectral
polynomial instead of the product of factors, its own one-sided stencils,
its own jets, basis and eps ladder), so that agreement with the producer
is a check and not a re-run. This file does NOT call any function of
m5_32_r20_0_class.py for a check; it imports that module ONLY to read the
gradient it ships (compared against this file's own exact gradient) and
reads its JSON for the numbers under audit. The stack modules
(m5_21_3_a_4d.py as B3, m5_32_lagrangian.py, m5_32_r19_entrants.py) are
imported for the certified energies (B3.e_parts, EN.block_reads), which
are the objects the claims are about.

EQUATIONS FIRST
---------------
Field M real symmetric 4x4, eta = diag(-1, 1, 1, 1), N = M eta. The stack
(code branch s = -1, g = 8, delta = 0.3): B3.vac4 = diag(-s g, 1, delta, 0)
= diag(8, 1, 0.3, 0), so N has the spectrum q = (s g, 1, delta, 0)
= (-8, 1, 0.3, 0) and B3.c4_of = (sum_i q_i^p)_{p=1..4}.
    V4     = w sum_{p=1..4} (tr N^p - C_p)^2,  C_p = sum_i q_i^p
    V_spec = gamma tr[P(N)^2],  P(x) = prod_i (x - q_i)

Claim 1 (Hessians on sym4 at M_vac). For a basis B_a of the 10-dim sym4
space, H_ab is the coefficient of t s in V(M_vac + t B_a + s B_b)
(twice the coefficient of t^2 when a = b), computed EXACTLY in sympy with
rational roots (delta = 3/10). Closed forms to compare against:
    H_V4 (diagonal directions)     = D (2 w J^T J) D,  J_pi = p q_i^(p-1),  D = eta
    H_V_spec (diagonal directions) = 2 gamma diag(P'(q_i)^2),  P'(q_i) = prod_{j != i} (q_i - q_j)
and zero on the six off-diagonal directions and on the cross block. The
diagonal V4 form follows from n_i = eta_i m_i (d/dm_i = eta_i d/dn_i) and
from tr N^p - C_p = 0 at the vacuum (only the gradient-times-gradient term
survives); the V_spec form from V_spec(diag) = gamma sum_i P(q_i + d_i)^2.
The off-diagonal directions are the Lorentz conjugation directions
dM = K M + M K^T with K = eta A, A antisymmetric (K eta + eta K^T = 0), on
which both potentials are constant (N -> O N O^-1); they are nonzero
directions iff the four N eigenvalues are distinct.
The author's four numbers 371866.88 / 48.02 / 5.229 / 11.52 are tested
against 2 P'(q_i)^2 on (8, 1, 0.3, 0) and on (-8, 1, 0.3, 0); and
P'(q_i)^2 under the full reflection q -> -q (invariant, each factor flips
sign, three factors) and under the single flip q_0 -> -q_0 (not invariant,
which is what makes the two root sets distinguishable).

Claim 2 (gradient of V_spec). d tr[P(N)^2] = 2 tr[P(N) P'(N) dN]
= 2 tr[eta P P' dM], so with dE = sum_ab G_ab dM_ab on symmetric dM,
G = 2 gamma sym4[(eta P P')^T]. Tested by the exact univariate route: at a
rational M0 the coefficient of t in tr[P(N(M0 + t B_a))^2] along the 10
basis directions (E_ii, and E_ij + E_ji giving 2 G_ij), against the
formula evaluated exactly in sympy (rational matrices). The mutants: eta
dropped, eta on the right (P P' eta), the factor 2 dropped, P' replaced by
P (all must FAIL); the transpose omitted is tested too and is NOT a
failing mutant, because sym4(X^T) = sym4(X). A complex-step check of the
same formula on a random numpy field closes the numeric side, and the
shipped vspec_energy_grad of the class script is compared with this
file's own numpy formula (Horner on the expanded polynomial).

Claim 3 (gamma). gamma = E_V4 / (h^3 sum_cells tr[P(N)^2]) on the record
field data/m5_32_r19_1/I1_un_single_d0_n32_g8.npz (n 32, L 48, s = -1,
g = 8, delta = 0.3). E_V4 recomputed here from the per-cell EIGENVALUES of
N (power sums of eigenvalues, np.linalg.eigvals) and compared with
B3.e_parts; sum tr[P(N)^2] recomputed from the eigenvalues
(tr P(N)^2 = sum_i P(lambda_i)^2 for any matrix) and from Horner. The
record's total 13.8349 is rebuilt as 4 h^3 sum I1 + V4 with this file's
own fwd / bwd stencils. The s branch of the record is read off the field
itself (the far-field M_00) and off V4 under the other branch's C_p.

Claim 4 (flat-vacuum jets). E_I1 = EN.block_reads(M, cfg, "I1")["E_curv_I1"]
on n = 8, L = 12 with M = M_vac + eps X(x): F_ij = A_i eta A_j - A_j eta A_i
is bilinear in the jets and A(M_vac) = 0, so E_I1 = eps^4 E_I1[X] EXACTLY
for every shape X (linear or not in x): the exponent must be 4 to
roundoff, single-axis jets must be exactly zero (one axis, no pair), and a
two-axis pair x B_a + y B_b is zero iff B_a eta B_b - B_b eta B_a = 0
(zero up to the entrants' roundoff floor: a same-matrix random pair
returns about 1e-29 where this file's direct F is exactly 0). Own
random symmetric jets, own eps ladder (0.3, 0.1, 0.03, 0.01), own quadratic
and mixed shapes; E_I1 also rebuilt with this file's own stencils. The
class script's 12 pairs (its JSON) are re-classified by the commutator
criterion. Qualifier recorded: the second variation of L_cert about the
flat vacuum is the mass term -H_V4 (nonzero on the four spectral
directions), zero only in the derivative sector.

Claim 5 (pair laws on record). Each number of the class JSON's stage (e)
is searched in the ledger text section by section, the six-point law is
refitted from the R19-3 table of the task record, the R3 lambda = 0 rows
are read from the R19-1 JSON, the R19-1 (Gam) tail is refitted as
A + B / d (R^2) and the undressed g 8 pair's E_int(d) trend is read from
the R19-1 JSON.

Claim 6 (ledger section 6.9 rows 1 and 2). Row 1's "(8, 1, 0.3, 0): the
spectrum of M eta in our embedding (M_00 = -g)" against the code branch;
row 2's V4 eigenvalues 0.059 / 1.50 / 12.6 / 4.23e6 against the exact
eigenvalues of 2 J^T J on both root sets (the ratio test: a constant
factor 2 on (8, 1, 0.3, 0), no constant factor on (-8, 1, 0.3, 0)).

Every printed line is PASS/FAIL on a number that can go either way.
Out: ../data/m5_32_r20_0_audit.json
Usage: python3 m5_32_r20_0_audit.py
"""
from __future__ import annotations

import importlib.util
import json
import os
import re
import sys
import time

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
FINDINGS = os.path.join(HERE, "..", "findings")
TASKS = os.path.join(HERE, "..", "tasks")
OUT = os.path.join(DATA, "m5_32_r20_0_audit.json")
CLASS_JSON = os.path.join(DATA, "m5_32_r20_0_class.json")
LEDGER = os.path.join(FINDINGS, "m5_32_candidate_ledger.md")
TASK_MD = os.path.join(TASKS, "m5_32_task_details.md")
RECORD_SINGLE = os.path.join(DATA, "m5_32_r19_1", "I1_un_single_d0_n32_g8.npz")
R19_1_JSON = os.path.join(DATA, "m5_32_r19_1_pair.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    argv = sys.argv
    sys.argv = [argv[0]]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


EN = _load("m5_32_r19_entrants", "m5_32_r19_entrants.py")   # block_reads: the certified curvature energy under test
LAG = EN.LAG
B3 = LAG.B3                                                   # e_parts, base_cfg, vac4, c4_of, W1: the certified stack
CLS = _load("m5_32_r20_0_class", "m5_32_r20_0_class.py")     # ONLY its shipped gradient is read, for comparison

E4 = np.array([-1.0, 1.0, 1.0, 1.0])
ETA = np.diag(E4)
ETA_S = sp.diag(-1, 1, 1, 1)
G8, DELTA = 8.0, 0.3
T0 = time.time()
LINES = {}


def log(*a):
    print(f"[{time.time() - T0:6.1f}s]", *a, flush=True)


def line(key, ok, detail):
    LINES[key] = {"pass": bool(ok), "detail": detail}
    print(f"{'PASS' if ok else 'FAIL'} {key}: {detail}", flush=True)


def rel(a, b):
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))


# ================= own exact (sympy) layer =================
def sym_basis_exact():
    """10 sym4 directions: E_ii (4), then (E_ij + E_ji) / sqrt 2 (6, lexicographic), exact."""
    B = []
    for i in range(4):
        M = sp.zeros(4, 4); M[i, i] = 1; B.append(M)
    for i in range(4):
        for j in range(i + 1, 4):
            M = sp.zeros(4, 4); M[i, j] = M[j, i] = sp.sqrt(2) / 2; B.append(M)
    return B


def V4_exact(M, q, w=1):
    N = M * ETA_S
    Np = sp.eye(4); V = 0
    for p in range(1, 5):
        Np = Np * N
        V += w * (Np.trace() - sum(qi ** p for qi in q)) ** 2
    return V


def Vspec_exact(M, q, gamma=1):
    N = M * ETA_S
    P = sp.eye(4)
    for qi in q:
        P = P * (N - qi * sp.eye(4))
    return gamma * (P * P).trace()


def exact_hessian(pot, Mvac, q):
    """H_ab = coefficient of t s in pot(Mvac + t B_a + s B_b); 2 x coefficient of t^2 on the diagonal."""
    t, s = sp.symbols("t s")
    B = sym_basis_exact()
    H = sp.zeros(10, 10)
    for a in range(10):
        for b in range(a, 10):
            V = sp.expand(pot(Mvac + t * B[a] + s * B[b], q))
            poly = sp.Poly(V, t, s)
            if a == b:
                H[a, a] = 2 * poly.coeff_monomial(t ** 2)
            else:
                H[a, b] = H[b, a] = poly.coeff_monomial(t * s)
    return H


def closed_V4(q, w=1):
    J = sp.Matrix(4, 4, lambda p, i: (p + 1) * q[i] ** p)
    return ETA_S * (2 * w * J.T * J) * ETA_S


def Pprime_exact(q):
    return [sp.prod([q[i] - q[j] for j in range(4) if j != i]) for i in range(4)]


# ================= own numpy layer =================
def poly_coeffs(q):
    """coefficients (highest first) of P(x) = prod (x - q_i) and of P'(x)."""
    c = np.array([1.0])
    for qi in q:
        c = np.convolve(c, [1.0, -qi])
    dc = c[:-1] * np.arange(len(c) - 1, 0, -1)
    return c, dc


def horner(N, c):
    """Horner evaluation of a matrix polynomial on a (..., 4, 4) stack."""
    I = np.broadcast_to(np.eye(4), N.shape)
    R = c[0] * I
    for ck in c[1:]:
        R = R @ N + ck * I
    return R


def own_vspec_energy(M, q, gamma, h3):
    c, _ = poly_coeffs(q)
    P = horner(M @ ETA, c)
    return gamma * h3 * np.sum(np.einsum("...ij,...ji->...", P, P))


def own_vspec_grad(M, q, gamma, h3, variant="correct"):
    c, dc = poly_coeffs(q)
    N = M @ ETA
    P, Pp = horner(N, c), horner(N, dc)
    X = P @ Pp
    if variant == "correct":
        G = 2.0 * (ETA @ X).swapaxes(-1, -2)
    elif variant == "no_transpose":
        G = 2.0 * (ETA @ X)
    elif variant == "eta_dropped":
        G = 2.0 * X.swapaxes(-1, -2)
    elif variant == "eta_right":
        G = 2.0 * (X @ ETA).swapaxes(-1, -2)
    elif variant == "factor_1":
        G = 1.0 * (ETA @ X).swapaxes(-1, -2)
    elif variant == "Pp_as_P":
        G = 2.0 * (ETA @ (P @ P)).swapaxes(-1, -2)
    else:
        raise ValueError(variant)
    return gamma * h3 * 0.5 * (G + G.swapaxes(-1, -2))


def own_d1(f, ax, h, side):
    out = np.zeros_like(f)
    sl_hi = [slice(None)] * f.ndim; sl_lo = [slice(None)] * f.ndim
    sl_hi[ax] = slice(1, None); sl_lo[ax] = slice(0, -1)
    diff = (f[tuple(sl_hi)] - f[tuple(sl_lo)]) / h
    tgt = [slice(None)] * f.ndim
    tgt[ax] = slice(0, -1) if side == "fwd" else slice(1, None)
    out[tuple(tgt)] = diff
    return out


def own_E_I1(M, h):
    """4 h^3 sum_cells (1/2)(fwd + bwd) sum_{i<j} <F_ij, F_ij>_eta with <F, G>_eta = sum_ab eta_a eta_b F_ab G_ab."""
    tot = 0.0
    for side in ("fwd", "bwd"):
        A = [own_d1(M, ax, h, side) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                tot += 0.5 * np.einsum("a,b,...ab,...ab->", E4, E4, F, F, optimize=True)
    return 4.0 * h ** 3 * tot


def own_V4_eig(M, q, w, h3):
    """V4 from the per-cell eigenvalues of N (power sums of eigenvalues, no matrix powers)."""
    lam = np.linalg.eigvals(M @ ETA)
    V = 0.0
    for p in range(1, 5):
        V = V + (np.sum(lam ** p, axis=-1).real - sum(qi ** p for qi in q)) ** 2
    return w * h3 * float(np.sum(V)), float(np.max(np.abs(lam.imag)))


# ================= claim 1 =================
def audit_1(cls):
    out = {}
    quoted = np.array([371866.88, 48.02, 5.229, 11.52])
    third = sp.Rational(3, 10)
    for label, s, key in (("stack_s_minus_1", -1, "stack_branch_s_minus_1"), ("report_s_plus_1", 1, "report_roots_s_plus_1")):
        cfg = B3.base_cfg(s=float(s), g=G8, n=2, L=2.0, delta=DELTA)
        Mv_code = B3.vac4(cfg)
        q_code = np.diag(Mv_code @ ETA)
        c4_code = np.array(B3.c4_of(cfg))
        q = [sp.Integer(s) * 8, sp.Integer(1), third, sp.Integer(0)]
        Mvac = sp.diag(-sp.Integer(s) * 8, 1, third, 0)
        H4 = exact_hessian(V4_exact, Mvac, q)
        Hs = exact_hessian(Vspec_exact, Mvac, q)
        H4c = closed_V4(q)
        Pp = Pprime_exact(q)
        Hsc = sp.diag(*[2 * v ** 2 for v in Pp])
        dev4 = (H4[:4, :4] - H4c).applyfunc(sp.simplify)
        devs = (Hs[:4, :4] - Hsc).applyfunc(sp.simplify)
        off4 = max(abs(float(H4[a, b])) for a in range(10) for b in range(10) if a >= 4 or b >= 4)
        offs = max(abs(float(Hs[a, b])) for a in range(10) for b in range(10) if a >= 4 or b >= 4)
        eig4 = sorted(np.linalg.eigvalsh(np.array(H4[:4, :4], dtype=float)))
        eigs_own = sorted(float(v) for v in Hs.diagonal()[:4])
        r = {"q_code": q_code.tolist(), "M_vac_code": np.diag(Mv_code).tolist(), "c4_code": c4_code.tolist(),
             "c4_from_q": [float(sum(qi ** p for qi in q_code)) for p in range(1, 5)],
             "V4_exact_diag_block_minus_closed_form_is_zero": bool(dev4 == sp.zeros(4, 4)),
             "V4_offdiag_and_cross_max_abs_exact": off4,
             "V4_diag_block_eigs_own": eig4, "V4_diag_block_eigs_script": cls["a"][key]["V4_w1"]["diag_block_eigs"],
             "V4_H00_exact": float(H4[0, 0]),
             "V4_rank": int(sp.Matrix(H4[:4, :4]).rank()),
             "Vspec_exact_diag_block_minus_closed_form_is_zero": bool(devs == sp.zeros(4, 4)),
             "Vspec_offdiag_and_cross_max_abs_exact": offs,
             "Vspec_diag_own": [float(v) for v in Hs.diagonal()[:4]], "Vspec_diag_script": cls["a"][key]["V_spec_gamma1"]["diag_block_diag"],
             "Pprime_own": [float(v) for v in Pp], "Pprime_script": cls["a"][key]["V_spec_gamma1"]["P_prime_at_roots"]}
        r["quoted_max_rel_dev"] = float(np.max(np.abs(np.array(r["Vspec_diag_own"]) - quoted) / quoted))
        out[label] = r
        ok_code = (np.allclose(c4_code, r["c4_from_q"]) and np.allclose(q_code, [s * 8.0, 1.0, 0.3, 0.0]))
        line(f"A1_{label}_code_vacuum_and_C_p_are_the_power_sums_of_q", ok_code,
             f"vac4 diag {np.diag(Mv_code).tolist()}, N spectrum {q_code.tolist()}, c4_of {c4_code.tolist()} vs power sums {np.round(r['c4_from_q'], 6).tolist()}")
        line(f"A1_{label}_V4_exact_hessian_equals_D_2JtJ_D_and_zero_off_the_diagonal_block",
             r["V4_exact_diag_block_minus_closed_form_is_zero"] and off4 == 0.0 and r["V4_rank"] == 4,
             f"exact diag-block deviation zero {r['V4_exact_diag_block_minus_closed_form_is_zero']}, off/cross max {off4}, rank {r['V4_rank']}, H_00 {r['V4_H00_exact']:.1f}")
        line(f"A1_{label}_V4_eigs_match_script_1e-8", rel(eig4, r["V4_diag_block_eigs_script"]) < 1e-8,
             f"own {[f'{v:.6g}' for v in eig4]} vs script {[f'{v:.6g}' for v in r['V4_diag_block_eigs_script']]}, rel {rel(eig4, r['V4_diag_block_eigs_script']):.1e}")
        line(f"A1_{label}_Vspec_exact_hessian_is_2Pprime2_diag_and_zero_off_the_diagonal_block",
             r["Vspec_exact_diag_block_minus_closed_form_is_zero"] and offs == 0.0,
             f"exact diag-block deviation zero {r['Vspec_exact_diag_block_minus_closed_form_is_zero']}, off/cross max {offs}, diag {[f'{v:.6g}' for v in r['Vspec_diag_own']]}")
        line(f"A1_{label}_Vspec_diag_and_Pprime_match_script_1e-8",
             rel(r["Vspec_diag_own"], r["Vspec_diag_script"]) < 1e-8 and rel(r["Pprime_own"], r["Pprime_script"]) < 1e-12,
             f"rel diag {rel(r['Vspec_diag_own'], r['Vspec_diag_script']):.1e}, rel P' {rel(r['Pprime_own'], r['Pprime_script']):.1e}")
    dp, dm = out["report_s_plus_1"]["quoted_max_rel_dev"], out["stack_s_minus_1"]["quoted_max_rel_dev"]
    line("A2_author_numbers_reproduce_on_8_1_0.3_0_within_1e-4_and_NOT_on_-8_1_0.3_0", dp < 1e-4 and dm > 0.5,
         f"max rel dev on (8, 1, 0.3, 0) {dp:.1e}; on (-8, 1, 0.3, 0) {dm:.2f} (own (-8) numbers {[f'{v:.5g}' for v in out['stack_s_minus_1']['Vspec_diag_own']]})")
    ledger_minus = np.array([714251, 79.4, 6.08, 11.52])
    dl = float(np.max(np.abs(np.array(out["stack_s_minus_1"]["Vspec_diag_own"]) - ledger_minus) / ledger_minus))
    line("A2b_ledger_6.9_numbers_for_-g_1_delta_0_reproduce_within_1e-3", dl < 1e-3, f"714251 / 79.4 / 6.08 / 11.52 vs own, max rel dev {dl:.1e}")
    # the reflection note: full reflection invariant, single-root flip not
    rng = np.random.default_rng(20)
    worst_full, worst_single = 0.0, 0.0
    for _ in range(50):
        qq = rng.standard_normal(4) * 3
        pp = lambda q: np.array([np.prod([q[i] - q[j] for j in range(4) if j != i]) for i in range(4)])  # noqa: E731
        worst_full = max(worst_full, rel(pp(-qq) ** 2, pp(qq) ** 2))
        q1 = qq.copy(); q1[0] = -q1[0]
        worst_single = max(worst_single, rel(pp(q1) ** 2, pp(qq) ** 2))
    out["reflection"] = {"full_reflection_max_rel_change": worst_full, "single_root_flip_min_over_50_max_rel_change": worst_single}
    line("A3_Pprime2_invariant_under_full_q_to_minus_q_but_not_under_a_single_root_flip", worst_full < 1e-12 and worst_single > 0.1,
         f"full reflection max rel change {worst_full:.1e}; single flip q_0 -> -q_0 rel change {worst_single:.2f}: the (8) vs (-8) sets are a SINGLE flip, so the invariance note does not identify them")
    return out


# ================= claim 2 =================
def audit_2(cls):
    out = {}
    rng = np.random.default_rng(2020)
    third = sp.Rational(3, 10)
    q = [sp.Integer(-8), sp.Integer(1), third, sp.Integer(0)]
    # a rational symmetric M0 near the vacuum
    R = sp.zeros(4, 4)
    for i in range(4):
        for j in range(i, 4):
            R[i, j] = R[j, i] = sp.Rational(int(rng.integers(-9, 10)), 7)
    M0 = sp.diag(8, 1, third, 0) + R
    t = sp.symbols("t")
    unit = []
    for i in range(4):
        B = sp.zeros(4, 4); B[i, i] = 1; unit.append((i, i, B))
    for i in range(4):
        for j in range(i + 1, 4):
            B = sp.zeros(4, 4); B[i, j] = B[j, i] = 1; unit.append((i, j, B))
    G_exact = sp.zeros(4, 4)
    for i, j, B in unit:
        V = sp.expand(Vspec_exact(M0 + t * B, q))
        c1 = sp.Poly(V, t).coeff_monomial(t)
        if i == j:
            G_exact[i, i] = c1
        else:
            G_exact[i, j] = G_exact[j, i] = c1 / 2          # dE = sum_ab G_ab dM_ab counts (i, j) and (j, i)
    # the formula, exactly
    N0 = M0 * ETA_S
    P = sp.eye(4); Pp = sp.zeros(4, 4)
    for qi in q:
        P = P * (N0 - qi * sp.eye(4))
    for k in range(4):
        Rk = sp.eye(4)
        for i in range(4):
            if i != k:
                Rk = Rk * (N0 - q[i] * sp.eye(4))
        Pp += Rk
    X = P * Pp
    def symm(Y):
        return (Y + Y.T) / 2
    variants = {"correct": 2 * symm((ETA_S * X).T), "no_transpose": 2 * symm(ETA_S * X), "eta_dropped": 2 * symm(X.T),
                "eta_right": 2 * symm((X * ETA_S).T), "factor_1": symm((ETA_S * X).T), "Pp_as_P": 2 * symm((ETA_S * P * P).T)}
    scale = max(abs(float(v)) for v in G_exact)
    devs = {k: float(max(abs(float(v)) for v in (G_exact - V))) / scale for k, V in variants.items()}
    out["exact_rational_point"] = {"M0": [[str(M0[i, j]) for j in range(4)] for i in range(4)], "G_scale": scale, "rel_dev_by_variant": devs}
    line("B1_exact_gradient_equals_2sym4_eta_P_Pprime_T_exactly", devs["correct"] == 0.0, f"exact rel deviation {devs['correct']} (scale {scale:.3e})")
    line("B2_mutants_eta_dropped_eta_right_factor_1_Pp_as_P_all_FAIL_by_more_than_1e-2",
         all(devs[k] > 1e-2 for k in ("eta_dropped", "eta_right", "factor_1", "Pp_as_P")),
         "rel deviations " + ", ".join(f"{k} {devs[k]:.3f}" for k in ("eta_dropped", "eta_right", "factor_1", "Pp_as_P")))
    line("B3_transpose_omitted_is_NOT_a_failing_mutant_sym4_absorbs_it", devs["no_transpose"] == 0.0,
         f"rel deviation with the transpose omitted {devs['no_transpose']} (sym4(X^T) = sym4(X)): the class docstring's transpose is redundant, not load-bearing")
    # numeric: complex step on this file's own energy, on a random field; and the shipped gradient vs this file's formula
    cfg = B3.base_cfg(s=-1.0, g=G8, n=5, L=7.5, delta=DELTA)
    h3 = cfg["h"] ** 3
    qn = (-8.0, 1.0, 0.3, 0.0)
    M = B3.vac4(cfg)[None, None, None] + 0.3 * B3.sym4(rng.standard_normal((5, 5, 5, 4, 4)))
    G = own_vspec_grad(M, qn, 1.0, h3)
    worst = 0.0
    for _ in range(8):
        D = B3.sym4(rng.standard_normal(M.shape))
        dd = float(np.sum(G * D))
        cs = np.imag(own_vspec_energy(M + 1e-20j * D, qn, 1.0, h3)) / 1e-20
        worst = max(worst, abs(cs - dd) / abs(dd))
    Gcls = CLS.vspec_energy_grad(M, cfg, qn, 1.0)[1]
    Ecls = CLS.vspec_energy_grad(M, cfg, qn, 1.0, need_grad=False)[0]
    out["numeric"] = {"complex_step_worst_rel": worst, "shipped_grad_vs_own_rel": rel(Gcls, G), "shipped_E_vs_own_rel": abs(Ecls - own_vspec_energy(M, qn, 1.0, h3)) / abs(Ecls)}
    line("B4_own_formula_complex_step_1e-12_on_a_random_field", worst < 1e-12, f"worst rel {worst:.1e} over 8 directions, n 5")
    line("B5_shipped_vspec_energy_grad_matches_own_horner_formula_1e-12",
         out["numeric"]["shipped_grad_vs_own_rel"] < 1e-12 and out["numeric"]["shipped_E_vs_own_rel"] < 1e-12,
         f"grad rel {out['numeric']['shipped_grad_vs_own_rel']:.1e}, E rel {out['numeric']['shipped_E_vs_own_rel']:.1e}")
    return out


# ================= claim 3 =================
def audit_3(cls):
    out = {}
    cfg = B3.base_cfg(s=-1.0, g=G8, n=32, L=48.0, delta=DELTA)
    h = cfg["h"]; h3 = h ** 3
    q = (-8.0, 1.0, 0.3, 0.0)
    M = np.load(RECORD_SINGLE)["M"]
    _, ev = B3.e_parts(M, cfg)
    ev = float(ev)
    V4_eig, max_imag = own_V4_eig(M, q, B3.W1, h3)
    lam = np.linalg.eigvals(M @ ETA)
    c, _ = poly_coeffs(q)
    S_eig = h3 * float(np.sum(np.polyval(c, lam) ** 2).real)
    S_hor = own_vspec_energy(M, q, 1.0, h3)
    gamma_own = ev / S_eig
    E_I1_own = own_E_I1(M, h)
    E_tot_own = E_I1_own + ev
    br = EN.block_reads(M, cfg, "I1")
    # the s branch of the record: far-field M_00 and V4 under the other branch's targets
    X, Y, Z = B3.coords(32, h)
    far = np.sqrt(X * X + Y * Y + Z * Z) > 20.0
    m00_far = float(np.mean(M[far][:, 0, 0]))
    V4_other, _ = own_V4_eig(M, (8.0, 1.0, 0.3, 0.0), B3.W1, h3)
    out.update({"E_V4_stack": ev, "E_V4_own_eigenvalue_route": V4_eig, "max_imag_eigenvalue": max_imag, "E_V4_stack_expected": 0.3069567721158226,
                "sum_h3_trP2_eig": S_eig, "sum_h3_trP2_horner": S_hor, "sum_h3_trP2_script": cls["c"]["sum_h3_trP2"],
                "gamma_own": gamma_own, "gamma_script": cls["c"]["gamma"], "gamma_quoted": 4.742198e-4,
                "E_I1_own_stencils": E_I1_own, "E_I1_block_reads": br["E_curv_I1"], "E_total_own": E_tot_own, "E_total_record": 13.834938211360924,
                "far_field_M00_mean_r_gt_20": m00_far, "V4_under_s_plus_1_targets": V4_other})
    line("C1_E_V4_by_eigenvalue_power_sums_equals_e_parts_1e-10_and_the_quoted_0.3069567721158226",
         abs(V4_eig - ev) / ev < 1e-10 and abs(ev - 0.3069567721158226) < 1e-12,
         f"own {V4_eig:.16g}, e_parts {ev:.16g}, max |Im lambda| {max_imag:.1e}")
    line("C2_sum_h3_trP2_by_eigenvalues_and_by_horner_agree_1e-10_and_match_script",
         abs(S_eig - S_hor) / S_hor < 1e-10 and abs(S_eig - cls["c"]["sum_h3_trP2"]) / S_eig < 1e-10,
         f"eig {S_eig:.10f}, horner {S_hor:.10f}, script {cls['c']['sum_h3_trP2']:.10f}")
    line("C3_gamma_4.742198e-4_within_1e-6_and_matches_script_1e-10",
         abs(gamma_own - 4.742198e-4) / 4.742198e-4 < 1e-6 and abs(gamma_own - cls["c"]["gamma"]) / gamma_own < 1e-10,
         f"own gamma {gamma_own:.9e}, script {cls['c']['gamma']:.9e}")
    line("C4_record_total_13.834938_rebuilt_as_own_I1_plus_V4_1e-9_and_own_I1_equals_block_reads",
         abs(E_tot_own - 13.834938211360924) / 13.834938211360924 < 1e-9 and abs(E_I1_own - br["E_curv_I1"]) / E_I1_own < 1e-10,
         f"own E_total {E_tot_own:.12f}, own E_I1 {E_I1_own:.12f}, block_reads {br['E_curv_I1']:.12f}")
    line("C5_record_field_lives_on_the_s_minus_1_branch_M00_far_is_plus_8_and_V4_under_8_1_0.3_0_targets_is_huge",
         abs(m00_far - 8.0) < 1e-6 and V4_other > 1e3 * ev,
         f"far-field M_00 {m00_far:.8f} (N_00 = {-m00_far:.3f}); V4 with C_p from (8, 1, 0.3, 0): {V4_other:.4g} vs the record's {ev:.4g}")
    return out


# ================= claim 4 =================
def audit_4(cls):
    out = {}
    n, L = 8, 12.0
    cfg = B3.base_cfg(s=-1.0, g=G8, n=n, L=L, delta=DELTA)
    h = cfg["h"]
    X = B3.coords(n, h)
    Mv = np.broadcast_to(B3.vac4(cfg), (n, n, n, 4, 4)).copy()
    eps = np.array([0.3, 0.1, 0.03, 0.01])
    rng = np.random.default_rng(404)

    def E_of(shape):
        return np.array([EN.block_reads(Mv + e * shape, cfg, "I1")["E_curv_I1"] for e in eps])

    def expo_dev(Es):
        return float(np.max(np.abs(np.log(Es[:-1] / Es[1:]) / np.log(eps[:-1] / eps[1:]) - 4.0)))

    def ratio_spread(Es):
        r = Es / eps ** 4
        return float((r.max() - r.min()) / r.max())
    # single-axis, own random symmetric jets (12) + the 10 unit ones per axis
    rand_syms = [B3.sym4(rng.standard_normal((4, 4))) for _ in range(12)]
    units = [np.array(b, dtype=float) for b in [np.array(m.tolist(), dtype=float) for m in sym_basis_exact()]]
    single_max = 0.0
    for ax in range(3):
        for B in rand_syms + units:
            Es = E_of(X[ax][..., None, None] * B)
            single_max = max(single_max, float(np.max(np.abs(Es))))
    out["single_axis_max_abs_E"] = single_max
    line("D1_single_axis_jets_66_shapes_give_exactly_zero", single_max == 0.0, f"max |E_curv_I1| over 3 axes x 22 jets x 4 eps = {single_max}")
    # two-axis pairs: zero iff [B_a, B_b]_eta = 0
    pairs = []
    n_pred_zero, n_pred_nz, agree = 0, 0, 0
    for k in range(24):
        ax1, ax2 = rng.choice(3, size=2, replace=False)
        if k < 12:
            Ba, Bb = units[rng.integers(10)], units[rng.integers(10)]
        else:
            Ba, Bb = rand_syms[rng.integers(12)], rand_syms[rng.integers(12)]
        comm = Ba @ ETA @ Bb - Bb @ ETA @ Ba
        pred_zero = bool(np.max(np.abs(comm)) == 0.0)
        Es = E_of(X[ax1][..., None, None] * Ba + X[ax2][..., None, None] * Bb)
        # zero up to the entrants' roundoff floor: a same-matrix random pair returns E ~ 1e-29 (1e-33 of the O(1) pairs),
        # the fwd_bwd summation order not cancelling exactly where this file's direct F = 0 does
        is_zero = bool(np.max(np.abs(Es)) < 1e-20)
        row = {"axes": [int(ax1), int(ax2)], "comm_zero": pred_zero, "E_zero": is_zero, "E": Es.tolist(), "max_abs_E": float(np.max(np.abs(Es)))}
        if not is_zero:
            row["exponent_dev"] = expo_dev(Es); row["eps4_ratio_spread"] = ratio_spread(Es)
        pairs.append(row)
        n_pred_zero += pred_zero; n_pred_nz += (not pred_zero); agree += (pred_zero == is_zero)
    out["two_axis_pairs"] = pairs
    worst_pair = max([r.get("exponent_dev", 0.0) for r in pairs] + [0.0])
    floor = max([r["max_abs_E"] for r in pairs if r["E_zero"]] + [0.0])
    smallest_nz = min([r["max_abs_E"] for r in pairs if not r["E_zero"]] + [np.inf])
    out["two_axis_zero_floor_max_abs_E"] = floor
    out["two_axis_smallest_nonzero_max_abs_E"] = float(smallest_nz)
    line("D2_two_axis_pairs_zero_iff_eta_commutator_zero_24_own_pairs_both_classes_present_nonzero_ones_exactly_quartic",
         agree == 24 and n_pred_zero > 0 and n_pred_nz > 0 and worst_pair < 1e-9 and floor < 1e-20 * smallest_nz,
         f"agreement {agree}/24 (predicted zero {n_pred_zero}, nonzero {n_pred_nz}); worst |exponent - 4| {worst_pair:.1e}; roundoff floor of the zero class {floor:.1e} vs smallest nonzero {smallest_nz:.2e}")
    # the class script's 12 pairs re-classified
    cls_pairs = cls["d"]["two_axis_pairs"]
    ok_cls = 0
    for r in cls_pairs:
        a, b = r["sym_indices"]
        comm = units[a] @ ETA @ units[b] - units[b] @ ETA @ units[a]
        pred_zero = bool(np.max(np.abs(comm)) == 0.0)
        ok_cls += (pred_zero == (r["E_curv"][2] == 0.0))
    out["class_pairs_reclassified_agree"] = ok_cls
    line("D3_class_script_12_pairs_9_zero_3_quartic_match_the_commutator_criterion", ok_cls == 12 and cls["d"]["pairs_identically_zero"] == 9 and cls["d"]["pairs_quartic"] == 3,
         f"{ok_cls}/12 agree; script counts zero {cls['d']['pairs_identically_zero']}, quartic {cls['d']['pairs_quartic']}")
    # random 30-jet combination, quadratic shapes, mixed shape: all exactly quartic; own stencil agrees with block_reads
    shapes = {}
    C = rng.standard_normal((3, 10))
    S = sum(C[ax, a] * X[ax][..., None, None] * units[a] for ax in range(3) for a in range(10))
    shapes["linear_30_jets"] = S
    Q = sum(rng.standard_normal() * (X[i] * X[j])[..., None, None] * rand_syms[rng.integers(12)] for i in range(3) for j in range(i, 3))
    shapes["quadratic_in_x"] = Q
    shapes["mixed_linear_plus_quadratic_plus_cubic"] = S + Q + (X[0] * X[1] * X[2])[..., None, None] * rand_syms[0]
    shapes["smooth_gaussian_bump"] = np.exp(-(X[0] ** 2 + X[1] ** 2 + X[2] ** 2) / 8.0)[..., None, None] * rand_syms[1] + np.sin(X[2])[..., None, None] * rand_syms[2]
    rows = {}
    worst_dev, worst_spread, worst_stencil = 0.0, 0.0, 0.0
    for name, sh in shapes.items():
        Es = E_of(sh)
        own = np.array([own_E_I1(Mv + e * sh, h) for e in eps])
        rows[name] = {"E": Es.tolist(), "exponent_dev": expo_dev(Es), "eps4_ratio_spread": ratio_spread(Es), "own_stencil_rel": rel(own, Es)}
        worst_dev = max(worst_dev, rows[name]["exponent_dev"]); worst_spread = max(worst_spread, rows[name]["eps4_ratio_spread"]); worst_stencil = max(worst_stencil, rows[name]["own_stencil_rel"])
    out["shapes"] = rows
    line("D4_random_30_jet_quadratic_mixed_and_gaussian_shapes_all_exactly_quartic_E_over_eps4_constant_1e-9",
         worst_dev < 1e-9 and worst_spread < 1e-9, f"worst |exponent - 4| {worst_dev:.1e}, worst spread of E / eps^4 {worst_spread:.1e} over 4 shapes x 4 eps")
    line("D5_own_stencil_I1_equals_block_reads_on_every_shape_1e-10", worst_stencil < 1e-10, f"worst rel {worst_stencil:.1e}")
    # the qualifier: the second variation of L_cert about the flat vacuum is the mass term, not zero
    lam4 = np.linalg.eigvalsh(np.array(closed_V4([sp.Integer(-8), sp.Integer(1), sp.Rational(3, 10), sp.Integer(0)], sp.nsimplify(B3.W1)), dtype=float))
    out["V4_mass_term_eigs_at_W1"] = lam4.tolist()
    line("D6_qualifier_second_variation_of_L_cert_is_not_zero_the_V4_mass_term_is_positive_on_four_directions", bool(np.all(lam4 > 0)),
         f"H_V4 eigenvalues at W1: {[f'{v:.3g}' for v in lam4]}; the class statement 'zero on the 30-jet space' holds for the derivative sector only")
    return out


# ================= claim 5 =================
def _sections(text):
    """(heading, body) list split at ## / ### headings."""
    parts = re.split(r"^(#{2,3} .*)$", text, flags=re.M)
    secs = []
    for i in range(1, len(parts), 2):
        secs.append((parts[i].strip(), parts[i + 1]))
    return secs


def audit_5(cls):
    out = {}
    text = open(LEDGER, encoding="utf-8").read()
    secs = _sections(text)
    def where(sub):
        return [h[:60] for h, b in secs if sub in b]
    strings = {"string_law": "-18.66 + 0.935 d", "r3_3869": "3869", "r3_2784": "2784", "one_over_d_fit": "fits `1 / d`",
               "best_exponent_1": "best exponent 1", "coulomb_identity": "Coulomb identity", "object_independent": "object-independent"}
    found = {k: where(v) for k, v in strings.items()}
    out["ledger_hits_by_section"] = found
    sec3 = [b for h, b in secs if h.startswith("## 3.")][0]
    sec68 = [b for h, b in secs if h.startswith("### 6.8")][0]
    sec69 = [b for h, b in secs if h.startswith("### 6.9")][0]
    line("E1_string_law_-18.66_+_0.935_d_appears_in_ledger_6.8_and_6.9", "-18.66 + 0.935 d" in sec68 and "-18.66 + 0.935 d" in sec69, f"sections holding it: {found['string_law']}")
    in3 = ("3869" in sec3) and ("2784" in sec3)
    rows_e = cls["e"]["rows"]
    line("E2_R3_numbers_3869_and_2784_are_in_ledger_6.8_and_6.9_not_in_section_3_and_the_class_row_cites_6.8",
         ("3869" in sec68 and "2784" in sec68 and "3869" in sec69 and "2784" in sec69) and not in3 and "6.8" in rows_e[1]["source"] and "section 3)" not in rows_e[1]["source"],
         f"section 3 heading: '{[h for h, b in secs if h.startswith('## 3.')][0][:50]}', holds 3869 {'3869' in sec3}, 2784 {'2784' in sec3}; sections holding 3869: {found['r3_3869']}")
    # the six-point law refit from the task record's R19-3 table
    d6 = np.array([6.0, 9.0, 12.0, 15.0, 18.0, 24.0]); e6 = np.array([-12.59, -10.53, -7.61, -5.54, -0.86, 3.67])
    b, a = np.polyfit(d6, e6, 1)
    rms = float(np.sqrt(np.mean((e6 - (a + b * d6)) ** 2)))
    out["six_point_refit"] = {"d": d6.tolist(), "E_int": e6.tolist(), "intercept": float(a), "slope": float(b), "rms": rms}
    line("E3_six_point_law_refits_to_-18.66_+_0.935_d_rms_0.59_from_the_task_record_points", abs(a + 18.66) < 0.03 and abs(b - 0.935) < 0.003 and abs(rms - 0.59) < 0.03,
         f"refit E_int = {a:.2f} + {b:.3f} d, rms {rms:.2f}")
    d_range_ok = ("d 6, 9, 12, 15, 18, 24" in rows_e[0]["law"]) and ("d 12 to 24, no bend" not in rows_e[0]["law"])
    line("E4_the_class_row_states_the_six_points_as_d_6_9_12_15_18_24", d_range_ok,
         f"the six points are d {d6.astype(int).tolist()} (R18-3's four at d 12 to 24 plus R19-3's d 6 and 9); class row law: '{rows_e[0]['law']}' (the first version said 'six points, d 12 to 24', refuted 2026-09-13 and fixed)")
    # the R3 lambda = 0 rows and the R19-1 numbers from the R19-1 JSON
    r19 = json.load(open(R19_1_JSON))
    r3 = r19["results"]["R3_lam0_dr_n32_g32"]["rows"]
    i1dr = r19["results"]["I1_dr_n32_g32"]["rows"]
    out["R3_lam0_dressed_g32"] = {k: v["E_int"] for k, v in r3.items()}
    line("E5_R3_lambda0_like_pair_g32_E_int_+3869_d10_and_+2784_d24_and_monotone_fall", abs(r3["d10"]["E_int"] - 3869) < 1 and abs(r3["d24"]["E_int"] - 2784) < 1
         and r3["d10"]["E_int"] > r3["d14"]["E_int"] > r3["d24"]["E_int"] and abs(i1dr["d10"]["E_int"] - r3["d10"]["E_int"]) < 1e-3,
         f"R3 lambda = 0 rows: {[(k, round(v['E_int'], 1)) for k, v in sorted(r3.items())]}; R19-1's own I1 dressed rows reproduce them")
    gam = r19["results"]["Gam_dr_n32_g32_ext2"]["rows"]
    dd = np.array([10.0, 14.0, 18.0, 24.0]); ee = np.array([gam[f"d{int(x)}"]["E_int"] for x in dd])
    def r2(basis):
        A = np.column_stack([np.ones_like(dd), basis]); coef, *_ = np.linalg.lstsq(A, ee, rcond=None)
        res = ee - A @ coef; return coef, 1.0 - np.sum(res ** 2) / np.sum((ee - ee.mean()) ** 2)
    (c1, r2_1), (c5, r2_5), (cl, r2_l) = r2(1.0 / dd), r2(1.0 / dd ** 5), r2(np.log(dd) / dd)
    out["Gam_g32_ext2_tail_fits"] = {"E_int": ee.tolist(), "A_plus_B_over_d": {"A": float(c1[0]), "B": float(c1[1]), "R2": float(r2_1)},
                                     "A_plus_B_over_d5": {"R2": float(r2_5)}, "A_plus_B_lnd_over_d": {"R2": float(r2_l)}}
    line("E6_R19-1_Gam_g32_tail_at_10500_steps_fits_1_over_d_with_R2_0.993_positive_B_and_beats_1_over_d5", abs(r2_1 - 0.993) < 0.002 and c1[1] > 0 and r2_1 > r2_5,
         f"A + B/d: B {c1[1]:.0f}, R^2 {r2_1:.4f}; 1/d^5 R^2 {r2_5:.3f}; ln d / d R^2 {r2_l:.4f}; E_int falls {ee[0]:.0f} -> {ee[-1]:.0f} (repulsive force)")
    un = r19["results"]["I1_un_n32_g8"]
    eu = [un["rows"][f"d{k}"]["E_int"] for k in (12, 18, 24, 30)]
    out["I1_undressed_g8"] = {"E_int_d12_18_24_30": eu, "record_outcome": un.get("outcome")}
    rising = all(eu[i + 1] > eu[i] for i in range(3))
    line("E7_undressed_I1_pair_g8_E_int_RISES_with_d_and_the_class_row_says_so", rising and "RISING" in rows_e[3]["law"] and rows_e[3]["form"] != "repulsive",
         f"E_int at d 12 / 18 / 24 / 30 = {[round(v, 2) for v in eu]}: positive but RISING with d (dE_int/dd > 0, a string); the record's outcome string: '{un.get('outcome')}'; class row form: '{rows_e[3]['form']}' (the first version said 'repulsive', refuted 2026-09-13 and fixed)")
    return out


# ================= claim 6 =================
def audit_6(cls, a1):
    out = {}
    text = open(LEDGER, encoding="utf-8").read()
    sec69 = [b for h, b in _sections(text) if h.startswith("### 6.9")][0]
    row1 = "M_00 = -g" in sec69 and "(8, 1, 0.3, 0)" in sec69 and "share the convention" in sec69
    out["row1_text_found"] = row1
    code_q = a1["stack_s_minus_1"]["q_code"]
    line("F1_ledger_6.9_row1_says_M00_=_-g_and_the_stack_spectrum_is_8_1_0.3_0_which_the_s_minus_1_code_branch_CONTRADICTS",
         row1 and code_q == [-8.0, 1.0, 0.3, 0.0] and a1["stack_s_minus_1"]["M_vac_code"] == [8.0, 1.0, 0.3, 0.0],
         f"row 1 text present {row1}; code s = -1: vac4 {a1['stack_s_minus_1']['M_vac_code']} (M_00 = +8), N spectrum {code_q}; (8, 1, 0.3, 0) is the s = +1 branch's spectrum (vac4 {a1['report_s_plus_1']['M_vac_code']})")
    sec11 = [b for h, b in _sections(text) if h.startswith("### 1.1")][0]
    wrong_form = "M_vac = diag(-g, 1, delta, 0) (s = -1" in sec11
    corrected = "diag(-sg, 1, delta, 0) = diag(+g, 1, delta, 0)" in sec11 and "spectrum (-g, 1, delta, 0)" in sec11 and "C_p = (s g)^p" in sec11
    out["sec_1_1_vacuum_line"] = {"wrong_form_present": wrong_form, "corrected_form_present": corrected}
    line("F1b_ledger_1.1_vacuum_line_carries_M_00_=_+g_on_s_=_-1_consistent_with_its_own_C_p_=_(sg)^p_and_vac4",
         corrected and not wrong_form,
         f"corrected form present {corrected}, the first version's 'M_vac = diag(-g, 1, delta, 0) (s = -1 ...)' present {wrong_form} (it contradicted C_p = (s g)^p + 1 + delta^p and vac4, which put N_00 = s g = -8, i.e. M_00 = +8 on s = -1; refuted 2026-09-13 and corrected)")
    quoted = np.array([0.059, 1.50, 12.6, 4.23e6])
    eig_plus = np.array(a1["report_s_plus_1"]["V4_diag_block_eigs_own"])
    eig_minus = np.array(a1["stack_s_minus_1"]["V4_diag_block_eigs_own"])
    rp, rm = eig_plus / quoted, eig_minus / quoted
    out["row2_ratio_true_over_quoted"] = {"(8,1,0.3,0)": rp.tolist(), "(-8,1,0.3,0)": rm.tolist()}
    row2 = "0.059 / 1.50 / 12.6 / 4.23e6" in sec69
    line("F2_ledger_6.9_row2_V4_eigenvalues_are_HALF_the_exact_2JtJ_eigenvalues_on_8_1_0.3_0_and_no_constant_factor_on_-8_1_0.3_0",
         row2 and np.all(np.abs(rp - 2.0) < 0.02) and (rm.max() - rm.min()) > 0.3,
         f"row 2 text present {row2}; true / quoted on (8, 1, 0.3, 0) = {[f'{v:.3f}' for v in rp]}; on (-8, 1, 0.3, 0) = {[f'{v:.2f}' for v in rm]}; exact on (8, 1, 0.3, 0): {[f'{v:.4g}' for v in eig_plus]}, on (-8, 1, 0.3, 0): {[f'{v:.4g}' for v in eig_minus]}")
    # the stiffness ratios quoted in row 2 (7e7 for V4, 7e4 for V_spec) on each root set
    out["stiffness_ratios"] = {"V4_(8)": float(eig_plus[-1] / eig_plus[0]), "V4_(-8)": float(eig_minus[-1] / eig_minus[0]),
                               "Vspec_(8)": float(max(a1["report_s_plus_1"]["Vspec_diag_own"]) / min(a1["report_s_plus_1"]["Vspec_diag_own"])),
                               "Vspec_(-8)": float(max(a1["stack_s_minus_1"]["Vspec_diag_own"]) / min(a1["stack_s_minus_1"]["Vspec_diag_own"]))}
    sr = out["stiffness_ratios"]
    line("F3_row2_stiffness_ratios_7e7_and_7e4_hold_on_8_1_0.3_0_and_differ_on_the_stack_roots",
         abs(np.log10(sr["V4_(8)"]) - np.log10(7e7)) < 0.1 and abs(np.log10(sr["Vspec_(8)"]) - np.log10(7e4)) < 0.1 and sr["V4_(-8)"] < 0.9 * sr["V4_(8)"],
         f"V4 ratio on (8) {sr['V4_(8)']:.2e}, on (-8) {sr['V4_(-8)']:.2e}; V_spec ratio on (8) {sr['Vspec_(8)']:.2e}, on (-8) {sr['Vspec_(-8)']:.2e}")
    return out


# ================= main =================
def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, float) and not np.isfinite(o):
        return str(o)
    return o


def main():
    cls = json.load(open(CLASS_JSON))
    out = {"task": "M5.32 R20-0 audit", "class_json_utc": cls.get("utc"), "convention_under_audit": cls.get("convention")}
    log("claim 1: exact Hessians")
    out["claim1_hessians"] = audit_1(cls)
    log("claim 2: exact gradient")
    out["claim2_gradient"] = audit_2(cls)
    log("claim 3: gamma")
    out["claim3_gamma"] = audit_3(cls)
    log("claim 4: flat-vacuum jets")
    out["claim4_jets"] = audit_4(cls)
    log("claim 5: pair laws on record")
    out["claim5_pair_laws"] = audit_5(cls)
    log("claim 6: ledger 6.9 rows 1 and 2")
    out["claim6_ledger_6_9"] = audit_6(cls, out["claim1_hessians"])
    out["lines"] = LINES
    out["runtime_s"] = time.time() - T0
    with open(OUT, "w") as f:
        json.dump(_jsonable(out), f, indent=1)
    npass = sum(1 for v in LINES.values() if v["pass"])
    print(f"{npass}/{len(LINES)} PASS, runtime {out['runtime_s']:.1f}s, wrote {OUT}")
    return out


if __name__ == "__main__":
    main()
