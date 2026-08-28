"""M5.32 R7 arm (a): the class-C4 two-derivative term K_T at FORM level.

Reuses (imports, never modifies) the R0 registry m5_32_lagrangian.py, the
R1 extended registry m5_32_terms_ext.py (h_cov_np, timelike_eig_sym,
_transform), the R1 audit channels (m5_32_r1_audit_symbolic.py), the R2
form-level machinery m5_32_r2_a_formmap.py (all_channels, jets_batch,
Q_forms, min_eig) and m5_32_r2_audit_formmap.py (channel_set, boost_to_
frame-free h_of), the R4 audit family m5_32_r4_audit_clock.py (LocFamily,
window) and its n = 64, L = 96 box ladder (the static gain G(R)).

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1,1,1,1), jets A_mu = d_mu M
(raw entries CONTRAVARIANT, the R0 rule), A_0 = omega a0 on a static
background. u(M) = the timelike unit eigenvector of N = M eta (u^T eta u
= -1), h_cov = eta + 2 (eta u)(eta u)^T (the R1 flip metric, covariant).

THE TERM (class C4, the lower-order kinetic term with the time row
isolated by the field-dependent metric):
    K_T = 1/2 sum_mu eta^{mu mu} [ tr(h_cov A_mu h_cov A_mu)
                                   - tr(eta A_mu eta A_mu) ]            (1)
u-frame (u = e0, h_cov = 1): tr(A A) - tr(eta A eta A) = sum_ab (1 - eta_a
eta_b) A_ab A_ba = 4 sum_j A_{0j} A_{j0} = 4 sum_j A_{0j}^2 for SYMMETRIC
jets (every physical jet d_mu M is symmetric; on an ANTISYMMETRIC probe
a0 the sign flips: the R1/R2 clock channels coms(gamma, M) with a boost
part and the x_boost probes are such probes, flagged per channel), so
    K_T = 2 sum_mu eta^{mu mu} sum_j (A_mu)_{0j}^2                      (2)
       = 2 sum_i sum_j (A_i)_{0j}^2  -  2 omega^2 sum_j (a0)_{0j}^2
       (time-row entries of the jets in the frame where u = e0);
K_T = 0 iff the u-frame time row of every A_mu vanishes. Both traces are
Lorentz scalars (h_cov, eta covariant internal metrics, the eta^{mu mu}
weight on the derivative pair), so K_T is SO(1,3)-invariant by
construction; dropping the eta^{mu mu} weights or replacing h_cov by the
Frobenius identity (the M5.21.16 variant-A read) breaks it.

Candidate action and Legendre read (the R2 formmap convention):
    L = -4 [(1 - lambda) I1 + lambda I1_h] - c2 K_T - V4,   c2 > 0.
For a term I(omega) = A + B omega + C omega^2 its energy is
H_I = C omega^2 - A; the action's energy is sum_k c_k (C_k omega^2 - A_k).
For K_T (eq. 2): A_KT = 2 sum_i sum_j (A_i)_{0j}^2 >= 0 (spatial part),
C_KT = -2 sum_j (a0)_{0j}^2 <= 0 (time part, eta^00 = -1), B_KT = 0, so
    E[-c2 K_T] = c2 [ 2 sum_i sum_j (A_i)_{0j}^2 + 2 omega^2 sum_j (a0)_{0j}^2 ]
which is >= 0 pointwise for c2 > 0: the omega^2 coefficient of the ENERGY
(the R2 "H2") contributed by -c2 K_T is K0 := 2 c2 sum_j (a0)_{0j}^2, a
CONSTANT in the spatial jet (like Pgrad), not a form. A constant cannot
make an indefinite homogeneous form PSD (scale the jet x -> s x), so c2
does not replace lambda >= 1/2; it never hurts (K0 >= 0).

Derrick: K_T carries two derivatives, so under the spatial dilation
b(r) -> b(r / mu) of a boost dressing its energy scales as mu^{3-2} =
mu^1 (the R6 table), i.e. E_KT ∝ R at fixed amp for a wide dressing on
the hedgehog background (the dressing enters through the boost-generator
gradient ~ amp / R times the background entries, squared, times R^3).
Static block-diagonal fields (uniform time row, u = e0): K_T = 0 exactly,
so the certified 3x3 hedgehog, the 3-lepton census and the Coulomb
anchors are untouched by construction.

Localization estimate: with the R4 audit's envelope static gain G(R) =
min_amp E_stat(amp, R) - E_stat(0) at lambda = 0.75, 1 (n = 64, L = 96)
and E_KT(R) ∝ R on the same family, the total c2 E_KT(R) + G(R) has an
interior minimum R*(c2) iff the gain grows slower than R^1 (power law
G = -a R^s, s < 1: R* = (a s / (c2 k))^(1/(1-s)), k = E_KT / R).

Out: ../data/m5_32_r7_kt_form.json
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time

import numpy as np
import sympy as sp
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r7_kt_form.json")
T0 = time.time()

sys.argv = [sys.argv[0]]                       # the imported mains parse argv

R2F = __import__("importlib.util").util
def _load(name, fname):
    spec = R2F.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = R2F.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


FM = _load("m5_32_r2_a_formmap", "m5_32_r2_a_formmap.py")      # channels, jets, Q_forms
L0, EXT, AUD = FM.L0, FM.EXT, FM.AUD
R2A = _load("m5_32_r2_audit_formmap", "m5_32_r2_audit_formmap.py")
R4 = _load("m5_32_r4_audit_clock", "m5_32_r4_audit_clock.py")
R2L, B3, B8 = R4.R2, R4.B3, R4.B8
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
XI = sp.diag(-1, 1, 1, 1)
G, DELTA, S = 32.0, 0.3, -1.0
P = L0.default_params(s=S, g=G, delta=DELTA)
RES = {"arm": "R7.a form-level K_T (class C4)",
       "reused": ["m5_32_lagrangian.py", "m5_32_terms_ext.py", "m5_32_r1_audit_symbolic.py",
                  "m5_32_r2_a_formmap.py", "m5_32_r2_audit_formmap.py", "m5_32_r4_audit_clock.py",
                  "m5_32_r4_audit_clock_box_n64_L96.json"],
       "conventions": {
           "eta": "diag(-1,1,1,1); jets A_mu = d_mu M raw contravariant; u = timelike unit eigenvector of M eta",
           "K_T": "1/2 sum_mu eta^{mu mu} [tr(h A_mu h A_mu) - tr(eta A_mu eta A_mu)], h = eta + 2 (eta u)(eta u)^T",
           "action": "L = -4[(1-lambda) I1 + lambda I1_h] - c2 K_T - V4, c2 > 0",
           "Legendre": "I(omega) = A + B omega + C omega^2 -> H_I = C omega^2 - A; energy of the action = sum_k c_k (C_k omega^2 - A_k); "
                       "H2 = the omega^2 coefficient of the energy = sum_k c_k C_k (the R2 formmap convention, c_I1 = -4)",
           "K_T energy": "E[-c2 K_T] = c2 (2 sum_i sum_j (A_i)_{0j}^2 + 2 omega^2 sum_j (a0)_{0j}^2) in the u-frame; K0 = 2 c2 sum_j (a0)_{0j}^2 = its H2 contribution (a constant)",
           "toy point": {"g": G, "delta": DELTA, "s": S}},
       "claims": {}}


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


# ================= the companion registry entry =================
def kt_density_np(A, M, p, h=None, weights=(-1.0, 1.0, 1.0, 1.0), h_mode="cov"):
    """eq. (1) per cell. A (4, ..., 4, 4), M (..., 4, 4). h_mode: 'cov'
    (registry h_cov), 'delta' (mutant: Frobenius identity), 'eta' (control:
    identically zero)."""
    if h is None:
        if h_mode == "cov":
            h = EXT.h_cov_np(M)
        elif h_mode == "delta":
            h = np.broadcast_to(np.eye(4), M.shape)
        else:
            h = np.broadcast_to(ETA, M.shape)
    tot = 0.0
    for mu in range(4):
        Am = A[mu]
        th = np.einsum("...ab,...bc,...cd,...da->...", h, Am, h, Am)
        te = np.einsum("ab,...bc,cd,...da->...", ETA, Am, ETA, Am)
        tot = tot + 0.5 * weights[mu] * (th - te)
    return tot


def kt_density_sym(F, M, p):
    """eq. (1) exact: jets F['jets'] (4 sympy matrices), M rational."""
    A = F["jets"]
    u0 = EXT.timelike_eig_sym(M)[0]
    hu = XI * u0
    h = XI + 2 * hu * hu.T
    tot = 0
    for mu in range(4):
        tot += sp.Rational(1, 2) * XI[mu, mu] * (
            sp.trace(h * A[mu] * h * A[mu]) - sp.trace(XI * A[mu] * XI * A[mu]))
    return sp.simplify(tot)


class KTTerm(EXT.TermX):
    """registry shape of m5_32_terms_ext.TermX (duck-typed for
    m5_32_lagrangian.term_lagrangian / omega_decompose)."""

    def __init__(self):
        super().__init__(
            "K_T",
            "K_T = 1/2 sum_mu eta^{mu mu} [tr(h_cov A_mu h_cov A_mu) - "
            "tr(eta A_mu eta A_mu)], A_mu = d_mu M (raw contravariant entries), "
            "h_cov = eta + 2 (eta u)(eta u)^T, u = timelike unit eigenvector of "
            "M eta (u^T eta u = -1); u-frame: 2 sum_mu eta^{mu mu} sum_j (A_mu)_{0j}^2 "
            "(class C4, two derivatives, the time-row kinetic term)",
            kt_density_sym, np_fn=kt_density_np, kind="kinetic_C4")
        self.definition_hash = hashlib.sha256(self.definition.encode()).hexdigest()


KT = KTTerm()
REGISTRY_R7 = {"K_T": KT}


# ================= 1. covariance + mutants + sympy =================
def stage_covariance(rng):
    A, M = L0._random_jets(rng, 48, P)
    out = {"n_jets": 48, "per_transform": {}, "mutants": {}}
    worst = 0.0
    d0 = KT.density(A, M, P)
    sc = max(np.max(np.abs(d0)), 1e-300)
    for kind in ("boost", "rotation"):
        for k in range(4):
            Lm = L0._lorentz(rng, kind, scale=0.5)
            Ap, Mp = EXT._transform(Lm, A, M)
            d1 = KT.density(Ap, Mp, P)
            dr = float(np.max(np.abs(d1 - d0)) / sc)
            out["per_transform"][f"{kind}_{k}"] = dr
            worst = max(worst, dr)
    # a combined SO(1,3) x SO(3) element, larger rapidity
    Lm = L0._lorentz(rng, "boost", 1.0) @ L0._lorentz(rng, "rotation", 1.0)
    Ap, Mp = EXT._transform(Lm, A, M)
    dr = float(np.max(np.abs(KT.density(Ap, Mp, P) - d0)) / sc)
    out["per_transform"]["boost1.0_x_rot1.0"] = dr
    worst = max(worst, dr)
    out["worst_drift"] = worst
    # mutants (must fail): drop the eta^{mu mu} weights; h -> delta
    for nm, kw in (("no_weights", {"weights": (1.0, 1.0, 1.0, 1.0)}),
                   ("h_to_delta", {"h_mode": "delta"})):
        m0 = kt_density_np(A, M, P, **kw)
        m1 = kt_density_np(Ap, Mp, P, **kw)
        out["mutants"][nm] = float(np.max(np.abs(m1 - m0)) / max(np.max(np.abs(m0)), 1e-300))
    out["control_h_to_eta_is_zero"] = float(np.max(np.abs(kt_density_np(A, M, P, h_mode="eta"))))
    # sympy vs numpy on rational jets at a rationally boosted vacuum
    Lr = EXT.rational_lorentz()
    rows = []
    for trial in range(2):
        def rsym():
            X = sp.zeros(4, 4)
            for i in range(4):
                for j in range(i, 4):
                    X[i, j] = X[j, i] = sp.Rational(int(rng.integers(-9, 10)), int(rng.integers(1, 7)))
            return X
        Mmu = [rsym() for _ in range(4)]
        D = sp.diag(32, 1, sp.Rational(3, 10), 0)
        Mpt = Lr * D * Lr.T
        vs = float(sp.nsimplify(KT.sympy({"jets": Mmu}, Mpt, P)))
        A_np = np.array([[[float(Mmu[m][i, j]) for j in range(4)] for i in range(4)] for m in range(4)])[:, None]
        M_np = np.array([[float(Mpt[i, j]) for j in range(4)] for i in range(4)])[None]
        vn = float(KT.density(A_np, M_np, P)[0])
        # u-frame check: pull back by Lr^-1 and read the time rows (eq. 2)
        Li = np.linalg.inv(np.array(Lr.tolist(), dtype=float))
        Ab, Mb = EXT._transform(Li, A_np, M_np)
        assert np.allclose(Mb[0], np.diag([32, 1, 0.3, 0]), atol=1e-12)
        uf = float(sum(np.diag(ETA)[mu] * 2.0 * np.sum(Ab[mu, 0, 0, 1:] ** 2) for mu in range(4)))
        rows.append({"sympy": vs, "numpy": vn, "u_frame_eq2": uf,
                     "rel_sympy_numpy": L0._rel(vn, vs), "rel_eq2": L0._rel(uf, vs)})
    out["sympy_vs_numpy"] = rows
    out["worst_sympy_numpy"] = max(r["rel_sympy_numpy"] for r in rows)
    out["worst_eq2"] = max(r["rel_eq2"] for r in rows)
    # EXACT covariance at large rapidity (the float drift of the combined transform above is a
    # cancellation floor: intermediate traces ~ 1e5-1e6 cancel to a density ~ 10): rational jets,
    # rational vacuum, the rational boost Lr^4 (rapidity 4 arccosh(5/4) = 2.77) composed with the
    # rational rotation; K_T(transformed) - K_T(original) must be exactly 0 in sympy
    L4 = Lr ** 4
    assert L4.T * XI * L4 == XI
    Mmu = [rsym() for _ in range(4)]
    D = sp.diag(32, 1, sp.Rational(3, 10), 0)
    Lm4 = sp.Matrix(L4)
    L4invT = (Lm4.inv()).T
    Mmu_t = [sum((L4invT[mu, nu] * (Lm4 * Mmu[nu] * Lm4.T) for nu in range(4)), sp.zeros(4, 4)) for mu in range(4)]
    v0 = sp.nsimplify(KT.sympy({"jets": Mmu}, D, P))
    v1 = sp.nsimplify(KT.sympy({"jets": Mmu_t}, Lm4 * D * Lm4.T, P))
    A_np = np.array([[[float(Mmu[m][i, j]) for j in range(4)] for i in range(4)] for m in range(4)])[:, None]
    M_np = np.array([[float(D[i, j]) for j in range(4)] for i in range(4)])[None]
    L4n = np.array(L4.tolist(), dtype=float)
    Ap4, Mp4 = EXT._transform(L4n, A_np, M_np)
    h4 = EXT.h_cov_np(Mp4)
    inter = float(np.max(np.abs(np.einsum("...ab,...bc,...cd,...da->...", h4, Ap4[0], h4, Ap4[0]))))
    out["exact_large_rapidity"] = {"rapidity": float(4 * np.arccosh(1.25)), "K_T_original": float(v0),
                                   "K_T_transformed_minus_original_exact": str(sp.simplify(v1 - v0)),
                                   "float_drift_same_transform": float(abs(KT.density(Ap4, Mp4, P)[0] - float(v0)) / abs(float(v0))),
                                   "max_transformed_field_entry": float(np.max(np.abs(Mp4))),
                                   "intermediate_trace_scale": inter,
                                   "note": "float64 cancellation floor ~ 1e-16 x intermediate / density"}
    out["worst_drift_standard_transforms"] = max(v for k, v in out["per_transform"].items() if not k.startswith("boost1.0"))
    out["combined_rapidity_transform_drift"] = out["per_transform"]["boost1.0_x_rot1.0"]
    ok = out["worst_drift_standard_transforms"] <= 1e-10 and sp.simplify(v1 - v0) == 0 \
        and out["worst_sympy_numpy"] <= 1e-12 and out["worst_eq2"] <= 1e-12 \
        and min(out["mutants"].values()) > 1e-3
    out["pass"] = bool(ok)
    log(f"covariance: standard drift {out['worst_drift_standard_transforms']:.2e}, combined rapidity-1 {worst:.2e}, exact large-rapidity residual {out['exact_large_rapidity']['K_T_transformed_minus_original_exact']} (float {out['exact_large_rapidity']['float_drift_same_transform']:.1e}, intermediates {out['exact_large_rapidity']['intermediate_trace_scale']:.2g}); mutants {out['mutants']}; sympy {out['worst_sympy_numpy']:.2e}; eq2 {out['worst_eq2']:.2e}")
    return out


# ================= 2. static identity + the R4 family =================
def random_block_diag_field(rng, n, L, cfg):
    """M(x) = diag(g, M3(x)) with a smooth random 3x3 symmetric M3(x)
    (uniform time row: u = e0 everywhere; the M eta spectrum stays real
    with the timelike eigenvector e0 as long as g dominates)."""
    X, Y, Z = B3.coords(n, cfg["h"])
    M = np.zeros((n, n, n, 4, 4))
    M[..., 0, 0] = -cfg["sg"]
    base = np.diag([1.0, DELTA, 0.0])
    for a in range(3):
        for b in range(a, 3):
            f = 0.0
            for _ in range(3):
                k = rng.normal(size=3) * 2 * np.pi / L * 2
                f = f + rng.normal() * 0.3 * np.cos(k[0] * X + k[1] * Y + k[2] * Z + rng.uniform(0, 2 * np.pi))
            M[..., 1 + a, 1 + b] = base[a, b] + f
            M[..., 1 + b, 1 + a] = M[..., 1 + a, 1 + b]
    return M


def kt_lattice(M, cfg, a0=None, omega=0.0, h=None):
    """h^3 sum of the K_T density over the certified stencil (the
    Lagrangian read; = A_KT for a0 = None). Optional analytic h."""
    tot = 0.0
    for A, wt in L0.lattice_jets(M, cfg, a0, omega):
        tot = tot + wt * np.sum(kt_density_np(A, M, P, h=h))
    return float(cfg["h"] ** 3 * tot)


def stage_static(rng):
    n, L = 16, 48.0
    cfg = B3.base_cfg(s=S, g=G, delta=DELTA, n=n, L=L)
    vals = []
    for k in range(20):
        M = random_block_diag_field(rng, n, L, cfg)
        # reference scale: the same field's I1 density magnitude
        e_u, _ = B3.e_parts(M, cfg)
        vals.append({"E_KT": kt_lattice(M, cfg), "E_u_scale": float(e_u)})
    worst = max(abs(v["E_KT"]) / max(v["E_u_scale"], 1e-300) for v in vals)
    out = {"n_fields": 20, "n": n, "L": L, "fields": vals, "worst_rel_to_E_u": worst,
           "pass": bool(worst <= 1e-14)}
    # a spatially varying time row must NOT vanish (the identity is not trivial)
    M = random_block_diag_field(rng, n, L, cfg)
    X, Y, Z = B3.coords(n, cfg["h"])
    M[..., 0, 1] = M[..., 1, 0] = 0.5 * np.sin(2 * np.pi * X / L)
    out["control_time_row_field_E_KT"] = kt_lattice(M, cfg)
    log(f"static identity: worst |E_KT| / E_u = {worst:.2e} (20 block-diagonal fields); control with a time row: {out['control_time_row_field_E_KT']:.4g}")
    return out


def stage_family():
    fam = R4.LocFamily(32, 48.0)
    cfg, h3 = fam.cfg, fam.cfg["h"] ** 3
    out = {"n": 32, "L": 48.0, "kind": "exp2", "amps": [0.005, 0.01, 0.02, 0.04], "R_grid": [6.0, 9.0, 12.0, 18.0, 24.0]}
    # hedgehog (amp = 0): K_T = 0
    out["E_KT_hedgehog_amp0"] = kt_lattice(fam.Mb, cfg)
    # the clock direction of the hedgehog has no time row: K0 = 0
    out["hedgehog_a0_time_row_max"] = float(np.max(np.abs(fam.a0[..., 0, 1:])))
    out["hedgehog_C_KT_omega2"] = L0.omega_decompose(KT, fam.Mb, cfg, P, fam.a0)[2]
    tab = {}
    for R in out["R_grid"]:
        for amp in out["amps"]:
            Qb = R2L.qb_field(cfg, fam.b_of(amp, R, "exp2"))
            Md = B3.sym4(R2L.conj(Qb, fam.Mb))
            a0d = B3.sym4(R2L.conj(Qb, fam.a0))
            hQ = R2L.h_from_Q(Qb)
            e_reg = kt_lattice(Md, cfg)                       # registry h (eigen-solved u)
            e_ana = kt_lattice(Md, cfg, h=hQ)                 # analytic h = Qb^-T Qb^-1
            # eq. (2) direct: pull the jets back to the u-frame with Qb^-1
            Qi = np.linalg.inv(Qb)
            e2 = 0.0
            for A, wt in L0.lattice_jets(Md, cfg):
                for mu in range(1, 4):
                    Au = np.einsum("...ab,...bc,...dc->...ad", Qi, A[mu], Qi)
                    e2 += wt * 2.0 * np.sum(Au[..., 0, 1:] ** 2)
            e2 *= h3
            C = L0.omega_decompose(KT, Md, cfg, P, a0d)[2]    # omega^2 coefficient (Lagrangian)
            a0u = np.einsum("...ab,...bc,...dc->...ad", Qi, a0d, Qi)
            K0_direct = -2.0 * h3 * float(np.sum(a0u[..., 0, 1:] ** 2))
            tab[f"R_{R:g}_amp_{amp:g}"] = {"R": R, "amp": amp, "E_KT": e_reg, "E_KT_analytic_h": e_ana,
                                          "E_KT_eq2_uframe": e2, "C_KT_omega2": C, "C_KT_eq2": K0_direct}
        log(f"family R={R:g}: E_KT(amp) " + " ".join(f"{tab[f'R_{R:g}_amp_{a:g}']['E_KT']:.4g}" for a in out["amps"]))
    out["table"] = tab
    out["worst_registry_vs_analytic_h"] = max(abs(v["E_KT"] - v["E_KT_analytic_h"]) / abs(v["E_KT"]) for v in tab.values())
    out["worst_vs_eq2"] = max(abs(v["E_KT"] - v["E_KT_eq2_uframe"]) / abs(v["E_KT"]) for v in tab.values())
    out["worst_abs_C_over_E_KT"] = max(abs(v["C_KT_omega2"]) / v["E_KT"] for v in tab.values())
    out["all_positive"] = bool(all(v["E_KT"] > 0 for v in tab.values()))
    out["all_C_nonpositive"] = bool(all(v["C_KT_omega2"] <= 1e-10 * v["E_KT"] for v in tab.values()))
    out["C_KT_zero_on_family"] = bool(out["worst_abs_C_over_E_KT"] <= 1e-10)
    out["finding_C_zero"] = ("the omega^2 coefficient of K_T VANISHES on the whole R4 boost-dressing family: the dressed clock "
                             "a0 = Qb a0_unit Qb^T pulled back to the u-frame (u = Qb e0) is the block-diagonal rotation flow, "
                             "no time row; so -c2 K_T adds NO omega^2 floor and no kinetic inertia there, only the static "
                             "localizer c2 E_KT(amp, R) > 0")
    # amp^2 law at R = 12
    e = np.array([tab[f"R_12_amp_{a:g}"]["E_KT"] for a in out["amps"]])
    out["amp_exponent_R12"] = float(np.polyfit(np.log(out["amps"]), np.log(e), 1)[0])
    # R law at amp = 0.02 (fit on R <= 18 where the box does not cut the window)
    Rs = np.array(out["R_grid"]); eR = np.array([tab[f"R_{R:g}_amp_0.02"]["E_KT"] for R in Rs])
    out["R_exponent_amp0.02_R6to18"] = float(np.polyfit(np.log(Rs[:4]), np.log(eR[:4]), 1)[0])
    out["R_exponent_amp0.02_R9to24"] = float(np.polyfit(np.log(Rs[1:]), np.log(eR[1:]), 1)[0])
    # Derrick dilation at (amp, R) = (0.02, 6): b(r) -> b(r / mu)
    mus = [0.5, 0.7, 1.0, 1.4, 2.0]
    ed = []
    for mu in mus:
        Qb = R2L.qb_field(cfg, fam.b_of(0.02, 6.0, "exp2", mu=mu))
        ed.append(kt_lattice(B3.sym4(R2L.conj(Qb, fam.Mb)), cfg))
    out["derrick"] = {"amp": 0.02, "R": 6.0, "mu": mus, "E_KT": ed,
                      "exponent_loglog": float(np.polyfit(np.log(mus), np.log(ed), 1)[0]),
                      "expected": 1.0}
    log(f"family: amp exponent {out['amp_exponent_R12']:.3f}; R exponent {out['R_exponent_amp0.02_R6to18']:.3f} (6-18) "
        f"{out['R_exponent_amp0.02_R9to24']:.3f} (9-24); Derrick mu exponent {out['derrick']['exponent_loglog']:.3f}; "
        f"registry-vs-analytic h {out['worst_registry_vs_analytic_h']:.1e}; eq2 {out['worst_vs_eq2']:.1e}; |C_KT|/E_KT <= {out['worst_abs_C_over_E_KT']:.1e}")
    return out


# ================= 3. channels, the constant K0, the LP =================
def kt_C_batch(a0, X, M):
    """the omega^2 coefficient of K_T on the jet18 batch (exact quadratic)."""
    N = X.shape[0]
    Mb = np.broadcast_to(M, (N, 4, 4)).copy()
    A = FM.jets_batch(X, M)
    vals = {}
    for w in (1.0, -1.0, 0.0):
        A2 = A.copy(); A2[0] = w * a0
        vals[w] = kt_density_np(A2, Mb, P)
    return 0.5 * (vals[1.0] + vals[-1.0]) - vals[0.0]


def kt_Q_form(a0, M, n=18):
    Xs = [np.zeros(n)]
    for i in range(n):
        x = np.zeros(n); x[i] = 1.0; Xs.append(x)
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            x = np.zeros(n); x[i] = 1.0; x[j] = 1.0; Xs.append(x); pairs.append((i, j))
    c = kt_C_batch(a0, np.array(Xs), M)
    c0 = float(c[0]); d = c[1:1 + n] - c[0]
    Q = np.diag(d)
    for m, (i, j) in enumerate(pairs):
        Q[i, j] = Q[j, i] = 0.5 * (c[1 + n + m] - c[0] - d[i] - d[j])
    return Q, c0


def uframe_time_row(a0, M):
    """(time row, time column) of a0 in the frame where u = e0 (Lambda u = e0)."""
    u, _ = R2A.timelike_u(M)
    if u[0] < 0:
        u = -u
    Lam = R2L.boost_to_frame(u)
    au = Lam @ a0 @ Lam.T
    return au[0, 1:], au[1:, 0]


def stage_channels(rng):
    Mv = np.diag([G, 1.0, DELTA, 0.0])
    chans = {nm: (a0, Mv) for nm, a0 in FM.all_channels(Mv).items()}          # the 19 R2 channels
    ext = R2A.channel_set(G, DELTA, rng, extended=True)
    extras = {f"aud_{nm}": v for nm, v in ext.items() if nm not in chans}        # the audit's extras
    chans.update(extras)
    out = {"n_channels": len(chans), "n_R2": 19, "n_audit_extras": len(extras), "per_channel": {}}
    Qs, H0 = {}, {}
    worstQ = worst_id = 0.0
    for nm, (a0, M) in chans.items():
        Q, c0 = kt_Q_form(a0, M)
        Qn = float(np.max(np.abs(Q)))
        tr, tc = uframe_time_row(a0, M)
        sym = float(np.max(np.abs(a0 - a0.T)) / max(np.max(np.abs(a0)), 1e-300))
        c0_pred = -2.0 * float(np.sum(tr * tc))          # general: -2 sum_j a0_0j a0_j0 (= -2 sum a0_0j^2 symmetric)
        idr = abs(c0 - c0_pred) / max(abs(c0_pred), 1e-12)
        worstQ = max(worstQ, Qn); worst_id = max(worst_id, idr)
        # the lambda-family forms on the same channel (Lagrangian C forms)
        Qf, c0f = FM.Q_forms(a0, M, P)
        Qs[nm] = Qf
        out["per_channel"][nm] = {
            "C_KT_omega2_constant": c0, "C_KT_predicted_-2sum_a0_0j^2": c0_pred, "identity_rel": idr,
            "Q_KT_form_max_abs": Qn, "u_frame_time_row_norm": float(np.linalg.norm(tr)),
            "K0_energy_per_c2": -c0, "a0_antisymmetry_rel": sym,
            "a0_symmetric": bool(sym < 1e-12),
            "sign_check": ("no time row (K0 = 0)" if abs(c0) <= 1e-10 else
                           ("boost-like symmetric (K0 > 0)" if -c0 > 0 else
                            "ANTISYMMETRIC probe (K0 < 0: not a physical jet direction, d_0 M is symmetric)")),
            "min_eig_H2_lambda": {f"{lam:g}": FM.min_eig(-4.0 * ((1 - lam) * Qf["I1"] + lam * Qf["I1_h"]))
                                  for lam in (0.0, 0.25, 0.5, 0.75, 1.0)}}
    out["worst_Q_KT_form_norm"] = worstQ
    out["worst_identity_rel"] = worst_id
    out["claim_constant"] = ("the omega^2 coefficient of K_T is jet-INDEPENDENT on every channel (Q_KT = 0 to machine precision) and equals "
                             "-2 sum_j a0_0j a0_j0 in the u-frame: -2 sum_j a0_0j^2 <= 0 on every SYMMETRIC a0 (physical), "
                             "+2 sum_j a0_0j^2 on the antisymmetric probes (clock_t, clock_tr, x_boost_k_probe, aud_clock_random_*: "
                             "coms(gamma, M) is antisymmetric in its boost part)")
    out["symmetric_channels"] = [nm for nm, v in out["per_channel"].items() if v["a0_symmetric"]]
    out["antisymmetric_probes"] = [nm for nm, v in out["per_channel"].items() if not v["a0_symmetric"]]
    log(f"channels: {len(chans)}; Q_KT form norm {worstQ:.2e}; identity {worst_id:.2e}")
    for nm, v in out["per_channel"].items():
        log(f"   {nm:28s} K0/c2 = {v['K0_energy_per_c2']:10.4g}   |a0 time row| = {v['u_frame_time_row_norm']:.4g}   {'sym ' if v['a0_symmetric'] else 'ANTI'} {v['sign_check']}")

    # ---- the (lambda, c2) PSD set ----
    def H2(nm, lam):
        return -4.0 * ((1 - lam) * Qs[nm]["I1"] + lam * Qs[nm]["I1_h"])

    def K0(nm, c2):
        return c2 * out["per_channel"][nm]["K0_energy_per_c2"]

    lam_grid = [0.0, 0.25, 0.4, 0.49, 0.5, 0.51, 0.75, 1.0]
    c2_grid = [0.0, 0.1, 1.0, 10.0, 100.0, 1000.0]
    scales = [1.0, 10.0, 100.0, 1000.0]
    psd = {}
    for lam in lam_grid:
        row = {}
        for c2 in c2_grid:
            # min over channels and over |x| = s of x^T H2 x + c2 K0 (eigen-exact on the sphere)
            worst_by_scale = {}
            for s in scales:
                w = min(FM.min_eig(H2(nm, lam)) * s * s + K0(nm, c2) for nm in out["symmetric_channels"])
                worst_by_scale[f"{s:g}"] = float(w)
            form_psd = all(FM.min_eig(H2(nm, lam)) >= -1e-9 * max(1.0, abs(np.linalg.eigvalsh(H2(nm, lam))[-1])) for nm in chans)
            row[f"c2_{c2:g}"] = {"min_energy_omega2_by_jet_scale": worst_by_scale,
                                 "form_PSD_all_channels": bool(form_psd),
                                 "nonneg_at_all_scales": bool(min(worst_by_scale.values()) >= -1e-9)}
        psd[f"lam_{lam:g}"] = row
    out["psd_set"] = psd
    # exact lambda* (bisection on the min eig over all channels), c2-independent
    def all_psd(lam):
        return all(FM.min_eig(H2(nm, lam)) >= -1e-9 * max(1.0, abs(np.linalg.eigvalsh(H2(nm, lam))[-1])) for nm in chans)
    lo, hi = 0.0, 1.0
    assert not all_psd(lo) and all_psd(hi)
    for _ in range(50):
        m = 0.5 * (lo + hi)
        if all_psd(m): hi = m
        else: lo = m
    out["lambda_star_form"] = hi
    # at lambda = 0: the jet scale at which c2 K0 stops covering the negative form, per boost channel
    fail = {}
    for nm in chans:
        me = FM.min_eig(H2(nm, 0.0))
        if me < -1e-9 and out["per_channel"][nm]["K0_energy_per_c2"] > 0:
            fail[nm] = {"min_eig_H2_lam0": me, "K0_per_c2": out["per_channel"][nm]["K0_energy_per_c2"],
                        "jet_scale_where_c2_fails": {f"c2_{c2:g}": float(np.sqrt(K0(nm, c2) / -me)) for c2 in (0.1, 1.0, 10.0, 100.0)}}
    out["lambda0_c2_failure_scales"] = fail
    out["c2_never_hurts"] = bool(all(out["per_channel"][nm]["K0_energy_per_c2"] >= -1e-12 for nm in out["symmetric_channels"]))
    out["c2_never_hurts_scope"] = "every SYMMETRIC a0 channel (physical jets); the antisymmetric probes carry K0 < 0 by the sign flip and are listed separately"
    out["c2_alone_stabilizes"] = bool(psd["lam_0"]["c2_1000"]["nonneg_at_all_scales"])
    out["verdict"] = ("c2 K_T adds a jet-independent constant K0 >= 0 to the omega^2 energy coefficient (never hurts); "
                      "it cannot make the indefinite lambda < 1/2 form PSD (fails at jet scale sqrt(c2 K0 / |min eig|)); "
                      f"lambda* = {hi:.6f} with or without c2: c2 does NOT replace lambda >= 1/2")
    log(f"LP: lambda* = {hi:.6f}; c2 never hurts = {out['c2_never_hurts']}; c2 alone stabilizes = {out['c2_alone_stabilizes']}")
    return out


# ================= 4. the certified hedgehog =================
def stage_hedgehog():
    out = {}
    # (a) the stored certified 3x3 end state (M5.21.11, s = +1, g = 8), embedded
    cfg = B3.base_cfg(s=1.0)
    p1 = L0.default_params(s=1.0)
    M4 = B3.embed34(L0.load_stored3(), cfg)
    e0 = L0.certified_energy(M4, cfg, p1)
    ekt = kt_lattice(M4, cfg)
    out["stored3"] = {"file": L0.STORED3_NPZ, "E_cert": float(e0), "E_KT": ekt,
                      "E_cert_plus_c2KT_minus_E_cert": {f"c2_{c2:g}": float(c2 * ekt) for c2 in (0.1, 1.0, 10.0)},
                      "max_time_row": float(np.max(np.abs(M4[..., 0, 1:])))}
    # (b) the R4 family hedgehog (the analytic rotation-only dressing at g = 32)
    fam = R4.LocFamily(32, 48.0)
    eu, ev = B3.e_parts(fam.Mb, fam.cfg)
    kk = B3.kin_of(fam.Mb, fam.a0, fam.cfg)
    out["r4_hedgehog_g32"] = {"E_u": float(eu), "E_V": float(ev), "kin": float(kk),
                              "E_KT": kt_lattice(fam.Mb, fam.cfg),
                              "C_KT_omega2_on_clock": L0.omega_decompose(KT, fam.Mb, fam.cfg, P, fam.a0)[2]}
    out["statement"] = ("K_T = 0 and its omega^2 coefficient = 0 on every block-diagonal static field with its clock flow "
                        "(uniform time row, u = e0): the hedgehog static energy, its Derrick balance, the 3-lepton census "
                        "and the Coulomb anchors (all 3x3-sector measurements) are untouched by -c2 K_T for ANY c2, by construction")
    out["pass"] = bool(abs(ekt) <= 1e-20 * max(abs(e0), 1.0) and abs(out["r4_hedgehog_g32"]["E_KT"]) <= 1e-20 * max(eu, 1.0))
    log(f"hedgehog: stored3 E_cert {e0:.4f} E_KT {ekt:.2e}; R4 g=32 E_u {eu:.4f} E_KT {out['r4_hedgehog_g32']['E_KT']:.2e}")
    return out


# ================= 5. the c2 localization estimate =================
def stage_localize():
    with open(os.path.join(DATA, "m5_32_r4_audit_clock_box_n64_L96.json")) as f:
        box = json.load(f)
    lad = box["kinds"]["exp2"]["ladder"]
    Rs = np.array(box["R_grid"])
    amps = np.array(box["amps"])
    fam = R4.LocFamily(64, 96.0)
    amp_ref = 0.02
    ekt = []
    for R in Rs:
        Qb = R2L.qb_field(fam.cfg, fam.b_of(amp_ref, R, "exp2"))
        ekt.append(kt_lattice(B3.sym4(R2L.conj(Qb, fam.Mb)), fam.cfg))
        log(f"localize: n64 L96 E_KT(R={R:g}, amp={amp_ref}) = {ekt[-1]:.5g}")
    ekt = np.array(ekt)
    kR = ekt / amp_ref ** 2                               # E_KT(amp, R) = amp^2 kR(R) (amp^2 law, stage 2)
    fitk = np.polyfit(np.log(Rs[:4]), np.log(ekt[:4]), 1)
    out = {"source_gain": "m5_32_r4_audit_clock_box_n64_L96.json (exp2 ladder)", "n": 64, "L": 96.0,
           "R_grid": Rs.tolist(), "amp_ref": amp_ref, "E_KT_amp_ref": ekt.tolist(),
           "E_KT_R_exponent_R6to18": float(fitk[0]), "E_KT_R_exponent_R6to48": float(np.polyfit(np.log(Rs), np.log(ekt), 1)[0]),
           "k_lin": float(np.mean(ekt[:4] / Rs[:4])), "per_lambda": {}}
    c2_grid = np.geomspace(0.01, 10.0, 13)
    for lam in (0.75, 1.0):
        Es = {}
        for R in Rs:
            r = lad[f"R_{R:g}"]
            e0 = np.array(r["E_stat_lam0"]); e1 = np.array(r["E_stat_lam1"])
            Es[R] = (1 - lam) * e0 + lam * e1
        E_ref = Es[Rs[0]][0]
        gain = np.array([np.min(Es[R]) - E_ref for R in Rs])
        amp_star = np.array([amps[int(np.argmin(Es[R]))] for R in Rs])
        neg = gain < 0
        sfit = np.polyfit(np.log(Rs[neg]), np.log(-gain[neg]), 1)
        a_g, s_g = float(np.exp(sfit[1])), float(sfit[0])
        amp_use = float(np.median(amp_star))
        row = {"amp_star_per_R": amp_star.tolist(), "amp_used_for_E_KT": amp_use, "envelope_gain": gain.tolist(),
               "gain_powerlaw": {"a": a_g, "s": s_g, "form": "G(R) = -a R^s"}, "scan": {}}
        for c2 in c2_grid:
            # (a) 1D in R at fixed amp: c2 E_KT(R; amp*) + G(R) on the grid
            tot = c2 * kR * amp_use ** 2 + gain
            i = int(np.argmin(tot))
            # continuum power-law prediction: R* = (a s / (c2 k))^(1/(1-s)) if s < 1
            kk = float(np.mean(kR[:4] * amp_use ** 2 / Rs[:4]))   # E_KT(R; amp*) = kk R
            Rstar_pl = float((a_g * s_g / (c2 * kk)) ** (1.0 / (1.0 - s_g))) if s_g < 1.0 else None
            # (b) joint (amp, R): min over the amp grid of E_stat(amp, R) - E_ref + c2 amp^2 kR(R), then over R
            EJ = np.array([np.min(Es[R] - E_ref + c2 * amps ** 2 * kR[j]) for j, R in enumerate(Rs)])
            aJ = np.array([amps[int(np.argmin(Es[R] - E_ref + c2 * amps ** 2 * kR[j]))] for j, R in enumerate(Rs)])
            ij = int(np.argmin(EJ))
            row["scan"][f"c2_{c2:.4g}"] = {
                "1D_total_on_grid": tot.tolist(), "1D_R_star_grid": float(Rs[i]),
                "1D_interior": bool(0 < i < len(Rs) - 1), "1D_undressed_wins": bool(tot[i] >= 0),
                "powerlaw_R_star": Rstar_pl,
                "joint_min_over_amp": EJ.tolist(), "joint_amp_star": aJ.tolist(), "joint_R_star_grid": float(Rs[ij]),
                "joint_interior": bool(0 < ij < len(Rs) - 1), "joint_undressed_wins": bool(EJ[ij] >= -1e-12)}
        out["per_lambda"][f"lam_{lam:g}"] = row
        log(f"localize lam {lam:g}: gain {np.round(gain, 2).tolist()} slope {s_g:.3f} amp* {amp_use}")
        for c2 in c2_grid:
            r = row["scan"][f"c2_{c2:.4g}"]
            log(f"   c2 {c2:8.4g}: 1D R* {r['1D_R_star_grid']:4g} ({'interior' if r['1D_interior'] else 'wall/edge'}"
                f"{', undressed wins' if r['1D_undressed_wins'] else ''}), power-law R* "
                f"{'none (s >= 1)' if r['powerlaw_R_star'] is None else f'{r['powerlaw_R_star']:.3g}'}; "
                f"joint R* {r['joint_R_star_grid']:4g} amp* {r['joint_amp_star'][int(np.argmin(r['joint_min_over_amp']))]:.4g}"
                f"{' (undressed wins)' if r['joint_undressed_wins'] else ''}")
    return out


def main():
    rng = np.random.default_rng(20260828 + 7)
    RES["term"] = {"name": KT.name, "definition": KT.definition, "hash12": KT.hash,
                   "definition_sha256": KT.definition_hash, "kind": KT.kind}
    RES["1_covariance"] = stage_covariance(rng)
    RES["2_static_identity"] = stage_static(rng)
    RES["2_family"] = stage_family()
    RES["3_channels_lp"] = stage_channels(rng)
    RES["4_hedgehog"] = stage_hedgehog()
    RES["5_localize"] = stage_localize()
    c = RES["claims"]
    c["G4_covariance"] = RES["1_covariance"]["pass"]
    c["static_identity"] = RES["2_static_identity"]["pass"] and RES["4_hedgehog"]["pass"]
    c["K_T_positive_on_family_and_C_nonpositive"] = RES["2_family"]["all_positive"] and RES["2_family"]["all_C_nonpositive"]
    c["C_KT_zero_on_R4_family"] = RES["2_family"]["C_KT_zero_on_family"]
    c["Derrick_exponent_1"] = abs(RES["2_family"]["derrick"]["exponent_loglog"] - 1.0) < 0.1
    c["omega2_coefficient_is_constant"] = RES["3_channels_lp"]["worst_Q_KT_form_norm"] < 1e-9
    c["c2_never_hurts"] = RES["3_channels_lp"]["c2_never_hurts"]
    c["c2_replaces_lambda_half"] = RES["3_channels_lp"]["c2_alone_stabilizes"]
    RES["runtime_s"] = round(time.time() - T0, 1)
    with open(OUT, "w") as f:
        json.dump(RES, f, indent=1, default=float)
    log(f"wrote {OUT}; claims {c}")


if __name__ == "__main__":
    main()
