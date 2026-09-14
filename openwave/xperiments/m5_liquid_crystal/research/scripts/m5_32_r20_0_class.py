"""M5.32 R20-0: the form level of the R20 ladder (ledger section 6.9): the
certified potential V4 and the author's degree-8 spectral potential V_spec
side by side, before any solve.

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1), N = M eta. The
certified stack (m5_21_3_a_4d.py, code branch s = -1): the vacuum
M_vac = diag(-sg, 1, delta, 0) with sg = s g, so on s = -1, g = 8 the
matrix is diag(8, 1, 0.3, 0) and the N-spectrum is q = (sg, 1, delta, 0)
= (-8, 1, 0.3, 0) (the ledger's convention line, section 1.1).
    V4     = w sum_{p=1..4} (tr N^p - C_p)^2,   C_p = sum_i q_i^p,  w = W1
    V_spec = gamma tr[P(N)^2],   P(x) = prod_i (x - q_i)    (report section 370)
Both vanish exactly on the target spectrum (the four power sums fix the
multiset of four eigenvalues; P(N) = 0 iff N is diagonalizable with the
target spectrum), both are degree 8 in M.
(a) Hessians on the 10-dimensional sym4 space at M_vac: the six off-diagonal
    directions are conjugation directions (dM = eta A M - M A eta, A
    antisymmetric, nonzero because the four N eigenvalues are distinct and
    m_0 + m_j != 0), so both Hessians vanish there; on the four diagonal
    directions
        H_V4     = D (2 w J^T J) D,  J_pi = p q_i^(p-1), D = eta   (full rank, not diagonal)
        H_V_spec = 2 gamma diag(P'(q_i)^2)
    computed here by Richardson-extrapolated central differences of the
    exact gradients and compared with the closed forms, on the stack's
    roots (-8, 1, 0.3, 0) and on the report's roots (8, 1, 0.3, 0).
(b) The exact lattice gradient of V_spec: d tr[P(N)^2] = 2 tr[P P' dN]
    = 2 tr[eta P(N) P'(N) dM], so dV_spec/dM = 2 gamma sym4[eta P P'] (the
    eta placement is load-bearing; sym4 absorbs the transpose),
    with P'(N) = sum_k prod_{i != k} (N - q_i). Gated by complex step
    (exact for a polynomial) and by a 4-point real stencil.
(c) gamma fixed by matching the potential energy of the record's relaxed
    undressed single (R19-1, I1_un_single_d0_n32_g8.npz, E 13.83) under
    V_spec to its V4 energy: gamma = E_V4 / (h^3 sum_cells tr[P(N)^2]).
(d) The flat-vacuum quadratic form of the certified curvature energy on
    the 30 jets (3 directions x 10 sym components): E_I1[M_vac + eps x_a B]
    scales as eps^4 (F_ij = [A_i, A_j]_eta is bilinear in the jets), so the
    quadratic form is zero on every jet; measured by the scaling exponent.
(e) The measured pair laws on record (no run), for the README's
    "quartic kinetic term gives a linear potential" line.
(f) The trace direction (Maciej's 2026-09-13 note): the curvature energy
    is invariant under N -> N + c(x) 1, so tr N is fixed by the potential
    alone; the trace curvature of each potential at the vacuum.

Out: ../data/m5_32_r20_0_class.json, ../plots/m5_32_r20_0_hessians.png
Usage: python3 m5_32_r20_0_class.py
"""
from __future__ import annotations

import os
os.environ["OMP_NUM_THREADS"] = "1"

import importlib.util  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
OUT_JSON = os.path.join(DATA, "m5_32_r20_0_class.json")
OUT_PNG = os.path.join(PLOTS, "m5_32_r20_0_hessians.png")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    argv = sys.argv
    sys.argv = [argv[0]]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


EN = _load("m5_32_r19_entrants", "m5_32_r19_entrants.py")
LAG = EN.LAG
B3 = EN.B3
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
W1 = B3.W1
G_MAIN, DELTA = 8.0, 0.3
RECORD_SINGLE = os.path.join(DATA, "m5_32_r19_1", "I1_un_single_d0_n32_g8.npz")
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ================= the two potentials (per cell, any target spectrum) =================
def roots_of(cfg, degenerate=False):
    """the N-spectrum the stack targets: (sg, 1, delta, 0), or (sg, 1, delta, delta)."""
    sg, d = cfg["sg"], cfg["delta"]
    return (sg, 1.0, d, d) if degenerate else (sg, 1.0, d, 0.0)


def v4_energy_grad(M, cfg, q, w=W1, need_grad=True):
    """V4 = w sum_p (tr N^p - C_p)^2 with C_p = sum_i q_i^p; h^3-weighted total and sym4 gradient."""
    h3 = cfg["h"] ** 3
    Me = M @ ETA
    pows = [np.broadcast_to(np.eye(4), M.shape).copy()]
    for k in range(1, 4):
        pows.append(pows[-1] @ Me)
    t = [np.einsum("...kk->...", P @ Me) for P in pows]
    cp = [sum(qi ** p for qi in q) for p in range(1, 5)]
    E = h3 * w * np.sum(sum((t[k] - cp[k]) ** 2 for k in range(4)))
    E = E if np.iscomplexobj(E) else float(E)
    if not need_grad:
        return E, None
    GV = np.zeros_like(M)
    for k in range(1, 5):
        coef = 2.0 * w * (t[k - 1] - cp[k - 1]) * k
        X = ETA @ pows[k - 1]
        GV += coef[..., None, None] * X.swapaxes(-1, -2)
    return E, h3 * B3.sym4(GV)


def _poly_mats(N, q):
    I = np.broadcast_to(np.eye(4), N.shape)
    facs = [N - qi * I for qi in q]
    P = facs[0] @ facs[1] @ facs[2] @ facs[3]
    Pp = np.zeros_like(N)
    for k in range(4):
        R = I.copy()
        for i in range(4):
            if i != k:
                R = R @ facs[i]
        Pp += R
    return P, Pp


def vspec_density(M, q):
    """tr[P(N)^2] per cell."""
    P, _ = _poly_mats(M @ ETA, q)
    return np.einsum("...ij,...ji->...", P, P)


def vspec_energy_grad(M, cfg, q, gamma, need_grad=True):
    """V_spec = gamma h^3 sum_cells tr[P(N)^2]; gradient 2 gamma h^3 sym4[(eta P P')^T]."""
    h3 = cfg["h"] ** 3
    P, Pp = _poly_mats(M @ ETA, q)
    E = gamma * h3 * np.sum(np.einsum("...ij,...ji->...", P, P))
    E = E if np.iscomplexobj(E) else float(E)
    if not need_grad:
        return E, None
    G = 2.0 * (ETA @ (P @ Pp)).swapaxes(-1, -2)
    return E, gamma * h3 * B3.sym4(G)


# ================= (a) the Hessians on sym4 =================
def sym_basis():
    B = []
    for i in range(4):
        E = np.zeros((4, 4)); E[i, i] = 1.0; B.append(E)
    for i in range(4):
        for j in range(i + 1, 4):
            E = np.zeros((4, 4)); E[i, j] = E[j, i] = 1.0 / np.sqrt(2.0); B.append(E)
    return B


def hessian_fd(grad_fn, M0, basis, h=1e-3):
    """H_ab = d/dx_b (G . B_a) by Richardson-extrapolated central differences (h, h/2)."""
    def col(hh):
        H = np.zeros((len(basis), len(basis)))
        for b, Bb in enumerate(basis):
            Gp = grad_fn(M0 + hh * Bb)
            Gm = grad_fn(M0 - hh * Bb)
            for a, Ba in enumerate(basis):
                H[a, b] = float(np.sum((Gp - Gm) * Ba)) / (2.0 * hh)
        return H
    H1, H2 = col(h), col(h / 2.0)
    return (4.0 * H2 - H1) / 3.0


def closed_forms(q, w, gamma):
    q = np.array(q, dtype=float)
    J = np.array([[p * qi ** (p - 1) for qi in q] for p in range(1, 5)])
    D = np.diag([-1.0, 1.0, 1.0, 1.0])            # n_i = eta_i m_i: the M-coordinate Hessian is D (2 w J^T J) D
    H4 = D @ (2.0 * w * J.T @ J) @ D
    Pp = np.array([np.prod([qi - qj for j, qj in enumerate(q) if j != i]) for i, qi in enumerate(q)])
    Hs = 2.0 * gamma * Pp ** 2
    return H4, Hs, Pp


def stage_a():
    out = {}
    basis = sym_basis()
    for label, s in (("stack_branch_s_minus_1", -1.0), ("report_roots_s_plus_1", 1.0)):
        cfg = B3.base_cfg(s=s, g=G_MAIN, n=2, L=2.0, delta=DELTA)
        q = roots_of(cfg)
        M0 = B3.vac4(cfg)[None, None, None]
        cfg1 = dict(cfg); cfg1["h"] = 1.0
        g4 = lambda M: v4_energy_grad(M, cfg1, q, 1.0)[1][0, 0, 0]          # noqa: E731  (w = 1)
        gs = lambda M: vspec_energy_grad(M, cfg1, q, 1.0)[1][0, 0, 0]       # noqa: E731  (gamma = 1)
        H4 = hessian_fd(g4, M0, basis)
        Hs = hessian_fd(gs, M0, basis)
        H4c, Hsc, Pp = closed_forms(q, 1.0, 1.0)
        r = {"roots_N_spectrum": list(map(float, q)), "M_vac_diag": list(map(float, np.diag(M0[0, 0, 0]))),
             "sym4_basis": "4 diagonal E_ii, then 6 off-diagonal (E_ij + E_ji) / sqrt 2, (i, j) lexicographic",
             "V4_w1": {"hessian_eigs": sorted(map(float, np.linalg.eigvalsh(H4))),
                       "diag_block_eigs": sorted(map(float, np.linalg.eigvalsh(H4[:4, :4]))),
                       "diag_block_closed_form_eigs": sorted(map(float, np.linalg.eigvalsh(H4c))),
                       "diag_block_max_abs_dev_from_closed_form": float(np.max(np.abs(H4[:4, :4] - H4c))),
                       "offdiag_block_max_abs": float(np.max(np.abs(H4[4:, 4:]))),
                       "cross_block_max_abs": float(np.max(np.abs(H4[:4, 4:]))),
                       "diag_block_matrix": H4[:4, :4].tolist(),
                       "with_W1_eigs": sorted(map(float, W1 * np.linalg.eigvalsh(H4c)))},
             "V_spec_gamma1": {"hessian_eigs": sorted(map(float, np.linalg.eigvalsh(Hs))),
                               "diag_block_diag": [float(Hs[i, i]) for i in range(4)],
                               "diag_block_closed_form": list(map(float, Hsc)),
                               "P_prime_at_roots": list(map(float, Pp)),
                               "diag_block_max_abs_dev_from_closed_form": float(np.max(np.abs(Hs[:4, :4] - np.diag(Hsc)))),
                               "diag_block_offdiag_max_abs": float(np.max(np.abs(Hs[:4, :4] - np.diag(np.diag(Hs[:4, :4]))))),
                               "offdiag_block_max_abs": float(np.max(np.abs(Hs[4:, 4:]))),
                               "cross_block_max_abs": float(np.max(np.abs(Hs[:4, 4:])))}}
        out[label] = r
        log(f"(a) {label}: q {q}  V_spec diag {np.round(Hsc, 3)}  V4 eigs {np.round(np.linalg.eigvalsh(H4c), 4)}")
    out["report_curvatures_quoted"] = [371866.88, 48.02, 5.229, 11.52]
    Hs_plus = np.array(out["report_roots_s_plus_1"]["V_spec_gamma1"]["diag_block_closed_form"])
    Hs_minus = np.array(out["stack_branch_s_minus_1"]["V_spec_gamma1"]["diag_block_closed_form"])
    quoted = np.array(out["report_curvatures_quoted"])
    out["report_curvatures_match_roots"] = {
        "(8, 1, 0.3, 0)_max_rel_dev": float(np.max(np.abs(Hs_plus - quoted) / quoted)),
        "(-8, 1, 0.3, 0)_max_rel_dev": float(np.max(np.abs(Hs_minus - quoted) / quoted)),
        "note": "(8, 1, 0.3, 0) and (-8, 1, 0.3, 0) differ by a single-root flip, under which P'(q_i)^2 changes (the two sets are distinguishable); "
                "the full reflection q -> -q leaves P'(q_i)^2 invariant, so the report's numbers are equally those of (-8, -1, -0.3, 0): "
                "the g eigenvalue of N sharing the sign of the 1 and delta eigenvalues, which on this stack is the s = +1 branch, not the certified s = -1 branch"}
    return out


# ================= (b) the gradient gates =================
def stage_b(seed=7):
    rng = np.random.default_rng(seed)
    cfg = B3.base_cfg(s=-1.0, g=G_MAIN, n=6, L=9.0, delta=DELTA)
    q = roots_of(cfg)
    M = B3.vac4(cfg)[None, None, None] + 0.3 * B3.sym4(rng.standard_normal((6, 6, 6, 4, 4)))
    out = {"n": 6, "amp": 0.3, "ndir": 6}
    for name, fn in (("V_spec", lambda X, ng=True: vspec_energy_grad(X, cfg, q, 1.0, ng)),
                     ("V4_targets", lambda X, ng=True: v4_energy_grad(X, cfg, q, W1, ng))):
        E, G = fn(M)
        worst_cs, worst_fd = 0.0, 0.0
        for _ in range(6):
            D = B3.sym4(rng.standard_normal(M.shape))
            dd = float(np.sum(G * D))
            # complex step (exact for a polynomial)
            hcs = 1e-20
            Ecs = fn(M + 1j * hcs * D, False)[0]
            cs = float(np.imag(Ecs)) / hcs if np.iscomplexobj(Ecs) else float("nan")
            # 4-point real stencil
            e = 1e-3
            f = [fn(M + k * e * D, False)[0] for k in (-2, -1, 1, 2)]
            fd = (f[0] - 8 * f[1] + 8 * f[2] - f[3]) / (12 * e)
            sc = max(abs(dd), 1e-300)
            worst_cs = max(worst_cs, abs(cs - dd) / sc)
            worst_fd = max(worst_fd, abs(fd - dd) / sc)
        out[name] = {"E": float(np.real(E)), "complex_step_max_rel_err": worst_cs, "stencil4_max_rel_err": worst_fd,
                     "PASS_complex_step_1e-12": bool(worst_cs < 1e-12), "PASS_stencil4_1e-8": bool(worst_fd < 1e-8)}
        log(f"(b) {name}: complex-step rel err {worst_cs:.2e}, 4-point stencil rel err {worst_fd:.2e}")
    # the V4 wrapper reproduces the stack's own V4 (B3.e_parts / LAG.v4_grad_np) on the same field
    Mr = np.real(M)
    _, ev = B3.e_parts(Mr, cfg)
    p = LAG.default_params(s=-1.0, g=G_MAIN)
    gv = cfg["h"] ** 3 * B3.sym4(LAG.v4_grad_np(Mr, p))
    E4, G4 = v4_energy_grad(Mr, cfg, q, W1)
    out["V4_wrapper_vs_stack"] = {"E_rel": float(abs(E4 - ev) / abs(ev)), "G_max_abs_dev": float(np.max(np.abs(G4 - gv))),
                                  "PASS": bool(abs(E4 - ev) / abs(ev) < 1e-13 and np.max(np.abs(G4 - gv)) < 1e-12 * np.max(np.abs(gv)))}
    log(f"(b) V4 wrapper vs stack: E rel {out['V4_wrapper_vs_stack']['E_rel']:.2e}, G dev {out['V4_wrapper_vs_stack']['G_max_abs_dev']:.2e}")
    return out


# ================= (c) gamma =================
def stage_c():
    cfg = B3.base_cfg(s=-1.0, g=G_MAIN, n=32, L=48.0, delta=DELTA)
    q = roots_of(cfg)
    M = np.load(RECORD_SINGLE)["M"]
    _, ev = B3.e_parts(M, cfg)
    E4, _ = v4_energy_grad(M, cfg, q, W1, need_grad=False)
    S = cfg["h"] ** 3 * float(np.sum(vspec_density(M, q)))
    gamma = float(ev) / S
    Es, _ = vspec_energy_grad(M, cfg, q, gamma, need_grad=False)
    out = {"record_field": os.path.relpath(RECORD_SINGLE, DATA), "record_E_total_R19": 13.834938211360924,
           "E_V4_stack": float(ev), "E_V4_wrapper": float(E4), "sum_h3_trP2": S, "gamma": gamma,
           "E_V_spec_at_gamma": float(Es), "gamma_check_rel": float(abs(Es - ev) / abs(ev))}
    # the curvature ratios the match implies (the softest and stiffest spectral curvatures under each potential)
    H4c, Hsc, _ = closed_forms(q, W1, gamma)
    out["curvatures_at_fixed_weights"] = {"V4_W1_eigs": sorted(map(float, np.linalg.eigvalsh(H4c))),
                                          "V_spec_gamma_diag": list(map(float, Hsc)),
                                          "roots": list(map(float, q))}
    # where the record field's V_spec density sits (core vs far) against V4's
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    d4 = cfg["h"] ** 3 * W1 * np.sum([(t - c) ** 2 for t, c in zip(LAG.v4_traces_np(M), [sum(qi ** p for qi in q) for p in range(1, 5)])], axis=0)
    ds = cfg["h"] ** 3 * gamma * vspec_density(M, q)
    for lab, msk in (("core_r_lt_4", r < 4.0), ("mid_4_12", (r >= 4.0) & (r < 12.0)), ("far_r_ge_12", r >= 12.0)):
        out[f"split_{lab}"] = {"V4": float(np.sum(d4[msk])), "V_spec": float(np.sum(ds[msk]))}
    log(f"(c) gamma {gamma:.6e} (E_V4 {float(ev):.6f}, sum h3 trP2 {S:.6e})")
    return out


# ================= (d) the flat-vacuum jets =================
def stage_d():
    cfg = B3.base_cfg(s=-1.0, g=G_MAIN, n=8, L=12.0, delta=DELTA)
    basis = sym_basis()
    X = B3.coords(cfg["n"], cfg["h"])
    Mv = np.broadcast_to(B3.vac4(cfg), (8, 8, 8, 4, 4)).copy()
    eps = (1e-1, 1e-2, 1e-3)
    rows = []
    worst_exp_dev = 0.0
    for ax in range(3):
        for a, Ba in enumerate(basis):
            Es = []
            for e in eps:
                M = Mv + e * X[ax][..., None, None] * Ba
                Es.append(EN.block_reads(M, cfg, "I1")["E_curv_I1"])
            # scaling exponent between the two smallest eps (E ~ eps^4 expected)
            expo = float(np.log(Es[1] / Es[2]) / np.log(eps[1] / eps[2])) if Es[2] > 0 else float("nan")
            quad = float(Es[2] / eps[2] ** 2)
            rows.append({"axis": ax, "sym_index": a, "E_curv": Es, "exponent": expo, "E_over_eps2_at_1e-3": quad})
            if np.isfinite(expo):
                worst_exp_dev = max(worst_exp_dev, abs(expo - 4.0))
    # two-axis pairs (the first nonvanishing commutator needs jets on two axes)
    rng = np.random.default_rng(3)
    pairs = []
    for _ in range(12):
        ax1, ax2 = rng.choice(3, size=2, replace=False)
        a, b = rng.integers(0, 10, size=2)
        Es2 = []
        for e in eps:
            M = Mv + e * X[ax1][..., None, None] * basis[a] + e * X[ax2][..., None, None] * basis[b]
            Es2.append(EN.block_reads(M, cfg, "I1")["E_curv_I1"])
        ex2 = float(np.log(Es2[1] / Es2[2]) / np.log(eps[1] / eps[2])) if Es2[2] > 0 else float("nan")
        pairs.append({"axes": [int(ax1), int(ax2)], "sym_indices": [int(a), int(b)], "E_curv": Es2, "exponent": ex2})
    C = rng.standard_normal((3, 10))
    Es = []
    for e in eps:
        M = Mv.copy()
        for ax in range(3):
            for a, Ba in enumerate(basis):
                M = M + e * C[ax, a] * X[ax][..., None, None] * Ba
        Es.append(EN.block_reads(M, cfg, "I1")["E_curv_I1"])
    expo = float(np.log(Es[1] / Es[2]) / np.log(eps[1] / eps[2]))
    n_quartic = sum(1 for r in rows if np.isfinite(r["exponent"]) and abs(r["exponent"] - 4.0) < 0.05)
    n_zero = sum(1 for r in rows if not np.isfinite(r["exponent"]) or r["E_curv"][2] == 0.0)
    n_pair_quartic = sum(1 for r in pairs if np.isfinite(r["exponent"]) and abs(r["exponent"] - 4.0) < 0.05)
    n_pair_zero = sum(1 for r in pairs if not np.isfinite(r["exponent"]))
    out = {"n": 8, "eps": list(eps), "jets": rows, "two_axis_pairs": pairs, "pairs_quartic": n_pair_quartic, "pairs_identically_zero": n_pair_zero,
           "random_combination": {"E_curv": Es, "exponent": expo},
           "jets_quartic": n_quartic, "jets_identically_zero": n_zero, "max_abs_exponent_dev_from_4": worst_exp_dev,
           "quadratic_form_max_E_over_eps2": float(max([r["E_over_eps2_at_1e-3"] for r in rows] + [Es[2] / eps[2] ** 2])),
           "statement": "derivative sector only (the second variation of L_cert about the flat vacuum is the V4 mass term, the four spectral curvatures at W1): a single jet gives F = 0 identically (one axis, no commutator); two-axis pairs and the random 30-jet combination are quartic in eps: the quadratic fluctuation form of L_cert is zero on the 30-jet space, so no propagator, TT or otherwise, exists at quadratic order about the flat vacuum (the R14-0 fact)"}
    log(f"(d) single jets: {n_zero} of 30 identically zero; two-axis pairs: {n_pair_quartic} quartic, {n_pair_zero} zero of 12; random combination exponent {expo:.3f}")
    return out


# ================= (e) the pair laws on record =================
def stage_e():
    return {"rows": [
        {"object": "the like pair on the degenerate vacuum (g, 1, delta, delta), g 8", "law": "E_int = -18.66 + 0.935 d (six points, d 6, 9, 12, 15, 18, 24; rms 0.59, no bend)", "source": "R18-3 (d 12 to 24) + R19-3 (d 6, 9) (ledger sections 6.7, 6.8)", "form": "linear (a string)"},
        {"object": "the certified split-vacuum like pair at g 32", "law": "E falls from +3869 (d 10) to +2784 (d 24)", "source": "R3 arm (ii), lambda = 0 rows (ledger sections 6.8 and 6.9; R19-1 JSON R3_lam0_dr_n32_g32)", "form": "monotone fall, repulsive"},
        {"object": "the boost-dressed (Gam) pair at g 8, n 32 L 48", "law": "E_int(d) positive; the tail fits 1/d (R^2 0.993) and ln d/d, not 1/d^5", "source": "R19-1 (ledger section 6.8)", "form": "repulsive, a 1/d tail"},
        {"object": "the undressed I1 pair at g 8 (the static part)", "law": "E_int positive and RISING with d (+17.2 to +20.1 over d 12 to 30), object-independent (the Coulomb identity)", "source": "R19-1 (ledger section 6.8; outcome STATIC_PART, attractive in the string sense, not a Newton read)", "form": "string-like rise"}],
        "readme_statement": "a quartic kinetic term gives a linear, not 1/r, potential (round501; the script is not online)",
        "reading": "on this stack the linear law is the degenerate-vacuum string (a line defect, not a propagator); on the split vacuum the g 32 like pair falls with d while the undressed g 8 pair rises with d (string-like), and the dressed pair has a 1/d tail; the quadratic fluctuation form about the flat vacuum is zero (stage d), so any propagator on L_cert is background-generated, as the README says"}


# ================= (f) the trace direction (Maciej, #186 2026-09-13 15:57 UTC: "the trace is a constraint, not a field") =================
def stage_f(seed=11):
    """N -> N + c(x) 1 leaves the certified curvature energy unchanged (the commutator kills the identity), so along
    dM = c(x) eta the static energy is potential only; the potential curvature there is 2 w |J 1|^2 for V4 and
    2 gamma sum_i P'(q_i)^2 for V_spec (P(q_i + c) = c P'(q_i) + O(c^2))."""
    rng = np.random.default_rng(seed)
    cfg = B3.base_cfg(s=-1.0, g=G_MAIN, n=8, L=12.0, delta=DELTA)
    q = roots_of(cfg)
    M = B3.vac4(cfg)[None, None, None] + 0.2 * B3.sym4(rng.standard_normal((8, 8, 8, 4, 4)))
    c = 0.3 * rng.standard_normal((8, 8, 8))
    Mc = M + c[..., None, None] * ETA                     # N + c 1 = (M + c eta) eta
    E1 = EN.block_reads(M, cfg, "I1")["E_curv_I1"]
    E2 = EN.block_reads(Mc, cfg, "I1")["E_curv_I1"]
    _, G, _ = EN.energy_grad(M, cfg, "I1", c=1.0, p=LAG.default_params(s=-1.0, g=G_MAIN))
    gv = cfg["h"] ** 3 * B3.sym4(LAG.v4_grad_np(M, LAG.default_params(s=-1.0, g=G_MAIN)))
    Gc = G - gv                                            # the curvature gradient alone
    proj = float(np.max(np.abs(np.einsum("...ab,ab->...", Gc, ETA))))   # its component along eta per cell
    # the potential curvatures along the trace direction at the vacuum
    J1 = np.array([p_ * sum(qi ** (p_ - 1) for qi in q) for p_ in range(1, 5)])
    Pp = np.array([np.prod([qi - qj for j, qj in enumerate(q) if j != i]) for i, qi in enumerate(q)])
    out = {"n": 8, "E_curv_before": float(E1), "E_curv_after_trace_shift": float(E2), "rel_change": float(abs(E2 - E1) / abs(E1)),
           "curvature_gradient_component_along_eta_max": proj, "curvature_gradient_scale": float(np.max(np.abs(Gc))),
           "trace_curvature_V4_w1": float(2.0 * np.dot(J1, J1)), "trace_curvature_V4_W1": float(2.0 * W1 * np.dot(J1, J1)),
           "trace_curvature_V_spec_gamma1": float(2.0 * np.sum(Pp ** 2)),
           "statement": "the curvature energy is exactly invariant under N -> N + c(x) 1 and its gradient has no component along eta, so the trace of N is fixed by the potential alone (an algebraic Euler-Lagrange relation, tr N = C_1 for V4 at the vacuum); both R20 potentials are stiff there"}
    log(f"(f) trace shift: E_curv rel change {out['rel_change']:.1e}, gradient along eta {proj:.1e} (scale {out['curvature_gradient_scale']:.2e}); trace curvature V4 (W1) {out['trace_curvature_V4_W1']:.1f}, V_spec (gamma 1) {out['trace_curvature_V_spec_gamma1']:.0f}")
    return out


def plot(res):
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axs[0]
    for lab, key, mk in (("stack roots (-8, 1, 0.3, 0)", "stack_branch_s_minus_1", "o"), ("report roots (8, 1, 0.3, 0)", "report_roots_s_plus_1", "s")):
        r = res["a"][key]
        ax.semilogy(range(4), r["V4_w1"]["diag_block_closed_form_eigs"], mk + "-", label=f"V4 (w = 1), {lab}")
        ax.semilogy(range(4), sorted(r["V_spec_gamma1"]["diag_block_closed_form"]), mk + "--", label=f"V_spec (gamma = 1), {lab}")
    ax.set_xlabel("spectral curvature, sorted"); ax.set_ylabel("Hessian eigenvalue"); ax.set_title("(a) the four spectral curvatures at the vacuum")
    ax.legend(fontsize=7)
    ax = axs[1]
    eps = res["d"]["eps"]
    for r in res["d"]["two_axis_pairs"]:
        if r["E_curv"][2] > 0:
            ax.loglog(eps, r["E_curv"], "-", color="gray", alpha=0.5, lw=0.8)
    ax.loglog(eps, res["d"]["random_combination"]["E_curv"], "ko-", label=f"random jet combination, exponent {res['d']['random_combination']['exponent']:.2f}")
    e = np.array(eps)
    ax.loglog(e, res["d"]["random_combination"]["E_curv"][0] * (e / e[0]) ** 4, "r:", label="eps^4 reference")
    ax.set_xlabel("jet amplitude eps"); ax.set_ylabel("E_curv(I1), n 8 L 12"); ax.set_title("(d) the flat-vacuum jets: quartic, no quadratic form")
    ax.legend(fontsize=7)
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=130)
    log(f"plot {OUT_PNG}")


def main():
    res = {"task": "M5.32 R20-0", "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "convention": "code branch s = -1: M_vac = diag(8, 1, 0.3, 0), N = M eta spectrum (-8, 1, 0.3, 0), W1 = %.12g" % W1}
    res["a"] = stage_a()
    res["b"] = stage_b()
    res["c"] = stage_c()
    res["d"] = stage_d()
    res["e"] = stage_e()
    res["f"] = stage_f()
    os.makedirs(DATA, exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1)
    plot(res)
    log(f"wrote {OUT_JSON}")
    return res


if __name__ == "__main__":
    main()
