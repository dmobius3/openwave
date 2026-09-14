"""M5.32 R21-0 adversarial audit: an independent attempt to refute the seven claim groups (a) to (g)
of the R21-0 form-level script (m5_32_r21_0_form.py and its data/m5_32_r21_0_form.json), written
from the definitions with this file's own methods: its own one-sided stencils and curvature
density, its own V4 (matrix-power traces) and V_spec (Horner) densities, complex-step and
Richardson derivatives of its OWN energy against the stack's gradient, its own restriction and a
face-loss identity for the n64 -> n32 seed, the closed-form factorization of the section 711.3
texture with a continuum quadrature for the -1/s constant, exact sympy valley slopes with an
implicit-function numerical cross-check, its own boost exponential (checked against scipy expm),
its own per-cell rescale with a c + m s^2 fit of the norm-preserving second difference, and its
own solid-angle degree by preimage counting on a BFS-oriented cube surface. The stack modules are
imported ONLY for the certified energy (m5_32_r20_1_axes energy_grad / energy_parts / density /
seed_axes, m5_21_3_a_4d as B3 for coords / pin_shell / sym4, m5_32_r3_ii_pair boost_at / conj),
the objects the claims are about; every READ is reimplemented here. This file calls no function
of m5_32_r21_0_form.py; it reads its JSON for the numbers under audit.

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4 per cell, eta = diag(-1, 1, 1, 1), N = M eta; code branch s = -1,
M_vac = diag(8, 1, delta, 0), N-spectrum q = (-8, 1, delta, 0), delta = 0.3 unless stated.
    E = 4 h^3 sum_br (1/2) sum_cells sum_{i<j} <F_ij, F_ij>_eta + V,  F_ij = A_i eta A_j - A_j eta A_i,
    A_i = d_i M (fwd and bwd one-sided branches, the missing end zero),
    <F, G>_eta = sum_ab eta_a eta_b F_ab G_ab,
    V4 = W1 h^3 sum_cells sum_{p=1..4} (tr N^p - C_p)^2, C_p = sum_i q_i^p,
    V_spec = gamma h^3 sum_cells tr[P(N)^2], P(x) = prod_i (x - q_i).
(a) dE = sum_ab G_ab dM_ab on symmetric dM (the stack's G), so a packed off-diagonal variable
    x_ab = M_ab = M_ba sees 2 G_ab. Tested with the complex step on THIS file's energy,
    dE/dt = Im E(M + i t D) / t (t = 1e-20, no cancellation), in directions supported on all
    cells, on the pinned shell only, on single face / edge / corner / interior entries, and with
    a Richardson central difference of the stack's own energy. The pin: the outer ceil(1.6/h) = 2
    planes per face, n^3 - (n - 4)^3 = 10816 of 32768 cells at n 32.
(b) x_i(n) = (i - (n - 1)/2) h: x_{i+16}(64) = x_i(32) exactly at h 1.5. The restricted field's
    n32 energy differs from the n64 density summed over the central cube ONLY by the lost
    one-sided branch at the cube faces: on the n32 box the bwd branch has A_ax = 0 on the low
    face and the fwd branch has A_ax = 0 on the high face along ax, so E_n32(restricted)
    = sum_cube e64 with those A_ax zeroed (an identity to roundoff, tested).
(c) M_s = M_vac + eps [E_00 f + (E_01 + E_10) q], f, q functions of u = x / s: A_i =
    eps (E_00 d_i f + B d_i q) with B = E_01 + E_10, so
        F_ij = eps^2 (d_i f d_j q - d_j f d_i q) F_0,  F_0 = E_00 eta B - B eta E_00 = -E_01 + E_10,
        E_curv = 4 h^3 sum_br (1/2) sum_cells eps^4 |grad f x grad q|^2 <F_0, F_0>_eta,
        <F_0, F_0>_eta = -2 (= tr(F_0 eta F_0^T eta)), tr(F_0 eta F_0 eta) = +2, tr(F_0 F_0^T) = +2.
    So E_curv is EXACTLY eps^4 times a shape integral, independent of delta and of g, zero for
    any two radial profiles (grad f parallel to grad q), and in the continuum
        E_curv(s) -> -8 eps^4 I / s,  I = integral |grad_u f x grad_u q|^2 d^3u,
    since |grad f x grad q|^2 scales as s^-4 and the volume as s^3. With the compact bump
    b(u) = exp(1 - 1/(1 - u^2)), b'(u) = -2 u b(u) / (1 - u^2)^2:
        (ii)  f = b, q = b u_x:      I_2 = (8 pi / 3) int_0^1 b^2 b'^2 u^2 du
        (iii) f = b u_x, q = b u_y:  I_3 = int_0^1 [(8 pi / 3) b'^2 b^2 u^4 + (16 pi / 3) b^3 b' u^3 + 4 pi b^4 u^2] du
    (angular integrals of sin^2 theta = 8 pi / 3 and of 1 = 4 pi), both also checked by a 3D
    midpoint quadrature. Two different radial profiles (b, b^2) have grad f parallel to grad q in
    the continuum (I = 0) but NOT on the lattice, where the one-sided differences of two different
    profiles are not parallel: the lattice E_curv of that control is measured against the p-wave's. V4 on the texture depends on the 2x2 block (0, 1) only, so it is also
    delta-independent: tr N^p - C_p = tr N_2x2^p - ((-8)^p + 1).
(d) Read from data/m5_32_r20_1_axes.json: c(g) = d2E/ds2 at s 0.02 on S_1 at g 8 (9000 steps)
    and g 32; ratio against (31/7)^2 and 16; the g 16 separation (15/7)^2 / 4 - 1. The two raw
    curvatures are also recomputed here with this file's own energy along the R3 boost.
(e) phi(nu, Delta) = V(nu + Delta, nu, nu, nu); on a valley branch (phi_nu = 0, phi_nunu > 0)
    nu'(Delta) = -phi_nuDelta / phi_nunu (implicit function theorem), cross-checked by re-solving
    phi_nu = 0 at Delta +- 1e-4. At a fully symmetric point of any permutation-symmetric V with
    V_ii = A, V_ij = B (i != j): phi_nuDelta = A + 3 B, phi_nunu = 4 A + 12 B, ratio 1/4 (derived
    here on a generic symmetric polynomial with symbolic coefficients). For V_spec at Delta -9 the
    valley is exact, (nu + Delta, nu, nu, nu) = (-8, 1, 1, 1), and
        nu' = -P'(-8)^2 / (P'(-8)^2 + 3 P'(1)^2).
(f) Q(x) = exp(s K_x), K_x = b*(r) n.K (n.K symmetric, (n.K)^3 = n.K, so exp = I + sinh K +
    (cosh - 1) K^2, checked against scipy expm), Q eta Q^T = eta, M_s = Q M Q^T,
    N_s = Q N Q^-1 (spectrum invariant). d/ds tr M_s^2 = 4 tr(M^2 K) (zero on block-diagonal M),
    d2/ds2 tr M_s^2 = 8 tr(K^2 M^2) + 8 tr(K M K M) (from M_s'' = K^2 M + 2 K M K + M K^2).
    Norm-preserving path M_s f_s, f_s = sqrt(tr M^2 / tr M_s^2) = 1 + c2 s^2 + O(s^4),
    c2 = -(d2/ds2 tr M_s^2) / (4 tr M^2) per cell (even in s: Q(-s) = eta Q(s) eta and the energy
    is eta-conjugation invariant). Expanding E about M:
        E_n(s) - E_0 = (1/2) c_raw s^2 + c2 s^2 <G, M> + (1/2) c2^2 s^4 H_V(M, M) + O(s^4),
    so the central second difference along the norm-preserving path is c + m s^2 with
    c = c_raw + 2 sum_cells c2 <G, M>: the s -> 0 limit is the RAW curvature up to the residual
    gradient term, and the s 0.02 number is that limit plus a quartic paid by the potential's
    stiffness along the scale direction, H_V4(M, M) = h^3 W1 sum_p [2 (p t_p)^2 + 2 (t_p - C_p)
    p (p - 1) t_p], H_Vspec(M, M) = gamma h^3 sum_i 2 [(l_i P'(l_i))^2 + P(l_i) l_i^2 P''(l_i)].
    A two-point fit D(s) = c0 + c1 s^2 on s 0.005 / 0.01 is validated at s 0.02 / 0.04 / 0.05 and c0
    is compared with c_raw; the script's least squares on s 0.01 / 0.02 / 0.05 is reproduced against
    its s_to_0 keys. Since D(s) is a CENTRAL SECOND DIFFERENCE, an even path E_0 + (1/2) c0 s^2 + k s^4
    gives D(s) = c0 + 2 k s^2, so the quartic wall is k = (c1_normpres - c1_raw) / 2 (the script's
    /12 would be right for the true second derivative 12 k s^2, which is not what was fitted); k is
    also measured directly as [E_n(s) - E_raw(s) - s^2 sum c2 <G, M>] / s^4 and as the s^4 scaling
    of the potential's excess dV(s) along the rescaled path (ratios 16 per doubling of s). The dip
    E_min = -c0^2 / (16 k) at s_min = sqrt(|c0| / 4 k) is evaluated DIRECTLY (E_n(s) - E_0 at the
    script's s_min and at this file's) on three rows.
(g) tr(M_3x3^2) = sum lam_i^2 (vacuum 1 + delta^2 = 1.09) on the axis cells rho < 0.75 h,
    4 < |z| < 22; the eigenvalue gaps there; the degree of the oriented winding eigenvector on
    the cube surface c0 +- 6 cells (r 9) by preimage counting: for 24 random directions v the
    signed number of oriented spherical triangles containing v (all three det(p, q, v),
    det(q, r, v), det(r, p, v) of one sign, the triangle on v's hemisphere), the surface oriented
    by a BFS over the surface graph.
    V_spec = gamma sum_i P(l_i)^2 with P(0) = 0 vanishes on the all-zero spatial block, so a
    melted axis is a V_spec vacuum; the V_spec and V4 densities on the melted cells are compared.

Every printed line is PASS/FAIL on a number that can go either way.
Out: ../data/m5_32_r21_0_audit.json
Usage: python3 m5_32_r21_0_audit.py   (about 3 minutes single-threaded; the n32 / n64 fields of data/m5_32_r20_1/)
"""

from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "1"

import importlib.util  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import heapq  # noqa: E402
from collections import deque  # noqa: E402

import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402
from scipy.integrate import quad  # noqa: E402
from scipy.linalg import expm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r21_0_audit.json")
FORM_JSON = os.path.join(DATA, "m5_32_r21_0_form.json")
R20_JSON = os.path.join(DATA, "m5_32_r20_1_axes.json")
R20_0_JSON = os.path.join(DATA, "m5_32_r20_0_class.json")
R20_NPZ = os.path.join(DATA, "m5_32_r20_1")
BSTAR_JSON = os.path.join(DATA, "m5_21_14_minimize.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    argv = sys.argv
    sys.argv = [argv[0]]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


R20 = _load(
    "m5_32_r20_1_axes", "m5_32_r20_1_axes.py"
)  # energy_grad, energy_parts, density, seed_axes: the certified stack
B3 = R20.B3  # coords, pin_shell, sym4, W1
R3 = R20.R3  # boost_at, conj: the record's dressing

E4 = np.array([-1.0, 1.0, 1.0, 1.0])
ETA = np.diag(E4)
W1 = float(B3.W1)
DELTA = 0.3
LAM = {"S1": (1.0, DELTA, 0.0), "Sd": (DELTA, 1.0, 0.0), "S0": (0.0, 1.0, DELTA)}
WRANK = {"S1": 2, "Sd": 1, "S0": 0}
T0 = time.time()
LINES = {}


def log(*a):
    print(f"[{time.time() - T0:6.1f}s]", *a, flush=True)


def line(key, ok, detail):
    LINES[key] = {"pass": bool(ok), "detail": detail}
    print(f"{'PASS' if ok else 'FAIL'} {key}: {detail}", flush=True)


def rel(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-300))


# ================= own energy layer =================
def own_coords(n, h):
    x = (np.arange(n) - 0.5 * (n - 1)) * h
    return np.meshgrid(x, x, x, indexing="ij")


def own_pin(n, h, depth=1.6):
    w = int(np.ceil(depth / h))
    P = np.ones((n, n, n), dtype=bool)
    P[w : n - w, w : n - w, w : n - w] = False
    return P


def own_d1(f, ax, h, side):
    out = np.zeros_like(f)
    hi = [slice(None)] * f.ndim
    lo = [slice(None)] * f.ndim
    hi[ax] = slice(1, None)
    lo[ax] = slice(0, -1)
    diff = (f[tuple(hi)] - f[tuple(lo)]) / h
    tgt = [slice(None)] * f.ndim
    tgt[ax] = slice(0, -1) if side == "fwd" else slice(1, None)
    out[tuple(tgt)] = diff
    return out


def own_curv_density_from_jets(Aby, h):
    """Aby = {side: [A_x, A_y, A_z]}: 4 h^3 sum_br (1/2) sum_{i<j} <F_ij, F_ij>_eta per cell."""
    e = None
    for side, A in Aby.items():
        for i in range(3):
            for j in range(i + 1, 3):
                F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                term = 0.5 * np.einsum("a,b,...ab,...ab->...", E4, E4, F, F, optimize=True)
                e = term if e is None else e + term
    return 4.0 * h**3 * e


def own_curv_density(M, h):
    return own_curv_density_from_jets(
        {s: [own_d1(M, ax, h, s) for ax in range(3)] for s in ("fwd", "bwd")}, h
    )


def own_traces(M):
    N = M @ ETA
    P = N
    t = [np.einsum("...ii->...", P)]
    for _ in range(3):
        P = P @ N
        t.append(np.einsum("...ii->...", P))
    return t


def own_v4_density(M, q, w, h3):
    t = own_traces(M)
    cp = [sum(qi**p for qi in q) for p in range(1, 5)]
    return h3 * w * sum((t[k] - cp[k]) ** 2 for k in range(4))


def poly_coeffs(q):
    c = np.array([1.0])
    for qi in q:
        c = np.convolve(c, [1.0, -qi])
    return c


def horner(N, c):
    I = np.broadcast_to(np.eye(4), N.shape)
    R = c[0] * I
    for ck in c[1:]:
        R = R @ N + ck * I
    return R


def own_vspec_density(M, q, gamma, h3):
    P = horner(M @ ETA, poly_coeffs(q))
    return gamma * h3 * np.einsum("...ij,...ji->...", P, P)


def own_energy(M, cfg, potk="v4std", gamma=None):
    """(E_curv, V) as sums of this file's densities; complex-safe (no float casts)."""
    h = cfg["h"]
    h3 = h**3
    q = (cfg["sg"], 1.0, cfg["delta"], 0.0)
    ec = np.sum(own_curv_density(M, h))
    ev = (
        np.sum(own_v4_density(M, q, W1, h3))
        if potk == "v4std"
        else np.sum(own_vspec_density(M, q, gamma, h3))
    )
    return ec, ev


def own_density(M, cfg, potk="v4std", gamma=None):
    h = cfg["h"]
    h3 = h**3
    q = (cfg["sg"], 1.0, cfg["delta"], 0.0)
    e = own_curv_density(M, h)
    return e + (
        own_v4_density(M, q, W1, h3) if potk == "v4std" else own_vspec_density(M, q, gamma, h3)
    )


def gamma_r20():
    with open(R20_0_JSON) as f:
        return float(json.load(f)["c"]["gamma"])


# ================= (a) the free boundary =================
def audit_a(FJ):
    out = {}
    cfg = R20.cfg_of(8, 12.0, 8.0)
    p = R20.params_of(8.0)
    h = cfg["h"]
    n = 8
    rng = np.random.default_rng(2111)
    M = R20.seed_axes(cfg, LAM["S1"]) + 0.06 * B3.sym4(rng.standard_normal((n, n, n, 4, 4)))
    M[..., 0, 1:] = 0.0
    M[..., 1:, 0] = 0.0
    E_stack, G, info = R20.energy_grad(M, cfg, p, None)
    ec, ev = own_energy(M, cfg)
    out["E_stack"] = float(E_stack)
    out["E_own"] = float(ec + ev)
    line(
        "A0_own_energy_equals_stack_energy_grad_E_on_the_n8_test_field_1e-12",
        abs(ec + ev - E_stack) / abs(E_stack) < 1e-12,
        f"own {float(ec + ev):.12f} vs stack {float(E_stack):.12f}",
    )
    pin = own_pin(n, h)

    def cstep(D, base=None, Gb=None):
        base = M if base is None else base
        Gb = G if Gb is None else Gb
        dd = float(np.sum(Gb * D))
        ec_, ev_ = own_energy(base + 1e-20j * D, cfg)
        return dd, float(np.imag(ec_ + ev_) / 1e-20)

    # a second base with M_0i != 0 (the gradient in the 0i directions vanishes at any block-diagonal field by the eta parity M_0i -> -M_0i)
    Mb = M.copy()
    Rb = 0.03 * B3.sym4(rng.standard_normal(M.shape))
    Mb[..., 0, 1:] += Rb[..., 0, 1:]
    Mb[..., 1:, 0] += Rb[..., 1:, 0]
    Eb, Gb, _ = R20.energy_grad(Mb, cfg, p, None)
    ecb, evb = own_energy(Mb, cfg)
    out["E_stack_M0i_base"] = float(Eb)
    out["E_own_M0i_base"] = float(ecb + evb)

    def blockdiag(D):
        D = B3.sym4(D)
        D[..., 0, 1:] = 0.0
        D[..., 1:, 0] = 0.0
        return D

    fam = {}
    for lab, nd, maskfn, zero0i in (
        ("all_cells", 6, lambda: np.ones((n, n, n), bool), True),
        ("pinned_shell_only", 6, lambda: pin, True),
        ("interior_only", 3, lambda: ~pin, True),
        ("all_cells_M0i_entries_only", 3, lambda: np.ones((n, n, n), bool), False),
        ("all_cells_all_entries_at_the_M0i_base", 3, lambda: np.ones((n, n, n), bool), None),
    ):
        worst, dds = 0.0, []
        for _ in range(nd):
            D = rng.standard_normal(M.shape) * maskfn()[..., None, None]
            if zero0i is True:
                D = blockdiag(D)
            elif zero0i is False:
                D = B3.sym4(D)
                D[..., 1:, 1:] = 0.0
                D[..., 0, 0] = 0.0
            else:
                D = B3.sym4(D)
            dd, cs = cstep(D) if zero0i is True else cstep(D, Mb, Gb)
            worst = max(worst, abs(cs - dd) / max(abs(dd), 1e-300))
            dds.append(dd)
        fam[lab] = {
            "n_dirs": nd,
            "worst_rel": worst,
            "cells": int(maskfn().sum()),
            "typical_abs_dd": float(np.mean(np.abs(dds))),
        }
    out["complex_step_families"] = fam
    line(
        "A1_stack_gradient_equals_complex_step_of_own_energy_1e-9_in_random_directions_on_ALL_cells_and_on_the_pinned_shell_only",
        fam["all_cells"]["worst_rel"] < 1e-9
        and fam["pinned_shell_only"]["worst_rel"] < 1e-9
        and fam["interior_only"]["worst_rel"] < 1e-9,
        f"worst rel: all cells {fam['all_cells']['worst_rel']:.1e} (6 dirs), pinned shell only {fam['pinned_shell_only']['worst_rel']:.1e} (6 dirs, {fam['pinned_shell_only']['cells']} cells), interior {fam['interior_only']['worst_rel']:.1e}",
    )
    line(
        "A1b_stack_gradient_exact_also_in_the_M0i_directions_and_in_general_directions_at_a_base_with_M0i_ne_0_1e-9",
        fam["all_cells_M0i_entries_only"]["worst_rel"] < 1e-9
        and fam["all_cells_all_entries_at_the_M0i_base"]["worst_rel"] < 1e-9
        and abs(out["E_stack_M0i_base"] - out["E_own_M0i_base"]) / abs(out["E_stack_M0i_base"])
        < 1e-12,
        f"worst rel: 0i-only dirs {fam['all_cells_M0i_entries_only']['worst_rel']:.1e}, general dirs {fam['all_cells_all_entries_at_the_M0i_base']['worst_rel']:.1e} at a base with |M_0i| up to {float(np.max(np.abs(Mb[..., 0, 1:]))):.3f} (own E vs stack rel {abs(out['E_stack_M0i_base'] - out['E_own_M0i_base']) / abs(out['E_stack_M0i_base']):.1e}); informational: the boost-dressed fields of (f) carry M_0i != 0",
    )
    cells = {
        "face": (0, 3, 4),
        "edge": (0, 0, 5),
        "corner": (7, 7, 7),
        "face_hi": (3, 7, 2),
        "interior": (3, 4, 4),
    }
    per = {}
    for kind, (i, j, k) in cells.items():
        worst_abs, worst_rel, gmax = 0.0, 0.0, 0.0
        for a, b in ((0, 0), (1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3)):
            D = np.zeros_like(M)
            D[i, j, k, a, b] = 1.0
            D[i, j, k, b, a] = 1.0
            dd, cs = cstep(D)
            expect = G[i, j, k, a, b] if a == b else 2.0 * G[i, j, k, a, b]
            worst_abs = max(worst_abs, abs(cs - dd), abs(dd - expect))
            worst_rel = max(worst_rel, abs(cs - dd) / max(abs(dd), 1e-6))
            gmax = max(gmax, abs(dd))
        per[kind] = {
            "cell": [i, j, k],
            "pinned": bool(pin[i, j, k]),
            "worst_abs": worst_abs,
            "worst_rel_floor_1e-6": worst_rel,
            "max_abs_entry_derivative": gmax,
        }
    out["per_entry"] = per
    line(
        "A2_single_entry_derivatives_at_face_edge_corner_cells_match_complex_step_1e-9_and_the_packing_rule_G_diag_2G_off",
        all(v["worst_rel_floor_1e-6"] < 1e-9 for v in per.values())
        and all(v["worst_abs"] < 1e-7 for v in per.values())
        and per["face"]["pinned"]
        and per["corner"]["pinned"]
        and not per["interior"]["pinned"],
        "; ".join(
            f"{k} {v['cell']} pinned {v['pinned']}: rel {v['worst_rel_floor_1e-6']:.1e}, |dE/dx| up to {v['max_abs_entry_derivative']:.1f}"
            for k, v in per.items()
        ),
    )
    # Richardson on the stack's own energy (need_grad False)
    worst_r = 0.0
    for _ in range(3):
        D = blockdiag(rng.standard_normal(M.shape))
        dd = float(np.sum(G * D))

        def cd(e):
            Ep = R20.energy_grad(M + e * D, cfg, p, None, need_grad=False)[0]
            Em = R20.energy_grad(M - e * D, cfg, p, None, need_grad=False)[0]
            return (Ep - Em) / (2 * e)

        rich = (4.0 * cd(5e-4) - cd(1e-3)) / 3.0
        worst_r = max(worst_r, abs(rich - dd) / abs(dd))
    out["richardson_stack_energy_worst_rel"] = worst_r
    line(
        "A3_richardson_central_difference_of_the_STACK_energy_matches_its_gradient_1e-9_in_all_cell_directions",
        worst_r < 1e-9,
        f"worst rel {worst_r:.1e} (3 dirs, steps 1e-3 and 5e-4)",
    )
    # the pin
    pin32 = own_pin(32, 1.5)
    own_count = 32**3 - (32 - 4) ** 3
    fa = FJ["a_free_boundary_gate"]
    out["pin_n32"] = {
        "own_count": int(pin32.sum()),
        "formula_n3_minus_n_minus_4_cubed": own_count,
        "own_frac": float(pin32.mean()),
        "script": fa["pin_n32"],
        "stack_pin_shell_agrees": bool(np.array_equal(pin32, B3.pin_shell(32, 1.5))),
    }
    line(
        "A4_pin_covers_10816_of_32768_n32_cells_33.0_percent_two_planes_per_face",
        int(pin32.sum()) == 10816 == own_count == fa["pin_n32"]["pinned_cells"]
        and abs(fa["pin_n32"]["frac"] - 10816 / 32768) < 1e-12
        and out["pin_n32"]["stack_pin_shell_agrees"],
        f"own {int(pin32.sum())} = 32^3 - 28^3 = {own_count}, frac {pin32.mean():.6f}; script {fa['pin_n32']['pinned_cells']} ({fa['pin_n32']['frac']:.6f}), depth cells {fa['pin_n32']['depth_cells']}",
    )
    line(
        "A5_script_stencil_numbers_below_1e-10_on_all_cells_and_interior",
        fa["all_cells_faces_included"]["stencil4_rel_e0.001"] < 1e-10
        and fa["interior_only"]["stencil4_rel_e0.001"] < 1e-10,
        f"script all cells {fa['all_cells_faces_included']['stencil4_rel_e0.001']:.1e}, interior {fa['interior_only']['stencil4_rel_e0.001']:.1e}; own complex step {fam['all_cells']['worst_rel']:.1e}",
    )
    return out


# ================= (b) the far seed =================
def audit_b(FJ, RJ):
    out = {}
    x64 = (np.arange(64) - 31.5) * 1.5
    x32 = (np.arange(32) - 15.5) * 1.5
    dev = float(np.max(np.abs(x64[16:48] - x32)))
    out["coords"] = {
        "max_abs_dev": dev,
        "x64_16": float(x64[16]),
        "x32_0": float(x32[0]),
        "stack_coords_dev": float(
            np.max(np.abs(B3.coords(64, 1.5)[0][16:48, 0, 0] - B3.coords(32, 1.5)[0][:, 0, 0]))
        ),
    }
    line(
        "B1_central_32_cells_of_n64_L96_coincide_with_the_n32_L48_cells_exactly_offset_16",
        dev == 0.0
        and out["coords"]["stack_coords_dev"] == 0.0
        and FJ["b_far_seed"]["coord_coincidence_max_abs"] == 0.0,
        f"own max |x64[16:48] - x32| = {dev}, first central x {x64[16]}; stack coords dev {out['coords']['stack_coords_dev']}; script {FJ['b_far_seed']['coord_coincidence_max_abs']}",
    )
    cfg32 = R20.cfg_of(32, 48.0, 8.0)
    cfg64 = R20.cfg_of(64, 96.0, 8.0)
    h = 1.5
    h3 = h**3
    q = (-8.0, 1.0, DELTA, 0.0)
    pin = own_pin(32, h)
    worst = {
        "E_total": 0.0,
        "E_curv": 0.0,
        "V": 0.0,
        "E_pin_shell": 0.0,
        "cube": 0.0,
        "n64_total": 0.0,
        "face_identity": 0.0,
    }
    for obj in ("S1", "Sd", "S0"):
        M64 = np.load(os.path.join(R20_NPZ, f"{obj}_v4std_n64_g8.npz"))["M"]
        Mr = M64[16:48, 16:48, 16:48]
        e32c = own_curv_density(Mr, h)
        e32v = own_v4_density(Mr, q, W1, h3)
        E32 = {"E_curv": float(e32c.sum()), "V": float(e32v.sum())}
        E32["E_total"] = E32["E_curv"] + E32["V"]
        E32["E_pin_shell"] = float((e32c + e32v)[pin].sum())
        # the n64 density, own, and the face-loss identity
        jets = {s: [own_d1(M64, ax, h, s) for ax in range(3)] for s in ("fwd", "bwd")}
        e64 = own_curv_density_from_jets(jets, h) + own_v4_density(M64, q, W1, h3)
        cube = float(e64[16:48, 16:48, 16:48].sum())
        jr = {s: [A[16:48, 16:48, 16:48].copy() for A in jets[s]] for s in jets}
        for ax in range(3):
            sl = [slice(None)] * 3
            sl[ax] = 0
            jr["bwd"][ax][tuple(sl)] = 0.0
            sl[ax] = 31
            jr["fwd"][ax][tuple(sl)] = 0.0
        e_mod = own_curv_density_from_jets(jr, h)
        face_id = abs(float(e_mod.sum()) - E32["E_curv"]) / E32["E_curv"]
        loss = cube - E32["E_total"]
        del jets, jr, e_mod
        fr = FJ["b_far_seed"][obj]["restricted_n32"]
        r64 = RJ["rows"][f"{obj}_v4std_n64_g8"]["end_reads"]["energy"]["E_total"]
        rec = {
            "own_restricted_n32": E32,
            "script_restricted_n32": {k: fr[k] for k in ("E_total", "E_curv", "V", "E_pin_shell")},
            "own_n64_total": float(e64.sum()),
            "R20_n64_E_total": r64,
            "own_n64_cube_sum": cube,
            "script_n64_cube_sum": FJ["b_far_seed"][obj]["n64_density_inside_the_central_32_cube"],
            "face_loss_cube_minus_restricted": loss,
            "face_identity_rel": face_id,
            "shell_frac_of_restricted": E32["E_pin_shell"] / E32["E_total"],
            "R20_B_seed_9000_E_total": FJ["b_far_seed"][obj]["R20_B_seed_n32_9000"]["E_total"],
            "R20_B_seed_9000_shell": FJ["b_far_seed"][obj]["R20_B_seed_n32_9000"]["E_pin_shell"],
        }
        for k in ("E_total", "E_curv", "V", "E_pin_shell"):
            worst[k] = max(worst[k], abs(E32[k] - fr[k]) / abs(fr[k]))
        worst["cube"] = max(worst["cube"], abs(cube - rec["script_n64_cube_sum"]) / cube)
        worst["n64_total"] = max(worst["n64_total"], abs(rec["own_n64_total"] - r64) / r64)
        worst["face_identity"] = max(worst["face_identity"], face_id)
        out[obj] = rec
        log(
            f"(b) {obj}: own restricted E {E32['E_total']:.6f} (shell {E32['E_pin_shell']:.6f}); cube sum {cube:.6f}, loss {loss:.6f}; face identity {face_id:.1e}"
        )
    out["worst_rel"] = worst
    line(
        "B2_own_n32_energy_of_the_restricted_field_matches_script_E_total_E_curv_V_1e-10_for_S1_Sd_S0",
        max(worst[k] for k in ("E_total", "E_curv", "V")) < 1e-10,
        f"worst rel E_total {worst['E_total']:.1e}, E_curv {worst['E_curv']:.1e}, V {worst['V']:.1e}; own {[round(out[o]['own_restricted_n32']['E_total'], 4) for o in ('S1', 'Sd', 'S0')]}",
    )
    line(
        "B3_own_pinned_shell_energy_of_the_restricted_field_matches_script_1e-10",
        worst["E_pin_shell"] < 1e-10,
        f"worst rel {worst['E_pin_shell']:.1e}; shell energies {[round(out[o]['own_restricted_n32']['E_pin_shell'], 4) for o in ('S1', 'Sd', 'S0')]} of {[round(out[o]['own_restricted_n32']['E_total'], 3) for o in ('S1', 'Sd', 'S0')]} (the R20 B_seed 9000 rows: {[round(out[o]['R20_B_seed_9000_shell'], 2) for o in ('S1', 'Sd', 'S0')]} of {[round(out[o]['R20_B_seed_9000_E_total'], 2) for o in ('S1', 'Sd', 'S0')]})",
    )
    line(
        "B4_own_n64_density_matches_the_R20_n64_total_1e-10_and_its_central_cube_sum_matches_script_1e-10",
        worst["n64_total"] < 1e-10 and worst["cube"] < 1e-10,
        f"worst rel n64 total {worst['n64_total']:.1e}, cube sum {worst['cube']:.1e}; cube sums {[round(out[o]['own_n64_cube_sum'], 4) for o in ('S1', 'Sd', 'S0')]}",
    )
    line(
        "B5_restricted_n32_energy_equals_the_n64_cube_sum_with_the_lost_one_sided_branch_zeroed_identity_1e-12",
        worst["face_identity"] < 1e-12,
        f"identity rel dev {worst['face_identity']:.1e}; the loss cube - restricted = {[round(out[o]['face_loss_cube_minus_restricted'], 5) for o in ('S1', 'Sd', 'S0')]} (positive: the cube sum exceeds the n32 energy by the missing half-branch at the faces)",
    )
    return out


# ================= (c) the 711.3 texture =================
def bump(u):
    out = np.zeros_like(u)
    m = np.abs(u) < 1.0
    out[m] = np.exp(1.0 - 1.0 / (1.0 - u[m] ** 2))
    return out


def dbump(u):
    out = np.zeros_like(u)
    m = np.abs(u) < 1.0
    out[m] = np.exp(1.0 - 1.0 / (1.0 - u[m] ** 2)) * (-2.0 * u[m] / (1.0 - u[m] ** 2) ** 2)
    return out


def texture_fields(kind, s, X, Y, Z):
    r = np.sqrt(X * X + Y * Y + Z * Z)
    u, ux, uy = r / s, X / s, Y / s
    if kind == "i":
        return bump(u), bump(u)
    if kind == "ii":
        return bump(u), bump(u) * ux
    if kind == "iii":
        return bump(u) * ux, bump(u) * uy
    if kind == "two_radial":
        return bump(u), bump(u) ** 2
    raise ValueError(kind)


def own_texture_curv(f, q, h, pairing):
    """4 h^3 sum_br (1/2) sum_cells |grad f x grad q|^2 pairing, own one-sided stencils on the scalar profiles (eps = 1)."""
    tot = 0.0
    for side in ("fwd", "bwd"):
        gf = [own_d1(f, ax, h, side) for ax in range(3)]
        gq = [own_d1(q, ax, h, side) for ax in range(3)]
        cross = sum(
            (gf[i] * gq[j] - gf[j] * gq[i]) ** 2 for i in range(3) for j in range(i + 1, 3)
        )
        tot += 0.5 * float(cross.sum())
    return 4.0 * h**3 * pairing * tot


def texture_M(kind, s, eps, delta, X, Y, Z):
    f, q = texture_fields(kind, s, X, Y, Z)
    M = np.zeros(X.shape + (4, 4))
    M[..., 0, 0] = 8.0 + eps * f
    M[..., 1, 1] = 1.0
    M[..., 2, 2] = delta
    M[..., 0, 1] = eps * q
    M[..., 1, 0] = eps * q
    return M


def audit_c(FJ):
    out = {}
    fc = FJ["c_texture_711_3"]
    # the jet identity, own
    A = np.zeros((4, 4))
    A[0, 0] = 1.0
    Bm = np.zeros((4, 4))
    Bm[0, 1] = Bm[1, 0] = 1.0
    F0 = A @ ETA @ Bm - Bm @ ETA @ A
    F0_expected = np.zeros((4, 4))
    F0_expected[0, 1] = -1.0
    F0_expected[1, 0] = 1.0
    pair = float(np.einsum("a,b,ab,ab->", E4, E4, F0, F0))
    trT = float(np.trace(F0 @ ETA @ F0.T @ ETA))
    trN = float(np.trace(F0 @ ETA @ F0 @ ETA))
    trF = float(np.trace(F0 @ F0.T))
    jet = fc["jet"]
    key_auth = "tr_F0_eta_F0T_eta_the_authors_-2"
    out["jet"] = {
        "F0_is_minus_E01_plus_E10": bool(np.array_equal(F0, F0_expected)),
        "pairing_eta": pair,
        "tr_F0_eta_F0T_eta": trT,
        "tr_F0_eta_F0_eta": trN,
        "tr_F0_F0T": trF,
        "script": {k: v for k, v in jet.items() if k != "note"},
        "script_has_transpose_key": key_auth in jet,
    }
    line(
        "C1_jet_identity_F0_=_-E01_+_E10_pairing_-2_=_tr_F0_eta_F0T_eta_and_+2_without_the_transpose_and_+2_frobenius_matching_script",
        out["jet"]["F0_is_minus_E01_plus_E10"]
        and pair == -2.0
        and trT == -2.0
        and trN == 2.0
        and trF == 2.0
        and key_auth in jet
        and jet[key_auth] == -2.0
        and jet.get("tr_F0_eta_F0_eta_no_transpose") == 2.0
        and jet.get("tr_F0_F0T_frobenius") == 2.0
        and jet.get("stack_inner_eta_F0_F0") == -2.0,
        f"own: pairing {pair}, tr(F0 eta F0^T eta) {trT}, tr(F0 eta F0 eta) {trN}, tr(F0 F0^T) {trF}; script keys present {key_auth in jet}: {out['jet']['script']}",
    )
    # continuum constants, own 1D quadratures and a 3D midpoint check
    b = lambda u: np.exp(1.0 - 1.0 / (1.0 - u * u)) if abs(u) < 1 else 0.0  # noqa: E731
    db = lambda u: b(u) * (-2.0 * u / (1.0 - u * u) ** 2) if abs(u) < 1 else 0.0  # noqa: E731
    I2, _ = quad(lambda u: (8 * np.pi / 3) * b(u) ** 2 * db(u) ** 2 * u**2, 0, 1, limit=200)
    I3, _ = quad(
        lambda u: (8 * np.pi / 3) * db(u) ** 2 * b(u) ** 2 * u**4
        + (16 * np.pi / 3) * b(u) ** 3 * db(u) * u**3
        + 4 * np.pi * b(u) ** 4 * u**2,
        0,
        1,
        limit=200,
    )
    ng = 120
    g1 = (np.arange(ng) + 0.5) / ng * 2.0 - 1.0
    GX, GY, GZ = np.meshgrid(g1, g1, g1, indexing="ij")
    U = np.sqrt(GX * GX + GY * GY + GZ * GZ)
    Us = np.where(U < 1e-12, 1e-12, U)
    B, dB = bump(U), dbump(U)
    uh = [GX / Us, GY / Us, GZ / Us]

    def grad_of(kind_f):
        if kind_f == "b":
            return [dB * uh[k] for k in range(3)]
        if kind_f == "bx":
            return [dB * GX * uh[0] + B, dB * GX * uh[1], dB * GX * uh[2]]
        if kind_f == "by":
            return [dB * GY * uh[0], dB * GY * uh[1] + B, dB * GY * uh[2]]

    def cross2(gf, gq):
        return sum((gf[i] * gq[j] - gf[j] * gq[i]) ** 2 for i in range(3) for j in range(i + 1, 3))

    dv = (2.0 / ng) ** 3
    I2_3d = float(cross2(grad_of("b"), grad_of("bx")).sum() * dv)
    I3_3d = float(cross2(grad_of("bx"), grad_of("by")).sum() * dv)
    I_two_radial_3d = float(
        cross2(grad_of("b"), [2.0 * B * gk for gk in grad_of("b")]).sum() * dv
    )  # f = b, q = b^2: grad q = 2 b grad b, parallel to grad f
    del GX, GY, GZ, U, Us, B, dB, uh
    out["continuum"] = {
        "I2_quad": I2,
        "I2_3d_midpoint_120": I2_3d,
        "I3_quad": I3,
        "I3_3d_midpoint_120": I3_3d,
        "I_two_radial_b_b2_3d": I_two_radial_3d,
        "law": "E_curv(s) s / eps^4 -> -8 I",
    }
    line(
        "C2_continuum_shape_integrals_by_1D_quadrature_and_3D_midpoint_agree_1e-3",
        abs(I2 - I2_3d) / I2 < 1e-3 and abs(I3 - I3_3d) / I3 < 1e-3,
        f"I2 {I2:.6f} (3D {I2_3d:.6f}), I3 {I3:.6f} (3D {I3_3d:.6f}); continuum constants -8 I: {-8 * I2:.5f} (ii), {-8 * I3:.5f} (iii)",
    )
    # own lattice closed form on the n32 L48 box
    cfg = R20.cfg_of(32, 48.0, 8.0)
    p = R20.params_of(8.0)
    h = cfg["h"]
    X, Y, Z = own_coords(32, h)
    S_LAD = (12.0, 9.0, 6.0, 4.5, 3.0)
    own = {}
    for kind in ("i", "ii", "iii", "two_radial"):
        own[kind] = [own_texture_curv(*texture_fields(kind, s, X, Y, Z), h, pair) for s in S_LAD]
    out["own_lattice_E_curv_eps1"] = {
        k: dict(zip([f"s{s:g}" for s in S_LAD], v)) for k, v in own.items()
    }
    line(
        "C3_identity_control_f_eq_q_gives_EXACTLY_zero_in_the_closed_form_on_the_lattice",
        all(v == 0.0 for v in own["i"]),
        f"(i) E_curv at s {list(S_LAD)}: {own['i']}",
    )
    ratio_tr = [own["two_radial"][k] / own["ii"][k] for k in range(5)]
    line(
        "C3b_FINDING_two_DIFFERENT_radial_profiles_b_and_b2_give_a_NONZERO_negative_lattice_E_curv_above_10_percent_of_the_pwave_value_at_s_12_though_their_continuum_integral_is_zero",
        all(v < 0.0 for v in own["two_radial"])
        and abs(ratio_tr[0]) > 0.1
        and abs(I_two_radial_3d) < 1e-12,
        f"(b, b^2) E_curv at s {list(S_LAD)}: {[f'{v:.4f}' for v in own['two_radial']]} = {[f'{100 * v:.0f}' for v in ratio_tr]} percent of the p-wave (ii) values; continuum integral {I_two_radial_3d:.1e}: the form script's 'any two radial profiles' statement holds in the continuum only, on this lattice the one-sided differences of two different radial profiles are not parallel (an O(h/s) stencil artifact)",
    )
    # stack calls: energy_parts on own textures
    checks = [
        ("ii", 6.0, 1.0, 0.3),
        ("ii", 12.0, 2.0, 0.3),
        ("iii", 4.5, 1.0, 0.3),
        ("ii", 6.0, 0.5, 0.3),
        ("ii", 6.0, 4.0, 0.3),
        ("i", 6.0, 1.0, 0.3),
        ("two_radial", 6.0, 1.0, 0.3),
        ("ii", 6.0, 1.0, 1.0),
    ]
    stack = {}
    for kind, s, eps, delta in checks:
        cfgd = (
            R20.cfg_of(32, 48.0, 8.0)
            if delta == 0.3
            else B3.base_cfg(s=-1.0, g=8.0, n=32, L=48.0, delta=delta)
        )
        pd = (
            R20.params_of(8.0)
            if delta == 0.3
            else R20.LAG.default_params(s=-1.0, g=8.0, delta=delta)
        )
        Ms = texture_M(kind, s, eps, delta, X, Y, Z)
        ep = R20.energy_parts(Ms, cfgd, pd, None)
        own_c = eps**4 * own_texture_curv(*texture_fields(kind, s, X, Y, Z), h, pair)
        own_v = float(np.sum(own_v4_density(Ms, (-8.0, 1.0, delta, 0.0), W1, h**3)))
        stack[f"{kind}_s{s:g}_eps{eps:g}_d{delta:g}"] = {
            "stack_E_curv": ep["E_curv"],
            "own_closed_form": own_c,
            "stack_V4": ep["V"],
            "own_V4": own_v,
        }
    out["stack_vs_closed_form"] = stack
    nz = [k for k in stack if not k.startswith("i_") and not k.startswith("two_radial")]
    w_c = max(
        abs(stack[k]["stack_E_curv"] - stack[k]["own_closed_form"])
        / abs(stack[k]["own_closed_form"])
        for k in nz
    )
    w_v = max(
        abs(stack[k]["stack_V4"] - stack[k]["own_V4"]) / abs(stack[k]["own_V4"]) for k in stack
    )
    line(
        "C4_stack_E_curv_on_the_textures_equals_the_closed_form_eps4_|grad_f_x_grad_q|2_pairing_1e-10_and_stack_V4_equals_own_V4_1e-10",
        w_c < 1e-10 and w_v < 1e-10,
        f"worst rel E_curv {w_c:.1e} over {len(nz)} (kind, s, eps) points, V4 {w_v:.1e} over {len(stack)}",
    )
    zero_floor = abs(stack["i_s6_eps1_d0.3"]["stack_E_curv"])
    tr = stack["two_radial_s6_eps1_d0.3"]
    line(
        "C5_stack_E_curv_on_the_identity_control_is_zero_to_roundoff_below_1e-20_of_the_pwave_value_while_the_stack_confirms_the_two_radial_lattice_value_1e-10",
        zero_floor < 1e-20 * abs(stack["ii_s6_eps1_d0.3"]["stack_E_curv"])
        and abs(tr["stack_E_curv"] - tr["own_closed_form"]) < 1e-10 * abs(tr["own_closed_form"]),
        f"stack (i) {stack['i_s6_eps1_d0.3']['stack_E_curv']:.1e} vs (ii) {stack['ii_s6_eps1_d0.3']['stack_E_curv']:.4f}; stack (b, b^2) at s 6 {tr['stack_E_curv']:.6f} vs own closed form {tr['own_closed_form']:.6f}",
    )
    q4 = [stack[f"ii_s6_eps{e:g}_d0.3"]["stack_E_curv"] / e**4 for e in (0.5, 1.0, 4.0)]
    line(
        "C6_exact_quartic_stack_E_curv_over_eps4_at_eps_0.5_1_4_constant_to_1e-13_by_the_bilinear_F_argument",
        (max(q4) - min(q4)) / abs(np.mean(q4)) < 1e-13,
        f"E_curv / eps^4 at s 6: {[f'{v:.15f}' for v in q4]}, spread {(max(q4) - min(q4)) / abs(np.mean(q4)):.1e}",
    )
    d03, d1 = stack["ii_s6_eps1_d0.3"], stack["ii_s6_eps1_d1"]
    line(
        "C7_delta_independence_is_structural_E_curv_and_V4_identical_at_delta_0.3_and_1_to_1e-13",
        abs(d03["stack_E_curv"] - d1["stack_E_curv"]) < 1e-13 * abs(d03["stack_E_curv"])
        and abs(d03["stack_V4"] - d1["stack_V4"]) < 1e-13 * abs(d03["stack_V4"]),
        f"E_curv {d03['stack_E_curv']:.12f} vs {d1['stack_E_curv']:.12f}; V4 {d03['stack_V4']:.6f} vs {d1['stack_V4']:.6f} (the texture lives on the (0, 1) block; the (2, 2) entry cancels in tr N^p - C_p)",
    )
    # the -1/s law: own lattice E s at eps 1 vs the script's and vs the continuum
    for kind, tag, Ic in (
        ("ii", "delta0.3_ii_f_radial_q_pwave", I2),
        ("iii", "delta0.3_iii_f_px_q_py", I3),
    ):
        cs_own = np.array(own[kind]) * np.array(S_LAD)
        cs_scr = np.array(fc[tag]["eps1"]["E_curv_times_s"])
        above = cs_own[:4]
        spread = float((above.max() - above.min()) / abs(above.mean()))
        floor_dev = float(cs_own[4] / above.mean() - 1.0)
        cont = -8.0 * Ic
        devs = cs_own / cont - 1.0
        rich = (4.0 * cs_own[0] - cs_own[2]) / 3.0  # s 12 and s 6: h/s ratio 2
        out[f"law_{kind}"] = {
            "own_E_curv_times_s_eps1": cs_own.tolist(),
            "script": cs_scr.tolist(),
            "spread_above_floor": spread,
            "floor_s3_rel": floor_dev,
            "continuum_-8I": cont,
            "dev_from_continuum": devs.tolist(),
            "richardson_s12_s6": rich,
            "script_outcome": fc[tag]["outcome"],
        }
        gate = 0.11 if kind == "ii" else 0.30
        line(
            f"C8_{kind}_own_lattice_E_curv_times_s_matches_script_1e-8_all_negative_spread_above_floor_within_{int(100 * gate)}_percent_and_s3_low",
            rel(cs_own, cs_scr) < 1e-8
            and np.all(cs_own < 0)
            and spread < gate
            and floor_dev < -0.2,
            f"own {np.round(cs_own, 5).tolist()} vs script {np.round(cs_scr, 5).tolist()} (rel {rel(cs_own, cs_scr):.1e}); spread over s 12..4.5 {100 * spread:.1f} percent; s 3 {100 * floor_dev:.1f} percent from the mean above; script outcome {fc[tag]['outcome']}",
        )
        line(
            f"C9a_{kind}_continuum_-1_over_s_constant_-8I_from_own_quadrature_reproduces_the_s12_lattice_value_within_5_percent_negative_so_the_texture_lowers_E_curv_in_the_continuum",
            abs(devs[0]) < 0.05 and cont < 0,
            f"-8 I = {cont:.5f}; lattice E s / (-8 I) - 1 at s {list(S_LAD)}: {[f'{100 * v:+.2f}' for v in devs]} percent; Richardson (s 12, 6) {rich:.5f}",
        )
        line(
            f"C9b_{kind}_lattice_deviation_from_the_continuum_shrinks_monotonically_as_s_grows",
            bool(np.all(np.diff(np.abs(devs)) > 0)),
            f"|dev| at s 12, 9, 6, 4.5, 3: {[f'{100 * abs(v):.2f}' for v in devs]} percent"
            + (
                ""
                if np.all(np.diff(np.abs(devs)) > 0)
                else ": the lattice sequence crosses the continuum value (stencil artifacts of both signs, cf. C3b), so the 1 percent agreement at s 12 is partly accidental"
            ),
        )
    ii = fc["delta0.3_i_f_eq_q_radial"]
    out["identity_control_labels_in_script"] = {
        "outcome": ii["outcome"],
        "eps1_all_negative": ii["eps1"]["all_negative"],
        "eps1_E_curv_times_s": ii["eps1"]["E_curv_times_s"],
        "E_curv_at_s6": ii["eps4_law_at_s6"],
    }
    line(
        "C10_script_JSON_labels_the_identity_control_by_the_sign_of_roundoff_all_negative_True_on_1e-28_values_and_outcome_TEXTURE_MIXED",
        ii["eps1"]["all_negative"] is True
        and ii["outcome"].startswith("TEXTURE_MIXED")
        and max(abs(v) for v in ii["eps4_law_at_s6"].values()) < 1e-20,
        f"script (i): all_negative {ii['eps1']['all_negative']}, outcome '{ii['outcome'][:14]}', E_curv/eps^4 at s 6 {[f'{v:.1e}' for v in ii['eps4_law_at_s6'].values()]}: a sign read on roundoff; the correct label for (i) is E_curv = 0 (the closed form is exactly 0)",
    )
    # the amplitude condition and the fit's a coefficient at eps 4 (the V4-dominated fit)
    fit_a = {
        e: fc["delta0.3_ii_f_radial_q_pwave"][f"eps{e:g}"]["fit_-a/s+c*s3"]["a"]
        for e in (0.5, 1.0, 2.0, 4.0)
    }
    a_pred = {e: 8.0 * I2 * e**4 for e in (0.5, 1.0, 2.0, 4.0)}
    out["fit_a_vs_closed_form"] = {"script_fit_a": fit_a, "closed_form_8I_eps4": a_pred}
    line(
        "C11_script_fit_coefficient_a_of_-a/s+c_s3_is_NOT_the_curvature_amplitude_8I_eps4_it_is_V4_dominated_and_goes_negative_at_eps_4",
        fit_a[4.0] < 0 and all(abs(fit_a[e] / a_pred[e]) > 10 for e in (0.5, 1.0, 2.0)),
        f"script a: {{{', '.join(f'eps {e:g}: {v:.4g}' for e, v in fit_a.items())}}} vs the closed-form curvature amplitude 8 I eps^4: {{{', '.join(f'eps {e:g}: {v:.4g}' for e, v in a_pred.items())}}}: the -a/s term of the fit absorbs V4 residuals (V4 ~ 1e6 vs E_curv ~ 1)",
    )
    return out


# ================= (d) the two-roots pre-registration =================
def audit_d(FJ, RJ, fres):
    out = {}
    fr = lambda tag: RJ["rows"][tag]["end_reads"]["frame"]["s0.02"]["d2E_ds2"]  # noqa: E731
    c8, c8b, c32 = fr("S1_v4std_n32_g8_x4500"), fr("S1_v4std_n32_g8"), fr("S1_v4std_n32_g32")
    fd = FJ["d_two_roots_preregistration"]
    own = {
        "c_g8_9000": c8,
        "c_g8_4500": c8b,
        "c_g32": c32,
        "ratio": c32 / c8,
        "ratio_4500": c32 / c8b,
        "law_gm1": (31.0 / 7.0) ** 2,
        "law_g2": 16.0,
        "g16_sep_percent": 100.0 * ((15.0 / 7.0) ** 2 / 4.0 - 1.0),
        "own_raw_c_g8": fres["S1_v4std_n32_g8_x4500"]["own_raw_d2E_s0.02"],
        "own_raw_c_g32": fres["S1_v4std_n32_g32"]["own_raw_d2E_s0.02"],
    }
    out.update(own)
    line(
        "D1_R20_frame_curvatures_ratio_19.11_and_the_laws_19.61_16.0_and_the_g16_separation_14.8_percent_match_script_1e-10",
        abs(own["ratio"] - fd["ratio_g32_over_g8"]) < 1e-10
        and abs(own["law_gm1"] - 19.612244897959187) < 1e-12
        and abs(own["law_gm1"] - fd["law_(g-1)^2_ratio"]) < 1e-12
        and fd["law_g^2_ratio"] == 16.0
        and abs(own["g16_sep_percent"] - fd["g16_law_separation_percent"]) < 1e-10
        and abs(own["ratio"] - 19.11) < 0.01,
        f"c(g8, 9000) {c8:.4f}, c(g32) {c32:.3f}, ratio {own['ratio']:.4f}; (31/7)^2 {own['law_gm1']:.4f}, g^2 16; g 16 separation {own['g16_sep_percent']:.3f} percent",
    )
    line(
        "D2_own_raw_boost_curvatures_recomputed_with_own_energy_reproduce_the_R20_frame_numbers_1e-6",
        abs(own["own_raw_c_g8"] - c8) / abs(c8) < 1e-6
        and abs(own["own_raw_c_g32"] - c32) / abs(c32) < 1e-6,
        f"own g8 {own['own_raw_c_g8']:.5f} vs {c8:.5f}; own g32 {own['own_raw_c_g32']:.4f} vs {c32:.4f}",
    )
    mid = 0.5 * (16.0 + own["law_gm1"])
    line(
        "D3_the_ratio_is_convergence_dependent_the_4500_step_row_gives_17.11_which_sits_between_the_laws_while_the_9000_step_row_sits_within_3_percent_of_(g-1)^2",
        abs(own["ratio_4500"] - 17.108) < 0.01
        and 16.0 < own["ratio_4500"] < own["law_gm1"]
        and abs(own["ratio"] / own["law_gm1"] - 1) < 0.03
        and abs(own["ratio_4500"] - mid) < 1.0,
        f"ratio at 4500 steps {own['ratio_4500']:.3f} (midpoint of the laws {mid:.2f}), at 9000 steps {own['ratio']:.3f}: both FALLING rows, the g 8 curvature moved {100 * (c8 / c8b - 1):+.1f} percent between them",
    )
    return out


# ================= (e) the valley slopes =================
def audit_e(FJ):
    out = {}
    nu, D = sp.symbols("nu Delta", real=True)
    q = [sp.Integer(-8), sp.Integer(1), sp.Rational(3, 10), sp.Integer(0)]
    lam = [nu + D, nu, nu, nu]
    C = [sum(qi**p for qi in q) for p in range(1, 5)]
    V4 = sum(
        (sum(l**p for l in lam) - C[p - 1]) ** 2 for p in range(1, 5)
    )  # W1 cancels in the ratio
    P = lambda x: sp.prod([x - qi for qi in q])  # noqa: E731
    Vs = sum(P(l) ** 2 for l in lam)

    def branches(V, D0):
        Vn, Vnn, VnD = sp.diff(V, nu), sp.diff(V, nu, 2), sp.diff(V, nu, D)
        pol = sp.Poly(sp.expand(Vn.subs(D, D0)), nu)
        res = []
        for r in pol.nroots(n=40, maxsteps=200):
            if abs(sp.im(r)) > 1e-20:
                continue
            r = sp.re(r)
            vnn = Vnn.subs({nu: r, D: D0}).evalf(30)
            if vnn > 0:
                slope = float((-VnD.subs({nu: r, D: D0}) / vnn).evalf(30))
                # implicit-function check: re-solve at D0 +- dh from r by Newton on Vn
                dh = sp.Rational(1, 10**4)

                def resolve(Dv):
                    x = sp.Float(r, 40)
                    f_ = sp.lambdify(nu, Vn.subs(D, Dv), "mpmath")
                    fp_ = sp.lambdify(nu, Vnn.subs(D, Dv), "mpmath")
                    for _ in range(60):
                        x = x - f_(x) / fp_(x)
                    return x

                num = float((resolve(D0 + dh) - resolve(D0 - dh)) / (2 * dh))
                res.append({"nu0": float(r), "slope_formula": slope, "slope_implicit_fd": num})
        return sorted(res, key=lambda d: d["nu0"])

    own = {
        "V4_D0": branches(V4, 0),
        "V4_D-9": branches(V4, -9),
        "Vspec_D0": branches(Vs, 0),
        "Vspec_D-9": branches(Vs, -9),
    }
    out["own"] = own
    fe = FJ["e_valley_slope"]
    scr = {
        "V4_D0": fe["V4_D0"],
        "V4_D-9": fe["V4_D-9_hierarchy"],
        "Vspec_D0": fe["Vspec_D0"],
        "Vspec_D-9": fe["Vspec_D-9_hierarchy"],
    }
    worst_fd, worst_scr = 0.0, 0.0
    for k in own:
        for b in own[k]:
            worst_fd = max(worst_fd, abs(b["slope_formula"] - b["slope_implicit_fd"]))
        s_sorted = sorted(scr[k], key=lambda t: t[0])
        if len(s_sorted) != len(own[k]):
            worst_scr = np.inf
        else:
            for b, (n0s, sls) in zip(own[k], s_sorted):
                worst_scr = max(worst_scr, abs(b["nu0"] - n0s), abs(b["slope_formula"] - sls))
    out["worst_formula_vs_implicit_fd"] = worst_fd
    out["worst_vs_script"] = worst_scr
    line(
        "E1_own_sympy_valley_slopes_equal_the_implicit_function_re_solve_1e-6_on_every_branch",
        worst_fd < 1e-6,
        f"worst |formula - re-solve| {worst_fd:.1e} over {sum(len(v) for v in own.values())} branches",
    )
    line(
        "E2_branch_counts_nu0_and_slopes_match_script_1e-8_V4_D0_-1/4_x2_V4_D-9_-1.0054_and_+0.049_Vspec_D0_-1/4_x4",
        worst_scr < 1e-8
        and [round(b["slope_formula"], 6) for b in own["V4_D0"]] == [-0.25, -0.25]
        and abs(own["V4_D-9"][0]["slope_formula"] + 1.0054186) < 1e-6
        and abs(own["V4_D-9"][1]["slope_formula"] - 0.0490591) < 1e-6
        and all(abs(b["slope_formula"] + 0.25) < 1e-12 for b in own["Vspec_D0"])
        and len(own["Vspec_D0"]) == 4,
        "; ".join(
            f"{k}: {[(round(b['nu0'], 5), round(b['slope_formula'], 6)) for b in v]}"
            for k, v in own.items()
        )
        + f"; worst vs script {worst_scr:.1e}",
    )
    # V_spec at Delta -9: the exact valley and the closed form
    Pn = lambda x: np.prod([x - qi for qi in (-8.0, 1.0, 0.3, 0.0)])  # noqa: E731
    Pp = lambda x: sum(
        np.prod([x - qj for j, qj in enumerate((-8.0, 1.0, 0.3, 0.0)) if j != i]) for i in range(4)
    )  # noqa: E731
    closed = -Pp(-8.0) ** 2 / (Pp(-8.0) ** 2 + 3 * Pp(1.0) ** 2)
    out["Vspec_D-9_closed_form"] = {
        "nu0_exact": 1.0,
        "slope": closed,
        "Pprime_-8": Pp(-8.0),
        "Pprime_1": Pp(1.0),
    }
    line(
        "E3_Vspec_D-9_single_branch_at_the_EXACT_valley_-8_1_1_1_with_slope_-Pprime(-8)^2/(Pprime(-8)^2+3Pprime(1)^2)_=_-0.99967_matching_script",
        len(own["Vspec_D-9"]) == 1
        and abs(own["Vspec_D-9"][0]["nu0"] - 1.0) < 1e-12
        and abs(own["Vspec_D-9"][0]["slope_formula"] - closed) < 1e-10
        and abs(fe["Vspec_D-9_hierarchy"][0][1] - closed) < 1e-8
        and abs(Pn(1.0)) == 0.0
        and abs(Pn(-8.0)) == 0.0,
        f"closed form {closed:.10f} (P'(-8) {Pp(-8.0):.1f}, P'(1) {Pp(1.0):.2f}); own branches {[(round(b['nu0'], 8), round(b['slope_formula'], 8)) for b in own['Vspec_D-9']]}; script {fe['Vspec_D-9_hierarchy']}",
    )
    # the symmetric-point identity on a generic permutation-symmetric polynomial with symbolic coefficients
    l = sp.symbols("l1:5", real=True)
    a = sp.symbols("a1:5", real=True)
    c = sp.symbols("c1:5", real=True)
    m = sp.symbols("m1:7", real=True)
    ps = [sum(li**p for li in l) for p in range(1, 5)]
    Vg = sum(a[p] * (ps[p] - c[p]) ** 2 for p in range(4)) + sum(
        m[k] * (ps[i] - c[i]) * (ps[j] - c[j])
        for k, (i, j) in enumerate([(i, j) for i in range(4) for j in range(i + 1, 4)])
    )
    x = sp.symbols("x", real=True)
    sub = {li: x for li in l}
    A_ = sp.diff(Vg, l[0], 2).subs(sub)
    B_ = sp.diff(Vg, l[0], l[1]).subs(sub)
    phi = Vg.subs({l[0]: nu + D, l[1]: nu, l[2]: nu, l[3]: nu})
    phi_nD = sp.diff(phi, nu, D).subs({D: 0, nu: x})
    phi_nn = sp.diff(phi, nu, 2).subs({D: 0, nu: x})
    id1 = sp.simplify(sp.expand(phi_nD - (A_ + 3 * B_)))
    id2 = sp.simplify(sp.expand(phi_nn - (4 * A_ + 12 * B_)))
    # a numeric spot check with random coefficients
    rng = np.random.default_rng(5)
    vals = {s_: float(v) for s_, v in zip(a + c + m, rng.standard_normal(14))}
    ratio_num = float((phi_nD / phi_nn).subs(vals).subs(x, 0.7))
    out["symmetric_point_identity"] = {
        "phi_nD_minus_A_plus_3B": str(id1),
        "phi_nn_minus_4A_plus_12B": str(id2),
        "numeric_ratio_random_coefficients": ratio_num,
        "script": fe["symmetric_point_identity_(A+3B)/(4A+12B)"],
    }
    line(
        "E4_symmetric_point_identity_phi_nD_=_A+3B_and_phi_nn_=_4A+12B_hold_identically_for_a_generic_symmetric_polynomial_so_the_ratio_is_1/4",
        id1 == 0
        and id2 == 0
        and abs(ratio_num - 0.25) < 1e-12
        and fe["symmetric_point_identity_(A+3B)/(4A+12B)"] == "1/4",
        f"identities: {id1}, {id2}; numeric ratio with 14 random coefficients {ratio_num:.15f}; script '{fe['symmetric_point_identity_(A+3B)/(4A+12B)']}'",
    )
    return out


# ================= (f) the norm along the boost =================
def own_Khat(n, h):
    X, Y, Z = own_coords(n, h)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    nx, ny, nz = X / R, Y / R, Z / R
    K = np.zeros(X.shape + (4, 4))
    for i, a in enumerate((nx, ny, nz)):
        K[..., 0, 1 + i] = a
        K[..., 1 + i, 0] = a
    return K, R


def own_Q(Khat, bl):
    K2 = Khat @ Khat
    I = np.broadcast_to(np.eye(4), Khat.shape)
    return I + np.sinh(bl)[..., None, None] * Khat + (np.cosh(bl) - 1.0)[..., None, None] * K2


def audit_f(FJ):
    out = {}
    with open(BSTAR_JSON) as f:
        rec = json.load(f)
    rs, bstar = np.array(rec["rs"]), np.array(rec["b_star"])
    gamma = gamma_r20()
    n, h = 32, 1.5
    Khat, R = own_Khat(n, h)
    bl1 = np.interp(R, rs, bstar)
    Kx = Khat * bl1[..., None, None]
    # the exponential: closed form vs expm on random cells; (n.K)^3 = n.K
    rng = np.random.default_rng(77)
    idx = [tuple(rng.integers(0, n, 3)) for _ in range(40)]
    worst_expm = max(
        float(np.max(np.abs(own_Q(Khat[i], np.array(0.37 * bl1[i]))[()] - expm(0.37 * Kx[i]))))
        for i in idx
    )
    cube_dev = float(np.max(np.abs(Khat @ Khat @ Khat - Khat)))
    out["exponential"] = {"worst_abs_vs_expm_40_cells": worst_expm, "max_abs_K3_minus_K": cube_dev}
    line(
        "F1_own_closed_form_exponential_matches_scipy_expm_1e-13_on_40_cells_and_(n.K)^3_=_n.K",
        worst_expm < 1e-13 and cube_dev < 1e-13,
        f"vs expm {worst_expm:.1e}, |K^3 - K| {cube_dev:.1e}",
    )
    Q05 = own_Q(Khat, 0.05 * bl1)
    lor = float(np.max(np.abs(Q05 @ ETA @ Q05.swapaxes(-1, -2) - ETA)))
    Qs, _ = R3.boost_at(R20.cfg_of(32, 48.0, 8.0), 0.0, 0.05)
    out["Q_vs_stack_boost_at_s0.05"] = float(np.max(np.abs(Qs - Q05)))
    line(
        "F2_Q_eta_QT_=_eta_1e-13_at_s0.05_and_own_Q_equals_the_stack_boost_at_1e-13",
        lor < 1e-13 and out["Q_vs_stack_boost_at_s0.05"] < 1e-13,
        f"|Q eta Q^T - eta| {lor:.1e}; |own Q - boost_at| {out['Q_vs_stack_boost_at_s0.05']:.1e}",
    )
    rows = [
        ("S1_v4std_n32_g8_x4500", "v4std", 8.0),
        ("Sd_v4std_n32_g8_x4500", "v4std", 8.0),
        ("S0_v4std_n32_g8_x4500", "v4std", 8.0),
        ("S1_vspec_n32_g8", "vspec", 8.0),
        ("Sd_vspec_n32_g8", "vspec", 8.0),
        ("S0_vspec_n32_g8", "vspec", 8.0),
        ("S1_v4std_n32_g32", "v4std", 32.0),
    ]
    ff = FJ["f_norm_along_boost"]
    SV = (0.005, 0.01, 0.02, 0.04, 0.05)
    res = {}
    c_poly = poly_coeffs((-8.0, 1.0, 0.3, 0.0))
    c_p1 = np.polyder(c_poly)
    c_p2 = np.polyder(c_p1)
    DIP_ROWS = ("S1_v4std_n32_g8_x4500", "S1_vspec_n32_g8", "S1_v4std_n32_g32")

    def even_lsq(svals, dvals):
        """least-squares fit D(s) = c0 + c1 s^2 (the script's procedure, own implementation)."""
        A_ = np.column_stack([np.ones(len(svals)), np.asarray(svals) ** 2])
        coef, *_ = np.linalg.lstsq(A_, np.asarray(dvals), rcond=None)
        return float(coef[0]), float(coef[1])

    for tag, potk, g in rows:
        cfg = R20.cfg_of(32, 48.0, g)
        p = R20.params_of(g)
        pot = R20.pot_of(potk, cfg, gamma)
        q = (-g, 1.0, DELTA, 0.0)
        h3 = h**3
        M = np.load(os.path.join(R20_NPZ, tag + ".npz"))["M"]
        r = {"potential": potk, "g": g, "max_abs_M0i": float(np.max(np.abs(M[..., 0, 1:])))}
        M2 = M @ M
        n0 = np.einsum("...ii->...", M2)
        d1c = 4.0 * np.einsum("...ii->...", M2 @ Kx)
        d2c = 8.0 * np.einsum("...ii->...", Kx @ Kx @ M2) + 8.0 * np.einsum(
            "...ii->...", Kx @ M @ Kx @ M
        )

        def norm_at(s):
            Q = own_Q(Khat, s * bl1)
            Ms = Q @ M @ Q.swapaxes(-1, -2)
            return np.einsum("...ij,...ji->...", Ms, Ms), Ms

        hh = 0.02  # a second difference at 1e-3 is roundoff-limited near 1e-6 (the script's 4e-6); Richardson on 0.02 / 0.04 removes the O(h^2) term and keeps the roundoff near 1e-9
        nA, MsA = norm_at(hh)
        nB, _ = norm_at(-hh)
        nC, _ = norm_at(2 * hh)
        nD, _ = norm_at(-2 * hh)
        asym = float(np.max(np.abs(MsA - MsA.swapaxes(-1, -2))))
        d1_h, d1_2h = (nA - nB) / (2 * hh), (nC - nD) / (4 * hh)
        d2_h, d2_2h = (nA + nB - 2 * n0) / hh**2, (nC + nD - 2 * n0) / (2 * hh) ** 2
        d1f, d2f = (4 * d1_h - d1_2h) / 3, (4 * d2_h - d2_2h) / 3
        r["norm"] = {
            "d1_closed_max_abs": float(np.max(np.abs(d1c))),
            "d1_richardson_max_abs": float(np.max(np.abs(d1f))),
            "d2_closed_max": float(np.max(d2c)),
            "d2_closed_min": float(np.min(d2c)),
            "d2_richardson_minus_closed_max_abs": float(np.max(np.abs(d2f - d2c))),
            "Ms_asymmetry": asym,
            "script_d2_diff": ff[tag]["d2_max_abs_diff"],
            "script_d2_max": ff[tag]["d2_max_closed"],
        }
        # the spectrum: own eigenvalues and own trace invariants at s 0.05
        Q = own_Q(Khat, 0.05 * bl1)
        Ms = Q @ M @ Q.swapaxes(-1, -2)
        e0 = np.sort(np.linalg.eigvals(M @ ETA).real, axis=-1)
        e1 = np.sort(np.linalg.eigvals(Ms @ ETA).real, axis=-1)
        t0, t1 = own_traces(M), own_traces(Ms)
        r["spectrum_drift_s0.05"] = {
            "eigenvalues_max_abs": float(np.max(np.abs(e1 - e0))),
            "trace_invariants_max_rel": float(
                max(np.max(np.abs(t1[k] - t0[k])) / np.max(np.abs(t0[k])) for k in range(4))
            ),
            "script": ff[tag]["spectrum_drift_s0.05"],
        }
        # energies, own
        ec0, ev0 = own_energy(M, cfg, potk, gamma)
        E0 = float(ec0 + ev0)
        sp_ = R20.energy_parts(M, cfg, p, pot)
        r["E0_own_vs_stack_rel"] = max(
            abs(float(ec0) - sp_["E_curv"]) / sp_["E_curv"],
            abs(float(ev0) - sp_["V"]) / abs(sp_["V"]),
        )

        def path_energies(sv):
            """(E_raw(+s), E_raw(-s), E_n(+s), E_n(-s)) with their (curv, V) parts, own dressing and own rescale."""
            Qp = own_Q(Khat, sv * bl1)
            Qm = own_Q(Khat, -sv * bl1)
            Mp = Qp @ M @ Qp.swapaxes(-1, -2)
            Mm = Qm @ M @ Qm.swapaxes(-1, -2)
            fp = np.sqrt(n0 / np.einsum("...ij,...ji->...", Mp, Mp))
            fm = np.sqrt(n0 / np.einsum("...ij,...ji->...", Mm, Mm))
            parts = [
                own_energy(Mp, cfg, potk, gamma),
                own_energy(Mm, cfg, potk, gamma),
                own_energy(Mp * fp[..., None, None], cfg, potk, gamma),
                own_energy(Mm * fm[..., None, None], cfg, potk, gamma),
            ]
            return [(float(a), float(b)) for a, b in parts], float(min(fp.min(), fm.min()) - 1.0)

        paths = {}
        for sv in SV:
            (rp, rm, npp, npm), fmin = path_energies(sv)
            raw = (sum(rp) + sum(rm) - 2 * E0) / sv**2
            npres = (sum(npp) + sum(npm) - 2 * E0) / sv**2
            key = f"s{sv:g}"
            scr = ff[tag]["paths"].get(key, {})
            paths[key] = {
                "raw": raw,
                "raw_dV_max": float(max(abs(rp[1] - float(ev0)), abs(rm[1] - float(ev0)))),
                "normpres": npres,
                "normpres_curv": (npp[0] + npm[0] - 2 * float(ec0)) / sv**2,
                "normpres_V": (npp[1] + npm[1] - 2 * float(ev0)) / sv**2,
                "dE_raw_plus": sum(rp) - E0,
                "dE_normpres_plus": sum(npp) - E0,
                "dV_normpres_plus": npp[1] - float(ev0),
                "dEcurv_normpres_minus_raw_plus": npp[0] - rp[0],
                "f_minus_1_min": fmin,
                "even_check_E_plus_minus_E_minus": float(abs(sum(rp) - sum(rm))),
                "script_raw": scr.get("raw_d2E_ds2"),
                "script_normpres": scr.get("normpres_d2E_ds2"),
            }
        r["paths"] = paths
        r["own_raw_d2E_s0.02"] = paths["s0.02"]["raw"]
        # (1) the small-step two-point fit D(s) = c0 + c1 s^2 on s 0.005 / 0.01, validated at 0.02 / 0.04 / 0.05
        Dn = {sv: paths[f"s{sv:g}"]["normpres"] for sv in SV}
        Dr = {sv: paths[f"s{sv:g}"]["raw"] for sv in SV}
        c1_n = (Dn[0.01] - Dn[0.005]) / (0.01**2 - 0.005**2)
        c0_n = Dn[0.005] - c1_n * 0.005**2
        c1_r = (Dr[0.01] - Dr[0.005]) / (0.01**2 - 0.005**2)
        c0_r = Dr[0.005] - c1_r * 0.005**2
        valid = {
            f"s{sv:g}": abs(c0_n + c1_n * sv**2 - Dn[sv]) / abs(Dn[sv] - c0_n)
            for sv in (0.02, 0.04, 0.05)
        }
        # (2) the script's procedure: least squares on s 0.01 / 0.02 / 0.05, against its s_to_0 keys
        ls_n = even_lsq([0.01, 0.02, 0.05], [Dn[0.01], Dn[0.02], Dn[0.05]])
        ls_r = even_lsq([0.01, 0.02, 0.05], [Dr[0.01], Dr[0.02], Dr[0.05]])
        s2z = ff[tag].get("s_to_0", {})
        # (3) the closed-form pieces: the residual-gradient shift and the potential's scale stiffness
        _, G, _ = R20.energy_grad(M, cfg, p, pot)
        c2 = -d2c / (4.0 * n0)
        GM = np.einsum("...ab,...ab->...", G, M)
        shift_half = float(np.sum(c2 * GM))
        if potk == "v4std":
            t = own_traces(M)
            cp = [sum(qi**k for qi in q) for k in range(1, 5)]
            HMM = (
                h3
                * W1
                * sum(
                    2.0 * (k * t[k - 1]) ** 2
                    + 2.0 * (t[k - 1] - cp[k - 1]) * k * (k - 1) * t[k - 1]
                    for k in range(1, 5)
                )
            )
        else:
            lam = np.linalg.eigvals(M @ ETA).real
            HMM = (
                gamma
                * h3
                * np.sum(
                    2.0
                    * (
                        (lam * np.polyval(c_p1, lam)) ** 2
                        + np.polyval(c_poly, lam) * lam**2 * np.polyval(c_p2, lam)
                    ),
                    axis=-1,
                )
            )
        k_closed = 0.5 * float(np.sum(c2**2 * HMM))
        # (4) the quartic wall measured directly: the excess of the rescaled path over the raw path, minus the s^2 shift, over s^4
        k_direct = {
            f"s{sv:g}": (
                paths[f"s{sv:g}"]["dE_normpres_plus"]
                - paths[f"s{sv:g}"]["dE_raw_plus"]
                - shift_half * sv**2
            )
            / sv**4
            for sv in (0.01, 0.02, 0.04)
        }
        dV = {sv: paths[f"s{sv:g}"]["dV_normpres_plus"] for sv in SV}
        k_own = k_direct["s0.02"]
        r["quartic_structure"] = {
            "c0_normpres_2pt": c0_n,
            "c1_normpres_2pt": c1_n,
            "c0_raw_2pt": c0_r,
            "c1_raw_2pt": c1_r,
            "two_point_validation_rel_of_rise": valid,
            "c0_normpres_lsq": ls_n[0],
            "c1_normpres_lsq": ls_n[1],
            "c0_raw_lsq": ls_r[0],
            "c1_raw_lsq": ls_r[1],
            "script_s_to_0": s2z,
            "raw_s0.02": Dr[0.02],
            "c0_minus_raw_rel": (c0_n - Dr[0.02]) / abs(Dr[0.02]),
            "shift_2_sum_c2_GM": 2.0 * shift_half,
            "c0_minus_raw_minus_shift_rel": (c0_n - Dr[0.02] - 2.0 * shift_half) / abs(Dr[0.02]),
            "k_from_c1_half": 0.5 * (c1_n - c1_r),
            "k_direct_s4_scaling": k_direct,
            "k_closed_form_half_sum_c2sq_HVMM": k_closed,
            "k_script_c1_over_12": s2z.get("quartic_wall_k"),
            "dV_normpres_plus": {f"s{sv:g}": dV[sv] for sv in SV},
            "dV_ratio_0.02_over_0.01": dV[0.02] / dV[0.01],
            "dV_ratio_0.04_over_0.02": dV[0.04] / dV[0.02],
            "dEcurv_normpres_minus_raw_over_dV_at_s0.02": paths["s0.02"][
                "dEcurv_normpres_minus_raw_plus"
            ]
            / dV[0.02],
            "max_abs_G": float(np.max(np.abs(G))),
        }
        # (5) the dip: E_n(s) - E_0 evaluated directly at the script's s_min and at this file's s_min = sqrt(|c0| / 4k)
        if tag in DIP_ROWS:
            s_own = float(np.sqrt(abs(c0_n) / (4.0 * k_own)))
            E_own = -(c0_n**2) / (16.0 * k_own)
            s_scr = s2z.get("dip_s_min")
            E_scr = s2z.get("dip_E_min")
            grid = [0.5 * s_own, s_own, 1.5 * s_own, 2.0 * s_own] + ([s_scr] if s_scr else [])
            direct = {}
            for sv in grid:
                (rp, rm, npp, npm), _ = path_energies(sv)
                direct[f"{sv:.6f}"] = {
                    "dE_normpres": sum(npp) - E0,
                    "dE_raw": sum(rp) - E0,
                    "dV_normpres": npp[1] - float(ev0),
                }
            r["dip"] = {
                "s_min_own": s_own,
                "E_min_own_formula": E_own,
                "direct_at_s_min_own": direct[f"{s_own:.6f}"]["dE_normpres"],
                "s_min_script": s_scr,
                "E_min_script": E_scr,
                "direct_at_s_min_script": direct[f"{s_scr:.6f}"]["dE_normpres"] if s_scr else None,
                "grid": direct,
            }

        # labels, own, at each s and in the limit
        def label(npv, rawv):
            if np.sign(npv) != np.sign(rawv):
                return "NORM_LIFTS"
            if abs(npv - rawv) <= 0.1 * abs(rawv):
                return "NORM_SILENT"
            return "NORM_LOWERS" if npv < rawv else "NORM_RAISES"

        r["labels_own"] = {f"s{sv:g}": label(Dn[sv], Dr[sv]) for sv in SV}
        r["label_script"] = ff[tag]["outcome"]
        r["label_s_to_0_limit"] = label(c0_n, c0_r)
        res[tag] = r
        qs_ = r["quartic_structure"]
        log(
            f"(f) {tag}: raw {Dr[0.02]:+.3f} (script {paths['s0.02']['script_raw']:+.3f}); normpres s 0.005/0.01/0.02/0.04/0.05 {Dn[0.005]:+.2f} / {Dn[0.01]:+.2f} / {Dn[0.02]:+.2f} / {Dn[0.04]:+.2f} / {Dn[0.05]:+.2f}; "
            f"2pt c0 {c0_n:+.3f} c1 {c1_n:.4g}; lsq c0 {ls_n[0]:+.3f} c1 {ls_n[1]:.4g} (script {s2z.get('c0_normpres')}, {s2z.get('c1_normpres')}); k direct {k_own:.4g} closed {k_closed:.4g} script {s2z.get('quartic_wall_k')}; "
            f"dV ratios {qs_['dV_ratio_0.02_over_0.01']:.2f} / {qs_['dV_ratio_0.04_over_0.02']:.2f}; labels {r['labels_own']}"
            + (
                f"; dip own s {r['dip']['s_min_own']:.4f} E {r['dip']['direct_at_s_min_own']:+.4g}, at script s {r['dip']['s_min_script']:.4f}: {r['dip']['direct_at_s_min_script']:+.4g} (script {r['dip']['E_min_script']:+.4g})"
                if "dip" in r
                else ""
            )
        )
    out["rows"] = res
    tags = [t for t, _, _ in rows]
    w = lambda key: max(res[t]["norm"][key] for t in tags)  # noqa: E731
    line(
        "F3_d/ds_trM2_is_zero_on_every_block_diagonal_end_field_and_the_second_derivative_closed_form_matches_own_Richardson_differences_1e-7",
        w("d1_closed_max_abs") == 0.0
        and w("d1_richardson_max_abs") < 1e-12
        and w("d2_richardson_minus_closed_max_abs") < 1e-7
        and w("Ms_asymmetry") < 1e-13
        and all(res[t]["max_abs_M0i"] == 0.0 for t in tags),
        f"max |d1| closed {w('d1_closed_max_abs'):.1e}, Richardson {w('d1_richardson_max_abs'):.1e}; |d2 Richardson - closed| {w('d2_richardson_minus_closed_max_abs'):.1e} (script's plain-difference {max(res[t]['norm']['script_d2_diff'] for t in tags):.1e}); d2 max {[round(res[t]['norm']['d2_closed_max'], 4) for t in tags]} vs script {[round(res[t]['norm']['script_d2_max'], 4) for t in tags]}",
    )
    line(
        "F4_N_spectrum_invariant_under_the_dressing_eigenvalue_drift_below_1e-12_and_trace_invariants_below_1e-13_relative_at_s0.05",
        max(res[t]["spectrum_drift_s0.05"]["eigenvalues_max_abs"] for t in tags) < 1e-12
        and max(res[t]["spectrum_drift_s0.05"]["trace_invariants_max_rel"] for t in tags) < 1e-13,
        f"eigenvalue drift max {max(res[t]['spectrum_drift_s0.05']['eigenvalues_max_abs'] for t in tags):.1e}, trace-invariant rel max {max(res[t]['spectrum_drift_s0.05']['trace_invariants_max_rel'] for t in tags):.1e}; raw dV max {max(res[t]['paths']['s0.02']['raw_dV_max'] for t in tags):.1e}",
    )
    S3 = (0.01, 0.02, 0.05)
    wr = max(
        abs(res[t]["paths"][f"s{s:g}"]["raw"] - res[t]["paths"][f"s{s:g}"]["script_raw"])
        / abs(res[t]["paths"][f"s{s:g}"]["script_raw"])
        for t in tags
        for s in S3
    )
    wn = max(
        abs(res[t]["paths"][f"s{s:g}"]["normpres"] - res[t]["paths"][f"s{s:g}"]["script_normpres"])
        / abs(res[t]["paths"][f"s{s:g}"]["script_normpres"])
        for t in tags
        for s in S3
    )
    we = max(res[t]["E0_own_vs_stack_rel"] for t in tags)
    line(
        "F5_own_energies_reproduce_the_script_raw_and_norm_preserving_second_differences_1e-6_on_all_7_rows_x_3_steps_and_E0_1e-10",
        wr < 1e-6 and wn < 1e-6 and we < 1e-10,
        f"worst rel raw {wr:.1e}, norm-preserving {wn:.1e}, E0 parts {we:.1e}; own s 0.02: "
        + ", ".join(
            f"{t.split('_n32')[0]} raw {res[t]['paths']['s0.02']['raw']:+.1f} np {res[t]['paths']['s0.02']['normpres']:+.1f}"
            for t in tags
        ),
    )
    expected = {
        "S1_v4std_n32_g8_x4500": "NORM_RAISES",
        "Sd_v4std_n32_g8_x4500": "NORM_LIFTS",
        "S0_v4std_n32_g8_x4500": "NORM_RAISES",
        "S1_vspec_n32_g8": "NORM_SILENT",
        "Sd_vspec_n32_g8": "NORM_SILENT",
        "S0_vspec_n32_g8": "NORM_SILENT",
        "S1_v4std_n32_g32": "NORM_LIFTS",
    }
    line(
        "F6_own_labels_at_s0.02_give_the_brief's_RAISES_LIFTS_RAISES_SILENT_x3_LIFTS_and_the_script_now_labels_every_row_NORM_SILENT_at_s_to_0",
        all(res[t]["labels_own"]["s0.02"] == expected[t] for t in tags)
        and all(res[t]["label_script"].startswith("NORM_SILENT") for t in tags),
        "; ".join(
            f"{t.split('_n32')[0]}: s0.02 {res[t]['labels_own']['s0.02']}, script '{res[t]['label_script']}'"
            for t in tags
        ),
    )
    qs = {t: res[t]["quartic_structure"] for t in tags}
    line(
        "F7_the_norm_preserving_second_difference_is_c0_+_c1_s2_the_s0.02_0.04_0.05_points_predicted_from_s0.005_0.01_within_2_percent_of_their_rise_on_every_row",
        all(v < 0.02 for t in tags for v in qs[t]["two_point_validation_rel_of_rise"].values()),
        "; ".join(
            f"{t.split('_n32')[0]}: c0 {qs[t]['c0_normpres_2pt']:+.2f}, c1 {qs[t]['c1_normpres_2pt']:.3g}, validation {max(qs[t]['two_point_validation_rel_of_rise'].values()):.1e}"
            for t in tags
        ),
    )
    line(
        "F8_the_s_to_0_limit_c0_of_the_norm_preserving_path_equals_the_RAW_curvature_within_1_percent_on_every_row_NORM_SILENT_in_the_limit_matching_the_script's_c0_within_0.5_percent_except_g32",
        all(abs(qs[t]["c0_minus_raw_rel"]) < 0.01 for t in tags)
        and all(res[t]["label_s_to_0_limit"] == "NORM_SILENT" for t in tags)
        and all(
            abs(qs[t]["c0_normpres_2pt"] - qs[t]["script_s_to_0"]["c0_normpres"])
            / abs(qs[t]["c0_normpres_2pt"])
            < 0.005
            for t in tags
            if t != "S1_v4std_n32_g32"
        ),
        "; ".join(
            f"{t.split('_n32')[0]}: c0 {qs[t]['c0_normpres_2pt']:+.2f} vs raw {qs[t]['raw_s0.02']:+.2f} ({100 * qs[t]['c0_minus_raw_rel']:+.2f} percent), script c0 {qs[t]['script_s_to_0'].get('c0_normpres', float('nan')):+.2f}"
            for t in tags
        ),
    )
    line(
        "F9_own_least_squares_on_s0.01_0.02_0.05_reproduces_the_script's_c0_c1_raw_and_normpres_1e-4_and_the_g32_lsq_c0_is_5_percent_off_the_raw_where_the_two_point_small_step_fit_is_within_1_percent",
        all(
            abs(qs[t]["c0_normpres_lsq"] - qs[t]["script_s_to_0"]["c0_normpres"])
            / abs(qs[t]["c0_normpres_lsq"])
            < 1e-4
            and abs(qs[t]["c1_normpres_lsq"] - qs[t]["script_s_to_0"]["c1_normpres"])
            / abs(qs[t]["c1_normpres_lsq"])
            < 1e-4
            and abs(qs[t]["c0_raw_lsq"] - qs[t]["script_s_to_0"]["c0_raw"])
            / abs(qs[t]["c0_raw_lsq"])
            < 1e-4
            for t in tags
        )
        and abs(
            qs["S1_v4std_n32_g32"]["c0_normpres_lsq"] / qs["S1_v4std_n32_g32"]["raw_s0.02"] - 1
        )
        > 0.03
        and abs(qs["S1_v4std_n32_g32"]["c0_minus_raw_rel"]) < 0.01,
        "; ".join(
            f"{t.split('_n32')[0]}: lsq c0 {qs[t]['c0_normpres_lsq']:+.2f} (script {qs[t]['script_s_to_0'].get('c0_normpres', float('nan')):+.2f}), c1 {qs[t]['c1_normpres_lsq']:.4g} (script {qs[t]['script_s_to_0'].get('c1_normpres', float('nan')):.4g})"
            for t in tags
        ),
    )
    line(
        "F10_the_potential's_excess_along_the_rescaled_path_scales_as_s4_dV(0.02)/dV(0.01)_and_dV(0.04)/dV(0.02)_within_5_percent_of_16_on_every_row",
        all(
            abs(qs[t]["dV_ratio_0.02_over_0.01"] - 16) < 0.8
            and abs(qs[t]["dV_ratio_0.04_over_0.02"] - 16) < 0.8
            for t in tags
        ),
        "; ".join(
            f"{t.split('_n32')[0]}: {qs[t]['dV_ratio_0.02_over_0.01']:.2f} / {qs[t]['dV_ratio_0.04_over_0.02']:.2f} (dV at s0.02 {qs[t]['dV_normpres_plus']['s0.02']:.3g}, the curvature part's change is {100 * qs[t]['dEcurv_normpres_minus_raw_over_dV_at_s0.02']:+.2f} percent of it)"
            for t in tags
        ),
    )
    line(
        "F11_the_quartic_wall_k_measured_directly_equals_(c1_normpres_-_c1_raw)/2_within_5_percent_and_the_closed_form_half_sum_c2sq_H_V(M,M)_within_35_percent",
        all(
            abs(qs[t]["k_direct_s4_scaling"]["s0.02"] / qs[t]["k_from_c1_half"] - 1) < 0.05
            and 0.65
            < qs[t]["k_direct_s4_scaling"]["s0.02"] / qs[t]["k_closed_form_half_sum_c2sq_HVMM"]
            < 1.35
            for t in tags
        ),
        "; ".join(
            f"{t.split('_n32')[0]}: k direct {qs[t]['k_direct_s4_scaling']['s0.02']:.4g} (s0.01 {qs[t]['k_direct_s4_scaling']['s0.01']:.4g}, s0.04 {qs[t]['k_direct_s4_scaling']['s0.04']:.4g}), (c1n - c1r)/2 {qs[t]['k_from_c1_half']:.4g}, closed {qs[t]['k_closed_form_half_sum_c2sq_HVMM']:.4g}"
            for t in tags
        ),
    )
    ratio_k = {
        t: qs[t]["k_direct_s4_scaling"]["s0.02"] / qs[t]["k_script_c1_over_12"] for t in tags
    }
    line(
        "F12_script's_quartic_wall_k_=_(c1_normpres_-_c1_raw)/12_equals_the_directly_measured_s4_coefficient_within_10_percent",
        all(abs(v - 1) < 0.1 for v in ratio_k.values()),
        "; ".join(f"{t.split('_n32')[0]}: direct / script = {v:.3f}" for t, v in ratio_k.items())
        + " (a central second difference of k s^4 is 2 k s^2, so k = c1 / 2; the true second derivative 12 k s^2 is not what was fitted)",
    )
    dips = {t: res[t]["dip"] for t in DIP_ROWS}
    line(
        "F13_own_dip_E_n(s_min)_-_E0_=_-c0^2/16k_at_s_min_=_sqrt(|c0|/4k)_is_negative_and_reproduced_by_direct_evaluation_within_15_percent_on_S1_V4_S1_Vspec_g32",
        all(
            d["direct_at_s_min_own"] < 0
            and abs(d["direct_at_s_min_own"] / d["E_min_own_formula"] - 1) < 0.15
            for d in dips.values()
        ),
        "; ".join(
            f"{t.split('_n32')[0]}: s_min {d['s_min_own']:.4f}, formula {d['E_min_own_formula']:+.4g}, direct {d['direct_at_s_min_own']:+.4g}"
            for t, d in dips.items()
        ),
    )
    line(
        "F14_script's_dip_numbers_hold_E_n(s)_-_E0_at_its_s_min_within_15_percent_of_its_E_min",
        all(
            d["direct_at_s_min_script"] is not None
            and abs(d["direct_at_s_min_script"] / d["E_min_script"] - 1) < 0.15
            for d in dips.values()
        ),
        "; ".join(
            f"{t.split('_n32')[0]}: script s_min {d['s_min_script']:.4f} E_min {d['E_min_script']:+.4g}, direct E_n - E0 there {d['direct_at_s_min_script']:+.4g}"
            for t, d in dips.items()
        )
        + "; own s_min = script s_min / sqrt(6) and own E_min = script E_min / 6 if k is 6 times the script's",
    )
    changes = {t: [res[t]["labels_own"][f"s{s:g}"] for s in SV] for t in tags}
    n_change = sum(1 for t in tags if len(set(changes[t])) > 1)
    line(
        "F15_the_finite_step_labels_are_step_artifacts_at_least_4_of_7_rows_change_label_across_s0.005_to_s0.05",
        n_change >= 4,
        "; ".join(f"{t.split('_n32')[0]}: {' / '.join(changes[t])}" for t in tags),
    )
    return out


# ================= (g) the escape-route calibration =================
def surface_cells(n, c0, k):
    lo, hi = c0 - k, c0 + k
    cells = set()
    for i in range(lo, hi + 1):
        for j in range(lo, hi + 1):
            for face in range(6):
                ax, side = face // 2, face % 2
                idx = [i, j]
                idx.insert(ax, hi if side == 0 else lo)
                cells.add(tuple(idx))
    return sorted(cells), lo, hi


def own_degree(vec_field, c0, k, rng, greedy=True):
    """the degree of the oriented line field on the cube surface c0 +- k by preimage counting over 24 random directions;
    orientation over the SURFACE graph only, either a plain BFS or Prim's greedy walk that always extends along the
    strongest link |v . v'| (so near-degenerate cells cannot flip a region); a conflict = an edge whose oriented vectors
    disagree in sign (a nonzero count means the line field is not orientable on the surface and has no integer degree).
    The sign of the degree is a convention (the start cell's sign): |degree| is the invariant."""
    cells, lo, hi = surface_cells(vec_field.shape[0], c0, k)
    cset = set(cells)
    v = {c: vec_field[c] / np.linalg.norm(vec_field[c]) for c in cells}
    sign = {cells[0]: 1.0}

    def nbrs(c):
        for ax in range(3):
            for dd in (-1, 1):
                nb = list(c)
                nb[ax] += dd
                nb = tuple(nb)
                if nb in cset:
                    yield nb

    if greedy:
        heap = []
        for nb in nbrs(cells[0]):
            heapq.heappush(heap, (-abs(float(np.dot(v[cells[0]], v[nb]))), cells[0], nb))
        while heap:
            _, c, nb = heapq.heappop(heap)
            if nb in sign:
                continue
            sign[nb] = sign[c] * (1.0 if np.dot(v[c], v[nb]) >= 0 else -1.0)
            for nb2 in nbrs(nb):
                if nb2 not in sign:
                    heapq.heappush(heap, (-abs(float(np.dot(v[nb], v[nb2]))), nb, nb2))
    else:
        dq = deque([cells[0]])
        while dq:
            c = dq.popleft()
            for nb in nbrs(c):
                if nb not in sign:
                    sign[nb] = sign[c] * (1.0 if np.dot(v[c], v[nb]) >= 0 else -1.0)
                    dq.append(nb)
    vo = {c: sign[c] * v[c] for c in cells}
    conflicts = 0
    for c in cells:
        for ax in range(3):
            nb = list(c)
            nb[ax] += 1
            nb = tuple(nb)
            if nb in cset and np.dot(vo[c], vo[nb]) < 0:
                conflicts += 1
    tris = []
    for ax in range(3):
        for side, plane in ((0, hi), (1, lo)):
            o1, o2 = [a for a in range(3) if a != ax]
            for i in range(lo, hi):
                for j in range(lo, hi):

                    def at(a, b_):
                        idx = [0, 0, 0]
                        idx[ax] = plane
                        idx[o1] = a
                        idx[o2] = b_
                        return vo[tuple(idx)]

                    a_, b_, c_, d_ = at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)
                    # outward orientation: (o1, o2) is right-handed with +ax outward for ax = x (y, z) and ax = z (x, y), left-handed for ax = y (x, z)
                    flip = (side == 1) ^ (ax == 1)
                    for pp, qq, rr in ((a_, b_, c_), (a_, c_, d_)):
                        tris.append((pp, rr, qq) if flip else (pp, qq, rr))
    P_ = np.array([t[0] for t in tris])
    Q_ = np.array([t[1] for t in tris])
    R_ = np.array([t[2] for t in tris])
    degs = []
    for _ in range(24):
        vdir = rng.standard_normal(3)
        vdir /= np.linalg.norm(vdir)
        s1 = np.sign(np.cross(P_, Q_) @ vdir)
        s2 = np.sign(np.cross(Q_, R_) @ vdir)
        s3 = np.sign(np.cross(R_, P_) @ vdir)
        inside = (
            (s1 == s2) & (s2 == s3) & ((P_ + Q_ + R_) @ vdir > 0)
        )  # the same hemisphere: excludes the antipodal triangle, whose three dets are all negative
        degs.append(int(np.sum(s1[inside])))
    return float(np.median(degs)), int(max(degs) - min(degs)), conflicts, len(cells)


def audit_g(FJ):
    out = {}
    gamma = gamma_r20()
    cfg = R20.cfg_of(32, 48.0, 8.0)
    n, h, L = 32, 1.5, 48.0
    h3 = h**3
    X, Y, Z = own_coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rho = np.sqrt(X * X + Y * Y)
    axis = (rho < 0.75 * h) & (np.abs(Z) > 4.0) & (np.abs(Z) < 0.5 * L - 2.0)
    vac = 1.0 + DELTA**2
    q = (-8.0, 1.0, DELTA, 0.0)
    c_poly = poly_coeffs(q)
    rng = np.random.default_rng(909)
    fg = FJ["g_escape_route_calibration"]
    rows = [
        ("Sd_vspec_n32_g8", "Sd", "vspec"),
        ("S0_vspec_n32_g8", "S0", "vspec"),
        ("S1_vspec_n32_g8", "S1", "vspec"),
        ("Sd_v4std_n32_g8_x4500", "Sd", "v4std"),
    ]
    res = {}
    for tag, obj, potk in rows:
        M = np.load(os.path.join(R20_NPZ, tag + ".npz"))["M"]
        M3 = M[..., 1:, 1:]
        lam, vec = np.linalg.eigh(M3)
        nrm = np.sum(lam**2, axis=-1)
        nrm_direct = np.einsum("...ij,...ji->...", M3, M3)
        gaps = np.minimum(lam[..., 1] - lam[..., 0], lam[..., 2] - lam[..., 1])
        deg, spread, conflicts, ncells = own_degree(
            vec[..., :, WRANK[obj]], n // 2, int(round(9.0 / h)), rng, greedy=True
        )
        deg_b, spread_b, conflicts_b, _ = own_degree(
            vec[..., :, WRANK[obj]], n // 2, int(round(9.0 / h)), rng, greedy=False
        )
        seed = R20.seed_axes(cfg, LAM[obj])
        _, vs = np.linalg.eigh(seed[..., 1:, 1:])
        deg_seed, spread_seed, conf_seed, _ = own_degree(
            vs[..., :, WRANK[obj]], n // 2, 6, rng, greedy=True
        )
        # the densities on the melted axis cells (norm below 10 percent of the vacuum), own
        melted = axis & (nrm < 0.1 * vac)
        ev = own_vspec_density(M, q, gamma, h3)
        e4 = own_v4_density(M, q, W1, h3)
        rec = {
            "potential": potk,
            "degree_own": deg,
            "degree_spread_24_dirs": spread,
            "orientation_conflicts": conflicts,
            "surface_cells": ncells,
            "degree_script": fg[tag]["degree_r9_solid_angle"],
            "bfs_orientation": {"degree": deg_b, "spread": spread_b, "conflicts": conflicts_b},
            "orientable_on_surface": conflicts == 0,
            "seed_degree_own": deg_seed,
            "seed_degree_script": fg[tag]["seed"]["degree_r9"],
            "seed_conflicts": conf_seed,
            "axis_cells": int(axis.sum()),
            "axis_norm_min_rel": float(nrm[axis].min() / vac),
            "axis_norm_mean": float(nrm[axis].mean()),
            "axis_crossings_gap_lt_0.02": int(np.sum(gaps[axis] < 0.02)),
            "axis_gap_min": float(gaps[axis].min()),
            "norm_eig_vs_direct_max_abs": float(np.max(np.abs(nrm - nrm_direct))),
            "melted_axis_cells": int(melted.sum()),
            "script": {
                k: fg[tag][k]
                for k in (
                    "axis_norm_min_rel",
                    "axis_norm_mean",
                    "axis_crossings_gap_lt_0.02",
                    "axis_gap_min",
                    "winding_lost",
                    "route",
                )
            },
        }
        if melted.sum() > 0:
            rec["melted_cells"] = {
                "Vspec_density_max": float(ev[melted].max()),
                "V4_density_min": float(e4[melted].min()),
                "Vspec_over_V4_max": float(np.max(ev[melted] / e4[melted])),
                "Vspec_density_field_max": float(ev.max()),
                "lam_at_min_norm": lam[axis][np.argmin(nrm[axis])].tolist(),
            }
        lost = (conflicts > 0) or abs(deg) < 0.5
        rec["winding_lost_own"] = lost
        rec["route_own"] = (
            None
            if not lost
            else (
                "ESCAPE_BY_EXCHANGE"
                if (rec["axis_norm_min_rel"] >= 0.9 and rec["axis_crossings_gap_lt_0.02"] > 0)
                else (
                    "ESCAPE_BY_MELTING"
                    if rec["axis_norm_min_rel"] < 0.7
                    else "ESCAPE_ROUTE_UNDECIDED"
                )
            )
        )
        res[tag] = rec
        log(
            f"(g) {tag}: degree own {deg:+.1f} (spread {spread}, conflicts {conflicts}) script {fg[tag]['degree_r9_solid_angle']:+.3f}; axis norm min/vac {rec['axis_norm_min_rel']:.4f}, crossings {rec['axis_crossings_gap_lt_0.02']}, melted cells {rec['melted_axis_cells']}"
        )
    out["rows"] = res
    tags = [t for t, _, _ in rows]
    orientable = [t for t in tags if res[t]["orientable_on_surface"]]
    line(
        "G1_own_preimage_counting_|degree|_on_the_r9_cube_reproduces_the_script_0_0_1_1_with_zero_spread_where_the_line_field_is_orientable_and_the_seeds_|1|",
        all(
            abs(abs(res[t]["degree_own"]) - abs(round(res[t]["degree_script"]))) < 1e-12
            and res[t]["degree_spread_24_dirs"] == 0
            for t in orientable
        )
        and all(
            abs(abs(res[t]["seed_degree_own"]) - 1.0) < 1e-12 and res[t]["seed_conflicts"] == 0
            for t in tags
        )
        and all(
            abs(res[t]["degree_script"]) < 0.5 for t in tags if not res[t]["orientable_on_surface"]
        )
        and len(orientable) >= 3,
        "; ".join(
            f"{t.split('_n32')[0]}: own {res[t]['degree_own']:+.1f} (greedy conflicts {res[t]['orientation_conflicts']}, BFS {res[t]['bfs_orientation']['conflicts']}, spread {res[t]['degree_spread_24_dirs']}) script {res[t]['degree_script']:+.3f}, seed own {res[t]['seed_degree_own']:+.0f}"
            for t in tags
        )
        + "; the sign is the start cell's convention (opposite to the script's center-rooted tree), |degree| is the invariant",
    )
    line(
        "G1b_the_Sd_Vspec_winding_eigenvector_is_NOT_orientable_on_the_r9_surface_conflicts_under_both_orientations_so_its_degree_is_undefined_rather_than_0",
        res["Sd_vspec_n32_g8"]["orientation_conflicts"] > 0
        and res["Sd_vspec_n32_g8"]["bfs_orientation"]["conflicts"] > 0
        and all(res[t]["orientation_conflicts"] == 0 for t in tags if t != "Sd_vspec_n32_g8"),
        f"Sd V_spec rank-1 eigenvector: greedy conflicts {res['Sd_vspec_n32_g8']['orientation_conflicts']}, BFS conflicts {res['Sd_vspec_n32_g8']['bfs_orientation']['conflicts']}, preimage spread {res['Sd_vspec_n32_g8']['degree_spread_24_dirs']} (the middle eigenvalue's gap min on the axis {res['Sd_vspec_n32_g8']['axis_gap_min']:.4f}); the other three rows orient with 0 conflicts",
    )
    w = max(
        max(
            abs(res[t]["axis_norm_min_rel"] - res[t]["script"]["axis_norm_min_rel"]),
            abs(res[t]["axis_norm_mean"] - res[t]["script"]["axis_norm_mean"]),
            abs(res[t]["axis_gap_min"] - res[t]["script"]["axis_gap_min"]),
        )
        for t in tags
    )
    ok_cross = all(
        res[t]["axis_crossings_gap_lt_0.02"] == res[t]["script"]["axis_crossings_gap_lt_0.02"]
        for t in tags
    )
    line(
        "G2_own_axis_norm_min_0.002_0.000_0.055_0.60_of_1.09_gaps_and_crossings_26_4_58_26_match_script_1e-10",
        w < 1e-10
        and ok_cross
        and max(res[t]["norm_eig_vs_direct_max_abs"] for t in tags) < 1e-12
        and abs(res["Sd_vspec_n32_g8"]["axis_norm_min_rel"] - 0.002) < 5e-4
        and res["S0_vspec_n32_g8"]["axis_norm_min_rel"] < 5e-4
        and abs(res["S1_vspec_n32_g8"]["axis_norm_min_rel"] - 0.055) < 1e-3
        and abs(res["Sd_v4std_n32_g8_x4500"]["axis_norm_min_rel"] - 0.60) < 1e-2,
        "; ".join(
            f"{t.split('_n32')[0]}: min/vac {res[t]['axis_norm_min_rel']:.5f}, crossings {res[t]['axis_crossings_gap_lt_0.02']}, gap min {res[t]['axis_gap_min']:.4f}"
            for t in tags
        )
        + f"; worst |own - script| {w:.1e}; {res[tags[0]]['axis_cells']} axis cells",
    )
    line(
        "G3_route_rule_reproduces_MELTING_MELTING_None_None_and_winding_lost_only_on_the_two_Vspec_escapes",
        all(
            res[t]["route_own"] == res[t]["script"]["route"]
            and res[t]["winding_lost_own"] == res[t]["script"]["winding_lost"]
            for t in tags
        )
        and [res[t]["route_own"] for t in tags]
        == ["ESCAPE_BY_MELTING", "ESCAPE_BY_MELTING", None, None],
        "; ".join(
            f"{t.split('_n32')[0]}: lost {res[t]['winding_lost_own']} -> {res[t]['route_own']}"
            for t in tags
        ),
    )
    P0 = float(np.polyval(c_poly, 0.0))
    m1, m2 = res["Sd_vspec_n32_g8"]["melted_cells"], res["S0_vspec_n32_g8"]["melted_cells"]
    line(
        "G4_P(0)_=_0_so_the_all_zero_spatial_block_is_a_Vspec_vacuum_and_the_melted_axis_cells_carry_a_Vspec_density_below_1_percent_of_their_V4_density",
        P0 == 0.0 and m1["Vspec_over_V4_max"] < 0.01 and m2["Vspec_over_V4_max"] < 0.01,
        f"P(0) = {P0}; Sd: {res['Sd_vspec_n32_g8']['melted_axis_cells']} melted cells, V_spec density max {m1['Vspec_density_max']:.2e} vs V4 density min {m1['V4_density_min']:.2e} (ratio max {m1['Vspec_over_V4_max']:.1e}); S0: {res['S0_vspec_n32_g8']['melted_axis_cells']} cells, {m2['Vspec_density_max']:.2e} vs {m2['V4_density_min']:.2e} (ratio {m2['Vspec_over_V4_max']:.1e}); lam at the Sd minimum {np.round(m1['lam_at_min_norm'], 4).tolist()}",
    )
    line(
        "G5_CALIBRATION_CAVEAT_the_MELTING_norm_criterion_below_70_percent_is_ALSO_met_by_both_fields_that_KEPT_the_winding",
        res["S1_vspec_n32_g8"]["axis_norm_min_rel"] < 0.7
        and res["Sd_v4std_n32_g8_x4500"]["axis_norm_min_rel"] < 0.7
        and abs(res["S1_vspec_n32_g8"]["degree_own"]) == 1.0
        and abs(res["Sd_v4std_n32_g8_x4500"]["degree_own"]) == 1.0,
        f"S1 V_spec: degree {res['S1_vspec_n32_g8']['degree_own']:+.0f} kept, axis norm min {res['S1_vspec_n32_g8']['axis_norm_min_rel']:.3f} of the vacuum; Sd V4: degree {res['Sd_v4std_n32_g8_x4500']['degree_own']:+.0f} kept, {res['Sd_v4std_n32_g8_x4500']['axis_norm_min_rel']:.3f}: the axis norm does not discriminate escape from survival, only the degree does",
    )
    out["vacuum_norm"] = {
        "1_plus_delta2": vac,
        "shell_spectrum_0.36_0.38_0.98_norm": 0.36**2 + 0.38**2 + 0.98**2,
    }
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
    with open(FORM_JSON) as f:
        FJ = json.load(f)
    with open(R20_JSON) as f:
        RJ = json.load(f)
    out = {"task": "M5.32 R21-0 audit", "form_json_wall_s": FJ.get("wall_s")}
    log("(a) the free boundary")
    out["a"] = audit_a(FJ)
    log("(b) the far seed")
    out["b"] = audit_b(FJ, RJ)
    log("(c) the 711.3 texture")
    out["c"] = audit_c(FJ)
    log("(e) the valley slopes")
    out["e"] = audit_e(FJ)
    log("(f) the norm along the boost")
    out["f"] = audit_f(FJ)
    log("(d) the two-roots pre-registration")
    out["d"] = audit_d(FJ, RJ, out["f"]["rows"])
    log("(g) the escape-route calibration")
    out["g"] = audit_g(FJ)
    out["lines"] = LINES
    out["runtime_s"] = time.time() - T0
    with open(OUT, "w") as f:
        json.dump(_jsonable(out), f, indent=1)
    npass = sum(1 for v in LINES.values() if v["pass"])
    print(f"{npass}/{len(LINES)} PASS, runtime {out['runtime_s']:.1f}s, wrote {OUT}")
    return out


if __name__ == "__main__":
    main()
