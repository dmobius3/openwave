"""M5.32 R21-0: the form level of the R21 ladder (ledger section 6.11, amended on the author's
2026-09-14 14:45 UTC reply), one script, seven parts (a) to (g), audited by a fresh agent.

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1), N = M eta; code branch s = -1,
M_vac = diag(8, 1, delta, 0), delta = 0.3 (delta = 1 for the author's vacuum in (c)). The
static energy E = 4 h^3 sum_br wt sum_cells I1(A) + V4, A_i = d_i M on the sym stencil
(forward and backward branches, one-sided at the faces), V4 = W1 sum_p (tr N^p - C_p)^2
(m5_32_r20_1_axes.energy_grad, consumed read-only; its density, m5_32_r20_1_axes.density).
(a) The free-boundary gate. The pinned shell B3.pin_shell (depth 1.6, the outer cells on
    every face) is the only imposed condition in R3 to R20; dropping its mask leaves the same
    energy and the same exact adjoint gradient B3.d1_adj, so the free boundary IS the discrete
    functional's natural boundary. Checked: the packed gradient over ALL cells (the 10 free
    entries per cell, dE/dx_ab = G_ab on the diagonal and 2 G_ab off it) against a 4-point
    stencil in random directions on an n8 box at three step sizes; the per-entry stencil at a
    face, an edge, a corner and an interior cell; the gradient the pin had zeroed on the
    smooth S_d seed (the cells and its size).
(b) The B_far seed. The n64 L96 box's central 32^3 cells coincide with the n32 L48 cells at
    the same h 1.5 ((i - 31.5) h against (i - 15.5) h, the offset 16 cells exact): checked to
    0.0; the restricted field's shell values are the n64 field's bitwise; the restricted
    field's n32 energy and its shell energy against the R20 B_seed row's (S_d: 45.5 of 63.6).
(c) The report's section 711.3 texture on the certified action:
        M_s(x) = M_vac + eps [E_00 f(x / s) + (E_01 + E_10) q(x / s)]
    on the n32 L48 lattice, s in {12, 9, 6, 4.5, 3}, eps in {0.5, 1, 2, 4}, at the delta 0.3
    vacuum and at the author's delta 1 vacuum (-8, 1, 1, 0), with the compact bump
    b(u) = exp(1 - 1 / (1 - |u|^2)) for |u| < 1 (0 outside, sup 1), u = x / s, in three pairs:
        (i)   f = q = b            the identity control: every jet A_i = eps (E_00 + E_01 + E_10)
                                   d_i f is one matrix times a scalar, so F_ij = 0 exactly and
                                   E_curv = 0 to roundoff (labeled TEXTURE_ZERO); in the
                                   continuum the same holds for any two radial profiles
                                   (d_i f d_j q - d_j f d_i q = 0 for f = f(r), q = q(r)), but
                                   on the lattice two different radial profiles keep a
                                   residue (the audit: (b, b^2) gives 27 to 59 percent of the
                                   p-wave E_curv, the one-sided stencil's cross terms): the
                                   texture needs f and q with independent gradients
        (ii)  f = b, q = b u_x     the p-wave pairing (F_ij = eps^2 (d_i f d_j q - d_j f d_i q) F_0)
        (iii) f = b u_x, q = b u_y two p-waves
    The jet identity
    F_0 = A eta B - B eta A = -E_01 + E_10 for A = E_00, B = E_01 + E_10; the record's
    tr(F_0 eta F_0^T eta) = -2 (the transpose is in the record's text) is the stack's pairing
    <F_0, F_0>_eta = sum eta_a eta_b (F_0)_ab^2 = -2; without the transpose the trace is +2,
    and the Frobenius pairing tr(F_0 F_0^T) = +2 is the form-level control (this stack has no
    Frobenius branch). Per (s, eps): E_curv and V4 separately, E_curv s (constant iff the
    -1 / s law), E_curv / eps^4 (constant iff the eps^4 law), V4 / s^3; the fit
    E(s) = -a eps^4 / s + c s^3 (no sigma term on this stack) per eps; the author's amplitude
    condition eps (sup f + 2 sup q) < |g - 1| (eps < 7/3 here) and the locus gate per row.
    TEXTURE_UNBOUNDED iff E_curv < 0 on the texture and |E_curv| s is constant within 30
    percent over the s ladder above the lattice floor (s 12, 9, 6, 4.5, at fixed eps: the
    -1 / s law) at either vacuum; TEXTURE_BOUNDED iff E_curv > 0 on every (s, eps); the
    lattice floor at s = 2h = 3 reported beside the law (its deviation from the mean above).
(d) The two-roots pre-registration from the R20 record alone: the frame curvatures at s 0.02
    on S_1 (g 8, n32, 9000 steps) and S_1 (g 32, n32), their ratio, against (31/7)^2 = 19.61
    for a (g - 1)^2 law and (32/8)^2 = 16.0 for a g^2 law, and the g 16 predictions of each
    (the row R21-3 adds): 225/49 = 4.59 and 4.0 times the g 8 curvature.
(e) The valley slope nu' = -V_nD / V_nn on (nu + Delta, nu, nu, nu): the audited functions of
    m5_32_r21_plan_check.py rerun here (V4 at Delta 0, -9, 8; the separable control -1/4; the
    rank-weighted control -a / (a + 3b) as a right-derivative; the symmetric-point identity).
(f) The norm along the boost on the R20 n32 end fields (S_1, S_d, S_0 under V4 at 9000 steps
    and under V_spec; S_1 at g 32): the dressing M -> Q M Q^T with Q(x) = exp(s b*(r) n.K)
    (m5_32_r3_ii_pair.boost_at, the R20 frame read), pointwise K_x = b*(r) n.K, so the closed
    forms hold per cell: d/ds tr(M^2) = 4 tr(M^2 K_x), d2/ds2 tr(M^2) = 8 tr(K_x^2 M^2)
    + 8 tr(K_x M K_x M); checked against central differences of tr(M_s^2) at s = 1e-3 (max
    over cells); the N-spectrum drift (max over cells and s); then the curvature along the
    norm-preserving path, M_s -> M_s sqrt(tr M^2 / tr M_s^2) per cell, d2E / ds2 at s 0.01 /
    0.02 / 0.05 against the raw boost curvature: NORM_SILENT iff within 10 percent of the raw
    (at s 0.02), NORM_LIFTS iff the sign changes, NORM_LOWERS iff more negative by more than
    10 percent, and NORM_RAISES (added at the read, 2026-09-14 17:05 UTC: the pre-registered
    three cases did not cover it) iff less negative by more than 10 percent with the sign
    kept; the potential's change along the rescaled path reported (it is no longer silent
    there: the per-cell rescale moves the N-spectrum off the vacuum, and V4 pays). Added at
    the read (17:10 UTC, corrected on the audit at 17:45 UTC): the label is taken at s -> 0
    (each path's central-difference curvature D(s) fitted c0 + c1 s^2 on s 0.005 and 0.01,
    validated at 0.02 and 0.05), because the rescale factor is 1 + O(s^2) and the potential's
    cost along the path is a QUARTIC wall k s^4, not a curvature: D(s) of an even path
    E0 + c0 s^2 / 2 + k s^4 is c0 + 2 k s^2, so k = (c1_np - c1_raw) / 2 (the first version
    divided by 12, the audit refuted it by the direct s^4 scaling and the closed form
    (1/2) sum c2^2 H_V(M, M): six times larger k, six times shallower dips); the wall's dip
    (s_min = sqrt(|c0| / 4k), E_min = -c0^2 / 16k, evaluated directly on the path as the
    check) is the bounded remainder of the escape; the s 0.02 label is kept beside it.
(g) The escape-route calibration on the two R20 V_spec escapes (S_d, S_0 under V_spec, n32,
    CONVERGED, the winding lost) and on the electron that kept it (S_1 under V_spec): the
    R21 norm read (m5_32_r21_1_runs.extra_reads: tr(M_3x3^2) along the former polar axis
    against the vacuum's 1.09, the eigenvalue gaps there, the fractions inside r 12), and the
    R21 route rule applied (ESCAPE_BY_EXCHANGE iff the axis norm stays within 10 percent and
    the eigenvalues cross, ESCAPE_BY_MELTING iff it falls below 70 percent, else UNDECIDED),
    so that the vocabulary is calibrated on a known escape before R21-1 reads its own rows.
    Note that tr(M_3x3^2) = sum lam_i^2 depends on the eigenvalues only: every seed carries
    the vacuum's 1.09 outside its core by construction, and a V_spec end field at the shell
    spectrum (0.36, 0.38, 0.98) carries 1.235.

Run: python3 m5_32_r21_0_form.py   (about 5 minutes; the n32 fields from data/m5_32_r20_1/)
Out: ../data/m5_32_r21_0_form.json, ../plots/m5_32_r21_0_texture.png
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
OUT = os.path.join(DATA, "m5_32_r21_0_form.json")
R20_NPZ = os.path.join(DATA, "m5_32_r20_1")
R20_JSON = os.path.join(DATA, "m5_32_r20_1_axes.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    argv = sys.argv
    sys.argv = [argv[0]]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


RUNS = _load("m5_32_r21_1_runs", "m5_32_r21_1_runs.py")
PC = _load("m5_32_r21_plan_check", "m5_32_r21_plan_check.py")
R20, R0, EN, B3, LAG, R3 = RUNS.R20, RUNS.R0, RUNS.EN, RUNS.B3, RUNS.LAG, RUNS.R3
RB = R3.RB
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ================= (a) =================
def part_a():
    cfg = RUNS.cfg_of(8, 12.0, 8.0, 0.3)
    p = RUNS.params_of(8.0, 0.3)
    rng = np.random.default_rng(7)
    M = R20.seed_axes(cfg, RUNS.lam_of("Sd", 0.3))
    Ms = M + 0.05 * B3.sym4(rng.standard_normal(M.shape))
    Ms[..., 0, 1:] = 0.0
    Ms[..., 1:, 0] = 0.0
    out = {
        "box": "n8 L12 (h 1.5)",
        "field": "the S_d seed plus a 0.05 symmetric random perturbation on the spatial block",
    }
    E, G, info = R20.energy_grad(Ms, cfg, p, None)
    out["E"] = float(E)
    for pinned, lab in ((False, "all_cells_faces_included"), (True, "interior_only")):
        mask = RUNS.free_mask(cfg, pinned)
        x0 = RUNS.pack(Ms, mask)

        def fun(x):
            Mx = RUNS.unpack(x, Ms, mask)
            E, G, _ = R20.energy_grad(Mx, cfg, p, None)
            return E, (G[mask][:, RUNS.IU[0], RUNS.IU[1]] * RUNS.OFF).ravel()

        E0, g0 = fun(x0)
        rec = {"cells": int(mask.sum())}
        for e in (1e-3, 1e-4):
            w = 0.0
            for _ in range(4):
                d = rng.standard_normal(x0.shape)
                dd = float(g0 @ d)
                f = [fun(x0 + k * e * d)[0] for k in (-2, -1, 1, 2)]
                fd = (f[0] - 8 * f[1] + 8 * f[2] - f[3]) / (12 * e)
                w = max(w, abs(fd - dd) / abs(dd))
            rec[f"stencil4_rel_e{e:g}"] = w
        out[lab] = rec
    cells = {"face": (0, 3, 4), "edge": (0, 0, 5), "corner": (7, 7, 7), "interior": (3, 4, 4)}
    per = {}
    for kind, (i, j, k) in cells.items():
        worst_abs, worst_rel, gmax = 0.0, 0.0, 0.0
        for a, b in ((1, 1), (1, 2), (2, 3), (0, 0), (3, 3)):
            D = np.zeros_like(Ms)
            D[i, j, k, a, b] = 1.0
            D[i, j, k, b, a] = 1.0
            dd = float(np.sum(G * D))
            e = 1e-3
            f = [
                R20.energy_grad(Ms + s * e * D, cfg, p, None, need_grad=False)[0]
                for s in (-2, -1, 1, 2)
            ]
            fd = (f[0] - 8 * f[1] + 8 * f[2] - f[3]) / (12 * e)
            worst_abs = max(worst_abs, abs(fd - dd))
            worst_rel = max(worst_rel, abs(fd - dd) / max(abs(dd), 1e-3))
            gmax = max(gmax, abs(dd))
        per[kind] = {
            "max_abs_err": worst_abs,
            "max_rel_err_floored_1e-3": worst_rel,
            "max_abs_entry_gradient": gmax,
        }
    out["per_entry_stencil"] = per
    pin = B3.pin_shell(8, 1.5)
    out["pin_on_perturbed_field"] = {
        "pinned_cells": int(pin.sum()),
        "max_abs_G_on_pinned": float(np.max(np.abs(G[pin]))),
        "max_abs_G_interior": float(np.max(np.abs(G[~pin]))),
    }
    E2, G2, _ = R20.energy_grad(M, cfg, p, None)
    out["pin_on_smooth_seed"] = {
        "max_abs_G_on_pinned": float(np.max(np.abs(G2[pin]))),
        "max_abs_G_interior": float(np.max(np.abs(G2[~pin]))),
    }
    cfg32 = RUNS.cfg_of(32, 48.0, 8.0, 0.3)
    pin32 = B3.pin_shell(32, 1.5)
    out["pin_n32"] = {
        "pinned_cells": int(pin32.sum()),
        "of": 32**3,
        "frac": float(pin32.mean()),
        "depth_cells": int(np.ceil(1.6 / 1.5)),
    }
    out["free_boundary_is_the_same_call"] = (
        "energy_grad takes no boundary argument: B_free drops the mask on the descent direction only (m5_32_r21_1_runs.descend, free_mask)"
    )
    log(
        f"(a) all-cells packed stencil {out['all_cells_faces_included']}; per entry {per['face']['max_abs_err']:.1e} (face) {per['corner']['max_abs_err']:.1e} (corner); pin n32 {out['pin_n32']['frac']:.3f} of cells"
    )
    return out


# ================= (b) =================
def part_b():
    out = {}
    X64 = B3.coords(64, 1.5)[0][:, 0, 0]
    X32 = B3.coords(32, 1.5)[0][:, 0, 0]
    out["coord_coincidence_max_abs"] = float(np.max(np.abs(X64[16:48] - X32)))
    out["coord_formula"] = {
        "x64_first_central": float(X64[16]),
        "x32_first": float(X32[0]),
        "offset_cells": 16,
    }
    cfg32 = RUNS.cfg_of(32, 48.0, 8.0, 0.3)
    p = RUNS.params_of(8.0, 0.3)
    with open(R20_JSON) as f:
        J = json.load(f)
    pin = B3.pin_shell(32, 1.5)
    for obj in ("S1", "Sd", "S0"):
        M64 = np.load(os.path.join(R20_NPZ, f"{obj}_v4std_n64_g8.npz"))["M"]
        Mr = RUNS.restrict_center(M64, 32)
        rec = {
            "bitwise_restriction": bool(np.array_equal(Mr, M64[16:48, 16:48, 16:48])),
            "shell_values_equal_n64_bitwise": bool(
                np.array_equal(Mr[pin], M64[16:48, 16:48, 16:48][pin])
            ),
        }
        ep = R20.energy_parts(Mr, cfg32, p, None)
        rad = R20.radial_profile(Mr, cfg32, None)
        rec["restricted_n32"] = {
            "E_total": ep["E_total"],
            "E_curv": ep["E_curv"],
            "V": ep["V"],
            "E_pin_shell": rad["E_pin_shell"],
            "E_lt_12": rad["E_lt_R"]["12"],
        }
        r20 = J["rows"][f"{obj}_v4std_n32_g8_x4500"]["end_reads"]
        rec["R20_B_seed_n32_9000"] = {
            "E_total": r20["energy"]["E_total"],
            "E_pin_shell": r20["radial"]["E_pin_shell"],
            "E_lt_12": r20["radial"]["E_lt_R"]["12"],
        }
        r64 = J["rows"][f"{obj}_v4std_n64_g8"]["end_reads"]
        rec["R20_n64_4500"] = {
            "E_total": r64["energy"]["E_total"],
            "E_lt_12": r64["radial"]["E_lt_R"]["12"],
            "E_lt_24": r64["radial"]["E_lt_R"]["24"],
        }
        # the shell of the restricted field carries the n64 field's r 22.4 to 24 values: its own energy, and the n64 field's energy inside the same 24-box
        c = 16
        e64 = R20.density(M64, RUNS.cfg_of(64, 96.0, 8.0, 0.3), None)
        rec["n64_density_inside_the_central_32_cube"] = float(
            np.sum(e64[c : c + 32, c : c + 32, c : c + 32])
        )
        seed = R20.seed_axes(cfg32, RUNS.lam_of(obj, 0.3))
        rec["seed_n32"] = {
            "E_total": R20.energy_parts(seed, cfg32, p, None)["E_total"],
            "E_pin_shell": R20.radial_profile(seed, cfg32, None)["E_pin_shell"],
        }
        rec["max_abs_diff_restricted_vs_seed_on_shell"] = float(
            np.max(np.abs(Mr[pin] - seed[pin]))
        )
        out[obj] = rec
        log(
            f"(b) {obj}: restricted E {ep['E_total']:.3f} (shell {rad['E_pin_shell']:.3f}) vs B_seed 9000 {rec['R20_B_seed_n32_9000']['E_total']:.3f} (shell {rec['R20_B_seed_n32_9000']['E_pin_shell']:.3f}); n64 density in the central cube {rec['n64_density_inside_the_central_32_cube']:.3f}"
        )
    return out


# ================= (c) =================
def bump(u):
    out = np.zeros_like(u)
    m = np.abs(u) < 1.0
    out[m] = np.exp(1.0 - 1.0 / (1.0 - u[m] ** 2))
    return out


def part_c():
    out = {}
    A = np.zeros((4, 4))
    A[0, 0] = 1.0
    Bm = np.zeros((4, 4))
    Bm[0, 1] = Bm[1, 0] = 1.0
    F0 = A @ ETA @ Bm - Bm @ ETA @ A
    out["jet"] = {
        "F0": F0.tolist(),
        "tr_F0_eta_F0T_eta_the_authors_-2": float(np.trace(F0 @ ETA @ F0.T @ ETA)),
        "tr_F0_eta_F0_eta_no_transpose": float(np.trace(F0 @ ETA @ F0 @ ETA)),
        "tr_F0_F0T_frobenius": float(np.trace(F0 @ F0.T)),
        "stack_inner_eta_F0_F0": float(B3.inner_eta(F0[None], F0[None])[0]),
        "note": "the record's -2 carries a transpose (tr(F_0 eta F_0^T eta), the superscript T in the pasted text); it is the stack's pairing <F, F>_eta = sum eta_a eta_b F_ab^2 = -2 on F_0 = -E_01 + E_10; without the transpose the trace is +2, as is the Frobenius pairing",
    }
    S_LAD = (12.0, 9.0, 6.0, 4.5, 3.0)
    EPS = (0.5, 1.0, 2.0, 4.0)
    out["amplitude_condition_eps_max"] = 7.0 / 3.0
    TEXTURES = {
        "i_f_eq_q_radial": (lambda u, ux, uy: bump(u), lambda u, ux, uy: bump(u)),
        "ii_f_radial_q_pwave": (lambda u, ux, uy: bump(u), lambda u, ux, uy: bump(u) * ux),
        "iii_f_px_q_py": (lambda u, ux, uy: bump(u) * ux, lambda u, ux, uy: bump(u) * uy),
    }
    for delta in (0.3, 1.0):
        for tname, (ff, qq) in TEXTURES.items():
            cfg = RUNS.cfg_of(32, 48.0, 8.0, delta)
            p = RUNS.params_of(8.0, delta)
            X, Y, Z = B3.coords(32, 1.5)
            r = np.sqrt(X * X + Y * Y + Z * Z)
            Mvac = np.zeros((32, 32, 32, 4, 4))
            Mvac[..., 0, 0], Mvac[..., 1, 1], Mvac[..., 2, 2] = 8.0, 1.0, delta
            E_vac = R20.energy_parts(Mvac, cfg, p, None)
            rec = {"vacuum": [8.0, 1.0, delta, 0.0], "texture": tname, "E_vac": E_vac, "rows": []}
            f1 = ff(r / 6.0, X / 6.0, Y / 6.0)
            q1 = qq(r / 6.0, X / 6.0, Y / 6.0)
            rec["sup_f"], rec["sup_q"] = float(np.max(np.abs(f1))), float(np.max(np.abs(q1)))
            rec["amplitude_condition_eps_max"] = 7.0 / (rec["sup_f"] + 2.0 * rec["sup_q"])
            for eps in EPS:
                for s in S_LAD:
                    f = ff(r / s, X / s, Y / s)
                    q = qq(r / s, X / s, Y / s)
                    Ms = Mvac.copy()
                    Ms[..., 0, 0] += eps * f
                    Ms[..., 0, 1] += eps * q
                    Ms[..., 1, 0] += eps * q
                    Ew, Gw, info = R20.energy_grad(Ms, cfg, p, None, need_grad=False)
                    ep = R20.energy_parts(Ms, cfg, p, None)
                    row = {
                        "eps": eps,
                        "s": s,
                        "E_curv": ep["E_curv"],
                        "V4": ep["V"] - E_vac["V"],
                        "E_total_minus_vac": ep["E_total"] - E_vac["E_total"],
                        "E_curv_times_s": ep["E_curv"] * s,
                        "E_curv_over_eps4": ep["E_curv"] / eps**4,
                        "V4_over_s3": (ep["V"] - E_vac["V"]) / s**3,
                        "locus_ok": bool(info["ok"]),
                        "min_gap": float(info["min_gap"]),
                        "cells_in_bump": int(np.sum(r < s)),
                        "max_abs_M0i": float(np.max(np.abs(Ms[..., 0, 1:]))),
                    }
                    rec["rows"].append(row)
                rows = [q for q in rec["rows"] if q["eps"] == eps]
                cs = np.array([q["E_curv_times_s"] for q in rows])
                sv = np.array([q["s"] for q in rows])
                Et = np.array([q["E_total_minus_vac"] for q in rows])
                Aeq = np.stack([-1.0 / sv, sv**3], axis=1)
                coef, *_ = np.linalg.lstsq(Aeq, Et, rcond=None)
                rec[f"eps{eps:g}"] = {
                    "E_curv_times_s": cs.tolist(),
                    "spread_rel": float((cs.max() - cs.min()) / max(abs(cs.mean()), 1e-300)),
                    "all_negative": bool(np.all(np.array([q["E_curv"] for q in rows]) < 0)),
                    "fit_-a/s+c*s3": {
                        "a": float(coef[0]),
                        "c": float(coef[1]),
                        "max_rel_resid": float(
                            np.max(np.abs(Aeq @ coef - Et) / np.maximum(np.abs(Et), 1e-12))
                        ),
                    },
                }
            ec = np.array([q["E_curv"] for q in rec["rows"]])
            law = [rec[f"eps{e:g}"] for e in EPS]
            # the -1/s law is tested on s >= 3h (12, 9, 6, 4.5); the lattice floor s = 2h = 3 is reported beside it (the ledger: "reported either way")
            for e in EPS:
                cs = np.array(
                    [q["E_curv_times_s"] for q in rec["rows"] if q["eps"] == e and q["s"] >= 4.5]
                )
                rec[f"eps{e:g}"]["spread_rel_above_floor"] = float(
                    (cs.max() - cs.min()) / max(abs(cs.mean()), 1e-300)
                )
                cf = [q["E_curv_times_s"] for q in rec["rows"] if q["eps"] == e and q["s"] == 3.0][
                    0
                ]
                rec[f"eps{e:g}"]["floor_s3_rel_to_mean_above"] = float(cf / cs.mean() - 1.0)
            unb = any(l["all_negative"] and l["spread_rel_above_floor"] < 0.3 for l in law)
            if np.max(np.abs(ec)) < 1e-20:
                rec["outcome"] = "TEXTURE_ZERO (F_ij = 0 identically; the identity control)"
            else:
                rec["outcome"] = (
                    "TEXTURE_UNBOUNDED"
                    if unb
                    else (
                        "TEXTURE_BOUNDED"
                        if np.all(ec > 0)
                        else "TEXTURE_MIXED (E_curv negative somewhere, the -1/s law not within 30 percent)"
                    )
                )
            rec["E_curv_at_lattice_floor_s3"] = {
                f"eps{q['eps']:g}": q["E_curv"] for q in rec["rows"] if q["s"] == 3.0
            }
            # the eps^4 law across eps at fixed s
            rec["eps4_law_at_s6"] = {
                f"eps{q['eps']:g}": q["E_curv_over_eps4"] for q in rec["rows"] if q["s"] == 6.0
            }
            out[f"delta{delta:g}_{tname}"] = rec
            log(
                f"(c) delta {delta} {tname}: {rec['outcome']}; E_curv s at eps 1: {[round(v, 4) for v in rec['eps1']['E_curv_times_s']]}; eps^4 law at s 6: {rec['eps4_law_at_s6']}; eps max {rec['amplitude_condition_eps_max']:.2f}"
            )
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, delta in zip(axs, (0.3, 1.0)):
        rec = out[f"delta{delta:g}_ii_f_radial_q_pwave"]
        for eps in EPS:
            rows = [q for q in rec["rows"] if q["eps"] == eps]
            ax.plot(
                [q["s"] for q in rows],
                [q["E_curv"] for q in rows],
                "o-",
                label=f"E_curv, eps {eps:g}",
            )
            ax.plot(
                [q["s"] for q in rows],
                [q["V4"] for q in rows],
                "s--",
                ms=3,
                label=f"V4 - V4_vac, eps {eps:g}",
            )
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xlabel("s (bump radius)")
        ax.set_ylabel("energy")
        ax.set_yscale("symlog", linthresh=1e-2)
        ax.set_title(
            f"section 711.3 texture (ii) on 4 I1 + V4, vacuum (8, 1, {delta:g}, 0): {rec['outcome']}",
            fontsize=8,
        )
        ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r21_0_texture.png"), dpi=130)
    return out


# ================= (d) =================
def part_d():
    with open(R20_JSON) as f:
        J = json.load(f)
    r8 = J["rows"]["S1_v4std_n32_g8_x4500"]["end_reads"]["frame"]
    r8b = J["rows"]["S1_v4std_n32_g8"]["end_reads"]["frame"]
    r32 = J["rows"]["S1_v4std_n32_g32"]["end_reads"]["frame"]
    c8, c8b, c32 = r8["s0.02"]["d2E_ds2"], r8b["s0.02"]["d2E_ds2"], r32["s0.02"]["d2E_ds2"]
    out = {
        "c_g8_9000": c8,
        "c_g8_4500": c8b,
        "c_g32": c32,
        "ratio_g32_over_g8": c32 / c8,
        "ratio_range": [c32 / c8, c32 / c8b],
        "law_(g-1)^2_ratio": (31.0 / 7.0) ** 2,
        "law_g^2_ratio": 16.0,
        "g16_prediction_(g-1)^2_times_c_g8": c8 * (15.0 / 7.0) ** 2,
        "g16_prediction_g^2_times_c_g8": c8 * 4.0,
        "g16_law_separation_percent": 100.0 * ((15.0 / 7.0) ** 2 / 4.0 - 1.0),
        "escape_at_(8,0.3)_-2(delta-1)^2(g-1)^2": -2.0 * (0.3 - 1.0) ** 2 * (8.0 - 1.0) ** 2,
        "pre_registration": "R21-3 fits c (g - g_0)^2 on g in {1, 2, 4, 8, 16, 32} with the root free; FRAME_TWO_ROOTS iff g_0 within 0.3 of 1 and delta_0 within 0.15 of 1 with the product fitting all eight rows to 20 percent",
    }
    log(
        f"(d) ratio g32 / g8 = {c32 / c8:.2f} (4500-step row {c32 / c8b:.2f}); (g-1)^2 19.61, g^2 16.0; the g 16 row separates the laws by {out['g16_law_separation_percent']:.1f} percent"
    )
    return out


# ================= (e) =================
def part_e():
    import sympy as sp

    res = PC.part1()
    res["note"] = (
        "the audited functions of m5_32_r21_plan_check.py rerun (m5_32_r21_plan_check_audit.py, 27 checks); W1 cancels in every nu' (a ratio)"
    )
    # added on the author's 16:50 UTC reply (18437206): the same functions on V_spec = gamma sum_i P(l_i)^2, P(x) = prod (x - q_i),
    # q = (-8, 1, 3/10, 0): separable, so -1/4 at Delta 0 by the theorem; the hierarchy value is the new number (gamma cancels)
    nu, D = sp.symbols("nu Delta", real=True)
    P = lambda x: sp.prod([x - q for q in PC.SPEC_VAC])  # noqa: E731
    Vspec = sum(P(l) ** 2 for l in (nu + D, nu, nu, nu))
    res["Vspec_D0"] = PC.valley_slopes(Vspec, nu, D, 0)
    res["Vspec_D-9_hierarchy"] = PC.valley_slopes(Vspec, nu, D, -9)
    res["Vspec_D8"] = PC.valley_slopes(Vspec, nu, D, 8)
    for k in ("Vspec_D0", "Vspec_D-9_hierarchy", "Vspec_D8"):
        log(f"(e) {k}: {res[k]}")
    return res


# ================= (f) =================
def kx_of(cfg, s):
    """the per-cell generator K_x = b*(r) n.K of the R3 dressing at unit scale, from boost_at's construction."""
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    R = np.sqrt(X * X + Y * Y + Z * Z)
    nx, ny, nz = X / R, Y / R, Z / R
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    rs, bstar = RB.bstar_record()
    bl = np.interp(R.ravel(), rs, bstar).reshape(R.shape)
    return K * bl[..., None, None]


def part_f():
    with open(R20_JSON) as f:
        J = json.load(f)
    gamma = R20.gamma_r20_0()
    tags = [
        ("S1_v4std_n32_g8_x4500", "v4std", 8.0),
        ("Sd_v4std_n32_g8_x4500", "v4std", 8.0),
        ("S0_v4std_n32_g8_x4500", "v4std", 8.0),
        ("S1_vspec_n32_g8", "vspec", 8.0),
        ("Sd_vspec_n32_g8", "vspec", 8.0),
        ("S0_vspec_n32_g8", "vspec", 8.0),
        ("S1_v4std_n32_g32", "v4std", 32.0),
    ]
    out = {}
    for tag, potk, g in tags:
        cfg = R20.cfg_of(32, 48.0, g)
        p = R20.params_of(g)
        pot = R20.pot_of(potk, cfg, gamma)
        M = np.load(os.path.join(R20_NPZ, tag + ".npz"))["M"]
        Kx = kx_of(cfg, 1.0)
        M2 = M @ M
        n0 = np.einsum("...ii->...", M2)
        d1c = 4.0 * np.einsum("...ii->...", M2 @ Kx)
        d2c = 8.0 * np.einsum("...ii->...", Kx @ Kx @ M2) + 8.0 * np.einsum(
            "...ii->...", Kx @ M @ Kx @ M
        )
        hs = 1e-3
        Qp, _ = R3.boost_at(cfg, 0.0, hs)
        Qm, _ = R3.boost_at(cfg, 0.0, -hs)
        Mp, Mm = R3.conj(Qp, M), R3.conj(Qm, M)
        npl = np.einsum("...ij,...ji->...", Mp, Mp)
        nmi = np.einsum("...ij,...ji->...", Mm, Mm)
        d1f = (npl - nmi) / (2 * hs)
        d2f = (npl - 2 * n0 + nmi) / hs**2
        spec = lambda X: np.sort(np.linalg.eigvals(X @ ETA).real, axis=-1)  # noqa: E731
        s0 = spec(M)
        rec = {
            "potential": potk,
            "g": g,
            "norm_vac": float(64.0 + 1.0 + 0.09) if g == 8.0 else float(g * g + 1.09),
            "d1_max_abs_closed": float(np.max(np.abs(d1c))),
            "d1_max_abs_fd": float(np.max(np.abs(d1f))),
            "d1_max_abs_diff": float(np.max(np.abs(d1c - d1f))),
            "d2_max_closed": float(np.max(d2c)),
            "d2_max_abs_diff": float(np.max(np.abs(d2c - d2f))),
            "d2_min_closed": float(np.min(d2c)),
            "max_abs_M0i": float(np.max(np.abs(M[..., 0, 1:]))),
            "spectrum_drift_s0.05": float(
                np.max(np.abs(spec(R3.conj(R3.boost_at(cfg, 0.0, 0.05)[0], M)) - s0))
            ),
        }
        E0 = R20.energy_parts(M, cfg, p, pot)
        rec["E_0"] = E0
        rec["paths"] = {}
        for sv in (0.005, 0.01, 0.02, 0.05):
            Qp, _ = R3.boost_at(cfg, 0.0, sv)
            Qm, _ = R3.boost_at(cfg, 0.0, -sv)
            Mp, Mm = R3.conj(Qp, M), R3.conj(Qm, M)
            Ep, Em = R20.energy_parts(Mp, cfg, p, pot), R20.energy_parts(Mm, cfg, p, pot)
            raw = (Ep["E_total"] + Em["E_total"] - 2 * E0["E_total"]) / sv**2
            fp = np.sqrt(n0 / np.einsum("...ij,...ji->...", Mp, Mp))[..., None, None]
            fm = np.sqrt(n0 / np.einsum("...ij,...ji->...", Mm, Mm))[..., None, None]
            Mpn, Mmn = Mp * fp, Mm * fm
            Epn, Emn = R20.energy_parts(Mpn, cfg, p, pot), R20.energy_parts(Mmn, cfg, p, pot)
            npres = (Epn["E_total"] + Emn["E_total"] - 2 * E0["E_total"]) / sv**2
            rec["paths"][f"s{sv:g}"] = {
                "raw_d2E_ds2": float(raw),
                "raw_dV": float(max(abs(Ep["V"] - E0["V"]), abs(Em["V"] - E0["V"]))),
                "normpres_d2E_ds2": float(npres),
                "normpres_d2Ecurv_ds2": float(
                    (Epn["E_curv"] + Emn["E_curv"] - 2 * E0["E_curv"]) / sv**2
                ),
                "normpres_d2V_ds2": float((Epn["V"] + Emn["V"] - 2 * E0["V"]) / sv**2),
                "rescale_factor_min": float(min(fp.min(), fm.min())),
                "rescale_factor_max": float(max(fp.max(), fm.max())),
                "norm_restored_max_abs_rel": float(
                    np.max(np.abs(np.einsum("...ij,...ji->...", Mpn, Mpn) / n0 - 1.0))
                ),
            }

        def label(np2, raw2):
            if np.sign(np2) != np.sign(raw2):
                return "NORM_LIFTS"
            if abs(np2 - raw2) <= 0.1 * abs(raw2):
                return "NORM_SILENT"
            return (
                "NORM_LOWERS"
                if np2 < raw2
                else "NORM_RAISES (less negative by more than 10 percent, the sign kept)"
            )

        raw2, np2 = rec["paths"]["s0.02"]["raw_d2E_ds2"], rec["paths"]["s0.02"]["normpres_d2E_ds2"]
        rec["outcome_at_s0.02"] = label(np2, raw2)
        # the s -> 0 limit (added at the read, 17:10 UTC): the central-difference curvature of each path fitted c0 + c1 s^2 on
        # s 0.01 / 0.02 / 0.05; the rescale factor deviates from 1 at O(s^2), so the potential's cost along the norm-preserving
        # path is O(s^4): a quartic wall k s^4 with k = (c1_np - c1_raw) / 12, not a curvature; the dip it leaves at
        # s_min = sqrt(|c0| / 4k), E_min = -c0^2 / 16k (for c0 < 0, k > 0)
        # the fitted quantity D(s) = [E(s) + E(-s) - 2 E(0)] / s^2 of an even path E0 + c0 s^2 / 2 + k s^4 is c0 + 2 k s^2, so
        # k = (c1_np - c1_raw) / 2 (the audit refuted the first version's / 12, which would be right for the true second
        # derivative, not for the central difference); c0 from the two smallest steps (the audit: the s 0.05 point biases a
        # least-squares c0 where c1 s^2 is large, the g 32 row); the direct s^4 check: (E_np(s) - E_raw(s)) / s^4 at two s
        sv = np.array([0.005, 0.01])
        A = np.stack([np.ones(2), sv**2], axis=1)
        c_raw = np.array([rec["paths"][f"s{x:g}"]["raw_d2E_ds2"] for x in sv])
        c_np = np.array([rec["paths"][f"s{x:g}"]["normpres_d2E_ds2"] for x in sv])
        cr = np.linalg.solve(A, c_raw)
        cn = np.linalg.solve(A, c_np)
        k = (cn[1] - cr[1]) / 2.0
        s_all = np.array([0.005, 0.01, 0.02, 0.05])
        pred_raw = cr[0] + cr[1] * s_all**2
        pred_np = cn[0] + cn[1] * s_all**2
        excess = {
            f"s{x:g}": float(
                rec["paths"][f"s{x:g}"]["normpres_d2E_ds2"]
                - rec["paths"][f"s{x:g}"]["raw_d2E_ds2"]
            )
            for x in s_all
        }
        k_direct = {f"s{x:g}": float(excess[f"s{x:g}"] / (2.0 * x**2)) for x in s_all}
        rec["s_to_0"] = {
            "c0_raw": float(cr[0]),
            "c1_raw": float(cr[1]),
            "c0_normpres": float(cn[0]),
            "c1_normpres": float(cn[1]),
            "fit_points_s": [0.005, 0.01],
            "quartic_wall_k": float(k),
            "k_direct_from_excess_over_2s2": k_direct,
            "validation_at_s0.02_s0.05_abs_dev_raw": [
                float(abs(pred_raw[2] - rec["paths"]["s0.02"]["raw_d2E_ds2"])),
                float(abs(pred_raw[3] - rec["paths"]["s0.05"]["raw_d2E_ds2"])),
            ],
            "validation_at_s0.02_s0.05_abs_dev_normpres": [
                float(abs(pred_np[2] - rec["paths"]["s0.02"]["normpres_d2E_ds2"])),
                float(abs(pred_np[3] - rec["paths"]["s0.05"]["normpres_d2E_ds2"])),
            ],
        }
        if cn[0] < 0 and k > 0:
            smin = float(np.sqrt(-cn[0] / (4.0 * k)))
            rec["s_to_0"]["dip_s_min"] = smin
            rec["s_to_0"]["dip_E_min"] = float(-cn[0] ** 2 / (16.0 * k))
            # the dip evaluated directly on the path at s_min (the audit's check: the first version's dips were positive there)
            Qd, _ = R3.boost_at(cfg, 0.0, smin)
            Md = R3.conj(Qd, M)
            Md = Md * np.sqrt(n0 / np.einsum("...ij,...ji->...", Md, Md))[..., None, None]
            rec["s_to_0"]["dip_E_direct_at_s_min"] = float(
                R20.energy_parts(Md, cfg, p, pot)["E_total"] - E0["E_total"]
            )
        rec["outcome"] = label(cn[0], cr[0]) + " (at s -> 0)"
        out[tag] = rec
        log(
            f"(f) {tag}: d2 closed vs fd {rec['d2_max_abs_diff']:.1e}; at s 0.02 raw {raw2:+.1f} norm-preserving {np2:+.1f} (V {rec['paths']['s0.02']['normpres_d2V_ds2']:+.1f}): {rec['outcome_at_s0.02']}; s -> 0: raw {cr[0]:+.1f} norm-preserving {cn[0]:+.1f}, wall k {k:.3g}, dip {rec['s_to_0'].get('dip_E_min', float('nan')):.3g} at s {rec['s_to_0'].get('dip_s_min', float('nan')):.3g}: {rec['outcome']}"
        )
    return out


# ================= (g) =================
def part_g():
    gamma = R20.gamma_r20_0()
    out = {}
    cfg = R20.cfg_of(32, 48.0, 8.0)
    for tag, obj in (
        ("Sd_vspec_n32_g8", "Sd"),
        ("S0_vspec_n32_g8", "S0"),
        ("S1_vspec_n32_g8", "S1"),
        ("Sd_v4std_n32_g8_x4500", "Sd"),
    ):
        potk = "vspec" if "vspec" in tag else "v4std"
        pot = R20.pot_of(potk, cfg, gamma)
        M = np.load(os.path.join(R20_NPZ, tag + ".npz"))["M"]
        ex = RUNS.extra_reads(M, cfg, pot, True)
        seed = R20.seed_axes(cfg, RUNS.lam_of(obj, 0.3))
        exs = RUNS.extra_reads(seed, cfg, pot, True)
        wr = RUNS.WINDING_RANK[obj]
        deg = ex["solid_angle"][f"rank{wr}"]["degree"]["9"]
        nr = ex["norm3"]
        lost = abs(deg) < 0.5
        route = None
        if lost:
            route = (
                "ESCAPE_BY_EXCHANGE"
                if (nr["axis_norm_min_rel"] >= 0.9 and nr["axis_crossings_gap_lt_0.02"] > 0)
                else (
                    "ESCAPE_BY_MELTING"
                    if nr["axis_norm_min_rel"] < 0.7
                    else "ESCAPE_ROUTE_UNDECIDED"
                )
            )
        prof = nr["axis_profile"]
        out[tag] = {
            "potential": potk,
            "degree_r9_solid_angle": deg,
            "winding_lost": lost,
            "route": route,
            "axis_norm_min_rel": nr["axis_norm_min_rel"],
            "axis_norm_mean": nr["axis_norm_mean"],
            "axis_crossings_gap_lt_0.02": nr["axis_crossings_gap_lt_0.02"],
            "axis_gap_min": nr["axis_gap_min"],
            "shells": nr["shells"],
            "frac_within_10pct_r_lt_12": nr["frac_within_10pct_r_lt_12"],
            "frac_below_70pct_r_lt_12": nr["frac_below_70pct_r_lt_12"],
            "core_norm_center": nr["core_norm_center"],
            "axis_profile_norm": [q["norm"] for q in prof],
            "axis_profile_z": [q["z"] for q in prof],
            "axis_profile_lam": [q["lam"] for q in prof],
            "seed": {
                "axis_norm_min_rel": exs["norm3"]["axis_norm_min_rel"],
                "axis_gap_min": exs["norm3"]["axis_gap_min"],
                "degree_r9": exs["solid_angle"][f"rank{wr}"]["degree"]["9"],
                "core_norm_center": exs["norm3"]["core_norm_center"],
            },
        }
        log(
            f"(g) {tag}: degree r9 {deg:+.3f} lost {lost} -> {route}; axis norm min / vac {nr['axis_norm_min_rel']:.3f}, mean {nr['axis_norm_mean']:.3f}, crossings {nr['axis_crossings_gap_lt_0.02']}, gap min {nr['axis_gap_min']:.3f}"
        )
    out["calibration_note"] = (
        "the vocabulary as pre-registered: EXCHANGE needs the axis norm within 10 percent of 1.09 with a crossing (gap < 0.02); "
        "MELTING needs it below 70 percent; the V_spec escapes sit at the shell spectrum (0.36, 0.38, 0.98), norm 1.235, 13 percent above the vacuum"
    )
    return out


if __name__ == "__main__":
    os.makedirs(PLOTS, exist_ok=True)
    out = {
        "task": "M5.32 R21-0",
        "a_free_boundary_gate": part_a(),
        "b_far_seed": part_b(),
        "c_texture_711_3": part_c(),
        "d_two_roots_preregistration": part_d(),
        "e_valley_slope": part_e(),
        "f_norm_along_boost": part_f(),
        "g_escape_route_calibration": part_g(),
        "wall_s": round(time.time() - T0, 1),
    }
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    log(f"wrote {os.path.relpath(OUT)} in {out['wall_s']} s")
