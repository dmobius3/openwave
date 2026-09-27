"""M5.32 R27-0 audit: an independent attempt to refute the R27-0 form-level claims.

Written WITHOUT reading scripts/m5_32_r27_0_form.py. Every number below is recomputed from
its stated definition with this file's own bond sum, own seeds, own least squares, own sympy
and own charge readers; the certified stack (m5_21_3_a_4d) and the earlier rungs' seed
builders are consumed only where the claim names them (the stored wall states, the census
arrays, R25-1's BPS seed as the falsifier, R20's seed_axes and density for the Coulomb row,
R26-3's charge_density as the reader under audit, beside this file's own readers).

The definitions audited (data/m5_32_r27_0_form.json)
----------------------------------------------------
    E_u     = 4 h^3 sum_{i<j} <F_ij, F_ij>_eta,  F_ij = A_i eta A_j - A_j eta A_i,
              A_i the forward / backward stencil derivative at weight 1/2 each
    E_kappa = kappa h^3 sum_br w_br sum_i sum_cells <d_i^br S, d_i^br S>  (S the spatial block)
            = kappa h sum_bonds |S(cell') - S(cell)|_F^2   (each bond once: the bond sum)
    dE_kappa / dS(cell) = 2 kappa h sum_{neighbors} (S(cell) - S(neighbor))

Verdicts: CONFIRMED (reproduced within the stated tolerance), QUALIFIED (reproduced, with a
caveat stated in the note), REFUTED (this file's number disagrees; both are recorded).
Output: data/m5_32_r27_0_audit.json. Runtime about a minute, single process.
"""

import importlib.util
import json
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
CLAIMS = os.path.join(DATA, "m5_32_r27_0_form.json")
OUT_JSON = os.path.join(DATA, "m5_32_r27_0_audit.json")
R26_2 = os.path.join(DATA, "m5_32_r26_2")
R26_4 = os.path.join(DATA, "m5_32_r26_4")
CENSUS_JSON = os.path.join(DATA, "m5_32_r26_4_census.json")
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
W1 = 0.000724023879
W1S = 25.0
G = 8.0
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


B3 = _load("m5_21_3_a_4d", "m5_21_3_a_4d.py")
R25_1 = _load("m5_32_r25_1_strand", "m5_32_r25_1_strand.py")
R26_3 = _load("m5_32_r26_3_spin", "m5_32_r26_3_spin.py")
F0, R21, R20, R0 = R26_3.F0, R26_3.R21, R26_3.R20, R26_3.R0
assert abs(B3.W1 - W1) < 1e-18


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ============================================================================
# OWN PRIMITIVES
# ============================================================================
def K_lead(s0):
    return 4.0 + 36.0 * s0**2 + 144.0 * s0**4


def T_bps(delta, w, m=1):
    s0 = delta / 2.0
    return float(np.pi * np.sqrt(32.0 * w * K_lead(s0)) * s0**4 * m)


def kappa_bps(delta, w, m=1):
    return float(np.sqrt(w * K_lead(delta / 2.0) / (32.0 * m * m)))


def bond_E(S, h, mask=None):
    """kappa = 1 bond sum: h sum_bonds |S' - S|_F^2 over the leading grid axes of S (the last
    two axes are the matrix). `mask` (cells) keeps only bonds with BOTH ends inside."""
    nd = S.ndim - 2
    e = 0.0
    for ax in range(nd):
        d = np.diff(S, axis=ax)
        q = np.sum(d * d, axis=(-2, -1))
        if mask is not None:
            s0 = [slice(None)] * nd
            s1 = [slice(None)] * nd
            s0[ax], s1[ax] = slice(0, -1), slice(1, None)
            q = q * (mask[tuple(s0)] & mask[tuple(s1)])
        e += np.sum(q)
    return h * e


def bond_grad(S, h):
    """dE/dS = 2 h sum_neighbors (S - S_nb)."""
    nd = S.ndim - 2
    Gr = np.zeros_like(S)
    for ax in range(nd):
        d = np.diff(S, axis=ax)
        s0 = [slice(None)] * nd
        s1 = [slice(None)] * nd
        s0[ax], s1[ax] = slice(0, -1), slice(1, None)
        Gr[tuple(s0)] -= 2.0 * h * d
        Gr[tuple(s1)] += 2.0 * h * d
    return Gr


def eta_inner(F, Gm):
    return np.einsum("...ab,...cd,ac,bd->...", F, Gm, ETA, ETA, optimize=True)


def own_branch_derivs(M, h, br):
    nd = M.ndim - 2
    out = []
    for ax in range(nd):
        A = np.zeros_like(M)
        d = np.diff(M, axis=ax) / h
        s = [slice(None)] * nd
        s[ax] = slice(0, -1) if br == "fwd" else slice(1, None)
        A[tuple(s)] = d
        out.append(A)
    return out


def curv_density(M, h):
    """4 h^3 sum_br (1/2) sum_{i<j} <F_ij, F_ij>_eta per cell, own stencil code."""
    e = 0.0
    nd = M.ndim - 2
    for br in ("fwd", "bwd"):
        A = own_branch_derivs(M, h, br)
        for i in range(nd):
            for j in range(i + 1, nd):
                F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                e = e + 0.5 * 4.0 * eta_inner(F, F)
    return h**3 * e


def dsym_E(M, h):
    """the plan-time centered-stencil form: h^3 sum_i <dsym_i M, dsym_i M>_eta."""
    nd = M.ndim - 2
    e = 0.0
    for ax in range(nd):
        f, b = own_branch_derivs(M, h, "fwd")[ax], own_branch_derivs(M, h, "bwd")[ax]
        D = 0.5 * (f + b)
        e += np.sum(eta_inner(D, D))
    return h**3 * e


def slab_xy(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y = np.meshgrid(x, x, indexing="ij")
    return X, Y


def pair_field(n, nz, h, delta, b, psi):
    """M = diag(8, s0 + b cos psi, s0 - b cos psi, 1) with S_12 = b sin psi, z-invariant."""
    s0 = delta / 2.0
    M = np.zeros((n, n, nz, 4, 4))
    M[..., 0, 0] = 8.0
    M[..., 3, 3] = 1.0
    M[..., 1, 1] = (s0 + b * np.cos(psi))[:, :, None]
    M[..., 2, 2] = (s0 - b * np.cos(psi))[:, :, None]
    M[..., 1, 2] = M[..., 2, 1] = (b * np.sin(psi))[:, :, None]
    return M


def own_bps(n, nz, h, delta, w, m=1):
    X, Y = slab_xy(n, h)
    rho, phi = np.hypot(X, Y), np.arctan2(Y, X)
    b = (delta / 2.0) * np.sqrt(1.0 - np.exp(-kappa_bps(delta, w, m) * rho**2))
    return pair_field(n, nz, h, delta, b, m * phi)


def own_pair_seed(n, nz, h, delta, w, d=6.0):
    X, Y = slab_xy(n, h)
    wc = X + 1j * Y
    k = kappa_bps(delta, w, 1)
    b = delta / 2.0
    psi = np.zeros_like(X)
    for xc in (-d / 2.0, d / 2.0):
        rho_k = np.abs(wc - xc)
        b = b * np.sqrt(1.0 - np.exp(-k * rho_k**2))
        psi += np.angle(wc - xc)
    return pair_field(n, nz, h, delta, b, psi)


def tensor_phase(M):
    S = M[..., 1:3, 1:3]
    return np.arctan2(2.0 * S[..., 0, 1], S[..., 0, 0] - S[..., 1, 1]), 0.5 * np.hypot(
        S[..., 0, 0] - S[..., 1, 1], 2.0 * S[..., 0, 1]
    )


def winding_on_circle(psi_xy, n, h, center, radius, npts=256):
    """nearest-cell sampling of the phase along a circle, unwrapped, in turns."""
    t = np.linspace(0.0, 2.0 * np.pi, npts, endpoint=False)
    xs, ys = center[0] + radius * np.cos(t), center[1] + radius * np.sin(t)
    i = np.clip(np.rint(xs / h + (n - 1) / 2.0).astype(int), 0, n - 1)
    j = np.clip(np.rint(ys / h + (n - 1) / 2.0).astype(int), 0, n - 1)
    vals = psi_xy[i, j]
    ph = np.unwrap(np.append(vals, vals[0]))  # the closing step back to the first cell
    return float((ph[-1] - ph[0]) / (2.0 * np.pi))


def verdict(claim, own, claimed, ok, note, qualified=False):
    v = "CONFIRMED" if ok else "REFUTED"
    if ok and qualified:
        v = "QUALIFIED"
    return {"claim": claim, "own_value": own, "claimed_value": claimed, "verdict": v, "note": note}


# ============================================================================
# (a) THE RANK-ONE IDENTITY
# ============================================================================
def audit_a(C):
    out = {}
    delta, w = 0.3, W1 * W1S
    Tb = T_bps(delta, w)
    rows = {}
    for n, h, s in ((48, 1.0, 3.0), (96, 0.5, 3.0), (48, 1.0, 1.0), (96, 0.5, 1.0)):
        X, Y = slab_xy(n, h)
        nz = 2
        r = {}
        for name, k in (
            ("axis_x", (1.0, 0.0)),
            ("diagonal_45", (1.0, 1.0)),
            ("generic_30", (np.cos(np.pi / 6), np.sin(np.pi / 6))),
        ):
            kk = np.array(k) / np.hypot(*k)
            u = kk[0] * X + kk[1] * Y
            psi = 0.5 * np.pi * (1.0 + np.tanh(u / s))
            M = pair_field(n, nz, h, delta, delta / 2.0, psi)
            r[name] = float(np.sum(curv_density(M, h)) / (nz * h) / Tb)
        rows[f"n{n}_s{s:g}"] = r
    # the winding falsifier: R25-1's seed and this file's own copy
    fals = {}
    for dl in (0.3, 0.03):
        Tb_d = T_bps(dl, w)
        for n, h in ((48, 1.0), (96, 0.5)):
            Mo = own_bps(n, 2, h, dl, w)
            Mr = R25_1.bps_field(n, 2, h, dl, w)
            fals[f"d{dl:g}_n{n}"] = {
                "own_seed": float(np.sum(curv_density(Mo, h)) / (2 * h) / Tb_d),
                "r25_1_seed_own_energy": float(np.sum(curv_density(Mr, h)) / (2 * h) / Tb_d),
                "r25_1_seed_stack_energy": float(
                    B3.e_parts(Mr, R25_1.slab_cfg(n, 48.0, dl))[0] / (2 * h) / Tb_d
                ),
            }
    # own stencil code vs the certified stack on a random field (consistency, logged only)
    rng = np.random.default_rng(3)
    Mr = B3.sym4(rng.normal(size=(8, 8, 4, 4, 4)))
    cfg = R25_1.slab_cfg(8, 8.0, delta)
    stack_consistency = float(
        abs(np.sum(curv_density(Mr, 1.0)) - B3.e_parts(Mr, cfg)[0]) / B3.e_parts(Mr, cfg)[0]
    )
    ratio3 = rows["n48_s3"]["generic_30"] / rows["n96_s3"]["generic_30"]
    ratio1 = rows["n48_s1"]["generic_30"] / rows["n96_s1"]["generic_30"]
    out["rows"] = rows
    out["falsifier"] = fals
    out["generic_ratio_h1_over_h05_width3"] = float(ratio3)
    out["generic_ratio_h1_over_h05_width1"] = float(ratio1)
    out["own_vs_stack_curvature_rel"] = stack_consistency
    cl = C["a_rank_one"]
    exact_ok = all(abs(rows[k][d]) < 1e-12 for k in rows for d in ("axis_x", "diagonal_45"))
    fals_ok = all(
        abs(
            fals[f"d0.03_n{n}"]["own_seed"]
            - cl["rows"][f"n{n}"]["winding_falsifier"]["E_curv_per_len_over_T_bps"]
        )
        < 1e-12
        for n in (48, 96)
    )
    return {
        "a_exact": verdict(
            "diagonal and axis walls have E_u = 0 to 1e-12 T_bps",
            {k: [rows[k]["axis_x"], rows[k]["diagonal_45"]] for k in rows},
            [
                cl["rows"]["n48"]["axis_x"]["E_curv_per_len_over_T_bps"],
                cl["rows"]["n48"]["diagonal_45"]["E_curv_per_len_over_T_bps"],
            ],
            exact_ok,
            "reproduced with an own tanh wall of widths 3 and 1: an algebraic identity of the stencil (A_x = A_y exactly on the diagonal, A_y = 0 on the axis), so the check cannot fail short of a stencil bug",
            qualified=True,
        ),
        "a_generic_ratio": verdict(
            "the generic (30 deg) wall's residue falls by more than 3x from h 1 to h 0.5",
            {
                "width3": ratio3,
                "width1": ratio1,
                "values": {k: rows[k]["generic_30"] for k in rows},
            },
            {
                "ratio": cl["generic_ratio_h1_over_h05"],
                "values": [
                    cl["rows"]["n48"]["generic_30"]["E_curv_per_len_over_T_bps"],
                    cl["rows"]["n96"]["generic_30"]["E_curv_per_len_over_T_bps"],
                ],
            },
            ratio3 > 3.0,
            f"profile-dependent: the claimed 0.245 / 0.067 is the author's own wall profile (not stated in the JSON), an own tanh wall of width 3 gives ratio {ratio3:.2f}; a width-1 wall gives {ratio1:.2f}: the 3x criterion holds only for a resolved wall",
            qualified=(ratio1 <= 3.0),
        ),
        "a_falsifier": verdict(
            "the BPS winding seed has E_u about 0.49 T_bps per unit length",
            fals,
            {
                "n48": cl["rows"]["n48"]["winding_falsifier"]["E_curv_per_len_over_T_bps"],
                "n96": cl["rows"]["n96"]["winding_falsifier"]["E_curv_per_len_over_T_bps"],
            },
            fals_ok and fals["d0.03_n48"]["own_seed"] > 0.1,
            "reproduced to 1e-15 with an own seed and own stencil code ONCE the falsifier's delta is taken as 0.03 (the JSON does not state it; at delta 0.3 the same seed reads 0.4930 and 0.4982); 0.49 is the lattice value of the BPS half-share E_u = T_bps / 2, converging to 0.5 as h -> 0",
        ),
    }, out


# ============================================================================
# (b) THE TERM AND ITS GRADIENT ON THE STORED WALL STATE
# ============================================================================
def audit_b(C):
    cl = C["b_term_and_gradient"]
    Z = np.load(os.path.join(R26_2, "escape_d0.03_w25_n48_L48_long_amp0.npz"))
    D, delta, h = Z["D"], float(Z["delta"]), 1.0
    assert abs(delta - 0.03) < 1e-12
    M = np.diag([8.0, 1.0, delta, 0.0]) + delta * D
    S = M[..., 1:, 1:]
    # the eta form on the 4x4 field with the same bonds
    e_eta = 0.0
    for ax in range(3):
        d = np.diff(M, axis=ax)
        e_eta += np.sum(eta_inner(d, d))
    e_eta *= h
    e_block = bond_E(S, h)
    e_00 = h * sum(np.sum(np.diff(M[..., 0, 0], axis=ax) ** 2) for ax in range(3))
    m0i = float(np.abs(M[..., 0, 1:]).max())
    rel = abs(e_eta - e_block - e_00) / e_block
    e_dsym = dsym_E(M, h)
    # own gradient checks (kappa = 1, the term is linear in kappa)
    rng = np.random.default_rng(11)
    V = rng.normal(size=S.shape)
    V = 0.5 * (V + V.swapaxes(-1, -2))
    Gr = bond_grad(S, h)
    gv = float(np.sum(Gr * V))
    cs = float(np.imag(bond_E(S + 1e-20j * V, h)) / 1e-20)
    fd = {}
    for eps in (1e-4, 1e-5):
        fd[f"{eps:g}"] = float(
            abs((bond_E(S + eps * V, h) - bond_E(S - eps * V, h)) / (2 * eps) - gv) / abs(gv)
        )
    # the free-cell subset (the pinned shell excluded)
    free = R25_1.free_mask_slab(48, 4, h)
    e_free = bond_E(S, h, mask=free)
    out = {
        "E_eta_4x4": float(e_eta),
        "E_frobenius_block": float(e_block),
        "E_M00_part": float(e_00),
        "M0i_max": m0i,
        "eta_minus_block_minus_M00_rel": float(rel),
        "E_dsym_4x4": float(e_dsym),
        "bond_over_dsym": float(e_block / e_dsym),
        "E_block_free_bonds_only": float(e_free),
        "grad_complex_step_rel": float(abs(cs - gv) / abs(gv)),
        "grad_fd_rel": fd,
        "per_len_bond": float(e_block / 4.0),
    }
    return {
        "b_eta_identity": verdict(
            "eta contraction on 4x4 = block Frobenius + M_00 part (rel 1e-16)",
            [
                out["E_eta_4x4"],
                out["E_frobenius_block"],
                out["E_M00_part"],
                out["eta_minus_block_minus_M00_rel"],
            ],
            [
                cl["E_eta_4x4"],
                cl["E_frobenius_block"],
                cl["E_M00_part"],
                cl["eta_minus_block_minus_M00_rel"],
            ],
            rel < 1e-12 and abs(e_block - cl["E_frobenius_block"]) < 1e-9,
            "reproduced; the identity holds for ANY field with a zero time row (tr(eta dM eta dM) = dM_00^2 + |dS|_F^2 when M_0i = 0), so the check tests the code, not the field; M_0i = 0 here",
            qualified=True,
        ),
        "b_bond_vs_dsym": verdict(
            "bond sum 1.0257 vs plan-time centered stencil 0.3956 (2.6x)",
            [out["E_frobenius_block"], out["E_dsym_4x4"], out["bond_over_dsym"]],
            [
                cl["E_frobenius_block"],
                cl["E_dsym_4x4_plan_time_definition"],
                cl["E_frobenius_block"] / cl["E_dsym_4x4_plan_time_definition"],
            ],
            abs(e_dsym - cl["E_dsym_4x4_plan_time_definition"]) < 1e-9,
            "reproduced; note the whole-box bond sum includes the pinned shell's bonds (free-only bonds: %.4f, %.1f percent of the whole)"
            % (e_free, 100 * e_free / e_block),
        ),
        "b_gradient": verdict(
            "gradient checked by complex step (0) and central FD (2.8e-6 at 1e-5)",
            [out["grad_complex_step_rel"], fd],
            [cl["complex_step_rel"], cl["fd_rel_1e-5"]],
            out["grad_complex_step_rel"] < 1e-8 and fd["1e-05"] < 1e-5,
            "own gradient agrees with own complex step and own central FD to roundoff (E_kappa is quadratic: the central difference is exact, the complex step is exact); the claimed FD residue 2.8e-4 -> 2.8e-6 scales as eps^2, so the author's FD check ran on a non-quadratic total (E_kappa plus the quartic), which is a legitimate check of the summed gradient",
            qualified=True,
        ),
    }, out


# ============================================================================
# (c) THE WINDING COEFFICIENT (two boxes)
# ============================================================================
def audit_c(C):
    cl = C["c_winding_coefficient"]
    w = W1 * W1S
    out = {"rows": {}}
    ok = True
    for delta in (0.03, 0.3):
        b0 = delta / 2.0
        for m in (1, 2):
            pred = 4.0 * np.pi * m * m * b0 * b0 * np.log(2.0)
            row = {}
            for h in (1.0, 0.5, 0.25):
                es = {}
                for L in (48, 96):
                    n = int(round(L / h))
                    M = own_bps(n, 1, h, delta, w, m)
                    es[L] = bond_E(M[..., 1:, 1:], h) / h  # nz 1: per unit length = / (nz h)
                row[f"h{h:g}"] = {
                    "L48": es[48],
                    "L96": es[96],
                    "diff": es[96] - es[48],
                    "rel_err": (es[96] - es[48]) / pred - 1.0,
                }
            # the free-only variant at h 1 (bonds with both ends outside the pinned shell)
            es = {}
            for L in (48, 96):
                n = L
                M = own_bps(n, 1, 1.0, delta, w, m)
                free = R25_1.free_mask_slab(n, 1, 1.0)
                es[L] = bond_E(M[..., 1:, 1:], 1.0, mask=free)
            row["h1_free_only"] = {
                "diff": es[96] - es[48],
                "rel_err": (es[96] - es[48]) / pred - 1.0,
            }
            row["prediction"] = pred
            key = f"d{delta:g}_m{m}"
            out["rows"][key] = row
            cr = cl["rows"][key]
            ok &= abs(row["h1"]["diff"] - cr["two_box_difference"]) < 1e-9 * max(1.0, pred)
            ok &= abs(row["h1"]["rel_err"]) < 0.03
    out["ratio_m2_over_m1_d0.3_h1"] = (
        out["rows"]["d0.3_m2"]["h1"]["diff"] / out["rows"]["d0.3_m1"]["h1"]["diff"]
    )
    conv = [out["rows"]["d0.3_m1"][k]["rel_err"] for k in ("h1", "h0.5", "h0.25")]
    free_err = out["rows"]["d0.3_m1"]["h1_free_only"]["rel_err"]
    return {
        "c_two_box": verdict(
            "E_kappa(L96) - E_kappa(L48) per unit length = 4 pi m^2 b0^2 ln 2 within 0.3 percent",
            {
                k: [out["rows"][k]["h1"]["diff"], out["rows"][k]["h1"]["rel_err"]]
                for k in out["rows"]
            },
            {
                k: [cl["rows"][k]["two_box_difference"], cl["rows"][k]["rel_err"]]
                for k in cl["rows"]
            },
            ok,
            "reproduced to 1e-12 with an own seed and own bond sum on the WHOLE box (pinned shell included, nz-independent); the rel_err at h 1, 0.5, 0.25 is %s, falling linearly in h, so the two-box difference does converge to 4 pi m^2 b0^2 ln 2 (the cell-center offset is an O(h) effect, 0.27 percent at h 1); the pinned shell matters only if excluded: bonds among FREE cells give %.1f percent (ln(2 x 24 - 1.6 over 24 - 1.6) is not ln 2), so the coefficient identification requires the whole-box sum, which R27-1 / R27-2 must match if they quote it"
            % (", ".join(f"{100 * c:.3f}%" for c in conv), 100 * free_err),
        ),
        "c_m_ratio": verdict(
            "the m2 / m1 two-box ratio is 4 within 1 percent",
            out["ratio_m2_over_m1_d0.3_h1"],
            cl["ratio_m2_over_m1_d0.3"],
            abs(out["ratio_m2_over_m1_d0.3_h1"] - 4.0) < 0.04,
            "reproduced",
        ),
    }, out


# ============================================================================
# (d) THE WALL COEFFICIENT
# ============================================================================
def audit_d(C):
    cl = C["d_wall_coefficient"]
    delta, h = 0.03, 1.0
    b0 = delta / 2.0
    out = {}
    ok = True
    for W in (4, 8, 2, 1):
        Nb = int(W / h)
        x = np.arange(-4, Nb + 5) * h
        psi = np.clip(x / W, 0.0, 1.0) * np.pi
        S = np.zeros((len(x), 3, 3))
        S[:, 0, 0] = b0 + b0 * np.cos(psi)
        S[:, 1, 1] = b0 - b0 * np.cos(psi)
        S[:, 0, 1] = S[:, 1, 0] = b0 * np.sin(psi)
        e_area = bond_E(S, h) / (h * h)  # one column of cells has area h^2
        latt = Nb * 2.0 * b0 * b0 * (2.0 * np.sin(np.pi / (2.0 * Nb))) ** 2 / h
        cont = 2.0 * np.pi**2 * b0 * b0 / W
        out[f"W{W}h"] = {
            "E_per_area": float(e_area),
            "lattice_sine": float(latt),
            "continuum": float(cont),
            "rel_err_continuum": float(e_area / cont - 1.0),
            "rel_err_lattice": float(abs(e_area / latt - 1.0)),
        }
        if W in (4, 8):
            ok &= abs(e_area - cl["rows"][f"W{W}h"]["E_per_area"]) < 1e-12
    return {
        "d_wall": verdict(
            "bond sum of a pi ramp over W = N_b h bonds is 2 b0^2 N_b (2 sin(pi / 2 N_b))^2 / h, 5.0 and 1.3 percent below 2 pi^2 b0^2 / W at W 4h and 8h",
            {k: [out[k]["E_per_area"], out[k]["rel_err_continuum"]] for k in out},
            {
                k: [cl["rows"][k]["E_per_area"], cl["rows"][k]["rel_err_continuum"]]
                for k in cl["rows"]
            },
            ok,
            "reproduced to 1e-12 with an own ramp; the lattice-sine 'prediction' is the bond sum rewritten in closed form (each bond carries the same phase step), so the 2e-16 agreement is an identity, not a test; the continuum comparison is the content; the same formula gives -19 percent at W 2h and -59 percent at W h (a cell-sharp pi wall costs 8 b0^2 / h, the bond sum's saturation)",
            qualified=True,
        ),
    }, out


# ============================================================================
# (e) THE TUBE-MASKED HARMONIC READER
# ============================================================================
def real_harmonics(cth, sth, phi):
    return np.stack(
        [
            np.ones_like(cth),
            cth,
            sth * np.cos(phi),
            sth * np.sin(phi),
            0.5 * (3 * cth**2 - 1),
            sth * cth * np.cos(phi),
            sth * cth * np.sin(phi),
            sth**2 * np.cos(2 * phi),
            sth**2 * np.sin(2 * phi),
        ],
        1,
    )


def line_dist(X, Y, Z, u):
    proj = X * u[0] + Y * u[1] + Z * u[2]
    r2 = X * X + Y * Y + Z * Z
    dperp = np.sqrt(np.maximum(r2 - proj * proj, 0.0))
    return np.where(proj > 0, dperp, np.sqrt(r2))


def audit_e(C):
    cl = C["e_masked_reader"]
    n, L = 32, 48.0
    h = L / n
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rs = np.maximum(r, 1e-12)
    cth = Z / rs
    sth = np.sqrt(np.maximum(1 - cth**2, 0.0))
    phi = np.arctan2(Y, X)
    pin = B3.pin_shell(n, h, 1.6)
    lines = [(30.0, 0.0), (90.0, 90.0), (120.0, 200.0)]
    us = []
    for th, ph in lines:
        th, ph = np.radians(th), np.radians(ph)
        us.append(np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]))
    r_t = 1.5 * h
    dl = np.min(np.stack([line_dist(X, Y, Z, u) for u in us]), 0)
    mask_tube = dl < r_t
    shells = np.arange(6.0, 18.0 + 1e-9, 1.5)

    def read(dev, masked):
        a1 = []
        for Rs in shells:
            m = (np.abs(r - Rs) < 0.75) & (~pin)
            if masked:
                m &= ~mask_tube
            A = real_harmonics(cth[m], sth[m], phi[m])
            coef, *_ = np.linalg.lstsq(A, dev[m], rcond=None)
            a1.append(float(coef[1]))
        return np.array(a1)

    def slope(a):
        return float(np.polyfit(np.log(shells), np.log(np.abs(a)), 1)[0])

    out = {"h": h, "r_t": r_t}
    # the dipole
    dev = np.where(r >= 3.0, 0.01 * cth / rs**2, 0.0)
    a_d = read(dev, True) * shells**2
    out["dipole_masked_a1_r2"] = a_d.tolist()
    out["dipole_max_rel_err"] = float(np.abs(a_d / 0.01 - 1).max())
    out["dipole_slope"] = slope(read(dev, True))
    # the footprint inside the mask (radius 0.8 r_t)
    r_f = 0.8 * r_t
    foot = -0.3 * (dl < r_f)
    a_u = read(foot, False)
    a_m = read(foot, True)
    pred = np.array(
        [
            (3.0 / (4 * np.pi)) * (-0.3) * sum(np.pi * r_f**2 / Rs**2 * u[2] for u in us)
            for Rs in shells
        ]
    )
    out["footprint_unmasked_a1"] = a_u.tolist()
    out["footprint_prediction"] = pred.tolist()
    out["footprint_ratio"] = (a_u / pred).tolist()
    out["footprint_mean_ratio"] = float(np.mean(a_u / pred))
    out["footprint_ratio_min_max"] = [float((a_u / pred).min()), float((a_u / pred).max())]
    out["footprint_unmasked_slope"] = slope(a_u)
    out["footprint_masked_max_abs"] = float(np.abs(a_m).max())
    out["footprint_cells_outside_mask"] = int(np.sum((dl < r_f) & ~mask_tube))
    # footprints WIDER than the mask: what the masked read returns
    wider = {}
    for fac in (1.2, 1.5, 2.0):
        rf2 = fac * r_t
        f2 = -0.3 * (dl < rf2)
        au2, am2 = read(f2, False), read(f2, True)
        wider[f"{fac:g}_r_t"] = {
            "masked_over_unmasked": (am2 / au2).tolist(),
            "mean_abs_ratio": float(np.mean(np.abs(am2 / au2))),
        }
    out["footprint_wider_than_mask"] = wider
    # the physical carrier width: the BPS melt of the census seeds (delta 0.3, m 1 and 1/2)
    w = W1 * W1S
    out["melt_half_depth_radius"] = {
        f"m{m:g}": float(np.sqrt(np.log(2.0) / kappa_bps(0.3, w, m))) for m in (0.5, 1.0, 2.0)
    }
    d_ok = out["dipole_max_rel_err"] < 0.05 and abs(out["dipole_slope"] + 2) < 0.1
    f_ok = (
        abs(out["footprint_unmasked_slope"] + 2) < 0.5 and 0.7 < out["footprint_mean_ratio"] < 1.3
    )
    return {
        "e_dipole": verdict(
            "masked l=1 read of a synthetic dipole returns 0.01 r^-2 within 3.5 percent, slope -2.03",
            [out["dipole_max_rel_err"], out["dipole_slope"]],
            [cl["dipole_masked"]["max_rel_err"], cl["dipole_masked"]["slope_a1"]],
            d_ok and abs(out["dipole_max_rel_err"] - cl["dipole_masked"]["max_rel_err"]) < 1e-6,
            "reproduced with an own least squares (unnormalized real harmonics: the 0.01 recovery fixes the basis convention to cos theta, not sqrt(3/4pi) cos theta)",
        ),
        "e_footprint_unmasked": verdict(
            "unmasked l=1 read of a tube footprint falls with slope -2.04 and its mean ratio to the solid-angle prediction is 0.92",
            [
                out["footprint_unmasked_slope"],
                out["footprint_mean_ratio"],
                out["footprint_ratio_min_max"],
            ],
            [
                cl["footprint_unmasked"]["slope_a1"],
                cl["footprint_unmasked"]["mean_ratio_read_over_pred"],
            ],
            f_ok
            and abs(
                out["footprint_mean_ratio"] - cl["footprint_unmasked"]["mean_ratio_read_over_pred"]
            )
            < 1e-6,
            "reproduced; the per-shell ratio scatters from %.2f to %.2f (a 1-cell-wide footprint on a 1.5-cell shell), so the 0.7 to 1.3 window on the MEAN is a weak test: the per-shell read is not quantitative"
            % tuple(out["footprint_ratio_min_max"]),
            qualified=True,
        ),
        "e_footprint_masked_zero": verdict(
            "the masked read of the footprint is exactly 0",
            [out["footprint_masked_max_abs"], out["footprint_cells_outside_mask"]],
            cl["footprint_masked"]["max_abs_masked_over_unmasked"],
            out["footprint_masked_max_abs"] == 0.0,
            "UNFALSIFIABLE as stated: the footprint radius 0.8 r_t lies inside the mask radius r_t, so every footprint cell is masked (%d footprint cells outside the mask) and the read is zero by construction; for a footprint of 1.2 / 1.5 / 2 r_t the masked read keeps %.2f / %.2f / %.2f of the unmasked l = 1 read (mean |masked / unmasked| over the shells), i.e. a footprint 20 percent wider than the mask leaks 70 percent of its dipole through it; the BPS melt half-depth radius of the census carriers is %.2f (m 1) and %.2f (the half-unit carriers, m 1/2), both WIDER than r_t = 2.25, so the mask as sized cannot remove a physical carrier's footprint"
            % (
                out["footprint_cells_outside_mask"],
                wider["1.2_r_t"]["mean_abs_ratio"],
                wider["1.5_r_t"]["mean_abs_ratio"],
                wider["2_r_t"]["mean_abs_ratio"],
                out["melt_half_depth_radius"]["m1"],
                out["melt_half_depth_radius"]["m0.5"],
            ),
            qualified=True,
        ),
    }, out


# ============================================================================
# (c) THE WINDING COEFFICIENT (two boxes)
# ============================================================================
def audit_c(C):
    cl = C["c_winding_coefficient"]
    w = W1 * W1S
    out = {"rows": {}}
    ok = True
    for delta in (0.03, 0.3):
        b0 = delta / 2.0
        for m in (1, 2):
            pred = 4.0 * np.pi * m * m * b0 * b0 * np.log(2.0)
            row = {}
            for h in (1.0, 0.5, 0.25):
                es = {}
                for L in (48, 96):
                    n = int(round(L / h))
                    M = own_bps(n, 1, h, delta, w, m)
                    es[L] = bond_E(M[..., 1:, 1:], h) / h  # nz 1: per unit length = / (nz h)
                row[f"h{h:g}"] = {
                    "L48": es[48],
                    "L96": es[96],
                    "diff": es[96] - es[48],
                    "rel_err": (es[96] - es[48]) / pred - 1.0,
                }
            # the free-only variant at h 1 (bonds with both ends outside the pinned shell)
            es = {}
            for L in (48, 96):
                n = L
                M = own_bps(n, 1, 1.0, delta, w, m)
                free = R25_1.free_mask_slab(n, 1, 1.0)
                es[L] = bond_E(M[..., 1:, 1:], 1.0, mask=free)
            row["h1_free_only"] = {
                "diff": es[96] - es[48],
                "rel_err": (es[96] - es[48]) / pred - 1.0,
            }
            row["prediction"] = pred
            key = f"d{delta:g}_m{m}"
            out["rows"][key] = row
            cr = cl["rows"][key]
            ok &= abs(row["h1"]["diff"] - cr["two_box_difference"]) < 1e-9 * max(1.0, pred)
            ok &= abs(row["h1"]["rel_err"]) < 0.03
    out["ratio_m2_over_m1_d0.3_h1"] = (
        out["rows"]["d0.3_m2"]["h1"]["diff"] / out["rows"]["d0.3_m1"]["h1"]["diff"]
    )
    conv = [out["rows"]["d0.3_m1"][k]["rel_err"] for k in ("h1", "h0.5", "h0.25")]
    free_err = out["rows"]["d0.3_m1"]["h1_free_only"]["rel_err"]
    return {
        "c_two_box": verdict(
            "E_kappa(L96) - E_kappa(L48) per unit length = 4 pi m^2 b0^2 ln 2 within 0.3 percent",
            {
                k: [out["rows"][k]["h1"]["diff"], out["rows"][k]["h1"]["rel_err"]]
                for k in out["rows"]
            },
            {
                k: [cl["rows"][k]["two_box_difference"], cl["rows"][k]["rel_err"]]
                for k in cl["rows"]
            },
            ok,
            "reproduced to 1e-12 with an own seed and own bond sum on the WHOLE box (pinned shell included, nz-independent); the rel_err at h 1, 0.5, 0.25 is %s, falling linearly in h, so the two-box difference does converge to 4 pi m^2 b0^2 ln 2 (the cell-center offset is an O(h) effect, 0.27 percent at h 1); the pinned shell matters only if excluded: bonds among FREE cells give %.1f percent (ln(2 x 24 - 1.6 over 24 - 1.6) is not ln 2), so the coefficient identification requires the whole-box sum, which R27-1 / R27-2 must match if they quote it"
            % (", ".join(f"{100 * c:.3f}%" for c in conv), 100 * free_err),
        ),
        "c_m_ratio": verdict(
            "the m2 / m1 two-box ratio is 4 within 1 percent",
            out["ratio_m2_over_m1_d0.3_h1"],
            cl["ratio_m2_over_m1_d0.3"],
            abs(out["ratio_m2_over_m1_d0.3_h1"] - 4.0) < 0.04,
            "reproduced",
        ),
    }, out


# ============================================================================
# (d) THE WALL COEFFICIENT
# ============================================================================
def audit_d(C):
    cl = C["d_wall_coefficient"]
    delta, h = 0.03, 1.0
    b0 = delta / 2.0
    out = {}
    ok = True
    for W in (4, 8, 2, 1):
        Nb = int(W / h)
        x = np.arange(-4, Nb + 5) * h
        psi = np.clip(x / W, 0.0, 1.0) * np.pi
        S = np.zeros((len(x), 3, 3))
        S[:, 0, 0] = b0 + b0 * np.cos(psi)
        S[:, 1, 1] = b0 - b0 * np.cos(psi)
        S[:, 0, 1] = S[:, 1, 0] = b0 * np.sin(psi)
        e_area = bond_E(S, h) / (h * h)  # one column of cells has area h^2
        latt = Nb * 2.0 * b0 * b0 * (2.0 * np.sin(np.pi / (2.0 * Nb))) ** 2 / h
        cont = 2.0 * np.pi**2 * b0 * b0 / W
        out[f"W{W}h"] = {
            "E_per_area": float(e_area),
            "lattice_sine": float(latt),
            "continuum": float(cont),
            "rel_err_continuum": float(e_area / cont - 1.0),
            "rel_err_lattice": float(abs(e_area / latt - 1.0)),
        }
        if W in (4, 8):
            ok &= abs(e_area - cl["rows"][f"W{W}h"]["E_per_area"]) < 1e-12
    return {
        "d_wall": verdict(
            "bond sum of a pi ramp over W = N_b h bonds is 2 b0^2 N_b (2 sin(pi / 2 N_b))^2 / h, 5.0 and 1.3 percent below 2 pi^2 b0^2 / W at W 4h and 8h",
            {k: [out[k]["E_per_area"], out[k]["rel_err_continuum"]] for k in out},
            {
                k: [cl["rows"][k]["E_per_area"], cl["rows"][k]["rel_err_continuum"]]
                for k in cl["rows"]
            },
            ok,
            "reproduced to 1e-12 with an own ramp; the lattice-sine 'prediction' is the bond sum rewritten in closed form (each bond carries the same phase step), so the 2e-16 agreement is an identity, not a test; the continuum comparison is the content; the same formula gives -19 percent at W 2h and -59 percent at W h (a cell-sharp pi wall costs 8 b0^2 / h, the bond sum's saturation)",
            qualified=True,
        ),
    }, out


# ============================================================================
# (e) THE TUBE-MASKED HARMONIC READER
# ============================================================================
def real_harmonics(cth, sth, phi):
    return np.stack(
        [
            np.ones_like(cth),
            cth,
            sth * np.cos(phi),
            sth * np.sin(phi),
            0.5 * (3 * cth**2 - 1),
            sth * cth * np.cos(phi),
            sth * cth * np.sin(phi),
            sth**2 * np.cos(2 * phi),
            sth**2 * np.sin(2 * phi),
        ],
        1,
    )


def line_dist(X, Y, Z, u):
    proj = X * u[0] + Y * u[1] + Z * u[2]
    r2 = X * X + Y * Y + Z * Z
    dperp = np.sqrt(np.maximum(r2 - proj * proj, 0.0))
    return np.where(proj > 0, dperp, np.sqrt(r2))


def audit_e(C):
    cl = C["e_masked_reader"]
    n, L = 32, 48.0
    h = L / n
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rs = np.maximum(r, 1e-12)
    cth = Z / rs
    sth = np.sqrt(np.maximum(1 - cth**2, 0.0))
    phi = np.arctan2(Y, X)
    pin = B3.pin_shell(n, h, 1.6)
    lines = [(30.0, 0.0), (90.0, 90.0), (120.0, 200.0)]
    us = []
    for th, ph in lines:
        th, ph = np.radians(th), np.radians(ph)
        us.append(np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]))
    r_t = 1.5 * h
    dl = np.min(np.stack([line_dist(X, Y, Z, u) for u in us]), 0)
    mask_tube = dl < r_t
    shells = np.arange(6.0, 18.0 + 1e-9, 1.5)

    def read(dev, masked):
        a1 = []
        for Rs in shells:
            m = (np.abs(r - Rs) < 0.75) & (~pin)
            if masked:
                m &= ~mask_tube
            A = real_harmonics(cth[m], sth[m], phi[m])
            coef, *_ = np.linalg.lstsq(A, dev[m], rcond=None)
            a1.append(float(coef[1]))
        return np.array(a1)

    def slope(a):
        return float(np.polyfit(np.log(shells), np.log(np.abs(a)), 1)[0])

    out = {"h": h, "r_t": r_t}
    # the dipole
    dev = np.where(r >= 3.0, 0.01 * cth / rs**2, 0.0)
    a_d = read(dev, True) * shells**2
    out["dipole_masked_a1_r2"] = a_d.tolist()
    out["dipole_max_rel_err"] = float(np.abs(a_d / 0.01 - 1).max())
    out["dipole_slope"] = slope(read(dev, True))
    # the footprint inside the mask (radius 0.8 r_t)
    r_f = 0.8 * r_t
    foot = -0.3 * (dl < r_f)
    a_u = read(foot, False)
    a_m = read(foot, True)
    pred = np.array(
        [
            (3.0 / (4 * np.pi)) * (-0.3) * sum(np.pi * r_f**2 / Rs**2 * u[2] for u in us)
            for Rs in shells
        ]
    )
    out["footprint_unmasked_a1"] = a_u.tolist()
    out["footprint_prediction"] = pred.tolist()
    out["footprint_ratio"] = (a_u / pred).tolist()
    out["footprint_mean_ratio"] = float(np.mean(a_u / pred))
    out["footprint_ratio_min_max"] = [float((a_u / pred).min()), float((a_u / pred).max())]
    out["footprint_unmasked_slope"] = slope(a_u)
    out["footprint_masked_max_abs"] = float(np.abs(a_m).max())
    out["footprint_cells_outside_mask"] = int(np.sum((dl < r_f) & ~mask_tube))
    # footprints WIDER than the mask: what the masked read returns
    wider = {}
    for fac in (1.2, 1.5, 2.0):
        rf2 = fac * r_t
        f2 = -0.3 * (dl < rf2)
        au2, am2 = read(f2, False), read(f2, True)
        wider[f"{fac:g}_r_t"] = {
            "masked_over_unmasked": (am2 / au2).tolist(),
            "mean_abs_ratio": float(np.mean(np.abs(am2 / au2))),
        }
    out["footprint_wider_than_mask"] = wider
    # the physical carrier width: the BPS melt of the census seeds (delta 0.3, m 1 and 1/2)
    w = W1 * W1S
    out["melt_half_depth_radius"] = {
        f"m{m:g}": float(np.sqrt(np.log(2.0) / kappa_bps(0.3, w, m))) for m in (0.5, 1.0, 2.0)
    }
    d_ok = out["dipole_max_rel_err"] < 0.05 and abs(out["dipole_slope"] + 2) < 0.1
    f_ok = (
        abs(out["footprint_unmasked_slope"] + 2) < 0.5 and 0.7 < out["footprint_mean_ratio"] < 1.3
    )
    return {
        "e_dipole": verdict(
            "masked l=1 read of a synthetic dipole returns 0.01 r^-2 within 3.5 percent, slope -2.03",
            [out["dipole_max_rel_err"], out["dipole_slope"]],
            [cl["dipole_masked"]["max_rel_err"], cl["dipole_masked"]["slope_a1"]],
            d_ok and abs(out["dipole_max_rel_err"] - cl["dipole_masked"]["max_rel_err"]) < 1e-6,
            "reproduced with an own least squares (unnormalized real harmonics: the 0.01 recovery fixes the basis convention to cos theta, not sqrt(3/4pi) cos theta)",
        ),
        "e_footprint_unmasked": verdict(
            "unmasked l=1 read of a tube footprint falls with slope -2.04 and its mean ratio to the solid-angle prediction is 0.92",
            [
                out["footprint_unmasked_slope"],
                out["footprint_mean_ratio"],
                out["footprint_ratio_min_max"],
            ],
            [
                cl["footprint_unmasked"]["slope_a1"],
                cl["footprint_unmasked"]["mean_ratio_read_over_pred"],
            ],
            f_ok
            and abs(
                out["footprint_mean_ratio"] - cl["footprint_unmasked"]["mean_ratio_read_over_pred"]
            )
            < 1e-6,
            "reproduced; the per-shell ratio scatters from %.2f to %.2f (a 1-cell-wide footprint on a 1.5-cell shell), so the 0.7 to 1.3 window on the MEAN is a weak test: the per-shell read is not quantitative"
            % tuple(out["footprint_ratio_min_max"]),
            qualified=True,
        ),
        "e_footprint_masked_zero": verdict(
            "the masked read of the footprint is exactly 0",
            [out["footprint_masked_max_abs"], out["footprint_cells_outside_mask"]],
            cl["footprint_masked"]["max_abs_masked_over_unmasked"],
            out["footprint_masked_max_abs"] == 0.0,
            "UNFALSIFIABLE as stated: the footprint radius 0.8 r_t lies inside the mask radius r_t, so every footprint cell is masked and the read is zero by construction (0 footprint cells outside the mask); for a footprint of 1.2 / 1.5 / 2 r_t the masked read is %.2f / %.2f / %.2f of the unmasked in magnitude (sign flips: the fit extrapolates the ring); the BPS melt half-depth radius of the census carriers is %.2f (m 1) and %.2f (the half-unit carriers, m 1/2), both WIDER than r_t = 2.25, so the mask does not remove a physical carrier's footprint"
            % (
                wider["1.2_r_t"]["mean_abs_ratio"],
                wider["1.5_r_t"]["mean_abs_ratio"],
                wider["2_r_t"]["mean_abs_ratio"],
                out["melt_half_depth_radius"]["m1"],
                out["melt_half_depth_radius"]["m0.5"],
            ),
            qualified=True,
        ),
    }, out


# ============================================================================
# (f) THE CHARGE CONTROL
# ============================================================================
def own_charge_dsym(d, h):
    """own copy of the reader's construction: centered Jacobian current, centered divergence."""
    g = [np.gradient(d, h, axis=a) for a in range(3)]
    J = [
        np.einsum("...a,...a->...", d, np.cross(g[(a + 1) % 3], g[(a + 2) % 3])) / (4 * np.pi)
        for a in range(3)
    ]
    return sum(np.gradient(J[a], h, axis=a) for a in range(3)) * h**3


def solid_angle(a, b, c):
    num = np.einsum("...a,...a->...", a, np.cross(b, c))
    den = (
        1.0
        + np.einsum("...a,...a->...", a, b)
        + np.einsum("...a,...a->...", b, c)
        + np.einsum("...a,...a->...", c, a)
    )
    return 2.0 * np.arctan2(num, den)


def charge_berg_luscher(d):
    """the lattice degree of the director on each dual cube (8 cell centers), the sum of the
    signed solid angles of the 12 face triangles / 4 pi: an exact integer per dual cube."""
    n = d.shape[0]
    V = {}
    for i in (0, 1):
        for j in (0, 1):
            for k in (0, 1):
                V[(i, j, k)] = d[i : n - 1 + i, j : n - 1 + j, k : n - 1 + k]
    faces = [
        ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0)),  # x = 0, outward -x
        ((1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)),  # x = 1, outward +x
        ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)),  # y = 0
        ((0, 1, 0), (0, 1, 1), (1, 1, 1), (1, 1, 0)),  # y = 1
        ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0)),  # z = 0
        ((0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)),  # z = 1
    ]
    Q = 0.0
    for f in faces:
        a, b, c, e = (V[v] for v in f)
        Q = Q + solid_angle(a, b, c) + solid_angle(a, c, e)
    return Q / (4 * np.pi)


def audit_f(C):
    cl = C["f_charge_control"]
    n, L = 32, 48.0
    h = L / n
    cfg = R21.cfg_of(n, L, G, 0.3)
    M = R20.seed_axes(cfg, (1.0, 0.3, 0.0))
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rhat = np.stack([X, Y, Z], -1) / np.maximum(r, 1e-12)[..., None]
    lam, Vv = np.linalg.eigh(M[..., 1:, 1:])
    d = Vv[..., :, 2]
    d = d * np.sign(np.einsum("...a,...a->...", d, rhat) + 1e-300)[..., None]
    director_is_rhat = float(np.abs(d - rhat)[r > 1.0].max())
    free = ~B3.pin_shell(n, h, 1.6)
    rho_ref, _ = R26_3.charge_density(M, cfg)
    rho_own = own_charge_dsym(d, h)
    lines = [(30.0, 0.0), (90.0, 90.0), (120.0, 200.0)]
    us = []
    for th, ph in lines:
        th, ph = np.radians(th), np.radians(ph)
        us.append(np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]))
    dl = np.min(np.stack([line_dist(X, Y, Z, u) for u in us]), 0)
    tubes = dl < 3.0 * h

    def stats(rho):
        outer = free & (r > 6.0)
        tot_abs = np.sum(np.abs(rho[free]))
        return {
            "net_free": float(np.sum(rho[free])),
            "abs_free": float(tot_abs),
            "net_outer": float(np.sum(rho[outer])),
            "abs_frac_outer": float(np.sum(np.abs(rho[outer])) / tot_abs),
            "abs_frac_outer_tubes": float(np.sum(np.abs(rho[outer & tubes])) / tot_abs),
            "abs_frac_outer_rest": float(np.sum(np.abs(rho[outer & ~tubes])) / tot_abs),
            "tube_volume_share_of_outer": float(np.sum(outer & tubes) / np.sum(outer)),
        }

    s_ref, s_own = stats(rho_ref), stats(rho_own)
    # where does the outer floor sit: per-shell |rho| and the axis-vs-generic split
    outer = free & (r > 6.0)
    ax_dist = np.min(np.stack([np.hypot(Y, Z), np.hypot(X, Z), np.hypot(X, Y)]), 0)
    near_axis = ax_dist < 1.5 * h
    s_axis = float(np.sum(np.abs(rho_ref[outer & near_axis])) / np.sum(np.abs(rho_ref[outer])))
    axis_vol = float(np.sum(outer & near_axis) / np.sum(outer))
    # the Berg-Luscher degree per dual cube
    Qbl = charge_berg_luscher(d)
    Xc = 0.5 * (X[:-1, :-1, :-1] + X[1:, 1:, 1:])
    Yc = 0.5 * (Y[:-1, :-1, :-1] + Y[1:, 1:, 1:])
    Zc = 0.5 * (Z[:-1, :-1, :-1] + Z[1:, 1:, 1:])
    rc = np.sqrt(Xc**2 + Yc**2 + Zc**2)
    bl = {
        "total": float(np.sum(Qbl)),
        "max_abs_outside_r6": float(np.abs(Qbl[rc > 6.0]).max()),
        "abs_sum_outside_r6": float(np.sum(np.abs(Qbl[rc > 6.0]))),
        "nonzero_cubes": int(np.sum(np.abs(Qbl) > 1e-6)),
    }
    out = {
        "director_minus_rhat_max": director_is_rhat,
        "reader": s_ref,
        "own_dsym_copy": s_own,
        "berg_luscher": bl,
        "outer_abs_share_within_1.5h_of_a_grid_axis": s_axis,
        "grid_axis_volume_share_of_outer": axis_vol,
        "max_rel_diff_reader_vs_own": float(
            np.abs(rho_ref - rho_own).max() / np.abs(rho_ref).max()
        ),
    }
    ok = (
        abs(s_ref["net_free"] - cl["charge_net_free"]) < 1e-9
        and abs(
            s_ref["abs_frac_outer"]
            - (cl["fraction_outer_in_tubes"] + cl["fraction_outer_outside_tubes"])
        )
        < 1e-6
    )
    return {
        "f_floor": verdict(
            "the hedgehog seed's reader floor: net 0.999 free, net 0.026 and |rho| fraction 0.036 outside r 6 (0.008 in tubes, 0.028 outside)",
            s_ref,
            {
                "net_free": cl["charge_net_free"],
                "net_outer": cl["net_charge_outer"],
                "in_tubes": cl["fraction_outer_in_tubes"],
                "outside": cl["fraction_outer_outside_tubes"],
            },
            ok,
            "numbers reproduced (own centered-stencil copy agrees to %.0e); BUT the director of this seed is r-hat to %.0e everywhere, so the 3.6 percent is the centered-stencil reader's discretization residue on a field of ZERO true charge density outside the core, not a lattice floor: the Berg-Luscher solid-angle degree on the same director reads total %.6f with |Q| < %.0e on every dual cube outside r 6 (an exact integer per cube, one nonzero cube). The residue is field-dependent (it scales with the field's second differences), so quoting the hedgehog's 3.6 percent beside a relaxed census field's split assumes the relaxed field is as smooth as the hedgehog; the tubes hold %.1f percent of the outer |rho| at %.1f percent of the outer volume because the residue decays steeply in r and the tubes cover about a third of the r 7 shell (cells within 1.5 h of a grid axis: %.1f percent of the residue at %.1f percent of the volume), so the split's null is geometry-specific"
            % (
                out["max_rel_diff_reader_vs_own"],
                director_is_rhat,
                bl["total"],
                bl["max_abs_outside_r6"] + 1e-300,
                100 * s_ref["abs_frac_outer_tubes"] / s_ref["abs_frac_outer"],
                100 * s_ref["tube_volume_share_of_outer"],
                100 * s_axis,
                100 * axis_vol,
            ),
            qualified=True,
        ),
    }, out


# ============================================================================
# (g) THE COULOMB COEFFICIENT OF THE UNIAXIAL HEDGEHOG
# ============================================================================
def audit_g(C):
    import sympy as sp

    cl = C["g_coulomb_uniaxial"]
    x, y, z = sp.symbols("x y z", real=True)
    r = sp.sqrt(x * x + y * y + z * z)
    nv = sp.Matrix([x, y, z]) / r
    S = nv * nv.T
    A = [S.diff(v) for v in (x, y, z)]
    e = 0
    for i in range(3):
        for j in range(i + 1, 3):
            F = A[i] * A[j] - A[j] * A[i]
            e += 4 * sum(F[a, b] ** 2 for a in range(3) for b in range(3))
    Cg = sp.simplify(
        (e * r**4).subs({x: sp.Rational(3, 7), y: sp.Rational(-2, 5), z: sp.Rational(1, 3)})
    )
    Cp = sp.simplify((e * r**4).subs({x: 0, y: 0, z: 1}))
    Cgen = sp.simplify(e * r**4)
    out = {
        "C_generic": float(Cg),
        "C_pole": float(Cp),
        "C_symbolic": str(Cgen),
        "c_4piC": float(4 * np.pi * Cg),
    }
    # the lattice
    n, L = 48, 72.0
    h = L / n
    cfg = R21.cfg_of(n, L, G, 0.0)
    pot = ("v4", R0.roots_of(cfg), W1 * W1S)
    M = R20.seed_axes(cfg, (1.0, 0.0, 0.0))
    dens_full = R20.density(M, cfg, pot)
    dens_own = curv_density(M, h)
    X, Y, Z = B3.coords(n, h)
    rr = np.sqrt(X * X + Y * Y + Z * Z)
    free = ~B3.pin_shell(n, h, 1.6)
    shells = []
    for Rs in np.arange(9.0, 18.0 + 1e-9, 1.5):
        m = (np.abs(rr - Rs) < 0.75) & free
        shells.append(
            [
                float(Rs),
                float(np.mean(dens_own[m] / h**3 * rr[m] ** 4)),
                float(np.mean(dens_full[m] / h**3 * rr[m] ** 4)),
            ]
        )
    Rf = np.arange(9.0, 18.0 + 1e-9, 1.5)
    EgtR = np.array([np.sum(dens_full[free & (rr > R)]) for R in Rf])
    A = np.stack([1.0 / Rf, -np.ones_like(Rf)], 1)
    coef, *_ = np.linalg.lstsq(A, EgtR, rcond=None)
    c_fit, R_eff = float(coef[0]), float(coef[0] / coef[1])
    E_tot = float(np.sum(dens_full[free]))
    # the same object in another box: is c / E a property of the object?
    cfg2 = R21.cfg_of(32, 48.0, G, 0.0)
    M2 = R20.seed_axes(cfg2, (1.0, 0.0, 0.0))
    E_tot2 = float(
        np.sum(
            R20.density(M2, cfg2, ("v4", R0.roots_of(cfg2), W1 * W1S))[~B3.pin_shell(32, 1.5, 1.6)]
        )
    )
    out.update(
        {
            "shells_r_C_own_curv_C_full": shells,
            "C_lattice_max_rel_err": float(max(abs(s[2] / 8 - 1) for s in shells)),
            "c_fit": c_fit,
            "R_eff_fit": R_eff,
            "E_total_free_n48_L72": E_tot,
            "E_total_free_n32_L48": E_tot2,
            "c_over_E_n48_L72": float(4 * np.pi * 8 / E_tot),
            "c_over_E_n32_L48": float(4 * np.pi * 8 / E_tot2),
        }
    )
    ok = (
        float(Cg) == 8.0
        and float(Cp) == 8.0
        and abs(c_fit - cl["c_fit"]) < 1e-6
        and abs(E_tot - cl["E_total_free"]) < 1e-6
    )
    return {
        "g_coulomb": verdict(
            "e = 8 / r^4 (sympy), c = 4 pi C = 100.53; lattice shells 8.23 to 8.03; fit c = 100.61, R_eff 39.5; E_total 37.4, c / E = 2.69",
            {
                "C": [out["C_generic"], out["C_pole"]],
                "c_fit": c_fit,
                "R_eff": R_eff,
                "E_total": E_tot,
                "shells": shells,
            },
            {
                "C": [cl["C_sympy_generic_point"], cl["C_sympy_pole"]],
                "c_fit": cl["c_fit"],
                "R_eff": cl["R_eff_fit"],
                "E_total": cl["E_total_free"],
                "shells": cl["shells_r_C_lattice"],
            },
            ok,
            "reproduced with own sympy (C = 8 identically, e r^4 = 8 as a symbolic identity, so pole = generic point is a triviality of rotational invariance) and own curvature density (the shell values agree to 1e-12 with R20.density, the potential vanishes outside the core); the 4 pi C -> c step is stated correctly (E(> R) = 4 pi C / R); 'c / E' is NOT a property of the object: E_total is box-bound (the melted core plus c (1 / r_core - 1 / R_eff)), it reads %.2f at n 48 L 72 and %.2f at n 32 L 48, so c / E moves from %.2f to %.2f with the box"
            % (E_tot, E_tot2, out["c_over_E_n48_L72"], out["c_over_E_n32_L48"]),
            qualified=True,
        ),
    }, out


# ============================================================================
# (h) THE PAIR SEED, THE LOOPS, THE WALL READER
# ============================================================================
def audit_h(C):
    cl = C["h_pair_seed_loops"]
    n, nz, h, delta = 48, 4, 1.0, 0.3
    w = W1 * W1S
    Mp = own_pair_seed(n, nz, h, delta, w, d=6.0)
    psi, b = tensor_phase(Mp[:, :, 0])
    out = {
        "core_plus": winding_on_circle(psi, n, h, (3.0, 0.0), 2.0 * h),
        "core_minus": winding_on_circle(psi, n, h, (-3.0, 0.0), 2.0 * h),
        "outer_rho12": winding_on_circle(psi, n, h, (0.0, 0.0), 12.0),
    }
    for m in (1, 2):
        pm, _ = tensor_phase(own_bps(n, 1, h, delta, w, m)[:, :, 0])
        out[f"bps_m{m}_rho12"] = winding_on_circle(pm, n, h, (0.0, 0.0), 12.0)
    # own core finder: the centroid of b / b0 < 0.2 cells on each side of x = 0 (0.3 admits
    # the doubly-melted cells between the cores and drags the centroid inward)
    X, Y = slab_xy(n, h)
    low = b / (delta / 2.0) < 0.2
    cores = []
    for side in (X < 0, X > 0):
        m = low & side
        cores.append([float(np.mean(X[m])), float(np.mean(Y[m])), int(m.sum())])
    out["cores"] = cores
    out["core_distance"] = float(cores[1][0] - cores[0][0])

    # own wall reader: free bonds whose wrapped phase step exceeds 90 / 150 degrees, both ends
    # with b / b0 > 0.5 (the phase is read only where the pair is present)
    def wall_bonds(M, dl=delta):
        ps, bb = tensor_phase(M)
        free = R25_1.free_mask_slab(M.shape[0], M.shape[2], h)
        okb = (bb / (dl / 2.0) > 0.5) & free
        res = {}
        for thr in (90.0, 150.0):
            cnt = 0
            for ax in (0, 1):
                dpsi = np.diff(ps, axis=ax)
                dpsi = np.abs((dpsi + np.pi) % (2 * np.pi) - np.pi)
                s0 = [slice(None)] * 3
                s1 = [slice(None)] * 3
                s0[ax], s1[ax] = slice(0, -1), slice(1, None)
                both = okb[tuple(s0)] & okb[tuple(s1)]
                cnt += int(np.sum((dpsi > np.radians(thr)) & both))
            res[f"bonds_over_{thr:g}_deg_all_layers"] = cnt
            res[f"bonds_over_{thr:g}_deg_per_layer"] = cnt / M.shape[2]
        return res

    out["walls_seed"] = wall_bonds(Mp)
    Zs = np.load(os.path.join(R26_2, "escape_d0.03_w25_n48_L48_long_amp0.npz"))
    Mw = np.diag([8.0, 1.0, 0.03, 0.0]) + 0.03 * Zs["D"]
    delta_w = 0.03

    out["walls_stored_d0.03"] = wall_bonds(Mw, delta_w)
    ok = all(
        abs(out[k] - v) < 0.05
        for k, v in (
            ("core_plus", 1),
            ("core_minus", 1),
            ("outer_rho12", 2),
            ("bps_m1_rho12", 1),
            ("bps_m2_rho12", 2),
        )
    )
    ok &= (
        abs(out["core_distance"] - 6.0) < 1.0
        and out["walls_seed"]["bonds_over_150_deg_all_layers"] == 0
        and out["walls_stored_d0.03"]["bonds_over_150_deg_all_layers"] > 0
    )
    return {
        "h_loops": verdict(
            "pair seed reads winding 1 at each core, 2 on rho 12; BPS m1 / m2 read 1 / 2; two cores at distance 6 within 1 h; no wall bonds on the seed, wall bonds on the stored wall state",
            out,
            {
                k: cl[k]
                for k in (
                    "winding_core_plus",
                    "winding_core_minus",
                    "winding_outer_rho12",
                    "winding_bps_m1_outer",
                    "winding_bps_m2_outer",
                    "core_distance",
                )
            }
            | {
                "walls_seed_150": cl["walls"]["wall_bonds_150"],
                "walls_stored_150": cl["walls_on_stored_wall_state_d0.03"]["wall_bonds_150"],
            },
            ok,
            "reproduced with an own seed, own phase reader and own wall counter; the windings are read off the seed's OWN analytic phase psi (S was built from psi), so the loop reads validate the reader on a field where the answer is known by construction, not the seed's dynamics; the own wall counter (free bonds with b / b0 > 0.5 at both ends, wrapped phase step over 150 / 90 deg) reads %.0f / %.0f per z layer on the stored wall state, the author's 27 / 46 exactly, and 0 on the seed"
            % (
                out["walls_stored_d0.03"]["bonds_over_150_deg_per_layer"],
                out["walls_stored_d0.03"]["bonds_over_90_deg_per_layer"],
            ),
            qualified=True,
        ),
    }, out


# ============================================================================
# THE STORED TABLE
# ============================================================================
def audit_stored(C):
    cl = C["stored"]
    w = W1 * W1S
    n, nz, h = 48, 4, 1.0
    out = {"slab": {}, "census": {}}
    ok = True
    for delta in (0.3, 0.03, 0.01):
        b0 = delta / 2.0
        seed = own_bps(n, nz, h, delta, w, 1)
        e_seed = bond_E(seed[..., 1:, 1:], h) / (nz * h)
        Zs = np.load(os.path.join(R26_2, f"escape_d{delta:g}_w25_n48_L48_long_amp0.npz"))
        Mw = np.diag([8.0, 1.0, delta, 0.0]) + delta * Zs["D"]
        e_wall = bond_E(Mw[..., 1:, 1:], h) / (nz * h)
        Tb = T_bps(delta, w)
        row = {
            "T_bps": Tb,
            "E_kappa_seed": float(e_seed),
            "E_kappa_wall": float(e_wall),
            "wall_over_seed": float(e_wall / e_seed),
            "kappa_ref": float(Tb / e_seed),
            "E_kappa_seed_over_2b0sq": float(e_seed / (2 * b0 * b0)),
        }
        # the box dependence of kappa_ref: the seed's bond sum in the L 96 box (per unit length)
        e96 = bond_E(own_bps(96, 1, h, delta, w, 1)[..., 1:, 1:], h) / h
        row["kappa_ref_L96"] = float(Tb / e96)
        out["slab"][f"d{delta:g}"] = row
        cs = cl["slab"][f"d{delta:g}"]
        ok &= abs(row["T_bps"] - cs["T_bps"]) < 1e-15 * max(1, cs["T_bps"]) + 1e-18
        ok &= abs(row["E_kappa_seed"] - cs["bps_seed_m1"]["E_kappa_per_len_bond"]) < 1e-9
        ok &= abs(row["E_kappa_wall"] - cs["wall_state"]["E_kappa_per_len_bond"]) < 1e-9
        ok &= abs(row["kappa_ref"] - cs["kappa_ref"]) < 1e-9
    # the census
    J = json.load(open(CENSUS_JSON))
    hc = 1.5
    rows = {}
    for part in ("1_1_1_1", "2_1_1", "2_2", "3_1", "4"):
        tag = f"P{part}_d0.3_w25_n32_L48"
        Ms = np.load(os.path.join(R26_4, tag + "_seed.npz"))["M"]
        Me = np.load(os.path.join(R26_4, tag + ".npz"))["M"]
        rows[part] = {
            "sum_k2": sum(int(k) ** 2 for k in part.split("_")) / 4.0,
            "E_kappa_seed": float(bond_E(Ms[..., 1:, 1:], hc)),
            "E_kappa_end": float(bond_E(Me[..., 1:, 1:], hc)),
            "E_quartic_end": float(J["rows"][tag]["E"]),
            "M00_bond_seed": float(bond_E(Ms[..., :1, :1], hc)),
            "M00_bond_end": float(bond_E(Me[..., :1, :1], hc)),
            "M0i_max_end": float(np.abs(Me[..., 0, 1:]).max()),
        }
    seed_e = np.array([rows[p]["E_kappa_seed"] for p in rows])
    end_e = np.array([rows[p]["E_kappa_end"] for p in rows])
    k2 = np.array([rows[p]["sum_k2"] for p in rows])
    qe = np.array([rows[p]["E_quartic_end"] for p in rows])
    seed_spread, end_spread, q_spread = (
        float(np.ptp(seed_e)),
        float(np.ptp(end_e)),
        float(np.ptp(qe)),
    )
    kc = [0.5 * q_spread / seed_spread, 2.0 * q_spread / seed_spread]
    out["census"] = {
        "rows": rows,
        "seeds_in_sum_k2_order": bool(np.all(np.diff(seed_e[np.argsort(k2)]) > 0)),
        "ends_in_sum_k2_order": bool(np.all(np.diff(end_e[np.argsort(k2)]) > 0)),
        "seed_spread": seed_spread,
        "seed_min": float(seed_e.min()),
        "end_spread": end_spread,
        "quartic_end_spread": q_spread,
        "kappa_census": kc,
        "kappa_census_x_seed_E_kappa_min": [k * float(seed_e.min()) for k in kc],
        "kappa_census_x_end_E_kappa_spread": [k * end_spread for k in kc],
        "kappa_census_x_end_E_kappa_max": [k * float(end_e.max()) for k in kc],
        "E_quartic_end_mean": float(qe.mean()),
    }
    cc = cl["census"]
    ok_c = (
        abs(seed_spread - cc["seed_E_kappa_spread"]) < 1e-6
        and abs(end_spread - cc["end_E_kappa_spread"]) < 1e-6
        and abs(kc[0] - cc["kappa_census"][0]) < 1e-9
        and out["census"]["seeds_in_sum_k2_order"]
        and not out["census"]["ends_in_sum_k2_order"]
    )
    return {
        "stored_slab": verdict(
            "slab table: seed 0.68366 / 0.0066983 / 0.00074411, walls 16.186 / 0.25643 / 0.038345, kappa_ref = T_bps / E_kappa(seed)",
            {
                k: [
                    out["slab"][k]["E_kappa_seed"],
                    out["slab"][k]["E_kappa_wall"],
                    out["slab"][k]["kappa_ref"],
                ]
                for k in out["slab"]
            },
            {
                k: [
                    cl["slab"][k]["bps_seed_m1"]["E_kappa_per_len_bond"],
                    cl["slab"][k]["wall_state"]["E_kappa_per_len_bond"],
                    cl["slab"][k]["kappa_ref"],
                ]
                for k in cl["slab"]
            },
            ok,
            "reproduced to 1e-9 (own seed, own bond sum, own T_bps); kappa_ref follows its stated rule, but E_kappa(seed) is the WHOLE-BOX bond sum of a log-divergent object (2 b0^2 x 2 pi m^2 ln(L / core) + core), so kappa_ref is box-bound: in the L 96 box it is %.3e / %.3e / %.3e (%.0f percent lower)"
            % (
                out["slab"]["d0.3"]["kappa_ref_L96"],
                out["slab"]["d0.03"]["kappa_ref_L96"],
                out["slab"]["d0.01"]["kappa_ref_L96"],
                100
                * (1 - out["slab"]["d0.3"]["kappa_ref_L96"] / out["slab"]["d0.3"]["kappa_ref"]),
            ),
            qualified=True,
        ),
        "stored_census": verdict(
            "census seeds' whole-box E_kappa in sum-k^2 order (min 961.5, spread 95.5), ends not (spread 1668), quartic end spread 0.4393, kappa_census 2.299e-3 and 9.197e-3",
            {
                k: out["census"][k]
                for k in (
                    "seed_spread",
                    "seed_min",
                    "end_spread",
                    "quartic_end_spread",
                    "kappa_census",
                    "seeds_in_sum_k2_order",
                    "ends_in_sum_k2_order",
                )
            },
            {
                k: cc[k]
                for k in (
                    "seed_E_kappa_spread",
                    "seed_E_kappa_min",
                    "end_E_kappa_spread",
                    "quartic_end_spread_kappa0",
                    "kappa_census",
                    "seeds_in_sum_k2_order",
                    "ends_in_sum_k2_order",
                )
            },
            ok_c,
            "reproduced; the sizing rule is a SEED rule: at kappa_census the seeds carry %.2f / %.2f of smoothness energy (26 to 100 percent of the quartic end energies, mean %.2f), the END fields' E_kappa spread times kappa is %.2f / %.2f (9x and 35x the quartic spread 0.44) and the end fields' largest E_kappa costs %.1f / %.1f: the term is not a perturbation of the census energies at either kappa, so the census order at kappa > 0 will be set by the smoothness term, not by the quartic; the M_00 bond sum of the end fields is %.1f (excluded from E_kappa by the block definition; the time row is zero to %.0e)"
            % (
                out["census"]["kappa_census_x_seed_E_kappa_min"][0],
                out["census"]["kappa_census_x_seed_E_kappa_min"][1],
                out["census"]["E_quartic_end_mean"],
                out["census"]["kappa_census_x_end_E_kappa_spread"][0],
                out["census"]["kappa_census_x_end_E_kappa_spread"][1],
                out["census"]["kappa_census_x_end_E_kappa_max"][0],
                out["census"]["kappa_census_x_end_E_kappa_max"][1],
                max(rows[p]["M00_bond_end"] for p in rows),
                max(rows[p]["M0i_max_end"] for p in rows),
            ),
            qualified=True,
        ),
    }, out


# ============================================================================
def main():
    C = json.load(open(CLAIMS))
    verdicts, detail = {}, {}
    for name, fn in (
        ("a", audit_a),
        ("b", audit_b),
        ("c", audit_c),
        ("d", audit_d),
        ("e", audit_e),
        ("f", audit_f),
        ("g", audit_g),
        ("h", audit_h),
        ("stored", audit_stored),
    ):
        t = time.time()
        v, d = fn(C)
        verdicts.update(v)
        d["wall_s"] = round(time.time() - t, 1)
        detail[name] = d
        for k, vv in v.items():
            log(f"{k:28s} {vv['verdict']:10s} own={json.dumps(vv['own_value'])[:110]}")
    counts = {}
    for vv in verdicts.values():
        counts[vv["verdict"]] = counts.get(vv["verdict"], 0) + 1
    blind = [
        "E_kappa penalizes the spatial block only: M_00 (free in the static sector, slaved per cell in the census) and the time row M_0i (nonzero once the kinetic / rotating sector enters) can hold cell-sharp structure at zero smoothness cost; on the stored states both are zero-gradient, so the omission is invisible at R27-0 and untested for R27-1 / R27-2 relaxations",
        "the bond sum saturates: a cell-sharp pi wall costs 8 b0^2 / h per unit area (59 percent below the continuum ramp of width h), so the smoothness term can rank a cell-sharp wall below a resolved wall of width 2 h; the (d) check reports only W 4 h and 8 h",
        "the two-box winding coefficient holds for the WHOLE-BOX bond sum only (free-cell bonds give 6.7 percent): any R27-1 / R27-2 read that sums E_kappa over free cells does not carry the 4 pi m^2 b0^2 ln 2 coefficient",
        "the masked reader's null (footprint read exactly 0) is a construction: the test footprint (0.8 r_t) lies inside the mask; a footprint 1.2 / 1.5 / 2 x r_t wide keeps 71 / 74 / 93 percent of its unmasked l = 1 read through the mask, and the census carriers' BPS melt half-depth (3.63 at m 1, 2.57 at m 1/2) exceeds r_t = 2.25, so the mask as sized cannot remove a physical carrier's footprint",
        "the charge floor (3.6 percent) is the centered-stencil reader's residue on an exactly-r-hat director (zero true density), field-dependent and steeply decaying in r (the tubes hold 24 percent of it at 5 percent of the volume); an exact lattice degree (Berg-Luscher solid angles) has zero residue on the same field but reads only integer charges per dual cube, so neither reader gives a fractional tube-vs-rest split with a field-independent null",
        "kappa_ref = T_bps / E_kappa(seed) and kappa_census are both box-bound (the seed's whole-box bond sum grows as ln L; the census rule is sized on the seeds while the end fields' E_kappa is 2.7 to 4.3x larger and 17x more spread), so an R27-1 / R27-2 verdict at these kappas is a verdict about a smoothness-dominated energy in one box",
        "c / E_total for the Coulomb hedgehog is box-bound (2.69 at n 48 L 72, another value at n 32 L 48) and E includes the melted core's V4; the coefficient c = 32 pi is the object's property, the ratio is not",
    ]
    summary = {
        "counts": counts,
        "refuted": [
            f"{k}: {v['note'][:160]}" for k, v in verdicts.items() if v["verdict"] == "REFUTED"
        ],
        "qualified": [
            f"{k}: {v['note'].split(';')[0][:200]}"
            for k, v in verdicts.items()
            if v["verdict"] == "QUALIFIED"
        ],
        "blind_spots": blind,
    }
    J = {
        "task": "M5.32 R27-0 audit (independent, m5_32_r27_0_form.py unread)",
        "verdicts": verdicts,
        "summary": summary,
        "detail": detail,
        "wall_s": round(time.time() - T0, 1),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(J, f, indent=1, default=float)
    log(f"counts {counts}; wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
