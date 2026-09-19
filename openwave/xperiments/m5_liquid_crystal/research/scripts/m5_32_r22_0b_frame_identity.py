"""M5.32 R22-0 (b): the Faber identity extended to a full frame field (the
biaxial vacuum), and the R21 polished fields read against it.

EQUATIONS FIRST
---------------
Spatial block M_sp = R Lambda R^T, Lambda = diag(l_1, l_2, l_3) constant,
R(x) in SO(3), body rates w_i = R^T d_i R (antisymmetric, vector W_i).
Then d_i M = R [w_i, Lambda] R^T, a symmetric off-diagonal matrix with
entries (w_i)_ab (l_b - l_a), and for two such matrices the commutator is
antisymmetric with vector u_i x u_j, u_i^a = g_a W_i^a, g_a = the gap
OPPOSITE to axis a (g_1 = l_2 - l_3, g_2 = l_1 - l_3, g_3 = l_1 - l_2, signs
immaterial). Hence, in the stack normalization of R22-0 (a):
    4 sum_{i<j} <F_ij, F_ij>_eta = 8 sum_{i<j} |u_i x u_j|^2
                                 = 8 sum_a G_a^2 sum_{i<j} ((W_i x W_j)^a)^2,
    G_a = prod_{b != a} (l_a - l_b)      (the two gaps ADJACENT to axis a)
and by Maurer-Cartan (W_i x W_j)^a = rho_ij[e_a] = e_a . (d_i e_a x d_j e_a),
the topological current of eigenvector a. So the certified static quartic of
a frame field is THREE dual-Maxwell energies, one per eigenvector, weighted
by G_a^2. Uniaxial (1, delta, delta): G = ((1 - delta)^2, 0, 0), the author's
identity. Biaxial (1, delta, 0) at delta 0.3: G^2 = (0.49, 0.0441, 0.09).
Sign-free lattice form of the current (no orientation needed, so it is
defined on non-orientable textures): with P_a = e_a e_a^T,
    rho_ij[e_a]^2 = (1/2) tr(F F^T),  F = [d_i P_a, d_j P_a].

Checks:
    A  exact jets: random R, W_i, Lambda; the stack density formula against
       the closed form, and the projector form against ((W_i x W_j)^a)^2.
       Negative control: the closed form with the OPPOSITE-gap-squared
       weights (g_a^2 in place of G_a) must fail (median residual over 0.1).
    B  the R21 polished fields (n32 L48, the three axes on the three
       boundaries at W1, plus the W1 x 25 electron): on the shells r 8 to 16
       the measured curvature density against the frame prediction built
       from the LOCAL eigenvalues and projectors, cells with a gap under
       GAP_MIN excluded and counted, and the polar cone abs(cos theta) > 0.8
       left out (the seeds' frame carries a line defect on the polar axis,
       where one cell turns the frame by order one and no stencil applies).
       Calibration rows: the pure seeds at h 1.5 and h 0.75, exact frame
       fields with constant eigenvalues, where the ratio must go to 1. Reads: the ratio measured / predicted
       per shell, the share of each eigenvector's current, the share of the
       curvature that sits in the 0i block (not covered by the identity).
       No PASS gate: it is a read of what carries the R21 tails.
       Audit notes (C6.2, finding 3): the prediction is built with the same
       fwd / bwd stencil as the measurement; with the central stencil the
       r 12 to 14 ratio is 1.101 (robust) and the r 8 to 10 ratio moves from
       1.24 to 1.37 (not robust). e r^4 is not a plateau on the W1 x 25
       electron (3.38 / 3.87 / 4.24 / 4.52 on the four shells) because the
       local eigenvalues are (0.05, 0.27, ~1), not (0, 0.3, 1): the local
       G_top^2 runs 0.36 to 0.51. The supportable wording is a Coulomb weight
       of 0.49 within about 15 percent, against 0.2401 for (1 - delta)^4.

Regenerate: python m5_32_r22_0b_frame_identity.py   (about 1 min)
"""

import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r22_0b_frame_identity.json")
R21_NPZ = os.path.join(DATA, "m5_32_r21_1")
GAP_MIN = 0.05
COS_MAX = 0.8  # shells are read off the polar axis, where the seeds carry their line defect
SHELLS = ((6.0, 8.0), (8.0, 10.0), (10.0, 12.0), (12.0, 14.0), (14.0, 16.0))


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


R20 = _load("m5_32_r20_1_axes", "m5_32_r20_1_axes.py")
B3 = R20.B3
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ================= A: exact jets =================
def hat(v):
    return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0.0]])


def rand_rot(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def check_a(samples=2000, seed=1):
    rng = np.random.default_rng(seed)
    worst, worst_proj, negs = 0.0, 0.0, []
    for _ in range(samples):
        R = rand_rot(rng)
        lam = rng.uniform(-1.0, 2.0, size=3)
        L3 = np.diag(lam)
        W = [rng.normal(size=3) for _ in range(3)]
        dM = []
        for Wi in W:
            D = np.zeros((4, 4))
            D[1:, 1:] = R @ (hat(Wi) @ L3 - L3 @ hat(Wi)) @ R.T
            dM.append(D)
        lhs = 0.0
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(dM[i], dM[j])
                lhs += 4.0 * float(B3.inner_eta(F, F))
        Gw = np.array([np.prod([lam[a] - lam[b] for b in range(3) if b != a]) for a in range(3)])
        gw = np.array([lam[1] - lam[2], lam[0] - lam[2], lam[0] - lam[1]])
        cur2 = np.zeros(3)
        for i in range(3):
            for j in range(i + 1, 3):
                cur2 += np.cross(W[i], W[j]) ** 2
        rhs = 8.0 * float(np.sum(Gw**2 * cur2))
        rhs_neg = 8.0 * float(np.sum(gw**4 * cur2))
        worst = max(worst, abs(lhs - rhs) / max(abs(rhs), 1e-300))
        negs.append(abs(lhs - rhs_neg) / max(abs(rhs_neg), 1e-300))
        # projector form of the current, axis by axis
        for a in range(3):
            Ea = np.zeros((3, 3))
            Ea[a, a] = 1.0
            dP = [R @ (hat(Wi) @ Ea - Ea @ hat(Wi)) @ R.T for Wi in W]
            c2 = 0.0
            for i in range(3):
                for j in range(i + 1, 3):
                    F = dP[i] @ dP[j] - dP[j] @ dP[i]
                    c2 += 0.5 * np.sum(F * F)
            worst_proj = max(worst_proj, abs(c2 - cur2[a]) / max(cur2[a], 1e-300))
    neg = float(np.median(negs))
    log(
        f"A: closed form max rel residual {worst:.2e}; projector current {worst_proj:.2e}; "
        f"negative control (opposite-gap weights) median residual {neg:.2e}"
    )
    return {
        "samples": samples,
        "max_rel_residual": worst,
        "projector_form_max_rel": worst_proj,
        "negative_control_median_rel_residual": float(neg),
        "PASS": bool(worst < 1e-9 and worst_proj < 1e-9 and neg > 0.1),
        "G2_biaxial_1_0.3_0": [0.49, 0.0441, 0.09],
        "G2_uniaxial_1_0.3_0.3": [0.2401, 0.0, 0.0],
    }


# ================= B: the R21 fields =================
def curv_density(M, cfg):
    """stack curvature density per cell (not h^3-weighted)."""
    e = np.zeros(M.shape[:3])
    for br, (A, wt) in B3.a_fields(M, cfg).items():
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                e += wt * 4.0 * B3.inner_eta(F, F)
    return e


def frame_prediction(M3, cfg):
    """8 sum_a G_a(x)^2 sum_{i<j} rho_ij[e_a]^2 with local eigenvalues, the stack's fwd/bwd branches."""
    lam, vec = np.linalg.eigh(M3)  # ascending: index 2 is the largest
    h = cfg["h"]
    Gw = np.stack(
        [
            np.prod(np.stack([lam[..., a] - lam[..., b] for b in range(3) if b != a], -1), -1)
            for a in range(3)
        ],
        -1,
    )
    gaps = np.minimum(lam[..., 1] - lam[..., 0], lam[..., 2] - lam[..., 1])
    per_axis = []
    for a in range(3):
        P = vec[..., :, a, None] * vec[..., None, :, a]
        c2 = np.zeros(M3.shape[:3])
        for br, wt in B3.branches("sym"):
            dP = [B3.d1(P, ax, h, br) for ax in range(3)]
            for i in range(3):
                for j in range(i + 1, 3):
                    F = dP[i] @ dP[j] - dP[j] @ dP[i]
                    c2 += wt * 0.5 * np.einsum("...ab,...ab->...", F, F)
        per_axis.append(8.0 * Gw[..., a] ** 2 * c2)
    return per_axis, lam, gaps


def read_field(tag, M=None):
    if M is None:
        path = os.path.join(R21_NPZ, tag + ".npz")
        if not os.path.exists(path):
            return None
        M = np.load(path)["M"]
    n = M.shape[0]
    cfg = B3.base_cfg(s=-1.0, g=float(abs(M[0, 0, 0, 0, 0])), n=n, L=48.0, delta=0.3)
    h = cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    e_full = curv_density(M, cfg)
    Ms = M.copy()
    Ms[..., 0, 1:] = 0.0
    Ms[..., 1:, 0] = 0.0
    Ms[..., 0, 0] = Ms[0, 0, 0, 0, 0]
    e_sp = curv_density(Ms, cfg)
    per_axis, lam, gaps = frame_prediction(M[..., 1:, 1:], cfg)
    e_pred = per_axis[0] + per_axis[1] + per_axis[2]
    rows = []
    for lo, hi in SHELLS:
        sh = (r >= lo) & (r < hi) & (np.abs(Z) < COS_MAX * r)
        ok = sh & (gaps > GAP_MIN)
        S = lambda a: float(a[ok].sum())  # noqa: E731
        pred = S(e_pred)
        rows.append(
            {
                "shell": [lo, hi],
                "cells": int(sh.sum()),
                "cells_gap_excluded": int((sh & ~ok).sum()),
                "E_full_r4": float((e_full[ok] * r[ok] ** 4).mean()) if ok.any() else None,
                "ratio_spatial_over_pred": S(e_sp) / pred if pred > 0 else None,
                "share_0i_block": 1.0 - S(e_sp) / S(e_full) if S(e_full) > 0 else None,
                "share_axis_small_mid_large": [
                    S(p) / pred if pred > 0 else None for p in per_axis
                ],
                "lam_mean_small_mid_large": (
                    [float(lam[..., a][ok].mean()) for a in range(3)] if ok.any() else None
                ),
            }
        )
    return {"tag": tag, "g": cfg["g"] if "g" in cfg else None, "rows": rows}


def check_b():
    seeds = []
    tags = [
        f"{o}_{b}_g8_d0.3_w1_n32" for o in ("S1", "Sd", "S0") for b in ("Bseed", "Bfar", "Bfree")
    ]
    tags += ["S1_Bseed_g8_d0.3_w25_n32", "S1_Bfree_g8_d0.3_w25_n32"]
    out = []
    # calibration: the pure seeds (exact frame fields with the vacuum eigenvalues outside the core),
    # at the R21 spacing and at half of it, so the reader's own lattice error is on the record
    for n_ in (32, 64):
        cfg = R20.cfg_of(n_, 48.0, 8.0)
        for o in ("S1", "Sd", "S0"):
            seeds.append((f"SEED_{o}_n{n_}", R20.seed_axes(cfg, R20.OBJECTS[o])))
    for t in [s_[0] for s_ in seeds] + tags:
        rd = read_field(t, dict(seeds).get(t))
        if rd is None:
            log(f"B {t}: missing")
            continue
        out.append(rd)
        for row in rd["rows"]:
            if row["shell"][0] in (8.0, 12.0):
                sh = row["share_axis_small_mid_large"]
                log(
                    f"B {t:28s} r {row['shell'][0]:4.0f}-{row['shell'][1]:<4.0f} spatial/pred "
                    f"{row['ratio_spatial_over_pred']:.3f}  0i share {row['share_0i_block']:+.3f}  axis shares "
                    f"{sh[0]:.2f}/{sh[1]:.2f}/{sh[2]:.2f}  excl {row['cells_gap_excluded']}/{row['cells']}  "
                    f"<e r^4> {row['E_full_r4']:.3f}"
                )
    return out


def main():
    res = {"a": check_a(), "b": check_b(), "gap_min": GAP_MIN}
    res["runtime_s"] = time.time() - T0
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1)
    log(f"PASS A {res['a']['PASS']} -> {os.path.basename(OUT_JSON)}")


if __name__ == "__main__":
    main()
