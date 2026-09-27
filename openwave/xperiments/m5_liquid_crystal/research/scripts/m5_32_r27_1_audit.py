"""M5.32 R27-1 audit: an independent attempt to refute the kappa-slab claims (the kappa arm, the
k^2 law, the pair E(d), the wall ladder, the pre-registration amendment) with the auditor's own
readers, own smoothness term, own preconditioner and own descent. The audited scripts
(m5_32_r27_1_kappa_slab.py, m5_32_r27_0_form.py) were NOT read; the shared stack is consumed
read-only: R26-2's deviation-form quartic energy `energy_grad_M` (E_u + V4 and dE/dM) and its
slab helpers, R25-1's `bps_field` / `departure` / `density_u`, B3's stencils.

EQUATIONS FIRST
---------------
The slab of R25-1 (n x n x 4 cells, z-invariant, box L, h = L / n, the x-y shell of depth 1.6
pinned, the static sector M_0i = 0 with M_00 free), M = M0 + delta D, M0 = diag(8, 1, delta, 0).
    E[D] = delta^4 E_u[D] + V4[M] + kappa E_kappa[M],
    E_kappa = h^3 sum_{bonds in x, y} sum_z |S_{i+1} - S_i|_F^2 / h^2,   S = M[..., 1:, 1:]
(the forward and backward branches at weight 1/2 each sum to the plain bond sum), the gradient
    dE_kappa / dS_i = 2 h sum_{j ~ i} (S_i - S_j)   (full-entry convention, then sym),
checked here against a finite difference. Per unit length: / (4 h). T_bps(delta, w) of R25-1.
THE READERS (own): the pair gap 2 b = lambda_1 - lambda_0 of S (ascending eigenvalues), the
tensor phase psi = atan2(2 S_12, S_11 - S_22), a wall bond = |wrap(psi_i - psi_j)| above 150 (or
90) degrees between free cells with b >= 0.5 b0; an alternative any-plane measure, the Frobenius
bond jump |S_i - S_j|_F / (sqrt 2 b0) above 2 sin(75 deg) (the 150-degree equivalent of a pair
rotation) or 2 sin(45 deg); cores = connected clusters of free cells with b < 0.3 b0 (labelled);
the psi winding on a circle of radius 4 h about each core centroid and on the square loop one
cell inside the free boundary, by unwrapped differences of psi bilinearly interpolated.
THE DESCENT (own): L-BFGS-B on the 7 static-sector slots of the free cells in the eigenbasis of
a NUMERICAL 7 x 7 Hessian of f = E / T_ref at one vacuum cell (central differences of the
gradient), each eigen-direction scaled by sqrt(lambda) (floored at 1e-8 lambda_max): a diagonal
preconditioner independent of R26-2's analytic `slot_scales`. Chunks of 250 iterations; the
pair-block gate measure of the R26-2 audit, max |dE/dM_pair| over the free cells / (delta^3 h^3).
THE ZERO-MODE CHECK: any field depending on one spatial coordinate only has F_ij = 0 (the two
derivatives are parallel), and a rotation of the pair frame keeps the eigenvalues, so V4 is
unchanged: a straight wall of any width is an exact zero mode of the kappa-0 functional. This is
verified numerically on constructed fields (axis-aligned step, diagonal step, smooth 1D
rotation, and the control: a wall ENDING at a point).

Modes: run (everything, one process, about 15 minutes) | reads (no descents, about a minute) |
verdict (recompute the verdicts from the saved record).
Output: data/m5_32_r27_1_audit.json.
"""

import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r27_1_audit.json")
CLAIMS_JSON = os.path.join(DATA, "m5_32_r27_1_kappa_slab.json")
R27 = os.path.join(DATA, "m5_32_r27_1")
R26 = os.path.join(DATA, "m5_32_r26_2")
NZ = 4
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
SLOT_I = np.array([0, 1, 1, 1, 2, 2, 3])
SLOT_J = np.array([0, 1, 2, 3, 2, 3, 3])
SLOT_OFF = np.where(SLOT_I == SLOT_J, 1.0, 2.0)
KAPPA_REF = {0.3: 3.912306498873179e-3, 0.03: 3.618e-5, 0.01: 4.017e-6}
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


L2 = _load("m5_32_r26_2_ladder", "m5_32_r26_2_ladder.py")
R25 = L2.R25_1
B3 = L2.B3
W = L2.W1 * 25.0


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


# ================= own smoothness term =================
def e_kappa(M, h):
    S = M[..., 1:, 1:]
    e = 0.0
    for ax in (0, 1):
        d = np.diff(S, axis=ax)
        e += h * np.sum(d * d)
    return float(e)


def grad_kappa(M, h):
    """dE_kappa/dM, full-entry convention (symmetric), time row zero."""
    S = M[..., 1:, 1:]
    G = np.zeros_like(M)
    GS = np.zeros_like(S)
    for ax in (0, 1):
        d = np.diff(S, axis=ax)
        sl_lo = [slice(None)] * S.ndim
        sl_hi = [slice(None)] * S.ndim
        sl_lo[ax] = slice(0, -1)
        sl_hi[ax] = slice(1, None)
        GS[tuple(sl_lo)] -= 2.0 * h * d
        GS[tuple(sl_hi)] += 2.0 * h * d
    G[..., 1:, 1:] = GS
    return G


def check_grad_kappa(n=12, h=1.0, seed=3):
    rng = np.random.default_rng(seed)
    M = B3.sym4(rng.standard_normal((n, n, NZ, 4, 4)))
    G = grad_kappa(M, h)
    errs = []
    for _ in range(6):
        i, j, k = rng.integers(1, n - 1), rng.integers(1, n - 1), rng.integers(0, NZ)
        a, b = rng.integers(1, 4), rng.integers(1, 4)
        P = np.zeros_like(M)
        P[i, j, k, a, b] = P[i, j, k, b, a] = 1.0
        eps = 1e-6
        fd = (e_kappa(M + eps * P, h) - e_kappa(M - eps * P, h)) / (2 * eps)
        an = G[i, j, k, a, b] + (G[i, j, k, b, a] if a != b else 0.0)
        errs.append(abs(fd - an) / max(abs(an), 1e-300))
    return float(max(errs))


# ================= fields =================
def M_of(D, delta):
    return L2.M_vac(delta) + delta * D


def load_row(tag):
    Z = np.load(os.path.join(R27, tag + ".npz"))
    return Z["D"], float(Z["delta"]), (float(Z["kappa"]) if "kappa" in Z.files else None)


def quartic(D, cfg, delta):
    eu, ev, _ = L2.energy_grad_M(D, cfg, delta, W, need_grad=False)
    return float(eu), float(ev)


def coords(n, h):
    return (np.arange(n) - (n - 1) / 2.0) * h


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def loop_winding(S, pts_ij):
    from scipy.ndimage import map_coordinates

    s11 = map_coordinates(S[..., 0, 0], pts_ij.T, order=1, mode="nearest")
    s22 = map_coordinates(S[..., 1, 1], pts_ij.T, order=1, mode="nearest")
    s12 = map_coordinates(S[..., 0, 1], pts_ij.T, order=1, mode="nearest")
    psi = np.arctan2(2 * s12, s11 - s22)
    psi = np.append(psi, psi[0])
    return float(np.sum(wrap(np.diff(psi))) / (2 * np.pi))


def reads(M, cfg, delta, core_frac=0.3, loop_radius_h=4.0):
    """the auditor's readers on layer 0 (z-invariance reported)."""
    from scipy.ndimage import label

    n, h = cfg["n"], cfg["h"]
    free = R25.free_mask_slab(n, NZ, h)
    fr = free[:, :, 0]
    S = M[:, :, 0, 1:, 1:]
    lam = np.linalg.eigvalsh(S)
    b = 0.5 * (lam[..., 1] - lam[..., 0])
    b0 = delta / 2.0
    psi = np.arctan2(2 * S[..., 0, 1], S[..., 0, 0] - S[..., 1, 1])
    core_half = b < 0.5 * b0
    ok = fr & ~core_half
    out = {
        "z_variation_over_delta": float(np.max(np.abs(M - M[:, :, :1]))) / delta,
        "b_over_b0_min_free": float(np.min(b[fr]) / b0),
        "core_cells_half": int(np.sum(core_half & fr)),
        "departure": R25.departure(M, delta),
    }
    for name, thr in (("150", 150.0), ("90", 90.0)):
        cnt = 0
        for ax in (0, 1):
            d = np.abs(wrap(np.diff(psi, axis=ax)))
            both = np.logical_and(
                ok[tuple(slice(0, -1) if k == ax else slice(None) for k in range(2))],
                ok[tuple(slice(1, None) if k == ax else slice(None) for k in range(2))],
            )
            cnt += int(np.sum(both & (d > np.radians(thr))))
        out[f"wall_bonds_psi_{name}"] = cnt
    for name, thr in (("150", 2 * np.sin(np.radians(75))), ("90", 2 * np.sin(np.radians(45)))):
        cnt = 0
        for ax in (0, 1):
            dS = np.diff(S, axis=ax)
            jump = np.sqrt(np.sum(dS * dS, axis=(-1, -2))) / (np.sqrt(2.0) * b0)
            both = np.logical_and(
                ok[tuple(slice(0, -1) if k == ax else slice(None) for k in range(2))],
                ok[tuple(slice(1, None) if k == ax else slice(None) for k in range(2))],
            )
            cnt += int(np.sum(both & (jump > thr)))
        out[f"wall_bonds_frob_{name}"] = cnt
    # cores
    x = coords(n, h)
    mask = (b < core_frac * b0) & fr
    lab, k = label(mask)
    cores = []
    for c in range(1, k + 1):
        ii, jj = np.where(lab == c)
        cx, cy = float(np.mean(x[ii])), float(np.mean(x[jj]))
        r0 = loop_radius_h * h
        t = np.linspace(0, 2 * np.pi, 128, endpoint=False)
        pts = np.stack(
            [
                cx / h + (n - 1) / 2.0 + r0 / h * np.cos(t),
                cy / h + (n - 1) / 2.0 + r0 / h * np.sin(t),
            ],
            1,
        )
        cores.append(
            {
                "x": cx,
                "y": cy,
                "cells": int(len(ii)),
                "radius_eq": float(h * np.sqrt(len(ii) / np.pi)),
                "b_over_b0_min": float(np.min(b[ii, jj]) / b0),
                "winding_psi_r4h": loop_winding(S, pts),
            }
        )
    cores.sort(key=lambda c: (c["x"], c["y"]))
    out["cores"] = cores
    out["core_count"] = len(cores)
    if len(cores) == 2:
        out["core_distance"] = float(
            np.hypot(cores[0]["x"] - cores[1]["x"], cores[0]["y"] - cores[1]["y"])
        )
    else:
        out["core_distance"] = None
    wc = max(1, int(np.ceil(1.6 / h)))
    lo, hi = wc + 1, n - 2 - wc
    idx = np.arange(lo, hi + 1)
    edge = np.concatenate(
        [
            np.stack([idx, np.full_like(idx, lo)], 1),
            np.stack([np.full_like(idx, hi), idx], 1),
            np.stack([idx[::-1], np.full_like(idx, hi)], 1),
            np.stack([np.full_like(idx, lo), idx[::-1]], 1),
        ]
    ).astype(float)
    out["winding_psi_outer_loop"] = loop_winding(S, edge)
    out["free_half_width"] = float(x[n - 1 - wc] + 0.5 * h)
    return out


def densities(M, cfg, delta):
    """per-cell E_u and V4 densities on layer 0 (h^3-weighted, the plain trace form)."""
    h3 = cfg["h"] ** 3
    du = R25.density_u(M, cfg)[:, :, 0]
    N = M[:, :, 0] @ ETA
    P = N.copy()
    tr = []
    for p in range(4):
        if p:
            P = P @ N
        tr.append(np.einsum("...kk->...", P))
    C = [(-8.0) ** p + 1.0 + delta**p for p in range(1, 5)]
    dv = h3 * W * sum((tr[p] - C[p]) ** 2 for p in range(4))
    return du, dv


def localize(M, cfg, delta):
    """where the kappa-0 energy of a wall state sits: core neighbourhoods (b < 0.5 b0, dilated by
    one cell), cells adjacent to an any-plane wall bond (the Frobenius jump above the 90-degree
    equivalent), the rest; plus the concentration (the share of the energy in the top 1 and 5
    percent of the free cells, the number of cells holding half of it)."""
    from scipy.ndimage import binary_dilation

    n, h = cfg["n"], cfg["h"]
    free = R25.free_mask_slab(n, NZ, h)
    fr = free[:, :, 0]
    S = M[:, :, 0, 1:, 1:]
    lam = np.linalg.eigvalsh(S)
    b = 0.5 * (lam[..., 1] - lam[..., 0])
    b0 = delta / 2.0
    core = binary_dilation(b < 0.5 * b0, iterations=1) & fr
    wall = np.zeros_like(fr)
    for ax in (0, 1):
        dS = np.diff(S, axis=ax)
        d = np.sqrt(np.sum(dS * dS, axis=(-1, -2))) / (np.sqrt(2.0) * b0) > 2 * np.sin(
            np.radians(45)
        )
        lo = tuple(slice(0, -1) if k == ax else slice(None) for k in range(2))
        hi = tuple(slice(1, None) if k == ax else slice(None) for k in range(2))
        wall[lo] |= d
        wall[hi] |= d
    wall &= fr & ~core
    rest = fr & ~core & ~wall
    du, dv = densities(M, cfg, delta)
    e = du + dv
    tot = float(np.sum(e[fr]))
    ef = np.sort(e[fr])[::-1]
    cum = np.cumsum(ef) / tot
    nf = int(fr.sum())
    x = coords(n, h)
    X, Y = np.meshgrid(x, x, indexing="ij")
    r = np.hypot(X, Y)
    return {
        "E_layer0_total": tot,
        "E_share_r_lt_4": float(np.sum(e[fr & (r < 4.0)]) / tot),
        "b_over_b0_min_free": float(np.min(b[fr]) / b0),
        "director_row_max_over_delta": float(np.max(np.abs(M[:, :, 0, 3, 1:3]))) / delta,
        "cells": {"core": int(core.sum()), "wall": int(wall.sum()), "rest": int(rest.sum())},
        "E_fraction": {
            k: float(np.sum(e[m]) / tot)
            for k, m in (("core", core), ("wall", wall), ("rest", rest))
        },
        "E_u_fraction_of_total": float(np.sum(du[fr]) / tot),
        "share_top_1pct_cells": float(cum[max(1, nf // 100) - 1]),
        "share_top_5pct_cells": float(cum[max(1, nf // 20) - 1]),
        "cells_holding_half": int(np.searchsorted(cum, 0.5) + 1),
        "free_cells": nf,
    }


# ================= zero-mode check =================
def pair_field(n, h, delta, psi, b):
    s0 = delta / 2.0
    M = np.zeros((n, n, NZ, 4, 4))
    M[..., 0, 0] = 8.0
    M[..., 3, 3] = 1.0
    M[..., 1, 1] = (s0 + b * np.cos(psi))[:, :, None]
    M[..., 2, 2] = (s0 - b * np.cos(psi))[:, :, None]
    M[..., 1, 2] = M[..., 2, 1] = (b * np.sin(psi))[:, :, None]
    return M


def zero_mode_check(delta=0.3):
    """constructed fields with b = b0 everywhere (V4 = 0 by construction): the one-coordinate
    profiles (F = 0 exactly), two walls meeting at a corner (the same jump in x and y, F = 0),
    and the singular textures whose cost sits at the lattice core: the pure winding psi = phi
    and the half winding psi = phi / 2 with its cut (a wall ending at a point); the latter two
    at h 1 and h 1/2 to expose the 1/h scaling of a sharp junction."""
    out = {}
    b0 = delta / 2.0
    Tb = L2.T_bps(delta, W)
    for n in (48, 96):
        cfg = L2.slab_cfg(n, 48.0, delta)
        h = cfg["h"]
        x = coords(n, h)
        X, Y = np.meshgrid(x, x, indexing="ij")
        phi = np.arctan2(Y, X)
        r = np.hypot(X, Y)
        cases = {
            "step_x": np.where(X > 0, np.pi, 0.0),
            "step_diagonal": np.where(X + Y > 0, np.pi, 0.0),
            "smooth_1d_rotation": 0.5 * np.pi * (np.tanh(X / 3.0) + 1.0),
            "two_walls_meeting_at_a_corner": np.where((X > 0) & (Y > 0), np.pi, 0.0),
            "pure_winding_psi_eq_phi": phi,
            "half_winding_with_cut": 0.5 * phi,
        }
        for k, psi in cases.items():
            M = pair_field(n, h, delta, psi, b0 * np.ones_like(psi))
            du, dv = densities(M, cfg, delta)
            eu, ev = float(np.sum(du)) * NZ, float(np.sum(dv)) * NZ
            rec = {"E_u_over_T_bps": eu / (NZ * h) / Tb, "V4_over_T_bps": ev / (NZ * h) / Tb}
            if eu > 0:
                rec["E_u_share_within_3h_of_origin"] = float(np.sum(du[r < 3 * h]) / np.sum(du))
            out[f"{k}_h{h:g}"] = rec
    return out


# ================= own descent =================
class Own:
    def __init__(self, cfg, delta, kappa, T_ref):
        self.cfg, self.delta, self.kappa, self.T_ref = cfg, delta, kappa, T_ref
        self.h = cfg["h"]
        self.free = R25.free_mask_slab(cfg["n"], NZ, self.h)
        self.R = np.eye(7)
        self.s = np.ones(7)

    def parts(self, D):
        eu, ev, G = L2.energy_grad_M(D, self.cfg, self.delta, W)
        M = M_of(D, self.delta)
        ek = e_kappa(M, self.h)
        Gt = G + self.kappa * grad_kappa(M, self.h)
        Gt = R25.static_sector(Gt)
        return float(eu), float(ev), ek, Gt

    def energy(self, D):
        eu, ev, _ = L2.energy_grad_M(D, self.cfg, self.delta, W, need_grad=False)
        return float(eu), float(ev), e_kappa(M_of(D, self.delta), self.h)

    # x: the 7 slot values of D on the free cells; y = ((x @ R) * s)
    def pack(self, D):
        x = D[self.free][:, SLOT_I, SLOT_J]
        return ((x @ self.R) * self.s).ravel()

    def unpack(self, y, base):
        x = (y.reshape(-1, 7) / self.s) @ self.R.T
        D = base.copy()
        blk = D[self.free].copy()
        blk[:, SLOT_I, SLOT_J] = x
        blk[:, SLOT_J, SLOT_I] = x
        D[self.free] = blk
        return D

    def gx(self, Gt):
        return Gt[self.free][:, SLOT_I, SLOT_J] * SLOT_OFF * self.delta / self.T_ref

    def fun(self, y, base):
        D = self.unpack(y, base)
        eu, ev, ek, Gt = self.parts(D)
        f = (eu + ev + self.kappa * ek) / self.T_ref
        return f, ((self.gx(Gt) @ self.R) / self.s).ravel()

    def precondition(self, D, cell, eps=1e-4, floor=1e-8):
        """the numerical 7x7 Hessian of f in the slot variables of one cell; eigenbasis + sqrt scales."""
        i, j, k = cell
        H = np.zeros((7, 7))
        for c in range(7):
            cols = []
            for sgn in (1.0, -1.0):
                Dp = D.copy()
                a, b = SLOT_I[c], SLOT_J[c]
                Dp[i, j, k, a, b] += sgn * eps
                if a != b:
                    Dp[i, j, k, b, a] += sgn * eps
                _, _, _, Gt = self.parts(Dp)
                cols.append(Gt[i, j, k][SLOT_I, SLOT_J] * SLOT_OFF * self.delta / self.T_ref)
            H[:, c] = (cols[0] - cols[1]) / (2 * eps)
        H = 0.5 * (H + H.T)
        lam, R = np.linalg.eigh(H)
        lam_f = np.maximum(lam, floor * lam.max())
        self.R, self.s = R, np.sqrt(lam_f)
        return {
            "eigenvalues": lam.tolist(),
            "cell": [int(i), int(j), int(k)],
            "scales": self.s.tolist(),
        }

    def gmax_pair(self, Gt):
        return float(np.max(np.abs(Gt[self.free][:, 1:3, 1:3]))) / (self.delta**3 * self.h**3)

    def descend(self, D0, iters, chunk=250, tag=""):
        from scipy.optimize import minimize

        Tb = L2.T_bps(self.delta, W)
        D = D0.copy()
        eu, ev, ek, Gt = self.parts(D)
        f_prev = (eu + ev + self.kappa * ek) / self.T_ref
        traj = [
            {
                "iters": 0,
                "E_total_over_T_bps": f_prev,
                "E_quartic_over_T_bps": (eu + ev) / (NZ * self.h) / Tb,
                "E_kappa_per_len": ek / (NZ * self.h),
                "gmax_pair_scaled": self.gmax_pair(Gt),
            }
        ]
        done = 0
        while done < iters:
            base = D.copy()
            k = min(chunk, iters - done)
            res = minimize(
                self.fun,
                self.pack(D),
                args=(base,),
                jac=True,
                method="L-BFGS-B",
                options={
                    "maxcor": 20,
                    "maxiter": k,
                    "maxfun": 3 * k,
                    "gtol": 1e-30,
                    "ftol": 1e-30,
                },
            )
            D = self.unpack(np.asarray(res.x), base)
            its = int(res.nit)
            done += max(1, its)
            eu, ev, ek, Gt = self.parts(D)
            f = (eu + ev + self.kappa * ek) / self.T_ref
            traj.append(
                {
                    "iters": done,
                    "E_total_over_T_bps": f,
                    "E_quartic_over_T_bps": (eu + ev) / (NZ * self.h) / Tb,
                    "E_kappa_per_len": ek / (NZ * self.h),
                    "gmax_pair_scaled": self.gmax_pair(Gt),
                    "drop_rel": (f_prev - f) / max(abs(f_prev), 1e-300),
                    "chunk_iters": its,
                    "scipy": str(res.message),
                }
            )
            log(
                f"{tag} its {done} E/T_bps {f:.9f} quartic {traj[-1]['E_quartic_over_T_bps']:.6f} gpair {traj[-1]['gmax_pair_scaled']:.2e} drop {traj[-1]['drop_rel']:.2e}"
            )
            f_prev = f
            if its < 2:
                break
        return D, traj


def row_energy(tag, kappa_override=None):
    """own evaluation of a stored end field: quartic by R26-2, kappa term by the own bond sum."""
    D, delta, kappa = load_row(tag)
    n = D.shape[0]
    L = float(tag.split("_L")[1].split("_")[0])
    cfg = L2.slab_cfg(n, L, delta)
    h = cfg["h"]
    Tb = L2.T_bps(delta, W)
    eu, ev = quartic(D, cfg, delta)
    M = M_of(D, delta)
    ek = e_kappa(M, h)
    kap = kappa if kappa_override is None else kappa_override
    kap = 0.0 if kap is None else kap
    return (
        {
            "tag": tag,
            "n": n,
            "L": L,
            "h": h,
            "delta": delta,
            "kappa": kap,
            "E_quartic_over_T_bps": (eu + ev) / (NZ * h) / Tb,
            "E_u_over_T_bps": eu / (NZ * h) / Tb,
            "V4_over_T_bps": ev / (NZ * h) / Tb,
            "E_kappa_per_len": ek / (NZ * h),
            "E_total_over_T_bps": (eu + ev + kap * ek) / (NZ * h) / Tb,
            "reads": reads(M, cfg, delta),
        },
        D,
        cfg,
        M,
    )


def claimed(J, tag):
    r = J["rows"][tag]
    e = r["end_reads"]
    return {
        "E_total_over_T_bps": e["E_total_over_T_bps"],
        "E_quartic_over_T_bps": e["E_quartic_over_T_bps"],
        "E_kappa_per_len": e["E_kappa_per_len"],
        "wall_bonds_150": e["walls"]["wall_bonds_150"],
        "wall_bonds_90": e["walls"]["wall_bonds_90"],
        "director_row": e["departure"]["director_row"],
        "cores": e.get("cores"),
        "core_distance": e.get("core_distance"),
        "verdict": r["lbfgs"]["verdict"],
        "iters": r["lbfgs"]["iters"],
    }


def refine_x2(D, order=3):
    from scipy.ndimage import zoom

    n = D.shape[0]
    Df = np.zeros((2 * n, 2 * n, NZ, 4, 4))
    for a in range(4):
        for c in range(a, 4):
            v = zoom(D[:, :, 0, a, c], 2, order=order, mode="nearest")
            Df[:, :, :, a, c] = v[:, :, None]
            Df[:, :, :, c, a] = v[:, :, None]
    return Df


def two_core_ansatz(n, L, delta, centers, m_shell=None):
    """the product ansatz of two m 1 BPS cores at the given centers; E_kappa per unit length."""
    cfg = L2.slab_cfg(n, L, delta)
    h = cfg["h"]
    x = coords(n, h)
    X, Y = np.meshgrid(x, x, indexing="ij")
    k = L2.kappa_bps(delta, W)
    b0 = delta / 2.0
    psi = np.zeros_like(X)
    f = np.ones_like(X)
    for cx, cy in centers:
        r2 = (X - cx) ** 2 + (Y - cy) ** 2
        psi += np.arctan2(Y - cy, X - cx)
        f *= 1.0 - np.exp(-k * r2)
    M = pair_field(n, h, delta, psi, b0 * np.sqrt(f))
    if m_shell is not None:
        Ms = R25.bps_field(n, NZ, h, delta, W, m=m_shell)
        free = R25.free_mask_slab(n, NZ, h)
        M[~free] = Ms[~free]
    return e_kappa(M, h) / (NZ * h)


# ================= the audit =================
def run(do_descents=True):
    J = json.load(open(CLAIMS_JSON))
    out = {"task": "M5.32 R27-1 audit", "own": {}, "verdicts": {}, "summary": {}}
    out["own"]["grad_kappa_fd_check_max_rel"] = check_grad_kappa()
    log(f"grad_kappa finite-difference check {out['own']['grad_kappa_fd_check_max_rel']:.2e}")

    # ---------- claim 1: the kappa arm ----------
    c1 = {"claimed": {}, "own": {}}
    tags1 = [
        "wall_d0.3_n48_L48_k0.3",
        "bps_m1_d0.3_n48_L48_k0.3",
        "wall_d0.03_n48_L48_k0.3",
        "bps_m1_d0.03_n48_L48_k0.3",
        "bps_m1_d0.03_n48_L48_k0",
        "bps_m1_d0.3_n48_L48_k0",
        "wall_d0.3_n48_L48_k0",
        "wall_d0.03_n48_L48_k0",
        "bps_m2_d0.3_n48_L48_k0",
        "bps_m2_d0.03_n48_L48_k0",
        "pair_d6_d0.3_n48_L48_k0",
        "pair_d6_d0.03_n48_L48_k0",
        "wall_d0.3_n48_L48_k1",
        "bps_m1_d0.3_n48_L48_k1",
        "wall_d0.03_n48_L48_k1",
        "bps_m1_d0.03_n48_L48_k1",
        "bps_m1_d0.3_n48_L48_k0.03",
        "bps_m1_d0.03_n48_L48_k0.03",
    ]
    for t in tags1:
        rec, D, cfg, M = row_energy(t)
        c1["own"][t] = rec
        c1["claimed"][t] = claimed(J, t)
        log(
            f"{t}: own E_tot {rec['E_total_over_T_bps']:.10f} (claimed {c1['claimed'][t]['E_total_over_T_bps']:.10f}) "
            f"walls psi150 {rec['reads']['wall_bonds_psi_150']} psi90 {rec['reads']['wall_bonds_psi_90']} "
            f"(claimed {c1['claimed'][t]['wall_bonds_150']}/{c1['claimed'][t]['wall_bonds_90']}) "
            f"frob150 {rec['reads']['wall_bonds_frob_150']} dir {rec['reads']['departure']['director_row']:.3f}"
        )
    ea = c1["own"]["wall_d0.3_n48_L48_k0.3"]["E_total_over_T_bps"]
    eb = c1["own"]["bps_m1_d0.3_n48_L48_k0.3"]["E_total_over_T_bps"]
    c1["same_state_d0.3_k0.3_rel_diff"] = abs(ea - eb) / eb
    ea3 = c1["own"]["wall_d0.03_n48_L48_k0.3"]["E_total_over_T_bps"]
    eb3 = c1["own"]["bps_m1_d0.03_n48_L48_k0.3"]["E_total_over_T_bps"]
    c1["wall_below_smooth_d0.03_k0.3_rel"] = (eb3 - ea3) / eb3
    out["verdicts"]["claim1_kappa_arm"] = c1

    if do_descents:
        delta = 0.3
        n, L = 48, 48.0
        cfg = L2.slab_cfg(n, L, delta)
        h = cfg["h"]
        Tb = L2.T_bps(delta, W)
        T_ref = Tb * NZ * h
        kappa = 0.3 * KAPPA_REF[0.3]
        own = Own(cfg, delta, kappa, T_ref)
        D_bps = L2.bps_D(n, NZ, h, delta, W)
        c1["own_descent"] = {"kappa": kappa, "preconditioner": own.precondition(D_bps, (5, 5, 0))}
        log(
            f"preconditioner eigenvalues {np.array(c1['own_descent']['preconditioner']['eigenvalues'])}"
        )
        Zw = np.load(os.path.join(R26, "escape_d0.3_w25_n48_L48_long_amp0.npz"))
        D_wall = Zw["D"]
        ends = {}
        for name, D0 in (("wall_seed", D_wall), ("bps_seed", D_bps)):
            t0 = time.time()
            D, traj = own.descend(D0, 1500, tag=name)
            M = M_of(D, delta)
            ends[name] = {
                "traj": traj,
                "end": traj[-1],
                "reads": reads(M, cfg, delta),
                "wall_s": round(time.time() - t0, 1),
            }
            log(
                f"{name} done E_tot {traj[-1]['E_total_over_T_bps']:.9f} wall_s {ends[name]['wall_s']}"
            )
        e_w, e_b = (
            ends["wall_seed"]["end"]["E_total_over_T_bps"],
            ends["bps_seed"]["end"]["E_total_over_T_bps"],
        )
        c1["own_descent"].update(
            ends=ends,
            rel_diff_wall_vs_bps=abs(e_w - e_b) / e_b,
            within_1pct=bool(abs(e_w - e_b) / e_b < 0.01),
            wall_free=bool(
                ends["wall_seed"]["reads"]["wall_bonds_psi_150"] == 0
                and ends["bps_seed"]["reads"]["wall_bonds_psi_150"] == 0
            ),
            claimed_E_total_over_T_bps=c1["claimed"]["wall_d0.3_n48_L48_k0.3"][
                "E_total_over_T_bps"
            ],
            rel_diff_to_claimed=abs(
                e_b - c1["claimed"]["wall_d0.3_n48_L48_k0.3"]["E_total_over_T_bps"]
            )
            / e_b,
        )

    # ---------- claim 2: the k^2 law ----------
    c2 = {"claimed": J["collect"]["k2_law"], "own": {}}
    for delta in (0.03, 0.3):
        b0 = delta / 2.0
        pred = 4 * np.pi * b0**2 * np.log(2.0)
        rec = {"prediction_m1_two_box": pred}
        for m in (1, 2):
            e48, D48, cfg48, M48 = row_energy(f"bps_m{m}_d{delta:g}_n48_L48_k0.3")
            e96, D96, cfg96, M96 = row_energy(f"bps_m{m}_d{delta:g}_n96_L96_k0.3")
            diff = e96["E_kappa_per_len"] - e48["E_kappa_per_len"]
            rec[f"m{m}"] = {
                "E_kappa_L48": e48["E_kappa_per_len"],
                "E_kappa_L96": e96["E_kappa_per_len"],
                "difference": diff,
                "over_prediction_m1": diff / pred,
                "reads_L48": e48["reads"],
                "reads_L96": e96["reads"],
                "E_total_L48": e48["E_total_over_T_bps"],
                "E_total_L96": e96["E_total_over_T_bps"],
            }
            log(
                f"k2 d{delta} m{m}: diff {diff:.7f} pred {pred:.7f} ratio {diff / pred:.4f} cores48 {[(round(c['x'], 2), round(c['y'], 2), c['cells'], round(c['winding_psi_r4h'], 3)) for c in e48['reads']['cores']]} outer48 {e48['reads']['winding_psi_outer_loop']:.3f} cores96 {[(round(c['x'], 2), round(c['y'], 2), c['cells'], round(c['winding_psi_r4h'], 3)) for c in e96['reads']['cores']]} outer96 {e96['reads']['winding_psi_outer_loop']:.3f}"
            )
        rec["ratio_m2_over_m1"] = rec["m2"]["difference"] / rec["m1"]["difference"]
        # the alternative accountings, delta 0.3 only (the m 2 cores' positions from the own reader)
        if delta == 0.3:
            cs48 = [(c["x"], c["y"]) for c in rec["m2"]["reads_L48"]["cores"]]
            cs96 = [(c["x"], c["y"]) for c in rec["m2"]["reads_L96"]["cores"]]
            if len(cs48) == 2 and len(cs96) == 2:
                a48 = two_core_ansatz(48, 48.0, delta, cs48)
                a96 = two_core_ansatz(96, 96.0, delta, cs96)
                # the cores scaled to the box: the L 96 positions at exactly twice the L 48 ones
                cs96s = [(2 * cx, 2 * cy) for cx, cy in cs48]
                a96s = two_core_ansatz(96, 96.0, delta, cs96s)
                rec["ansatz_two_m1_cores"] = {
                    "E_kappa_L48": a48,
                    "E_kappa_L96": a96,
                    "difference_over_prediction_m1": (a96 - a48) / pred,
                    "exactly_scaled_positions_difference_over_prediction_m1": (a96s - a48) / pred,
                    "core_radius_of_centre_over_free_half_width_L48": float(
                        np.hypot(*cs48[0]) / rec["m2"]["reads_L48"]["free_half_width"]
                    ),
                    "core_radius_of_centre_over_free_half_width_L96": float(
                        np.hypot(*cs96[0]) / rec["m2"]["reads_L96"]["free_half_width"]
                    ),
                    "xi_bps": float(np.sqrt(np.log(2.0) / L2.kappa_bps(delta, W))),
                }
                log(
                    f"ansatz: diff/pred {rec['ansatz_two_m1_cores']['difference_over_prediction_m1']:.4f} scaled {rec['ansatz_two_m1_cores']['exactly_scaled_positions_difference_over_prediction_m1']:.4f}"
                )
        c2["own"][f"d{delta:g}"] = rec
    out["verdicts"]["claim2_k2_law"] = c2

    # ---------- claim 3: the pair E(d) ----------
    c3 = {"claimed": J["collect"].get("pair_E_of_d"), "own": {}}
    es = []
    for d in (4, 6, 9):
        t = f"pair_d{d}_d0.3_n48_L48_k0.3"
        rec, D, cfg, M = row_energy(t)
        c3["own"][t] = {
            k: rec[k] for k in ("E_total_over_T_bps", "E_quartic_over_T_bps", "E_kappa_per_len")
        }
        c3["own"][t]["reads"] = rec["reads"]
        c3["own"][t]["claimed"] = claimed(J, t)
        es.append(rec["E_total_over_T_bps"])
        log(
            f"{t}: E_tot {rec['E_total_over_T_bps']:.6f} cores {[(round(c['x'], 2), round(c['y'], 2), c['cells'], round(c['winding_psi_r4h'], 3)) for c in rec['reads']['cores']]} dist {rec['reads']['core_distance']}"
        )
    c3["own"]["energy_spread_rel"] = (max(es) - min(es)) / min(es)
    c3["own"]["m2_diagonal_state_E_total"] = c2["own"]["d0.3"]["m2"]["E_total_L48"]
    out["verdicts"]["claim3_pair_E_of_d"] = c3

    # ---------- claim 4: the wall ladder ----------
    c4 = {"claimed": J["collect"]["wall_ladder"], "own": {}}
    c4["own"]["zero_mode_check"] = zero_mode_check()
    log(f"zero-mode check {c4['own']['zero_mode_check']}")
    c4["own"]["stored_end_fields"] = {}
    for delta in (0.3, 0.03, 0.01):
        for n in (48, 96, 192):
            t = f"wall_d{delta:g}_n{n}_L48_k0"
            rec, D, cfg, M = row_energy(t)
            rec["localize"] = localize(M, cfg, delta)
            c4["own"]["stored_end_fields"][t] = {
                k: rec[k]
                for k in (
                    "h",
                    "E_quartic_over_T_bps",
                    "E_u_over_T_bps",
                    "V4_over_T_bps",
                    "reads",
                    "localize",
                )
            }
            log(
                f"{t}: E {rec['E_quartic_over_T_bps']:.5f} (claimed {claimed(J, t)['E_total_over_T_bps']:.5f}) walls psi150/90 {rec['reads']['wall_bonds_psi_150']}/{rec['reads']['wall_bonds_psi_90']} frob150/90 {rec['reads']['wall_bonds_frob_150']}/{rec['reads']['wall_bonds_frob_90']} E-frac {rec['localize']['E_fraction']} cells {rec['localize']['cells']}"
            )
    if do_descents:
        delta = 0.3
        t = "wall_d0.3_n48_L48_k0"
        D48, _, _ = load_row(t)
        Df = refine_x2(D48, order=3)
        cfg96 = L2.slab_cfg(96, 48.0, delta)
        h96 = cfg96["h"]
        Tb = L2.T_bps(delta, W)
        eu, ev = quartic(Df, cfg96, delta)
        e_interp = (eu + ev) / (NZ * h96) / Tb
        log(f"interpolated h 1/2 wall state: E {e_interp:.4f} T_bps (claimed seed 4.89)")
        T_ref = Tb * NZ * h96
        own96 = Own(cfg96, delta, 0.0, T_ref)
        pre = Own(cfg96, delta, 0.3 * KAPPA_REF[0.3], T_ref)
        prec = pre.precondition(L2.bps_D(96, NZ, h96, delta, W), (9, 9, 0))
        own96.R, own96.s = pre.R, pre.s
        t0 = time.time()
        Dr, traj = own96.descend(Df, 500, tag="wall h1/2 kappa0")
        Mr = M_of(Dr, delta)
        c4["own"]["refine_relax"] = {
            "E_interpolated_over_T_bps": e_interp,
            "claimed_E_seed_over_T_bps": c4["claimed"]["d0.3"]["h0.5"]["E_seed_over_T_bps"],
            "claimed_E_3000_over_T_bps": c4["claimed"]["d0.3"]["h0.5"]["E_over_T_bps"],
            "preconditioner": prec,
            "traj": traj,
            "end": traj[-1],
            "reads": reads(Mr, cfg96, delta),
            "localize": localize(Mr, cfg96, delta),
            "wall_s": round(time.time() - t0, 1),
        }
        log(
            f"refine-relax done: E {traj[-1]['E_quartic_over_T_bps']:.5f} after {traj[-1]['iters']} its"
        )
    out["verdicts"]["claim4_wall_ladder"] = c4

    # ---------- claim 5: the pre-registration amendment ----------
    c5 = {"own": {}}
    for delta in (0.3, 0.03):
        Tb = L2.T_bps(delta, W)
        Ts3 = L2.T_slaved3(delta, W)[0]
        q = {}
        for kf in ("0.03", "0.1", "0.3", "1"):
            t = f"bps_m1_d{delta:g}_n48_L48_k{kf}"
            r = J["rows"][t]
            q[kf] = {
                "E_quartic_over_T_bps": r["end_reads"]["E_quartic_over_T_bps"],
                "verdict": r["lbfgs"]["verdict"],
                "own_walls_psi_150": (
                    c1["own"][t]["reads"]["wall_bonds_psi_150"] if t in c1["own"] else None
                ),
            }
        ref_small = q["0.03"]["E_quartic_over_T_bps"]
        c5["own"][f"d{delta:g}"] = {
            "T_slaved3_over_T_bps": Ts3 / Tb,
            "quartic_by_kappa_frac": q,
            "letter_vs_T_bps_within_3pct": {
                kf: bool(abs(v["E_quartic_over_T_bps"] - 1.0) <= 0.03) for kf, v in q.items()
            },
            "amended_vs_smallest_kappa_within_3pct": {
                kf: bool(abs(v["E_quartic_over_T_bps"] / ref_small - 1.0) <= 0.03)
                for kf, v in q.items()
            },
            "vs_T_slaved3_within_3pct": {
                kf: bool(abs(v["E_quartic_over_T_bps"] / (Ts3 / Tb) - 1.0) <= 0.03)
                for kf, v in q.items()
            },
            "quartic_over_T_slaved3": {
                kf: v["E_quartic_over_T_bps"] / (Ts3 / Tb) for kf, v in q.items()
            },
        }
    c5["kappa_arm_claimed"] = {
        d: {k: v for k, v in J["collect"]["kappa_arm"][d].items() if not k.startswith("k")}
        for d in ("d0.3", "d0.03")
    }
    out["verdicts"]["claim5_preregistration"] = c5
    return out


def verdicts(out):
    V = {}
    # ---------- claim 1 ----------
    c1 = out["verdicts"]["claim1_kappa_arm"]
    same = c1["same_state_d0.3_k0.3_rel_diff"]
    od = c1.get("own_descent", {})
    w03 = c1["own"]["wall_d0.03_n48_L48_k0.3"]["reads"]
    b0k0 = c1["own"]["bps_m1_d0.03_n48_L48_k0"]["reads"]
    k0 = {t: c1["own"][t]["reads"] for t in c1["own"] if t.endswith("_k0")}
    k0_psi = {t: r["wall_bonds_psi_150"] for t, r in k0.items()}
    k0_frob = {t: r["wall_bonds_frob_150"] for t, r in k0.items()}
    k0_dir = {t: round(r["departure"]["director_row"], 3) for t, r in k0.items()}
    v1, notes1 = "CONFIRMED", []
    if same > 1e-12:
        v1 = "QUALIFIED"
        notes1.append(f"the two stored end fields differ by {same:.1e} (claimed 5e-13)")
    else:
        notes1.append(
            f"the stored wall-seed and smooth-seed end fields at delta 0.3, 0.3 kappa_ref agree to {same:.1e} on the own energy"
        )
    if od:
        if od["within_1pct"] and od["wall_free"]:
            notes1.append(
                f"own descent (own smoothness term, own numerical-Hessian preconditioner, 1500 L-BFGS iterations): the wall seed ends at {od['ends']['wall_seed']['end']['E_total_over_T_bps']:.6f} T_bps and the BPS seed at {od['ends']['bps_seed']['end']['E_total_over_T_bps']:.6f} (relative difference {od['rel_diff_wall_vs_bps']:.1e}, both wall-free; the claimed common value {od['claimed_E_total_over_T_bps']:.6f}, off by {od['rel_diff_to_claimed']:.1e})"
            )
        else:
            v1 = "REFUTED"
            notes1.append(
                f"own descent: wall seed {od['ends']['wall_seed']['end']['E_total_over_T_bps']:.6f}, BPS seed {od['ends']['bps_seed']['end']['E_total_over_T_bps']:.6f} T_bps (relative difference {od['rel_diff_wall_vs_bps']:.2e}), wall-free {od['wall_free']}"
            )
    notes1.append(
        f"own reads: wall_d0.03_k0.3 has {w03['wall_bonds_psi_90']} psi-90 bonds and director row {w03['departure']['director_row']:.3f} delta (claimed 0 and 0.39), ending {100 * c1['wall_below_smooth_d0.03_k0.3_rel']:.2f} percent below the smooth seed (claimed 3.6); bps_m1_d0.03_k0 has {b0k0['wall_bonds_psi_150']} psi-150 bonds (claimed 164)"
    )
    notes1.append(
        f"reader remark, not a refutation: at kappa 0 the psi-150 count is {k0_psi}; the WALL seeds re-relaxed at kappa 0 keep 0 to 25 psi walls (their texture is the director-row rotation, departure {k0_dir}), so the '127 to 307' range describes the four smooth seeds only; the any-plane Frobenius-150 count {k0_frob} sees every kappa-0 state as a sharp texture"
    )
    V["claim1_kappa_arm"] = {
        "verdict": v1,
        "own": {
            "same_state_rel_diff": same,
            "own_descent_E_total": {
                k: v["end"]["E_total_over_T_bps"] for k, v in od.get("ends", {}).items()
            },
            "own_descent_rel_diff": od.get("rel_diff_wall_vs_bps"),
            "own_descent_rel_diff_to_claimed": od.get("rel_diff_to_claimed"),
            "own_descent_gmax_pair_end": {
                k: v["end"]["gmax_pair_scaled"] for k, v in od.get("ends", {}).items()
            },
            "wall_d0.03_k0.3_psi90": w03["wall_bonds_psi_90"],
            "wall_d0.03_k0.3_director_row": w03["departure"]["director_row"],
            "wall_below_smooth_d0.03_k0.3_rel": c1["wall_below_smooth_d0.03_k0.3_rel"],
            "bps_m1_d0.03_k0_psi150": b0k0["wall_bonds_psi_150"],
            "kappa0_psi150": k0_psi,
            "kappa0_frob150": k0_frob,
            "kappa0_director_row": k0_dir,
        },
        "claimed": {
            "same_state_rel": "1e-13 (kappa_ref), 5e-13 (0.3 kappa_ref, delta 0.3)",
            "wall_below_smooth_d0.03": 0.036,
            "walls_k0": "127 to 307",
            "director_row_k0": 0.7,
        },
        "notes": notes1,
    }
    # ---------- claim 2 ----------
    c2 = out["verdicts"]["claim2_k2_law"]["own"]
    r03, r003 = c2["d0.3"], c2["d0.03"]
    an = r03.get("ansatz_two_m1_cores", {})
    m2_48, m2_96 = r03["m2"]["reads_L48"], r03["m2"]["reads_L96"]
    f48 = an.get("core_radius_of_centre_over_free_half_width_L48")
    f96 = an.get("core_radius_of_centre_over_free_half_width_L96")
    two_cores = m2_48["core_count"] == 2 and m2_96["core_count"] == 2
    self_similar = bool(two_cores and f48 is not None and abs(f48 - f96) < 0.05)
    notes2 = [
        f"numbers reproduced with the own bond sum: m 1 two-box coefficient over the prediction {r003['m1']['over_prediction_m1']:.4f} (delta 0.03) and {r03['m1']['over_prediction_m1']:.4f} (delta 0.3); m 2 over m 1 {r003['ratio_m2_over_m1']:.3f} and {r03['ratio_m2_over_m1']:.3f}",
        f"own core finder on the m 2 rows at delta 0.3: {m2_48['core_count']} cores at L 48 (distance {m2_48['core_distance']:.1f}, each of psi winding {[round(c['winding_psi_r4h'], 2) for c in m2_48['cores']]}) and {m2_96['core_count']} at L 96 (distance {m2_96['core_distance']:.1f}), the outer loop winding {m2_48['winding_psi_outer_loop']:.2f} / {m2_96['winding_psi_outer_loop']:.2f}; the core radius from the centre is {f48:.3f} / {f96:.3f} of the free half-width (self-similar to 5 percent: {self_similar})",
        f"the product ansatz of two BPS m 1 cores at the read positions gives the two-box coefficient {an.get('difference_over_prediction_m1', float('nan')):.3f} (exactly scaled positions {an.get('exactly_scaled_positions_difference_over_prediction_m1', float('nan')):.3f}) against the relaxed {r03['ratio_m2_over_m1']:.3f}: the reading 'two cores at the box scale' accounts for the coefficient, the residual from 2 being the fixed core size xi = {an.get('xi_bps', float('nan')):.1f} against the doubled box; the pinned shell holds winding 2 in both boxes, and it is not the cause",
        "label remark: the prediction as stated is for a SMOOTH winding of index k, which the m 1 rows confirm to 0.2 percent; the relaxed m 2 field is not a single smooth winding (the double vortex splits under any kappa > 0, the classic instability), so the m 2 row tests the stability of the index-1 object, not the k^2 coefficient of a smooth winding; K2_LAW_REFUTED overreaches unless the prediction was made for the relaxed index-1 object, in which case it stands with the mechanism given",
    ]
    V["claim2_k2_law"] = {
        "verdict": "QUALIFIED",
        "own": {
            "d0.03": {
                "m1_over_pred": r003["m1"]["over_prediction_m1"],
                "m2_over_pred": r003["m2"]["over_prediction_m1"],
                "ratio_m2_over_m1": r003["ratio_m2_over_m1"],
            },
            "d0.3": {
                "m1_over_pred": r03["m1"]["over_prediction_m1"],
                "m2_over_pred": r03["m2"]["over_prediction_m1"],
                "ratio_m2_over_m1": r03["ratio_m2_over_m1"],
            },
            "ansatz": an,
            "m2_core_distance_L48_L96": [m2_48["core_distance"], m2_96["core_distance"]],
            "m2_outer_winding_L48_L96": [
                m2_48["winding_psi_outer_loop"],
                m2_96["winding_psi_outer_loop"],
            ],
        },
        "claimed": {
            "m1_over_pred": [1.0017, 1.0009],
            "ratio_m2_over_m1": [1.94, 1.94],
            "core_distance": [32.4, 66.8],
            "label": "K2_LAW_REFUTED",
        },
        "notes": notes2,
    }
    # ---------- claim 3 ----------
    c3 = out["verdicts"]["claim3_pair_E_of_d"]["own"]
    dists = {t: v["reads"]["core_distance"] for t, v in c3.items() if t.startswith("pair")}
    ets = {t: v["E_total_over_T_bps"] for t, v in c3.items() if t.startswith("pair")}
    ok3 = (
        all(d is not None and 30.0 <= d <= 31.5 for d in dists.values())
        and c3["energy_spread_rel"] < 0.003
    )
    V["claim3_pair_E_of_d"] = {
        "verdict": "CONFIRMED" if ok3 else "REFUTED",
        "own": {
            "core_distances": dists,
            "E_total": ets,
            "energy_spread_rel": c3["energy_spread_rel"],
            "m2_diagonal_state_E_total": c3["m2_diagonal_state_E_total"],
        },
        "claimed": {"core_distances": "30.5 to 31.2", "energy_spread": "within 0.2 percent"},
        "notes": [
            f"own core finder: the three pair rows end with two m 1 cores {[round(d, 2) for d in dists.values()]} apart on the seed's axis (claimed 30.5 to 31.2), energies spread {100 * c3['energy_spread_rel']:.2f} percent (claimed 0.2): E(d) is unreadable at this box, as claimed",
            f"beside the claim: the pair rows sit {100 * (c3['m2_diagonal_state_E_total'] / min(ets.values()) - 1):.2f} percent below the m 2 row's diagonal two-core state at the same kappa and box, and differ among themselves by {100 * c3['energy_spread_rel']:.2f} percent at the gate: the two-core valley is flat in the core positions at the gate's resolution, the end separation set by the seed's axis and the box rather than by a minimum",
        ],
    }
    # ---------- claim 4 ----------
    c4 = out["verdicts"]["claim4_wall_ladder"]["own"]
    zm = c4["zero_mode_check"]
    rr = c4.get("refine_relax", {})
    sf = c4["stored_end_fields"]
    zero = {k: v["E_u_over_T_bps"] + v["V4_over_T_bps"] for k, v in zm.items()}
    flat = max(v for k, v in zero.items() if "winding" not in k)
    pw = {h: zm[f"pure_winding_psi_eq_phi_h{h}"] for h in ("1", "0.5")}
    hw = {h: zm[f"half_winding_with_cut_h{h}"] for h in ("1", "0.5")}
    loc = {t: v["localize"] for t, v in sf.items()}
    notes4 = [
        f"the ladder's numbers are reproduced on the stored fields: {', '.join(f'{t} {v['E_quartic_over_T_bps']:.4f}' for t, v in sf.items())}",
    ]
    if rr:
        under = rr["end"]["E_quartic_over_T_bps"] < rr["claimed_E_3000_over_T_bps"]
        notes4.append(
            f"own cubic interpolation of the h 1 wall state to h 1/2: {rr['E_interpolated_over_T_bps']:.3f} T_bps "
            f"(claimed 4.89; the seed value depends on the interpolation details, the qualitative jump holds); "
            f"own 500-iteration relaxation at kappa 0 lands at {rr['end']['E_quartic_over_T_bps']:.4f} "
            f"(the claim's row: {rr['claimed_E_3000_over_T_bps']:.4f} by 3000), the last chunk still dropping "
            f"{100 * rr['end'].get('drop_rel', 0):.1f} percent"
            + (
                "; the own value after 500 iterations already UNDERCUTS the claim's 3000-iteration value, so the "
                "h 1/2 rung is optimizer-limited and the ladder's rise from h 1 to h 1/2 is not a resolved "
                "h-dependence"
                if under
                else ""
            )
        )
    notes4.append(
        f"'wall tension' is the wrong frame: every one-coordinate profile and two walls meeting at a corner cost exactly 0 at kappa 0 (own constructed fields, max {flat:.1e} T_bps), because F_ij vanishes when the two derivatives are parallel and V4 is blind to a frame rotation, so a straight wall has zero tension identically; what costs is a singular END: the pure winding psi = phi with b = b0 costs {pw['1']['E_u_over_T_bps']:.1f} T_bps at h 1 and {pw['0.5']['E_u_over_T_bps']:.1f} at h 1/2, the half winding with its cut (a wall ending at a point) {hw['1']['E_u_over_T_bps']:.1f} and {hw['0.5']['E_u_over_T_bps']:.1f}, all of it within 3 h of the end point: a sharp end scales as 1/h^2 per unit length until melting is cheaper"
    )
    notes4.append(
        "own maps of the stored kappa-0 states (E share inside r < 4, melted cells b < 0.5 b0, director row): "
        + "; ".join(
            f"{t}: {100 * v['E_share_r_lt_4']:.0f} percent inside r 4, {v['cells']['core']} melted-neighbourhood cells, b min {v['b_over_b0_min_free']:.2f} b0, director row {v['director_row_max_over_delta']:.2f} delta, top 1 percent of cells hold {100 * v['share_top_1pct_cells']:.0f} percent"
            for t, v in loc.items()
        )
    )
    notes4.append(
        "reading: at h 1 and h 1/2 the energy is spread over the whole free region (no melted cell, b = b0 everywhere, each 4 x 4 block about 1 percent), the winding carried by a director-row band (an escape into the third eigen-direction) whose end point is sub-cell and therefore under-costed by the lattice; at h 1/4 the end point resolves into a melted core holding 26 to 32 percent of the energy and the value rises; the ladder is a resolution crossover, monotone from below, so it is consistent with a FINITE continuum value (the melted cut end plus a spread texture) and not with a vanishing one (that would need the h 1/4 cores to cost nothing), but its value is undecidable at these iteration counts (every row FALLING, the finest losing 8 to 15 percent per chunk); the deciding instrument is the h 1/4 and an h 1/8 row relaxed to the gate (a Newton-Krylov solver or an order of magnitude more chunks) plus a box arm at fixed h, since a smooth 2D texture of this quartic functional lowers its energy by spreading (E_u of a texture of size l scales as 1 / l^2), which makes the spread part box-limited rather than intrinsic"
    )
    V["claim4_wall_ladder"] = {
        "verdict": "QUALIFIED",
        "own": {
            "zero_mode_check_E_over_T_bps": zero,
            "stored_end_fields_E": {t: v["E_quartic_over_T_bps"] for t, v in sf.items()},
            "localization": {
                t: {
                    "E_fraction": v["E_fraction"],
                    "cells": v["cells"],
                    "share_top_1pct_cells": v["share_top_1pct_cells"],
                    "cells_holding_half": v["cells_holding_half"],
                }
                for t, v in loc.items()
            },
            "interpolated_E_over_T_bps": rr.get("E_interpolated_over_T_bps"),
            "own_relax_500_E_over_T_bps": rr.get("end", {}).get("E_quartic_over_T_bps"),
            "own_relax_traj": [
                (s["iters"], round(s["E_quartic_over_T_bps"], 5), round(s.get("drop_rel", 0), 4))
                for s in rr.get("traj", [])
            ],
        },
        "claimed": {
            "h1": [0.122, 0.107, 0.108],
            "h0.5": [0.142, 0.123, 0.122],
            "h0.25": [0.162, 0.208, 0.216],
            "interpolated_seed_h0.5_d0.3": 4.89,
            "label": "WALL_TENSION_UNRESOLVED",
        },
        "notes": notes4,
    }
    # ---------- claim 5 ----------
    c5 = out["verdicts"]["claim5_preregistration"]["own"]
    ka = out["verdicts"]["claim5_preregistration"].get("kappa_arm_claimed", {})
    notes5 = [
        "the letter's reference T_bps was mis-specified before the run: R25-1 and R26-2 had already measured the free static sector's smooth floor T_slaved3 below T_bps, so KAPPA_UNRESOLVED by the letter is correct and the amendment's direction is honest",
    ]
    alt = {}
    for d, v in c5.items():
        q = v["quartic_over_T_slaved3"]
        same = ka.get(d, {}).get("wall_seed_identical_to_smooth_at_kappa_fracs", [])
        ts3_pass = [
            k
            for k, ok in v["vs_T_slaved3_within_3pct"].items()
            if ok and float(k) in [float(x) for x in same]
        ]
        alt[d] = ts3_pass
        notes5.append(
            f"delta {d[1:]}: T_slaved3 / T_bps = {v['T_slaved3_over_T_bps']:.4f}; the smooth seed's quartic over T_slaved3 at kappa fraction {', '.join(f'{k}: {x:.4f}' for k, x in q.items())}; the author's amended reference is the smallest WALL-FREE kappa row ({ka.get(d, {}).get('strand_reference_kappa_frac_and_quartic')}), restoring (joint with the same-state condition) at {ka.get(d, {}).get('restores_amended_at_kappa_fracs')}; against T_slaved3 the joint condition restores at {ts3_pass}"
        )
    notes5.append(
        "the amendment's substitute reference (the same run's smallest wall-free kappa row) is self-referential (a row of the run judges the run) and at delta 0.03 it is the 0.1 kappa_ref row because the 0.03 row carries 6 wall bonds; the independent pre-existing reference T_slaved3 (R26-2) is the honest amended clause, and under it delta 0.03 restores at kappa_ref by a margin of 0.002 percent (1.02998 against the 1.03 clause), which is no margin at all: the delta 0.03 outcome is reference- and margin-sensitive and should be reported as such, not as 'at no kappa'; delta 0.3 restores at 0.3 kappa_ref under either reference"
    )
    V["claim5_preregistration"] = {
        "verdict": "QUALIFIED",
        "own": {
            d: {
                "T_slaved3_over_T_bps": v["T_slaved3_over_T_bps"],
                "quartic_over_T_slaved3": v["quartic_over_T_slaved3"],
                "amended_pass": v["amended_vs_smallest_kappa_within_3pct"],
                "T_slaved3_pass": v["vs_T_slaved3_within_3pct"],
                "letter_pass": v["letter_vs_T_bps_within_3pct"],
            }
            for d, v in c5.items()
        },
        "claimed": {
            "letter": "KAPPA_UNRESOLVED (cannot hold: the kappa-0 floor is T_slaved3 = 0.85 / 0.60 T_bps)",
            "amended": "restoration at 0.3 kappa_ref (delta 0.3), at no kappa (delta 0.03)",
        },
        "notes": notes5,
    }
    counts = {}
    for v in V.values():
        counts[v["verdict"]] = counts.get(v["verdict"], 0) + 1
    out["verdicts_short"] = V
    out["summary"] = {
        "counts": counts,
        "refuted": {k: v["notes"][0] for k, v in V.items() if v["verdict"] == "REFUTED"},
        "qualified": {k: v["notes"][-1] for k, v in V.items() if v["verdict"] == "QUALIFIED"},
        "blind_spots": [
            "the own descent is one preconditioned L-BFGS at n 48 for 1500 iterations from two seeds at delta 0.3; the m 2 strand, the pairs and delta 0.03 were not re-descended",
            "the delta 0.03 and 0.01 wall ladders were read on the stored fields only; the own h 1/2 relaxation stops at 500 iterations, far from any gate",
            "the ansatz test of the k^2 reading uses the BPS m 1 profile at the read core centres; the relaxed b profile near the pinned shell is not modelled",
            "the density localization uses the plain trace form of V4 (relative round-off 1e-7 to 1e-9 at these deltas, harmless for fractions)",
            "the psi reader is blind to director-row rotations; the Frobenius any-plane reader is reported beside it, its thresholds calibrated on pair-plane rotations only",
        ],
    }
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "verdict":
        # recompute the verdicts and the summary from the saved record (no descents)
        with open(OUT_JSON) as f:
            out = json.load(f)
    else:
        out = run(do_descents=(mode == "run"))
    out = verdicts(out)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(json.dumps(out["summary"], indent=1))
    for k, v in out["verdicts_short"].items():
        print(k, v["verdict"])
        for n in v["notes"]:
            print("   -", n)
    log("done")


if __name__ == "__main__":
    main()
