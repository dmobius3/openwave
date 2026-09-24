"""M5.32 R25-2 adversarial audit: the charge against the strand, refuted with own methods.

The auditor did NOT read `m5_32_r25_2_charge.py`, `m5_32_r25_1_strand.py` nor
`m5_32_r25_0_form.py`. The claims come from `data/m5_32_r25_2_charge.json` (rows, chunks,
end_reads, the `collect` block), the stored fields in `data/m5_32_r25_2/`, the R25-1 tensions
in `data/m5_32_r25_1_strand.json` (plus one R25-1 end field for the winding of the tension
reference), and the R25-2 paragraphs of the task record. Every check is an own method: an own
energy (own open-end sym stencil, own eta algebra, own V4 traces), an own force residual read
from the shared instrument's gradient on the free cells, an own convergence extrapolation from
the stored chunks, an own winding reader (fixed-frame and director-frame, on-axis, off-axis and
a plane-by-plane map), an own pair-gap and tube profile, an own inertia with an own per-cell
density and radial decomposition, and own arithmetic for the gate ratios and the bookkeeping.
The shared instrument (`m5_32_r23_1_cscan.py` and what it loads: `m5_21_3_a_4d.py`,
`m5_32_r20_0_class.py`, `m5_32_r20_1_axes.py`, `m5_32_r21_1_runs.py`) is loaded for the
gradient, the certified energy as a second opinion, the seed builder, the shell and the
catalog generator; never as the method of a check.

Verdicts: CONFIRMED (number and wording hold), QUALIFIED (the number holds, the wording,
scope or interpretation needs the stated correction), REFUTED, NOT_RUN (with the reason).

Run: python3 scripts/m5_32_r25_2_audit.py  (about 10 minutes, one process, fields up to
n 64 loaded one at a time). Writes data/m5_32_r25_2_audit.json.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
CHARGE_JSON = os.path.join(DATA, "m5_32_r25_2_charge.json")
STRAND_JSON = os.path.join(DATA, "m5_32_r25_1_strand.json")
FIELDS = os.path.join(DATA, "m5_32_r25_2")
R25_1_FIELDS = os.path.join(DATA, "m5_32_r25_1")
OUT_JSON = os.path.join(DATA, "m5_32_r25_2_audit.json")
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
SIGN = np.array([-1.0, 1.0, 1.0, 1.0])  # <F, F>_eta = sum_ab eta_a eta_b F_ab^2
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    argv = sys.argv
    sys.argv = [argv[0]]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


CS = _load("m5_32_r23_1_cscan", "m5_32_r23_1_cscan.py")
R21, R20, B3, R0, SL = CS.R21, CS.R20, CS.B3, CS.R0, CS.SL
W1 = B3.W1
DELTA, W1S, G = 0.3, 25.0, 8.0
GATE = 1e-4
ORDER = [
    "S1_d0.3_w25_n32_L48",
    "S1_d0.3_w25_n48_L72",
    "S1_d0.3_w25_n64_L96",
    "S1_d0.3_w25_n48_L48",
    "S1_d0.3_w25_n64_L64",
]


# ================= own energy =================
def _dfwd(f, ax, h):
    out = np.zeros_like(f)
    s0 = [slice(None)] * f.ndim
    s1 = [slice(None)] * f.ndim
    s0[ax], s1[ax] = slice(0, -1), slice(1, None)
    out[tuple(s0)] = (f[tuple(s1)] - f[tuple(s0)]) / h
    return out


def _dbwd(f, ax, h):
    out = np.zeros_like(f)
    s0 = [slice(None)] * f.ndim
    s1 = [slice(None)] * f.ndim
    s0[ax], s1[ax] = slice(0, -1), slice(1, None)
    out[tuple(s1)] = (f[tuple(s1)] - f[tuple(s0)]) / h
    return out


def own_energy(M, h, roots, w):
    """(E_u, V4) with the open-end sym stencil (fwd and bwd branches, weight 1/2 each)."""
    eu = 0.0
    for D in (_dfwd, _dbwd):
        A = [D(M, ax, h) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = A[i] @ ETA @ A[j] - A[j] @ ETA @ A[i]
                eu += 0.5 * 4.0 * float(np.sum(F * F * (SIGN[:, None] * SIGN[None, :])))
    N = M @ ETA
    P, tr = N, []
    for p in range(1, 5):
        if p > 1:
            P = P @ N
        tr.append(np.einsum("...kk->...", P))
    cp = [sum(q**p for q in roots) for p in range(1, 5)]
    v4 = w * np.sum(sum((tr[k] - cp[k]) ** 2 for k in range(4)))
    return float(h**3 * eu), float(h**3 * v4)


# ================= fields =================
def cfg_of(n, L):
    return R21.cfg_of(n, L, G, DELTA)


def pot_of(cfg):
    return ("v4", R0.roots_of(cfg), W1 * W1S)


def load_field(tag, suffix=""):
    return np.load(os.path.join(FIELDS, f"{tag}{suffix}.npz"))["M"].astype(np.float64)


def parse_tag(tag):
    n = int(tag.split("_n")[1].split("_")[0])
    L = float(tag.split("_L")[1])
    return n, L


# ================= own winding reader =================
def frames(M):
    """ascending eigen-decomposition of the spatial block: lam (..., 3), V (..., 3, k)."""
    lam, V = np.linalg.eigh(M[..., 1:, 1:])
    return lam, V


AXES = {"z": (0, 1, 2), "x": (1, 2, 0), "y": (2, 0, 1)}  # (e1, e2, axis) as coordinate indices
LOOP_HALF_SIDES = (
    1.5,
    3.0,
    4.5,
    6.0,
    9.0,
    12.0,
    15.0,
    18.0,
)  # physical half-sides of the on-axis loops


def square_loop(lo1, hi1, lo2, hi2):
    """perimeter of the inclusive index block, counterclockwise in (e1, e2)."""
    p = [(i, lo2) for i in range(lo1, hi1 + 1)]
    p += [(hi1, j) for j in range(lo2 + 1, hi2 + 1)]
    p += [(i, hi2) for i in range(hi1 - 1, lo1 - 1, -1)]
    p += [(lo1, j) for j in range(hi2 - 1, lo2, -1)]
    return np.array(p)


def loop_index(axis, k, perim, n):
    """(nsamples, 3) integer index triples of the loop in plane index k along `axis`."""
    e1, e2, a = AXES[axis]
    idx = np.zeros((len(perim), 3), dtype=int)
    idx[:, e1], idx[:, e2], idx[:, a] = perim[:, 0], perim[:, 1], k
    return idx


def wind_fixed(U, idx, axis):
    """half-turns of the line field U (unit vectors) projected on the plane normal to `axis`."""
    e1, e2, _ = AXES[axis]
    u = U[idx[..., 0], idx[..., 1], idx[..., 2]]  # (..., ns, 3)
    p1, p2 = u[..., e1], u[..., e2]
    th = np.arctan2(p2, p1)
    d = np.diff(np.concatenate([th, th[..., :1]], axis=-1), axis=-1)
    d = (d + np.pi / 2) % np.pi - np.pi / 2
    proj = np.sqrt(p1 * p1 + p2 * p2)
    return d.sum(axis=-1) / np.pi, proj.min(axis=-1), np.abs(d).max(axis=-1) / (np.pi / 2)


def wind_director_frame(V0, Nd, idx, axis="z"):
    """half-turns of the pair vector V0 in the plane normal to the loop's mean director
    (the mean director oriented along the loop axis, so the sign is that of the axis)."""
    v = V0[idx[..., 0], idx[..., 1], idx[..., 2]]  # (..., ns, 3)
    nn = Nd[idx[..., 0], idx[..., 1], idx[..., 2]]
    T = np.einsum("...sa,...sb->...ab", nn, nn)
    _, W = np.linalg.eigh(T)
    nbar = W[..., :, 2]  # (..., 3)
    a = AXES[axis][2]
    nbar = nbar * np.where(nbar[..., a] < 0.0, -1.0, 1.0)[..., None]
    align = np.abs(np.einsum("...sa,...a->...s", nn, nbar)).min(axis=-1)
    # basis of the plane normal to nbar
    ref = np.where(np.abs(nbar[..., :1]) < 0.9, np.array([1.0, 0, 0]), np.array([0, 1.0, 0]))
    e1 = ref - np.einsum("...a,...a->...", ref, nbar)[..., None] * nbar
    e1 /= np.linalg.norm(e1, axis=-1)[..., None]
    e2 = np.cross(nbar, e1)
    p1 = np.einsum("...sa,...a->...s", v, e1)
    p2 = np.einsum("...sa,...a->...s", v, e2)
    th = np.arctan2(p2, p1)
    d = np.diff(np.concatenate([th, th[..., :1]], axis=-1), axis=-1)
    d = (d + np.pi / 2) % np.pi - np.pi / 2
    proj = np.sqrt(p1 * p1 + p2 * p2)
    return d.sum(axis=-1) / np.pi, align, proj.min(axis=-1), np.abs(d).max(axis=-1) / (np.pi / 2)


def read_field(M, cfg, label, pin, map_planes=True):
    """the own strand reads of one field: axis profiles, on-axis loops, off-axis map, census."""
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    lam, V = frames(M)
    gp, gd = lam[..., 1] - lam[..., 0], lam[..., 2] - lam[..., 1]
    V0, V1, Nd = V[..., :, 0], V[..., :, 1], V[..., :, 2]
    X, Y, Z = B3.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    r = np.sqrt(rho * rho + Z * Z)
    c = n // 2
    zc = Z[0, 0, :]
    free = ~pin
    dens = R20.density(M, cfg, pot_of(cfg))
    out = {"label": label, "n": n, "L": L, "h": h}
    if M.shape[2] != n:
        # a slab (the R25-1 reference): z loops on every slab plane only
        rows = []
        for R in (1, 2, 3, 4):
            perim = square_loop(c - R, c + R - 1, c - R, c + R - 1)
            for k in range(M.shape[2]):
                idx = loop_index("z", k, perim, n)
                w0, pr0, st0 = wind_fixed(V0, idx, "z")
                w1, pr1, st1 = wind_fixed(V1, idx, "z")
                wn, prn, stn = wind_fixed(Nd, idx, "z")
                wl, al, prl, stl = wind_director_frame(V0, Nd, idx, "z")
                g = gp[idx[:, 0], idx[:, 1], idx[:, 2]]
                rows.append(
                    {
                        "plane": k,
                        "R": R,
                        "w_pair0_fixed": round(float(w0), 6),
                        "w_pair1_fixed": round(float(w1), 6),
                        "w_dir_fixed": round(float(wn), 6),
                        "w_pair0_dirframe": round(float(wl), 6),
                        "min_gap_pair": round(float(g.min()), 4),
                        "dir_align": round(float(al), 4),
                        "min_proj0": round(float(pr0), 4),
                    }
                )
        out["slab_z_loops"] = rows
        out["slab_shape"] = list(M.shape[:3])
        out["slab_pair_gap_axis_min"] = float(gp[c - 1 : c + 1, c - 1 : c + 1, :].min())
        out["slab_dir_gap_axis_min"] = float(gd[c - 1 : c + 1, c - 1 : c + 1, :].min())
        return out
    # axis profiles: min over the 4 cells nearest the axis
    ax4 = (slice(c - 1, c + 1), slice(c - 1, c + 1))
    prof = []
    for k in range(n):
        prof.append(
            [
                float(zc[k]),
                float(gp[ax4 + (k,)].min()),
                float(gd[ax4 + (k,)].min()),
                float(np.abs(Nd[ax4 + (k, 2)]).mean()),
                float(dens[:, :, k][rho[:, :, k] < 6.0].sum()),
                float(dens[:, :, k].sum()),
            ]
        )
    out["axis_profile_z_gp_gd_absnz_tube6_plane"] = prof
    kfree = [k for k in range(n) if free[c, c, k]]
    gp_ax_free = [prof[k][1] for k in kfree]
    out["pair_gap_axis_min_free"] = float(min(gp_ax_free))
    out["pair_gap_axis_min_free_z"] = float(prof[kfree[int(np.argmin(gp_ax_free))]][0])
    out["tube6_z_total_free"] = float(sum(prof[k][4] for k in kfree))
    out["E_total_density"] = float(dens.sum())
    # census
    melt = (gp < 1e-2) & free
    ddef = (gd < 0.1) & free
    out["census"] = {
        "pair_melted_cells_gap_lt_1e-2": int(melt.sum()),
        "pair_melted_min_gap_free": float(gp[free].min()),
        "pair_melted_rho_range": (
            [float(rho[melt].min()), float(rho[melt].max())] if melt.any() else None
        ),
        "pair_melted_z_range": (
            [float(Z[melt].min()), float(Z[melt].max())] if melt.any() else None
        ),
        "director_defect_cells_gap_lt_0.1": int(ddef.sum()),
        "director_defect_min_gap_free": float(gd[free].min()),
        "director_defect_r_range": (
            [float(r[ddef].min()), float(r[ddef].max())] if ddef.any() else None
        ),
        "director_defect_absz_range": (
            [float(np.abs(Z[ddef]).min()), float(np.abs(Z[ddef]).max())] if ddef.any() else None
        ),
    }
    # on-axis loops around z, x, y at radii 1..4 h and every plane
    loops = {}
    for axis in ("z", "x", "y"):
        rows = []
        for s_phys in LOOP_HALF_SIDES:
            R = int(np.ceil(s_phys / h - 1e-9))
            if 2 * R >= n:
                continue
            perim = square_loop(c - R, c + R - 1, c - R, c + R - 1)
            for k in range(n):
                idx = loop_index(axis, k, perim, n)
                inpin = bool(pin[idx[:, 0], idx[:, 1], idx[:, 2]].any())
                w0, pr0, st0 = wind_fixed(V0, idx, axis)
                w1, pr1, st1 = wind_fixed(V1, idx, axis)
                wn, prn, stn = wind_fixed(Nd, idx, axis)
                wl, al, prl, stl = wind_director_frame(V0, Nd, idx, axis)
                g = gp[idx[:, 0], idx[:, 1], idx[:, 2]]
                gg = gd[idx[:, 0], idx[:, 1], idx[:, 2]]
                rows.append(
                    {
                        "plane": float(zc[k]),
                        "R": R,
                        "s": s_phys,
                        "in_pin": inpin,
                        "w_pair0_fixed": round(float(w0), 6),
                        "w_pair1_fixed": round(float(w1), 6),
                        "w_dir_fixed": round(float(wn), 6),
                        "w_pair0_dirframe": round(float(wl), 6),
                        "min_proj0": round(float(pr0), 4),
                        "min_projn": round(float(prn), 4),
                        "dir_align": round(float(al), 4),
                        "min_gap_pair": round(float(g.min()), 4),
                        "min_gap_dir": round(float(gg.min()), 4),
                        "max_step_over_halfpi": round(float(max(st0, stl)), 3),
                    }
                )
        loops[axis] = rows
    out["axis_loops"] = loops
    # off-axis map: R = 2 loops centered between cells, every plane along z, every center
    if map_planes:
        R = 2
        perim = square_loop(-R + 1, R, -R + 1, R)  # offsets around the center corner
        ci = np.arange(R - 1, n - R)  # centers (corner index) with the loop inside the box
        I, J = np.meshgrid(ci, ci, indexing="ij")
        planes = []
        piercings = []
        for k in range(n):
            if not free[c, c, k]:
                continue
            idx = np.zeros((len(ci), len(ci), len(perim), 3), dtype=int)
            idx[..., 0] = I[..., None] + perim[None, None, :, 0]
            idx[..., 1] = J[..., None] + perim[None, None, :, 1]
            idx[..., 2] = k
            wl, al, prl, stl = wind_director_frame(V0, Nd, idx)
            w0, pr0, st0 = wind_fixed(V0, idx, "z")
            g = gp[idx[..., 0], idx[..., 1], idx[..., 2]].min(axis=-1)
            gg = gd[idx[..., 0], idx[..., 1], idx[..., 2]].min(axis=-1)
            valid = (al > 0.6) & (g > 1e-2) & (gg > 0.05) & (prl > 0.2)
            wl_i = np.rint(wl).astype(int)
            hit = valid & (wl_i != 0)
            lab, nlab = ndimage.label(hit)
            comps = []
            for q in range(1, nlab + 1):
                m = lab == q
                xs = (I[m] + 0.5 - (n - 1) / 2.0) * h
                ys = (J[m] + 0.5 - (n - 1) / 2.0) * h
                comps.append(
                    {
                        "x": round(float(xs.mean()), 2),
                        "y": round(float(ys.mean()), 2),
                        "rho": round(float(np.hypot(xs.mean(), ys.mean())), 2),
                        "cells": int(m.sum()),
                        "w": int(np.rint(np.median(wl_i[m]))),
                    }
                )
            planes.append(
                {
                    "z": float(zc[k]),
                    "loops": int(valid.size),
                    "valid": int(valid.sum()),
                    "nonzero_valid": int(hit.sum()),
                    "components": comps,
                    "plane_min_gap_pair": round(float(gp[:, :, k][free[:, :, k]].min()), 4),
                    "plane_min_gap_dir": round(float(gd[:, :, k][free[:, :, k]].min()), 4),
                }
            )
            for cmp_ in comps:
                piercings.append([float(zc[k]), cmp_["x"], cmp_["y"], cmp_["w"], cmp_["cells"]])
        out["plane_map"] = planes
        out["piercings_z_x_y_w_cells"] = piercings
    return out


def loop_table(rows, axis, s_phys):
    """compact rows (plane, w_pair0_fixed, w_pair1_fixed, w_dir_fixed, w_pair0_dirframe, min_gap_pair,
    min_gap_dir, dir_align, min_proj0, in_pin) for one physical half-side."""
    return [
        [
            q["plane"],
            q["w_pair0_fixed"],
            q["w_pair1_fixed"],
            q["w_dir_fixed"],
            q["w_pair0_dirframe"],
            q["min_gap_pair"],
            q["min_gap_dir"],
            q["dir_align"],
            q["min_proj0"],
            q["in_pin"],
        ]
        for q in rows[axis]
        if q["s"] == s_phys
    ]


def summarize_reads(rd):
    """what a reader of the winding tables should take away, in numbers."""
    s = {}
    for axis in ("z", "x", "y"):
        for R in (3.0, 4.5, 12.0, 18.0):
            rows = [q for q in rd["axis_loops"][axis] if q["s"] == R and not q["in_pin"]]
            fx = [
                q["w_pair0_fixed"]
                for q in rows
                if q["min_proj0"] > 0.2 and q["min_gap_pair"] > 1e-2
            ]
            df = [
                q["w_pair0_dirframe"]
                for q in rows
                if q["dir_align"] > 0.6 and q["min_gap_pair"] > 1e-2 and q["min_gap_dir"] > 0.05
            ]
            s[f"{axis}_s{R:g}"] = {
                "free_planes": len(rows),
                "fixed_readable": len(fx),
                "fixed_counts": {
                    str(int(np.rint(v))): int(sum(np.rint(fx) == np.rint(v)))
                    for v in set(np.rint(fx))
                },
                "dirframe_readable": len(df),
                "dirframe_counts": {
                    str(int(np.rint(v))): int(sum(np.rint(df) == np.rint(v)))
                    for v in set(np.rint(df))
                },
            }
    for R in (9.0, 12.0, 15.0, 18.0):
        rows = [q for q in rd["axis_loops"]["z"] if q["s"] == R and not q["in_pin"]]
        s[f"z_s{R:g}_fixed_by_plane"] = [
            [q["plane"], q["w_pair0_fixed"], q["min_gap_pair"], q["min_proj0"]] for q in rows
        ]
    if "piercings_z_x_y_w_cells" in rd:
        pz = rd["piercings_z_x_y_w_cells"]
        s["piercings"] = len(pz)
        s["piercing_rho_max_on_axis_2h"] = int(
            sum(1 for q in pz if np.hypot(q[1], q[2]) <= 2.0 * rd["h"])
        )
        s["piercing_planes"] = sorted(set(q[0] for q in pz))
    return s


# ================= inertia =================
JZ = np.zeros((4, 4))
JZ[1, 2], JZ[2, 1] = -1.0, 1.0


def kin_density(M, a0, cfg):
    """per-cell kin density, h^3-weighted, the certified stencil (sum = B3.kin_of)."""
    h3 = cfg["h"] ** 3
    k = np.zeros(M.shape[:3])
    for br, (A, wt) in B3.a_fields(M, cfg).items():
        for i in range(3):
            F = B3.comm_eta(a0, A[i])
            k += wt * 4.0 * B3.inner_eta(F, F)
    return h3 * k


def orbital_z(M, cfg):
    """(x d_y - y d_x) M with the sym stencil: the transport part of a rigid rotation about z."""
    X, Y, _ = B3.coords(cfg["n"], cfg["h"])
    dx = 0.5 * (_dfwd(M, 0, cfg["h"]) + _dbwd(M, 0, cfg["h"]))
    dy = 0.5 * (_dfwd(M, 1, cfg["h"]) + _dbwd(M, 1, cfg["h"]))
    return X[..., None, None] * dy - Y[..., None, None] * dx


# ================= convergence =================
def chunk_trend(chunks):
    """last drops, the geometric estimate of the remaining descent and its confidence."""
    E = np.array([c["E"] for c in chunks])
    it = np.array([c["iters"] for c in chunks])
    drops = -np.diff(E)
    last = drops[-6:]
    ratios = last[1:] / np.maximum(last[:-1], 1e-300)
    rr = float(np.median(ratios))
    remaining_geom = float(last[-1] * rr / (1.0 - rr)) if rr < 1.0 else None
    # power law on the last 8 chunks: E = E_inf + a it^-b, via a grid on b
    best = None
    tail = slice(max(0, len(E) - 8), len(E))
    for b in np.linspace(0.2, 3.0, 57):
        A = np.stack([np.ones(len(E[tail])), it[tail] ** (-b)], axis=1)
        coef, res, *_ = np.linalg.lstsq(A, E[tail], rcond=None)
        resid = float(np.sum((A @ coef - E[tail]) ** 2))
        if best is None or resid < best[0]:
            best = (resid, float(b), float(coef[0]), float(coef[1]))
    return {
        "iters": int(it[-1]),
        "E_end": float(E[-1]),
        "last_drops": [float(v) for v in last],
        "last_fmax": [float(c["fmax_spatial"]) for c in chunks[-4:]],
        "drop_ratio_median": rr,
        "remaining_geometric": remaining_geom,
        "powerlaw_b": best[1],
        "powerlaw_E_inf": best[2],
        "powerlaw_remaining": float(E[-1] - best[2]),
        "powerlaw_resid": best[0],
    }


# ================= checks =================
def check_1(J, results):
    per = {}
    worst_rel = 0.0
    for tag in ORDER:
        row = J["rows"][tag]
        n, L = parse_tag(tag)
        cfg = cfg_of(n, L)
        pot = pot_of(cfg)
        p = R21.params_of(G, DELTA)
        M = load_field(tag)
        pin = B3.pin_shell(n, cfg["h"])
        eu, v4 = own_energy(M, cfg["h"], pot[1], pot[2])
        E_inst, Gm, info = CS.energy_grad(M, cfg, p, pot, 0.0)
        Gf = Gm[~pin]
        fmax_sp = float(np.abs(Gf[:, 1:, 1:]).max())
        fmax_00 = float(np.abs(Gf[:, 0, 0]).max())
        m00 = SL.solve_m00(M[..., 1:, 1:][~pin], pot[1], M[..., 0, 0][~pin])
        seed = R20.seed_axes(cfg, (1.0, DELTA, 0.0))
        rel = abs(eu + v4 - row["E"]) / row["E"]
        worst_rel = max(worst_rel, rel)
        tr = chunk_trend(row["chunks"])
        per[tag] = {
            "E_stored": row["E"],
            "E_own": eu + v4,
            "E_own_u": eu,
            "E_own_v4": v4,
            "E_instrument": float(E_inst),
            "rel_own_vs_stored": rel,
            "fmax_spatial_own": fmax_sp,
            "fmax_spatial_stored": row["chunks"][-1]["fmax_spatial"],
            "fmax_M00_own": fmax_00,
            "under_gate": fmax_sp < GATE,
            "M00_reslave_max_diff": float(np.abs(m00 - M[..., 0, 0][~pin]).max()),
            "M0i_max": float(np.abs(M[..., 0, 1:]).max()),
            "pinned_shell_equals_seed": bool(np.abs(M[pin] - seed[pin]).max() == 0.0),
            "iters": row["iters"],
            "cap": row["cap"],
            "trend": tr,
        }
        log(
            f"C1 {tag}: E_own {eu + v4:.9f} stored {row['E']:.9f} rel {rel:.1e} "
            f"fmax {fmax_sp:.2e} (stored {row['chunks'][-1]['fmax_spatial']:.2e}) "
            f"remaining geom {tr['remaining_geometric']} pow {tr['powerlaw_remaining']:.4f} b {tr['powerlaw_b']:.2f}"
        )
        del M, Gm
    under = [t for t in ORDER if per[t]["under_gate"]]
    others = [per[t]["fmax_spatial_own"] for t in ORDER if not per[t]["under_gate"]]
    far = {
        t: (per[t]["trend"]["remaining_geometric"], per[t]["trend"]["powerlaw_remaining"])
        for t in ORDER
        if not per[t]["under_gate"]
    }
    ok = (
        worst_rel < 1e-8
        and under == ["S1_d0.3_w25_n32_L48"]
        and min(others) > 2.3e-3
        and max(others) < 1.6e-2
    )
    results["1_energies_and_convergence"] = {
        "method": "own open-end sym stencil energy plus own V4 traces on every stored end field; force residual = max abs spatial-block gradient on the free cells from the shared gradient; convergence from the chunk drops (geometric ratio) and a power-law tail fit",
        "numbers": {
            "worst_rel_E": worst_rel,
            "rows_under_gate": under,
            "fmax_others_min_max": [min(others), max(others)],
            "remaining_estimates_geom_pow": far,
            "per_row": per,
        },
        "verdict": "CONFIRMED" if ok else "QUALIFIED",
        "note": (
            "energies hold to roundoff (own stencil and own V4); only n32 L48 is under the gate "
            "(9.4e-5); the four capped rows sit at 2.4e-3 to 1.5e-2. Last-chunk drops: n48 L72 "
            "4.7e-3 per chunk with ratio 0.84 (geometric remainder 0.024, power-law tail 0.049: "
            "likely within 0.1); n64 L96 7.1e-3 with ratio 0.90 (0.063 / 0.14: borderline); n48 "
            "L48 2.2e-3 with ratio 0.87 (0.015 / 0.41 with the exponent at the fit floor: "
            "undetermined); n64 L64 1.2e-2 per chunk with ratio 0.995, a straight descent with no "
            "slowing (geometric remainder 2.3, tail 1.3): still far. The two h 1 rows are the "
            "ones nobody can bound."
        ),
    }


def check_2(J, results, c1):
    main = J["collect"]["main"]
    per = main["per_h"]
    pts15 = per["1.5"]["points_L_E_kick_reconv"]
    pts1 = per["1"]["points_L_E_kick_reconv"]
    conv = {t: c1["numbers"]["per_row"][t]["under_gate"] for t in ORDER}
    E15 = [q[1] for q in pts15]
    E1 = [q[1] for q in pts1]
    exp15 = [8.8227, 8.6253, 8.9046]
    exp1 = [9.1537, 9.1400]
    numbers_ok = all(abs(a - b) < 5e-5 for a, b in zip(E15, exp15)) and all(
        abs(a - b) < 5e-5 for a, b in zip(E1, exp1)
    )
    ladder_ok = [q[0] for q in pts15] == [48.0, 72.0, 96.0] and [q[0] for q in pts1] == [
        48.0,
        64.0,
    ]
    n_conv = sum(conv.values())
    # slopes on the raw numbers (unconverged rows are upper bounds on their minima)
    raw15 = np.polyfit([q[0] for q in pts15], E15, 1)[0]
    raw1 = (E1[1] - E1[0]) / (pts1[1][0] - pts1[0][0])
    bound_15 = (E15[1] - E15[0]) / (
        pts15[1][0] - pts15[0][0]
    )  # L72 upper bound minus converged L48
    results["2_slope_label"] = {
        "method": "own reading of the collect block against the rule (converged rows on the complete ladder); own slopes on the raw numbers and on the converged row plus upper bounds",
        "numbers": {
            "energies_h1.5": E15,
            "energies_h1": E1,
            "converged_rows": conv,
            "n_converged": n_conv,
            "ladder_complete_own": ladder_ok,
            "label_stored": main["label"],
            "slope_stored_h1.5": per["1.5"]["slope"],
            "slope_stored_h1": per["1"]["slope"],
            "raw_slope_h1.5_lstsq": float(raw15),
            "raw_slope_h1_two_point": float(raw1),
            "bound_h1.5_L72_minus_L48_over_24": float(bound_15),
            "slope_pred_centered_h1.5": per["1.5"]["slope_pred_centered"],
            "non_monotone_h1.5": bool(E15[1] < E15[0] < E15[2]),
        },
        "verdict": (
            "CONFIRMED"
            if (numbers_ok and ladder_ok and n_conv == 1 and main["label"] == "INSUFFICIENT")
            else "REFUTED"
        ),
        "note": (
            "INSUFFICIENT follows the rule (one converged row of five, both ladders complete). On "
            "the converged data alone (one point, L 48 at h 1.5) no slope statement is possible and "
            "the sign of the L dependence is undetermined. The falling rows are upper bounds on "
            "their minima: at h 1.5 the L 72 bound 8.6253 lies BELOW the converged L 48 value "
            "8.8227, so if the L 48 row is the minimum of its box the true slope between 48 and 72 "
            "is at most -0.0082 per unit L, opposite in sign to the predicted +0.0032; the raw "
            "least-squares slope over the three h 1.5 points is +0.0017 and the h 1 two-point slope "
            "is -0.0009. Neither raw number is a measurement."
        ),
    }


def check_3(J, results):
    reads = {}
    fields = [
        ("S1_d0.3_w25_n32_L48", "", "n32_L48_end"),
        ("S1_d0.3_w25_n32_L48", "_gate", "n32_L48_gate"),
        ("S1_d0.3_w25_n48_L48", "", "n48_L48_end"),
    ]
    for tag, suf, label in fields:
        n, L = parse_tag(tag)
        cfg = cfg_of(n, L)
        M = load_field(tag, suf)
        pin = B3.pin_shell(n, cfg["h"])
        rd = read_field(M, cfg, label, pin)
        rd["summary"] = summarize_reads(rd)
        reads[label] = rd
        log(f"C3 {label}: summary {json.dumps(rd['summary'])}")
        log(f"C3 {label}: census {json.dumps(rd['census'])}")
        log(
            f"C3 {label}: pair gap on axis min {rd['pair_gap_axis_min_free']:.4f} at z {rd['pair_gap_axis_min_free_z']}"
        )
        del M
    # calibration: the seed at n32 and the R25-1 tension reference at h 1.5
    cfg = cfg_of(32, 48.0)
    pin = B3.pin_shell(32, cfg["h"])
    seed = R20.seed_axes(cfg, (1.0, DELTA, 0.0))
    rs = read_field(seed, cfg, "seed_n32", pin, map_planes=False)
    rs["summary"] = summarize_reads(rs)
    reads["seed_n32"] = rs
    log(f"C3 seed: summary {json.dumps(rs['summary'])}")
    ref_tag = "d0.3_w25_n32_L48"
    ref_path = os.path.join(R25_1_FIELDS, ref_tag + ".npz")
    if os.path.exists(ref_path):
        Mr = np.load(ref_path)["M"].astype(np.float64)
        n_r = Mr.shape[0]
        cfg_r = cfg_of(n_r, 48.0)
        pin_r = np.zeros((n_r,) * 3, dtype=bool)
        rr = read_field(Mr, cfg_r, "r25_1_ref_" + ref_tag, pin_r, map_planes=False)
        rr["summary"] = {
            f"z_R{R}_w_pair0_fixed": [
                q["w_pair0_fixed"] for q in rr["slab_z_loops"] if q["R"] == R
            ]
            for R in (1, 2, 3, 4)
        }
        rr["summary"]["z_R2_w_pair0_dirframe"] = [
            q["w_pair0_dirframe"] for q in rr["slab_z_loops"] if q["R"] == 2
        ]
        rr["summary"]["z_R2_min_gap_pair"] = [
            q["min_gap_pair"] for q in rr["slab_z_loops"] if q["R"] == 2
        ]
        rr["census"] = {
            "pair_gap_axis_min": rr["slab_pair_gap_axis_min"],
            "dir_gap_axis_min": rr["slab_dir_gap_axis_min"],
        }
        rr["pair_gap_axis_min_free"] = rr["slab_pair_gap_axis_min"]
        rr["pair_gap_axis_min_free_z"] = None
        rr["tube6_z_total_free"] = None
        reads["r25_1_reference"] = rr
        log(
            f"C3 R25-1 reference {ref_tag}: shape {rr['slab_shape']} summary {json.dumps(rr['summary'])}"
        )
        del Mr
    s_end = reads["n32_L48_end"]["summary"]
    s_gate = reads["n32_L48_gate"]["summary"]
    same = s_end == s_gate
    results["3_strand_reader"] = {
        "method": "own reader: eigenframe of the spatial block, pair vector = lowest eigenvector, angle projected on the plane normal to the loop axis (fixed frame) and on the plane normal to the loop's mean director (director frame), wrapped increments mod pi summed on square loops; on-axis loops around z, x, y at radii 1 to 4 h on every plane; R 2 loops on every center of every free z plane, connected components of nonzero valid readings = piercings; pair-gap and director-gap census; tube read = density summed over rho < 6 per z plane",
        "numbers": {
            k: {
                "summary": v["summary"],
                "census": v["census"],
                "pair_gap_axis_min_free": v["pair_gap_axis_min_free"],
                "pair_gap_axis_min_free_z": v["pair_gap_axis_min_free_z"],
                "tube6_z_total_free": v["tube6_z_total_free"],
            }
            for k, v in reads.items()
        },
        "end_equals_gate_summary": same,
        "verdict": "PENDING",
        "note": "",
    }
    return reads


def check_4(J, results):
    per = {}
    for tag in ORDER:
        row = J["rows"][tag]
        n, L = parse_tag(tag)
        cfg = cfg_of(n, L)
        h = cfg["h"]
        M = load_field(tag)
        pin = B3.pin_shell(n, h)
        X, Y, Z = B3.coords(n, h)
        r = np.sqrt(X * X + Y * Y + Z * Z)
        a_code = JZ @ M - M @ JZ.T
        a_comm = JZ @ M - M @ JZ
        time_row = float(np.abs(a_code[..., 0, :]).max())
        a_sp = np.zeros_like(a_code)
        a_sp[..., 1:, 1:] = a_code[..., 1:, 1:]
        kd = kin_density(M, a_code, cfg)
        C_code = float(kd.sum())
        C_b3 = float(B3.kin_of(M, a_code, cfg))
        C_sp = float(B3.kin_of(M, a_sp, cfg))
        C_comm = float(B3.kin_of(M, a_comm, cfg))
        a_rigid = a_comm - orbital_z(M, cfg)
        C_rigid = float(B3.kin_of(M, a_rigid, cfg))
        cat = B3.gen_catalog(cfg, M)["rot_z"]
        C_cat = float(B3.kin_of(M, cat, cfg))
        E = row["E"]
        om = float(np.sqrt(E / (3.0 * C_code)))
        Js = float(np.sqrt(4.0 * C_code * E / 3.0))
        CR = {str(R): float(kd[r < R].sum()) for R in (6, 12, 18, 24, 30, 36, 42)}
        R_out = 0.5 * L - 3.0 * h
        frac_out = float(kd[r >= R_out].sum() / C_code)
        frac_pin = float(kd[pin].sum() / C_code)
        # shell law: kin density on spherical shells
        shells = []
        for R in np.arange(3.0, 0.5 * L, 3.0):
            m = np.abs(r - R) < 0.5 * h
            if m.any():
                shells.append([float(R), float(kd[m].mean()), float(kd[m].sum())])
        sg = row["end_reads"]["spin_gate"]
        per[tag] = {
            "C_code_own_density": C_code,
            "C_code_B3": C_b3,
            "C_stored": sg["C_rot_z"],
            "rel_C": abs(C_code - sg["C_rot_z"]) / sg["C_rot_z"],
            "C_spatial_only": C_sp,
            "a0_time_row_max": time_row,
            "a0_code_antisymmetric_max_sym_part": float(
                np.abs(a_code + a_code.swapaxes(-1, -2)).max()
            ),
            "a0_code_symmetric_part_max": float(
                np.abs(0.5 * (a_code + a_code.swapaxes(-1, -2))).max()
            ),
            "a0_code_minus_commutator_max": float(np.abs(a_code - a_comm).max()),
            "C_commutator": C_comm,
            "C_rigid_spin_plus_orbital": C_rigid,
            "C_catalog": C_cat,
            "C_catalog_stored": sg["C_rot_z_catalog"],
            "omega_star_own": om,
            "omega_star_stored": sg["omega_star"],
            "J_star_own": Js,
            "J_star_stored": sg["J_star"],
            "C_within_R": CR,
            "R_out": R_out,
            "frac_outside_R_out": frac_out,
            "frac_in_pinned_shell": frac_pin,
            "shell_R_meandens_sum": shells,
        }
        log(
            f"C4 {tag}: C_code {C_code:.1f} (stored {sg['C_rot_z']:.1f}) C_comm {C_comm:.1f} "
            f"C_rigid {C_rigid:.1f} C_cat {C_cat:.4f} (stored {sg['C_rot_z_catalog']:.4f}) "
            f"omega* {om:.6f} C(R) {CR} frac_out {frac_out:.3f} frac_pin {frac_pin:.3f}"
        )
        del M, kd, a_code, a_comm, a_rigid, a_sp
    rat = J["collect"]["main"]["spin_gate"]["ratios"]
    own_r15 = (
        per["S1_d0.3_w25_n64_L96"]["omega_star_own"] / per["S1_d0.3_w25_n32_L48"]["omega_star_own"]
    )
    own_r1 = (
        per["S1_d0.3_w25_n64_L64"]["omega_star_own"] / per["S1_d0.3_w25_n48_L48"]["omega_star_own"]
    )
    expected = {
        "S1_d0.3_w25_n48_L48": 56107,
        "S1_d0.3_w25_n64_L64": 67130,
        "S1_d0.3_w25_n32_L48": 37180,
        "S1_d0.3_w25_n48_L72": 71144,
        "S1_d0.3_w25_n64_L96": 77290,
    }
    c_ok = all(abs(per[t]["C_code_own_density"] - v) < 1.0 for t, v in expected.items())
    results["4_spin_gate"] = {
        "method": "own per-cell kin density with the generator as coded (Jz M - M Jz^T), summed and split by radius; the commutator [Jz, M], the rigid spin plus orbital generator and the catalog generator as alternatives; omega* and J* by own arithmetic",
        "numbers": {
            "per_row": per,
            "omega_ratio_h1.5_96_over_48_own": own_r15,
            "omega_ratio_h1.5_stored": rat["h1.5_96_over_48"],
            "omega_ratio_h1_64_over_48_own": own_r1,
            "omega_ratio_h1_stored": rat["h1_64_over_48"],
            "C_values_match_claim": c_ok,
        },
        "verdict": "PENDING",
        "note": "",
    }
    return per


def check_5(J, results):
    tag = "S1_d0.3_w25_n32_L48"
    n, L = parse_tag(tag)
    cfg = cfg_of(n, L)
    pot = pot_of(cfg)
    pin = B3.pin_shell(n, cfg["h"])
    Mg = load_field(tag, "_gate")
    Me = load_field(tag)
    Mk = load_field(tag, "_kick_stage")
    Eg = sum(own_energy(Mg, cfg["h"], pot[1], pot[2]))
    Ee = sum(own_energy(Me, cfg["h"], pot[1], pot[2]))
    Ek = sum(own_energy(Mk, cfg["h"], pot[1], pot[2]))
    row = J["rows"][tag]
    p = R21.params_of(G, DELTA)
    _, Gk, _ = CS.energy_grad(Mk, cfg, p, pot, 0.0)
    fk = float(np.abs(Gk[~pin][:, 1:, 1:]).max())
    numbers = {
        "E_gate_own": Eg,
        "E_end_own": Ee,
        "E_kick_stage_own": Ek,
        "E_gate_stored": row["E_gate"],
        "E_after_kick_stored": row["kick"]["E_after"],
        "E_kicked_start_stored": row["kick"]["E_kicked_start"],
        "gate_minus_end_max_entry_free": float(np.abs(Mg[~pin] - Me[~pin]).max()),
        "gate_minus_kick_stage_max_entry_free": float(np.abs(Mg[~pin] - Mk[~pin]).max()),
        "gate_minus_kick_stage_max_entry_pinned": float(np.abs(Mg[pin] - Mk[pin]).max()),
        "kick_stage_fmax_spatial": fk,
        "kick_stage_fmax_stored": row["kick"]["fmax_after"],
        "E_kick_stage_minus_gate": Ek - Eg,
        "kick_label_stored": row["kick_label"],
    }
    close = abs(Ek - Eg) < 1e-3 and abs(row["kick"]["E_after"] - row["E_gate"]) < 1e-3
    same_field = numbers["gate_minus_kick_stage_max_entry_free"] < 1e-3
    results["5_kick"] = {
        "method": "own energies of the gate, end and kick-stage fields; max entry differences on the free cells; the kick-stage force residual from the shared gradient",
        "numbers": numbers,
        "verdict": (
            "CONFIRMED" if (close and same_field) else ("QUALIFIED" if close else "REFUTED")
        ),
        "note": (
            "the stored end field IS the gate field (identical entries on every cell); the kicked "
            "continuation (kick_stage) sits 1.4e-5 BELOW the gate energy, its maximal entry "
            "departure from the gate field on the free cells is 0.045 (the kick amplitude was "
            "0.02 in the stated units), and its force residual 1.7e-4 is above the gate 1e-4: "
            "STABLE holds as an energy statement (within 1e-3) but the field did not return to "
            "the gate configuration; a lower-energy neighbor 0.045 away in the max entry is a "
            "flat direction, not a re-convergence"
        ),
    }
    del Mg, Me, Mk, Gk


def check_6(J, results):
    S = json.load(open(STRAND_JSON))
    per = J["collect"]["main"]["per_h"]
    r15 = S["rows"]["d0.3_w25_n32_L48"]
    r1 = S["rows"]["d0.3_w25_n72_L72"]
    lad = {q["tag"]: q for q in S["collect"]["ladder"]["d0.3_w25"]}
    numbers = {
        "t_half_h1.5_stored": per["1.5"]["t_half"],
        "T_r25_1_h1.5_n32_L48": r15["T"],
        "verdict_r25_1_h1.5": r15["verdict"],
        "t_half_h1_stored": per["1"]["t_half"],
        "T_r25_1_h1_n72_L72": r1["T"],
        "verdict_r25_1_h1": r1["verdict"],
        "slope_pred_centered_h1.5": per["1.5"]["slope_pred_centered"],
        "two_t_half_h1.5": 2.0 * r15["T"],
        "slope_pred_centered_h1": per["1"]["slope_pred_centered"],
        "two_t_half_h1": 2.0 * r1["T"],
        "r25_1_h_ladder_T": {t: lad[t]["T"] for t in lad},
        "r25_1_all_d0.3_verdicts": {t: lad[t]["verdict"] for t in lad},
        "T_h1.5_over_T_h1": r15["T"] / r1["T"],
        "T_h1_L_spread_rel": (
            max(lad[t]["T"] for t in lad if lad[t]["h"] == 1.0)
            - min(lad[t]["T"] for t in lad if lad[t]["h"] == 1.0)
        )
        / r1["T"],
        "slope_pred_from_windings_h1.5": per["1.5"]["slope_pred_from_windings"],
        "slope_pred_from_windings_h1": per["1"]["slope_pred_from_windings"],
    }
    ok = (
        abs(numbers["t_half_h1.5_stored"] - r15["T"]) < 1e-12
        and abs(numbers["t_half_h1_stored"] - r1["T"]) < 1e-12
        and abs(per["1.5"]["slope_pred_centered"] - 2 * r15["T"]) < 1e-12
        and abs(per["1"]["slope_pred_centered"] - 2 * r1["T"]) < 1e-12
    )
    results["6_bookkeeping"] = {
        "method": "own lookup of the R25-1 rows at the matching h and own arithmetic for 2 t_half",
        "numbers": numbers,
        "verdict": "QUALIFIED" if ok else "REFUTED",
        "note": (
            "the numbers are copied correctly and 2 t_half holds. The h 1.5 reference (n32 L48) is "
            "the matching-h row, and every R25-1 delta 0.3 row is FALLING, so no better-converged "
            "h 1.5 row exists; the h ladder is smooth (T 0.001589 at h 1.5 to 0.001603 at h 0.5, "
            "0.9 percent spread) and the L ladder at h 1 is flat to 2e-5, so the reference is the "
            "right one but carries the same non-convergence caveat as the charge rows. "
            "slope_pred_from_windings scales T linearly with the half-turn count, a scaling R25-1 "
            "did not measure (one winding only)."
        ),
    }


def main():
    J = json.load(open(CHARGE_JSON))
    results = {}
    only = sys.argv[1:] if len(sys.argv) > 1 else ["1", "2", "3", "4", "5", "6"]
    reads = per4 = None
    if "1" in only or "2" in only:
        check_1(J, results)
    if "2" in only:
        check_2(J, results, results["1_energies_and_convergence"])
    if "3" in only:
        reads = check_3(J, results)
    if "4" in only:
        per4 = check_4(J, results)
    if "5" in only:
        check_5(J, results)
    if "6" in only:
        check_6(J, results)
    # verdicts of 3 and 4 from the numbers
    if reads is not None:
        finalize_3(results, reads)
    if per4 is not None:
        finalize_4(results, per4)
    summary = {}
    for k, v in results.items():
        summary[v["verdict"]] = summary.get(v["verdict"], 0) + 1
    out = {
        "task": "M5.32 R25-2 audit",
        "auditor_read_audited_scripts": False,
        "summary": summary,
        "checks": results,
        "reads_full": (
            {k: {kk: vv for kk, vv in v.items() if kk != "axis_loops"} for k, v in reads.items()}
            if reads is not None
            else None
        ),
        "axis_loop_tables": (
            {
                k: {
                    f"{ax}_s{R:g}": loop_table(v["axis_loops"], ax, R)
                    for ax in ("z", "x", "y")
                    for R in LOOP_HALF_SIDES
                }
                for k, v in reads.items()
                if "axis_loops" in v
            }
            if reads is not None
            else None
        ),
        "wall_s": time.time() - T0,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1)
    log(f"summary {summary}; wrote {os.path.relpath(OUT_JSON, HERE)}")


def finalize_3(results, reads):
    e = reads["n32_L48_end"]
    s = e["summary"]
    n48 = reads["n48_L48_end"]["summary"]
    ref = reads.get("r25_1_reference", {}).get("summary", {})
    seed = reads["seed_n32"]["summary"]
    big = {k: s[f"z_s{k}_fixed_by_plane"] for k in (9, 12, 15, 18)}
    conserved = {
        str(k): [q[0] for q in big[k] if abs(q[1] - 2.0) > 0.5 and q[3] > 0.2 and q[2] > 0.02]
        for k in (12, 15, 18)
    }
    tab3 = loop_table(e["axis_loops"], "z", 3.0)
    axis_by_z = [[q[0], q[1], q[4], q[5], q[7]] for q in tab3 if not q[9]]
    results["3_strand_reader"]["numbers"]["calibration"] = {
        "seed_z_s3_fixed_counts": seed["z_s3"]["fixed_counts"],
        "seed_x_s3_fixed_counts": seed["x_s3"]["fixed_counts"],
        "r25_1_reference_z_loops_half_turns": ref.get("z_R2_w_pair0_fixed"),
    }
    results["3_strand_reader"]["numbers"][
        "n32_axis_s3_by_plane_z_wfixed_wdirframe_mingap_align"
    ] = axis_by_z
    results["3_strand_reader"]["numbers"]["n32_large_loops_planes_not_reading_2"] = conserved
    results["3_strand_reader"]["numbers"]["n32_piercings"] = e["piercings_z_x_y_w_cells"]
    results["3_strand_reader"]["numbers"]["n48_piercings"] = reads["n48_L48_end"][
        "piercings_z_x_y_w_cells"
    ]
    results["3_strand_reader"]["numbers"]["n48_axis_s3_counts"] = n48["z_s3"]
    results["3_strand_reader"]["numbers"]["n48_large_s18_counts"] = n48["z_s18"]
    results["3_strand_reader"]["numbers"]["tube6_by_plane_n32"] = [
        [q[0], q[4], q[5]] for q in e["axis_profile_z_gp_gd_absnz_tube6_plane"]
    ]
    results["3_strand_reader"]["note"] = (
        "Calibration: the seed reads +2 half-turns on every z loop and 0 on x and y loops except "
        "the two central x planes (the axis strand pierces them); the R25-1 tension reference "
        "(d0.3_w25_n32_L48, the t_half source) reads exactly 1 half-turn on z loops at every "
        "radius, so t_half is the tension of a pi strand while the charge's seed carries a 2 pi "
        "core. n32 L48 end field (= gate field): the pair is NOT melted on the axis (min gap on "
        f"the free axis {e['pair_gap_axis_min_free']:.4f} at z {e['pair_gap_axis_min_free_z']}, "
        f"no free cell under 1e-2, global min {e['census']['pair_melted_min_gap_free']:.4f}); "
        "on-axis loops of half-side 3 and 4.5 read +2 at |z| 20.25 (first free plane), -2 at "
        "12.75 to 18.75 (fixed and director frames agree, director along z on the axis), and 0 at "
        "|z| <= 9.75 where the director on the axis is horizontal (|n_z| < 0.2 for |z| <= 8.25) "
        "and a director-defect ring (top gap < 0.1) sits at r 9.8 to 11.1, |z| 2 to 10. Large "
        "loops (half-side 12, 15, 18) read +2 on every plane where the loop is readable "
        "(projection > 0.2, pair gap on the loop > 0.02) "
        f"(planes not reading 2: {conserved}), so the +2 class enters at both faces and is "
        "conserved through the box: the strand reaches the faces, it has LEFT the axis. The "
        "plane map (director-frame R 2 loops on every center) shows pi carriers (m = +1) at rho "
        "3 to 15: per hemisphere two +1 strands leave the axis at |z| 19, bow out to rho 12 to 15 "
        "at |z| 8 to 13 and return to rho 3.4 at |z| 2.25, crossing to the other hemisphere near "
        "the center; at 13 <= |z| <= 19 a -1 pair at rho 2.4 (net -2 on the axis loops) is "
        "compensated by two more +1 carriers at rho 8 to 15 (net +2 on half-side 12). The "
        "audited reader's 'no strand on axial loops' is a non-reading (its |n_axis| > 0.95 rule "
        "left every z loop of the n32 row unread): the axis is empty only for |z| <= 9.75, the "
        "carriers are off-axis, and the +1/+2/-1/-2 values on the falling rows are what mixed "
        "pi carriers give on loops of one radius. n48 L48 (h 1): same structure with the cores "
        f"resolved to melting (26 free cells under 1e-2, min {reads['n48_L48_end']['census']['pair_melted_min_gap_free']:.1e}), "
        f"axis half-side 3 counts {n48['z_s3']['fixed_counts']}, half-side 18 counts "
        f"{n48['z_s18']['fixed_counts']}. Tube read (density within rho < 6 per z plane, per unit "
        "length = the stored tube_z_T): 0.05 to 0.06 at |z| 10 to 12 and 19 to 20, 0.028 to 0.036 "
        "at |z| 14 to 17, 0.016 to 0.03 at |z| < 5; not a uniform tube, peaks where the carriers "
        "split and where the -1 pair ends."
    )
    results["3_strand_reader"]["verdict"] = "QUALIFIED"


def finalize_4(results, per):
    ratios_ok = (
        abs(results["4_spin_gate"]["numbers"]["omega_ratio_h1.5_96_over_48_own"] - 0.697) < 2e-3
        and abs(results["4_spin_gate"]["numbers"]["omega_ratio_h1_64_over_48_own"] - 0.914) < 2e-3
    )
    fr = {t: per[t]["frac_outside_R_out"] for t in per}
    laws = {}
    for t in per:
        sh = np.array(per[t]["shell_R_meandens_sum"])
        m = (sh[:, 0] >= 12.0) & (sh[:, 0] <= 0.5 * parse_tag(t)[1] - 6.0)
        if m.sum() >= 3:
            b = np.polyfit(np.log(sh[m, 0]), np.log(sh[m, 1]), 1)[0]
            laws[t] = float(b)
    c18 = {t: per[t]["C_within_R"]["18"] for t in per}
    results["4_spin_gate"]["numbers"]["frac_outside_R_out"] = fr
    results["4_spin_gate"]["numbers"]["shell_density_exponent_R12_to_edge"] = laws
    results["4_spin_gate"]["numbers"]["C_within_18_by_row"] = c18
    results["4_spin_gate"]["verdict"] = "QUALIFIED" if ratios_ok else "REFUTED"
    n32 = per["S1_d0.3_w25_n32_L48"]
    n96 = per["S1_d0.3_w25_n64_L96"]
    results["4_spin_gate"]["note"] = (
        "C_rot_z, omega* and J* reproduce to roundoff by an own per-cell density; the generator "
        "time row is zero and the spatial-only kin is identical. Two qualifications. (1) The "
        "generator as coded, Jz M - M Jz^T = Jz M + M Jz, is ANTISYMMETRIC (its symmetric part "
        f"is {n32['a0_code_symmetric_part_max']:.1e}), not the rotation [Jz, M] (symmetric, "
        f"max difference {n32['a0_code_minus_commutator_max']:.2f}); with the commutator C is "
        f"{n32['C_commutator']:.0f} instead of {n32['C_code_own_density']:.0f} on n32 L48, and "
        f"the rigid spin plus orbital generator gives {n32['C_rigid_spin_plus_orbital']:.0f}. "
        "The catalog value is the same generator, envelope-weighted (renv 10) and unit-normalized, "
        "which is why it is 4 orders smaller. (2) The inertia is a halo integral: on n32 L48 "
        f"C(6) = {n32['C_within_R']['6']:.0f}, C(12) = {n32['C_within_R']['12']:.0f}, C(18) = "
        f"{n32['C_within_R']['18']:.0f}, C(24) = {n32['C_within_R']['24']:.0f} of "
        f"{n32['C_code_own_density']:.0f}, fraction outside L/2 - 3h = {fr['S1_d0.3_w25_n32_L48']:.3f}, "
        f"pinned shell {n32['frac_in_pinned_shell']:.3f}; shell density exponents "
        f"{ {k: round(v, 2) for k, v in laws.items()} } (r 12 to the edge; a density falling as "
        "r^-3 to r^-3.8 makes C(R) converge, but slowly: not a surface law, a tail that a 96 box "
        f"has not exhausted); on n64 L96 C(18) = {n96['C_within_R']['18']:.0f} against "
        f"{n96['C_code_own_density']:.0f} total, so 56 percent of the L 96 inertia sits beyond r 18; "
        "the inner part is not box-independent either (C(18) by row: "
        f"{ {k[-7:]: round(v) for k, v in c18.items()} }, a 21 percent change from L 48 to L 72 at "
        "h 1.5 on unconverged rows). The growth with L is the halo's tail plus the unconverged "
        "inner field, not the pinned shell (1 to 3 percent). "
        "SPIN_GATE_BOX_DEPENDENT is therefore a statement about the box and the generator on "
        "unconverged rows, not about the charge; the pre-registered reading 'the carrier is "
        "elsewhere' is the box's halo."
    )


if __name__ == "__main__":
    main()
