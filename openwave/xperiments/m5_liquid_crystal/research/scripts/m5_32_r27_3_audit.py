"""M5.32 R27-3 AUDIT: an independent attempt to refute the four reads of `m5_32_r27_3_reads.py`.

The audited script was NOT read. Every quantity below is recomputed from the stored fields with
this file's own code, on top of the shared stack modules (B3 stencils and coordinates, the R20
certified energy density, the R26-0 partition reader, which is allowed as an already-audited
instrument, and the R26-3 charge density, used only as a cross-check against this file's own
Jacobian).

(a) the tube-masked l = 1 moment of the pair-gap deviation: own least squares (nine real
    harmonics), own mask from the reader's carriers, four floors (the claimed rms / sqrt(N), the
    normal-equation standard error, a permutation null, a bootstrap), the model-extension shift
    (l <= 4), a solid-angle projection on interpolated spheres, and a second shell width.
(b) the charge in the lines: own Jacobian (np.gradient), own tubes, own sphere-degree route
    (solid-angle triangles), the |rho| radial and per-volume distribution, and its concentration
    where the director is ill-defined.
(c) the carriers against the box: the reader at r 9 and 18, and an independent read of where the
    pair gap closes (lattice clusters and interpolated-sphere clusters).
(d) the Coulomb calibration: own sympy for C on M = n n^T, own lattice sums, own E(> R), density
    exponents (9-18, 12-24, 24-33), c(R) with the author's R_max and with the true free radius,
    and the L 72 versus L 96 density comparison at equal r.

Single process, no multiprocessing. Writes data/m5_32_r27_3_audit.json.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy import ndimage
from scipy.ndimage import map_coordinates
from scipy.special import sph_harm_y

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r27_3_audit.json")
CLAIMS_JSON = os.path.join(DATA, "m5_32_r27_3_reads.json")
T0 = time.time()
G = 8.0
W1S = 25.0
DELTA = 0.3
RNG = np.random.default_rng(20260927)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


F0 = _load("m5_32_r26_0_form", "m5_32_r26_0_form.py")
R26_3 = _load("m5_32_r26_3_spin", "m5_32_r26_3_spin.py")
R21, R20, B3, R0, W1 = F0.R21, F0.R20, F0.B3, F0.R0, F0.W1


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def cfg_pot(n, L, delta=DELTA, w1s=W1S):
    cfg = R21.cfg_of(n, L, G, delta)
    pot = ("v4", R0.roots_of(cfg), W1 * w1s)
    return cfg, pot


S1_ROWS = {
    "S1_d0.3_w25_n32_L48": (32, 48.0),
    "S1_d0.3_w25_n48_L72": (48, 72.0),
    "S1_d0.3_w25_n64_L96": (64, 96.0),
    "S1_d0.3_w25_n48_L48": (48, 48.0),
    "S1_d0.3_w25_n64_L64": (64, 64.0),
}
BRANCHES = ["1_1_1_1", "3_1", "2_2", "2_1_1", "4"]


def row_path(name):
    if name.startswith("S1_"):
        return os.path.join(DATA, "m5_32_r25_2", name + ".npz")
    return os.path.join(DATA, "m5_32_r26_4", name + ".npz")


def branch_row(part):
    return f"P{part}_d0.3_w25_n32_L48"


def load_row(name):
    if name.startswith("S1_"):
        n, L = S1_ROWS[name]
    else:
        n, L = 32, 48.0
    M = np.load(row_path(name))["M"].astype(np.float64)
    cfg, pot = cfg_pot(n, L)
    return M, cfg, pot


def geometry(cfg):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    cth = Z / np.maximum(r, 1e-12)
    phi = np.arctan2(Y, X)
    pin = B3.pin_shell(n, h)
    return X, Y, Z, r, cth, phi, pin


def spectral(M, delta=DELTA):
    lam = np.linalg.eigvalsh(M[..., 1:, 1:])
    gap_pair = lam[..., 1] - lam[..., 0]
    gap_dir = lam[..., 2] - lam[..., 1]
    return gap_pair - delta, gap_pair, gap_dir


# ============================================================================
# (a) the tube-masked moment
# ============================================================================
def real_sh(cth, phi, lmax):
    """real harmonics up to lmax, column 0 = 1, column 1 = cos(theta) exactly (the claim's a_1
    convention), the rest any spanning set of the remaining (l, m)."""
    th = np.arccos(np.clip(cth, -1.0, 1.0))
    cols = [np.ones_like(cth), cth]
    for l in range(1, lmax + 1):
        for m in range(-l, l + 1):
            if l == 1 and m == 0:
                continue
            if m == 0:
                cols.append(np.real(sph_harm_y(l, 0, th, phi)))
            elif m > 0:
                cols.append(np.sqrt(2.0) * np.real(sph_harm_y(l, m, th, phi)))
            else:
                cols.append(np.sqrt(2.0) * np.imag(sph_harm_y(l, -m, th, phi)))
    return np.stack(cols, 1)


def lsq_a1(vals, cth, phi, lmax=2):
    A = real_sh(cth, phi, lmax)
    coef, *_ = np.linalg.lstsq(A, vals, rcond=None)
    res = vals - A @ coef
    rms = float(np.sqrt(np.mean(res**2)))
    N = len(vals)
    dof = max(N - A.shape[1], 1)
    cov = np.linalg.pinv(A.T @ A) * (np.sum(res**2) / dof)
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    # the l = 1 equatorial components are columns 2 and 3 in the ordering (l=1, m=-1) then (1, +1)
    # with this file's normalization; convert them to the plain sin(theta) sin(phi) / cos(phi)
    # coefficients: sqrt(2) Re/Im Y_1^1 = -sqrt(3 / 8pi) sqrt(2) sin(theta) {cos, sin}(phi)
    k = -np.sqrt(3.0 / (8.0 * np.pi)) * np.sqrt(2.0)
    return {
        "a0": float(coef[0]),
        "a1": float(coef[1]),
        "c1_y": float(coef[2] * k),
        "b1_x": float(coef[3] * k),
        "l1_norm": float(np.sqrt(coef[1] ** 2 + (coef[2] * k) ** 2 + (coef[3] * k) ** 2)),
        "m0_content": float(
            coef[1] ** 2 / max(coef[1] ** 2 + (coef[2] * k) ** 2 + (coef[3] * k) ** 2, 1e-300)
        ),
        "rms": rms,
        "N": int(N),
        "floor_claimed": float(rms / np.sqrt(N)),
        "se_a1": float(se[1]),
        "se_b1_x": float(se[3] * abs(k)),
        "se_c1_y": float(se[2] * abs(k)),
    }


def perm_boot_floor(vals, cth, phi, nrep=200):
    A = real_sh(cth, phi, 2)
    N = len(vals)
    a1p, a1b = [], []
    for _ in range(nrep):
        v = RNG.permutation(vals)
        c, *_ = np.linalg.lstsq(A, v, rcond=None)
        a1p.append(c[1])
        idx = RNG.integers(0, N, N)
        c, *_ = np.linalg.lstsq(A[idx], vals[idx], rcond=None)
        a1b.append(c[1])
    return float(np.std(a1p)), float(np.std(a1b))


def carriers_of(M, cfg, radii):
    """the reader's carriers per sphere as cartesian points (theta finite only)."""
    out = {}
    for R in radii:
        if R > 0.5 * cfg["L"] - 1.6:
            continue
        rd = F0.partition_reader(M, cfg, R, DELTA)
        pts = []
        for c in rd["carriers"]:
            th = c["theta"]
            if th is None or not np.isfinite(th):
                continue
            ph = c["phi"]
            pts.append(
                (
                    R * np.sin(th) * np.cos(ph),
                    R * np.sin(th) * np.sin(ph),
                    R * np.cos(th),
                    c["half_units"],
                    c["rho_from_axis"],
                    c["z"],
                )
            )
        out[R] = {"partition": rd["partition"], "pts": pts}
    return out


def build_mask(X, Y, Z, r, carr, rad, r_core=6.0):
    """the union of balls of radius `rad` around every carrier point on every sphere, plus tubes
    of radius `rad` inside r < r_core along the half-lines through the r 9 carriers."""
    mask = np.zeros_like(r, dtype=bool)
    for R, d in carr.items():
        for p in d["pts"]:
            mask |= (X - p[0]) ** 2 + (Y - p[1]) ** 2 + (Z - p[2]) ** 2 < rad**2
    if 9.0 in carr:
        for p in carr[9.0]["pts"]:
            u = np.array(p[:3]) / np.linalg.norm(p[:3])
            proj = X * u[0] + Y * u[1] + Z * u[2]
            perp2 = r * r - proj * proj
            mask |= (proj > 0) & (perp2 < rad**2) & (r < r_core)
    return mask


def project_a1_sphere(dev, cfg, R, nth=90, nph=180):
    """a_1 = 3 / (4 pi) * integral dev cos(theta) dOmega on the trilinearly interpolated sphere:
    the proper P_1 projection, no least squares, no lattice-shell weighting."""
    n, h = cfg["n"], cfg["h"]
    th = (np.arange(nth) + 0.5) * np.pi / nth
    ph = np.arange(nph) * 2 * np.pi / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    pts = R * np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
    coords = np.stack([(pts[..., a] / h + (n - 1) / 2.0).ravel() for a in range(3)])
    v = map_coordinates(dev, coords, order=1, mode="nearest").reshape(TH.shape)
    dO = np.sin(TH) * (np.pi / nth) * (2 * np.pi / nph)
    a1 = 3.0 / (4 * np.pi) * np.sum(v * np.cos(TH) * dO)
    b1 = 3.0 / (4 * np.pi) * np.sum(v * np.sin(TH) * np.cos(PH) * dO)
    c1 = 3.0 / (4 * np.pi) * np.sum(v * np.sin(TH) * np.sin(PH) * dO)
    a0 = np.sum(v * dO) / (4 * np.pi)
    return float(a1), float(b1), float(c1), float(a0)


def slope_loglog(rs, vals, lo, hi):
    rs, vals = np.asarray(rs), np.asarray(vals)
    m = (rs >= lo) & (rs <= hi) & (np.abs(vals) > 0)
    if m.sum() < 3:
        return None
    return float(np.polyfit(np.log(rs[m]), np.log(np.abs(vals[m])), 1)[0])


def audit_a(part, claims):
    name = branch_row(part)
    M, cfg, pot = load_row(name)
    X, Y, Z, r, cth, phi, pin = geometry(cfg)
    dev, gap_pair, gap_dir = spectral(M)
    radii = list(np.arange(3.0, 21.0 + 1e-9, 1.5))
    carr8 = carriers_of(M, cfg, [3.0, 4.5, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0])
    mask8 = build_mask(X, Y, Z, r, carr8, 3.0 * cfg["h"])
    carr13 = carriers_of(M, cfg, radii)
    mask13 = build_mask(X, Y, Z, r, carr13, 3.0 * cfg["h"])
    out = {"row": name, "carriers_r9": carr8[9.0], "shells": []}
    for Rs in radii:
        sh = {"r": float(Rs)}
        for tag, msk, width in (
            ("unmasked_w1.5", None, 1.5),
            ("unmasked_w1.0", None, 1.0),
            ("balls8_w1.5", mask8, 1.5),
            ("balls13_w1.5", mask13, 1.5),
        ):
            m = (np.abs(r - Rs) < 0.5 * width) & (~pin)
            if msk is not None:
                m &= ~msk
            if m.sum() < 30:
                sh[tag] = {"N": int(m.sum()), "skipped": True}
                continue
            f = lsq_a1(dev[m], cth[m], phi[m], 2)
            pf, bf = perm_boot_floor(dev[m], cth[m], phi[m])
            f["floor_perm"] = pf
            f["floor_boot"] = bf
            if m.sum() > 120:
                f4 = lsq_a1(dev[m], cth[m], phi[m], 4)
                f["a1_lmax4"] = f4["a1"]
                f["shift_l2_to_l4"] = float(f4["a1"] - f["a1"])
            if m.sum() > 300:
                f6 = lsq_a1(dev[m], cth[m], phi[m], 6)
                f["a1_lmax6"] = f6["a1"]
            sh[tag] = f
        a1p, b1p, c1p, a0p = project_a1_sphere(dev, cfg, Rs)
        sh["projection_sphere"] = {"a1": a1p, "b1_x": b1p, "c1_y": c1p, "a0": a0p}
        a1m, _, _, _ = project_a1_sphere(dev, cfg, Rs - 0.5)
        a1q, _, _, _ = project_a1_sphere(dev, cfg, Rs + 0.5)
        sh["projection_sphere_pm0.5"] = [a1m, a1q]
        out["shells"].append(sh)
    rs = [s["r"] for s in out["shells"]]
    for tag in ("unmasked_w1.5", "unmasked_w1.0", "balls8_w1.5", "balls13_w1.5"):
        vals = [s[tag].get("a1", 0.0) for s in out["shells"]]
        out[f"slope_a1_{tag}_6_18"] = slope_loglog(rs, vals, 6.0, 18.0)
        out[f"signs_a1_{tag}"] = [int(np.sign(v)) for v in vals]
    out["slope_a1_projection_6_18"] = slope_loglog(
        rs, [s["projection_sphere"]["a1"] for s in out["shells"]], 6.0, 18.0
    )
    # the pre-registered first condition on the 3h balls mask: masked |a1| < 3 floor on all
    # shells 6..18, with each floor
    win = [s for s in out["shells"] if 6.0 <= s["r"] <= 18.0 and "a1" in s["balls8_w1.5"]]
    for fl in ("floor_claimed", "se_a1", "floor_perm", "floor_boot"):
        out[f"masked_below_3x_{fl}_all_shells_6_18"] = bool(
            all(abs(s["balls8_w1.5"]["a1"]) < 3 * s["balls8_w1.5"][fl] for s in win)
        )
    out["equatorial_l1_over_se_unmasked_6_18"] = [
        float(
            np.hypot(s["unmasked_w1.5"]["b1_x"], s["unmasked_w1.5"]["c1_y"])
            / max(np.hypot(s["unmasked_w1.5"]["se_b1_x"], s["unmasked_w1.5"]["se_c1_y"]), 1e-300)
        )
        for s in out["shells"]
        if 6.0 <= s["r"] <= 18.0
    ]
    # the claimed numbers for the same shells
    cl = claims["a_moment"][part]["tube_3h"]
    out["claimed_unmasked_a1"] = {str(s["r"]): s["a1"] for s in cl["unmasked"]["shells"]}
    out["claimed_balls_a1"] = {str(s["r"]): s["a1"] for s in cl["balls"]["shells"]}
    out["claimed_floor"] = {
        str(s["r"]): s["residual_rms"] / np.sqrt(s["cells"]) for s in cl["balls"]["shells"]
    }
    out["claimed_slope_unmasked"] = cl["unmasked"]["slope_a1"]
    return out


# ============================================================================
# (b) the charge in the lines
# ============================================================================
def director_outward(M, X, Y, Z, r):
    lam, V = np.linalg.eigh(M[..., 1:, 1:])
    d = V[..., :, 2]
    rhat = np.stack([X, Y, Z], -1) / np.maximum(r, 1e-12)[..., None]
    al = np.einsum("...a,...a->...", d, rhat)
    sgn = np.where(al < 0, -1.0, 1.0)
    return d * sgn[..., None], np.abs(al), lam[..., 2] - lam[..., 1]


def own_rho(d, h):
    """rho = div J, J_k = (1 / 8 pi) eps_kij d . (d_i d x d_j d), central differences (np.gradient),
    h^3-weighted so that the sum is the charge."""
    dd = [np.gradient(d, h, axis=a) for a in range(3)]
    J = [
        np.einsum("...a,...a->...", d, np.cross(dd[1], dd[2])) / (4 * np.pi),
        np.einsum("...a,...a->...", d, np.cross(dd[2], dd[0])) / (4 * np.pi),
        np.einsum("...a,...a->...", d, np.cross(dd[0], dd[1])) / (4 * np.pi),
    ]
    rho = sum(np.gradient(J[k], h, axis=k) for k in range(3))
    return rho * h**3


def sphere_degree(M, cfg, R, nth=96, nph=192):
    """the outward-oriented director's degree on the sphere R by signed spherical-triangle areas
    of the trilinearly interpolated block (own implementation of the solid-angle route)."""
    n, h = cfg["n"], cfg["h"]
    th = (np.arange(nth) + 0.5) * np.pi / nth
    ph = np.arange(nph) * 2 * np.pi / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    rhat = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
    pts = R * rhat
    coords = np.stack([(pts[..., a] / h + (n - 1) / 2.0).ravel() for a in range(3)])
    S = np.zeros(TH.shape + (3, 3))
    for a in range(3):
        for b in range(a, 3):
            v = map_coordinates(M[..., 1 + a, 1 + b], coords, order=1, mode="nearest")
            S[..., a, b] = v.reshape(TH.shape)
            S[..., b, a] = S[..., a, b]
    lam, V = np.linalg.eigh(S)
    d = V[..., :, 2]
    al = np.einsum("...a,...a->...", d, rhat)
    d = d * np.where(al < 0, -1.0, 1.0)[..., None]

    def tri(a, b, c):
        num = np.einsum("...a,...a->...", a, np.cross(b, c))
        den = (
            1.0
            + np.einsum("...a,...a->...", a, b)
            + np.einsum("...a,...a->...", b, c)
            + np.einsum("...a,...a->...", c, a)
        )
        return 2.0 * np.arctan2(num, den)

    # grid quads split in two triangles; orientation: (theta, phi) frame is clockwise from
    # outside, so the outward degree is minus the sum
    a, b = d[:-1, :], np.roll(d, -1, axis=1)[:-1, :]
    c, e = np.roll(d, -1, axis=1)[1:, :], d[1:, :]
    area = np.sum(tri(a, b, c)) + np.sum(tri(a, c, e))
    north = d[0].mean(axis=0)
    north /= np.linalg.norm(north)
    south = d[-1].mean(axis=0)
    south /= np.linalg.norm(south)
    r0 = d[0]
    area += np.sum(tri(np.broadcast_to(north, r0.shape), r0, np.roll(r0, -1, axis=0)))
    r1 = d[-1]
    area += np.sum(tri(r1, np.broadcast_to(south, r1.shape), np.roll(r1, -1, axis=0)))
    return {
        "R": float(R),
        "degree": float(-area / (4 * np.pi)),
        "min_abs_align": float(np.abs(al).min()),
        "min_gap_dir": float((lam[..., 2] - lam[..., 1]).min()),
        "frac_align_below_0.3": float(np.mean(np.abs(al) < 0.3)),
    }


def lift_3d(d):
    """sign-propagate the director into a continuous vector field (a line field lifts iff no
    disclination line threads the region): along axis 2 on the first line, axis 1 on the first
    plane, axis 0 through the volume. Returns the lifted field, the frustrated bonds (neighbor
    pairs still antiparallel after the lift) and the bond count."""
    n = d.shape[0]
    s = np.ones(d.shape[:3])
    dots = np.einsum("...a,...a->...", d[0, 0, :-1], d[0, 0, 1:])
    s[0, 0, 1:] = np.cumprod(np.where(dots < 0, -1.0, 1.0))
    dots = np.einsum("...a,...a->...", d[0, :-1, :], d[0, 1:, :])
    s[0, 1:, :] = s[0, 0:1, :] * np.cumprod(np.where(dots < 0, -1.0, 1.0), axis=0)
    dots = np.einsum("...a,...a->...", d[:-1], d[1:])
    s[1:] = s[0:1] * np.cumprod(np.where(dots < 0, -1.0, 1.0), axis=0)
    dl = d * s[..., None]
    fr = np.zeros(d.shape[:3], dtype=bool)
    for ax in range(3):
        a = np.take(dl, range(0, n - 1), axis=ax)
        b = np.take(dl, range(1, n), axis=ax)
        bad = np.einsum("...a,...a->...", a, b) < 0
        sl = [slice(None)] * 3
        sl[ax] = slice(0, n - 1)
        fr[tuple(sl)] |= bad
        sl[ax] = slice(1, n)
        fr[tuple(sl)] |= bad
    return dl, fr, int(np.sum(fr)), 3 * n * n * (n - 1)


def sphere_degree_lift(M, cfg, R, nth=96, nph=192):
    """the degree of the director on the sphere R under a CONTINUOUS lift (sign propagation along
    phi on the first ring, then ring to ring), the frustrated edges counted; |degree| is the
    lift-independent quantity (the global sign is a convention)."""
    n, h = cfg["n"], cfg["h"]
    th = (np.arange(nth) + 0.5) * np.pi / nth
    ph = np.arange(nph) * 2 * np.pi / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    rhat = np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
    pts = R * rhat
    coords = np.stack([(pts[..., a] / h + (n - 1) / 2.0).ravel() for a in range(3)])
    S = np.zeros(TH.shape + (3, 3))
    for a in range(3):
        for b in range(a, 3):
            v = map_coordinates(M[..., 1 + a, 1 + b], coords, order=1, mode="nearest")
            S[..., a, b] = v.reshape(TH.shape)
            S[..., b, a] = S[..., a, b]
    lam, V = np.linalg.eigh(S)
    d = V[..., :, 2]
    s = np.ones(TH.shape)
    dots = np.einsum("...a,...a->...", d[0, :-1], d[0, 1:])
    s[0, 1:] = np.cumprod(np.where(dots < 0, -1.0, 1.0))
    dots = np.einsum("...a,...a->...", d[:-1, :], d[1:, :])
    s[1:, :] = s[0:1, :] * np.cumprod(np.where(dots < 0, -1.0, 1.0), axis=0)
    d = d * s[..., None]
    fr = int(np.sum(np.einsum("...a,...a->...", d, np.roll(d, -1, axis=1)) < 0)) + int(
        np.sum(np.einsum("...a,...a->...", d[:-1], d[1:]) < 0)
    )

    def tri(a, b, c):
        num = np.einsum("...a,...a->...", a, np.cross(b, c))
        den = (
            1.0
            + np.einsum("...a,...a->...", a, b)
            + np.einsum("...a,...a->...", b, c)
            + np.einsum("...a,...a->...", c, a)
        )
        return 2.0 * np.arctan2(num, den)

    a, b = d[:-1, :], np.roll(d, -1, axis=1)[:-1, :]
    c, e = np.roll(d, -1, axis=1)[1:, :], d[1:, :]
    area = np.sum(tri(a, b, c)) + np.sum(tri(a, c, e))
    north = d[0].mean(axis=0)
    north /= max(np.linalg.norm(north), 1e-300)
    south = d[-1].mean(axis=0)
    south /= max(np.linalg.norm(south), 1e-300)
    r0 = d[0]
    area += np.sum(tri(np.broadcast_to(north, r0.shape), r0, np.roll(r0, -1, axis=0)))
    r1 = d[-1]
    area += np.sum(tri(r1, np.broadcast_to(south, r1.shape), np.roll(r1, -1, axis=0)))
    return {
        "R": float(R),
        "abs_degree_lift": float(abs(area / (4 * np.pi))),
        "frustrated_edges": fr,
        "min_gap_dir": float((lam[..., 2] - lam[..., 1]).min()),
    }


def tube_mask(X, Y, Z, r, pts9, r_t, r_core=6.0):
    """cells outside the core within r_t of the half-lines from the origin through the r 9
    carriers (own construction)."""
    tubes = np.zeros_like(r, dtype=bool)
    for p in pts9:
        u = np.array(p[:3]) / np.linalg.norm(p[:3])
        proj = X * u[0] + Y * u[1] + Z * u[2]
        perp2 = np.maximum(r * r - proj * proj, 0.0)
        tubes |= (proj > 0) & (perp2 < r_t**2)
    return tubes & (r >= r_core)


def audit_b(name, claims):
    M, cfg, pot = load_row(name)
    n, h = cfg["n"], cfg["h"]
    X, Y, Z, r, cth, phi, pin = geometry(cfg)
    d, al, gap_dir = director_outward(M, X, Y, Z, r)
    rho = own_rho(d, h)
    rho_ref, _ = R26_3.charge_density(M, cfg)
    free = ~pin
    carr9 = carriers_of(M, cfg, [9.0])[9.0]
    r_t = 3.0 * h
    tubes = tube_mask(X, Y, Z, r, carr9["pts"], r_t)
    core = r < 6.0
    A = np.abs(rho)
    tot_abs = float(np.sum(A[free]))
    out = {
        "row": name,
        "n_L_h": [n, cfg["L"], h],
        "r_t": r_t,
        "n_lines": len(carr9["pts"]),
        "net_free": float(np.sum(rho[free])),
        "abs_free": tot_abs,
        "abs_free_R26_3": float(np.sum(np.abs(rho_ref[free]))),
        "net_outer": float(np.sum(rho[free & ~core])),
        "net_outer_in_tubes": float(np.sum(rho[free & tubes])),
        "frac_abs_core": float(np.sum(A[free & core]) / tot_abs),
        "frac_abs_tubes": float(np.sum(A[free & tubes]) / tot_abs),
        "frac_abs_rest": float(np.sum(A[free & ~core & ~tubes]) / tot_abs),
        "r2_abs_weighted": float(np.sum(A[free] * r[free] ** 2) / tot_abs),
        "corr_own_vs_R26_3": float(np.corrcoef(rho[free].ravel(), rho_ref[free].ravel())[0, 1]),
    }
    # where does |rho| sit: per shell, its sum, its per-volume mean, the fraction where the
    # director is ill-defined (gap_dir < 0.1 or |align| < 0.3)
    shells = []
    for Rs in np.arange(1.5, 0.5 * cfg["L"], 3.0):
        m = free & (np.abs(r - Rs) < 1.5)
        if m.sum() == 0:
            continue
        shells.append(
            {
                "r": float(Rs),
                "abs_sum": float(np.sum(A[m])),
                "net_sum": float(np.sum(rho[m])),
                "abs_per_cell": float(np.mean(A[m])),
                "cells": int(m.sum()),
            }
        )
    out["shells"] = shells
    ill = (gap_dir < 0.1) | (al < 0.3)
    out["frac_abs_in_ill_defined_director_cells"] = float(np.sum(A[free & ill]) / tot_abs)
    out["frac_cells_ill_defined"] = float(np.mean(ill[free]))
    out["frac_abs_in_cells_align_below_0.5"] = float(np.sum(A[free & (al < 0.5)]) / tot_abs)
    # the continuous-lift route: a lifted vector field (frustration counted), its own Jacobian
    d_un = d * np.where(d[..., 2] < 0, -1.0, 1.0)[..., None]  # any fixed convention first
    dl, fr_mask, n_fr, n_bonds = lift_3d(d_un)
    rho_l = own_rho(dl, h)
    Al = np.abs(rho_l)
    tot_l = float(np.sum(Al[free]))
    near_fr = ndimage.binary_dilation(fr_mask, iterations=2)
    out["lift3d"] = {
        "frustrated_cells": n_fr,
        "frustrated_fraction_of_cells": float(n_fr / dl.shape[0] ** 3),
        "frustrated_cells_in_free": int(np.sum(fr_mask & free)),
        "min_gap_dir_at_frustrated": float(gap_dir[fr_mask].min()) if n_fr else None,
        "net_free": float(np.sum(rho_l[free])),
        "abs_free": tot_l,
        "net_inside_r6": float(np.sum(rho_l[free & core])),
        "frac_abs_core": float(np.sum(Al[free & core]) / tot_l),
        "frac_abs_tubes": float(np.sum(Al[free & tubes]) / tot_l),
        "frac_abs_rest": float(np.sum(Al[free & ~core & ~tubes]) / tot_l),
        "frac_abs_within_2_cells_of_frustration": float(np.sum(Al[free & near_fr]) / tot_l),
        "r2_abs_weighted": float(np.sum(Al[free] * r[free] ** 2) / tot_l),
        "abs_per_cell_by_shell": [
            [float(Rs), float(np.mean(Al[free & (np.abs(r - Rs) < 1.5)]))]
            for Rs in np.arange(1.5, 0.5 * cfg["L"], 3.0)
            if np.any(free & (np.abs(r - Rs) < 1.5))
        ],
        "net_inside_R": [
            [float(R), float(np.sum(rho_l[free & (r < R)]))]
            for R in (3.0, 6.0, 9.0, 12.0, 15.0, 18.0)
        ],
    }
    out["sphere_degree_lift"] = [
        sphere_degree_lift(M, cfg, R)
        for R in (3.0, 4.5, 6.0, 7.5, 9.0, 10.5, 12.0, 15.0, 18.0, 21.0)
        if R <= 0.5 * cfg["L"] - 2 * h
    ]
    # the sphere-degree route
    degs = []
    for R in (3.0, 4.5, 6.0, 7.5, 9.0, 10.5, 12.0, 15.0, 18.0, 21.0):
        if R > 0.5 * cfg["L"] - 2 * h:
            continue
        degs.append(sphere_degree(M, cfg, R))
    out["sphere_degree"] = degs
    # the Jacobian net charge inside the same spheres (own rho), for the cross-check
    out["jacobian_net_inside"] = [
        {"R": dg["R"], "net": float(np.sum(rho[free & (r < dg["R"])]))} for dg in degs
    ]
    cl = claims["b_charge_split"]["rows"][name]
    out["claimed"] = {
        k: cl[k]
        for k in (
            "fraction_inside_core",
            "fraction_outer_in_tubes",
            "fraction_outer_outside_tubes",
            "net_charge_outer",
            "charge_abs_total_free",
            "r2_abs_weighted_all",
        )
    }
    return out


def audit_b_control():
    """own hedgehog control: M = n n^T (delta 0) on n 32 L 48; the charge sits inside r 6."""
    cfg, pot = cfg_pot(32, 48.0, delta=0.0)
    X, Y, Z, r, cth, phi, pin = geometry(cfg)
    M = hedgehog_field(cfg)
    d, al, gap_dir = director_outward(M, X, Y, Z, r)
    rho = own_rho(d, cfg["h"])
    free = ~pin
    A = np.abs(rho)
    return {
        "net_free": float(np.sum(rho[free])),
        "abs_free": float(np.sum(A[free])),
        "net_inside_r6": float(np.sum(rho[free & (r < 6)])),
        "abs_outside_r6": float(np.sum(A[free & (r >= 6)])),
        "sphere_degree_r6": sphere_degree(M, cfg, 6.0)["degree"],
    }


# ============================================================================
# (c) the carriers against the box
# ============================================================================
def closure_clusters_lattice(gap_pair, X, Y, Z, r, R=9.0, width=1.5, thr=0.1 * DELTA):
    m = (np.abs(r - R) < 0.5 * width) & (gap_pair < thr)
    lab, k = ndimage.label(m, structure=np.ones((3, 3, 3)))
    out = []
    for i in range(1, k + 1):
        sel = lab == i
        x, y, z = X[sel].mean(), Y[sel].mean(), Z[sel].mean()
        out.append(
            {
                "cells": int(sel.sum()),
                "rho": float(np.hypot(x, y)),
                "z": float(z),
                "phi_deg": float(np.degrees(np.arctan2(y, x)) % 360),
            }
        )
    out.sort(key=lambda c: -c["cells"])
    return out


def closure_clusters_sphere(M, cfg, R=9.0, thr=0.1 * DELTA, nth=128, nph=256):
    """the pair gap on the interpolated sphere; connected regions below thr (phi wraps)."""
    n, h = cfg["n"], cfg["h"]
    th = (np.arange(nth) + 0.5) * np.pi / nth
    ph = np.arange(nph) * 2 * np.pi / nph
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    pts = R * np.stack([np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)], -1)
    coords = np.stack([(pts[..., a] / h + (n - 1) / 2.0).ravel() for a in range(3)])
    S = np.zeros(TH.shape + (3, 3))
    for a in range(3):
        for b in range(a, 3):
            v = map_coordinates(M[..., 1 + a, 1 + b], coords, order=1, mode="nearest")
            S[..., a, b] = v.reshape(TH.shape)
            S[..., b, a] = S[..., a, b]
    lam = np.linalg.eigvalsh(S)
    gp = lam[..., 1] - lam[..., 0]
    m = gp < thr
    # wrap phi by tiling three times, labeling, and keeping the middle copy's labels
    m3 = np.concatenate([m, m, m], axis=1)
    lab3, k = ndimage.label(m3, structure=np.ones((3, 3)))
    lab = lab3[:, nph : 2 * nph]
    out = []
    for i in np.unique(lab):
        if i == 0:
            continue
        sel = lab == i
        w = np.sin(TH[sel])
        x = np.sum(w * pts[..., 0][sel]) / w.sum()
        y = np.sum(w * pts[..., 1][sel]) / w.sum()
        z = np.sum(w * pts[..., 2][sel]) / w.sum()
        out.append(
            {
                "solid_angle_frac": float(
                    np.sum(w) * (np.pi / nth) * (2 * np.pi / nph) / (4 * np.pi)
                ),
                "rho": float(np.hypot(x, y) * R / max(np.sqrt(x * x + y * y + z * z), 1e-12)),
                "z": float(z * R / max(np.sqrt(x * x + y * y + z * z), 1e-12)),
                "phi_deg": float(np.degrees(np.arctan2(y, x)) % 360),
                "min_gap_pair": float(gp[sel].min()),
            }
        )
    out.sort(key=lambda c: -c["solid_angle_frac"])
    return out, float(gp.min())


def audit_c(name):
    M, cfg, pot = load_row(name)
    X, Y, Z, r, cth, phi, pin = geometry(cfg)
    dev, gap_pair, gap_dir = spectral(M)
    out = {"row": name, "n_L_h": [cfg["n"], cfg["L"], cfg["h"]]}
    for R in (9.0, 18.0):
        rd = F0.partition_reader(M, cfg, R, DELTA)
        out[f"reader_R{R:g}"] = {
            "partition": rd["partition"],
            "n_carriers": len(rd["carriers"]),
            "rho_hu_z": [[c["rho_from_axis"], c["half_units"], c["z"]] for c in rd["carriers"]],
            "unreadable_plaquettes": rd.get("unreadable_plaquettes"),
            "flagged": rd.get("cap_ring_flagged"),
        }
    out["closure_lattice_r9"] = closure_clusters_lattice(gap_pair, X, Y, Z, r)
    cs, gmin = closure_clusters_sphere(M, cfg, 9.0)
    out["closure_sphere_r9"] = cs
    out["min_gap_pair_on_sphere_r9"] = gmin
    cs18, gmin18 = closure_clusters_sphere(M, cfg, 18.0)
    out["closure_sphere_r18_count"] = len(cs18)
    out["closure_sphere_r18_top4"] = cs18[:4]
    out["min_gap_pair_on_sphere_r18"] = gmin18
    return out


# ============================================================================
# (d) the Coulomb calibration
# ============================================================================
def sympy_C():
    import sympy as sp

    x, y, z = sp.symbols("x y z", real=True, positive=True)
    rr = sp.sqrt(x * x + y * y + z * z)
    nv = sp.Matrix([x, y, z]) / rr
    Mm = nv * nv.T
    A = [Mm.diff(v) for v in (x, y, z)]
    tot = 0
    for i in range(3):
        for j in range(i + 1, 3):
            F = A[i] * A[j] - A[j] * A[i]
            tot += sum(f**2 for f in F)
    e_r4 = sp.simplify(4 * tot * rr**4)
    num = float((4 * tot * rr**4).subs({x: 0.3, y: -1.1, z: 0.7}).evalf())
    return {"C_symbolic": str(e_r4), "C_numeric_check": num}


def hedgehog_field(cfg):
    X, Y, Z, r, cth, phi, pin = geometry(cfg)
    nh = np.stack([X, Y, Z], -1) / np.maximum(r, 1e-12)[..., None]
    M = np.zeros(r.shape + (4, 4))
    M[..., 0, 0] = -cfg["sg"]
    M[..., 1:, 1:] = nh[..., :, None] * nh[..., None, :]
    return M


def energy_tail(e, r, pin, cfg, Rs):
    """E(> R) on all free cells with r > R (the author's convention) and on the spherical window
    R < r < R_free (this file's), with the true free radius R_free = L/2 - wc h."""
    wc = max(1, int(np.ceil(1.6 / cfg["h"])))
    R_free = 0.5 * cfg["L"] - wc * cfg["h"]
    R_max_author = 0.5 * cfg["L"] - 1.6
    free = ~pin
    rows = []
    for R in Rs:
        E_all = float(np.sum(e[free & (r > R)]))
        E_sph = float(np.sum(e[free & (r > R) & (r < R_free)]))
        rows.append(
            {
                "R": float(R),
                "E_gt_R_all": E_all,
                "E_gt_R_sphere": E_sph,
                "c_author": (
                    float(E_all / (1.0 / R - 1.0 / R_max_author)) if R < R_max_author else None
                ),
                "c_sphere": float(E_sph / (1.0 / R - 1.0 / R_free)) if R < R_free else None,
            }
        )
    return rows, R_free, R_max_author


def shell_density(e, r, pin, cfg, Rs, width=1.5):
    h3 = cfg["h"] ** 3
    out = []
    for R in Rs:
        m = (~pin) & (np.abs(r - R) < 0.5 * width)
        if m.sum() < 10:
            continue
        out.append(
            [
                float(R),
                float(np.mean(e[m]) / h3),
                int(m.sum()),
                float(np.mean(e[m] * r[m] ** 4) / h3),
            ]
        )
    return out


def local_exponents(sd):
    rs = np.array([s[0] for s in sd])
    es = np.array([s[1] for s in sd])
    return [
        [
            float(rs[i]),
            float(rs[i + 1]),
            float(np.log(es[i + 1] / es[i]) / np.log(rs[i + 1] / rs[i])),
        ]
        for i in range(len(rs) - 1)
    ]


def audit_d_hedgehog():
    cfg, pot = cfg_pot(48, 72.0, delta=0.0)
    X, Y, Z, r, cth, phi, pin = geometry(cfg)
    M = hedgehog_field(cfg)
    e = R20.density(M, cfg, pot)
    Rs = list(np.arange(6.0, 30.0 + 1e-9, 1.5))
    tail, R_free, R_ma = energy_tail(e, r, pin, cfg, Rs)
    sd = shell_density(e, r, pin, cfg, Rs)
    # own fit E(> R) = c / R + b over 9..18, c = slope on 1/R
    sel = [t for t in tail if 9.0 <= t["R"] <= 18.0]
    xs = np.array([1.0 / t["R"] for t in sel])
    ys = np.array([t["E_gt_R_all"] for t in sel])
    c, b = np.polyfit(xs, ys, 1)
    ys2 = np.array([t["E_gt_R_sphere"] for t in sel])
    c2, b2 = np.polyfit(xs, ys2, 1)
    # the potential must vanish on the hedgehog (roots sg, 1, 0, 0)
    e_curv_only = R20.density(M, cfg, ("v4", pot[1], 0.0))
    return {
        "cfg": [48, 72.0, cfg["h"]],
        "E_total_free": float(np.sum(e[~pin])),
        "V4_share": float(1.0 - np.sum(e_curv_only[~pin]) / np.sum(e[~pin])),
        "c_fit_all_9_18": float(c),
        "R_eff_all": float(-c / b) if b != 0 else None,
        "c_fit_sphere_9_18": float(c2),
        "R_eff_sphere": float(-c2 / b2) if b2 != 0 else None,
        "c_theory_32pi": float(32 * np.pi),
        "C_lattice_per_shell": [[s[0], s[1] * s[0] ** 4] for s in sd],
        "C_lattice_per_shell_cellwise_e_r4": [[s[0], s[3]] for s in sd],
        "local_exponents": local_exponents(sd),
        "tail": tail,
        "R_free": R_free,
        "R_max_author": R_ma,
    }


def audit_d_row(name, claims):
    M, cfg, pot = load_row(name)
    X, Y, Z, r, cth, phi, pin = geometry(cfg)
    e = R20.density(M, cfg, pot)
    Rs = list(np.arange(6.0, 0.5 * cfg["L"] + 1e-9, 1.5))
    tail, R_free, R_ma = energy_tail(e, r, pin, cfg, Rs)
    sd = shell_density(e, r, pin, cfg, Rs)
    rs = [s[0] for s in sd]
    es = [s[1] for s in sd]
    out = {
        "row": name,
        "n_L_h": [cfg["n"], cfg["L"], cfg["h"]],
        "E_total_free": float(np.sum(e[~pin])),
        "R_free": R_free,
        "R_max_author": R_ma,
        "tail": tail,
        "shell_density": sd,
        "local_exponents": local_exponents(sd),
        "exp_9_18": slope_loglog(rs, es, 9.0, 18.0),
        "exp_12_24": slope_loglog(rs, es, 12.0, 24.0),
        "exp_13.5_24": slope_loglog(rs, es, 13.5, 24.0),
        "exp_24_33": slope_loglog(rs, es, 24.0, 33.0),
        "r_of_density_peak": float(rs[int(np.argmax(es))]),
    }
    ca = {t["R"]: t["c_author"] for t in tail if t["c_author"]}
    cs = {t["R"]: t["c_sphere"] for t in tail if t["c_sphere"]}
    out["c_author_9_18_drift"] = float(ca[18.0] / ca[9.0]) if 18.0 in ca else None
    out["c_sphere_9_18_drift"] = float(cs[18.0] / cs[9.0]) if 18.0 in cs else None
    if 30.0 in ca:
        out["c_author_24_30"] = [ca.get(R) for R in (24.0, 25.5, 27.0, 28.5, 30.0)]
        out["c_sphere_24_30"] = [cs.get(R) for R in (24.0, 25.5, 27.0, 28.5, 30.0)]
    cl = claims["d_coulomb"]["biaxial_rows"][name]
    out["claimed"] = {
        "shell_density_exponent_on_9_18": cl["shell_density_exponent_on_9_18"],
        "c_drift_9_to_18": cl["c_drift_9_to_18"],
        "E_total_free": cl["E_total_free"],
    }
    return out


# ============================================================================
def main():
    with open(CLAIMS_JSON) as f:
        claims = json.load(f)
    J = {"task": "M5.32 R27-3 audit (independent)", "date": "2026-09-27"}

    log("(d) sympy")
    J["d_sympy"] = sympy_C()
    log(f"  C r^4 = {J['d_sympy']}")
    log("(d) hedgehog lattice n48 L72")
    J["d_hedgehog"] = audit_d_hedgehog()
    hh = J["d_hedgehog"]
    log(
        f"  c_fit_all {hh['c_fit_all_9_18']:.2f} (R_eff {hh['R_eff_all']:.1f}), c_fit_sphere "
        f"{hh['c_fit_sphere_9_18']:.2f}, V4 share {hh['V4_share']:.2e}, "
        f"C per shell {[round(c[1], 3) for c in hh['C_lattice_per_shell'][:4]]}"
    )
    J["d_rows"] = {}
    for name in S1_ROWS:
        log(f"(d) row {name}")
        J["d_rows"][name] = audit_d_row(name, claims)
        o = J["d_rows"][name]
        log(
            f"  E {o['E_total_free']:.4f} exp9-18 {o['exp_9_18']:.2f} exp12-24 "
            f"{o['exp_12_24']} exp24-33 {o['exp_24_33']} peak r {o['r_of_density_peak']} "
            f"c drift author {o['c_author_9_18_drift']:.2f} sphere {o['c_sphere_9_18_drift']:.2f}"
        )

    J["c_rows"] = {}
    for name in list(S1_ROWS) + [branch_row("1_1_1_1")]:
        log(f"(c) row {name}")
        J["c_rows"][name] = audit_c(name)
        o = J["c_rows"][name]
        log(
            f"  reader R9 {o['reader_R9']['partition']} rho {[round(c[0], 2) if c[0] is not None else None for c in o['reader_R9']['rho_hu_z']]}"
            f" | R18 n {o['reader_R18']['n_carriers']} part {o['reader_R18']['partition']}"
        )
        log(
            f"  closure lattice r9: {[(c['cells'], round(c['rho'], 2), round(c['z'], 2)) for c in o['closure_lattice_r9']]}"
        )
        log(
            f"  closure sphere r9: {[(round(c['solid_angle_frac'], 4), round(c['rho'], 2), round(c['z'], 2)) for c in o['closure_sphere_r9']]}"
            f" min gap {o['min_gap_pair_on_sphere_r9']:.4f}; r18 regions {o['closure_sphere_r18_count']}"
        )

    log("(b) control")
    J["b_control"] = audit_b_control()
    log(f"  {J['b_control']}")
    J["b_rows"] = {}
    for name in list(S1_ROWS) + [branch_row(p) for p in BRANCHES]:
        log(f"(b) row {name}")
        J["b_rows"][name] = audit_b(name, claims)
        o = J["b_rows"][name]
        log(
            f"  net {o['net_free']:.4f} abs {o['abs_free']:.3f} (R26-3 {o['abs_free_R26_3']:.3f}, corr {o['corr_own_vs_R26_3']:.3f})"
            f" core {o['frac_abs_core']:.4f} tubes {o['frac_abs_tubes']:.4f} rest {o['frac_abs_rest']:.4f}"
            f" net_outer {o['net_outer']:.4f} r2 {o['r2_abs_weighted']:.1f} ill {o['frac_abs_in_ill_defined_director_cells']:.3f}"
        )
        log(
            f"  degree(R) outward lift: {[(d['R'], round(d['degree'], 3), round(d['min_abs_align'], 3)) for d in o['sphere_degree']]}"
        )
        log(
            f"  |degree|(R) continuous lift (frustrated edges, min gap_dir): {[(d['R'], round(d['abs_degree_lift'], 3), d['frustrated_edges'], round(d['min_gap_dir'], 3)) for d in o['sphere_degree_lift']]}"
        )
        l3 = o["lift3d"]
        log(
            f"  lift3d: frustrated cells {l3['frustrated_cells']} (free {l3['frustrated_cells_in_free']}, min gap_dir there {l3['min_gap_dir_at_frustrated']})"
            f" net {l3['net_free']:.4f} abs {l3['abs_free']:.3f} net_in_r6 {l3['net_inside_r6']:.4f} core {l3['frac_abs_core']:.4f} tubes {l3['frac_abs_tubes']:.4f} rest {l3['frac_abs_rest']:.4f}"
            f" near-frustration {l3['frac_abs_within_2_cells_of_frustration']:.3f} r2 {l3['r2_abs_weighted']:.1f} net_in_R {[(x[0], round(x[1], 3)) for x in l3['net_inside_R']]}"
        )
        log(
            f"  abs per cell by shell: {[(s['r'], round(s['abs_per_cell'], 5)) for s in o['shells']]}"
        )

    J["a_branches"] = {}
    for part in ("1_1_1_1", "3_1", "2_2"):
        log(f"(a) branch {part}")
        J["a_branches"][part] = audit_a(part, claims)
        o = J["a_branches"][part]
        for s in o["shells"]:
            u = s["unmasked_w1.5"]
            b = s["balls8_w1.5"]
            if "a1" in u and "a1" in b:
                log(
                    f"  r {s['r']:5.1f} unm a1 {u['a1']:+.2e} (claim {o['claimed_unmasked_a1'].get(str(s['r']), float('nan')):+.2e})"
                    f" l4 {u.get('a1_lmax4', float('nan')):+.2e} proj {s['projection_sphere']['a1']:+.2e}"
                    f" | balls a1 {b['a1']:+.2e} (claim {o['claimed_balls_a1'].get(str(s['r']), float('nan')):+.2e})"
                    f" floors: claimed {b['floor_claimed']:.1e} se {b['se_a1']:.1e} perm {b['floor_perm']:.1e} boot {b['floor_boot']:.1e}"
                    f" shift {b.get('shift_l2_to_l4', float('nan')):+.1e} | eq l1 {np.hypot(u['b1_x'], u['c1_y']):.3f}"
                )
        log(
            f"  slopes: unm {o['slope_a1_unmasked_w1.5_6_18']} (claim {o['claimed_slope_unmasked']}) proj {o['slope_a1_projection_6_18']} signs {o['signs_a1_unmasked_w1.5']}"
        )
        log(
            f"  masked below 3x floor (claimed/se/perm/boot): {[o[k] for k in o if k.startswith('masked_below')]}"
        )

    J["verdicts"], J["summary"] = verdicts(J, claims)
    with open(OUT_JSON, "w") as f:
        json.dump(J, f, indent=1)
    log(f"wrote {OUT_JSON}")
    print(json.dumps(J["summary"], indent=1))


def verdicts(J, claims):
    V = []

    def add(claim, own, claimed, verdict, note):
        V.append(
            {
                "claim": claim,
                "own_value": own,
                "claimed_value": claimed,
                "verdict": verdict,
                "note": note,
            }
        )

    def r1(x, k=3):
        return None if x is None else float(f"{x:.{k}g}")

    # ---------------- (a)
    a = J["a_branches"]["1_1_1_1"]
    sh = {s["r"]: s for s in a["shells"]}
    R5 = (6.0, 9.0, 12.0, 15.0, 18.0)
    own_u = {r: r1(sh[r]["unmasked_w1.5"]["a1"]) for r in R5}
    cl_u = {r: r1(a["claimed_unmasked_a1"][str(r)]) for r in R5}
    rel = max(
        abs(sh[r]["unmasked_w1.5"]["a1"] - a["claimed_unmasked_a1"][str(r)])
        / max(abs(a["claimed_unmasked_a1"][str(r)]), 1e-12)
        for r in R5
    )
    add(
        "(a) {1,1,1,1} unmasked a_1: 2.1e-4 (r 6) to 4e-5 (r 18), positive at 6-7.5, negative at 9-15, positive from 16.5",
        own_u,
        cl_u,
        "CONFIRMED" if rel < 0.02 else "REFUTED",
        f"own nine-harmonic least squares on the same shells reproduces every value (max relative difference {rel:.1e}); sign pattern {a['signs_a1_unmasked_w1.5']} (shells 3 to 21)",
    )
    proj = {r: r1(sh[r]["projection_sphere"]["a1"]) for r in R5}
    shift = {r: r1(sh[r]["unmasked_w1.5"].get("shift_l2_to_l4"), 2) for r in R5}
    add(
        "(a) {1,1,1,1} unmasked a_1 by an independent method (solid-angle P_1 projection on the trilinearly interpolated sphere, no least squares, no lattice shell)",
        {
            "projection": proj,
            "l2_to_l4_shift_of_lsq": shift,
            "lsq_shell_width_1.0": {r: r1(sh[r]["unmasked_w1.0"]["a1"]) for r in R5},
        },
        cl_u,
        "CONFIRMED",
        "the projection agrees with the shell least squares to 10-20 percent at r >= 7.5 (r 6 differs by 2x, the shell straddles the {4} carrier at r 6), the l <= 4 model extension moves the unmasked a_1 by under 2 percent, and the 1.0-wide shell gives the same values: the unmasked axial coefficient is a DETERMINED number of order 1e-4, not noise (this matters for the floor claim below)",
    )
    add(
        "(a) {1,1,1,1} log-log slope -1.98 of |a_1| on 6-18",
        {
            "lsq_w1.5": r1(a["slope_a1_unmasked_w1.5_6_18"]),
            "lsq_w1.0": r1(a["slope_a1_unmasked_w1.0_6_18"]),
            "projection": r1(a["slope_a1_projection_6_18"]),
        },
        r1(a["claimed_slope_unmasked"]),
        "QUALIFIED",
        "the number reproduces exactly, but a_1 crosses zero twice inside the window (between 7.5 and 9, between 15 and 16.5), so a log-log fit of |a_1| is not the exponent of anything: the 1.0-wide shell and the projection give different slopes for the same field",
    )
    win = [sh[r] for r in (6.0, 7.5, 9.0, 10.5, 12.0, 13.5, 15.0, 16.5, 18.0)]
    own_floor = {
        s["r"]: {
            "claimed_conv": r1(s["balls8_w1.5"]["floor_claimed"], 2),
            "se": r1(s["balls8_w1.5"]["se_a1"], 2),
            "perm": r1(s["balls8_w1.5"]["floor_perm"], 2),
            "boot": r1(s["balls8_w1.5"]["floor_boot"], 2),
            "l4_shift_masked": r1(s["balls8_w1.5"].get("shift_l2_to_l4"), 2),
        }
        for s in win
    }
    cl_floor = {s["r"]: r1(a["claimed_floor"][str(s["r"])], 2) for s in win}
    add(
        "(a) the floor rms / sqrt(N) (3.7e-3 to 5.4e-3) is a sound noise estimate for a_1",
        own_floor,
        cl_floor,
        "QUALIFIED",
        "the residual rms is the l >= 3 angular structure of the four carriers, not noise; the normal-equation standard error and the bootstrap are 1.7x the claimed floor (sqrt(3), the cos^2 mean), the permutation null 2x; on the z-symmetric {2,2} field the same least squares returns a_1 = 1e-8, and on {1,1,1,1} the unmasked a_1 is reproduced by the projection to 10-20 percent and moves under 2 percent under l <= 4, so its real uncertainty is ~2e-5, 100x below the floor; the floor is an ALIASING bound (how much unfitted structure can leak into a_1 when the shell is incomplete), which is the right yardstick for the MASKED fits (their l <= 4 shift is 1e-3 to 3e-3, the size of the masked a_1 itself) and a 100x overstatement for the unmasked fit; the conclusion (no axial l = 1 on {1,1,1,1} beyond order 1e-4) is unchanged under every floor",
    )
    add(
        "(a) {1,1,1,1} masked (3h balls) a_1 ~ 1e-3, noisy, below 3 x floor on every shell 6-18 (MOMENT_FOOTPRINT first condition)",
        {
            "balls_own_reader_8_spheres": {s["r"]: r1(s["balls8_w1.5"]["a1"]) for s in win},
            "balls_own_13_spheres": {s["r"]: r1(s["balls13_w1.5"].get("a1")) for s in win},
        },
        {s["r"]: r1(a["claimed_balls_a1"][str(s["r"])]) for s in win},
        "CONFIRMED",
        f"an own reader run gives the same carriers, so the eight-sphere mask reproduces the claimed values exactly; a thirteen-sphere mask (own variant) gives different values of the same 1e-3 size, which is the noisiness claimed; the condition holds under the claimed floor {a['masked_below_3x_floor_claimed_all_shells_6_18']}, the standard error {a['masked_below_3x_se_a1_all_shells_6_18']}, the permutation null {a['masked_below_3x_floor_perm_all_shells_6_18']}, the bootstrap {a['masked_below_3x_floor_boot_all_shells_6_18']}",
    )
    eq = {
        s["r"]: [
            r1(s["unmasked_w1.5"]["b1_x"], 3),
            r1(s["unmasked_w1.5"]["c1_y"], 3),
            r1(s["projection_sphere"]["b1_x"], 3),
            r1(s["projection_sphere"]["c1_y"], 3),
        ]
        for s in win
    }
    add(
        "(a) substance: no l = 1 signal above the reader floor on {1,1,1,1}, masked or unmasked",
        {
            "b1_x_c1_y_lsq_then_projection_by_shell": eq,
            "equatorial_over_se_6_18": [
                round(x, 1) for x in a["equatorial_l1_over_se_unmasked_6_18"]
            ],
        },
        "no l = 1 signal",
        "QUALIFIED",
        "true for the m = 0 (axial) component only: the equatorial l = 1 vector (b_1, c_1) is 0.02 to 0.06 on r 6-7.5 and 13.5-21 (near zero, 0.004-0.013, only on r 9-12), with b_1 = -c_1 to 1 percent on every shell (a dipole along (1, -1, 0), the reflection symmetry of the S1 seed) and reproduced by the projection, so it is a determined moment of the field, not aliasing (3-6 residual-based sigma where it is large); the branch carries a non-axial l = 1 moment that an axial read cannot see; the pre-registered rule reads a_1 only, so the LABEL is unaffected",
    )
    a3 = J["a_branches"]["3_1"]
    sh3 = {s["r"]: s for s in a3["shells"]}
    add(
        "(a) {3,1}: unmasked a_1 rises from -5e-3 (r 9) to 0.036 (r 15-21); masked 2e-3 to 2.4e-2, above the floor beyond r 12; m0 content 0.15-0.83",
        {
            "unmasked": {
                r: r1(sh3[r]["unmasked_w1.5"]["a1"]) for r in (9.0, 12.0, 15.0, 18.0, 21.0)
            },
            "projection": {
                r: r1(sh3[r]["projection_sphere"]["a1"]) for r in (9.0, 12.0, 15.0, 18.0, 21.0)
            },
            "balls": {r: r1(sh3[r]["balls8_w1.5"]["a1"]) for r in (12.0, 15.0, 18.0)},
            "balls_over_se": {
                r: round(sh3[r]["balls8_w1.5"]["a1"] / sh3[r]["balls8_w1.5"]["se_a1"], 1)
                for r in (12.0, 15.0, 18.0)
            },
            "m0_content_unmasked": {
                r: round(sh3[r]["unmasked_w1.5"]["m0_content"], 2)
                for r in (13.5, 15.0, 18.0, 21.0)
            },
        },
        {
            "unmasked": {
                r: r1(a3["claimed_unmasked_a1"][str(r)]) for r in (9.0, 12.0, 15.0, 18.0, 21.0)
            },
            "balls": {r: r1(a3["claimed_balls_a1"][str(r)]) for r in (12.0, 15.0, 18.0)},
        },
        "CONFIRMED",
        "reproduced exactly; the projection route confirms the unmasked axial moment (0.035 at r 15), so the {3,1} rise is a field property; the masked values sit at 3 standard errors beyond r 12 (6x the claimed floor), a weaker significance than the claimed floor suggests but still above it",
    )
    a2 = J["a_branches"]["2_2"]
    sh2 = {s["r"]: s for s in a2["shells"]}
    add(
        "(a) {2,2}: unmasked a_1 ~ 1e-8 (z-symmetric), masked ~ 1e-3 (the mask's asymmetry)",
        {
            "unmasked_max_abs": r1(max(abs(sh2[r]["unmasked_w1.5"]["a1"]) for r in sh2), 2),
            "balls": {r: r1(sh2[r]["balls8_w1.5"]["a1"]) for r in (9.0, 15.0, 18.0)},
        },
        {
            "unmasked": "1e-8",
            "balls": {r: r1(a2["claimed_balls_a1"][str(r)]) for r in (9.0, 15.0, 18.0)},
        },
        "CONFIRMED",
        "reproduced; this row shows the least squares has a 1e-8 numerical floor and that a masked a_1 of 1e-3 is the mask's footprint, not the field's",
    )
    # ---------------- (b)
    b = J["b_rows"]
    cl_b = claims["b_charge_split"]["rows"]
    add(
        "(b) |rho| fractions: core 0.001-0.024, tubes 0.001-0.027, rest 0.95-0.998 (S1 rows + certified branches); {2,1,1} tubes 0.14 core 0.04; {4} tubes 0.21 core 0.10",
        {
            k: {
                "core": r1(v["frac_abs_core"], 2),
                "tubes": r1(v["frac_abs_tubes"], 2),
                "rest": r1(v["frac_abs_rest"], 3),
                "core_lift": r1(v["lift3d"]["frac_abs_core"], 2),
                "tubes_lift": r1(v["lift3d"]["frac_abs_tubes"], 2),
            }
            for k, v in b.items()
        },
        {
            k: {
                "core": r1(cl_b[k]["fraction_inside_core"], 2),
                "tubes": r1(cl_b[k]["fraction_outer_in_tubes"], 2),
                "rest": r1(cl_b[k]["fraction_outer_outside_tubes"], 3),
            }
            for k in b
        },
        "CONFIRMED",
        "own Jacobian (np.gradient central differences, correlation 1.000 with the R26-3 density) and own tubes from an own reader run reproduce the split; the split is the same under a continuous (sign-propagated) lift of the director, so it does not depend on the outward orientation",
    )
    add(
        "(b) net charge outside r 6 is 0.97-0.999 of a unit; the hedgehog seed control has 0.964 inside r 6; the reader floor abs 0.036 outside r 6",
        {
            "net_outer_outward_lift": {k: r1(v["net_outer"], 3) for k, v in b.items()},
            "net_inside_r6_continuous_lift": {
                k: r1(v["lift3d"]["net_inside_r6"], 2) for k, v in b.items()
            },
            "own_hedgehog_control": {k: r1(v, 3) for k, v in J["b_control"].items()},
        },
        {k: r1(cl_b[k]["net_charge_outer"], 3) for k in b},
        "CONFIRMED",
        "reproduced by the outward-lift Jacobian, by the continuous-lift Jacobian (|net inside r 6| <= 0.013 on the S1 rows and certified branches, 0.07 on {4}), and by the lift-independent sphere degree (0 on r 3 with no frustrated edge on every row); the own uniaxial hedgehog control gives 0.973 inside r 6 and 0.037 |rho| outside",
    )
    r2k = ["S1_d0.3_w25_n32_L48", "S1_d0.3_w25_n48_L72", "S1_d0.3_w25_n64_L96"]
    add(
        "(b) r^2 moments of |rho|: 205 (L 48), 347 (L 72), 394 (L 96) at h 1.5",
        {
            k: {
                "outward": r1(b[k]["r2_abs_weighted"], 4),
                "continuous_lift": r1(b[k]["lift3d"]["r2_abs_weighted"], 4),
                "abs_total": r1(b[k]["abs_free"], 3),
            }
            for k in r2k
        },
        {k: r1(cl_b[k]["r2_abs_weighted_all"], 4) for k in r2k},
        "QUALIFIED",
        "the numbers reproduce (and the continuous lift gives the same to 2 percent), but they are moments of |rho| whose free-volume total is 3.6 (L 48), 6.8 (L 72), 7.0 (L 96) units against a net of 1, growing with the box and with resolution at fixed L (3.6 at h 1.5 to 5.6 at h 1 on L 48): this |rho| is the Jacobian of a LINE field near its disclinations (22-80 percent of it lies within two cells of a frustrated bond of the continuous lift), not a charge density whose moment converges",
    )
    add(
        "(b) the charge outside the core is a property of the field (the director not radial inside r 9, R26-4), not an artifact of the discrete Jacobian",
        {
            "abs_degree_continuous_lift_by_sphere_(R,|deg|,frustrated_edges,min_gap_dir)": {
                k: [
                    [
                        d["R"],
                        round(d["abs_degree_lift"], 3),
                        d["frustrated_edges"],
                        round(d["min_gap_dir"], 3),
                    ]
                    for d in v["sphere_degree_lift"]
                ]
                for k, v in b.items()
            },
            "outward_lift_degree_min_abs_align": {
                k: [
                    [d["R"], round(d["degree"], 3), round(d["min_abs_align"], 3)]
                    for d in v["sphere_degree"]
                ]
                for k, v in b.items()
            },
            "lift3d_frustrated_free_cells": {
                k: v["lift3d"]["frustrated_cells_in_free"] for k, v in b.items()
            },
            "frac_abs_rho_within_2_cells_of_frustration": {
                k: r1(v["lift3d"]["frac_abs_within_2_cells_of_frustration"], 2)
                for k, v in b.items()
            },
        },
        "property of the field",
        "QUALIFIED",
        "the LOCATION of the net charge is a property of the field: under a continuous lift the director's degree is exactly 0 on r 3 with zero frustrated edges on every row (r 4.5 too at h 1.5, and out to r 10.5 on {3,1}: the core is a topologically trivial, well-defined line field), and exactly 1 from r 12 on the S1 rows at h 1.5 (10.5 at h 1) and from r 15-18 on the branches; but between those radii disclination lines pierce every sphere (director gap < 0.01, 16-495 frustrated edges), so the director is not a vector field there and a charge DENSITY is not defined: the outward-lift degrees inside r 12 (-5, -3, 3, 0.9, 1.99, 2.97) are meaningless (min |align| = 0 on every such sphere) and the Jacobian's |rho| there is the finite-difference response to the disclinations, 1-6 percent of all cells being frustrated; the split core/tubes/rest is nevertheless lift-independent, so 'the charge is not in the core and not in the tubes' stands while 'the |rho| distribution' should not be read as where charge sits",
    )
    # ---------------- (c)
    c = J["c_rows"]
    add(
        "(c) on r 9 the S1 rows read four half-unit carriers at rho 7.7-8.4, z +-3.2 to 4.6, box-independent; {1,1,1,1} at rho 4.4-8.2",
        {
            "reader_rho": {
                k: [round(x[0], 2) for x in v["reader_R9"]["rho_hu_z"] if x[0] is not None]
                for k, v in c.items()
            },
            "reader_partition": {k: v["reader_R9"]["partition"] for k, v in c.items()},
            "gap_closure_regions_on_sphere_(rho,z,phi)": {
                k: [
                    [round(x["rho"], 2), round(x["z"], 2), round(x["phi_deg"])]
                    for x in v["closure_sphere_r9"]
                ]
                for k, v in c.items()
            },
        },
        "rho 7.7-8.4 (S1), 4.4-8.2 ({1,1,1,1})",
        "CONFIRMED",
        "the reader reproduces on every row; the independent gap-closure read (regions of pair gap < 0.1 delta on the interpolated r 9 sphere) finds the melted cores at rho 7.9-8.4, z +-3.2 to 4.3 on every S1 row (and at rho 4.6-4.7 / 8.2 on {1,1,1,1}), so the carrier positions are a field property; caveat: at h 1.5 the closure read also finds four extra melted spots at rho 6.9-7.4 (z +-5.1 to 5.9) that carry no winding and are absent at h 1, a resolution-dependent feature the reader does not report",
    )
    add(
        "(c) on r 18 the reader returns 14-26 fragments of +-1 and -2",
        {
            "reader_fragments": {k: v["reader_R18"]["n_carriers"] for k, v in c.items()},
            "gap_closure_regions_r18": {k: v["closure_sphere_r18_count"] for k, v in c.items()},
        },
        "14 to 26",
        "CONFIRMED",
        "14, 26, 24, 14, 16 on the S1 rows (8 on {1,1,1,1}); the closure read finds 18-54 disjoint melted regions on the r 18 sphere of the S1 rows, so the fragmentation is in the field, not only in the reader",
    )
    # ---------------- (d)
    hh = J["d_hedgehog"]
    Cs = {c[0]: round(c[1], 3) for c in hh["C_lattice_per_shell_cellwise_e_r4"]}
    add(
        "(d) C = 8 (sympy) on M = n n^T, c = 32 pi = 100.53, lattice fit 100.61 (R_eff 39.5) on n 48 L 72",
        {
            "C_sympy": J["d_sympy"]["C_symbolic"],
            "c_fit_all_cells": r1(hh["c_fit_all_9_18"], 5),
            "R_eff": r1(hh["R_eff_all"], 3),
            "c_fit_spherical_window": r1(hh["c_fit_sphere_9_18"], 5),
            "C_lattice_cellwise_r9_r18": [Cs[9.0], Cs[18.0]],
            "local_exponent_range": [
                round(min(e[2] for e in hh["local_exponents"]), 2),
                round(max(e[2] for e in hh["local_exponents"]), 2),
            ],
            "V4_share": hh["V4_share"],
        },
        {
            "C": 8.0,
            "c_4piC": 100.53,
            "c_fit": 100.61,
            "R_eff": 39.5,
            "C_lattice_r9_r18": [8.23, 8.03],
        },
        "CONFIRMED",
        "own sympy gives 8 exactly; own lattice fit reproduces c and R_eff; the cellwise e r^4 reads 8.15 at r 9 falling to 8.03 at r 18 (the stencil's 2 percent bias), the local exponent is -4.0 to -4.2 on every shell; the potential vanishes identically on the hedgehog",
    )
    d = J["d_rows"]
    add(
        "(d) biaxial shell-density exponent on 9-18: -1.0 to -1.2 at h 1.5, -1.8 to -1.9 at h 1 (a Coulomb exterior has -4, the R23 halo -2)",
        {
            k: {
                "exp_9_18": r1(v["exp_9_18"], 3),
                "peak_r": v["r_of_density_peak"],
                "local_exp_pairs_9_to_18": [
                    [e[0], round(e[2], 2)] for e in v["local_exponents"] if 9.0 <= e[0] < 18.0
                ],
                "exp_13.5_24": r1(v["exp_13.5_24"], 3),
            }
            for k, v in d.items()
        },
        {k: r1(v["claimed"]["shell_density_exponent_on_9_18"], 3) for k, v in d.items()},
        "QUALIFIED",
        "the numbers reproduce, but the density is not a power law on 9-18: at h 1.5 it RISES from r 9 to a peak at 10.5-12 (local exponent +0.9 to +1.2) and falls beyond with local exponents -1.5 to -3.5, at h 1 the peak sits at r 9 and the window is all fall, so the h 1.5 versus h 1 difference (-1.0 versus -1.9) is the peak's resolution shift, not two exponents; fitted beyond the peak (13.5-24) every row gives -2.9 to -3.3, and beyond r 24 on the two large boxes the local exponent is near -4",
    )
    L96 = d["S1_d0.3_w25_n64_L96"]
    L72 = d["S1_d0.3_w25_n48_L72"]
    add(
        "(d) c(R) = E(> R) / (1/R - 1/R_max) drifts by 1.5-2.4 across 9-18 (109 to 247 on L 48; 89 to 135 on L 96), then flat at 139 on R 24-30",
        {
            "drift_author_convention": {k: r1(v["c_author_9_18_drift"], 3) for k, v in d.items()},
            "drift_spherical_window_true_R_free": {
                k: r1(v["c_sphere_9_18_drift"], 3) for k, v in d.items()
            },
            "L96_c_24_30_author": [round(x, 1) for x in L96["c_author_24_30"]],
            "L96_c_24_30_sphere": [round(x, 1) for x in L96["c_sphere_24_30"]],
            "L72_c_24_30_author": [round(x, 1) for x in L72["c_author_24_30"]],
            "L72_c_24_30_sphere": [round(x, 1) for x in L72["c_sphere_24_30"]],
            "R_max_author_vs_true_free": {
                k: [v["R_max_author"], v["R_free"]] for k, v in d.items()
            },
        },
        {k: r1(v["claimed"]["c_drift_9_to_18"], 3) for k, v in d.items()},
        "CONFIRMED",
        "the drift reproduces with the author's R_max = L/2 - 1.6 and E(> R) over the whole free cube, and persists (1.5-2.1) with the true free radius L/2 - 2h and a spherical window; the L 96 flatness at 139 (author convention) becomes a 4 percent fall from 132 to 127 in the spherical window, and the L 72 box has its own flat-to-4-percent window at 172-179 on the same radii: the flat value is box-dependent",
    )
    s72 = {s[0]: s[1] for s in L72["shell_density"]}
    s96 = {s[0]: s[1] for s in L96["shell_density"]}
    ratio = {r: round(s72[r] / s96[r], 3) for r in (9.0, 12.0, 18.0, 24.0, 27.0, 30.0, 33.0)}
    add(
        "(d) 'no calibration readable' on the biaxial rows (c read on the uniaxial hedgehog only)",
        {
            "L72_over_L96_density_at_equal_r": ratio,
            "L96_exp_24_33": r1(L96["exp_24_33"], 3),
            "L72_exp_24_33": r1(L72["exp_24_33"], 3),
            "L96_c_24_30_over_32pi": [round(x / (32 * np.pi), 2) for x in L96["c_author_24_30"]],
            "L72_c_24_30_over_32pi": [round(x / (32 * np.pi), 2) for x in L72["c_author_24_30"]],
        },
        "no calibration readable",
        "CONFIRMED",
        "fair, and for a stronger reason than the drift: beyond r 24 both large boxes show a local density exponent near -4 (L 72: -3.8, L 96: -4.3 on 24-33) but their densities differ by 15-36 percent at equal r and their would-be c by 30 percent (172-179 versus 127-140, i.e. 1.7-1.8 versus 1.3-1.4 times 32 pi), so the exponent window exists but its amplitude is set by the box; a Coulomb-like far zone is not excluded, it is unread until a larger box converges the amplitude",
    )
    n_conf = sum(v["verdict"] == "CONFIRMED" for v in V)
    n_qual = sum(v["verdict"] == "QUALIFIED" for v in V)
    n_ref = sum(v["verdict"] == "REFUTED" for v in V)
    summary = {
        "counts": {"CONFIRMED": n_conf, "QUALIFIED": n_qual, "REFUTED": n_ref, "total": len(V)},
        "refuted": [
            v["claim"] + " :: " + v["note"].split(";")[0] for v in V if v["verdict"] == "REFUTED"
        ],
        "qualified": [
            v["claim"] + " :: " + v["note"].split(";")[0] for v in V if v["verdict"] == "QUALIFIED"
        ],
        "blind_spots": [
            "the partition reader is used as an instrument (allowed, R26-0 audited); the carrier POSITIONS were cross-checked by the gap-closure read, the half-unit VALUES were not re-derived",
            "the sphere-degree and gap-closure routes interpolate the same lattice field trilinearly; the disclination structure between r 6 and 12 is read at one resolution per row and its h dependence (the extra windingless closures at h 1.5) is noted, not resolved",
            "the 3D continuous lift propagates along straight lines, so frustrated cells are where the lift path meets a disclination, which may overcount the disclination volume; the 2D sphere lifts are path-independent on the unfrustrated spheres",
            "the charge density on the capped branches {2,1,1} and {4} is read on fields the census did not certify",
            "no box larger than L 96 exists, so whether the far-zone exponent near -4 on 24-33 survives box convergence is untested",
            "the floors are per-shell statistics; a joint test across shells (a_1 correlated in r) was not built; the equatorial l = 1 significance is quoted against the same residual-based error that the audit calls an overstatement",
        ],
    }
    return V, summary


if __name__ == "__main__":
    main()
