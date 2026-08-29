#!/usr/bin/env python3
"""M5.32 R12 ADVERSARIAL AUDIT of the charged disclination RING.

Independent re-implementation (AI_HYGIENE.md par. 1). This file does NOT
import, read or execute m5_32_r12_a_ring.py; the seed, the winding reader,
the far-surface degree, the clock generator and the sub-grid cord locator
are all rebuilt here from the equations below. The certified stack
m5_21_3_a_4d.py (INS4) is imported for the SHARED instruments only
(e_parts, kin_of, fire, pin_shell), because those are the objects under
test, not the objects being re-derived.

=========================== EQUATIONS FIRST ===========================

(0) FIELD AND VACUUM
    M(x) is a symmetric 4x4 field on an n^3 grid, x_i = (i - (n-1)/2) h,
    h = L/n. eta = diag(-1, 1, 1, 1). s = -1, g = 8, delta = 0.3, so
        vac4 = diag(-s g, 1, delta, 0) = diag(8, 1, 0.3, 0).
    A_i = d_i M / h (sym = mean of fwd and bwd one-sided stencils).

(1) STATIC ENERGY (INS4.e_parts)
    E_u = h^3 . 4 sum_{i<j} sum_cells <[A_i, A_j]_eta, [A_i, A_j]_eta>_eta
    with [A,B]_eta = A eta B - B eta A and <F,G>_eta = tr(eta F eta G^T).
    V4  = h^3 W1 sum_cells sum_{p=1..4} (tr((M eta)^p) - c_p)^2,
    c_p = (s g)^p + 1 + delta^p, W1 = 7.24023879e-4.
    E_stat = E_u + V4.

(2) THE RING SEED (rebuilt here from the engine1_seeds docstring)
    rho = sqrt(x^2 + y^2), rho_hat = (x, y, 0)/rho
    psi = 1/2 [ atan2(-z, rho - a) + atan2(-z, rho + a) ] + pi/2
    n   = sin(psi) rho_hat + cos(psi) z_hat
    e_phi = smoothstep(rho/rho_c) . (-y, x, 0)/rho, then projected
            perpendicular to n (NOT renormalized: it melts on the axis)
    e_theta = e_phi x n
    d_c = sqrt((rho - a)^2 + z^2),  smelt = smoothstep(d_c / w_c)
    d_iso = (1 + delta)/3
    (d0, d1, d2) = d_iso + smelt . ((1, delta, 0) - d_iso)
    M_sp = d0 n n^T + d1 e_theta e_theta^T + d2 e_phi e_phi^T
    M    = blockdiag(-s g, M_sp)
    with a the ring radius, w_c = rho_c = 3.0 physical units.
    a = 0 reduces EXACTLY to the radial point hedgehog (psi = pi/2 - beta).

(3) CLOCK GENERATOR (generalized, per-cell)
    n1(x) = leading eigenvector of the spatial 3x3 block of M.
    J(x): spatial block J_ab = eps_abc n1_c, time row/column zero.
    a0 = J M - M J   (= dM/dt for a rigid rotation about n1; symmetric,
    since J is antisymmetric and M symmetric). Defined up to the
    per-cell sign of n1, which cancels in kin (pointwise quadratic).
    TAPERED: a0 <- w(r) a0, w = 1 for r <= r0, linear to 0 at r1;
    the record convention is (r0, r1) = (12, 15).

(4) KINETIC FORM (INS4.kin_of)
    kin(M; a0) = h^3 . 4 sum_i sum_cells <[a0, A_i]_eta, [a0, A_i]_eta>_eta
    E_kin(omega) = omega^2 kin. Per-cell density used here for the shell
    profiles is the same summand before the cell sum.

(5) FIXED-J ENSEMBLE (M5.21.15 algebra)
    E(omega) = E_stat + omega^2 kin;  J = dE/domega = 2 omega kin
    => omega*(J) = J / (2 kin),  E_J = E_stat + J^2 / (4 kin).

(6) MERIDIONAL WINDING (my reader, eigenvector form -- NOT the record's
    atan2(2 M_xz, M_xx - M_zz) form, which is also computed for contrast)
    On a circle of radius R about the cord point (rho = a, z = 0) in the
    y = 0 half-plane, M is trilinearly interpolated, n1 read by eigh, and
    alpha = atan2(n1_z, n1_x). The DIRECTOR is a line field, so the
    single-valued object is 2 alpha; accumulate d(2 alpha) wrapped into
    (-pi, pi] and set q = (total) / (4 pi). q = 1/2 is a half winding.

(7) FAR-SURFACE DEGREE (my reader, sphere + Jacobian; NOT the record's
    cube + solid-angle-of-quads reader, which is also run for contrast)
    On a sphere of radius R, v(theta, phi) = leading eigenvector sign-fixed
    by sign(v . r_hat) (legitimate only while |v . r_hat| stays away from
    0, which is reported), and
        deg = (1/4 pi) INT v . (d_theta v x d_phi v) dtheta dphi.

(8) SUB-GRID CORD LOCATOR (the C2 attack)
    Three independent scalars along the +x ray at y = z = 0:
      gap12 = l1 - l2, gap13 = l1 - l3, and l1 itself (l1 >= l2 >= l3 the
      spatial eigenvalues). Each has an extremum at the melted cord.
      (a) PARABOLA: exact 3-point vertex fit through the discrete
          extremum on the grid radii (no x-interpolation).
      (b) FINE: trilinear interpolation on a 0.01-spaced ray.
    Reported for the seed and for the it1500 / it3000 relaxed fields, so
    the DRIFT is measured at ~0.01 resolution instead of the 1.5 grid step.

============================== THE CHECKS =============================
C1 instrument : half winding on MY circle radii + MY far-surface degree
                + a0 identity against B8.a0_unit on the hedgehog.
C2 statics    : sub-grid cord drift, E_u/V4 recomputed, energy still
                falling at 3000, ring-vs-hedgehog lead at MATCHED and at
                LONGER hedgehog relaxation, far-field pin integrity,
                unwinding-vs-motion.
C3 clock      : rigid/tapered kin recomputed, n48-vs-n32 interior
                identity (tautology test), tapered shell localization,
                hedgehog tapered comparison.
C4 fixed J    : seed ladder EXTENDED to a = 0 (the hedgehog limit) and
                w_c varied to 2 and 4; taper-artifact test with a wider
                taper (16, 20).

================================ OUTPUTS =============================
../data/m5_32_r12_audit.json   (all numbers + verdict fields)
Run from scripts/:  python3 m5_32_r12_audit.py <stage>
  stages: c1 | c2 | c3 | c4 | probe | rerelax | verdicts
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
CKPT = os.path.join(HERE, "..", "checkpoints", "m5_32_r12")
CKPT10 = os.path.join(HERE, "..", "checkpoints", "m5_32_r10")
OUT = os.path.join(DATA, "m5_32_r12_audit.json")
SCRATCH = os.path.join(HERE, "..", "checkpoints", "m5_32_r12_audit")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


INS4 = _load("ins4", os.path.join(HERE, "m5_21_3_a_4d.py"))
B8 = _load("b8", os.path.join(HERE, "m5_21_8_b_lattice.py"))
EAUD = _load("eaud", os.path.join(HERE, "m5_22_e_audit.py"))

DELTA = 0.3


# ============================ (2) the seed ============================
def smoothstep(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def ring_seed(cfg, a, w_c=3.0, rho_c=None):
    """The charged disclination ring, rebuilt from the eqn (2) block."""
    rho_c = w_c if rho_c is None else rho_c
    X, Y, Z = INS4.coords(cfg["n"], cfg["h"])
    d = cfg["delta"]
    rho = np.sqrt(X * X + Y * Y)
    ainv = 1.0 / np.sqrt(rho * rho + 1e-12)
    rh = np.stack([X * ainv, Y * ainv, np.zeros_like(X)], axis=-1)
    az = np.stack([-Y * ainv, X * ainv, np.zeros_like(X)], axis=-1)
    psi = 0.5 * (np.arctan2(-Z, rho - a) + np.arctan2(-Z, rho + a)) \
        + 0.5 * np.pi
    sp, cp = np.sin(psi), np.cos(psi)
    nv = np.stack([sp * rh[..., 0], sp * rh[..., 1], cp], axis=-1)
    ephi = az * smoothstep(rho / rho_c)[..., None]
    ephi = ephi - np.sum(ephi * nv, axis=-1)[..., None] * nv
    eth = np.cross(ephi, nv)
    dc = np.sqrt((rho - a) ** 2 + Z * Z)
    sm = smoothstep(dc / w_c)
    d_iso = (1.0 + d) / 3.0
    d0 = d_iso + sm * (1.0 - d_iso)
    d1 = d_iso + sm * (d - d_iso)
    d2 = d_iso + sm * (0.0 - d_iso)
    Msp = (d0[..., None, None] * nv[..., :, None] * nv[..., None, :]
           + d1[..., None, None] * eth[..., :, None] * eth[..., None, :]
           + d2[..., None, None] * ephi[..., :, None] * ephi[..., None, :])
    M = np.zeros(X.shape + (4, 4))
    M[..., 1:, 1:] = Msp
    M[..., 0, 0] = -cfg["sg"]
    return INS4.sym4(M)


# ====================== (3) the clock generator =======================
def leading_vec(M):
    _, V = np.linalg.eigh(M[..., 1:4, 1:4])
    return V[..., :, 2]


def a0_general(M, taper=None, cfg=None):
    """a0 = J M - M J with J_ab = eps_abc n1_c in the spatial block."""
    v = leading_vec(M)
    J = np.zeros(M.shape)
    n1, n2, n3 = v[..., 0], v[..., 1], v[..., 2]
    # eps_abc n_c on the spatial block (rows/cols 1..3)
    J[..., 1, 2], J[..., 2, 1] = n3, -n3
    J[..., 1, 3], J[..., 3, 1] = -n2, n2
    J[..., 2, 3], J[..., 3, 2] = n1, -n1
    a0 = J @ M - M @ J
    if taper is not None:
        X, Y, Z = INS4.coords(cfg["n"], cfg["h"])
        r = np.sqrt(X * X + Y * Y + Z * Z)
        r0, r1 = taper
        w = np.clip((r1 - r) / (r1 - r0), 0.0, 1.0)
        a0 = a0 * w[..., None, None]
    return a0


def kin_density(M, a0, cfg):
    """per-cell kin summand (eqn 4 before the cell sum), h^3-weighted."""
    h3 = cfg["h"] ** 3
    dens = np.zeros(M.shape[:3])
    for br, (A, wt) in INS4.a_fields(M, cfg).items():
        for i in range(3):
            F = INS4.comm_eta(a0, A[i])
            dens = dens + wt * 4.0 * INS4.inner_eta(F, F)
    return h3 * dens


# ========================= interpolation ==============================
def trilerp(M, cfg, pts):
    """M at physical points pts (..., 3); nearest-clamped at the border."""
    n, h = cfg["n"], cfg["h"]
    t = np.asarray(pts) / h + (n - 1) / 2.0
    i0 = np.floor(t).astype(int)
    f = t - i0
    i0 = np.clip(i0, 0, n - 2)
    f = np.clip(f, 0.0, 1.0)
    out = np.zeros(t.shape[:-1] + (4, 4))
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                w = ((f[..., 0] if dx else 1 - f[..., 0])
                     * (f[..., 1] if dy else 1 - f[..., 1])
                     * (f[..., 2] if dz else 1 - f[..., 2]))
                out += w[..., None, None] * M[i0[..., 0] + dx,
                                              i0[..., 1] + dy,
                                              i0[..., 2] + dz]
    return out


# ==================== (6) my meridional winding =======================
def winding_eig(M, cfg, a, R, npt=1440):
    """q from the LEADING EIGENVECTOR angle (my reader)."""
    t = np.linspace(0.0, 2.0 * np.pi, npt, endpoint=False)
    pts = np.stack([a + R * np.cos(t), np.zeros_like(t), R * np.sin(t)],
                   axis=-1)
    Mi = trilerp(M, cfg, pts)
    v = leading_vec(Mi)
    al = np.arctan2(v[..., 2], v[..., 0])
    two = 2.0 * al
    d = np.diff(np.concatenate([two, two[:1]]))
    d = (d + np.pi) % (2.0 * np.pi) - np.pi
    return float(np.sum(d) / (4.0 * np.pi))


def winding_atan2M(M, cfg, a, R, npt=1440):
    """the record's form: theta = 1/2 atan2(2 M_xz, M_xx - M_zz)."""
    t = np.linspace(0.0, 2.0 * np.pi, npt, endpoint=False)
    pts = np.stack([a + R * np.cos(t), np.zeros_like(t), R * np.sin(t)],
                   axis=-1)
    Mi = trilerp(M, cfg, pts)
    two = np.arctan2(2.0 * Mi[..., 1, 3], Mi[..., 1, 1] - Mi[..., 3, 3])
    d = np.diff(np.concatenate([two, two[:1]]))
    d = (d + np.pi) % (2.0 * np.pi) - np.pi
    return float(np.sum(d) / (4.0 * np.pi))


def aniso_on_circle(M, cfg, a, R, npt=360):
    t = np.linspace(0.0, 2.0 * np.pi, npt, endpoint=False)
    pts = np.stack([a + R * np.cos(t), np.zeros_like(t), R * np.sin(t)],
                   axis=-1)
    lam = np.linalg.eigvalsh(trilerp(M, cfg, pts)[..., 1:4, 1:4])
    return float(np.mean(lam[..., 2] - lam[..., 0])), \
        float(np.min(lam[..., 2] - lam[..., 0]))


# ==================== (7) my far-surface degree =======================
def degree_sphere(M, cfg, R, nth=181, nph=360):
    th = np.linspace(1e-4, np.pi - 1e-4, nth)
    ph = np.linspace(0.0, 2.0 * np.pi, nph, endpoint=False)
    T, P = np.meshgrid(th, ph, indexing="ij")
    rhat = np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P),
                     np.cos(T)], axis=-1)
    Mi = trilerp(M, cfg, R * rhat)
    v = leading_vec(Mi)
    dot = np.sum(v * rhat, axis=-1)
    v = v * np.sign(dot)[..., None]
    dth = th[1] - th[0]
    dph = ph[1] - ph[0]
    dvt = np.gradient(v, dth, axis=0)
    dvp = (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1)) / (2 * dph)
    jac = np.sum(v * np.cross(dvt, dvp), axis=-1)
    deg = float(np.sum(jac) * dth * dph / (4.0 * np.pi))
    return deg, float(np.min(np.abs(dot)))


# =============== (8) sub-grid cord locator on the +x ray ==============
def ray_scalars(M, cfg, xs):
    pts = np.stack([xs, np.zeros_like(xs), np.zeros_like(xs)], axis=-1)
    lam = np.linalg.eigvalsh(trilerp(M, cfg, pts)[..., 1:4, 1:4])
    l3, l2, l1 = lam[..., 0], lam[..., 1], lam[..., 2]
    return {"gap12": l1 - l2, "gap13": l1 - l3, "l1": l1}


def parabola_min(xs, ys):
    k = int(np.argmin(ys))
    if k == 0 or k == len(ys) - 1:
        return float(xs[k]), int(k), False
    a, b, c = ys[k - 1], ys[k], ys[k + 1]
    den = (c - 2 * b + a)
    if abs(den) < 1e-300:
        return float(xs[k]), int(k), False
    dx = xs[1] - xs[0]
    return float(xs[k] - 0.5 * dx * (c - a) / den), int(k), True


def cord_locate(M, cfg, lo, hi):
    """(a) parabola on grid radii, (b) fine 0.01 trilinear scan."""
    n, h = cfg["n"], cfg["h"]
    gx = (np.arange(n) - (n - 1) / 2.0) * h
    gx = gx[(gx >= lo) & (gx <= hi)]
    sg = ray_scalars(M, cfg, gx)
    fx = np.arange(lo, hi + 1e-9, 0.01)
    sf = ray_scalars(M, cfg, fx)
    out = {}
    for key in ("gap12", "gap13", "l1"):
        xp, kg, ok = parabola_min(gx, sg[key])
        out[key + "_grid_argmin"] = float(gx[int(np.argmin(sg[key]))])
        out[key + "_parab"] = xp if ok else None
        out[key + "_fine"] = float(fx[int(np.argmin(sf[key]))])
    return out


# ============================== stages ================================
def cfg_of(n, L):
    return INS4.base_cfg(s=-1.0, g=8.0, n=n, L=L, delta=DELTA)


def load_ck(name, n):
    A = np.load(os.path.join(CKPT, name))
    assert A.shape == (n, n, n, 4, 4), A.shape
    return A


def st(M, cfg):
    eu, ev = INS4.e_parts(M, cfg)
    return float(eu), float(ev)


def stage_c1():
    t0 = time.time()
    res = {}
    cfg = cfg_of(32, 48.0)
    # --- a0 identity on the m = 0 hedgehog ---
    Mh = B8.dressed(cfg, 0.0)
    a0b = B8.a0_unit(cfg, 0.0)
    a0g = a0_general(Mh)
    s = np.sign(np.sum(a0b * a0g, axis=(-1, -2)))
    s[s == 0] = 1.0
    a0s = a0g * s[..., None, None]
    num = np.max(np.abs(a0s - a0b))
    den = np.max(np.abs(a0b))
    res["a0_identity"] = {
        "max_abs_diff_after_cellsign": float(num),
        "rel": float(num / den),
        "absdiff_of_magnitudes_rel": float(
            np.max(np.abs(np.abs(a0g) - np.abs(a0b))) / den),
        "kin_B8": float(INS4.kin_of(Mh, a0b, cfg)),
        "kin_mine": float(INS4.kin_of(Mh, a0g, cfg)),
        "kin_R10_record_rigid": 426.5070121483972,
    }
    # self-test of my degree reader on the hedgehog
    dh, mdh = degree_sphere(Mh, cfg, 15.0)
    res["degree_selftest_hedgehog"] = {"deg_sphere": dh, "min_abs_dot": mdh}
    # --- the ring seeds ---
    for a in (6.0, 9.0):
        M = ring_seed(cfg, a)
        eu, ev = st(M, cfg)
        mine = {f"R={R}": winding_eig(M, cfg, a, R)
                for R in (1.875, 2.625, 3.375, 4.5)}
        rec = {f"R={R}": winding_atan2M(M, cfg, a, R)
               for R in (2.25, 3.0, 3.75)}
        deg, mdot = degree_sphere(M, cfg, 15.0)
        degc = EAUD.read_charge_from_M(M[..., 1:4, 1:4], 5, cfg["n"] - 6)
        res[f"seed_a{a:g}"] = {
            "E_u": eu, "V4": ev, "E_stat": eu + ev,
            "winding_mine_eigvec": mine,
            "winding_record_form": rec,
            "degree_sphere_R15": deg, "min_abs_dot": mdot,
            "degree_record_cube_reader": [float(degc[0]), int(degc[1])],
        }
        print(f"a={a}: E_u {eu:.4f} V4 {ev:.4f} deg {deg:+.4f} "
              f"wind {mine}", flush=True)
    res["runtime_s"] = round(time.time() - t0, 1)
    return res


CASES = [("n32_L48_a6", 32, 48.0, 6.0), ("n32_L48_a9", 32, 48.0, 9.0),
         ("n48_L72_a6", 48, 72.0, 6.0)]


def stage_c2():
    t0 = time.time()
    res = {}
    for tag, n, L, a in CASES:
        cfg = cfg_of(n, L)
        rec = {}
        fields = {"seed": ring_seed(cfg, a)}
        for it in (1500, 3000):
            fields[f"it{it}"] = load_ck(f"{tag}_it{it}.npy", n)
        lo, hi = max(1.0, a - 4.5), a + 4.5
        for key, M in fields.items():
            eu, ev = st(M, cfg)
            loc = cord_locate(M, cfg, lo, hi)
            wnd = {f"R={R}": winding_eig(M, cfg, a, R)
                   for R in (1.875, 2.625, 3.375, 4.5)}
            an = {f"R={R}": aniso_on_circle(M, cfg, a, R)
                  for R in (1.875, 3.375, 4.5)}
            # far-field pin integrity: spectrum vs vac4 outside r > 0.6 Lhalf
            X, Y, Z = INS4.coords(n, cfg["h"])
            r = np.sqrt(X * X + Y * Y + Z * Z)
            far = r > 18.0          # fixed PHYSICAL radius, both boxes
            lam = np.linalg.eigvalsh(M[far][..., 1:4, 1:4])
            vv = np.array([0.0, DELTA, 1.0])
            rec[key] = {
                "E_u": eu, "V4": ev, "E_stat": eu + ev,
                "cord": loc, "winding": wnd, "aniso_circle": an,
                "far_max_dev_from_vac_spectrum":
                    float(np.max(np.abs(lam - vv))),
                "far_ncells": int(far.sum()),
                "deg_sphere_R15": degree_sphere(M, cfg, 15.0)[0],
            }
            print(f"{tag}/{key}: E_u {eu:.4f} V4 {ev:.4f} "
                  f"cord {loc['gap12_fine']:.3f}/{loc['gap13_fine']:.3f}",
                  flush=True)
        # pinned-cell integrity: relaxed must equal seed on the pin shell
        P = INS4.pin_shell(n, cfg["h"], 1.6)
        rec["pin_max_dev_it3000_vs_seed"] = float(
            np.max(np.abs(fields["it3000"][P] - fields["seed"][P])))
        rec["pinned_frac"] = float(P.mean())
        res[tag] = rec
    # --- the hedgehog comparison, recomputed by me, at 3 depths ---
    cfg = cfg_of(32, 48.0)
    hh = {}
    for it in (3000, 6000, 12000):
        p = os.path.join(CKPT10, f"relax_g8_n32_L48_it{it}.npy")
        if not os.path.exists(p):
            continue
        Mh = np.load(p)
        eu, ev = st(Mh, cfg)
        hh[f"it{it}"] = {"E_u": eu, "V4": ev, "E_stat": eu + ev}
        print(f"hedgehog it{it}: E_u {eu:.4f} V4 {ev:.4f}", flush=True)
    res["hedgehog_R10"] = hh
    res["runtime_s"] = round(time.time() - t0, 1)
    return res


def shell_profile(dens, cfg, edges):
    X, Y, Z = INS4.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    return [float(dens[(r >= a) & (r < b)].sum())
            for a, b in zip(edges[:-1], edges[1:])]


def stage_c3():
    t0 = time.time()
    res = {}
    edges = list(np.arange(0.0, 36.1, 3.0))
    fields = {}
    for tag, n, L, a in CASES:
        cfg = cfg_of(n, L)
        M = load_ck(f"{tag}_it3000.npy", n)
        fields[tag] = (cfg, M)
        ar = a0_general(M)
        at = a0_general(M, taper=(12.0, 15.0), cfg=cfg)
        aw = a0_general(M, taper=(16.0, 20.0), cfg=cfg)
        dr = kin_density(M, ar, cfg)
        dt = kin_density(M, at, cfg)
        res[tag] = {
            "kin_rigid": float(INS4.kin_of(M, ar, cfg)),
            "kin_taper_12_15": float(INS4.kin_of(M, at, cfg)),
            "kin_taper_16_20": float(INS4.kin_of(M, aw, cfg)),
            "shell_edges": edges,
            "shell_rigid": shell_profile(dr, cfg, edges),
            "shell_taper": shell_profile(dt, cfg, edges),
            "taper_sum_check": float(dt.sum()),
        }
        print(f"{tag}: rigid {res[tag]['kin_rigid']:.3f} "
              f"taper {res[tag]['kin_taper_12_15']:.4f}", flush=True)
    # --- tautology test: n48 interior vs n32, same h, commensurate grids ---
    cfg32, M32 = fields["n32_L48_a6"]
    cfg48, M48 = fields["n48_L72_a6"]
    x32 = (np.arange(32) - 15.5) * cfg32["h"]
    x48 = (np.arange(48) - 23.5) * cfg48["h"]
    idx = np.array([int(np.argmin(np.abs(x48 - v))) for v in x32])
    assert np.max(np.abs(x48[idx] - x32)) < 1e-9
    sub = M48[np.ix_(idx, idx, idx)]
    X, Y, Z = INS4.coords(32, cfg32["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    for lab, msk in (("r<12", r < 12.0), ("r<15", r < 15.0),
                     ("all_n32_box", np.ones_like(r, bool))):
        res.setdefault("n48_vs_n32", {})[lab] = {
            "max_abs_dM": float(np.max(np.abs(sub[msk] - M32[msk]))),
            "rms_dM": float(np.sqrt(np.mean((sub[msk] - M32[msk]) ** 2))),
            "ncells": int(msk.sum()),
        }
    # the same comparison on the SEEDS (is the interior identical a priori?)
    s32 = ring_seed(cfg32, 6.0)
    s48 = ring_seed(cfg48, 6.0)[np.ix_(idx, idx, idx)]
    res["n48_vs_n32"]["seed_r<15"] = {
        "max_abs_dM": float(np.max(np.abs(s48[r < 15] - s32[r < 15])))}
    # --- hedgehog tapered comparison ---
    Mh = np.load(os.path.join(CKPT10, "relax_g8_n32_L48_it3000.npy"))
    ah_r = a0_general(Mh)
    ah_t = a0_general(Mh, taper=(12.0, 15.0), cfg=cfg32)
    dh_t = kin_density(Mh, ah_t, cfg32)
    res["hedgehog_it3000"] = {
        "kin_rigid": float(INS4.kin_of(Mh, ah_r, cfg32)),
        "kin_taper_12_15": float(INS4.kin_of(Mh, ah_t, cfg32)),
        "shell_taper": shell_profile(dh_t, cfg32, edges),
        "R10_record_kin_relaxed": 351.16985121697314,
    }
    print("hedgehog:", res["hedgehog_it3000"]["kin_rigid"],
          res["hedgehog_it3000"]["kin_taper_12_15"], flush=True)
    res["runtime_s"] = round(time.time() - t0, 1)
    return res


def stage_c4():
    t0 = time.time()
    cfg = cfg_of(32, 48.0)
    ladder = [0.0, 1.5, 3.0, 4.5, 6.0, 7.5, 9.0, 12.0, 15.0]
    res = {"ladder_a": ladder, "w_c_runs": {}}
    for wc in (2.0, 3.0, 4.0):
        rows = []
        for a in ladder:
            M = ring_seed(cfg, a, w_c=wc)
            eu, ev = st(M, cfg)
            kr = float(INS4.kin_of(M, a0_general(M), cfg))
            kt = float(INS4.kin_of(
                M, a0_general(M, taper=(12.0, 15.0), cfg=cfg), cfg))
            kw = float(INS4.kin_of(
                M, a0_general(M, taper=(16.0, 20.0), cfg=cfg), cfg))
            rows.append({"a": a, "E_u": eu, "V4": ev, "E_stat": eu + ev,
                         "kin_rigid": kr, "kin_tap_12_15": kt,
                         "kin_tap_16_20": kw})
            print(f"wc={wc} a={a:5.1f} E_stat {eu + ev:8.4f} "
                  f"kin_tap {kt:8.3f} kin_w {kw:8.3f} kin_rig {kr:9.3f}",
                  flush=True)
        res["w_c_runs"][f"wc={wc:g}"] = rows
    # --- fixed-J curves (eqn 5) on each kin convention ---
    EJ = {}
    for wc, rows in res["w_c_runs"].items():
        for kkey in ("kin_tap_12_15", "kin_tap_16_20", "kin_rigid"):
            for J in (50.0, 200.0, 800.0):
                vals = [r["E_stat"] + J * J / (4.0 * r[kkey]) for r in rows]
                k = int(np.argmin(vals))
                EJ[f"{wc}|{kkey}|J={J:g}"] = {
                    "E_J": [float(v) for v in vals],
                    "argmin_a": ladder[k],
                    "interior": bool(0 < k < len(ladder) - 1),
                    "omega_star": float(J / (2.0 * rows[k][kkey])),
                }
    res["fixed_J"] = EJ
    res["runtime_s"] = round(time.time() - t0, 1)
    return res


def stage_rerelax():
    """independent re-run of the R10 protocol on MY seed, a = 6, 1500 it."""
    t0 = time.time()
    cfg = cfg_of(32, 48.0)
    M0 = ring_seed(cfg, 6.0)
    free = ~INS4.pin_shell(cfg["n"], cfg["h"], 1.6)
    M, info = INS4.fire(M0, cfg, free, 1500, a0=None, omega=0.0,
                        dt0=0.01, dt_max=0.1, tag="audit_a6")
    eu, ev = st(M, cfg)
    ref = load_ck("n32_L48_a6_it1500.npy", 32)
    os.makedirs(SCRATCH, exist_ok=True)
    np.save(os.path.join(SCRATCH, "mine_a6_it1500.npy"), M)
    out = {"E_u": eu, "V4": ev, "stop": info["stop"],
           "max_abs_dM_vs_producer": float(np.max(np.abs(M - ref))),
           "rms_dM_vs_producer": float(np.sqrt(np.mean((M - ref) ** 2))),
           "max_abs_M": float(np.max(np.abs(M))),
           "cord": cord_locate(M, cfg, 1.5, 10.5),
           "runtime_s": round(time.time() - t0, 1)}
    print(json.dumps(out, indent=1), flush=True)
    return out


def stage_probe():
    """P1 blindness of the grid-argmin cord reader; P2 far field on FREE
    cells only; P3 cord-drift extrapolation."""
    t0 = time.time()
    cfg = cfg_of(32, 48.0)
    # P1: which TRUE cord radii a all read back as the same grid argmin?
    p1 = []
    for a in np.arange(4.0, 7.01, 0.125):
        M = ring_seed(cfg, float(a))
        loc = cord_locate(M, cfg, 1.5, 10.5)
        p1.append({"a_true": float(a),
                   "grid_argmin": loc["gap12_grid_argmin"],
                   "parab": loc["gap12_parab"],
                   "fine": loc["gap12_fine"]})
    blind = [r["a_true"] for r in p1 if r["grid_argmin"] == 5.25]
    # P2: far field on FREE cells only (the pin is exact by construction)
    p2 = {}
    for tag, n, L, a in CASES:
        c = cfg_of(n, L)
        X, Y, Z = INS4.coords(n, c["h"])
        r = np.sqrt(X * X + Y * Y + Z * Z)
        free = ~INS4.pin_shell(n, c["h"], 1.6)
        vv = np.array([0.0, DELTA, 1.0])
        row = {}
        for key in ("seed", "it3000"):
            M = ring_seed(c, a) if key == "seed" else \
                load_ck(f"{tag}_it3000.npy", n)
            for lo in (12.0, 15.0, 18.0):
                msk = free & (r > lo)
                lam = np.linalg.eigvalsh(M[msk][..., 1:4, 1:4])
                row[f"{key}_r>{lo:g}_maxdev"] = float(
                    np.max(np.abs(lam - vv)))
                row[f"{key}_r>{lo:g}_meandev"] = float(
                    np.mean(np.abs(lam - vv)))
        p2[tag] = row
    # P3: geometric extrapolation of the cord drift
    with open(OUT) as f:
        c2 = json.load(f).get("c2", {})
    p3 = {}
    for tag in ("n32_L48_a6", "n32_L48_a9"):
        if tag not in c2:
            continue
        s = c2[tag]["seed"]["cord"]["gap12_fine"]
        a1 = c2[tag]["it1500"]["cord"]["gap12_fine"]
        a2 = c2[tag]["it3000"]["cord"]["gap12_fine"]
        d1v, d2v = a1 - s, a2 - a1
        ratio = d2v / d1v if d1v else float("nan")
        p3[tag] = {"seed": s, "it1500": a1, "it3000": a2,
                   "drift_0_1500": d1v, "drift_1500_3000": d2v,
                   "ratio": float(ratio),
                   "geometric_asymptote": float(
                       a2 + d2v * ratio / (1 - ratio)) if ratio < 1 else None,
                   "pct_shrink_at_3000": float(100 * (a2 - s) / s)}
    out = {"P1_grid_argmin_blindness": {
        "rows": p1, "a_true_reading_5.25": blind,
        "blind_window": [min(blind), max(blind)] if blind else None},
        "P2_far_field_free_cells": p2, "P3_cord_extrapolation": p3,
        "runtime_s": round(time.time() - t0, 1)}
    print(json.dumps({"blind_window": out["P1_grid_argmin_blindness"]
                      ["blind_window"], "P3": p3}, indent=1), flush=True)
    return out


def stage_verdicts():
    """Assemble the verdict table from the stored stages (run LAST)."""
    with open(OUT) as f:
        d = json.load(f)
    c2, c3, c4 = d["c2"], d["c3"], d["c4"]
    pr = d["probe"]
    wc3 = {r["a"]: r for r in c4["w_c_runs"]["wc=3"]}
    astar = {w: c4["fixed_J"][f"{w}|kin_tap_12_15|J=50"]["argmin_a"]
             for w in ("wc=2", "wc=3", "wc=4")}
    estat_argmin = {w: min(c4["w_c_runs"][w], key=lambda r: r["E_stat"])["a"]
                    for w in ("wc=2", "wc=3", "wc=4")}
    v = {
        "C1_instrument": {
            "verdict": "CONFIRMED",
            "note": "my eigenvector winding reader on radii "
                    "1.875/2.625/3.375/4.5 gives q = 0.5000 on both seeds; "
                    "my sphere-Jacobian degree at R = 15 gives +0.9999 "
                    "(hedgehog self-test +0.9999); the record cube reader "
                    "gives +1.0000 with 0 lift conflicts; a0 identity rel "
                    "6.667e-09, kin 426.50701 vs 426.50702.",
            "my_numbers": {
                "a0_rel": d["c1"]["a0_identity"]["rel"],
                "kin_mine": d["c1"]["a0_identity"]["kin_mine"],
                "deg_a6": d["c1"]["seed_a6"]["degree_sphere_R15"],
                "deg_a9": d["c1"]["seed_a9"]["degree_sphere_R15"]}},
        "C2_statics": {
            "verdict": "REFUTED on 'does not shrink'; energies CONFIRMED",
            "refuting_number":
                "sub-grid cord radius (0.01 resolution, three locators) "
                "a = 6: 6.000 -> 5.760 -> 5.610 (-6.50 % at it3000, still "
                "moving -0.15 per 1500 it); a = 9: 9.000 -> 8.820 -> 8.710 "
                "(-3.22 %). The producer's grid argmin reads 5.25 for ANY "
                "true cord in [4.5, 6.0] (probe P1), so it cannot resolve a "
                "shrink of up to 25 %.",
            "energies_reproduced": True,
            "hedgehog_lead_is_a_stopping_artifact":
                "ring a6 E_u 11.060 at it3000 vs hedgehog 13.540 at it3000, "
                "but the same hedgehog reaches 11.035 at it6000 and 9.047 at "
                "it12000; both states are unconverged at 3000.",
            "cord": pr["P3_cord_extrapolation"],
            "hedgehog": c2["hedgehog_R10"]},
        "C3_clock": {
            "verdict": "QUALIFIED",
            "tapered_equality_is_a_tautology":
                "n48 and n32 relaxed interiors agree to max|dM| = 6.7e-07 on "
                "r < 15 (seeds identical to 0.0 there), so the tapered kin "
                "agreeing to 1e-06 is the same field measured twice.",
            "does_not_localize_on_the_cord":
                "tapered kin shells (a = 6): 0.37 / 7.55 / 22.76 / 36.27 / "
                "12.16 over r = 0-3..12-15; the peak is the 9-12 shell, not "
                "the cord at r = 6, and the hedgehog profile has the same "
                "shape (1.38 / 11.94 / 27.36 / 37.57 / 12.48).",
            "cutoff_dependence":
                "taper (16, 20) instead of (12, 15) gives 136.35 not 79.11, "
                "so the tapered kin is a cutoff number, not a converged one.",
            "hedgehog_reference_mismatch":
                "the record's 115.5 uses the ANALYTIC a0_unit on the relaxed "
                "hedgehog; with the SAME generalized a0 the ring must use, "
                "the hedgehog gives rigid 319.13 / tapered 90.74, so the "
                "ring's tapered deficit is 79.11/90.74 = -12.8 %, not -31 %.",
            "numbers": {k: {kk: c3[k][kk] for kk in c3[k]
                            if kk.startswith("kin")}
                        for k in ("n32_L48_a6", "n32_L48_a9", "n48_L72_a6",
                                  "hedgehog_it3000")},
            "n48_vs_n32": c3["n48_vs_n32"]},
        "C4_fixed_J": {
            "verdict": "QUALIFIED (arithmetic CONFIRMED, interpretation "
                       "REFUTED)",
            "ladder_reproduced": {a: wc3[a]["E_stat"] for a in
                                  (3.0, 4.5, 6.0, 9.0, 12.0)},
            "w_c_moves_the_minimum":
                f"E_stat argmin over a moves with the seed's cord width: "
                f"{estat_argmin}; a*(J=50) moves {astar}. The interior "
                f"minimum is a property of the fixed w_c = 3, not of the "
                f"functional (a* ~ 2.2 w_c).",
            "J200_interior_min_exists":
                "extending the ladder below a = 3 (to a = 1.5 and the exact "
                "hedgehog limit a = 0) gives argmin a = 1.5 at J = 200 "
                "(E_J 108.47 vs 112.35 at a = 3 and 117.96 at a = 0), so "
                "'no interior minimum for J = 200' is a ladder-truncation "
                "artifact; only J = 800 sits on the a = 0 endpoint.",
            "omega_star_is_cutoff_dependent":
                "omega* at the minimum: 0.268 (taper 12-15), 0.176 (taper "
                "16-20), 0.083 (rigid). kin is IR-divergent (flat ~40 per "
                "3-unit shell to the box wall), so omega*(J fixed) -> 0 as "
                "L grows and the radiation-window test is not converged.",
            "a0_hedgehog_limit_check":
                "a = 0 IS the radial hedgehog seed; its E_stat 44.26 is the "
                "highest on the ladder, so the E_stat interior minimum does "
                "survive the extension (this half of C4 stands)."},
        "reproduction_of_the_producer_protocol": {
            "verdict": "BIT-IDENTICAL",
            "note": "my independently built a = 6 seed + the R10 protocol "
                    "(FIRE, pin 1.6, dt0 0.01, dt_max 0.1, 1500 it) "
                    "reproduces the producer's n32_L48_a6_it1500.npy with "
                    "max|dM| = 0.0 exactly.",
            "max_abs_dM": d["rerelax"]["max_abs_dM_vs_producer"]},
        "far_field_qualifier":
            "the far field is NOT the exact vacuum spectrum even in the "
            "seed: max |lambda - (0, 0.3, 1)| = 0.2754 on free cells at "
            "r > 18 (the e_phi axis melt inside rho < rho_c runs the whole "
            "z axis). Relaxation leaves the max at 0.2667 but raises the "
            "MEAN deviation 1.46e-04 -> 5.07e-04. The pin cells themselves "
            "are exact (max dev 0.0).",
    }
    return v


def merge(stage, payload):
    old = {}
    if os.path.exists(OUT):
        with open(OUT) as f:
            old = json.load(f)
    old[stage] = payload
    old["_meta"] = {"script": "m5_32_r12_audit.py",
                    "python": sys.version.split()[0],
                    "numpy": np.__version__}
    with open(OUT, "w") as f:
        json.dump(old, f, indent=1)
    print(f"wrote {OUT} [{stage}]", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "c1"
    fn = {"c1": stage_c1, "c2": stage_c2, "c3": stage_c3,
          "c4": stage_c4, "rerelax": stage_rerelax,
          "probe": stage_probe, "verdicts": stage_verdicts}[mode]
    merge(mode, fn())
