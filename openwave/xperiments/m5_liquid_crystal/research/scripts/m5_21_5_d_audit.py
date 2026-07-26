"""M5.21.5 D: INDEPENDENT ADVERSARIAL AUDIT of the mu/S/g claims
(C1-C6) of scripts/m5_21_5_a_mu.py + data/m5_21_5_bridge.json.

Independence contract: every physics quantity recomputed here uses
OWN implementations (own eigenframe gauge code, own Rodrigues, own
stencils, own mu integrator, own tangent construction, own EID seed
build). The ONLY reuse from the primary instrument stack is
m5_21_3_a_4d.kin_of / e_parts / base_cfg / embed34, permitted as the
certified functional of record; the tangents fed to kin_of are built
here and cross-checked by finite-differencing the rotation flow.
Deliberate probe deviations from the record pipeline: central (not
forward) dtn stencil, dphi = 3e-4 (not 1e-4), dphi2 variant 2e-3,
pole-fallback threshold variant, 1e-6 relative noise injection.

Out: data/m5_21_5_audit.json
"""
from __future__ import annotations

import importlib.util
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_s4 = importlib.util.spec_from_file_location(
    "ins4", os.path.join(HERE, "m5_21_3_a_4d.py"))
INS4 = importlib.util.module_from_spec(_s4)
_s4.loader.exec_module(INS4)

ETA4 = np.diag([-1.0, 1.0, 1.0, 1.0])


# ================= own geometry =====================================
def grid(n, h):
    c = (n - 1) / 2.0
    ax = (np.arange(n) - c) * h
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    pos = np.stack([X, Y, Z], axis=-1)
    return pos, np.sqrt(X * X + Y * Y + Z * Z)


def dcen(f, ax, h):
    return (np.roll(f, -1, axis=ax) - np.roll(f, 1, axis=ax)) / (2 * h)


def sp(M):
    return M if M.shape[-1] == 3 else M[..., 1:4, 1:4]


def env_of(r, renv=10.0):
    return np.exp(-((r / renv) ** 4))


# ================= own eigenframe (same definition, own code) =======
def frame(Msp, pos, pole_tol=0.1):
    lam, V = np.linalg.eigh(Msp)
    v1, v2 = V[..., :, 2].copy(), V[..., :, 1].copy()
    r = np.linalg.norm(pos, axis=-1, keepdims=True)
    rhat = pos / np.maximum(r, 1e-12)
    s1 = np.sign(np.einsum("...a,...a->...", v1, rhat))
    s1[s1 == 0] = 1.0
    v1 *= s1[..., None]
    # polar reference p = normalize((v1.z) v1 - z); algebraically the
    # same direction as the record's normalize(z x v1) x v1
    z = np.array([0.0, 0.0, 1.0])
    p = v1 * v1[..., 2:3] - z
    pn = np.linalg.norm(p, axis=-1)
    pole = pn < pole_tol                # |z x v1| = |p| identically
    if pole.any():
        x = np.array([1.0, 0.0, 0.0])
        p2 = v1 * v1[..., 0:1] - x
        p[pole] = p2[pole]
    p /= np.maximum(np.linalg.norm(p, axis=-1, keepdims=True), 1e-300)
    s2 = np.sign(np.einsum("...a,...a->...", v2, p))
    s2[s2 == 0] = 1.0
    v2 *= s2[..., None]
    v3 = np.cross(v1, v2)
    v3 /= np.maximum(np.linalg.norm(v3, axis=-1, keepdims=True),
                     1e-300)
    return v1, v2, v3, lam[..., 2] - lam[..., 1], lam[..., 1] - lam[..., 0]


def rot_about(axis, ang):
    """Per-voxel Rodrigues, own assembly via the rotation formula
    v' = v cos + (a x v) sin + a (a.v)(1 - cos), built as a matrix."""
    c, s = np.cos(ang), np.sin(ang)
    a = axis
    R = np.empty(a.shape[:-1] + (3, 3))
    for i in range(3):
        for j in range(3):
            R[..., i, j] = a[..., i] * a[..., j] * (1 - c)
        R[..., i, i] += c
    R[..., 0, 1] -= s * a[..., 2]
    R[..., 1, 0] += s * a[..., 2]
    R[..., 0, 2] += s * a[..., 1]
    R[..., 2, 0] -= s * a[..., 1]
    R[..., 1, 2] -= s * a[..., 0]
    R[..., 2, 1] += s * a[..., 0]
    return R


def rot_state(M, axis, ang):
    R = rot_about(axis, ang)
    if M.shape[-1] == 3:
        return np.einsum("...ab,...bc,...dc->...ad", R, M, R)
    R4 = np.zeros(M.shape[:-2] + (4, 4))
    R4[..., 0, 0] = 1.0
    R4[..., 1:4, 1:4] = R
    return np.einsum("...ab,...bc,...dc->...ad", R4, M, R4)


def align(nn, nref):
    s = np.sign(np.einsum("...a,...a->...", nn, nref))
    s[s == 0] = 1.0
    return nn * s[..., None]


def mermin(nh, h):
    dn = [dcen(nh, ax, h) for ax in range(3)]
    F12 = np.einsum("...a,...a->...", nh, np.cross(dn[0], dn[1]))
    F13 = np.einsum("...a,...a->...", nh, np.cross(dn[0], dn[2]))
    F23 = np.einsum("...a,...a->...", nh, np.cross(dn[1], dn[2]))
    return np.stack([F23, -F13, F12], axis=-1), dn


def efield(nh, dtn, dn):
    return np.stack([np.einsum("...a,...a->...", nh,
                               np.cross(dtn, dn[ax]))
                     for ax in range(3)], axis=-1)


def curl3(V, h):
    return np.stack(
        [dcen(V[..., 2], 1, h) - dcen(V[..., 1], 2, h),
         dcen(V[..., 0], 2, h) - dcen(V[..., 2], 0, h),
         dcen(V[..., 1], 0, h) - dcen(V[..., 0], 1, h)], axis=-1)


# ================= own mu pipeline ==================================
def mu_audit(M, h, channel, mask, w=None, dphi=3e-4, dphi2=1e-3,
             pole_tol=0.1, rwin=(4, 6, 8, 10, 12, 14, 16, 100)):
    """Own read of the SAME definition. Probe deviations: central
    dtn stencil, dphi = 3e-4 default. Returns diagnostics fields."""
    n = M.shape[0]
    pos, r = grid(n, h)
    v1, v2, v3, g12, g23 = frame(sp(M), pos, pole_tol)
    ax = {"twist": v1, "tilt12": v3, "tilt13": -v2}.get(channel)
    if ax is None and channel == "globalz":
        ax = np.zeros_like(v1)
        ax[..., 2] = 1.0
    wf = 1.0 if w is None else w

    def nof(phi, ref):
        if phi == 0.0:
            return v1
        Mp = rot_state(M, ax, phi * wf)
        nb, _, _, _, _ = frame(sp(Mp), pos, pole_tol)
        return align(nb, ref)

    n0 = v1
    B0, dn0 = mermin(n0, h)
    dtn0 = (nof(dphi, n0) - nof(-dphi, n0)) / (2 * dphi)
    E0 = efield(n0, dtn0, dn0)
    n1 = nof(dphi2, n0)
    B1, dn1 = mermin(n1, h)
    dtn1 = (align(nof(dphi2 + dphi, n0), n1)
            - align(nof(dphi2 - dphi, n0), n1)) / (2 * dphi)
    E1 = efield(n1, dtn1, dn1)
    jcb = curl3(B0, h)
    jck = -(E1 - E0) / dphi2
    cck = 0.5 * np.cross(pos, jck) * h ** 3      # clock integrand
    ccb = 0.5 * np.cross(pos, jcb) * h ** 3
    rows = {}
    for R in rwin:
        m = mask & (r < R)
        mu_ck = cck[m].sum(axis=0)
        mu_cb = ccb[m].sum(axis=0)
        rows[f"R{R}"] = {
            "mu_clock": float(np.linalg.norm(mu_ck)),
            "mu_curlB": float(np.linalg.norm(mu_cb)),
            "mu_total": float(np.linalg.norm(mu_ck + mu_cb)),
            "vox": int(m.sum())}
    mu_ck = cck[mask].sum(axis=0)
    mu_cb = ccb[mask].sum(axis=0)
    return {"channel": channel,
            "mu_clock": float(np.linalg.norm(mu_ck)),
            "mu_clock_vec": [float(x) for x in mu_ck],
            "mu_curlB": float(np.linalg.norm(mu_cb)),
            "mu_total": float(np.linalg.norm(mu_ck + mu_cb)),
            "radial": rows,
            "_cck": cck, "_mask": mask, "_r": r, "_pos": pos,
            "_v1": v1}


def masks(M, h, margin=2, gap_min=0.02, pole_tol=0.1):
    n = M.shape[0]
    pos, _ = grid(n, h)
    _, _, _, g12, g23 = frame(sp(M), pos, pole_tol)
    inter = np.zeros((n, n, n), bool)
    inter[margin:-margin, margin:-margin, margin:-margin] = True
    return inter & (g12 > gap_min), inter & (g12 > gap_min) & (
        g23 > gap_min), inter


# ================= own tangent + S ==================================
def gen4_of(axis):
    G = np.zeros(axis.shape[:-1] + (4, 4))
    a = axis
    G[..., 1, 2], G[..., 1, 3] = -a[..., 2], a[..., 1]
    G[..., 2, 1], G[..., 2, 3] = a[..., 2], -a[..., 0]
    G[..., 3, 1], G[..., 3, 2] = -a[..., 1], a[..., 0]
    return G


def tangent(M4, axis):
    G = gen4_of(axis)
    return G @ M4 - M4 @ G


def s_audit(M4, cfg4, channel, pole_tol=0.1):
    n = M4.shape[0]
    pos, r = grid(n, cfg4["h"])
    v1, v2, v3, _, _ = frame(sp(M4), pos, pole_tol)
    ax = {"twist": v1, "tilt12": v3, "tilt13": -v2}[channel]
    tg = tangent(M4, ax)
    # finite-difference gate on the tangent construction
    e = 1e-5
    tg_fd = (rot_state(M4, ax, e) - rot_state(M4, ax, -e)) / (2 * e)
    fd_err = float(np.abs(tg_fd - tg).max() / max(np.abs(tg).max(),
                                                  1e-300))
    w = env_of(r)[..., None, None]
    kin_r = float(INS4.kin_of(M4, tg, cfg4))
    kin_e = float(INS4.kin_of(M4, w * tg, cfg4))
    return {"S_rigid": 2 * kin_r, "S_env": 2 * kin_e,
            "tangent_fd_relerr": fd_err}


def kin_density(M4, a0, cfg4):
    """Own per-cell kin integrand (sym stencil), to locate where S
    lives; total must match INS4.kin_of."""
    h = cfg4["h"]
    dens = np.zeros(M4.shape[:3])
    for st, wt in (("fwd", 0.5), ("bwd", 0.5)):
        for axd in range(3):
            if st == "fwd":
                A = (np.roll(M4, -1, axis=axd) - M4) / h
                sl = [slice(None)] * 3
                sl[axd] = -1
                A[tuple(sl)] = 0.0
            else:
                A = (M4 - np.roll(M4, 1, axis=axd)) / h
                sl = [slice(None)] * 3
                sl[axd] = 0
                A[tuple(sl)] = 0.0
            F = a0 @ ETA4 @ A - A @ ETA4 @ a0
            dens += wt * 4.0 * np.einsum(
                "...ab,...cd,ac,bd->...", F, F, ETA4, ETA4,
                optimize=True)
    return dens * h ** 3


# ================= C1: the algebra ==================================
def c1_algebra():
    out = {}
    rng = np.random.default_rng(7)
    # hedgehog closed form from exact analytic derivatives
    x = rng.normal(size=(4000, 3)) * 3 + 0.5
    r = np.linalg.norm(x, axis=-1, keepdims=True)
    nh = x / r
    Ident = np.eye(3)
    dn = (Ident[None] - nh[..., :, None] * nh[..., None, :]) / r[
        ..., None]                       # dn[p, i, j] = d_i n_j
    errF = errB = 0.0
    F = np.zeros((x.shape[0], 3, 3))
    for i in range(3):
        for j in range(3):
            F[:, i, j] = np.einsum(
                "pa,pa->p", nh, np.cross(dn[:, i, :], dn[:, j, :]))
    eps = np.zeros((3, 3, 3))
    eps[0, 1, 2] = eps[1, 2, 0] = eps[2, 0, 1] = 1
    eps[0, 2, 1] = eps[2, 1, 0] = eps[1, 0, 2] = -1
    Fth = np.einsum("ijk,pk->pij", eps, nh) / r[..., None] ** 2
    errF = float(np.abs(F - Fth).max() / np.abs(Fth).max())
    B = 0.5 * np.einsum("kij,pij->pk", eps, F)
    errB = float(np.abs(np.sum(B * B, axis=-1) * r[:, 0] ** 4
                        - 1).max())
    out["hedgehog"] = {
        "F_eq_eps_n_over_r2_relerr": errF,
        "B2_r4_minus_1_max": errB,
        "flux_exact": "B_r = 1/r^2 so the sphere integral = 4 pi "
                      "identically",
        "sum_Fij2_ilej_r4": float(np.abs(
            0.5 * np.sum(F * F, axis=(1, 2)) * r[:, 0] ** 4
            - 1).max())}
    # numeric flux via lattice divergence (own stencil, ball sum)
    n, h = 48, 0.5
    pos, r3 = grid(n, h)
    nh3 = pos / np.maximum(r3[..., None], 1e-12)
    Bl = nh3 / np.maximum(r3[..., None], 1e-12) ** 2
    divB = sum(dcen(Bl[..., a], a, h) for a in range(3))
    annulus = (r3 > 0.9) & (r3 < 8)      # excludes the singular core
    # monopole: flux(R) - flux(0.9) = 0, so the annulus div-sum ~ 0
    out["hedgehog"]["annulus_divB_sum_over_4pi"] = float(
        (divB[annulus].sum() * h ** 3) / (4 * np.pi))
    shell = (r3 > 6.0) & (r3 < 6.0 + h)
    Br = np.einsum("...a,...a->...", Bl, nh3)
    # exact surface law: Br * r^2 = 1 identically -> flux = 4 pi
    out["hedgehog"]["Br_r2_mean_shell"] = float(
        (Br[shell] * r3[shell] ** 2).mean())
    # unit chain: Gaussian cgs numbers
    alpha = 7.2973525693e-3
    hbar, c = 1.054571817e-27, 2.99792458e10
    me = 9.1093837015e-28
    e_g = np.sqrt(alpha * hbar * c)
    rows = []
    for (c2p, ell) in ((0.211, 3.7e-12), (0.05, 1.9e-11)):
        C_E = alpha * hbar * c / (64 * np.pi * c2p * ell)
        # far field check: u_phys(r) = 8 c2p C_E ell / r^4 target
        u_coef = 8 * c2p * C_E * ell
        q = 8 * np.sqrt(np.pi * c2p * C_E * ell)
        lam = 8 * np.sqrt(np.pi * c2p * C_E / ell ** 3)
        mu_coef = lam * ell ** 3 / (4 * np.pi)   # mu_phys per mu_lat
        # independent rederivation of mu_coef from Maxwell factors:
        # Gaussian Ampere J = (c/4pi)(curl B - dE/dt / c); with
        # B_phys = lam B_lat, curl_phys = curl_lat/ell, t_phys =
        # (ell/c) t_lat (c_lat = 1): J_phys = (c lam/(4pi ell)) j_lat
        # mu = (1/2c) int r x J dV = (1/2c)(c lam/(4pi ell))
        #      * (ell)(ell^3) * [int r_lat x j_lat dV_lat]
        #    = lam ell^3/(8pi) * (2 mu_lat)
        mu_coef_ind = (1 / (2 * c)) * (c * lam / (4 * np.pi * ell)) \
            * ell * ell ** 3 * 2.0
        mu_lat, S_lat, = 3.3776756, 37.1822136
        E_lat = me * c ** 2 / C_E        # the mass anchor
        mu_phys = mu_coef * mu_lat
        S_phys = C_E * (ell / c) * S_lat
        g_dim = 2 * me * c * mu_phys / (e_g * S_phys)
        g_formula = mu_lat * E_lat / (2 * np.pi * S_lat)
        k = g_dim * S_lat / (2 * mu_lat * E_lat)
        rows.append({
            "c2p": c2p, "l_cm": ell,
            "u_coef_over_alpha_hbar_c_8pi": float(
                u_coef / (alpha * hbar * c / (8 * np.pi))),
            "q_over_e": float(q / e_g),
            "lam_l2_over_q": float(lam * ell ** 2 / q),
            "mu_coef_indep_ratio": float(mu_coef_ind / mu_coef),
            "g_dimensional": float(g_dim),
            "g_formula_muE_2piS": float(g_formula),
            "g_ratio": float(g_dim / g_formula),
            "k_derived": float(k),
            "k_target_1_over_4pi": float(1 / (4 * np.pi))})
    out["unit_chain"] = rows
    # rigid rotor S = dE/domega at omega = 1 equals 2 kin
    cfg = INS4.base_cfg(s=-1.0, n=8, L=12.0)
    rng2 = np.random.default_rng(3)
    Mt = rng2.normal(size=(8, 8, 8, 4, 4))
    Mt = 0.5 * (Mt + np.swapaxes(Mt, -1, -2))
    at = rng2.normal(size=(8, 8, 8, 4, 4))
    at = 0.5 * (at + np.swapaxes(at, -1, -2))
    k1 = float(INS4.kin_of(Mt, at, cfg))
    k2 = float(INS4.kin_of(Mt, 2 * at, cfg))
    ee = 1e-4
    dEdo = (float(INS4.kin_of(Mt, (1 + ee) * at, cfg))
            - float(INS4.kin_of(Mt, (1 - ee) * at, cfg))) / (2 * ee)
    out["rigid_rotor"] = {"kin_2a_over_4kin_a": float(k2 / (4 * k1)),
                          "dE_domega_at_1_over_2kin": float(
                              dEdo / (2 * k1))}
    return out


# ================= EID seed (own build) =============================
def eid_seed(n=24, box=6.0, delta=0.3, rc=0.8, rhoc=0.8):
    xs = np.linspace(-box, box, n)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rho = np.sqrt(X * X + Y * Y)
    rhat = np.stack([X, Y, Z], -1) / np.maximum(r[..., None], 1e-15)
    azim = np.stack([-Y, X, np.zeros_like(X)], -1) / np.sqrt(
        rho ** 2 + 1e-12)[..., None]
    s = np.clip(rho / rhoc, 0.0, 1.0)
    shrink = (s * s * (3.0 - 2.0 * s))[..., None]
    ephi = azim * shrink
    ephi = ephi - np.einsum("...i,...i->...", ephi,
                            rhat)[..., None] * rhat
    etheta = np.cross(ephi, rhat)
    M3 = (rhat[..., :, None] * rhat[..., None, :]
          + delta * etheta[..., :, None] * etheta[..., None, :])
    inter = np.zeros((n, n, n), bool)
    inter[2:-2, 2:-2, 2:-2] = True
    act = inter & (r > 2 * rc) & (rho > rhoc)
    return M3, h, act, inter, r, rho


# ================= diagnostics ======================================
def structure(res, label):
    cck, mask, r, pos = res["_cck"], res["_mask"], res["_r"], res[
        "_pos"]
    mu = np.array(res["mu_clock_vec"])
    mun = np.linalg.norm(mu)
    mhat = mu / max(mun, 1e-300)
    proj = cck @ mhat
    pm = np.where(mask, proj, 0.0)
    net = float(pm.sum())
    pplus = float(pm[pm > 0].sum())
    pminus = float(-pm[pm < 0].sum())
    edges = [0, 2, 4, 6, 8, 10, 12, 14, 16, 1e9]
    shells = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = mask & (r >= a) & (r < b)
        shells.append({"r": [a, min(b, 999)],
                       "dmu": float(pm[m].sum()),
                       "vox": int(m.sum())})
    flat = np.abs(pm[mask])
    order = np.argsort(flat)[::-1]
    tot = flat.sum()
    cum = np.cumsum(flat[order])
    n50 = int(np.searchsorted(cum, 0.5 * tot) + 1)
    sgn = pm[mask][order]
    top10_net = float(sgn[:10].sum())
    top1pc = int(max(1, round(0.01 * mask.sum())))
    top1pc_net = float(sgn[:top1pc].sum())
    hemi = {}
    for axn, axi in (("x", 0), ("y", 1), ("z", 2)):
        c = pos[..., axi]
        hemi[axn] = [float(pm[mask & (c > 0)].sum()),
                     float(pm[mask & (c < 0)].sum())]
    octs = []
    for sx in (1, -1):
        for sy in (1, -1):
            for sz in (1, -1):
                m = mask & (sx * pos[..., 0] > 0) & (
                    sy * pos[..., 1] > 0) & (sz * pos[..., 2] > 0)
                octs.append(float(pm[m].sum()))
    return {"label": label, "net": net, "gross": pplus + pminus,
            "cancel_net_over_gross": float(
                net / max(pplus + pminus, 1e-300)),
            "pos_sum": pplus, "neg_sum": pminus,
            "mu_hat": [float(x) for x in mhat],
            "angle_to_z_deg": float(np.degrees(np.arccos(
                np.clip(abs(mhat[2]), 0, 1)))),
            "shells": shells,
            "hemispheres_pos_neg": hemi, "octants": octs,
            "vox_for_50pct_of_gross": n50,
            "n_masked": int(mask.sum()),
            "top10_net_over_net": float(top10_net / net),
            "top1pct_net_over_net": float(top1pc_net / net)}


# ================= main =============================================
def main():
    out = {}
    print("== C1 algebra ==", flush=True)
    out["C1"] = c1_algebra()
    print(json.dumps(out["C1"]["unit_chain"][0], indent=1)[:600],
          flush=True)

    # ---- load states
    def load(p):
        z = np.load(os.path.join(DATA, p))
        return z["M"].astype(np.float64), float(
            z["h"]) if "h" in z else 1.5

    Mfj, hfj = load("m5_21_9_fixedj_conj_om0.2_end.npz")
    Mt32, ht32 = load("m5_21_6_end_t32_A.npz")
    Mt24, ht24 = load("m5_21_5_end_t24_A.npz")

    # ---- C2: fjom0.2 envelope-localized tilt12 clock moment
    print("== C2 fjom0.2 ==", flush=True)
    n = Mfj.shape[0]
    _, r = grid(n, hfj)
    w = env_of(r)
    mdir, mfr, minter = masks(Mfj, hfj)
    base = mu_audit(Mfj, hfj, "tilt12", mfr, w=w)
    var_d2 = mu_audit(Mfj, hfj, "tilt12", mfr, w=w, dphi2=2e-3)
    var_pole = mu_audit(Mfj, hfj, "tilt12", mfr, w=w, pole_tol=0.05)
    rng = np.random.default_rng(11)
    scale = np.abs(Mfj).max()
    P = rng.normal(size=Mfj.shape) * 1e-6 * scale
    P = 0.5 * (P + np.swapaxes(P, -1, -2))
    noisy = mu_audit(Mfj + P, hfj, "tilt12",
                     masks(Mfj + P, hfj)[1], w=w)
    twist = mu_audit(Mfj, hfj, "twist", mdir, w=w)
    out["C2"] = {
        "record": 3.3776756140804354,
        "mine_dphi3e-4_central": base["mu_clock"],
        "mine_vec": base["mu_clock_vec"],
        "var_dphi2_2e-3": var_d2["mu_clock"],
        "var_pole_tol_0.05": var_pole["mu_clock"],
        "noise_1e-6_rel": noisy["mu_clock"],
        "radial": {k: v["mu_clock"]
                   for k, v in base["radial"].items()},
        "twist_mu_clock": twist["mu_clock"],
        "mu_curlB": base["mu_curlB"], "mu_total": base["mu_total"]}
    print(" mu_clock mine", base["mu_clock"], "record 3.3777",
          flush=True)

    # ---- C3: t32A suppression + where the moments live
    print("== C3 t32A ==", flush=True)
    _, r32 = grid(Mt32.shape[0], ht32)
    w32 = env_of(r32)
    mdir32, mfr32, _ = masks(Mt32, ht32)
    t32 = mu_audit(Mt32, ht32, "tilt12", mfr32, w=w32)
    out["C3"] = {
        "record": 0.018245417515429206,
        "mine": t32["mu_clock"],
        "mine_vec": t32["mu_clock_vec"],
        "ratio_fj_over_t32": base["mu_clock"] / t32["mu_clock"],
        "radial": {k: v["mu_clock"]
                   for k, v in t32["radial"].items()}}
    out["structure_fjom"] = structure(base, "fjom0.2 tilt12 clock")
    out["structure_t32A"] = structure(t32, "t32A tilt12 clock")
    # direction vs the state's own structure axes (fjom)
    t13 = mu_audit(Mfj, hfj, "tilt13", mfr, w=w)
    v1fj = base["_v1"]
    M2 = np.einsum("pa,pb->ab", v1fj[mfr], v1fj[mfr]) / mfr.sum()
    lam2, V2 = np.linalg.eigh(M2)
    mh12 = np.array(base["mu_clock_vec"])
    mh12 /= np.linalg.norm(mh12)
    mh13 = np.array(t13["mu_clock_vec"])
    mh13 /= np.linalg.norm(mh13)
    axes_ang = {f"eig{i}_lam{lam2[i]:.3f}": float(np.degrees(
        np.arccos(np.clip(abs(mh12 @ V2[:, i]), 0, 1))))
        for i in range(3)}
    om_dirs = {}
    for tg in ("fjom0.5", "fjom1"):
        try:
            dj = json.load(open(os.path.join(
                DATA, f"m5_21_5_mu_{tg}.json")))
            vv = np.array(dj["mu_tilt12"]["radial"]["R100"][
                "mu_clock_vec"])
            om_dirs[tg] = float(np.degrees(np.arccos(np.clip(
                abs(vv @ mh12) / np.linalg.norm(vv), 0, 1))))
        except FileNotFoundError:
            pass
    out["structure_direction"] = {
        "tilt12_vs_tilt13_deg": float(np.degrees(np.arccos(
            np.clip(abs(mh12 @ mh13), 0, 1)))),
        "tilt12_mu_vs_nematic_axes_deg": axes_ang,
        "nematic_eigvals": [float(x) for x in lam2],
        "angle_to_fjom0.2_dir_deg": om_dirs,
        "tilt12_angle_to_z_deg": out["structure_fjom"][
            "angle_to_z_deg"]}
    print(" mu_clock mine", t32["mu_clock"], "ratio",
          out["C3"]["ratio_fj_over_t32"], flush=True)

    # ---- C4: analytic EID seed, own build, rigid tilt12
    print("== C4 EID seed ==", flush=True)
    M3, he, act, inter, re_, rho = eid_seed()
    # NOTE own grid centers == linspace centers (both symmetric)
    _, mfr_e, _ = masks(M3, he)
    sphere = inter & (re_ > 1.6)         # sphere cut, NO cylinder
    rw = (2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10, 100)
    a_act = mu_audit(M3, he, "tilt12", act, w=None, rwin=rw)
    a_sph = mu_audit(M3, he, "tilt12", sphere, w=None, rwin=rw)
    a_gap = mu_audit(M3, he, "tilt12", mfr_e, w=None, rwin=rw)
    a_raw = mu_audit(M3, he, "tilt12", inter, w=None, rwin=rw)
    # power-law fit of mu(<R) over COMPLETE shells (R <= 6 = box
    # half-width): R^3 growth means box-truncated, not converged
    Rs = np.array([2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6])
    mus = np.array([a_act["radial"][f"R{g}" if g != int(g)
                                    else f"R{int(g)}"]["mu_clock"]
                    for g in Rs])
    slope = float(np.polyfit(np.log(Rs), np.log(mus), 1)[0])
    out["C4"] = {
        "record": 0.2209, "h": he,
        "mine_act_mask": a_act["mu_clock"],
        "mine_act_mu_total": a_act["mu_total"],
        "mine_sphere_cut_only": a_sph["mu_clock"],
        "mine_gap_mask": a_gap["mu_clock"],
        "mine_full_interior_nogap": a_raw["mu_clock"],
        "radial_act": {k: v["mu_clock"]
                       for k, v in a_act["radial"].items()},
        "mu_growth_exponent_R2.5_to_6": slope,
        "frac_from_corner_region_r_gt_6": float(
            1 - a_act["radial"]["R6"]["mu_clock"]
            / a_act["mu_clock"]),
        "note": "own seed build: exact r_hat (the record's rc "
                "regularization is nullified by its own "
                "renormalization), same ePhi shrink + Gram-Schmidt"}
    print(" act", a_act["mu_clock"], "sphere", a_sph["mu_clock"],
          "raw", a_raw["mu_clock"], "growth exp", slope, flush=True)

    # ---- C5: S values, own tangents
    print("== C5 S ==", flush=True)
    cfg32 = INS4.base_cfg(s=-1.0, n=32, L=48.0)
    cfg24 = INS4.base_cfg(s=-1.0, n=24, L=36.0)
    M4t32 = INS4.embed34(Mt32, cfg32)
    M4t24 = INS4.embed34(Mt24, cfg24)
    s_fj_tw = s_audit(Mfj, cfg32, "twist")
    s_fj_t12 = s_audit(Mfj, cfg32, "tilt12")
    s_t32_tw = s_audit(M4t32, cfg32, "twist")
    s_t32_t12 = s_audit(M4t32, cfg32, "tilt12")
    s_t24_tw = s_audit(M4t24, cfg24, "twist")
    s_t24_t12 = s_audit(M4t24, cfg24, "tilt12")
    # vacuum commutator: [G_z, diag(8, 1, 0.3, 0)]
    Dv = np.diag([8.0, 1.0, 0.3, 0.0])
    Gz = gen4_of(np.array([0.0, 0.0, 1.0]))
    comm = Gz @ Dv - Dv @ Gz
    # IR growth: synthetic hedgehog-texture vacuum at 3 box sizes
    grow = {}
    for nn in (24, 32, 48):
        cfgN = INS4.base_cfg(s=-1.0, n=nn, L=nn * 1.5)
        posN, rN = grid(nn, 1.5)
        rhatN = posN / np.maximum(rN[..., None], 1e-12)
        aziN = np.stack([-posN[..., 1], posN[..., 0],
                         np.zeros_like(rN)], -1)
        rhoN = np.maximum(np.linalg.norm(aziN, axis=-1,
                                         keepdims=True), 1e-12)
        aziN = aziN / rhoN
        sN = np.clip(rhoN[..., 0] / 0.8, 0, 1)
        aziN *= (sN * sN * (3 - 2 * sN))[..., None]
        aziN -= np.einsum("...i,...i->...", aziN,
                          rhatN)[..., None] * rhatN
        thN = np.cross(aziN, rhatN)
        M3N = (rhatN[..., :, None] * rhatN[..., None, :]
               + 0.3 * thN[..., :, None] * thN[..., None, :])
        M4N = INS4.embed34(M3N, cfgN)
        v1N, _, _, _, _ = frame(M3N, posN)
        tgN = tangent(M4N, v1N)
        grow[f"n{nn}"] = {"L": nn * 1.5,
                          "S_rigid_twist": float(
                              2 * INS4.kin_of(M4N, tgN, cfgN))}
    Ls = np.array([grow[k]["L"] for k in grow])
    Ss = np.array([grow[k]["S_rigid_twist"] for k in grow])
    expo = float(np.polyfit(np.log(Ls), np.log(Ss), 1)[0])
    # where does the actual t32A rigid S live: own kin density
    pos32, _ = grid(32, ht32)
    v1_32, _, _, _, _ = frame(sp(M4t32), pos32)
    tg32 = tangent(M4t32, v1_32)
    dens = kin_density(M4t32, tg32, cfg32)
    dtot = float(dens.sum())
    shellS = []
    for a, b in zip([0, 4, 8, 12, 16, 20], [4, 8, 12, 16, 20, 99]):
        m = (r32 >= a) & (r32 < b)
        shellS.append({"r": [a, b], "dS2": float(2 * dens[m].sum())})
    out["C5"] = {
        "fjom0.2": {"S_env_twist": s_fj_tw["S_env"],
                    "S_rigid_twist": s_fj_tw["S_rigid"],
                    "S_env_tilt12": s_fj_t12["S_env"],
                    "rec_env_tw": 37.18221360010632,
                    "rec_rigid_tw": 576.8219558899685,
                    "rec_env_t12": 109.92492994146926,
                    "fd_err": s_fj_tw["tangent_fd_relerr"]},
        "t32A": {"S_env_twist": s_t32_tw["S_env"],
                 "S_rigid_twist": s_t32_tw["S_rigid"],
                 "S_env_tilt12": s_t32_t12["S_env"],
                 "rec_env_tw": 16.20668362638778,
                 "rec_rigid_tw": 440.39169514213893,
                 "rec_env_t12": 93.5599804556743,
                 "fd_err": s_t32_tw["tangent_fd_relerr"]},
        "t24A": {"S_env_twist": s_t24_tw["S_env"],
                 "S_rigid_twist": s_t24_tw["S_rigid"],
                 "S_env_tilt12": s_t24_t12["S_env"],
                 "rec_env_tw": 13.113772707475016,
                 "rec_rigid_tw": 255.46131925958989,
                 "rec_env_t12": 22.380584453209718,
                 "fd_err": s_t24_tw["tangent_fd_relerr"]},
        "vacuum_comm_fro": float(np.sqrt((comm * comm).sum())),
        "vacuum_comm_theory_sqrt2_x_0.7": float(np.sqrt(2) * 0.7),
        "ir_growth": grow, "ir_exponent_vs_L": expo,
        "kin_density_total_2x_vs_kinof": float(
            2 * dtot / s_t32_tw["S_rigid"]),
        "t32A_rigid_S_shells": shellS}
    print(" S_env tw fj", s_fj_tw["S_env"], "t32", s_t32_tw["S_env"],
          "IR exponent", expo, flush=True)

    # ---- C6: bridge arithmetic + the 2.0023 question
    print("== C6 bridge ==", flush=True)
    br = json.load(open(os.path.join(DATA, "m5_21_5_bridge.json")))
    checks, gs = [], []
    for tag, st in br["states"].items():
        for pname, p in st["pairings"].items():
            g = p["mu_lat"] * p["E_lat"] / (2 * np.pi * p["S_lat"])
            kneed = 2.0023 * p["S_lat"] / (
                2 * p["mu_lat"] * p["E_lat"])
            checks.append({
                "state": tag, "pairing": pname,
                "g_recomputed": float(g),
                "g_recorded": p["g_first_principles"],
                "match": bool(abs(g / p["g_first_principles"] - 1)
                              < 1e-12),
                "k_needed_recomputed": float(kneed),
                "k_needed_ratio_mine_over_rec": float(
                    kneed / p["k_needed_for_2.0023"])})
            gs.append((tag, pname, float(g)))
    gs.sort(key=lambda t: abs(t[2] / 2.0023 - 1))
    # my own fully-independent g for the two headline states
    g_fj = base["mu_clock"] * 6.740891714696357 / (
        2 * np.pi * s_fj_tw["S_env"])
    g_t32 = t32["mu_clock"] * 4.767566808218357 / (
        2 * np.pi * s_t32_tw["S_env"])
    eu, ev = INS4.e_parts(Mfj, cfg32)
    knr = [c["k_needed_ratio_mine_over_rec"] for c in checks]
    # the recorded column equals g_meas*S/(mu E) with g_meas =
    # 2.00231930436, i.e. the g = k mu E/S convention: 2x the value
    # consistent with the file's own k_derived (g = 2 k mu E/S)
    g_meas = 2.00231930436
    out["C6"] = {
        "E_fjom_recomputed_4d": float(eu + ev),
        "E_fjom_record": 6.740891714696357,
        "rows_checked": len(checks),
        "all_g_match": bool(all(c["match"] for c in checks)),
        "k_needed_ratio_min_max": [float(min(knr)),
                                   float(max(knr))],
        "k_needed_ratio_explained": float(0.5 * 2.0023 / g_meas),
        "closest_to_2.0023": [
            {"state": t, "pairing": p, "g": g,
             "rel_dev_from_2.0023": g / 2.0023 - 1}
            for t, p, g in gs[:4]],
        "my_g_fjom_mixed": float(g_fj),
        "my_g_t32A_mixed": float(g_t32),
        "k_derived_1_over_4pi": float(1 / (4 * np.pi))}
    print(" closest g:", gs[0], flush=True)

    with open(os.path.join(DATA, "m5_21_5_audit.json"), "w") as f:
        json.dump({k: v for k, v in out.items()}, f, indent=1,
                  default=lambda o: None)
    print("saved data/m5_21_5_audit.json", flush=True)


if __name__ == "__main__":
    main()
