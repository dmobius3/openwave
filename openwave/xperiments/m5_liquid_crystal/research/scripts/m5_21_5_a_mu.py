"""M5.21.5 G1: the mu/S instrument on verified-L states (EID-C
mechanism, per-voxel eigenframe channels).

Director route (Mermin-Ho, the m5_8_2r construction generalized off
the analytic frame):
    n(x)   = leading eigenvector of the spatial 3x3 block,
             sign-aligned to +r_hat (headless director; flip stats
             recorded; near-degenerate cells masked by eigengap)
    F_ij   = n . (d_i n x d_j n)          B_k = 1/2 eps_kij F_ij
    d_t n  from the channel flow M(phi) = R4 M R4^T (Rodrigues
             per-voxel about the channel axis, unit angular rate)
    E_i    = n . (d_t n x d_i n)
    j      = curl B - dE/dt               (displacement term, EID)
    mu     = 1/2 sum_mask r x j h^3       (+ mu(<R) radial windows)

Channels per state (axes from the local eigenframe v1 >= v2 >= v3,
right-handed v3 = v1 x v2):
    twist   axis v1  (rotation about the major axis; the m5_21_9
                      conjugation-clock direction; expect mu = 0)
    tilt12  axis v3  (rotates v1 -> v2; the EID-C tilt channel)
    tilt13  axis -v2 (rotates v1 -> v3; contrast)
    globalz axis z   (rigid internal rotation; control)

S side (the verified-L clock inertia, PHYSICAL-RATE convention:
tangent a = [G_loc, M] at unit angular rate, NO envelope, NO unit-
normalization): S = 2 * kin_of(M4, a, cfg4) so E_kin = omega^2 kin
and S = dE_kin/domega at omega = 1. The m5_21_9 carried-J convention
(envelope + unit-Frobenius a0) is ALSO read on the fixed-J states
with the conversion ratio recorded (the kin-convention flag made
quantitative).

Modes:
  gates                      P0: closed-form hedgehog + EID repro
  read tag=<t> path=<npz>    G1/G2 read on a state (3x3 or 4x4)
Out: data/m5_21_5_gates.json / data/m5_21_5_mu_<tag>.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import types

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_s4 = importlib.util.spec_from_file_location(
    "ins4", os.path.join(HERE, "m5_21_3_a_4d.py"))
INS4 = importlib.util.module_from_spec(_s4)
_s4.loader.exec_module(INS4)


# ================= EID shim (P0 gate only) =========================
def load_eid():
    """The m5_8_2r imports live under retired sandbox package paths;
    register the flat scripts/ files under those dotted names."""
    base = "openwave.xperiments.m5_liquid_crystal.research"
    chain = [
        (f"{base}.sandbox_v6.m5_6_2a_biaxial_hedgehog",
         "m5_6_2a_biaxial_hedgehog.py"),
        (f"{base}.sandbox_v8.m5_8_2a_4d_hamiltonian",
         "m5_8_2a_4d_hamiltonian.py"),
        (f"{base}.sandbox_v8.m5_8_2c1_full_evolution",
         "m5_8_2c1_full_evolution.py"),
        (f"{base}.sandbox_v8.m5_8_2cb_taichi_constrained",
         "m5_8_2cb_taichi_constrained.py"),
        (f"{base}.sandbox_vn.m5_8_2h_omega_attractor",
         "m5_8_2h_omega_attractor.py"),
    ]
    for pkg in (base.rsplit(".", 3)[0],):
        pass
    parents = set()
    for name, _ in chain:
        parts = name.split(".")
        for i in range(1, len(parts)):
            parents.add(".".join(parts[:i]))
    for p in sorted(parents):
        if p not in sys.modules:
            m = types.ModuleType(p)
            m.__path__ = []
            sys.modules[p] = m
    for name, fn in chain:
        if name in sys.modules:
            continue
        sp = importlib.util.spec_from_file_location(
            name, os.path.join(HERE, fn))
        mod = importlib.util.module_from_spec(sp)
        sys.modules[name] = mod
        sp.loader.exec_module(mod)
    sp = importlib.util.spec_from_file_location(
        "eid", os.path.join(HERE, "m5_8_2r_electron_id.py"))
    eid = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(eid)
    return eid


# ================= geometry helpers ================================
def cgrid(n, h):
    c = (n - 1) / 2.0
    idx = (np.indices((n, n, n)).transpose(1, 2, 3, 0) - c) * h
    return idx  # (n,n,n,3) positions


def central(f, ax, h):
    return (np.roll(f, -1, axis=ax) - np.roll(f, 1, axis=ax)) / (2 * h)


def frame_of(Msp, pos):
    """(v1, v2, v3, gap12, gap23): leading eigenframe of the spatial
    block. v1 sign-aligned to +r_hat (headless director). v2 (the
    middle/delta eigenvector) is per-voxel sign-ambiguous from eigh;
    an incoherent v2 decoheres every frame-based tangent field, so
    it is GAUGE-FIXED: aligned to the azimuthal reference
    t_hat = normalize(z_hat x v1) (fallback x_hat x v1 near the
    poles where v1 || z). Right-handed v3 = v1 x v2. The polar-line
    gauge defect is recorded via pole_frac."""
    lam, V = np.linalg.eigh(Msp)
    v1 = V[..., :, 2].copy()
    v2 = V[..., :, 1].copy()
    r = np.linalg.norm(pos, axis=-1, keepdims=True)
    rhat = np.divide(pos, np.maximum(r, 1e-12))
    dot = np.einsum("...a,...a->...", v1, rhat)
    flip = dot < 0
    v1[flip] = -v1[flip]
    zhat = np.zeros_like(v1)
    zhat[..., 2] = 1.0
    t = np.cross(zhat, v1)
    tn = np.linalg.norm(t, axis=-1)
    pole = tn < 0.1
    if pole.any():
        xhat = np.zeros_like(v1)
        xhat[..., 0] = 1.0
        t2 = np.cross(xhat, v1)
        t[pole] = t2[pole]
    t = t / np.maximum(np.linalg.norm(t, axis=-1, keepdims=True),
                       1e-300)
    p = np.cross(t, v1)  # polar-like transverse reference (the EID
    # analytic frame's col2 = eTheta is polar: gauge matches it)
    d2 = np.einsum("...a,...a->...", v2, p)
    v2[d2 < 0] = -v2[d2 < 0]
    v3 = np.cross(v1, v2)
    nz = np.linalg.norm(v3, axis=-1, keepdims=True)
    v3 = np.divide(v3, np.maximum(nz, 1e-300))
    gap12 = lam[..., 2] - lam[..., 1]
    gap23 = lam[..., 1] - lam[..., 0]
    stats = {"flip_frac": float(np.mean(flip)),
             "amb_frac": float(np.mean(np.abs(dot) < 0.05)),
             "pole_frac": float(np.mean(pole))}
    return v1, v2, v3, gap12, gap23, stats


def rodrigues(axis, phi):
    """(...,3) unit axis, scalar OR per-voxel phi -> (...,3,3)
    rotation (per-voxel phi carries the envelope-localized clock)."""
    a = axis
    K = np.zeros(a.shape[:-1] + (3, 3))
    K[..., 0, 1], K[..., 0, 2] = -a[..., 2], a[..., 1]
    K[..., 1, 0], K[..., 1, 2] = a[..., 2], -a[..., 0]
    K[..., 2, 0], K[..., 2, 1] = -a[..., 1], a[..., 0]
    eye = np.broadcast_to(np.eye(3), K.shape)
    if np.ndim(phi):
        phi = np.asarray(phi)[..., None, None]
    return eye + np.sin(phi) * K + (1 - np.cos(phi)) * (K @ K)


def advance(M, axis, phi):
    """M(phi) = R M R^T on the spatial block; 4x4 states rotate the
    (0,i) row/column consistently via R4 = blockdiag(1, R3)."""
    R3 = rodrigues(axis, phi)
    if M.shape[-1] == 3:
        return R3 @ M @ np.swapaxes(R3, -1, -2)
    R4 = np.zeros(M.shape)
    R4[..., 0, 0] = 1.0
    R4[..., 1:4, 1:4] = R3
    return R4 @ M @ np.swapaxes(R4, -1, -2)


def spatial(M):
    return M if M.shape[-1] == 3 else M[..., 1:4, 1:4]


def director(M, pos):
    v1, _, _, g12, g23, st = frame_of(spatial(M), pos)
    return v1, g12, g23, st


def align_to(n_new, n_ref):
    dot = np.einsum("...a,...a->...", n_new, n_ref)
    out = n_new.copy()
    out[dot < 0] = -out[dot < 0]
    return out


def mermin_B(n_hat, h):
    dn = [central(n_hat, ax, h) for ax in range(3)]
    F12 = np.einsum("...a,...a->...", n_hat, np.cross(dn[0], dn[1]))
    F13 = np.einsum("...a,...a->...", n_hat, np.cross(dn[0], dn[2]))
    F23 = np.einsum("...a,...a->...", n_hat, np.cross(dn[1], dn[2]))
    return np.stack([F23, -F13, F12], axis=-1), dn


def E_of(n_hat, dtn, dn):
    return np.stack(
        [np.einsum("...a,...a->...", n_hat,
                   np.cross(dtn, dn[ax])) for ax in range(3)],
        axis=-1)


def curl(V, h):
    return np.stack(
        [central(V[..., 2], 1, h) - central(V[..., 1], 2, h),
         central(V[..., 0], 2, h) - central(V[..., 2], 0, h),
         central(V[..., 1], 0, h) - central(V[..., 0], 1, h)],
        axis=-1)


# ================= the mu read =====================================
def axis_of(channel, v1, v2, v3, n):
    if channel == "twist":
        return v1
    if channel == "tilt12":
        return v3
    if channel == "tilt13":
        return -v2
    if channel == "globalz":
        z = np.zeros_like(v1)
        z[..., 2] = 1.0
        return z
    raise ValueError(channel)


def mu_read(M, h, channel, mask, dphi=1e-4, dphi2=1e-3,
            rwin=(4, 6, 8, 10, 12, 14, 16, 100), w=None):
    """w = None: rigid unit-rate flow (EID-comparable; measured to
    be pin-layer/tail-dominated on relaxed pinned states). w given:
    the envelope-localized clock (per-voxel angle phi*w, the
    m5_21_9 state-of-record flow shape) -> the mu of the clock the
    state actually carries; pairs consistently with S_env."""
    n = M.shape[0]
    pos = cgrid(n, h)
    v1, v2, v3, g12, g23, st = frame_of(spatial(M), pos)
    ax = axis_of(channel, v1, v2, v3, n)
    wf = 1.0 if w is None else w

    def nb(phi):
        if phi == 0.0:
            return v1
        Mp = advance(M, ax, phi * wf)
        np_, _, _, _ = director(Mp, pos)
        return align_to(np_, v1)

    n0 = v1
    B0, dn0 = mermin_B(n0, h)
    dtn = (nb(dphi) - n0) / dphi
    E0 = E_of(n0, dtn, dn0)
    n1 = nb(dphi2)
    B1, dn1 = mermin_B(n1, h)
    dtn1 = (align_to(nb(dphi2 + dphi), n1) - n1) / dphi
    E1 = E_of(n1, dtn1, dn1)
    j = curl(B0, h) - (E1 - E0) / dphi2
    jcb = curl(B0, h)
    r = np.linalg.norm(pos, axis=-1)
    rows = {}
    for R in rwin:
        m = mask & (r < R)
        mu = 0.5 * np.cross(pos, j)[m].sum(axis=0) * h ** 3
        mucb = 0.5 * np.cross(pos, jcb)[m].sum(axis=0) * h ** 3
        muck = mu - mucb  # the clock-induced part (E-term only)
        rows[f"R{R}"] = {"mu": [float(x) for x in mu],
                         "mu_norm": float(np.linalg.norm(mu)),
                         "mu_curlB_vec": [float(x) for x in mucb],
                         "mu_curlB_only": float(np.linalg.norm(mucb)),
                         "mu_clock_vec": [float(x) for x in muck],
                         "mu_clock": float(np.linalg.norm(muck)),
                         "vox": int(m.sum())}
    full = rows[f"R{rwin[-1]}"]
    return {"channel": channel, "mu_norm": full["mu_norm"],
            "mu_vec": full["mu"], "mu_curlB_only":
            full["mu_curlB_only"], "mu_clock": full["mu_clock"],
            "radial": rows,
            "dtn_max": float(np.abs(dtn[mask]).max()),
            "frame_stats": st}


def masks_of(M, h, margin=2, gap_min=0.02):
    """Two masks: director-grade (gap12: v1 well-defined; twist +
    globalz channels) and frame-grade (also gap23: v2/v3 well-
    defined off the uniaxial-escape rod cores; tilt channels)."""
    n = M.shape[0]
    pos = cgrid(n, h)
    _, g12, g23, _ = director(M, pos)
    inter = np.zeros((n, n, n), bool)
    inter[margin:-margin, margin:-margin, margin:-margin] = True
    good12 = g12 > gap_min
    good23 = g23 > gap_min
    m_dir = inter & good12
    m_frame = inter & good12 & good23
    return m_dir, m_frame, {
        "masked_frac_dir": float(1 - good12[inter].mean()),
        "masked_frac_frame":
            float(1 - (good12 & good23)[inter].mean()),
        "gap_min": gap_min}


# ================= the S read ======================================
def S_read(M4, cfg4, channel="twist"):
    """S = 2 kin at unit angular rate, physical tangent [G_loc, M].
    TWO variants: rigid (no envelope: IR-extensive, the M5.21.8
    constant-omega pathology, reported as diagnostic) and env (the
    m5_21_9 state-of-record clock localization w(r), renv = 10:
    w = 1 through the mu-active core, truncates only the IR tail:
    the PRIMARY S)."""
    n = M4.shape[0]
    pos = cgrid(n, cfg4["h"])
    v1, v2, v3, _, _, _ = frame_of(spatial(M4), pos)
    ax = axis_of(channel, v1, v2, v3, n)
    G = np.zeros(M4.shape)
    a = ax
    G[..., 1, 2], G[..., 1, 3] = -a[..., 2], a[..., 1]
    G[..., 2, 1], G[..., 2, 3] = a[..., 2], -a[..., 0]
    G[..., 3, 1], G[..., 3, 2] = -a[..., 1], a[..., 0]
    tang = G @ M4 - M4 @ G
    w = INS4.envelope(cfg4)[..., None, None]
    kin_r = INS4.kin_of(M4, tang, cfg4)
    kin_e = INS4.kin_of(M4, w * tang, cfg4)
    return {"channel": channel,
            "kin_rigid": float(kin_r),
            "S_rigid": float(2.0 * kin_r),
            "kin_env": float(kin_e),
            "S_env": float(2.0 * kin_e),
            "S_phys_rate": float(2.0 * kin_e),
            "tang_fro": float(np.sqrt(np.sum(tang * tang)))}


# ================= gates (P0) ======================================
def gates():
    out = {}
    # (i) closed form: analytic n = r_hat on a 24^3 grid, h = 1.5:
    # sum_{i<j} F_ij^2 = 1/r^4 and |B| flux = 4pi (unit monopole)
    n, h = 24, 1.5
    pos = cgrid(n, h)
    r = np.linalg.norm(pos, axis=-1)
    nhat = pos / np.maximum(r[..., None], 1e-12)
    B, dn = mermin_B(nhat, h)
    F2 = np.sum(B * B, axis=-1)      # = sum_{i<j} F_ij^2
    band = (r > 6) & (r < 12)
    ratio = F2[band] * r[band] ** 4
    out["closed_form"] = {
        "r4F2_mean": float(ratio.mean()),
        "r4F2_p95_dev": float(np.percentile(np.abs(ratio - 1), 95)),
        "target": 1.0}
    # flux through a coordinate shell (sum of B_r over a thin band)
    shell = (r > 8) & (r < 8 + h)
    Br = np.einsum("...a,...a->...", B, nhat)
    flux = float((Br[shell] / (r[shell] ** 0)).mean()
                 * 4 * np.pi * 8 ** 2)
    out["closed_form"]["flux_vs_4pi"] = flux / (4 * np.pi)
    # (ii) EID reproduction on their analytic 24^3 seed
    eid = load_eid()
    gr, W, M, Mth, act = eid.grid_and_seed(24, plane=(1, 2), phi=0.0,
                                           b=0.0)
    heid = gr["h"]
    # their instrument numbers on this seed (tilt mu, twist L_int)
    res = {}
    n0 = eid.director_of(W)
    _, Wy, _, _, _ = eid.grid_and_seed(24, plane=(1, 2), phi=1e-4,
                                       b=0.0)
    E0B0 = eid.EB_of(n0, (eid.director_of(Wy) - n0) / 1e-4, heid)
    _, Wz, _, _, _ = eid.grid_and_seed(24, plane=(1, 2), phi=1e-3,
                                       b=0.0)
    n1 = eid.director_of(Wz)
    _, Wz2, _, _, _ = eid.grid_and_seed(24, plane=(1, 2),
                                        phi=1e-3 + 1e-4, b=0.0)
    E1B1 = eid.EB_of(n1, (eid.director_of(Wz2) - n1) / 1e-4, heid)
    mu_theirs, _, _ = eid.mu_of(gr, E0B0, E1B1, 1e-3, act)
    res["mu_tilt_theirs"] = float(np.linalg.norm(mu_theirs))
    # MY pipeline on the SAME M (their act mask for comparability)
    mine = mu_read(M, heid, "tilt12", act)
    res["mu_tilt12_mine"] = mine["mu_norm"]
    res["record_2r"] = 0.2209
    # twist channel EM-silence on the analytic seed, my pipeline
    res["mu_twist_mine"] = mu_read(M, heid, "twist", act)["mu_norm"]
    out["eid_repro"] = res
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, "m5_21_5_gates.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1), flush=True)
    return out


# ================= state reads =====================================
def load_state(path):
    z = np.load(path)
    M = z["M"].astype(np.float64)
    h = float(z["h"]) if "h" in z else 1.5
    return M, h


def e_of_record(tag, path):
    """The state's OWN functional-of-record energy: 3x3 census
    states carry the T2 3x3 E_end in their row JSON; 4x4 states use
    the 4D e_parts. Mixing conventions breaks the m_e anchor."""
    base = os.path.basename(path)
    row = None
    for cand in (base.replace("_end_", "_row_").replace(".npz",
                                                        ".json"),):
        p = os.path.join(DATA, cand)
        if os.path.exists(p):
            row = json.load(open(p))
            break
    return row


def read_state(tag, path, channels=("twist", "tilt12", "tilt13",
                                    "globalz")):
    M, h = load_state(path)
    n = M.shape[0]
    is4 = M.shape[-1] == 4
    m_dir, m_frame, mstat = masks_of(M, h)
    out = {"tag": tag, "path": os.path.basename(path), "n": n,
           "h": h, "is4d": bool(is4), "mask": mstat}
    cfg4 = INS4.base_cfg(s=-1.0, n=n, L=n * h)
    if is4:
        eu, ev = INS4.e_parts(M, cfg4)
        out["E_lat"] = float(eu + ev)  # 4D functional of record
        out["E_u"], out["E_v"] = float(eu), float(ev)
        M4 = M
    else:
        M4 = INS4.embed34(M, cfg4)
        eu, ev = INS4.e_parts(M4, cfg4)
        out["E_lat_4dconv"] = float(eu + ev)  # sensitivity only
        row = e_of_record(tag, path)
        if row is not None:
            out["E_lat"] = float(row["E_end"])  # T2 3x3 of record
            out["E_row_src"] = row["tag"]
        else:
            out["E_lat"] = float(eu + ev)
            out["E_row_src"] = "4dconv_fallback"
    w_env = INS4.envelope(cfg4)
    for ch in channels:
        mask = m_frame if ch.startswith("tilt") else m_dir
        out[f"mu_{ch}"] = mu_read(M, h, ch, mask, w=w_env)
        out[f"S_{ch}"] = S_read(M4, cfg4, ch)
        print(tag, ch, "mu_env", out[f"mu_{ch}"]["mu_norm"],
              "mu_clock_env", out[f"mu_{ch}"]["mu_clock"],
              "S_env", out[f"S_{ch}"]["S_env"], flush=True)
    # the rigid tilt12 read kept as the EID-comparable diagnostic
    out["mu_tilt12_rigid"] = mu_read(M, h, "tilt12", m_frame)
    # the m5_21_9 carried-J convention on fixed-J states
    if is4:
        w = INS4.envelope(cfg4)[..., None, None]
        pos = cgrid(n, h)
        v1, v2, v3, _, _, _ = frame_of(spatial(M4), pos)
        G = np.zeros(M4.shape)
        a = v1
        G[..., 1, 2], G[..., 1, 3] = -a[..., 2], a[..., 1]
        G[..., 2, 1], G[..., 2, 3] = a[..., 2], -a[..., 0]
        G[..., 3, 1], G[..., 3, 2] = -a[..., 1], a[..., 0]
        A = w * (G @ M4 - M4 @ G)
        norm = float(np.sqrt(np.sum(A * A)))
        a0u = A / max(norm, 1e-300)
        kin_u = INS4.kin_of(M4, a0u, cfg4)
        out["fixedj_conv"] = {"kin_unitnorm": float(kin_u),
                              "a0_env_fro": norm}
    with open(os.path.join(DATA, f"m5_21_5_mu_{tag}.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print("saved", tag, flush=True)
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "gates"
    if mode == "gates":
        gates()
    elif mode == "read":
        kw = dict(a.split("=", 1) for a in sys.argv[2:])
        read_state(kw["tag"], os.path.join(DATA, kw["path"]))
    else:
        raise SystemExit(f"unknown mode {mode}")
