"""M5.21.4 instrument: 2-defect seeds, signed charge, E(d), capture dynamics.

Consumes the M5.21.2b certified 3D instrument (INS: T2 term, sym stencil,
FIRE) and the M5.21.6 dynamics stack (DEC: leapfrog, sponge, e_tot,
loop_read) by import; new physics code only.

Seeds (all exactly on the (1, delta, 0) spectrum manifold, isotropic-core
blended like INS.seed3):
  antipair   alpha = atan2(rho, z-zt) + pi - atan2(rho, z-zb)
             (radial hedgehog + mirror hedgehog; far field uniform)
  samepair   inverse-stereographic product ansatz u =
             tan(th_t/2)*tan(th_b/2)*e^{2 i phi} with a smooth +z escape
             blend on the inter-core axis tube (integer winding escapes);
             far field = the charge-2 texture
  single / mirror  one-defect limits (regression + charge-sign gates)

Signed charge (the M5.21.5 emergent-monopole convention): B_i = 1/2
eps_ijk n.(d_j n x d_k n) (Mermin-Ho); Q = (1/4 pi) closed-surface sum of
B.dS over a lattice cube. Relative signs need ONE global orientation of
the director lift: analytic seeds are smooth by construction; eigh-derived
states are oriented by continuity sweeps (conflict count reported).

Phases (CLI): gates | ladder | evolve | collect
  ladder: equal-depth healed E(d) reads, both sectors, n=32 pinned
  evolve: the antipair capture run, n=48 free, leapfrog + sponge
"""

import importlib.util
import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")

_s1 = importlib.util.spec_from_file_location(
    "ins", os.path.join(HERE, "m5_21_2b_a_instrument.py"))
INS = importlib.util.module_from_spec(_s1)
_s1.loader.exec_module(INS)
_s2 = importlib.util.spec_from_file_location(
    "dec", os.path.join(HERE, "m5_21_6_a_decay.py"))
DEC = importlib.util.module_from_spec(_s2)
_s2.loader.exec_module(DEC)


# ================= seeds =================
def _nhat_from_alpha(n, h, alpha):
    X, Y, Z = INS.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    rhos = np.where(rho < 1e-12, 1e-12, rho)
    cphi, sphi = X / rhos, Y / rhos
    near = rho < 1e-9
    if np.any(near):
        cphi[near], sphi[near] = 1.0, 0.0
    sa, ca = np.sin(alpha), np.cos(alpha)
    return np.stack([sa * cphi, sa * sphi, ca], axis=-1)


def _tensor_from_nhat(n, h, delta, nhat, centers, r_c=4.0):
    """M = n n^T + delta phihat phihat^T, isotropic-blended at each core."""
    X, Y, Z = INS.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    rhos = np.where(rho < 1e-12, 1e-12, rho)
    phihat = np.stack([-Y / rhos, X / rhos, np.zeros_like(Z)], axis=-1)
    near = rho < 1e-9
    if np.any(near):
        phihat[near] = np.array([0.0, 1.0, 0.0])
    # phihat is not exactly orthogonal to nhat only where the escape blend
    # tilted nhat off the meridional plane; re-orthogonalize there
    dot = np.einsum("...a,...a->...", phihat, nhat)[..., None]
    ph = phihat - dot * nhat
    ph = ph / np.maximum(np.linalg.norm(ph, axis=-1)[..., None], 1e-300)
    S = (nhat[..., :, None] * nhat[..., None, :]
         + delta * ph[..., :, None] * ph[..., None, :])
    a = (1.0 + delta) / 3.0
    w = np.ones_like(rho)
    for zc in centers:
        di = np.sqrt(X * X + Y * Y + (Z - zc) ** 2)
        w = w * (1.0 - np.exp(-((di / r_c) ** 2)))
    return w[..., None, None] * S + (1.0 - w[..., None, None]) * (a * np.eye(3))


def seed_pair(cfg, kind, d):
    """kind: 'anti' | 'same' | 'single' | 'mirror'. d = separation."""
    n, h, delta = cfg["n"], cfg["h"], cfg["delta"]
    X, Y, Z = INS.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    zt, zb = +d / 2.0, -d / 2.0
    tht = np.arctan2(rho, Z - zt)
    thb = np.arctan2(rho, Z - zb)
    if kind == "single":
        return _tensor_from_nhat(
            n, h, delta, _nhat_from_alpha(n, h, np.arctan2(rho, Z)), [0.0])
    if kind == "mirror":
        return _tensor_from_nhat(
            n, h, delta, _nhat_from_alpha(n, h, np.pi - np.arctan2(rho, Z)),
            [0.0])
    if kind == "anti":
        return _tensor_from_nhat(
            n, h, delta, _nhat_from_alpha(n, h, tht + np.pi - thb), [zt, zb])
    if kind == "same":
        tt = np.clip(np.tan(0.5 * tht), 0.0, 1e6)
        tb = np.clip(np.tan(0.5 * thb), 0.0, 1e6)
        umag = np.clip(tt * tb, 0.0, 1e12)
        rhos = np.where(rho < 1e-12, 1e-12, rho)
        c2p = (X * X - Y * Y) / (rhos * rhos)      # cos 2phi
        s2p = -(2.0 * X * Y) / (rhos * rhos)       # -sin 2phi: +2 sector
        den = 1.0 + umag * umag
        nhat = np.stack([2.0 * umag * c2p / den, 2.0 * umag * s2p / den,
                         (1.0 - umag * umag) / den], axis=-1)
        # smooth +z escape on the inter-core axis tube (integer winding)
        r_t = 3.0
        zed = 0.5 * (np.tanh((Z - zb) / 2.0) - np.tanh((Z - zt) / 2.0))
        eta = np.exp(-((rho / r_t) ** 2)) * zed
        nhat = nhat + eta[..., None] * np.array([0.0, 0.0, 1.0])
        nhat = nhat / np.maximum(
            np.linalg.norm(nhat, axis=-1)[..., None], 1e-300)
        return _tensor_from_nhat(n, h, delta, nhat, [zt, zb])
    raise ValueError(kind)


# ================= signed charge (Mermin-Ho flux) =================
def orient_v1(M):
    """leading eigenvector of the 3x3 field, oriented by continuity:
    axis-sweep first guess, then vectorized majority-vote fix-up passes
    (flip any voxel anti-aligned with the majority of its 6 neighbors)
    until stable; returns (nhat, n_conflicts). Conflicts concentrate at
    degenerate cores; the flux cubes read faces away from them."""
    _, vec = np.linalg.eigh(M)
    v = vec[..., :, 2].copy()
    sgn = np.ones(v.shape[:3])
    for ax in range(3):
        dots = np.einsum("...a,...a->...",
                         np.roll(v * sgn[..., None], 1, axis=ax),
                         v * sgn[..., None])
        flip = np.where(dots < 0.0, -1.0, 1.0)
        idx = [slice(None)] * 3
        idx[ax] = 0
        flip[tuple(idx)] = 1.0
        sgn = sgn * np.cumprod(flip, axis=ax)
    for _ in range(60):
        vo = v * sgn[..., None]
        vote = np.zeros(v.shape[:3])
        for ax in range(3):
            up = np.einsum("...a,...a->...", np.roll(vo, 1, ax), vo)
            dn = np.einsum("...a,...a->...", np.roll(vo, -1, ax), vo)
            sl = [slice(None)] * 3
            sl[ax] = 0
            up[tuple(sl)] = 0.0
            sl[ax] = -1
            dn[tuple(sl)] = 0.0
            vote = vote + up + dn
        bad = vote < 0.0
        if not bad.any():
            break
        sgn = np.where(bad, -sgn, sgn)
    vo = v * sgn[..., None]
    ncf = 0
    for ax in range(3):
        dd = np.einsum("...a,...a->...", np.roll(vo, 1, ax), vo)
        sl = [slice(None)] * 3
        sl[ax] = slice(1, None)
        ncf += int((dd[tuple(sl)] < 0).sum())
    return vo, ncf


def mermin_B(nhat, h):
    dn = [_central(nhat, ax, h) for ax in range(3)]
    F12 = np.einsum("...a,...a->...", nhat, np.cross(dn[0], dn[1]))
    F13 = np.einsum("...a,...a->...", nhat, np.cross(dn[0], dn[2]))
    F23 = np.einsum("...a,...a->...", nhat, np.cross(dn[1], dn[2]))
    return np.stack([F23, -F13, F12], axis=-1)


def _central(f, ax, h):
    out = (np.roll(f, -1, axis=ax) - np.roll(f, 1, axis=ax)) / (2.0 * h)
    sl0 = [slice(None)] * f.ndim
    sln = [slice(None)] * f.ndim
    sl0[ax], sln[ax] = 0, -1
    out[tuple(sl0)] = out[tuple(sln)] = 0.0
    return out


def cube_flux(B, cfg, center_z, half):
    """Q = (1/4pi) sum_faces B.dS over the lattice cube |x|,|y|<=half,
    |z-center_z|<=half (physical units)."""
    n, h = cfg["n"], cfg["h"]
    c = (n - 1) / 2.0
    iz = int(round(center_z / h + c))
    k = int(round(half / h))
    i0, i1 = int(round(c)) - k, int(round(c)) + k
    z0, z1 = iz - k, iz + k
    if min(i0, z0) < 1 or max(i1, z1) > n - 2:
        return float("nan")
    s = 0.0
    s += B[i1, i0:i1 + 1, z0:z1 + 1, 0].sum() - \
        B[i0, i0:i1 + 1, z0:z1 + 1, 0].sum()
    s += B[i0:i1 + 1, i1, z0:z1 + 1, 1].sum() - \
        B[i0:i1 + 1, i0, z0:z1 + 1, 1].sum()
    s += B[i0:i1 + 1, i0:i1 + 1, z1, 2].sum() - \
        B[i0:i1 + 1, i0:i1 + 1, z0, 2].sum()
    return float(s * h * h / (4.0 * np.pi))


def charge_suite(M, cfg, d, analytic_nhat=None):
    """per-core + far cube fluxes at two sizes each."""
    if analytic_nhat is not None:
        nhat, ncf = analytic_nhat, 0
    else:
        nhat, ncf = orient_v1(M)
    B = mermin_B(nhat, cfg["h"])
    zt, zb = +d / 2.0, -d / 2.0
    hc = max(2.0 * cfg["h"], 0.28 * d)
    far0 = 0.5 * cfg["L"] - 4.0 * cfg["h"]
    out = {"n_conflicts": ncf}
    for lab, zc, halves in [
            ("top", zt, (hc, hc + cfg["h"])),
            ("bot", zb, (hc, hc + cfg["h"])),
            ("far", 0.0, (far0 - 2 * cfg["h"], far0))]:
        out[lab] = [cube_flux(B, cfg, zc, hv) for hv in halves]
    return out


# ================= core tracking =================
def core_zs(M, cfg, with_gaps=False):
    """two core positions from the on-axis biaxiality dip (min eigen-gap
    along the axis column), sub-voxel by parabolic fit; search window
    excludes 4 voxels at each boundary (pin-shell artifacts). The dip
    VALUE distinguishes a real core (gap -> 0) from no-core (bulk gap)."""
    n, h = cfg["n"], cfg["h"]
    c = (n - 1) // 2
    col = M[c, c, :, :, :]
    lam = np.linalg.eigvalsh(col)
    gap = np.minimum(lam[:, 2] - lam[:, 1], lam[:, 1] - lam[:, 0])
    zax = (np.arange(n) - (n - 1) / 2.0) * h
    mid = n // 2
    zs, gv = [], []
    for sl in [slice(mid, n - 4), slice(4, mid)]:
        i = np.argmin(gap[sl]) + sl.start
        y0, y1, y2 = gap[i - 1], gap[i], gap[i + 1]
        den = (y0 - 2 * y1 + y2)
        off = 0.5 * (y0 - y2) / den if abs(den) > 1e-30 else 0.0
        zs.append(float(zax[i] + np.clip(off, -1, 1) * h))
        gv.append(float(y1))
    return (zs, gv) if with_gaps else zs  # [z_top, z_bot]


# ================= drivers =================
def heal(cfg, kind, d, it):
    M0 = seed_pair(cfg, kind, d)
    if cfg["bc"] == "pinned":
        free = ~INS.pin_shell(cfg["n"], cfg["h"])
    else:
        free = np.ones((cfg["n"],) * 3, dtype=bool)
    M, states, info = INS.fire(M0, cfg, free, it,
                               log_every=max(it // 4, 1),
                               tag=f"{kind}_d{d:g}")
    return M0, M, info


def ladder(kinds=("anti", "same"), ds=(12.0, 18.0, 24.0), it=1500, n=32):
    cfg = INS.base_cfg(term="T2", n=n, L=1.5 * n, bc="pinned",
                       stencil="sym")
    rows = []
    for kind in kinds:
        for d in ds:
            t0 = time.time()
            M0, M, info = heal(cfg, kind, d, it)
            e_u, e_d_, e_v = INS.e_parts(M, cfg)
            q = charge_suite(M, cfg, d)
            zs, gv = core_zs(M, cfg, with_gaps=True)
            row = {"kind": kind, "d": d, "it": it, "n": n,
                   "E": float(e_u + e_d_ + e_v), "E_u": float(e_u),
                   "stop": info["stop"], "cores": zs, "core_gaps": gv,
                   "charge": q, "wall_s": time.time() - t0}
            rows.append(row)
            print(json.dumps(row), flush=True)
            np.savez_compressed(os.path.join(
                DATA, f"m5_21_4_lad_{kind}_d{d:g}_n{n}_it{it}.npz"),
                M=M.astype(np.float32))
            with open(os.path.join(
                    DATA, f"m5_21_4_ladder_it{it}.json"), "w") as f:
                json.dump(rows, f, indent=1)
    # the single-defect reference + far-coefficient self-calibration
    _, Ms, _ = heal(cfg, "single", 0.0, it)
    e_single = float(sum(INS.e_parts(Ms, cfg)))
    r = DEC.rgrid(cfg)
    dens = _edens(Ms, cfg)
    sel = (r > 9.0) & (r < 0.5 * cfg["L"] - 4 * cfg["h"])
    c2fit = float(np.median(dens[sel] * r[sel] ** 4) / 8.0)
    with open(os.path.join(
            DATA, f"m5_21_4_ladder_it{it}.json"), "w") as f:
        json.dump({"rows": rows, "E_single": e_single,
                   "c2_selfcal": c2fit,
                   "coulomb_pred_coeff_64pi_c2": 64.0 * np.pi * c2fit},
                  f, indent=1)
    print(f"single E {e_single:.6f}  c2_selfcal {c2fit:.6e}", flush=True)


def _edens(M, cfg):
    """per-voxel energy density matching the INS functional split
    (curvature + potential; sym = average of fwd/bwd densities)."""
    h, st = cfg["h"], cfg["stencil"]
    stl = ("fwd", "bwd") if st == "sym" else (st,)
    e = 0.0
    for s in stl:
        dMs = [INS.d1(M, ax, h, s) for ax in range(3)]
        for a in range(3):
            for b in range(a + 1, 3):
                C = dMs[a] @ dMs[b] - dMs[b] @ dMs[a]
                e = e + (4.0 / len(stl)) * np.einsum(
                    "...ab,...ab->...", C, C)
    return e + INS.v_density(M, cfg)


def evolve(d0=24.0, n=48, steps=4000, dt=0.025, snap_every=100,
           heal_it=800, gam_in=0.02, tag="cap"):
    cfg = INS.base_cfg(term="T2", n=n, L=1.5 * n, bc="free", stencil="sym")
    free = np.ones((n, n, n), dtype=bool)
    M0 = seed_pair(cfg, "anti", d0)
    M, _, _ = INS.fire(M0, cfg, free, heal_it, log_every=200,
                       tag=f"{tag}_heal")
    V = np.zeros_like(M)
    gam = DEC.sponge(cfg) + gam_in
    fr = free[..., None, None].astype(float)
    rows, snaps = [], [{"it": 0, "M": M.astype(np.float32)}]
    E0, KE0 = DEC.e_tot(M, V, cfg)
    absorbed = 0.0
    t0 = time.time()
    for it in range(1, steps + 1):
        Eb, KEb = DEC.e_tot(M, V, cfg)
        M, V = DEC.leap_step(M, V, cfg, fr, gam, dt)
        Ea, KEa = DEC.e_tot(M, V, cfg)
        absorbed += (Eb + KEb) - (Ea + KEa)
        if not np.isfinite(Ea):
            print("NON-FINITE, stop", flush=True)
            break
        if it % snap_every == 0:
            zs = core_zs(M, cfg)
            q = charge_suite(M, cfg, max(zs[0] - zs[1], 4 * cfg["h"]))
            rows.append({"it": it, "t": it * dt, "E": Ea, "KE": KEa,
                         "absorbed": absorbed,
                         "ledger": Ea + KEa + absorbed - (E0 + KE0),
                         "cores": zs, "dsep": zs[0] - zs[1],
                         "charge": {k: q[k] for k in ("top", "bot", "far")},
                         "conflicts": q["n_conflicts"]})
            print(json.dumps(rows[-1]), flush=True)
            if it % (4 * snap_every) == 0:
                snaps.append({"it": it, "M": M.astype(np.float32)})
            with open(os.path.join(
                    DATA, f"m5_21_4_ev_{tag}_rows.json"), "w") as f:
                json.dump({"cfg_n": n, "d0": d0, "dt": dt,
                           "E0": E0, "KE0": KE0, "rows": rows}, f, indent=1)
    np.savez_compressed(
        os.path.join(DATA, f"m5_21_4_ev_{tag}.npz"),
        **{f"M_it{s['it']}": s["M"] for s in snaps})
    print(f"evolve done {time.time() - t0:.0f}s", flush=True)


# ================= gates =================
def gates():
    out = {}
    cfg = INS.base_cfg(term="T2", n=32, L=48.0, bc="pinned", stencil="sym")
    n, h, delta = cfg["n"], cfg["h"], cfg["delta"]
    # GS1 single-defect regression: seed_pair('single') == INS.seed3 A
    Ma = seed_pair(cfg, "single", 0.0)
    Mb = INS.seed3(n, h, delta, "A")
    out["GS1_single_vs_seed3"] = float(np.abs(Ma - Mb).max())
    # GS2 spectrum manifold: pair seeds outside cores sit on (1,delta,0)
    r = DEC.rgrid(cfg)
    X, Y, Z = INS.coords(n, h)
    for kind in ("anti", "same"):
        M = seed_pair(cfg, kind, 18.0)
        dt_ = np.sqrt(X * X + Y * Y + (Z - 9.0) ** 2)
        db_ = np.sqrt(X * X + Y * Y + (Z + 9.0) ** 2)
        away = (dt_ > 10.0) & (db_ > 10.0)
        lam = np.linalg.eigvalsh(M[away])
        tgt = np.array([0.0, delta, 1.0])
        out[f"GS2_{kind}_spec_dev"] = float(
            np.abs(np.sort(lam, axis=-1) - tgt).max())
    # GQ charge table (analytic nhat for single/mirror; eigh for pairs)
    rows = {}
    for kind, d in [("single", 0.0), ("mirror", 0.0),
                    ("anti", 18.0), ("same", 18.0)]:
        M = seed_pair(cfg, kind, d)
        q = charge_suite(M, cfg, d if d > 0 else 18.0)
        rows[kind] = q
    out["GQ"] = rows
    # GD dynamics regression on the healed antipair (short)
    cfgf = INS.base_cfg(term="T2", n=32, L=48.0, bc="free", stencil="sym")
    free = np.ones((32, 32, 32), dtype=bool)
    M0 = seed_pair(cfgf, "anti", 18.0)
    M, _, _ = INS.fire(M0, cfgf, free, 300, log_every=300, tag="GD_heal")
    V = np.zeros_like(M)
    fr = free[..., None, None].astype(float)
    E0, K0 = DEC.e_tot(M, V, cfgf)
    # small kick so KE is nonzero, then gamma=0 conservation
    rng = np.random.default_rng(214)
    V = 0.01 * rng.normal(size=V.shape)
    V = 0.5 * (V + np.swapaxes(V, -1, -2)) * fr
    E0, K0 = DEC.e_tot(M, V, cfgf)
    Mg, Vg = M.copy(), V.copy()
    for _ in range(200):
        Mg, Vg = DEC.leap_step(Mg, Vg, cfgf, fr, 0.0, 0.025)
    E1, K1 = DEC.e_tot(Mg, Vg, cfgf)
    out["GD_econs_rel"] = abs((E1 + K1) - (E0 + K0)) / max(E0 + K0, 1e-12)
    gam = DEC.sponge(cfgf)
    Ms, Vs = M.copy(), V.copy()
    tots, ok = [], True
    for i in range(200):
        Ms, Vs = DEC.leap_step(Ms, Vs, cfgf, fr, gam, 0.025)
        if i % 20 == 0:
            e, k = DEC.e_tot(Ms, Vs, cfgf)
            tots.append(e + k)
    ok = all(tots[i + 1] <= tots[i] + 1e-9 for i in range(len(tots) - 1))
    out["GD_sponge_monotone"] = bool(ok)
    with open(os.path.join(DATA, "m5_21_4_gates.json"), "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(json.dumps(out, indent=1, default=float), flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "gates"
    kv = {}
    for a in sys.argv[2:]:
        k, v = a.split("=")
        kv[k] = float(v) if "." in v or v.lstrip("-").isdigit() else v
    if cmd == "gates":
        gates()
    elif cmd == "ladder":
        ladder(it=int(kv.get("it", 1500)), n=int(kv.get("n", 32)))
    elif cmd == "evolve":
        evolve(d0=float(kv.get("d0", 24.0)), n=int(kv.get("n", 48)),
               steps=int(kv.get("steps", 4000)),
               gam_in=float(kv.get("gam", 0.02)),
               heal_it=int(kv.get("heal", 800)),
               tag=str(kv.get("tag", "cap")))
    else:
        raise SystemExit(f"unknown cmd {cmd}")
