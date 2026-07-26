"""M5.21.6: the 3D mu/tau decay instrument.

Consumes the M5.21.2b certified 3D stack (T2 eigenvalue-penalty term,
sym stencil; gates green there) by IMPORT: seeds, energy, gradient,
FIRE are the audited originals. New physics added here:

    relax   T2 endpoint regeneration under m5_21_6_ naming
            (B/C at n=32 pinned; A/B/C at n=48 FREE, the honest
            cross-sector arena per the census N>=48 protection read)
    gates   new-code gates (cap 3 tries, pre-registered):
            GL1 leapfrog gamma=0 conserves E_tot = E + KE (<1% / 200 steps)
            GL2 sponge on -> E_tot monotone down
            GK  kick operators keep M exactly symmetric + pinned shell
                untouched
    p1      Arena 1 kick ladder (pinned, in-sector only): K1 random
            smooth core kick eps in {0.05,0.15,0.4} RMS-relative;
            K2 core-twist rotation kick theta in {30,60,90} deg about
            z with the seed envelope; FIRE 4000; classify returned /
            moved (E + core spectrum vs the start endpoint)
    evolve  damped wave evolution M_tt = -grad E / h^3 - gamma(r) M_t
            (leapfrog kick-drift-kick), interior gamma = 0, absorbing
            sponge ramp in the outer shell; per-snap reads:
            E(r<10) / mid / sponge-absorbed cumulative, core
            spectrum-deviation vs frame-rotation (rotation-vs-melt),
            emitted-pulse radial energy profile
    loops   biaxial-core bookkeeping at snaps: per-cell min eigen-gap
            < thr mask -> 26-connected components; compact components
            disjoint from the edge = closed-loop candidates; COUNT
            (the two-loop conjecture read, thr sensitivity {0.06,
            0.09, 0.15})
    collect merge rows -> ../data/m5_21_6_all.json

Energy convention: E = h^3 sum_cells (u + V) as in the instrument;
KE = 0.5 h^3 sum |M_t|^2, so M_tt = -grad/h^3 with grad = dE/dM.
"""
from __future__ import annotations

import glob
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

_spec = importlib.util.spec_from_file_location(
    "ins", os.path.join(HERE, "m5_21_2b_a_instrument.py"))
INS = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(INS)


def rgrid(cfg):
    X, Y, Z = INS.coords(cfg["n"], cfg["h"])
    return np.sqrt(X * X + Y * Y + Z * Z)


# ================= relax (t-named wrapper over the instrument) =====
def relax(cfg):
    tag = cfg["tag"]
    M0 = INS.make_seed(cfg)
    n = cfg["n"]
    free = ~INS.pin_shell(n, cfg["h"]) if cfg["bc"] == "pinned" else \
        np.ones((n, n, n), dtype=bool)
    e0 = INS.e_parts(M0, cfg)
    snaps = (250, 500, 1000, 2000, 4000, 6000, 8000, 10000) \
        if cfg["snaps"] else ()
    M, states, info = INS.fire(M0, cfg, free, max_iter=cfg["maxit"],
                               log_every=500, snaps=snaps, tag=tag)
    e_u, e_d, e_v = INS.e_parts(M, cfg)
    row = {k: cfg[k] for k in ("seed", "term", "stencil", "eps", "n",
                               "L", "h", "delta", "bc", "maxit")}
    row.update({
        "tag": tag, "E_end": float(e_u + e_d + e_v),
        "E_u": float(e_u), "E_d": float(e_d), "E_v": float(e_v),
        "E_seed": float(sum(e0)),
        "r_half": INS.r_half(M, cfg),
        "retention": INS.retention(M, cfg),
        "min_gap_end": INS.min_gap(M),
        "stop": info["stop"], "trace": info["trace"][-6:],
        "wall_s": info["wall_s"]})
    os.makedirs(DATA, exist_ok=True)
    snap_arrays = {f"M_it{st['it']}": st["M"].astype(np.float32)
                   for st in states if st["it"] > 0}
    np.savez_compressed(
        os.path.join(DATA, f"m5_21_6_end_{tag}.npz"),
        M=M.astype(np.float32), delta=cfg["delta"], h=cfg["h"],
        n=cfg["n"], bc=cfg["bc"], term=cfg["term"],
        maxit=cfg["maxit"], **snap_arrays)
    with open(os.path.join(DATA, f"m5_21_6_row_{tag}.json"), "w") as f:
        json.dump(row, f, indent=1)
    print(json.dumps({k: row[k] for k in
                      ("tag", "E_end", "E_u", "E_v", "r_half",
                       "stop", "wall_s")}))
    return row


def load_end(tag):
    z = np.load(os.path.join(DATA, f"m5_21_6_end_{tag}.npz"))
    return z


# ================= kicks =================
def kick_random(M, cfg, eps, rng_seed=216, r_k=8.0):
    from scipy import ndimage
    rng = np.random.default_rng(rng_seed)
    n = cfg["n"]
    K = rng.standard_normal((n, n, n, 3, 3))
    K = ndimage.gaussian_filter(K, sigma=(2, 2, 2, 0, 0))
    K = 0.5 * (K + np.swapaxes(K, -1, -2))
    w = np.exp(-((rgrid(cfg) / r_k) ** 2))[..., None, None]
    K = K * w
    iso = (1.0 + cfg["delta"]) / 3.0 * np.eye(3)
    scale = np.sqrt(np.mean((M - iso) ** 2))
    K *= eps * scale / max(np.sqrt(np.mean(K ** 2)), 1e-300)
    out = M + K
    return 0.5 * (out + np.swapaxes(out, -1, -2))


def kick_rotate(M, cfg, theta_deg, r_k=8.0):
    """core-twist: rotate M by R(z, theta*w(r)) M R^T, w = gaussian."""
    th = np.deg2rad(theta_deg) * np.exp(-((rgrid(cfg) / r_k) ** 2))
    c, s = np.cos(th), np.sin(th)
    n = cfg["n"]
    R = np.zeros((n, n, n, 3, 3))
    R[..., 0, 0] = c
    R[..., 0, 1] = -s
    R[..., 1, 0] = s
    R[..., 1, 1] = c
    R[..., 2, 2] = 1.0
    out = np.einsum("...ij,...jk,...lk->...il", R, M, R)
    return 0.5 * (out + np.swapaxes(out, -1, -2))


def spec_core(M, cfg, r_c=8.0):
    lam = np.linalg.eigvalsh(M)
    m = rgrid(cfg) < r_c
    return np.sort(lam[m].mean(axis=0))[::-1]


# ================= gates (new code only; cap 3 tries) ==============
def leap_step(M, V, cfg, free, gam, dt):
    h3 = cfg["h"] ** 3
    F = -INS.grad(M, cfg) / h3
    V = (V + 0.5 * dt * F) * free
    V = V / (1.0 + 0.5 * dt * gam)
    M = M + dt * V
    F = -INS.grad(M, cfg) / h3
    V = (V + 0.5 * dt * F) * free
    V = V / (1.0 + 0.5 * dt * gam)
    return M, V


def sponge(cfg, g_max=0.5, r0=None, r1=None):
    r = rgrid(cfg)
    half = cfg["L"] / 2.0
    r0 = 0.65 * half if r0 is None else r0
    r1 = 0.98 * half if r1 is None else r1
    s = np.clip((r - r0) / max(r1 - r0, 1e-9), 0.0, 1.0)
    return (g_max * s * s)[..., None, None]


def e_tot(M, V, cfg):
    e_u, e_d, e_v = INS.e_parts(M, cfg)
    ke = 0.5 * cfg["h"] ** 3 * float(np.sum(V * V))
    return float(e_u + e_d + e_v), ke


def gates():
    cfg = INS.base_cfg(term="T2", n=32, L=48.0, bc="pinned")
    z = np.load(os.path.join(DATA, "m5_21_2b_end_i2_A_T2.npz"))
    M0 = z["M"].astype(float)
    free = (~INS.pin_shell(cfg["n"], cfg["h"]))[..., None, None] \
        .astype(float)
    out = {}
    # GK: the kick operators AS APPLIED (masked form) keep M exactly
    # symmetric and leave the pinned shell exactly untouched
    Mk1 = M0 + (kick_random(M0, cfg, 0.15) - M0) * free
    Mk2 = M0 + (kick_rotate(M0, cfg, 60.0) - M0) * free
    shell = 1.0 - free
    out["GK_sym_max"] = float(max(
        np.max(np.abs(Mk1 - np.swapaxes(Mk1, -1, -2))),
        np.max(np.abs(Mk2 - np.swapaxes(Mk2, -1, -2)))))
    out["GK_shell_leak"] = float(max(
        np.max(np.abs((Mk1 - M0) * shell)),
        np.max(np.abs((Mk2 - M0) * shell))))
    # GL1: gamma=0 conservation from a kicked state, dt ladder;
    # production dt = the rung with drift < 2.5e-3 over 200 steps
    ladder = {}
    dt_ok, noise = None, None
    for dt in (0.05, 0.025, 0.0125):
        M, V = Mk1.copy(), np.zeros_like(M0)
        E0, K0 = e_tot(M, V, cfg)
        worst, step_rise, bad = 0.0, 0.0, False
        prev = E0 + K0
        for it in range(200):
            M, V = leap_step(M, V, cfg, free, 0.0, dt)
            E, K = e_tot(M, V, cfg)
            if not np.isfinite(E + K):
                bad = True
                break
            worst = max(worst, abs((E + K) - (E0 + K0))
                        / max(abs(E0 + K0), 1e-300))
            step_rise = max(step_rise, (E + K) - prev)
            prev = E + K
        ladder[str(dt)] = None if bad else worst
        if not bad and worst < 2.5e-3 and dt_ok is None:
            dt_ok, noise = dt, step_rise
    out["GL1_ladder"] = ladder
    out["GL1_dt"] = dt_ok
    out["GL1_drift"] = None if dt_ok is None else ladder[str(dt_ok)]
    # GL2 (physical form): with the sponge on, (i) all finite,
    # (ii) net dissipation vs start, (iii) ends BELOW the
    # conservative run's E_tot, (iv) no single-step rise above 2x
    # the measured gamma=0 integrator noise floor
    gam = sponge(cfg)
    dtp = dt_ok or 0.0125
    M, V = Mk1.copy(), np.zeros_like(M0)
    E0s = sum(e_tot(M, V, cfg))
    Mc, Vc = Mk1.copy(), np.zeros_like(M0)
    prev, rise_max = E0s, 0.0
    fin = True
    for it in range(200):
        M, V = leap_step(M, V, cfg, free, gam, dtp)
        Mc, Vc = leap_step(Mc, Vc, cfg, free, 0.0, dtp)
        cur = sum(e_tot(M, V, cfg))
        if not np.isfinite(cur):
            fin = False
            break
        rise_max = max(rise_max, cur - prev)
        prev = cur
    E_end_sp = prev
    E_end_cons = sum(e_tot(Mc, Vc, cfg))
    out["GL2_E0"] = E0s
    out["GL2_E_end_sponge"] = E_end_sp
    out["GL2_E_end_conservative"] = E_end_cons
    out["GL2_step_rise_max"] = rise_max
    out["GL2_noise_floor"] = noise
    gl2 = bool(fin and E_end_sp < E0s and E_end_sp < E_end_cons
               and rise_max <= 2.0 * (noise or 0.0) + 1e-12)
    out["GL2_pass"] = gl2
    out["PASS"] = bool(out["GK_sym_max"] == 0.0
                       and out["GK_shell_leak"] == 0.0
                       and dt_ok is not None and gl2)
    with open(os.path.join(DATA, "m5_21_6_gates.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out))
    return out


# ================= P1: kick ladder (Arena 1, pinned) ===============
def p1(start_tag, maxit=4000):
    z = load_end(start_tag)
    cfg = INS.base_cfg(term="T2", n=int(z["n"]), L=48.0, bc="pinned")
    M0 = z["M"].astype(float)
    free = ~INS.pin_shell(cfg["n"], cfg["h"])
    e_u, e_d, e_v = INS.e_parts(M0, cfg)
    E_start = float(e_u + e_d + e_v)
    s0 = spec_core(M0, cfg)
    rows = []
    kicks = ([("K1", eps) for eps in (0.05, 0.15, 0.4)]
             + [("K2", th) for th in (30.0, 60.0, 90.0)])
    for fam, amp in kicks:
        Mk = kick_random(M0, cfg, amp) if fam == "K1" else \
            kick_rotate(M0, cfg, amp)
        Mk = M0 + (Mk - M0) * free[..., None, None]
        tagk = f"{start_tag}_{fam}_{amp:g}"
        M, _, info = INS.fire(Mk, cfg, free, max_iter=maxit,
                              log_every=1000, tag=tagk)
        e_u, e_d, e_v = INS.e_parts(M, cfg)
        E = float(e_u + e_d + e_v)
        s1 = spec_core(M, cfg)
        dE = (E - E_start) / max(abs(E_start), 1e-300)
        ds = float(np.max(np.abs(s1 - s0)))
        verdict = "returned" if (abs(dE) < 0.01 and ds < 0.02) else \
            ("moved_lower" if dE < -0.01 else "moved")
        ek_u, ek_d, ek_v = INS.e_parts(Mk, cfg)
        rows.append({"start": start_tag, "family": fam, "amp": amp,
                     "E_kicked": float(ek_u + ek_d + ek_v),
                     "E_end": E, "dE_rel": float(dE),
                     "spec_shift": ds, "verdict": verdict,
                     "stop": info["stop"],
                     "wall_s": info["wall_s"]})
        print(json.dumps(rows[-1]), flush=True)
        np.savez_compressed(
            os.path.join(DATA, f"m5_21_6_end_{tagk}.npz"),
            M=M.astype(np.float32), delta=cfg["delta"], h=cfg["h"],
            n=cfg["n"], bc="pinned", term="T2", maxit=maxit)
    with open(os.path.join(DATA,
                           f"m5_21_6_p1_{start_tag}.json"), "w") as f:
        json.dump({"E_start": E_start,
                   "spec_core_start": [float(x) for x in s0],
                   "rows": rows}, f, indent=1)
    return rows


# ================= evolve (P3, dynamics-grade) =====================
def evolve(start_tag, steps=4000, dt=None, snap_every=250,
           out_tag=None, kick=None, snap=None):
    z = load_end(start_tag)
    n = int(z["n"])
    bc = str(z["bc"])
    cfg = INS.base_cfg(term="T2", n=n, L=48.0, bc=bc)
    M = z[snap or "M"].astype(float)
    if kick:
        fam, amp = kick.split(":")
        M = kick_random(M, cfg, float(amp)) if fam == "K1" else \
            kick_rotate(M, cfg, float(amp))
    free_b = ~INS.pin_shell(n, cfg["h"]) if bc == "pinned" else \
        np.ones((n, n, n), dtype=bool)
    free = free_b[..., None, None].astype(float)
    if bc == "pinned":
        M_base = z[snap or "M"].astype(float)
        M = M_base + (M - M_base) * free
    gam = sponge(cfg)
    if dt is None:
        gts = json.load(open(os.path.join(DATA,
                                          "m5_21_6_gates.json")))
        dt = gts["GL1_dt"]
    V = np.zeros_like(M)
    r = rgrid(cfg)
    half = cfg["L"] / 2.0
    m_core = r < 10.0
    m_mid = (r >= 10.0) & (r < 0.65 * half)
    lam0 = np.linalg.eigvalsh(M)
    vecs0 = np.linalg.eigh(M)[1][..., -1]
    vac = np.sort(np.array([1.0, cfg["delta"], 0.0]))
    tag = out_tag or f"ev_{start_tag}"
    hist, snaps = [], {}
    absorbed = 0.0
    h3 = cfg["h"] ** 3
    t0 = time.time()
    for it in range(1, steps + 1):
        absorbed += float(np.sum(gam * V * V)) * h3 * dt
        M, V = leap_step(M, V, cfg, free, gam, dt)
        if it % snap_every == 0 or it == steps:
            e_u, e_d, e_v = INS.e_parts(M, cfg)
            dens = INS.v_density(M, cfg)
            lam = np.linalg.eigvalsh(M)
            spec_dev = float(np.mean(
                np.abs(lam[m_core] - vac[None, :])))
            v1 = np.linalg.eigh(M)[1][..., -1]
            dots = np.abs(np.sum(v1 * vecs0, axis=-1))
            rot = float(np.mean(np.arccos(np.clip(
                dots[m_core], 0, 1))) * 180 / np.pi)
            ke = 0.5 * h3 * float(np.sum(V * V))
            row = {"it": it, "t": it * dt,
                   "E": float(e_u + e_d + e_v), "KE": ke,
                   "absorbed": absorbed,
                   "spec_dev_core": spec_dev,
                   "rot_core_deg": rot}
            hist.append(row)
            print(f"  {tag} it {it:6d} E {row['E']:10.5f} "
                  f"KE {ke:9.5f} abs {absorbed:9.5f} "
                  f"specdev {spec_dev:.4f} rot {rot:6.2f}deg "
                  f"[{time.time() - t0:.0f}s]", flush=True)
            if it % (snap_every * 4) == 0 or it == steps:
                snaps[f"M_it{it}"] = M.astype(np.float32)
    np.savez_compressed(
        os.path.join(DATA, f"m5_21_6_ev_{tag}.npz"),
        M=M.astype(np.float32), delta=cfg["delta"], h=cfg["h"], n=n,
        bc=bc, dt=dt, term="T2", maxit=steps, **snaps)
    with open(os.path.join(DATA, f"m5_21_6_ev_{tag}.json"), "w") as f:
        json.dump({"start": start_tag, "kick": kick, "dt": dt,
                   "steps": steps, "hist": hist,
                   "wall_s": time.time() - t0}, f, indent=1)
    return hist


# ================= loops (topology bookkeeping) ====================
def loop_read(M, cfg, thr):
    from scipy import ndimage
    lam = np.linalg.eigvalsh(M)
    gap = np.minimum(lam[..., 1] - lam[..., 0],
                     lam[..., 2] - lam[..., 1])
    r = rgrid(cfg)
    half = cfg["L"] / 2.0
    mask = (gap < thr) & (r < 0.65 * half)
    lab, ncomp = ndimage.label(mask, structure=np.ones((3, 3, 3)))
    comps = []
    for k in range(1, ncomp + 1):
        sel = lab == k
        size = int(np.sum(sel))
        if size < 3:
            continue
        rc = float(np.mean(r[sel]))
        touches = bool(np.any(r[sel] > 0.62 * half))
        comps.append({"size": size, "r_centroid": rc,
                      "touches_edge_zone": touches})
    comps.sort(key=lambda c: -c["size"])
    return {"thr": thr, "n_components": len(comps),
            "n_compact": sum(1 for c in comps
                             if not c["touches_edge_zone"]),
            "components": comps[:12]}


def loops(npz_name):
    z = np.load(os.path.join(DATA, npz_name))
    n = int(z["n"])
    bc = str(z["bc"]) if "bc" in z else "pinned"
    cfg = INS.base_cfg(term="T2", n=n, L=48.0, bc=bc)
    out = {"file": npz_name, "snaps": {}}
    keys = ["M"] + [k for k in z.files if k.startswith("M_it")]
    for k in keys:
        M = z[k].astype(float)
        out["snaps"][k] = [loop_read(M, cfg, thr)
                           for thr in (0.06, 0.09, 0.15)]
        r9 = out["snaps"][k][1]
        print(f"  {npz_name}:{k} thr0.09 comps "
              f"{r9['n_components']} compact {r9['n_compact']}",
              flush=True)
    on = npz_name.replace(".npz", "_loops.json") \
        .replace("m5_21_6_", "m5_21_6_lp_")
    with open(os.path.join(DATA, on), "w") as f:
        json.dump(out, f, indent=1)
    return out


# ================= collect =================
def collect():
    rows = []
    for p in sorted(glob.glob(os.path.join(DATA,
                                           "m5_21_6_row_*.json"))):
        rows.append(json.load(open(p)))
    allp = {"rows": rows}
    for extra in ("gates", ):
        p = os.path.join(DATA, f"m5_21_6_{extra}.json")
        if os.path.exists(p):
            allp[extra] = json.load(open(p))
    for p in sorted(glob.glob(os.path.join(DATA,
                                           "m5_21_6_p1_*.json"))):
        allp[os.path.basename(p)[:-5]] = json.load(open(p))
    with open(os.path.join(DATA, "m5_21_6_all.json"), "w") as f:
        json.dump(allp, f, indent=1)
    print(f"collected {len(rows)} rows")


# ================= CLI =================
if __name__ == "__main__":
    argv = sys.argv[1:]
    mode = argv[0] if argv else "status"
    if mode == "relax":
        kw = INS.parse_kv(argv[1:])
        cfg = INS.base_cfg(term="T2", **kw)
        relax(cfg)
    elif mode == "gates":
        gates()
    elif mode == "p1":
        p1(argv[1], maxit=int(argv[2]) if len(argv) > 2 else 4000)
    elif mode == "evolve":
        kw = dict(a.split("=", 1) for a in argv[2:])
        evolve(argv[1],
               steps=int(kw.get("steps", 4000)),
               dt=float(kw["dt"]) if "dt" in kw else None,
               snap_every=int(kw.get("snap_every", 250)),
               out_tag=kw.get("out_tag"),
               kick=kw.get("kick"),
               snap=kw.get("snap"))
    elif mode == "loops":
        loops(argv[1])
    elif mode == "collect":
        collect()
    else:
        for p in sorted(glob.glob(os.path.join(
                DATA, "m5_21_6_row_*.json"))):
            r = json.load(open(p))
            print(f"{r['tag']:24s} E {r['E_end']:10.4f} "
                  f"stop {r['stop']}")
