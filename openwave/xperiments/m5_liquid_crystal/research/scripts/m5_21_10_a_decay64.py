"""M5.21.10: the decay-grade extension at the n = 64 free arena.

Thin wrapper over the certified M5.21.6 instrument (imported, not
copied): the physics (leap_step, sponge, kicks, loop_read, spec_core,
FIRE via the M5.21.2b stack) is the audited original. What this file
fixes/adds:

    cfg from npz    D.evolve/p1 hardcode L = 48.0; the f64 arena is
                    n = 64, L = 64, h = 1.0, so every mode here
                    rebuilds cfg as L = n * h from the npz record
    gates           GL1 gamma=0 E_tot conservation dt ladder and GL2
                    sponge monotone drain, run AT n = 64 free on a
                    kicked seed state (the M5.21.6 gates certified
                    n = 32 pinned only); GK kick symmetry at n = 64
    evolve          damped-wave dynamics identical to D.evolve but
                    with the corrected cfg, m5_21_10_ naming, and
                    snapshots every snap_every*2 (kinematics grade)
    loops           the M5.21.6 loop census (26-connected biaxial
                    cores, thr ladder {0.06, 0.09, 0.15}) over ALL
                    snapshots of an evolution npz, with full 3D
                    centroids added per component
    kin             the NEW release-kinematics read: match components
                    across snapshots (nearest-centroid), report per-
                    track velocity vector, speed, ejection direction

Convention notes: energies as in the instrument (E = h^3 sum(u + V));
KE = 0.5 h^3 sum |M_t|^2. Free arena: no pinned shell; the sponge
(outer-shell damping ramp) is the only boundary treatment.
"""
from __future__ import annotations

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

_spec = importlib.util.spec_from_file_location(
    "dec", os.path.join(HERE, "m5_21_6_a_decay.py"))
D = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(D)
INS = D.INS


def cfg_from_npz(z):
    n, h = int(z["n"]), float(z["h"])
    return INS.base_cfg(term="T2", n=n, L=n * h, bc=str(z["bc"]))


def load_end(tag):
    return np.load(os.path.join(DATA, f"m5_21_6_end_{tag}.npz"))


# ================= gates at n = 64 free ============================
def gates():
    cfg = INS.base_cfg(term="T2", seed="A", n=64, L=64.0, bc="free")
    assert cfg["h"] == 1.0
    M0 = INS.make_seed(cfg)
    free = np.ones((64, 64, 64, 1, 1), dtype=float)
    out = {"arena": {"n": 64, "L": 64.0, "h": 1.0, "bc": "free"}}
    # GK: kick operators keep M exactly symmetric at n = 64
    Mk1 = D.kick_random(M0, cfg, 0.15)
    Mk2 = D.kick_rotate(M0, cfg, 60.0)
    out["GK_sym_max"] = float(max(
        np.max(np.abs(Mk1 - np.swapaxes(Mk1, -1, -2))),
        np.max(np.abs(Mk2 - np.swapaxes(Mk2, -1, -2)))))
    # GL1: gamma=0 conservation from the kicked seed, dt ladder
    ladder, dt_ok = {}, None
    for dt in (0.05, 0.025, 0.0125):
        M, V = Mk1.copy(), np.zeros_like(M0)
        E0, K0 = D.e_tot(M, V, cfg)
        worst, bad = 0.0, False
        prev = E0 + K0
        for it in range(200):
            M, V = D.leap_step(M, V, cfg, free, 0.0, dt)
            E, K = D.e_tot(M, V, cfg)
            if not np.isfinite(E + K):
                bad = True
                break
            worst = max(worst, abs(E + K - (E0 + K0)) /
                        max(abs(E0 + K0), 1e-300))
            prev = E + K
        ladder[f"dt{dt:g}"] = None if bad else float(worst)
        if not bad and worst < 2.5e-3 and dt_ok is None:
            dt_ok = dt
        print(f"GL1 dt={dt:g} drift={ladder[f'dt{dt:g}']}", flush=True)
    out["GL1_ladder"] = ladder
    out["GL1_dt"] = dt_ok
    # GL2: sponge on -> E_tot + absorbed constant, E_tot monotone down
    gam = D.sponge(cfg)
    M, V = Mk1.copy(), np.zeros_like(M0)
    h3 = cfg["h"] ** 3
    absorbed, prev, mono = 0.0, None, True
    for it in range(200):
        absorbed += float(np.sum(gam * V * V)) * h3 * 0.025
        M, V = D.leap_step(M, V, cfg, free, gam, 0.025)
        if it % 20 == 19:
            E, K = D.e_tot(M, V, cfg)
            tot = E + K
            if prev is not None and tot > prev + 1e-9:
                mono = False
            prev = tot
    out["GL2_monotone"] = bool(mono)
    out["GL2_absorbed_200"] = float(absorbed)
    with open(os.path.join(DATA, "m5_21_10_gates.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out))


# ================= evolve (corrected cfg, kinematics grade) ========
def evolve(start_tag, steps=6000, dt=None, snap_every=200,
           out_tag=None, snap=None):
    z = load_end(start_tag)
    cfg = cfg_from_npz(z)
    n = cfg["n"]
    assert cfg["bc"] == "free", "f64 arena is free"
    M = z[snap or "M"].astype(float)
    free = np.ones((n, n, n, 1, 1), dtype=float)
    gam = D.sponge(cfg)
    if dt is None:
        dt = json.load(open(os.path.join(
            DATA, "m5_21_10_gates.json")))["GL1_dt"]
    V = np.zeros_like(M)
    r = D.rgrid(cfg)
    m_core = r < 10.0
    vecs0 = np.linalg.eigh(M)[1][..., -1]
    vac = np.sort(np.array([1.0, cfg["delta"], 0.0]))
    tag = out_tag or f"ev_{start_tag}"
    hist, snaps = [], {}
    absorbed = 0.0
    h3 = cfg["h"] ** 3
    t0 = time.time()
    for it in range(1, steps + 1):
        absorbed += float(np.sum(gam * V * V)) * h3 * dt
        M, V = D.leap_step(M, V, cfg, free, gam, dt)
        if it % snap_every == 0 or it == steps:
            e_u, e_d, e_v = INS.e_parts(M, cfg)
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
            if it % (snap_every * 2) == 0 or it == steps:
                snaps[f"M_it{it}"] = M.astype(np.float32)
    np.savez_compressed(
        os.path.join(DATA, f"m5_21_10_ev_{tag}.npz"),
        M=M.astype(np.float32), delta=cfg["delta"], h=cfg["h"], n=n,
        bc=cfg["bc"], dt=dt, term="T2", maxit=steps, **snaps)
    with open(os.path.join(DATA,
                           f"m5_21_10_ev_{tag}.json"), "w") as f:
        json.dump({"start": start_tag, "dt": dt, "steps": steps,
                   "snap_every": snap_every, "hist": hist,
                   "wall_s": time.time() - t0}, f, indent=1)
    return hist


# ================= loops over all snapshots ========================
def loop_read_xyz(M, cfg, thr):
    """D.loop_read + full 3D centroids per component."""
    from scipy import ndimage
    lam = np.linalg.eigvalsh(M)
    gap = np.minimum(lam[..., 1] - lam[..., 0],
                     lam[..., 2] - lam[..., 1])
    r = D.rgrid(cfg)
    half = cfg["L"] / 2.0
    X, Y, Z = INS.coords(cfg["n"], cfg["h"])
    mask = (gap < thr) & (r < 0.65 * half)
    lab, ncomp = ndimage.label(mask, structure=np.ones((3, 3, 3)))
    comps = []
    for k in range(1, ncomp + 1):
        sel = lab == k
        size = int(np.sum(sel))
        if size < 3:
            continue
        comps.append({
            "size": size,
            "r_centroid": float(np.mean(r[sel])),
            "xyz": [float(np.mean(X[sel])), float(np.mean(Y[sel])),
                    float(np.mean(Z[sel]))],
            "touches_edge_zone": bool(np.any(r[sel] > 0.62 * half))})
    comps.sort(key=lambda c: -c["size"])
    return {"thr": thr, "n_components": len(comps),
            "n_compact": sum(1 for c in comps
                             if not c["touches_edge_zone"]),
            "components": comps[:12]}


def loops(ev_tag):
    z = np.load(os.path.join(DATA, f"m5_21_10_ev_{ev_tag}.npz"))
    cfg = cfg_from_npz(z)
    keys = sorted([k for k in z.files if k.startswith("M_it")],
                  key=lambda k: int(k[4:]))
    dt = float(z["dt"])
    out = {"ev": ev_tag, "dt": dt, "snaps": {}}
    for k in keys:
        M = z[k].astype(float)
        out["snaps"][k] = {
            "t": int(k[4:]) * dt,
            "thr": [loop_read_xyz(M, cfg, thr)
                    for thr in (0.06, 0.09, 0.15)]}
        n9 = out["snaps"][k]["thr"][1]
        print(f"{ev_tag} {k} t={out['snaps'][k]['t']:.1f} "
              f"thr0.09: {n9['n_components']} comps "
              f"({n9['n_compact']} compact)", flush=True)
    with open(os.path.join(DATA,
                           f"m5_21_10_lp_{ev_tag}.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


# ================= kin: release kinematics =========================
def kin(ev_tag, thr_idx=1, match_d=6.0):
    """component tracks across snapshots at one thr rung
    (default 0.09): greedy nearest-centroid matching within match_d;
    per-track velocity by finite difference on the centroid."""
    lp = json.load(open(os.path.join(
        DATA, f"m5_21_10_lp_{ev_tag}.json")))
    keys = sorted(lp["snaps"].keys(), key=lambda k: int(k[4:]))
    tracks = []           # each: list of (t, xyz, size, touches)
    for k in keys:
        s = lp["snaps"][k]
        comps = s["thr"][thr_idx]["components"]
        t = s["t"]
        used = set()
        for tr in tracks:
            _, last_xyz, _, _ = tr[-1]
            best, bd = None, match_d
            for i, c in enumerate(comps):
                if i in used:
                    continue
                d = float(np.linalg.norm(
                    np.array(c["xyz"]) - np.array(last_xyz)))
                if d < bd:
                    best, bd = i, d
            if best is not None:
                c = comps[best]
                used.add(best)
                tr.append((t, c["xyz"], c["size"],
                           c["touches_edge_zone"]))
        for i, c in enumerate(comps):
            if i not in used:
                tracks.append([(t, c["xyz"], c["size"],
                                c["touches_edge_zone"])])
    rows = []
    for tr in tracks:
        if len(tr) < 3:
            continue
        ts = np.array([p[0] for p in tr])
        xyz = np.array([p[1] for p in tr])
        v = (xyz[-1] - xyz[0]) / max(ts[-1] - ts[0], 1e-300)
        sp = float(np.linalg.norm(v))
        rows.append({
            "t_first": float(ts[0]), "t_last": float(ts[-1]),
            "n_snaps": len(tr),
            "size_first": tr[0][2], "size_last": tr[-1][2],
            "xyz_first": [float(x) for x in xyz[0]],
            "xyz_last": [float(x) for x in xyz[-1]],
            "v": [float(x) for x in v], "speed": sp,
            "dir": [float(x) for x in v / max(sp, 1e-300)],
            "r_first": float(np.linalg.norm(xyz[0])),
            "r_last": float(np.linalg.norm(xyz[-1])),
            "departing": bool(np.linalg.norm(xyz[-1]) >
                              np.linalg.norm(xyz[0]) + 2.0),
            "ever_touches_edge_zone": bool(any(p[3] for p in tr))})
    rows.sort(key=lambda r: -r["speed"])
    out = {"ev": ev_tag, "thr": lp["snaps"][keys[0]]["thr"][thr_idx]
           ["thr"], "match_d": match_d, "n_tracks": len(rows),
           "tracks": rows}
    with open(os.path.join(DATA,
                           f"m5_21_10_kin_{ev_tag}.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"ev": ev_tag, "n_tracks": len(rows),
                      "departing": sum(1 for r in rows
                                       if r["departing"])}))
    return out


# ================= CLI =================
if __name__ == "__main__":
    argv = sys.argv[1:]
    mode = argv[0] if argv else "gates"
    if mode == "gates":
        gates()
    elif mode == "evolve":
        kw = dict(a.split("=", 1) for a in argv[2:])
        evolve(argv[1],
               steps=int(kw.get("steps", 6000)),
               dt=float(kw["dt"]) if "dt" in kw else None,
               snap_every=int(kw.get("snap_every", 200)),
               out_tag=kw.get("out_tag"),
               snap=kw.get("snap"))
    elif mode == "loops":
        loops(argv[1])
    elif mode == "kin":
        kw = dict(a.split("=", 1) for a in argv[2:])
        kin(argv[1], thr_idx=int(kw.get("thr_idx", 1)),
            match_d=float(kw.get("match_d", 6.0)))
    else:
        print("modes: gates | evolve <tag> | loops <ev_tag> | "
              "kin <ev_tag>")
