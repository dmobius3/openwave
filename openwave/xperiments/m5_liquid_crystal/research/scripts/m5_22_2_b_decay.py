"""M5.22.2: beta decay of the kicked neutral states (rung 4).

The M5.21.6 energy-injection kick protocol (K1 random smooth core
kick, K2 core-twist rotation; GL/GK gates verified there) transplanted
to the two neutral ring-antiring states, against the author's expected
channels (2026-08-02 reply, m5_22_convo.md):

    neutron:   n -> p + e + nu (the neutrino "here topological
               vortex knot")
    bineutron: -> p + n + e, or two neutrons

Targets (the pre-registered set):
    P-1_plane_sc6_n32_pinned_d0.3          census neutron-analog
                                           (the ring-antiring pair)
    d2_s-0.5_s-0.5_a2_sc6_n32_d0.3 (m5_22_1)  the pp-control cousin
                                           (second neutral basin)

Per run: endpoint -> kick (masked to the free interior, shell exact)
-> damped leapfrog evolve (DT = 0.025 GL1-verified; sponge r0frac 0.8
so rings at r ~ 16 sit undamped) -> FIRE to the product basin.

Reads: energy budget (E_start, E_injected, absorbed-by-sponge, E_end);
ring ledger + per-fragment slab charges (KICK.endpoint_reads); the
BOTH-VARIANT div E moments (m5_22_2_a_dive: full = the arXiv
2108.07896 dual curvature, curv = long-axis field-line curvature).
Classification per endpoint: returned / rearranged / decayed
(ledger vs the start endpoint) with products by ring count + charge.
A non-decay is a reportable outcome, not a failure (the author's own
caveat: toy parameters + missing angular momenta).

Modes:
  decay   src=m5_22 tag=... kick=K1:0.4 [steps=3000 maxit=8000
          seed=216 r0frac=0.8]
  stage1  the pre-registered ladder: both targets x
          (K1:0.4, K1:0.8, K2:90), sequential, checkpoint per run
  spectrum src=... tag=... kick=K1:0.4 seeds=8 [steps=3000 maxit=6000]
          seed ensemble for the emitted-energy distribution
Out: ../data/m5_22_2_end_<tag>__<kick>_s<seed>.npz + row json
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
CKPT = os.path.join(HERE, "..", "checkpoints", "m5_22_2_progress.md")


def _load(name, fname):
    sp = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


KICK = _load("kick", "m5_22_1_a_kick.py")
DIVE = _load("dive", "m5_22_2_a_dive.py")
M56 = _load("m56", "m5_21_6_a_decay.py")
CEN = KICK.CEN
INS = KICK.INS
PAIR = KICK.PAIR
W2_T2 = KICK.W2_T2
DT = 0.025


def _ckpt(line):
    with open(CKPT, "a") as f:
        f.write(line.rstrip() + "\n")


def kick_ring(M, cfg, eps, ring_rho, ring_z, rng_seed=216, r_k=4.0):
    """K3: ring-localized random smooth kick. The K1 envelope
    exp(-(r/8)^2) is ORIGIN-centered, so the rings at r ~ 17 receive
    weight ~ 0.01 (measured in stage1): the decay targets were barely
    kicked. K3 centers the envelope on ONE ring torus,
    d_ring = sqrt((rho - ring_rho)^2 + (z - ring_z)^2): the
    convert-one-neutron probe."""
    from scipy import ndimage
    rng = np.random.default_rng(rng_seed)
    n = cfg["n"]
    K = rng.standard_normal((n, n, n, 3, 3))
    K = ndimage.gaussian_filter(K, sigma=(2, 2, 2, 0, 0))
    K = 0.5 * (K + np.swapaxes(K, -1, -2))
    X, Y, Z = INS.coords(n, cfg["h"])
    rho = np.sqrt(X * X + Y * Y)
    d_ring = np.sqrt((rho - ring_rho) ** 2 + (Z - ring_z) ** 2)
    w = np.exp(-((d_ring / r_k) ** 2))[..., None, None]
    K = K * w
    iso = (1.0 + cfg["delta"]) / 3.0 * np.eye(3)
    scale = np.sqrt(np.mean((M - iso) ** 2))
    K *= eps * scale / max(np.sqrt(np.mean(K ** 2)), 1e-300)
    out = M + K
    return 0.5 * (out + np.swapaxes(out, -1, -2))


def top_ring(M, cfg):
    """(rho_mean, z_centroid) of the upper (z > 0) ring at thr0.15."""
    rr = KICK.ring_read(M, cfg, 0.15)
    ups = [r for r in rr["rings"] if r["z_centroid"] > 0]
    if not ups:
        raise RuntimeError("no upper ring found at thr0.15")
    r = max(ups, key=lambda e: e["voxels"])
    return r["rho_mean"], r["z_centroid"]


def apply_kick(M0, cfg, kick_spec, seed):
    fam, amp = kick_spec.split(":")
    if fam == "K1":
        Mk = M56.kick_random(M0, cfg, float(amp), rng_seed=seed)
    elif fam == "K3":
        ring_rho, ring_z = top_ring(M0, cfg)
        Mk = kick_ring(M0, cfg, float(amp), ring_rho, ring_z,
                       rng_seed=seed)
    else:
        Mk = M56.kick_rotate(M0, cfg, float(amp))
    free = (~INS.pin_shell(cfg["n"], cfg["h"]))[..., None, None] \
        .astype(float)
    Mk = M0 + (Mk - M0) * free
    return 0.5 * (Mk + np.swapaxes(Mk, -1, -2))


def decay_run(src, tag, kick_spec, steps=3000, maxit=8000, seed=216,
              r0frac=0.8, light=False):
    z = np.load(os.path.join(DATA, f"{src}_end_{tag}.npz"))
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(z["n"]), delta=float(z["delta"]),
                       bc="pinned")
    n, h = cfg["n"], cfg["h"]
    h3 = h ** 3
    M0 = z["M"].astype(np.float64)
    e0 = float(sum(INS.e_parts(M0, cfg)))
    Mk = apply_kick(M0, cfg, kick_spec, seed)
    ek = float(sum(INS.e_parts(Mk, cfg)))
    free_b = ~INS.pin_shell(n, h)
    free = free_b[..., None, None].astype(float)
    gam = KICK.sponge(cfg, r0frac=r0frac)
    out_tag = f"{tag}__{kick_spec.replace(':', '')}_s{seed}"
    t0 = time.time()
    M, V = Mk, np.zeros_like(Mk)
    absorbed = 0.0
    hist = []
    for it in range(1, steps + 1):
        absorbed += float(np.sum(gam * V * V)) * h3 * DT
        M, V = KICK.leap_step(M, V, cfg, free, gam, DT)
        if it % 500 == 0 or it == steps:
            e_now = float(sum(INS.e_parts(M, cfg)))
            ke = 0.5 * h3 * float(np.sum(V * V))
            rr = KICK.ring_read(M, cfg, 0.09)
            row = {"it": it, "E": e_now, "KE": ke,
                   "absorbed": absorbed,
                   "n_rings_thr0.09": rr["n_rings"]}
            hist.append(row)
            print(f"  {out_tag} it {it:5d} E {e_now:9.4f} "
                  f"KE {ke:8.4f} abs {absorbed:8.4f} "
                  f"rings {rr['n_rings']} "
                  f"[{time.time() - t0:.0f}s]", flush=True)
    M, _, info = INS.fire(M, cfg, free_b, max_iter=maxit,
                          log_every=2000, tag=out_tag)
    row = {"src": src, "start_tag": tag, "kick": kick_spec,
           "seed": seed, "steps": steps, "maxit": maxit,
           "r0frac": r0frac, "dt": DT,
           "E_start": e0, "E_kicked": ek, "E_injected": ek - e0,
           "absorbed": absorbed, "stop": info["stop"],
           "evolve_hist": hist, "n": n, "delta": cfg["delta"]}
    row.update(KICK.endpoint_reads(M, cfg))
    row["dive"] = DIVE.both_variants(M, cfg)
    row["E_returned_delta"] = row["E_end"] - e0
    row["wall_s"] = time.time() - t0
    if not light:
        np.savez_compressed(
            os.path.join(DATA, f"m5_22_2_end_{out_tag}.npz"),
            M=M.astype(np.float32), M0=M0.astype(np.float32),
            delta=cfg["delta"], h=h, n=n)
    with open(os.path.join(DATA, f"m5_22_2_row_{out_tag}.json"),
              "w") as f:
        json.dump(row, f, indent=1)
    r9 = row["rings"]["thr0.09"]
    summary = {"out_tag": out_tag, "E_start": e0,
               "E_injected": round(ek - e0, 4),
               "E_end": row["E_end"],
               "dE_vs_start": round(row["E_returned_delta"], 4),
               "absorbed": round(absorbed, 4),
               "rings": r9["n_rings"],
               "ring_zs": [round(r["z_centroid"], 1)
                           for r in r9["rings"]],
               "q_far_full": row["dive"]["full"]["q_flux"]["half18"],
               "stop": info["stop"],
               "wall_s": round(row["wall_s"])}
    print(json.dumps(summary))
    return row, summary


STAGE1 = [
    ("m5_22", "P-1_plane_sc6_n32_pinned_d0.3"),
    ("m5_22_1", "d2_s-0.5_s-0.5_a2_sc6_n32_d0.3"),
]
# widened after the smoke test measured the injection scale
# (K1:0.4 injects E ~ 295 on E_start ~ 12.7, 23x): the ladder must
# bracket the decay barrier from below, not only blast the state
KICKS1 = ["K1:0.05", "K1:0.15", "K1:0.4", "K2:30", "K2:90"]


def stage1(_kw):
    _ckpt(f"\n## stage1 launched {time.strftime('%Y-%m-%d %H:%M')}")
    for src, tag in STAGE1:
        for ks in KICKS1:
            try:
                _, summary = decay_run(src, tag, ks)
                _ckpt(f"- ✅ `{summary['out_tag']}`: E_inj "
                      f"{summary['E_injected']}, dE_vs_start "
                      f"{summary['dE_vs_start']}, rings "
                      f"{summary['rings']} at {summary['ring_zs']}, "
                      f"absorbed {summary['absorbed']}, "
                      f"stop {summary['stop']}")
            except Exception as exc:  # keep the ladder going
                _ckpt(f"- ❌ {tag} {ks}: {exc!r}")
                print(f"FAILED {tag} {ks}: {exc!r}", flush=True)
    _ckpt(f"- stage1 complete {time.strftime('%Y-%m-%d %H:%M')}")
    print("STAGE1 COMPLETE", flush=True)


def extend_run(out_tag, maxit=12000):
    """converge an unfinished m5_22_2 endpoint by more FIRE."""
    z = np.load(os.path.join(DATA, f"m5_22_2_end_{out_tag}.npz"))
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(z["n"]), delta=float(z["delta"]),
                       bc="pinned")
    M = z["M"].astype(np.float64)
    free_b = ~INS.pin_shell(cfg["n"], cfg["h"])
    t0 = time.time()
    M, _, info = INS.fire(M, cfg, free_b, max_iter=maxit,
                          log_every=2000, tag=out_tag + "_ext")
    row = {"branch": "extend", "parent": out_tag, "maxit": maxit,
           "stop": info["stop"], "n": cfg["n"],
           "delta": cfg["delta"]}
    row.update(KICK.endpoint_reads(M, cfg))
    row["dive"] = DIVE.both_variants(M, cfg)
    row["wall_s"] = time.time() - t0
    np.savez_compressed(
        os.path.join(DATA, f"m5_22_2_end_{out_tag}_ext.npz"),
        M=M.astype(np.float32), delta=cfg["delta"], h=cfg["h"],
        n=cfg["n"])
    with open(os.path.join(DATA,
                           f"m5_22_2_row_{out_tag}_ext.json"),
              "w") as f:
        json.dump(row, f, indent=1)
    r15 = row["rings"]["thr0.15"]
    summary = {"out_tag": out_tag + "_ext", "E_end": row["E_end"],
               "stop": info["stop"], "rings": r15["n_rings"],
               "ring_zs": [round(r["z_centroid"], 1)
                           for r in r15["rings"]],
               "q_upper": round(row["charges"]["q_upper_half"], 3),
               "q_lower": round(row["charges"]["q_lower_half"], 3),
               "wall_s": round(row["wall_s"])}
    print(json.dumps(summary))
    return row, summary


def stage2(_kw):
    """the exhaustion follow-ups: converge the two K1:0.4 blasts,
    then the ring-localized K3 ladder on both targets."""
    _ckpt(f"\n## stage2 launched {time.strftime('%Y-%m-%d %H:%M')}")
    for out_tag in [
            "P-1_plane_sc6_n32_pinned_d0.3__K10.4_s216",
            "d2_s-0.5_s-0.5_a2_sc6_n32_d0.3__K10.4_s216"]:
        try:
            _, s = extend_run(out_tag)
            _ckpt(f"- ✅ `{s['out_tag']}`: E_end {s['E_end']:.4f}, "
                  f"rings {s['rings']} at {s['ring_zs']}, q_upper "
                  f"{s['q_upper']}, q_lower {s['q_lower']}, "
                  f"stop {s['stop']}")
        except Exception as exc:
            _ckpt(f"- ❌ extend {out_tag}: {exc!r}")
            print(f"FAILED extend {out_tag}: {exc!r}", flush=True)
    for src, tag in STAGE1:
        # K3 injects ~7x K1 at equal eps (the torus envelope sits on
        # the high-gradient region; K3:0.4 -> E_inj ~ 671 measured),
        # so the ladder reaches DOWN to 0.02 to bracket the barrier
        for ks in ["K3:0.02", "K3:0.05", "K3:0.15", "K3:0.4"]:
            try:
                _, s = decay_run(src, tag, ks)
                _ckpt(f"- ✅ `{s['out_tag']}`: E_inj "
                      f"{s['E_injected']}, dE_vs_start "
                      f"{s['dE_vs_start']}, rings {s['rings']} at "
                      f"{s['ring_zs']}, absorbed {s['absorbed']}, "
                      f"stop {s['stop']}")
            except Exception as exc:
                _ckpt(f"- ❌ {tag} {ks}: {exc!r}")
                print(f"FAILED {tag} {ks}: {exc!r}", flush=True)
    _ckpt(f"- stage2 complete {time.strftime('%Y-%m-%d %H:%M')}")
    print("STAGE2 COMPLETE", flush=True)


def spectrum(kw):
    src = kw.get("src", "m5_22")
    tag = kw["tag"]
    ks = kw.get("kick", "K1:0.4")
    n_seeds = int(kw.get("seeds", 8))
    steps = int(kw.get("steps", 3000))
    maxit = int(kw.get("maxit", 6000))
    rows = []
    _ckpt(f"\n## spectrum {tag} {ks} x{n_seeds} "
          f"{time.strftime('%Y-%m-%d %H:%M')}")
    for k in range(n_seeds):
        seed = 1000 + 37 * k
        _, summary = decay_run(src, tag, ks, steps=steps,
                               maxit=maxit, seed=seed, light=True)
        rows.append(summary)
        _ckpt(f"- seed {seed}: E_end {summary['E_end']:.4f}, "
              f"absorbed {summary['absorbed']}, "
              f"rings {summary['rings']}")
    with open(os.path.join(
            DATA, f"m5_22_2_spectrum_{tag}__"
            f"{ks.replace(':', '')}.json"), "w") as f:
        json.dump(rows, f, indent=1)
    print("SPECTRUM COMPLETE", flush=True)


if __name__ == "__main__":
    ARGV = sys.argv[1:]
    mode = ARGV[0]
    kw = dict(a.split("=", 1) for a in ARGV[1:])
    if mode == "decay":
        decay_run(kw.get("src", "m5_22"), kw["tag"], kw["kick"],
                  steps=int(kw.get("steps", 3000)),
                  maxit=int(kw.get("maxit", 8000)),
                  seed=int(kw.get("seed", 216)),
                  r0frac=float(kw.get("r0frac", 0.8)))
    elif mode == "stage1":
        stage1(kw)
    elif mode == "stage2":
        stage2(kw)
    elif mode == "extend":
        extend_run(kw["tag"], maxit=int(kw.get("maxit", 12000)))
    elif mode == "spectrum":
        spectrum(kw)
