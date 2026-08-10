"""M5.22.4 main course: the omega-twist ladder on the baryon states.

The author's directive (2026-08-02, re-endorsed 2026-08-06): baryon
minima "should not be static, but very dynamical", suspected "mainly
by twists of the long axis", approximated as "twists of constant
frequency omega". The M5.21.3 machinery transplanted verbatim (the
audited 4D stack: eta = diag(-1,1,1,1), trace-target V4, sym stencil,
FIRE): E*(omega) = min_M [E_static(M) + omega^2 kin(M; a0)], a0 =
clock_local (rotation about the local long axis = the named twist).

The four targets (M5.22 census + M5.22.1), all n = 32, L = 48,
delta = 0.3, relaxed under the census T2 term; the 4D lift relaxes
under the trace-target V4 (the M5.21.3 precedent), so P1 carries a
RING-SURVIVAL read (census charge profile + core ledger on the
spatial block) - a deformed or lost ring is a first-class outcome,
not a silent assumption.

Envelope: renv = 18 (NOT the M5.21.3 electron's 10): the ring cores
sit at r ~ 16-17 (the M5.22.2 K1 lesson: an origin-centered envelope
gives the rings ~1% weight). Logged in the task planning.

Modes:
  p1 <key> [s] [maxit]   static 4D lift + ring survival (default s=+1,
                         maxit 6000)
  p2 <key>               generator catalog: kin signs + twist channel
  p3 <key> [gen] [omegas] omega-ladder, warm-started (default
                         clock_local, 0.05,0.1,0.2,0.4,0.8, 3000/rung)
  ctrl <key> [maxit]     matched-depth static control (default 15000)
  collect                merge rows -> m5_22_4_all.json
Out: ../data/m5_22_4_row_*.json + m5_22_4_p1_<key>[_s-1].npz
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def _load(name, fname):
    sp = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


G4 = _load("g4", "m5_21_3_a_4d.py")
CEN = _load("cen", "m5_22_b_census.py")

RENV = 18.0

STATES = {
    "prot": ("m5_22", "P-0.5_plane_sc6_n32_pinned_d0.3",
             "census proton-analog"),
    "neut": ("m5_22", "P-1_plane_sc6_n32_pinned_d0.3",
             "census neutron-analog (ring-antiring)"),
    "d2": ("m5_22_1", "d2_s-0.5_s-0.5_a2_sc6_n32_d0.3",
           "second neutral basin (pp-control cousin)"),
    "deut": ("m5_22_1", "dn_-0.5at-2_+0.5at+0_-0.5at+2_n32_d0.3",
             "deuteron candidate"),
}


def load3(key):
    src, tag, note = STATES[key]
    z = np.load(os.path.join(DATA, f"{src}_end_{tag}.npz"))
    return z["M"].astype(np.float64), tag, note


def ring_reads(M4, cfg4):
    """census instruments on the spatial block: charge profile
    (q_far + core/shell) + core ledger (ring components)."""
    M3 = M4[..., 1:4, 1:4]
    c3 = CEN.INS.base_cfg(term="T2", stencil="sym", eps=0.0,
                          w2=CEN.W2_T2, n=cfg4["n"],
                          delta=cfg4["delta"], bc="pinned")
    prof = CEN.charge_profile(M3, c3)
    led = CEN.core_ledger(M3, c3)
    return {"q_far": prof["q_far"],
            "profile_tail": prof["profile"][-3:],
            "min_gap": led["min_gap"],
            "core_ledger": led}


def p1(key, s=1.0, maxit=6000):
    M3, tag, note = load3(key)
    sfx = "" if s > 0 else "_s-1"
    cfg = G4.base_cfg(s=float(s), tag=f"p1_{key}{sfx}")
    M0 = G4.embed34(M3, cfg)
    cfg["renv"] = RENV
    free = ~G4.pin_shell(cfg["n"], cfg["h"])
    e0u, e0v = G4.e_parts(M0, cfg)
    ring0 = ring_reads(M0, cfg)
    M, info = G4.fire(M0, cfg, free, int(maxit), tag=cfg["tag"])
    eu, ev = G4.e_parts(M, cfg)
    ring1 = ring_reads(M, cfg)
    offb = float(np.max(np.abs(M[..., 0, 1:4])))
    row = {"tag": cfg["tag"], "state": tag, "note": note,
           "s": cfg["s"], "E_seed": float(e0u + e0v),
           "E_end": float(eu + ev), "E_u": float(eu),
           "E_v": float(ev), "offblock_max": offb,
           "ring_seed": ring0, "ring_end": ring1,
           "stop": info["stop"], "trace": info["trace"][-3:],
           "wall_s": info["wall_s"]}
    np.savez_compressed(
        os.path.join(DATA, f"m5_22_4_p1_{key}{sfx}.npz"),
        M=M.astype(np.float32), s=cfg["s"], delta=cfg["delta"],
        h=cfg["h"], n=cfg["n"])
    with open(os.path.join(DATA,
                           f"m5_22_4_row_p1_{key}{sfx}.json"),
              "w") as f:
        json.dump(row, f, indent=1)
    print(json.dumps({k: row[k] for k in
                      ("tag", "E_seed", "E_end", "offblock_max",
                       "stop")}
                     | {"q_far_seed": ring0["q_far"],
                        "q_far_end": ring1["q_far"]}))
    return row


def load_p1(key, sfx=""):
    Z = np.load(os.path.join(DATA, f"m5_22_4_p1_{key}{sfx}.npz"))
    cfg = G4.base_cfg(s=float(Z["s"]), tag=f"{key}{sfx}")
    cfg["renv"] = RENV
    return Z["M"].astype(np.float64), cfg


def p2(key, sfx=""):
    M, cfg = load_p1(key, sfx)
    a0s = G4.gen_catalog(cfg, M)
    rows = {}
    for nm, a0 in a0s.items():
        k = G4.kin_of(M, a0, cfg)
        tw = G4.twist_read(M, a0, cfg)
        rows[nm] = {"kin": float(k), "twist": tw}
        print(f"  {key}{sfx} {nm:12s} kin {k:+.6e} "
              f"b {tw['b']:+.3e} k* {tw['k_star']:+.4f} "
              f"dE(k*) {tw['dE_at_kstar']:+.3e}", flush=True)
    out = {"tag": f"p2_{key}{sfx}", "rows": rows,
           "min_kin": float(min(r["kin"] for r in rows.values())),
           "argmin": min(rows, key=lambda nm: rows[nm]["kin"])}
    with open(os.path.join(DATA, f"m5_22_4_row_p2_{key}{sfx}.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"tag": out["tag"], "min_kin": out["min_kin"],
                      "argmin": out["argmin"]}))
    return out


def p3(key, gen_name="clock_local",
       omegas=(0.05, 0.1, 0.2, 0.4, 0.8), maxit=3000):
    M, cfg = load_p1(key)
    free = ~G4.pin_shell(cfg["n"], cfg["h"])
    E_stat = float(G4.e_total(M, cfg))
    ladder = [{"omega": 0.0, "E": E_stat, "kin": None,
               "stop": "ref"}]
    dive_floor = E_stat - 10.0 * abs(E_stat) - 100.0
    for om in omegas:
        a0 = G4.gen_catalog(cfg, M)[gen_name]
        k0 = float(G4.kin_of(M, a0, cfg))
        tg = f"p3_{key}_{gen_name}_w{om:g}"
        M, info = G4.fire(M, cfg, free, int(maxit), a0=a0,
                          omega=om, tag=tg, dive_floor=dive_floor)
        eu, ev = G4.e_parts(M, cfg)
        kin_end = float(G4.kin_of(M, a0, cfg))
        E = float(eu + ev + om ** 2 * kin_end)
        ring = ring_reads(M, cfg)
        ladder.append({"omega": om, "E": E, "E_u": float(eu),
                       "E_v": float(ev), "kin": kin_end,
                       "kin_pre": k0, "q_far": ring["q_far"],
                       "stop": info["stop"]})
        print(json.dumps(ladder[-1]), flush=True)
        if info["stop"] in ("dive", "non-finite"):
            break
    out = {"tag": f"p3_{key}_{gen_name}", "E_static": E_stat,
           "ladder": ladder}
    with open(os.path.join(
            DATA, f"m5_22_4_row_p3_{key}_{gen_name}.json"),
            "w") as f:
        json.dump(out, f, indent=1)
    np.savez_compressed(
        os.path.join(DATA, f"m5_22_4_p3_{key}_{gen_name}.npz"),
        M=M.astype(np.float32), s=cfg["s"], delta=cfg["delta"],
        h=cfg["h"], n=cfg["n"])
    return out


def ctrl(key, maxit=15000):
    M, cfg = load_p1(key)
    free = ~G4.pin_shell(cfg["n"], cfg["h"])
    M2, info = G4.fire(M, cfg, free, int(maxit),
                       tag=f"ctrl_{key}")
    eu, ev = G4.e_parts(M2, cfg)
    ring = ring_reads(M2, cfg)
    row = {"tag": f"ctrl_{key}", "maxit": int(maxit),
           "E_end": float(eu + ev), "E_u": float(eu),
           "E_v": float(ev), "q_far": ring["q_far"],
           "stop": info["stop"], "trace": info["trace"][-3:]}
    with open(os.path.join(DATA, f"m5_22_4_row_ctrl_{key}.json"),
              "w") as f:
        json.dump(row, f, indent=1)
    print(json.dumps({k: row[k] for k in ("tag", "E_end", "stop",
                                          "q_far")}))
    return row


def collect():
    import glob
    rows = {}
    for p in sorted(glob.glob(
            os.path.join(DATA, "m5_22_4_row_*.json"))):
        key = os.path.basename(p)[len("m5_22_4_row_"):-len(".json")]
        with open(p) as f:
            rows[key] = json.load(f)
    with open(os.path.join(DATA, "m5_22_4_all.json"), "w") as f:
        json.dump(rows, f, indent=1)
    print(f"collected {len(rows)} rows")


if __name__ == "__main__":
    ARGV = sys.argv[1:]
    mode = ARGV[0]
    if mode == "p1":
        p1(ARGV[1], *(float(a) for a in ARGV[2:3]),
           **({"maxit": int(ARGV[3])} if len(ARGV) > 3 else {}))
    elif mode == "p2":
        p2(ARGV[1], *(a for a in ARGV[2:3]))
    elif mode == "p3":
        om = tuple(float(x) for x in ARGV[3].split(",")) if \
            len(ARGV) > 3 else (0.05, 0.1, 0.2, 0.4, 0.8)
        gen = ARGV[2] if len(ARGV) > 2 else "clock_local"
        p3(ARGV[1], gen, om)
    elif mode == "ctrl":
        ctrl(ARGV[1], *(int(a) for a in ARGV[2:3]))
    elif mode == "collect":
        collect()
