"""M5.21.11 branch-ID readers: the frozen § 5 triple per ladder endpoint.

Framework § 5 (FROZEN): branch identity at every rung = the
(charge class, line census, core class) triple matching the branch's
rung-0.30 triple; topology only, ENERGY IS NEVER USED. A read with
active gap flags is an EXCLUDED rung (§ 4), not a guessed identity.

Instruments (all pre-existing, imported unchanged):
    charge class   surface-flux of the oriented long axis
                   (m5_21_4_a_pair.orient_v1 -> mermin_B -> cube_flux
                   at halves 12 AND 18; class = flux rounded to the
                   nearest half-integer, both halves must agree)
    line census    m5_23_2_tracer.trace (count + closure verdicts,
                   min_size = 2 voxels; summary = (n_lines, sorted
                   verdict multiset))
    core class     which eigenvalue pair (12 vs 23) holds the smaller
                   mean gap over the potential-excess core set
                   (v_density > 0.5 max)
    gap flags      per-contour min eigen-gap (m5_21_2b_b_split
                   contour_winding, band 1, planes +-4/+-8,
                   rho_c in {6, 10, 14, 18}); flag bar delta/6 and the
                   endpoint min-gap guard max(delta/60, 1e-3)
                   (pre-run interpretation pins, checkpoint doc)

Modes: all | one tag=t11lad_A_n48_d0.3
Out: ../data/m5_21_11_read_<tag>.json (idempotent per endpoint)
     + merged ../data/m5_21_11_reads.json (all mode)
"""
from __future__ import annotations

import glob
import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PAIR = _load("pair", "m5_21_4_a_pair.py")
TRACER = _load("tracer", "m5_23_2_tracer.py")
INS = _load("ins2b", "m5_21_2b_a_instrument.py")
SPLIT = _load("bsplit", "m5_21_2b_b_split.py")

W2_PIN = 0.002758100
FLAG_PLANES = (-8.0, -4.0, 4.0, 8.0)
FLAG_RHOS = (6.0, 10.0, 14.0, 18.0)


def charge_read(M, n, h):
    nhat, ncf = PAIR.orient_v1(M)
    B = PAIR.mermin_B(nhat, h)
    cfg = {"n": n, "h": h}
    out = {"orient_conflicts": int(ncf)}
    classes = []
    for half in (12.0, 18.0):
        q = PAIR.cube_flux(B, cfg, 0.0, half)
        out[f"q_flux_{int(half)}"] = q
        classes.append(None if not np.isfinite(q)
                       else round(2.0 * q) / 2.0)
    out["class_12"], out["class_18"] = classes
    out["charge_class"] = classes[0] if classes[0] == classes[1] \
        else None
    return out


def census_read(M, h):
    T = TRACER.trace(M, h=h, min_size=2)
    lines = [{k: ln[k] for k in ("n_vox", "verdict", "extent_h")}
             for ln in T["lines"]]
    return {"bulk_split": T["bulk_split"], "warn": T["warn"],
            "n_lines": len(lines),
            "verdicts": sorted(ln["verdict"] for ln in lines),
            "lines": lines}


def core_read(M, n, h, delta):
    cfg = INS.base_cfg(seed="A", term="T2", stencil="sym", eps=0.0,
                       n=n, L=n * h, delta=delta, bc="pinned",
                       w2=W2_PIN)
    vd = INS.v_density(M, cfg)
    core = vd > 0.5 * vd.max()
    lam = np.linalg.eigvalsh(M)
    g12 = float(np.mean(lam[core][:, 1] - lam[core][:, 0]))
    g23 = float(np.mean(lam[core][:, 2] - lam[core][:, 1]))
    return {"gap12_core": g12, "gap23_core": g23,
            "core_class": "12" if g12 <= g23 else "23",
            "n_core_vox": int(core.sum())}


def flag_read(M, h, delta):
    bar = delta / 6.0
    rows, active = [], 0
    for z in FLAG_PLANES:
        n = M.shape[0]
        k = int(round((n - 1) / 2.0 + z / h))
        P = M[:, :, max(0, min(n - 1, k))]
        for rho in FLAG_RHOS:
            w, gmin = SPLIT.contour_winding(P, h, rho, 1)
            fl = bool(gmin < bar)
            active += int(fl)
            rows.append({"z": z, "rho": rho, "half_units": w,
                         "min_gap": round(gmin, 5), "flag": fl})
    return {"bar": bar, "n_flagged": active,
            "n_contours": len(rows), "contours": rows}


def read_one(tag):
    Z = np.load(os.path.join(DATA, f"m5_21_11_end_{tag}.npz"))
    M = Z["M"].astype(np.float64)
    n, h, delta = int(Z["n"]), float(Z["h"]), float(Z["delta"])
    row = json.load(open(os.path.join(DATA,
                                      f"m5_21_11_row_{tag}.json")))
    guard_bar = max(delta / 60.0, 1e-3)
    out = {"tag": tag, "branch": str(Z["branch"]), "n": n,
           "delta": delta,
           "charge": charge_read(M, n, h),
           "census": census_read(M, h),
           "core": core_read(M, n, h, delta),
           "flags": flag_read(M, h, delta),
           "gap_guard": {"bar": guard_bar,
                         "min_gap_end": row["min_gap_end"],
                         "silent": bool(row["min_gap_end"] > guard_bar
                                        and row["stop"]
                                        != "non-finite")}}
    with open(os.path.join(DATA, f"m5_21_11_read_{tag}.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"tag": tag,
                      "charge": out["charge"]["charge_class"],
                      "n_lines": out["census"]["n_lines"],
                      "verdicts": out["census"]["verdicts"],
                      "core": out["core"]["core_class"],
                      "flags": out["flags"]["n_flagged"],
                      "guard_silent": out["gap_guard"]["silent"]}),
          flush=True)
    return out


def all_mode():
    outs = {}
    for p in sorted(glob.glob(
            os.path.join(DATA, "m5_21_11_end_t11lad_*.npz"))):
        tag = os.path.basename(p)[len("m5_21_11_end_"):-len(".npz")]
        rp = os.path.join(DATA, f"m5_21_11_read_{tag}.json")
        if os.path.exists(rp):
            outs[tag] = json.load(open(rp))
            continue
        if not os.path.exists(os.path.join(
                DATA, f"m5_21_11_row_{tag}.json")):
            continue
        outs[tag] = read_one(tag)
    with open(os.path.join(DATA, "m5_21_11_reads.json"), "w") as f:
        json.dump(outs, f, indent=1)
    print(f"collected {len(outs)} reads")
    return outs


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    kw = dict(a.split("=", 1) for a in sys.argv[2:])
    if mode == "one":
        read_one(kw["tag"])
    else:
        all_mode()
