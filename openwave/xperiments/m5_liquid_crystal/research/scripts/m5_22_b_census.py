"""M5.22 census driver: rotate + relax the author seeds, endpoint instruments.

Consumes by import: INS = the certified 3D stack (m5_21_2b_a_instrument:
T2 term w2 = 0.002758100, sym stencil, eps = 0, L = 48, FIRE, pin shell),
PAIR = the signed-charge instruments (m5_21_4_a_pair: orient_v1, mermin_B,
cube_flux), SEEDS = the seed factory (m5_22_a_seeds, GATE-A validated).

Config of record: bc = pinned (far-field topology held at the seed's),
delta = 0.3, census rung n = 32 (maxit 8000), confirmation n = 48
(maxit 12000). Box L = 48 fixed.

Endpoint instruments (per relaxed state):
  E parts / virial / r_half           (INS)
  total charge Q_far + PROFILE Q(half) over growing centered cubes of
    the oriented Mermin-Ho flux -> the core/shell read. The director
    lift's GLOBAL sign is a gauge (orient_v1 continuity fix): |Q| is the
    charge class; WITHIN a state the profile's relative signs are
    meaningful (positive core / negative shell)
  core ledger: connected components (scipy.ndimage.label, 3^3
    connectivity) of the low-eigengap set at two thresholds, pin shell
    excluded; per component: voxels, mean/max rho, z extent, min gap
  symmetry read: relative L2 deviation of M from its 90-degree
    azimuthal rotation about z (cylindrical-symmetry survival under
    full-3D relaxation; the seed's own deviation printed as baseline)
  seed 2D winding q2d carried in the row (pre-rotation diagnostic)

Modes:
  relax fam=N s=-0.5 conv=plane scale=6 n=32 maxit=8000 bc=pinned
        [delta=0.3 rc=2.0 tag=...]
  instruments tag=<saved tag>       (re-read a saved endpoint npz)
  collect                           (merge rows -> m5_22_census.json)
Out: ../data/m5_22_row_<tag>.json + m5_22_end_<tag>.npz per run.
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
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")


def _load(name, fname):
    sp = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


INS = _load("ins", "m5_21_2b_a_instrument.py")
PAIR = _load("pair", "m5_21_4_a_pair.py")
SEEDS = _load("seeds", "m5_22_a_seeds.py")

W2_T2 = 0.002758100          # the certified T2 normalization (2b note)


# ================= endpoint instruments =================
def charge_profile(M, cfg):
    nhat, ncf = PAIR.orient_v1(M)
    B = PAIR.mermin_B(nhat, cfg["h"])
    far0 = 0.5 * cfg["L"] - 4.0 * cfg["h"]
    halves = list(np.arange(3.0, far0 + 1e-9, 1.5 * cfg["h"]))
    if halves[-1] < far0 - 1e-9:
        halves.append(far0)
    prof = [(float(hv), PAIR.cube_flux(B, cfg, 0.0, hv))
            for hv in halves]
    qs = [q for _, q in prof if np.isfinite(q)]
    q_far = float(np.mean(qs[-2:])) if len(qs) >= 2 else float("nan")
    return {"n_conflicts": int(ncf), "q_far": q_far,
            "profile": [[hv, q] for hv, q in prof]}


def core_ledger(M, cfg, thresholds=(0.05, 0.10)):
    n, h = cfg["n"], cfg["h"]
    lam = np.linalg.eigvalsh(M)
    gap = np.minimum(lam[..., 1] - lam[..., 0],
                     lam[..., 2] - lam[..., 1])
    X, Y, Z = INS.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    interior = ~INS.pin_shell(n, h)
    out = {"min_gap": float(gap[interior].min())}
    for th in thresholds:
        mask = (gap < th) & interior
        lab, ncomp = ndimage.label(mask,
                                   structure=np.ones((3, 3, 3), int))
        comps = []
        for k in range(1, ncomp + 1):
            sel = lab == k
            comps.append({
                "voxels": int(sel.sum()),
                "rho_mean": float(rho[sel].mean()),
                "rho_max": float(rho[sel].max()),
                "z_min": float(Z[sel].min()),
                "z_max": float(Z[sel].max()),
                "gap_min": float(gap[sel].min())})
        comps.sort(key=lambda c: -c["voxels"])
        out[f"th{th:g}"] = {"n_components": ncomp, "components": comps[:8]}
    return out


def sym_read(M):
    """relative deviation from the 90-degree rotation about z:
    M'(x) = R M(R^{-1} x) R^T, R = z-rotation by +90 deg. On the grid
    (X, Y) -> rot90 of axes (0, 1)."""
    R = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    Mg = np.rot90(M, k=1, axes=(0, 1))
    Mr = np.einsum("ab,...bc,dc->...ad", R, Mg, R)
    num = float(np.sqrt(np.sum((M - Mr) ** 2)))
    den = float(np.sqrt(np.sum(M ** 2)))
    return num / max(den, 1e-300)


def instruments(M, M0, cfg, fam, s, extra=None):
    e_u, e_d, e_v = INS.e_parts(M, cfg)
    row = {
        "E_end": float(e_u + e_d + e_v), "E_u": float(e_u),
        "E_d": float(e_d), "E_v": float(e_v),
        "E_seed": float(sum(INS.e_parts(M0, cfg))),
        "virial_resid": float((-e_u + e_d + 3 * e_v)
                              / max(e_u + e_d + e_v, 1e-300)),
        "r_half": INS.r_half(M, cfg),
        "min_gap_end": INS.min_gap(M),
        "consistency": INS.consistency(M, cfg),
        "charge": charge_profile(M, cfg),
        "charge_seed": charge_profile(M0, cfg),
        "ledger": core_ledger(M, cfg),
        "sym_dev_end": sym_read(M),
        "sym_dev_seed": sym_read(M0),
        "q2d_seed": SEEDS.q2d(fam, s),
    }
    if extra:
        row.update(extra)
    return row


# ================= the relax driver =================
def relax(kw):
    fam = kw.pop("fam")
    s = float(kw.pop("s"))
    conv = kw.pop("conv", "plane")
    scale = float(kw.pop("scale", 6.0))
    r_c = float(kw.pop("rc", 2.0))
    qshift = float(kw.pop("qshift", 0.0))
    tag = kw.pop("tag", "")
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       **{k: v for k, v in kw.items()})
    if not tag:
        tag = (f"{fam}{s:+g}_{conv}_sc{scale:g}_n{cfg['n']}"
               f"_{cfg['bc']}_d{cfg['delta']:g}")
        if qshift:
            tag += f"_q{qshift:g}"
    M0 = SEEDS.seed_3d(fam, s, cfg, conv=conv, scale=scale, r_c=r_c,
                       qshift=qshift)
    free = ~INS.pin_shell(cfg["n"], cfg["h"]) if cfg["bc"] == "pinned" \
        else np.ones((cfg["n"],) * 3, dtype=bool)
    t0 = time.time()
    M, states, info = INS.fire(M0, cfg, free, max_iter=cfg["maxit"],
                               log_every=500, tag=tag)
    row = {"fam": fam, "s": s, "conv": conv, "scale": scale, "r_c": r_c,
           "qshift": qshift, "tag": tag}
    row.update({k: cfg[k] for k in ("n", "L", "h", "delta", "bc",
                                    "maxit", "w2")})
    row.update(instruments(M, M0, cfg, fam, s,
                           extra={"stop": info["stop"],
                                  "trace": info["trace"][-4:],
                                  "wall_s": time.time() - t0}))
    os.makedirs(DATA, exist_ok=True)
    np.savez_compressed(os.path.join(DATA, f"m5_22_end_{tag}.npz"),
                        M=M.astype(np.float32),
                        M0=M0.astype(np.float32),
                        delta=cfg["delta"], h=cfg["h"], n=cfg["n"],
                        fam=fam, s=s, conv=conv, scale=scale)
    with open(os.path.join(DATA, f"m5_22_row_{tag}.json"), "w") as f:
        json.dump(row, f, indent=1)
    print(json.dumps({k: row[k] for k in
                      ("tag", "E_end", "E_u", "virial_resid", "r_half",
                       "stop", "sym_dev_end", "min_gap_end")}
                     | {"q_far": row["charge"]["q_far"],
                        "q_far_seed": row["charge_seed"]["q_far"],
                        "q2d": row["q2d_seed"],
                        "xratio":
                        row["consistency"]["xstencil_ratio"]}))
    return row


def extend(kw):
    """continue a saved endpoint: mode extend tag=... maxit=...
    [noise=0.0 newtag=...]; noise > 0 = the perturbation-return probe
    (uniform noise of that amplitude on free cells before relaxing)."""
    tag = kw.pop("tag")
    noise = float(kw.pop("noise", 0.0))
    newtag = kw.pop("newtag", tag + ("_prt" if noise else "_ext"))
    z = np.load(os.path.join(DATA, f"m5_22_end_{tag}.npz"))
    with open(os.path.join(DATA, f"m5_22_row_{tag}.json")) as f:
        old = json.load(f)
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(z["n"]), delta=float(z["delta"]),
                       bc=old["bc"], maxit=int(kw.get("maxit", 8000)))
    M0 = z["M"].astype(np.float64)
    free = ~INS.pin_shell(cfg["n"], cfg["h"]) if cfg["bc"] == "pinned" \
        else np.ones((cfg["n"],) * 3, dtype=bool)
    if noise:
        rng = np.random.default_rng(7)
        pert = rng.normal(size=M0.shape) * noise
        pert = 0.5 * (pert + pert.swapaxes(-1, -2))
        M0 = M0 + pert * free[..., None, None]
    t0 = time.time()
    M, states, info = INS.fire(M0, cfg, free, max_iter=cfg["maxit"],
                               log_every=500, tag=newtag)
    fam, s = old["fam"], old["s"]
    row = {k: old.get(k, 0.0 if k == "qshift" else None)
           for k in ("fam", "s", "conv", "scale", "r_c", "qshift",
                     "n", "L", "h", "delta", "bc", "w2")}
    row.update({"tag": newtag, "maxit": cfg["maxit"], "parent": tag,
                "noise": noise})
    row.update(instruments(M, M0, cfg, fam, s,
                           extra={"stop": info["stop"],
                                  "trace": info["trace"][-4:],
                                  "wall_s": time.time() - t0}))
    np.savez_compressed(os.path.join(DATA, f"m5_22_end_{newtag}.npz"),
                        M=M.astype(np.float32),
                        M0=M0.astype(np.float32),
                        delta=cfg["delta"], h=cfg["h"], n=cfg["n"],
                        fam=fam, s=s, conv=old["conv"],
                        scale=old["scale"])
    with open(os.path.join(DATA, f"m5_22_row_{newtag}.json"),
              "w") as f:
        json.dump(row, f, indent=1)
    print(json.dumps({"tag": newtag, "E_end": row["E_end"],
                      "q_far": row["charge"]["q_far"],
                      "stop": row["stop"],
                      "fmax": row["trace"][-1]["fmax"]}))
    return row


def rerun_instruments(tag):
    z = np.load(os.path.join(DATA, f"m5_22_end_{tag}.npz"))
    with open(os.path.join(DATA, f"m5_22_row_{tag}.json")) as f:
        row = json.load(f)
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(z["n"]), delta=float(z["delta"]),
                       bc=row["bc"])
    row.update(instruments(z["M"].astype(np.float64),
                           z["M0"].astype(np.float64), cfg,
                           str(z["fam"]), float(z["s"])))
    with open(os.path.join(DATA, f"m5_22_row_{tag}.json"), "w") as f:
        json.dump(row, f, indent=1)
    print(f"instruments refreshed for {tag}")


def collect():
    import glob
    rows = {}
    for p in sorted(glob.glob(os.path.join(DATA, "m5_22_row_*.json"))):
        key = os.path.basename(p)[len("m5_22_row_"):-len(".json")]
        with open(p) as f:
            rows[key] = json.load(f)
    out = os.path.join(DATA, "m5_22_census.json")
    with open(out, "w") as f:
        json.dump(rows, f, indent=1)
    print(f"collected {len(rows)} rows -> {os.path.basename(out)}")


def parse_kv(argv):
    kw = {}
    casts = {"n": int, "L": float, "delta": float, "maxit": int,
             "scale": float, "rc": float, "s": float, "w2": float,
             "qshift": float}
    for a in argv:
        k, v = a.split("=", 1)
        kw[k] = casts[k](v) if k in casts else v
    return kw


if __name__ == "__main__":
    ARGV = sys.argv[1:]
    mode = ARGV[0] if ARGV else "collect"
    if mode == "relax":
        relax(parse_kv(ARGV[1:]))
    elif mode == "extend":
        extend(parse_kv(ARGV[1:]))
    elif mode == "instruments":
        rerun_instruments(parse_kv(ARGV[1:])["tag"])
    elif mode == "collect":
        collect()
