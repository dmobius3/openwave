"""M5.21.9 P0-P2: the delta < 0 lattice rung, core-profile
sensitivity, and the fixed-J build, on the certified 4D stack.

P0  mcurve/relax at delta < 0: the bounded family from the negative
    delta verify (findings/m5_21_9_negdelta_note.md). The +delta
    rigid family was lattice-singular under FIRE (M5.21.8 § 5);
    the question is whether the bounded family relaxes.
P1  core-profile sensitivity: blend the field toward its azimuthal
    average inside a core radius r_c (3 profiles); measure whether
    E(m), kin, and relax survival depend on the choice (the one
    author-gated regularization left; fork (a)).
P2  fixed-J build: minimize E_stat(M) + J^2/(4 kin(M)) at FIXED J
    (constraint-carried isorotation; omega* = J/(2 kin) at the
    minimum). Implemented as FIRE with the adaptive coefficient
    -J^2/(4 kin^2) on kin_grad, refreshed every refresh_every
    iterations. Hold read: convergence + the post-relax kick return.

Modes:
    mcurve dl=<v> [g= n= L= m0= m1= nm=]
    relax  dl=<v> m=<v> [maxit= dt0= prof= rc=]
    fixedj dl=<v> m=<v> J=<v> [maxit= refresh=]
Out: data/m5_21_9_lat_*.json
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

_s3 = importlib.util.spec_from_file_location(
    "ins4", os.path.join(HERE, "m5_21_3_a_4d.py"))
INS4 = importlib.util.module_from_spec(_s3)
_s3.loader.exec_module(INS4)
_s8 = importlib.util.spec_from_file_location(
    "lat8", os.path.join(HERE, "m5_21_8_b_lattice.py"))
LAT8 = importlib.util.module_from_spec(_s8)
_s8.loader.exec_module(LAT8)


def save(tag, obj):
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, f"m5_21_9_lat_{tag}.json"),
              "w") as f:
        json.dump(obj, f, indent=1)
    print("saved", tag, flush=True)


def core_profile(M, cfg, rc, kind="tanh"):
    """blend M toward its azimuthal average inside rho < ~rc
    (the vortex-axis regularization; kind picks the transition)."""
    if rc <= 0:
        return M
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS4.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    # azimuthal average via 90-degree rotational symmetrization of
    # the axis neighborhood is overkill; the axis limit of the
    # dressed family is the phi-average, computed cheaply by
    # averaging the 4 grid quadrant rotations (C4 subgroup).
    Mq = M.copy()
    for k in (1, 2, 3):
        Mr = np.rot90(M, k=k, axes=(0, 1))
        # rotate tensor indices by the same C4 element
        c, s = np.cos(k * np.pi / 2), np.sin(k * np.pi / 2)
        Rz = np.array([[1, 0, 0, 0], [0, c, -s, 0],
                       [0, s, c, 0], [0, 0, 0, 1.0]])
        Mq = Mq + np.einsum("ab,...bc,dc->...ad", Rz, Mr, Rz)
    Mq = Mq / 4.0
    x = rho / rc
    if kind == "tanh":
        w = 0.5 * (1 + np.tanh(4.0 * (x - 1.0)))
    elif kind == "quad":
        w = np.clip(x * x, 0.0, 1.0)
    else:  # hard
        w = (x >= 1.0).astype(float)
    W = w[..., None, None]
    return INS4.sym4(W * M + (1 - W) * Mq)


def mcurve(dl, g=8.0, n=32, L=48.0, m0=-0.35, m1=0.35, nm=29):
    cfg = INS4.base_cfg(s=-1.0, g=g, n=n, L=L, delta=dl)
    ms = np.linspace(m0, m1, nm)
    rows = []
    for m in ms:
        M = LAT8.dressed(cfg, m)
        e_u, e_v = INS4.e_parts(M, cfg)
        a0 = LAT8.a0_unit(cfg, m)
        k = INS4.kin_of(M, a0, cfg)
        rows.append({"m": float(m), "E_u": float(e_u),
                     "E_v": float(e_v), "kin": float(k)})
        print(f"  dl={dl:+.3g} m {m:+.3f} E_u {e_u:10.4f} "
              f"E_v {e_v:.2e} kin {k:+.4f}", flush=True)
    Es = [r["E_u"] + r["E_v"] for r in rows]
    i = int(np.argmin(Es))
    mfit = ms[i]
    if 0 < i < len(ms) - 1:
        a, b, c = Es[i - 1], Es[i], Es[i + 1]
        mfit = ms[i] - 0.5 * (ms[1] - ms[0]) * (c - a) / (c - 2 * b + a)
    out = {"dl": dl, "g": g, "n": n, "L": L, "rows": rows,
           "m_star_lattice": float(mfit),
           "m_star_his": float(0.5 * np.log((1 + g) / (g - 1))),
           "E_min": float(min(Es))}
    save(f"mcurve_dl{dl:+g}_g{g:g}_n{n}", out)
    return out


def relax(dl, m, maxit=3000, dt0=0.02, dt_max=0.2, prof="none",
          rc=0.0):
    cfg = INS4.base_cfg(s=-1.0, g=8.0, delta=dl)
    M0 = LAT8.dressed(cfg, m)
    if prof != "none":
        M0 = core_profile(M0, cfg, rc, kind=prof)
    free = ~INS4.pin_shell(cfg["n"], cfg["h"])
    tag = f"relax_dl{dl:+g}_m{m:+.3f}_{prof}{rc:g}"
    M, info = INS4.fire(M0, cfg, free, max_iter=maxit,
                        log_every=250, tag=tag, dt0=dt0,
                        dt_max=dt_max)
    e0, e1 = INS4.e_parts(M0, cfg), INS4.e_parts(M, cfg)
    a0 = LAT8.a0_unit(cfg, m)
    dv = float(np.sqrt(np.sum((M - M0) ** 2))
               / max(np.sqrt(np.sum((M0 - INS4.vac4(cfg)) ** 2)),
                     1e-300))
    out = {"dl": dl, "m": m, "prof": prof, "rc": rc,
           "E0": float(sum(e0)), "E_end": float(sum(e1)),
           "E_u_end": float(e1[0]), "E_v_end": float(e1[1]),
           "rel_move": dv,
           "kin_end": float(INS4.kin_of(M, a0, cfg)),
           "stop": info["stop"], "finite": bool(np.isfinite(
               np.sum(M))), }
    save(tag, out)
    np.savez_compressed(os.path.join(
        DATA, f"m5_21_9_lat_{tag}_end.npz"), M=M.astype(np.float32))
    return out


def fixedj(dl, m, J, maxit=4000, refresh=200, dt0=0.02,
           dt_max=0.2, prof="none", rc=0.0):
    """fixed-J relax: FIRE on E_stat + J^2/(4 kin), coefficient on
    kin_grad refreshed every `refresh` iters (adaptive constraint)."""
    cfg = INS4.base_cfg(s=-1.0, g=8.0, delta=dl)
    M0 = LAT8.dressed(cfg, m)
    if prof != "none":
        M0 = core_profile(M0, cfg, rc, kind=prof)
    a0 = LAT8.a0_unit(cfg, m)
    free = ~INS4.pin_shell(cfg["n"], cfg["h"])
    M = M0.copy()
    hist = []
    stop = "maxit"
    # Fixed-J gradient: d/dM [E_stat + J^2/(4 kin)]
    #   = grad_stat - (J/(2 kin))^2 * kin_grad = grad_stat
    #     - omega*^2 * kin_grad.
    # fire() adds +omega^2 * kin_grad, so run it with kin_grad
    # sign-wrapped (module-level patch; fire resolves the name at
    # call time), refreshing omega* = J/(2 kin) every `refresh`
    # iterations (the adaptive constraint carry).
    _orig_kin_grad = INS4.kin_grad
    INS4.kin_grad = lambda Mx, a0x, cfgx, da0_of=None: \
        -_orig_kin_grad(Mx, a0x, cfgx, da0_of)
    try:
        for outer in range(0, maxit, refresh):
            kin = INS4.kin_of(M, a0, cfg)
            if kin <= 0:
                stop = f"kin_nonpositive({kin:.3g})"
                break
            om = J / (2.0 * kin)
            Mn, info = INS4.fire(M, cfg, free, max_iter=refresh,
                                 a0=a0, omega=float(om),
                                 log_every=refresh, tag="", dt0=dt0,
                                 dt_max=dt_max)
            M = Mn
            if not np.isfinite(np.sum(M)):
                stop = "non_finite"
                break
            e_u, e_v = INS4.e_parts(M, cfg)
            row = {"iter": outer + refresh, "kin": float(kin),
                   "omega_star": float(om),
                   "E_rot": float(J * J / (4.0 * kin)),
                   "E_u": float(e_u), "E_v": float(e_v),
                   "fire_stop": info["stop"]}
            hist.append(row)
            print(row, flush=True)
            if info["stop"] not in ("maxit", "max_iter"):
                stop = info["stop"]
                break
    finally:
        INS4.kin_grad = _orig_kin_grad
    kin_f = INS4.kin_of(M, a0, cfg) if np.isfinite(np.sum(M)) else \
        float("nan")
    out = {"dl": dl, "m": m, "J": J, "prof": prof, "rc": rc,
           "stop": stop, "hist": hist,
           "kin_final": float(kin_f),
           "omega_star_final": float(J / (2 * kin_f))
           if np.isfinite(kin_f) and kin_f > 0 else None}
    tag = f"fixedj_dl{dl:+g}_m{m:+.3f}_J{J:g}"
    save(tag, out)
    if np.isfinite(np.sum(M)):
        np.savez_compressed(os.path.join(
            DATA, f"m5_21_9_lat_{tag}_end.npz"),
            M=M.astype(np.float32))
    return out


if __name__ == "__main__":
    mode = sys.argv[1]
    kw = dict(a.split("=", 1) for a in sys.argv[2:])
    if mode == "mcurve":
        mcurve(float(kw["dl"]), g=float(kw.get("g", 8)),
               n=int(kw.get("n", 32)), L=float(kw.get("L", 48)),
               m0=float(kw.get("m0", -0.35)),
               m1=float(kw.get("m1", 0.35)),
               nm=int(kw.get("nm", 29)))
    elif mode == "relax":
        relax(float(kw["dl"]), float(kw["m"]),
              maxit=int(kw.get("maxit", 3000)),
              dt0=float(kw.get("dt0", 0.02)),
              prof=kw.get("prof", "none"),
              rc=float(kw.get("rc", 0)))
    elif mode == "fixedj":
        fixedj(float(kw["dl"]), float(kw["m"]), float(kw["J"]),
               maxit=int(kw.get("maxit", 4000)),
               refresh=int(kw.get("refresh", 200)),
               prof=kw.get("prof", "none"),
               rc=float(kw.get("rc", 0)))
