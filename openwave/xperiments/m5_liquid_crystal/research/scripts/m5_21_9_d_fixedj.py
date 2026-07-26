"""M5.21.9 P2 (stop-string fix: fire reports max_iter): the fixed-J isorotation electron on the certified
M5.21.3 4D state (the author's s = -1 vacuum branch).

Arena decision (deviation logged 2026-07-20 11:57, checkpoint): the
boost-dressed family is V4-flat (E_v = 0 on the whole family), so
free FIRE descent runs away through the boost channel at BOTH signs
of delta (P0b). The certified stable 4D statics is the M5.21.3
embedded-electron endpoint, where V4 fires off-spectrum: the fixed-J
state is built THERE.

Construction: constraint-carried J on the clock flow
    a0 = clock_local(M)   (gen_catalog; unit Frobenius norm;
                           the conjugation-convention ZBW clock)
    E(J; M) = E_stat(M) + J^2 / (4 kin(M)),  omega* = J / (2 kin)
FIRE descent with kin_grad sign-wrapped (fire adds +omega^2 kin_grad;
fixed-J needs -omega*^2 kin_grad), omega* AND a0 refreshed from the
current M every `refresh` iterations (adaptive carry; a0 frozen
within a window, the kin_grad frozen-a0 convention).

Reads per J rung:
    convergence (fire stop + E trace), omega*_final, kin drift,
    core spectrum (center eigenvalues) vs the start state,
    twist_read on start vs end (the axial-twist linear response b:
    the statics-grade field-coupling channel for the Larmor read)

Mode:  build om=<omega_target> [maxit= refresh=]
       (J = 2 kin0 om: the target sets the initial omega*)
Out:   data/m5_21_9_fixedj_om<val>.json + _end.npz
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


def spec_core(M, cfg):
    n = cfg["n"]
    c = n // 2
    lam = np.linalg.eigvalsh(M[c, c, c])
    return [float(x) for x in lam]


def a0_conj(cfg, M):
    """the CONJUGATION-tangent clock flow (the canonical convention,
    m5_theory_canonical § 6 anti-recipe row: the physical inertia is
    the conjugation value; gen_catalog's GM - MG^T probe flow is
    antisymmetric and exits the symmetric configuration space):
    a0 = w [W, M] with W = rotation about the local leading spatial
    eigenvector, unit Frobenius norm."""
    w = INS4.envelope(cfg)[..., None, None]
    lam, V = np.linalg.eigh(M[..., 1:4, 1:4])
    vhat = V[..., :, 2]
    W = np.zeros(vhat.shape[:-1] + (4, 4))
    n1, n2, n3 = vhat[..., 0], vhat[..., 1], vhat[..., 2]
    W[..., 1, 2], W[..., 1, 3] = -n3, n2
    W[..., 2, 1], W[..., 2, 3] = n3, -n1
    W[..., 3, 1], W[..., 3, 2] = -n2, n1
    A = w * (W @ M - M @ W)
    return A / max(np.sqrt(np.sum(A * A)), 1e-300)


def build(om_target, maxit=1500, refresh=300, tagbase="p1_s-1",
          conv="conj"):
    M, cfg = INS4.load_p1(tagbase)
    free = ~INS4.pin_shell(cfg["n"], cfg["h"])
    a0 = a0_conj(cfg, M) if conv == "conj" else \
        INS4.gen_catalog(cfg, M)["clock_local"]
    kin0 = INS4.kin_of(M, a0, cfg)
    J = 2.0 * kin0 * om_target
    e0 = INS4.e_parts(M, cfg)
    tw0 = INS4.twist_read(M, a0, cfg)
    start = {"kin0": float(kin0), "J": float(J),
             "om_target": float(om_target),
             "E_u0": float(e0[0]), "E_v0": float(e0[1]),
             "spec_core0": spec_core(M, cfg),
             "twist0": {k: tw0[k] for k in ("b", "d", "k_star")}}
    print("start", start, flush=True)
    if kin0 <= 0:
        raise SystemExit("kin0 nonpositive; construction undefined")

    _orig = INS4.kin_grad
    INS4.kin_grad = lambda Mx, a0x, cfgx, da0_of=None: \
        -_orig(Mx, a0x, cfgx, da0_of)
    hist, stop = [], "maxit"
    try:
        for outer in range(0, maxit, refresh):
            a0 = a0_conj(cfg, M) if conv == "conj" else \
                INS4.gen_catalog(cfg, M)["clock_local"]
            kin = INS4.kin_of(M, a0, cfg)
            if kin <= 0:
                stop = f"kin_nonpositive({kin:.4g})"
                break
            om = J / (2.0 * kin)
            Mn, info = INS4.fire(M, cfg, free, max_iter=refresh,
                                 a0=a0, omega=float(om),
                                 log_every=refresh, tag="",
                                 dt0=0.02, dt_max=0.2)
            M = Mn
            if not np.isfinite(np.sum(M)):
                stop = "non_finite"
                break
            eu, ev = INS4.e_parts(M, cfg)
            row = {"iter": outer + refresh, "kin": float(kin),
                   "omega_star": float(om),
                   "E_rot": float(J * J / (4 * kin)),
                   "E_u": float(eu), "E_v": float(ev),
                   "fire_stop": info["stop"]}
            hist.append(row)
            print(row, flush=True)
            if info["stop"] not in ("maxit", "max_iter"):
                stop = info["stop"]
                break
    finally:
        INS4.kin_grad = _orig

    out = {"start": start, "hist": hist, "stop": stop}
    if np.isfinite(np.sum(M)):
        a0f = a0_conj(cfg, M) if conv == "conj" else \
            INS4.gen_catalog(cfg, M)["clock_local"]
        kinf = INS4.kin_of(M, a0f, cfg)
        ef = INS4.e_parts(M, cfg)
        twf = INS4.twist_read(M, a0f, cfg)
        out["final"] = {
            "kin_final": float(kinf),
            "omega_star_final": float(J / (2 * kinf))
            if kinf > 0 else None,
            "E_u": float(ef[0]), "E_v": float(ef[1]),
            "spec_core": spec_core(M, cfg),
            "twist": {k: twf[k] for k in ("b", "d", "k_star")},
            "rel_move_from_start": float(
                np.sqrt(np.sum((M - INS4.load_p1(tagbase)[0]) ** 2))
                / np.sqrt(np.sum(INS4.load_p1(tagbase)[0] ** 2)))}
        print("final", out["final"], flush=True)
        np.savez_compressed(
            os.path.join(DATA, f"m5_21_9_fixedj_{conv}_om"
                         f"{om_target:g}_end.npz"),
            M=M.astype(np.float32))
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA,
                           f"m5_21_9_fixedj_{conv}_om"
                           f"{om_target:g}.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print("saved", f"fixedj_om{om_target:g}", flush=True)
    return out


if __name__ == "__main__":
    kw = dict(a.split("=", 1) for a in sys.argv[2:]) \
        if len(sys.argv) > 2 else {}
    build(float(kw.get("om", sys.argv[1].split("=")[-1]
                       if len(sys.argv) > 1 and "=" in sys.argv[1]
                       else 0.2)),
          maxit=int(kw.get("maxit", 1500)),
          refresh=int(kw.get("refresh", 300)),
          conv=kw.get("conv", "conj"))
