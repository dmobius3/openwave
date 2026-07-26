"""M5.21.9 P3: the Larmor precession read on the fixed-J electron
(4x4 dynamics; the author's acceptance protocol).

The author's protocol (2026-07-19 05:09): "introduce constant
external magnetic field by constant field derivative in its
direction, and temporal field derivative - should lead to electron
precession with frequency proportional to magnetic field strength".

OUR IMPLEMENTATION CHOICES (author-gated conventions, flagged in the
note; the twist_read linear response is exactly 0 at the symmetric
state, so the read must be dynamical):
    background   M_B(x) = eps * (x T_y - y T_x), T_i = the symmetric
                 tilt generator (E_{0i} + E_{i0}): both spatial
                 derivatives are constant tilts, so the commutator
                 F_xy picks up a CONSTANT z-oriented block ~ eps^2
                 (the author's "constant field derivative in its
                 direction"); the measured field strength B_meas =
                 the F_xy EM-block norm, reported instead of eps
    dynamics     M_tt = -grad E / h^3 - gamma(r) M_t, leapfrog
                 kick-drift-kick (the M5.21.6 pattern lifted to 4x4),
                 interior gamma = 0, absorbing outer ramp; boundary
                 shell pinned to vacuum + M_B (the background is
                 boundary-maintained)
    init         M(0) = fixed-J state + M_B, M_t(0) = omega* a0
                 (the isorotation velocity; "temporal field
                 derivative" per the protocol)
    observable   the J-vector track: J_k(t) = <M_t, a_rot_k>_eta
                 projections on the three GLOBAL rotation flows
                 (gen_catalog rot_x/rot_z + rot_y built here);
                 precession = the azimuth of (J_x, J_y) advancing
                 linearly; Omega vs B_meas over >= 3 eps rungs

Gates (cap 3 tries, pre-registered):
    GL4a  eps = 0, gamma = 0: E_tot = E + KE drift < 1% / 200 steps
    GL4b  eps = 0: the J-vector holds direction (no self-precession
          above the numerical floor)
Modes:
    gates            run GL4a/GL4b
    run eps=<v> [steps= dt= snap=]   one precession run
    ladder           eps in {0.5e-3, 1e-3, 2e-3} runs + the fit
Out: data/m5_21_9_larmor_*.json
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

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def tilt(i):
    T = np.zeros((4, 4))
    T[0, i], T[i, 0] = 1.0, 1.0
    return T


def m_background(cfg, eps, rev=False):
    # B flips with the COMMUTATOR orientation, not with eps (F ~ eps^2):
    # rev=True swaps the relative sign, F_xy -> -F_xy (the -B mirror)
    X, Y, Z = INS4.coords(cfg["n"], cfg["h"])
    sgn = 1.0 if rev else -1.0
    return eps * (X[..., None, None] * tilt(2)[None, None, None]
                  + sgn * Y[..., None, None] * tilt(1)[None, None, None])


def b_measured(cfg, eps):
    """the constant F_xy of the background alone: eps^2 [T_y, -T_x]
    eta-commutator block norm."""
    Ty, Tx = tilt(2), tilt(1)
    F = eps * Ty @ ETA @ (-eps * Tx) - (-eps * Tx) @ ETA @ (eps * Ty)
    return float(np.sqrt(np.sum(F * F)))


def rot_flows(cfg, M):
    """unit-Frobenius global-rotation flows a_rot_{x,y,z} of M."""
    Js = {}
    for nm, (a, b) in (("x", (2, 3)), ("y", (3, 1)), ("z", (1, 2))):
        G = np.zeros((4, 4))
        G[a, b], G[b, a] = -1.0, 1.0
        A = G @ M - M @ G.T
        Js[nm] = A / max(np.sqrt(np.sum(A * A)), 1e-300)
    return Js


def j_vector(Mt, flows):
    return {k: float(np.sum(Mt * A)) for k, A in flows.items()}


def gamma_of(cfg, g_max=0.5, r0=None, r1=None):
    X, Y, Z = INS4.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    r0 = r0 or 0.70 * cfg["L"] / 2
    r1 = r1 or 0.95 * cfg["L"] / 2
    t = np.clip((r - r0) / max(r1 - r0, 1e-9), 0.0, 1.0)
    return g_max * t * t


def leap(M, Mt, cfg, free, gam, dt):
    """the M5.21.6 leap_step lifted to 4x4: velocity masked to the
    free cells on EVERY kick (pinned shell carries no KE), implicit
    gamma damping."""
    h3 = cfg["h"] ** 3
    fw = free[..., None, None]
    gw = gam[..., None, None]
    F = -INS4.grad(M, cfg) / h3
    Mt = (Mt + 0.5 * dt * F) * fw
    Mt = Mt / (1.0 + 0.5 * dt * gw)
    M = INS4.sym4(M + dt * Mt)
    F = -INS4.grad(M, cfg) / h3
    Mt = (Mt + 0.5 * dt * F) * fw
    Mt = Mt / (1.0 + 0.5 * dt * gw)
    return M, Mt


def e_tot(M, Mt, cfg):
    eu, ev = INS4.e_parts(M, cfg)
    ke = 0.5 * cfg["h"] ** 3 * float(np.sum(Mt * Mt))
    return float(eu + ev), ke


def load_state(om_target=0.2):
    Z = np.load(os.path.join(DATA,
                             f"m5_21_9_fixedj_om{om_target:g}_end.npz"))
    M = Z["M"].astype(np.float64)
    with open(os.path.join(DATA,
                           f"m5_21_9_fixedj_om{om_target:g}.json")) \
            as f:
        rec = json.load(f)
    _, cfg = INS4.load_p1("p1_s-1")
    return M, cfg, rec


def setup(eps, om_target=0.2, env=False, rev=False, jsign=1.0):
    # jsign = -1 flips the CARRIED J (Mt(0) = -omega* a0): the audit's
    # J-flip discriminator: a genuine Larmor rate is J-sign-even,
    # the drift and relaxation channels flip with J
    M, cfg, rec = load_state(om_target)
    a0 = INS4.gen_catalog(cfg, M)["clock_local"]
    om = jsign * rec["final"]["omega_star_final"]
    MB = m_background(cfg, eps, rev=rev)
    if env:
        # localized field: constant over the core (envelope ~ 1 for
        # r < ~7, renv = 10), vacuum at the boundary; kills the
        # linear-growth self-energy blowup seen at eps >= 0.03
        MB = MB * INS4.envelope(cfg)[..., None, None]
    Mtot = INS4.sym4(M + MB)
    free = ~INS4.pin_shell(cfg["n"], cfg["h"])
    Mt = om * a0 * free[..., None, None]
    gam = gamma_of(cfg)
    return Mtot, Mt, cfg, free, gam, om


def gates(steps=200, dt=0.02):
    M, Mt, cfg, free, gam, om = setup(0.0)
    gam0 = np.zeros_like(gam)
    E0, K0 = e_tot(M, Mt, cfg)
    flows = rot_flows(cfg, M)
    j0 = j_vector(Mt, flows)
    for k in range(steps):
        M, Mt = leap(M, Mt, cfg, free, gam0, dt)
    E1, K1 = e_tot(M, Mt, cfg)
    j1 = j_vector(Mt, rot_flows(cfg, M))
    drift = abs((E1 + K1) - (E0 + K0)) / max(abs(E0 + K0), 1e-300)
    out = {"steps": steps, "dt": dt,
           "E0": E0, "KE0": K0, "E1": E1, "KE1": K1,
           "etot_drift_rel": float(drift),
           "GL4a_pass": bool(drift < 0.01),
           "j0": j0, "j1": j1}
    print(json.dumps(out, indent=1), flush=True)
    with open(os.path.join(DATA, "m5_21_9_larmor_gates.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    return out


def run(eps, steps=2000, dt=0.02, snap=50, om_target=0.2,
        tag=None, env=False, rev=False, jsign=1.0):
    M, Mt, cfg, free, gam, om = setup(eps, om_target, env=env,
                                      rev=rev, jsign=jsign)
    Bm = b_measured(cfg, eps)
    rows = []
    for k in range(steps + 1):
        if k % snap == 0:
            flows = rot_flows(cfg, M)
            j = j_vector(Mt, flows)
            E, KE = e_tot(M, Mt, cfg)
            phi = float(np.arctan2(j["y"], j["x"]))
            rows.append({"t": k * dt, "J": j, "phi_xy": phi,
                         "E": E, "KE": KE})
            print(f"  t {k * dt:7.2f} Jz {j['z']:+.4e} "
                  f"Jxy ({j['x']:+.3e},{j['y']:+.3e}) "
                  f"phi {phi:+.3f} E {E:.4f}", flush=True)
        if k < steps:
            M, Mt = leap(M, Mt, cfg, free, gam, dt)
    out = {"eps": eps, "B_meas": Bm, "om_start": om, "dt": dt,
           "steps": steps, "rows": rows}
    tg = tag or f"eps{eps:g}"
    with open(os.path.join(DATA, f"m5_21_9_larmor_{tg}.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print("saved larmor", tg, flush=True)
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "gates"
    kw = dict(a.split("=", 1) for a in sys.argv[2:])
    if mode == "gates":
        gates(steps=int(kw.get("steps", 200)),
              dt=float(kw.get("dt", 0.02)))
    elif mode == "run":
        run(float(kw["eps"]), steps=int(kw.get("steps", 2000)),
            dt=float(kw.get("dt", 0.02)),
            snap=int(kw.get("snap", 50)),
            env=bool(int(kw.get("env", 0))),
            rev=bool(int(kw.get("rev", 0))),
            jsign=float(kw.get("jsign", 1)),
            tag=kw.get("tag"))
    elif mode == "ladder":
        for eps in (0.0, 2e-3, 5e-3, 1e-2):
            run(eps, steps=int(kw.get("steps", 2000)),
                dt=float(kw.get("dt", 0.02)),
                snap=int(kw.get("snap", 50)))
