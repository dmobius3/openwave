"""M5.32 R20-1 / R20-2: the paper's three-axis hedgehog (the three-lepton
mechanism, report sections 363 to 364) solved on the certified static sector
4 I1 + V4 and on 4 I1 + V_spec (the author's degree-8 spectral potential),
with the swapped-assignment, degenerate-vacuum and g 32 controls
(ledger section 6.9).

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1), N = M eta. Code
branch s = -1: M_vac = diag(8, 1, delta, 0), the N-spectrum
q = (-8, 1, delta, 0), delta = 0.3, g = 8 (g = 32 control). The static
energy of a field (the R3 / R19 instrument, m5_32_r19_entrants.energy_grad):
    E[M] = 4 h^3 sum_br wt sum_cells I1(A) + V[M],  A_i = d_i M (sym stencil)
    V = V4     = W1 sum_p (tr N^p - C_p)^2, C_p = sum_i q_i^p     (R20-1)
    V = V4^dd  = the same with q = (-8, 1, delta, delta)         (the degenerate control)
    V = V_spec = gamma h^3 sum_cells tr[P(N)^2], P(x) = prod (x - q_i),
                 gamma from R20-0 (c)                             (R20-2)
Seeds (the record's builder generalized, m5_21_4_a_pair._tensor_from_nhat):
    M3(x) = lam_r n n^T + lam_phi phi phi^T + lam_theta theta theta^T
    n radial (the hedgehog), phi azimuthal, theta = phi x n (polar; phi and
    theta are undefined on the polar axis, where the record sets phi = y),
    isotropic blend at the core w = 1 - exp(-(r / r_c)^2), r_c = 4, toward
    a I with a = (lam_r + lam_phi + lam_theta) / 3; embedded with M_00 = 8.
    S_1 (1, delta, 0)  the record's electron   [identical to seed_pair 'single', gated]
    S_d (delta, 1, 0)  the delta-axis object
    S_0 (0, 1, delta)  the 0-axis object
    S_k' the transverse swap of each: (1, 0, delta), (delta, 0, 1), (0, delta, 1)
    on the degenerate vacuum (1, delta, delta): S_1^dd (1, delta, delta),
    S_d^dd (delta, 1, delta), S_0^dd (delta, delta, 1) (the transverse swap
    of S_d^dd: the two must coincide there, the mechanism's own null)
Relaxation: R3's FIRE with energy-monotone backtracking (dt0 0.02, dt_max
0.2, alpha 0.1, dt halved on rejection, dt_min 1e-7), the pinned Dirichlet
shell B3.pin_shell depth 1.6 at the seed values, STEPS_ACC = 4500 accepted
steps (three times R3's), IT_CAP 9000; a resumable checkpoint (M, v, the
FIRE state, the trace) every CKPT_EVERY accepted steps, picked up by a
rerun of the same stage. Kill rules as R19 (RUNAWAY on the
time row, DIVERGED, LOCUS-HIT).
The convergence gate (pre-registered): CONVERGED iff the relative energy
drift over the last third of the accepted steps is below 1e-3 AND max |G|
at the end is two decades below the seed's; else FALLING, and the `extend`
stage doubles the budget once from the saved end field; the label is
carried into every quoted number.
Reads per object and box: E (total, curvature, potential); the degree of
EACH eigenvector of the 3x3 block (rank 0 / 1 / 2 by eigenvalue) on the
lattice cubes of half-width 6, 9, 12 (and the far cube), by the record's
Mermin flux after orientation by continuity (the escape read; the number
of orientation conflicts is the line-defect read); on the polar axis
(the four columns nearest x = y = 0, 4 < |z| < L/2 - 2): the eigenvalue
gaps and the biaxiality, and the energy in a tube of radius 3 along z
minus the same tube along x and y, per unit length (the string read);
the core ball r < 4: gaps, biaxiality, energy. At collect: E(L) on the
ladder n32 L48 / n48 L72 / n64 L96 (h 1.5), its slope, and the outcomes
    AXIS_ESCAPES    the winding eigenvector's |degree| on the r 9 and r 12
                    cubes fell below 0.5 at the end (seed: 1)
    AXIS_STRING     E grows with L: the L72 -> L96 increment is positive
                    and above half the L48 -> L72 increment, and above
                    0.01 per unit length
    AXIS_CONVERGED  the L72 -> L96 increment is below 1e-2 |E| or below a
                    quarter of the L48 -> L72 increment
    AXIS_BOX_LIMITED E falls with L on both rungs and the virial E_curv / V
                    on the largest box exceeds 3 (the scaling force says
                    expand: the box, not the field, sets the size; added
                    2026-09-13 after the author's section 579 remark,
                    before any ladder increment was read)
    THREE_AXES_DISTINCT / AXES_DEGENERATE on the three converged energies,
                    "equal" meaning within three times the ladder's
                    resolution (the largest L72 -> L96 increment)
Koide Q = sum E / (sum sqrt E)^2 and the ratios are reads, not gates.
The Derrick reads (per row, from the end field): the virial ratio
E_curv / V (3 at a quartic + potential equilibrium, E = a/R + b R^3), the
direct dilation scan E(lambda), lambda 0.9 to 1.1 by cubic interpolation
(linear as the check), its slope and parabolic minimum, the half-energy
radius and R_* = r (E_curv / 3 V)^(1/4).

STAGES (python3 m5_32_r20_1_axes.py STAGE [--workers W] [--mem-budget GB]):
    smoke    the seed identity gate, the energy-wrapper gates, one tiny job
    all      the 22 jobs (9 main, 7 controls, 6 V_spec), largest first
    n64      the three n64 main rows only (the n48 rung was dropped on 2026-09-13)
    extend   +STEPS_ACC from the end field of every FALLING row
    ddladder the degenerate-vacuum trio at n48 L72 (added 2026-09-13)
    phase2   extend (n32 and n48 only; n64 carried at its label) + ddladder in one pool
    collect  tables, outcomes, plots (merges every side-pool JSON)
Out: ../data/m5_32_r20_1_axes.json (partials after every job),
     ../data/m5_32_r20_1/*.npz (local), ../plots/m5_32_r20_1_*.png
"""
from __future__ import annotations

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import argparse  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa: E402
import multiprocessing as mp  # noqa: E402

import numpy as np  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
OUT_JSON = os.path.join(DATA, "m5_32_r20_1_axes.json")
MAIN_JSON = OUT_JSON                     # the main batch's file; a side pool redirects OUT_JSON, never MAIN_JSON
OUT_NPZ = os.path.join(DATA, "m5_32_r20_1")
R20_0_JSON = os.path.join(DATA, "m5_32_r20_0_class.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    argv = sys.argv
    sys.argv = [argv[0]]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


R0 = _load("m5_32_r20_0_class", "m5_32_r20_0_class.py")
EN = R0.EN
B3 = R0.B3
LAG = R0.LAG
PAIR = _load("m5_21_4_a_pair", "m5_21_4_a_pair.py")
R3 = _load("m5_32_r3_ii_pair", "m5_32_r3_ii_pair.py")       # the record's boost dressing (frame_reads)
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
W1 = B3.W1
DELTA = 0.3
G_MAIN, G_CTRL = 8.0, 32.0
LADDER = ((32, 48.0), (48, 72.0), (64, 96.0))
STEPS_ACC = 4500
IT_CAP = 9000
RUNAWAY_FACTOR = 3.0
DIVE_FLOOR = -1e6
R_CORE = 4.0
R_TUBE = 3.0
T0 = time.time()

OBJECTS = {"S1": (1.0, DELTA, 0.0), "Sd": (DELTA, 1.0, 0.0), "S0": (0.0, 1.0, DELTA),
           "S1p": (1.0, 0.0, DELTA), "Sdp": (DELTA, 0.0, 1.0), "S0p": (0.0, DELTA, 1.0),
           "S1dd": (1.0, DELTA, DELTA), "Sddd": (DELTA, 1.0, DELTA), "S0dd": (DELTA, DELTA, 1.0)}
WINDING_RANK = {"S1": 2, "Sd": 1, "S0": 0, "S1p": 2, "Sdp": 1, "S0p": 0, "S1dd": 2, "Sddd": None, "S0dd": None}


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def cfg_of(n, L, g):
    return B3.base_cfg(s=-1.0, g=g, n=n, L=float(L), delta=DELTA)


def params_of(g):
    return LAG.default_params(s=-1.0, g=g)


def gamma_r20_0():
    with open(R20_0_JSON) as f:
        return float(json.load(f)["c"]["gamma"])


# ================= seeds =================
def seed_axes(cfg, lam, r_c=R_CORE):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    nhat = PAIR._nhat_from_alpha(n, h, np.arctan2(rho, Z))
    rhos = np.where(rho < 1e-12, 1e-12, rho)
    phihat = np.stack([-Y / rhos, X / rhos, np.zeros_like(Z)], axis=-1)
    near = rho < 1e-9
    if np.any(near):
        phihat[near] = np.array([0.0, 1.0, 0.0])
    dot = np.einsum("...a,...a->...", phihat, nhat)[..., None]
    ph = phihat - dot * nhat
    ph = ph / np.maximum(np.linalg.norm(ph, axis=-1)[..., None], 1e-300)
    th = np.cross(ph, nhat)
    lr, lp, lt = lam
    S = (lr * nhat[..., :, None] * nhat[..., None, :] + lp * ph[..., :, None] * ph[..., None, :]
         + lt * th[..., :, None] * th[..., None, :])
    a = (lr + lp + lt) / 3.0
    r = np.sqrt(X * X + Y * Y + Z * Z)
    w = 1.0 - np.exp(-((r / r_c) ** 2))
    M3 = w[..., None, None] * S + (1.0 - w[..., None, None]) * (a * np.eye(3))
    return B3.embed34(M3, cfg)


# ================= the energy with the potential hook =================
def pot_of(kind, cfg, gamma=None):
    """None: the stack's V4 as EN computes it. ('v4dd'): the degenerate targets. ('vspec'): the author's potential."""
    if kind == "v4std":
        return None
    if kind == "v4dd":
        return ("v4", R0.roots_of(cfg, degenerate=True), W1)
    if kind == "vspec":
        return ("vspec", R0.roots_of(cfg), gamma)
    raise ValueError(kind)


def energy_grad(M, cfg, p, pot, need_grad=True):
    E, G, info = EN.energy_grad(M, cfg, "I1", c=1.0, p=p, need_grad=need_grad)
    if pot is None or not np.isfinite(E):
        return E, G, info
    q_std = R0.roots_of(cfg)
    ev, gv = R0.v4_energy_grad(M, cfg, q_std, W1, need_grad)
    if pot[0] == "v4":
        ep, gp = R0.v4_energy_grad(M, cfg, pot[1], pot[2], need_grad)
    else:
        ep, gp = R0.vspec_energy_grad(M, cfg, pot[1], pot[2], need_grad)
    E = E - ev + ep
    if need_grad:
        G = G - gv + gp
    return E, G, info


def energy_parts(M, cfg, p, pot):
    r = EN.block_reads(M, cfg, "I1")
    out = {"E_curv": r["E_curv_I1"], "V4_std": r["V4"], "min_gap": r["min_gap"]}
    if pot is None:
        out["V"] = r["V4"]
    elif pot[0] == "v4":
        out["V"] = R0.v4_energy_grad(M, cfg, pot[1], pot[2], need_grad=False)[0]
    else:
        out["V"] = R0.vspec_energy_grad(M, cfg, pot[1], pot[2], need_grad=False)[0]
    out["E_total"] = out["E_curv"] + out["V"]
    return out


def density(M, cfg, pot):
    """per-cell energy density (curvature + potential), h^3-weighted, the certified stencil."""
    h3 = cfg["h"] ** 3
    e = np.zeros(M.shape[:3])
    for br, (A, wt) in B3.a_fields(M, cfg).items():
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                e += wt * 4.0 * B3.inner_eta(F, F)
    e = h3 * e
    if pot is None:
        q, w = R0.roots_of(cfg), W1
        t = LAG.v4_traces_np(M)
        cp = [sum(qi ** k for qi in q) for k in range(1, 5)]
        e = e + h3 * w * sum((t[k] - cp[k]) ** 2 for k in range(4))
    elif pot[0] == "v4":
        t = LAG.v4_traces_np(M)
        cp = [sum(qi ** k for qi in pot[1]) for k in range(1, 5)]
        e = e + h3 * pot[2] * sum((t[k] - cp[k]) ** 2 for k in range(4))
    else:
        e = e + h3 * pot[2] * R0.vspec_density(M, pot[1])
    return e


# ================= the descent (R3 / R19 FIRE with the potential hook) =================
CKPT_EVERY = 500      # accepted steps between resumable checkpoints (the 2026-09-13 kills lost 3.7 h of n64 descent)


def _ckpt_save(path, M, v, state):
    np.savez_compressed(path + ".npz", M=M, v=v)
    with open(path + ".json", "w") as f:
        json.dump(state, f)


def _ckpt_load(path):
    if not (os.path.exists(path + ".npz") and os.path.exists(path + ".json")):
        return None
    Z = np.load(path + ".npz")
    with open(path + ".json") as f:
        st = json.load(f)
    return Z["M"], Z["v"], st


def descend(M0, cfg, p, pot, steps_acc, it_cap, tag, log_every=100, dt0=0.02, dt_max=0.2, ckpt_path=None):
    free = (~B3.pin_shell(cfg["n"], cfg["h"]))[..., None, None].astype(float)
    M = M0.copy()
    E0, G, info = energy_grad(M, cfg, p, pot)
    m0i_seed = float(np.max(np.abs(M0[..., 0, 1:])))
    out = {"E0": float(E0), "steps_acc_budget": steps_acc, "it_cap": it_cap,
           "pin": "B3.pin_shell depth 1.6 (Dirichlet at the seed values)",
           "fire": {"dt0": dt0, "dt_max": dt_max, "alpha0": 0.1, "dt_min": 1e-7},
           "max_abs_M0i_seed": m0i_seed, "trace": [], "resumed_from_ckpt": None}
    if not np.isfinite(E0) or G is None:
        out.update(stop="DIVERGED (seed energy undefined)", verdict="DIVERGED", steps_run=0, accepted=0, E_end=float("nan"))
        return M, out
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    dt_min = 1e-7
    F = -G * free
    E_prev = E0
    stop = "budget"
    n_rej, n_rej_locus, n_acc = 0, 0, 0
    fmax0 = float(np.max(np.abs(F)))
    out["fmax_seed"] = fmax0
    fmax = fmax0
    it = 0
    runaway = None
    E = E0
    ck = _ckpt_load(ckpt_path) if ckpt_path else None
    if ck is not None:
        M, v, st = ck
        E, G, info = energy_grad(M, cfg, p, pot)
        F = -G * free
        it, n_acc, n_rej, n_rej_locus = st["it"], st["acc"], st["rej"], st["rej_locus"]
        dt, alpha, n_up, E_prev, fmax = st["dt"], st["alpha"], st["n_up"], st["E_prev"], float(np.max(np.abs(F)))
        out["trace"] = st["trace"]
        out["resumed_from_ckpt"] = {"it": it, "acc": n_acc, "E": float(E)}
        log(f"{tag} RESUMED from checkpoint at it {it} acc {n_acc} E {E:.6f}")
    while it < it_cap and n_acc < steps_acc:
        it += 1
        P = float(np.sum(F * v))
        if P > 0.0:
            n_up += 1
            vn = np.sqrt(np.sum(v * v))
            fn = np.sqrt(np.sum(F * F))
            v = (1 - alpha) * v + alpha * (F / max(fn, 1e-300)) * vn
            if n_up > 5:
                dt = min(dt * 1.1, dt_max)
                alpha *= 0.99
        else:
            v[:] = 0.0
            alpha, n_up = 0.1, 0
        v_try = v + dt * F
        M_try = M + dt * v_try
        E, G, info = energy_grad(M_try, cfg, p, pot)
        locus_loss = not info["ok"]
        reject = locus_loss or not np.isfinite(E) or E > E_prev + 1e-12 * max(abs(E_prev), 1.0)
        if reject:
            n_rej += 1
            if locus_loss:
                n_rej_locus += 1
            dt *= 0.5
            v[:] = 0.0
            alpha, n_up = 0.1, 0
            if dt < dt_min:
                stop = "LOCUS-HIT" if locus_loss else "STALLED (dt collapsed, no descent direction accepted)"
                break
            continue
        n_acc += 1
        M, v, E_prev = M_try, v_try, E
        F = -G * free
        fmax = float(np.max(np.abs(F)))
        m0i = info["max_abs_M0i"]
        if m0i_seed > 0 and m0i > RUNAWAY_FACTOR * m0i_seed:
            runaway = {"step": it, "accepted": n_acc, "max_abs_M0i": m0i, "E": float(E)}
            stop = "RUNAWAY"
            break
        if E < DIVE_FLOOR:
            stop = "DIVERGED (dive floor)"
            break
        if it % log_every == 0 or n_acc == steps_acc:
            row = {"it": it, "acc": n_acc, "E": float(E), "fmax": fmax, "dt": dt, "min_gap": info["min_gap"]}
            out["trace"].append(row)
            log(f"{tag} it {it:5d} acc {n_acc:5d} E {E:14.6f} fmax {fmax:.3e} dt {dt:.2e} rej {n_rej}")
        if ckpt_path and n_acc % CKPT_EVERY == 0 and n_acc < steps_acc:
            _ckpt_save(ckpt_path, M, v, {"it": it, "acc": n_acc, "rej": n_rej, "rej_locus": n_rej_locus, "dt": dt, "alpha": alpha,
                                         "n_up": n_up, "E_prev": float(E_prev), "trace": out["trace"]})
    if stop == "budget" and n_acc < steps_acc:
        stop = f"IT_CAP ({it_cap} iterations before {steps_acc} accepted)"
    out.update({"stop": stop, "steps_run": it, "accepted": n_acc, "rejected": n_rej, "rejected_locus": n_rej_locus,
                "dt_final": dt, "fmax_end": fmax, "E_end": float(E_prev), "E_drop": float(E0 - E_prev), "runaway": runaway})
    if ckpt_path:
        for ext in (".npz", ".json"):
            if os.path.exists(ckpt_path + ext):
                os.remove(ckpt_path + ext)
    tr = out["trace"]
    q = [r for r in tr if r["acc"] >= (2.0 / 3.0) * n_acc] if n_acc else []
    if stop.startswith("budget") or stop.startswith("IT_CAP"):
        if len(q) >= 2:
            dE = q[-1]["E"] - q[0]["E"]
            out["last_third_dE"] = float(dE)
            out["last_third_rel"] = float(abs(dE) / max(abs(q[-1]["E"]), 1.0))
            out["fmax_decades"] = float(np.log10(fmax0 / max(fmax, 1e-300)))
            conv = out["last_third_rel"] < 1e-3 and out["fmax_decades"] >= 2.0
            out["verdict"] = "CONVERGED" if conv else ("FALLING" if dE < 0 else "RISING")
            out["gate"] = {"last_third_rel_lt_1e-3": bool(out["last_third_rel"] < 1e-3), "fmax_fell_2_decades": bool(out["fmax_decades"] >= 2.0)}
        else:
            out["verdict"] = "budget (trace too short)"
    else:
        out["verdict"] = stop
    return M, out


# ================= reads =================
def orient_field(v):
    """orientation by continuity of a unit-vector field (m5_21_4_a_pair.orient_v1 on a given field)."""
    v = v.copy()
    sgn = np.ones(v.shape[:3])
    for ax in range(3):
        dots = np.einsum("...a,...a->...", np.roll(v * sgn[..., None], 1, axis=ax), v * sgn[..., None])
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


def halves_of(cfg):
    far = 0.5 * cfg["L"] - 4.0 * cfg["h"]
    hs = [6.0, 9.0, 12.0] + [x for x in (18.0, 24.0, 30.0, 36.0, 42.0) if x <= far]
    return [x for x in hs if x <= far]


def biaxiality(M3):
    """beta^2 = 1 - 6 (tr S^3)^2 / (tr S^2)^3 on the traceless part (0 uniaxial, 1 maximally biaxial)."""
    S = M3 - (np.einsum("...kk->...", M3) / 3.0)[..., None, None] * np.eye(3)
    t2 = np.einsum("...ij,...ji->...", S, S)
    t3 = np.einsum("...ij,...jk,...ki->...", S, S, S)
    return 1.0 - 6.0 * t3 ** 2 / np.maximum(t2 ** 3, 1e-300)


def reads(M, cfg, p, pot):
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rho = np.sqrt(X * X + Y * Y)
    M3 = M[..., 1:, 1:]
    lam, vec = np.linalg.eigh(M3)                     # ascending: rank 0 (0), rank 1 (delta), rank 2 (1)
    out = {"energy": energy_parts(M, cfg, p, pot)}
    # the escape read: the degree of each eigenvector on the cubes
    deg = {}
    for k in range(3):
        vo, ncf = orient_field(vec[..., :, k])
        Bm = PAIR.mermin_B(vo, h)
        deg[f"rank{k}"] = {"conflicts": ncf, "flux": {f"{hv:g}": PAIR.cube_flux(Bm, cfg, 0.0, hv) for hv in halves_of(cfg)}}
    out["degree"] = deg
    # eigenvalues on the shells (the identity of each rank there)
    sh = {}
    for hv in (6.0, 9.0, 12.0):
        m = (r > hv - 1.0) & (r < hv + 1.0)
        sh[f"{hv:g}"] = {"lam_mean": [float(np.mean(lam[m][:, k])) for k in range(3)],
                         "gap_min": [float(np.min(lam[m][:, k + 1] - lam[m][:, k])) for k in range(2)],
                         "biax_mean": float(np.mean(biaxiality(M3[m])))}
    out["shells"] = sh
    # the string read on the polar axis
    ax = (rho < 0.75 * h) & (np.abs(Z) > R_CORE) & (np.abs(Z) < 0.5 * L - 2.0)
    out["axis"] = {"cells": int(ax.sum()), "gap_min": [float(np.min(lam[ax][:, k + 1] - lam[ax][:, k])) for k in range(2)],
                   "gap_mean": [float(np.mean(lam[ax][:, k + 1] - lam[ax][:, k])) for k in range(2)],
                   "biax_mean": float(np.mean(biaxiality(M3[ax]))), "biax_max": float(np.max(biaxiality(M3[ax])))}
    e = density(M, cfg, pot)
    tube = {}
    for lab, perp, along in (("z", rho, Z), ("x", np.sqrt(Y * Y + Z * Z), X), ("y", np.sqrt(X * X + Z * Z), Y)):
        m = (perp < R_TUBE) & (np.abs(along) > R_CORE) & (np.abs(along) < 0.5 * L - 2.0)
        tube[lab] = {"E": float(np.sum(e[m])), "length": float(2.0 * (0.5 * L - 2.0 - R_CORE))}
    out["tube"] = tube
    out["string_tension_read"] = float((tube["z"]["E"] - 0.5 * (tube["x"]["E"] + tube["y"]["E"])) / tube["z"]["length"])
    # the core ball
    cb = r < R_CORE
    out["core"] = {"E": float(np.sum(e[cb])), "gap_min": [float(np.min(lam[cb][:, k + 1] - lam[cb][:, k])) for k in range(2)],
                   "biax_max": float(np.max(biaxiality(M3[cb]))), "lam_center": [float(x) for x in lam[r == r.min()][0]]}
    out["E_density_total"] = float(np.sum(e))
    out["E_far_shell_r_gt_12"] = float(np.sum(e[r > 12.0]))
    out["derrick"] = derrick_reads(M, cfg, p, pot)
    out["frame"] = frame_reads(M, cfg, p, pot)
    out["radial"] = radial_profile(M, cfg, pot)
    return out


# ================= the Derrick reads (added 2026-09-13 17:xx UTC after the author's section 579 remark, before any ladder increment was read) =================
def rescale_field(M, cfg, lam, order=3):
    """M_lam(x) = M(x / lam): the object dilated by lam about the box center (values beyond the grid: the nearest edge value)."""
    n = cfg["n"]
    c = (n - 1) / 2.0
    src = (np.arange(n) - c) / lam + c
    X, Y, Z = np.meshgrid(src, src, src, indexing="ij")
    coords = np.stack([X.ravel(), Y.ravel(), Z.ravel()])
    out = np.empty_like(M)
    for a in range(4):
        for b in range(a, 4):
            v = map_coordinates(M[..., a, b], coords, order=order, mode="nearest").reshape(n, n, n)
            out[..., a, b] = v
            out[..., b, a] = v
    return out


def derrick_reads(M, cfg, p, pot):
    """the scaling structure of a static quartic + potential object: E(lam) = a / lam + b lam^3 for a pure dilation,
    so at an equilibrium E_curv = 3 V (the virial) and R_* = R (E_curv / 3 V)^(1/4); measured here directly by
    dilating the end field (cubic interpolation, the linear one as the check) and by the virial ratio."""
    E1 = energy_parts(M, cfg, p, pot)
    lams = (0.9, 0.95, 1.05, 1.1)
    out = {"E_1": E1, "virial_E_curv_over_V": float(E1["E_curv"] / E1["V"]) if E1["V"] > 0 else float("nan")}
    for order in (3, 1):
        Es = {}
        for lam in lams:
            Es[f"{lam:g}"] = energy_parts(rescale_field(M, cfg, lam, order), cfg, p, pot)
        t = lambda k: Es[k]["E_total"]                                   # noqa: E731
        dE = (-t("1.1") + 8.0 * t("1.05") - 8.0 * t("0.95") + t("0.9")) / (12.0 * 0.05)
        d2E = (-t("1.1") + 16.0 * t("1.05") - 30.0 * E1["E_total"] + 16.0 * t("0.95") - t("0.9")) / (12.0 * 0.05 ** 2)
        out[f"order{order}"] = {"E_of_lambda": Es, "dE_dlambda_at_1": float(dE), "d2E_dlambda2_at_1": float(d2E),
                                "lambda_star_parabolic": float(1.0 - dE / d2E) if d2E > 0 else float("nan"),
                                "scaling_force": "EXPAND" if dE < 0 else "CONTRACT"}
    # the core radius: the half-energy radius of the density and the radius where the largest spectral gap reaches half its far value
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z).ravel()
    e = density(M, cfg, pot).ravel()
    o = np.argsort(r)
    cum = np.cumsum(e[o])
    r_half = float(r[o][np.searchsorted(cum, 0.5 * cum[-1])])
    out["r_half_energy"] = r_half
    v = out["virial_E_curv_over_V"]
    out["R_star_from_virial"] = float(r_half * (v / 3.0) ** 0.25) if np.isfinite(v) else float("nan")
    return out


# ================= the frame and radial reads (added 2026-09-13 19:xx UTC after the author's #186 comment 18424648) =================
FRAME_S = (0.02, 0.05)
PROFILE_R = (3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 30.0, 36.0, 42.0)


def frame_reads(M, cfg, p, pot):
    """the undressed minimum's curvature along the record's boost dressing M -> Q M Q^T, Q = exp(s b*(r) n . K)
    (m5_32_r3_ii_pair.boost_at about the center): an eta-orthogonal conjugation, so both potentials are invariant
    (the author's item 4) and only the curvature energy moves; d2E/ds2 < 0 means the undressed minimum is a saddle
    in the full sector along this boost-rotation direction (the author's item 1, R19's RUNAWAY rows)."""
    E0 = energy_parts(M, cfg, p, pot)
    out = {"E_0": E0, "s": list(FRAME_S)}
    for sv in FRAME_S:
        Ep = energy_parts(R3.conj(R3.boost_at(cfg, 0.0, sv)[0], M), cfg, p, pot)
        Em = energy_parts(R3.conj(R3.boost_at(cfg, 0.0, -sv)[0], M), cfg, p, pot)
        out[f"s{sv:g}"] = {"E_plus": Ep["E_total"], "E_minus": Em["E_minus"] if "E_minus" in Em else Em["E_total"],
                           "dV_plus": Ep["V"] - E0["V"], "dV_minus": Em["V"] - E0["V"],
                           "d2E_ds2": float((Ep["E_total"] + Em["E_total"] - 2.0 * E0["E_total"]) / sv ** 2),
                           "dE_ds": float((Ep["E_total"] - Em["E_total"]) / (2.0 * sv)),
                           "max_abs_M0i_plus": float(np.max(np.abs(R3.conj(R3.boost_at(cfg, 0.0, sv)[0], M)[..., 0, 1:])))}
    out["saddle_along_boost_dressing"] = bool(out[f"s{FRAME_S[0]:g}"]["d2E_ds2"] < 0.0)
    return out


def radial_profile(M, cfg, pot):
    """the cumulative energy E(< R) on the shells PROFILE_R (the author's matching-radius request): the smallest R at which
    E(< R) is box-independent across the ladder is read at collect."""
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    e = density(M, cfg, pot)
    pin = B3.pin_shell(n, h)
    out = {"E_total": float(np.sum(e)), "E_pin_shell": float(np.sum(e[pin])),
           "E_lt_R": {f"{R:g}": float(np.sum(e[r < R])) for R in PROFILE_R if R <= 0.5 * L}}
    return out


# ================= jobs =================
def job_tag(obj, potk, n, g):
    return f"{obj}_{potk}_n{n}_g{g:g}"


def job_list(stage, gamma):
    jobs = []
    if stage in ("all", "main"):
        for n, L in reversed(LADDER):
            for obj in ("S1", "Sd", "S0"):
                jobs.append(dict(obj=obj, potk="v4std", n=n, L=L, g=G_MAIN, steps=STEPS_ACC))
    if stage in ("all", "vspec"):
        for n, L in reversed(LADDER[:2]):
            for obj in ("S1", "Sd", "S0"):
                jobs.append(dict(obj=obj, potk="vspec", n=n, L=L, g=G_MAIN, steps=STEPS_ACC, gamma=gamma))
    if stage in ("all", "controls"):
        n, L = LADDER[0]
        for obj in ("S1p", "Sdp", "S0p"):
            jobs.append(dict(obj=obj, potk="v4std", n=n, L=L, g=G_MAIN, steps=STEPS_ACC))
        for obj in ("S1dd", "Sddd", "S0dd"):
            jobs.append(dict(obj=obj, potk="v4dd", n=n, L=L, g=G_MAIN, steps=STEPS_ACC))
        jobs.append(dict(obj="S1", potk="v4std", n=n, L=L, g=G_CTRL, steps=STEPS_ACC))
    # largest boxes first (the long pole)
    jobs.sort(key=lambda j: -j["n"])
    return jobs


def run_job(j):
    t0 = time.time()
    cfg = cfg_of(j["n"], j["L"], j["g"])
    p = params_of(j["g"])
    pot = pot_of(j["potk"], cfg, j.get("gamma"))
    resume = j.get("resume")
    tag = job_tag(j["obj"], j["potk"], j["n"], j["g"]) + (f"_x{j['steps']}" if resume else "")
    row = dict(j)
    row.update({"tag": tag, "h": cfg["h"], "lam": OBJECTS[j["obj"]], "roots": list(R0.roots_of(cfg, degenerate=(j["potk"] == "v4dd")))})
    try:
        if resume:
            M0 = np.load(os.path.join(OUT_NPZ, resume + ".npz"))["M"]
        else:
            M0 = seed_axes(cfg, OBJECTS[j["obj"]])
        row["seed_reads"] = reads(M0, cfg, p, pot)
        os.makedirs(OUT_NPZ, exist_ok=True)
        M, des = descend(M0, cfg, p, pot, j["steps"], IT_CAP if not resume else 2 * j["steps"], tag,
                         ckpt_path=os.path.join(OUT_NPZ, f"{tag}_ckpt"))
        row["descent"] = des
        finite = bool(np.all(np.isfinite(M))) and np.isfinite(des["E_end"])
        row["status"] = "OK" if finite and (des["stop"].startswith("budget") or des["stop"].startswith("IT_CAP")) else des["stop"].split(" ")[0]
        if finite:
            row["end_reads"] = reads(M, cfg, p, pot)
            row["E"] = row["end_reads"]["energy"]["E_total"]
            os.makedirs(OUT_NPZ, exist_ok=True)
            np.savez_compressed(os.path.join(OUT_NPZ, f"{tag}.npz"), M=M.astype(np.float64))
        else:
            row["E"] = None
    except Exception as e:                                # noqa: BLE001
        row["status"] = "DIVERGED"
        row["stop"] = f"exception: {e!r}"
        row["E"] = None
    row["wall_s"] = round(time.time() - t0, 1)
    log(f"DONE {tag} status {row['status']} verdict {row.get('descent', {}).get('verdict')} E {row.get('E')} wall {row['wall_s']}")
    return row


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R20-1 / R20-2", "rows": {}}


def load_json_all():
    """the main JSON plus every side-pool JSON (m5_32_r20_1_axes_<suffix>.json), rows merged (a side pool runs its own
    process pool and writes its own file, so two pools never race on one read-modify-write)."""
    if os.path.exists(MAIN_JSON):
        with open(MAIN_JSON) as f:
            J = json.load(f)
    else:
        J = {"task": "M5.32 R20-1 / R20-2", "rows": {}}
    base = os.path.basename(MAIN_JSON)[:-5]
    for f in sorted(os.listdir(DATA)):
        if f.startswith(base + "_") and f.endswith(".json") and f != base + "_smoke.json":
            with open(os.path.join(DATA, f)) as fh:
                J["rows"].update(json.load(fh).get("rows", {}))
    return J


def save_json(J):
    os.makedirs(DATA, exist_ok=True)
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1)
    os.replace(tmp, OUT_JSON)


MEM_GB = {32: 0.7, 48: 2.0, 64: 3.5}      # measured peak RSS of one energy_grad call per box (2026-09-13)


def run_pool(jobs, workers, label, mem_budget=12.0):
    """at most `workers` jobs at once AND at most `mem_budget` GB of summed per-box peaks (the 2026-09-13 kill:
    12 workers with 3 n64 + 6 n48 + 3 n32 = 23.5 GB on a machine with 12 GB free)."""
    rows = load_json_all()["rows"]
    todo = [j for j in jobs if job_tag(j["obj"], j["potk"], j["n"], j["g"]) + (f"_x{j['steps']}" if j.get("resume") else "") not in rows]
    log(f"{label}: {len(todo)} jobs, {workers} workers, memory budget {mem_budget} GB (done already: {len(jobs) - len(todo)})")
    t0 = time.time()
    ctx = mp.get_context("spawn")
    pending = list(todo)
    running = {}
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as ex:
        while pending or running:
            used = sum(MEM_GB[j["n"]] for j in running.values())
            for j in list(pending):
                if len(running) >= workers:
                    break
                if used + MEM_GB[j["n"]] <= mem_budget + 1e-9:
                    fut = ex.submit(run_job, j)
                    running[fut] = j
                    pending.remove(j)
                    used += MEM_GB[j["n"]]
                    log(f"submit {job_tag(j['obj'], j['potk'], j['n'], j['g'])} (running {len(running)}, memory {used:.1f} GB, pending {len(pending)})")
            if not running:
                break
            done = next(as_completed(list(running)))
            row = done.result()
            running.pop(done)
            J = load_json()
            J["rows"][row["tag"]] = row
            J[f"{label}_wall_s"] = round(time.time() - t0, 1)
            save_json(J)
    log(f"{label} done in {time.time() - t0:.0f} s")


def extend_jobs(max_n=None):
    """max_n: the largest box extended (the 2026-09-13 20:xx UTC decision: the n64 rows are carried at their FALLING
    label, the n32 and n48 rows extended; the label travels with every quoted number)."""
    J = load_json_all()
    jobs = []
    for tag, r in J["rows"].items():
        if max_n is not None and r.get("n", 0) > max_n:
            continue
        if r.get("status") == "OK" and r.get("descent", {}).get("verdict") == "FALLING" and not r.get("resume") and not any(k.startswith(tag + "_x") for k in J["rows"]):
            j = {k: r[k] for k in ("obj", "potk", "n", "L", "g", "steps")}
            if "gamma" in r:
                j["gamma"] = r["gamma"]
            j["resume"] = tag
            j["steps"] = r["steps"]
            jobs.append(j)
    return jobs


def ddladder_jobs():
    """the degenerate-vacuum trio on the second rung (added 2026-09-13 after the report's sections 555 to 556 put the
    neutrino three-axis mechanism on the (1, delta, delta) vacuum)."""
    n, L = LADDER[1]
    return [dict(obj=obj, potk="v4dd", n=n, L=L, g=G_MAIN, steps=STEPS_ACC) for obj in ("S1dd", "Sddd", "S0dd")]


def stage_extend(workers, mem_budget=12.0):
    jobs = extend_jobs()
    jobs.sort(key=lambda j: -j["n"])
    run_pool(jobs, workers, "extend", mem_budget)


def stage_phase2(workers, mem_budget=12.0):
    jobs = extend_jobs(max_n=48) + ddladder_jobs()
    jobs.sort(key=lambda j: -j["n"])
    run_pool(jobs, workers, "phase2", mem_budget)


# ================= smoke =================
def smoke():
    cfg = cfg_of(16, 24.0, G_MAIN)
    p = params_of(G_MAIN)
    out = {}
    # the seed identity: S_1 is the record's single
    Ma = seed_axes(cfg, OBJECTS["S1"])
    Mb = B3.embed34(PAIR.seed_pair(cfg, "single", 0.0), cfg)
    out["seed_identity_S1_max_abs"] = float(np.max(np.abs(Ma - Mb)))
    # the transverse swap has the same eigenvalue multiset, different field
    Mc = seed_axes(cfg, OBJECTS["S1p"])
    lam_a = np.sort(np.linalg.eigvalsh(Ma[..., 1:, 1:]), axis=-1)
    lam_c = np.sort(np.linalg.eigvalsh(Mc[..., 1:, 1:]), axis=-1)
    out["swap_same_spectrum_max_abs"] = float(np.max(np.abs(lam_a - lam_c)))
    out["swap_field_max_abs_diff"] = float(np.max(np.abs(Ma - Mc)))
    # the energy wrapper: the v4-with-standard-targets hook reproduces the stack bitwise
    E1, G1, _ = energy_grad(Ma, cfg, p, None)
    E2, G2, _ = energy_grad(Ma, cfg, p, ("v4", R0.roots_of(cfg), W1))
    out["wrapper_identity"] = {"E_rel": float(abs(E1 - E2) / abs(E1)), "G_max_abs": float(np.max(np.abs(G1 - G2))), "G_scale": float(np.max(np.abs(G1)))}
    # the density sums to the energy
    for potk in ("v4std", "vspec", "v4dd"):
        pot = pot_of(potk, cfg, gamma_r20_0())
        e = density(Ma, cfg, pot)
        ep = energy_parts(Ma, cfg, p, pot)
        Ew, _, _ = energy_grad(Ma, cfg, p, pot, need_grad=False)
        out[f"density_vs_energy_{potk}"] = {"sum_density": float(np.sum(e)), "E_total_parts": ep["E_total"], "E_wrapper": float(Ew),
                                            "rel": float(abs(np.sum(e) - Ew) / abs(Ew))}
    # the wrapper's gradient against complex step for vspec and v4dd on a random field
    rng = np.random.default_rng(1)
    Mr = Ma + 0.05 * B3.sym4(rng.standard_normal(Ma.shape))
    Mr[..., 0, 1:] = 0.0; Mr[..., 1:, 0] = 0.0
    for potk in ("vspec", "v4dd"):
        pot = pot_of(potk, cfg, gamma_r20_0())
        E, G, _ = energy_grad(Mr, cfg, p, pot)
        worst = 0.0
        for _ in range(3):
            D = B3.sym4(rng.standard_normal(Mr.shape))
            D[..., 0, :] = 0.0; D[..., :, 0] = 0.0
            dd = float(np.sum(G * D))
            e = 1e-4
            f = [energy_grad(Mr + k * e * D, cfg, p, pot, need_grad=False)[0] for k in (-2, -1, 1, 2)]
            fd = (f[0] - 8 * f[1] + 8 * f[2] - f[3]) / (12 * e)
            worst = max(worst, abs(fd - dd) / max(abs(dd), 1e-300))
        out[f"wrapper_grad_stencil4_rel_{potk}"] = worst
    # one tiny job end to end under each potential
    for potk in ("v4std", "vspec"):
        j = dict(obj="Sd", potk=potk, n=16, L=24.0, g=G_MAIN, steps=20, gamma=gamma_r20_0())
        cfg16 = cfg_of(16, 24.0, G_MAIN)
        pot = pot_of(potk, cfg16, j["gamma"])
        M0 = seed_axes(cfg16, OBJECTS["Sd"])
        M, des = descend(M0, cfg16, p, pot, 20, 60, f"smoke_{potk}", log_every=5)
        rd = reads(M, cfg16, p, pot)
        out[f"tiny_{potk}"] = {"E0": des["E0"], "E_end": des["E_end"], "accepted": des["accepted"], "stop": des["stop"],
                               "degree_rank1_flux6": rd["degree"]["rank1"]["flux"]["6"], "string_read": rd["string_tension_read"]}
    for k, v in out.items():
        log(f"smoke {k}: {v}")
    with open(os.path.join(DATA, "m5_32_r20_1_smoke.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


# ================= collect =================
def outcome_axis(obj, rows_by_n, potk):
    """rows_by_n: {n: row (the latest extension)}; returns the per-axis outcome and its reads."""
    Ls = {32: 48.0, 48: 72.0, 64: 96.0}
    Es = {n: r["E"] for n, r in rows_by_n.items() if r.get("E") is not None}
    labels = {n: r["descent"]["verdict"] for n, r in rows_by_n.items() if "descent" in r}
    out = {"E": Es, "labels": labels, "outcomes": []}
    wr = WINDING_RANK.get(obj)
    if wr is not None:
        esc = {}
        for n, r in rows_by_n.items():
            if "end_reads" in r:
                fl = r["end_reads"]["degree"][f"rank{wr}"]["flux"]
                fl0 = r["seed_reads"]["degree"][f"rank{wr}"]["flux"]
                esc[n] = {"seed": [fl0.get("9"), fl0.get("12")], "end": [fl.get("9"), fl.get("12")]}
        out["winding_degree"] = esc
        if any(all(abs(v) < 0.5 for v in e["end"] if v is not None and np.isfinite(v)) and any(abs(v) > 0.5 for v in e["seed"] if v is not None) for e in esc.values()):
            out["outcomes"].append("AXIS_ESCAPES")
    if all(n in Es for n in (32, 48, 64)):
        d1, d2 = Es[48] - Es[32], Es[64] - Es[48]
        s2 = d2 / (Ls[64] - Ls[48])
        out["increments"] = {"L48_to_L72": d1, "L72_to_L96": d2, "slope_per_unit_L_hi": s2}
        vir64 = rows_by_n[64].get("end_reads", {}).get("derrick", {}).get("virial_E_curv_over_V", float("nan"))
        out["virial_largest_box"] = vir64
        if d2 > 0 and d2 > 0.5 * d1 and s2 > 0.01:
            out["outcomes"].append("AXIS_STRING")
        elif abs(d2) < 1e-2 * abs(Es[64]) or abs(d2) < 0.25 * abs(d1):
            out["outcomes"].append("AXIS_CONVERGED")
        elif d1 < 0 and d2 < 0 and np.isfinite(vir64) and vir64 > 3.0:
            out["outcomes"].append("AXIS_BOX_LIMITED")
        else:
            out["outcomes"].append("E(L)_UNDECIDED")
    elif all(n in Es for n in (32, 64)):
        # the two-rung rule (the n48 rung was dropped on 2026-09-13 20:xx UTC, before any n64 row was read): one increment over
        # L 48 -> 96; the same thresholds, the resolution being that increment itself
        d = Es[64] - Es[32]
        sl = d / (Ls[64] - Ls[32])
        vir64 = rows_by_n[64].get("end_reads", {}).get("derrick", {}).get("virial_E_curv_over_V", float("nan"))
        out["virial_largest_box"] = vir64
        out["increments"] = {"L48_to_L96": d, "slope_per_unit_L": sl, "two_rung_rule": True}
        if d > 0 and sl > 0.01:
            out["outcomes"].append("AXIS_STRING")
        elif abs(d) < 1e-2 * abs(Es[64]):
            out["outcomes"].append("AXIS_CONVERGED")
        elif d < 0 and np.isfinite(vir64) and vir64 > 3.0:
            out["outcomes"].append("AXIS_BOX_LIMITED")
        else:
            out["outcomes"].append("E(L)_UNDECIDED")
    elif all(n in Es for n in (32, 48)):
        out["increments"] = {"L48_to_L72": Es[48] - Es[32]}
    # the label travels: an E(L) outcome read off rows that are not CONVERGED is provisional (the increment then carries the
    # convergence state as much as the physics; the tube read and the Derrick scan are the discriminators)
    not_conv = [n for n, lab in labels.items() if lab != "CONVERGED"]
    if not_conv:
        out["outcomes"] = [o + (f" (at {','.join(labels[n] for n in sorted(not_conv))} on n{','.join(str(n) for n in sorted(not_conv))})" if o.startswith("AXIS_") and o != "AXIS_ESCAPES" else o) for o in out["outcomes"]]
    out["string_tension_read"] = {n: r["end_reads"]["string_tension_read"] for n, r in rows_by_n.items() if "end_reads" in r}
    prof = {n: r["end_reads"]["radial"]["E_lt_R"] for n, r in rows_by_n.items() if "end_reads" in r and "radial" in r["end_reads"]}
    if len(prof) >= 2:
        common = sorted(set.intersection(*[set(v) for v in prof.values()]), key=float)
        table = {R: {n: prof[n][R] for n in prof} for R in common}
        spread = {R: float(max(v.values()) - min(v.values())) for R, v in table.items()}
        rel = {R: spread[R] / max(abs(np.mean(list(table[R].values()))), 1e-300) for R in common}
        match = [R for R in common if rel[R] <= 0.02]
        out["matching_radius"] = {"E_lt_R_by_box": table, "spread": spread, "rel_spread": rel,
                                  "R_match_2pct": (float(match[0]) if match else None),
                                  "note": "the smallest R at which E(< R) agrees across the boxes to 2 percent; None = the interior energy is itself box-dependent at every R"}
    out["frame"] = {n: {"d2E_ds2_s0.02": r["end_reads"]["frame"]["s0.02"]["d2E_ds2"], "saddle": r["end_reads"]["frame"]["saddle_along_boost_dressing"]}
                    for n, r in rows_by_n.items() if "end_reads" in r and "frame" in r["end_reads"]}
    return out


def backfill_derrick(J):
    """rows relaxed before the Derrick / frame / radial reads existed: compute them from the saved end field (no rerun)."""
    changed = False
    for tag, r in J["rows"].items():
        if r.get("status") != "OK" or "end_reads" not in r:
            continue
        missing = [k for k in ("derrick", "frame", "radial") if k not in r["end_reads"]]
        if not missing:
            continue
        f = os.path.join(OUT_NPZ, tag + ".npz")
        if not os.path.exists(f):
            continue
        cfg = cfg_of(r["n"], r["L"], r["g"])
        p = params_of(r["g"])
        pot = pot_of(r["potk"], cfg, r.get("gamma"))
        M = np.load(f)["M"]
        if "derrick" in missing:
            r["end_reads"]["derrick"] = derrick_reads(M, cfg, p, pot)
            r["end_reads"]["derrick"]["backfilled_from_end_field"] = True
        if "frame" in missing:
            r["end_reads"]["frame"] = frame_reads(M, cfg, p, pot)
            r["end_reads"]["frame"]["backfilled_from_end_field"] = True
        if "radial" in missing:
            r["end_reads"]["radial"] = radial_profile(M, cfg, pot)
            r["end_reads"]["radial"]["backfilled_from_end_field"] = True
        d = r["end_reads"]["derrick"]
        log(f"backfill {tag} {missing}: virial {d['virial_E_curv_over_V']:.2f}, dE/dlam {d['order3']['dE_dlambda_at_1']:+.3f}, boost d2E/ds2 {r['end_reads']['frame']['s0.02']['d2E_ds2']:+.3e}")
        changed = True
    return changed


def collect():
    J = load_json_all()
    backfill_derrick(J)
    rows = J["rows"]

    def latest(obj, potk, n, g):
        base = job_tag(obj, potk, n, g)
        cands = [t for t in rows if t == base or t.startswith(base + "_x")]
        if not cands:
            return None
        return rows[max(cands, key=lambda t: rows[t].get("steps", 0) * (2 if rows[t].get("resume") else 1))]
    res = {"axes": {}}
    for potk in ("v4std", "vspec"):
        per = {}
        for obj in ("S1", "Sd", "S0"):
            rb = {n: r for n, _ in LADDER for r in [latest(obj, potk, n, G_MAIN)] if r is not None}
            per[obj] = outcome_axis(obj, rb, potk)
        # the triple
        conv = {o: per[o]["E"][max(per[o]["E"])] for o in per if per[o]["E"]}
        trip = {"E_largest_box": conv, "largest_box_n": {o: max(per[o]["E"]) for o in per if per[o]["E"]}}
        if len(conv) == 3 and all(v is not None for v in conv.values()):
            incs = [abs(per[o].get("increments", {}).get("L72_to_L96", per[o].get("increments", {}).get("L48_to_L96", np.nan))) for o in per]
            resol = float(np.nanmax(incs)) if np.any(np.isfinite(incs)) else float("nan")
            Es = np.array([conv[o] for o in ("S1", "Sd", "S0")])
            diffs = [abs(Es[i] - Es[j]) for i in range(3) for j in range(i + 1, 3)]
            trip["resolution"] = resol
            trip["pairwise_abs_diffs"] = diffs
            trip["outcome"] = "THREE_AXES_DISTINCT" if np.isfinite(resol) and min(diffs) > 3.0 * resol else ("AXES_DEGENERATE" if np.isfinite(resol) else "LADDER_INCOMPLETE")
            if np.all(Es > 0):
                srt = np.sort(Es)
                trip["koide_Q"] = float(np.sum(Es) / np.sum(np.sqrt(Es)) ** 2)
                trip["ratios_sorted"] = [float(x / srt[0]) for x in srt]
        res["axes"][potk] = {"per_axis": per, "triple": trip}
    # controls
    ctrl = {}
    for obj, potk, g in (("S1p", "v4std", G_MAIN), ("Sdp", "v4std", G_MAIN), ("S0p", "v4std", G_MAIN),
                         ("S1dd", "v4dd", G_MAIN), ("Sddd", "v4dd", G_MAIN), ("S0dd", "v4dd", G_MAIN), ("S1", "v4std", G_CTRL)):
        r = latest(obj, potk, 32, g)
        if r is not None:
            ctrl[job_tag(obj, potk, 32, g)] = {"E": r.get("E"), "verdict": r.get("descent", {}).get("verdict"), "status": r.get("status"),
                                               "string_read": r.get("end_reads", {}).get("string_tension_read"),
                                               "degrees_end": {k: v["flux"] for k, v in r.get("end_reads", {}).get("degree", {}).items()}}
    if "Sddd_v4dd_n32_g8" in ctrl and "S0dd_v4dd_n32_g8" in ctrl and ctrl["Sddd_v4dd_n32_g8"]["E"] and ctrl["S0dd_v4dd_n32_g8"]["E"]:
        a, b = ctrl["Sddd_v4dd_n32_g8"]["E"], ctrl["S0dd_v4dd_n32_g8"]["E"]
        ctrl["degenerate_null_rel"] = float(abs(a - b) / max(abs(a), 1e-300))
    if "S1_v4std_n32_g32" in ctrl:
        ctrl["S1_g32_vs_R3_record_18.970"] = ctrl["S1_v4std_n32_g32"]["E"]
    res["controls"] = ctrl
    # the swap pairs at n32 (main vs swapped)
    sw = {}
    for a, b in (("S1", "S1p"), ("Sd", "Sdp"), ("S0", "S0p")):
        ra, rb = latest(a, "v4std", 32, G_MAIN), latest(b, "v4std", 32, G_MAIN)
        if ra and rb and ra.get("E") is not None and rb.get("E") is not None:
            sw[f"{a}_vs_{b}"] = {"E": [ra["E"], rb["E"]], "diff": rb["E"] - ra["E"]}
    res["swap_pairs_n32"] = sw
    res["collected_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    J["results"] = res
    save_json(J)                       # the merged rows + results land in the main JSON
    plots(J)
    print(json.dumps(res, indent=1))
    return res


def plots(J):
    rows = J["rows"]
    os.makedirs(PLOTS, exist_ok=True)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, potk in zip(axs, ("v4std", "vspec")):
        for obj, mk in (("S1", "o"), ("Sd", "s"), ("S0", "^")):
            Ls, Es = [], []
            for n, L in LADDER:
                cands = [t for t in rows if t.startswith(job_tag(obj, potk, n, G_MAIN))]
                if cands:
                    t = max(cands, key=len)
                    if rows[t].get("E") is not None:
                        Ls.append(L); Es.append(rows[t]["E"])
            if Ls:
                ax.plot(Ls, Es, mk + "-", label=f"{obj} {OBJECTS[obj]}")
        ax.set_xlabel("box L (h = 1.5)"); ax.set_ylabel("E (static, end of descent)")
        ax.set_title(f"the three axes on 4 I1 + {'V4' if potk == 'v4std' else 'V_spec'}" + (" (n32 extended to 9000 steps, n64 at 4500, all FALLING)" if potk == "v4std" else " (n32, CONVERGED)"), fontsize=9)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r20_1_E_of_L.png"), dpi=130)
    # the descents
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for t, r in rows.items():
        tr = r.get("descent", {}).get("trace", [])
        if tr and r.get("n") in (32, 48, 64) and r.get("potk") == "v4std" and r.get("g") == G_MAIN and r["obj"] in ("S1", "Sd", "S0"):
            off = STEPS_ACC if r.get("resume") else 0            # an extension continues the base row's count
            ax.plot([q["acc"] + off for q in tr], [q["E"] for q in tr], lw=0.9, label=t)
    ax.set_xlabel("accepted steps (an extension continues its base row)"); ax.set_ylabel("E"); ax.set_title("descents, 4 I1 + V4 (main)")
    ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r20_1_descents.png"), dpi=130)
    # the radial profiles E(< R) per object and box (the matching-radius read)
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.0), sharey=False)
    for ax, obj in zip(axs, ("S1", "Sd", "S0")):
        for n, L in LADDER:
            cands = [t for t in rows if t.startswith(job_tag(obj, "v4std", n, G_MAIN))]
            for t in sorted(cands):
                rr = rows[t]
                if "end_reads" in rr and "radial" in rr["end_reads"]:
                    pr = rr["end_reads"]["radial"]["E_lt_R"]
                    Rs = sorted(pr, key=float)
                    ax.plot([float(R) for R in Rs], [pr[R] for R in Rs], "o-", ms=3, label=f"{t} (E {rr['E']:.2f})")
        ax.set_xlabel("R"); ax.set_ylabel("E(< R)"); ax.set_title(f"{obj} {OBJECTS[obj]}: the energy inside R, per box", fontsize=9)
        ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r20_1_radial.png"), dpi=130)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=("smoke", "all", "main", "n64", "vspec", "controls", "extend", "ddladder", "phase2", "collect"))
    ap.add_argument("--json-suffix", default="", help="a side pool writes m5_32_r20_1_axes_<suffix>.json (merged at collect)")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--mem-budget", type=float, default=12.0, help="GB, the sum of per-box peak RSS allowed at once")
    a = ap.parse_args()
    if a.json_suffix:
        global OUT_JSON
        OUT_JSON = os.path.join(DATA, f"m5_32_r20_1_axes_{a.json_suffix}.json")
    if a.stage == "smoke":
        smoke()
    elif a.stage == "extend":
        stage_extend(a.workers, a.mem_budget)
    elif a.stage == "phase2":
        stage_phase2(a.workers, a.mem_budget)
    elif a.stage == "ddladder":
        run_pool(ddladder_jobs(), a.workers, "ddladder", a.mem_budget)
    elif a.stage == "n64":
        # the 2026-09-13 20:xx UTC decision: the n48 rung dropped; only the three n64 main rows (resumed from their checkpoints)
        run_pool([j for j in job_list("main", gamma_r20_0()) if j["n"] == 64], a.workers, "n64", a.mem_budget)
    elif a.stage == "collect":
        collect()
    else:
        run_pool(job_list(a.stage, gamma_r20_0()), a.workers, a.stage, a.mem_budget)


if __name__ == "__main__":
    main()
