"""M5.32 R3 arm (ii): RELAXED two-defect boost-dressed pairs on the pinned
lattice, the Newton instrument for the covariant lambda-family.

EQUATIONS FIRST
---------------
Field M(x): real symmetric 4x4 per cell, eta = diag(-1, 1, 1, 1), jets
A_i = d_i M on the certified sym stencil (1/2 (fwd + bwd), density per
branch), F_ij = A_i eta A_j - A_j eta A_i. The candidate (R2, audited
vacuum-stable for lambda >= 1/2, static sector = I1 exactly):
    L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4
    E_lambda[M] = 4 h^3 sum_br wt sum_cells sum_{i<j} q_lambda(F_ij) + V4
    q_lambda(F) = (1 - lambda) tr(eta F eta F^T) + lambda tr(h F h F^T)
    h(M) = eta + 2 (eta u)(eta u)^T,  u the timelike unit eigenvector of
           M eta (u^T eta u = -1)
    V4 = w sum_p (tr((M eta)^p) - C_p)^2   (conjugation invariant)
The energy and its exact gradient are IMPORTED from the R2 arm (b)
instrument (m5_32_r2_b_bounded.energy_grad: the I1 part is the certified
B3.grad, the I1_h part carries the eigenframe piece through h(M); FD-gated
at 4.2e-7 there and re-gated here in the `gate` stage on the PAIR seeds).

Seeds (the M5.21.14 dressed electron, twice):
  3x3 part: the certified two-center composition of m5_21_4_a_pair
  (seed_pair: 'same' = the inverse-stereographic product ansatz with the
  +z escape tube, 'anti' = the mirror-hedgehog sum, 'single' = one
  hedgehog), isotropic-blended cores (r_c = 4), on the (1, delta, 0)
  spectrum manifold; embedded with M_00 = g (B3.embed34), so the seed's
  spectrum of M eta is the vacuum (-g, 1, delta, 0) away from the cores.
  4x4 part: per center c at z_c the boost dressing of the record,
      Q_c(x) = exp(b*(r_c) n_c . K) = I + sinh(b*) K_c + (cosh(b*) - 1) K2_c,
      r_c = |x - z_c e_z|,  n_c = (x - z_c e_z)/r_c,
  with b*(r) the M5.21.14 record profile (m5_21_14_minimize.json, rs /
  b_star, np.interp, held at its last value beyond r = 24). The PAIR
  dressing is the ORDERED PRODUCT of the two per-cell boosts,
      Q(x) = Q_top(x) Q_bot(x),   M = Q M_emb Q^T   (sym4),
  the single is Q_0 alone. Q is a per-cell Lorentz transform, so V4 is
  EXACTLY the undressed seed's V4 in every cell (conjugation invariance);
  the noncommutativity |Q_top Q_bot - Q_bot Q_top| and the energy of the
  reversed order are reported (the `seed` stage). Because b* does NOT
  decay (b*(24) = +0.025), the two dressings overlap everywhere and the
  pair's far field carries the rapidity ~2 b*; this is the record's
  profile, used as is (never tuned), and its d-independent bulk offset is
  absorbed by the constant A of the fits.

Relaxation (equal-depth heal): the arm (b) FIRE with energy-monotone
backtracking (dt0 0.02, dt_max 0.2, alpha 0.1, dt halved on a rejected
step, dt_min 1e-7), the pinned Dirichlet shell B3.pin_shell depth 1.6
(the record's pin depth) held at the SEED values (vacuum-pinned frame
boundaries: the boundary cells keep their seed dressing), a FIXED budget
of STEPS_ACC accepted steps (iteration cap IT_CAP). The dressing amplitude
is FREE (the point of the arm); it is recorded per center after the heal
as max |M_0i| inside a ball of radius BALL_R about each core, and on the
whole grid. Kill rules (never fitted): RUNAWAY when the whole-grid
max |M_0i| exceeds RUNAWAY_FACTOR x its seed value (the step recorded);
DIVERGED when E is non-finite or below the dive floor; LOCUS-HIT when the
flip metric loses its eigenframe (arm (b) rule).

Reads:
    E_int(d) = E(pair) - 2 E(single)          (same protocol, same box)
    static part      = E_int of the UNDRESSED pair (lambda-independent:
                       on block-diagonal fields I1_h = I1 and the eigen-
                       frame piece vanishes; run at lambda = 0, the
                       identity checked at lambda = 1 in the null control)
    dressing part    = E_int(dressed) - E_int(undressed)   (G2's observable)
    block split      = per field, q_lambda(F) with the time row / column
                       of F zeroed (the 3x3 sector) vs the remainder
    force sign: attraction <=> E_int INCREASES with d (F = -dE/dd < 0)
    fits over d: A + B/d (B > 0 repulsive) and A + B/d + C ln(d)/d
Controls (pre-registered): (a) lambda = 0 calibration (repulsive dressing
part or the runaway); (b) vacuum null (undressed pair == dressed pair at
amplitude 0, to roundoff; lambda = 1 == lambda = 0 on undressed fields);
(c) the anti pair attracts in its static part; (d) the mutation: the
dressing-part sign flips between lambda = 0 and lambda = 1; (e) the
reference-matched Dirichlet Poisson comparison (the M9.2 DST-I solver,
import I6) with two Gaussian sources of the dressing's rms width in the
same box, read as forces (image correction = box / free-space ratio).

STAGES (python3 m5_32_r3_ii_pair.py STAGE [--workers W] [--ladder]):
    gate      identities + the FD gradient gate on the PAIR seeds
    seed      composition report + the seed-level (UNRELAXED) E_int table
    relax     the heal jobs (n = 32, L = 48; --ladder adds n = 48, L = 72)
    collect   tables, fits, controls, the Poisson reference, plots
Out: ../data/m5_32_r3_pair.json (partials after every relaxation),
     ../data/m5_32_r3_ii/*.npz (local), ../plots/m5_32_r3_pair_*.png
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
from concurrent.futures import ProcessPoolExecutor, as_completed  # noqa
import multiprocessing as mp  # noqa: E402

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
OUT_JSON = os.path.join(DATA, "m5_32_r3_pair.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r3_ii")
M92_PATH = os.path.join(HERE, "..", "..", "..", "m9_emergent_gravity",
                        "research", "scripts", "m9_2_newton_limit.py")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


RB = _load("m5_32_r2_b_bounded", os.path.join(HERE, "m5_32_r2_b_bounded.py"))
PAIR = _load("m5_21_4_a_pair", os.path.join(HERE, "m5_21_4_a_pair.py"))
B3 = RB.B3
ETA = B3.ETA

G_MAIN, S_MAIN, DELTA = 32.0, -1.0, 0.3
LAMBDAS = (0.0, 0.75, 1.0)
DS = (10.0, 14.0, 18.0, 24.0)
N_MAIN, L_MAIN = 32, 48.0
N_LAD, L_LAD = 48, 72.0
DS_LAD = (14.0, 24.0)
STEPS_ACC = 1500
IT_CAP = 3000
BALL_R = 5.0
RUNAWAY_FACTOR = 3.0
DIVE_FLOOR = -1e6
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def cfg_of(n, L):
    return RB.cfg_of(n, L)


# ================= seeds =================
def centers_of(kind, d):
    return [0.0] if kind == "single" else [+d / 2.0, -d / 2.0]


def boost_at(cfg, zc, scale=1.0):
    """Q_c(x) = exp(scale b*(r_c) n_c . K) per cell (the arm (b) formula
    about the center z_c e_z)."""
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    Zc = Z - zc
    R = np.sqrt(X * X + Y * Y + Zc * Zc)
    nx, ny, nz = X / R, Y / R, Zc / R
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, bb in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * bb
    rs, bstar = RB.bstar_record()
    bl = scale * np.interp(R.ravel(), rs, bstar).reshape(R.shape)
    Q = (np.eye(4)[None, None, None] + np.sinh(bl)[..., None, None] * K
         + (np.cosh(bl) - 1.0)[..., None, None] * K2)
    return Q, bl


def conj(Q, M):
    return B3.sym4(np.einsum("...ab,...bc,...dc->...ad", Q, M, Q))


def seed_field(cfg, kind, d, scale, order="top_bot"):
    """(M4 seed, meta). scale = 0 is the undressed seed exactly."""
    M3 = PAIR.seed_pair(cfg, kind, d)
    M4 = B3.embed34(M3, cfg)
    meta = {"kind": kind, "d": d, "scale": scale, "order": order}
    if scale == 0.0:
        return M4, meta
    cs = centers_of(kind, d)
    Qs = [boost_at(cfg, zc, scale)[0] for zc in cs]
    if len(Qs) == 1:
        Q = Qs[0]
    else:
        Q = Qs[0] @ Qs[1] if order == "top_bot" else Qs[1] @ Qs[0]
    return conj(Q, M4), meta


def composition_report(cfg, kind, d, scale=1.0):
    """V4 exactness, noncommutativity, reversed-order energy."""
    M3 = PAIR.seed_pair(cfg, kind, d)
    Mu = B3.embed34(M3, cfg)
    Md, _ = seed_field(cfg, kind, d, scale, "top_bot")
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    h3 = cfg["h"] ** 3
    _, evu = B3.e_parts(Mu, cfg)
    _, evd = B3.e_parts(Md, cfg)
    out = {"V4_undressed": float(evu), "V4_dressed": float(evd),
           "V4_abs_diff": float(abs(evu - evd)),
           "max_abs_M0i_seed": float(np.max(np.abs(Md[..., 0, 1:])))}
    # per-cell V4 density split: overlap (between the cores) vs outside
    Me = Md @ ETA
    P = Me
    t = []
    for p in range(4):
        if p:
            P = P @ Me
        t.append(np.einsum("...kk->...", P))
    cp = B3.c4_of(cfg)
    vd = h3 * B3.W1 * sum((t[p] - cp[p]) ** 2 for p in range(4))
    if kind != "single":
        between = (np.abs(Z) < d / 2.0) & (np.sqrt(X * X + Y * Y) < d / 2.0)
        core = np.zeros_like(between)
        for zc in centers_of(kind, d):
            core |= np.sqrt(X * X + Y * Y + (Z - zc) ** 2) < 4.0
        out["V4_between_cores_excl_core_r4"] = float(np.sum(vd[between & ~core]))
        out["V4_outside_between"] = float(np.sum(vd[~between]))
        out["V4_core_balls_r4"] = float(np.sum(vd[core]))
        Qt = boost_at(cfg, +d / 2.0, scale)[0]
        Qb = boost_at(cfg, -d / 2.0, scale)[0]
        C = Qt @ Qb - Qb @ Qt
        out["noncommutativity_max_frob"] = float(
            np.max(np.sqrt(np.sum(C * C, axis=(-1, -2)))))
        Mr = conj(Qb @ Qt, Mu)
        out["max_abs_field_diff_reversed_order"] = float(np.max(np.abs(Md - Mr)))
        for lam in LAMBDAS:
            e1 = RB.energy_grad(Md, cfg, lam)[0]
            e2 = RB.energy_grad(Mr, cfg, lam)[0]
            out[f"E_lam{lam:g}_top_bot"] = float(e1)
            out[f"E_lam{lam:g}_bot_top"] = float(e2)
    return out


# ================= reads =================
def block_reads(M, cfg, lam):
    """E_stat_lambda split into the 3x3 sector (F with its time row and
    column zeroed) and the remainder (time-row / dressing part)."""
    h3, h = cfg["h"] ** 3, cfg["h"]
    u0, *_ = RB.tl_eig(M)
    hh = RB.h_of(u0)
    tot = spa = 0.0
    for br, wt in B3.branches(cfg["stencil"]):
        A = [B3.d1(M, ax, h, br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                Fs = F.copy()
                Fs[..., 0, :] = 0.0
                Fs[..., :, 0] = 0.0
                for FF, acc in ((F, "tot"), (Fs, "spa")):
                    q = 0.0
                    if lam != 1.0:
                        q = q + (1 - lam) * np.sum(RB.q_eta(FF))
                    if lam != 0.0:
                        q = q + lam * np.sum(RB.q_h(FF, hh))
                    if acc == "tot":
                        tot += wt * q
                    else:
                        spa += wt * q
    _, ev = B3.e_parts(M, cfg)
    return {"E_total": float(4 * h3 * tot + ev), "E_curv": float(4 * h3 * tot),
            "E_curv_3x3": float(4 * h3 * spa),
            "E_curv_timerow": float(4 * h3 * (tot - spa)), "V4": float(ev)}


def dressing_amplitude(M, cfg, kind, d):
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    m0i = np.sqrt(np.sum(M[..., 0, 1:] ** 2, axis=-1))
    out = {"grid_max_abs_M0i": float(np.max(np.abs(M[..., 0, 1:]))),
           "grid_max_norm_M0i": float(np.max(m0i))}
    for lab, zc in zip(("top", "bot"), centers_of(kind, d)):
        ball = np.sqrt(X * X + Y * Y + (Z - zc) ** 2) <= BALL_R
        out[f"{lab}_ball_max_norm_M0i"] = float(np.max(m0i[ball]))
        out[f"{lab}_ball_mean_norm_M0i"] = float(np.mean(m0i[ball]))
    return out


def amp_trend(a_seed, a_end):
    r = a_end / max(a_seed, 1e-300)
    if a_seed == 0.0:
        return "zero (undressed)" if a_end < 1e-10 else "appeared"
    if r < 0.05:
        return "vanished"
    if r < 0.9:
        return "shrank"
    if r > 1.1:
        return "grew"
    return "held"


# ================= descent (the arm (b) loop, seed-agnostic) =================
def descend(M0, cfg, lam, steps_acc, it_cap, tag, log_every=100,
            dt0=0.02, dt_max=0.2):
    free = (~B3.pin_shell(cfg["n"], cfg["h"]))[..., None, None].astype(float)
    M = M0.copy()
    E0, G, info = RB.energy_grad(M, cfg, lam)
    m0i_seed = float(np.max(np.abs(M0[..., 0, 1:])))
    out = {"E0": float(E0), "steps_acc_budget": steps_acc, "it_cap": it_cap,
           "pin": "B3.pin_shell depth 1.6 (Dirichlet at the seed values)",
           "fire": {"dt0": dt0, "dt_max": dt_max, "alpha0": 0.1,
                    "dt_min": 1e-7},
           "max_abs_M0i_seed": m0i_seed, "trace": []}
    if not np.isfinite(E0) or G is None:
        out.update(stop="DIVERGED (seed energy undefined)", verdict="DIVERGED",
                   steps_run=0, accepted=0)
        return M, out
    v = np.zeros_like(M)
    dt, alpha, n_up = dt0, 0.1, 0
    dt_min = 1e-7
    F = -G * free
    E_prev = E0
    stop = "budget"
    n_rej, n_rej_locus, n_acc = 0, 0, 0
    fmax = float(np.max(np.abs(F)))
    it = 0
    runaway = None
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
            alpha = 0.1
            n_up = 0
        v_try = v + dt * F
        M_try = M + dt * v_try
        E, G, info = RB.energy_grad(M_try, cfg, lam)
        locus_loss = not info["ok"]
        reject = locus_loss or not np.isfinite(E) or \
            E > E_prev + 1e-12 * max(abs(E_prev), 1.0)
        if reject:
            n_rej += 1
            if locus_loss:
                n_rej_locus += 1
            dt *= 0.5
            v[:] = 0.0
            alpha, n_up = 0.1, 0
            if dt < dt_min:
                stop = "LOCUS-HIT" if locus_loss else \
                    "STALLED (dt collapsed, no descent direction accepted)"
                break
            continue
        n_acc += 1
        M, v, E_prev = M_try, v_try, E
        F = -G * free
        fmax = float(np.max(np.abs(F)))
        m0i = info["max_abs_M0i"]
        if m0i_seed > 0 and m0i > RUNAWAY_FACTOR * m0i_seed:
            runaway = {"step": it, "accepted": n_acc, "max_abs_M0i": m0i,
                       "E": float(E)}
            stop = "RUNAWAY"
            break
        if E < DIVE_FLOOR:
            stop = "DIVERGED (dive floor)"
            break
        if it % log_every == 0 or n_acc == steps_acc:
            row = {"it": it, "acc": n_acc, "E": float(E), "fmax": fmax,
                   "dt": dt, "min_gap": info["min_gap"], "max_abs_M0i": m0i}
            out["trace"].append(row)
            log(f"{tag} it {it:5d} acc {n_acc:5d} E {E:14.5f} fmax {fmax:.3e}"
                f" dt {dt:.2e} max|M0i| {m0i:.4g} rej {n_rej}")
    if stop == "budget" and n_acc < steps_acc:
        stop = f"IT_CAP ({it_cap} iterations before {steps_acc} accepted)"
    out.update({"stop": stop, "steps_run": it, "accepted": n_acc,
                "rejected": n_rej, "rejected_locus": n_rej_locus,
                "dt_final": dt, "fmax_end": fmax, "E_end": float(E_prev),
                "E_drop": float(E0 - E_prev), "runaway": runaway,
                "max_abs_M0i_end": float(np.max(np.abs(M[..., 0, 1:])))})
    tr = out["trace"]
    q = [r for r in tr if r["acc"] >= 0.75 * n_acc] if n_acc else []
    if stop.startswith("budget") or stop.startswith("IT_CAP"):
        if len(q) >= 2:
            dE = q[-1]["E"] - q[0]["E"]
            out["last_quarter_dE"] = float(dE)
            out["last_quarter_rel"] = float(abs(dE) / max(abs(q[-1]["E"]), 1.0))
            out["verdict"] = ("PLATEAU" if out["last_quarter_rel"] <= 1e-3
                              else "FALLING (still descending at the budget)"
                              if dE < 0 else "RISING")
        else:
            out["verdict"] = "budget (trace too short)"
    else:
        out["verdict"] = stop
    return M, out


# ================= jobs =================
def job_tag(lam, kind, d, scale, n):
    return f"lam{lam:g}_{'dr' if scale else 'un'}{scale:g}_{kind}_d{d:g}_n{n}"


def run_job(args):
    lam, kind, d, scale, n, L, steps = args
    t0 = time.time()
    cfg = cfg_of(n, L)
    tag = job_tag(lam, kind, d, scale, n)
    row = {"lam": lam, "kind": kind, "d": d, "scale": scale, "n": n, "L": L,
           "h": cfg["h"], "steps_acc_budget": steps, "tag": tag,
           "dressed": bool(scale != 0.0)}
    try:
        M0, meta = seed_field(cfg, kind, d, scale)
        row["seed_reads"] = block_reads(M0, cfg, lam)
        row["seed_amp"] = dressing_amplitude(M0, cfg, kind, d)
        M, des = descend(M0, cfg, lam, steps, IT_CAP, tag)
        row["descent"] = des
        finite = bool(np.all(np.isfinite(M))) and np.isfinite(des["E_end"])
        row["status"] = ("OK" if finite and des["stop"].startswith("budget")
                         else des["stop"].split(" ")[0] if not des["stop"].startswith("budget")
                         else "OK")
        if finite:
            row["end_reads"] = block_reads(M, cfg, lam)
            row["E"] = row["end_reads"]["E_total"]
            row["end_amp"] = dressing_amplitude(M, cfg, kind, d)
            row["amp_trend"] = {
                lab: amp_trend(row["seed_amp"].get(f"{lab}_ball_max_norm_M0i", 0.0),
                               row["end_amp"].get(f"{lab}_ball_max_norm_M0i", 0.0))
                for lab in ("top", "bot") if f"{lab}_ball_max_norm_M0i" in row["end_amp"]}
            row["amp_trend"]["grid"] = amp_trend(row["seed_amp"]["grid_max_norm_M0i"],
                                                 row["end_amp"]["grid_max_norm_M0i"])
            try:
                dq = d if d > 0 else 18.0
                row["charge"] = PAIR.charge_suite(M[..., 1:, 1:], cfg, dq)
            except Exception as e:                        # noqa: BLE001
                row["charge"] = f"charge suite failed: {e!r}"
            os.makedirs(OUT_NPZ, exist_ok=True)
            np.savez_compressed(os.path.join(OUT_NPZ, f"{tag}.npz"),
                                M=M.astype(np.float64))
        else:
            row["E"] = None
    except Exception as e:                                # noqa: BLE001
        row["status"] = "DIVERGED"
        row["stop"] = f"exception: {e!r}"
        row["E"] = None
    row["wall_s"] = round(time.time() - t0, 1)
    log(f"DONE {tag} status {row['status']} E {row.get('E')} wall {row['wall_s']}")
    return row


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R3 arm (ii): relaxed boost-dressed pairs",
            "candidate": "L_lambda = -4[(1-lambda) I1 + lambda I1_h] - V4",
            "point": {"g": G_MAIN, "s": S_MAIN, "delta": DELTA},
            "lambdas": list(LAMBDAS), "ds": list(DS),
            "protocol": {"steps_acc": STEPS_ACC, "it_cap": IT_CAP,
                         "pin_depth": 1.6, "ball_r": BALL_R,
                         "runaway_factor": RUNAWAY_FACTOR,
                         "dive_floor": DIVE_FLOOR}}


def save_json(J):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1, default=float)
    os.replace(tmp, OUT_JSON)


# ================= GATE =================
def fd_gate(M, cfg, lam, rng, ndir=2, eps=1e-3):
    E0, G, info = RB.energy_grad(M, cfg, lam)
    rows = []
    for _ in range(ndir):
        D = B3.sym4(rng.standard_normal(M.shape))
        D /= np.sqrt(np.sum(D * D))
        gd = float(np.sum(G * D))

        def e_at(t):
            return RB.energy_grad(M + t * D, cfg, lam)[0]

        def fd4(e):
            return (8 * (e_at(e) - e_at(-e)) - (e_at(2 * e) - e_at(-2 * e))) / (12 * e)
        f1, f2 = fd4(eps), fd4(eps / 2)
        rich = (16 * f2 - f1) / 15.0
        rows.append({"g_dot_D": gd, "fd_richardson": rich,
                     "rel_err": float(abs(rich - gd) / max(abs(gd), 1e-300))})
    return {"E": float(E0), "min_gap": info["min_gap"], "dirs": rows}


def stage_gate():
    J = load_json()
    rng = np.random.default_rng(3230)
    out = {"checks": {}, "fd": {}}
    chk = out["checks"]
    cfg16 = cfg_of(16, 24.0)
    cfg32 = cfg_of(N_MAIN, L_MAIN)
    # identities on the undressed PAIR seed: lambda = 0 == B3, lambda = 1 == lambda = 0
    for cfg, lab in ((cfg16, "n16"), (cfg32, "n32")):
        Mu, _ = seed_field(cfg, "same", 10.0 if lab == "n16" else 18.0, 0.0)
        e0, g0, _ = RB.energy_grad(Mu, cfg, 0.0)
        e1, g1, _ = RB.energy_grad(Mu, cfg, 1.0)
        eb, gb = B3.e_total(Mu, cfg), B3.grad(Mu, cfg)
        chk[f"{lab}_undressed_E_lam0_vs_B3"] = float(abs(e0 - eb) / abs(eb))
        chk[f"{lab}_undressed_grad_lam0_vs_B3"] = float(np.max(np.abs(g0 - gb)) / np.max(np.abs(gb)))
        chk[f"{lab}_undressed_E_lam1_vs_lam0"] = float(abs(e1 - e0) / abs(e0))
        chk[f"{lab}_undressed_grad_lam1_vs_lam0"] = float(np.max(np.abs(g1 - g0)) / np.max(np.abs(g0)))
        # vacuum null at the seed: dressed with scale 0 == undressed exactly
        Mz, _ = seed_field(cfg, "same", 10.0 if lab == "n16" else 18.0, 0.0)
        chk[f"{lab}_scale0_field_diff"] = float(np.max(np.abs(Mz - Mu)))
    # the FD gradient gate on the DRESSED pair seeds
    worst = 0.0
    for lam in (0.75, 1.0):
        for kind, d in (("same", 10.0), ("anti", 10.0)):
            Md, _ = seed_field(cfg16, kind, d, 1.0)
            r = fd_gate(Md, cfg16, lam, rng, ndir=3)
            out["fd"][f"n16_lam{lam:g}_{kind}_d{d:g}"] = r
            worst = max(worst, max(x["rel_err"] for x in r["dirs"]))
    Md, _ = seed_field(cfg32, "same", 18.0, 1.0)
    r = fd_gate(Md, cfg32, 1.0, rng, ndir=2)
    out["fd"]["n32_lam1_same_d18"] = r
    worst = max(worst, max(x["rel_err"] for x in r["dirs"]))
    out["fd_worst_rel_err"] = float(worst)
    out["fd_gate_arm_b_record"] = 4.2422932739602266e-07
    out["gate_pass"] = bool(worst <= 1e-6 and all(v <= 1e-9 for v in chk.values()))
    J["gate"] = out
    save_json(J)
    log(f"GATE worst FD rel {worst:.3e}; checks " + ", ".join(f"{k} {v:.1e}" for k, v in chk.items()))
    log(f"GATE PASS = {out['gate_pass']}")
    return out


# ================= SEED =================
def stage_seed():
    J = load_json()
    cfg = cfg_of(N_MAIN, L_MAIN)
    out = {"n": N_MAIN, "L": L_MAIN, "h": cfg["h"], "composition": {},
           "seed_energies": {}, "seed_Eint": {}}
    rs, bstar = RB.bstar_record()
    out["bstar_profile"] = {"r_min": float(rs[0]), "r_max": float(rs[-1]),
                            "b_at_rmin": float(bstar[0]), "b_max": float(bstar.max()),
                            "r_of_bmax": float(rs[bstar.argmax()]),
                            "b_min": float(bstar.min()), "b_at_rmax": float(bstar[-1]),
                            "decays": bool(abs(bstar[-1]) < 0.1 * bstar.max())}
    for kind in ("single", "same", "anti"):
        for d in ((0.0,) if kind == "single" else DS):
            out["composition"][f"{kind}_d{d:g}"] = composition_report(cfg, kind, d)
            log(f"SEED comp {kind} d{d:g}: " + json.dumps(out["composition"][f"{kind}_d{d:g}"])[:300])
    # seed-level energies per lambda (the UNRELAXED read)
    for scale in (0.0, 1.0):
        for kind in ("single", "same", "anti"):
            for d in ((0.0,) if kind == "single" else DS):
                M, _ = seed_field(cfg, kind, d, scale)
                for lam in LAMBDAS:
                    key = job_tag(lam, kind, d, scale, N_MAIN)
                    out["seed_energies"][key] = block_reads(M, cfg, lam)
    for lam in LAMBDAS:
        for scale in (0.0, 1.0):
            es = out["seed_energies"][job_tag(lam, "single", 0.0, scale, N_MAIN)]
            for kind in ("same", "anti"):
                for d in DS:
                    ep = out["seed_energies"][job_tag(lam, kind, d, scale, N_MAIN)]
                    out["seed_Eint"][job_tag(lam, kind, d, scale, N_MAIN)] = {
                        k: ep[k] - 2 * es[k] for k in ep}
    J["seed"] = out
    save_json(J)
    for lam in LAMBDAS:
        for kind in ("same", "anti"):
            line = " ".join(
                f"d{d:g}: un {out['seed_Eint'][job_tag(lam, kind, d, 0.0, N_MAIN)]['E_total']:+.4f}"
                f" dr {out['seed_Eint'][job_tag(lam, kind, d, 1.0, N_MAIN)]['E_total']:+.4f}"
                for d in DS)
            log(f"SEED E_int lam {lam:g} {kind}: {line}")
    return out


# ================= RELAX =================
LADDER_KINDS = tuple(os.environ.get("R3_LADDER_KINDS", "same,anti").split(","))


def job_list(ladder=False):
    jobs = []
    # undressed (static part), lambda = 0, shared by every lambda
    jobs.append((0.0, "single", 0.0, 0.0, N_MAIN, L_MAIN, STEPS_ACC))
    for kind in ("same", "anti"):
        for d in DS:
            jobs.append((0.0, kind, d, 0.0, N_MAIN, L_MAIN, STEPS_ACC))
    # the vacuum-null / identity controls at lambda = 1 on undressed fields
    jobs.append((1.0, "same", 18.0, 0.0, N_MAIN, L_MAIN, STEPS_ACC))
    jobs.append((1.0, "single", 0.0, 0.0, N_MAIN, L_MAIN, STEPS_ACC))
    # dressed, every lambda
    for lam in LAMBDAS:
        jobs.append((lam, "single", 0.0, 1.0, N_MAIN, L_MAIN, STEPS_ACC))
        for kind in ("same", "anti"):
            for d in DS:
                jobs.append((lam, kind, d, 1.0, N_MAIN, L_MAIN, STEPS_ACC))
    if ladder:
        for scale in (0.0, 1.0):
            jobs.append((1.0, "single", 0.0, scale, N_LAD, L_LAD, STEPS_ACC))
            for kind in LADDER_KINDS:
                for d in DS_LAD:
                    jobs.append((1.0, kind, d, scale, N_LAD, L_LAD, STEPS_ACC))
    return jobs


def stage_relax(workers=6, ladder=False, only_ladder=False):
    J = load_json()
    rows = {r["tag"]: r for r in J.get("rows", [])}
    jobs = [j for j in job_list(ladder) if job_tag(j[0], j[1], j[2], j[3], j[4]) not in rows]
    if only_ladder:
        jobs = [j for j in jobs if j[4] == N_LAD]
    log(f"RELAX {len(jobs)} jobs on {workers} workers (done already: {len(rows)})")
    ctx = mp.get_context("spawn")
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as ex:
        futs = {ex.submit(run_job, j): j for j in jobs}
        for fut in as_completed(futs):
            j = futs[fut]
            try:
                r = fut.result()
            except Exception as e:                        # noqa: BLE001
                r = {"lam": j[0], "kind": j[1], "d": j[2], "scale": j[3], "n": j[4],
                     "L": j[5], "tag": job_tag(j[0], j[1], j[2], j[3], j[4]),
                     "status": "DIVERGED", "stop": f"exception: {e!r}", "E": None}
            rows[r["tag"]] = r
            J = load_json()
            J["rows"] = list(rows.values())
            J["relax_wall_s"] = round(time.time() - t0, 1)
            save_json(J)
            log(f"RELAX [{len(rows)}] saved {r['tag']}")
    return rows


# ================= COLLECT =================
def fit_inv(ds, es):
    ds, es = np.asarray(ds, float), np.asarray(es, float)
    X = np.stack([np.ones_like(ds), 1.0 / ds], axis=1)
    c, *_ = np.linalg.lstsq(X, es, rcond=None)
    pred = X @ c
    ss = np.sum((es - es.mean()) ** 2)
    r2 = 1.0 - np.sum((es - pred) ** 2) / ss if ss > 0 else float("nan")
    return {"A": float(c[0]), "B": float(c[1]), "R2": float(r2)}


def fit_log(ds, es):
    ds, es = np.asarray(ds, float), np.asarray(es, float)
    X = np.stack([np.ones_like(ds), 1.0 / ds, np.log(ds) / ds], axis=1)
    c, *_ = np.linalg.lstsq(X, es, rcond=None)
    pred = X @ c
    ss = np.sum((es - es.mean()) ** 2)
    r2 = 1.0 - np.sum((es - pred) ** 2) / ss if ss > 0 else float("nan")
    return {"A": float(c[0]), "B": float(c[1]), "C": float(c[2]), "R2": float(r2)}


def force_read(ds, es):
    """attraction <=> E_int increases with d over the outer window."""
    if len(ds) < 2:
        return {"sign": "n/a"}
    o = np.argsort(ds)
    ds, es = np.asarray(ds)[o], np.asarray(es)[o]
    dE = es[-1] - es[-2]
    return {"outer_window": [float(ds[-2]), float(ds[-1])],
            "dEint_dd_outer": float(dE / (ds[-1] - ds[-2])),
            "sign": "ATTRACTIVE" if dE > 0 else "REPULSIVE" if dE < 0 else "flat",
            "monotone": "increasing" if np.all(np.diff(es) > 0) else
            "decreasing" if np.all(np.diff(es) < 0) else "non-monotone"}


def poisson_reference(n, L, ds, sigma, label):
    """Two Gaussian sources (unit mass, width sigma) at z = +-d/2 in the
    same box, Dirichlet Phi = 0 on the boundary (the M9.2 DST-I solver);
    E_ref(d) = -int rho_1 Phi_2 (Newton sign: attractive, E increases with
    d), F_z on the top source = -int rho_1 dPhi_2/dz; free-space values
    for the same Gaussians (erf smearing) and the image ratio."""
    M92 = _load("m9_2_newton_limit", M92_PATH)
    cfg = cfg_of(n, L)
    X, Y, Z = B3.coords(n, cfg["h"])
    h = cfg["h"]
    out = {"n": n, "L": L, "h": h, "sigma": sigma, "rows": []}
    from math import erf

    def gauss(zc):
        r2 = X * X + Y * Y + (Z - zc) ** 2
        g = np.exp(-r2 / (2 * sigma ** 2))
        return g / (np.sum(g) * h ** 3)
    for d in ds:
        r1, r2_ = gauss(+d / 2), gauss(-d / 2)
        phi2 = np.zeros_like(r2_)
        phi2[1:-1, 1:-1, 1:-1] = M92._dst_poisson((4 * np.pi * r2_)[1:-1, 1:-1, 1:-1], h)
        # Phi_2 solves lap Phi = 4 pi rho: Phi = -1/r in free space
        E_box = float(np.sum(r1 * phi2) * h ** 3)
        dphi = (np.roll(phi2, -1, axis=2) - np.roll(phi2, 1, axis=2)) / (2 * h)
        F_box = float(-np.sum(r1 * dphi) * h ** 3)
        # free space for two Gaussians of width sigma: E = -erf(d / (2 sigma)) / d
        s2 = sigma * np.sqrt(2.0)
        E_free = -erf(d / s2) / d
        F_free = -(erf(d / s2) / d ** 2 - 2 / (np.sqrt(np.pi) * s2) * np.exp(-(d / s2) ** 2) / d)
        out["rows"].append({"d": d, "E_box": E_box, "E_free": E_free,
                            "F_box": F_box, "F_free": F_free,
                            "image_ratio_E": E_box / E_free,
                            "image_ratio_F": F_box / F_free})
    es = [r["E_box"] for r in out["rows"]]
    out["fit_inv_box"] = fit_inv(ds, es)
    out["fit_inv_free"] = fit_inv(ds, [r["E_free"] for r in out["rows"]])
    out["force_read_box"] = force_read(ds, es)
    out["label"] = label
    return out


def stage_collect():
    J = load_json()
    rows = {r["tag"]: r for r in J.get("rows", [])}
    res = {"tables": {}, "fits": {}, "controls": {}}

    def E_of(lam, kind, d, scale, n):
        r = rows.get(job_tag(lam, kind, d, scale, n))
        if r is None or r.get("E") is None:
            return None, r
        return r["E"], r

    # ---- tables per box
    for n, L, dlist in ((N_MAIN, L_MAIN, DS), (N_LAD, L_LAD, DS_LAD)):
        lam_list = LAMBDAS if n == N_MAIN else (1.0,)
        # static part: the undressed pair, lambda = 0 at n32; lambda = 1 in the ladder
        lam_un = 0.0 if n == N_MAIN else 1.0
        Es_un, r_su = E_of(lam_un, "single", 0.0, 0.0, n)
        tab = {}
        for lam in lam_list:
            Es_dr, r_sd = E_of(lam, "single", 0.0, 1.0, n)
            for kind in ("same", "anti"):
                for d in dlist:
                    Ep_un, r_pu = E_of(lam_un, kind, d, 0.0, n)
                    Ep_dr, r_pd = E_of(lam, kind, d, 1.0, n)
                    key = f"lam{lam:g}_{kind}_d{d:g}_n{n}"
                    row = {"lam": lam, "kind": kind, "d": d, "n": n, "L": L,
                           "h": L / n}
                    if Ep_un is not None and Es_un is not None:
                        row["Eint_static_undressed"] = Ep_un - 2 * Es_un
                        row["static_status"] = r_pu["status"]
                        row["static_verdict"] = r_pu["descent"]["verdict"]
                        row["static_steps"] = [r_pu["descent"]["accepted"], r_pu["descent"]["steps_run"]]
                        row["static_fmax"] = r_pu["descent"]["fmax_end"]
                    if Ep_dr is not None and Es_dr is not None:
                        row["Eint_dressed_total"] = Ep_dr - 2 * Es_dr
                        row["Eint_dressed_3x3_block"] = (r_pd["end_reads"]["E_curv_3x3"] + r_pd["end_reads"]["V4"]
                                                         - 2 * (r_sd["end_reads"]["E_curv_3x3"] + r_sd["end_reads"]["V4"]))
                        row["Eint_dressed_timerow_block"] = (r_pd["end_reads"]["E_curv_timerow"]
                                                             - 2 * r_sd["end_reads"]["E_curv_timerow"])
                        row["dressed_status"] = r_pd["status"]
                        row["dressed_verdict"] = r_pd["descent"]["verdict"]
                        row["dressed_steps"] = [r_pd["descent"]["accepted"], r_pd["descent"]["steps_run"]]
                        row["dressed_fmax"] = r_pd["descent"]["fmax_end"]
                        row["dressed_E_drop"] = r_pd["descent"]["E_drop"]
                        row["single_dressed_status"] = r_sd["status"]
                        row["single_dressed_verdict"] = r_sd["descent"]["verdict"]
                        row["amp_end"] = {k: v for k, v in r_pd["end_amp"].items()}
                        row["amp_seed"] = {k: v for k, v in r_pd["seed_amp"].items()}
                        row["amp_trend"] = r_pd["amp_trend"]
                        row["single_amp_end_grid"] = r_sd["end_amp"]["grid_max_norm_M0i"]
                        row["single_amp_trend"] = r_sd["amp_trend"]
                        row["runaway"] = r_pd["descent"].get("runaway")
                        if "Eint_static_undressed" in row:
                            row["Eint_dressing_part"] = row["Eint_dressed_total"] - row["Eint_static_undressed"]
                    else:
                        row["dressed_status"] = (r_pd or {}).get("status", "missing")
                        row["single_dressed_status"] = (r_sd or {}).get("status", "missing")
                    tab[key] = row
        res["tables"][f"n{n}"] = tab
        # ---- fits and force reads
        fits = {}
        for lam in lam_list:
            for kind in ("same", "anti"):
                for q in ("Eint_static_undressed", "Eint_dressed_total", "Eint_dressing_part",
                          "Eint_dressed_timerow_block", "Eint_dressed_3x3_block"):
                    pts = [(d, tab[f"lam{lam:g}_{kind}_d{d:g}_n{n}"][q]) for d in dlist
                           if q in tab[f"lam{lam:g}_{kind}_d{d:g}_n{n}"]]
                    if len(pts) < 2:
                        continue
                    ds_, es_ = zip(*pts)
                    f = {"points": [list(p) for p in pts], "force": force_read(ds_, es_)}
                    if len(pts) >= 3:
                        f["fit_inv"] = fit_inv(ds_, es_)
                    if len(pts) >= 4:
                        f["fit_log"] = fit_log(ds_, es_)
                    fits[f"lam{lam:g}_{kind}_{q}"] = f
        res["fits"][f"n{n}"] = fits
    # ---- controls
    ctr = res["controls"]
    tab = res["tables"]["n32"]
    # (a) calibration at lambda = 0
    cal = {}
    for kind in ("same", "anti"):
        rr = [tab[f"lam0_{kind}_d{d:g}_n32"] for d in DS]
        cal[kind] = {"dressing_part": {f"d{r['d']:g}": r.get("Eint_dressing_part") for r in rr},
                     "runaway": {f"d{r['d']:g}": r.get("runaway") for r in rr},
                     "status": {f"d{r['d']:g}": r.get("dressed_status") for r in rr},
                     "force": res["fits"]["n32"].get(f"lam0_{kind}_Eint_dressing_part", {}).get("force")}
    single0 = rows.get(job_tag(0.0, "single", 0.0, 1.0, N_MAIN), {})
    cal["single_dressed_lam0"] = {"status": single0.get("status"),
                                  "verdict": single0.get("descent", {}).get("verdict"),
                                  "E0": single0.get("descent", {}).get("E0"),
                                  "E_end": single0.get("descent", {}).get("E_end"),
                                  "runaway": single0.get("descent", {}).get("runaway"),
                                  "amp_trend": single0.get("amp_trend")}
    f0 = cal["same"]["force"]
    cal["verdict"] = ("RUNAWAY shown" if any(v for v in cal["same"]["runaway"].values())
                      else "dressing part REPULSIVE (calibration met)" if f0 and f0["sign"] == "REPULSIVE"
                      else "dressing part ATTRACTIVE at lambda = 0 (calibration NOT met as stated)"
                      if f0 else "incomplete")
    ctr["a_calibration_lam0"] = cal
    # (b) vacuum null
    nul = {}
    a = rows.get(job_tag(0.0, "same", 18.0, 0.0, N_MAIN))
    b = rows.get(job_tag(1.0, "same", 18.0, 0.0, N_MAIN))
    s0 = rows.get(job_tag(0.0, "single", 0.0, 0.0, N_MAIN))
    s1 = rows.get(job_tag(1.0, "single", 0.0, 0.0, N_MAIN))
    if a and b and a.get("E") is not None and b.get("E") is not None:
        nul["same_d18_E_lam0"] = a["E"]
        nul["same_d18_E_lam1"] = b["E"]
        nul["abs_diff"] = abs(a["E"] - b["E"])
        try:
            Ma = np.load(os.path.join(OUT_NPZ, a["tag"] + ".npz"))["M"]
            Mb = np.load(os.path.join(OUT_NPZ, b["tag"] + ".npz"))["M"]
            nul["field_max_abs_diff"] = float(np.max(np.abs(Ma - Mb)))
            nul["max_abs_M0i_lam1_end"] = float(np.max(np.abs(Mb[..., 0, 1:])))
        except Exception as e:                            # noqa: BLE001
            nul["field_max_abs_diff"] = f"unavailable: {e!r}"
    if s0 and s1 and s0.get("E") is not None and s1.get("E") is not None:
        nul["single_E_lam0"] = s0["E"]
        nul["single_E_lam1"] = s1["E"]
        nul["single_abs_diff"] = abs(s0["E"] - s1["E"])
    nul["seed_scale0_field_diff"] = J.get("gate", {}).get("checks", {}).get("n32_scale0_field_diff")
    nul["note"] = ("a dressed seed at amplitude 0 IS the undressed seed (same function, "
                   "scale = 0 branch), so the null is the lambda = 1 vs lambda = 0 identity "
                   "on undressed fields after the heal")
    nul["pass_1e-8"] = bool(nul.get("abs_diff", 1) <= 1e-8 and nul.get("single_abs_diff", 1) <= 1e-8)
    ctr["b_vacuum_null"] = nul
    # (c) the anti pair attracts in its static part
    fa = res["fits"]["n32"].get("lam0_anti_Eint_static_undressed", {})
    fs = res["fits"]["n32"].get("lam0_same_Eint_static_undressed", {})
    ctr["c_static_coulomb"] = {
        "anti": fa.get("force"), "anti_fit_inv": fa.get("fit_inv"),
        "same": fs.get("force"), "same_fit_inv": fs.get("fit_inv"),
        "same_Eint_values": fs.get("points"),
        "R1b_like_charge_control_it120": "E_int(same) ~ 5.05..5.13, nearly flat (INS 3x3 stack, T2 potential, it = 120)",
        "pass": bool(fa.get("force", {}).get("sign") == "ATTRACTIVE")}
    # (d) the mutation: dressing-part sign flip lambda = 0 -> lambda = 1
    mut = {}
    for kind in ("same", "anti"):
        per_d = {}
        for d in DS:
            r0 = tab[f"lam0_{kind}_d{d:g}_n32"].get("Eint_dressing_part")
            r1 = tab[f"lam1_{kind}_d{d:g}_n32"].get("Eint_dressing_part")
            r75 = tab[f"lam0.75_{kind}_d{d:g}_n32"].get("Eint_dressing_part")
            per_d[f"d{d:g}"] = {"lam0": r0, "lam0.75": r75, "lam1": r1,
                                "value_sign_flips_0_to_1": (None if r0 is None or r1 is None
                                                            else bool(np.sign(r0) != np.sign(r1)))}
        f0 = res["fits"]["n32"].get(f"lam0_{kind}_Eint_dressing_part", {}).get("force", {})
        f1 = res["fits"]["n32"].get(f"lam1_{kind}_Eint_dressing_part", {}).get("force", {})
        mut[kind] = {"per_d": per_d, "force_sign_lam0": f0.get("sign"),
                     "force_sign_lam1": f1.get("sign"),
                     "force_sign_flips": bool(f0.get("sign") and f1.get("sign")
                                              and f0["sign"] != f1["sign"]
                                              and "flat" not in (f0["sign"], f1["sign"]))}
    ctr["d_mutation"] = mut
    # (e) the Poisson reference: sigma = rms radius of |M_0i|^2 of the relaxed dressed single
    sig = {}
    for lam in LAMBDAS:
        r = rows.get(job_tag(lam, "single", 0.0, 1.0, N_MAIN))
        if r and r.get("E") is not None:
            try:
                Mz = np.load(os.path.join(OUT_NPZ, r["tag"] + ".npz"))["M"]
                cfg = cfg_of(N_MAIN, L_MAIN)
                X, Y, Z = B3.coords(N_MAIN, cfg["h"])
                w = np.sum(Mz[..., 0, 1:] ** 2, axis=-1)
                sig[f"lam{lam:g}"] = float(np.sqrt(np.sum((X * X + Y * Y + Z * Z) * w) / np.sum(w)))
            except Exception:                              # noqa: BLE001
                pass
    Ms, _ = seed_field(cfg_of(N_MAIN, L_MAIN), "single", 0.0, 1.0)
    X, Y, Z = B3.coords(N_MAIN, L_MAIN / N_MAIN)
    w = np.sum(Ms[..., 0, 1:] ** 2, axis=-1)
    sig["seed"] = float(np.sqrt(np.sum((X * X + Y * Y + Z * Z) * w) / np.sum(w)))
    pois = {"sigma_rms_of_M0i_sq": sig}
    sig_use = sig.get("lam1", sig["seed"])
    pois["n32_L48_sigma_rms"] = poisson_reference(N_MAIN, L_MAIN, DS, sig_use, "rms width of the dressing")
    pois["n32_L48_sigma_2h"] = poisson_reference(N_MAIN, L_MAIN, DS, 2 * L_MAIN / N_MAIN, "point-like (2h)")
    pois["n48_L72_sigma_rms"] = poisson_reference(N_LAD, L_LAD, DS, sig_use, "rms width of the dressing, ladder box")
    pois["n48_L72_sigma_2h"] = poisson_reference(N_LAD, L_LAD, DS, 2 * L_LAD / N_LAD, "point-like (2h), ladder box")
    ctr["e_poisson_reference"] = pois
    # ---- job summary
    res["job_summary"] = {t: {"status": r.get("status"), "E": r.get("E"),
                              "verdict": r.get("descent", {}).get("verdict"),
                              "accepted": r.get("descent", {}).get("accepted"),
                              "steps_run": r.get("descent", {}).get("steps_run"),
                              "fmax_end": r.get("descent", {}).get("fmax_end"),
                              "E_drop": r.get("descent", {}).get("E_drop"),
                              "runaway": r.get("descent", {}).get("runaway"),
                              "wall_s": r.get("wall_s")} for t, r in sorted(rows.items())}
    J["results"] = res
    J["collected_utc"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
    save_json(J)
    plot(J)
    # terminal tables
    for n in ("n32", "n48"):
        for k, r in res["tables"][n].items():
            log(f"{k}: static {r.get('Eint_static_undressed')} dressed {r.get('Eint_dressed_total')} "
                f"dressing_part {r.get('Eint_dressing_part')} [{r.get('dressed_verdict')}] "
                f"amp grid seed {r.get('amp_seed', {}).get('grid_max_norm_M0i')} -> end "
                f"{r.get('amp_end', {}).get('grid_max_norm_M0i')} trend {r.get('amp_trend')}")
    for n in ("n32", "n48"):
        for k, f in res["fits"][n].items():
            log(f"FIT {n} {k}: force {f['force'].get('sign')} inv {f.get('fit_inv')} log {f.get('fit_log')}")
    log(f"CONTROLS: a {ctr['a_calibration_lam0']['verdict']}; b {ctr['b_vacuum_null']}; "
        f"c anti {ctr['c_static_coulomb']['anti']}; d {json.dumps(ctr['d_mutation'])[:600]}")
    return res


def plot(J):
    res = J["results"]
    os.makedirs(PLOTS, exist_ok=True)
    tab = res["tables"]["n32"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax, q, ttl in zip(axes, ("Eint_static_undressed", "Eint_dressed_total", "Eint_dressing_part"),
                          ("static (undressed) E_int", "dressed E_int", "dressing part = dressed - undressed")):
        for lam, c in zip(LAMBDAS, ("k", "tab:orange", "tab:red")):
            for kind, mk in (("same", "o-"), ("anti", "s--")):
                pts = [(d, tab[f"lam{lam:g}_{kind}_d{d:g}_n32"].get(q)) for d in DS]
                pts = [p for p in pts if p[1] is not None]
                if pts:
                    ax.plot(*zip(*pts), mk, color=c, label=f"lam {lam:g} {kind}")
        ax.set_title(ttl, fontsize=10)
        ax.set_xlabel("d")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("E_int = E(pair) - 2 E(single)  (n32 L48 h1.5)")
    axes[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r3_pair_eint.png"), dpi=130)
    plt.close(fig)
    # dressing amplitudes
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for lam, c in zip(LAMBDAS, ("k", "tab:orange", "tab:red")):
        for kind, mk in (("same", "o-"), ("anti", "s--")):
            pts = [(d, tab[f"lam{lam:g}_{kind}_d{d:g}_n32"].get("amp_end", {}).get("top_ball_max_norm_M0i"))
                   for d in DS]
            pts = [p for p in pts if p[1] is not None]
            if pts:
                ax[0].plot(*zip(*pts), mk, color=c, label=f"lam {lam:g} {kind}")
            pts = [(d, tab[f"lam{lam:g}_{kind}_d{d:g}_n32"].get("amp_end", {}).get("grid_max_norm_M0i"))
                   for d in DS]
            pts = [p for p in pts if p[1] is not None]
            if pts:
                ax[1].plot(*zip(*pts), mk, color=c, label=f"lam {lam:g} {kind}")
    s = tab[f"lam1_same_d18_n32"].get("amp_seed", {})
    for a, k, t in ((ax[0], "top_ball_max_norm_M0i", "core ball r <= 5: max |M_0i| after heal"),
                    (ax[1], "grid_max_norm_M0i", "whole grid: max |M_0i| after heal")):
        if k in s:
            a.axhline(s[k], color="gray", ls=":", label="seed")
        a.set_title(t, fontsize=10)
        a.set_xlabel("d")
        a.grid(alpha=0.3)
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r3_pair_dressing.png"), dpi=130)
    plt.close(fig)
    # descent traces
    rows = {r["tag"]: r for r in J.get("rows", [])}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax, lam in zip(axes, LAMBDAS):
        for tag, r in rows.items():
            if r.get("lam") != lam or r.get("n") != N_MAIN or "descent" not in r:
                continue
            tr = r["descent"]["trace"]
            if not tr:
                continue
            ax.plot([t["acc"] for t in tr], [t["E"] - r["descent"]["E0"] for t in tr],
                    "-" if r["dressed"] else ":", lw=1,
                    label=f"{r['kind']} d{r['d']:g} {'dr' if r['dressed'] else 'un'}")
        ax.set_title(f"lambda = {lam:g}: E - E0 vs accepted steps", fontsize=10)
        ax.set_xlabel("accepted steps")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=5, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r3_pair_traces.png"), dpi=130)
    plt.close(fig)
    # poisson
    pois = res["controls"]["e_poisson_reference"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for key, mk in (("n32_L48_sigma_rms", "o-"), ("n32_L48_sigma_2h", "o--"),
                    ("n48_L72_sigma_rms", "s-"), ("n48_L72_sigma_2h", "s--")):
        p = pois[key]
        ax[0].plot([r["d"] for r in p["rows"]], [r["E_box"] for r in p["rows"]], mk,
                   label=f"{key} box (sigma {p['sigma']:.2f})")
        ax[0].plot([r["d"] for r in p["rows"]], [r["E_free"] for r in p["rows"]], "x:", color="gray")
        ax[1].plot([r["d"] for r in p["rows"]], [r["image_ratio_F"] for r in p["rows"]], mk, label=key)
    ax[0].set_title("Poisson reference E_ref(d) (box vs free x)", fontsize=10)
    ax[1].set_title("image correction F_box / F_free", fontsize=10)
    for a in ax:
        a.set_xlabel("d")
        a.grid(alpha=0.3)
        a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r3_pair_poisson.png"), dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=("gate", "seed", "relax", "collect"))
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--ladder", action="store_true")
    ap.add_argument("--only-ladder", action="store_true")
    a = ap.parse_args()
    if a.stage == "gate":
        out = stage_gate()
        sys.exit(0 if out["gate_pass"] else 1)
    elif a.stage == "seed":
        stage_seed()
    elif a.stage == "relax":
        stage_relax(a.workers, a.ladder or a.only_ladder, a.only_ladder)
    else:
        stage_collect()


if __name__ == "__main__":
    main()
