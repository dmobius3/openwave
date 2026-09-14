"""M5.32 R21-1 to R21-4: the three-axis hedgehog under three boundaries at
equal residual, the (g, delta) scan of the frame saddle, and the W1 ladder
(ledger section 6.11, amended on the author's 2026-09-14 14:45 UTC reply).

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1), N = M eta. Code
branch s = -1: M_vac = diag(g, 1, delta, 0), the N-spectrum (-g, 1, delta, 0);
the main rows at g = 8, delta = 0.3. The static energy is the R20 instrument's
(m5_32_r20_1_axes.energy_grad, consumed read-only):
    E[M] = 4 h^3 sum_br wt sum_cells I1(A) + V[M],  A_i = d_i M (sym stencil)
    V = V4 = w sum_p (tr N^p - C_p)^2, C_p = sum_i q_i^p, w = W1 (R21-1, R21-3)
                                                  or W1 x 5, W1 x 25 (R21-4)
The objects: S_1 (1, delta, 0), S_d (delta, 1, 0), S_0 (0, 1, delta), the R20
seeds (m5_32_r20_1_axes.seed_axes). The three boundaries, n32 L48 (h 1.5):
    B_seed  the R20 boundary: the Dirichlet shell (B3.pin_shell, depth 1.6) at
            the seed's values, the seed's line defect frozen at half-width 24
    B_far   the Dirichlet shell at the relaxed far field: the R20 n64 L96 end
            field restricted to its central 32^3 cells (the same h, the cell
            centers coincide exactly: (i - 31.5) h = (i + 16 - 15.5 - 32) h),
            the shell then holds the values a bigger box chose at half-width 24
    B_free  no pin: every cell relaxes; the stack's derivative B3.d1 is
            one-sided at the faces (forward / backward branches), so the
            discrete functional carries its natural boundary and the exact
            adjoint gradient B3.d1_adj is unchanged (R21-0 (a) gates the face
            cells' gradient, which the pin had zeroed since R3)
The descent: R3's FIRE (the R20 copy with the pin optional), the R20 kill
rules, a resumable checkpoint every CKPT_EVERY accepted steps; then the
R21-2 polish inside the same job: scipy L-BFGS (memory 20) on the free
entries (the 10 upper-triangle entries per free cell; the gradient of the
symmetric parametrization is G_aa on the diagonal and 2 G_ab off it), the
gate POLISHED iff max |G| on the free cells < 1e-3 (a decade under R20's
best FIRE end) within POLISH_ITERS iterations, else the residual is
recorded and the row is labeled FALLING at that max |G|. Every comparison
across boundaries uses the polished fields; the FIRE-end field is kept and
read too (the equal-budget row, R20's bar).
Reads per row (the R20 reads on the polished field: E and its parts, the
Mermin degree of each eigenvector on the r 6 / 9 / 12 cubes, the shells'
eigenvalues, the polar-axis gaps and biaxiality, the tube read and the
string tension, the core ball, the Derrick reads with R_* and the virial,
the frame read at s 0.02 / 0.05, E(< R) on the shells and the pinned-shell
energy), plus the R21 reads:
    solid-angle degree  the exact degree of each oriented eigenvector map on
                        the r 6 / 9 / 12 cubes (the R20 audit's reader,
                        consumed read-only; the audit corrected the Mermin
                        flux on S_d's near-degenerate axis)
    frame3              d2E / ds2 along the R3 boost dressing at s 0.01 /
                        0.02 / 0.05 with the even extrapolation
                        c(s) = c0 + c1 s^2 to s -> 0 (R21-3's number)
    tail slope          the slope of the mean cell density against r on the
                        shells r 8 to 16 (Mikulski's read; its sign)
    face layer          the energy on the cells the pin used to hold (the
                        outer 1.6 of the box), pinned or not
    norm3               tr(M_3x3^2) pointwise against the vacuum's
                        1 + delta^2 (1.09): the profile along the former
                        polar axis with the eigenvalue gaps there, the
                        r 6 / 9 / 12 shells, the fractions within 10 percent
                        of and below 70 percent of the vacuum inside r 12
                        (the escape-route read, calibrated in R21-0 (g))
Outcomes at collect (pre-registered, ledger section 6.11):
    per axis    INTERIOR_SET iff E(< 12) agrees to 2 percent across the three
                boundaries and the solid-angle degree of the winding
                eigenvector stays 1 (within 0.1) on the r 9 cube under all
                three; BOUNDARY_SET iff the polished total under B_free or
                B_far differs from B_seed's by more than 5 percent of the
                total, or the winding is lost under any boundary; else
                BOUNDARY_UNDECIDED (with the residuals)
    the string  STRING_ESCAPES iff under B_free the tube read in the flat
                middle falls below S_1's 0.02 or the degree is lost;
                STRING_HELD iff the flat-middle tension stays within
                30 percent of B_seed's 0.12 with the degree kept; the route
                of any lost winding: ESCAPE_BY_EXCHANGE iff the axis norm
                stays within 10 percent of the vacuum's where the winding
                was lost and the eigenvalues cross there (gap < 0.02),
                ESCAPE_BY_MELTING iff that norm falls by more than 30 percent
    the box     per axis under B_free: AXIS_BOX_LIMITED iff dE / d lambda < 0
                and R_* > L / 4; AXIS_LOCALIZED iff R_* < L / 6 with the
                virial within 30 percent of 3
    the triple  THREE_AXES_DISTINCT iff the three polished B_free energies
                differ pairwise by more than three times the largest
                boundary difference of any axis; else AXES_DEGENERATE;
                Koide Q a read
    the W1 ladder  SCALE_IS_THE_POTENTIALS iff at W1 x 25 the polished S_1 row
                has R_* < L / 6, the virial within 30 percent of 3 and
                dE / d lambda within 10 percent of zero relative to R20's;
                SCALE_IS_THE_BOXS iff R_* stays above L / 4 or the virial
                above 6; else SCALE_UNDECIDED
    the frame   the fits c (g - g_0)^2 on the six g rows (g 1, 2, 4, 8, 16,
                32 at delta 0.3; g 8 and g 32 the polished R20 rows) and
                c' (delta - delta_0)^2 on the three delta rows (0.3, 0.6,
                1.0 at g 8), on the s -> 0 curvatures, least squares in
                log |c|; FRAME_TWO_ROOTS iff g_0 within 0.3 of 1 and delta_0
                within 0.15 of 1 and the product c (g - g_0)^2
                (delta - delta_0)^2 fits all eight rows within 20 percent
                (relative to each row, floored at 20 percent of the median
                |c| so that a row at a root cannot fail by 0 / 0);
                FRAME_ONE_ROOT iff exactly one root sits there (g_0 within
                0.3 of 0 is the pure g^2 law); FRAME_NOT_PRODUCT iff the
                product fit misses 20 percent on any row

STAGES (python3 m5_32_r21_1_runs.py STAGE [--workers W] [--mem-budget GB] [--json-suffix S]):
    smoke    the gates: the free-boundary gradient at face / edge cells, the
             n64 restriction, the polish parametrization, one tiny job
    main     the 15 n32 rows (R21-1 six, R21-4 three plus the optional B_free
             copy of the W1 x 25 S_1 row, R21-3 six) and the four polish-only
             R20 rows (S_1, S_d, S_0 under B_seed from their 9000-step end
             fields; S_1 at g 32), longest first
    n64      the optional S_d B_free row at n64 from the R20 n64 end field,
             2000 accepted steps, a 200-iteration polish
    collect  tables, outcomes, plots (merges every side-pool JSON)
Out: ../data/m5_32_r21_1_runs.json (partials after every job),
     ../data/m5_32_r21_1/*.npz (local), ../plots/m5_32_r21_1_*.png
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
from scipy.optimize import minimize  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
OUT_JSON = os.path.join(DATA, "m5_32_r21_1_runs.json")
MAIN_JSON = OUT_JSON
OUT_NPZ = os.path.join(DATA, "m5_32_r21_1")
R20_NPZ = os.path.join(DATA, "m5_32_r20_1")
R20_JSON = os.path.join(DATA, "m5_32_r20_1_axes.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    argv = sys.argv
    sys.argv = [argv[0]]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


R20 = _load("m5_32_r20_1_axes", "m5_32_r20_1_axes.py")
R20A = _load(
    "m5_32_r20_1_audit", "m5_32_r20_1_audit.py"
)  # orient_tree + solid_angle_degree, read-only
R0, EN, B3, LAG, PAIR, R3 = R20.R0, R20.EN, R20.B3, R20.LAG, R20.PAIR, R20.R3
W1 = B3.W1
G_MAIN, DELTA_MAIN = 8.0, 0.3
STEPS_LONG, STEPS_SCAN = 4500, 3000
POLISH_ITERS, POLISH_GATE, POLISH_MEM = 3000, 1e-3, 20
POLISH_ITERS_N64 = 200
CKPT_EVERY = 500
RUNAWAY_FACTOR, DIVE_FLOOR = 3.0, -1e6
FRAME_S3 = (0.01, 0.02, 0.05)
R_CORE, R_TUBE = R20.R_CORE, R20.R_TUBE
MEM_GB = R20.MEM_GB
IU = np.triu_indices(4)
OFF = np.where(IU[0] != IU[1], 2.0, 1.0)
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def lam_of(obj, delta):
    return {"S1": (1.0, delta, 0.0), "Sd": (delta, 1.0, 0.0), "S0": (0.0, 1.0, delta)}[obj]


WINDING_RANK = {"S1": 2, "Sd": 1, "S0": 0}


def cfg_of(n, L, g, delta):
    return B3.base_cfg(s=-1.0, g=g, n=n, L=float(L), delta=delta)


def params_of(g, delta):
    return LAG.default_params(s=-1.0, g=g, delta=delta)


def pot_of(cfg, w1s):
    """the stack's V4 (None) or the same targets at the scaled weight (the R20 'v4' hook)."""
    return None if w1s == 1.0 else ("v4", R0.roots_of(cfg), W1 * w1s)


def job_tag(j):
    return f"{j['obj']}_B{j['bnd']}_g{j['g']:g}_d{j['delta']:g}_w{j['w1s']:g}_n{j['n']}"


def restrict_center(M_big, n_small):
    """the central n_small^3 cells of an n_big^3 field at the same h (the cell centers coincide exactly)."""
    nb = M_big.shape[0]
    c = (nb - n_small) // 2
    return M_big[c : c + n_small, c : c + n_small, c : c + n_small].copy()


def seed_of(j, cfg):
    sf = j.get("seed_from")
    if sf == "n64restrict":
        M64 = np.load(os.path.join(R20_NPZ, f"{j['obj']}_v4std_n64_g8.npz"))["M"]
        return (
            restrict_center(M64, cfg["n"]),
            f"R20 {j['obj']}_v4std_n64_g8 (4500 steps, FALLING) restricted to the central 32^3",
        )
    if sf == "r20_x4500":
        return (
            np.load(os.path.join(R20_NPZ, f"{j['obj']}_v4std_n32_g8_x4500.npz"))["M"],
            f"R20 {j['obj']}_v4std_n32_g8_x4500 (9000 steps, FALLING)",
        )
    if sf == "r20_g32":
        return (
            np.load(os.path.join(R20_NPZ, "S1_v4std_n32_g32.npz"))["M"],
            "R20 S1_v4std_n32_g32 (4500 steps, CONVERGED)",
        )
    if sf == "r20_n64":
        return (
            np.load(os.path.join(R20_NPZ, f"{j['obj']}_v4std_n64_g8.npz"))["M"],
            f"R20 {j['obj']}_v4std_n64_g8 (4500 steps, FALLING)",
        )
    return (
        R20.seed_axes(cfg, lam_of(j["obj"], j["delta"])),
        f"seed_axes {lam_of(j['obj'], j['delta'])}",
    )


# ================= the descent (the R20 FIRE copied with the pin optional) =================
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


def free_mask(cfg, pinned):
    n, h = cfg["n"], cfg["h"]
    return (~B3.pin_shell(n, h)) if pinned else np.ones((n, n, n), dtype=bool)


def descend(
    M0,
    cfg,
    p,
    pot,
    steps_acc,
    it_cap,
    tag,
    pinned,
    log_every=100,
    dt0=0.02,
    dt_max=0.2,
    ckpt_path=None,
):
    free = free_mask(cfg, pinned)[..., None, None].astype(float)
    M = M0.copy()
    E0, G, info = R20.energy_grad(M, cfg, p, pot)
    m0i_seed = float(np.max(np.abs(M0[..., 0, 1:])))
    out = {
        "E0": float(E0),
        "steps_acc_budget": steps_acc,
        "it_cap": it_cap,
        "pin": (
            "B3.pin_shell depth 1.6 (Dirichlet at the seed values)"
            if pinned
            else "none (B_free: every cell relaxes, the faces one-sided)"
        ),
        "fire": {"dt0": dt0, "dt_max": dt_max, "alpha0": 0.1, "dt_min": 1e-7},
        "max_abs_M0i_seed": m0i_seed,
        "trace": [],
        "resumed_from_ckpt": None,
    }
    if not np.isfinite(E0) or G is None:
        out.update(
            stop="DIVERGED (seed energy undefined)",
            verdict="DIVERGED",
            steps_run=0,
            accepted=0,
            E_end=float("nan"),
        )
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
        E, G, info = R20.energy_grad(M, cfg, p, pot)
        F = -G * free
        it, n_acc, n_rej, n_rej_locus = st["it"], st["acc"], st["rej"], st["rej_locus"]
        dt, alpha, n_up, E_prev, fmax = (
            st["dt"],
            st["alpha"],
            st["n_up"],
            st["E_prev"],
            float(np.max(np.abs(F))),
        )
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
        E, G, info = R20.energy_grad(M_try, cfg, p, pot)
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
                stop = (
                    "LOCUS-HIT"
                    if locus_loss
                    else "STALLED (dt collapsed, no descent direction accepted)"
                )
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
            row = {
                "it": it,
                "acc": n_acc,
                "E": float(E),
                "fmax": fmax,
                "dt": dt,
                "min_gap": info["min_gap"],
            }
            out["trace"].append(row)
            log(
                f"{tag} it {it:5d} acc {n_acc:5d} E {E:14.6f} fmax {fmax:.3e} dt {dt:.2e} rej {n_rej}"
            )
        if ckpt_path and n_acc % CKPT_EVERY == 0 and n_acc < steps_acc:
            _ckpt_save(
                ckpt_path,
                M,
                v,
                {
                    "it": it,
                    "acc": n_acc,
                    "rej": n_rej,
                    "rej_locus": n_rej_locus,
                    "dt": dt,
                    "alpha": alpha,
                    "n_up": n_up,
                    "E_prev": float(E_prev),
                    "trace": out["trace"],
                },
            )
    if stop == "budget" and n_acc < steps_acc:
        stop = f"IT_CAP ({it_cap} iterations before {steps_acc} accepted)"
    out.update(
        {
            "stop": stop,
            "steps_run": it,
            "accepted": n_acc,
            "rejected": n_rej,
            "rejected_locus": n_rej_locus,
            "dt_final": dt,
            "fmax_end": fmax,
            "E_end": float(E_prev),
            "E_drop": float(E0 - E_prev),
            "runaway": runaway,
        }
    )
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
            out["gate"] = {
                "last_third_rel_lt_1e-3": bool(out["last_third_rel"] < 1e-3),
                "fmax_fell_2_decades": bool(out["fmax_decades"] >= 2.0),
            }
        else:
            out["verdict"] = "budget (trace too short)"
    else:
        out["verdict"] = stop
    return M, out


# ================= the R21-2 polish: L-BFGS on the free entries =================
class _Polished(Exception):
    pass


def pack(M, mask):
    return M[mask][:, IU[0], IU[1]].ravel().copy()


def unpack(x, M_base, mask):
    S = x.reshape(-1, 10)
    F = np.zeros((S.shape[0], 4, 4))
    F[:, IU[0], IU[1]] = S
    F[:, IU[1], IU[0]] = S
    Mx = M_base.copy()
    Mx[mask] = F
    return Mx


def polish(M, cfg, p, pot, pinned, tag, max_iter=POLISH_ITERS, gate=POLISH_GATE, log_every=50):
    """scipy L-BFGS (memory POLISH_MEM) on the 10 free entries per free cell; the gate on the stack's own
    max |G| over the free cells (the FIRE residual measure), read from the accepted iterates through a cache
    of the evaluations (no extra energy call). Returns the polished field and the record."""
    mask = free_mask(cfg, pinned)
    x0 = pack(M, mask)
    cache = {}
    st = {"n_eval": 0, "it": 0, "trace": [], "x_stop": None, "t0": time.time()}

    def fun(x):
        Mx = unpack(x, M, mask)
        E, G, info = R20.energy_grad(Mx, cfg, p, pot)
        st["n_eval"] += 1
        if G is None or not np.isfinite(E) or not info["ok"]:
            return 1e30, np.zeros_like(x)
        Gm = G[mask]
        fmax = float(np.max(np.abs(Gm)))
        key = hash(x.tobytes())
        cache[key] = (float(E), fmax, float(info["max_abs_M0i"]), float(info["min_gap"]))
        if len(cache) > 64:
            for k in list(cache)[:32]:
                cache.pop(k, None)
        return float(E), (Gm[:, IU[0], IU[1]] * OFF).ravel()

    def cb(intermediate_result):
        res = intermediate_result  # scipy passes the OptimizeResult when the parameter carries this name
        st["it"] += 1
        rec = cache.get(hash(np.asarray(res.x).tobytes()))
        if rec is None:
            return
        E, fmax, m0i, gap = rec
        if st["it"] % log_every == 0 or st["it"] == 1:
            st["trace"].append(
                {"it": st["it"], "n_eval": st["n_eval"], "E": E, "fmax": fmax, "min_gap": gap}
            )
            log(
                f"{tag} polish it {st['it']:4d} eval {st['n_eval']:4d} E {E:14.6f} fmax {fmax:.3e}"
            )
        if fmax < gate:
            st["x_stop"] = np.asarray(res.x).copy()
            st["trace"].append(
                {"it": st["it"], "n_eval": st["n_eval"], "E": E, "fmax": fmax, "min_gap": gap}
            )
            raise _Polished()

    E_in, G_in, _ = R20.energy_grad(M, cfg, p, pot)
    fmax_in = float(np.max(np.abs(G_in[mask])))
    out = {
        "method": f"scipy L-BFGS-B, maxcor {POLISH_MEM}, on the free entries",
        "free_cells": int(mask.sum()),
        "pinned": pinned,
        "max_iter": max_iter,
        "gate_max_abs_G": gate,
        "E_in": float(E_in),
        "fmax_in": fmax_in,
        "trace": [],
    }
    if fmax_in < gate:
        out.update(
            E_end=float(E_in),
            fmax_end=fmax_in,
            iters=0,
            n_eval=1,
            verdict="POLISHED",
            note="already under the gate",
        )
        return M.copy(), out
    x_end = None
    try:
        res = minimize(
            fun,
            x0,
            jac=True,
            method="L-BFGS-B",
            callback=cb,
            options={
                "maxcor": POLISH_MEM,
                "maxiter": max_iter,
                "maxfun": 3 * max_iter,
                "gtol": 1e-14,
                "ftol": 1e-16,
            },
        )
        x_end = np.asarray(res.x)
        out["scipy_message"] = str(res.message)
    except _Polished:
        x_end = st["x_stop"]
        out["scipy_message"] = "stopped at the gate by the callback"
    Mp = unpack(x_end, M, mask)
    E_end, G_end, info = R20.energy_grad(Mp, cfg, p, pot)
    fmax_end = float(np.max(np.abs(G_end[mask])))
    out.update(
        {
            "trace": st["trace"],
            "iters": st["it"],
            "n_eval": st["n_eval"],
            "E_end": float(E_end),
            "fmax_end": fmax_end,
            "E_drop": float(E_in - E_end),
            "max_abs_M0i_end": float(info["max_abs_M0i"]),
            "min_gap_end": float(info["min_gap"]),
            "wall_s": round(time.time() - st["t0"], 1),
        }
    )
    out["verdict"] = "POLISHED" if fmax_end < gate else f"FALLING (at max |G| = {fmax_end:.2e})"
    log(
        f"{tag} polish {out['verdict']}: E {E_in:.6f} -> {E_end:.6f}, fmax {fmax_in:.2e} -> {fmax_end:.2e}, {st['it']} its, {st['n_eval']} evals, {out['wall_s']} s"
    )
    return Mp, out


# ================= the R21 reads =================
def frame_reads3(M, cfg, p, pot):
    """d2E / ds2 along the R3 boost dressing at s in FRAME_S3 and the even extrapolation c(s) = c0 + c1 s^2 to s -> 0."""
    E0 = R20.energy_parts(M, cfg, p, pot)
    out = {"E_0": E0["E_total"], "s": list(FRAME_S3), "d2E_ds2": {}, "dV": {}}
    for sv in FRAME_S3:
        Ep = R20.energy_parts(R3.conj(R3.boost_at(cfg, 0.0, sv)[0], M), cfg, p, pot)
        Em = R20.energy_parts(R3.conj(R3.boost_at(cfg, 0.0, -sv)[0], M), cfg, p, pot)
        out["d2E_ds2"][f"{sv:g}"] = float(
            (Ep["E_total"] + Em["E_total"] - 2.0 * E0["E_total"]) / sv**2
        )
        out["dV"][f"{sv:g}"] = float(max(abs(Ep["V"] - E0["V"]), abs(Em["V"] - E0["V"])))
    s = np.array(FRAME_S3)
    c = np.array([out["d2E_ds2"][f"{sv:g}"] for sv in FRAME_S3])
    A = np.stack([np.ones_like(s), s**2], axis=1)
    coef, *_ = np.linalg.lstsq(A, c, rcond=None)
    out["c0_extrapolated"] = float(coef[0])
    out["c1_s2_coefficient"] = float(coef[1])
    out["c0_two_point_0.01_0.02"] = float((4.0 * c[0] - c[1]) / 3.0)
    out["fit_max_abs_resid"] = float(np.max(np.abs(A @ coef - c)))
    out["saddle_along_boost_dressing"] = bool(coef[0] < 0.0)
    return out


def _outer_plane(n):
    P = np.zeros((n, n, n), dtype=bool)
    for ax in range(3):
        sl = [slice(None)] * 3
        sl[ax] = 0
        P[tuple(sl)] = True
        sl[ax] = n - 1
        P[tuple(sl)] = True
    return P


def extra_reads(M, cfg, pot, pinned):
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rho = np.sqrt(X * X + Y * Y)
    M3 = M[..., 1:, 1:]
    lam, vec = np.linalg.eigh(M3)
    e = R20.density(M, cfg, pot)
    pin = B3.pin_shell(n, h)
    out = {}
    out["face_layer"] = {
        "E": float(np.sum(e[pin])),
        "cells": int(pin.sum()),
        "frac_of_total": float(np.sum(e[pin]) / max(np.sum(e), 1e-300)),
        "pinned": pinned,
        "E_outer_plane": float(np.sum(e[_outer_plane(n)])),
    }
    # the tail slope on the shells r 8 to 16: the mean cell density per shell of width h, and the shell energy
    edges = np.arange(8.0, 16.0 + 0.5 * h, h)
    rc, dens, esh = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (r >= a) & (r < b)
        if m.sum() > 0:
            rc.append(0.5 * (a + b))
            dens.append(float(np.mean(e[m])))
            esh.append(float(np.sum(e[m])))
    rc, dens, esh = np.array(rc), np.array(dens), np.array(esh)
    out["tail"] = {
        "r": rc.tolist(),
        "mean_cell_density": dens.tolist(),
        "shell_energy": esh.tolist(),
        "slope_density": float(np.polyfit(rc, dens, 1)[0]),
        "slope_shell_energy": float(np.polyfit(rc, esh, 1)[0]),
        "slope_log_density": float(np.polyfit(rc, np.log(np.maximum(dens, 1e-300)), 1)[0]),
    }
    # the spatial-block norm against the vacuum's
    nrm = np.einsum("...ij,...ji->...", M3, M3)
    vac = 1.0 + cfg["delta"] ** 2
    ax = (rho < 0.75 * h) & (np.abs(Z) > R_CORE) & (np.abs(Z) < 0.5 * L - 2.0)
    zs = np.unique(np.round(Z[ax], 6))
    prof = []
    for zv in zs:
        m = ax & (np.abs(Z - zv) < 1e-6)
        prof.append(
            {
                "z": float(zv),
                "norm": float(np.mean(nrm[m])),
                "gap01": float(np.mean(lam[m][:, 1] - lam[m][:, 0])),
                "gap12": float(np.mean(lam[m][:, 2] - lam[m][:, 1])),
                "lam": [float(v) for v in np.mean(lam[m], axis=0)],
            }
        )
    gmin = np.minimum(lam[ax][:, 1] - lam[ax][:, 0], lam[ax][:, 2] - lam[ax][:, 1])
    ball = r < 12.0
    out["norm3"] = {
        "vacuum": vac,
        "axis_profile": prof,
        "axis_norm_min": float(np.min(nrm[ax])),
        "axis_norm_mean": float(np.mean(nrm[ax])),
        "axis_norm_min_rel": float(np.min(nrm[ax]) / vac),
        "axis_cells": int(ax.sum()),
        "axis_crossings_gap_lt_0.02": int(np.sum(gmin < 0.02)),
        "axis_gap_min": float(np.min(gmin)),
        "shells": {
            f"{hv:g}": {
                "mean": float(np.mean(nrm[(r > hv - 1) & (r < hv + 1)])),
                "min": float(np.min(nrm[(r > hv - 1) & (r < hv + 1)])),
                "max": float(np.max(nrm[(r > hv - 1) & (r < hv + 1)])),
            }
            for hv in (6.0, 9.0, 12.0)
        },
        "frac_within_10pct_r_lt_12": float(np.mean(np.abs(nrm[ball] - vac) <= 0.1 * vac)),
        "frac_below_70pct_r_lt_12": float(np.mean(nrm[ball] < 0.7 * vac)),
        "core_norm_center": float(nrm[r == r.min()][0]),
    }
    # the solid-angle degree of each oriented eigenvector map (the R20 audit's exact reader)
    c0 = n // 2
    sa = {}
    for k in range(3):
        vo, ncf = R20A.orient_tree(vec[..., :, k])
        sa[f"rank{k}"] = {
            "conflicts": int(ncf),
            "degree": {
                f"{hv:g}": float(R20A.solid_angle_degree(vo, c0, int(round(hv / h))))
                for hv in (6.0, 9.0, 12.0)
                if c0 + int(round(hv / h)) < n
            },
        }
    out["solid_angle"] = sa
    return out


def reads_all(M, cfg, p, pot, pinned, frame3=True):
    out = R20.reads(M, cfg, p, pot)
    out.update(extra_reads(M, cfg, pot, pinned))
    if frame3:
        out["frame3"] = frame_reads3(M, cfg, p, pot)
    return out


# ================= jobs =================
def jobs_main():
    J = []
    for obj in ("S1", "Sd", "S0"):
        J.append(
            dict(
                obj=obj,
                bnd="far",
                g=G_MAIN,
                delta=DELTA_MAIN,
                w1s=1.0,
                n=32,
                L=48.0,
                steps=STEPS_LONG,
                seed_from="n64restrict",
                stage="R21-1",
            )
        )
        J.append(
            dict(
                obj=obj,
                bnd="free",
                g=G_MAIN,
                delta=DELTA_MAIN,
                w1s=1.0,
                n=32,
                L=48.0,
                steps=STEPS_LONG,
                seed_from="r20_x4500",
                stage="R21-1",
            )
        )
    J.append(
        dict(
            obj="S1",
            bnd="seed",
            g=G_MAIN,
            delta=DELTA_MAIN,
            w1s=5.0,
            n=32,
            L=48.0,
            steps=STEPS_LONG,
            stage="R21-4",
        )
    )
    J.append(
        dict(
            obj="S1",
            bnd="seed",
            g=G_MAIN,
            delta=DELTA_MAIN,
            w1s=25.0,
            n=32,
            L=48.0,
            steps=STEPS_LONG,
            stage="R21-4",
        )
    )
    J.append(
        dict(
            obj="Sd",
            bnd="seed",
            g=G_MAIN,
            delta=DELTA_MAIN,
            w1s=25.0,
            n=32,
            L=48.0,
            steps=STEPS_LONG,
            stage="R21-4",
        )
    )
    J.append(
        dict(
            obj="S1",
            bnd="free",
            g=G_MAIN,
            delta=DELTA_MAIN,
            w1s=25.0,
            n=32,
            L=48.0,
            steps=STEPS_LONG,
            stage="R21-4 (optional)",
        )
    )
    for g in (1.0, 2.0, 4.0, 16.0):
        J.append(
            dict(
                obj="S1",
                bnd="seed",
                g=g,
                delta=DELTA_MAIN,
                w1s=1.0,
                n=32,
                L=48.0,
                steps=STEPS_SCAN,
                stage="R21-3",
            )
        )
    for d in (0.6, 1.0):
        J.append(
            dict(
                obj="S1",
                bnd="seed",
                g=G_MAIN,
                delta=d,
                w1s=1.0,
                n=32,
                L=48.0,
                steps=STEPS_SCAN,
                stage="R21-3",
            )
        )
    for obj in ("S1", "Sd", "S0"):
        J.append(
            dict(
                obj=obj,
                bnd="seed",
                g=G_MAIN,
                delta=DELTA_MAIN,
                w1s=1.0,
                n=32,
                L=48.0,
                steps=0,
                seed_from="r20_x4500",
                stage="R21-2 (the R20 rows)",
            )
        )
    J.append(
        dict(
            obj="S1",
            bnd="seed",
            g=32.0,
            delta=DELTA_MAIN,
            w1s=1.0,
            n=32,
            L=48.0,
            steps=0,
            seed_from="r20_g32",
            stage="R21-2 (the R20 g 32 row)",
        )
    )
    J.sort(key=lambda j: -j["steps"])
    return J


def jobs_n64():
    return [
        dict(
            obj="Sd",
            bnd="free",
            g=G_MAIN,
            delta=DELTA_MAIN,
            w1s=1.0,
            n=64,
            L=96.0,
            steps=2000,
            seed_from="r20_n64",
            stage="R21-1 (optional n64)",
            polish_iters=POLISH_ITERS_N64,
        )
    ]


def run_job(j):
    t0 = time.time()
    cfg = cfg_of(j["n"], j["L"], j["g"], j["delta"])
    p = params_of(j["g"], j["delta"])
    pot = pot_of(cfg, j["w1s"])
    pinned = j["bnd"] != "free"
    tag = job_tag(j)
    row = dict(j)
    row.update(
        {
            "tag": tag,
            "h": cfg["h"],
            "lam": lam_of(j["obj"], j["delta"]),
            "roots": list(R0.roots_of(cfg)),
            "W1_eff": W1 * j["w1s"],
            "pinned": pinned,
        }
    )
    try:
        M0, seed_note = seed_of(j, cfg)
        row["seed_note"] = seed_note
        row["seed_reads"] = reads_all(M0, cfg, p, pot, pinned, frame3=False)
        os.makedirs(OUT_NPZ, exist_ok=True)
        if j["steps"] > 0:
            M, des = descend(
                M0,
                cfg,
                p,
                pot,
                j["steps"],
                2 * j["steps"],
                tag,
                pinned,
                ckpt_path=os.path.join(OUT_NPZ, f"{tag}_ckpt"),
            )
            row["descent"] = des
            finite = bool(np.all(np.isfinite(M))) and np.isfinite(des["E_end"])
            row["status"] = (
                "OK"
                if finite
                and (des["stop"].startswith("budget") or des["stop"].startswith("IT_CAP"))
                else des["stop"].split(" ")[0]
            )
            if finite:
                row["fire_reads"] = reads_all(M, cfg, p, pot, pinned, frame3=False)
                row["E_fire"] = row["fire_reads"]["energy"]["E_total"]
                np.savez_compressed(
                    os.path.join(OUT_NPZ, f"{tag}_fire.npz"), M=M.astype(np.float64)
                )
        else:
            M, finite = M0, True
            row["status"] = "OK"
            row["descent"] = {"verdict": "none (polish-only row)", "accepted": 0}
        if finite:
            Mp, pol = polish(
                M, cfg, p, pot, pinned, tag, max_iter=j.get("polish_iters", POLISH_ITERS)
            )
            row["polish"] = pol
            row["end_reads"] = reads_all(Mp, cfg, p, pot, pinned, frame3=(j["n"] == 32))
            row["E"] = row["end_reads"]["energy"]["E_total"]
            row["label"] = pol["verdict"]
            np.savez_compressed(os.path.join(OUT_NPZ, f"{tag}.npz"), M=Mp.astype(np.float64))
        else:
            row["E"] = None
    except Exception as e:  # noqa: BLE001
        import traceback

        row["status"] = "DIVERGED"
        row["stop"] = f"exception: {e!r}"
        row["traceback"] = traceback.format_exc()
        row["E"] = None
    row["wall_s"] = round(time.time() - t0, 1)
    log(
        f"DONE {tag} status {row['status']} label {row.get('label')} E {row.get('E')} wall {row['wall_s']}"
    )
    return row


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R21-1 to R21-4", "rows": {}}


def load_json_all():
    if os.path.exists(MAIN_JSON):
        with open(MAIN_JSON) as f:
            J = json.load(f)
    else:
        J = {"task": "M5.32 R21-1 to R21-4", "rows": {}}
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


def run_pool(jobs, workers, label, mem_budget=12.0):
    rows = load_json_all()["rows"]
    todo = [j for j in jobs if job_tag(j) not in rows]
    log(
        f"{label}: {len(todo)} jobs, {workers} workers, memory budget {mem_budget} GB (done already: {len(jobs) - len(todo)})"
    )
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
                    log(
                        f"submit {job_tag(j)} (running {len(running)}, memory {used:.1f} GB, pending {len(pending)})"
                    )
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


# ================= smoke =================
def smoke():
    out = {}
    # (a) the free-boundary gradient at face and edge cells against a 4-point stencil (n8, the cells the pin zeroed)
    cfg = cfg_of(8, 12.0, G_MAIN, DELTA_MAIN)
    p = params_of(G_MAIN, DELTA_MAIN)
    rng = np.random.default_rng(3)
    M = R20.seed_axes(cfg, lam_of("Sd", DELTA_MAIN))
    M = M + 0.05 * B3.sym4(rng.standard_normal(M.shape))
    M[..., 0, 1:] = 0.0
    M[..., 1:, 0] = 0.0
    E, G, info = R20.energy_grad(M, cfg, p, None)
    worst = {"face": 0.0, "edge": 0.0, "corner": 0.0, "interior": 0.0}
    cells = {"face": (0, 3, 4), "edge": (0, 0, 5), "corner": (7, 7, 7), "interior": (3, 4, 4)}
    for kind, (i, j, k) in cells.items():
        for a, b in ((1, 1), (1, 2), (2, 3), (0, 0)):
            D = np.zeros_like(M)
            D[i, j, k, a, b] = 1.0
            D[i, j, k, b, a] = 1.0
            dd = float(np.sum(G * D))
            e = 1e-4
            f = [
                R20.energy_grad(M + s * e * D, cfg, p, None, need_grad=False)[0]
                for s in (-2, -1, 1, 2)
            ]
            fd = (f[0] - 8 * f[1] + 8 * f[2] - f[3]) / (12 * e)
            worst[kind] = max(worst[kind], abs(fd - dd) / max(abs(dd), 1e-6))
    out["a_free_gradient_stencil4_rel"] = worst
    out["a_pin_zeroes_faces"] = {
        "pinned_cells": int(B3.pin_shell(8, 1.5).sum()),
        "max_abs_G_on_pinned": float(np.max(np.abs(G[B3.pin_shell(8, 1.5)]))),
    }
    # (b) the restriction: coordinates coincide, values bitwise
    X64 = B3.coords(64, 1.5)[0][:, 0, 0]
    X32 = B3.coords(32, 1.5)[0][:, 0, 0]
    out["b_coord_coincidence_max_abs"] = float(np.max(np.abs(X64[16:48] - X32)))
    M64 = np.load(os.path.join(R20_NPZ, "Sd_v4std_n64_g8.npz"))["M"]
    Mr = restrict_center(M64, 32)
    out["b_restriction_bitwise"] = bool(np.array_equal(Mr, M64[16:48, 16:48, 16:48]))
    cfg32 = cfg_of(32, 48.0, G_MAIN, DELTA_MAIN)
    ep = R20.energy_parts(Mr, cfg32, p, None)
    rad = R20.radial_profile(Mr, cfg32, None)
    out["b_restricted_Sd_energy_n32"] = {
        "E_total": ep["E_total"],
        "E_pin_shell": rad["E_pin_shell"],
    }
    # the polish parametrization: the packed gradient against a stencil on the packed function
    out["M_symmetry_max_abs"] = float(np.max(np.abs(M - np.swapaxes(M, -1, -2))))
    for pinned_k, lab in ((True, "interior_only"), (False, "all_cells_faces_included")):
        mask = free_mask(cfg, pinned_k)
        x0 = pack(M, mask)

        def fun(x):
            Mx = unpack(x, M, mask)
            E, G, _ = R20.energy_grad(Mx, cfg, p, None)
            return E, (G[mask][:, IU[0], IU[1]] * OFF).ravel()

        out[f"polish_unpack_roundtrip_max_abs_{lab}"] = float(
            np.max(np.abs(unpack(x0, M, mask) - M))
        )
        E0, g0 = fun(x0)
        wp = {}
        for e in (1e-3, 1e-4, 1e-5):
            w = 0.0
            for _ in range(3):
                d = rng.standard_normal(x0.shape)
                dd = float(g0 @ d)
                f = [fun(x0 + s * e * d)[0] for s in (-2, -1, 1, 2)]
                fd = (f[0] - 8 * f[1] + 8 * f[2] - f[3]) / (12 * e)
                w = max(w, abs(fd - dd) / max(abs(dd), 1e-300))
            wp[f"{e:g}"] = w
        out[f"a_packed_gradient_stencil4_rel_{lab}"] = wp
    # one tiny job end to end: FIRE 30 steps then the polish (gate loose so that the callback path is exercised)
    cfg16 = cfg_of(16, 24.0, G_MAIN, DELTA_MAIN)
    for bnd in ("seed", "free"):
        M0 = R20.seed_axes(cfg16, lam_of("Sd", DELTA_MAIN))
        Mf, des = descend(M0, cfg16, p, None, 30, 90, f"smoke_{bnd}", bnd != "free", log_every=10)
        Mp, pol = polish(
            Mf,
            cfg16,
            p,
            None,
            bnd != "free",
            f"smoke_{bnd}",
            max_iter=60,
            gate=(2.5 if bnd == "seed" else 0.05),
            log_every=10,
        )
        rd = reads_all(Mp, cfg16, p, None, bnd != "free")
        out[f"tiny_{bnd}"] = {
            "E0": des["E0"],
            "E_fire": des["E_end"],
            "E_polished": pol["E_end"],
            "fmax_fire": des["fmax_end"],
            "fmax_polished": pol["fmax_end"],
            "polish_verdict": pol["verdict"],
            "polish_iters": pol["iters"],
            "face_layer_E": rd["face_layer"]["E"],
            "solid_angle_rank1_r6": rd["solid_angle"]["rank1"]["degree"]["6"],
            "frame3_c0": rd["frame3"]["c0_extrapolated"],
            "axis_norm_min_rel": rd["norm3"]["axis_norm_min_rel"],
        }
    # g 1 accepted by the stack's locus gate on the S_1 seed
    cfg_g1 = cfg_of(16, 24.0, 1.0, DELTA_MAIN)
    E1, G1, info1 = R20.energy_grad(
        R20.seed_axes(cfg_g1, lam_of("S1", DELTA_MAIN)), cfg_g1, params_of(1.0, DELTA_MAIN), None
    )
    out["g1_seed"] = {"E": float(E1), "ok": info1["ok"], "min_gap": info1["min_gap"]}
    cfg_d1 = cfg_of(16, 24.0, G_MAIN, 1.0)
    Ed, Gd, infod = R20.energy_grad(
        R20.seed_axes(cfg_d1, lam_of("S1", 1.0)), cfg_d1, params_of(G_MAIN, 1.0), None
    )
    out["delta1_seed"] = {
        "E": float(Ed),
        "ok": infod["ok"],
        "min_gap": infod["min_gap"],
        "roots": list(R0.roots_of(cfg_d1)),
    }
    for k, v in out.items():
        log(f"smoke {k}: {v}")
    with open(os.path.join(DATA, "m5_32_r21_1_smoke.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


# ================= collect =================
def orient_on_surface(v, c0, k):
    """orient the unit-vector field ON the closed cube surface only (half-width k cells about c0): the maximum spanning
    tree of |v . v'| over the surface's nearest-neighbor links (scipy), the signs propagated from the root along the
    tree, then the conflicts (v . v' < 0) counted over every surface link. Zero conflicts iff the surface field is
    orientable (the degree is then defined); the 3D orientation of the R20 audit can leave a sign cut crossing the
    surface even when the surface field is orientable, so the check is done on the surface itself.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import minimum_spanning_tree, breadth_first_order

    n = v.shape[0]
    idx = np.arange(n)
    I, Jg, K = np.meshgrid(idx, idx, idx, indexing="ij")
    cheb = np.maximum.reduce([np.abs(I - c0), np.abs(Jg - c0), np.abs(K - c0)])
    on = cheb == k
    cells = np.argwhere(on)
    lab = -np.ones((n, n, n), dtype=int)
    lab[on] = np.arange(len(cells))
    rows, cols, w = [], [], []
    for ax in range(3):
        for d in (1,):
            nb = cells.copy()
            nb[:, ax] += d
            ok = nb[:, ax] < n
            ok &= on[nb[ok][:, 0], nb[ok][:, 1], nb[ok][:, 2]] if ok.any() else ok
            a = lab[cells[ok][:, 0], cells[ok][:, 1], cells[ok][:, 2]]
            b = lab[nb[ok][:, 0], nb[ok][:, 1], nb[ok][:, 2]]
            dots = np.einsum(
                "ij,ij->i",
                v[cells[ok][:, 0], cells[ok][:, 1], cells[ok][:, 2]],
                v[nb[ok][:, 0], nb[ok][:, 1], nb[ok][:, 2]],
            )
            rows += list(a)
            cols += list(b)
            w += list(1.0 - np.abs(dots) + 1e-9)
    m = len(cells)
    G = coo_matrix((w, (rows, cols)), shape=(m, m)).tocsr()
    T = minimum_spanning_tree(G).tocoo()
    adj = {}
    for a, b in zip(T.row, T.col):
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    sign = np.zeros(m)
    order, pred = breadth_first_order(T + T.T, 0, directed=False, return_predecessors=True)
    vs = v[cells[:, 0], cells[:, 1], cells[:, 2]]
    sign[order[0]] = 1.0
    for node in order[1:]:
        pr = pred[node]
        sign[node] = sign[pr] * (1.0 if np.dot(vs[pr], vs[node]) >= 0 else -1.0)
    vo = v.copy()
    vo[cells[:, 0], cells[:, 1], cells[:, 2]] = vs * sign[:, None]
    dots = np.einsum(
        "ij,ij->i",
        vo[cells[rows][:, 0], cells[rows][:, 1], cells[rows][:, 2]],
        vo[cells[cols][:, 0], cells[cols][:, 1], cells[cols][:, 2]],
    )
    ncf = int(np.sum(dots < 0.0))
    return vo, ncf, int(len(order) == m)


def degree_surface_reads(M, cfg):
    """the solid-angle degree of each eigenvector on the r 6 / 9 / 12 cubes WITH the surface-orientability check
    (added at collect on 2026-09-14 22:40 UTC after the R21-1 audit: the R20 audit's reader orients in 3D and returns
    an integer on a non-orientable surface too, where the degree is undefined and a disclination pierces the cube;
    here the field is oriented on the surface itself, the way the audit did it)."""
    n, h = cfg["n"], cfg["h"]
    c0 = n // 2
    lam, vec = np.linalg.eigh(M[..., 1:, 1:])
    out = {}
    for r in range(3):
        rec = {}
        for hv in (6.0, 9.0, 12.0):
            k = int(round(hv / h))
            if c0 + k >= n:
                continue
            vo, sc, connected = orient_on_surface(vec[..., :, r], c0, k)
            deg = float(R20A.solid_angle_degree(vo, c0, k))
            if sc == 0:
                reading = (
                    "kept"
                    if abs(abs(deg) - 1.0) < 0.1
                    else ("lost" if abs(deg) < 0.5 else "fractional")
                )
            else:
                reading = "undefined (non-orientable: a disclination pierces the cube)"
            rec[f"{hv:g}"] = {
                "surface_conflicts": sc,
                "surface_connected": connected,
                "solid_angle": deg,
                "degree": (deg if sc == 0 else None),
                "reading": reading,
            }
        out[f"rank{r}"] = rec
    return out


def backfill_degree_surface(J):
    changed = False
    for tag, r in J["rows"].items():
        if r.get("status") != "OK" or "end_reads" not in r or "degree_surface" in r["end_reads"]:
            continue
        f = os.path.join(OUT_NPZ, tag + ".npz")
        if not os.path.exists(f):
            continue
        cfg = cfg_of(r["n"], r["L"], r["g"], r["delta"])
        r["end_reads"]["degree_surface"] = degree_surface_reads(np.load(f)["M"], cfg)
        ff = os.path.join(OUT_NPZ, tag + "_fire.npz")
        if os.path.exists(ff) and "fire_reads" in r:
            r["fire_reads"]["degree_surface"] = degree_surface_reads(np.load(ff)["M"], cfg)
        wr = WINDING_RANK[r["obj"]]
        log(
            f"degree_surface {tag}: winding rank {wr} on r 9: {r['end_reads']['degree_surface'][f'rank{wr}'].get('9')}"
        )
        changed = True
    return changed


def _wind9(r, obj):
    """(reading, degree or None, surface conflicts) of the winding eigenvector on the r 9 cube, surface-checked."""
    q = r["end_reads"]["degree_surface"][f"rank{WINDING_RANK[obj]}"]["9"]
    return [q["reading"], q["degree"], q["surface_conflicts"]]


def _polished_E(r):
    return r.get("E")


def _deg9(r, obj):
    wr = WINDING_RANK[obj]
    return r["end_reads"]["solid_angle"][f"rank{wr}"]["degree"]["9"]


def _fit_root_law(xs, cs):
    """least squares of log |c| = log A + 2 log |x - x0| over x0 (a 1D scan then a polish); returns (A, x0, max rel miss)."""
    xs, cs = np.asarray(xs, float), np.asarray(cs, float)
    ac = np.abs(cs)

    def miss(x0):
        w = np.log(np.maximum(np.abs(xs - x0), 1e-9))
        logA = np.mean(np.log(np.maximum(ac, 1e-300)) - 2 * w)
        pred = np.exp(logA) * (xs - x0) ** 2
        return logA, float(np.max(np.abs(pred - ac) / np.maximum(ac, 0.2 * np.median(ac))))

    grid = np.linspace(min(xs) - 2.0, max(xs) + 2.0, 4001)
    best = min(grid, key=lambda x0: miss(x0)[1])
    from scipy.optimize import minimize_scalar

    rs = minimize_scalar(lambda x0: miss(x0)[1], bracket=(best - 0.05, best, best + 0.05))
    x0 = float(rs.x) if rs.success and miss(rs.x)[1] <= miss(best)[1] else float(best)
    logA, mx = miss(x0)
    return float(np.exp(logA)), x0, mx


def collect():
    J = load_json_all()
    backfill_degree_surface(J)
    rows = J["rows"]
    with open(R20_JSON) as f:
        R20J = json.load(f)
    res = {"per_axis": {}, "string": {}, "box": {}, "triple": {}, "w1_ladder": {}, "frame": {}}
    Lbox = 48.0

    def row(obj, bnd, g=G_MAIN, d=DELTA_MAIN, w=1.0, n=32):
        return rows.get(job_tag(dict(obj=obj, bnd=bnd, g=g, delta=d, w1s=w, n=n)))

    # per axis: the boundary triple at equal residual
    bdiff = {}
    for obj in ("S1", "Sd", "S0"):
        tri = {b: row(obj, b) for b in ("seed", "far", "free")}
        rec = {
            "E": {},
            "label": {},
            "fmax_end": {},
            "E_lt_12": {},
            "degree_r9": {},
            "E_fire": {},
            "shell_E": {},
            "face_layer_E": {},
            "tail_slope": {},
            "tension": {},
        }
        for b, r in tri.items():
            if r is None or r.get("E") is None:
                continue
            rec["E"][b] = r["E"]
            rec["label"][b] = r["label"]
            rec["fmax_end"][b] = r["polish"]["fmax_end"]
            rec.setdefault("winding_r9", {})[b] = _wind9(r, obj)
            rec.setdefault("fmax_fire_end", {})[b] = r.get("descent", {}).get("fmax_end")
            rec.setdefault("degree_r9_by_rank", {})[b] = {
                k: v["degree"]["9"] for k, v in r["end_reads"]["solid_angle"].items()
            }
            rec.setdefault("mermin_r9_winding_rank", {})[b] = r["end_reads"]["degree"][
                f"rank{WINDING_RANK[obj]}"
            ]["flux"]["9"]
            rec.setdefault("conflicts_winding_rank", {})[b] = r["end_reads"]["solid_angle"][
                f"rank{WINDING_RANK[obj]}"
            ]["conflicts"]
            rec.setdefault("shell_lam_r9", {})[b] = r["end_reads"]["shells"]["9"]["lam_mean"]
            rec.setdefault("virial", {})[b] = r["end_reads"]["derrick"]["virial_E_curv_over_V"]
            rec.setdefault("R_star", {})[b] = r["end_reads"]["derrick"]["R_star_from_virial"]
            rec.setdefault("dE_dlambda", {})[b] = r["end_reads"]["derrick"]["order3"][
                "dE_dlambda_at_1"
            ]
            rec.setdefault("polish_iters", {})[b] = r["polish"].get("iters")
            rec.setdefault("E_curv_V", {})[b] = [
                r["end_reads"]["energy"]["E_curv"],
                r["end_reads"]["energy"]["V"],
            ]
            rec["E_lt_12"][b] = r["end_reads"]["radial"]["E_lt_R"]["12"]
            rec["degree_r9"][b] = _deg9(r, obj)
            rec["E_fire"][b] = r.get("E_fire")
            rec["shell_E"][b] = r["end_reads"]["radial"]["E_pin_shell"]
            rec["face_layer_E"][b] = r["end_reads"]["face_layer"]["E"]
            rec["tail_slope"][b] = r["end_reads"]["tail"]["slope_density"]
            rec["tension"][b] = r["end_reads"]["string_tension_read"]
        outc = []
        if len(rec["E"]) == 3:
            Es = rec["E"]
            e12 = rec["E_lt_12"]
            rel12 = (max(e12.values()) - min(e12.values())) / max(
                abs(np.mean(list(e12.values()))), 1e-300
            )
            wind_all = all(w[0] == "kept" for w in rec["winding_r9"].values())
            wind_lost = any(
                w[0] != "kept" for w in rec["winding_r9"].values()
            )  # lost, or undefined
            rec["residual_spread_max_over_min"] = float(
                max(rec["fmax_end"].values()) / min(rec["fmax_end"].values())
            )
            reldiff = max(abs(Es["far"] - Es["seed"]), abs(Es["free"] - Es["seed"])) / max(
                abs(Es["seed"]), 1e-300
            )
            rec["rel_spread_E_lt_12"] = float(rel12)
            rec["rel_boundary_diff_total"] = float(reldiff)
            rec["max_boundary_abs_diff"] = float(
                max(
                    abs(Es["far"] - Es["seed"]),
                    abs(Es["free"] - Es["seed"]),
                    abs(Es["free"] - Es["far"]),
                )
            )
            bdiff[obj] = rec["max_boundary_abs_diff"]
            if rel12 <= 0.02 and wind_all:
                outc.append("INTERIOR_SET")
            if reldiff > 0.05 or wind_lost:
                outc.append("BOUNDARY_SET")
            if not outc:
                outc.append("BOUNDARY_UNDECIDED")
            unpol = [b for b, lab in rec["label"].items() if lab != "POLISHED"]
            if unpol:
                outc = [
                    o
                    + f" (at max |G| {min(rec['fmax_end'].values()):.1e} to {max(rec['fmax_end'].values()):.1e}, not one residual)"
                    for o in outc
                ]
        rec["outcomes"] = outc
        res["per_axis"][obj] = rec
    # the string
    sd = {b: row("Sd", b) for b in ("seed", "far", "free")}
    if sd["free"] is not None and sd["free"].get("E") is not None:
        rf = sd["free"]
        rs_ = sd["seed"]
        t_free = rf["end_reads"]["string_tension_read"]
        t_seed = (
            rs_["end_reads"]["string_tension_read"] if rs_ and rs_.get("E") is not None else None
        )
        wreading, deg, wsc = _wind9(rf, "Sd")
        deg = deg if deg is not None else 0.0
        rec = {
            "tension_free": t_free,
            "tension_seed": t_seed,
            "degree_r9_free": deg,
            "winding_r9_free": [wreading, wsc],
            "winding_r9_by_boundary": {
                b: _wind9(sd[b], "Sd")
                for b in ("seed", "far", "free")
                if sd[b] is not None and sd[b].get("E") is not None
            },
            "face_layer_E_free": rf["end_reads"]["face_layer"]["E"],
            "axis_norm_min_rel_free": rf["end_reads"]["norm3"]["axis_norm_min_rel"],
            "axis_crossings_free": rf["end_reads"]["norm3"]["axis_crossings_gap_lt_0.02"],
            "axis_gap_min_free": rf["end_reads"]["norm3"]["axis_gap_min"],
        }
        outc = []
        lost = wreading != "kept"
        if t_free < 0.02 or lost:
            outc.append("STRING_ESCAPES (by the tension; the winding " + wreading + ")")
        elif (
            t_seed is not None and abs(t_free - t_seed) <= 0.3 * abs(t_seed) and wreading == "kept"
        ):
            outc.append("STRING_HELD")
        else:
            outc.append("STRING_UNDECIDED")
        rec["outcomes"] = outc
        res["string"] = rec
    # the escape route on any axis and boundary where the winding was lost
    routes = {}
    for obj in ("S1", "Sd", "S0"):
        for b in ("seed", "far", "free"):
            r = row(obj, b)
            if r is None or r.get("E") is None:
                continue
            wreading, wdeg, wsc = _wind9(r, obj)
            if wreading != "kept":
                nr = r["end_reads"]["norm3"]
                route = (
                    "ESCAPE_BY_EXCHANGE"
                    if (nr["axis_norm_min_rel"] >= 0.9 and nr["axis_crossings_gap_lt_0.02"] > 0)
                    else (
                        "ESCAPE_BY_MELTING"
                        if nr["axis_norm_min_rel"] < 0.7
                        else "ESCAPE_ROUTE_UNDECIDED"
                    )
                )
                # added at the read (2026-09-14 19:30 UTC): the winding can move to another rank when the eigenvalues cross (the
                # exchange route read directly on the degrees, beside the norm rule the plan pre-registered)
                ds = r["end_reads"]["degree_surface"]
                degs = {
                    k: [v["9"]["reading"], v["9"]["degree"], v["9"]["surface_conflicts"]]
                    for k, v in ds.items()
                }
                carrier = [
                    k
                    for k, v in degs.items()
                    if v[0] == "kept" and k != f"rank{WINDING_RANK[obj]}"
                ]
                routes[f"{obj}_B{b}"] = {
                    "route_by_norm_rule": route,
                    "winding_r9": [wreading, wdeg, wsc],
                    "reading": (
                        f"a disclination pierces the r 9 cube (the winding eigenvector non-orientable on the surface, {wsc} surface conflicts): a disclination texture (the non-orientability read; its biaxiality not read)"
                        if wreading.startswith("undefined")
                        else (
                            ("the winding sits on " + ",".join(carrier))
                            if carrier
                            else "no eigenvector carries a winding on the r 9 cube"
                        )
                    ),
                    "axis_norm_min_rel": nr["axis_norm_min_rel"],
                    "axis_crossings": nr["axis_crossings_gap_lt_0.02"],
                    "degree_r9_by_rank_surface_checked": degs,
                    "winding_transferred_to": carrier,
                    "route_by_rank_transfer": (
                        (
                            "EXCHANGE_BY_RANK_TRANSFER (the winding sits on "
                            + ",".join(carrier)
                            + ")"
                        )
                        if carrier
                        else "no rank carries the winding"
                    ),
                }
    res["escape_routes"] = routes
    # the box, per axis under B_free
    for obj in ("S1", "Sd", "S0"):
        r = row(obj, "free")
        if r is None or r.get("E") is None:
            continue
        d = r["end_reads"]["derrick"]
        en = r["end_reads"]["energy"]
        rec = {
            "dE_dlambda_scan": d["order3"]["dE_dlambda_at_1"],
            "dE_dlambda_identity": float(-en["E_curv"] + 3.0 * en["V"]),
            "R_star": d["R_star_from_virial"],
            "virial": d["virial_E_curv_over_V"],
            "r_half": d["r_half_energy"],
            "note": "the identity dE / d lambda = -E_curv + 3 V is the derivative of a pure dilation of a field at the vacuum outside the dilated region; on these box-filling fields (not at the vacuum anywhere) the scan and the identity measure different things and disagree in sign, the interior carrying 82 to 97 percent of the scan's negative slope (the R21-1 audit, C5e, C5g, C5h)",
        }
        same_sign = np.sign(rec["dE_dlambda_scan"]) == np.sign(rec["dE_dlambda_identity"])
        if (
            rec["dE_dlambda_identity"] < 0
            and rec["dE_dlambda_scan"] < 0
            and rec["R_star"] > Lbox / 4
        ):
            rec["outcome"] = "AXIS_BOX_LIMITED"
        elif rec["R_star"] < Lbox / 6 and abs(rec["virial"] - 3.0) <= 0.9:
            rec["outcome"] = "AXIS_LOCALIZED"
        elif not same_sign:
            rec["outcome"] = (
                "BOX_UNDECIDED (the dilation scan and the Derrick identity disagree in sign: the field fills the box; the virial below 3 says the potential's share exceeds the balance)"
            )
        else:
            rec["outcome"] = "BOX_UNDECIDED"
        res["box"][obj] = rec
    # the triple on the polished B_free energies
    Ef = {
        obj: row(obj, "free")["E"]
        for obj in ("S1", "Sd", "S0")
        if row(obj, "free") and row(obj, "free").get("E") is not None
    }
    if len(Ef) == 3 and bdiff:
        Es = np.array([Ef[o] for o in ("S1", "Sd", "S0")])
        diffs = [abs(Es[i] - Es[j]) for i in range(3) for j in range(i + 1, 3)]
        resol = max(bdiff.values())
        trip = {
            "E_free": Ef,
            "pairwise_abs_diffs": diffs,
            "resolution_3x_max_boundary_diff": 3.0 * resol,
            "outcome": "THREE_AXES_DISTINCT" if min(diffs) > 3.0 * resol else "AXES_DEGENERATE",
        }
        if np.all(Es > 0):
            srt = np.sort(Es)
            trip["koide_Q"] = float(np.sum(Es) / np.sum(np.sqrt(Es)) ** 2)
            trip["ratios_sorted"] = [float(x / srt[0]) for x in srt]
        res["triple"] = trip
    # the W1 ladder
    lad = {}
    r20s1 = row("S1", "seed")
    base_dE = (
        r20s1["end_reads"]["derrick"]["order3"]["dE_dlambda_at_1"]
        if r20s1 and r20s1.get("E") is not None
        else None
    )
    for obj, w, b in (
        ("S1", 1.0, "seed"),
        ("S1", 5.0, "seed"),
        ("S1", 25.0, "seed"),
        ("S1", 25.0, "free"),
        ("Sd", 1.0, "seed"),
        ("Sd", 25.0, "seed"),
    ):
        r = row(obj, b, w=w)
        if r is None or r.get("E") is None:
            continue
        d = r["end_reads"]["derrick"]
        lad[f"{obj}_w{w:g}_B{b}"] = {
            "E": r["E"],
            "label": r["label"],
            "R_star": d["R_star_from_virial"],
            "virial": d["virial_E_curv_over_V"],
            "dE_dlambda": d["order3"]["dE_dlambda_at_1"],
            "dE_dlambda_identity": float(
                -r["end_reads"]["energy"]["E_curv"] + 3.0 * r["end_reads"]["energy"]["V"]
            ),
            "winding_r9": _wind9(r, obj),
            "r_half": d["r_half_energy"],
            "E_lt_12": r["end_reads"]["radial"]["E_lt_R"]["12"],
            "shell_share": r["end_reads"]["radial"]["E_pin_shell"] / r["E"],
            "E_curv": r["end_reads"]["energy"]["E_curv"],
            "V": r["end_reads"]["energy"]["V"],
            "degree_r9": _deg9(r, obj),
        }
    k25 = "S1_w25_Bseed"
    if k25 in lad:
        q = lad[k25]
        if (
            q["R_star"] < Lbox / 6
            and abs(q["virial"] - 3.0) <= 0.9
            and base_dE is not None
            and abs(q["dE_dlambda"]) <= 0.1 * abs(base_dE)
        ):
            lad["outcome"] = "SCALE_IS_THE_POTENTIALS"
        elif q["R_star"] > Lbox / 4 or q["virial"] > 6.0:
            lad["outcome"] = "SCALE_IS_THE_BOXS"
        else:
            lad["outcome"] = "SCALE_UNDECIDED"
        lad["R20_S1_dE_dlambda_reference"] = base_dE
    res["w1_ladder"] = lad
    # the frame saddle: the s -> 0 curvatures on the polished rows
    gpts, dpts = {}, {}
    for g in (1.0, 2.0, 4.0, 8.0, 16.0, 32.0):
        r = row("S1", "seed", g=g)
        if r is not None and r.get("E") is not None and "frame3" in r["end_reads"]:
            gpts[g] = {
                "c0": r["end_reads"]["frame3"]["c0_extrapolated"],
                "c_0.02": r["end_reads"]["frame3"]["d2E_ds2"]["0.02"],
                "label": r["label"],
                "E": r["E"],
                "dV_max": max(r["end_reads"]["frame3"]["dV"].values()),
                "winding_r9": _wind9(r, "S1"),
            }
    for d in (0.3, 0.6, 1.0):
        r = row("S1", "seed", d=d)
        if r is not None and r.get("E") is not None and "frame3" in r["end_reads"]:
            dpts[d] = {
                "c0": r["end_reads"]["frame3"]["c0_extrapolated"],
                "c_0.02": r["end_reads"]["frame3"]["d2E_ds2"]["0.02"],
                "label": r["label"],
                "E": r["E"],
                "winding_r9": _wind9(r, "S1"),
            }
    fr = {
        "g_rows": {f"{g:g}": v for g, v in gpts.items()},
        "delta_rows": {f"{d:g}": v for d, v in dpts.items()},
    }
    if len(gpts) >= 4:
        A, g0, mg = _fit_root_law(list(gpts), [v["c0"] for v in gpts.values()])
        fr["fit_g"] = {"A": A, "g0": g0, "max_rel_miss": mg, "law": "c = A (g - g0)^2"}
        fr["ratio_g32_over_g8"] = (
            gpts[32.0]["c0"] / gpts[8.0]["c0"]
            if 32.0 in gpts and 8.0 in gpts and gpts[8.0]["c0"] != 0
            else None
        )
    if len(dpts) >= 3:
        Ad, d0, md = _fit_root_law(list(dpts), [v["c0"] for v in dpts.values()])
        fr["fit_delta"] = {
            "A": Ad,
            "delta0": d0,
            "max_rel_miss": md,
            "law": "c = A (delta - delta0)^2",
        }
    if "fit_g" in fr and "fit_delta" in fr:
        g0, d0 = fr["fit_g"]["g0"], fr["fit_delta"]["delta0"]
        xs = [(g, 0.3, gpts[g]["c0"]) for g in gpts] + [
            (8.0, d, dpts[d]["c0"]) for d in dpts if d != 0.3
        ]
        ac = np.array([abs(c) for _, _, c in xs])
        basis = np.array([(g - g0) ** 2 * (d - d0) ** 2 for g, d, _ in xs])
        C = float(
            np.exp(np.mean(np.log(np.maximum(ac, 1e-300)) - np.log(np.maximum(basis, 1e-300))))
        )
        pred = C * basis
        miss = np.abs(pred - ac) / np.maximum(ac, 0.2 * np.median(ac))
        fr["product_fit"] = {
            "C": C,
            "g0": g0,
            "delta0": d0,
            "rows": [
                {"g": g, "delta": d, "c0": c, "pred": float(-pr), "rel_miss": float(m)}
                for (g, d, c), pr, m in zip(xs, pred, miss)
            ],
            "max_rel_miss": float(np.max(miss)),
        }
        g_at_1, d_at_1, g_at_0 = abs(g0 - 1.0) <= 0.3, abs(d0 - 1.0) <= 0.15, abs(g0) <= 0.3
        if fr["product_fit"]["max_rel_miss"] > 0.2:
            fr["outcome"] = "FRAME_NOT_PRODUCT"
        elif g_at_1 and d_at_1:
            fr["outcome"] = "FRAME_TWO_ROOTS"
        elif (g_at_1 or d_at_1) or g_at_0:
            fr["outcome"] = "FRAME_ONE_ROOT" + (" (the pure g^2 law)" if g_at_0 else "")
        else:
            fr["outcome"] = "FRAME_ROOTS_ELSEWHERE"
        unpol = [f"g{g:g}" for g in gpts if gpts[g]["label"] != "POLISHED"] + [
            f"d{d:g}" for d in dpts if dpts[d]["label"] != "POLISHED"
        ]
        if unpol:
            fr["outcome"] += f" (rows not POLISHED: {','.join(unpol)})"
    res["frame"] = fr
    # the optional n64 row
    r64 = row("Sd", "free", n=64)
    if r64 is not None and r64.get("E") is not None:
        res["n64_string"] = {
            "E": r64["E"],
            "E_fire": r64.get("E_fire"),
            "label": r64["label"],
            "tension": r64["end_reads"]["string_tension_read"],
            "winding_r9": _wind9(r64, "Sd"),
            "face_layer_E": r64["end_reads"]["face_layer"]["E"],
            "axis_norm_min_rel": r64["end_reads"]["norm3"]["axis_norm_min_rel"],
        }
    # the R20 reference energies at the FIRE bar
    res["R20_reference"] = {
        t: R20J["rows"][t]["E"]
        for t in (
            "S1_v4std_n32_g8_x4500",
            "Sd_v4std_n32_g8_x4500",
            "S0_v4std_n32_g8_x4500",
            "S1_v4std_n32_g32",
        )
        if t in R20J["rows"]
    }
    res["collected_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    J["results"] = res
    save_json(J)
    plots(J)
    print(json.dumps(res, indent=1, default=str))
    return res


def plots(J):
    rows = J["rows"]
    os.makedirs(PLOTS, exist_ok=True)

    def row(obj, bnd, g=G_MAIN, d=DELTA_MAIN, w=1.0, n=32):
        return rows.get(job_tag(dict(obj=obj, bnd=bnd, g=g, delta=d, w1s=w, n=n)))

    # 1: the polished energies per boundary and axis, with E(< 12)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    xb = {"seed": 0, "far": 1, "free": 2}
    for obj, mk in (("S1", "o"), ("Sd", "s"), ("S0", "^")):
        xs, Es, e12 = [], [], []
        for b in ("seed", "far", "free"):
            r = row(obj, b)
            if r and r.get("E") is not None:
                xs.append(xb[b])
                Es.append(r["E"])
                e12.append(r["end_reads"]["radial"]["E_lt_R"]["12"])
        if xs:
            axs[0].plot(xs, Es, mk + "-", label=f"{obj} {lam_of(obj, DELTA_MAIN)}")
            axs[1].plot(xs, e12, mk + "-", label=obj)
    for ax, t in zip(axs, ("E polished (max |G| < 1e-3)", "E(< 12) polished")):
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["B_seed", "B_far", "B_free"])
        ax.set_ylabel(t)
        ax.legend(fontsize=8)
    axs[0].set_title(
        "R21-1 / R21-2: the three axes under three boundaries at equal residual, n32 L48",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r21_1_boundaries.png"), dpi=130)
    # 2: the descents and polishes
    fig, ax = plt.subplots(figsize=(8, 4.4))
    for t, r in rows.items():
        if r.get("n") != 32 or r.get("stage", "").startswith("R21-3"):
            continue
        tr = r.get("descent", {}).get("trace", [])
        if tr:
            ax.plot([q["acc"] for q in tr], [q["E"] for q in tr], lw=0.9, label=t)
        pt = r.get("polish", {}).get("trace", [])
        if pt:
            off = r.get("descent", {}).get("accepted", 0)
            ax.plot([off + q["it"] for q in pt], [q["E"] for q in pt], lw=0.9, ls="--")
    ax.set_xlabel("accepted FIRE steps, then L-BFGS iterations (dashed)")
    ax.set_ylabel("E")
    ax.set_title("R21-1 / R21-4 descents and polishes", fontsize=9)
    ax.legend(fontsize=5)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r21_1_descents.png"), dpi=130)
    # 3: E(< R) per boundary per axis
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.0))
    for ax, obj in zip(axs, ("S1", "Sd", "S0")):
        for b in ("seed", "far", "free"):
            r = row(obj, b)
            if r and r.get("E") is not None:
                pr = r["end_reads"]["radial"]["E_lt_R"]
                Rs = sorted(pr, key=float)
                ax.plot(
                    [float(R) for R in Rs],
                    [pr[R] for R in Rs],
                    "o-",
                    ms=3,
                    label=f"B_{b} (E {r['E']:.2f}, {r['label'][:8]})",
                )
        ax.set_xlabel("R")
        ax.set_ylabel("E(< R)")
        ax.set_title(f"{obj}: the energy inside R per boundary", fontsize=9)
        ax.legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r21_1_radial.png"), dpi=130)
    # 4: the spatial-block norm along the former polar axis
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.0))
    for ax, obj in zip(axs, ("S1", "Sd", "S0")):
        for b in ("seed", "far", "free"):
            r = row(obj, b)
            if r and r.get("E") is not None:
                pr = r["end_reads"]["norm3"]["axis_profile"]
                ax.plot([q["z"] for q in pr], [q["norm"] for q in pr], "-", label=f"B_{b}")
        ax.axhline(1.0 + DELTA_MAIN**2, color="k", lw=0.6, ls=":")
        ax.set_xlabel("z on the polar axis")
        ax.set_ylabel("tr(M_3x3^2)")
        ax.set_title(f"{obj}: the block norm along the axis (vacuum 1.09 dotted)", fontsize=9)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r21_1_axis_norm.png"), dpi=130)
    # 5: the frame curvature against g and delta
    fr = J.get("results", {}).get("frame", {})
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
    gr = fr.get("g_rows", {})
    if gr:
        gs = sorted(float(k) for k in gr)
        axs[0].plot(
            gs, [-gr[f"{g:g}"]["c0"] for g in gs], "o", label="-c0 (s -> 0), polished rows"
        )
        if "fit_g" in fr:
            gg = np.linspace(0.5, 34, 300)
            axs[0].plot(
                gg,
                fr["fit_g"]["A"] * (gg - fr["fit_g"]["g0"]) ** 2,
                "-",
                lw=0.8,
                label=f"A (g - {fr['fit_g']['g0']:.2f})^2",
            )
        axs[0].set_xscale("log")
        axs[0].set_yscale("log")
        axs[0].set_xlabel("g (delta 0.3)")
        axs[0].set_ylabel("-d2E/ds2 along the boost")
        axs[0].legend(fontsize=8)
    dr = fr.get("delta_rows", {})
    if dr:
        ds = sorted(float(k) for k in dr)
        axs[1].plot(ds, [-dr[f"{d:g}"]["c0"] for d in ds], "s", label="-c0 (s -> 0)")
        if "fit_delta" in fr:
            dd = np.linspace(0.2, 1.1, 200)
            axs[1].plot(
                dd,
                fr["fit_delta"]["A"] * (dd - fr["fit_delta"]["delta0"]) ** 2,
                "-",
                lw=0.8,
                label=f"A (delta - {fr['fit_delta']['delta0']:.2f})^2",
            )
        axs[1].set_xlabel("delta (g 8)")
        axs[1].set_ylabel("-d2E/ds2")
        axs[1].legend(fontsize=8)
    axs[0].set_title("R21-3: the frame saddle's curvature, the roots fitted free", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_32_r21_1_frame.png"), dpi=130)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=("smoke", "main", "n64", "collect"))
    ap.add_argument(
        "--json-suffix",
        default="",
        help="a side pool writes m5_32_r21_1_runs_<suffix>.json (merged at collect)",
    )
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--mem-budget", type=float, default=12.0)
    ap.add_argument(
        "--only", default="", help="comma-separated job tags to run (a subset of the stage)"
    )
    a = ap.parse_args()
    if a.json_suffix:
        global OUT_JSON
        OUT_JSON = os.path.join(DATA, f"m5_32_r21_1_runs_{a.json_suffix}.json")
    if a.stage == "smoke":
        smoke()
    elif a.stage == "collect":
        collect()
    else:
        jobs = jobs_main() if a.stage == "main" else jobs_n64()
        if a.only:
            keep = set(a.only.split(","))
            jobs = [j for j in jobs if job_tag(j) in keep]
        run_pool(jobs, a.workers, a.stage, a.mem_budget)


if __name__ == "__main__":
    main()
