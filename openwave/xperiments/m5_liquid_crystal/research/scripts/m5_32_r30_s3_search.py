"""M5.32 R30: S3, the global search of the Lorentz-safe vortex nucleon on the author's instrument (the search kit of
2026-10-26, S3), run with the author's scripts unchanged. The driver: R30-0 (the controls), R30-1 (a) the differential
evolution per structural cell, (b) the multistart sequential LP, and the consolidation of the evaluation history.
The front, the fine-grid replay, the LP prices and the gate word are `m5_32_r30_s3_front.py`.

EQUATIONS FIRST (the author's instrument, consumed read-only from a scratch copy whose hashes are checked at R30-0 (a))
---------------
State: the ten coefficients of V48's `ell` mode, theta = (ln f, ln c_total, ln K, ln b, ln kappa, ln tau, ln L_s,
ln d_G, ln eta, ell), mapped by `V48_covariant_omega_refit.params` to `vortex_v3.par_v3(..., drag = 1.0)`, the comoving
clock; the structural cell enters through `vortex_v3`'s option keys `npow` (the power n of the strain invariant) and
`Gpow` (the lock power of G = z^2 + d_G (z^2 - 1)^Gpow), with `tscale` the start scale the author's V52 uses for n = 4.
Score: `vortex_v3.evaluate(p, N = 320, X = 11.0)` (the coarse grid of V48), cold from the v2.6 profile the kit ships;
    E15 = max_r |100 (O_r / REF_r - 1)| over the fifteen rows of `vortex_v3.ROWS` against `vortex_radial.REFS`
    (the frozen isoscalar target r2_ES = 0.596962357 fm^2);
    the heavy-limit omega proxy  P = pi^3 c_d m_omega^2 / f^2,  c_d = (1 - eta) c_total,  m_omega = 783 MeV;
    the CFD phase gate  max |gate_dev(p)| (the tree I = 0 S-wave phase minus the CFD phase, 400 to 800 MeV);
    the photon bound  K_gamma = 1/2 - e^2 ell^2 max[(z^2)']^2 / (4 pi^4 c_d).
Merit (the DE objective) = E15 + RHO max(0, W - P) + PEN_GATE max(0, gate - 5.02) + PEN_KG max(0, 0.05 - K_gamma),
    RHO = 10 and W = 5.8 as V48; the two hard penalties keep the DE inside V48's constraints; a failed solve scores
    FAIL_MERIT and is counted. Feasible = gate <= 5.02 (V48's acceptance), K_gamma >= 0.05, eta <= 0.9 (the box),
    P >= 5.8. The finite-range g_omega^2 / 4 pi (V48b's radial Fourier transform) is computed at the front only.
Box (the declared space; a statement about this box and nothing outside it): f [84, 101] MeV; c_total, b, kappa, L_s
    within factors 3, 4, 5, 4 of the shipped V48 state; K [5, 60]; tau [3, 5000]; d_G [0.05, 1.2]; eta [0.02, 0.9];
    ell [-0.15, 0.20]. Cells: n in {3, 4} x Gpow in {2, 4}; tscale 1 for n = 3 and the author's V52 value for n = 4.
Search: SciPy differential evolution (best1bin, 15 members per dimension, Latin-hypercube start with the shipped V48
    state injected, and the matching V52 state where one exists, deferred updating, no polish, fixed seed), stopped
    when the merit spread of the lowest 90 percent of the population stays below 0.01 points for 20 generations or at
    the cap; the population is checkpointed every generation and a rerun resumes from it. Then the author's `slp()`
    (elastic sequential LP) imported unchanged from `V48_covariant_omega_refit.py`, from the best DE members and from
    Latin-hypercube starts (the W ladder 3.0 -> 4.5 -> 5.8 of V48 for a cold start, W = 5.8 for a DE member).

Modes (python3 m5_32_r30_s3_search.py MODE ..., with the scratch copy of the kit in $M5_32_S3_KIT):
    controls                       R30-0 (a) hashes, (b) the shipped V48 state coarse and fine, (c) V55 prices, (d) timing
                                   and the solve-failure rate on 200 Latin-hypercube points
    de CELL [--workers W] [--maxiter G]   the differential evolution of one cell (resumable)
    slp [--workers W]              the multistart SLP in every cell that has a DE record
    consolidate                    the history jsonl files -> data/m5_32_r30_s3/history.npz (local) + the summary JSON
One thread per worker (the four thread variables are set to 1 before NumPy is imported); the pool size is the
machine's, never written here.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import argparse
import contextlib
import glob
import hashlib
import io
import json
import multiprocessing as mp
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_DIR = os.path.join(DATA, "m5_32_r30_s3")
OUT_JSON = os.path.join(DATA, "m5_32_r30_s3.json")
THEORY = os.path.join(HERE, "..", "..", "theory", "duda_2026-10-27_ldgs_scorecard_v1_2", "package", "LdGS_vacuum_sectors")
KIT = os.environ.get("M5_32_S3_KIT")
KIT_FILES = ["vortex_radial.py", "vortex_v3.py", "V48_covariant_omega_refit.py", "V48b_readout.py", "V54b_replay_fine.py",
             "V55_shadow_prices.py", "V24_alpha_particle.py", "V20_v26_profile.npz", "V44_comoving_clock_refit.json",
             "V48_covariant_omega_refit_ell.json", "V52_refit_n4.json", "V52_refit_lock.json"]
T0 = time.time()

NG, XG = 320, 11.0                   # V48's coarse grid
W, RHO = 5.8, 10.0                   # V48's omega gate and merit weight
GATE_ACC, KG_MIN, ETA_MAX = 5.02, 0.05, 0.9
PEN_GATE, PEN_KG, FAIL_MERIT = 20.0, 200.0, 1e3
LOGP = ["f", "c_total", "K", "b", "kappa", "tau", "Ls", "dG", "eta"]
FREE = LOGP + ["ell"]
CELLS = {"n3g2": dict(npow=3, Gpow=2, tscale=1.0), "n4g2": dict(npow=4, Gpow=2, tscale=0.022858),
         "n3g4": dict(npow=3, Gpow=4, tscale=1.0), "n4g4": dict(npow=4, Gpow=4, tscale=0.022858)}
MAIN_CELL = "n3g2"
POP_PER_DIM, CAP_MAIN, CAP_VARIANT = 15, 200, 100
CONV_SPREAD, CONV_GENS, CONV_FRACTION = 0.01, 20, 0.9
SEED = 20261029


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


# ================= the author's instrument, loaded from the scratch copy =================
def _exec_part(fname, cut, name):
    """the function part of an author script that runs at module level (V48b, V55): the source up to `cut`, executed
    in its own namespace (the author's code unchanged; the cut point is recorded in the controls)."""
    src = open(os.path.join(KIT, fname)).read()
    head = src.split(cut)[0]
    ns = {"__name__": name}
    exec(compile(head, fname, "exec"), ns)
    return ns, len(head.splitlines())


def load_kit():
    if not KIT or not os.path.isdir(KIT):
        raise SystemExit("set M5_32_S3_KIT to the scratch copy of LdGS_vacuum_sectors (never the theory copy: the scripts write into their cwd)")
    os.chdir(KIT)
    if KIT not in sys.path:
        sys.path.insert(0, KIT)
    import vortex_radial, vortex_v3                                     # noqa: E402
    argv = sys.argv; sys.argv = ["V48_covariant_omega_refit.py", "ell"]  # V48 reads its mode from argv at import
    try:
        import V48_covariant_omega_refit as V48                         # noqa: E402
    finally:
        sys.argv = argv
    v48b, n48b = _exec_part("V48b_readout.py", "\nout = {}\n", "V48b_part")
    v55, n55 = _exec_part("V55_shadow_prices.py", "\nimport sys\ntry:\n", "V55_part")
    return vortex_radial, vortex_v3, V48, v48b, v55, {"V48b_readout.py": n48b, "V55_shadow_prices.py": n55}


VR, V3, V48, V48B, V55, CUTS = load_kit()
assert V48.MODE == "ell" and V48.FREE == FREE
SHIPPED = json.load(open(os.path.join(KIT, "V48_covariant_omega_refit_ell.json")))["5.8"]
V52 = {t: json.load(open(os.path.join(KIT, f"V52_refit_{t}.json")))["5.8"] for t in ("n4", "lock")}


def theta_of(v):
    return np.array([np.log(v[k]) for k in LOGP] + [v.get("ell", 0.0)])


THETA_V48 = theta_of(SHIPPED["params"])
BOX = {"f": (np.log(84.0), np.log(101.0)), "c_total": (np.log(SHIPPED["params"]["c_total"] / 3), np.log(SHIPPED["params"]["c_total"] * 3)),
       "K": (np.log(5.0), np.log(60.0)), "b": (np.log(SHIPPED["params"]["b"] / 4), np.log(SHIPPED["params"]["b"] * 4)),
       "kappa": (np.log(SHIPPED["params"]["kappa"] / 5), np.log(SHIPPED["params"]["kappa"] * 5)), "tau": (np.log(3.0), np.log(5000.0)),
       "Ls": (np.log(SHIPPED["params"]["Ls"] / 4), np.log(SHIPPED["params"]["Ls"] * 4)), "dG": (np.log(0.05), np.log(1.2)),
       "eta": (np.log(0.02), np.log(ETA_MAX)), "ell": (-0.15, 0.20)}
BOUNDS = [BOX[k] for k in FREE]
LO, HI = np.array([b[0] for b in BOUNDS]), np.array([b[1] for b in BOUNDS])


def in_box(theta):
    return bool(np.all(theta >= LO - 1e-12) and np.all(theta <= HI + 1e-12))


def params_cell(theta, cell):
    p = V48.params(theta)
    p.update(CELLS[cell])
    return p


# ================= one evaluation =================
def score_point(theta, cell, warm=None, source="de"):
    theta = np.asarray(theta, float)
    p = params_cell(theta, cell)
    t0 = time.time()
    row = {"theta": theta.tolist(), "cell": cell, "source": source, "ok": False, "merit": FAIL_MERIT, "fail": None}
    try:
        if warm is None:
            o = V3.evaluate(p, N=NG, X=XG)
        else:
            o = V3.evaluate(p, F0=(warm[0], warm[1]), z0=(warm[0], warm[2]), N=NG, X=XG)
    except Exception as ex:                                            # a solve that raises (singular Hessian, overflow)
        o, row["fail"] = None, repr(ex)[:120]
    row["wall"] = time.time() - t0
    if o is None or not np.isfinite(o["M"]) or not np.all(np.isfinite(o["errors"])):
        row["fail"] = row["fail"] or "solver returned None or a non-finite value"
        return row, None
    gatev = VR.gate_dev(p); gate = float(np.max(np.abs(gatev)))
    x, z = o["profile"]["x"], o["profile"]["z"]
    dz2 = np.diff(z ** 2) / np.diff(x)
    kg = float(0.5 - V48.E2 * p["ell"] ** 2 * np.max(dz2 ** 2) / (4 * np.pi ** 4 * p["c"]))
    proxy = float(np.pi ** 3 * p["c"] * 783.0 ** 2 / p["f"] ** 2)
    e15 = float(o["E15"])
    merit = e15 + RHO * max(0.0, W - proxy) + PEN_GATE * max(0.0, gate - GATE_ACC) + PEN_KG * max(0.0, KG_MIN - kg)
    row.update(ok=True, merit=float(merit), E15=e15, errors=[float(v) for v in o["errors"]], gate=gate, Kgamma=kg, proxy=proxy,
               feasible=bool(gate <= GATE_ACC and kg >= KG_MIN and p["eta"] <= ETA_MAX + 1e-12 and proxy >= W),
               M=float(o["M"]), z0=float(o["z0"]), Delta=float(o["Delta"]), gA=float(o["gA"]), rEV2=float(o["rEV2"]), rn2=float(o["rn2"]),
               rMV2=float(o["rMV2"]), E_tr=float(o["E_tr"]), eta=float(p["eta"]), nit=int(o.get("nit", -1)), converged=bool(o.get("converged", True)))
    return row, o


_HIST = None


def _hist_file(cell):
    global _HIST
    if _HIST is None or not _HIST.name.endswith(f"hist_{cell}_{os.getpid()}.jsonl"):
        os.makedirs(OUT_DIR, exist_ok=True)
        _HIST = open(os.path.join(OUT_DIR, f"hist_{cell}_{os.getpid()}.jsonl"), "a")
    return _HIST


def de_objective(theta, cell):
    row, _ = score_point(theta, cell, source="de")
    f = _hist_file(cell); f.write(json.dumps(row) + "\n"); f.flush()
    return row["merit"]


# ================= R30-0, the controls =================
def sha(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def fine_eval(theta, cell):
    """V54b's continuation: 500 / 12 -> 800 / 18 -> 1200 / 24, cold at the first stage; returns the final evaluate() dict."""
    p = params_cell(theta, cell)
    prev, o = None, None
    for N, X in ((500, 12.0), (800, 18.0), (1200, 24.0)):
        o = V3.evaluate(p, N=N, X=X) if prev is None else V3.evaluate(p, F0=(prev["x"], prev["F"]), z0=(prev["x"], prev["z"]), N=N, X=X)
        if o is None:
            return None
        prev = o["profile"]
    return o


def lhs(n, seed):
    from scipy.stats import qmc
    u = qmc.LatinHypercube(d=len(FREE), seed=seed).random(n)
    return LO + u * (HI - LO)


def _ctrl_point(args):
    theta, cell = args
    row, _ = score_point(theta, cell, source="lhs")
    return row


def controls(workers):
    out = {"kit": "a scratch copy of LdGS_vacuum_sectors outside the repo (hashed against the theory/ copy, check a)", "cut_lines_of_exec_parts": CUTS, "coarse_grid": {"N": NG, "X": XG}, "box_theta": {k: list(v) for k, v in BOX.items()},
           "cells": CELLS, "pop_per_dim": POP_PER_DIM, "caps": {"main": CAP_MAIN, "variant": CAP_VARIANT},
           "convergence": {"spread": CONV_SPREAD, "generations": CONV_GENS, "fraction": CONV_FRACTION}, "seed": SEED}
    # (a) the hashes
    hashes = {}
    for f in KIT_FILES:
        a, b = sha(os.path.join(KIT, f)), sha(os.path.join(THEORY, f))
        hashes[f] = {"sha256": a, "same_as_theory_copy": a == b}
    out["a_hashes"] = hashes; out["a_pass"] = all(h["same_as_theory_copy"] for h in hashes.values())
    log(f"(a) hashes: {'all equal' if out['a_pass'] else 'DIFFER'}")
    # the box holds the record states
    states = {"V48_shipped": SHIPPED["params"], "V52_n4": V52["n4"]["params"], "V52_lock": V52["lock"]["params"]}
    for nm, P in (("C9", V3.C9), ("O10", V3.O10)):
        states[nm] = {k: P[k] for k in LOGP} | {"ell": P["ell"], "c_total": P["c_total"]}
    out["box_holds"] = {k: in_box(theta_of(v)) for k, v in states.items()}
    # (b) the shipped V48 state, coarse and fine
    t0 = time.time(); row, o = score_point(THETA_V48, MAIN_CELL, source="ctrl")
    t_cold = time.time() - t0
    t0 = time.time(); row_w, _ = score_point(THETA_V48, MAIN_CELL, warm=(o["profile"]["x"], o["profile"]["F"], o["profile"]["z"]), source="ctrl")
    t_warm = time.time() - t0
    fo = fine_eval(THETA_V48, MAIN_CELL)
    out["b_coarse"] = {k: row[k] for k in ("E15", "proxy", "gate", "Kgamma", "M", "Delta", "gA", "rEV2", "rn2", "rMV2", "z0", "feasible")}
    out["b_fine"] = {"E15": fo["E15"], "gA": fo["gA"], "rMV2": fo["rMV2"], "rEV2": fo["rEV2"], "rn2": fo["rn2"], "M": fo["M"], "Delta": fo["Delta"]}
    out["b_expected"] = {"coarse_E15": SHIPPED["E15"], "coarse_proxy": SHIPPED["proxy"], "fine_E15": 9.1997, "fine_gA": 1.1583, "fine_rMV2": 0.471}
    out["b_pass"] = bool(abs(row["E15"] / SHIPPED["E15"] - 1) < 1e-4 and abs(row["proxy"] / SHIPPED["proxy"] - 1) < 1e-4
                         and abs(fo["E15"] / 9.1997 - 1) < 1e-4 and abs(fo["gA"] / 1.1583 - 1) < 1e-4 and abs(fo["rMV2"] / 0.471 - 1) < 2e-3)
    log(f"(b) coarse E15 {row['E15']:.4f} proxy {row['proxy']:.3f}; fine E15 {fo['E15']:.4f} gA {fo['gA']:.4f} rMV2 {fo['rMV2']:.4f}; pass {out['b_pass']}")
    # (c) the V55 prices at the shipped state
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        pr = V55["prices"](SHIPPED["params"], {}, "V48")
    lam = pr["0.01"]["lam"]
    top = sorted(lam.items(), key=lambda kv: -kv[1])[:4]
    out["c_prices"] = {"trust_0.01": pr["0.01"], "n_active": pr["n_active"], "rank_active": pr["rank_active"], "top": top, "stdout": buf.getvalue()[-1500:]}
    out["c_expected"] = {"Delta": 0.374, "gA": 0.285, "n_active": 9, "rank_active": 9}
    out["c_pass"] = bool(abs(lam["Delta"] - 0.374) < 0.01 and abs(lam["gA"] - 0.285) < 0.01 and pr["n_active"] == 9 and pr["rank_active"] == 9)
    log(f"(c) prices Delta {lam['Delta']:.3f} gA {lam['gA']:.3f} active {pr['n_active']} rank {pr['rank_active']}; pass {out['c_pass']}")
    # (d) timing and the failure rate on 200 LHS points
    pts = lhs(200, SEED + 7)
    ctx = mp.get_context("spawn")
    t0 = time.time()
    with ctx.Pool(workers) as pool:
        rows = pool.map(_ctrl_point, [(th, MAIN_CELL) for th in pts])
    wall = time.time() - t0
    ok = [r for r in rows if r["ok"]]
    out["d_timing"] = {"cold_s_V48_state": t_cold, "warm_s_V48_state": t_warm, "lhs_points": 200, "lhs_wall_s_pool": wall, "workers": workers,
                       "lhs_mean_wall_s_per_point": float(np.mean([r["wall"] for r in rows])), "lhs_failures": 200 - len(ok), "lhs_failure_rate": (200 - len(ok)) / 200,
                       "lhs_feasible": int(sum(r.get("feasible", False) for r in ok)), "lhs_best_E15": float(min(r["E15"] for r in ok)) if ok else None,
                       "lhs_best_merit": float(min(r["merit"] for r in ok)) if ok else None,
                       "fail_reasons": sorted({r["fail"] for r in rows if not r["ok"]})[:10]}
    log(f"(d) cold {t_cold:.2f} s, warm {t_warm:.2f} s; 200 LHS points in {wall:.1f} s on the pool, failures {200 - len(ok)}, feasible {out['d_timing']['lhs_feasible']}")
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "hist_lhs_controls.jsonl"), "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
    out["wall_s"] = time.time() - T0
    prev = json.load(open(OUT_JSON)) if os.path.exists(OUT_JSON) else {}
    prev["controls"] = out
    json.dump(prev, open(OUT_JSON, "w"), indent=1, default=float)
    return out


# ================= R30-1 (a), the differential evolution =================
class DEState:
    def __init__(self, cell, cap):
        self.cell, self.cap = cell, cap
        self.gen, self.streak, self.t0 = 0, 0, time.time()
        self.path = os.path.join(OUT_DIR, f"de_{cell}_population.npz")
        self.logf = open(os.path.join(OUT_DIR, f"de_{cell}.log"), "a")
        try:
            sys.path.insert(0, os.path.expanduser("~/Library/Application Support/UPSSentinel")); import ups_guard; ups_guard.register(); self.guard = ups_guard
        except Exception:
            self.guard = None
        self.stop_reason = None

    def __call__(self, intermediate_result):
        self.gen += 1
        E = np.sort(intermediate_result.population_energies)
        k = max(1, int(CONV_FRACTION * len(E)) - 1)
        spread = float(E[k] - E[0])
        self.streak = self.streak + 1 if spread < CONV_SPREAD else 0
        np.savez(self.path, population=intermediate_result.population, energies=intermediate_result.population_energies, generation=self.gen,
                 best_x=intermediate_result.x, best_fun=intermediate_result.fun, streak=self.streak)
        line = f"{self.cell} gen {self.gen} best {E[0]:.4f} p90 {E[k]:.4f} worst {E[-1]:.1f} spread {spread:.4f} streak {self.streak} n_fail {int(np.sum(E >= FAIL_MERIT))} t {time.time() - self.t0:.0f}s"
        self.logf.write(line + "\n"); self.logf.flush(); log(line)
        if self.streak >= CONV_GENS:
            self.stop_reason = f"converged: spread < {CONV_SPREAD} for {CONV_GENS} generations"; return True
        if self.guard is not None and self.guard.wrap_up():
            self.stop_reason = "UPS wrap-up flag"; return True
        return False


def run_de(cell, workers, cap):
    from scipy.optimize import differential_evolution
    os.makedirs(OUT_DIR, exist_ok=True)
    st = DEState(cell, cap)
    npop = POP_PER_DIM * len(FREE)
    resumed = False
    if os.path.exists(st.path):
        z = np.load(st.path); init = z["population"]; st.gen = int(z["generation"]); st.streak = int(z["streak"]); resumed = True
        log(f"{cell}: resuming from generation {st.gen} (population {init.shape})")
    else:
        init = lhs(npop, SEED + {"n3g2": 1, "n4g2": 2, "n3g4": 3, "n4g4": 4}[cell])
        init[0] = THETA_V48
        if cell == "n4g2":
            init[1] = theta_of(V52["n4"]["params"])
        if cell == "n3g4":
            init[1] = theta_of(V52["lock"]["params"])
    remaining = max(1, cap - st.gen)
    t0 = time.time()
    res = differential_evolution(de_objective, BOUNDS, args=(cell,), strategy="best1bin", maxiter=remaining, popsize=POP_PER_DIM, tol=0.0, atol=0.0,
                                 mutation=(0.5, 1.0), recombination=0.7, seed=SEED, callback=st, polish=False, init=init, updating="deferred", workers=workers)
    rec = {"cell": cell, "options": CELLS[cell], "workers": workers, "cap": cap, "generations": st.gen, "resumed": resumed, "nfev_this_run": int(res.nfev),
           "wall_s_this_run": time.time() - t0, "stop_reason": st.stop_reason or str(res.message), "best_merit": float(res.fun), "best_theta": [float(v) for v in res.x],
           "converged": st.stop_reason is not None and st.stop_reason.startswith("converged"), "population_spread_final": None}
    E = np.sort(res.population_energies); k = max(1, int(CONV_FRACTION * len(E)) - 1)
    rec["population_spread_final"] = float(E[k] - E[0]); rec["population_best_10"] = [float(v) for v in E[:10]]
    brow, _ = score_point(res.x, cell, source="de-best")
    rec["best_row"] = brow
    prev = json.load(open(OUT_JSON)) if os.path.exists(OUT_JSON) else {}
    prev.setdefault("de", {})[cell] = rec
    json.dump(prev, open(OUT_JSON, "w"), indent=1, default=float)
    log(f"{cell}: DE done, {rec['stop_reason']}, best merit {res.fun:.4f} (E15 {brow.get('E15')}), {res.nfev} evaluations this run, {rec['wall_s_this_run']:.0f} s")
    return rec


# ================= R30-1 (b), the multistart sequential LP =================
def _slp_worker(args):
    cell, name, theta0, ladder = args
    theta0 = np.asarray(theta0, float)
    orig = V48.params

    def patched(th, _c=cell, _o=orig):                            # the cell enters through par_v3's option keys; slp() itself unchanged
        pp = _o(th); pp.update(CELLS[_c]); return pp
    V48.params = patched
    buf = open(os.path.join(OUT_DIR, f"slp_{cell}_{name}.log"), "w")
    t0 = time.time(); rec = {"cell": cell, "start": name, "theta0": theta0.tolist(), "ladder": ladder, "ok": False}
    try:
        with contextlib.redirect_stdout(buf):
            ev = V48.Ev()
            theta = theta0.copy()
            if ladder:
                for Wl in (3.0, 4.5, 5.8):
                    theta, o = V48.slp(ev, theta, Wl, 45 if Wl == 5.8 else 20, f"{cell} {name} W={Wl}")
            else:
                theta, o = V48.slp(ev, theta, W, 45, f"{cell} {name} W={W}")
        row, _ = score_point(theta, cell, source="slp-end")           # the endpoint re-scored cold, comparable with the history
        rec.update(ok=True, theta=theta.tolist(), E15_slp=float(np.max(np.abs(o["errors"]))), proxy_slp=float(V48.proxy(theta)), endpoint=row)
    except Exception as ex:
        rec["error"] = repr(ex)[:200]
    finally:
        V48.params = orig; buf.close()
    rec["wall_s"] = time.time() - t0
    return rec


def run_slp(workers):
    prev = json.load(open(OUT_JSON))
    jobs = []
    for cell, d in prev.get("de", {}).items():
        z = np.load(os.path.join(OUT_DIR, f"de_{cell}_population.npz"))
        order = np.argsort(z["energies"])
        nbest = 6 if cell == MAIN_CELL else 3
        for i, k in enumerate(order[:nbest]):
            jobs.append((cell, f"de{i}", z["population"][k], False))
        if cell == MAIN_CELL:
            for i, th in enumerate(lhs(18, SEED + 99)):
                jobs.append((cell, f"lhs{i}", th, True))
    jobs.append((MAIN_CELL, "v48state", THETA_V48, False))          # the known-answer control: the SLP from the shipped state stays at 9.20
    log(f"slp: {len(jobs)} starts on {workers} workers")
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers) as pool:
        recs = pool.map(_slp_worker, jobs)
    with open(os.path.join(OUT_DIR, "hist_slp_endpoints.jsonl"), "w") as f:
        for r in recs:
            f.write(json.dumps(r, default=float) + "\n")
    ok = [r for r in recs if r["ok"] and r["endpoint"]["ok"]]
    summ = {"starts": len(jobs), "ok": len(ok), "wall_s_max": max(r["wall_s"] for r in recs),
            "endpoints": [{"cell": r["cell"], "start": r["start"], "E15": r["endpoint"]["E15"], "merit": r["endpoint"]["merit"], "proxy": r["endpoint"]["proxy"],
                           "feasible": r["endpoint"]["feasible"], "M": r["endpoint"]["M"], "z0": r["endpoint"]["z0"]} for r in ok],
            "errors": [{"cell": r["cell"], "start": r["start"], "error": r.get("error")} for r in recs if not r["ok"]]}
    # the plateau statement: distinct endpoints at 1e-3 in the log-coefficients (feasible, within 1 point of the best)
    feas = [r for r in ok if r["endpoint"]["feasible"]]
    if feas:
        best = min(r["endpoint"]["E15"] for r in feas)
        near = [np.array(r["theta"]) for r in feas if r["endpoint"]["E15"] < best + 1.0]
        distinct = []
        for th in near:
            if not any(np.max(np.abs(th - d)) < 1e-3 for d in distinct):
                distinct.append(th)
        summ["plateau"] = {"best_E15": best, "within_1_point": len(near), "distinct_at_1e-3": len(distinct),
                           "E15_spread_within_1_point": [float(r["endpoint"]["E15"]) for r in feas if r["endpoint"]["E15"] < best + 1.0]}
    prev["slp"] = summ
    json.dump(prev, open(OUT_JSON, "w"), indent=1, default=float)
    log(f"slp done: {len(ok)} / {len(jobs)} endpoints; " + (f"best feasible E15 {summ['plateau']['best_E15']:.4f}, {summ['plateau']['distinct_at_1e-3']} distinct within 1 point" if feas else "no feasible endpoint"))
    return summ


# ================= the history =================
def consolidate():
    rows = []
    for f in sorted(glob.glob(os.path.join(OUT_DIR, "hist_*.jsonl"))):
        if f.endswith("hist_slp_endpoints.jsonl"):
            continue
        for line in open(f):
            rows.append(json.loads(line))
    n = len(rows); d = len(FREE)
    theta = np.full((n, d), np.nan); err = np.full((n, 15), np.nan)
    cols = {k: np.full(n, np.nan) for k in ("merit", "E15", "gate", "Kgamma", "proxy", "M", "z0", "Delta", "gA", "rEV2", "rn2", "rMV2", "E_tr", "eta", "wall")}
    cell = np.array([r["cell"] for r in rows]); src = np.array([r["source"] for r in rows])
    ok = np.array([r["ok"] for r in rows]); feas = np.array([bool(r.get("feasible", False)) for r in rows])
    for i, r in enumerate(rows):
        theta[i] = r["theta"]
        for k in cols:
            if k in r:
                cols[k][i] = r[k]
        if r["ok"]:
            err[i] = r["errors"]
    os.makedirs(OUT_DIR, exist_ok=True)
    np.savez_compressed(os.path.join(OUT_DIR, "history.npz"), theta=theta, errors=err, cell=cell, source=src, ok=ok, feasible=feas, **cols)
    summ = {"evaluations": n, "solved": int(ok.sum()), "feasible": int(feas.sum()), "per_cell": {}}
    for c in CELLS:
        m = cell == c
        if m.any():
            mo = m & ok; mf = m & feas
            summ["per_cell"][c] = {"evaluations": int(m.sum()), "solved": int(mo.sum()), "feasible": int(mf.sum()),
                                   "best_feasible_E15_coarse": float(np.nanmin(cols["E15"][mf])) if mf.any() else None,
                                   "best_merit": float(np.nanmin(cols["merit"][mo])) if mo.any() else None,
                                   "mean_wall_s": float(np.nanmean(cols["wall"][m]))}
    prev = json.load(open(OUT_JSON)) if os.path.exists(OUT_JSON) else {}
    prev["history"] = summ
    json.dump(prev, open(OUT_JSON, "w"), indent=1, default=float)
    log(f"history: {n} evaluations, {int(ok.sum())} solved, {int(feas.sum())} feasible")
    return summ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode"); ap.add_argument("cell", nargs="?")
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--maxiter", type=int, default=None)
    a = ap.parse_args()
    if a.mode == "controls":
        controls(a.workers)
    elif a.mode == "de":
        cap = a.maxiter or (CAP_MAIN if a.cell == MAIN_CELL else CAP_VARIANT)
        run_de(a.cell, a.workers, cap)
    elif a.mode == "slp":
        run_slp(a.workers)
    elif a.mode == "consolidate":
        consolidate()
    else:
        raise SystemExit(f"unknown mode {a.mode}")


if __name__ == "__main__":
    main()
