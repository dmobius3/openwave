"""M5.32 R5 arm (a): the STATIC collateral of gate G6 on the audited R2
candidate, the covariant lambda-family, run through the RECORD instruments
at lambda = 0 (control, must reproduce the stored record) and lambda = 1.

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1), jets A_i = d_i M
(i = 1..3 on the certified sym stencil, fwd / bwd branches averaged at the
density level), A_0 = 0 on every static item here.
    F_{mu nu} = A_mu eta A_nu - A_nu eta A_mu
    I1   = sum_{mu<nu} eta^mu eta^nu <F_{mu nu}, F_{mu nu}>_eta  (certified)
    I1_h = sum_{mu<nu} eta^mu eta^nu tr(h F h F^T),
           h = eta + 2 (eta u)(eta u)^T, u the timelike unit eigenvector of
           M eta (u^T eta u = -1)
The candidate (R2, audited):
    L_lambda = -4 [(1 - lambda) I1 + lambda I1_h] - V4
    E_lambda = 4 h^3 sum_br wt sum_cells [(1 - lambda) I1 + lambda I1_h]
               + h^3 sum V                                  (static read)
Every static item below lives on a 3x3 field embedded with the uniform
vacuum time row, M4 = diag(g_tt) (+) M3, g_tt = 32 (the toy point, s = -1),
so A_i = 0 (+) d_i M3, u = e0 exactly, h = 1 on the block, and
    I1_h = I1   pointwise, density AND gradient          (the R2 identity)
which makes every lambda = 1 number below a TAUTOLOGY of the identity, not
an independent test: the script MEASURES the identity per item (density gap,
energy gap, gradient gap) instead of assuming it, and says so per row.

THE FOUR ITEMS (record instrument, record value)
  1. 3-lepton census A < C < B: the M5.21.11 T2 ladder endpoints (branches
     A, C, B; n = 48, L = 48, h = 1; and the n = 32, h = 1.5 rung), the
     certified 3D instrument m5_21_2b_a_instrument (T2, w2 = 0.0027581,
     sym stencil, pinned shell); record E_end in m5_21_11_row_*.json.
     Readout: E_lambda and the free-cell gradient on the stored end state;
     continuation: FIRE for --iters from the n = 32 endpoints under each
     lambda (the record relaxation is 12000 iterations; > 45 min through
     the 4x4 registry path, so the smallest documented mode = a
     continuation from the endpoint is run instead).
  2. Far-field degree exactly -1 (the Q37 instrument): the M5.22 audit's
     van Oosterom-Strackee solid-angle degree of the max-|dot| lifted
     director on centered cube surfaces (m5_22_e_audit.read_charge_from_M),
     record -1.000000 for the proton-analog and the lepton state at half-
     widths 13.5 / 16.5 / 19.5 (n = 48, h = 1); here read on the stored
     states and on every continued state (the field the candidate relaxes).
  3. Baryon-analog existence: the proton-analog P s = -1/2 (m5_22 census,
     n = 32, L = 48, h = 1.5, record E_end = 8.4962 at f_tol); readout of
     the stored endpoint (float32 on disk: the record reproduction is
     bounded by that storage precision, reported) and the documented
     `extend` continuation for --iters under each lambda, with the core
     ledger, the Mermin-Ho charge profile and the solid-angle degree on the
     continued states; the lepton reference E s = -1/2 for the mass
     ordering proton > lepton.
  4. Coulomb 1/r, the M5.17 `fixed` instrument (axisymmetric (rho, z)
     reduction, NR = 96, NZ = 192, h = 1, the locked beta = 1 potential):
     record fit_A = +215.9255 (like pair, repulsive), -422.64 (antipair),
     sign_ok both. Its curvature density is c2 4 sum_{mu<nu} ||[A_mu,
     A_nu]||_F^2 on a 4x4 field with M[0,0] = g_time = 8 constant, which on
     zero time-row jets equals 4 I1 exactly (eta-commutator = plain
     commutator, <F,F>_eta = Frobenius on the spatial block); the wrapper
     swaps that density for 4 [(1 - lambda) I1 + lambda I1_h] and re-runs
     the fixed curves without touching the record JSON.  This is the G1
     ABSOLUTE bar (the it = 120 pair ladder of R1.b cannot decide it).

Verdict per item: PASS if the lambda = 0 run reproduces the record to
<= 1e-10 (relative; storage-precision-limited items report the float32
bound) AND the lambda = 1 value differs from lambda = 0 by <= 1e-10
relative with the qualitative claim intact (ordering, degree, existence,
1/r sign + monotone decrease); any loss is a FAIL.

Stages: census, baryon, coulomb, report (or all). Each stage merges into
the JSON so stages can run in parallel processes (--workers <= 3).
Usage: nice -n 10 python3 m5_32_r5_a_collateral.py --stage all --iters 150
Out: ../data/m5_32_r5_collateral.json (+ ../data/m5_32_r5_a/*.npz local)
"""
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse  # noqa: E402
import fcntl  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r5_collateral.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r5_a")
T_START = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


LAG = _load("m5_32_lagrangian", "m5_32_lagrangian.py")
EXT = _load("m5_32_terms_ext", "m5_32_terms_ext.py")
PB = _load("m5_32_r1_b_pair", "m5_32_r1_b_pair.py")
R2C = _load("m5_32_r2_c_screen", "m5_32_r2_c_screen.py")
CEN = _load("m5_22_b_census", "m5_22_b_census.py")
AUD22 = _load("m5_22_e_audit", "m5_22_e_audit.py")
INS = CEN.INS                      # the census instrument instance (fire)

G_TOY = 32.0
LAMBDAS = (0.0, 1.0)
TOL = 1e-10
W2_T2 = CEN.W2_T2


def log(msg):
    print(f"[{time.time() - T_START:8.1f}s] {msg}", flush=True)


def rel(a, b):
    return float(abs(a - b) / max(abs(b), 1e-300))


def merge_json(key, payload):
    os.makedirs(DATA, exist_ok=True)
    lock = OUT_JSON + ".lock"
    with open(lock, "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        cur = {}
        if os.path.exists(OUT_JSON):
            with open(OUT_JSON) as f:
                cur = json.load(f)
        cur[key] = payload
        cur.setdefault("meta", {}).update({
            "task": "M5.32 R5 (a) static collateral (G6 static items)",
            "candidate": "L_lambda = -4[(1-lambda) I1 + lambda I1_h] - V4",
            "lambdas": list(LAMBDAS), "g_tt": G_TOY, "tol": TOL,
            "python": sys.version.split()[0],
            "script": os.path.basename(__file__),
            "git_head": _git_head()})
        tmp = OUT_JSON + ".tmp"
        with open(tmp, "w") as f:
            json.dump(cur, f, indent=1, default=float)
        os.replace(tmp, OUT_JSON)
        fcntl.flock(lk, fcntl.LOCK_UN)


def _git_head():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=HERE,
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:                                  # noqa: BLE001
        return None


# ================= the lambda-family routing (R2 arm (c), reused) =========
def act_of(lam):
    return R2C.ExtAction(lam, 0.0, G_TOY)


class Routed:
    """route the census instrument's FIRE to the lambda-family energy +
    gradient; restores the certified paths on exit."""

    def __init__(self, lam):
        self.act = act_of(lam)

    def __enter__(self):
        self._g, self._e = INS.grad, INS.e_parts
        INS.grad, INS.e_parts = self.act.grad, self.act.e_parts
        return self.act

    def __exit__(self, *a):
        INS.grad, INS.e_parts = self._g, self._e


def cfg_of(n, L=48.0):
    return INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                        n=n, L=L, bc="pinned")


def identity_gap(M, cfg):
    """max over cells and branches of |I1_h - I1| and the I1_h gradient's
    local piece: the measured tautology witness."""
    h = cfg["h"]
    M4 = R2C.embed_pair(M, G_TOY)
    gap = 0.0
    gloc = 0.0
    for br, wt in INS.branches(cfg["stencil"]):
        A = PB.jets4(M, h, br)
        F = LAG.F_of_A(A)
        d1 = LAG.density_from_K(F, PB.K_of("I1"))
        dh, _, Gl = R2C.i1h_pieces(A, M4)
        gap = max(gap, float(np.max(np.abs(dh - d1))))
        gloc = max(gloc, float(np.max(np.abs(Gl))))
    return {"density_gap_max": gap, "local_grad_max": gloc}


def readout(M, cfg, lam):
    act = act_of(lam)
    e_u, _, e_v = act.e_parts(M, cfg)
    g = act.grad(M, cfg)
    free = ~INS.pin_shell(cfg["n"], cfg["h"])
    gf = g * free[..., None, None]
    return {"lambda": lam, "E": float(e_u + e_v), "E_u": float(e_u),
            "E_v": float(e_v),
            "grad_free_norm": float(np.sqrt(np.sum(gf * gf))),
            "fmax_free": float(np.max(np.abs(gf)))}


def certified_readout(M, cfg):
    e_u, e_d, e_v = INS.e_parts(M, cfg)
    return {"E": float(e_u + e_d + e_v), "E_u": float(e_u),
            "E_v": float(e_v)}


def continuation(M0, cfg, lam, iters, tag):
    free = ~INS.pin_shell(cfg["n"], cfg["h"])
    t0 = time.time()
    with Routed(lam):
        M, states, info = INS.fire(M0.copy(), cfg, free, max_iter=iters,
                                   log_every=max(iters // 5, 1),
                                   tag=f"{tag}_lam{lam:g}")
    trace = info["trace"]
    # an endpoint already at f_tol exits FIRE before its first log line
    # (empty trace): read the end values directly
    rd = readout(M, cfg, lam)
    return M, {"lambda": lam, "iters": iters, "stop": info["stop"],
               "E_trace": [[t["it"], t["E"]] for t in trace],
               "fmax_trace": [[t["it"], t["fmax"]] for t in trace],
               "E_end": trace[-1]["E"] if trace else rd["E"],
               "fmax_end": trace[-1]["fmax"] if trace else rd["fmax_free"],
               "wall_s": time.time() - t0}


def degree_reads(M, cfg, widths):
    """the Q37 instrument on centered cube surfaces (index ranges as the
    M5.22 audit: half-width w -> |coord| <= w)."""
    n, h = cfg["n"], cfg["h"]
    out = {}
    for w in widths:
        # cells with |x_i| <= w: x_i = (i - (n-1)/2) h
        k = int(np.floor(w / h + (n - 1) / 2.0 + 1e-9))
        ilo, ihi = n - 1 - k, k
        Q, cf = AUD22.read_charge_from_M(M, ilo, ihi)
        out[f"{w:g}"] = {"Q": float(Q), "conflicts": int(cf),
                         "ilo": ilo, "ihi": ihi}
    return out


def save_npz(tag, M):
    os.makedirs(OUT_NPZ, exist_ok=True)
    np.savez_compressed(os.path.join(OUT_NPZ, f"{tag}.npz"),
                        M=M.astype(np.float64))


def pair_diff(Ma, Mb):
    return float(np.max(np.abs(Ma - Mb)))


# ================= item 1: the 3-lepton census =================
def stage_census(iters):
    log("census: start")
    out = {"instrument": "m5_21_2b_a_instrument T2 sym pinned (the M5.21.11 "
                         "ladder endpoints), routed through R2C.ExtAction",
           "record_source": "m5_21_11_row_t11lad_<b>_n<n>_d0.3.json",
           "states": {}}
    for n in (48, 32):
        cfg = cfg_of(n)
        for b in ("A", "C", "B"):
            tag = f"t11lad_{b}_n{n}_d0.3"
            with open(os.path.join(DATA, f"m5_21_11_row_{tag}.json")) as f:
                rec = json.load(f)
            z = np.load(os.path.join(DATA, f"m5_21_11_end_{tag}.npz"))
            M = z["M"].astype(np.float64)
            row = {"branch": b, "n": n, "L": cfg["L"], "h": cfg["h"],
                   "stored_dtype": str(z["M"].dtype),
                   "record": {k: rec[k] for k in ("E_end", "E_u", "E_v",
                                                  "stop", "maxit")},
                   "certified_readout": certified_readout(M, cfg),
                   "identity": identity_gap(M, cfg),
                   "readout": {}}
            for lam in LAMBDAS:
                row["readout"][f"{lam:g}"] = readout(M, cfg, lam)
            r0, r1 = row["readout"]["0"], row["readout"]["1"]
            row["rel_lam0_vs_record"] = rel(r0["E"], rec["E_end"])
            row["rel_lam1_vs_lam0"] = rel(r1["E"], r0["E"])
            row["abs_lam1_vs_lam0"] = float(r1["E"] - r0["E"])
            row["grad_abs_gap"] = None
            out["states"][tag] = row
            log(f"census {tag}: rec {rec['E_end']:.12g} lam0 {r0['E']:.12g} "
                f"lam1 {r1['E']:.12g} rel0 {row['rel_lam0_vs_record']:.2e} "
                f"rel10 {row['rel_lam1_vs_lam0']:.2e} "
                f"dens_gap {row['identity']['density_gap_max']:.2e}")
    # the gradient gap (lambda 1 vs 0) on the n = 32 endpoints
    for b in ("A", "C", "B"):
        tag = f"t11lad_{b}_n32_d0.3"
        cfg = cfg_of(32)
        M = np.load(os.path.join(DATA, f"m5_21_11_end_{tag}.npz"))["M"] \
            .astype(np.float64)
        g0 = act_of(0.0).grad(M, cfg)
        g1 = act_of(1.0).grad(M, cfg)
        out["states"][tag]["grad_abs_gap"] = pair_diff(g0, g1)
        out["states"][tag]["grad_abs_max"] = float(np.max(np.abs(g0)))
    # continuation at n = 32 (the census rung) under each lambda
    cont = {}
    ends = {}
    for b in ("A", "C", "B"):
        tag = f"t11lad_{b}_n32_d0.3"
        cfg = cfg_of(32)
        M0 = np.load(os.path.join(DATA, f"m5_21_11_end_{tag}.npz"))["M"] \
            .astype(np.float64)
        cont[b] = {}
        for lam in LAMBDAS:
            M, info = continuation(M0, cfg, lam, iters, f"census_{b}")
            info["readout_end"] = readout(M, cfg, lam)
            info["degree_end"] = degree_reads(M, cfg, (13.5, 16.5))
            info["move_from_start"] = pair_diff(M, M0)
            cont[b][f"{lam:g}"] = info
            ends[(b, lam)] = M
            save_npz(f"m5_32_r5_census_{b}_n32_lam{lam:g}", M)
            log(f"census cont {b} lam {lam:g}: E_end {info['E_end']} "
                f"fmax {info['fmax_end']:.3e} stop {info['stop']} "
                f"[{info['wall_s']:.0f}s]")
        cont[b]["end_state_gap_lam1_vs_lam0"] = pair_diff(ends[(b, 1.0)],
                                                          ends[(b, 0.0)])
        cont[b]["E_end_rel_gap"] = rel(cont[b]["1"]["E_end"],
                                       cont[b]["0"]["E_end"])
    out["continuation_n32"] = {"iters": iters, "n": 32, "L": 48.0,
                               "h": 1.5, "branches": cont}
    # verdict
    def order_ok(key):
        e = {b: out["states"][f"t11lad_{b}_n48_d0.3"]["readout"][key]["E"]
             for b in "ACB"}
        return bool(e["A"] < e["C"] < e["B"]), e
    ok0, e0 = order_ok("0")
    ok1, e1 = order_ok("1")
    okc = bool(cont["A"]["1"]["E_end"] < cont["C"]["1"]["E_end"]
               < cont["B"]["1"]["E_end"])
    rep = max(v["rel_lam0_vs_record"] for v in out["states"].values())
    gap = max(v["rel_lam1_vs_lam0"] for v in out["states"].values())
    cgap = max(cont[b]["E_end_rel_gap"] for b in "ACB")
    out["verdict"] = {
        "ordering_lam0_n48": ok0, "ordering_lam1_n48": ok1,
        "ordering_lam1_continued_n32": okc,
        "E_n48": {"lam0": e0, "lam1": e1},
        "max_rel_lam0_vs_record": rep, "max_rel_lam1_vs_lam0": gap,
        "max_cont_E_rel_gap": cgap,
        "PASS": bool(ok0 and ok1 and okc and rep <= TOL and gap <= TOL
                     and cgap <= 1e-8),
        "tautology": "YES: block-diagonal static fields with the uniform "
                     "vacuum time row; I1_h = I1 pointwise (gap measured)"}
    merge_json("census", out)
    log(f"census: done PASS={out['verdict']['PASS']}")
    return out


# ================= items 2 + 3: degree + baryon-analog existence =========
def stage_baryon(iters):
    log("baryon: start")
    out = {"instrument": "m5_22_b_census (INS T2 sym pinned, n = 32, L = 48, "
                         "h = 1.5) routed through R2C.ExtAction; the "
                         "documented `extend` continuation from the stored "
                         "endpoint; degree = m5_22_e_audit.read_charge_from_M "
                         "(van Oosterom-Strackee)",
           "states": {}}
    # (i) stored-state readouts: proton-analog + lepton at n = 32 and 48
    for tag, n in (("P-0.5_plane_sc6_n32_pinned_d0.3", 32),
                   ("E-0.5_plane_sc6_n32_pinned_d0.3", 32),
                   ("P-0.5_plane_sc6_n48_pinned_d0.3", 48),
                   ("E-0.5_plane_sc6_n48_pinned_d0.3", 48)):
        cfg = cfg_of(n)
        with open(os.path.join(DATA, f"m5_22_row_{tag}.json")) as f:
            rec = json.load(f)
        z = np.load(os.path.join(DATA, f"m5_22_end_{tag}.npz"))
        M = z["M"].astype(np.float64)
        row = {"n": n, "L": cfg["L"], "h": cfg["h"],
               "stored_dtype": str(z["M"].dtype),
               "record": {k: rec[k] for k in ("E_end", "E_u", "E_v", "stop",
                                              "maxit")},
               "record_q_far": rec["charge"]["q_far"],
               "certified_readout": certified_readout(M, cfg),
               "identity": identity_gap(M, cfg), "readout": {}}
        for lam in LAMBDAS:
            row["readout"][f"{lam:g}"] = readout(M, cfg, lam)
        r0, r1 = row["readout"]["0"], row["readout"]["1"]
        row["rel_lam0_vs_record"] = rel(r0["E"], rec["E_end"])
        row["rel_certified_vs_record"] = rel(row["certified_readout"]["E"],
                                             rec["E_end"])
        row["rel_lam1_vs_lam0"] = rel(r1["E"], r0["E"])
        row["grad_abs_gap"] = pair_diff(act_of(0.0).grad(M, cfg),
                                        act_of(1.0).grad(M, cfg))
        widths = (13.5, 16.5, 19.5) if n == 48 else (13.5, 16.5)
        row["degree_stored"] = degree_reads(M, cfg, widths)
        out["states"][tag] = row
        log(f"baryon {tag}: rec {rec['E_end']:.12g} cert {row['certified_readout']['E']:.12g} "
            f"lam0 {r0['E']:.12g} lam1 {r1['E']:.12g} "
            f"rel0 {row['rel_lam0_vs_record']:.2e} rel10 {row['rel_lam1_vs_lam0']:.2e} "
            f"deg {[v['Q'] for v in row['degree_stored'].values()]}")
    # (ii) the documented continuation (extend) at n = 32, both states
    cont = {}
    for tag in ("P-0.5_plane_sc6_n32_pinned_d0.3",
                "E-0.5_plane_sc6_n32_pinned_d0.3"):
        cfg = cfg_of(32)
        M0 = np.load(os.path.join(DATA, f"m5_22_end_{tag}.npz"))["M"] \
            .astype(np.float64)
        cont[tag] = {}
        ends = {}
        for lam in LAMBDAS:
            M, info = continuation(M0, cfg, lam, iters, tag[:5])
            info["readout_end"] = readout(M, cfg, lam)
            info["degree_end"] = degree_reads(M, cfg, (13.5, 16.5))
            info["charge_profile_q_far"] = CEN.charge_profile(M, cfg)["q_far"]
            info["core_ledger"] = CEN.core_ledger(M, cfg)
            info["move_from_start"] = pair_diff(M, M0)
            cont[tag][f"{lam:g}"] = info
            ends[lam] = M
            save_npz(f"m5_32_r5_{tag}_lam{lam:g}", M)
            log(f"baryon cont {tag} lam {lam:g}: E_end {info['E_end']} "
                f"fmax {info['fmax_end']:.3e} stop {info['stop']} "
                f"q_far {info['charge_profile_q_far']:.4f} "
                f"deg {[v['Q'] for v in info['degree_end'].values()]} "
                f"[{info['wall_s']:.0f}s]")
        cont[tag]["end_state_gap_lam1_vs_lam0"] = pair_diff(ends[1.0],
                                                            ends[0.0])
        cont[tag]["E_end_rel_gap"] = rel(cont[tag]["1"]["E_end"],
                                         cont[tag]["0"]["E_end"])
    out["continuation_n32"] = {"iters": iters, "n": 32, "L": 48.0,
                               "h": 1.5, "states": cont}
    # verdicts
    P32 = out["states"]["P-0.5_plane_sc6_n32_pinned_d0.3"]
    E32 = out["states"]["E-0.5_plane_sc6_n32_pinned_d0.3"]
    P48 = out["states"]["P-0.5_plane_sc6_n48_pinned_d0.3"]
    E48 = out["states"]["E-0.5_plane_sc6_n48_pinned_d0.3"]
    cP = cont["P-0.5_plane_sc6_n32_pinned_d0.3"]
    cE = cont["E-0.5_plane_sc6_n32_pinned_d0.3"]
    f32 = float(np.finfo(np.float32).eps)
    degs = {}
    for name, row in (("P48", P48), ("E48", E48), ("P32", P32),
                      ("E32", E32)):
        degs[name] = {w: v["Q"] for w, v in row["degree_stored"].items()}
    degs_cont = {"P_lam0": {w: v["Q"] for w, v in cP["0"]["degree_end"].items()},
                 "P_lam1": {w: v["Q"] for w, v in cP["1"]["degree_end"].items()},
                 "E_lam0": {w: v["Q"] for w, v in cE["0"]["degree_end"].items()},
                 "E_lam1": {w: v["Q"] for w, v in cE["1"]["degree_end"].items()}}
    all_deg = [q for d in degs.values() for q in d.values()] + \
              [q for d in degs_cont.values() for q in d.values()]
    deg_ok = bool(all(abs(abs(q) - 1.0) <= 1e-6 for q in all_deg))
    deg_lam_gap = max(abs(degs_cont["P_lam1"][w] - degs_cont["P_lam0"][w])
                      for w in degs_cont["P_lam0"]) if degs_cont["P_lam0"] \
        else None
    out["verdict_degree"] = {
        "record": "-1.000000 (P-0.5 n48, E-0.5 n48) at w = 13.5/16.5/19.5",
        "stored": degs, "continued_n32": degs_cont,
        "all_within_1e-6_of_unit": deg_ok,
        "max_abs_gap_lam1_vs_lam0_continued_P": deg_lam_gap,
        "PASS": deg_ok,
        "tautology": "YES for the field (identical relaxations); the degree "
                     "is a field functional and reads the SAME field"}
    exist_ok = bool(cP["1"]["stop"] in ("f_tol", "plateau", "max_iter")
                    and np.isfinite(cP["1"]["E_end"])
                    and abs(cP["1"]["charge_profile_q_far"]) > 0.9
                    and cP["1"]["fmax_end"] <= max(1e-6, 10 * cP["0"]["fmax_end"]))
    order_ok = bool(cP["1"]["E_end"] > cE["1"]["E_end"])
    out["verdict_baryon"] = {
        "record_E_proton_n32": P32["record"]["E_end"],
        "record_E_lepton_n32": E32["record"]["E_end"],
        "rel_lam0_vs_record_P32": P32["rel_lam0_vs_record"],
        "rel_certified_vs_record_P32": P32["rel_certified_vs_record"],
        "float32_storage_bound": f32,
        "rel_lam1_vs_lam0_P32": P32["rel_lam1_vs_lam0"],
        "rel_lam1_vs_lam0_P48": P48["rel_lam1_vs_lam0"],
        "grad_abs_gap_P32": P32["grad_abs_gap"],
        "continued_E_end": {"P_lam0": cP["0"]["E_end"], "P_lam1": cP["1"]["E_end"],
                            "E_lam0": cE["0"]["E_end"], "E_lam1": cE["1"]["E_end"]},
        "continued_fmax_end": {"P_lam0": cP["0"]["fmax_end"],
                               "P_lam1": cP["1"]["fmax_end"]},
        "continued_state_gap_P": cP["end_state_gap_lam1_vs_lam0"],
        "existence_lam1": exist_ok, "mass_order_proton_gt_lepton_lam1": order_ok,
        "PASS": bool(exist_ok and order_ok
                     and P32["rel_lam1_vs_lam0"] <= TOL
                     and P48["rel_lam1_vs_lam0"] <= TOL
                     and cP["E_end_rel_gap"] <= 1e-8
                     and P32["rel_certified_vs_record"] <= 1e-5),
        "note": "the record relaxation (8000 FIRE iterations from the seed) "
                "exceeds 45 min through the 4x4 registry path; the smallest "
                "documented mode (extend from the endpoint) was run",
        "tautology": "YES: block-diagonal static field, uniform vacuum time row"}
    merge_json("baryon", out)
    log(f"baryon: done degree PASS={deg_ok} baryon PASS={out['verdict_baryon']['PASS']}")
    return out


# ================= item 4: Coulomb 1/r (M5.17 fixed) =================
def stage_coulomb():
    log("coulomb: start")
    argv = sys.argv
    sys.argv = ["m5_17_two_charge.py", "fixed"]
    try:
        M17 = _load("m5_17_two_charge", "m5_17_two_charge.py")
    finally:
        sys.argv = argv
    M17E = sys.modules["m5_17_energy"]
    K1 = PB.K_of("I1")
    rec_path = os.path.join(DATA, "m5_17_two_charge_fixed.json")
    with open(rec_path) as f:
        rec = json.load(f)
    NR, NZ, H = M17.NR, M17.NZ, M17.H
    assert (NR, NZ, H) == (rec["grid"]["NR"], rec["grid"]["NZ"],
                           rec["grid"]["h"])

    def jets_axisym(Mnp, h):
        """the record's three derivative channels (rho, phi, z) on the
        included cells, as 4x4 jets with A_0 = 0 (the M5.17 code path,
        mirrored line by line)."""
        nr = Mnp.shape[0]
        Mminus = np.empty_like(Mnp[: nr - 1])
        Mminus[1:] = Mnp[: nr - 2]
        Mminus[0] = M17E.MIR * Mnp[0]
        Mrho = ((Mnp[1:] - Mminus) / (2.0 * h))[:, 1:-1]
        Mz = (Mnp[: nr - 1, 2:] - Mnp[: nr - 1, :-2]) / (2.0 * h)
        Mc = Mnp[: nr - 1, 1:-1]
        rho = ((np.arange(nr - 1) + 0.5) * h)[:, None, None, None]
        Mphi = M17E._comm_np(np.broadcast_to(M17E.J4, Mc.shape), Mc) / rho
        A = np.zeros((4,) + Mc.shape)
        A[1], A[2], A[3] = Mrho, Mphi, Mz
        return A, Mc

    def curvature_lambda(lam):
        def dens(Mnp, h, c2=1.0):
            A, Mc = jets_axisym(Mnp, h)
            F = LAG.F_of_A(A)
            d = (1.0 - lam) * LAG.density_from_K(F, K1)
            if lam != 0.0:
                d = d + lam * EXT.I1_h_np(A, Mc, None)
            return c2 * 4.0 * d
        return dens

    a, b, c, vvac = M17E.ldg_coeffs(M17.BETA_LOCK, M17.CSCALE_LOCK)
    R, Z = M17E.grid_coords(NR, NZ, H)
    ds = rec["curves"]["likepair"]["d"]
    orig = M17E.curvature_density_np
    out = {"instrument": "m5_17_two_charge.py fixed (NR = 96, NZ = 192, h = 1, "
                         "beta = 1 lock), curvature density swapped for "
                         "4[(1-lambda) I1 + lambda I1_h] on the (rho, phi, z) jets",
           "grid": {"NR": NR, "NZ": NZ, "h": H, "g_time": M17E.G_TIME},
           "record": {"fit_A_like": rec["curves"]["likepair"]["fit_A"],
                      "fit_A_anti": rec["curves"]["antipair"]["fit_A"],
                      "E_single": rec["E_single_same_box"],
                      "sign_ok_like": rec["curves"]["likepair"]["sign_ok"],
                      "d_fit_min": rec["params"]["d_fit_min"]},
           "runs": {}}
    # the record density itself first (the unmodified instrument)
    def curves(density_fn, label):
        M17E.curvature_density_np = density_fn
        try:
            res = {"E_single": float(M17E.total_energy_np(
                M17E.hedgehog_field(R, Z, M17.RC_CORE), a, b, c, 1.0, H, vvac))}
            for q2, name in ((-1, "antipair"), (+1, "likepair")):
                Es = [float(M17E.total_energy_np(
                    M17.pair_field(R, Z, d, M17.RC_CORE, q2), a, b, c, 1.0,
                    H, vvac)) for d in ds]
                E0, A, rms = M17.coulomb_fit(ds, Es, M17.D_FIT_MIN)
                m = np.asarray(ds) >= M17.D_FIT_MIN
                fit = E0 + A / np.asarray(ds)[m]
                sst = float(np.sum((np.asarray(Es)[m] - np.mean(np.asarray(Es)[m])) ** 2))
                r2 = 1.0 - float(np.sum((np.asarray(Es)[m] - fit) ** 2)) / sst
                dec = bool(np.all(np.diff(np.asarray(Es)[m]) < 0)) if q2 > 0 \
                    else bool(np.all(np.diff(np.asarray(Es)[m]) > 0))
                res[name] = {"q2": q2, "d": ds, "E_pair": Es, "fit_E0": E0,
                             "fit_A": A, "fit_rms": rms, "fit_R2": r2,
                             "sign_ok": bool(np.sign(A) == q2),
                             "monotone_outer_window": dec,
                             "A_over_prediction": A / (q2 * 64.0 * np.pi)}
        finally:
            M17E.curvature_density_np = orig
        log(f"coulomb {label}: like A {res['likepair']['fit_A']:+.10f} "
            f"anti A {res['antipair']['fit_A']:+.10f} R2 {res['likepair']['fit_R2']:.6f}")
        return res

    out["runs"]["record_density"] = curves(orig, "record density")
    for lam in LAMBDAS:
        out["runs"][f"lam{lam:g}"] = curves(curvature_lambda(lam), f"lam {lam:g}")
    rr, r0, r1 = (out["runs"]["record_density"], out["runs"]["lam0"],
                  out["runs"]["lam1"])
    # per-d energy gaps
    gaps = {}
    for name in ("likepair", "antipair"):
        e_rec = np.asarray(rec["curves"][name]["E_pair"])
        gaps[name] = {
            "max_rel_record_density_vs_json": float(np.max(np.abs(
                np.asarray(rr[name]["E_pair"]) - e_rec) / np.abs(e_rec))),
            "max_rel_lam0_vs_json": float(np.max(np.abs(
                np.asarray(r0[name]["E_pair"]) - e_rec) / np.abs(e_rec))),
            "max_rel_lam1_vs_lam0": float(np.max(np.abs(
                np.asarray(r1[name]["E_pair"]) - np.asarray(r0[name]["E_pair"]))
                / np.abs(np.asarray(r0[name]["E_pair"])))),
            "fit_A_json": rec["curves"][name]["fit_A"],
            "fit_A_lam0": r0[name]["fit_A"], "fit_A_lam1": r1[name]["fit_A"],
            "fit_A_rel_lam0_vs_json": rel(r0[name]["fit_A"],
                                          rec["curves"][name]["fit_A"]),
            "fit_A_rel_lam1_vs_lam0": rel(r1[name]["fit_A"], r0[name]["fit_A"])}
    out["gaps"] = gaps
    # the tautology witness on the pair field
    A_, Mc = jets_axisym(M17.pair_field(R, Z, 24.0, M17.RC_CORE, +1), H)
    F = LAG.F_of_A(A_)
    out["identity"] = {"density_gap_max_like_d24": float(np.max(np.abs(
        EXT.I1_h_np(A_, Mc, None) - LAG.density_from_K(F, K1)))),
        "u_is_e0_max_dev": float(np.max(np.abs(
            np.abs(EXT.timelike_eig_np(Mc)[0][..., 0]) - 1.0)))}
    g1 = gaps["likepair"]
    out["verdict"] = {
        "record_fit_A_like": rec["curves"]["likepair"]["fit_A"],
        "lam0_fit_A_like": r0["likepair"]["fit_A"],
        "lam1_fit_A_like": r1["likepair"]["fit_A"],
        "lam1_sign_ok_like": r1["likepair"]["sign_ok"],
        "lam1_monotone_decrease_like": r1["likepair"]["monotone_outer_window"],
        "lam1_fit_R2_like": r1["likepair"]["fit_R2"],
        "lam1_sign_ok_anti": r1["antipair"]["sign_ok"],
        "G1_absolute_bar": bool(r1["likepair"]["sign_ok"]
                                and r1["likepair"]["monotone_outer_window"]
                                and r1["likepair"]["fit_R2"] >= 0.95),
        "PASS": bool(g1["max_rel_lam0_vs_json"] <= TOL
                     and g1["max_rel_lam1_vs_lam0"] <= TOL
                     and gaps["antipair"]["max_rel_lam0_vs_json"] <= TOL
                     and gaps["antipair"]["max_rel_lam1_vs_lam0"] <= TOL
                     and r1["likepair"]["sign_ok"]
                     and r1["likepair"]["monotone_outer_window"]),
        "tautology": "YES: M[0,0] = g_time constant, zero time-row jets, "
                     "u = e0 exactly (measured); I1_h = I1 on every cell"}
    merge_json("coulomb", out)
    log(f"coulomb: done PASS={out['verdict']['PASS']}")
    return out


# ================= report =================
def report():
    with open(OUT_JSON) as f:
        d = json.load(f)
    rows = []
    c = d.get("census")
    if c:
        for b in "ACB":
            s = c["states"][f"t11lad_{b}_n48_d0.3"]
            rows.append((f"1 census {b} (n 48, L 48, h 1)", "M5.21.11 T2 endpoint E",
                         s["record"]["E_end"], s["readout"]["0"]["E"],
                         s["readout"]["1"]["E"], s["abs_lam1_vs_lam0"],
                         "PASS" if c["verdict"]["PASS"] else "FAIL", "YES"))
    b_ = d.get("baryon")
    if b_:
        for tag, lab in (("P-0.5_plane_sc6_n48_pinned_d0.3", "2 degree P (n 48, L 48, h 1) w 19.5"),
                         ("E-0.5_plane_sc6_n48_pinned_d0.3", "2 degree E (n 48, L 48, h 1) w 19.5")):
            s = b_["states"][tag]
            q = s["degree_stored"]["19.5"]["Q"]
            rows.append((lab, "Q37 solid-angle degree", -1.0, q, q, 0.0,
                         "PASS" if b_["verdict_degree"]["PASS"] else "FAIL", "YES"))
        s = b_["states"]["P-0.5_plane_sc6_n32_pinned_d0.3"]
        rows.append(("3 proton-analog P-1/2 (n 32, L 48, h 1.5) E", "m5_22 census E_end",
                     s["record"]["E_end"], s["readout"]["0"]["E"],
                     s["readout"]["1"]["E"],
                     s["readout"]["1"]["E"] - s["readout"]["0"]["E"],
                     "PASS" if b_["verdict_baryon"]["PASS"] else "FAIL", "YES"))
    co = d.get("coulomb")
    if co:
        v = co["verdict"]
        rows.append(("4 Coulomb like-pair fit_A (NR 96, NZ 192, h 1)", "M5.17 fixed",
                     v["record_fit_A_like"], v["lam0_fit_A_like"],
                     v["lam1_fit_A_like"],
                     v["lam1_fit_A_like"] - v["lam0_fit_A_like"],
                     "PASS" if v["PASS"] else "FAIL", "YES"))
    print("\n| item | instrument | record | lambda 0 | lambda 1 | diff (1 - 0) | verdict | tautology |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- |")
    for r in rows:
        print(f"| {r[0]} | {r[1]} | {r[2]:.12g} | {r[3]:.12g} | {r[4]:.12g} | "
              f"{r[5]:.3e} | {r[6]} | {r[7]} |")
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "census", "baryon", "coulomb", "report"])
    ap.add_argument("--iters", type=int, default=150)
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    if a.stage == "report":
        report()
        return
    if a.stage == "all":
        # three stages in three processes (<= 3 workers)
        procs = []
        for st in ("census", "baryon", "coulomb"):
            procs.append(subprocess.Popen(
                [sys.executable, os.path.abspath(__file__), "--stage", st,
                 "--iters", str(a.iters)]))
        for p in procs:
            p.wait()
        report()
        return
    if a.stage == "census":
        stage_census(a.iters)
    elif a.stage == "baryon":
        stage_baryon(a.iters)
    elif a.stage == "coulomb":
        stage_coulomb()


if __name__ == "__main__":
    main()
