"""M5.32 R1.b: the PHYSICAL Coulomb gate (G1, pair half) on the static sector.

The certified two-hedgehog instrument (m5_21_4_a_pair.py: seeds, FIRE loop,
pinned shell, charge suite, core tracker) is re-run with the STATIC action
modified term by term through the M5.32 registry, the pair RELAXED under
the modified action, and the like-charge / anti-charge E(d) trends read.

EQUATIONS FIRST (what is minimized)
-----------------------------------
Field: M(x) real symmetric 3x3 on the n^3 lattice, h = L/n, pinned shell
(depth 1.6) held at the seed. Static 3x3 embed into the 4x4 registry:
    A_0 = 0,   A_i = 0 (+) d_i M   (i = 1..3, the certified fwd/bwd sym
    stencil, each branch's density evaluated and averaged),
    F_ij = A_i eta A_j - A_j eta A_i   (eta restricted to the block = 1,
    so F_ij = [d_i M, d_j M], the certified 3x3 commutator).
Energy relaxed (per term k in {I2, I6, I4} and coefficient c):
    E[M] = 4 h^3 sum_branch wt sum_cells ( I1(F) + c I_k(F) )
           + h^3 sum_cells V_T2(M),
    I1 = (1/2) F_{mu nu a b} F^{mu nu a b}   (the certified action; at
         c = 0 this is EXACTLY the instrument's E_u = 4 h^3 sum_{i<j}
         tr(C_ij C_ij^T), the control below proves it to 1e-10),
    I2 = F_{mu nu a b} F^{a b mu nu},  I6 = R^2 with R = F_{mu nu}^{nu mu},
    I4 = R_{nu a} R^{nu a}   (the redundancy check),
    V_T2 = w2 sum_lambda (lambda(M) - lambda_vac)^2, lambda_vac = (1, delta, 0)
         (the instrument's certified T2 eigenvalue potential; the registry's
         4x4 V4 is NOT used here: the pair instrument is the certified 3x3
         T2 stack and the control demands its record).
The c-grid is in units where I1 has coefficient 1 (the certified
CERTIFIED_COEFFS convention L = -4 I1 - V, E = 4 I1 + V, is carried as the
overall factor 4 h^3 on the curvature block).

Static identities (audited at R0, m5_32_r0_audit.json C3) with
G_mn = (1/4) eps_ijm eps_kln F_ij^kl, a = tr(G G^T), b = tr(G^2),
c_ = (tr G)^2:
    I1 = 2a,  I2 = 4b,  I6 = 4c_,  I3 = a + b,  I4 = a + c_,  I5 = b + c_
so on this sector I4 = I1/2 + I6/4 exactly (the redundancy check), and the
curvature block is bounded below iff
    I1 + c I2 :  2a + 4c b   >= 0 for all G  <=>  |c| <= 1/2  (b ranges
                 over [-a, a]: b = -a on antisymmetric G, b = a on symmetric G)
    I1 + c I6 :  2a + 4c c_  >= 0 for all G  <=>  c >= -1/6   (c_ <= 3a)
    I1 + c I4 :  (2 + c) a + c c_ >= 0       <=>  c >= -1/2
(these are predictions; the lattice run records DIVERGED where they bite).

Reads per (k, c): E_single (the single hedgehog relaxed under the same
action), E_same(d), E_anti(d), E_int = E_pair - 2 E_single, the fit
E_int(d) = A + B/d (B > 0: E_int falls with d, force -dE/dd > 0 =
repulsive), R^2 over the record window d = 12, 18, 24 and the extended
window with d = 30; the Mermin-Ho flux degrees (charge_suite: single far
cube, pair top / bot / far cubes); the single hedgehog half-radius
(INS.r_half, the potential-density median radius); the integrated
(a, b, c_) content and I3, I4, I5 of the relaxed single.

CONTROL: c = 0 reproduces m5_21_4_ladder_it120.json (E, E_u, E_single)
to 1e-10 at it = 120, n = 32, L = 48.

Gate (pre-registered G1, physical reading): like charges repel (B_same > 0)
with R^2 >= 0.95, anti charges attract (B_anti < 0), degree +-1 held
(single far flux within 0.1 of +-1, same-pair far flux within 0.2 of 2),
the hedgehog exists (E_single finite and positive, r_half finite). The
sign pattern RELATIVE to the c = 0 control is reported as its own column
(the record's own like-charge read is the M5.21.4 string-confinement
form, not a clean 1/d, so the absolute R^2 bar and the relative bar are
both printed).

Usage:
    python3 m5_32_r1_b_pair.py control            (the c = 0 reproduction)
    python3 m5_32_r1_b_pair.py grid [workers=12]  (the full grid, parallel)
    python3 m5_32_r1_b_pair.py plot               (re-plot from the json)

Out: ../data/m5_32_r1_pair.json, ../data/m5_32_r1_b/*.npz (local),
     ../plots/m5_32_r1_pair_gate.png
"""
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

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
OUT_JSON = os.path.join(DATA, "m5_32_r1_pair.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r1_b")
RECORD = os.path.join(DATA, "m5_21_4_ladder_it120.json")
PLOT = os.path.join(PLOTS, "m5_32_r1_pair_gate.png")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


PAIR = _load("m5_21_4_a_pair", "m5_21_4_a_pair.py")
INS = PAIR.INS                      # the instrument whose fire() we reuse
LAG = _load("m5_32_lagrangian", "m5_32_lagrangian.py")

TERMS = ("I2", "I6", "I4")
CGRID = (-1.0, -0.5, -0.25, 0.25, 0.5, 1.0)
DS = (12.0, 18.0, 24.0, 30.0)
IT, N, L = 120, 32, 48.0
BOUND_C = {"I2": (-0.5, 0.5), "I6": (-1.0 / 6.0, np.inf),
           "I4": (-0.5, np.inf)}   # analytic bounded-below windows

_K = {}


def K_of(name):
    if name not in _K:
        _K[name] = LAG.REGISTRY[name]._K()
    return _K[name]


def cfg_of():
    return INS.base_cfg(term="T2", n=N, L=L, bc="pinned", stencil="sym")


# ================= the modified static action =================
def jets4(M, h, br):
    """the static 3x3 embed: A (4, n, n, n, 4, 4), A_0 = 0, A_i = 0 (+) d_i M."""
    A = np.zeros((4,) + M.shape[:3] + (4, 4), dtype=M.dtype)
    for ax in range(3):
        A[1 + ax, ..., 1:, 1:] = INS.d1(M, ax, h, br)
    return A


def curv_integrals(M, cfg, names):
    """h^3 sum_branch wt sum_cells I_name(F) for each name (NO factor 4)."""
    h = cfg["h"]
    out = {k: 0.0 for k in names}
    for br, wt in INS.branches(cfg["stencil"]):
        F = LAG.F_of_A(jets4(M, h, br))
        for k in names:
            out[k] += wt * float(np.sum(LAG.density_from_K(F, K_of(k))))
    return {k: h ** 3 * v for k, v in out.items()}


class ModifiedAction:
    """E = 4 h^3 sum (I1 + c I_k) + V_T2; grad by the registry's exact
    K-matrix adjoint chained through the certified stencil adjoints."""

    def __init__(self, term, c):
        self.term, self.c = term, float(c)
        self.K = K_of("I1") + (self.c * K_of(term) if c != 0.0 else 0.0)

    def e_parts(self, M, cfg, st=None):
        st = st or cfg["stencil"]
        h = cfg["h"]
        e_u = 0.0
        for br, wt in INS.branches(st):
            F = LAG.F_of_A(jets4(M, h, br))
            e_u += wt * 4.0 * float(np.sum(LAG.density_from_K(F, self.K)))
        e_v = float(np.sum(INS.v_density(M, cfg)))
        return h ** 3 * e_u, 0.0, h ** 3 * e_v

    def grad(self, M, cfg):
        h = cfg["h"]
        G = np.zeros_like(M)
        for br, wt in INS.branches(cfg["stencil"]):
            A = jets4(M, h, br)
            F = LAG.F_of_A(A)
            dA = LAG.dA_from_W(LAG.dW_from_K(F, self.K), A)
            for ax in range(3):
                G += wt * 4.0 * INS.d1_adj(dA[1 + ax][..., 1:, 1:],
                                           ax, h, br)
        G = 0.5 * (G + G.swapaxes(-1, -2))
        G += INS.v_grad(M, cfg)
        return (h ** 3) * G

    def install(self):
        """route INS.fire's energy + gradient to this action (the FIRE
        loop, the pin mask, the seeds, the logging stay the certified
        code paths)."""
        INS.grad = self.grad
        INS.e_parts = self.e_parts


# ================= diagnostics =================
def fit_inv_d(ds, e):
    ds, e = np.asarray(ds, float), np.asarray(e, float)
    X = np.stack([np.ones_like(ds), 1.0 / ds], axis=1)
    coef, *_ = np.linalg.lstsq(X, e, rcond=None)
    res = e - X @ coef
    sst = float(np.sum((e - e.mean()) ** 2))
    r2 = 1.0 - float(np.sum(res ** 2)) / sst if sst > 0 else float("nan")
    return float(coef[0]), float(coef[1]), r2


def run_point(term, c, kind, d, save=True):
    """one relaxation: returns the row dict (JSON-safe)."""
    t0 = time.time()
    cfg = cfg_of()
    act = ModifiedAction(term, c)
    act.install()
    tag = f"{term}_c{c:+g}_{kind}_d{d:g}"
    try:
        M0, M, info = PAIR.heal(cfg, kind, d, IT)
    except FloatingPointError:
        return {"term": term, "c": c, "kind": kind, "d": d, "status":
                "DIVERGED", "stop": "exception", "wall_s": time.time() - t0}
    e_u, _, e_v = act.e_parts(M, cfg)
    E = e_u + e_v
    row = {"term": term, "c": c, "kind": kind, "d": d, "it": IT, "n": N,
           "L": L, "E": float(E), "E_curv": float(e_u), "E_v": float(e_v),
           "stop": info["stop"], "trace": info["trace"],
           "wall_s": time.time() - t0}
    finite = bool(np.isfinite(E)) and np.all(np.isfinite(M)) \
        and info["stop"] != "non-finite"
    row["status"] = "OK" if finite else "DIVERGED"
    if finite:
        # the trace's energy monotonicity (a bounded-below sanity read)
        Es = [t["E"] for t in info["trace"]]
        row["E_trace_min"] = float(min(Es))
        row["negative_curvature_energy"] = bool(e_u < 0.0)
        ints = curv_integrals(M, cfg, ("I1", "I2", "I3", "I4", "I5", "I6"))
        row["integrals"] = ints
        row["abc"] = {"a": ints["I1"] / 2.0, "b": ints["I2"] / 4.0,
                      "c_": ints["I6"] / 4.0}
        row["redundancy"] = {
            "I3_minus_(a+b)": ints["I3"] - (row["abc"]["a"] + row["abc"]["b"]),
            "I4_minus_(a+c_)": ints["I4"] - (row["abc"]["a"] + row["abc"]["c_"]),
            "I5_minus_(b+c_)": ints["I5"] - (row["abc"]["b"] + row["abc"]["c_"])}
        dq = d if d > 0 else 18.0
        q = PAIR.charge_suite(M, cfg, dq)
        row["charge"] = q
        if kind == "single":
            row["r_half"] = float(INS.r_half(M, cfg))
            lam = np.linalg.eigvalsh(M)
            row["min_gap_center"] = float(np.min(
                np.minimum(lam[..., 2] - lam[..., 1],
                           lam[..., 1] - lam[..., 0])))
        else:
            zs, gv = PAIR.core_zs(M, cfg, with_gaps=True)
            row["cores"], row["core_gaps"] = zs, gv
        if save:
            os.makedirs(OUT_NPZ, exist_ok=True)
            np.savez_compressed(os.path.join(OUT_NPZ, f"{tag}.npz"),
                                M=M.astype(np.float32))
    print(json.dumps({k: v for k, v in row.items() if k != "trace"}),
          flush=True)
    return row


def _worker(args):
    return run_point(*args)


# ================= assembly =================
def assemble(rows, control_ref=None):
    """per (term, c): E_single, E_int(d) both kinds, fits, verdict."""
    by = {}
    for r in rows:
        by.setdefault((r["term"], r["c"]), {})[(r["kind"], r["d"])] = r
    out = []
    for (term, c), pts in sorted(by.items()):
        s = pts.get(("single", 0.0))
        lo, hi = BOUND_C[term]
        rec = {"term": term, "c": c,
               "analytic_bounded": bool(lo - 1e-12 <= c <= hi + 1e-12)}
        statuses = [p["status"] for p in pts.values()]
        rec["n_runs"] = len(pts)
        rec["n_diverged"] = int(sum(st == "DIVERGED" for st in statuses))
        if s is None or s["status"] != "OK":
            rec["verdict"] = "DIVERGED (single)"
            rec["failure_mode"] = "divergence"
            out.append(rec)
            continue
        rec["E_single"] = s["E"]
        rec["E_single_curv"] = s["E_curv"]
        rec["r_half"] = s["r_half"]
        rec["single_far_flux"] = s["charge"]["far"]
        rec["abc_single"] = s["abc"]
        rec["integrals_single"] = s["integrals"]
        rec["redundancy_single"] = s["redundancy"]
        rec["exists"] = bool(np.isfinite(s["E"]) and s["E"] > 0
                             and np.isfinite(s["r_half"]))
        rec["negative_curvature_energy_single"] = s[
            "negative_curvature_energy"]
        deg_single_ok = abs(abs(s["charge"]["far"][1]) - 1.0) <= 0.1
        kinds = {}
        for kind in ("same", "anti"):
            E, E_int, dd, flux, st = [], [], [], {}, {}
            for d in DS:
                p = pts.get((kind, d))
                if p is None:
                    continue
                st[str(d)] = p["status"]
                if p["status"] != "OK":
                    continue
                dd.append(d)
                E.append(p["E"])
                E_int.append(p["E"] - 2.0 * s["E"])
                flux[str(d)] = {k: p["charge"][k] for k in
                                ("top", "bot", "far")}
            k = {"d": dd, "E": E, "E_int": E_int, "status": st,
                 "flux": flux}
            if len(dd) >= 3:
                sel = [i for i, d in enumerate(dd) if d <= 24.0]
                if len(sel) >= 3:
                    A, B, r2 = fit_inv_d([dd[i] for i in sel],
                                         [E_int[i] for i in sel])
                    k["fit_record_window"] = {"A": A, "B": B, "R2": r2}
                A, B, r2 = fit_inv_d(dd, E_int)
                k["fit_all"] = {"A": A, "B": B, "R2": r2}
                fw = k.get("fit_record_window", k["fit_all"])
                k["B_sign"] = "repulsive" if fw["B"] > 0 else "attractive"
                k["dEint_dd_sign_pairs"] = [
                    float(np.sign(E_int[i + 1] - E_int[i]))
                    for i in range(len(E_int) - 1)]
            kinds[kind] = k
        rec["kinds"] = kinds
        same, anti = kinds.get("same", {}), kinds.get("anti", {})
        rec["like_sign"] = same.get("B_sign", "n/a")
        rec["anti_sign"] = anti.get("B_sign", "n/a")
        rec["like_B"] = same.get("fit_record_window", {}).get("B")
        rec["like_R2"] = same.get("fit_record_window", {}).get("R2")
        rec["like_B_all"] = same.get("fit_all", {}).get("B")
        rec["like_R2_all"] = same.get("fit_all", {}).get("R2")
        rec["anti_B"] = anti.get("fit_record_window", {}).get("B")
        rec["anti_R2"] = anti.get("fit_record_window", {}).get("R2")
        same_far = [v["far"][1] for v in same.get("flux", {}).values()]
        rec["same_far_flux"] = same_far
        rec["degree_ok"] = bool(deg_single_ok and same_far
                                and all(abs(f - 2.0) <= 0.2
                                        for f in same_far))
        modes = []
        if rec["n_diverged"]:
            modes.append("divergence")
        if not rec["exists"]:
            modes.append("no hedgehog")
        if not rec["degree_ok"]:
            modes.append("degree loss")
        if rec["like_sign"] == "attractive":
            modes.append("like-charge sign reversal")
        if rec["anti_sign"] == "repulsive":
            modes.append("anti-charge sign reversal")
        if rec["like_R2"] is not None and rec["like_R2"] < 0.95:
            modes.append("1/d R2 < 0.95")
        gate_abs = (rec["like_sign"] == "repulsive"
                    and rec["anti_sign"] == "attractive"
                    and rec["degree_ok"] and rec["exists"]
                    and rec["n_diverged"] == 0
                    and (rec["like_R2"] or 0.0) >= 0.95)
        rec["gate_absolute"] = bool(gate_abs)
        if control_ref is not None:
            rec["same_sign_pattern_as_control"] = bool(
                rec["like_sign"] == control_ref["like_sign"]
                and rec["anti_sign"] == control_ref["anti_sign"]
                and rec["n_diverged"] == 0 and rec["exists"]
                and rec["degree_ok"])
        rec["failure_modes"] = modes
        rec["verdict"] = "PASS" if gate_abs else "FAIL"
        out.append(rec)
    return out


def write_json(payload):
    with open(OUT_JSON, "w") as f:
        json.dump(payload, f, indent=1, default=float)


# ================= phases =================
def control():
    """c = 0 through the registry path vs the record, plus the gradient
    identity vs the certified INS.grad on the seeds."""
    rec = json.load(open(RECORD))
    cfg = cfg_of()
    grad_cert, eparts_cert = INS.grad, INS.e_parts
    act = ModifiedAction("I2", 0.0)
    out = {"gradient_identity": {}, "energy_identity": {}, "record": {}}
    for kind, d in [("single", 0.0), ("anti", 18.0), ("same", 18.0)]:
        M0 = PAIR.seed_pair(cfg, kind, d)
        g0 = grad_cert(M0, cfg)
        g1 = act.grad(M0, cfg)
        out["gradient_identity"][f"{kind}_d{d:g}"] = {
            "max_abs_diff": float(np.max(np.abs(g0 - g1))),
            "max_abs_grad": float(np.max(np.abs(g0)))}
        e0, e1 = eparts_cert(M0, cfg), act.e_parts(M0, cfg)
        out["energy_identity"][f"{kind}_d{d:g}"] = {
            "E_u_cert": e0[0], "E_u_reg": e1[0],
            "rel": LAG._rel(e1[0], e0[0])}
    rows = []
    for r in rec["rows"]:
        row = run_point("I2", 0.0, r["kind"], r["d"], save=False)
        rows.append(row)
        out["record"][f"{r['kind']}_d{r['d']:g}"] = {
            "E_record": r["E"], "E_here": row["E"],
            "abs_diff": abs(r["E"] - row["E"]),
            "E_u_record": r["E_u"], "E_u_here": row["E_curv"],
            "abs_diff_u": abs(r["E_u"] - row["E_curv"])}
    srow = run_point("I2", 0.0, "single", 0.0, save=False)
    out["record"]["single"] = {"E_record": rec["E_single"],
                               "E_here": srow["E"],
                               "abs_diff": abs(rec["E_single"] - srow["E"])}
    worst = max(v["abs_diff"] for v in out["record"].values())
    out["control_pass_1e-10"] = bool(worst <= 1e-10)
    out["worst_abs_diff"] = worst
    INS.grad, INS.e_parts = grad_cert, eparts_cert
    print(json.dumps(out, indent=1, default=float), flush=True)
    return out, rows + [srow]


def grid(workers=12):
    t0 = time.time()
    ctrl, ctrl_rows = control()
    jobs = [("I2", 0.0, "single", 0.0)]      # a stored c = 0 single (npz)
    jobs += [("I2", 0.0, k, d) for k in ("same", "anti") for d in DS]
    for term in TERMS:
        for c in CGRID:
            jobs.append((term, c, "single", 0.0))
            for kind in ("same", "anti"):
                for d in DS:
                    jobs.append((term, c, kind, d))
    print(f"grid: {len(jobs)} relaxations on {workers} workers", flush=True)
    rows, done = [], 0
    ctx = mp.get_context("spawn")
    payload = {"control": ctrl, "grid": {"terms": TERMS, "c": CGRID,
                                         "d": DS, "it": IT, "n": N, "L": L},
               "analytic_bounds": BOUND_C, "rows": rows, "summary": []}
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as ex:
        futs = {ex.submit(_worker, j): j for j in jobs}
        for fut in as_completed(futs):
            j = futs[fut]
            try:
                rows.append(fut.result())
            except Exception as e:      # a crashed worker = DIVERGED row
                rows.append({"term": j[0], "c": j[1], "kind": j[2],
                             "d": j[3], "status": "DIVERGED",
                             "stop": f"exception: {e!r}"})
            done += 1
            if done % 10 == 0 or done == len(jobs):
                payload["elapsed_s"] = time.time() - t0
                payload["summary"] = summarize(rows)
                write_json(payload)
                print(f"[{done}/{len(jobs)} {time.time() - t0:.0f}s]",
                      flush=True)
    payload["elapsed_s"] = time.time() - t0
    payload["summary"] = summarize(rows)
    payload["windows"] = windows(payload["summary"])
    write_json(payload)
    plot(payload)
    print(f"grid done {time.time() - t0:.0f}s", flush=True)


def summarize(rows):
    # the c = 0 control rows are stored under term I2; they serve as the
    # control reference for every term
    ctrl = [r for r in rows if r["c"] == 0.0]
    ctrl_rec = assemble(ctrl)
    ref = ctrl_rec[0] if ctrl_rec else None
    out = assemble([r for r in rows if r["c"] != 0.0], control_ref=ref)
    if ref is not None:
        ref["term"] = "control"
        ref["same_sign_pattern_as_control"] = True
        out = ctrl_rec + out
    return out


def windows(summary):
    """the physical Coulomb window per term: the c values where the
    absolute gate holds, and where the relative (control-pattern) bar
    holds; the failure mode elsewhere."""
    w = {}
    for term in TERMS:
        recs = [r for r in summary if r["term"] == term]
        w[term] = {
            "absolute_pass_c": [r["c"] for r in recs if r.get("verdict")
                                == "PASS"],
            "relative_pass_c": [r["c"] for r in recs
                                if r.get("same_sign_pattern_as_control")],
            "fail": {str(r["c"]): r.get("failure_modes",
                                        [r.get("failure_mode", "")])
                     for r in recs if r.get("verdict") != "PASS"}}
    return w


def plot(payload=None):
    if payload is None:
        payload = json.load(open(OUT_JSON))
    summ = payload["summary"]
    ctrl = [r for r in summ if r["term"] == "control"]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), sharex=True)
    cmap = plt.get_cmap("coolwarm")
    for j, term in enumerate(TERMS):
        for i, kind in enumerate(("same", "anti")):
            ax = axes[i, j]
            if ctrl:
                k = ctrl[0]["kinds"].get(kind, {})
                ax.plot(k.get("d", []), k.get("E_int", []), "k-o", lw=2,
                        label="record (c = 0)")
            for r in summ:
                if r["term"] != term or "kinds" not in r:
                    continue
                k = r["kinds"].get(kind, {})
                if not k.get("d"):
                    continue
                col = cmap(0.5 + 0.5 * r["c"] / max(CGRID))
                ls = "-" if r.get("verdict") == "PASS" else "--"
                ax.plot(k["d"], k["E_int"], ls, marker="s", color=col,
                        label=f"c = {r['c']:+g} ({r.get('verdict')})")
            div = [r["c"] for r in summ if r["term"] == term
                   and r.get("n_diverged", 0)]
            ax.set_title(f"I1 + c {term}: {kind}-charge E_int(d)"
                         + (f"\nDIVERGED c = {div}" if div else ""),
                         fontsize=9)
            ax.set_xlabel("d")
            ax.set_ylabel("E_pair - 2 E_single")
            ax.grid(alpha=0.3)
            ax.legend(fontsize=7)
    fig.suptitle("M5.32 R1.b: the physical Coulomb gate on the pair "
                 f"instrument (n = {N}, L = {L:g}, it = {IT}; "
                 "solid = gate PASS, dashed = FAIL)")
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(PLOT, dpi=130)
    print(f"plot -> {PLOT}", flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "control"
    kv = {}
    for a in sys.argv[2:]:
        k, v = a.split("=")
        kv[k] = v
    if cmd == "control":
        control()
    elif cmd == "grid":
        grid(workers=int(kv.get("workers", 12)))
    elif cmd == "plot":
        plot()
    elif cmd == "assemble":          # re-summarize from the stored rows
        payload = json.load(open(OUT_JSON))
        payload["summary"] = summarize(payload["rows"])
        payload["windows"] = windows(payload["summary"])
        write_json(payload)
        plot(payload)
    else:
        raise SystemExit(f"unknown cmd {cmd}")
