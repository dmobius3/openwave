"""M5.32 R23-1 / R23-2: the halo radius against a second-order stiffness on the
degenerate pair, and the Coulomb-anchored delta.

EQUATIONS FIRST
---------------
N = M eta, the vacuum spectrum q = (-g, 1, delta, delta). V4 = w sum_p (tr N^p -
C_p)^2 is a sum of squares, so every first derivative in the invariants vanishes
on the vacuum and the split delta +/- eps of the degenerate pair costs eps^4
(R22-0 (c)). Inside V(tr N, ..., tr N^4) the spectrum stays critical only if
the first derivatives V_p satisfy sum_p p V_p x^(p-1) = 0 at the three distinct
eigenvalues, which fixes them up to one factor (R23-0 (a)):
    L = sum_p a_p (tr N^p - C_p),   sum_p p a_p x^(p-1) = (x + g)(x - 1)(x - delta),
    V = V4 - c L,   c >= 0   (the plan post writes V4 + c L with c < 0; its
                              "abs(c)" is this script's c),
    on the split:  -c L = c (g + delta)(1 - delta) eps^2 - c eps^4 / 2   (exact).
The static energy is E[M] = certified quartic + V4 - c L on the R22-1 lattice.

The instrument is R22-1s (m5_32_r22_1s_slaved.py): on a block-diagonal field
M_00 enters through the potential alone, so it is solved per cell by Newton on
    w sum_p (u^p + s_p - C_p)^2 - c sum_p a_p (u^p + s_p - C_p),  u = -M_00,
and L-BFGS descends on the six spatial entries; by the envelope theorem the
reduced gradient is the spatial block of the full gradient at the solved M_00.
With c = 0 every code path is the R22-1s one (the term is skipped, not
multiplied by zero), which is the known-answer control of the smoke.

Rows start from the analytic seed (m5_32_r22_1_cores.seed_core), not from an
R22 end field, so every box and spacing starts the same way; one row starts from
the box-filling R22-1s end field instead, so one c is approached from both sides.

Gate (relative, R23-2 needs it): max |G_spatial| < 1e-3 k(delta) / k(0.3) read
inside a chunk (never on a chunk's first iterate) AND the last-chunk drop under
1e-4 |E|. The slope fit of collect() uses a second, separate condition: the
half-energy radius drifts under 2 percent over the last 1000 iterations.

Reads per chunk: E and its parts (curvature, V4, -c L), the half-energy radius,
the shell profile of eps = (lam_mid - lam_small) / 2 on a 0.75 grid, the sorted
eigenvalues on the R22 shells, the surface-oriented degree of the top
eigenvector on the r 6 / 9 / 12 cubes (R21's reader).

Modes: smoke | run [workers] | extra [workers] | extra2 [workers] | collect | plot | jobs
Regenerate: smoke about 6 min; run about 8 h on 12 workers (16 jobs, the three
n 48 rows the long pole); collect seconds. Arrays: data/m5_32_r23_1/*.npz (local).
"""

import importlib.util
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r23_1_cscan.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r23_1")
R22S_NPZ = os.path.join(DATA, "m5_32_r22_1s")
R22_NPZ = os.path.join(DATA, "m5_32_r22_1")
R22S_JSON = os.path.join(DATA, "m5_32_r22_1s_slaved.json")
G = 8.0
CHUNK = 250
MAX_ITER = {32: 6000, 48: 4000}
GATE, GATE_DROP = 1e-3, 1e-4
DRIFT_WINDOW, DRIFT_MAX = 1000, 0.02
ALPHA = 1.0 / 137.035999
DELTA_PHYS = 1.0 - (ALPHA / (64.0 * np.pi)) ** 0.25  # 0.92238
K_RATIO = (0.7 / (1.0 - DELTA_PHYS)) ** 4  # 6615
EPS_R = np.arange(3.0, 34.6, 0.75)
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


SL = _load("m5_32_r22_1s_slaved", "m5_32_r22_1s_slaved.py")
CORES = _load("m5_32_r22_1_cores", "m5_32_r22_1_cores.py")
R21, R20, B3, R0 = SL.R21, SL.R20, SL.B3, SL.R0
W1 = B3.W1
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
IU3, OFF3 = SL.IU3, SL.OFF3


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


# ================= the linear member =================
def lin_coeffs(roots):
    """a_p, p = 1..4, with sum_p p a_p x^(p-1) = prod over the DISTINCT roots of (x - q)."""
    distinct = sorted(set(float(q) for q in roots))
    if len(distinct) != 3:
        raise ValueError("the linear member needs exactly three distinct vacuum eigenvalues")
    poly = np.poly(distinct)[::-1]  # ascending: coefficient of x^(p-1) at index p-1
    return [float(poly[p - 1]) / p for p in range(1, 5)]


def lin_energy_grad(M, cfg, roots, c, need_grad=True):
    """-c L, h^3-weighted total and the sym4 gradient (the conventions of R0.v4_energy_grad)."""
    h3 = cfg["h"] ** 3
    a = lin_coeffs(roots)
    Me = M @ ETA
    pows = [np.broadcast_to(np.eye(4), M.shape).copy()]
    for _ in range(1, 4):
        pows.append(pows[-1] @ Me)
    t = [np.einsum("...kk->...", P @ Me) for P in pows]
    cp = [sum(q**p for q in roots) for p in range(1, 5)]
    E = -c * h3 * np.sum(sum(a[k] * (t[k] - cp[k]) for k in range(4)))
    E = E if np.iscomplexobj(E) else float(E)
    if not need_grad:
        return E, None
    GV = np.zeros_like(M)
    for k in range(1, 5):
        GV += (-c * a[k - 1] * k) * (ETA @ pows[k - 1]).swapaxes(-1, -2)
    return E, h3 * B3.sym4(GV)


def lin_density(M, cfg, roots, c):
    a = lin_coeffs(roots)
    t = R20.LAG.v4_traces_np(M)
    cp = [sum(q**p for q in roots) for p in range(1, 5)]
    return -c * cfg["h"] ** 3 * sum(a[k] * (t[k] - cp[k]) for k in range(4))


def energy_grad(M, cfg, p, pot, c, need_grad=True):
    E, Gm, info = R20.energy_grad(M, cfg, p, pot, need_grad)
    if c == 0.0 or not np.isfinite(E):
        return E, Gm, info
    el, gl = lin_energy_grad(M, cfg, pot[1], c, need_grad)
    return E + el, (Gm + gl if need_grad else None), info


def solve_m00(S, roots, w, c, m_start, iters=14):
    """per cell argmin over M_00 of w sum_p r_p^2 - c sum_p a_p r_p, r_p = (-M_00)^p + tr S^p - C_p."""
    if c == 0.0:
        return SL.solve_m00(S, roots, m_start)
    a = lin_coeffs(roots)
    cp = [sum(q**k for q in roots) for k in range(1, 5)]
    s, P = [], S
    for k in range(1, 5):
        if k > 1:
            P = P @ S
        s.append(np.einsum("...kk->...", P))
    x = -m_start.copy()
    for _ in range(iters):
        r = [x**k + s[k - 1] - cp[k - 1] for k in range(1, 5)]
        d = [k * x ** (k - 1) for k in range(1, 5)]
        dd = [k * (k - 1) * x ** (k - 2) if k > 1 else 0.0 * x for k in range(1, 5)]
        g = sum((2.0 * w * r[k] - c * a[k]) * d[k] for k in range(4))
        H = sum(2.0 * w * d[k] ** 2 + (2.0 * w * r[k] - c * a[k]) * dd[k] for k in range(4))
        x = x - g / H
    return -x


class Reduced:
    def __init__(self, M, cfg, p, pot, c, pinned):
        self.cfg, self.p, self.pot, self.c = cfg, p, pot, c
        self.mask = R21.free_mask(cfg, pinned)
        self.M = M.copy()
        self.last = None

    def build(self, x):
        M = self.M.copy()
        S = M[..., 1:, 1:].copy()
        blk = np.zeros((int(self.mask.sum()), 3, 3), dtype=x.dtype)
        blk[:, IU3[0], IU3[1]] = x.reshape(-1, 6)
        blk = blk + blk.swapaxes(-1, -2) - np.einsum("...ii->...i", blk)[..., None] * np.eye(3)
        S[self.mask] = blk
        M[..., 1:, 1:] = S
        m00 = M[..., 0, 0].copy()
        m00[self.mask] = solve_m00(S[self.mask], self.pot[1], self.pot[2], self.c, m00[self.mask])
        M[..., 0, 0] = m00
        return M

    def pack(self, M):
        return M[..., 1:, 1:][self.mask][:, IU3[0], IU3[1]].ravel().copy()

    def fun(self, x):
        M = self.build(x)
        E, Gm, info = energy_grad(M, self.cfg, self.p, self.pot, self.c)
        if Gm is None or not np.isfinite(E) or not info["ok"]:
            return 1e30, np.zeros_like(x)
        Gf = Gm[self.mask]
        gs = Gf[:, 1:, 1:][:, IU3[0], IU3[1]] * OFF3
        self.last = (
            float(E),
            float(np.max(np.abs(Gf[:, 1:, 1:]))),
            float(np.max(np.abs(Gf[:, 0, 0]))),
        )
        self.M[..., 0, 0] = M[..., 0, 0]
        return float(E), gs.ravel()


# ================= reads =================
def chunk_reads(M, cfg, p, pot, c):
    parts = R20.energy_parts(M, cfg, p, pot)
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    e = R20.density(M, cfg, pot)
    E_lin = 0.0
    if c != 0.0:
        el = lin_density(M, cfg, pot[1], c)
        e = e + el
        E_lin = float(el.sum())
    e = e.ravel()
    o = np.argsort(r.ravel())
    cum = np.cumsum(e[o])
    lam = np.linalg.eigvalsh(M[..., 1:, 1:])
    eps = 0.5 * (lam[..., 1] - lam[..., 0])
    prof = []
    for R in EPS_R:
        if R > 0.5 * cfg["L"] - 2.0 * h:
            break
        sh = np.abs(r - R) < 0.5 * h
        if sh.any():
            prof.append([float(R), float(eps[sh].mean()), float(lam[..., 2][sh].mean())])
    shells = []
    for R in SL.SHELL_R:
        sh = np.abs(r - R) < 0.75 * h
        if sh.any():
            shells.append([R] + [float(lam[..., a][sh].mean()) for a in range(3)])
    deg = R21.degree_surface_reads(M, cfg)["rank2"]
    return {
        "E": parts["E_total"] + E_lin,
        "E_curv": parts["E_curv"],
        "V4": parts["V"],
        "E_lin": E_lin,
        "r_half": float(r.ravel()[o][np.searchsorted(cum, 0.5 * cum[-1])]),
        "e_min_cell": float(e.min()),
        "eps_profile_r_eps_top": prof,
        "degree_top": {k: [v["degree"], v["surface_conflicts"]] for k, v in deg.items()},
        "shells_r_small_mid_top": shells,
        "M00_range": [float(M[..., 0, 0].min()), float(M[..., 0, 0].max())],
    }


def r_eps_of(prof, p_tail=0.6180339887498949, r_ref=4.5):
    """the radius where eps r^p_tail falls to half of its value at r_ref (linear interpolation)."""
    R = np.array([q[0] for q in prof])
    y = np.array([q[1] for q in prof]) * R**p_tail
    y0 = float(np.interp(r_ref, R, y))
    if y0 <= 0:
        return None
    k = np.where((R > r_ref) & (y < 0.5 * y0))[0]
    if len(k) == 0:
        return None  # never falls to half inside the read window
    i = k[0]
    return float(R[i - 1] + (0.5 * y0 - y[i - 1]) * (R[i] - R[i - 1]) / (y[i] - y[i - 1]))


# ================= jobs =================
def job_tag(j):
    src = "" if j.get("src") is None else "_from_r22s"
    return (
        f"{j['seed']}_pin_d{j['delta']:.4g}_w{j['w1s']:.4g}_c{j['c']:g}_n{j['n']}_L{j['L']:g}{src}"
    )


def jobs_all():
    J = []
    base = {"delta": 0.3, "w1s": 25.0, "seed": "rad", "src": None}
    for c in (3e-5, 1e-4):  # the box check first: the long pole
        J.append(dict(base, part="R23-1 box", c=c, n=48, L=72.0))
    J.append(dict(base, part="R23-1 spacing", c=3e-4, n=48, L=48.0))
    for c in (0.0, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3):
        J.append(dict(base, part="R23-1 scan", c=c, n=32, L=48.0))
    J.append(
        dict(
            base,
            part="R23-1 two-sided",
            c=3e-4,
            n=32,
            L=48.0,
            src="rad_pin_d0.3_w25_n32_L48_x9000_s",
        )
    )
    for c in (1e-4, 1e-3):
        J.append(dict(base, part="R23-1 seed", seed="bia", c=c, n=32, L=48.0))
    # R23-2: c chosen so that m2 / 2K equals the delta 0.3, c 3e-4 row (the same halo equation)
    d = float(DELTA_PHYS)
    c_phys = (
        3e-4 * ((G + 0.3) * 0.7 / (8.0 * 0.49)) / ((G + d) * (1.0 - d) / (8.0 * (1.0 - d) ** 2))
    )
    # row B runs at c = 0 only: R23-0 (c) shows the vacuum unstable there for c > 7.5e-7 (deviation 1)
    for w1s, nm, cl in ((25.0, "A", (0.0, c_phys)), (25.0 / K_RATIO, "B", (0.0,))):
        for c in cl:
            J.append(
                {
                    "part": f"R23-2 row {nm}",
                    "delta": d,
                    "w1s": w1s,
                    "seed": "rad",
                    "src": None,
                    "c": float(f"{c:.4g}"),
                    "n": 32,
                    "L": 48.0,
                }
            )
    return J


def jobs_extra():
    """added at EXECUTE (deviation 7): the L 48 box limits the halo for c <= 3e-4, so the scan is repeated on
    n 48 L 72 at the three larger c values; and the c 3e-3 row, which met the gate at 1000 iterations with its
    half-energy radius still moving 4 percent per 500 iterations, is continued past the gate (the drift rule).
    """
    base = {"delta": 0.3, "w1s": 25.0, "seed": "rad", "src": None}
    J = [dict(base, part="R23-1 box", c=c, n=48, L=72.0) for c in (3e-4, 1e-3, 3e-3)]
    J.append(dict(base, part="R23-1 scan", c=3e-3, n=32, L=48.0, min_iters=3000, force=True))
    return J


def jobs_extra2():
    """added at EXECUTE on the author's 2026-09-20 01:26 UTC reply (deviation 9): the ladder extended upward to
    seven points over three decades (c 1e-2 and 3e-2, both under c_crit 0.086), and the three highest c values
    repeated at h 1.0 (n 48 L 48), since the halo approaches the core there and 3 h is the lattice floor.
    """
    base = {"delta": 0.3, "w1s": 25.0, "seed": "rad", "src": None}
    J = [dict(base, part="R23-1 scan, upper ladder", c=c, n=32, L=48.0) for c in (1e-2, 3e-2)]
    J += [
        dict(base, part="R23-1 spacing, upper ladder", c=c, n=48, L=48.0)
        for c in (3e-3, 1e-2, 3e-2)
    ]
    return J


def gate_of(delta):
    return GATE * ((1.0 - delta) / 0.7) ** 4


def run_job(j):
    from scipy.optimize import minimize

    t0 = time.time()
    tag = job_tag(j)
    ug = SL.ups_guard()
    cfg = R21.cfg_of(j["n"], j["L"], G, j["delta"])
    p = R21.params_of(G, j["delta"])
    pot = ("v4", R0.roots_of(cfg, degenerate=True), W1 * j["w1s"])
    c, gate = j["c"], gate_of(j["delta"])
    os.makedirs(OUT_NPZ, exist_ok=True)
    stage = os.path.join(OUT_NPZ, tag + "_stage.npz")
    row = dict(j, tag=tag, gate=gate)
    try:
        if os.path.exists(stage):
            Zs = np.load(stage, allow_pickle=True)
            M, done, chunks = Zs["M"], int(Zs["done"]), json.loads(str(Zs["chunks"]))
        else:
            if j.get("src"):
                M = np.load(os.path.join(R22S_NPZ, j["src"] + ".npz"))["M"]
            else:
                M = CORES.seed_core(cfg, j["seed"], j["delta"])
            if np.abs(M[..., 0, 1:]).max() > 0:
                raise RuntimeError("the start field is not block-diagonal")
            done, chunks = 0, [
                dict(chunk_reads(M, cfg, p, pot, c), iters=0, note="the start field")
            ]
        verdict = "FALLING"
        while done < MAX_ITER[j["n"]]:
            red = Reduced(M, cfg, p, pot, c, True)
            st = {"it": 0, "gate_hit": None}

            def cb(xk, red=red, st=st):
                st["it"] += 1
                if st["it"] >= 5 and red.last is not None and red.last[1] < gate:
                    st["gate_hit"] = st["it"]

            res = minimize(
                red.fun,
                red.pack(M),
                jac=True,
                method="L-BFGS-B",
                callback=cb,
                options={
                    "maxcor": 20,
                    "maxiter": CHUNK,
                    "maxfun": 3 * CHUNK,
                    "gtol": 1e-14,
                    "ftol": 1e-16,
                },
            )
            M = red.build(np.asarray(res.x))
            done += max(1, st["it"])
            _, Gm, _ = energy_grad(M, cfg, p, pot, c)
            Gf = Gm[red.mask]
            rec = dict(
                chunk_reads(M, cfg, p, pot, c),
                iters=done,
                fmax_spatial=float(np.max(np.abs(Gf[:, 1:, 1:]))),
                fmax_M00=float(np.max(np.abs(Gf[:, 0, 0]))),
                gate_hit_inside_chunk=st["gate_hit"],
                chunk_iters=st["it"],
                scipy=str(res.message),
            )
            rec["drop"] = chunks[-1]["E"] - rec["E"]
            chunks.append(rec)
            np.savez_compressed(stage, M=M, done=done, chunks=json.dumps(chunks))
            log(
                f"{tag} its {done} E {rec['E']:.6g} drop {rec['drop']:.2e} fmax_sp {rec['fmax_spatial']:.2e} "
                f"(gate {gate:.1e}) r_half {rec['r_half']:.2f} r_eps {r_eps_of(rec['eps_profile_r_eps_top'])} "
                f"deg {rec['degree_top']}"
            )
            if (
                rec["fmax_spatial"] < gate
                and 0 <= rec["drop"] < GATE_DROP * abs(rec["E"])
                and st["it"] >= 5
            ):
                verdict = "AT_GATE"
                if done >= j.get("min_iters", 0):
                    break
            if st["it"] < 5 and rec["drop"] <= 0:
                verdict = "LINE_SEARCH_STALL"
                break
            if ug is not None and ug.wrap_up():
                verdict = "UPS wrap-up (resumable)"
                break
        row.update(chunks=chunks, iters=done, label=verdict, E=chunks[-1]["E"], status="OK")
        np.savez_compressed(os.path.join(OUT_NPZ, tag + ".npz"), M=M)
    except Exception as e:  # noqa: BLE001
        import traceback

        row.update(status="FAILED", stop=repr(e), traceback=traceback.format_exc())
    row["wall_s"] = round(time.time() - t0, 1)
    log(f"DONE {tag} {row.get('label')} E {row.get('E')} wall {row['wall_s']}")
    return row


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R23-1 / R23-2", "rows": {}}


EXTRA_JSON = OUT_JSON.replace(".json", "_extra.json")


def load_all_rows():
    """the main pool's rows, overridden by the extra pool's rows (its own file, so two pools never race)."""
    rows = dict(load_json()["rows"])
    for fn in (EXTRA_JSON, EXTRA_JSON.replace("_extra.json", "_extra2.json")):
        if os.path.exists(fn):
            with open(fn) as f:
                rows.update(json.load(f)["rows"])
    return rows


def run_pool(jobs, workers):
    rows = load_json()["rows"]
    pending = [j for j in jobs if j.get("force") or rows.get(job_tag(j), {}).get("status") != "OK"]
    log(f"pool: {len(pending)} jobs, {workers} workers")
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        futs = [ex.submit(run_job, j) for j in pending]
        for fut in as_completed(futs):
            row = fut.result()
            Jn = load_json()
            Jn["rows"][row["tag"]] = row
            tmp = OUT_JSON + ".tmp"
            with open(tmp, "w") as f:
                json.dump(Jn, f, indent=1)
            os.replace(tmp, OUT_JSON)
    log("pool done")


# ================= smoke =================
def smoke():
    from scipy.optimize import minimize

    out = {}
    cfg = R21.cfg_of(32, 48.0, G, 0.3)
    p = R21.params_of(G, 0.3)
    roots = R0.roots_of(cfg, degenerate=True)
    pot = ("v4", roots, W1 * 25.0)
    c = 3e-4
    out["a_p"] = lin_coeffs(roots)
    # (1) -c L on a uniform split field: the closed form c (g + delta)(1 - delta) eps^2 - c eps^4 / 2
    e_ = 0.05
    Mu = np.broadcast_to(np.diag([G, 1.0, 0.3 + e_, 0.3 - e_]), (4, 4, 4, 4, 4)).copy()
    cfg4 = dict(cfg, n=4)
    got = lin_energy_grad(Mu, cfg4, roots, c, False)[0] / (64 * cfg["h"] ** 3)
    want = c * (G + 0.3) * 0.7 * e_**2 - c * e_**4 / 2
    out["split_closed_form_rel"] = abs(got - want) / want
    # (2) the gradient of -c L against complex step, full 4x4 random symmetric field
    rng = np.random.default_rng(3)
    Mr = Mu + 0.05 * B3.sym4(rng.normal(size=Mu.shape))
    _, gl = lin_energy_grad(Mr, cfg4, roots, c)
    cs = []
    for _ in range(4):
        V = B3.sym4(rng.normal(size=Mu.shape))
        num = np.imag(lin_energy_grad(Mr + 1e-30j * V, cfg4, roots, c, False)[0]) / 1e-30
        cs.append(abs(num - float(np.sum(gl * V))) / abs(num))
    out["lin_grad_complex_step_rel"] = cs
    # (3) the slaved solve and the reduced gradient with c != 0, on an R22 end field
    M = np.load(os.path.join(R22_NPZ, "bia_pin_d0.3_w25_n32_L48.npz"))["M"]
    red = Reduced(M, cfg, p, pot, c, True)
    x0 = red.pack(M)
    E0, G0, _ = energy_grad(M, cfg, p, pot, c)
    E1, g1 = red.fun(x0)
    out["E_unslaved"], out["E_slaved"] = float(E0), E1
    out["max_G00_unslaved"] = float(np.max(np.abs(G0[red.mask][:, 0, 0])))
    out["max_G00_slaved"] = red.last[2]
    fd = []
    for _ in range(3):
        v = rng.normal(size=x0.shape)
        v /= np.linalg.norm(v)
        t = 1e-4
        num = (red.fun(x0 + t * v)[0] - red.fun(x0 - t * v)[0]) / (2 * t)
        fd.append(abs(num - float(g1 @ v)) / max(abs(num), 1e-300))
    out["reduced_grad_fd_rel"] = fd
    # (4) the known-answer control: c = 0, one chunk from the R22-1s source, against the R22-1s record
    src = "rad_pin_d0.3_w25_n32_L48_x9000"
    M = np.load(os.path.join(R22_NPZ, src + ".npz"))["M"]
    red0 = Reduced(M, cfg, p, pot, 0.0, True)
    res = minimize(
        red0.fun,
        red0.pack(M),
        jac=True,
        method="L-BFGS-B",
        options={
            "maxcor": 20,
            "maxiter": CHUNK,
            "maxfun": 3 * CHUNK,
            "gtol": 1e-14,
            "ftol": 1e-16,
        },
    )
    E_c0 = chunk_reads(red0.build(np.asarray(res.x)), cfg, p, pot, 0.0)["E"]
    with open(R22S_JSON) as f:
        ref = json.load(f)["rows"][src + "_s"]["chunks"][1]["E"]
    out["c0_first_chunk_E"], out["r22_1s_first_chunk_E"] = E_c0, ref
    out["c0_control_abs"] = abs(E_c0 - ref)
    # (5) the seed is block-diagonal and the pinned shell is untouched
    Ms = CORES.seed_core(cfg, "rad", 0.3)
    reds = Reduced(Ms, cfg, p, pot, c, True)
    Mb = reds.build(reds.pack(Ms))
    out["seed_block_diagonal"] = bool(np.abs(Ms[..., 0, 1:]).max() == 0)
    out["pinned_shell_unchanged"] = bool(np.array_equal(Mb[~reds.mask], Ms[~reds.mask]))
    out["PASS"] = bool(
        out["split_closed_form_rel"] < 1e-10
        and max(cs) < 1e-10
        and E1 <= E0 + 1e-12
        and out["max_G00_slaved"] < 1e-6 * max(out["max_G00_unslaved"], 1e-300) + 1e-9
        and max(fd) < 1e-5
        and out["c0_control_abs"] < 1e-10
        and out["seed_block_diagonal"]
        and out["pinned_shell_unchanged"]
    )
    print(json.dumps(out, indent=1))
    with open(OUT_JSON.replace(".json", "_smoke.json"), "w") as f:
        json.dump(out, f, indent=1)


# ================= post-audit reads (added after m5_32_r23_1_audit.py; NOT pre-registered) =================
def r_half_interp(tag, row):
    """the half-energy radius with the cumulative energy interpolated between distinct cell radii (the stored
    r_half is the radius of the crossing cell; audit check C7: one step of it is 4 percent at r 7.5).
    """
    f = os.path.join(OUT_NPZ, tag + ".npz")
    if not os.path.exists(f):
        return None
    M = np.load(f)["M"]
    cfg = R21.cfg_of(row["n"], row["L"], G, row["delta"])
    pot = ("v4", R0.roots_of(cfg, degenerate=True), W1 * row["w1s"])
    e = R20.density(M, cfg, pot)
    if row["c"] != 0.0:
        e = e + lin_density(M, cfg, pot[1], row["c"])
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z).ravel()
    ru, inv = np.unique(np.round(r, 9), return_inverse=True)
    cum = np.cumsum(np.bincount(inv, weights=e.ravel()))
    return float(np.interp(0.5 * cum[-1], cum, ru))


def cutoff_fits(prof, L):
    """one-scale fits of the shell profile eps(r) on 5 <= r <= L/2 - 6, in log eps, amplitude free:
    K: sqrt(r) K_nu(r^2 / 2 R^2), nu = sqrt(5)/4 (the decaying solution of the 1/r^2-stiffness equation);
    Y: exp(-r / R) / r (a constant-stiffness Klein-Gordon halo)."""
    from scipy.optimize import minimize_scalar
    from scipy.special import kve

    P = np.array([q[:2] for q in prof])
    k = (P[:, 0] >= 5.0) & (P[:, 0] <= 0.5 * L - 6.0) & (P[:, 1] > 0)
    k &= P[:, 1] > 1e-3 * np.interp(
        5.0, P[:, 0], P[:, 1]
    )  # above the floor of the eigenvalue solver's noise
    r, y = P[k, 0], np.log(P[k, 1])
    nu = np.sqrt(5.0) / 4.0
    shapes = {
        "K_nu": lambda R: 0.5 * np.log(r)
        + np.log(kve(nu, r * r / (2 * R * R)))
        - r * r / (2 * R * R),
        "yukawa": lambda R: -r / R - np.log(r),
    }
    out = {}
    for nm, fn in shapes.items():

        def res(lr, fn=fn):
            m = fn(np.exp(lr))
            return float(np.sqrt(np.mean((y - m - np.mean(y - m)) ** 2)))

        o = minimize_scalar(res, bounds=(np.log(1.0), np.log(400.0)), method="bounded")
        out[nm] = {"R": float(np.exp(o.x)), "rms_log_residual": float(o.fun)}
    m = -0.6180339887498949 * np.log(r)
    out["free_power_law"] = {
        "rms_log_residual": float(np.sqrt(np.mean((y - m - np.mean(y - m)) ** 2)))
    }
    return out


# ================= collect =================
def row_summary(row):
    ch = row["chunks"]
    last = ch[-1]
    win = [q for q in ch if q["iters"] >= last["iters"] - DRIFT_WINDOW and q["iters"] > 0]
    rh = np.array([q["r_half"] for q in win])
    re_ = [r_eps_of(q["eps_profile_r_eps_top"]) for q in win]
    out = {
        "part": row["part"],
        "label": row["label"],
        "iters": row["iters"],
        "E": last["E"],
        "E_curv": last["E_curv"],
        "V4": last["V4"],
        "E_lin": last["E_lin"],
        "fmax_spatial": last.get("fmax_spatial"),
        "gate": row["gate"],
        "last_drop_over_E": last.get("drop", 0) / abs(last["E"]),
        "r_half": last["r_half"],
        "r_half_drift": float((rh.max() - rh.min()) / rh.mean()),
        "r_eps": re_[-1],
        "r_eps_drift": (
            None if any(v is None for v in re_) else float((max(re_) - min(re_)) / np.mean(re_))
        ),
        "degree_top": last["degree_top"],
        "interior_spectrum_r1.5_3_4.5": last["shells_r_small_mid_top"][:3],
        "e_min_cell": last["e_min_cell"],
        "wall_s": row["wall_s"],
    }
    if (
        row["label"] == "LINE_SEARCH_STALL"
        and out["fmax_spatial"] is not None
        and out["fmax_spatial"] < row["gate"]
    ):
        out["label"] = "BELOW_GATE (the line search stalled under the gate)"
    out["CONVERGED_FOR_FIT"] = bool(len(win) >= 4 and out["r_half_drift"] < DRIFT_MAX)
    return out


def slope(cs, rs):
    return float(np.polyfit(np.log(cs), np.log(rs), 1)[0]) if len(cs) >= 3 else None


def collect():
    rows = {t: r for t, r in load_all_rows().items() if r.get("status") == "OK"}
    S = {t: row_summary(r) for t, r in rows.items()}
    res = {"rows": S}
    tag = lambda c, n=32, L=48.0, seed="rad": job_tag(
        {"seed": seed, "delta": 0.3, "w1s": 25.0, "c": c, "n": n, "L": L, "src": None}
    )  # noqa: E731
    box = {}
    for c in (3e-5, 1e-4):
        a_, b_ = S.get(tag(c)), S.get(tag(c, 48, 72.0))
        if a_ and b_:
            box[str(c)] = {
                k: [
                    a_[k],
                    b_[k],
                    None if (a_[k] is None or b_[k] is None) else abs(a_[k] - b_[k]) / b_[k],
                ]
                for k in ("r_half", "r_eps")
            }
            box[str(c)]["both_converged"] = a_["CONVERGED_FOR_FIT"] and b_["CONVERGED_FOR_FIT"]
    res["box_check_L48_vs_L72"] = box
    sp_a, sp_b = S.get(tag(3e-4)), S.get(tag(3e-4, 48, 48.0))
    if sp_a and sp_b:
        res["spacing_check_c3e-4"] = {k: [sp_a[k], sp_b[k]] for k in ("r_half", "r_eps", "E")}
    two = S.get(tag(3e-4) + "_from_r22s")
    if sp_a and two:
        res["two_sided_c3e-4"] = {k: [sp_a[k], two[k]] for k in ("r_half", "r_eps", "E")}
    for c in (1e-4, 1e-3):
        a_, b_ = S.get(tag(c)), S.get(tag(c, seed="bia"))
        if a_ and b_:
            res[f"seed_check_c{c:g}"] = {k: [a_[k], b_[k]] for k in ("r_half", "r_eps", "E")}
    # the fit: converged rows; a c passes the box check if the largest box-checked c at or below it passed
    passed = [float(c) for c, v in box.items() if v["both_converged"] and v["r_half"][2] < 0.10]
    fit = {}
    for read in ("r_half", "r_eps"):
        pts = [
            (c, S[tag(c)][read])
            for c in (3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2)
            if tag(c) in S
            and S[tag(c)]["CONVERGED_FOR_FIT"]
            and S[tag(c)][read] is not None
            and any(cp <= c for cp in passed)
        ]
        fit[read] = {
            "points": pts,
            "slope": slope([q[0] for q in pts], [q[1] for q in pts]),
            "slope_lowest_three": slope([q[0] for q in pts[:3]], [q[1] for q in pts[:3]]),
        }
    res["fit"] = fit
    res["all_rows_r_half_vs_c"] = [
        (c, S[tag(c)]["r_half"], S[tag(c)]["r_eps"], S[tag(c)]["CONVERGED_FOR_FIT"])
        for c in (0.0, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2)
        if tag(c) in S
    ]
    sl = [fit[k]["slope"] for k in ("r_half", "r_eps")]
    seed_core_r = 4.4
    if all(
        abs(S[tag(c)]["r_half"] - seed_core_r) < 1.5 * 1.5
        for c in (3e-5, 1e-4, 3e-4, 1e-3, 3e-3)
        if tag(c) in S
    ):
        lab = "HALO_LATTICE_PINNED"
    elif any(v is None for v in sl) or len(fit["r_half"]["points"]) < 3:
        lab = "HALO_BOX_LIMITED"
    elif all(-0.6 <= v <= -0.4 for v in sl):
        lab = "HALO_HALF"
    elif all(-0.3 <= v <= -0.2 for v in sl):
        lab = "HALO_QUARTER"
    else:
        lab = "HALO_OTHER"
    res["label_R23_1"] = lab
    # post-audit reads, per box: not pre-registered, reported beside the label
    post = {}
    for n, L in ((32, 48.0), (48, 72.0), (48, 48.0)):
        tab = []
        for c in (0.0, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2):
            t = tag(c, n, L)
            if t not in S:
                continue
            fits = cutoff_fits(rows[t]["chunks"][-1]["eps_profile_r_eps_top"], L)
            tab.append(
                {
                    "c": c,
                    "label": S[t]["label"],
                    "iters": S[t]["iters"],
                    "r_half": S[t]["r_half"],
                    "r_half_interp": r_half_interp(t, rows[t]),
                    "r_eps": S[t]["r_eps"],
                    **fits,
                }
            )
        if not tab:
            continue
        post[f"n{n}_L{L:g}"] = {"rows": tab}
        cut = lambda R_: (
            R_ if R_ < 0.5 * L else None
        )  # noqa: E731  a fitted scale beyond the half-box is "not cut"
        for key, get in (
            ("r_eps", lambda q: q["r_eps"]),
            ("K_nu_R", lambda q: cut(q["K_nu"]["R"])),
            ("yukawa_R", lambda q: cut(q["yukawa"]["R"])),
            ("r_half_interp", lambda q: q["r_half_interp"]),
        ):
            for cmin in (1e-4, 3e-4):
                pts = [(q["c"], get(q)) for q in tab if q["c"] >= cmin and get(q) is not None]
                if len(pts) >= 3:
                    post[f"n{n}_L{L:g}"][f"slope_{key}_c>={cmin:g}"] = slope(
                        [a_ for a_, _ in pts], [b_ for _, b_ in pts]
                    )
    res["post_audit_reads_NOT_preregistered"] = post
    with open(OUT_JSON.replace(".json", "_collect.json"), "w") as f:
        json.dump(res, f, indent=1)
    for t, v in S.items():
        print(
            f"{t:58s} {v['label']:10s} its {v['iters']:5d} E {v['E']:.5g} r_half {v['r_half']:.2f} (drift {v['r_half_drift']:.3f}) "
            f"r_eps {v['r_eps']} fmax {v['fmax_spatial']:.1e} fit-ok {v['CONVERGED_FOR_FIT']}"
        )
    print(json.dumps({k: v for k, v in res.items() if k != "rows"}, indent=1, default=str))


# ================= plot =================
def plot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = {t: r for t, r in load_all_rows().items() if r.get("status") == "OK"}
    with open(os.path.join(DATA, "m5_32_r23_0_form.json")) as f:
        pred = {q["c"]: q["reader_radius_f1"] for q in json.load(f)["f"]["rows"]}
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
    cs = (3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2)
    for n, L, mk, nm in (
        (32, 48.0, "o", "n 32 L 48"),
        (48, 72.0, "s", "n 48 L 72"),
        (48, 48.0, "^", "n 48 L 48"),
    ):
        pts = []
        for c in cs:
            t = job_tag(
                {"seed": "rad", "delta": 0.3, "w1s": 25.0, "c": c, "n": n, "L": L, "src": None}
            )
            if t in rows:
                S = row_summary(rows[t])
                pts.append((c, S["r_half"], S["r_eps"], S["CONVERGED_FOR_FIT"]))
        if pts:
            ax[0].loglog(
                [q[0] for q in pts],
                [q[1] for q in pts],
                mk + "-",
                label=f"half-energy radius, {nm}",
            )
            pe = [q for q in pts if q[2] is not None]
            ax[0].loglog(
                [q[0] for q in pe],
                [q[2] for q in pe],
                mk + "--",
                mfc="none",
                label=f"split-profile radius, {nm}",
            )
    cc = np.array(cs)
    cp = [c for c in cs if c in pred]
    ax[0].loglog(
        cp,
        [pred[c] for c in cp],
        "k:",
        label="split-profile reader on the c^(-1/4) solution (R23-0 f)",
    )
    ax[0].loglog(cc, 12.0 * (cc / 3e-4) ** -0.5, color="0.6", lw=1, label="slope -1/2")
    ax[0].loglog(cc, 12.0 * (cc / 3e-4) ** -0.25, color="0.6", lw=1, ls="-.", label="slope -1/4")
    ax[0].set_xlabel("c"), ax[0].set_ylabel("radius"), ax[0].set_title(
        "halo radius against the stiffness c"
    )
    ax[0].legend(fontsize=7)
    for c in (0.0,) + cs:
        t = job_tag(
            {"seed": "rad", "delta": 0.3, "w1s": 25.0, "c": c, "n": 32, "L": 48.0, "src": None}
        )
        if t in rows:
            pr = np.array(rows[t]["chunks"][-1]["eps_profile_r_eps_top"])
            ax[1].loglog(pr[:, 0], np.maximum(pr[:, 1], 1e-8), label=f"c = {c:g}")
    rr = np.array([3.0, 22.0])
    ax[1].loglog(rr, 0.12 * (rr / 3.0) ** -0.618, "k:", label="r^(-0.618)")
    ax[1].set_xlabel("r"), ax[1].set_ylabel("shell mean of eps = (lam_mid - lam_small) / 2")
    ax[1].set_title("the split of the pair, n 32 L 48"), ax[1].legend(fontsize=7)
    for c in (0.0,) + cs:
        t = job_tag(
            {"seed": "rad", "delta": 0.3, "w1s": 25.0, "c": c, "n": 32, "L": 48.0, "src": None}
        )
        if t in rows:
            ch = rows[t]["chunks"]
            ax[2].plot([q["iters"] for q in ch], [q["r_half"] for q in ch], label=f"c = {c:g}")
    ax[2].set_xlabel("iterations"), ax[2].set_ylabel("half-energy radius"), ax[2].set_title(
        "the drift read"
    )
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    out = os.path.join(HERE, "..", "plots", "m5_32_r23_1_cscan.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        smoke()
    elif mode == "run":
        run_pool(jobs_all(), int(sys.argv[2]) if len(sys.argv) > 2 else 12)
    elif mode == "extra":
        global OUT_JSON
        OUT_JSON = EXTRA_JSON
        run_pool(jobs_extra(), int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    elif mode == "extra2":
        OUT_JSON = EXTRA_JSON.replace("_extra.json", "_extra2.json")
        run_pool(jobs_extra2(), int(sys.argv[2]) if len(sys.argv) > 2 else 5)
    elif mode == "collect":
        collect()
    elif mode == "plot":
        plot()
    elif mode == "jobs":
        for j in jobs_all():
            print(job_tag(j), j["part"])


if __name__ == "__main__":
    main()
