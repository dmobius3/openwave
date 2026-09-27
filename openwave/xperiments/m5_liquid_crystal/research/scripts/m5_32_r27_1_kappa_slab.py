"""M5.32 R27-1: the strand under the smoothness term kappa |dM|^2 on the slab (the kappa arm,
the k^2 law on two boxes, the half-strand pair's E(d)) and the wall tension's continuum limit
(the stored wall states re-relaxed at h 1, 1/2, 1/4).

EQUATIONS FIRST
---------------
The slab is R26-2's: the deviation form M = M0 + delta D, M0 = diag(8, 1, delta, 0), nz 4
z-invariant layers, the x-y shell of depth 1.6 pinned at the seed, the FULL static sector
(M_0i = 0, the director row free: the sector the walls live in), the energy
    E[D] = delta^4 E_u[D] + V4[M0 + delta D] + kappa E_kappa[M0 + delta D],
E_u the certified commutator curvature, V4 by the telescoped brackets (w = W1 x 25), and the
smoothness term of R27-0 (the bond sum on the spatial block, its gradient checked there):
    E_kappa = h^3 sum_br w_br sum_i <d_i S, d_i S>.
The optimizer is R26-2's slot-preconditioned L-BFGS (the diagonal slots rescaled by their V4
stiffness, the variables the seven static-sector entries per free cell), the one that
resolves the wall class (R26-2: the plain variables stay on the strand, the preconditioned
ones find the walls), in chunks of 500 to the cap. THE GATE is the R26-2 audit's: the
pair-slot gradient in scaled units, max |dE/dM|_{pair block, free cells} / (delta^3 h^3),
under GATE_PAIR, with the last chunk's relative energy drop under GATE_DROP; a row at the cap
without both is FALLING; a line-search stall is STALL.
THE KAPPA UNIT: kappa_ref = T_bps / E_kappa(BPS seed, per unit length) per delta (R27-0's
stored table: 3.912e-3 at delta 0.3, 3.618e-5 at delta 0.03, 4.017e-6 at delta 0.01), the kappa
at which the smooth strand's smoothness energy equals its Bogomolny tension.
THE KAPPA ARM: delta 0.03 and 0.3, h 1, n 48, L 48; seeds: the BPS strand (m 1), the stored
wall state (R26-2's escape row at the same delta), the full strand (m 2), the half-strand pair
at d 6 h; kappa in {0, 0.03, 0.1, 0.3, 1} kappa_ref (40 rows). Under the ungauged term two
like half strands repel: E_pair - E_full = -8 pi kappa b0^2 ln(d / xi) (the XY energy with
J = 2 kappa b0^2); the pair rows at d 4, 6, 9 h at kappa 0.3 kappa_ref read the slope.
THE k^2 LAW: E_kappa(L 96) - E_kappa(L 48) = 4 pi m^2 b0^2 ln 2 per unit length on the RELAXED
rows at kappa 0.3 kappa_ref, m 1 and 2 (the two-box difference isolates the log coefficient;
a single box mixes in the core constant, R27-0 c: the m 2 to m 1 ratio on one box is 3.0).
THE WALL LADDER: the stored wall states (delta 0.3, 0.03, 0.01, h 1, T / T_bps 0.13 to 0.16 at
R26-2's cap) re-relaxed at kappa 0 at h 1 (n 48), interpolated (cubic, z-invariant) to h 1/2
(n 96) and h 1/4 (n 192) and re-relaxed there, each to the gate or the cap: a continuum
texture keeps its energy, a lattice artifact loses it under refinement (R26-2's audit:
the walls RE-FORM at h/2).
READS per row: E_quartic / T_bps (the curvature plus V4), E_kappa (kappa = 1 units) and
kappa E_kappa / T_bps, the wall bonds (tensor-phase jumps over 150 and 90 degrees between
free cells outside the melted cores), the wall components, the wall width, the walls' b / b0,
the director-row departure, the core b / b0, the cores' positions and distance, the outer
loop winding.

PRE-REGISTERED LABELS (the record's R27 PLANNING section)
--------------------------------------------------------
    WALL_TENSION_VANISHES   at every delta the re-relaxed wall state's E / T_bps falls by more
                            than 30 percent per halving from h 1 to 0.5 to 0.25 (both ratios
                            below 0.7) with the rows at the gate
    WALL_TENSION_FINITE     the h 0.5 and h 0.25 values agree within 10 percent and sit above
                            0.05 T_bps at every delta
    WALL_TENSION_UNRESOLVED otherwise (a row not at the gate, or the deltas disagree)
    KAPPA_RESTORES_FLOOR    at both deltas some kappa <= kappa_ref makes the wall seed end with
                            zero walls and within 1 percent of the smooth seed's end energy at
                            the same kappa, the smooth seed's quartic part within 3 percent of
                            T_bps
    KAPPA_WALLS_PERSIST     at kappa_ref the wall seed still ends more than 5 percent below the
                            smooth seed with walls present, at either delta
    KAPPA_UNRESOLVED        otherwise
    K2_LAW_CONFIRMED        the two-box difference ratio m 2 / m 1 within 10 percent of 4 at both
                            deltas, the m 1 value within 10 percent of 4 pi b0^2 ln 2
    K2_LAW_REFUTED          off by more than 25 percent at both deltas with the rows at the gate
    K2_LAW_UNRESOLVED       otherwise
    the pair's E(d) slope is reported against -8 pi kappa b0^2, never labeled.

Modes: smoke | run [workers] | collect | jobs. Output: data/m5_32_r27_1_kappa_slab.json,
arrays in data/m5_32_r27_1/ (local, kept). Regenerate: run 10 about 1.5 h (the n 192 rows the
long pole); collect seconds; smoke a minute.
"""

import importlib.util
import json
import multiprocessing as mp
import os
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r27_1_kappa_slab.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r27_1")
R27_0_JSON = os.path.join(DATA, "m5_32_r27_0_form.json")
R26_2_DIR = os.path.join(DATA, "m5_32_r26_2")
NZ = 4
W1S = 25.0
KAPPA_FRACS = (0.0, 0.03, 0.1, 0.3, 1.0)
K2_FRAC = 0.3
CAP = {16: 20, 48: 4000, 96: 3000, 192: 1500}
LBFGS_CHUNK = 500
GATE_PAIR = 1e-3
GATE_DROP = 1e-5
KICK_EXCESS = 0.0  # the seeds are relaxed as given (a kick would move the wall seed off its class)
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


F = _load("m5_32_r27_0_form", "m5_32_r27_0_form.py")
R2, R25_1, B3, W1 = F.R2, F.R25_1, F.B3, F.W1


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def kappa_refs():
    J = json.load(open(R27_0_JSON))
    return {k[1:]: v["kappa_ref"] for k, v in J["stored"]["slab"].items()}


# ============================================================================
# the objective and the descent
# ============================================================================
class ScaledK:
    """R26-2's scaled objective plus kappa E_kappa, the full static sector."""

    def __init__(self, cfg, delta, w, free, T_ref, kappa):
        self.cfg, self.delta, self.w, self.free, self.T_ref, self.kappa = (
            cfg,
            delta,
            w,
            free,
            T_ref,
            kappa,
        )
        self.fr = free[..., None, None].astype(float)
        self.M0 = R2.M_vac(delta)

    def parts(self, D):
        eu, ev, _ = R2.energy_grad_M(D, self.cfg, self.delta, self.w, need_grad=False)
        ek = F.e_kappa(self.M0 + self.delta * D, self.cfg)
        return eu, ev, ek

    def fun_grad(self, D):
        eu, ev, G = R2.energy_grad_M(D, self.cfg, self.delta, self.w)
        M = self.M0 + self.delta * D
        ek = F.e_kappa(M, self.cfg)
        if self.kappa != 0.0:
            G = R25_1.static_sector(G + self.kappa * F.grad_kappa(M, self.cfg))
        return (eu + ev + self.kappa * ek) / self.T_ref, (self.delta / self.T_ref) * G * self.fr

    def grad_M(self, D):
        eu, ev, G = R2.energy_grad_M(D, self.cfg, self.delta, self.w)
        M = self.M0 + self.delta * D
        if self.kappa != 0.0:
            G = R25_1.static_sector(G + self.kappa * F.grad_kappa(M, self.cfg))
        return G

    def energy(self, D):
        eu, ev, ek = self.parts(D)
        return eu + ev + self.kappa * ek


def gmax_pair(S, D):
    """the R26-2 audit's measure: the pair-block gradient in scaled units."""
    G = S.grad_M(D)
    Gf = G[S.free]
    return float(np.max(np.abs(Gf[:, 1:3, 1:3]))) / (S.delta**3 * S.cfg["h"] ** 3)


def lbfgs(S, D0, max_iter, sc, chunk=LBFGS_CHUNK, tag=""):
    from scipy.optimize import minimize

    D = D0.copy()
    done, chunks, verdict = 0, [], "FALLING"
    f_prev = S.fun_grad(D)[0]
    f_start = f_prev
    gp = gmax_pair(S, D)
    chunks.append({"iters": 0, "f": f_prev, "gmax_pair": gp, "note": "the seed"})
    while done < max_iter:
        base = D.copy()

        def fun(x, base=base):
            Dx = R2.unpack(x, base, S.free, sc)
            f, g = S.fun_grad(Dx)
            return f, R2.grad_pack(g, S.free, sc)

        k = min(chunk, max_iter - done)
        res = minimize(
            fun,
            R2.pack(D, S.free, sc),
            jac=True,
            method="L-BFGS-B",
            options={"maxcor": 20, "maxiter": k, "maxfun": 3 * k, "gtol": 1e-16, "ftol": 1e-18},
        )
        D = R2.unpack(np.asarray(res.x), base, S.free, sc)
        its = int(res.nit)
        done += max(1, its)
        f = S.fun_grad(D)[0]
        gp = gmax_pair(S, D)
        drop = (f_prev - f) / max(abs(f), 1e-300)
        chunks.append(
            {
                "iters": done,
                "f": f,
                "gmax_pair": gp,
                "drop_rel": drop,
                "chunk_iters": its,
                "scipy": str(res.message),
            }
        )
        log(f"{tag} its {done} f {f:.8e} gpair {gp:.2e} drop {drop:.1e}")
        if gp < GATE_PAIR and 0.0 <= drop < GATE_DROP:
            verdict = "AT_GATE"
            break
        if its < 3:
            verdict = "STALL" if gp < 10.0 * GATE_PAIR else "LINE_SEARCH_STALL"
            break
        f_prev = f
    return D, {
        "f_start": f_start,
        "f_end": f,
        "gmax_pair_end": gp,
        "iters": done,
        "verdict": verdict,
        "chunks": chunks,
    }


# ============================================================================
# the seeds
# ============================================================================
def wall_state(delta):
    Z = np.load(os.path.join(R26_2_DIR, f"escape_d{delta:g}_w25_n48_L48_long_amp0.npz"))
    return Z["D"]


def refine_D(D, factor, order=3):
    """cubic interpolation of the z-invariant deviation to a grid `factor` times finer."""
    from scipy.ndimage import zoom

    n = D.shape[0]
    Dz = D[:, :, 0]
    Df = np.zeros((n * factor, n * factor, NZ, 4, 4))
    for a in range(4):
        for c in range(a, 4):
            v = zoom(Dz[:, :, a, c], factor, order=order, mode="nearest")
            Df[:, :, :, a, c] = v[:, :, None]
            Df[:, :, :, c, a] = v[:, :, None]
    return Df


def seed_of(j, cfg):
    n, h, delta, w = j["n"], cfg["h"], j["delta"], W1 * W1S
    s = j["seed"]
    if s == "bps_m1":
        return R2.bps_D(n, NZ, h, delta, w), "the BPS strand, m 1"
    if s == "bps_m2":
        return (
            R25_1.bps_field(n, NZ, h, delta, w, m=2) - R2.M_vac(delta)
        ) / delta, "the full strand, m 2"
    if s.startswith("pair_d"):
        d = float(s[6:]) * h
        return F.pair_D(n, NZ, h, delta, w, d), f"two like half strands at d {d:g}"
    if s == "wall":
        D = wall_state(delta)
        factor = n // D.shape[0]
        if factor == 1:
            return D, "the stored wall state (R26-2 escape row)"
        return refine_D(D, factor), f"the stored wall state interpolated x{factor} (cubic)"
    raise ValueError(s)


# ============================================================================
# the jobs
# ============================================================================
def job_tag(j):
    return f"{j['seed']}_d{j['delta']:g}_n{j['n']}_L{j['L']:g}_k{j['kappa_frac']:g}"


def jobs_all():
    jobs = []
    for delta in (0.03, 0.3):
        for seed in ("bps_m1", "wall", "bps_m2", "pair_d6"):
            for kf in KAPPA_FRACS:
                jobs.append(dict(arm="kappa", seed=seed, delta=delta, n=48, L=48.0, kappa_frac=kf))
        for seed in ("bps_m1", "bps_m2"):
            jobs.append(dict(arm="k2", seed=seed, delta=delta, n=96, L=96.0, kappa_frac=K2_FRAC))
        for seed in ("pair_d4", "pair_d9"):
            jobs.append(dict(arm="pair", seed=seed, delta=delta, n=48, L=48.0, kappa_frac=K2_FRAC))
    for delta in (0.3, 0.03, 0.01):
        for n in (48, 96, 192):
            j = dict(arm="wall", seed="wall", delta=delta, n=n, L=48.0, kappa_frac=0.0)
            if job_tag(j) not in {job_tag(q) for q in jobs}:
                jobs.append(j)
    return jobs


def reads(D, cfg, delta, w, free, S):
    M = R2.M_vac(delta) + delta * D
    h = cfg["h"]
    eu, ev, ek = S.parts(D)
    Tb = R2.T_bps(delta, w)
    per = NZ * h
    out = {
        "E_u_per_len": eu / per,
        "V4_per_len": ev / per,
        "E_quartic_per_len": (eu + ev) / per,
        "E_quartic_over_T_bps": (eu + ev) / per / Tb,
        "E_kappa_per_len": ek / per,
        "kappaE_kappa_over_T_bps": S.kappa * ek / per / Tb,
        "E_total_over_T_bps": (eu + ev + S.kappa * ek) / per / Tb,
        "walls": F.wall_reads(M, delta, cfg, free),
        "departure": R25_1.departure(M, delta),
        "cores": F.core_positions(M, delta, h, free),
        "winding_outer_rho12": F.loop_winding(M, 0.0, 0.0, min(12.0, 0.5 * cfg["L"] - 4.0 * h), h),
        "winding_outer_rho20": F.loop_winding(M, 0.0, 0.0, min(20.0, 0.5 * cfg["L"] - 4.0 * h), h),
        "b_over_b0_min": F.wall_reads(M, delta, cfg, free)["b_over_b0_min_free"],
    }
    c = out["cores"]
    out["core_distance"] = (
        float(np.hypot(c[0][0] - c[1][0], c[0][1] - c[1][1])) if len(c) == 2 else None
    )
    out["core_count"] = len(c)
    return out


def run_job(j):
    t0 = time.time()
    tag = job_tag(j)
    n, L, delta, w = j["n"], j["L"], j["delta"], W1 * W1S
    cfg = R2.slab_cfg(n, L, delta)
    h = cfg["h"]
    row = dict(j, tag=tag, h=h, w=w)
    try:
        krefs = kappa_refs()
        kref = krefs[f"{delta:g}"]
        kappa = j["kappa_frac"] * kref
        free = R2.free_mask_slab(n, NZ, h)
        Tb = R2.T_bps(delta, w)
        T_ref = Tb * NZ * h
        D0, note = seed_of(j, cfg)
        row.update(
            kappa_ref=kref,
            kappa=kappa,
            T_bps=Tb,
            T_ref=T_ref,
            seed_note=note,
            cap=CAP[n],
            gate_pair=GATE_PAIR,
            gate_drop=GATE_DROP,
        )
        S = ScaledK(cfg, delta, w, free, T_ref, kappa)
        row["seed_reads"] = reads(D0, cfg, delta, w, free, S)
        sc, ks = R2.slot_scales(cfg, delta, w, T_ref)
        row["preconditioner"] = {
            "scales_e_u_v_t": [float(sc[0]), float(sc[1]), float(sc[4]), float(sc[6])],
            **ks,
        }
        os.makedirs(OUT_NPZ, exist_ok=True)
        D, lb = lbfgs(S, D0, CAP[n], sc, tag=tag)
        row["lbfgs"] = lb
        row["verdict"] = lb["verdict"]
        row["end_reads"] = reads(D, cfg, delta, w, free, S)
        row["status"] = "OK" if np.all(np.isfinite(D)) else "NONFINITE"
        np.savez_compressed(
            os.path.join(OUT_NPZ, tag + ".npz"), D=D.astype(np.float64), delta=delta, kappa=kappa
        )
    except Exception as e:  # noqa: BLE001
        import traceback

        row.update(status="FAILED", stop=repr(e), traceback=traceback.format_exc())
    row["wall_s"] = round(time.time() - t0, 1)
    er = row.get("end_reads") or {}
    log(
        f"DONE {tag} {row.get('status')} {row.get('verdict')} Eq/Tb {er.get('E_quartic_over_T_bps')} "
        f"Etot/Tb {er.get('E_total_over_T_bps')} walls {(er.get('walls') or {}).get('wall_bonds_150')} wall {row['wall_s']}"
    )
    return row


def continue_job(tag, extra):
    """RUN-TIME ADDITION (15:05 UTC): continue a saved row from its end field for `extra` more
    iterations (the h 1/4 wall rows were far from the gate at their 1500-iteration cap); the
    row's chunks are appended, its end reads and verdict replaced, the continuation logged."""
    t0 = time.time()
    J = load_json()
    r = J["rows"][tag]
    n, L, delta, w = r["n"], r["L"], r["delta"], W1 * W1S
    cfg = R2.slab_cfg(n, L, delta)
    h = cfg["h"]
    free = R2.free_mask_slab(n, NZ, h)
    S = ScaledK(cfg, delta, w, free, r["T_ref"], r["kappa"])
    D0 = np.load(os.path.join(OUT_NPZ, tag + ".npz"))["D"]
    sc, _ = R2.slot_scales(cfg, delta, w, r["T_ref"])
    D, lb = lbfgs(S, D0, extra, sc, tag=tag + "_cont")
    prev = r["lbfgs"]
    off = prev["iters"]
    for c in lb["chunks"][1:]:
        c["iters"] += off
    r["lbfgs"] = {
        "f_start": prev["f_start"],
        "f_end": lb["f_end"],
        "gmax_pair_end": lb["gmax_pair_end"],
        "iters": off + lb["iters"],
        "verdict": lb["verdict"],
        "chunks": prev["chunks"] + lb["chunks"][1:],
    }
    r["verdict"] = lb["verdict"]
    r["continued"] = (r.get("continued") or []) + [
        {"from_iters": off, "extra": extra, "wall_s": round(time.time() - t0, 1)}
    ]
    r["end_reads_before_continuation"] = r.get("end_reads_before_continuation") or r["end_reads"]
    r["end_reads"] = reads(D, cfg, delta, w, free, S)
    for k in (
        "gmax_all_scaled",
        "gmax_pair_scaled",
        "gmax_director_row_scaled",
        "gmax_diag_scaled",
        "walls_90_end",
    ):
        r.pop(k, None)
    np.savez_compressed(
        os.path.join(OUT_NPZ, tag + ".npz"), D=D.astype(np.float64), delta=delta, kappa=r["kappa"]
    )
    Jn = load_json()
    Jn["rows"][tag] = r
    save_json(Jn)
    log(
        f"DONE {tag} continued to {r['lbfgs']['iters']} {r['verdict']} Eq/Tb {r['end_reads']['E_quartic_over_T_bps']}"
    )
    return r


def continue_pool(tags, extra, workers):
    log(f"continue pool: {len(tags)} rows, {workers} workers, +{extra} iterations")
    with ProcessPoolExecutor(max_workers=int(workers), mp_context=mp.get_context("spawn")) as ex:
        futs = [ex.submit(continue_job, t, extra) for t in tags]
        for fut in as_completed(futs):
            fut.result()
    log("pool done")


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R27-1", "rows": {}, "collect": {}}


def save_json(J):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1, default=str)
    os.replace(tmp, OUT_JSON)


def run_pool(jobs, workers):
    workers = min(int(workers), 12)
    rows = load_json()["rows"]
    pending = [j for j in jobs if rows.get(job_tag(j), {}).get("status") != "OK"]
    pending.sort(key=lambda j: -j["n"])
    log(f"pool: {len(pending)} jobs, {workers} workers")
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        futs = [ex.submit(run_job, j) for j in pending]
        for fut in as_completed(futs):
            row = fut.result()
            Jn = load_json()
            Jn["rows"][row["tag"]] = row
            save_json(Jn)
    log("pool done")


# ============================================================================
# collect
# ============================================================================
def _row(rows, seed, delta, n, L, kf):
    return rows.get(job_tag(dict(seed=seed, delta=delta, n=n, L=L, kappa_frac=kf)))


def post_reads(rows):
    """RUN-TIME ADDITION (14:10 UTC deviation): the gradient over ALL static slots in the audit's
    scaled units, read from the saved end field, beside the pair-slot measure the gate used."""
    w = W1 * W1S
    for tag, r in rows.items():
        if r.get("status") != "OK" or "gmax_all_scaled" in r:
            continue
        f = os.path.join(OUT_NPZ, tag + ".npz")
        if not os.path.exists(f):
            continue
        Z = np.load(f)
        D, delta = Z["D"], float(Z["delta"])
        cfg = R2.slab_cfg(r["n"], r["L"], delta)
        free = R2.free_mask_slab(r["n"], NZ, cfg["h"])
        S = ScaledK(cfg, delta, w, free, r["T_ref"], r["kappa"])
        G = S.grad_M(D)[free]
        sc = delta**3 * cfg["h"] ** 3
        r["gmax_all_scaled"] = float(np.max(np.abs(G))) / sc
        r["gmax_pair_scaled"] = float(np.max(np.abs(G[:, 1:3, 1:3]))) / sc
        r["gmax_director_row_scaled"] = float(np.max(np.abs(G[:, 3, 1:3]))) / sc
        r["gmax_diag_scaled"] = (
            float(max(np.max(np.abs(G[:, 0, 0])), np.max(np.abs(G[:, 3, 3])))) / sc
        )
        r["walls_90_end"] = r["end_reads"]["walls"]["wall_bonds_90"]
    return rows


def collect(rows=None):
    J = load_json()
    rows = J["rows"] if rows is None else rows
    if rows is J["rows"]:
        post_reads(rows)
    w = W1 * W1S
    out = {}

    def ok(r):
        return r is not None and r.get("status") == "OK"

    def gate(r):
        # RUN-TIME DEFINITION (14:50 UTC): a STALL row (scipy's line search stopping) whose
        # gradient over ALL static slots is under 1e-3 in the audit's scaled units is converged
        # (the STALL verdict was this rung's own addition, not pre-registered); such rows count
        # as at the gate, flagged STALL_CONVERGED in the tables
        if not ok(r):
            return False
        if r["verdict"] == "AT_GATE":
            return True
        return r["verdict"] == "STALL" and (r.get("gmax_all_scaled") or 1.0) < 1e-3

    def verdict_of(r):
        if r["verdict"] == "STALL" and gate(r):
            return "STALL_CONVERGED"
        return r["verdict"]

    # ---- the wall ladder ----
    lad = {}
    for delta in (0.3, 0.03, 0.01):
        rec = {}
        for n in (48, 96, 192):
            r = _row(rows, "wall", delta, n, 48.0, 0.0)
            rec[f"h{48.0 / n:g}"] = (
                {
                    "E_over_T_bps": r["end_reads"]["E_quartic_over_T_bps"],
                    "E_seed_over_T_bps": r["seed_reads"]["E_quartic_over_T_bps"],
                    "verdict": verdict_of(r),
                    "iters": r["lbfgs"]["iters"],
                    "wall_bonds_150": r["end_reads"]["walls"]["wall_bonds_150"],
                    "wall_components": r["end_reads"]["walls"]["wall_components"],
                    "wall_width_h": r["end_reads"]["walls"]["wall_width_h_median"],
                    "wall_bonds_90": r["end_reads"]["walls"]["wall_bonds_90"],
                    "director_row": r["end_reads"]["departure"]["director_row"],
                    "gmax_all_scaled": r.get("gmax_all_scaled"),
                    "gmax_pair_scaled_end": r["lbfgs"]["gmax_pair_end"],
                    "last_drop_rel": r["lbfgs"]["chunks"][-1].get("drop_rel"),
                }
                if ok(r)
                else None
            )
        vals = [rec[f"h{x:g}"] for x in (1.0, 0.5, 0.25)]
        if all(v is not None for v in vals):
            e1, e2, e3 = (v["E_over_T_bps"] for v in vals)
            rec["ratio_h05_over_h1"] = e2 / e1
            rec["ratio_h025_over_h05"] = e3 / e2
            rec["h05_h025_agree_10pct"] = abs(e2 - e3) / max(abs(e3), 1e-300) < 0.10
            rec["all_at_gate"] = all(v["verdict"] in ("AT_GATE", "STALL_CONVERGED") for v in vals)
        lad[f"d{delta:g}"] = rec
    complete = all("ratio_h05_over_h1" in lad[k] for k in lad)
    if complete:
        van = all(
            lad[k]["ratio_h05_over_h1"] < 0.7
            and lad[k]["ratio_h025_over_h05"] < 0.7
            and lad[k]["all_at_gate"]
            for k in lad
        )
        fin = all(
            lad[k]["h05_h025_agree_10pct"] and lad[k]["h0.25"]["E_over_T_bps"] > 0.05 for k in lad
        )
        lad["label"] = (
            "WALL_TENSION_VANISHES"
            if van
            else ("WALL_TENSION_FINITE" if fin else "WALL_TENSION_UNRESOLVED")
        )
        lad["gate_note"] = (
            "FINITE by the rule does not require the gate; the verdicts are beside each value"
        )
    else:
        lad["label"] = "WALL_TENSION_UNRESOLVED (rows missing)"
    out["wall_ladder"] = lad

    # ---- the kappa arm ----
    arm = {}
    restores_per_delta = {}
    persist_any = False
    for delta in (0.03, 0.3):
        tab = {}
        for kf in KAPPA_FRACS:
            ent = {}
            for seed in ("bps_m1", "wall", "bps_m2", "pair_d6"):
                r = _row(rows, seed, delta, 48, 48.0, kf)
                ent[seed] = (
                    {
                        "E_total_over_T_bps": r["end_reads"]["E_total_over_T_bps"],
                        "E_quartic_over_T_bps": r["end_reads"]["E_quartic_over_T_bps"],
                        "kappaE_kappa_over_T_bps": r["end_reads"]["kappaE_kappa_over_T_bps"],
                        "E_kappa_per_len": r["end_reads"]["E_kappa_per_len"],
                        "wall_bonds_150": r["end_reads"]["walls"]["wall_bonds_150"],
                        "wall_components": r["end_reads"]["walls"]["wall_components"],
                        "verdict": r["verdict"],
                        "iters": r["lbfgs"]["iters"],
                        "core_count": r["end_reads"]["core_count"],
                        "core_distance": r["end_reads"]["core_distance"],
                        "director_row": r["end_reads"]["departure"]["director_row"],
                        "b_over_b0_min": r["end_reads"]["b_over_b0_min"],
                        "wall_bonds_90": r["end_reads"]["walls"]["wall_bonds_90"],
                        "wall_width_h": r["end_reads"]["walls"]["wall_width_h_median"],
                        "gmax_all_scaled": r.get("gmax_all_scaled"),
                        "gmax_pair_scaled_end": r["lbfgs"]["gmax_pair_end"],
                    }
                    if ok(r)
                    else None
                )
            for sd in ("bps_m1", "wall", "bps_m2", "pair_d6"):
                r = _row(rows, sd, delta, 48, 48.0, kf)
                if ent.get(sd) and ok(r):
                    ent[sd]["verdict"] = verdict_of(r)
            if ent["bps_m1"] and ent["wall"]:
                ent["wall_minus_smooth_rel"] = (
                    ent["wall"]["E_total_over_T_bps"] - ent["bps_m1"]["E_total_over_T_bps"]
                ) / max(abs(ent["bps_m1"]["E_total_over_T_bps"]), 1e-300)
                ent["wall_seed_ends_without_walls"] = ent["wall"]["wall_bonds_150"] == 0
                ent["smooth_quartic_within_3pct_of_T_bps"] = (
                    abs(ent["bps_m1"]["E_quartic_over_T_bps"] - 1.0) < 0.03
                )
            if ent["bps_m2"] and ent["pair_d6"]:
                ent["pair_minus_full_over_T_bps"] = (
                    ent["pair_d6"]["E_total_over_T_bps"] - ent["bps_m2"]["E_total_over_T_bps"]
                )
            tab[f"k{kf:g}"] = ent
        restores = any(
            tab[f"k{kf:g}"].get("wall_seed_ends_without_walls")
            and abs(tab[f"k{kf:g}"].get("wall_minus_smooth_rel", 9)) < 0.01
            and tab[f"k{kf:g}"].get("smooth_quartic_within_3pct_of_T_bps")
            for kf in KAPPA_FRACS
            if kf <= 1.0
        )
        restores_per_delta[f"d{delta:g}"] = restores
        # PRE-REGISTRATION FLAW (found at collect, 14:50 UTC): the clause "the smooth seed's
        # quartic part within 3 percent of T_bps" cannot hold, since the strand's kappa-0 value is
        # T_slaved3 = 0.85 T_bps (delta 0.03) and 0.60 T_bps (delta 0.3) (R25, R26); the amended
        # clause below takes the smooth seed's wall-free quartic at the smallest kappa as the
        # strand's reference and is REPORTED beside the letter, never substituted for it
        ref = None
        for kf in KAPPA_FRACS[1:]:
            q = tab[f"k{kf:g}"].get("bps_m1")
            if q and q["wall_bonds_90"] == 0:
                ref = (kf, q["E_quartic_over_T_bps"])
                break
        tab["strand_reference_kappa_frac_and_quartic"] = ref
        for kf in KAPPA_FRACS:
            e = tab[f"k{kf:g}"]
            if ref and e.get("bps_m1"):
                e["smooth_quartic_over_strand_reference"] = (
                    e["bps_m1"]["E_quartic_over_T_bps"] / ref[1]
                )
                e["restores_amended"] = bool(
                    e.get("wall_seed_ends_without_walls")
                    and e["wall"]["wall_bonds_90"] == 0
                    and abs(e.get("wall_minus_smooth_rel", 9)) < 0.01
                    and abs(e["smooth_quartic_over_strand_reference"] - 1.0) < 0.03
                )
        tab["restores_amended_at_kappa_fracs"] = [
            kf for kf in KAPPA_FRACS if tab[f"k{kf:g}"].get("restores_amended")
        ]
        tab["wall_seed_identical_to_smooth_at_kappa_fracs"] = [
            kf for kf in KAPPA_FRACS if abs(tab[f"k{kf:g}"].get("wall_minus_smooth_rel", 9)) < 1e-6
        ]
        e1 = tab["k1"]
        if e1.get("wall") and e1.get("bps_m1"):
            if e1["wall_minus_smooth_rel"] < -0.05 and e1["wall"]["wall_bonds_150"] > 0:
                persist_any = True
        arm[f"d{delta:g}"] = tab
    have_all = all(
        _row(rows, s, d, 48, 48.0, kf) is not None
        for d in (0.03, 0.3)
        for s in ("bps_m1", "wall")
        for kf in KAPPA_FRACS
    )
    if have_all:
        if all(restores_per_delta.values()):
            arm["label"] = "KAPPA_RESTORES_FLOOR"
        elif persist_any:
            arm["label"] = "KAPPA_WALLS_PERSIST"
        else:
            arm["label"] = "KAPPA_UNRESOLVED"
    else:
        arm["label"] = "KAPPA_UNRESOLVED (rows missing)"
    arm["restores_per_delta"] = restores_per_delta
    out["kappa_arm"] = arm

    # ---- the k^2 law ----
    k2 = {}
    conf, refut = [], []
    for delta in (0.03, 0.3):
        b0 = delta / 2.0
        rec = {}
        diffs = {}
        gated = True
        for m in (1, 2):
            r48 = _row(rows, f"bps_m{m}", delta, 48, 48.0, K2_FRAC)
            r96 = _row(rows, f"bps_m{m}", delta, 96, 96.0, K2_FRAC)
            if ok(r48) and ok(r96):
                diffs[m] = (
                    r96["end_reads"]["E_kappa_per_len"] - r48["end_reads"]["E_kappa_per_len"]
                )
                gated = gated and gate(r48) and gate(r96)
                rec[f"m{m}"] = {
                    "E_kappa_L48": r48["end_reads"]["E_kappa_per_len"],
                    "E_kappa_L96": r96["end_reads"]["E_kappa_per_len"],
                    "difference": diffs[m],
                    "verdicts": [verdict_of(r48), verdict_of(r96)],
                    "core_count_end_L48_L96": [
                        r48["end_reads"]["core_count"],
                        r96["end_reads"]["core_count"],
                    ],
                    "core_distance_end_L48_L96": [
                        r48["end_reads"]["core_distance"],
                        r96["end_reads"]["core_distance"],
                    ],
                    "seed_difference": r96["seed_reads"]["E_kappa_per_len"]
                    - r48["seed_reads"]["E_kappa_per_len"],
                }
        if 1 in diffs and 2 in diffs:
            pred1 = 4.0 * np.pi * b0**2 * np.log(2.0)
            rec["ratio_m2_over_m1"] = diffs[2] / diffs[1]
            rec["m1_over_prediction"] = diffs[1] / pred1
            rec["prediction_m1"] = pred1
            rec["all_at_gate"] = gated
            c = (
                abs(rec["ratio_m2_over_m1"] - 4.0) < 0.4
                and abs(rec["m1_over_prediction"] - 1.0) < 0.10
            )
            rf = (
                abs(rec["ratio_m2_over_m1"] - 4.0) > 1.0
                or abs(rec["m1_over_prediction"] - 1.0) > 0.25
            ) and gated
            conf.append(c)
            refut.append(rf)
        k2[f"d{delta:g}"] = rec
    if len(conf) == 2:
        k2["label"] = (
            "K2_LAW_CONFIRMED"
            if all(conf)
            else ("K2_LAW_REFUTED" if all(refut) else "K2_LAW_UNRESOLVED")
        )
    else:
        k2["label"] = "K2_LAW_UNRESOLVED (rows missing)"
    out["k2_law"] = k2

    # ---- the pair E(d) ----
    pr = {}
    for delta in (0.03, 0.3):
        b0 = delta / 2.0
        kref = kappa_refs()[f"{delta:g}"]
        kappa = K2_FRAC * kref
        Tb = R2.T_bps(delta, w)
        pts = []
        for d in (4, 6, 9):
            r = _row(rows, f"pair_d{d}", delta, 48, 48.0, K2_FRAC)
            if ok(r):
                pts.append(
                    {
                        "d_h": d,
                        "E_total_per_len": r["end_reads"]["E_total_over_T_bps"] * Tb,
                        "E_total_over_T_bps": r["end_reads"]["E_total_over_T_bps"],
                        "core_distance_end": r["end_reads"]["core_distance"],
                        "core_count_end": r["end_reads"]["core_count"],
                        "verdict": r["verdict"],
                    }
                )
        rec = {
            "points": pts,
            "prediction_slope_dE_dlnd": -8.0 * np.pi * kappa * b0**2,
            "kappa": kappa,
        }
        good = [q for q in pts if q["core_distance_end"] and q["core_count_end"] == 2]
        if len(good) >= 2:
            x = np.log([q["core_distance_end"] for q in good])
            y = [q["E_total_per_len"] for q in good]
            rec["slope_dE_dlnd_measured"] = float(np.polyfit(x, y, 1)[0])
            rec["slope_over_prediction"] = (
                rec["slope_dE_dlnd_measured"] / rec["prediction_slope_dE_dlnd"]
            )
        full = _row(rows, "bps_m2", delta, 48, 48.0, K2_FRAC)
        if ok(full):
            rec["E_full_strand_per_len"] = full["end_reads"]["E_total_over_T_bps"] * Tb
            rec["E_full_strand_over_T_bps"] = full["end_reads"]["E_total_over_T_bps"]
        pr[f"d{delta:g}"] = rec
    pr["note"] = (
        "reported against -8 pi kappa b0^2, never labeled (the pinned shell differs per d)"
    )
    out["pair_E_of_d"] = pr

    out["labels"] = {
        "wall_ladder": lad["label"],
        "kappa_arm": arm["label"],
        "k2_law": k2["label"],
    }
    out["rows_total"] = len(rows)
    out["rows_ok"] = sum(1 for r in rows.values() if r.get("status") == "OK")
    out["rows_at_gate"] = sum(1 for r in rows.values() if r.get("verdict") == "AT_GATE")
    J["collect"] = out
    save_json(J)
    return out


# ============================================================================
# smoke
# ============================================================================
def smoke():
    global OUT_JSON, OUT_NPZ
    OUT_JSON = OUT_JSON.replace(".json", "_smoke.json")
    OUT_NPZ = OUT_NPZ + "_smoke"
    res = {"rows": {}}
    for seed in ("bps_m1", "bps_m2", "pair_d4", "wall"):
        j = dict(
            arm="smoke",
            seed=seed,
            delta=0.03,
            n=48 if seed == "wall" else 16,
            L=48.0,
            kappa_frac=0.3,
        )
        r = run_job(j)
        res["rows"][r["tag"]] = {
            k: r.get(k) for k in ("status", "verdict", "seed_note", "stop", "traceback")
        }
        res["rows"][r["tag"]]["end"] = (
            {
                k: r["end_reads"][k]
                for k in ("E_quartic_over_T_bps", "E_total_over_T_bps", "core_count")
            }
            if r.get("end_reads")
            else None
        )
        # the wall row: n 48 would take minutes, cap it to one chunk of 20
    # a synthetic collect: rows with the label conditions met
    rows = {}
    for delta in (0.3, 0.03, 0.01):
        for n, e in ((48, 0.14), (96, 0.09), (192, 0.05)):
            j = dict(seed="wall", delta=delta, n=n, L=48.0, kappa_frac=0.0)
            er = {
                "E_quartic_over_T_bps": e,
                "E_total_over_T_bps": e,
                "kappaE_kappa_over_T_bps": 0.0,
                "E_kappa_per_len": 0.1,
                "walls": {"wall_bonds_150": 3, "wall_components": 1, "wall_width_h_median": 2.0},
                "departure": {"director_row": 0.1},
                "core_count": 0,
                "core_distance": None,
                "b_over_b0_min": 0.9,
            }
            rows[job_tag(j)] = dict(
                j,
                status="OK",
                verdict="AT_GATE",
                lbfgs={"iters": 1},
                seed_reads=dict(er),
                end_reads=er,
            )
    c = collect(rows)
    res["synthetic_wall_label"] = c["labels"]["wall_ladder"]
    res["PASS"] = (
        all(v["status"] == "OK" for v in res["rows"].values())
        and res["synthetic_wall_label"] == "WALL_TENSION_VANISHES"
    )
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "rows"}, indent=1))
    for t, v in res["rows"].items():
        print(t, v["status"], v["verdict"], v["end"], (v.get("stop") or "")[:200])
    return res


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        CAP[48] = 20
        smoke()
    elif mode == "run":
        run_pool(jobs_all(), sys.argv[2] if len(sys.argv) > 2 else 10)
        print(json.dumps(collect()["labels"], indent=1))
    elif mode == "collect":
        print(json.dumps(collect(), indent=1, default=str))
    elif mode == "continue_wall":
        extra = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
        tags = [f"wall_d{d:g}_n192_L48_k0" for d in (0.3, 0.03, 0.01)]
        continue_pool(tags, extra, 3)
        print(json.dumps(collect()["labels"], indent=1))
    elif mode == "jobs":
        for j in jobs_all():
            print(job_tag(j))
        print(len(jobs_all()), "jobs")
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
