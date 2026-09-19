"""M5.32 R22-1s: the R22-1 rows continued with the time entry slaved (added at
EXECUTE on the R22-1 audit, deviation 11).

EQUATIONS FIRST
---------------
The R22-1 static energy, E[M] = certified quartic + V4 with the targets
(-g, 1, delta, delta). On a block-diagonal field (M_0i = 0) the 1 x 1 time
block drops out of every commutator, so M_00 enters E through V4 alone:
    V4 = w sum_p (u^p + s_p - C_p)^2,   u = N_00 = -M_00,  s_p = tr S^p,
S the 3 x 3 spatial block. The V4 Hessian along M_00 is 8.46e6 w against 3.09 w
and 35.0 w for the two other massive directions (R22-0 (c)), and the R22-1
audit (m5_32_r22_1_audit.py, Q2.3, Q3, Q4.2) showed that this one stiff entry
rate-limits the 10-entry L-BFGS polish by a factor 44 to 53 and that the
max |G| gate was read on it.
The reduced functional: E_red[S] = min over M_00 of E[S, M_00], per cell, by
Newton on the quartic-in-u polynomial (warm-started), the pinned cells held.
By the envelope theorem dE_red / dS is the spatial block of the production
gradient at the slaved M_00. The minima of E_red are the minima of E on the
block-diagonal sector, which every R22-1 field is in (M_0i = 0 exactly, and
the static gradient has no time-row component).

Descent: scipy L-BFGS on the 6 spatial entries per free cell, from the R22-1
end fields, in resumable chunks of CHUNK iterations. Gate: max |G_spatial| on
the free cells < 1e-3 measured on the accepted iterate INSIDE a chunk, never
on the first iterate of a fresh chunk, AND the energy drop of the last full
chunk under 1e-4; the slaving residual max |G_00| is recorded beside it.

Reads per chunk: E and its parts, the virial E_curv / V, the half-energy
radius, the surface-oriented solid-angle degree of the top eigenvector on the
r 6 / 9 / 12 cubes with the conflicting-link count (R21's reader), the sorted
eigenvalues on the shells.

Smoke gates: (1) max |G_00| after slaving under 1e-9 of its unslaved value
scale; (2) E at the slaved field not above E before slaving; (3) the reduced
gradient against a central difference of E_red in random directions.

Modes: smoke | run [workers] | collect
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
OUT_JSON = os.path.join(DATA, "m5_32_r22_1s_slaved.json")
SRC_NPZ = os.path.join(DATA, "m5_32_r22_1")
OUT_NPZ = os.path.join(DATA, "m5_32_r22_1s")
G = 8.0
CHUNK = 250
MAX_ITER = 3000
GATE, GATE_DROP = 1e-3, 1e-4
IU3 = np.triu_indices(3)
OFF3 = np.where(IU3[0] != IU3[1], 2.0, 1.0)
SHELL_R = (1.5, 3.0, 4.5, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0)
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


R21 = _load("m5_32_r21_1_runs", "m5_32_r21_1_runs.py")
R20, B3, R0 = R21.R20, R21.B3, R21.R0
W1 = B3.W1


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def ups_guard():
    try:
        sys.path.insert(0, os.path.expanduser("~/Library/Application Support/UPSSentinel"))
        import ups_guard as ug

        ug.register()
        return ug
    except Exception:  # noqa: BLE001
        return None


# ================= the slaving =================
def solve_m00(S, roots, m_start, iters=12):
    """per cell argmin over M_00 of sum_p ((-M_00)^p + tr S^p - C_p)^2, Newton from m_start."""
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
        g = sum(2.0 * r[k] * d[k] for k in range(4))
        H = sum(2.0 * d[k] ** 2 + 2.0 * r[k] * dd[k] for k in range(4))
        x = x - g / H
    return -x


class Reduced:
    def __init__(self, M, cfg, p, pot, pinned):
        self.cfg, self.p, self.pot = cfg, p, pot
        self.mask = R21.free_mask(cfg, pinned)
        self.M = M.copy()
        self.last = None

    def build(self, x):
        M = self.M.copy()
        S = M[..., 1:, 1:].copy()
        blk = np.zeros((int(self.mask.sum()), 3, 3))
        blk[:, IU3[0], IU3[1]] = x.reshape(-1, 6)
        blk = blk + blk.swapaxes(-1, -2) - np.einsum("...ii->...i", blk)[..., None] * np.eye(3)
        S[self.mask] = blk
        M[..., 1:, 1:] = S
        m00 = M[..., 0, 0].copy()
        m00[self.mask] = solve_m00(S[self.mask], self.pot[1], m00[self.mask])
        M[..., 0, 0] = m00
        return M

    def pack(self, M):
        return M[..., 1:, 1:][self.mask][:, IU3[0], IU3[1]].ravel().copy()

    def fun(self, x):
        M = self.build(x)
        E, Gm, info = R20.energy_grad(M, self.cfg, self.p, self.pot)
        if Gm is None or not np.isfinite(E) or not info["ok"]:
            return 1e30, np.zeros_like(x)
        Gf = Gm[self.mask]
        gs = Gf[:, 1:, 1:][:, IU3[0], IU3[1]] * OFF3
        self.last = (
            float(E),
            float(np.max(np.abs(Gf[:, 1:, 1:]))),
            float(np.max(np.abs(Gf[:, 0, 0]))),
        )
        self.M[..., 0, 0] = M[..., 0, 0]  # warm start of the next Newton solve
        return float(E), gs.ravel()


# ================= reads =================
def chunk_reads(M, cfg, p, pot):
    parts = R20.energy_parts(M, cfg, p, pot)
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    e = R20.density(M, cfg, pot).ravel()
    o = np.argsort(r.ravel())
    cum = np.cumsum(e[o])
    lam = np.linalg.eigvalsh(M[..., 1:, 1:])
    shells = []
    for R in SHELL_R:
        sh = np.abs(r - R) < 0.75 * h
        if sh.any():
            shells.append([R] + [float(lam[..., a][sh].mean()) for a in range(3)])
    deg = R21.degree_surface_reads(M, cfg)["rank2"]
    return {
        "E": parts["E_total"],
        "E_curv": parts["E_curv"],
        "V": parts["V"],
        "virial": parts["E_curv"] / parts["V"] if parts["V"] > 0 else None,
        "r_half": float(r.ravel()[o][np.searchsorted(cum, 0.5 * cum[-1])]),
        "E_inside_r12": float(e[r.ravel() < 12.0].sum()),
        "degree_top": {k: [v["degree"], v["surface_conflicts"]] for k, v in deg.items()},
        "shells_r_small_mid_top": shells,
        "M00_range": [float(M[..., 0, 0].min()), float(M[..., 0, 0].max())],
    }


# ================= jobs =================
def job_tag(j):
    return j["src"] + "_s"


def jobs_all():
    J = []
    for s_ in ("rad", "bia", "perm", "obl"):
        J.append({"src": f"{s_}_pin_d0.3_w25_n32_L48_x9000", "w1s": 25.0, "bnd": "pin"})
    for s_ in ("rad", "bia"):
        J.append({"src": f"{s_}_free_d0.3_w25_n32_L48", "w1s": 25.0, "bnd": "free"})
        J.append({"src": f"{s_}_pin_d0.3_w1_n32_L48", "w1s": 1.0, "bnd": "pin"})
        J.append({"src": f"{s_}_free_d0.3_w1_n32_L48", "w1s": 1.0, "bnd": "free"})
    return J


def run_job(j):
    from scipy.optimize import minimize

    t0 = time.time()
    tag = job_tag(j)
    ug = ups_guard()
    cfg = R21.cfg_of(32, 48.0, G, 0.3)
    p = R21.params_of(G, 0.3)
    pot = ("v4", R0.roots_of(cfg, degenerate=True), W1 * j["w1s"])
    pinned = j["bnd"] == "pin"
    os.makedirs(OUT_NPZ, exist_ok=True)
    stage = os.path.join(OUT_NPZ, tag + "_stage.npz")
    row = dict(j, tag=tag)
    try:
        if os.path.exists(stage):
            Zs = np.load(stage, allow_pickle=True)
            M, done, chunks = Zs["M"], int(Zs["done"]), json.loads(str(Zs["chunks"]))
        else:
            M = np.load(os.path.join(SRC_NPZ, j["src"] + ".npz"))["M"]
            if np.abs(M[..., 0, 1:]).max() > 0:
                raise RuntimeError("the source field is not block-diagonal")
            done, chunks = 0, [
                dict(chunk_reads(M, cfg, p, pot), iters=0, note="the R22-1 end field, unslaved")
            ]
        verdict = "FALLING"
        while done < MAX_ITER:
            red = Reduced(M, cfg, p, pot, pinned)
            st = {"it": 0, "gate_hit": None}

            def cb(xk, red=red, st=st):
                st["it"] += 1
                if st["it"] >= 5 and red.last is not None and red.last[1] < GATE:
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
            E_end, Gm, _ = R20.energy_grad(M, cfg, p, pot)
            Gf = Gm[red.mask]
            rec = dict(
                chunk_reads(M, cfg, p, pot),
                iters=done,
                fmax_spatial=float(np.max(np.abs(Gf[:, 1:, 1:]))),
                fmax_M00=float(np.max(np.abs(Gf[:, 0, 0]))),
                gate_hit_inside_chunk=st["gate_hit"],
                scipy=str(res.message),
            )
            rec["drop"] = chunks[-1]["E"] - rec["E"]
            chunks.append(rec)
            np.savez_compressed(stage, M=M, done=done, chunks=json.dumps(chunks))
            log(
                f"{tag} its {done} E {rec['E']:.6f} drop {rec['drop']:.2e} fmax_sp {rec['fmax_spatial']:.2e} "
                f"virial {rec['virial']:.2f} r_half {rec['r_half']:.1f} deg {rec['degree_top']}"
            )
            if rec["fmax_spatial"] < GATE and 0 <= rec["drop"] < GATE_DROP and st["it"] >= 5:
                verdict = "AT_GATE"
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
    return {"task": "M5.32 R22-1s", "rows": {}}


def run_pool(jobs, workers):
    rows = load_json()["rows"]
    pending = [j for j in jobs if job_tag(j) not in rows]
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


def smoke():
    cfg = R21.cfg_of(32, 48.0, G, 0.3)
    p = R21.params_of(G, 0.3)
    pot = ("v4", R0.roots_of(cfg, degenerate=True), W1 * 25.0)
    M = np.load(os.path.join(SRC_NPZ, "bia_pin_d0.3_w25_n32_L48.npz"))["M"]
    out = {}
    red = Reduced(M, cfg, p, pot, True)
    x0 = red.pack(M)
    E0, G0, _ = R20.energy_grad(M, cfg, p, pot)
    E1, g1 = red.fun(x0)
    out["E_unslaved"], out["E_slaved"] = float(E0), E1
    out["max_G00_unslaved"] = float(np.max(np.abs(G0[red.mask][:, 0, 0])))
    out["max_G00_slaved"] = red.last[2]
    rng = np.random.default_rng(5)
    fd = []
    for _ in range(3):
        v = rng.normal(size=x0.shape)
        v /= np.linalg.norm(v)
        t = 1e-4
        num = (red.fun(x0 + t * v)[0] - red.fun(x0 - t * v)[0]) / (2 * t)
        fd.append(abs(num - float(g1 @ v)) / max(abs(num), 1e-300))
    out["grad_fd_rel"] = fd
    Mb = red.build(x0)
    out["pinned_shell_unchanged"] = bool(np.array_equal(Mb[~red.mask], M[~red.mask]))
    out["PASS"] = bool(
        E1 <= E0 + 1e-12
        and out["max_G00_slaved"] < 1e-6 * max(out["max_G00_unslaved"], 1e-300) + 1e-9
        and max(fd) < 1e-5
        and out["pinned_shell_unchanged"]
    )
    print(json.dumps(out, indent=1))
    with open(OUT_JSON.replace(".json", "_smoke.json"), "w") as f:
        json.dump(out, f, indent=1)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        smoke()
    elif mode == "run":
        run_pool(jobs_all(), int(sys.argv[2]) if len(sys.argv) > 2 else 10)


if __name__ == "__main__":
    main()
