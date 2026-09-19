"""M5.32 R22-1: inequivalent one-charge core seeds under ONE uniaxial exterior
(the protocol of the lepton paper, Zenodo 22801032, sections 2.1 and 2.4, on
the full 4x4 field with the R21 polish as the instrument).

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1), N = M eta, code
branch s = -1, g = 8. Vacuum: the EXACTLY uniaxial N-spectrum
(-g, 1, delta, delta) (the working note, Zenodo 22788604, section 4 (iv)).
    E[M] = 4 h^3 sum_br wt sum_cells sum_{i<j} <F_ij, F_ij>_eta + V4[M]
    V4   = w sum_p (tr N^p - C_p)^2,  C_p from (-g, 1, delta, delta)
    w    = W1 x 25 (main) or W1 x 1 (control)
the R20 / R21 instrument (m5_32_r20_1_axes.energy_grad with the 'v4' hook on
the degenerate targets, consumed read-only). By R22-0 (c) V4 is quartic-flat
along the split of the degenerate pair (eps^4 coefficient w (4 + 36 delta^2 +
144 delta^4)), so a biaxial core costs no quadratic potential.

The common exterior: the one-charge hedgehog of the eigenvalue-1 director,
    S_ext = delta I + (1 - delta) r-hat r-hat^T.
The seeds differ only in the interior, blend u(r) = 1 - exp(-(r / r_c)^2),
M_sp = u S_ext + (1 - u) C, r_c = 4 (the difference from the common exterior
is exp(-(r / r_c)^2): 1e-4 at r 12, 1e-14 on the pinned shell):
    rad   C = a I, a = (1 + 2 delta) / 3            the radially melted core
    big   the same with r_c = 8                      a wider melted core
    obl   C = diag((1 + delta) / 2, (1 + delta) / 2, delta)   oblate about z
    bia   C = diag(1, delta + e, delta - e), e = min(delta, 0.3)   biaxial
    perm  C = v S_in + (1 - v) a I, S_in = delta I + (1 - delta) phi-hat phi-hat^T,
          v = 1 - exp(-(rho / 2)^2): the eigenvalue 1 moved from the radial
          to the azimuthal axis inside the core (the distinguished core axis
          permuted, section 2.1 item 2), melted on the polar line
On the cubic lattice the permutations of a FIXED core frame over x, y, z are
lattice symmetries of a spherical exterior, so they are one seed, not three;
'perm' is the inequivalent permutation (radial against tangential).

Boundaries: 'pin' the Dirichlet shell (B3.pin_shell, depth 1.6) at the exact
exterior; 'free' no pin (the R21 B_free: the faces one-sided).
Descent: FIRE (R21.descend, FIRE_STEPS accepted steps) then the R21 L-BFGS
polish in chunks of CHUNK iterations, a resumable field saved after each
chunk, gate max |G| < 1e-3 on the free cells; the attained residual is
recorded per row.

Reads per endpoint: the R20 / R21 reads (energy and parts, Derrick, shells,
biaxiality), plus here: the sorted spatial eigenvalues and the biaxiality
(lambda_mid - lambda_small) on the shells r 1.5 to 21, the degree of the
top eigenvector oriented along r-hat through the r 9 cube (with the count of
cells where the orientation is unreliable, |e . r-hat| < 0.2), and at collect
the RMS field difference between endpoints, minimized over the 48 cubic
images M(x) -> R M(R^T x) R^T and, separately, over the 12 of them that are
exact symmetries of the stack's functional (see rms_mod_cubic).

Pre-registered outcomes (the plan post): ONE_BRANCH, BRANCH_BOUNDARY_HELD,
BRANCH_LATTICE_ONLY, BRANCHES_SURVIVE.

Modes: smoke | run [workers] | collect
"""

import importlib.util
import itertools
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r22_1_cores.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r22_1")
G = 8.0
R_C = 4.0
FIRE_STEPS = 600
CHUNK = 500
POLISH_GATE = 1e-3
MEM_GB = {16: 0.2, 32: 0.8, 48: 2.2, 64: 4.0}
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


# ================= seeds =================
def seed_core(cfg, kind, delta):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rh = np.stack([X, Y, Z], axis=-1) / np.maximum(r, 1e-300)[..., None]
    I3 = np.eye(3)
    S_ext = delta * I3 + (1.0 - delta) * rh[..., :, None] * rh[..., None, :]
    a = (1.0 + 2.0 * delta) / 3.0
    r_c = 2.0 * R_C if kind == "big" else R_C
    u = (1.0 - np.exp(-((r / r_c) ** 2)))[..., None, None]
    if kind in ("rad", "big"):
        C = a * I3
    elif kind == "obl":
        C = np.diag([(1.0 + delta) / 2.0, (1.0 + delta) / 2.0, delta])
    elif kind == "bia":
        e = min(delta, 0.3)
        C = np.diag([1.0, delta + e, delta - e])
    elif kind == "perm":
        rho = np.sqrt(X * X + Y * Y)
        ph = np.stack([-Y, X, np.zeros_like(Z)], axis=-1) / np.maximum(rho, 1e-300)[..., None]
        S_in = delta * I3 + (1.0 - delta) * ph[..., :, None] * ph[..., None, :]
        v = (1.0 - np.exp(-((rho / 2.0) ** 2)))[..., None, None]
        C = v * S_in + (1.0 - v) * a * I3
    else:
        raise ValueError(kind)
    return B3.embed34(u * S_ext + (1.0 - u) * C, cfg)


# ================= reads =================
def shell_reads(M, cfg):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    lam, vec = np.linalg.eigh(M[..., 1:, 1:])
    rows = []
    for R in SHELL_R:
        sh = np.abs(r - R) < 0.75 * h
        if not sh.any():
            continue
        rows.append(
            {
                "r": R,
                "cells": int(sh.sum()),
                "lam_small": float(lam[..., 0][sh].mean()),
                "lam_mid": float(lam[..., 1][sh].mean()),
                "lam_top": float(lam[..., 2][sh].mean()),
                "biax_mean": float((lam[..., 1] - lam[..., 0])[sh].mean()),
                "biax_max": float((lam[..., 1] - lam[..., 0])[sh].max()),
                "M0i_max": float(np.abs(M[..., 0, 1:][sh]).max()),
                "M00_dev": float(np.abs(M[..., 0, 0][sh] - M[0, 0, 0, 0, 0]).max()),
            }
        )
    # the degree of the top eigenvector, oriented along r-hat, through the r 9 cube
    e1 = vec[..., :, 2]
    rh = np.stack([X, Y, Z], axis=-1) / np.maximum(r, 1e-300)[..., None]
    dot = np.einsum("...a,...a->...", e1, rh)
    e1 = e1 * np.where(dot < 0, -1.0, 1.0)[..., None]
    dn = [np.gradient(e1, h, axis=ax) for ax in range(3)]
    cr = lambda a_, b_: np.einsum("...a,...a->...", e1, np.cross(dn[a_], dn[b_]))  # noqa: E731
    E = np.stack([cr(1, 2), cr(2, 0), cr(0, 1)], axis=-1)
    c = (np.arange(n) - (n - 1) / 2.0) * h
    idx = np.where(np.abs(c) < 9.0)[0]
    tot, unreliable = 0.0, 0
    for ax in range(3):
        o = [b for b in range(3) if b != ax]
        for sgn, i_face in ((+1, idx[-1]), (-1, idx[0])):
            sl = [None] * 3
            sl[ax] = i_face
            sl[o[0]], sl[o[1]] = idx[:, None], idx[None, :]
            sl2 = list(sl)
            sl2[ax] = i_face + sgn
            tot += (
                sgn * float(np.sum(0.5 * (E[tuple(sl) + (ax,)] + E[tuple(sl2) + (ax,)]))) * h * h
            )
            unreliable += int(np.sum(np.abs(dot[tuple(sl)]) < 0.2))
    return {
        "shells": rows,
        "degree_top_r9": tot / (4 * np.pi),
        "unreliable_orientation_cells_r9": unreliable,
        "biax_max_anywhere": float((lam[..., 1] - lam[..., 0]).max()),
        "lam_top_min": float(lam[..., 2].min()),
    }


def cubic_group():
    out = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1.0, -1.0), repeat=3):
            R = np.zeros((3, 3))
            for i, (pi, s) in enumerate(zip(perm, signs)):
                R[i, pi] = s
            out.append((perm, signs, R))
    return out


def transform(M3, perm, signs, R):
    """(R . M3)(x) = R M3(R^T x) R^T on the cell-centered cube."""
    # spatial part: new[x] = old[R^T x]; R maps old axis perm[i] (times sign) to new axis i
    F = np.transpose(M3, axes=[perm[0], perm[1], perm[2], 3, 4])
    for i, s in enumerate(signs):
        if s < 0:
            F = np.flip(F, axis=i)
    return np.einsum("ab,...bc,dc->...ad", R, F, R)


def rms_mod_cubic(Ma, Mb, mask):
    """min RMS over the 48 cubic images of Mb, and over the 12 that are EXACT symmetries of the stack's
    functional (axis permutations times the full inversion: the fwd / bwd average takes all three axes
    forward or all three backward, so a reflection of one or two axes maps it to a mixed branch that the
    average does not contain; measured in the smoke, 36 of the 48 change E_curv at the 1e-3 level).
    """
    best, best12 = None, None
    A = Ma[..., 1:, 1:]
    for perm, signs, R in cubic_group():
        T = transform(Mb[..., 1:, 1:], perm, signs, R)
        v = float(np.sqrt(np.mean(np.sum((A - T) ** 2, axis=(-1, -2))[mask])))
        if best is None or v < best[0]:
            best = (v, perm, signs)
        if len(set(signs)) == 1 and (best12 is None or v < best12[0]):
            best12 = (v, perm, signs)
    return {
        "rms": best[0],
        "perm": list(best[1]),
        "signs": list(best[2]),
        "rms_exact12": best12[0],
        "perm_exact12": list(best12[1]),
        "signs_exact12": list(best12[2]),
    }


# ================= the job =================
def job_tag(j):
    base = f"{j['seed']}_{j['bnd']}_d{j['delta']:g}_w{j['w1s']:g}_n{j['n']}_L{j['L']:g}"
    return base + (f"_x{j['polish']}" if j.get("ext") else "")


def jobs_ext():
    """The extension rows (added at EXECUTE, 2026-09-19 11:20 EDT, on the interim read: after 3000
    iterations every W1 x 25 endpoint was still falling by 0.08 to 0.24 per 500 iterations, the four
    energies ordered by how biaxial the seed's core was, the profiles converging on one another).
    Each row continues the base row's polished field to 9000 iterations in all (the stage file of
    the base row is copied, the base row is left as it is)."""
    base = {
        "delta": 0.3,
        "n": 32,
        "L": 48.0,
        "fire": FIRE_STEPS,
        "polish": 9000,
        "w1s": 25.0,
        "ext": True,
    }
    return [dict(base, seed=s_, bnd="pin") for s_ in ("rad", "bia", "perm", "obl")]


def run_job(j):
    t0 = time.time()
    tag = job_tag(j)
    ug = ups_guard()
    cfg = R21.cfg_of(j["n"], j["L"], G, j["delta"])
    p = R21.params_of(G, j["delta"])
    pot = ("v4", R0.roots_of(cfg, degenerate=True), W1 * j["w1s"])
    pinned = j["bnd"] == "pin"
    row = dict(j)
    row.update(
        {"tag": tag, "h": cfg["h"], "roots": list(pot[1]), "W1_eff": pot[2], "pinned": pinned}
    )
    os.makedirs(OUT_NPZ, exist_ok=True)
    stage = os.path.join(OUT_NPZ, tag + "_stage.npz")
    if j.get("ext") and not os.path.exists(stage):
        import shutil

        shutil.copyfile(os.path.join(OUT_NPZ, job_tag(dict(j, ext=False)) + "_stage.npz"), stage)
    try:
        M0 = seed_core(cfg, j["seed"], j["delta"])
        row["seed_E"] = R20.energy_parts(M0, cfg, p, pot)
        row["seed_shells"] = shell_reads(M0, cfg)
        done_it, chunks = 0, []
        if os.path.exists(stage):
            Zs = np.load(stage, allow_pickle=True)
            M, done_it = Zs["M"], int(Zs["done_it"])
            chunks = json.loads(str(Zs["chunks"]))
            row["resumed_at_polish_iter"] = done_it
        else:
            M, des = R21.descend(
                M0,
                cfg,
                p,
                pot,
                j["fire"],
                2 * j["fire"],
                tag,
                pinned,
                ckpt_path=os.path.join(OUT_NPZ, tag + "_fireckpt"),
            )
            row["descent"] = {k: v for k, v in des.items() if k != "trace"}
            if not (np.all(np.isfinite(M)) and np.isfinite(des["E_end"])):
                raise RuntimeError(f"FIRE ended non-finite: {des.get('stop')}")
            row["E_fire"] = float(des["E_end"])
            np.savez_compressed(stage, M=M, done_it=0, chunks=json.dumps([]))
        verdict = None
        while done_it < j["polish"]:
            k = min(CHUNK, j["polish"] - done_it)
            M, pol = R21.polish(
                M, cfg, p, pot, pinned, tag, max_iter=k, gate=POLISH_GATE, log_every=100
            )
            done_it += max(1, pol["iters"])
            chunks.append(
                {
                    "iters": pol["iters"],
                    "E_end": pol["E_end"],
                    "fmax_end": pol["fmax_end"],
                    "min_gap_end": pol.get("min_gap_end"),
                    "max_abs_M0i_end": pol.get("max_abs_M0i_end"),
                }
            )
            np.savez_compressed(stage, M=M, done_it=done_it, chunks=json.dumps(chunks))
            verdict = pol["verdict"]
            if verdict == "POLISHED":
                break
            if ug is not None and ug.wrap_up():
                verdict = "UPS wrap-up (resumable)"
                break
        row["polish_chunks"] = chunks
        row["polish_iters"] = done_it
        row["label"] = verdict
        row["fmax_end"] = chunks[-1]["fmax_end"] if chunks else None
        row["end_E"] = R20.energy_parts(M, cfg, p, pot)
        row["E"] = row["end_E"]["E_total"]
        row["end_shells"] = shell_reads(M, cfg)
        try:
            row["derrick"] = R20.derrick_reads(M, cfg, p, pot)
        except Exception as e:  # noqa: BLE001
            row["derrick"] = {"error": repr(e)}
        np.savez_compressed(os.path.join(OUT_NPZ, tag + ".npz"), M=M.astype(np.float64))
        row["status"] = "OK"
    except Exception as e:  # noqa: BLE001
        import traceback

        row.update(status="FAILED", stop=repr(e), traceback=traceback.format_exc(), E=None)
    row["wall_s"] = round(time.time() - t0, 1)
    log(
        f"DONE {tag} status {row['status']} label {row.get('label')} E {row.get('E')} wall {row['wall_s']}"
    )
    return row


def jobs_all():
    J = []
    base = {"delta": 0.3, "n": 32, "L": 48.0, "fire": FIRE_STEPS, "polish": 3000}
    for seed in ("rad", "bia", "obl", "perm", "big"):
        for bnd in ("pin", "free"):
            J.append(dict(base, seed=seed, bnd=bnd, w1s=25.0))
    for seed in ("rad", "bia"):
        for bnd in ("pin", "free"):
            J.append(dict(base, seed=seed, bnd=bnd, w1s=1.0))
    for seed in ("rad", "bia"):  # refinement at the fixed box
        J.append(dict(base, seed=seed, bnd="pin", w1s=25.0, n=48, polish=2000))
    J.append(dict(base, seed="rad", bnd="pin", w1s=25.0, delta=0.89))  # side row (cut line 2)
    J.append(
        dict(base, seed="rad", bnd="pin", w1s=25.0, n=64, L=96.0, polish=1000)
    )  # box row (cut line 3)
    return J


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R22-1", "rows": {}}


def save_json(J):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1)
    os.replace(tmp, OUT_JSON)


def run_pool(jobs, workers, mem_budget=16.0):
    rows = load_json()["rows"]
    pending = [j for j in jobs if job_tag(j) not in rows]
    log(f"pool: {len(pending)} jobs, {workers} workers")
    running = {}
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        while pending or running:
            used = sum(MEM_GB[j["n"]] for j in running.values())
            for j in list(pending):
                if len(running) >= workers:
                    break
                if used + MEM_GB[j["n"]] <= mem_budget:
                    running[ex.submit(run_job, j)] = j
                    pending.remove(j)
                    used += MEM_GB[j["n"]]
                    log(f"submit {job_tag(j)}")
            done = next(as_completed(list(running)))
            row = done.result()
            running.pop(done)
            Jn = load_json()
            Jn["rows"][row["tag"]] = row
            save_json(Jn)
    log("pool done")


def smoke():
    global OUT_NPZ
    OUT_NPZ = os.path.join(DATA, "m5_32_r22_1", "smoke")
    out = {}
    cfg = R21.cfg_of(32, 48.0, G, 0.3)
    p = R21.params_of(G, 0.3)
    pot = ("v4", R0.roots_of(cfg, degenerate=True), W1 * 25.0)
    for seed in ("rad", "big", "obl", "bia", "perm"):
        M0 = seed_core(cfg, seed, 0.3)
        E = R20.energy_parts(M0, cfg, p, pot)
        sh = shell_reads(M0, cfg)
        ext = float(np.abs(M0 - seed_core(cfg, "rad", 0.3))[B3.pin_shell(32, cfg["h"])].max())
        out[seed] = {
            "E_total": E["E_total"],
            "E_curv": E["E_curv"],
            "V": E["V"],
            "degree_top_r9": sh["degree_top_r9"],
            "biax_max": sh["biax_max_anywhere"],
            "shell_diff_from_rad_on_pin": ext,
        }
        log(f"smoke seed {seed}: {out[seed]}")
    # cubic symmetry reader: a rotated copy must come back at rms 0
    Ma = seed_core(cfg, "bia", 0.3)
    perm, signs, R = cubic_group()[17]
    Mb = Ma.copy()
    Mb[..., 1:, 1:] = transform(Ma[..., 1:, 1:], perm, signs, R)
    mask = np.ones(Ma.shape[:3], dtype=bool)
    out["rms_rotated_copy"] = rms_mod_cubic(Ma, Mb, mask)
    out["rms_direct_rotated_copy"] = float(np.sqrt(np.mean(np.sum((Ma - Mb) ** 2, axis=(-1, -2)))))
    out["rms_rad_vs_bia"] = rms_mod_cubic(seed_core(cfg, "rad", 0.3), Ma, mask)
    # a short end-to-end row
    row = run_job(
        {
            "seed": "bia",
            "bnd": "free",
            "delta": 0.3,
            "w1s": 25.0,
            "n": 16,
            "L": 24.0,
            "fire": 20,
            "polish": 20,
        }
    )
    out["short_row"] = {
        k: row.get(k) for k in ("status", "label", "E", "E_fire", "fmax_end", "wall_s", "stop")
    }
    print(json.dumps(out, indent=1, default=str))


def collect():
    """Pre-registered reading (written before any endpoint existed, 2026-09-19 09:20 EDT):
    two endpoints of one (boundary, w, n, L, delta) cell are THE SAME BRANCH iff
        abs(dE) < max(1 percent of E, 3 x the larger last-chunk energy drop)   and
        RMS over r < 12 of the field difference, minimized over the exact 12 symmetries,
        is under 0.25 x the same RMS between their SEEDS.
    ONE_BRANCH: every seed pair is the same branch on both boundaries.
    BRANCH_BOUNDARY_HELD: distinct under the pin, the same branch under the free boundary.
    BRANCH_LATTICE_ONLY: distinct at n 32, the same branch at n 48.
    BRANCHES_SURVIVE: distinct on both boundaries and at n 48."""
    J = load_json()
    ext = OUT_JSON.replace(".json", "_ext.json")
    if os.path.exists(ext):
        with open(ext) as f:
            for t, r in json.load(f)["rows"].items():
                J["rows"][t] = dict(r, seed=r["seed"] + "_x")
    rows = {t: r for t, r in J["rows"].items() if r.get("status") == "OK"}
    cells = {}
    for t, r in rows.items():
        cells.setdefault((r["bnd"], r["w1s"], r["n"], r["L"], r["delta"]), []).append(r)
    out = {"cells": {}}
    for key, rs in sorted(cells.items(), key=lambda kv: str(kv[0])):
        cfg = R21.cfg_of(key[2], key[3], G, key[4])
        X, Y, Z = B3.coords(cfg["n"], cfg["h"])
        mask = np.sqrt(X * X + Y * Y + Z * Z) < 12.0
        name = f"{key[0]}_w{key[1]:g}_n{key[2]}_L{key[3]:g}_d{key[4]:g}"
        cell = {"endpoints": {}, "pairs": []}
        fields = {}
        for r in rs:
            fields[r["seed"]] = (
                np.load(os.path.join(OUT_NPZ, r["tag"] + ".npz"))["M"],
                seed_core(cfg, r["seed"].replace("_x", ""), key[4]),
            )
            ch = r["polish_chunks"]
            drop = abs(ch[-2]["E_end"] - ch[-1]["E_end"]) if len(ch) >= 2 else None
            sh = {s_["r"]: s_ for s_ in r["end_shells"]["shells"]}
            cell["endpoints"][r["seed"]] = {
                "E": r["E"],
                "E_curv": r["end_E"]["E_curv"],
                "V": r["end_E"]["V"],
                "E_seed": r["seed_E"]["E_total"],
                "fmax_end": r["fmax_end"],
                "label": r["label"],
                "polish_iters": r["polish_iters"],
                "last_chunk_drop": drop,
                "degree_top_r9": r["end_shells"]["degree_top_r9"],
                "biax_max": r["end_shells"]["biax_max_anywhere"],
                "lam_top_min": r["end_shells"]["lam_top_min"],
                "biax_r1.5_r3_r6": [sh.get(q, {}).get("biax_mean") for q in (1.5, 3.0, 6.0)],
                "virial_E_curv_over_V": (r.get("derrick") or {}).get("virial_E_curv_over_V"),
                "dE_dlambda": ((r.get("derrick") or {}).get("order3") or {}).get(
                    "dE_dlambda_at_1"
                ),
                "r_half_energy": (r.get("derrick") or {}).get("r_half_energy"),
                "chunk_E": [c["E_end"] for c in ch],
                "wall_s": r["wall_s"],
            }
        seeds = sorted(fields)
        for a, b in itertools.combinations(seeds, 2):
            ea, eb = cell["endpoints"][a], cell["endpoints"][b]
            end = rms_mod_cubic(fields[a][0], fields[b][0], mask)
            sd = rms_mod_cubic(fields[a][1], fields[b][1], mask)
            drops = [x for x in (ea["last_chunk_drop"], eb["last_chunk_drop"]) if x is not None]
            tolE = max(0.01 * abs(ea["E"]), 3 * max(drops) if drops else 0.0)
            same = bool(
                abs(ea["E"] - eb["E"]) < tolE and end["rms_exact12"] < 0.25 * sd["rms_exact12"]
            )
            cell["pairs"].append(
                {
                    "a": a,
                    "b": b,
                    "dE": ea["E"] - eb["E"],
                    "tol_E": tolE,
                    "rms_end_exact12": end["rms_exact12"],
                    "rms_end_48": end["rms"],
                    "rms_seed_exact12": sd["rms_exact12"],
                    "same_branch": same,
                }
            )
        cell["all_same"] = bool(cell["pairs"]) and all(p_["same_branch"] for p_ in cell["pairs"])
        out["cells"][name] = cell
        print(name)
        for s_, e in cell["endpoints"].items():
            print(
                f"  {s_:5s} E {e['E']:.5f} (curv {e['E_curv']:.4f} V {e['V']:.4f}; seed {e['E_seed']:.3f}) fmax {e['fmax_end']:.2e} "
                f"its {e['polish_iters']} drop {e['last_chunk_drop']} deg {e['degree_top_r9']:.3f} biax max {e['biax_max']:.3f} "
                f"virial {e['virial_E_curv_over_V']} dE/dlam {e['dE_dlambda']} r_half {e['r_half_energy']}"
            )
        for p_ in cell["pairs"]:
            print(
                f"  {p_['a']:5s} vs {p_['b']:5s} dE {p_['dE']:+.5f} (tol {p_['tol_E']:.5f}) rms end {p_['rms_end_exact12']:.4f} "
                f"(48: {p_['rms_end_48']:.4f}) seed {p_['rms_seed_exact12']:.4f} same {p_['same_branch']}"
            )
    with open(OUT_JSON.replace(".json", "_collect.json"), "w") as f:
        json.dump(out, f, indent=1, default=str)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        smoke()
    elif mode == "run":
        run_pool(jobs_all(), int(sys.argv[2]) if len(sys.argv) > 2 else 8)
    elif mode == "ext":
        global OUT_JSON
        OUT_JSON = OUT_JSON.replace(".json", "_ext.json")
        run_pool(jobs_ext(), int(sys.argv[2]) if len(sys.argv) > 2 else 4)
    elif mode == "collect":
        collect()


if __name__ == "__main__":
    main()
