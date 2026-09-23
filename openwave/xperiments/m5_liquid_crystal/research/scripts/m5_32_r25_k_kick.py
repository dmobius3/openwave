"""M5.32 R25-K: the kick-and-continue pass owed from R24 over the R22-1 and R23-1 rows.

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1), N = M eta. The rows live on
the UNIAXIAL exterior (the N-spectrum (-g, 1, delta, delta), g = 8): the R20 energy
    E[M] = 4 h^3 sum_br wt sum_cells I1(A) + V[M],  A_i = d_i M (sym stencil),
    V = w sum_p (tr N^p - C_p)^2 - c L,  C_p = sum_i q_i^p, w = W1 x w1s,
with c = 0 on the R22-1 rows and the row's own c on the R23-1 rows (L the unique
linear invariant of R23-0). Nothing in the energy is new here: every number comes
from the row's OWN instrument, consumed read-only.

THE OBJECTS
-----------
The 17 base rows of data/m5_32_r22_1_cores.json (fields data/m5_32_r22_1/<tag>.npz;
the 10-entry L-BFGS polish of R21 with the pin optional) and the 16 rows of
data/m5_32_r23_1_cscan.json (fields data/m5_32_r23_1/<tag>.npz; the slaved M_00
and the L-BFGS on the six spatial entries of R23-1, the pinned shell). The R24
audit found that the fresh symmetric-seed rows of these families converged inside
a 12-operation symmetric subspace with the core on a lattice vertex, so their end
states may be symmetry-protected saddles; no energy of theirs is a minimum until a
kick that breaks the symmetry has been continued and found to return.

THE PROTOCOL (per row, two continuations from the stored end field, 500 iterations
each with the row's own instrument)
    CONTROL  the stored field continued unkicked (the drift the instrument still
             has at this budget, which the kicked run is compared against)
    KICKED   0.02 x N(0, 1) added to the six upper-triangle entries of the spatial
             3x3 block on the FREE cells only (symmetrized; the generator seeded by
             the CRC-32 of the tag), the field kept block-diagonal, M_00 re-slaved
             per free cell on the R23-1 rows (their instrument raises otherwise),
             then continued
Reads per continuation: E before (recomputed from the stored field), E after, the
max abs gradient on the free cells' spatial block after, the core position before
and after (the energy-density-weighted centroid over the cells whose density
exceeds half its maximum) and the kicked core's shift against the control in
units of h.

PRE-REGISTERED LABELS (per row)
    STABLE        abs(E_kick - E_ctrl) < 1e-5 max(1, abs(E)) after the continuation
                  AND the kicked core moved under 0.2 h relative to the control
    SADDLE        E_kick < E_ctrl - 1e-4 max(1, abs(E)): the kick found a lower
                  state (its energy and the core shift recorded)
    UNRESOLVED    anything else (both energies recorded)
    INSUFFICIENT  the row's status is not OK or its field file is absent (not run)
The R22 and R23 record claims that rest on a SADDLE row are re-labeled in the R25
section of the task record; the continuation clause enters the gate of every
later rung.

Modes: smoke | run [workers] | collect. Output: data/m5_32_r25_k_kick.json (rows
keyed by the source tag, atomic saves, resumable: a row with status OK is skipped),
the two end fields per row in data/m5_32_r25_k/<tag>_ctrl.npz and <tag>_kick.npz
(local, gitignored). The smoke reads no stored field (synthetic n 16 fields only).
Regenerate: run 12 about 2 to 3 hours (66 continuations, the n 48 rows the long
pole); collect seconds.
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
OUT_JSON = os.path.join(DATA, "m5_32_r25_k_kick.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r25_k")
SRC = {
    "r22_1": (os.path.join(DATA, "m5_32_r22_1_cores.json"), os.path.join(DATA, "m5_32_r22_1")),
    "r23_1": (os.path.join(DATA, "m5_32_r23_1_cscan.json"), os.path.join(DATA, "m5_32_r23_1")),
}
G = 8.0
ITERS = 500
KICK_AMP = 0.02
TOL_STABLE, TOL_SADDLE, CORE_TOL_H = 1e-5, 1e-4, 0.2
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


CS = _load("m5_32_r23_1_cscan", "m5_32_r23_1_cscan.py")
R21, R20, R0, B3, W1 = CS.R21, CS.R20, CS.R0, CS.B3, CS.W1
IU3 = CS.IU3


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


# ================= the rows =================
def source_rows():
    """the 33 source rows in a fixed order: (family, tag, source row)."""
    out = []
    for fam in ("r22_1", "r23_1"):
        jpath, _ = SRC[fam]
        if not os.path.exists(jpath):
            continue
        with open(jpath) as f:
            rows = json.load(f)["rows"]
        for tag in sorted(rows):
            out.append((fam, tag, rows[tag]))
    return out


def setup_of(fam, r):
    cfg = R21.cfg_of(r["n"], r["L"], G, r["delta"])
    p = R21.params_of(G, r["delta"])
    pot = ("v4", R0.roots_of(cfg, degenerate=True), W1 * r["w1s"])
    pinned = True if fam == "r23_1" else (r.get("bnd") == "pin")
    c = float(r.get("c", 0.0)) if fam == "r23_1" else 0.0
    return cfg, p, pot, pinned, c


def energy_of(fam, M, cfg, p, pot, c):
    if fam == "r23_1":
        return float(CS.energy_grad(M, cfg, p, pot, c, need_grad=False)[0])
    return float(R20.energy_parts(M, cfg, p, pot)["E_total"])


def fmax_of(fam, M, cfg, p, pot, c, mask):
    if fam == "r23_1":
        _, Gm, _ = CS.energy_grad(M, cfg, p, pot, c)
    else:
        _, Gm, _ = R20.energy_grad(M, cfg, p, pot)
    return float(np.max(np.abs(Gm[mask][:, 1:, 1:])))


def core_of(M, cfg, pot):
    """the energy-density-weighted centroid over the cells above half the peak density."""
    e = R20.density(M, cfg, pot)
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    sel = e > 0.5 * e.max()
    wts = e[sel]
    return [float(np.sum(wts * C[sel]) / np.sum(wts)) for C in (X, Y, Z)]


def kick_field(M, mask, amp, seed, reslave, roots, w, c):
    """the six spatial entries kicked on the free cells, block-diagonal kept, M_00 re-slaved if asked."""
    rng = np.random.default_rng(seed)
    Mk = M.copy()
    S = Mk[..., 1:, 1:].copy()
    nf = int(mask.sum())
    noise = amp * rng.standard_normal((nf, 6))
    blk = np.zeros((nf, 3, 3))
    blk[:, IU3[0], IU3[1]] = noise
    blk = blk + blk.swapaxes(-1, -2) - np.einsum("...ii->...i", blk)[..., None] * np.eye(3)
    S[mask] = S[mask] + blk
    Mk[..., 1:, 1:] = S
    Mk[..., 0, 1:] = 0.0
    Mk[..., 1:, 0] = 0.0
    if reslave:
        m00 = Mk[..., 0, 0].copy()
        m00[mask] = CS.solve_m00(S[mask], roots, w, c, m00[mask])
        Mk[..., 0, 0] = m00
    return Mk


def continue_r23(M, cfg, p, pot, c, iters, tag):
    """CS.run_job's descent, consumed as it is: chunks of CS.CHUNK L-BFGS-B iterations on the
    reduced (slaved M_00) field, until `iters` iterations are done."""
    from scipy.optimize import minimize

    done = 0
    while done < iters:
        red = CS.Reduced(M, cfg, p, pot, c, True)
        st = {"it": 0}

        def cb(xk, st=st):
            st["it"] += 1

        res = minimize(
            red.fun,
            red.pack(M),
            jac=True,
            method="L-BFGS-B",
            callback=cb,
            options={
                "maxcor": 20,
                "maxiter": min(CS.CHUNK, iters - done),
                "maxfun": 3 * min(CS.CHUNK, iters - done),
                "gtol": 1e-14,
                "ftol": 1e-16,
            },
        )
        M = red.build(np.asarray(res.x))
        done += max(1, st["it"])
        log(
            f"{tag} r23 chunk done at {done} its, E {red.last[0] if red.last else float('nan'):.6f}"
        )
        if st["it"] < 5:
            break  # a line-search stall: the instrument's own stop
    return M, done


def continue_row(fam, M, cfg, p, pot, c, pinned, iters, tag):
    if fam == "r23_1":
        return continue_r23(M, cfg, p, pot, c, iters, tag)
    Mp, pol = R21.polish(M, cfg, p, pot, pinned, tag, max_iter=iters, gate=1e-3, log_every=100)
    return Mp, int(pol["iters"])


def classify(E_ref, E_ctrl, E_kick, core_shift_h):
    scale = max(1.0, abs(E_ref))
    if E_kick < E_ctrl - TOL_SADDLE * scale:
        return "SADDLE"
    if abs(E_kick - E_ctrl) < TOL_STABLE * scale and core_shift_h < CORE_TOL_H:
        return "STABLE"
    return "UNRESOLVED"


def run_job(job):
    fam, tag, r = job
    t0 = time.time()
    row = {
        "tag": tag,
        "family": fam,
        "n": r["n"],
        "L": r["L"],
        "delta": r["delta"],
        "w1s": r["w1s"],
        "c": float(r.get("c", 0.0)) if fam == "r23_1" else 0.0,
        "source_label": r.get("label"),
        "kick_amp": KICK_AMP,
        "kick_seed": int(zlib.crc32(tag.encode())),
    }
    fpath = os.path.join(SRC[fam][1], tag + ".npz")
    if r.get("status") != "OK" or not os.path.exists(fpath):
        row.update(
            status="OK",
            label="INSUFFICIENT",
            why=f"source status {r.get('status')}, field present {os.path.exists(fpath)}",
            wall_s=0.0,
        )
        log(f"DONE {tag} INSUFFICIENT")
        return row
    try:
        cfg, p, pot, pinned, c = setup_of(fam, r)
        row.update(h=cfg["h"], pinned=pinned, iters=ITERS)
        mask = R21.free_mask(cfg, pinned)
        M0 = np.load(fpath)["M"].astype(np.float64)
        row["E_before"] = energy_of(fam, M0, cfg, p, pot, c)
        row["core_before"] = core_of(M0, cfg, pot)
        os.makedirs(OUT_NPZ, exist_ok=True)
        out = {}
        for kind in ("control", "kicked"):
            Ms = (
                M0
                if kind == "control"
                else kick_field(
                    M0, mask, KICK_AMP, row["kick_seed"], fam == "r23_1", pot[1], pot[2], c
                )
            )
            E_start = energy_of(fam, Ms, cfg, p, pot, c)
            Me, its = continue_row(fam, Ms, cfg, p, pot, c, pinned, ITERS, f"{tag}:{kind}")
            rec = {
                "E_start": E_start,
                "E_after": energy_of(fam, Me, cfg, p, pot, c),
                "fmax_after": fmax_of(fam, Me, cfg, p, pot, c, mask),
                "core": core_of(Me, cfg, pot),
                "iters_done": its,
            }
            rec["drop"] = E_start - rec["E_after"]
            np.savez_compressed(
                os.path.join(OUT_NPZ, f"{tag}_{'ctrl' if kind == 'control' else 'kick'}.npz"),
                M=Me.astype(np.float64),
            )
            out[kind] = rec
        shift = float(
            np.linalg.norm(np.array(out["kicked"]["core"]) - np.array(out["control"]["core"]))
            / cfg["h"]
        )
        out["kicked"]["core_shift_h"] = shift
        row.update(out)
        row["label"] = classify(
            row["E_before"], out["control"]["E_after"], out["kicked"]["E_after"], shift
        )
        row["status"] = "OK"
    except Exception as e:  # noqa: BLE001
        import traceback

        row.update(status="FAILED", stop=repr(e), traceback=traceback.format_exc())
    row["wall_s"] = round(time.time() - t0, 1)
    log(
        f"DONE {tag} status {row['status']} label {row.get('label')} "
        f"E_ctrl {row.get('control', {}).get('E_after')} E_kick {row.get('kicked', {}).get('E_after')} "
        f"wall {row['wall_s']}"
    )
    return row


# ================= json + pool =================
def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R25-K", "rows": {}}


def save_json(J):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1)
    os.replace(tmp, OUT_JSON)


def run_pool(workers):
    workers = min(int(workers), 12)
    rows = load_json()["rows"]
    pending = [j for j in source_rows() if rows.get(j[1], {}).get("status") != "OK"]
    log(f"pool: {len(pending)} jobs, {workers} workers")
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        futs = [ex.submit(run_job, j) for j in pending]
        for fut in as_completed(futs):
            row = fut.result()
            Jn = load_json()
            Jn["rows"][row["tag"]] = row
            save_json(Jn)
    log("pool done")


# ================= collect =================
def collect():
    J = load_json()
    rows = J["rows"]
    lines = [
        "| tag | family | label | E_before | E_ctrl | E_kick | core shift (h) |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]
    counts = {}
    for fam, tag, _ in source_rows():
        r = rows.get(tag)
        if r is None:
            continue
        lab = r.get("label") if r.get("status") == "OK" else r.get("status")
        counts[lab] = counts.get(lab, 0) + 1
        fmt = lambda v: "" if v is None else f"{v:.6f}"  # noqa: E731
        lines.append(
            f"| `{tag}` | {fam} | {lab} | {fmt(r.get('E_before'))} "
            f"| {fmt(r.get('control', {}).get('E_after'))} | {fmt(r.get('kicked', {}).get('E_after'))} "
            f"| {'' if r.get('kicked') is None else f'{r['kicked'].get('core_shift_h', float('nan')):.2f}'} |"
        )
    table = "\n".join(lines)
    print(table)
    print("counts:", json.dumps(counts))
    J["collect"] = {"counts": counts, "table": table, "rows_done": len(rows)}
    save_json(J)


# ================= smoke =================
def smoke():
    """synthetic only: no stored field is read."""
    global OUT_NPZ
    import shutil
    import tempfile

    out = {"classify": {}, "wiring": {}}
    out["classify"]["stable"] = classify(5.0, 5.000001, 5.000002, 0.05)
    out["classify"]["saddle"] = classify(5.0, 5.0, 4.99, 1.5)
    out["classify"]["unresolved"] = classify(5.0, 5.0, 5.0, 0.8)
    ok = (
        out["classify"]["stable"] == "STABLE"
        and out["classify"]["saddle"] == "SADDLE"
        and out["classify"]["unresolved"] == "UNRESOLVED"
    )
    src = source_rows()
    tags = [t for _, t, _ in src]
    fields = sum(os.path.exists(os.path.join(SRC[f][1], t + ".npz")) for f, t, _ in src)
    out["wiring"] = {
        "source_json_present": {k: os.path.exists(v[0]) for k, v in SRC.items()},
        "rows": len(tags),
        "tags_unique": len(set(tags)) == len(tags),
        "source_fields_present": int(fields),
    }
    ok = ok and out["wiring"]["rows"] == 33 and out["wiring"]["tags_unique"]
    # the two instrument paths on a synthetic n 16 field, two iterations each
    OUT_NPZ = tempfile.mkdtemp(prefix="r25_k_smoke_", dir=DATA)
    chunk0 = CS.CHUNK
    CS.CHUNK = 2
    try:
        cfg = R21.cfg_of(16, 24.0, G, 0.3)
        M = CS.CORES.seed_core(cfg, "rad", 0.3)
        for fam, r in (
            ("r22_1", {"n": 16, "L": 24.0, "delta": 0.3, "w1s": 25.0, "bnd": "pin"}),
            ("r23_1", {"n": 16, "L": 24.0, "delta": 0.3, "w1s": 25.0, "c": 1e-3}),
        ):
            cfg, p, pot, pinned, c = setup_of(fam, r)
            mask = R21.free_mask(cfg, pinned)
            E0 = energy_of(fam, M, cfg, p, pot, c)
            Mk = kick_field(M, mask, KICK_AMP, 7, fam == "r23_1", pot[1], pot[2], c)
            Ek = energy_of(fam, Mk, cfg, p, pot, c)
            Me, its = continue_row(fam, Mk, cfg, p, pot, c, pinned, 2, f"smoke:{fam}")
            Ee = energy_of(fam, Me, cfg, p, pot, c)
            rec = {
                "E0": E0,
                "E_kicked": Ek,
                "E_after_2_its": Ee,
                "iters": its,
                "block_diagonal_kept": bool(np.abs(Mk[..., 0, 1:]).max() == 0.0),
                "pinned_cells_untouched": bool(np.array_equal(Mk[~mask], M[~mask])),
                "core_shift_h": float(
                    np.linalg.norm(
                        np.array(core_of(Me, cfg, pot)) - np.array(core_of(M, cfg, pot))
                    )
                    / cfg["h"]
                ),
                "fmax_after": fmax_of(fam, Me, cfg, p, pot, c, mask),
            }
            rec["kick_raised_E"] = Ek > E0
            rec["descent_lowered_E"] = Ee < Ek
            out[fam] = rec
            ok = ok and rec["block_diagonal_kept"] and rec["pinned_cells_untouched"]
            ok = ok and rec["kick_raised_E"] and rec["descent_lowered_E"]
    finally:
        CS.CHUNK = chunk0
        shutil.rmtree(OUT_NPZ, ignore_errors=True)
    out["PASS"] = bool(ok)
    with open(OUT_JSON.replace(".json", "_smoke.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        smoke()
    elif mode == "run":
        run_pool(int(sys.argv[2]) if len(sys.argv) > 2 else 12)
    elif mode == "collect":
        collect()
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
