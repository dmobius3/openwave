"""M5.32 R24-2: the threshold of the biaxial halo between c 3e-3 and 1e-2.

R23-1 found the split amplitude at r 6 falling 12.7 times (h 1.5) and 26.7
times (h 1.0) between c 3e-3 and 1e-2, where a cutoff moving through r 6 allows
2.5 at most. The 2026-09-20 reply reads this as the cutoff R_c = beta^(-1/4)
crossing the core near c 5e-3. R24-0 (e3) finds no linear term in the split on
the hedgehog exterior and a positive second variation, so the far field cannot
start the halo: it starts in the core, and what kind of transition that is, is
what this ladder reads.

The instrument is R23-1's, unchanged (m5_32_r23_1_cscan.py imported read-only:
the slaved M_00, L-BFGS on the six spatial entries, chunks of 250, the gate
fmax_spatial < 1e-3 read inside a chunk and a last-chunk drop under 1e-4 abs(E)).
Pinned exterior, delta 0.3, W1 x 25, V = V4 - c L.

Jobs (16), all from the analytic `rad` seed unless said: on n 48 L 48 (h 1.0)
c = 4e-3, 5e-3, 6e-3, 7e-3, 8e-3; on n 32 L 48 (h 1.5) the same five and
3.5e-3, 4.5e-3, 5.5e-3, 6.5e-3; and two DOWNWARD rows on n 32, c 5e-3 and 7e-3,
started from the stored R23 c 3e-3 end field (the halo state), which with the
seed rows at the same c read hysteresis.

Reads per row: eps(6) and eps(4.5) (shell mean, R24-1's profile), the floor
eps_floor(h) = eps(6) of the R23 c 1e-2 row of the same spacing, the energy
inside r 6, the degree of the top eigenvector (R21's reader, stored by the
instrument).

Pre-registered labels (the R23 rows at c 3e-3 and 1e-2 join the ladder):
    THRESHOLD_FIRST_ORDER   a downward row and the seed row at the same c end in different
                            states (eps(6) apart by over 3 times, both at the gate), or the
                            last ladder point above 3 floors still holds 0.8 of eps(6) at
                            c 3e-3 (a square-root law sampled at these steps allows 0.71)
    THRESHOLD_CONTINUOUS    no such jump, and eps(6)^2 - eps_floor^2 is linear in c on at
                            least 3 ladder points below its zero (R^2 >= 0.98), the zero
                            c* agreeing between the two spacings within 20 percent
    THRESHOLD_SMOOTH        no jump and no linear law: Q(c) of R24-1 (eps(6) over the master
                            prediction) within a factor 2 across the ladder, the fall
                            explained by the cutoff alone
    UNRESOLVED              none of these, or fewer than 3 gate rows per spacing
The reply's estimate (the cutoff crossing the core) is compared as c* against
4.9e-3 to 6.7e-3; it is a number to compare, not a label.

Modes: smoke | wire | run [workers] | collect | jobs. Order at the run: run, then
m5_32_r24_1_collapse.py run (it reads these rows too), then collect (it reads Q(c) back). Output: data/m5_32_r24_2_threshold.json
(its own file, its own folder data/m5_32_r24_2/ for the local arrays).
Regenerate: run 12 about 3 h (the n 48 rows the long pole); collect seconds.
"""

import importlib.util
import json
import multiprocessing as mp
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r24_2_threshold.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r24_2")
C_LADDER = (4e-3, 5e-3, 6e-3, 7e-3, 8e-3)
C_LADDER_32_EXTRA = (3.5e-3, 4.5e-3, 5.5e-3, 6.5e-3)
SRC_TAG = "rad_pin_d0.3_w25_c0.003_n32_L48"


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


CS = _load("m5_32_r23_1_cscan", "m5_32_r23_1_cscan.py")
COL = _load("m5_32_r24_1_collapse", "m5_32_r24_1_collapse.py")
R23_NPZ = CS.OUT_NPZ


def jobs_all():
    base = {"delta": 0.3, "w1s": 25.0, "seed": "rad", "src": None}
    J = [
        dict(base, part="R24-2 seed h 1.0", c=c, n=48, L=48.0) for c in C_LADDER
    ]  # the long pole first
    J += [
        dict(base, part="R24-2 seed h 1.5", c=c, n=32, L=48.0)
        for c in C_LADDER + C_LADDER_32_EXTRA
    ]
    J += [
        dict(base, part="R24-2 downward h 1.5", c=c, n=32, L=48.0, down=True) for c in (5e-3, 7e-3)
    ]
    return J


def jobs_ext():
    """the extension logged at EXECUTE: the first ladder left no h 1.0 point between c 3e-3 and 4e-3, where the
    whole fall happens, and its two downward rows sat far above the fall. Fresh rows inside the fall at both
    spacings, and started rows from both sides of it (from the R23 halo state, from this rung's residue state).
    """
    base = {"delta": 0.3, "w1s": 25.0, "seed": "rad", "src": None}
    J = [
        dict(base, part="R24-2 ext seed h 1.0", c=c, n=48, L=48.0)
        for c in (3.25e-3, 3.5e-3, 3.75e-3)
    ]
    J += [
        dict(base, part="R24-2 ext seed h 1.5", c=c, n=32, L=48.0)
        for c in (3.25e-3, 3.75e-3, 4.25e-3)
    ]
    halo = lambda n: (
        "r23",
        f"rad_pin_d0.3_w25_c0.003_n{n}_L48",
        "the R23 c 3e-3 end field",
    )  # noqa: E731
    J += [
        dict(
            base,
            part="R24-2 ext from the halo",
            c=4e-3,
            n=n,
            L=48.0,
            down=True,
            start=("rad_down",) + halo(n),
        )
        for n in (48, 32)
    ]
    J += [
        dict(
            base,
            part="R24-2 ext from the residue",
            c=c,
            n=n,
            L=48.0,
            down=True,
            start=(
                "rad_up",
                "r24",
                f"rad_pin_d0.3_w25_c{src}_n{n}_L48",
                f"this rung's c {src} end field",
            ),
        )
        for c, n, src in ((3.5e-3, 48, "0.004"), (3.5e-3, 32, "0.005"), (4e-3, 32, "0.005"))
    ]
    return J


def jobs_ext2():
    """the second extension logged at EXECUTE: the first one found two gate states at c 4e-3 on both spacings (the
    halo-started one lower in energy), so the fresh seed's fall is a basin boundary of the seed and not the end of
    the halo branch. The halo branch followed upward from its c 4e-3 state, to bracket where it ends.
    """
    base = {"delta": 0.3, "w1s": 25.0, "seed": "rad", "src": None}
    J = []
    for n, cs in ((48, (4.5e-3, 5e-3, 6e-3)), (32, (4.25e-3, 4.5e-3, 4.75e-3))):
        src = f"rad_down_pin_d0.3_w25_c0.004_n{n}_L48"
        J += [
            dict(
                base,
                part="R24-2 ext2 along the halo branch",
                c=c,
                n=n,
                L=48.0,
                down=True,
                start=("rad_down", "r24", src, "this rung's halo-started c 4e-3 end field"),
            )
            for c in cs
        ]
    return J


def run_job(j):
    """R23-1's run_job on this rung's folder; a downward row is pre-staged from the stored halo state."""
    CS.OUT_NPZ = OUT_NPZ
    os.makedirs(OUT_NPZ, exist_ok=True)
    jj = dict(j)
    if j.get("down"):
        seed, folder, src_tag, note = j.get("start") or (
            "rad_down",
            "r23",
            SRC_TAG,
            "the R23 c 3e-3 end field",
        )
        jj.pop("start", None)
        jj["seed"] = seed  # only the tag changes: the stage file below carries the start field
        stage = os.path.join(OUT_NPZ, CS.job_tag(jj) + "_stage.npz")
        if not os.path.exists(stage):
            M = np.load(os.path.join(R23_NPZ if folder == "r23" else OUT_NPZ, src_tag + ".npz"))[
                "M"
            ]
            cfg = CS.R21.cfg_of(j["n"], j["L"], CS.G, j["delta"])
            p = CS.R21.params_of(CS.G, j["delta"])
            pot = ("v4", CS.R0.roots_of(cfg, degenerate=True), CS.W1 * j["w1s"])
            start = dict(
                CS.chunk_reads(M, cfg, p, pot, j["c"]),
                iters=0,
                note="the start field: " + note,
            )
            np.savez_compressed(stage, M=M, done=0, chunks=json.dumps([start]))
    row = CS.run_job(jj)
    row["down"] = bool(j.get("down"))
    if j.get("start"):
        row["start_from"] = j["start"][2]
    return row


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R24-2", "rows": {}}


def save_json(J):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1)
    os.replace(tmp, OUT_JSON)


def run_pool(workers, jobs=None):
    rows = load_json()["rows"]
    pending = [j for j in (jobs or jobs_all()) if rows.get(_row_tag(j), {}).get("status") != "OK"]
    CS.log(f"pool: {len(pending)} jobs, {workers} workers")
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        futs = [ex.submit(run_job, j) for j in pending]
        for fut in as_completed(futs):
            row = fut.result()
            Jn = load_json()
            Jn["rows"][row["tag"]] = row
            save_json(Jn)
    CS.log("pool done")


def _row_tag(j):
    if not j.get("down"):
        return CS.job_tag(j)
    jj = {k: v for k, v in j.items() if k != "start"}
    return CS.job_tag(dict(jj, seed=j["start"][0] if j.get("start") else "rad_down"))


# ================= collect =================
def eps_at(tag, folder, n, L, radii=(4.5, 6.0)):
    f = os.path.join(folder, tag + ".npz")
    if not os.path.exists(f):
        return None
    M = np.load(f)["M"]
    cfg = CS.R21.cfg_of(n, L, CS.G, 0.3)
    X, Y, Z = CS.B3.coords(cfg["n"], cfg["h"])
    r = np.sqrt(X * X + Y * Y + Z * Z)
    prof = COL.shell_profile(COL.eps_of_field(M), r, cfg["h"])
    return [float(np.interp(R, prof[:, 0], prof[:, 1])) for R in radii]


def classify(ladder):
    """ladder: {h: [ {c, eps6, at_gate, down} ... ]}; the pre-registered labels of the docstring."""
    out = {"per_spacing": {}}
    jump = False
    cstar = {}
    enough = True
    for h, pts in ladder.items():
        seed = sorted((q for q in pts if not q["down"] and q["at_gate"]), key=lambda q: q["c"])
        floor = next((q["eps6"] for q in seed if abs(q["c"] - 1e-2) < 1e-12), None)
        rec = {"floor_eps6_at_c_1e-2": floor, "gate_rows": len(seed)}
        if len(seed) < 3 or floor is None:
            enough = False
            out["per_spacing"][h] = rec
            continue
        rec["adjacent_ratios"] = [a["eps6"] / b["eps6"] for a, b in zip(seed[:-1], seed[1:])]
        ref = next((q["eps6"] for q in seed if abs(q["c"] - 3e-3) < 1e-12), None)
        above = [q for q in seed if q["eps6"] > 3.0 * floor]
        if ref is not None and above and above[-1] is not seed[-1]:
            rec["last_point_above_3_floors_over_eps6_at_3e-3"] = above[-1]["eps6"] / ref
            if above[-1]["eps6"] / ref >= 0.8:
                jump = True
        for d in (q for q in pts if q["down"] and q["at_gate"]):
            s = next((q for q in seed if abs(q["c"] - d["c"]) < 1e-12), None)
            if s is not None:
                hyst = max(d["eps6"], s["eps6"]) / min(d["eps6"], s["eps6"])
                rec.setdefault("hysteresis_ratio", {})[f"{d['c']:g}"] = hyst
                if hyst > 3.0:
                    jump = True
        y = np.array([q["eps6"] ** 2 - floor**2 for q in seed])
        c = np.array([q["c"] for q in seed])
        k = y > 0.05 * y.max()
        if k.sum() >= 3:
            a, b = np.polyfit(c[k], y[k], 1)
            fit = a * c[k] + b
            var = np.sum((y[k] - y[k].mean()) ** 2)
            r2 = 1 - np.sum((y[k] - fit) ** 2) / var if var > 0 else 0.0
            rec.update(
                linear_points=int(k.sum()), r2=float(r2), c_star=float(-b / a) if a < 0 else None
            )
            if r2 >= 0.98 and a < 0:
                cstar[h] = -b / a
        out["per_spacing"][h] = rec
    if not enough:
        out["label"] = "UNRESOLVED"
    elif jump:
        out["label"] = "THRESHOLD_FIRST_ORDER"
    elif len(cstar) == len(ladder) and max(cstar.values()) / min(cstar.values()) < 1.2:
        out["label"] = "THRESHOLD_CONTINUOUS"
        out["c_star"] = {str(h): float(v) for h, v in cstar.items()}
    else:
        out["label"] = "SMOOTH_OR_UNRESOLVED"  # collect() decides with Q(c) of R24-1
    out["reply_estimate_c_star"] = [4.9e-3, 6.7e-3]
    return out


def collect():
    J = load_json()
    r23 = CS.load_all_rows()
    ladder = {}
    table = []
    src = [(OUT_NPZ, t, q) for t, q in J["rows"].items()]
    for c in (3e-3, 1e-2):
        for n in (32, 48):
            t = CS.job_tag(
                {"seed": "rad", "delta": 0.3, "w1s": 25.0, "c": c, "n": n, "L": 48.0, "src": None}
            )
            if t in r23:
                src.append((R23_NPZ, t, r23[t]))
    for folder, t, q in src:
        if q.get("status") != "OK":
            table.append({"tag": t, "status": q.get("status"), "stop": q.get("stop")})
            continue
        e = eps_at(t, folder, q["n"], q["L"])
        last = q["chunks"][-1]
        gate = last.get("fmax_spatial") is not None and last["fmax_spatial"] < q["gate"]
        h = q["L"] / q["n"]
        rec = {
            "tag": t,
            "c": q["c"],
            "h": h,
            "down": bool(q.get("down")),
            "label": q["label"],
            "at_gate": bool(gate),
            "iters": q["iters"],
            "E": last["E"],
            "fmax_spatial": last.get("fmax_spatial"),
            "eps4.5": e[0] if e else None,
            "eps6": e[1] if e else None,
            "r_half": last["r_half"],
            "degree_top": last["degree_top"],
            "wall_s": q.get("wall_s"),
        }
        table.append(rec)
        if e:
            ladder.setdefault(f"{h:g}", []).append(rec)
    thr = classify(ladder)
    if thr["label"] == "SMOOTH_OR_UNRESOLVED":
        thr["label"] = "UNRESOLVED"
        if os.path.exists(
            COL.OUT_JSON
        ):  # R24-1 run AFTER this pool: its supporting read carries Q(c)
            with open(COL.OUT_JSON) as f:
                sup = json.load(f).get("supporting", [])
            gate_tags = {q["tag"] for q in table if q.get("at_gate") and not q.get("down")}
            for A in sorted({q["A"] for q in sup}):
                Q = [
                    q["Q"]
                    for q in sup
                    if q["A"] == A and q["tag"] in gate_tags and 3e-3 <= q["c"] <= 1e-2
                ]
                if len(Q) >= 3:
                    thr.setdefault("Q_max_over_min", {})[f"A_{A:g}"] = float(max(Q) / min(Q))
            if thr.get("Q_max_over_min") and min(thr["Q_max_over_min"].values()) <= 2.0:
                thr["label"] = "THRESHOLD_SMOOTH"
    J["collect"] = {
        "table": sorted(table, key=lambda q: (q.get("h", 0), q.get("c", 0))),
        "threshold": thr,
    }
    save_json(J)
    print(json.dumps(J["collect"], indent=1))


# ================= smoke: the classifier on synthetic ladders, and the job wiring =================
def smoke():
    cs = (3e-3, 3.5e-3, 4e-3, 4.5e-3, 5e-3, 5.5e-3, 6e-3, 6.5e-3, 7e-3, 8e-3, 1e-2)
    floor = 0.004

    def lad(f, down=None):
        pts = [
            {"c": c, "eps6": float(np.hypot(f(c), floor)), "at_gate": True, "down": False}
            for c in cs
        ]
        for c, v in (down or {}).items():
            pts.append({"c": c, "eps6": v, "at_gate": True, "down": True})
        return pts

    pitch = lambda c: 0.09 * np.sqrt(max(0.0, 1 - c / 6.2e-3) / (1 - 3e-3 / 6.2e-3))  # noqa: E731
    step = lambda c: 0.09 if c < 5.5e-3 else 0.0  # noqa: E731
    out = {
        "pitchfork": classify({"1.5": lad(pitch), "1": lad(pitch)}),
        "jump": classify({"1.5": lad(step), "1": lad(step)}),
        "hysteresis": classify({"1.5": lad(pitch, {5e-3: 0.085, 7e-3: 0.08}), "1": lad(pitch)}),
        "smooth": classify(
            {h: lad(lambda c: 0.09 * np.exp(-(c - 3e-3) / 1.5e-3)) for h in ("1.5", "1")}
        ),
    }
    ok = (
        out["pitchfork"]["label"] == "THRESHOLD_CONTINUOUS"
        and abs(out["pitchfork"]["c_star"]["1.5"] - 6.2e-3) < 2e-4
        and out["jump"]["label"] == "THRESHOLD_FIRST_ORDER"
        and out["hysteresis"]["label"] == "THRESHOLD_FIRST_ORDER"
        and out["smooth"]["label"] == "SMOOTH_OR_UNRESOLVED"
    )
    tags = [_row_tag(j) for j in jobs_all()]
    wiring = {
        "jobs": len(tags),
        "tags_unique": len(set(tags)) == len(tags),
        "source_field_present": os.path.exists(os.path.join(R23_NPZ, SRC_TAG + ".npz")),
        "no_tag_collides_with_r23": not any(
            os.path.exists(os.path.join(R23_NPZ, t + ".npz")) for t in tags
        ),
        "c_below_c_crit_0.0863": max(C_LADDER) < 0.0863,
    }
    out["wiring"] = wiring
    out["PASS"] = bool(ok and all(v for k, v in wiring.items() if k != "jobs"))
    with open(OUT_JSON.replace(".json", "_smoke.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


def wire():
    """the job wiring end to end on a throwaway folder: the downward row's pre-stage and three iterations."""
    import shutil
    import tempfile

    global OUT_NPZ
    OUT_NPZ = tempfile.mkdtemp(prefix="r24_2_wire_", dir=DATA)
    try:
        CS.CHUNK, CS.MAX_ITER = 3, {32: 3, 48: 3}
        row = run_job(dict(jobs_all()[-1]))
        rec = {k: row.get(k) for k in ("tag", "status", "label", "iters", "E", "down", "stop")}
        rec["start_note"] = row["chunks"][0].get("note") if row.get("chunks") else None
        rec["E_start_minus_E_end"] = (
            row["chunks"][0]["E"] - row["chunks"][-1]["E"] if row.get("chunks") else None
        )
        print(json.dumps(rec, indent=1))
    finally:
        shutil.rmtree(OUT_NPZ)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        smoke()
    elif mode == "run":
        run_pool(int(sys.argv[2]) if len(sys.argv) > 2 else 12)
    elif mode == "run_ext":
        run_pool(int(sys.argv[2]) if len(sys.argv) > 2 else 12, jobs_ext())
    elif mode == "run_ext2":
        run_pool(int(sys.argv[2]) if len(sys.argv) > 2 else 12, jobs_ext2())
    elif mode == "wire":
        wire()
    elif mode == "collect":
        collect()
    elif mode == "jobs":
        for j in jobs_all() + jobs_ext() + jobs_ext2():
            print(_row_tag(j), j["part"])


if __name__ == "__main__":
    main()
