"""M5.32 R27-2: the index-partition census under the smoothness term kappa |dM|^2 (the
author's sum k^2 ordering prediction, 2026-09-25 19:04 UTC) and the free-shell twins.

EQUATIONS FIRST
---------------
The stack and the descent are R26-4's (R25-2's reduced L-BFGS-B at c = 0 with M_00 slaved per
cell to V4, the pinned shell of depth 1.6, the gate fmax_spatial < 1e-4 with the last-chunk
drop under 1e-5 abs(E), the post-gate kick, the STABLE / SADDLE / UNRESOLVED / FALLING
labels), with the energy wrapped:
    E = 4 h^3 sum_{i<j} <F_ij, F_ij>_eta + V4 + kappa E_kappa,
    E_kappa = h^3 sum_br w_br sum_i <d_i S, d_i S>   (R27-0, the spatial block, the bond sum),
the gradient of the term added to the reduced descent's (R27-0 check b); the M_00 slot stays
slaved to V4 alone (the term does not act on it). Delta 0.3, w = W1 x 25, n 32, L 48 (h 1.5),
the five seeds of R26-0.seed_partition ({4}, {3,1}, {2,2}, {2,1,1}, {1,1,1,1}) gated by the
partition reader before relaxation.
THE TWO KAPPAS (R27-0's stored table, the pre-registered rule): the seeds' E_kappa spread
across the family (95.5, the sum k^2 part; the hedgehog's common part 962) times kappa equals
0.5 and 2 times the kappa-0 quartic spread of the R26-4 rows (0.439): kappa 2.299e-3 and
9.197e-3. The author's prediction: the smooth winding costs 16 pi kappa k^2 b0^2 ln(L / xi)
per unit length on the slab, so the partitions order by sum k^2 = 1.0 ({1,1,1,1}) < 1.5
({2,1,1}) < 2.0 ({2,2}) < 2.5 ({3,1}) < 4.0 ({4}); the kappa-0 order of R26-4 was {2,1,1} <
{3,1} < {1,1,1,1} < {2,2} < {4} (two rows capped). The seeds sit in the sum k^2 order by
construction (R27-0), the kappa-0 ends do not: the test is the RELAXED rows under kappa.
THE FREE-SHELL TWINS: {3,1} (the certified kappa-0 winner) and {1,1,1,1} at kappa 0 with the
shell free (R21's unpinned mask: every cell relaxes, the faces one-sided), the same descent
and readers, against their pinned twins of R26-4 (energy, the interior partition, the
carriers' positions).
READS per row: E_quartic, E_kappa (kappa = 1 units), the total, the last-chunk drop as the
error, the interior partition on the spheres r 3 to 21, the carriers' distance from the axis,
beta^2, the gap tail, the physical generator, the virial (R26-4's own_reads).

PRE-REGISTERED LABELS (the record's R27 PLANNING section)
--------------------------------------------------------
    ORDER_K2_CONFIRMED    at both kappa at least four of five rows certified (at the gate,
                          STABLE under the kick or SADDLE reconverged) and every certified pair
                          ordered by sum k^2 beyond two errors
    ORDER_K2_REFUTED      at both kappa a certified pair inverted beyond two errors
    ORDER_K2_UNRESOLVED   otherwise
    CENSUS_INSUFFICIENT   at a kappa with fewer than four rows at the gate (R26's rule)
    FREE_SHELL_SAME_CLASS the twin's interior partition on r 9 equals its pinned twin's
    FREE_SHELL_REORGANIZES otherwise; the energies and the carriers' positions reported

Modes: smoke | run [workers] | collect | jobs. Output: data/m5_32_r27_2_kappa_census.json,
arrays in data/m5_32_r27_2/ (local, kept). Regenerate: run 12 about 3 h; collect seconds.
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
OUT_JSON = os.path.join(DATA, "m5_32_r27_2_kappa_census.json")
OUT_NPZ = os.path.join(DATA, "m5_32_r27_2")
R27_0_JSON = os.path.join(DATA, "m5_32_r27_0_form.json")
R26_4_JSON = os.path.join(DATA, "m5_32_r26_4_census.json")
T0 = time.time()
DELTA = 0.3
W1S = 25.0
CAP = {16: 4, 32: 8000}
ORDER = ("1_1_1_1", "2_1_1", "2_2", "3_1", "4")
SUMK2 = {"1_1_1_1": 1.0, "2_1_1": 1.5, "2_2": 2.0, "3_1": 2.5, "4": 4.0}
TWINS = ("3_1", "1_1_1_1")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


F = _load("m5_32_r27_0_form", "m5_32_r27_0_form.py")
R4, F0, R25, R21, R20, B3, R0, W1 = F.R4, F.F0, F.R25, F.R21, F.R20, F.B3, F.R0, F.W1
CS = R25.CS

# ---- the energy wrapped: kappa E_kappa added to the reduced descent's energy and gradient ----
_ORIG_ENERGY_GRAD = CS.energy_grad
CS._R27_KAPPA = 0.0


def energy_grad_kappa(M, cfg, p, pot, c, need_grad=True):
    E, Gm, info = _ORIG_ENERGY_GRAD(M, cfg, p, pot, c, need_grad)
    k = CS._R27_KAPPA
    if k == 0.0 or not np.isfinite(E):
        return E, Gm, info
    ek = F.e_kappa(M, cfg)
    if need_grad and Gm is not None:
        Gm = Gm + k * F.grad_kappa(M, cfg)
    return E + k * ek, Gm, info


CS.energy_grad = energy_grad_kappa
R25.OUT_NPZ = OUT_NPZ


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def kappas():
    J = json.load(open(R27_0_JSON))
    return [float(x) for x in J["stored"]["census"]["kappa_census"]]


def jobs_all(n=32, L=48.0):
    ka, kb = kappas()
    jobs = []
    for k in (ka, kb):
        for part in ORDER:
            jobs.append(
                dict(part=part, delta=DELTA, w1s=W1S, n=n, L=L, cap=CAP[n], kappa=k, pinned=True)
            )
    for part in TWINS:
        jobs.append(
            dict(part=part, delta=DELTA, w1s=W1S, n=n, L=L, cap=CAP[n], kappa=0.0, pinned=False)
        )
    return jobs


def job_tag(j):
    return (
        f"P{j['part']}_d{j['delta']:g}_w{j['w1s']:g}_n{j['n']}_L{j['L']:g}_k{j['kappa']:.3e}"
        + ("" if j["pinned"] else "_free")
    )


def descend_chunks(M, cfg, p, pot, mask, tag, cap, chunks, done, stage, pinned, n_chunks=None):
    """R25-2's chunked descent with the pin as an argument (the free-shell twins)."""
    from scipy.optimize import minimize

    verdict = "FALLING"
    k_run = 0
    while done < cap:
        red = CS.Reduced(M, cfg, p, pot, 0.0, pinned)
        st = {"it": 0, "gate_hit": None}

        def cb(xk, red=red, st=st):
            st["it"] += 1
            if st["it"] >= 5 and red.last is not None and red.last[1] < R25.GATE:
                st["gate_hit"] = st["it"]

        res = minimize(
            red.fun,
            red.pack(M),
            jac=True,
            method="L-BFGS-B",
            callback=cb,
            options={
                "maxcor": 20,
                "maxiter": CS.CHUNK,
                "maxfun": 3 * CS.CHUNK,
                "gtol": 1e-14,
                "ftol": 1e-16,
            },
        )
        M = red.build(np.asarray(res.x))
        done += max(1, st["it"])
        k_run += 1
        rec = R25.chunk_light(M, cfg, p, pot, mask, done)
        rec["gate_hit_inside_chunk"] = st["gate_hit"]
        rec["chunk_iters"] = st["it"]
        rec["scipy"] = str(res.message)
        rec["drop"] = chunks[-1]["E"] - rec["E"]
        chunks.append(rec)
        np.savez_compressed(stage, M=M, done=done, chunks=json.dumps(chunks))
        log(
            f"{tag} its {done} E {rec['E']:.7g} drop {rec['drop']:.2e} fmax_sp {rec['fmax_spatial']:.2e} (gate {R25.GATE:.0e})"
        )
        if (
            rec["fmax_spatial"] < R25.GATE
            and 0 <= rec["drop"] < R25.GATE_DROP * abs(rec["E"])
            and st["it"] >= 5
        ):
            verdict = "AT_GATE"
            break
        if st["it"] < 5 and rec["drop"] <= 0:
            verdict = "LINE_SEARCH_STALL"
            break
        if n_chunks is not None and k_run >= n_chunks:
            verdict = "CHUNKS_DONE"
            break
    return M, done, chunks, verdict


def energy_split(M, cfg, p, pot, kappa):
    parts = R20.energy_parts(M, cfg, p, pot)
    ek = F.e_kappa(M, cfg)
    return {
        "E_quartic": float(parts["E_total"]),
        "E_curv": float(parts["E_curv"]),
        "V4": float(parts["V"]),
        "E_kappa": ek,
        "kappaE_kappa": kappa * ek,
        "E_total": float(parts["E_total"]) + kappa * ek,
    }


def run_job(j):
    t0 = time.time()
    tag = job_tag(j)
    CS._R27_KAPPA = float(j["kappa"])
    ug = CS.SL.ups_guard()
    cfg, p, pot = R4.cfg_pot(j)
    mask = R21.free_mask(cfg, j["pinned"])
    os.makedirs(OUT_NPZ, exist_ok=True)
    stage = os.path.join(OUT_NPZ, tag + "_stage.npz")
    row = dict(
        j,
        tag=tag,
        h=cfg["h"],
        roots=list(pot[1]),
        W1_eff=pot[2],
        gate=R25.GATE,
        sum_k2=SUMK2[j["part"]],
    )
    try:
        if os.path.exists(stage):
            Zs = np.load(stage, allow_pickle=True)
            M, done, chunks = Zs["M"], int(Zs["done"]), json.loads(str(Zs["chunks"]))
            row["resumed_at"] = done
            _, intended = F0.partition_of(j["part"])
            row["seed_gate"] = {"intended": intended, "resumed": True, "PASS": True}
        else:
            M, intended = F0.seed_partition(cfg, j["delta"], j["part"], w=pot[2])
            if np.abs(M[..., 0, 1:]).max() > 0:
                raise RuntimeError("the start field is not block-diagonal")
            row["seed_gate"] = R4.seed_gate(M, cfg, j["delta"], intended)
            if not row["seed_gate"]["PASS"]:
                row.update(status="SEED_GATE_FAIL", E=None, kick_label="NOT_RUN")
                row["wall_s"] = round(time.time() - t0, 1)
                return row
            np.savez_compressed(os.path.join(OUT_NPZ, tag + "_seed.npz"), M=M.astype(np.float64))
            done = 0
            chunks = [dict(R25.chunk_light(M, cfg, p, pot, mask, 0), note="the start field")]
            np.savez_compressed(stage, M=M, done=0, chunks=json.dumps(chunks))
        row["seed_split"] = energy_split(M, cfg, p, pot, j["kappa"])
        row["E_seed"] = chunks[0]["E"]
        M, done, chunks, verdict = descend_chunks(
            M, cfg, p, pot, mask, tag, j["cap"], chunks, done, stage, j["pinned"]
        )
        row.update(chunks=chunks, iters=done, gate_label=verdict)
        E_gate = chunks[-1]["E"]
        row["E_gate"] = E_gate
        row["E_err"] = abs(chunks[-1]["drop"]) if len(chunks) > 1 else None
        np.savez_compressed(os.path.join(OUT_NPZ, tag + "_gate.npz"), M=M.astype(np.float64))
        if ug is not None and ug.wrap_up():
            row["stop"] = "UPS wrap-up before the kick (resumable)"
        if verdict == "AT_GATE":
            Mk, seed, amp_eff, E_ref = R25.kick_field(
                M, mask, pot, tag, energy=lambda X: R25.spatial_fmax(X, cfg, p, pot, mask)[0]
            )
            E_k0 = R25.spatial_fmax(Mk, cfg, p, pot, mask)[0]
            kchunks = [dict(R25.chunk_light(Mk, cfg, p, pot, mask, 0), note="the kicked field")]
            kstage = os.path.join(OUT_NPZ, tag + "_kick_stage.npz")
            Mk, kdone, kchunks, kverdict = descend_chunks(
                Mk,
                cfg,
                p,
                pot,
                mask,
                tag + "_kick",
                10**9,
                kchunks,
                0,
                kstage,
                j["pinned"],
                R25.KICK_CHUNKS,
            )
            E_k = kchunks[-1]["E"]
            row["kick"] = {
                "amp_effective": amp_eff,
                "E_ref_reslaved": E_ref,
                "seed": seed,
                "E_kicked_start": E_k0,
                "E_after": E_k,
                "iters": kdone,
                "chunks": kchunks,
                "verdict": kverdict,
                "fmax_after": kchunks[-1]["fmax_spatial"],
            }
            if E_k < E_gate - R25.KICK_LOWER:
                row["kick_label"] = "SADDLE"
                row["E"] = E_k
                M = Mk
                row["reconverged"] = bool(kchunks[-1]["fmax_spatial"] < R25.GATE)
                row["E_err"] = abs(kchunks[-1]["drop"])
            elif E_k - E_gate < R25.KICK_RETURN * max(1.0, abs(E_gate)):
                row["kick_label"] = "STABLE"
                row["E"] = E_gate
                row["reconverged"] = True
            else:
                row["kick_label"] = "UNRESOLVED"
                row["E"] = E_gate
                row["reconverged"] = False
        else:
            row["kick_label"] = "FALLING"
            row["E"] = E_gate
            row["reconverged"] = False
        row["end_split"] = energy_split(M, cfg, p, pot, j["kappa"])
        row["own_reads"] = R4.own_reads(M, cfg, p, pot, j["delta"], j["part"])
        np.savez_compressed(os.path.join(OUT_NPZ, tag + ".npz"), M=M.astype(np.float64))
        row["status"] = "OK"
    except Exception as e:  # noqa: BLE001
        import traceback

        row.update(status="FAILED", stop=repr(e), traceback=traceback.format_exc(), E=None)
    row["wall_s"] = round(time.time() - t0, 1)
    log(
        f"DONE {tag} status {row['status']} gate {row.get('gate_label')} kick {row.get('kick_label')}"
        f" E {row.get('E')} split {row.get('end_split')} interior {(row.get('own_reads') or {}).get('interior_partition')}"
        f" wall {row['wall_s']}"
    )
    return row


def load_json():
    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as f:
            return json.load(f)
    return {"task": "M5.32 R27-2", "rows": {}, "collect": {}}


def save_json(J):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(J, f, indent=1, default=str)
    os.replace(tmp, OUT_JSON)


def run_pool(jobs, workers):
    workers = min(int(workers), 12)
    rows = load_json()["rows"]
    pending = [j for j in jobs if rows.get(job_tag(j), {}).get("status") != "OK"]
    log(f"pool: {len(pending)} jobs, {workers} workers")
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        futs = [ex.submit(run_job, j) for j in pending]
        for fut in as_completed(futs):
            row = fut.result()
            Jn = load_json()
            Jn["rows"][row["tag"]] = row
            save_json(Jn)
    log("pool done")


def certified(r):
    return (
        r.get("status") == "OK"
        and r.get("gate_label") == "AT_GATE"
        and (
            r.get("kick_label") == "STABLE"
            or (r.get("kick_label") == "SADDLE" and r.get("reconverged"))
        )
    )


def collect(rows=None):
    J = load_json()
    rows = J["rows"] if rows is None else rows
    prev = json.load(open(R26_4_JSON))["rows"] if os.path.exists(R26_4_JSON) else {}
    out = {"per_kappa": {}, "twins": {}}
    ka, kb = kappas()
    labels_k = {}
    for k in (ka, kb):
        tab = {}
        for part in ORDER:
            r = rows.get(
                job_tag(dict(part=part, delta=DELTA, w1s=W1S, n=32, L=48.0, kappa=k, pinned=True))
            )
            if r is None:
                tab[part] = None
                continue
            sp = r.get("end_split") or {}
            orr = r.get("own_reads") or {}
            tab[part] = {
                "sum_k2": SUMK2[part],
                "status": r.get("status"),
                "gate_label": r.get("gate_label"),
                "kick_label": r.get("kick_label"),
                "certified": certified(r),
                "iters": r.get("iters"),
                "E_total": sp.get("E_total"),
                "E_quartic": sp.get("E_quartic"),
                "E_kappa": sp.get("E_kappa"),
                "kappaE_kappa": sp.get("kappaE_kappa"),
                "E_err": r.get("E_err"),
                "E_seed_total": (r.get("seed_split") or {}).get("E_total"),
                "interior_partition": orr.get("interior_partition"),
                "largest_readable_sphere": orr.get("largest_readable_sphere"),
                "r9_partition": orr.get("r9_partition"),
                "carriers_r9": ((orr.get("spheres") or {}).get("R9") or {}).get(
                    "carriers_hu_rho_z_phi"
                ),
                "beta2": orr.get("biaxiality"),
                "virial": (orr.get("virial") or {}).get("E_u_over_3V"),
                "C_rigid": (orr.get("physical_generator") or {}).get("rigid"),
                "E_quartic_kappa0_R26_4": (prev.get(f"P{part}_d0.3_w25_n32_L48") or {}).get("E"),
            }
        cert = [q for q in ORDER if tab.get(q) and tab[q]["certified"]]
        gated = [q for q in ORDER if tab.get(q) and tab[q]["gate_label"] == "AT_GATE"]
        rec = {"rows": tab, "certified": cert, "at_gate": gated}
        if len(gated) < 4:
            rec["label"] = "CENSUS_INSUFFICIENT"
            labels_k[f"k{k:.3e}"] = "CENSUS_INSUFFICIENT"
        else:
            order_total = sorted(cert, key=lambda q: tab[q]["E_total"])
            order_quartic = sorted(cert, key=lambda q: tab[q]["E_quartic"])
            rec["certified_order_by_E_total"] = order_total
            rec["certified_order_by_E_quartic"] = order_quartic
            rec["sum_k2_order_of_certified"] = sorted(cert, key=lambda q: SUMK2[q])
            pairs_ok, pairs_inv, pairs_unres = [], [], []
            for i, a in enumerate(cert):
                for b in cert[i + 1 :]:
                    if SUMK2[a] == SUMK2[b]:
                        continue
                    lo, hi = (a, b) if SUMK2[a] < SUMK2[b] else (b, a)
                    dE = tab[hi]["E_total"] - tab[lo]["E_total"]
                    err = 2.0 * ((tab[lo]["E_err"] or 0.0) + (tab[hi]["E_err"] or 0.0))
                    entry = {
                        "lower_sum_k2": lo,
                        "higher_sum_k2": hi,
                        "E_hi_minus_E_lo": dE,
                        "two_errors": err,
                    }
                    if dE > err:
                        pairs_ok.append(entry)
                    elif dE < -err:
                        pairs_inv.append(entry)
                    else:
                        pairs_unres.append(entry)
            rec["pairs_in_order"] = pairs_ok
            rec["pairs_inverted"] = pairs_inv
            rec["pairs_unresolved"] = pairs_unres
            if len(cert) >= 4 and not pairs_inv and not pairs_unres:
                rec["label"] = "ORDER_K2_CONFIRMED_AT_THIS_KAPPA"
            elif pairs_inv:
                rec["label"] = "ORDER_K2_INVERTED_AT_THIS_KAPPA"
            else:
                rec["label"] = "ORDER_K2_UNRESOLVED_AT_THIS_KAPPA"
            labels_k[f"k{k:.3e}"] = rec["label"]
        out["per_kappa"][f"k{k:.3e}"] = rec
    vals = list(labels_k.values())
    if len(vals) == 2 and all(v == "ORDER_K2_CONFIRMED_AT_THIS_KAPPA" for v in vals):
        out["label"] = "ORDER_K2_CONFIRMED"
    elif len(vals) == 2 and all(v == "ORDER_K2_INVERTED_AT_THIS_KAPPA" for v in vals):
        out["label"] = "ORDER_K2_REFUTED"
    elif any(v == "CENSUS_INSUFFICIENT" for v in vals):
        out["label"] = "CENSUS_INSUFFICIENT"
    else:
        out["label"] = "ORDER_K2_UNRESOLVED"
    # the twins
    for part in TWINS:
        r = rows.get(
            job_tag(dict(part=part, delta=DELTA, w1s=W1S, n=32, L=48.0, kappa=0.0, pinned=False))
        )
        pr = prev.get(f"P{part}_d0.3_w25_n32_L48") or {}
        if r is None:
            out["twins"][part] = None
            continue
        orr = r.get("own_reads") or {}
        por = pr.get("own_reads") or {}
        rec = {
            "status": r.get("status"),
            "gate_label": r.get("gate_label"),
            "kick_label": r.get("kick_label"),
            "iters": r.get("iters"),
            "E_free": r.get("E"),
            "E_pinned_R26_4": pr.get("E"),
            "r9_partition_free": orr.get("r9_partition"),
            "r9_partition_pinned": por.get("r9_partition"),
            "interior_partition_free": orr.get("interior_partition"),
            "interior_partition_pinned": por.get("interior_partition"),
            "carriers_r9_free": ((orr.get("spheres") or {}).get("R9") or {}).get(
                "carriers_hu_rho_z_phi"
            ),
            "carriers_r9_pinned": ((por.get("spheres") or {}).get("R9") or {}).get(
                "carriers_hu_rho_z_phi"
            ),
            "carriers_r18_free": ((orr.get("spheres") or {}).get("R18") or {}).get(
                "carriers_hu_rho_z_phi"
            ),
            "boundary_partition_free_r21": ((orr.get("spheres") or {}).get("R21") or {}).get(
                "partition"
            ),
        }
        if r.get("status") == "OK":
            rec["label"] = (
                "FREE_SHELL_SAME_CLASS"
                if rec["r9_partition_free"] == rec["r9_partition_pinned"]
                and rec["r9_partition_free"] is not None
                else "FREE_SHELL_REORGANIZES"
            )
        out["twins"][part] = rec
    out["kappas"] = [ka, kb]
    J["collect"] = out
    save_json(J)
    return out


def smoke():
    global OUT_JSON, OUT_NPZ
    OUT_JSON = OUT_JSON.replace(".json", "_smoke.json")
    OUT_NPZ = OUT_NPZ + "_smoke"
    R25.OUT_NPZ = OUT_NPZ
    ka, kb = kappas()
    res = {"rows": {}}
    for j in (
        dict(part="3_1", delta=DELTA, w1s=W1S, n=16, L=48.0, cap=4, kappa=ka, pinned=True),
        dict(part="1_1_1_1", delta=DELTA, w1s=W1S, n=16, L=48.0, cap=4, kappa=0.0, pinned=False),
    ):
        r = run_job(j)
        res["rows"][r["tag"]] = {
            k: r.get(k) for k in ("status", "gate_label", "kick_label", "E", "end_split", "stop")
        }
        if r.get("traceback"):
            print(r["traceback"][-1500:])
    # the wrapper: the total energy in the descent equals quartic + kappa E_kappa
    cfg, p, pot = R4.cfg_pot(dict(n=16, L=48.0, delta=DELTA, w1s=W1S))
    M, _ = F0.seed_partition(cfg, DELTA, "2_2", w=pot[2])
    CS._R27_KAPPA = ka
    E_wrapped = CS.energy_grad(M, cfg, p, pot, 0.0, need_grad=False)[0]
    sp = energy_split(M, cfg, p, pot, ka)
    res["wrapper_rel_err"] = abs(E_wrapped - sp["E_total"]) / abs(sp["E_total"])
    # the wrapped gradient by finite differences on one spatial entry
    E0, G0, _ = CS.energy_grad(M, cfg, p, pot, 0.0)
    i = (8, 8, 8)
    eps = 1e-6
    Mp = M.copy()
    Mp[i][1, 2] += eps
    Mp[i][2, 1] += eps
    Mm = M.copy()
    Mm[i][1, 2] -= eps
    Mm[i][2, 1] -= eps
    fd = (
        CS.energy_grad(Mp, cfg, p, pot, 0.0, need_grad=False)[0]
        - CS.energy_grad(Mm, cfg, p, pot, 0.0, need_grad=False)[0]
    ) / (2 * eps)
    res["wrapper_grad_fd_rel_err"] = abs(fd - 2.0 * G0[i][1, 2]) / max(abs(fd), 1e-300)
    # a synthetic collect
    rows = {}
    for k in (ka, kb):
        for q, e in zip(ORDER, (8.0, 8.1, 8.2, 8.3, 8.4)):
            j = dict(part=q, delta=DELTA, w1s=W1S, n=32, L=48.0, kappa=k, pinned=True)
            rows[job_tag(j)] = dict(
                j,
                status="OK",
                gate_label="AT_GATE",
                kick_label="STABLE",
                iters=1,
                E=e,
                E_err=1e-4,
                end_split={
                    "E_total": e,
                    "E_quartic": e - 0.1,
                    "E_kappa": 0.1 / k,
                    "kappaE_kappa": 0.1,
                },
                seed_split={"E_total": e + 1},
                own_reads={"interior_partition": [4], "r9_partition": [4], "spheres": {}},
            )
    c = collect(rows)
    res["synthetic_label"] = c["label"]
    res["PASS"] = (
        all(v["status"] == "OK" for v in res["rows"].values())
        and res["wrapper_rel_err"] < 1e-12
        and res["wrapper_grad_fd_rel_err"] < 1e-5
        and res["synthetic_label"] == "ORDER_K2_CONFIRMED"
    )
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "rows"}, indent=1))
    for t, v in res["rows"].items():
        print(
            t, v["status"], v["gate_label"], v["kick_label"], v["E"], (v.get("stop") or "")[:200]
        )
    return res


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        smoke()
    elif mode == "run":
        run_pool(jobs_all(), sys.argv[2] if len(sys.argv) > 2 else 12)
        c = collect()
        print(
            json.dumps(
                {
                    "label": c["label"],
                    "twins": {k: (v or {}).get("label") for k, v in c["twins"].items()},
                },
                indent=1,
            )
        )
    elif mode == "collect":
        print(json.dumps(collect(), indent=1, default=str))
    elif mode == "jobs":
        for j in jobs_all():
            print(job_tag(j))
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
