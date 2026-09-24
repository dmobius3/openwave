"""M5.32 R25-K adversarial audit: the kick-and-continue census (scripts/m5_32_r25_k_kick.py).

INDEPENDENCE: the auditor did NOT read the audited script. Every check below is built from the
census DATA only (data/m5_32_r25_k_kick.json, the end fields under data/m5_32_r25_k/, the source
fields under data/m5_32_r22_1/ and data/m5_32_r23_1/) plus the shared instruments the two families
were built with (R21.polish, R20.energy_parts, CS.energy_grad, CS.solve_m00, B3.pin_shell), with
the energies cross-checked by an auditor-written formula (B3 curvature + own V4 + own linear term).

Checks (claim ids C1..C7, one dict each in the output JSON):
  C1 labels recomputed from the stored numbers AND from the stored end fields (own energies);
  C2 the 2 percent kick sizing, with the r23_1 reference re-slaved by the auditor;
  C3 block-diagonality of the end fields and the untouched pinned shells;
  C4 the r22_1 family: controls still descending, kicked arm never below;
  C5 the r23_1 family: the 7 STABLE + 8 UNRESOLVED ranges, the n48 controls;
  C6 own kick-and-continue on three r22_1 n32 rows (own seed, own sizing, R21.polish gate 0)
     plus one 10 percent kick;
  C7 own energy-density centroid, the core shift on the STABLE rows.

Run:  python3 scripts/m5_32_r25_k_audit.py [--iters 120] [--static-only]
Writes data/m5_32_r25_k_audit.json. Two worker subprocesses at most (a pool shares the machine);
each arm runs as `--arm tag arm seed frac iters outfile` in its own interpreter.
"""

import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
KICK_JSON = os.path.join(DATA, "m5_32_r25_k_kick.json")
KICK_NPZ = os.path.join(DATA, "m5_32_r25_k")
SRC_NPZ = {"r22_1": os.path.join(DATA, "m5_32_r22_1"), "r23_1": os.path.join(DATA, "m5_32_r23_1")}
OUT_JSON = os.path.join(DATA, "m5_32_r25_k_audit.json")
G = 8.0
TOL_SADDLE, TOL_STABLE, CORE_MAX_H = 1e-4, 1e-5, 0.2
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    argv = sys.argv
    sys.argv = [argv[0]]
    spec.loader.exec_module(mod)
    sys.argv = argv
    return mod


R21 = _load("m5_32_r21_1_runs", "m5_32_r21_1_runs.py")
CS = _load("m5_32_r23_1_cscan", "m5_32_r23_1_cscan.py")
R20, R0, B3 = R21.R20, R21.R0, R21.B3
W1 = B3.W1


# ================= setup per row =================
def setup(row):
    cfg = R21.cfg_of(row["n"], row["L"], G, row["delta"])
    p = R21.params_of(G, row["delta"])
    roots = R0.roots_of(cfg, degenerate=True)
    pot = ("v4", roots, W1 * row["w1s"])
    pinned = bool(row["pinned"]) if row["family"] == "r22_1" else True
    mask = R21.free_mask(cfg, pinned)
    return cfg, p, pot, float(row.get("c", 0.0) or 0.0), pinned, mask


def ref_energy(M, row, cfg, p, pot, c):
    """the family's own energy: R20.energy_parts (r22_1) or CS.energy_grad at c (r23_1)."""
    if row["family"] == "r22_1":
        return float(R20.energy_parts(M, cfg, p, pot)["E_total"])
    return float(CS.energy_grad(M, cfg, p, pot, c, need_grad=False)[0])


def own_energy(M, cfg, roots, w, c):
    """auditor-written total: B3 curvature + w sum_p (tr N^p - C_p)^2 - c sum_p a_p (tr N^p - C_p),
    a_p = b_(p-1)/p with prod over the distinct roots of (x - q) = sum b_k x^k."""
    h3 = cfg["h"] ** 3
    e_curv = float(B3.e_parts(M, cfg)[0])
    Me = M @ ETA
    P, t = None, []
    for k in range(4):
        P = Me if k == 0 else P @ Me
        t.append(np.einsum("...ii->...", P))
    cp = [sum(float(q) ** k for q in roots) for k in range(1, 5)]
    v4 = h3 * w * float(np.sum(sum((t[k] - cp[k]) ** 2 for k in range(4))))
    lin = 0.0
    if c != 0.0:
        b = np.poly(sorted(set(float(q) for q in roots)))[::-1]  # ascending powers
        a = [float(b[k]) / (k + 1) for k in range(4)]
        lin = -c * h3 * float(np.sum(sum(a[k] * (t[k] - cp[k]) for k in range(4))))
    return e_curv + v4 + lin


def label_of(E_before, E_ctrl, E_kick, core_shift_h):
    s = max(1.0, abs(E_before))
    if E_kick < E_ctrl - TOL_SADDLE * s:
        return "SADDLE"
    if abs(E_kick - E_ctrl) < TOL_STABLE * s and core_shift_h < CORE_MAX_H:
        return "STABLE"
    return "UNRESOLVED"


def density(M, cfg, pot, c):
    e = R20.density(M, cfg, pot)
    if c != 0.0:
        e = e + CS.lin_density(M, cfg, pot[1], c)
    return e


def centroid(e, cfg):
    """energy-density centroid over the cells above half the peak (own read)."""
    n, h = cfg["n"], cfg["h"]
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X = np.meshgrid(x, x, x, indexing="ij")
    sel = e >= 0.5 * e.max()
    w = e[sel]
    return np.array([float(np.sum(w * Xi[sel]) / np.sum(w)) for Xi in X])


def reslave(M, cfg, pot, c, mask):
    """the source with M_00 re-slaved on the free cells (the r23_1 family's reduced variable)."""
    Ms = M.copy()
    S = Ms[..., 1:, 1:]
    m00 = Ms[..., 0, 0].copy()
    m00[mask] = CS.solve_m00(S[mask], pot[1], pot[2], c, m00[mask])
    Ms[..., 0, 0] = m00
    return Ms


# ================= the static checks (C1..C5, C7) =================
def static_checks(rows):
    per = {}
    for tag, row in rows.items():
        r = {"family": row["family"], "label_stored": row["label"], "n": row.get("n")}
        per[tag] = r
        src = os.path.join(SRC_NPZ[row["family"]], tag + ".npz")
        ctrl_f = os.path.join(KICK_NPZ, tag + "_ctrl.npz")
        kick_f = os.path.join(KICK_NPZ, tag + "_kick.npz")
        r["files_present"] = [os.path.exists(f) for f in (src, ctrl_f, kick_f)]
        if row["label"] == "INSUFFICIENT" or "control" not in row:
            r["why_stored"] = row.get("why")
            r["stage_file_present"] = os.path.exists(
                os.path.join(SRC_NPZ[row["family"]], tag + "_stage.npz")
            )
            continue
        cfg, p, pot, c, pinned, mask = setup(row)
        ctl, kck = row["control"], row["kicked"]
        s = max(1.0, abs(row["E_before"]))
        # C1a: the label from the stored numbers
        r["label_from_json"] = label_of(
            row["E_before"], ctl["E_after"], kck["E_after"], kck["core_shift_h"]
        )
        r["rel_kick_minus_ctrl_json"] = (kck["E_after"] - ctl["E_after"]) / s
        r["ctrl_drop_json"] = ctl["E_start"] - ctl["E_after"]
        r["ctrl_drop_rel_E"] = r["ctrl_drop_json"] / abs(row["E_before"])
        r["ctrl_drop_rel_s"] = r["ctrl_drop_json"] / s
        r["kick_excess_json"] = kck["E_start"] - row["E_before"]
        r["kick_excess_over_target"] = r["kick_excess_json"] / (0.02 * s)
        if not all(r["files_present"]):
            r["note"] = "a field file is missing; the field checks are skipped"
            continue
        M0 = np.load(src)["M"]
        Mc = np.load(ctrl_f)["M"]
        Mk = np.load(kick_f)["M"]
        # C1b: own energies of the end fields
        E0, Ec, Ek = (ref_energy(M, row, cfg, p, pot, c) for M in (M0, Mc, Mk))
        O0, Oc, Ok = (own_energy(M, cfg, pot[1], pot[2], c) for M in (M0, Mc, Mk))
        r["E_src_ref"], r["E_ctrl_ref"], r["E_kick_ref"] = E0, Ec, Ek
        r["own_vs_ref_rel"] = max(
            abs(O0 - E0) / max(1, abs(E0)),
            abs(Oc - Ec) / max(1, abs(Ec)),
            abs(Ok - Ek) / max(1, abs(Ek)),
        )
        r["E_before_vs_src_rel"] = abs(row["E_before"] - E0) / max(1, abs(E0))
        r["ctrl_after_vs_field_rel"] = abs(ctl["E_after"] - Ec) / max(1, abs(Ec))
        r["kick_after_vs_field_rel"] = abs(kck["E_after"] - Ek) / max(1, abs(Ek))
        # C7: own centroid and core shift
        ec, ek = density(Mc, cfg, pot, c), density(Mk, cfg, pot, c)
        cc, ck = centroid(ec, cfg), centroid(ek, cfg)
        r["core_ctrl_own"], r["core_kick_own"] = cc.tolist(), ck.tolist()
        r["core_shift_h_own"] = float(np.linalg.norm(ck - cc) / cfg["h"])
        r["core_shift_h_stored"] = kck["core_shift_h"]
        r["core_ctrl_vs_stored_h"] = float(np.linalg.norm(cc - np.array(ctl["core"])) / cfg["h"])
        # the same centroid on the density WITHOUT the linear member (a definition probe)
        cc2 = centroid(R20.density(Mc, cfg, pot), cfg)
        ck2 = centroid(R20.density(Mk, cfg, pot), cfg)
        r["core_shift_h_own_nolin"] = float(np.linalg.norm(ck2 - cc2) / cfg["h"])
        r["core_ctrl_nolin_vs_stored_h"] = float(
            np.linalg.norm(cc2 - np.array(ctl["core"])) / cfg["h"]
        )
        r["label_from_fields"] = label_of(E0, Ec, Ek, r["core_shift_h_own"])
        # C3: block-diagonal + pinned shell
        r["max_abs_M0i"] = [float(np.abs(M[..., 0, 1:]).max()) for M in (M0, Mc, Mk)]
        shell = B3.pin_shell(cfg["n"], cfg["h"])
        r["pinned"] = pinned
        r["shell_max_dev"] = [float(np.abs(M[shell] - M0[shell]).max()) for M in (Mc, Mk)]
        r["shell_free_cells_moved"] = [
            float(np.abs(M[~shell] - M0[~shell]).max()) for M in (Mc, Mk)
        ]
        # C2 (r23_1): the re-slaved reference
        if row["family"] == "r23_1":
            Ms = reslave(M0, cfg, pot, c, mask)
            Er = ref_energy(Ms, row, cfg, p, pot, c)
            r["E_src_reslaved"] = Er
            r["reslave_shift_rel"] = (Er - E0) / max(1, abs(E0))
            r["kick_excess_over_target_reslaved"] = (kck["E_start"] - Er) / (
                0.02 * max(1, abs(Er))
            )
            r["m00_reslave_max_change"] = float(np.abs(Ms[..., 0, 0] - M0[..., 0, 0]).max())
        log(
            f"{tag:46s} {row['label']:10s} json->{r['label_from_json']:10s} fields->{r['label_from_fields']:10s} "
            f"dE_after {r['ctrl_after_vs_field_rel']:.1e}/{r['kick_after_vs_field_rel']:.1e} own/ref {r['own_vs_ref_rel']:.1e} "
            f"core {r['core_shift_h_own']:.3f} vs {kck['core_shift_h']:.3f} M0i {max(r['max_abs_M0i']):.1e} shell {max(r['shell_max_dev']):.1e}"
        )
    return per


# ================= C6: the auditor's own kick-and-continue =================
def kick_field(M, cfg, p, pot, mask, seed, frac):
    """a symmetric white kick on the spatial 3x3 block of the free cells, sized so the energy sits
    frac x max(1, |E0|) above the source (secant iteration on the amplitude); M_0i stays 0, M_00 fixed.
    """
    rng = np.random.default_rng(seed)
    N = rng.standard_normal(M.shape[:3] + (3, 3))
    N = 0.5 * (N + N.swapaxes(-1, -2))
    N[~mask] = 0.0
    E0, _, info0 = R20.energy_grad(M, cfg, p, pot, need_grad=False)
    target = frac * max(1.0, abs(E0))

    def excess(a):
        Mk = M.copy()
        Mk[..., 1:, 1:] += a * N
        E, _, info = R20.energy_grad(Mk, cfg, p, pot, need_grad=False)
        return (E - E0) if (np.isfinite(E) and info["ok"]) else np.inf, Mk

    a, hist = 1e-3, []
    for _ in range(25):
        f, Mk = excess(a)
        hist.append((a, f))
        if not np.isfinite(f):
            a *= 0.5
            continue
        if abs(f / target - 1.0) < 2e-3:
            break
        a *= float(np.sqrt(target / max(f, 1e-300)))
    return Mk, {
        "E0": float(E0),
        "excess": float(f),
        "target": target,
        "amp": a,
        "n_size": len(hist),
    }


def arm_job(args):
    tag, row, arm, seed, frac, iters = args
    cfg, p, pot, c, pinned, mask = setup(row)
    M0 = np.load(os.path.join(SRC_NPZ[row["family"]], tag + ".npz"))["M"]
    out = {"tag": tag, "arm": arm, "frac": frac, "seed": seed, "iters_req": iters}
    t0 = time.time()
    if arm == "ctrl":
        Ms = M0
    else:
        Ms, siz = kick_field(M0, cfg, p, pot, mask, seed, frac)
        out["sizing"] = siz
        out["max_abs_M0i_kicked"] = float(np.abs(Ms[..., 0, 1:]).max())
    Mend, rec = R21.polish(
        Ms, cfg, p, pot, pinned, f"audit-{tag}-{arm}", max_iter=iters, gate=0.0, log_every=10
    )
    out.update(
        E_start=rec["E_in"],
        E_end=rec["E_end"],
        fmax_end=rec["fmax_end"],
        iters_done=rec["iters"],
        n_eval=rec["n_eval"],
        scipy_message=rec.get("scipy_message"),
        trace=[(t["it"], t["E"]) for t in rec["trace"]],
        wall_s=round(time.time() - t0, 1),
    )
    e = density(Mend, cfg, pot, c)
    out["core_end"] = centroid(e, cfg).tolist()
    out["shell_max_dev"] = float(
        np.abs(Mend[B3.pin_shell(cfg["n"], cfg["h"])] - M0[B3.pin_shell(cfg["n"], cfg["h"])]).max()
    )
    return out


def dynamic_checks(rows, iters, workers=2):
    picks = ["perm_free_d0.3_w25_n32_L48", "bia_pin_d0.3_w25_n32_L48", "rad_pin_d0.3_w1_n32_L48"]
    seeds = {t: 20260923 + i for i, t in enumerate(picks)}
    jobs = []
    for t in picks:
        jobs.append((t, rows[t], "ctrl", seeds[t], 0.02, iters))
        jobs.append((t, rows[t], "kick", seeds[t], 0.02, iters))
    jobs.append((picks[0], rows[picks[0]], "kick10", seeds[picks[0]] + 7, 0.10, iters))
    tmp = tempfile.mkdtemp(prefix="r25k_audit_")
    pending, running, res = list(jobs), [], []
    while pending or running:
        while pending and len(running) < workers:
            j = pending.pop(0)
            outf = os.path.join(tmp, f"{j[0]}_{j[2]}.json")
            cmd = [
                sys.executable,
                os.path.abspath(__file__),
                "--arm",
                j[0],
                j[2],
                str(j[3]),
                str(j[4]),
                str(j[5]),
                outf,
            ]
            running.append((subprocess.Popen(cmd), outf, j))
        time.sleep(2.0)
        still = []
        for pr, outf, j in running:
            if pr.poll() is None:
                still.append((pr, outf, j))
            elif pr.returncode == 0 and os.path.exists(outf):
                res.append(json.load(open(outf)))
                log(
                    f"arm done {j[0]} {j[2]} E_end {res[-1]['E_end']:.6f} in {res[-1]['wall_s']} s"
                )
            else:
                res.append({"tag": j[0], "arm": j[2], "failed": True, "returncode": pr.returncode})
                log(f"arm FAILED {j[0]} {j[2]} rc {pr.returncode}")
        running = still
    shutil.rmtree(tmp, ignore_errors=True)
    by = {}
    for r in res:
        by.setdefault(r["tag"], {})[r["arm"]] = r
    return by


# ================= verdicts =================
def verdict_static(rows, per):
    checks = []
    ok_rows = {t: r for t, r in per.items() if "label_from_fields" in r}
    r22 = {t: r for t, r in ok_rows.items() if r["family"] == "r22_1"}
    r23 = {t: r for t, r in ok_rows.items() if r["family"] == "r23_1"}
    # ---- C1
    counts = {}
    for r in per.values():
        counts[r["label_stored"]] = counts.get(r["label_stored"], 0) + 1
    mism_json = [t for t, r in ok_rows.items() if r["label_from_json"] != r["label_stored"]]
    mism_field = [t for t, r in ok_rows.items() if r["label_from_fields"] != r["label_stored"]]
    e_dis = {
        t: (r["ctrl_after_vs_field_rel"], r["kick_after_vs_field_rel"], r["E_before_vs_src_rel"])
        for t, r in ok_rows.items()
        if max(r["ctrl_after_vs_field_rel"], r["kick_after_vs_field_rel"]) > 1e-8
    }
    own_max = max(r["own_vs_ref_rel"] for r in ok_rows.values())
    e_before_max = max(r["E_before_vs_src_rel"] for r in ok_rows.values())
    insuff = [t for t, r in per.items() if r["label_stored"] == "INSUFFICIENT"]
    c1_ok = (
        counts == {"UNRESOLVED": 25, "STABLE": 7, "INSUFFICIENT": 1}
        and not mism_json
        and not mism_field
        and not e_dis
    )
    checks.append(
        {
            "id": "C1",
            "claim": "labels follow the stated rule; counts 25 / 7 / 1 / 0; stored E_after = own energy of the end field",
            "method": "relabel every row from the JSON numbers and from own energies of the stored end fields (R20.energy_parts for r22_1, CS.energy_grad at c for r23_1, cross-checked by an auditor-written formula)",
            "numbers": {
                "counts_stored": counts,
                "rows_with_fields": len(ok_rows),
                "label_mismatch_from_json": mism_json,
                "label_mismatch_from_fields": mism_field,
                "E_after_disagree_gt_1e-8": e_dis,
                "max_E_after_rel_disagreement": max(
                    max(r["ctrl_after_vs_field_rel"], r["kick_after_vs_field_rel"])
                    for r in ok_rows.values()
                ),
                "max_E_before_vs_source_field_rel": e_before_max,
                "own_formula_vs_family_energy_max_rel": own_max,
                "insufficient_rows": {t: per[t] for t in insuff},
            },
            "verdict": "CONFIRMED" if c1_ok else "REFUTED",
        }
    )
    # ---- C2
    ex22 = {t: r["kick_excess_over_target"] for t, r in r22.items()}
    ex23_raw = {t: r["kick_excess_over_target"] for t, r in r23.items()}
    ex23_res = {t: r["kick_excess_over_target_reslaved"] for t, r in r23.items()}
    resl = {t: r["reslave_shift_rel"] for t, r in r23.items()}
    tol = 0.01
    c2_22 = all(abs(v - 1) < tol for v in ex22.values())
    c2_23 = all(abs(v - 1) < tol for v in ex23_res.values())
    c2_23raw = all(abs(v - 1) < tol for v in ex23_raw.values())
    checks.append(
        {
            "id": "C2",
            "claim": "kicked start energy = source energy + 0.02 max(1, |E_ref|), E_ref re-slaved for r23_1",
            "method": "(kicked.E_start - E_ref) / (0.02 max(1,|E_ref|)) per row, E_ref = E_before for r22_1, own re-slaved source energy (CS.solve_m00 on the free cells) for r23_1",
            "numbers": {
                "r22_1_excess_over_target_min_max": [min(ex22.values()), max(ex22.values())],
                "r23_1_excess_over_target_vs_E_before_min_max": [
                    min(ex23_raw.values()),
                    max(ex23_raw.values()),
                ],
                "r23_1_excess_over_target_vs_reslaved_min_max": [
                    min(ex23_res.values()),
                    max(ex23_res.values()),
                ],
                "r23_1_reslave_shift_rel_min_max": [min(resl.values()), max(resl.values())],
                "r23_1_E_before_is_reslaved": all(
                    abs(rows[t]["E_before"] - r["E_src_reslaved"])
                    < 1e-9 * max(1, abs(r["E_src_reslaved"]))
                    for t, r in r23.items()
                ),
                "r23_1_E_before_is_raw_source": all(
                    r["E_before_vs_src_rel"] < 1e-9 for r in r23.values()
                ),
                "rows_off_target_gt_1pct": {
                    t: {
                        "excess_over_target": v,
                        "E_before": rows[t]["E_before"],
                        "kick_excess_abs": rows[t]["kicked"]["E_start"] - rows[t]["E_before"],
                        "w1s": rows[t]["w1s"],
                    }
                    for t, v in {**ex22, **ex23_res}.items()
                    if abs(v - 1) >= tol
                },
                "rows_with_kick_gt_10x_E": {
                    t: (rows[t]["kicked"]["E_start"] - rows[t]["E_before"])
                    / abs(rows[t]["E_before"])
                    for t in ok_rows
                    if (rows[t]["kicked"]["E_start"] - rows[t]["E_before"])
                    > 10 * abs(rows[t]["E_before"])
                },
                "per_row_r22_1": ex22,
                "per_row_r23_1_reslaved": ex23_res,
            },
            "verdict": (
                "CONFIRMED"
                if (c2_22 and c2_23)
                else ("QUALIFIED" if (c2_22 and c2_23raw) else "REFUTED")
            ),
            "reading": "the r23_1 sources are already M_00-slaved (re-slaving moves nothing), so the raw and re-slaved references coincide; a row whose excess misses the 2 percent target by more than 1 percent is listed in rows_off_target_gt_1pct; rows_with_kick_gt_10x_E lists where the max(1, |E|) floor made the kick many times the state's own energy",
        }
    )
    # ---- C3
    m0i = max(max(r["max_abs_M0i"]) for r in ok_rows.values())
    shell_pinned = {t: max(r["shell_max_dev"]) for t, r in ok_rows.items() if r["pinned"]}
    shell_free = {t: max(r["shell_max_dev"]) for t, r in ok_rows.items() if not r["pinned"]}
    checks.append(
        {
            "id": "C3",
            "claim": "end fields block-diagonal (M_0i = 0) and the pinned shells untouched",
            "method": "max |M_0i| over source, ctrl, kick fields; max |M_end - M_src| on B3.pin_shell cells, split by pinned / free rows",
            "numbers": {
                "max_abs_M0i_all_fields": m0i,
                "pinned_rows_shell_max_dev": max(shell_pinned.values()),
                "free_rows_shell_max_dev_min_max": (
                    [min(shell_free.values()), max(shell_free.values())] if shell_free else None
                ),
                "n_pinned_rows": len(shell_pinned),
                "n_free_rows": len(shell_free),
            },
            "verdict": (
                "CONFIRMED" if (m0i == 0.0 and max(shell_pinned.values()) <= 1e-15) else "REFUTED"
            ),
        }
    )
    # ---- C4
    drops = {t: r["ctrl_drop_json"] for t, r in r22.items()}
    drel = {t: r["ctrl_drop_rel_E"] for t, r in r22.items()}
    kabove = {t: r["rel_kick_minus_ctrl_json"] for t, r in r22.items()}
    kabove_E = {
        t: (rows[t]["kicked"]["E_after"] - rows[t]["control"]["E_after"])
        / abs(rows[t]["E_before"])
        for t in r22
    }
    all_unres = all(r["label_stored"] == "UNRESOLVED" for r in r22.values())
    never_below = all(v > 0 for v in kabove.values())
    stated_abs = 0.0046 <= min(drops.values()) and max(drops.values()) <= 0.196 + 5e-4
    stated_rel = 0.002 <= min(drel.values()) and max(drel.values()) <= 0.04
    stated_k = 1e-4 <= min(kabove.values()) and max(kabove.values()) <= 1.4e-2 + 5e-4
    checks.append(
        {
            "id": "C4",
            "claim": "every r22_1 row UNRESOLVED; control drop 0.2 to 4 percent of E (0.0046 to 0.196 abs); kicked end 1e-4 to 1.4e-2 relative above control, never below",
            "method": "from the stored numbers: ctrl.E_start - ctrl.E_after, its ratio to |E_before| and to max(1,|E_before|), (E_kick - E_ctrl)/max(1,|E_before|)",
            "numbers": {
                "n_rows": len(r22),
                "all_unresolved": all_unres,
                "kick_never_below_ctrl": never_below,
                "ctrl_drop_abs_min_max": [min(drops.values()), max(drops.values())],
                "ctrl_drop_abs_argmin": min(drops, key=drops.get),
                "ctrl_drop_rel_E_min_max": [min(drel.values()), max(drel.values())],
                "ctrl_drop_rel_E_argmax": max(drel, key=drel.get),
                "ctrl_drop_rel_E_argmin": min(drel, key=drel.get),
                "kick_above_ctrl_rel_s_min_max": [min(kabove.values()), max(kabove.values())],
                "kick_above_ctrl_rel_s_argmin": min(kabove, key=kabove.get),
                "kick_above_ctrl_rel_E_min_max": [min(kabove_E.values()), max(kabove_E.values())],
                "stated_ranges_hold": {"abs": stated_abs, "rel": stated_rel, "kick": stated_k},
                "per_row": {
                    t: {"drop": drops[t], "drop_rel_E": drel[t], "kick_above_rel_s": kabove[t]}
                    for t in r22
                },
            },
            "verdict": (
                "CONFIRMED"
                if (all_unres and never_below and stated_abs and stated_rel and stated_k)
                else ("QUALIFIED" if (all_unres and never_below) else "REFUTED")
            ),
        }
    )
    # ---- C5
    stab = {
        t: abs(r["rel_kick_minus_ctrl_json"])
        for t, r in r23.items()
        if r["label_stored"] == "STABLE"
    }
    unres = {
        t: r["rel_kick_minus_ctrl_json"]
        for t, r in r23.items()
        if r["label_stored"] == "UNRESOLVED"
    }
    n48 = {t: (r["label_stored"], r["ctrl_drop_json"]) for t, r in r23.items() if r["n"] == 48}
    below = [t for t, r in r23.items() if r["rel_kick_minus_ctrl_json"] < 0]
    c5_ok = (
        len(stab) == 7
        and len(unres) == 8
        and 1e-13 <= min(stab.values()) <= 2e-13
        and 6.5e-6 <= max(stab.values()) <= 7e-6
        and 1.2e-5 <= min(unres.values()) <= 1.4e-5
        and 3.6e-4 <= max(unres.values()) <= 3.8e-4
        and all(lab == "UNRESOLVED" for lab, _ in n48.values())
        and all(5e-3 <= d <= 2e-2 for _, d in n48.values())
        and not below
    )
    checks.append(
        {
            "id": "C5",
            "claim": "r23_1: 7 STABLE (1.7e-13 to 6.8e-6 rel), 8 UNRESOLVED (1.3e-5 to 3.7e-4 above control), n48 rows UNRESOLVED with controls dropping ~1e-2, none below control",
            "method": "from the stored numbers, (E_kick - E_ctrl)/max(1,|E_before|) per r23_1 row, the n48 control drops",
            "numbers": {
                "n_stable": len(stab),
                "n_unresolved": len(unres),
                "stable_rel_min_max": [min(stab.values()), max(stab.values())] if stab else None,
                "unresolved_rel_min_max": (
                    [min(unres.values()), max(unres.values())] if unres else None
                ),
                "n48_rows": n48,
                "rows_below_control": below,
                "stable_tol_over_abs_E": {
                    t: TOL_STABLE * max(1.0, abs(rows[t]["E_before"])) / abs(rows[t]["E_before"])
                    for t in stab
                },
                "stable_kick_return_over_abs_E": {
                    t: abs(rows[t]["kicked"]["E_after"] - rows[t]["control"]["E_after"])
                    / abs(rows[t]["E_before"])
                    for t in stab
                },
                "stable_ctrl_iters_done": {t: rows[t]["control"]["iters_done"] for t in stab},
                "stable_kick_iters_done": {t: rows[t]["kicked"]["iters_done"] for t in stab},
                "note_rel_kick_minus_ctrl_can_be_negative_tiny": {
                    t: r["rel_kick_minus_ctrl_json"]
                    for t, r in r23.items()
                    if r["rel_kick_minus_ctrl_json"] < 0
                },
            },
            "verdict": (
                "CONFIRMED"
                if (
                    c5_ok
                    and all(
                        TOL_STABLE * max(1.0, abs(rows[t]["E_before"]))
                        < 1e-2 * abs(rows[t]["E_before"])
                        for t in stab
                    )
                )
                else ("QUALIFIED" if not below else "REFUTED")
            ),
            "reading": "the label tolerance 1e-5 x max(1, |E|) is an ABSOLUTE 1e-5 on the rows with |E| < 1; on the delta 0.9224 rows (E ~ 1e-4 to 6e-4) that is percent-level of the state's own energy, so STABLE there is a statement at the absolute scale (stable_tol_over_abs_E, stable_kick_return_over_abs_E)",
        }
    )
    # ---- C7
    shifts = {t: (r["core_shift_h_own"], r["core_shift_h_stored"]) for t, r in ok_rows.items()}
    stable_shift = {t: v for t, v in shifts.items() if per[t]["label_stored"] == "STABLE"}
    max_dev = max(abs(a - b) for a, b in shifts.values())
    ctrl_dev = max(r["core_ctrl_vs_stored_h"] for r in ok_rows.values())
    checks.append(
        {
            "id": "C7",
            "claim": "core shift (half-peak energy-density centroid) under 0.2 h on every STABLE row",
            "method": "own density (R20.density + own linear density) and own half-peak centroid on the stored ctrl and kick end fields, shift in units of h",
            "numbers": {
                "stable_rows_shift_own_vs_stored": stable_shift,
                "max_stable_shift_own_h": max(v[0] for v in stable_shift.values()),
                "max_abs_dev_own_vs_stored_all_rows_h": max_dev,
                "max_ctrl_centroid_dev_vs_stored_h": ctrl_dev,
                "max_abs_dev_nolin_vs_stored_all_rows_h": max(
                    abs(r["core_shift_h_own_nolin"] - r["core_shift_h_stored"])
                    for r in ok_rows.values()
                ),
                "max_ctrl_centroid_nolin_dev_vs_stored_h": max(
                    r["core_ctrl_nolin_vs_stored_h"] for r in ok_rows.values()
                ),
                "stable_rows_shift_nolin_h": {
                    t: per[t]["core_shift_h_own_nolin"] for t in stable_shift
                },
                "all_rows_shift_own_h": {t: v[0] for t, v in shifts.items()},
            },
            "verdict": (
                "CONFIRMED"
                if (max(v[0] for v in stable_shift.values()) < CORE_MAX_H and max_dev < 0.02)
                else (
                    "QUALIFIED"
                    if max(v[0] for v in stable_shift.values()) < CORE_MAX_H
                    else "REFUTED"
                )
            ),
            "reading": "max_abs_dev_nolin_vs_stored_all_rows_h = 0 means the stored core is the centroid of the curvature + V4 density WITHOUT the linear member (the c term); the auditor's centroid on the full density differs on the c != 0 rows but stays under 0.2 h on every STABLE row",
        }
    )
    return checks


def verdict_dynamic(rows, dyn):
    per, any_below, any_below_trace = {}, False, False
    for tag, arms in dyn.items():
        arms = {a: r for a, r in arms.items() if not r.get("failed")}
        c, k = arms.get("ctrl"), arms.get("kick")
        s = max(1.0, abs(rows[tag]["E_before"]))
        d = {"E_before": rows[tag]["E_before"]}
        if c and k:
            d["ctrl"] = {
                kk: c[kk]
                for kk in (
                    "E_start",
                    "E_end",
                    "fmax_end",
                    "iters_done",
                    "n_eval",
                    "wall_s",
                    "scipy_message",
                )
            }
            d["kick"] = {
                kk: k[kk]
                for kk in (
                    "E_start",
                    "E_end",
                    "fmax_end",
                    "iters_done",
                    "n_eval",
                    "wall_s",
                    "scipy_message",
                    "sizing",
                    "max_abs_M0i_kicked",
                    "shell_max_dev",
                )
            }
            d["kick_excess_over_target"] = k["sizing"]["excess"] / k["sizing"]["target"]
            d["ctrl_drop"] = c["E_start"] - c["E_end"]
            d["kick_minus_ctrl_rel_s"] = (k["E_end"] - c["E_end"]) / s
            core_shift = float(
                np.linalg.norm(np.array(k["core_end"]) - np.array(c["core_end"])) / rows[tag]["h"]
            )
            d["core_shift_h"] = core_shift
            d["label_own"] = label_of(rows[tag]["E_before"], c["E_end"], k["E_end"], core_shift)
            ct, kt = dict(c["trace"]), dict(k["trace"])
            common = sorted(set(ct) & set(kt))
            d["trace_kick_minus_ctrl_rel_s"] = [(it, (kt[it] - ct[it]) / s) for it in common]
            d["kick_ever_below_ctrl_in_trace"] = any(kt[it] < ct[it] for it in common)
            any_below |= k["E_end"] < c["E_end"] - TOL_SADDLE * s
            any_below_trace |= d["kick_ever_below_ctrl_in_trace"]
        b = arms.get("kick10")
        if b and c:
            d["kick10"] = {
                kk: b[kk]
                for kk in (
                    "E_start",
                    "E_end",
                    "fmax_end",
                    "iters_done",
                    "n_eval",
                    "wall_s",
                    "sizing",
                    "shell_max_dev",
                )
            }
            d["kick10_excess_over_target"] = b["sizing"]["excess"] / b["sizing"]["target"]
            d["kick10_minus_ctrl_rel_s"] = (b["E_end"] - c["E_end"]) / s
            cs10 = float(
                np.linalg.norm(np.array(b["core_end"]) - np.array(c["core_end"])) / rows[tag]["h"]
            )
            d["kick10_core_shift_h"] = cs10
            d["kick10_label_own"] = label_of(rows[tag]["E_before"], c["E_end"], b["E_end"], cs10)
            d["kick10_drop_vs_own_start"] = b["E_start"] - b["E_end"]
            any_below |= b["E_end"] < c["E_end"] - TOL_SADDLE * s
        per[tag] = d
    labels = {t: d.get("label_own") for t, d in per.items()}
    return {
        "id": "C6",
        "claim": "no kicked descent found a state below its control (no SADDLE); the saddle hypothesis unsupported at 2 percent / 500 its; the R22-1 rows cannot be called minima",
        "method": "own white symmetric kick on the spatial block (own seed, secant-sized to 2 percent of max(1,|E|)), both arms continued by R21.polish (L-BFGS, gate 0) for the same short budget; one 10 percent kick on the first row",
        "numbers": {
            "per_row": per,
            "own_labels_2pct": labels,
            "any_kick_end_below_ctrl_by_saddle_tol": any_below,
            "any_kick_below_ctrl_at_a_logged_iteration": any_below_trace,
        },
        "verdict": "CONFIRMED" if not any_below else "REFUTED",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iters", type=int, default=120)
    ap.add_argument("--static-only", action="store_true")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--arm", nargs=6, metavar=("TAG", "ARM", "SEED", "FRAC", "ITERS", "OUT"))
    a = ap.parse_args()
    if a.arm:
        tag, arm, seed, frac, iters, outf = a.arm
        rows = json.load(open(KICK_JSON))["rows"]
        r = arm_job((tag, rows[tag], arm, int(seed), float(frac), int(iters)))
        with open(outf, "w") as f:
            json.dump(r, f, default=float)
        return
    J = json.load(open(KICK_JSON))
    rows = J["rows"]
    log(f"{len(rows)} rows; stored counts {J['collect']['counts']}")
    per = static_checks(rows)
    checks = verdict_static(rows, per)
    out = {
        "task": "M5.32 R25-K audit",
        "independence": "the auditor did not read scripts/m5_32_r25_k_kick.py; checks built from the census data and the shared family instruments only",
        "checks": checks,
        "per_row": per,
    }

    def finish(out):
        c = {}
        for ch in out["checks"]:
            c[ch["verdict"]] = c.get(ch["verdict"], 0) + 1
        out["summary"] = c
        out["wall_s"] = round(time.time() - T0, 1)
        with open(OUT_JSON, "w") as f:
            json.dump(out, f, indent=1, default=float)
        log(f"summary {c}; wrote {OUT_JSON}")

    if a.static_only:
        prev = json.load(open(OUT_JSON)) if os.path.exists(OUT_JSON) else {}
        c6 = [c for c in prev.get("checks", []) if c["id"] == "C6" and c["verdict"] != "NOT_RUN"]
        if c6:
            out["checks"].insert(5, c6[0])
            out["dynamic_raw"] = prev.get("dynamic_raw")
        else:
            out["checks"].append({"id": "C6", "verdict": "NOT_RUN", "reason": "--static-only"})
        finish(out)
        return
    finish(out)
    log(f"C6: own kick-and-continue, {a.iters} iterations per arm, {a.workers} workers")
    dyn = dynamic_checks(rows, a.iters, a.workers)
    out["checks"].insert(5, verdict_dynamic(rows, dyn))
    out["dynamic_raw"] = dyn
    finish(out)
    for ch in out["checks"]:
        log(f"{ch['id']} {ch['verdict']}")


if __name__ == "__main__":
    main()
