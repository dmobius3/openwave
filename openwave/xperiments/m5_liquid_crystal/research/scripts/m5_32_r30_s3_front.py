"""M5.32 R30: S3, the observable (R30-2): the fine-grid front over the search history, the finite-range omega coupling,
the merged 4He, the LP prices at each cell's minimum, the Pareto fronts, the pre-registered gate word and the plot.
Runs after `m5_32_r30_s3_search.py` (controls, de, slp); consumes the author's instrument through that module.

EQUATIONS FIRST
---------------
Candidates: over the feasible coarse history (gate <= 5.02, K_gamma >= 0.05, eta <= 0.9, proxy >= 5.8), per cell the
    E15 minimum, the two-dimensional non-dominated sets of (E15, |r2_EV - 0.823|), (E15, |r_n2 + 0.1161|),
    (E15, |r2_MV - 0.73|) and (E15, -proxy), the ten lowest E15 states, and every feasible SLP endpoint within 1.0 point
    of the best; states deduplicated at 1e-6 in theta.
Fine grid: V54b's continuation 500 / 12 -> 800 / 18 -> 1200 / 24 (cold at the first stage), the fifteen errors and E15
    on the fine grid, the held-out radii, M and z0 (the branch), the gate and K_gamma on the fine profile.
Finite-range omega coupling (V48b `omega_finite`, unchanged): rho_B(r) = -sin^2 F F' / (2 pi^2 r^2), its radial Fourier
    transform rho(k), the factor  fac = int k^2 rho^2 m^2 / (k^2 + m^2) dk / int k^2 rho^2 dk  (m = m_omega / hbar c),
    g_omega^2 / 4 pi (finite) = proxy / fac.
Merged 4He (V48b `merged`, unchanged): the rational-map cube of V43 at B = 1 and B = 4 with the covariant sextic folded
    through the map moments, E4 / 4 E1 and 4 M_N - E4 against the measured binding 28.3 MeV.
LP prices (V55 `prices`, unchanged) at each cell's fine-grid minimum: the dual multipliers lambda_r of the local minimax
    LP; a row is ACTIVE when lambda_r > 1e-3 (V55's own print threshold).
Words (the task doc, § R30-2, pre-registered):
    THE CLAIM FALSIFIED  a feasible state with g_omega^2 / 4 pi (finite) >= 7.9, fine E15 < 8.0 and r2_EV > 0, r_n2 < 0,
                         r2_MV > 0 (a state at 8.0 with a wrong-sign radius is reported, not counted);
    THE CLAIM CONFIRMED  the feasible fine-grid minimum over the four cells within 0.5 of 9.2 with Delta and gA active;
    NEITHER              otherwise;
    the wording limit: no outcome is a lower bound; "no counterexample found over the declared space".
Out: data/m5_32_r30_s3.json (front, prices, gate), plots/m5_32_r30_s3_front.png; the per-candidate fine rows in
     data/m5_32_r30_s3/front_candidates.json (tracked, small).
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import argparse
import contextlib
import io
import json
import multiprocessing as mp
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import m5_32_r30_s3_search as S                                          # noqa: E402  (loads the kit from $M5_32_S3_KIT)

PLOT = os.path.join(HERE, "..", "plots", "m5_32_r30_s3_front.png")
HELD = {"rEV2": 0.823, "rn2": -0.1161, "rMV2": 0.73}
HE4_BINDING = 28.3
GW_GATE, E15_FALSIFY, CONFIRM_CENTER, CONFIRM_HALF, ACTIVE = 7.9, 8.0, 9.2, 0.5, 1e-3
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:8.1f}s] {msg}", flush=True)


def pareto(points):
    """indices of the non-dominated points of a list of (a, b) pairs, both minimised."""
    idx = sorted(range(len(points)), key=lambda i: (points[i][0], points[i][1]))
    out, best_b = [], np.inf
    for i in idx:
        if points[i][1] < best_b:
            out.append(i); best_b = points[i][1]
    return out


def candidates(hist, slp_rows, cap_front=25):
    ok = hist["ok"] & hist["feasible"]
    th = hist["theta"]; cell = hist["cell"]
    sel = {}

    def add(i, why):
        key = (str(cell[i]), tuple(np.round(th[i], 6)))
        sel.setdefault(key, {"theta": th[i].tolist(), "cell": str(cell[i]), "why": [], "coarse": {k: float(hist[k][i]) for k in ("E15", "merit", "proxy", "gate", "Kgamma", "M", "z0", "rEV2", "rn2", "rMV2", "gA", "Delta")}})
        sel[key]["why"].append(why)
    for c in S.CELLS:
        m = np.where(ok & (cell == c))[0]
        if len(m):
            add(int(m[np.argmin(hist["E15"][m])]), f"min E15 {c}")
    fi = np.where(ok)[0]
    if len(fi):
        for i in fi[np.argsort(hist["E15"][fi])[:10]]:
            add(int(i), "top10 E15")
        for name, pts in (("rEV2", [(hist["E15"][i], abs(hist["rEV2"][i] - HELD["rEV2"])) for i in fi]),
                          ("rn2", [(hist["E15"][i], abs(hist["rn2"][i] - HELD["rn2"])) for i in fi]),
                          ("rMV2", [(hist["E15"][i], abs(hist["rMV2"][i] - HELD["rMV2"])) for i in fi]),
                          ("proxy", [(hist["E15"][i], -hist["proxy"][i]) for i in fi])):
            front = pareto(pts)
            front = sorted(front, key=lambda j: pts[j][0])[:cap_front]
            for j in front:
                add(int(fi[j]), f"front {name}")
    for r in slp_rows:
        if r.get("ok") and r["endpoint"]["ok"] and r["endpoint"]["feasible"]:
            key = (r["cell"], tuple(np.round(np.array(r["theta"]), 6)))
            sel.setdefault(key, {"theta": r["theta"], "cell": r["cell"], "why": [], "coarse": {k: float(r["endpoint"][k]) for k in ("E15", "merit", "proxy", "gate", "Kgamma", "M", "z0", "rEV2", "rn2", "rMV2", "gA", "Delta")}})
            sel[key]["why"].append(f"slp {r['start']}")
    return list(sel.values())


def _fine_worker(cand):
    theta = np.array(cand["theta"]); cell = cand["cell"]
    t0 = time.time(); out = dict(cand); out["fine_ok"] = False
    try:
        p = S.params_cell(theta, cell)
        o = S.fine_eval(theta, cell)
        if o is None:
            out["fine_fail"] = "continuation returned None"; return out
        x, F, z = o["profile"]["x"], o["profile"]["F"], o["profile"]["z"]
        proxy, fac, gw = S.V48B["omega_finite"](p, x, F)
        ratio, B4 = S.V48B["merged"](p, (x, F, z), o["M"])
        gatev = S.VR.gate_dev(p); gate = float(np.max(np.abs(gatev)))
        dz2 = np.diff(z ** 2) / np.diff(x)
        kg = float(0.5 - S.V48.E2 * p["ell"] ** 2 * np.max(dz2 ** 2) / (4 * np.pi ** 4 * p["c"]))
        fine = {"E15": float(o["E15"]), "errors": [float(v) for v in o["errors"]], "M": float(o["M"]), "z0": float(o["z0"]), "Delta": float(o["Delta"]), "gA": float(o["gA"]),
                "rEV2": float(o["rEV2"]), "rn2": float(o["rn2"]), "rMV2": float(o["rMV2"]), "gate": gate, "Kgamma": kg, "proxy": float(proxy), "omega_factor": float(fac),
                "g2w_4pi_finite": float(gw), "E4_over_4E1": float(ratio), "merged_4He_binding_MeV": float(B4), "E_tr": float(o["E_tr"]), "eta": float(p["eta"]),
                "rES2_consistent_target_error": float(100 * (o["rES2"] / (0.84087 ** 2 - 0.1161) - 1))}
        fine["feasible"] = bool(gate <= S.GATE_ACC and kg >= S.KG_MIN and p["eta"] <= S.ETA_MAX + 1e-12 and proxy >= S.W)
        fine["signs_right"] = bool(fine["rEV2"] > 0 and fine["rn2"] < 0 and fine["rMV2"] > 0)
        fine["gw_gate"] = bool(gw >= GW_GATE)
        out["fine"] = fine; out["fine_ok"] = True
    except Exception as ex:
        out["fine_fail"] = repr(ex)[:200]
    out["fine_wall_s"] = time.time() - t0
    return out


def prices_at(theta, cell, tag):
    v = {k: float(np.exp(t)) for k, t in zip(S.LOGP, theta[:9])}; v["ell"] = float(theta[9])
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        pr = S.V55["prices"](v, dict(S.CELLS[cell]), tag)
    lam = pr["0.01"]["lam"]
    return {"trust_0.01": pr["0.01"], "n_active": pr["n_active"], "rank_active": pr["rank_active"], "lam": lam,
            "Delta_active": bool(lam["Delta"] > ACTIVE), "gA_active": bool(lam["gA"] > ACTIVE), "stdout_tail": buf.getvalue()[-600:]}


def run(workers):
    S.consolidate()
    rec = json.load(open(S.OUT_JSON))
    hist = dict(np.load(os.path.join(S.OUT_DIR, "history.npz"), allow_pickle=False))
    slp_rows = [json.loads(l) for l in open(os.path.join(S.OUT_DIR, "hist_slp_endpoints.jsonl"))] if os.path.exists(os.path.join(S.OUT_DIR, "hist_slp_endpoints.jsonl")) else []
    cands = candidates(hist, slp_rows)
    log(f"{len(cands)} candidate states for the fine grid")
    ctx = mp.get_context("spawn")
    with ctx.Pool(workers) as pool:
        rows = pool.map(_fine_worker, cands)
    fine = [r for r in rows if r["fine_ok"]]
    json.dump(rows, open(os.path.join(S.OUT_DIR, "front_candidates.json"), "w"), indent=1, default=float)
    feas = [r for r in fine if r["fine"]["feasible"]]
    out = {"candidates": len(cands), "fine_ok": len(fine), "fine_feasible": len(feas), "fine_failures": [r.get("fine_fail") for r in rows if not r["fine_ok"]][:10]}
    # the per-cell fine minimum and its prices
    out["cell_minimum"] = {}
    for c in S.CELLS:
        fc = [r for r in feas if r["cell"] == c]
        if not fc:
            out["cell_minimum"][c] = None; continue
        b = min(fc, key=lambda r: r["fine"]["E15"])
        pr = prices_at(np.array(b["theta"]), c, c)
        out["cell_minimum"][c] = {"theta": b["theta"], "params": {k: float(np.exp(t)) for k, t in zip(S.LOGP, b["theta"][:9])} | {"ell": float(b["theta"][9])},
                                  "fine": b["fine"], "coarse_E15": b["coarse"]["E15"], "why": b["why"], "prices": pr}
        log(f"{c}: fine minimum E15 {b['fine']['E15']:.4f} (coarse {b['coarse']['E15']:.4f}), gw finite {b['fine']['g2w_4pi_finite']:.2f}, radii {b['fine']['rEV2']:.3f} {b['fine']['rn2']:.3f} {b['fine']['rMV2']:.3f}, prices Delta {pr['lam']['Delta']:.3f} gA {pr['lam']['gA']:.3f}")
    # the fronts on the fine grid (feasible candidates)
    fronts = {}
    for name, key in (("rEV2", lambda r: abs(r["fine"]["rEV2"] - HELD["rEV2"])), ("rn2", lambda r: abs(r["fine"]["rn2"] - HELD["rn2"])),
                      ("rMV2", lambda r: abs(r["fine"]["rMV2"] - HELD["rMV2"])), ("g2w_4pi_finite", lambda r: -r["fine"]["g2w_4pi_finite"]),
                      ("merged_4He", lambda r: abs(r["fine"]["merged_4He_binding_MeV"] - HE4_BINDING))):
        pts = [(r["fine"]["E15"], key(r)) for r in feas]
        idx = pareto(pts) if pts else []
        fronts[name] = [{"cell": feas[i]["cell"], "E15": feas[i]["fine"]["E15"], "value": feas[i]["fine"][name] if name != "merged_4He" else feas[i]["fine"]["merged_4He_binding_MeV"],
                         "M": feas[i]["fine"]["M"], "z0": feas[i]["fine"]["z0"], "theta": feas[i]["theta"]} for i in sorted(idx, key=lambda i: pts[i][0])]
    out["fronts_fine"] = fronts
    # the words
    falsifiers = [r for r in feas if r["fine"]["gw_gate"] and r["fine"]["E15"] < E15_FALSIFY and r["fine"]["signs_right"]]
    near_misses = [r for r in fine if r["fine"]["E15"] < E15_FALSIFY and not (r["fine"]["feasible"] and r["fine"]["gw_gate"] and r["fine"]["signs_right"])]
    gmin = min(feas, key=lambda r: r["fine"]["E15"]) if feas else None
    word, detail = "NEITHER", ""
    if falsifiers:
        b = min(falsifiers, key=lambda r: r["fine"]["E15"])
        word = "THE CLAIM FALSIFIED"; detail = f"a feasible state at fine E15 {b['fine']['E15']:.3f} with g_omega^2/4pi (finite) {b['fine']['g2w_4pi_finite']:.2f} and the three held-out signs right, cell {b['cell']}"
    elif gmin is not None:
        pr = out["cell_minimum"][gmin["cell"]]["prices"]
        within = abs(gmin["fine"]["E15"] - CONFIRM_CENTER) <= CONFIRM_HALF
        if within and pr["Delta_active"] and pr["gA_active"]:
            word = "THE CLAIM CONFIRMED"; detail = f"the feasible fine-grid minimum {gmin['fine']['E15']:.3f} (cell {gmin['cell']}) lies within {CONFIRM_HALF} of {CONFIRM_CENTER} with Delta ({pr['lam']['Delta']:.3f}) and gA ({pr['lam']['gA']:.3f}) active"
        else:
            detail = f"the feasible fine-grid minimum {gmin['fine']['E15']:.3f} (cell {gmin['cell']}); within 0.5 of 9.2: {within}; Delta active {pr['Delta_active']} ({pr['lam']['Delta']:.3f}), gA active {pr['gA_active']} ({pr['lam']['gA']:.3f})"
    else:
        detail = "no feasible fine-grid state"
    out["word"] = word; out["word_detail"] = detail
    out["global_minimum_fine"] = {"cell": gmin["cell"], "fine": gmin["fine"], "theta": gmin["theta"], "why": gmin["why"]} if gmin else None
    out["falsifiers"] = [{"cell": r["cell"], "E15": r["fine"]["E15"], "gw": r["fine"]["g2w_4pi_finite"], "theta": r["theta"]} for r in falsifiers]
    out["below_8_but_not_counted"] = [{"cell": r["cell"], "E15": r["fine"]["E15"], "feasible": r["fine"]["feasible"], "gw_finite": r["fine"]["g2w_4pi_finite"], "signs_right": r["fine"]["signs_right"],
                                       "radii": [r["fine"]["rEV2"], r["fine"]["rn2"], r["fine"]["rMV2"]], "proxy": r["fine"]["proxy"], "gate": r["fine"]["gate"], "theta": r["theta"]} for r in near_misses]
    out["wording_limit"] = ("no outcome here is a lower bound on E15: the strongest wording is 'no counterexample found over the declared space', "
                            "with the box, the cells, the evaluation count and the multistart spread stated beside it")
    # the V48 state on the fine grid for reference, and the coarse-vs-fine shift over the candidates
    shifts = [r["fine"]["E15"] - r["coarse"]["E15"] for r in fine]
    out["coarse_to_fine_E15_shift"] = {"mean": float(np.mean(shifts)), "min": float(np.min(shifts)), "max": float(np.max(shifts))} if shifts else None
    rec["front"] = out
    json.dump(rec, open(S.OUT_JSON, "w"), indent=1, default=float)
    plot(hist, fine, feas, rec)
    log(f"word: {word}; {detail}")
    return out


def plot(hist, fine, feas, rec):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ok = hist["ok"] & hist["feasible"]; cell = hist["cell"]
    colors = {"n3g2": "C0", "n4g2": "C1", "n3g4": "C2", "n4g4": "C3"}
    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    panels = [("rEV2", "r2_EV (fm2), target 0.823 (axis clipped to 0..3)", (0, 3)), ("rn2", "r_n2 (fm2), target -0.1161 (clipped to -1..0.3)", (-1, 0.3)),
              ("rMV2", "r2_MV (fm2), comparator 0.73 (clipped to 0..1.5)", (0, 1.5)), ("proxy", "g_omega^2/4pi: proxy (history), finite (front); clipped to 0..30", (0, 30))]
    v48 = rec.get("controls", {}).get("b_fine", {}); v48g = {"rEV2": v48.get("rEV2"), "rn2": v48.get("rn2"), "rMV2": v48.get("rMV2"), "proxy": 7.913, "E15": v48.get("E15")}
    for a, (k, lab, xl) in zip(ax.flat[:4], panels):
        a.set_xlim(*xl)
        for c in S.CELLS:
            m = ok & (cell == c)
            if m.any():
                a.plot(hist[k][m], hist["E15"][m], ".", ms=2, alpha=0.25, color=colors[c], label=f"{c} coarse ({int(m.sum())})")
        for r in feas:
            v = r["fine"]["g2w_4pi_finite"] if k == "proxy" else r["fine"][k]
            a.plot([v], [r["fine"]["E15"]], "o", ms=4, mfc="none", color=colors[r["cell"]])
        if k in HELD:
            a.axvline(HELD[k], color="k", lw=0.6, ls="--")
        if k == "proxy":
            a.axvline(S.W, color="gray", lw=0.6, ls=":"); a.axvline(GW_GATE, color="k", lw=0.6, ls="--")
        if v48g.get(k) is not None and v48g.get("E15") is not None:
            a.plot([v48g[k]], [v48g["E15"]], "k*", ms=11, mfc="none", label="shipped V48 state (fine)" if k == "rEV2" else None)
        a.axhline(E15_FALSIFY, color="r", lw=0.6, ls="--"); a.axhline(CONFIRM_CENTER, color="g", lw=0.6, ls=":")
        a.set_xlabel(lab); a.set_ylabel("E15 (percent)"); a.set_ylim(min(7.0, min([r["fine"]["E15"] for r in feas] + [9.0]) - 0.5), 16)
    ax[0, 0].legend(fontsize=7)
    a = ax[1, 1]
    for r in feas:
        a.plot([r["fine"]["merged_4He_binding_MeV"]], [r["fine"]["E15"]], "o", ms=4, color=colors[r["cell"]])
    a.axvline(HE4_BINDING, color="k", lw=0.6, ls="--"); a.axvline(1128.1, color="gray", lw=0.6, ls=":"); a.set_xlabel("merged 4He binding 4 M_N - E4 (MeV); measured 28.3 (dashed), shipped V48 1128 (dotted)"); a.set_ylabel("E15 fine (percent)")
    a = ax[1, 2]
    for c in S.CELLS:
        lf = os.path.join(S.OUT_DIR, f"de_{c}.log")
        if os.path.exists(lf):
            best = [float(l.split("best ")[1].split()[0]) for l in open(lf) if " best " in l]
            a.plot(range(1, len(best) + 1), best, "-", color=colors[c], label=c)
    a.set_xlabel("DE generation"); a.set_ylabel("best merit"); a.legend(fontsize=7); a.set_ylim(7, 12)
    fr = rec.get("front", {})
    fig.suptitle(f"M5.32 R30 (S3): the fine-grid front on the author's instrument; word: {fr.get('word')}", fontsize=10)
    fig.tight_layout(); os.makedirs(os.path.dirname(PLOT), exist_ok=True); fig.savefig(PLOT, dpi=130); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", nargs="?", default="front"); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if a.mode == "plot":                                                  # the figure alone, from the stored candidates and history
        rec = json.load(open(S.OUT_JSON)); hist = dict(np.load(os.path.join(S.OUT_DIR, "history.npz"), allow_pickle=False))
        rows = json.load(open(os.path.join(S.OUT_DIR, "front_candidates.json"))); fine = [r for r in rows if r["fine_ok"]]
        plot(hist, fine, [r for r in fine if r["fine"]["feasible"]], rec); log("plot written"); return
    run(a.workers)


if __name__ == "__main__":
    main()
