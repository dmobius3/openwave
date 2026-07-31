"""M5.22 ranking + reads: the census table, the selection principle, plots.

Consumes ../data/m5_22_census.json (m5_22_b_census.py collect). The
selection principle (author, 2026-07-28 13:20): identity is ASSIGNED by
ranking, per MEASURED charge class of the relaxed state: proton-analog =
the lightest |Q| = 1 baryon-family state, neutron-analog = the lightest
Q = 0 baryon-family state; heavier states = candidate excited baryons;
the E family is the charged-lepton reference (reads (v)-(vi)).

Charge class = round(|q_far|); class confidence = |q_far| - round.
Dedup within a class: pairwise relative field distance (same-n npz) +
relative energy gap; states within (dist < 0.08, dE/E < 0.01) merge.

Modes: rank [conv=plane scale=6 n=32 bc=pinned delta=0.3] | panel
Out: ../data/m5_22_ranking.json, ../plots/m5_22_census_panel.png,
     ../plots/m5_22_profiles.png
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")


def load(sel):
    with open(os.path.join(DATA, "m5_22_census.json")) as f:
        J = json.load(f)
    rows = []
    for k, r in J.items():
        if r.get("qshift", 0.0) or r["tag"].startswith("q38"):
            continue
        if (r["conv"] == sel["conv"] and r["scale"] == sel["scale"]
                and r["n"] == sel["n"] and r["bc"] == sel["bc"]
                and abs(r["delta"] - sel["delta"]) < 1e-12
                and not r["tag"].startswith("shakedown")):
            rows.append(r)
    return rows


def field_dist(tag_a, tag_b):
    try:
        A = np.load(os.path.join(DATA, f"m5_22_end_{tag_a}.npz"))["M"]
        B = np.load(os.path.join(DATA, f"m5_22_end_{tag_b}.npz"))["M"]
    except FileNotFoundError:
        return float("nan")
    num = float(np.sqrt(np.sum((A - B) ** 2)))
    den = float(np.sqrt(np.sqrt(np.sum(A ** 2)) * np.sqrt(np.sum(B ** 2))))
    return num / max(den, 1e-300)


def rank(sel):
    rows = load(sel)
    for r in rows:
        q = r["charge"]["q_far"]
        r["q_abs"] = abs(q)
        r["q_class"] = int(round(abs(q)))
        r["q_conf"] = abs(abs(q) - round(abs(q)))
        # the citation gate: the selection principle operates only on
        # states the instrument certifies: discretization-consistent
        # (I1 bar, 2b note) AND stationary (f_tol, or a residual force
        # below 1e-5 at max_iter). Gated-out states stay listed as
        # not-citable (still descending / under-resolved), never ranked.
        fmax_last = r["trace"][-1]["fmax"] if r.get("trace") else 1e9
        r["citable"] = bool(
            r["consistency"]["xstencil_ratio"] <= 1.5
            and (r["stop"] == "f_tol" or fmax_last <= 1e-5)
            and r["q_conf"] < 0.15)
    classes = {}
    for r in sorted(rows, key=lambda r: r["E_end"]):
        classes.setdefault(r["q_class"], []).append(r)
    # dedup within class
    for qc, rs in classes.items():
        for i, r in enumerate(rs):
            r["dup_of"] = None
            for p in rs[:i]:
                if p["dup_of"]:
                    continue
                dE = abs(r["E_end"] - p["E_end"]) / max(p["E_end"], 1e-300)
                fd = field_dist(r["tag"], p["tag"])
                if dE < 0.01 and fd < 0.08:
                    r["dup_of"] = p["tag"]
                    break
    baryon = {qc: [r for r in rs if r["fam"] in ("N", "P")
                   and r["citable"]]
              for qc, rs in classes.items()}
    lepton = [r for rs in classes.values() for r in rs
              if r["fam"] == "E" and r["citable"]]
    out = {"selection": sel, "n_rows": len(rows),
           "not_citable": [{k: r[k] for k in ("tag", "E_end", "stop")}
                           | {"xr": r["consistency"]["xstencil_ratio"],
                              "q_far": r["charge"]["q_far"]}
                           for r in rows if not r["citable"]]}
    p_cand = next((r for r in baryon.get(1, []) if not r["dup_of"]), None)
    n_cand = next((r for r in baryon.get(0, []) if not r["dup_of"]), None)
    lep = min(lepton, key=lambda r: r["E_end"]) if lepton else None
    out["proton_analog"] = p_cand["tag"] if p_cand else None
    out["neutron_analog"] = n_cand["tag"] if n_cand else None
    out["lepton_ref"] = lep["tag"] if lep else None
    reads = {}
    if p_cand and n_cand:
        reads["ii_mass_ordering"] = {
            "E_p": p_cand["E_end"], "E_n": n_cand["E_end"],
            "neutron_heavier": bool(n_cand["E_end"] > p_cand["E_end"]),
            "ratio": n_cand["E_end"] / p_cand["E_end"]}
    if n_cand:
        prof = n_cand["charge"]["profile"]
        qs = [p[1] for p in prof if np.isfinite(p[1])]
        reads["iii_core_shell"] = {
            "profile": prof,
            "q_inner_peak": (max(qs, key=abs) if qs else None),
            "q_far": n_cand["charge"]["q_far"]}
    if p_cand and lep:
        reads["v_same_charge_hierarchy"] = {
            "E_proton": p_cand["E_end"], "E_lepton": lep["E_end"],
            "proton_heavier": bool(p_cand["E_end"] > lep["E_end"]),
            "ratio": p_cand["E_end"] / max(lep["E_end"], 1e-300)}
    if p_cand and lep:
        reads["vi_q37_topology"] = {
            "q_far_proton": p_cand["charge"]["q_far"],
            "q_far_lepton": lep["charge"]["q_far"],
            "abs_degree_equal": bool(
                abs(abs(p_cand["charge"]["q_far"])
                    - abs(lep["charge"]["q_far"])) < 0.15),
            "ledger_proton": p_cand["ledger"]["th0.05"]["n_components"],
            "ledger_lepton": lep["ledger"]["th0.05"]["n_components"]}
    out["reads"] = reads
    out["ladder"] = {
        str(qc): [{k: r[k] for k in ("tag", "fam", "s", "E_end",
                                     "q_abs", "q_conf", "stop",
                                     "sym_dev_end", "dup_of")}
                  for r in rs]
        for qc, rs in classes.items()}
    with open(os.path.join(DATA, "m5_22_ranking.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: out[k] for k in ("proton_analog",
                                          "neutron_analog",
                                          "lepton_ref")}, indent=1))
    print(json.dumps(out["reads"], indent=1, default=str))
    return out


def panel(sel):
    rows = load(sel)
    with open(os.path.join(DATA, "m5_22_ranking.json")) as f:
        RK = json.load(f)
    fig = plt.figure(figsize=(15, 9))
    # (a) energy ladder per charge class
    ax = fig.add_subplot(2, 2, 1)
    cols = {"E": "tab:green", "N": "tab:blue", "P": "tab:red"}
    for r in rows:
        qc = int(round(abs(r["charge"]["q_far"])))
        x = qc + {"E": -0.18, "N": 0.0, "P": 0.18}[r["fam"]]
        ax.scatter([x], [r["E_end"]], c=cols[r["fam"]], s=40)
        ax.annotate(f"{r['fam']}{r['s']:+g}", (x, r["E_end"]),
                    fontsize=7, xytext=(4, 2),
                    textcoords="offset points")
    ax.set_xlabel("|Q| class (far cube flux, rounded)")
    ax.set_ylabel("E_end")
    ax.set_title("(a) the census ladder: relaxed energy per charge class"
                 "\n(green = lepton family, blue = central-pi, "
                 "red = fractional)")
    # (b) charge profiles of the headline states
    ax = fig.add_subplot(2, 2, 2)
    for key, lab in (("proton_analog", "proton-analog"),
                     ("neutron_analog", "neutron-analog"),
                     ("lepton_ref", "lepton ref")):
        tag = RK.get(key)
        if not tag:
            continue
        r = next(x for x in rows if x["tag"] == tag)
        prof = np.array([[p[0], p[1]] for p in r["charge"]["profile"]
                         if np.isfinite(p[1])])
        ax.plot(prof[:, 0], prof[:, 1], "-o", ms=3, label=f"{lab} ({tag})")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("cube half-width")
    ax.set_ylabel("Q(cube) [Mermin-Ho / 4pi]")
    ax.set_title("(b) charge profile: core/shell structure")
    ax.legend(fontsize=7)
    # (c) energy vs |s| per family
    ax = fig.add_subplot(2, 2, 3)
    for fam in ("N", "P", "E"):
        rs = sorted([r for r in rows if r["fam"] == fam],
                    key=lambda r: r["s"])
        if rs:
            ax.plot([r["s"] for r in rs], [r["E_end"] for r in rs],
                    "-o", color=cols[fam], label=fam)
    ax.set_xlabel("seed sign s")
    ax.set_ylabel("E_end")
    ax.set_title("(c) relaxed energy vs seed sign")
    ax.legend()
    # (d) seed charge vs relaxed charge
    ax = fig.add_subplot(2, 2, 4)
    for r in rows:
        ax.scatter([abs(r["q2d_seed"])], [r["q_abs"] if "q_abs" in r
                                          else abs(r["charge"]["q_far"])],
                   c=cols[r["fam"]], s=40)
    lim = max(3.2, max(abs(r["q2d_seed"]) for r in rows) + 0.3)
    ax.plot([0, lim], [0, lim], "k--", lw=0.6)
    ax.set_xlabel("|q2d| of the seed cross-section (author instrument)")
    ax.set_ylabel("|Q_far| of the relaxed 3D state")
    ax.set_title("(d) does the charge class survive rotation + "
                 "relaxation?")
    fig.tight_layout()
    p = os.path.join(PLOTS, "m5_22_census_panel.png")
    fig.savefig(p, dpi=110)
    print(f"saved {p}")


SEL = {"conv": "plane", "scale": 6.0, "n": 32, "bc": "pinned",
       "delta": 0.3}

if __name__ == "__main__":
    kw = {}
    for a in sys.argv[2:]:
        k, v = a.split("=", 1)
        kw[k] = float(v) if k in ("scale", "delta") else \
            (int(v) if k == "n" else v)
    sel = SEL | kw
    mode = sys.argv[1] if len(sys.argv) > 1 else "rank"
    if mode == "rank":
        rank(sel)
    elif mode == "panel":
        panel(sel)
