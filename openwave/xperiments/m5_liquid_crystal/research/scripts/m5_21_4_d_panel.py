"""M5.21.4 panel: seed-level E_int(d) + Coulomb prediction | the
same-sector string fit | the capture energy ledger | the charge
dissolution. Out: ../plots/m5_21_4_panel.png"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")


def main(ev_tag="cap"):
    sl = json.load(open(os.path.join(DATA, "m5_21_4_seed_ladder.json")))
    ev = json.load(open(os.path.join(DATA,
                                     f"m5_21_4_ev_{ev_tag}_rows.json")))
    post = json.load(open(os.path.join(DATA,
                                       f"m5_21_4_{ev_tag}_post.json")))
    er = ev["rows"]

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    a = ax[0, 0]
    for nk, mk in (("n32", "o"), ("n48", "s")):
        blk = sl[nk]
        ref = blk["E_single"] + blk["E_mirror"]
        rr = sorted([r for r in blk["rows"] if r["kind"] == "anti"],
                    key=lambda r: r["d"])
        ds = np.array([r["d"] for r in rr])
        Ei = np.array([r["E"] for r in rr]) - ref
        a.plot(ds, Ei, mk + "-", label=f"anti E_int ({nk})")
        dd = np.linspace(ds.min(), ds.max(), 50)
        a.plot(dd, -64 * np.pi * blk["c2_seed"] / dd, "--", lw=1,
               label=f"pred −64πc2/d ({nk}, c2={blk['c2_seed']:.3f})")
    a.axhline(0, color="k", lw=0.5)
    a.set_xlabel("d")
    a.set_ylabel("E_int(d)")
    a.set_title("antipair seed-level binding vs the derived Coulomb form")
    a.legend(fontsize=7)

    b = ax[0, 1]
    for nk, mk in (("n32", "o"), ("n48", "s")):
        blk = sl[nk]
        ref = blk["E_single"] + blk["E_mirror"]
        rr = sorted([r for r in blk["rows"] if r["kind"] == "same"],
                    key=lambda r: r["d"])
        ds = np.array([r["d"] for r in rr])
        Ei = np.array([r["E"] for r in rr]) - ref
        sfit = np.polyfit(ds, Ei, 1)
        b.plot(ds, Ei, mk, label=f"same E_int ({nk})")
        b.plot(ds, np.polyval(sfit, ds), "-", lw=1,
               label=f"linear fit slope {sfit[0]:.2f} ({nk})")
    b.set_xlabel("d")
    b.set_ylabel("E_int(d)")
    b.set_title("same-sector: the string term (linear in d)")
    b.legend(fontsize=7)

    c = ax[1, 0]
    t = [r["t"] for r in er]
    c.plot(t, [r["E"] for r in er], label="E")
    c.plot(t, [r["KE"] for r in er], label="KE")
    c.plot(t, [r["absorbed"] for r in er], label="absorbed (by difference)")
    c.set_xlabel("t")
    c.set_title(f"the {ev_tag} run energy ledger (gate: γ=0 cons. 4.6e-6)")
    c.legend(fontsize=8)

    d = ax[1, 1]
    tp = [r["t"] for r in post]
    d.plot(tp, [r["Q_top_fix"] for r in post], "-o", ms=3,
           label="Q(top cube, fixed at +12)")
    d.plot(tp, [r["Q_bot_fix"] for r in post], "-o", ms=3,
           label="Q(bot cube, fixed at −12)")
    d.plot(tp, [r["Q_far"] for r in post], "-", label="Q(far)")
    d.axhline(0, color="k", lw=0.5)
    d.set_xlabel("t")
    d.set_title("charge dissolution (fixed cubes; seed exact ±1, audit)")
    d.legend(fontsize=8)

    fig.suptitle("M5.21.4: two-defect interaction + antipair capture "
                 "(T2 sym verified-L stack)")
    fig.tight_layout()
    out = os.path.join(PLOTS, "m5_21_4_panel.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
