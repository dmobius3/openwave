"""M5.21.15 panel figure, rebuildable from the saved JSONs.

Six panels: (1) the free coupled envelope E*(omega) per channel with
the concavity certificate; (2) kin_tot(b*(omega)) per channel (sign =
the favorability read); (3) the R-cut sign map at the channel-minimal
kin; (4) the fixed-J money curve E(omega)|_J with the interior
minimum; (5) the E(J) constrained envelope with dE/dJ = omega*;
(6) beta*(r) profiles (free at omega_max vs fixed-J mid).

Regen: python3 m5_21_15_d_panel.py
Out:   ../plots/m5_21_15_panel.png
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")


def main():
    with open(os.path.join(DATA, "m5_21_15_coupled.json")) as f:
        cp = json.load(f)
    with open(os.path.join(DATA, "m5_21_15_fixedj.json")) as f:
        fj = json.load(f)
    blk = cp["main_s-1_g32"]

    fig, ax = plt.subplots(2, 3, figsize=(15.5, 8.5))

    a = ax[0, 0]
    for ch, mk in (("clock", "o"), ("rot_z", "s"), ("boost_z", "^")):
        lad = blk[f"scan_{ch}"]["ladder"]
        om = [r["omega"] for r in lad]
        E = [r["E_env"] for r in lad]
        a.plot(om, E, mk + "-", label=f"{ch} "
               f"({blk[f'scan_{ch}']['monotone']})")
    a.axhline(0, color="k", lw=0.5)
    a.set_yscale("symlog", linthresh=100)
    a.set_xlabel("omega (family units)")
    a.set_ylabel("E* - E_base (corr + omega^2 kin_tot)")
    a.set_title("A2: the FREE coupled envelope (concave in omega^2,\n"
                "no interior minimum: the theorem, measured)")
    a.legend(fontsize=8)

    a = ax[0, 1]
    for ch, mk in (("clock", "o"), ("rot_z", "s"), ("boost_z", "^")):
        lad = blk[f"scan_{ch}"]["ladder"]
        om = [r["omega"] for r in lad]
        kt = [r["kin_tot"] for r in lad]
        a.plot(om, kt, mk + "-", label=ch)
    a.axhline(0, color="k", lw=0.5)
    a.set_yscale("symlog", linthresh=10)
    a.set_xlabel("omega")
    a.set_ylabel("kin_tot(b*(omega))")
    a.set_title("the omega^2 coefficient at the free minimizer")
    a.legend(fontsize=8)

    a = ax[0, 2]
    chans = ("clock", "rot_z", "boost_z")
    rc_labels = None
    for i, ch in enumerate(chans):
        ct = blk[f"minkin_{ch}"]["cut_table"]
        if rc_labels is None:
            rc_labels = list(ct.keys())
        kt = [ct[k]["kin_base"] + ct[k]["kin_corr"]
              for k in rc_labels]
        a.plot(range(len(rc_labels)), kt, "o-", label=f"{ch} "
               f"(min kin_tot = "
               f"{blk[f'minkin_{ch}']['kin_tot_min']:.3g})")
    a.axhline(0, color="k", lw=0.5)
    a.set_xticks(range(len(rc_labels)), rc_labels)
    a.set_yscale("symlog", linthresh=10)
    a.set_ylabel("kin_tot at channel-minimal b")
    a.set_title("A3: the R-cut sign map at min-kin dressing")
    a.legend(fontsize=8)

    a = ax[1, 0]
    with open(os.path.join(DATA, "m5_21_15_fom_narrow.json")) as f:
        fn = json.load(f)
    rows = [r for r in fn["rows"] if r["feasible"]]
    om = [r["omega"] for r in rows]
    E = [r["E_total"] for r in rows]
    a.plot(om, E, "o-", color="tab:red",
           label=f"J = {fn['J']:.1f}, guard {fn['bound']:g}")
    inf_rows = [r for r in fn["rows"] if not r["feasible"]]
    if inf_rows:
        a.plot([r["omega"] for r in inf_rows],
               [r["E_total"] for r in inf_rows], "x", color="gray",
               label="infeasible (kin target out of family reach)")
    if "omega_min" in fn:
        a.axvline(fn["omega_min"], ls=":", color="tab:red")
    und = fn["J"] / (2.0 * fn["kin_base"])
    a.axvline(und, ls="--", color="gray", lw=0.8)
    a.axhline(0, color="k", lw=0.5)
    a.set_xlabel("omega")
    a.set_ylabel("E_total (E_base + corr + omega^2 kin)")
    a.set_title("A4: E(omega) at FIXED J, narrow guard: the interior\n"
                "minimum at positive energy and nonzero omega")
    a.legend(fontsize=8)

    a = ax[1, 1]
    ej = fj["EJ"]
    Js = [r["J"] for r in ej]
    Es = [r["E_total"] for r in ej]
    a.plot(Js, Es, "o-", label="E(J), dressed constrained envelope")
    a2 = a.twinx()
    a2.plot(Js, [r["dEdJ_over_omega_star"] for r in ej], "s--",
            color="tab:green", label="dE/dJ / omega*")
    a2.axhline(1.0, color="tab:green", lw=0.5, ls=":")
    a2.set_ylabel("dE/dJ / omega*", color="tab:green")
    a.set_xlabel("J (family units)")
    a.set_ylabel("E_total")
    a.set_title("the E(J) envelope + clock thermodynamics "
                "dE/dJ = omega*")
    a.legend(fontsize=8, loc="upper left")

    a = ax[1, 2]
    rs = np.geomspace(0.15, 24.0, 72)
    lad = blk["scan_clock"]["ladder"]
    g = blk["g"]
    RHOS = np.geomspace(0.5, 16.0, 9)

    def b_of(avec, r):
        val = avec[0] * np.tanh(r / 2.0)
        for k, rho in enumerate(RHOS):
            val = val + avec[k + 1] * (r / rho) \
                * np.exp(-((r / rho) ** 2))
        return val

    a.plot(rs, g * b_of(np.array(lad[-1]["avec"]), rs),
           label=f"free b* at omega = {lad[-1]['omega']:g} (clock)")
    mid = [r for r in fj["EJ"]
           if r["omega_target_undressed"] == 0.5][0]
    a.plot(rs, g * b_of(np.array(mid["avec"]), rs), "--",
           label=f"fixed-J b* (J = {mid['J']:.1f})")
    a.axhline(0, color="k", lw=0.5)
    a.set_xlabel("r")
    a.set_ylabel("beta = g b(r)")
    a.set_title("the minimizing dressings (free vs fixed-J)")
    a.legend(fontsize=8)

    fig.suptitle(
        "M5.21.15: free minimization has no interior omega-minimum "
        "(concave envelope); the minimum EXISTS at fixed J "
        "(s = -1, g = 32, exact dressed functional)")
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(os.path.join(PLOTS, "m5_21_15_panel.png"), dpi=110)
    print("panel written")


if __name__ == "__main__":
    main()
