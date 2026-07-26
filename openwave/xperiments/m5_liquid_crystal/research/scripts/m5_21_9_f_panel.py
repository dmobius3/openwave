"""M5.21.9 analysis + summary panel.

(a) P0: E(m) + kin(m) of the delta = -0.3 dressed family
(b) P2/P4: the fixed-J family: E_tot(J) three rungs + the Legendre
    check dE/dJ vs omega* (finite differences between rungs)
(c) P3: the Larmor read: phi_xy(t) per eps rung + Omega vs B_meas

Run:  python3 m5_21_9_f_panel.py
Out:  ../plots/m5_21_9_panel.png + printed fit tables
      (rerunnable as rungs land; skips missing files)
"""
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")


def jload(name):
    p = os.path.join(DATA, name)
    return json.load(open(p)) if os.path.exists(p) else None


fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.0))

# (a) P0 m-curve
mc = jload("m5_21_9_lat_mcurve_dl-0.3_g8_n32.json")
if mc:
    ms = [r["m"] for r in mc["rows"]]
    Es = [r["E_u"] for r in mc["rows"]]
    ks = [r["kin"] for r in mc["rows"]]
    ax[0].plot(ms, Es, "o-", ms=3, color="tab:blue",
               label="E_u(m), delta = -0.3")
    ax[0].axhline(0, color="k", lw=0.8)
    for v in (0.1508, -0.1508):
        ax[0].axvline(v, color="gray", ls=":", lw=1)
    ax[0].set_ylim(-40, 400)
    axk = ax[0].twinx()
    axk.plot(ms, ks, "s--", ms=2, color="tab:orange", alpha=0.7,
             label="kin(m)")
    axk.axhline(0, color="tab:orange", lw=0.6, alpha=0.5)
    axk.set_ylabel("kin (clock inertia)", color="tab:orange")
    ax[0].set_xlabel("boost rapidity m")
    ax[0].set_ylabel("E_u", color="tab:blue")
    ax[0].set_title("(a) the delta < 0 dressed family: twin minima\n"
                    "at |m| ~ 0.15 (analytic 0.1508, gray), kin > 0\n"
                    "there; E_u(m*) = -6.5 SUB-VACUUM (flag)")
    ax[0].legend(loc="upper center", fontsize=8)

# (b) the fixed-J family + Legendre
rungs = []
for f in sorted(glob.glob(os.path.join(DATA,
                                       "m5_21_9_fixedj_conj_om*.json"))):
    r = json.load(open(f))
    if "final" in r:
        J = r["start"]["J"]
        om = r["final"]["omega_star_final"]
        E = r["final"]["E_u"] + r["final"]["E_v"] \
            + J * J / (4 * r["final"]["kin_final"])
        rungs.append({"J": J, "om": om, "E": E})
rungs.sort(key=lambda x: x["J"])
if len(rungs) >= 2:
    Js = np.array([r["J"] for r in rungs])
    Es = np.array([r["E"] for r in rungs])
    oms = np.array([r["om"] for r in rungs])
    ax[1].plot(Js, Es, "o-", color="tab:blue")
    for r in rungs:
        ax[1].annotate(f"omega*={r['om']:.3f}", (r["J"], r["E"]),
                       fontsize=8, textcoords="offset points",
                       xytext=(6, -4))
    ax[1].set_xlabel("J (constraint-carried)")
    ax[1].set_ylabel("E_tot(J) = E_stat + J^2/(4 kin)")
    rows = []
    for i in range(len(rungs) - 1):
        dEdJ = (Es[i + 1] - Es[i]) / (Js[i + 1] - Js[i])
        om_mid = 0.5 * (oms[i + 1] + oms[i])
        rows.append({"dEdJ": float(dEdJ), "om_mid": float(om_mid),
                     "ratio": float(dEdJ / om_mid)})
        print("Legendre", rows[-1], flush=True)
    ax[1].set_title("(b) the fixed-J family: E(J) with the\n"
                    "Legendre check dE/dJ = omega* "
                    + ", ".join(f"{r['ratio']:.3f}x" for r in rows))
else:
    ax[1].set_title("(b) fixed-J rungs pending")

# (c) Larmor
lad = []
for pat in ("m5_21_9_larmor_eps*.json", "m5_21_9_larmor_env*.json"):
    for f in sorted(glob.glob(os.path.join(DATA, pat))):
        r = json.load(open(f))
        r["_tag"] = os.path.basename(f).replace(
            "m5_21_9_larmor_", "").replace(".json", "")
        lad.append(r)
BROKEN = {"eps0.03", "eps0.1", "env0.1"}  # arena blew up mid-run:
# the finite prefix is pre-blowup transient, not a precession read
fit_rows = []
for r in lad:
    if r["_tag"] in BROKEN:
        print("larmor", {"tag": r["_tag"],
                         "skipped": "arena blowup (kept as the "
                                    "saturation record only)"},
              flush=True)
        continue
    t = np.array([x["t"] for x in r["rows"]])
    phi = np.array([x["phi_xy"] for x in r["rows"]])
    ok = np.isfinite(phi)
    if ok.sum() < 4:
        print("larmor", {"tag": r["_tag"], "skipped": "non-finite"},
              flush=True)
        continue
    phi = np.unwrap(phi[ok])
    tt = t[ok]
    om_fit = float(np.polyfit(tt, phi, 1)[0])
    fit_rows.append({"tag": r["_tag"], "eps": r["eps"],
                     "B": r["B_meas"], "Omega": om_fit})
    ax[2].plot(tt, phi - phi[0], "-",
               label=f"{r['_tag']} (B={r['B_meas']:.2e})")
    print("larmor", fit_rows[-1], flush=True)
if fit_rows:
    ax[2].set_xlabel("t")
    ax[2].set_ylabel("phi_xy(t) - phi_xy(0)")
    Bs = np.array([x["B"] for x in fit_rows])
    Oms = np.array([x["Omega"] for x in fit_rows])
    if (Bs > 0).sum() >= 2:
        sl = np.polyfit(Bs[Bs >= 0], Oms[Bs >= 0], 1)
        ax[2].set_title("(c) the Larmor read: J-azimuth vs t;\n"
                        f"Omega vs B slope {sl[0]:+.3e} "
                        f"(intercept {sl[1]:+.2e})")
        print("Omega_vs_B", {"slope": float(sl[0]),
                             "intercept": float(sl[1])}, flush=True)
    else:
        ax[2].set_title("(c) Larmor rungs pending")
    ax[2].legend(fontsize=8)
else:
    ax[2].set_title("(c) Larmor rungs pending")

fig.suptitle("M5.21.9: the fixed-J isorotation electron: the "
             "delta < 0 rung, the constraint-carried clock, the "
             "Larmor acceptance")
fig.tight_layout()
out = os.path.join(PLOTS, "m5_21_9_panel.png")
fig.savefig(out, dpi=110)
print("wrote", out)
