"""M5.21.8 panel: the analytical-hedgehog verification summary.

(a) lattice E_u(m) (the author's family on the audited 4D stack) + fine inset
(b) kin(m): the clock-flow kinetic sign structure
(c) the omega^2-coefficient sign map of the author's own Hm (region map)
(d) the g-ladder: m*_lattice vs the author's (1/2) ln((1+g)/(g-1))

Run:  python3 m5_21_8_d_panel.py
Out:  ../plots/m5_21_8_panel.png
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")


def load(tag):
    return json.load(open(os.path.join(DATA, f"m5_21_8_{tag}.json")))


fig, ax = plt.subplots(2, 2, figsize=(13, 10))

# (a) E_u(m)
c = load("lat_mcurve_s-1_g8_n32")
ms = [r["m"] for r in c["rows"]]
Es = [r["E_u"] for r in c["rows"]]
ax[0, 0].semilogy(ms, Es, "o-", ms=3, label="lattice E_u(m), g=8")
f = load("lat_fine_g8")
ax[0, 0].semilogy([r["m"] for r in f["rows"]],
                  [r["E_u"] for r in f["rows"]], ".-", ms=2,
                  label="fine grid")
f24 = load("lat_fine_g8_n24")
ax[0, 0].semilogy([r["m"] for r in f24["rows"]],
                  [r["E_u"] for r in f24["rows"]], ".--", ms=2,
                  label="fine, n=24 (h-robustness)")
for v, lab in ((0.12566, "the author's m* (+)"), (-0.12566, "the author's m* (-)")):
    ax[0, 0].axvline(v, color="gray", ls=":", lw=1)
ax[0, 0].set_xlabel("boost-hedgehog rapidity m")
ax[0, 0].set_ylabel("E_u (E_V = 0 exactly on the family)")
ax[0, 0].set_title("(a) the author's dressed family on the audited stack: "
                   "twin finite minima near the author's m*")
ax[0, 0].legend(fontsize=8)
ax[0, 0].grid(alpha=0.3)

# (b) kin(m)
ks = [r["kin"] for r in c["rows"]]
ax[0, 1].plot(ms, ks, "o-", ms=3)
ax[0, 1].axhline(0, color="k", lw=1)
ax[0, 1].axvline(c["m_star_lattice"], color="tab:red", ls=":",
                 label=f"E-minimum m* = {c['m_star_lattice']:+.3f}")
ax[0, 1].axvline(-c["m_star_lattice"], color="tab:red", ls=":")
ax[0, 1].set_xlabel("m")
ax[0, 1].set_ylabel("kin (E_kin = omega^2 kin, eta-reading)")
ax[0, 1].set_title("(b) the clock-flow kinetic: POSITIVE at the "
                   "minima (omega* = 0), negative only past "
                   "|m| ~ 0.125")
ax[0, 1].legend(fontsize=8)
ax[0, 1].grid(alpha=0.3)

# (c) region map
p0 = load("p0")
m = p0["C3_region_map"]
sign = np.array(m["sign"])
ax[1, 0].imshow(sign, aspect="auto", cmap="RdBu", vmin=-1, vmax=1,
                origin="lower")
ax[1, 0].set_xticks(range(len(m["dls"])),
                    [f"{d:g}" for d in m["dls"]], rotation=45)
ax[1, 0].set_yticks(range(len(m["gs"])),
                    [f"{g:g}" for g in m["gs"]])
ax[1, 0].set_xlabel("delta")
ax[1, 0].set_ylabel("g")
ax[1, 0].set_title("(c) sign of the author's Hm omega^2 coefficient: "
                   "NEGATIVE (red) for all g > 1, delta < 2\n"
                   "incl. (1e10, 1e-10); positive only outside the "
                   "construction")
# (d) g-ladder
gl = load("lat_gladder")
gs = [r["g"] for r in gl]
fine = {8.0: load("lat_fine_g8")["m_star_lattice"],
        16.0: load("lat_fine_g16")["m_star_lattice"],
        32.0: load("lat_fine_g32")["m_star_lattice"],
        64.0: load("lat_fine_g64")["m_star_lattice"]}
mlat = [abs(fine.get(r["g"], r["m_star_lattice"])) for r in gl]
mhis = [r["m_star_his"] for r in gl]
ax[1, 1].loglog(gs, mhis, "o-", label="the author's m* = 0.5 ln((1+g)/(g-1))")
ax[1, 1].loglog(gs, mlat, "s-", label="lattice m* (fine grids)")
for g, a, b in zip(gs, mlat, mhis):
    ax[1, 1].annotate(f"{a / b:.2f}", (g, a), fontsize=7,
                      textcoords="offset points", xytext=(4, -10))
ax[1, 1].set_xlabel("g")
ax[1, 1].set_ylabel("m*")
ax[1, 1].set_title("(d) the gravitational-mass law: lattice tracks "
                   "the author's 1/g formula at ~0.83x across "
                   "the ladder")
ax[1, 1].legend(fontsize=8)
ax[1, 1].grid(alpha=0.3, which="both")

fig.suptitle("M5.21.8: the analytical twisting massive hedgehog, "
             "verified (m-sector) and refuted-as-stated (omega "
             "free minimization)")
fig.tight_layout()
out = os.path.join(PLOTS, "m5_21_8_panel.png")
fig.savefig(out, dpi=110)
print("wrote", out)
