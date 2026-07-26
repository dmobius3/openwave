"""M5.21.3 summary panel: saddle spectrum, kin sign table, the
omega-ladders, P1 descent traces.

Run after `m5_21_3_a_4d.py collect`:  python3 m5_21_3_d_panel.py
Out: ../plots/m5_21_3_panel.png
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

with open(os.path.join(DATA, "m5_21_3_all.json")) as f:
    R = json.load(f)

fig, ax = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle("M5.21.3: the 4D lift (N=32, sym stencil, V4 trace-"
             "target, both g branches)", fontsize=12)

# (0,0) saddle spectrum
a = ax[0, 0]
for s, c in (("+1", "tab:blue"), ("-1", "tab:red")):
    h = R.get(f"hess_p1_s{s}")
    if h:
        a.plot(sorted(h["curv_timemix"]), "o-", color=c, ms=4,
               label=f"time-mixing, s={s}")
        a.plot(sorted(h["curv_spatial_ctrl"]), "s--", color=c, ms=4,
               alpha=0.5, label=f"spatial ctrl, s={s}")
a.axhline(0, color="k", lw=0.8)
a.set_title("directional curvatures at the static endpoint\n"
            "(ALL 24 time-mixing directions NEGATIVE = saddle)")
a.set_xlabel("direction (sorted)"); a.set_ylabel("d2E per unit norm2")
a.legend(fontsize=7)

# (0,1) kin sign table
a = ax[0, 1]
p2 = R.get("p2_p1_s+1")
if p2:
    names = list(p2["rows"].keys())
    x = np.arange(len(names))
    for dx, (s, c) in ((-.18, ("+1", "tab:blue")),
                       (.18, ("-1", "tab:red"))):
        rr = R.get(f"p2_p1_s{s}")
        if rr:
            vals = [rr["rows"][nm]["kin"] for nm in names]
            a.bar(x + dx, vals, width=0.34, color=c, label=f"s={s}")
    a.axhline(0, color="k", lw=0.8)
    a.set_xticks(x, names, rotation=25, fontsize=8)
    a.set_title("kin(a0): E_kin = omega^2 * kin\n"
                "(rotations ALL positive; boosts negative)")
    a.legend(fontsize=8)

# (1,0) omega ladders
a = ax[1, 0]
for s, c in (("+1", "tab:blue"), ("-1", "tab:red")):
    r = R.get(f"p3_p1_s{s}_boost_z")
    if r:
        om = [d["omega"] for d in r["ladder"]]
        E = [d["E"] for d in r["ladder"]]
        a.plot(om, E, "o-", color=c, label=f"E*(omega), s={s}")
        st = [d["stop"] for d in r["ladder"]]
        for o, e, t in zip(om, E, st):
            if t in ("dive", "non-finite"):
                a.plot([o], [e], "x", color="k", ms=12, mew=3)
a.set_title("the omega-ladder, boost_z (re-relaxed per rung;\n"
            "x = dive-detector fired)")
a.set_xlabel("omega"); a.set_ylabel("E* = min_M [E_stat + w^2 kin]")
a.legend(fontsize=8)

# (1,1) P1 descent traces
a = ax[1, 1]
for s, c in (("+1", "tab:blue"), ("-1", "tab:red")):
    r = R.get(f"p1_s{s}")
    if r:
        tr = r["trace"]
        a.plot([t["it"] for t in tr], [t["E"] for t in tr], "o-",
               color=c, label=f"s={s} (tail)")
a.set_title("P1 static descent (trace tail; offblock stays 0.0)")
a.set_xlabel("iteration"); a.set_ylabel("E")
a.legend(fontsize=8)

fig.tight_layout()
out = os.path.join(PLOTS, "m5_21_3_panel.png")
fig.savefig(out, dpi=130)
print("wrote", out)
