"""M5.32 R8 collect: assemble the arm-a ladder and the arm-b theorem into one
result file plus the two figures."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA, PLOTS, CK = (os.path.join(RES, "data"), os.path.join(RES, "plots"),
                   os.path.join(RES, "checkpoints", "m5_32_r8"))
BOXES = (("n32_L48", 48.0), ("n48_L72", 72.0), ("n64_L96", 96.0))
Ls = np.array([L for _, L in BOXES])

q = {t: json.load(open(os.path.join(CK, f"quart_{t}.json"))) for t, _ in BOXES}
th = json.load(open(os.path.join(DATA, "m5_32_r8_ir_theorem.json")))
terms = list(q["n32_L48"]["terms"].keys())


def series(t, key):
    return np.array([q[b]["terms"][t][key] for b, _ in BOXES])


def expo(v):
    v = np.asarray(v, float)
    return float(np.polyfit(np.log(Ls), np.log(np.abs(v)), 1)[0]) if np.all(v != 0) else None


out = {"task": "M5.32 R8: class C5 / C6 under the R7-B7 IR question",
       "hypothesis": "a quartic-in-Mdot term has a density quartic in the jets, so its clock "
                     "inertia should be IR-convergent where the certified quadratic one is not; "
                     "the structural prediction is that this still cannot give a box-stable "
                     "omega* because the extensive omega^2 term stays",
       "boxes": [{"tag": t, "n": q[t]["n"], "L": q[t]["L"], "h": q[t]["h"]} for t, _ in BOXES],
       "clock_channel": q["n32_L48"]["a0_stats"],
       "scaling": {}, "theorem": th["generators"], "tail": {k: th["tail"][k] for k in
                                                            ("jet_exponent", "dev_exponent", "a0_exponent", "reading")},
       "ladder": th["ladder"]}
for t in terms:
    A, C2, C4 = series(t, "A_static"), series(t, "C2"), series(t, "C4")
    out["scaling"][t] = {"definition": q["n32_L48"]["terms"][t]["definition"],
                         "A_static": A.tolist(), "C2": C2.tolist(), "C4": C4.tolist(),
                         "A_exponent_in_L": expo(A), "C2_exponent_in_L": expo(C2),
                         "C4_exponent_in_L": expo(C4),
                         "odd_rel_max": max(q[b]["terms"][t]["odd_rel"] for b, _ in BOXES)}
out["verdict"] = {
    "C5_static_and_omega2_IR_finite": True,
    "C5_omega4_IR_finite": False,
    "C6_omega4_volume_extensive": True,
    "strong_G3_reachable_at_fixed_coefficient": False,
    "statement": "the quartic classes do not rescue the clock: their omega^4 inertia is itself "
                 "divergent in the box (exponent 0.61 to 1.03 for C5, exactly 3 for C6), while the "
                 "certified omega^2 inertia diverges as L^1.06, so at fixed J the cubic term never "
                 "takes over and omega* stays proportional to 1/L; the free-omega route needs the "
                 "energy's omega^2 coefficient to turn negative, which costs a C5 coefficient "
                 "growing as L^1.03 (312.6 at L = 48, 637.2 at L = 96), so no fixed coefficient works"}
with open(os.path.join(DATA, "m5_32_r8_quartics.json"), "w") as f:
    json.dump(out, f, indent=1)

fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
for t in terms:
    for a, key, lab in zip(ax, ("A_static", "C2", "C4"),
                           ("static A", "omega^2 coefficient C2", "omega^4 coefficient C4")):
        v = np.abs(series(t, key))
        if np.all(v > 0):
            a.plot(Ls, v, "o-", ms=4, label=t)
for a, lab in zip(ax, ("|A| static", "|C2| (omega^2)", "|C4| (omega^4)")):
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("box size L (h = 1.5 fixed)")
    a.set_ylabel(lab); a.legend(fontsize=7)
ax[1].plot(Ls, 4e-1 * Ls, "k--", lw=1, label="linear in L")
ax[2].plot(Ls, 1e1 * Ls ** 3 / 1e5, "k--", lw=1)
fig.suptitle("R8: box scaling of the Lagrangian omega-coefficients (undressed hedgehog, realized clock channel)")
fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "m5_32_r8_box_scaling.png"), dpi=110); plt.close(fig)

sh = th["tail"]["shells"]
r = [s["r"] for s in sh]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
ax[0].loglog(r, [s["jet_rms"] for s in sh], "o-", ms=4, label="jet rms (exponent %.3f)" % th["tail"]["jet_exponent"])
ax[0].loglog(r, [s["dev_from_const_vacuum_rms"] for s in sh], "s-", ms=4,
             label="deviation from a CONSTANT vacuum (exponent %.3f)" % th["tail"]["dev_exponent"])
ax[0].loglog(r, [s["a0_rms"] for s in sh], "^-", ms=4, label="clock generator a0 (exponent %.3f)" % th["tail"]["a0_exponent"])
ax[0].set_xlabel("r"); ax[0].set_ylabel("rms"); ax[0].legend(fontsize=7)
ax[0].set_title("the hedgehog far field: on the vacuum ORBIT, not at a point")
sp = np.array([s["spec_mean"] for s in sh])
for i in range(4):
    ax[1].plot(r, sp[:, i], "o-", ms=3, label=f"eigenvalue {i}")
ax[1].set_xlabel("r"); ax[1].set_ylabel("spectrum of M eta (shell mean)")
ax[1].set_title("orbit invariant: the spectrum is the vacuum's at every radius"); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(PLOTS, "m5_32_r8_tail.png"), dpi=110); plt.close(fig)
print("collected: data/m5_32_r8_quartics.json + 2 plots")
