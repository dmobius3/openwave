"""M5.21.1 phase P4: the (g, delta) scaling ladder toward his physical
regime (Duda 2026-07-15: delta ~ 1e-10, g ~ 1e10, signs open; paper
anchors delta^2 ~ hbar c, g^4 ~ 1e38; "need practical approximations").

Physical values are unreachable directly on any grid: this phase
charts the TOY-REGIME scaling laws of the P1-P3 observables so every
claim carries its regime qualifier + extrapolation direction. One
knob at a time (blindspot 6): the delta-ladder at g = 8, the g-ladder
at delta = 0.3, plus two corner spot checks; sign s = +1 throughout
(P0: block-diagonal statics is sign-mirror-equivalent).

OBSERVABLES per (g, delta) point (64 x 128 relax, 4000 iters)
    E_total; core spectrum (a measured vs (1 + delta)/3, spread);
    containment r50 / r90; melt r_half; axis transverse split
    (C1 read); the k = 0 twist gap omega_true(0) (the KG-fit m at
    k = 0: 1 complex-step HV, the clock-gap proxy); the vacuum
    ladder max omega (analytic, free).

FIT  log-log slopes of each observable along each ladder + the
     honest "does not extrapolate" flag when a fit fails R2 > 0.98.

FILM the delta = 0.5 relax sequence renders in BOTH templates (the
     strongest-biax representative; the film standard for this
     script's evolving sequences).

Run:  python m5_21_1_d_scaling.py [iters]
Out:  ../data/m5_21_1_d_scaling.json,
      ../plots/m5_21_1_d_scaling.png / _film_basic.png /
      _film_thermal.png
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ARGV = sys.argv[1:]
sys.argv = sys.argv[:1]

from m5_17_energy import cell_weights, grid_coords                  # noqa: E402
from m5_20_2_a_eom import (G_T, H, WSCALE, grad_static_4,           # noqa: E402
                           u_eta_density, v4_density)
from m5_20_3_a_constraint import e_static_c, t_total_c              # noqa: E402
from m5_21_a_snap import eig_fields, spectral_amplitude             # noqa: E402
from m5_21_c_clockrun import local_gen                              # noqa: E402
from m5_21_1_c_4d import (FREE42, NR2, NZ2, PIN2, RHO42, W42,       # noqa: E402
                          core64, fire64, hv_cs2, q_meridional64,
                          seed64)
from m5_film import film_strip                                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")

DELTA_LADDER = (0.1, 0.2, 0.3, 0.4, 0.5)     # at g = 8
G_LADDER = (4.0, 8.0, 16.0, 32.0)            # at delta = 0.3
CORNERS = ((16.0, 0.1), (4.0, 0.5))          # spot checks
FILM_POINT = (8.0, 0.5)


def eig_sector_max_omega(delta, g):
    lam = np.array([g, 1.0, delta, 0.0])
    sgn = np.array([-1.0, 1.0, 1.0, 1.0])
    J = np.array([[p * lam[i] ** (p - 1) * sgn[i] for i in range(4)]
                  for p in range(1, 5)])
    ev = np.linalg.eigvalsh(2.0 * WSCALE * (J.T @ J))
    return float(np.sqrt(max(ev.max(), 0.0)))


def obs_point(M, delta, g):
    R, Z = grid_coords(NR2, NZ2, H)
    r = np.sqrt(R ** 2 + Z ** 2)
    w = cell_weights(NR2, NZ2, H)
    dens = (u_eta_density(M, H) + v4_density(M, WSCALE, delta, g)) * w
    Rc, Zc = R[: NR2 - 1, 1:-1], Z[: NR2 - 1, 1:-1]
    rc = np.sqrt(Rc ** 2 + Zc ** 2)
    E_tot = float(dens.sum())
    rb = np.arange(1.0, 60.0, 1.0)
    csum = np.array([dens[rc <= x].sum() for x in rb])
    def r_at(frac):
        i = int(np.searchsorted(csum, frac * csum[-1]))
        return float(rb[min(i, len(rb) - 1)])
    A = spectral_amplitude(M, (g, 1.0, delta, 0.0))
    prof = np.array([np.median(A[(r >= x - 0.5) & (r < x + 0.5)])
                     for x in rb])
    below = np.where(prof < 0.5 * prof[0])[0]
    r_half = float(rb[below[0]]) if below.size else float("nan")
    # axis transverse split (C1)
    S = M[:2, :, 1:4, 1:4]
    lam, V = np.linalg.eigh(0.5 * (S + np.swapaxes(S, -1, -2)))
    zdot = np.abs(V[..., 2, :])
    axk = np.argmax(zdot, axis=-1)
    idx = np.arange(3)[None, None, :]
    lam_t = lam[idx != axk[..., None]].reshape(2, NZ2, 2)
    dax = np.abs(lam_t[..., 1] - lam_t[..., 0]).mean(axis=0)
    z = (np.arange(NZ2) - NZ2 / 2 + 0.5) * H
    mid = (np.abs(z) > 5.0) & (np.abs(z) < 50.0)
    # k = 0 twist gap under the true kinetic metric
    _, Vv = eig_fields(M)
    W = local_gen(Vv[..., :, 0])
    env = np.exp(-((r / 20.0) ** 4))[..., None, None]
    v = env * (W @ M - M @ W) * FREE42
    T2 = float(t_total_c(M, v).real)
    hv = hv_cs2(M, v, g, W42, RHO42, FREE42)
    D2E = float(np.sum(v * hv))
    core = core64(M)
    return {"E_total": E_tot, "r50": r_at(0.5), "r90": r_at(0.9),
            "r_half": r_half, "q": q_meridional64(M),
            "core_a": core["core_mean"],
            "core_spread": core["core_spread"],
            "a_pred": (1.0 + delta) / 3.0,
            "axis_split_mid": float(np.median(dax[mid])),
            "w_twist0_true": float(np.sqrt(D2E / (2 * T2)))
            if T2 > 0 and D2E > 0 else None,
            "T2_twist_positive": bool(T2 > 0),
            "vac_ladder_max_omega": eig_sector_max_omega(delta, g)}


def run_point(g, delta, iters, film=False):
    M0 = seed64(delta, g)
    snaps = (0, 500, 1500, 3000, iters) if film else ()
    # the stiff g-mode tightens the stable FIRE step (the first run
    # went non-finite at g = 32 with the flat dt0 = 0.02)
    dt0 = 0.02 * min(1.0, 8.0 / g)
    M, states, rel = fire64(
        M0, lambda MM: float(e_static_c(MM, WSCALE, delta, g=g).real),
        lambda MM: grad_static_4(MM, WSCALE, delta, g=g, w4=W42,
                                 rho4=RHO42),
        max_iter=iters, dt0=dt0, dt_max=10.0 * dt0, log_every=2000,
        snaps=snaps, tag=f"g{g}_d{delta}")
    if not np.isfinite(M).all():
        row = {"g": g, "delta": delta, "diverged": True,
               "fire_tail": rel["trace"][-1]}
        print(f"[P4 g={g} d={delta}] DIVERGED (non-finite state)",
              flush=True)
        return row
    row = {"g": g, "delta": delta, **obs_point(M, delta, g),
           "fire_tail": rel["trace"][-1]}
    print(f"[P4 g={g} d={delta}] E {row['E_total']:.4f} "
          f"a {row['core_a']:.4f} (pred {row['a_pred']:.4f}) "
          f"r50/r90 {row['r50']:.0f}/{row['r90']:.0f} r_half "
          f"{row['r_half']} axis {row['axis_split_mid']:.4f} "
          f"w_twist0 {row['w_twist0_true']}", flush=True)
    if film:
        for st in states:
            st["vac_spec"] = (g, 1.0, delta, 0.0)
        film_strip(states,
                   os.path.join(PLOTS, "m5_21_1_d_film_basic.png"),
                   template="basic", delta=delta, h=H, g=g,
                   wscale=WSCALE,
                   suptitle=f"M5.21.1-P4 relax at (g, delta) = "
                            f"({g}, {delta})")
        film_strip(states,
                   os.path.join(PLOTS, "m5_21_1_d_film_thermal.png"),
                   template="thermal", delta=delta, h=H, g=g,
                   suptitle=f"M5.21.1-P4 relax at (g, delta) = "
                            f"({g}, {delta}) (thermal)", step=6)
    return row


def fit_ladder(rows, xkey, ykey):
    x = np.array([r.get(xkey, np.nan) for r in rows], dtype=float)
    y = np.array([r.get(ykey) if r.get(ykey) is not None else np.nan
                  for r in rows], dtype=float)
    ok = np.isfinite(y) & (y > 0) & np.isfinite(x) & (x > 0)
    if ok.sum() < 3:
        return {"slope": None, "R2": None,
                "extrapolates": False, "n": int(ok.sum())}
    lx, ly = np.log(x[ok]), np.log(y[ok])
    A = np.vstack([np.ones(ok.sum()), lx]).T
    coef, *_ = np.linalg.lstsq(A, ly, rcond=None)
    pred = A @ coef
    r2 = 1.0 - (np.sum((ly - pred) ** 2)
                / max(np.sum((ly - ly.mean()) ** 2), 1e-300))
    return {"slope": float(coef[1]), "R2": float(r2),
            "extrapolates": bool(r2 > 0.98), "n": int(ok.sum())}


OBS_KEYS = ("E_total", "r50", "r90", "r_half", "core_spread",
            "axis_split_mid", "w_twist0_true", "vac_ladder_max_omega")


def main(iters=4000):
    os.makedirs(DATA, exist_ok=True)
    res = {"task": "M5.21.1", "phase": "P4", "grid": [NR2, NZ2, H],
           "iters": iters, "sign": "+1",
           "physical_regime_note":
           "delta ~ 1e-10, g ~ 1e10 (his 2026-07-15 scales; paper "
           "anchors delta^2 ~ hbar c, g^4 ~ 1e38); toy regime here, "
           "extrapolation by the fitted slopes only where "
           "extrapolates = true"}
    rows_d, rows_g = [], []
    for d in DELTA_LADDER:
        film = (8.0, d) == FILM_POINT
        rows_d.append(run_point(8.0, d, iters, film=film))
    for g in G_LADDER:
        if g == 8.0:
            rows_g.append(rows_d[DELTA_LADDER.index(0.3)])
        else:
            rows_g.append(run_point(g, 0.3, iters))
    corners = [run_point(g, d, iters) for (g, d) in CORNERS]
    res["delta_ladder_g8"] = rows_d
    res["g_ladder_d0p3"] = rows_g
    res["corners"] = corners
    res["fits"] = {
        "vs_delta": {k: fit_ladder(rows_d, "delta", k)
                     for k in OBS_KEYS},
        "vs_g": {k: fit_ladder(rows_g, "g", k) for k in OBS_KEYS}}
    for lad, ff in res["fits"].items():
        for k, f in ff.items():
            print(f"[P4 fit {lad} {k}] slope {f['slope']} R2 {f['R2']}"
                  f" extrapolates={f['extrapolates']}", flush=True)
    # figure
    fig, axes = plt.subplots(2, 3, figsize=(15, 7.5))
    ds = [r["delta"] for r in rows_d]
    gs = [r["g"] for r in rows_g]
    axes[0, 0].plot(ds, [r["core_a"] for r in rows_d], "o-",
                    label="a measured")
    axes[0, 0].plot(ds, [r["a_pred"] for r in rows_d], "--",
                    label="(1+delta)/3")
    axes[0, 0].set_xlabel("delta")
    axes[0, 0].set_title("core mean a(delta)", fontsize=9)
    axes[0, 0].legend(fontsize=7)
    axes[0, 1].plot(ds, [r["core_spread"] for r in rows_d], "o-",
                    label="spread(delta)")
    axes[0, 1].plot(ds, [r["axis_split_mid"] for r in rows_d], "s-",
                    label="axis split (C1)")
    axes[0, 1].plot(ds, ds, ":", lw=0.8, label="delta (combed ref)")
    axes[0, 1].set_xlabel("delta")
    axes[0, 1].legend(fontsize=7)
    axes[0, 2].plot(ds, [r["r50"] for r in rows_d], "o-", label="r50")
    axes[0, 2].plot(ds, [r["r90"] for r in rows_d], "s-", label="r90")
    axes[0, 2].plot(ds, [r["r_half"] for r in rows_d], "v-",
                    label="r_half")
    axes[0, 2].set_xlabel("delta")
    axes[0, 2].legend(fontsize=7)
    for ax, ykey, ttl in ((axes[1, 0], "E_total", "E_total(g)"),
                          (axes[1, 1], "w_twist0_true",
                           "twist gap omega(g)"),
                          (axes[1, 2], "r50", "r50 / r90 (g)")):
        y = [r.get(ykey) if r.get(ykey) is not None else np.nan
             for r in rows_g]
        ax.loglog(gs, y, "o-", label=ykey)
        if ykey == "r50":
            ax.loglog(gs, [r.get("r90", np.nan) for r in rows_g], "s-",
                      label="r90")
        ax.set_xlabel("g")
        ax.set_title(ttl, fontsize=9)
        ax.legend(fontsize=7)
    fig.suptitle("M5.21.1-P4 (g, delta) scaling ladders "
                 "(64x128, sign +1)", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_21_1_d_scaling.png"), dpi=130)
    plt.close(fig)
    with open(os.path.join(DATA, "m5_21_1_d_scaling.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    print("wrote data/m5_21_1_d_scaling.json "
          "+ plots/m5_21_1_d_scaling.png")
    return res


if __name__ == "__main__":
    main(int(ARGV[0]) if ARGV else 4000)
