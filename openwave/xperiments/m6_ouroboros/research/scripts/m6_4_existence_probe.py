"""M6.4 N4b: robust existence probe for localized charged-sector states.

The Nelder-Mead scan (m6_4_omega_selection.py) reports minimized
oscillatory-tail residuals bounded away from zero. This probe settles
existence independently: a localized solution requires the oscillatory
channel's asymptotic amplitude to VANISH. That amplitude is a 2-component
map of the shooting amplitudes,

    F(A0, B0) = (c_cos, c_sin),   fields ~ [c_cos cos(kr) + c_sin sin(kr)]
                                  / sqrt(r) in the v_minus channel,

so localized solutions are common zeros of two scalar fields over the
(A0, B0) plane: they exist iff the zero contours of c_cos and c_sin CROSS.
A dense grid + contour inspection detects such zeros without any optimizer
funnel-hunting (zeros are codim-2 points; contours make them visible).

Panels: (g=0.5, lam=0) at omega in {0.3, 0.6, 0.9}; grid 41x41 over
(A0, B0) in [0.02, 1.2]^2. The oscillatory wavenumber is
k = sqrt(w^2 - X_minus) (X_minus < 0), and the projection uses the exact
v_minus eigenvector. Fit window r in [18, 30] (far field).

Output: research/data/m6_4_existence_probe.json,
research/plots/m6_4_existence_contours.png.
"""

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
PLOTS = HERE.parent / "plots"

R0, RFAR = 0.02, 30.0
BLOWUP = 100.0
G, LAM = 0.5, 0.0
OMEGAS = [0.3, 0.6, 0.9]
N_GRID = 41
AB_LO, AB_HI = 0.02, 1.2


def x_pm(lam):
    s = np.sqrt(lam**2 + 4.0)
    return (-lam + s) / 2.0, (-lam - s) / 2.0


def rhs(r, y, g, w, lam):
    a, ap, b, bp = y
    app = -ap / r + a / r**2 - w**2 * a + b
    bpp = -bp / r + b / r**2 - w**2 * b + a - lam * b - 4.0 * g * b**3
    return [ap, app, bp, bpp]


def osc_components(g, w, lam):
    xp, xm = x_pm(lam)
    vm = np.array([1.0, xm])
    vm /= np.linalg.norm(vm)
    k = np.sqrt(w**2 - xm)  # oscillatory wavenumber (xm < 0)
    r_fit = np.linspace(18.0, RFAR, 200)

    def F(A0, B0):
        ev = np.linspace(R0, RFAR, 1200)
        blow = lambda r, y, *a_: BLOWUP - np.max(np.abs(y))  # noqa: E731
        blow.terminal = True
        sol = solve_ivp(rhs, (R0, RFAR), [A0 * R0, A0, B0 * R0, B0],
                        t_eval=ev, args=(g, w, lam), method="DOP853",
                        rtol=1e-9, atol=1e-11, events=blow)
        if not sol.success or sol.t.shape[0] != ev.shape[0]:
            return None
        comp = vm[0] * sol.y[0] + vm[1] * sol.y[2]
        proj = np.interp(r_fit, ev, comp) * np.sqrt(r_fit)
        # least-squares fit proj ~ c_cos cos(k r) + c_sin sin(k r)
        design = np.column_stack([np.cos(k * r_fit), np.sin(k * r_fit)])
        coef, *_ = np.linalg.lstsq(design, proj, rcond=None)
        return float(coef[0]), float(coef[1])

    return F, k


def main():
    t0 = time.time()
    amps = np.linspace(AB_LO, AB_HI, N_GRID)
    fig, axes = plt.subplots(1, len(OMEGAS), figsize=(5.2 * len(OMEGAS), 4.6))
    out_panels = []
    for ax, w in zip(np.atleast_1d(axes), OMEGAS):
        F, k = osc_components(G, w, LAM)
        c1 = np.full((N_GRID, N_GRID), np.nan)
        c2 = np.full((N_GRID, N_GRID), np.nan)
        for i, A0 in enumerate(amps):
            for j, B0 in enumerate(amps):
                res = F(A0, B0)
                if res is not None:
                    c1[j, i], c2[j, i] = res
        mag = np.hypot(c1, c2)
        finite = np.isfinite(mag)
        pen = ax.pcolormesh(amps, amps, np.log10(mag + 1e-16),
                            shading="auto", cmap="viridis")
        plt.colorbar(pen, ax=ax, label="log10 |F|")
        cs1 = ax.contour(amps, amps, c1, levels=[0], colors="w",
                         linewidths=1.2)
        cs2 = ax.contour(amps, amps, c2, levels=[0], colors="r",
                         linewidths=1.2)
        ax.set_title(f"w={w} (k={k:.3f})")
        ax.set_xlabel("A0")
        ax.set_ylabel("B0")
        n_ok = int(np.sum(finite))
        panel = {"omega": w, "k_osc": round(float(k), 4),
                 "grid_points_ok": n_ok,
                 "grid_points_diverged": int(N_GRID * N_GRID - n_ok),
                 "min_abs_F": float(np.nanmin(mag)),
                 "median_abs_F": float(np.nanmedian(mag))}
        # crude crossing detector: any grid cell where both components
        # change sign among its corners
        cross = 0
        for i in range(N_GRID - 1):
            for j in range(N_GRID - 1):
                q1 = c1[j:j + 2, i:i + 2]
                q2 = c2[j:j + 2, i:i + 2]
                if (np.all(np.isfinite(q1)) and np.all(np.isfinite(q2))
                        and np.nanmin(q1) < 0 < np.nanmax(q1)
                        and np.nanmin(q2) < 0 < np.nanmax(q2)):
                    cross += 1
        panel["cells_with_both_sign_changes"] = cross
        out_panels.append(panel)
        print(panel, flush=True)
        del cs1, cs2
    fig.suptitle(f"M6.4 existence probe, frozen-spec system, g={G} "
                 f"lam={LAM}: localized solution = white and red contours "
                 "crossing")
    fig.tight_layout()
    fig.savefig(PLOTS / "m6_4_existence_contours.png", dpi=140)
    plt.close(fig)
    out = {"task": "M6.4 N4b existence probe", "g": G, "lam": LAM,
           "grid": {"n": N_GRID, "lo": AB_LO, "hi": AB_HI},
           "panels": out_panels,
           "runtime_s": round(time.time() - t0, 1)}
    (DATA / "m6_4_existence_probe.json").write_text(
        json.dumps(out, indent=2, default=str))
    print(f"done in {out['runtime_s']} s")


if __name__ == "__main__":
    main()
