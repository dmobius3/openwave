"""M6.4 N4: the omega-selection hunt (canonical OQ3).

Question: does ANY mechanism in the frozen term set make the lepton omega
values (electron 1, mu 11.0/12.82, tau 40.7/50.0, pion 15.0) discrete?

System (the frozen v11-class reduction = May z20044392 (3.3)-(3.4) = the Lean
ODE of z20866581 section 1; correct Klein-Gordon +w^2 sign):

    alpha'' + alpha'/r - alpha/r^2 + w^2*alpha = beta
    beta''  + beta'/r  - beta/r^2  + w^2*beta  = alpha - lam*beta - 4g*beta^3

Far-field channel structure (analytic, verified numerically here): linearized
about the vacuum, Lap1(a,b) = M(a,b) with M = [[-w^2, 1], [1, -w^2-lam]];
eigenvalues mu2_pm = X_pm - w^2, X_pm = (-lam +/- sqrt(lam^2+4))/2.
X_minus < 0 always -> one OSCILLATORY channel at every (w, lam).
X_plus in (0, 1] for lam >= 0 -> the decaying channel exists iff w^2 < X_plus
(the May record's own printed threshold; w_max(lam) <= 1 for lam >= 0).

True localization therefore requires the solution to have ZERO amplitude in
the oscillatory channel: 2 asymptotic conditions on the 2 shooting DOF
(A0, B0) -> at fixed (g, lam, w), localized solutions sit at ISOLATED
amplitude pairs. The omega-discreteness question is whether those solutions
exist only at isolated w (mechanism!) or across a continuum of w in the
window (clean negative: the term set selects AMPLITUDES, never frequencies).

Method: minimize the oscillatory-channel tail amplitude over (A0, B0) at each
w on a fine grid across the window, with continuation seeding; a solution
"exists" if the minimized residual is < TOL_EXIST and the profile's decaying
tail tracks kappa(w) = sqrt(X_plus - w^2). Two panels: (g=0.5, lam=0) and
(g=1.0, lam=1.0). Pre-registered: discreteness = isolated w with existence,
separated by gaps of non-existence resolved by the grid; a continuum of
existence, or uniform non-existence, is a clean negative.

Also states the flipped-sign (June z20866581 integrated) system's structure:
both far-field channels decay for every w (kappa^2 eigenvalues of
[[1+w^2, -1/2], [-1/2, 1+lam+w^2]] are positive for all w, lam >= 0), so
EVERY omega admits decaying tails by construction: no selection mechanism is
possible in the system that produced the mu/tau ladder.

Outputs: research/data/m6_4_omega_selection.json,
research/plots/m6_4_omega_window.png, m6_4_omega_existence.png.
"""

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
PLOTS = HERE.parent / "plots"

R0, RFAR = 0.02, 30.0
TOL_EXIST = 1e-3  # minimized oscillatory-channel envelope amplitude
BLOWUP = 100.0

LADDER = {"electron (v11, w=1)": 1.0, "muon May16": 11.0, "muon May23": 12.82,
          "pion": 15.0, "tau May16": 40.7, "tau May23": 50.0}


def x_pm(lam):
    s = np.sqrt(lam**2 + 4.0)
    return (-lam + s) / 2.0, (-lam - s) / 2.0


def eigvecs(lam):
    """Unit eigenvectors of M = [[0,1],[1,-lam]] shifted by -w^2*I (shift
    does not change eigenvectors)."""
    xp, xm = x_pm(lam)
    vp = np.array([1.0, xp])
    vm = np.array([1.0, xm])
    return vp / np.linalg.norm(vp), vm / np.linalg.norm(vm)


def rhs(r, y, g, w, lam):
    a, ap, b, bp = y
    app = -ap / r + a / r**2 - w**2 * a + b
    bpp = -bp / r + b / r**2 - w**2 * b + a - lam * b - 4.0 * g * b**3
    return [ap, app, bp, bpp]


def shoot(g, w, lam, A0, B0, rfar=RFAR):
    y0 = [A0 * R0, A0, B0 * R0, B0]
    ev = np.linspace(R0, rfar, 1500)
    blow = lambda r, y, *a_: BLOWUP - np.max(np.abs(y))  # noqa: E731
    blow.terminal = True
    sol = solve_ivp(rhs, (R0, rfar), y0, t_eval=ev, args=(g, w, lam),
                    method="DOP853", rtol=1e-9, atol=1e-11, events=blow)
    if not sol.success or sol.t.shape[0] != ev.shape[0]:
        return None
    return sol.t, sol.y


def osc_amplitude(r_grid, y, lam):
    """Envelope amplitude of the oscillatory channel over the outer tail,
    2D-normalized (fields ~ amp/sqrt(r) in that channel)."""
    _, vm = eigvecs(lam)
    tail = r_grid > 0.6 * r_grid[-1]
    comp = vm[0] * y[0][tail] + vm[1] * y[2][tail]
    return float(np.max(np.abs(comp) * np.sqrt(r_grid[tail])))


def dec_amplitude(r_grid, y, lam):
    vp, _ = eigvecs(lam)
    tail = r_grid > 0.6 * r_grid[-1]
    comp = vp[0] * y[0][tail] + vp[1] * y[2][tail]
    return float(np.max(np.abs(comp)))


def existence_at(g, w, lam, seed):
    """Minimize the oscillatory-channel amplitude over (A0, B0). Returns
    (residual, (A0, B0), kappa_fit) or (inf, seed, None) if all shots blow."""

    def objective(p):
        out = shoot(g, w, lam, p[0], p[1])
        if out is None:
            return 1e3
        return osc_amplitude(*out, lam=lam)

    res = minimize(objective, seed, method="Nelder-Mead",
                   options=dict(xatol=1e-8, fatol=1e-10, maxiter=400))
    out = shoot(g, w, lam, res.x[0], res.x[1])
    kappa_fit = None
    if out is not None:
        r_grid, y = out
        vp, _ = eigvecs(lam)
        comp = np.abs(vp[0] * y[0] + vp[1] * y[2])
        sel = (r_grid > 8.0) & (r_grid < 20.0) & (comp > 1e-14)
        if np.sum(sel) > 20:
            kappa_fit = -float(np.polyfit(
                r_grid[sel], np.log(comp[sel] * np.sqrt(r_grid[sel])), 1)[0])
        # degenerate-minimum guard: the zero solution kills every channel
        if dec_amplitude(r_grid, y, lam) < 1e-10:
            return np.inf, tuple(res.x), None
    return float(res.fun), tuple(res.x), kappa_fit


def scan_panel(g, lam, n_w=24):
    xp, _ = x_pm(lam)
    w_max = float(np.sqrt(xp))
    w_grid = np.linspace(0.15 * w_max, 0.985 * w_max, n_w)
    seed = (0.3, 0.3)
    rows = []
    for w in w_grid:
        resid, amps, kappa_fit = existence_at(g, w, lam, seed)
        kappa_pred = float(np.sqrt(max(xp - w**2, 0.0)))
        rows.append(dict(w=float(w), residual=resid,
                         A0=amps[0], B0=amps[1],
                         kappa_fit=kappa_fit, kappa_pred=kappa_pred,
                         exists=bool(resid < TOL_EXIST)))
        if resid < TOL_EXIST:
            seed = amps  # continuation
        print(f"  g={g} lam={lam} w={w:.3f}: resid={resid:.2e} "
              f"exists={resid < TOL_EXIST} kfit={kappa_fit} "
              f"kpred={kappa_pred:.3f}", flush=True)
    return dict(g=g, lam=lam, w_max=w_max, rows=rows)


def june_flipped_structure():
    """Both kappa^2 eigenvalues of the June system's far-field matrix,
    sampled across the ladder omegas: positive everywhere."""
    out = {}
    for name, w in LADDER.items():
        for lam in (0.0, 0.5, 1.0):
            mat = np.array([[1 + w**2, -0.5], [-0.5, 1 + lam + w**2]])
            ev = np.linalg.eigvalsh(mat)
            out[f"{name}, lam={lam}"] = [round(float(e), 4) for e in ev]
    return out


def main():
    t0 = time.time()
    lam_axis = np.linspace(0, 2, 100)
    w_max_axis = np.sqrt([x_pm(la)[0] for la in lam_axis])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(lam_axis, w_max_axis, lw=2,
            label=r"$\omega_{max}(\lambda)=\sqrt{X_+}$ (frozen spec)")
    ax.fill_between(lam_axis, 0, w_max_axis, alpha=0.2,
                    label="localization window")
    for name, w in LADDER.items():
        ax.axhline(w, ls=":", lw=1, color="crimson")
        ax.annotate(f"{name}", (1.02, w), fontsize=7,
                    va="bottom", color="crimson")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$\omega$ (log)")
    ax.set_title("The claimed omega ladder vs the frozen-spec\n"
                 "localization window (electron sits ON the boundary)")
    ax.legend(loc="center right")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "m6_4_omega_window.png", dpi=140)
    plt.close(fig)

    panels = [scan_panel(0.5, 0.0), scan_panel(1.0, 1.0)]

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=False)
    for ax, panel in zip(axes, panels):
        wv = [r["w"] for r in panel["rows"]]
        rv = [r["residual"] for r in panel["rows"]]
        ax.semilogy(wv, rv, "o-", ms=4)
        ax.axhline(TOL_EXIST, color="k", ls="--", lw=1,
                   label=f"existence tol {TOL_EXIST}")
        ax.set_title(f"g={panel['g']}, lam={panel['lam']}: minimized "
                     "oscillatory-channel amplitude vs omega "
                     f"(window edge {panel['w_max']:.3f})")
        ax.set_xlabel("omega")
        ax.set_ylabel("residual (log)")
        ax.legend()
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOTS / "m6_4_omega_existence.png", dpi=140)
    plt.close(fig)

    out = {
        "task": "M6.4 N4: omega-selection hunt (OQ3)",
        "window": {"formula": "w_max = sqrt((-lam+sqrt(lam^2+4))/2)",
                   "w_max_lam0": 1.0,
                   "note": "w_max <= 1 for all lam >= 0; the record's own "
                           "grids use lam in [0,2]"},
        "ladder_vs_window": {
            name: {"w": w, "in_window_any_lam_0_2": bool(w < 1.0)}
            for name, w in LADDER.items()},
        "panels": panels,
        "june_flipped_sign_kappa2": june_flipped_structure(),
        "runtime_s": round(time.time() - t0, 1),
    }
    (DATA / "m6_4_omega_selection.json").write_text(
        json.dumps(out, indent=2, default=str))
    n_exist = [sum(r["exists"] for r in p["rows"]) for p in panels]
    print(f"existence counts per panel (of {len(panels[0]['rows'])} omegas):",
          n_exist)
    print(f"done in {out['runtime_s']} s")


if __name__ == "__main__":
    main()
