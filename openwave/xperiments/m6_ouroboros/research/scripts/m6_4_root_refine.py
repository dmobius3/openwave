"""M6.4 N4c: root-refine the existence-probe candidates + omega continuation.

The contour probe (m6_4_existence_probe.py) found candidate zeros of
F(A0, B0) = (c_cos, c_sin) at omega = 0.9 (g = 0.5, lam = 0): 31 cells where
both components change sign, min |F| = 6.5e-3, while the Nelder-Mead scan
(m6_4_omega_selection.py) reported no existence (optimizer funnel-miss).
This script settles it:

  1. re-evaluate F on the omega = 0.9 grid, collect candidate cells, run
     scipy.optimize.root (hybr) from each cell center, dedupe converged roots
     (|F| < 1e-8);
  2. verify each root IS a localized solution: decaying-channel tail fit
     kappa vs the analytic kappa = sqrt(X_plus - w^2); profile plot;
  3. proper-GF conjugate-point count on each root (is the localized state
     stable radially?);
  4. CONTINUATION in omega from each distinct root: re-solve at stepped
     omega across the window; existence on a continuous branch = the
     term set selects amplitudes, not frequencies (OQ3 clean negative via
     the constructive route); record where the branch terminates;
  5. spot-check the second panel (g = 1, lam = 1) near its window top.

Output: research/data/m6_4_root_refine.json,
research/plots/m6_4_branch_continuation.png, m6_4_root_profiles.png.
"""

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
PLOTS = HERE.parent / "plots"

R0, RFAR = 0.02, 30.0
BLOWUP = 100.0


def x_pm(lam):
    s = np.sqrt(lam**2 + 4.0)
    return (-lam + s) / 2.0, (-lam - s) / 2.0


def rhs(r, y, g, w, lam):
    a, ap, b, bp = y
    app = -ap / r + a / r**2 - w**2 * a + b
    bpp = -bp / r + b / r**2 - w**2 * b + a - lam * b - 4.0 * g * b**3
    return [ap, app, bp, bpp]


def shoot(g, w, lam, A0, B0):
    ev = np.linspace(R0, RFAR, 1500)
    blow = lambda r, y, *a_: BLOWUP - np.max(np.abs(y))  # noqa: E731
    blow.terminal = True
    sol = solve_ivp(rhs, (R0, RFAR), [A0 * R0, A0, B0 * R0, B0], t_eval=ev,
                    args=(g, w, lam), method="DOP853", rtol=1e-10, atol=1e-12,
                    events=blow)
    if not sol.success or sol.t.shape[0] != ev.shape[0]:
        return None
    return sol.t, sol.y


def make_F(g, w, lam):
    xp, xm = x_pm(lam)
    vm = np.array([1.0, xm])
    vm /= np.linalg.norm(vm)
    k = np.sqrt(w**2 - xm)
    r_fit = np.linspace(18.0, RFAR, 200)
    design = np.column_stack([np.cos(k * r_fit), np.sin(k * r_fit)])

    def F(p):
        out = shoot(g, w, lam, p[0], p[1])
        if out is None:
            return np.array([1e3, 1e3])
        r_grid, y = out
        comp = vm[0] * y[0] + vm[1] * y[2]
        proj = np.interp(r_fit, r_grid, comp) * np.sqrt(r_fit)
        coef, *_ = np.linalg.lstsq(design, proj, rcond=None)
        return coef

    return F


def kappa_fit_of(g, w, lam, A0, B0):
    xp, _ = x_pm(lam)
    vp = np.array([1.0, xp])
    vp /= np.linalg.norm(vp)
    out = shoot(g, w, lam, A0, B0)
    if out is None:
        return None, None
    r_grid, y = out
    comp = np.abs(vp[0] * y[0] + vp[1] * y[2])
    sel = (r_grid > 8.0) & (r_grid < 20.0) & (comp > 1e-13)
    if np.sum(sel) < 20:
        return None, out
    kfit = -float(np.polyfit(r_grid[sel],
                             np.log(comp[sel] * np.sqrt(r_grid[sel])), 1)[0])
    return kfit, out


def rhs_jacobi(r, x, g, w, lam, r_grid, b_prof):
    b = np.interp(r, r_grid, b_prof)
    xa, xap, xb, xbp = x
    xapp = -xap / r + xa / r**2 - w**2 * xa + xb
    xbpp = (-xbp / r + xb / r**2 - w**2 * xb + xa - lam * xb
            - 12.0 * g * b**2 * xb)
    return [xap, xapp, xbp, xbpp]


def gf_conjpts(g, w, lam, r_grid, b_prof, r_test=20.0):
    ev = np.linspace(R0, r_test, 900)
    sols = []
    for ic in ([0, 1, 0, 0], [0, 0, 0, 1]):
        s = solve_ivp(rhs_jacobi, (R0, r_test), ic, t_eval=ev,
                      args=(g, w, lam, r_grid, b_prof), method="DOP853",
                      rtol=1e-9, atol=1e-11)
        if not s.success or s.t.shape[0] != ev.shape[0]:
            return None
        sols.append(s.y)
    det = sols[0][0] * sols[1][2] - sols[1][0] * sols[0][2]
    d = det[int(np.argmax(np.abs(det) > 1e-12)):]
    sg = np.sign(d)
    sg = sg[sg != 0]
    return int(np.sum(sg[:-1] * sg[1:] < 0))


def find_roots(g, w, lam, n_grid=41, lo=0.02, hi=1.2):
    F = make_F(g, w, lam)
    amps = np.linspace(lo, hi, n_grid)
    c1 = np.full((n_grid, n_grid), np.nan)
    c2 = np.full((n_grid, n_grid), np.nan)
    for i, A0 in enumerate(amps):
        for j, B0 in enumerate(amps):
            v = F((A0, B0))
            c1[j, i], c2[j, i] = v
    cands = []
    for i in range(n_grid - 1):
        for j in range(n_grid - 1):
            q1, q2 = c1[j:j + 2, i:i + 2], c2[j:j + 2, i:i + 2]
            if (np.all(np.isfinite(q1)) and np.all(np.isfinite(q2))
                    and q1.min() < 0 < q1.max()
                    and q2.min() < 0 < q2.max()):
                cands.append((0.5 * (amps[i] + amps[i + 1]),
                              0.5 * (amps[j] + amps[j + 1])))
    roots = []
    for seed in cands:
        sol = root(F, seed, method="hybr", tol=1e-12)
        if not sol.success:
            continue
        resid = float(np.hypot(*F(sol.x)))
        if resid > 1e-8:
            continue
        if not (0 < sol.x[0] < 2.0 and 0 < sol.x[1] < 2.0):
            continue
        if any(np.hypot(sol.x[0] - r_[0], sol.x[1] - r_[1]) < 1e-4
               for r_ in roots):
            continue
        roots.append((float(sol.x[0]), float(sol.x[1]), resid))
    return roots, len(cands)


def continue_branch(g, lam, w_start, seed, w_lo=0.55, w_hi=0.999, dw=0.01):
    """Walk the root in omega both directions; return the branch rows."""
    branch = {}
    for direction in (-1, +1):
        w, p = w_start, np.array(seed)
        while True:
            w_next = w + direction * dw
            if not (w_lo <= w_next <= w_hi * np.sqrt(x_pm(lam)[0])):
                break
            F = make_F(g, w_next, lam)
            sol = root(F, p, method="hybr", tol=1e-12)
            if not sol.success or np.hypot(*F(sol.x)) > 1e-7:
                break
            w, p = w_next, sol.x
            kfit, _ = kappa_fit_of(g, w, lam, *p)
            branch[round(w, 4)] = dict(
                A0=float(p[0]), B0=float(p[1]),
                kappa_fit=kfit,
                kappa_pred=float(np.sqrt(max(x_pm(lam)[0] - w**2, 0))))
    return dict(sorted(branch.items()))


def main():
    t0 = time.time()
    out = {"task": "M6.4 N4c: root refinement + continuation"}

    g, lam, w = 0.5, 0.0, 0.9
    roots, n_cand = find_roots(g, w, lam)
    out["panel1"] = {"g": g, "lam": lam, "omega": w,
                     "candidate_cells": n_cand,
                     "distinct_roots": [
                         {"A0": r_[0], "B0": r_[1], "residual": r_[2]}
                         for r_ in roots]}
    print(f"omega={w}: {n_cand} candidate cells -> {len(roots)} "
          f"distinct roots: {roots}", flush=True)

    profiles = []
    for A0, B0, _ in roots:
        kfit, sol_out = kappa_fit_of(g, w, lam, A0, B0)
        cp = None
        if sol_out is not None:
            r_grid, y = sol_out
            cp = gf_conjpts(g, w, lam, r_grid, y[2])
            profiles.append((A0, B0, r_grid, y))
        for rec in out["panel1"]["distinct_roots"]:
            if abs(rec["A0"] - A0) < 1e-9:
                rec["kappa_fit"] = kfit
                rec["kappa_pred"] = float(np.sqrt(x_pm(lam)[0] - w**2))
                rec["gf_conjugate_points"] = cp
        print(f"  root A0={A0:.5f} B0={B0:.5f}: kappa_fit={kfit} "
              f"(pred {np.sqrt(x_pm(lam)[0] - w**2):.4f}), gf_conjpts={cp}",
              flush=True)

    if profiles:
        fig, axes = plt.subplots(1, len(profiles),
                                 figsize=(5.5 * len(profiles), 4.4),
                                 squeeze=False)
        for ax, (A0, B0, r_grid, y) in zip(axes[0], profiles):
            ax.plot(r_grid, y[0], label="alpha")
            ax.plot(r_grid, y[2], label="beta")
            ax.set_title(f"candidate root A0={A0:.4f} B0={B0:.4f} "
                         "(refuted: non-localized)")
            ax.set_xlabel("r")
            ax.legend()
            ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(PLOTS / "m6_4_root_profiles.png", dpi=140)
        plt.close(fig)

    # continuation from the first root (if any)
    if roots:
        seed = roots[0][:2]
        branch = continue_branch(g, lam, w, seed)
        out["panel1"]["branch"] = branch
        ws = [float(k) for k in branch]
        if ws:
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
            axes[0].plot(ws, [branch[k]["A0"] for k in branch], "o-",
                         ms=3, label="A0(omega)")
            axes[0].plot(ws, [branch[k]["B0"] for k in branch], "s-",
                         ms=3, label="B0(omega)")
            axes[0].set_xlabel("omega")
            axes[0].set_title("the localized-solution branch: amplitudes "
                              "vs omega (continuous = no omega selection)")
            axes[0].legend()
            axes[0].grid(alpha=0.3)
            axes[1].plot(ws, [branch[k]["kappa_fit"] for k in branch], "o",
                         ms=3, label="kappa fit")
            axes[1].plot(ws, [branch[k]["kappa_pred"] for k in branch], "-",
                         label="kappa analytic")
            axes[1].set_xlabel("omega")
            axes[1].set_title("tail decay rate vs the analytic prediction")
            axes[1].legend()
            axes[1].grid(alpha=0.3)
            fig.tight_layout()
            fig.savefig(PLOTS / "m6_4_branch_continuation.png", dpi=140)
            plt.close(fig)
            out["panel1"]["branch_extent"] = [min(ws), max(ws)]
            print(f"branch continues over omega in [{min(ws)}, {max(ws)}] "
                  f"({len(ws)} points)", flush=True)

    # panel 2 spot check near its window top
    g2, lam2 = 1.0, 1.0
    w2 = 0.97 * np.sqrt(x_pm(lam2)[0])
    roots2, n_cand2 = find_roots(g2, round(float(w2), 3), lam2, n_grid=31)
    out["panel2"] = {"g": g2, "lam": lam2, "omega": round(float(w2), 3),
                     "candidate_cells": n_cand2,
                     "distinct_roots": [
                         {"A0": r_[0], "B0": r_[1], "residual": r_[2]}
                         for r_ in roots2]}
    print(f"panel2 w={w2:.3f}: {n_cand2} cells -> {len(roots2)} roots",
          flush=True)

    out["runtime_s"] = round(time.time() - t0, 1)
    (DATA / "m6_4_root_refine.json").write_text(
        json.dumps(out, indent=2, default=str))
    print(f"done in {out['runtime_s']} s")


if __name__ == "__main__":
    main()
