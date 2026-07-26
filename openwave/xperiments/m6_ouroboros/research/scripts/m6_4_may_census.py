"""M6.4 D1-backing: reproduce the May z20044392 census (62 vs internal 42).

The May record's ODE (3.3)-(3.4) is EXACTLY the Lean-formalized system of
z20866581 section 1 (correct Klein-Gordon +w^2 sign, centrifugal -1/r^2,
cubic in beta only):

    alpha'' + alpha'/r - alpha/r^2 + w^2*alpha = beta
    beta''  + beta'/r  - beta/r^2  + w^2*beta  = alpha - lam*beta - 4g*beta^3

Printed grid: g in [0.5,5.0], w in [0.3,1.5], A0,B0 in [0.1,1.0], lam in
[0,2], "covering 1280 combinations on a regular grid". 1280 = 4*4*4*4*5;
reconstruction (documented, evenly spaced with endpoints):
    g in {0.5, 2.0, 3.5, 5.0}; w in {0.3, 0.7, 1.1, 1.5};
    A0, B0 in {0.1, 0.4, 0.7, 1.0}; lam in {0, 0.5, 1.0, 1.5, 2.0}.

Printed criteria (4.2): (i) |a|+|b| < 0.1 for all r >= 8; (ii) <= 4 radial
nodes of a or b; (iii) GF conjugate-point test on H' before r_max = 8.
Shooting from the regular core solution a ~ A0*r, b ~ B0*r (RK45, rtol 1e-8,
atol 1e-10 as printed).

Also verifies the printed far-field threshold analytically + the
representative case (g=0.5, w=0.9, lam=0): kappa = sqrt(X+ - w^2) ~ 0.436,
X+ = (-lam + sqrt(lam^2+4))/2.

Outputs: research/data/m6_4_may_census.json, research/plots/m6_4_may_*.png.
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

R0, RMAX = 0.02, 10.0  # integrate past r=8 so criterion (i) has a window
RTOL, ATOL = 1e-8, 1e-10  # as printed (section 4.1)
BLOWUP = 100.0

G_GRID = [0.5, 2.0, 3.5, 5.0]
W_GRID = [0.3, 0.7, 1.1, 1.5]
AB_VALS = [0.1, 0.4, 0.7, 1.0]
LAM_GRID = [0.0, 0.5, 1.0, 1.5, 2.0]


def x_plus(lam):
    return (-lam + np.sqrt(lam**2 + 4.0)) / 2.0


def rhs_may(r, y, g, w, lam):
    a, ap, b, bp = y
    app = -ap / r + a / r**2 - w**2 * a + b
    bpp = -bp / r + b / r**2 - w**2 * b + a - lam * b - 4.0 * g * b**3
    return [ap, app, bp, bpp]


def solve_may(g, w, lam, A0, B0, rmax=RMAX, method="RK45"):
    y0 = [A0 * R0, A0, B0 * R0, B0]  # regular branch: a ~ A0*r
    ev = np.linspace(R0, rmax, 1000)
    blowup = lambda r, y, *a_: BLOWUP - np.max(np.abs(y))  # noqa: E731
    blowup.terminal = True
    sol = solve_ivp(rhs_may, (R0, rmax), y0, t_eval=ev, args=(g, w, lam),
                    method=method, rtol=RTOL, atol=ATOL, events=blowup)
    if not sol.success or sol.t.shape[0] != ev.shape[0]:
        return None
    return sol.t, sol.y[0], sol.y[2]


def count_nodes(prof):
    s = np.sign(prof)
    s = s[np.abs(prof) > 1e-9]
    return int(np.sum(s[:-1] * s[1:] < 0))


def rhs_jacobi_may(r, x, g, w, lam, r_grid, b_prof):
    b = np.interp(r, r_grid, b_prof)
    xa, xap, xb, xbp = x
    xapp = -xap / r + xa / r**2 - w**2 * xa + xb
    xbpp = (-xbp / r + xb / r**2 - w**2 * xb + xa - lam * xb
            - 12.0 * g * b**2 * xb)
    return [xap, xapp, xbp, xbpp]


def gf_conjugate_points(r_grid, b_prof, g, w, lam, r_test=8.0):
    """Proper GF: two solutions vanishing at r0 (derivative ICs = identity);
    conjugate point = zero of the 2x2 value determinant in (r0, r_test]."""
    ev = np.linspace(R0, r_test, 800)
    sols = []
    for ic in ([0, 1, 0, 0], [0, 0, 0, 1]):
        s = solve_ivp(rhs_jacobi_may, (R0, r_test), ic, t_eval=ev,
                      args=(g, w, lam, r_grid, b_prof),
                      method="RK45", rtol=1e-8, atol=1e-10)
        if not s.success or s.t.shape[0] != ev.shape[0]:
            return None
        sols.append(s.y)
    det = sols[0][0] * sols[1][2] - sols[1][0] * sols[0][2]
    mask = np.abs(det) > 1e-12
    first = int(np.argmax(mask))
    d = det[first:]
    s = np.sign(d)
    s = s[s != 0]
    return int(np.sum(s[:-1] * s[1:] < 0))


def run_census():
    rows = []
    for g in G_GRID:
        for w in W_GRID:
            for lam in LAM_GRID:
                for A0 in AB_VALS:
                    for B0 in AB_VALS:
                        row = dict(g=g, w=w, lam=lam, A0=A0, B0=B0,
                                   below_threshold=bool(w**2 < x_plus(lam)))
                        bg = solve_may(g, w, lam, A0, B0)
                        if bg is None:
                            row["diverged"] = True
                            rows.append(row)
                            continue
                        r_grid, a_prof, b_prof = bg
                        row["diverged"] = False
                        tail = slice(np.searchsorted(r_grid, 8.0), None)
                        row["crit_i_loc"] = bool(
                            np.all(np.abs(a_prof[tail]) + np.abs(b_prof[tail])
                                   < 0.1))
                        row["crit_ii_nodes"] = bool(
                            count_nodes(a_prof) <= 4
                            and count_nodes(b_prof) <= 4)
                        if row["crit_i_loc"] and row["crit_ii_nodes"]:
                            cp = gf_conjugate_points(r_grid, b_prof, g, w, lam)
                            row["gf_conjpts"] = cp
                            row["crit_iii_gf"] = bool(cp == 0) \
                                if cp is not None else None
                        rows.append(row)
    return rows


def representative_case():
    """g=0.5, w=0.9, lam=0: printed kappa ~ 0.436, range ~ 2.3."""
    g, w, lam = 0.5, 0.9, 0.0
    kappa_pred = float(np.sqrt(x_plus(lam) - w**2))
    # small-amplitude probe: in the linear regime the decaying channel decays
    # at kappa_pred; fit the tail of the best grid solution
    best, best_tail = None, np.inf
    for A0 in AB_VALS:
        for B0 in AB_VALS:
            bg = solve_may(g, w, lam, A0, B0)
            if bg is None:
                continue
            r_grid, a_prof, b_prof = bg
            tail_mag = np.max(np.abs(a_prof[-100:]) + np.abs(b_prof[-100:]))
            if tail_mag < best_tail:
                best_tail, best = tail_mag, (A0, B0, r_grid, a_prof, b_prof)
    out = {"kappa_printed": 0.436, "kappa_analytic": round(kappa_pred, 4)}
    if best is not None:
        A0, B0, r_grid, a_prof, b_prof = best
        out["best_grid_shot"] = {"A0": A0, "B0": B0,
                                 "tail_mag_at_10": float(best_tail)}
        m = np.abs(b_prof) + np.abs(a_prof)
        sel = (r_grid > 5.0) & (m > 1e-12)
        if np.sum(sel) > 10:
            slope = np.polyfit(r_grid[sel], np.log(m[sel] * np.sqrt(
                r_grid[sel])), 1)[0]
            out["best_grid_shot"]["tail_decay_fit"] = round(float(-slope), 4)
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.semilogy(r_grid, np.abs(a_prof) + 1e-16, label="|alpha|")
        ax.semilogy(r_grid, np.abs(b_prof) + 1e-16, label="|beta|")
        ref = m[np.searchsorted(r_grid, 5.0)]
        ax.semilogy(r_grid, ref * np.exp(-kappa_pred * (r_grid - 5.0)),
                    "k--", label=f"exp(-{kappa_pred:.3f} r) (printed 0.436)")
        ax.set_title("May representative case g=0.5 w=0.9 lam=0, "
                     f"best grid shot A0={A0} B0={B0}")
        ax.set_xlabel("r")
        ax.legend()
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(PLOTS / "m6_4_may_representative.png", dpi=140)
        plt.close(fig)
    return out


def main():
    t0 = time.time()
    thr = {str(lam): round(float(np.sqrt(x_plus(lam))), 4)
           for lam in LAM_GRID}
    print("localization thresholds w_max(lam) =", thr, flush=True)
    rows = run_census()
    ok = [r for r in rows if not r["diverged"]]
    passed = [r for r in ok if r.get("crit_i_loc") and r.get("crit_ii_nodes")
              and r.get("crit_iii_gf")]
    summary = {
        "n": len(rows),
        "diverged": sum(r["diverged"] for r in rows),
        "crit_i_pass": sum(bool(r.get("crit_i_loc")) for r in ok),
        "crit_i_and_ii_pass": sum(
            bool(r.get("crit_i_loc") and r.get("crit_ii_nodes"))
            for r in ok),
        "all_three_pass": len(passed),
        "passers_below_threshold": sum(
            r["below_threshold"] for r in passed),
        "published_counts": {"sec_5_1": 62, "sec_5_4_internal": 42},
        "thresholds_w_max": thr,
        "grid_reconstruction": "4g x 4w x 4A0 x 4B0 x 5lam = 1280 evenly "
                               "spaced incl. endpoints (factorization "
                               "unprinted in the record)",
    }
    # where do passers sit in omega?
    if passed:
        summary["passer_omegas"] = sorted({r["w"] for r in passed})
        summary["passer_examples"] = [
            {k: r[k] for k in ("g", "w", "lam", "A0", "B0")}
            for r in passed[:10]]
    rep = representative_case()
    out = {"task": "M6.4 D1: May z20044392 census reconstruction",
           "summary": summary, "representative_case": rep,
           "runtime_s": round(time.time() - t0, 1)}
    (DATA / "m6_4_may_census.json").write_text(
        json.dumps(out, indent=2, default=str))
    np.savez_compressed(DATA / "m6_4_may_census_rows.npz",
                        rows=json.dumps(rows, default=str))
    print(json.dumps(summary, indent=2, default=str))
    print("representative:", json.dumps(rep, indent=2, default=str))
    print(f"done in {out['runtime_s']} s")


if __name__ == "__main__":
    main()
