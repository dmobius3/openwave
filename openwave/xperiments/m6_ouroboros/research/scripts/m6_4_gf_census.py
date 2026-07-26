"""M6.4 N1/N2: rerun the Zenodo 20866581 v2 Gelfand-Fomin census as printed.

The record's "attached" script (chaoiton_gf_verification.py) does not exist on
either Zenodo version (v1 = chaoiton_theorem.lean.txt, v2 = the docx), so this
is a reconstruction from the docx's printed recipe (research/data/
m6_4_record_v2_fulltext.txt), with every under-specified fork enumerated and
reported (pre-registered: no fork is resolved by matching the published counts).

Printed system (docx section 2.1), c = 1:
    alpha'' + (1/r) alpha' - m2a*alpha + (c/2)*beta + 4g*alpha^3 = 0
    beta''  + (1/r) beta'  - m2b*beta  + (c/2)*alpha + 4g*beta^3  = 0
    m2a = 1 + w^2,  m2b = 1 + lam + w^2

Printed test (section 2.3): Jacobi system with Q11 = m2a - 12g*a^2,
Q22 = m2b - 12g*b^2, Q12 = -c/4; two solutions from IC1 = (1,0,0,0),
IC2 = (0,0,1,0); W(r) = xa1*xb2 - xa2*xb1; a zero crossing = conjugate point.

Forks enumerated (the docx under-specifies the amplitude convention):
    T: alpha(rmax) = A*K0(ka*rmax)      (A = tail multiplier), inward
    V: alpha(rmax) = A                  (A = boundary value),   inward
    S: alpha(r0) = A*r0, alpha'(r0) = A (A = core slope, Lean A0), outward

Extra columns beyond the printed recipe (reported separately, never merged):
    - corrected Jacobi cross term Q12 = -c/2 (the true linearization of 2.1)
    - proper GF conjugate-point determinant (solutions VANISHING at r0,
      derivative ICs = identity; zeros of det after r0)
    - far-field channel analysis (both kappa^2 > 0 for all grid points: the
      flipped-sign system localizes structurally at every omega)
    - window sensitivity: printed-W crossings at rmax = 10 vs 15

Grid (section 3.1): g in {0.3,0.5,0.8,1.0,2.0}, w in {1.0,1.2,1.5,1.8,2.0,2.2},
lam in {0,0.5,1.0}, (A,B) in {(.1,.1),(.1,.2),(.2,.1),(.2,.2)} -> 360.
Published: 340 localized, 319 GF-pass, L/Q in [0.57, 2.60],
best L/Q = 0.567 at g=0.3, w=2.2, lam=0.5, A=0.1, B=0.2.

Outputs: research/data/m6_4_gf_census.json (+ per-convention npz),
research/plots/m6_4_census_*.png. Runtime: ~2 min single process.
"""

import json
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import k0, k1

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
PLOTS = HERE.parent / "plots"

C = 1.0
R0, RMAX = 0.02, 10.0
RTOL, ATOL = 1e-10, 1e-12
BLOWUP = 100.0  # |field| beyond this = diverged (reconstruction guard: fields
# of physical interest are O(1); the guard only CLASSIFIES rows as diverged,
# it never alters numbers of surviving rows. Set low so saturated inward
# integrations terminate instead of grinding at machine-small step sizes.

G_GRID = [0.3, 0.5, 0.8, 1.0, 2.0]
W_GRID = [1.0, 1.2, 1.5, 1.8, 2.0, 2.2]
LAM_GRID = [0.0, 0.5, 1.0]
AB_GRID = [(0.1, 0.1), (0.1, 0.2), (0.2, 0.1), (0.2, 0.2)]

BEST_POINT = dict(g=0.3, w=2.2, lam=0.5, A=0.1, B=0.2)  # published L/Q = 0.567


def rhs_background(r, y, g, m2a, m2b):
    a, ap, b, bp = y
    app = -ap / r + m2a * a - (C / 2.0) * b - 4.0 * g * a**3
    bpp = -bp / r + m2b * b - (C / 2.0) * a - 4.0 * g * b**3
    return [ap, app, bp, bpp]


def solve_background(conv, g, w, lam, A, B, rmax=RMAX):
    """Return (r_grid, a, b) on an even grid r0..rmax, or None if diverged."""
    m2a, m2b = 1.0 + w**2, 1.0 + lam + w**2
    ka, kb = np.sqrt(m2a), np.sqrt(m2b)
    r_eval = np.linspace(R0, rmax, 800)
    blowup = lambda r, y, *a_: BLOWUP - np.max(np.abs(y))  # noqa: E731
    blowup.terminal = True
    if conv == "S":  # outward from the core, Lean A0/B0 = slopes at 0
        y0 = [A * R0, A, B * R0, B]
        span, ev = (R0, rmax), r_eval
    else:  # inward from the K0 tail
        if conv == "T":
            a_m, ap_m = A * k0(ka * rmax), -A * ka * k1(ka * rmax)
            b_m, bp_m = B * k0(kb * rmax), -B * kb * k1(kb * rmax)
        else:  # V: boundary VALUE = A, slope from the K0 log-derivative
            a_m, ap_m = A, -A * ka * k1(ka * rmax) / k0(ka * rmax)
            b_m, bp_m = B, -B * kb * k1(kb * rmax) / k0(kb * rmax)
        y0 = [a_m, ap_m, b_m, bp_m]
        span, ev = (rmax, R0), r_eval[::-1]
    sol = solve_ivp(rhs_background, span, y0, t_eval=ev, args=(g, m2a, m2b),
                    method="DOP853", rtol=RTOL, atol=ATOL, events=blowup)
    if not sol.success or sol.t.shape[0] != r_eval.shape[0]:
        return None
    if conv == "S":
        return sol.t, sol.y[0], sol.y[2]
    return sol.t[::-1], sol.y[0][::-1], sol.y[2][::-1]


def rhs_jacobi(r, x, g, m2a, m2b, q12, r_grid, a_prof, b_prof):
    a = np.interp(r, r_grid, a_prof)
    b = np.interp(r, r_grid, b_prof)
    q11 = m2a - 12.0 * g * a**2
    q22 = m2b - 12.0 * g * b**2
    xa, xap, xb, xbp = x
    return [xap, -xap / r + q11 * xa + q12 * xb,
            xbp, -xbp / r + q12 * xa + q22 * xb]


def count_crossings(w_vals):
    s = np.sign(w_vals)
    s = s[s != 0]
    return int(np.sum(s[:-1] * s[1:] < 0))


def jacobi_tests(r_grid, a_prof, b_prof, g, w, lam, rmax=RMAX):
    """Printed W-test (Q12=-c/4), corrected (Q12=-c/2), and proper GF det."""
    m2a, m2b = 1.0 + w**2, 1.0 + lam + w**2
    out = {}
    r_eval = np.linspace(R0, rmax, 800)
    for tag, q12 in (("printed", -C / 4.0), ("linearized", -C / 2.0)):
        sols = []
        for ic in ([1, 0, 0, 0], [0, 0, 1, 0]):
            s = solve_ivp(rhs_jacobi, (R0, rmax), ic, t_eval=r_eval,
                          args=(g, m2a, m2b, q12, r_grid, a_prof, b_prof),
                          method="DOP853", rtol=1e-9, atol=1e-11)
            if not s.success:
                sols = None
                break
            sols.append(s.y)
        if sols is None:
            out[f"w_cross_{tag}"] = None
            continue
        wv = sols[0][0] * sols[1][2] - sols[1][0] * sols[0][2]
        out[f"w_cross_{tag}"] = count_crossings(wv)
        # proper GF: solutions VANISHING at r0, derivative ICs = identity
        sols_gf = []
        for ic in ([0, 1, 0, 0], [0, 0, 0, 1]):
            s = solve_ivp(rhs_jacobi, (R0, rmax), ic, t_eval=r_eval,
                          args=(g, m2a, m2b, q12, r_grid, a_prof, b_prof),
                          method="DOP853", rtol=1e-9, atol=1e-11)
            sols_gf.append(s.y if s.success else None)
        if all(s is not None for s in sols_gf):
            det = sols_gf[0][0] * sols_gf[1][2] - sols_gf[1][0] * sols_gf[0][2]
            # skip the trivial zero at r0 itself: start counting after det
            # first leaves 0 (det ~ (r-r0)^2 near r0)
            mask = np.abs(det) > 1e-12
            first = np.argmax(mask)
            out[f"gf_conjpts_{tag}"] = count_crossings(det[first:])
        else:
            out[f"gf_conjpts_{tag}"] = None
    return out


def observables(r_grid, a_prof, b_prof):
    q = 2 * np.pi * np.trapezoid(r_grid * b_prof**2, r_grid)
    ell = 2 * np.pi * np.trapezoid(r_grid * a_prof * b_prof, r_grid)
    return q, ell


def run_census(conv):
    rows = []
    for g in G_GRID:
        for w in W_GRID:
            for lam in LAM_GRID:
                for A, B in AB_GRID:
                    row = dict(conv=conv, g=g, w=w, lam=lam, A=A, B=B)
                    bg = solve_background(conv, g, w, lam, A, B)
                    if bg is None:
                        row.update(diverged=True)
                        rows.append(row)
                        continue
                    r_grid, a_prof, b_prof = bg
                    row["diverged"] = False
                    row["abs_at_rmax"] = float(
                        max(abs(a_prof[-1]), abs(b_prof[-1])))
                    row["loc_printed"] = bool(row["abs_at_rmax"] < 0.05)
                    row["core_amp"] = float(
                        max(np.max(np.abs(a_prof)), np.max(np.abs(b_prof))))
                    q, ell = observables(r_grid, a_prof, b_prof)
                    row["Q"], row["L"] = float(q), float(ell)
                    row["LoverQ"] = float(ell / q) if abs(q) > 1e-300 else None
                    row.update(jacobi_tests(r_grid, a_prof, b_prof, g, w, lam))
                    # window sensitivity on the printed W-test
                    bg15 = solve_background(conv, g, w, lam, A, B, rmax=15.0)
                    if bg15 is not None:
                        j15 = jacobi_tests(*bg15, g, w, lam, rmax=15.0)
                        row["w_cross_printed_rmax15"] = j15["w_cross_printed"]
                    rows.append(row)
    return rows


def summarize(rows):
    ok = [r for r in rows if not r["diverged"]]
    s = {
        "n": len(rows),
        "diverged": sum(r["diverged"] for r in rows),
        "loc_printed_pass": sum(r.get("loc_printed", False) for r in ok),
        "gf_printed_pass": sum(
            1 for r in ok if r.get("w_cross_printed") == 0),
        "gf_linearized_pass": sum(
            1 for r in ok if r.get("w_cross_linearized") == 0),
        "gf_proper_pass_printedQ": sum(
            1 for r in ok if r.get("gf_conjpts_printed") == 0),
        "gf_proper_pass_linQ": sum(
            1 for r in ok if r.get("gf_conjpts_linearized") == 0),
        "both_printed_pass": sum(
            1 for r in ok
            if r.get("loc_printed") and r.get("w_cross_printed") == 0),
        "window_flips_printed": sum(
            1 for r in ok
            if r.get("w_cross_printed_rmax15") is not None
            and (r["w_cross_printed"] == 0)
            != (r["w_cross_printed_rmax15"] == 0)),
    }
    lq = [r["LoverQ"] for r in ok if r.get("LoverQ") is not None]
    if lq:
        s["LoverQ_min"], s["LoverQ_max"] = float(min(lq)), float(max(lq))
    best = [r for r in ok
            if (r["g"], r["w"], r["lam"], r["A"], r["B"])
            == tuple(BEST_POINT.values())]
    s["best_point"] = best[0] if best else None
    return s


def main():
    t0 = time.time()
    all_rows, summaries = {}, {}
    for conv in ("T", "V", "S"):
        print(f"== convention {conv} ==", flush=True)
        rows = run_census(conv)
        all_rows[conv] = rows
        summaries[conv] = summarize(rows)
        print(json.dumps(summaries[conv], indent=2, default=str), flush=True)

    # fingerprint plot: background at the published best point, per convention
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax, conv in zip(axes, ("T", "V", "S")):
        bg = solve_background(conv, **BEST_POINT)
        if bg is not None:
            r_grid, a_prof, b_prof = bg
            ax.plot(r_grid, a_prof, label="alpha(r)")
            ax.plot(r_grid, b_prof, label="beta(r)")
        bp = summaries[conv]["best_point"]
        lq = bp.get("LoverQ") if bp else None
        lq_s = f"{lq:.4f}" if lq is not None else "n/a"
        ax.set_title(f"conv {conv}: L/Q = {lq_s} (published 0.567)")
        ax.set_xlabel("r")
        ax.legend()
        ax.grid(alpha=0.3)
    fig.suptitle("M6.4 June-census best point g=0.3 w=2.2 lam=0.5 A=0.1 B=0.2")
    fig.tight_layout()
    fig.savefig(PLOTS / "m6_4_census_best_point.png", dpi=140)
    plt.close(fig)

    out = {
        "task": "M6.4 N1/N2 June (z20866581 v2) census rerun",
        "published": {"n": 360, "localized": 340, "gf_pass": 319,
                      "LoverQ_range": [0.57, 2.60],
                      "best_LoverQ": 0.567, "best_point": BEST_POINT},
        "summaries": summaries,
        "runtime_s": round(time.time() - t0, 1),
    }
    (DATA / "m6_4_gf_census.json").write_text(
        json.dumps(out, indent=2, default=str))
    for conv, rows in all_rows.items():
        np.savez_compressed(DATA / f"m6_4_census_rows_{conv}.npz",
                            rows=json.dumps(rows, default=str))
    print(f"done in {out['runtime_s']} s -> m6_4_gf_census.json")


if __name__ == "__main__":
    main()
