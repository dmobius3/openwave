"""M6.4 ADVERSARIAL AUDIT: claims A5 + A6 (June z20866581 census fingerprint
and the W-test vs Gelfand-Fomin construction).

Auditor-owned rerun, methodologically independent of the primary
(m6_4_gf_census.py):
  - integrator LSODA (primary: DOP853), tolerances 1e-9/1e-11;
  - background stored on a 2500-point grid with CubicSpline interpolation in
    the Jacobi RHS (primary: 800-point linear np.interp);
  - the two Jacobi solutions integrated as ONE 8-dim system;
  - observables by Simpson quadrature (primary: trapezoid);
  - my own crossing counters.

NOTE: the primary's m6_4_gf_census.json did not exist at audit time (the
conv-S leg of the primary run had not completed); the conv-T counts audited
here are those in research/data/m6_4_census_run.log and the method note.

A5 audits, under the tail-multiplier convention (start at r_max = 10 with
a = A*K0(ka r), a' = -A*ka*K1(ka r), integrate inward to r0 = 0.02):
  - printed localization check max(|a(rmax)|,|b(rmax)|) < 0.05: expected
    360/360 trivially (published: 340);
  - printed W-test (Q12 = -c/4, value ICs): primary 318/360, published 319;
  - best-point L/Q at (g=0.3, w=2.2, lam=0.5, A=0.1, B=0.2): primary 0.471,
    published 0.567; L/Q range (published [0.57, 2.60]);
  - window flip on a 24-combo stratified subsample rerun at r_max = 15.

A6 audits: the printed W construction (value ICs (1,0,0,0)/(0,0,1,0)) vs the
Gelfand-Fomin conjugate-point construction (solutions VANISHING at r0 with
independent derivative ICs (0,1,0,0)/(0,0,0,1); zeros of the value
determinant). GF (Calculus of Variations, 1963, Ch. 5 section 29) defines a
conjugate point r~ by a nontrivial Jacobi solution with xi(r0) = 0 = xi(r~);
detecting these REQUIRES the fundamental set vanishing at r0. Value-IC
solutions instead detect zeros of a determinant tied to a Neumann-type
condition at r0 -- a disconjugacy statement for a different boundary problem.
The two must disagree on some backgrounds; this script counts and lists the
disagreements across the full 360-combo grid.

Output: research/data/m6_4_audit_june_census.json
"""

import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import k0, k1

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"

C = 1.0
R0 = 0.02
BLOW = 100.0
G_GRID = [0.3, 0.5, 0.8, 1.0, 2.0]
W_GRID = [1.0, 1.2, 1.5, 1.8, 2.0, 2.2]
LAM_GRID = [0.0, 0.5, 1.0]
AB_GRID = [(0.1, 0.1), (0.1, 0.2), (0.2, 0.1), (0.2, 0.2)]
BEST = (0.3, 2.2, 0.5, 0.1, 0.2)


def rhs_bg(r, y, g, m2a, m2b):
    """June system (my transcription of z20866581 section 2.1), c = 1."""
    a, ap, b, bp = y
    return [ap, m2a * a - (C / 2) * b - 4 * g * a**3 - ap / r,
            bp, m2b * b - (C / 2) * a - 4 * g * b**3 - bp / r]


def background(g, w, lam, A, B, rmax):
    m2a, m2b = 1 + w**2, 1 + lam + w**2
    ka, kb = np.sqrt(m2a), np.sqrt(m2b)
    y0 = [A * k0(ka * rmax), -A * ka * k1(ka * rmax),
          B * k0(kb * rmax), -B * kb * k1(kb * rmax)]
    ev = np.linspace(rmax, R0, 2500)
    hit = lambda r, y, *a_: BLOW - np.max(np.abs(y))  # noqa: E731
    hit.terminal = True
    sol = solve_ivp(rhs_bg, (rmax, R0), y0, t_eval=ev, args=(g, m2a, m2b),
                    method="LSODA", rtol=1e-9, atol=1e-11, events=hit)
    if not sol.success or sol.t.shape[0] != ev.shape[0]:
        return None
    return sol.t[::-1], sol.y[0][::-1], sol.y[2][::-1]


def jac_rhs(r, x, g, m2a, m2b, q12, spl_a, spl_b):
    av, bv = spl_a(r), spl_b(r)
    q11 = m2a - 12 * g * av**2
    q22 = m2b - 12 * g * bv**2
    out = np.empty(8)
    for j in (0, 1):
        xa, xap, xb, xbp = x[4 * j: 4 * j + 4]
        out[4 * j] = xap
        out[4 * j + 1] = q11 * xa + q12 * xb - xap / r
        out[4 * j + 2] = xbp
        out[4 * j + 3] = q12 * xa + q22 * xb - xbp / r
    return out


def crossings(v, skip_below=None):
    if skip_below is not None:
        v = v[np.abs(v) > skip_below]
    s = np.sign(v)
    s = s[s != 0]
    return int(np.sum(s[:-1] * s[1:] < 0))


def stability_tests(rr, a, b, g, w, lam, rmax):
    m2a, m2b = 1 + w**2, 1 + lam + w**2
    spl_a, spl_b = CubicSpline(rr, a), CubicSpline(rr, b)
    ev = np.linspace(R0, rmax, 2000)
    out = {}
    # printed W-test: value ICs
    ic = np.zeros(8)
    ic[0] = 1.0  # sol 1 = (1,0,0,0)
    ic[6] = 1.0  # sol 2 = (0,0,1,0)
    s = solve_ivp(jac_rhs, (R0, rmax), ic, t_eval=ev,
                  args=(g, m2a, m2b, -C / 4, spl_a, spl_b),
                  method="LSODA", rtol=1e-9, atol=1e-11)
    if s.success and s.t.shape[0] == ev.shape[0]:
        wr = s.y[0] * s.y[6] - s.y[4] * s.y[2]
        out["w_cross"] = crossings(wr)
    else:
        out["w_cross"] = None
    # proper GF: vanishing ICs
    ic = np.zeros(8)
    ic[1] = 1.0  # sol 1 = (0,1,0,0)
    ic[7] = 1.0  # sol 2 = (0,0,0,1)
    s = solve_ivp(jac_rhs, (R0, rmax), ic, t_eval=ev,
                  args=(g, m2a, m2b, -C / 4, spl_a, spl_b),
                  method="LSODA", rtol=1e-9, atol=1e-11)
    if s.success and s.t.shape[0] == ev.shape[0]:
        det = s.y[0] * s.y[6] - s.y[4] * s.y[2]
        start = np.searchsorted(ev, R0 + 0.1)
        out["gf_conj"] = crossings(det[start:], skip_below=1e-14)
    else:
        out["gf_conj"] = None
    return out


def main():
    t0 = time.time()
    rows = []
    for g in G_GRID:
        for w in W_GRID:
            for lam in LAM_GRID:
                for A, B in AB_GRID:
                    row = dict(g=g, w=w, lam=lam, A=A, B=B)
                    bg = background(g, w, lam, A, B, 10.0)
                    if bg is None:
                        row["diverged"] = True
                        rows.append(row)
                        continue
                    rr, a, b = bg
                    row["diverged"] = False
                    row["abs_rmax"] = float(max(abs(a[-1]), abs(b[-1])))
                    row["loc_printed"] = bool(row["abs_rmax"] < 0.05)
                    q = 2 * np.pi * simpson(rr * b**2, x=rr)
                    ell = 2 * np.pi * simpson(rr * a * b, x=rr)
                    row["LoverQ"] = float(ell / q) if abs(q) > 0 else None
                    row.update(stability_tests(rr, a, b, g, w, lam, 10.0))
                    rows.append(row)
    ok = [r for r in rows if not r["diverged"]]

    # window-flip subsample: every 15th combo -> 24 stratified combos
    flips = []
    for r in rows[::15]:
        if r["diverged"] or r["w_cross"] is None:
            continue
        bg15 = background(r["g"], r["w"], r["lam"], r["A"], r["B"], 15.0)
        if bg15 is None:
            continue
        t15 = stability_tests(*bg15, r["g"], r["w"], r["lam"], 15.0)
        flips.append({
            "combo": {k: r[k] for k in ("g", "w", "lam", "A", "B")},
            "w_cross_rmax10": r["w_cross"],
            "w_cross_rmax15": t15["w_cross"],
            "flip": (r["w_cross"] == 0) != (t15["w_cross"] == 0)
            if t15["w_cross"] is not None else None})

    # best-point L/Q convergence study (tolerance + integrator sweep):
    # the conv-T background core is tolerance-sensitive, so pin the
    # converged value before comparing to primary 0.471 / published 0.567.
    conv = []
    for method, npts, tol in (("LSODA", 2500, 1e-9), ("LSODA", 8000, 1e-11),
                              ("Radau", 8000, 1e-10),
                              ("DOP853", 8000, 1e-12)):
        g, w, lam, A, B = BEST
        m2a, m2b = 1 + w**2, 1 + lam + w**2
        ka, kb = np.sqrt(m2a), np.sqrt(m2b)
        y0 = [A * k0(ka * 10.0), -A * ka * k1(ka * 10.0),
              B * k0(kb * 10.0), -B * kb * k1(kb * 10.0)]
        ev = np.linspace(10.0, R0, npts)
        s = solve_ivp(rhs_bg, (10.0, R0), y0, t_eval=ev, args=(g, m2a, m2b),
                      method=method, rtol=tol, atol=tol * 1e-2)
        rr2 = s.t[::-1]
        a2, b2 = s.y[0][::-1], s.y[2][::-1]
        q2 = simpson(rr2 * b2**2, x=rr2)
        e2 = simpson(rr2 * a2 * b2, x=rr2)
        conv.append({"method": method, "n": npts, "rtol": tol,
                     "LoverQ": float(e2 / q2),
                     "a_r0": float(a2[0]), "b_r0": float(b2[0])})

    disagree = [r for r in ok
                if r["w_cross"] is not None and r["gf_conj"] is not None
                and (r["w_cross"] == 0) != (r["gf_conj"] == 0)]
    best_row = [r for r in ok
                if (r["g"], r["w"], r["lam"], r["A"], r["B"]) == BEST][0]
    lq = [r["LoverQ"] for r in ok if r["LoverQ"] is not None]
    summary = {
        "n": len(rows),
        "diverged": sum(r["diverged"] for r in rows),
        "loc_printed_pass": sum(bool(r.get("loc_printed")) for r in ok),
        "w_test_pass": sum(r["w_cross"] == 0 for r in ok
                           if r["w_cross"] is not None),
        "gf_proper_pass": sum(r["gf_conj"] == 0 for r in ok
                              if r["gf_conj"] is not None),
        "published": {"localized": 340, "gf_pass": 319,
                      "LoverQ_range": [0.57, 2.60], "best_LoverQ": 0.567},
        "primary_convT": {"loc": 360, "w_pass": 318, "gf_proper": 338,
                          "best_LoverQ": 0.471},
        "best_point_LoverQ": best_row["LoverQ"],
        "best_point_LoverQ_convergence": conv,
        "best_point_w_cross": best_row["w_cross"],
        "LoverQ_min": float(min(lq)), "LoverQ_max": float(max(lq)),
        "LoverQ_negative_count": sum(v < 0 for v in lq),
        "window_flip_subsample_n": len(flips),
        "window_flips": sum(bool(f["flip"]) for f in flips),
        "flip_details": flips,
        "w_vs_gf_disagreements": len(disagree),
        "disagreement_examples": [
            {k: r[k] for k in ("g", "w", "lam", "A", "B",
                               "w_cross", "gf_conj")}
            for r in disagree[:10]],
    }
    out = {"task": "M6.4 audit A5/A6: June census independent rerun (conv T)",
           "summary": summary, "rows": rows,
           "runtime_s": round(time.time() - t0, 1)}
    (DATA / "m6_4_audit_june_census.json").write_text(
        json.dumps(out, indent=2, default=str))
    s2 = {k: v for k, v in summary.items() if k != "flip_details"}
    print(json.dumps(s2, indent=2, default=str))
    print(f"runtime {out['runtime_s']} s")


if __name__ == "__main__":
    main()
