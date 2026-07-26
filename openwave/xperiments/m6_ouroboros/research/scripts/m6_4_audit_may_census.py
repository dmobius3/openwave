"""M6.4 ADVERSARIAL AUDIT: claim A4 (May z20044392 census irreproducibility).

Auditor-owned reconstruction, methodologically independent of the primary
(m6_4_may_census.py):
  - integrator: LSODA (Adams/BDF) instead of the primary's RK45, run at BOTH
    the printed tolerances (rtol 1e-8, atol 1e-10) and tighter (1e-10/1e-12);
    plus a hand-rolled fixed-step classical RK4 cross-check on the printed
    representative combo (no scipy in that path);
  - Jacobi/GF machinery written from my own linearization of the May system,
    integrated as a single 8-dim system (both vanishing-IC solutions at once)
    on a cubic-spline background (primary: separate 4-dim solves on linear
    interpolation);
  - my own node counting and criteria implementation.

Subsample (>= 100 combos, self-chosen):
  - the full A0 = B0 = 0.1 slice: 4g x 4w x 5lam = 80 combos (this contains
    all 26 crit-i near-passers the primary reports);
  - a below-threshold block: w in {0.3, 0.7}, lam in {0, 1, 2}, (A0, B0) in
    {0.1, 0.4, 1.0}^2 minus overlaps -> the localized-window region;
  - the printed representative "stable" combo (g=0.5, w=1.5, A0=B0=0.1,
    lam=0) plus 8 random +-5% perturbations of it (the paper's section 5.3
    claims a >= 5% stable neighborhood around every passer).

Claim audited: under the printed criteria (i) |a|+|b| < 0.1 for r >= 8,
(ii) <= 4 radial nodes, (iii) GF conjugate-point test, the count of full
passers is 0 (published: 62 in section 5.1, 42 in section 5.4).

Output: research/data/m6_4_audit_may_census.json
"""

import itertools
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"

R0, RMAX, RTEST = 0.02, 10.0, 8.0
BLOW = 100.0


def rhs(r, y, g, w, lam):
    """May system (my transcription of z20044392 (3.3)-(3.4))."""
    a, ap, b, bp = y
    return [ap,
            b - w**2 * a + a / r**2 - ap / r,
            bp,
            a - lam * b - 4 * g * b**3 - w**2 * b + b / r**2 - bp / r]


def integrate(g, w, lam, A0, B0, rtol, atol):
    ev = np.linspace(R0, RMAX, 2000)
    hit = lambda r, y, *a_: BLOW - np.max(np.abs(y))  # noqa: E731
    hit.terminal = True
    sol = solve_ivp(rhs, (R0, RMAX), [A0 * R0, A0, B0 * R0, B0],
                    t_eval=ev, args=(g, w, lam), method="LSODA",
                    rtol=rtol, atol=atol, events=hit)
    if not sol.success or sol.t.shape[0] != ev.shape[0]:
        return None
    return sol.t, sol.y


def rk4_integrate(g, w, lam, A0, B0, h=5e-4):
    """Scipy-free fixed-step RK4 cross-check."""
    n = int((RMAX - R0) / h)
    rr = R0 + h * np.arange(n + 1)
    y = np.array([A0 * R0, A0, B0 * R0, B0], dtype=float)
    out = np.empty((n + 1, 4))
    out[0] = y
    f = lambda r, y_: np.array(rhs(r, y_, g, w, lam))  # noqa: E731
    for i in range(n):
        r_i = rr[i]
        k1 = f(r_i, y)
        k2 = f(r_i + h / 2, y + h / 2 * k1)
        k3 = f(r_i + h / 2, y + h / 2 * k2)
        k4 = f(r_i + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = y
        if np.max(np.abs(y)) > BLOW:
            return None
    return rr, out.T


def my_nodes(rr, prof):
    """Sign changes of the profile where it is meaningfully nonzero."""
    m = prof[np.abs(prof) > 1e-8]
    s = np.sign(m)
    return int(np.sum(s[:-1] * s[1:] < 0))


def crit_i(rr, a, b):
    sel = rr >= 8.0
    return bool(np.all(np.abs(a[sel]) + np.abs(b[sel]) < 0.1))


def jac_rhs(r, x, g, w, lam, spl_b):
    """My linearization of the May system about (a, b): the a-perturbation
    equation is unchanged in form; the b equation gains -12 g b^2 xb.
    x packs BOTH vanishing-IC solutions: (xa1, xa1', xb1, xb1', xa2, ...)."""
    bv = spl_b(r)
    out = np.empty(8)
    for j in (0, 1):
        xa, xap, xb, xbp = x[4 * j: 4 * j + 4]
        out[4 * j] = xap
        out[4 * j + 1] = xb - w**2 * xa + xa / r**2 - xap / r
        out[4 * j + 2] = xbp
        out[4 * j + 3] = (xa - lam * xb - 12 * g * bv**2 * xb
                          - w**2 * xb + xb / r**2 - xbp / r)
    return out


def gf_conj_points(rr, b, g, w, lam):
    """Proper GF determinant: solutions vanishing at r0 with derivative ICs
    = identity; count sign changes of det after it leaves the trivial zero."""
    spl_b = CubicSpline(rr, b)
    ev = np.linspace(R0, RTEST, 1600)
    ic = np.zeros(8)
    ic[1] = 1.0   # solution 1: xa'(r0) = 1
    ic[7] = 1.0   # solution 2: xb'(r0) = 1
    sol = solve_ivp(jac_rhs, (R0, RTEST), ic, t_eval=ev,
                    args=(g, w, lam, spl_b), method="LSODA",
                    rtol=1e-9, atol=1e-11)
    if not sol.success or sol.t.shape[0] != ev.shape[0]:
        return None
    det = sol.y[0] * sol.y[6] - sol.y[4] * sol.y[2]
    # skip the (r - r0)^2 departure zone near r0
    start = np.searchsorted(ev, R0 + 0.1)
    d = det[start:]
    s = np.sign(d[np.abs(d) > 1e-14])
    return int(np.sum(s[:-1] * s[1:] < 0))


def assess(g, w, lam, A0, B0, rtol=1e-8, atol=1e-10):
    row = dict(g=g, w=w, lam=lam, A0=round(A0, 4), B0=round(B0, 4))
    out = integrate(g, w, lam, A0, B0, rtol, atol)
    if out is None:
        row["diverged"] = True
        return row
    rr, y = out
    a, b = y[0], y[2]
    row["diverged"] = False
    row["crit_i"] = crit_i(rr, a, b)
    row["crit_ii"] = bool(my_nodes(rr, a) <= 4 and my_nodes(rr, b) <= 4)
    if row["crit_i"] and row["crit_ii"]:
        cp = gf_conj_points(rr, b, g, w, lam)
        row["gf_conj_points"] = cp
        row["crit_iii"] = bool(cp == 0) if cp is not None else None
        row["passes_all"] = bool(row["crit_iii"]) \
            if cp is not None else None
    else:
        row["passes_all"] = False
    return row


def main():
    t0 = time.time()
    combos = []
    # block 1: full small-amplitude slice (contains all 26 primary
    # near-passers)
    for g in (0.5, 2.0, 3.5, 5.0):
        for w in (0.3, 0.7, 1.1, 1.5):
            for lam in (0.0, 0.5, 1.0, 1.5, 2.0):
                combos.append((g, w, lam, 0.1, 0.1))
    # block 2: below-threshold block with larger amplitudes
    for w in (0.3, 0.7):
        for lam in (0.0, 1.0, 2.0):
            for A0, B0 in itertools.product((0.4, 1.0), repeat=2):
                for g in (0.5, 3.5):
                    combos.append((g, w, lam, A0, B0))
    # block 3: the printed representative + 8 random 5% perturbations
    rep = (0.5, 1.5, 0.0, 0.1, 0.1)
    combos.append(rep)
    rng = np.random.default_rng(20260721)
    for _ in range(8):
        f = 1 + 0.05 * rng.uniform(-1, 1, size=5)
        combos.append((rep[0] * f[0], rep[1] * f[1],
                       max(rep[2], 1e-6) * f[2], rep[3] * f[3],
                       rep[4] * f[4]))
    combos = list(dict.fromkeys(combos))

    rows = [assess(*cbo) for cbo in combos]
    # tolerance robustness: rerun every crit-i passer at tighter tolerance
    tighter = [assess(r["g"], r["w"], r["lam"], r["A0"], r["B0"],
                      rtol=1e-10, atol=1e-12)
               for r in rows if r.get("crit_i")]

    # scipy-free RK4 cross-check on the printed representative combo
    rk4 = rk4_integrate(*rep)
    rk4_out = None
    if rk4 is not None:
        rr, y = rk4
        a, b = y[0], y[2]
        rk4_out = {
            "crit_i": crit_i(rr, a, b),
            "nodes_a": my_nodes(rr, a), "nodes_b": my_nodes(rr, b),
            "tail_abs_max_r_ge_8": float(np.max(
                np.abs(a[rr >= 8]) + np.abs(b[rr >= 8]))),
            "gf_conj_points_on_rk4_background": gf_conj_points(
                rr, b, rep[0], rep[1], rep[2]),
        }

    ok = [r for r in rows if not r["diverged"]]
    crit_i_rows = [r for r in ok if r.get("crit_i")]
    passers = [r for r in ok if r.get("passes_all")]
    summary = {
        "n_combos": len(rows),
        "diverged": sum(r["diverged"] for r in rows),
        "crit_i_pass": len(crit_i_rows),
        "full_passers": len(passers),
        "passer_list": passers,
        "crit_i_conj_point_range": [
            min(r["gf_conj_points"] for r in crit_i_rows
                if r.get("gf_conj_points") is not None),
            max(r["gf_conj_points"] for r in crit_i_rows
                if r.get("gf_conj_points") is not None)]
        if crit_i_rows else None,
        "crit_i_combos": [
            {k: r[k] for k in ("g", "w", "lam", "A0", "B0",
                               "gf_conj_points")}
            for r in crit_i_rows],
        "tighter_tol_full_passers": sum(
            bool(r.get("passes_all")) for r in tighter),
        "tighter_tol_conj_point_range": [
            min(r["gf_conj_points"] for r in tighter
                if r.get("gf_conj_points") is not None),
            max(r["gf_conj_points"] for r in tighter
                if r.get("gf_conj_points") is not None)]
        if tighter else None,
        "representative_rk4_crosscheck": rk4_out,
    }
    # published-count grep (verbatim, with line numbers)
    txt = (DATA / "m6_4_record_may_fulltext.txt").read_text().splitlines()
    summary["printed_62_lines"] = [
        i + 1 for i, ln in enumerate(txt) if "62" in ln and "famil" in ln
        or "62 yield stable" in ln or "62 stable chaoiton" in ln]
    summary["printed_42_lines"] = [
        i + 1 for i, ln in enumerate(txt) if "42 stable chaoiton" in ln]
    summary["printed_62_quotes"] = [txt[i - 1].strip()
                                    for i in summary["printed_62_lines"]]
    summary["printed_42_quotes"] = [txt[i - 1].strip()
                                    for i in summary["printed_42_lines"]]
    out = {"task": "M6.4 audit A4: May census independent subsample",
           "summary": summary,
           "rows": rows,
           "runtime_s": round(time.time() - t0, 1)}
    (DATA / "m6_4_audit_may_census.json").write_text(
        json.dumps(out, indent=2, default=str))
    print(json.dumps(summary, indent=2, default=str))
    print(f"runtime {out['runtime_s']} s")


if __name__ == "__main__":
    main()
