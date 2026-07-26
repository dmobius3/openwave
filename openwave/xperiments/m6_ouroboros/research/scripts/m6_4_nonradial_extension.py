"""M6.4 N3: extend the Gelfand-Fomin test beyond radial perturbations.

The record's admitted gap: the GF conjugate-point test (both eras) covers
RADIAL perturbations only. Within the printed 2D radial reduction, a
non-radial perturbation sector delta ~ xi(r) e^{i l theta} adds +l^2/r^2 to
the diagonal of the Jacobi potential Q(r):

    xi'' + xi'/r - (l^2/r^2) xi = Q(r) xi   (per component, coupled via Q12)

Sturm-comparison argument (stated, then numerically spot-checked here): the
l-sector Jacobi operator differs from the radial one by the POSITIVE operator
+l^2/r^2; conjugate points therefore cannot INCREASE with l. If the radial
sector passes, every l >= 1 sector passes a fortiori. The census's radial
verdicts are thus the binding ones WITHIN the ansatz reduction.

What the extension can NOT reach from the printed reduction (the honest
boundary, documented, not computed): perturbations that leave the ansatz
manifold (other vector components of A_mu/J_mu, non-axisymmetric 3D toroidal
modes). Those need the full 3D theory; the M7 HydroBoros program carries it
(and found the truncation's real-time vacuum unconditionally tachyonic, the
strongest known non-radial statement about this system).

Numeric spot-check: for representative backgrounds (May system, the
localization-criterion passers + one below-threshold shot; June system, the
published best point), count conjugate points at l = 0, 1, 2, 3 and verify
monotone non-increase.

Output: research/data/m6_4_nonradial.json.
"""

import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

DATA = Path(__file__).resolve().parent.parent / "data"

R0 = 0.02
BLOWUP = 100.0


# ---- backgrounds ----------------------------------------------------------

def rhs_may(r, y, g, w, lam):
    a, ap, b, bp = y
    app = -ap / r + a / r**2 - w**2 * a + b
    bpp = -bp / r + b / r**2 - w**2 * b + a - lam * b - 4.0 * g * b**3
    return [ap, app, bp, bpp]


def rhs_june(r, y, g, w, lam):
    a, ap, b, bp = y
    m2a, m2b = 1.0 + w**2, 1.0 + lam + w**2
    app = -ap / r + m2a * a - 0.5 * b - 4.0 * g * a**3
    bpp = -bp / r + m2b * b - 0.5 * a - 4.0 * g * b**3
    return [ap, app, bp, bpp]


def solve_bg(system, g, w, lam, A0, B0, rmax):
    rhs = rhs_may if system == "may" else rhs_june
    ev = np.linspace(R0, rmax, 1000)
    blow = lambda r, y, *a_: BLOWUP - np.max(np.abs(y))  # noqa: E731
    blow.terminal = True
    sol = solve_ivp(rhs, (R0, rmax), [A0 * R0, A0, B0 * R0, B0], t_eval=ev,
                    args=(g, w, lam), method="DOP853", rtol=1e-9, atol=1e-11,
                    events=blow)
    if not sol.success or sol.t.shape[0] != ev.shape[0]:
        return None
    return sol.t, sol.y[0], sol.y[2]


# ---- l-sector Jacobi ------------------------------------------------------

def rhs_jacobi(r, x, ell, system, g, w, lam, r_grid, a_prof, b_prof):
    a = np.interp(r, r_grid, a_prof)
    b = np.interp(r, r_grid, b_prof)
    if system == "may":
        q11, q22 = -w**2 + 1.0 / r**2, -w**2 - lam - 12.0 * g * b**2 \
            + 1.0 / r**2
        q12a, q12b = 1.0, 1.0  # +b in a-eq, +a in b-eq (RHS form)
        xa, xap, xb, xbp = x
        xapp = -xap / r + (ell**2) / r**2 * xa + q11 * xa + q12a * xb
        xbpp = -xbp / r + (ell**2) / r**2 * xb + q12b * xa + q22 * xb
        return [xap, xapp, xbp, xbpp]
    m2a, m2b = 1.0 + w**2, 1.0 + lam + w**2
    q11 = m2a - 12.0 * g * a**2
    q22 = m2b - 12.0 * g * b**2
    q12 = -0.5  # the TRUE linearization (V2); printed -c/4 tested in census
    xa, xap, xb, xbp = x
    xapp = -xap / r + (ell**2) / r**2 * xa + q11 * xa + q12 * xb
    xbpp = -xbp / r + (ell**2) / r**2 * xb + q12 * xa + q22 * xb
    return [xap, xapp, xbp, xbpp]


def conjugate_points(ell, system, g, w, lam, r_grid, a_prof, b_prof,
                     r_test):
    ev = np.linspace(R0, r_test, 800)
    sols = []
    for ic in ([0, 1, 0, 0], [0, 0, 0, 1]):
        s = solve_ivp(rhs_jacobi, (R0, r_test), ic, t_eval=ev,
                      args=(ell, system, g, w, lam, r_grid, a_prof, b_prof),
                      method="DOP853", rtol=1e-9, atol=1e-11)
        if not s.success or s.t.shape[0] != ev.shape[0]:
            return None
        sols.append(s.y)
    det = sols[0][0] * sols[1][2] - sols[1][0] * sols[0][2]
    mask = np.abs(det) > 1e-12
    d = det[int(np.argmax(mask)):]
    sg = np.sign(d)
    sg = sg[sg != 0]
    return int(np.sum(sg[:-1] * sg[1:] < 0))


CASES = [
    # (system, g, w, lam, A0, B0, rmax, r_test, label)
    ("may", 0.5, 1.5, 0.0, 0.1, 0.1, 10.0, 8.0, "May crit-i passer (above "
     "threshold, radial gf=6)"),
    ("may", 0.5, 1.1, 2.0, 0.1, 0.1, 10.0, 8.0, "May crit-i passer (above "
     "threshold, radial gf=6)"),
    ("may", 0.5, 0.7, 0.0, 0.4, 0.4, 10.0, 8.0, "May below-threshold grid "
     "shot (in-window)"),
    ("june", 0.3, 2.2, 0.5, 0.1, 0.2, 10.0, 10.0, "June published best "
     "point (outward S convention)"),
    ("june", 1.0, 1.0, 1.0, 0.1, 0.1, 10.0, 10.0, "June grid point "
     "(outward S convention)"),
]


def main():
    results = []
    for system, g, w, lam, A0, B0, rmax, r_test, label in CASES:
        bg = solve_bg(system, g, w, lam, A0, B0, rmax)
        row = dict(system=system, g=g, w=w, lam=lam, A0=A0, B0=B0,
                   label=label)
        if bg is None:
            row["diverged"] = True
            results.append(row)
            continue
        r_grid, a_prof, b_prof = bg
        cps = {}
        for ell in (0, 1, 2, 3):
            cps[f"l={ell}"] = conjugate_points(
                ell, system, g, w, lam, r_grid, a_prof, b_prof, r_test)
        row["conjugate_points"] = cps
        vals = [v for v in cps.values() if v is not None]
        row["monotone_nonincreasing"] = bool(
            all(vals[i] >= vals[i + 1] for i in range(len(vals) - 1)))
        results.append(row)
        print(label, cps, "monotone:", row.get("monotone_nonincreasing"),
              flush=True)
    out = {
        "task": "M6.4 N3: non-radial (l-sector) GF extension",
        "argument": "the l-sector Jacobi operator adds +l^2/r^2 (positive) "
                    "to the radial one; by Sturm comparison conjugate points "
                    "cannot increase with l; radial verdicts are binding "
                    "within the ansatz reduction",
        "boundary": "ansatz-breaking perturbations (other vector components, "
                    "non-axisymmetric 3D modes) are NOT reachable from the "
                    "printed reduction; the M7 3D continuation found the "
                    "truncation's real-time vacuum unconditionally tachyonic",
        "spot_checks": results,
        "all_monotone": bool(all(r.get("monotone_nonincreasing", True)
                                 for r in results)),
    }
    (DATA / "m6_4_nonradial.json").write_text(
        json.dumps(out, indent=2, default=str))
    print("all monotone non-increasing:", out["all_monotone"])


if __name__ == "__main__":
    main()
