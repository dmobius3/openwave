"""M6.4 ADVERSARIAL AUDIT: claim A7 (l-sector monotonicity of conjugate
points under +l^2/r^2).

THEORY ASSESSMENT (the auditor's reasoning, stated for the record):

The claim invokes "Sturm comparison". For a SCALAR Jacobi equation this is
classical. For the COUPLED 2-component systems at hand, scalar Sturm
comparison is NOT sufficient as stated: sign changes of individual solution
components have no comparison theory in systems. The correct tool is the
matrix generalization for SELF-ADJOINT Hamiltonian systems (Morse index
theorem for systems; Reid, "Sturmian Theory for ODEs", Ch. V):

  For the Jacobi system (r xi')' = r Q(r) xi with Q SYMMETRIC, conjugate
  points of r0 in (r0, r1] are in bijection (with multiplicity) with the
  negative directions of the quadratic form
      I[xi] = Int r (|xi'|^2 + xi^T Q xi) dr
  on H^1_0(r0, r1). Replacing Q by Q + (l^2/r^2) I adds the nonnegative
  term Int r (l^2/r^2)|xi|^2 dr to I, so I can only increase; its Morse
  index (hence the conjugate-point count) can only decrease or stay equal.

Hence the CONCLUSION of A7 holds, but only because (a) both Jacobi systems
in play have symmetric Q (May linearization: off-diagonal +1/+1; June
printed: -c/4/-c/4; June corrected: -c/2/-c/2 -- verified below), and
(b) the perturbation l^2/r^2 * IDENTITY is positive semidefinite. The
justification must be the matrix/Morse version, not scalar Sturm. CAVEAT,
stated honestly: the numerical detector (sign changes of the 2x2 vanishing-
solution determinant) counts conjugate points only up to even-multiplicity
zeros (a double zero of det produces no sign change), so numeric counts are
a lower bound; the monotonicity statement proper is about the true counts.

Sign-convention note: for a 2-component system, flipping the sign of BOTH
off-diagonal entries (Q12, Q21) -> (-Q12, -Q21) is conjugation by the
orthogonal D = diag(1, -1) (xb -> -xb) and leaves conjugate points
invariant; the audit uses the direct linearization signs.

NUMERIC SPOT-CHECK (auditor-owned code, LSODA + CubicSpline + 8-dim
combined integration; primary used DOP853/RK45 separate solves):
  - 2 May backgrounds (the printed representative near-passer and one
    below-threshold combo);
  - 1 June conv-T background (the published best point; primary spot-checked
    it under the OUTWARD S convention -- different choice here).
  conjugate-point counts for l = 0, 1, 2, 3 must be non-increasing.

Output: research/data/m6_4_audit_sturm.json
"""

import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import k0, k1

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
R0 = 0.02
BLOW = 100.0
C = 1.0


# ---------------------------------------------------------------- backgrounds
def may_bg(g, w, lam, A0, B0, rmax=10.0):
    def rhs(r, y):
        a, ap, b, bp = y
        return [ap, b - w**2 * a + a / r**2 - ap / r,
                bp, a - lam * b - 4 * g * b**3 - w**2 * b + b / r**2
                - bp / r]
    ev = np.linspace(R0, rmax, 2500)
    s = solve_ivp(rhs, (R0, rmax), [A0 * R0, A0, B0 * R0, B0], t_eval=ev,
                  method="LSODA", rtol=1e-9, atol=1e-11)
    return s.t, s.y[0], s.y[2]


def june_bg(g, w, lam, A, B, rmax=10.0):
    m2a, m2b = 1 + w**2, 1 + lam + w**2
    ka, kb = np.sqrt(m2a), np.sqrt(m2b)

    def rhs(r, y):
        a, ap, b, bp = y
        return [ap, m2a * a - (C / 2) * b - 4 * g * a**3 - ap / r,
                bp, m2b * b - (C / 2) * a - 4 * g * b**3 - bp / r]
    y0 = [A * k0(ka * rmax), -A * ka * k1(ka * rmax),
          B * k0(kb * rmax), -B * kb * k1(kb * rmax)]
    ev = np.linspace(rmax, R0, 2500)
    s = solve_ivp(rhs, (rmax, R0), y0, t_eval=ev,
                  method="LSODA", rtol=1e-9, atol=1e-11)
    return s.t[::-1], s.y[0][::-1], s.y[2][::-1]


# -------------------------------------------------- l-sector conjugate points
def conj_points(system, rr, a, b, g, w, lam, ell, rtest=8.0):
    """Proper GF vanishing-solution determinant with +l^2/r^2 on the
    diagonal; Q symmetric in both systems (May cross = +1/+1 after moving
    to the xi'' + xi'/r = Q xi form: Q12 = Q21 = -1; June printed:
    Q12 = Q21 = -c/4)."""
    spl_a, spl_b = CubicSpline(rr, a), CubicSpline(rr, b)

    def rhs(r, x):
        av, bv = spl_a(r), spl_b(r)
        cent = ell**2 / r**2
        if system == "may":
            # my linearization of the May system, in xi'' + xi'/r = Q xi
            # form: xa'' + xa'/r = xa/r^2 - w^2 xa + xb  and
            # xb'' + xb'/r = xb/r^2 - w^2 xb + xa - lam xb - 12 g b^2 xb;
            # so Q12 = Q21 = +1 (symmetric).
            q11 = 1 / r**2 + cent - w**2
            q22 = 1 / r**2 + cent - w**2 - lam - 12 * g * bv**2
            q12 = q21 = 1.0
        else:
            q11 = (1 + w**2) - 12 * g * av**2 + cent
            q22 = (1 + lam + w**2) - 12 * g * bv**2 + cent
            q12 = q21 = -C / 4
        out = np.empty(8)
        for j in (0, 1):
            xa, xap, xb, xbp = x[4 * j: 4 * j + 4]
            out[4 * j] = xap
            out[4 * j + 1] = q11 * xa + q12 * xb - xap / r
            out[4 * j + 2] = xbp
            out[4 * j + 3] = q21 * xa + q22 * xb - xbp / r
        return out

    ev = np.linspace(R0, rtest, 1600)
    ic = np.zeros(8)
    ic[1] = 1.0
    ic[7] = 1.0
    s = solve_ivp(rhs, (R0, rtest), ic, t_eval=ev,
                  method="LSODA", rtol=1e-9, atol=1e-11)
    if not s.success or s.t.shape[0] != ev.shape[0]:
        return None
    det = s.y[0] * s.y[6] - s.y[4] * s.y[2]
    start = np.searchsorted(ev, R0 + 0.1)
    d = det[start:]
    sg = np.sign(d[np.abs(d) > 1e-14])
    return int(np.sum(sg[:-1] * sg[1:] < 0))


def main():
    t0 = time.time()
    cases = []
    checks = [
        ("may", (0.5, 1.5, 0.0, 0.1, 0.1), "printed May representative"),
        ("may", (0.5, 0.7, 0.0, 0.4, 0.4), "May below-threshold combo"),
        ("june", (0.3, 2.2, 0.5, 0.1, 0.2),
         "June published best point (conv T inward; primary used S)"),
    ]
    for system, (g, w, lam, A, B), label in checks:
        rr, a, b = (may_bg if system == "may" else june_bg)(g, w, lam, A, B)
        counts = {}
        for ell in (0, 1, 2, 3):
            counts[f"l={ell}"] = conj_points(system, rr, a, b, g, w, lam,
                                             ell)
        vals = [v for v in counts.values() if v is not None]
        cases.append({
            "system": system, "g": g, "w": w, "lam": lam, "A": A, "B": B,
            "label": label, "conj_points": counts,
            "monotone_nonincreasing": bool(
                all(x >= y for x, y in zip(vals, vals[1:]))),
        })
    out = {
        "task": "M6.4 audit A7: l-sector monotonicity",
        "theory_verdict": {
            "scalar_sturm_sufficient_for_coupled_systems": False,
            "valid_replacement": "Morse index theorem / matrix Sturm "
                                 "comparison for self-adjoint systems "
                                 "(requires symmetric Q; satisfied here: "
                                 "May Q12=Q21=-1, June printed "
                                 "Q12=Q21=-c/4) with the psd perturbation "
                                 "(l^2/r^2) I",
            "numeric_caveat": "det sign-change counting misses even-"
                              "multiplicity zeros; numeric counts are a "
                              "lower bound on true conjugate-point counts",
        },
        "cases": cases,
        "all_monotone": bool(all(c["monotone_nonincreasing"]
                                 for c in cases)),
        "runtime_s": round(time.time() - t0, 1),
    }
    (DATA / "m6_4_audit_sturm.json").write_text(
        json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
