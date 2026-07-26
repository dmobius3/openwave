"""ADVERSARIAL AUDIT M6.2 — independent numerics for B5, B6, B8 + overturn scans.

Deliberately different from the task's pipeline:
 - integrator: LSODA at rtol 1e-11 / atol 1e-13 (task: RK45 1e-9/1e-11 max_step 0.02)
 - quadrature: numpy.trapz on the native dense grid (task: midpoint rectangles)
 - grid: 4001 points (task: 800)
Cross-check: also run the task's exact config to verify their certificate numbers.
"""
import json
import numpy as np
from scipy.integrate import solve_ivp

G, LAM, W, A0, B0 = 1.0, 1.0, 1.0, 0.1, 0.1
R_INNER = 0.02


def rhs(r, y, g=G, w=W, lam=LAM):
    a, da, b, db = y
    d2a = -da / r + a / r**2 - w**2 * a + b
    d2b = -db / r + b / r**2 - w**2 * b + a - lam * b - 4 * g * b**3
    return [da, d2a, db, d2b]


def solve(r_max, method="LSODA", rtol=1e-11, atol=1e-13, n=4001, max_step=np.inf,
          A0_=A0, B0_=B0):
    r_eval = np.linspace(R_INNER, r_max, n)
    y0 = [A0_ * R_INNER, A0_, B0_ * R_INNER, B0_]
    sol = solve_ivp(rhs, [R_INNER, r_max], y0, t_eval=r_eval, method=method,
                    rtol=rtol, atol=atol, max_step=max_step)
    assert sol.success
    return sol.t, sol.y[0], sol.y[2], sol.y[1], sol.y[3]


def observables(r, a, b, da, db):
    """trapezoid quadrature on the native grid, analytic derivatives from the
    solver state (not finite differences)."""
    tz = lambda f: float(np.trapz(f * r, r))
    H_p = tz(da**2 + db**2 + a**2 / r**2 + b**2 / r**2 + 4 * G * b**4)
    Q_p = tz(b**2)
    quad = (W**2 * a**2 + da**2 + a**2 / r**2 + W**2 * b**2 + db**2 + b**2 / r**2)
    H_d = 0.25 * tz(quad) - 0.5 * tz(a * b) + 0.375 * G * tz(b**4)
    Q_d = W * tz(a**2 + b**2)
    Q_J = W * tz(b**2)
    H_dp = 0.25 * tz(quad) + 0.375 * G * tz(b**4)  # printed phases: <J.A>=0
    sgn = int(np.sum(np.sign(b[:-1]) != np.sign(b[1:])))
    return dict(H_p=H_p, Q_p=Q_p, HQ_p=H_p / Q_p, H_d=H_d, Q_d=Q_d,
                HQ_d=H_d / Q_d, HQ_J=H_d / Q_J, HQ_dp_Qd=H_dp / Q_d,
                HQ_dp_QJ=H_dp / Q_J, tail=float(abs(a[-1]) + abs(b[-1])),
                beta_signchg=sgn)


print("=" * 72)
print("B5: r_max = 12, MY pipeline (LSODA 1e-11/1e-13, trapz, 4001 pts)")
obs = observables(*solve(12.0))
for k, v in obs.items():
    print(f"  {k:10s} = {v:.6f}" if isinstance(v, float) else f"  {k:10s} = {v}")

print("\nB5 cross-check with RK45 tight (1e-12/1e-14), no max_step cap:")
obs_rk = observables(*solve(12.0, method="RK45", rtol=1e-12, atol=1e-14))
print(f"  HQ_p = {obs_rk['HQ_p']:.6f}  HQ_d = {obs_rk['HQ_d']:.6f}")

print("\nB5 gap arithmetic (claimed 91.53 / 91.577 / 68.441):")
for label, val in [("HQ_d vs 1.6875", abs(obs['HQ_d'] - 1.6875) / 1.6875 * 100),
                   ("HQ_d vs 1.6969", abs(obs['HQ_d'] - 1.6969) / 1.6969 * 100),
                   ("HQ_J vs 1.6875", abs(obs['HQ_J'] - 1.6875) / 1.6875 * 100),
                   ("claimed JSON 0.142929 vs 1.6875",
                    abs(0.142929 - 1.6875) / 1.6875 * 100),
                   ("claimed JSON 0.532557 vs 1.6875",
                    abs(0.532557 - 1.6875) / 1.6875 * 100)]:
    print(f"  {label}: {val:.3f}%")

print("\nB5 replication of the TASK's exact pipeline (RK45 1e-9/1e-11, "
      "max_step 0.02, 800 pts, midpoint):")
r_eval = np.linspace(R_INNER, 12.0, 800)
y0 = [A0 * R_INNER, A0, B0 * R_INNER, B0]
s = solve_ivp(rhs, [R_INNER, 12.0], y0, t_eval=r_eval, method="RK45",
              rtol=1e-9, atol=1e-11, max_step=0.02)
r, a, b = s.t, s.y[0], s.y[2]
dr = np.diff(r); rm = 0.5 * (r[:-1] + r[1:])
am = 0.5 * (a[:-1] + a[1:]); bm = 0.5 * (b[:-1] + b[1:])
dam = np.diff(a) / dr; dbm = np.diff(b) / dr
I = lambda f: float(np.sum(f * rm * dr))
H_p_task = I(dam**2 + dbm**2 + am**2 / rm**2 + bm**2 / rm**2 + 4 * G * bm**4)
Q_p_task = I(bm**2)
print(f"  H_p = {H_p_task:.6f} (claimed 0.49555), Q_p = {Q_p_task:.6f} "
      f"(claimed 0.293404), H/Q = {H_p_task / Q_p_task:.6f} (claimed 1.688971)")

print("\n" + "=" * 72)
print("B6: R_phys closure: hbar_c * H_p / m_e")
hbar_c, m_e = 197.3269804, 0.511
print(f"  with task H_p 0.49555 : R_phys = {hbar_c * 0.49555 / m_e:.2f} fm "
      f"(v11 prints 191 fm)")
print(f"  with MY H_p {obs['H_p']:.5f}: R_phys = {hbar_c * obs['H_p'] / m_e:.2f} fm")
print(f"  inverse: H_code needed for exactly 191 fm = "
      f"{191.0 * m_e / hbar_c:.5f}")

print("\n" + "=" * 72)
print("B8: window dependence, MY pipeline, r_max = 12 / 24 / 48")
rows = {}
for rmx in (12.0, 24.0, 48.0):
    o = observables(*solve(rmx, n=int(400 * rmx) + 1))
    rows[rmx] = o
    print(f"  r_max={rmx:5.0f}: H_p={o['H_p']:.4f}  Q_p={o['Q_p']:.4f}  "
          f"HQ_p={o['HQ_p']:.4f}  HQ_d={o['HQ_d']:.4f}  tail={o['tail']:.4f}  "
          f"beta sign changes={o['beta_signchg']}")
hqs = [rows[k]["HQ_p"] for k in rows]
print(f"  HQ_p drift over windows: max-min = {max(hqs) - min(hqs):.4f} "
      f"({(max(hqs) - min(hqs)) / 1.689 * 100:.1f}% of 1.689)")

print("\nOVERTURN (ii): is the tail numerical artifact? LSODA 1e-11 vs RK45 1e-9 "
      "task settings at r_max=48:")
o_tight = rows[48.0]
r_eval = np.linspace(R_INNER, 48.0, 3200)
s2 = solve_ivp(rhs, [R_INNER, 48.0], y0, t_eval=r_eval, method="RK45",
               rtol=1e-9, atol=1e-11, max_step=0.02)
r2, a2, b2 = s2.t, s2.y[0], s2.y[2]
sgn2 = int(np.sum(np.sign(b2[:-1]) != np.sign(b2[1:])))
print(f"  task-settings RK45: beta(48) = {b2[-1]:+.6f}, sign changes {sgn2}")
i_t = np.searchsorted(r_eval, 47.999)
print(f"  tight LSODA:        beta(48) = {solve(48.0, n=int(400*48)+1)[2][-1]:+.6f}, "
      f"sign changes {o_tight['beta_signchg']}")

print("\nOVERTURN (ii)b: far-field linear analysis (analytic):")
M = np.array([[-W**2, 1.0], [1.0, -(W**2 + LAM)]])
ev = np.linalg.eigvals(M)
print(f"  lap v = M v far field, eigenvalues of M = {np.sort(ev)}")
print("  both NEGATIVE => both channels oscillatory Bessel (rho^-1/2 envelope):")
print("  NO exponentially decaying branch exists at w=1, lam=1 for ANY A0,B0;")
print("  integrals of f^2 rho drho grow ~linearly. Localization impossible.")
print(f"  condition for one decaying channel: w^4 + lam*w^2 < 1 -> w < "
      f"{np.sqrt((np.sqrt(LAM**2 + 4) - LAM) / 2):.4f} at lam=1")

print("\nOVERTURN (iii): coarse A0/B0 scan near 0.1 for a truly decaying solution")
print("  (envelope test: min over scan of |beta| max on rho in [40,48])")
best = None
for A0_ in (0.05, 0.08, 0.1, 0.12, 0.15):
    for B0_ in (0.05, 0.08, 0.1, 0.12, 0.15):
        r3, a3, b3, _, _ = solve(48.0, n=1921, A0_=A0_, B0_=B0_)
        mask = r3 > 40
        env = float(np.max(np.abs(b3[mask])) + np.max(np.abs(a3[mask])))
        sgnc = int(np.sum(np.sign(b3[:-1]) != np.sign(b3[1:])))
        if best is None or env < best[2]:
            best = (A0_, B0_, env, sgnc)
        print(f"  A0={A0_:.2f} B0={B0_:.2f}: far-envelope={env:.4f}, "
              f"beta sign changes={sgnc}")
print(f"  best far-envelope: A0={best[0]}, B0={best[1]}, env={best[2]:.4f} "
      f"(a true bound state would be ~0); sign changes {best[3]}")
