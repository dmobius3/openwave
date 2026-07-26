"""M6.4 ADVERSARIAL AUDIT: claim A8 (OQ3 non-existence in the May system).

Claim: in the MAY system, panels (g=0.5, lam=0) and (g=1, lam=1), no
localized solution exists for (A0, B0) in (0, 1.2]^2 at any in-window omega;
the 3 candidate zeros at omega = 0.9 are aliasing artifacts (fields stay
O(0.65-1.49) at r = 25 vs ~3e-4 expected for kappa = 0.436 decay).

Auditor-owned methods (all different from the primary's):
  1. CANDIDATE PROFILES: integrate the exact candidate amplitudes
     (0.8297241423223741, 0.14884440333386637) and
     (1.0190475618320765, 0.09022087372298898) at (g=0.5, lam=0, w=0.9)
     with LSODA at two tolerances AND a scipy-free fixed-step RK4;
     report |a|+|b| at r = 25 vs the K1(kappa r) localization expectation.
  2. EXACT-BESSEL CHANNEL AMPLITUDE (no fitting; primary used cos/sin
     least squares and an envelope max): in the linear far field the
     oscillatory channel u = v-.(a,b) satisfies the order-1 Bessel equation
     with argument k r, k = sqrt(w^2 - X_-); from (u(R), u'(R)) alone,
     solve the 2x2 system for (cJ, cY) in u = cJ J1(kR) + cY Y1(kR); the
     asymptotic envelope is amp = sqrt(cJ^2 + cY^2) * sqrt(2/(pi k)).
     Evaluate at R = 12, 18, 24: a genuine localized solution needs
     amp ~ 0 CONSISTENTLY across R; artifacts are R-inconsistent.
     Scan a 15x15 (A0, B0) grid over (0.02, 1.2]^2 at 3 omegas per panel
     (incl. w = 0.9 for panel 1 and the near-window-top w = 0.763 for
     panel 2) and record min over the grid.
  3. DIRECT BVP (solve_bvp, a solver family the primary never used): on
     [r0, 20], boundary conditions = regularity at r0 (a' - a/r = 0,
     b' - b/r = 0), pure-decay at R (v- component = 0 AND
     v+' + kappa v+ = 0), i.e. exactly the 4 conditions of a localized
     state, PLUS a fixed-amplitude normalization a'(r0) = A0* with omega
     as an UNKNOWN PARAMETER (solve_bvp(p=[w])): if a localized state
     existed at any nearby omega, this converges to it. Seeded at 3 omegas
     per panel with a decaying-bump initial guess. Convergence is then
     validated against the leftover 5th condition (v-'(R) = 0) and against
     criterion (1)-style re-integration; non-convergence or
     zero-collapse = no localized state found.

Output: research/data/m6_4_audit_oq3.json
"""

import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.special import j1, jvp, k1, kvp, y1, yvp

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
R0 = 0.02
BLOW = 100.0


def x_pm(lam):
    s = np.sqrt(lam**2 + 4.0)
    return (-lam + s) / 2.0, (-lam - s) / 2.0


def channel_vectors(lam):
    """Unit eigenvectors of [[0, 1], [1, -lam]] for X_+ (decay) and X_-
    (oscillation); (a, b) = va * v+ + vb * v-."""
    xp, xm = x_pm(lam)
    vp = np.array([1.0, xp])
    vm = np.array([1.0, xm])
    return vp / np.linalg.norm(vp), vm / np.linalg.norm(vm)


def rhs_may(r, y, g, w, lam):
    a, ap, b, bp = y
    return [ap, b - w**2 * a + a / r**2 - ap / r,
            bp, a - lam * b - 4 * g * b**3 - w**2 * b + b / r**2 - bp / r]


def shoot(g, w, lam, A0, B0, rfar, rtol=1e-10, atol=1e-12):
    ev = np.linspace(R0, rfar, 3000)
    hit = lambda r, y, *a_: BLOW - np.max(np.abs(y))  # noqa: E731
    hit.terminal = True
    s = solve_ivp(rhs_may, (R0, rfar), [A0 * R0, A0, B0 * R0, B0],
                  t_eval=ev, args=(g, w, lam), method="LSODA",
                  rtol=rtol, atol=atol, events=hit)
    if not s.success or s.t.shape[0] != ev.shape[0]:
        return None
    return s.t, s.y


def rk4_shoot(g, w, lam, A0, B0, rfar, h=2e-4):
    n = int((rfar - R0) / h)
    y = np.array([A0 * R0, A0, B0 * R0, B0])
    f = lambda r, y_: np.array(rhs_may(r, y_, g, w, lam))  # noqa: E731
    r = R0
    keep_r, keep_y = [r], [y.copy()]
    stride = max(1, n // 3000)
    for i in range(n):
        k1_ = f(r, y)
        k2_ = f(r + h / 2, y + h / 2 * k1_)
        k3_ = f(r + h / 2, y + h / 2 * k2_)
        k4_ = f(r + h, y + h * k3_)
        y = y + h / 6 * (k1_ + 2 * k2_ + 2 * k3_ + k4_)
        r += h
        if np.max(np.abs(y)) > BLOW:
            return None
        if (i + 1) % stride == 0:
            keep_r.append(r)
            keep_y.append(y.copy())
    return np.array(keep_r), np.array(keep_y).T


def bessel_osc_amplitude(rr, y, w, lam, R):
    """Exact-Bessel oscillatory-channel amplitude from point values at R."""
    _, xm = x_pm(lam)
    _, vm = channel_vectors(lam)
    k = np.sqrt(w**2 - xm)
    i = np.searchsorted(rr, R)
    if i >= rr.shape[0] - 1:
        return None
    u = vm[0] * y[0][i] + vm[1] * y[2][i]
    up = vm[0] * y[1][i] + vm[1] * y[3][i]
    x = k * rr[i]
    A_ = np.array([[j1(x), y1(x)], [k * jvp(1, x), k * yvp(1, x)]])
    cJ, cY = np.linalg.solve(A_, [u, up])
    return float(np.hypot(cJ, cY) * np.sqrt(2 / (np.pi * k)))


def candidate_profiles():
    out = []
    g, lam, w = 0.5, 0.0, 0.9
    xp, _ = x_pm(lam)
    kappa = np.sqrt(xp - w**2)
    for A0, B0 in ((0.8297241423223741, 0.14884440333386637),
                   (1.0190475618320765, 0.09022087372298898)):
        row = {"A0": A0, "B0": B0, "kappa_analytic": float(kappa)}
        sol = shoot(g, w, lam, A0, B0, 30.0)
        if sol is None:
            row["lsoda"] = "diverged"
        else:
            rr, y = sol
            mag = np.abs(y[0]) + np.abs(y[2])
            i5, i25 = np.searchsorted(rr, 5.0), np.searchsorted(rr, 25.0)
            row["lsoda"] = {
                "mag_r5": float(mag[i5]), "mag_r25": float(mag[i25]),
                "expected_r25_if_localized": float(
                    mag[i5] * (k1(kappa * 25) * 1) / k1(kappa * 5)),
                "bessel_amp_R12": bessel_osc_amplitude(rr, y, w, lam, 12),
                "bessel_amp_R18": bessel_osc_amplitude(rr, y, w, lam, 18),
                "bessel_amp_R24": bessel_osc_amplitude(rr, y, w, lam, 24),
            }
        rk = rk4_shoot(g, w, lam, A0, B0, 30.0)
        if rk is None:
            row["rk4"] = "diverged"
        else:
            rr, y = rk
            mag = np.abs(y[0]) + np.abs(y[2])
            i25 = np.searchsorted(rr, 25.0)
            row["rk4"] = {"mag_r25": float(mag[i25])}
        out.append(row)
    return out


def panel_scan(g, lam, omegas, n=15):
    xp, _ = x_pm(lam)
    res = []
    amps = np.linspace(0.02, 1.2, n)
    for w in omegas:
        best = None
        n_small_tail = 0
        n_diverged = 0
        for A0 in amps:
            for B0 in amps:
                sol = shoot(g, w, lam, A0, B0, 26.0, rtol=1e-9, atol=1e-11)
                if sol is None:
                    n_diverged += 1
                    continue
                rr, y = sol
                sel = rr >= 23.0
                t24 = float(np.max(np.abs(y[0][sel]) + np.abs(y[2][sel])))
                if t24 < 1e-3:
                    n_small_tail += 1
                amp = bessel_osc_amplitude(rr, y, w, lam, 18.0)
                if amp is not None and (best is None or amp < best[0]):
                    best = (amp, float(A0), float(B0), t24)
        res.append({
            "w": float(w),
            "kappa_pred": float(np.sqrt(max(xp - w**2, 0))),
            "grid": f"{n}x{n} over (0.02,1.2]^2",
            "diverged": n_diverged,
            "tail_mag_below_1e-3_at_r24": n_small_tail,
            "min_bessel_osc_amp_R18": best[0] if best else None,
            "argmin": best[1:3] if best else None,
            "tail_mag_at_argmin": best[3] if best else None,
        })
    return res


def bvp_hunt(g, lam, w_seeds, a0_norm=0.3, Rb=20.0):
    """solve_bvp with omega as unknown parameter; 5 BCs for 4 ODES + 1
    parameter. Returns per-seed outcome."""
    xp, xm = x_pm(lam)
    vp, vm = channel_vectors(lam)

    def fun(r, y, p):
        w = p[0]
        a, ap, b, bp = y
        return np.vstack([
            ap, b - w**2 * a + a / r**2 - ap / r,
            bp, a - lam * b - 4 * g * b**3 - w**2 * b + b / r**2 - bp / r])

    def bc(ya, yb, p):
        w = p[0]
        kappa = np.sqrt(max(xp - w**2, 1e-9))
        um = vm[0] * yb[0] + vm[1] * yb[2]
        upv = vp[0] * yb[0] + vp[1] * yb[2]
        upd = vp[0] * yb[1] + vp[1] * yb[3]
        return np.array([
            ya[1] - ya[0] / R0,          # a regular
            ya[3] - ya[2] / R0,          # b regular
            ya[1] - a0_norm,             # amplitude normalization a'(r0)
            um,                          # oscillatory channel value = 0
            upd + kappa * upv,           # decaying-channel log-derivative
        ])

    out = []
    for w0 in w_seeds:
        kappa0 = np.sqrt(max(xp - w0**2, 0.05))
        rmesh = np.linspace(R0, Rb, 400)
        bump = rmesh * np.exp(-kappa0 * rmesh)
        y0 = np.vstack([a0_norm * bump,
                        a0_norm * (1 - kappa0 * rmesh) * np.exp(
                            -kappa0 * rmesh),
                        0.5 * a0_norm * bump,
                        0.5 * a0_norm * (1 - kappa0 * rmesh) * np.exp(
                            -kappa0 * rmesh)])
        try:
            sol = solve_bvp(fun, bc, rmesh, y0, p=[w0], max_nodes=200000,
                            tol=1e-8)
        except Exception as e:  # noqa: BLE001
            out.append({"w_seed": w0, "status": f"exception: {e}"})
            continue
        row = {"w_seed": w0, "status": sol.status,
               "message": str(sol.message), "converged": bool(sol.success)}
        if sol.success:
            w_f = float(sol.p[0])
            row["w_found"] = w_f
            row["max_field"] = float(np.max(np.abs(sol.y[[0, 2]])))
            # leftover 5th condition: v- derivative at Rb must ALSO vanish
            umd = vm[0] * sol.y[1][-1] + vm[1] * sol.y[3][-1]
            row["leftover_vminus_deriv_at_Rb"] = float(umd)
            # re-integration check of the "solution" by shooting
            A0r = float(sol.y[1][0])
            B0r = float(sol.y[3][0])
            reint = shoot(g, w_f, lam, A0r, B0r, 25.0)
            if reint is None:
                row["reintegration"] = "diverged"
            else:
                rr, y = reint
                sel = rr >= 23.0
                row["reintegration_tail_mag_r24"] = float(
                    np.max(np.abs(y[0][sel]) + np.abs(y[2][sel])))
        out.append(row)
    return out


def main():
    t0 = time.time()
    out = {"task": "M6.4 audit A8: OQ3 existence (May system)"}
    out["part1_candidate_profiles"] = candidate_profiles()
    out["part2_panel1"] = panel_scan(0.5, 0.0, (0.45, 0.7, 0.9))
    out["part2_panel2"] = panel_scan(1.0, 1.0, (0.4, 0.6, 0.763))
    out["part3_bvp_panel1"] = bvp_hunt(0.5, 0.0, (0.45, 0.7, 0.9))
    out["part3_bvp_panel2"] = bvp_hunt(1.0, 1.0, (0.4, 0.6, 0.763))
    out["runtime_s"] = round(time.time() - t0, 1)
    (DATA / "m6_4_audit_oq3.json").write_text(
        json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
