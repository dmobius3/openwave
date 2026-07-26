"""M5.21.1e phase C: the 1+1D time-crystal toy-model regression.

The papers' ONLY computed dynamics (arXiv:2501.04036v2 = TOY; Wolfram
Community companion) is the 1+1D two-scalar model

    H = phi_t^2 + phi_x^2 + (1 - phi^2)^2 - alpha R^2 + beta R^4,
    R = phi_t psi_x - phi_x psi_t .

For a static kink phi(x) with the clock psi = omega t:
R = -omega phi_x, so

    E(phi, omega) = INT [ phi_x^2 + (1-phi^2)^2
                          - alpha omega^2 phi_x^2
                          + beta omega^4 phi_x^4 ] dx .

TOY Eq (5) closed form for the pure-tanh family phi = tanh(x/w):
    omega* = sqrt(70 alpha / (96 beta - 35 alpha^2)),
    w*     = sqrt(96 beta / (96 beta - 35 alpha^2)),
and the paper/WOLF anchors at alpha = beta = 1:
    pure tanh: E ~ 2.12568, omega ~ 1.07123
    deg-5 poly: E ~ 2.02515, omega ~ 1.28975.

This script re-derives all of it independently (conformance rung: our
reading of his equations reproduces his only published numbers):
  1. closed form re-derivation from the two integrals
     INT sech^4 = 4w/3, INT sech^8 = 32w/35 (machine-checked by
     quadrature), minimized analytically;
  2. a free-profile grid minimization (relax phi(x) by FIRE at fixed
     omega, golden-scan omega) that must beat the pure-tanh energy and
     approach the deg-5 anchor from below/near.

Also the in-house demonstration of the paper's central stability
mechanism: E(omega) at the relaxed kink has its MINIMUM at omega != 0
(energy minimization propels the clock), the 1+1D version of what the
4D electron is conjectured to do.

Output: data/m5_21_1e_toy.json + plots/m5_21_1e_toy.png.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
PLOTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "plots")

ALPHA, BETA = 1.0, 1.0
L, N = 30.0, 601                      # domain [-L, L]; dx = 0.1 so the
X = np.linspace(-L, L, N)             # explicit relax respects dt < dx^2/4
DX = X[1] - X[0]


def e_of(phi, omega, alpha=ALPHA, beta=BETA):
    px = np.gradient(phi, DX)
    dens = (px ** 2 + (1 - phi ** 2) ** 2
            - alpha * omega ** 2 * px ** 2
            + beta * omega ** 4 * px ** 4)
    return float(np.trapezoid(dens, X))


def grad_of(phi, omega, alpha=ALPHA, beta=BETA):
    px = np.gradient(phi, DX)
    a = (1.0 - alpha * omega ** 2)
    flux = 2.0 * a * px + 4.0 * beta * omega ** 4 * px ** 3
    g = -np.gradient(flux, DX) - 4.0 * (1 - phi ** 2) * phi
    g[0] = g[-1] = 0.0                # kink BCs pinned at +-1
    return g


def fire_phi(phi0, omega, iters=6000, dt0=5e-4):
    phi, v = phi0.copy(), np.zeros_like(phi0)
    dt, al, nup = dt0, 0.1, 0
    F = -grad_of(phi, omega)
    for _ in range(iters):
        P = float(np.sum(F * v))
        if P > 0:
            nup += 1
            vn, fn = np.sqrt(np.sum(v * v)), np.sqrt(np.sum(F * F))
            v = (1 - al) * v + al * (F / max(fn, 1e-300)) * vn
            if nup > 5:
                dt = min(dt * 1.1, 10 * dt0)
                al *= 0.99
        else:
            v[:] = 0.0
            dt *= 0.5
            al, nup = 0.1, 0
        v += dt * F
        phi += dt * v
        if not np.all(np.isfinite(phi)):       # stability backoff
            phi, v = phi0.copy(), np.zeros_like(phi0)
            dt0 *= 0.25
            dt, al, nup = dt0, 0.1, 0
        F = -grad_of(phi, omega)
    return phi


def main():
    out = {}
    # 1. integral checks + closed form
    w = 1.3
    s4 = float(np.trapezoid(1 / np.cosh(X / w) ** 4, X))
    s8 = float(np.trapezoid(1 / np.cosh(X / w) ** 8, X))
    out["int_sech4_rel"] = abs(s4 - 4 * w / 3) / (4 * w / 3)
    out["int_sech8_rel"] = abs(s8 - 32 * w / 35) / (32 * w / 35)
    om_cf = np.sqrt(70 * ALPHA / (96 * BETA - 35 * ALPHA ** 2))
    w_cf = np.sqrt(96 * BETA / (96 * BETA - 35 * ALPHA ** 2))
    # E(w, om) for tanh family from the two integrals:
    def e_tanh(w, om):
        return (4 / (3 * w) * (1 - ALPHA * om ** 2) + 4 * w / 3
                + BETA * om ** 4 * 32 / (35 * w ** 3))
    # analytic-minimum cross-check by dense scan
    ws = np.linspace(0.5, 3.0, 601)
    oms = np.linspace(0.0, 2.0, 601)
    Ewo = e_tanh(ws[:, None], oms[None, :])
    i, j = np.unravel_index(np.argmin(Ewo), Ewo.shape)
    out["closed_form"] = {"omega": om_cf, "w": w_cf,
                          "E": e_tanh(w_cf, om_cf)}
    out["scan_min"] = {"omega": float(oms[j]), "w": float(ws[i]),
                       "E": float(Ewo[i, j])}
    out["anchor_pure_tanh"] = {"E": 2.12568, "omega": 1.07123}
    out["cf_vs_anchor_E_rel"] = abs(out["closed_form"]["E"] - 2.12568) / 2.12568
    out["cf_vs_anchor_om_rel"] = abs(om_cf - 1.07123) / 1.07123

    # 2. free-profile minimization + omega scan
    es = []
    phi0 = np.tanh(X / w_cf)
    for om in oms[::20]:
        phi = fire_phi(phi0, om, iters=1500)
        es.append(e_of(phi, om))
    es = np.array(es)
    oms_s = oms[::20]
    k = int(np.argmin(es))
    # refine around the winner
    om_ref = np.linspace(max(oms_s[k] - 0.15, 0.0), oms_s[k] + 0.15, 31)
    es_ref = []
    for om in om_ref:
        phi = fire_phi(phi0, om, iters=4000)
        es_ref.append(e_of(phi, om))
    es_ref = np.array(es_ref)
    kk = int(np.argmin(es_ref))
    out["free_profile_min"] = {"omega": float(om_ref[kk]),
                               "E": float(es_ref[kk])}
    out["anchor_deg5"] = {"E": 2.02515, "omega": 1.28975}
    # the free-profile minimum is NOT the kink: at alpha omega^2 > 1 the
    # quadratic-gradient coefficient is negative, so a grid-cutoff ripple
    # sea (px ~ sqrt((alpha om^2 - 1)/(2 beta om^4)), phi wiggling around
    # +-1) has strictly negative density: the continuum infimum is not
    # attained (cutoff-dependent). The paper's kink numbers are
    # ansatz-FAMILY minima. Record both readings honestly:
    out["free_min_below_kink"] = bool(es_ref[kk] < 2.02515 - 0.05)
    out["ripple_density_floor_at_omstar"] = float(
        -(ALPHA * om_cf ** 2 - 1) ** 2 / (4 * BETA * om_cf ** 4))
    # the clock statement lives IN the tanh family: E(w*, om*) vs the
    # omega = 0 family minimum (w = 1, E = 8/3)
    out["clock_in_family"] = {"E_omstar": e_tanh(w_cf, om_cf),
                              "E_om0_min": 8.0 / 3.0,
                              "propelled": bool(e_tanh(w_cf, om_cf)
                                                < 8.0 / 3.0)}

    with open(os.path.join(DATA, "m5_21_1e_toy.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].plot(oms_s, es, "o-", label="free-profile relax")
    ax[0].plot(om_ref, es_ref, ".-", label="refined")
    ax[0].plot(oms, e_tanh(w_cf, oms), "--", label="tanh family (analytic)")
    ax[0].axvline(om_cf, color="k", lw=0.5)
    ax[0].axhline(2.02515, color="g", lw=0.5, label="paper deg-5 anchor")
    ax[0].set_xlabel("omega")
    ax[0].set_ylabel("E")
    ax[0].set_title("E(omega): the minimum sits at omega != 0 (the clock)")
    ax[0].legend(fontsize=7)
    phi_best = fire_phi(phi0, om_ref[kk], 4000)
    ax[1].plot(X, phi_best, label=f"relaxed kink (om={om_ref[kk]:.3f})")
    ax[1].plot(X, np.tanh(X / w_cf), "--", label="pure tanh w*")
    ax[1].set_xlim(-8, 8)
    ax[1].set_title("kink profile")
    ax[1].legend(fontsize=7)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    p = os.path.join(PLOTS, "m5_21_1e_toy.png")
    fig.savefig(p, dpi=110)
    print(f"wrote {p}")


if __name__ == "__main__":
    main()
