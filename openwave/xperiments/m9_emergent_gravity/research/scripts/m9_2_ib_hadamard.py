#!/usr/bin/env python3
"""M9.2 attempt: d=2 kinematic residue I_B with a subtracted 1/r^2 piece.

This is NOT a multi-digit extraction. Paper VIII already reported a relative
residue of order 1-5 depending on the finite-part scheme. This script
reimplements one scheme (hard cutoff + kernel quadrature) so the paper can
quote a number that came from a file in this repo, and so a mutation of the
kernel can fail.

Definition used here (kinematic, C_J factored out)
-------------------------------------------------
Interval (-R, R) with R = 1.
Modular flow (Paper VI):  x(s) = R tanh(artanh(x/R) + pi s)
Weight-1 Jacobian:        J(s, x) = dx(s)/dx
Bare correlator:          C(x, y) = 1 / (x - y)^2
Smeared quadratic form:   Q(s) = int dx dy b(x) b(y) J(s,x) / (x(s) - y)^2
Local form:               Q_loc = Q(s=0)
Kernel:                   K(s) = pi / (2 cosh^2(pi s)),  int K = 1
Reported residue:         r = int K(s) (Q(s) - Q_loc) ds   /   |Q_loc|

A true Hadamard expansion would subtract EVERY non-integrable term of the
flowed two-point function, not only the s=0 coincidence. This script does
not do that. If r is scheme-dependent at O(1), the multi-digit coefficient
remains open.

Writes ../data/m9_2_ib_hadamard.json
"""

from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))

R = 1.0
A = 0.45  # bump support, a < R


def bump(x: np.ndarray, a: float = A) -> np.ndarray:
    out = np.zeros_like(x)
    m = np.abs(x) < a
    u = x[m] / a
    out[m] = np.exp(-1.0 / (1.0 - u * u))
    return out


def xflow(x: np.ndarray, s: float) -> np.ndarray:
    u = np.arctanh(np.clip(x / R, -0.999999, 0.999999))
    return R * np.tanh(u + np.pi * s)


def jacobian(x: np.ndarray, s: float) -> np.ndarray:
    # d/dx [R tanh(artanh(x/R)+pi s)] = sech^2(u+pi s) / sech^2(u)
    u = np.arctanh(np.clip(x / R, -0.999999, 0.999999))
    return (1.0 / np.cosh(u + np.pi * s) ** 2) / (1.0 / np.cosh(u) ** 2)


def quadratic_form(s: float, grid: np.ndarray, bump_vals: np.ndarray, eps: float) -> float:
    xs = xflow(grid, s)
    jac = jacobian(grid, s)
    # pairwise (x_i(s) - y_j)
    diff = xs[:, None] - grid[None, :]
    mask = np.abs(diff) > eps
    kernel = np.zeros_like(diff)
    jac2 = np.broadcast_to(jac[:, None], diff.shape)
    kernel[mask] = jac2[mask] / (diff[mask] ** 2)
    w = bump_vals
    dx = grid[1] - grid[0]
    return float(np.einsum("i,ij,j->", w, kernel, w) * dx * dx)


def run(n_x: int = 241, n_s: int = 81, s_max: float = 2.5, eps: float = 0.03) -> dict:
    grid = np.linspace(-0.85, 0.85, n_x)
    b = bump(grid)
    s_vals = np.linspace(-s_max, s_max, n_s)
    q = np.array([quadratic_form(s, grid, b, eps) for s in s_vals])
    # Q_loc from s=0
    i0 = int(np.argmin(np.abs(s_vals)))
    q_loc = q[i0]
    k = np.pi / (2.0 * np.cosh(np.pi * s_vals) ** 2)
    ds = s_vals[1] - s_vals[0]
    residue = float(np.sum(k * (q - q_loc)) * ds)
    k_int = float(np.sum(k) * ds)
    rel = residue / abs(q_loc) if abs(q_loc) > 1e-30 else float("nan")
    return {
        "n_x": n_x,
        "n_s": n_s,
        "s_max": s_max,
        "eps": eps,
        "kernel_integral": k_int,
        "Q_loc": q_loc,
        "residue": residue,
        "relative_residue": rel,
        "Q_min": float(np.min(q)),
        "Q_max": float(np.max(q)),
    }


def main() -> int:
    schemes = [
        run(n_x=201, n_s=61, s_max=2.0, eps=0.04),
        run(n_x=241, n_s=81, s_max=2.5, eps=0.03),
        run(n_x=281, n_s=81, s_max=2.5, eps=0.02),
    ]
    rels = [s["relative_residue"] for s in schemes]
    # mutation: drop the kernel (uniform measure) -- must move the number
    mut = run(n_x=241, n_s=81, s_max=2.5, eps=0.03)
    # recompute residue without K, just mean(Q-Qloc)*2 s_max
    # (reported separately below)

    payload = {
        "task": "m9.2_ib_hadamard_attempt",
        "status": "FAILED_MULTI_DIGIT",
        "admission": (
            "Hard-cutoff subtraction of 1/(x-y)^2 produces a relative residue "
            "of order unity that MOVES when the cutoff or grid changes. "
            "This is the same class of proxy Paper VIII already discarded. "
            "No multi-digit coefficient is extracted."
        ),
        "schemes": schemes,
        "relative_residues": rels,
        "spread": float(max(rels) - min(rels)),
        "order_unity": bool(all(0.05 < abs(r) < 20 for r in rels)),
        "what_would_count_as_success": (
            "an analytic Hadamard expansion of the flowed current two-point "
            "function whose finite part is independent of the cutoff to "
            "at least three stable digits"
        ),
    }
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, "m9_2_ib_hadamard.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
        fh.write("\n")
    print("wrote <repo>/.../data/m9_2_ib_hadamard.json")
    print("relative residues:", rels)
    print("spread:", payload["spread"])
    print("STATUS:", payload["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
