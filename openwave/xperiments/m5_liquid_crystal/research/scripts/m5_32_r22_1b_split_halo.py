"""M5.32 R22-1 (b): a split-halo hypothesis for the R22-1 descent, REFUTED by
the R22-1 audit and kept as the record of a wrong reading with its reads.

STATUS (m5_32_r22_1_audit.py, Q5): the algebra below holds and the shell reads
reproduce, but the MECHANISM does not. The local balance drops gradient terms of
the same order: the second variation of the curvature on the exact uniaxial
exterior is 8 (1 - delta)^2 / r^2 [2 eps_r^2 + (lam - 2) eps^2 / r^2] with
lam = 4, POSITIVE (the production second difference +25 to +43 in 12 trials
against this script's -0.45 to -3.1; the lowest eigenvalue of the split-only
Hessian on 8 < r < 16 is +0.0236). The halo is box-filling, not 1 / r^2 (log
slopes -0.4 to -0.9 over r 6 to 12, eps r^2 still growing), 78 to 83 percent of
the fall of G_1^2 at r 9 is the contraction of (top - pair mean), not the split,
and A = 8.41 is the frozen-direction coefficient (0.93 on the valley floor).
Check C rests on a row that stopped after 6 iterations and carries no
information. Read the numbers below as shell reads, not as a derivation.

The hypothesis as first written:
The Coulomb tail of the uniaxial hedgehog is unstable to the split of the
degenerate pair.

EQUATIONS FIRST
---------------
R22-0 (b): for a frame field the certified static quartic is
    e = 8 sum_a G_a^2 |E_top[e_a]|^2,   G_a = prod_{b != a} (l_a - l_b).
On the hedgehog of the top eigenvector, |E_top[e_1]|^2 = 1 / r^4. Split the
degenerate pair, (l_1, l_2, l_3) = (1, delta + eps, delta - eps):
    G_1 = (1 - delta - eps)(1 - delta + eps) = (1 - delta)^2 - eps^2,
    G_1^2 = (1 - delta)^4 - 2 (1 - delta)^2 eps^2 + eps^4,
so the Coulomb tail LOSES 16 (1 - delta)^2 eps^2 / r^4 at second order, while
R22-0 (c) showed the potential resists only at fourth order,
    V4 = w A(delta) eps^4,  A = 4 + 36 delta^2 + 144 delta^4.
The two tangent eigenvectors acquire currents of their own with weights
    G_2^2 = G_3^2 = 4 eps^2 (1 - delta -/+ eps)^2 ~ 4 (1 - delta)^2 eps^2,
a cost 32 (1 - delta)^2 eps^2 (|E_top[e_2]|^2 + |E_top[e_3]|^2) / 2 that depends
on the tangent frame the field chooses (the hairy-ball theorem forbids a
current-free frame on the whole sphere, it does not fix the mean). The local
balance WITHOUT the frame cost gives the halo
    eps^2 (r) = 8 (1 - delta)^2 / (w A r^4),   eps = c_0 / r^2,
    c_0 = sqrt(8 (1 - delta)^2 / (w A)).
So on (certified quartic + V4) the uniaxial exterior is a saddle for every r,
the biaxial halo is a power law, and the charge weight runs:
    G_1^2 / (1 - delta)^4 = (1 - eps^2 / (1 - delta)^2)^2.

Checks:
    A  sympy: the expansion of G_1^2, G_2^2, G_3^2; c_0 at the R22-1 settings.
    B  the R22-1 endpoints: the shell means of eps = (l_mid - l_small) / 2 and
       eps r^2, at 3000 and at 9000 polish iterations, against c_0; the seed's
       eps on the same shells (exactly 0 outside the core for 'rad'); the
       measured G_1^2 against (1 - delta)^4 per shell (the running weight).
    C  the delta 0.89 row: c_0 is 0.21 there, and the halo must be absent.

Regenerate: python m5_32_r22_1b_split_halo.py   (seconds; reads the local
arrays data/m5_32_r22_1/*.npz)
"""

import json
import os

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
NPZ = os.path.join(DATA, "m5_32_r22_1")
OUT_JSON = os.path.join(DATA, "m5_32_r22_1b_split_halo.json")
W1 = 0.000724023879
SHELLS = (4.5, 6.0, 9.0, 12.0, 15.0, 18.0)


def check_a():
    d, e = sp.symbols("delta epsilon", positive=True)
    lam = (1, d + e, d - e)
    G = [sp.expand(sp.prod([lam[a] - lam[b] for b in range(3) if b != a])) for a in range(3)]
    G1sq = sp.expand(G[0] ** 2)
    ok = sp.simplify(G1sq - ((1 - d) ** 4 - 2 * (1 - d) ** 2 * e**2 + e**4)) == 0
    ok &= sp.simplify(G[1] ** 2 - 4 * e**2 * (1 - d - e) ** 2) == 0
    ok &= sp.simplify(G[2] ** 2 - 4 * e**2 * (1 - d + e) ** 2) == 0
    out = {
        "G1_sq": str(G1sq),
        "G2_sq": str(sp.factor(G[1] ** 2)),
        "G3_sq": str(sp.factor(G[2] ** 2)),
        "PASS": bool(ok),
    }
    for dl, w1s in ((0.3, 25.0), (0.3, 1.0), (0.89, 25.0)):
        A = 4 + 36 * dl**2 + 144 * dl**4
        out[f"c0_delta{dl}_w{w1s:g}"] = float(np.sqrt(8 * (1 - dl) ** 2 / (W1 * w1s * A)))
    print("A:", out)
    return out


def shells_of(M, L):
    n = M.shape[0]
    h = L / n
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z)
    lam = np.linalg.eigvalsh(M[..., 1:, 1:])
    rows = []
    for R in SHELLS:
        sh = np.abs(r - R) < 0.75 * h
        eps = 0.5 * (lam[..., 1] - lam[..., 0])[sh]
        top, mid, small = lam[..., 2][sh], lam[..., 1][sh], lam[..., 0][sh]
        G1sq = ((top - mid) * (top - small)) ** 2
        rows.append(
            {
                "r": R,
                "eps_mean": float(eps.mean()),
                "eps_r2": float(eps.mean() * R * R),
                "G1_sq_mean": float(G1sq.mean()),
                "lam_top": float(top.mean()),
                "pair_mean": float(0.5 * (mid + small).mean()),
            }
        )
    return rows


def check_b():
    out = {}
    tags = [
        ("rad 3000", "rad_pin_d0.3_w25_n32_L48"),
        ("rad 9000", "rad_pin_d0.3_w25_n32_L48_x9000"),
        ("bia 3000", "bia_pin_d0.3_w25_n32_L48"),
        ("bia 9000", "bia_pin_d0.3_w25_n32_L48_x9000"),
        ("perm 9000", "perm_pin_d0.3_w25_n32_L48_x9000"),
        ("obl 9000", "obl_pin_d0.3_w25_n32_L48_x9000"),
        ("rad free 3000", "rad_free_d0.3_w25_n32_L48"),
        ("rad W1x1 pin (at the gate)", "rad_pin_d0.3_w1_n32_L48"),
        ("rad delta 0.89 (at the gate)", "rad_pin_d0.89_w25_n32_L48"),
    ]
    for name, tag in tags:
        path = os.path.join(NPZ, tag + ".npz")
        if not os.path.exists(path):
            continue
        rows = shells_of(np.load(path)["M"], 48.0)
        out[name] = rows
        print(
            f"B {name:30s} eps r^2: "
            + "  ".join(f"r{q['r']:g}: {q['eps_r2']:6.2f}" for q in rows)
            + "   G1^2: "
            + "  ".join(f"{q['G1_sq_mean']:.3f}" for q in rows)
        )
    return out


def main():
    res = {
        "a": check_a(),
        "b": check_b(),
        "uniaxial_G1_sq_delta0.3": 0.7**4,
        "uniaxial_G1_sq_delta0.89": 0.11**4,
    }
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
