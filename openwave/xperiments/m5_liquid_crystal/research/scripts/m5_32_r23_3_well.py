"""M5.32 R23-3 (stretch): what makes the well of report 008's frozen-profile
ladder, read by contraction and by block on the report's committed N = 32 field.

EQUATIONS FIRST (report 008's route 2, as re-typed in m5_32_r22_3_omega_e.py)
---------------
    E(omega) - E(0) = gamma (-2 omega^2 C1 + omega^4 C2),
    C1 = h^3 sum i1s k1,  C2 = h^3 sum k1^2,
    i1s = 4 sum_{i<j} <[A_i, A_j]_eta, .>_X,   k1 = 4 sum_i <[a0, A_i]_eta, .>_X,
    <F, F>_X = F_ab X_ac X_bd F_cd.
The report contracts with X = G = eta - 2 q eta, q the Lagrange projector on the
time eigenvalue; on a block-diagonal field G is the Euclidean metric, so the
report's contraction IS the positive one. The well exists iff C1 > 0.
F is antisymmetric. Split F into its boost block (the 0i entries) and its
rotation block (the ij entries): <F, F>_X = <Fb, Fb> + <Fr, Fr> + 2 <Fb, Fr>.
On a block-diagonal static field every A_i is block-diagonal, so [A_i, A_j]_eta
has no 0i entries (i1s is pure rotation block), and the frozen tangent
a0 = env (W M + M W^T), W the boost-x generator, has only 0i entries, so
[a0, A_i]_eta has only 0i entries (k1 is pure boost block). Then
    X = Euclidean:  k1 = +2 |F_0i|^2 >= 0, i1s >= 0  =>  C1 >= 0, a well for ANY field;
    X = eta:        k1 = -2 |F_0i|^2 <= 0, i1s >= 0  =>  C1 <= 0, no well.
This script measures how far the committed field is from that statement.

Checks: (1) the control, C1 / C2 / omega_E of the record under X = G;
(2) the three contractions G, Frobenius (X = 1), eta, each with the block split of
i1s and k1 and with C1 split by k1's blocks, the pieces summing to the full C1;
(3) pointwise signs: the fraction of cells with i1s < 0 or k1 < 0 per contraction;
(4) the size of the field's M_0i and of G - 1.

AUDIT NOTES (m5_32_r23_0_audit.py, checks B1 to B3; every number reproduced to 1e-15)
- Positivity of i1s and k1 under a positive contraction holds for ANY antisymmetric F
  (<F, F>_X = |R F R^T|^2 for X = R^T R); the block structure is needed only for the
  sign under eta.
- "No well under eta" holds on an exactly block-diagonal field only: white-noise M_0i
  of amplitude 1e-2 on the committed field moves C1(eta) from -9.07e-6 to +4.1e-5
  (it stays negative, -8.7e-6, at 1e-3).
- The well is uninformative as an existence statement: C1 is a sum of products of
  non-negative densities, so any non-uniform field has one (three block-diagonal
  white-noise fields with no soliton return omega_E about 25; the uniform vacuum
  returns C1 = 0). Existence is fixed by the form (i1s - omega^2 k1)^2 and by the
  sign of the contraction on the boost block; only the magnitude could inform, and
  it depends on the a0 normalization, the envelope and h^(-3/2) (R22-3).

Labels (fixed in the plan post): WELL_NEGATIVE_BLOCK if C1 <= 0 under the
positive contraction, WELL_SURVIVES otherwise.

Usage: python m5_32_r23_3_well.py <path to report 008's results/L_ladder>   (seconds)
"""

import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, "..", "data", "m5_32_r23_3_well.json")
spec = importlib.util.spec_from_file_location(
    "r22_3", os.path.join(HERE, "m5_32_r22_3_omega_e.py")
)
R = importlib.util.module_from_spec(spec)
spec.loader.exec_module(R)
ETA = R.ETA


def blocks(F):
    Fb = np.zeros_like(F)
    Fb[..., 0, :], Fb[..., :, 0] = F[..., 0, :], F[..., :, 0]
    return Fb, F - Fb


def inner2(F1, F2, X):
    return np.einsum("...ab,...ac,...bd,...cd->...", F1, X, X, F2)


def densities_split(M, a0, h, delta, metric):
    if metric == "G":
        X = R.G_of(M, delta)
    elif metric == "frobenius":
        X = np.broadcast_to(np.eye(4), M.shape)
    else:
        X = np.broadcast_to(ETA, M.shape)
    out = {k: 0.0 for k in ("i1s_b", "i1s_r", "i1s_x", "k1_b", "k1_r", "k1_x")}
    for st in ("fwd", "bwd"):
        A = [R.d1(M, ax, st, h) for ax in range(3)]
        for i in range(3):
            Fb, Fr = blocks(R.comm(a0, A[i]))
            out["k1_b"] = out["k1_b"] + 2.0 * inner2(Fb, Fb, X)
            out["k1_r"] = out["k1_r"] + 2.0 * inner2(Fr, Fr, X)
            out["k1_x"] = out["k1_x"] + 4.0 * inner2(Fb, Fr, X)
            for j in range(i + 1, 3):
                Fb, Fr = blocks(R.comm(A[i], A[j]))
                out["i1s_b"] = out["i1s_b"] + 2.0 * inner2(Fb, Fb, X)
                out["i1s_r"] = out["i1s_r"] + 2.0 * inner2(Fr, Fr, X)
                out["i1s_x"] = out["i1s_x"] + 4.0 * inner2(Fb, Fr, X)
    return out


def main():
    ext = sys.argv[1]
    h = 1.5
    M = R.pinned_field(
        np.load(os.path.join(ext, "M_G_polished_N32.npz"))["M"],
        np.load(os.path.join(ext, "seeds", "m5_21_2b_end_A_T2_sym_e0_n32_d0.3_pinned.npz"))[
            "M"
        ].astype(np.float64),
    )
    a0, _ = R.tangent(M, h)
    res = {"field": "report 008, M_G_polished_N32 with its pinned shell", "rows": {}}
    ctrl = R.omega_e(M, h)
    res["control"] = {k: ctrl[k] for k in ("C1", "C2", "omega_E")}
    res["control_PASS"] = bool(abs(ctrl["omega_E"] - R.RECORD["omega_E"]) < 1e-7)
    Gm = R.G_of(M, R.DELTA)
    res["max_abs_M0i"] = float(np.abs(M[..., 0, 1:]).max())
    res["max_abs_G_minus_identity"] = float(np.abs(Gm - np.eye(4)).max())
    res["max_abs_a0_spatial_block"] = float(np.abs(a0[..., 1:, 1:]).max())
    res["max_abs_a0_00"] = float(np.abs(a0[..., 0, 0]).max())
    for metric in ("G", "frobenius", "eta"):
        d = densities_split(M, a0, h, R.DELTA, metric)
        i1s = d["i1s_b"] + d["i1s_r"] + d["i1s_x"]
        k1 = d["k1_b"] + d["k1_r"] + d["k1_x"]
        C1 = h**3 * float(np.sum(i1s * k1))
        C2 = h**3 * float(np.sum(k1 * k1))
        row = {
            "C1": C1,
            "C2": C2,
            "omega_E": float(np.sqrt(C1 / C2)) if C1 > 0 else None,
            "well": bool(C1 > 0),
            "C1_by_k1_block": {
                b: h**3 * float(np.sum(i1s * d[f"k1_{b}"])) for b in ("b", "r", "x")
            },
            "i1s_total_by_block": {
                b: h**3 * float(np.sum(d[f"i1s_{b}"])) for b in ("b", "r", "x")
            },
            "k1_total_by_block": {b: h**3 * float(np.sum(d[f"k1_{b}"])) for b in ("b", "r", "x")},
            "fraction_cells_i1s_negative": float(np.mean(i1s < 0)),
            "fraction_cells_k1_negative": float(np.mean(k1 < 0)),
            "most_negative_i1s_over_max": float(i1s.min() / i1s.max()),
        }
        row["blocks_sum_to_C1_rel"] = abs(sum(row["C1_by_k1_block"].values()) - C1) / abs(C1)
        res["rows"][metric] = row
        print(metric, json.dumps(row, indent=1), flush=True)
    g = res["rows"]["G"]
    res["G_matches_control_rel"] = abs(g["C1"] / ctrl["C1"] - 1)
    res["label"] = (
        "WELL_SURVIVES"
        if res["rows"]["frobenius"]["C1"] > 0 and g["C1"] > 0
        else "WELL_NEGATIVE_BLOCK"
    )
    res["PASS"] = bool(
        res["control_PASS"]
        and res["G_matches_control_rel"] < 1e-10
        and all(r["blocks_sum_to_C1_rel"] < 1e-10 for r in res["rows"].values())
    )
    print({k: v for k, v in res.items() if k != "rows"})
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
