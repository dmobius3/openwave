"""M5.32 R22-3 (stretch): the frozen-profile clock frequency omega_E of
Mikulski's report 008 (the (I_1^G)^2 well), read on the compact W1 x 25
electron, with the report's own N = 32 field as the known-answer control.

EQUATIONS FIRST (report 008, route 2: verify_L_ladder_energies.py, re-typed
here with the spacing h as a parameter; nothing else changed)
---------------
eta = diag(-1, 1, 1, 1), x = eta M, the Lagrange-projector Euclideanizer
    q = x (x - 1) (x - delta) / (g (g - 1) (g - delta)),   G = eta - 2 q eta
<F, F>_G = F_ab G_ac G_bd F_cd. One-sided differences, fwd / bwd averaged:
    i1s(x) = 4 sum_{i<j} <[A_i, A_j]_eta, .>_G       the static I_1^G density
    k1(x)  = 4 sum_i <[a0, A_i]_eta, .>_G            the time density at omega 1
    E(omega) = E_stat + gamma h^3 sum (i1s - omega^2 k1)^2
so E(omega) - E(0) = gamma (-2 omega^2 C1 + omega^4 C2) with
    C1 = h^3 sum i1s k1,   C2 = h^3 sum k1^2,   omega_E = sqrt(C1 / C2)
(the report's chain records carry C1 and C2 without gamma; matched to 1e-6).

Branch translation. Report 008 sits on the branch M_00 = -g (eta M has the time
eigenvalue +g); the R20 to R22 stack sits on the code branch s = -1, M_00 = +g.
For a block-diagonal field (M_0i = 0 exactly) the CURVATURE is the same on both
branches (the 1 x 1 time block drops out of every commutator; audited to 0.0),
and C1, C2, omega_E contain no potential, so the field is carried to the
report's branch by M_00 -> -M_00 before the read. A field with any M_0i above
1e-12 is refused. Audit correction (m5_32_r22_0_audit.py, C9.2b): V4 is NOT
branch-invariant where M_00 departs from g (the cross term of the p = 1, 3
traces flips sign): on the R21 W1 x 25 fields V4 differs by 16 to 18 percent,
so the translated field is not a stationary point of the other branch's static
energy. The read is a frozen-profile read on a borrowed field, nothing more.

What the number is (audit C9.4, C9.5). a0 is normalized by the cell sum of a^2
with no h^3, so omega_E scales as h^(-3/2), and it grows about linearly with
the envelope radius (the report's own control: 0.222 to 3.49 from radius 6 to
no envelope). omega_E is therefore not a continuum quantity; this script reports
it beside (i) the value rescaled to h 1.5, omega_E (h / 1.5)^(3/2), and (ii) the
RATIO to the report's control at the same envelope radius, for radii 6 to 16:
the ratio is the stable read.
Frozen tangent: a0 = env (W M + M W^T) normalized to unit Frobenius norm,
W the boost-x generator, env = exp(-(r / 10)^4) (a physical radius).
gamma = 70.61005836040782 (the report's value, fixed across boxes).

Control (must pass before any other number is read): the report's committed
polished N = 32 field with its pinned shell returns the recorded
C1 9.075435624180151e-06, C2 8.564773350793511e-05, omega_E 0.3255185953964089.

Rows: the R21 W1 x 25 electron (the biaxial vacuum) on the seed pin and on the
free boundary; the R22-1 endpoints (the uniaxial vacuum) where present. The
rungs are NOT relaxed and gamma is NOT refit (the plan post's statement).

Usage: python m5_32_r22_3_omega_e.py <path to report 008's results/L_ladder>
"""

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r22_3_omega_e.json")
SG, DELTA = 8.0, 0.3
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
GAMMA = 70.61005836040782
RECORD = {"C1": 9.075435624180151e-06, "C2": 8.564773350793511e-05, "omega_E": 0.3255185953964089}
T0 = time.time()


def d1(f, ax, st, h):
    out = np.zeros_like(f)
    sl = [slice(None)] * f.ndim
    lo, hi = [slice(None)] * f.ndim, [slice(None)] * f.ndim
    lo[ax], hi[ax] = slice(0, -1), slice(1, None)
    sl[ax] = slice(0, -1) if st == "fwd" else slice(1, None)
    out[tuple(sl)] = (f[tuple(hi)] - f[tuple(lo)]) / h
    return out


def comm(A, B):
    return A @ ETA @ B - B @ ETA @ A


def G_of(M, delta):
    x = np.einsum("ab,...bc->...ac", ETA, M)
    I4 = np.broadcast_to(np.eye(4), M.shape)
    q = (x @ (x - I4) @ (x - delta * I4)) / (SG * (SG - 1) * (SG - delta))
    return ETA - 2.0 * q @ ETA


def inner_pc(F, X):
    return np.einsum("...ab,...ac,...bd,...cd->...", F, X, X, F)


def densities(M, a0, h, delta):
    Gm = G_of(M, delta)
    i1s, k1 = 0.0, 0.0
    for st in ("fwd", "bwd"):
        A = [d1(M, ax, st, h) for ax in range(3)]
        for i in range(3):
            k1 = k1 + 0.5 * 4.0 * inner_pc(comm(a0, A[i]), Gm)
            for j in range(i + 1, 3):
                i1s = i1s + 0.5 * 4.0 * inner_pc(comm(A[i], A[j]), Gm)
    return i1s, k1


def tangent(M, h, renv=10.0):
    N = M.shape[0]
    x = (np.arange(N) - (N - 1) / 2.0) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z)
    env = np.exp(-((r / renv) ** 4))
    W = np.zeros((4, 4))
    W[0, 1] = W[1, 0] = 1.0
    a = env[..., None, None] * (
        np.einsum("ab,...bc->...ac", W, M) + np.einsum("...ab,cb->...ac", M, W)
    )
    return a / np.linalg.norm(a), r


def omega_scan(M, h, delta=DELTA, radii=(6.0, 8.0, 10.0, 12.0, 16.0)):
    out = {}
    for R in radii:
        a0, _ = tangent(M, h, R)
        i1s, k1 = densities(M, a0, h, delta)
        c1, c2 = float(np.sum(i1s * k1)), float(np.sum(k1 * k1))
        out[str(R)] = float(np.sqrt(c1 / c2) * (h / 1.5) ** 1.5) if c1 > 0 and c2 > 0 else None
    return out


def omega_e(M, h, delta=DELTA):
    a0, r = tangent(M, h)
    i1s, k1 = densities(M, a0, h, delta)
    C1 = h**3 * float(np.sum(i1s * k1))
    C2 = h**3 * float(np.sum(k1 * k1))
    k1s = np.sort(k1.ravel())[::-1]
    half = int(np.searchsorted(np.cumsum(k1s), 0.5 * k1s.sum())) + 1
    return {
        "C1": C1,
        "C2": C2,
        "omega_E": float(np.sqrt(C1 / C2)) if C1 > 0 and C2 > 0 else None,
        "dE_on_the_report_rungs": {
            str(w): GAMMA * (-2 * w**2 * C1 + w**4 * C2) for w in (0.1, 0.2, 0.35)
        },
        "depth_at_omega_E": float(GAMMA * C1**2 / C2) if C2 > 0 else None,
        "I1G_static": float(h**3 * np.sum(i1s)),
        "k1_half_mass_cells": half,
        "k1_mean_radius": float(np.sum(k1 * r) / np.sum(k1)),
        "i1s_mean_radius": float(np.sum(i1s * r) / np.sum(i1s)),
        "i1s_fraction_outside_r16": float(np.sum(i1s[r > 16.0]) / np.sum(i1s)),
    }


def pinned_field(M_raw, seed3):
    N = M_raw.shape[0]
    wc = max(1, int(np.ceil(1.6 / 1.5)))
    mask = np.zeros((N, N, N), dtype=bool)
    for ax in range(3):
        sl = [slice(None)] * 3
        sl[ax] = slice(0, wc)
        mask[tuple(sl)] = True
        sl[ax] = slice(N - wc, N)
        mask[tuple(sl)] = True
    seed4 = np.zeros((N, N, N, 4, 4))
    seed4[..., 1:, 1:] = seed3
    seed4[..., 0, 0] = -SG
    Ms = 0.5 * (M_raw + np.swapaxes(M_raw, -1, -2))
    return np.where(mask[..., None, None], seed4, Ms)


def main():
    res = {"gamma": GAMMA, "record": RECORD, "rows": {}}
    ext = sys.argv[1] if len(sys.argv) > 1 else None
    if ext:
        ctrl = {}
        for N in (32, 48):
            fp = os.path.join(ext, f"M_G_polished_N{N}.npz")
            fs = os.path.join(ext, "seeds", f"m5_21_2b_end_A_T2_sym_e0_n{N}_d0.3_pinned.npz")
            if not (os.path.exists(fp) and os.path.exists(fs)):
                continue
            M = pinned_field(np.load(fp)["M"], np.load(fs)["M"].astype(np.float64))
            ctrl[str(N)] = omega_e(M, 1.5)
            ctrl[str(N)]["scan_h15"] = omega_scan(M, 1.5)
            print(f"control N {N}: {ctrl[str(N)]}", flush=True)
        res["control"] = ctrl
        c = ctrl.get("32")
        res["control_PASS"] = bool(
            c
            and abs(c["C1"] / RECORD["C1"] - 1) < 1e-6
            and abs(c["C2"] / RECORD["C2"] - 1) < 1e-6
            and abs(c["omega_E"] - RECORD["omega_E"]) < 1e-7
        )
        print(
            "control PASS:",
            res["control_PASS"],
            "(the report's appendix: 0.326 at N 32, 0.280 at N 48)",
            flush=True,
        )
    fields = []
    r21 = os.path.join(DATA, "m5_32_r21_1")
    for t in (
        "S1_Bseed_g8_d0.3_w25_n32",
        "S1_Bfree_g8_d0.3_w25_n32",
        "S1_Bseed_g8_d0.3_w1_n32",
        "S1_Bfar_g8_d0.3_w1_n32",
        "S1_Bfree_g8_d0.3_w1_n32",
    ):
        fields.append(("R21 " + t, os.path.join(r21, t + ".npz"), 48.0, 0.3))
    r22 = os.path.join(DATA, "m5_32_r22_1")
    if os.path.isdir(r22):
        for f in sorted(os.listdir(r22)):
            if f.endswith(".npz") and "_stage" not in f and "ckpt" not in f:
                L = float(f.split("_L")[1].split(".npz")[0].split("_")[0])
                dl = float(f.split("_d")[1].split("_")[0])
                fields.append(("R22-1 " + f[:-4], os.path.join(r22, f), L, dl))
    for name, path, L, dl in fields:
        if not os.path.exists(path):
            continue
        M = np.load(path)["M"].copy()
        if M[0, 0, 0, 0, 0] > 0:
            if np.abs(M[..., 0, 1:]).max() > 1e-12:
                res["rows"][name] = {
                    "refused": f"not block-diagonal: max |M_0i| {np.abs(M[..., 0, 1:]).max():.2e}"
                }
                print(name, res["rows"][name], flush=True)
                continue
            M[..., 0, 0] *= -1.0
        h = L / M.shape[0]
        res["rows"][name] = dict(omega_e(M, h, dl), h=h, n=int(M.shape[0]), L=L, delta=dl)
        r = res["rows"][name]
        r["omega_E_at_h15"] = r["omega_E"] * (h / 1.5) ** 1.5 if r["omega_E"] else None
        r["scan_h15"] = omega_scan(M, h, dl)
        c32 = res.get("control", {}).get("32", {}).get("scan_h15")
        if c32:
            r["ratio_to_control_by_radius"] = {
                k_: (v / c32[k_] if v and c32.get(k_) else None) for k_, v in r["scan_h15"].items()
            }
            print(
                "   ratio to the report's N 32 control by envelope radius:",
                {
                    k_: (round(v, 2) if v else None)
                    for k_, v in r["ratio_to_control_by_radius"].items()
                },
                flush=True,
            )
        print(
            f"{name:44s} omega_E {r['omega_E']}  C1 {r['C1']:.4e}  C2 {r['C2']:.4e}  <r>_k1 {r['k1_mean_radius']:.2f}  "
            f"i1s outside r16 {r['i1s_fraction_outside_r16']:.3f}",
            flush=True,
        )
    res["runtime_s"] = time.time() - T0
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1)


if __name__ == "__main__":
    main()
