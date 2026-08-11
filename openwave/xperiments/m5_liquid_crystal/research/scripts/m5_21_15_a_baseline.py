"""M5.21.15 A1: the undressed retro-gate.

Three parts, all seed-free (the M5.21.3 relaxed endpoints were deleted
under the pre-2026-07-20 dataset rule; the row JSONs survive and are
the baseline record):

  R  RECORD gate: re-assert the M5.21.3 claims from the saved rows
     (m5_21_3_all.json): every rotation-channel kin > 0, boosts < 0,
     the p3 ladders monotone decreasing with every rung depth-bounded
     (stop == max_iter). NEW diagnostic: the envelope-concavity test.
     E*(u), u = omega^2, is min_M of affine functions of u, hence
     CONCAVE in u with non-increasing slopes; the measured ladder
     slopes INCREASE, certifying the rungs as under-converged upper
     bounds (contained-not-converged), the honest read of the
     author's "it stopped minimization process".
  G  instrument gates on the current stack, seed-free variants of the
     M5.21.3 gates: G0 static + G0k kinetic complex-step gradients
     (random fields), G1 SO(1,3) invariance + negative control
     (random smooth field), G2 vacuum == 0 both branches,
     G3 3D-embed regression (random 3x3 field).
  Q  the quadraticity premise on the analytic family (the concavity
     theorem's hypothesis): for the rotating orbit M(t), F_0i is
     exactly linear in omega, so E_kin = omega^2 * kin with kin
     t-independent along the orbit. Measured: kin at t = 0 vs
     t = 0.7 on the analytic family (g = 32, s = -1), and the
     finite-difference orbit velocity vs the analytic a0.

Out: ../data/m5_21_15_baseline.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_s3 = importlib.util.spec_from_file_location(
    "b3", os.path.join(HERE, "m5_21_3_a_4d.py"))
B3 = importlib.util.module_from_spec(_s3)
_s3.loader.exec_module(B3)

_s14 = importlib.util.spec_from_file_location(
    "c14", os.path.join(HERE, "m5_21_14_c_minimize.py"))
C14 = importlib.util.module_from_spec(_s14)
_s14.loader.exec_module(C14)


def stage_record():
    with open(os.path.join(DATA, "m5_21_3_all.json")) as f:
        all3 = json.load(f)
    out = {}
    for s in ("+1", "-1"):
        rows = all3[f"p2_p1_s{s}"]["rows"]
        rot = {nm: r["kin"] for nm, r in rows.items()
               if not nm.startswith("boost")}
        boo = {nm: r["kin"] for nm, r in rows.items()
               if nm.startswith("boost")}
        lad = all3[f"p3_p1_s{s}_boost_z"]["ladder"]
        om = np.array([r["omega"] for r in lad])
        E = np.array([r["E"] for r in lad])
        u = om ** 2
        slopes = np.diff(E) / np.diff(u)
        out[f"s{s}"] = {
            "rot_kin_all_positive": bool(min(rot.values()) > 0),
            "rot_kin": rot,
            "boost_kin_all_negative": bool(max(boo.values()) < 0),
            "boost_kin": boo,
            "ladder_monotone_decreasing": bool(np.all(np.diff(E) < 0)),
            "ladder_all_depth_bounded": bool(all(
                r["stop"] == "max_iter" for r in lad[1:])),
            "u_slopes": slopes.tolist(),
            "u_slopes_increasing": bool(np.all(np.diff(slopes) > 0)),
            "concavity_violated": bool(np.any(np.diff(slopes) > 0)),
        }
        out[f"s{s}"]["reading"] = (
            "slopes in u = omega^2 INCREASE: a true envelope is "
            "concave in u (min of affine functions), so the measured "
            "ladder is an under-converged upper bound per rung, not "
            "the envelope; consistent with every stop == max_iter")
    out["pass"] = all(
        out[f"s{s}"]["rot_kin_all_positive"]
        and out[f"s{s}"]["boost_kin_all_negative"]
        and out[f"s{s}"]["ladder_monotone_decreasing"]
        and out[f"s{s}"]["ladder_all_depth_bounded"]
        for s in ("+1", "-1"))
    return out


def stage_gates():
    rng = np.random.default_rng(2115)
    out = {}
    cfg = B3.base_cfg(n=10, L=15.0)
    M = B3.sym4(rng.normal(size=(10, 10, 10, 4, 4)))
    G = B3.grad(M, cfg)
    errs = []
    for _ in range(4):
        V = B3.sym4(rng.normal(size=M.shape))
        de_an = float(np.sum(G * V))
        t = 1e-30
        de = B3.e_total(M + 1j * t * V, cfg).imag / t
        errs.append(abs(de - de_an) / max(abs(de), 1e-300))
    out["g0_static"] = float(np.max(errs))
    a0 = rng.normal(size=M.shape)
    a0 = a0 / np.sqrt(np.sum(a0 * a0))
    Gk = B3.kin_grad(M, a0, cfg)
    errs = []
    for _ in range(4):
        V = B3.sym4(rng.normal(size=M.shape))
        de_an = float(np.sum(Gk * V))
        t = 1e-30
        de = B3.kin_of(M + 1j * t * V, a0, cfg).imag / t
        errs.append(abs(de - de_an) / max(abs(de), 1e-300))
    out["g0_kin"] = float(np.max(errs))
    # G1 seed-free: random SMOOTH field (band-limited) instead of the
    # deleted 2b endpoint; SO(1,3) invariance is field-agnostic
    from scipy.linalg import expm
    from scipy.ndimage import gaussian_filter
    cfgB = B3.base_cfg(n=14, L=21.0)
    MB = np.stack([[gaussian_filter(
        rng.normal(size=(14, 14, 14)), 2.0)
        for _ in range(4)] for _ in range(4)], axis=-1)
    MB = MB.reshape(14, 14, 14, 4, 4)
    # vacuum + smooth perturbation: keeps the V4 traces near target so
    # the negative control is not diluted by a huge invariant offset
    MB = B3.vac4(cfgB)[None, None, None] + 0.5 * B3.sym4(MB)
    E0 = B3.e_total(MB, cfgB)
    Gm = np.zeros((4, 4))
    Gm[0, 1] = Gm[1, 0] = 0.11
    Gm[2, 3], Gm[3, 2] = -0.23, 0.23
    L = expm(Gm)
    assert np.allclose(L.T @ B3.ETA @ L, B3.ETA, atol=1e-12)
    ML = np.einsum("ab,...bc,dc->...ad", L, MB, L)
    out["g1_so13"] = abs(float(B3.e_total(ML, cfgB)) - E0) / abs(E0)
    Lb = L + 0.05 * rng.normal(size=(4, 4))
    MLb = np.einsum("ab,...bc,dc->...ad", Lb, MB, Lb)
    out["g1_negctrl"] = abs(float(B3.e_total(MLb, cfgB)) - E0) / abs(E0)
    for s in (1.0, -1.0):
        c2 = B3.base_cfg(s=s, n=8, L=12.0)
        Mv = np.zeros((8, 8, 8, 4, 4))
        Mv[:] = B3.vac4(c2)
        eu, ev = B3.e_parts(Mv, c2)
        out[f"g2_vac_s{int(s):+d}"] = [float(eu), float(ev)]
    # G3 seed-free: random smooth 3x3 field embedded block-diagonally
    cfg3 = B3.base_cfg(n=16, L=24.0)
    M3f = np.stack([[gaussian_filter(
        rng.normal(size=(16, 16, 16)), 2.0)
        for _ in range(3)] for _ in range(3)], axis=-1)
    M3f = M3f.reshape(16, 16, 16, 3, 3)
    M3f = 0.5 * (M3f + M3f.swapaxes(-1, -2))
    M4 = B3.embed34(M3f, cfg3)
    eu4, _ = B3.e_parts(M4, cfg3)
    h = cfg3["h"]
    e3 = 0.0
    for br, wt in B3.branches("sym"):
        A = [B3.d1(M3f, ax, h, br) for ax in range(3)]
        for i in range(3):
            for j in range(i + 1, 3):
                C = A[i] @ A[j] - A[j] @ A[i]
                e3 += wt * 4.0 * np.sum(
                    np.einsum("...kl,...kl->...", C, C))
    e3 *= h ** 3
    out["g3_regression"] = abs(float(eu4) - float(e3)) \
        / max(abs(float(e3)), 1e-300)
    out["pass"] = bool(
        out["g0_static"] < 5e-9 and out["g0_kin"] < 5e-9
        and out["g1_so13"] < 1e-9 and out["g1_negctrl"] > 1e-6
        and max(abs(x) for v in (out["g2_vac_s+1"], out["g2_vac_s-1"])
                for x in v) < 1e-16
        and out["g3_regression"] < 1e-12)
    return out


def stage_quad():
    """the concavity theorem's premise on the analytic family:
    kin is t-independent along the clock orbit, and the analytic a0
    matches the finite-difference orbit velocity."""
    grid = C14.make_grid(48, 8, 16)
    g = 32.0
    P, w = grid["P"], grid["wvol"]
    out = {}
    kins = {}
    for t in (0.0, 0.7):
        h = 1e-4
        A = []
        for ax in range(3):
            e = np.zeros(3)
            e[ax] = 1.0
            d = (8.0 * (C14.m4h_batch(P + h * e, g, t=t)
                        - C14.m4h_batch(P - h * e, g, t=t))
                 - (C14.m4h_batch(P + 2 * h * e, g, t=t)
                    - C14.m4h_batch(P - 2 * h * e, g, t=t))) / (12 * h)
            A.append(d)
        a0 = (C14.m4h_batch(P, g, t=t + 1e-5)
              - C14.m4h_batch(P, g, t=t - 1e-5)) / 2e-5
        kins[f"t{t:g}"] = float(np.sum(w * C14.dens_k_batch(a0, A)))
    out["kin_t0"] = kins["t0"]
    out["kin_t0.7"] = kins["t0.7"]
    out["kin_t_invariance_rel"] = abs(kins["t0"] - kins["t0.7"]) \
        / max(abs(kins["t0"]), 1e-300)
    # analytic a0 vs finite-orbit velocity at larger dt: linearity of
    # dM/dt in the orbit parameter (F_0i linear in omega premise)
    a0_ref = C14.a0h_batch(P, g)
    for dt in (0.01, 0.1):
        a0_fd = (C14.m4h_batch(P, g, t=dt)
                 - C14.m4h_batch(P, g, t=-dt)) / (2 * dt)
        out[f"a0_fd_rel_dt{dt:g}"] = float(
            np.max(np.abs(a0_fd - a0_ref))
            / max(np.max(np.abs(a0_ref)), 1e-300))
    out["pass"] = bool(out["kin_t_invariance_rel"] < 1e-6)
    out["reading"] = (
        "kin is t-invariant along the clock orbit and dM/dt is the "
        "fixed a0 scaled by omega, so E(omega) at fixed config is "
        "EXACTLY E_stat + omega^2*kin (no higher omega terms): the "
        "premise of the envelope-concavity theorem; the small "
        "a0_fd deviations at dt = 0.1 are the o(dt^2) central-"
        "difference curvature of the orbit, not an omega^3 energy "
        "term (the energy is exactly quadratic in the VELOCITY)")
    return out


def main():
    t0 = time.time()
    out = {"RECORD": stage_record()}
    print(json.dumps({"RECORD_pass": out["RECORD"]["pass"]},
                     indent=None), flush=True)
    out["GATES"] = stage_gates()
    print(json.dumps({"GATES": {k: v for k, v in out["GATES"].items()
                                if k != "pass"},
                      "GATES_pass": out["GATES"]["pass"]}), flush=True)
    out["QUAD"] = stage_quad()
    print(json.dumps({"QUAD": {k: v for k, v in out["QUAD"].items()
                               if k not in ("pass", "reading")},
                      "QUAD_pass": out["QUAD"]["pass"]}), flush=True)
    out["all_green"] = bool(out["RECORD"]["pass"]
                            and out["GATES"]["pass"]
                            and out["QUAD"]["pass"])
    out["deviation"] = (
        "the M5.21.3 relaxed endpoints (npz) were deleted under the "
        "pre-2026-07-20 dataset rule, so A1 runs seed-free: the "
        "RECORD gate re-asserts the saved rows, the instrument gates "
        "use random smooth fields (field-agnostic checks), and the "
        "premise check runs on the analytic family")
    out["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(DATA, "m5_21_15_baseline.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"all_green": out["all_green"],
                      "runtime_s": out["runtime_s"]}))


if __name__ == "__main__":
    main()
