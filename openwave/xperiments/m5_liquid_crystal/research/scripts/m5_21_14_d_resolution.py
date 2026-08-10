"""M5.21.14 D: the resolution + instrument-agreement addendum for the
failed lattice_E gate.

The c-script minimizer concentrated the static gain in an OSCILLATORY
core structure (alternating bumps to +-1.6 at rho = 0.5-2.8 over a
0.116 plateau): 58% of E_corr below r = 1.3 (the n = 32 lattice's
smallest radius), 87% below r = 3. Two instrument-comparison hygiene
facts drive this addendum's design:
  - a NON-DECAYING plateau makes every lattice read corner-dominated
    (the cube minus the R = 24 sphere is half the box volume), so all
    comparisons run on COMPACT profiles (cos-taper on r in [14, 20]);
  - a resolvable profile must be RE-MINIMIZED in its own subfamily,
    not obtained by zeroing tuned components.

Stages:
  R1  core-mass split of the continuum E_corr density at b* + the
      bulk kin slopes (base vs dressing, per unit R, mid-radius
      window): the bulk-flip criterion in numbers
  R2  instrument agreement: the re-minimized COMPACT RESOLVABLE
      profile (plateau + bumps rho >= 4, tapered) evaluated by the
      continuum quadrature AND the certified lattice (n = 32):
      gates at 10%
  R3  the n-ladder (32/48/64) on the compact FULL b*: lattice E_corr
      descending toward the continuum number as h shrinks (the
      UV-regulator reading of the BND stage)

Out: ../data/m5_21_14_resolution.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import time

import numpy as np
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_specc = importlib.util.spec_from_file_location(
    "cmin", os.path.join(HERE, "m5_21_14_c_minimize.py"))
C = importlib.util.module_from_spec(_specc)
_specc.loader.exec_module(C)
INS4, B8 = C.INS4, C.B8


def taper(r):
    t = np.ones_like(r)
    t[r >= 20.0] = 0.0
    band = (r > 14.0) & (r < 20.0)
    t[band] = 0.5 * (1 + np.cos(np.pi * (r[band] - 14.0) / 6.0))
    return t


RES_RHOS = [4.3666, 6.7420, 10.4097, 16.0]  # the resolvable bumps


def b_res_of(p, r):
    val = p[0] * np.tanh(r / 2.0)
    for k, rho in enumerate(RES_RHOS):
        val = val + p[k + 1] * (r / rho) * np.exp(-((r / rho) ** 2))
    return val * taper(r)


def lattice_eval(bfun, n):
    cfg = INS4.base_cfg(s=C.S_SIGN, g=C.G_MAIN, n=n, L=48.0)
    h = cfg["h"]
    X, Y, Z = INS4.coords(n, h)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    bl = bfun(R.ravel()).reshape(R.shape)
    nx, ny, nz = X / R, Y / R, Z / R
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, bb in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * bb
    Qb = (np.eye(4)[None, None, None] + np.sinh(bl)[..., None, None]
          * K + (np.cosh(bl) - 1.0)[..., None, None] * K2)
    Mb = B8.dressed(cfg, 0.0)
    a0b = B8.a0_unit(cfg, 0.0)
    Md = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad", Qb, Mb, Qb))
    a0d = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad",
                              Qb, a0b, Qb))
    eu0, ev0 = INS4.e_parts(Mb, cfg)
    eud, evd = INS4.e_parts(Md, cfg)
    k0 = INS4.kin_of(Mb, a0b, cfg)
    kd = INS4.kin_of(Md, a0d, cfg)
    return {"n": n, "h": h, "min_r": float(R.min()),
            "E_corr": float(eud + evd - eu0 - ev0),
            "kin_corr": float(kd - k0), "kin_base": float(k0)}


def lattice_masked(bfun, n, rmin):
    """cellwise instrument densities, summed over r > rmin only
    (replicates e_parts / kin_of cellwise, both stencil branches)."""
    cfg = INS4.base_cfg(s=C.S_SIGN, g=C.G_MAIN, n=n, L=48.0)
    h = cfg["h"]
    X, Y, Z = INS4.coords(n, h)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    bl = bfun(R.ravel()).reshape(R.shape)
    nx, ny, nz = X / R, Y / R, Z / R
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, bb in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * bb
    Qb = (np.eye(4)[None, None, None] + np.sinh(bl)[..., None, None]
          * K + (np.cosh(bl) - 1.0)[..., None, None] * K2)
    Mb = B8.dressed(cfg, 0.0)
    a0b = B8.a0_unit(cfg, 0.0)
    Md = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad", Qb, Mb, Qb))
    a0d = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad",
                              Qb, a0b, Qb))
    mask = (R > rmin).astype(float)

    def cell_dens(Mfield, a0field):
        du = np.zeros_like(R)
        dk = np.zeros_like(R)
        for br, wt in (("fwd", 0.5), ("bwd", 0.5)):
            A = [INS4.d1(Mfield, ax, h, br) for ax in range(3)]
            for i in range(3):
                for j in range(i + 1, 3):
                    F = INS4.comm_eta(A[i], A[j])
                    du += wt * 4.0 * INS4.inner_eta(F, F)
            for i in range(3):
                F = INS4.comm_eta(a0field, A[i])
                dk += wt * 4.0 * INS4.inner_eta(F, F)
        return du, dk

    du_d, dk_d = cell_dens(Md, a0d)
    du_b, dk_b = cell_dens(Mb, a0b)
    return {"n": n, "rmin": rmin,
            "E_corr_masked": float(h ** 3 * np.sum(mask
                                                   * (du_d - du_b))),
            "kin_corr_masked": float(h ** 3 * np.sum(mask
                                                     * (dk_d - dk_b)))}


def main():
    t0 = time.time()
    np.seterr(over="ignore", invalid="ignore")
    with open(os.path.join(DATA, "m5_21_14_minimize.json")) as f:
        M = json.load(f)
    avec = np.array(M["avec"])
    bfull_c = lambda r: C.b_of(avec, r) * taper(r)
    out = {"avec": avec.tolist(),
           "profile_read": {"plateau": float(avec[0]),
                            "bump_amps": avec[1:].round(4).tolist(),
                            "bump_rhos": C.RHOS.round(3).tolist()}}

    grid = C.make_grid(72, 12, 24)
    ec = C.ExactCorr(grid, C.G_MAIN)

    # R1: core-mass split + the bulk kin slopes at the plateau level
    du, dk = ec.densities(lambda r: C.b_of(avec, r))
    r = ec.r_c
    w = grid["wvol"]
    split = {}
    for rc in (0.75, 1.3, 3.0):
        w_in = np.where(r < rc, w, 0.0)
        split[f"r_lt_{rc:g}"] = {
            "E_corr": float(np.sum(w_in * du)),
            "kin_corr": float(np.sum(w_in * dk))}
    # radial ledgers: d(kin)/dR in the mid window (plateau physics)
    rs = grid["rs"]
    kin_base_r, kin_corr_r = [], []
    for i, rr in enumerate(rs):
        sel = np.abs(r - rr) < 1e-9
        kin_base_r.append(float(np.sum(w[sel] * ec.dk_base[sel]))
                          / max(np.gradient(rs)[i], 1e-12))
        kin_corr_r.append(float(np.sum(w[sel] * dk[sel]))
                          / max(np.gradient(rs)[i], 1e-12))
    win = (rs > 10.0) & (rs < 20.0)
    slope_base = float(np.mean(np.array(kin_base_r)[win]))
    slope_corr = float(np.mean(np.array(kin_corr_r)[win]))
    out["R1"] = {
        "continuum_full": {
            "E_corr": float(np.sum(w * du)),
            "kin_corr": float(np.sum(w * dk))},
        "core_mass": split,
        "bulk_kin_slopes_per_unit_R": {
            "base": slope_base, "dressing": slope_corr,
            "net": slope_base + slope_corr,
            "window": "r in (10, 20)",
            "reading": ("net < 0 means the dressed constant-omega "
                        "kin ledger DESCENDS with box radius: the "
                        "flip is bulk, not core")}}
    print(json.dumps({"R1": out["R1"]}), flush=True)

    # R2: re-minimized compact resolvable profile, both instruments
    grid_opt = C.make_grid(48, 8, 16)
    ec_opt = C.ExactCorr(grid_opt, C.G_MAIN)
    best = None
    for start in ([0.116, 0.0, 0.0, 0.0, 0.0], [0.03, 0, 0, 0, 0],
                  [0.1, 0.1, -0.1, 0.05, 0.0]):
        res = minimize(
            lambda p: ec_opt.e_corr(lambda rr: b_res_of(p, rr)),
            np.array(start), method="L-BFGS-B",
            options={"maxiter": 120, "ftol": 1e-13, "eps": 1e-5})
        print(json.dumps({"R2_start": start[:2],
                          "E_coarse": float(res.fun)}), flush=True)
        if best is None or res.fun < best.fun:
            best = res
    pres = np.asarray(best.x)
    e_res, k_res = ec.both(lambda rr: b_res_of(pres, rr))
    lat_res = lattice_eval(lambda rr: b_res_of(pres, rr), 32)
    out["R2"] = {"params": pres.tolist(),
                 "continuum": {"E_corr": e_res, "kin_corr": k_res},
                 "lattice_n32": lat_res,
                 "E_rel_dev": abs(lat_res["E_corr"] - e_res)
                 / abs(e_res),
                 "kin_rel_dev": abs(lat_res["kin_corr"] - k_res)
                 / abs(k_res)}
    out["R2"]["gate"] = bool(out["R2"]["E_rel_dev"] < 0.10
                             and out["R2"]["kin_rel_dev"] < 0.10)
    print(json.dumps({"R2": {k: out["R2"][k] for k in
                             ("continuum", "lattice_n32",
                              "E_rel_dev", "kin_rel_dev", "gate")}}),
          flush=True)

    # R3: n-ladder on the compact FULL profile
    e_cf, k_cf = ec.both(bfull_c)
    out["R3"] = {"continuum_target": {"E_corr": e_cf,
                                      "kin_corr": k_cf},
                 "ladder": [lattice_eval(bfull_c, n)
                            for n in (32, 48, 64)]}
    print(json.dumps({"R3": out["R3"]}), flush=True)

    # R4: instrument agreement where BOTH resolve (r > 3), on the
    # resolvable profile: cellwise-masked lattice vs cut continuum
    du_res, dk_res = ec.densities(lambda rr: b_res_of(pres, rr))
    r_c = ec.r_c
    w_out = np.where(r_c > 3.0, grid["wvol"], 0.0)
    cont_out = {"E_corr": float(np.sum(w_out * du_res)),
                "kin_corr": float(np.sum(w_out * dk_res))}
    lat_out = lattice_masked(lambda rr: b_res_of(pres, rr), 32, 3.0)
    out["R4"] = {"continuum_r_gt_3": cont_out,
                 "lattice_masked_n32": lat_out,
                 "E_rel_dev": abs(lat_out["E_corr_masked"]
                                  - cont_out["E_corr"])
                 / max(abs(cont_out["E_corr"]), 1e-12),
                 "kin_rel_dev": abs(lat_out["kin_corr_masked"]
                                    - cont_out["kin_corr"])
                 / max(abs(cont_out["kin_corr"]), 1e-12)}
    out["R4"]["gate"] = bool(out["R4"]["E_rel_dev"] < 0.10
                             and out["R4"]["kin_rel_dev"] < 0.10)
    print(json.dumps({"R4": out["R4"]}), flush=True)

    out["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(DATA, "m5_21_14_resolution.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"R2_gate": out["R2"]["gate"],
                      "R4_gate": out["R4"]["gate"],
                      "runtime_s": out["runtime_s"]}))


if __name__ == "__main__":
    main()
