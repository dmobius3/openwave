"""M5.21.15 A4: the fixed-J bridge: the constrained minimum at
nonzero omega, dressed (clock channel, s = -1, g = 32).

The envelope-concavity theorem (task doc) forbids a FREE interior
omega-minimum; the constraint that breaks it is fixed angular
momentum J = dE/domega = 2 omega kin (the author's "for electron
analog WITH angular momentum"). Three parts:

  EJ   the constrained envelope E(J) = E_base
       + min_b [E_corr(b) + J^2 / (4 kin_tot(b))], kin_tot > 0
       (the guard is structural here: at fixed J the J^2/(4 kin)
       term PENALIZES small kin, so the minimization pushes kin UP,
       away from the runaway channel). Per J: b*, kin_tot,
       omega* = J/(2 kin_tot), E > 0 check, and the clock
       thermodynamics dE/dJ = omega* (numeric derivative).
  FOM  the money curve: E(omega) at fixed J traced by the penalty-
       constrained family (kin_tot(b) = J/(2 omega)), the measured
       "energy minimum at positive energy and nonzero omega".
       Infeasible rungs (target kin outside the family's reach)
       reported as such, never silently dropped.
  LAT  certified-lattice cross-check of b*(J_mid) on the n = 32,
       L = 48 box (the M5.21.14 pattern): E_corr and kin_corr
       relative deviations < 10 percent gates.

Units: family units throughout (a0 un-normalized, the M5.21.14
convention). Out: ../data/m5_21_15_fixedj.json
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

_sb = importlib.util.spec_from_file_location(
    "b15", os.path.join(HERE, "m5_21_15_b_coupled.py"))
B15 = importlib.util.module_from_spec(_sb)
_sb.loader.exec_module(B15)
C14 = B15.C14

G_MAIN = 32.0
KIN_FLOOR = 1.0                     # hard floor: kin_tot > this
OMEGA_TARGETS = (0.1, 0.2, 0.5, 1.0, 1.5)


def ej_point(cc, cc_full, J, warm):
    """min over b of E_corr + J^2/(4 kin_tot), kin_tot floored."""

    def obj(a):
        e, kc = cc.chan_kins(B15.bfun_of(a), ("clock",))
        kt = cc.kin_base["clock"] + kc["clock"]
        if kt <= KIN_FLOOR:
            # smooth barrier: continuous at kt = KIN_FLOOR, pushes back
            return e + J * J / (4.0 * KIN_FLOOR) \
                + 1e3 * (KIN_FLOOR - kt) ** 2
        return e + J * J / (4.0 * kt)

    # E_corr and kin_corr are EVEN in b (measured: the +-plateau probe
    # coincide), so avec = 0 is a stationary point: zero-start alone
    # is blind, the +-0.01 plateau starts break the symmetry
    best = None
    for st in (warm, [0.01] + [0.0] * 9, [-0.01] + [0.0] * 9,
               np.zeros(10)):
        res = minimize(obj, np.asarray(st), method="L-BFGS-B",
                       bounds=[(-B15.BOUND, B15.BOUND)] * 10,
                       options={"maxiter": 150, "ftol": 1e-13,
                                "eps": 1e-5})
        if best is None or res.fun < best.fun:
            best = res
    avec = np.asarray(best.x)
    e_c, kc = cc_full.chan_kins(B15.bfun_of(avec), ("clock",))
    kt = cc_full.kin_base["clock"] + kc["clock"]
    return {"J": J, "E_corr": e_c, "kin_tot": kt,
            "E_rot": J * J / (4.0 * kt),
            "E_over_base": e_c + J * J / (4.0 * kt),
            "omega_star": J / (2.0 * kt),
            "avec": avec.tolist(),
            "at_bound": bool(np.any(np.abs(avec)
                                    > 0.999 * B15.BOUND))}


def fom_curve(cc, cc_full, J, omegas, warm):
    """E(omega) at fixed J via the penalty-constrained family."""
    rows = []
    avec = np.asarray(warm, dtype=float)
    for om in omegas:
        K_t = J / (2.0 * om)
        mu = 1e4

        def obj(a):
            e, kc = cc.chan_kins(B15.bfun_of(a), ("clock",))
            kt = cc.kin_base["clock"] + kc["clock"]
            return e + mu * ((kt - K_t) / max(abs(K_t), 1.0)) ** 2

        best = None
        for st in (avec, [0.01] + [0.0] * 9, [-0.01] + [0.0] * 9):
            res = minimize(obj, np.asarray(st), method="L-BFGS-B",
                           bounds=[(-B15.BOUND, B15.BOUND)] * 10,
                           options={"maxiter": 150, "ftol": 1e-13,
                                    "eps": 1e-5})
            if best is None or res.fun < best.fun:
                best = res
        avec = np.asarray(best.x)
        e_c, kc = cc_full.chan_kins(B15.bfun_of(avec), ("clock",))
        kt = cc_full.kin_base["clock"] + kc["clock"]
        mismatch = abs(kt - K_t) / max(abs(K_t), 1.0)
        rows.append({"omega": om, "kin_target": K_t,
                     "kin_tot": kt, "feasible": bool(mismatch < 0.01),
                     "mismatch_rel": mismatch,
                     "E_corr": e_c,
                     "E_over_base": e_c + om * om * kt,
                     "avec": avec.tolist()})
        print(json.dumps({"FOM": {"omega": om,
                                  "E_over_base":
                                      rows[-1]["E_over_base"],
                                  "feasible":
                                      rows[-1]["feasible"]}}),
              flush=True)
    feas = [r for r in rows if r["feasible"]]
    res = {"J": J, "rows": rows, "n_feasible": len(feas)}
    if len(feas) >= 3:
        E = np.array([r["E_over_base"] for r in feas])
        om = np.array([r["omega"] for r in feas])
        i = int(np.argmin(E))
        res["interior_minimum"] = bool(0 < i < len(feas) - 1)
        res["omega_min"] = float(om[i])
        res["E_min_over_base"] = float(E[i])
    return res


def lattice_check(bstar_rs, bstar, e_corr_cont, kin_corr_cont):
    INS4, B8 = C14.INS4, C14.B8
    cfg = INS4.base_cfg(s=-1.0, g=G_MAIN, n=32, L=48.0)
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS4.coords(n, h)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    bl = np.interp(R.ravel(), bstar_rs, bstar).reshape(R.shape)
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
    k0l = INS4.kin_of(Mb, a0b, cfg)
    kdl = INS4.kin_of(Md, a0d, cfg)
    lat = {"E_corr_lattice": float(eud + evd - eu0 - ev0),
           "kin_corr_lattice": float(kdl - k0l),
           "kin_base_lattice": float(k0l),
           "E_rel_dev": float(abs((eud + evd - eu0 - ev0)
                                  - e_corr_cont)
                              / max(abs(e_corr_cont), 1e-300)),
           "kin_rel_dev": float(abs((kdl - k0l) - kin_corr_cont)
                                / max(abs(kin_corr_cont), 1e-300))}
    lat["gates"] = {"E": bool(lat["E_rel_dev"] < 0.10),
                    "kin": bool(lat["kin_rel_dev"] < 0.10)}
    return lat


def main():
    t0 = time.time()
    C14.S_SIGN = -1.0
    grid_opt = C14.make_grid(48, 8, 16)
    grid_full = C14.make_grid(72, 12, 24)
    cc = B15.ChanCorr(grid_opt, G_MAIN)
    cc_full = B15.ChanCorr(grid_full, G_MAIN)
    kb = cc_full.kin_base["clock"]
    E_base = float(np.sum(grid_full["wvol"] * cc_full.ec.du_base))
    out = {"kin_base_clock": kb, "E_base_u": E_base,
           "kin_floor": KIN_FLOOR}
    print(json.dumps({"kin_base_clock": kb, "E_base_u": E_base}),
          flush=True)

    with open(os.path.join(DATA, "m5_21_14_minimize.json")) as f:
        warm14 = json.load(f)["avec"]

    # EJ: the constrained envelope
    ej = []
    warm = warm14
    for om_t in OMEGA_TARGETS:
        J = 2.0 * kb * om_t
        row = ej_point(cc, cc_full, J, warm)
        row["omega_target_undressed"] = om_t
        row["E_total"] = E_base + row["E_over_base"]
        row["E_positive"] = bool(row["E_total"] > 0)
        warm = row["avec"]
        ej.append(row)
        print(json.dumps({"EJ": {k: row[k] for k in
                                 ("J", "E_total", "omega_star",
                                  "kin_tot", "E_positive",
                                  "at_bound")}}), flush=True)
    # clock thermodynamics: dE/dJ vs omega* (central differences)
    Js = np.array([r["J"] for r in ej])
    Es = np.array([r["E_over_base"] for r in ej])
    dEdJ = np.gradient(Es, Js)
    for r, d in zip(ej, dEdJ):
        r["dEdJ_numeric"] = float(d)
        r["dEdJ_over_omega_star"] = float(d / r["omega_star"])
    out["EJ"] = ej

    # FOM: the fixed-J energy-vs-omega curve at the middle J
    J_mid = 2.0 * kb * 0.5
    omegas = (0.15, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1.0, 1.25, 1.5)
    warm_mid = next(r["avec"] for r in ej
                    if r["omega_target_undressed"] == 0.5)
    out["FOM"] = fom_curve(cc, cc_full, J_mid, omegas, warm_mid)
    out["FOM"]["E_base_u"] = E_base

    # LAT: certified-lattice cross-check at the mid-J minimizer
    avec_mid = np.asarray(warm_mid)
    rs = grid_full["rs"]
    bstar = C14.b_of(avec_mid, rs)
    e_c, kc = cc_full.chan_kins(B15.bfun_of(avec_mid), ("clock",))
    out["LAT"] = lattice_check(rs, bstar, e_c, kc["clock"])
    print(json.dumps({"LAT": out["LAT"]}), flush=True)

    out["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(DATA, "m5_21_15_fixedj.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"done": True, "runtime_s": out["runtime_s"]}))


if __name__ == "__main__":
    main()
