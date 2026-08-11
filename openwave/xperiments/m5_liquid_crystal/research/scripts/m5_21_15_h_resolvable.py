"""M5.21.15 addendum: the lattice-resolvable fixed-J point.

The narrow-guard minimizer concentrates dressing action below the
h = 1.5 lattice scale, failing the cross-check gates (the M5.21.14
subgrid-well phenomenon). This run restricts the family to the
RESOLVABLE scales (plateau + bumps rho >= 2.83 only; the small-rho
coefficients pinned to zero), re-minimizes E(J) at J = 2*kin_base*0.5
within guard 0.02, and lattice-checks THAT profile at n = 32:
does a positive-energy fixed-J minimum survive on scales the
certified box can see, and do the gates then pass?

Out: ../data/m5_21_15_resolvable.json
"""
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
_sc = importlib.util.spec_from_file_location(
    "c15", os.path.join(HERE, "m5_21_15_c_fixedj.py"))
C15 = importlib.util.module_from_spec(_sc)
_sc.loader.exec_module(C15)
C14 = B15.C14

BOUND_N = 0.02
KIN_FLOOR = 1.0
FREE_IDX = [0, 5, 6, 7, 8, 9]     # plateau + bumps rho 2.83..16


def expand(sub):
    a = np.zeros(10)
    a[FREE_IDX] = sub
    return a


def main():
    t0 = time.time()
    C14.S_SIGN = -1.0
    grid_opt = C14.make_grid(48, 8, 16)
    grid_full = C14.make_grid(72, 12, 24)
    cc = B15.ChanCorr(grid_opt, 32.0)
    cc_full = B15.ChanCorr(grid_full, 32.0)
    kb = cc_full.kin_base["clock"]
    E_base = float(np.sum(grid_full["wvol"] * cc_full.ec.du_base))
    J = 2.0 * kb * 0.5

    def obj(sub):
        a = expand(sub)
        e, kc = cc.chan_kins(B15.bfun_of(a), ("clock",))
        kt = cc.kin_base["clock"] + kc["clock"]
        if kt <= KIN_FLOOR:
            return e + J * J / (4.0 * KIN_FLOOR) \
                + 1e3 * (KIN_FLOOR - kt) ** 2
        return e + J * J / (4.0 * kt)

    best = None
    for st in ([0.01, 0, 0, 0, 0, 0], [-0.01, 0, 0, 0, 0, 0],
               [0.01] * 6, [-0.01] * 6):
        res = minimize(obj, np.asarray(st, float), method="L-BFGS-B",
                       bounds=[(-BOUND_N, BOUND_N)] * 6,
                       options={"maxiter": 200, "ftol": 1e-14,
                                "eps": 2e-6})
        if best is None or res.fun < best.fun:
            best = res
    avec = expand(np.asarray(best.x))
    e_c, kc = cc_full.chan_kins(B15.bfun_of(avec), ("clock",))
    kt = cc_full.kin_base["clock"] + kc["clock"]
    out = {"J": J, "bound": BOUND_N, "free_idx": FREE_IDX,
           "E_base_u": E_base, "kin_base": kb,
           "avec": avec.tolist(), "E_corr": e_c,
           "kin_tot": kt, "omega_star": J / (2.0 * kt),
           "E_rot": J * J / (4.0 * kt),
           "E_total": E_base + e_c + J * J / (4.0 * kt)}
    out["E_positive"] = bool(out["E_total"] > 0)
    print(json.dumps({k: out[k] for k in
                      ("E_total", "omega_star", "kin_tot",
                       "E_corr", "E_positive")}), flush=True)
    rs = grid_full["rs"]
    bstar = C14.b_of(avec, rs)
    out["LAT_n32"] = C15.lattice_check(rs, bstar, e_c, kc["clock"])
    print(json.dumps({"LAT_n32": out["LAT_n32"]}), flush=True)
    out["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(DATA, "m5_21_15_resolvable.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"done": True, "runtime_s": out["runtime_s"]}))


if __name__ == "__main__":
    main()
