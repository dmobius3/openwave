"""M5.21.15 addendum: the guard-width ladder for the fixed-J minimum.

The A4 first point showed the fixed-J minimum rides the deep
dressing well (E_corr ~ -4600 at trust region BOUND = 2), making
E_total NEGATIVE: the energy sign of the constrained minimum is set
by the PROVISIONAL guard's width, not by physics (the Q25 fork,
the M5.21.14 unbounded-below structure). This ladder measures the
verdict flip: E(J) minimized within trust regions of increasing
width; report E_total sign, omega*, kin_tot, well depth per rung.

Out: ../data/m5_21_15_guard.json
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
KIN_FLOOR = 1.0
BOUNDS = (0.01, 0.02, 0.05, 0.15, 0.5, 2.0)
OMEGA_TARGETS = (0.2, 0.5)


def ej_min(cc, cc_full, J, bound):
    def obj(a):
        e, kc = cc.chan_kins(B15.bfun_of(a), ("clock",))
        kt = cc.kin_base["clock"] + kc["clock"]
        if kt <= KIN_FLOOR:
            return e + J * J / (4.0 * KIN_FLOOR) \
                + 1e3 * (KIN_FLOOR - kt) ** 2
        return e + J * J / (4.0 * kt)

    best = None
    sc = 0.5 * bound
    for st in ([sc] + [0.0] * 9, [-sc] + [0.0] * 9,
               [sc / 5] + [0.0] * 9):
        res = minimize(obj, np.asarray(st), method="L-BFGS-B",
                       bounds=[(-bound, bound)] * 10,
                       options={"maxiter": 150, "ftol": 1e-13,
                                "eps": min(1e-5, bound / 50)})
        if best is None or res.fun < best.fun:
            best = res
    avec = np.asarray(best.x)
    e_c, kc = cc_full.chan_kins(B15.bfun_of(avec), ("clock",))
    kt = cc_full.kin_base["clock"] + kc["clock"]
    return avec, e_c, kt


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
           "bounds": list(BOUNDS), "rows": []}
    for om_t in OMEGA_TARGETS:
        J = 2.0 * kb * om_t
        for bd in BOUNDS:
            avec, e_c, kt = ej_min(cc, cc_full, J, bd)
            row = {"J": J, "omega_target_undressed": om_t,
                   "bound": bd, "E_corr": e_c, "kin_tot": kt,
                   "E_rot": J * J / (4.0 * kt),
                   "E_total": E_base + e_c + J * J / (4.0 * kt),
                   "omega_star": J / (2.0 * kt),
                   "at_bound": bool(np.any(np.abs(avec)
                                           > 0.999 * bd)),
                   "avec": avec.tolist()}
            row["E_positive"] = bool(row["E_total"] > 0)
            out["rows"].append(row)
            print(json.dumps({k: row[k] for k in
                              ("J", "bound", "E_total", "omega_star",
                               "E_positive", "at_bound")}),
                  flush=True)
            with open(os.path.join(DATA, "m5_21_15_guard.json"),
                      "w") as f:
                json.dump(out, f, indent=1)
    out["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(DATA, "m5_21_15_guard.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"done": True, "runtime_s": out["runtime_s"]}))


if __name__ == "__main__":
    main()
