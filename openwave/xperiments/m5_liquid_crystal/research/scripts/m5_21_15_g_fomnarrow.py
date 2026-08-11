"""M5.21.15 addendum: the narrow-guard money curve.

The full-guard FOM run met its kin constraint only to 2-9 percent
(the deep well rewards mismatch), so the quantitative fixed-J curve
is re-measured in the NARROW-guard regime (bound = 0.02, the
positive-energy regime from the guard ladder): J = 2*kin_base*0.5,
omega swept across the feasible window [J/(2*kin_base), J/(2*kin_min)],
penalty mu = 1e6, feasibility threshold 2 percent. The EJ minimum at
this guard (omega* = 0.592, E_total = +119) should appear as an
INTERIOR minimum of the curve at positive energy: the author's asked
object, measured.

Also: the lattice cross-check re-run on the NARROW-guard EJ
minimizer (features plateau-scale, resolvable at h = 1.5, unlike the
subgrid deep well that failed the full-guard E gate, cf. the
M5.21.14 resolution ladder), at n = 32 and n = 48.

Out: ../data/m5_21_15_fom_narrow.json
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
_sc = importlib.util.spec_from_file_location(
    "c15", os.path.join(HERE, "m5_21_15_c_fixedj.py"))
C15 = importlib.util.module_from_spec(_sc)
_sc.loader.exec_module(C15)
C14 = B15.C14

G_MAIN = 32.0
BOUND_N = 0.02
MU = 1e6
FEAS = 0.02
OMEGAS = (0.48, 0.50, 0.53, 0.56, 0.59, 0.62, 0.66, 0.70)


def main():
    t0 = time.time()
    C14.S_SIGN = -1.0
    grid_opt = C14.make_grid(48, 8, 16)
    grid_full = C14.make_grid(72, 12, 24)
    cc = B15.ChanCorr(grid_opt, G_MAIN)
    cc_full = B15.ChanCorr(grid_full, G_MAIN)
    kb = cc_full.kin_base["clock"]
    E_base = float(np.sum(grid_full["wvol"] * cc_full.ec.du_base))
    J = 2.0 * kb * 0.5

    with open(os.path.join(DATA, "m5_21_15_guard.json")) as f:
        gd = json.load(f)
    warm = next(r["avec"] for r in gd["rows"]
                if r["bound"] == BOUND_N and
                r["omega_target_undressed"] == 0.5)

    out = {"J": J, "bound": BOUND_N, "mu": MU, "E_base_u": E_base,
           "kin_base": kb, "rows": []}
    avec = np.asarray(warm, dtype=float)
    for om in OMEGAS:
        K_t = J / (2.0 * om)

        def obj(a):
            e, kc = cc.chan_kins(B15.bfun_of(a), ("clock",))
            kt = cc.kin_base["clock"] + kc["clock"]
            return e + MU * ((kt - K_t) / max(abs(K_t), 1.0)) ** 2

        best = None
        for st in (avec, [0.01] + [0.0] * 9, [-0.01] + [0.0] * 9,
                   np.zeros(10)):
            res = minimize(obj, np.asarray(st), method="L-BFGS-B",
                           bounds=[(-BOUND_N, BOUND_N)] * 10,
                           options={"maxiter": 200, "ftol": 1e-14,
                                    "eps": 2e-6})
            if best is None or res.fun < best.fun:
                best = res
        avec = np.asarray(best.x)
        e_c, kc = cc_full.chan_kins(B15.bfun_of(avec), ("clock",))
        kt = cc_full.kin_base["clock"] + kc["clock"]
        mism = abs(kt - K_t) / max(abs(K_t), 1.0)
        row = {"omega": om, "kin_target": K_t, "kin_tot": kt,
               "mismatch_rel": mism, "feasible": bool(mism < FEAS),
               "E_corr": e_c,
               "E_total": E_base + e_c + om * om * kt,
               "avec": avec.tolist()}
        out["rows"].append(row)
        print(json.dumps({k: round(row[k], 4) if
                          isinstance(row[k], float) else row[k]
                          for k in ("omega", "kin_tot", "E_total",
                                    "feasible", "mismatch_rel")}),
              flush=True)
    feas = [r for r in out["rows"] if r["feasible"]]
    if len(feas) >= 3:
        E = np.array([r["E_total"] for r in feas])
        om = np.array([r["omega"] for r in feas])
        i = int(np.argmin(E))
        out["interior_minimum"] = bool(0 < i < len(feas) - 1)
        out["omega_min"] = float(om[i])
        out["E_min_total"] = float(E[i])
        out["E_min_positive"] = bool(E[i] > 0)
    # lattice cross-check on the narrow-guard EJ minimizer
    rs = grid_full["rs"]
    bstar = C14.b_of(np.asarray(warm), rs)
    e_c, kc = cc_full.chan_kins(B15.bfun_of(np.asarray(warm)),
                                ("clock",))
    out["LAT_narrow_n32"] = C15.lattice_check(rs, bstar, e_c,
                                              kc["clock"])
    print(json.dumps({"LAT_narrow_n32": out["LAT_narrow_n32"]}),
          flush=True)
    out["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(DATA, "m5_21_15_fom_narrow.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"done": True, "runtime_s": out["runtime_s"]}))


if __name__ == "__main__":
    main()
