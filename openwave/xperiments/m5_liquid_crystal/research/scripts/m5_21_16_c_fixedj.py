"""M5.21.16 C: boundedness without the guard + the fixed-J electron
under the FLIP functional.

  WIDE   minimize the static dressing correction E_corr over the 10-dim
         smooth b(r) family at WIDE bounds (|a| <= 0.5, vs the M5.21.14
         guard bracket [0.02, 0.05]), eta vs flip, multi-start: does the
         flip close the unbounded-below channel (guard moot)?
  FIXJ   the fixed-J dressed electron under flip: E(amp; J) =
         E_stat_base + E_corr_flip(amp) + J^2 / (4 kin_flip(amp)),
         amp ladder x J ladder; locate the optimal dressing, omega* =
         J / (2 kin), and the energy sign. Under eta the dressing LOWERS
         kin (raising omega*, opening the runaway); under flip it RAISES
         kin, so dressing trades positive static cost against kinetic
         relief: a finite interior optimum is possible.

Family/conventions: C14 (m5_21_14_c_minimize.py) quadrature machinery,
g = 32, s = -1, delta = 0.3. The base family is rotation-only, so base
energies are metric-independent (measured in arm B).

Out: ../data/m5_21_16_fixedj.json
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

_s14 = importlib.util.spec_from_file_location(
    "c14", os.path.join(HERE, "m5_21_14_c_minimize.py"))
C14 = importlib.util.module_from_spec(_s14)
_s14.loader.exec_module(C14)

_sb = importlib.util.spec_from_file_location(
    "b16", os.path.join(HERE, "m5_21_16_b_field.py"))
B16 = importlib.util.module_from_spec(_sb)
_sb.loader.exec_module(B16)


def corr_pair(ec, w, bfun):
    """(E_corr, kin_corr) under (eta, flip) for one dressing."""
    A = ec._A(bfun)
    Qb = C14.qb_from(ec.K_c, ec.K2_c, bfun(ec.r_c))
    a0d = Qb @ ec.a0_base @ np.swapaxes(Qb, -1, -2)
    A_base = ec._A(lambda r: np.zeros_like(r))
    out = {}
    out["eta"] = (
        float(np.sum(w * (C14.dens_u_batch(A) - ec.du_base))),
        float(np.sum(w * (C14.dens_k_batch(a0d, A) - ec.dk_base))))
    du_bf = B16.dens_u_flip(A_base)
    dk_bf = B16.dens_k_flip(ec.a0_base, A_base)
    out["flip"] = (
        float(np.sum(w * (B16.dens_u_flip(A) - du_bf))),
        float(np.sum(w * (B16.dens_k_flip(a0d, A) - dk_bf))))
    return out


def stage_wide(ec, w, rng):
    """multi-start minimization of E_corr over the smooth family."""
    out = {}
    for metric in ("eta", "flip"):
        dens_u = C14.dens_u_batch if metric == "eta" else B16.dens_u_flip
        base = (ec.du_base if metric == "eta"
                else B16.dens_u_flip(ec._A(lambda r: np.zeros_like(r))))

        def obj(avec):
            A = ec._A(lambda r: C14.b_of(avec, r))
            return float(np.sum(w * (dens_u(A) - base)))

        best = None
        starts = []
        for k in range(6):
            a0 = 0.05 * rng.standard_normal(10)
            res = minimize(obj, a0, method="Powell",
                           bounds=[(-0.5, 0.5)] * 10,
                           options={"maxfev": 1500, "xtol": 1e-4})
            starts.append(float(res.fun))
            if best is None or res.fun < best[0]:
                best = (float(res.fun), res.x.tolist())
        out[metric] = {"best_E_corr": best[0], "starts": starts,
                       "best_avec": best[1]}
    out["eta_unbounded_signature"] = bool(
        out["eta"]["best_E_corr"] < -500.0)
    out["flip_bounded"] = bool(out["flip"]["best_E_corr"] > -100.0)
    return out


def stage_fixj(ec, w):
    A_base = ec._A(lambda r: np.zeros_like(r))
    e_stat_base = float(np.sum(w * ec.du_base))
    kin_base = float(np.sum(w * B16.dens_k_flip(ec.a0_base, A_base)))
    amps = [0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4]
    rows = []
    for amp in amps:
        def bf(r, a=amp):
            return a * np.tanh(r / 2.0)
        if amp == 0.0:
            ec_f, kc_f = 0.0, 0.0
        else:
            pair = corr_pair(ec, w, bf)
            ec_f, kc_f = pair["flip"]
        rows.append({"amp": amp, "E_corr_flip": ec_f,
                     "kin_flip": kin_base + kc_f})
    Js = [50.0, 200.0, 800.0]
    table = {}
    for J in Js:
        col = []
        for r in rows:
            kin = r["kin_flip"]
            E = e_stat_base + r["E_corr_flip"] + J * J / (4.0 * kin)
            col.append({"amp": r["amp"], "E_total": E,
                        "omega_star": J / (2.0 * kin)})
        Es = [c["E_total"] for c in col]
        i0 = int(np.argmin(Es))
        table[f"J_{J:g}"] = {
            "ladder": col, "opt_amp": col[i0]["amp"],
            "opt_E": col[i0]["E_total"],
            "opt_omega": col[i0]["omega_star"],
            "interior_optimum": bool(0 < i0 < len(col) - 1),
            "E_positive": bool(col[i0]["E_total"] > 0)}
    return {"e_stat_base": e_stat_base, "kin_base_flip": kin_base,
            "amp_rows": rows, "byJ": table}


def main():
    t0 = time.time()
    rng = np.random.default_rng(21163)
    grid = C14.make_grid(36, 6, 12)
    ec = C14.ExactCorr(grid, 32.0)
    w = grid["wvol"]
    out = {"WIDE": stage_wide(ec, w, rng)}
    print(json.dumps({"WIDE": {m: {k: v for k, v in out["WIDE"][m].items()
                                   if k != "best_avec"}
                               for m in ("eta", "flip")},
                      "eta_unbounded_signature":
                          out["WIDE"]["eta_unbounded_signature"],
                      "flip_bounded": out["WIDE"]["flip_bounded"]},
                     indent=1), flush=True)
    out["FIXJ"] = stage_fixj(ec, w)
    print(json.dumps({"FIXJ_byJ": {k: {kk: vv for kk, vv in v.items()
                                       if kk != "ladder"}
                                   for k, v in out["FIXJ"]["byJ"].items()},
                      "kin_base_flip": out["FIXJ"]["kin_base_flip"],
                      "e_stat_base": out["FIXJ"]["e_stat_base"]},
                     indent=1), flush=True)
    out["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(DATA, "m5_21_16_fixedj.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"runtime_s": out["runtime_s"]}))


if __name__ == "__main__":
    main()
