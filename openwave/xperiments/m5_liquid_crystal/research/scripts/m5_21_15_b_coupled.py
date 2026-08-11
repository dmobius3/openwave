"""M5.21.15 A2 + A3: the coupled dressed omega-scan and the
angular-momentum (channel) read on the analytic vacuum-hedgehog
family, exact functional at finite g (the M5.21.14 machinery).

A2 (the free coupled envelope): E*(omega) = min over the guarded
smooth b(r) family (10-dim, trust region |a_k| <= BOUND) of
    E_corr(b) + omega^2 * (kin_base + kin_corr(b))
per channel, omega swept THROUGH AND PAST the M5.21.3 probed range
(0..2.0), warm-started along the ladder. The envelope-concavity
theorem (task doc) predicts NO interior omega-minimum: measured
slopes in u = omega^2 must be non-increasing (concave), so the scan
either sits at omega = 0 (kin_tot > 0) or runs away (kin_tot < 0).
The scan MEASURES this rather than assuming it. Trust-region hits
are flagged (runaway-within-family, the M5.21.14 guard discipline).

A3 (the channel read): per channel (the internal clock a0h, the
global rotation Jz, the global boost Kz), minimize kin_tot(b) alone
within the family: does the dressing turn any ROTATION-sector
omega^2 coefficient negative? J = dE/domega = 2 omega kin_tot
reported on every candidate.

Dressed-channel readings: the internal clock uses the M5.21.14
reading (the core clock ticks inside the dressing, a0d = Qb a0 Qb^T);
global channels use the whole-config reading (a0d = G Md + Md G^T,
Md = Qb M Qb^T). Frequencies are in family units (a0 un-normalized,
matching M5.21.14; the M5.21.3 catalog normalized a0, so magnitudes
are not directly comparable across the two records, signs are).

R-axis: post-hoc cuts of every found b* at R = 8, 12, 16, 24
(the M5.21.14 bulk-flip threshold sits at the certified L = 48 box).

Out: ../data/m5_21_15_coupled.json
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

BOUND = 2.0                       # trust region per family coefficient
OMEGAS = (0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0)
RCUTS = (8.0, 12.0, 16.0, 24.0)

JZ = np.zeros((4, 4)); JZ[1, 2], JZ[2, 1] = -1.0, 1.0
KZ = np.zeros((4, 4)); KZ[0, 3] = KZ[3, 0] = 1.0
CHANNELS = ("clock", "rot_z", "boost_z")


class ChanCorr:
    """channel-resolved corrections on top of C14.ExactCorr."""

    def __init__(self, grid, g):
        self.ec = C14.ExactCorr(grid, g)
        self.grid, self.g = grid, g
        P = grid["P"]
        self.M_base = C14.m4h_batch(P, g)
        self.A_base = self.ec._A(lambda r: np.zeros_like(r))
        self.a0_chan = {"clock": self.ec.a0_base,
                        "rot_z": JZ @ self.M_base + self.M_base @ JZ.T,
                        "boost_z": KZ @ self.M_base
                        + self.M_base @ KZ.T}
        self.kin_base = {
            nm: float(np.sum(grid["wvol"]
                             * C14.dens_k_batch(a0, self.A_base)))
            for nm, a0 in self.a0_chan.items()}

    def _dressed(self, bfun):
        Qb = C14.qb_from(self.ec.K_c, self.ec.K2_c,
                         bfun(self.ec.r_c))
        A = self.ec._A(bfun)
        return Qb, A

    def chan_kins(self, bfun, channels=CHANNELS):
        """(E_corr, {channel: kin_corr}) at b."""
        Qb, A = self._dressed(bfun)
        w = self.grid["wvol"]
        du = C14.dens_u_batch(A) - self.ec.du_base
        e_corr = float(np.sum(w * du))
        Md = Qb @ self.M_base @ np.swapaxes(Qb, -1, -2)
        kc = {}
        for nm in channels:
            if nm == "clock":
                a0d = Qb @ self.ec.a0_base @ np.swapaxes(Qb, -1, -2)
            else:
                G = JZ if nm == "rot_z" else KZ
                a0d = G @ Md + Md @ G.T
            dk = C14.dens_k_batch(a0d, A) - C14.dens_k_batch(
                self.a0_chan[nm], self.A_base)
            kc[nm] = float(np.sum(w * dk))
        return e_corr, kc

    def cut_table(self, bfun, channel):
        """per-R-cut (E_corr, kin_corr, kin_base) for one channel."""
        Qb, A = self._dressed(bfun)
        du = C14.dens_u_batch(A) - self.ec.du_base
        Md = Qb @ self.M_base @ np.swapaxes(Qb, -1, -2)
        if channel == "clock":
            a0d = Qb @ self.ec.a0_base @ np.swapaxes(Qb, -1, -2)
        else:
            G = JZ if channel == "rot_z" else KZ
            a0d = G @ Md + Md @ G.T
        dkd = C14.dens_k_batch(a0d, A)
        dkb = C14.dens_k_batch(self.a0_chan[channel], self.A_base)
        out = {}
        for rc in RCUTS:
            w = np.where(self.ec.r_c <= rc, self.grid["wvol"], 0.0)
            out[f"R{rc:g}"] = {
                "E_corr": float(np.sum(w * du)),
                "kin_corr": float(np.sum(w * (dkd - dkb))),
                "kin_base": float(np.sum(w * dkb))}
        return out


def bfun_of(avec):
    a = np.asarray(avec)
    return lambda r: C14.b_of(a, r)


def scan_channel(cc, cc_full, channel, warm, tag):
    """the free coupled omega-ladder for one channel."""
    kb = cc.kin_base[channel]
    ladder = []
    avec = np.array(warm, dtype=float)
    for om in OMEGAS:
        u = om * om

        def obj(a):
            e, kc = cc.chan_kins(bfun_of(a), (channel,))
            return e + u * (kb + kc[channel])

        # E_corr / kin_corr are EVEN in b, so avec = 0 is stationary:
        # the +-0.01 plateau starts break the symmetry (zero-start
        # alone is blind)
        starts = [avec, [0.01] + [0.0] * 9, [-0.01] + [0.0] * 9]
        best = None
        for st in starts:
            res = minimize(obj, np.asarray(st), method="L-BFGS-B",
                           bounds=[(-BOUND, BOUND)] * 10,
                           options={"maxiter": 120, "ftol": 1e-13,
                                    "eps": 1e-5})
            if best is None or res.fun < best.fun:
                best = res
        avec = np.asarray(best.x)
        e_c, kc = cc_full.chan_kins(bfun_of(avec), (channel,))
        kb_full = cc_full.kin_base[channel]
        kin_tot = kb_full + kc[channel]
        at_bound = bool(np.any(np.abs(avec) > 0.999 * BOUND))
        row = {"omega": om, "E_corr": e_c,
               "kin_corr": kc[channel], "kin_base": kb_full,
               "kin_tot": kin_tot,
               "E_env": e_c + u * kin_tot,
               "J": 2.0 * om * kin_tot,
               "avec": avec.tolist(), "at_bound": at_bound,
               "converged": bool(best.success)}
        ladder.append(row)
        print(json.dumps({"tag": tag, "omega": om,
                          "E_env": row["E_env"],
                          "kin_tot": kin_tot,
                          "at_bound": at_bound}), flush=True)
    E = np.array([r["E_env"] for r in ladder])
    us = np.array([r["omega"] ** 2 for r in ladder])
    slopes = np.diff(E) / np.diff(us)
    interior_min = bool(np.any((E[1:-1] < E[:-2])
                               & (E[1:-1] < E[2:])))
    return {"ladder": ladder, "u_slopes": slopes.tolist(),
            "concave_certified": bool(np.all(np.diff(slopes) <= 1e-9)),
            "interior_minimum_found": interior_min,
            "monotone": ("decreasing" if np.all(np.diff(E) < 0)
                         else "increasing" if np.all(np.diff(E) > 0)
                         else "mixed")}


def min_kin_channel(cc, cc_full, channel):
    """A3: minimize kin_tot(b) alone within the family."""
    kb = cc.kin_base[channel]

    def obj(a):
        _, kc = cc.chan_kins(bfun_of(a), (channel,))
        return kb + kc[channel]

    best = None
    for st in ([0.01] + [0.0] * 9, [0.026] + [0.0] * 9,
               [-0.026] + [0.0] * 9):
        res = minimize(obj, np.asarray(st), method="L-BFGS-B",
                       bounds=[(-BOUND, BOUND)] * 10,
                       options={"maxiter": 120, "ftol": 1e-13,
                                "eps": 1e-5})
        if best is None or res.fun < best.fun:
            best = res
    avec = np.asarray(best.x)
    e_c, kc = cc_full.chan_kins(bfun_of(avec), (channel,))
    kt = cc_full.kin_base[channel] + kc[channel]
    return {"kin_tot_min": kt, "E_corr_there": e_c,
            "kin_base": cc_full.kin_base[channel],
            "avec": avec.tolist(),
            "at_bound": bool(np.any(np.abs(avec) > 0.999 * BOUND)),
            "turns_negative": bool(kt < 0),
            "cut_table": cc_full.cut_table(bfun_of(avec), channel)}


def main():
    t0 = time.time()
    out = {"BOUND": BOUND, "omegas": list(OMEGAS)}
    grid_opt = C14.make_grid(48, 8, 16)
    grid_full = C14.make_grid(72, 12, 24)

    # warm start from the M5.21.14 minimizer
    with open(os.path.join(DATA, "m5_21_14_minimize.json")) as f:
        warm14 = json.load(f)["avec"]

    for s, g, chans, tag in ((-1.0, 32.0, CHANNELS, "main_s-1_g32"),
                             (+1.0, 32.0, ("clock",), "s+1_g32"),
                             (-1.0, 8.0, ("clock",), "s-1_g8")):
        C14.S_SIGN = s
        cc = ChanCorr(grid_opt, g)
        cc_full = ChanCorr(grid_full, g)
        blk = {"s": s, "g": g,
               "kin_base": cc_full.kin_base}
        for ch in chans:
            warm = warm14 if (s < 0 and g == 32.0) else [0.0] * 10
            blk[f"scan_{ch}"] = scan_channel(cc, cc_full, ch, warm,
                                             f"{tag}_{ch}")
            blk[f"minkin_{ch}"] = min_kin_channel(cc, cc_full, ch)
            print(json.dumps({"tag": tag, "minkin": ch,
                              "kin_tot_min":
                                  blk[f"minkin_{ch}"]["kin_tot_min"],
                              "turns_negative":
                                  blk[f"minkin_{ch}"]
                                  ["turns_negative"]}), flush=True)
        out[tag] = blk
        with open(os.path.join(DATA, "m5_21_15_coupled.json"),
                  "w") as f:
            json.dump(out, f, indent=1)   # eager checkpoint per block

    out["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(DATA, "m5_21_15_coupled.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"done": True, "runtime_s": out["runtime_s"]}))


if __name__ == "__main__":
    C14.S_SIGN = -1.0
    main()
