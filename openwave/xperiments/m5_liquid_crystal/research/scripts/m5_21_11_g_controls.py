"""M5.21.11 g-arm instrument controls (run before trusting the F4 verdict).

Three controls validating m5_21_11_d_garm.py before its F4 read is
accepted (AI_HYGIENE: a terminal criterion from a fresh wrapper is a
hypothesis until the wrapper is verified against non-wrapper ground
truth):

C1 vacuum null      dressing the CONSTANT vacuum d4 must have its
                    minimum AT m = 0 with zero gain (no texture, no
                    spontaneous dressing) and E(m) > 0 away from it
C2 field identity   my Qb-conjugation pipeline applied to the
                    m5_21_8 textured base (Qh d4 Qh^T) must equal
                    m5_21_8_b_lattice.dressed(cfg, m) field-for-field
C3 record match     evaluating the m5_21_8 gladder grid through this
                    stack must reproduce the RECORDED E_min values
                    (data/m5_21_8_lat_gladder.json) at their grid
                    points: 1.1439 (g=8), 1.8534 (g=16), 2.1807 (g=32)

Also records the family-gain read off the m5_21_8 record itself:
E(0) = 62.85 vs E_min ~ 2: family gain ~ -61 FLAT across g = 8-32,
the same q ~ 0 behavior the endpoint g-arm measures (the two
measurements corroborate; neither knew about the other).

Out: ../data/m5_21_11_garm_controls.json
"""
from __future__ import annotations

import importlib.util
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


GM = _load("gm", "m5_21_11_d_garm.py")
L8 = _load("lat8", "m5_21_8_b_lattice.py")
INS4 = GM.INS4

RECORD = {8.0: 1.143928509786395, 16.0: 1.8534164079499864,
          32.0: 2.180748207153326}          # m5_21_8_lat_gladder.json


def main():
    out = {}
    # C1 vacuum null
    c1 = []
    for g in (8.0, 32.0):
        cfg = INS4.base_cfg(s=-1.0, g=g, n=32, L=48.0, delta=0.3)
        n = cfg["n"]
        M4 = np.zeros((n, n, n, 4, 4))
        M4[:] = INS4.vac4(cfg)
        mh = float(np.arctanh(1.0 / g))
        ms = np.linspace(-2.5 * mh, 2.5 * mh, 21)
        Es = [GM.e_dressed(M4, cfg, m) for m in ms]
        i = int(np.argmin(Es))
        c1.append({"g": g, "m_min": float(ms[i]),
                   "E_min": float(Es[i]),
                   "E_edge": float(Es[0]),
                   "pass": bool(abs(ms[i]) < 1e-12
                                and Es[i] < 1e-10
                                and Es[0] > 1.0)})
    out["C1_vacuum_null"] = c1

    # C2 field identity + C3 record match
    c2, c3, fam = [], [], []
    for g in (8.0, 16.0, 32.0):
        cfg = INS4.base_cfg(s=-1.0, g=g, n=32, L=48.0)
        base = L8.dressed(cfg, 0.0)
        m_probe = -0.025
        Mt = L8.dressed(cfg, m_probe)
        Qb = GM.qb_field(cfg, m_probe)
        Mm = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad",
                                 Qb, base, Qb))
        dmax = float(np.abs(Mt - Mm).max())
        c2.append({"g": g, "m": m_probe, "field_diff_max": dmax,
                   "pass": bool(dmax < 1e-10)})
        ms = np.linspace(-0.35, 0.35, 29)
        Es = [float(sum(INS4.e_parts(L8.dressed(cfg, m), cfg)))
              for m in ms]
        emin = float(min(Es))
        c3.append({"g": g, "E_min_reproduced": emin,
                   "E_min_recorded": RECORD[g],
                   "pass": bool(abs(emin - RECORD[g])
                                / RECORD[g] < 1e-6)})
        E0 = float(sum(INS4.e_parts(base, cfg)))
        fam.append({"g": g, "E0_family": E0,
                    "family_gain": emin - E0})
    out["C2_field_identity"] = c2
    out["C3_record_match"] = c3
    out["family_gain_from_record"] = fam
    out["all_pass"] = bool(all(r["pass"] for r in c1 + c2 + c3))
    with open(os.path.join(DATA, "m5_21_11_garm_controls.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    main()
