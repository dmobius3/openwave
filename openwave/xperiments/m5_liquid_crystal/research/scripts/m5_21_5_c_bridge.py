"""M5.21.5 G3: the first-principles g bridge from the M5.16 Coulomb
anchor (replaces the structurally-motivated K = 4/alpha).

THE DERIVATION (every factor traced; Gaussian units):

1. The emergent EM sector. Director n, Mermin-Ho curvature
   F_ij = n.(d_i n x d_j n), B_k = 1/2 eps_kij F_ij. For the unit
   hedgehog n = r_hat: F_ij = eps_ijk n_k / r^2, so
   |B_lat|^2 = sum_{i<j} F_ij^2 = 1/r^4 and B_lat = r_hat/r^2
   (unit monopole, flux 4 pi). [P0 gate (i) checks this on-lattice.]

2. The Coulomb anchor (M5.16). The far-field curvature energy
   density is 8 c2/r^4 with c2 = alpha hbar c / 64 pi fixed by
   matching the exterior integral to the EM self-energy
   alpha hbar c / 2r. In lattice units u_lat = 8 c2' |B_lat|^2 with
   the physical map E_phys = C_E E_lat (C_E from the m_e anchor:
   C_E = m_e c^2 / E_lat(electron state)) and x_phys = l x_lat.
   Matching u_phys = C_E u_lat / l^3 to alpha hbar c/(8 pi r_phys^4)
   gives  c2' C_E l = alpha hbar c / 64 pi.        (*)

3. The field normalization is then DERIVED, not assumed:
   u_phys = B_phys^2 / 8 pi  =>  B_phys = lambda B_lat with
   lambda = 8 sqrt(pi c2' C_E / l^3). The emergent monopole charge:
   q = lambda l^2 = 8 sqrt(pi c2' C_E l) = sqrt(alpha hbar c) = e
   EXACTLY by (*): the relaxed hedgehog carries the electron charge
   with zero freedom (the alpha consistency is exact by anchoring).

4. mu. j_phys = (c/4 pi) curl B_phys (+ the displacement term, both
   carried in the lattice j); mu_phys = (1/2c) int r x j_phys d3x.
   Unit-tracing (r_phys = l r_lat, curl_phys = curl_lat / l):
   mu_phys = (lambda l^3 / 4 pi) mu_lat with mu_lat = the EID
   integral 1/2 sum r x j_lat h^3. In Bohr magnetons:
   mu_phys / mu_B = mu_lat (lambda l^3 / 4 pi)(2 m_e c / e hbar)
                  = mu_lat (q l /(2 pi)) (m_e c / e hbar)  [q=e]
                  = mu_lat (l / lambda_C) / (2 pi),
   lambda_C = hbar / m_e c.

5. S. The verified-L clock inertia kin (physical-rate tangent),
   E_kin = omega^2 kin, S_lat = 2 kin (rigid-rotor identity
   J = dE/domega at omega = 1). Physical action: S_phys =
   C_E tau0 S_lat with tau0 = l / c (ASSUMPTION: lattice signal
   speed c_lat = 1, the leapfrog-certified convention), so
   S_phys / hbar = (l / lambda_C) S_lat / E_lat.

6. g. g = (mu_phys / mu_B) / (S_phys / hbar)
       = [mu_lat (l/lambda_C)/(2 pi)] / [(l/lambda_C) S_lat/E_lat]
       = mu_lat E_lat / (2 pi S_lat).
   THE LENGTH UNIT CANCELS: g is computed from three measured
   lattice numbers with NO free factor. The residual assumptions,
   stated: (i) c_lat = 1 enters S_phys; (ii) the Gaussian-units
   energy matching defines B_phys (the E||B dyon duality of the
   hedgehog is carried as in EID); (iii) E_lat anchored to m_e c^2.

Legacy contrast: K = 4/alpha gave g_legacy = (mu/S) * 548
(structurally motivated, never derived). Both assemblies are
reported side by side, plus the "k needed" factor that WOULD close
g = 2.0023 under the first-principles form.

Usage: python3 m5_21_5_c_bridge.py [tags...]
Consumes data/m5_21_5_mu_<tag>.json; writes data/m5_21_5_bridge.json
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

ALPHA_INV = 137.035999
G_MEASURED = 2.0023193


def one_g(mu, S, E):
    if S <= 0 or mu is None:
        return {"mu_lat": mu, "S_lat": S, "note": "S<=0"}
    g_fp = mu * E / (2 * np.pi * S)
    return {
        "mu_lat": mu, "S_lat": S, "E_lat": E,
        "g_first_principles": g_fp,
        "g_vs_measured": g_fp / G_MEASURED,
        "g_legacy_4_over_alpha": (mu / S) * 4.0 * ALPHA_INV,
        "k_needed_for_2.0023": G_MEASURED * S / (2.0 * mu * E)
        if mu > 0 else None,
        "k_derived": 1.0 / (4 * np.pi)}


def assemble(tag):
    with open(os.path.join(DATA, f"m5_21_5_mu_{tag}.json")) as f:
        d = json.load(f)
    E = d["E_lat"]
    rows = {}
    # PRIMARY: clock-induced mu (E-term), state-of-record S (env)
    for pairing, muk, sk in (
            ("mixed", "mu_tilt12", "S_twist"),
            ("matched", "mu_tilt12", "S_tilt12"),
            ("tilt13_mixed", "mu_tilt13", "S_twist")):
        rows[pairing] = one_g(d[muk]["mu_clock"], d[sk]["S_env"], E)
    # sensitivity: total mu (incl. the static curl-B part)
    rows["mixed_mu_total"] = one_g(d["mu_tilt12"]["mu_norm"],
                                   d["S_twist"]["S_env"], E)
    # sensitivity: rigid S (IR-extensive diagnostic)
    rows["mixed_S_rigid"] = one_g(d["mu_tilt12"]["mu_clock"],
                                  d["S_twist"]["S_rigid"], E)
    # sensitivity: the 4D-convention E on 3x3 states
    if "E_lat_4dconv" in d:
        rows["mixed_E_4dconv"] = one_g(d["mu_tilt12"]["mu_clock"],
                                       d["S_twist"]["S_env"],
                                       d["E_lat_4dconv"])
    out = {"tag": tag, "n": d["n"], "is4d": d["is4d"],
           "E_lat": E, "E_src": d.get("E_row_src", "4d_native"),
           "pairings": rows,
           "mu_static_curlB": d["mu_twist"]["mu_curlB_only"],
           "mu_native_twist_clock": d["mu_twist"]["mu_clock"]}
    # the carried J of the fixed-J rung (unit-norm convention), for
    # the kin-convention quantification
    if tag.startswith("fjom"):
        om = tag[4:]
        p = os.path.join(DATA, f"m5_21_9_fixedj_conj_om{om}.json")
        if os.path.exists(p):
            fj = json.load(open(p))
            out["J_carried_unitnorm"] = fj["start"]["J"]
            out["om_target"] = fj["start"]["om_target"]
    return out


def main(tags):
    out = {"tags": tags, "states": {},
           "notes": {
               "g_formula": "g = mu_lat * E_lat / (2 pi S_lat)",
               "assumptions": ["c_lat = 1 (enters S_phys)",
                               "Gaussian energy matching (dyon "
                               "duality as in EID)",
                               "E_lat anchored to m_e c^2"],
               "charge_consistency":
                   "q = 8 sqrt(pi c2' C_E l) = e exact by the "
                   "Coulomb anchor (*): no lattice read needed for "
                   "g (l cancels)"}}
    for t in tags:
        try:
            out["states"][t] = assemble(t)
            print(json.dumps(out["states"][t], indent=1), flush=True)
        except FileNotFoundError:
            print(f"skip {t}: mu json missing", flush=True)
    with open(os.path.join(DATA, "m5_21_5_bridge.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("saved bridge", flush=True)


if __name__ == "__main__":
    tags = sys.argv[1:] or ["fjom0.2", "fjom0.5", "fjom1",
                            "t32A", "t24A", "t48A"]
    main(tags)
