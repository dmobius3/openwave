"""M5.21.16 A: symbolic reproduction of the author's 2026-08-13 notebook
(theory/duda_fmunu_4d_hamiltonian_imaginary_2026_08_13.pdf) + the bridge
lemma to the certified lattice functional.

The notebook's conventions (reproduced exactly):
  Gamma_mu = 4x4 generator with boost components tG_mu^j symmetric in
  row/col 0 and rotation components G_mu^j in the so(3) spatial block;
  xi = diag(-1,1,1,1); coms(A,B) = A xi B - B xi A; d = diag(g,1,delta,0);
  M_mu = coms(Gamma_mu, d); F_munu = coms(M_mu, M_nu);
  H = sum_{mu,nu} [ (F_munu)_23^2 + (F_munu)_24^2 + (F_munu)_34^2
                  - (F_munu)_12^2 - (F_munu)_13^2 - (F_munu)_14^2 ]
  (1-indexed internal entries; signature weighting on INTERNAL indices),
  expanded as a parabola in omega = G_0^1, then g -> 1/delta at leading
  order (the notebook's Series[..., {delta, 0, 0}]).

Checks:
  N1  baseline omega^2 coefficient == -2 sum of six (tG_i^{2,3})^2
  N2  variant A (imaginary tG) omega^2 coefficient == +2 (same six)
  N3  variant B (imaginary g + conjugated second F) coefficient ==
      2 (2 tG_0^2^2 + 2 tG_0^3^2 + the six)   [positive-definite]
  N4  the FLIP (all-plus internal contraction, tr(F F^T)) coefficient,
      compared against N2/N3 (is the plain flip == variant A/B at
      leading order?)
  B1  bridge: H == sum_{mu<nu} <F,F>_eta with <F,G>_eta = tr(eta F eta G^T)
      (the certified-lattice inner product), exactly, at the density level
  B2  bridge: H_flip == sum_{mu<nu} tr(F F^T), exactly

Out: ../data/m5_21_16_symbolic.json
"""
from __future__ import annotations

import json
import os
import time

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

# 24 real curvature symbols: G[mu][j] rotations, T[mu][j] boosts
G = [[sp.Symbol(f"G{m}_{j}", real=True) for j in (1, 2, 3)]
     for m in range(4)]
T = [[sp.Symbol(f"T{m}_{j}", real=True) for j in (1, 2, 3)]
     for m in range(4)]
g, delta = sp.symbols("g delta", positive=True)
XI = sp.diag(-1, 1, 1, 1)


def gamma_mu(m, boost_style):
    """the notebook's Gamma_mu. boost_style: 'real' | 'imag'."""
    t1, t2, t3 = T[m]
    if boost_style == "imag":
        row0 = [sp.I * t1, sp.I * t2, sp.I * t3]
        col0 = [-sp.I * t1, -sp.I * t2, -sp.I * t3]
    else:
        row0 = [t1, t2, t3]
        col0 = [t1, t2, t3]
    r1, r2, r3 = G[m]
    Gm = sp.zeros(4, 4)
    for j in range(3):
        Gm[0, 1 + j] = row0[j]
        Gm[1 + j, 0] = col0[j]
    Gm[1, 2], Gm[1, 3] = -r3, r2
    Gm[2, 1], Gm[2, 3] = r3, -r1
    Gm[3, 1], Gm[3, 2] = -r2, r1
    return Gm


def coms(A, B):
    return A * XI * B - B * XI * A


def build_F(boost_style, dvec, conj_second=False):
    d = sp.diag(*dvec)
    M = [coms(gamma_mu(m, boost_style), d) for m in range(4)]
    F = {}
    for mu in range(4):
        for nu in range(4):
            B = M[nu].conjugate() if conj_second else M[nu]
            F[(mu, nu)] = coms(M[mu], B)
    return F


def h_of(F, weights):
    """weights = (w_spatial, w_time): notebook H has (+1, -1); flip (+1, +1).
    Internal entries (0-indexed): spatial (1,2),(1,3),(2,3); time (0,1),(0,2),(0,3)."""
    ws, wt = weights
    H = 0
    for Fm in F.values():
        H += ws * (Fm[1, 2] ** 2 + Fm[1, 3] ** 2 + Fm[2, 3] ** 2)
        H += wt * (Fm[0, 1] ** 2 + Fm[0, 2] ** 2 + Fm[0, 3] ** 2)
    return sp.expand(H)


def omega2_leading(H):
    """coefficient of omega^2 (omega = G[0][0] i.e. G0_1), then g -> 1/delta
    leading order (through delta^0)."""
    om = G[0][0]
    c2 = sp.expand(H).coeff(om, 2)
    c2 = sp.expand(c2.subs(g, 1 / delta))
    # Laurent series through delta^0
    lead = sp.series(c2, delta, 0, 1).removeO()
    return sp.expand(lead)


SIX = sum(T[m][j] ** 2 for m in (1, 2, 3) for j in (1, 2))  # tG_i^{2,3}
T0 = T[0][1] ** 2 + T[0][2] ** 2                            # tG_0^{2,3}


def main():
    t0 = time.time()
    out = {}
    dvec = (g, 1, delta, 0)

    # N1 baseline
    Fb = build_F("real", dvec)
    Hb = h_of(Fb, (1, -1))
    n1 = omega2_leading(Hb)
    out["N1_baseline_lead"] = str(n1)
    out["N1_matches_minus2six"] = bool(
        sp.simplify(n1 - (-2 * SIX)) == 0)

    # N2 variant A: imaginary tilde-Gamma
    Fa = build_F("imag", dvec)
    Ha = h_of(Fa, (1, -1))
    n2 = omega2_leading(Ha)
    out["N2_variantA_lead"] = str(n2)
    out["N2_matches_plus2six"] = bool(sp.simplify(n2 - (2 * SIX)) == 0)

    # N3 variant B: imaginary g + conjugated second factor
    Fc = build_F("real", (sp.I * g, 1, delta, 0), conj_second=True)
    Hc = h_of(Fc, (1, -1))
    n3 = omega2_leading(Hc)
    out["N3_variantB_lead"] = str(n3)
    tgtB = 2 * (2 * T0 + SIX)
    out["N3_matches_2_2T0_plus_six"] = bool(sp.simplify(n3 - tgtB) == 0)

    # N4 the plain FLIP: all-plus internal weights on the baseline build
    Hf = h_of(Fb, (1, 1))
    n4 = omega2_leading(Hf)
    out["N4_flip_lead"] = str(n4)
    out["N4_equals_variantA"] = bool(sp.simplify(n4 - n2) == 0)
    out["N4_equals_variantB"] = bool(sp.simplify(n4 - n3) == 0)
    out["N4_minus_variantB"] = str(sp.simplify(n4 - n3))

    # B1/B2 bridge lemmas at the density level (exact, no series)
    def inner_eta(F1, F2):
        return sp.trace(XI * F1 * XI * F2.T)

    Hb_eta = 0
    Hf_frob = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            Fm = Fb[(mu, nu)]
            Hb_eta += inner_eta(Fm, Fm)
            Hf_frob += sp.trace(Fm * Fm.T)
    out["B1_H_equals_sum_inner_eta"] = bool(
        sp.expand(Hb - Hb_eta) == 0)
    out["B2_Hflip_equals_sum_frobenius"] = bool(
        sp.expand(Hf - Hf_frob) == 0)

    # N4_equals_variantA / B are structural FINDINGS (which mechanism the
    # plain flip implements), not pass criteria: excluded from all_pass
    finding_keys = {"N4_equals_variantA", "N4_equals_variantB"}
    out["all_pass"] = all(v for k, v in out.items()
                          if isinstance(v, bool) and k not in finding_keys)
    out["runtime_s"] = round(time.time() - t0, 1)
    os.makedirs(DATA, exist_ok=True)
    with open(os.path.join(DATA, "m5_21_16_symbolic.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items()
                      if isinstance(v, (bool, float))}, indent=1))


if __name__ == "__main__":
    main()
