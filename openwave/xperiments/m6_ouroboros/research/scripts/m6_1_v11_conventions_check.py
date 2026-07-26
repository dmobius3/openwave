"""M6.1 arm (a): machine checks behind the v11 convention sheet.

Verifies, from the LoE v11 paper's own printed content (Zenodo 20357670),
every arithmetic and symbolic claim the convention sheet
(research/m6_1_v11_convention_sheet.md) makes:

  C1  e_nat = sqrt(4 pi alpha) reproduces the paper's 0.30282, and the
      physical target H/Q = m_e / e_nat = 0.511 / 0.30282 = 1.6875.
  C2  Gap arithmetic between every pair of the four circulating H/Q values
      (1.6875 physical, 1.6890 paper scan, 1.6918 reproduction, 1.6969
      relabeled target): which printed "gap" labels are arithmetically
      consistent with which target.
  C3  The Section 8 Step-3 "predicted" charge Q = 0.3011 equals
      m_e / 1.6969 (and 0.511 / 1.6890 = 0.30255): the charge prediction
      is the H/Q discrepancy restated, not an independent number.
  C4  The R_phys = 191 fm scale-setting: which H_code value it implies
      (H_code = 0.511 * R_phys / hbar_c), and the numerical identity
      hbar_c / 191 fm = 1.033 MeV (the value later published as m_J).
  C5  The printed Section 5.1 ansatz A = r_hat x grad(phi(r)) cos(wt) is
      identically zero for radial phi (symbolic).
  C6  The printed Euler-Lagrange pair (2.1)/(2.2) is NOT the EL system of
      the printed Lagrangian -F^2 - G^2 + J.A - g(J.J)^2 under either
      metric signature; the unique reconstruction preserving the printed
      interaction signs is -1/4 F^2 - 1/4 G^2 + J.A - g(J.J)^2 under
      mostly-plus signature (symbolic, residual-verified).

Headless; writes data/m6_1_v11_conventions.json. Runtime: ~1 min.
"""

import json
from pathlib import Path

import numpy as np
import sympy as sp

OUT = Path(__file__).resolve().parent.parent / "data" / "m6_1_v11_conventions.json"
report = {}

# ---------------------------------------------------------------- C1: e_nat
ALPHA_INV = 137.035999084  # CODATA fine-structure constant, 1/alpha
M_E = 0.51099895  # MeV
e_nat = float(np.sqrt(4 * np.pi / ALPHA_INV))
target_phys = M_E / e_nat
report["C1_e_nat"] = {
    "e_nat_sqrt_4pi_alpha": round(e_nat, 6),
    "paper_value": 0.30282,
    "match": abs(e_nat - 0.30282) < 5e-5,
    "target_phys_me_over_enat": round(target_phys, 6),
    "paper_target": 1.6875,
    "match_target": abs(target_phys - 1.6875) < 5e-4,
}

# ------------------------------------------------- C2: the four-value table
vals = {"1.6875_physical": 1.6875, "1.6890_scan": 1.6890,
        "1.6918_reproduction": 1.6918, "1.6969_relabeled": 1.6969}
gaps = {}
for na, va in vals.items():
    for nb, vb in vals.items():
        if na < nb:
            gaps[f"{na}_vs_{nb}"] = round(abs(va - vb) / vb * 100, 3)
report["C2_gap_matrix_pct"] = gaps
report["C2_verdicts"] = {
    "scan_1.6890_vs_physical": f"{gaps['1.6875_physical_vs_1.6890_scan']}% (paper prints 0.09%)",
    "repro_1.6918_vs_relabeled_1.6969": f"{gaps['1.6918_reproduction_vs_1.6969_relabeled']}% (paper 5.1 prints 0.30%)",
    "repro_1.6969_vs_physical_1.6875": (
        f"{round(abs(1.6969 - 1.6875) / 1.6875 * 100, 3)}% "
        "(paper Table 2 prints this pairing as 0.30%: arithmetically inconsistent; "
        "0.30% only holds against the model-internal 1.6969)"),
}

# --------------------------------------------- C3: the Q_predicted identity
q_pred_from_1p6969 = M_E / 1.6969
q_pred_from_1p6890 = M_E / 1.6890
report["C3_Q_predicted_identity"] = {
    "me_over_1.6969": round(q_pred_from_1p6969, 5),
    "paper_Q_predicted": 0.3011,
    "identity_holds": abs(q_pred_from_1p6969 - 0.3011) < 5e-4,
    "me_over_1.6890": round(q_pred_from_1p6890, 5),
    "gap_vs_e_nat_pct_1.6969": round(abs(q_pred_from_1p6969 - e_nat) / e_nat * 100, 3),
    "gap_vs_e_nat_pct_1.6890": round(abs(q_pred_from_1p6890 - e_nat) / e_nat * 100, 3),
    "verdict": ("Q_predicted = m_e/(H/Q): its 0.56% gap vs e_nat is the H/Q gap vs the "
                "physical target restated, not an independent prediction"),
}

# ------------------------------------------------- C4: R_phys scale-setting
HBAR_C = 197.3269804  # MeV fm
h_code_implied = M_E * 191.0 / HBAR_C
unit_energy_191 = HBAR_C / 191.0
report["C4_R_phys"] = {
    "H_code_implied_by_191fm": round(h_code_implied, 5),
    "printed_anywhere_in_v11": False,
    "unit_energy_hbar_c_over_191fm_MeV": round(unit_energy_191, 5),
    "july_papers_m_J_MeV": 1.033,
    "numeric_identity": abs(unit_energy_191 - 1.033) < 5e-3,
    "note": ("191 fm requires H_code = 0.4946 (equivalently Q_code = 0.2928 at "
             "H/Q = 1.6890), neither printed in v11; hbar_c/191 fm = 1.0331 MeV "
             "coincides numerically with the m_J = 1.033 MeV of the July 5/8 papers"),
}

# ------------------------------------- C5: the printed ansatz is degenerate
x, y, z = sp.symbols("x y z", real=True)
r = sp.sqrt(x**2 + y**2 + z**2)
phi = sp.Function("phi")
grad_phi = sp.Matrix([sp.diff(phi(r), c) for c in (x, y, z)])
r_hat = sp.Matrix([x, y, z]) / r
cross = r_hat.cross(grad_phi)
cross_zero = all(sp.simplify(c) == 0 for c in cross)
report["C5_ansatz_degenerate"] = {
    "r_hat_cross_grad_phi_of_r_is_zero": bool(cross_zero),
    "verdict": ("the Section 5.1 production ansatz as printed vanishes identically "
                "for radial phi(r); the actual field configuration exists only in "
                "the benchmark code" if cross_zero else "NOT degenerate: recheck"),
}

# ------------------------------- C6: EL pair vs Lagrangian normalization
# General Lagrangian  L = -a F^2 - a G^2 + c J.A - q (J.J)^2, fields as
# upper-index components, box = d_t^2 - lap as printed in v11 Table 1.
# For each metric signature (the paper never states one), derive the EL
# equations symbolically, verify the reduced-form identity by residual,
# then solve in Lorenz gauge for the (a, c, q) reproducing the printed pair
#   (2.1) box A^nu = J^nu     (2.2) box J^nu = A^nu - 4 g (J.J) J^nu.
t = sp.symbols("t", real=True)
coords = (t, x, y, z)
box = lambda f: (sp.diff(f, t, 2) - sp.diff(f, x, 2)
                 - sp.diff(f, y, 2) - sp.diff(f, z, 2))
c6 = {}
for sig_name, diag_sig in [("mostly_minus(+,-,-,-)", (1, -1, -1, -1)),
                           ("mostly_plus(-,+,+,+)", (-1, 1, 1, 1))]:
    eta = sp.diag(*diag_sig)
    A = [sp.Function(f"A{i}")(*coords) for i in range(4)]  # upper-index comps
    J = [sp.Function(f"J{i}")(*coords) for i in range(4)]
    a_c, c_c, q_c = sp.symbols("a c q", positive=True)
    lower = lambda V: [sum(eta[m, n] * V[n] for n in range(4)) for m in range(4)]
    d = lambda mu, f: sp.diff(f, coords[mu])
    A_lo, J_lo = lower(A), lower(J)
    F_lo = [[d(m, A_lo[n]) - d(n, A_lo[m]) for n in range(4)] for m in range(4)]
    G_lo = [[d(m, J_lo[n]) - d(n, J_lo[m]) for n in range(4)] for m in range(4)]
    raise2 = lambda T: [[sum(eta[m, mm] * eta[n, nn] * T[mm][nn]
                             for mm in range(4) for nn in range(4))
                         for n in range(4)] for m in range(4)]
    F_up, G_up = raise2(F_lo), raise2(G_lo)
    F2 = sum(F_up[m][n] * F_lo[m][n] for m in range(4) for n in range(4))
    G2 = sum(G_up[m][n] * G_lo[m][n] for m in range(4) for n in range(4))
    JA = sum(J[m] * A_lo[m] for m in range(4))
    JJ = sum(J[m] * J_lo[m] for m in range(4))
    L = -a_c * F2 - a_c * G2 + c_c * JA - q_c * JJ**2
    nu = 1  # representative spatial component
    el_A = sum(sp.diff(sp.diff(L, sp.diff(A[nu], coords[m])), coords[m])
               for m in range(4)) - sp.diff(L, A[nu])
    el_J = sum(sp.diff(sp.diff(L, sp.diff(J[nu], coords[m])), coords[m])
               for m in range(4)) - sp.diff(L, J[nu])
    divA = sum(d(m, A[m]) for m in range(4))
    divJ = sum(d(m, J[m]) for m in range(4))
    s2 = 1 if diag_sig[0] == 1 else -1  # source sign, verified by residual
    idA = sp.simplify(el_A - (4 * a_c * (box(A[nu]) + sp.diff(divA, x))
                              + s2 * c_c * J[nu])) == 0
    idJ = sp.simplify(el_J - (4 * a_c * (box(J[nu]) + sp.diff(divJ, x))
                              + s2 * c_c * A[nu] - s2 * 4 * q_c * JJ * J[nu])) == 0
    # Lorenz gauge (div = 0):  box A = -(s2)(c/4a) J
    #                          box J = -(s2)(c/4a) A + (s2)(q/a)(J.J) J
    # printed pair needs:  -(s2) c/4a = +1   and   (s2) q/a = -4g
    g_val = sp.symbols("g", positive=True)
    a_s, c_s, q_s = sp.symbols("a_s c_s q_s")
    sol = sp.solve([sp.Eq(-s2 * c_s / (4 * a_s), 1),
                    sp.Eq(s2 * q_s / a_s, -4 * g_val)], [c_s, q_s], dict=True)[0]
    c6[sig_name] = {
        "identity_A_verified": bool(idA),
        "identity_J_verified": bool(idJ),
        "coeffs_reproducing_printed_pair_at_a_quarter": {
            "a": "1/4 (kinetic normalization; printed: 1)",
            "c": str(sp.simplify(sol[c_s].subs(a_s, sp.Rational(1, 4)))) + " (printed: +1)",
            "q": str(sp.simplify(sol[q_s].subs(a_s, sp.Rational(1, 4)))) + " (printed: g)",
        },
    }
report["C6_EL_normalization"] = {
    "per_signature": c6,
    "el_of_printed_lagrangian_in_lorenz_gauge": {
        "mostly_plus": "box A = J/4, box J = A/4 - g(J.J)J",
        "mostly_minus": "box A = -J/4, box J = -A/4 + g(J.J)J",
    },
    "verdict": ("the printed EL pair (2.1)/(2.2) is NOT the EL system of the "
                "printed Lagrangian -F^2 - G^2 + J.A - g(J.J)^2 under either "
                "signature; the unique reconstruction preserving the printed "
                "interaction signs is L_ref = -1/4 F^2 - 1/4 G^2 + J.A - "
                "g(J.J)^2 under mostly-plus signature (-,+,+,+); under "
                "(+,-,-,-) the coupling and quartic signs must flip as well; "
                "v11 states no metric signature; the dropped 1/4 cannot be "
                "absorbed by rescaling L because the quartic coefficient -4g "
                "in (2.2) pins the scale"),
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
print(f"\nwritten: {OUT}")
