"""M6.1 arm (b): day-scale characterization of Ouroboros+Eli+Fable v4.

Verifies what is checkable in the v4 spec (Zenodo 21447590) as printed,
without touching the unspecified constraint C[phi]:

  V1  Euler-Lagrange equation of the printed Lagrangian (sympy): the EL
      system cannot be closed because C[phi] is never written; the open
      equation is  box phi + m^2 phi + (lam/6) phi^3 + g rho = lam dC/dphi.
  V2  The printed Hamiltonian density follows from the printed L by
      Legendre transform (sympy residual check).
  V3  Boundedness below of V_eff = 1/2 m^2 phi^2 + (lam/4!) phi^4 + g rho phi
      for lam > 0 and bounded rho: analytic minimum of the tilted quartic
      (cubic stationarity), numeric global minimum scan.
  V4  Dimensional analysis: [phi] = 1, [lam] = 0, [g] = 0 given [rho_n] = 3;
      E_int = g phi A carries dimension 1 (energy).  Note: a dimensionless
      quartic coupling is strictly renormalizable, NOT superrenormalizable
      by power counting (v11's quartic had [g] = -4); v4 Section 1 carries
      the legacy superrenormalizability claim without a v4 derivation.
  V5  Range arithmetic: hbar_c / m_phi = 428.9 fm ("~400 fm" as printed);
      a Yukawa of that range deviates from 1/r by only 0.23% at 1 fm and
      2.3% at 10 fm, so the massive-vs-massless distinction the constraint
      C[phi] is invoked to repair is a <= 2.3% effect across the whole
      r = 1-10 fm window the paper quotes; the constraint has leverage
      only at r >~ 400 fm.
  V6  A-linearity genericity: with a Woods-Saxon rho_n normalized to A and
      any phi with range >> R_nucleus, E_int = int g phi rho d3x tracks
      g phi(0) A to within ~2% for A = 1..208.  The A-linear ENERGY follows
      from the normalization int rho = A plus long-rangedness alone; it is
      equally true of a standard per-nucleon scalar Yukawa coupling.
  V7  The in-medium field value is unconstrained by the paper: the printed
      E_int ~= g phi A leaves phi free; the self-consistent tadpole minimum
      at nuclear density is phi_min = -(6 g rho / lam)^(1/3) (quartic
      regime), which depends on the numerically unpinned lam.
  V8  The mass-migration table: 0.460 / 0.618 / 1.033 MeV <-> 428.9 /
      319.3 / 191.0 fm.

Also documents (text facts, no computation): the lambda symbol collision
(quartic coupling AND constraint multiplier in one Lagrangian line), the
duplicated paragraphs and "[INSERT FIGURE 1 HERE]" drafting artifacts, and
that C[phi] appears only as prose + a one-line surface integral with the
undefined symbol delta_Sigma.

Headless; writes data/m6_1_v4_characterization.json and
plots/m6_1_v4_veff_range.png. Runtime: seconds.
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent.parent
OUT_JSON = HERE / "data" / "m6_1_v4_characterization.json"
OUT_PLOT = HERE / "plots" / "m6_1_v4_veff_range.png"
report = {}

HBAR_C = 197.3269804  # MeV fm
M_PHI = 0.460  # MeV (v4)
G_J = 0.0054  # per nucleon (July 5/8 papers; v4 leaves g unpinned)
RHO0_FM3 = 0.16  # fm^-3 (v4 Woods-Saxon central density)
RHO0_MEV3 = RHO0_FM3 * HBAR_C**3  # nucleon density in MeV^3

# ------------------------------------------------------- V1/V2: EL + Legendre
t, x, y, z = sp.symbols("t x y z", real=True)
m, lam, g = sp.symbols("m lambda g", positive=True)
phi = sp.Function("phi")(t, x, y, z)
rho = sp.Function("rho")(x, y, z)
V = sp.Rational(1, 2) * m**2 * phi**2 + lam / 24 * phi**4
dphi = [sp.diff(phi, c) for c in (t, x, y, z)]
L = (sp.Rational(1, 2) * (dphi[0]**2 - dphi[1]**2 - dphi[2]**2 - dphi[3]**2)
     - V - g * phi * rho)  # the lam*C[phi] term is unwritable: C never given
el = sum(sp.diff(sp.diff(L, sp.diff(phi, c)), c) for c in (t, x, y, z)) - sp.diff(L, phi)
box_phi = (sp.diff(phi, t, 2) - sp.diff(phi, x, 2)
           - sp.diff(phi, y, 2) - sp.diff(phi, z, 2))
el_expected = box_phi + m**2 * phi + lam / 6 * phi**3 + g * rho
report["V1_euler_lagrange"] = {
    "el_matches_box_plus_m2_plus_cubic_plus_source": sp.simplify(el - el_expected) == 0,
    "open_equation": "box phi + m^2 phi + (lam/6) phi^3 + g rho = lam dC/dphi",
    "closable": False,
    "reason": "C[phi] is never written out in v4; the EL system cannot be closed",
}
pi_sym = sp.diff(phi, t)
H_legendre = sp.expand(pi_sym * sp.diff(phi, t) - L)
H_printed = (sp.Rational(1, 2) * pi_sym**2
             + sp.Rational(1, 2) * (dphi[1]**2 + dphi[2]**2 + dphi[3]**2)
             + sp.Rational(1, 2) * m**2 * phi**2 + lam / 24 * phi**4 + g * phi * rho)
report["V2_hamiltonian"] = {
    "printed_H_matches_legendre_transform": sp.simplify(H_legendre - H_printed) == 0,
    "caveat": "verified for L WITHOUT the lam*C[phi] term (unwritable); v4 prints H without it too",
}

# --------------------------------------------------------- V3: boundedness
phi_s = sp.symbols("phi_r", real=True)
veff = (sp.Rational(1, 2) * m**2 * phi_s**2 + lam / 24 * phi_s**4
        + g * phi_s * sp.symbols("rho_c", positive=True))
lims = {str(pt): sp.limit(veff, phi_s, pt) for pt in (sp.oo, -sp.oo)}
mins = {}
for lam_v in (0.1, 1.0, 10.0):
    grid = np.linspace(-200, 200, 400001)
    v = 0.5 * M_PHI**2 * grid**2 + lam_v / 24 * grid**4 + G_J * RHO0_MEV3 * grid
    mins[f"lam_{lam_v}"] = {"phi_min_MeV": round(float(grid[v.argmin()]), 3),
                            "Veff_min_MeV4": round(float(v.min()), 1),
                            "analytic_cubic_root_MeV": round(float(-(6 * G_J * RHO0_MEV3 / lam_v) ** (1 / 3)), 3)}
report["V3_boundedness"] = {
    "limits_at_pm_infinity": {k: str(v) for k, v in lims.items()},
    "bounded_below": all(str(v) == "oo" for v in lims.values()),
    "assumptions": "lam > 0; rho_n bounded with compact support (nuclear profile)",
    "numeric_minima_at_rho0_gJ": mins,
    "verdict": "v4's boundedness claim VERIFIED under its stated assumptions",
}

# --------------------------------------------------- V4: dimensional analysis
report["V4_dimensions"] = {
    "phi": 1, "m_phi": 1, "lambda_quartic": 0, "g_given_rho_dim3": 0,
    "E_int_g_phi_A": 1,
    "consistent": True,
    "renormalizability_note": ("dimensionless quartic -> strictly renormalizable, "
                               "NOT superrenormalizable by power counting; v11's "
                               "quartic coupling had [g] = -4; v4 Section 1 carries "
                               "the legacy superrenormalizability claim with no v4 "
                               "derivation"),
}

# ------------------------------------------------------- V5: range arithmetic
yukawa_range = HBAR_C / M_PHI
dev = lambda r_fm: 1.0 - float(np.exp(-r_fm / yukawa_range))
report["V5_range"] = {
    "hbar_c_over_m_phi_fm": round(yukawa_range, 1),
    "paper_says": "~400 fm",
    "yukawa_vs_1_over_r_deviation_pct": {"r_1fm": round(dev(1) * 100, 3),
                                         "r_10fm": round(dev(10) * 100, 2),
                                         "r_100fm": round(dev(100) * 100, 1),
                                         "r_400fm": round(dev(400) * 100, 1)},
    "verdict": ("a 428.9 fm Yukawa IS 1/r to within 2.3% over the whole "
                "r = 1-10 fm window the paper quotes; the mass/range paradox "
                "the constraint C[phi] is invoked to resolve is a percent-level "
                "effect below ~40 fm and only becomes order-one at r >~ 400 fm"),
}

# ------------------------------------------------- V6: A-linearity genericity
def woods_saxon_Eint(A):
    """E_int = int g phi rho d3x with rho Woods-Saxon normalized to A."""
    R_A = 1.2 * A ** (1 / 3)  # fm
    a_ws = 0.5  # fm surface diffuseness
    r_grid = np.linspace(1e-4, R_A + 12 * a_ws, 40000)
    rho_r = 1.0 / (1.0 + np.exp((r_grid - R_A) / a_ws))
    norm = np.trapezoid(4 * np.pi * r_grid**2 * rho_r, r_grid)
    rho_r *= A / norm
    phi_r = np.exp(-r_grid / yukawa_range)  # ambient long-range field, phi(0)=1
    e_int = np.trapezoid(4 * np.pi * r_grid**2 * G_J * phi_r * rho_r, r_grid)
    return e_int, G_J * 1.0 * A  # vs g*phi(0)*A

rows = {}
for A in (1, 2, 40, 85, 184, 208):
    e_num, e_lin = woods_saxon_Eint(A)
    rows[f"A_{A}"] = {"E_int_numeric": round(e_num, 5), "g_phi0_A": round(e_lin, 5),
                      "ratio": round(e_num / e_lin, 4)}
report["V6_A_linearity"] = {
    "rows": rows,
    "max_deviation_pct": round(max(abs(1 - v["ratio"]) for v in rows.values()) * 100, 2),
    "verdict": ("E_int proportional to A holds to ~2% for ANY field with range >> "
                "R_nucleus once rho_n is normalized to A: the v4 'derivation' is "
                "the normalization identity int rho = A plus long-rangedness, and "
                "the A-linear ENERGY is equally a property of a standard "
                "per-nucleon scalar Yukawa coupling; distinguishing v4 substance "
                "lives in cross-sections, which v4 defers ('ongoing work')"),
}

# -------------------------------------------- V7: the unpinned in-medium phi
report["V7_in_medium_phi"] = {
    "tadpole_minimum_formula": "phi_min = -(6 g rho / lam)^(1/3) in the quartic regime",
    "phi_min_MeV_at_lam_1": round(-(6 * G_J * RHO0_MEV3 / 1.0) ** (1 / 3), 2),
    "lam_dependence": "phi_min ~ lam^(-1/3); lam is not pinned numerically in v4",
    "verdict": ("the phi entering E_int ~= g phi A is fixed by nothing in the "
                "paper: no equation, no value; with the July g_J and lam = 1 the "
                "self-consistent in-medium value would be ~ -34 MeV, 74x m_phi"),
}

# ------------------------------------------------------- V8: mass migration
report["V8_mass_migration"] = {
    f"{mv}_MeV": f"{HBAR_C / mv:.1f} fm"
    for mv in (0.460, 0.618, 1.033)
}

# ------------------------------------------------------------- text facts
report["T1_lambda_collision"] = ("lambda is the quartic self-coupling in V(phi) AND "
                                 "the multiplier of C[phi] in the same Lagrangian line")
report["T2_C_phi_unspecified"] = ("C[phi] appears only as prose + the one-line surface "
                                  "integral I = int phi delta_Sigma - Sigma_n with "
                                  "delta_Sigma undefined; no functional form anywhere")
report["T3_drafting_artifacts"] = ("duplicated paragraphs (Sections 4.1/5), '[INSERT "
                                   "FIGURE 1 HERE]' placeholders, 'In Word, insert each "
                                   "figure image from your /OuroborosEli_figures folder'")

# ------------------------------------------------------------------- plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
grid = np.linspace(-75, 75, 3000)
for lam_v, style in ((1.0, "-"), (10.0, "--")):
    for rho_v, color, lab in ((0.0, "tab:blue", "vacuum"),
                              (RHO0_MEV3, "tab:red", "nuclear rho_0")):
        v = 0.5 * M_PHI**2 * grid**2 + lam_v / 24 * grid**4 + G_J * rho_v * grid
        ax1.plot(grid, v, style, color=color, label=f"lam={lam_v:g}, {lab}", lw=1.4)
        if rho_v > 0:
            ax1.plot(grid[v.argmin()], v.min(), "o", color=color, ms=5)
ax1.set_xlabel("phi [MeV]"); ax1.set_ylabel("V_eff [MeV^4]")
ax1.set_title("V4 effective potential: bounded below, tilted in medium\n"
              f"(m_phi = {M_PHI} MeV, g = {G_J}, rho_0 = 0.16 fm^-3; dots = in-medium minima)")
ax1.legend(fontsize=8); ax1.grid(alpha=0.3)
ax1.set_ylim(-2.4e5, 2.5e5)

r_fm = np.logspace(-1, 3.5, 500)
ax2.semilogx(r_fm, 100 * (1 - np.exp(-r_fm / yukawa_range)), lw=1.6)
for rv, lab in ((1, "1 fm"), (10, "10 fm"), (400, "400 fm")):
    ax2.axvline(rv, color="gray", ls=":", lw=0.8)
    ax2.annotate(lab, (rv, 60), rotation=90, fontsize=8, color="gray")
ax2.set_xlabel("r [fm]"); ax2.set_ylabel("Yukawa deviation from 1/r [%]")
ax2.set_title("How different is the 428.9 fm Yukawa from 1/r?\n"
              "(the C[phi] 'mass/range paradox' is percent-level below ~40 fm)")
ax2.grid(alpha=0.3, which="both")
fig.tight_layout()
OUT_PLOT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_PLOT, dpi=150)

OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
OUT_JSON.write_text(json.dumps(report, indent=2))
print(json.dumps(report, indent=2))
print(f"\nwritten: {OUT_JSON}\nplot: {OUT_PLOT}")
