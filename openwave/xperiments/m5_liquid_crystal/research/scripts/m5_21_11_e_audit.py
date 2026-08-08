"""M5.21.11 adversarial audit: independent numeric checks against
findings/m5_21_11_framework.md section 1 derivations and the section
2/3/6 arithmetic. Written by the audit agent with its own algebra;
the instrument module is imported ONLY to evaluate the functional of
record (so the checks bind to the implementation, not to the doc).

Checks:
  A1 Derrick: E_u ~ s^-1 and E_V ~ s^3 under spatial dilation
  A2 far field: commutator identity ([G,D])_ab = G_ab (d_b - d_a);
     dM = O [Gamma, D] O^T on a path; E_u(delta) for an O(x)-texture
     is EXACTLY a quartic polynomial in delta, finite at delta = 0,
     with a nonzero linear coefficient
  A3 line-core competition exponents (and the point-core variant),
     plus the dimensional-analysis premises
  A4 artanh identity and series
  A5 evenness counterexample: two even families with the SAME
     m*(g) = artanh(1/g) but gain ~ m*^2 vs gain ~ m*^4; shows the
     doc's step is not a derivation
  A6 g-freeness of the 3x3 instrument (config keys + source scan)
  A7 ladder statistics arithmetic (points, params, df)
  A8 internal numbers: a(0.05), L/a, delta_phys^0.4, nu refit from
     the quoted radii, a*(delta_phys) vs the box, F2 false-alarm rate

Output: ../data/m5_21_11_audit.json
ASCII only; no em dashes.
"""
from __future__ import annotations

import importlib.util
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
INSTR_PATH = os.path.join(HERE, "m5_21_2b_a_instrument.py")

spec = importlib.util.spec_from_file_location("instr", INSTR_PATH)
instr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(instr)

rng = np.random.default_rng(7)
R: dict = {}


def smooth_field(n, h, shape, n_modes=5, amp=0.4):
    """sum of low-frequency sinusoids: smooth on the grid."""
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    L = n * h
    out = np.zeros((n, n, n) + shape)
    for _ in range(n_modes):
        k = rng.integers(-2, 3, size=3) * (2 * np.pi / L)
        ph = rng.uniform(0, 2 * np.pi)
        c = amp * rng.normal(size=shape)
        base = np.sin(k[0] * X + k[1] * Y + k[2] * Z + ph)
        out += base.reshape(base.shape + (1,) * len(shape)) * c
    return out


def rodrigues(w):
    """SO(3) exponential of the axis-angle field w (..., 3)."""
    th = np.linalg.norm(w, axis=-1, keepdims=True)
    th = np.where(th < 1e-12, 1e-12, th)
    k = w / th
    th = th[..., 0]
    K = np.zeros(w.shape[:-1] + (3, 3))
    K[..., 0, 1] = -k[..., 2]
    K[..., 0, 2] = k[..., 1]
    K[..., 1, 0] = k[..., 2]
    K[..., 1, 2] = -k[..., 0]
    K[..., 2, 0] = -k[..., 1]
    K[..., 2, 1] = k[..., 0]
    st = np.sin(th)[..., None, None]
    ct = (1.0 - np.cos(th))[..., None, None]
    return np.eye(3) + st * K + ct * (K @ K)


# ---------------- A1: Derrick scaling exponents ----------------
cfg = instr.base_cfg(term="T2", eps=0.0, n=20, L=10.0, delta=0.3,
                     w2=1.0)
Mr = smooth_field(20, cfg["h"], (3, 3))
Mr = 0.5 * (Mr + Mr.swapaxes(-1, -2)) + np.diag([1.0, 0.3, 0.0])
s = 2.0
eu1, _, ev1 = instr.e_parts(Mr, cfg, "sym")
cfg2 = dict(cfg)
cfg2["h"] = cfg["h"] * s
cfg2["L"] = cfg["L"] * s
eu2, _, ev2 = instr.e_parts(Mr, cfg2, "sym")
exp_u = float(np.log(eu2 / eu1) / np.log(s))
exp_v = float(np.log(ev2 / ev1) / np.log(s))
# stationarity algebra: d/ds (Eu/s + s^3 Ev) at s=1 -> -Eu + 3Ev = 0
R["A1_derrick"] = {
    "exp_u_measured": exp_u, "exp_u_expected": -1.0,
    "exp_v_measured": exp_v, "exp_v_expected": 3.0,
    "algebra": "d/ds(Eu/s + s^3 Ev)|1 = -Eu + 3Ev = 0 => Eu = 3Ev; "
               "E = Eu + Ev = 4Ev = (4/3)Eu",
    "pass": abs(exp_u + 1.0) < 1e-10 and abs(exp_v - 3.0) < 1e-10,
    "caveat": "exact only for a boundary-free dilation family; the "
              "pinned finite box breaks it (doc measures +0.034 "
              "residual and gates on 0.05, consistent)"}

# ---------------- A2: far-field analyticity ----------------
# commutator identity, exact
G = rng.normal(size=(3, 3))
G = G - G.T
dvec = rng.normal(size=3)
Dm = np.diag(dvec)
comm = G @ Dm - Dm @ G
expect = G * (dvec[None, :] - dvec[:, None])
id_err = float(np.max(np.abs(comm - expect)))

# dM = O [Gamma, D] O^T on a one-parameter path (central difference)
A = rng.normal(size=(3, 3))
A = A - A.T
avec = np.array([A[2, 1], A[0, 2], A[1, 0]])
O0 = rodrigues(rng.normal(size=(1, 3)))[0]
D0 = np.diag([1.0, 0.3, 0.0])
t = 1e-6


def M_of(tt):
    Ot = O0 @ rodrigues((tt * avec)[None, :])[0]
    return Ot @ D0 @ Ot.T


dM_num = (M_of(t) - M_of(-t)) / (2 * t)
dM_an = O0 @ (A @ D0 - D0 @ A) @ O0.T
path_err = float(np.max(np.abs(dM_num - dM_an))
                 / np.max(np.abs(dM_an)))

# E_u(delta) for a pure O(x) D O(x)^T texture: exact quartic in delta
n2, L2 = 20, 10.0
cfg_ff = instr.base_cfg(term="T2", eps=0.0, n=n2, L=L2, delta=0.3,
                        w2=1.0)
wax = smooth_field(n2, cfg_ff["h"], (3,), n_modes=5, amp=0.6)
Ofield = rodrigues(wax)
deltas = np.linspace(0.0, 0.4, 9)
eus, evs = [], []
for d in deltas:
    Dd = np.diag([1.0, float(d), 0.0])
    Mff = np.einsum("...ab,bc,...dc->...ad", Ofield, Dd, Ofield)
    c = dict(cfg_ff)
    c["delta"] = float(d)
    eu, _, ev = instr.e_parts(Mff, c, "sym")
    eus.append(float(eu))
    evs.append(float(ev))
coef = np.polyfit(deltas, eus, 4)
resid = float(np.max(np.abs(np.polyval(coef, deltas) - eus))
              / np.max(np.abs(eus)))
p0, p1 = float(coef[-1]), float(coef[-2])
R["A2_far_field"] = {
    "commutator_identity_max_err": id_err,
    "dM_path_rel_err": path_err,
    "E_u_at_delta0": eus[0],
    "quartic_fit_rel_resid": resid,
    "c1_linear_over_const": p1 / p0,
    "poly_coeffs_desc": [float(x) for x in coef],
    "E_V_max_over_deltas": float(np.max(np.abs(evs))),
    "pass": (id_err < 1e-12 and path_err < 1e-6 and resid < 1e-9
             and eus[0] > 0.0 and abs(p1 / p0) > 1e-3),
    "note": "M is linear in delta, so lattice E_u is EXACTLY a "
            "quartic polynomial in delta; finite nonzero delta->0 "
            "limit; linear coefficient generically nonzero. E_V = 0 "
            "for the pure texture (eigenvalues equal targets)."}

# ---------------- A3: core competition exponents ----------------
s_exp, ku, kv, wgt = 0.8, 1.7, 0.9, 1.3
dgrid = np.logspace(-6, -1, 6)


def min_scan(f):
    a = np.logspace(-3, 5, 40001)
    e = f(a)
    i = int(np.argmin(e))
    a2 = np.linspace(a[max(i - 1, 0)], a[min(i + 1, len(a) - 1)],
                     20001)
    e2 = f(a2)
    j = int(np.argmin(e2))
    return a2[j], e2[j]


astars, estars, astp, estp = [], [], [], []
for d in dgrid:
    pen = kv * wgt * d ** s_exp
    aL, eL = min_scan(lambda a: ku / a ** 2 + pen * a ** 2)
    aP, eP = min_scan(lambda a: ku / a + pen * a ** 3)
    astars.append(aL)
    estars.append(eL)
    astp.append(aP)
    estp.append(eP)
ld = np.log(dgrid)
sl_a_line = float(np.polyfit(ld, np.log(astars), 1)[0])
sl_e_line = float(np.polyfit(ld, np.log(estars), 1)[0])
sl_a_pt = float(np.polyfit(ld, np.log(astp), 1)[0])
sl_e_pt = float(np.polyfit(ld, np.log(estp), 1)[0])
R["A3_core_competition"] = {
    "s_exp_used": s_exp,
    "line_a_slope": sl_a_line, "line_a_expected": -s_exp / 4,
    "line_e_slope": sl_e_line, "line_e_expected": s_exp / 2,
    "point_a_slope": sl_a_pt, "point_a_expected": -s_exp / 4,
    "point_e_slope": sl_e_pt, "point_e_expected": s_exp / 4,
    "theta_eq_2nu_line": "theta = s/2 = 2*(s/4) = 2*nu confirmed",
    "premises": "quartic-in-gradient density (Delta/a)^4 over a "
                "cross-section a^2 gives ku/a^2 per unit line "
                "length; over a volume a^3 gives ku/a for a point; "
                "T2 penalty w*delta^s over a^2 (line) or a^3 "
                "(point): both premises dimensionally correct",
    "pass": (abs(sl_a_line + s_exp / 4) < 1e-3
             and abs(sl_e_line - s_exp / 2) < 1e-3
             and abs(sl_a_pt + s_exp / 4) < 1e-3
             and abs(sl_e_pt - s_exp / 4) < 1e-3)}

# ---------------- A4: artanh identity ----------------
gg = np.array([2.0, 8.0, 16.0, 32.0, 64.0, 1e3, 1e6])
lhs = 0.5 * np.log((gg + 1) / (gg - 1))
rhs = np.arctanh(1.0 / gg)
ser = 1.0 / gg + 1.0 / (3 * gg ** 3)
# float note: at g = 1e6 the ln form loses ~5 digits to
# cancellation; the identity itself is exact. Gate on g <= 64.
rel = np.abs(lhs - rhs) / np.abs(rhs)
R["A4_artanh"] = {
    "identity_max_rel_err_g_le_64": float(np.max(rel[:5])),
    "identity_rel_err_g1e6_float_cancellation": float(rel[6]),
    "series_rel_err_at_g8": float(abs(rhs[1] - ser[1]) / rhs[1]),
    "series_next_term_over_at_g8": float((1.0 / (5 * 8.0 ** 5))
                                         / rhs[1]),
    "pass": bool(np.max(rel[:5]) < 1e-14)}

# ---------------- A5: evenness does not pin the gain power ------
# Two even families, both with minima exactly at +-m*(g),
# m*(g) = artanh(1/g):
#   family A (stiffening quartic): E = -m^2 + m^4/(2 m*^2)
#       gain = -m*^2/2                (kappa fixed, doc's reading)
#   family B (softening quadratic, pitchfork): E = -2 m*^2 m^2 + m^4
#       gain = -m*^4                  (kappa NOT fixed; g^-4 law)
gg5 = np.array([64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0])
ms = np.arctanh(1.0 / gg5)
mgrid = np.linspace(-0.2, 0.2, 400001)
gainA, gainB, locA, locB = [], [], [], []
for m0 in ms:
    EA = -mgrid ** 2 + mgrid ** 4 / (2 * m0 ** 2)
    EB = -2 * m0 ** 2 * mgrid ** 2 + mgrid ** 4
    iA, iB = int(np.argmin(EA)), int(np.argmin(EB))
    locA.append(abs(abs(mgrid[iA]) / m0 - 1.0))
    locB.append(abs(abs(mgrid[iB]) / m0 - 1.0))
    gainA.append(EA[iA])
    gainB.append(EB[iB])
slA = float(np.polyfit(np.log(gg5), np.log(-np.array(gainA)), 1)[0])
slB = float(np.polyfit(np.log(gg5), np.log(-np.array(gainB)), 1)[0])
R["A5_evenness_gain"] = {
    "minima_at_mstar_rel_err_max": float(max(max(locA), max(locB))),
    "gain_slope_family_A_vs_g": slA, "expected_A": -2.0,
    "gain_slope_family_B_vs_g": slB, "expected_B": -4.0,
    "verdict": "evenness + minima at m*(g)->0 do NOT imply gain = "
               "-kappa*m*^2 with fixed kappa: family B (quadratic "
               "coefficient softening, a pitchfork) satisfies the "
               "same empirical inputs and gives gain ~ -m*^4 ~ "
               "g^-4. The doc's O(m*^2) gain law is an extra "
               "assumption (fixed E''(0) < 0 as g varies), not a "
               "consequence; it must come from the measured g-arm.",
    "consequence": "separability SURVIVES either way (g^-4 is even "
                   "smaller than g^-2 at g ~ 1e10) but F4 as "
                   "written could terminally fail route (b) on a "
                   "healthy instrument if the true law is g^-4",
    "pass_arithmetic": abs(slA + 2) < 0.02 and abs(slB + 4) < 0.02}

# ---------------- A6: g-freeness of the 3x3 instrument ----------
cfg_keys = sorted(instr.base_cfg().keys())
src = open(INSTR_PATH).read()
R["A6_g_freeness"] = {
    "base_cfg_keys": cfg_keys,
    "has_g_key": "g" in cfg_keys,
    "source_mentions_artanh_or_boostparam": ("artanh" in src
                                             or "g_phys" in src),
    "pass": "g" not in cfg_keys and "artanh" not in src,
    "note": "T2 potential = w2 * sum_k (lam_k - v_k)^2 with "
            "v = sorted(1, delta, 0); curvature term parameter-"
            "free; no g anywhere in the instrument"}

# ---------------- A7: ladder statistics ----------------
points = 8 * 3
holdout = 2 * 3
fitpts = points - holdout
params = 1 + 3 * 3
df = fitpts - params
df_fallback = fitpts - 3 * 2
R["A7_statistics"] = {
    "points": points, "holdout": holdout, "fit_points": fitpts,
    "params_shared_theta": params, "df": df, "df_claimed": 8,
    "df_after_degeneracy_fallback": df_fallback,
    "refinement_runs": 3 * 3 * 3, "beyond_production": 27 - 9,
    "pass": df == 8 and fitpts == 18 and params == 10,
    "determinism_gap": "the theta 68pct interval METHOD (delta-"
                       "method vs profile) and whether the fit "
                       "enforces theta <= 1 are not pinned; both "
                       "affect whether the interval 'includes 1'"}

# ---------------- A8: internal numbers ----------------
a005 = 2.5 * (0.05 / 0.3) ** (-0.2)
rad_d = np.array([0.3, 0.2, 0.1])
rad_a = np.array([2.5, 2.9, 3.4])
nu_refit = -float(np.polyfit(np.log(rad_d), np.log(rad_a), 1)[0])
a005_refit = 2.5 * (0.05 / 0.3) ** (-nu_refit)
a_phys = 2.5 * (1e-10 / 0.3) ** (-0.2)
p = 0.05
f2_false_alarm = 1.0 - (1 - p) ** 6 - 6 * p * (1 - p) ** 5
R["A8_internal_numbers"] = {
    "a_005_nu02": float(a005), "doc_value": 3.6,
    "L_over_a_nu02": float(48 / a005), "doc_value_La": 13,
    "nu_refit_from_quoted_radii": nu_refit,
    "a_005_nu_refit": float(a005_refit),
    "L_over_a_nu_refit": float(48 / a005_refit),
    "delta_phys_pow_04": float((1e-10) ** 0.4),
    "artanh_sq_at_gphys": float(np.arctanh(1e-10) ** 2),
    "a_star_at_delta_phys_nu02": float(a_phys),
    "box_L": 48.0,
    "g_from_1e38_quarter": float(1e38 ** 0.25),
    "F2_false_alarm_prob_2of6_at_95pct": float(f2_false_alarm),
    "pass_arithmetic": abs(a005 - 3.578) < 0.01
                       and abs((1e-10) ** 0.4 - 1e-4) < 1e-18,
    "catches": [
        "nu point estimate from the quoted radii is 0.27, not 0.2 "
        "(inside the stated 0.2 +/- 0.1 band, but the central "
        "value the doc propagates is low)",
        "a*(delta_phys) ~ 2.0e2 >> L = 48: the physical point is "
        "8+ decades of delta beyond the measured window and the "
        "core does not fit any feasible box there; the prediction "
        "is pure functional-form extrapolation (doc partly "
        "concedes via 'numerically the E_inf limit')",
        "F2 has a ~3.3 percent false-alarm probability from the "
        ">=2-of-6-at-2sigma clause alone if intervals are "
        "calibrated"]}

os.makedirs(DATA, exist_ok=True)
with open(os.path.join(DATA, "m5_21_11_audit.json"), "w") as f:
    json.dump(R, f, indent=1)
print(json.dumps(R, indent=1))
