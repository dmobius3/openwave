"""M5.32 R29-3: the three cosmology add-ons (the author's post of 2026-10-22), at the form level.

EQUATIONS FIRST
---------------
Vacuum N = M_vac eta with eigenvalues lambda = (-g, 1, delta, 0) (toy) and (-g, 1, delta, delta) (R23); the four
eigenvalue modes are the diagonal directions diag_a (mode a scales lambda_a).
(1) The gap map. The author's combination c = 16 pi Delta^4 has no definition of Delta on this side (the reply left
    it as a placeholder), so the add-on reports, along each eigenvalue mode a, the logarithmic derivatives of every
    eigenvalue gap, d ln |lambda_i - lambda_j| / d ln lambda_a (the only directions the map can be built on once Delta
    is given as a function of the gaps), and beside them the author's own d ln(M_N / m_pi) / d ln g chain numbers,
    labelled as the author's (3.76 time-axis only, 3.29 all four, 0 projector form).
(2) The energy density of the degenerate point: V4 at four equal eigenvalues, minimised over the common value,
        rho_deg = min_lambda w sum_p (4 lambda^p - C_p)^2,
    in model units (energy per h^3 cell is V4 x h^3 on the lattice; the density is V4 itself), converted with the
    OpenWave anchor (1 unit = 0.0305 MeV from M5.8 P1; the length unit from the same record's r0, see the json) and
    reported as not expressible in the author's anchor (R29-0 (f)).
(3) The sign of the boosts' kinetic norm per term: the coefficient of eps_dot^2 for X = boost0j about M_vac, from
    R29-0 (a) for -4 I1, K_T, K_cov, K_h, K_spec, and for the torsion row the omega^2 coefficient of the record boost
    direction dE = E0 K_j (E -> E Lambda^T) at the scan points, from the A2 block of R29-0's torsion_qep.
Run from this directory: python3 m5_32_r29_3_addons.py -> ../data/m5_32_r29_3_addons.json (seconds)
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r29_3_addons.json")
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
E_ETA = np.diag(ETA)
W1 = 0.000724023879
T0 = time.time()
out = {"versions": {"python": sys.version.split()[0], "numpy": np.__version__}}

VACUA = {"toy": dict(lam=np.array([-32.0, 1.0, 0.3, 0.0]), w=W1, g=32.0, delta=0.3),
         "R23": dict(lam=np.array([-8.0, 1.0, 0.3, 0.3]), w=25.0 * W1, g=8.0, delta=0.3)}

# ---------------- (1) the gap derivatives per eigenvalue mode ----------------
sec1 = {}
for vn, v in VACUA.items():
    lam = v["lam"]
    rows = {}
    for a in range(4):
        if lam[a] == 0.0:
            rows[f"mode {a} (lambda_a = 0)"] = "a logarithmic derivative in lambda_a is undefined at lambda_a = 0; the linear derivative d|lambda_i - lambda_j| / d lambda_a is +-1 for the pairs containing a, 0 otherwise"
            continue
        r = {}
        for i in range(4):
            for j in range(i + 1, 4):
                gap = lam[i] - lam[j]
                if abs(gap) < 1e-15:
                    r[f"gap {i}{j}"] = "degenerate pair: the gap is 0"
                    continue
                d = (1.0 if i == a else 0.0) - (1.0 if j == a else 0.0)
                r[f"gap {i}{j}"] = float(d * lam[a] / gap)
        rows[f"mode {a} (lambda_a = {lam[a]:g})"] = r
    sec1[vn] = rows
sec1["the author's chain (labelled theirs, not recomputed)"] = {"d ln(M_N/m_pi)/d ln g": {"time-axis only": 3.76, "all four eigenvalues": 3.29, "projector form": 0.0},
                                                               "source": "discussioncomment-19208635 and the Skyrmion paper v2 section 6.3; the alpha map c = 16 pi Delta^4 waits on the definition of Delta (unfilled in the reply of 2026-10-28)"}
out["1_gap_derivatives_per_eigenvalue_mode"] = sec1

# ---------------- (2) the energy density of the degenerate point ----------------
sec2 = {}
for vn, v in VACUA.items():
    lam, w = v["lam"], v["w"]
    C = np.array([np.sum(lam ** p) for p in range(1, 5)])
    f = lambda x: w * np.sum((4.0 * x ** np.arange(1, 5) - C) ** 2)
    # scan then refine: the function is a polynomial in x with several local minima
    xs = np.linspace(-60.0, 60.0, 120001)
    vals = np.array([f(x) for x in xs])
    i0 = int(np.argmin(vals))
    res = minimize_scalar(f, bounds=(xs[max(0, i0 - 3)], xs[min(len(xs) - 1, i0 + 3)]), method="bounded", options={"xatol": 1e-12})
    locmins = []
    for i in range(1, len(xs) - 1):
        if vals[i] < vals[i - 1] and vals[i] < vals[i + 1]:
            locmins.append((float(xs[i]), float(vals[i])))
    sec2[vn] = {"C_p": C.tolist(), "w": w, "lambda_star": float(res.x), "rho_deg_model_units": float(res.fun),
                "local_minima_(lambda, V4)": locmins[:8],
                "rho_deg_MeV_per_unit_volume_openwave_anchor": float(res.fun) * 0.0305,
                "note": "1 model energy unit = 0.0305 MeV (M5.8 P1); the volume unit is the lattice length unit cubed, whose fm value follows from the record's r0 once its lattice value is read off the M5.8 record (the json carries the model-unit number as the primary result)",
                "author_anchor": "not expressible from the records (R29-0 (f)): reported as model units only; the author's rho_c^(1/4) >~ 300 MeV is their bound, compared, not adopted"}
out["2_energy_density_of_the_degenerate_point"] = sec2

# ---------------- (3) the sign of the boosts' kinetic norm per term ----------------
sec3 = {}
for vn, v in VACUA.items():
    M = np.diag(-v["lam"] * np.array([1.0, 1.0, 1.0, 1.0]))      # M = N eta^-1: diag(g, 1, delta, ...)
    M = np.diag([v["g"], 1.0, v["delta"], abs(v["lam"][3])])
    rows = {}
    for j in range(1, 4):
        X = np.zeros((4, 4)); X[0, j] = X[j, 0] = 1.0
        h = np.eye(4)                                                 # h_cov at the vacuum
        qT = float(np.trace(h @ X @ h @ X) - np.trace(ETA @ X @ ETA @ X))
        qcov = float(np.trace(ETA @ X @ ETA @ X)); qh = float(np.trace(X @ X))
        N = M @ ETA; qs = 0.0; P = np.eye(4)
        for p in range(1, 5):
            qs += 0.5 * (p * np.trace(P @ X @ ETA)) ** 2; P = P @ N
        rows[f"boost0{j}"] = {"-4 I1": 0.0, "-c2 K_T (c2 = 1)": 0.5 * qT, "-(kappa/2) K_cov (kappa = 1)": 0.5 * qcov,
                              "-(kappa/2) K_h (kappa = 1)": 0.5 * qh, "K_spec (kappa_p = 1)": float(qs)}
    sec3[vn] = {"coefficient_of_epsdot2_in_L": rows,
                "signs": {"-4 I1": "0 (no quadratic norm)", "K_T": "+ (the boosts' only term)", "K_cov": "-", "K_h": "+", "K_spec": "0 (eigenvalue modes only)"}}
# the torsion row: the omega^2 coefficient of the record boost direction dE = E0 K_j at the scan points (uniform in space: k = 0)
src = open(os.path.join(HERE, "m5_32_r29_0_form.py")).read()
i = src.index("def torsion_qep("); j = src.find("\n\n\n", i)
code = src[i:j]
# expose A2: cut the function at the linearisation and return (A2, idx)
code = code[: code.index("    # Euler-Lagrange of")] + "    return A2, A1k, L0, idx\n"
def _sym_basis():
    out_ = []
    for a in range(4):
        X = np.zeros((4, 4)); X[a, a] = 1.0; out_.append(X)
    for jj in range(1, 4):
        X = np.zeros((4, 4)); X[0, jj] = X[jj, 0] = 1.0; out_.append(X)
    for ii, jj in [(1, 2), (1, 3), (2, 3)]:
        X = np.zeros((4, 4)); X[ii, jj] = X[jj, ii] = 1.0; out_.append(X)
    return np.array(out_)
ns = {"np": np, "ETA": ETA, "E_ETA": E_ETA, "BX": _sym_basis()}
exec(code, ns)
tors = {}
g_t, d_t = 32.0, 0.3
E0 = np.diag([np.sqrt(g_t), 1.0, np.sqrt(d_t), 0.0])
for pname, (a1, a2, a3) in {"TEGR (1/4, 1/2, -1)": (0.25, 0.5, -1.0), "TEGR healthy sign (-1/4, -1/2, 1)": (-0.25, -0.5, 1.0),
                            "a1 only (1, 0, 0)": (1.0, 0.0, 0.0), "a1 only healthy (-1/4, 0, 0)": (-0.25, 0.0, 0.0)}.items():
    A2, A1k, L0, idx = ns["torsion_qep"](g_t, d_t, a1, a2, a3, 0.0, "rank3_pseudoinverse", (0.0, 0.0, 1.0), 0.0, np.zeros((10, 10)))
    r = {}
    for jb in range(1, 4):
        K = np.zeros((4, 4)); K[0, jb] = K[jb, 0] = 1.0
        dE = E0 @ K                                                    # E -> E Lambda^T, Lambda = exp(theta K): dE = E0 K^T = E0 K
        x = np.array([dE[a, m] for (a, m) in idx])
        r[f"record boost 0{jb}"] = {"omega2_coefficient_x^T A2 x": float(x @ A2 @ x), "|dE|^2": float(x @ x),
                                   "sign": int(np.sign(x @ A2 @ x)) if abs(x @ A2 @ x) > 1e-12 else 0}
    tors[pname] = r
sec3["torsion_locked_variant_(record boost direction dE = E0 K_j, k = 0)"] = tors
out["3_sign_of_the_boosts_kinetic_norm"] = sec3
out["wall_s"] = round(time.time() - T0, 1)
json.dump(out, open(OUT, "w"), indent=1, default=float)
print(json.dumps({"rho_deg": {k: v["rho_deg_model_units"] for k, v in sec2.items()}, "lambda_star": {k: v["lambda_star"] for k, v in sec2.items()},
                  "torsion_boost_signs": {p: {k: v["sign"] for k, v in r.items()} for p, r in tors.items()}}, indent=1))
print(f"wrote {os.path.relpath(OUT, HERE)}")
