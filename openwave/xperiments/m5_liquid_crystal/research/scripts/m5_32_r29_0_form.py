"""M5.32 R29-0: the form-level checks of the D8 gap run (sympy and small arrays, minutes).

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1), jets A_mu = d_mu M, N = M eta.
The certified action and the quadratic terms of the S2 list, every symbol as in the record:
    L_cert = -4 I1 - V4,  I1 = sum_{mu<nu} eta^mu eta^nu <F_mu nu, F_mu nu>_eta,  F = [A_mu, A_nu]_eta,
    V4 = w sum_p (tr N^p - C_p)^2
    K_T   = 1/2 sum_mu eta^mu [tr(h A_mu h A_mu) - tr(eta A_mu eta A_mu)],  h = eta + 2 (eta u)(eta u)^T   (R7)
    K_cov = sum_mu eta^mu tr(eta A_mu eta A_mu)
    K_h   = sum_mu eta^mu tr(h A_mu h A_mu) = K_cov + 2 K_T
    K_spec= sum_p (kappa_p / 2) d_mu S_p d^mu S_p,  S_p = tr N^p                         (the author's V10)
    I1_h  = sum_{mu<nu} eta^mu eta^nu <F, F>_h                                            (R4's h-metric curvature)
About the uniform vacuum F is bilinear in two derivative directions, so -4 I1 and I1_h give no quadratic
norm to any direction (M5.34 R1v); the quadratic pencil about M_vac is set by the listed term and V4.

(a) SIGN CONVENTION. Healthy = positive coefficient of eps_dot^2 in L for M = M_vac + eps(t) X
(= positive kinetic energy). For a term entering L as  -(kappa/2) Q(A)  with Q = sum_mu eta^mu q(A_mu),
the coefficient is +(kappa/2) q(X). The R1v table listed q(X) itself (K_cov, K_h) and q_T(X); the author's
V3 census lists +(1/2) q(X) at kappa = 1 and +(1/2) q_T at c2 = 1: ours = -2 x the author's for K_cov and K_h, -1 x for K_T.

(b) GENERATORS on the relaxed hedgehog (the stored R23-1 c = 0 field, n 32, L 48, h 1.5, G = 8,
roots (-8, 1, 0.3, 0.3), w = 25 W1): the thirteen directions
    translations  D = d_k M (sym stencil),  orbital rotations  D = (x x grad)_k M,
    internal rotations  D = [J_k, M],  internal boosts  D = K_k M + M K_k,  (L + I) rotations,
    the clock's time translation (zero on a static field),
each restricted to the free cells and normalised, S2(D) = d^2 E_stat[M + eps D] / d eps^2 by the Richardson
pair of central differences of the full gradient; the Gram rank of the span.  The K-signature census: the
omega^2 weight of the certified action is LOCAL in D,
    K(D) = 4 h^3 sum_x sum_br wt sum_i <[D, A_i]_eta, [D, A_i]_eta>_eta,
so per cell it is a 10 x 10 quadratic form k_x on the symmetric basis; its eigenvalues count the K-positive,
K-negative and K-null directions per cell (static-stiff = K-null).  The trace mode D = sigma(x) eta has
K(D) = 0 exactly ([eta, X]_eta = 0) and S2(D) from V4 (the author's V1 prediction).

(c) VACUUM PENCILS. For M = M_vac + eps X e^{i(k.x - omega t)} and a term L = -(kappa/2) Q(A) - V4:
    (k^2 kappa Q_X + H) X = omega^2 kappa Q_X X,   H = Hess V4 (10 x 10),
solved as a generalised eigenproblem per k; norm sign = sign(kappa X^T Q X); helicity content = the
Frobenius weight of X in the |h| = 0, 1, 2 eigenspaces of the rotation about k on the spatial indices
(the toy vacuum is anisotropic, so it is a content, not a quantum number).
The torsion row runs in the author's locked variant (Skyrmion paper v2, App. A, eq. A.1): field E^a_mu,
M = E^T E, E0 = diag(sqrt g, 1, sqrt delta, 0), T^a_{lm} = d_l E^a_m - d_m E^a_l,
    L_tors = a1 T_{a lm} T^{a lm} + a2 T_{a lm} T^{m l a} + a3 T^l_{a l} T_n^{a n},
a1 closed with delta_ab and eta; a2, a3 need the internal index soldered to a spacetime one, done with the
rank-three pseudo-inverse E0^+ = diag(1/sqrt g, 1, 1/sqrt delta, 0) (OURS, pre-registered; the null leg
enters a1 only). The plane wave gives a quadratic eigenproblem in omega (the omega k cross terms of T_0i),
solved by linearisation; the Krein norm X^+ (2 omega A + omega-linear part) X / (2 omega) is the mode's norm.
The lock V_K is modelled as mu_K^2 |P_fibre dE|^2 on the three non-pion internal rotations (the so(3) of the
internal indices 1..3 acting as dE = w E0), at mu_K = 0 and mu_K large (OURS; App. A gives V_K its role, not
its form). The second soldering (ii) drops the null leg from every invariant; the norm signs of (i) and (ii)
are compared at the TEGR point (1/4, 1/2, -1).

(e) Q65. The author's dressed internal coefficient (spin-hydrodynamics v2, eq. 16)
    C = Z [s_i^2 + s_j^2 - 4 Z s_i^2 s_j^2 / (Z (s_i^2 + s_j^2) + Y)]
expanded in 1/g at s_i^2 = g, s_j^2 = lambda (a spatial eigenvalue of M_vac), against M5.21.14's
T1_kin = -8 sum_i |Mdot3 v_i|^2 (g-flat, partner-independent, non-positive).

(f) UNITS. The OpenWave anchor (M5.8 P1: H_static 16.74 <-> 0.511 MeV, r0 = 2.2132 fm) and the author's
(the Skyrme calibration on M_N and M_Delta), both as stated in the records; whether each is a number.

Run from this directory:  python3 m5_32_r29_0_form.py  ->  ../data/m5_32_r29_0_form.json
The author's V3 census is read from the local theory folder when present (gitignored; reported absent otherwise).
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
import scipy.linalg as sla
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r29_0_form.json")
AUTHOR_V3 = os.environ.get(
    "AUTHOR_V3",
    os.path.join(HERE, "..", "..", "theory", "duda_2026-10-27_ldgs_scorecard_v1_2", "package",
                 "LdGS_vacuum_sectors", "V3_mode_census.json"),
)
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
E_ETA = np.diag(ETA)
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


C23 = _load("m5_32_r23_1_cscan", "m5_32_r23_1_cscan.py")
R21, R20, R0, W1, G23 = C23.R21, C23.R20, C23.R0, C23.W1, C23.G
B3, LAG = R21.B3, R21.R20.LAG
out: dict = {"versions": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": __import__("scipy").__version__,
                          "sympy": sp.__version__}}


# ================= bases and small forms =================
def sym_basis():
    out_ = []
    for a in range(4):
        X = np.zeros((4, 4)); X[a, a] = 1.0; out_.append((f"diag{a}", X))
    for j in range(1, 4):
        X = np.zeros((4, 4)); X[0, j] = X[j, 0] = 1.0; out_.append((f"boost0{j}", X))
    for i, j in [(1, 2), (1, 3), (2, 3)]:
        X = np.zeros((4, 4)); X[i, j] = X[j, i] = 1.0; out_.append((f"rot{i}{j}", X))
    return out_


BASIS = sym_basis()
BN = [b[0] for b in BASIS]
BX = np.array([b[1] for b in BASIS])                     # (10, 4, 4)


def M_vac_of(spec):
    return np.diag(list(spec))


def h_of_vac(M):
    """h_cov at a vacuum whose timelike unit eigenvector is e0: h = eta + 2 e0 e0^T = identity."""
    w, V = np.linalg.eig(M @ ETA)
    best = None
    for k in range(4):
        v = np.real(V[:, k]); n = v @ ETA @ v
        if n < 0 and (best is None or n < best[0]):
            best = (n, v / np.sqrt(-n))
    u = best[1]; eu = ETA @ u
    return ETA + 2.0 * np.outer(eu, eu)


def q_T(X, h):
    return float(np.trace(h @ X @ h @ X) - np.trace(ETA @ X @ ETA @ X))


def q_cov(X):
    return float(np.trace(ETA @ X @ ETA @ X))


def q_h(X, h):
    return float(np.trace(h @ X @ h @ X))


def q_spec(X, M, kappas=(1.0, 1.0, 1.0, 1.0)):
    """sum_p (kappa_p/2) (p tr(N^{p-1} X eta))^2 : the spectral term's quadratic form on one jet."""
    N = M @ ETA
    tot, P = 0.0, np.eye(4)
    for p in range(1, 5):
        tot += 0.5 * kappas[p - 1] * (p * np.trace(P @ X @ ETA)) ** 2
        P = P @ N
    return float(tot)


def quad_matrix(q):
    """the 10 x 10 symmetric matrix of a quadratic form q on the symmetric basis (polarisation)."""
    Q = np.zeros((10, 10))
    for a in range(10):
        for b in range(a, 10):
            Q[a, b] = Q[b, a] = 0.5 * (q(BX[a] + BX[b]) - q(BX[a]) - q(BX[b]))
    return Q


def v4_of(M, roots, w):
    N = M @ ETA
    cp = [sum(q ** p for p in [k] for q in roots) for k in range(1, 5)]
    tot, P = 0.0, np.eye(4)
    for p in range(1, 5):
        P = P @ N
        tot += (np.trace(P) - cp[p - 1]) ** 2
    return float(w * tot)


def grad_sp(M):
    """the ten-vectors v_p[a] = dS_p/dM along X_a = p tr(N^{p-1} X_a eta), p = 1..4."""
    N = M @ ETA
    V = np.zeros((4, 10)); P = np.eye(4)
    for p in range(1, 5):
        for a in range(10):
            V[p - 1, a] = p * np.trace(P @ BX[a] @ ETA)
        P = P @ N
    return V


def hess_v4(M, roots, w):
    """EXACT at a vacuum (S_p = C_p): Hess V4 = 2 w sum_p v_p v_p^T (the second-derivative term of S_p carries the
    vanishing residual). Rank <= 4: the six O(1,3) directions are exactly flat. A plain finite difference at eps 1e-4
    reads the quartic term of the flat directions as a mass (9e-4 on the toy boosts), which is why this is exact."""
    V = grad_sp(M)
    return 2.0 * w * (V.T @ V)


def hess_v4_fd_check(M, roots, w, X, eps=1e-3):
    """Richardson pair of the 5-point second derivative along X: the control of the exact Hessian."""
    f = lambda e: v4_of(M + e * X, roots, w)
    d2 = lambda e: (-f(2 * e) + 16 * f(e) - 30 * f(0.0) + 16 * f(-e) - f(-2 * e)) / (12 * e * e)
    a, b = d2(eps), d2(eps / 2)
    return float((16 * b - a) / 15.0)


# ================= (a) the sign convention on both vacua =================
log("(a) sign convention")
VACUA = {
    "toy (32, 1, 0.3, 0), w 7.24e-4": (dict(roots=(-32.0, 1.0, 0.3, 0.0), M=np.diag([32.0, 1.0, 0.3, 0.0]), w=W1)),
    "R23 (8, 1, 0.3, 0.3), w 25 W1": (dict(roots=(-8.0, 1.0, 0.3, 0.3), M=np.diag([8.0, 1.0, 0.3, 0.3]), w=25.0 * W1)),
}
sec_a = {}
for vname, v in VACUA.items():
    M, h = v["M"], h_of_vac(v["M"])
    rows = {}
    for name, X in BASIS:
        rows[name] = {
            "q_T": q_T(X, h), "q_cov": q_cov(X), "q_h": q_h(X, h), "q_spec(kappa_p=1)": q_spec(X, M),
            "coef_epsdot2_in_L: -c2 K_T (c2=1)": 0.5 * q_T(X, h),
            "coef_epsdot2_in_L: -(kappa/2) K_cov (kappa=1)": 0.5 * q_cov(X),
            "coef_epsdot2_in_L: -(kappa/2) K_h (kappa=1)": 0.5 * q_h(X, h),
            "coef_epsdot2_in_L: +K_spec (kappa_p=1)": q_spec(X, M),
        }
    sec_a[vname] = {"h_at_vacuum_is_identity": bool(np.allclose(h, np.eye(4))), "per_direction": rows}
# the author's census, when present
auth = None
if os.path.exists(AUTHOR_V3):
    try:
        auth = json.load(open(AUTHOR_V3))
    except Exception as e:  # noqa: BLE001
        auth = {"error": repr(e)}
cmp = {"author_census_found": auth is not None}
if auth and "error" not in auth:
    # the author's lists are in the order of their sym_basis; their keys as shipped
    try:
        kin = auth["kinetic_coefficient_of_epsdot^2"]
        a_names = auth["directions"]
        ours = sec_a["toy (32, 1, 0.3, 0), w 7.24e-4"]["per_direction"]
        pairs = {}
        # the author lists the coefficient of eps_dot^2 in L: +1/2 q(X) for K_cov and K_h (kappa = 1), +1/2 q_T for K_T (c2 = 1);
        # the R1v table lists eta^00 q(X) = -q(X) for K_cov / K_h and 1/2 eta^00 q_T = -1/2 q_T for K_T: ours = -2 x (K_cov, K_h), -1 x (K_T)
        for key_auth, key_ours in (("K_cov (kappa=1)", "q_cov"), ("K_h (kappa=1)", "q_h"), ("K_T (c2=1)", "q_T")):
            vals = kin[key_auth]
            diffs = [abs(0.5 * ours[n][key_ours] - vals[i]) for i, n in enumerate(a_names) if n in ours]
            pairs[key_auth] = {"author = +1/2 q(X); R1v = eta^00 q(X) (K_cov, K_h) or 1/2 eta^00 q_T (K_T)": True,
                               "factor_R1v_to_author": -2.0 if key_ours != "q_T" else -1.0,
                               "max_abs_diff_author_vs_half_q": max(diffs), "n": len(diffs)}
        pairs["L_cert(-4 I1)"] = {"author_all_zero": bool(all(abs(v) < 1e-14 for v in kin["L_cert(-4 I1)"]))}
        cmp["per_term"] = pairs
        cmp["healthy_under_author_K_h_all_ten"] = bool(all(0.5 * ours[n]["q_h"] > 0 for n in BN))
        cmp["K_cov_boosts_negative"] = bool(all(0.5 * ours[n]["q_cov"] < 0 for n in BN if n.startswith("boost")))
    except Exception as e:  # noqa: BLE001
        cmp["error"] = repr(e)
sec_a["comparison_with_the_author_V3_census"] = cmp
sec_a["reading"] = ("one convention for the run: healthy = positive coefficient of eps_dot^2 in L; the R1v table's q(X) "
                    "is -2 x the author's census for K_cov and K_h and -1 x for K_T; K_h at kappa = 1 gives every one "
                    "of the ten directions a positive coefficient (the author's 'all ten healthy'), K_T only the boosts, "
                    "K_spec only the four eigenvalue modes (diag), -4 I1 none")
sec_a["fails_if"] = "a direction whose health differs between the R1v table and the author's census beyond the stated factor"
out["a_sign_convention"] = sec_a

# ================= (b) generators, the K census, the trace mode on the stored hedgehog =================
log("(b) generators on the stored n 32 c = 0 hedgehog")
TAG = "rad_pin_d0.3_w25_c0_n32_L48"
Mh = np.load(os.path.join(C23.OUT_NPZ, TAG + ".npz"))["M"]
n = Mh.shape[0]
cfg = R21.cfg_of(n, 48.0, G23, 0.3)
p23 = R21.params_of(G23, 0.3)
pot = ("v4", R0.roots_of(cfg, degenerate=True), W1 * 25.0)
h = cfg["h"]
free = (~B3.pin_shell(n, h)).astype(float)[..., None, None]
X, Y, Z = B3.coords(n, h)


def egrad(M):
    E, Gm, info = C23.energy_grad(M, cfg, p23, pot, 0.0)
    return E, Gm


E0, G0 = egrad(Mh)
rec_b = {"field": TAG, "E": E0, "max|G| spatial free": float(np.abs(G0[..., 1:, 1:] * free).max()),
         "max|G| M00": float(np.abs(G0[..., 0, 0]).max()), "max|G| time-row": float(np.abs(G0[..., 0, 1:]).max())}


def S2_of(D, eps=1e-3):
    """D . Hess . D by the Richardson pair of central differences of the gradient (4 gradients)."""
    def d1(e):
        _, Gp = egrad(Mh + e * D); _, Gm = egrad(Mh - e * D)
        return float(np.sum(D * (Gp - Gm)) / (2 * e))
    a, b = d1(eps), d1(eps / 2)
    return float((4 * b - a) / 3.0), float(abs(a - b) / max(abs(b), 1e-300))


def norm(D):
    return D / np.sqrt(np.sum(D * D))


def d_central(M, ax):
    return 0.5 * (B3.d1(M, ax, h, "fwd") + B3.d1(M, ax, h, "bwd"))


gens = {}
A = [d_central(Mh, ax) for ax in range(3)]
for k in range(3):
    gens[f"translation_{k}"] = A[k] * free
coords = [X, Y, Z]
Jm = []
for k in range(3):
    J = np.zeros((4, 4)); i, j = [(2, 3), (3, 1), (1, 2)][k]
    J[i, j], J[j, i] = -1.0, 1.0
    Jm.append(J)
    i0, j0 = [(1, 2), (2, 0), (0, 1)][k]                   # (x x grad)_k = x_i d_j - x_j d_i (cyclic)
    gens[f"orbital_rotation_{k}"] = (coords[i0][..., None, None] * A[j0] - coords[j0][..., None, None] * A[i0]) * free
    gens[f"internal_rotation_{k}"] = (J @ Mh - Mh @ J) * free
    Kb = np.zeros((4, 4)); Kb[0, k + 1] = Kb[k + 1, 0] = 1.0
    gens[f"internal_boost_{k}"] = (Kb @ Mh + Mh @ Kb) * free
for k in range(3):
    gens[f"(L-I)_rotation_{k}"] = gens[f"orbital_rotation_{k}"] - gens[f"internal_rotation_{k}"]   # annihilates a uniaxial hedgehog (sign checked on the seed below)
gens["clock_time_translation"] = np.zeros_like(Mh)        # static field
# the control: on the analytic rad seed the orbital and internal rotations about one axis are the SAME direction
Ms = C23.CORES.seed_core(cfg, "rad", 0.3)
As = [d_central(Ms, ax) for ax in range(3)]
seed_ctrl = {}
for k in range(3):
    i0, j0 = [(1, 2), (2, 0), (0, 1)][k]
    Jk = np.zeros((4, 4)); i, j = [(2, 3), (3, 1), (1, 2)][k]; Jk[i, j], Jk[j, i] = -1.0, 1.0
    Ls = (coords[i0][..., None, None] * As[j0] - coords[j0][..., None, None] * As[i0]) * free
    Is = (Jk @ Ms - Ms @ Jk) * free
    seed_ctrl[f"axis_{k}"] = {"cos(L, I)": float(np.sum(Ls * Is) / np.sqrt(np.sum(Ls * Ls) * np.sum(Is * Is))),
                              "|L-I|/|L|": float(np.sqrt(np.sum((Ls - Is) ** 2) / np.sum(Ls * Ls)))}
lam_h = np.linalg.eigvalsh(Mh[..., 1:, 1:]); eps_b = 0.5 * (lam_h[..., 1] - lam_h[..., 0])
lam_s = np.linalg.eigvalsh(Ms[..., 1:, 1:]); eps_s = 0.5 * (lam_s[..., 1] - lam_s[..., 0])
# scale: random smooth directions
rng = np.random.default_rng(29)
r = np.sqrt(X * X + Y * Y + Z * Z)
env = np.exp(-(r / 9.0) ** 2)
kmax = 3 * 2 * np.pi / 48.0
smooth = []
for _ in range(8):
    D = np.zeros(Mh.shape)
    for m in range(6):
        kk = rng.uniform(-kmax, kmax, 3); ph = rng.uniform(0, 2 * np.pi)
        S = B3.sym4(rng.standard_normal((4, 4)))
        D += np.cos(kk[0] * X + kk[1] * Y + kk[2] * Z + ph)[..., None, None] * S
    smooth.append(norm(D * env[..., None, None] * free))
s2_smooth = [S2_of(D) for D in smooth]
scale = float(np.median([abs(v[0]) for v in s2_smooth]))
rec_b["S2_scale_median_over_8_random_smooth_unit_directions"] = scale
rec_b["S2_random_smooth"] = s2_smooth
gen_rows = {}
for name, D in gens.items():
    nn = float(np.sqrt(np.sum(D * D)))
    if nn < 1e-14:
        gen_rows[name] = {"norm": nn, "S2": 0.0, "richardson_rel_err": 0.0, "S2_over_scale": 0.0, "zero_direction": True}
        continue
    s2, err = S2_of(D / nn)
    gen_rows[name] = {"norm": nn, "S2": s2, "richardson_rel_err": err, "S2_over_scale": s2 / scale}
rec_b["generators"] = gen_rows
# the span rank: Gram matrix of the 13 normalised directions (the zero one excluded)
names13 = [k for k in gens if k != "clock_time_translation" and not k.startswith("(L-I)")]
Gv = np.array([norm(gens[k]).ravel() for k in names13])
sv = np.linalg.svd(Gv @ Gv.T, compute_uv=False)
rec_b["span"] = {"directions": names13, "gram_singular_values": sv.tolist(),
                 "rank_at_1e-6": int(np.sum(sv > 1e-6 * sv[0])), "rank_at_1e-3": int(np.sum(sv > 1e-3 * sv[0]))}
# (L+I) coincidence: cosine between orbital and -internal rotation about the same axis
shell = (r > 6.0) & (r < 9.0)
rec_b["L_vs_I_on_the_relaxed_field"] = {f"axis_{k}": {"cos(orbital, internal)": float(np.sum(norm(gens[f'orbital_rotation_{k}']) * norm(gens[f'internal_rotation_{k}']))),
                                  "|L-I| / |L|": float(np.sqrt(np.sum(gens[f'(L-I)_rotation_{k}'] ** 2) / np.sum(gens[f'orbital_rotation_{k}'] ** 2)))}
                     for k in range(3)}
rec_b["L_vs_I_on_the_analytic_seed_(control)"] = seed_ctrl
rec_b["biaxiality"] = {"relaxed field, eps = (lam_mid - lam_small)/2 on the shell r 6 to 9": {"mean": float(eps_b[shell].mean()), "max": float(eps_b[shell].max())},
                       "seed, same read": {"mean": float(eps_s[shell].mean()), "max": float(eps_s[shell].max())},
                       "reading": "the relaxed c = 0 field is biaxial (the record's transverse half-strings), so a spatial rotation and the internal rotation about the same axis are DIFFERENT directions on it (L = I only on the uniaxial seed): the generator span is twelve on this field, not nine; both rotation families are zero modes of the action to the lattice tolerance"}


# the K weight (local) and its per-cell census
def K_of(D):
    tot = 0.0
    for br, wt in B3.branches(cfg["stencil"]):
        Ab = [B3.d1(Mh, ax, h, br) for ax in range(3)]
        for i in range(3):
            F = B3.comm_eta(D, Ab[i])
            tot += wt * np.sum(B3.inner_eta(F, F))
    return float(4.0 * h ** 3 * tot)


def K_cell_matrix():
    """k_x[a, b] per cell: K on the basis pairs (polarisation), vectorised over the grid."""
    Ab = {br: [B3.d1(Mh, ax, h, br) for ax in range(3)] for br, _ in B3.branches(cfg["stencil"])}

    def kc(Xc):
        tot = np.zeros(Mh.shape[:3])
        for br, wt in B3.branches(cfg["stencil"]):
            for i in range(3):
                F = B3.comm_eta(Xc, Ab[br][i])
                tot += wt * B3.inner_eta(F, F)
        return 4.0 * h ** 3 * tot
    kx = np.zeros(Mh.shape[:3] + (10, 10))
    diag = [kc(np.broadcast_to(BX[a], Mh.shape)) for a in range(10)]
    for a in range(10):
        kx[..., a, a] = diag[a]
        for b in range(a + 1, 10):
            kx[..., a, b] = kx[..., b, a] = 0.5 * (kc(np.broadcast_to(BX[a] + BX[b], Mh.shape)) - diag[a] - diag[b])
    return kx


log("(b) K census per cell")
kx = K_cell_matrix()
ev = np.linalg.eigvalsh(kx)                               # (n, n, n, 10)
kscale = float(np.max(np.abs(ev)))
tol = 1e-8 * kscale
freeb = free[..., 0, 0] > 0
pos = int(np.sum((ev > tol) & freeb[..., None])); neg = int(np.sum((ev < -tol) & freeb[..., None])); nul = int(np.sum((np.abs(ev) <= tol) & freeb[..., None]))
rec_b["K_census_free_cells"] = {
    "K_scale_max_abs_cell_eigenvalue": kscale, "null_tolerance": tol,
    "directions_positive": pos, "directions_negative": neg, "directions_null": nul, "total": int(np.sum(freeb)) * 10,
    "per_cell_signature_histogram": {f"{int(a)}+/{int(b)}-/{int(c)}0": int(cnt) for (a, b, c), cnt in
                                      zip(*np.unique(np.stack([np.sum(ev > tol, -1), np.sum(ev < -tol, -1), np.sum(np.abs(ev) <= tol, -1)], -1)[freeb].reshape(-1, 3), axis=0, return_counts=True))},
}
# trace mode
sig = np.cos(0.3 * X + 0.2 * Y) * env * free[..., 0, 0]        # scalar field (n, n, n)
Dtr = norm(sig[..., None, None] * ETA[None, None, None])
s2tr, errtr = S2_of(Dtr)
rec_b["trace_mode_sigma_eta"] = {"K": K_of(Dtr), "K_over_K(random_smooth)": K_of(Dtr) / max(1e-300, abs(K_of(smooth[0]))),
                                 "S2": s2tr, "richardson_rel_err": errtr, "S2_over_scale": s2tr / scale,
                                 "reading": "K = 0 to round-off (F blind to eta shifts, the author's V1); S2 from V4 with its sign as measured: a K-null, non-symmetry direction, STATIC-STIFF by the pre-registered words"}
rec_b["K_random_smooth"] = [K_of(D) for D in smooth]
rec_b["fails_if"] = ("a generator direction with S2 above tolerance is reported, not projected out; a span rank other than "
                     "nine on the static pinned field is reported with the coincidence; a trace mode with K above 1e-10 of the random scale")
out["b_generators_and_K_census"] = rec_b

# ================= (c) the vacuum pencils per term =================
log("(c) vacuum pencils")


def helicity_projectors(khat):
    """orthonormal eigenbasis of the rotation about khat on the symmetric 4x4 (spatial indices), grouped by |h|."""
    kx_, ky_, kz_ = khat
    J = np.zeros((4, 4)); J[2, 3], J[3, 2] = -kx_, kx_; J[3, 1], J[1, 3] = -ky_, ky_; J[1, 2], J[2, 1] = -kz_, kz_
    L = np.zeros((10, 10))
    for b in range(10):
        dX = J @ BX[b] - BX[b] @ J                        # the generator action on X
        for a in range(10):
            L[a, b] = np.sum(BX[a] * dX) / np.sum(BX[a] * BX[a])
    sw = np.sqrt(np.array([np.sum(b * b) for b in BX]))    # orthonormal coordinates y = sw x
    Lo = np.diag(sw) @ L @ np.diag(1.0 / sw)
    w, V = np.linalg.eig(Lo)
    groups = {0: [], 1: [], 2: []}
    for i in range(10):
        hval = int(round(abs(w[i].imag)))
        groups[hval].append(V[:, i])
    P = {}
    for hval, vs in groups.items():
        if vs:
            Vm = np.array(vs).T
            Qm, _ = np.linalg.qr(Vm)
            P[hval] = np.diag(1.0 / sw) @ (Qm @ Qm.conj().T) @ np.diag(sw)
        else:
            P[hval] = np.zeros((10, 10))
    return P


def norm_vec10(x):
    # Frobenius weights of the basis (off-diagonal basis matrices have norm^2 = 2)
    wts = np.array([np.sum(b * b) for b in BX])
    return float(np.real(np.sum(wts * np.abs(x) ** 2)))


def content(x, P):
    tot = norm_vec10(x)
    return {f"|h|={hv}": norm_vec10(P[hv] @ x) / max(tot, 1e-300) for hv in (0, 1, 2)}


def pencil_rows(Q, H, kappa, kvals, khat, label):
    """(k^2 kappa Q + H) X = omega^2 kappa Q X per k; returns the finite modes with norm sign and helicity content."""
    P = helicity_projectors(khat)
    rows = []
    for kk in kvals:
        Amat = kk * kk * kappa * Q + H
        Bmat = kappa * Q
        w, V = sla.eig(Amat, Bmat)
        for i in range(10):
            if not np.isfinite(w[i]) or abs(w[i]) > 1e12:
                continue
            x = V[:, i]
            nrm = float(np.real(np.conj(x) @ Bmat @ x))
            rows.append({"k": kk, "omega2": complex(w[i]).real, "omega2_imag": complex(w[i]).imag,
                         "norm_sign": int(np.sign(nrm)) if abs(nrm) > 1e-12 else 0,
                         "helicity_content": content(x, P), "dominant_basis": BN[int(np.argmax(np.abs(x)))]})
    n_static = 10 - int(np.sum(np.abs(np.linalg.eigvalsh(Q)) > 1e-12))
    return {"label": label, "static_directions_(Q-null)": n_static, "modes": rows}


sec_c = {}
kvals = [0.0, 0.25, 0.5, 1.0]
for vname, v in VACUA.items():
    M, roots, w, hm = v["M"], v["roots"], v["w"], h_of_vac(v["M"])
    H = hess_v4(M, roots, w)
    Hev = np.linalg.eigvalsh(H)
    terms = {
        "K_T (c2 = 1)": quad_matrix(lambda X: q_T(X, hm)),
        "K_cov (kappa = 1)": quad_matrix(q_cov),
        "K_h (kappa = 1)": quad_matrix(lambda X: q_h(X, hm)),
        "K_spec (kappa_p = 1)": 2.0 * quad_matrix(lambda X: q_spec(X, M)),     # enters L with +, so kappa Q/2 = q_spec
    }
    fd_boost = hess_v4_fd_check(M, roots, w, BX[4]); fd_diag1 = hess_v4_fd_check(M, roots, w, BX[1])
    rows = {"HessV4_eigenvalues": Hev.tolist(), "HessV4_flat_directions_(exact_rank_4)": int(np.sum(np.abs(Hev) < 1e-12 * Hev[-1])),
            "FD_Richardson_check": {"boost01 (exact 0)": fd_boost, "diag1 exact": float(H[1, 1]), "diag1 FD": fd_diag1}}
    for lab, Q in terms.items():
        for khat_name, khat in (("k||x1", (1.0, 0, 0)), ("k||x3", (0, 0, 1.0))):
            rows[f"{lab}, {khat_name}"] = pencil_rows(Q, H, 1.0, kvals, khat, lab)
    rows["-4 I1 and I1_h"] = "no quadratic norm about the uniform vacuum (F bilinear in two derivative directions, M5.34 R1v): every direction static-stiff under them alone; reported, not scanned"
    sec_c[vname] = rows
out["c_vacuum_pencils_M_field"] = sec_c


# ---- the torsion row in the locked variant ----
def torsion_qep(g, delta, a1, a2, a3, mu_K, soldering, khat, kk, H_M):
    """the 16-dim quadratic eigenproblem  (omega^2 A2 + omega A1 + A0) x = 0  for the plane wave, with
    L_tors + V4(E^T E) + V_K; returns the finite omega^2 with the Krein norm sign and helicity content."""
    sig = np.array([np.sqrt(g), 1.0, np.sqrt(delta), 0.0])
    E0 = np.diag(sig)
    Einv = np.diag([1.0 / s if s > 0 else 0.0 for s in sig])            # rank-three pseudo-inverse
    legs = [0, 1, 2, 3] if soldering == "rank3_pseudoinverse" else [0, 1, 2]   # (ii) drops the null leg everywhere
    # basis of dE: 16 matrices e_{a mu}
    idx = [(a, m) for a in range(4) for m in range(4)]
    nb = len(idx)
    # T_{a l m}(dE) for a plane wave with k_mu = (k0 = -omega, k): T = i (k_l dE_{a m} - k_m dE_{a l}).
    # Quadratic forms split by powers of omega: build the coefficient matrices numerically from the bilinear forms.
    kvec = np.array([0.0, khat[0] * kk, khat[1] * kk, khat[2] * kk])

    def Tof(dE, k4):
        T = np.zeros((4, 4, 4))
        for a in range(4):
            for l in range(4):
                for m in range(4):
                    T[a, l, m] = k4[l] * dE[a, m] - k4[m] * dE[a, l]
        return T

    def invariants(T):
        # a1: delta_ab on a, eta on (l, m)
        i1 = sum(E_ETA[l] * E_ETA[m] * T[a, l, m] ** 2 for a in legs for l in range(4) for m in range(4))
        # soldered T^rho_{lm} = Einv[rho, a] T^a_{lm} (Einv diagonal: rho = a on the non-null legs)
        Ts = np.zeros((4, 4, 4))
        for a in range(4):
            if Einv[a, a] != 0.0 and a in legs:
                Ts[a] = Einv[a, a] * T[a]
        # a2: T_{rho l m} T^{m l rho}: all indices spacetime, raised with eta
        i2 = sum(E_ETA[r] * E_ETA[l] * E_ETA[m] * Ts[r, l, m] * Ts[m, l, r] for r in range(4) for l in range(4) for m in range(4))
        # a3: T^l_{r l} T_n^{r n} = (sum_l eta^l T_{l r l})(sum_n eta^n T_{n r n})  with the first index lowered by eta
        v = np.array([sum(E_ETA[l] * Ts[l, r, l] for l in range(4)) for r in range(4)])
        i3 = sum(E_ETA[r] * v[r] * v[r] for r in range(4))
        return a1 * i1 + a2 * i2 + a3 * i3

    # L_tors is quadratic in (k0, k): L = sum over omega powers; extract A2 (omega^2), A1 (omega), A0 (k only) by polarisation
    def Lmat(k4):
        Lm = np.zeros((nb, nb))
        Ts = []
        for (a, m) in idx:
            dE = np.zeros((4, 4)); dE[a, m] = 1.0
            Ts.append(Tof(dE, k4))
        vals = [invariants(T) for T in Ts]
        for i in range(nb):
            Lm[i, i] = vals[i]
            for j in range(i + 1, nb):
                Lm[i, j] = Lm[j, i] = 0.5 * (invariants(Ts[i] + Ts[j]) - vals[i] - vals[j])
        return Lm
    L0 = Lmat(kvec)                                   # omega = 0 part
    kp = kvec.copy(); kp[0] = 1.0; Lp = Lmat(kp)      # k0 = +1
    km = kvec.copy(); km[0] = -1.0; Lm_ = Lmat(km)    # k0 = -1
    A2 = 0.5 * (Lp + Lm_) - L0                         # coefficient of k0^2
    A1k = 0.5 * (Lp - Lm_)                             # coefficient of k0 (k0 = -omega)
    # the sign: the plane-wave action density is  -(L_tors) ... we take L = + L_tors - V4 - V_K as the author writes the scan
    # (a_i scanned over both signs covers the overall sign). Mass: V4(E^T E) Hessian through J: dM = E0^T dE + dE^T E0.
    Jm = np.zeros((10, nb))
    for j, (a, m) in enumerate(idx):
        dE = np.zeros((4, 4)); dE[a, m] = 1.0
        dM = E0.T @ dE + dE.T @ E0
        for i in range(10):
            Jm[i, j] = np.sum(BX[i] * dM) / np.sum(BX[i] * BX[i])
    H_E = Jm.T @ H_M @ Jm
    # the lock: the three non-pion internal rotations dE = w E0, w in so(3) on internal indices 1..3
    fib = []
    for (i, j) in [(1, 2), (1, 3), (2, 3)]:
        wmat = np.zeros((4, 4)); wmat[i, j], wmat[j, i] = 1.0, -1.0
        dE = wmat @ E0
        vecd = np.array([dE[a, m] for (a, m) in idx])
        if np.linalg.norm(vecd) > 0:
            fib.append(vecd / np.linalg.norm(vecd))
    Pf = np.zeros((nb, nb))
    if fib:
        Fm = np.array(fib).T; Qf, _ = np.linalg.qr(Fm); Pf = Qf @ Qf.T
    # Euler-Lagrange of  L = x^T (omega^2 A2 - omega A1k + A0(k)) x / 1 - x^T (H_E/2 + mu_K^2 Pf) x  ... quadratic action
    # the plane-wave action density is a Hermitian form; the stationary condition gives the QEP
    #   (omega^2 A2 - omega A1k + L0 - H_E/2 - mu_K^2 Pf) x = 0     (A1k enters with k0 = -omega)
    C0 = L0 - 0.5 * H_E - mu_K * mu_K * Pf
    C1 = -A1k
    C2 = A2
    # linearisation: [[0, I], [-C0, -C1]] y = omega [[I, 0], [0, C2]] y
    Z = np.zeros((nb, nb)); I = np.eye(nb)
    Abig = np.block([[Z, I], [-C0, -C1]]); Bbig = np.block([[I, Z], [Z, C2]])
    # The pencil is SINGULAR wherever the quadratic action has a gauge symmetry (Maxwell at the a1-only point, the local
    # Lorentz and diffeomorphism directions at the TEGR point): a direction annihilated at every omega makes QZ return
    # indeterminate (NaN) eigenvalues and spoils the degenerate physical blocks. Treatment: solve twice with a tiny
    # fixed random symmetric perturbation of C0 (eps1, eps2) and keep only the eigenvalues that are stable under it
    # (the regular part moves O(eps); the singular part's eigenvalues are arbitrary and move O(sqrt eps)).
    # (a perturbation of one size at two magnitudes does NOT discriminate: on the singular subspace the equation
    # reduces to eps R x = 0, whose roots are independent of eps; two DIFFERENT random R at one eps do)
    rng_ = np.random.default_rng(2026_10_28)
    sols = []
    for _ in range(2):
        Rreg = rng_.standard_normal((nb, nb)); Rreg = Rreg @ Rreg.T / nb
        A_ = np.block([[Z, I], [-(C0 + 1e-8 * Rreg), -C1]])
        w_, V_ = sla.eig(A_, Bbig)
        sols.append((w_, V_))
    w1, V1 = sols[0]; w2, V2 = sols[1]
    keep = []
    for i in range(len(w2)):
        if not np.isfinite(w2[i]) or abs(w2[i]) > 1e8:
            continue
        j = int(np.argmin(np.abs(w1 - w2[i])))
        if np.isfinite(w1[j]) and abs(w1[j] - w2[i]) < 1e-4 * (1.0 + abs(w2[i])):
            keep.append(i)
    w = w2[keep]; V = V2[:, keep]
    n_dropped_unstable = int(np.sum(np.isfinite(w2) & (np.abs(w2) <= 1e8))) - len(keep)
    P = helicity_projectors(khat)
    rows = []
    for i in range(len(w)):
        om = w[i]
        if not np.isfinite(om) or abs(om) > 1e8:
            continue
        x = V[:nb, i]
        if np.linalg.norm(x) < 1e-12:
            continue
        x = x / np.linalg.norm(x)
        if np.real(om) < -1e-12 and abs(np.imag(om)) < 1e-9:
            continue                                      # one of each +-omega pair (the positive root); complex ones kept
        krein = np.real(np.conj(x) @ (2 * om * C2 + C1) @ x)
        sig_ = krein * (np.sign(np.real(om)) if abs(np.real(om)) > 1e-12 else 1.0)
        rows.append({"omega": [float(np.real(om)), float(np.imag(om))], "omega2": float(np.real(om * om)),
                     "krein_norm": float(sig_), "norm_sign": int(np.sign(sig_)) if abs(sig_) > 1e-10 else 0,
                     "fibre_weight": float(np.real(np.conj(x) @ Pf @ x)),
                     "helicity_content": hel16(x, khat),
                     "dominant_component": f"dE[{idx[int(np.argmax(np.abs(x)))][0]},{idx[int(np.argmax(np.abs(x)))][1]}]"})
    rows.sort(key=lambda r: r["omega2"])
    return {"n_finite_modes": len(rows), "n_unstable_under_regularisation_(gauge_or_singular)": n_dropped_unstable, "modes": rows}


def hel16(x, khat):
    """helicity content of dE (16 components) under the diagonal rotation about khat on both indices."""
    kx_, ky_, kz_ = khat
    J = np.zeros((4, 4)); J[2, 3], J[3, 2] = -kx_, kx_; J[3, 1], J[1, 3] = -ky_, ky_; J[1, 2], J[2, 1] = -kz_, kz_
    idx_ = [(a, m) for a in range(4) for m in range(4)]
    L = np.zeros((16, 16))
    for jb, (a, m) in enumerate(idx_):
        dE = np.zeros((4, 4)); dE[a, m] = 1.0
        d = J @ dE - dE @ J
        for ia, (a2, m2) in enumerate(idx_):
            L[ia, jb] = d[a2, m2]
    w, V = np.linalg.eig(L)
    out_ = {}
    for hv in (0, 1, 2):
        cols = [V[:, i] for i in range(16) if int(round(abs(w[i].imag))) == hv]
        if cols:
            Q, _ = np.linalg.qr(np.array(cols).T)
            out_[f"|h|={hv}"] = float(np.real(np.conj(x) @ (Q @ Q.conj().T) @ x) / max(1e-300, np.real(np.conj(x) @ x)))
        else:
            out_[f"|h|={hv}"] = 0.0
    return out_


tors = {}
ctrl = {}
for kk in (0.5, 1.0):
    ctrl[f"Maxwell control k {kk:g}"] = torsion_qep(1.0, 1.0, -0.25, 0.0, 0.0, 0.0, "rank3_pseudoinverse", (0.0, 0.0, 1.0), kk, np.zeros((10, 10)))
    m = ctrl[f"Maxwell control k {kk:g}"]["modes"]
    real_ = [q for q in m if abs(q["omega"][1]) < 1e-9]
    ctrl[f"Maxwell control k {kk:g}"]["summary"] = {"real_modes": len(real_), "complex_modes": len(m) - len(real_),
        "modes_at_omega2_eq_k2_positive_norm": sum(1 for q in real_ if abs(q["omega2"] - kk * kk) < 1e-4 and q["norm_sign"] > 0),
        "expected": "8 (two transverse modes per leg, four legs: a1 closes with delta_ab so the null leg carries a1 too); no complex mode; the longitudinal and temporal components are gauge and must NOT appear as stable eigenvalues", "fails_if": "a count other than 8, or any complex mode"}
g_t, d_t = 32.0, 0.3
H_M_toy = hess_v4(np.diag([g_t, 1.0, d_t, 0.0]), (-g_t, 1.0, d_t, 0.0), W1)
TEGR = (0.25, 0.5, -1.0)
scan_pts = {"TEGR (1/4, 1/2, -1)": TEGR, "a1 only (1, 0, 0)": (1.0, 0.0, 0.0), "(1/4, 1/2, 0)": (0.25, 0.5, 0.0),
            "(1/4, 0, -1)": (0.25, 0.0, -1.0), "minus TEGR": (-0.25, -0.5, 1.0)}
for pname, (a1, a2, a3) in scan_pts.items():
    for mu_K in (0.0, 30.0):
        for sold in ("rank3_pseudoinverse", "null_leg_dropped"):
            for kk in (0.5, 1.0):
                key = f"{pname} | mu_K {mu_K:g} | {sold} | k {kk:g} along x3"
                try:
                    tors[key] = torsion_qep(g_t, d_t, a1, a2, a3, mu_K, sold, (0.0, 0.0, 1.0), kk, H_M_toy)
                except Exception as e:  # noqa: BLE001
                    tors[key] = {"error": repr(e)}
# the agreement test at the TEGR point: the norm-sign multiset of (i) vs (ii)
agree = {}
for mu_K in (0.0, 30.0):
    for kk in (0.5, 1.0):
        ki = f"TEGR (1/4, 1/2, -1) | mu_K {mu_K:g} | rank3_pseudoinverse | k {kk:g} along x3"
        kj = f"TEGR (1/4, 1/2, -1) | mu_K {mu_K:g} | null_leg_dropped | k {kk:g} along x3"
        if "modes" in tors.get(ki, {}) and "modes" in tors.get(kj, {}):
            si = sorted((round(m["omega2"], 6), m["norm_sign"]) for m in tors[ki]["modes"] if m["omega"][1] == 0.0 or abs(m["omega"][1]) < 1e-9)
            sj = sorted((round(m["omega2"], 6), m["norm_sign"]) for m in tors[kj]["modes"] if m["omega"][1] == 0.0 or abs(m["omega"][1]) < 1e-9)
            agree[f"mu_K {mu_K:g}, k {kk:g}"] = {"(i) real modes (omega2, sign)": si, "(ii) real modes (omega2, sign)": sj,
                                                  "negative_norm_count_(i)": sum(1 for _, s in si if s < 0), "negative_norm_count_(ii)": sum(1 for _, s in sj if s < 0)}
out["c_torsion_locked_variant"] = {
    "setup": "E0 = diag(sqrt 32, 1, sqrt 0.3, 0); L = a1 T.T(delta_ab) + a2 T.T^(mla) + a3 (tr T)^2 - V4(E^T E) - mu_K^2 |P_fibre dE|^2; QEP by linearisation; Krein norm sign",
    "control": ctrl, "scan": tors, "agreement_i_vs_ii_at_TEGR": agree,
    "fails_if": "the (i) and (ii) solderings disagree on the sign of a propagating mode's norm at the TEGR point (then that mode is frame-dependent at the degenerate leg); the mu_K = 0 and mu_K large readings disagree on the ghost count (then the row is V_K-dependent)",
}

# ================= (e) Q65 matching =================
log("(e) Q65 sympy expansion")
Zs, Ys, g_, lam, z = sp.symbols("Z Y g lambda z", real=True)
C = Zs * (g_ + lam - 4 * Zs * g_ * lam / (Zs * (g_ + lam) + Ys))
ser = sp.series(C, g_, sp.oo, 3).removeO()
ser_Zneg = sp.series(C.subs(Zs, -z / g_), g_, sp.oo, 3).removeO()
# C >= Z (s_i^2 - s_j^2)^2 / (s_i^2 + s_j^2) at Y -> 0 (the free-record limit is the minimum over Y >= 0 for Z > 0)
si2, sj2 = sp.symbols("si2 sj2", positive=True)
Cgen = Zs * (si2 + sj2 - 4 * Zs * si2 * sj2 / (Zs * (si2 + sj2) + Ys))
dC_dY = sp.simplify(sp.diff(Cgen, Ys))
C_Y0 = sp.simplify(Cgen.subs(Ys, 0))
out["e_Q65_matching"] = {
    "C(Z, Y, g, lambda)": str(C),
    "series_in_1/g_(Z, Y fixed)": str(sp.simplify(ser)),
    "series_with_Z = -z/g": str(sp.expand(ser_Zneg)),
    "dC/dY": str(dC_dY),
    "C_at_Y=0_(free_record)": str(sp.factor(C_Y0)),
    "T1_kin (M5.21.14 § 1.3)": "-8 sum_i |Mdot3 v_i|^2 : non-positive, g-flat at leading order, the same coefficient for every spatial partner i",
    "reading": ("for Z > 0 and Y >= 0 the author's coefficient is non-negative (its minimum over Y is the free-record value "
                "Z (s_i^2 - s_j^2)^2 / (s_i^2 + s_j^2) >= 0) and its leading term Z g grows with g, while T1_kin is non-positive and g-flat: "
                "with a positive internal inertia the two forms part at the SIGN (the Lorentzian time row against a Euclidean internal pair). "
                "The only identification that reproduces a g-flat, partner-independent, negative leading term is Z = -z/g (a negative, "
                "1/g-suppressed internal inertia): then C = -z + lambda z (3 z + Y) / ((z - Y) g) + O(1/g^2) for Y < z (a pole at Y = z), the first correction partner-dependent. "
                "Form level: 'a different object' unless the author's Z is allowed negative and O(1/g); the step where the forms part is the sign of Z"),
    "fails_if": "the series or the sign bound is wrong (sympy)",
}

# ================= (f) units =================
out["f_units"] = {
    "openwave_anchor": {"energy": "M5.8 P1: H_static 16.74 lattice units <-> 0.511 MeV, 1 unit = 0.0305 MeV (m5_8_2x_findings.md)",
                        "length": "r0 = 2.2132 fm from the same record", "expressible": True},
    "author_anchor": {"statement": "the Skyrme-sector calibration on M_N and M_Delta (Skyrmion paper v2: f_pi, e from the nucleon and Delta masses; m_N = 938.9 MeV named at its magnetic-moment check)",
                      "expressible_for_the_4x4_record_field": False,
                      "reason": "no record on this side maps a model unit of the 4 x 4 record energy to MeV through M_N and M_Delta; the author's calibration lives in the internal Skyrme sector (f_pi, e), not in the record action's V4 scale"},
    "reading": "the R29-3 energy density reports model units and the OpenWave anchor; the author's anchor is reported as not expressible from the records (the pre-registered fallback of R29-0 (f))",
}

out["wall_s"] = round(time.time() - T0, 1)
json.dump(out, open(OUT, "w"), indent=1, default=float)
log(f"wrote {os.path.relpath(OUT, HERE)}")
