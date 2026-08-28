"""M5.32 R1.a (symbolic arm): the extended catalog on three backgrounds,
the external no-go under both Coulomb gates, the Q54 column.

Consumes m5_32_lagrangian.py (I1..I6, controls) and m5_32_terms_ext.py
(E1..Er, I1_h, J1, J2, Pgrad); never modifies either.

EQUATIONS FIRST
---------------
omega decomposition (every registered term, Lagrangian read):
    I(omega) = A + B omega + C omega^2,   H_I = omega dI/domega - I = C omega^2 - A
For an action L = sum_k c_k I_k the Hamiltonian omega^2 weight on a channel is
    H_2 = sum_k c_k C_k          (the certified point: c_I1 = -4, H_2 = -4 C_I1)
so a channel is vacuum-stable iff H_2 >= 0 there; in the term's own read
(the certified sign c = -1 per unit) the kinetic weight is kin_k = -C_k and
kin_k < 0 (C_k > 0) means the term drives the runaway on that channel.

Backgrounds:
  (alpha-sym) the affine clock witness: the notebook jets M_mu = [Gamma_mu, d]_xi
      at the vacuum point d = diag(g, 1, delta, 0) (s = -1 branch), the 18
      spatial generator components (Gamma_i, Gamma~_i) free symbols, the
      time jet A_0 = omega a0 with a0 one of the six generator channels
      (rot_1 = clock_local: rotation in the (2,3) plane about the leading
      spatial eigenvector; rot_2, rot_3; boost_1..3 in the CHAN reading
      [K, d], the notebook "real" style). A, B, C are polynomials in the
      spatial symbols; C is a quadratic form (its 18x18 matrix at the toy
      point g = 32, delta = 3/10 is the sign object).
  (alpha-lat) the same witness on the lattice: the M5.21.16 CHAN field
      (the analytic hedgehog at g = 32, s = -1, n = 32, L = 48) with the
      gen_catalog channels clock_local, boost_z, boost_x (a0 = w [G, M]) and
      the Lorentz variants boost_z_lor, boost_x_lor (a0 = w (K M + M K)).
  (beta) the boost-dressed twisting hedgehog of the M5.21.14 analytic
      family at the record b* (m5_21_14_minimize.json avec), toy point,
      on the make_grid(48, 8, 16) radius set; a0 = the dressed twist
      generator (Qb a0_base Qb^T); the undressed family (b = 0) as control.
  (gamma) vacuum perturbations: M(x) = Q(x) d Q(x)^T with Q a Lorentz
      boost of rapidity eps tanh(r/2) along the radial direction
      (boost_radial), a rotation by eps tanh(r/2) about the radial axis
      (rot_radial), a single-generator scalar profile (boost_z_scalar,
      F_ij == 0 identically since all A_i are parallel), and an abelian
      eigenvalue mix (eig_mix, F_ij == 0 since diagonal jets commute);
      time channels a0 = K_z d + d K_z (boost_t, Lorentz), J_1 d - d J_1
      (rot_t = clock), E_00 (eig_t). Orders = log-log slopes in eps.

Off-shell gate (static 3x3 sector, audited relations I1 = 2a, I2 = 4b,
I6 = 4c, 2I3 = I1 + I2/2, 4I4 = 2I1 + I6, 4I5 = I2 + I6):
    a-content = 2c1 + c3 + c4,  b-content = 4c2 + c3 + c5,  c-content = c4 + c5 + 4c6
    static-null basis (b = c = 0, modulo I1):
      N1 = 4 I3 - I2 - 2 I1,  N2 = 4 I4 - I6 - 2 I1,  N3 = 4 I5 - I2 - I6

Physical gate: c1 = -4 fixed, (c2..c6, eps coefficients) free; require
H_2^chan = sum c_k Q_k^chan  positive semidefinite (as a form in the 18
spatial symbols) on every generator channel; solved as an LMI by cutting
planes over an LP (scipy linprog).

Q54: omega* = (J - B)/(2C) at J = 0 -> -B/(2C) on (beta); Legendre H =
C w^2 - A verified symbolically by sp.diff on two terms.

Out: ../data/m5_32_r1_symbolic.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
import sympy as sp
from scipy.linalg import expm
from scipy.optimize import linprog

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


X = _load("m5_32_terms_ext", "m5_32_terms_ext.py")
L0 = X.L0
B3 = L0.B3
B16 = _load("m5_21_16_b_field", "m5_21_16_b_field.py")
C14 = B16.C14

ETA = L0.ETA
T0 = time.time()
OUT = {"arm": "R1.a symbolic", "toy_point": {"g": 32.0, "delta": 0.3, "s": -1.0}}


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def terms_catalog():
    """ordered catalog: base six, the two controls, eps basis, C2 ext."""
    d = {}
    for k in ("I1", "I2", "I3", "I4", "I5", "I6", "I1_frob", "I3_mixed_eta"):
        d[k] = L0.REGISTRY[k]
    d.update(X.REGISTRY_EXT)
    return d


TERMS = terms_catalog()
BASE6 = ["I1", "I2", "I3", "I4", "I5", "I6"]
EPS = [k for k in TERMS if k.startswith("E")]
CONTROLS = {"I1_frob", "I3_mixed_eta"}
STATIC_NULL = {"N1": {"I3": 4, "I2": -1, "I1": -2},
               "N2": {"I4": 4, "I6": -1, "I1": -2},
               "N3": {"I5": 4, "I2": -1, "I6": -1}}


def parity_of(T):
    return getattr(T, "parity", "even")


# ============ A1/A2: the catalog table ============
def stage_catalog():
    rng = np.random.default_rng(3220)
    p = L0.default_params(s=-1.0, g=32.0)
    A, M = L0._random_jets(rng, 64, p)
    out = {}
    # covariance + parity drift
    trans = [L0._lorentz(rng, k) for k in ("boost", "boost", "boost",
                                           "rotation", "rotation", "rotation")]
    Pref = np.diag([1.0, -1.0, 1.0, 1.0])
    for nm, T in TERMS.items():
        d0 = T.density(A, M, p)
        sc = max(np.max(np.abs(d0)), 1e-300)
        drift = 0.0
        for Lm in trans:
            Ap, Mp = X._transform(Lm, A, M)
            drift = max(drift, float(np.max(np.abs(T.density(Ap, Mp, p) - d0)) / sc))
        Ap, Mp = X._transform(Pref, A, M)
        d1 = T.density(Ap, Mp, p)
        flip = float(np.max(np.abs(d1 + d0)) / sc)
        keep = float(np.max(np.abs(d1 - d0)) / sc)
        out[nm] = {"definition": T.definition, "hash": T.hash,
                   "covariance_drift_SO13_SO3": drift,
                   "parity_x_reflection": ("odd (flips)" if flip <= 1e-10
                                           else "even (keeps)" if keep <= 1e-10
                                           else "neither"),
                   "control": nm in CONTROLS}
    # rank on generic jets (with random M so the M-dependent terms count)
    names = [n for n in TERMS if n not in CONTROLS]
    V = np.array([TERMS[n].density(A, M, p) for n in names]).T
    V = V / np.max(np.abs(V), axis=0)
    sv = np.linalg.svd(V, compute_uv=False)
    rank = int(np.sum(sv > 1e-9))
    # rank contribution: greedy in catalog order
    contrib = {}
    cur = []
    for n in names:
        cand = cur + [n]
        r = np.linalg.matrix_rank(V[:, [names.index(c) for c in cand]], tol=1e-9)
        contrib[n] = int(r - len(cur))
        if r == len(cand):
            cur = cand
    for n in names:
        out[n]["rank_contribution_generic"] = contrib[n]
    # static 3x3 sector: block-diagonal static fields, fit on (I1, I2, I6)
    from scipy.ndimage import gaussian_filter
    cfg = B3.base_cfg(n=12, L=18.0, s=-1.0, g=32.0)
    vals = {nm: [] for nm in TERMS}
    for trial in range(6):
        M3 = np.stack([[gaussian_filter(rng.normal(size=(12,) * 3), 1.5)
                        for _ in range(3)] for _ in range(3)], axis=-1)
        M3 = M3.reshape(12, 12, 12, 3, 3)
        M3 = 0.5 * (M3 + M3.swapaxes(-1, -2)) + np.diag([1.0, 0.3, 0.0])
        M4 = B3.embed34(M3, cfg)
        for nm, T in TERMS.items():
            vals[nm].append(float(L0.term_lagrangian(T, M4, cfg, p)))
    Bm = np.array([vals["I1"], vals["I2"], vals["I6"]]).T
    Bm = Bm @ np.diag([0.5, 0.25, 0.25])            # -> (a, b, c) columns
    sc = np.max(np.abs(Bm))
    for nm in TERMS:
        v = np.array(vals[nm])
        if np.max(np.abs(v)) <= 1e-11 * sc:
            out[nm]["static_abc_content"] = [0.0, 0.0, 0.0]
            out[nm]["static_status"] = "vanishes identically on the static sector"
            continue
        c, *_ = np.linalg.lstsq(Bm, v, rcond=None)
        res = float(np.max(np.abs(Bm @ c - v)) / max(np.max(np.abs(v)), 1e-300))
        out[nm]["static_abc_content"] = [float(round(x, 10)) for x in c]
        out[nm]["static_fit_residual"] = res
        out[nm]["static_status"] = ("= (a, b, c) content above" if res <= 1e-8
                                    else "NOT in span{a, b, c} (field-dependent)")
    OUT["catalog"] = {"terms": out, "rank_generic_noncontrol": rank,
                      "singular_values": sv.tolist(),
                      "names_ranked": names,
                      "eps_report": {k: v for k, v in X.EPS_REPORT.items()
                                     if k != "reductions"}}
    log(f"catalog: rank {rank} on {len(names)} non-control terms; "
        f"eps rank {X.EPS_REPORT['rank_eps']}")


# ============ (alpha-sym): the symbolic witness ============
G = L0.G_SYM
TT = L0.T_SYM
g_s, delta_s = L0.g_s, L0.delta_s
OM = sp.Symbol("omega", real=True)
SPATIAL = [G[m][j] for m in (1, 2, 3) for j in range(3)] + \
          [TT[m][j] for m in (1, 2, 3) for j in range(3)]
CHANNELS_SYM = {"rot_1_clock": G[0][0], "rot_2": G[0][1], "rot_3": G[0][2],
                "boost_1": TT[0][0], "boost_2": TT[0][1], "boost_3": TT[0][2]}
TIME_SYMS = list(G[0]) + list(TT[0])
TOY = {g_s: sp.Integer(32), delta_s: sp.Rational(3, 10)}


def sym_expressions():
    Mmu, d = L0.notebook_jets()
    F = L0.F_of_jets_sym(Mmu)
    F["jets"] = Mmu
    p = {"s": sp.Integer(-1), "g": g_s, "delta": delta_s, "w": 1}
    ex = {}
    for nm, T in TERMS.items():
        t0 = time.time()
        ex[nm] = sp.expand(T.sympy(F, d, p))
        log(f"  sympy {nm}: {len(ex[nm].args)} terms, {time.time() - t0:.1f}s")
    return ex


def abc_of(expr, chan_sym):
    sub = {s: 0 for s in TIME_SYMS}
    sub[chan_sym] = OM
    e = sp.expand(expr.subs(sub))
    return e.coeff(OM, 0), e.coeff(OM, 1), e.coeff(OM, 2)


def form_matrix(C, const=None):
    """18x18 matrix of the quadratic part of C in SPATIAL at the toy point;
    a degree-0 part (Pgrad on boost channels) is returned through const."""
    Ct = sp.expand(C.subs(TOY))
    n = len(SPATIAL)
    Q = np.zeros((n, n))
    P = sp.Poly(Ct, *SPATIAL) if Ct != 0 else None
    if P is None:
        return Q
    for mon, coef in P.terms():
        idx = [i for i, e in enumerate(mon) for _ in range(e)]
        if sum(mon) == 0:
            if const is not None:
                const.append(float(coef))
            continue
        if sum(mon) != 2:
            raise ValueError(f"C has a degree-{sum(mon)} part in the spatial symbols")
        i, j = idx
        if i == j:
            Q[i, i] += float(coef)
        else:
            Q[i, j] += float(coef) / 2
            Q[j, i] += float(coef) / 2
    return Q


def definiteness(Q, tol=1e-9):
    w = np.linalg.eigvalsh(Q)
    sc = max(np.max(np.abs(w)), 1e-300)
    if np.max(np.abs(w)) <= 1e-12:
        return "zero", w
    if w.min() >= -tol * sc:
        return "PSD (C >= 0)", w
    if w.max() <= tol * sc:
        return "NSD (C <= 0)", w
    return "indefinite", w


def stage_alpha_sym():
    ex = sym_expressions()
    OUT["alpha_sym"] = {"definition": (
        "notebook jets M_mu = [Gamma_mu, d]_xi at d = diag(g,1,delta,0) "
        "(s = -1), spatial symbols free, A_0 = omega a0 per channel; C is a "
        "quadratic form in the 18 spatial generator components; sign object "
        "= its 18x18 matrix at g = 32, delta = 3/10"), "terms": {}}
    forms = {}
    for nm, e in ex.items():
        row = {}
        for ch, cs in CHANNELS_SYM.items():
            A_, B_, C_ = abc_of(e, cs)
            const = []
            Q = form_matrix(C_, const)
            forms[(nm, ch)] = Q
            dfn, w = definiteness(Q)
            if const:
                dfn = f"constant {const[0]:+.6g} (order eps^0) + quadratic part {dfn}"
            Bq = sp.expand(B_.subs(TOY))
            row[ch] = {"A_zero": bool(sp.expand(A_) == 0),
                       "B_zero_symbolic": bool(sp.expand(B_) == 0),
                       "B_leading_order_delta": str(sp.series(
                           sp.expand(B_.subs(g_s, 1 / delta_s)), delta_s, 0, 1
                       ).removeO()) if sp.expand(B_) != 0 else "0",
                       "C_definiteness": dfn,
                       "C_eig_min": float(w.min()), "C_eig_max": float(w.max()),
                       "C_trace": float(np.trace(Q)),
                       "C_leading_order_delta": str(sp.expand(sp.series(
                           sp.expand(C_.subs(g_s, 1 / delta_s)), delta_s, 0, 1
                       ).removeO())) if sp.expand(C_) != 0 else "0",
                       "sign_C": ("+" if dfn.startswith("PSD") else
                                  "-" if dfn.startswith("NSD") else
                                  "0" if dfn == "zero" else "+-")}
        OUT["alpha_sym"]["terms"][nm] = row
        log(f"  alpha-sym {nm}: " + ", ".join(
            f"{ch}:{row[ch]['sign_C']}" for ch in CHANNELS_SYM))
    # the static-null subspace on this witness (symbolic, g and delta free)
    sub = {}
    for k, combo in STATIC_NULL.items():
        e = sum(c * ex[n] for n, c in combo.items())
        rowk = {}
        for ch, cs in CHANNELS_SYM.items():
            A_, B_, C_ = abc_of(e, cs)
            rowk[ch] = {"A_zero": bool(sp.expand(A_) == 0),
                        "B_zero": bool(sp.expand(B_) == 0),
                        "C_zero": bool(sp.expand(C_) == 0),
                        "C_expr_short": str(sp.factor(C_))[:400]}
        sub[k] = rowk
        sub[k]["whole_expression_zero"] = bool(sp.expand(e) == 0)
        log(f"  static-null {k}: whole=0? {sub[k]['whole_expression_zero']}; "
            + ", ".join(f"{ch}: C=0 {rowk[ch]['C_zero']}" for ch in CHANNELS_SYM))
    OUT["alpha_sym"]["static_null_subspace"] = sub
    OUT["alpha_sym"]["static_null_basis"] = STATIC_NULL
    return ex, forms


# ============ (alpha-lat): the lattice witness ============
def stage_alpha_lat():
    cfgc = B3.base_cfg(s=-1.0, g=32.0, n=32, L=48.0)
    pc = L0.default_params(s=-1.0, g=32.0)
    Mc = B16.lattice_family_M(cfgc, 32.0)
    a0s = B3.gen_catalog(cfgc, Mc)
    w = B3.envelope(cfgc)[..., None, None]
    Kz = np.zeros((4, 4)); Kz[0, 3] = Kz[3, 0] = 1.0
    Kx = np.zeros((4, 4)); Kx[0, 1] = Kx[1, 0] = 1.0
    for nm, Km in (("boost_z_lor", Kz), ("boost_x_lor", Kx)):
        a0 = w * (Km @ Mc + Mc @ Km)
        a0s[nm] = a0 / np.sqrt(np.sum(a0 * a0))
    # the SYMMETRIC clock tangent [J_local, M] (gen_catalog's clock_local is
    # w (J M - M J^T) = w (J M + M J), an ANTISYMMETRIC matrix: record fact)
    lam, V = np.linalg.eigh(Mc[..., 1:4, 1:4])
    vh = V[..., :, 2]
    Wl = np.zeros(Mc.shape)
    n1, n2, n3 = vh[..., 0], vh[..., 1], vh[..., 2]
    Wl[..., 1, 2], Wl[..., 1, 3] = -n3, n2
    Wl[..., 2, 1], Wl[..., 2, 3] = n3, -n1
    Wl[..., 3, 1], Wl[..., 3, 2] = -n2, n1
    a0 = w * (Wl @ Mc - Mc @ Wl)
    a0s["clock_local_sym"] = a0 / np.sqrt(np.sum(a0 * a0))
    asym = {k: float(np.sum((a + a.swapaxes(-1, -2)) ** 2) / np.sum(a * a))
            for k, a in a0s.items()}
    chans = ["clock_local", "clock_local_sym", "boost_z", "boost_x",
             "boost_z_lor", "boost_x_lor"]
    res = {"definition": (
        "M5.21.16 CHAN field (analytic hedgehog, g = 32, s = -1, n = 32, "
        "L = 48, sym stencil); a0 = envelope * [G, M] normalized (gen_catalog) "
        "for clock_local / boost_z / boost_x; a0 = envelope * (K M + M K) "
        "normalized for the Lorentz variants; clock_local_sym = envelope * "
        "[J_local, M] normalized (the symmetric tangent)"),
        "a0_symmetry_|a+aT|^2/|a|^2 (0 = antisymmetric, 2 = symmetric)": asym,
        "terms": {}}
    for nm, T in TERMS.items():
        row = {}
        for ch in chans:
            A_, B_, C_ = L0.omega_decompose(T, Mc, cfgc, pc, a0s[ch])
            row[ch] = {"A": float(A_), "B": float(B_), "C": float(C_),
                       "sign_C": "+" if C_ > 1e-12 * max(abs(A_), 1e-300) + 1e-14
                       else "-" if C_ < -1e-12 * max(abs(A_), 1e-300) - 1e-14
                       else "0"}
        res["terms"][nm] = row
        log(f"  alpha-lat {nm}: " + ", ".join(
            f"{ch}:C={row[ch]['C']:+.4e}" for ch in chans))
    # certified cross-check: -4 C_I1 == kin (record)
    with open(os.path.join(DATA, "m5_21_16_field.json")) as f:
        f16 = json.load(f)
    chk = {ch: {"minus4C_I1": -4 * res["terms"]["I1"][ch]["C"],
                "kin_eta_stored": f16["CHAN"]["rows"][ch]["kin_eta"]}
           for ch in ("clock_local", "boost_z", "boost_x")}
    res["certified_check"] = chk
    # static-null combos + 3x3 matrix (rows channels, cols N1..N3) of C
    mat = {}
    for k, combo in STATIC_NULL.items():
        mat[k] = {ch: {q: sum(c * res["terms"][n][ch][q] for n, c in combo.items())
                       for q in ("A", "B", "C")} for ch in chans}
    res["static_null_subspace"] = mat
    M3 = np.array([[mat[k][ch]["C"] for k in STATIC_NULL]
                   for ch in ("clock_local", "boost_z", "boost_x")])
    scale = max(abs(res["terms"]["I1"][ch]["C"]) for ch in chans)
    sv = np.linalg.svd(M3, compute_uv=False)
    res["static_null_C_matrix_rows_clock_boostz_boostx"] = M3.tolist()
    res["static_null_C_matrix_singular_values"] = sv.tolist()
    res["static_null_C_matrix_rel_scale_C_I1"] = float(np.max(np.abs(M3)) / scale)
    res["static_null_C_matrix_rank_rel1e-8"] = int(np.sum(sv > 1e-8 * scale))
    OUT["alpha_lat"] = res
    return res


# ============ numeric full-jet witness at the vacuum point ============
def stage_alpha_num(n_jets=20):
    rng = np.random.default_rng(3221)
    p = L0.default_params(s=-1.0, g=32.0)
    d = np.diag([32.0, 1.0, 0.3, 0.0])
    gens = {}
    J1 = np.zeros((4, 4)); J1[2, 3], J1[3, 2] = -1.0, 1.0
    J3 = np.zeros((4, 4)); J3[1, 2], J3[2, 1] = -1.0, 1.0
    Kz = np.zeros((4, 4)); Kz[0, 3] = Kz[3, 0] = 1.0
    Kx = np.zeros((4, 4)); Kx[0, 1] = Kx[1, 0] = 1.0
    gens["rot_1_clock_chan"] = J1 @ d - d @ J1.T      # antisymmetric (gen_catalog)
    gens["rot_1_clock_sym"] = J1 @ d - d @ J1         # the orbit tangent
    gens["rot_3_chan"] = J3 @ d - d @ J3.T
    gens["boost_3_chan"] = Kz @ d - d @ Kz.T
    gens["boost_1_chan"] = Kx @ d - d @ Kx.T
    gens["boost_3_lor"] = Kz @ d + d @ Kz.T
    gens["boost_1_lor"] = Kx @ d + d @ Kx.T
    gens["eig_t"] = np.diag([1.0, 0, 0, 0])
    gens["random_sym"] = B3.sym4(rng.normal(size=(4, 4)))
    Asp = B3.sym4(rng.normal(size=(3, n_jets, 4, 4)))     # full spatial jets
    M = np.broadcast_to(d, (n_jets, 4, 4)).copy()
    res = {"definition": (
        "vacuum point M = d, 20 random FULL symmetric spatial jets (eigenvalue "
        "directions included), A_0 = omega a0 per channel; A/B/C per jet"),
        "terms": {}}
    for nm, T in TERMS.items():
        row = {}
        for ch, a0 in gens.items():
            def lag(om):
                A = np.zeros((4, n_jets, 4, 4))
                A[1:] = Asp
                A[0] = om * a0
                return T.density(A, M, p)
            l0, lp, lm = lag(0.0), lag(1.0), lag(-1.0)
            A_, B_, C_ = l0, 0.5 * (lp - lm), 0.5 * (lp + lm) - l0
            row[ch] = {"A": A_.tolist(), "B": B_.tolist(), "C": C_.tolist()}
        res["terms"][nm] = row
    # static-null combos: max |C| over jets, relative to I1's
    sub = {}
    for k, combo in STATIC_NULL.items():
        sub[k] = {}
        for ch in gens:
            Cc = sum(c * np.array(res["terms"][n][ch]["C"]) for n, c in combo.items())
            Ac = sum(c * np.array(res["terms"][n][ch]["A"]) for n, c in combo.items())
            Bc = sum(c * np.array(res["terms"][n][ch]["B"]) for n, c in combo.items())
            sc = max(np.max(np.abs(res["terms"]["I1"][ch]["C"])), 1e-300)
            sub[k][ch] = {"max_abs_C_rel_I1": float(np.max(np.abs(Cc)) / sc),
                          "max_abs_A": float(np.max(np.abs(Ac))),
                          "max_abs_B": float(np.max(np.abs(Bc))),
                          "C_signs": "".join("+" if x > 1e-9 * sc else "-"
                                             if x < -1e-9 * sc else "0" for x in Cc)}
    res["static_null_subspace"] = sub
    OUT["alpha_num"] = res
    for k in sub:
        log(f"  alpha-num {k}: " + ", ".join(
            f"{ch}:{sub[k][ch]['max_abs_C_rel_I1']:.2e}" for ch in gens))


# ============ (beta): the M5.21.14 dressed family ============
def stage_beta():
    with open(os.path.join(DATA, "m5_21_14_minimize.json")) as f:
        rec = json.load(f)
    avec = np.array(rec["avec"])
    grid = C14.make_grid(48, 8, 16)
    ec = C14.ExactCorr(grid, C14.G_MAIN)
    p = L0.default_params(s=-1.0, g=32.0)
    P = grid["P"]
    wvol = grid["wvol"]
    res = {"definition": (
        "M5.21.14 analytic boost-dressed twisting hedgehog, toy point (g = 32, "
        "delta = 0.3, s = -1), radius set make_grid(48, 8, 16) (r in "
        "[0.15, 24]), b*(r) = b_of(avec_record, r); spatial jets by the "
        "record's 4-point Richardson; a0 = Qb a0_base Qb^T (the dressed "
        "twist generator); control column = undressed b = 0"),
        "avec_record": avec.tolist(), "terms": {}}
    for lab, bfun in (("b_star", lambda r: C14.b_of(avec, r)),
                      ("b_zero", lambda r: np.zeros_like(r))):
        Asp = ec._A(bfun)
        Qb = C14.qb_from(ec.K_c, ec.K2_c, bfun(ec.r_c))
        M4 = Qb @ C14.m4h_batch(P, C14.G_MAIN) @ np.swapaxes(Qb, -1, -2)
        a0 = Qb @ ec.a0_base @ np.swapaxes(Qb, -1, -2)
        for nm, T in TERMS.items():
            def lag(om):
                A = np.zeros((4,) + M4.shape)
                for i in range(3):
                    A[1 + i] = Asp[i]
                A[0] = om * a0
                return float(np.sum(wvol * T.density(A, M4, p)))
            l0, lp, lm = lag(0.0), lag(1.0), lag(-1.0)
            A_, B_, C_ = l0, 0.5 * (lp - lm), 0.5 * (lp + lm) - l0
            om_star = -B_ / (2 * C_) if abs(C_) > 1e-14 * max(abs(A_), 1.0) else None
            res["terms"].setdefault(nm, {})[lab] = {
                "A": A_, "B": B_, "C": C_,
                "sign_C": "+" if C_ > 0 else "-" if C_ < 0 else "0",
                "omega_star_J0": om_star}
        log(f"  beta {lab}: done")
    # certified check: I1 at b* vs the record (E_corr, kin_corr) is on a
    # different grid (72,12,24); here the 4 h^3 factor -> 4 * (A, -C)
    res["certified_read"] = {
        "4A_I1_bstar_minus_4A_I1_b0 (E_corr on this grid)":
            4 * (res["terms"]["I1"]["b_star"]["A"] - res["terms"]["I1"]["b_zero"]["A"]),
        "record_E_corr_grid_72_12_24": rec["verdicts"]["E_corr_at_bstar"],
        "-4C_I1_bstar_plus_4C_I1_b0 (kin_corr on this grid)":
            -4 * (res["terms"]["I1"]["b_star"]["C"] - res["terms"]["I1"]["b_zero"]["C"]),
        "record_kin_corr_grid_72_12_24": rec["verdicts"]["kin_corr_at_bstar"]}
    OUT["beta"] = res


# ============ (gamma): vacuum perturbations ============
def _boost_radial(P, amp):
    K, K2, r = C14.kgeom(P)
    return C14.qb_from(K, K2, amp * np.tanh(r / 2.0))


def _rot_radial(P, amp):
    r = np.linalg.norm(P, axis=1)
    n = P / r[:, None]
    N = P.shape[0]
    J = np.zeros((N, 4, 4))
    J[:, 1, 2], J[:, 2, 1] = -n[:, 2], n[:, 2]
    J[:, 1, 3], J[:, 3, 1] = n[:, 1], -n[:, 1]
    J[:, 2, 3], J[:, 3, 2] = -n[:, 0], n[:, 0]
    th = amp * np.tanh(r / 2.0)
    J2 = J @ J
    return (np.eye(4)[None] + np.sin(th)[:, None, None] * J
            + (1 - np.cos(th))[:, None, None] * J2)


def _profile(P):
    return (np.sin(0.7 * P[:, 0] + 0.3) * np.cos(0.4 * P[:, 1])
            + 0.5 * np.sin(0.5 * P[:, 2] - 0.2))


def gamma_field(kind, P, eps, d):
    if kind == "boost_radial":
        Q = _boost_radial(P, eps)
        return Q @ d[None] @ np.swapaxes(Q, -1, -2)
    if kind == "rot_radial":
        Q = _rot_radial(P, eps)
        return Q @ d[None] @ np.swapaxes(Q, -1, -2)
    if kind == "boost_z_scalar":
        Kz = np.zeros((4, 4)); Kz[0, 3] = Kz[3, 0] = 1.0
        b = eps * _profile(P)
        Q = np.array([expm(bb * Kz) for bb in b])
        return Q @ d[None] @ np.swapaxes(Q, -1, -2)
    if kind == "eig_mix":
        f1 = _profile(P)
        f2 = _profile(P[:, ::-1] * 1.3)
        M = np.broadcast_to(d, (P.shape[0], 4, 4)).copy()
        M[:, 0, 0] += eps * f1
        M[:, 1, 1] += eps * f2
        return M
    raise ValueError(kind)


def gamma_jets(kind, P, eps, d, h=1e-3):
    A = np.zeros((4, P.shape[0], 4, 4))
    for ax in range(3):
        acc = 0.0
        for k, w in C14.RICH_SHIFTS:
            Q = P.copy()
            Q[:, ax] += k * h
            acc = acc + w * gamma_field(kind, Q, eps, d)
        A[1 + ax] = acc / (12.0 * h)
    return A, gamma_field(kind, P, eps, d)


def stage_gamma():
    rng = np.random.default_rng(3222)
    p = L0.default_params(s=-1.0, g=32.0)
    d = np.diag([32.0, 1.0, 0.3, 0.0])
    N = 300
    P = rng.normal(size=(N, 3))
    P = P / np.linalg.norm(P, axis=1)[:, None] * (0.5 + 5.5 * rng.random(N))[:, None]
    Kz = np.zeros((4, 4)); Kz[0, 3] = Kz[3, 0] = 1.0
    J1 = np.zeros((4, 4)); J1[2, 3], J1[3, 2] = -1.0, 1.0
    tchan = {"boost_t_lor": lambda M: Kz @ M + M @ Kz,
             "boost_t_chan": lambda M: Kz @ M - M @ Kz,
             "rot_t_clock": lambda M: J1 @ M - M @ J1.T,
             "eig_t": lambda M: np.broadcast_to(np.diag([1.0, 0, 0, 0]), M.shape)}
    epss = [0.2, 0.1, 0.05, 0.025]
    res = {"definition": (
        "300 random points in 0.5 < r < 6; spatial channels boost_radial "
        "(Lorentz boost, rapidity eps tanh(r/2) along r-hat), rot_radial "
        "(rotation by eps tanh(r/2) about r-hat), boost_z_scalar (fixed K_z, "
        "scalar profile), eig_mix (abelian eigenvalue profiles); time "
        "channels a0(M) = K_z M + M K_z (Lorentz), K_z M - M K_z (CHAN "
        "reading), J_1 M - M J_1^T (clock), E_00; A, B, C = mean density; "
        "order = log2 ratio between successive eps halvings (last pair)"),
        "eps": epss, "terms": {}}
    for kind in ("boost_radial", "rot_radial", "boost_z_scalar", "eig_mix"):
        fields = [gamma_jets(kind, P, e, d) for e in epss]
        scale_ref = None
        for nm, T in TERMS.items():
            row = {}
            for tn, tf in tchan.items():
                As, Bs, Cs = [], [], []
                for (A, M) in fields:
                    a0 = tf(M)
                    a0 = a0 / np.sqrt(np.mean(np.sum(a0 * a0, axis=(-1, -2))))
                    def lag(om):
                        AA = A.copy()
                        AA[0] = om * a0
                        return float(np.mean(T.density(AA, M, p)))
                    l0, lp, lm = lag(0.0), lag(1.0), lag(-1.0)
                    As.append(l0); Bs.append(0.5 * (lp - lm))
                    Cs.append(0.5 * (lp + lm) - l0)

                if scale_ref is None:
                    scale_ref = max(abs(x) for x in As + Cs)   # I1 first

                def order(v):
                    v = np.abs(np.array(v))
                    if v[0] <= 1e-8 * scale_ref or v[-1] <= 1e-13 * v[0]:
                        return "0 (vanishes)"
                    return f"{np.log2(v[-2] / v[-1]):.2f}"
                row[tn] = {"A": As, "B": Bs, "C": Cs,
                           "order_A": order(As), "order_B": order(Bs),
                           "order_C": order(Cs),
                           "sign_C_smallest_eps": "+" if Cs[-1] > 1e-8 * scale_ref
                           else "-" if Cs[-1] < -1e-8 * scale_ref else "0",
                           "sign_A_smallest_eps": "+" if As[-1] > 1e-8 * scale_ref
                           else "-" if As[-1] < -1e-8 * scale_ref else "0"}
            res["terms"].setdefault(nm, {})[kind] = row
        log(f"  gamma {kind}: " + ", ".join(
            f"{nm}:A^{res['terms'][nm][kind]['boost_t_lor']['order_A']}"
            for nm in ("I1", "I2", "Pgrad")))
    OUT["gamma"] = res


# ============ B2: the physical gate (LMI) ============
def stage_physical(forms):
    free = [n for n in BASE6 if n != "I1"] + EPS
    chans = list(CHANNELS_SYM)
    Qs = {ch: {n: forms[(n, ch)] for n in ["I1"] + free} for ch in chans}
    n = len(SPATIAL)
    # cutting-plane LP: maximize t s.t. v^T (sum c_k Q_k) v >= t for cuts,
    # |c_k| <= 100; c_I1 = -4 fixed
    cuts = {ch: [] for ch in chans}
    rng = np.random.default_rng(3223)
    for ch in chans:
        for _ in range(40):
            v = rng.normal(size=n)
            cuts[ch].append(v / np.linalg.norm(v))
    hist = []
    c_best = None
    for it in range(200):
        Aub, bub = [], []
        for ch in chans:
            for v in cuts[ch]:
                base = -4.0 * v @ Qs[ch]["I1"] @ v
                row = [-(v @ Qs[ch][k] @ v) for k in free] + [1.0]
                Aub.append(row)
                bub.append(base)
        cobj = np.zeros(len(free) + 1); cobj[-1] = -1.0
        bounds = [(-100, 100)] * len(free) + [(None, None)]
        lp = linprog(cobj, A_ub=np.array(Aub), b_ub=np.array(bub),
                     bounds=bounds, method="highs")
        c = lp.x[:-1]
        t = lp.x[-1]
        worst = 0.0
        for ch in chans:
            Qtot = -4.0 * Qs[ch]["I1"] + sum(ck * Qs[ch][k] for ck, k in zip(c, free))
            w, V = np.linalg.eigh(Qtot)
            worst = min(worst, w[0])
            if w[0] < -1e-9:
                cuts[ch].append(V[:, 0])
        hist.append({"it": it, "lp_t": float(t), "min_eig_over_channels": float(worst)})
        if t < -1e-9:
            verdict = ("EMPTY: the LP over the accumulated necessary cuts has "
                       f"optimum t = {t:.4g} < 0 (no coefficient vector with "
                       "|c_k| <= 100 makes every channel form PSD)")
            c_best = c
            break
        if worst >= -1e-9:
            verdict = "NON-EMPTY: a feasible coefficient vector found"
            c_best = c
            break
    else:
        verdict = "UNDECIDED after 200 cutting-plane iterations"
    res = {"free_coefficients": free, "c_I1_fixed": -4.0,
           "channels": chans, "verdict": verdict, "iterations": len(hist),
           "history_tail": hist[-5:]}
    if c_best is not None:
        cc = dict(zip(free, [float(x) for x in c_best]))
        cc["I1"] = -4.0
        res["coefficients_at_stop"] = cc
        res["static_abc_content_at_stop"] = {
            "a": 2 * cc["I1"] + cc["I3"] + cc["I4"],
            "b": 4 * cc["I2"] + cc["I3"] + cc["I5"],
            "c": cc["I4"] + cc["I5"] + 4 * cc["I6"]}
        res["channel_min_eig_at_stop"] = {}
        for ch in chans:
            Qtot = -4.0 * Qs[ch]["I1"] + sum(cc[k] * Qs[ch][k] for k in free)
            res["channel_min_eig_at_stop"][ch] = float(np.linalg.eigvalsh(Qtot)[0])
    # --- certificate: the common null space of the FREE forms per channel,
    # and I1's form on it (a direction where no coefficient can act)
    cert = {}
    for ch in chans:
        S = np.concatenate([Qs[ch][k] for k in free], axis=0)     # stack
        u, sv, vt = np.linalg.svd(S)
        null = vt[np.sum(sv > 1e-9 * sv[0]):]                     # rows
        if null.shape[0] == 0:
            cert[ch] = {"common_null_dim": 0}
            continue
        Qn = null @ Qs[ch]["I1"] @ null.T
        w = np.linalg.eigvalsh(Qn)
        cert[ch] = {"common_null_dim": int(null.shape[0]),
                    "I1_form_on_null_eig_min": float(w.min()),
                    "I1_form_on_null_eig_max": float(w.max()),
                    "H2_certified_on_null_range": [float(-4 * w.max()),
                                                   float(-4 * w.min())],
                    "verdict": ("CERTIFICATE: I1's C is positive on a spatial "
                                "direction where every free term's C vanishes, "
                                "so H_2 = -4 C_I1 < 0 there for ANY c2..c6, eps"
                                if w.max() > 1e-9 else "no certificate here")}
        # name the direction: which symbols carry it
        i = int(np.argmax(w))
        vec = np.linalg.eigh(Qn)[1][:, i] @ null
        cert[ch]["direction_top_symbols"] = [
            (str(SPATIAL[j]), float(round(vec[j], 4)))
            for j in np.argsort(-np.abs(vec))[:6]]
    res["certificate_common_null"] = cert
    # --- LP variants: wider box; c1 free with the a-content held <= -8
    def lp_variant(c1_fixed, box, a_content_max=None, extra_chans=None):
        chs = extra_chans or chans
        cuts2 = {ch: list(cuts[ch]) for ch in chs}
        names = (["I1"] if c1_fixed is None else []) + free
        for it in range(300):
            Aub, bub = [], []
            for ch in chs:
                for v in cuts2[ch]:
                    # sum_k c_k vQ_k v + c1 vQ_I1 v >= t  <=>
                    # -sum_k c_k vQ_k v + t <= c1 vQ_I1 v
                    base = 0.0 if c1_fixed is None else c1_fixed * v @ Qs[ch]["I1"] @ v
                    row = [-(v @ Qs[ch][k] @ v) for k in names] + [1.0]
                    Aub.append(row); bub.append(base)
            if a_content_max is not None:
                # 2 c1 + c3 + c4 <= a_content_max
                row = [0.0] * (len(names) + 1)
                for k, coef in (("I1", 2.0), ("I3", 1.0), ("I4", 1.0)):
                    if k in names:
                        row[names.index(k)] = coef
                Aub.append(row); bub.append(a_content_max - (0.0 if c1_fixed is None
                                                             else 2 * c1_fixed))
            cobj = np.zeros(len(names) + 1); cobj[-1] = -1.0
            bounds = [(-box, box)] * len(names) + [(None, None)]
            lp = linprog(cobj, A_ub=np.array(Aub), b_ub=np.array(bub),
                         bounds=bounds, method="highs")
            if lp.status != 0:
                return {"status": f"LP status {lp.status}", "it": it}
            c = lp.x[:-1]; t = lp.x[-1]
            worst = 0.0
            for ch in chs:
                Qtot = (0.0 if c1_fixed is None else -c1_fixed) * Qs[ch]["I1"] * 0 \
                    + sum(ck * Qs[ch][k] for ck, k in zip(c, names))
                if c1_fixed is not None:
                    Qtot = Qtot + c1_fixed * Qs[ch]["I1"]
                w, V = np.linalg.eigh(Qtot)
                worst = min(worst, w[0])
                if w[0] < -1e-9:
                    cuts2[ch].append(V[:, 0])
            if t < -1e-9:
                return {"status": "EMPTY", "lp_t": float(t), "it": it,
                        "c_at_stop": dict(zip(names, [float(x) for x in c]))}
            if worst >= -1e-9:
                return {"status": "NON-EMPTY", "lp_t": float(t), "it": it,
                        "c": dict(zip(names, [float(x) for x in c]))}
        return {"status": "UNDECIDED", "it": it}
    # --- Farkas certificate WITHOUT box bounds: y >= 0 over cuts with
    # sum_i y_i (vQ_k v)_i = 0 for every free k, sum y = 1, and
    # sum_i y_i c1 (vQ_I1 v)_i < 0: then no coefficient vector of any size
    # makes every cut nonnegative (a positive combination of necessary
    # conditions contradicts itself).
    def farkas(chs, c1_fixed=-4.0):
        n = len(SPATIAL)
        dirs = [np.eye(n)[i] for i in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                dirs.append((np.eye(n)[i] + np.eye(n)[j]) / np.sqrt(2))
                dirs.append((np.eye(n)[i] - np.eye(n)[j]) / np.sqrt(2))
        rows, bs, labels = [], [], []
        for ch in chs:
            for v in dirs + list(cuts[ch]):
                rows.append([v @ Qs[ch][k] @ v for k in free])
                bs.append(c1_fixed * (v @ Qs[ch]["I1"] @ v))
                top = np.argsort(-np.abs(v))[:2]
                labels.append((ch, " ".join(f"{v[j]:+.2f}*{SPATIAL[j]}" for j in top
                                            if abs(v[j]) > 1e-6)))
        Rm = np.array(rows); bv = np.array(bs)
        m = len(bs)
        # minimize b.y  s.t. R^T y = 0, sum y = 1, y >= 0
        Aeq = np.vstack([Rm.T, np.ones((1, m))])
        beq = np.zeros(len(free) + 1); beq[-1] = 1.0
        lp = linprog(bv, A_eq=Aeq, b_eq=beq, bounds=[(0, None)] * m, method="highs")
        if lp.status != 0:
            return {"status": f"LP status {lp.status}"}
        y = lp.x
        sup = np.where(y > 1e-9)[0]
        return {"min_b_dot_y": float(lp.fun),
                "verdict": ("EMPTY for ALL coefficient magnitudes (Farkas "
                            "certificate: the listed cuts combine to 0 >= "
                            f"{-lp.fun:.4g} > 0)" if lp.fun < -1e-9 else
                            "no unbounded certificate (a feasible point may exist)"),
                "support": [{"channel": labels[i][0], "direction": labels[i][1],
                             "y": float(y[i]),
                             "row_free_terms": dict(zip(free, [float(r) for r in Rm[i]])),
                             "c1_part": float(bv[i])} for i in sup],
                "residual_Rt_y": float(np.max(np.abs(Rm.T @ y)))}
    res["farkas"] = {"all six channels": farkas(chans),
                     "boost_3 only": farkas(["boost_3"]),
                     "clock only": farkas(["rot_1_clock"]),
                     "boost_3 + boost_1": farkas(["boost_3", "boost_1"])}
    for k, v in res["farkas"].items():
        log(f"farkas {k}: {v.get('verdict')}  support {len(v.get('support', []))}")
    res["variants"] = {
        "c1=-4, box 1e4": lp_variant(-4.0, 1e4),
        "c1 free, box 1e4, a-content <= -8": lp_variant(None, 1e4, -8.0),
        "c1 free, box 1e4, no static constraint (homogeneous cone)":
            lp_variant(None, 1e4),
        "c1=-4, box 1e4, clock + boost_3 only": lp_variant(
            -4.0, 1e4, extra_chans=["rot_1_clock", "boost_3"]),
        "c1=-4, box 1e4, clock only": lp_variant(-4.0, 1e4, extra_chans=["rot_1_clock"]),
        "c1=-4, box 1e4, boost_3 only": lp_variant(-4.0, 1e4, extra_chans=["boost_3"])}
    log(f"physical gate variants: " + "; ".join(
        f"{k}: {v['status']}" for k, v in res["variants"].items()))
    # the per-channel structure: which monomials carry each term's C
    diag_only = {}
    for ch in chans:
        offd = max(float(np.max(np.abs(Qs[ch][k] - np.diag(np.diag(Qs[ch][k])))))
                   / max(np.max(np.abs(Qs[ch][k])), 1e-300) for k in ["I1"] + free
                   if np.max(np.abs(Qs[ch][k])) > 0)
        diag_only[ch] = offd
    res["max_offdiagonal_fraction_per_channel"] = diag_only
    # scalar version: trace weights per channel (the number-level table)
    res["trace_weights"] = {ch: {k: float(np.trace(Qs[ch][k])) for k in ["I1"] + free}
                            for ch in chans}
    # the boost channels alone (the reversal question) and the clock alone
    OUT["physical_gate"] = res
    log(f"physical gate: {verdict}")


# ============ B3: the Legendre check ============
def stage_legendre(ex):
    res = {}
    for nm in ("I2", "I4", "E1", "Pgrad"):
        e = ex[nm]
        sub = {s: 0 for s in TIME_SYMS}
        sub[G[0][0]] = OM
        Lg = sp.expand(e.subs(sub))
        H = sp.expand(OM * sp.diff(Lg, OM) - Lg)
        A_, B_, C_ = Lg.coeff(OM, 0), Lg.coeff(OM, 1), Lg.coeff(OM, 2)
        res[nm] = {"L_degree_in_omega": int(sp.degree(Lg, OM)) if Lg != 0 else 0,
                   "H_equals_C_w2_minus_A": bool(sp.expand(H - (C_ * OM ** 2 - A_)) == 0),
                   "B_nonzero": bool(sp.expand(B_) != 0)}
    OUT["legendre"] = res
    log(f"legendre: {res}")


def main():
    stage_catalog()
    ex, forms = stage_alpha_sym()
    stage_alpha_num()
    stage_alpha_lat()
    stage_beta()
    stage_gamma()
    stage_physical(forms)
    stage_legendre(ex)
    OUT["runtime_s"] = round(time.time() - T0, 1)
    with open(os.path.join(DATA, "m5_32_r1_symbolic.json"), "w") as f:
        json.dump(OUT, f, indent=1)
    log("written m5_32_r1_symbolic.json")


if __name__ == "__main__":
    main()
