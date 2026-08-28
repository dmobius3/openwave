"""M5.32 R1 arm (c) INDEPENDENT ADVERSARIAL AUDIT of the term screen.

Built WITHOUT reading the producer's script or JSON. Every background,
channel and probe is rebuilt from the record scripts (m5_21_16_b_field.py,
m5_21_14_c_minimize.py, m5_21_3_a_4d.py); every term density is
re-implemented here by explicit einsum contractions (NOT the registry's
256x256 K matrices). The registry is imported only for the stored 3D
non-vacuum field loader and a cross-check row.

EQUATIONS FIRST
---------------
Field M(x) real symmetric 4x4, eta = diag(-1, 1, 1, 1). Jets A_mu = d_mu M
(mu = 1..3 the certified sym stencil, 1/2(fwd + bwd) with the density
averaged per branch); A_0 = omega a0.
    F_{mu nu} = A_mu eta A_nu - A_nu eta A_mu            F[mu, nu, a, b]
Contraction rule (derivative pair: eta^{mu nu}; internal pair: eta_{ab};
mixed derivative-internal pair: delta):
    I1 = 1/2 sum eta^mu eta^nu eta_a eta_b F[mu nu a b]^2
    I2 = sum F[mu nu a b] F[a b mu nu]
    I3 = sum eta^mu eta_b F[mu nu a b] F[mu a nu b]
    R[nu, a] = sum_mu F[mu nu a mu]
    I4 = sum eta^nu eta_a R[nu a]^2      I5 = sum R[nu a] R[a nu]
    I6 = (sum_nu R[nu nu])^2
    I1_frob = I1 with eta_a eta_b -> 1 (the M5.21.16 flip; control)
Every term is exactly quadratic in omega: I(omega) = A + B omega + C omega^2
with C = (I(+1) + I(-1))/2 - I(0). Energy convention (the producer's):
    L' = -4 (I1 + c T) - V4,  E = -4 (A_I1 + c A_T) ... static part
    K(c) = C_E(I1) + c C_E(T),  C_E = -4 C   (energy kinetic coefficient)
    static dressing energy of a term: dE_T = 4 h^3 sum (dens_T - dens_T,base)
"stable" = K(c) >= 0 on every channel; "reversal" = K(c) > 0 on both
boost channels. For I1, C_E equals the certified B3.kin_of (checked).

CHANNELS. The record catalog (m5_21_3_a_4d.gen_catalog) builds
    a0_probe = w(x) (G M - M G^T) / ||.||_F
which is ANTISYMMETRIC for symmetric M (m5_theory_canonical.md: "a probe,
not a motion of symmetric M"); the physical conjugation tangent is
    a0_tan = w(x) (G M + M G^T) / ||.||_F.
Both forms are screened here, plus EXTRA channels: a radial boost twist
K(x) = n_hat(x) . K_vec (spatially varying boost direction), a combined
boost+rotation generator Kz + Jz, a random smooth symmetric a0, and the
clock on the stored NON-VACUUM 3D field (M5.21.11 end state, embedded).

S3 family (m5_21_14 wide plateau dressing): b(r) = amp tanh(r/2),
Q_b = 1 + sinh(b) K + (cosh(b) - 1) K2 with K the radial boost, the field
M_d = Q_b M_h Q_b^T on the analytic hedgehog M_h (C14.m4h_batch, g = 32).
Large-amplitude structure: Q_b ~ e^{|b|}, A ~ e^{2b}, F ~ e^{4b},
dens ~ e^{8b}: the dressing energy is a Laurent polynomial in e^{amp};
its leading coefficient is read from a fit E(amp) = sum_k c_k e^{2 k amp}.

S4 bumps (n = 48, L = 72, rho): static boost bump M = Q(eps f n_z) d Q^T
on the vacuum d = diag(g, 1, delta, 0) (s = -1) with f = exp(-|x-x0|^2 /
(2 rho^2)); omega bump M = R(eps f Jz) d R^T with a0 = Jz M - M Jz,
omega = 1. Overlap(d) = E(bump1 + bump2) - E(bump1) - E(bump2).

S5: A'_mu = (Lambda^{-T})_mu^nu Lambda A_nu Lambda^T, M' = Lambda M
Lambda^T on random jets; conjugation-only control drops the Lambda^{-T}.

Out: ../data/m5_32_r1_audit_screen.json
Run: python3 m5_32_r1_audit_screen.py [--fast]
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.linalg import expm
from scipy.ndimage import gaussian_filter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(DATA, "m5_32_r1_audit_screen.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


B3 = _load("b3", "m5_21_3_a_4d.py")
C14 = _load("c14", "m5_21_14_c_minimize.py")
REG = _load("reg", "m5_32_lagrangian.py")

ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
ETAD = np.diag(ETA)
TERMS = ["I1", "I2", "I3", "I4", "I5", "I6"]
ALL = TERMS + ["I1_frob"]
T0 = time.time()
FAST = "--fast" in sys.argv


def log(**kw):
    kw["t_s"] = round(time.time() - T0, 1)
    print(json.dumps(kw), flush=True)


# ================= own densities =================
def F_of(A):
    """A (4, N, 4, 4) -> F (N, 4, 4, 4, 4) = F[mu nu a b]."""
    AE = A @ ETA
    P = np.einsum("mnab,knbc->nmkac", AE, A, optimize=True)
    return P - P.swapaxes(1, 2)


def densities(A):
    """dict term -> density per cell, own einsum contractions."""
    F = F_of(A)
    out = {}
    out["I1"] = 0.5 * np.einsum("m,n,a,b,xmnab,xmnab->x", ETAD, ETAD, ETAD,
                                ETAD, F, F, optimize=True)
    out["I1_frob"] = 0.5 * np.einsum("m,n,xmnab,xmnab->x", ETAD, ETAD, F, F,
                                     optimize=True)
    out["I2"] = np.einsum("xmnab,xabmn->x", F, F, optimize=True)
    out["I3"] = np.einsum("m,b,xmnab,xmanb->x", ETAD, ETAD, F, F,
                          optimize=True)
    R = np.einsum("xmnam->xna", F)
    out["I4"] = np.einsum("n,a,xna,xna->x", ETAD, ETAD, R, R, optimize=True)
    out["I5"] = np.einsum("xna,xan->x", R, R, optimize=True)
    out["I6"] = np.einsum("xnn->x", R) ** 2
    return out


def lattice_sums(M, h, a0=None, chunk=None):
    """h^3-weighted sums of each term density, sym stencil (per-branch
    density averaged); a0 given -> A_0 = a0 (omega = 1)."""
    n = M.shape[0]
    tot = {k: 0.0 for k in ALL}
    chunk = chunk or n
    for br in ("fwd", "bwd"):
        Ax = [B3.d1(M, ax, h, br) for ax in range(3)]
        for i0 in range(0, n, chunk):
            sl = slice(i0, i0 + chunk)
            A = np.zeros((4, (min(i0 + chunk, n) - i0) * n * n, 4, 4))
            for ax in range(3):
                A[1 + ax] = Ax[ax][sl].reshape(-1, 4, 4)
            if a0 is not None:
                A[0] = a0[sl].reshape(-1, 4, 4)
            d = densities(A)
            for k in ALL:
                tot[k] += 0.5 * float(np.sum(d[k]))
    return {k: h ** 3 * v for k, v in tot.items()}


def omega_C(M, h, a0, chunk=None):
    """C_E = -4 C per term, C the omega^2 coefficient."""
    l0 = lattice_sums(M, h, None, chunk)
    lp = lattice_sums(M, h, a0, chunk)
    lm = lattice_sums(M, h, -a0, chunk)
    return {k: -4.0 * (0.5 * (lp[k] + lm[k]) - l0[k]) for k in ALL}


def static_E(M, h, chunk=None):
    """4 h^3 sum dens per term (the static energy read, coef -4)."""
    s = lattice_sums(M, h, None, chunk)
    return {k: 4.0 * v for k, v in s.items()}


# ================= backgrounds =================
def hedgehog_M(cfg, g):
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    P = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    P[:, 0] += 1e-9
    n = cfg["n"]
    return C14.m4h_batch(P, g).reshape(n, n, n, 4, 4)


JZ = np.zeros((4, 4)); JZ[1, 2], JZ[2, 1] = -1.0, 1.0
JX = np.zeros((4, 4)); JX[2, 3], JX[3, 2] = -1.0, 1.0
KZ = np.zeros((4, 4)); KZ[0, 3] = KZ[3, 0] = 1.0
KX = np.zeros((4, 4)); KX[0, 1] = KX[1, 0] = 1.0


def local_rot(vhat):
    W = np.zeros(vhat.shape[:-1] + (4, 4))
    n1, n2, n3 = vhat[..., 0], vhat[..., 1], vhat[..., 2]
    W[..., 1, 2], W[..., 1, 3] = -n3, n2
    W[..., 2, 1], W[..., 2, 3] = n3, -n1
    W[..., 3, 1], W[..., 3, 2] = -n2, n1
    return W


def local_boost(vhat):
    K = np.zeros(vhat.shape[:-1] + (4, 4))
    K[..., 0, 1:] = vhat
    K[..., 1:, 0] = vhat
    return K


def unit(a):
    return a / max(np.sqrt(np.sum(a * a)), 1e-300)


def channels(cfg, M, rng, extras=True):
    """(name -> a0) for probe (G M - M G^T) and tangent (G M + M G^T)."""
    w = B3.envelope(cfg)[..., None, None]
    lam, V = np.linalg.eigh(M[..., 1:4, 1:4])
    gens = {"clock_local": local_rot(V[..., :, 2]),
            "plane_1d": local_rot(V[..., :, 0]),
            "rot_z": np.broadcast_to(JZ, M.shape),
            "rot_x": np.broadcast_to(JX, M.shape),
            "boost_z": np.broadcast_to(KZ, M.shape),
            "boost_x": np.broadcast_to(KX, M.shape)}
    if extras:
        X, Y, Z = B3.coords(cfg["n"], cfg["h"])
        R = np.sqrt(X * X + Y * Y + Z * Z) + 1e-12
        nhat = np.stack([X / R, Y / R, Z / R], axis=-1)
        gens["x_boost_twist_radial"] = local_boost(nhat)
        gens["x_boost_twist_azim"] = local_boost(
            np.stack([-Y / R, X / R, np.zeros_like(Z)], axis=-1))
        gens["x_boost_plus_rot_z"] = np.broadcast_to(KZ + JZ, M.shape)
        gens["x_boost_z_rot_x"] = np.broadcast_to(KZ + JX, M.shape)
    out = {}
    for nm, G in gens.items():
        GT = G.swapaxes(-1, -2)
        out[nm + "|probe"] = unit(w * (G @ M - M @ GT))
        out[nm + "|tan"] = unit(w * (G @ M + M @ GT))
    if extras:
        n = cfg["n"]
        rnd = np.stack([[gaussian_filter(rng.normal(size=(n,) * 3), 2.0)
                         for _ in range(4)] for _ in range(4)], axis=-1)
        rnd = B3.sym4(rnd.reshape(n, n, n, 4, 4))
        out["x_random_smooth|tan"] = unit(w * rnd)
        out["x_random_smooth_antisym|probe"] = unit(
            w * (rnd @ ETA - ETA @ rnd))      # antisymmetric probe control
    return out


# ================= S1 + S2 =================
def windows(table, chan_names, terms=("I2", "I3", "I4", "I5", "I6")):
    """stable-with-reversal c windows from a C_E table."""
    out = {}
    for T in terms:
        lo, hi = -np.inf, np.inf
        binding = {}
        for ch in chan_names:
            k1, kt = table[ch]["I1"], table[ch][T]
            boost = ch.split("|")[0].startswith("boost")
            # K(c) = k1 + c kt >= 0 (> 0 for boosts)
            if abs(kt) < 1e-14:
                if k1 < 0 or (boost and k1 <= 0):
                    lo, hi = np.inf, -np.inf
                    binding[ch] = "empty (kt = 0, k1 <= 0)"
                continue
            c0 = -k1 / kt
            if kt > 0:
                if c0 > lo:
                    lo = c0; binding["lo"] = ch
            else:
                if c0 < hi:
                    hi = c0; binding["hi"] = ch
        out[T] = {"c_lo": float(lo), "c_hi": float(hi),
                  "nonempty": bool(lo <= hi), "binding": binding}
    return out


def stage_s1_s2(rng):
    cfg = B3.base_cfg(s=-1.0, g=32.0, n=32, L=48.0, delta=0.3)
    M = hedgehog_M(cfg, 32.0)
    chans = channels(cfg, M, rng)
    table = {}
    for nm, a0 in chans.items():
        ce = omega_C(M, cfg["h"], a0)
        ce["kin_of_B3"] = float(B3.kin_of(M, a0, cfg))
        table[nm] = ce
        log(S1=nm, **{k: round(v, 5) for k, v in ce.items()})
    # exact-identity probe for I2 = I6 = 0: generic pointwise jets
    ident = {}
    Ag = B3.sym4(rng.normal(size=(4, 200, 4, 4)))
    d = densities(Ag)
    ident["generic_symmetric_jets_maxabs"] = {
        k: float(np.max(np.abs(d[k]))) for k in ALL}
    # a0 = G M -/+ M G^T with constant Lorentz G on generic M, generic A_i
    Mg = B3.vac4(cfg)[None] + 0.5 * B3.sym4(rng.normal(size=(200, 4, 4)))
    for gname, G in (("Jz", JZ), ("Kz", KZ)):
        for form, sgn in (("probe", -1.0), ("tan", +1.0)):
            a0 = G @ Mg + sgn * Mg @ G.T
            A = B3.sym4(rng.normal(size=(4, 200, 4, 4)))
            A[0] = a0
            A0 = A.copy(); A0[0] = 0.0
            Am = A.copy(); Am[0] = -a0
            dp, dm, d0 = densities(A), densities(Am), densities(A0)
            ident[f"{gname}|{form}"] = {
                k: float(np.max(np.abs(0.5 * (dp[k] + dm[k]) - d0[k])))
                for k in ALL}
    # windows
    cat_probe = [c for c in table if c.endswith("|probe")
                 and not c.startswith("x_")]
    cat_tan = [c for c in table if c.endswith("|tan")
               and not c.startswith("x_")]
    all_probe = [c for c in table if c.endswith("|probe")]
    all_tan = [c for c in table if c.endswith("|tan")]
    win = {"catalog_probe": windows(table, cat_probe),
           "catalog_probe_plus_extras": windows(table, all_probe),
           "catalog_tangent": windows(table, cat_tan),
           "catalog_tangent_plus_extras": windows(table, all_tan),
           "everything": windows(table, list(table))}
    # stored non-vacuum 3D field (M5.21.11 end state), clock channels
    cfg3 = B3.base_cfg(s=1.0)
    M3 = B3.embed34(REG.load_stored3(), cfg3)
    ch3 = channels(cfg3, M3, rng, extras=False)
    table3 = {}
    for nm in ("clock_local|probe", "clock_local|tan", "rot_z|probe",
               "rot_z|tan", "boost_z|probe", "boost_z|tan",
               "boost_x|tan"):
        ce = omega_C(M3, cfg3["h"], ch3[nm])
        ce["kin_of_B3"] = float(B3.kin_of(M3, ch3[nm], cfg3))
        table3[nm] = ce
        log(S1_stored3=nm, **{k: round(v, 5) for k, v in ce.items()})
    return {"toy_point": {k: cfg[k] for k in ("s", "g", "n", "L", "h",
                                                 "delta", "stencil")},
            "C_E_table": table, "identity_probe": ident,
            "windows": win, "stored3_nonvacuum_C_E": table3}


# ================= S3 =================
def dressed_lattice(cfg, g, bfun):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    P = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    P[:, 0] += 1e-9
    M = C14.m4h_batch(P, g)
    K, K2, r = C14.kgeom(P)
    Qb = C14.qb_from(K, K2, bfun(r))
    Md = Qb @ M @ np.swapaxes(Qb, -1, -2)
    return Md.reshape(n, n, n, 4, 4), M.reshape(n, n, n, 4, 4)


def stage_s3():
    out = {"ladder": {}}
    ns = (64, 96) if FAST else (64, 96, 128)
    amps = (0.05, 0.1, 0.2, 0.4)
    for n in ns:
        cfg = B3.base_cfg(s=-1.0, g=32.0, n=n, L=24.0, delta=0.3)
        chunk = max(1, int(2.0e5 // (n * n)))
        Md0, M0 = dressed_lattice(cfg, 32.0, lambda r: np.zeros_like(r))
        E0 = static_E(M0, cfg["h"], chunk)
        row = {"h": cfg["h"], "E_base": E0, "amps": {}}
        for amp in amps:
            Md, _ = dressed_lattice(cfg, 32.0,
                                    lambda r, a=amp: a * np.tanh(r / 2.0))
            E = static_E(Md, cfg["h"], chunk)
            dE = {k: E[k] - E0[k] for k in ALL}
            dE["I1+c3(-1.9)"] = dE["I1"] - 1.9 * dE["I3"]
            dE["I1+c5(-2.5)"] = dE["I1"] - 2.5 * dE["I5"]
            dE["I1+c4(-1.07)"] = dE["I1"] - 1.07 * dE["I4"]
            row["amps"][f"{amp:g}"] = dE
            log(S3=f"n{n}", amp=amp, **{k: round(v, 2) for k, v in dE.items()
                                        if k in ("I1", "I3", "I5", "I1+c3(-1.9)",
                                                 "I1+c5(-2.5)")})
        out["ladder"][f"n{n}"] = row
    # quadrature (record instrument) large-amplitude fit of the Laurent
    # leading coefficient, per term
    grid = C14.make_grid(72, 12, 24)
    P, w = grid["P"], grid["wvol"]
    K, K2, r = C14.kgeom(P)
    hfd = 1e-4
    base_sets = {}
    for ax in range(3):
        e = np.zeros(3); e[ax] = 1.0
        for k, _ in C14.RICH_SHIFTS:
            Q = P + k * hfd * e
            base_sets[(ax, k)] = (C14.m4h_batch(Q, 32.0),) + C14.kgeom(Q)

    def quad_E(bfun):
        A = np.zeros((4, P.shape[0], 4, 4))
        for ax in range(3):
            acc = 0.0
            for k, wt in C14.RICH_SHIFTS:
                M4, Kq, K2q, rq = base_sets[(ax, k)]
                Qb = C14.qb_from(Kq, K2q, bfun(rq))
                acc = acc + wt * (Qb @ M4 @ np.swapaxes(Qb, -1, -2))
            A[1 + ax] = acc / (12.0 * hfd)
        d = densities(A)
        return {k: 4.0 * float(np.sum(w * d[k])) for k in ALL}

    Eq0 = quad_E(lambda r: np.zeros_like(r))
    amps_fit = [0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]
    quad = {}
    for amp in amps_fit:
        E = quad_E(lambda r, a=amp: a * np.tanh(r / 2.0))
        quad[f"{amp:g}"] = {k: E[k] - Eq0[k] for k in ALL}
    out["quadrature_dE"] = quad
    # leading coefficient: E(amp) ~ c8 e^{8 amp} at large amp; fit on the
    # last points the local log-slope and the coefficient
    fit = {}
    big = [a for a in amps_fit if a >= 1.5]
    for k in ALL + ["I1+c3(-1.9)", "I1+c5(-2.5)"]:
        def val(a):
            if k.startswith("I1+c3"):
                return quad[f"{a:g}"]["I1"] - 1.9 * quad[f"{a:g}"]["I3"]
            if k.startswith("I1+c5"):
                return quad[f"{a:g}"]["I1"] - 2.5 * quad[f"{a:g}"]["I5"]
            return quad[f"{a:g}"][k]
        vals = np.array([val(a) for a in big])
        slopes = np.diff(np.log(np.abs(vals))) / np.diff(big)
        fit[k] = {"sign_at_amp3": float(np.sign(vals[-1])),
                  "log_slope_pairs": slopes.tolist(),
                  "c8_estimate_at_amp3": float(vals[-1] * np.exp(-8 * 3.0)),
                  "values_large_amp": vals.tolist()}
    out["laurent_leading_fit"] = {"amps_used": big, "per_term": fit}
    return out


# ================= S4 =================
def bump_field(cfg, g, centers, kind, eps, rho, senses=None):
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    d = np.diag([g, 1.0, 0.3, 0.0])            # s = -1 vacuum
    G = {"boost": KZ, "boost_x": KX, "rot": JZ}[kind]
    ang = np.zeros((n, n, n))
    for i, c in enumerate(centers):
        s = 1.0 if senses is None else senses[i]
        ang += s * eps * np.exp(-((X - c[0]) ** 2 + (Y - c[1]) ** 2
                                  + (Z - c[2]) ** 2) / (2 * rho * rho))
    Gm = ang[..., None, None] * G
    Q = np.zeros(Gm.shape)
    if kind.startswith("boost"):
        Q[:] = (np.eye(4) + np.sinh(ang)[..., None, None] * G
                + (np.cosh(ang) - 1.0)[..., None, None] * (G @ G))
    else:
        Q[:] = (np.eye(4) + np.sin(ang)[..., None, None] * JZ
                + (1 - np.cos(ang))[..., None, None] * (JZ @ JZ))
    M = Q @ d @ Q.swapaxes(-1, -2)
    return M, X, Y, Z


def stage_s4():
    cfg = B3.base_cfg(s=-1.0, g=32.0, n=48, L=72.0, delta=0.3)
    h = cfg["h"]
    out = {}
    # static boost bump kernel order + two-bump overlap
    stat = {}
    for eps in (0.05, 0.1):
        M, *_ = bump_field(cfg, 32.0, [(0, 0, 0)], "boost", eps, 3.0)
        stat[f"eps{eps:g}"] = static_E(M, h)
    stat["ratio_E(2eps)/E(eps)"] = {k: stat["eps0.1"][k] / stat["eps0.05"][k]
                                    if stat["eps0.05"][k] != 0 else None
                                    for k in ALL}
    ov = {}
    for dsep in (10.0, 15.0):
        M2, *_ = bump_field(cfg, 32.0, [(0, 0, -dsep / 2), (0, 0, dsep / 2)],
                            "boost", 0.1, 3.0)
        M1a, *_ = bump_field(cfg, 32.0, [(0, 0, -dsep / 2)], "boost", 0.1, 3.0)
        M1b, *_ = bump_field(cfg, 32.0, [(0, 0, dsep / 2)], "boost", 0.1, 3.0)
        E2, Ea, Eb = static_E(M2, h), static_E(M1a, h), static_E(M1b, h)
        ov[f"d{dsep:g}"] = {k: E2[k] - Ea[k] - Eb[k] for k in ALL}
        ov[f"d{dsep:g}_single"] = Ea
    out["static_boost_bump"] = {"kernel": stat, "overlap": ov}
    log(S4="static", ratio=stat["ratio_E(2eps)/E(eps)"],
        overlap_d10={k: round(v, 4) for k, v in ov["d10"].items()})

    # omega bump: a0 = Jz M - M Jz, omega = 1
    def a0_of(M, X, Y, Z, senses, centers):
        a0 = JZ @ M - M @ JZ
        if senses is not None:
            # sign mask by nearest center
            dist = np.stack([np.sqrt((X - c[0]) ** 2 + (Y - c[1]) ** 2
                                     + (Z - c[2]) ** 2) for c in centers], 0)
            idx = np.argmin(dist, axis=0)
            sgn = np.array(senses)[idx]
            a0 = sgn[..., None, None] * a0
        return a0

    om = {"kernel": {}}
    for eps in (0.1, 0.2):
        M, X, Y, Z = bump_field(cfg, 32.0, [(0, 0, 0)], "rot", eps, 3.0)
        om["kernel"][f"eps{eps:g}"] = omega_C(M, h, a0_of(M, X, Y, Z, None,
                                                          None))
    om["kernel"]["ratio_K(2eps)/K(eps)"] = {
        k: om["kernel"]["eps0.2"][k] / om["kernel"]["eps0.1"][k]
        if abs(om["kernel"]["eps0.1"][k]) > 1e-12 else None for k in ALL}
    log(S4="omega_kernel", ratio=om["kernel"]["ratio_K(2eps)/K(eps)"],
        single={k: round(v, 4) for k, v in om["kernel"]["eps0.2"].items()})
    rob = {}
    for dsep in (10.0, 15.0):
        for rho in (3.0, 4.5):
            for label, senses, angs in (("like", None, (1, 1)),
                                        ("a0_flipped_on_bump2", (1, -1),
                                         (1, 1)),
                                        ("angle_flipped_on_bump2", None,
                                         (1, -1))):
                cs = [(0, 0, -dsep / 2), (0, 0, dsep / 2)]
                M2, X, Y, Z = bump_field(cfg, 32.0, cs, "rot", 0.2, rho,
                                         senses=angs)
                K2 = omega_C(M2, h, a0_of(M2, X, Y, Z, senses, cs))
                Ma, *_ = bump_field(cfg, 32.0, [cs[0]], "rot", 0.2, rho)
                Ka = omega_C(Ma, h, a0_of(Ma, X, Y, Z, None, None))
                Mb, *_ = bump_field(cfg, 32.0, [cs[1]], "rot", 0.2 * angs[1],
                                    rho)
                a0b = a0_of(Mb, X, Y, Z, None, None)
                if senses is not None:
                    a0b = -a0b
                Kb = omega_C(Mb, h, a0b)
                key = f"d{dsep:g}_rho{rho:g}_{label}"
                rob[key] = {"overlap": {k: K2[k] - Ka[k] - Kb[k] for k in ALL},
                            "single_a": Ka}
                log(S4=key, overlap={k: round(v, 4) for k, v in
                                     rob[key]["overlap"].items()})
    om["robustness"] = rob
    out["omega_bump"] = om
    return out


# ================= S5 =================
def stage_s5(rng):
    cfg = B3.base_cfg(s=-1.0, g=32.0)
    N = 128
    M = B3.vac4(cfg)[None] + 0.5 * B3.sym4(rng.normal(size=(N, 4, 4)))
    A = B3.sym4(rng.normal(size=(4, N, 4, 4)))
    d0 = densities(A)
    out = {"full_jet": {}, "conjugation_only": {}}
    for kind in ("boost", "rotation", "mixed"):
        for k in range(3):
            if kind == "mixed":
                L = REG._lorentz(rng, "boost", 0.5) @ REG._lorentz(
                    rng, "rotation", 0.5)
            else:
                L = REG._lorentz(rng, kind, 0.5)
            LiT = np.linalg.inv(L).T
            Ac = np.einsum("ab,mxbc,dc->mxad", L, A, L)
            Af = np.einsum("mn,nxab->mxab", LiT, Ac)
            for tag, Ap in (("full_jet", Af), ("conjugation_only", Ac)):
                d1 = densities(Ap)
                out[tag][f"{kind}_{k}"] = {
                    t: float(np.max(np.abs(d1[t] - d0[t]))
                             / np.max(np.abs(d0[t]))) for t in ALL}
    out["worst_full_jet_covariant"] = max(
        v[t] for v in out["full_jet"].values() for t in TERMS)
    out["min_conjugation_only_mixed_terms"] = min(
        v[t] for v in out["conjugation_only"].values()
        for t in ("I2", "I3", "I4", "I5", "I6"))
    out["conjugation_only_I1_worst"] = max(
        v["I1"] for v in out["conjugation_only"].values())
    log(S5=out["worst_full_jet_covariant"],
        conj_only_min_mixed=out["min_conjugation_only_mixed_terms"],
        conj_only_I1=out["conjugation_only_I1_worst"])
    # registry cross-check of the own densities on the same jets
    p = REG.default_params(s=-1.0, g=32.0)
    xc = {}
    for t in ALL:
        dr = REG.REGISTRY[t].density(A, M, p)
        xc[t] = float(np.max(np.abs(dr - d0[t])) / np.max(np.abs(d0[t])))
    out["own_vs_registry_density_rel"] = xc
    log(S5_registry_xcheck=xc)
    return out


def main():
    rng = np.random.default_rng(3212)
    res = {"fast": FAST}
    res["S5"] = stage_s5(rng)
    res["S1_S2"] = stage_s1_s2(rng)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    res["S4"] = stage_s4()
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    res["S3"] = stage_s3()
    res["runtime_s"] = round(time.time() - T0, 1)
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    log(done=res["runtime_s"])


if __name__ == "__main__" and "--s4b" not in sys.argv:
    main()


# ================= S4b: bumps with two independent matrix directions ====
def stage_s4b():
    """(i) radial boost bump Q(eps f n_hat . K) on the vacuum (quartic
    check); (ii) omega bump = x-boost bump carrying the Jz clock (Kx and Jz share
    index 1; a z-boost commutes with Jz and gives F_0i = 0 exactly)
    a0 = Jz M - M Jz (F_0i ~ eps: the quadratic kernel and the omega
    overlap of I2, I5, I6); senses: a0 -> -a0 on bump 2, boost sign
    flipped on bump 2."""
    cfg = B3.base_cfg(s=-1.0, g=32.0, n=48, L=72.0, delta=0.3)
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    d = np.diag([32.0, 1.0, 0.3, 0.0])

    def radial_bump(centers, eps, rho):
        M = np.broadcast_to(d, (n, n, n, 4, 4)).copy()
        for c in centers:
            dx, dy, dz = X - c[0], Y - c[1], Z - c[2]
            r = np.sqrt(dx * dx + dy * dy + dz * dz) + 1e-12
            f = eps * np.exp(-r * r / (2 * rho * rho))
            nh = np.stack([dx / r, dy / r, dz / r], -1)
            K = local_boost(nh)
            K2 = np.zeros_like(K); K2[..., 0, 0] = 1.0
            K2[..., 1:, 1:] = nh[..., :, None] * nh[..., None, :]
            Q = (np.eye(4) + np.sinh(f)[..., None, None] * K
                 + (np.cosh(f) - 1.0)[..., None, None] * K2)
            M = Q @ M @ Q.swapaxes(-1, -2)
        return M

    out = {"radial_static": {}}
    for eps in (0.05, 0.1):
        out["radial_static"][f"eps{eps:g}"] = static_E(
            radial_bump([(0, 0, 0)], eps, 3.0), h)
    out["radial_static"]["ratio_E(2eps)/E(eps)"] = {
        k: out["radial_static"]["eps0.1"][k] / out["radial_static"]["eps0.05"][k]
        if abs(out["radial_static"]["eps0.05"][k]) > 1e-14 else None
        for k in ALL}
    for dsep in (10.0, 15.0):
        cs = [(0, 0, -dsep / 2), (0, 0, dsep / 2)]
        E2 = static_E(radial_bump(cs, 0.1, 3.0), h)
        Ea = static_E(radial_bump([cs[0]], 0.1, 3.0), h)
        Eb = static_E(radial_bump([cs[1]], 0.1, 3.0), h)
        out["radial_static"][f"overlap_d{dsep:g}"] = {
            k: E2[k] - Ea[k] - Eb[k] for k in ALL}
        out["radial_static"][f"single_d{dsep:g}"] = Ea
    log(S4b="radial_static", ratio=out["radial_static"]["ratio_E(2eps)/E(eps)"],
        ov10={k: round(v, 4) for k, v in
              out["radial_static"]["overlap_d10"].items()})

    def a0_masked(M, senses, centers):
        a0 = JZ @ M - M @ JZ
        if senses is not None:
            dist = np.stack([np.sqrt((X - c[0]) ** 2 + (Y - c[1]) ** 2
                                     + (Z - c[2]) ** 2) for c in centers], 0)
            sgn = np.array(senses)[np.argmin(dist, axis=0)]
            a0 = sgn[..., None, None] * a0
        return a0

    om = {"kernel": {}}
    for eps in (0.1, 0.2):
        M, *_ = bump_field(cfg, 32.0, [(0, 0, 0)], "boost_x", eps, 3.0)
        om["kernel"][f"eps{eps:g}"] = omega_C(M, h, a0_masked(M, None, None))
    om["kernel"]["ratio_K(2eps)/K(eps)"] = {
        k: om["kernel"]["eps0.2"][k] / om["kernel"]["eps0.1"][k]
        if abs(om["kernel"]["eps0.1"][k]) > 1e-12 else None for k in ALL}
    log(S4b="omega_kernel_boostbump_Jzclock",
        ratio=om["kernel"]["ratio_K(2eps)/K(eps)"],
        single={k: round(v, 4) for k, v in om["kernel"]["eps0.2"].items()})
    rob = {}
    for dsep in (10.0, 15.0):
        for rho in (3.0, 4.5):
            for label, senses, sg in (("like", None, (1, 1)),
                                      ("a0_flipped_on_bump2", (1, -1), (1, 1)),
                                      ("boost_flipped_on_bump2", None, (1, -1))):
                cs = [(0, 0, -dsep / 2), (0, 0, dsep / 2)]
                M2, *_ = bump_field(cfg, 32.0, cs, "boost_x", 0.2, rho, senses=sg)
                K2 = omega_C(M2, h, a0_masked(M2, senses, cs))
                Ma, *_ = bump_field(cfg, 32.0, [cs[0]], "boost_x", 0.2, rho)
                Ka = omega_C(Ma, h, a0_masked(Ma, None, None))
                Mb, *_ = bump_field(cfg, 32.0, [cs[1]], "boost_x", 0.2 * sg[1], rho)
                a0b = a0_masked(Mb, None, None)
                if senses is not None:
                    a0b = -a0b
                Kb = omega_C(Mb, h, a0b)
                key = f"d{dsep:g}_rho{rho:g}_{label}"
                rob[key] = {"overlap": {k: K2[k] - Ka[k] - Kb[k] for k in ALL},
                            "overlap_over_single": {
                                k: (K2[k] - Ka[k] - Kb[k]) / Ka[k]
                                if abs(Ka[k]) > 1e-12 else None for k in ALL},
                            "single_a": Ka}
                log(S4b=key, overlap={k: round(v, 5) for k, v in
                                      rob[key]["overlap"].items()})
    om["robustness"] = rob
    out["omega_bump_boost_Jz"] = om
    return out


if __name__ == "__main__" and "--s4b" in sys.argv:
    r = stage_s4b()
    with open(os.path.join(DATA, "m5_32_r1_audit_screen_s4b.json"), "w") as f:
        json.dump(r, f, indent=1)
    log(done_s4b=round(time.time() - T0, 1))
