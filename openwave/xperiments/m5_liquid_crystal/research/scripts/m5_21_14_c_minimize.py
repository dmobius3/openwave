"""M5.21.14 S2/S3: variational b(r) on the analytic family: the exact
functional at finite g, the leading-order boundedness diagnosis, the
pre-registered sign verdicts, and the certified-lattice cross-check.

The derived leading-order (beta = g*b) corrections (S1, verified):
  T1_static[beta] = 4 sum_{i<j} [ 2<[G_i,G_j], v_j v_i^T - v_i v_j^T>
                    + |v_j v_i^T - v_i v_j^T|^2 - 2|G_i v_j - G_j v_i|^2 ]
  T1_kin[beta]    = -8 sum_i |Mdot3 v_i|^2          (omega^2 multiplier)
with v_i = d_i(beta(r) nhat).

MEASURED STRUCTURE (first pass, recorded as the DIAG stage): the free
minimization of T1_static ALONE is UNBOUNDED BELOW: for pure-radial
v_i (grid-scale oscillation: large beta', small beta) the positive
quartic |W|^2 vanishes identically while the negative
-2|G_i v_j - G_j v_i|^2 channel grows as beta'^2. The leading term is
ill-posed as a standalone objective; its stabilizers are the
g-suppressed higher orders (or a constraint). Stages:
  DIAG  reproduce + record the runaway (leading order, grid family)
  BND   probe the EXACT functional at g = 32 on sawtooth families
        (does finite g restore a floor?)
  MIN   minimize the EXACT dressing correction at g = 32 within a
        smooth 10-dim radial family (plateau + bumps, multi-start),
        the author's "minimization of general b(r)" read; coarse
        quadrature for the search, full quadrature for the numbers
  then verdicts A/B/C + R-sensitivity + rigid-vs-variational +
  g-flatness + the certified-lattice cross-check (e_parts / kin_of,
  n = 32, L = 48)

Family/conventions: the analytic vacuum-hedgehog family (M5.21.8
builder conventions, s = -1, delta = 0.3, t = 0 evaluation).

Out: ../data/m5_21_14_minimize.json + ../plots/m5_21_14_panel.png
"""
from __future__ import annotations

import importlib.util
import json
import os
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")

_spec8 = importlib.util.spec_from_file_location(
    "b8", os.path.join(HERE, "m5_21_8_b_lattice.py"))
B8 = importlib.util.module_from_spec(_spec8)
_spec8.loader.exec_module(B8)
INS4 = B8.INS4

ETA4 = np.diag([-1.0, 1.0, 1.0, 1.0])
DELTA = 0.3
S_SIGN = -1.0
G_MAIN = 32.0

G1 = np.zeros((4, 4)); G1[2, 3] = -1.0; G1[3, 2] = 1.0
G2 = np.zeros((4, 4)); G2[1, 3] = 1.0; G2[3, 1] = -1.0
G3 = np.zeros((4, 4)); G3[1, 2] = -1.0; G3[2, 1] = 1.0


def rotb(G, A):
    G2m = G @ G
    return (np.eye(4)[None] + np.sin(A)[:, None, None] * G[None]
            + (1 - np.cos(A))[:, None, None] * G2m[None])


def m4h_batch(P, g, t=0.0):
    """full 4x4 analytic base Qh d4 Qh^T at points P (N,3)."""
    X, Y, Z = P[:, 0], P[:, 1], P[:, 2]
    rho = np.sqrt(X * X + Y * Y)
    phi = np.arctan2(Y, X)
    th = -np.arctan2(Z, rho)
    Qh = rotb(G3, phi) @ rotb(G2, th) @ rotb(G1, t * np.ones_like(phi))
    d4 = np.diag([-S_SIGN * g, 1.0, DELTA, 0.0])
    return Qh @ d4[None] @ np.swapaxes(Qh, -1, -2)


def a0h_batch(P, g, dt=1e-5):
    return (m4h_batch(P, g, t=dt) - m4h_batch(P, g, t=-dt)) / (2 * dt)


def kgeom(P):
    """K, K2, r for the boost at points P."""
    r = np.linalg.norm(P, axis=1)
    n = P / r[:, None]
    N = P.shape[0]
    K = np.zeros((N, 4, 4))
    K[:, 0, 1:] = n
    K[:, 1:, 0] = n
    K2 = np.zeros((N, 4, 4))
    K2[:, 0, 0] = 1.0
    K2[:, 1:, 1:] = n[:, :, None] * n[:, None, :]
    return K, K2, r


def qb_from(K, K2, b):
    return (np.eye(4)[None] + np.sinh(b)[:, None, None] * K
            + (np.cosh(b) - 1.0)[:, None, None] * K2)


def dens_u_batch(A):
    tot = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = (A[i] @ ETA4 @ A[j] - A[j] @ ETA4 @ A[i])
            EF = ETA4 @ F @ ETA4
            tot = tot + 4.0 * np.einsum("nab,nab->n", EF, F)
    return tot


def dens_k_batch(a0, A):
    tot = 0.0
    for i in range(3):
        F = a0 @ ETA4 @ A[i] - A[i] @ ETA4 @ a0
        EF = ETA4 @ F @ ETA4
        tot = tot + 4.0 * np.einsum("nab,nab->n", EF, F)
    return tot


# ---------------- quadrature grid ----------------
def make_grid(nr, ndir_u, ndir_phi, rmax=24.0):
    rs = np.geomspace(0.15, rmax, nr)
    u, wu = np.polynomial.legendre.leggauss(ndir_u)
    phis = (np.arange(ndir_phi) + 0.5) * 2 * np.pi / ndir_phi
    U, PH = np.meshgrid(u, phis, indexing="ij")
    WD = (np.repeat(wu[:, None], ndir_phi, 1)
          * (2 * np.pi / ndir_phi)).ravel()
    st = np.sqrt(1 - U ** 2)
    dirs = np.stack([st * np.cos(PH), st * np.sin(PH), U],
                    -1).reshape(-1, 3)
    P = (rs[:, None, None] * dirs[None]).reshape(-1, 3)
    wvol = (np.gradient(rs) * rs ** 2)[:, None] * WD[None]
    return {"rs": rs, "dirs": dirs, "P": P, "wvol": wvol.ravel()}


# ---------------- b(r) families ----------------
RHOS = np.geomspace(0.5, 16.0, 9)


def b_of(avec, r):
    """smooth family: plateau tanh(r/2) + 9 bumps."""
    val = avec[0] * np.tanh(r / 2.0)
    for k, rho in enumerate(RHOS):
        val = val + avec[k + 1] * (r / rho) * np.exp(-((r / rho) ** 2))
    return val


def saw_of(A, lam, r):
    return A * np.sin(np.pi * r / lam) * (r <= 8.0)


# ---------------- the exact correction functionals ----------------
RICH_SHIFTS = [(+1, 8.0), (-1, -8.0), (+2, -1.0), (-2, 1.0)]


class ExactCorr:
    """E_u and kin corrections of the b(r)-dressed analytic family.
    All b-independent pieces (the shifted base fields and the boost
    geometry) are precomputed; an evaluation only builds Qb."""

    def __init__(self, grid, g, h=1e-4):
        self.grid, self.g, self.h = grid, g, h
        P = grid["P"]
        self.sets = {}
        for ax in range(3):
            e = np.zeros(3)
            e[ax] = 1.0
            for k, _ in RICH_SHIFTS:
                Q = P + k * h * e
                self.sets[(ax, k)] = (Q, m4h_batch(Q, g)) + kgeom(Q)
        self.K_c, self.K2_c, self.r_c = kgeom(P)
        self.a0_base = a0h_batch(P, g)
        A_base = self._A(lambda r: np.zeros_like(r))
        self.du_base = dens_u_batch(A_base)
        self.dk_base = dens_k_batch(self.a0_base, A_base)

    def _A(self, bfun):
        A = []
        for ax in range(3):
            acc = 0.0
            for k, w in RICH_SHIFTS:
                Q, M4, K, K2, r = self.sets[(ax, k)]
                Qb = qb_from(K, K2, bfun(r))
                acc = acc + w * (Qb @ M4 @ np.swapaxes(Qb, -1, -2))
            A.append(acc / (12.0 * self.h))
        return A

    def densities(self, bfun):
        A = self._A(bfun)
        du = dens_u_batch(A)
        Qb = qb_from(self.K_c, self.K2_c, bfun(self.r_c))
        a0d = Qb @ self.a0_base @ np.swapaxes(Qb, -1, -2)
        dk = dens_k_batch(a0d, A)
        return du - self.du_base, dk - self.dk_base

    def e_corr(self, bfun):
        A = self._A(bfun)
        return float(np.sum(self.grid["wvol"]
                            * (dens_u_batch(A) - self.du_base)))

    def both(self, bfun):
        du, dk = self.densities(bfun)
        w = self.grid["wvol"]
        return float(np.sum(w * du)), float(np.sum(w * dk))

    def cut(self, bfun, rcut):
        du, dk = self.densities(bfun)
        w = np.where(self.r_c <= rcut, self.grid["wvol"], 0.0)
        return (float(np.sum(w * du)), float(np.sum(w * dk)),
                float(np.sum(w * self.dk_base)))


# ---------------- DIAG: the leading-order runaway (record) ----------
def stage_diag():
    grid = make_grid(40, 8, 16)
    P, w = grid["P"], grid["wvol"]

    def pad0(M):
        return np.pad(M, ((0, 0), (1, 0), (1, 0)))

    h = 1e-4
    G = []
    for ax in range(3):
        e = np.zeros(3)
        e[ax] = 1.0
        d = (8.0 * (m4h_batch(P + h * e, 1.0)
                    - m4h_batch(P - h * e, 1.0))
             - (m4h_batch(P + 2 * h * e, 1.0)
                - m4h_batch(P - 2 * h * e, 1.0))) / (12.0 * h)
        G.append(d[:, 1:, 1:])
    rs = grid["rs"]
    r = np.linalg.norm(P, axis=1)
    n = P / r[:, None]

    def t1_of(beta_nodes):
        beta = np.interp(r, rs, beta_nodes)
        bp = np.interp(r, rs, np.gradient(beta_nodes, rs))
        V = []
        for ax in range(3):
            e = np.zeros(3)
            e[ax] = 1.0
            V.append(bp[:, None] * n[:, ax][:, None] * n
                     + (beta / r)[:, None]
                     * (e[None] - n[:, ax][:, None] * n))
        tot = 0.0
        for i in range(3):
            for j in range(i + 1, 3):
                C = G[i] @ G[j] - G[j] @ G[i]
                W = (V[j][:, :, None] * V[i][:, None, :]
                     - V[i][:, :, None] * V[j][:, None, :])
                ww = (np.einsum("nab,nb->na", G[i], V[j])
                      - np.einsum("nab,nb->na", G[j], V[i]))
                tot = tot + 4.0 * (
                    2.0 * np.einsum("nab,nab->n", C, W)
                    + np.einsum("nab,nab->n", W, W)
                    - 2.0 * np.einsum("na,na->n", ww, ww))
        return float(np.sum(w * tot))

    b0 = 0.5 * rs * np.exp(-rs / 4.0)
    res = minimize(lambda v: t1_of(np.concatenate([[0.0], v])),
                   b0[1:], method="L-BFGS-B",
                   options={"maxiter": 120, "ftol": 1e-16})
    prof = np.concatenate([[0.0], res.x])
    flips = int(np.sum(np.sign(prof[1:]) != np.sign(prof[:-1])))
    return {"t1_reached": float(res.fun), "sign_flips": flips,
            "max_abs_beta": float(np.max(np.abs(prof))),
            "note": ("bounded only by the iteration cap; the full "
                     "first-pass record reached -1.04e7 with a "
                     "grid-scale +-11 sawtooth (12 sign flips)"),
            "reading": ("T1_static alone is UNBOUNDED BELOW: the "
                        "pure-radial UV channel kills the quartic "
                        "|W|^2 while -2|G_i v_j - G_j v_i|^2 grows "
                        "as beta'^2; the leading term needs its "
                        "g-suppressed stabilizers or a constraint")}


# ---------------- stages ----------------
def stage_bnd(ec):
    out = []
    for lam in (2.0, 1.0, 0.5):
        row = {"lam": lam, "amps": [], "E": []}
        for A in (0.02, 0.05, 0.1, 0.2, 0.4, 0.8):
            e = ec.e_corr(lambda r, A=A, lam=lam: saw_of(A, lam, r))
            row["amps"].append(A)
            row["E"].append(e)
        row["floor_found"] = bool(np.argmin(row["E"])
                                  < len(row["E"]) - 1)
        out.append(row)
        print(json.dumps({"BND": row}), flush=True)
    return out


def stage_min(ec_opt):
    best = None
    for start in ([0.026] + [0.0] * 9, [0.0] * 10,
                  [0.01] + [0.02] * 9):
        res = minimize(
            lambda a: ec_opt.e_corr(lambda r: b_of(np.asarray(a), r)),
            np.array(start), method="L-BFGS-B",
            options={"maxiter": 120, "ftol": 1e-13, "eps": 1e-5})
        print(json.dumps({"MIN_start": start[:2],
                          "E_coarse": float(res.fun)}), flush=True)
        if best is None or res.fun < best.fun:
            best = res
    return np.asarray(best.x)


def main():
    t0 = time.time()
    out = {"DIAG": stage_diag()}
    print(json.dumps({"DIAG": out["DIAG"]}), flush=True)

    grid_opt = make_grid(48, 8, 16)
    grid_full = make_grid(72, 12, 24)
    ec_opt = ExactCorr(grid_opt, G_MAIN)
    ec = ExactCorr(grid_full, G_MAIN)
    out["BND"] = stage_bnd(ec)

    avec = stage_min(ec_opt)
    e_corr, k_corr = ec.both(lambda r: b_of(avec, r))
    print(json.dumps({"MIN_full": {"E_corr": e_corr,
                                   "kin_corr": k_corr}}), flush=True)
    rs = grid_full["rs"]
    bstar = b_of(avec, rs)

    # rigid constant-b scan on the same (full) functional
    consts = np.linspace(0.005, 0.06, 12)
    e_c = [ec.e_corr(lambda r, c=c: np.full_like(r, c))
           for c in consts]
    ic = int(np.argmin(e_c))
    rigid = {"b_const_best": float(consts[ic]),
             "E_rigid": float(e_c[ic]),
             "consts": consts.tolist(), "E_scan": e_c}
    print(json.dumps({"rigid": {k: rigid[k] for k in
                                ("b_const_best", "E_rigid")}}),
          flush=True)

    # R-sensitivity (kin base is ~R-extensive)
    sens = {}
    for rcut in (8.0, 12.0, 16.0, 24.0):
        du, dk, k0 = ec.cut(lambda r: b_of(avec, r), rcut)
        sens[f"R{rcut:g}"] = {"E_corr": du, "kin_corr": dk,
                              "kin_base": k0, "kin_total": k0 + dk}
    print(json.dumps({"R_sensitivity": sens}), flush=True)

    # g-flatness: beta* = g b* carried across g
    beta_star = G_MAIN * bstar
    flat = {}
    for g in (8.0, 16.0):
        ecg = ExactCorr(grid_full, g)
        eg, kg = ecg.both(lambda r: np.interp(r, rs, beta_star) / g)
        flat[f"g{g:g}"] = {"E_corr": eg, "kin_corr": kg}
    flat[f"g{G_MAIN:g}"] = {"E_corr": e_corr, "kin_corr": k_corr}
    print(json.dumps({"g_flatness": flat}), flush=True)

    # lattice cross-check (certified instruments)
    cfg = INS4.base_cfg(s=S_SIGN, g=G_MAIN, n=32, L=48.0)
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS4.coords(n, h)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    bl = np.interp(R.ravel(), rs, bstar).reshape(R.shape)
    nx, ny, nz = X / R, Y / R, Z / R
    K = np.zeros(X.shape + (4, 4))
    K[..., 0, 1], K[..., 0, 2], K[..., 0, 3] = nx, ny, nz
    K[..., 1, 0], K[..., 2, 0], K[..., 3, 0] = nx, ny, nz
    K2 = np.zeros_like(K)
    K2[..., 0, 0] = 1.0
    for i, a in enumerate((nx, ny, nz)):
        for j, bb in enumerate((nx, ny, nz)):
            K2[..., 1 + i, 1 + j] = a * bb
    Qb = (np.eye(4)[None, None, None] + np.sinh(bl)[..., None, None]
          * K + (np.cosh(bl) - 1.0)[..., None, None] * K2)
    Mb = B8.dressed(cfg, 0.0)
    a0b = B8.a0_unit(cfg, 0.0)
    Md = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad", Qb, Mb, Qb))
    a0d = INS4.sym4(np.einsum("...ab,...bc,...dc->...ad",
                              Qb, a0b, Qb))
    eu0, ev0 = INS4.e_parts(Mb, cfg)
    eud, evd = INS4.e_parts(Md, cfg)
    k0l = INS4.kin_of(Mb, a0b, cfg)
    kdl = INS4.kin_of(Md, a0d, cfg)
    lat = {"E_corr_lattice": float(eud + evd - eu0 - ev0),
           "kin_corr_lattice": float(kdl - k0l),
           "kin_base_lattice": float(k0l),
           "E_v_dressed": float(evd),
           "E_rel_dev": float(abs((eud + evd - eu0 - ev0) - e_corr)
                              / abs(e_corr)),
           "kin_rel_dev": float(abs((kdl - k0l) - k_corr)
                                / abs(k_corr))}
    print(json.dumps({"lattice": lat}), flush=True)

    verdicts = {
        "A_variational_gain_negative": bool(e_corr < 0),
        "E_corr_at_bstar": e_corr,
        "E_rigid_best": rigid["E_rigid"],
        "variational_vs_rigid": e_corr - rigid["E_rigid"],
        "B_kin_term_negative_derived": True,
        "kin_corr_at_bstar": k_corr,
        "B_kin_corr_negative_measured": bool(k_corr < 0),
        "C_kin_total_by_R": {k: bool(v["kin_total"] < 0)
                             for k, v in sens.items()},
        "T1_leading_order_unbounded": True}
    print(json.dumps({"verdicts": verdicts}), flush=True)

    out.update({"rs": rs.tolist(), "b_star": bstar.tolist(),
                "beta_star": beta_star.tolist(),
                "avec": avec.tolist(), "rigid": rigid,
                "R_sensitivity": sens, "g_flatness": flat,
                "lattice_crosscheck": lat, "verdicts": verdicts,
                "gates": {"lattice_E": bool(lat["E_rel_dev"] < 0.10),
                          "lattice_kin":
                              bool(lat["kin_rel_dev"] < 0.10)},
                "runtime_s": round(time.time() - t0, 1)})
    out["all_green"] = all(out["gates"].values())
    with open(os.path.join(DATA, "m5_21_14_minimize.json"), "w") as f:
        json.dump(out, f, indent=1)

    panel(out)
    print(json.dumps({"gates": out["gates"],
                      "all_green": out["all_green"],
                      "runtime_s": out["runtime_s"]}))


def panel(out):
    """the figure, rebuildable from the saved JSONs (mode: panel)."""
    rs = np.array(out["rs"])
    beta_star = np.array(out["beta_star"])
    rigid = out["rigid"]
    sens = out["R_sensitivity"]
    with open(os.path.join(DATA, "m5_21_14_verify.json")) as f:
        ver = json.load(f)
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))
    a = ax[0, 0]
    gs = [1e2, 1e3, 1e4]
    for row in ver["V1"]["rows"]:
        a.loglog(gs, row["errs_u"], "o-", alpha=0.6)
        a.loglog(gs, row["errs_k"], "s--", alpha=0.6)
    a.loglog(gs, [0.5 / g for g in gs], "k:", label="~1/g")
    a.set_xlabel("g")
    a.set_ylabel("rel err vs T1")
    a.set_title("V1: pointwise convergence (o static, s kin)")
    a.legend()
    a = ax[0, 1]
    for c in ver["V2"]["curves"]:
        m = np.array(c["m"]); E = np.array(c["E"])
        a.plot(m, E - E.min(), label=f"g={c['g']:g}")
        mh = np.arctanh(1.0 / c["g"])
        a.axvline(mh, ls=":", alpha=0.4)
        a.axvline(-mh, ls=":", alpha=0.4)
    a.set_xlabel("m"); a.set_ylabel("E - Emin")
    a.set_title("V2: lattice E(m); dotted = artanh(1/g)")
    a.legend()
    a = ax[0, 2]
    for c in ver["V2"]["curves"]:
        m = np.array(c["m"]); E = np.array(c["E"])
        sel = m > 0
        iz = int(np.argmin(np.abs(m)))
        a.plot(c["g"] * m[sel], E[sel] - E[iz], label=f"g={c['g']:g}")
    a.set_xlabel("beta = g*m"); a.set_ylabel("E(m) - E(0)")
    a.set_title("the beta-collapse (flat-gain retrodiction)")
    a.legend()
    a = ax[1, 0]
    a.plot(rs, beta_star, label="beta* = g b* (variational)")
    a.axhline(G_MAIN * rigid["b_const_best"], ls="--",
              label=f"rigid best beta = "
                    f"{G_MAIN * rigid['b_const_best']:.2f}")
    a.set_xlabel("r"); a.set_ylabel("beta")
    a.set_title("S2: the smooth minimizer (g = 32)")
    a.legend()
    a = ax[1, 1]
    for row in out["BND"]:
        a.plot(row["amps"], row["E"], "o-",
               label=f"sawtooth lam={row['lam']:g}")
    a.axhline(0, color="k", lw=0.5)
    a.set_xlabel("amplitude A"); a.set_ylabel("exact E_corr")
    a.set_xscale("log")
    a.set_yscale("symlog", linthresh=100)
    a.set_title("BND: exact functional, UV probe (g = 32)")
    a.legend()
    a = ax[1, 2]
    rcuts = [8.0, 12.0, 16.0, 24.0]
    a.plot(rcuts, [sens[f"R{r:g}"]["kin_base"] for r in rcuts],
           "o-", label="kin base (~R)")
    a.plot(rcuts, [sens[f"R{r:g}"]["kin_corr"] for r in rcuts],
           "s-", label="kin corr (dressing)")
    a.plot(rcuts, [sens[f"R{r:g}"]["kin_total"] for r in rcuts],
           "^-", label="kin total")
    a.axhline(0, color="k", lw=0.5)
    a.set_xlabel("box radius R"); a.set_ylabel("kin")
    a.set_title("C: the oscillation ledger vs box size")
    a.legend()
    fig.suptitle("M5.21.14: the (1/g) dressing term: derived, "
                 "verified, minimized (analytic family, delta = 0.3)")
    fig.tight_layout()
    os.makedirs(PLOTS, exist_ok=True)
    fig.savefig(os.path.join(PLOTS, "m5_21_14_panel.png"), dpi=110)


if __name__ == "__main__":
    import sys
    if sys.argv[1:] == ["panel"]:
        with open(os.path.join(DATA, "m5_21_14_minimize.json")) as f:
            panel(json.load(f))
    else:
        main()
