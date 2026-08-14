"""M5.21.16 B: the FLIP functional on the certified lattice stack.

FLIP = keep comm_eta brackets and V4 unchanged; replace the curvature
contraction <F,F>_eta = tr(eta F eta F^T) by the Frobenius contraction
tr(F F^T). Licensed by the arm-A bridge lemma: this reverses exactly the
(tilde-Gamma)^2 boost-row contributions (the author's variant A at
leading order).

Stages:
  IDENT  static 3x3-embedded field: flip == eta identically (the charge/
         Coulomb sector is untouched by construction)
  INV    invariance table: E_u under SO(3) rotation, SO(1,3) boost, and
         compact SO(4)-style (0,1)-plane rotation conjugations, for eta
         and flip functionals (random smooth 4x4 field)
  CHAN   per-channel kin on the analytic boost-hedgehog family evaluated
         on the certified lattice (n = 32, L = 48, g = 32, s = -1):
         eta vs flip for the M5.21.3 generator catalog
  DRESS  the M5.21.14 dressed family (quadrature): kin correction and
         the sawtooth static-runaway probe, eta vs flip

Out: ../data/m5_21_16_field.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import time

import numpy as np
from scipy.linalg import expm
from scipy.ndimage import gaussian_filter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

_s3 = importlib.util.spec_from_file_location(
    "b3", os.path.join(HERE, "m5_21_3_a_4d.py"))
B3 = importlib.util.module_from_spec(_s3)
_s3.loader.exec_module(B3)

_s14 = importlib.util.spec_from_file_location(
    "c14", os.path.join(HERE, "m5_21_14_c_minimize.py"))
C14 = importlib.util.module_from_spec(_s14)
_s14.loader.exec_module(C14)

ETA = B3.ETA


def inner_frob(F, G):
    return np.einsum("...ab,...ab->...", F, G, optimize=True)


def e_u_of(M, cfg, contraction):
    """spatial curvature energy with a chosen contraction."""
    inner = B3.inner_eta if contraction == "eta" else inner_frob
    h3 = cfg["h"] ** 3
    e_u = 0.0
    for br, (A, wt) in B3.a_fields(M, cfg).items():
        for i in range(3):
            for j in range(i + 1, 3):
                F = B3.comm_eta(A[i], A[j])
                e_u = e_u + wt * 4.0 * np.sum(inner(F, F))
    return h3 * e_u


def kin_of(M, a0, cfg, contraction):
    inner = B3.inner_eta if contraction == "eta" else inner_frob
    h3 = cfg["h"] ** 3
    k = 0.0
    for br, (A, wt) in B3.a_fields(M, cfg).items():
        for i in range(3):
            F = B3.comm_eta(a0, A[i])
            k = k + wt * 4.0 * np.sum(inner(F, F))
    return h3 * k


def dens_u_flip(A):
    tot = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = (A[i] @ C14.ETA4 @ A[j] - A[j] @ C14.ETA4 @ A[i])
            tot = tot + 4.0 * np.einsum("nab,nab->n", F, F)
    return tot


def dens_k_flip(a0, A):
    tot = 0.0
    for i in range(3):
        F = a0 @ C14.ETA4 @ A[i] - A[i] @ C14.ETA4 @ a0
        tot = tot + 4.0 * np.einsum("nab,nab->n", F, F)
    return tot


def lattice_family_M(cfg, g):
    """the C14 analytic boost-hedgehog family on the lattice grid."""
    X, Y, Z = B3.coords(cfg["n"], cfg["h"])
    P = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    # avoid the axis singularity at rho = 0 by an epsilon shift
    P[:, 0] += 1e-9
    M = C14.m4h_batch(P, g)
    n = cfg["n"]
    return M.reshape(n, n, n, 4, 4)


def stage_ident(rng):
    cfg = B3.base_cfg(n=16, L=24.0)
    M3 = np.stack([[gaussian_filter(rng.normal(size=(16,) * 3), 2.0)
                    for _ in range(3)] for _ in range(3)], axis=-1)
    M3 = M3.reshape(16, 16, 16, 3, 3)
    M3 = 0.5 * (M3 + M3.swapaxes(-1, -2))
    M4 = B3.embed34(M3, cfg)
    ue = e_u_of(M4, cfg, "eta")
    uf = e_u_of(M4, cfg, "flip")
    rel = abs(ue - uf) / max(abs(ue), 1e-300)
    return {"u_eta": float(ue), "u_flip": float(uf),
            "rel_diff": float(rel), "pass": bool(rel < 1e-13)}


def stage_inv(rng):
    cfg = B3.base_cfg(n=14, L=21.0)
    MB = np.stack([[gaussian_filter(rng.normal(size=(14,) * 3), 2.0)
                    for _ in range(4)] for _ in range(4)], axis=-1)
    MB = MB.reshape(14, 14, 14, 4, 4)
    MB = B3.vac4(cfg)[None, None, None] + 0.5 * B3.sym4(MB)
    out = {}
    gens = {}
    Grot = np.zeros((4, 4)); Grot[1, 2], Grot[2, 1] = -0.3, 0.3
    gens["so3_rot"] = Grot
    Gboo = np.zeros((4, 4)); Gboo[0, 1] = Gboo[1, 0] = 0.25
    gens["so13_boost"] = Gboo
    Gcmp = np.zeros((4, 4)); Gcmp[0, 1], Gcmp[1, 0] = -0.25, 0.25
    gens["so4_compact"] = Gcmp
    for contraction in ("eta", "flip"):
        E0 = e_u_of(MB, cfg, contraction)
        row = {}
        for nm, Gm in gens.items():
            L = expm(Gm)
            ML = np.einsum("ab,...bc,dc->...ad", L, MB, L)
            row[nm] = float(abs(e_u_of(ML, cfg, contraction) - E0)
                            / max(abs(E0), 1e-300))
        out[contraction] = row
    out["pass"] = bool(
        out["eta"]["so3_rot"] < 1e-9 and out["eta"]["so13_boost"] < 1e-9
        and out["flip"]["so3_rot"] < 1e-9
        and out["flip"]["so13_boost"] > 1e-6)
    out["reading"] = (
        "eta: SO(1,3)-invariant (certified); flip: SO(3) kept, boost "
        "invariance BROKEN (the signature surgery costs Lorentz "
        "invariance of the energy; the compact (0,1) rotation is not a "
        "symmetry of either since the bracket keeps eta)")
    return out


def stage_chan():
    cfg = B3.base_cfg(s=-1.0, g=32.0, n=32, L=48.0)
    M = lattice_family_M(cfg, 32.0)
    a0s = B3.gen_catalog(cfg, M)
    rows = {}
    for nm, a0 in a0s.items():
        ke = kin_of(M, a0, cfg, "eta")
        kf = kin_of(M, a0, cfg, "flip")
        rows[nm] = {"kin_eta": float(ke), "kin_flip": float(kf)}
    boosts = {nm: r for nm, r in rows.items() if nm.startswith("boost")}
    out = {"rows": rows,
           "boost_kin_eta_all_negative": bool(
               max(r["kin_eta"] for r in boosts.values()) < 0),
           "boost_kin_flip_all_positive": bool(
               min(r["kin_flip"] for r in boosts.values()) > 0),
           "all_kin_flip_nonnegative": bool(
               min(r["kin_flip"] for r in rows.values()) >= 0)}
    out["pass"] = bool(out["boost_kin_eta_all_negative"])
    return out


def stage_dress():
    grid = C14.make_grid(40, 6, 12)
    ec = C14.ExactCorr(grid, 32.0)
    w = grid["wvol"]
    out = {}
    # base (undressed) kin density integral, eta vs flip
    A_base = ec._A(lambda r: np.zeros_like(r))
    out["kin_base_eta"] = float(np.sum(w * C14.dens_k_batch(
        ec.a0_base, A_base)))
    out["kin_base_flip"] = float(np.sum(w * dens_k_flip(
        ec.a0_base, A_base)))
    # guard-level smooth dressing (plateau amplitude ladder): static
    # correction E_corr eta vs flip: does the negative well close?
    ladder = {}
    for amp in (0.02, 0.05, 0.1, 0.2, 0.4):
        def bf(r, a=amp):
            return a * np.tanh(r / 2.0)
        A = ec._A(bf)
        du_eta = float(np.sum(w * (C14.dens_u_batch(A) - ec.du_base)))
        du_base_flip = dens_u_flip(A_base)
        du_flip = float(np.sum(w * (dens_u_flip(A) - du_base_flip)))
        Qb = C14.qb_from(ec.K_c, ec.K2_c, bf(ec.r_c))
        a0d = Qb @ ec.a0_base @ np.swapaxes(Qb, -1, -2)
        dk_eta = float(np.sum(w * (C14.dens_k_batch(a0d, A)
                                   - ec.dk_base)))
        dk_flip = float(np.sum(w * (dens_k_flip(a0d, A)
                                    - dens_k_flip(ec.a0_base, A_base))))
        ladder[f"amp_{amp:g}"] = {
            "E_corr_eta": du_eta, "E_corr_flip": du_flip,
            "kin_corr_eta": dk_eta, "kin_corr_flip": dk_flip}
    out["plateau_ladder"] = ladder
    # sawtooth static-runaway probe (the M5.21.14 BND family)
    saw = {}
    for A_amp in (0.05, 0.1, 0.2, 0.4):
        def bf(r, a=A_amp):
            return C14.saw_of(a, 1.0, r)
        A = ec._A(bf)
        saw[f"A_{A_amp:g}"] = {
            "E_corr_eta": float(np.sum(w * (C14.dens_u_batch(A)
                                            - ec.du_base))),
            "E_corr_flip": float(np.sum(w * (dens_u_flip(A)
                                             - dens_u_flip(A_base))))}
    out["sawtooth_ladder"] = saw
    e_eta = [v["E_corr_eta"] for v in saw.values()]
    e_flip = [v["E_corr_flip"] for v in saw.values()]
    out["saw_eta_goes_negative"] = bool(min(e_eta) < 0)
    out["saw_flip_stays_positive"] = bool(min(e_flip) > 0)
    return out


def main():
    t0 = time.time()
    rng = np.random.default_rng(21160)
    out = {"IDENT": stage_ident(rng)}
    print(json.dumps({"IDENT": out["IDENT"]}), flush=True)
    out["INV"] = stage_inv(rng)
    print(json.dumps({"INV": out["INV"]}), flush=True)
    out["CHAN"] = stage_chan()
    print(json.dumps({"CHAN": out["CHAN"]}, indent=1), flush=True)
    out["DRESS"] = stage_dress()
    print(json.dumps({"DRESS": out["DRESS"]}, indent=1), flush=True)
    out["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(DATA, "m5_21_16_field.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({"runtime_s": out["runtime_s"]}))


if __name__ == "__main__":
    main()
