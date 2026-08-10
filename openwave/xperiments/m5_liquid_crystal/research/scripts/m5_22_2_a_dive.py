"""M5.22.2 opening move: the div E electric instrument (both variants).

The author's 2026-08-02 correction (m5_22_convo.md): the M5.22.1 note's
charge density rho = div B / 4pi (B = oriented Mermin-Ho flux of the
long axis) "seems magnetic charge"; the ELECTRIC charge needs div E,
"for E being e.g. curvature of long axis here, or full as in
arXiv 2108.07896".

Two variants, both built and cross-checked (series exhaustion rule):

  E_curv  the literal reading: E = curvature of the long-axis field
          lines, E_i = (nhat . grad) nhat_i. Lift-invariant (n -> -n
          leaves it unchanged), so it needs no orientation fix-up.
  E_full  the paper form (arXiv 2108.07896 eq 3-4): Gamma_mu =
          n x d_mu n, R_munu = Gamma_mu x Gamma_nu = [n.(d_mu n x
          d_nu n)] n, and the DUAL tensor maps the space-space
          curvature to the ELECTRIC sector. The dual vector
          E_k = eps_kij n.(d_i n x d_j n) / 2 is numerically the
          m5_21_4 oriented Mermin-Ho vector: the SAME array the
          M5.22.1 moments called "B". For the ideal hedgehog
          n = rhat it is EXACTLY rhat/r^2 (Coulomb), which the
          calib mode verifies. The correction therefore lands as:
          the old instrument's numbers belong to the ELECTRIC
          sector IF the paper's dual mapping is adopted, and the
          discriminator between the two readings is measured, not
          assumed.

Charges are reported three ways per variant: q_vol = sum div E /4pi
over the interior (pin shell excluded), q_flux(half) = surface flux
through centered cubes (the div read underreads ~11% at ring cores,
note section 5b), and the moments p_z / Q2_zz / Q2_xx of rho.

Convention: signed values are MATHEMATICAL topological orientation;
the electric reading NEGATES them (the author's 2026-07-30 hedgehog
convention, census note section 7). E_curv carries its own sign
meaning (reported raw).

Modes:
  calib                      analytic hedgehog + escaped-ring controls
  read tag=... src=m5_22     both-variant charges + moments, one state
  all                        the pre-registered target set (see TARGETS)
Out: ../data/m5_22_2_dive_<mode>.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def _load(name, fname):
    sp = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


CEN = _load("cen", "m5_22_b_census.py")
INS = CEN.INS
PAIR = CEN.PAIR
W2_T2 = CEN.W2_T2


# ================= the two E fields =================
def e_full(nhat, h):
    """arXiv 2108.07896 eq 3-4 dual curvature vector = the oriented
    Mermin-Ho vector (the array the M5.22.1 moments called B)."""
    return PAIR.mermin_B(nhat, h)


def e_curv(nhat, h):
    """field-line curvature of the long axis: E_i = (n.grad)n_i.
    Lift-invariant: n -> -n gives (+n.grad)(+n). Central differences,
    one-cell boundary zeroed (matches PAIR._central)."""
    E = np.zeros_like(nhat)
    for j in range(3):
        dj = PAIR._central(nhat, j, h)
        E += nhat[..., j:j + 1] * dj
    return E


def divergence(F, h):
    div = np.zeros(F.shape[:3])
    for ax in range(3):
        d = np.zeros(F.shape[:3])
        sl = [slice(None)] * 3
        sp, sm = list(sl), list(sl)
        sl[ax] = slice(1, -1)
        sp[ax] = slice(2, None)
        sm[ax] = slice(None, -2)
        d[tuple(sl)] = (F[tuple(sp) + (ax,)]
                        - F[tuple(sm) + (ax,)]) / (2 * h)
        div += d
    return div


def variant_reads(E, cfg, label):
    """q_vol + flux ladder + moments of rho = div E / 4pi."""
    n, h = cfg["n"], cfg["h"]
    rho = divergence(E, h) / (4.0 * np.pi)
    X, Y, Z = INS.coords(n, h)
    interior = ~INS.pin_shell(n, h)
    r2 = X * X + Y * Y + Z * Z
    w = rho * interior
    h3 = h ** 3
    halves = [6.0, 12.0, 18.0]
    flux = {f"half{hh:g}": float(PAIR.cube_flux(E, cfg, 0.0, hh))
            for hh in halves}
    return {"variant": label,
            "q_vol": float(w.sum() * h3),
            "q_flux": flux,
            "p_z": float((w * Z).sum() * h3),
            "Q2_zz": float((w * (3 * Z * Z - r2)).sum() * h3),
            "Q2_xx": float((w * (3 * X * X - r2)).sum() * h3)}


def both_variants(M, cfg):
    nhat, ncf = PAIR.orient_v1(M)
    out = {"orient_conflicts": int(ncf)}
    out["full"] = variant_reads(e_full(nhat, cfg["h"]), cfg, "full")
    out["curv"] = variant_reads(e_curv(nhat, cfg["h"]), cfg, "curv")
    return out


# ================= calibration =================
def calib(_kw):
    """controls with KNOWN electric content:
    h1  analytic hedgehog n = rhat (charge +1 topological): E_full
        must reproduce Coulomb rhat/r^2 pointwise and q = +1;
        E_curv is identically 0 there (straight field lines).
    h2  the escaped ring texture is NOT constructed analytically;
        instead the relaxed census lepton E+0.5 and proton-analog
        P-0.5 endpoints (known electric -/+1 magnitude, |q_far| = 1
        by the audit exact degrees) are read with both variants."""
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=32, delta=0.3, bc="pinned")
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rs = np.where(r < 1e-9, 1e-9, r)
    nhat = np.stack([X / rs, Y / rs, Z / rs], axis=-1)
    Ef = e_full(nhat, h)
    Ec = e_curv(nhat, h)
    coul = nhat / (rs * rs)[..., None]
    mid = (r > 4.0) & (r < 16.0)
    rel = (np.linalg.norm(Ef - coul, axis=-1)
           / np.linalg.norm(coul, axis=-1))[mid]
    out = {"hedgehog": {
        "full_q_flux_half12": float(PAIR.cube_flux(Ef, cfg, 0.0, 12.0)),
        "full_vs_coulomb_rel_median": float(np.median(rel)),
        "full_vs_coulomb_rel_p90": float(np.percentile(rel, 90)),
        "curv_max_abs": float(np.abs(Ec).max()),
        "curv_q_flux_half12": float(PAIR.cube_flux(Ec, cfg, 0.0, 12.0))}}
    for tag, known in [("P-0.5_plane_sc6_n32_pinned_d0.3", -1.0),
                       ("E+0.5_plane_sc6_n32_pinned_d0.3", -1.0)]:
        z = np.load(os.path.join(DATA, f"m5_22_end_{tag}.npz"))
        c2 = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                          n=int(z["n"]), delta=float(z["delta"]),
                          bc="pinned")
        r2 = both_variants(z["M"].astype(np.float64), c2)
        r2["known_topological_q"] = known
        out[tag] = r2
    with open(os.path.join(DATA, "m5_22_2_dive_calib.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
    return out


# ================= state reads =================
TARGETS = [
    # (src, tag, note)
    ("m5_22", "P-0.5_plane_sc6_n32_pinned_d0.3", "proton-analog n32"),
    ("m5_22", "P-0.5_plane_sc6_n48_pinned_d0.3_ext", "proton-analog n48"),
    ("m5_22", "P-1_plane_sc6_n32_pinned_d0.3", "neutron-analog n32"),
    ("m5_22", "P-1_plane_sc6_n48_pinned_d0.3", "neutron-analog n48"),
    ("m5_22", "E+0.5_plane_sc6_n32_pinned_d0.3", "lepton ref n32"),
    ("m5_22_1", "d2_s-0.5_s-0.5_a2_sc6_n32_d0.3", "pp-control cousin"),
    ("m5_22_1", "dn_-0.5at-2_+0.5at+0_-0.5at+2_n32_d0.3",
     "deuteron candidate n32"),
    ("m5_22_1", "dn_-0.5at-2_+0.5at+0_-0.5at+2_n48_d0.3_ext",
     "deuteron candidate n48"),
]


def read_one(src, tag):
    z = np.load(os.path.join(DATA, f"{src}_end_{tag}.npz"))
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=int(z["n"]), delta=float(z["delta"]),
                       bc="pinned")
    return both_variants(z["M"].astype(np.float64), cfg)


def all_mode(_kw):
    out = {}
    for src, tag, note in TARGETS:
        row = read_one(src, tag)
        row["note"] = note
        out[tag] = row
        print(json.dumps({"tag": tag, "note": note,
                          "full_q_vol": row["full"]["q_vol"],
                          "full_Q2_zz": row["full"]["Q2_zz"],
                          "curv_q_vol": row["curv"]["q_vol"],
                          "curv_Q2_zz": row["curv"]["Q2_zz"]}))
    with open(os.path.join(DATA, "m5_22_2_dive_all.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


# ================= the frame-axis fork =================
def orient_axis(M, idx):
    """oriented eigenvector #idx (ascending eigh order; 2 = long
    axis): remap the eigenvalues so axis idx leads, then reuse the
    audited orient_v1 continuity machinery unchanged."""
    lam, vec = np.linalg.eigh(M)
    w = np.zeros_like(lam)
    order = [i for i in range(3) if i != idx] + [idx]
    for rank, i in enumerate(order):
        w[..., i] = float(rank)
    M2 = np.einsum("...ai,...i,...bi->...ab", vec, w, vec)
    return PAIR.orient_v1(M2)


def axes_mode(_kw):
    """which frame axis carries the quantized Gauss charge: the
    e_full flux per eigenframe axis (0 = short, 1 = middle, 2 = long)
    on the calibration + target states. The arXiv 2108.07896 eq 16
    biaxial decomposition question (R^ee vs R^gg sectors) made
    machine-checkable."""
    picks = [("m5_22", "P-0.5_plane_sc6_n32_pinned_d0.3"),
             ("m5_22", "P-1_plane_sc6_n32_pinned_d0.3"),
             ("m5_22", "E+0.5_plane_sc6_n32_pinned_d0.3"),
             ("m5_22_1", "dn_-0.5at-2_+0.5at+0_-0.5at+2_n32_d0.3")]
    out = {}
    for src, tag in picks:
        z = np.load(os.path.join(DATA, f"{src}_end_{tag}.npz"))
        cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0,
                           w2=W2_T2, n=int(z["n"]),
                           delta=float(z["delta"]), bc="pinned")
        M = z["M"].astype(np.float64)
        row = {}
        for idx, name in [(2, "long"), (1, "middle"), (0, "short")]:
            nhat, ncf = orient_axis(M, idx)
            E = e_full(nhat, cfg["h"])
            row[name] = {
                "q_flux_half12": float(
                    PAIR.cube_flux(E, cfg, 0.0, 12.0)),
                "q_flux_half18": float(
                    PAIR.cube_flux(E, cfg, 0.0, 18.0)),
                "orient_conflicts": int(ncf)}
        out[tag] = row
        print(json.dumps({"tag": tag} | {
            k: round(v["q_flux_half18"], 3)
            for k, v in row.items()}))
    with open(os.path.join(DATA, "m5_22_2_dive_axes.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    return out


if __name__ == "__main__":
    ARGV = sys.argv[1:]
    mode = ARGV[0]
    kw = dict(a.split("=", 1) for a in ARGV[1:])
    if mode == "calib":
        calib(kw)
    elif mode == "read":
        row = read_one(kw.get("src", "m5_22"), kw["tag"])
        print(json.dumps(row, indent=1))
    elif mode == "all":
        all_mode(kw)
    elif mode == "axes":
        axes_mode(kw)
