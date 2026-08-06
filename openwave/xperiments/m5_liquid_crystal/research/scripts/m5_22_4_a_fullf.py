"""M5.22.4 opening add-on: the FULL-F electric instrument.

The author's 2026-08-06 definition (m5_22_convo.md): the basic
electric instrument is the curvature of the longest axis only
("corresponding to Faber's n field", = the M5.22.2 calibrated
dual-curvature read); the FULL one uses the spatial coordinates of
the F tensor of arXiv 2108.07896, E = (F_23, F_31, F_12), "this way
also containing contributions from other eigenvalues, eigenvectors".

The literal build (paper eq 5-8, static 3x3 case):
  O(x)    the oriented full eigenframe of M (columns e1 short,
          e2 middle, e3 long; right-handed by construction)
  Gamma_i = O^T d_i O, antisymmetric; the rotation vector
  Gvec_i  = ((Gamma_i)_32, (Gamma_i)_13, (Gamma_i)_21)   (eq 7)
  Rvec_ij = Gvec_i x Gvec_j                              (eq 3/8)
  E_k     = the space-space (dual-pair) components: the stacking
            mirrors mermin_B exactly: E = (R_23, -R_13, R_12).

Rvec_ij is a vector in the INTERNAL frame. Scalar reads run as a
measured fork (exhaustion rule):
  comp3     internal long-axis component. DERIVED IDENTITY (verified
            numerically here, not assumed): (Gvec_i x Gvec_j)_3
            = e3.(d_i e3 x d_j e3), i.e. the full F CONTAINS the
            basic instrument as its long-axis component. Proof:
            (Gvec_i)_1 (Gvec_j)_2 - (1<->2) with (Gamma)_32 =
            e3.d e2 = -e2.d e3 expands to the e1/e2 determinant of
            d e3, which is e3.(d_i e3 x d_j e3) since d e3 _|_ e3.
  comp1/2   internal short/middle components = by the same identity
            the Mermin-Ho densities of e1/e2: the objects whose flux
            the M5.22.2 axis fork already measured (<= 0.01/0.34).
  norm3     sign(comp3) * ||Rvec||: the all-contributions magnitude
            read, labeled extra arm.

Frame note: e3 comes from the audited orient_v1 continuity fix; e1
from the same machinery via eigenvalue remap (m5_22_2_a_dive
orient_axis); e2 = e3 x e1 completes the right-handed frame, so O is
orientation-consistent wherever both fixes are.

Modes:
  calib     analytic hedgehog (+ known-charge states): identity diff
            + flux table per read
  all       the four M5.22.4 targets + the two calibration states
Out: ../data/m5_22_4_fullf_<mode>.json
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


DIVE = _load("dive", "m5_22_2_a_dive.py")
INS = DIVE.INS
PAIR = DIVE.PAIR
W2_T2 = DIVE.W2_T2


# ================= the full frame =================
def full_frame(M):
    """oriented right-handed eigenframe: e3 = orient_v1 long axis,
    e1 = orient_axis short axis, e2 = e3 x e1. Returns (O, ncf3, ncf1)
    with O[..., :, a] = e_{a+1}."""
    e3, ncf3 = PAIR.orient_v1(M)
    e1, ncf1 = DIVE.orient_axis(M, 0)
    # exact orthonormalization guard (eigh vectors are orthonormal to
    # float eps; remove the residual e3-component then renormalize)
    e1 = e1 - np.einsum("...a,...a->...", e1, e3)[..., None] * e3
    e1 = e1 / np.linalg.norm(e1, axis=-1, keepdims=True)
    e2 = np.cross(e3, e1)
    O = np.stack([e1, e2, e3], axis=-1)
    return O, int(ncf3), int(ncf1)


def gamma_vecs(O, h):
    """Gvec_i = ((Gamma_i)_32, (Gamma_i)_13, (Gamma_i)_21) with
    Gamma_i = O^T d_i O (central differences, boundary zeroed to
    match PAIR._central)."""
    out = []
    for ax in range(3):
        dO = PAIR._central(O, ax, h)
        Gm = np.einsum("...ka,...kb->...ab", O, dO)
        out.append(np.stack([Gm[..., 2, 1], Gm[..., 0, 2],
                             Gm[..., 1, 0]], axis=-1))
    return out


def full_F(M, h):
    """the three internal-vector curvatures R_23, R_13, R_12 and the
    assembled per-read E fields (dict read -> (n,n,n,3) E array,
    stacking mirroring mermin_B: E = (R_23, -R_13, R_12))."""
    O, ncf3, ncf1 = full_frame(M)
    G = gamma_vecs(O, h)
    R23 = np.cross(G[1], G[2])
    R13 = np.cross(G[0], G[2])
    R12 = np.cross(G[0], G[1])

    def assemble(comp):
        return np.stack([R23[..., comp], -R13[..., comp],
                         R12[..., comp]], axis=-1)

    E3 = assemble(2)
    E2 = assemble(1)
    E1 = assemble(0)
    nrm = np.stack([np.linalg.norm(R23, axis=-1),
                    np.linalg.norm(R13, axis=-1),
                    np.linalg.norm(R12, axis=-1)], axis=-1)
    En = np.sign(E3) * nrm
    return {"comp3": E3, "comp2": E2, "comp1": E1, "norm3": En}, \
        {"orient_conflicts_long": ncf3, "orient_conflicts_short": ncf1}


def reads(M, cfg):
    """per-read flux ladder + the identity diff vs the basic
    instrument (mermin_B of the oriented long axis)."""
    h = cfg["h"]
    fields, meta = full_F(M, h)
    e3, _ = PAIR.orient_v1(M)
    basic = PAIR.mermin_B(e3, h)
    interior = ~INS.pin_shell(cfg["n"], h)
    d = (fields["comp3"] - basic)[interior]
    b = basic[interior]
    out = {"meta": meta,
           "identity_diff": {
               "max_abs": float(np.abs(d).max()),
               "rms": float(np.sqrt((d * d).mean())),
               "basic_rms": float(np.sqrt((b * b).mean()))}}
    for name, E in fields.items():
        out[name] = {f"half{hh:g}": float(
            PAIR.cube_flux(E, cfg, 0.0, hh))
            for hh in (6.0, 12.0, 18.0)}
    out["basic"] = {f"half{hh:g}": float(
        PAIR.cube_flux(basic, cfg, 0.0, hh))
        for hh in (6.0, 12.0, 18.0)}
    return out


def cfg_for(z):
    return INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                        n=int(z["n"]), delta=float(z["delta"]),
                        bc="pinned")


# ================= modes =================
def calib(_kw):
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=32, delta=0.3, bc="pinned")
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rs = np.where(r < 1e-9, 1e-9, r)
    nhat = np.stack([X / rs, Y / rs, Z / rs], axis=-1)
    # uniaxial hedgehog as a NEARLY-DEGENERATE biaxial M: long axis
    # rhat, weak transverse split so the full frame is defined
    t1 = np.stack([-Y, X, np.zeros_like(X)], axis=-1)
    t1n = np.linalg.norm(t1, axis=-1, keepdims=True)
    t1 = np.where(t1n > 1e-9, t1 / np.where(t1n < 1e-9, 1.0, t1n),
                  np.stack([np.ones_like(X), np.zeros_like(X),
                            np.zeros_like(X)], axis=-1))
    M = (np.einsum("...a,...b->...ab", nhat, nhat)
         + 0.3 * np.einsum("...a,...b->...ab", t1, t1))
    out = {"hedgehog": reads(M, cfg)}
    for tag, known in [("P-0.5_plane_sc6_n32_pinned_d0.3", -1.0),
                       ("E+0.5_plane_sc6_n32_pinned_d0.3", -1.0)]:
        z = np.load(os.path.join(DATA, f"m5_22_end_{tag}.npz"))
        row = reads(z["M"].astype(np.float64), cfg_for(z))
        row["known_topological_q"] = known
        out[tag] = row
        print(json.dumps({"tag": tag,
                          "comp3_half12": row["comp3"]["half12"],
                          "basic_half12": row["basic"]["half12"],
                          "iden_max": row["identity_diff"]["max_abs"]}))
    with open(os.path.join(DATA, "m5_22_4_fullf_calib.json"),
              "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out["hedgehog"], indent=1))
    return out


TARGETS = [
    ("m5_22", "P-0.5_plane_sc6_n32_pinned_d0.3", "proton-analog"),
    ("m5_22", "P-1_plane_sc6_n32_pinned_d0.3", "neutron-analog"),
    ("m5_22_1", "d2_s-0.5_s-0.5_a2_sc6_n32_d0.3", "d2 neutral basin"),
    ("m5_22_1", "dn_-0.5at-2_+0.5at+0_-0.5at+2_n32_d0.3",
     "deuteron candidate"),
]


def all_mode(_kw):
    out = {}
    for src, tag, note in TARGETS:
        z = np.load(os.path.join(DATA, f"{src}_end_{tag}.npz"))
        row = reads(z["M"].astype(np.float64), cfg_for(z))
        row["note"] = note
        out[tag] = row
        print(json.dumps({"tag": tag, "note": note,
                          "comp3_half18": row["comp3"]["half18"],
                          "basic_half18": row["basic"]["half18"],
                          "comp2_half18": row["comp2"]["half18"],
                          "comp1_half18": row["comp1"]["half18"],
                          "norm3_half18": row["norm3"]["half18"],
                          "iden_max": row["identity_diff"]["max_abs"]}))
    with open(os.path.join(DATA, "m5_22_4_fullf_all.json"), "w") as f:
        json.dump(out, f, indent=1)
    return out


if __name__ == "__main__":
    ARGV = sys.argv[1:]
    mode = ARGV[0]
    kw = dict(a.split("=", 1) for a in ARGV[1:])
    if mode == "calib":
        calib(kw)
    elif mode == "all":
        all_mode(kw)
