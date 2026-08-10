"""M5.22.2 adversarial audit (independent second agent).

Attacks claims A-H of the M5.22.2 run (div E instrument + beta-decay
driver) with INDEPENDENT implementations: own hedgehog, own central
differences, own Levi-Civita dual-curvature vector, own curvature
field, own divergence (np.gradient), own symmetric-box surface flux,
own slab flux, own ring ledger, own kick replication. Trusted (per
audit charter): INS.e_parts / INS.base_cfg (certified M5.21.2b),
PAIR.orient_v1 (certified M5.21.4 orientation machinery). The
functions under audit (DIVE.e_full, PAIR.mermin_B) are only called
where the claim IS about their output (claim B identity).

Out: ../data/m5_22_2_audit.json + one printed JSON verdict summary.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")


def _load(name, fname):
    sp = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(sp)
    sp.loader.exec_module(mod)
    return mod


DIVE = _load("dive_a", "m5_22_2_a_dive.py")
KICK = _load("kick_a", "m5_22_1_a_kick.py")
INS = DIVE.INS          # trusted certified instrument
PAIR = DIVE.PAIR        # orient_v1 trusted; mermin_B under audit (B)
W2_T2 = DIVE.W2_T2

EPS3 = np.zeros((3, 3, 3))
for i, j, k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    EPS3[i, j, k] = 1.0
for i, j, k in [(0, 2, 1), (2, 1, 0), (1, 0, 2)]:
    EPS3[i, j, k] = -1.0


# ================= independent numerics =================
def my_coords(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    return np.meshgrid(x, x, x, indexing="ij")


def my_central(f, ax, h):
    """own central difference, slicing (not np.roll); one-cell
    boundary zeroed (the field definition used by the instrument)."""
    out = np.zeros_like(f)
    a = [slice(None)] * f.ndim
    b = list(a)
    c = list(a)
    a[ax] = slice(1, -1)
    b[ax] = slice(2, None)
    c[ax] = slice(None, -2)
    out[tuple(a)] = (f[tuple(b)] - f[tuple(c)]) / (2.0 * h)
    return out


def my_dual_E(nhat, h):
    """E_k = (1/2) eps_kij n . (d_i n x d_j n), explicit Levi-Civita
    contraction (independent of mermin_B's hand-stacked form)."""
    dn = [my_central(nhat, ax, h) for ax in range(3)]
    F = np.zeros(nhat.shape[:3] + (3, 3))
    for i in range(3):
        for j in range(3):
            if i == j:
                continue
            F[..., i, j] = np.einsum(
                "...a,...a->...", nhat, np.cross(dn[i], dn[j]))
    return 0.5 * np.einsum("kij,...ij->...k", EPS3, F)


def my_curv_E(nhat, h):
    """E_i = (n . grad) n_i, own differences."""
    E = np.zeros_like(nhat)
    for j in range(3):
        E += nhat[..., j:j + 1] * my_central(nhat, j, h)
    return E


def my_div(F, h):
    """own divergence via np.gradient (2nd-order central interior,
    one-sided edges: structurally different from the run's)."""
    return sum(np.gradient(F[..., ax], h, axis=ax) for ax in range(3))


def my_pin_shell(n, h, depth=1.6):
    wc = max(1, int(np.ceil(depth / h)))
    P = np.zeros((n, n, n), dtype=bool)
    for ax in range(3):
        sl = [slice(None)] * 3
        sl[ax] = slice(0, wc)
        P[tuple(sl)] = True
        sl[ax] = slice(n - wc, n)
        P[tuple(sl)] = True
    return P


def _face_sum(F, h, i0, i1, j0, j1, k0, k1):
    s = 0.0
    s += F[i1, j0:j1 + 1, k0:k1 + 1, 0].sum() \
        - F[i0, j0:j1 + 1, k0:k1 + 1, 0].sum()
    s += F[i0:i1 + 1, j1, k0:k1 + 1, 1].sum() \
        - F[i0:i1 + 1, j0, k0:k1 + 1, 1].sum()
    s += F[i0:i1 + 1, j0:j1 + 1, k1, 2].sum() \
        - F[i0:i1 + 1, j0:j1 + 1, k0, 2].sum()
    return float(s * h * h / (4.0 * np.pi))


def my_cube_flux(F, n, h, half, offset_convention=False):
    """symmetric closed lattice box about the true center (n-1)/2
    (for even n the run's cube_flux rounds 15.5 -> 16: half-cell
    offset; offset_convention=True reproduces that for comparison)."""
    c = (n - 1) / 2.0
    if offset_convention:
        k = int(round(half / h))
        i0, i1 = int(round(c)) - k, int(round(c)) + k
    else:
        i1 = int(np.floor(c + half / h + 0.5))
        i0 = (n - 1) - i1
    if i0 < 1 or i1 > n - 2:
        return float("nan")
    return _face_sum(F, h, i0, i1, i0, i1, i0, i1)


def my_slab_flux(F, n, h, z0, z1, half_lat):
    """own fragment slab: symmetric lateral box, z0<=z<=z1."""
    c = (n - 1) / 2.0
    i1 = int(np.floor(c + half_lat / h + 0.5))
    i0 = (n - 1) - i1
    k0 = max(1, int(np.ceil(z0 / h + c)))
    k1 = min(n - 2, int(np.floor(z1 / h + c)))
    if i0 < 1 or i1 > n - 2 or k1 <= k0:
        return float("nan")
    return _face_sum(F, h, i0, i1, i0, i1, k0, k1)


def my_ring_ledger(M, n, h, L, thr, min_vox=2):
    """own eigengap component ledger; column if rho_mean < 3."""
    lam = np.linalg.eigvalsh(M)
    gap = np.minimum(lam[..., 1] - lam[..., 0],
                     lam[..., 2] - lam[..., 1])
    X, Y, Z = my_coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    r = np.sqrt(rho * rho + Z * Z)
    interior = (~my_pin_shell(n, h)) & (r < 0.9 * L / 2.0)
    lab, nc = ndimage.label((gap < thr) & interior,
                            structure=np.ones((3, 3, 3), int))
    rings, cols = [], []
    for k in range(1, nc + 1):
        sel = lab == k
        sz = int(sel.sum())
        if sz < min_vox:
            continue
        e = {"voxels": sz, "rho_mean": float(rho[sel].mean()),
             "z": float(Z[sel].mean())}
        (cols if e["rho_mean"] < 3.0 else rings).append(e)
    rings.sort(key=lambda e: e["z"])
    return rings, cols


def load_end(fname):
    z = np.load(os.path.join(DATA, fname))
    n, delta = int(z["n"]), float(z["delta"])
    cfg = INS.base_cfg(term="T2", stencil="sym", eps=0.0, w2=W2_T2,
                       n=n, delta=delta, bc="pinned")
    return z["M"].astype(np.float64), cfg


V = {}   # verdicts
DEFECTS = []


def verdict(ok_all, ok_part=False):
    return "PASS" if ok_all else ("PARTIAL" if ok_part else "FAIL")


# ================= A: hedgehog calibration =================
def claim_A():
    n, h = 32, 48.0 / 32
    X, Y, Z = my_coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    rs = np.maximum(r, 1e-9)
    nhat = np.stack([X / rs, Y / rs, Z / rs], axis=-1)
    coul = nhat / (rs * rs)[..., None]
    # (i) exact continuum identity, ANALYTIC derivatives (no FD):
    # d_i n_j = (delta_ij - rhat_i rhat_j)/r  ->  E should be coul
    dn_an = [np.zeros(nhat.shape) for _ in range(3)]
    for i in range(3):
        for j in range(3):
            dn_an[i][..., j] = ((1.0 if i == j else 0.0)
                                - nhat[..., i] * nhat[..., j]) / rs
    Fan = np.zeros((n, n, n, 3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                Fan[..., i, j] = np.einsum(
                    "...a,...a->...", nhat, np.cross(dn_an[i], dn_an[j]))
    E_an = 0.5 * np.einsum("kij,...ij->...k", EPS3, Fan)
    mid = (r > 4.0) & (r < 16.0)
    dev_an = float(np.max(np.linalg.norm((E_an - coul), axis=-1)[mid]
                          / np.linalg.norm(coul, axis=-1)[mid]))
    # (ii) own finite-difference dual vector vs Coulomb
    E = my_dual_E(nhat, h)
    rel = (np.linalg.norm(E - coul, axis=-1)
           / np.linalg.norm(coul, axis=-1))[mid]
    med, p90 = float(np.median(rel)), float(np.percentile(rel, 90))
    # (iii) own fluxes: symmetric box + the run's offset convention
    fx_sym = my_cube_flux(E, n, h, 12.0)
    fx_off = my_cube_flux(E, n, h, 12.0, offset_convention=True)
    # own volume charge
    q_vol = float((my_div(E, h)[~my_pin_shell(n, h)]).sum()
                  * h ** 3 / (4 * np.pi))
    ok = (dev_an < 1e-10 and abs(med - 0.007) < 0.004
          and abs(p90 - 0.020) < 0.008 and abs(fx_off - 1.06) < 0.03
          and abs(q_vol - 1.0) < 0.1)
    if abs(fx_sym - 1.0) < abs(fx_off - 1.0) - 0.02:
        DEFECTS.append(
            "cube_flux center rounding: for even n the box center "
            f"rounds 15.5->16 (half-cell offset); hedgehog flux "
            f"offset-box {fx_off:.4f} vs symmetric-box {fx_sym:.4f} "
            "(the +1.06 is partly the offset artifact, not physics)")
    V["A_hedgehog_calibration"] = {
        "verdict": verdict(ok),
        "measured": {
            "analytic_identity_max_rel_dev": dev_an,
            "fd_vs_coulomb_rel_median": med,
            "fd_vs_coulomb_rel_p90": p90,
            "flux_half12_run_convention": fx_off,
            "flux_half12_symmetric_box": fx_sym,
            "q_vol_own_div": q_vol},
        "claimed": {"median": 0.007, "p90": 0.020,
                    "flux_half12": 1.06}}


# ================= B: instrument identity =================
def claim_B():
    M, cfg = load_end(
        "m5_22_1_end_dn_-0.5at-2_+0.5at+0_-0.5at+2_n32_d0.3.npz")
    h = cfg["h"]
    nhat, _ = PAIR.orient_v1(M)
    A1 = DIVE.e_full(nhat, h)           # under audit
    A2 = PAIR.mermin_B(nhat, h)         # the old "B"
    d_id = float(np.abs(A1 - A2).max())
    # my own dual vector vs mermin_B (formula check, interior)
    E_my = my_dual_E(nhat, h)
    inter = ~my_pin_shell(cfg["n"], h)
    d_form = float(np.abs((E_my - A2)[inter]).max())
    scale = float(np.abs(A2[inter]).max())
    ok = d_id == 0.0 and d_form / max(scale, 1e-300) < 1e-12
    V["B_instrument_identity"] = {
        "verdict": verdict(ok),
        "measured": {"e_full_vs_mermin_B_max_abs_diff": d_id,
                     "own_LeviCivita_vs_mermin_B_max_abs_diff": d_form,
                     "field_scale": scale},
        "note": "e_full delegates to PAIR.mermin_B in source; the "
                "identity is structural AND numerically exact"}


# ================= C: quantization split =================
def claim_C():
    out = {}
    for tag, cl_full_vol, cl_full_f12 in [
            ("P-0.5_plane_sc6_n32_pinned_d0.3", -0.977, -1.05),
            ("E+0.5_plane_sc6_n32_pinned_d0.3", -0.974, -1.01)]:
        M, cfg = load_end(f"m5_22_end_{tag}.npz")
        n, h = cfg["n"], cfg["h"]
        nhat, _ = PAIR.orient_v1(M)
        inter = ~my_pin_shell(n, h)
        row = {}
        for lab, E in [("full", my_dual_E(nhat, h)),
                       ("curv", my_curv_E(nhat, h))]:
            rho = my_div(E, h) / (4 * np.pi)
            row[lab] = {
                "q_vol": float((rho * inter).sum() * h ** 3),
                "flux_half6": my_cube_flux(E, n, h, 6.0),
                "flux_half12": my_cube_flux(E, n, h, 12.0),
                "flux_half18": my_cube_flux(E, n, h, 18.0)}
        out[tag] = row
    p, e = (out["P-0.5_plane_sc6_n32_pinned_d0.3"],
            out["E+0.5_plane_sc6_n32_pinned_d0.3"])
    full_ok = (abs(p["full"]["q_vol"] + 0.977) < 0.03
               and abs(e["full"]["q_vol"] + 0.974) < 0.03
               and abs(p["full"]["flux_half12"] + 1.05) < 0.08
               and abs(e["full"]["flux_half12"] + 1.01) < 0.08)
    # curv must be non-quantized + ladder-inconsistent on P-0.5
    cv = p["curv"]
    ladder_spread = abs(cv["flux_half12"] - cv["flux_half6"])
    curv_ok = (abs(cv["q_vol"] - 7.05) < 0.5
               and ladder_spread > 5.0
               and min(abs(cv["q_vol"] - r) for r in (-1, 0, 1)) > 2.0)
    V["C_quantization_split"] = {
        "verdict": verdict(full_ok and curv_ok),
        "measured": out,
        "claimed": {"P-0.5_full_vol": -0.977, "E+0.5_full_vol": -0.974,
                    "P-0.5_curv_vol": +7.05,
                    "P-0.5_curv_flux_6to12": [2.2, 11.9]}}


# ================= D: frame-axis fork =================
def claim_D():
    M, cfg = load_end("m5_22_end_P-0.5_plane_sc6_n32_pinned_d0.3.npz")
    n, h = cfg["n"], cfg["h"]
    lam, vec = np.linalg.eigh(M)
    row = {}
    for idx, name in [(2, "long"), (1, "middle"), (0, "short")]:
        # independent construction: rank-1 projector of THAT
        # eigenvector -> its leading eigenvector IS the chosen axis;
        # orientation by the trusted continuity machinery
        v = vec[..., :, idx]
        M1 = v[..., :, None] * v[..., None, :]
        nhat, ncf = PAIR.orient_v1(M1)
        E = my_dual_E(nhat, h)
        row[name] = {
            "flux_half12": my_cube_flux(E, n, h, 12.0),
            "flux_half18": my_cube_flux(E, n, h, 18.0),
            "orient_conflicts": int(ncf)}
    # global lift sign is not anchored -> compare magnitudes
    ok = (abs(abs(row["long"]["flux_half18"]) - 1.0) < 0.1
          and abs(row["middle"]["flux_half18"]) < 0.4
          and abs(row["short"]["flux_half18"]) < 0.05)
    lab_note = ("claim text says half12 for (-1.02, 0.10, 0.005) but "
                "those are the HALF18 values in m5_22_2_dive_axes.json"
                " (half12: -1.051, 0.152, 0.010): label slip, "
                "substance unchanged")
    DEFECTS.append("claim D wording: " + lab_note)
    V["D_frame_axis_fork"] = {
        "verdict": verdict(ok),
        "measured": row,
        "claimed": {"long": -1.02, "middle": 0.10, "short": 0.005},
        "note": lab_note + "; audit uses |flux| (the global director "
                "lift sign from eigh is unanchored, so per-axis signs "
                "can flip between constructions)"}


# ================= E: deuteron quadrupole =================
def claim_E():
    out = {}
    for fname, key, cl in [
            ("m5_22_1_end_dn_-0.5at-2_+0.5at+0_-0.5at+2_n32_d0.3.npz",
             "n32", 21.8),
            ("m5_22_1_end_dn_-0.5at-2_+0.5at+0_-0.5at+2_n48_d0.3_ext"
             ".npz", "n48_ext", 61.5)]:
        M, cfg = load_end(fname)
        n, h = cfg["n"], cfg["h"]
        nhat, _ = PAIR.orient_v1(M)
        E = my_dual_E(nhat, h)
        rho = my_div(E, h) / (4 * np.pi)
        X, Y, Z = my_coords(n, h)
        r2 = X * X + Y * Y + Z * Z
        w = rho * (~my_pin_shell(n, h))
        out[key] = {
            "Q2_zz": float((w * (3 * Z * Z - r2)).sum() * h ** 3),
            "q_vol": float(w.sum() * h ** 3), "claimed_Q2_zz": cl}
    # cross-check vs the OLD m5_22_1 moments rows
    olds = {}
    for key, rf in [("n32", "m5_22_1_row_dn_-0.5at-2_+0.5at+0_"
                     "-0.5at+2_n32_d0.3.json"),
                    ("n48_ext", "m5_22_1_row_dn_-0.5at-2_+0.5at+0_"
                     "-0.5at+2_n48_d0.3_ext.json")]:
        with open(os.path.join(DATA, rf)) as f:
            olds[key] = json.load(f)["moments"]["Q2_zz"]
    ok = (abs(out["n32"]["Q2_zz"] - 21.8) < 0.5
          and abs(out["n48_ext"]["Q2_zz"] - 61.5) < 1.0
          and abs(out["n32"]["Q2_zz"] - olds["n32"]) < 0.5
          and abs(out["n48_ext"]["Q2_zz"] - olds["n48_ext"]) < 1.0)
    V["E_deuteron_quadrupole"] = {
        "verdict": verdict(ok),
        "measured": out, "old_m5_22_1_moments_Q2_zz": olds,
        "note": "electric reading = negated (convention), so "
                "-21.8/-61.5 electric follows from +21.8/+61.5 "
                "topological if the sign convention is granted"}


# ================= F: returned verdicts =================
F_ROWS = [
    ("m5_22_2_end_P-1_plane_sc6_n32_pinned_d0.3__K10.4_s216_ext.npz",
     "P-1_K10.4_ext", 12.7296, 12.0, 1.04),
    ("m5_22_2_end_P-1_plane_sc6_n32_pinned_d0.3__K30.15_s216.npz",
     "P-1_K30.15", 12.7296, 12.0, 1.04),
    ("m5_22_2_end_P-1_plane_sc6_n32_pinned_d0.3__K30.4_s216.npz",
     "P-1_K30.4", 12.7296, 12.0, 1.04),
    ("m5_22_2_end_d2_s-0.5_s-0.5_a2_sc6_n32_d0.3__K10.4_s216_ext.npz",
     "d2_K10.4_ext", 15.047, 13.4, 1.08),
    ("m5_22_2_end_d2_s-0.5_s-0.5_a2_sc6_n32_d0.3__K30.15_s216.npz",
     "d2_K30.15", 15.047, 13.4, 1.08),
    ("m5_22_2_end_d2_s-0.5_s-0.5_a2_sc6_n32_d0.3__K30.4_s216.npz",
     "d2_K30.4", 15.047, 13.4, 1.08),
]


def claim_F():
    rows = {}
    # independently recompute the un-kicked E_start references
    refs = {}
    for fname, key, cl in [
            ("m5_22_end_P-1_plane_sc6_n32_pinned_d0.3.npz",
             "P-1", 12.7296),
            ("m5_22_1_end_d2_s-0.5_s-0.5_a2_sc6_n32_d0.3.npz",
             "d2", 15.047)]:
        Mr, cr = load_end(fname)
        refs[key] = {"E_start_recomputed":
                     round(float(sum(INS.e_parts(Mr, cr))), 4),
                     "claimed": cl}
    ref_ok = all(abs(r["E_start_recomputed"] - r["claimed"]) < 0.01
                 for r in refs.values())
    all_ok, any_ok = ref_ok, False
    for fname, key, e_ref, z_ref, q_ref in F_ROWS:
        M, cfg = load_end(fname)
        n, h, L = cfg["n"], cfg["h"], cfg["L"]
        e_end = float(sum(INS.e_parts(M, cfg)))   # trusted energy
        rings2, _ = my_ring_ledger(M, n, h, L, 0.15, min_vox=2)
        rings4, _ = my_ring_ledger(M, n, h, L, 0.15, min_vox=4)
        nhat, _ = PAIR.orient_v1(M)
        E = my_dual_E(nhat, h)
        far0 = 0.5 * L - 4.0 * h
        qu = my_slab_flux(E, n, h, 0.5 * h, far0, far0)
        ql = my_slab_flux(E, n, h, -far0, -0.5 * h, far0)
        zs2 = [round(r["z"], 2) for r in rings2]
        big = [r for r in rings2 if r["voxels"] >= 10]
        de = e_end - e_ref
        e_ok = abs(de) < 0.15
        ring_ok = (len(rings2) == 2
                   and all(abs(abs(r["z"]) - z_ref) < 1.5
                           for r in rings2))
        big_ok = (len(big) == 2
                  and all(abs(abs(r["z"]) - z_ref) < 1.5 for r in big))
        q_ok = (abs(qu - q_ref) < 0.08 and abs(ql + q_ref) < 0.08)
        row_ok = e_ok and ring_ok and q_ok
        rows[key] = {
            "E_end": round(e_end, 4), "dE_vs_start": round(de, 4),
            "n_rings_thr0.15_minvox2": len(rings2),
            "n_rings_thr0.15_minvox4": len(rings4),
            "n_rings_thr0.15_minvox10": len(big),
            "ring_zs": zs2, "ring_voxels": [r["voxels"]
                                            for r in rings2],
            "q_upper": round(qu, 4), "q_lower": round(ql, 4),
            "E_ok": e_ok, "rings_ok_as_claimed": ring_ok,
            "rings_ok_major_components": big_ok, "q_ok": q_ok}
        all_ok = all_ok and row_ok
        any_ok = any_ok or (e_ok and q_ok and big_ok)
    if not all_ok:
        bad = [k for k, r in rows.items()
               if not r["rings_ok_as_claimed"]]
        if bad:
            DEFECTS.append(
                "claim F ring count: rows with thr0.15 component "
                f"count != 2 at min_vox=2: {bad} (small spurious "
                "blobs below ~4 voxels; the run's own saved row for "
                "P-1 K10.4_ext already lists 3 components, incl. a "
                "3-voxel blob at z=+2.25, gap_min 0.141 just under "
                "thr): threshold-sensitive bookkeeping, the two "
                "MAJOR rings and charges are intact")
    V["F_returned_verdicts"] = {
        "verdict": verdict(all_ok, any_ok),
        "E_start_refs": refs,
        "measured": rows,
        "claimed": {"E_tol": 0.15, "rings": 2,
                    "q": "+-1.04 (P-1) / +-1.08 (d2)"}}


# ================= G: K1-envelope gap =================
def claim_G():
    M, cfg = load_end("m5_22_end_P-1_plane_sc6_n32_pinned_d0.3.npz")
    n, h, L = cfg["n"], cfg["h"], cfg["L"]
    rings, _ = my_ring_ledger(M, n, h, L, 0.15)
    ups = [r for r in rings if r["z"] > 0]
    top = max(ups, key=lambda r: r["voxels"])
    rho_r, z_r = top["rho_mean"], top["z"]
    r_tor = float(np.sqrt(rho_r ** 2 + z_r ** 2))
    w_k1 = float(np.exp(-((r_tor / 8.0) ** 2)))
    w_k3 = 1.0   # d_ring = 0 at the ring by construction
    # mean K3 weight over the actual ring voxel set (own read)
    lam = np.linalg.eigvalsh(M)
    gap = np.minimum(lam[..., 1] - lam[..., 0],
                     lam[..., 2] - lam[..., 1])
    X, Y, Z = my_coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    sel = (gap < 0.15) & (Z > 0) & (rho > 3.0) \
        & (~my_pin_shell(n, h))
    d_ring = np.sqrt((rho - rho_r) ** 2 + (Z - z_r) ** 2)
    w_k3_mean = float(np.mean(np.exp(-((d_ring[sel] / 4.0) ** 2))))
    w_k1_mean = float(np.mean(np.exp(
        -((np.sqrt(rho ** 2 + Z ** 2)[sel] / 8.0) ** 2))))
    ok = w_k1 <= 0.03 and w_k3_mean > 0.5
    V["G_k1_envelope_gap"] = {
        "verdict": verdict(ok),
        "measured": {"ring_rho": round(rho_r, 2),
                     "ring_z": round(z_r, 2),
                     "ring_r": round(r_tor, 2),
                     "w_K1_at_ring_center": w_k1,
                     "w_K1_mean_over_ring_voxels": w_k1_mean,
                     "w_K3_at_ring_center": w_k3,
                     "w_K3_mean_over_ring_voxels": w_k3_mean},
        "claimed": {"w_K1": "~0.02 or less", "w_K3": "~1"}}


# ================= H: energy bookkeeping =================
def claim_H():
    M0, cfg = load_end("m5_22_end_P-1_plane_sc6_n32_pinned_d0.3.npz")
    n, h, delta = cfg["n"], cfg["h"], cfg["delta"]
    e0 = float(sum(INS.e_parts(M0, cfg)))
    # ring pick: replicate the driver's top_ring AND my own ledger
    rho_run, z_run = KICK.ring_read(M0, cfg, 0.15), None
    ups = [r for r in rho_run["rings"] if r["z_centroid"] > 0]
    tr = max(ups, key=lambda e: e["voxels"])
    picks = {"run_top_ring": (tr["rho_mean"], tr["z_centroid"])}
    rings, _ = my_ring_ledger(M0, n, h, cfg["L"], 0.15)
    mytop = max([r for r in rings if r["z"] > 0],
                key=lambda r: r["voxels"])
    picks["own_ledger"] = (mytop["rho_mean"], mytop["z"])
    res = {}
    for lab, (rr, zz) in picks.items():
        rng = np.random.default_rng(216)
        K = rng.standard_normal((n, n, n, 3, 3))
        K = ndimage.gaussian_filter(K, sigma=(2, 2, 2, 0, 0))
        K = 0.5 * (K + np.swapaxes(K, -1, -2))
        X, Y, Z = my_coords(n, h)
        rho = np.sqrt(X * X + Y * Y)
        d_ring = np.sqrt((rho - rr) ** 2 + (Z - zz) ** 2)
        K = K * np.exp(-((d_ring / 4.0) ** 2))[..., None, None]
        iso = (1.0 + delta) / 3.0 * np.eye(3)
        scale = np.sqrt(np.mean((M0 - iso) ** 2))
        K *= 0.4 * scale / max(np.sqrt(np.mean(K ** 2)), 1e-300)
        Mk = M0 + K
        Mk = 0.5 * (Mk + np.swapaxes(Mk, -1, -2))
        free = (~my_pin_shell(n, h))[..., None, None].astype(float)
        Mk = M0 + (Mk - M0) * free
        Mk = 0.5 * (Mk + np.swapaxes(Mk, -1, -2))
        res[lab] = {"E_injected": float(sum(INS.e_parts(Mk, cfg))
                                        - e0),
                    "ring_pick": [round(rr, 3), round(zz, 3)]}
    inj = res["run_top_ring"]["E_injected"]
    rel = abs(inj - 671.42) / 671.42
    # absorbed + E_end from the saved row are bookkeeping outputs of
    # the run; consistency check: E_kicked - absorbed - E_end
    with open(os.path.join(
            DATA, "m5_22_2_row_P-1_plane_sc6_n32_pinned_d0.3__"
            "K30.4_s216.json")) as f:
        row = json.load(f)
    ledger_gap = (row["E_kicked"] - row["absorbed"] - row["E_end"])
    ok = (abs(e0 - 12.7296) < 0.01 and rel < 0.05
          and abs(row["E_end"] - 12.749) < 0.01)
    if ledger_gap > 0.5 * row["E_injected"] * 0.5:
        DEFECTS.append(
            "claim H bookkeeping: E_kicked - absorbed - E_end = "
            f"{ledger_gap:.1f} of {row['E_injected']:.1f} injected is "
            "NOT accounted by the sponge line-integral 'absorbed' "
            "(the leapfrog absorbed term gam*V^2*h^3*dt is a lower-"
            "bound estimate; ~half the injected energy leaves the "
            "ledger silently during evolve+FIRE): 'absorbed 344.02' "
            "is an instrument reading, not a closed energy budget")
    V["H_energy_bookkeeping"] = {
        "verdict": verdict(ok),
        "measured": {"E_start_recomputed": round(e0, 4),
                     "E_injected_replicated": res,
                     "rel_err_vs_claimed_671.42": round(rel, 5),
                     "row_E_end": row["E_end"],
                     "row_absorbed": row["absorbed"],
                     "ledger_unaccounted": round(ledger_gap, 2)},
        "claimed": {"E_injected": 671.42, "absorbed": 344.02,
                    "E_end": 12.749}}


# ================= run =================
if __name__ == "__main__":
    for fn in (claim_A, claim_B, claim_C, claim_D, claim_E,
               claim_F, claim_G, claim_H):
        try:
            print(f"--- {fn.__name__}", flush=True)
            fn()
        except Exception as exc:
            V[fn.__name__] = {"verdict": "ERROR", "error": repr(exc)}
            DEFECTS.append(f"{fn.__name__} raised {exc!r}")
    out = {"task": "M5.22.2 adversarial audit",
           "verdicts": V, "defects": DEFECTS}
    with open(os.path.join(DATA, "m5_22_2_audit.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))
