"""M5.32 R22-1s adversarial audit: an independent attempt to refute the claims made on the
slaved continuation rows (m5_32_r22_1s_slaved.py, data/m5_32_r22_1s_slaved.json).

EQUATIONS FIRST
---------------
Energy under audit: E[M] = certified quartic (production energy_parts) + V4,
    V4 = w h^3 sum_cells sum_{p=1..4} (tr N^p - C_p)^2,  N = M eta,  C_p = sum_roots q^p,
roots (-8, 1, 0.3, 0.3). The audited instrument solves M_00 per cell and descends on the six
spatial entries.

Own reads (none of them taken from the audited script):
  (a) own V4 per cell and own d V4 / d M_00, plus a brute-force scan of the per-cell polynomial
      in u = -M_00 to test that the slaved M_00 is the GLOBAL per-cell minimum;
  (b) own surface degree reader: the top eigenvector on the closed surface of a cell cube,
      oriented on the surface by a maximum spanning tree of |v . v'|, conflicting links counted
      over every surface link, global sign fixed OUTWARD, solid angle summed on two different
      triangulations of every quad;
  (c) own plaquette map: a lattice plaquette is pierced by a half-integer disclination of the
      top eigenvector when the product of the four link dot products around it is negative
      (sign-convention free); elementary cubes with no pierced face carry an integer degree;
  (d) rotation-invariant profiles: sorted eigenvalues per radial bin (mean and spread), the
      trace, the biaxiality (mid - small), E(<R), own V per unit volume;
  (e) direct field comparison modulo the 48 cubic symmetries (self-tested on the exact exterior).

Writes data/m5_32_r22_1s_audit.json and prints S<id>.<n> PASS|FAIL lines.
"""

import importlib.util
import itertools
import json
import os
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import breadth_first_order, minimum_spanning_tree

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
ROWS_JSON = os.path.join(DATA, "m5_32_r22_1s_slaved.json")
NPZ_DIR = os.path.join(DATA, "m5_32_r22_1s")
SRC_DIR = os.path.join(DATA, "m5_32_r22_1")
OUT_JSON = os.path.join(DATA, "m5_32_r22_1s_audit.json")
T0 = time.time()
G_TIME, DELTA = 8.0, 0.3
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
PIN25 = [f"{s}_pin_d0.3_w25_n32_L48_x9000_s" for s in ("rad", "bia", "perm", "obl")]
FREE25 = [f"{s}_free_d0.3_w25_n32_L48_s" for s in ("rad", "bia")]
PIN1 = [f"{s}_pin_d0.3_w1_n32_L48_s" for s in ("rad", "bia")]
FREE1 = [f"{s}_free_d0.3_w1_n32_L48_s" for s in ("rad", "bia")]
ALL = PIN25 + FREE25 + PIN1 + FREE1
CUBES = (("r6", 12, 19), ("r9", 10, 21), ("r12", 8, 23))
CUBES_WIDE = CUBES + (("r15", 6, 25), ("r18", 4, 27), ("r21", 2, 29))
CUBES_PROD = (("p6", 12, 20), ("p9", 10, 22), ("p12", 8, 24))  # the audited reader's placement
CHECKS = []


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def clean(v):
    if isinstance(v, dict):
        return {str(k): clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [clean(x) for x in v]
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    if isinstance(v, np.bool_):
        return bool(v)
    if isinstance(v, np.ndarray):
        return clean(v.tolist())
    return v


def check(cid, claim, method, value, expected, ok):
    rec = {
        "id": cid,
        "claim": claim,
        "method": method,
        "value": clean(value),
        "expected": expected,
        "verdict": "PASS" if ok else "FAIL",
    }
    CHECKS.append(rec)
    txt = json.dumps(rec["value"])
    print(f"{cid} {rec['verdict']} {claim} :: {txt[:600]}", flush=True)


# ================= own reads =================
def own_v4_density(M, roots, w):
    """own V4 per cell WITHOUT the volume factor: w sum_p (tr N^p - C_p)^2."""
    N = M @ ETA
    P = N.copy()
    v = np.zeros(M.shape[:3])
    for k in range(1, 5):
        if k > 1:
            P = P @ N
        v += (np.trace(P, axis1=-2, axis2=-1) - sum(q**k for q in roots)) ** 2
    return w * v


def own_dv_du(M, roots):
    """own d/du of sum_p (u^p + s_p - C_p)^2 at u = -M_00 (block-diagonal field), per cell."""
    S = M[..., 1:, 1:]
    u = -M[..., 0, 0]
    P = S.copy()
    g = np.zeros(M.shape[:3])
    for k in range(1, 5):
        if k > 1:
            P = P @ S
        res = u**k + np.trace(P, axis1=-2, axis2=-1) - sum(q**k for q in roots)
        g += 2.0 * res * k * u ** (k - 1)
    return g


def global_min_gap(M, roots, mask):
    """brute-force scan of the per-cell polynomial in u: value at the slaved u minus the scanned
    minimum over u in [-14, 6] (positive = the slaved root is NOT the global minimum)."""
    S = M[..., 1:, 1:][mask]
    u0 = -M[..., 0, 0][mask]
    s, P = [], S.copy()
    for k in range(1, 5):
        if k > 1:
            P = P @ S
        s.append(np.trace(P, axis1=-2, axis2=-1) - sum(q**k for q in roots))
    ug = np.linspace(-14.0, 6.0, 4001)
    worst = 0.0
    for a in range(0, len(u0), 2048):
        sl = slice(a, a + 2048)
        val = sum((ug[None, :] ** k + s[k - 1][sl, None]) ** 2 for k in range(1, 5))
        at = sum((u0[sl] ** k + s[k - 1][sl]) ** 2 for k in range(1, 5))
        worst = max(worst, float(np.max(at - val.min(axis=1))))
    return worst


def reslave(M, roots, mask, iters=40):
    """own Newton solve of the per-cell M_00 on the masked cells (warm start = current M_00)."""
    S = M[..., 1:, 1:][mask]
    s, P = [], S.copy()
    for k in range(1, 5):
        if k > 1:
            P = P @ S
        s.append(np.trace(P, axis1=-2, axis2=-1) - sum(q**k for q in roots))
    u = -M[..., 0, 0][mask]
    for _ in range(iters):
        g = sum(2.0 * (u**k + s[k - 1]) * k * u ** (k - 1) for k in range(1, 5))
        H = sum(
            2.0 * (k * u ** (k - 1)) ** 2
            + (2.0 * (u**k + s[k - 1]) * k * (k - 1) * u ** (k - 2) if k > 1 else 0.0)
            for k in range(1, 5)
        )
        u = u - g / H
    out = M.copy()
    m00 = out[..., 0, 0].copy()
    m00[mask] = -u
    out[..., 0, 0] = m00
    return out


def squeeze_split(M, mask, t):
    """move the two lower spatial eigenvalues toward their mean by the fraction t on the masked
    cells (t = 1: uniaxial about the top eigenvector; t < 0: more biaxial), eigenvectors kept."""
    lam, vec = np.linalg.eigh(M[..., 1:, 1:])
    mean = 0.5 * (lam[..., 0] + lam[..., 1])
    lam2 = lam.copy()
    lam2[..., 0] = lam[..., 0] + t * (mean - lam[..., 0])
    lam2[..., 1] = lam[..., 1] + t * (mean - lam[..., 1])
    S2 = np.einsum("...ak,...k,...bk->...ab", vec, lam2, vec)
    out = M.copy()
    blk = out[..., 1:, 1:].copy()
    blk[mask] = S2[mask]
    out[..., 1:, 1:] = blk
    return out


def solid_angle(p, q, r):
    num = np.einsum("ia,ia->i", p, np.cross(q, r))
    den = (
        1.0
        + np.einsum("ia,ia->i", p, q)
        + np.einsum("ia,ia->i", q, r)
        + np.einsum("ia,ia->i", r, p)
    )
    return 2.0 * np.arctan2(num, den)


def cube_degree(v, lo, hi, n=32):
    """own surface reader on the cell cube [lo, hi]^3 (both ends included)."""
    quads = []
    rng = np.arange(lo, hi)
    A, B = np.meshgrid(rng, rng, indexing="ij")
    A, B = A.ravel(), B.ravel()
    for ax in range(3):
        ua, va = (ax + 1) % 3, (ax + 2) % 3
        for side in (hi, lo):
            q = np.zeros((len(A), 4, 3), dtype=int)
            for c, (da, db) in enumerate(((0, 0), (1, 0), (1, 1), (0, 1))):
                q[:, c, ax] = side
                q[:, c, ua] = A + da
                q[:, c, va] = B + db
            if side == lo:
                q = q[:, ::-1]
            quads.append(q)
    quads = np.concatenate(quads)
    lin = (quads[..., 0] * n + quads[..., 1]) * n + quads[..., 2]
    uniq, inv = np.unique(lin, return_inverse=True)
    Q = inv.reshape(lin.shape)
    cells = np.stack(np.unravel_index(uniq, (n, n, n)), axis=-1)
    vs = v[cells[:, 0], cells[:, 1], cells[:, 2]]
    e = np.concatenate([Q[:, [c, (c + 1) % 4]] for c in range(4)])
    e = np.unique(np.sort(e, axis=1), axis=0)
    dots = np.einsum("ia,ia->i", vs[e[:, 0]], vs[e[:, 1]])
    m = len(uniq)
    Wt = coo_matrix((1.0 - np.abs(dots) + 1e-9, (e[:, 0], e[:, 1])), shape=(m, m)).tocsr()
    T = minimum_spanning_tree(Wt)
    order, pred = breadth_first_order(T + T.T, 0, directed=False, return_predecessors=True)
    sign = np.zeros(m)
    sign[order[0]] = 1.0
    for b in order[1:]:
        a = pred[b]
        sign[b] = sign[a] * (1.0 if np.dot(vs[a], vs[b]) >= 0.0 else -1.0)
    vo = vs * sign[:, None]
    conflicts = int(np.sum(np.einsum("ia,ia->i", vo[e[:, 0]], vo[e[:, 1]]) < 0.0))
    ctr = 0.5 * (lo + hi)
    if np.sum(vo * (cells - ctr)) < 0.0:
        vo = -vo
    degs = []
    for tri in (((0, 1, 2), (0, 2, 3)), ((0, 1, 3), (1, 2, 3))):
        om = 0.0
        for a, b, c in tri:
            om += float(np.sum(solid_angle(vo[Q[:, a]], vo[Q[:, b]], vo[Q[:, c]])))
        degs.append(om / (4.0 * np.pi))
    return {"deg_tri_a": degs[0], "deg_tri_b": degs[1], "conflicts": conflicts, "links": len(e)}


def plaquette_map(v, h, gap=None):
    """pierced plaquettes of the director v: centers (physical coordinates) of every lattice
    plaquette whose four link dot products multiply to a negative number."""
    n = v.shape[0]
    x = (np.arange(n) - (n - 1) / 2.0) * h
    out, diag, weak = [], [], []
    pierced_faces = {}
    for ax in range(3):
        ua, va = (ax + 1) % 3, (ax + 2) % 3

        def sh(du, dv, ua=ua, va=va):
            sl = [slice(None)] * 3
            sl[ua] = slice(du, n - 1 + du)
            sl[va] = slice(dv, n - 1 + dv)
            return v[tuple(sl)]

        c00, c10, c11, c01 = sh(0, 0), sh(1, 0), sh(1, 1), sh(0, 1)
        d = (
            np.einsum("...a,...a->...", c00, c10)
            * np.einsum("...a,...a->...", c10, c11)
            * np.einsum("...a,...a->...", c11, c01)
            * np.einsum("...a,...a->...", c01, c00)
        )
        P = d < 0.0
        pierced_faces[ax] = P
        idx = np.argwhere(P).astype(float)
        if len(idx) and gap is not None:

            def shg(du, dv, ua=ua, va=va):
                sl = [slice(None)] * 3
                sl[ua] = slice(du, n - 1 + du)
                sl[va] = slice(dv, n - 1 + dv)
                return gap[tuple(sl)]

            gmin = np.minimum.reduce([shg(0, 0), shg(1, 0), shg(1, 1), shg(0, 1)])
            diag.append(gmin[P])
            dmin = np.minimum.reduce(
                [
                    np.abs(np.einsum("...a,...a->...", c00, c10)),
                    np.abs(np.einsum("...a,...a->...", c10, c11)),
                    np.abs(np.einsum("...a,...a->...", c11, c01)),
                    np.abs(np.einsum("...a,...a->...", c01, c00)),
                ]
            )
            weak.append(dmin[P])
        if len(idx):
            pos = np.zeros_like(idx)
            for a in range(3):
                pos[:, a] = (idx[:, a] - (n - 1) / 2.0 + (0.5 if a in (ua, va) else 0.0)) * h
            out.append(pos)
    pos = np.concatenate(out) if out else np.zeros((0, 3))
    gaps = np.concatenate(diag) if diag else np.zeros(0)
    wk = np.concatenate(weak) if weak else np.zeros(0)
    return pos, pierced_faces, gaps, wk


def point_defects(v, pierced, h):
    """integer degree of every elementary cube none of whose six faces is pierced."""
    n = v.shape[0]
    m = n - 1
    bad = np.zeros((m, m, m), dtype=bool)
    for ax in range(3):
        P = pierced[ax]
        sl_lo, sl_hi = [slice(0, m)] * 3, [slice(0, m)] * 3
        sl_lo[ax], sl_hi[ax] = slice(0, m), slice(1, m + 1)
        bad |= P[tuple(sl_lo)] | P[tuple(sl_hi)]

    def c(dx, dy, dz):
        return v[dx : m + dx, dy : m + dy, dz : m + dz].reshape(-1, 3)

    def sg(a, b):
        return np.where(np.einsum("ia,ia->i", a, b) >= 0.0, 1.0, -1.0)

    V = {k: c(*k) for k in itertools.product((0, 1), repeat=3)}
    s = {(0, 0, 0): np.ones(m**3)}
    for k, par in (
        ((1, 0, 0), (0, 0, 0)),
        ((0, 1, 0), (0, 0, 0)),
        ((0, 0, 1), (0, 0, 0)),
        ((1, 1, 0), (1, 0, 0)),
        ((0, 1, 1), (0, 1, 0)),
        ((1, 0, 1), (0, 0, 1)),
        ((1, 1, 1), (1, 1, 0)),
    ):
        s[k] = s[par] * sg(V[par], V[k])
    Vo = {k: V[k] * s[k][:, None] for k in V}
    om = np.zeros(m**3)
    for ax in range(3):
        ua, va = (ax + 1) % 3, (ax + 2) % 3
        for side in (1, 0):
            q = []
            for da, db in ((0, 0), (1, 0), (1, 1), (0, 1)):
                k = [0, 0, 0]
                k[ax], k[ua], k[va] = side, da, db
                q.append(Vo[tuple(k)])
            if side == 0:
                q = q[::-1]
            om += solid_angle(q[0], q[1], q[2]) + solid_angle(q[0], q[2], q[3])
    deg = np.rint(om / (4.0 * np.pi)).reshape(m, m, m)
    deg[bad] = 0
    idx = np.argwhere(deg != 0)
    pos = (idx - (n - 1) / 2.0 + 0.5) * h
    return pos, int(bad.sum())


def link_stats(v):
    """|v . v'| over every nearest-neighbor link: fractions under cos 45 deg and under 0.2."""
    n = v.shape[0]
    d = []
    for ax in range(3):
        a = np.take(v, range(0, n - 1), axis=ax)
        b = np.take(v, range(1, n), axis=ax)
        d.append(np.abs(np.einsum("...a,...a->...", a, b)).ravel())
    d = np.concatenate(d)
    return {
        "frac_links_over_45_deg": float((d < np.sqrt(0.5)).mean()),
        "frac_links_abs_dot_under_0.2": float((d < 0.2).mean()),
        "min_abs_dot": float(d.min()),
    }


def cubic_group():
    out = []
    for perm in itertools.permutations(range(3)):
        for flips in itertools.product((1, -1), repeat=3):
            out.append((perm, flips))
    return out


def apply_sym(S, perm, flips):
    B = S.transpose(tuple(perm) + (3, 4))
    B = B[..., list(perm), :][..., :, list(perm)]
    f = np.array(flips, dtype=float)
    for a in range(3):
        if flips[a] < 0:
            B = np.flip(B, axis=a)
    return B * f[:, None] * f[None, :]


def radial_profiles(M, e_cell, v_cell, r, h):
    lam = np.linalg.eigvalsh(M[..., 1:, 1:])
    edges = np.arange(0.0, 24.0 + 1e-9, 3.0)
    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        sh = (r >= a) & (r < b)
        rows.append(
            {
                "r_lo": a,
                "r_hi": b,
                "cells": int(sh.sum()),
                "lam_mean": [float(lam[..., k][sh].mean()) for k in range(3)],
                "lam_std": [float(lam[..., k][sh].std()) for k in range(3)],
                "trace_mean": float(lam[sh].sum(axis=-1).mean()),
                "trace_std": float(lam[sh].sum(axis=-1).std()),
                "p2_mean": float((lam[sh] ** 2).sum(axis=-1).mean()),
                "biax_mean": float((lam[..., 1] - lam[..., 0])[sh].mean()),
                "biax_std": float((lam[..., 1] - lam[..., 0])[sh].std()),
                "v_per_volume": float(v_cell[sh].mean()),
                "e_per_volume": float(e_cell[sh].mean() / h**3),
            }
        )
    Rs = np.arange(3.0, 24.0 + 1e-9, 3.0)
    e_in = [float(e_cell[r < R].sum()) for R in Rs]
    return lam, rows, Rs.tolist(), e_in


def implied_split(lam):
    """the (mid - small) split that holds tr S = 1.6 and tr S^2 = 1.18 (the vacuum values of
    the first two power sums) at each cell's measured top eigenvalue, against the actual one."""
    top = lam[:, 2]
    arg = 2.0 * (1.18 - top**2) - (1.6 - top) ** 2
    imp = np.sqrt(np.maximum(arg, 0.0))
    act = lam[:, 1] - lam[:, 0]
    p2 = (lam**2).sum(axis=1)
    return {
        "p2_mean": float(p2.mean()),
        "p2_std": float(p2.std()),
        "actual_split_mean": float(act.mean()),
        "implied_split_mean": float(imp.mean()),
        "corr_actual_vs_implied": float(np.corrcoef(act, imp)[0, 1]),
        "rms_actual_minus_implied": float(np.sqrt(np.mean((act - imp) ** 2))),
    }


def r_half_of(e_cell, r):
    o = np.argsort(r.ravel())
    cum = np.cumsum(e_cell.ravel()[o])
    return float(r.ravel()[o][np.searchsorted(cum, 0.5 * cum[-1])])


# ================= main =================
def main():
    R21 = _load("m5_32_r21_1_runs", "m5_32_r21_1_runs.py")
    R20, B3, R0 = R21.R20, R21.B3, R21.R0
    CORES = _load("m5_32_r22_1_cores_for_audit", "m5_32_r22_1_cores.py")
    cfg = R21.cfg_of(32, 48.0, G_TIME, DELTA)
    p = R21.params_of(G_TIME, DELTA)
    roots = R0.roots_of(cfg, degenerate=True)
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = B3.coords(n, h)
    r = np.sqrt(X * X + Y * Y + Z * Z)
    cheb = np.maximum.reduce([np.abs(X), np.abs(Y), np.abs(Z)])
    rows = json.load(open(ROWS_JSON))["rows"]
    pin = B3.pin_shell(n, h)
    seed = CORES.seed_core(cfg, "rad", DELTA)
    rh = np.stack([X, Y, Z], axis=-1) / r[..., None]
    own_ext = np.zeros((n, n, n, 4, 4))
    own_ext[..., 1:, 1:] = DELTA * np.eye(3) + (1.0 - DELTA) * rh[..., :, None] * rh[..., None, :]
    own_ext[..., 0, 0] = G_TIME
    log(f"roots {list(roots)} n {n} h {h} pinned cells {int(pin.sum())} W1 {B3.W1}")

    F = {}
    for tag in ALL:
        w = B3.W1 * rows[tag]["w1s"]
        pot = ("v4", roots, w)
        M = np.load(os.path.join(NPZ_DIR, tag + ".npz"))["M"]
        parts = R20.energy_parts(M, cfg, p, pot)
        E, Gm, info = R20.energy_grad(M, cfg, p, pot)
        e_cell = R20.density(M, cfg, pot)
        v_cell = own_v4_density(M, roots, w)
        mask = R21.free_mask(cfg, rows[tag]["bnd"] == "pin")
        lam, vec = np.linalg.eigh(M[..., 1:, 1:])
        F[tag] = dict(
            M=M,
            w=w,
            pot=pot,
            parts=parts,
            E=float(E),
            G=Gm,
            e=e_cell,
            v=v_cell,
            mask=mask,
            lam=lam,
            top=vec[..., :, 2],
            row=rows[tag],
        )
        log(f"loaded {tag} E {parts['E_total']:.6f}")

    # ---------------- S1 ----------------
    val, ok = {}, True
    for tag in ALL:
        f = F[tag]
        d = abs(f["parts"]["E_total"] - f["row"]["E"])
        d2 = abs(f["parts"]["E_total"] - f["row"]["chunks"][-1]["E"])
        dv = abs(float(f["v"].sum()) * h**3 - f["parts"]["V"])
        val[tag] = {"E": f["parts"]["E_total"], "abs_diff_json": d, "own_V_minus_prod_V": dv}
        ok &= d < 1e-10 and d2 < 1e-10 and dv < 1e-9 * max(1.0, f["parts"]["V"])
    check(
        "S1.1",
        "every end field reproduces its JSON energy; own V4 equals production V",
        "production energy_parts on the npz; own per-cell V4 summed",
        val,
        "abs diff under 1e-10",
        ok,
    )

    val = {tag: float(np.abs(F[tag]["M"][..., 0, 1:]).max()) for tag in ALL}
    val_sym = {tag: float(np.abs(F[tag]["M"] - F[tag]["M"].swapaxes(-1, -2)).max()) for tag in ALL}
    check(
        "S1.2",
        "M_0i = 0 and M symmetric on every end field",
        "max abs",
        {"M0i": val, "asym": val_sym},
        "M_0i exactly 0; asymmetry under 1e-15",
        all(x == 0.0 for x in val.values()) and all(x < 1e-15 for x in val_sym.values()),
    )

    val, ok = {}, True
    for tag in PIN25 + PIN1:
        M = F[tag]["M"]
        bit = bool(np.array_equal(M[pin], seed[pin]))
        own_seed = CORES.seed_core(cfg, tag.split("_")[0], DELTA)
        val[tag] = {
            "bitwise_equal_seed_core_rad": bit,
            "bitwise_equal_seed_core_of_own_kind": bool(np.array_equal(M[pin], own_seed[pin])),
            "max_abs_vs_seed_core_rad": float(np.abs(M[pin] - seed[pin]).max()),
            "max_abs_vs_own_formula": float(np.abs(M[pin] - own_ext[pin]).max()),
            "M00_on_shell": [float(M[pin][:, 0, 0].min()), float(M[pin][:, 0, 0].max())],
        }
        ok &= bit and val[tag]["max_abs_vs_own_formula"] < 1e-12
    check(
        "S1.3",
        "pinned shell is bitwise the exact exterior, M_00 = 8",
        "array_equal against seed_core(cfg,'rad',0.3) on pin_shell; own formula beside it",
        val,
        "bitwise equal; own formula within 1e-14",
        ok,
    )

    val, ok = {}, True
    for tag in ALL:
        f = F[tag]
        Gf = f["G"][f["mask"]]
        own00 = f["w"] * h**3 * own_dv_du(f["M"], roots)[f["mask"]]
        gap = global_min_gap(f["M"], roots, f["mask"])
        val[tag] = {
            "max_G00_prod": float(np.abs(Gf[:, 0, 0]).max()),
            "max_G00_own": float(np.abs(own00).max()),
            "max_G0i": float(np.abs(Gf[:, 0, 1:]).max()),
            "max_G_spatial": float(np.abs(Gf[:, 1:, 1:]).max()),
            "max_G_all10": float(np.abs(Gf).max()),
            "json_fmax_spatial": f["row"]["chunks"][-1]["fmax_spatial"],
            "slaved_minus_global_min_per_cell": gap,
        }
        ok &= (
            val[tag]["max_G00_prod"] < 1e-8
            and val[tag]["max_G00_own"] < 1e-8
            and val[tag]["max_G0i"] < 1e-12
            and abs(val[tag]["max_G_spatial"] - val[tag]["json_fmax_spatial"]) < 1e-9
            and gap < 1e-9
        )
    check(
        "S1.4",
        "slaving residual negligible, slaved M_00 is the global per-cell minimum, no "
        "hidden gradient component",
        "production energy_grad all 10 entries on free cells; own "
        "dV/dM_00; brute-force scan of the per-cell polynomial over u in [-14, 6]",
        val,
        "max|G_00| < 1e-8, max|G_0i| < 1e-12, max|G| = the JSON spatial value, gap < 1e-9",
        ok,
    )

    # ---------------- own degree reads on every field ----------------
    DEG, PLQ = {}, {}
    for tag in ALL:
        f = F[tag]
        DEG[tag] = {nm: cube_degree(f["top"], lo, hi) for nm, lo, hi in CUBES_WIDE + CUBES_PROD}
        gap_tm = f["lam"][..., 2] - f["lam"][..., 1]
        pos, pierced, pg, wk = plaquette_map(f["top"], h, gap_tm)
        pts, nbad = point_defects(f["top"], pierced, h)
        rec = {
            "pierced_plaquettes": len(pos),
            "cubes_with_pierced_face": nbad,
            "links": link_stats(f["top"]),
        }
        if len(wk):
            rec["frac_pierced_with_weakest_link_abs_dot_under_0.2"] = float((wk < 0.2).mean())
        if len(pos):
            rr = np.linalg.norm(pos, axis=1)
            ch = np.abs(pos).max(axis=1)
            ctr = pos.mean(axis=0)
            evals, evecs = np.linalg.eigh(np.cov((pos - ctr).T))
            ev = np.sqrt(np.maximum(evals, 0.0))
            rec.update(
                normal_of_smallest_extent=evecs[:, 0].tolist(),
                corner_gap_top_mid_at_pierced_min_median_max=[
                    float(pg.min()),
                    float(np.median(pg)),
                    float(pg.max()),
                ],
                r_min=float(rr.min()),
                r_mean=float(rr.mean()),
                r_max=float(rr.max()),
                cheb_max=float(ch.max()),
                cheb_min=float(ch.min()),
                n_inside_cheb_12=int((ch < 12.0).sum()),
                centroid=ctr.tolist(),
                extent_rms_sorted=ev.tolist(),
            )
        rec["point_defects"] = [
            {"pos": q.tolist(), "r": float(np.linalg.norm(q))} for q in pts[:12]
        ]
        rec["n_point_defects"] = len(pts)
        gap = f["lam"][..., 2] - f["lam"][..., 1]
        rec["min_top_mid_gap"] = float(gap.min())
        rec["min_top_mid_gap_r"] = float(r.ravel()[np.argmin(gap)])
        PLQ[tag] = rec
        log(
            f"degree {tag} "
            + " ".join(f"{k}:{v['deg_tri_a']:+.2f}/{v['conflicts']}" for k, v in DEG[tag].items())
            + f" pierced {rec['pierced_plaquettes']} points {rec['n_point_defects']}"
        )

    # ---------------- S2 ----------------
    val, ok = {}, True
    claimed = {
        FREE25[0]: (5.74, 0.20),
        FREE25[1]: (5.19, 0.44),
        FREE1[0]: (2.22, 0.10),
        FREE1[1]: (2.22, 0.06),
    }
    for tag in FREE25 + FREE1:
        ch = F[tag]["row"]["chunks"]
        src = np.load(os.path.join(SRC_DIR, F[tag]["row"]["src"] + ".npz"))["M"]
        e0 = R20.energy_parts(src, cfg, p, F[tag]["pot"])["E_total"]
        val[tag] = {
            "E_start_own": e0,
            "E_end_own": F[tag]["parts"]["E_total"],
            "last_drop": ch[-1]["drop"],
            "claimed": claimed[tag],
        }
        ok &= abs(e0 - claimed[tag][0]) < 0.006 and ch[-1]["drop"] > 1e-3
        ok &= abs(F[tag]["parts"]["E_total"] - claimed[tag][1]) < 0.006
    check(
        "S2.1",
        "free rows: start and end energies as claimed, still falling",
        "production energy_parts on source and end fields",
        val,
        "within 0.006; drop > 1e-3",
        ok,
    )

    val, ok = {}, True
    for tag in FREE25 + FREE1:
        d = DEG[tag]
        val[tag] = {
            k: [round(v["deg_tri_a"], 6), round(v["deg_tri_b"], 6), v["conflicts"]]
            for k, v in d.items()
        }
        for k in ("r6", "r9"):
            ok &= d[k]["conflicts"] == 0 and abs(d[k]["deg_tri_a"]) < 1e-6
    check(
        "S2.2",
        "free rows: top-eigenvector degree is 0 on the r 6 and r 9 cubes",
        "own surface reader, symmetric cubes, [deg tri a, deg tri b, conflicting links]",
        val,
        "0 conflicts and degree 0 on r6 and r9 in all four free rows",
        ok,
    )

    val = {}
    for tag in FREE25 + FREE1:
        hist = []
        for c in F[tag]["row"]["chunks"]:
            s = ""
            for k in ("6", "9", "12"):
                dg, cf = c["degree_top"][k]
                s += "N" if cf else ("1" if abs(abs(dg) - 1) < 0.1 else "0")
            hist.append([c["iters"], s])
        val[tag] = hist
    check(
        "S2.3",
        "degree history (JSON chunks): per chunk the r6 r9 r12 state, 1 = charge, "
        "N = non-orientable, 0 = no charge",
        "read of the audited chunk records",
        val,
        "informational: inside-out order (r6 opens first, r6 clears first) = a ring leaving",
        True,
    )

    val, ok = {}, True
    for tag in FREE25 + FREE1:
        f = F[tag]
        Et = float(f["e"].sum())
        inner = cheb < 12.0
        top = f["top"][inner]
        Qm = np.einsum("ia,ib->ab", top, top) / len(top)
        axis = np.linalg.eigh(Qm)[1][:, -1]
        ang = np.degrees(np.arccos(np.clip(np.abs(top @ axis), 0.0, 1.0)))
        lam_in = f["lam"][inner]
        val[tag] = {
            "angle_to_mean_axis_deg_median_p90": [
                float(np.median(ang)),
                float(np.percentile(ang, 90)),
            ],
            "frac_cells_more_than_30_deg_off_axis": float(np.mean(ang > 30.0)),
            "frac_E_within_3_cells_of_a_face": float(f["e"][cheb > 24.0 - 4.5].sum() / Et),
            "frac_E_within_6_cells_of_a_face": float(f["e"][cheb > 24.0 - 9.0].sum() / Et),
            "frac_E_cheb_lt_12": float(f["e"][inner].sum() / Et),
            "director_order_cheb_lt_12": float(np.linalg.eigvalsh(Qm)[-1]),
            "lam_mean_cheb_lt_12": lam_in.mean(axis=0).tolist(),
            "lam_std_cheb_lt_12": lam_in.std(axis=0).tolist(),
            "plaquettes": PLQ[tag],
        }
        ok &= val[tag]["director_order_cheb_lt_12"] > 0.95
    check(
        "S2.4",
        "free rows: interior nearly uniform; remaining energy and winding near faces",
        "energy fractions by distance to the nearest face; top eigenvalue of <n n^T> over the "
        "inner cube (1 = uniform); pierced-plaquette map",
        val,
        "director order above 0.95 in the inner cube on all four",
        ok,
    )

    # ---------------- S3 ----------------
    val, ok = {}, True
    claimedE = {PIN25[0]: 3.540, PIN25[1]: 3.430, PIN25[2]: 3.495, PIN25[3]: 3.454}
    for tag in PIN25:
        f = F[tag]
        Et = float(f["e"].sum())
        val[tag] = {
            "E": f["parts"]["E_total"],
            "density_sum": Et,
            "last_drop": f["row"]["chunks"][-1]["drop"],
            "iters": f["row"]["iters"],
            "r_half": r_half_of(f["e"], r),
            "E_inside_r12": float(f["e"][r < 12.0].sum()),
            "frac_inside_r12": float(f["e"][r < 12.0].sum() / Et),
        }
        ok &= abs(val[tag]["E"] - claimedE[tag]) < 6e-4
        ok &= 15.3 < val[tag]["r_half"] < 16.7 and 0.735 < val[tag]["E_inside_r12"] < 0.925
        ok &= 2.5e-4 < val[tag]["last_drop"] < 0.0115
    Es = [val[t]["E"] for t in PIN25]
    val["spread_percent_of_min"] = 100.0 * (max(Es) - min(Es)) / min(Es)
    val["seed_rad_reference"] = {
        "r_half": r_half_of(R20.density(seed, cfg, F[PIN25[0]]["pot"]), r),
        "E": R20.energy_parts(seed, cfg, p, F[PIN25[0]]["pot"])["E_total"],
    }
    ok &= val["spread_percent_of_min"] < 3.25
    check(
        "S3.1",
        "pinned W1 x 25: energies, drops, r_half 15.4 to 16.6, E inside r 12 0.74 to "
        "0.92 (ABSOLUTE energy, see frac_inside_r12), spread within 3.2 percent",
        "production energy_parts and density; own radial sort",
        val,
        "as claimed",
        ok,
    )

    val, ok = {}, True
    for tag in PIN25:
        d = DEG[tag]
        val[tag] = {
            k: [round(v["deg_tri_a"], 6), round(v["deg_tri_b"], 6), v["conflicts"]]
            for k, v in d.items()
        }
        val[tag]["plaquettes"] = PLQ[tag]
        ok &= d["r12"]["conflicts"] == 0 and abs(d["r12"]["deg_tri_a"] - 1.0) < 1e-6
        ok &= abs(d["r12"]["deg_tri_b"] - 1.0) < 1e-6
        ok &= d["p12"]["conflicts"] == 0 and abs(d["p12"]["deg_tri_a"] - 1.0) < 1e-6
    check(
        "S3.2",
        "pinned W1 x 25: degree exactly 1 on the r 12 cube on all four",
        "own surface reader (outward sign, two triangulations), centered cubes r6 to r21 and "
        "the audited placement p6 p9 p12; plaquette map of the disclination loop",
        val,
        "r12 and p12: 0 conflicts, degree +1 on both triangulations",
        ok,
    )
    val, ok = {}, True
    for tag in PIN25:
        d = DEG[tag]
        val[tag] = {
            k: [d[k]["conflicts"], round(d[k]["deg_tri_a"], 6)] for k in ("r6", "r9", "p6", "p9")
        }
        val[tag]["json_last_chunk"] = F[tag]["row"]["chunks"][-1]["degree_top"]
        ok &= all(d[k]["conflicts"] > 0 for k in ("r6", "r9", "p6", "p9"))
    check(
        "S3.2b",
        "pinned W1 x 25: top eigenvector non-orientable on the r 6 and r 9 cubes on " "all four",
        "own reader: [conflicting links, solid-angle degree] on the centered cubes "
        "(r) and on the audited placement (p)",
        val,
        "conflicts > 0 on r6, r9, p6, p9 of all four",
        ok,
    )

    def same_object(tags, sid):
        prof = {}
        for tag in tags:
            f = F[tag]
            _, prows, Rs, e_in = radial_profiles(f["M"], f["e"], f["v"], r, h)
            prof[tag] = (prows, Rs, e_in, float(f["e"].sum()))
        lam_d = max(
            max(
                abs(prof[a][0][i]["lam_mean"][k] - prof[b][0][i]["lam_mean"][k])
                for i in range(len(prof[a][0]))
                for k in range(3)
            )
            for a in tags
            for b in tags
        )
        lam_where = max(
            (
                (
                    abs(prof[a][0][i]["lam_mean"][k] - prof[b][0][i]["lam_mean"][k]),
                    prof[a][0][i]["r_lo"],
                    k,
                )
                for a in tags
                for b in tags
                for i in range(len(prof[a][0]))
                for k in range(3)
            )
        )
        biax_d = max(
            max(
                abs(prof[a][0][i]["biax_mean"] - prof[b][0][i]["biax_mean"])
                for i in range(len(prof[a][0]))
            )
            for a in tags
            for b in tags
        )
        ein_d = max(
            max(abs(prof[a][2][i] - prof[b][2][i]) for i in range(len(prof[a][2])))
            for a in tags
            for b in tags
        )
        einf_d = max(
            max(
                abs(prof[a][2][i] / prof[a][3] - prof[b][2][i] / prof[b][3])
                for i in range(len(prof[a][2]))
            )
            for a in tags
            for b in tags
        )
        outer = [i for i in range(len(prof[tags[0]][0])) if prof[tags[0]][0][i]["r_lo"] >= 6.0]
        lam_d_outer = max(
            abs(prof[a][0][i]["lam_mean"][k] - prof[b][0][i]["lam_mean"][k])
            for a in tags
            for b in tags
            for i in outer
            for k in range(3)
        )
        biax_d_outer = max(
            abs(prof[a][0][i]["biax_mean"] - prof[b][0][i]["biax_mean"])
            for a in tags
            for b in tags
            for i in outer
        )
        group = cubic_group()
        S_ext = own_ext[..., 1:, 1:]
        self_test = max(float(np.abs(apply_sym(S_ext, *g) - S_ext).max()) for g in group)
        direct = {}
        for ia, a in enumerate(tags):
            for b in tags[ia + 1 :]:
                Sa, Sb = F[a]["M"][..., 1:, 1:], F[b]["M"][..., 1:, 1:]
                best = min(
                    float(np.sqrt(np.mean(np.sum((Sa - apply_sym(Sb, *g)) ** 2, axis=(-1, -2)))))
                    for g in group
                )
                ident = float(np.sqrt(np.mean(np.sum((Sa - Sb) ** 2, axis=(-1, -2)))))
                direct[a[:3] + "_vs_" + b[:3]] = {"identity": ident, "best_of_48": best}
        scale = {
            t[:3]: float(
                np.sqrt(np.mean(np.sum((F[t]["M"][..., 1:, 1:] - S_ext) ** 2, axis=(-1, -2))))
            )
            for t in tags
        }
        return {
            "max_diff_shell_mean_eigenvalue": lam_d,
            "where_[diff, r_lo, eig_index]": list(lam_where),
            "max_diff_shell_mean_biaxiality": biax_d,
            "max_diff_shell_mean_eigenvalue_r_ge_6": lam_d_outer,
            "max_diff_shell_mean_biaxiality_r_ge_6": biax_d_outer,
            "max_diff_E_inside_R_absolute": ein_d,
            "max_diff_E_inside_R_fraction": einf_d,
            "cubic_group_self_test_on_exterior": self_test,
            "rms_frobenius_distance_per_cell": direct,
            "rms_distance_of_each_field_from_the_exterior": scale,
            "disclination_map": {
                t[:3]: {
                    k: PLQ[t].get(k)
                    for k in (
                        "pierced_plaquettes",
                        "r_min",
                        "r_mean",
                        "r_max",
                        "extent_rms_sorted",
                        "centroid",
                        "normal_of_smallest_extent",
                        "n_point_defects",
                        "corner_gap_top_mid_at_pierced_min_median_max",
                    )
                }
                for t in tags
            },
            "profiles": {t[:3]: prof[t][0] for t in tags},
            "E_inside_R": {
                t[:3]: dict(zip([str(x) for x in prof[t][1]], prof[t][2])) for t in tags
            },
        }, sid

    so, _ = same_object(PIN25, "S3.3")
    ratio = max(v["best_of_48"] for v in so["rms_frobenius_distance_per_cell"].values()) / min(
        so["rms_distance_of_each_field_from_the_exterior"].values()
    )
    so["worst_pair_distance_over_smallest_distance_from_exterior"] = ratio
    check(
        "S3.3",
        "pinned W1 x 25: are the four end fields one object (rotation-invariant reads "
        "and direct comparison modulo the 48 cubic symmetries)",
        "radial-bin profiles (3.0 wide), E(<R), disclination map, best-of-48 rms distance",
        so,
        "ONE object if invariant reads agree within 0.05 in eigenvalue and 0.05 in E(<R) "
        "fraction and the best-of-48 distance is under 0.3 of the distance from the exterior",
        so["max_diff_shell_mean_eigenvalue"] < 0.05
        and so["max_diff_E_inside_R_fraction"] < 0.05
        and ratio < 0.3,
    )

    # ---------------- S4 ----------------
    val, ok = {}, True
    claimed4 = {PIN1[0]: 2.0167, PIN1[1]: 2.0157}
    for tag in PIN1:
        f = F[tag]
        Gf = f["G"][f["mask"]]
        val[tag] = {
            "E": f["parts"]["E_total"],
            "last_drop": f["row"]["chunks"][-1]["drop"],
            "max_G_spatial": float(np.abs(Gf[:, 1:, 1:]).max()),
            "r_half": r_half_of(f["e"], r),
            "own_degree_[a, b, conflicts]": {
                k: [round(v["deg_tri_a"], 6), round(v["deg_tri_b"], 6), v["conflicts"]]
                for k, v in DEG[tag].items()
            },
            "json_degree_last_chunk": f["row"]["chunks"][-1]["degree_top"],
            "plaquettes": PLQ[tag],
        }
        ok &= abs(val[tag]["E"] - claimed4[tag]) < 6e-5
        ok &= 1.5e-4 < val[tag]["last_drop"] < 5.5e-4 and 19.0 < val[tag]["r_half"] < 19.6
        ok &= 2.3e-4 < val[tag]["max_G_spatial"] < 1.0e-3
        for k in ("r6", "r9", "r12"):
            ok &= DEG[tag][k]["conflicts"] == 0 and abs(DEG[tag][k]["deg_tri_a"] - 1.0) < 1e-6
    Es = [val[t]["E"] for t in PIN1]
    val["apart_percent"] = 100.0 * abs(Es[0] - Es[1]) / min(Es)
    ok &= val["apart_percent"] < 0.055
    check(
        "S4.1",
        "pinned W1: energies, drops, gradient, r_half 19.3; with the OUTWARD sign "
        "convention the degree is +1 on all three cubes of both rows with 0 conflicts (so the "
        "JSON sign flips are the reader's arbitrary root sign, not physics)",
        "production energy; own reader with the global sign fixed outward",
        val,
        "as claimed",
        ok,
    )

    so4, _ = same_object(PIN1, "S4.2")
    ratio4 = max(v["best_of_48"] for v in so4["rms_frobenius_distance_per_cell"].values()) / min(
        so4["rms_distance_of_each_field_from_the_exterior"].values()
    )
    so4["worst_pair_distance_over_smallest_distance_from_exterior"] = ratio4
    check(
        "S4.2",
        "pinned W1: rad and bia are one object",
        "same reads as S3.3",
        so4,
        "ONE object by the same criteria as S3.3 (0.05, 0.05, ratio under 0.3)",
        so4["max_diff_shell_mean_eigenvalue"] < 0.05
        and so4["max_diff_E_inside_R_fraction"] < 0.05
        and ratio4 < 0.3,
    )

    # ---------------- S5 ----------------
    val, ok = {}, True
    for tag in PIN25 + PIN1 + FREE25 + FREE1:
        f = F[tag]
        band = (r >= 6.0) & (r <= 18.0)
        lam = f["lam"][band]
        Eb = float(f["e"][band].sum())
        Vb = float(f["v"][band].sum()) * h**3
        val[tag] = {
            "lam_mean_r6_18": lam.mean(axis=0).tolist(),
            "lam_std_r6_18": lam.std(axis=0).tolist(),
            "lam_min_max_r6_18": [lam.min(axis=0).tolist(), lam.max(axis=0).tolist()],
            "trace_mean_std": [float(lam.sum(axis=1).mean()), float(lam.sum(axis=1).std())],
            "p2_mean": float((lam**2).sum(axis=1).mean()),
            "split_mid_minus_small_mean_std": [
                float((lam[:, 1] - lam[:, 0]).mean()),
                float((lam[:, 1] - lam[:, 0]).std()),
            ],
            "frac_cells_split_below_0.1": float(np.mean((lam[:, 1] - lam[:, 0]) < 0.1)),
            "V_per_volume_r6_18": float(f["v"][band].mean()),
            "V_per_volume_over_w": float(f["v"][band].mean() / f["w"]),
            "V_band_over_E_band": Vb / Eb,
            "V_total": f["parts"]["V"],
            "E_curv_total": f["parts"]["E_curv"],
            "V_band_over_V_total": Vb / f["parts"]["V"],
        }
    for tag in PIN25:
        m = val[tag]["lam_mean_r6_18"]
        ok &= abs(m[0] - 0.17) < 0.03 and abs(m[1] - 0.45) < 0.03 and abs(m[2] - 0.98) < 0.03
    for tag in PIN1:
        m = val[tag]["lam_mean_r6_18"]
        ok &= abs(m[0] - 0.20) < 0.03 and abs(m[1] - 0.47) < 0.03 and abs(m[2] - 0.96) < 0.03
    check(
        "S5.1",
        "pinned rows: interior spectrum about (0.17, 0.45, 0.98) at W1 x 25 and "
        "(0.20, 0.47, 0.96) at W1 on r 6 to 18, trace near 1.6 (free rows listed beside them)",
        "own eigenvalues over all cells with 6 <= r <= 18: mean, spread, extremes; own V4 per "
        "unit volume",
        val,
        "means within 0.03 of the claimed triple",
        ok,
    )

    tr = {t: val[t]["trace_mean_std"] for t in PIN25 + PIN1}
    prof_tr = {
        t: [
            [pr["r_lo"], round(pr["trace_mean"], 4), round(pr["trace_std"], 4)]
            for pr in radial_profiles(F[t]["M"], F[t]["e"], F[t]["v"], r, h)[1]
        ]
        for t in PIN1 + PIN25[:1]
    }
    check(
        "S5.1b",
        "pinned rows: the trace is kept near 1.6 on r 6 to 18",
        "own trace of the spatial block: [mean, std] on the band, and per 3.0 wide shell "
        "[r_lo, mean, std] for the two W1 rows and one W1 x 25 row",
        {"band": tr, "per_shell": prof_tr},
        "band mean within 0.01 of 1.6 and std under 0.01 on all six",
        all(abs(v[0] - 1.6) < 0.01 and v[1] < 0.01 for v in tr.values()),
    )

    val3, ok3 = {}, True
    for tag in PIN25 + PIN1:
        f = F[tag]
        base = f["parts"]
        rec = {}
        for t in (-0.25, 0.25, 0.5, 1.0):
            Mt = reslave(squeeze_split(f["M"], f["mask"], t), roots, f["mask"])
            pt = R20.energy_parts(Mt, cfg, p, f["pot"])
            rec[f"t={t:g}"] = {
                "dE": pt["E_total"] - base["E_total"],
                "dE_curv": pt["E_curv"] - base["E_curv"],
                "dV": pt["V"] - base["V"],
            }
        val3[tag] = rec
        ok3 &= rec["t=0.25"]["dE"] > 0.0 and rec["t=-0.25"]["dE"] > 0.0
    check(
        "S5.3",
        "is the biaxial split held by the energy (E rises when the two lower "
        "eigenvalues are moved toward OR away from each other, eigenvectors kept, M_00 "
        "re-slaved), or is it a residual the descent has not removed yet",
        "own squeeze of the split by the fraction t on every free cell, own Newton re-slave, "
        "production energy_parts; t = 1 is the uniaxial projection",
        val3,
        "dE > 0 at t = +0.25 and at t = -0.25 on all six pinned rows",
        ok3,
    )

    aniso = {}
    for tag in PIN25 + PIN1:
        f = F[tag]
        band = (r >= 6.0) & (r <= 18.0)
        sp = (f["lam"][..., 1] - f["lam"][..., 0])[band]
        rhb = rh[band]
        Tm = np.einsum("i,ia,ib->ab", sp, rhb, rhb) / sp.sum()
        Tev, Tvec = np.linalg.eigh(Tm)
        cosn = np.abs(rhb @ Tvec[:, 0])
        aniso[tag] = {
            "eig_of_split_weighted_rhat_rhat (isotropic = 1/3 each)": Tev.tolist(),
            "axis_of_least_split": Tvec[:, 0].tolist(),
            "split_mean_by_abs_cos_to_that_axis_[0-1/3, 1/3-2/3, 2/3-1]": [
                float(sp[(cosn >= a) & (cosn < b + 1e-12)].mean())
                for a, b in ((0.0, 1 / 3), (1 / 3, 2 / 3), (2 / 3, 1.0))
            ],
        }
    check(
        "S5.4",
        "is the interior split uniform over the sphere (the claim says nearly uniform)",
        "angular distribution of (mid - small) on r 6 to 18",
        aniso,
        "informational: PASS if the three angular bins agree within 25 percent",
        all(
            max(v[k]) / min(v[k]) < 1.25
            for v in aniso.values()
            for k in v
            if k.startswith("split_mean_by")
        ),
    )

    s25 = np.mean([val[t]["split_mid_minus_small_mean_std"][0] for t in PIN25])
    s1 = np.mean([val[t]["split_mid_minus_small_mean_std"][0] for t in PIN1])
    v25 = np.mean([val[t]["V_per_volume_over_w"] for t in PIN25])
    v1 = np.mean([val[t]["V_per_volume_over_w"] for t in PIN1])
    sf = np.mean([val[t]["split_mid_minus_small_mean_std"][0] for t in FREE25 + FREE1])
    val2 = {
        "split_mean_W1x25": float(s25),
        "split_mean_W1": float(s1),
        "split_ratio_25_over_1": float(s25 / s1),
        "V_per_volume_over_w_W1x25": float(v25),
        "V_per_volume_over_w_W1": float(v1),
        "V_per_volume_W1x25": float(v25 * 25 * B3.W1),
        "V_per_volume_W1": float(v1 * B3.W1),
        "split_mean_free_rows": float(sf),
        "p2_band_mean_std_and_split_implied_by_vacuum_p1_p2": {
            t: implied_split(F[t]["lam"][(r >= 6.0) & (r <= 18.0)]) for t in PIN25 + PIN1
        },
        "per_shell_split": {
            t[:3]
            + ("_w25" if "w25" in t else "_w1"): [
                [
                    pr["r_lo"],
                    round(pr["biax_mean"], 4),
                    round(pr["biax_std"], 4),
                    round(pr["trace_mean"], 4),
                    round(pr["v_per_volume"] / F[t]["w"], 6),
                ]
                for pr in radial_profiles(F[t]["M"], F[t]["e"], F[t]["v"], r, h)[1]
            ]
            for t in PIN25 + PIN1
        },
    }
    check(
        "S5.2",
        "the split of the degenerate pair is about the same at w and 25 w",
        "mean (mid - small) on r 6 to 18 per scale; V per unit volume divided by w (the pure "
        "spectral misfit) per scale; per shell [r_lo, split mean, split std, trace, V/w]",
        val2,
        "split ratio within 10 percent of 1",
        abs(s25 / s1 - 1.0) < 0.10,
    )

    # ---------------- S6 ----------------
    free_lines = {
        t: (PLQ[t]["pierced_plaquettes"], PLQ[t].get("n_inside_cheb_12", 0))
        for t in FREE25 + FREE1
    }
    val = {
        "free_rows_[pierced plaquettes, of which inside cheb 12]": free_lines,
        "free_rows_E_end": {t: F[t]["parts"]["E_total"] for t in FREE25 + FREE1},
        "free_rows_last_drop": {t: F[t]["row"]["chunks"][-1]["drop"] for t in FREE25 + FREE1},
        "pinned_w25_spread_percent": 100.0
        * (
            max(F[t]["parts"]["E_total"] for t in PIN25)
            / min(F[t]["parts"]["E_total"] for t in PIN25)
            - 1
        ),
        "pinned_w25_last_drops": {t: F[t]["row"]["chunks"][-1]["drop"] for t in PIN25},
        "pinned_w25_same_object_reads": {
            k: so[k]
            for k in (
                "max_diff_shell_mean_eigenvalue",
                "max_diff_shell_mean_biaxiality",
                "max_diff_E_inside_R_fraction",
                "rms_frobenius_distance_per_cell",
                "rms_distance_of_each_field_from_the_exterior",
            )
        },
    }
    ok = all(v[1] == 0 for v in free_lines.values())
    check(
        "S6.1",
        "free rows: no winding (no pierced plaquette, no point defect) left inside the "
        "inner cube cheb < 12",
        "plaquette map",
        val,
        "0 pierced plaquettes inside cheb 12 on all four free rows",
        ok,
    )

    tot = {
        "PASS": sum(c["verdict"] == "PASS" for c in CHECKS),
        "FAIL": sum(c["verdict"] == "FAIL" for c in CHECKS),
    }
    with open(OUT_JSON, "w") as fjs:
        json.dump(
            {
                "task": "M5.32 R22-1s audit",
                "totals": tot,
                "checks": CHECKS,
                "wall_s": round(time.time() - T0, 1),
            },
            fjs,
            indent=1,
        )
    print(f"TOTAL PASS {tot['PASS']} FAIL {tot['FAIL']} wall {time.time() - T0:.0f}s")


if __name__ == "__main__":
    main()
