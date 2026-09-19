"""M5.32 R22-2 adversarial audit: an independent attempt to refute the pair-gate claims.

Audited: m5_32_r22_2_pair.py, its rows (data/m5_32_r22_2_pair.json), its collect
(data/m5_32_r22_2_pair_collect.json) and the end fields data/m5_32_r22_2/<tag>.npz.

Nothing is imported from the audited script or from the stack. Own pieces:
    energy      4 (1 - delta)^4 h^3 sum_cells sum_{i<j} |[D_i P, D_j P]|^2, P = n n^T,
                D_i one-sided forward or backward (zero at the last / first cell), the two
                energies averaged; the central read uses second-order central differences
    degree      Van Oosterom-Strackee solid angles of the 12 triangles of every elementary
                lattice cube (n oriented as stored); a core degree is a sum over a cube
    G_D         (a) a 2D sine series with the exact sinh factor along the pair axis,
                (b) the image series of the grounded cube summed in neutral cells
    E_pred      the image series field (not the lattice DST of the audited script)
    Neumann     the image series with all-positive images, for the neutral unlike pair

Checks print as  P<id>.<n> PASS|FAIL ...  and land in data/m5_32_r22_2_audit.json.
A FAIL means the claim as worded was broken; the value field carries the counter-number.

Run: python m5_32_r22_2_audit.py      (single thread, about 1 minute, under 2 GB)
"""

import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
NPZ = os.path.join(DATA, "m5_32_r22_2")
ROWS_JSON = os.path.join(DATA, "m5_32_r22_2_pair.json")
COLLECT_JSON = os.path.join(DATA, "m5_32_r22_2_pair_collect.json")
OUT_JSON = os.path.join(DATA, "m5_32_r22_2_audit.json")

DELTA = 0.3
K = 8.0 * (1.0 - DELTA) ** 4
R_C = 2.5
SLOPE = 8.0 * np.pi * K
BOXES_H1 = ((48, 48, (6, 8, 10, 12)), (64, 64, (8, 12)))
BOX_H15 = (32, 48, (6, 8, 10, 12))
IMG = {"like": np.array([-1.0, -1.0, 1.0]), "unlike": np.array([-1.0, 1.0, 1.0])}
Q2 = {"ref": 0.0, "like": 1.0, "unlike": -1.0}
T0 = time.time()
CHECKS = []
_FIELDS = {}


def check(cid, claim, method, value, expected, ok):
    verdict = "PASS" if ok else "FAIL"
    CHECKS.append(
        {
            "id": cid,
            "claim": claim,
            "method": method,
            "value": value,
            "expected": expected,
            "verdict": verdict,
        }
    )
    print(
        f"{cid} {verdict} [{time.time() - T0:6.1f}s] {claim} | value {value} | "
        f"expected {expected}",
        flush=True,
    )


def rnd(v, p=4):
    if isinstance(v, dict):
        return {k: rnd(x, p) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [rnd(x, p) for x in v]
    if isinstance(v, (float, np.floating)):
        return float(np.round(v, p))
    if isinstance(v, (np.integer,)):
        return int(v)
    return v


def tag_of(kind, d, n_, big_l):
    return f"{kind}_d{d:g}_n{n_}_L{big_l:g}"


def field(kind, d, n_, big_l):
    t = tag_of(kind, d, n_, big_l)
    if t not in _FIELDS:
        _FIELDS[t] = np.load(os.path.join(NPZ, t + ".npz"))["n"]
    return _FIELDS[t]


def grid(n_, big_l):
    h = big_l / n_
    x = (np.arange(n_) - (n_ - 1) / 2.0) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    return h, x, X, Y, Z


# ================= own energy =================
def one_sided(P, ax, h, mode):
    D = np.zeros_like(P)
    hi = [slice(None)] * P.ndim
    lo = [slice(None)] * P.ndim
    hi[ax], lo[ax] = slice(1, None), slice(0, -1)
    diff = (P[tuple(hi)] - P[tuple(lo)]) / h
    D[tuple(lo if mode == "f" else hi)] = diff
    return D


def density(n, h, mode):
    P = n[..., :, None] * n[..., None, :]
    if mode == "c":
        D = [np.gradient(P, h, axis=ax) for ax in range(3)]
    else:
        D = [one_sided(P, ax, h, mode) for ax in range(3)]
    rho = np.zeros(n.shape[:3])
    for i in range(3):
        for j in range(i + 1, 3):
            C = np.einsum("...ab,...bc->...ac", D[i], D[j])
            C = C - np.swapaxes(C, -1, -2)
            rho += np.einsum("...ab,...ab->...", C, C)
    return 4.0 * (1.0 - DELTA) ** 4 * rho


def stack_density(n, h):
    return 0.5 * (density(n, h, "f") + density(n, h, "b"))


def etop_central(n, h):
    dn = [np.gradient(n, h, axis=ax) for ax in range(3)]

    def tr(a, b):
        return np.einsum("...a,...a->...", n, np.cross(dn[a], dn[b]))

    return np.stack([tr(1, 2), tr(2, 0), tr(0, 1)], axis=-1)


# ================= own degree reader =================
def tri_omega(a, b, c):
    num = np.einsum("...i,...i->...", a, np.cross(b, c))
    den = (
        1.0
        + np.einsum("...i,...i->...", a, b)
        + np.einsum("...i,...i->...", b, c)
        + np.einsum("...i,...i->...", a, c)
    )
    return 2.0 * np.arctan2(num, den)


CUBE_FACES = (
    ((1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)),
    ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0)),
    ((0, 1, 0), (0, 1, 1), (1, 1, 1), (1, 1, 0)),
    ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)),
    ((0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)),
    ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0)),
)


def cube_degrees(n):
    """degree of every elementary cube (corners i .. i+1), outward orientation."""
    N = n.shape[0]

    def corner(q):
        return n[q[0] : N - 1 + q[0], q[1] : N - 1 + q[1], q[2] : N - 1 + q[2]]

    tot = np.zeros((N - 1,) * 3)
    for f in CUBE_FACES:
        p = [corner(q) for q in f]
        tot += tri_omega(p[0], p[1], p[2]) + tri_omega(p[0], p[2], p[3])
    return tot / (4.0 * np.pi)


def link_angles(n, pin):
    """angles of nearest-neighbor links whose two cells are both free; degrees."""
    out = []
    for ax in range(3):
        hi = [slice(None)] * 3
        lo = [slice(None)] * 3
        hi[ax], lo[ax] = slice(1, None), slice(0, -1)
        dot = np.einsum("...i,...i->...", n[tuple(hi)], n[tuple(lo)])
        ang = np.degrees(np.arccos(np.clip(dot, -1.0, 1.0)))
        both_free = (~pin[tuple(hi)]) & (~pin[tuple(lo)])
        out.append(np.where(both_free, ang, 0.0))
    return out


# ================= own Green functions =================
def gd_series(d, a, M=600):
    """Dirichlet cube [-a, a]^3, source (d/2, 0, 0), field point (-d/2, 0, 0)."""
    big = 2.0 * a
    m = np.arange(1, M + 1)
    s = np.sin(m * np.pi / 2.0) ** 2
    g = np.pi * np.sqrt(m[:, None] ** 2 + m[None, :] ** 2) / big
    xl, xg = -d / 2.0 + a, d / 2.0 + a
    t = (
        np.exp(-g * (xg - xl))
        * (1 - np.exp(-2 * g * xl))
        * (1 - np.exp(-2 * g * (big - xg)))
        / (2 * g * (1 - np.exp(-2 * g * big)))
    )
    return float(4.0 * np.pi * (4.0 / big**2) * np.sum(s[:, None] * s[None, :] * t))


def axis_images(c, a, M, neumann=False):
    r = np.arange(-M, M + 1)
    pos = np.concatenate([4 * a * r + c, 4 * a * r + 2 * a - c])
    sg = np.concatenate([np.ones(len(r)), (1.0 if neumann else -1.0) * np.ones(len(r))])
    return pos, sg


def phi_regular(at, charges, a, M, neumann):
    """potential at `at` of every image of every charge, the coincident one left out."""
    tot = 0.0
    for q, src in charges:
        A = [axis_images(c, a, M, neumann) for c in src]
        DX = (at[0] - A[0][0])[:, None, None]
        DY = (at[1] - A[1][0])[None, :, None]
        DZ = (at[2] - A[2][0])[None, None, :]
        R = np.sqrt(DX**2 + DY**2 + DZ**2)
        S = A[0][1][:, None, None] * A[1][1][None, :, None] * A[2][1][None, None, :]
        with np.errstate(divide="ignore"):
            tot += q * float(np.where(R > 1e-9, S / R, 0.0).sum())
    return tot


def field_dirichlet(pts, src, a, M=2, chunk=1500):
    A = [axis_images(c, a, M) for c in src]
    P = np.stack(np.meshgrid(A[0][0], A[1][0], A[2][0], indexing="ij"), -1).reshape(-1, 3)
    S = (A[0][1][:, None, None] * A[1][1][None, :, None] * A[2][1][None, None, :]).ravel()
    E = np.zeros_like(pts)
    for i in range(0, len(pts), chunk):
        D = pts[i : i + chunk, None, :] - P[None, :, :]
        r3 = np.einsum("pqi,pqi->pq", D, D) ** 1.5
        E[i : i + chunk] = np.einsum("pqi,pq->pi", D, S[None, :] / r3)
    return E


def pin_even(d):
    return (
        -2.0
        * K
        * sum(
            (l / (l + 1.0)) * (4 * np.pi / (2 * l + 1)) * R_C ** (2 * l + 1) / d ** (2 * l + 2)
            for l in range(1, 40)
        )
    )


# ================= analytic seed (own formula) for the stencil test =================
def seed_field(X, Y, Z, d, kind):
    def w(cx):
        x = X - cx
        r = np.sqrt(x * x + Y * Y + Z * Z)
        return (x + 1j * Y) / np.maximum(r + Z, 1e-12)

    w1, w2 = w(d / 2.0), w(-d / 2.0)
    ww = {"ref": w1, "like": w1 * w2, "unlike": w1 * np.conj(w2)}[kind]
    a = np.abs(ww) ** 2
    return np.stack([2 * ww.real, 2 * ww.imag, 1 - a], -1) / (1 + a)[..., None]


def seed_exact_density(X, Y, Z, d, kind, e=1e-4):
    n = seed_field(X, Y, Z, d, kind)
    dn = []
    for a in range(3):
        sh = [e * (a == 0), e * (a == 1), e * (a == 2)]
        up = seed_field(X + sh[0], Y + sh[1], Z + sh[2], d, kind)
        dw = seed_field(X - sh[0], Y - sh[1], Z - sh[2], d, kind)
        dn.append((up - dw) / (2 * e))
    E = np.stack(
        [
            np.einsum("...a,...a->...", n, np.cross(dn[1], dn[2])),
            np.einsum("...a,...a->...", n, np.cross(dn[2], dn[0])),
            np.einsum("...a,...a->...", n, np.cross(dn[0], dn[1])),
        ],
        -1,
    )
    return K * (E**2).sum(-1)


# =====================================================================================
def p1(rows):
    claimed = {
        (48, 6): (6.503, -6.499),
        (48, 8): (4.372, -4.488),
        (48, 10): (3.160, -3.226),
        (48, 12): (2.369, -2.402),
        (64, 8): (4.805, -4.902),
        (64, 12): (2.771, -2.799),
    }
    U, worst_u, worst_row = {}, 0.0, 0.0
    for n_, big_l, ds in BOXES_H1 + (BOX_H15,):
        h = big_l / n_
        for d in ds:
            E = {}
            for kind in ("ref", "like", "unlike"):
                n = field(kind, d, n_, big_l)
                es = float(stack_density(n, h).sum() * h**3)
                ec = float(density(n, h, "c").sum() * h**3)
                E[kind] = (es, ec)
                row = rows[tag_of(kind, d, n_, big_l)]
                worst_row = max(
                    worst_row, abs(es - row["E"]), abs(ec - row["end_reads"]["E_central"])
                )
            U[(n_, d)] = {
                "like": E["like"][0] - 2 * E["ref"][0],
                "unlike": E["unlike"][0] - 2 * E["ref"][0],
                "like_c": E["like"][1] - 2 * E["ref"][1],
                "unlike_c": E["unlike"][1] - 2 * E["ref"][1],
                "E": E,
            }
            if (n_, d) in claimed:
                worst_u = max(
                    worst_u,
                    abs(U[(n_, d)]["like"] - claimed[(n_, d)][0]),
                    abs(U[(n_, d)]["unlike"] - claimed[(n_, d)][1]),
                )
    val = {f"n{k[0]}_d{k[1]}": [v["like"], v["unlike"]] for k, v in U.items() if k in claimed}
    check(
        "P1.1",
        "U(++) and U(+-) at h 1.0, six rows, own stack functional",
        "own one-sided commutator energy on P = n n^T, fwd and bwd averaged",
        rnd(val),
        "claimed values to 1e-3",
        worst_u < 1.5e-3,
    )
    check(
        "P1.2",
        "every one of the 30 row energies (stack and central) matches the JSON",
        "max |E_own - E_row| over 30 rows, both stencils",
        float(worst_row),
        "< 1e-8",
        worst_row < 1e-8,
    )
    return U


def p2(rows):
    unit_err, pin_err, self_err, counts, offsets = 0.0, 0.0, 0.0, {}, {}
    for n_, big_l, ds in BOXES_H1 + (BOX_H15,):
        h, x, X, Y, Z = grid(n_, big_l)
        for d in ds:
            c1, c2 = d / 2.0, -d / 2.0
            r1 = np.sqrt((X - c1) ** 2 + Y**2 + Z**2)
            r2 = np.sqrt((X - c2) ** 2 + Y**2 + Z**2)
            rh1 = np.stack([X - c1, Y, Z], -1) / r1[..., None]
            rh2 = np.stack([X - c2, Y, Z], -1) / r2[..., None]
            b1, b2 = r1 < R_C, r2 < R_C
            counts[f"n{n_}_d{d}"] = [int(b1.sum()), int(b2.sum())]
            offsets[f"n{n_}_d{d}"] = float((X[b1] - c1).mean())
            ball_self = {}
            for kind in ("ref", "like", "unlike"):
                n = field(kind, d, n_, big_l)
                unit_err = max(unit_err, float(np.abs(np.linalg.norm(n, axis=-1) - 1).max()))
                pin_err = max(pin_err, float(np.abs(n[b1] - rh1[b1]).max()))
                if kind != "ref":
                    pin_err = max(pin_err, float(np.abs(n[b2] - rh2[b2] * IMG[kind]).max()))
                # energy of the stencils that touch pinned cells only
                rf, rb = density(n, h, "f"), density(n, h, "b")
                for name, ball in (("b1", b1), ("b2", b2)):
                    if kind == "ref" and name == "b2":
                        continue
                    mf, mb = ball.copy(), ball.copy()
                    for ax in range(3):
                        mf &= np.roll(ball, -1, axis=ax)
                        mb &= np.roll(ball, 1, axis=ax)
                    ball_self[(kind, name)] = 0.5 * float(rf[mf].sum() + rb[mb].sum()) * h**3
            ref_self = ball_self[("ref", "b1")]
            for key, v in ball_self.items():
                self_err = max(self_err, abs(v - ref_self) / ref_self)
    check(
        "P2.1",
        "n is unit in every end field",
        "max ||n| - 1| over 30 fields",
        unit_err,
        "< 1e-12",
        unit_err < 1e-12,
    )
    check(
        "P2.2",
        "ball 1 = r-hat about +d/2, ball 2 = the stated constant O(3) image of r-hat "
        "about -d/2, in the END fields",
        "max |n - formula| over pinned cells, 30 fields",
        pin_err,
        "< 1e-12",
        pin_err < 1e-12,
    )
    check(
        "P2.3",
        "the pinned-only stencil energy of each ball of each pair equals the "
        "reference ball (self-energies cancel exactly)",
        "stack density summed over stencils whose four cells are all pinned; max rel diff",
        self_err,
        "< 1e-12",
        self_err < 1e-12,
    )
    bad = {k: v for k, v in offsets.items() if abs(v) > 1e-9}
    check(
        "P2.4",
        "the pinned ball is centered on the core in every row (needed for a "
        "like-for-like h sequence)",
        "centroid offset of ball 1 along x; cells per ball",
        {"offset_rows": rnd(bad), "cells_per_ball": counts},
        "offset 0 in every row",
        len(bad) == 0,
    )
    return counts


def p3(rows):
    summary, ok_deg, extra, flips = {}, True, 0, 0
    tube = {}
    smooth_h1 = 0.0
    for n_, big_l, ds in BOXES_H1 + (BOX_H15,):
        h, x, X, Y, Z = grid(n_, big_l)
        xc = x[:-1] + h / 2.0
        XC, YC, ZC = np.meshgrid(xc, xc, xc, indexing="ij")
        for d in ds:
            r1 = np.sqrt((X - d / 2.0) ** 2 + Y**2 + Z**2)
            r2 = np.sqrt((X + d / 2.0) ** 2 + Y**2 + Z**2)
            half = R_C + 2.0 * h
            in1 = (np.abs(XC - d / 2.0) < half) & (np.abs(YC) < half) & (np.abs(ZC) < half)
            in2 = (np.abs(XC + d / 2.0) < half) & (np.abs(YC) < half) & (np.abs(ZC) < half)
            for kind in ("ref", "like", "unlike"):
                n = field(kind, d, n_, big_l)
                pin = (r1 < R_C) | ((r2 < R_C) & (kind != "ref"))
                q = cube_degrees(n)
                deg = [float(q[in1].sum()), float(q[in2].sum()), float(q.sum())]
                n_def = int((np.abs(q) > 0.5).sum())
                want = {"ref": [1, 0, 1], "like": [1, 1, 2], "unlike": [1, -1, 0]}[kind]
                good = all(abs(a - b) < 1e-6 for a, b in zip(deg, want))
                n_extra = n_def - (1 if kind == "ref" else 2)
                extra += abs(n_extra)
                ang = link_angles(n, pin)
                over90 = int(sum((a > 90).sum() for a in ang))
                over60 = int(sum((a > 60).sum() for a in ang))
                over45 = int(sum((a > 45).sum() for a in ang))
                amax = float(max(a.max() for a in ang))
                flips += over90
                if h == 1.0:
                    ok_deg &= good
                    smooth_h1 = max(smooth_h1, amax)
                t = tag_of(kind, d, n_, big_l)
                summary[t] = {
                    "deg": rnd(deg, 6),
                    "extra_defects": n_extra,
                    "links_over_90": over90,
                    "links_over_60": over60,
                    "max_free_link_deg": round(amax, 2),
                }
                if over45:
                    ax_i = int(np.argmax([a.max() for a in ang]))
                    idx = np.argwhere(ang[1] > 60)
                    tube[t] = {
                        "links_over_45": over45,
                        "links_over_60": over60,
                        "max_deg": round(amax, 2),
                        "x_span_of_y_links": (
                            [float(x[idx[:, 0]].min()), float(x[idx[:, 0]].max())]
                            if len(idx)
                            else None
                        ),
                        "yz_index_span": (
                            [
                                int(idx[:, 1].min()),
                                int(idx[:, 1].max()),
                                int(idx[:, 2].min()),
                                int(idx[:, 2].max()),
                            ]
                            if len(idx)
                            else None
                        ),
                        "worst_axis": ax_i,
                    }
    check(
        "P3.1",
        "degrees on every h 1.0 row: like (+1, +1, 2), unlike (+1, -1, 0), ref +1",
        "solid-angle sum over triangulated elementary cubes, summed per core cube and box",
        {k: v["deg"] for k, v in summary.items() if "_n32_" not in k},
        "integers as claimed",
        ok_deg,
    )
    check(
        "P3.2",
        "no extra point defect in any of the 30 end fields",
        "count of elementary cubes with |degree| > 0.5 beyond the pinned cores",
        extra,
        0,
        extra == 0,
    )
    check(
        "P3.3",
        "no flip tube: no free link with neighboring directors more than 90 degrees "
        "apart, 30 rows",
        "nearest-neighbor angles, both cells free",
        flips,
        0,
        flips == 0,
    )
    check(
        "P3.4",
        "h 1.0 rows are smooth (max free-link angle under 45 degrees)",
        "same scan, h 1.0 rows",
        round(smooth_h1, 2),
        "< 45",
        smooth_h1 < 45.0,
    )
    check(
        "P3.5",
        "the h 1.5 unlike d 8 and d 10 anomalies are only resolution-limited "
        "(no tube-like object between the cores)",
        "free links over 45 and 60 degrees: count, max angle, x span and (y, z) index span "
        "of the y links over 60",
        tube,
        "no row with free links over 45 degrees",
        len(tube) == 0,
    )
    return summary


def p4(U):
    claimed = {
        (24, 6): 6.309,
        (24, 8): 4.313,
        (24, 10): 3.126,
        (24, 12): 2.345,
        (32, 8): 4.732,
        (32, 12): 2.739,
    }
    G, worst, img_err = {}, 0.0, 0.0
    for (a, d), c in claimed.items():
        g = gd_series(d, a)
        gi = phi_regular((-d / 2.0, 0, 0), [(1.0, (d / 2.0, 0, 0))], a, 30, False)
        G[(a, d)] = g
        img_err = max(img_err, abs(g - gi) / g)
        worst = max(worst, abs(SLOPE * g - c))
    check(
        "P4.1",
        "8 pi k G_D = 6.309/4.313/3.126/2.345 (half-width 24), 4.732/2.739 (32)",
        "2D sine series with exact sinh factor; image series as second method",
        {
            "values": rnd({f"a{a}_d{d}": SLOPE * g for (a, d), g in G.items()}),
            "series_vs_images_rel": img_err,
        },
        "claimed to 1e-3",
        worst < 1e-3 and img_err < 1e-6,
    )

    # wall components of E_top on the end fields, against the Dirichlet image field
    wall, zone, tq = {}, {}, {}
    for n_, big_l, ds in BOXES_H1:
        h, x, X, Y, Z = grid(n_, big_l)
        pts = np.stack([X, Y, Z], -1).reshape(-1, 3)
        r_in = np.zeros((n_,) * 3, bool)
        r_in[2:-2, 2:-2, 2:-2] = True
        outer = np.ones((n_,) * 3, bool)
        outer[6:-6, 6:-6, 6:-6] = False
        for d in ds:
            E1 = field_dirichlet(pts, (d / 2.0, 0, 0), big_l / 2.0).reshape(n_, n_, n_, 3)
            E2 = -E1[::-1, ::-1, ::-1]
            r1 = np.sqrt((X - d / 2.0) ** 2 + Y**2 + Z**2)
            r2 = np.sqrt((X + d / 2.0) ** 2 + Y**2 + Z**2)
            for kind in ("ref", "like", "unlike"):
                E = etop_central(field(kind, d, n_, big_l), h)
                Ep = E1 + Q2[kind] * E2
                m = r_in & (r1 > R_C + 2 * h) & ((r2 > R_C + 2 * h) | (kind == "ref"))
                sg = float(np.sign((E[m] * Ep[m]).sum()))
                dE = sg * E - Ep
                T = K * float((dE[m] ** 2).sum()) * h**3
                C = K * float((Ep[m] ** 2).sum()) * h**3
                corr = float(
                    (sg * E[m] * Ep[m]).sum() / np.sqrt((E[m] ** 2).sum() * (Ep[m] ** 2).sum())
                )
                slab = m & (np.abs(X) < 2.0) & (np.minimum(r1, r2) > 6.5)
                key = f"{kind}_n{n_}_d{d}"
                tq[key] = {
                    "T": T,
                    "C": C,
                    "T_over_C": T / C,
                    "corr": corr,
                    "T_midplane_slab": K * float((dE[slab] ** 2).sum()) * h**3,
                }
                zone[key] = K * float((dE[outer] ** 2).sum()) * h**3
                fr = {}
                for j in (2, 3):
                    en = et = pn = pt = 0.0
                    for ax in range(3):
                        for idx in (j, n_ - 1 - j):
                            sl = [slice(j, n_ - j)] * 3
                            sl[ax] = idx
                            F, Fp = E[tuple(sl)], Ep[tuple(sl)]
                            en += float((F[..., ax] ** 2).sum())
                            et += float((F**2).sum() - (F[..., ax] ** 2).sum())
                            pn += float((Fp[..., ax] ** 2).sum())
                            pt += float((Fp**2).sum() - (Fp[..., ax] ** 2).sum())
                    fr[f"j{j}"] = [et / (en + et), pt / (pn + pt)]
                wall[key] = fr
    dev = {
        k: max(abs(v["j2"][0] - v["j2"][1]), abs(v["j3"][0] - v["j3"][1])) for k, v in wall.items()
    }
    dev_rl = max(v for k, v in dev.items() if not k.startswith("unlike"))
    dev_un = max(v for k, v in dev.items() if k.startswith("unlike"))
    check(
        "P4.2",
        "ref and like: E_top at the free faces is normal as the grounded wall predicts",
        "tangential energy fraction on planes 2 and 3 cells inside the faces, data against "
        "the Dirichlet image field; max abs difference over ref and like rows",
        round(dev_rl, 4),
        "< 0.01",
        dev_rl < 0.01,
    )
    check(
        "P4.3",
        "unlike: E_top at the free faces is normal as the grounded wall predicts",
        "same read, unlike rows: [data, Dirichlet] tangential fractions; k |E - E_pred|^2 "
        "summed over the outer 6 layers",
        {
            "fractions": rnd({k: v for k, v in wall.items() if k.startswith("unlike")}),
            "outer_zone_deviation_energy": rnd({k: v for k, v in zone.items()}),
        },
        "max abs difference < 0.01",
        dev_un < 0.01,
    )

    # Neumann analog for the neutral pair: W_N(pair) - 2 W_D(ref)
    comp = {}
    for n_, big_l, ds in BOXES_H1:
        a = big_l / 2.0
        for d in ds:
            x1, x2 = (d / 2.0, 0, 0), (-d / 2.0, 0, 0)
            psi_d = phi_regular(x1, [(1.0, x1)], a, 30, False)
            phi_n = phi_regular(x1, [(1.0, x1), (-1.0, x2)], a, 30, True)
            u_n = 4 * np.pi * K * (2 * phi_n - 2 * psi_d)
            u_d = -SLOPE * G[(a, d)]
            comp[f"n{n_}_d{d}"] = {
                "U_unlike": U[(n_, d)]["unlike"],
                "U_dirichlet": u_d,
                "U_dirichlet_plus_pin": u_d + pin_even(d),
                "U_neumann": u_n,
                "U_free": -SLOPE / d,
            }
    rms = {
        m: float(np.sqrt(np.mean([(v["U_unlike"] - v[m]) ** 2 for v in comp.values()])))
        for m in ("U_dirichlet", "U_dirichlet_plus_pin", "U_neumann", "U_free")
    }
    check(
        "P4.4",
        "the grounded-wall (Dirichlet) cube is the right yardstick for U(+-), not "
        "Neumann and not free space",
        "rms of U(+-) minus each prediction over the six h 1.0 "
        "rows; Neumann = W_N(pair) - 2 W_D(ref) by all-positive images",
        {"rms": rnd(rms), "rows": rnd(comp)},
        "Dirichlet rms the smallest by a factor > 3",
        rms["U_dirichlet"] * 3 < min(rms["U_neumann"], rms["U_free"]),
    )
    # box to box at fixed d: cancels the near-core lattice error
    b2b = {}
    for d in (8, 12):
        pred = SLOPE * (G[(32.0, d)] - G[(24.0, d)])
        b2b[f"d{d}"] = {
            "pred": pred,
            "like": (U[(64, d)]["like"] - U[(48, d)]["like"]) / pred,
            "unlike": -(U[(64, d)]["unlike"] - U[(48, d)]["unlike"]) / pred,
        }
    worst_b = max(abs(v[k2] - 1) for v in b2b.values() for k2 in ("like", "unlike"))
    check(
        "P4.5",
        "the box-to-box shift of U at fixed d equals 8 pi k (G_D(32) - G_D(24))",
        "ratio measured / predicted, d 8 and 12, both pair types",
        rnd(b2b),
        "within 5 percent",
        worst_b < 0.05,
    )
    # bound-type prediction: U(+-) >= -8 pi k G_D + pin if the reference attains its bound
    viol = {k: v["U_unlike"] - v["U_dirichlet_plus_pin"] for k, v in comp.items()}
    check(
        "P4.6",
        "8 pi k G_D is a lower bound: U(+-) is not below -(8 pi k G_D) + pin",
        "U(+-) - (-(8 pi k G_D) + pin) per h 1.0 row (negative = the lattice pair sits below "
        "the continuum bound: a direct measure of the lattice error)",
        rnd(viol),
        ">= 0 on every row",
        min(viol.values()) >= 0,
    )
    return G, tq, viol


def p5(U, collect):
    claimed = {6: (4.88, -6.49), 8: (3.71, -4.55), 10: (2.82, -3.20), 12: (2.16, -2.38)}
    worst = max(
        max(abs(U[(48, d)]["like_c"] - c[0]), abs(U[(48, d)]["unlike_c"] - c[1]))
        for d, c in claimed.items()
    )
    check(
        "P5.1",
        "central-stencil U of the same end fields: U(+-) -6.49/-4.55/-3.20/-2.38, "
        "U(++) +4.88/+3.71/+2.82/+2.16",
        "own central-difference commutator energy",
        rnd({d: [U[(48, d)]["like_c"], U[(48, d)]["unlike_c"]] for d in claimed}),
        "claimed to 0.01",
        worst < 0.011,
    )
    # where the like excess (stack minus central) sits
    n_, big_l = 48, 48
    h, x, X, Y, Z = grid(n_, big_l)
    where = {}
    for d in (6, 8, 10, 12):
        ref = field("ref", d, n_, big_l)
        rs, rc = stack_density(ref, h), density(ref, h, "c")
        rs2, rc2 = rs[::-1, ::-1, ::-1], rc[::-1, ::-1, ::-1]
        r1 = np.sqrt((X - d / 2.0) ** 2 + Y**2 + Z**2)
        r2 = np.sqrt((X + d / 2.0) ** 2 + Y**2 + Z**2)
        rm = np.minimum(r1, r2)
        for kind in ("like", "unlike"):
            n = field(kind, d, n_, big_l)
            D = (stack_density(n, h) - rs - rs2) - (density(n, h, "c") - rc - rc2)
            D *= h**3
            tot = float(D.sum())
            far = rm >= 4.5
            edge = (np.abs(X) > 21.9) | (np.abs(Y) > 21.9) | (np.abs(Z) > 21.9)
            where[f"{kind}_d{d}"] = {
                "total": tot,
                "cores_r_lt_4.5": float(D[~far].sum()),
                "fold_slab_absx_lt_1.5_far": float(D[far & (np.abs(X) < 1.5)].sum()),
                "walls_outer_2_layers": float(D[edge].sum()),
                "bulk_rest": float(D[far & (np.abs(X) >= 1.5) & ~edge].sum()),
                "z_below": float(D[Z < 0].sum()),
                "z_above": float(D[Z > 0].sum()),
            }
    zasym = max(
        abs(v["z_below"] - v["z_above"]) / max(abs(v["total"]), 1e-9)
        for k2, v in where.items()
        if k2.startswith("like")
    )
    check(
        "P5.2",
        "the like excess does not sit on the stereographic string (z < 0 half equals "
        "z > 0 half)",
        "stack-minus-central interaction density, halves in z",
        round(zasym, 4),
        "< 0.02 of the total",
        zasym < 0.02,
    )
    shares = {
        k2: {
            "cores": v["cores_r_lt_4.5"] / v["total"],
            "fold": v["fold_slab_absx_lt_1.5_far"] / v["total"],
            "walls": v["walls_outer_2_layers"] / v["total"],
            "bulk": v["bulk_rest"] / v["total"],
        }
        for k2, v in where.items()
        if k2.startswith("like")
    }
    frac_bulk = max(v["bulk"] for v in shares.values())
    check(
        "P5.3",
        "the like excess is localized (cores, midplane fold, walls) rather than "
        "spread through the bulk",
        "shares of the like stack-minus-central difference: "
        "cores r < 4.5, fold slab |x| < 1.5, outer 2 layers, bulk = the rest",
        {"shares": rnd(shares, 3), "map": rnd(where)},
        "bulk share < 0.5 on every d",
        frac_bulk < 0.5,
    )
    # which stencil reads a smooth degree-2 texture better: the analytic seed, known density
    seed = {}
    for h_, n2 in ((1.0, 48), (0.5, 96)):
        _, x2, X2, Y2, Z2 = grid(n2, 48)
        inner = np.zeros(X2.shape, bool)
        inner[2:-2, 2:-2, 2:-2] = True
        for d in (6, 12):
            r1 = np.sqrt((X2 - d / 2.0) ** 2 + Y2**2 + Z2**2)
            r2 = np.sqrt((X2 + d / 2.0) ** 2 + Y2**2 + Z2**2)
            m = inner & (r1 > 4.0) & (r2 > 4.0)
            res = {}
            for kind in ("ref", "like", "unlike"):
                n = seed_field(X2, Y2, Z2, d, kind)
                ex = float(seed_exact_density(X2, Y2, Z2, d, kind)[m].sum()) * h_**3
                st = float(stack_density(n, h_)[m].sum()) * h_**3
                ce = float(density(n, h_, "c")[m].sum()) * h_**3
                res[kind] = (ex, st, ce)
            for kd in ("like", "unlike"):
                u = [res[kd][i] - 2 * res["ref"][i] for i in range(3)]
                seed[f"h{h_}_d{d}_{kd}"] = {
                    "exact": u[0],
                    "stack_err": u[1] - u[0],
                    "central_err": u[2] - u[0],
                }
            seed[f"h{h_}_d{d}_single_rel"] = {
                kd: [
                    (res[kd][1] - res[kd][0]) / res[kd][0],
                    (res[kd][2] - res[kd][0]) / res[kd][0],
                ]
                for kd in res
            }
            del n
    like_ratio = min(
        abs(seed[f"h1.0_d{d}_like"]["central_err"]) / abs(seed[f"h1.0_d{d}_like"]["stack_err"])
        for d in (6, 12)
    )
    sigma = max(
        abs(seed[f"h1.0_d{d}_{kd}"]["stack_err"]) for d in (6, 12) for kd in ("like", "unlike")
    )
    check(
        "P5.4",
        "the stack number is the better estimate of the continuum like-pair energy",
        "analytic product seed (known k |E_top|^2) outside r 4 of both cores: error of the "
        "pair-minus-two-singles energy per stencil at h 1.0 and h 0.5; [stack, central] "
        "relative errors of each single energy",
        {
            "min_central_over_stack_error_like_h1": round(like_ratio, 2),
            "max_stack_abs_err_h1": round(sigma, 4),
            "table": rnd(seed),
        },
        "central error at least 3 times the stack error",
        like_ratio >= 3.0,
    )
    return where, seed, sigma


def p6(tq):
    claimed_like, claimed_unlike = (0.14, 0.21), (0.07, 0.09)
    like = [v["T"] for k2, v in tq.items() if k2.startswith("like")]
    unlike = [v["T"] for k2, v in tq.items() if k2.startswith("unlike")]
    ok = (
        claimed_like[0] - 0.005 <= min(like)
        and max(like) <= claimed_like[1] + 0.005
        and claimed_unlike[0] - 0.005 <= min(unlike)
        and max(unlike) <= claimed_unlike[1] + 0.005
    )
    check(
        "P6.1",
        "E_transverse is 0.07 to 0.09 (unlike) and 0.14 to 0.21 (like) at h 1.0",
        "own image-series E_pred, own central E_top, same mask",
        rnd({"like": [min(like), max(like)], "unlike": [min(unlike), max(unlike)]}),
        "inside the claimed ranges",
        ok,
    )
    floor = 2.0 * np.mean([v["T"] for k2, v in tq.items() if k2.startswith("ref")])
    r_un = [u / floor for u in unlike]
    r_li = [v / floor for v in like]
    check(
        "P6.2",
        "unlike realizes the box Coulomb field, like keeps a small transverse part: "
        "measured against the floor of the read (two times the single-core value)",
        "T / (2 T_ref) per row",
        rnd({"floor": floor, "unlike": r_un, "like": r_li}),
        "unlike within 0.75 to 1.25 of the floor, like above 1.5",
        0.75 <= min(r_un) and max(r_un) <= 1.25 and min(r_li) > 1.5,
    )
    toc = {k2: v["T_over_C"] for k2, v in tq.items()}
    un_toc = [v for k2, v in toc.items() if k2.startswith("unlike")]
    li_toc = [v for k2, v in toc.items() if k2.startswith("like")]
    check(
        "P6.3",
        "the same statement in relative terms: unlike is closer to its Coulomb field " "than like",
        "T / (Coulomb energy in the same mask) and the pointwise correlation",
        rnd({"T_over_C": toc, "corr": {k2: v["corr"] for k2, v in tq.items()}}, 5),
        "max unlike T/C < min like T/C",
        max(un_toc) < min(li_toc),
    )
    return floor


def p7(rows, U, viol):
    rem = {}
    for n_, big_l, ds in BOXES_H1:
        for d in ds:
            for kind in ("ref", "like", "unlike"):
                row = rows[tag_of(kind, d, n_, big_l)]
                tr = np.array(row["trace"])
                it_end = tr[-1, 0]

                def e_at(it):
                    return float(tr[np.argmin(np.abs(tr[:, 0] - it)), 1])

                d1 = e_at(it_end - 1000) - e_at(it_end - 500)
                d2 = e_at(it_end - 500) - row["E"]
                rho = min(max(d2 / d1, 0.0), 0.95) if d1 > 0 else 0.95
                geo = d2 * rho / (1 - rho)
                pw = d2 * (it_end - 500) / 500.0  # E - E_inf = B / it
                rem[tag_of(kind, d, n_, big_l)] = {
                    "geometric": geo,
                    "power_law_p1": pw,
                    "residual_end": row["residual_end"],
                    "label": row["label"][:7],
                }
    prop = {}
    for n_, big_l, ds in BOXES_H1:
        for d in ds:
            g = {kd: rem[tag_of(kd, d, n_, big_l)] for kd in ("ref", "like", "unlike")}
            prop[f"n{n_}_d{d}"] = {
                m: [-(g["like"][m] + g["unlike"][m]) / 2.0, 2.0 * g["ref"][m]]
                for m in ("geometric", "power_law_p1")
            }
    worst_relax = max(
        max(abs(v["power_law_p1"][0]), abs(v["power_law_p1"][1])) for v in prop.values()
    )
    at_gate = sum(1 for v in rem.values() if v["label"] == "AT_GATE")
    check(
        "P7.1",
        "relaxation error of U_grav is small (the rows stopped at the iteration cap)",
        "last 1000 iterations of each trace: geometric tail and the slower 1 / it tail, "
        "propagated as [from pairs, from 2 x ref]",
        {
            "U_grav_shift_bounds": rnd(prop, 5),
            "rows_at_gate_of_18": at_gate,
            "max_residual_end": max(v["residual_end"] for v in rem.values()),
        },
        "max |shift| < 0.01",
        worst_relax < 0.01,
    )
    remainder, ug = {}, {}
    for d in (6, 8, 10, 12):
        ug[d] = 0.5 * (U[(48, d)]["like"] + U[(48, d)]["unlike"])
        remainder[d] = ug[d] - pin_even(d)
    worst_rem = max(abs(v) for v in remainder.values())
    check(
        "P7.2",
        "U_grav is zero within the pin artifact plus the relaxation error",
        "U_grav - pin multipole sum per d (L48), against the relaxation bound of P7.1",
        rnd({"U_grav": ug, "remainder": remainder, "relaxation_bound": worst_relax}, 4),
        "|remainder| <= relaxation bound",
        worst_rem <= worst_relax,
    )
    # the lattice error: three independent handles
    cen = {d: 0.5 * (U[(48, d)]["like_c"] + U[(48, d)]["unlike_c"]) for d in (6, 8, 10, 12)}
    h15 = {d: 0.5 * (U[(32, d)]["like"] + U[(32, d)]["unlike"]) for d in (6, 12)}
    rich = {d: ug[d] + (ug[d] - h15[d]) / (1.5**2 - 1.0) for d in (6, 12)}
    margin = {
        d: U[(48, d)]["like"] - (SLOPE * gd_series(d, 24.0) + pin_even(d)) for d in (6, 8, 10, 12)
    }
    deficit = {d: viol[f"n48_d{d}"] for d in (6, 8, 10, 12)}
    bound = {
        d: max(abs(margin[d]), abs(deficit[d]), abs(rich[d] - ug[d]) if d in rich else 0.0)
        for d in (6, 8, 10, 12)
    }
    bound_rel = {d: bound[d] / (SLOPE * gd_series(d, 24.0)) for d in bound}
    check(
        "P7.3",
        "U_grav at h 1.0 is resolved: its value exceeds its own lattice dependence",
        "U_grav by central stencil; U_grav at h 1.5 (centered-ball rows d 6, 12) and the h^2 "
        "extrapolation; the like margin over and the unlike deficit under the continuum "
        "bound, whose half-sum is U_grav - pin",
        rnd(
            {
                "U_grav_stack_h1": ug,
                "U_grav_central_h1": cen,
                "U_grav_stack_h15": h15,
                "U_grav_h2_extrapolated": rich,
                "like_margin_over_bound": margin,
                "unlike_deficit_under_bound": deficit,
                "upper_bound_abs": bound,
                "upper_bound_over_U_box": bound_rel,
            }
        ),
        "|U_grav| > lattice dependence on every d",
        all(abs(ug[d]) > bound[d] for d in ug),
    )
    return bound


def p8(U, G, collect, sigma):
    fits = {}
    for n_, a, ds in ((48, 24.0, (6, 8, 10, 12)), (64, 32.0, (8, 12))):
        gd = np.array([G[(a, d)] for d in ds])
        inv = 1.0 / np.array(ds, float)
        for key in ("like", "unlike"):
            u = np.array([U[(n_, d)][key] for d in ds])
            s, c = np.polyfit(gd, u, 1)
            s1, c1 = np.polyfit(inv, u, 1)
            fits[f"n{n_}_{key}"] = {
                "slope_vs_G_D": s / SLOPE,
                "intercept_vs_G_D": c,
                "slope_vs_inv_d": s1 / SLOPE,
                "intercept_vs_inv_d": c1,
                "slope_through_origin_vs_G_D": float((u * gd).sum() / (gd**2).sum() / SLOPE),
            }
    cb = collect["boxes"]["n48_L48"]
    dif = max(
        abs(fits["n48_like"]["slope_vs_G_D"] - cb["fit_U_like"]["slope_over_8pik_vs_G_D"]),
        abs(fits["n48_unlike"]["slope_vs_inv_d"] - cb["fit_U_unlike"]["slope_over_8pik_vs_inv_d"]),
        abs(fits["n48_like"]["intercept_vs_G_D"] - cb["fit_U_like"]["intercept_vs_G_D"]),
    )
    check(
        "P8.1",
        "the collect fits reproduce (own G_D, own U)",
        "np.polyfit, L48 and L64",
        rnd(fits),
        "collect values to 2e-3",
        dif < 2e-3,
    )
    box_aware = [abs(abs(fits[k2]["slope_vs_G_D"]) - 1.0) for k2 in fits]
    check(
        "P8.2",
        "the slope is 8 pi k within about 3 percent in the box-aware fit too (the "
        "claim quotes the slope from the 1/d fit and the intercept from the G_D fit)",
        "|slope / 8 pi k| - 1 against G_D, like and unlike, both boxes",
        rnd(box_aware),
        "all <= 0.035",
        max(box_aware) <= 0.035,
    )
    odd = {
        f"n{n_}_d{d}": 0.5 * (U[(n_, d)]["like"] - U[(n_, d)]["unlike"]) / (SLOPE * G[(a, d)])
        for n_, a, ds in ((48, 24.0, (6, 8, 10, 12)), (64, 32.0, (8, 12)))
        for d in ds
    }
    check(
        "P8.3",
        "pointwise: the charge-odd energy U_odd equals 8 pi k G_D within 3.5 percent "
        "on all six h 1.0 points",
        "U_odd / (8 pi k G_D)",
        rnd(odd),
        "0.965 to 1.035",
        all(0.965 <= v <= 1.035 for v in odd.values()),
    )
    # can the data tell 8 pi k from a slope 5 percent off? zero intercept, pin removed
    chi = {}
    for s in (0.95, 1.00, 1.05):
        c2 = 0.0
        for d in (6, 8, 10, 12):
            ub = SLOPE * G[(24.0, d)]
            c2 += ((U[(48, d)]["like"] - pin_even(d) - s * ub) / sigma) ** 2
            c2 += ((U[(48, d)]["unlike"] - pin_even(d) + s * ub) / sigma) ** 2
        chi[f"s{s:.2f}"] = c2
    check(
        "P8.4",
        "the data distinguish slope 1.00 from slope 1.05 (times 8 pi k)",
        "chi^2 over the 8 L48 points, zero intercept, pin removed, sigma = the largest "
        "stack truncation error of the seed test (P5.4)",
        rnd({"sigma": sigma, "chi2": chi}, 3),
        "|chi2(1.05) - chi2(1.00)| > 4",
        abs(chi["s1.05"] - chi["s1.00"]) > 4.0,
    )
    check(
        "P8.5",
        "the data distinguish slope 1.00 from slope 0.95",
        "same chi^2",
        rnd(chi, 3),
        "chi2(0.95) - chi2(1.00) > 4",
        chi["s0.95"] - chi["s1.00"] > 4.0,
    )
    return fits


def refinement():
    tags = [tag_of(kd, 8, 64, 48) for kd in ("ref", "like", "unlike")]
    if not all(os.path.exists(os.path.join(NPZ, t + ".npz")) for t in tags):
        print("R.0 refinement triple *_d8_n64_L48 not present at audit time: not used")
        return None
    h = 48.0 / 64
    E = {}
    for kd in ("ref", "like", "unlike"):
        n = field(kd, 8, 64, 48)
        E[kd] = (float(stack_density(n, h).sum() * h**3), float(density(n, h, "c").sum() * h**3))
    out = {
        "U_like": E["like"][0] - 2 * E["ref"][0],
        "U_unlike": E["unlike"][0] - 2 * E["ref"][0],
        "U_like_central": E["like"][1] - 2 * E["ref"][1],
        "U_unlike_central": E["unlike"][1] - 2 * E["ref"][1],
    }
    out["U_grav"] = 0.5 * (out["U_like"] + out["U_unlike"])
    print("R.0 refinement triple at h 0.75 (d 8, L48):", rnd(out))
    return out


def main():
    with open(ROWS_JSON) as f:
        rows = json.load(f)["rows"]
    with open(COLLECT_JSON) as f:
        collect = json.load(f)
    U = p1(rows)
    p2(rows)
    p3(rows)
    G, tq, viol = p4(U)
    where, seed, sigma = p5(U, collect)
    p6(tq)
    bound = p7(rows, U, viol)
    p8(U, G, collect, sigma)
    ref = refinement()
    n_pass = sum(1 for c in CHECKS if c["verdict"] == "PASS")
    out = {
        "task": "M5.32 R22-2 adversarial audit",
        "checks": CHECKS,
        "totals": {"PASS": n_pass, "FAIL": len(CHECKS) - n_pass, "checks": len(CHECKS)},
        "U_grav_upper_bound_abs": rnd({str(k2): v for k2, v in bound.items()}),
        "refinement_h075": rnd(ref) if ref else None,
        "runtime_s": round(time.time() - T0, 1),
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1, default=float)
    print(
        f"TOTAL PASS {n_pass} FAIL {len(CHECKS) - n_pass} of {len(CHECKS)} "
        f"({time.time() - T0:.0f}s)"
    )


if __name__ == "__main__":
    main()
