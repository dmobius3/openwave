"""M5.32 R22-0 (a): the Faber identity on the certified stencils, and the
single-hedgehog self-energy (the known-answer control of R22-2).

EQUATIONS FIRST
---------------
Degenerate-pair vacuum, spatial block M_sp = delta I + (1 - delta) n n^T,
n in S^2, M_00 = +g constant (code branch s = -1: N_00 = -g; the R20 / R21 stack).

Claim under test (the author's 2026-09-15 post on the coordination thread,
section II): with d_i M = (1 - delta) (d_i n n^T + n d_i n^T) and n . d_i n = 0,
    F_ij = [d_i M, d_j M] = (1 - delta)^2 rho_ij [n]_x,
    rho_ij = n . (d_i n x d_j n),  E_top^k = (1/2) eps_kij rho_ij,
so tr(F_ij F_ij^T) = 2 (1 - delta)^4 rho_ij^2.

Two normalizations of the same density:
    author   sum over ordered (i, j) of <F_ij, F_ij>  = 4 (1 - delta)^4 |E_top|^2
    stack    4 sum over i < j of <F_ij, F_ij>_eta     = 8 (1 - delta)^4 |E_top|^2
The stack (m5_21_3_a_4d.e_parts, consumed read-only through
m5_32_r20_1_axes.density) carries the second, so in the author's
e^2 = 16 pi c4 (1 - delta)^4 the stack sits at c4 = 2: the single-hedgehog
self-energy outside r_c is 32 pi (1 - delta)^4 / r_c and the pair slope is
-/+ 64 pi (1 - delta)^4. Check A measures the factor, it does not assume it.

Checks (each can fail; the negative control shows it):
    A  pointwise, exact jets: random n, a_i perpendicular to n; the stack
       density formula against 8 (1 - delta)^4 |E_top|^2, to 1e-12.
       Negative control: a_i NOT perpendicular to n must break the identity.
    B  lattice residual on a smooth analytic director field (the author's
       stereographic plane-wave field, same seed), at three spacings, for
       the stack's production density (fwd / bwd averaged) and for the
       central stencil of the author's script; reference = the ANALYTIC
       8 (1 - delta)^4 |E_top|^2 from the analytic derivatives of n.
       Reported: relative L2 and max-norm residual, the fitted order.
       The author's field carries max |d n| h = 1.3 at n 40 (a director turn
       of more than a radian per cell near the zeros of w), so it is kept as
       the pre-asymptotic record (b_author_field) and the gate is read on the
       same field stretched by 8 (kscale 0.125), which the lattice resolves.
       PASS iff the L2 order of the stack density is within 1.5 to 2.5.
       Audit notes (m5_32_r22_0_audit.py, C1.3, C2.4, C3.4, C3.5): off the
       constraint the quartic is 8 (1 - delta)^4 (|E_top|^2 + 4 sum_{i<j}
       |gamma_j t_i - gamma_i t_j|^2), gamma_i = n . d_i n, and a one-sided
       difference carries gamma_i = -/+ (h / 2) |d_i n|^2. The fwd / bwd energy
       average cancels the O(h) term but keeps a positive-definite
       (h^2 / 4) |F_1|^2 that does not vanish where E_top does; the central
       stencil has no such term. So the ratio of the two error constants is
       field-dependent: 142 on this field (mean |E_top|^2 / mean |d n|^4 =
       0.016), 11 on the auditor's field (0.10), 0.3 to 0.5 on the hedgehog,
       where the central stencil is the worse one. The post's 6.5 and 1.7
       percent do come out of the author's field, at n 160 and n 320 (L2
       against the lattice right-hand side), not at the script's n 40.
    C  single hedgehog n = r-hat in the R21 box (L 48) at h 1.5, 1.0, 0.75:
       lattice energy on the shell r_c < r < R against the closed form
       32 pi (1 - delta)^4 (1 / r_c - 1 / R), and against the same-cell
       sum of the analytic density (separates stencil error from mask error).

Regenerate: python m5_32_r22_0a_faber_identity.py   (about 2 min, under 3 GB)
"""

import importlib.util
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r22_0a_faber_identity.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


R20 = _load("m5_32_r20_1_axes", "m5_32_r20_1_axes.py")
B3 = R20.B3
R0 = R20.R0
G = 8.0
DELTA = 0.3
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ================= A: pointwise =================
def stack_density_jets(n, a, delta):
    """4 sum_{i<j} <F_ij, F_ij>_eta on exact jets d_i M = (1-delta)(a_i n^T + n a_i^T)."""
    dM = []
    for ai in a:
        D = np.zeros((4, 4))
        D[1:, 1:] = (1.0 - delta) * (np.outer(ai, n) + np.outer(n, ai))
        dM.append(D)
    e = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            F = B3.comm_eta(dM[i], dM[j])
            e += 4.0 * float(B3.inner_eta(F, F))
    return e


def etop2_jets(n, a):
    rho = [n @ np.cross(a[1], a[2]), n @ np.cross(a[2], a[0]), n @ np.cross(a[0], a[1])]
    return float(sum(r * r for r in rho))


def check_a(samples=2000, seed=0):
    rng = np.random.default_rng(seed)
    worst, worst_neg, ratios = 0.0, np.inf, []
    for _ in range(samples):
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        delta = rng.uniform(0.05, 0.95)
        a = [rng.normal(size=3) for _ in range(3)]
        a_perp = [ai - n * (n @ ai) for ai in a]
        lhs = stack_density_jets(n, a_perp, delta)
        e2 = etop2_jets(n, a_perp)
        rhs = 8.0 * (1.0 - delta) ** 4 * e2
        worst = max(worst, abs(lhs - rhs) / max(abs(rhs), 1e-300))
        ratios.append(lhs / ((1.0 - delta) ** 4 * e2))
        # negative control: jets with a radial component, same formula
        lhs_n = stack_density_jets(n, a, delta)
        rhs_n = 8.0 * (1.0 - delta) ** 4 * etop2_jets(n, a)
        worst_neg = min(worst_neg, abs(lhs_n - rhs_n) / max(abs(rhs_n), 1e-300))
    out = {
        "samples": samples,
        "max_rel_residual": worst,
        "factor_min": float(np.min(ratios)),
        "factor_max": float(np.max(ratios)),
        "negative_control_median_note": "a_i with a component along n",
        "negative_control_min_rel_residual": float(worst_neg),
        "PASS": bool(worst < 1e-10),
    }
    log(
        f"A: max rel residual {worst:.2e}; factor in [{out['factor_min']:.12f}, {out['factor_max']:.12f}]"
        f" (stack normalization 8, the author's 4); negative control min residual {worst_neg:.2e}"
    )
    return out


# ================= B: lattice residual against the analytic field =================
def modes(seed=0, kscale=1.0):
    """The author's faber_identity.py field: six plane waves, the same generator order.
    kscale < 1 stretches the same field so that the lattice resolves it."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(6):
        c = rng.normal() + 1j * rng.normal()
        k = rng.normal(size=3) * 0.6 * kscale
        out.append((c, k))
    return out


def n_and_grad(X, Y, Z, md):
    """n (…,3) and analytic d_i n, list of three (…,3), from w = sum c exp(i k.x)."""
    w = np.zeros(X.shape, dtype=complex)
    dw = [np.zeros(X.shape, dtype=complex) for _ in range(3)]
    for c, k in md:
        ph = c * np.exp(1j * (k[0] * X + k[1] * Y + k[2] * Z))
        w += ph
        for i in range(3):
            dw[i] += 1j * k[i] * ph
    u, v = w.real, w.imag
    s = 1.0 + u * u + v * v
    n = np.stack([2 * u / s, 2 * v / s, (2.0 - s) / s], axis=-1)
    # partials of n with respect to u and v
    n_u = np.stack([2 / s - 4 * u * u / s**2, -4 * u * v / s**2, -4 * u / s**2], axis=-1)
    n_v = np.stack([-4 * u * v / s**2, 2 / s - 4 * v * v / s**2, -4 * v / s**2], axis=-1)
    dn = [n_u * dw[i].real[..., None] + n_v * dw[i].imag[..., None] for i in range(3)]
    return n, dn


def etop2_analytic(n, dn):
    rho = [
        np.einsum("...a,...a->...", n, np.cross(dn[1], dn[2])),
        np.einsum("...a,...a->...", n, np.cross(dn[2], dn[0])),
        np.einsum("...a,...a->...", n, np.cross(dn[0], dn[1])),
    ]
    return sum(r * r for r in rho)


def m_of_n(n, cfg, delta):
    M3 = delta * np.eye(3) + (1.0 - delta) * n[..., :, None] * n[..., None, :]
    return B3.embed34(M3, cfg)


def central_density(M, h):
    """The author's stencil: central differences, 4 sum_{i<j} <F,F>_eta (stack normalization)."""
    D = []
    for ax in range(3):
        d = np.zeros_like(M)
        sl_c = [slice(None)] * M.ndim
        sl_p = list(sl_c)
        sl_m = list(sl_c)
        sl_c[ax], sl_p[ax], sl_m[ax] = slice(1, -1), slice(2, None), slice(0, -2)
        d[tuple(sl_c)] = (M[tuple(sl_p)] - M[tuple(sl_m)]) / (2 * h)
        D.append(d)
    e = np.zeros(M.shape[:3])
    for i in range(3):
        for j in range(i + 1, 3):
            F = B3.comm_eta(D[i], D[j])
            e += 4.0 * B3.inner_eta(F, F)
    return e


def check_b(L=10.0, ns=(32, 64, 128), delta=DELTA, kscale=1.0):
    md = modes(0, kscale)
    rows = []
    for n_ in ns:
        cfg = B3.base_cfg(s=-1.0, g=G, n=n_, L=float(L), delta=delta)
        h = cfg["h"]
        X, Y, Z = B3.coords(n_, h)
        n, dn = n_and_grad(X, Y, Z, md)
        ref = 8.0 * (1.0 - delta) ** 4 * etop2_analytic(n, dn)
        gh = float(max(np.linalg.norm(d, axis=-1).max() for d in dn) * h)
        del dn
        M = m_of_n(n, cfg, delta)
        del n
        e_stack = R20.density(M, cfg, ("v4", R0.roots_of(cfg, degenerate=True), 0.0)) / h**3
        e_cent = central_density(M, h)
        del M
        k = max(2, int(round(0.08 * n_)))  # the same physical margin at every spacing
        inner = (slice(k, -k),) * 3
        r, es, ec = ref[inner], e_stack[inner], e_cent[inner]
        row = {
            "n": n_,
            "h": h,
            "max_grad_n_times_h": gh,
            "stack_rel_L2": float(np.sqrt(np.sum((es - r) ** 2) / np.sum(r**2))),
            "stack_rel_max": float(np.abs(es - r).max() / np.abs(r).max()),
            "stack_rel_sum": float((es.sum() - r.sum()) / r.sum()),
            "central_rel_L2": float(np.sqrt(np.sum((ec - r) ** 2) / np.sum(r**2))),
            "central_rel_max": float(np.abs(ec - r).max() / np.abs(r).max()),
            "central_rel_sum": float((ec.sum() - r.sum()) / r.sum()),
        }
        rows.append(row)
        log(
            f"B kscale {kscale} n {n_:4d} h {h:.4f} |dn|h {gh:.3f}: stack L2 {row['stack_rel_L2']:.3e} max {row['stack_rel_max']:.3e} "
            f"sum {row['stack_rel_sum']:+.3e} | central L2 {row['central_rel_L2']:.3e} "
            f"max {row['central_rel_max']:.3e} sum {row['central_rel_sum']:+.3e}"
        )
    hs = np.log([r["h"] for r in rows])
    order = {
        key: float(np.polyfit(hs, np.log([abs(r[key]) for r in rows]), 1)[0])
        for key in (
            "stack_rel_L2",
            "stack_rel_max",
            "stack_rel_sum",
            "central_rel_L2",
            "central_rel_max",
            "central_rel_sum",
        )
    }
    log("B orders: " + ", ".join(f"{k} {v:.2f}" for k, v in order.items()))
    return {
        "L": L,
        "kscale": kscale,
        "rows": rows,
        "order": order,
        "PASS": bool(1.5 <= order["stack_rel_L2"] <= 2.5),
    }


def check_b_author_setting(delta=DELTA):
    """The posted script's own setting (n 40, L 10, periodic roll, margin 2), its own norm."""
    n_, L = 40, 10.0
    h = L / n_
    a = (np.arange(n_) - n_ / 2 + 0.5) * h
    X, Y, Z = np.meshgrid(a, a, a, indexing="ij")
    rows = {}
    for label, nn in (("n40", 40), ("n80", 80)):
        h = L / nn
        a = (np.arange(nn) - nn / 2 + 0.5) * h
        X, Y, Z = np.meshgrid(a, a, a, indexing="ij")
        n, dn = n_and_grad(X, Y, Z, modes(0))
        cfg = B3.base_cfg(s=-1.0, g=G, n=nn, L=float(L), delta=delta)
        M = m_of_n(n, cfg, delta)
        lhs = central_density(M, h)
        # the author's right-hand side: E_top from central differences of n
        dnl = [np.gradient(n, h, axis=ax) for ax in range(3)]
        rhs_lat = 8.0 * (1.0 - delta) ** 4 * etop2_analytic(n, dnl)
        rhs_an = 8.0 * (1.0 - delta) ** 4 * etop2_analytic(n, dn)
        inner = (slice(2, -2),) * 3
        rows[label] = {
            "max_norm_vs_lattice_rhs": float(
                np.abs(lhs[inner] - rhs_lat[inner]).max() / np.abs(rhs_lat[inner]).max()
            ),
            "L2_vs_lattice_rhs": float(
                np.sqrt(np.sum((lhs[inner] - rhs_lat[inner]) ** 2) / np.sum(rhs_lat[inner] ** 2))
            ),
            "max_norm_vs_analytic": float(
                np.abs(lhs[inner] - rhs_an[inner]).max() / np.abs(rhs_an[inner]).max()
            ),
            "L2_vs_analytic": float(
                np.sqrt(np.sum((lhs[inner] - rhs_an[inner]) ** 2) / np.sum(rhs_an[inner] ** 2))
            ),
            "max_grad_n_times_h": float(max(np.linalg.norm(d, axis=-1).max() for d in dn) * h),
        }
        log(
            f"B' author setting {label}: "
            + ", ".join(f"{k} {v:.3e}" for k, v in rows[label].items())
        )
    return rows


# ================= C: the single hedgehog =================
def check_c(L=48.0, ns=(32, 48, 64), r_c=4.0, delta=DELTA):
    rows = []
    for n_ in ns:
        cfg = B3.base_cfg(s=-1.0, g=G, n=n_, L=float(L), delta=delta)
        h = cfg["h"]
        X, Y, Z = B3.coords(n_, h)
        r = np.sqrt(X * X + Y * Y + Z * Z)
        n = np.stack([X, Y, Z], axis=-1) / np.maximum(r, 1e-300)[..., None]
        M = m_of_n(n, cfg, delta)
        e = R20.density(M, cfg, ("v4", R0.roots_of(cfg, degenerate=True), 0.0))  # h^3-weighted
        R = L / 2.0 - 2.0 * h
        mask = (r > r_c) & (r < R)
        E_lat = float(e[mask].sum())
        E_cells = float((8.0 * (1.0 - delta) ** 4 / r[mask] ** 4).sum() * h**3)
        E_exact = 32.0 * np.pi * (1.0 - delta) ** 4 * (1.0 / r_c - 1.0 / R)
        rows.append(
            {
                "n": n_,
                "h": h,
                "r_c": r_c,
                "R": R,
                "E_lattice": E_lat,
                "E_same_cells_analytic": E_cells,
                "E_closed_form": E_exact,
                "stencil_rel": (E_lat - E_cells) / E_cells,
                "mask_rel": (E_cells - E_exact) / E_exact,
                "total_rel": (E_lat - E_exact) / E_exact,
            }
        )
        log(
            f"C n {n_} h {h:.3f}: E_lat {E_lat:.4f} same-cells {E_cells:.4f} closed {E_exact:.4f} | "
            f"stencil {rows[-1]['stencil_rel']:+.3%} mask {rows[-1]['mask_rel']:+.3%} total {rows[-1]['total_rel']:+.3%}"
        )
    ok = (
        abs(rows[-1]["stencil_rel"]) < abs(rows[0]["stencil_rel"])
        and abs(rows[-1]["stencil_rel"]) < 0.02
    )
    return {
        "rows": rows,
        "PASS": bool(ok),
        "gate": "the stencil residual falls with h and is under 2 percent at h 0.75",
    }


def main():
    res = {
        "a": check_a(),
        "b_author_setting": check_b_author_setting(),
        "b_author_field": check_b(kscale=1.0),
        "b": check_b(kscale=0.125),
        "c": check_c(),
    }
    res["constants"] = {
        "delta": DELTA,
        "one_minus_delta_4": (1 - DELTA) ** 4,
        "pair_slope_author_norm": 32 * np.pi * (1 - DELTA) ** 4,
        "pair_slope_stack_norm": 64 * np.pi * (1 - DELTA) ** 4,
    }
    res["runtime_s"] = time.time() - T0
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1)
    log(
        f"PASS A {res['a']['PASS']}  B {res['b']['PASS']}  C {res['c']['PASS']}  -> {os.path.basename(OUT_JSON)}"
    )


if __name__ == "__main__":
    main()
