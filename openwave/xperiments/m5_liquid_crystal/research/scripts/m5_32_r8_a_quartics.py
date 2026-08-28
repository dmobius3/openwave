"""M5.32 R8 arm a: class C5 / C6 (curvature^4, saturation, non-commutator
quartics) under the R7-B7 IR question.

Pre-registered hypothesis (before any number here):
  the omega^2 clock inertia is IR-divergent because its density is QUADRATIC in
  the jets while the hedgehog's orientation tail falls as 1/r (density r^-2,
  integral linear in L).  A quartic-in-Mdot term has a density QUARTIC in the
  jets (r^-4, convergent), so C5 / C6 is the only class whose clock inertia can
  be box-stable.  The structural prediction to test is that this does NOT save
  the fixed-J clock: with J = 2 C2 omega + 4 C4 omega^3, an extensive C2 and a
  finite C4 still give omega* -> J / (2 C2) ~ 1/L, and a negative C2 (forbidden
  by vacuum stability) would make omega*^2 = -C2 / (2 C4) GROW with L.

Stages:
  quartics <boxes>   the term build, covariance, the omega decomposition to
                     degree 4 on the REALIZED clock channel (the R7-N8 lesson),
                     the radial decay of each omega order's density, and the
                     box ladder of C2 and C4
  fixedj             the quartic-corrected fixed-J relation, solved per box
  collect            assemble + plots

Conventions: the run's locked table (eta = diag(-1,1,1,1), index 0 time, raw
contravariant internal entries, M -> L M L^T, certified sym stencil, h = L/n).
"""
import importlib.util, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA, PLOTS, CK = (os.path.join(RES, "data"), os.path.join(RES, "plots"),
                   os.path.join(RES, "checkpoints", "m5_32_r8"))
os.makedirs(CK, exist_ok=True)
OUT = os.path.join(DATA, "m5_32_r8_quartics.json")
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


RB = _load("m5_32_r2_b_bounded", "m5_32_r2_b_bounded.py")
B3, B8, LAG = RB.B3, RB.B8, RB.LAG
ETA = B3.ETA
G, S, DELTA = RB.G_MAIN, RB.S_MAIN, RB.DELTA
BOXES = ((32, 48.0), (48, 72.0), (64, 96.0))
OMEGAS = (0.25, 0.5)          # 5 samples with 0 and the negatives: exact to degree 4
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def dump(name, obj):
    with open(os.path.join(CK, name), "w") as f:
        json.dump(obj, f, indent=1)


def rd(name):
    p = os.path.join(CK, name)
    return json.load(open(p)) if os.path.exists(p) else None


# ---------------- quartic densities (per cell, from the jets) -------------
def d_I1(A):
    return LAG.density_from_K(LAG.F_of_A(A), LAG.REGISTRY["I1"]._K())


def d_I4(A):
    return LAG.density_from_K(LAG.F_of_A(A), LAG.REGISTRY["I4"]._K())


def _trAA(A):
    """t[mu,nu] per cell = tr(A_mu eta A_nu eta), the covariant jet gram."""
    Ae = np.einsum("m...ab,bc->m...ac", A, ETA)
    return np.einsum("m...ab,n...ba->mn...", Ae, Ae)


def d_C6a(A):
    """[sum_mu eta^mumu tr(A_mu eta A_mu eta)]^2."""
    t = _trAA(A)
    s = sum(ETA[m, m] * t[m, m] for m in range(4))
    return s * s


def d_C6b(A):
    """sum_{mu nu} eta^mumu eta^nunu [tr(A_mu eta A_nu eta)]^2."""
    t = _trAA(A)
    out = 0.0
    for m in range(4):
        for n in range(4):
            out = out + ETA[m, m] * ETA[n, n] * t[m, n] ** 2
    return out


def d_Fpair(A):
    """sum_{mu<nu, rho<sigma} eta-weighted <F_munu, F_rhosigma>_eta ^2.
    F is stored with the spatial indices FIRST: F[..., mu, nu, a, b]."""
    F = LAG.F_of_A(A)
    Fe = np.einsum("...ab,bc->...ac", F, ETA)      # F eta on the internal pair
    out = 0.0
    for m in range(4):
        for n in range(m + 1, 4):
            X = np.einsum("ab,...bc->...ac", ETA, Fe[..., m, n, :, :])
            for r in range(4):
                for s in range(r + 1, 4):
                    br = np.einsum("...ab,...ab->...", X, F[..., r, s, :, :])
                    w = ETA[m, m] * ETA[n, n] * ETA[r, r] * ETA[s, s]
                    out = out + w * br * br
    return out


def d_BI(A, b2=1.0e4):
    """Born-Infeld saturation b^2 (sqrt(1 + 2 d_I1 / b^2) - 1); the quartic
    piece is -d_I1^2 / (2 b^2)."""
    x = 1.0 + 2.0 * d_I1(A) / b2
    ok = x > 0
    out = np.zeros_like(x)
    out[ok] = b2 * (np.sqrt(x[ok]) - 1.0)
    return out, int((~ok).sum())


QUARTICS = {
    "Q_I1sq":  ("(F_abcd F^abcd)^2, the square of the certified curvature density (C5)",
                lambda A: d_I1(A) ** 2),
    "Q_I4sq":  ("(R_ac R^ac)^2, the square of the mixed-trace density (C5)",
                lambda A: d_I4(A) ** 2),
    "Q_Fpair": ("sum_{mu<nu,rho<sigma} <F_munu, F_rhosigma>_eta^2, a quartic that is not a square of a scalar (C5)",
                d_Fpair),
    "Q_C6a":   ("[sum_mu eta^mumu tr(d_mu M eta d^mu M eta)]^2 (C6)", d_C6a),
    "Q_C6b":   ("sum_{mu nu} eta^mumu eta^nunu [tr(d_mu M eta d_nu M eta)]^2 (C6)", d_C6b),
    "Q_BI":    ("Born-Infeld b^2 (sqrt(1 + 2 I1/b^2) - 1), b^2 = 1e4 (C5 saturation)",
                lambda A: d_BI(A)[0]),
    "I1":      ("the certified curvature (control, quadratic in omega)", d_I1),
}


def jets(M, cfg, a0, omega, branch):
    A = np.zeros((4,) + M.shape)
    for ax in range(3):
        A[1 + ax] = B3.d1(M, ax, cfg["h"], branch)
    A[0] = omega * a0
    return A


def density_omega(fn, M, cfg, a0, omega):
    """stencil-averaged density per cell at this omega."""
    out = None
    for br, wt in B3.branches(cfg["stencil"]):
        d = fn(jets(M, cfg, a0, omega, br))
        out = wt * d if out is None else out + wt * d
    return out


def omega_poly(fn, M, cfg, a0, omegas=OMEGAS):
    """exact degree-4 decomposition of the DENSITY: returns the per-cell
    coefficient arrays (c0, c2, c4) and the odd residual scale."""
    w1, w2 = omegas
    d0 = density_omega(fn, M, cfg, a0, 0.0)
    dp1 = density_omega(fn, M, cfg, a0, w1)
    dm1 = density_omega(fn, M, cfg, a0, -w1)
    dp2 = density_omega(fn, M, cfg, a0, w2)
    dm2 = density_omega(fn, M, cfg, a0, -w2)
    e1 = 0.5 * (dp1 + dm1) - d0          # c2 w1^2 + c4 w1^4
    e2 = 0.5 * (dp2 + dm2) - d0
    det = w1 ** 2 * w2 ** 4 - w2 ** 2 * w1 ** 4
    c2 = (e1 * w2 ** 4 - e2 * w1 ** 4) / det
    c4 = (e2 * w1 ** 2 - e1 * w2 ** 2) / det
    odd = 0.5 * (dp1 - dm1)
    return d0, c2, c4, odd


def radial_profile(dens, cfg, nbin=24):
    """shell-integrated |density| per unit r, and the fitted decay exponent of
    the density itself over the outer half."""
    r = cfg["R"] if "R" in cfg else None
    if r is None:
        r = RB.boost_geom(cfg)[0]
    rmax = float(cfg["h"]) * (cfg["n"] // 2)
    edges = np.linspace(0.0, rmax, nbin + 1)
    mid, shell, mean = [], [], []
    h3 = cfg["h"] ** 3
    for i in range(nbin):
        m = (r >= edges[i]) & (r < edges[i + 1])
        if not m.any():
            continue
        mid.append(0.5 * (edges[i] + edges[i + 1]))
        shell.append(float(np.sum(dens[m]) * h3 / (edges[i + 1] - edges[i])))
        mean.append(float(np.mean(np.abs(dens[m]))))
    mid, shell, mean = np.array(mid), np.array(shell), np.array(mean)
    sel = (mid > 0.35 * rmax) & (mid < 0.95 * rmax) & (mean > 0)
    exp_d = float(np.polyfit(np.log(mid[sel]), np.log(mean[sel]), 1)[0]) if sel.sum() > 3 else None
    return {"r": mid.tolist(), "shell_per_unit_r": shell.tolist(),
            "mean_abs_density": mean.tolist(), "density_exponent_outer": exp_d}


def stage_quartics(n, L):
    tag = f"n{n}_L{L:g}"
    have = rd(f"quart_{tag}.json")
    if have is not None:
        log(f"{tag}: restored")
        return have
    cfg = RB.cfg_of(n, L)
    M = B8.dressed(cfg, 0.0)                 # the undressed hedgehog
    a0 = B8.a0_unit(cfg, 0.0)                # the REALIZED clock channel
    h3 = cfg["h"] ** 3
    out = {"n": n, "L": L, "h": cfg["h"], "terms": {}}
    out["a0_stats"] = {"max_abs": float(np.abs(a0).max()),
                       "symmetric_rel": float(np.abs(a0 - np.swapaxes(a0, -1, -2)).max()
                                              / max(np.abs(a0).max(), 1e-300)),
                       "block_diagonal_time_row_max": float(np.abs(a0[..., 0, 1:]).max())}
    for name, (defn, fn) in QUARTICS.items():
        t0 = time.time()
        d0, c2, c4, odd = omega_poly(fn, M, cfg, a0)
        rec = {"definition": defn,
               "A_static": float(np.sum(d0) * h3),
               "C2": float(np.sum(c2) * h3),
               "C4": float(np.sum(c4) * h3),
               "odd_rel": float(np.abs(np.sum(odd) * h3) /
                                max(abs(float(np.sum(d0) * h3)), 1e-300)),
               "radial_static": radial_profile(d0, cfg),
               "radial_C2": radial_profile(c2, cfg),
               "radial_C4": radial_profile(c4, cfg),
               "runtime_s": round(time.time() - t0, 1)}
        out["terms"][name] = rec
        log(f"{tag} {name}: A {rec['A_static']:.6g} C2 {rec['C2']:.6g} "
            f"C4 {rec['C4']:.6g} | density exponents static "
            f"{rec['radial_static']['density_exponent_outer']} C2 "
            f"{rec['radial_C2']['density_exponent_outer']} C4 "
            f"{rec['radial_C4']['density_exponent_outer']} ({rec['runtime_s']}s)")
        dump(f"quart_{tag}.json", out)
    dump(f"quart_{tag}.json", out)
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "quartics"
    kw = dict(a.split("=", 1) for a in sys.argv[2:])
    if mode == "quartics":
        stage_quartics(int(kw.get("n", 32)), float(kw.get("L", 48)))
