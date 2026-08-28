"""M5.32 R9 arm b: is a string-free hedgehog possible at all?

The R8 audit found the author's ansatz M = Qh d4 Qh^T discontinuous along the z
axis.  Arm a showed the C5 quartics are IR-convergent once that line is excised,
which reopens the class, so the run now needs an instrument WITHOUT the defect.
This arm asks whether one can exist.

The structure of the ansatz answers it.  d4 = diag(-s g, 1, delta, 0) has three
DISTINCT spatial eigenvalues, and the frame Qh = R3(phi) R2(theta) carries the
eigenvalue-1 eigenvector to n-hat.  Fixing M then also requires placing the
PERPENDICULAR pair (delta, 0) in the tangent plane of the sphere, which is a
continuous unit tangent field on S^2: impossible by the hairy-ball theorem,
unless the perpendicular pair is DEGENERATE.  The pair is degenerate exactly at
delta = 0, and G1 (the clock generator) rotates that same pair, so at delta = 0
the clock flow vanishes identically.

Prediction, written before the numbers: a hedgehog is string-free if and only if
it has no clock.  Stages:

  1 equivalence  at delta = 0, M equals a manifestly frame-free expression in
                 n-hat n-hat^T (so the string is gone, exactly, not numerically)
  2 ladder       the phi-ring spread of M near the axis, the clock inertia, and
                 the axis-tube share of each C5 quartic, all as functions of
                 delta: they must vanish together
  3 consequence  the C5 coefficients measured at each delta, so the class's
                 "clean" numbers can be read off against the clock they cost
"""
import importlib.util, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA, CK = os.path.join(RES, "data"), os.path.join(RES, "checkpoints", "m5_32_r9")
os.makedirs(CK, exist_ok=True)
OUT = os.path.join(DATA, "m5_32_r9_string.json")
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


QA = _load("m5_32_r8_a_quartics", "m5_32_r8_a_quartics.py")
RB, B3, B8 = QA.RB, QA.B3, QA.B8
INS4 = B8.INS4
DELTAS = (0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0)
C5 = ("Q_I1sq", "Q_Fpair", "Q_I4sq")


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def cfg_delta(n, L, delta):
    c = RB.cfg_of(n, L)
    c = dict(c)
    c["delta"] = float(delta)
    return c


def hedgehog(cfg):
    """the author's ansatz at m = 0, t = 0, with cfg's delta."""
    return B8.dressed(cfg, 0.0)


def a0_of(cfg):
    return B8.a0_unit(cfg, 0.0)


# ---------------- stage 1: the frame-free identity at delta = 0 ------------
def stage_equivalence(n=32, L=48.0):
    cfg = cfg_delta(n, L, 0.0)
    M = hedgehog(cfg)
    X, Y, Z = INS4.coords(n, cfg["h"])
    R = np.sqrt(X * X + Y * Y + Z * Z)
    R = np.where(R == 0, 1.0, R)
    nx, ny, nz = X / R, Y / R, Z / R
    d4 = INS4.vac4(cfg)
    lam_par, lam_perp = float(d4[1, 1]), float(d4[2, 2])   # 1 and delta( = 0)
    # frame-free: M = diag(t, 0,0,0) + lam_perp (I3 - nn^T) + lam_par nn^T
    Mf = np.zeros_like(M)
    Mf[..., 0, 0] = d4[0, 0]
    nn = np.stack([nx, ny, nz], axis=-1)
    P = nn[..., :, None] * nn[..., None, :]
    Mf[..., 1:, 1:] = lam_perp * np.eye(3) + (lam_par - lam_perp) * P
    dev = float(np.abs(M - Mf).max())
    scale = float(np.abs(M).max())
    # the same comparison at delta != 0 must FAIL
    cfg2 = cfg_delta(n, L, 0.3)
    M2 = hedgehog(cfg2)
    d42 = INS4.vac4(cfg2)
    Mf2 = np.zeros_like(M2)
    Mf2[..., 0, 0] = d42[0, 0]
    Mf2[..., 1:, 1:] = float(d42[2, 2]) * np.eye(3) + (float(d42[1, 1]) - float(d42[2, 2])) * P
    dev2 = float(np.abs(M2 - Mf2).max())
    return {"delta_0_max_abs_deviation_from_frame_free_form": dev,
            "field_scale": scale, "relative": dev / scale,
            "delta_0.3_same_comparison": dev2,
            "form": "M = diag(d00,0,0,0) + lam_perp (I3 - n n^T) + (lam_par - lam_perp) n n^T",
            "reading": "at delta = 0 the ansatz IS a function of n n^T alone, so it is smooth away "
                       "from the origin and carries no string; at delta = 0.3 the same form fails, "
                       "because the third eigenvalue needs a tangent frame the sphere does not admit"}



# ---------------- stage 1b: the continuum ring (the string, measured) ------
def M_continuum(x, y, z, g, delta, s):
    """the ansatz evaluated OFF-lattice, so the ring can be shrunk freely."""
    def rot(Gm, a):
        G2m = Gm @ Gm
        return (np.eye(4)[None] + np.sin(a)[:, None, None] * Gm[None]
                + (1 - np.cos(a))[:, None, None] * G2m[None])
    rho = np.sqrt(x * x + y * y)
    phi, th = np.arctan2(y, x), -np.arctan2(z, rho)
    Qh = np.einsum("iab,ibc->iac", rot(B8.G3, phi), rot(B8.G2, th))
    d4 = np.diag([-s * g, 1.0, delta, 0.0])
    return np.einsum("iab,bc,idc->iad", Qh, d4, Qh)


def stage_continuum_ring(z=12.0, rhos=(1e-2, 1e-4, 1e-6), nphi=64):
    """spread of M around a shrinking phi ring about the z axis. A smooth field
    gives a spread proportional to the ring radius; a frame string leaves it
    finite and radius-independent."""
    out = {"z": z, "rhos": list(rhos), "rows": {}}
    for delta in DELTAS:
        row = []
        for rho in rhos:
            phi = np.linspace(0, 2 * np.pi, nphi, endpoint=False)
            M = M_continuum(rho * np.cos(phi), rho * np.sin(phi),
                            np.full_like(phi, z), RB.G_MAIN, delta, RB.S_MAIN)
            row.append(float(np.abs(M - M.mean(axis=0)).max()))
        out["rows"][f"delta_{delta:g}"] = {
            "spread_per_rho": row,
            "radius_ratio_1e-6_over_1e-2": (row[-1] / row[0]) if row[0] > 0 else None,
            "string": bool(row[0] > 0 and row[-1] / row[0] > 0.5),
            "note": "the ring shrinks by 1e4 between the first and last entry: a SMOOTH field's "
                    "spread shrinks with it (ratio 1e-4), a frame string leaves the ratio at 1",
            "spread_over_delta": (row[1] / delta) if delta > 0 else None}
        log(f"ring delta {delta:g}: spread {row} "
            f"({'STRING' if out['rows'][f'delta_{delta:g}']['string'] else 'smooth'})")
    return out


# ---------------- stage 2: the ladder in delta ----------------------------
def phi_ring_spread(M, cfg, n_rings=2, nz=(6, 12, 18)):
    """spread of M around a phi ring at the smallest available cylindrical radii,
    per z. A smooth field has this going to zero as the ring shrinks; a frame
    string leaves it finite."""
    n, h = cfg["n"], cfg["h"]
    X, Y, Z = INS4.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    pos = rho[rho > 1e-12]
    if pos.size == 0:
        return {}
    rmin = float(pos.min())
    out = {}
    for z in nz:
        for k in range(1, n_rings + 1):
            target = k * rmin
            sel = (np.abs(np.abs(Z) - z) < 0.51 * h) & (np.abs(rho - target) < 0.1 * rmin)
            if sel.sum() < 3:
                continue
            vals = M[sel]
            out[f"z_{z}_rho_{target:g}"] = float(np.abs(vals - vals.mean(axis=0)).max())
    return out


def tube_share(fn, M, cfg, a0, rho_cut_cells=2.0):
    d0, c2, c4, _ = QA.omega_poly(fn, M, cfg, a0)
    n, h = cfg["n"], cfg["h"]
    X, Y, _ = INS4.coords(n, h)
    rho = np.sqrt(X * X + Y * Y)
    m = rho < rho_cut_cells * h
    tot4, tot2 = float(np.sum(c4)), float(np.sum(c2))
    return {"C4": tot4 * h ** 3, "C2": tot2 * h ** 3,
            "C4_tube_share": float(np.sum(c4[m]) / tot4) if tot4 else None,
            "C2_tube_share": float(np.sum(c2[m]) / tot2) if tot2 else None}


def stage_ladder(n=32, L=48.0):
    rows = {}
    for delta in DELTAS:
        cfg = cfg_delta(n, L, delta)
        M, a0 = hedgehog(cfg), a0_of(cfg)
        h3 = cfg["h"] ** 3
        _, c2_i1, _, _ = QA.omega_poly(QA.d_I1, M, cfg, a0)
        kin = -4.0 * float(np.sum(c2_i1) * h3)
        rec = {"delta": delta, "kin": kin,
               "a0_max_abs": float(np.abs(a0).max()),
               "phi_ring_spread_LATTICE_CONFOUNDED": phi_ring_spread(M, cfg),
               "quartics": {t: tube_share(QA.QUARTICS[t][1], M, cfg, a0) for t in C5}}
        rows[f"delta_{delta:g}"] = rec
        sp = rec["phi_ring_spread_LATTICE_CONFOUNDED"]
        log(f"delta {delta:g}: kin {kin:10.4f}  |a0| {rec['a0_max_abs']:.4f}  "
            f"phi-ring spread max {max(sp.values() or [0.0]):.5f}  "
            f"Q_I1sq C4 {rec['quartics']['Q_I1sq']['C4']:.5g} "
            f"(tube {rec['quartics']['Q_I1sq']['C4_tube_share']})")
    # the joint statement: fit spread ~ delta^p and kin ~ delta^q
    ds = np.array([d for d in DELTAS if d > 0])
    sp = np.array([max(rows[f"delta_{d:g}"]["phi_ring_spread_LATTICE_CONFOUNDED"].values() or [0.0])
                   for d in ds])
    kn = np.array([rows[f"delta_{d:g}"]["kin"] for d in ds])
    kn0 = rows["delta_0"]["kin"]
    fit = {"spread_exponent_in_delta": float(np.polyfit(np.log(ds), np.log(sp), 1)[0]),
           "kin_exponent_in_delta": float(np.polyfit(np.log(ds), np.log(np.abs(kn - kn0)), 1)[0]),
           "kin_at_delta_0": kn0,
           "spread_at_delta_0": max(rows["delta_0"]["phi_ring_spread_LATTICE_CONFOUNDED"].values() or [0.0]),
           "lattice_ring_caveat": "the lattice ring is CONFOUNDED: its smallest radius is one cell, "
                                  "where the smooth hedgehog itself varies with phi, so it reads "
                                  "0.137 even at delta = 0 where the field is provably frame-free. "
                                  "The continuum ring is the clean instrument"}
    return {"rows": rows, "fit": fit}


if __name__ == "__main__":
    out = {"task": "M5.32 R9 arm b: a string-free hedgehog and the clock are mutually exclusive",
           "theorem": "M = Qh d4 Qh^T with three distinct spatial eigenvalues requires, beyond the "
                      "hedgehog direction n-hat, a continuous placement of the perpendicular pair in "
                      "the tangent plane of S^2, i.e. a nowhere-zero tangent vector field on the "
                      "sphere, which does not exist (hairy ball). The obstruction disappears exactly "
                      "when the perpendicular pair is degenerate, delta = 0, and the clock generator "
                      "G1 rotates that same pair, so a string-free hedgehog has a0 = 0 identically."}
    out["equivalence"] = stage_equivalence()
    log(f"equivalence: delta 0 deviation {out['equivalence']['relative']:.3e} relative; "
        f"delta 0.3 same form deviates {out['equivalence']['delta_0.3_same_comparison']:.4f}")
    out["continuum_ring"] = stage_continuum_ring()
    out["ladder"] = stage_ladder()
    out["runtime_s"] = round(time.time() - T0, 1)
    json.dump(out, open(OUT, "w"), indent=1)
    log(f"written {OUT}")
