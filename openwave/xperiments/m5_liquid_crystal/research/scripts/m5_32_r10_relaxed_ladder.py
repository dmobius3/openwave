"""M5.32 R10: does the clock inertia stay extensive on the RELAXED, core-resolved field?

R9's audit resolved the ansatz's axis line into a finite physical core and found
the clock survives with an h-convergent inertia (351.17 at h = 1.5, 351.14 at
h = 0.75, against 426.51 and 445.22 rigid).  One question decides whether the
run's whole obstruction is physics or an artifact of the rigid ansatz:

    is the RELAXED inertia still LINEAR IN THE BOX?

PRE-REGISTERED PREDICTION (written before any number here): core resolution is a
LOCAL repair, radius about 4, and R7 measured only 6 % of the inertia inside
r = 6.  So the extensive SLOPE should be unchanged and only the INTERCEPT should
move.  Two independent tests, because they fail differently:

  A  SHELL TEST (cheap, no convergence confound): relax one box and compare the
     kin shell profile rigid vs relaxed.  If the shells beyond the core are
     unchanged, the slope is unchanged by construction.  PREDICT: shells beyond
     r = 8 agree to a few percent, and the whole 75.3 drop sits inside r ~ 6.
  B  LADDER (direct, but confounded): relax at L = 48 and 72 under an identical
     protocol and fit the slope.  PREDICT kin_relaxed(72) = 553.7 +- 15.
     CAVEAT registered in advance: FIRE stops at max_iter, not convergence, and
     relaxation time grows with the box, so the larger box is LESS relaxed under
     an identical budget.  Test B is therefore read only WITH its iteration
     ladder (1500 / 3000 / 6000), never from a single endpoint.

PROTOCOL (the R9 auditor's, verbatim, so the 351.17 is reproduced before anything
is built on it): cfg = base_cfg(s = -1, g = 8, delta = 0.3) at h = 1.5;
M0 = the m = 0 ansatz; free = ~pin_shell(n, h) at the default depth 1.6;
fire(dt0 = 0.01, dt_max = 0.1, max_iter = 3000), a0 = None and omega = 0 so the
relaxed functional is pure E_static; kin measured on the relaxed field with the
FROZEN ansatz a0 computed once from M0.

WHY g = 8 (the auditor's finding, carried as a scope limit): E_u and kin are
EXACTLY g-independent on this family (62.8517443315 and 426.5070121484 at both
g = 8 and g = 32, ten digits) because M_00 is a spatial constant and M_0i = 0, so
g cancels from the jets.  V4 is NOT g-independent and is the whole stiffness:
4096x larger at g = 32.  The core melt trades u-energy against V4, so the
EXISTENCE of a resolved core and a surviving clock is established at g = 8, while
the core radius and 351.17 are g = 8 numbers.  Stage `gcore` measures how much
the g = 32 core tightens.
"""
import importlib.util, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA, CK = os.path.join(RES, "data"), os.path.join(RES, "checkpoints", "m5_32_r10")
os.makedirs(CK, exist_ok=True)
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


QA = _load("m5_32_r8_a_quartics", "m5_32_r8_a_quartics.py")
RB, B3, B8 = QA.RB, QA.B3, QA.B8
INS4 = B8.INS4


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def cfg_of(n, L, g=8.0, delta=0.3):
    return INS4.base_cfg(s=-1.0, g=g, n=n, L=L, delta=delta)


def kin_c2(M, a0, cfg):
    """-4 x the omega^2 coefficient of I1: the measure R7 and R8 used."""
    _, c2, _, _ = QA.omega_poly(QA.d_I1, M, cfg, a0)
    return -4.0 * float(np.sum(c2) * cfg["h"] ** 3)


def kin_shells(M, a0, cfg, dr=3.0):
    _, c2, _, _ = QA.omega_poly(QA.d_I1, M, cfg, a0)
    dens = -4.0 * c2
    r = RB.boost_geom(cfg)[0]
    h3 = cfg["h"] ** 3
    edges = np.arange(0.0, cfg["L"] / 2 + 1e-9, dr)
    out = []
    for i in range(len(edges) - 1):
        m = (r >= edges[i]) & (r < edges[i + 1])
        out.append({"r_lo": float(edges[i]), "r_hi": float(edges[i + 1]),
                    "kin_shell": float(np.sum(dens[m]) * h3)})
    return out


def relax(n, L, g=8.0, maxit=3000, dt0=0.01, dt_max=0.1, plateau=None, tag=""):
    key = f"g{g:g}_n{n}_L{L:g}_it{maxit}{tag}"
    p = os.path.join(CK, f"relax_{key}.json")
    npy = os.path.join(CK, f"relax_{key}.npy")
    if os.path.exists(p) and os.path.exists(npy):
        log(f"{key}: restored")
        return json.load(open(p)), np.load(npy)
    cfg = cfg_of(n, L, g)
    M0, a0 = B8.dressed(cfg, 0.0), B8.a0_unit(cfg, 0.0)
    free = ~INS4.pin_shell(n, cfg["h"])
    kw = dict(max_iter=maxit, log_every=500, tag=f"m5_32_r10_{key}", dt0=dt0, dt_max=dt_max)
    if plateau is not None:
        kw["plateau"] = plateau
    t0 = time.time()
    M, info = INS4.fire(M0, cfg, free, **kw)
    e0, e1 = INS4.e_parts(M0, cfg), INS4.e_parts(M, cfg)
    rec = {"g": g, "n": n, "L": L, "h": cfg["h"], "maxit": maxit, "dt0": dt0, "dt_max": dt_max,
           "stop": info.get("stop"), "pinned_frac": float((~free).mean()),
           "E_u_start": float(e0[0]), "E_u_end": float(e1[0]),
           "V4_start": float(e0[1]), "V4_end": float(e1[1]),
           "kin_rigid_INS4": float(INS4.kin_of(M0, a0, cfg)),
           "kin_relaxed_INS4": float(INS4.kin_of(M, a0, cfg)),
           "kin_rigid_C2": kin_c2(M0, a0, cfg), "kin_relaxed_C2": kin_c2(M, a0, cfg),
           "rel_move": float(np.sqrt(np.sum((M - M0) ** 2)) /
                             max(np.sqrt(np.sum((M0 - INS4.vac4(cfg)) ** 2)), 1e-300)),
           "runtime_s": round(time.time() - t0, 1)}
    np.save(npy, M)
    json.dump(rec, open(p, "w"), indent=1)
    log(f"{key}: stop {rec['stop']} in {rec['runtime_s']}s | E_u {rec['E_u_start']:.4f} -> "
        f"{rec['E_u_end']:.4f} | V4 {rec['V4_end']:.5f} | kin {rec['kin_rigid_INS4']:.4f} -> "
        f"{rec['kin_relaxed_INS4']:.4f}")
    return rec, M


def stage_reproduce():
    rec, M = relax(32, 48.0)
    ref = {"E_u_end": 13.540076362700246, "V4_end": 0.25966593821613715,
           "kin_relaxed": 351.1698512170}
    dev = {k: abs(rec[{"E_u_end": "E_u_end", "V4_end": "V4_end",
                       "kin_relaxed": "kin_relaxed_INS4"}[k]] / v - 1.0)
           for k, v in ref.items()}
    log(f"reproduction of the R9 audit endpoint: relative deviations {dev}")
    return {"record": rec, "audit_reference": ref, "relative_deviation": dev,
            "reproduced": bool(max(dev.values()) < 1e-3)}


def stage_shell():
    rec, M = relax(32, 48.0)
    cfg = cfg_of(32, 48.0)
    M0, a0 = B8.dressed(cfg, 0.0), B8.a0_unit(cfg, 0.0)
    sr, sx = kin_shells(M0, a0, cfg), kin_shells(M, a0, cfg)
    rows = []
    for a, b in zip(sr, sx):
        d = b["kin_shell"] - a["kin_shell"]
        rows.append({"r_lo": a["r_lo"], "r_hi": a["r_hi"], "rigid": a["kin_shell"],
                     "relaxed": b["kin_shell"], "delta": d,
                     "rel": d / a["kin_shell"] if a["kin_shell"] else None})
    tot = sum(r["delta"] for r in rows)
    inner = sum(r["delta"] for r in rows if r["r_hi"] <= 6.0)
    log("shell test (kin per 3-unit shell, rigid -> relaxed):")
    for r in rows:
        log(f"   r {r['r_lo']:5.1f}-{r['r_hi']:5.1f}: {r['rigid']:9.3f} -> {r['relaxed']:9.3f} "
            f"({r['rel'] * 100 if r['rel'] is not None else float('nan'):+7.2f} %)")
    log(f"total change {tot:.3f}; inside r = 6: {inner:.3f} ({100 * inner / tot:.1f} %)")
    return {"shells": rows, "total_change": tot, "inside_r6": inner,
            "fraction_inside_r6": inner / tot if tot else None}


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "reproduce"
    out = {}
    if mode in ("reproduce", "all"):
        out["reproduce"] = stage_reproduce()
    if mode in ("shell", "all"):
        out["shell"] = stage_shell()
    if mode == "ladder":
        for n, L in ((48, 72.0),):
            for it in (1500, 3000, 6000):
                relax(n, L, maxit=it)
    if mode == "gcore":
        relax(32, 48.0, g=32.0, maxit=6000, dt0=1.5e-4, dt_max=1.5e-3,
              plateau=(10 ** 9, 0.0), tag="_g32")
    if out:
        p = os.path.join(DATA, f"m5_32_r10_{mode}.json")
        json.dump(out, open(p, "w"), indent=1)
        log(f"written {p}")
