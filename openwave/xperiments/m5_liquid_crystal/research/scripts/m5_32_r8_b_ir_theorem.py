"""M5.32 R8 arm b: why the clock inertia is extensive, and what could fix it.

Three questions, each answered by a measurement or an exact statement:

 1 GENERATORS.  The clock flow is a0 = [G, M].  Its far field is [G, M_vac],
   which vanishes only if G commutes with M_vac.  Enumerate every rotation and
   boost generator of SO(1,3) acting on the internal index and report which (if
   any) annihilate the toy vacuum diag(-s g, 1, delta, 0) and what the vacuum
   spectrum would have to be for one to exist.

 2 THE TAIL.  Measure the hedgehog's far field directly: is it a CONSTANT
   vacuum (jets decaying faster than any power) or an angle-dependent point of
   the vacuum orbit (jets ~ 1/r)?  Report the decay exponent of the jets, of
   the deviation from the nearest constant vacuum, and the spectrum of M eta in
   the far field (which is orbit-invariant by the R6 theorem).

 3 THE COEFFICIENT LADDER.  For the strong clock (a free interior minimum in
   omega) the energy's omega^2 coefficient must be negative.  With the
   certified inertia extensive and every quartic's static and omega^2 parts
   IR-finite, measure the coefficient c5 a C5 term would need to flip the sign,
   per box, and its growth with L.  A coefficient that must grow with the box
   is not a coefficient.

Conventions: the run's locked table.
"""
import importlib.util, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA, CK = os.path.join(RES, "data"), os.path.join(RES, "checkpoints", "m5_32_r8")
OUT = os.path.join(DATA, "m5_32_r8_ir_theorem.json")
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


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def generators():
    """the six so(1,3) generators acting on the internal index, in the raw
    contravariant convention M -> L M L^T, so the flow is d/dt M = [X, M]_eta
    with X the generator matrix; a rotation in the (a,b) plane and a boost in
    the (0,a) plane."""
    out = {}
    for a in range(1, 4):
        for b in range(a + 1, 4):
            X = np.zeros((4, 4))
            X[a, b], X[b, a] = 1.0, -1.0
            out[f"rot_{a}{b}"] = X
    for a in range(1, 4):
        X = np.zeros((4, 4))
        X[0, a], X[a, 0] = 1.0, 1.0
        out[f"boost_{a}"] = X
    return out


def stage_generators():
    d = np.diag([-S * G, 1.0, DELTA, 0.0])
    rows = {}
    for name, X in generators().items():
        a0 = X @ d - d @ X.T            # the flow of M -> (1 + t X) M (1 + t X)^T
        rows[name] = {"a0_vac_max_abs": float(np.abs(a0).max()),
                      "annihilates_vacuum": bool(np.abs(a0).max() < 1e-14),
                      "a0_vac": a0.tolist()}
    # which vacua WOULD admit an annihilating rotation
    spec = [-S * G, 1.0, DELTA, 0.0]
    pairs = {}
    for a in range(4):
        for b in range(a + 1, 4):
            pairs[f"{a}{b}"] = {"eigenvalues": [spec[a], spec[b]],
                                "degenerate": bool(abs(spec[a] - spec[b]) < 1e-12),
                                "gap": float(abs(spec[a] - spec[b]))}
    return {"vacuum": spec, "generators": rows, "eigenvalue_pairs": pairs,
            "statement": "a rotation generator in the (a,b) plane annihilates a diagonal vacuum "
                         "if and only if its two eigenvalues are equal; the toy spectrum "
                         "(g, 1, delta, 0) is non-degenerate, so NO generator annihilates it and "
                         "the clock flow is a nonzero constant at infinity for every choice",
            "min_gap": float(min(v["gap"] for v in pairs.values()))}


def stage_tail(n=64, L=96.0):
    cfg = RB.cfg_of(n, L)
    M = B8.dressed(cfg, 0.0)
    a0 = B8.a0_unit(cfg, 0.0)
    r = RB.boost_geom(cfg)[0]
    h = cfg["h"]
    A = np.zeros((4,) + M.shape)
    for ax in range(3):
        A[1 + ax] = B3.d1(M, ax, h, "sym")
    jet = np.sqrt(sum(np.sum(A[1 + ax] ** 2, axis=(-1, -2)) for ax in range(3)))
    dvac = np.diag([-S * G, 1.0, DELTA, 0.0])
    dev = np.sqrt(np.sum((M - dvac) ** 2, axis=(-1, -2)))
    # the orbit-invariant spectrum of M eta in the far field
    Me = np.einsum("...ab,bc->...ac", M, ETA)
    edges = np.arange(6.0, L / 2 + 1e-9, 3.0)
    rows = []
    for i in range(len(edges) - 1):
        m = (r >= edges[i]) & (r < edges[i + 1])
        if not m.any():
            continue
        ev = np.sort(np.real(np.linalg.eigvals(Me[m][:: max(1, m.sum() // 200)])), axis=-1)
        rows.append({"r": 0.5 * float(edges[i] + edges[i + 1]),
                     "jet_rms": float(np.sqrt(np.mean(jet[m] ** 2))),
                     "dev_from_const_vacuum_rms": float(np.sqrt(np.mean(dev[m] ** 2))),
                     "a0_rms": float(np.sqrt(np.mean(np.sum(a0[m] ** 2, axis=(-1, -2))))),
                     "spec_mean": np.mean(ev, axis=0).tolist(),
                     "spec_std": np.std(ev, axis=0).tolist()})
    rr = np.array([x["r"] for x in rows])
    sel = rr > 0.3 * L / 2

    def ex(key):
        v = np.array([x[key] for x in rows])
        return float(np.polyfit(np.log(rr[sel]), np.log(v[sel]), 1)[0])
    return {"n": n, "L": L, "shells": rows,
            "jet_exponent": ex("jet_rms"),
            "dev_exponent": ex("dev_from_const_vacuum_rms"),
            "a0_exponent": ex("a0_rms"),
            "reading": "jet ~ r^-1 with the deviation from a CONSTANT vacuum NOT decaying is the "
                       "signature of an angle-dependent point of the vacuum orbit (the degree is "
                       "carried at infinity); a0 flat confirms the generator does not decay"}


def stage_ladder():
    """the c5 a C5 term would need to flip the energy's omega^2 sign, per box."""
    q = {t: json.load(open(os.path.join(CK, f"quart_{t}.json")))
         for t in ("n32_L48", "n48_L72", "n64_L96")}
    Ls = [48.0, 72.0, 96.0]
    out = {"note": "energy omega^2 coefficient H2 = sum_k c_k C2_k with c_I1 = -4; a C5 term with "
                   "c5 > 0 contributes c5 C2_q < 0, so the flip needs c5 >= 4 |C2_I1| / |C2_q|",
           "per_term": {}}
    for t in ("Q_I1sq", "Q_Fpair", "Q_I4sq"):
        need, w2, w4 = [], [], []
        for b, L in zip(("n32_L48", "n48_L72", "n64_L96"), Ls):
            c2i1 = q[b]["terms"]["I1"]["C2"]
            c2q = q[b]["terms"][t]["C2"]
            c4q = q[b]["terms"][t]["C4"]
            need.append(4.0 * abs(c2i1) / abs(c2q))
            w2.append(-4.0 * c2i1)          # the certified inertia (positive)
            w4.append(c4q)
        e = float(np.polyfit(np.log(Ls), np.log(need), 1)[0])
        # the omega at which the quartic term (c5 = 1) would match the certified inertia
        wmatch = [float(np.sqrt(k / (2.0 * c))) if c > 0 else None for k, c in zip(w2, w4)]
        out["per_term"][t] = {"c5_required_to_flip_H2": need, "exponent_in_L": e,
                              "certified_inertia_H2": w2, "quartic_C4": w4,
                              "omega_where_quartic_matches_at_c5_1": wmatch,
                              "radiation_window_omega_max": 0.786}
        log(f"{t}: c5 needed {need[0]:.1f} -> {need[2]:.1f} (L exponent {e:.3f}); "
            f"omega match at c5 = 1: {wmatch}")
    return out


if __name__ == "__main__":
    out = {"task": "M5.32 R8 arm b: the IR obstruction, its generators, its tail, its coefficient ladder",
           "conventions": {"vacuum": f"diag(-s g, 1, delta, 0), s = {S}, g = {G}, delta = {DELTA}",
                           "clock flow": "a0 = X M + M X^T for a generator X (the raw contravariant law)",
                           "stencil": "certified sym, h = L / n"}}
    out["generators"] = stage_generators()
    log("generators done")
    out["tail"] = stage_tail()
    log(f"tail done: jet exponent {out['tail']['jet_exponent']:.3f}, "
        f"deviation exponent {out['tail']['dev_exponent']:.3f}, a0 exponent {out['tail']['a0_exponent']:.3f}")
    out["ladder"] = stage_ladder()
    out["runtime_s"] = round(time.time() - T0, 1)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    log(f"written {OUT}")
