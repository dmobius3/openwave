"""M5.32 R9 arm a: excise the coordinate string and re-decide C5.

The R8 audit found that the author's ansatz M = Qh d4 Qh^T is discontinuous
along the whole z axis (the Euler-angle frame), and that 73 to 98 % of every C5
quartic coefficient sits in the lattice columns beside that line, so the C5
numbers of R8 have no continuum limit.  This arm removes a tube of FIXED
PHYSICAL radius (h-independent, so the excision is the same region on every
lattice) and asks the two questions that decide the class:

  Q1  with the string out, is the C5 omega^4 inertia h-STABLE (a continuum
      quantity at all) and is its L exponent below 1 (IR-convergent)?
  Q2  does the certified inertia's linear growth (R7-B7) survive the same
      excision, on both the L ladder and the h ladder?

A tube is not a physical regulator, it is a diagnostic: a quantity that depends
on the excision radius is a property of the string, and a quantity that does not
is a property of the field.
"""
import importlib.util, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA, CK = os.path.join(RES, "data"), os.path.join(RES, "checkpoints", "m5_32_r9")
os.makedirs(CK, exist_ok=True)
OUT = os.path.join(DATA, "m5_32_r9_tube.json")
T0 = time.time()


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


QA = _load("m5_32_r8_a_quartics", "m5_32_r8_a_quartics.py")
RB, B3, B8 = QA.RB, QA.B3, QA.B8
RHO_CUTS = (0.0, 1.5, 3.0, 4.5, 6.0)
BOXES = ((32, 48.0), (48, 72.0), (64, 96.0), (64, 48.0), (96, 72.0))
TERMS = ("I1", "Q_I1sq", "Q_I4sq", "Q_Fpair", "Q_C6a")


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def rho_of(cfg):
    n, h = cfg["n"], cfg["h"]
    ax = np.arange(n) * h - h * (n // 2)
    X, Y, _ = np.meshgrid(ax, ax, ax, indexing="ij")
    return np.sqrt(X * X + Y * Y)


def run_box(n, L):
    tag = f"n{n}_L{L:g}"
    p = os.path.join(CK, f"tube_{tag}.json")
    if os.path.exists(p):
        log(f"{tag}: restored")
        return json.load(open(p))
    cfg = RB.cfg_of(n, L)
    M, a0 = B8.dressed(cfg, 0.0), B8.a0_unit(cfg, 0.0)
    rho, h3 = rho_of(cfg), cfg["h"] ** 3
    out = {"n": n, "L": L, "h": cfg["h"], "rho_cuts": list(RHO_CUTS), "terms": {}}
    for t in TERMS:
        fn = QA.QUARTICS[t][1]
        d0, c2, c4, _ = QA.omega_poly(fn, M, cfg, a0)
        rec = {}
        for rc in RHO_CUTS:
            m = rho >= rc
            rec[f"rho_ge_{rc:g}"] = {"A": float(np.sum(d0[m]) * h3),
                                     "C2": float(np.sum(c2[m]) * h3),
                                     "C4": float(np.sum(c4[m]) * h3),
                                     "cells_kept_frac": float(m.mean())}
        out["terms"][t] = rec
        log(f"{tag} {t}: C4 full {rec['rho_ge_0']['C4']:.5g} -> "
            f"rho>=3 {rec['rho_ge_3']['C4']:.5g}; C2 full {rec['rho_ge_0']['C2']:.5g} -> "
            f"rho>=3 {rec['rho_ge_3']['C2']:.5g}")
    json.dump(out, open(p, "w"), indent=1)
    return out


def analyze(res):
    """L exponents at fixed h = 1.5, and h exponents at fixed L."""
    L_ladder = [("n32_L48", 48.0), ("n48_L72", 72.0), ("n64_L96", 96.0)]
    h_pairs = [("n32_L48", "n64_L48", 48.0), ("n48_L72", "n96_L72", 72.0)]
    out = {"L_exponent_at_h_1.5": {}, "h_exponent_at_fixed_L": {}}
    for t in TERMS:
        out["L_exponent_at_h_1.5"][t] = {}
        for rc in RHO_CUTS:
            k = f"rho_ge_{rc:g}"
            Ls = np.array([L for _, L in L_ladder])
            for key in ("A", "C2", "C4"):
                v = np.array([res[b]["terms"][t][k][key] for b, _ in L_ladder])
                e = (float(np.polyfit(np.log(Ls), np.log(np.abs(v)), 1)[0])
                     if np.all(v != 0) else None)
                out["L_exponent_at_h_1.5"][t].setdefault(k, {})[key] = e
        out["h_exponent_at_fixed_L"][t] = {}
        for rc in RHO_CUTS:
            k = f"rho_ge_{rc:g}"
            for key in ("A", "C2", "C4"):
                es = []
                for a, b, _ in h_pairs:
                    va, vb = res[a]["terms"][t][k][key], res[b]["terms"][t][k][key]
                    ha, hb = res[a]["h"], res[b]["h"]
                    if va != 0 and vb != 0 and np.sign(va) == np.sign(vb):
                        es.append(float(np.log(abs(vb / va)) / np.log(hb / ha)))
                out["h_exponent_at_fixed_L"][t].setdefault(k, {})[key] = (
                    float(np.mean(es)) if es else None)
    return out


if __name__ == "__main__":
    res = {}
    for n, L in BOXES:
        res[f"n{n}_L{L:g}"] = run_box(n, L)
    an = analyze(res)
    out = {"task": "M5.32 R9 arm a: the coordinate string excised, C5 re-decided",
           "excision": "cells with cylindrical radius rho < rho_cut about the z axis are dropped; "
                       "rho_cut is a FIXED PHYSICAL length, identical on every lattice",
           "boxes": res, "analysis": an,
           "reading": "an exponent that moves with rho_cut belongs to the string; one that does not "
                      "belongs to the field. h exponent near 0 = a continuum quantity. "
                      "L exponent: 0 = IR-finite, 1 = linear, 3 = volume"}
    json.dump(out, open(OUT, "w"), indent=1)
    log(f"written {OUT}")
    print()
    print(f"{'term':>9} {'cut':>5} | {'C2 L-exp':>9} {'C2 h-exp':>9} | {'C4 L-exp':>9} {'C4 h-exp':>9}")
    for t in TERMS:
        for rc in RHO_CUTS:
            k = f"rho_ge_{rc:g}"
            a = an["L_exponent_at_h_1.5"][t][k]
            b = an["h_exponent_at_fixed_L"][t][k]
            f = lambda x: "   None" if x is None else f"{x:9.3f}"
            print(f"{t:>9} {rc:5.1f} | {f(a['C2'])} {f(b['C2'])} | {f(a['C4'])} {f(b['C4'])}")
