"""M5.32 R22-0 (c): the flatness of the certified potential on the
degenerate-pair vacuum (the working note, Zenodo 22788604, section 8.6).

EQUATIONS FIRST
---------------
N = M eta, vacuum N-spectrum q = (-g, 1, delta, delta) (code branch s = -1).
    V4    = w sum_{p=1..4} (tr N^p - C_p)^2,   C_p = sum_i q_i^p
    Vspec = gamma tr[P(N)^2],                  P(x) = prod_i (x - q_i)
Claim under test: along the split of the degenerate pair,
(delta, delta) -> (delta + eps, delta - eps), V4 = A(delta) eps^4 + O(eps^6)
with no quadratic term, so a biaxial core costs no quadratic potential.
Derivation: tr N^p shifts by (delta + eps)^p + (delta - eps)^p - 2 delta^p
= p (p - 1) delta^(p - 2) eps^2 + O(eps^4), hence
    A(delta) = w sum_p [p (p - 1) delta^(p - 2)]^2 = w (4 + 36 delta^2 + 144 delta^4).
Audit note (m5_32_r22_1_audit.py, Q5.2): this A is the FROZEN-direction
coefficient (the top eigenvalue, the pair mean and N_00 held). With those three
massive directions relaxed along the split the valley-floor coefficient is 0.93
at delta 0.3 (the shifts per eps^2: top -0.70, pair mean +0.157).

Checks:
    A  sympy: the series of V4 along the split (the eps^2 coefficient must
       vanish, the eps^4 coefficient must equal A(delta)); the same series
       along three control directions that must NOT be flat: the common
       shift of the pair (delta + eps, delta + eps), the shift of the
       eigenvalue 1, and the shift of g. Negative control built in.
    B  the full 10 x 10 Hessian of V4 and of Vspec on the vacuum (finite
       differences of the production gradient on a uniform field): the
       count of zero modes. Expected for V4: 5 broken generators of the
       frame (the 6 of so(1,3) minus the unbroken rotation of the
       degenerate pair) plus the 2 split directions of the pair
       (eps on the diagonal and its off-diagonal partner) = 7 of 10.
    C  the lattice scan: the production v4_energy_grad on a uniform field
       at eps = 0.2 ... 0.0125, the fitted power and the coefficient.
    D  Vspec along the same directions: a double root of P at delta makes
       every shift of the pair quartic, the common shift included.

Regenerate: python m5_32_r22_0c_v4_flatness.py   (seconds)
"""

import importlib.util
import json
import os
import sys
import time

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r22_0c_v4_flatness.json")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


R20 = _load("m5_32_r20_1_axes", "m5_32_r20_1_axes.py")
B3, R0 = R20.B3, R20.R0
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
G, DELTA = 8.0, 0.3
ZERO_TOL = (
    1e-3  # the finite-difference noise of the Hessian is about 1e-5 (eigenvalues run to 8e6)
)
T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def check_a():
    g, d, e = sp.symbols("g delta epsilon", real=True)
    q = [-g, 1, d, d]
    cp = [sum(qi**p for qi in q) for p in range(1, 5)]
    dirs = {
        "split (d+e, d-e)": [-g, 1, d + e, d - e],
        "common (d+e, d+e)": [-g, 1, d + e, d + e],
        "axis-1 shift": [-g, 1 + e, d, d],
        "g shift": [-g + e, 1, d, d],
    }
    out, sym = {}, {}
    for name, lam in dirs.items():
        V = sum((sum(li**p for li in lam) - cp[p - 1]) ** 2 for p in range(1, 5))
        ser = sp.series(sp.expand(V), e, 0, 5).removeO()
        c2 = sp.simplify(ser.coeff(e, 2))
        c4 = sp.simplify(ser.coeff(e, 4))
        sym[name] = (c2, c4)
        out[name] = {
            "eps2": str(sp.factor(c2)),
            "eps4": str(sp.factor(c4)),
            "eps2_at_g8_d0.3": float(c2.subs({g: 8, d: sp.Rational(3, 10)})),
            "eps4_at_g8_d0.3": float(c4.subs({g: 8, d: sp.Rational(3, 10)})),
        }
        log(
            f"A {name:20s} eps^2: {out[name]['eps2_at_g8_d0.3']:.6g}   eps^4: {out[name]['eps4_at_g8_d0.3']:.6g}"
        )
    A_closed = 4 + 36 * d**2 + 144 * d**4
    ok = (
        sp.simplify(sym["split (d+e, d-e)"][0]) == 0
        and sp.simplify(sym["split (d+e, d-e)"][1] - A_closed) == 0
        and all(out[k]["eps2_at_g8_d0.3"] > 1e-6 for k in out if not k.startswith("split"))
    )
    out["A_delta_over_w"] = str(A_closed)
    out["A_at_d0.3_over_w"] = float(A_closed.subs(d, sp.Rational(3, 10)))
    out["PASS"] = bool(ok)
    return out


def uniform_grad(fn, M0):
    """gradient per unit volume of a potential on the uniform field M0 (a 3^3 block, h = 1)."""
    cfg = B3.base_cfg(s=-1.0, g=G, n=3, L=3.0, delta=DELTA)
    M = np.broadcast_to(M0, (3, 3, 3, 4, 4)).copy()
    return fn(M, cfg)[1][1, 1, 1]


def hessian(fn, M0, t=1e-4):
    basis = R0.sym_basis()
    H = np.zeros((len(basis), len(basis)))
    for b, Eb in enumerate(basis):
        gp = uniform_grad(fn, M0 + t * Eb)
        gm = uniform_grad(fn, M0 - t * Eb)
        dG = (gp - gm) / (2 * t)
        for a, Ea in enumerate(basis):
            H[a, b] = np.sum(dG * Ea)
    return 0.5 * (H + H.T)


def check_b():
    cfg = B3.base_cfg(s=-1.0, g=G, n=3, L=3.0, delta=DELTA)
    out = {}
    for vac, deg in (("degenerate (g,1,d,d)", True), ("biaxial (g,1,d,0)", False)):
        q = R0.roots_of(cfg, degenerate=deg)
        # N = M eta has spectrum q, M diagonal: M_00 = -q_0 (eta_00 = -1), M_ii = q_i
        M0 = np.diag([-q[0], q[1], q[2], q[3]])
        for pot, fn in (
            ("V4", lambda M, c, q=q: R0.v4_energy_grad(M, c, q, 1.0)),
            ("Vspec", lambda M, c, q=q: R0.vspec_energy_grad(M, c, q, 1.0)),
        ):
            ev = np.linalg.eigvalsh(hessian(fn, M0))
            scale = np.abs(ev).max()
            zero = int(np.sum(np.abs(ev) < ZERO_TOL))
            out[f"{pot} on {vac}"] = {
                "eigs": [float(x) for x in ev],
                "zero_modes": zero,
                "negative": int(np.sum(ev < -ZERO_TOL)),
            }
            log(
                f"B {pot:6s} on {vac:22s}: zero modes {zero}/10, negative {out[f'{pot} on {vac}']['negative']}, "
                f"nonzero {[round(float(x), 3) for x in ev if abs(x) >= ZERO_TOL]}"
            )
    return out


def check_c():
    cfg = B3.base_cfg(s=-1.0, g=G, n=3, L=3.0, delta=DELTA)
    q = R0.roots_of(cfg, degenerate=True)
    eps = np.array([0.2, 0.1, 0.05, 0.025, 0.0125])
    rows = {}
    for pot, fn in (
        ("V4", lambda M: R0.v4_energy_grad(M, cfg, q, 1.0, need_grad=False)[0]),
        ("Vspec", lambda M: R0.vspec_energy_grad(M, cfg, q, 1.0, need_grad=False)[0]),
    ):
        for name, dvec in (("split", (0, 0, 1, -1)), ("common", (0, 0, 1, 1))):
            v = []
            for e in eps:
                M0 = np.diag([-q[0], q[1] + e * dvec[1], q[2] + e * dvec[2], q[3] + e * dvec[3]])
                M = np.broadcast_to(M0, (3, 3, 3, 4, 4)).copy()
                v.append(fn(M) / (27 * cfg["h"] ** 3))
            v = np.array(v)
            power = float(np.polyfit(np.log(eps[-3:]), np.log(v[-3:]), 1)[0])
            rows[f"{pot} {name}"] = {
                "eps": eps.tolist(),
                "V_per_volume": v.tolist(),
                "power": power,
                "coef_smallest_eps": float(v[-1] / eps[-1] ** round(power)),
            }
            log(
                f"C {pot:6s} {name:7s}: power {power:.3f}, coefficient {rows[f'{pot} {name}']['coef_smallest_eps']:.5g}"
            )
    A = 4 + 36 * DELTA**2 + 144 * DELTA**4
    ok = (
        abs(rows["V4 split"]["power"] - 4) < 0.05
        and abs(rows["V4 split"]["coef_smallest_eps"] / A - 1) < 0.01
        and abs(rows["V4 common"]["power"] - 2) < 0.05
    )
    return {"rows": rows, "A_closed_over_w": A, "PASS": bool(ok)}


def main():
    res = {
        "a": check_a(),
        "b": check_b(),
        "c": check_c(),
        "W1": float(B3.W1),
        "A_times_W1": float(B3.W1 * (4 + 36 * DELTA**2 + 144 * DELTA**4)),
        "A_times_25W1": float(25 * B3.W1 * (4 + 36 * DELTA**2 + 144 * DELTA**4)),
    }
    res["runtime_s"] = time.time() - T0
    with open(OUT_JSON, "w") as f:
        json.dump(res, f, indent=1)
    log(f"PASS A {res['a']['PASS']}  C {res['c']['PASS']} -> {os.path.basename(OUT_JSON)}")


if __name__ == "__main__":
    main()
