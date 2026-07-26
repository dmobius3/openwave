"""M5.21.8 P0/P1: independent verification of the author's analytical
twisting massive hedgehog (theory/duda_2026-07-19_3p1D_analytical_
twisting_massive_hedgehog.pdf).

INDEPENDENT PIPELINE (no transcription of his intermediate blobs):
closed-form group elements ->
    Qb = exp(m r-hat . boosts)        I + sinh(m) K + (cosh(m)-1) K^2
    Qh = Rz(phi) Ry(theta) Rx(omega t) per-generator Rodrigues
    M  = Qb Qh d Qh^T Qb^T,  d = diag(g, 1, delta, 0)
derivatives by high-accuracy central differences (Richardson h, h/2),
F_ij = dM_i xi dM_j - dM_j xi dM_i  (xi = diag(-1,1,1,1)),
his Hamiltonian split per entry:
    H = sum_ij [F(1,2)^2+F(1,3)^2+F(2,3)^2 - F(0,1)^2-F(0,2)^2-F(0,3)^2]
(0-based; time row negative), evaluated on the y = 0 plane
x = r cos(TH), z = r sin(TH); his normalizations replicated:
    Ha  = time average * 2  (his Integrate * omega/Pi)
    E(Delta) = 2 Int_0^{pi/2-Delta} Ha cos(TH) r^2 dTH   (cone removed)

CHECKS (pre-registered in the task PLAN):
    C1 m*: numeric minimum over m of the small-delta energy vs his
        m* = +/- (1/2) ln((1+g)/(g-1)), several g, cone-robust
    C2 Hm: d E / d(omega^2) from the pipeline vs the transcribed Hm
        formula (the dropped cone term is omega-independent)
    C3 omega-boundedness: sign of dE/d(omega^2) at (g=8, delta=0.3)
        and the (g, delta) region map of the transcribed coefficient
    C4 cone divergence: E(Delta) ~ log as Delta -> 0, coefficient
        scaling ~ delta^2
    C5 constant-omega r-divergence: r^2 * density -> const * omega^2
    C6 the finite branch: min_omega Hm vs his closed form where the
        omega^2 coefficient is positive
    P1 bridge: per-(i,j) his component-split vs our eta-trace
        tr(eta F eta F^T) on the same F

Out: data/m5_21_8_p0.json (verdict rows) + printed tables.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
XI = np.diag([-1.0, 1.0, 1.0, 1.0])

# rotation generators (his Gamma_mu layout, rows/cols = t,x,y,z)
G1 = np.zeros((4, 4)); G1[2, 3] = -1.0; G1[3, 2] = 1.0
G2 = np.zeros((4, 4)); G2[1, 3] = 1.0; G2[3, 1] = -1.0
G3 = np.zeros((4, 4)); G3[1, 2] = -1.0; G3[2, 1] = 1.0


def rot(G, a):
    return np.eye(4) + np.sin(a) * G + (1.0 - np.cos(a)) * (G @ G)


def boost_hedge(m, x, y, z):
    R = np.sqrt(x * x + y * y + z * z)
    n = np.array([x, y, z]) / R
    K = np.zeros((4, 4))
    K[0, 1:] = n
    K[1:, 0] = n
    return np.eye(4) + np.sinh(m) * K + (np.cosh(m) - 1.0) * (K @ K)


def Mfield(t, x, y, z, g, dl, m, om):
    """his ansatz, closed form."""
    phi = np.arctan2(y, x)
    th = -np.arctan2(z, np.sqrt(x * x + y * y))
    Qh = rot(G3, phi) @ rot(G2, th) @ rot(G1, om * t)
    Qb = boost_hedge(m, x, y, z)
    Q = Qb @ Qh
    d = np.diag([g, 1.0, dl, 0.0])
    return Q @ d @ Q.T


def dM_all(t, x, y, z, g, dl, m, om, h=1e-5):
    """central differences with one Richardson step (h, h/2)."""
    out = []
    for k, v in enumerate((t, x, y, z)):
        args = [t, x, y, z]

        def d_at(hh):
            a_p = list(args); a_p[k] = v + hh
            a_m = list(args); a_m[k] = v - hh
            return (Mfield(*a_p, g, dl, m, om)
                    - Mfield(*a_m, g, dl, m, om)) / (2 * hh)
        d1, d2 = d_at(h), d_at(h / 2)
        out.append((4.0 * d2 - d1) / 3.0)
    return out


def H_split(t, r, TH, g, dl, m, om):
    """his Hamiltonian density on the y=0 plane + the per-block
    bridge pieces."""
    x, z = r * np.cos(TH), r * np.sin(TH)
    dM = dM_all(t, x, 1e-12, z, g, dl, m, om)
    Hs = Ht = 0.0
    his_ss = his_tt = eta_ss = eta_tt = 0.0
    for i in range(4):
        for j in range(4):
            F = dM[i] @ XI @ dM[j] - dM[j] @ XI @ dM[i]
            sp = F[1, 2] ** 2 + F[1, 3] ** 2 + F[2, 3] ** 2
            tm = F[0, 1] ** 2 + F[0, 2] ** 2 + F[0, 3] ** 2
            Hs += sp
            Ht += tm
            et = np.trace(XI @ F @ XI @ F.T)
            if i > 0 and j > 0:
                his_ss += sp - tm
                eta_ss += et
            elif i != j:
                his_tt += sp - tm
                eta_tt += et
    return {"H": Hs - Ht,
            "his_ss": his_ss, "his_tt": his_tt,
            "eta_ss": eta_ss, "eta_tt": eta_tt}


def Ha_avg(r, TH, g, dl, m, om, nt=24):
    """his Ha = (omega/pi) * integral over one period = 2 * mean."""
    ts = (np.arange(nt) + 0.5) * (2 * np.pi / max(om, 1e-12)) / nt
    if om == 0.0:
        ts = np.array([0.0])
    return 2.0 * np.mean([H_split(t, r, TH, g, dl, m, om)["H"]
                          for t in ts])


def E_cone(g, dl, m, om, r=1.0, Delta=0.15, nTH=48):
    """2 Int_0^{pi/2-Delta} Ha cos(TH) r^2 dTH (Gauss-Legendre)."""
    a, b = 0.0, np.pi / 2 - Delta
    xs, ws = np.polynomial.legendre.leggauss(nTH)
    TH = 0.5 * (b - a) * xs + 0.5 * (b + a)
    W = 0.5 * (b - a) * ws
    vals = [Ha_avg(r, th, g, dl, m, om) * np.cos(th) * r * r
            for th in TH]
    return 2.0 * float(np.sum(W * np.array(vals)))


# ---- his transcribed downstream formulas (for cross-check only) ----
def Hm_formula(g, dl, r, om):
    """his Out[54] (at m = log((-1+g)/(1+g))/2)."""
    w = -1.0 + 2.0 * r * r * om * om
    inner = (-2.0 * g * (-1.0 + dl) * dl * w + dl * dl * w
             + g * g * (2.0 + (-2.0 + dl) * dl * w))
    return 4.0 * dl * dl * inner / ((-1.0 + g) ** 2 * r * r)


def om2_coeff(g, dl):
    """coefficient of 2 r^2 om^2 in Hm's inner bracket * 4 dl^2/(g-1)^2:
    sign decides omega-boundedness."""
    T = (-2.0 * g * (-1.0 + dl) * dl + dl * dl
         + g * g * (-2.0 + dl) * dl)
    return 8.0 * dl * dl * T / ((-1.0 + g) ** 2)


def fin_formula(g, dl, r0):
    """his Out[51] finite-branch value."""
    return ((8.0 * g * g * dl * dl
             + 8.0 * (-1.0 + g) * g * dl ** 3
             - 4.0 * (-1.0 + g) ** 2 * dl ** 4)
            / ((-1.0 + g) ** 2 * r0))


def mstar(g):
    return 0.5 * np.log((1.0 + g) / (g - 1.0))


# ================= checks =================
def c1():
    rows = []
    for g in (3.0, 8.0, 50.0):
        for Delta in (0.1, 0.2):
            ms = np.linspace(-0.5, 0.5, 41)
            dl = 1e-4
            Es = [E_cone(g, dl, m, 0.3, Delta=Delta) for m in ms]
            i = int(np.argmin(Es))
            m_fit = ms[i]
            if 0 < i < len(ms) - 1:
                a, b, c = Es[i - 1], Es[i], Es[i + 1]
                m_fit = ms[i] - 0.5 * (ms[1] - ms[0]) * (c - a) / \
                    (c - 2 * b + a)
            rows.append({"g": g, "Delta": Delta,
                         "m_min_num": float(m_fit),
                         "m_star_abs": float(abs(mstar(g))),
                         "match_abs": bool(abs(abs(m_fit)
                                               - abs(mstar(g)))
                                           < 0.01)})
            print(rows[-1], flush=True)
    return rows


def c2c3():
    rows = []
    for (g, dl) in ((8.0, 0.3), (8.0, 0.05), (3.0, 0.3), (0.5, 0.3),
                    (100.0, 1e-3)):
        m = -abs(mstar(g)) if g > 1 else 0.5 * np.log((1 + g)
                                                      / (1 - g + 1e-99))
        if g < 1:
            m = 0.0  # his m* branch is complex for g<1; formula probe only
        r = 1.3
        om1, om2 = 0.2, 0.4
        dE = (E_cone(g, dl, m, om2, r=r) - E_cone(g, dl, m, om1, r=r)) \
            / (om2 ** 2 - om1 ** 2)
        dHm = (Hm_formula(g, dl, r, om2) - Hm_formula(g, dl, r, om1)) \
            / (om2 ** 2 - om1 ** 2) * r * r * 2.0
        # pipeline E integrates Ha*cos*r^2 over TH; Hm is the density
        # AFTER his angular reduction: compare signs + report both
        rows.append({"g": g, "dl": dl, "m_used": float(m),
                     "dE_dom2_pipeline": float(dE),
                     "dHm_dom2_x2r2": float(dHm),
                     "om2_coeff_formula": float(om2_coeff(g, dl)),
                     "same_sign": bool(np.sign(dE)
                                       == np.sign(om2_coeff(g, dl)))})
        print(rows[-1], flush=True)
    return rows


def c3map():
    gs = np.array([0.2, 0.5, 0.9, 1.5, 3.0, 8.0, 30.0, 1e3, 1e10])
    dls = np.array([1e-10, 1e-3, 0.05, 0.3, 1.0, 1.9, 2.5, 3.0])
    grid = [[float(np.sign(om2_coeff(g, d))) for d in dls]
            for g in gs]
    print("om2-coeff sign map (rows g, cols delta):")
    for g, row in zip(gs, grid):
        print(f"  g={g:8g}: {row}")
    return {"gs": gs.tolist(), "dls": dls.tolist(), "sign": grid}


def c4():
    rows = []
    for dl in (0.3, 0.15):
        Ds = (0.4, 0.2, 0.1, 0.05)
        Es = [E_cone(8.0, dl, -abs(mstar(8.0)), 0.3, Delta=D)
              for D in Ds]
        d32 = (Es[2] - Es[1]) / (np.log(Ds[1]) - np.log(Ds[2]))
        d43 = (Es[3] - Es[2]) / (np.log(Ds[2]) - np.log(Ds[3]))
        rows.append({"dl": dl, "E_at_Delta": dict(zip(map(str, Ds),
                                                      map(float, Es))),
                     "log_slope_mid": float(d32),
                     "log_slope_last": float(d43)})
        print(rows[-1], flush=True)
    ratio = rows[1]["log_slope_last"] / max(rows[0]["log_slope_last"],
                                            1e-300)
    rows.append({"slope_ratio_dl_half": float(ratio),
                 "delta2_scaling_expect": 0.25})
    print(rows[-1], flush=True)
    return rows


def c5():
    g, dl, om = 8.0, 0.3, 0.3
    m = -abs(mstar(g))
    rows = []
    for r in (2.0, 4.0, 8.0, 16.0):
        v = E_cone(g, dl, m, om, r=r) - E_cone(g, dl, m, 0.0, r=r)
        rows.append({"r": r, "E_omega_part": float(v)})
        print(rows[-1], flush=True)
    return rows


def c6():
    rows = []
    for (g, dl) in ((8.0, 2.5), (3.0, 2.5), (0.5, 0.3)):
        if om2_coeff(g, dl) <= 0:
            rows.append({"g": g, "dl": dl, "skip": "coeff negative"})
            continue
        oms = np.linspace(0.0, 2.0, 81)
        vals = [Hm_formula(g, dl, 1.0, om) for om in oms]
        rows.append({"g": g, "dl": dl,
                     "min_num": float(np.min(vals)),
                     "at_om": float(oms[int(np.argmin(vals))]),
                     "his_fin_formula": float(fin_formula(g, dl, 1.0)),
                     "hm_at_0": float(Hm_formula(g, dl, 1.0, 0.0))})
        print(rows[-1], flush=True)
    return rows


def p1():
    rows = []
    rng = np.random.default_rng(218)
    for _ in range(6):
        g, dl = 8.0, 0.3
        m, om = rng.uniform(-0.3, 0.3), rng.uniform(0.05, 0.6)
        t = rng.uniform(0, 5.0)
        r, TH = rng.uniform(0.8, 3.0), rng.uniform(0.1, 1.2)
        s = H_split(t, r, TH, g, dl, m, om)
        rows.append({k: float(v) for k, v in s.items()}
                    | {"m": float(m), "om": float(om), "r": float(r)})
        print({k: round(v, 6) if isinstance(v, float) else v
               for k, v in rows[-1].items()}, flush=True)
    return rows


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    out = {}
    if mode in ("all", "c1"):
        out["C1_mstar"] = c1()
    if mode in ("all", "c2"):
        out["C2_C3_omega_slope"] = c2c3()
        out["C3_region_map"] = c3map()
    if mode in ("all", "c4"):
        out["C4_cone_divergence"] = c4()
    if mode in ("all", "c5"):
        out["C5_r_divergence"] = c5()
    if mode in ("all", "c6"):
        out["C6_finite_branch"] = c6()
    if mode in ("all", "p1"):
        out["P1_bridge_samples"] = p1()
    os.makedirs(DATA, exist_ok=True)
    tag = "p0" if mode == "all" else f"p0_{mode}"
    with open(os.path.join(DATA, f"m5_21_8_{tag}.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("saved", tag)
