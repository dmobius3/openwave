"""M5.21.1 phase P0: the (-g)^p spec correction, BOTH signs of g, gated.

THE CORRECTION (Duda 2026-07-15, tasks/m5_20_convo.md): the potential
targets read C_p = (-g)^p + 1 + delta^p with the SIGN of g left open,
where the current verified build (m5_20_2_a_eom.c4_of) is
C_p = g^p + 1 + delta^p. Both readings are reachable through the
existing spec stack by the SIGNED g value
    sg = s * G_T,  s in {+1, -1}:
    C_p(s)   = sg^p + 1 + delta^p        (c4_of(delta, g = sg))
    M_vac(s) = diag(-sg, 1, delta, 0)    (eta M_vac spectrum
                                          (sg, 1, delta, 0))
so s = +1 IS the verified build and s = -1 IS the corrected reading.

THE PRE-REGISTERED ANALYTIC CLAIM (measured below, not assumed):
u_eta never sees g. For exactly time-BLOCK-DIAGONAL fields (M_0i = 0)
the time row enters V4 only through mu = -M_00 (the eta-spectrum time
eigenvalue), and
    sum_p (mu^p - sg^p)^2   under   (mu, s) -> (-mu, -s)
is EXACT (the (-1)^p squares away), so the statics of block-diagonal
seeds is sign-MIRROR-EQUIVALENT: identical spatial observables up to
roundoff-seeded divergence (block-diagonality is preserved exactly by
the gradient flow: products/commutators of block-diagonal matrices
stay block-diagonal). The sign can only bite through TIME-MIXING
(M_0i != 0) textures beyond quadratic order: the 4D sector (P2), not
the statics (P1).

OUTCOME (2026-07-16, the run's own SP4 gate): the claim above is
REFUTED. The residual is r_p = a_p + s_p with a_p = the time-row
deviation term and s_p = the SPATIAL trace deviation; under the
mirror map a_p -> (-1)^p a_p while s_p is unchanged, so the ODD-p
CROSS TERMS 2 a_p s_p do NOT square away wherever both sectors
deviate simultaneously (the defect core). Measured: same-state
landscape gap E_-(mirror) - E_+ = +0.1296 (0.94%); matched-snapshot
E gap 1.1% (hedgehog) / 2.8% (loop); the anchors hold (q diff
<= 1e-4, r_half equal, core lams <= 2e-3); block-diagonality exactly
preserved (max offblock 0.0). The sign is a GENUINE statics-level
knob at the ~1% E level; s = +1 sits LOWER on both backgrounds.
SP5: the boost-texture cross-sign difference onsets at QUARTIC order
(fitted power 4.00 vacuum / 3.98 hedgehog).

GATES (pre-registered)
    SP1  target pinning per sign: V4(branch rep) = 0 exactly
         (wscale = 1) for all 4 branch reps x s = +-1 x delta grid
    SP2  branch census per sign: the timelike-eigenvector-carried
         eta-eigenvalue is an orbit invariant (checked under random
         eta-orthogonal conjugation; V4 stays 0 along the orbit);
         4 distinct labels => 4 disjoint vacuum branches per sign
         (the M5.18 structure, re-enumerated per sign)
    SP3  gap map per sign: analytic signed eigenvalue-sector Hessian
         2 w Jt^T Jt, Jt_pi = p lam_i^{p-1} eta_ii on
         lam = (sg, 1, delta, 0), 6 exact boost/rotation flats;
         FD cross-check <= 1e-6 rel; ladders REPORTED per sign
    SP4  statics regression per sign (the anchors): FIRE relax of the
         M5.21-B electron hedgehog seed (+ a short loop-seed relax)
         under s = +1 vs s = -1: E within 1e-6 rel, q within 5e-3,
         core spectrum within 1e-4 abs, r_half equal bins
    SP5  where the sign DOES bite: boost-texture energy curves
         E_s(eps) on the vacuum field and on the relaxed hedgehog;
         the cross-sign difference D(eps) = dE_+(eps) - dE_-(eps)
         fitted for its leading power in eps: quantifies "the
         correction is a 4D-sector question, not a statics one"

Run:  python m5_21_1_a_spec.py gates | regress [iters] | all [iters]
Out:  ../data/m5_21_1_a_spec.json,
      ../data/m5_21_1_a_relaxed_sp.npz / _sm.npz (P1 warm starts),
      ../plots/m5_21_1_a_spec.png
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ARGV = sys.argv[1:]
sys.argv = sys.argv[:1]

from scipy.linalg import expm                                       # noqa: E402

from m5_17_energy import cell_weights, grid_coords                  # noqa: E402
from m5_20_2_a_eom import (ETA, G_T, H, NR, NZ, WSCALE, c4_of,      # noqa: E402
                           grad_static_4, sym_basis4,
                           total_energy_4, vac4, vdens4_point)
from m5_20_1_b_seeds import loop_field_biax, winding_measure_biax   # noqa: E402
from m5_20_4_a_bvp import gen                                       # noqa: E402
from m5_21_a_snap import spectral_amplitude                         # noqa: E402
from m5_21_b_electron import (core_spectrum, fire_relax,            # noqa: E402
                              meridional_charge)

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
DELTA = 0.3
SIGNS = {"sp": +1.0, "sm": -1.0}          # sp = verified build, sm = (-g)^p
BRANCHES = ("g_timelike", "one_timelike", "delta_timelike",
            "zero_timelike")


# ================= signed spec helpers =================
def make_egf_s(delta, sg):
    """energy + gradient closures at the signed g (the make_egf_4
    pattern with the g knob exposed)."""
    w4 = cell_weights(NR, NZ, H)[..., None, None]
    rho4 = ((np.arange(NR - 1) + 0.5) * H)[:, None, None, None]

    def e_fn(MM):
        return total_energy_4(MM, WSCALE, delta, g=sg)

    def g_fn(MM):
        return grad_static_4(MM, WSCALE, delta, g=sg, w4=w4, rho4=rho4)
    return e_fn, g_fn


def electron_seed_s(delta, sg, r_c=4.0):
    """the M5.21-B biaxial-hedgehog + 3-equal-core seed with the SIGNED
    time row M_00 = -sg (the branch representative per sign)."""
    R, Z = grid_coords(NR, NZ, H)
    r = np.sqrt(R ** 2 + Z ** 2)
    rs = np.where(r < 1e-12, 1e-12, r)
    n = np.stack([R / rs, np.zeros_like(R), Z / rs], axis=-1)
    m = np.broadcast_to(np.array([0.0, 1.0, 0.0]), n.shape)
    S = (n[..., :, None] * n[..., None, :]
         + delta * m[..., :, None] * m[..., None, :])
    a = (1.0 + delta) / 3.0
    core = a * np.eye(3)
    w = (1.0 - np.exp(-((r / r_c) ** 2)))[..., None, None]
    M = np.zeros((NR, NZ, 4, 4))
    M[..., 1:4, 1:4] = w * S + (1.0 - w) * core
    M[..., 0, 0] = -sg
    return M


def loop_seed_s(delta, sg, pairing="pair_1d", R0=17.0, q=0.5):
    R, Z = grid_coords(NR, NZ, H)
    M = loop_field_biax(R, Z, R0, q, delta, pairing)
    M[..., 0, 0] = -sg
    return M


def r_half_s(Mnp, delta, sg):
    """the melt half-radius against the SIGNED vacuum eta-spectrum."""
    R, Z = grid_coords(NR, NZ, H)
    r = np.sqrt(R ** 2 + Z ** 2)
    A = spectral_amplitude(Mnp, (sg, 1.0, delta, 0.0))
    rb = np.arange(0.5, 40.0, 1.0)
    prof = np.array([np.median(A[(r >= x - 0.5) & (r < x + 0.5)])
                     for x in rb])
    below = np.where(prof < 0.5 * prof[0])[0]
    return float(rb[below[0]]) if below.size else float("nan")


# ================= SP1 + SP2: pinning + branch census =================
def sp1_pinning():
    out = {}
    for tag, s in SIGNS.items():
        for delta in (0.1, 0.3, 0.5):
            for br in BRANCHES:
                v = vdens4_point(vac4(delta, g=s * G_T, branch=br),
                                 delta, 1.0, g=s * G_T)
                out[f"{tag}_d{delta}_{br}"] = v
    ok = all(abs(v) < 1e-18 for v in out.values())
    return ok, out


def _timelike_label(M):
    """the eta-eigenvalue carried by the timelike eigenvector of eta M
    (the orbit-invariant branch label). Generalized eigenproblem
    M v = lam eta v; causal character = sign(v^T eta v)."""
    lam, V = np.linalg.eig(ETA @ M)
    labels = []
    for k in range(4):
        v = np.real(V[:, k])
        nrm = float(v @ ETA @ v)
        if nrm < -1e-10:
            labels.append(float(np.real(lam[k])))
    return sorted(labels)


def sp2_branch_census(seed=7):
    rng = np.random.default_rng(seed)
    out = {}
    ok = True
    for tag, s in SIGNS.items():
        sg = s * G_T
        rows = []
        for br in BRANCHES:
            M0 = vac4(DELTA, g=sg, branch=br)
            lab0 = _timelike_label(M0)
            # random eta-orthogonal conjugation: L = exp(A), A in so(1,3)
            drift, v_orbit = 0.0, 0.0
            for _ in range(3):
                A = sum(c * gen(n) for c, n in zip(
                    rng.normal(scale=0.4, size=6),
                    ("J12", "J13", "J23", "K1", "K2", "K3")))
                L = expm(A)
                Mo = L @ M0 @ L.T
                v_orbit = max(v_orbit,
                              abs(vdens4_point(Mo, DELTA, 1.0, g=sg)))
                lab = _timelike_label(Mo)
                if len(lab) == len(lab0):
                    drift = max(drift, float(np.max(np.abs(
                        np.array(lab) - np.array(lab0)))))
                else:
                    drift = np.inf
            rows.append({"branch": br, "label": lab0,
                         "label_drift_under_orbit": drift,
                         "V4_on_orbit_max": v_orbit})
        labels = [tuple(np.round(r["label"], 9)) for r in rows]
        distinct = len(set(labels)) == len(labels)
        ok = ok and distinct and all(
            r["label_drift_under_orbit"] < 1e-7
            and r["V4_on_orbit_max"] < 1e-16 for r in rows)
        out[tag] = {"rows": rows, "n_branches_distinct": len(set(labels)),
                    "disjoint": distinct}
    return ok, out


# ================= SP3: signed gap map =================
def eig_sector_s(delta, sg, wscale=WSCALE):
    """the signed analytic eigenvalue-sector Hessian 2 w Jt^T Jt on
    lam = (sg, 1, delta, 0), sgn = eta_ii (the m5_20_2_a_eom
    analytic_eig_sector with the g sign exposed)."""
    lam = np.array([sg, 1.0, delta, 0.0])
    sgn = np.array([-1.0, 1.0, 1.0, 1.0])
    J = np.array([[p * lam[i] ** (p - 1) * sgn[i] for i in range(4)]
                  for p in range(1, 5)])
    return 2.0 * wscale * (J.T @ J)


def hessian4_s(delta, sg, wscale=WSCALE, eps=1e-5):
    M0 = vac4(delta, g=sg)
    basis = sym_basis4()
    Hm = np.zeros((10, 10))
    for a in range(10):
        for b in range(10):
            def v(pa, pb):
                return vdens4_point(M0 + eps * (pa * basis[a]
                                                + pb * basis[b]),
                                    delta, wscale, g=sg)
            Hm[a, b] = (v(1, 1) - v(1, -1) - v(-1, 1) + v(-1, -1)) \
                / (4 * eps ** 2)
    return 0.5 * (Hm + Hm.T)


def sp3_gapmap():
    rows = []
    worst = 0.0
    for delta in (0.0, 0.1, 0.3, 0.5):
        row = {"delta": delta}
        for tag, s in SIGNS.items():
            sg = s * G_T
            Ha = np.zeros((10, 10))
            Ha[:4, :4] = eig_sector_s(delta, sg)
            ev = np.linalg.eigvalsh(Ha)
            Hn = hessian4_s(delta, sg)
            rel = float(np.max(np.abs(Hn - Ha))
                        / max(np.max(np.abs(Ha)), 1e-300))
            worst = max(worst, rel)
            row[tag] = {"eigs": ev.tolist(),
                        "omegas": np.sqrt(np.maximum(ev, 0.0)).tolist(),
                        "n_zero": int(np.sum(np.abs(ev) < 1e-12)),
                        "n_negative": int(np.sum(ev < -1e-12)),
                        "FD_crosscheck_rel": rel}
        omp = np.array(row["sp"]["omegas"])
        omm = np.array(row["sm"]["omegas"])
        nz = omp > 1e-9
        row["ladder_rel_diff_sm_vs_sp"] = float(np.max(
            np.abs(omm[nz] - omp[nz]) / omp[nz])) if nz.any() else 0.0
        rows.append(row)
    return worst < 1e-6, rows, worst


# ================= SP4: statics regression =================
def _hedge_obs(M, delta, sg, e_fn):
    return {"E": e_fn(M), "core": core_spectrum(M),
            "q_meridional": meridional_charge(M),
            "r_half": r_half_s(M, delta, sg)}


def sp4_regress(iters=4000, loop_iters=1200):
    out = {}
    states = {}
    for tag, s in SIGNS.items():
        sg = s * G_T
        e_fn, g_fn = make_egf_s(DELTA, sg)
        M0 = electron_seed_s(DELTA, sg)
        print(f"[SP4 {tag}] hedgehog relax ({iters} iters, sg = {sg:+.1f})",
              flush=True)
        Mr, rel = fire_relax(M0, e_fn, g_fn, max_iter=iters)
        states[tag] = Mr
        out[tag] = {"hedgehog": _hedge_obs(Mr, DELTA, sg, e_fn),
                    "fire_tail": rel["trace"][-1],
                    "block_diag_max_offblock": float(
                        np.max(np.abs(Mr[..., 0, 1:])))}
        np.savez_compressed(
            os.path.join(DATA, f"m5_21_1_a_relaxed_{tag}.npz"), M=Mr)
        # loop-seed short regression
        Ml0 = loop_seed_s(DELTA, sg)
        print(f"[SP4 {tag}] loop relax ({loop_iters} iters)", flush=True)
        Mlr, rel_l = fire_relax(Ml0, e_fn, g_fn, max_iter=loop_iters,
                                log_every=400)
        ql, _ = winding_measure_biax(Mlr, NR, NZ, H, 17.0, 0.0)
        out[tag]["loop"] = {"E": e_fn(Mlr),
                            "q_ring": None if not np.isfinite(ql)
                            else float(ql),
                            "fire_tail": rel_l["trace"][-1]}
    hp, hm = out["sp"]["hedgehog"], out["sm"]["hedgehog"]
    lp, lm = out["sp"]["loop"], out["sm"]["loop"]
    cmp_ = {
        "E_rel_diff": abs(hp["E"] - hm["E"]) / abs(hp["E"]),
        "q_abs_diff": abs(hp["q_meridional"] - hm["q_meridional"]),
        "core_lams_abs_diff": float(np.max(np.abs(
            np.array(hp["core"]["core_lams"])
            - np.array(hm["core"]["core_lams"])))),
        "r_half_equal": bool(hp["r_half"] == hm["r_half"]),
        "loop_E_rel_diff": abs(lp["E"] - lm["E"]) / abs(lp["E"]),
        "loop_q_abs_diff": None
        if lp["q_ring"] is None or lm["q_ring"] is None
        else abs(lp["q_ring"] - lm["q_ring"]),
    }
    ok = (cmp_["E_rel_diff"] < 1e-6 and cmp_["q_abs_diff"] < 5e-3
          and cmp_["core_lams_abs_diff"] < 1e-4 and cmp_["r_half_equal"]
          and cmp_["loop_E_rel_diff"] < 1e-6)
    out["compare"] = cmp_
    return ok, out, states


# ================= SP5: where the sign bites =================
def _boost_bump(shape, center, sig, eps):
    R, Z = grid_coords(NR, NZ, H)
    bump = eps * np.exp(-((R - center[0]) ** 2 + (Z - center[1]) ** 2)
                        / (2 * sig ** 2))
    B = np.zeros(shape)
    B[..., 0, 1] = bump
    B[..., 1, 0] = bump
    return B


def sp5_boost_curves(states=None,
                     eps_grid=(0.02, 0.05, 0.1, 0.2)):
    out = {}
    for bg_tag in ("vacuum", "hedgehog"):
        rows = []
        for eps in eps_grid:
            row = {"eps": eps}
            for tag, s in SIGNS.items():
                sg = s * G_T
                e_fn, _ = make_egf_s(DELTA, sg)
                if bg_tag == "vacuum":
                    M = np.zeros((NR, NZ, 4, 4))
                    M[..., :, :] = vac4(DELTA, g=sg)
                    center = (40.0, 0.0)
                elif states is not None and tag in states:
                    M = states[tag]
                    center = (2.0, 0.0)
                else:
                    row = None
                    break
                E0 = e_fn(M)
                Ep = e_fn(M + _boost_bump(M.shape, center, 3.0, eps))
                Em = e_fn(M + _boost_bump(M.shape, center, 3.0, -eps))
                row[tag] = {"dE_plus": Ep - E0, "dE_minus": Em - E0,
                            "odd_part": 0.5 * ((Ep - E0) - (Em - E0))}
            if row is None:
                break
            if "sp" in row and "sm" in row:
                row["cross_sign_diff"] = (row["sp"]["dE_plus"]
                                          - row["sm"]["dE_plus"])
            rows.append(row)
        if rows:
            eps_a = np.array([r["eps"] for r in rows])
            d_a = np.array([abs(r["cross_sign_diff"]) for r in rows])
            pw = None
            if np.all(d_a > 0):
                pw = float(np.polyfit(np.log(eps_a), np.log(d_a), 1)[0])
            out[bg_tag] = {"rows": rows, "cross_diff_leading_power": pw}
    return out


# ================= plot =================
def make_plot(res):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.0))
    rows = res["SP3"]["rows"]
    ds = [r["delta"] for r in rows]
    for tag, mk in (("sp", "o-"), ("sm", "s--")):
        om = np.array([r[tag]["omegas"] for r in rows])
        for k in range(10):
            axes[0].plot(ds, om[:, k], mk, ms=3, lw=1,
                         color=f"C{k % 10}", alpha=0.8)
    axes[0].set_xlabel("delta")
    axes[0].set_ylabel("omega")
    axes[0].set_title("vacuum mass ladder per sign\n"
                      "(circles s = +1, squares s = -1)", fontsize=9)
    for bg_tag, ax in (("vacuum", axes[1]), ("hedgehog", axes[2])):
        if bg_tag not in res["SP5"]:
            continue
        rr = res["SP5"][bg_tag]["rows"]
        eps = [r["eps"] for r in rr]
        ax.loglog(eps, [abs(r["cross_sign_diff"]) for r in rr], "o-",
                  label="|E_+ - E_-|(eps)")
        ax.loglog(eps, [abs(r["sp"]["odd_part"]) for r in rr], "s--",
                  label="odd part, s = +1")
        ax.loglog(eps, [abs(r["sm"]["odd_part"]) for r in rr], "v--",
                  label="odd part, s = -1")
        pw = res["SP5"][bg_tag]["cross_diff_leading_power"]
        ax.set_title(f"boost-texture E curves on {bg_tag}\n"
                     f"cross-sign leading power "
                     f"{'n/a' if pw is None else f'{pw:.2f}'}", fontsize=9)
        ax.set_xlabel("eps")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOTS, "m5_21_1_a_spec.png"), dpi=130)
    plt.close(fig)


def main(which="all", iters=4000):
    os.makedirs(DATA, exist_ok=True)
    res = {"task": "M5.21.1", "phase": "P0", "delta": DELTA, "G_T": G_T,
           "wscale": WSCALE, "grid": [NR, NZ, H],
           "c4_sp": list(c4_of(DELTA, G_T)),
           "c4_sm": list(c4_of(DELTA, -G_T))}
    print(f"[P0] C_p(s=+1) = {np.round(res['c4_sp'], 4)}")
    print(f"[P0] C_p(s=-1) = {np.round(res['c4_sm'], 4)}")
    states = None
    if which in ("gates", "all"):
        ok1, d1 = sp1_pinning()
        res["SP1"] = {"ok": bool(ok1),
                      "max_abs": float(max(abs(v) for v in d1.values()))}
        print(f"[SP1] {'PASS' if ok1 else 'FAIL'} max |V4(rep)| "
              f"{res['SP1']['max_abs']:.2e}")
        ok2, d2 = sp2_branch_census()
        res["SP2"] = {"ok": bool(ok2), **d2}
        for tag in SIGNS:
            labs = [r["label"] for r in d2[tag]["rows"]]
            print(f"[SP2 {tag}] {'PASS' if d2[tag]['disjoint'] else 'FAIL'}"
                  f" timelike labels per branch: {labs}")
        ok3, rows3, worst3 = sp3_gapmap()
        res["SP3"] = {"ok": bool(ok3), "rel_worst": worst3, "rows": rows3}
        print(f"[SP3] {'PASS' if ok3 else 'FAIL'} FD worst {worst3:.2e}")
        for r in rows3:
            print(f"  d = {r['delta']:.1f}: ladder rel diff sm vs sp = "
                  f"{r['ladder_rel_diff_sm_vs_sp']:.3e} "
                  f"(zeros {r['sp']['n_zero']}/{r['sm']['n_zero']}, "
                  f"neg {r['sp']['n_negative']}/{r['sm']['n_negative']})")
    if which in ("regress", "all"):
        ok4, d4, states = sp4_regress(iters=iters)
        res["SP4"] = {"ok": bool(ok4), **d4}
        print(f"[SP4] {'PASS' if ok4 else 'FAIL'} "
              + json.dumps(d4["compare"], default=float))
        res["SP5"] = sp5_boost_curves(states)
        for bg, dd in res["SP5"].items():
            print(f"[SP5 {bg}] cross-sign leading power: "
                  f"{dd['cross_diff_leading_power']} | "
                  f"|D(0.1)| = "
                  + str([abs(r['cross_sign_diff'])
                         for r in dd['rows'] if r['eps'] == 0.1]))
    if "SP3" in res and "SP5" in res:
        make_plot(res)
        print("wrote plots/m5_21_1_a_spec.png")
    with open(os.path.join(DATA, "m5_21_1_a_spec.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    print("wrote data/m5_21_1_a_spec.json")
    return res


if __name__ == "__main__":
    w = ARGV[0] if ARGV else "all"
    it = int(ARGV[1]) if len(ARGV) > 1 else 4000
    main(w, it)
