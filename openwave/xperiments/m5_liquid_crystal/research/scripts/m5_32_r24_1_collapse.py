"""M5.32 R24-1: the collapse read of the biaxial halo (no power law in it).

The linear halo around the hedgehog obeys eps'' = A eps / r^2 + beta r^2 eps,
beta = (g + delta)(1 - delta) c / 2K, K = 8 (1 - delta)^2. Its decaying solution
is sqrt(r) K_nu(rho), rho = sqrt(beta) r^2 / 2, nu = sqrt(1 + 4A) / 4, so

    eps(r) = Amp * r^(1/2 - 2 nu) * m_nu(rho),   m_nu(rho) = rho^nu K_nu(rho) / (2^(nu-1) Gamma(nu)),

one master curve for every c, one free amplitude per row (the core's boundary
condition), no fitted exponent and no fitted radius. Two values of A are read,
both fixed before any field was opened:

    A = 1   the form of the 2026-09-20 reply (nu = sqrt 5 / 4, tail r^-0.618),
            with the master curve corrected to rho^nu K_nu (R24-0 b);
    A = 3   the R24-0 (e3) second variation on the lowest spin-2 mode
            (nu = sqrt 13 / 4, tail r^-1.303).

The read uses the stored R23 end fields (and the R24-2 threshold rows when they
exist); nothing is relaxed here. The shell profile is the mean of
eps = (lam_mid - lam_small) / 2 over cells binned by radius (bin 0.5 h), with
the shell RMS beside it.

Pre-registered window of a row: 6 <= r <= L/2 - 6, eps above 1e-3 of its r 6
value (the eigen-solver floor), and LINEAR cells only: the quartic force of V4
over the quadratic force of c, 0.0472 eps^2 / c, under 0.5. A row qualifies with
at least 6 shells in its window and a gate label (AT_GATE or BELOW_GATE).

Pre-registered labels, per A:
    COLLAPSE           at least 3 qualifying rows on at least 2 values of c, the pooled
                       RMS of log(eps / prediction) at most 0.10 with the radius scale s = 1
    COLLAPSE_SHIFTED   fails at s = 1, passes with ONE global s in [0.80, 1.25]
                       (rho -> rho / s^2), and the per-row best s spread under 10 percent
    NO_COLLAPSE        otherwise
    INSUFFICIENT       fewer than 3 qualifying rows, or one value of c
and overall: MASTER_A1, MASTER_A3 (one passes, or both pass and the pooled RMS
differ by over 0.02), BOTH, NEITHER.

Known going in (disclosed): the R23 tail fits with a free radius per row gave
RMS 0.03, 0.07, 0.15 and radii 1.09 to 1.18 of beta^(-1/4) at A = 1, so
COLLAPSE at s = 1 and A = 1 is expected to fail.

Supporting read near the threshold: Q(c) = eps(6) over the master prediction at
r 6 carrying the pooled amplitude of the collapsed rows of the same spacing.

Modes: smoke | run | plot. Output: data/m5_32_r24_1_collapse.json, plots/m5_32_r24_1_collapse.png.
"""

import importlib.util
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, kve

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_JSON = os.path.join(DATA, "m5_32_r24_1_collapse.json")
R24_2_JSON = os.path.join(DATA, "m5_32_r24_2_threshold.json")
R24_2_NPZ = os.path.join(DATA, "m5_32_r24_2")
G_, D_ = 8.0, 0.3
K_EL = 8 * (1 - D_) ** 2
M2 = (G_ + D_) * (1 - D_)
A_LIST = (1.0, 3.0)
R_IN, WALL, FLOOR, LIN_MAX, MIN_SHELLS = 6.0, 6.0, 1e-3, 0.5, 6
QUARTIC_OVER_QUADRATIC = 0.0472  # 4 * 0.137 / (2 * 5.81): the force ratio is this times eps^2 / c
RMS_MAX, S_RANGE, S_SPREAD = 0.10, (0.80, 1.25), 0.10


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def beta_of(c):
    return M2 * c / (2 * K_EL)


def nu_of(A):
    return np.sqrt(1 + 4 * A) / 4


def log_master(rho, nu):
    """log of rho^nu K_nu(rho) / (2^(nu - 1) Gamma(nu)); -> 0 as rho -> 0."""
    return nu * np.log(rho) + np.log(kve(nu, rho)) - rho - ((nu - 1) * np.log(2.0) + gammaln(nu))


def log_pred(r, c, A, s=1.0):
    nu = nu_of(A)
    rho = np.sqrt(beta_of(c)) * r * r / (2 * s * s)
    return (0.5 - 2 * nu) * np.log(r) + log_master(rho, nu)


# ================= the profile =================
def shell_profile(eps, r, h):
    """mean, RMS and count of eps over cells binned by radius, bin width h / 2."""
    k = np.floor(r.ravel() / (0.5 * h)).astype(int)
    n = np.bincount(k)
    s1 = np.bincount(k, weights=eps.ravel())
    s2 = np.bincount(k, weights=eps.ravel() ** 2)
    rs = np.bincount(k, weights=r.ravel())
    ok = n > 0
    return np.stack([rs[ok] / n[ok], s1[ok] / n[ok], np.sqrt(s2[ok] / n[ok]), n[ok]], 1)


def eps_of_field(M):
    lam = np.linalg.eigvalsh(M[..., 1:, 1:])
    return 0.5 * (lam[..., 1] - lam[..., 0])


def window(prof, c, L, col=1):
    r, e = prof[:, 0], prof[:, col]
    e6 = float(np.interp(R_IN, r, e))
    k = (r >= R_IN) & (r <= 0.5 * L - WALL) & (e > FLOOR * e6) & (e > 0)
    k &= QUARTIC_OVER_QUADRATIC * e * e / c < LIN_MAX
    return k


# ================= the fits =================
def row_resid(row, A, s, col=1):
    r, e = row["prof"][row["win"], 0], row["prof"][row["win"], col]
    d = np.log(e) - log_pred(r, row["c"], A, s)
    return d - d.mean(), float(d.mean())


def pooled_rms(rows, A, s):
    d = np.concatenate([row_resid(q, A, s)[0] for q in rows])
    return float(np.sqrt(np.mean(d * d)))


def best_s(rows, A):
    o = minimize_scalar(
        lambda ls: pooled_rms(rows, A, np.exp(ls)),
        bounds=(np.log(0.5), np.log(2.0)),
        method="bounded",
    )
    return float(np.exp(o.x)), float(o.fun)


def read_A(rows, A):
    out = {"A": A, "nu": float(nu_of(A)), "tail_exponent": float(0.5 - 2 * nu_of(A))}
    cs = sorted({q["c"] for q in rows})
    if len(rows) < 3 or len(cs) < 2:
        return dict(out, label="INSUFFICIENT", qualifying_rows=len(rows), c_values=cs)
    out["pooled_rms_s1"] = pooled_rms(rows, A, 1.0)
    s, rms = best_s(rows, A)
    out["global_s"], out["pooled_rms_global_s"] = s, rms
    per = []
    for q in rows:
        sr, rr = best_s([q], A)
        d1, amp = row_resid(q, A, 1.0)
        per.append(
            {
                "tag": q["tag"],
                "c": q["c"],
                "shells": int(q["win"].sum()),
                "rms_s1": float(np.sqrt(np.mean(d1 * d1))),
                "log_amplitude_s1": amp,
                "s_row": sr,
                "rms_s_row": rr,
            }
        )
    out["rows"] = per
    sr = np.array([p["s_row"] for p in per])
    out["s_row_spread"] = float((sr.max() - sr.min()) / sr.mean())
    lc = np.log([p["c"] for p in per])
    out["slope_log_s_row_vs_log_c"] = (
        float(np.polyfit(lc, np.log(sr), 1)[0]) if len(cs) >= 2 else None
    )
    if out["pooled_rms_s1"] <= RMS_MAX:
        out["label"] = "COLLAPSE"
    elif rms <= RMS_MAX and S_RANGE[0] <= s <= S_RANGE[1] and out["s_row_spread"] < S_SPREAD:
        out["label"] = "COLLAPSE_SHIFTED"
    else:
        out["label"] = "NO_COLLAPSE"
    return out


def free_A(rows):
    """not a label: the A that minimizes the pooled RMS, at s = 1 and at the global s."""
    res = {}
    for nm, shifted in (("s1", False), ("global_s", True)):
        f = lambda A: best_s(rows, A)[1] if shifted else pooled_rms(rows, A, 1.0)  # noqa: E731
        o = minimize_scalar(f, bounds=(-0.2, 8.0), method="bounded")
        res[nm] = {"A": float(o.x), "pooled_rms": float(o.fun)}
    return res


def overall(reads):
    ok = {q["A"]: q["label"] in ("COLLAPSE", "COLLAPSE_SHIFTED") for q in reads}
    if any(q["label"] == "INSUFFICIENT" for q in reads):
        return "INSUFFICIENT"
    if ok[1.0] and ok[3.0]:
        r1, r3 = (min(q["pooled_rms_s1"], q["pooled_rms_global_s"]) for q in reads)
        if abs(r1 - r3) <= 0.02:
            return "BOTH"
        return "MASTER_A1" if r1 < r3 else "MASTER_A3"
    if ok[1.0]:
        return "MASTER_A1"
    if ok[3.0]:
        return "MASTER_A3"
    return "NEITHER"


def peel(row, A, s, amp, col=1):
    """walking inward from the window: the last shell that stays within 0.2 in log of the prediction."""
    r, e = row["prof"][:, 0], row["prof"][:, col]
    k0 = np.where(row["win"])[0]
    if len(k0) == 0:
        return None
    i = k0[0]
    while (
        i > 0
        and e[i - 1] > 0
        and abs(np.log(e[i - 1]) - log_pred(r[i - 1], row["c"], A, s) - amp) <= 0.2
    ):
        i -= 1
    return {
        "r_peel": float(r[i]),
        "rho_peel": float(np.sqrt(beta_of(row["c"])) * r[i] ** 2 / (2 * s * s)),
    }


def analyse(rows):
    good = [q for q in rows if q["qualifies"]]
    reads = [read_A(good, A) for A in A_LIST]
    out = {"reads": reads, "overall": overall(reads)}
    if len(good) >= 3:
        out["free_A_not_a_label"] = free_A(good)
    out["rms_profile_reads"] = []
    for A in A_LIST:  # the same labels on the shell RMS instead of the shell mean
        alt = [dict(q, win=window(q["prof"], q["c"], q["L"], col=2)) for q in good]
        alt = [q for q in alt if q["win"].sum() >= MIN_SHELLS]
        if len(alt) >= 3:
            d = np.concatenate([row_resid(q, A, 1.0, col=2)[0] for q in alt])
            out["rms_profile_reads"].append(
                {"A": A, "pooled_rms_s1": float(np.sqrt(np.mean(d * d))), "rows": len(alt)}
            )
    return out


# ================= rows =================
def load_rows():
    CS = _load("m5_32_r23_1_cscan", "m5_32_r23_1_cscan.py")
    meta = [(CS.OUT_NPZ, t, q) for t, q in CS.load_all_rows().items()]
    if os.path.exists(R24_2_JSON):
        with open(R24_2_JSON) as f:
            meta += [(R24_2_NPZ, t, q) for t, q in json.load(f)["rows"].items()]
    rows = []
    for folder, tag, q in meta:
        f = os.path.join(folder, tag + ".npz")
        if (
            q.get("status") != "OK"
            or abs(q["delta"] - 0.3) > 1e-9
            or q["c"] <= 0
            or not os.path.exists(f)
        ):
            continue
        cfg = CS.R21.cfg_of(q["n"], q["L"], G_, q["delta"])
        X, Y, Z = CS.B3.coords(cfg["n"], cfg["h"])
        r = np.sqrt(X * X + Y * Y + Z * Z)
        prof = shell_profile(eps_of_field(np.load(f)["M"]), r, cfg["h"])
        last = q["chunks"][-1]
        gate = last.get("fmax_spatial") is not None and last["fmax_spatial"] < q["gate"]
        win = window(prof, q["c"], q["L"])
        rows.append(
            {
                "tag": tag,
                "c": q["c"],
                "n": q["n"],
                "L": q["L"],
                "h": cfg["h"],
                "label": q["label"],
                "at_gate": bool(gate),
                "prof": prof,
                "win": win,
                "qualifies": bool(gate and win.sum() >= MIN_SHELLS),
            }
        )
    return rows


def run():
    rows = load_rows()
    out = {
        "task": "M5.32 R24-1",
        "constants": {"R_IN": R_IN, "WALL": WALL, "FLOOR": FLOOR, "LIN_MAX": LIN_MAX},
    }
    out["all_rows"] = analyse(rows)
    for h in sorted({q["h"] for q in rows}):
        out[f"h_{h:g}"] = analyse([q for q in rows if q["h"] == h])
    # the supporting read and the peel-off, with the winning A (or both when undecided)
    sup = []
    for A in A_LIST:
        rd = next(q for q in out["all_rows"]["reads"] if q["A"] == A)
        if "rows" not in rd:
            continue
        s = 1.0 if rd["label"] == "COLLAPSE" else rd["global_s"]
        for h in sorted({q["h"] for q in rows}):
            amps = [row_resid(q, A, s)[1] for q in rows if q["qualifies"] and q["h"] == h]
            if not amps:
                continue
            amp = float(np.mean(amps))
            for q in (x for x in rows if x["h"] == h):
                e6 = float(np.interp(6.0, q["prof"][:, 0], q["prof"][:, 1]))
                rec = {"A": A, "s": s, "tag": q["tag"], "c": q["c"], "h": h, "eps6": e6}
                rec["Q"] = float(e6 / np.exp(log_pred(6.0, q["c"], A, s) + amp))
                if q["qualifies"]:
                    rec["peel"] = peel(q, A, s, row_resid(q, A, s)[1])
                    e = q["prof"][:, 1]
                    k = np.where(
                        (q["prof"][:, 0] >= 3.0) & (QUARTIC_OVER_QUADRATIC * e * e / q["c"] < 1.0)
                    )[0]
                    rec["r_crossover_force_ratio_1"] = (
                        float(q["prof"][k[0], 0]) if len(k) else None
                    )
                sup.append(rec)
    out["supporting"] = sup
    out["profiles"] = {
        q["tag"]: {
            "c": q["c"],
            "n": q["n"],
            "L": q["L"],
            "label": q["label"],
            "at_gate": q["at_gate"],
            "qualifies": q["qualifies"],
            "window_shells": int(q["win"].sum()),
            "r_mean_rms_count": np.round(q["prof"], 8).tolist(),
        }
        for q in rows
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: out[k] for k in out if k not in ("profiles", "supporting")}, indent=1))


# ================= smoke: synthetic fields only =================
def smoke():
    """the instrument on synthetic fields with a known answer; no stored field is opened."""
    from scipy.integrate import solve_ivp
    from scipy.special import kv

    rng = np.random.default_rng(5)
    n, L = 32, 48.0
    h = L / n
    ax = (np.arange(n) - (n - 1) / 2) * h
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    r = np.sqrt(X * X + Y * Y + Z * Z)
    ct = Z / np.maximum(r, 1e-9)
    ang = 1 - ct * ct  # the modulus of an l 2, m 0 spin-2 mode

    def synth(c, A, s, plateau, through_matrix=False):
        beta = beta_of(c) / s**4
        nu = nu_of(A)
        f = lambda q: np.sqrt(q) * kv(nu, np.sqrt(beta) * q * q / 2)  # noqa: E731
        ro = 5.5 * beta**-0.25
        sol = solve_ivp(
            lambda q, y: [y[1], A * y[0] / q**2 + beta * q**2 * y[0]],
            [ro, 2.0],
            [f(ro), (f(ro * 1.00001) - f(ro * 0.99999)) / (2e-5 * ro)],
            rtol=1e-10,
            atol=1e-300,
            dense_output=True,
        )
        rr = np.clip(r, 2.0, ro)
        e = sol.sol(rr.ravel())[0].reshape(r.shape) * (r < ro)
        e = e / sol.sol(6.0)[0] * 0.08
        e = np.minimum(e, plateau) * ang * np.exp(0.03 * rng.normal(size=r.shape))
        if through_matrix:
            M = np.zeros(r.shape + (4, 4))
            M[..., 0, 0] = G_
            M[..., 1, 1], M[..., 2, 2], M[..., 3, 3] = 1.0, D_ + e, D_ - e
            e = eps_of_field(M)
        return shell_profile(e, r, h)

    cases = {}
    for nm, A_true, s_true in (("A3_s1", 3.0, 1.0), ("A3_s1.1", 3.0, 1.1), ("A1_s1.1", 1.0, 1.1)):
        rows = []
        for c in (3e-4, 1e-3, 3e-3):
            prof = synth(c, A_true, s_true, plateau=0.2, through_matrix=(c == 1e-3))
            win = window(prof, c, L)
            rows.append(
                {
                    "tag": f"synthetic_c{c:g}",
                    "c": c,
                    "n": n,
                    "L": L,
                    "h": h,
                    "prof": prof,
                    "win": win,
                    "qualifies": bool(win.sum() >= MIN_SHELLS),
                }
            )
        res = analyse(rows)
        cases[nm] = {
            "truth": {"A": A_true, "s": s_true},
            "overall": res["overall"],
            "labels": {str(q["A"]): q["label"] for q in res["reads"]},
            "pooled_rms_s1": {str(q["A"]): q.get("pooled_rms_s1") for q in res["reads"]},
            "global_s": {str(q["A"]): q.get("global_s") for q in res["reads"]},
            "pooled_rms_global_s": {
                str(q["A"]): q.get("pooled_rms_global_s") for q in res["reads"]
            },
            "free_A": res.get("free_A_not_a_label"),
            "window_shells": [int(q["win"].sum()) for q in rows],
        }
    # a Yukawa halo must not collapse
    rows = []
    for c in (3e-4, 1e-3, 3e-3):
        R = beta_of(c) ** -0.25
        e = 0.08 * np.exp(-(r - 6.0) / R) * 6.0 / np.maximum(r, 1.0) * ang
        prof = shell_profile(e, r, h)
        win = window(prof, c, L)
        rows.append(
            {
                "tag": f"yukawa_c{c:g}",
                "c": c,
                "n": n,
                "L": L,
                "h": h,
                "prof": prof,
                "win": win,
                "qualifies": bool(win.sum() >= MIN_SHELLS),
            }
        )
    res = analyse(rows)
    cases["yukawa"] = {
        "overall": res["overall"],
        "labels": {str(q["A"]): q["label"] for q in res["reads"]},
        "pooled_rms_s1": {str(q["A"]): q.get("pooled_rms_s1") for q in res["reads"]},
    }
    ok = (
        cases["A3_s1"]["overall"] == "MASTER_A3"
        and cases["A3_s1"]["labels"]["3.0"] == "COLLAPSE"
        and cases["A3_s1.1"]["labels"]["3.0"] == "COLLAPSE_SHIFTED"
        and abs(cases["A3_s1.1"]["global_s"]["3.0"] - 1.1) < 0.03
        and cases["A1_s1.1"]["labels"]["1.0"] == "COLLAPSE_SHIFTED"
        and cases["yukawa"]["overall"] == "NEITHER"
    )
    cases["PASS"] = bool(ok)
    with open(OUT_JSON.replace(".json", "_smoke.json"), "w") as f:
        json.dump(cases, f, indent=1)
    print(json.dumps(cases, indent=1))


def plot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with open(OUT_JSON) as f:
        J = json.load(f)
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2))
    for ax, A in zip(axs, A_LIST):
        nu = nu_of(A)
        rd = next(q for q in J["all_rows"]["reads"] if q["A"] == A)
        amps = {p["tag"]: p["log_amplitude_s1"] for p in rd.get("rows", [])}
        rho = np.logspace(-2.5, 1.2, 300)
        ax.loglog(
            rho,
            np.exp(log_master(rho, nu)),
            "k-",
            lw=2,
            label=r"$\rho^\nu K_\nu(\rho)$, normalized",
        )
        for tag, q in sorted(J["profiles"].items(), key=lambda kv_: kv_[1]["c"]):
            if tag not in amps:
                continue
            P = np.array(q["r_mean_rms_count"])
            k = (P[:, 0] >= 3.0) & (P[:, 0] <= 0.5 * q["L"] - 3.0) & (P[:, 1] > 0)
            x = np.sqrt(beta_of(q["c"])) * P[k, 0] ** 2 / 2
            y = P[k, 1] * P[k, 0] ** (2 * nu - 0.5) / np.exp(amps[tag])
            ax.loglog(x, y, ".", ms=4, label=f"c {q['c']:g} n {q['n']} L {q['L']:g}")
        ax.set_xlabel(r"$\rho = \sqrt{\beta}\, r^2 / 2$")
        ax.set_ylabel(r"$\varepsilon\, r^{2\nu - 1/2}$ / row amplitude")
        ax.set_title(f"A = {A:g} (nu = {nu:.3f}): {rd['label']}")
        ax.set_ylim(1e-3, 30)
        ax.legend(fontsize=7)
    fig.tight_layout()
    out = os.path.join(HERE, "..", "plots", "m5_32_r24_1_collapse.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    {"smoke": smoke, "run": run, "plot": plot}[sys.argv[1] if len(sys.argv) > 1 else "smoke"]()
