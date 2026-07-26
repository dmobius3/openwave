"""M5.21.6 INDEPENDENT ADVERSARIAL AUDIT.

Recomputes the M5.21.6 claims with an OWN energy implementation
(no import of m5_21_2b_a_instrument.py / m5_21_6_a_decay.py; those
were only READ to fix conventions). Energy convention audited against
the instrument source:

    A_i^fwd = (M(x + h e_i) - M(x)) / h   [derivative plane at the far
              edge NOT computed -> left 0, per the instrument d1();
              NOTE: the audit brief said "np.roll periodic wrap", that
              is WRONG vs the actual code -> replicated the code]
    A_i^bwd = (M(x) - M(x - h e_i)) / h   [first plane left 0]
    u  = 4 * sum_{i<j} tr(C_ij^T C_ij),  C_ij = [A_i, A_j]
    E_u(sym) = h^3 * 0.5 * (sum u[fwd] + sum u[bwd])
    V  = W2 * sum_i (lam_i - v_i)^2, lam ascending,
         v = sorted([1, delta, 0]),  W2 = 0.000724023879
    E  = E_u + h^3 * sum V          (eps = 0, no Dirichlet)

Outputs ../data/m5_21_6_audit.json + a summary table on stdout.
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless (no plots produced, kept per brief)
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

W2 = 0.000724023879


# ================= OWN energy implementation =================
def d_one(F, ax, h, side):
    """fwd: D[i]=(F[i+1]-F[i])/h, last plane 0.
    bwd: D[i]=(F[i]-F[i-1])/h, first plane 0."""
    D = np.zeros_like(F)
    n = F.shape[ax]
    sl = [slice(None)] * F.ndim
    lo = list(sl)
    hi = list(sl)
    lo[ax] = slice(0, n - 1)
    hi[ax] = slice(1, n)
    diff = (F[tuple(hi)] - F[tuple(lo)]) / h
    if side == "fwd":
        D[tuple(lo)] = diff
    else:
        D[tuple(hi)] = diff
    return D


def u_sum(M, h, side):
    A = [d_one(M, ax, h, side) for ax in range(3)]
    tot = 0.0
    for i in range(3):
        for j in range(i + 1, 3):
            C = A[i] @ A[j] - A[j] @ A[i]
            tot += 4.0 * float(np.sum(C * C))  # sum C_kl^2 = tr(C^T C)
    return tot


def energy_parts(M, h, delta=0.3, w2=W2):
    M = np.asarray(M, dtype=np.float64)
    e_u = h ** 3 * 0.5 * (u_sum(M, h, "fwd") + u_sum(M, h, "bwd"))
    lam = np.linalg.eigvalsh(M)
    vac = np.sort(np.array([1.0, float(delta), 0.0]))
    e_v = h ** 3 * w2 * float(np.sum((lam - vac) ** 2))
    return e_u, e_v


def energy(M, h, delta=0.3, w2=W2):
    e_u, e_v = energy_parts(M, h, delta, w2)
    return e_u + e_v


def load(name):
    return np.load(os.path.join(DATA, name))


def jload(name):
    with open(os.path.join(DATA, name)) as f:
        return json.load(f)


def cell_r(n, h):
    x = (np.arange(n) - (n - 1) / 2.0) * h
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    return np.sqrt(X * X + Y * Y + Z * Z), X, Y, Z


def core_spec(M, n, h, r_c):
    r, _, _, _ = cell_r(n, h)
    lam = np.linalg.eigvalsh(np.asarray(M, dtype=np.float64))
    return np.sort(lam[r < r_c].mean(axis=0))[::-1]


def components(M, n, h, thr=0.09, r_in=15.6, r_edge=14.88, conn=26):
    """own topology read: min eigen-gap mask -> labeled components."""
    lam = np.linalg.eigvalsh(np.asarray(M, dtype=np.float64))
    gap = np.minimum(lam[..., 1] - lam[..., 0],
                     lam[..., 2] - lam[..., 1])
    r, X, Y, Z = cell_r(n, h)
    mask = (gap < thr) & (r < r_in)
    st = np.ones((3, 3, 3)) if conn == 26 else \
        ndimage.generate_binary_structure(3, 1)
    lab, nc = ndimage.label(mask, structure=st)
    rho = np.sqrt(X * X + Y * Y)
    out = []
    for k in range(1, nc + 1):
        sel = lab == k
        out.append({"size": int(sel.sum()),
                    "r_max": float(r[sel].max()),
                    "rho_c": float(rho[sel].mean()),
                    "z_c": float(Z[sel].mean()),
                    "compact": bool(r[sel].max() <= r_edge)})
    out.sort(key=lambda c: -c["size"])
    return out


RESULTS = []


def add(claim, method, computed, expected, verdict, note=""):
    RESULTS.append({"claim": claim, "method": method,
                    "computed": computed, "expected": expected,
                    "verdict": verdict, "note": note})


def main():
    # ---------------- CLAIM 1: pinned levels ----------------
    levels = [("m5_21_6_end_p32_B.npz", 55.599677, "p32_B"),
              ("m5_21_6_end_p32_C.npz", 15.934288, "p32_C")]
    vals = {}
    for fn, exp, tag in levels:
        z = load(fn)
        E = energy(z["M"], float(z["h"]), float(z["delta"]))
        vals[tag] = E
        rel = abs(E - exp) / abs(exp)
        add("1 pinned level " + tag,
            "own energy (sym stencil, T2 potential) on stored float32 M",
            round(E, 6), exp,
            "CONFIRMED" if rel < 2e-4 else "REFUTED",
            f"rel diff {rel:.2e}")
    # i2_A_T2: the stated W2 does NOT reproduce 6.604177; the 2b row
    # json records w2 = 0.0027581 for this run. Audit both weights.
    zA = load("m5_21_2b_end_i2_A_T2.npz")
    rowA = jload("m5_21_2b_row_i2_A_T2.json")
    w2_run = float(rowA["w2"])
    eu_A, ev_A_stated = energy_parts(zA["M"], float(zA["h"]),
                                     float(zA["delta"]), W2)
    E_A_stated = eu_A + ev_A_stated
    E_A_runw2 = energy(zA["M"], float(zA["h"]), float(zA["delta"]),
                       w2=w2_run)
    vals["i2_A_T2"] = E_A_runw2
    rel_run = abs(E_A_runw2 - 6.604177) / 6.604177
    rel_eu = abs(eu_A - rowA["E_u"]) / abs(rowA["E_u"])
    add("1 pinned level i2_A_T2",
        "own energy, stated W2 AND the run-recorded w2=0.0027581",
        f"E(W2 stated)={E_A_stated:.6f}, "
        f"E(run w2)={E_A_runw2:.6f}, own E_u={eu_A:.6f}",
        6.604177,
        "PARTIAL" if (rel_run < 2e-4 and rel_eu < 2e-4) else "REFUTED",
        "claim as stated REFUTED: with W2=0.000724023879 the file "
        "scores 5.318035, not 6.604177. The 2b run json records "
        "w2=0.0027581 (the n=32 calib match); under THAT weight the "
        f"claimed value reproduces (rel {rel_run:.1e}) and own E_u "
        f"matches the row E_u=4.860236 (rel {rel_eu:.1e}). So the "
        "number is genuine but the A level uses a ~3.8x heavier "
        "potential weight than B and C -> the three levels are NOT "
        "on one common energy scale")
    order_ok = (E_A_runw2 < vals["p32_C"] < vals["p32_B"]
                and E_A_stated < vals["p32_C"] < vals["p32_B"])
    add("1 ordering A < C < B", "compare recomputed levels",
        f"A={E_A_runw2:.4f} (or {E_A_stated:.4f} at stated W2) "
        f"C={vals['p32_C']:.4f} B={vals['p32_B']:.4f}", "A < C < B",
        "CONFIRMED" if order_ok else "REFUTED",
        "ordering holds under BOTH weight conventions")

    # ---------------- CLAIM 2: kick returns ----------------
    all_de = []
    for jf in ("m5_21_6_p1_p32_B.json", "m5_21_6_p1_p32_C.json"):
        J = jload(jf)
        all_de += [r["dE_rel"] for r in J["rows"]]
    n_runs = len(all_de)
    mx = max(abs(d) for d in all_de)
    add("2 JSON: 12 runs |dE_rel| < 0.002",
        "read both p1 JSONs",
        f"{n_runs} runs, max |dE_rel| = {mx:.2e}",
        "12 runs, all < 0.002",
        "CONFIRMED" if (n_runs == 12 and mx < 0.002) else "REFUTED")
    for tag, jf, start_fn in (
            ("p32_B", "m5_21_6_p1_p32_B.json", "m5_21_6_end_p32_B.npz"),
            ("p32_C", "m5_21_6_p1_p32_C.json", "m5_21_6_end_p32_C.npz")):
        J = jload(jf)
        row = [r for r in J["rows"]
               if r["family"] == "K1" and r["amp"] == 0.4][0]
        z = load(f"m5_21_6_end_{tag}_K1_0.4.npz")
        E_end = energy(z["M"], float(z["h"]), float(z["delta"]))
        zs = load(start_fn)
        E_start = energy(zs["M"], float(zs["h"]), float(zs["delta"]))
        rel_json = abs(E_end - row["E_end"]) / abs(row["E_end"])
        de = (E_end - E_start) / E_start
        ok = rel_json < 2e-4 and abs(de) < 0.002
        add(f"2 recompute {tag}_K1_0.4 endpoint",
            "own energy on kicked-relaxed endpoint + start",
            f"E_end={E_end:.6f} (json {row['E_end']:.6f}), "
            f"own dE_rel={de:.2e}",
            f"json dE_rel={row['dE_rel']:.2e}, |dE_rel|<0.002",
            "CONFIRMED" if ok else "REFUTED",
            f"rel diff vs json E_end {rel_json:.2e}")

    # ---------------- CLAIM 3: free-arena endpoints ----------------
    f48 = {}
    for tag, expE in (("A", 1.946), ("B", 0.853), ("C", 2.099)):
        z = load(f"m5_21_6_end_f48_{tag}.npz")
        f48[tag] = z
        E = energy(z["M"], float(z["h"]), float(z["delta"]))
        rel = abs(E - expE) / abs(expE)
        add(f"3 f48_{tag} endpoint E",
            "own energy, n=48 h=1.0", round(E, 6), expE,
            "CONFIRMED" if rel < 1e-3 else "REFUTED",
            f"rel diff {rel:.2e} (claim quoted to 3 decimals)")
    spec_exp = {"A": (0.77, 0.49, 0.04), "B": (1.00, 0.22, 0.07),
                "C": (0.77, 0.51, 0.02)}
    specs = {}
    for tag in "ABC":
        s = core_spec(f48[tag]["M"], 48, 1.0, r_c=6.0)
        specs[tag] = s
        dev = max(abs(a - b) for a, b in zip(s, spec_exp[tag]))
        add(f"3 f48_{tag} core spectrum (r<6)",
            "own eigvalsh mean over r<6",
            [round(float(x), 3) for x in s], list(spec_exp[tag]),
            "CONFIRMED" if dev < 0.05 else "REFUTED",
            f"max dev vs claim {dev:.3f}")
    ac = max(abs(a - c) for a, c in zip(specs["A"], specs["C"]))
    bv = [round(float(abs(x - v)), 3)
          for x, v in zip(specs["B"], (1.0, 0.3, 0.0))]
    add("3 A == C within ~0.03; B near vacuum",
        "compare own spectra",
        f"max|A-C|={ac:.3f}; |B-vac|={bv}",
        "|A-C| <~0.03; B near (1, 0.3, 0)",
        "CONFIRMED" if ac < 0.04 else "REFUTED",
        "B mid-eigenvalue sits 0.08 below vacuum 0.3 (claim itself "
        "says 'near', shift noted)")

    # ---------------- CLAIM 4: topology ----------------
    counts26 = {}
    counts6 = {}
    for tag, exp_n in (("A", 2), ("B", 0), ("C", 2)):
        comp = components(f48[tag]["M"], 48, 1.0, conn=26)
        comp6 = components(f48[tag]["M"], 48, 1.0, conn=6)
        nc = sum(1 for c in comp if c["compact"])
        nc3 = sum(1 for c in comp if c["compact"] and c["size"] >= 3)
        counts26[tag] = nc
        counts6[tag] = sum(1 for c in comp6 if c["compact"])
        add(f"4 f48_{tag} compact components (26-conn)",
            "own gap<0.09 & r<15.6 mask + scipy label",
            f"{nc} compact (sizes "
            f"{[c['size'] for c in comp if c['compact']]}); "
            f"{nc3} with size>=3", exp_n,
            "CONFIRMED" if nc == exp_n else
            ("PARTIAL" if nc3 == exp_n else "REFUTED"),
            f"total components {len(comp)}")
    add("4 robustness: 6-connectivity",
        "same mask, 6-conn label",
        f"A={counts6['A']} B={counts6['B']} C={counts6['C']} "
        f"(26-conn: A={counts26['A']} B={counts26['B']} "
        f"C={counts26['C']})",
        "report whether counts change",
        "CONFIRMED" if counts6 == counts26 else "PARTIAL",
        "counts unchanged under 6-conn" if counts6 == counts26
        else "counts CHANGE under 6-conn")
    comp10k = components(f48["C"]["M_it10000"], 48, 1.0, conn=26)
    eq = [c for c in comp10k
          if abs(c["z_c"]) < 3.0 and c["rho_c"] > 13.0]
    eq3 = [c for c in eq if c["size"] >= 3]
    sizes3 = sorted(c["size"] for c in eq3)
    ok_eq = (len(eq3) == 2 and sizes3[0] == sizes3[-1]
             and 5 <= sizes3[0] <= 12
             and all(abs(c["rho_c"] - 15.2) < 1.0
                     and abs(c["z_c"]) < 1.0 for c in eq3)) \
        if eq3 else False
    add("4d f48_C it10000 equatorial arcs",
        "own components; select rho_c>13, |z_c|<3; size>=3 "
        "(the instrument's loop_read convention)",
        f"{len(eq3)} arcs size>=3: sizes {sizes3}, "
        f"rho_c {[round(c['rho_c'], 1) for c in eq3]}, "
        f"z_c {[round(c['z_c'], 2) for c in eq3]}; raw count in "
        f"band {len(eq)} (extra sizes "
        f"{sorted(c['size'] for c in eq if c['size'] < 3)})",
        "2 arcs, equal size ~8, rho~15.2, z~0",
        "PARTIAL" if ok_eq else "REFUTED",
        "the two size-8 equal arcs EXIST exactly as claimed "
        "(rho_c 15.2, z_c 0.0), but the raw 26-conn count in the "
        "band is 4: two additional 2-voxel specks at rho 15.6, "
        "z +-0.5. The 'exactly two' phrasing holds only under the "
        "size>=3 filter the source pipeline (loop_read) applies")

    # ---------------- CLAIM 5: evolve ledger ----------------
    ev = jload("m5_21_6_ev_C_free.json")
    fin = ev["hist"][-1]
    ledger = fin["E"] + fin["KE"] + fin["absorbed"]
    E_start_own = energy(f48["C"]["M_it1000"], 1.0, 0.3)
    ok5a = abs(ledger - 5.064) < 0.01 and \
        abs(ledger - E_start_own) < 0.01
    add("5 ledger E+KE+absorbed == start E",
        "sum final JSON row; own E of f48_C M_it1000",
        f"sum={ledger:.4f}; own start E={E_start_own:.4f}",
        "5.064 within 0.01",
        "CONFIRMED" if ok5a else "REFUTED",
        f"leak {ledger - E_start_own:+.2e}")
    ok5b = (abs(fin["E"] - 3.5816) < 5e-4
            and abs(fin["KE"] - 0.3337) < 5e-4
            and abs(fin["absorbed"] - 1.1491) < 5e-4)
    add("5 final-row values",
        "read JSON final row",
        f"E={fin['E']:.4f} KE={fin['KE']:.4f} "
        f"abs={fin['absorbed']:.4f}",
        "E=3.5816 KE=0.3337 absorbed=1.1491",
        "CONFIRMED" if ok5b else "REFUTED")
    zev = load("m5_21_6_ev_C_free.npz")
    E_ev_own = energy(zev["M"], float(zev["h"]), float(zev["delta"]))
    rel = abs(E_ev_own - fin["E"]) / abs(fin["E"])
    add("5 evolve endpoint static E",
        "own energy on ev_C_free.npz M",
        round(E_ev_own, 6), round(fin["E"], 6),
        "CONFIRMED" if rel < 2e-4 else "REFUTED",
        f"rel diff {rel:.2e}")

    # ---------------- CLAIM 6: rotation vs melt ----------------
    rots = [r["rot_core_deg"] for r in ev["hist"]]
    sdevs = [r["spec_dev_core"] for r in ev["hist"]]
    ok6a = (abs(rots[0] - 0.16) < 0.02 and abs(rots[-1] - 4.27) < 0.05
            and min(sdevs) >= 0.14 and max(sdevs) <= 0.19)
    add("6 JSON: rot 0.16 -> ~4.27, spec_dev in [0.14, 0.19]",
        "read ev JSON history",
        f"rot {rots[0]:.2f} -> {rots[-1]:.2f} "
        f"(peak {max(rots):.2f}); spec_dev "
        f"[{min(sdevs):.3f}, {max(sdevs):.3f}]",
        "0.16 -> 4.27; [0.14, 0.19]",
        "CONFIRMED" if ok6a else "REFUTED",
        "rot peaks 4.61 at it2750 then settles 4.27 "
        "(claim 'grows to ~4.27' is the endpoint, not monotone)")
    Ms = np.asarray(f48["C"]["M_it1000"], dtype=np.float64)
    Me = np.asarray(zev["M"], dtype=np.float64)
    v0 = np.linalg.eigh(Ms)[1][..., -1]
    v1 = np.linalg.eigh(Me)[1][..., -1]
    r, _, _, _ = cell_r(48, 1.0)
    dots = np.abs(np.sum(v0 * v1, axis=-1))
    rot_own = float(np.degrees(np.mean(
        np.arccos(np.clip(dots[r < 10.0], 0.0, 1.0)))))
    add("6 independent final rotation (r<10)",
        "own eigh leading vectors, mean arccos|dot|",
        round(rot_own, 3), "~4.27 within 0.5",
        "CONFIRMED" if abs(rot_own - 4.27) < 0.5 else "REFUTED")

    # ---------------- write + print ----------------
    out = os.path.join(DATA, "m5_21_6_audit.json")
    with open(out, "w") as f:
        json.dump(RESULTS, f, indent=1)
    wc = max(len(r["claim"]) for r in RESULTS)
    print(f"\n{'CLAIM':<{wc}}  VERDICT    COMPUTED")
    print("-" * (wc + 60))
    for r in RESULTS:
        print(f"{r['claim']:<{wc}}  {r['verdict']:<9}  "
              f"{r['computed']}")
    nref = sum(1 for r in RESULTS if r["verdict"] == "REFUTED")
    npar = sum(1 for r in RESULTS if r["verdict"] == "PARTIAL")
    print(f"\n{len(RESULTS)} checks: "
          f"{len(RESULTS) - nref - npar} CONFIRMED, "
          f"{npar} PARTIAL, {nref} REFUTED -> {out}")


if __name__ == "__main__":
    main()
