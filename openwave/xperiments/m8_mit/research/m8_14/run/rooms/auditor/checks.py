"""Checks.  Every check is a function returning (ok, detail).  Each is run on
the real computation and then on at least one planted defect, which must make
it FAIL; a planted defect that does not fire is itself reported as a failure.
Writes results_checks.json."""
import math

import numpy as np

import fem1d_sectors as B
import fem2d_q2 as A
from common import WIDTHS, NEIG, load_json, richardson, save_json

LOG = []


def report(name, ok, detail, planted=False):
    tag = ("PASS" if ok else "FAIL") if not planted else ("FIRED" if not ok else "NOT-FIRED")
    line = "%s  %s%s: %s" % (tag, "[planted] " if planted else "", name, detail)
    print(line, flush=True)
    LOG.append({"check": name, "planted": planted, "ok": bool(ok), "detail": detail})


# ---- C1: 1D assembly / quadrature identities ------------------------------
def c1(N=16, M=8, W=1.0, nquad=A.NQUAD):
    Ay, By, My = A.y_matrices(N, A.G, nquad)
    Kw, Mw = A.w_matrices(W, M)
    e1 = abs(My.sum() - 2.0)                  # int_0^pi |cos y| dy = 2
    e2 = np.abs(Ay @ np.ones(Ay.shape[0])).max()  # constants in kernel
    e3 = np.abs(Kw @ np.ones(Kw.shape[0])).max()
    e4 = abs(Mw.sum() - 2 * W)
    worst = max(e1, e2, e3, e4)
    return worst < 1e-12, "max identity defect %.2e (tol 1e-12)" % worst


# ---- C2: eigen-residual at own level --------------------------------------
def c2(W=1.7, N=32, perturb=0.0):
    r = A.solve(W, N, return_vectors=True)
    Kr, Mr, X = r["Kr"], r["Mr"], r["vecs"]
    worst = 0.0
    for i in range(NEIG):
        lam = r["vals"][i] * (1 + perturb)
        x = X[:, i]
        worst = max(worst, np.linalg.norm(Kr @ x - lam * (Mr @ x)) / np.linalg.norm(Kr @ x))
    return worst < 1e-9, "max relative residual %.2e (tol 1e-9)" % worst


# ---- C3: seam / boundary data of reconstructed eigenfunctions --------------
def c3(W=1.4, N=16, seam_sign=-1.0, reflect=True):
    r = A.solve(W, N, seam_sign=seam_sign, reflect=reflect, return_vectors=True)
    M = r["M"]; nw = 2 * M + 1; ny = 4 * N + 1
    worst = 0.0
    for i in range(NEIG):
        u = (r["P"] @ r["vecs"][:, i]).reshape(ny, nw)
        sc = np.abs(u).max()
        # spec: phi(pi, -w) = -phi(0, w); node w_j has -w_j = w_{2M-j}
        seam = np.abs(u[4 * N, ::-1] + u[0, :]).max() / sc
        bdry = max(np.abs(u[:, 0]).max(), np.abs(u[:, -1]).max(), np.abs(u[2 * N, :]).max()) / sc
        worst = max(worst, seam, bdry)
    return worst < 1e-12, "max seam/boundary violation %.2e (tol 1e-12)" % worst, r["vals"]


# ---- C4: asymptotic-regime test on a level sequence ------------------------
def c4_seq(seq):
    d = np.diff(seq)
    mono = all(d[i] * d[i + 1] > 0 and abs(d[i + 1]) < abs(d[i]) for i in range(len(d) - 1))
    p1 = richardson(seq)["p"]; p0 = richardson(seq[:-1])["p"]
    stable = p1 is not None and p0 is not None and abs(p1 - p0) < 0.1
    return mono and stable, "monotone=%s p_last=%s p_prev=%s" % (mono, p1, p0)


# ---- C5: cross-method agreement within the stated errors -------------------
def c5(a_vals, a_err, b_vals, b_err):
    worst = 0.0
    for x, ex, y, ey in zip(a_vals, a_err, b_vals, b_err):
        worst = max(worst, abs(x - y) / (ex + ey + 1e-12 * abs(y)))
    return worst <= 1.0, "max |diff|/(errA+errB) = %.3g (must be <= 1)" % worst


# ---- C6: sector cutoff (adding the next sector leaves lowest six unchanged)
def c6(W, low, used, N=8):
    k = max(used) + 1
    _, extra = B.sector_eigs(W, k, N)
    merged = sorted(list(low) + list(extra))[:NEIG]
    worst = max(abs(x - y) for x, y in zip(merged, low))
    return worst == 0.0, "next sector k=%d changes lowest six by %.2e" % (k, worst)


def main():
    RA = load_json("results_A.json")["widths"]
    RB = load_json("results_B.json")["widths"]
    RC = load_json("results_C.json")
    extra = {}

    ok, d = c1(); report("C1 assembly identities", ok, d)
    ok, d = c1(nquad=1); report("C1 assembly identities, 1-point quadrature", ok, d, planted=True)

    ok, d = c2(); report("C2 eigen-residual W=17/10 N=32", ok, d)
    ok, d = c2(perturb=1e-6); report("C2 eigen-residual with eigenvalue shifted 1e-6", ok, d, planted=True)

    ok, d, v = c3(); report("C3 seam phi(pi,-w)=-phi(0,w), zero on w=+-W and fiber", ok, d)
    extra["C3_correct_vals_W7/5_N16"] = [float(x) for x in v]
    ok, d, v = c3(seam_sign=+1.0); report("C3 with seam sign flipped (+)", ok, d, planted=True)
    extra["planted_seam_plus_vals_W7/5_N16"] = [float(x) for x in v]
    print("      eigenvalues with flipped seam sign:", " ".join("%.12f" % x for x in v))
    ok, d, v = c3(reflect=False); report("C3 with seam reflection w->-w omitted", ok, d, planted=True)
    extra["planted_noreflect_vals_W7/5_N16"] = [float(x) for x in v]
    print("      eigenvalues without reflection: ", " ".join("%.12f" % x for x in v))

    for label, _ in WIDTHS:
        for i in range(NEIG):
            seq = [lv["lowest"][i] for lv in RA[label]["levels"]]
            ok, d = c4_seq(seq); report("C4 A asymptotic regime W=%s eig%d" % (label, i + 1), ok, d)
    seq = [lv["lowest"][0] for lv in RA["17/10"]["levels"]]
    bad = [x + (1e-6 if j % 2 else -1e-6) for j, x in enumerate(seq)]
    ok, d = c4_seq(bad); report("C4 on W=17/10 bottom sequence with +-1e-6 alternating noise", ok, d, planted=True)

    for label, W in WIDTHS:
        a = [r["extrap"] for r in RA[label]["richardson"]]
        ea = [r["err_combined"] for r in RA[label]["richardson"]]
        b, eb = RB[label]["best"], RB[label]["err"]
        ok, d = c5(a, ea, b, eb); report("C5 A vs B, W=%s" % label, ok, d)
        c = RC[label]["lowest"]
        ok, d = c5(b, eb, c, [0.0] * NEIG); report("C5' B vs closed form, W=%s" % label, ok, d)
        ok, d = c5(a, ea, c, [0.0] * NEIG); report("C5' A vs closed form, W=%s" % label, ok, d)
    # planted: A with the seam replaced by zero data on y=0 and y=pi, run
    # through the same level/Richardson pipeline (N = 32, 64, 128), vs B.
    # (Seam-sign flip and omitted reflection are NOT used here: C3 shows they
    # leave the eigenvalues unchanged, so they cannot be caught by C5.)
    # INFO (not a planted defect): seam_sign=0 puts zero data on y=pi but
    # leaves y=0 free (natural condition).  Found to be isospectral, since the
    # surface is symmetric about the seam line; recorded, not counted.
    seqs = [A.solve(1.4, n, seam_sign=0.0)["vals"] for n in (32, 64, 128)]
    rr = [richardson([s[i] for s in seqs]) for i in range(NEIG)]
    ok, d = c5([r["extrap"] for r in rr], [r["err"] for r in rr], RB["7/5"]["best"], RB["7/5"]["err"])
    print("INFO  seam replaced by zero data on y=pi only (y=0 free), W=7/5, extrapolated "
          "N=32,64,128, vs B: %s -> %s" % (d, "indistinguishable" if ok else "different"))
    extra["info_seam_cut_levels_W7/5_N32_64_128"] = [[float(x) for x in s] for s in seqs]
    extra["info_seam_cut_extrap_W7/5"] = [r["extrap"] for r in rr]
    # planted: width misread as the full width 2W (A solved at W/2), same pipeline
    seqs = [A.solve(0.7, n)["vals"] for n in (32, 64, 128)]
    rr = [richardson([s[i] for s in seqs]) for i in range(NEIG)]
    ok, d = c5([r["extrap"] for r in rr], [r["err"] for r in rr], RB["7/5"]["best"], RB["7/5"]["err"])
    report("C5 A(width misread: solved at W/2, extrapolated N=32,64,128) vs B, W=7/5", ok, d, planted=True)
    extra["planted_halfwidth_extrap_W7/5"] = [r["extrap"] for r in rr]
    # planted: B with mu_k = (k pi / W)^2 instead of (k pi/(2W))^2
    lowbad, _ = B.spectrum(1.4, 8, mu_scale=4.0)
    ok, d = c5([p[0] for p in lowbad], [1e-9] * NEIG, RB["7/5"]["best"], RB["7/5"]["err"])
    report("C5 B(wrong transverse eigenvalue) vs B, W=7/5", ok, d, planted=True)
    extra["planted_wrongmu_vals_W7/5_N8"] = [p[0] for p in lowbad]

    for label, W in WIDTHS:
        lv = RB[label]["levels"][3]  # N=8 level
        ok, d = c6(W, lv["lowest"], lv["sectors"], N=lv["N"]); report("C6 sector cutoff W=%s" % label, ok, d)
    low1, _ = B.spectrum(math.pi / 2, 8)
    _, only1 = B.sector_eigs(math.pi / 2, 1, 8)
    ok, d = c6(math.pi / 2, list(only1), [1], N=8)
    report("C6 sector cutoff with only k=1 kept, W=pi/2", ok, d, planted=True)

    real = [x for x in LOG if not x["planted"]]
    planted = [x for x in LOG if x["planted"]]
    print("SUMMARY real checks: %d PASS, %d FAIL; planted defects: %d FIRED, %d NOT-FIRED" % (
        sum(x["ok"] for x in real), sum(not x["ok"] for x in real),
        sum(not x["ok"] for x in planted), sum(x["ok"] for x in planted)))
    save_json("results_checks.json", {"log": LOG, "extra_numbers": extra})


if __name__ == "__main__":
    main()
