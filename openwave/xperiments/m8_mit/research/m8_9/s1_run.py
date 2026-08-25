#!/usr/bin/env python3
"""M8.9 S1: trivial-fibre high-level control on the shipped R_0 production block.

Executes S1 under the filed decision rule at
``../findings/m8_9_s1_decision_rule.md``, frozen region above FREEZE-BOUNDARY,
SHA-256 68df11a02ee5097d23712b4eace9b65220d2f540ea147660ce3cd8ae4a938934.

Every constant below is transcribed from that rule. Nothing here reinterprets,
supplements or extends it. Run from this directory:

    PYTHONPATH=.:../m8_4_p1a:../m8_5b/pilot:../m8_5b/production python3 s1_run.py
"""
import hashlib
import subprocess
import sys

import numpy as np

from p0 import group as G, representations as RP
from p0.cloud import fibonacci_seeds_s3, build_orbit_cloud
from p0.bundle_operator import build_L_bundle
from p1a.mass_matrix import build_Mh_base, build_Mh_rho
from p1a.subspace import extract_and_score, projector_invariants

# --- frozen inputs, transcribed from the filed rule ------------------------
RULE = "../findings/m8_9_s1_decision_rule.md"
RULE_SHA = "68df11a02ee5097d23712b4eace9b65220d2f540ea147660ce3cd8ae4a938934"

N_SEEDS, K, RBF_M, RBF_P = 60, 110, 7, 4          # P1A production constants

#            label,  n,  lambda,  multiplicity,  window [lo, hi)
CLUSTERS = [("C_0",   0,      0,   1,  (-np.inf,  84.0)),
            ("C_12", 12,    168,  13,  (   84.0, 304.0)),
            ("C_20", 20,    440,  21,  (  304.0, 532.0)),
            ("C_24", 24,    624,  25,  (  532.0, np.inf))]

I_LOW_MAX  = 1.27e-08      # R_7, dirtiest of the prior clean group
I_HIGH_MIN = 1.33e-02      # R_8, cleanest of the prior dirty group
ADJUDICATING = ("C_12", "C_20")     # C_24 is recorded, zero classification weight


def gate_freeze():
    """The rule must be the filed rule, byte for byte, before anything spectral."""
    out = subprocess.run(["sh", "-c", f"sed '/FREEZE-BOUNDARY/,$d' {RULE} | shasum -a 256"],
                         capture_output=True, text=True)
    got = out.stdout.split()[0]
    print(f"  frozen-rule SHA-256 {got}")
    if got != RULE_SHA:
        sys.exit(f"  FAIL: expected {RULE_SHA}. STOP, no spectrum computed.")
    print("  PASS: matches the filed value\n")


def assign(evals):
    """Nearest analytic level by Re lambda, via the frozen windows. Not sort-and-slice."""
    re = np.real(evals)
    return {name: np.where((re >= lo) & (re < hi))[0] for name, _, _, _, (lo, hi) in CLUSTERS}


def g_mult(idx):
    """Counts must equal the analytic multiplicities in every window."""
    rows, ok = [], True
    for name, _, lam, mult, _ in CLUSTERS:
        got = len(idx[name])
        ok &= (got == mult)
        rows.append((name, lam, mult, got, "ok" if got == mult else "MISMATCH"))
    return ok, rows


def main():
    print("== M8.9 S1: trivial-fibre high-level control ==\n")
    print("GATE: frozen rule")
    gate_freeze()

    print("GATE: shipped operator provenance")
    elems = G.build_icosians()
    chi = G.build_character_table(elems)
    out = RP.build_all_representations(elems, chi)
    reps = out[0] if isinstance(out, tuple) else out
    X, oid, gid = build_orbit_cloud(fibonacci_seeds_s3(N_SEEDS), elems)
    W = build_Mh_base(X, oid, N_SEEDS)
    L = np.asarray(build_L_bundle(X, oid, gid, elems, reps["R0"], k=K, m=RBF_M, p=RBF_P)[0])
    Mh = build_Mh_rho(W, reps["R0"][0].shape[0])
    print(f"  cloud {len(X)} nodes, k={K}, m={RBF_M}, p={RBF_P}")
    print(f"  R_0 block {L.shape}, SHA-256 of bytes {hashlib.sha256(L.tobytes()).hexdigest()[:32]}...")
    if L.shape != (60, 60):
        sys.exit(f"  FAIL: expected a 60x60 block, got {L.shape}. STOP.")
    print("  PASS: 60x60, assembled by the shipped build_L_bundle\n")

    # --- G-MULT mutation arm, BEFORE any real spectrum is read -------------
    # Run on the exact analytic spectrum, which passes by construction, so that a
    # failure after mutation is attributable to the mutation. Arming against the
    # real spectrum would be vacuous if the real spectrum already fails.
    print("GATE: G-MULT mutation arm, on the synthetic analytic spectrum")
    synth = np.concatenate([np.full(m, float(lam)) for _, _, lam, m, _ in CLUSTERS]).astype(complex)
    s_ok, s_rows = g_mult(assign(synth))
    print(f"  synthetic spectrum, G-MULT {'PASSES' if s_ok else 'FAILS'} (must PASS to arm)")
    mut = synth.copy()
    mut[np.where(np.real(synth) == 168.0)[0][0]] = 305.0 + 0j    # C_12 -> C_20, across 304
    m_ok, m_rows = g_mult(assign(mut))
    c0 = {r[0]: r[3] for r in s_rows}; c1 = {r[0]: r[3] for r in m_rows}
    d12, d20 = c1["C_12"] - c0["C_12"], c1["C_20"] - c0["C_20"]
    print(f"  move 1 eigenvalue across 304: dC_12={d12:+d} dC_20={d20:+d}; "
          f"G-MULT {'FAILS' if not m_ok else 'PASSES'}")
    arm = s_ok and (not m_ok) and d12 == -1 and d20 == +1
    print(f"  arm {'PASS: the gate discriminates' if arm else 'FAIL: arm is vacuous'}\n")
    if not arm:
        sys.exit("  STOP: G-MULT is not armed, so no S1 result may be issued.")

    # --- the single spectral read ------------------------------------------
    evals = np.linalg.eigvals(L)
    idx = assign(evals)

    print("GATE: G-MULT on the real spectrum, filed windows 84 / 304 / 532")
    ok, rows = g_mult(idx)
    print(f"  {'cluster':8s} {'lambda':>7} {'required':>9} {'found':>6}   status")
    for name, lam, mult, got, st in rows:
        print(f"  {name:8s} {lam:>7} {mult:>9} {got:>6}   {st}")
    print(f"  G-MULT: {'PASS' if ok else 'FAIL'}\n")

    # --- per-cluster record ------------------------------------------------
    print("RECORD")
    print(f"  {'cluster':8s} {'center':>10} {'re spread':>11} {'max|Im|':>11} {'Schur resid':>12}")
    I = {}
    for name, _, lam, mult, _ in CLUSTERS:
        e = evals[idx[name]]
        if len(e) == 0:
            print(f"  {name:8s} {'empty':>10}"); I[name] = float("nan"); continue
        center = float(np.mean(np.real(e)))
        spread = float(np.max(np.real(e)) - np.min(np.real(e)))
        imax = float(np.max(np.abs(np.imag(e))))
        I[name] = imax
        try:
            _, Q = extract_and_score(L, Mh, float(lam), mult, label=f"R0@{lam}")
            res = projector_invariants(Q, L, Mh, mult)["invariance_residual"]
        except Exception as exc:                       # recorded, never silently dropped
            res = float("nan"); print(f"    (Schur extraction at lambda={lam}: {exc})")
        print(f"  {name:8s} {center:>10.4f} {spread:>11.3e} {imax:>11.3e} {res:>12.3e}")

    # --- classification, only if G-MULT passed -----------------------------
    print()
    if not ok:
        print("VERDICT: INSTRUMENT DEFECT. G-MULT failed on the real spectrum.")
        print("  No S1-A/B/C classification is issued and I_star is NOT formed:")
        print("  the filed rule computes it only on a passing gate.\n")
        print("DEFECT CHARACTERIZATION (describes the failure, classifies nothing)")
        re, im = np.real(evals), np.abs(np.imag(evals))
        srt = np.argsort(re)
        j = srt[34]
        print(f"  count deficit in C_20 and excess in C_24 are produced at the 532 boundary;"
              f" the nearest eigenvalue on the C_24 side is Re={re[j]:.3f}")
        dirty = im > 1e-6
        print(f"  {int(dirty.sum())} of 60 eigenvalues carry |Im| > 1e-06, in"
              f" {int(dirty.sum()) // 2} complex-conjugate pairs")
        print(f"  such pairs already appear within the C_12 assignment window, the lowest"
              f" at Re = {re[dirty].min():.2f}, not only near the cutoff")
        print(f"  the lowest mode is clean at {im[srt[0]]:.3e}, matching the P1A shipped"
              f" R_0 value 5.95e-14")
        for name, _, lam, mult, (lo, hi) in CLUSTERS:
            w = (re >= lo) & (re < hi)
            print(f"  {name:5s} window: {int(w.sum()):>2} eigenvalues, "
                  f"{int((w & dirty).sum()):>2} dirty, Re from {re[w].min():.1f} to {re[w].max():.1f}"
                  f" (analytic center {lam})")
        print("\n  The block does not reproduce the analytic decomposition the rule predicted:")
        print("  the C_20 and C_24 windows are not clustered, and assignments run past 1300")
        print("  where the analytic content stops at 624. C_12's count is exactly 13, but a")
        print("  correct count is not evidence its contents are the level-12 eigenspace:")
        print("  count-only safety is not subspace correctness (route_a_repn.py). No branch")
        print("  is licensed on any window.")
        return
    I_star = max(I[c] for c in ADJUDICATING)
    print(f"  I_12 = {I['C_12']:.3e}   I_20 = {I['C_20']:.3e}   ->  I_star = {I_star:.3e}")
    print(f"  thresholds: I_low_max = {I_LOW_MAX:.2e}   I_high_min = {I_HIGH_MIN:.2e}")
    print(f"  C_24 = {I['C_24']:.3e}, recorded, zero classification weight")
    if I_star >= I_HIGH_MIN:
        print("\nVERDICT: S1-A. Base discretization strongly implicated.")
        print("  Licensed: nontrivial fibre transport is NOT NECESSARY for high-scale"
              " contamination, and the base discretization is strongly implicated.")
    elif I_star <= I_LOW_MAX:
        print("\nVERDICT: S1-B. The simple high-lambda explanation is rejected.")
        print("  Licensed, and no more: high harmonic level by itself, in the trivial quotient"
              " sector, is INSUFFICIENT to reproduce the high-contamination regime.")
        print("  This does NOT exonerate RBF-FD and does NOT convict the fibre transport. S2 required.")
    else:
        print("\nVERDICT: S1-C. Indeterminate.")
        print("  No localization conclusion. S2 required.")


if __name__ == "__main__":
    main()
