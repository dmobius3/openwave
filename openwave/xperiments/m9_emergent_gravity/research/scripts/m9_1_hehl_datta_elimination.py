#!/usr/bin/env python3
"""M9.1 solver: extract the Hehl-Datta coefficient from Palatini + Dirac.

The contraction returns a float. Comparison to 3/16 happens only after
extraction (task C1). Mutations (C5) must fail.

Writes ../data/m9_1_hehl_datta.json (repo-relative paths only).
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
sys.path.insert(0, HERE)

from hehl_datta import (  # noqa: E402
    check_clifford,
    contact_lagrangian,
    mutated_ratio,
    random_spinor,
    signature_mostly_minus,
    signature_mostly_plus,
)

TARGET = 3.0 / 16.0
TOL = 1e-12
N_SAMPLES = 8
SEED = 20260815
KAPPA = 1.0


def _eval_signature(sig, rng: np.random.Generator) -> dict:
    cliff = check_clifford(sig)
    ratios = []
    spin_res = []
    mixed_res = []
    q_alt_res = []
    jsqs = []
    factors = []
    for _ in range(N_SAMPLES):
        psi = random_spinor(rng)
        out = contact_lagrangian(sig, psi, kappa=KAPPA)
        ratios.append(out["ratio"])
        spin_res.append(out["spin_residual"])
        mixed_res.append(out["mixed_residual"])
        q_alt_res.append(out["Q_alt_residual"])
        jsqs.append(out["j5_sq"])
        factors.append(out["spin_dual_factor"])
    ratios = np.array(ratios, dtype=np.float64)
    psi0 = random_spinor(rng)
    mut_palatini = mutated_ratio(
        sig, psi0, palatini_factor=1.0, source="dirac", kappa=KAPPA
    )
    mut_spin = mutated_ratio(
        sig, psi0, palatini_factor=0.5, source="plus_quarter_dual", kappa=KAPPA
    )
    return {
        "signature": sig.name,
        "clifford": cliff,
        "n_samples": N_SAMPLES,
        "ratios": ratios.tolist(),
        "ratio_mean": float(np.mean(ratios)),
        "ratio_std": float(np.std(ratios)),
        "ratio_max_abs_dev_from_mean": float(np.max(np.abs(ratios - np.mean(ratios)))),
        "spin_residual_max": float(np.max(spin_res)),
        "spin_dual_factor_mean": float(np.mean(factors)),
        "mixed_residual_max": float(np.max(mixed_res)),
        "Q_alt_residual_max": float(np.max(q_alt_res)),
        "j5_sq_samples": jsqs,
        "mutation_double_palatini_ratio": float(mut_palatini),
        "mutation_flip_spin_sign_ratio": float(mut_spin),
    }


def _verdicts(minus: dict, plus: dict) -> dict:
    def near(x: float, y: float) -> bool:
        return abs(x - y) < TOL

    c1_minus = near(minus["ratio_mean"], TARGET) and minus["ratio_max_abs_dev_from_mean"] < TOL
    c1_plus = near(plus["ratio_mean"], TARGET) and plus["ratio_max_abs_dev_from_mean"] < TOL
    c2 = minus["spin_residual_max"] < TOL and plus["spin_residual_max"] < TOL
    c3 = minus["mixed_residual_max"] < TOL and plus["mixed_residual_max"] < TOL
    c4 = abs(minus["ratio_mean"] - plus["ratio_mean"]) < TOL
    # C5: double Palatini (1/kappa instead of 1/(2kappa)) should give 3/8
    #     flip spin sign should give -3/16
    c5_pal_ok = (
        abs(minus["mutation_double_palatini_ratio"] - TARGET) > 0.05
        and abs(plus["mutation_double_palatini_ratio"] - TARGET) > 0.05
    )
    c5_sgn_ok = (
        abs(minus["mutation_flip_spin_sign_ratio"] - TARGET) > 0.05
        and abs(plus["mutation_flip_spin_sign_ratio"] - TARGET) > 0.05
    )
    return {
        "C1_PRIMARY": {
            "pass": bool(c1_minus and c1_plus),
            "minus_ratio": minus["ratio_mean"],
            "plus_ratio": plus["ratio_mean"],
            "target": TARGET,
        },
        "C2_spin_identification": {
            "pass": bool(c2),
            "minus_residual": minus["spin_residual_max"],
            "plus_residual": plus["spin_residual_max"],
            "minus_dual_factor": minus["spin_dual_factor_mean"],
            "plus_dual_factor": plus["spin_dual_factor_mean"],
            "paper_factor": -0.25,
        },
        "C3_no_mixed_symmetry": {
            "pass": bool(c3),
            "minus_residual": minus["mixed_residual_max"],
            "plus_residual": plus["mixed_residual_max"],
        },
        "C4_signature_independence": {
            "pass": bool(c4),
            "delta": abs(minus["ratio_mean"] - plus["ratio_mean"]),
        },
        "C5_mutations_fail": {
            "pass": bool(c5_pal_ok and c5_sgn_ok),
            "double_palatini_minus": minus["mutation_double_palatini_ratio"],
            "double_palatini_plus": plus["mutation_double_palatini_ratio"],
            "flip_spin_minus": minus["mutation_flip_spin_sign_ratio"],
            "flip_spin_plus": plus["mutation_flip_spin_sign_ratio"],
        },
    }


def main() -> int:
    rng = np.random.default_rng(SEED)
    minus = _eval_signature(signature_mostly_minus(), rng)
    plus = _eval_signature(signature_mostly_plus(), rng)
    verdicts = _verdicts(minus, plus)
    payload = {
        "task": "m9.1",
        "seed": SEED,
        "kappa": KAPPA,
        "n_samples_per_signature": N_SAMPLES,
        "tolerance": TOL,
        "target_ratio_compared_after_extraction": TARGET,
        "mostly_minus": minus,
        "mostly_plus": plus,
        "verdicts": verdicts,
        "gate_pass": bool(verdicts["C1_PRIMARY"]["pass"]),
    }
    os.makedirs(DATA, exist_ok=True)
    out_path = os.path.join(DATA, "m9_1_hehl_datta.json")
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
        fh.write("\n")
    print(f"wrote <repo>/openwave/xperiments/m9_emergent_gravity/research/data/m9_1_hehl_datta.json")
    print(f"mostly_minus ratio mean = {minus['ratio_mean']!r}")
    print(f"mostly_plus  ratio mean = {plus['ratio_mean']!r}")
    print(f"target (post-extraction comparison) = {TARGET!r}")
    for key, val in verdicts.items():
        flag = "PASS" if val["pass"] else "FAIL"
        print(f"  {key}: {flag}  {val}")
    print(f"GATE: {'PASS' if payload['gate_pass'] else 'FAIL'}")
    return 0 if payload["gate_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
