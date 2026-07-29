#!/usr/bin/env python3
"""
m8_1_1_coexact_audit.py  --  ADVERSARIAL AUDIT of solverA/m8_1_1_coexact_solver.py

Deliberately different computational path from the audited script:

  audited script (A)                     this audit (B)
  -------------------------------------  ---------------------------------------------
  float64 matrices, element identity      exact 2x2 matrices over GF(p), p = k*L+1 with
  decided by |dM| < 1e-8                  L = lcm(1..24, 60); identity is dict equality,
                                          NO tolerance anywhere in the group construction
  i, sqrt2, sqrt5, phi as float           i = zeta_4, sqrt2 = zeta_8+zeta_8^-1,
                                          sqrt5 = Gauss sum over zeta_5, phi exact in GF(p)
  characters via Burnside-Dixon class-     characters via ISOTYPIC SPLITTING of V_a by the
  multiplication matrices + mpmath        class sums acting on the module (exact GF(p)
  50-digit floating eigenvectors          eigenspaces), multiplicity m from the EXACT
                                          commutant dimension, then a Dixon Fourier lift
                                          to exact values in Z[zeta_M].  No a_ijk, no eig.
  branching by 50-digit inner products    branching two ways: (i) module multiplicities from
                                          the isotypic split, (ii) exact Z[zeta_M] inner
                                          products.  Both must agree as integers.
  invariant dims by SVD rank, tol 1e-8    invariant dims by EXACT GF(p) rank of the full
                                          unreduced averaging projector + exact idempotency
  a <= 14, m <= 12, 18 groups             a <= 32, m <= 30, 34 groups (C_1..C_20,
                                          BD_2..BD_12, 2T, 2O, 2I) + a boundary probe to
                                          BD_16 / C_32

Run:  python3 m8_1_1_coexact_audit.py      (writes m8_1_1_coexact_audit.json)

Modules part1/part2/core/driver/proj live beside this file and hold the implementation;
sweep.py, compare.py, compare2.py, patterns.py, patterns2.py, probe.py, mutate.py,
runmuts.py, runmuts2.py, finalchecks.py are the individual audit stages.  This file
re-runs every stage and assembles the single results document.
"""
import json, subprocess, sys, os, time

HERE = os.path.dirname(os.path.abspath(__file__))
STAGES = [
    ("group_construction_and_orders", "part1.py",     {}),
    ("exact_arithmetic_selftests",    "part2.py",     {}),
    ("main_sweep_34_groups",          "sweep.py",     {}),
    ("compare_naive_signature_match", "compare.py",   {}),
    ("compare_rigorous_relabelling",  "compare2.py",  {}),
    ("pattern_attack",                "patterns.py",  {}),
    ("general_rule_check",            "patterns2.py", {}),
    ("truncation_boundary_probe",     "probe.py",     {"AUDIT_LMAX": "64"}),
    ("mutation_tests_round1",         "runmuts.py",   {}),
    ("mutation_tests_round2",         "runmuts2.py",  {}),
    ("mutation_tests_round3",         "runmuts3.py",  {}),
    ("independent_float_crosschecks", "finalchecks.py", {}),
    ("verdict_table",                 "mkverdicts.py",{}),
    ("defect_list",                   "mkdefects.py", {}),
]

def main():
    logs = {}
    for tag, script, env in STAGES:
        e = dict(os.environ); e.update(env)
        t = time.time()
        r = subprocess.run([sys.executable, os.path.join(HERE, script)],
                           capture_output=True, text=True, cwd=HERE, env=e)
        logs[tag] = {"script": script, "returncode": r.returncode,
                     "seconds": round(time.time() - t, 1),
                     "stdout_tail": r.stdout[-4000:], "stderr_tail": r.stderr[-1500:]}
        print("[%-34s] rc=%d  %.1fs" % (tag, r.returncode, time.time() - t))

    L = lambda f: json.load(open(os.path.join(HERE, f)))
    doc = {
        "audit_of": "solverA/m8_1_1_coexact_solver.py + solverA/m8_1_1_coexact.json",
        "method": json.loads(json.dumps({
            "field": "GF(p), p = 10708457761 = 2*lcm(1..24,60)+1; all roots of unity, "
                     "sqrt2, sqrt5, phi exact in GF(p). Boundary probe uses "
                     "L = lcm(1..64), p = 50837476016421357502988548801.",
            "group_construction": "exact closure under right multiplication, identity by "
                                  "dict equality (no tolerance, no epsilon)",
            "irreducible_characters": "isotypic splitting of V_a by class sums acting on the "
                                      "module (exact GF(p) eigenspaces), multiplicities from "
                                      "the exact commutant dimension, Dixon Fourier lift to "
                                      "exact Z[zeta_M]; class-multiplication matrices never built",
            "inner_products": "exact in Z[zeta_M], reduced mod Phi_M, Fraction results",
            "invariant_dimensions": "three exact routes: character sum convention A, "
                                    "convention B, and the branching formula "
                                    "(m-1)*mult(tau* in V_m)+(m+1)*mult(tau* in V_(m-2)); "
                                    "plus exact GF(p) rank of the full unreduced projector",
            "ranges": {"symmetric_powers_a": [0, 32], "levels_m": [2, 30],
                       "groups": "C_1..C_20, BD_2..BD_12, 2T, 2O, 2I (34 groups)"},
        })),
        "stage_logs": logs,
        "my_results_per_group": L("_sweep_raw.json"),
        "exact_projector_rank_checks": L("_projchecks.json"),
        "comparison_naive_signature_match": L("_compare.json"),
        "comparison_rigorous_relabelling": L("_compare2.json"),
        "pattern_analysis": L("_patterns.json"),
        "mutation_tests": L("_mutations.json") + L("_mutations2.json") + L("_mutations3.json"),
        "defects_found": L("_defects.json"),
        "verdicts": L("_verdicts.json"),
    }
    with open(os.path.join(HERE, "m8_1_1_coexact_audit.json"), "w") as f:
        json.dump(doc, f, indent=1)
    print("wrote m8_1_1_coexact_audit.json (%d bytes)"
          % os.path.getsize(os.path.join(HERE, "m8_1_1_coexact_audit.json")))

if __name__ == "__main__":
    main()
