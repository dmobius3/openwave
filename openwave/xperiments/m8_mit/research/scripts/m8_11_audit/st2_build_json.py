"""Write audit_stage2.json: task-1 labels (with the comparison counts from out2/compare.json),
task-2 reproduction result, task-3 mutation outcomes (out2/mutations.json), task-5 grades,
and the consulted-material manifest (last key)."""
import json

cmp_ = json.load(open("out2/compare.json"))["summary"]
mut = json.load(open("out2/mutations.json"))
thm = json.load(open("out2/theorem_checks.json"))
exm = json.load(open("out2/exact_m.json"))

RULE_EXACT = "sympy: simplify(radsimp(mine - theirs)) == 0 (theirs squared where it reports a norm and I report a squared norm)"
items = {
    "0": ("CONFIRMED", "== on integers/booleans; exact sympy equality on |q|^2 and <3 3;3 -3|6 0>"),
    "1": ("CONFIRMED", "== on all 38 dim Hom entries (n = 0..18, both sectors) and the DERIVED declaration; "
                       "intertwiner residuals: both below the common 1e-40 requirement at 50 digits"),
    "2": ("CONFIRMED", RULE_EXACT + " for Q (12) and r6hat (6); == on the multiple-of-Phi and stationarity flags"),
    "3": ("CONFIRMED", "== on the block chi-subspace dimension and on the chi-subspace dimensions of V_J, J = 0..9 (60); "
                       "stabilizer orders; and st2_item3.py: the solver's closed-form generators applied with MY D-matrices "
                       "fix [u] with the claimed characters (residual <= 7.8e-16), generate groups of orders 48, 24, 10, 12, "
                       "and the character is a homomorphism on each whole group (<= 2.6e-15). Float-level on both sides "
                       "(neither worker certified the finite stabilizers exactly)."),
    "4": ("CONFIRMED", RULE_EXACT + " for every ||Pi_n xi||^2/g^2 (level 6: their explicit 0 vs my 'not listed, 0 by definition') "
                       "and for every ||Pi_n N||^2; == on the character flag"),
    "5": ("CONFIRMED", RULE_EXACT + " for along-Phi, perp norm (squared), e_t / i e_t / tau_x / tau_y components and remainders"),
    "6": ("CONFIRMED", RULE_EXACT + " for the forms and the r6hat second derivatives, including cross terms (same polarized reading)"),
    "7": ("CONFIRMED", RULE_EXACT + " for ||kappa||^2, components, remainders and the perp equation without kappa; "
                       "== on free-direction dimensions; with-kappa perp equation: mine exact 0, theirs '0 (residual < 1e-120 at 120 digits)'"),
    "8": ("CONFIRMED", RULE_EXACT + " for lambda4/g^2 with and without kappa; sign: both 'negative for every g != 0' "
                       "(their field is a literal in s99_assemble.py, but their values and argument support it)"),
    "9": ("PARTIAL", "set of computed zeros in items 4-8 mapped onto a common vocabulary: both set differences empty; "
                     "unresolved: 0 = 0. PARTIAL because one of their reasons is incomplete as stated (the with-kappa perpendicular "
                     "equation is 'zero by construction', which needs the solvability / Noether argument they give only in item 10) "
                     "and their machine check of reasons (s09) cannot fail (mutations m1, m2, m5)."),
    "10": ("CONFIRMED", "== on existence (yes for all six), orbit dimension in Phi-perp (2,2,3,3,3,3), rank 12 - orbit, gaps 72/32; "
                        "arguments compared in task 4 (both correct)"),
    "11": ("CONFIRMED", RULE_EXACT + " for r6hat and Re<Phi, DN e_t> in both sectors; == on stationarity (no) and Pi6N parallel (no); "
                        "their |grad|^2 = 13225/92928 equals the sum of squares of my two nonzero gradient components"),
    "12": ("CONFIRMED", "each side below its own declared requirement: mine 8.5e-15 < 1e-10 (float64), theirs 1.5e-29 < 1e-25 (30 digits); "
                        "levels 20, 22 empty on both sides"),
}
labels = {}
for k, (lab, rule) in items.items():
    s = cmp_.get(k, {})
    labels["item" + k] = {"label": lab, "comparison": rule, "values_compared": s.get("compared"), "values_agreeing": s.get("agree")}

grades = {
    "T.a": {"grade": "ESTABLISHED", "note": "levels from item 1 (none below 6; nearest 10 / 8; gaps 72 / 32). T.a's 'nearest level' "
                                            "wording alone does not exclude levels 0, 2 in the 4-dim sector; item 1 does."},
    "T.b": {"grade": "ESTABLISHED", "note": "F commutes with rotations and N(chi psi) = chi N(psi); xi in X_H checked exactly on my xi "
                                            "(U1-U4, U6 exact generators; U5 exact weight argument)."},
    "T.c": {"grade": "ESTABLISHED", "note": "holds with 'tangential part' taken in the moving tangent space at theta (the theta-gradient "
                                            "of the reduced functional); with the fixed projection Pi_T the a^3 division needs lambda "
                                            "substituted from the radial equation first; same Jacobian. The 1/4 comes from dQ = 4 Re<N, E>, "
                                            "the first-order form of the bridge."},
    "T.d": {"grade": "ESTABLISHED", "note": "IFT in (lambda, theta) plus the uniqueness of the range solution; uniqueness only within X_H."},
    "T.e": {"grade": "ESTABLISHED", "note": "needs the reparametrization a_T.c = |k| -> a = <Phi, psi> (analytic, odd, a(1 + O(a^4))); "
                                            "T's ansatz allows a block tilt at every odd order, the worklist's section 2 ansatz does not "
                                            "at order a^5; they agree through a^3."},
    "E1": {"grade": "ESTABLISHED", "note": "m = 1 checked exactly for exact generators at U1-U4 (st2_exact_m.py)."},
    "E2": {"grade": "ESTABLISHED", "note": "all stated facts checked (orbit direction = -(2/5)sqrt39 * i e_t exactly; no continuous "
                                           "symmetry of K_H at U6; chart actions; m = 2 exact at U6, by weights at U5). The T.c "
                                           "nondegeneracy that E2 uses but does not state: Hess_T Q = -224/165, -42/55 (U5); "
                                           "det 5770240/7913763, 202860/879307 (U6)."},
    "bridge": {"grade": "ESTABLISHED", "note": "derived; checked exactly in 12 cases, Q'' from an N-only quartic interpolation = 4 x form, "
                                               "cross terms by polarization."},
    "lambda4_lemma": {"grade": "DEFECT", "note": "the formula and the kappa-independence hold (2X + conj X checked exactly). The sign step "
                                                 "leaves its argument unsupplied and leaves a necessary hypothesis unstated: every non-block "
                                                 "level of the sector lies above 6 (mu_n > 0). T.a does not imply it for the 4-dim sector. "
                                                 "Rescued: item 1 gives mu_n > 0; Q != 1 => |Phi|^2 not constant => sum_{n!=6} "
                                                 "||Pi_n N||^2 = int|Phi|^6 - Q^2 > 0 (Cauchy-Schwarz). Minor: the sum must exclude n = 6."},
    "S_lemma": {"grade": "ESTABLISHED", "note": "section-level action derived: (S psi)(g) = conj(psi(r^-1 g)) R in the P-trick form "
                                                "(R_{m,m'} = (-1)^m delta_{m,-m'}), fibres x -> conj(D^J(r)) Theta x, columns b -> "
                                                "(-1)^b Theta(col -b); keeps the sector because conj(P) = R P R exactly. On the chart it "
                                                "acts as u(z) -> u(conj z) with no phase. All consequences hold exactly on my objects. "
                                                "Two remarks: the phase freedom must be taken in the lift (or in eta, which Phi and tau_x "
                                                "share), not by rephasing u alone; 'hence v_y = 0' also uses L_T(tau_y, tau_y) = 56/429, "
                                                "21/286 != 0, i.e. the T.c hypothesis."},
}

out = {
    "task1_labels": labels,
    "task2_reproduction": {"reproduces_own_results": True,
                           "detail": "./py solver_rerun/run_all.py from a clean copy: 503 PASS, 0 FAIL; results.json identical "
                                     "(recursive diff: 0 differences); all 11 results/*.json byte-identical; log identical line by "
                                     "line (964 lines) after masking timings"},
    "task3_mutations": {k: {"n_PASS": v["n_PASS"], "n_FAIL": v["n_FAIL"], "lines": v["lines"][:6]} for k, v in mut.items()},
    "task5_grades": grades,
    "task5_supporting_checks": {"theorem_checks": {"total": len(thm["log"]), "pass": sum(l["ok"] for l in thm["log"]),
                                                   "mutants_run": sum(l["mutant_caught"] is not None for l in thm["log"]),
                                                   "mutants_caught": sum(l["mutant_caught"] is True for l in thm["log"])},
                                "exact_m": exm, "bridge_values": thm["bridge"], "lambda4_values": thm["lambda4"]},
    "consulted_material_manifest": {
        "read_in_this_room": ["BRIEF_stage2.md", "BRIEF.md", "worklist.md (== solver_work/worklist.md byte for byte)", "METHOD.md",
                              "AUDIT_STAGE1.md", "audit_results.json", "THEOREM.md",
                              "solver_work/RETURN.md", "solver_work/results.json", "solver_work/run_all.log",
                              "solver_work/run_all.py", "solver_work/common.py", "solver_work/s00_group.py",
                              "solver_work/s01_levels.py", "solver_work/s02_su2_exact.py", "solver_work/s03_stabilizers.py",
                              "solver_work/s04_pipeline.py", "solver_work/s09_zeros.py", "solver_work/s10_existence.py",
                              "solver_work/s12_quadrature.py", "solver_work/s13_identification_B.py",
                              "solver_work/s14_negative_controls.py", "solver_work/s99_assemble.py",
                              "my stage-1 code: su2.py, engine.py, s_main.py, mq.py, numD.py, linalg_mq.py, s03_stabilizers.py",
                              "my stage-1 outputs: out/main_all.pkl, out/main_all.json, out/group.pkl, out/item3.json, out/item01.json",
                              "py, .room.sb (the room interpreter wrapper and sandbox profile)"],
        "outside_the_room": "nothing (no web, papers, books or other files)",
        "software": "Python 3.12 standard library, numpy, mpmath, sympy (the solver's scripts also use sympy.physics.quantum.cg "
                    "as their own cross-check; I did not call it)",
        "item1_declaration_stage1": "DERIVED (unchanged)",
    },
}
json.dump(out, open("audit_stage2.json", "w"), indent=1)
print("wrote audit_stage2.json")
