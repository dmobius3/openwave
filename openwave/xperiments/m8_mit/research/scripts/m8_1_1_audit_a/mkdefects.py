import json
D=[
{"id":"D1","severity":"robustness only, changes no reported number",
 "where":"m8_1_1_coexact_solver.py, isotypic-block irrep extraction, `assert rk == ds`",
 "finding":"The assertion checks only the RANK of the isotypic projector, so it cannot tell sigma "
           "from sigma*. Mutation M2 removed the complex conjugate in "
           "`cs = conj(chartab[s][...])`, which projects onto the sigma* block instead: the "
           "assertion stayed green and the script completed normally. The wrong-representation "
           "matrices were caught only by `irrep_matrix_character_error`, which jumped from "
           "5.6e-16 to 1.7320508 for 2T. That quantity is written into the JSON but is never "
           "compared to any threshold anywhere in the script.",
 "effect_on_numbers":"none: in the real run irrep_matrix_character_error = 1.4e-15 at worst, "
                     "and my independent recomputation reproduces every reported number."},
{"id":"D2","severity":"wording only",
 "where":"T3, `agree: (va == vb == rk)` and the reported pair char_tau / char_tau_conj",
 "finding":"The spec asked for mu computed 'BOTH with chi_tau and with its complex conjugate "
           "(the tau versus tau* convention)'. Because chi_(E_m) is real for every m and every "
           "Gamma (E_m is a sum of restricted SU(2) symmetric powers, which are self-dual), the "
           "two conventions are FORCED to be equal, so `va == vb` is a check that cannot fail for "
           "the E_m actually used. Mutation M6 swapped the two conventions outright: every "
           "reported number was unchanged and no check went red. Mutation M5 replaced chi_(E_m) "
           "by a non-self-dual class function and the check DID go red, confirming that the only "
           "thing it tests is the reality of chi_(E_m), not the tau/tau* convention.",
 "effect_on_numbers":"none. Convention A is the one the spec's definition of mu (invariants of "
                     "E_m tensor V_tau) calls for, and that is the one `first_nonzero_m` uses, "
                     "so q, q^2 and e are computed from the right side."},
{"id":"D3","severity":"wording only",
 "where":"`G.closed = True` with the comment 'closure is verified by idx() succeeding'",
 "finding":"`closed_under_multiplication` in the output JSON is a hard-coded literal True; there "
           "is no code path that can set it False. Mutation M4 dropped one element from the "
           "closed group: the script raised RuntimeError('element not found, dist 0.707107') "
           "rather than reporting closed=False. The property IS enforced (by the exception), but "
           "the JSON field is not evidence of it.",
 "effect_on_numbers":"none"},
{"id":"D4","severity":"robustness only",
 "where":"irrep_matrix_character_error / _homomorphism_error / _unitarity_error / "
         "sym_power_max_homomorphism_error",
 "finding":"All are computed, stored and serialized, but never asserted. Mutation M13 multiplied "
           "the extracted irrep matrices by i: irrep_char_err went to 4.24 and irrep_hom_err to "
           "1.41 and the script still completed with every boolean check green. In addition the "
           "homomorphism errors are sampled on a stride (range(0, N, N//6) and N//5), not over all "
           "pairs, so a homomorphism failure confined to unsampled pairs would go unseen. My audit "
           "verifies the homomorphism exactly for all pairs over GF(p).",
 "effect_on_numbers":"none; the real values are ~1e-15"},
{"id":"D5","severity":"scope limit, honestly reported as null but worth naming",
 "where":"JMAX = 14 and MMIN, MMAX = 2, 12",
 "finding":"For the 18 groups the spec lists these bounds have headroom (largest least_a and "
           "largest e are both 7, at 2I). One step outside that list the bounds bite: I ran agent "
           "A's own analyse_group on wider members of the same families. C_26 -> T6 e is null for "
           "the irreducible at d = 13. BD_13 -> e null at d = 13. C_30 and BD_15 -> T1 least_a "
           "null at d = 15. BD_16 -> hard failure, RuntimeError('BD_16: no multiplicity-one host "
           "for sigma 10'), because no V_a with a <= 14 contains that 2-dim irreducible. The "
           "results for the spec's own groups are unaffected; the point is that the ranges are "
           "not stated as a limitation and the transition to null (and then to a crash) is one "
           "group away.",
 "effect_on_numbers":"none within the spec's 18 groups"},
{"id":"D6","severity":"wording only",
 "where":"T2 block, `parity_rule_holds`",
 "finding":"For C_1, C_3, C_5, C_7, C_9 (the groups without -I) the JSON carries "
           "`parity_rule_holds: false` with n_violations of 7, 20, 29, 36, 39. The adjacent field "
           "`rule_applicable: false` says the spec never asked for it there, but a reader "
           "grepping the output for a false parity flag sees five apparent failures that are not "
           "failures. My run reproduces the same violation counts.",
 "effect_on_numbers":"none"},
{"id":"D7","severity":"limitation of one check, not a defect in the result",
 "where":"`gram_matches = (Gram_exact == M M^T)`",
 "finding":"This is the strongest exact check in the script and it does fire: mutation M1 "
           "(corrupting one multiplicity) turned it red. But M M^T is invariant under ANY "
           "permutation of the columns of M, i.e. under relabelling the irreducibles, so the "
           "check by construction cannot validate which sigma is which, and therefore cannot "
           "validate the distance vector, the T1 least_a assignment or the A entries. Mutation "
           "M14 swapped two columns and gram_matches stayed green. In agent A's full pipeline a "
           "relabelling is caught downstream (M14b, on 2O, tripped `assert rk == ds`), so this is "
           "a scope note about the check, not a hole in the result.",
 "effect_on_numbers":"none"},
{"id":"D8","severity":"spec-compliance gap, minor",
 "where":"character table serialization, `exact_rational`",
 "finding":"The spec says 'Exact rationals as strings like \"73/144\", never as floats'. 534 of "
           "834 character-table entries carry an exact rational. The other 300 are null: 256 are "
           "genuinely complex and 44 are real but irrational (sqrt2 in BD_4, the golden ratio in "
           "BD_5 / BD_6 / 2I). For those 44 the only record in the JSON is a 30-digit decimal, so "
           "nothing in the file certifies that -1.4142135623730950488 is exactly -sqrt2. This "
           "audit's JSON gives every character value as exact Z[zeta_M] coefficients so the "
           "algebraic identity is checkable.",
 "effect_on_numbers":"none; the values are correct"},
]
json.dump(D,open("_defects.json","w"),indent=1)
for d in D: print("%s [%s] %s"%(d["id"],d["severity"],d["where"]))
