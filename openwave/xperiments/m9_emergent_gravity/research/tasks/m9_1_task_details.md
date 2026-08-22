# M9.1: CERTIFICATION GATE, independent elimination of Einstein-Cartan torsion

> Proposed column (ID pending maintainer assignment). Spec of record:
> [`../m9_theory_canonical.md`](../m9_theory_canonical.md) section 2. Source of
> the claim: McGwier Paper III equation (6) and the 1971 Hehl-Datta identity it
> invokes. This task certifies the *effective-interaction algebra*, not the
> holographic first-law argument.

## TASK PLANNING (2026-08-15)

### Scope

Independent algebraic verification of the one interaction the Emergent Gravity
program places beyond the Standard Model plus Einstein gravity: after solving
the Cartan equation for a Dirac field and substituting back, the Palatini plus
Hermitian Dirac action reduces to Riemannian Einstein-Hilbert plus

\[
\mathcal{L}_{\mathrm{HD}} = -\frac{3\kappa}{16}\, J_5^\mu J_{5\mu},
\qquad
J_5^\mu = \bar\psi\gamma^5\gamma^\mu\psi,
\qquad
\kappa = 8\pi G.
\]

This is the column's certification gate. Either outcome is a result. A match
does **not** certify Faulkner-Guica-Hartman-Myers-Van Raamsdonk, the finite-ball
Condition NL, the residue \(I_B\), or any de Sitter claim. Those stay untested.

### The independence protocol

| Role | Sees | Does not treat as input |
| --- | --- | --- |
| Designer (this session) | the author's Paper III and Final Status table | n/a |
| Solver | the Palatini 4-form, the Hermitian Dirac action, \(\kappa=8\pi G\), \(\varepsilon^{0123}=+1\), and the two metric signatures | the rational \(3/16\) is not coded into the contraction; it is compared only after the coefficient is extracted |
| Auditor | the same action, a *different* gamma representation and a *different* index contraction | the solver's intermediate tensors; it may read the solver's extracted number only to try to refute it |

No-search rule: the extracted coefficient is reported first; comparison to
\(3/16\) is a post-step. If the number lands elsewhere, that is the result.

### Pre-registered claims under test (fixed BEFORE the scripts ran)

Tolerance: machine residual \(\lvert\Delta\rvert < 10^{-12}\) on dimensionless
ratios, across both metric signatures.

| ID | Claim | Pass | Fail |
| --- | --- | --- | --- |
| C1 (PRIMARY, the gate) | On-shell elimination of algebraic torsion from Palatini + Hermitian Dirac yields \(\mathcal{L}_{\mathrm{int}} / (-\kappa J_5\cdot J_5) = 3/16\) | extracted ratio equals \(3/16\) within \(10^{-12}\) in both signatures | any other simple rational (\(0\), \(3/8\), \(3/32\), \(1/2\), \(-3/16\)) |
| C2 | The canonical Dirac spin tensor is totally antisymmetric and equals \(s^{\lambda\mu\nu} = -\frac{1}{4}\varepsilon^{\lambda\mu\nu\rho} J_{5\rho}\) | the two 24-component tensors agree within \(10^{-12}\) relative to \(\lVert J_5\rVert\) | a vector-trace piece, a mixed-symmetry piece, or a different rational prefactor |
| C3 | Minimal Dirac produces no mixed-symmetry spin source (algebraic input to the author's modular selection rule, not the holographic theorem) | the Young projection of \(s\) onto mixed-symmetry vanishes at \(10^{-12}\) | a nonzero mixed-symmetry component |
| C4 | The gate ratio is independent of metric signature once \(J_5\cdot J_5\) is formed with the same \(\eta\) | mostly-plus and mostly-minus ratios agree to \(10^{-12}\) | signature-dependent coefficient |
| C5 (mutation) | The C1 check can fail: doubling Palatini \(\frac{1}{2\kappa}\to\frac{1}{\kappa}\), or flipping \(s\to +\frac{1}{4}\varepsilon J_5\), must move the ratio off \(3/16\) by the predicted factor | each mutation fails C1 in the predicted direction | a mutation still reports \(3/16\) (the check is tautological) |

### Definition of done

| # | Item |
| --- | --- |
| 1 | Solver writes the extracted coefficient, the five claim residuals, and mutation outcomes to `data/m9_1_hehl_datta.json` |
| 2 | Adversarial auditor recomputes with its own method and files per-claim verdicts in `data/m9_1_audit.json` |
| 3 | Method note `findings/m9_1_hehl_datta_note.md` (equations first, equation-to-code map, audit record) |
| 4 | Briefing and application draft name this task as the proposed first cell, without moving any `MODELS.md` icon |

### Blindspot pass

| Blindspot | Mitigation |
| --- | --- |
| Palatini 4-form vs \(\frac{1}{2\kappa}eR\) off by 2 | derive both conversions in code; C5 mutates the factor |
| \(\gamma^5\) definition (\(i\gamma^0\gamma^1\gamma^2\gamma^3\) vs \(\gamma^0\gamma^1\gamma^2\gamma^3\)) flips \(\varepsilon\) | lock \(\gamma^5 := i\gamma^0\gamma^1\gamma^2\gamma^3\) in mostly-minus and the Clifford-consistent analogue in mostly-plus; auditor uses a different representation |
| Hermitian vs non-Hermitian Dirac doubles the spin | solver uses the Hermitian kinetic term written in the canonical; auditor repeats from \(\delta S/\delta\omega\) |
| \(J_5\cdot J_5\) sign convention (mostly-plus timelike vs spacelike) | report the dimensionless ratio \(\mathcal{L}/(-\kappa J\cdot J)\), not a signed energy story |
| Completing the square drops the Dirac-linear piece or the Palatini-quadratic piece | both terms are evaluated on-shell and off-shell (random \(K\)) as a consistency check |
| Treating a literature \(3/16\) as an input | the contraction returns a float; \(3/16\) appears only in the comparison block |

### Sub-experiments

| ID | What | Artifact |
| --- | --- | --- |
| S1 | Solver: gamma-matrix variation + Riemann quadratic + on-shell substitution | `scripts/hehl_datta.py`, `scripts/m9_1_hehl_datta_elimination.py` |
| S2 | Auditor: Hodge dual + Levi-Civita identities + 4-form Palatini, different representation | `scripts/m9_1_audit_hehl_datta.py` |
| S3 | Method note + briefing/application | `findings/m9_1_hehl_datta_note.md`, `__M9_model_briefing.md`, `APPLICATION.md` |

### Not computed

FGHMV equivalence, Condition NL, \(I_B\), Hehl-Datta *from relative entropy*,
de Sitter, nonlinear Einstein-Cartan, collider or astrophysical bounds, the
high-density bounce, any `MODELS.md` gravity-row measurement on a lattice.

## DEVIATIONS LOG

Cartan was pre-registered as the algebraic map \(T=\kappa s\). The solver
instead finds \(K\) by stationarity of \(\mathcal{L}_g+\mathcal{L}_D\)
(finite-difference Hessian). The printed \(T=\kappa s\) map is not used.
Reason: C2 failed (measured \(s=-\frac12(*J_5)\)), so inserting
\(T=\kappa s\) with the paper's \(s\) would have silently mixed two
normalizations. Stationarity does not need that map. The gate (C1) is
unchanged.

## FINDINGS

Full record: [`../findings/m9_1_hehl_datta_note.md`](../findings/m9_1_hehl_datta_note.md).

| ID | Finding |
| --- | --- |
| F1 | GATE PASSES. C1: \(r=3/16\) in both signatures, scatter \(<10^{-15}\). Independent complete-the-square audit (Weyl Clifford, no solver import) reproduces \(0.18750000000000006\) |
| F2 | C2 FAILS as written. Dirac and Weyl both give \(s=-\frac12(*J_5)\), not the paper's \(-\frac14\). Convention: \(2\delta L/\delta\omega\) vs \(\delta L/\delta\omega\). Not folded into a silent redefinition |
| F3 | C3 and C4 pass. Minimal Dirac spin is a 3-form. The ratio is signature-independent in the solver; the auditor left C4 QUALIFIED (one Clifford only) |
| F4 | C5 passes. Double Palatini \(\to 3/32\). Paper \(\alpha\) used as if it were this \(s\) \(\to 3/64\) |
| F5 | `palatini_quadratic_alt` was \(-Q\) from an index swap. Fixed: residual \(7\times 10^{-15}\). Not used in \(r\). \(3/16\) unchanged |

Artifacts: `scripts/hehl_datta.py`, `m9_1_hehl_datta_elimination.py`,
`m9_1_audit_hehl_datta.py`; `data/m9_1_hehl_datta.json`, `data/m9_1_audit.json`;
`findings/m9_1_hehl_datta_note.md`.
