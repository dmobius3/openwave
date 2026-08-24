# Solver-versus-discretization localization test

> Run 2026-08-24 as a standalone diagnostic whose outcome CANNOT alter P1A. P1A is closed as
> filed. Reproducer: `localization_test.py` beside this note. This answers one causal question with no fit, no threshold, no envelope, and no sector
> rulings.

## The question

P1A.4b's corrected floor did not reclassify anything, leaving `R_1` and `R_3` at about `1e-10`
above a `5.4e-12` floor. A commissioner hypothesis proposed that the missing ingredient was the
`O(N)` accumulation of an `N x N` eigensolve, since

```text
R_1: |Im| = 1.30e-10   vs   N eps ||L|| = 9.2e-11    ratio 1.4
R_3: |Im| = 2.76e-10   vs   N eps ||L|| = 1.4e-10    ratio 2.0
```

**That hypothesis was wrong to act on and is now disproved.** It substituted a worst-case bound
for a MEASURED a posteriori residual of about `3e-12`, on the grounds that it happened to cover
the two sectors in question. It also borrowed the cluster spectral-projector norm to bound
INDIVIDUAL eigenvalue error, which P1A.3 never licensed: a nearly degenerate cluster can have a
beautifully conditioned invariant subspace while the individual eigenvalues inside it are far
more sensitive.

## The experiment

Hold the assembled floating-point `L_h` byte-for-byte fixed. Change ONLY the arithmetic precision
of the eigensolve. No other variable moves.

| sector | `N` | cluster `λ` | float64 (~16 dps) | mpmath 30 dps | mpmath 50 dps |
| --- | --- | --- | --- | --- | --- |
| `R_1` | 40 | 3 | `1.021e-10` | `1.021e-10` | `1.021e-10` |
| `R_3` | 60 | 8 | `1.942e-10` | `1.942e-10` | `1.942e-10` |
| `R_2` | 40 | 63 | `2.222e-01` | `2.222e-01` | `2.222e-01` |

**Identical to four significant figures across three precisions spanning sixteen to fifty digits.**

## The conclusion

The imaginary parts are **not eigensolver arithmetic**. They are eigenvalues of the discrete
matrix `L_h` to that scale. Tripling the working precision moves them not at all.

Three consequences follow, and none of them is discretionary:

1. **The `N eps ‖L‖` agreement was a coincidence.** Two ratios near 1.4 and 2.0 on a quantity
   that does not depend on precision at all. Curve fit, not mechanism.

2. **No floor repair can legitimately reclassify `R_1` and `R_3`.** They are genuinely
   ABOVE_FLOOR, because their imaginary parts are properties of the matrix rather than of the
   solver. Arm A must not classify them away regardless of what empirical formula happens to
   cover them. **P1A.4b's global FAIL is therefore correct and is not an artefact of a defective
   floor.** A fourth Arm-A pass would have been chasing something that is not there.

3. **The source lies upstream of the eigensolver, in the assembly or the discretization.** `L_h`
   has genuinely complex eigenvalues, running from about `1e-10` on the first cluster to order
   `10` at mid and high modes, consistent with the measured `eps_SA` of roughly 0.45 being
   dominated by the high end. The continuum operator is self-adjoint and this discrete one is
   not. **The precision test does not say WHICH upstream stage is responsible**: it separates
   solver from matrix, not equivariant assembly from the underlying RBF-FD discretization. That
   attribution is a further question and is not answered here.

## What this means for M8.4

The contamination qualification failed globally under a rule frozen before the fit, and this test
shows the failure is substantive rather than procedural. On the frozen outcome rules that is
estimator-qualification failure and **no scientific label**, which was always a legitimate
preregistered result.

Continuing on this substrate for M8.4 scoring would require a structural discretization
correction, independently justified, that makes `L_h` self-adjoint in the physical metric to a
degree the target observable can tolerate. That is a different piece of work from anything P1A
was chartered to do, and it should not be started as a repair to P1A.

What survives untouched: P0 PASS, P1A.0 to P1A.3 PASS with the invariant-subspace estimator
qualified, the fixed-`k` refinement family PASS, and the cloud admissibility gate PASS and
demonstrably discriminating.
