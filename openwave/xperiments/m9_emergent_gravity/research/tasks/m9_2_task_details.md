# M9.2: PRE-REGISTRATION, Newton limit (attractive \(1/r^2\))

> Proposed column (ID pending maintainer assignment). Spec of record:
> [`../m9_theory_canonical.md`](../m9_theory_canonical.md). This page **locks
> the gates on 2026-08-15, before any Newton script exists**. Later changes go
> in a dated addendum at the end, never in-place. Official M9.2 is this
> Newton-limit task. The files `m9_2_ib_hadamard.py`, `m9_3_ib_analytic.py`,
> `m9_4_ib_hadamard_complete.py`, and `m9_5_ec_symplectic.py` are campaign
> scripts for \(I_B\) and the axial obstruction; they are **not** this task.

## TASK PLANNING (2026-08-15)

### Scope

Script-verify the Newtonian limit of the *metric* sector of the NSM action:
linearized Einstein-Hilbert about Minkowski, sourced by a static compact
rest-mass density, produces an attractive \(1/r^2\) acceleration. This is the
[`MODELS.md`](../../../../../MODELS.md) simplest test for
**Gravity: Newton limit (GEM)** read as "attractive \(1/r^2\) between masses".
The GEM *route* (M5 boost-tilt) is **not** this column's mechanism and is
out of scope.

A match does **not** certify Faulkner-Guica-Hartman-Myers-Van Raamsdonk,
Condition NL, \(I_B\), de Sitter, nonlinear Einstein-Cartan, or the modular
selection of Einstein-Cartan. Those stay untested. Torsion is out of scope:
it is algebraic, vacuum-vanishing, and does not propagate, so a static
spinless source has \(\omega=\mathring{\omega}(e)\).

This task does **not** move a `MODELS.md` icon until a script, an audit, and
a method note exist and the pre-registered gates below have been scored.

### Why this is inherited Einstein, and why a script still earns a cell

The NSM gravitational sector *is* Palatini Einstein-Cartan. In vacuum, or
for a spinless source, that is Einstein-Hilbert. The Newton limit is
therefore inherited, not discovered. M8 left the same inheritance at 🚧
because it was not separately checked. This task is the check: a runnable
artifact that can fail.

### Pre-registered claims under test (fixed BEFORE any script is written)

Signature lock for this task: mostly minus,
\(\eta=\mathrm{diag}(-1,+1,+1,+1)\). Weak-field metric and Poisson equation:

\[
\mathrm{d}s^2
=
-(1+2\Phi)\,\mathrm{d}t^2
+
(1-2\Phi)\,\delta_{ij}\,\mathrm{d}x^i\mathrm{d}x^j,
\qquad
\nabla^2\Phi = 4\pi G\rho,
\qquad
\mathbf{a} = -\nabla\Phi.
\]

A single static compact source of mass
\(M=\int\rho\,\mathrm{d}^3x>0\), centered at the origin. The test observable
is the radial acceleration of a static probe at far-field points that do not
overlap the source. Implementation lock: a 3-dimensional discrete Poisson
solver on a cubic grid (not the continuum Green's function written in by
hand, and not a 2-dimensional Laplacian). \(G\) is a coded positive
constant. No holographic input.

Tolerance: relative residual \(\lvert\Delta\rvert<0.05\) on dimensionless
ratios at the three pre-registered probe radii, on at least two grids.

| ID | Claim | Pass | Fail |
| --- | --- | --- | --- |
| C1 (PRIMARY, the cell) | Far-field \(\lvert\mathbf{a}\rvert\) is attractive \(GM/r^2\) | at three probe radii \(r\in\{0.30 L,\,0.35 L,\,0.40 L\}\) on a box of half-width \(L\), with the source compact compared to \(0.30 L\): (i) \(\mathbf{a}\cdot\hat{r}<0\); (ii) \(\bigl\lvert \lvert\mathbf{a}\rvert r^2/(GM)-1\bigr\rvert<0.05\); (iii) a log-log fit of \(\lvert\mathbf{a}\rvert\) vs \(r\) on those three points has slope \(\alpha\) with \(\lvert\alpha+2\rvert<0.08\) | repulsive, or any of (ii)-(iii) missed on the finer of two grids |
| C2 | The potential itself is \(\Phi=-GM/r\) in the same far field | \(\bigl\lvert \Phi r/(GM)+1\bigr\rvert<0.05\) at the same three radii | a different power, the wrong sign, or a residual above the cut |
| C3 (vacuum) | \(\rho=0\) produces no Newtonian field | \(\max\lvert\mathbf{a}\rvert\) and \(\max\lvert\Phi\rvert\) on the probe set are \(<10^{-8}\) times the C1 values at the same points | a ghost \(1/r^2\) or a nonzero constant field |
| C4 (mutation, 2-d) | Replacing the 3-d Laplacian by a 2-d Laplacian in the same plane must fail C1 | the 2-d run reports \(\lvert\alpha+1\rvert<0.15\) (logarithmic potential, \(1/r\) force) and fails C1 (ii) | the 2-d mutation still reports C1 PASS |
| C5 (mutation, \(G\)) | \(G\to -G\) must reverse attraction | C1 (i) flips sign; C1 (ii) still holds with \(\lvert G\rvert\) | the flipped-\(G\) run still reports inward acceleration |

### Definition of done

| # | Item |
| --- | --- |
| 1 | Solver writes C1-C5 residuals to `data/m9_2_newton.json`. Suggested name: `scripts/m9_2_newton_limit.py`. That file does not exist at lock time |
| 2 | Adversarial auditor recomputes with a different Poisson discretization or a different source shape and files per-claim verdicts in `data/m9_2_audit_newton.json` |
| 3 | Method note `findings/m9_2_newton_note.md` (equations first, equation-to-code map, audit record) |
| 4 | A `MODELS.md` Newton-limit cell is proposed only after 1-3, and only at the weight the gates support. This lock page does not move the cell |

### Blindspot pass

| Blindspot | Mitigation |
| --- | --- |
| Hand-inserting \(\Phi=-GM/r\) and calling it a measurement | C1-C2 must come from a discrete Poisson solve. The continuum Green's function is a check, not the primary observable |
| Periodic images on a cubic box fake the exponent | use Dirichlet \(\Phi=0\) on the boundary, keep probes well inside, and require two grid spacings |
| Self-force / source overlap | probes locked at \(r\ge 0.30 L\); source rms radius locked \(<0.05 L\) |
| Scoring GEM (M5 boost-tilt) as this cell | GEM is out of scope. The source is static and spinless |
| Scoring holography as this cell | no entanglement, no ball, no FGHMV input |
| Scoring torsion / HD as this cell | spinless source; \(\omega=\mathring{\omega}\) |
| Moving `MODELS.md` from this lock page | definition of done item 4. No icon without the script |

### Sub-experiments (not started)

| ID | What | Artifact |
| --- | --- | --- |
| S1 | 3-d Poisson + geodesic acceleration, C1-C3 | `scripts/m9_2_newton_limit.py` (absent at lock) |
| S2 | Mutations C4-C5 | same driver, flagged runs |
| S3 | Second-method audit | `scripts/m9_2_audit_newton.py` (absent at lock) |
| S4 | Method note | `findings/m9_2_newton_note.md` (absent at lock) |

### Not computed (and not to be smuggled in)

FGHMV equivalence, Condition NL, \(I_B\), Hehl-Datta from relative entropy,
de Sitter, nonlinear Einstein-Cartan, light bending, time dilation, \(\Lambda\),
GEM, collider or astrophysical bounds, any particle-row cell.

## DEVIATIONS LOG

None to the locked gates. Implementation note 2026-08-15:
grids \(n=65\) and \(n=97\) (both resolve \(R=0.04L\)).
A trial \(n=41\) put the whole source in one cell (rms \(0\))
and is not used. Acceleration is the discrete gradient
interpolated to the probe, not a derivative of interpolated
\(\Phi\).

## FINDINGS

Run 2026-08-15. C1 PRIMARY **PASS** (finer residuals
\(1.29\%\), \(1.04\%\), \(1.07\%\); \(\alpha=-2.008\)).
C2 **FAIL** (Dirichlet images, \(26\)--\(35\%\)).
C3 PASS. C4 PASS. C5 PASS. Auditor C1 CONFIRMED, C2
REFUTED. Note:
[`../findings/m9_2_newton_note.md`](../findings/m9_2_newton_note.md).
Status: **RUN. C1 PASS. C2 FAIL. No MODELS.md edit.**
