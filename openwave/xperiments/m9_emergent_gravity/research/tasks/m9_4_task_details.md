# M9.4: Minimal axial UV deformation (not a selected completion)

> Proposed column. Spec of record:
> [`../m9_theory_canonical.md`](../m9_theory_canonical.md). Paper 11 refused
> to invent a CFT to fill Q4. Paper 14 Thm. obs forces a *deformation* if
> axial CFT Fisher is to match a local bulk. This task tests that
> deformation only. Official M9.4 is this note. The file
> `m9_4_ib_hadamard_complete.py` is an \(I_B\) campaign script, not this
> task.

## TASK PLANNING (2026-08-15)

### Scope

Split Q4.

- **Q4a.** Select a microscopic pair (CFT, compactification, asymptotically
  safe fixed point) that produces the installed SM and a renormalizable
  metric sector. Still open by construction. Out of scope.
- **Q4b.** The minimal local bulk deformation that evades Paper 14
  Thm. obs in the axial channel and reduces to Hehl-Datta in the IR.

Q4b Lagrangian, axial sector only, tree level:

\[
\mathcal{L}_k
=
A\bigl(1+k^2/M^2\bigr)\lvert S(k)\rvert^2
+
B\, S(-k)\cdot J_5(k),
\]

with \(A,B\) locked so that \(M\to\infty\) recovers the M9.1 ratio
\(3/16\). Then

\[
r(k)
:=
\frac{\mathcal{L}_{\mathrm{eff}}(k)}{-\kappa\,J_5\cdot J_5}
=
\frac{3/16}{1+k^2/M^2}.
\]

A match of C1-C4 does **not** select the SM, does not renormalize
Einstein-Hilbert, and does not provide a CFT. Final Status "UV
completion" stays `[O]` for Q4a.

### Pre-registered claims

| ID | Claim | Pass | Fail |
| --- | --- | --- | --- |
| C1 | \(r(0)=3/16\), and \(M\to\infty\) recovers \(3/16\) at finite \(k\) | residuals \(<10^{-12}\) | any other IR rational |
| C2 | \(r(M)=3/32\); \(r(k)\) strictly decreases in \(k^2\) | both | a contact (flat \(r\)) |
| C3 | position-space kernel is 4d Yukawa, not \(\delta^{(4)}\) | 4d radial ODE residual \(<5\cdot10^{-4}\) | a numerical delta |
| C4 | \(M\to 0\) sends \(r(k\to\mathrm{fixed}\neq 0)\to 0\) | true | a residual contact |
| C5 (mutation) | locking the IR target to \(3/8\) leaves \(r(0)=3/8\) | the check can fail | tautology |

### Definition of done

| # | Item |
| --- | --- |
| 1 | Solver writes C1-C5 to `data/m9_4_uv_axial.json` |
| 2 | Auditor, no solver import, files `data/m9_4_audit_uv.json` |
| 3 | Method note and Paper 15. No `MODELS.md` edit |
| 4 | Final Status UV label for Q4a is not moved to `[P]` |

### Not computed

A CFT with \(SU(3)\times SU(2)\times U(1)\). A string compactification.
Asymptotic safety. Metric renormalizability. de Sitter. Loop corrections.
Ghost-freedom of a general Poincaré-gauge Lagrangian (literature only).

## DEVIATIONS LOG

None.

## FINDINGS

Full record:
[`../findings/m9_4_uv_deformation_note.md`](../findings/m9_4_uv_deformation_note.md).
Q4b tree-level matching holds. Q4a remains `[O]`.
