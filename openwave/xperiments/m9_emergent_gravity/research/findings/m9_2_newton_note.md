# M9.2: inherited Newton limit --- force passes, isolated \(\Phi\) does not

> Locked 2026-08-15, run 2026-08-15. Palatini Einstein-Cartan
> with a spinless source *is* Einstein-Hilbert. This note
> checks the Newton limit of that inheritance. It does not
> discover gravity from entanglement.

> **Archive pointer (2026-08-22).** The campaign records this document cites (M9.3 to M9.73: notes, scripts, data, latex, APPLICATION) are not in the tree; they are parked in PR #441's branch history at `a5640709` and return as per-task PRs against registered roadmap rows. Their links below were reduced to plain text at the merge trim.

## Equations

Mostly-minus weak field (lock):

\[
\mathrm{d}s^2=-(1+2\Phi)\,\mathrm{d}t^2+(1-2\Phi)\,\delta_{ij}\,\mathrm{d}x^i\mathrm{d}x^j,
\qquad
\nabla^2\Phi=4\pi G\rho,
\qquad
\mathbf{a}=-\nabla\Phi.
\]

The last line is the slow-motion geodesic of that metric.
Cube \([-L,L]^3\), Dirichlet \(\Phi=0\) on \(\partial\).
Seven-point Laplacian, inverted by DST-I (exact for that
stencil). Discrete \(a_x=-(\Phi_{i+1}-\Phi_{i-1})/(2h)\).

Source: uniform ball, \(R=0.04L\), \(M=1\), \(G=1\),
rms \(0.029L<0.05L\). Grids \(n=65\) and \(n=97\).
Probes \(r\in\{0.30L,0.35L,0.40L\}\) on the \(x\)-axis.

## Verdicts

Finer grid \(n=97\). Isolated Newton
\(\lvert\mathbf{a}\rvert=GM/r^2\), \(\Phi=-GM/r\).

| Gate | Result | Numbers |
| --- | --- | --- |
| C1 PRIMARY | **PASS** | attractive; residuals \(1.29\%\), \(1.04\%\), \(1.07\%\); slope \(\alpha=-2.008\) |
| C2 | **FAIL** | \(\lvert\Phi r/(GM)+1\rvert=0.260\), \(0.305\), \(0.350\) |
| C3 vacuum | PASS | \(\Phi=\mathbf{a}=0\) |
| C4 2-d mutation | PASS | \(\alpha=-0.988\), C1(ii) fails |
| C5 \(G\to-G\) | PASS | \(\mathbf{a}\) flips; C1(ii) holds with \(\lvert G\rvert\) |

Coarse \(n=65\): C1 residuals \(3.0\%\), \(2.1\%\), \(1.9\%\),
\(\alpha=-2.037\). C2 still \(\sim 26\)--\(35\%\).

Auditor, Gaussian source, own gradient, \(n=61\):
C1 **CONFIRMED**, C2 **REFUTED**, C3--C5 **CONFIRMED**.

`NEWTON_INHERITED_PASS` on C1 only. C2 fails because a
Dirichlet cube is not isolated space: images pull \(\Phi\)
toward zero by a quarter at these probe radii. That is the
locked box, not a coding error. The force (a derivative)
still matches \(GM/r^2\) to the \(5\%\) cut.

This is not FGHMV, not entanglement gravity, not de Sitter,
not GEM, not Hehl-Datta. No `MODELS.md` column is added:
the official M9 ID is still unassigned, and C2 did not pass.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(\nabla^2\Phi=4\pi G\rho\), DST-I | `scripts/m9_2_newton_limit.py` (`poisson_3d`, `_dst_poisson`) |
| \(\mathbf{a}=-\nabla_h\Phi\), C1--C5 | same file, `score_grid`, `main` |
| Gaussian adversary | `scripts/m9_2_audit_newton.py` |

Paper: `../latex/40_Inherited_Newton.tex`.
Lock: [`../tasks/m9_2_task_details.md`](../tasks/m9_2_task_details.md).
