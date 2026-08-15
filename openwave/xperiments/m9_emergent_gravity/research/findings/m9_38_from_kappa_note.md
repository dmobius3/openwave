# M9.38: \(\kappa\) reads the mass; inherited Newton does \(1/r^2\)

> Paper 47 used \(\kappa\) to weigh and locate. This note
> feeds that mass into the M9.2 DST Poisson solver.

## Equations

Even well-inside split (seed 38):

\[
\kappa=\mathrm{median}(\delta S/P_{\mathrm{flat}}).
\]

Odd half: \(M_{\mathrm{hat}}=\mathrm{median}(\delta S/\kappa)\),
compared to \(M_{\mathrm{enc}}=\mathrm{median}(P_{\mathrm{flat}})\).
Location: intersection of balls with \(\delta S>0.70\max\delta S\).

Inherited Poisson on a Dirichlet cube of half-width \(L=1\),
\(n=65\), \(G=1\), compact blob of mass \(M_{\mathrm{hat}}\)
at the centre:

\[
\nabla^2\Phi=4\pi G\rho,\qquad
\mathbf{a}=-\nabla_h\Phi.
\]

C_newt (same C1 as M9.2): attractive; \(\lvert |a|r^2/(GM_{\mathrm{hat}})-1\rvert<0.05\)
at \(r=0.30L,0.35L,0.40L\); log-log slope \(\lvert\alpha+2\rvert<0.08\).

C_newt is inherited Einstein sourced by an entanglement
mass. It is not a derivation of Poisson.

## Verdicts

\(N=12\), \(R=3\), source \((6,6,6)\). \(\kappa=0.968\).

| Gate | Result |
| --- | --- |
| C_vac | PASS \(\rho=0.974\) |
| C_mass | PASS, \(\lvert M_{\mathrm{hat}}/M_{\mathrm{enc}}-1\rvert=1.7\times 10^{-8}\) |
| C_loc | PASS, \(\hat x=(6,6,6)\) |
| C_newt PRIMARY | PASS. \(a_r<0\); residuals \(3.0\%\), \(2.1\%\), \(1.9\%\); slope \(-2.037\) |

Auditor \(N=10\), source \((4,5,5)\), own \(M_{\mathrm{hat}}\),
own \(n=65\) Poisson: C_mass **CONFIRMED** (\(0.18\%\)).
C_newt **CONFIRMED** (same residuals, same slope).

The C1 residuals are identical on the two masses because
Poisson is linear: \(a\propto M\), so
\(\lvert |a|r^2/(GM)-1\rvert\) is independent of
\(M_{\mathrm{hat}}\). C_newt would pass for any compact
positive mass at the cube centre. That is M9.2, reused.

`KAPPA_TO_NEWTON`. *computed.* Entanglement supplies
\(M_{\mathrm{hat}}\). Einstein still does \(1/r^2\).
Not a derived Poisson. Not \(8\pi G\) from \(\kappa\).
Not \(1/4G\). Not FGHMV. Not de Sitter. No `MODELS.md`.

## Equation-to-code

| Object | Where |
| --- | --- |
| Split, \(\kappa\), \(M_{\mathrm{hat}}\), DST C1 | `scripts/m9_38_from_kappa.py` |
| Adversary | `scripts/m9_38_audit_from.py` |

Paper: [`../latex/48_From_Kappa.tex`](../latex/48_From_Kappa.tex).
