# M9.45: first-law Gauss force, no Poisson solver

> Papers 48--54 solved Dirichlet Poisson to get \(\mathbf{a}\).
> Spherical Gauss is \(a(R)=-GM(R)/R^2\). This note uses
> \(M(R)=\delta S(R)/\kappa\) and never calls DST.

## Equations

\[
a(R)=-\frac{G}{R^2}\frac{\delta S(R)}{\kappa},\qquad G=1,
\]

with lattice \(R\). Same rule on two states.

STAR: open hop, compact packet. \(\kappa\) from the
smallest enclosing ball. SEA: periodic band-edge
transfer. \(\kappa\) from \(R=2\).

## Verdicts

\(N=12\). STAR encloses at \(R=3\). Fit \(R=3,4,5\):
slope \(-1.997\). C_star **PASS**.

SEA, \(R=2,3,4,5\): slope \(1.266\). Closer to \(+1\)
than to \(-2\). C_sea **PASS**.

C_split **PASS**. Auditor \(N=10\): star \(-1.889\),
sea \(1.237\), both **CONFIRMED**.

Gauss is inherited Newtonian integral calculus. The
input is the first-law mass. A star then falls as
\(1/R^2\). A uniform sea then grows as \(R\). Not a
derived Einstein equation. Not a de Sitter dual.

`GAUSS_TWO_LAWS`. *computed.*

## Equation-to-code

| Object | Where |
| --- | --- |
| Star, sea, Gauss \(a(R)\) | `scripts/m9_45_gauss_force.py` |
| Adversary | `scripts/m9_45_audit_gauss.py` |

Paper: [`../latex/55_Gauss_Force.tex`](../latex/55_Gauss_Force.tex).
