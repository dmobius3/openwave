# M9.46: the Gauss slope interpolates. The force is dust, not \(\Lambda\)

> Paper 54 called the uniform sea a Newtonian \(\Lambda\)
> signature. Interior Newton of positive \(\rho\) is
> \(a=-(4\pi G\rho/3)r\) (inward). de Sitter / \(\Lambda>0\)
> is \(a=+(\Lambda/3)r\) (outward). Every allowed occupation
> swap on the Fermi sea *raises* \(\langle H\rangle\). There
> is no negative-\(\delta e\) sea in this instrument.

## Equations

First-law Gauss of Paper 55, \(G=1\):

\[
a(R)=-\frac{1}{R^2}\frac{\delta S(R)}{\kappa}.
\]

Open hop, \(\sigma\in\{1,2,4,8\}\). PBC band-edge as the
flat endpoint. Slope of \(\lvert a\rvert\) vs \(R\) on
\(R=2,3,4,5\).

## Verdicts

\(N=12\). All \(a(R)<0\).

| \(\sigma\) | slope |
| --- | --- |
| \(1\) | \(-1.868\) |
| \(2\) | \(-0.903\) |
| \(4\) | \(+0.404\) |
| \(8\) | \(+1.001\) |
| PBC sea | \(+1.266\) |

C_star, C_mono, C_wide, C_sea, C_in **PASS**.
Auditor \(N=10\): \(-1.837\to +0.601\), inward
**CONFIRMED**, monotone **CONFIRMED**.

`DUST_NOT_LAMBDA`. *computed.* The first-law Gauss
slope interpolates from a star to a sea. The force is
inward at every width. Paper 54's "Newtonian \(\Lambda\)"
is withdrawn. This is the interior of positive energy,
not de Sitter.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(\sigma\) scan, sea, sign | `scripts/m9_46_sigma.py` |
| Adversary | `scripts/m9_46_audit_sigma.py` |

Paper: [`../latex/56_Dust_Not_Lambda.tex`](../latex/56_Dust_Not_Lambda.tex).
