# M9.41: stars plateau; an extended source does not

> Compact packets (Papers 47--50) are stars: \(\delta S\)
> saturates once the ball holds the packet. A wide packet
> is an extended source. The first law should keep growing.
> Nested volume vs area is collinear (Paper 42). Linear
> \(a\propto r\) (Newtonian \(\Lambda\)) is a Poisson
> question; this occupation transfer is not a uniform fluid.

## Equations

Same occupation transfer. Wide \(\sigma=8\), compact
\(\sigma=1\). Source-centered balls \(R=2,3,4,5\).
\(P_{\mathrm{flat}}=\sum_{B}\delta e\). Growth
\(\delta S(5)/\delta S(2)\).

## Verdicts

\(N=12\). First law tracks enclosed energy on the wide
source (\(\rho=0.999\)). Wide growth \(16.1\). Compact
growth \(1.13\) (plateau). C_fl / C_grow **PASS**.
Compact control **PASS**.

Volume vs area: \(\rho_V=0.996\), \(\rho_A=1.000\).
Nested balls are collinear. Discriminator withdrawn.

Wide \(\delta e\) is not uniform (\(\mathrm{std}/\mathrm{mean}=0.35\)).
Fill-cube \(a_r\) changes sign. C_lin **FAIL**. Auditor
C_lin **REFUTED**. This is not Newtonian \(\Lambda\).

`EXTENDED_FIRST_LAW`. *computed.* The first law can tell
a star from an extended source without Poisson. It did
not produce de Sitter.

## Equation-to-code

| Object | Where |
| --- | --- |
| Wide/compact scan, Poisson diagnostic | `scripts/m9_41_uniform.py` |
| Adversary | `scripts/m9_41_audit_uniform.py` |

Paper: [`../latex/51_Extended.tex`](../latex/51_Extended.tex).
