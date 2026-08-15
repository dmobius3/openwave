# M9.40: two real packets, one enclosing ball, midpoint cancels

> Fewer assumptions than Papers 48--49: no \(M_{\mathrm{hat}}\),
> no point mass, no 60-ball median. Mass from one
> source-centered first-law ball. \(\rho\) is the actual
> pair \(\delta e\). Still inherited Poisson.

## Equations

Orthonormal two-source occupation transfer (M9.34).
\(\kappa\) from packet \(A\), source-centered \(R=3\):

\[
\kappa=\delta S_A/P_{\mathrm{flat},A},\qquad
M_{\mathrm{FL}}=\delta S_{AB}/\kappa.
\]

Pair ball at the midpoint. \(R=5\) on \(N=12\):
\(R=4\) encloses only \(82\%\) of \(\sum\delta e\); Gauss
needs the packets, not just the two site-centres.

Inherited DST Poisson, \(G=1\). Midpoint at the origin.
Separation \(0.20L\). Exterior compared to two-point
Coulomb of \((M_A,x_A)\) and \((M_B,x_B)\).

\(a_{AB}=a_A+a_B\) is linear Poisson. Sanity only.

## Verdicts

\(N=12\), \(A=(3,6,6)\), \(B=(9,6,6)\). \(C\) in \([0,1]\).
Energy additive to \(5\times 10^{-8}\).

| Gate | Result |
| --- | --- |
| C_read PRIMARY | PASS, \(\lvert M_{\mathrm{FL}}/M_{AB}-1\rvert=0.84\%\) (enclose \(98.8\%\)) |
| C_mid PRIMARY | PASS, \(\lvert a(0)\rvert/a_{\mathrm{char}}=3.7\times 10^{-5}\) |
| C_ext | PASS, Poisson vs Coulomb \(0.09\%\), \(0.57\%\), \(0.81\%\) |
| C_super sanity | \(8.6\times 10^{-8}\) (not a discovery) |

Auditor \(N=10\): C_read **CONFIRMED** (\(0.38\%\)).
C_mid **CONFIRMED** (\(4.4\times 10^{-4}\)).
C_ext **CONFIRMED** (\(0.8\)--\(1.4\%\)).

`PAIR_NEWTON`. *computed.* One enclosing ball and the
one-mass \(\kappa\) read the pair. Two real densities
cancel in the middle and match two-point Newton outside.
Not a derived Poisson. Not \(8\pi G\). Not de Sitter.
No `MODELS.md`.

## Equation-to-code

| Object | Where |
| --- | --- |
| Pair, \(\kappa\), Poisson, gates | `scripts/m9_40_pair.py` |
| Adversary | `scripts/m9_40_audit_pair.py` |

Paper: [`../latex/50_Pair_Newton.tex`](../latex/50_Pair_Newton.tex).
