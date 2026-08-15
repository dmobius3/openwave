# M9.32: hop-area is an energy proxy, not a second term

> Jacobson / FGHMV is \(\delta S=\eta\,\delta A+\beta\,\delta E\).
> Papers 39--41 saw \(\delta S\sim\delta A\) with no constant
> \(\eta\). This note asks whether \(\delta A\) is independent
> of \(P_{\mathrm{CHM}}\).

## Equations

Two variations, same \(512\) balls. Face area
\(A=\sum_{\partial B}1/t^2\). Energy predictor
\(P=\sum_B(R^2-r^2)\,\delta e\).

Matter: occupation transfer, hops fixed, \(\delta A\equiv 0\).
Geometry: hop conformal bump, new ground state.

Partial correlation (only scored if \(\lvert\rho(A,P)\rvert\le 0.90\))

\[
\rho(\delta S,\delta A\mid P)
=\frac{\rho_{SA}-\rho_{SP}\rho_{AP}}{\sqrt{(1-\rho_{SP}^2)(1-\rho_{AP}^2)}}.
\]

If \(\lvert\rho_{AP}\rvert>0.90\), \(A\) and \(P\) are the same
bump. The partial is then a degeneracy artifact and is not
cashed as a second term.

## Verdicts

Instrument: matter \(\rho(K_{\mathrm{vac}})=0.999\), geometry
\(1-8\times 10^{-9}\). Matter: \(\max\lvert\delta A\rvert=0\),
\(\rho(\delta S,P)=0.9998\).

| Geometry object | Value |
| --- | --- |
| \(\rho(\delta A,P)\) | \(0.996\) |
| \(\rho(\delta S,\delta A)\) | \(0.918\) |
| \(\rho(\delta S,P)\) | \(0.948\) |
| raw partial \(\rho(\delta S,\delta A\mid P)\) | \(-0.987\) (not scored) |
| \(R(P)\) / \(R(A)\) / \(R(A,P)\) | \(0.317\) / \(0.397\) / \(0.050\) |

`AREA_IS_ENERGY_PROXY`. *computed.* The two-column residual
drop is what a duplicate predictor does. Matter-only already
gives a first law with \(\delta A=0\).

Auditor \(N=10\): \(\rho(A,P)=0.988\), independence
**REFUTED**.

Papers 39--41's area correlation is the energy first law
in geometric clothing. Not a two-term Clausius law. Not
\(8\pi G\). Not Einstein.

## Equation-to-code

| Object | Where |
| --- | --- |
| Both variations, partial, gates | `scripts/m9_32_two_term.py` |
| Adversary | `scripts/m9_32_audit_term.py` |

Paper: [`../latex/42_Area_Is_Energy.tex`](../latex/42_Area_Is_Energy.tex).
