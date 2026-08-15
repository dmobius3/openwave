# M9.39: Newton sees \(\sum\delta e\), not \(M_{\mathrm{hat}}\)

> Paper 48 put a point of mass \(M_{\mathrm{hat}}\) into Poisson.
> This note sources Poisson with the actual site-energy map \(\delta e\).

## Equations

Same occupation transfer as Papers 47--48. Site energy
\(\delta e_i=\sum_j H_{ij}(C_{ij}-C^{(0)}_{ij})\).
\(M_{\mathrm{global}}=\sum_i\delta e_i\).
\(M_{\mathrm{hat}}=\mathrm{median}(\delta S/\kappa)\) on the
odd well-inside split.

Embed \(\delta e\) as \(\rho\) on the M9.2 DST Poisson cube.
Far-field lock: 3 lattice sites map to \(0.05L\). Probes
\(0.30L,0.35L,0.40L\). Residual of a mass \(M\):

\[
\bigl\lvert |a|r^2/(GM)-1\bigr\rvert.
\]

C_hat / C_global: residual \(<0.05\), attractive,
\(\lvert\alpha+2\rvert<0.08\). C_which: smaller mean residual.

Near-field embedding (lattice fills the cube) is a
diagnostic: those probes sit inside the packet.

## Verdicts

\(N=12\), \(R=3\), source \((6,6,6)\). Point-mass control
still PASS (Paper 48 reproduced).

The packet is compact. Cumulative \(\delta e\): \(90\%\)
inside 1 site, \(99.7\%\) inside 2, essentially all inside 3.
Negative \(\delta e\) is \(10^{-6}\) dust. There is no distant
tail. \(M_{\mathrm{hat}}/M_{\mathrm{global}}-1=-28\%\) is
offset-ball bookkeeping, not leaked energy.

| Test | Result |
| --- | --- |
| Near field | not Coulomb (slope \(-0.065\)). Probes inside the packet. |
| Far C_global | PASS. Residuals \(2.9\%\), \(2.0\%\), \(1.9\%\). Slope \(-2.035\) |
| Far C_hat | FAIL. Residuals \(43\%\), \(41\%\), \(41\%\) |
| C_which PRIMARY | \(M_{\mathrm{global}}\) |

Auditor \(N=10\): C_global **CONFIRMED** (\(2.9\%\), \(2.0\%\),
\(1.9\%\)). C_hat **CONFIRMED FAIL** (\(29\%\), \(28\%\),
\(28\%\)). C_which **CONFIRMED**.

`TOTAL_ENERGY_IS_SOURCE`. *computed.* Real \(\delta e\)
sources inherited Newton with mass \(\sum\delta e\).
\(\kappa\) reads a ball-population functional, not the
far-field mass. Not a derived Poisson. Not \(8\pi G\).
Not \(1/4G\). Not FGHMV. Not de Sitter. No `MODELS.md`.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(\delta e\), near/far embed, gates | `scripts/m9_39_tail.py` |
| Adversary | `scripts/m9_39_audit_tail.py` |

Paper: [`../latex/49_Which_Mass.tex`](../latex/49_Which_Mass.tex).
