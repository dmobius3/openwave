# M9.30: area correlates, Clausius does not --- not gravity

> Papers 37--38 bookkeep \(\delta S\) against energy at fixed
> hops. That cannot be gravity: hop-length area does not
> move. This note changes the hops.

## Equations

Conformal metric bump on the hop graph

\[
t_{ij}(\varepsilon)=-\Bigl(1+\varepsilon\frac{\Phi_i+\Phi_j}{2}\Bigr),
\qquad
\Phi_i=\exp\bigl(-r_i^2/(2\sigma^2)\bigr).
\]

Lattice area of a ball (independent of \(S\))

\[
A=\sum_{\langle ij\rangle\in\partial B}\frac{1}{\lvert t_{ij}\rvert}.
\]

Ground state of \(H(\varepsilon)\). Instrument
\(\delta S\stackrel{?}{=}\mathrm{Tr}(K_{\mathrm{vac}}\Delta C)\).
Clausius would be a *single* \(\eta\) with \(\delta S=\eta\,\delta A\).

Matter-only identity: hops fixed \(\Rightarrow\delta A\equiv 0\).

## Verdicts

\(N=12\), \(R=2\), \(512\) balls, \(\sigma=2\), \(\varepsilon=0.02\).
Occupancy stable \(864/864/864\).

| Gate | Lock | Result |
| --- | --- | --- |
| C_vac | \(\lvert\rho(\delta S,\mathrm{Tr}K_{\mathrm{vac}}\Delta C)\rvert>0.95\) | PASS \(1-8\times 10^{-9}\) |
| C0 | \(\max\|\delta S\|>10^{-6}\), \(\max\|\delta A\|>10^{-4}\) | PASS \(0.008\), \(0.91\) |
| C1 | Pearson \(\delta S(\varepsilon),\delta S(2\varepsilon)>0.95\) | PASS \(0.99999\) |
| C_area | \(\lvert\rho(\delta S,\delta A)\rvert>0.80\) | PASS \(0.918\) |
| C_eta PRIMARY | IQR\((\delta S/\delta A)/\|\mathrm{med}\|<0.35\) | **FAIL** \(3.17\) |
| C_pred | \(R_{\mathrm{area}}<R_{\mathrm{CHM}}\) | FAIL \(0.397>0.317\) |

Energy still beats area. \(\eta_{\mathrm{median}}=-1.7\times 10^{-3}\)
is not \(1/4G\).

Auditor \(N=10\), off-center \(\Phi\): C_vac **CONFIRMED**,
C_area **CONFIRMED**, C_eta **REFUTED** (rel IQR \(4.96\)).

`AREA_CORRELATES_NOT_CLAUSIUS`. *computed.* A metric
variation makes \(\delta S\) and \(\delta A\) move together.
It does not make a universal Clausius factor. It is not
Einstein. It is not a `MODELS.md` cell. Inherited Newton
(M9.2) is a different task and is still un-run.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(\Phi\), \(t(\varepsilon)\), \(A\), gates | `scripts/m9_30_area.py` |
| Adversary | `scripts/m9_30_audit_area.py` |

Paper: [`../latex/39_Area_Not_Clausius.tex`](../latex/39_Area_Not_Clausius.tex).
