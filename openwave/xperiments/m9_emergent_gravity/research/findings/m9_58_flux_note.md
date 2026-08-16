# M9.58: flux law \(g=-M/A\). The hop has no \(1/r\)

> Entanglement supplies \(M(B)=\sum_{i\in B}\delta e\).
> Stokes plus isotropy is the only step that produces
> a radial field: \(g_R=-M/A\). Paper 61 used
> \(-M/R^2\). This run uses \(A\). Two-packet
> \(E_{\mathrm{int}}(d)\) is computed, not assumed.

## Equations

If a local \(\mathbf{g}\) exists and
\(\oint_{\partial B}\mathbf{g}\cdot d\mathbf{A}=-M(B)\),
isotropy on a source-centered ball forces

\[
g_R=-\frac{M}{A},\qquad
A=\text{outgoing NN count}.
\]

That is not a theorem of the hop. The hop energy
of two orthonormal packets:

\[
E_{\mathrm{int}}(d)=E_{AB}(d)-E_A-E_B.
\]

C_flat PRIMARY: \(\max_d\lvert E_{\mathrm{int}}\rvert/(\lvert E_A\rvert+\lvert E_B\rvert)<0.02\).

## Verdicts

\(N=12\), \(\alpha=0.02\). \(A/R^2\approx 19.3\)--\(19.5\).

| State | slope \(-M/A\) | slope \(-M/R^2\) | fit |
| --- | --- | --- | --- |
| star \(R=3,4,5\) | \(-1.9998\) | \(-1.9979\) | enclose |
| sea \(R=2,3,4,5\) | \(+0.992\) | \(+0.967\) | all \(R\) |

C_star, C_sea, C_class **PASS**.

| \(d\) | \(E_{\mathrm{int}}\) | rel |
| --- | --- | --- |
| \(2\) | \(-1.05\times 10^{-2}\) | \(0.0282\) |
| \(3\) | \(-1.89\times 10^{-3}\) | \(0.0051\) |
| \(4\) | \(-1.16\times 10^{-4}\) | \(0.0003\) |
| \(5\) | \(-2.38\times 10^{-6}\) | \(6\times 10^{-6}\) |
| \(6\) | \(-1.72\times 10^{-8}\) | \(5\times 10^{-8}\) |

C_flat **FAIL** at \(d=2\) (overlap, \(\sigma=1\)).
The tail is not \(1/r\) (successive ratios \(5.6,16,48\)).
Auditor \(N=10\): slopes \(-1.867\), \(+1.061\)
**CONFIRMED**. C_flat **REFUTED** (\(d=2\) rel \(0.026\)).

`FLUX_GAUSS_OVERLAP`. *computed.* The force law
used from here is \(g=-M/A\). It is Stokes plus
the measured mass. It is not in \(H\). Not Einstein.

## Equation-to-code

| Object | Where |
| --- | --- |
| Flux \(a=-M/A\), pair \(E_{\mathrm{int}}\) | `scripts/m9_58_flux.py` |
| Adversary | `scripts/m9_58_audit_flux.py` |

Paper: [`../latex/68_Flux.tex`](../latex/68_Flux.tex).
