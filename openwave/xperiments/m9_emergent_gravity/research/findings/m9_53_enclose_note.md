# M9.53: \(\delta S/S_{\mathrm{global}}\) tracks \(P_{\mathrm{flat}}/M_{\mathrm{global}}\)

> Paper 60: on one enclosing ball, \(\delta S=S_{\mathrm{global}}=2h(\alpha)\).
> Paper 62: on 216 balls, \(\delta S\) tracks \(P_{\mathrm{flat}}\).
> This run asks whether they are the same fraction.

## Equations

\[
f_S=\frac{\delta S}{S_{\mathrm{global}}},\qquad
f_E=\frac{P_{\mathrm{flat}}}{M_{\mathrm{global}}},\qquad
S_{\mathrm{global}}=2h(\alpha),\qquad
M_{\mathrm{global}}=\sum_i\delta e_i.
\]

C_rho PRIMARY: \(\rho(f_S,f_E)>0.95\) at every \(\alpha\).
C_rms: \(\mathrm{RMS}(f_S-f_E)<0.10\).
C_enc: on source-inside balls, median \(\lvert f-1\rvert<0.05\).

## Verdicts

\(N=12\), \(R=3\), \(216\) balls, \(120\) source-inside.

| \(\alpha\) | \(\rho\) | RMS | med \(\lvert f_S-1\rvert\) | centered \(f_S\) | centered \(f_E\) |
| --- | --- | --- | --- | --- | --- |
| \(0.005\) | \(0.991\) | \(0.098\) | \(0.391\) | \(0.996\) | \(0.999\) |
| \(0.02\) | \(0.996\) | \(0.062\) | \(0.336\) | \(0.998\) | \(0.999\) |
| \(0.08\) | \(0.999\) | \(0.028\) | \(0.285\) | \(0.998\) | \(0.999\) |

C_rho **PASS**. C_rms **PASS**. C_enc **FAIL**.
The pre-registered well-inside set is leak: a
source site inside \(R=3\) does not enclose a
\(\sigma=1\) packet. Offset-binned, both
fractions fall together (Paper 49's \(28\%\)
gap again). On the centered ball both recover
\(1\) to \(<0.5\%\).

Auditor \(N=10\), own source, own \(\alpha\):
\(\rho=0.998\) at \(0.01\) and \(0.06\). C_rho
**CONFIRMED**. C_rms **CONFIRMED**. C_enc
**REFUTED**. Centered \(f_S=0.993,0.996\).

`ENCLOSURE_FRACTION_QUALIFIED`. *computed.*
Finite \(\delta S\) is mixing entropy times how
much of the packet sits in the ball. Gravity
is still \(\sum\delta e\) plus inherited Gauss.
Not Clausius, not \(8\pi G\), not de Sitter.

## Equation-to-code

| Object | Where |
| --- | --- |
| 216 balls, three \(\alpha\) | `scripts/m9_53_enclose.py` |
| Offset bins (same construction) | `data/m9_53_enclose_bins.json` |
| Adversary | `scripts/m9_53_audit_enclose.py` |

Paper: [`../latex/63_Enclosure.tex`](../latex/63_Enclosure.tex).
