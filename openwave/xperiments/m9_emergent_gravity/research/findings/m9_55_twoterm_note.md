# M9.55: grow \(=2.5\) is not a universal area term

> Paper 64: sea-transfer \(\delta S\) is extensive.
> Slab grow is \(2.5\), not \(3\). A two-term
> \(\delta S=aV+bA\) would be the remaining
> Clausius/horizon handle. In-sample fit is not
> enough: the term must predict a held-out shape.

## Equations

\[
\delta S \stackrel{?}{=} aV+bA
\qquad\text{(no intercept).}
\]

Volume-only: \(\delta S=cV\).
C_loo PRIMARY: leave one family out; \(\rho>0.90\)
and rel RMS \(<0.25\) on the held-out family.

## Verdicts

\(N=12\). Fourteen regions: balls, cubes, slabs, rods.

| \(\alpha\) | \(a\) | \(b\) | rel RMS two | rel RMS vol |
| --- | --- | --- | --- | --- |
| \(0.01\) | \(8.04\times 10^{-5}\) | \(-2.0\times 10^{-7}\) | \(0.081\) | \(0.081\) |
| \(0.02\) | \(1.42\times 10^{-4}\) | \(3.8\times 10^{-6}\) | \(0.073\) | \(0.074\) |
| \(0.04\) | \(2.48\times 10^{-4}\) | \(1.5\times 10^{-5}\) | \(0.066\) | \(0.069\) |

C_fit **PASS**. C_gain **FAIL** (area does not
earn its keep). C_loo **FAIL**: held-out rods
rel RMS \(\approx 5.2\); \(b\) flips sign when
slabs are held out. Within-family \(\rho\) stays
high because size is monotone. That is not
accuracy.

Auditor \(N=10\): C_fit **CONFIRMED**. C_gain
**REFUTED**. C_loo **REFUTED** (rods \(\approx 5.5\)).

`TWO_TERM_FAIL`. *computed.* Paper 64's leftover
is shape junk, not a horizon piece. Volume
already explains the sea transfer. Gravity is
still \(\sum\delta e\) plus inherited Gauss.
Not Clausius. Not \(1/4G\).

## Equation-to-code

| Object | Where |
| --- | --- |
| Four families, LOO | `scripts/m9_55_twoterm.py` |
| Adversary | `scripts/m9_55_audit_twoterm.py` |

Paper: [`../latex/65_Two_Term_Fail.tex`](../latex/65_Two_Term_Fail.tex).
