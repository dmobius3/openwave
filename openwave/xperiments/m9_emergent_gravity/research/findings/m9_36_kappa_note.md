# M9.36: a universal Gauss first-law constant

> Paper 45: well inside, \(\delta S\) tracks enclosed energy.
> This note asks whether the ratio is one number.

## Equations

Well-inside balls (every source in \(B\)):

\[
\kappa=\frac{\delta S}{P_{\mathrm{flat}}},\qquad
P_{\mathrm{flat}}=\sum_{i\in B}\delta e_i.
\]

Config 1: one packet. Config 2: two packets. Predict
\(\delta S^{(2)}\stackrel{?}{=}\kappa_1 P_{\mathrm{flat}}^{(2)}\).

## Verdicts

\(N=12\), \(R=3\). Instrument holds.

| Config | \(n\) | med \(\kappa\) | rel IQR | \(\rho(\delta S,P)\) |
| --- | --- | --- | --- | --- |
| one mass | \(120\) | \(0.968\) | \(0.084\) | \(0.999\) |
| two masses | \(69\) | \(0.984\) | \(0.099\) | \(0.999\) |

\(\lvert\kappa_1-\kappa_2\rvert/\mathrm{mean}=0.017\). Cross
prediction Pearson \(0.999\).

C_univ **PASS**. C_pred **PASS**. Auditor \(N=10\):
\(\kappa=0.972\) vs \(0.992\) (rel \(0.021\)), **CONFIRMED**.

`KAPPA_UNIVERSAL`. *computed.* One lattice constant, to
about two percent, takes enclosed energy to \(\delta S\).
It is not \(1/4G\). It is not Einstein. It is the Gauss
first law of this free fermion, with a reusable \(\kappa\).

## Equation-to-code

| Object | Where |
| --- | --- |
| Configs, \(\kappa\), gates | `scripts/m9_36_kappa.py` |
| Adversary | `scripts/m9_36_audit_kappa.py` |

Paper: [`../latex/46_Kappa.tex`](../latex/46_Kappa.tex).
