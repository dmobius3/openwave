# M9.57: a thermal scan is not de Sitter

> Paper 66: \(T=0\) transfer is dust. This run heats
> the periodic Fermi gas, \(\mu=0\), same virial.
> Source is the excess over \(T=0\) (Paper 58).

## Equations

\[
n_k(T)=\frac{1}{e^{\varepsilon_k/T}+1},\qquad
r(T)=\frac{P(T)-P(0)}{E(T)-E(0)}.
\]

C_lambda PRIMARY: \(\lvert r+1\rvert<0.25\) at any
kept \(T\). Forbidden: \(p=-E/V\).

## Verdicts

\(N=12\), \(T\in\{0.1,0.25,0.5,1,2,4\}\).
Raw \(P(0)/E(0)=-0.882\) is not a source.

| \(T\) | \(r=\delta P/\delta E\) |
| --- | --- |
| \(0.1\) | \(-57.8\) |
| \(0.25\) | \(-5.57\) |
| \(0.5\) | \(-0.959\) |
| \(1\) | \(-0.162\) |
| \(2\) | \(+0.018\) |
| \(4\) | \(+0.069\) |

C_lambda **trips** at \(T=0.5\). C_pos **FAIL**.
C_nolam **FAIL**. The ratio *runs*. It is not a
constant \(-1\). Low \(T\) is dominated by
\(\varepsilon\approx 0\) modes with \(w=\pi\)
(\(n=1\to\frac12\), \(\delta E\approx 0\),
\(\delta P\) large).

Auditor \(N=10\), own \(T\in\{0.2,0.8,3\}\):
\(r=+0.719,+0.233,+0.169\). All positive.
C_lambda **REFUTED**. C_pos **CONFIRMED**.

`THERMAL_NOT_LAMBDA`. *computed.* A one-lattice
crossing is not \(\Lambda\). Gravity remains
\(\sum\delta e\) plus inherited Gauss.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(N=12\) \(T\) scan | `scripts/m9_57_thermal.py` |
| Adversary | `scripts/m9_57_audit_thermal.py` |

Paper: [`../latex/67_Thermal.tex`](../latex/67_Thermal.tex).
