# M9.56: kinetic stress of the transfer is dust. Not \(\Lambda\)

> Papers 54--65: mass is \(\sum\delta e\); Gauss is inward.
> de Sitter needs \(p/\rho=-1\). This run measures the
> hop-dispersion virial, not \(p=-E/V\) (forbidden:
> that identity fakes \(\Lambda\) at \(\mu=0\)).

## Equations

\[
\varepsilon(\mathbf{k})=-2\sum_\mu\cos k_\mu,\qquad
w(\mathbf{k})=\frac23\sum_\mu k_\mu\sin k_\mu,
\]

\(\mathbf{k}\) in the first Brillouin zone.

\[
r=\frac{\sum_m\delta n_m\,w_m}{\sum_m\delta n_m\,\varepsilon_m}=\frac{\delta P}{\delta E}.
\]

STAR: open product modes, \(k_\mu=\pi n_\mu/(N+1)\).
SEA: \(L=0\), \(R=(\pi,\pi,\pi)\). C_lambda PRIMARY:
\(\lvert r_{\mathrm{sea}}+1\rvert<0.25\).

## Verdicts

\(N=12\). \(\alpha\in\{0.01,0.02,0.04\}\). \(r\) does
not run.

| State | \(r=\delta P/\delta E\) |
| --- | --- |
| sea transfer | \(0\) (to \(10^{-16}\)) |
| star packet | \(0.163\) |
| vacuum (not a source) | \(P/E=-0.881\) |

C_e **PASS** (mode sum matches site energy,
\(4\times 10^{-15}\)). C_lambda **FAIL**.
C_dust **PASS**. C_rad **FAIL**. C_hold **PASS**.
Auditor \(N=10\): sea \(0\), star \(0.170\).
C_lambda **REFUTED**. C_dust **CONFIRMED**.

Band edges have vanishing virial (\(\sin 0=\sin\pi=0\)).
They carry energy and no kinetic pressure. That is
dust. Vacuum \(P/E\approx-0.88\) is the filled sea;
Paper 58 subtracts it. \(p=-E/V\) is not used.

`STRESS_DUST`. *computed.* This transfer cannot
source de Sitter. Gravity remains \(\sum\delta e\)
plus inherited Gauss.

## Equation-to-code

| Object | Where |
| --- | --- |
| Virial, star, sea, vacuum | `scripts/m9_56_stress.py` |
| Adversary | `scripts/m9_56_audit_stress.py` |

Paper: [`../latex/66_Stress_Dust.tex`](../latex/66_Stress_Dust.tex).
