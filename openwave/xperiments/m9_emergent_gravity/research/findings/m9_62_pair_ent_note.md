# M9.62: pair \(S_{\mathrm{global}}=2h(\alpha_A)+2h(\alpha_B)\). Enclosure tracks. \(\nabla\cdot\mathbf g\) is not Poisson.

> Rank-4 occupation of the orthonormal two-source
> state. Exact open-hop basis. Peschel by
> `mpmath eigsy` at \(dps=80\). No hop LAPACK.

## Equations

If \(\{L_A,L_B\}\) and \(\{R_A,R_B\}\) are orthonormal
and occupied \(\perp\) empty,

\[
\mathrm{spec}_{\mathrm{move}}(C)=\{\alpha_A,1-\alpha_A,\alpha_B,1-\alpha_B\},
\qquad
S_{\mathrm{global}}=2h(\alpha_A)+2h(\alpha_B).
\]

\[
f_S=\delta S/S_{\mathrm{global}},\qquad
f_E=M_{AB}/M_{\mathrm{tot}}.
\]

\(\mathbf g=-(M/A)\hat n\) as in Papers 69--70.
Lattice \(\nabla\cdot\mathbf g\) on the \(R=2\) field.

## Verdicts

Gram error \(8\times 10^{-81}\). \(C\) off-diagonal
\(2\times 10^{-84}\). \(S_4=S_{\mathrm{theory}}\)
to \(10^{-80}\). *derived* then *computed*.

\(N=12\), \(64\) two-sided \(R=3\) balls.
Peschel \(123\times 123\), \(dps=80\).

| Gate | Solver | Auditor \(N=11\) |
| --- | --- | --- |
| C_rho \(\rho(f_S,f_E)\) | **PASS** \(0.957\) | **CONFIRMED** \(0.980\) |
| C_rms | **PASS** \(0.043\) | **CONFIRMED** \(0.063\) |
| C_both \(\lvert f-1\rvert\) both-inside | **FAIL** \(0.73,0.62\) (\(n=1\)) | **REFUTED** \(0.75,0.61\) |
| C_src \(\lvert\nabla\cdot g\rvert\) source vs mid | **PASS** ratio \(3.83\) | **REFUTED** mid larger |

`PAIR_ENT_TRACKS`. *computed.* Finite pair
\(\delta S\) tracks enclosed pair mass. Absolute
\(f=1\) on a ball that only *touches* both sites
is leak (Paper 63 again). Discrete
\(\nabla\cdot\mathbf g\) is not a local Newton
source law: signs mix, and the auditor's
midpoint is larger than the near-packet values.
\(\mathbf g\) remains a ball probe, not a
derived Poisson field. Not Einstein.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(4\times 4\) \(C\), Peschel, \(\nabla\cdot g\) | `scripts/m9_62_pair_ent.py` |
| Adversary \(N=11\) | `scripts/m9_62_audit_pair.py` |

Paper: [`../latex/72_Pair_Entropy.tex`](../latex/72_Pair_Entropy.tex).
