# M9.42: the interior null is \(M/r^2\), not the centre of mass

> Paper 50 cancelled at the midpoint because the masses
> were equal. That is a symmetry. Between two attractors
> the force balances at \(M_A/r_A^2=M_B/r_B^2\), closer
> to the *lighter* mass.

## Equations

Unequal occupation transfers \(\alpha_A=0.02\),
\(\alpha_B=0.05\). \(M_A,M_B=\sum\delta e\).

\[
x_{\mathrm{cm}}=\frac{M_A x_A+M_B x_B}{M_A+M_B},\qquad
\frac{M_A}{(x-x_A)^2}=\frac{M_B}{(x_B-x)^2}.
\]

\(x_{\mathrm{null}}\) is the interpolated zero of \(a_x\)
on the segment. Inherited DST Poisson.

A first draft asked for the centre of mass. That was
wrong. The gates above are the corrected pre-register.

## Verdicts

\(N=12\), \(A=(2,6,6)\), \(B=(10,6,6)\), \(M_B/M_A=2.46\).

| Gate | Result |
| --- | --- |
| C_uneq | PASS, ratio \(2.46\) |
| C_force PRIMARY | PASS, \(\lvert x_{\mathrm{null}}-x_{\mathrm{force}}\rvert=0.0025L\) |
| C_notcm | PASS, \(x_{\mathrm{null}}\) is \(0.16L\) from the CM |
| C_side | PASS, both on the light-mass side |

Auditor \(N=10\): C_side **CONFIRMED**. C_notcm
**CONFIRMED**. C_force **REFUTED** (\(0.15L\)): the
auditor packets are fatter relative to the gap, so the
zero sits in the light near field. The side is Newton;
the millimetre is not.

`INVSQ_NULL`. *computed.* Unequal real densities cancel
closer to the lighter mass, at the inverse-square point,
not at the centre of mass. Not a derived Poisson.

## Equation-to-code

| Object | Where |
| --- | --- |
| Unequal pair, null, gates | `scripts/m9_42_bary.py` |
| Adversary | `scripts/m9_42_audit_bary.py` |

Paper: [`../latex/52_Invsq_Null.tex`](../latex/52_Invsq_Null.tex).
