# M9.63: a ball that holds the pair has \(f_S=f_E=1\)

> Paper 72 tracked \(\delta S\) vs \(M_{AB}\) on \(R=3\),
> but those balls only touch the source sites.
> Paper 50: \(R=5\) is the enclosing radius.
> Three balls, exact Peschel, \(dps=80\).

## Equations

\[
S_{\mathrm{global}}=2h(\alpha_A)+2h(\alpha_B),\qquad
f_S=\delta S/S_{\mathrm{global}},\qquad
f_E=M_{AB}/M_{\mathrm{tot}}.
\]

ENC: midpoint \(R=5\). ENC_OFF: one site toward \(A\).
MISS: centred on \(A\), \(B\) site outside.

## Verdicts

\(N=12\). \(515\)-site blocks. \(S_4\) matches theory
to \(10^{-80}\).

| Ball | \(f_S\) | \(f_E\) | both sites |
| --- | --- | --- | --- |
| ENC \((6,6,6)\) \(R=5\) | \(0.9907\) | \(0.9881\) | yes |
| ENC_OFF \((5,6,6)\) \(R=5\) | \(0.886\) | \(0.881\) | yes |
| MISS \((3,6,6)\) \(R=5\) | \(0.382\) | \(0.386\) | no |

C_sg, C_fe, C_fs PRIMARY, C_miss **PASS**.
C_off **FAIL** (\(\lvert f-1\rvert=0.114>0.08\)):
both *sites* inside is not both *packets* inside.
The heavier tail leaks. \(f_S\) still tracks \(f_E\).

Auditor \(N=11\), perpendicular offset:
ENC \(0.994,0.994\); ENC_OFF \(0.990,0.991\);
MISS \(0.306,0.308\). All **CONFIRMED**.

`PAIR_ENC_IDENTITY`. *computed.* When the ball
holds the pair mass, finite \(\delta S\) is the
pair mixing entropy. A site-inside offset can
still leak. Not Einstein.

## Equation-to-code

| Object | Where |
| --- | --- |
| Three balls, `eigsy` \(dps=80\) | `scripts/m9_63_enclose_pair.py` |
| Adversary \(N=11\) | `scripts/m9_63_audit_enclose.py` |

Paper: [`../latex/73_Pair_Enclose.tex`](../latex/73_Pair_Enclose.tex).
