# M9.49: the complement is not a cosmological horizon

> Paper 20 (metric SdS): \(T\mathrm{d}S+\mathrm{d}M=0\), so
> \(\mathrm{d}S/\mathrm{d}M<0\) on the cosmological horizon.
> Ordinary balls have \(\delta S>0\) when enclosed energy
> rises. This note asks the complement.

## Equations

Fermi sea \(C_0\) is pure: \(S(B)=S(B^c)\). Occupation
transfer makes \(C_1\) mixed,
\(S_{\mathrm{global}}=2h(\alpha)\). Balls \(R=3,4,5\)
enclose the packet. Measure \(\delta S(B)\) and
\(\delta S(B^c)\).

C_comp: \(\delta S(B^c)<0\) at every \(R\).

## Verdicts

\(N=12\). Vacuum purity holds (\(\mathrm{rel}\sim 10^{-9}\)).
\(\delta S(B)\approx 0.196>0\) at every \(R\) (ordinary
first law). \(S_{\mathrm{global}}=0.196\): the new
entropy lives in \(B\).

| \(R\) | \(\delta S(B)\) | \(\delta S(B^c)\) |
| --- | --- | --- |
| \(3\) | \(+0.196\) | \(+4.7\times 10^{-5}\) |
| \(4\) | \(+0.196\) | \(-9.6\times 10^{-6}\) |
| \(5\) | \(+0.196\) | \(-4.6\times 10^{-6}\) |

C_comp **FAIL** (R=3 is plus). The complement moves
four orders of magnitude less than the ball. That is
not \(T\mathrm{d}S+\mathrm{d}M=0\).

Auditor \(N=10\): \(\delta S(B^c)\) also \(\sim 10^{-5}\)
and negative. C_comp **CONFIRMED** on the tiny minus.
The *magnitude* is the same story: not a horizon.

`COMPLEMENT_PLUS`. *computed.* The cosmological minus
sign of Paper 20 is metric SdS. This lattice complement
does not carry it at leading order.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(B\), \(B^c\), gates | `scripts/m9_49_complement.py` |
| Adversary | `scripts/m9_49_audit_comp.py` |

Paper: [`../latex/59_Complement.tex`](../latex/59_Complement.tex).
