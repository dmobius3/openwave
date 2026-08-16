# M9.60: pair direction from the exact open-hop basis

> Paper 69 used float64 LAPACK for \(H\). The open
> cube hop is separable. This run drops the eigensolve.

## Equations

\[
k_n=\frac{\pi n}{N+1},\qquad
\varepsilon_n=-2\cos k_n,\qquad
\psi_j^{(n)}=\sqrt{\frac{2}{N+1}}\sin\bigl(k_n(j+1)\bigr).
\]

Occupied iff \(\varepsilon_x+\varepsilon_y+\varepsilon_z<0\).
Same pair field as Paper 69:

\[
\hat n=\frac{\nabla_c M_{AB}}{\lvert\nabla_c M_{AB}\rvert},\qquad
\mathbf g=-\frac{M_{AB}}{A}\,\hat n.
\]

`mpmath` \(dps=80\). No LAPACK.

## Verdicts

\(N=12\), \(864\) occupied modes. \(M_B/M_A=2\).
62 leaking centres.

| Gate | Exact basis | Paper 69 (LAPACK) |
| --- | --- | --- |
| median \(\angle(\hat n,\mathbf g_N)\) | \(10.746^\circ\) | \(10.746^\circ\) |
| median \(\angle(\hat n,\mathrm{CM})\) | \(62.173^\circ\) | \(62.173^\circ\) |
| \(\lvert\nabla M\rvert\) null / CM | \(0.012/0.080\) | \(0.012/0.080\) |
| equal-pair midpoint ratio | \(2.50\times 10^{-4}\) | \(2.50\times 10^{-4}\) |
| additivity | \(6.15\times 10^{-8}\) | \(6.15\times 10^{-8}\) |

C_dir, C_notcm, C_null, C_mid, C_cert **PASS**.
C_add **FAIL** at the \(10^{-12}\) pre-register:
the \(6\times 10^{-8}\) survives the exact basis.
It is orthonormalization of overlapping packets,
not rounding.

Auditor \(N=11\), own pair, \(3\) axis sites:
C_dir **CONFIRMED** (\(10.92^\circ\)). C_notcm
**CONFIRMED** (\(10.92^\circ<39.7^\circ\)).
C_null **CONFIRMED**.

`EXACT_PAIR_DIRECTION`. *derived* (basis) /
*computed* (angles). Paper 69 is physics.
The pair field is \(\mathbf g=-(M/A)\hat n\).
Not Einstein.

## Equation-to-code

| Object | Where |
| --- | --- |
| Closed-form 1d modes, pair \(M_{AB}\) | `scripts/m9_60_exact.py` |
| Adversary \(N=11\) | `scripts/m9_60_audit_exact.py` |

Paper: [`../latex/70_Exact_Basis.tex`](../latex/70_Exact_Basis.tex).
