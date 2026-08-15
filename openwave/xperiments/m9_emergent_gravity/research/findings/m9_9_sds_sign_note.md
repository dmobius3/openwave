# M9.9 method note: cosmological minus sign derived from SdS

> [P]: Einstein+\(\Lambda\) \(\Rightarrow\) \(T\,\mathrm{d}S+\mathrm{d}M=0\)
> on the cosmological horizon. Not FGHMV. Not a prize claim.
> Paper: [`../latex/20_Derived_Sign_Einstein_Lambda.tex`](../latex/20_Derived_Sign_Einstein_Lambda.tex).

## VERDICT

On Schwarzschild--de Sitter, \(f=1-2M/r-r^2/\ell^2\), the cosmological
root satisfies \(\mathrm{d}r_c/\mathrm{d}M|_{M=0}=-1\) (*derived*,
auditor CONFIRMED identically). Then \(T\,\mathrm{d}S=-\mathrm{d}M\).
The minus sign is an identity of Einstein+\(\Lambda\), not an AdS
import. Finite-\(M\) residual \(\sim 10^{-6}\). Without \(\Lambda\)
there is no cosmological root.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(f(r)\), \(r_c(M)\) | `scripts/m9_9_sds_sign.py` `f`, `r_cosmo` |
| Implicit \(\mathrm{d}r/\mathrm{d}M\) | `scripts/m9_9_audit_sds.py` |

## Tags

| Statement | Tag |
| --- | --- |
| \(\mathrm{d}r_c/\mathrm{d}M\|_0=-1\), \(T\mathrm{d}S+\mathrm{d}M=0\) | *proved* |
| Local Clausius \(\Rightarrow\) Einstein+\(\Lambda\) | *cited* (Jacobson 1995) |
| Vacuum is dS | *proved* (Paper 19) |
| FGHMV for balls in dS | *not claimed* |
| Value of \(\Lambda\) | *not claimed* |
