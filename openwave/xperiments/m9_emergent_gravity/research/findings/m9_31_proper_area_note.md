# M9.31: proper area does not save Clausius

> Paper 39 summed cut *lengths*. A surface should scale as
> \(\ell^2\). At linear order those two measures are
> proportional, and \(\eta\) still scatters.

## Equations

Same conformal hops as Paper 39. Three cut measures:

\[
A_{\mathrm{len}}=\sum_{\partial B}\frac{1}{\lvert t\rvert},
\qquad
A_{\mathrm{face}}=\sum_{\partial B}\frac{1}{t^2},
\qquad
A_{\mathrm{wf}}=\sum_{\partial B}\bigl(1-\varepsilon\Phi_{\mathrm{avg}}\bigr).
\]

With \(t=-(1+\varepsilon\Phi_{\mathrm{avg}})\),

\[
\delta A_{\mathrm{face}}=2\,\delta A_{\mathrm{len}}+O(\varepsilon^2)
=\mathrm{const}\times\delta A_{\mathrm{wf}}+O(\varepsilon^2).
\]

Relative IQR of \(\delta S/\delta A\) is scale-invariant, so
the three \(\eta\) tests are the same test. That is derived,
then confirmed on the lattice.

## Verdicts

\(N=12\), \(512\) balls. Instrument holds
(\(\rho(K_{\mathrm{vac}})=1-8\times 10^{-9}\)).

| Measure | \(\rho(\delta S,\delta A)\) | rel IQR \(\eta\) | C_eta |
| --- | --- | --- | --- |
| length (Paper 39) | \(0.918\) | \(3.173\) | FAIL |
| face PRIMARY | \(0.918\) | \(3.174\) | **FAIL** |
| weak-field | \(0.919\) | \(3.172\) | FAIL |

Energy still beats every geometric predictor
(\(R_{\mathrm{CHM}}=0.317<R_{\mathrm{face}}=0.397\)).

Auditor \(N=10\): C_vac **CONFIRMED**, C_area **CONFIRMED**,
C_eta **REFUTED** (rel IQR \(4.96\)).

`FACE_CORRELATES_NOT_CLAUSIUS`. The excuse that Paper 39
used the wrong dimension is **refuted**. Not \(1/4G\).
Not Einstein. Not a `MODELS.md` cell.

## Equation-to-code

| Object | Where |
| --- | --- |
| Three areas, gates | `scripts/m9_31_proper_area.py` |
| Adversary, face only | `scripts/m9_31_audit_area.py` |

Paper: [`../latex/41_Proper_Area.tex`](../latex/41_Proper_Area.tex).
