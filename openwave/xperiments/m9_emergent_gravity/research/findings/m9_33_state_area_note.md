# M9.33: cut-correlator area responds at fixed \(H\) --- still an energy proxy

> Paper 42 killed hop-area as a second term. At fixed hops the
> only geometric object that can move is built from \(C\).

## Equations

Occupation transfer, \(H\) fixed. For each ball \(B\),

\[
A_{\mathrm{bond}}=\sum_{\langle ij\rangle\in\partial B}\lvert C_{ij}\rvert,
\qquad
A_{\mathrm{pur}}=\sum_k 4\lambda_k(1-\lambda_k).
\]

\(A_{\mathrm{pur}}\) is a function of the same spectrum as
\(S\). \(A_{\mathrm{bond}}\) is a cut observable of the
state, not of \(H\).

## Verdicts

\(N=12\), \(512\) balls. Instrument holds
(\(\rho(K_{\mathrm{vac}})=0.999\), \(\rho_\varepsilon=0.9995\)).

| Object | \(\rho\) with \(\delta S\) | \(\rho\) with \(P_{\mathrm{CHM}}\) |
| --- | --- | --- |
| \(A_{\mathrm{pur}}\) | \(0.99998\) (is \(S\)) | --- |
| \(A_{\mathrm{bond}}\) | \(-0.924\) | \(-0.925\) |
| \(P_{\mathrm{CHM}}\) | \(0.9998\) | \(1\) |

C_indepP **FAIL** (\(\lvert\rho\rvert=0.925>0.90\)).
C_eta **FAIL** (rel IQR \(0.87\)).
C_track PASS. C_indepS PASS (not \(S\)).

`BOND_IS_ENERGY_PROXY`. *computed.* The cut correlator
moves at fixed \(H\), unlike hop-area. It is not entropy.
It is the energy first law on the cut.

Auditor \(N=10\): independence QUALIFIED
(\(\rho=-0.818<0.90\), **CONFIRMED**). C_eta **REFUTED**
(rel IQR \(1.21\)).

Not a constant \(\eta\). Not \(8\pi G\). Not Einstein.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(A_{\mathrm{bond}}\), \(A_{\mathrm{pur}}\), gates | `scripts/m9_33_state_area.py` |
| Adversary | `scripts/m9_33_audit_state.py` |

Paper: [`../latex/43_State_Area.tex`](../latex/43_State_Area.tex).
