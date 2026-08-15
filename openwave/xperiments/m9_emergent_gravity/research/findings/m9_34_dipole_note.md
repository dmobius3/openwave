# M9.34: two masses --- CHM on all balls, Gauss not killed inside

> A charge dipole is impossible on the Fermi-sea vacuum:
> every allowed occupation swap raises \(\langle H\rangle\).
> Two positive PH packets are two masses. Balls that contain
> both have nearly constant enclosed energy; the CHM weight
> still varies.

## Equations

Orthonormal occupied \(L_\pm\) and unoccupied \(R_\pm\),

\[
C=C_0+\alpha\bigl(\Delta_++\Delta_-\bigr),\qquad
\Delta=\lvert R\rangle\langle R\rvert-\lvert L\rangle\langle L\rvert.
\]

Sources \((5,6,6)\) and \((7,6,6)\). Same \(P_{\mathrm{CHM}}\),
\(P_{\mathrm{flat}}\) as Paper 37. Subset: both sources
inside the ball (\(n=11\)).

## Verdicts

Instrument holds. \(C\) in \([0,1]\). \(\rho(K_{\mathrm{vac}})=0.998\).

| Set | \(R_{\mathrm{CHM}}\) | \(R_{\mathrm{flat}}\) | C2 |
| --- | --- | --- | --- |
| all \(512\) | \(0.044\) | \(0.159\) | PASS |
| both-inside \(11\) | \(0.107\) | \(0.115\) | PASS (gap \(0.008\)) |

On the subset, \(\mathrm{std}(P_{\mathrm{flat}})=0.034\),
\(\mathrm{std}(P_{\mathrm{CHM}})=0.098\). Flat is the
near-constant Gauss predictor. \(\delta S\) still tracks
that small enclosed-energy variation (\(\rho_{\mathrm{flat}}=0.993\)).

Auditor \(N=10\), \(n_{\mathrm{both}}=11\): all-ball C2
**CONFIRMED** (\(0.054<0.158\)). Both-inside C2b **REFUTED**
(\(0.151>0.102\), flat wins).

`PAIR_CHM_BEATS_GAUSS` on the pre-registered all-ball gate.
The subset that can kill Gauss does **not** robustly select
CHM. The all-ball win is mostly the pair entering or leaving
the ball (Paper 37 again). Not \(8\pi G\). Not Einstein.

## Equation-to-code

| Object | Where |
| --- | --- |
| Pair, gates, subset | `scripts/m9_34_dipole.py` |
| Adversary | `scripts/m9_34_audit_dipole.py` |

Paper: [`../latex/44_Two_Source.tex`](../latex/44_Two_Source.tex).
