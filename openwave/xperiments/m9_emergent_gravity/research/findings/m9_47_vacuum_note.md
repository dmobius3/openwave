# M9.47: the Fermi-sea vacuum is area-law. Not \(\Lambda\)

> Paper 56: occupation transfer cannot make negative
> \(\delta e\). The unperturbed sea has \(E_{\mathrm{vac}}<0\).
> If that energy were a first-law volume source, Gauss
> would give outward \(a\propto r\). Vacuum entanglement
> is not a volume source.

## Equations

No transfer. \(C_0\) is the occupied projector.
\(e_i=\sum_j H_{ij}(C_0)_{ij}\), \(E_{\mathrm{vac}}=\sum_i e_i\).
Peschel \(S\) of source-centered balls \(R=2,3,4,5\).
\(A\) = outgoing NN count. \(V\) = site count.
Shape diagnostic \(a_{\mathrm{try}}=-S/R^2\).

## Verdicts

\(N=12\). Open and periodic.

| Vacuum | \(E_{\mathrm{vac}}\) | \(\rho(S,A)\) | \(\rho(S,V)\) | \(S_5/S_2\) | \(A_5/A_2\) | \(V_5/V_2\) | slope \(a_{\mathrm{try}}\) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| open | \(-1657\) | \(0.9999\) | \(0.9957\) | \(8.11\) | \(6.23\) | \(15.61\) | \(0.263\) |
| PBC | \(-1723\) | \(0.9997\) | \(0.9949\) | \(8.08\) | \(6.23\) | \(15.61\) | \(0.262\) |

C_neg, C_area, C_ratio, C_notds **PASS**.
Auditor \(N=10\): \(E_{\mathrm{vac}}<0\) **CONFIRMED**.
Area over volume **CONFIRMED**.

Nested \(A\) and \(V\) are still collinear
(\(\rho>0.99\) both). The growth ratio is the cleaner
cut: \(S\) tracks area, not volume. \(a_{\mathrm{try}}\)
is nearly constant (area / \(R^2\)), not \(a\propto r\).

`VACUUM_AREA_NOT_LAMBDA`. *computed.* The sea's energy
is negative and (on the torus) uniform. That is the
usual cutoff Fermi-sea energy, \(N\)-dependent
(\(-1657\) vs auditor \(-950\)). It is not a measured
\(\Lambda\). Vacuum \(S\) is an area law. The first
law does not promote \(E_{\mathrm{vac}}\) to a
cosmological constant.

## Equation-to-code

| Object | Where |
| --- | --- |
| Open/PBC vacuum, gates | `scripts/m9_47_vacuum.py` |
| Adversary | `scripts/m9_47_audit_vacuum.py` |

Paper: [`../latex/57_Vacuum_Area.tex`](../latex/57_Vacuum_Area.tex).
