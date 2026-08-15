# M9.54: enclosure of \(2h(\alpha)\) survives the sea

> Paper 63: on a compact packet, \(f_S\) tracks \(f_E\).
> Nested balls cannot decide the sea (\(A,V\) collinear,
> Paper 42). Vacuum \(S\) is area-law (Paper 57).
> This run uses the PBC band-edge transfer and slabs,
> whose area does not grow with thickness.

## Equations

Periodic hop. \(L=\) ground state, \(R=\) top state.

\[
f_S=\frac{\delta S}{2h(\alpha)},\qquad
f_E=\frac{P_{\mathrm{flat}}}{M_{\mathrm{global}}}=\frac{V}{N^3}.
\]

C_slab PRIMARY: \(\mathrm{grow}=\delta S(t{=}3)/\delta S(t{=}1)\)
closer to \(3\) than to \(1\), and \(\mathrm{grow}>1.5\).

## Verdicts

\(N=12\). Nine regions: balls \(R=2,3,4\), cubes
side \(2,3,4\), slabs \(t=1,2,3\). Slab area is
\(288\) at every \(t\).

| \(\alpha\) | grow | \(\rho(f_S,f_E)\) | RMS | \(\rho(\delta S,V)\) | \(\rho(\delta S,A)\) |
| --- | --- | --- | --- | --- | --- |
| \(0.01\) | \(2.59\) | \(0.994\) | \(0.031\) | \(0.994\) | \(0.867\) |
| \(0.02\) | \(2.55\) | \(0.994\) | \(0.036\) | \(0.994\) | \(0.872\) |
| \(0.04\) | \(2.51\) | \(0.995\) | \(0.043\) | \(0.995\) | \(0.877\) |

C_slab, C_rho, C_rms **PASS**. Full-system
\(\Delta S=2h(\alpha)\) to \(10^{-10}\). Grow is
\(2.5\), not \(3\): slabs sit above \(f_E\)
(residual, not a second law). Auditor \(N=10\):
grow \(2.56,2.47\), \(\rho=0.996\). All
**CONFIRMED**.

`SEA_ENCLOSURE`. *computed.* Paper 63 is not a
packet artifact. Sea-transfer \(\delta S\) tracks
enclosed energy, not area. Vacuum \(S\) (Paper 57)
is a different object. Gravity is still
\(\sum\delta e\) plus inherited Gauss. Not Clausius.

## Equation-to-code

| Object | Where |
| --- | --- |
| Sea transfer, slabs, cubes, balls | `scripts/m9_54_sea_enclose.py` |
| Adversary | `scripts/m9_54_audit_sea.py` |

Paper: [`../latex/64_Sea_Enclosure.tex`](../latex/64_Sea_Enclosure.tex).
