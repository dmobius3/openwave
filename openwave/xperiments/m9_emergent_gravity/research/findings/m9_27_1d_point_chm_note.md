# M9.27: 1d point-source calibration --- instrument fails where CHM is a theorem

> Paper 35: in 3d, flat / enclosed energy beats CHM.
> That is physics only if the same instrument recovers CHM in 1d,
> where Casini--Huerta--Myers is a theorem for \(K\) itself.
> It does not.

## Equations

Massless open chain, hop \(-1\). Interval \(I\) of length \(L\).
CHM weight on the interval (affine to \(x(L-x)\))

\[
w_i=\Bigl(\frac{L}{2}\Bigr)^2-x_i^2.
\]

Site energy and the two predictors (same operators as Paper 35)

\[
e_i=\sum_j H_{ij}C_{ij},\qquad
P_{\mathrm{CHM}}=\sum_{i\in I}w_i\,\delta e_i,\qquad
P_{\mathrm{flat}}=\sum_{i\in I}\delta e_i.
\]

Peschel entropy and modular Hamiltonian of a block

\[
S=-\mathrm{Tr}\bigl[C\log C+(1-C)\log(1-C)\bigr],\qquad
K=\log\frac{1-C}{C}.
\]

Exact differential: \(\mathrm{d}S=\mathrm{Tr}(K\,\mathrm{d}C)\).
Diagnostics, not gates:

\[
\delta S\;\stackrel{?}{=}\;\mathrm{Tr}(K_{\mathrm{vac}}\Delta C)
\qquad\text{vs}\qquad
\delta S\;\stackrel{?}{=}\;\mathrm{Tr}(K_{\mathrm{mid}}\Delta C).
\]

The vacuum first law that CHM would turn into \(\int w\,\delta\langle T_{00}\rangle\)
is the left-hand side.

## Pre-registered gates

Chain \(N=200\), \(L=16\), source site \(100\), \(\varepsilon=0.05\)
and \(0.10\). \(E<0\) occupancy; half-fill only if \(n_{\mathrm{occ}}\)
flips (it did not: \(100/100/100\)). \(185\) intervals.

| Gate | Lock | Result |
| --- | --- | --- |
| C0 | \(\max\|\delta S\|>10^{-6}\) | PASS \(3.30\times 10^{-4}\) |
| C1 | Pearson\(\delta S(\varepsilon),\delta S(2\varepsilon)>0.95\) | PASS \(1-4\times 10^{-9}\) |
| C2 PRIMARY | \(R(\delta S,P_{\mathrm{CHM}})<R(\delta S,P_{\mathrm{flat}})\) | **FAIL** \(0.999>0.961\) |
| C4 | \(\lvert\rho(\delta S,P_{\mathrm{CHM}})\rvert>0.60\) | **FAIL** \(-0.054\) |

\(\rho_{\mathrm{flat}}=-0.276\). \(\delta e\) is localized
(\(\mathrm{rms}\) width \(0.92\) sites). Verdict
`1D_INSTRUMENT_OR_FLAT`.

## Diagnostics (not gates)

| Object | Solver | Auditor \(N=160\), \(L=12\) |
| --- | --- | --- |
| \(\rho(\delta S,P_{\mathrm{CHM}})\) | \(-0.054\) | \(-0.060\) |
| \(\rho(\delta S,\mathrm{Tr}(K_{\mathrm{vac}}\Delta C))\) | \(-0.025\) | \(-0.016\) |
| \(\rho(\delta S,\mathrm{Tr}(K_{\mathrm{mid}}\Delta C))\) | \(0.999999\) | \(0.999996\) |
| rel.\(\lvert\delta S-\mathrm{Tr}(K_{\mathrm{mid}}\Delta C)\rvert\) | \(0.14\%\) | \(0.27\%\) |

Five other local \(\delta e\) identifications (fixed \(H_0\),
bond-symmetric, \(\varepsilon n_i\) only, kinetic bonds) all
fail C2 the same way. The failure is not the one-sided
assignment \(e_i=\sum_j H_{ij}C_{ij}\).

## Adversarial audit

Own script, no import, different size and source
(`m9_27_audit_1d.py`). C2 **REFUTED**. Vacuum first law
**REFUTED**. Midpoint identity **CONFIRMED**.

## What this is

Paper 25 still stands: the *shape of \(K\)* on hops is CHM.
This note is about a different object: \(\delta S\) of a
sliding interval versus a local functional of \(\delta e\),
after a point Hamiltonian source.

That comparison already fails in 1d. \(\delta S\) is the
path integral of \(K(C)\), not \(\mathrm{Tr}(K_{\mathrm{vac}}\Delta C)\).
A point potential inside the interval rotates the entanglement
spectrum enough that the vacuum modular energy is not \(\delta S\).
Paper 35 therefore cannot be cashed as ``3d / non-CFT physics
selects Gauss's law instead of Einstein.'' It is the same
instrument, failing the theorem case.

Not \(8\pi G\). Not de Sitter. Not a claim that CHM is false.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(H\), \(C\), \(e_i\), \(S\), \(K\) | `scripts/m9_27_1d_point_chm.py` (`hop_H`, `occupy`, `site_energy`, `peschel_s`, `peschel_k`) |
| \(P_{\mathrm{CHM}}\), \(P_{\mathrm{flat}}\), gates | same file, `main` |
| Independent sizes | `scripts/m9_27_audit_1d.py` |

Paper: [`../latex/36_1d_Instrument.tex`](../latex/36_1d_Instrument.tex).
