# M9.26: point source separates the kernels --- flat wins

> Paper 34 could not tell CHM from flat (smooth Gaussian).
> A single-site potential can. The enclosed-energy kernel wins.

## Setup

Same \(512\) balls of radius \(2\) on \(N=12\). \(V=\varepsilon\)
at \((6,6,6)\) only. \(n_{\mathrm{occ}}\) stable.

## Verdicts

| Gate | Solver | Auditor \(N=10\) |
| --- | --- | --- |
| C1 linear | PASS \(1-8\times 10^{-10}\) | --- |
| C2 CHM \(<\) flat | **FAIL** \(0.885>0.559\) | **REFUTED** \(0.895>0.639\) |
| C4 \(\lvert\rho_{\mathrm{CHM}}\rvert>0.60\) | FAIL \(0.465\) | FAIL \(0.446\) |
| \(\rho_{\mathrm{flat}}\) | \(-0.829\) | --- |

`POINT_KERNEL_NOT_CHM`. When the source is sharp enough to
discriminate, \(\delta S\) tracks **whether the energy sits in
the ball** (flat / Gauss), not the CHM boost weight
\((R^2-r^2)\).

That is a linear functional of enclosed energy. It is not the
FGHMV / Jacobson kernel. It is not linearized Einstein.

Scripts: `m9_26_point_source.py`, `m9_26_audit_point.py`.
Paper: [`../latex/35_Point_Source_Flat.tex`](../latex/35_Point_Source_Flat.tex).
