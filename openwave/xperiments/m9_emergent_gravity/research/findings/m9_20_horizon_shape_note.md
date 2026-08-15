# M9.20: horizon first law, CHM vs linear vs flat

> Paper 27: CHM beat flat on the cut. This run asks whether the
> *parabola* is selected, with occupancy held fixed.

## Gates

Surface \(r_{\mathrm{mid}}\ge 0.75R\). Drop hops that flip
\(n_{\mathrm{occ}}\). Three kernels on the kept set: CHM
\(R^2-r^2\), linear \(R-r\), flat. Score \(\delta S\) vs
\(\mathrm{Tr}(K\Delta C)\).

| ID | Claim | Solver \(N=12,R=4\) |
| --- | --- | --- |
| C0 | 1d \(\lvert\rho_{\mathrm{CHM}}\rvert>0.70\) | FAIL, \(0.591\) |
| C1 | 3d \(\lvert\rho\rvert>0.60\) | FAIL, \(0.589\) |
| C2 | \(R_{\mathrm{CHM}}<0.70\) | FAIL, \(0.808\) |
| C3 | \(R_{\mathrm{CHM}}<R_{\mathrm{flat}}\) | **PASS**, \(0.808<0.994\) |
| C4 | \(R_{\mathrm{CHM}}<R_{\mathrm{lin}}\) | **tie**, \(0.8082<0.8099\) |

\(n=294\), \(n_{\mathrm{flip}}=0\). Occupancy is stable on this
even cube.

## Auditor \(N=14,R=3\)

\(n_{\mathrm{flip}}=142\), only \(8\) hops kept. Sample is too
small. C4 **FAIL** (linear \(0.9589\) vs CHM \(0.9594\)). C3
still CHM \(<\) flat. Not a C4 confirmation.

## Verdict

`HORIZON_CHM_NOT_UNIQUE`. C3 replicates Paper 27 on a new radius
with no occupancy flips: a radial CHM-type weight beats a flat
kernel on the cut. C4 does not select the parabola over a linear
envelope. Tracking floors still fail.

*computed.* Not \(\eta=1/4G\). Not Einstein. Not the Bloch guess.

## Equation-to-code

`scripts/m9_20_horizon_shape.py`, `scripts/m9_20_audit_horizon.py`.

Paper: [`../latex/29_Horizon_Not_Unique.tex`](../latex/29_Horizon_Not_Unique.tex).
