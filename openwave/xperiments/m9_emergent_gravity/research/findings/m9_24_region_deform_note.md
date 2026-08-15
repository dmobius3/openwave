# M9.24: deform the region at fixed \(H\)

> One vacuum \(C\). Vary only the region. Jacobson-style, not
> another hop C4.

## Equations

Peschel \(S\) of \(C\) restricted to region \(A\). Lattice area
\(A_{\mathrm{cut}}\) = bonds leaving \(A\). Fit on balls only:

\[
S=\alpha A_{\mathrm{cut}}+\beta.
\]

C2: cubes have larger RMS about that line than the balls do.

## Verdicts

\(\alpha=0.245\) (same as Paper 24). Shifted \(R=4\) ball:
\(\Delta S/S=6\times 10^{-4}\). C1 PASS.

| Family | RMS about the ball line |
| --- | --- |
| balls \(R=2,3,4,5\) | \(0.27\) |
| cubes \(L=3,4,5,6\) | \(5.67\) |
| taxicab \(t=2,3,4,5\) | \(9.46\) |

C2 PASS. C3 PASS. Auditor \(N=14\): \(0.29\) vs \(3.93\),
CONFIRMED.

Same cut area, different shape: ball \(R=2\) and taxicab \(t=2\)
both have \(A_{\mathrm{cut}}=78\), but \(S=13.60\) vs \(12.00\).
Area law alone cannot do that.

`SHAPE_MATTERS`. *computed.* Not \(\eta=1/4G\). Not Einstein.
The finite piece of \(S\) knows the region's shape. The ball is
the tight family; cubes and diamonds miss the ball line.

Scripts: `m9_24_region_deform.py`, `m9_24_audit_deform.py`.
Paper: [`../latex/33_Region_Shape.tex`](../latex/33_Region_Shape.tex).
