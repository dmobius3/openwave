# M9.22: C4 scored at \(R=5\) with half-filling

> M9.21 left C4 unscored because \(E<0\) filling flickered.
> Occupy the lowest \(V/2\) states. Occupancy cannot flip.

## Result

`C4_TIE`. Scored. \(n=30\), \(n_{\mathrm{flip}}=0\).

| Kernel | \(R_{\mathrm{shape}}\) | \(\rho\) |
| --- | --- | --- |
| CHM \(R^2-r^2\) | \(0.99951\) | \(-0.031\) |
| linear \(R-r\) | \(0.99945\) | \(-0.033\) |
| flat | \(0.99542\) | \(+0.096\) |

Gap CHM vs linear: \(5.8\times 10^{-5}\) (lock: \(<0.005\) is a tie).
Neither tracks \(\delta S\) (\(\lvert\rho\rvert\sim 0\)). Flat is
slightly better. C3 FAIL. C0/C1/C2 FAIL.

Auditor \(N=16\), \(R=5\), \(15\) hops, half-fill, \(n_{\mathrm{flip}}=0\):
\(R_{\mathrm{CHM}}=0.889\), \(R_{\mathrm{lin}}=0.892\),
\(R_{\mathrm{flat}}=0.997\). Same story: C4 is a dead heat;
C3 beats flat; tracking \(\lvert\rho\rvert=0.46<0.60\).

## Verdict

C4 is **scored** at \(R=5\). The parabola is **not selected**
over a linear weight. First-law tracking collapses under
half-filling (the sea is no longer the \(E<0\) vacuum of
Papers 26--29). Paper 29 remains the last *tracking* C4
(a \(0.002\) tie at \(R=4\), \(E<0\), \(294\) hops).

Not Planck. Not Einstein.

Scripts: `m9_22_halffill_horizon.py`, `m9_22_audit_horizon.py`.
Paper: [`../latex/31_C4_Scored_HalfFill.tex`](../latex/31_C4_Scored_HalfFill.tex).
