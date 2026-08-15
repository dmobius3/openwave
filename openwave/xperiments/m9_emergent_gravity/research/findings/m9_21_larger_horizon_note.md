# M9.21: larger-horizon first law --- instrument dies

> Tried \(R=5\) so CHM vs linear could split. The occupancy-stable
> hop probe does not survive. Not a Planck-scale result.

## What was locked

\(N=14\), \(R=5\) (margin 1). Surface subsample: three
\(r_{\mathrm{mid}}\) bins, 10 lex-first bonds each (\(30\) hops).
Same C3/C4 as M9.20. First \(\varepsilon=0.01\), then
\(\varepsilon=0.002\) after \(27/30\) occupancy flips.

## Result

| \(\varepsilon\) | kept | flipped | C4 |
| --- | --- | --- | --- |
| \(0.01\) | \(3\) | \(27\) | \(n=3\), not scored |
| \(0.002\) | \(5\) | \(25\) | linear slightly better on \(n=5\) |

\(n_{\mathrm{occ}}(0)=1370\) on a \(14^3=2744\) grid. The sea is
one or two states from half filling. A single surface hop of
size \(0.002\) crosses the Fermi level on the locked lex sample.

`LARGER_HORIZON_INSTRUMENT_REJECTED`. C4 at \(R=5\) is **not
scored**. Paper 29 (\(R=4\), \(294\) hops, \(n_{\mathrm{flip}}=0\),
C4 tie) remains the last valid uniqueness test.

This is not \(a\to\ell_P\). It is a more degenerate Fermi sea
on a larger ball.

## Equation-to-code

`scripts/m9_21_larger_horizon.py`.

Paper: [`../latex/30_Larger_Horizon_Instrument.tex`](../latex/30_Larger_Horizon_Instrument.tex).
