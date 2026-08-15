# M9.14: A1 on a 3+1D diamond waist

> Same logical A1 as M9.11, now an area law. Physical theory is
> \(3+1\)D. Not \(\eta=1/4G\). Not mean-zero curvature.

## Equations

Peschel entropy of the ball \(A\) at the diamond waist:

\[
S=-\mathrm{Tr}\bigl(C_A\log C_A+(1-C_A)\log(1-C_A)\bigr).
\]

Lattice area \(A_{\mathrm{cut}}\) is the number of hopping bonds that
leave the ball. Least-squares ansatz, locked:

\[
S(R)=\alpha\,A_{\mathrm{cut}}(R)+\beta.
\]

Unsubtracted sea energy density on the \(N^3\) grid (diagnostic,
not a gate):

\[
\varepsilon=\frac{1}{N^3}\sum_{E<0}E.
\]

## Pre-registered gates

| ID | Claim | Pass |
| --- | --- | --- |
| C1 | \(\alpha(0)>0\) | yes, \(0.245\) |
| C2 | \(\lvert\alpha(m)-\alpha(0)\rvert/\alpha(0)<0.20\) for \(mR_{\max}\le 0.5\) | yes, max drift \(4.3\%\) |
| C3 | area-law RMSE \(<\) volume-law RMSE at \(m=0\) | yes, \(0.27<2.96\) |
| C4 | at least three radii | yes, \(R=2,3,4,5\) |

## Verdicts

| Object | Value | Tag |
| --- | --- | --- |
| \(\alpha(0)\), \(N=16\) | \(0.245\) | *computed* |
| UV relative drifts | \(0.24\%\), \(1.4\%\), \(4.3\%\) | *computed* |
| auditor \(N=14\), \(\alpha(0)\), rel | \(0.245\), \(2.1\%\) | CONFIRMED |
| adversary \(N=13\), \(\alpha(0)\), worst UV rel | \(0.245\), \(5.9\%\) | C2 CONFIRMED |
| IR diagnostic \(mR_{\max}=2\) | rel \(22.5\%\) (outside the UV window) | *computed* |
| \(\varepsilon(m=0)\) | \(-0.970\) per site | *computed* |

`A1_DIAMOND_4D_PASS` on the **locked** observable \(A_{\mathrm{cut}}\).
The leading piece is a cut-bond area law and is IR-mass stable in
the locked UV window.

Adversary catch, QUALIFIED not REFUTED: C3 is
observable-dependent. Boundary-site count and \(4\pi R^2\) can
lose to a volume fit if \(R=1\) is included. Absolute \(\alpha\)
is proxy-dependent (\(0.245\) cut / \(0.35\) geometric /
\(0.48\) boundary). Relative UV drift is not. The locked area is
cut bonds. Other proxies were not the gate.

## Foam / mean curvature (not a result)

Author target: the quantum foam has mean-zero curvature when the
program is done. That is a destination, not a theorem of this run.

- If "mean zero" is an **IR** statement, it is \(\Lambda=0\)
  (Minkowski). That contradicts "entanglement gravity is de Sitter"
  unless a second scale is named. *unresolved / author ontology.*
- If it is only a statement about **UV fluctuations** \(\langle\delta R\rangle=0\)
  around a smooth background, it does not select \(\Lambda\). *conjectured.*
- This lattice vacuum has \(\varepsilon\simeq -1\) per site. Feeding
  that unsubtracted \(\langle T_{00}\rangle\) to Einstein is the
  cosmological-constant problem (Planck-scale curvature), not
  mean-zero foam. *computed.*
- Normal-ordering \(\langle:T:\rangle=0\) inserts mean-zero stress
  by definition. That is not a derivation. *derived as circularity.*

A1 does not measure curvature and does not finish the foam.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(H\), ball, \(S\), \(\alpha\), \(\varepsilon\) | `scripts/m9_14_A1_diamond_4d.py` |
| Independent \(\alpha\) | `scripts/m9_14_audit_A1.py` `alpha_of` |

## What this is not

Not \(\eta=1/4G\). Not quantum foam. Not mean-zero curvature.
Not Einstein in vacuum. Not a value of \(\Lambda\). Not FGHMV
in de Sitter. Not victory.

Scripts: `m9_14_A1_diamond_4d.py`, `m9_14_audit_A1.py`.
Data: `data/m9_14_A1_diamond_4d.json`, `data/m9_14_audit_A1.json`.
Paper: [`../latex/24_A1_4d_Diamond_Area.tex`](../latex/24_A1_4d_Diamond_Area.tex).
