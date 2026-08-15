# M9.43: \(\kappa\) survives \(3+1\)D mass

> Papers 46--52 measured \(\kappa\) on a massless cube hop.
> The diamond waist of Papers 23--25 is the geodesic ball
> of a \(3+1\)D staggered fermion. This note turns the
> staggered mass on and asks whether \(\kappa\) moves.

## Equations

Spatial Hamiltonian of a \(3+1\)D staggered fermion
(Paper 23):

\[
H_{ii}=m\,(-1)^{x+y+z},\qquad H_{\langle ij\rangle}=-1.
\]

Same occupation transfer as Paper 37, fixed \(H\).
Source-centered balls. \(R_{\mathrm{enc}}\) is the
smallest \(R\) with \(P_{\mathrm{flat}}/M_{\mathrm{global}}>0.95\).

\[
\kappa(m)=\delta S(R_{\mathrm{enc}})/P_{\mathrm{flat}}(R_{\mathrm{enc}}).
\]

## Verdicts

\(N=12\), packet \((6,6,6)\). Still a star (growth
\(1.13\to 1.12\)). \(R_{\mathrm{enc}}=3\) at every \(m\).
\(\kappa(R=4)/\kappa(R=3)\) agrees to \(0.03\%\).

| \(m\) | \(\kappa\) | \(\lvert\kappa/\kappa(0)-1\rvert\) |
| --- | --- | --- |
| \(0\) | \(1.0473\) | --- |
| \(0.25\) | \(1.0458\) | \(0.14\%\) |
| \(0.50\) | \(1.0410\) | \(0.60\%\) |

C_univ **PASS**. Auditor \(N=10\), \(m=0\) vs \(0.40\):
\(0.975\) vs \(0.977\) (\(0.21\%\)), **CONFIRMED**.

Single-ball \(\lvert\delta S-\mathrm{Tr}(K_{\mathrm{vac}}\Delta C)\rvert/\lvert\delta S\rvert\)
grows with \(R\). That is not Paper 37's many-ball Pearson.
Recorded, not a gate.

`KAPPA_SURVIVES_MASS`. *computed.* The first-law constant
is not a massless-cube artifact. Not a continuum proof.
Not \(8\pi G\). Not de Sitter.

## Equation-to-code

| Object | Where |
| --- | --- |
| Staggered \(H\), scan, gates | `scripts/m9_43_diamond_kappa.py` |
| Adversary | `scripts/m9_43_audit_diamond.py` |

Paper: [`../latex/53_Diamond_Kappa.tex`](../latex/53_Diamond_Kappa.tex).
