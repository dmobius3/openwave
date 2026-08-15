# M9.44: a valid fermion state with uniform \(\delta e\); \(a\propto r\)

> Paper 51's wide Gaussian was not a uniform fluid. On a
> periodic cube the band-edge hop eigenstates are unique
> and flat: \(k=0\) and \((\pi,\pi,\pi)\), \(\lvert\psi\rvert^2=1/V\).

## Equations

Periodic hop. \(L=\) ground state \(E=-6\), \(R=\) top
state \(E=+6\). Occupation transfer at fixed \(H\):

\[
\delta e_i=\alpha(E_R-E_L)/V=12\alpha/V.
\]

Balls use the toroidal (min-image) metric. Continuum
image of the measured-uniform \(\delta e\): constant
\(\rho\) on the Dirichlet interior. Inherited Poisson.

## Verdicts

\(N=12\), \(\alpha=0.02\). \(E_L=-6\), \(E_R=+6\).
\(\mathrm{std}(\delta e)/\mathrm{mean}=7\times 10^{-13}\).

| Gate | Result |
| --- | --- |
| C_unif | PASS |
| C_fl | PASS \(\rho(\delta S,P_{\mathrm{flat}})=0.999\) |
| C_grow | PASS, factor \(20.6\) |
| C_dens | PASS, \(P/V\) IQR \(\sim 0\), \(P/A\) IQR \(0.36\) |
| C_lin PRIMARY | PASS, slope \(1.029\), \(a_r<0\) |
| C_invsq | FAIL as required (\(\lvert\alpha+2\rvert=3.03\)) |

Auditor \(N=10\): C_unif **CONFIRMED**. C_lin **CONFIRMED**
(same slope: Poisson of uniform \(\rho\)). C_fl **REFUTED**
(\(\rho=0.047\)): wrapping \(R=4\) balls on \(N=10\) are
most of the torus. The volume law is not robust there.

`UNIFORM_NEWTON_LAMBDA`. *computed.* A valid fermion
state has uniform energy. The first law is extensive on
the primary lattice. Inherited Newton of that density is
\(a\propto r\), not \(1/r^2\). That is the Newtonian
\(\Lambda\) signature of this state. Not a derived
Poisson. Not a de Sitter dual. No `MODELS.md`.

## Equation-to-code

| Object | Where |
| --- | --- |
| PBC edges, first law, Poisson | `scripts/m9_44_uniform_pbc.py` |
| Adversary | `scripts/m9_44_audit_pbc.py` |

Paper: [`../latex/54_Newton_Lambda.tex`](../latex/54_Newton_Lambda.tex).
