# M9.15: CHM shape of the modular hopping kernel

> A2 said \(K\) is local. This test asks whether that local kernel
> is the Casini--Huerta--Myers envelope. It is.

## Equations

CHM / Bisognano--Wichmann (cited):

\[
K=2\pi\int_A w(x)\,T_{00}(x),\qquad
w=R^2-r^2
\quad\text{(ball)},\qquad
w=(L/2)^2-x^2
\quad\text{(interval)}.
\]

\(T_{00}\) contains the hop. The lattice statement is therefore
about nearest-neighbour entries of the Peschel kernel
\(K=\log((1-C_A)/C_A)\):

\[
K_{ij}\;\text{tracks}\;w_{ij}=\tfrac12(w_i+w_j)
\quad\text{on spatial nearest neighbours}.
\]

Hops are \(-1\), so the Pearson correlation \(\rho(K_{ij},w_{ij})\)
must be **negative**. Residual after \(K_{ij}\approx a\,w_{ij}+b\):

\[
R_{\mathrm{shape}}=\frac{\lVert K_{\mathrm{nn}}-(aw+b)\rVert}{\lVert K_{\mathrm{nn}}-\langle K_{\mathrm{nn}}\rangle\rVert}.
\]

A first instrument (diagonal \(K_{ii}\) vs \(w\)) failed the 1d
theorem and was rejected. \(r^2\) is affine to \(-w\) and cannot
mutate the test.

## Pre-registered gates

| ID | Claim | Pass |
| --- | --- | --- |
| C0 | 1d massless, \(\rho<-0.70\) | yes, \(-0.942\) |
| C1 | 3+1D ball \(R=5\), \(m=0\), \(\rho<-0.60\) | yes, \(-0.987\) |
| C2 | \(R_{\mathrm{shape}}<0.50\) at \(m=0\) | yes, \(0.162\) |
| C3 | permuted \(w\) has \(\lvert\rho\rvert<0.30\) | yes, \(-0.076\) |
| C4 | C1 at \(mR=0.5\) | yes, \(-0.989\) |

## Verdicts

| Object | Value | Tag |
| --- | --- | --- |
| 1d \(\rho\), \(R_{\mathrm{shape}}\) | \(-0.942\), \(0.334\) | *computed* |
| 3+1D \(m=0\) \(\rho\), \(R_{\mathrm{shape}}\) | \(-0.987\), \(0.162\) | *computed* |
| auditor 1d \(N=200\), \(L=24\) | \(\rho=-0.984\) | CONFIRMED |
| auditor 3d \(N=14\), \(R=4\) | \(\rho=-0.987\), \(R_{\mathrm{shape}}=0.161\) | CONFIRMED |

`CHM_SHAPE_PASS` on the locked gates. The local modular hop
tracks a bulk-peaked radial envelope on this free \(3+1\)D
fermion. Combined with Papers 23--24 this is the lattice
content of Jacobson's input: a local \(K\) with a geometric
weight, and an area-law \(S\).

Adversary catch, QUALIFIED not REFUTED: C1/C2 hold on
independent sizes (\(\rho\in[-0.97,-0.999]\),
\(R_{\mathrm{shape}}\le 0.27\)). The locked quadratic is
\emph{not unique} on small balls: a linear \(R-r\) weight can
beat \(R^2-r^2\) at \(R=3\). The parabola wins at \(R=4,5\).
The hop tracks a CHM-\emph{type} envelope, not a uniquely
identified quadratic. Hop \(=+1\) flips the sign of \(\rho\)
and is not the staggered Hamiltonian of this column.

Conditional identification, not a measurement: if
\(S=\alpha A_{\mathrm{cut}}=A_{\mathrm{phys}}/(4G)\) and
\(A_{\mathrm{cut}}=A_{\mathrm{phys}}/a^2\), then
\(G=a^2/(4\alpha)\). With Paper 24's \(\alpha=0.245\),
\(G/a^2\approx 1.02\). That sets the cutoff equal to a Planck
length. It does not measure \(G\).

## Equation-to-code

| Object | Where |
| --- | --- |
| 1d / 3d \(K\), \(\rho\), \(R_{\mathrm{shape}}\) | `scripts/m9_15_chm_shape.py` |
| Independent sizes | `scripts/m9_15_audit_chm.py` |

## What this is not

Not \(\eta=1/4G\) as a derived number. Not mean-zero foam.
Not a selection of \(\Lambda\). Not FGHMV in de Sitter.
Not the Standard Model. Not a continuum \(a\to 0\) proof.

Scripts: `m9_15_chm_shape.py`, `m9_15_audit_chm.py`.
Data: `data/m9_15_chm_shape.json`, `data/m9_15_audit_chm.json`.
Paper: [`../latex/25_CHM_Shape_Diamond.tex`](../latex/25_CHM_Shape_Diamond.tex).
