# M9.13: A2 on a 4d-continuum diamond waist

> Same locked gate as M9.10: \(R(m)/R(0)<2\) for \(0<mL\le 8\).
> Physical theory is \(3+1\)D. The waist is a geodesic ball.
> Two lattice spacings at matched \(mL\). Not \(a\to 0\).

## Equations

Spatial Hamiltonian of a free staggered fermion in three-space
(the equal-time slice of a \(3+1\)D field):

\[
H_{ii}=m\,(-1)^{x+y+z},\qquad
H_{\langle ij\rangle}=-1
\]

on an \(N^3\) open grid. The causal-diamond waist is the integer ball

\[
(x-c)^2+(y-c)^2+(z-c)^2\le R^2,\qquad c=N/2,\qquad L=2R.
\]

Peschel modular Hamiltonian on the ball, from the occupied projector
\(C=\sum_{E<0}v_E v_E^\top\) restricted to the ball:

\[
K=\log\frac{1-C_A}{C_A}.
\]

Local ansatz: \(K_{\mathrm{loc}}\) keeps only the diagonal and the six
spatial nearest neighbours inside the ball. Remainder

\[
R(m)=\frac{\lVert K-K_{\mathrm{loc}}\rVert_F}{\lVert K\rVert_F}.
\]

Continuum identities *cited*, not derived here. Rindler wedge in
Minkowski, any dimension including \(3+1\) (Bisognano--Wichmann 1975):

\[
K=2\pi\int_{x^1>0}x^1\,T_{00}(x)\,d^3x.
\]

CFT ball (Casini--Huerta--Myers 2011):

\[
K=2\pi\int_{|x|<R}\frac{R^2-r^2}{2R}\,T_{00}(x)\,d^3x.
\]

Both are local integrals of \(T_{00}\). That is the continuum meaning
of A2 for those regions. This script tests the *massive* free-fermion
deformation of the ball, which is not a CFT.

## Pre-registered gates

| ID | Claim | Pass |
| --- | --- | --- |
| C1 | \(0<R(0)<1\) on each spacing | both |
| C2 | \(R(m)/R(0)<2\) for \(0<mL\le 8\) on each spacing | both |
| C_CONT | C2 on \((N,R)=(12,3)\) and \((16,5)\) at matched \(mL\in\{3,6\}\) | both |
| C4 | \(C_A\) spectrum in \((0,1)\); \(K\) Hermitian | both |

## Verdicts

| Spacing | \(n_{\mathrm{ball}}\) | \(R(0)\) | ratios at \(mL=3,6\) | C2 |
| --- | --- | --- | --- | --- |
| \(N=12\), \(R=3\) | 123 | 0.132 | 1.36, 1.33 | PASS |
| \(N=16\), \(R=5\) | 515 | 0.158 | 1.13, 1.19 | PASS |
| auditor \(N=14\), \(R=4\) | --- | 0.150 | 1.19, 1.22 | CONFIRMED |
| adversary \(N=13\), \(R=3\) | 123 | 0.139 | worst \(1.29\) at \(mL=4\) | CONFIRMED |
| adversary \(N=15\), \(R=4\) | 257 | 0.157 | worst \(1.16\) at \(mL=5\) | CONFIRMED |

`A2_DIAMOND_4D_PASS`. The finer spacing has the *smaller* ratios
(*computed*, not a continuum proof). Author ontology that the field
modes live in \(3+1\)D is the domain of the test, not a result.

Adversary catch, not a C2 fail: \(R(m)\) is **not monotone**. It
rises from \(m=0\), peaks near \(mL\sim 4\)--\(5\), then falls.
``Mass localizes \(K\)'' is false on this stencil until
\(mL\gtrsim 5\). The factor-of-two cut is doing work; the peak
stayed at \(1.29\).

## Equation-to-code

| Object | Where |
| --- | --- |
| \(H\), ball, \(K\), \(R\) | `scripts/m9_13_A2_diamond_4d.py` `staggered_H_3d`, `ball_sites`, `run_resolution` |
| Independent \(R\) | `scripts/m9_13_audit_A2.py` `diamond_R` |

## Tags

| Statement | Tag |
| --- | --- |
| C2 on this diamond waist, two spacings + auditor | *computed* (threshold 2, locked) |
| Ratios decrease from coarse to fine | *computed* (two points) |
| BW / CHM locality for wedges / CFT balls in \(3+1\)D | *cited* |
| \(a\to 0\) continuum theorem for massive Dirac | *not claimed* |
| A2 for the Standard Model | *unresolved* |
| Entanglement selects de Sitter / a value of \(\Lambda\) | *not claimed* |
| FGHMV in de Sitter | *not claimed* (Paper 17) |

## What this is not

Not a 4-spatial-dimensional lattice. Not the Standard Model. Not
\(\eta=1/4G\). Not Jacobson 2016 as \([P]\). Not FGHMV in de Sitter.
A two-spacing check is not \(a\to 0\).

Scripts: `m9_13_A2_diamond_4d.py`, `m9_13_audit_A2.py`.
Data: `data/m9_13_A2_diamond_4d.json`, `data/m9_13_audit_A2.json`.
Paper: [`../latex/23_A2_4d_Continuum_Diamond.tex`](../latex/23_A2_4d_Continuum_Diamond.tex).
