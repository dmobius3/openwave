# M9.19: \(S^3\) with a curvature axis

> Geometry of a guess. Wick sign-flip is *derived*. Virtual
> particles living there is *conjectured*, not a habitat.

## Equations

Round \(S^3\) of radius \(\rho\), hyperspherical:

\[
\mathrm{d}s^2=\rho^2\bigl(\mathrm{d}\chi^2+\sin^2\chi\,\mathrm{d}\Omega_2^2\bigr),
\qquad\chi\in[0,\pi].
\]

Ricci scalar from the metric (Christoffel, this script):

\[
\mathcal{R}=\frac{6}{\rho^2}.
\]

Imaginary radius \(\rho=i\ell\):

\[
\mathcal{R}\mapsto\frac{6}{(i\ell)^2}=-\frac{6}{\ell^2}.
\]

That is \(H^3\) (Euclidean AdS\(_3\)). One 3-sphere family, two
signs of curvature. *derived.*

Equator \(\chi=\pi/2\) (\(\sin\chi=1\)) is a round \(S^2\) of
radius \(\rho\). That is the same 2-sphere topology as the
diamond waist of Papers 23--27. *derived* as an induced metric,
not as a derivation of those papers.

Embedding coordinate \(X_4=\cos\chi\). Haar mean on \(S^3\):

\[
\langle\cos\chi\rangle
=\frac{\int_0^\pi\cos\chi\sin^2\chi\,\mathrm{d}\chi}
{\int_0^\pi\sin^2\chi\,\mathrm{d}\chi}=0.
\]

*derived.* If one *identifies* \(X_4\) with curvature, the
uniform mean is zero. The identification is author ontology,
not a theorem. It is the mean-zero foam *given* that map.

Isometries: \(\dim\mathfrak{so}(4)=6\),
\(\dim\mathfrak{so}(1,4)=10\), \(\dim\mathfrak{so}(2,4)=15\).
This \(S^3\) does not give the CHM net of balls in dS\(_4\).
Paper 17 is unchanged.

## Gates

| ID | Claim | Result |
| --- | --- | --- |
| C1 | \(\mathcal{R}=6/\rho^2\) from the metric | PASS |
| C2 | \(\rho\to i\ell\) flips the sign | PASS, auditor CONFIRMED |
| C3 | equator is round \(S^2\) | PASS |
| C4 | Haar \(\langle\cos\chi\rangle=0\) | PASS, auditor CONFIRMED |

`S3_WICK_IDENTITIES_PASS`.

## Tags

| Statement | Tag |
| --- | --- |
| Ricci, Wick sign, equator \(S^2\), Haar mean | *derived* / *computed* |
| Virtual particles live on this \(S^3\) | *guess* |
| \(X_4\) *is* curvature | *guess* |
| The two caps are dS / AdS *spacetime* | *not claimed* (this is 3-geometry) |
| FGHMV in dS / selected \(\Lambda\) | *not claimed* |

## Equation-to-code

| Object | Where |
| --- | --- |
| Ricci from \(g\) | `scripts/m9_19_s3_curvature_axis.py` `ricci_scalar_s3` |
| Wick + Haar | same; audit `m9_19_audit_s3.py` |

Paper: [`../latex/28_S3_Curvature_Axis.tex`](../latex/28_S3_Curvature_Axis.tex).
