# M9.59: two-mass direction from \(\nabla_c M_{AB}\)

> Paper 68 fixed \(\lvert g\rvert=M/A\). The direction was
> \(\hat r\). Enclosure selects one direction: the
> steepest rise of the *pair* mass. The state is the
> orthonormal two-source transfer (M9.40).

## Equations

\[
\hat n(c)=\frac{\nabla_c M_{AB}(c)}{\lvert\nabla_c M_{AB}(c)\rvert},\qquad
\mathbf g(c)=-\frac{M_{AB}(c)}{A(c)}\,\hat n(c).
\]

\(\nabla_c\) is the two-sided lattice difference of
the mass map (exact on this graph). Newton comparison,
not a hop force:

\[
\mathbf g_N(c)=M_A\frac{\mathbf A-\mathbf c}{\lvert\mathbf A-\mathbf c\rvert^3}
+M_B\frac{\mathbf B-\mathbf c}{\lvert\mathbf B-\mathbf c\rvert^3}.
\]

eigh is float64 LAPACK. Masses, gradients, Newton
vectors, angles: `mpmath` \(dps=50\). Floor
\(\lvert\nabla M\rvert>10^{-6}M_{AB}\).

## Verdicts

Unequal pair \(N=12\), \(\alpha_A=0.02\),
\(\alpha_B=0.04\). \(M_B/M_A=2.00\).
Additivity \(\lvert M_{AB}-M_A-M_B\rvert/M_{AB}=6\times 10^{-8}\).
62 leaking centres.

| Comparison | median angle |
| --- | --- |
| \(\hat n\) vs \(\mathbf g_N\) | \(10.75^\circ\) |
| \(\hat n\) vs CM | \(62.17^\circ\) |

\(\lvert\nabla M\rvert\) at the \(1/r^2\) null \(0.012\),
at the CM \(0.080\). Equal pair: midpoint
\(\lvert\nabla M\rvert/\mathrm{median}=2.5\times 10^{-4}\).

C_dir, C_notcm, C_null, C_mid, C_add **PASS**.

Auditor \(N=10\): C_dir **CONFIRMED**
(\(2.55^\circ\)). C_notcm **CONFIRMED**
(\(2.55^\circ<3.38^\circ\)). C_null **REFUTED**
(only two axis sites with a two-sided gradient).

`PAIR_DIRECTION`. *computed.* The field used on
an entangled pair is \(g=-(M/A)\hat n\) with
\(\hat n=\nabla_c M_{AB}\). It tracks \(M/r^2\)
superposition, not the CM. Not Einstein.

## Equation-to-code

| Object | Where |
| --- | --- |
| Pair \(M_{AB}\), lattice \(\nabla_c\), mpmath gates | `scripts/m9_59_direction.py` |
| Adversary | `scripts/m9_59_audit_direction.py` |

Paper: [`../latex/69_Direction.tex`](../latex/69_Direction.tex).
