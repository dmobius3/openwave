# M9.61: pair masses add. Pair \(\mathbf g\)s do not.

> Papers 69--70: \(\mathbf g=-(M/A)\hat n\),
> \(\hat n=\nabla_c M/\lvert\nabla_c M\rvert\).
> For two packets that field is not linear:

\[
\mathbf g_{AB}\propto M_{AB}\frac{\nabla M_A+\nabla M_B}{\lvert\nabla M_A+\nabla M_B\rvert},
\quad
\mathbf g_A+\mathbf g_B\propto
M_A\hat n_A+M_B\hat n_B.
\]

They agree iff the two gradients are parallel.

Exact open-hop basis. `mpmath` \(dps=80\). No LAPACK.

## Verdicts

\(N=12\), unequal pair, \(64\) centres with all
three gradients.

| Gate | Solver | Auditor \(N=11\) |
| --- | --- | --- |
| C_map \(\lvert M_{AB}-M_A-M_B\rvert\) | **PASS** \(1.9\times 10^{-7}\) | **CONFIRMED** \(6\times 10^{-7}\) |
| C_ang PRIMARY \(\angle(\mathbf g_{AB},\mathbf g_A+\mathbf g_B)\) | **PASS** \(5.43^\circ\) | **CONFIRMED** \(3.15^\circ\) |
| C_rel relative vector residual | **FAIL** \(0.255\) | **REFUTED** \(0.454\) |
| C_far far centres, \(\lvert c-A\rvert,\lvert c-B\rvert\ge 3\) | **PASS** \(0.008^\circ\) | **CONFIRMED** \(4\times 10^{-6\circ}\) |

`PAIR_NOT_LINEAR`. *computed* (exact basis).
The mass maps add. The \(\mathbf g\)s do not,
except far away where the gradients are
parallel. To use two packets: form
\(M=M_A+M_B\), then one \(\hat n=\nabla_c M\).
Do not add the one-body fields.
Not Einstein.

## Equation-to-code

| Object | Where |
| --- | --- |
| Exact pair vs one-body \(\mathbf g\) | `scripts/m9_61_super.py` |
| Adversary \(N=11\) | `scripts/m9_61_audit_super.py` |

Paper: [`../latex/71_Superpose.tex`](../latex/71_Superpose.tex).
