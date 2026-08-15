# M9.5 method note: Q4a as holographic-pair selection

> Equations first. Q4a in this program is not "invent a CFT". It is
> the monograph's question: the NSM is bookkeeping for a holographic
> pair; is that pair selected by the certified principles? Task:
> [`../tasks/m9_5_task_details.md`](../tasks/m9_5_task_details.md).
> Paper: [`../latex/16_Q4a_Holographic_Pair.tex`](../latex/16_Q4a_Holographic_Pair.tex).
> No `MODELS.md` cell. Final Status "UV completion" stays `[O]`
> (existence). Selection-uniqueness is answered in the negative.

## VERDICT (one paragraph)

Two readings. **Reading B** (the Standard Model *is* the boundary
CFT) fails necessary conditions: one-loop
\((b_1,b_2,b_3)=(41/10,-19/6,-7)\), all nonzero (*computed*,
auditor CONFIRMED from a different species sum), and \(118\) light
on-shell degrees of freedom before EWSB, not a large-\(N\) gap
(*computed*; Heemskerk-Penedones-Polchinski-Sully *cited*).
**Reading A** (SM fields are bulk; the CFT is something else) is
the monograph. The certified first law uses only \(T_{\mu\nu}\).
\(G_{\mathrm{SM}}\), \(n_g\), and the Yukawas can be replaced
without leaving the class of pairs the principles allow. The pair
is underdetermined. Existence of some pair remains `[O]`. No CFT
is named.

## 1. Objects

Casini-Huerta-Myers ball modular Hamiltonian:

\[
H_B
=
2\pi\int_{\lvert x\rvert<R}
\frac{R^2-\lvert x\rvert^2}{2R}\,T_{00}.
\]

One-loop convention \(\beta_i=b_i g_i^3/(16\pi^2)\), GUT-normalized
\(U(1)\), \(Q=T_3+Y/2\):

\[
b
=
-\frac{11}{3}C_2(G)
+\frac{2}{3}T_F
+\frac{1}{3}T_S.
\]

NSM matter ledger (installed, not derived):
\(G_{\mathrm{SM}}\), \(n_g=3\), Yukawas, \(\theta_{\mathrm{QCD}}\),
Higgs potential, \(\Lambda\).

## 2. Equation-to-code map

| Object | Function | File |
| --- | --- | --- |
| SM \(b_i\) from Dynkin sums | `sm_b_coefficients` | `scripts/m9_5_q4a_pair.py` |
| On-shell dof before EWSB | `sm_on_shell_dof` | same |
| \(\mathcal{N}=4\) one-loop \(b\) | `n4_sym_b_su_n` | same |
| Species-table recompute | `b_from_species` | `scripts/m9_5_audit_q4a.py` |

## 3. Results after methods

| ID | Result | Status |
| --- | --- | --- |
| C1 | \(b_1=41/10\), \(b_2=-19/6\), \(b_3=-7\), all nonzero | PASS, auditor CONFIRMED |
| C2 | dof \(=24+4+90=118<200\) | PASS (count) |
| C3 | CHM / FGHMV lists contain \(T_{\mu\nu}\), \(C_T\), \(R\); not flavor data | PASS (structural) |
| C4 | Paper II scope is Young symmetry, not \(G_{\mathrm{SM}}\) | PASS (structural) |
| C5 | \(\mathcal{N}=4\) \(b=0\) | PASS, auditor CONFIRMED |

Verdict string: `Q4A_SELECTION_ANSWERED_NEGATIVE_EXISTENCE_OPEN`.

C3 and C4 are not numerical. They are readings of cited formulas.
The auditor does not pretend to re-derive CHM.

## 4. What this is not

- A holographic pair. A compactification. A selected CFT.
- A claim that \(\mathcal{N}=4\) is the UV. The mutation only shows
  the CFT-necessary check can fail.
- A metric ultraviolet completion. That is still Q4a existence, `[O]`.
- A `MODELS.md` cell. A move of Final Status "UV completion" to `[P]`.

## 5. Adversarial audit

Auditor: `m9_5_audit_q4a.py`. No solver import. Sums hypercharge and
color indices generation-by-generation from an explicit species
table. Recovers the same three fractions. Recovers \(\mathcal{N}=4\)
\(b=0\) from \(T(\mathrm{adj})=C_2\).

| Claim | Auditor |
| --- | --- |
| C1, C5 | CONFIRMED |
| C3 | CONFIRMED as structural |
| Pair constructed | NOT_CLAIMED |
| Existence | OPEN |

## 6. Epistemic tags

| Statement | Tag |
| --- | --- |
| SM \((b_1,b_2,b_3)\) | *computed* (exact fractions) |
| \(118\) light dof | *computed* (species count) |
| SM is not a CFT | *derived* from the \(b_i\) |
| SM is not a holographic CFT | *cited* (HPPS 2009) plus the count |
| CHM / FGHMV independent of flavor | *derived* from the cited formulas |
| Pair uniqueness | *derived* (underdetermination) |
| Existence of some pair | *unresolved* |
| Maldacena / Witten / CHM / FGHMV / HPPS | *cited* (real identifiers) |
