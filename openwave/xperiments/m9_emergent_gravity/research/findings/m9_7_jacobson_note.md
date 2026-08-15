# M9.7 method note: Jacobson is not a [P] substitute

> [P] only. Paper:
> [`../latex/18_Jacobson_Not_P_Substitute.tex`](../latex/18_Jacobson_Not_P_Substitute.tex).
> Task: [`../tasks/m9_7_task_details.md`](../tasks/m9_7_task_details.md).

## VERDICT (one paragraph)

Jacobson 1995, under the hypotheses he states, implies Einstein
with free \(\Lambda\) and without torsion (*derived* from the null
lemma + Bianchi; *cited* PRL 75, 1260). The null lemma is an
identity: auditor kernel \((-1,1,1,1,0,0,0)\) after normalizing
\(S_{00}=-1\). Jacobson 2016's conformal half requires a CFT;
\(b_3=-7\) so it does not apply to NSM matter (*computed*). The
2016 nonconformal half is a conjecture in Jacobson's own sentence
and is not labeled [P]. Paper IX requires Einstein-Cartan plus
Hehl-Datta. Therefore Jacobson is **not** a [P] substitute for Q2.

## 1. Objects

\[
T_{ab}k^ak^b=\frac{\hbar\eta}{2\pi}R_{ab}k^ak^b
\quad(\forall\text{ null }k)
\;\Longrightarrow\;
G_{ab}+\Lambda g_{ab}=\frac{2\pi}{\hbar\eta}T_{ab},
\quad
G=\frac{1}{4\hbar\eta}.
\]

\(\Lambda\) is not fixed. No \(J_5\).

## 2. Equation-to-code map

| Object | Function | File |
| --- | --- | --- |
| \(S=f\eta\) has \(S_{kk}=0\) | `c1_null_lemma` | `scripts/m9_7_jacobson.py` |
| Polarization kernel of \(S_{kk}=0\) | `null_lemma_sympy` | `scripts/m9_7_audit_jacobson.py` |
| SM \(b_3\) | `sm_not_cft` / `sm_b3` | both |

## 3. Results

| ID | Result | Status |
| --- | --- | --- |
| C1 | kernel dim 1, spanned by \(\eta\) | PASS / CONFIRMED |
| C2 | \(\Lambda\) free | PASS |
| C3 | no spin, no \(J_5\) | PASS |
| C4 | conservation is load-bearing | PASS |
| C5 | \(b_3=-7\) | PASS / CONFIRMED |

Verdict string: `JACOBSON_NOT_A_P_SUBSTITUTE`.

## 4. Epistemic tags ([P] only)

| Statement | Tag |
| --- | --- |
| Null lemma | *proved* |
| 1995 \(\Rightarrow\) Einstein, \(\Lambda\) free | *cited* + *derived* (algebra) |
| 1995 \(\Rightarrow\) HD or EC | *proved false* |
| 2016 conformal half applies to SM | *proved false* |
| 2016 nonconformal | not [P] (Jacobson's word: conjecture) |
| Jacobson substitutes for Q2 at the Paper IX bar | *proved false* |
