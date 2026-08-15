# M9.6 method note: de Sitter decided at the FGHMV standard

> Equations first. The \pP{} is an obstruction, not a dual. Task:
> [`../tasks/m9_6_task_details.md`](../tasks/m9_6_task_details.md).
> Paper: [`../latex/17_deSitter_FGHMV_Obstruction.tex`](../latex/17_deSitter_FGHMV_Obstruction.tex).
> Jacobson is not used. No `MODELS.md` cell.

## VERDICT (one paragraph)

FGHMV-standard closure of the cosmological domain is **obstructed**.
Gibbons-Hawking gives \(S=3\pi/(G\Lambda)\) and
\(\partial S/\partial\Lambda=-3\pi/(G\Lambda^2)<0\) (*derived*,
auditor CONFIRMED identically in sympy). Casini-Huerta-Myers has a
positive \(T_{00}\) kernel, so FGHMV sends \(+\delta\langle T_{00}\rangle\)
to \(+\delta S\). The cosmological-horizon first law sends added
Killing energy to \(-\delta S\) (Banihashemi-Jacobson-Svesko-Visser
*cited*). Those signs conflict. Independently,
\(\dim\mathfrak{so}(1,4)=10<15=\dim\mathfrak{so}(2,4)\) (*computed*):
de Sitter isometries do not contain the CHM family of arbitrary
balls. Einstein+\(\Lambda\) from a cosmological CFT is **not**
labeled \pP{}. The AdS theorems stand.

## 1. Objects

\[
S_{\mathrm{GH}}=\frac{\pi\ell^2}{G}=\frac{3\pi}{G\Lambda},
\qquad
\Lambda=\frac{3}{\ell^2},
\qquad
\frac{\partial S_{\mathrm{GH}}}{\partial\Lambda}
=
-\frac{3\pi}{G\Lambda^2}.
\]

\[
H_B
=
2\pi\int_{\lvert x\rvert<R}
\frac{R^2-\lvert x\rvert^2}{2R}\,T_{00},
\qquad
\dim\mathfrak{so}(p,q)=\frac{(p+q)(p+q-1)}{2}.
\]

## 2. Equation-to-code map

| Object | Function | File |
| --- | --- | --- |
| \(S=\pi\ell^2/G\) | `gh_entropy` | `scripts/m9_6_ds_sign.py` |
| \(\dim\mathfrak{so}(p,q)\) | `so_dim` | same |
| sympy \(S(\Lambda)\), \(\partial S/\partial\Lambda\) | `main` | `scripts/m9_6_audit_ds.py` |
| combinatorial \(\dim\mathfrak{so}(n)\) | `so_dim_count` | same |

## 3. Results after methods

| ID | Result | Status |
| --- | --- | --- |
| C1 | \(S\) matches \(3\pi/(G\Lambda)\); \(\partial S/\partial\Lambda<0\) | PASS, auditor CONFIRMED (identity) |
| C2 | \(10\) and \(15\) | PASS, auditor CONFIRMED |
| C3 | CHM sign \(+1\) | PASS (cited formula) |
| C4 | C1 and C3 opposite | PASS |
| C5 | AdS branch \(\partial S/\partial\Lambda>0\) | PASS |

Verdict string: `FGHMV_STANDARD_DS_CLOSURE_OBSTRUCTED`.

## 4. What is \pP{}, and what is not

| Statement | Tag |
| --- | --- |
| \(S=3\pi/(G\Lambda)\), \(\partial S/\partial\Lambda<0\) | *derived* / *proved* |
| \(\dim\mathfrak{so}(1,4)=10\), \(\dim\mathfrak{so}(2,4)=15\) | *proved* (counting) |
| Copy of FGHMV onto the cosmological horizon | *proved false* (sign) |
| CHM net of balls from dS isometries | *proved false* (dimension) |
| Static-patch first-law minus sign | *cited* (JHEP 01 (2023) 054) |
| Einstein+\(\Lambda\) from a cosmological CFT | *unresolved* as a positive |
| Jacobson equilibrium as dS closure | *conjectured*; not used |
| A dS holographic dual | *not claimed* |

## 5. Adversarial audit

Auditor: `m9_6_audit_ds.py`. No solver import. Recovers
\(S=3\pi/(G\Lambda)\) and \(\partial S/\partial\Lambda=-3\pi/(G\Lambda^2)\)
as sympy identities. Counts antisymmetric generators of \(5\times5\)
and \(6\times6\) matrices.

| Claim | Auditor |
| --- | --- |
| C1, C2, C4 | CONFIRMED |
| Dual / Einstein+\(\Lambda\) from a CFT | NOT_CLAIMED |
