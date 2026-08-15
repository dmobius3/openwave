# M9.16: first law with the CHM kernel (effective gates)

> Designed to fail if the CHM envelope is not doing Clausius work.
> It failed that gate. A local kernel is enough; the envelope is not
> selected.

## What was discarded

Bond-local \(\chi=\Delta S/\Delta t\) vs \(w\) has 1d \(\rho\approx 0.28\).
\(\Delta C\) is nonlocal. That instrument is rejected.

An independent adversary, own code, never saw the locked
\(\mathrm{Tr}(K_{\mathrm{CHM}}\Delta C)\) test, attacked \(\chi\)
anyway and **REFUTED** it: 1d \(\rho\le 0.42\) on four lattices;
3d \(N=13\) \(\rho=0.21\); on the clean even grid \(\lvert\chi\rvert\)
is *largest at the cut* (surface, not bulk-peaked CHM). Exterior
\(\lvert\chi\rvert\) exceeds interior. That is area-law surface
response, the opposite Clausius weight. It confirms the discard.
It does not audit Paper 26's locked observable.

A permutation of \(w\) is not a mutation: the intercept of
\(K_{ij}=aw_{ij}+b\) dominates, so every NN kernel looks the same.

## Equations

Vacuum Peschel kernel \(K=\log((1-C)/C)\). Fit a CHM hop kernel

\[
(K_{\mathrm{CHM}})_{ij}=a\,w_{ij}+b,\qquad
w_{ij}=\tfrac12(w_i+w_j),\qquad
w=R^2-r^2
\]

on nearest neighbours only. Null kernel: \((K_{\mathrm{flat}})_{ij}=\langle K_{\mathrm{NN}}\rangle\)
on every hop. Perturb one hop, \(t\to t(1+\varepsilon)\), and compare

\[
\delta S=S_\varepsilon-S_0
\qquad\text{to}\qquad
\delta S_K=\mathrm{Tr}\bigl(K(C_\varepsilon-C_0)\bigr)
\]

for \(K=K_{\mathrm{CHM}}\) and \(K=K_{\mathrm{flat}}\).
\(\mathrm{Tr}(K_{\mathrm{exact}}\delta C)\) is an identity and is
not used as a gate.

## Pre-registered gates

| ID | Claim | Result |
| --- | --- | --- |
| C0 | 1d \(\rho(\delta S,\delta S_{\mathrm{CHM}})>0.80\) | PASS, \(0.918\) |
| C1 | 3d all NN, \(\rho>0.60\) | PASS, \(0.850\) |
| C2 | \(R_{\mathrm{shape}}<0.70\) | PASS, \(0.527\) |
| C3 | \(R_{\mathrm{shape}}(\mathrm{CHM})<R_{\mathrm{shape}}(\mathrm{flat})\) | **FAIL** |

## Verdicts

| Object | CHM | flat | Tag |
| --- | --- | --- | --- |
| 1d \(\rho\) / \(R_{\mathrm{shape}}\) | \(0.918\) / \(0.397\) | \(0.981\) / \(0.192\) | *computed* |
| 3d \(282\) hops \(\rho\) / \(R_{\mathrm{shape}}\) | \(0.850\) / \(0.527\) | \(0.912\) / \(0.411\) | *computed* |
| auditor 1d \(L=20\) \(\rho\) | \(0.775\) (below C0 \(0.80\)) | --- | QUALIFIED |
| auditor 3d \(N=10\) | \(0.797\) / \(0.604\) | \(0.865\) / \(0.501\) | C3 CONFIRMED |

`FIRST_LAW_FAIL` on C3. The first law **tracks a local hop kernel**.
The CHM *envelope is not selected*: a shapeless local \(K\) predicts
\(\delta S\) better. Paper 25 remains a statement about the shape of
\(K\), not about Clausius.

That is the effective result. Locality (A2) is the first-law input.
The geometric weight is extra structure that this probe does not
need.

## Equation-to-code

| Object | Where |
| --- | --- |
| \(K_{\mathrm{CHM}}\), \(K_{\mathrm{flat}}\), \(\delta S\) | `scripts/m9_16_first_law.py` |
| Independent sizes | `scripts/m9_16_audit_first_law.py` |

## What this is not

Not a refutation of CHM as a continuum theorem. Not a refutation of
Paper 25 ( \(K_{\mathrm{NN}}\) still tracks \(w\) ). Not
\(\eta=1/4G\). Not Einstein. Not foam. Not de Sitter.

Scripts: `m9_16_first_law.py`, `m9_16_audit_first_law.py`.
Data: `data/m9_16_first_law.json`, `data/m9_16_audit_first_law.json`.
Paper: [`../latex/26_First_Law_Local_Not_CHM.tex`](../latex/26_First_Law_Local_Not_CHM.tex).
