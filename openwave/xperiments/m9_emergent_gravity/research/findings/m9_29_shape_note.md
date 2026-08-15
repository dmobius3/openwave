# M9.29: not a sphere --- a guess, then a measurement of this probe

> The author guessed that the first-law kernel needs a
> *shape*, not a sphere. That stays a guess about gravity.
> This note only asks whether Paper 37's CHM win is
> ball-only on the same fixed-\(H\) occupation transfer.

## Equations

Same state as Paper 37: \(H\) fixed,

\[
C(\alpha)=C_0+\alpha\bigl(\lvert R\rangle\langle R\rvert
-\lvert L\rangle\langle L\rvert\bigr).
\]

Regions at every legal center: Euclidean balls \(R=2\)
(\(33\) sites), cubes of side \(3\) (\(27\)), taxicab
\(t=2\) (\(25\)). Predictors \(P=\sum_A w\,\delta e\) with

\[
w_{\mathrm{flat}}=1,\qquad
w_{\mathrm{export}}=r_{\max}^2-r^2
\]

(\(w_{\mathrm{export}}\) *is* CHM on a ball; off a ball it
is an illegal export). Shape-native, itself a guess:

\[
w_{\mathrm{cube}}=\prod_{\mu=x,y,z}\Bigl(\bigl(\tfrac{s}{2}\bigr)^2-d_\mu^2\Bigr),
\qquad
w_{\mathrm{taxi}}=t^2-\ell_1^2.
\]

Instrument: \(\delta S\stackrel{?}{=}\mathrm{Tr}(K_{\mathrm{vac}}\Delta C)\).

## Verdicts

Control: balls reproduce Paper 37. \(C_{\mathrm{vac}}=0.999\).
\(R_{\mathrm{CHM}}=0.018<R_{\mathrm{flat}}=0.102\).

| Family | \(C_{\mathrm{vac}}\) | \(R_{\mathrm{export}}\) | \(R_{\mathrm{native}}\) | \(R_{\mathrm{flat}}\) | C2e | C2n | C2x |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ball \(512\) | PASS | \(0.018\) | \(0.018\) | \(0.102\) | PASS | PASS | \(=\) |
| cube \(1000\) | PASS | \(0.047\) | \(0.022\) | \(0.063\) | PASS | PASS | PASS |
| taxi \(512\) | PASS | \(0.049\) | \(0.045\) | \(0.138\) | PASS | PASS | PASS |

Auditor \(N=10\): ball control **CONFIRMED**. Cube export
C2e **REFUTED** (\(0.087>0.077\)). Cube native C2n and C2x
**CONFIRMED**.

`SHAPE_NATIVE_MEASURED_STILL_A_GUESS`. *computed* for this
probe. The exported ball kernel is not robust on cubes.
A cube product weight fits cubes better than CHM or flat.
That does not derive a cube modular Hamiltonian, and it
does not make the author's guess a theorem of gravity.

Not \(8\pi G\). Not de Sitter. Not ``CHM off balls''.

## Equation-to-code

| Object | Where |
| --- | --- |
| Families, kernels, gates | `scripts/m9_29_shape.py` |
| Adversary, \(N=10\) | `scripts/m9_29_audit_shape.py` |

Paper: [`../latex/38_Shape_Not_Sphere.tex`](../latex/38_Shape_Not_Sphere.tex).
