# M9.4 method note: axial UV deformation, not a selected completion

> Equations first. SuperGrok does not repeal Paper 11: inventing a
> CFT to fill Q4 is a fabrication. This note constructs the unique
> *quadratic* deformation forced by Paper 14 Thm. obs and checks it.
> Task: [`../tasks/m9_4_task_details.md`](../tasks/m9_4_task_details.md).
> Paper: [`../latex/15_UV_Deformation_Axial.tex`](../latex/15_UV_Deformation_Axial.tex).
> No `MODELS.md` cell. Final Status Q4a stays `[O]`.

## VERDICT (one paragraph)

Q4 splits. **Q4a** (a selected holographic pair, a string
compactification, or a metric UV completion that produces the
installed SM) is still open by construction and is **not invented
here**. **Q4b** is the unique local quadratic axial action that
(i) recovers the M9.1 ratio \(3/16\) as \(M\to\infty\) and
(ii) is nonlocal at finite \(M\). Completing the square and locking
the infrared coefficient gives

\[
r(k)=\frac{3/16}{1+k^2/M^2}.
\]

That algebra is *derived*. The solver and a no-import auditor both
report C1-C5 PASS (*computed*). The position-space kernel is the
4d Yukawa function; the auditor's radial Laplacian of
\((M/(4\pi^2 r))K_1(Mr)\) is identically zero. This evades the
contact obstruction by *changing the theory*. It is not Einstein-Cartan
at finite \(M\), and it is not a selected ultraviolet completion.

## 1. Objects

Infrared NSM axial sector (M9.1, Hehl-Datta 1971):

\[
\mathcal{L}_{\mathrm{HD}}
=
-\frac{3\kappa}{16}\,J_5\cdot J_5,
\qquad
\kappa=8\pi G,
\qquad
r_{\mathrm{IR}}=3/16.
\]

Ultraviolet deformation, axial channel only, tree level:

\[
\mathcal{L}_k
=
A\bigl(1+k^2/M^2\bigr)\lvert S(k)\rvert^2
+
B\,S(-k)\cdot J_5(k).
\]

On-shell,

\[
\mathcal{L}_{\mathrm{eff}}(k)
=
-\frac{B^2}{4A(1+k^2/M^2)}\,J_5\cdot J_5.
\]

Infrared lock \(B^2/(4A)=3\kappa/16\) is imposed so that M9.1 is
not given away. Without the lock the infrared coefficient is
\(B^2/(4A\kappa)\), a free parameter (auditor, symbolic).

Euclidean 4d kernel of \(1/(k^2+M^2)\):

\[
G(x)=\frac{M}{4\pi^2 r}\,K_1(Mr),
\qquad
(-\square+M^2)G=\delta^{(4)}.
\]

## 2. Equation-to-code map

Task-scoped files; `blob/main` after merge.

| Object | Function | File |
| --- | --- | --- |
| \(r(k)=(3/16)/(1+k^2/M^2)\) | `r_of_k2` | `scripts/m9_4_uv_axial.py` |
| \(G=(M/(4\pi^2 r))K_1(Mr)\) | `yukawa_4d` | same |
| Finite-difference 4d radial ODE (solver check) | `main` C3 | same |
| Symbolic complete-the-square | `complete_the_square` | `scripts/m9_4_audit_uv.py` |
| Symbolic \((-\square+M^2)G=0\) for \(r>0\) | `yukawa_identity_residual` | same |

## 3. Results after methods

Tolerance: \(10^{-12}\) on dimensionless ratios unless noted.

| ID | Result | Status |
| --- | --- | --- |
| C1 | \(r(0)=0.1875=3/16\); \(M=10^{12}\) recovers \(3/16\) at all sampled \(k\) | PASS, auditor CONFIRMED |
| C2 | \(r(M)=0.09375=3/32\); \(r\) strictly decreasing in \(k^2\) | PASS, auditor CONFIRMED |
| C3 | solver 4d ODE relative residual \(8.76\times 10^{-7}\); auditor symbolic ODE \(=0\) | PASS / CONFIRMED |
| C4 | \(M\to 0\) at fixed \(k\neq 0\) gives \(r=0\) (sympy limit) | PASS, auditor CONFIRMED |
| C5 | unlocked infrared target \(3/8\) stays \(3/8\); lock is necessary | PASS |

Verdict string: `IR_MATCHING_HOLDS_DEFORMATION_NOT_SELECTION`.

## 4. What this is not

- A CFT. A compactification. Asymptotic safety.
- A metric ultraviolet completion. Einstein-Hilbert is still
  non-renormalizable.
- A proof that \(S_\mu\) *is* axial torsion. That reading is
  `[C]`. The algebra does not care.
- A laboratory torsion wave. Compatibility with the program
  requires \(M\gtrsim M_{\mathrm{Pl}}\) (consistency, not a
  derivation).
- A `MODELS.md` cell. A move of Final Status "UV completion" to
  `[P]`.

## 5. Adversarial audit

Auditor: `m9_4_audit_uv.py`. No solver import. Completes the square
in sympy; checks the Yukawa ODE by differentiating
`besselk(1, M r)` symbolically. First FD attempt on \(G''\) was
discarded: it reported a \(6\times 10^{-6}\) residual from step
size, which is not a physics fail. The symbolic Laplacian is
identically zero.

| Claim | Auditor |
| --- | --- |
| C1-C4 | CONFIRMED |
| Q4a selected UV | NOT_CLAIMED |

## 6. Epistemic tags

| Statement | Tag |
| --- | --- |
| Unique quadratic \(r(k)\) given the infrared lock | *derived* |
| C1-C5 numbers | *computed* (tolerance as above) |
| \(S_\mu\) is axial torsion | *conjectured* (identification) |
| \(M\gtrsim M_{\mathrm{Pl}}\) | *consistency requirement*, not derived |
| Selected UV of the NSM | *unresolved* (Q4a) |
| Sezgin-van Nieuwenhuizen ghost-free PGT | *cited* (Phys. Rev. D 21 (1980) 3269), not re-derived |
| Witten dictionary | *cited* (ATMP 2 (1998) 253), not re-derived |
