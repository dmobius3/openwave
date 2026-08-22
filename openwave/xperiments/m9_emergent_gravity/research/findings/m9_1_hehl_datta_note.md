# M9.1 method note: Hehl-Datta elimination

> Certification gate for the proposed Emergent Gravity / NSM column. Equations
> first. The reader should be able to audit the coefficient by reading this page
> and clicking the equation-to-code map. Task record:
> [`../tasks/m9_1_task_details.md`](../tasks/m9_1_task_details.md). Code links
> are `blob/main` paths for the intended landing; they resolve after the first
> PR merges.

## VERDICT (one paragraph)

The on-shell four-fermion Lagrangian produced by Palatini Einstein-Cartan plus
a Hermitian Dirac field satisfies

\[
\frac{\mathcal{L}_{\mathrm{int}}}{-\kappa\, J_5\cdot J_5} = \frac{3}{16}
\]

in both metric signatures, residual \(<10^{-15}\) (*computed*, 8 random spinors
each, seed 20260815). That is C1, the gate, and it **passes**. The paper's
printed spin identification \(s^{\lambda\mu\nu}=-\frac14\varepsilon^{\lambda\mu\nu\rho}J_{5\rho}\)
does **not** match the spin tensor defined by \(\delta S_D/\delta\omega\) in
this note: the measured factor is \(-\frac12\), not \(-\frac14\) (C2 **fails**
as pre-registered). The mismatch is a \(2\delta L/\delta\omega\) versus
\(\delta L/\delta\omega\) convention. C1 does not use that printed factor: the
contorsion is fixed by stationarity of the action, not by inserting
\(T=\kappa s\). Doubling the Palatini normalization moves the ratio to
\(3/32\); a fake source \(s=+\frac14(*J_5)\) moves it to \(3/64\). The check
can fail. This does **not** certify FGHMV, Condition NL, \(I_B\), or de Sitter.

## 1. Objects and conventions

```text
kappa = 8 pi G
eps_{0123} = +1
eta = diag(+1,-1,-1,-1)   or   diag(-1,+1,+1,+1)

e^a = coframe,  omega^{ab} = independent Lorentz connection
T^a = de^a + omega^a_b /\\ e^b
R^{ab} = d omega^{ab} + omega^a_c /\\ omega^{cb}

gamma_{ab} = (1/2) [gamma_a, gamma_b]
gamma^5 (mostly minus) = i gamma^0 gamma^1 gamma^2 gamma^3
J5^mu = psibar gamma^5 gamma^mu psi
```

Hermitian Dirac kinetic term:

\[
\mathcal{L}_D
=
\frac{i}{2}\,e\,\bar\psi\gamma^\mu\overleftrightarrow{D}_\mu\psi
-
e\,m\,\bar\psi\psi,
\qquad
D_\mu=\partial_\mu+\tfrac14\omega_\mu^{ab}\gamma_{ab}.
\]

The \(\omega\)-linear piece is

\[
\mathcal{L}_D[\omega]
=
\frac{i}{8}\,e\,\omega_\mu^{ab}\,
\bar\psi\{\gamma^\mu,\gamma_{ab}\}\psi.
\]

**Definition (this note).** The canonical spin used below is

\[
s^\mu{}_{ab}
:=
\frac{i}{4}\,\bar\psi\{\gamma^\mu,\gamma_{ab}\}\psi,
\]

so that \(\mathcal{L}_D[\omega]=\frac12 e\, s^\mu{}_{ab}\,\omega_\mu^{ab}\).
This is Hehl's \(\tau=2\delta L/\delta\omega\) normalization, not
\(\delta L/\delta\omega\).

Palatini:

\[
S_g
=
\frac{1}{4\kappa}
\int\varepsilon_{abcd}\,e^a\wedge e^b\wedge R^{cd}
=
\int\frac{e\,R}{2\kappa}.
\]

Decompose \(\omega=\mathring{\omega}+K\). The quadratic remainder of the
curvature scalar is

\[
Q
=
\eta^{\sigma\nu}
\bigl(
K^\mu{}_{\lambda\mu}K^\lambda{}_{\sigma\nu}
-
K^\mu{}_{\lambda\nu}K^\lambda{}_{\sigma\mu}
\bigr),
\qquad
\mathcal{L}_g[K]=\frac{Q}{2\kappa}.
\]

Totally antisymmetric contorsion is parameterized by a 4-vector,

\[
K_{abc}=\varepsilon_{abcr}\,v^r
\qquad(\varepsilon_{0123}=+1,\text{ component contraction}).
\]

**On-shell rule (this note).** \(v\) is the stationary point of
\(\mathcal{L}_g[K]+\mathcal{L}_D[K]\). The action is exactly quadratic in
\(K\), so the Hessian is constant and the critical point is unique.

**Reported observable.**

\[
r
:=
\frac{\mathcal{L}_{\mathrm{int}}(v_\star)}{-\kappa\,J_5^\mu J_{5\mu}}.
\]

The rational \(3/16\) is not an input to the contraction. It is compared after
extraction.

## 2. Equation-to-code map

Physics lives in
[`hehl_datta.py`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m9_emergent_gravity/research/scripts/hehl_datta.py).
The driver only samples spinors, compares, and mutates.

| Term | Function | Link |
| --- | --- | --- |
| Clifford + \(\gamma^5\) | `signature_mostly_minus`, `signature_mostly_plus`, `check_clifford` | [hehl_datta.py L73-L98](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m9_emergent_gravity/research/scripts/hehl_datta.py#L73-L98) |
| \(\bar\psi=\psi^\dagger A\) | `adjoint` | [L160-L166](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m9_emergent_gravity/research/scripts/hehl_datta.py#L160-L166) |
| \(J_5^\mu\) | `axial_current` | [L169-L175](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m9_emergent_gravity/research/scripts/hehl_datta.py#L169-L175) |
| \(s^\mu{}_{ab}\) | `spin_tensor_from_dirac` | [L178-L193](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m9_emergent_gravity/research/scripts/hehl_datta.py#L178-L193) |
| \((*J)^{\lambda\mu\nu}\) | `dual_of_vector` | [L206-L222](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m9_emergent_gravity/research/scripts/hehl_datta.py#L206-L222) |
| Palatini \(Q\) | `palatini_quadratic` | [L230-L248](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m9_emergent_gravity/research/scripts/hehl_datta.py#L230-L248) |
| \(\mathcal{L}(K)\) | `lagrangian_of_k` | [L271-L279](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m9_emergent_gravity/research/scripts/hehl_datta.py#L271-L279) |
| \(v_\star\) | `stationary_vector` | [L282-L318](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m9_emergent_gravity/research/scripts/hehl_datta.py#L282-L318) |
| ratio \(r\) | `contact_lagrangian` | [L353-L404](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m9_emergent_gravity/research/scripts/hehl_datta.py#L353-L404) |
| C5 mutations | `mutated_ratio` | [L407-L470](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m9_emergent_gravity/research/scripts/hehl_datta.py#L407-L470) |
| comparison / JSON | driver | [m9_1_hehl_datta_elimination.py](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m9_emergent_gravity/research/scripts/m9_1_hehl_datta_elimination.py) |

`palatini_quadratic_alt` is the same \(Q\) by explicit index
sums. A previous contraction swapped the first two indices of
the second factor and returned \(-Q\) for totally antisymmetric
\(K\) (residuals \(O(10)\)). That spectator is fixed. Residual
now \(7.1\times 10^{-15}\) (minus) and \(4.4\times 10^{-15}\)
(plus). It is not used in \(r\). The ratio \(3/16\) is
unchanged.

## 3. Pre-registered gates and results

Seed 20260815, \(\kappa=1\), 8 spinors per signature, tolerance \(10^{-12}\).
Raw table: [`../data/m9_1_hehl_datta.json`](../data/m9_1_hehl_datta.json).

| ID | Gate | Result | Status |
| --- | --- | --- | --- |
| C1 | \(r=3/16\) both signatures | minus \(0.1874999999999999\), plus \(0.1875\), sample scatter \(<10^{-15}\) | PASS (*computed*) |
| C2 | \(s=-\frac14(*J_5)\) | measured dual factor \(-1/2\) both signatures; residual vs paper \(0.25\,\lVert J_5\rVert\) | FAIL (*computed*) |
| C3 | \(s\) totally antisymmetric | mixed-symmetry residual \(<3\times10^{-16}\) | PASS (*computed*) |
| C4 | \(r\) signature-independent | \(\Delta r=1.1\times10^{-16}\) | PASS (*computed*) |
| C5 | mutations leave \(3/16\) | double Palatini \(\to 3/32=0.09375\); fake \(s=+\frac14(*J_5)\) \(\to 3/64=0.046875\) | PASS (*computed*) |

C1 is the column gate. It passes. C2 is a real discrepancy with Paper III's
printed factor and is left as a fail, not silently redefined.

## 4. What C2 failing does and does not mean

Under the definition \(s^\mu{}_{ab}=\frac{i}{4}\bar\psi\{\gamma^\mu,\gamma_{ab}\}\psi\),

\[
s^{\lambda\mu\nu}=-\tfrac12\,(*J_5)^{\lambda\mu\nu}.
\]

Paper III writes \(s=-\frac14\varepsilon J_5\). That is the factor one gets by
taking \(s=\delta L/\delta\omega\) instead of \(2\delta L/\delta\omega\), or by
absorbing a 2 into the Cartan map \(T=\kappa s\). Because this note **does not
insert** \(T=\kappa s\), and instead solves \(\delta(\mathcal{L}_g+\mathcal{L}_D)=0\),
C1 is insensitive to that labeling. A later holographic task that *does* use
the printed \(s=-\frac14\varepsilon J_5\) inside Condition NL must say which
normalization it inherited.

## 5. Scope honesty (not computed)

- Faulkner-Guica-Hartman-Myers-Van Raamsdonk equivalence.
- Condition NL, absorbability, the residue \(I_B\).
- Hehl-Datta *from relative entropy* (Papers V, VII).
- de Sitter / cosmology (Paper IX).
- Nonlinear Einstein-Cartan.
- Lattice Newton \(1/r^2\), light bending, \(\Lambda\).
- Collider or astrophysical bounds; the high-density bounce.
- Any `MODELS.md` icon.

## 6. Minimal inspection set

1. This note, sections 1-3.
2. [`hehl_datta.py`](../scripts/hehl_datta.py) (the functional).
3. [`m9_1_hehl_datta.json`](../data/m9_1_hehl_datta.json).
4. The auditor script and JSON, section 10 (last).

## 7. Regeneration

From the repository root:

```bash
python3 openwave/xperiments/m9_emergent_gravity/research/scripts/m9_1_hehl_datta_elimination.py
python3 openwave/xperiments/m9_emergent_gravity/research/scripts/m9_1_audit_hehl_datta.py
```

No GPU. Deterministic given the seed.

## 10. Independent adversarial audit

A second script,
[`m9_1_audit_hehl_datta.py`](../scripts/m9_1_audit_hehl_datta.py), was
forbidden to import the solver. It uses (i) the Levi-Civita identity
\(\varepsilon_{lmnr}\varepsilon^{lmns}=6\delta_r{}^s\) in both numpy and
sympy, (ii) a direct 4-form evaluation of Palatini on a random curvature,
(iii) complete-the-square on \(K=*v\), and (iv) a **Weyl / chiral** Clifford
in sympy to measure \(\alpha\) in \(s=\alpha(*J_5)\). Raw table:
[`../data/m9_1_audit.json`](../data/m9_1_audit.json).

| Check | Auditor number | Verdict |
| --- | --- | --- |
| \(\varepsilon_{lmn r}\varepsilon^{lmn s}\) | 6 (numpy and sympy) | identity holds |
| Palatini 4-form / \(R\) | \(0.5000000000000001\) | CONFIRMED \(eR/2\) |
| \(Q/v^2\) on \(K=*v\) | \(6\) | used in the square |
| Weyl \(\alpha\) | \(-1/2\) | C2 **REFUTED** |
| \(r=9\alpha^2/(2\,Q/v^2)\) with Weyl \(\alpha\) | \(0.18750000000000006\) | C1 **CONFIRMED** |
| same formula with paper \(\alpha=-1/4\) | \(0.046875=3/64\) | would have killed C1; it is the wrong \(\alpha\) for this \(s\) |
| doubled \(Q/v^2\) | \(0.09375=3/32\) | C5 **CONFIRMED** (the check can fail) |
| C3 mixed residual | \(8.9\times10^{-16}\) | CONFIRMED |
| C4 | not re-run on a second Clifford | QUALIFIED |

C1 survives a method that never formed a finite-difference Hessian and never
imported the solver. C2 is independently refuted: two representations (Dirac
and Weyl) give \(\alpha=-1/2\). Using the paper's printed \(-1/4\) *inside
this \(s\) definition* produces \(3/64\), which is exactly the solver's C5
fake-source mutation. That is a convention trap, not a C1 failure.

A separate subagent was also launched with a refute-only brief and no access
to the solver. If that agent files a conflicting JSON, the conflict is itself
a finding and this section must be updated before any outward send.

