# M5.32 method note: the autonomous Lagrangian hunt (rungs R0 to R10)

> **Status: the task is PAUSED, not finished.** This note is built to the
> [`METHOD_NOTE.md`](../../../../../dev_docs/METHOD_NOTE.md) standard: the reader must be able to
> audit every number below by reading, without trusting the run and without reverse-engineering
> Python. It is written for the model author and for the maintainer re-reading it later, who are the
> same reader. Section 8 carries the open questions.
>
> Task record: [`tasks/m5_32_task_details.md`](../tasks/m5_32_task_details.md) (plan, RUNG LOG, the
> pause records). Machine ledger: [`data/m5_32_ledger.json`](../data/m5_32_ledger.json).
> Code links resolve once the task's branch is merged to `main`.

## 1. The physics, before any result

### 1.1 Field, metric, conventions

```text
M(x)            a real 4x4 matrix field on a periodic cubic lattice
eta             diag(-1, +1, +1, +1); index 0 is time, as a derivative index AND as
                the internal row of M
A_mu            d_mu M, the jets, with RAW CONTRAVARIANT internal entries
Lorentz action  M -> L M L^T   (raw entries contravariant)
                under M -> L^-T M L^-1 the roles swap; the two agree on M_cov = eta M eta,
                so any covariant-metric object must be converted before mixing
vacuum          M_vac = diag(-s g, 1, delta, 0);  toy point s = -1, g = 32, delta = 0.3
lattice         h = L / n, certified symmetric stencil: the density is formed per
                forward and backward branch, then weight-averaged
```

Contraction rule, locked at R0 and audited: a derivative-derivative index pair contracts with
`eta`, an internal-internal pair with `eta`, and a MIXED derivative-internal pair with `delta`.
The all-`eta` reading is not covariant (measured boost drift 32.2, the R0 audit's figure) and is retained only as the
control term `I3_mixed_eta`.

### 1.2 The certified action

```text
F_mu nu      = A_mu eta A_nu  -  A_nu eta A_mu            (curvature, quadratic in the jets)
<F, G>_eta   = tr( eta F eta G^T )                        (the bracket)
I1           = sum_{mu < nu} eta^mu mu eta^nu nu <F_mu nu, F_mu nu>_eta
             = (1/2) F_abcd F^abcd
V4           = w sum_{p = 1..4} ( tr((M eta)^p) - C_p )^2 ,  C_p = (s g)^p + 1 + delta^p
w            = 7.24023879e-4
L_cert       = -4 I1  -  V4                               (CERTIFIED_COEFFS)
```

### 1.3 The clock, and the two energies

The clock is a one-parameter internal rotation applied in the body frame of the hedgehog ansatz.
Its tangent at `t = 0` is the field `a0`, and the time jet is `A_0 = omega a0`:

```text
a0(x)        = Qh(x) ( G1 d4 + d4 G1^T ) Qh(x)^T          (the co-moving flow)
G1           the (2,3)-plane rotation generator: G1[2,3] = -1, G1[3,2] = +1
Qh(x)        = R3(phi) R2(theta), the Euler frame carrying the eigenvalue-1
             eigenvector to n-hat = x / |x|
```

Every term in the registry is a polynomial in `omega`, so the Lagrangian read and the Hamiltonian
(energy) read are related term by term by a Legendre transform:

```text
quadratic terms    I(omega) = A + B omega + C omega^2
                   H_I      = C omega^2 - A
quartic terms      I(omega) = A + C2 omega^2 + C4 omega^4
                   J        = dI/domega = 2 C2 omega + 4 C4 omega^3
                   H_I      = C2 omega^2 + 3 C4 omega^4 - A
lattice energy     E_cert   = 4 (U + omega^2 T) + V4      (the factor 4 is |CERTIFIED_COEFFS[I1]|)
clock inertia      kin      = -4 x (the omega^2 coefficient of I1) = INS4.kin_of(M, a0, cfg)
                            the two measures agree on the rigid ansatz to ten digits
fixed-J energy     E_J      = E_stat + J^2 / (4 kin)  ,  omega* = J / (2 kin)
```

`B` (the `omega`-odd piece) is zero for `I1` and nonzero for the mixed contractions `I2` to `I6`;
it shifts the fixed-J relation to `omega* = (J - B) / (2 C)` but leaves `H` even, so it never
creates a free minimum. That Legendre argument is re-verified symbolically per term.

### 1.4 The candidate families tested

```text
lambda-family (class C2, rung R2)
    L_lambda = -4 [ (1 - lambda) I1 + lambda I1_h ] - V4
    h_cov    = eta + 2 (eta u)(eta u)^T ,  u the timelike unit eigenvector of M eta, u^T eta u = -1
    I1_h     the bracket with eta -> h_cov on the INTERNAL pair only

K_T (class C4, rung R7)
    K_T      = (1/2) sum_mu eta^mu mu [ tr(h A_mu h A_mu) - tr(eta A_mu eta A_mu) ]
             = 2 sum_mu eta^mu mu sum_j (A_mu)_{0j}^2      in the u-frame
    L        = -4 [ (1 - lambda) I1 + lambda I1_h ] - c2 K_T - V4 ,  c2 > 0

quartics (classes C5 and C6, rung R8)
    Q_I1sq   = (I1 density)^2
    Q_I4sq   = (I4 density)^2 ,  I4 = R_ac R^ac ,  R[nu, a] = sum_mu F[mu, nu, a, mu]
    Q_Fpair  = sum_{mu<nu, rho<sigma} eta-weighted <F_mu nu, F_rho sigma>_eta^2
    Q_C6a    = [ sum_mu eta^mu mu tr(A_mu eta A_mu eta) ]^2
    Q_C6b    = sum_{mu nu} eta^mu mu eta^nu nu [ tr(A_mu eta A_nu eta) ]^2
    Q_BI     = b^2 ( sqrt(1 + 2 I1 / b^2) - 1 ) ,  b^2 = 1e4      (not polynomial in omega)
```

### 1.5 The ansatz, the relaxation, and the degree

```text
ansatz       M = Q d4 Q^T ,  Q = Qb Qh ,  d4 = diag(-s g, 1, delta, 0)
             Qb = I + sinh(m) K + (cosh(m) - 1) K2, the boost dressing built from n-hat
             at m = 0 (used throughout R7 to R10) Q = Qh, the Euler frame alone
relaxation   FIRE on E_static only (a0 = None, omega = 0), boundary shell pinned at the ansatz
             by ~pin_shell(n, h) with default depth 1.6
degree       read_charge_from_M on the SPATIAL 3x3 block of M (every caller passes
             M[..., 1:4, 1:4]): eigh, take V[..., -1] (the LEADING eigenvector, n-hat), lift
             its sign field over a centered cube surface, integrate the RP^2 degree
```

The last line is stated here in full because rung R10 turns on it: the instrument reads ONE
eigenvector of the spatial block, not the full order parameter, and the time row is invisible to
it. The sign of that degree is a lift convention, so it is reported as `+-1` throughout.

### 1.6 The gate and class vocabulary used below

The run pre-registered seven gates a candidate had to pass, and worked through term classes in a
fixed order. Both are defined in full in the task record
([`tasks/m5_32_task_details.md`](../tasks/m5_32_task_details.md), sections 4 and 8); the short
form, so the tables below can be read without it:

| Gate | Statement |
| --- | --- |
| G1 | Coulomb kept: like charges repel with a 1/r trend, the static sector unchanged |
| G2 | Newton reversed: two boost-dressed defects attract |
| G3 | Electron clock: a finite nonzero omega* at positive energy, the vacuum preferring omega = 0 |
| G4 | Lorentz covariance of the energy functional, verified numerically |
| G5 | Bounded below with no guard on every runaway family |
| G6 | Collateral: the certified positives survive (census ordering, protection, the fixed-J clock) |
| G7 | Parsimony and robustness: few terms, O(1) coefficients, every sign holding over a factor >= 4 |

| Class | Content |
| --- | --- |
| C0, C1 | the author's two contractions, and the full quadratic basis of `F x F` |
| C2 | field-dependent internal metrics (`h_cov`), `M`-inserted contractions, projector currents |
| C3 | the potential axis: eigenvalue penalties, the LdG lift, dressing-sensitive `V` |
| C4 | the lower-order 2-derivative term (`K_T`) |
| C5, C6 | curvature^4 and saturation; non-commutator quartics |
| C7, C8 | higher-order timelike-current / Skyrme contractions; cross-model imports (never opened) |

## 2. Equation-to-code map

Base: `https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/`

| Object in section 1 | Function | File and lines |
| --- | --- | --- |
| `F_mu nu` from the jets | `F_of_A` | [`scripts/m5_32_lagrangian.py#L134`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py#L134-L141) |
| the bracket and every contraction pattern | `_K_from_pattern`, `density_from_K` | [`#L142`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py#L142-L182) |
| `I1` (sympy reference) | `I1_sym` | [`#L313`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py#L313-L321) |
| `I4 = R_ac R^ac`, the mixed trace | `I4_sym`, `R_readings_np` | [`#L355`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py#L355-L360), [`#L205`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py#L205-L219) |
| `V4` and its weight `w` | `V4_sym`, `v4_density_np`, `W1` | [`#L371`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py#L371-L382), [`#L236`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py#L236-L242), [`#L109`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py#L109) |
| `L_cert = -4 I1 - V4` | `CERTIFIED_COEFFS` | [`#L469`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py#L469) |
| `A_0 = omega a0`, the stencil branches | `lattice_jets` | [`#L479`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py#L479-L491) |
| `H_I = C omega^2 - A` | `term_hamiltonian` | [`#L500`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py#L500-L510) |
| `(A, B, C)` per term | `omega_decompose` | [`#L515`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py#L515-L522) |
| `K_T` density, both readings | `kt_density_np`, `kt_density_sym` | [`scripts/m5_32_r7_a_kt_form.py#L122`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r7_a_kt_form.py#L122-L154) |
| the u-frame time row | `uframe_time_row` | [`#L410`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r7_a_kt_form.py#L410-L419) |
| `E_stat` and `kin` under `L_lambda + c2 K_T` | `es_kin` | [`scripts/m5_32_r7_b_kt_lattice.py#L207`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r7_b_kt_lattice.py#L207-L213) |
| `E_J` minimized over the dressing | `min_over_amp`, `scan_R` | [`#L230`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r7_b_kt_lattice.py#L230-L275) |
| the six quartic densities | `d_I1`, `d_I4`, `d_Fpair`, `d_C6a`, `d_C6b`, `d_BI`, `QUARTICS` | [`scripts/m5_32_r8_a_quartics.py#L68`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r8_a_quartics.py#L68-L140) |
| the exact degree-4 `omega` extraction | `omega_poly` | [`#L158`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r8_a_quartics.py#L158-L175) |
| the generator enumeration `[X, M_vac]` | `generators`, `stage_generators` | [`scripts/m5_32_r8_b_ir_theorem.py#L53`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r8_b_ir_theorem.py#L53-L94) |
| the far-field tail measurement | `stage_tail` | [`#L95`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r8_b_ir_theorem.py#L95-L136) |
| the frame-free identity at `delta = 0` | `stage_equivalence` | [`scripts/m5_32_r9_b_string.py#L74`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r9_b_string.py#L74-L109) |
| the continuum ring (the string, measured) | `M_continuum`, `stage_continuum_ring` | [`#L110`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r9_b_string.py#L110-L147) |
| the fixed-physical-radius excision | `rho_of`, `run_box` | [`scripts/m5_32_r9_a_tube.py#L48`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r9_a_tube.py#L48-L82) |
| `kin` on a relaxed field, and its shells | `kin_c2`, `kin_shells` | [`scripts/m5_32_r10_relaxed_ladder.py#L72`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r10_relaxed_ladder.py#L72-L91) |
| the relaxation protocol | `relax` | [`#L92`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r10_relaxed_ladder.py#L92-L125) |
| FIRE, the pinned shell, `kin_of`, `e_parts` | `fire`, `pin_shell`, `kin_of`, `e_parts` | [`scripts/m5_21_3_a_4d.py#L327`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_21_3_a_4d.py#L327), [`#L109`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_21_3_a_4d.py#L109), [`#L274`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_21_3_a_4d.py#L274), [`#L179`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_21_3_a_4d.py#L179) |
| the ansatz and the clock tangent | `dressed`, `a0_unit` | [`scripts/m5_21_8_b_lattice.py#L56`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_21_8_b_lattice.py#L56-L87) |
| **the degree instrument** (the record's) | `read_charge_from_M` on the spatial 3x3 block | [`scripts/m5_22_e_audit.py#L192`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_22_e_audit.py#L192-L200), called through [`m5_32_r6_a_deltaladder.py#L156`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r6_a_deltaladder.py#L156-L169) |
| the R10 audit's own degree reader (a different lift and triangulation) | `directors`, `read_surface` | [`scripts/m5_32_r10_audit.py#L142`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r10_audit.py#L142), [`#L247`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r10_audit.py#L247) |
| the degree-0 vacuum-interior seed (the 78 % figure) | `unwound_seed` | [`scripts/m5_32_r10_audit.py#L260`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r10_audit.py#L260) |
| the melt paths, the clock taper, the `g = 32` probe (the R10 audit's typed results) | the `zero_barrier_robustness`, `clock_taper`, `scope_g32` blocks | [`scripts/m5_32_r10_audit.py#L1036`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r10_audit.py#L1036-L1092) |
| the note audit's independent re-derivation of all of the above | `path_energy`, `stage_barrier`, `stage_boundary`, `stage_taper`, `stage_g32` | [`scripts/m5_32_note_audit.py#L446`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_note_audit.py#L446-L708) |

## 3. The physics module

[`scripts/m5_32_lagrangian.py`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py)
is the single-purpose registry: each term has ONE definition string (hashed, so a term is never
silently re-tried under a different meaning), a sympy implementation on the notebook conventions,
and a numpy implementation on the certified stencil, plus per-term selftests. Every driver imports
it; no driver re-implements physics. Run `python3 scripts/m5_32_lagrangian.py --selftest` for the
17-line check, and `--mutant eta_time_row` for the negative control that must redden it.

## 4. Results, each with its pre-registered gate

| # | Result | Gate it was pre-registered against | Convergence evidence |
| --- | --- | --- | --- |
| R0 | The stack, the record and the author's 2026-08-17 Newton notebook all reproduce | selftests within 1e-3 of the record | 17/17 selftests; 10/10 record items at <= 2.2e-15; notebook fit `A = 863.733`, `B = 167.668`, sign `+`, to six digits |
| R1 | **The whole constant-coefficient current-order class is infeasible**: no coefficients make the energy's `omega^2` form PSD on every time channel with the boost weight reversed | a coefficient region existing under either Coulomb gate | Farkas / LP certificates at `(g, delta)` = (32, 0.3), (8, 0.3), (32, 0.1), with and without the parity-odd terms, on every channel alone, even at `c_I1 = +4` |
| R2 | The covariant flip family `L_lambda` is bounded below for `lambda >= 1/2` by a pointwise theorem wherever `M eta` has a real timelike eigenvector, keeps the static sector exactly, and gives `lambda* = 1/2` on every channel | G4, G5, the sector half of G1 | 0 negative densities in 27,560 random non-Lorentz samples; 36 channel x g cases; lattice probes bounded with no guard |
| R3 | G2 not met on any of three constructions | sign robust across 2 of 3 constructions, 2 boxes, both boundary types | ansatz repulsive at 348 reads; 44 relaxed pair heals; the cross-inertia undecidable at this resolution |
| R4 | The fixed-J minimizer runs to the box wall on every localized dressing family, `omega*` proportional to `1/L` | an interior `R*` with `omega*` stable across the domain ladder | 96/96 producer cases and every audit case; `omega* L` constant at 7.1 |
| R6 | **C3 orbit-blindness theorem**: any Lorentz-invariant derivative-free `V` is constant along a Lorentz dressing | a potential that localizes the dressing | variation <= 2e-7 on 50 dressed points and the whole R4 family up to rapidity 3; Euclidean controls O(1e4) |
| R7 | `K_T` moves the fixed-J minimizer off the box wall, but the audit found the interior minimum is the dressing switching off (0.43 % deep), and the term is exactly inert on the realized clock channel; G7 fails on the drift gate alone | interior `R*` with `omega*` drift <= 10 % over a c2 range >= factor 4 | interior at `c2` = 0.03 and 0.1 in both boxes; drift never below 0.301 against the 0.10 bar; the range half is MET on a dense ladder (factor 4.87) |
| R8 | C6's `omega^4` inertia is exactly VOLUME extensive and h-independent | an IR-convergent `omega^4` inertia | L exponent 3.0000 to 1e-13, ratio 8.000 over a factor 2 in L; h exponent -3.6e-14 |
| R9 | The ansatz carries a topologically protected biaxial disclination on the z axis; at `delta = 0` the field is exactly frame-free and the PERIODIC clock vanishes identically (a smooth radial-boost flow survives). The audit's headline: a relaxation resolves the line into a finite core (radius 3.98 at h = 1.5, 3.58 at h = 0.75) with the clock surviving at an h-convergent inertia 351.17 / 351.14, at `g = 8` | a string-free hedgehog with a nonzero clock | frame-free identity to 2.08e-17 relative; continuum ring spread exactly `delta / 2` and radius-INDEPENDENT over a 1e4 shrink, against a spread proportional to the radius at `delta = 0` |
| R10 | **No energy barrier protects the ansatz's degree, and the extensive inertia is not a property of an object** | the relaxed core-resolved soliton's inertia still extensive (the producer's registered prediction; the unwinding and boundary claims below are the AUDIT's findings and carry no pre-registered gate of their own) | from the UNRELAXED ansatz, a straight-line melt to the degree-0 state never rises above its start (energy 62.852 at the start, 14.794 at the end, \|Q\| 1 -> 0 inside each of five melt windows, 201 points each); the barrier from the RELAXED state was not computed, and the probes run from it rise 0.73 to 4.49; a degree-0 configuration with a vacuum interior out to `r = 15` still carries 78 % of the inertia (272.20 of 351.17, equal-budget unconverged comparison); a linear taper of the clock flow over `r = 12` to `15` leaves 32.9 % and removes the box dependence, as any compactly supported clock on a box-independent interior must |

### 4.1 The two results that reverse earlier readings

**The degree instrument is not an invariant of this order-parameter space.**
`read_charge_from_M` takes `V[..., -1]` of the spatial 3x3 block, the leading eigenvector alone, so
it measures an `RP^2` degree of one eigenvector. The stabilizer of `d4` in `SO(1,3)+` is the Klein
four-group (an `L` with `L^T eta L = eta` and `L d4 L^T = d4` must commute with `d4 eta`, whose
spectrum is distinct, so `L` is a diagonal sign matrix; `det = +1` and `L_00 = +1` leave four), so
`pi_1 = Q8` and `pi_2 = 0`. With `pi_2 = 0` and the eigenvalues frozen on the `SO(1,3)+` orbit of
`d4`, any CONTINUOUS map `S^2 -> OPS` has leading-eigenvector degree 0; the degree `+-1` reading is
possible only because the ansatz is discontinuous on the measurement surface. On the same three
surfaces the instrument calls conflict-free, the MIDDLE eigenvector admits no consistent sign lift
at all. What was measured about protection is narrower than a barrier: from the unrelaxed ansatz a
straight-line melt to degree 0 never rises above its start, while the relaxed state's own barrier
was not computed and FIRE holds `|Q| = 1` on every surface through 12000 iterations; the protection
argument is `pi_2 = 0`, not a measured barrier.

**The extensive clock inertia belongs to the boundary and the frozen-clock convention.**
`kin` is quadratic in `a0`, so the frozen (non-tapering) clock flow is only an UPPER BOUND on the
fixed-J inertia, and an upper bound that grows with `L` cannot establish that the inertia grows with
`L`. A linear taper of the clock flow to zero over `r = 12` to `15` leaves 32.9 % of it, and a
compactly supported clock on a box-independent interior is L-independent by construction.

## 5. Minimal inspection set (physics first, driver last)

| Order | Artifact | Why |
| --- | --- | --- |
| 1 | [`scripts/m5_32_lagrangian.py`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_lagrangian.py) | the functional: every term's definition, sympy and numpy side by side |
| 2 | [`scripts/m5_22_e_audit.py#L192`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_22_e_audit.py#L192-L200) | the degree instrument, because section 4.1 turns on what it measures |
| 3 | [`scripts/m5_32_note_audit.py`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_note_audit.py) | the independent instrument behind section 9: its own ansatz, energy, clock tangent and degree reader, so every load-bearing number here has a second implementation to read against |
| 4 | [`data/m5_32_ledger.json`](../data/m5_32_ledger.json) | every rung's hypothesis, claims and audited verdicts in machine form; the drivers are listed there and are the last thing to read |

## 6. What was NOT computed

| Not computed | Why it matters |
| --- | --- |
| Any candidate carried through the full G1 to G7 battery | no candidate survived far enough; the `lambda`-family died at G2 and G3, `K_T` at G3 and G7 |
| The classes C7 (higher-order timelike-current / Skyrme contractions) and C8 (cross-model imports) | never opened; R10's criterion says they cannot move G3, but that is an argument, not a measurement |
| A relaxed two-box ladder at the toy point `g = 32` | every relaxation here is at `g = 8`, where `V4` is 4096x softer; the `g = 32` probe reached `V4 = 0.00097` with no cell isotropized (the smallest top eigenvalue gap 0.616 against the 0.35 threshold that defines the melt front) but ended unconverged at `fmax` 5.14 |
| Whether a protected object exists in this space at all | `pi_1 = Q8` suggests a disclination LOOP rather than a point hedgehog; not tested |
| A converged relaxation anywhere | every FIRE run stops on `max_iter`; the audited slope decrements per doubling depend on which box pair is read (-0.267 / -0.403 on (24, 36), -0.159 / -0.130 on (36, 48)), so neither convergence nor decay is established |
| The unwinding barrier of the RELAXED state | only the melt from the unrelaxed ansatz was run; the probes from the relaxed state rise 0.73 to 4.49, and FIRE never moves the degree |
| A degree instrument that sees the time row | the instrument reads the spatial 3x3 block only, so anything the relaxation does to the time row is invisible to it |
| Converged comparisons behind the 78 % and 32.9 % figures | both are equal-budget, unconverged comparisons |
| The physical clock localization | which clock flow is the physical one is a convention question the run could not settle from inside |
| `J = hbar / 2` in program units | undefined in the record; never invented, so every fixed-J number is at an arbitrary `J` |
| The Coulomb pair half of G1 on the 4x4 stack | the like-charge static control fails on this stack (the string form), so the instrument could not decide it |

## 7. The adversarial audit record

Every rung was audited by an independent agent instructed to REFUTE, with its own implementation
(different stencil branch order, own amp grid, own minimizer, own densities) and forbidden from
reading the producer's scripts. The audits are the reason several headline claims below R7 no longer
stand as first written.

| Rung | Verdicts | What the audit changed |
| --- | --- | --- |
| R7 | 8 CONFIRMED, 5 QUALIFIED, 0 REFUTED | found that the LP channel list contains no channel built from the clock the model actually runs, so `c2` gives exactly zero help there; found the dressed-pair Coulomb anchors are not `c2`-independent |
| R8 | 4 CONFIRMED, 4 QUALIFIED, 2 REFUTED | found the ansatz's z-axis discontinuity and that 73 to 98 % of every C5 quartic coefficient sits beside it; refuted the producer's `c5` coefficient ladder as an h-artifact (`h^+2.99`) |
| R9 | 5 CONFIRMED, 2 QUALIFIED, 2 REFUTED | refuted the producer's exclusion theorem by RELAXING it: the line resolves into a finite core and the clock survives; established `pi_1 = Q8`, `pi_2 = 0` |
| R10 | 2 CONFIRMED, 3 REFUTED | refuted the degree's topological meaning, measured the unwinding barrier at exactly 0.0, showed the inertia belongs to the boundary, and scoped the whole core-melt effect to `g = 8` |

Producer errors caught and logged rather than buried: an off-center excision mask (built on a
cell-centered grid while the density lives on the certified offset grid); a spherical-shell
integration biased 18.5 % low because it discards the cube corners; a generator table built from
`X M - M X^T`, which is antisymmetric and not a tangent; a "converging to a nonzero value" claim
withdrawn by the producer on its own 12000-iteration point and sent to the auditor as a claim to
refute BEFORE that auditor ruled.

One certified-stack defect was found and is owed as a platform issue: `gen_catalog` normalizes `a0`
by `max(norm, 1e-300)`, so at `delta = 0` it returns a unit-norm noise field and reports a phantom
`kin` of 2.25. Any `delta -> 0` study routed through it sees a clock that is not there.

## 8. Open questions for the author

Each of these is author-gated in the strict sense: the run can measure around it but cannot settle
it from inside, because the answer is a statement of intent or of convention about the model.

| # | Question | Why it is author-gated |
| --- | --- | --- |
| 1 | Which clock localization is physical: the rigid co-moving flow, or one that decays away from the defect? | the run measured that the answer decides whether the clock inertia is extensive at all: tapering the flow at `r = 12` leaves 32.9 % of it and removes the box dependence entirely |
| 2 | Is the electron intended as a point hedgehog, given `pi_2 = 0` and `pi_1 = Q8` in this order-parameter space? | the protected objects here are lines, and a disclination loop is a different object |
| 3 | The intended reading of the mixed trace `R_ac` | exactly one independent mixed trace exists up to sign, and it is not symmetric, so `I4 != I5` |
| 4 | `J = hbar / 2` in program units | undefined in the record; every fixed-J number is at an arbitrary `J` without it |
| 5 | Whether a spectral function of `M` (a projector, `h_cov`) is admissible under the model's own boundary | `h_cov` is undefined past the degeneracy locus `t* = (g + 1) / 2`, where the spectrum of `M eta` goes complex |
| 6 | The earlier open questions carried from the M5.21 series | in the model's question tracker ([`m5_question_tracker.md`](../m5_question_tracker.md)), unchanged |

## 9. The pre-send audit of this note

Per the standard, an independent agent (a different model from the one that produced the run)
audited THIS DOCUMENT before it was sent: it re-derived the six load-bearing claims with its own
ansatz, stencil energy, clock tangent, kin density, degree reader (a different sign lift and
triangulation) and melt paths, traced every number in the note to its artifact, re-walked all rows
of the equation-to-code map on the working tree, and checked every section 1 equation against the
code term by term. Instrument: [`scripts/m5_32_note_audit.py`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_note_audit.py);
record: [`data/m5_32_note_audit.json`](../data/m5_32_note_audit.json).

| Load-bearing claim | Verdict | The auditor's own number |
| --- | --- | --- |
| the degree instrument reads one leading eigenvector | CONFIRMED, with a correction to the note's wording (it acts on the spatial 3x3 block; fed the 4x4 it raises) | `+-1` on three surfaces; the middle eigenvector admits no lift |
| stabilizer = Klein four-group, `pi_1 = Q8`, `pi_2 = 0` | CONFIRMED by its own derivation | stabilizer order 4; the SU(2) lift closes on a non-abelian group of order 8 containing `-1` |
| the unwinding barrier is exactly 0.0 | **QUALIFIED, and the note's lead claim was re-scoped on it**: the number reproduces exactly (five melt windows, 201 points each, endpoint 14.7940) but the path starts at the UNRELAXED ansatz; from the relaxed state the probes rise 0.73 to 4.49 and FIRE holds `\|Q\| = 1` through 12000 iterations | 0.0 from the rigid start; 0.725 / 3.037 / 4.490 from the relaxed states |
| a degree-0 vacuum-interior state carries 78 % of the inertia | CONFIRMED from its own seed | 272.204 / 351.170 = 77.5 % |
| a clock taper at `r = 12` leaves 32.9 %, L-independent | CONFIRMED, with the taper's shape made explicit (a linear ramp to zero over 12 to 15) | 115.385 / 351.170 = 32.86 %; L = 36 / 48 / 60 agree to 4 digits |
| at `g = 32` the relaxation leaves `V4 = 0.00097` and no melt | CONFIRMED | own `V4` 0.000970; smallest top gap 0.616 against the 0.35 threshold |

Beyond the six: all code-map anchors resolved to the named function; every section 1 equation
matched the code (`I1 = (1/2) F_abcd F^abcd` to 8e-15, `E_cert = 4(U + omega^2 T) + V4` to 2e-16,
both Legendre forms, `kin = -4 C(I1)`, the six quartics); the per-rung audit verdict counts matched
the audit records. Six blocking findings were applied before sending (the barrier re-scope above,
the instrument's input block, "monotone" replaced by "never above the start", the gate and class
legend added as section 1.6, two transcription errors corrected: the boost drift figure and a
mid-run `fmax` reported as an endpoint) and the fourteen minor findings were folded into the rows
they concerned. The audit did not edit this note; the producer applied its findings and this
section records them.
