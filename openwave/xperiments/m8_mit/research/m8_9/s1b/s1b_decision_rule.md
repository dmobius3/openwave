# M8.9 S1b decision rule: analytic-subspace action control

> **Status: FROZEN. S1b is NOT commissioned and this document does not authorize the run.** The
> frozen region is every byte above the boundary marker at the foot of this file, and is
> append-only from here.
> Seven adversarial design passes preceded this freeze; what they closed is recorded in the closing
> sections rather than tidied away, because several of the defects they caught were in the
> instrument's own gates.
>
> **The commissioned unit does IMPLEMENTATION QUALIFICATION ONLY.** Reproduce every gate and every
> mutation from this document, attack vacuity, and STOP on any contract contradiction. Do NOT
> redesign thresholds, and do NOT adjust anything in response to what the live target produces.
>
> Claim boundary, unchanged from M8.9: locate where the non-real spectrum enters. No repair, no
> nonlinear dynamics, no reopening of P1A. S2 remains untriggered and is not authorized by this
> document.

## Why the architecture changes

S1 closed `INSTRUMENT DEFECT`. Its architecture was: diagonalize the whole 60x60 block, sort the
spectrum into four windows, and infer which discrete eigenspace corresponds to which continuum
harmonic level. G-MULT failed, so that inference was never licensed. The assumption that died is
**that the finite block's global spectrum self-identifies its continuum levels**, and raising the
seed count does not remove that assumption, it only makes it more comfortable.

S1b removes it. The continuum eigenspaces are constructed INDEPENDENTLY of `L`, sampled on the
production cloud, and the shipped operator is then asked what it does to them. The continuum
statement is not "there ought to be eigenvalues near here", it is the far stronger

    L_cont Q_n = lambda_n Q_n,     lambda_n = n(n+2)

so the measurement is an action residual, not a cluster identification. No edge cluster, no Voronoi
assignment, and no dependence on `1 + 13 + 21 + 25 = 60` meaning anything about resolution.

**The trace result motivates this and does not adjudicate it.** `trace(L) = 32258.4` against the
four-level content `13(168) + 21(440) + 25(624) = 27024` is a window-free, basis-independent
statement that global spectral truncation failed, which is precisely why global spectral truncation
must not be the successor instrument. But the trace is nearly blind to what M8.9 is hunting: for a
real matrix the complex eigenvalues come in conjugate pairs and their imaginary parts cancel in the
trace EXACTLY. It is recorded as a required supporting diagnostic and adjudicates nothing.

**Scope of the trace comparator, frozen.** `27024` is the four-level analytic content of the
60-DIMENSIONAL block specifically. It is not a comparator at any other dimension. The trace check
is therefore recorded ONLY at the 60-seed baseline and is never reused at larger dimension, where
the block's analytic content is a different number and the four-level sum means nothing.

**Seed ladder, frozen now, before any singular value is seen.** G-SAMPLE is allowed to demand more
resolution, so the ladder cannot be chosen after seeing where sampling fails. It is

    60  ->  120  ->  180

with `k = 110`, `m = 7`, `p = 4` held at the P1A production values throughout, so seed count is the
only parameter that moves. 60 is the baseline that carries the trace comparator; the two rungs
double and triple it. S1b runs at 60 first. A rung above 60 is entered ONLY on a G-SAMPLE failure
at the rung below, never to improve a result that already qualified.

## The sampling realization map, written from the live contract

**There is no shipped map from the invariant basis to node values.**
`invariant_dim_and_basis` in `m8_5b/production/route_a_repn.py` returns
`(dimension, basis, conditioning_gap)` with the basis in the SU(2) symmetric-power COEFFICIENT
space. The only evaluation routine in the tree, `eval_monos` in `m8_5b/pilot/route_a_tuned.py`,
evaluates MONOMIALS, a different basis spanning a different space, since degree-`n` monomials on
S³ carry harmonics of degree `n`, `n-2`, `n-4` and so on. The two do not compose.

**Derived, then checked against the module, then machine-verified.** The module states the
contract, but this is new load-bearing code and the P0 precedent is derive-then-implement rather
than adopt-by-analogy. The derivation is short. On S³ isomorphic to SU(2) the degree-`n` harmonics
are the matrix entries of the `(n+1)`-dimensional irrep, so a function at level `n` is a linear
functional on those entries, `f(x) = c . vec(rho_n(x))`, an ordinary bilinear pairing. The pullback
by the two-sided action `x -> u x v` gives `rho_n(u x v) = rho_n(u) rho_n(x) rho_n(v)`, whose
column-major vec is `(rho_n(v)^T kron rho_n(u)) vec(rho_n(x)) = M_n vec(rho_n(x))`. Hence

    f(u x v) = c . (M_n vec(rho_n(x))) = (M_n^T c) . vec(rho_n(x))

so the COEFFICIENT action is `M_n^T` while the BASIS-VALUE action is `M_n`, and the pairing carries
no conjugation because it is a dot product and not an inner product. That derivation reproduces the
module's `REALIZATION` string and its `VEC_ORDER = "F"` exactly, which is a check on the derivation
rather than its source. G-REAL is the machine verification of it.

Frozen here, with the objects given distinct names so the contraction cannot drift:

    v_n(x)  =  sym_power(quat_to_su2(x), n).reshape(-1, order=VEC_ORDER)     VEC_ORDER = "F"

    C_n     =  invariant_dim_and_basis(pairs, n)[1]        shape ((n+1)^2, k), COLUMNS are the
                                                           basis vectors (`Vh.conj().T[:, dim-k:]`)

    V_n     =  the (N_nodes x (n+1)^2) matrix whose ROWS are v_n(x_i)^T

    F_n     =  V_n C_n                                     the raw sampled function values

**The contraction is a plain transpose, never a Hermitian conjugate**, for the reason the
derivation gives: the pairing is a dot product, not an inner product. So
`F_n^cover[i, j] = f_{c_j}(x_i)` with no conjugation anywhere. Inserting a `.conj()` here would be
silent and wrong, and would still produce a full-rank well-conditioned matrix.

**The wrong-realization mutation is the no-transpose form, NOT row-major.** The module header is
explicit that column-major `rho_n(v)^T kron rho_n(u)` and row-major `rho_n(u) kron rho_n(v)^T` are
the SAME operator in different vec conventions, and that the transpose on the right factor is
load-bearing rather than a vec artifact. Labelling a consistently implemented row-major convention
as wrong would be a false mutation. The genuine one is shipped: `no_transpose=True`, giving
`rho_n(u) kron rho_n(v)`, which is CHARACTER-EQUIVALENT, reproduces every invariant dimension
exactly, and fixes a maximally different subspace at `n >= 3`.

**THE n = 2 TRAP, and it is why the arm has a level floor.** The module records that at `n = 2`
the correct and no-transpose subspaces COINCIDE, so a regression that stops there passes under the
wrong realization. `verify_realization` sweeps `n >= 3` for exactly this reason. S1b's arm inherits
that floor. S1b's own targets are `n = 12, 20`, both safely above it, but the ARM must be run where
the two realizations demonstrably differ or it proves nothing.

**One caution the module raises that this rule must respect.** `coefficient_operator` records
narrowly that replacing `M_n` by `M_n^T` changed neither the invariant dimensions nor the resulting
invariant subspace at `n = 2..6` on the tuning groups, and that whether this is coincidence or a
theorem is NOT settled. So the coefficient transpose may not be separable by SUBSPACE comparison.
G-REAL therefore tests the POINTWISE law, which is the module's own stated evidence, and never
infers correctness from subspace agreement.

## Gates, in order. Every one must be able to fail

**G-REAL, pointwise realization. This is the load-bearing gate.** On the ACTUAL cover cloud, the
sampled invariant functions must satisfy deck invariance pointwise,

    f(gamma . x) = f(x)      for every frozen group element and every node

tested independently of `L`, in the manner of the existing `verify_realization`. Armed with the
shipped
`no_transpose=True` form, which is character-equivalent and dimension-preserving, swept at
`n >= 3` because the two coincide at `n = 2`.

Operational criterion, INHERITED from `verify_realization` rather than invented: pass requires
`worst_residual_correct < 1e-10` AND `best_residual_no_transpose > 1e-10`, against the recorded
separation of 1.1e-15 for the derived law and 4.9e+00 for the untransposed one. If G-REAL cannot
reproduce that separation, S1b stops and reports an instrument defect.

**The three stages are distinct objects and the gates attach to specific ones.** Writing them
separately is what stops a gate being satisfied by the step that follows it:

    C_n  --G-RANK-->  F_n^cover  --G-REAL-->  F_n^seed  --G-ALIGN, G-SAMPLE-->  Q_n  --> L_q

Only `F_n^seed` may feed G-SAMPLE, `Q_n`, `A_n`, `K_n` and `J_n`. `F_n^cover` exists to be tested
by G-REAL and for nothing else, since deck invariance is a statement about the cover and `L_q`
lives on the seeds.

**G-RANK, on `C_n`, before any sampling.** `invariant_dim_and_basis` at its shipped
`tol_rel = 1e-8` must return dimension exactly 13 at `n = 12` and exactly 21 at `n = 20`, and the
returned `gap` diagnostic is recorded with its `state` field.

**G-SUBSPACE, two independent constructions of the invariant subspace.** The subspace whose
compression is supposed to answer the science question is built TWICE, by algebraically different
routes on the same qualified representation matrices:

    C_n^svd  =  invariant_dim_and_basis(pairs, n)[1]        stacked-nullspace SVD, as shipped

    Pi_n     =  (1/|G|) sum_{gamma} coefficient_operator(u, v, n)      finite-group averaging
    C_n^avg  =  an orthonormal basis of ran(Pi_n)

Measured on the shipped code at the left action `pairs_left(elems)`, `|G| = 120`:

| n | invariant dim, SVD route | rank of `Pi_n` | max sin(principal angle) | SVD time | `Pi_n` time |
| --- | --- | --- | --- | --- | --- |
| 0 | 1 | 1 | 0.00e+00 | 0.0 s | 0.0 s |
| 3 | 0 | 0 | n/a, both empty | 0.1 s | 0.0 s |
| 12 | 13 | 13 | 2.58e-08 | 415 s | 0.1 s |
| 20 | 21 | 21 | 2.11e-08 | 3919 s | 0.3 s |

Three things follow, and none of them was assumed.

**`Pi_n` is idempotent but NOT Hermitian**, `‖Pi^2 - Pi‖` of 7.5e-09 at `n = 12` and 2.3e-07 at
`n = 20`, against `‖Pi - Pi^H‖` of 3.9e+01 and 4.5e+02. It is an OBLIQUE projector, because
`sym_power` in this basis is not unitary. The gate therefore orthonormalizes `ran(Pi_n)` and
compares SUBSPACES by principal angles. Treating `Pi_n` itself as an orthogonal projector would be
wrong.

**Rank must be taken at an ABSOLUTE cutoff.** A first version of this check used a relative cutoff
and reported rank 16 at `n = 3`, where the true invariant dimension is 0. `‖Pi_3‖_2 = 4.1e-18`, so
the matrix is numerically zero and a relative threshold divides by nothing. Frozen: rank of `Pi_n`
at absolute cutoff `1e-8`. This is the same defect as the ladder's zero overlap, in a different
place, and both are now guarded.

**Two angles, in two different spaces, and only one of them controls the operator-level
uncertainty.** The measured 2.58e-08 and 2.11e-08 are COEFFICIENT-space angles. Write

    theta_C(n) = theta_max( ran C_n^svd , ran C_n^avg )                in coefficient space
    theta_Q(n) = theta_max( ran Qtilde_n^svd , ran Qtilde_n^avg )      in the space Ltilde acts on

with `Qtilde = M_h^{1/2} Q`. These are NOT the same angle. The chain
`C -> F^seed -> W = M_h^{1/2} F^seed -> Qtilde` can amplify a coefficient-space disagreement, and
G-SAMPLE's `kappa(W_n^r) <= 1e6` says amplification is POSSIBLE without making the two angles
equal.
`K_floor`'s discrepancy term therefore uses `theta_Q`, the principal angle of the subspaces
ACTUALLY being compressed, since the derivation `‖Q_1 - Q_2 U‖_2 = 2 sin(theta/2)` is about those
subspaces.

**The 8.34e-05 and 6.82e-05 figures are PROVISIONAL**, computed from `theta_C` because `theta_Q`
cannot be known before the sampling runs. They establish only that the discrepancy term is likely
to dominate the `kappa` term at 7.18e-06 and the arithmetic term at 1.3e-11. The operative value is
whatever `theta_Q` measures.

**`theta_Q` is recorded and NOT gated, deliberately.** Gating it would need a threshold that cannot
be set before the run without seeing it, which is the failure this rule exists to avoid. It is
instead self-regulating: a large `theta_Q` enlarges `K_floor`, which makes `K` harder to resolve,
which is the correct consequence. Route disagreement in the `J` or `K` readings is caught
separately by rules 3 and 5 of the branch algorithm. If `theta_Q` is large enough that `K_floor`
swamps any plausible `K`, S1b reports that outcome rather than adjusting anything.

The G-SUBSPACE tilt arm is ALSO run in `Qtilde`-space and reported as a capability demonstration
that the comparison discriminates there, with no gate attached to its result.

**Gate and routing.** `G-SUBSPACE` passes when both routes return the same dimension, rank of
`Pi_n` is taken at ABSOLUTE cutoff `1e-8`, and `sin(theta_C(n)) <= 1e-6`. The gate is on
`theta_C`, which certifies that the two ALGEBRAIC constructions identify the same
representation-theoretic space; `theta_Q` is measured later and is not gated. `A_n`, `K_n` and
`J_n` are
then formed from BOTH `Q_n^svd` and `Q_n^avg`, and construction agreement is enforced by rules 3
and 5 of the branch algorithm, not by prose here.

**Armed, because `1e-6` was chosen after seeing 2.6e-08 and a threshold picked that way has to
prove it has power.** Two arms, both with analytically controlled answers.

*Subspace arm.* From the known-good `C_n`, take a unit vector `w` orthogonal to `ran(C_n)` and
replace one column by `c_0 cos(phi) + w sin(phi)`, then re-orthonormalize. The principal angles of
the result against the original are `{phi, 0, ..., 0}`, so `sin(theta_max) = sin(phi)` exactly.
Verified at `d = 169`, `k = 13`:

| frozen `sin(phi)` | measured `sin(theta_max)` | required |
| --- | --- | --- |
| 1e-04 | 1.000e-04 | RED |
| 1e-05 | 1.000e-05 | RED |
| 1e-06 | 9.997e-07 | green |
| 1e-07 | 9.541e-08 | green |

The `1e-04` tilt is the frozen arm and must fail the gate. The `1e-07` tilt is the frozen
green control and must pass. Together they show `1e-6` discriminates rather than merely sitting
38 times above what was observed.

*Rank arm.* At the absolute cutoff `1e-8`, a numerically-zero matrix at the `‖Pi_3‖ ~ 1e-17` scale
must return rank 0, and a known rank-`k` projector must return `k`. Verified: 0 and 13
respectively, where the RELATIVE cutoff returns 16 and 13, which is the defect that produced the
spurious `n = 3` rank in the first place.

The `1e-6` remains a declared engineering tolerance, fixed before any target quantity is computed,
but it is now a tolerance with demonstrated power rather than a number with headroom.

**Practical note, not a criterion.** The averaging route is roughly 1.3e+04 times faster at
`n = 20`, 0.3 s against 65 minutes. An earlier draft suggested dropping the SVD route if a seed
rung is entered. That was wrong: `invariant_dim_and_basis(pairs, n, ...)` takes only the group
pairs and the level, with no cloud and no seed count, so `C_n^svd` is SEED-INDEPENDENT. The
65-minute cost is paid once and the basis is REUSED at 120 and 180 seeds; only the sampling changes.
Neither route may be dropped, because the scientific branch now requires two-construction
agreement.

**G-ALIGN, the cover-to-seed restriction, before anything is orthonormalized.** `L_q` is the
QUOTIENT operator on `seed_orbits = sorted(plan.keys())`, so the restriction is frozen as

    node_of, mult, plan = orbit_stencils(X, oid, gid, elems, k=110)   # same k as the operator
    L_q, seed_orbits    = build_L_bundle(X, oid, gid, elems, rho, k=110, m=7, p=4)
    F_n^seed[a, :]      = F_n^cover[ plan[seed_orbits[a]][0] , : ]

**The gate is a structural predicate on the representative-index vector, not a numerical effect.**
Define

    r = ( plan[seed_orbits[a]][0] )_{a = 0 .. N-1}

and require, exactly:

    r == ( node_of[(o, 0)] )_{o in seed_orbits}       both verified true on the shipped code
    seed_orbits == sorted(plan.keys())                verified true on the shipped code
    F_n^seed == F_n^cover[r, :]                       by construction, asserted

**Armed structurally.** Change exactly ONE entry of `r` to a different valid cover-node index,
leave every other entry fixed, and require the alignment predicate to fail. Known-green parent,
mechanically red child.

An earlier draft armed this by permuting the rows of `F_n^seed` and requiring `D_n` to move by a
factor of 10. That is withdrawn, and the constant it introduced is deleted. There is no theorem
that every non-identity permutation raises the action residual: a swap of two numerically similar
representative rows could move `D_n` very little, or conceivably lower it, and a correctly aligned
implementation would then fail its own arm because the MUTATION lacked power rather than because
the gate was broken. G-ALIGN has an exact structural predicate and does not need a downstream
numerical effect to be testable.

**G-SAMPLE, on the RAW `F_n^seed`, before `M_h`-orthonormalization exists.** This ordering is the
whole
gate. Once `Q_n` has been successfully orthonormalized its conditioning is good BY CONSTRUCTION, so
a conditioning test applied to `Q_n` is satisfied by the step that produced it and proves nothing.
That would be the sixth member of the green-for-the-wrong-reason family this program has already
paid for five times.

Define, PER ROUTE, the object that is actually orthonormalized, and qualify each of them. One
route being well sampled does not qualify the other, so BOTH must pass:

    W_n^r = M_h^{1/2} F_n^{seed, r},     r in {svd, avg}

`M_h`-orthonormalizing `F_n^seed` IS ordinary orthonormalization of `W_n`, so both criteria are on
`W_n`'s singular values. Exact rank is unchanged by the invertible `M_h^{1/2}`, but NUMERICAL rank
at a relative cutoff is not, so qualifying `F_n^seed` and then orthonormalizing `W_n` would again
be qualifying the wrong object. Frozen criteria:

- numerical rank of `W_n^r` exactly `k` at relative cutoff `1e-8`, inheriting
  `invariant_dim_and_basis`'s own `tol_rel` so the two stages use one convention;
- `kappa(W_n^r) = sigma_max / sigma_min <= 1e6`.

Both criteria must hold for BOTH `r`. The precision ladder likewise runs independently on
`A_n^svd` and `A_n^avg`, which the branch algorithm already assumes.

The `1e6` is derived, not chosen for comfort: orthonormalization in float64 loses about
`log10(kappa)` digits, so `kappa = 1e6` leaves roughly ten, which is the margin the `K` and `J`
floors below need in order to sit far above their own arithmetic noise. If a 21-dimensional
analytic subspace collapses under sampling the run stops HERE, before `L` is applied, and that is
the only circumstance in which spending compute on more seeds has a principled purpose: qualifying
the sampling of a KNOWN subspace, not hoping a global spectrum looks more cluster-like.

**G-BASIS, basis invariance. Run separately on `Q_n^svd` and `Q_n^avg`, and BOTH must pass.**
Replace `Q_n^r -> Q_n^r U` with `U` the unitary factor of a QR
factorization of a Gaussian matrix drawn from a FROZEN seed, `numpy.random.default_rng(20260825)`,
so the check is deterministic and rerunnable. **Tested as covariance, not as relative scalar
agreement.** A relative tolerance is undefined
when the reference vanishes, which is exactly possible for `K_n` and `J_n`, whose continuum targets
are zero. So the check is on the objects and is norm-scaled:

    ‖A'_n - U^H A_n U‖_F  <=  1e-10 * ‖A_n‖_F
    ‖P'_n - P_n‖_F        <=  1e-10 * max(1, ‖P_n‖_F)

An earlier draft justified a bare `1e-10` by claiming `‖P_n‖ = 1` for the `M_h`-orthogonal
projector. That is wrong twice: an orthogonal projector of rank `k` has Frobenius norm `sqrt(k)`,
not 1, and an `M_h`-orthogonal projector need not even have Euclidean operator norm 1 in the
original coordinates. The bound above is norm-scaled instead, and it uses FROBENIUS
deliberately: G-BASIS tests aggregate covariance of the whole compression, which is a different
question from the spectral one `K_floor` adjudicates, so the norms differ on purpose rather than by
oversight. `‖K‖_2` and `J` are then covariant
consequences rather than separately compared scalars. Any dependence above that is a defect, not a
tolerance question.

**G-DISCRIM, and it exists because `K` and `J` are different predicates.** See below.

**G-WIRE, machinery control at `n = 0`.** The constant mode with `lambda_0 = 0` verifies that the
construct, sample, orthonormalize, apply-`L` pipeline is wired correctly end to end. Criterion,
INHERITED from the merged S1 record rather than invented: `‖L Q_0‖_{M_h,F} <= 1e-8`, against the
5.7e-10 measured there.

`max abs Im` on the constant mode is RECORDED as a regression against S1's 5.954e-14 and the
shipped P1A 5.95e-14, and is NOT gate evidence. For a real `L`, `M_h` and `Q_0` the constant
compression is a real 1x1 scalar, so its imaginary part is zero by algebra whatever the pipeline
did. The residual and its mutation are the parts of G-WIRE that carry evidence.

**Its arm, and the first version of this was vacuous.** An earlier draft armed G-WIRE by
permuting the sampled node order of `F_0`. `F_0` is the CONSTANT harmonic, so every row is the
same scalar and `P F_0 = F_0` for every permutation `P`. The mutation could not change the object
and would have been the sixth green-for-the-wrong-reason arm in this program. It is replaced by
one whose effect is analytically forced. After `M_h`-normalizing `Q_0`, apply

    L_mut = L + delta * P_0,     P_0 = Q_0 Q_0^H M_h,     delta = 1e-4 frozen

Since `Q_0^H M_h Q_0 = I`, this gives `L_mut Q_0 = L Q_0 + delta Q_0` EXACTLY, so against the
green criterion `‖L Q_0‖ <= 1e-8` the mutated residual satisfies `‖L_mut Q_0‖ >= 1e-4 - 1e-8`,
four orders above the gate. Known-green parent, analytically forced red child.

**It is explicitly NOT the calibration for `K`**: `A_0` is a 1x1 real scalar, so `K_0 = 0` by
algebra whatever the higher-dimensional pipeline did. Verified: an arbitrary real 1x1 block returns
`‖K‖_2 = 0.000e+00` while an arbitrary real 13x13 returns a nonzero value. Using `n = 0` to
certify the `K`
machinery would be green for the wrong reason, which is the failure family that has now cost this
program five separate arms.

## What is measured

> **Route indexing, governing from here on.** Unless expressly stated otherwise, every `Q_n`,
> `P_n`, `A_n`, `D_n`, `E_n`, `K_n`, `J_n`, every precision-ladder result and every diagnostic below
> is evaluated SEPARATELY for `r` in `{svd, avg}`. Singular notation in the definitions is for
> readability; there is no singular object. `K_floor(n)` is the one deliberate exception, being a
> common worst-case bar shared by both routes.

For `n` in `{12, 20}`, after `M_h`-orthonormalizing the independently sampled invariant basis, with
`P_n = Q_n Q_n^H M_h` and `A_n = Q_n^H M_h L Q_n`:

| quantity | definition | role |
| --- | --- | --- |
| `D_n` | `‖(L - lambda_n I) Q_n‖_{M_h,F} / ‖lambda_n Q_n‖_{M_h,F}` | supporting, non-adjudicating |
| `E_n` | `‖(I - P_n) L Q_n‖_{M_h,F} / ‖L Q_n‖_{M_h,F}` | supporting, non-adjudicating |
| `‖K_n‖_2` | `K_n = (A_n - A_n^H)/2`, OPERATOR 2-norm | **self-adjointness witness** |
| `J_n` | `max abs Im mu` over the spectrum of `A_n` | **non-real-spectrum witness** |
| `trace(L)` | against `13(168) + 21(440) + 25(624)` | supporting, non-adjudicating |

`D_n` and `E_n` measure ACCURACY, and there is no inherited accuracy scale for high-level action,
so no principled a priori threshold exists for them. They are recorded and they do not route.

## `K` and `J` are different predicates, and the rule must not collapse them

Both have an exact continuum target of zero. They are not the same statement.

    A = [[1, 1],
         [0, 2]]        ‖K‖_2 = 0.5,  J = 0,  eigenvalues 1 and 2, both real

So `K_n != 0` certifies only that the trivial-fibre discrete action already fails to be
`M_h`-self-adjoint on a known harmonic subspace. It does NOT certify that the trivial fibre
reproduces the non-real spectral phenomenon M8.9 is trying to localize. `J_n` is what asks the
M8.9 question.

**G-DISCRIM, the arm that proves the implementation knows the difference.** From a known-green
Hermitian
parent `A_green = lambda I`, apply two frozen mutations whose answers are analytically exact. These
are SYNTHETIC checks on exact matrices, so the criteria are machine allowances rather than physical
tolerances, and they are named rather than left implicit. The ANALYTIC answers are `J = 0`, `0`
and `2`; the GATE tests their numerical realization, because no eigensolver is obliged to return a
bitwise-zero imaginary part on every BLAS or LAPACK build:

    green parent:   ‖K‖_2 <= 100 eps ‖A_green‖_2      and   J <= 1e-12
    arm A:          | ‖K‖_2 - 0.5 | <= 1e-12          and   J <= 1e-12
    arm B:          | ‖K‖_2 - 2   | <= 1e-12          and   | J - 2 | <= 1e-12

| arm | mutation | required |
| --- | --- | --- |
| green parent | `lambda I` | analytic `‖K‖_2 = 0`, `J = 0` |
| A | `lambda I + N`, `N` nilpotent, `N[0,1] = 1` | analytic `‖K‖_2 = 0.5`, `J = 0` |
| B | `lambda I + S`, `S` real antisymmetric, `S[0,1] = -S[1,0] = 2` | analytic `‖K‖_2 = 2`, `J = 2` |

Arm A is the one that matters: it is non-Hermitian with a provably real spectrum, so an
implementation that equates "non-Hermitian" with "complex spectrum" fails it. Arm B confirms the
`J` path fires when a conjugate pair genuinely exists. Both targets are exact, so neither arm can
be satisfied by a tolerance chosen afterwards.

## The `K_n` floor: the right adjoint defect, and why `kappa` and this floor are ONE problem

**The parent anti-self-adjoint object is not `L - L^H`.** From `A_n = Q_n^H M_h L Q_n`,

    A_n - A_n^H = Q_n^H ( M_h L - L^H M_h ) Q_n

so the compressed object is the `M_h`-adjoint defect, and `L - L^H` is the right answer only when
`M_h = I`. Equivalently, writing `Q_n = M_h^{-1/2} Qt_n` with `Qt_n` ORDINARY-orthonormal,

    A_n = Qt_n^H Ltilde Qt_n,        Ltilde = M_h^{1/2} L M_h^{-1/2}

so `A_n` is the ordinary compression of `Ltilde` onto an ordinary orthonormal basis, and the scale
that propagates is `‖Ltilde - Ltilde^H‖_2`. Verified on the shipped 60-seed `R_0` operator:
`‖Q^H M_h L Q - Qt^H Ltilde Qt‖ = 7.6e-13` and `‖Q^H M_h Q - I‖ = 1.7e-15`.

An earlier draft used `‖L - L^H‖_2 = 3230.88`. The correct quantity is
`‖Ltilde - Ltilde^H‖_2 = 3232.03`. **The correction is required and, at 60 seeds, nearly inert**:
`M_h` here is very well conditioned, `kappa(M_h) = 1.209` with diagonal entries between 2.45e-03
and 2.97e-03, so the similarity barely moves the norm. It is written the right way anyway, because
the near-agreement is a property of THIS mass matrix at THIS seed count and carries no guarantee
for another rung of the ladder or for a nontrivial `rho` whose `M_h` is worse conditioned.

**Why this floor and `kappa` are one problem.** A subspace admitted at `kappa eps` carries an
uncertainty in WHICH subspace `Qt_n` represents, and because the parent is strongly
non-self-adjoint, a small error in `Qt_n` EXPOSES part of that defect inside the compression.
Numerically at 60 seeds:

| target | arithmetic `100 eps ‖A_n‖_2` | `kappa` term | discrepancy term, provisional | ratio, largest over arithmetic |
| --- | --- | --- | --- | --- |
| `n = 12`, `k = 13` | 3.73e-12 | 7.18e-06 | 8.34e-05 | 2.2e+07 |
| `n = 20`, `k = 21` | 9.77e-12 | 7.18e-06 | 6.82e-05 | 7.0e+06 |

Recomputed in the operator 2-norm, since that is the norm the rule adjudicates in. An earlier
draft tabulated `100 eps ‖A_n‖_F`, giving 1.34e-11 and 4.48e-11, which for `A ~ lambda I` is
`sqrt(k)` larger than the operative value. The history is worth keeping: the gap between the
arithmetic floor and the real uncertainty was five orders of magnitude when only the `kappa` term
was known, and is seven once the measured discrepancy term is included.

Freezing the arithmetic floor alone would assert a distinction five orders of magnitude below what
the admitted subspace precision supports. The floor therefore carries both terms:

**Every downstream object is ROUTE-INDEXED.** Two constructions means two of each:
`W_n^svd, W_n^avg`, `Q_n^svd, Q_n^avg`, `A_n^svd, A_n^avg`, `K_n^svd, K_n^avg`, `J^svd, J^avg`.
In general `‖A_n^svd‖_2 != ‖A_n^avg‖_2` and `kappa(W_n^svd) != kappa(W_n^avg)`, so a singular
`K_floor` would silently adjudicate the two routes against different budgets. One COMMON,
worst-case floor is used instead, so both routes face the same bar:

    B          = Ltilde - Ltilde^H,   r ranges over {svd, avg}
    K_floor(n) = max( 100 * eps * max_r ‖A_n^r‖_2 ,
                      10  * eps * max_r kappa(W_n^r) * ‖B‖_2 ,
                      2 * sin( theta_Q(n) / 2 ) * ‖B‖_2 )

**All three terms and the witness are now in the OPERATOR 2-norm.** An earlier draft adjudicated a
Frobenius witness against operator-norm uncertainty terms. For a `k x k` compression
`‖X‖_F <= sqrt(k) ‖X‖_2`, so that comparison understated the uncertainty by up to `sqrt(13) = 3.61`
at `n = 12` and `sqrt(21) = 4.58` at `n = 20`. The spectral norm is also the sharper question:
is there ANY direction in this analytic subspace with a resolved anti-self-adjoint component?
Frobenius aggregation is not what self-adjointness failure needs.

**The third term uses `2 sin(theta_max/2)`, which is the correct alignment distance and is also
tighter.** For subspaces at principal angle `theta_max`, the optimally aligned basis difference is
`‖Q_1 - Q_2 U‖_2 = 2 sin(theta_max/2)`, and since `K(Q) = (1/2) Q^H B Q` the induced change is
bounded by `2 sin(theta_max/2) ‖B‖_2`. Checked against directly measured perturbations at
`theta_max` of 1e-03, 1e-04 and 1e-05: worst measured-over-bound ratio 0.629 for this form, so it
holds, against 0.315 for the looser `2 sin(theta_max)` version an earlier draft used.

At 60 seeds the third term DOMINATES, 8.34e-05 at `n = 12` and 6.82e-05 at `n = 20`, against
7.18e-06 for the `kappa` term and 1.3e-11 for the arithmetic term. Without it, `‖K_n‖_2` could
clear the floor purely because two independent constructions of the subspace disagree, which is not
a property of the operator.

**It is a discrepancy ALLOWANCE, not a certified error bound.** What is measured is the
disagreement between two independent numerical constructions of the same exact invariant space. Two
methods agreeing does not prove the exact subspace lies between them, so this term is named an
independent-construction discrepancy allowance and never a bound on distance to the continuum
subspace. The stronger practical evidence is the branch-agreement requirement: both constructions
must independently yield the same scientific disposition.

A common floor makes branch rules 5 and 6 clean: both `‖K_n^svd‖_2` and `‖K_n^avg‖_2` above it is
ADJOINT-eligible, exactly one above is NO_LABEL, neither above is no resolved defect.

`‖K_n^r‖_2 > K_floor(n)` is a self-adjointness defect; at or below it, `K_n^r` is consistent with
zero
AT THE ACHIEVABLE SUBSPACE PRECISION, which is weaker than "consistent with zero" and is the honest
statement. `kappa <= 1e6` survives only as an admissibility bound in G-SAMPLE; the floor adapts to
the measured value rather than assuming it away.

**Two stability checks, and they measure different things.** Neither is a full re-derivation, and
the rule says so rather than overclaiming.

- **Orthonormalization stability, per route.** For each `r`, re-derive `Q_n^r` at 50 digits from
  the SAME binary `C_n^r` and `F_n^{seed, r}`, and record the movement in `K_n^r` and `J_n^r`. This
  bounds the arithmetic of `M_h`-orthonormalization ONLY. It does NOT bound uncertainty in the
  invariant basis or in the sampling, because both are held fixed.
- **Independent construction, which replaces an earlier and ill-defined proposal.** A previous
  draft said to perturb `C_n` "inside its own certified nullspace scale, using the `gap`
  diagnostic". That was undefined: `invariant_dim_and_basis` returns `gap = min(s_above) /
  max(s_below)`, a dimensionless SEPARATION RATIO with special states for the empty and exact
  cases, described in the source only as how cleanly the nullspace separates from the rest. It is
  not a norm in coefficient space and not an allowed perturbation amplitude. It is withdrawn.
  G-SUBSPACE above replaces it with a second, algebraically different construction of the same
  subspace, which measures the uncertainty instead of inferring a scale for it.

A genuine high-precision re-derivation of `C_n^svd` is still not proposed: it needs a 50-digit SVD
of a `(120 (n+1)^2) x (n+1)^2` stack, 52920 x 441 at `n = 20`, where the float64 SVD alone takes 65
minutes. G-SUBSPACE makes that unnecessary rather than merely standing in for it, because two
independent constructions bound the ambiguity directly.

`n = 0` cannot calibrate any of this, for the algebraic reason recorded under G-WIRE.

## The `J_n` floor, and the ladder is NOT threshold-free

The exact target is `J_n = 0`, so no empirical contamination envelope is used and P1A's fitted
envelope is explicitly not the model. But an earlier draft of this rule claimed the precision
ladder "needs NO threshold at all", and that was wrong: a 50-digit eigensolve is still a numerical
computation, and "persists unchanged" is not self-defining. Both statements below are frozen now.

**Precision ladder, primary, and it must cover FORMATION and not only the eigensolve.** P1A held
`L` byte-fixed and raised only the eigensolver precision, which was correct there because `L` was
itself the object under investigation. Here `A_n = Q_n^H M_h L Q_n` is a NEWLY FORMED measurement
object. Holding a float64 `A_n` fixed and eigensolving it at higher precision would prove only that
the complex pairs belong to the ROUNDED `A_n`, not that they survived its formation. For a
nonnormal or nearly defective compression, entry-level rounding during formation can manufacture
conjugate pairs outright.

So from the SAME fixed binary inputs `Q_n`, `M_h`, `L`, the matrix product is RE-FORMED at each
rung and then eigensolved at the matching precision:

    rung 1: form A_n in float64, eigensolve in float64        -> J(64)
    rung 2: re-form A_n at 30 dps, eigensolve at 30 dps       -> J(30)
    rung 3: re-form A_n at 50 dps, eigensolve at 50 dps       -> J(50)

The dps sequence is inherited from P1A's localization test rather than chosen here. The criterion
is stated as an ALGORITHM with ordered branches rather than as a predicate table, because table
predicates overlap: an earlier draft's `J(50) >= 0.5 J(64)` and `J(50) <= 1e-3 J(64)` are BOTH true
when `J(64) = J(30) = J(50) = 0`, since each reduces to `0 >= 0` and `0 <= 0`, so PERSISTS and
COLLAPSES both fired at exact zero. Ordered, first match wins:

    1.  if J(64) == 0 and J(30) == 0 and J(50) == 0:      -> COLLAPSES   (exact zero, no ratios)
    2.  if J(50) < 0.5 * J(30) or J(50) > 2 * J(30):      -> AMBIGUOUS   (tail not settled)
    3.  if J(50) >= 0.5 * J(64):                          -> PERSISTS
    4.  if J(50) <= 1e-3 * J(64):                         -> COLLAPSES
    5.  otherwise:                                        -> AMBIGUOUS

Rule 1 takes exact zero out of the ratio tests entirely. Rule 2 makes the 30-dps rung a voting
member rather than decoration: an earlier draft listed three rungs and compared only two. Every
AMBIGUOUS routes to `S1b-NO_LABEL`, never to a branch.

Armed by the same discriminating control P1A used: the companion matrix of `(x-1)^k`, whose exact
spectrum is real, must COLLAPSE under this ladder. A control that persists means the ladder is
mis-wired and S1b stops.

**Certified perturbation bound, secondary and demotable.** The eigensolver is backward stable, so
the computed spectrum is exact for `A_n + E`, and Bauer-Fike bounds the eigenvalue displacement by
`kappa(V) ‖E‖` with `V` the eigenvector matrix. Recorded as
`J_bf(n) = kappa(V) * 100 * eps * ‖A_n‖_2`. **This bound is demoted to a recorded diagnostic, with
no routing power, whenever `kappa(V) > 1e8`**, because near defectiveness it becomes loose enough
to admit anything and `V` may not be trustworthy at all. Which regime holds is decided by
`kappa(V)`, a property of `A_n`, and the demotion rule is frozen here rather than after seeing it.

## Branch logic. Five outcomes, distinguished rather than collapsed

**Precedence is global across both targets, evaluated in this order, first match wins.** Without
it the branches are not deterministic: `n = 12` collapsing with `K > K_floor` while `n = 20` is
AMBIGUOUS could be read as either ADJOINT or NO_LABEL. The frozen order:

Every quantity is computed on BOTH constructions from G-SUBSPACE, `Q_n^svd` and `Q_n^avg`, and
construction agreement is an explicit PREDICATE rather than prose sitting above the adjudicator.
An earlier draft stated the agreement requirement in G-SUBSPACE but omitted it from this table, so
the two contradicted each other: a case with `J` collapsed on both routes, `‖K‖_2 > K_floor` on the
SVD route and `<= K_floor` on the averaging route, was NO_LABEL by G-SUBSPACE and ADJOINT by the
table. Frozen, first match wins:

| # | condition | outcome |
| --- | --- | --- |
| 1 | any instrument gate fails | `S1b-DEFECT` |
| 2 | any target has an AMBIGUOUS `J` on EITHER construction | `S1b-NO_LABEL` |
| 3 | the two constructions disagree on the `J` reading for any target | `S1b-NO_LABEL` |
| 4 | any target is PERSISTS on BOTH constructions | `S1b-SPECTRAL` |
| 5 | all `J` collapsed or zero on both, and the constructions disagree on `‖K_n‖_2 > K_floor(n)` for any target | `S1b-NO_LABEL` |
| 6 | all `J` collapsed or zero on both, and some target has `‖K_n‖_2 > K_floor(n)` on BOTH | `S1b-ADJOINT` |
| 7 | else | `S1b-NULL` |

Rules 2 and 3 sit above rule 4 deliberately: an ambiguous or route-dependent second target must not
be hidden behind a clean-looking reading on the first. Rule 5 does the same for `K`.

**S1b-DEFECT.** G-REAL, G-RANK, G-SUBSPACE, G-ALIGN, G-SAMPLE, G-BASIS, G-DISCRIM or G-WIRE
fails. Instrument defect.
No `K` or `J` interpretation is issued.

**S1b-SPECTRAL.** `J_n` is QUALIFIED PERSISTENT for at least one of `n = 12, 20`, meaning it
returns PERSISTS from the ladder AND agrees across both independent invariant-subspace
constructions. Licensed: the trivial fibre alone can produce a non-real compressed action, so
nontrivial fibre transport is NOT NECESSARY for that phenomenon, and the base discretization or the
scalar quotient reduction is strongly implicated.

An earlier draft made this branch conditional on `J_n` exceeding "the certified floor". That phrase
had no referent once the Bauer-Fike bound became demotable to a non-routing diagnostic. The branch
predicate is the ladder plus subspace reproducibility, which is a stronger statement than clearing
a loose perturbation bound, and the Bauer-Fike value is recorded alongside without routing power.

**S1b-ADJOINT.** `‖K_n‖_2 > K_floor(n)` on both constructions, but `J_n` consistent with zero.
Licensed, and no more: the
trivial-fibre discrete action already has a self-adjointness defect on a known harmonic subspace,
and S1b did NOT reproduce the non-real spectral phenomenon. This is a bounded result and must not
be promoted to a localization.

**S1b-NO_LABEL.** Reached by rule 2, 3 or 5: an ambiguous ladder result on either construction, a
`J` reading that differs between constructions, or a `K` verdict that differs between them. S1b
issues NO scientific reading. This outcome exists so that
an intermediate ladder result is reported as undecided rather than argued into SPECTRAL or NULL,
which is precisely the failure P1A's fitted envelope produced at larger scale.

**S1b-NULL.** Both consistent with zero. Licensed, and no more: trivial-fibre high-level action
alone is insufficient. This does not exonerate RBF-FD and does not convict the fibre transport,
since the defect may be an interaction that vanishes at trivial fibre by construction.

**Routing to a possible reduction-mechanism comparison: BOTH `S1b-ADJOINT` and `S1b-NULL`.** An
earlier draft made NULL the only such branch. That was a policy choice presented as a mathematical
consequence, and it is withdrawn. `S1b-ADJOINT` says the trivial-fibre substrate carries a
self-adjointness defect without yet producing a non-real spectrum, which is exactly the substrate
an interaction with the twisted reduction could act on to produce one. Excluding it would rule out
the most plausible interaction hypothesis on no stated grounds. Both branches therefore leave M8.9
OPEN with the same disposition: a separately frozen interaction or reduction comparison MAY be
commissioned. Neither branch commissions one, and neither authorizes S2 as filed.

## What pass 7 closed, and what is still open

Both fixes have one root cause: introducing a second construction made several downstream objects
ambiguous while the notation still read as though there were one.

- **Two angles, in two spaces.** `theta_C` compares the coefficient-space subspaces and is what
  G-SUBSPACE gates at `1e-6`, certifying that the two algebraic constructions identify the same
  representation-theoretic space. `theta_Q` compares `ran Qtilde^svd` against `ran Qtilde^avg` in
  the space `Ltilde` acts on, and it is what `K_floor`'s discrepancy term must use, since the
  `2 sin(theta/2)` derivation is about the subspaces actually compressed. The chain
  `C -> F^seed -> W -> Qtilde` can amplify, and `kappa(W) <= 1e6` says amplification is possible
  without making the angles equal. **The 8.34e-05 and 6.82e-05 figures are now PROVISIONAL**: they
  were computed from `theta_C` and establish only that the discrepancy term is likely to dominate.
- **`theta_Q` is recorded and not gated.** Gating it would require a threshold set after seeing it.
  It is self-regulating instead, since a large `theta_Q` enlarges `K_floor`, and route disagreement
  is already caught by branch rules 3 and 5.
- **Everything downstream is route-indexed**, `W_n^r`, `A_n^r`, `K_n^r`, `J^r`, and `K_floor` is a
  single COMMON worst-case bar taking `max_r` over both routes. A singular floor would have
  adjudicated the two routes against different budgets, since `‖A^svd‖_2 != ‖A^avg‖_2` and
  `kappa(W^svd) != kappa(W^avg)` in general. G-SAMPLE now requires BOTH routes to pass rank and
  conditioning; one being well sampled does not qualify the other.
- **The last unnamed tolerance is named.** G-DISCRIM's green parent requires
  `‖K‖_2 <= 100 eps ‖A_green‖_2` and `J <= 1e-12`, with each arm's target to `1e-12` absolute.
  The analytic answers stay `0`, `0.5`, `2`; the gate tests their numerical realization, since no
  eigensolver is obliged to return a bitwise-zero imaginary part on every build.
- **The SVD route is reused, not dropped.** `invariant_dim_and_basis(pairs, n, ...)` takes only the
  group and the level, so `C_n^svd` is seed-independent; the 65-minute cost is paid once and the
  basis carries to 120 and 180 seeds unchanged. An earlier draft had this backwards, and dropping
  it would have destroyed the two-construction agreement the branch table now depends on.

Still open, all declared engineering choices:

1. **Constants `10`, `100` and the `2` in `K_floor`**, `1e-4` in the G-WIRE mutation, `1e-6` in
   G-SUBSPACE with its `1e-4` red arm and `1e-7` green control, the `1e-8` absolute rank cutoff,
   and `1e-12` for the G-DISCRIM arm targets.
2. **`kappa(W_n^r) <= 1e6`**, admissibility only.
3. **The seed ladder `60, 120, 180`**, preregistered, entered only on a G-SAMPLE failure.
4. **`theta_Q` has no gate.** A deliberate omission, not an oversight; it feeds the floor instead.
<!-- FREEZE-BOUNDARY -->

**Freeze record.** SHA-256 covers every byte ABOVE the boundary comment: `c44c603a8058ed8529e5bb0f42ec168b443ff22b050d3748f1f1de3537c7d297`

```bash
sed '/^<!-- FREEZE-BOUNDARY -->$/,$d' S1B_DECISION_RULE.md | shasum -a 256
```
