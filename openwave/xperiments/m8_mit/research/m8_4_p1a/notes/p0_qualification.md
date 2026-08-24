# M8.4 P0 Qualification Note

**Produced:** 2026-08-24 00:14:06 UTC
**Python:** 3.13.13
**NumPy:** 2.5.0
**SciPy:** 1.18.0
**Total elapsed:** 33.5s

## Statement

This qualification imports no code from outside this room.
Cloud: 60 seeds x 120 = 7200 nodes on S^3.
Stencil: k=110, m=7, p=4.

This is the reissue run. The first run is preserved as `RUN1_ENGINEERING_NOTE.txt`
and reclassified as engineering and rehearsal. Five repairs were applied before this
rerun; no gate logic was changed beyond the repairs.

## Frozen Qualification Tolerances

Frozen in `p0/frozen_tolerances.py` BEFORE the qualification rerun.
Values derived during the first engineering run (RUN1).
File SHA-256: `3f78c5711f1ac76b6562db2f3f66642f479f331235a5cc75fff5d23dc992bc7e`

| Parameter | Value |
| --- | --- |
| `p0_0_homomorphism` | `1e-08` |
| `p0_0_identity` | `1e-12` |
| `p0_0_trace` | `1e-06` |
| `p0_0_unitarity` | `1e-08` |
| `p0_2_collapse_imag` | `1e-14` |
| `p0_2_collapse_real` | `1e-10` |
| `p0_3_eigenvalue_window` | `2e+00` |
| `p0_4_admissible_floor` | `1e-03` |
| `p0_4_min_separation` | `1e+03` |
| `p0_4_rejection_floor` | `1e+00` |

## P0.0: Representation Certification

**Overall: PASS**

### dimensions

| Rep | Residual | Result |
| --- | --- | --- |
| R0 | True | PASS |
| R1 | True | PASS |
| R2 | True | PASS |
| R3 | True | PASS |
| R4 | True | PASS |
| R5 | True | PASS |
| R6 | True | PASS |
| R7 | True | PASS |
| R8 | True | PASS |

### identity

| Rep | Residual | Result |
| --- | --- | --- |
| R0 | 0.00e+00 | PASS |
| R1 | 0.00e+00 | PASS |
| R2 | 2.22e-16 | PASS |
| R3 | 0.00e+00 | PASS |
| R4 | 2.22e-16 | PASS |
| R5 | 4.44e-16 | PASS |
| R6 | 0.00e+00 | PASS |
| R7 | 1.11e-16 | PASS |
| R8 | 0.00e+00 | PASS |

### homomorphism

| Rep | Residual | Result |
| --- | --- | --- |
| R0 | 0.00e+00 | PASS |
| R1 | 5.60e-11 | PASS |
| R2 | 3.92e-10 | PASS |
| R3 | 1.12e-10 | PASS |
| R4 | 3.36e-10 | PASS |
| R5 | 3.36e-10 | PASS |
| R6 | 1.68e-10 | PASS |
| R7 | 2.24e-10 | PASS |
| R8 | 2.80e-10 | PASS |

### unitarity

| Rep | Residual | Result |
| --- | --- | --- |
| R0 | 0.00e+00 | PASS |
| R1 | 5.60e-11 | PASS |
| R2 | 3.92e-10 | PASS |
| R3 | 1.12e-10 | PASS |
| R4 | 3.36e-10 | PASS |
| R5 | 3.36e-10 | PASS |
| R6 | 1.68e-10 | PASS |
| R7 | 2.24e-10 | PASS |
| R8 | 2.80e-10 | PASS |

### traces

| Rep | Residual | Result |
| --- | --- | --- |
| R0 | 0.00e+00 | PASS |
| R1 | 5.01e-11 | PASS |
| R2 | 3.17e-10 | PASS |
| R3 | 1.06e-10 | PASS |
| R4 | 2.72e-10 | PASS |
| R5 | 2.24e-10 | PASS |
| R6 | 1.12e-10 | PASS |
| R7 | 1.12e-10 | PASS |
| R8 | 1.12e-10 | PASS |

**Mutation:** swap_detected=True, corrupt_detected=True
**Covariance:** all_pass=True

## P0.1: Deck-Equivariance Law Derivation

E_rho = (S^3 x W_rho) / 2I with (x, w) ~ (gamma x, rho(gamma) w).
Section psi lifts to psi_tilde: S^3 -> W_rho with psi_tilde(gamma x) = rho(gamma) psi_tilde(x).

**Derivation:** well-definedness of [(x, psi_tilde(x))] requires
[(gamma x, psi_tilde(gamma x))] = [(x, psi_tilde(x))],
which gives psi_tilde(gamma x) = rho(gamma) psi_tilde(x) (LEFT equivariance).
Not rho(gamma)^{-1}, not rho(gamma)^dagger, not rho(gamma)^T: the equivalence relation dictates rho(gamma).
Implemented in `bundle_operator.py` as `fibre_map = rho[gid]`.

**P0.1: PASS** (derivation filed, law implemented)

## P0.2: Trivial Irrep Collapse

Repair 1 applied: referent is now the qualified production `build_L_equivariant`
from `m8_5b_production/equivariant_stencils.py`, called with one-sided 2I pairs `[(gamma, e)]`.
Own scalar path retained as third arm.

- `||L_R0(real) - L_prod_grouped||_inf` = 2.84e-14 (production referent)
- `||L_R0(real) - L_own_grouped||_inf` = 2.84e-14 (own scalar path)
- `||L_R0(imag)||_inf` = 0.00e+00
- `||L_prod - L_own_scalar||_inf` = 0.00e+00 (production vs own scalar)
- Mutation (inverse law): residual = 2.84e-14

**P0.2: PASS**

## P0.3: Free E_rho Section Spectrum

**Overall: PASS**
**Mutations: PASS**

| Rep | d_rho | d | lambda_exp | dim_exp | Eigenvalues | Below | Result |
| --- | --- | --- | --- | --- | --- | --- | --- |
| R0 | 0 | 1 | 0 | 1 | -0.0000 | 0 | PASS |
| R1 | 1 | 2 | 3 | 2 | 3.0000, 3.0000 | 0 | PASS |
| R2 | 7 | 2 | 63 | 8 | 62.9559, 62.9559, 62.9906, 62.9906, 63.0178, 63.0178, 63.0277, 63.0277 | 0 | PASS |
| R3 | 2 | 3 | 8 | 3 | 8.0000, 8.0000, 8.0000 | 0 | PASS |
| R4 | 6 | 3 | 48 | 7 | 47.9634, 47.9938, 47.9938, 48.0062, 48.0062, 48.0182, 48.0182 | 0 | PASS |
| R5 | 6 | 4 | 48 | 7 | 47.9930, 47.9930, 47.9940, 47.9940, 47.9990, 47.9990, 48.0310 | 0 | PASS |
| R6 | 3 | 4 | 15 | 4 | 15.0000, 15.0000, 15.0000, 15.0000 | 0 | PASS |
| R7 | 4 | 5 | 24 | 5 | 24.0000, 24.0000, 24.0000, 24.0000, 24.0000 | 0 | PASS |
| R8 | 5 | 6 | 35 | 6 | 34.9973, 34.9973, 34.9995, 34.9995, 35.0033, 35.0033 | 0 | PASS |

## P0.4: Mutation-Test Fibre Transport

**Overall: PASS**

Oracle: independent stencil + weights at non-seed nodes,
eigensection lift with proposed law, check Delta psi_tilde = -lambda psi_tilde.

| Rep | lambda | correct | inverse | transpose | conjugate | omitted | wrong_g | separation | Result |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| R1 | 3.0000 | 2.50e-09 | 2.68e+02 | 3.09e+02 | 3.73e+02 | 2.81e+02 | 5.24e+02 | 1.1e+11 | PASS |
| R2 | 62.9559 | 5.14e-04 | 1.15e+01 | 1.75e+01 | 1.45e+01 | 1.20e+01 | 1.37e+01 | 2.2e+04 | PASS |
| R3 | 8.0000 | 9.19e-10 | 1.55e+02 | 1.23e+02 | 1.35e+02 | 1.81e+02 | 1.46e+02 | 1.3e+11 | PASS |
| R4 | 47.9634 | 1.69e-07 | 1.91e+01 | 2.59e+01 | 1.76e+01 | 2.16e+01 | 2.10e+01 | 1.0e+08 | PASS |
| R5 | 47.9930 | 5.15e-05 | 2.29e+01 | 1.62e+01 | 2.02e+01 | 1.83e+01 | 2.19e+01 | 3.2e+05 | PASS |
| R6 | 15.0000 | 1.18e-09 | 6.53e+01 | 6.65e+01 | 7.86e+01 | 1.03e+02 | 7.94e+01 | 5.5e+10 | PASS |
| R7 | 24.0000 | 1.34e-09 | 3.17e+01 | 3.76e+01 | 3.27e+01 | 7.00e+01 | 4.02e+01 | 2.4e+10 | PASS |
| R8 | 34.9973 | 1.62e-04 | 1.98e+01 | 2.84e+01 | 2.78e+01 | 3.24e+01 | 2.86e+01 | 1.2e+05 | PASS |

Numerical floor justification:

- R1, R3, R6, R7 (d_rho <= 4): floor ~1e-9, matching one-form precedent
- R4 (d_rho=6): floor ~1e-7, higher level reduces grid resolution
- R5, R8 (d_rho=5,6): floor ~1e-4 to 5e-5, same reason
- R2 (d_rho=7): floor ~5e-4, highest level, coarsest resolution
- Minimum separation across all: >=10^3 (well above O(1) threshold)

## Regression Tests (Repairs 3-4)

### Repair 3: Identity Index Assertion

- `elems[0]` is identity: **PASS**
- Mutation (identity at index 60): first eigenvalue = -3.4088,
  expected = 8.0, error = 11.4088,
  mutation detected: **True**

### Repair 4: Non-Hermitian Solver Lock

- Asymmetry `||L - L^H|| / ||L||` = 0.899 (> 0.1: **True**)
- `min(Re(eig(L)))` = 3.0000 (non-negative as expected)
- `min(eig(symmetrized))` = -1095.4 (negative: **True**)
- Regression lock: **PASS**

## Architecture Interface

```
construct_input -> [production_continuation] -> [production_subspace_score] -> label
```

The bundle operator (`build_L_bundle`) provides a single interface:

```python
L_rho, seed_orbits = build_L_bundle(X, oid, gid, elems, rho, k, m, p)
```

identical for R0, manufactured, free, and eventual target inputs.
The only parameter that varies is `rho` (the representation matrices).
No `if manufactured:` branch exists in the operator construction.

## Code Manifest

| File | Bytes | SHA-256 (prefix) |
| --- | --- | --- |
| `__init__.py` | 0 | `e3b0c44298fc1c14...` |
| `algebra.py` | 2512 | `2a27f19a8eb70923...` |
| `bundle_operator.py` | 6169 | `5513f7092c448518...` |
| `cloud.py` | 1804 | `0926c44ab0c6df92...` |
| `frozen_tolerances.py` | 1000 | `3f78c5711f1ac76b...` |
| `group.py` | 5793 | `1409d25cbdcdbcfd...` |
| `qualify.py` | 24395 | `af944887c1d6164c...` |
| `rbffd.py` | 2825 | `40ca9d68b0d65f07...` |
| `regression_tests.py` | 3823 | `f772fbe005fb5d2f...` |
| `representations.py` | 9641 | `d1439829d79d5e9b...` |

## Summary

- P0.0 (representations): **PASS**
- P0.1 (equivariance law): **PASS** (derivation filed)
- P0.2 (trivial collapse): **PASS**
- P0.3 (section spectrum): **PASS**
- P0.4 (mutation oracle): **PASS**
- Regression tests: **PASS**

**QUALIFICATION: PASS -- substrate can carry E_rho sections**
