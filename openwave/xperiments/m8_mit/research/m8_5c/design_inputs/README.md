# M8.5-C design inputs: the five verified scripts

Secured 2026-08-27 from volatile locations (`/tmp`, `~/Downloads`). These are the checkable
half of the chassis decision memo (frozen `44c664d1…`): the memo carries their results in its
Verified-inputs table; these files are the provenance. They ship with the M8.5-C protocol.

| script | verifies | key numbers |
| --- | --- | --- |
| `mode_count.py` | per-sector Galerkin mode counts from the 2I character table; self-checks: characters orthonormal, first occurrence `n = d_ρ` mult 1, all nine sectors | memo cost table; asymptotic `d_ρ N³/360` |
| `exact_quad_check.py` | degree-`4N` Hopf product rule returns the cubic Galerkin projection to rounding; node-drop mutation armed | `2.2e-15` vs `1.5e-01` (2N rule); Gram `3.1e-14`; cross-level `1.8e-15`; `sym_power` basis non-unitary |
| `cascade_quad_check.py` | the cascade monitor's band `N < n ≤ 3N` needs its OWN `6N`-exact rule; production `4N` aliases the band | `6.0e-15` vs `2.35e-01`; band-norm reading 1.4% off under 4N; node ratio 3.05x measured |
| `jacobian_check.py` | the fluctuation operator of `R(φ;ω) = c²Δφ − c1⟨φ,φ⟩φ + ω²φ` is REAL-linear on `R^{2m}` (a complex Jacobian drops the `⟨δφ,φ⟩φ` term, rel diff 1.000); `iφ` is an exact kernel vector INSIDE the scored eigenspace; gauge-breaking mutation goes red | kernel residual `2.5e-09`; mutation `1.3e-01` |
| `right_translation_check.py` | script 5 (per D-5): arms the seven theorem-to-code bridges of `../../findings/m8_5c_symmetry_derivation.md` plus two armed preflights (P0 group forensics with a non-icosian substitution mutation; P1 character table with R5 DERIVED in-room by column orthogonality and a unique sign solution, perturbed-R5 mutation), all routed through the same report() so the ledger records structure, never a traceback; each check has its own real green parent and real mutation: C1 right-action realization (l = 1..7, 12); C2 assembled-vs-analytic coefficient rep (kron-swap mutation); C3 cubic equivariance (corrected sides, mutation under the 6N rule so redness is the symmetry break, never aliasing); C4 complete-level necessity; C5 multiplicity one at every scored level, BOTH halves of step 6 (Molien levels {0,12,20,24,30,32,36} by character + Reynolds; the eight sectors with first occurrence DISCOVERED over a fixed scan n = 0..10 and compared to the pinned distances, so neither the distances nor R5's character row is a transcription); C6 2I-commutant census (dims = <chi_l,chi_l> = {1,1,1,1,1,2,2,4} at l = 1..7,12, commutative iff multiplicity-free; pi1+pi1 mutation shows dim 4 noncommutative); C7 level-diagonal spectrum vs the assembled right action (renamed per redline: the assembled object is the action, the spectrum is analytic; intra-level-broken mutation). Prints provenance: resolved module path, its SHA-256, numpy/scipy versions | all nine report lines armed, wall ~3 s; the 120-icosian group is explicit with UNIQUENESS, FULL 120x120 closure, and the order census {1:1,2:1,3:20,4:30,5:24,6:20,10:24} |

## Dependency and run commands

`exact_quad_check.py` and `cascade_quad_check.py` import `quat_to_su2` and `sym_power` from the
repo module `openwave/xperiments/m8_mit/research/m8_5b/pilot/route_a_nonabelian.py`:

```bash
PYTHONPATH=<repo>/openwave/xperiments/m8_mit/research/m8_5b/pilot python3 exact_quad_check.py
PYTHONPATH=<repo>/openwave/xperiments/m8_mit/research/m8_5b/pilot python3 cascade_quad_check.py
python3 mode_count.py
python3 jacobian_check.py
```

All four re-run PASS from this directory 2026-08-27; script 5 re-run PASS 2026-08-28 at its
final pinned form (nine report lines: P0, P1, C1
through C7; every parent green, every mutation red; exit 0; R5 derived, first occurrence
discovered).

## SHA-256 at securing

```
7f31bc305150c9d5e37ec312fe8fd00060583a87427feb7a56d6743d2fd66bfa  cascade_quad_check.py
51fa9b52bc9d3367fbb20baf084235da260c7f15b36fad1b78b31770e7f55ed7  exact_quad_check.py
93eeb001a555a086758d341e13be543ce373fe5a3585a177535a669db01918da  jacobian_check.py
794d7063bea01ad3637a96c33e4e01e7fd71acfc636eb2d443c4f3cecb6959a6  mode_count.py
752baa2de36dbaea9fe4108ec6df9121351b386453d01b06bd8111562de8f547  right_translation_check.py
```

Those are the values AT SECURING. `jacobian_check.py` moved after it, at review, to
`dce687ea793522be40c4d06d25df19a0dfb99c72309a19361e220077007eb9d8`: its two verdict labels
printed unconditionally, so the gauge-breaking arm could not fail. The shipped pin for all
five is the protocol's § 15 table.

Run command for script 5 (same dependency):

```bash
PYTHONPATH=<repo>/openwave/xperiments/m8_mit/research/m8_5b/pilot python3 right_translation_check.py
```

Development record worth keeping: the script's first two runs went RED on C2 and C3 because of
MY OWN convention errors in the coefficient representation (first `P^T ⊗ I` for a row-major
vec, then the missing contragredient: field coefficients transform by `I ⊗ P`, no transpose).
The derivation was never in doubt; the bridges were, twice, which is exactly what per-bridge
arming is for. Those conventions are now pinned by a check that fails if they drift. Round two (both
redlines): the C3 mutation had computed R_g N_w(psi) on BOTH sides, so its redness was
quadrature aliasing of a degree-(4N+1) integrand under the 4N rule, not the symmetry break;
C2's mutation and C4's parent were hard-coded constants; and the first C6 measured the SU(2)
commutant of Sym^n, which Schur makes 1 for every n, so it could not fail. All four were
rebuilt: C3's sides corrected and run under the 6N rule (fires at 1.2e+00), C2/C4 given real
arms, and C6 replaced by the 2I-commutant census, whose mutation (pi1+pi1, commutant M_2)
fires on noncommutativity.

SCOPE, two layers (do not conflate): green here = the repo SCALAR primitives realize the
derived right action. The W_rho-valued sector bases and the real Galerkin system do not exist
yet; they are qualified later by protocol gates 3 and 5, and this script gives them no
pre-emptive credit.

## Jacobian findings NOT in the frozen memo (postdate the freeze; protocol must encode them)

The fluctuation operator is real-linear on `R^{2(d_ρ+1)}`; `iφ` (the U(1) gauge mode of the
standing-wave ansatz) is an exact kernel vector sitting inside the scored eigenspace; structural
zeros are branch-dependent. Consequences for the protocol: a predicted-vs-measured zero-count
gate; scored splitting measured on the symmetry-quotiented complement; phase fixing or a bordered
system for Newton; prereg § 0's complexification is load-bearing for continuation.
