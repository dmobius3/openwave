# Convention refactor: index-3 -> index-0 ("flip the living engine, freeze the history")

Decision (the maintainer, 2026-06-21): align the M5 storage convention to Duda's `D = diag(g, 1, δ, 0)`
(time/g axis = array index 0, `eta = diag(-1,1,1,1)`). Flip the LIVING engine + active sandboxes;
FREEZE historical sandboxes as index-3 records. N0 Test B already proved the relabel is physics-
neutral; the golden master proves the WHOLE engine is.

Baseline commit (pre-flip): `8ea9633`.

## The flip rule (exact)

Cyclic right-shift of the array index by 1: every engine index `k -> (k+1) mod 4`.

| engine (index-3, now) | -> new (index-0, Duda) | meaning |
| --- | --- | --- |
| index 0 (eigenvalue 1, r̂ principal) | index 1 | EM/principal spatial axis |
| index 1 (eigenvalue δ, e_Θ) | index 2 | δ-twist axis |
| index 2 (eigenvalue 0, e_Φ) | index 3 | null axis |
| index 3 (eigenvalue g, time/boost) | index 0 | time/g axis (gets the eta minus) |

So `D = diag(1,δ,0,g) -> diag(g,1,δ,0)`. Spatial block `[0:3,0:3] -> [1:4,1:4]`; `g` `[3,3] -> [0,0]`;
the metric minus and every hardcoded `3` for the time axis `-> 0`; 4-vectors put the time component
at index 0 (was 3). `eta = diag(1,1,1,-1) -> diag(-1,1,1,1)`.

## Scope

| Layer | Action | Files |
| --- | --- | --- |
| Engine COMPUTE path (golden-master-gated) | FLIP | `medium.py`, `engine1_seeds.py`, `engine2_pde.py`, `engine3_observables.py` |
| Engine viz/launch (smoke-test) | FLIP | `engine4_render.py`, `instrumentation.py`, `_launcher.py`, `xparameters/` |
| Active sandboxes | FLIP | `sandbox_v9` (#200), `sandbox_v10` (neutrino: n0 only; n1/n2 already index-0) |
| Historical sandboxes | FREEZE (header label, no logic edit) | `sandbox_v1..v8`, `sandbox_vn` (default; confirm) |

## Procedure (golden-master gated)

1. ✅ Golden master CAPTURE (`golden_master.py`) , 57 invariants frozen (`golden_master_baseline.json`).
2. 🔶 Flip the 4 COMPUTE files consistently (a half-flip is inconsistent; do them as one unit).
3. Golden master CHECK (`python3 golden_master.py --check`) , MUST PASS (tol 1e-4). This both
   recompiles every kernel (catches shape/syntax errors) AND proves physics unchanged.
4. Flip render/launcher/instrumentation/xparameters; viz smoke-test (import + one headless frame).
5. Flip sandbox_v9 + sandbox_v10/n0; re-run #200 + n0 spot checks.
6. Freeze sandbox_v1..v8 (+ vn) with a one-line header label.
7. Standardize the human-readable notation everywhere to `diag(g,1,δ,0)` + the implementation note
   (the doc/audit pass); fix any contradicting docstrings (the m5_9_4 class).

## Gate

The flip is ACCEPTED only when `golden_master.py --check` PASSES (worst rel diff < 1e-4). If any
invariant moves, a site was flipped wrong , bisect by file.

## Compute-path flip , site checklist (DONE ✅)

| File | Sites flipped |
| --- | --- |
| `engine1_seeds.py` | `embed4` (spatial->[1:4], g->[0,0]); `seed_dressed_hedgehog_M` 4-vectors (time comp -> index 0) |
| `engine2_pde.py` | `rebuild_M_from_director` embed; `V_M`/`dV_M`/`compute_tstar`/`dV_M_dressed` spatial-block extraction [0:3]->[1:4] + writebacks; integrator time-freeze; `eta_twist_masked` + `eta_twist` ((α,3)->(α,0), range(1,4)); `signed_dot4` ((p==3)->(p==0)); `sample_v03_drift` + `sample_p03_drift` (d/p[a_,3]->[a_+1,0]); `evolve_M_4d` drift removal; `apply_p03_clamp` |
| `engine3_observables.py` | none (uses `pde.V_M` + Frobenius curvature, index-blind) |
| `medium.py` | none (comments only; the `[i,j,3]` hits are render quad-vertex indices, not matrix axes) |

## Status
Step 1 ✅ (golden master captured, 57 invariants).
Step 2 ✅ (4 compute files flipped) + Step 3 ✅ **GOLDEN MASTER PASS, worst rel diff 4.09e-07**
  => the index-3 -> index-0 flip is PROVEN physics-neutral. The constrained clock integrator
  (the riskiest 23-site path) matches to ~4e-7.
### Render-path cascades (the golden master did not originally cover these)

Found by careful reading, NOT by the compute suite. Both fixed + now covered by the extended
golden master:

| Cascade | Fix |
| --- | --- |
| `eigen_decompose` (engine2) passed the FULL 4×4 to `principal_director`/`eigvec_for` (which read [0:3,0:3]) | extract the spatial `[1:4,1:4]` block first |
| amplitude tracker `d_vac` (engine3) was `diag(δ,δ,1,g)` | -> `diag(g,δ,δ,1)` (g at index 0) |

Golden master EXTENDED to 66 invariants (added eigen_decompose director + spatial eigenvalues +
amplitude tracker). Re-verified via git-stash pre-flip baseline:

**GOLDEN MASTER PASS (66/66), worst rel diff 3.84e-07; the 9 render invariants bit-identical
(rel 0.0). `render.dir_probe_x=1.0` (hedgehog director = radial, unchanged).**

ENGINE FLIP COMPLETE + PROVEN (medium + engine1/2/3): compute path 4e-7, render/tracker bit-identical.

## Status , ALL CORE STEPS ✅

| Step | Status |
| --- | --- |
| 1. Golden master capture | ✅ 66 invariants (compute + render) |
| 2-3. Flip engine compute path (medium/engine1/2/3) + render cascades (eigen_decompose, d_vac) | ✅ |
| 3. Golden master CHECK | ✅ **PASS 66/66, worst rel 3.87e-07; render invariants bit-identical** |
| 4. engine4_render/instrumentation/_launcher | ✅ no matrix-index code sites (read derived fields) |
| 5. sandbox_v9 (#200) m5_9_4 + m5_9_5 (docstrings + comps) | ✅ smoke-tested: boost_GEM dressed=1.549e-6, undressed=0 exactly |
| 5. sandbox_v10/n0 bridge | ✅ updated to index-0, re-PASS (fidelity 3.4e-5, invariance bit-identical) |
| 6. Freeze sandbox_v1-v8 + vn | ✅ untouched (perfect provenance), documented in `FROZEN_SANDBOXES.md` |
| 7. Comment audit (engine + #200) | ✅ swept to index-0 (perl + manual); canonical note `CONVENTION.md` |

### Step 7b , research-doc notation sweep ✅ (2026-06-21)

| Doc | Action |
| --- | --- |
| `99_summary_report`, `5a_lagrangian_evolution`, `0b_M5_roadmap` | convention BANNER (dense `(α,3)` narratives describe index-3-era work accurately; banner = "read `(α,3)` as `(α,0)`") |
| `9a_lepton_mass_planning` | banner + spec lines fixed (g on index 0, eta=diag(-1,1,1,1)) |
| `10a_neutrino_oscillations` | banner (Duda quotes left verbatim, predate the flip) |
| `sandbox_v9/m5_9_lepton_mass_clock_findings` | banner (active #200) |
| `0a_background`, `1a_framework`, `m5_question_tracker` | no change , already use `diag(g,1,δ,0)` |
| `4a_/4c_convo` conversation logs | LEFT verbatim (historical record) |
| frozen `sandbox_v1-v8`, `sandbox_vn` | untouched , covered by `FROZEN_SANDBOXES.md` |

CONVENTION REFACTOR COMPLETE. Code flipped + proven (golden master 66/66), #200 re-runs unchanged,
n0 re-passes, sandboxes frozen, comments + docs swept. Canonical note: `CONVENTION.md`.
