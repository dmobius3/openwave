# M5.9 lepton mass spectrum (mu, tau): corrected build plan (#200)

> **Index convention (2026-06-21).** The M5 engine now stores the order parameter at INDEX-0: `D = diag(g, 1, δ, 0)`, `eta = diag(-1, 1, 1, 1)` (time/g axis = array index 0). The engine itself was flipped to match the notation this doc already adopts for the new code (proven physics-neutral: [`_convention_refactor/CONVENTION.md`](../_convention_refactor/CONVENTION.md)). Any legacy `(α,3)` / `index 3` reference below now means `(α,0)` / `index 0`.

## Status

Issue #200 is OPEN (board: Next). The first pass (M5.9.1 + M5.9.3) correctly DIAGNOSED the missing ingredient (a core-volume confiner / Higgs-like vacuum) but used a setup too simplified to compute the masses. Dr. Duda reviewed the posted result and corrected the direction. This doc captures his review, evaluates it, inventories what the real engine already provides, and scopes the full production engine build.

> **Execution status (2026-06-20): the full production engine build RAN, locally.** Sub-tasks 1-6 are complete on the production 4x4 engine. Scripts: `m5_9_4_engine_lepton.py` (the signed-energy ledger with the negative boost-GEM gravity + clock-oscillation contributions; minimization) and `m5_9_5_higgs_clock_spectrum.py` (the Higgs/confiner scale-selection, the #220 dynamical clock scan, the spectrum attempt). Results + the honest frontier map are written up in [`sandbox_v9/m5_9_lepton_mass_clock_findings.md`](../findings/m5_9_lepton_mass_clock_findings.md) § "M5.9.4-5: the full production engine 4x4 build". Step-by-step progress in [`sandbox_v9/checkpoints/`](../findings). All LOCAL: nothing posted to the issue or committed yet (the maintainer commits/PRs, then publishes the distilled writeup to #200). Headline: the negative gravity term is real and dominant; the Higgs/confiner does not robustly select a scale on the production V (open); the hierarchy origin stays input; the dynamical clock readout is stability-limited. A convincing documented framework + a precise open-problem map (matching mu/tau in one pass was never the bar).
>
> **Round 2 (2026-06-20): Dr. Duda corrected the PARAMETER REGIME and recommended a pivot. #200 PARKS.** Duda's reply (reviewing `sandbox_v9/m5_9_1_lepton_mass_law.py`, the original toy) identified the load-bearing error: our parameter regime is wrong by ~20 orders of magnitude. `delta` should be ~`1e-10` (the QED quantum-phase contribution, a tiny correction), not 0.3; `g` should be ~`1e10` (the rest-mass scale), not 8. Both the toy AND the production build used delta=0.30, g=8, which is why we get huge values and no hierarchy. He also flagged the potential needs a proper Landau-de Gennes TENSOR form (the production build already uses one; the exact form for the vortex-core regularization is the open piece). **He recommended doing less but more rigorously: pivot the first serious deliverable to NEUTRINO OSCILLATIONS** (topological vortices only, no point-like; the oscillation parameters are a guessing game today, so deriving them from the deeper theory is genuinely novel and article-worthy). The charged-lepton-mass target (#200) PARKS behind the shared blocker (the parameters g, delta + the regularizing potential), which Duda just began to pin (g~1e10, delta~1e-10). The neutrino pivot is scoped in **[`m5_10a_neutrino_oscillations.md`](m5_10a_neutrino_oscillations.md)** and filed as issue [#236](https://github.com/openwave-labs/openwave/issues/236) (board: In progress; Duda's round-2 reply captured verbatim in `10a`). #200 moved to **Next**; it resumes once g/delta/potential are locked.
>
> **Round 3 (2026-06-21): the blocker is an ACTIVE SEARCH, not a wait; eta docstring fixed.** Duda confirmed that the exact `g`, `delta`, and (worse) the potential are **surprisingly hard to find , even for him** ("you should start here, finding these parameters/details for agreements"). So #200's parked status is precise: it resumes only once the **parameter + potential search** (now the live work on [`10a`](m5_10a_neutrino_oscillations.md) / #236) yields a `(g, delta, potential)` that agrees with data , the charged-lepton case (point-like) reuses whatever that search finds. Duda also re-flagged the `eta` placement in `m5_9_4`: it was a **docstring** bug (the docstring wrote `diag(g,1,delta,0)` while the engine orders `diag(1,delta,0,g)` with g + the Minkowski minus both at index 3 , physics correct, notation misleading). **Fixed** in `m5_9_4`'s docstring (the CONVENTION note); the new neutrino build adopts Duda's `diag(g,1,delta,0)` + `eta=diag(-1,1,1,1)` index-0 convention. The neutrino mass-as-loop-length-density model is neutrino-specific; #200's mass approach (the `g`/`delta` QED-Lagrangian structure) is unchanged, still gated on the same search.

## Dr. Duda's review (verbatim, 2026-06-20, models-of-particles list)

```text
Dear Rodrigo,

I see lepton mass turn out completely incorrect ... but using 3x3 diag(1, delta, 0) and Faber's probably ansatz?

They definitely require full 4x4 field e.g. allowing these negative mass/energy contributions from gravity and oscillations.

Also Faber's ansatz is only approximation already for electron - serious calculations (e.g. https://www.mdpi.com/2218-1997/11/4/113) use it only as the starting point of energy minimization.

Generally for simulation the most difficult is regularization in the center of singularity, depending on Higgs-like potential deforming from this diag(g, 1, delta, 0) ... which details are still to be found.

Also I have no idea what parameters are used, starting from the basic g, delta ... I have briefly looked at various Python files there, but don't see them, and they are much too simple.

There are required serious simulations e.g. to convince people - if they are performed, we need the details: parameters, potential, minimizing energy configurations, etc.

Best,
Jarek
```

Reference he gives: Faber MTF regularization paper, MDPI Universe 11/4/113 (2025), <https://www.mdpi.com/2218-1997/11/4/113> (Faber's ansatz is the START of an energy minimization, not the answer).

## Evaluation of Duda's critique

| Duda's point | Assessment | Verdict |
| --- | --- | --- |
| "using 3x3 diag(1, delta, 0) ... they definitely require full 4x4 field, allowing negative mass/energy contributions from gravity and oscillations" | Correct, and the core error. Our 3x3 static field has no `g` (gravity/boost) axis and no clock oscillation, so it computed only the POSITIVE gradient + potential energy. The real rest mass is that MINUS the negative contributions: clock activation (the platform measured ~21% rest-energy reduction) and boost-dressing (GEM proportional to -(b.g)^2). A static 3x3 field structurally cannot see either. | ✅ valid, decisive |
| "Faber's ansatz is only an approximation already for the electron; serious calculations use it only as the starting point of energy minimization" | Correct. M5.9.1/M5.9.3 evaluated/scanned energy on the FIXED Faber profile `s(r)=r/sqrt(r^2+r0^2)`, never relaxing the field to its true minimum. So "E ~ Lambda^3" is a property of the fixed ansatz, not of the minimized configuration. | ✅ valid |
| "the most difficult is regularization in the center of singularity, depending on a Higgs-like potential deforming diag(g,1,delta,0) ... details still to be found" | Confirms our diagnosis (the confiner / Higgs term is the missing piece) but corrects the setting: it deforms the 4x4 `diag(g,1,delta,0)`, not the 3x3, and the exact form is genuinely unknown (an open search). | ✅ confirms + sharpens |
| "I have no idea what parameters are used (g, delta...); the Python files are much too simple" | Fair, with a clarification: our `sandbox_v9` scripts WERE simple (numpy toys, `delta=0.3` hardcoded, no `g`). But the production M5 engine is NOT simple, it is a full 4x4 Taichi field with documented parameters (see the inventory below). The fix is to run the build on the PRODUCTION engine and document every parameter, not to write more toy scripts. | ✅ valid (methodology) |

**Net.** We correctly diagnosed the missing ingredient, and Duda agrees that is the crux. But the setup was too simplified to compute the masses: wrong field dimension (3x3, no gravity/clock), wrong method (fixed ansatz, not minimization), and the toy scripts hid the real (documented) engine. The "E ~ Lambda^3" result is a toy artifact, not a physics claim. Duda is raising the bar to a full production engine simulation, not rejecting the direction.

## A unifying insight (#200 + #220)

The fix for #200 is the SAME machinery flagged as the residual for #220. Both need the dynamical 4x4 clock field: #200's masses need the clock oscillation + boost (the negative contributions); #220's clock-scaling needs the V-on dynamical readout. One full production engine build closes both. The clock and the mass are intertwined: the de Broglie clock energy is part of what sets the rest mass.

## What the production engine already provides (the inventory)

This directly answers "the Python files are much too simple": the toys were simple; the production engine is a full 4x4 Taichi stack. Verified 2026-06-20.

| Capability | Where | Notes |
| --- | --- | --- |
| 4x4 matrix field (MDIM=4) | `medium.py` (`TensorField`) | M is block-diag(`g` on index 0, spatial 3x3); `LC_G=8.0` (the boost/gravity axis), `LC_DELTA`/`BIAXIAL_DELTA` |
| Boost-dressed hedgehog seeder | `engine1_seeds.py: seed_dressed_hedgehog_M(tf, cx,cy,cz, r0_vox, rhoc_vox, delta, b_star, rw_vox, kick_theta)` | biaxial frame + Faber melt + core-localized boost mixing `e_Theta` with the time axis (the GEM dip = the negative gravity energy); writes the clock tangent `M_psi` and a clock-phase kick |
| 4D evolution (time axis live) | `engine2_pde.py: compute_curvature_flux_4d` (Minkowski eta=diag(-1,1,1,1)), `evolve_M_4d(...)` | the clock + boost evolve; drift guard + diagonal inertia `km` |
| Constrained integrator (validated) | `engine2_pde.py: flux_4d_constrained / update_P_4d / solve_constrained_4d / update_M_4d_constrained` | spectral-projection symplectic; f64-vs-f32 validated to 5e-15 |
| LdG potential | `engine2_pde.py: V_M / dV_M`; coeffs in `_launcher.py` | `V = a.Tr(M^2) - b.Tr(M^3) + c.(Tr M^2)^2` on the spatial block; `ldg_c=K.c^2/dx^4`, `ldg_a=-2.ldg_c.(1+delta^2)`, `ldg_b=0`; well minimum at `Tr(M^2)=1+delta^2`. KNOWN GAP: confines amplitude but does not fully stabilize biaxiality (needs an extra invariant, "Q7", already flagged to Duda) |
| Hamiltonian energy | `engine3_observables.py: compute_energyH_density_M` | `H = (1/2)\|\|Mdot\|\|^2 + c^2.4.sum\|\|[M_mu,M_nu]\|\|^2 + (V - v0)`; vacuum-subtracted |
| Clock-frequency tracker | `engine3_observables.py: update_trackers_M` | EMA of `\|\|Mdot\|\|_F` = the de Broglie clock rate |
| Director relaxation (gradient descent) | `engine2_pde.py: relax_director_step` + `AUTO_RELAX_STEPS` | descends the Frank energy on the director (biaxial seeds currently skip it) |
| Parameter block | `xparameters/_topo_biaxial1_von.py: TOPOLOGY_SEED` | `MODE, CENTER, R0_FRACTION, RHOC_VOXELS, BIAXIAL_DELTA, B_STAR, RW_FRACTION, CLOCK_KICK, LDG_STIFFNESS_K, AUTO_RELAX_STEPS, DT_SCALE_4D, INTEGRATOR_4D, KM_INERTIA_4D` + `UNIVERSE_EDGE, TARGET_VOXELS` |
| Headless template | `research/sandbox_v8/m5_8_1_headless_check.py` | full instantiate -> seed -> step -> read-energy path |

What must be ADDED for the lepton-mass build:

| Missing piece | What it is |
| --- | --- |
| GEM / boost-energy measurement | a kernel to extract the boost field `b(x)` from `M[a,3]` and compute `-(b.g)^2` per voxel, so the negative gravity contribution enters the rest mass |
| Clock-energy separation | decompose `(1/2)\|\|Mdot\|\|^2` into the ZBW-oscillation (kinetic) and the clock-activation floor (negative), so the rest mass reflects the ~21% reduction |
| True energy-minimization mode | currently only leapfrog dynamics + director relaxation; add gradient-flow / steepest-descent on the FULL 4x4 energy functional (Faber profile -> relax to the true minimum) |
| Higgs-potential search | refine the LdG `V` into a fixed-density confiner that (a) selects the core scale (the `r0^3` term the M5.9.3 toy showed is needed) and (b) stabilizes biaxiality (the Q7 invariant). The exact form is the open "details to be found" |
| Parameter-sweep + documentation | loop over `(g, delta, K, b*, rw, r0)`, save the minimizing configurations + the energy ledger (kinetic / curvature / V / GEM / clock), plot, so it is reproducible and convincing |

## TASK PLANNING (corrected #200 scope, draft for the issue)

### Goal

Compute the charged-lepton rest masses by ENERGY MINIMIZATION on the full 4x4 field `diag(g,1,delta,0)`, with a Higgs-like core regularization and the negative energy contributions from gravity (boost) and the clock (oscillation), documented to a standard that convinces working physicists.

### Approach (the corrected pipeline)

1. Run on the PRODUCTION 4x4 engine (`medium.TensorField` MDIM=4 + `seed_dressed_hedgehog_M` + the 4D / constrained integrator), not the numpy toy.
2. Field: `diag(g, 1, delta, 0)`, `g=LC_G` the boost/time axis, `delta=BIAXIAL_DELTA`, with the boost dressing `(b_star, rw)` that mixes `e_Theta` with the time axis (the GEM dip).
3. Energy: include the NEGATIVE contributions. Add the GEM kernel `-(b.g)^2` and the clock-activation energy separation, so the rest mass `H = (gradient + V) - (clock + boost)`. This is the load-bearing new physics the 3x3 missed.
4. Minimization: Faber profile as the INITIAL condition, then relax the full 4x4 field (gradient flow, or the constrained integrator with damping) to the true energy minimum. The minimized `H` = the rest mass.
5. Higgs potential: refine `V` into a fixed-density confiner that selects the core scale and stabilizes biaxiality (the Q7 invariant); sweep candidate forms (this is the open search).
6. The three leptons: as the three eigen-configurations of the biaxial structure; minimize each, read masses, compare ratios to `m_mu/m_e = 207`, `m_tau/m_e = 3477`. Pin the electron (`r0 = 2.213 fm -> 0.511 MeV`); mu/tau follow without re-tuning if it works.
7. Document every parameter (`g, delta, L, N`, the LdG coeffs / `K`, `b_star, rw, r0`, the integrator) + the minimizing configurations (saved field + plots) + the energy ledger.

### Definition of done

A documented, energy-minimizing 4x4 simulation that: (a) reproduces the electron rest energy by MINIMIZATION (not the fixed ansatz), (b) includes the negative boost + clock contributions in the rest mass, (c) searches the Higgs-potential form for one that selects discrete core scales, (d) reports the three lepton masses (or, honestly, how far off they are and what the potential search showed), (e) ships a full parameter block + the minimizing configurations, reproducible by a third party. Honest: matching mu/tau is NOT guaranteed in one pass; the convincing deliverable is the complete documented framework + the potential search, whether or not the spectrum lands.

### Sub-tasks (sequence)

1. Stand up the production 4x4 engine headless for one dressed hedgehog; reproduce the electron rest energy by minimization (relax, not the fixed ansatz); document all parameters. [foundation]
2. Add the negative-energy instrumentation: the GEM kernel `-(b.g)^2` + the clock-activation energy separation; verify the ~21% reduction and the boost dip enter the rest mass. [the key new physics]
3. Add a true energy-minimization mode (gradient flow / steepest descent on the full M functional). [method]
4. The Higgs-potential search: refine `V` to a fixed-density confiner selecting the core scale + stabilizing biaxiality (Q7); sweep forms. [the crux, open]
5. The three eigen-configurations: minimize each, read masses, compare to 207 / 3477. [the result]
6. Documentation + the convincing writeup (parameters, potential, minimizing configurations, energy ledger, plots). [Duda's explicit ask]

### Risks / honest unknowns

- The Higgs-potential form is genuinely unknown (Duda: "details still to be found"); the search may not converge to the spectrum in one effort.
- The biaxial vacuum needs an extra invariant (Q7) to be stable; may have to be added to `V`.
- Running the constrained 4D integrator to a true energy MINIMUM (vs stable dynamics) may need a damping/relaxation addition.
- Compute: the 4D Taichi engine at 64^3 is GPU-heavy; a minimization over a parameter sweep is a multi-session, possibly multi-hour effort.
- Absolute-scale calibration ties to the geometric-mean clock law (established as a scale-covariant principle); masses come out as ratios first, then anchored at the electron.

### Stage / tracking

OpenWave M5.9 sector, issue #200 (board: Next). Multi-session. Total invisibility upheld (public OpenWave physics only). The research body lands in `sandbox_v9/` (this plan + the new scripts + a findings doc); the issue gets the distilled writeup only after local work is solid and committed.

## Open questions for Dr. Duda (content bullets for the maintainer's reply, not a draft)

These would let us converge faster; phrased as points to raise, for the maintainer to put in his own words.

- Confirm the eigenstructure and the lepton-distinguishing knob: is it `diag(g, 1, delta, 0)` with `g` the time/boost axis (our `LC_G=8.0`), and do the three charged leptons correspond to the three spatial eigen-axes, or to a hierarchy in `g, delta`?
- The Higgs / core-regularization form: does he have a preferred LdG potential (beyond the Faber MDPI-2025 starting point) for the deformation of `diag(g,1,delta,0)`? Specifically, what invariant stabilizes the BIAXIAL vacuum (our `V` confines amplitude `Tr(M^2)=1+delta^2` but does not hold biaxiality, the "Q7" gap we already flagged)?
- The negative contributions: confirm the mass ledger is `m = (gradient + V) - (clock-activation + boost-GEM)`; is ~21% the right clock-activation magnitude, and how does the boost `b` differ per lepton?
- Parameters: what `g, delta` does he expect for the electron, and do mu/tau arise by varying `g`, `delta`, or the eigen-axis selection? (We will document whatever we run; we want his expected starting values.)
- The minimization protocol: confirm "Faber profile as the initial condition, then relax to the energy minimum" is the intended method; any preferred minimizer or convergence criterion?
- Validation order: agree we should first reproduce the ELECTRON by minimization (with the negative contributions included) as the gate, before attempting mu/tau?

---
_Local plan (research/9a_lepton_mass_planning.md; the research body lives in research/sandbox_v9/). Iterate here; publish the distilled version to #200 only after it is solid and committed._
