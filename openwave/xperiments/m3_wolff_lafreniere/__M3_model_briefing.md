# M3 Wolff-LaFreniere (WSM): Model Briefing

> **What M3 brings.** The Wave Structure of Matter made runnable: particles as standing-wave
> centers in an elastic wave medium, where force is an energy gradient (`F = −∇E`). It is the
> scalar engine that holds OpenWave's EWT validation record, demonstrating the standing-wave
> lock-in *mechanism* while honestly documenting where selectivity fails. It shares the
> "EWT (M4)" [`MODELS.md`](../../../MODELS.md) column with M4.

![M3_hero](../../../images/demo3.gif)

## Identity

| Field | Value |
| --- | --- |
| Model ID | M3 |
| Name | Wolff-LaFreniere (Wave Structure of Matter / WSM) |
| Author | Milo Wolff (Space Resonance Theory) + Gabriel LaFreniere (wave mechanics); formalized as Energy Wave Theory by Jeff Yee (with Dieter Hauger) |
| Extension | Łukasz Smoliński (BCC-lattice, zero-parameter EWT) |
| Relationship | the scalar WSM engine; the vector-PDE successor is M4 (see [`../m4_ewt/__M4_model_briefing.md`](../m4_ewt/__M4_model_briefing.md)) |
| Primary sources | Wolff, *Exploring the Physics of the Unknown Universe* (Ch. 12); LaFreniere wave-mechanics pages; Jeff Yee EWT papers |
| In-repo | `medium.py` + `particle.py` + `wave_engine.py` + `force_motion.py` + `_launcher.py`; `research/0_OVERVIEW.md`, `0_STATUS.md`, `0_ROADMAP.md`, `0a_equations.md` |

## Model Profile (what it brings, short form)

| Attribute | M3 |
| --- | --- |
| Substrate | elastic fluid-like spacetime medium; scalar longitudinal displacement field ψ |
| Particle | a wave-center (WC): a point that elastically disturbs the base wave; the standing-wave structure IS the particle |
| Charge | intended emergent (tetrahedral geometry + spin, L→T); currently imposed via `cos(source_offset)` (0 = electron, π = positron) |
| Stability | standing-wave lock-in: same-phase WCs settle into sinc energy wells at λ spacing |
| Governing eqn | Combined Wolff-LaFreniere `ψ(r,t) = A·[sin(ωt−kr) − sin(ωt)]/r`, 1/r energy-conserving falloff |
| Force | `F = −∇E` (energy gradient) |
| Free parameters | EWT wave constants (amplitude A, wavelength λ, density ρ) + per-script envelope / threshold choices |
| Anchor | ρ = 3.86e22 kg/m³, c, A₀ = 9.22e-19 m, f₀ = 1.05e25 Hz, λ₀ = 2.85e-17 m |

## Field Configuration of Particles

Standing demand of any particle model: *state the field configuration of each particle, and
say whether it uses topological vortices.* M3's answers (all "no topology": it is a scalar field, so there
is no winding number, no vortex, unlike M5's defects or M7's toroidal knots):

| Particle | Field configuration in M3 | Topological vortex? |
| --- | --- | --- |
| Electron neutrino | 1 wave-center (K = 1), symmetric, neutral | ❌ standing wave |
| Electron / positron | K = 10 "1-3-6 tetrahedron" of in-phase WCs (charged); positron = opposite nodes (λ/2 offset) | ❌ standing wave |
| Annihilation pair | two opposite-phase WCs (0 vs π), head-on → positronium wells | ❌ |
| Heavier families | neutral K = 1, 8, 20 (ν); charged K = 10, 28, 50 (e, μ, τ); magic numbers 2, 8, 10, 20, 28, 50 | ❌ |

## Implementation Status

The EWT family shares one MODELS.md column, "EWT (M4)": **0 ✅, 8 ⚠️, 3 ❌, 10 🚧** (of 21).
These are the results recorded on this scalar engine; the vector-PDE successor is M4.

| Sector | Status |
| --- | --- |
| Standing-wave lock-in | ⚠️ same-phase wells at λ spacing (mechanism shown) |
| Annihilation | ⚠️ opposite-phase pair annihilates with documented assists (0.5λ threshold, damping, velocity clamp) |
| Electron mass | ⚠️ lock-in shown; mass values from analytic EWT, not in-sim dynamics |
| Charge quantization | ❌ charge sign imposed via `cos(source_offset)`, not emergent (#203) |
| Lepton spectrum (K-selectivity) | ❌ all K = 2..10 equally stable at perfect placement; K = 10 breaks worst under perturbation (#201) |
| Electric / Coulomb force | ❌ sinc-envelope barriers block far-field attraction / repulsion (#202) |
| de Broglie clock · spin · quarks · gravity | 🚧 not yet |

Honest conclusion (`0_STATUS.md`): the scalar monochromatic model shows the *mechanism*
(lock-in) but not the *selectivity* (why K = 10, not K = 3). It needs missing physics
(variable λ, spin, non-linear terms). Deep dives:
[`research/0_OVERVIEW.md`](research/0_OVERVIEW.md),
[`research/0_STATUS.md`](research/0_STATUS.md) (honest blockers),
[`research/0_ROADMAP.md`](research/0_ROADMAP.md).

## Roadmap

| Next | What lands |
| --- | --- |
| Magic-number stability | stabilize K = 2, 8, 10 (the electron tetrahedron) under perturbation |
| Variable λ(r) | Yee-Hauger wavelength shells to fix the non-node pair distances in K = 10 |
| Non-linear ψ³ | Smoliński soliton / Lagrangian stabilizer |
| Flux-based Coulomb | `S = −c²·ψ·∇ψ` as the far-field force candidate; 3D flux integration |
| Vector successor | the L→T spin + far-field EM route moves to M4 (the vector engine) |

## Help Wanted

| You can contribute | How |
| --- | --- |
| A K-selectivity mechanism | make K = 10 uniquely stable under perturbation (variable λ, ψ³, volume-integrated ∇E) |
| A far-field Coulomb test | a flux-based force that survives the sinc envelope |
| An independent reproduction | re-run the lock-in / annihilation configs and confirm pass / fail |
| A documented negative | a runnable "this does not work, here is why" is as valuable as a positive |

Flow: open an issue or discussion → fork → branch → PR with a DCO sign-off
(`git commit -s`), under Apache 2.0. See [`../../../MODELS.md`](../../../MODELS.md)
§ Contributing, [`../../../ONBOARDING_MODELS.md`](../../../ONBOARDING_MODELS.md),
[`../../../CONTRIBUTING.md`](../../../CONTRIBUTING.md). Model discussion runs in the
Models-of-Particles group.

## Rich Context for Deep Reader

This is a top level documentation and orientation content. For additional context on this model, a detailed read in the /theory and /research folders is recommended, as well as the production files in this model root folder, that may contain the canonical full PDE implementation of the theory.
