# M2 Free Wave: Model Briefing

> **What M2 brings.** A stable GPU voxel-grid PDE wave substrate: an energy-wave freely
> propagates, reflects, and superposes on a 3D lattice (longitudinal + transverse), with no
> particles present. Like M1, M2 is a **wave-physics library to mine, not a
> particle-emergence model**, and it is not scored in the
> [`MODELS.md`](../../../MODELS.md) coverage matrix.

![M2_hero](../../../images/demo1.gif)

## Identity

| Field | Value |
| --- | --- |
| Model ID | M2 |
| Name | Free Wave |
| Author | OpenWave wave-dynamics engine (no single originating theorist named) |
| Physics anchor | Energy Wave Theory constants (Jeff Yee lineage); Wolff standing-wave design |
| Role | PDE wave substrate + primitives for the particle models, not a coverage column |
| Primary sources | `research/00_LEVEL1_PLAN.md`, `research/10_STABILITY_ANALYSIS.md`, `research/stable_waves.md` |
| In-repo | `medium.py` (voxel grid, `psiL`/`psiT`) + `wave_engine.py` + `_launcher.py` |

## Model Profile (what it brings, short form)

| Attribute | M2 |
| --- | --- |
| Substrate | cell-centered cubic voxel grid (~1e9 target voxels, f32); scalar displacement split into longitudinal `psiL` + transverse `psiT` |
| Free wave | an energy-wave propagating / reflecting / superposing with no wave-centers present |
| Dynamics | 3D wave eqn `∂²ψ/∂t² = c²∇²ψ`, leap-frog / Verlet, 6-pt (opt. 18-pt) Laplacian, Dirichlet reflective walls |
| Stability | escapes the LEVEL-0 numerical explosion via slow-motion on wave *speed* (not timestep); CFL `dt ≤ dx/(c√3)`, single step per frame |
| Charge / particle | "charge" = energy-charging the field (pulse injection), NOT electric charge; no particles |
| Free parameters | universe size / edge, target voxels (resolution), static boost (charge amplitude), sim speed; physics constants fixed |
| Anchor | EWT energy-wave constants (λ ≈ 28.5 am, f ≈ 1.05e25 Hz, ρ ≈ 3.86e22 kg/m³); analytic, no lab data |

## Wave Primitives Provided

M2 does not model particles, so in place of a per-particle field-configuration table it
supplies the PDE wave mechanics the higher models reuse:

| Primitive | What it is |
| --- | --- |
| charging (energy injection) | `charge_full` (radial outgoing), `charge_gaussian` (packet), `charge_oscillator_*` (continuous drivers) |
| propagation | `propagate_wave` (leap-frog PDE, L + T), `compute_laplacian` (6 / 18-pt) |
| damping / envelopes | `damp_full`, charge-envelope + damping-factor helpers |
| wave-center interactions | `interact_wc_*` (standing, spin up / down, lens, drain), experimental, NOT wired into the live loop |

Wave-centers and particle formation are Phases C-D of the plan (designed + experimental
code), not validated; the live loop today is charge + free propagation only.

## Implementation Status

Not a MODELS.md coverage column. M2 is self-described as NOT YET VALIDATED (architecture
complete, implementation in progress). What runs live (`python _launcher.py`): grid build +
energy charging + free-wave PDE propagation with reflective boundaries and flux-mesh
visualization; a 20,000-step stable-charging verification is on record. Wave-center /
particle dynamics are not running. The M7 plan lists M2 among the wave library (M1-M4) to
mine, "not a parent" model.

## Roadmap

| Next | What lands |
| --- | --- |
| Phase A / B (core + charging) | solver + energy-conservation + boundary tests; measure wave speed and λ |
| Phase C / D (wave centers + particles) | reflective voxels → standing waves → mass from trapped energy → force + Newtonian motion |
| Phase E / F (viz + analysis) | data export, plots, video |
| Validation targets (open) | wave speed ≈ c (±5-10%), energy drift < 1% over 1M steps, standing waves at `r = nλ/2` |

## Help Wanted

M2 is a shared wave library, not a scored model, so there are no coverage cells to fill.
Useful contributions:

| You can contribute | How |
| --- | --- |
| A new wave primitive or demo | a new source / boundary / driver, as an `xparameters/` preset |
| Wire the wave-center engine into the live loop | move the experimental `interact_wc_*` kernels from designed to running |
| Reuse into the particle models | use the PDE substrate as a seeder for M7 (the M7 plan calls M1-M4 the library to mine) |

Flow: open an issue or discussion → fork → branch → PR with a DCO sign-off
(`git commit -s`), under Apache 2.0. See [`../../../CONTRIBUTING.md`](../../../CONTRIBUTING.md)
and [`../../../MODELS.md`](../../../MODELS.md). Model discussion runs in the
Models-of-Particles group.

## Rich Context for Deep Reader

This is a top level documentation and orientation content. For additional context on this model, a detailed read in the /theory and /research folders is recommended, as well as the production files in this model root folder, that may contain the canonical full PDE implementation of the theory.
