# M1 Granule Motion: Model Briefing

> **What M1 brings.** A granule-medium wave-propagation sandbox: it shows how spherical
> longitudinal waves from harmonic sources superpose into interference and standing-wave
> patterns on a discrete BCC spacetime lattice. M1 is a **wave-physics library to mine,
> not a particle-emergence model**, and it is not scored in the
> [`MODELS.md`](../../../MODELS.md) coverage matrix.

![M1_hero](../../../images/demo2.gif)

## Identity

| Field | Value |
| --- | --- |
| Model ID | M1 |
| Name | Granule Motion |
| Author | OpenWave educational library (no single originating theorist named) |
| Physics anchor | Energy Wave Theory wave constants (Jeff Yee); Wolff-style in / out standing waves |
| Role | reusable wave primitives for the particle models (M3-M7), not a coverage column |
| Primary sources | EWT papers (Jeff Yee); Milo Wolff wave-structure-of-matter (in / out waves) |
| In-repo | `medium.py` (BCC lattice) + `wave_engine.py` + `_launcher.py`; `research/README.md`, `research/spherical_wave.md`, `research/phase_shift.md` |

## Model Profile (what it brings, short form)

| Attribute | M1 |
| --- | --- |
| Substrate | BCC granule lattice (fluid-like medium, sub-attometer granules; 2 per unit cell) |
| Granule | point-mass unit vibrating harmonically about its equilibrium |
| Wave | displacement field `x(t) = x_eq + A(r)·cos(ωt + φ)·dir`, phase `φ = -kr` |
| Wave types | standing · traveling · spherical · energy (base / in / out toggles) |
| Charge / particle | none (per-source phase offset is the only particle-like handle) |
| Free parameters | fixed EWT constants (c, λ, f, A, ρ); per-demo: source count / positions / phase, base-in-out toggles |
| Anchor | analytic PSHO (phase-synchronized harmonic oscillation), `v = fλ = c` by construction; theory only, no lab data |

## Wave Primitives Provided

M1 does not model particles, so in place of a per-particle field-configuration table it
supplies the wave mechanics the higher models reuse:

| Primitive | What it is |
| --- | --- |
| base_wave | background wave summed from the 8 universe-box vertices (constant amplitude) |
| in_wave | inward wave toward a wave-center (full amplitude) |
| out_wave | outward spherical wave with 1/r amplitude falloff |
| superposition engine | per-granule sum of all sources → constructive / destructive interference |
| energy-conservation primitives | 1/r falloff, singularity clamp `r_min = 1λ`, physical cap `A ≤ r` |

Wave-center → particle formation (constructive / destructive interference assembling a
particle) is aspirational here, not implemented; that work lives in the particle models.

## Implementation Status

Not a MODELS.md coverage column. M1 is a runnable Taichi GGUI application
(`python _launcher.py`) with 9 preset demos (granule motion, 3D spherical wave, standing
wave, wave interference, spacetime vibration, medium disturbance, golden-ratio spiral, and
two multi-source force demos). Parameters are correct *by construction* (the analytic PSHO
solution), not claim-backed validations. The M7 implementation plan classifies M1 as part
of the wave-physics library to mine for insights, explicitly "not a parent" model.

## Roadmap

| Next | What lands |
| --- | --- |
| Phase manipulation | richer multi-source interference phase + a particle-interaction phase (hooks only today) |
| Wave-center → particle | assembling particles from interference (aspirational; belongs to M3-M7) |
| Doc cleanup | `research/README.md` drift (legacy names, removed presets) |

## Help Wanted

M1 is a shared wave library, not a scored model, so there are no coverage cells to fill.
The value is clean, reusable wave mechanics. Useful contributions:

| You can contribute | How |
| --- | --- |
| A new wave primitive or demo | e.g. polarization, dispersion, a new source geometry, as an `xparameters/` preset |
| Reuse into the particle models | wire these seeders into M3-M7 (the M7 plan calls M1-M4 the library to mine) |
| A doc / layout cleanup | reconcile `research/README.md` with the current file tree |

Flow: open an issue or discussion → fork → branch → PR with a DCO sign-off
(`git commit -s`), under Apache 2.0. See [`../../../CONTRIBUTING.md`](../../../CONTRIBUTING.md)
and [`../../../MODELS.md`](../../../MODELS.md). Model discussion runs in the
Models-of-Particles group.

## Rich Context for Deep Reader

This is a top level documentation and orientation content. For additional context on this model, a detailed read in the /theory and /research folders is recommended, as well as the production files in this model root folder, that may contain the canonical full PDE implementation of the theory.
