# M5.23.3: Dynamic on-screen demo suite (charge/ZBW/spin/2-charge/gravity/EM)

> Task **M5.23.3** (M5 / Liquid-Crystal model; renumbered from M5.13 on 2026-07-23, the rendering-series consolidation). Status: **Backlog** · Gated by: headless gates (demo-tier, non-gating) · Roadmap: [`m5_roadmap.md`](../m5_roadmap.md)

This doc is the task's full record: planning + findings + future planning + documentation.

---

## Current detail (from the roadmap, 2026-07-02 migration)

The breathing time-crystal on screen: which would now SELF-START at kick = 0 per 2g. Demo value only; no UI troubleshooting on any critical path

---

## GitHub issue archive (#213)

> Migrated from OpenWave GitHub issue [#213](https://github.com/openwave-labs/openwave/issues/213) on 2026-07-02 (M5 tracking moving to this local roadmap). Title: "M5 dynamic demos (on-screen rendering): charge stable, ZBW clock, spin, 2-charges, gravity, EM field & waves". State at migration: OPEN. Labels: none.

### Issue body

## Goal

A suite of **on-screen / rendered M5 dynamic demos**: charge stable at evolve PDE, the ZBW lock + breathing time-crystal, spin/magnetism, two-charge dynamics, gravity, EM waves. Communication / demo value, built on the existing rendering stack. **Policy (NG-6, headless-first):** rendering gates nothing scientific (the physics is validated headless); these are for *seeing* and explaining the validated results, not for validating them.

## The rendering map: the time axis vs its shadows

M5.8 is a 4D tensor field (time = the 4th axis, via boosts). You cannot render the time *axis* as static geometry; you render its time-**evolution** (animation) or its 3D **shadows** (imprints). This sets what each demo is:

| Phenomenon | What it is in M5 | Visualization |
| --- | --- | --- |
| **Electric / charge** | `E = ∇·n̂` (director splay), the radial hedgehog diverging at the core like Coulomb charge | **pure static 3D geometry**: no time-projection (already done: `WAVE_MENU` 6 + E glyph) |
| **Magnetic / spin** | `B = ∇×n̂` (curl), zero for a static hedgehog, appears from the clock's circulation | a **shadow of the time axis**, but the spatial curl is 3D-renderable |
| **Gravity** | boost-tilt of the 4×4 time axis (GEM ∝ (b·g)²) | a **shadow of the time axis** |
| **The clock / ZBW / breathing** | the time-evolution of the defect's oscillation | **animate it**: the most direct view of the time-dynamics |

Electric/charge is the one field that is pure static 3D space geometry (the easiest, already built). Magnetism and gravity are 3D shadows of the time-axis coupling. The clock is best animated. Full plan: review [`m5_4b_rendering_features.md`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/4b_rendering_features.md) (its new "the time axis is not renderable, only its shadows" section + Parts 3-5).

## Demo checklist (maps to the roadmap NG backlog + `4b`)

| Demo | Home / status | Notes |
| --- | --- | --- |
| Charge stable @ evolvePDE / E-field | `4b` WAVE_MENU 6 + E glyph: ✅ mostly built | the static-3D charge geometry |
| **Clock + breathing time-crystal + ZBW + particle stability** (the headline GUI demo) | **NG-6** | the breathing time-crystal on screen, self-starts at kick = 0 (per M5.8.2g spontaneity) |
| Magnetic dipole / spin | **NG-8** | real circulating `B = ∇×n̂` from the live clock (physics half in EID-C); delete the VIZ.4 placeholder scaffolding |
| 2 charges, attraction / repulsion | **NG-9** | two-defect demo; needs the multi-center seeder + the charge-class (`π₁ = Q₈`) handling |
| Gravity | `4b § 4.7` (M5.8-4D, future) | render the GEM / gravitational-field shadow |
| EM field & waves | the outgoing tilt-wave (EM channels exist) | the divergence/curl decomposition propagating |
| Annihilation | **merged → #198** | the capture → breather → vacuum dynamics; rendering is a deliverable of #198, not duplicated here |

## Review (the rendering plan)

Review [`openwave/xperiments/m5_liquid_crystal/research/4b_rendering_features.md`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/4b_rendering_features.md): the full feature catalog, the observable→channel map, the glyph/flux-mesh assignments, the gauge-stable charge fix, the magnetic-dipole sample, the gravitational-field (future) section, and the M5.8 4×4 safety contract. The new "time axis vs its shadows" subsection (Part 3) frames which demos are static-3D / animation / shadow.

## Definition of done

- The clock / breathing / ZBW / particle-stability GUI demo renders on screen (NG-6).
- Charge (static 3D), magnetism/spin (NG-8), two-charge (NG-9), gravity (4b §4.7), EM-wave demos each render their observable per the map above.
- Placeholder scaffolding (VIZ.4 hardcoded moment, etc.) removed where the live physics now feeds the channel.
- A short `4b` "as-built" log entry per demo as it lands.

Model: M5. Demo-only (gates nothing, headless-first). Annihilation rendering tracked in #198.

### Issue comments

_No issue comments._
