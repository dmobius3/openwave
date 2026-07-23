# M5.23.6: Gravity / time-dilation viz + electron gravitational mass

> Task **M5.23.6** (M5 / Liquid-Crystal model; renumbered from M5.8.7 on 2026-07-23, the rendering-series consolidation). Status: **Backlog** · Gated by: M5.8 4D g-axis live · Roadmap: [`m5_roadmap.md`](../m5_roadmap.md)

This doc is the task's full record: planning + findings + future planning + documentation.

---

## Current detail (from the roadmap, 2026-07-02 migration)

- [ ] **M5.23.6: Gravity / time-dilation viz suite (REVISIT when the 4D `g`-axis is live; breadth, NOT headline-gating)**, once the time axis exists, build the renders specced in [`4b §4.7`](m5_4b_rendering_features.md): **(a)** `g(x)` scalar flux-mesh = the **gravity well** (sequential single-sign palette, deep = clock-slow); **(b)** `∇g` **PURPLE** vector glyph = the **gravitational pull** (monopole `1/r²`); **(c)** per-voxel **clock-rate `ω`** colormap / animation-speed = **time dilation** (the wristwatch-per-voxel). The "gravity = a *bend in time*" intuition reads as the *pattern* of clock-rates across voxels (same "read the arrangement, not one glyph" lesson as bend in the director glyphs). Physics-intuition home = `0c` Lesson 3 (gravity = time-tilt gradient + Wheeler's curved-time). ⚠️ Needs the boost-`g` axis (M5.8 promotion); the gravity sector is the framework's least-developed (design expectation, not yet verified). **Seed from 2a/2b (2026-06-05): the static GEM effect is now MEASURED**, boost dressing alone lowers the defect's energy (`A(b)`: 5.97→0.39→rise; the dressed defect IS the family's static ground state) ⇒ the gravity-well `g(x)` viz has a concrete, already-validated target: render the dressing profile `b·w(r)` + its energy dip, not just a design expectation.

---

## GitHub issue archive (#209)

> Migrated from OpenWave GitHub issue [#209](https://github.com/openwave-labs/openwave/issues/209) on 2026-07-02 (M5 tracking moving to this local roadmap). Title: "Gravity sector: electron gravitational mass from the boost/GEM coupling". State at migration: OPEN. Labels: help wanted.

### Issue body

## The target (raised by J. Duda, Models-of-Particles, 2026-06)

From the **preferred boosts** of the M5 defect, derive the **gravitational mass of the electron**. In M5, gravity enters via the **boost tilt of the 4×4 time axis** (GEM ∝ (b·g)², the negative "clock-fuel" block); the de Broglie clock's energy minimum is reached by minimizing over exactly this boost dressing (the negative boost-sector kinetic term that makes E(ω) have a minimum). So the boost/GEM sector is the natural place a gravitational-mass prediction would come from.

## Why it's timely

The electron's gravitational mass is **not** experimentally well-established: the classic free-fall measurement (Witteborn & Fairbank, [PRL 19, 1049 (1967)](https://link.aps.org/doi/10.1103/PhysRevLett.19.1049)) returned ~zero, later attributed to a screening artifact. A modern measurement may now be in reach, so a concrete model prediction (even a `m_grav / m_inertial` ratio) would be valuable and falsifiable.

## Current status

Partially validated: the boost-tilt coupling is measured (GEM ∝ (b·g)², exactly zero at zero boost; [`m5_8_2q_delta_scaling.py`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/sandbox_vn/m5_8_2q_delta_scaling.py)); the dynamical metric is not implemented. See the Gravity row of [`MODELS.md`](https://github.com/openwave-labs/openwave/blob/main/MODELS.md).

## Definition of done

- A derivation of the electron gravitational mass (or the `m_grav / m_inertial` ratio) from the boost/GEM sector.
- A clear statement of the resulting prediction and how it could be tested (e.g. against a modern free-fall / interferometric measurement).
- Honest scope: what the absence of a dynamical metric does and does not allow at this stage.

Model: M5 Liquid Crystal. Gated in practice on the absolute-scale calibration (the units → physical-scale issue) for a physical-number prediction.


---

## DUDA INPUT + REFERENCES (2026-06-17)

Context (Duda's "4 types of mass" framing): **gravitational mass is experimentally confirmed only for nucleons.** For the electron the classic 1967 Witteborn-Fairbank free-fall reading was ~zero (later attributed to a shielding/screening artifact); "generally we don't know much about gravitational and de Broglie mass." So a model prediction here is genuinely open territory, a place the model can predict where experiment cannot yet measure.
- Gravitational mass of the electron, measurement-attempt slides (WAG 2015): https://indico.cern.ch/event/361413/contributions/1776296/attachments/1137816/1628821/WAG2015.pdf
- Witteborn-Fairbank (1967): https://link.aps.org/doi/10.1103/PhysRevLett.19.1049

### Issue comments

_No issue comments._
