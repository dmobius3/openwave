# M5.23.5: Magnetic-dipole viz Stage 2 (real circulating B + placeholder cleanup)

> Task **M5.23.5** (M5 / Liquid-Crystal model; renumbered from M5.6.5f on 2026-07-23, the rendering-series consolidation). Status: **Backlog** · Gated by: M5.23.4 (was M5.6.5e) · Roadmap: [`m5_roadmap.md`](../m5_roadmap.md)

This doc is the task's full record: planning + findings + future planning + documentation.

---

## Current detail (from the roadmap, 2026-07-02 migration)

Real circulating B from the live clock + delete ALL the VIZ.4 placeholder scaffolding (incl. the hardcoded `+ẑ` moment). The physics half is being absorbed by ELECTRON-ID EID-C (the fixed-clock μ integration); what stays here is the VIZ cleanup

- [ ] **⤺ M5.23.5: magnetic-dipole viz, STAGE 2: real circulating B + PLACEHOLDER CLEANUP.** **Stage 1 ✅ DONE (VIZ.4, 2026-05-30):** the render path (bluered N/S via radial `(∇×n̂)·r̂`, B glyphs, YELLOW moment glyph) is built + validated against an analytic dipole placeholder. Stage 2, **at M5.8 when the twisting/spinning clock generates a REAL circulating B** (its magnetic moment appears naturally): point the same render at the live `∇×n̂` and **delete ALL the VIZ.4 placeholder scaffolding**, (1) `engine3_observables.fill_dipole_sample_B` kernel; (2) the `_viz_sample_dipole.py` xparameter; (3) the `DIPOLE_SAMPLE`/`DIPOLE_AXIS`/`DIPOLE_CENTER`/`DIPOLE_R0_VOX` state + loader keys; (4) **NEW FEATURE, real magnetic-moment viz:** the moment is currently **HARDCODED** (`update_moment_glyph` draws a YELLOW arrow along `m̂ = DIPOLE_AXIS = +ẑ`, a constant in `_viz_sample_dipole`, NOT computed). Implement the real thing: **compute μ from the field** (`m̂ ∝ ∫∇×n̂` over the defect, or the clock's spin axis), render the YELLOW moment glyph from THAT, and **delete the hardcoded `+ẑ`** placeholder; (5) wire `curl_radial`/`curl_center` (now dipole-only) to the **real defect center + auto-axis** (derive `m̂` from the net circulation direction). Keep `_curl_signed_proj` + the radial/axial coloring + the `moment_glyph_*` buffers, those are the production render path. Co-implementation partner of the **M5.23.4 two-defect demo** (below, a twisting/interacting defect generates the circulation). Physics link: validates the "magnetism = T-component of the defect's outgoing wave" picture (Phase-4 EM-drive leg). Moved here from M5.6 (Rodrigo 2026-05-27). Viz Note: be able to "see" magnetic field components in the director vectors render: BEND (sequence of directors) + FRANK_TWIST (director sticking out of sample place).
