# M5.23.7: WM6 gauge-stable signed charge view (winding density)

> Task **M5.23.7** (M5 / Liquid-Crystal model; renumbered from M5.6.5b on 2026-07-23, the rendering-series consolidation). Status: **Backlog** · Gated by: - · Roadmap: [`m5_roadmap.md`](../m5_roadmap.md). Measurement-side input since staging: the production guarded winding instrument (the M5.24 round-2 § 5.4 port); this task is the RENDER view of the signed charge.

This doc is the task's full record: planning + findings + future planning + documentation.

---

## Current detail (from the roadmap, 2026-07-02 migration)

Via the topological winding density (Brouwer degree), not eigenvector-sign band-aids; matters when two-defect ± charge becomes load-bearing

- [ ] **⤺ M5.23.7: gauge-stable SIGNED charge view (the residual; the glyph half is DONE in VIZ.1).** Context: the apolar `n̂ ≡ −n̂` sign-flip corrupted *two* render targets. **(a) Director glyph, ✅ FIXED in VIZ.1 (2026-05-30):** centered + barbless render is gauge-stable (`n̂→−n̂` swaps endpoints = identical), no more 180° slosh. **(b) Charge field WM6 (`∇·n̂`), STILL the open item:** its *local sign* flips under Evolve-PDE (`n̂→−n̂ ⇒ ∇·n̂→−∇·n̂` → spurious +/− patches). **DECISION (Rodrigo 2026-05-30): leave WM6 honest-but-flipping for now**, the charge *magnitude/location* is physical, and the charge-region **expansion is REAL physics (free-defect orientation dispersal, M5.7.2 / M5.6.5c, `V` confines amplitude not orientation), NOT the flip.** The `|∇·n̂|`-unsigned mode was tried + dropped (redundant with WM6, UI clutter). **When to address: M5.8** (sustained dynamic runs) or two-defect work where reliable ± charge *between* defects is load-bearing, via the gauge-invariant **topological winding density** (Brouwer degree, M5.1 `compute_winding_number`, which never flips), NOT band-aiding the eigenvector sign. Cheaper interim if ever needed sooner: per-step defect-relative re-pin (`n̂·r̂_defect > 0`, the `relax_director_step` `pin_signs` logic, currently seed-only). Detail: [`4b §4.4`](m5_4b_rendering_features.md) + [`4b §5.2`](m5_4b_rendering_features.md) #1b.
