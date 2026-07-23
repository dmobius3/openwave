# M5.14: 3D dynamical pair-annihilation (capture to breather to vacuum)

> Task **M5.14** (M5 / Liquid-Crystal model). Status: ❌ **RETIRED 2026-07-23** (user decision, the series reorganization review) · Roadmap: [`m5_roadmap.md`](../m5_roadmap.md)
>
> **Why retired**: the scope has been delivered at higher rigor by the verified-L era. The capture film + the measured mechanism (CONDUIT ANNIHILATION through the connecting axial line, charge ±1 → ≤ 0.005, energy ledger 0.1%) landed at [M5.21.4](m5_21_4_task_details.md) (2026-07-21); the energy-release-as-outgoing-waves read landed at the [M5.21.6](m5_21_6_task_details.md) dynamics arm (23% of the start energy radiated, ledger closed to 3 decimals); the localized oscillating remnant (breather candidate) was measured at [M5.20](m5_20_task_details.md). The 1+1D Sine-Gordon `m5_14_*` scripts/data stay the frozen validated record; GitHub #198 stays archived below.

This doc is the task's full record: planning + findings + future planning + documentation.

---

## Current detail (from the roadmap, 2026-07-02 migration)

_No prior roadmap row; scoped entirely by GitHub #198 (migrated below)._

---

## GitHub issue archive (#198)

> Migrated from OpenWave GitHub issue [#198](https://github.com/openwave-labs/openwave/issues/198) on 2026-07-02 (M5 tracking moving to this local roadmap). Title: "3D dynamical pair-annihilation: capture to breather to vacuum evolution". State at migration: OPEN. Labels: none.

### Issue body

## Goal

Run the **full 3D dynamical pair-annihilation** of a ± hedgehog pair in the M5 field, the capture → breather → vacuum-decay evolution, and confirm that topology cancellation releases the stored field energy as outgoing waves.

## Why now

On the antimatter / annihilation row the **accounting is already settled**, but it is accounting, not dynamics:

| Piece | Status |
| --- | --- |
| Charge ledger | ✅ single ±1 → Q ≈ ±1; the enclosing sphere of a ± pair → Q = 0 (additive) (`m5_8_2v_pair_annihilation_budget.py`) |
| Energy ledger | ✅ pair rest energy ≈ 2× H_static, balances |
| 1+1D principle trail | ✅ SG kink + antikink → breather → vacuum, Q = 0 throughout (`m5_14_sine_gordon_annihilation.py`) |
| **3D dynamical capture → breather → vacuum** | 🚧 open, this issue |

`m5_8_2v` is a budget ledger (it sums charges and energies of a placed configuration); it does not evolve the pair. The missing piece is an EOM-driven 3D run where the two defects actually approach, capture, ring down through a breather, and decay to vacuum.

## Scope

- Initialize a +1 / −1 hedgehog pair at separation `d` in the M5 matrix field.
- Evolve under the M5 equation of motion (the same integrator as the molten-clock runs).
- Track over time: topological charge `Q(t)` (should stay ≈ 0 for the enclosing volume and resolve the two unit charges into the cancellation), total energy `H(t)` partitioned into bound vs radiated, and the outgoing-wave flux through a far sphere.
- Confirm the released energy ≈ the pair rest energy (2× H_static) leaves as radiation, and that the late state is vacuum (or a decaying breather), matching the 1+1D trail in 3D.

## Acceptance

A single headless script + plot showing `Q(t)`, the bound/radiated energy split, and the far-sphere flux for at least one separation, demonstrating capture → breather → vacuum with the energy ledger closing dynamically.

Physics-only, headless. Builds directly on `m5_8_2v` and `m5_14_sine_gordon_annihilation.py`.


---

## RENDERING DELIVERABLE (the annihilation demo, 2026-06-17)

The **annihilation dynamic demo** lives here (not in the M5-demos umbrella): render the **capture → breather → vacuum** evolution of a ± pair on screen, with the charge ledger (Q: ±1 → 0) and the energy ledger (pair rest energy released as outgoing waves) shown live. Demo-only / headless-first (the physics is validated headless; the render is for *seeing* it). Use the existing rendering stack (`WAVE_MENU` charge + energy channels, the E/B glyphs); see [`m5_4b_rendering_features.md`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m5_liquid_crystal/research/4b_rendering_features.md). A short `4b` as-built note when it lands.


---

## Antimatter + annihilation — detail moved from the MODELS.md coverage cell (2026-06-19)

The MODELS.md "Antimatter + annihilation" cell (M5 column) was condensed to a summary; its prior full descriptive text is preserved here:

> Q → −Q under reflection; annihilation = topology-cancellation releasing stored field energy as outgoing waves (natural for topological soliton models, starting with sine-Gordon). Gauss law + (topological, Gauss-Bonnet) charge quantization permit charge destruction ONLY by annihilation with the opposite charge, releasing the integrated-Hamiltonian energy (the particle's mass); stopping the clock cannot destroy charge. Principle-level trail run: SG kink + antikink pass through in the integrable limit, and with radiation losses capture → breather → vacuum decay, Q = 0 throughout. On the production matrix the accounting is settled: the charge ledger is validated (single ±1 → Q ≈ ±1, the enclosing sphere of a ± pair → Q = 0, additive) and the energy ledger balances (pair rest energy ≈ 2× H_static, released as outgoing waves); the remaining piece is the full 3D dynamical capture → breather → vacuum evolution.

### Issue comments

_No issue comments._
