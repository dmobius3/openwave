# M5.23.4: Two biaxial hedgehogs interaction demo

> Task **M5.23.4** (M5 / Liquid-Crystal model; renumbered from M5.6.5e on 2026-07-23, the rendering-series consolidation). Status: **Backlog** · Gated by: multi-center biaxial seeder, V-on (M5.6.5c) · Roadmap: [`m5_roadmap.md`](../m5_roadmap.md)

This doc is the task's full record: planning + findings + future planning + documentation.

---

## Current detail (from the roadmap, 2026-07-02 migration)

Needs the multi-center seeder + the charge-class decision (`π₁ = Q₈` conjugacy classes: the discovery, not a parameter)

- [ ] **⤺ M5.23.4: two biaxial hedgehogs under Evolve PDE (interaction demo).** A `_biaxial2.py` xperiment seeding TWO biaxial hedgehogs interacting (attract / repel / scatter) under the live Eq.18 matrix leapfrog, the **dynamic** biaxial analog of the M5.1/M5.4 *static-relaxation* Coulomb test (`m5_3a_coulomb_visual_geometry.md`). **Prerequisites:** (i) a **multi-center** biaxial seeder (`seed_biaxial_hedgehog_M` is single-center, and biaxial **frames don't superpose like uniaxial directors**, the overlap needs frame interpolation or a relax pass); (ii) V-on (M5.6.5c). **Open question, not a foregone result:** biaxial defects are NOT simple ±1 Coulomb charges (richer quaternion-group π₁ classification). Co-implementation partner of the **5f** dipole viz below, a twisting/interacting defect generates the circulating B. **⚠️ Charge-sign convention (flagged 2026-05-30):** unlike the uniaxial seed (`_topo_uniaxial2.py` exposes per-defect `DEFECTS:[{"SIGN": ±1}]` → kernel `n̂=SIGN·r̂`, `engine1_seeds.py:439`), the current biaxial seeder `seed_biaxial_hedgehog_M` is single-center and exposes **NO charge-sign**, it hard-builds one radial frame `O=[r̂|e_Θ|e_Φ]`. The two-defect demo must therefore first **decide what "opposite charge" even means for a biaxial defect** (the order-parameter space is `SO(3)/D₂`, `π₁=Q₈`, a quaternion *class*, not a `±` bit): pick the conjugacy-class pair that gives an attract/annihilate analog, then add the second center + that class to the seeder. This is the discovery, not a parameter. (Teaching note mirrored in `0c §L4`.)
