# M5.8.9: G-2c-2 hardening (stiff stepper + engine quartic port)

> Task **M5.8.9** (M5 / Liquid-Crystal model). Status: ❌ **RETIRED 2026-07-23** (user decision, the series reorganization review) · Roadmap: [`m5_roadmap.md`](../m5_roadmap.md)
>
> **Why retired**: its scope (the `u + βu²` quartic flux port + the deep-floor stiff stepper) is tied to the M5.8-era quartic, which was retired with that era ([M5.21.1e](m5_21_1e_task_details.md): the one 1D transplant, removed from the canonical stack). The production need it anticipated was delivered by [M5.24](m5_24_task_details.md) (canonical verified-L stack ported, dimension-split certified dt, FIRE relaxer, sponge).

This doc is the task's full record: planning + findings + future planning + documentation.

---

## Current detail (from the roadmap, 2026-07-02 migration)

Stiffness-aware stepper (the deep-floor cascade), horizons ≫45 periods, the engine-side quartic port (the `(1+2β·u)` flux prefactor) when production-scale runs need it
