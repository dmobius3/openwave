# M5.7.3: Close Eq.23 spin-density cross-check (optional)

> Task **M5.7.3** (M5 / Liquid-Crystal model). Status: ❌ **RETIRED 2026-07-23** (user decision, the series reorganization review) · Roadmap: [`m5_roadmap.md`](../m5_roadmap.md)
>
> **Why retired**: it would test Close Eq.23 as an alternative to Duda Eq.18, but the dynamics of record has moved past both to the M5.18-verified purely-quartic Lagrangian (canonical stack, [M5.24](m5_24_task_details.md) production port); the M5.7 question it cross-checked is already covered by the 4D-clock explanation. Self-labeled "nice-to-have" since staging.

This doc is the task's full record: planning + findings + future planning + documentation.

---

## Current detail (from the roadmap, 2026-07-02 migration)

The spin-density-wave alternative; nice-to-have — Duda's 4D-clock explanation already covers the M5.7 results

- [ ] **⤺ OPTIONAL cross-check — Close's Eq.23 (Plan B residual).** Our M5.7 nulls are for **Duda's Eq.18**; **Close's Eq.23** (spin-density wave, `∇·s=0`) is a different equation. Read it from his published paper (`../theory/Equation-of-Everything.pdf`), implement + test the 3 candidate forms (the physical one preserves `∇·s=0` + bounded energy). *Nice-to-have* — Duda's 4D-clock explanation already covers the M5.7 results, so this is not a blocker.
