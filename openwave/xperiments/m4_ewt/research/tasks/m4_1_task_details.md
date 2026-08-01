# M4.1: Golden-angle K-selectivity + perturbation-robust stability

> Task **M4.1** (M4 / EWT model). Status: **BACKLOG**, no owner · Roadmap:
> [`../m4_roadmap.md`](../m4_roadmap.md) · Migrated from GitHub issue
> [#201](https://github.com/openwave-labs/openwave/issues/201) on 2026-08-01, when
> [T5](../../../../../dev_docs/tasks/t5_task_details.md) settled that tasks live in roadmaps and
> issues are reserved for platform defects. The issue body is archived below in full, so closing
> it loses nothing.

This doc is the task's full record: planning, then findings when someone runs it.

## PLANNING

### The two open problems

They are the standing blockers on the [`0_STATUS.md`](../../../m3_wolff_lafreniere/research/0_STATUS.md)
sheet the M3 and M4 columns share, and they are coupled: the same suspected missing physics
(non-linearity, spin, variable `λ(r)`) sits behind both.

| # | Problem | State |
| --- | --- | --- |
| 1 | K-selectivity | `K = 2..10` are equally stable at perfect placement, and `K = 10` breaks worst under perturbation. The energy landscape does not discriminate the electron configuration |
| 2 | Perturbation-robust stability | Combined Wolff-LaFreniere has shallow equilibria with no fine structure to prefer `K = 10` |

Root cause per the proposal: Combined W-L scales energy linearly in the number of wave centers,
which makes the landscape flat and every `K` degenerate.

### The proposal to be evaluated

Two interconnected geometric mechanisms, from Łukasz Smoliński (contributor). The claims are the
proposal, not a platform result:

| Mechanism | Claim |
| --- | --- |
| `r⁵` vs `r³` energy density | A soliton's energy scales as `E ∝ r⁵` while the compensating volume scales as `r³`; the `r²` disparity opens a quantized gap that isolates `K = 10`, and beyond it forces nested shells (the muon and tau generations, the "Onion Model") |
| Golden-angle self-organization | Wave centers self-organize by spherical phyllotaxis (≈ 137.5°) rather than a regular polyhedron; `K = 10` is the smallest `K` at which the pattern closes on the sphere. Spin matters: gyroscopic rigidity plus the magnetic deficit `ε_M = 1/(8π⁷)` and the lattice factor `g_v` carve the deep narrow well |

The contributor's later formalization of the non-linear stabilization terms (Gaussian density
profiles, quintic saturation, modular profiles, shell structures) is in-repo at
[`../M4_k_selectivity_Formalization.md`](../M4_k_selectivity_Formalization.md) and is the current
spec to implement against.

### The falsifiable test

Initialize `K = 2..12` wave centers in golden-angle configurations **with an initial angular
momentum**, and measure whether only `K = 10` forms a stable spherical standing wave that survives
perturbation. A first code contribution is merged:
[PR #205](https://github.com/openwave-labs/openwave/pull/205), the golden-angle / spherical
phyllotaxis test for both the M3 and M4 engines. The substrate it needs is in place:
[`../M4_engine_upgrade.md`](../M4_engine_upgrade.md) P0-P4 delivered the vector-field non-linear
solver.

### Pass/fail, as the proposal states it

Inherited from the issue and to be pre-registered in full by whoever runs it, before any numbers
exist:

| # | Question | What decides it |
| --- | --- | --- |
| 1 | Does the `r⁵` vs `r³` argument produce a deep narrow well only at `K = 10` when implemented, against today's flat linear landscape? | the measured energy landscape over `K = 2..12` |
| 2 | Does golden-angle plus angular-momentum initialization make `K = 10` uniquely survive perturbation? | perturbation survival across the same range |
| 3 | Are `ε_M = 1/(8π⁷)` and `g_v` derivable in-platform, or fitted? | an in-platform derivation, or an honest "fitted" label |
| 4 | Does the Onion Model recover the muon and tau as nested shells? | the shell energies against the lepton ratios |

A negative closes this task exactly like a positive one.

### What it would move

The [`MODELS.md`](../../../../../MODELS.md) M4 lepton-mass-spectrum cell (❌ honest negative,
K-selectivity not achieved) and, through it, the stability and magnetic-moment cells that name
`K = 10` electron stability as their prerequisite.

**Gated by**: an owner. The M4 column has no active author, and the program belongs to whoever
extends the model.

## GitHub issue archive (#201)

> Migrated from OpenWave GitHub issue
> [#201](https://github.com/openwave-labs/openwave/issues/201) on 2026-08-01. Title: "Evaluate EWT
> geometric resolution of Combined W-L K-selectivity & perturbation-robust stability". Opened
> 2026-06-15 by `xrodz`. State at migration: OPEN. Labels: `help wanted`. Body verbatim, with its
> headings demoted one level so they nest under this section.

### Summary

Proposal from **Łukasz Smoliński** (contributor; paper + code links below) for resolving the two
open problems on the EWT / Combined Wolff-LaFreniere (W-L) status sheet
([`m3_wolff_lafreniere/research/0_STATUS.md`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md)).
This issue captures the full reasoning so anyone can evaluate it openly and run the proposed test.
The claims below are the **proposal to be evaluated**, not yet validated in-platform.

A first code contribution exists: **PR #205** (golden-angle / spherical-phyllotaxis K-selectivity
test).

### The two open problems (from `0_STATUS.md`)

| # | Problem | Current state |
| --- | --- | --- |
| 1 | **K-selectivity** | All `K = 2..10` are equally stable at perfect placement; `K = 10` actually breaks worst under perturbation. The energy landscape does not discriminate `K = 10` from simpler geometries |
| 2 | **Perturbation-robust stability** | Combined W-L has shallow equilibria, no fine structure to discriminate `K = 10` from others |

Root cause per the proposal: Combined W-L assumes **linear** scaling of energy with the number of
wave centers, giving a flat landscape where all `K` are degenerate.

### Proposed EWT resolution (two interconnected geometric mechanisms)

#### Mechanism 1: the r⁵ vs r³ energy-density non-linearity

In EWT a soliton's energy scales as `E ∝ r⁵` (volumetric occupancy `r³` × amplitude `r` ×
frequency `1/r`), while the volume available for geometric compensation scales only as `V ∝ r³`.
The `r⁵/r³ = r²` disparity is a strong non-linearity: growing the soliton radius (hence `K`)
sharply raises the energy density.

- For the electron configuration the ratio `(r_e/r_ν)⁵ = 10¹⁰` creates a **quantized energy gap**
  that isolates `K = 10` from its neighbors.
- Beyond `K = 10` the `r²` term forces the system to shed excess energy into **nested shells**, the
  muon and tau generations (the "Onion Model").
- Without the `r⁵` vs `r³` imbalance, all `K` stay degenerate and minima stay shallow (exactly the
  Combined W-L behavior).
- The **magnetic deficit** `ε_M = 1/(8π⁷)` and the lattice coupling factor `g_v` add a
  `K`-dependent non-linear compression (spin-induced torque) that carves a **deep, narrow potential
  well** only at `K = 10` (and recursively via the Onion Model).

> Note from the author: in the paper `K_WC` (wave centers) is not the same as `K` (nodal metrics).

#### Mechanism 2: golden-angle (spherical phyllotaxis) self-organization

Rather than a regular polyhedral arrangement of point sources, the wave centers self-organize by
**spherical phyllotaxis**, the golden-angle distribution (~137.5°) that minimizes destructive
interference (the same principle as sunflower seed packing).

- `K = 10` is the **smallest K for which the pattern "closes" on the sphere**, forming a coherent
  standing wave.
- For `K > 10`, excess energy forces recursive shells (muon, tau) rather than a single overloaded
  core.
- **Spin matters**: the golden-angle configuration alone may not suffice. Geometric rotation of the
  whole wave-center ensemble provides gyroscopic rigidity against perturbation; the magnetic
  deficit represents the coupling between this rotation and the BCC-lattice stiffness, and the
  resulting magnetic torque compresses the soliton, deepening the well exactly at `K = 10`.

### Falsifiable test for OpenWave

Initialize `K = 2..12` wave centers in **golden-angle configurations with an initial angular
momentum**, and measure whether **only `K = 10`** forms a stable, spherical standing wave that
survives perturbation. The geometry predicts `K = 10` is the unique ground state by phyllotactic
necessity, not by assumption.

Started in **PR #205**. The author notes spin / initial angular momentum should be included to
capture the full stabilization mechanism.

### Links

| Resource | URL |
| --- | --- |
| Paper (v4.4.15) | <https://zenodo.org/records/20313808> |
| Calculation scripts | <https://zenodo.org/records/19398255> |
| Code PR (golden-angle test) | <https://github.com/openwave-labs/openwave/pull/205> |
| Open-problem source | [`m3_wolff_lafreniere/research/0_STATUS.md`](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md) |

Relevant paper sections: "The 1:100 Decadic Resonance Discovery" (near `r_e/r_ν = 100`); "The
1:10^10 Resonance as the Geometric Foundation of the Onion Model"; "Physical Origin of the r⁵
Scaling: Geometric Energy Density"; "The Recursive Lepton Hierarchy: Nodal Shell Resonance Model";
"Natural Emergence of Three Lepton Generations"; "The Geometric Unification of Lepton Properties";
"Geometric Derivation of the Neutrino Radius and the g_v Factor".

### What to evaluate

- Does the `r⁵` vs `r³` energy-density argument actually produce a deep, narrow well only at
  `K = 10` when implemented in the platform (vs the current flat, linear landscape)?
- Does the golden-angle + angular-momentum initialization make `K = 10` uniquely survive
  perturbation across `K = 2..12`?
- Are `ε_M = 1/(8π⁷)` and `g_v` reproducible / derivable in-platform, or fitted?
- Does the Onion Model recover the muon / tau generations as nested shells?

Constructive criticism explicitly welcomed by the author. Proposal raised by Łukasz Smoliński (with
Jeff Yee on the thread); captured here so the platform community can evaluate and extend it openly.

## DEVIATIONS LOG

(none)

## FINDINGS

(pending: the task has not been run)
