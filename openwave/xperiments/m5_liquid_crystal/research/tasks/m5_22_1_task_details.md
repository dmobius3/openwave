# M5.22.1: nuclei at scale (the far arc): larger nuclei, halo nuclei, alpha clusters, toward the table of nuclides

**Status**: 🚧 PLANNED STUB, far-gated. **Split out of [M5.22](m5_22_task_details.md) on 2026-07-27** (user decision): M5.22 keeps everything that runs on the hardware this project actually has (the toy 3×3 vortex-knot census), and this task carries everything that does not. The split exists so the resource question is stated honestly instead of buried inside a single task that is half-runnable: what is promised to the author is the toy ladder, and what needs compute we do not have is named here and priced on the [README § Resource Contributors](../../../../../README.md#resource-contributors) page. Records written before the split call this scope "M5.22 Phase 2". **Sequencing CONFIRMED by the author 2026-07-27** ([`m5_22_convo.md § 2026-07-27 afternoon`](m5_22_convo.md)): focus on the small nuclei first; the larger ones need separately prepared large simulations, to be thought about after the small ones come out right.

## Scope: the rungs that do not fit on current hardware

The ladder is the author's ([`m5_22_convo.md § 2026-07-22`](m5_22_convo.md), extended 2026-07-26); M5.22 runs its lower rungs, this task owns the upper ones.

| Rung | Target | Read / anchor |
| --- | --- | --- |
| A | Larger nuclei, the toy "periodic table" | qualitative agreement on electromagnetic moments across a range of knot sizes; the author's public framing question is exactly this ("various size knots, and various size nuclei - is there correspondence between them?") |
| B | Halo nuclei | stable halo configurations as a qualitative signature |
| C | Alpha clusters | C12, O16, compared against the Skyrmion models of nuclei (the author's named baseline literature) |
| D | Strange-baryon decay | decay releasing a pion or a kaon; the meson side of the same reconnection mechanism (pion = twist/reconnection of a vortex loop, kaon = a Möbius-like twisted loop with strangeness as the twist) |
| E | The table of nuclides | the author's own far target: "trying to recreate <https://en.wikipedia.org/wiki/Table_of_nuclides> even with such toy model" ([`m5_22_convo.md § 2026-07-26 ps`](m5_22_convo.md)) |
| F | Fusion and fission of knots | the experimental LC anchor: the Nature paper on fusion and fission of particle-like chiral nematic vortex knots, named by the author (X, 2026-07-18) as the style of nuclei model to build. Locate + file into `theory/` at PLAN |

## What gates it

| Gate | Content |
| --- | --- |
| The toy results | [M5.22](m5_22_task_details.md) has to land first: the lower rungs are the instrument validation for these |
| The Lagrangian details | pinned by the e + ν programs ([M5.21.3](m5_21_3_task_details.md) → [M5.21.6](m5_21_6_task_details.md)/[M5.20.8](m5_20_8_task_details.md)) |
| The realistic-parameter bridge | [Q33](../m5_question_tracker.md#q33-detail) / [M5.21.11](m5_21_11_task_details.md): any QUANTITATIVE fusion/fission or binding claim needs it. Without the bridge this arc stays qualitative like M5.22 |
| **Compute** | The distinguishing gate, and the reason for the split. The author flags supercomputer-scale simulation for this arc, and the project has neither the machine nor the budget. This is a named cost center on the [README § Resource Contributors](../../../../../README.md#resource-contributors) page: a cluster or HPC allocation, or funded compute, opens it |
| User "go" | as always |

## The scoping sub-task (the first real work here)

Before any physics, this arc needs a plan that can be handed to whoever brings the machine: lattice sizes and memory footprint per nuclide, solver choice and parallel decomposition, wall-clock estimates per rung, checkpoint/restart strategy, and the allocation routes worth applying to (national HPC calls, university partners, the author's own institutional contacts, sponsored cloud). That scoping document is itself the deliverable that makes a resource offer actionable, and it can be written before any allocation exists.

## The skeptic gauntlet carries forward

The pre-registered honest caveats from [M5.22](m5_22_task_details.md) apply to the claim language here too, and harder, because scale invites overreach: (i) what plays SU(3) confinement and asymptotic freedom, (ii) sub-barrier fusion needs tunneling (the Gamow factor) and classical knots deform rather than tunnel, so either the quantized theory or an honest scope statement, (iii) nucleon statistics ([Q34](../m5_question_tracker.md#q34-detail)).

## Cross-links

| Doc | Why |
| --- | --- |
| [M5.22](m5_22_task_details.md) | the toy census this arc continues; the charge-quantization frame and the run ladder live there |
| [`m5_22_convo.md`](m5_22_convo.md) | the author-channel record for the whole nuclear program |
| [M5.28](m5_28_task_details.md) | the composite-particles hunt (atomic orbitals), gated on the M5.22 toy results |
| [`../m5_roadmap.md`](../m5_roadmap.md) | Backlog row + the change-log entry for the split |
