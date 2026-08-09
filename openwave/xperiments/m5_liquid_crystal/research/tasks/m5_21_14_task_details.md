# M5.21.14: the (1/g) dressing term (symbolic 4×4 boost hedgehog)

> Task **M5.21.14** (M5 / Liquid-Crystal model). Status: 🚧 **PLANNED STUB** · Roadmap:
> [`../m5_roadmap.md`](../m5_roadmap.md) · Staged 2026-08-09 from the author's 1:1 reply to the
> [M5.21.11](m5_21_11_task_details.md) ladder close (decode:
> [`m5_21_convo.md § 2026-08-09`](m5_21_convo.md)). Series: M5.21.x, the electron hunt
> (lepton-sector instrument work).

This doc is the task's full record: planning, then findings at the run.

## PLANNING

### Why this task exists

The [M5.21.11](m5_21_11_task_details.md) ladder closed route (b) terminally: the 4D boost-dressing
correction is O(1), branch-dependent, and FLAT in g (measured on the ladder endpoints and
retro-read in the [M5.21.8](m5_21_8_task_details.md) family record), so no 3×3-only ladder can
reach physical-regime energies or ratios. The close-out named the required shape of any successor:
the dressing carried INSIDE the ladder. The author's 2026-08-09 reply supplies the concrete first
rung, symbolic rather than numeric, so it needs no lattice, no certified branches, and no new
relaxations.

### The author's recipe (2026-08-09, pinned before any computation)

| Step | Content |
| --- | --- |
| 1 | Work symbolically in the 4×4 case with at least radial dependence, assuming the spherical boost hedgehog `MatrixExp[b(r) {x, y, z} . boostgenerators]` |
| 2 | Find the FIRST NONTRIVIAL TERM in the (1/g) expansion |
| 3 | Include that term in the further-considered 3×3 case |
| The criterion (the author's, pre-stated) | the term "needs to have negative Hamiltonian contribution to get oscillations" |
| The instrument correction it carries | the rigid uniform-m read is too crude (the measured m\* "seems much too large, and the radius dependence seems more complicated"): start from minimization of a GENERAL at-least-radius-dependent b(r), not a single rapidity knob |

### Scope sketch (to firm at PLAN if picked)

| Arm | Content |
| --- | --- |
| S1 symbolic derivation | sympy: M = Qb D Qbᵀ with Qb = exp(b(r) r̂·K) (the general-b(r) form of the rigid `qb_field` instrument in [`../scripts/m5_21_11_d_garm.py`](../scripts/m5_21_11_d_garm.py)); derive H[b] under the verified quartic L, expand in 1/g, extract the leading nontrivial term as a functional of b(r) |
| S2 the sign verdict | minimize over b(r) (variational / Euler-Lagrange); the pre-registered read: does the leading term contribute NEGATIVE Hamiltonian density at the minimizer? Either answer is the deliverable |
| S3 the 3×3 handoff | write the term as an additive correction to the 3×3 functional (the dressing carried inside any future ladder) and state what it retrodicts for the flat-gain record (the M5.21.11 g-arm + the M5.21.8 family) and for the m\*-too-large read |
| Numerical sequel (OUT of scope here) | lattice implementation of the corrected 3×3 functional plus an instrument that certifies B/C (bigger boxes or a new stencil generation): a separate task if S1-S3 land |

### What this task does NOT do

| Non-goal | Why |
| --- | --- |
| Reopen route (b) | the [framework](../findings/m5_21_11_framework.md) is frozen and its terminal verdict stands; this task builds a NEW instrument generation, it does not re-analyze the ladder data |
| Quote any mass ratio | B/C remain uncertifiable at N = 48; a symbolic term changes nothing about the certification problem |
| Wait on the potential details | the author states the current eigenvalue potential is "just a first guess" ([Q25](../m5_question_tracker.md#q25-detail)); S1 runs on the verified quartic L + the T2 base of record, and the V-dependence of the term's sign is REPORTED, not assumed away |

**Gated by**: user pick (candidate NEXT alongside [M5.22.3](m5_22_3_task_details.md) /
[M5.22.5](m5_22_5_task_details.md)).

## DEVIATIONS LOG

(none)

## FINDINGS

(pending: the task has not been run)
