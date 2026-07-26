# M5.20.6: the breathing-orbit BVP: the profile-dynamic clock on the loop

**Status**: ❌ ARCHIVED 2026-07-16 (the loop-side reserve). Duda's 2026-07-15 Q24 reply ([`m5_20_convo.md § 2026-07-15`](m5_20_convo.md)) DEFERRED the BVP question ("I think about it - can elaborate, but generally it seems quite difficult"; "numerically it doesn't seem practical (?)") and REDIRECTED to the electron hedgehog: the folding table's electron-redirect row fires exactly as pre-registered. **This stub archives unrun as the loop-side reserve**; the program continues at [M5.21.1](m5_21_1_task_details.md) (user decision 2026-07-16: electron-first). Reopen only if his promised Q24 elaboration lands AND the user routes back to the loop.

> Original stub (2026-07-14, written for the CONFIRM branch of [Q24](../m5_question_tracker.md#q24-detail)), preserved below unmodified:

**Lineage**: [M5.20.5](m5_20_5_task_details.md) (the rigid level measured OUT: the directional block = the loop's static force is 99.9997% time-row, a sector rigid conjugation cannot trade energy with; a breathing profile can) · [`../findings/m5_20_5_method_note.md`](../findings/m5_20_5_method_note.md) · [M5.20.4](m5_20_4_task_details.md) (the balance roots + the φ-average correction) · the M5.12 phase-D time-periodic action instrument (the container precedent) · canonical [§ 2 rows 4/4b](../m5_theory_canonical.md).

## TASK PLANNING (2026-07-14 stub; finalize at PLAN when the answer lands)

### Scope

Extend the free-period least-action BVP from rigid conjugation orbits to PROFILE-DYNAMIC (breathing) orbits on the loop, and converge an actual extremal: the particle clock with a time-periodic profile. One harmonic first:

```text
M(x, t) = Lam(w t) [ M0(x) + A(x) cos(w t) + B(x) sin(w t) ] Lam(w t)^T
S(w; M0, A, B) = INT_0^{2pi/w} dt [ T_true - U ]     (free period, phi-averaged kinetic per time slice)
extremize over (M0, A, B, w);  the dS/dw = 0 identity (generalized virial) derived + gated at F0.
```

| Step | Content | Kill / survive (pre-registered) |
| --- | --- | --- |
| F0 | Instrument gates: time-quadrature action + all gradients complex-step gated; the time-integrand of the quartic L with a one-harmonic ansatz is band-limited (harmonics ≤ 4 in ωt: derive the exact bound, gate the quadrature count like the nphi lesson); **the rigid-limit regression gate**: A = B = 0 must reproduce the M5.20.5 rigid numbers exactly; the φ-average applied per time slice (the § 3 equivariance boundary); gauge fixing (time-translation phase) verified to remove the solver null | gates green before any physics |
| F1 | **THE CHANNEL GATE (cheap, decisive, minutes)**: a2x-style alignment in the ENLARGED space: does the breathing sector's gradient block couple to the loop's time-row static force? Measure the cos alignment + sector split at the seed; natural A-seed candidates: the time-row force direction itself + the leading negative-inertia eigenvector of K10 (the theory's own boost-sector mode as the clock DOF: "negative Hamiltonian terms propel") | cos bounded away from 0 = the channel is OPEN (proceed); cos ≈ 0 again = the one-harmonic class is blocked too: STOP EARLY, ship the characterization, the formulation question goes back sharpened |
| F2 | The one-harmonic extremal solve at 2-3 working points (rigid working points + the F1-informed seeds; warm/cold starts): the M5.20.5 solver suite (descent on \|gradS\|², Newton-Krylov, L-BFGS) on the enlarged unknowns (M0, A, B, ω) | bar (pre-registered, unchanged): rel residual ≤ 1e-3, q intact, the F0 virial identity as the certificate; converged = THE FIRST TRUE-L CLOCK; all-fail = the honest record + what blocks |
| F3 | The harmonic ladder: k = 2, 3 if k = 1 signals; convergence of the extremal in harmonic content | monotone residual/energy behavior across k |
| F4 | THE OBSERVABLES on the converged orbit: clock ω* vs the chirped ladder at the ring; the RADIUS-BREATHING amplitude (Duda's own oscillation picture, read directly); energy split + the virial certificate; containment r50/r90 vs baseline; winding + core spectrum on-orbit | measured on the converged orbit; partials labeled 🔶 at the best state |

NOT in scope: the electron hedgehog (M5.21.1 inherits the winner), any kinetic-term modification (unless his answer sanctions one: then the gated s2/qc EL library from M5.20.5 is the implementation layer, folded here at PLAN).

### The folding table (his answer routes the plan)

| Duda's Q24 answer | Route |
| --- | --- |
| Confirms BVP + breathing | run this plan as written (finalize working points at PLAN) |
| Adds a constraint surface / kinetic normalization | fold as F0 gates + the modified functional; the s2/qc EL machinery implements any sanctioned term |
| Names a rigid class we missed | insert an F1-style alignment probe of that class BEFORE the breathing solve (cheap); proceed by its verdict |
| Redirects to the electron | STOP: re-plan as M5.21.1-first (this stub archives as the loop-side reserve) |

### Definition of done (stub level)

| ✅ when | Bar |
| --- | --- |
| F1 verdict | the channel gate measured (open/blocked), whichever way |
| F2 verdict | a converged breathing extremal (bar above) with F4 observables, OR the pre-registered non-convergence record + characterization |
| Records | method note `../findings/m5_20_6_method_note.md` (equations first, rigid-limit regression gate shown); independent adversarial audit (cardinal rule); tracker Q24 + roadmap + checkpoint routing; field-state prints + film standard where states evolve |

### Blindspot pass (stub; re-run at PLAN)

| # | Unknown unknown surfaced | Fold |
| --- | --- | --- |
| 1 | Time-quadrature aliasing (the quartic L multiplies harmonics) | derive the exact band limit, GATE the quadrature count (the nphi lesson: measured, not assumed) |
| 2 | The breathing sector may re-open the ill-posed IVP modes as action SADDLE directions | extremize by residual (gradS → 0), never by naive S-minimization; the M5.20.5 solver suite carries over |
| 3 | Gauge freedom (time-translation + Λ phase) leaves solver null directions | fix a phase condition at F0 and verify the Hessian null count drops accordingly |
| 4 | The rigid-limit regression could silently break under the new container | A = B = 0 must reproduce M5.20.5's numbers to machine precision (F0, mandatory) |
| 5 | Warm-start bias toward the rigid basin | seed A from the time-row force direction AND the negative-inertia eigenvector (two starts, pre-registered) |
| 6 | Cost: each action eval = n_t × nphi kinetic evals | budget from the measured 0.09 s/grad_q2_avg; pick n_t from the gated band limit, not generosity |

### Research-body destinations

| Artifact | Destination |
| --- | --- |
| Scripts | `../scripts/m5_20_6_bvp.py` · `m5_20_6_plots.py` |
| Data / plots | `../data/m5_20_6_*.json` (+ npz ≤ 1 MB) · `../plots/m5_20_6_*.png` |
| Findings | `../findings/m5_20_6_method_note.md` |
| Records | this file (FINDINGS + TASK REVIEW) · tracker Q24 · `../checkpoints/m5_20_6_progress.md` (opens at go) |

### Gating

Roadmap `Gated By`: **Duda's Q24 answer (this stub assumes the CONFIRM branch: re-plan first if it differs)** + user "go". Not runnable before the answer: the whole point of the park is to spend the next run on HIS formulation, not another self-directed guess.

### Model + effort

Fable 5 high for the F0 container derivation (new functional + virial identity + gauge fixing) and F2 solver work; default for the mechanical gates; independent agent for the audit (cardinal rule).
