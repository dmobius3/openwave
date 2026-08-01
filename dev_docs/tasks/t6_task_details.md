# T6: SMT-certified discrete topological invariants (contributor offer)

> Roadmap row: [`../platform_roadmap.md`](../platform_roadmap.md). Status: 🔶 **DEFERRED
> 2026-08-01, on arrival**. Owner: unassigned. Filed from GitHub issue
> [#298](https://github.com/openwave-labs/openwave/issues/298) (opened 2026-07-16 by
> `chadbrewbaker`) when [T5](t5_task_details.md) settled that tasks live in roadmaps. The row
> exists to keep a promise made in that thread, not because work is queued: the maintainer reply
> said the fit would be analysed and followed up if a consumer landed, and a closed issue is not a
> place a promise survives.

## PLANNING

### What was offered

An outside contributor benchmarked SMT (z3) encodings against small lattice topology and offered
to implement them properly, with tests, if the platform has a use. The measured results, from the
thread:

| Result | Time |
| --- | --- |
| Per-plaquette winding bound `\|w\| ≤ 1` proved UNSAT on a `Z_k` clock discretization | 18 ms |
| Total vortex charge on a torus is zero, proved as UNSAT | 19 ms |
| Vortex/antivortex pair constructed at prescribed plaquettes | 28 ms |
| Charge-1 vortex built on an open grid, with the boundary carrying the charge | 16 ms |

The finding that travels further than the timings is about encoding, not solving: giving each edge
one difference variable, reused with opposite signs by its two adjacent plaquettes, makes
`Σ w_p = 0` syntactic and collapses a 15-second timeout to 19 ms of propagation. Stated generally,
the solver does not find topology; the encoding carries it, and an impossibility proof is cheap
exactly when the variable sharing mirrors the boundary operator.

A second arm used Burnside / cycle-index counting to predict symmetry reduction before the solver
ran (100 raw states → 5 orbits under `C₄ × Z₅`; allSAT enumeration 98 ms → 24 ms, the predicted
factor of 5), which makes the count and the solver a correctness check on each other.

### Where it could slot, and where it could not

The contributor's own framing is the honest one and is adopted here: SMT is credible for
**certified discrete initial conditions** and **impossibility checks on small lattices** ("no
configuration with this charge profile exists under these boundary conditions"), and is not a
competitor to the platform's PDE solvers. Two plausible consumers, neither of them scheduled:

| Candidate | Why it might fit |
| --- | --- |
| A certified seed for a defect run | Today's seeds are constructed and then measured; a certificate that a seed's charge profile is the only one admissible under its boundary conditions would be an input guarantee rather than an output check |
| A machine-checked discrete invariant in a task that turns on one | Where a winding or charge-sum argument is load-bearing for a verdict, a mechanical proof on the discrete lattice is a second route to it, independent of the measuring script |

### Why deferred

No consumer exists in the current program: the live work is continuous-field PDE runs, and nothing
queued turns on a discrete impossibility proof. Filing a task with no consumer would put a
contributor's effort behind a use that has not been established, which the thread explicitly
avoided ("no need to implement anything speculatively in the meantime").

**Re-open trigger**: a task whose verdict turns on a machine-checked discrete invariant, or that
wants a certified seed, plus the contributor's availability at that time. Whoever re-opens it
starts by asking whether the offer still stands.

### Blindspots

| Risk | Guard |
| --- | --- |
| The row becomes a permanent parking space, and the offer quietly expires | the re-open trigger names a concrete condition rather than a date, and the thread stays linked so the offer's own terms are readable |
| Small-lattice certificates get read as claims about the continuum | the scope line above is the contributor's own and stays attached: certified discrete initial conditions and impossibility checks, never a statement about the PDE limit |

## DEVIATIONS LOG

(none)

## FINDINGS

(pending: deferred on arrival, nothing run)
