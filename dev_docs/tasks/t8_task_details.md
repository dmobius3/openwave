# T8: Research-script performance note: preallocation and pool shape for NumPy lattice work

> Platform task **T8**. Status: 🚧 **STAGED 2026-10-27**, no go yet · Roadmap:
> [`../platform_roadmap.md`](../platform_roadmap.md) · Deliverable lands in
> [`PERFORMANCE_GUIDELINES.md § Python Performance`](../PERFORMANCE_GUIDELINES.md).

## PLANNING

**Scope.** A short, generic section for research scripts that run NumPy on lattice fields (the shape of the M5 ingestion and reproduction runs, where many short jobs each build a 96³ field, or a few long single-threaded jobs each solve one problem). Two practices and one sizing rule, each with a one-line reason and a check: (1) preallocate the step's arrays once and write into them with `out=` instead of building temporaries per step; (2) write periodic shifts as slice adds into the preallocated buffer instead of `np.roll` (which allocates a full copy per call); (3) size a pool by what bounds the job: a pool of many short jobs on large fields is bound by memory bandwidth and saturates at a small worker count, while a set of long single-threaded jobs scales with the core count, so pin `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and `VECLIB_MAXIMUM_THREADS` to 1 in each worker before NumPy is imported and measure the saturation point with a fixed-work batch on the machine at hand. Also: only fixed-work runs compare across machines or settings; a run that iterates to convergence does different work every time. No machine's numbers go into the section: the rules are portable, the numbers are not.

**Why now.** The M5.33 to M5.35 ingestion runs reran author packages as shipped, which is the right thing for a reproduction, and left the pool shape to the operator. A maintainer's fixed-work test of the same lattice math in a preallocated form (verified against the original to round-off) ran about twice as fast alone and more than twice as fast as a pool, while the pool's saturation point did not move; so the gain is in the code, and it is the kind of practice a contributor can apply to any research script without touching the physics. The note belongs with the existing NumPy guidance so it is found there.

**Definition of done.** The section written under `§ Python Performance` with the two practices, the sizing rule and the fixed-work rule; one small script of our own under `dev_docs/utils/` or the research scripts folder that demonstrates the preallocated step against the naive one on a random field and asserts agreement to 1e-12 (the demonstration, not a benchmark table); the section links the script; no machine-specific number anywhere in the diff; the doc checkers clean.

**Gating.** None beyond the user's go; the deliverable is a doc section and a 40-line script.

**Blindspot pass.** (i) A preallocated step changes the arithmetic order in some reductions, so the agreement check must be stated with its tolerance, not as exact equality. (ii) `np.roll` is sometimes the clearer code; the note says when the copy cost matters (per-step on large fields) and when it does not. (iii) Thread pinning helps pools and can hurt a single large BLAS job; the note states the rule per job shape, not globally. (iv) The author's shipped scripts are never edited for speed by us; the note is for scripts this repo writes.

**Consumers.** Any research script that builds lattice fields per step (the M5 kit checks, the reproduction comparators, the planned Hessian ladders); the instrument check each needs is the agreement assertion of the demonstration script.
