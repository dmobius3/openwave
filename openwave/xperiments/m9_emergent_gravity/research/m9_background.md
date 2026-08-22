# M9 background: what NSM has, what it lacks, and why the column exists

> **Purpose.** This is the gap map behind the M9 column: an honest statement of what
> the Emergent Gravity / New Standard Model program brings, what it is missing, and
> what a gravity-certification column on this platform is for. It was written by the
> maintainers at admission time (2026-08-15) from the author's application
> ([discussion #442](https://github.com/openwave-labs/openwave/discussions/442)), the
> first package ([PR #441](https://github.com/openwave-labs/openwave/pull/441)), and
> the maintainer verification run recorded in that PR's review. The author owns the
> science; corrections via PR or discussion are welcome.

> **Archive pointer (2026-08-22).** Rows below that cite M9.3 to M9.73 results describe the author's parked campaign records, not files in the tree; they live in PR #441's branch history at `a5640709` and return as per-task PRs.

## 1. What NSM is

NSM is an **effective-theory assembly, not an emergence program**: Palatini
Einstein-Cartan gravity (coframe \(e^a\), independent Lorentz connection
\(\omega^{ab}\)) plus the empirically installed Standard Model. Matter is input;
particles are SM field excitations, not defects of a lattice medium. The one
interaction beyond SM + Einstein gravity is the Hehl-Datta axial-axial contact
\(\mathcal{L}_{\mathrm{HD}} = -\frac{3\kappa}{16} J_5\cdot J_5\), with the
coefficient fixed by \(G\) and no adjustable knob. The program's distinctive claim,
that entanglement bookkeeping *selects* Einstein-Cartan, is a holographic argument
certified (by the author's papers) only in AdS at linear order, and is explicitly
not the first in-platform task.

That makes M9 the mirror image of the emergence columns: they build particles and
hope gravity follows (M5's GEM row is ⚠️ with the inter-mass force not yet
computed); M9 installs the particles and works only the gravity sector. The
specification of record is [`m9_theory_canonical.md`](m9_theory_canonical.md);
the paper series lives in the author's repo
([n4hy/New_Model_Emergent_Gravity](https://github.com/n4hy/New_Model_Emergent_Gravity)).

## 2. The gap map (what the M9 program must supply)

| Gap | The honest statement | What supplies it |
| --- | --- | --- |
| Particle rows out of reach by design | The SM is installed, so every particle, charge, and wave-emergence criterion stays 🚧 permanently for this column; scoring them would be a category error ([`m9_theory_canonical.md`](m9_theory_canonical.md) § 7) | Nothing. The column claims gravity rows only |
| The selection claim is not in-platform | "Entanglement selects Einstein-Cartan" rests on holographic certificates (FGHMV-class) that exist as papers, not as platform scripts | Later tasks, one certificate at a time, each needing its own runnable artifact before anything is cited as in-platform |
| Newton limit unproven here | Attractive \(1/r^2\) from the inherited Einstein equations is pre-registered (task M9.2) and not run | M9.2, on its locked gates |
| de Sitter / cosmology obstructed | Parked campaign record: the FGHMV-standard copy onto the cosmological horizon fails on sign and isometries (M9.6); Jacobson 1995/2016 is not a substitute (M9.7); Einstein+\(\Lambda\) from a cosmological CFT is not claimed | Open. Recorded as documented negatives, which is the shape this platform wants |
| \(I_B\) multi-digit coefficient | Parked campaign record: negative on three extraction routes; no source-independent multi-digit coefficient exists | Open, with the negatives on record |
| Second-order Einstein-Cartan | Parked campaign record: metric Einstein at second order is cited (FHHPRV 2017), not re-derived; the axial matching is obstructed | Open |

## 3. Evidence weight (the grading the column starts with)

| Layer | Weight | Why |
| --- | --- | --- |
| The M9.1 algebra gate | verified in-platform | maintainer-reproduced by two independent routes (the solver and a sympy second method); r = 3/16 in both signatures; the mutation checks move the ratio, so the gate can fail |
| The C2 spin-dual defect | a documented negative against the author's own paper | the paper's printed \(s = -\frac14(*J_5)\) measures \(-\frac12\); the author shipped the FAIL as a scored gate rather than adjusting the paper silently. Adopted as the column's honesty baseline |
| The campaign negatives (\(I_B\), axial EC, first-law kernel probes) | documented negatives | each carries a script and a finding note; none moves a cell, and none claims more than its own lattice or expansion reaches |
| The holographic selection argument | paper-tier, out of platform scope | cited, never consumed as an in-platform result until a script exists per certificate |
| The paper series | preprint-tier, same-day | all venues unrefereed at admission (Zenodo; the sampled DOIs resolve); the standing red flag is recorded in [`../theory/_CITATIONS.md`](../theory/_CITATIONS.md), same as the M8 precedent |

The practical rule for M9 tasks: **target the gravity certificates one runnable
artifact at a time, and let negatives stand as results.** The column's credibility
is its verification density, not its claim count.

## 4. Onboarding evaluation of record (2026-08-15)

| Check | Result |
| --- | --- |
| Provenance | author account genuine and active since 2010; the declared specification repo exists and is public (CC-BY-4.0); sampled Zenodo DOIs machine-verified resolving; DCO signed on every commit of the first package |
| [`ONBOARDING_MODELS.md`](../../../../ONBOARDING_MODELS.md) STEP 1 fit | gravity rows native, everything else honestly declared out of scope; falsifiers named with current bounds; free-parameter ledger disclosed (zero free parameters is explicitly NOT claimed) |
| Claim-to-artifact | the headline gate recomputed by the maintainers from the shipped scripts and data, agreement to 1e-15; second-method audit CONFIRMED / REFUTED per criterion, matching the author's own statements |
| Rigor culture | pre-registration before runs, mutation-tested gates, negatives led with unprompted: compatible with [`AI_HYGIENE.md`](../../../../AI_HYGIENE.md) from day one |
| Sequencing | the first package predated the application by an hour and carried a full research campaign; admission terms therefore fixed the working mode below |

Admission terms of record (from discussion #442): the rigor bar as demonstrated;
**one PR per concern**, draft while work continues, review on a quiescent head;
**the author creates and manages the tasks on [`m9_roadmap.md`](m9_roadmap.md)**,
with every script, data file and finding carrying the id of a registered roadmap
row; papers stay in the author's repo with
[`../theory/_CITATIONS.md`](../theory/_CITATIONS.md) as the tracked pointer;
review and merging follow maintainer schedule and resource availability.
