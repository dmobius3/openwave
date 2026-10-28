# M8.17 method note: the slot dynamics census run

Task: [`tasks/m8_17_task_details.md`](../tasks/m8_17_task_details.md), governing SHA-256 `667159c3db993396…`, pre-registered in [#662](https://github.com/openwave-labs/openwave/pull/662), run from 2026-10-27T19:50:47Z to 2026-10-28T16:18:02Z on the census run's code frozen at MIT `338b873`, under the maintainer's [go](https://github.com/openwave-labs/openwave/pull/662#issuecomment-6443460379) and the [freeze block](../m8_17/freeze_block.md) of [#673](https://github.com/openwave-labs/openwave/pull/673). The record is [`RECORD.md`](https://github.com/dmobius3/mode-identity-theory/blob/7b66f8a0e207888042e6bb4a9ab636f946cb2a81/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/RECORD.md) beside the code. The run completed, with no stop, and the reproduction gate held at every known point it compared. Of the 56 first-time readings, 42 are classified. Four rows meet an S4 completion condition: T2/R4 and T3/R5, with each change of type located, and T2/R5 and T3/R4, Elliptic throughout. Five end not settled by this run: T1/R4, T4/R4, T1/R5, T4/R5 and T5/R2. This note states the classification, the locator and the row outcomes as equations, maps each to its code, gives the results against the frozen terms, and says what the run does not verify. Adjudication is the maintainers'.

## 1. Equations first

The instrument is M8.16.1's, byte for byte, as [M8.16.1's method note](../findings/m8_16_1_method_note.md) and, through it, [M8.16's](../findings/m8_16_method_note.md) give it: the law, the discretization, the exact kernel, the bordered Newton, the seeds and continuation, the linear record, K1, the convergence rule with its band-edge guard, and branch identity between `N = 24` and 32. M8.17 scores no instrument. It classifies points, brackets and locates changes of type, and records each row's outcome.

### 1.1 The points

```text
grid:      ε_k = 0.25 · 2^(k/2)
windows:   R4: k = 0, …, 16, and ε = 72 (18 points);   R5: k = 0, …, 14 (15);   R2: k = 0, …, 18, and ε = 132 (20)
paths:     at each N ∈ {16, 24, 32}, from the case's seed at ε = 0.25, along the grid, each point continued from the last;
           a solve that fails or is not accepted ends that cutoff's path
```

The run computes T2/R4 to 22.63, T3/R4 to 32, T5/R2 to 132, and T1 and T4 in both slots over their whole windows. T2/R4's path goes on to 32 only where its type bracket is (22.63, 32), and here it was (16, 22.63). T2/R5 and T3/R5 rest on known outcomes.

### 1.2 Classifying a point

```text
here(ε):    for every N:   accepted(ε, N)   and   K1(ε, N)
            reading(ε) = converged             the convergence rule, on the three cutoffs' summaries
            identity(ε; 24, 32)                aligned distance ≤ 10⁻² ‖φ‖, and the two records equal
            type(ε) = the rule's type,  with  (n(L), n(D), K)  and  k_i⁻  at  N = 32
there(ε):   M8.15 decided ε;   type(ε) = its adjudicated reading
classified(ε)  ⇔  here(ε) or there(ε),  and where both hold, they agree in type, K and k_i⁻
```

### 1.3 Brackets and the locator

```text
type bracket:    consecutive grid points a < b, both classified, type(a) ≠ type(b)
index bracket:   a type bracket with n_N(a) and n_N(b) resolved and n_N(a) ≠ n_N(b), at N = 24 and at N = 32
bisection, at each N:
    m = (a + b)/2, solved by continuation from a and from b
    the two solutions agree:  aligned distance ≤ 10⁻² ‖φ‖,  n_N equal;   K1 holds;   n_N read by inertia
    keep the half whose ends differ in n_N
    complete when  b − a < (w/2) · a,   w = 10⁻³
located:
    the bisections at N = 24 and 32 complete, with no unsuccessful stop
    their intervals overlap, and their union U = [u⁻, u⁺] has  u⁺ − u⁻ < w · u⁻
    Δn_24 = Δn_32 ≠ 0
    guards at u⁻(1 − δ) and u⁺(1 + δ), δ = 0.01, then 0.05 once, each solved from its own end and classified
    on this instrument with that end's type
```

An unsuccessful stop is a missing or unreproduced end, a failed or unaccepted solve, K1 failing, the two solutions disagreeing, `n(L)` unresolved, or `n(L)` matching neither end. `N = 16` is bisected where its ends' `n(L)` is resolved and different, and reported otherwise; it does not decide.

### 1.4 The row outcomes

```text
Elliptic throughout:           every window point classified Elliptic
Each change of type located:   every point classified, or between a located change's guard readings;
                               every type bracket located
Not settled by this run:       otherwise
```

Whether a row meets an S4 completion condition is stated separately. A row not settled by this run is not closed under S6 here, and a point it does not classify is not offered to this instrument again without a ruling on #512.

### 1.5 The reproduction gate

```text
rule 1:   a point M8.16.1 read, reached by its own steps:   the whole record equal, timings aside
rule 2:   a point M8.16.1 withheld:                         accepted, ε reached, λ, m = ‖φ‖², floor, evaluations equal
rule 3:   T1's and T4's ε = 0.5, 1 and 2:                   type, indices, V₀'s dimension and ranks, chains,
                                                            Krein counts and the rule's verdict equal
```

The first point that does not reproduce stops the run, and nothing is repaired or rerun.

## 2. Equation-to-code map

The links are pinned to MIT `338b873`, the frozen code. The instrument's map is [M8.16's note](../findings/m8_16_method_note.md), § 2, and M8.16.1's driver and checks are mapped in [its note](../findings/m8_16_1_method_note.md), § 2. The new code is in `run/`: `census_lib.py`, the classification, the locator and the outcomes, with no solve; `census.py`, the run; and `write_guard.py` and `freeze_block.py`, the freeze.

| Equation | Code |
| --- | --- |
| § 1.1 the paths | [`path_points`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L81), [`Cutoff`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census.py#L78) |
| § 1.2 here(ε) | [`here`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L218), [`readings_of`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census.py#L133) |
| § 1.2 there(ε), and the precedence | [`there`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L248), [`classify`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L261) |
| § 1.3 the brackets | [`changes`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L302), [`nl_at`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L279), [`bracket_kind`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L289), [`end_kind`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L603) |
| § 1.3 the bisection | [`bisect`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L325), [`aligned`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census.py#L125) |
| § 1.3 the location and the guards | [`location`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L363), [`guard_offsets`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L393), [`guards`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L401), [`span`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L425) |
| § 1.3 one bracket, end to end | [`locate`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L441), [`locate`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census.py#L223) |
| § 1.4 the row outcomes | [`row_outcome`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L509), [`meets_s4`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L527), [`classify_rows`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L574), [`assemble`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L624) |
| § 1.5 the three rules | [`gate_full`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L171), [`gate_solve`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L182), [`gate_categorical`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L198), [`rule_for`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census_lib.py#L132), [`gate_row`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census.py#L199) |
| the run: its order, the stop and the resume | [`main`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census.py#L236), [`stop_run`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census.py#L188), [`resume`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census.py#L156) |
| the freeze: written once, then admitted | [`write`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/write_guard.py#L72), [`frozen`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/census.py#L58) |
| the freeze block | [`build`](https://github.com/dmobius3/mode-identity-theory/blob/338b87365fea6e2b31bb9bb861851b2d46b7d4b5/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/run/freeze_block.py#L39) |

## 3. The run

- **The environment.** Python 3.13.13, NumPy 2.5.0 and SciPy 1.18.0, on Accelerate's BLAS and LAPACK as each reports them, under macOS 15.8 (24H23) on one Apple M3, as `run/ENVIRONMENT.json` records it. It was read just before the launch and again after the run, with no difference.
- **The order of the freeze.** The code on MIT `main` at `338b873`, with its manifest pin, `2deed9aa564a6ccb…`; the terms merged in #662; the maintainer's go; `run/write_guard.py`'s `FROZEN.json` at 2026-10-27T15:20:25Z, with its pin, `bea8d27a537198b3…`; the freeze block, whose pull request was opened at 2026-10-27T16:41:23Z and merged at 2026-10-27T18:59:56Z; then the run. `census.py` refuses to start unless the freeze validates.
- **The launch.** One launch, on mains power with sleep held off. The command checked a read-only resume check against its reviewed SHA-256, `2c4762ef369f8421…`, ran it, and only then `run/census.py`, so that no relaunch could pass a stop the run's log had not saved. The check sat outside the pinned tree and changes nothing in it. The console log opens with the hash check's `OK` and the check's `PASS`, then holds the run's transcript line for line.
- **The run.** It started at 2026-10-27T19:50:47Z and wrote its record at 2026-10-28T16:18:02Z, 73635 s later. The seven rows ran in the terms' order, then the locator's two brackets. Its log records the manifest pin at the start, one start, no resume and no stop.
- **Deviations.** None in the run.

## 4. Results against the frozen terms

[`RECORD.md`](https://github.com/dmobius3/mode-identity-theory/blob/7b66f8a0e207888042e6bb4a9ab636f946cb2a81/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17/RECORD.md) gives every point with its classification, its instrument and this instrument's reading (§ 3), every bracket (§ 4) and each row's ending (§ 5).

| Row | Slot | Outcome | Meets an S4 completion condition | Classified | Open |
| --- | --- | --- | --- | --- | --- |
| T1/R4 | `R4` | Not settled by this run | no | 13 of 18 | 1.41, 32, 45.25, 64 and 72 |
| T4/R4 | `R4` | Not settled by this run | no | 13 of 18 | 4, 8, 11.31, 64, 72 and the bracket (2.83, 5.66) |
| T1/R5 | `R5` | Not settled by this run | no | 14 of 15 | 32 |
| T4/R5 | `R5` | Not settled by this run | no | 15 of 15 | the bracket (8, 11.31) |
| T2/R4 | `R4` | Each change of type located | yes | 18 of 18 | none |
| T2/R5 | `R5` | Elliptic throughout | yes | 15 of 15 | none |
| T3/R4 | `R4` | Elliptic throughout | yes | 18 of 18 | none |
| T3/R5 | `R5` | Each change of type located | yes | 14 of 15 | none |
| T5/R2 | `R2` | Not settled by this run | no | 17 of 20 | 45.25, 128 and 132 |

- **The reproduction gate** held at all 63 known points it compared: rule 1 at 47 points, rule 2 at T2/R4's 16 and 22.63, T3/R4's 32 and T5/R2's 45.25, and rule 3 at T1's and T4's 0.5, 1 and 2 in both slots, 12 points. The locator's two recomputed paths, 14 points on T2/R4 and 12 points on T4/R5, reproduced before their bisections began.
- **The first-time readings.** Of 56, 42 are classified on this instrument and 14 are not: T1/R4 at 1.41, 32, 45.25, 64 and 72; T4/R4 at 4, 8, 11.31, 64 and 72; T1/R5 at 32; and T5/R2 at 45.25, 128 and 132. The four points M8.16.1 withheld, read here for the first time: T2/R4 at 16, Elliptic, (3, 3, 0); T2/R4 at 22.63, Hyperbolic, (5, 3, 2); T3/R4 at 32, Elliptic, (9, 3, 6); and T5/R2 at 45.25, Unresolved.
- **T2/R4's change of type** is located: "a change of type located between the guard readings, with `n(L)` changing by +2 inside the interval at `N = 24` and 32". Its `n(L)` interval is [16.2718, 16.2783] at `N = 24` and 32 alike, 4.0×10⁻⁴ of its lower end. The lower guard, at 16.1091, reads Elliptic, (3, 3, 0), and the upper, at 16.4411, Hyperbolic, (5, 3, 2), each at `δ = 0.01`. `N = 16`'s interval, [16.1942, 16.2006], is reported and does not decide.
- **T4/R5's bracket,** (8, 11.31), Elliptic to Hyperbolic, an index bracket with `n(L)` 7 to 6 at `N = 24` and 32, is not located. Its bisections at `N = 24` and 32 each stopped at their tenth midpoint, `ε = 8.7216`, where K1 fails on both solutions, from either end. Their interval was then [8.7184, 8.7249], 7.4×10⁻⁴ of its lower end, not yet below `w/2`. `N = 16`'s bisection completed, at [8.8511, 8.8543], and does not decide.
- **T4/R4's change of type,** Elliptic at 2.83 to Hyperbolic at 5.66, lies across the unclassified grid point 4. It has no type bracket and is not located.
- **T3/R5's change** keeps M8.15.2's location, which accounts for `ε = 16`.
- **No disagreement** between the instruments at any point both classify.

## 5. Findings after the run

The readings in this section are post-run findings, not further tests. No outcome is assigned here.

### 5.1 Where this instrument leaves a point Unresolved

On the seven rows the run computes, this instrument leaves 25 points unclassified, counting the two beyond T3/R4's path that M8.16.1 read, and 11 of them are classified by M8.15's readings. Every one carries a set whose count differs between cutoffs or a group or set not contracting above the floor, the two kinds the terms' § STATED LIMITS anticipated from the grouping check of M8.16.1's Unresolved points. By kind, a point counted under each it carries: a set whose count differs between cutoffs, at 21; a group or set not contracting above the floor, at 8; a record's structure differing along the ladder, at 2; a group's `U` above a tenth of its distance to the imaginary axis, at 2; branch identity failing, at 2; K1 failing at `N = 16`, at 1. This run does not test whether these count changes are the closure window's, as the grouping check found at nine of M8.16.1's ten (the terms, § STATED LIMITS: WHAT IS ADDED), and it claims no type at any of them. `RECORD.md` § 3 gives every point's reasons as the run states them.

### 5.2 T4/R5's bracket

On T4/R5, `n(L)` falls from 7 at `ε = 8` to 6 at 11.31 at all three cutoffs. At `N = 24` and 32 the bisections kept the same half at every step, and every midpoint before the last held K1 on both solutions. At the last, `ε = 8.7216`, both solutions were accepted, with `n(L) = 7`, and K1 failed on both, at both cutoffs: the run's unsuccessful stop. As the run's records state it, at `N = 24` "the zero cluster is not separated: an eigenvalue lies in (margin/3, 3 margin)", with `outside_min` 3.69×10⁻⁶ and 3.70×10⁻⁶ against `margin` 2.03×10⁻⁶, and "V0 does not match its expectation", with a symmetry leak of 9.4×10⁻⁶ and 8.2×10⁻⁶; at `N = 32` the separation alone, `outside_min` 6.12×10⁻⁶ on both solutions against `margin` 2.72×10⁻⁶, with `V₀` matching. At the midpoint before, K1 held, with `outside_min` 2.48×10⁻⁵ and 2.53×10⁻⁵. `N = 16`'s bisection parted from theirs at their common midpoint `ε = 8.8284`, where it read `n(L) = 7` against their 6, and completed near 8.85. Why K1 fails at the last midpoint is OPEN. The run reports the stop and reruns nothing.

### 5.3 T5's upper window

The terms expected T5's 45.25, 128 and 132 to end Unresolved on this instrument (§ THE ROW OUTCOMES), and they do. Their causes are, at 45.25, a set whose count differs between cutoffs; at 128, a set whose count differs between cutoffs, a group's `U` above a tenth of its distance to the imaginary axis and branch identity failing; and at 132, a group or set not contracting above the floor, a record's structure differing along the ladder, a group's `U` above a tenth of its distance to the imaginary axis, branch identity failing and K1 failing at `N = 16`. The census page may consider closing them as instrument-limited under S6, where its condition is established. This run makes no S6 ruling.

## 6. Provenance

- **Terms:** `tasks/m8_17_task_details.md` at OpenWave `a4c7a270`, SHA-256 `667159c3db993396…`.
- **Code:** MIT `338b873`. Its manifest pins 35 files, `2deed9aa564a6ccb…`, and `SHA256SUMS`, regenerated at the record's commit, lists every file of the folder.
- **Freeze:** `run/FROZEN.json`, pin `bea8d27a537198b3…`, carried in full by [the freeze block](../m8_17/freeze_block.md).
- **The resume check:** `launch/resume_check.py` in the record, SHA-256 `2c4762ef369f8421…`, with its 22 arms.
- **Record:** MIT [`7b66f8a`](https://github.com/dmobius3/mode-identity-theory/tree/7b66f8a0e207888042e6bb4a9ab636f946cb2a81/files/framework/files/working/files/scripts/soft-slot-dynamics-m8-17). It holds `RECORD.md`, `run/out/`, the freeze's files, the resume check and the post-run files.
- **Reproduction:** `REPRODUCE.md` beside the code, and `RECORD.md` § 10.

## 7. What this run does not verify

- **A type at an Unresolved point.** This instrument does not classify a point its rule leaves Unresolved, and no registered test here tells a genuine degeneracy from the instrument's resolution, so S4's third ending is not reached (§ 5.1).
- **Why K1 fails in T4/R5's bisection.** The record shows where it failed. The cause is OPEN (§ 5.2).
- **The location's identification** with the genuine index crossing. T2/R4's location is an `n(L)` interval with two guard readings. The Hamiltonian–Krein Identification's Corollary I is a conditional test, and it is not applied here.
- **The persistence of any Elliptic stretch.** S5 needs a pre-registration of its own, outside this admission.
- **The calibration's range.** M8.16.1's stated limits carry over (the terms, § STATED LIMITS: WHAT IS ADDED): the margin was fitted in `R2`, `R3` and `R4`, so T1's and T4's classifications in `R5` lie outside the fitted slots. K1 held at every grid point of the seven computed paths but T5/R2's 132, where it failed at `N = 16`.
- **Other machines.** The run is on one machine, and bitwise agreement is not claimed.
- **The platform's criteria.** No `MODELS.md` cell moves, and M8.7's gate is unchanged.

## 8. Audit record

Before filing, the record and this note were checked in two rounds against the frozen terms, the frozen code and the run's files. The checks rebuilt `RECORD.md` from the files and reran its eight refusal arms. They reclassified every computed point from the raw case files with separate code, matching `CENSUS.json` at each point, with no disagreement at the 68 points both instruments classify. They also traced both locator records midpoint by midpoint; recounted the 56 first-time readings, the 63 points through the reproduction gate and the causes at the 25 Unresolved points; reran the pre-freeze tests on copies of the frozen folder; and rebuilt this note against clones of both repositories, each code link on its definition. They led to five changes:
- K1's failing parts at T4/R5's stopped midpoint, stated as the run's records give them;
- the calibration's range among this note's stated limits;
- `SHA256SUMS` regenerated over the published folder;
- the launch command's test on stubs, published beside the resume check;
- one clause on T2/R4's path in § 1.1.

All were adopted, and no reading changed.
