# Auditor brief

You audit a computation done in two independent blind rooms. You work in two stages. At stage 1 you receive this brief and the spec sheet, and nothing else. At stage 2, after your stage-1 commit, you also receive the rooms' returns with their scripts and data. Work offline. You are given no expected values, and you must not consult any document, repository or web page beyond the files given to you. Report every number you compute. Set `R = 1` in every number you report.

## Stage 1: your own computation, before you see any return

**V2. Recomputation.** Solve the spec sheet's T4 eigenproblem yourself, at all seven widths, by a different discretization: for example P2 elements, finite differences, or a spectral method in `(y, w)`. Run your own convergence study. Report your lowest eigenvalues, your convergence evidence and your error estimates.

**Commit** your scripts, one JSON file holding every number, and a table of your results. Only after this commit do you receive the returns.

## Stage 2: the returns

You now receive two returns, from Room A and Room B. Each is a report, the scripts that produced it, and one JSON file of its numbers.

**V1. Validity.** Rule each return VALID or INVALID, and give the reasons.
- **Explicit requirements.** A return is INVALID if it breaks any explicit requirement of the spec sheet, whatever the numerical effect. That means:
  - the geometry, the normal, the seam condition, the admissible class or the data of sections 2 and 3;
  - T4's protocol as stated there: P1 elements in `(y, w)`, the rectangle, the anti-periodic seam, the area weight, zero data on the arcs and on the collapsed fiber, the node grading `(πR/2)(k/N)²`, `N = 16` per side and 8 cells across, exactly four uniform refinements, `p` from the last three levels, the Richardson extrapolation from the last two, the error estimate `|λ_extrap − λ_finest|`, and all seven widths.

  Some violations, such as a wrong seam sign, can leave the numbers unchanged here. They still invalidate.
- **Sources.** A return is INVALID if there is any sign that the room used a source beyond the spec sheet.
- **Instrument defects** that break no explicit requirement: coding errors, an eigensolver not converged at its own level, or a mismatch between the script, the JSON and the report. Such a defect invalidates if it moves any reported bottom by a relative `10⁻³` or more. Otherwise it may be ruled immaterial, with the reason.

Locate each defect.

**V3. Comparison.** For each return and each width, compare the room's extrapolated bottom with your stage-1 value. Where they differ by more than a relative `10⁻³`, locate the defect:
- in the return, which is then INVALID;
- or in your own computation, which you then correct, recording that the correction was made after the returns were seen.

If you can locate neither, say so.

**G. Grading.** Grade each room's arguments for T1, T2, T3 and T5 step by step, one verdict per step:
- **ESTABLISHED:** the stated route holds.
- **ESTABLISHED, SUPPLIED:** the step holds once a missing part is supplied. Name the part and show it.
- **GAP:** the step cannot be established.
- **DEFECT:** a false statement, or a necessary hypothesis left unstated, even if the result can be rescued.

Then state what each argument establishes as a whole:
- **T1:** the quadratic form, and any boundary or vertex terms.
- **T2:** the classification for each sector, and the realization the admissible class selects.
- **T3:** the spectrum with its labels, and whether it is complete.
- **T5:** the sign, the index and the nullity at each width.

A passing argument is recorded as an audited argument, never as a verified or proven theorem.

**Replacement rooms.** An INVALID return is replaced by a fresh room, with a clean context and the spec sheet alone. It never sees the first return, your audit or your numbers. You audit the replacement's return by stage 2's rules, against your stage-1 values as committed.

## What to return

- **Stage 1:** your scripts, one JSON file of every number, and your table, committed before anything else is shared with you.
- **Stage 2:** a report with V1 for each return, V3 for each return and width, and G for each room and task, with any stage-1 corrections recorded as such.
