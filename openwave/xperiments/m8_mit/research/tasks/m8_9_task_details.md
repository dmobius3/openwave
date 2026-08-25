# M8.9: Discrete-spectrum source localization

> **Diagnostic only.** M8.9 locates where the non-real spectrum enters the M8.4 substrate. It does
> not repair anything, runs no nonlinear dynamics, and does not reopen M8.4 P1A, whose outcome is
> final. If M8.9 identifies and motivates an operator correction, that correction and its
> requalification take a further task identity, so a diagnosis cannot quietly become a repair.

## The question

M8.4 P1A established that the imaginary parts in the computed bundle spectra belong to the
assembled matrix rather than to the eigensolver: holding `L_h` byte-for-byte fixed, they are
unchanged from float64 through 50-digit arithmetic
([P1A closeout](../findings/m8_4_p1a_closeout.md)). It did not establish which upstream stage
produces them. M8.9 asks exactly that, and nothing else:

> does the non-real spectrum enter through the base RBF-FD discretization, through the equivariant
> quotient assembly, or through their interaction?

The shipped contamination table cannot answer it, because each nontrivial bundle was inspected only
at its own first-section level, so `rho`, `d_rho` and `lambda = d_rho(d_rho + 2)` move together
across every row. "High harmonic level is hard for this discretization" and "nontrivial fibre
transport manufactures the complex part" fit the same numbers.

## Stages

**S1, trivial-fibre high-level control.** The trivial bundle `R_0` supplies the missing cell, a
trivial fibre at high harmonic level. Its production block is exactly 60x60 at 60 seeds and the
2I-invariant harmonics fill it exactly, multiplicities 1, 13, 21, 25 at `lambda = 0, 168, 440, 624`.
P1A examined only the `lambda = 0` cluster. Governed by
[the filed S1 decision rule](../findings/m8_9_s1_decision_rule.md), frozen and hashed before any S1
spectrum was computed. S1 is a control, not a verdict: its job is to decide whether S2 is needed.

**S2, independent restriction, if S1 does not settle it.** Compare the production quotient
transport assembly against an independent subspace restriction of the cover-space operator, on the
same cloud with the same RBF weights. Scope and the `R_1` pre-commitment are in the same filed rule.

## Why a new identity rather than an M8.4 row

The [#468](https://github.com/openwave-labs/openwave/pull/468) merge commitment records that any
continuation is separately commissioned substrate work rather than a P1A repair, and that future M8
reviews enforce that boundary. A row numbered inside M8.4 would read as reopening the phase by
placement even though its text says otherwise. M8.4 keeps its `pilot BLOCKED` status unchanged.

## Boundaries

No M8.4 thresholds, no amplitude ladder, no label calibration, no nonlinear run on any bundle, no
target sector spent. M8.9 reads a spectrum that already exists in shipped code and classifies it
under a rule frozen in advance.
