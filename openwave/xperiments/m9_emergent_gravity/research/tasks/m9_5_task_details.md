# M9.5: Q4a in the proper context (holographic-pair selection)

> Spec of record: [`../m9_theory_canonical.md`](../m9_theory_canonical.md).
> Q4a is the monograph's question: the microscopic description is the
> boundary theory of a holographic pair; is that pair selected?
> Official M9.5 is this task. The file `m9_5_ec_symplectic.py` is the
> Paper 14 campaign script, not this task.

## TASK PLANNING (2026-08-15)

### Scope

Answer Q4a as a selection problem. Do not invent a CFT.

- **Selection.** Do CHM, FGHMV, and the Paper II modular rule single
  out a pair of which the NSM is the bookkeeping?
- **Existence.** Does at least one such pair exist?

Reading B (the SM *is* the boundary CFT) is tested by necessary
conditions for a holographic CFT. Reading A (SM fields are bulk;
the CFT is something else) is tested by underdetermination.

### Pre-registered claims

| ID | Claim | Pass | Fail |
| --- | --- | --- | --- |
| C1 | SM one-loop \(b_i\) equal \((41/10,-19/6,-7)\) and are all nonzero | exact fractions | any other triple, or a vanishing \(b_i\) |
| C2 | SM light on-shell dof (before EWSB) are \(O(10^2)\), not a large-\(N\) limit | total \(<200\) | a parametric \(N\to\infty\) count |
| C3 | CHM / FGHMV use \(T_{\mu\nu}\) only; \(G_{\mathrm{SM}}\), \(n_g\), Yukawas do not appear | independence table | a first-law formula that contains flavor data |
| C4 | Paper II selects Young symmetry of the spin source, not \(G_{\mathrm{SM}}\) | stated scope | a derivation of \(SU(3)\times SU(2)\times U(1)\) from modular flow |
| C5 (mutation) | \(\mathcal{N}=4\) SYM has one-loop \(b=0\) | exact zero | the CFT-necessary check cannot fail |

A pass on C1-C5 is a *negative* selection answer plus an open existence
line. It is not a pair.

### Definition of done

| # | Item |
| --- | --- |
| 1 | Solver writes C1-C5 to `data/m9_5_q4a_pair.json` |
| 2 | Auditor recomputes \(b_i\) from a species table, no solver import |
| 3 | Method note and Paper 16. No `MODELS.md`. Q4a existence stays `[O]` |

### Not computed

A CFT. A compactification. Anomaly-inflow matching for a chosen flavor
group. A string vacuum. de Sitter.

## DEVIATIONS LOG

None.

## FINDINGS

Full record:
[`../findings/m9_5_q4a_pair_note.md`](../findings/m9_5_q4a_pair_note.md).
Selection-uniqueness: negative. Existence: still `[O]`.
