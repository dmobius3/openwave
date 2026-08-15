# M9.6: de Sitter at the FGHMV standard (sign and isometries)

> Spec of record: [`../m9_theory_canonical.md`](../m9_theory_canonical.md).
> Paper IX checklist is the closure bar. This task decides Q2 at that
> bar. It does not invent a dual and does not promote Jacobson.

## TASK PLANNING (2026-08-15)

### Scope

Decide whether the AdS bookkeeping can be copied onto de Sitter at
the FGHMV / CHM standard. Two identities are in scope: the
Gibbons-Hawking entropy as a function of \(\Lambda\), and the
dimension of \(\mathfrak{so}(1,4)\) versus \(\mathfrak{so}(2,4)\).

Out of scope: a dS/CFT pair, Hehl-Datta on a horizon, Jacobson
equilibrium as a substitute \pP{}.

### Pre-registered claims

| ID | Claim | Pass | Fail |
| --- | --- | --- | --- |
| C1 | \(S=\pi\ell^2/G=3\pi/(G\Lambda)\) and \(\partial S/\partial\Lambda<0\) | identities hold | the derivative nonnegative |
| C2 | \(\dim\mathfrak{so}(1,4)=10\), \(\dim\mathfrak{so}(2,4)=15\) | exact integers | any other count |
| C3 | CHM kernel of \(T_{00}\) is positive | cited formula | a minus in CHM |
| C4 | C1 and C3 have opposite implications for \(\delta S/\delta E\) | opposite signs | the same sign |
| C5 (mutation) | AdS branch \(\Lambda=-3/\ell^2\) flips \(\partial S/\partial\Lambda\) | positive on that branch | tautology |

A pass is an *obstruction*, not a cosmological Einstein theorem.

### Definition of done

| # | Item |
| --- | --- |
| 1 | Solver writes C1-C5 to `data/m9_6_ds_sign.json` |
| 2 | Auditor recovers \(S(\Lambda)\) and the Lie-algebra dimensions with no solver import |
| 3 | Method note and Paper 17. No `MODELS.md`. No invented dual |

## DEVIATIONS LOG

None.

## FINDINGS

Full record:
[`../findings/m9_6_ds_closure_note.md`](../findings/m9_6_ds_closure_note.md).
FGHMV-standard Q2: obstructed. Einstein+\(\Lambda\) from a
cosmological CFT: still not \pP{}.
