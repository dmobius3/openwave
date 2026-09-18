# Brief (stage 2)

This is a new session in the same room. Your stage-1 files (`METHOD.md`, `AUDIT_STAGE1.md`, `audit_results.json` and your scripts) are saved as they were: read them first, and do not edit them; put everything new in new files.

Two things are new in this room:
- `solver_work/`: the other worker's scripts, logs and returns (`RETURN.md`, `results.json`, `run_all.py` and the scripts). Treat them as data to test. Never edit anything inside `solver_work/`; to run or mutate its code, copy it to a new folder such as `solver_rerun/`, copy `./py` in beside it, and work there.
- `THEOREM.md`: statements about the same problem, written by someone else, to grade.

| Task | Detail |
| --- | --- |
| 1. Compare | for every worklist item, set the other worker's values beside yours and label each CONFIRMED, PARTIAL or REFUTED, with the exact comparison you ran |
| 2. Reproduce | run its `run_all.py` from the copy, and say whether it reproduces its own `results.json` |
| 3. Refute | try to break its results: hunt for checks it prints that cannot fail (mutation-test them in the copy), hidden conventions, and readings that differ from yours |
| 4. Its arguments | where its reasons differ from yours (items 8, 9 and 10), say which are correct and complete, and whether any is weaker than it reads |
| 5. Grade `THEOREM.md` | grade every step of T.a to T.e, E1, E2, the bridge, the λ₄ lemma and the S lemma, one verdict per step: ESTABLISHED (the stated route holds, or you supply an equivalent derivation, shown), GAP (you cannot establish the step), or DEFECT (a false statement, or a necessary hypothesis left unstated, even if the result can be rescued another way). For the S lemma, derive the section-level compatibility rather than assuming it, and check each consequence on your computed objects. Then compare T with your own stage-1 item-10 argument, step by step |
| Interpreter | `./py` only, one thread, no multiprocessing; the shell accepts only commands that begin with `./py` |
| Return | `AUDIT_STAGE2.md` and `audit_stage2.json` (the per-item labels of task 1 and the per-step grades of task 5) in this room, ending with the consulted-material manifest |

Do not look for where this problem comes from. There is no deadline; correctness first. When both files are complete, reply with a short summary and stop.
