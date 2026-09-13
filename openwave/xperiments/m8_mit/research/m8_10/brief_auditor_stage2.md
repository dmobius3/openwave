# Brief (stage 2)

Your stage-1 files are saved as they were. Do not edit `METHOD.md`, `AUDIT_STAGE1.md` or `audit_results.json`; put everything new in new files.

The other worker's scripts, logs and returns are now in `solver_work/` (`RETURN.md`, `results.json`, `run_all.sh` and the scripts). Treat them as data to test. Never edit anything inside `solver_work/`; to run or mutate its code, copy it to a new folder such as `solver_rerun/`, copy `./py` in beside it, and work there.

| Task | Detail |
| --- | --- |
| 1. Compare | for every worklist item, set the other worker's values beside yours and label each CONFIRMED, PARTIAL or REFUTED, with the exact comparison you ran |
| 2. Reproduce | run its `run_all.sh` from the copy, and say whether it reproduces its own `results.json` |
| 3. Refute | try to break its results: hunt for checks it prints that cannot fail (mutation-test them in the copy), hidden conventions, and readings that differ from yours |
| 4. Arguments | where its reasons differ from yours (item 5's reasons for the zeros, item 7's sign argument), say which are correct and complete, and whether any is weaker than it reads |
| Interpreter | `./py` only, one thread, no multiprocessing |
| Return | `AUDIT_STAGE2.md` and `audit_stage2.json` in this room, ending with the consulted-material manifest |
