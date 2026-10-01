# Brief, stage 1

You are the auditor. Read `auditor_brief.md` and `spec_sheet.md` in this room. This session is **stage 1 only**: do V2, commit it, and stop. Nothing else will be given to you in this session.

| Rule | Detail |
| --- | --- |
| Tools | your shell accepts only commands that begin with `./py` (this room's Python interpreter: Python 3.12 with numpy 2.5.3, scipy 1.18.1, sympy 1.14.0 and mpmath 1.3.0); use Python for anything else, such as listing files. Create and edit files with the file tools or from Python, inside this room only |
| The room | this directory is your whole world, and that is a rule, not a description of the machine. A read outside it is recoverable if you report it; concealing one is not. There is no network |
| Interpreter | `./py` is pinned to one thread and low priority because other work shares the machine; do not use multiprocessing or threads |
| Your discretization | a different one from T4's P1 protocol, as `auditor_brief.md` V2 asks. Run your own convergence study and state your error estimate per width and how it was formed |
| No search | report every number you compute, including ones that look wrong, and tune nothing toward any value |
| Scripts | every number you report comes from a script saved in this room, with paths relative to the script, and one script `run_all.py` runs them all in order through `./py` |
| Checks | any line a script prints as PASS must be able to fail: say how you know it can, and show at least one planted defect firing |
| Readings | if anything is underdetermined, say so and state the reading you took |
| Commit | write `STAGE1.md` (method, convergence evidence, the table of results) and `stage1.json` with every number; it must include a block `"v2"` keyed by the width labels `"1/4"`, `"1/2"`, `"1"`, `"7/5"`, `"3/2"`, `"pi/2"`, `"17/10"`, each holding `"lowest"` (your best values of the lowest six eigenvalues), `"bottom"` and `"err"` (your error estimate on the bottom). These files are your stage-1 commit and are hashed when you stop |
| Manifest | end `STAGE1.md` and your final reply with the consulted-material manifest: every file you read, including anything that loaded without your asking, and every result you used from memory rather than derived (name it) |

Do not look for where this problem comes from. Print relative paths only. There is no deadline; correctness first. When `STAGE1.md` and `stage1.json` are complete, reply with a short summary and stop.
