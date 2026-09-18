# Brief

You are the solver. The file `worklist.md` in this room is the whole problem: a conventions extract and a numbered list of questions. Work every item, in order, starting with item 0.

| Rule | Detail |
| --- | --- |
| Tools | your shell accepts only commands that begin with `./py` (this room's Python interpreter); use Python for anything else, such as listing files. Create and edit files with the file tools or from Python, inside this room only |
| Interpreter | `./py` is pinned to one thread and low priority because another workload shares the machine; do not use multiprocessing or threads |
| Exactness | where an item asks for an exact value, prefer exact arithmetic (rationals, or algebraic numbers such as square roots of rationals). If you identify a value from high-precision numerics instead, state the precision, the denominator bound and the field, and repeat the identification at a second precision |
| Scripts | every number you report comes from a script saved in this room, with paths relative to the script, and one script `run_all.py` runs them all in order through `./py` |
| Checks | any line a script prints as PASS must be able to fail: say how you know it can |
| Readings | if an item is underdetermined, say so and state the reading you took |
| Arguments | items 8, 9 and 10 ask for arguments; give them in full in prose, not only as computed evidence |
| Return | write `RETURN.md` in this room, item by item: values, method, readings, and anything that looked wrong. Also write `results.json` with every reported value, exact values as strings |
| Manifest | end `RETURN.md` and your final reply with the consulted-material manifest, including item 1's DERIVED or RECOGNIZED declaration |

Do not look for where this problem comes from. There is no deadline; correctness first. When `RETURN.md` and `results.json` are complete, reply with a short summary and stop.
