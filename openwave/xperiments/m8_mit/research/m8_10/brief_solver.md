# Brief

You are the solver. The file `worklist.md` in this room is the whole problem: a conventions extract and a numbered list of questions. Work every item, in order, starting with item 0.

| Rule | Detail |
| --- | --- |
| Interpreter | run Python only as `./py` from inside this room. It is pinned to one thread and low priority because another workload shares the machine; do not use multiprocessing or threads |
| Exactness | where an item asks for an exact value, prefer exact arithmetic (rationals, or the golden ratio as an algebraic number). If you identify a rational from high-precision numerics instead, state the precision, the denominator bound, and repeat the identification at a second precision |
| Scripts | every number you report comes from a script saved in this room, with paths relative to the script |
| Checks | any line a script prints as PASS must be able to fail: say how you know it can |
| Readings | if an item is underdetermined, say so and state the reading you took |
| Return | write `RETURN.md` in this room, item by item: values, method, readings, and anything that looked wrong. Also write `results.json` with every reported value |
| Manifest | end `RETURN.md` and your reply with the consulted-material manifest, including item 1's DERIVED or RECOGNIZED declaration |

Do not look for where this problem comes from. There is no deadline; correctness first.
