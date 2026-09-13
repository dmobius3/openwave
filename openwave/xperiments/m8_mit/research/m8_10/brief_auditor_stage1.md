# Brief (stage 1)

You are an independent implementer. The file `worklist.md` in this room is a conventions extract and a numbered list of questions. Another worker is answering the same list separately; you will not see its work in this stage. Later you will receive it for comparison, and only after your own results are saved.

| Rule | Detail |
| --- | --- |
| Method first | before computing anything, write `METHOD.md` in this room describing the route you will take for each item. Do not edit it afterwards; if the route changes, append a dated note |
| Independence | build your own route. Construct sections as explicit functions on `SU(2)` (unit quaternions) and obtain projections by averaging over the elements of the finite group and by representation theory you implement yourself. Do not call a library's tabulated Clebsch-Gordan routine; if you need coupling coefficients, compute them yourself, and cross-check them by a second construction |
| Interpreter | run Python only as `./py` from inside this room. It is pinned to one thread and low priority because another workload shares the machine; do not use multiprocessing or threads |
| Exactness | where an item asks for an exact value, prefer exact arithmetic. If you identify a rational from high-precision numerics, state the precision and the denominator bound, and repeat at a second precision |
| Checks | mutation-test every line you print as PASS: change what it checks to something wrong and confirm it fails; record the mutation |
| Return | write `AUDIT_STAGE1.md` item by item (values, method, readings), and `audit_results.json` with every reported value |
| Manifest | end `AUDIT_STAGE1.md` and your reply with the consulted-material manifest, including item 1's DERIVED or RECOGNIZED declaration |

Do not look for where this problem comes from. There is no deadline; correctness first.
