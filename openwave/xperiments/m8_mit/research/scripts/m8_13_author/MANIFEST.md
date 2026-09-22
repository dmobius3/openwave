# M8.13 author package: landing manifest

The author's withheld instruments for [M8.13](../../tasks/m8_13_task_details.md), landing after the verdict as the pre-registration said they would, so the provenance comparison can run. The 13 carried files land byte-identical to pinned bytes, each pinned before the stage it governed: seven at [#583](https://github.com/openwave-labs/openwave/pull/583) before the go, and six at [#585](https://github.com/openwave-labs/openwave/pull/585) after stage 1 had begun and before any stage-2 file was handed over. `MANIFEST.md` is new landing metadata. With it the directory holds 14 files.

## Where each hash was recorded before this landing

| file | recorded before landing in | landed SHA-256 |
| --- | --- | --- |
| `m813_equality.py` | the table in the maintainer's review at [#583](https://github.com/openwave-labs/openwave/pull/583#pullrequestreview-5278860214), which the author cannot edit, and the amendment table in the task doc, merged at `42e65f79` as unchanged | `3577872c112cefe13413eab3a706e834314e635129f7a914f5488d43f3cc090d` |
| `m813_equality_log.txt` | the table in the maintainer's review at [#583](https://github.com/openwave-labs/openwave/pull/583#pullrequestreview-5278860214), which the author cannot edit, and the amendment table in the task doc, merged at `42e65f79` as unchanged | `a30b8198d329d2a038ea3253fce2253b007954db0581fa4b71aeac5d12de0354` |
| `prereg/m813_prereg_build.py` | the amendment table in the task doc, merged at `42e65f79`, superseding the review's `2a09b63d237fd00ce6268c1081602addc03308d047bd473e9a18bba21d16dee4` | `eb0a440887bb778739464d6fed7bcc350038170f024d7820840a8d1480b66dd6` |
| `prereg/m813_prereg_arms.py` | the amendment table in the task doc, merged at `42e65f79`, superseding the review's `077ca28b1c8cda35227a8709a81703031713050cc9204932311a66adf985acda` | `f068c3d27eac293f424df7069d8559c35e8ee8e224964b9e62938a7206c128b0` |
| `prereg/gate_inventory.txt` | the amendment table in the task doc, merged at `42e65f79`, superseding the review's `062fae9d62dba00760fb78a0b89401242e78acb05959e709ad138f4f5dcf7b53` | `71e9a8887b077c382af0290cfc191f2f455298fb2284326e378619e813b78a6c` |
| `prereg/m813_prereg_build_log.txt` | the body of [#585](https://github.com/openwave-labs/openwave/pull/585), the reply on #583, and the commit message of `f7a4d357`, superseding the review's `ea62519dc96599c4f237b93acdf4ed57d5059d43b91054e24a8df56898fe444d` | `c435104a501439e582cf17846d7abc0e5caec7084d8a55baab56ebe56a57aee8` |
| `prereg/m813_prereg_arms_log.txt` | the body of #585, the reply on #583, and the commit message of `f7a4d357`, superseding the review's `b1d281bf69b3059ffda5f861ba788d12f5205c60c6845d8e67069e548ee4da49` | `046b10f8d1fef0680015e5cbe39ef70f38ed106bc47d23718501f86d58a1c21c` |
| `prereg/m8_13_task_details_as_qualified.md` | the file table in the body of #585 | `ecbddd2cf33c64f633c02c1bca11fb5e47fce1cf40149061ffe2bd5770286e3b` |

The five files the #585 amendment superseded land under `prereg/superseded_at_585/`, at the hashes the #583 review pinned, so every hash that review published resolves to a file here:

| file | SHA-256 at #583 |
| --- | --- |
| `prereg/superseded_at_585/m813_prereg_build.py` | `2a09b63d237fd00ce6268c1081602addc03308d047bd473e9a18bba21d16dee4` |
| `prereg/superseded_at_585/m813_prereg_arms.py` | `077ca28b1c8cda35227a8709a81703031713050cc9204932311a66adf985acda` |
| `prereg/superseded_at_585/gate_inventory.txt` | `062fae9d62dba00760fb78a0b89401242e78acb05959e709ad138f4f5dcf7b53` |
| `prereg/superseded_at_585/m813_prereg_build_log.txt` | `ea62519dc96599c4f237b93acdf4ed57d5059d43b91054e24a8df56898fe444d` |
| `prereg/superseded_at_585/m813_prereg_arms_log.txt` | `b1d281bf69b3059ffda5f861ba788d12f5205c60c6845d8e67069e548ee4da49` |

## Why the task doc lands here

The gate reads the task doc and prints its hash, so its logs are bytes of one particular version: the amendment as the author qualified it, `ecbddd2cf33c64f633c02c1bca11fb5e47fce1cf40149061ffe2bd5770286e3b`, the file table of #585's body. That version is not in `main`. #585 was merged as a squash, and the maintainer's note on the amendment's position against the go was added on the branch before the merge, so `main` carries the amended task doc and not the qualified one. The two differ by that note alone. Run against the task doc now on `main`, the gate passes 95 of its 97 checks: the two failures are the required-headings check and the amendment-record check, both keyed to the heading the note renamed. That state is final by design rather than debt: the gate's job was to qualify the pre-registration before the run, which it did against the bytes the rooms were scored on, and those bytes are here. A gate kept passing against a document the run has since appended to would be measuring the wrong thing. The qualified bytes land here so that the regeneration below needs nothing outside `main`.

The commit those bytes came from, `f7a4d357`, survives at `refs/pull/585/head`, whose parent it is, and on the author's fork. Its message carries the two log hashes. Neither is needed to run the route below.

## Regeneration

From this folder, in a clone that has `main`:

```bash
D=$(mktemp -d)
R=openwave/xperiments/m8_mit/research
cp m813_equality.py m813_equality_log.txt "$D"/
cp prereg/m813_prereg_build.py prereg/m813_prereg_arms.py prereg/gate_inventory.txt "$D"/
cp prereg/m813_prereg_build_log.txt prereg/m813_prereg_arms_log.txt "$D"/
cp prereg/m8_13_task_details_as_qualified.md "$D"/m8_13_task_details.md
for f in worklist.md stage2a_grading.md stage2b_grading.md step3_text.md; do git show 42e65f79:$R/m8_13/$f > "$D"/$f; done
for f in S0_S3_MAXIMUM.md check_s3_maximum.py check_s3_maximum_log.txt; do git show 42e65f79:$R/scripts/m8_12_author/$f > "$D"/$f; done
git show 42e65f79:$R/m8_12/worklist.md > "$D"/m8_12_worklist.md
cd "$D" && PYTHONHASHSEED=0 python3 m813_prereg_build.py 2>&1 | diff - m813_prereg_build_log.txt \
        && PYTHONHASHSEED=0 python3 m813_prereg_arms.py 2>&1 | diff - m813_prereg_arms_log.txt
```

Both `diff` commands print nothing. The gate takes about 2 seconds and prints `ALL PASS` on 97 checks; the mutation suite takes about a minute and catches 66 of 66 arms against a green parent.

The five superseded instruments regenerate the same way against the documents as #583 merged them, at `1ffe3a5c`:

```bash
D=$(mktemp -d)
R=openwave/xperiments/m8_mit/research
cp m813_equality.py m813_equality_log.txt "$D"/
cp prereg/superseded_at_585/* "$D"/
git show 1ffe3a5c:$R/tasks/m8_13_task_details.md > "$D"/m8_13_task_details.md
for f in worklist.md stage2_grading.md n2_variant_step3.md; do git show 1ffe3a5c:$R/m8_13/$f > "$D"/$f; done
for f in S0_S3_MAXIMUM.md check_s3_maximum.py check_s3_maximum_log.txt; do git show 1ffe3a5c:$R/scripts/m8_12_author/$f > "$D"/$f; done
git show 1ffe3a5c:$R/m8_12/worklist.md > "$D"/m8_12_worklist.md
cd "$D" && PYTHONHASHSEED=0 python3 m813_prereg_build.py 2>&1 | diff - m813_prereg_build_log.txt \
        && PYTHONHASHSEED=0 python3 m813_prereg_arms.py 2>&1 | diff - m813_prereg_arms_log.txt
```

That gate prints `ALL PASS` on 83 checks and the suite catches 53 of 53.

Both routes were run this way before landing, in a fresh directory holding only these files and the documents the commands name. That holds on the author's machine: Python 3.13.13, sympy 1.14.0, numpy 2.5.0, mpmath 1.3.0, and pandoc 3.9.0.2, which the gate calls to render the documents as GitHub does. `PYTHONHASHSEED=0` makes the ordering in the logs hold by construction rather than by observation.

## What the instruments are

`m813_equality.py` determines step 3's equality set exactly: for each block the equality condition is a conic whose resultant with the ellipse is not zero, so Bezout bounds their common points by four, exact factoring over the field of the square root of 6 enumerates them, and whether the bound is the block's larger eigenvalue there is decided symbolically. It returns the three points and nothing else. It was withheld from the auditor at every stage, so the audit's route is its own; [#586](https://github.com/openwave-labs/openwave/pull/586) records four independent derivations of the same set, of which this is one.

`prereg/m813_prereg_build.py` gates the pre-registration itself: it fails closed on a missing input before reading anything, reruns `m813_equality.py` and requires its output to equal the log byte for byte, checks that every hash the task doc states is the hash of a file it names, pins each frozen value cell and traces its quantities to the argument or to the check, holds the stage-1 worklist against leakage with a list of forbidden values and a pinned count of every numeral, and pins its own check names as a multiset that fails closed when the pin is absent. `prereg/m813_prereg_arms.py` is its mutation suite: each arm plants one defect and requires the specific check it targets to fail, reporting a wrong check or a broken build as a failure.

## These bytes are exempt from formatting

Every carried file in this folder, which is every file except this landing manifest, is pinned by a hash published before the stage it governed, so a formatter run, a line-wrap pass or a lint fix would break the only thing tying the folder to the #583 and #585 records. Exclude it from any repository-wide formatting sweep.
