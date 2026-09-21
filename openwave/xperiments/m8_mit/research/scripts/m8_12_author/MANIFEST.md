# M8.12 author package: landing manifest

The author's derivation for [M8.12](../../tasks/m8_12_task_details.md), landing after the verdict as the task doc says it would, so the provenance comparison can run. Every file lands byte-identical to the bytes that produced the pre-registration's frozen values and the run-up record on [#571](https://github.com/openwave-labs/openwave/pull/571). With this manifest the directory holds 26 files.

## Where each hash was recorded before this landing

The chain starts at a hash already public: `prereg/s1_prereg_build_log.txt` was posted in #571, and its `(D)` block lists the SHA-256 of sixteen files, fourteen of which land here. The other two are the filed task doc and worklist, which the gate reads from its own folder and which live in the repository at their canonical paths, not duplicated here.

| file | recorded before landing in | landed SHA-256 |
| --- | --- | --- |
| `s1_L0.py` | the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `346b17e3cfaa8e59c5bdf2c6be844f62b6558feef729e54c3dfd6a63e311c979` |
| `s1_L0_log.txt` | the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `cf63a7b795bb28ce12ec7838eeca80520dfe77ef232fce1a07f294133afda353` |
| `s1_L1.py` | the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `30bf03da2b41161b7c4085ca77adcb98a42eafc53e5ea8bde047cde76b73c13a` |
| `s1_L1_log.txt` | the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `b1158cbcd473c0bbdf26e854437cab721ea191bc3481a1bb6189dad24085637b` |
| `s1_H.py` | the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `ffd44322ec22edf4ab61e585b0658b68c414ef65e30fb2d63d2156585d760735` |
| `s1_H_log.txt` | the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `1f55602e43c57028db2ea9f8f491ca2e85130d6da6170060ac1cd1c5a7737297` |
| `s1_control.py` | posted directly in #571, and the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `e1363de6f2f546d9bdef46c7c2029239aa1388865d870a37832ca58e68228f1f` |
| `s1_control_log.txt` | the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `b88e39261339584259c9a17b28dc45b97b3998f50f7e0edee9f56b9202eacdd3` |
| `S1_STEP1_L0.md` | the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `536f63eb15173d537b0aa8d9aa9eaa60fe182b27457c929662b02d9c470dab0f` |
| `S1_STEP2_L1.md` | the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `eb6ecdc070b7c396519d43b09f0987f2db06e1410540da266942fc33591ab186` |
| `S1_STEP3_H.md` | the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `33b149f8cc146833e694e94da39a426e8295cafd60aa94a1d7232a4c82d063c2` |
| `DISCLOSURE_INVENTORY.md` | the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `1ea164c35b1bcd9f4f9c89576082a183376776cbe0f3c9abbce7b9b19071d6c6` |
| `prereg/s1_prereg_build.py` | posted directly in #571, and the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `ece3771f197f8012e87214200509e0b2adf1ac81bbcca27a31ca8638784a5977` |
| `prereg/s1_prereg_arms.py` | posted directly in #571, and the (D) block of `prereg/s1_prereg_build_log.txt`, whose own hash was posted in #571 | `49b8c5261ebdaf6667b973cefc8c1f7e751b83b90d62dc3bb3539e3e7d2256c9` |
| `prereg/gate_inventory.txt` | posted directly in #571, and printed by the gate itself | `c5de76439fa5f9d16dc3c3de1f94b4029f2308569f973d0743f0bd088b374266` |
| `prereg/s1_prereg_build_log.txt` | posted directly in #571 | `229ef1303f6775b3dc19475b9d8cafe7f42d49a1c4dd8938078343848d0096eb` |
| `prereg/s1_prereg_arms_log.txt` | posted directly in #571 | `d3e467d5b4ebadfbeff20bb51029eb2b85ba6b4456bb1696ae8606d963c63374` |
| `S0_S3_MAXIMUM.md` | hard-coded in `prereg/s1_prereg_build.py`, which checks it on every run | `7c634a33fdb0ee05d5a934395a414f6346ba26e049b07f2619961792ab725498` |
| `check_s3_maximum_log.txt` | `DISCLOSURE_INVENTORY.md`, as a 16-digit prefix | `37a00d34168b23200a1f07694627dc2a93dd588b9156be5c55661deffab67b1d` |
| `translate_rho6_log.txt` | `DISCLOSURE_INVENTORY.md`, as a 16-digit prefix | `228e8a35c6c38d1c19c418b620bf491321eb25f235b7e517677207281a13be0a` |
| `check_extrema_log.txt` | `DISCLOSURE_INVENTORY.md`, as a 16-digit prefix | `56ad57c3d54e813d64c7db476acfaf3b009320018dfe15d0238d016a0b525817` |
| `check_s3_maximum.py` | first recorded here; it regenerates its pinned log byte-identically | `34054633f9602e864c12d017b035e6517ebbc4d8b0ba6c043dc1ef4adc764d3a` |
| `translate_rho6.py` | first recorded here; it regenerates its pinned log byte-identically | `f1f254fd0f3f03f90798038528fa040eeec667511fc9818b6c1129df8776b471` |
| `check_extrema.py` | first recorded here; it regenerates its pinned log byte-identically | `86348025d4063973f092d9629f69eed89cf01d3dc2cd72fb9161613c6ac3ac93` |
| `make_disclosure_inventory.py` | first recorded here; its output is pinned, and it cannot rerun here, see below | `fedf9245bb597551e9e26298c71a5dffaa0934cf8c9dcad3be2131eccaf2cfbb` |

Four scripts were never pinned before this landing. Each of the first three regenerates, byte for byte, a log whose hash was pinned through the inventory, which is the strongest provenance available for a script whose own hash was not recorded: the output it produced was, and it still produces exactly that output.

## Regeneration

In a fresh directory holding this folder's files, copy the filed documents into `prereg/` under the names the gate reads, [`../../tasks/m8_12_task_details.md`](../../tasks/m8_12_task_details.md) as `prereg/s1_task_details.md` and [`../../m8_12/worklist.md`](../../m8_12/worklist.md) as `prereg/worklist.md`. Then, each from its own folder:

| command | runtime | prints exactly |
| --- | --- | --- |
| `python3 s1_L0.py` | about 30 s | `s1_L0_log.txt` |
| `python3 s1_L1.py` | about 10 s | `s1_L1_log.txt` |
| `python3 s1_H.py` | about 15 s | `s1_H_log.txt` |
| `python3 s1_control.py` | about 20 s, re-running `s1_H.py` inside it | `s1_control_log.txt` |
| `python3 translate_rho6.py` | about 1 s | `translate_rho6_log.txt` |
| `python3 check_extrema.py` | about 40 s | `check_extrema_log.txt` |
| `python3 check_s3_maximum.py` | about 45 s | `check_s3_maximum_log.txt` |
| `python3 s1_prereg_build.py`, from `prereg/` | a few seconds | `s1_prereg_build_log.txt`, 190 of 190 against the documents at `427e8ce5…` and `c07f9bc6…` |
| `python3 s1_prereg_arms.py`, from `prereg/` | about 10 min | `s1_prereg_arms_log.txt`, 73 of 73 |

All nine were rerun this way before landing, in a directory holding only these files, and all nine match. That holds on the author's machine: Python 3.13.13, sympy 1.14.0, numpy 2.5.0, scipy 1.18.0, mpmath 1.3.0, and pandoc 3.9.0.2, which the gate calls to render both documents as GitHub does. The logs print some floating-point values, so another interpreter or library version can differ in their last digits; the exact values are the ones the scripts identify, not the printed floats.

`make_disclosure_inventory.py` cannot rerun here. It sweeps the author's local exploratory folders, which do not land, and it fails on any record it cannot classify, which is the property that made the inventory complete. Its output, `DISCLOSURE_INVENTORY.md`, is pinned above.

## The two forward items from #571

**The inventory is bound to the documents it pins.** `prereg/gate_inventory.txt` at `c5de76439fa5f9d16dc3c3de1f94b4029f2308569f973d0743f0bd088b374266` is the gate-name pin for the task doc at `427e8ce51b63442c975a0ed0cb4da16171e716210fe51f0a20b22d83ee92167a` and the worklist at `c07f9bc64d39ffc99317188e28a91378be18c0fc8902490bd44bed76bba4efdc`. A different pair of documents may legitimately run a different set of gates, so the inventory is only meaningful against this pair.

**The logs do not depend on the hash seed.** The gate's full output hashes identically under `PYTHONHASHSEED` 0, 1, 2 and 3 and under a random seed, so no ordering in any log derives from a set. Setting `PYTHONHASHSEED=0` makes that hold by construction rather than by observation, and costs nothing.

## One argument this package carries that the run did not settle

The blind run reproduced G2's maximum as an argument, but neither room nor the audit settled whether any orbit other than the hexagon attains it ([method note](../../findings/m8_12_method_note.md), § 5.4). That uniqueness is argued in `S0_S3_MAXIMUM.md` from the equality cases of its three bounds: step 2 carries two of them, both immediate, and step 3 carries the substantive one, the tracing through the three 2×2 blocks that shows `λ_max = 15/√6` only at the three permutations of `(2, −1, −1)/√6`, and hence `TrN̄² = 171/2` only on `span{|3,3⟩ₙ, |3,−3⟩ₙ}`. Step 4 combines them. `check_s3_maximum.py` corroborates step 3's equality set on a 20001-point sweep of the circle and checks step 4 on the computed objects; the exhaustiveness of step 3's tracing is the argument an audit most needs to check. It lands here so that the argument can be audited; until it is, G2's uniqueness clause rests on the author's argument and not on the blind run.

## What does not land

The S0 literature memo and the author's other exploratory records. The inventory names 57 records swept from the two exploratory folders, with their hashes. The S0 memo is not among them, because the sweep excluded `.md` and `.py` files, and six references to it survive in the landed notes and in the gate, which checks that the task doc names it among the firewalled sources. It is available on request.

## Naming

These files keep the `s1_` prefix of the author's working names rather than `m812_`. The filed task doc names `s1_L0.py`, `s1_L1.py`, `s1_H.py` and `s1_control.py` in its feasibility table, the posted hashes are for those names, and the scripts refer to one another by them, so renaming would break both the filed references and the pins. The folder carries the task id.
