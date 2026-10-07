# M8.15: the freeze block

The items the pre-registration fixes at the freeze (§ FIXED AT THE FREEZE of [the task](../tasks/m8_15_task_details.md)), recorded on 2026-10-07 before any target was run. Every hash below was computed from the bytes it names.

## The terms

- **SHA-256 of [`tasks/m8_15_task_details.md`](../tasks/m8_15_task_details.md):** `b21eed1a487dea27431ecc3336417fa558df3b9b848fcce6ba3def7dafcbe5e4`.
- The same hash is on the author's parent page, MIT [`soft-slot-branches.md`](https://github.com/dmobius3/mode-identity-theory/blob/ae1286ace4a6be47f5637425eb0c4ba58806e2b9/files/framework/files/working/files/soft-slot-branches.md) § VI, at `ae1286a` on `main`, pushed on 2026-10-07 before any target was run. The page also keeps the terms' first hash, taken the same day before a rendering fix to one table row that changes no term.

## The solver

- **The commit:** MIT [`6058f95`](https://github.com/dmobius3/mode-identity-theory/tree/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics) (`6058f95c3a893652d8fb1db19e7d2d14060ba881`), folder `files/framework/files/working/files/scripts/soft-slot-dynamics/`, unchanged at `ae1286a`. It carries the reproduction route [`REPRODUCE.md`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/REPRODUCE.md) and [`SHA256SUMS`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/SHA256SUMS), which pins every file in the folder.
- **The manifest,** [`m8_15_solver/out/MANIFEST.json`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/out/MANIFEST.json): the SHA-256 of every file the run executes, the outside-orbit check with its imports and point copy, and the frozen numbers. Paths are relative to `m8_15_solver/`.

  | file | SHA-256 |
  | --- | --- |
  | `geometry.py` | `788e10eb6973668e84370c13a8bd7b16dc9e712f6fa681c7386e15e9b1d1abff` |
  | `irreps.py` | `cbc7f4481706a08d7ca0229914980e773125c642da616e867580fe1620b5a149` |
  | `fem.py` | `209fba67eecf8c2b6ecaa432aafda3c6a161024d83d14a060809b81ca2e1061d` |
  | `nonlinear.py` | `f00ca41c1406d123dac4384cd2f7336410be49e851426db15cad79da3102a060` |
  | `standing.py` | `455cfc94f4bb10b97eaf0c426c86cb933e32fd7b199b00fc153bcad88fde5efb` |
  | `linearize.py` | `7cca177bf0c8efc344cc76d2cbcee94f0a96987ae9f3c0895ae64f2f705eed58` |
  | `integrate.py` | `3c137e773f7890478b2495d6b466dc167c2aeffc4a8d179a0b7775b9b6b58096` |
  | `symmetry.py` | `c6bbb178dd9103ab5b2f168c6df8cd923b4c57bd8432215b1b9d4a5bef5aae19` |
  | `pinned.py` | `d4d9606d620986ae562063d6275fcb7c78df569d725f0607b89a0237b56f5b31` |
  | `continuation.py` | `f0117f59de69376ac130017005d9d5ae9d7bccd356ba7ac018ead14777b6cbcb` |
  | `verdict.py` | `17dcd3571480abd99a58c6c75903e911e0a87343e9fe138318c9537538f29864` |
  | `persistence.py` | `49d9856961b39b0c5077951c70e94c9325f2eae83f6cb863053c574053a12946` |
  | `m815_common.py` | `dc48c8aa5da834dc3aaf43fffda060bf08c5ef83cfc19dfa431cd48f6805aa3a` |
  | `m815_controls.py` | `dfe12487e5a8c495aeb7d30121fcd63e5d2e70f7c7b0154f69fdda76d95b783c` |
  | `m815_branch.py` | `0693617814bc28dc2549ebc017710370505eeb3e61035de845a44c26e8891c7e` |
  | `m815_freeze.py` | `6a47bf4376aadd589243408b6d925c1910c649080dd1df85075fa5135320a866` |
  | `m8_15_math/outside3_krein.py` | `94ba627cfc455d20cae34ef761f191bcfbf89e170748545d30f72d3a8b7d89e5` |
  | `m8_15_math/outside3_krein.out` | `ad01a27c993793eac968ec7ccba25499f15eaf42b10ff58234ea7079313a42b1` |
  | `m8_15_math/analyse.py` | `971f1a755cc9b3b329552e299f2e6f09a8851121b3e975e30015ffac9242c5a0` |
  | `m8_15_math/core.py` | `3507004f47d3bcc248012335e4e20e70c47cae118e9aa2a87cf85f9cc62f86b3` |
  | `m8_15_math/outside3_p.npy` | `027a797c805be453c8ea0cc1ca8c9eeb4d276d7e765a7979e1daf09954e82f0d` |
  | `out/freeze_numbers.json` | `c46b102c5ded4a78667dfb0e87e1311d729b006cb1362c6aeb91cbe1b9fe369f` |
  | `out/orientations.json` | `21b1f19161058da22a074acdc12150c8650ce60c7cdc9a86bc6622fc452290a1` |

- **The environment** the controls ran in, as the manifest records it: Python 3.13.13 (conda-forge), numpy 2.5.0, scipy 1.18.0, BLAS from Accelerate, macOS-15.8-arm64-arm-64bit-Mach-O, one thread for every BLAS (`OMP`, `OPENBLAS`, `VECLIB` and `MKL` thread counts all 1), ARPACK starting vectors from the seed `271828`. The solver is deterministic under the frozen environment; bit-for-bit agreement on other machines is not promised.

## The random seed

`20261006`: a fresh `numpy.random.default_rng(20261006)` for each persistence test (§ THE NONLINEAR PERSISTENCE VERDICT).

## The setup of every case, at both levels

[`m8_15_solver/out/orientations.json`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/out/orientations.json), hashed in the manifest, holds each case's orientation `g₀`, its pinned subgroup as indices into the mesh's right 2I, its fibre maps and its generators, at levels 2 and 3. The fibre fit is the residual of the fibre maps' fit; the seed's `Q` is the germ's at that level, against the Result's (§ FROZEN VALUES).

| case | slot | order of `H` | fixed dimension (levels 2, 3) | fibre fit, worst | seed `Q` (levels 2, 3) | the Result's `Q` | relative difference (levels 2, 3) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| C1 | `R3` | 10 | 1216, 14528 | 1.1e-09 | 1.0000000, 1.0000000 | not tabled | not tabled |
| C2a | `R4` | 10 | 1216, 14528 | 2.7e-09 | 1.0007811, 1.0007772 | 1288/1287 | 4.1e-06, 2.0e-07 |
| C2b | `R2` | 10 | 816, 9696 | 3.1e-09 | 1.0069691, 1.0069917 | 144/143 | 2.4e-05, 1.3e-06 |
| C3 | `R4` | 24 | 496, 6032 | 1.2e-09 | 1.2238683, 1.2237806 | 175/143 | 7.5e-05, 3.6e-06 |
| T1 | `R4` | 10 | 1216, 14528 | 2.7e-09 | 1.0279289, 1.0279699 | 147/143 | 4.2e-05, 2.1e-06 |
| T1 | `R5` | 10 | 1616, 19360 | 2.7e-09 | 1.0157480, 1.0157350 | 581/572 | 1.3e-05, 7.7e-07 |
| T2 | `R4` | 12 | 1024, 12128 | 2.7e-09 | 1.1589908, 1.1590137 | 5831/5031 | 2.0e-05, 3.8e-07 |
| T2 | `R5` | 12 | 1350, 16146 | 2.7e-09 | 1.0894409, 1.0894456 | 609/559 | 4.2e-06, 1.9e-07 |
| T3 | `R4` | 12 | 1024, 12128 | 2.7e-09 | 1.3595589, 1.3597434 | 1750/1287 | 1.4e-04, 5.9e-06 |
| T3 | `R5` | 12 | 1350, 16146 | 2.7e-09 | 1.2023665, 1.2023616 | 2751/2288 | 5.3e-06, 1.2e-06 |
| T4 | `R4` | 4 | 3040, 36320 | 1.8e-10 | 1.1852096, 1.1852039 | 1.185203461338846… | 5.2e-06, 3.8e-07 |
| T4 | `R5` | 4 | 4032, 48384 | 2.6e-10 | 1.1041759, 1.1041770 | 1.104176947003101… | 9.3e-07, 1.9e-08 |
| T5 | `R2` | 20 | 408, 4848 | 3.2e-09 | 1.6917316, 1.6922787 | 22/13 | 3.4e-04, 1.7e-05 |

## The numbers

From [`m8_15_solver/out/freeze_numbers.json`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_solver/out/freeze_numbers.json), hashed in the manifest; § FIXED AT THE FREEZE states them.

- **Floors:** `r_f⁽³⁾ = 1.0397e-06`, `r_f⁽²⁾ = 3.9702e-05`; `τ₀ = τ_Re = 10 r_f` at each level.
- **Level stability:** `d₃ = 1.247e-04`, so `δ_lev = 0.01`.
- **`s₂`:** 6.7391e-04 (`R4`), 9.5472e-04 (`R5`), 2.1362e-03 (`R2`).
- **Persistence:** `η = η_ref = 0.001`, `T = 1000`; level-2 `dt` = 0.008197 (`R2`), 0.008158 (`R3`), 0.008187 (`R4`, `R5`).
- **The unresolved modes:** `T4/R4`, `τ = 0.00573686`: first resolved at `ε = 1.41421`; `T4/R5`, `τ = 0.00573686`: first resolved at `ε = 2`. Every other frozen mode is resolved at `ε_min = 0.25`.
- **The gates:** cluster, linear and dynamic, all pass.

## The third outside orbit's frozen frequencies and Krein signs

[`m8_15_math/outside3_krein.py`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_math/outside3_krein.py) reads the certified point from `branches_check.py` at MIT `d0de8ca` and, at 60 digits, prints for each slow frequency `ν = −τ²` as an eigenvalue of the square of the slow matrix; the three smallest singular values of that square minus `ν`, two of them zero to the working precision, so the pair spans an invariant plane, and the third bounded away; and the quadratic form `M` on that plane, definite, whose sign is the Krein sign. Its record, [`outside3_krein.out`](https://github.com/dmobius3/mode-identity-theory/blob/6058f95c3a893652d8fb1db19e7d2d14060ba881/files/framework/files/working/files/scripts/soft-slot-dynamics/m8_15_math/outside3_krein.out):

```
pinned point equals the bench copy: True
Q(u) = 1.185203461338846   nrm before = 1.0
nu = -1.02432941669733: tau = 1.01209160489; smallest singular values ['2.56e-61', '6.03e-61', '0.892']; M on the pair: ['51.333', '78.4332']
nu = -0.0965299801561095: tau = 0.310692742362; smallest singular values ['1.44e-62', '3.6e-62', '0.0296']; M on the pair: ['-29.0265', '-21.2184']
nu = -0.019391072983931: tau = 0.139251832964; smallest singular values ['1.09e-62', '8.3e-62', '0.0102']; M on the pair: ['-1.35013', '-7.29987']
nu = -3.29115537456027e-5: tau = 0.00573685922309; smallest singular values ['4.72e-63', '2.44e-62', '1.84e-5']; M on the pair: ['-0.00862601', '-0.0415785']
```

## The continuation algorithm

As written in the task: § THE SOLVER for the continuation in `m`, the accepted solutions, the fold, the secant and the Newton retries, and § THE AMPLITUDE LADDER for the grid, the branch end and the bisection. It is implemented in `m8_15_solver/continuation.py`, hashed in the manifest.

## The order

The solver landed on MIT `main` at `6058f95` and the terms-hash line at `ae1286a`, both on 2026-10-07. The pre-registration is this pull request's merge, and no target is solved before it. After the maintainer's go, the guard file `out/FROZEN.json`, which the solver requires before any target runs, is written with the terms' SHA-256.
