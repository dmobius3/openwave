# M5.29: MODELS.md criteria housekeeping: split the coarse rows, name the simplest passing test per row

**Status**: 🚧 STAGED 2026-07-27 (user decision, same day as the author proposal). Origin: the author's 2026-07-27 reply on the go-pack thread ([`m5_22_convo.md § 2026-07-27 afternoon`](m5_22_convo.md)), two structural asks about the [`MODELS.md`](../../../../../MODELS.md) score-board. Platform-wide housekeeping: the criteria rows are shared by ALL model columns, so every edit here touches the whole matrix, not just M5. Docs-only, no runs.

## The two asks

| Ask | Content | Source |
| --- | --- | --- |
| Split the coarse criteria | Some rows bundle results of very different difficulty, so a single icon is dishonest in both directions. Named examples: **Weak force** contains muon decay (the simplest example, measured for M5) and beta decay (much harder, needs a proper neutron); **Baryons (p, n)** contains why the neutron is heavier, the positive core / negative shell, exact masses, and beta decay | author, 2026-07-27 (summarized; 1:1 channel) |
| Every row names its simplest passing experiment | "There should be specified some simplest e.g. experiment used to test given property." For "chiral SU(2)" the author's own question is what the test would even be. For **Magnetic force**: Coulomb + Lorentz covariance jointly imply magnetism, but **Lorentz covariance is itself missing as a criterion**; the closest direct test is **Larmor precession** (electron magnetic dipole in an external field); ortho- vs para-positronium (annihilation lifetimes, 2 vs 3 photons) is the harder follow-on | same |

## Scope

| Piece | Content |
| --- | --- |
| Row-split proposal | A proposed new criteria set: at minimum Weak force → (muon decay \| neutron beta decay) and Baryons → (mass ordering + core/shell profile \| exact masses \| beta decay), plus a sweep of the other 20 rows for the same bundling smell (candidates on sight: Neutrinos, Gravity, Strong force). Each split row keeps a defensible icon PER MODEL from the evidence already linked, no new runs |
| The test column | A per-row "simplest passing test" entry (design decision: a new column vs a first line of the row label; weigh against the 55-word cell budget and the matrix width). Larmor precession lands as the magnetic-force test; the chiral-SU(2) row gets its test named or the row gets rephrased as testable |
| The Lorentz covariance row | New criterion row, missing today and load-bearing (it is what upgrades Coulomb into magnetism) |
| Score-board integrity | The summary counts (`Total criteria` = 22 today) and every per-model detail table re-derived to match the new row set; `python3 dev_docs/check_models_md.py` green (55-word cells, icon/status sync) |
| Editing rules | The standing MODELS.md rules apply: cells are state of the art, no names/dates/endorsements in cells, verifiable refs only. The author-proposal provenance lives HERE and in the convo record, not in the matrix |

## What this task is not

No physics runs, no new evidence, no icon upgrades. Splitting a row may DOWNGRADE the apparent coverage of a model (a bundled ⚠️ can become one ✅ sub-row and one 🚧 sub-row); that is the point, and the score-board counts move accordingly. Any run the new rows motivate (Larmor precession on the 3×3/4×4 stack) is staged separately.

**Gated by**: user "go" (docs-only; can interleave between physics tasks).

## Cross-links

| Doc | Why |
| --- | --- |
| [`m5_22_convo.md`](m5_22_convo.md) | The 2026-07-27 proposal record + routing |
| [`MODELS.md`](../../../../../MODELS.md) | The target doc |
| [M5.22](m5_22_task_details.md) | The thread this proposal arrived on; its results will feed the split Baryons rows |
| [`../m5_roadmap.md`](../m5_roadmap.md) | Backlog row + the change-log entry |
