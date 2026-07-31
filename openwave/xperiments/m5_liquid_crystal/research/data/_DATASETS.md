# Local-only datasets manifest

> AUTO-GENERATED, do not hand-edit the table: `python3 dev_docs/utils/gen_datasets_manifest.py ../data --write`

Heavy binary arrays in this folder are **local-only**: gitignored, never deleted (policy 2026-07-20, which supersedes the earlier "delete raw data > 1 MB" rule). They stay on the working machine so later tasks can consume them directly, and they stay OUT of the repo so clones stay light. What IS tracked in git and readable on GitHub: the summary `.json` / `.csv` / `.txt` in this same folder, the plots, and the scripts that rebuild everything here.

**Inventory**: 76 local-only files, 611.10 MB, in 6 task groups.

| Task group | Files | Size | Producing script(s) | Record (regen commands + context) |
| --- | --- | --- | --- | --- |
| `m5_21_10` | 3 | 293.47 MB | [`m5_21_10_a_decay64.py`](../scripts/m5_21_10_a_decay64.py) · [`m5_21_10_b_ring.py`](../scripts/m5_21_10_b_ring.py) · [`m5_21_10_c_panel.py`](../scripts/m5_21_10_c_panel.py) (+1 more) | [`m5_21_10_task_details.md`](../tasks/m5_21_10_task_details.md) |
| `m5_21_4` | 20 | 67.54 MB | [`m5_21_4_a_pair.py`](../scripts/m5_21_4_a_pair.py) · [`m5_21_4_audit_check.py`](../scripts/m5_21_4_audit_check.py) · [`m5_21_4_c_films.py`](../scripts/m5_21_4_c_films.py) (+2 more) | [`m5_21_4_task_details.md`](../tasks/m5_21_4_task_details.md) |
| `m5_21_5` | 3 | 3.13 MB | [`m5_21_5_a_mu.py`](../scripts/m5_21_5_a_mu.py) · [`m5_21_5_b_ladder.py`](../scripts/m5_21_5_b_ladder.py) · [`m5_21_5_c_bridge.py`](../scripts/m5_21_5_c_bridge.py) (+2 more) | [`m5_21_5_task_details.md`](../tasks/m5_21_5_task_details.md) |
| `m5_21_6` | 6 | 166.69 MB | [`m5_21_6_a_decay.py`](../scripts/m5_21_6_a_decay.py) · [`m5_21_6_audit_check.py`](../scripts/m5_21_6_audit_check.py) · [`m5_21_6_c_films.py`](../scripts/m5_21_6_c_films.py) (+1 more) | [`m5_21_6_task_details.md`](../tasks/m5_21_6_task_details.md) |
| `m5_21_9` | 8 | 7.18 MB | [`m5_21_9_a_audit.py`](../scripts/m5_21_9_a_audit.py) · [`m5_21_9_a_negdelta.py`](../scripts/m5_21_9_a_negdelta.py) · [`m5_21_9_b_audit.py`](../scripts/m5_21_9_b_audit.py) (+6 more) | [`m5_21_9_task_details.md`](../tasks/m5_21_9_task_details.md) |
| `m5_22` | 36 | 73.09 MB | [`m5_22_a_seeds.py`](../scripts/m5_22_a_seeds.py) · [`m5_22_b_census.py`](../scripts/m5_22_b_census.py) · [`m5_22_c_rank.py`](../scripts/m5_22_c_rank.py) (+2 more) | [`m5_22_task_details.md`](../tasks/m5_22_task_details.md) |

**Regeneration**: the exact command + runtime per dataset lives in the task record linked on its row (the task_details / findings doc), which is where the run configuration is already written down. Runs are deterministic from their fixed seeds and configs, so a regenerated array reproduces the original bit-for-bit at the stored precision.
