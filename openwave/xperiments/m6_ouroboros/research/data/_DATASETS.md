# Local-only datasets manifest

> AUTO-GENERATED, do not hand-edit the table: `python3 ~/Documents/source_code/XRODZ/dotfiles/claude_projects/NEPTUNYA-SABER/scripts/gen_datasets_manifest.py research/data --write`

Heavy binary arrays in this folder are **local-only**: gitignored, never deleted (policy 2026-07-20, which supersedes the earlier "delete raw data > 1 MB" rule). They stay on the working machine so later tasks can consume them directly, and they stay OUT of the repo so clones stay light. What IS tracked in git and readable on GitHub: the summary `.json` / `.csv` / `.txt` in this same folder, the plots, and the scripts that rebuild everything here.

**Inventory**: 4 local-only files, 0.07 MB, in 1 task groups.

| Task group | Files | Size | Producing script(s) | Record (regen commands + context) |
| --- | --- | --- | --- | --- |
| `m6_4` | 4 | 0.07 MB | [`m6_4_audit_algebra.py`](../scripts/m6_4_audit_algebra.py) · [`m6_4_audit_june_census.py`](../scripts/m6_4_audit_june_census.py) · [`m6_4_audit_may_census.py`](../scripts/m6_4_audit_may_census.py) (+11 more) | [`m6_4_task_details.md`](../tasks/m6_4_task_details.md) |

**Regeneration**: the exact command + runtime per dataset lives in the task record linked on its row (the task_details / findings doc), which is where the run configuration is already written down. Runs in both repos are deterministic from their fixed seeds and configs, so a regenerated array reproduces the original bit-for-bit at the stored precision.
