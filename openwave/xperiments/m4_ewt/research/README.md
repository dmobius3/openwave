# M4 research: the headless solver

> This folder holds M4's **headless** work: scripts that compute a number, the data
> and plots they produce, and the notes that record what the number means. It is the
> companion to the **rendered** production stack in the model root
> ([`../_launcher.py`](../_launcher.py), [`../wave_engine.py`](../wave_engine.py),
> [`../xparameters/`](../xparameters/)).
>
> The layout mirrors the other model columns ([M5](../../m5_liquid_crystal/research/),
> [M7](../../m7_hydroboros/research/), [M8](../../m8_mit/research/)).

## 0. The documents in this folder

| Doc | What it is |
| --- | --- |
| [`m4_agent_orientation.md`](m4_agent_orientation.md) | the AI-agent front door: the reading list, how tasks run, the platform contract, M4-specific cautions. Say "read the m4_agent_orientation.md" to bootstrap an agent on this column |
| [`m4_theory_canonical.md`](m4_theory_canonical.md) | the spec of record. A skeleton: its § 0 maps where M4's theory record currently lives, and its OPEN QUESTIONS table is seeded from the briefing |
| [`m4_roadmap.md`](m4_roadmap.md) | the program. A skeleton with no tasks by design; its `CONVENTIONS` section defines IDs, row shape, and how to create a task |
| [`M4_engine_upgrade.md`](M4_engine_upgrade.md) | the implementation record for the vector-PDE engine, phase by phase with validation evidence |
| [`M4_k_selectivity_Formalization.md`](M4_k_selectivity_Formalization.md) | the contributed nonlinear-stabilisation formalization |

The two skeletons are empty on purpose. M4's program belongs to the people extending the
model, and a maintainer inventing tasks or paraphrasing theory into them would make the
column less trustworthy, not more.

## 1. The two solvers

OpenWave runs the same field equations two ways. They answer different questions, and
the words for each are used interchangeably across the repo and in discussion:

| | Headless solver | Rendered solver |
| --- | --- | --- |
| Also called | sandbox, research, headless | production, launcher, rendering, GGUI |
| Lives in | this folder, `research/` | the model root: `_launcher.py`, `wave_engine.py`, `xparameters/` |
| How it runs | `python path/to/script.py`, unattended, loopable over a parameter grid | `openwave -x`, one interactive window at a time |
| Output | numbers, JSON/CSV, matplotlib PNGs, a written note | on-screen 3D render, glyphs, flux meshes, exportable GIF/PNG |
| Answers | "which configuration is selected, and by how much" | "what does this mechanism look like, and how does it respond when perturbed" |
| Good for | quantitative validation, reproducibility, the numbers behind [`MODELS.md`](../../../../MODELS.md), AI-agent runs, CI | seeing a mechanism, teaching, exploratory cause-and-effect once the physics is known |

Full description of the split: [`TUTORIAL.md § 3`](../../../../TUTORIAL.md#3-the-two-solvers-headless-sandbox-vs-rendered).

## 2. The order they run in

**Headless first. The rendering port is the last step, not the starting point.**

A parameter study belongs in the headless solver because that is the only place its
result survives the run: a script plus a parameter loop runs the whole grid unattended
and ends in a table anyone can regenerate. The same study driven through the rendered
launcher is one window per configuration, and its conclusion lives with whoever watched
them. Rendering an unvalidated dynamics also showcases nothing, since there is not yet a
statement about what the viewer is looking at.

This is the convention across the platform, not an M4 rule:

| Column | Where it is written |
| --- | --- |
| M5 reached its `_launcher.py` through headless rounds first, and holds "rendering gates nothing" as a standing convention | [`m5_visualization.md`](../../m5_liquid_crystal/research/m5_visualization.md) |
| M7: "Headless first ... rendering graduates once the electron is canonical ... No GUI / viz work before the physics is canonical" | [`m7_roadmap.md`](../../m7_hydroboros/research/m7_roadmap.md) |
| M8 gates its rendering port on validated dynamics outright: "rendering an unvalidated dynamics showcases nothing" | [`m8_roadmap.md`](../../m8_mit/research/m8_roadmap.md) (task M8.7) |
| The port path itself, written for AI agents, opens with "Do not start here" | [`m8_platform_pointers.md § 7`](../../m8_mit/research/m8_platform_pointers.md) |

When a result does graduate to the rendered solver, the launcher runs the **same kernels
the headless work validated**. Interactive demos show the physics of record, including
its instabilities.

## 3. Folder layout

| Folder | Holds | Tracked in git |
| --- | --- | --- |
| `scripts/` | headless scripts, one per study or a runner plus a module | yes |
| `data/` | summary JSON / CSV / TXT the scripts emit, plus [`_DATASETS.md`](data/_DATASETS.md) | yes, except heavy binary arrays (`.npz`, `.npy`, `.h5`, ...) which stay local |
| `plots/` | matplotlib PNGs | yes |
| `tasks/` | one `m4_<id>_task_details.md` per study: plan, method, findings | yes |
| `findings/` | method notes and standalone result write-ups | yes |
| `checkpoints/` | in-flight progress state | no, gitignored repo-wide |

Name artifacts with the task id as a prefix (`m4_<n>_...`, matching the roadmap ID) so a
script, its data and its plots stay grouped for good.

Nothing that opens a GGUI window belongs under `research/`, whatever else it does. The
standalone harmonic-motion demo [`granule_motion.py`](../granule_motion.py) sits in the
model root for that reason: it renders, so it is a rendered artifact.

## 4. What counts as a finished result

A run is not a result. A result is a written statement backed by a script anyone can
rerun. The reporting standard for anything substantive is
[`dev_docs/METHOD_NOTE.md`](../../../../dev_docs/METHOD_NOTE.md): equations first, a map
from each equation to the code that implements it, figures embedded inline, an explicit
list of what was not computed, and an adversarial pass that tries to refute the claim.
The working contract for AI-assisted work is
[`AI_HYGIENE.md`](../../../../AI_HYGIENE.md): model output is a draft until something
that is not a language model confirms it.

A documented negative meets this bar exactly as well as a positive. "This does not work,
here is the script, here is why" is a contribution, and it is the form that lets the next
person skip the same dead end.
