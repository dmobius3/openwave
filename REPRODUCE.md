# Reproducing OpenWave Results

Every claim published from this repository is backed by a runnable artifact: a script, a dataset manifest, or a research document with the command that produced it. This page is the front door from a **clean clone** to re-running the artifact behind any claim.

> **What this file is, and is not.** This is a STATIC convention guide: it documents the environment and the fixed traversal from a claim to its command. It never logs results. Result-level reproduction commands live in exactly one place, the per-task research docs, and stay there; nothing lands here when a validation closes. See the [maintenance contract](#5-maintenance-contract).

## 1. Clean-clone setup

```bash
# Python >= 3.12 required
git clone https://github.com/openwave-labs/openwave.git
cd openwave
pip install .          # reads dependencies from pyproject.toml

# Sanity check
python3 -c "import openwave; print('openwave imports OK')"
```

| Note | Detail |
| --- | --- |
| Development install (editable mode, DCO hook, LaTeX/FFmpeg extras) | [`CONTRIBUTING.md`](CONTRIBUTING.md) |
| Research scripts are headless | They write plots and data to files via matplotlib; no display or GUI is required to reproduce research results |
| Taichi kernels | Scripts that use Taichi select a backend automatically; CPU works everywhere, GPU accelerates where available |
| Interactive xperiments (GUI) | `openwave -x` launches the CLI xperiment selector; optional, not needed for reproduction |

## 2. From a claim to its command

The repository keeps reproduction information **single-source**: the command that regenerates a result is recorded once, in the research doc of the task that produced it. Everything else links there. The traversal:

| You start from | Route to the command |
| --- | --- |
| A [`MODELS.md`](MODELS.md) coverage-matrix cell | The cell links its backing script or research doc directly; scripts live under `openwave/xperiments/<model>/research/scripts/` |
| A findings note (`research/findings/`) | The note links its task doc (`research/tasks/<taskID>_task_details.md`); task docs record the exact regeneration command(s) for their artifacts at close |
| A data array (`.npz` and similar) | Heavy arrays are regenerable and not tracked in git; a `data/` folder that owns such arrays carries a `_DATASETS.md` manifest listing each array's regeneration script, mode, and runtime |
| A plot (`research/plots/` or embedded in a doc) | Plots are tracked; the producing script shares the same task-id prefix in `research/scripts/` |

The glue is the **task-id naming convention**: within a model, one task id prefixes the task doc, its scripts, its data, and its plots (for example `m5_21_2b_*`), so any artifact resolves to its task doc, and the task doc holds the commands.

A good first reproduction: pick any ✅ cell in [`MODELS.md`](MODELS.md) and run the script it links. Task docs state runtimes where they are non-trivial.

## 3. Per-model entry points

| Model | Briefing | Research index |
| --- | --- | --- |
| M4 EWT | [`__M4_model_briefing.md`](openwave/xperiments/m4_ewt/__M4_model_briefing.md) | [`research/`](openwave/xperiments/m4_ewt/research/) (legacy layout, predates the task-doc convention) |
| M5 Liquid Crystal | [`__M5_model_briefing.md`](openwave/xperiments/m5_liquid_crystal/__M5_model_briefing.md) | [`m5_roadmap.md`](openwave/xperiments/m5_liquid_crystal/research/m5_roadmap.md) |
| M6 Ouroboros | [`__M6_model_briefing.md`](openwave/xperiments/m6_ouroboros/__M6_model_briefing.md) | [`m6_roadmap.md`](openwave/xperiments/m6_ouroboros/research/m6_roadmap.md) |
| M7 HydroBoros | [`__M7_model_briefing.md`](openwave/xperiments/m7_hydroboros/__M7_model_briefing.md) | [`m7_roadmap.md`](openwave/xperiments/m7_hydroboros/research/m7_roadmap.md) |
| M8 MIT | [`__M8_model_briefing.md`](openwave/xperiments/m8_mit/__M8_model_briefing.md) | [`m8_roadmap.md`](openwave/xperiments/m8_mit/research/m8_roadmap.md) |

## 4. What "reproducible" means here

| Bar | Meaning |
| --- | --- |
| Script-backed | No claim in the shared benchmark rests on prose; every [`MODELS.md`](MODELS.md) cell links something runnable or a research doc that does |
| Deterministic where stated | Task docs state seeds and settings; descent-type runs regenerate their endpoints from the recorded commands |
| Negatives included | ❌ cells are results with scripts behind them, reproducible like any other |
| Honest icons | ✅ ⚠️ ❌ 🔶 🚧 status semantics are defined in [`MODELS.md`](MODELS.md) |

The AI-collaboration contract that keeps this auditable is [`AI_HYGIENE.md`](AI_HYGIENE.md); the reporting standard for methods is [`dev_docs/METHOD_NOTE.md`](dev_docs/METHOD_NOTE.md).

## 5. Maintenance contract

| This file changes when | This file NEVER changes when |
| --- | --- |
| The environment or setup steps change | A validation lands or a task closes (that lives in the task docs and [`MODELS.md`](MODELS.md)) |
| A model is added or retired (one row in § 3) | Results, numbers, or statuses update anywhere |
| The repo-wide conventions in § 2 themselves change | A new script, dataset, or plot is added under an existing task |

Keeping this file static is deliberate: the multi-model arena already maintains one cross-model surface ([`MODELS.md`](MODELS.md)), and reproduction commands stay single-source in the task docs. This page only teaches the traversal.

---

Cross-refs: [`README.md`](README.md) (what OpenWave is) · [`CONTRIBUTING.md`](CONTRIBUTING.md) (dev setup + PR flow) · [`MODELS.md`](MODELS.md) (the coverage matrix) · [`ONBOARDING_MODELS.md`](ONBOARDING_MODELS.md) (adding a model) · [`AI_HYGIENE.md`](AI_HYGIENE.md) (the AI-collaboration contract)
