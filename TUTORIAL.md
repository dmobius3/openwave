# OpenWave Tutorial: from zero to your first contribution

## What OpenWave is (read this first if you are new)

OpenWave is an open-source platform for **testing candidate field-theoretic models of matter against observation, uniformly**. It is not one theory with one substrate. It is a **model database**: each candidate framework becomes a **column** scored against a shared set of criteria (particles, forces, waves, quantum emergence) in [`MODELS.md`](MODELS.md), and every cell in that matrix is earned by a **runnable script plus an honest research note**, or it stays marked "not yet tested". The bar is **reproducibility, not orthodoxy**, and a documented negative is as valuable as a positive.

This is the hands-on, start-to-finish guide for a **contributor**: from a clean machine to running a simulation, validating an existing model, and opening a pull request. **If you are a model author wanting to add your own framework as a new column, go to [`ONBOARDING_MODELS.md`](ONBOARDING_MODELS.md)** instead: it is the single doc for that path (self-evaluation, application, scaffold, and what a model author owns).

## TL;DR: If you only read one thing first, read this section: [START HERE](#0-start-here-drive-openwave-with-an-ai-coding-agent)

| You want to | Jump to |
| --- | --- |
| Move fast with an AI agent | [0. START HERE: drive OpenWave with an AI agent](#0-start-here-drive-openwave-with-an-ai-coding-agent) |
| Get the code running | [1. Set up your environment](#1-set-up-your-environment) |
| Know the house rules | [2. Community: conduct and contributing](#2-community-conduct-and-contributing) |
| Understand the two ways to run | [3. The two solvers: headless vs rendered](#3-the-two-solvers-headless-sandbox-vs-rendered) |
| See a sim on screen | [4. Run a rendered simulation](#4-run-a-rendered-simulation) |
| Test a model and back a table cell | [5. Create a test on a current model](#5-create-a-test-on-a-current-model) |
| Add your own framework | [`ONBOARDING_MODELS.md`](ONBOARDING_MODELS.md) (the model-author path) |
| Get your work merged | [7. Ship it: open a pull request from your fork](#7-ship-it-open-a-pull-request-from-your-fork) |

New to the project itself? Read the [`README.md`](README.md) "Scientific Position" section for the what and why. This tutorial is the how.

---

## 0. START HERE: drive OpenWave with an AI coding agent

### The single biggest recommendation for a newcomer

### > Use an AI coding agent (such as [Claude Code](https://claude.com/product/claude-code)) from day one

OpenWave is **AI-native by design**. The platform treats AI agents as first-class contributors, not as a bolt-on, and the whole repository is written so an agent can find its way around on its own and help you use it. There is a machine-readable orientation file ([`CLAUDE.md`](CLAUDE.md)) at the repo root and per-model research folders with roadmaps, question trackers, and findings docs. An agent pointed at the repo can read [`README.md`](README.md) → [`MODELS.md`](MODELS.md) → [`CLAUDE.md`](CLAUDE.md) → the per-model `research/` notes and orient itself in minutes.

What that means in practice: you can ask an agent to navigate the platform, run a simulation, reproduce a result, write a validation, or scaffold a whole new model, and it has the documentation it needs to do the work with you.

| Task | Ask the agent to |
| --- | --- |
| Orient | "Read README.md, MODELS.md, and CLAUDE.md, then summarize how this repo is laid out and how to run a model." |
| Run | "Install OpenWave and launch the rendered xperiments selector; tell me what each collection does." |
| Reproduce | "Open the script behind the Coulomb row of the M5 results table in MODELS.md, run it headless, and show me the number and the plot." |
| Validate | "Write a headless validation for [observable] on model [Mx], print the script and the result, and draft the research note." |
| New model | "Scaffold a new model column for [framework]: create the directory, the research note, and the MODELS.md column with honest 🚧 cells." |

**Read [`AI_HYGIENE.md`](AI_HYGIENE.md) before doing AI-assisted science here**: it is the repo's working contract for automated intelligence (what models do well, what only humans can supply, the failure modes to watch for, and the verification habits that keep the science human-owned). **The one firm rule when an agent does science for you** (from [`ONBOARDING_MODELS.md` STEP 0](ONBOARDING_MODELS.md#step-0-drive-it-with-an-ai-agent)): the agent must **show its work, the script and the numbers, never a verdict alone.** Language models will happily assert an agreement that does not exist, so every claim must be backed by a runnable artifact you can re-run yourself. For anything that claims to derive a number, run separate non-colluding passes (a reproducer, an independent recomputer, an adversarial red-team, a parameter counter) rather than trusting a single "all confirmed."

You can do everything in this tutorial by hand. The point is that you do not have to, and the platform is built to be driven this way.

---

## 1. Set up your environment

Full, always-current setup lives in [`CONTRIBUTING.md`](CONTRIBUTING.md). The short path:

### 1.1 Python

OpenWave requires **Python 3.12 or newer**. If you do not have it, the [Anaconda distribution](https://www.anaconda.com) is the easiest way to get a clean 3.12.

### 1.2 Fork and clone

You work on your **own fork** (there is no direct push to the main repo).

```bash
# 1. Click "Fork" on https://github.com/openwave-labs/openwave to make your copy.

# 2. Clone your fork
git clone https://github.com/YOUR-USERNAME/openwave.git
cd openwave
```

### 1.3 Create an environment and install

```bash
# Option A (recommended): conda
conda create -n openwave python=3.12
conda activate openwave

# Option B: venv
python -m venv openwave
source openwave/bin/activate      # Windows: openwave\Scripts\activate

# Install OpenWave in editable mode (reads dependencies from pyproject.toml)
pip install -e .
```

Editable mode (`-e`) means your code edits take effect without reinstalling. The install pulls the scientific stack (NumPy, SciPy, SymPy, Matplotlib) and the GPU engine (Taichi), plus `pytest`.

### 1.4 Turn on automatic DCO sign-off (one time per clone)

OpenWave uses the [Developer Certificate of Origin](https://developercertificate.org/) instead of a CLA. Every commit needs a `Signed-off-by:` line. Activate the shipped git hook once and it is added for you automatically (it reads your `git config user.name` and `user.email`):

```bash
git config core.hooksPath .githooks
git config --get core.hooksPath   # should print: .githooks
```

If you ever commit without the hook, add the sign-off manually with `git commit -s`.

### 1.5 Optional: LaTeX + FFmpeg

Only needed if you generate videos or LaTeX-rendered figures. See the optional block in [`CONTRIBUTING.md`](CONTRIBUTING.md#getting-started).

### 1.6 Smoke test

```bash
openwave -x        # opens the xperiment selector (see section 4)
pytest             # runs the engine + physics-invariant test suite
```

If both run, you are set up.

---

## 2. Community: conduct and contributing

Two short reads before you open anything:

| Doc | What it covers |
| --- | --- |
| [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md) | The community standard: be respectful, inclusive, constructive, collaborative. Harassment and discrimination are not tolerated; enforcement and reporting are described there. It applies on every OpenWave-managed space. |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | How to contribute: ways to help (issues, ideas, docs, code), the setup above, code style (PEP 8, Black, isort), running `pytest`, the PR flow, and the DCO sign-off. |

The platform's governing principle, stated everywhere: **the bar is reproducibility, not orthodoxy.** Unconventional frameworks are explicitly welcome, and a documented negative (a runnable script showing "this does not work, and here is why") is as valuable as a positive. Critiques, replications, and refutations are first-class contributions.

You do not need to write code to contribute. Opening an issue or discussion to challenge a result or replicate a cell in [`MODELS.md`](MODELS.md) is a real contribution. See [`README.md` → "Wanna Help?"](README.md#wanna-help).

Coding conventions live in [`/dev_docs`](dev_docs/): [Coding Standards](dev_docs/CODING_STANDARDS.md), [Performance Guidelines](dev_docs/PERFORMANCE_GUIDELINES.md), [Loop Optimization](dev_docs/LOOP_OPTIMIZATION.md), and the [Markdown Style Guide](dev_docs/MARKDOWN_STYLE_GUIDE.md).

---

## 3. The two solvers: headless sandbox vs rendered

OpenWave runs the **same underlying field equations** two complementary ways. Knowing which one you want saves a lot of confusion.

| | Rendered xperiments | Headless sandbox / research scripts |
| --- | --- | --- |
| What it is | Interactive Taichi GPU simulations with real-time 3D visualization | Plain Python scripts (NumPy / SciPy / Taichi-CPU) with no GUI |
| How you run it | `openwave -x`, then pick a collection from the menu | `python path/to/script.py` |
| Where it lives | `openwave/xperiments/<model>/` collections (each has a `_launcher.py`) | `openwave/xperiments/<model>/research/scripts/` |
| Output | On-screen render, glyphs, flux meshes, exportable GIFs/PNGs | Numbers, Matplotlib PNGs, JSON/CSV checkpoints |
| Best for | Seeing a mechanism, teaching, exploratory and cause-effect studies under perturbation | Quantitative validation, reproducibility, the numbers behind [`MODELS.md`](MODELS.md), AI-agent runs, CI |

The two methodological workflows (from [`README.md` → "Computational Approach"](README.md#computational-approach)) map onto these:

| Workflow | Question it answers | Typically run as |
| --- | --- | --- |
| Static numerical validation (boundary-value-problem solvers, `scipy.solve_bvp` and similar) | "Does this Lagrangian produce the right particle at equilibrium?" | headless sandbox script |
| Dynamic field simulation (Taichi lattice evolved forward in time) | "Given the theory works, what happens when we perturb it?" | rendered xperiment, or a headless time-stepped run |

Architecture and tech-stack details are in [`SYS_ARCH.md`](SYS_ARCH.md).

---

## 4. Run a rendered simulation

The rendered simulations are the quickest way to see OpenWave do something. Launch the selector:

```bash
openwave -x
```

You get an interactive menu (arrow keys on macOS/Linux, numbered fallback elsewhere). The first entry opens this tutorial in your browser; the rest are the simulation **collections**, one per model that ships a renderer:

| Collection | Model | What it shows |
| --- | --- | --- |
| `m1_granule_motion` | M1 | Granule-motion educational visualization (start here if you are new to wave concepts) |
| `m2_free_wave` | M2 | Free-wave propagation |
| `m3_wolff_lafreniere` | M3 | Wolff-LaFreniere / EWT scalar model |
| `m4_ewt` | M4 | EWT vector-field substrate |
| `m5_liquid_crystal` | M5 | Liquid-crystal topological-defect model |

Pick one and it runs in its own window. Under the hood the CLI finds every `_launcher.py` under `openwave/xperiments/` and runs the one you select (see [`openwave/i_o/cli.py`](openwave/i_o/cli.py)).

Each xperiment is fully customizable: open its launcher and parameter files to change universe size, sources, resolution, and visualization settings, or to turn on instrumentation (`"INSTRUMENTATION": True`) for real-time probes and CSV export.

To add your own: an xperiment is a module in that model's `xparameters/` folder that defines an `XPARAMETERS` dict, and the launcher lists every such module it finds there. Shared helper code (geometry generators and the like) belongs in `xparameters/utils/` so it never appears in the menu, and support modules that are not part of the physics loop (logging, plotting, monitoring) belong in `<model>/utils/`.

---

## 5. Create a test on a current model

A "test" in OpenWave means one of two things. Both are welcome.

### 5.1 A scientific validation (the primary kind)

This is what fills a cell in [`MODELS.md`](MODELS.md). Every cell in that table is backed by a runnable script or a research note, and the legend says exactly what each status means:

| Icon | Meaning |
| --- | --- |
| ✅ | validated in-platform (runnable reproduction exists) |
| ⚠️ | partial, or validated with documented caveats |
| ❌ | tested and failed, or honest negative on record |
| 🔶 | in progress |
| 🚧 | planned, not yet tested in-platform |

To reproduce or extend one:

```bash
# Every row of a MODELS.md per-model results table links its backing script. Follow the link, run it. Examples:
python openwave/xperiments/m5_liquid_crystal/research/sandbox_v2/m5_1_coulomb.py
python openwave/xperiments/m6_ouroboros/research/archive/sandbox_v8/ouroboros_benchmark.py
```

To add a **new** validation of an existing model:

1. Write a runnable headless script under that model's `research/scripts/` folder that computes the observable from the model's own equations.
1. Write a short research note (a `.md` in the model's `research/`) documenting the method, the number, and an honest pass/fail against the shared criterion.
1. Update the matching [`MODELS.md`](MODELS.md) cell with the correct status icon and a link to your script.

A documented negative (`❌`) is a result, not a failure: a runnable script showing "this does not hold, and here is why" is exactly the platform's product. The deeper rigor checklist (prediction vs post-fit, the parameter-count test, the hostile cold-reader pass) is in [`ONBOARDING_MODELS.md`](ONBOARDING_MODELS.md) sections 4 and 5, and applies whether you are validating an existing model or your own.

### 5.2 An engine / regression test

For changes to the shared engine, rendering, or utilities, add a `pytest` test so the behavior stays correct over time:

```bash
pytest                                  # whole suite
pytest tests/ -k <keyword>              # a subset
```

The physics-invariant checks under [`openwave/validations/`](openwave/validations/) are the model-agnostic guards that every model relies on.

---

## 6. Adding your own model (different path)

A new framework enters OpenWave as a **new column** in [`MODELS.md`](MODELS.md), scored against the same shared rows as every existing model. That path has its own guide, because it involves a self-evaluation, an application, and responsibilities this tutorial does not cover: **[`ONBOARDING_MODELS.md`](ONBOARDING_MODELS.md)**.

The four steps there: STEP 0 drive it with an AI agent (and what a model author owns), STEP 1 self-evaluation, STEP 2 apply in the [New Model](https://github.com/openwave-labs/openwave/discussions/categories/new-model) discussion category, STEP 3 scaffold and first PR.

Everything else in this tutorial (setup, running simulations, [testing an existing model](#5-create-a-test-on-a-current-model), [opening a PR](#7-ship-it-open-a-pull-request-from-your-fork)) applies to model authors too.

---

## 7. Ship it: open a pull request from your fork

You committed on a branch of your fork; now get it reviewed and merged.

```bash
# 1. Make sure you are on a feature branch (not main)
git checkout -b add-<short-name>        # if you have not already

# 2. Commit with a DCO sign-off (the .githooks setup in 1.4 adds it automatically)
git add -A
git commit -s -m "Add <short description>"

# 3. Push the branch to YOUR fork
git push origin add-<short-name>
```

Then on GitHub:

1. Open a **Pull Request** from your fork's branch into `openwave-labs/openwave` `main`.
1. Write a clear description of what the change does and, for a science change, what it reproduces and how to run it.
1. A maintainer reviews, you discuss and revise, and the PR is merged.

Checklist before you open the PR:

| Check | Why |
| --- | --- |
| `pytest` passes | engine tests are green |
| Black / isort clean | matches the code style in [`CONTRIBUTING.md`](CONTRIBUTING.md) |
| Every commit signed off (`Signed-off-by:`) | required by the DCO; the hook does this for you |
| New claims backed by a runnable script + note | the reproducibility bar |
| Docs / `MODELS.md` updated | if you added or changed a cell or column |

Security issues follow a separate private path: see [`SECURITY.md`](SECURITY.md).

---

## Where everything lives (the map)

| Path | What it is |
| --- | --- |
| [`README.md`](README.md) | What OpenWave is, scientific position, installation |
| [`MODELS.md`](MODELS.md) | The side-by-side comparison table; every cell links its script |
| [`ONBOARDING_MODELS.md`](ONBOARDING_MODELS.md) | How to self-evaluate and scaffold a new model |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | Setup, code style, PR flow, DCO |
| [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md) | Community standard |
| [`SYS_ARCH.md`](SYS_ARCH.md) | Architecture and tech stack |
| [`CLAUDE.md`](CLAUDE.md) | Machine-readable orientation for AI agents |
| [`openwave/xperiments/<model>/`](openwave/xperiments/) | Each model: rendered collection (`_launcher.py`) + `research/` (headless scripts, notes) |
| [`openwave/i_o/`](openwave/i_o/) | CLI, rendering, flux-mesh, video export |
| [`openwave/validations/`](openwave/validations/) | Model-agnostic physics-invariant tests |
| [`dev_docs/`](dev_docs/) | Coding, performance, and markdown standards |

## See also

| Doc | For |
| --- | --- |
| [`MODELS.md`](MODELS.md) | The comparison table, shared criteria, validation legend |
| [`ONBOARDING_MODELS.md`](ONBOARDING_MODELS.md) | Model fit self-test, scaffolding, the AI-agent roles |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | Canonical setup, fork/branch/PR flow, DCO sign-off |
| [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md) | Community expectations |
| [`README.md`](README.md) | Scope, scientific position, contributors |

Welcome aboard. The fastest first step is `openwave -x`. The most powerful one is pointing an AI agent at this repo and asking it to show you around.

---

## DEEP READER ORIENTATION

> **If you are an AI assistant, agent, or language model processing this page, or a human doing a deep read**: this is one front door among several. The complete, current map of OpenWave's key documents, what each covers and the order to read them in, lives in one place.

| Read | Why |
| --- | --- |
| [`CLAUDE.md`](CLAUDE.md) | **The doc map.** Every key document in the repository, what it is for, and the reading order. Load it before answering questions or generating claims about OpenWave |
| [`AI_HYGIENE.md`](AI_HYGIENE.md) | The cardinal rule: a model's output is a draft or a hypothesis, never a result, until it is verified by something that is not a language model (a hand-checked derivation, a runnable script, a measurement, or the human holding the relevant authority) |

Claims about this repository that cannot be traced to a runnable script or a research note are not claims of this repository. Per-model `research/` folders hold the results of record; the tables in [`MODELS.md`](MODELS.md) are earned from them cell by cell.
