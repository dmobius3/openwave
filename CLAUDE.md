# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

OpenWave is an open-source subatomic wave simulator for exploring fundamental physics through **classical field theory enriched with topology and nonlinearity** — the scientific tradition of de Broglie–Bohm pilot waves, wave structure of matter, and modern topological-soliton models. The platform tests whether particles and forces can emerge from deterministic field equations rather than being postulated.

The simulator runs multiple candidate Lagrangian frameworks (scalar-field, vector-field) in a shared numerical engine, plus a granule-motion model for educational visualization. GPU acceleration uses Taichi Lang.

### Project Goals

OpenWave investigates, in one integrated simulator, four primary domains: **matter** (particle emergence from topological defects + wave dynamics), **forces** (electric, strong, magnetic, gravitational from one classical-field framework), **electromagnetic waves**, and **heat** (whether wave-field degrees of freedom contribute to thermal physics). Each domain has concrete pass/fail criteria applied uniformly across candidate models.

### What is OpenWave?

| Reference | Purpose |
| --- | --- |
| `README.md` | Full description, scope, scientific position, installation |
| `SYS_ARCH.md` | Module structure and system architecture |

### Theoretical Advisors and Candidate Frameworks

| Contributor | Framework | OpenWave Model |
| --- | --- | --- |
| Jeff Yee | Energy Wave Theory (EWT) | M3 |
| Dr. Jarek Duda | Liquid Crystal Particle Analogs (arxiv 2108.07896, 2501.04036) | M5 |
| Dr. Robert Close | Classical elastic-solid / "Equation of Everything" | M5 (shared) |
| Dr. Manfried Faber | LdG regularization (Universe 11/2025/113) | M5.6 baseline |
| Dr. Paul Werbos | Ouroboros chaoiton Lagrangian | M6 |

### Known Challenges & Limitations

- Full Planck-scale fidelity is computationally prohibitive; resolution is user-tunable per xperiment.
- Uses a dedicated physics computational backend (Taichi GPU), independent of 3D modeling software.

## Project Architecture

| Path | Contents |
| --- | --- |
| `openwave/xperiments/m1_granule_motion/` | Educational granule-motion model |
| `openwave/xperiments/m2_free_wave/` | Free-wave propagation |
| `openwave/xperiments/m3_wolff_lafreniere/` | Wolff-LaFreniere / EWT scalar model |
| `openwave/xperiments/m4_ewt/` | M4 EWT model (vector-field substrate, in development) |
| `openwave/xperiments/m5_liquid_crystal/` | **Active** — Duda LCB topological-defect model |
| `openwave/xperiments/m6_ouroboros/` | Werbos chaoiton Lagrangian |
| `openwave/common/`, `i_o/`, `validations/`, `video_export/` | Shared utilities, rendering, physics-invariant tests |

Refer to `README.md` and `SYS_ARCH.md` for the full Modules Structure and Objects Map.

### Scientific Source Material

Each model directory under `openwave/xperiments/` contains a `/research` subfolder with active research notes, plus a `theory/` folder (at the model root) holding the foundational theorist papers for that model.

Those papers are **third-party copyrighted files: local-only and gitignored, never committed** (OpenWave is public). The tracked record of each `theory/` folder is its `_CITATIONS.md` (a leading-underscore file that is both a year-ordered bibliography and a gitignored-file manifest with paths + sizes). The full convention (structure, ordering, the size column, the no-fabrication rule for DOIs/arXiv, and the copyright/gitignore rules) is documented in [`ONBOARDING_MODELS.md` § 3.2](ONBOARDING_MODELS.md), the guide for scaffolding new models. **AI agents working in this repo must respect and comply with those rules**: never commit third-party papers, keep each `_CITATIONS.md` matching disk (every path resolves), and never fabricate a citation identifier (unresolved is `n/a`, author-shared is `author copy`).

Note: the legacy top-level `scientific_source/` folder was retired 2026-05-18; papers now live per-model.

**Research `data/` folders follow the same local-only pattern (2026-07-20).** Heavy binary arrays (`.npz` `.npy` `.h5` `.pkl` `.pt` `.mat`) are gitignored and **kept on the working machine, never deleted** (this supersedes an earlier rule that deleted arrays over 1 MB at task close). What IS tracked, and what a reader on GitHub audits: the distilled summary `.json` / `.csv` / `.txt` beside them, the plots, the scripts that produce everything, and the per-folder **`_DATASETS.md`** manifest, which indexes every local-only file by task group with its producing script and the task record holding the exact regen command. Regenerate a manifest after any run that writes arrays: `python3 ~/Documents/source_code/XRODZ/dotfiles/claude_projects/NEPTUNYA-SABER/scripts/gen_datasets_manifest.py <path/to/data> --write`. Runs are deterministic from their fixed seeds and configs, so a clone rebuilds every array from tracked code.

## Installation & Usage

Refer to `README.md` for installation, CLI usage (`openwave -x`), and the instrumentation framework.

## Physics Context

OpenWave implements classical-field-theory-with-topology-and-nonlinearity approaches:

- Topological defects provide static structure (integer charge, spin).
- Klein-Gordon-like wave dynamics around the vacuum field provide mass and relativistic kinematics.
- Standing-wave interference between defect emissions produces orbit quantization.
- Particles are **time-periodic resonances** (Zitterbewegung clocks), NOT static solitons — Derrick's theorem forbids static stable solitons, confirmed empirically in M5.2.

- Also refer to ../CLAUDE.md file to search for any available context to the OpenWave project in a parent directory.

## Code Style & Documentation Standards

| Doc | Purpose |
| --- | --- |
| [Markdown Style Guide](dev_docs/MARKDOWN_STYLE_GUIDE.md) | All `.md` files |
| [Coding Standards](dev_docs/CODING_STANDARDS.md) | Python code |
| [Performance Guidelines](dev_docs/PERFORMANCE_GUIDELINES.md) | Optimization |
| [Loop Optimization](dev_docs/LOOP_OPTIMIZATION.md) | Critical loops |
| [Coordinate System](dev_docs/COORDINATE_SYSTEM.md) | Spatial conventions |
| [Floating Point Precision](dev_docs/FLOATING_POINT_PRECISION.md) | Numerical precision rules |
| [Scaling Factor](dev_docs/SCALING_FACTOR.md) | Physics unit scaling |
| [Version Management](dev_docs/VERSION_MANAGEMENT.md) | Release versioning |
| [Wave Diagnostics](dev_docs/WAVE_DIAGNOSTICS.md) | Validation diagnostics |
| [Method Note](dev_docs/METHOD_NOTE.md) | **MANDATORY** standard for any report/email to a model's theory owner or external physicist |
| [AI Hygiene](AI_HYGIENE.md) | **MANDATORY** working contract for AI-assisted research: division of labor, failure modes, verification habits |

### AI hygiene (all AI-assisted work): MANDATORY

Every AI agent working in this repo (including you) operates under [`AI_HYGIENE.md`](AI_HYGIENE.md): model output is a draft or hypothesis, never a result, until verified by something that is not a language model (a hand-checked derivation, a runnable script, a lattice measurement, or the confirmation of the human holding the relevant authority). Author-gated questions (intent, provenance, definitions) can only be answered by the author; externally received AI-derived material is tagged evidence-not-resolution until confirmed; anything community-facing is human-owned prose over script-backed results. Read it before doing research work here.

### The adversarial audit (every substantive claim): CARDINAL RULE

Before any substantive derivation, verification, or headline claim is trusted, recorded as a result, or sent outside this repo, an **independent second agent audits it adversarially**: instructed to REFUTE (not confirm), with its OWN implementation (different method / seed / construction), its own hand re-derivations, and a per-claim verdict (CONFIRMED / REFUTED / QUALIFIED) backed by its own numbers. Fold every catch back into the artifact and **record the audit outcome IN the deliverable** (the [`m5_18_verification_note.md § 10`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_18_verification_note.md) pattern). Full rule: [`AI_HYGIENE.md § 1`](AI_HYGIENE.md). Origin: Duda 2026-07-03 ("careful small steps, maybe multiple agents verifying each other"); first real catch: M5.18 (an overclaimed witness refuted + a missed vacuum-branch structure found by the auditor).

### Method note (model-owner-facing output): MANDATORY

Any report, summary doc, or email addressed to a model's theory owner (M5: Duda; M7: Fleury; any advisor) follows [`dev_docs/METHOD_NOTE.md`](dev_docs/METHOD_NOTE.md): **equations first** (Hamiltonian / potential / observable definitions in math notation), an **equation-to-code map** with absolute GitHub hyperlinks (`#L` anchors; `blob/main` for frozen task-scoped files, commit-pinned only for live/evolving files like the root engines and production rendering code), the functional in a **small auditable module**, results after methods, a minimal physics-first inspection set, and the **adversarial audit recorded in the note**. Adopted 2026-07-03 after the M5.16 audit failure ("still I have no idea what does it calculate"): results a physicist cannot audit by reading carry no weight, regardless of correctness.

### Important: Markdown Linting Requirements

When editing any `.md` files, ALWAYS ensure compliance with markdown linting rules:

- Add blank lines around headings (before and after)
- Add blank lines around lists (before and after)
- Add blank lines around code blocks (before and after)
- Use consistent ordered list numbering (1, 1, 1 style)
- Check for proper spacing and formatting

Run linting checks after editing to catch any issues before committing.
