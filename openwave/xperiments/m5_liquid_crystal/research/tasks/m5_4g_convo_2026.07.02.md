# Contribution record 2026-07-02: Rudolf Golubich on the SU(2) soliton lattice method (Faber model)

Correspondent: **Dr. Rudolf Golubich** (Neusiedl am See; [ORCID 0000-0003-1990-3894](https://orcid.org/0000-0003-1990-3894), [statistikmacher.de](https://www.statistikmacher.de/)), **co-author of Faber & Golubich [arXiv:2604.12021](https://arxiv.org/abs/2604.12021)** (the high-precision lattice SU(2) solitonic-dipole paper). Direct email to the maintainer (cc Manfried Faber), in reply to the maintainer's note about integrating the topological SU(2) soliton model into OpenWave. This record captures the **technical content only** (the model, the method, the calibration); personal/logistical details stay in the email. The Faber SU(2)/S³ soliton model is what M5.6 ports (Faber regularization) and what M5.11 P1 reproduced; Golubich is an authoritative source for it. Clean for the public repo.

> **Why this matters.** Golubich is a co-author of the model M5 consumes, and he hands us the exact **lattice energy-minimization recipe** (nonlinear CG + line search, the derivative scheme, the single free parameter `r₀`) that [`m5_16_task_details.md`](m5_16_task_details.md) is being built to run. It validates the M5.16 static-minimization fork and de-risks its numerics. Added to the OpenWave README contributors table (after Faber).

---

## 1. The model + the single-parameter result

The topological SU(2) / S³ soliton (Faber's Model of Topological Fermions, MTF): charged fermions are stable topological solitons, charge = winding number, the soliton held at fixed size (evades Derrick).

| Claim (Golubich) | Detail |
| --- | --- |
| **One free parameter** | `r₀` only |
| **The calibration** | `r₀ = (π/4) × the classical electron radius` reproduces **both the electron mass AND the fine-structure constant** |
| **QED cross-check** | also matches the **perturbative QED vacuum-polarization corrections** "with remarkable accuracy" |

This is exactly the anchor M5 already reproduced in the M5.11 P1 minimizer (Faber electron 511.00 keV at `r₀=2.2132 fm`, `α⁻¹→137.03`); Golubich confirms it is the authoritative, single-parameter result of the source model.

## 2. The method (the recipe for M5.16)

They have **not** simulated the full dynamical model, only **static energy minimization** ("which already poses substantial numerical challenges"):

| Element | Golubich's prescription |
| --- | --- |
| Objective | minimize the **total energy** on the field degrees of freedom; the minimum **directly gives the mass** |
| Minimizer | **nonlinear conjugate gradient**, **Polak-Ribière** update |
| Line search | **bracketing + golden-section** search along the gradient |
| Cost | large lattices take **several weeks** on their hardware (matches Duda's "Faber said weeks", `4e §1`) |
| The main challenge | **integrating the derivative calculations** into OpenWave, the difference quotients are in Faber's `derivatives.tex` note |

## 3. Charge and spin as topological invariants (analytic, not fitted)

The energy minimization gives the mass; **charge and spin fall out of the analytic solitonic solution**, not the numerics:

| Observable | Origin (Golubich) |
| --- | --- |
| **Charge** | the **sign of `nᵢ`**: how the field rotates along the radial axis moving out of the soliton core |
| **Spin** | the behaviour of **`α(r)`**: covering the **upper half of SU(2)-S² → +1/2**, the **lower half → −1/2** |
| **Link to electrodynamics** | the SU(2) field couples to ED via the **rotation axis `n̂`, which is constant along electric field lines** |

## 4. Where this lands (placements)

| Content | Home | Action |
| --- | --- | --- |
| The full CG method (Polak-Ribière + bracketing/golden-section), the `r₀` single-parameter calibration, the derivative-integration challenge | **[`m5_16_task_details.md`](m5_16_task_details.md)** | method note added: Golubich's recipe is the concrete build for the M5.16 minimizer |
| Charge = sign of `nᵢ`; spin = `α(r)` upper/lower SU(2)-S² half; SU(2)↔ED via `n̂` along E-lines | electron-ID / M5.11 P1 (Faber electron) | referenced; the clean topological readout to adopt |
| `r₀ = (π/4)·r_classical` → mass + α + QED vacuum polarization | Faber calibration (M5.6.3 / Q8 in [`../m5_question_tracker.md`](../m5_question_tracker.md)) | corroboration; authoritative source |

## 5. Shared files (local-only) + presentation

Three `.tex` files were shared, now kept **local-only** in [`../../theory/golubich_faber_su2/`](../../theory) and **gitignored** (private author sources, not redistributable on a public repo; manifest in [`../../theory/_CITATIONS.md`](../../theory/_CITATIONS.md) #51-53):

| File | What |
| --- | --- |
| `MTF.tex` | the Model of Topological Fermions source (the SU(2) model) |
| `dipole.tex` | Faber & Golubich [arXiv:2604.12021](https://arxiv.org/abs/2604.12021) source, the compact model + lattice-implementation intro |
| `derivatives.tex` | Faber's German note on the **difference quotients** for the derivative calculation (guides the M5.16 gradient) |

Golubich is preparing a presentation for the **ÖPG-CMD joint meeting 2026** (how electron mass, spin, and `α` emerge from the topology), shared in advance at [statistikmacher.de/DEV/Solitons](https://www.statistikmacher.de/DEV/Solitons/) (access credentials sent privately in the email, not recorded here).

## 6. Standing offer (the collaboration model)

Golubich considered integrating the soliton model into OpenWave himself but is time-constrained and unfamiliar with the codebase. The working arrangement (the maintainer's reply): rather than have him set up and learn the OpenWave environment now, he **sends insights by email**, and the maintainer ensures they are **processed into the repo** (this record is the first). This is the model for the multi-particle-model validation initiative OpenWave hosts.
