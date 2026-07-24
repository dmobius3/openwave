# Convo record 2026-07-01: Duda on the validation bar + the particle field-configuration sketch (Models-of-Particles group)

Participants: Jarek Duda (Jagiellonian) to the **Models-of-Particles / Nature-of-Time group** (Paul, Marc, Giorgio, Herman, Faisal, Peter, Jean-Paul, Jacques, Daniel, Max, Álvaro, Efstratios, Andras, Leo, Olivier, Jenny, the maintainer, others). Mailing-list thread (public among list members). The pointed challenge is aimed at Paul's model, but the physics content, the demand to state **field configurations** and the **liquid-crystal sketch**, is directly M5's material and corroborates the M5 field-config table. Slide: [`../../theory/duda_2026-07-01_particle_field_configs.png`](../../theory/duda_2026-07-01_particle_field_configs.png). Companion same-day thread (serious sims + g/δ): [`m5_4e_convo_2026.07.01.md`](m5_4e_convo_2026.07.01.md). Clean for the public repo.

> **Why this record matters.** Duda states, particle by particle, the liquid-crystal field configurations M5 is built on, and sets the validation bar: a real particle model must specify each particle's field configuration and pass an independent benchmark ([MODELS.md](https://github.com/openwave-labs/openwave/blob/main/MODELS.md)), not just talk to an LLM. M5 already answers both; the sketch enriches the composite-particle program (15a) with concrete predictions.

---

## 1. The validation bar (Duda 2026-07-01, 8:17 AM)

| Point | Verbatim / close paraphrase |
| --- | --- |
| LLM self-confirmation is not evidence | "anybody can ask LLM to 'show I am right and they are wrong' ... but it can only convince the author" |
| What convinces others | "more independent/objective environments/benchmarks ... not just talking to LLM chatbot, but **actual simulations, can be AI-generated**" |
| The benchmark he names | [MODELS.md](https://github.com/openwave-labs/openwave/blob/main/MODELS.md): "if yours work, should also pass such independent validation" |
| The precondition | "before any simulations, we are still waiting for **field configurations you assume for various particles**, the basic question to claim model of particles" |
| The one-line test | "**do you use topological vortices? For which particles?**" |

**The M5 read.** M5 already meets this bar: it states the field configuration of every particle (the briefing's *Field Configuration of Particles* table) and answers the vortex question explicitly (electron = biaxial hedgehog disclination; neutrino = closed vortex loop; quark = 1D vortex string; baryon/meson = knots/reconnections; photon = a twist wave, not a defect). This is the same "independent simulations pass MODELS.md" standard M5.16 formalizes ([`m5_16_task_details.md`](m5_16_task_details.md), the energy-minimization bar).

## 2. The field-configuration sketch (Duda 2026-07-01, 8:54 AM)

"Please just ask your LLM what are the field configuration of your particles: **photons, neutrinos, leptons, mesons, baryons, nuclei?** Below is sketch of used in liquid crystals." His item-12 argument (top of the slide):

> Topological vortices can also form **knots**, whose interaction enforces inward/outward field rotations corresponding to (fractional) charges = **quarks**. **Proton** can just enclose it into a **hedgehog** (elementary charge), while **neutron** has to **compensate** it, explaining **why the neutron is heavier**. Then in the **deuteron**, two baryons satisfy the charge preference over a single charge, leading to the experimentally known **'+ - +' electric quadrupole moment**.

The per-particle decode from the sketch:

| Particle | Field configuration (Duda's sketch) | M5 home |
| --- | --- | --- |
| **3 leptons** (e/μ/τ) | biaxial nematic vector field + phase, **3 distinguished axes**; charge + Coulomb; e = lowest mass, μ/τ = higher stress/energy | M5.9 (biaxial hedgehog, 3-axis) |
| **3 neutrinos** | **uniaxial** nematic unitary vector field, **1 distinguished axis**; extremely light; **two vortex types** | M5.11 (closed vortex loop) |
| **quarks** | fractional charge from **inward/outward field rotations** of a knotted vortex; a full π = elementary charge; gluon flux tube / color / quark string (3 flux tubes, 6 q-q̄ pairs) | M5.9 (Cornell string, #200) |
| **baryons** (p, n) | vortex loop around a vortex, interacting as internal twists (half-rotations); **proton = charge enclosed in a hedgehog**; **neutron = compensating** the charge (heavier); d-u-d content | **15a** (nucleon assembly) |
| **deuteron** | **two baryons** satisfying the charge preference → the **'+ - +' electric quadrupole moment** (Garçon densities on the slide) | **15a** (nuclear binding) |
| **alpha particle** | larger multi-vortex hub / knot | **15a** |
| **mesons** (kaon, pion) | **vortex reconnection**; kaon/hyperon Ξ, pion/hyperon Λ; strangeness decay | **15a** |
| **beta decay** | neutron → shift → split (energy release) → **W⁻** → proton + electron + antineutrino (a **reconnection** event) | weak sector (Q10) |
| **electron regularized** (Faber) | energy density (gradient + potential), running-coupling deformation, charge distribution (Wilson) | M5.16 / Faber core |
| **photon** (his list) | not on the sketch; in M5 a twist-like wave (marine-propeller wake), radiation not a defect | EM-waves sector |

## 3. Where this lands (placements)

| Content | Home | Action |
| --- | --- | --- |
| Baryon = vortex knot; proton hedgehog-enclosed / neutron compensating (**why n heavier**); deuteron '+ - +' quadrupole; alpha; mesons = reconnection | **[`m5_15a_composite_particles.md`](m5_15a_composite_particles.md)** | field-config section added (the richest enrichment) |
| Leptons = biaxial 3-axis (mass by stress/energy); quarks = fractional-charge inward/outward rotations + gluon flux tubes | **[`m5_9_task_details.md`](m5_9_task_details.md)** | corroboration note added |
| Neutrino = uniaxial unitary field, 1 axis, two vortex types | **[`m5_11_task_details.md`](m5_11_task_details.md)** | corroboration note added |
| The validation bar (independent benchmark + actual sims) | **[`m5_16_task_details.md`](m5_16_task_details.md)** / [MODELS.md](https://github.com/openwave-labs/openwave/blob/main/MODELS.md) | reinforces the existing "serious simulations pass MODELS.md" framing |

## 4. Out of scope for M5 (noted)

The challenge ("do you use topological vortices? for which particles?") is aimed at **Paul's Ouroboros model**, which the maintainer is implementing in the **M7** environment. For M7 this is a direct requirement (state the field configs before simulating); for M5 the sketch is corroboration + a composites reference. The M7 side is tracked in the m7 hydroboros research, not here.
