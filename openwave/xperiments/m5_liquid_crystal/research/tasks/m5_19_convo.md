# M5.19 convo record (the program thread "Regularized neutrino vortex-loop: runs + questions")

Per-task convo record (public group thread; models-of-particles CC'd, so tracked in-repo). Outbound context: the M5.19 batched ask sent 2026-07-10 evening with the method-note link ([`../findings/m5_19_method_note.md`](../findings/m5_19_method_note.md)); tracker registry [`../m5_question_tracker.md`](../m5_question_tracker.md).

## 2026-07-10: Duda's reply to the M5.19 batch: the radiation-mechanism directive (spec for M5.20)

**From**: Jarek Duda, 2026-07-10 ~20:03 EDT, reply-all on the program thread. **New CCs**: Marcin Misiaczek (neutrino-oscillation specialist), Jacek Jendrej (radiation from solitons), models-of-particles. Opening remark: "analysis now with Fable5 is much better than previous with Opus 4.8" (his words, addressed at convincing the group).

### Verbatim core (his answers, lightly trimmed)

> Sure, such neutrino as loop of topological vortex has mass/energy density per length, so static energy minimization could reduce radius R of this loop to zero. However, this is a dynamical process, total energy rather has to be conserved - for neutrino to lose energy this way, this energy would need to be transferred into some different field excitations. In other words, it would require radiation mechanism - which personally I cannot imagine without interaction with particles (e.g. beta decay/capture) - just flying through empty vacuum (as potential minimum). To investigate it, instead of just energy minimization, there is required real evolution e.g. with Euler-Lagrange equations for the assumed Lagrangian, searching if there is a way to radiate this energy as some different field excitations (?) If such radiation mechanism is not found (?), we can assume that total energy during flight of neutrino at least in vacuum is conserved.
>
> Regarding the question if in vortex cross-section there is 2D topological charge 1 or 1/2 [...] While 1/2-vortex should have lower energy, charge as hedgehog seems to topologically require 1-vortices in cross-section. Therefore, I suspect electric-charge based sources of neutrino (beta decay, leptons) should use 1-vortex, what also concerns their absorption e.g. by neutrino detectors. However, e.g. in extremely high temperature like just after Big Bang, from vacuum fluctuations there could be also formed 1/2-vortices and their loops as kind of half-neutrinos, which would be nearly impossible for us to observe, but contributing e.g. to mass as one of types for dark matter/energy.
>
> Q16 (protection): What protects the physical loop? - I suspect lack of radiation mechanism at least flying through empty vacuum.
>
> Q21 (the temporal boosts): the Lagrangian is still the same, Hamiltonian has negative energy contributions Gamma * tilde{Gamma} - by energy minimization preferring nonzero time derivatives (oscillations), which require also local boosts corresponding to gravity. Marcin Misiaczek suggested there might be gravitational effect for oscillations e.g. in MSW.
>
> Q20 (calibration): required spectral deformation in vortex center means positive energy contributions from potential. There would be also negative contribution from oscillations-gravity, but I suspect it should be smaller - total energy density per length should be positive. But if there is no radiation mechanism, energy density per length could be also zero or negative - just creations/absorption of such vortex loop is quite non-trivial, requires e.g. beta.
>
> Q22 (delta): I don't know in this moment.
>
> Seems the main question now is radiation mechanism - if topological vortex loop in empty vacuum could transfer this energy into some field excitation - using not just energy minimization, but real evolution e.g. with Euler-Lagrange for this Lagrangian. There is also usually huge kinetic energy from e.g. beta decay, but special relativity allows to boost to rest frame leaving oscillations as the only dynamics. However, this boost is extremely large, and introduces very strong time dilation - so maybe there is radiation mechanism (neutrino losing rest mass with time), but time dilation makes it nearly negligible - e.g. for minutes from our Sun, for millions of years from supernova explosion (?)

Three figures attached: the beta-decay quark-string split/reconnect cartoon (neutron → proton + electron + neutrino as a released vortex loop); the hedgehog electron cross-section (charge topological, requiring 1-vortices); the Hamiltonian matrix slide (the `−(Γ0¹Γ̃μⁱ − Γμ¹Γ̃0ⁱ)²` clock/gravity terms, "clock propulsion?" annotation).

### Decode + routing (per the unknowns discipline; evidence, not resolution, until measured)

| Item | Decode | Routing |
| --- | --- | --- |
| The M5.19 statics negative reframed | Gradient flow = the infinite-friction limit: it hands the field a free dissipation channel, so unwinding under it never tested his protection mechanism. Conservative dynamics has no such channel unless the loop genuinely radiates | Machine-checkable → [M5.20](m5_20_task_details.md) (the run) |
| Q16 protection | Author answer: protection = absence of a radiation mechanism in vacuum (dynamical, not topological). A hypothesis he himself flags with "(?)": test, don't assume | Tracker Q16 → RESOLVED (author answer on record); empirical verdict = M5.20 |
| Degree 1 vs 1/2 | Author answer: beta-sourced (physical) neutrinos = 1-vortex, topologically required by hedgehog charge; 1/2-loops = primordial half-neutrinos, dark-sector candidates. M5.19's 2:1 energetic preference for 1/2 is reinterpreted, not contradicted | Seed priority in M5.20: q=1 first, q=1/2 kept as the dark-sector class |
| Q21 temporal boosts | Equation-level: same Lagrangian, Hamiltonian negative `Γ·Γ̃` contributions prefer nonzero time derivatives (oscillations = clock), requiring local boosts = gravitational mass. Marcin: possible gravitational effect in MSW | Tracker Q21 → RESOLVED (author answer); residual = the Γ↔M dictionary for our stack (new tracker item Q23, author-gated, next batch); clock sector stays M5.20.2 |
| Q20 calibration | PARTIAL: he answered the energy-sign framing (core potential positive, oscillation-gravity negative smaller, per-length total positive), NOT the Dirichlet-term presence question as asked | Tracker Q20 stays OPEN (partial 2026-07-10); M5.20 runs both functional variants so no verdict hinges on it |
| Q22 sector δ | "I don't know in this moment" | Parked: author-gated, unanswered by choice; stays OPEN at low priority |
| Time-dilation remark | Radiation might exist but be near-negligible in lab frame (minutes from the Sun, Myr from supernovae): ⚠️ speculative, needs the clock/boost machinery | Deferred M5.20.2 |
| Social layer | Two specialists CC'd (Misiaczek: oscillations; Jendrej: soliton radiation) + models-of-particles; public quality endorsement of the current analysis workflow | Raises the method-note rigor bar for M5.20 (Jendrej reads radiation claims professionally) |
