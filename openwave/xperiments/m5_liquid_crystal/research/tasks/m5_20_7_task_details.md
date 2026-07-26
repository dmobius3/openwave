# M5.20.7: the neutrino-aging phenomenology check (Duda's 2026-07-15 request)

**Status**: 🚧 PARKED STUB (2026-07-16). His direct request at the Q24 reply ([`m5_20_convo.md § 2026-07-15`](m5_20_convo.md)): "Could you check if it fits various neutrino data?". **User decision 2026-07-16: runs AFTER the M5.21 series** (electron-first); this stub holds the scope so nothing is lost.

## The hypothesis to check (his words, decoded)

| Element | His statement (2026-07-15) | What it predicts |
| --- | --- | --- |
| Rest-frame radiation | "If neutrino as vortex loop is radiating (in rest frame)" (the M5.20.3 blowup read as physics rather than pathology, slowed by time dilation) | a finite rest-frame lifetime; lab-frame decay length ∝ γ = E/m |
| "Neutrino aging" | "experimentally it would be seen as ... (growing with distance) energy threshold below which we don't see neutrinos" | low-energy neutrinos deplete first over long baselines |
| The final flash | "there might be final flash for R->0 (ending such aging) we could observe" | a photon/particle burst at end-of-life |
| DAMA/LIBRA | "DAMA/LIBRA data with annual oscillations seems fitting such 'final neutrino flash' hypothesis" | the annual modulation re-read as flux-modulated final flashes rather than WIMP scattering |

## Planned shape (literature + order-of-magnitude; NO simulation, no GPU)

| Step | Content |
| --- | --- |
| N1 | Formalize: rest-frame lifetime τ₀ → lab decay length L(E) = γcτ₀; the survival-probability threshold curve vs baseline |
| N2 | The constraint gauntlet: SN1987A (168 kly baseline, ~10-30 MeV, arrival intact = the hardest bound), solar (8 lm, sub-MeV to ~15 MeV), atmospheric + accelerator (km-scale), reactor, IceCube astrophysical (Gly, TeV-PeV); extract the allowed (τ₀/m) window per dataset, if any |
| N3 | DAMA/LIBRA scoping: event rate, modulation phase/amplitude vs what a final-flash flux would need; consistency with the null results of the other direct-detection experiments under a flash (non-scattering) signature |
| N4 | Verdict table: fits / conflicts / needs-a-number-from-the-model (the loop's rest-frame radiated power is exactly what the parked loop-side sims would have to supply) |

Records at run time: `../findings/m5_20_7_note.md` + this file; sources logged to [`../../theory/_CITATIONS.md`](../../theory/_CITATIONS.md). Adversarial audit per the cardinal rule (claims here are data-vs-model statements, auditable from the cited numbers).

## Gating

Roadmap `Gated By`: **the M5.21 series (user decision 2026-07-16: electron-first)** + user "go". Related parked item, recorded in the convo only (no task): his Boötes-void active-formation question.
