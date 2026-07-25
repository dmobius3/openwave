# Checkpoint 04 - #200 full production engine build COMPLETE (local)

Go 2026-06-20 08:31 EDT. All sub-tasks 1-6 done on the production 4x4 engine. Nothing posted to
GitHub / committed (the maintainer's call). The distilled writeup goes to issue #200 only after he commits.

## Deliverables (all LOCAL in sandbox_v9)

| Artifact | What |
| --- | --- |
| `m5_9_4_engine_lepton.py` | production 4x4 engine + signed-energy ledger (boost-GEM negative gravity + clock-negative oscillation) + gradient-flow minimization + electron calibration; env-parameterized |
| `m5_9_5_higgs_clock_spectrum.py` | Higgs/confiner scale-selection (A) + #220 dynamical clock omega(r0) (B) + delta-scan spectrum (C) |
| `m5_9_lepton_mass_clock_findings.md` | § "M5.9.4-5: the full production engine 4x4 build (FIRST READ)" - full writeup, toy kept below as precursor |
| `m5_9_4_results.json` (47^3), `m5_9_4_results_64.json` (63^3), `m5_9_5_results.json` + 3 PNGs | data + plots, all < 70K |
| `9a_lepton_mass_planning` | Status updated: full production engine build RAN locally |
| `checkpoints/00-04` | step-by-step progress |

## What it answers (Duda's review)

| Duda's point | Delivered |
| --- | --- |
| need full 4x4 with negative gravity + oscillation contributions | signed ledger: boost_GEM real + DOMINANT (3.0e-05 at 47^3, 9.3e-05 at 63^3, ~half the curvature; undressed control EXACTLY 0). clock_negative real but small at b*=0.13. |
| Faber is only the START of minimization | gradient-flow relaxation; dressed seed sits near a local min (<1%) |
| Higgs core-regularization, details to be found | confiner ~r0^3 confirmed, but production LdG V does NOT robustly select a scale (energy collapses to grid floor) - the open problem pinned precisely |
| no idea what parameters; files too simple | full TOPOLOGY_SEED block dumped to JSON; production engine, not toys; env-reproducible |

## Honest open frontier (NOT closed, by design)
- the scale-selecting Higgs potential form (production V is V-dominated, not scale-covariant)
- the eigenvalue hierarchy origin (why 206.8 / 3477) stays Yukawa input - delta knob gives 158x non-monotonic, not the spectrum
- #220 dynamical omega(r0): blocked by (clock period >> CFL window) + (constrained integrator instability at long t); analytic scale-covariance principle stands, REFINED (production V non-covariant)

## Bar check
DoD said matching mu/tau in one pass was never guaranteed; the deliverable is the full documented
framework + the potential search, honest about what lands. Delivered exactly that with real
production-engine runs (47^3 + 63^3 confirmation), the load-bearing new physics (the signed ledger),
real minimization, and an honest frontier map. Did not fake a spectrum match.

## For the maintainer's Duda reply (content bullets ready in 9a_lepton_mass_planning § "Open questions for Dr. Duda")
The open-questions bullets are staged; after the maintainer commits + publishes the distilled #200 writeup,
the reply to Duda can point at the signed ledger (the negative gravity term) + the honest open frontier.
