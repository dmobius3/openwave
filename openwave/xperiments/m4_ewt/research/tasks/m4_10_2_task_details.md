# M4.10.2 Task Details

## Objective
Address the strength (G) clause of the `Gravity: Newton limit (GEM)`
criterion.

## The criterion, precisely

> "attractive 1/r² between masses, via the GEM route; the strength (G)
> credited only when it follows from the model's own mechanism, never
> fitted"

The clause does not say "do not use G". It says the strength must
**follow from the mechanism**, not be **fitted**. G must be a
consequence of the model, not an input to it.

## Background
M4.10 computes F_EMC and F_Newton with the same G_geom, so G cancels
identically on both sides. The gate confirms normalization consistency
but does not demonstrate that G follows from the mechanism: G enters
as an input to A and K_emc, not as a derived quantity of the force
itself. The strength clause therefore cannot be credited on M4.10 alone.

## Method
1. Express A_i and K_emc in M4.7 chain factors, so that G does not
   appear anywhere in the force calculation.
2. Compute F_grav = 4 pi K_emc A_1 A_2 / R^2 from these factors.
3. Compare with F_obs = G_CODATA M1 M2 / R^2.
4. Mutation tests on geometric inputs (K_WC, N_geom, A_pi, r_e).

The absence of G in step 1 is not the goal. It is the consequence of
the goal: a force whose strength follows from the model's own
mechanism. If the mechanism is correct, G is not needed as an input;
the strength is already fixed by the geometry.

## Acceptance
- No G appears as an input in the force calculation.
- F_grav is compared with F_obs within the CODATA uncertainty on G.
  Since F_grav = G_geom M_1 M_2 / R^2 identically, this is the G_geom
  accuracy statement (result: not met, 21.9x).
- The section 8 rows report the closed-form exponents at fixed chain
  outputs. They cannot show the absence of a fitted coupling; that rests
  on provenance, M4.7 deriving lambda_l.
- G recovered from F_grav agrees with the M4.7 chain (structural
  consistency, not a test).

## Dependencies
- `m4_7_ewt_emergence_engine.py` (M4.7 chain).

## Artifacts
- `research/scripts/m4_10_2_emergent_gravitational_force_geometric.py`
- `research/findings/m4_10_2_emergent_gravitational_force_geometric.md`
