# M4.10.2 Emergent Gravitational Force from Geometric Amplitude

> Extension artifact: Enhanced EWT, Łukasz Smoliński
> (manuscript v5.0.0, DOI 10.5281/zenodo.22540635)

## Criterion
`Gravity: Newton limit (GEM)` — strength (G) clause:
attractive 1/r² between masses, via the GEM route; the strength (G)
credited only when it follows from the model's own mechanism, never fitted.

## Status
Strength clause: circularity discharged by the M4.7 chain; accuracy not met,
`G_geom` residual 0.048169 %, 21.9x the CODATA uncertainty on `G`.

## What was computed

The gravitational force between the Sun and the Earth was computed from
geometric factors alone, with no G value entering the calculation.

The force kernel follows from the overlap of the two EMC density
deficits. Each soliton sources a monopole deficit δη_i(r) = −A_i/r
whose gradient is ∇(δη_i) = (A_i/r²) r̂. The angular integral of
∇(δη_1)·∇(δη_2) over the sphere of radius r gives exactly 4π/r²
for r ≥ R and zero for r < R (Newton's shell theorem analogue). The
radial integral over r ∈ [R, ∞) gives I(R) = 4π A_1 A_2 / R. The force
is −dI/dR = 4π A_1 A_2 / R², and the EMC field-energy normalization
K_emc turns the kernel into the physical force:

    F_grav = 4 π K_emc A_1 A_2 / R^2

The amplitudes and coupling are written in M4.7 chain factors:

    A_i   = (2 M_i r_e / m_e) * sqrt(X_eff)
            / (A_pi^4 N_geom^3 K_WC sqrt(N_nu_stat))

    K_emc = c^2 m_e A_pi^4 N_geom^3 K_WC sqrt(N_nu_eff) / (16 pi r_e)

Every factor on the right-hand side is either geometric (A_pi, N_geom,
K_WC, X_eff, N_nu_stat, N_nu_eff) or a dimensional anchor (r_e, m_e, c).
No G value appears in any of these expressions.

Derivation source: M4.10 artifact, sections "Angular Integral Correction
& Field Overlap Formulation" and "Physical Interaction Energy and Sign
Convention". Amplitude chain factors: M4.12 artifact.

## Result

    F_grav (geometric, no G) = 3.544208708919220e+22 N
    F_obs  (from G_CODATA)   = 3.542499577248244e+22 N
    relative difference      = 0.048246 %
    residual / CODATA unc.   = 21.9 x

With A_i and K_emc as written, A_i = 2 G_geom M_i / c^2 and
K_emc = c^4 / (16 pi G_geom), where G_geom is the M4.7 chain value
(`G_EWT` in `gravity_sector`). So F_grav = G_geom M_1 M_2 / R^2 holds as
an identity for every M_1, M_2 and R, and the comparison with F_obs is
the comparison of G_geom with G_CODATA. With the analytic kernel the
residual is 0.048169 %; the printed 0.048246 % adds the radial
quadrature's 7.7e-7.

The briefing's [reading of the strength clause](../../__M4_model_briefing.md#reading-the-two-gravity-criteria) splits it in two
halves. Circularity is discharged by M4.7 deriving lambda_l instead of
setting it to the Planck length. Accuracy asks for agreement within the
target's own uncertainty: at 21.9x the CODATA uncertainty on G, the
accuracy half is not met.

## Closed-form exponents (fixed chain outputs)

| Perturbation    | Chain fixed (section 8) | Through the chain |
|-----------------|-------------------------|-------------------|
| K_WC: 10 -> 9   | 11.11 %                 | 21.00 %           |
| N_geom: * 1.001 |  0.30 %                 |  1.19 %           |
| A_pi: * 1.01    |  3.90 %                 | not propagated    |
| r_e: * 1.01     |  1.00 %                 |  7.21 %           |

Section 8 rescales the final force by the exponent the closed form gives
for each factor, holding lambda_l, N_nu_stat, N_nu_eff and X_eff fixed.
It does not rerun the chain, so its rows cannot fail. At fixed chain
outputs each amplitude scales as r_e and K_emc as 1/r_e, so
`F_grav ~ K_emc * A_1 * A_2 ~ r_e`. Through the chain lambda_l also moves:
the third column recomputes lambda_l, N_nu_stat and X_eff from the
perturbed input, and there the force goes as r_e^7 and N_geom^-12. The
A_pi row is not propagated, since A_pi also enters alpha. The third
column comes from the reproduction in the [PR #602 review](https://github.com/openwave-labs/openwave/pull/602#pullrequestreview-5330379192).

Sensitivity does not separate a fitted constant from a derived one: both
move the result when perturbed. What separates them is where the
constant came from, which is the circularity half of the clause.

## Criterion mapping

| Clause | Evidence |
|--------|----------|
| attractive 1/r² between masses | F_geom = 4 π A_1 A_2 / R^2, R² dependence confirmed |
| via the GEM route | A_i, K_emc from BCC lattice geometry (M4.7 chain) |
| strength (G) follows from mechanism | Circularity discharged: M4.7 derives lambda_l. Accuracy not met: F_grav = G_geom M_1 M_2 / R^2 identically, residual 0.048169 %, 21.9x the CODATA uncertainty |
| never fitted | Rests on provenance (the circularity half), not on section 8, whose rows are closed-form exponents |

## Relation to M4.10

M4.10 and M4.12 compared G_geom with itself: G enters both A and K_emc
and cancels identically, so those gates confirm normalization
consistency. This artifact compares G_geom with G_CODATA, so unlike them
it can fail. It adds no factor beyond G_geom.

## Artifacts

- `research/scripts/m4_10_2_emergent_gravitational_force_geometric.py`

## Reference

Enhanced EWT manuscript, version 5.0.0 or later:
[DOI: 10.5281/zenodo.22540635](https://doi.org/10.5281/zenodo.22540635)

Relevant manuscript section:
- "Newtonian Force from Interacting EMC Deficits"

Relevant M4 artifacts:
- M4.10 — angular integral, sign convention, K_emc normalization
- M4.12 — amplitude written in M4.7 chain factors
- M4.7  — geometric primitives (A_pi, N_geom, K_WC, X_eff,
  N_nu_stat, N_nu_eff)
