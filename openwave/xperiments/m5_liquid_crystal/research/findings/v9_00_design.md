# #200 full production engine build - design + progress log (LOCAL, do not publish)

Go-time: 2026-06-20 08:31 EDT. Reset 12:10pm EDT (resume ping armed 16:15 UTC).
Keeping everything LOCAL per the maintainer: no GitHub status, no comments. Progress here.

## The corrected pipeline (from 9a_lepton_mass_planning, Duda's critique)

The toy (m5_9_1/2/3) was 3x3, fixed Faber ansatz, no `g` axis, no clock. The full production engine build
runs on the PRODUCTION 4x4 Taichi engine and adds the NEGATIVE energy contributions the 3x3
structurally could not see.

## Load-bearing physics decision (verified by reading the engine, 2026-06-20)

`engine3_observables.compute_energyH_density_M` computes the EUCLIDEAN Hamiltonian:
`H = 1/2\|\|Mdot\|\|^2 + c^2*4*sum\|\|[M_mu,M_nu]\|\|_F^2 + (V - v0)`. The curvature uses the
Frobenius norm => always POSITIVE => it cannot represent the negative gravity contribution.

The engine's DYNAMICS, however, are Minkowski-signed (`engine2_pde.signed_dot4`,
`eta_twist_masked`, eta = diag(1,1,1,-1)): the `(alpha,3)` / `(3,alpha)` components of any
matrix enter with a MINUS sign. So the physically-conserved curvature energy is
`c^2*4*sum signed_dot4([M_mu,M_nu],[M_mu,M_nu])`. When the boost dressing puts weight on the
`(alpha,3)` (time-axis) components, the SIGNED curvature is LOWER than the Euclidean one. That
gap IS the boost-GEM negative gravity contribution Duda asked for:

  GEM dip = H_euclid - H_signed = 2 * c^2 * 4 * sum_(alpha,3 comps) [M_mu,M_nu]^2  (>= 0, lowers E)

The new instrumentation = a SIGNED energy kernel mirroring `compute_energyH_density_M` but using
`signed_dot4` for the curvature, reporting kinetic / curv_euclid / curv_signed / V separately so
the ledger reads `m_rest = (gradient + V) - boost_GEM - clock_activation`.

## Sub-task plan + status

| # | Sub-task | Script | Status |
| --- | --- | --- | --- |
| 1 | Stand up production 4x4 engine headless; one dressed hedgehog; reproduce electron rest energy by minimization | m5_9_4_engine_lepton.py | 🔶 building |
| 2 | Negative-energy instrumentation: signed/GEM energy + clock-activation separation | m5_9_4 (same) | 🔶 building |
| 3 | True energy-minimization mode (gradient flow / damped dynamics on the full 4x4 functional) | m5_9_4 (same) | 🔶 building |
| 4 | Higgs-potential search: fixed-density confiner that selects the core scale + stabilizes biaxiality | m5_9_5_higgs_search.py | 🚧 next |
| 5 | Three eigen-configurations -> masses -> compare to 207 / 3477 | m5_9_6_lepton_spectrum.py | 🚧 next |
| 6 | Documentation + convincing writeup (params, potential, minimizing configs, energy ledger, plots) | findings doc | 🚧 next |

## Engine API (verified, from the headless template m5_8_1_headless_check.py)

- `tf = medium.TensorField([EDGE,EDGE,EDGE], target_voxels, [0.5,0.5,0.5], viz_stride=2)`; EDGE=1e-15 m.
- `seeds.seed_dressed_hedgehog_M(tf, cx,cy,cz, r0_vox, rhoc_vox, delta, b_star, rw_vox, kick_theta)`.
- `pde.compute_stable_mask(tf)`; `pde.compute_tstar(tf)`.
- leapfrog 4D: `compute_curvature_flux_4d` -> `sample_v03_drift` -> `evolve_M_4d(tf,c,dt,ldg_c,vm0,vm1,vm2,km)` -> `swap_matrix_buffers`.
- constrained: `flux_4d_constrained`->`update_P_4d`->`sample_p03_drift`->`apply_p03_clamp`->`solve_constrained_4d`->`update_M_4d_constrained`->swap; `init_P_4d(tf,kick)` first.
- energy: `obs.compute_energyH_density_M(tf, observables, c_amrs, dt_rs, a,b,c, v0, e_scale)` -> `observables.energyH_density_aJ`.
- `obs = medium.FieldObservables(tf.grid_size)`; `trackers = medium.Trackers(tf.grid_size)`.
- LdG: `ldg_c = K*c_amrs^2/dx_am^4`, `ldg_a = -2*ldg_c*(1+delta^2)`, `ldg_b=0`, `v0 = ldg_a*(1+d^2)+ldg_c*(1+d^2)^2`.
- Minimization mode used: set `M_prev <- M` each step before evolve => leapfrog reduces to pure
  gradient-descent `M_new = M + dt^2*Force/inertia` (overdamped gradient flow on the full 4x4 energy).
