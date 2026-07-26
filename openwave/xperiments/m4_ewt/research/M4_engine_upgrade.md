# M4 EWT Engine Upgrade Plan

| Field | Value |
| --- | --- |
| Tracking issue | [#203](https://github.com/openwave-labs/openwave/issues/203) (M4 substrate: vector-field, non-linear PDE solver) |
| Unblocks | [#201](https://github.com/openwave-labs/openwave/issues/201) (K-selectivity, needs the non-linearity) |
| Scope of this cut | Vector + non-linear single-file PDE engine. NO vector calculus, NO director field, NO glyphs. |

## 1. Goal

Replace M4's analytical **Weighted Partial Standing Wave** core (`compute_voxel_wave()`, a closed-form per-source sum recomputed every step) with a real **time-integrated, non-linear, vector PDE solver**, built as a **single engine file** in M2's style and borrowing M5's non-linear potential and seed-wave insights.

The field stays a single 3-component vector `psi`. The engine evolves it by leapfrog over **all** voxels, hosts a swappable non-linear term `V(psi)`, and re-imposes wave centers through selectable interaction modes. Trackers, energy, force/motion, particle layer, scalar flux mesh, and granule render are kept as they are today.

## 2. Current state to target state

| Aspect | M4 today | M4 after upgrade |
| --- | --- | --- |
| Core | Analytical WPSW: `compute_voxel_wave()` recomputed per step (`wave_engine.py:104-293`) | Time-integrated leapfrog PDE: `propagate_wave()` |
| Voxel coverage | Selective: `select_voxels()` populates 7 voxels/center (`wave_engine.py:31-102`) | All voxels (a PDE must run everywhere); `select_voxels()` deleted |
| Field | `displacement_am`, single buffer (`medium.py:115`) | `psi_am` + triple buffer `psi_prev_am`, `psi_new_am` |
| Dynamics | Linear, monochromatic, closed-form | Non-linear, time-evolving, swappable `V(psi)` |
| Wave modes | Radial/longitudinal only (transverse = TODO) | Single `psi` vector (no L/T split this cut) |
| Seeding | Sources hardcoded into the per-step sum | `seed_wave()` @engine, called once @launcher |
| File layout | Single `wave_engine.py` | Single `wave_engine.py` (kept; NOT split like M5) |

## 3. Scope decisions (from the checklist)

| Do | Do NOT |
| --- | --- |
| Rename `displacement_am` to `psi_am` | Implement a director field |
| Add triple buffer `psi_prev_am`, `psi_new_am` | Implement vector calculus (div/curl L/T split) |
| Keep the single-file engine (M2-style, not M5's 4-file split) | Split the engine into multiple `.py` scripts |
| Remove `select_voxels()`; run the PDE over all voxels | Add new field observables to TRACKERS |
| Replace WPSW with a full vector PDE (M2 structure, M5 insights) | Add charging/damping controls to the seed |
| Add `seed_wave()` @engine + call @launcher | Display total-universe-energy dashboard (M2-style) |
| Add a non-linear term `V` with a swappable `V(psi)` | Create 2 wave modes (longitudinal/transverse) |
| Add WC interaction modes: Dirichlet, Neumann, + a third | Add a vector flux mesh |
| Keep trackers, energy, buffer swap, granule render, force/motion, particle as-is | Add glyph rendering |

## 4. The new wave engine in detail

All of the following lives in the single `m4_ewt/wave_engine.py`. The Taichi `@ti.kernel` / `@ti.func` split mirrors M2.

### 4.1 Field and buffers (MEDIUM)

Rename and promote the field to a triple buffer in `medium.py` (currently `medium.py:115`):

```text
psi_new_am   = ti.Vector.field(3, ti.f32, shape=grid)   # psi at t+dt
psi_am       = ti.Vector.field(3, ti.f32, shape=grid)   # psi at t      (renamed from displacement_am)
psi_prev_am  = ti.Vector.field(3, ti.f32, shape=grid)   # psi at t-dt
```

Trackers stay exactly as today (`medium.py:450-460`): `amp_local_emarms_am`, `freq_local_cross_rHz`, `last_crossing`, `energy_local_aJ`, plus the three global scalars. No new observable fields (the MEDIUM rule: keep trackers, do not add field observables, do not add the director / div / curl fields).

### 4.2 Seed wave = the EWT base wave (NOT wave-center sourced)

In EWT the medium carries an **always-on base wave**: a background ground-state oscillation that fills the whole medium. Particles (wave centers) are localized disturbances riding *on top of* that base wave, not its source. So `seed_wave()` seeds the **base wave**, and it is sourced from the **domain (universe) center**, not from the WC positions. WC sourcing is a separate concern handled by the WC interaction modes in P3 (re-driving at WC positions). This corrects the initial P1 design, which incorrectly summed radial waves from the WC positions.

`seed_wave(wave_field, seed_mode, boost, dt_rs)` is called once from the launcher before the main loop. It writes `psi_am` and `psi_prev_am` (the latter sets the initial velocity), radial (longitudinal) displacement `psi = A * profile(r) * r_hat` about the domain center, interior only (Dirichlet `psi = 0` at the boundary). It takes **no** `wave_center` argument. Charging/damping helpers (`compute_charge_envelope`, `compute_damping_factor`) are intentionally omitted.

| Seed mode | Profile (r = distance from domain center) | `psi_prev` | M2 analog |
| --- | --- | --- | --- |
| 0 `gaussian` | `A * exp(-r^2 / 2 sigma^2) * r_hat`, `sigma = lambda/2` | `psi_prev = psi` (rest start) | `charge_gaussian` |
| 1 `radial` | `A * exp(-r^2 / 2 sigma^2) * cos(k r) * r_hat` | `psi_prev = psi` (rest start) | (cosine-modulated) |
| 2 `full` | `A * cos(omega t - k r) * r_hat`, domain-filling outgoing wave | `psi_prev` at `t = -dt` (outgoing velocity) | `charge_full` |

Mode 2 is the closest emulation of the always-on base wave (a domain-filling radial oscillation). Modes 0 and 1 are localized pulses, useful for clean propagation / energy-conservation tests. Charge / particle identity (the per-WC `offset` phase, `particle.py:70`) is no longer applied at seed time: it belongs to the WC re-driving in P3, since the base wave itself is particle-independent.

### 4.3 Laplacian (componentwise, no vector calculus)

Vectorize M2's 6-point face-neighbor stencil (`compute_laplacianL`, `wave_engine.py:527`) by applying it to each of the 3 components independently. This is a plain vector Laplacian, **not** a div/curl decomposition, so it respects "no vector calculus":

```text
lap(psi)[c] = ( psi[i+1,c] + psi[i-1,c]
              + psi[i,j+1,c] + psi[i,j-1,c]
              + psi[i,j,k+1,c] + psi[i,j,k-1,c]
              - 6 * psi[i,j,k,c] ) / dx^2      for c in {x, y, z}
```

`compute_laplacian()` is a `@ti.func` returning a 3-vector, called only on interior voxels (Dirichlet psi=0 on the outer boundary, M2 convention).

### 4.4 Leapfrog evolve with non-linear term

`propagate_wave()` is the single `@ti.kernel` (M2's structure). Per interior voxel:

```text
force   = c^2 * compute_laplacian(psi_am) - dV_psi(psi_am, params)
psi_new = 2 * psi_am - psi_prev_am + (dt^2) * force
```

The linear part is exactly M2 (`propagate_wave`, `wave_engine.py:602`, update at `:647`). The `- dt^2 * dV_psi` term is the M5 pattern (`evolve_M`, `engine2_pde.py:376`, where force = `c^2 div(G) - dV_M`). CFL bound (M2): `dt <= dx / (c * sqrt(3))`.

### 4.5 Non-linear potential `V(psi)`, swappable

M5's `V_M` / `dV_M` act on matrix invariants `Tr(M^2)`, `Tr(M^3)` (`engine2_pde.py:286-320`). For a single vector field the natural scalar invariant is `u = ||psi||^2 = psi . psi`. Provide a `@ti.func` pair the user can edit, with a `V_MODE` selector and runtime coefficients (so potentials swap without recompiling the rest):

| `V_MODE` | `V(psi)` | Restoring force `dV_psi` | Coeffs | Tests |
| --- | --- | --- | --- | --- |
| 0 `linear` | 0 | 0 | none | linear sanity, energy conservation |
| 1 `cubic_nls` | `(c1/4) u^2` | `c1 * u * psi` | `c1 = k` | Smolinski psi^3 focusing (#201) |
| 2 `saturating` | `(c1/4) u^2 - (c2/6) u^3` | `c1 u psi - c2 u^2 psi` | `c1=k, c2=q` | soliton stabilizer (focus + defocus) |
| 3 `double_well` | `-(c1/2) u + (c2/4) u^2` | `-c1 psi + c2 u psi` | `c1=a, c2=b` | symmetry-broken vacuum |

`dV_psi(psi, v_mode, c1, c2)` returns a 3-vector; the two coefficients are passed at kernel-call time (each mode uses at most two), matching M5's runtime-argument style. `V_psi(...)` returns the matching scalar potential for offline energy diagnostics. Mode 0 reproduces M2's linear wave exactly, the regression baseline. **Focusing sign:** `c1 < 0` is focusing (self-amplifying), `c1 > 0` defocusing; confirmed empirically (P2). The leapfrog enters the force as `- dt^2 * dV_psi`.

### 4.6 Trackers, energy, buffer swap (kept, M2-placed)

Compute these **inside** `propagate_wave()` after the update, the way M2 does (`wave_engine.py:700-748`), but using M4's existing fields and formulas (the checklist: "keep trackers and energy as they are now"):

| Quantity | Method | Source pattern |
| --- | --- | --- |
| Amplitude | EMA of `\|\|psi\|\|` into `amp_local_emarms_am` | M2 EMA, alpha = 0.005 |
| Frequency | zero-crossing of a chosen component into `freq_local_cross_rHz` (now measurable, was hardcoded) | M2 zero-crossing, alpha = 0.05 |
| Energy | keep M4's per-voxel `energy_local_aJ = rho * V * (f * A)^2` (`wave_engine.py:284-293`) | M4 as-is |
| Buffer swap | `psi_prev <- psi`, `psi <- psi_new`, at end of kernel | M2 `wave_engine.py:743-748` |
| Global avg | `sample_avg_trackers()` 3-plane sampling, kept | M4 `wave_engine.py:442-505` |

Note: `energy_local_aJ` stays the existing heuristic so `force_motion.py` (which reads it for `F = -grad E`) is untouched. Strict PDE energy-conservation checking (true Hamiltonian density) is an **offline validation** in `validations/`, not a new shipped observable.

### 4.7 Wave center interaction modes

After the leapfrog update (`psi_am` now current), re-impose each active wave center through a selectable mode in a ball of `radius` voxels around the WC (the WC center is a node; the outer Dirichlet boundary is untouched). Positions come from `particle.py`; motion still flows through `force_motion.py` unchanged. Launcher config: `WC_INTERACT_MODE` (0 free / 1 dirichlet / 2 neumann / 3 soft), `WC_BOOST`, `WC_RADIUS`, `WC_SIGMA`.

| Mode | Function | Action in the WC ball | Physics it tests |
| --- | --- | --- | --- |
| 1 Dirichlet | `interact_wc_dirichlet` | hard-pin `psi = A sin(omega t + off) r_hat` (and `psi_prev`) | standing/reflective driven antenna; energy builds at the WC |
| 2 Neumann | `interact_wc_neumann` | bounded outgoing pin `psi = A sin(omega t + off - k r) r_hat` | radiating flux source (net energy outward); lower peak than Dirichlet |
| 3 Soft | `interact_wc_soft` | additive `psi += A exp(-r^2/2 sigma^2) cos(omega t + off) r_hat` | continuous radiator with full back-reaction (does not overwrite psi) |

A fourth mode (`0`) is **free**: no re-drive (seed once, evolve freely), the purest #201 test of whether the non-linearity self-selects a stable soliton without forcing. All four are selectable for A/B testing.

**Neumann finding (P3):** a strict velocity-Neumann pin (fixing `(psi - psi_prev)/dt`) is an unbounded integrator at the CFL `dt` (the drive period is under-resolved) and blows up. Mode 2 is therefore implemented as a value-bounded **outgoing radiating pin** (the `-k r` spatial phase sends energy outward, the flux-source role) rather than a raw velocity pin. A true velocity-Neumann would need a sub-CFL `dt` (the oscillon-search route, now tracked in [#201](https://github.com/openwave-labs/openwave/issues/201)).

## 5. Per-module work breakdown

### MEDIUM (`medium.py`)

| Task | Detail |
| --- | --- |
| Rename | `displacement_am` to `psi_am` (`:115`) |
| Add buffers | `psi_prev_am`, `psi_new_am` (same type/shape) |
| Trackers | keep `medium.py:450-460` unchanged; add no observables |

### WAVE_ENGINE (`wave_engine.py`, stays single file)

| Task | Detail |
| --- | --- |
| Delete | `select_voxels()` (`:31-102`), `compute_voxel_wave()` (`:104-293`), `propagate_wave_neighbors()` (`:297`), `propagate_wave_full()` (`:340`) |
| Add | `seed_wave()` (4.2), `compute_laplacian()` vector (4.3), `propagate_wave()` leapfrog + non-linear + trackers + energy + swap (4.4-4.6) |
| Add | `V_psi()` / `dV_psi()` with `V_MODE` (4.5) |
| Add | `interact_wc_dirichlet` / `interact_wc_neumann` / `interact_wc_soft` (4.7) |
| Keep | `update_flux_mesh_values()` (`:512-672`), `sample_position_to_render()` (`:679-704`), `sample_avg_trackers()` (`:442-505`) |

### FORCE & MOTION / PARTICLE (`force_motion.py`, `particle.py`)

| Task | Detail |
| --- | --- |
| Keep as-is | `compute_force_vector()` (`F = -grad E` from `energy_local_aJ`), `integrate_motion_leapfrog()`, `detect_annihilation()`, `WaveCenter` |
| Integration check | new `propagate_wave()` must still populate `energy_local_aJ` every step so `F = -grad E` keeps working |

### LAUNCHER (`_launcher.py`)

| Task | Detail |
| --- | --- |
| Seed call | call `seed_wave()` once before the main loop (replaces hardcoded source setup) |
| Main loop | `propagate_wave()` then `interact_wc_*()` then `compute_force_motion()` then advance time (was `compute_wave_oscillation` -> `compute_force_motion`, `:607-611`) |
| Single psi | drop any longitudinal/transverse mode branching; one `psi` only |
| Dashboard | do NOT add M2's total-universe-energy display (`charge_level` / `energy_total`) |
| Review model specs | audit resolution and stability (below) |

Model-level specs to review (a PDE has requirements the analytical engine did not):

| Spec | Why it now matters |
| --- | --- |
| Voxels per wavelength | PDE needs >= ~10 vox/lambda to resolve the wave; analytical WPSW needed none |
| Timestep | enforce CFL `dt <= dx / (c sqrt(3))` |
| Wave speed `c`, base `lambda`, base `freq` | must be mutually consistent with `c = lambda * freq` and the grid |
| Boundary | confirm Dirichlet psi=0 (or absorbing) is acceptable for the box size |
| Num sources / positions / offsets | unchanged inputs, but seeded via `seed_wave()` now |

### VISUALIZATION

| Task | Detail |
| --- | --- |
| Flux mesh | keep scalar mesh (`update_flux_mesh_values`), do NOT add a vector mesh |
| Granule render | keep vector granule render (`sample_position_to_render`) |
| Glyphs | do NOT add glyph rendering |

## 6. Phased implementation

| Phase | Work | Validation gate |
| --- | --- | --- |
| ✅ P0 Rename + buffers | `psi_am` + `psi_prev_am` + `psi_new_am`; no behavior change | engine still runs identically on `psi_am` (headless smoke: 19³ grid, max \|ψ\| ≈ 19.4 am, energy populated, prev buffer 0.0) |
| ✅ P1 Linear vector PDE | deleted WPSW + `select_voxels`; added `seed_wave`, `compute_laplacian`, `propagate_wave` (V off), swap; wired launcher seed + CFL `dt_rs`/`c_amrs` | headless (40³, CFL²=0.30): no blowup, energy E/E0 ∈ [1.03, 1.07] flat over 80 steps, wavefront expands. GUI visual pass pending (the maintainer) |
| ✅ P2 Non-linear V (machinery) | added `V_psi`/`dV_psi` + `V_MODE` (4 modes, 2 coeffs); injected `-dt^2 dV`; launcher `V_MODE`/`V_C1`/`V_C2` (default off) | linear mode (`V_MODE=0`) == M2 baseline ✅; non-linear term active + sign-correct (focusing `c1<0` preserves peak vs defocusing `c1>0`) ✅. A stable soliton core was NOT obtained from a released-gaussian seed (pure cubic collapses in 3D; saturating quintic does not arrest it at the large CFL `dt`) → the soliton search is tracked in [#201](https://github.com/openwave-labs/openwave/issues/201) (oscillon route) |
| ✅ P3 WC interaction modes | `interact_wc_dirichlet` / `_neumann` (bounded radiating pin) / `_soft` + free mode; launcher `WC_INTERACT_MODE`/`WC_BOOST`/`WC_RADIUS`/`WC_SIGMA` | headless (zero base field, 1 off-center WC): free=quiet, all 3 drives sustain the WC and stay bounded; dirichlet reflective (peak 3.1) vs neumann radiating (more E_tot, peak 1.8) vs soft back-reacting. Strict velocity-Neumann blew up → replaced by the bounded radiating pin |
| ✅ P4 Launcher + viz + specs + cleanup | single psi, seed calls, no energy dashboard, no glyphs, scalar mesh + granule kept; removed dead `selected_voxels`; `TIMESTEP`→`SIM_SPEED`; dashboard shows dt/CFL². `SEED_`/`V_`/`WC_` configs promoted to a per-xperiment `engine` xparameter section (2026-06-18) | end-to-end ✅ (headless: seed→propagate→WC drive→force→motion; WC moved 2.9 vox, field bounded) |

The M4 substrate upgrade is P0-P4 (done). The downstream physics (single-oscillon search → golden-angle K-selectivity) is **research, tracked in OpenWave issue [#201](https://github.com/openwave-labs/openwave/issues/201)** (see the [scoping comment](https://github.com/openwave-labs/openwave/issues/201#issuecomment-4746612842)), not a phase of this engineering plan.

### Post-P4 cleanup backlog

All items are done. The xparameters migration (deferred at P4) was completed 2026-06-18 once the M4 substrate froze: the `SEED_`/`V_`/`WC_` configs now live in each xperiment's `engine` section, read by the launcher via an `ENGINE_DEFAULTS` fallback.

| Item | Where | Action |
| --- | --- | --- |
| ✅ Dead `selected_voxels` machinery | `medium.py` (`max_selected_voxels`, `selected_voxels`, `num_selected_voxels`) | removed (P4) |
| ✅ Vestigial `TIMESTEP` slider | `_launcher.py` | repurposed to a `SIM_SPEED` control (0.5-1.0, live `c_amrs`); dashboard shows `dt_rs` + CFL² (P4) |
| ✅ Seed config constants | `_launcher.py` → xparameter `engine` section (`SEED_MODE`, `SEED_BOOST`) | promoted to per-xperiment xparameters |
| ✅ Potential config constants | `_launcher.py` → xparameter `engine` section (`V_MODE`, `V_C1`, `V_C2`) | promoted to per-xperiment xparameters |
| ✅ WC interaction config constants | `_launcher.py` → xparameter `engine` section (`WC_INTERACT_MODE`, `WC_BOOST`, `WC_RADIUS`, `WC_SIGMA`) | promoted to per-xperiment xparameters |

## 7. Acceptance criteria

| Criterion | Target |
| --- | --- |
| Linear regression | `V_MODE=0` reproduces an M2-class free wave (dispersion-free propagation at `c`) |
| Stability | CFL-respecting runs show no blow-up over a long horizon |
| Energy (linear) | total energy drift small over the run (offline Hamiltonian check in `validations/`) |
| Interaction modes | all four WC modes (Dirichlet / Neumann / soft / free) run and are A/B-selectable |

## 8. Open questions and risks

| Item | Note |
| --- | --- |
| `V(psi)` coefficients + soliton (P2 finding) | machinery is wired and sign-correct, but no stable core formed from a released-gaussian seed: weak `c1` is ~linear, stronger `c1<0` collapses (3D cubic instability), and the saturating quintic did not arrest it at the large CFL `dt` (`dt²≈523` amplifies the stiff term faster than the quintic catches it). The soliton search is tracked in [#201](https://github.com/openwave-labs/openwave/issues/201) (oscillon route): a localized seed profile, a likely sub-CFL `dt`, and a mass term, a time-periodic oscillon rather than a static soliton. |
| Frequency tracker | zero-crossing on a vector field needs a chosen component or `\|\|psi\|\|` extremum convention; pick and document one |
| Boundary reflections | Dirichlet box may reflect radiated energy back onto WCs; if it pollutes #201, add a simple absorbing margin (still no vector calculus) |
| #202 still blocked | emergent Coulomb needs the L/T div/curl split this cut omits; it is a deliberate follow-up, not delivered here |
| Discipline | port the solver, not the model: keep EWT physics (wave centers, base + disturbance), do not drift into M5's Landau-de Gennes ontology |

## 9. References

| Source | Anchor |
| --- | --- |
| Upgrade proposal | issue [#203](https://github.com/openwave-labs/openwave/issues/203) |
| Open problems | [#201](https://github.com/openwave-labs/openwave/issues/201) K-selectivity, [#202](https://github.com/openwave-labs/openwave/issues/202) Coulomb |
| M4 current core | `m4_ewt/wave_engine.py:31-102` (`select_voxels`), `:104-293` (`compute_voxel_wave`) |
| M4 field | `m4_ewt/medium.py:115` (`displacement_am`), trackers `:450-460` |
| M2 PDE template | `m2_free_wave/wave_engine.py:527` (laplacian), `:602` (`propagate_wave`), `:743-748` (swap), `:30/93` (`charge_full`/`charge_gaussian`) |
| M5 non-linear + seed | `m5_liquid_crystal/engine2_pde.py:286-320` (`V_M`/`dV_M`), `:376` (leapfrog force), `engine1_seeds.py:80-146` (`seed_hedgehog_M`) |
| M5 rendering catalog | `m5_liquid_crystal/research/4b_rendering_features.md` (reviewed; glyphs / vector mesh intentionally NOT ported) |
