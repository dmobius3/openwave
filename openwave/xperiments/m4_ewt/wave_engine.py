"""
ENERGY-WAVE ENGINE

ON VECTOR-WAVE MODEL

Wave Physics Engine @spacetime module. Wave dynamics and motion.
"""

import taichi as ti
import numpy as np

from openwave.common import colormap, constants

# ================================================================
# Energy-Wave Oscillation Parameters
# ================================================================
base_amplitude_am = constants.EWAVE_AMPLITUDE / constants.ATTOMETER  # am, oscillation amplitude
base_wavelength = constants.EWAVE_LENGTH  # in meters
base_wavelength_am = constants.EWAVE_LENGTH / constants.ATTOMETER  # in attometers
base_frequency = constants.EWAVE_FREQUENCY  # in Hz
base_frequency_rHz = constants.EWAVE_FREQUENCY * constants.RONTOSECOND  # in rHz (1/rontosecond)
rho = constants.MEDIUM_DENSITY  # medium density for Gaussian energy calc (kg/m³)
rho_qgam = constants.MEDIUM_DENSITY_QGAM  # qg/am³, for energy computation in scaled units


# ================================================================
# WAVE PROPAGATION ENGINE  (vector PDE solver)
# ================================================================
# Time-integrated leapfrog vector wave equation:  ∂²ψ/∂t² = c²∇²ψ
# (the nonlinear restoring term −dV(ψ) is added in P2).
#
# A single 3-component field ψ (psi_am) is evolved over ALL voxels — a PDE has
# no "selective voxel" shortcut. The structure mirrors M2's free-wave engine,
# vectorized componentwise; seeding takes the multi-center superposition insight
# from M5. NO vector calculus (div/curl) and NO director field (see #203 scope).
# ================================================================


@ti.kernel
def seed_wave(
    wave_field: ti.template(),  # type: ignore
    seed_mode: ti.i32,  # type: ignore
    boost: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
):
    """
    Seed the medium's BASE WAVE: the always-on background oscillation of the medium
    at its ground state (EWT). This is NOT sourced from the wave centers — in EWT the
    base wave fills the whole medium and the wave centers are localized disturbances on
    top of it, re-driven separately by the WC interaction modes (P3). The base wave
    radiates from the DOMAIN (universe) center voxel as a radial (longitudinal)
    displacement ψ = A·profile(r)·r̂. Interior only (Dirichlet ψ=0 at the outer boundary).

    seed_mode:
        0 = gaussian : A·exp(-r²/2σ²)·r̂ (σ = λ/2) localized pulse, released from rest
        1 = radial   : A·exp(-r²/2σ²)·cos(k·r)·r̂ cosine-modulated pulse, released from rest
        2 = full     : A·cos(ω·t − k·r)·r̂ domain-filling outgoing radial wave (à la M2
                       charge_full); ψ_prev is set at t = −dt to give an initial outgoing
                       velocity. Closest emulation of the always-on base wave.

    Args:
        wave_field: WaveField instance (psi_am / psi_prev_am / psi_new_am)
        seed_mode: base-wave profile selector
        boost: amplitude multiplier
        dt_rs: timestep (rs), used by mode 2 to set the t = −dt previous buffer
    """
    # Domain (universe) center, in grid indices — the base-wave source
    cx = wave_field.nx * 0.5
    cy = wave_field.ny * 0.5
    cz = wave_field.nz * 0.5

    # Wave numbers in grid/scaled units; Gaussian width σ = λ/2 (grid units)
    wavelength_grid = base_wavelength * wave_field.scale_factor / wave_field.dx
    k_grid = 2.0 * ti.math.pi / wavelength_grid  # radians per grid index
    omega_rs = 2.0 * ti.math.pi * base_frequency_rHz / wave_field.scale_factor  # rad/rs
    sigma_grid = wavelength_grid / 2.0
    amp = base_amplitude_am * wave_field.scale_factor * boost

    # Seed ψ (and ψ_prev) over interior voxels (Dirichlet ψ=0 at edges)
    for i, j, k in ti.ndrange(
        (1, wave_field.nx - 1), (1, wave_field.ny - 1), (1, wave_field.nz - 1)
    ):
        # Radial direction from the DOMAIN CENTER to this voxel
        dir_vec = ti.Vector(
            [ti.cast(i, ti.f32) - cx, ti.cast(j, ti.f32) - cy, ti.cast(k, ti.f32) - cz]
        )
        r_grid = dir_vec.norm()
        r_hat = dir_vec / (r_grid + 1e-6)
        node = ti.select(r_grid < 0.5, 0.0, 1.0)  # node at the exact center (r̂ undefined)
        envelope = ti.exp(-(r_grid * r_grid) / (2.0 * sigma_grid * sigma_grid))

        # Default: gaussian pulse released from rest (mode 0)
        psi = amp * envelope * r_hat * node
        psi_prev = psi
        if seed_mode == 1:  # cosine-modulated radial pulse, released from rest
            psi = amp * envelope * ti.cos(k_grid * r_grid) * r_hat * node
            psi_prev = psi
        elif seed_mode == 2:  # domain-filling outgoing base wave (à la charge_full)
            psi = amp * ti.cos(-k_grid * r_grid) * r_hat * node
            psi_prev = amp * ti.cos(omega_rs * (-dt_rs) - k_grid * r_grid) * r_hat * node

        wave_field.psi_am[i, j, k] = psi
        wave_field.psi_prev_am[i, j, k] = psi_prev

    # Clear the next buffer
    for i, j, k in ti.ndrange(wave_field.nx, wave_field.ny, wave_field.nz):
        wave_field.psi_new_am[i, j, k] = ti.Vector([0.0, 0.0, 0.0])


@ti.func
def compute_laplacian(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """
    Vector Laplacian ∇²ψ at voxel [i,j,k], 6-point face-neighbor stencil applied
    componentwise (a plain vector Laplacian, NO div/curl decomposition).
        ∇²ψ ≈ (face_sum - 6·center) / dx²

    Args:
        wave_field: WaveField instance (psi_am)
        i, j, k: interior voxel indices (0 < i,j,k < n-1)

    Returns:
        3-vector Laplacian in units [am / am²] = [1/am]
    """
    face_sum = (
        wave_field.psi_am[i + 1, j, k]
        + wave_field.psi_am[i - 1, j, k]
        + wave_field.psi_am[i, j + 1, k]
        + wave_field.psi_am[i, j - 1, k]
        + wave_field.psi_am[i, j, k + 1]
        + wave_field.psi_am[i, j, k - 1]
    )
    center = wave_field.psi_am[i, j, k]
    return (face_sum - 6.0 * center) / (wave_field.dx_am**2)


# ================================================================
# NON-LINEAR POTENTIAL V(ψ)  (swappable; acts on u = ‖ψ‖², NO vector calculus)
# ================================================================
# The dynamics use dV_psi (the restoring force, entered as −dt²·dV in the leapfrog).
# V_psi is the matching scalar potential, kept for offline energy diagnostics.
# Coefficients (c1, c2) are passed at call time so potentials swap without
# recompiling the rest.
#
#   v_mode  0 linear        : V = 0                                      dV = 0
#   v_mode  1 cubic_nls     : V = (c1/4)·u²                              dV = c1·u·ψ            (c1 = k)
#   v_mode  2 saturating    : V = (c1/4)·u² − (c2/6)·u³                  dV = c1·u·ψ − c2·u²·ψ  (c1 = k, c2 = q)
#   v_mode  3 double_well   : V = −(c1/2)·u + (c2/4)·u²                  dV = −c1·ψ + c2·u·ψ    (c1 = a, c2 = b)
#   v_mode  4 deficit       : V = (c1/4)·u²·mod(r)                       dV = c1·u·ψ·mod(r)     (c1 = γ, mod = 1−ρ/ρ₀, core deficit, no wall)
#   v_mode  5 deficit+wall  : V = (c1/4)·u²·mod(r)                       dV = c1·u·ψ·mod(r)     (c1 = γ, mod = 1−ρ/ρ₀, core deficit + EMC wall)
#   v_mode  6 flat          : V = (c1/4)·u²·mod_flat(r)                  dV = c1·u·ψ·mod_flat(r) (c1 = γ, constant deficit inside, sigmoid transition)
#   v_mode  7 flat+wall     : V = (c1/4)·u²·mod_flat_wall(r)             dV = c1·u·ψ·mod_flat_wall(r) (c1 = γ, flat deficit + sigmoid + wall)
#   v_mode  9 gaussian      : V = (c1/4)·u²·mod_gauss(r)                 dV = c1·u·ψ·mod_gauss(r) (c1 = γ, smooth Gaussian profile, no wall)
#   v_mode 10 gaussian+satur: V = (c1/4)·u²·mod_gauss(r) − (c2/6)·u³·mod_gauss(r)  dV = c1·u·ψ·mod_gauss(r) − c2·u²·ψ·mod_gauss(r) (c1 = γ, c2 = q)
#
# Legend:
#   u = ‖ψ‖²               = squared field magnitude (scalar invariant)
#   mod(r)                 = spatial modulation factor (density deficit profile)
#   mod_gauss(r)           = deficit_depth * exp(-(r / R_soliton)^2)
#   mod_flat(r)            = 1 - rho_flat(r)/rho_0  (sigmoid transition)
#   mod_flat_wall(r)       = 1 - rho_flat_wall(r)/rho_0  (flat + Gaussian wall)
#
# Modes 4-7 and 9-10 use a radial EMC density profile ρ(r) to modulate the nonlinearity.
# The modulation factor is active only where the density is depleted relative to
# the statutory background ρ₀ = N_nu_stat (vacuum equilibrium density).
#
# Mode 10 is the RECOMMENDED mode for electron soliton stability:
#   - Gaussian profile creates a smooth potential well (no sharp gradients)
#   - Quintic saturation prevents amplitude blow-up (physical limiting)
#   - Vacuum pressure (PRESSURE_STRENGTH) is active (V_MODE >= 4)
#
# For a focusing (attractive) potential, use c1 < 0.
# For saturation strength, c2 > 0 (quintic term) — larger values = stronger limiting.
# Formalization in doc: https://doi.org/10.5281/zenodo.21557369
# ================================================================


@ti.func
def V_psi(
    psi,  # type: ignore
    v_mode: ti.i32,  # type: ignore
    c1: ti.f32,  # type: ignore
    c2: ti.f32,  # type: ignore
    rho_local: ti.f32,  # type: ignore  # for modes 4/5 only
):
    """Scalar non-linear potential V(ψ) on u = ‖ψ‖² (see table above). Offline diagnostics only."""
    u = psi.dot(psi)
    out = 0.0
    if v_mode == 1:
        # V = (c1/4)·u²
        out = 0.25 * c1 * u * u
    elif v_mode == 2:
        # V = (c1/4)·u² − (c2/6)·u³  [A4] fixed: now returns scalar, not vector
        out = 0.25 * c1 * u * u - (c2 / 6.0) * u * u * u
    elif v_mode == 3:
        # V = −(c1/2)·u + (c2/4)·u²
        out = -0.5 * c1 * u + 0.25 * c2 * u * u
    elif v_mode in (4, 5):
        # Modulation depends on radial position, so we return 0.0 here.
        # The physical effect is captured in dV_psi.
        out = 0.0
    elif v_mode == 10:
        # Gaussian profile + quintic saturation
        # Returns 0 here; physical effect is in dV_psi
        out = 0.0
    return out


# ================================================================
# EMC DENSITY PROFILE (NEW for v_mode 4 and 5)
# ================================================================
# Provides a radial density profile ρ(r)/ρ₀ for the EMC medium.
#   v_mode 4: core deficit only (ρ rises from 1−depth at r=0 to 1.0 at ∞)
#   v_mode 5: core deficit + Gaussian EMC wall peak at r = r_wall
# ================================================================


@ti.func
def emc_density_profile(
    r: ti.f32,  # type: ignore
    r_wall: ti.f32,  # type: ignore
    wall_height: ti.f32,  # type: ignore
    deficit_depth: ti.f32,  # type: ignore
    v_mode: ti.i32,  # type: ignore
) -> ti.f32:
    """
    Radial EMC density profile rho(r) / N_nu_stat.

    Args:
        r: distance from soliton centre (grid units)
        r_wall: radius of the Degraded EMC Wall (grid units)
        wall_height: peak density at the wall, relative to background (e.g. 1.2 = +20%)
        deficit_depth: core density reduction (e.g. 0.8 = 80% below background)
        v_mode: 4 = deficit only, 5 = deficit + wall

    Returns:
        Normalised EMC density ρ(r) / ρ₀
    """
    rho0 = 1.0  # statutory background (normalised)

    if v_mode == 4:
        # Core deficit, smooth asymptotic return to background
        rho = rho0 - deficit_depth * ti.exp(-r / (r_wall * 0.3))
    elif v_mode == 5:
        # Core deficit + Gaussian wall peak
        rho_core = rho0 - deficit_depth * ti.exp(-r / (r_wall * 0.3))
        bump = (wall_height - 1.0) * ti.exp(-((r - r_wall) ** 2) / (2.0 * (r_wall * 0.1) ** 2))
        rho = rho_core + bump
    else:
        rho = rho0  # fallback (no modulation)
    return rho


@ti.func
def emc_density_profile_flat_with_wall(
    r: ti.f32,
    r_soliton: ti.f32,  # radius of the flat bottom
    sigma: ti.f32,  # transition width (sigmoid)
    deficit_depth: ti.f32,  # deficit depth at the centre
    r_wall: ti.f32,  # wall radius (centre of the bump)
    wall_height: ti.f32,  # wall height (rho_wall - rho0)
    wall_sigma: ti.f32,  # wall width
) -> ti.f32:
    """
    Flat-bottom density profile with outer wall.

    r < r_soliton:        rho = 1 - deficit_depth (constant deficit)
    r_soliton < r < r_wall: sigmoidal transition to rho = 1
    r ≈ r_wall:           Gaussian wall (rho > 1)
    """
    rho0 = 1.0

    # 1. Flat bottom inside soliton
    rho_inner = rho0 - deficit_depth

    # 2. Sigmoid transition from inner to vacuum
    # Smoothstep: 0 at r = r_soliton, 1 at r = r_wall (approx)
    # Using tanh for smooth transition
    transition = 0.5 * (1.0 + ti.tanh((r - r_soliton) / sigma))

    # 3. Outer wall (Gaussian bump)
    wall = (wall_height - 1.0) * ti.exp(-((r - r_wall) ** 2) / (2.0 * wall_sigma**2))

    # Combine: start with inner, transition to rho0, add wall
    rho = rho_inner + (rho0 - rho_inner) * transition + wall

    return rho


# ================================================================
# RESTORING FORCE dV/dψ (NEW spatial arguments for modes 4/5)
# ================================================================


@ti.func
def emc_density_profile_flat(
    r: ti.f32,
    r_soliton: ti.f32,
    sigma: ti.f32,
    deficit_depth: ti.f32,
) -> ti.f32:
    """
    Flat-bottom density profile: constant deficit inside soliton,
    rapid rise to vacuum outside via sigmoid.

    For r << r_soliton: rho ≈ 1 - deficit_depth (deficit)
    For r >> r_soliton: rho ≈ 1 (vacuum)
    """
    rho0 = 1.0
    # Sigmoid transition: smooth but steep rise at r = r_soliton
    rho = rho0 - deficit_depth * (0.5 * (1.0 - ti.tanh((r - r_soliton) / sigma)))
    return rho


@ti.func
def dV_psi(
    psi,  # type: ignore
    v_mode: ti.i32,  # type: ignore
    c1: ti.f32,  # type: ignore
    c2: ti.f32,  # type: ignore
    rho_local: ti.f32,  # type: ignore  # for modes 4/5 only
    # Spatial context for radial density profile (modes 4 and 5)
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
    nx: ti.i32,  # type: ignore
    ny: ti.i32,  # type: ignore
    nz: ti.i32,  # type: ignore
    r_wall: ti.f32,  # type: ignore
    wall_height: ti.f32,  # type: ignore
    deficit_depth: ti.f32,  # type: ignore
    r_soliton: ti.f32,
    sigma: ti.f32,
):
    """Restoring force dV/dψ (3-vector) for V_psi; enters the leapfrog as −dt²·dV_psi.

    v_mode = 4: density-modulated cubic with core deficit (no EMC wall).
    v_mode = 5: density-modulated cubic with core deficit + EMC wall.
    v_mode = 6: flat-bottom density profile (sigmoid transition, no wall).
    v_mode = 7: flat-bottom density profile + sigmoid + EMC wall.
    v_mode = 9: Gaussian density profile (smooth, no wall).
    v_mode = 10: Gaussian density profile + quintic saturation (RECOMMENDED).
    """
    u = psi.dot(psi)
    out = ti.Vector([0.0, 0.0, 0.0])

    if v_mode == 1:
        # Pure cubic nonlinearity (focusing if c1 < 0, defocusing if c1 > 0)
        out = c1 * u * psi
    elif v_mode == 2:
        # Quintic saturation: cubic focusing + quintic limiting
        # For small u: cubic term dominates (focusing)
        # For large u: quintic term dominates (saturation, prevents blow-up)
        out = c1 * u * psi - c2 * u * u * psi
    elif v_mode == 3:
        # Double-well potential: symmetry-broken vacuum
        out = -c1 * psi + c2 * u * psi
    elif v_mode == 4 or v_mode == 5:
        # V_MODE=4: Exponential density deficit (core only, no wall)
        # V_MODE=5: Exponential density deficit + EMC wall (Gaussian bump)
        # Compute radial distance from domain centre
        cx = ti.cast(nx, ti.f32) * 0.5
        cy = ti.cast(ny, ti.f32) * 0.5
        cz = ti.cast(nz, ti.f32) * 0.5
        dx = ti.cast(i, ti.f32) - cx
        dy = ti.cast(j, ti.f32) - cy
        dz = ti.cast(k, ti.f32) - cz
        r = ti.sqrt(dx * dx + dy * dy + dz * dz)

        # Get normalised EMC density at this voxel
        rho_norm = emc_density_profile(r, r_wall, wall_height, deficit_depth, v_mode)

        # Modulation: nonlinearity active only where density is below background
        modulation = 1.0 - ti.min(rho_norm, 1.0)
        out = c1 * u * psi * modulation
    elif v_mode == 6:
        # V_MODE=6: Flat-bottom density profile (sigmoid transition, no wall)
        # rho(r) = 1 - deficit * 0.5 * (1 - tanh((r - R_soliton) / sigma))
        # Constant deficit inside soliton, smooth transition to vacuum outside
        cx = ti.cast(nx, ti.f32) * 0.5
        cy = ti.cast(ny, ti.f32) * 0.5
        cz = ti.cast(nz, ti.f32) * 0.5
        dx = ti.cast(i, ti.f32) - cx
        dy = ti.cast(j, ti.f32) - cy
        dz = ti.cast(k, ti.f32) - cz
        r = ti.sqrt(dx * dx + dy * dy + dz * dz)

        # Get flat-bottom density profile
        rho_norm = emc_density_profile_flat(r, r_soliton, sigma, deficit_depth)

        # Modulation: 1 - rho/rho0 (active only where rho < rho0)
        modulation = 1.0 - ti.min(rho_norm, 1.0)

        # Nonlinear term: c1 * u * psi * modulation
        out = c1 * u * psi * modulation
    elif v_mode == 7:
        # V_MODE=7: Flat-bottom density profile + sigmoid + Gaussian EMC wall
        # Combines constant deficit inside, sigmoid transition to vacuum,
        # and a Gaussian wall (density bump) at r = r_wall
        cx = ti.cast(nx, ti.f32) * 0.5
        cy = ti.cast(ny, ti.f32) * 0.5
        cz = ti.cast(nz, ti.f32) * 0.5
        dx = ti.cast(i, ti.f32) - cx
        dy = ti.cast(j, ti.f32) - cy
        dz = ti.cast(k, ti.f32) - cz
        r = ti.sqrt(dx * dx + dy * dy + dz * dz)

        # Wall parameters for V_MODE=7:
        # wall_sigma controls the width of the Gaussian wall.
        # For simplicity we use wall_sigma = r_wall * 0.1 (constant).
        wall_sigma = r_wall * 0.1

        rho_norm = emc_density_profile_flat_with_wall(
            r, r_soliton, sigma, deficit_depth, r_wall, wall_height, wall_sigma
        )
        modulation = 1.0 - ti.min(rho_norm, 1.0)
        out = c1 * u * psi * modulation
    elif v_mode == 9:
        # V_MODE=9: Gaussian density profile (smooth, no wall)
        # rho(r) = 1 - deficit_depth * exp(-(r / R_soliton)^2)
        # Smooth, monotonic transition from core deficit to vacuum
        cx = ti.cast(nx, ti.f32) * 0.5
        cy = ti.cast(ny, ti.f32) * 0.5
        cz = ti.cast(nz, ti.f32) * 0.5
        dx = ti.cast(i, ti.f32) - cx
        dy = ti.cast(j, ti.f32) - cy
        dz = ti.cast(k, ti.f32) - cz
        r = ti.sqrt(dx * dx + dy * dy + dz * dz)

        rho = 1.0 - deficit_depth * ti.exp(-((r / r_soliton) ** 2))
        modulation = 1.0 - ti.min(rho, 1.0)
        out = c1 * u * psi * modulation
    elif v_mode == 10:
        # ------------------------------------------------------------------
        # V_MODE=10: Gaussian density profile + quintic saturation
        #
        # Combines the smooth Gaussian profile (V_MODE=9) with the
        # saturation mechanism (V_MODE=2). This creates a potential well
        # (from the density deficit) while preventing amplitude blow-up
        # (from the quintic term).
        #
        # Formula: out = c1 * u * psi * modulation - c2 * u * u * psi * modulation
        # where modulation = deficit_depth * exp(-(r / R_soliton)^2)
        #
        # This is the most physically complete V_MODE, as it includes:
        # - Vacuum density profile (Gaussian well)
        # - Nonlinear saturation (quintic)
        # - Pressure force (via PRESSURE_STRENGTH in force_motion)
        #
        # For electron (K=10) stability, recommended parameters:
        #   c1 < 0 (focusing), c2 > 0 (saturation), R_soliton ~ 35
        #
        # FUTURE EXTENSION:
        # The current implementation uses a fixed Gaussian profile.
        # A future enhancement could make the profile selectable via an
        # additional parameter (e.g., PROFILE_MODE), allowing V_MODE=10
        # to use any of the density profiles defined in V_MODE=4-7:
        #   - exponential (deficit only)
        #   - exponential + EMC wall
        #   - flat-bottom (sigmoid)
        #   - flat-bottom + wall
        # This would enable systematic comparison of different vacuum
        # density profiles while retaining the saturation mechanism.
        # Implementation: replace the modulation calculation below with
        # a call to a profile function selected by PROFILE_MODE.
        # ------------------------------------------------------------------
        cx = ti.cast(nx, ti.f32) * 0.5
        cy = ti.cast(ny, ti.f32) * 0.5
        cz = ti.cast(nz, ti.f32) * 0.5
        dx = ti.cast(i, ti.f32) - cx
        dy = ti.cast(j, ti.f32) - cy
        dz = ti.cast(k, ti.f32) - cz
        r = ti.sqrt(dx * dx + dy * dy + dz * dz)

        # Gaussian modulation: strongest at center, decays to zero
        modulation = deficit_depth * ti.exp(-((r / r_soliton) ** 2))

        # Cubic focusing + quintic saturation, both modulated by density profile
        out = c1 * u * psi * modulation - c2 * u * u * psi * modulation
    return out


@ti.kernel
def propagate_wave(
    wave_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    elapsed_t_rs: ti.f32,  # type: ignore
    v_mode: ti.i32,  # type: ignore
    v_c1: ti.f32,  # type: ignore
    v_c2: ti.f32,  # type: ignore
    # New profile parameters for modes 4 and 5
    r_wall: ti.f32,  # type: ignore
    wall_height: ti.f32,  # type: ignore
    deficit_depth: ti.f32,  # type: ignore
    r_soliton: ti.f32,
    sigma: ti.f32,
):
    """
    Evolve ψ one step by leapfrog over all interior voxels (non-linear vector wave PDE).

    Continuous PDE: ∂²ψ/∂t² = c²∇²ψ − dV(ψ)
    Discrete form (Verlet / leapfrog):
        ψ_new = 2ψ - ψ_prev + (c·dt)²·∇²ψ − dt²·dV(ψ)
    With v_mode = 0 the dV term is zero and this is exactly the linear wave (M2 baseline).

    Boundary: Dirichlet ψ=0 at edges (boundary voxels skipped).
    Trackers (RMS amplitude, zero-crossing frequency) and per-voxel energy are
    updated in-kernel using M4's existing fields/formulas, then buffers are swapped.

    Args:
        wave_field: WaveField (psi_am / psi_prev_am / psi_new_am)
        trackers: WaveTrackers (amp / freq / energy fields)
        c_amrs: wave speed (am/rs)
        dt_rs: timestep (rs)
        elapsed_t_rs: elapsed simulation time (rs)
        v_mode: non-linear potential selector (0 linear, 1 cubic, 2 saturating, 3 double-well, 4 deficit, 5 deficit+wall)
        v_c1, v_c2: potential coefficients (see V_psi/dV_psi)
        r_wall: EMC wall radius (grid units), for modes 4/5
        wall_height: EMC wall peak density relative to background, for mode 5
        deficit_depth: core density reduction, for modes 4/5
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    c2dt2 = (c_amrs * dt_rs) ** 2
    dt2 = dt_rs * dt_rs  # scales the non-linear restoring term −dt²·dV(ψ)

    # Domain center: reference for the radial (longitudinal) zero-crossing scalar
    cx = ti.cast(nx, ti.f32) * 0.5
    cy = ti.cast(ny, ti.f32) * 0.5
    cz = ti.cast(nz, ti.f32) * 0.5

    # Update interior voxels (Dirichlet BC: ψ=0 at edges; 6-pt stencil needs a 1-cell buffer)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        laplacian = compute_laplacian(wave_field, i, j, k)
        psi_cur = wave_field.psi_am[i, j, k]

        # Leap-Frog update: ψ_new = 2ψ - ψ_prev + (c·dt)²·∇²ψ − dt²·dV(ψ)
        # For v_mode 4/5 we pass the spatial indices and profile parameters
        psi_new = (
            2.0 * psi_cur
            - wave_field.psi_prev_am[i, j, k]
            + c2dt2 * laplacian
            - dt2
            * dV_psi(
                psi_cur,
                v_mode,
                v_c1,
                v_c2,
                trackers.energy_local_aJ[i, j, k],
                i,
                j,
                k,
                nx,
                ny,
                nz,
                r_wall,
                wall_height,
                deficit_depth,
                r_soliton,
                sigma,
            )
        )
        wave_field.psi_new_am[i, j, k] = psi_new

        # WAVE-TRACKERS ============================================
        # RMS AMPLITUDE via EMA on |ψ|², with M4's unconditional decay (clears stale trails)
        disp2 = psi_new.norm() ** 2
        current_rms2 = trackers.amp_local_emarms_am[i, j, k] ** 2
        alpha_rms = 0.005  # EMA smoothing factor for RMS tracking
        new_amp = ti.sqrt(alpha_rms * disp2 + (1.0 - alpha_rms) * current_rms2)
        decay_factor = ti.cast(0.99, ti.f32)  # ~100 frames to ~37%
        trackers.amp_local_emarms_am[i, j, k] = new_amp * decay_factor

        # FREQUENCY via zero-crossing of the radial (longitudinal) projection s = ψ·r̂
        # Reference direction = from domain center (a signed scalar that oscillates everywhere)
        rdir = ti.Vector(
            [ti.cast(i, ti.f32) - cx, ti.cast(j, ti.f32) - cy, ti.cast(k, ti.f32) - cz]
        )
        r_hat = rdir / (rdir.norm() + 1e-6)
        s_prev = wave_field.psi_am[i, j, k].dot(r_hat)
        s_curr = psi_new.dot(r_hat)
        if s_prev < 0.0 and s_curr >= 0.0:  # positive-going zero crossing
            period_rs = elapsed_t_rs - trackers.last_crossing[i, j, k]
            if period_rs > dt_rs * 2.0:  # reject spurious crossings
                measured_freq = 1.0 / period_rs  # in rHz
                current_freq = trackers.freq_local_cross_rHz[i, j, k]
                alpha_freq = 0.05  # EMA smoothing factor for frequency
                trackers.freq_local_cross_rHz[i, j, k] = (
                    alpha_freq * measured_freq + (1.0 - alpha_freq) * current_freq
                )
            trackers.last_crossing[i, j, k] = elapsed_t_rs

        # LOCAL ENERGY PER VOXEL: E = ρ · V · (f · A)² in aJ (M4 formula kept; f = base frequency)
        # F = -∇E (force = negative energy gradient, computed in force_motion)
        amp_am = trackers.amp_local_emarms_am[i, j, k]
        trackers.energy_local_aJ[i, j, k] = (
            rho_qgam
            * wave_field.dx_am**3
            * (base_frequency_rHz * amp_am) ** 2
            * constants.INTERNAL_ENERGY_TO_AJ
        )

    # Swap time levels: prev ← current, current ← new
    for i, j, k in ti.ndrange(nx, ny, nz):
        wave_field.psi_prev_am[i, j, k] = wave_field.psi_am[i, j, k]
        wave_field.psi_am[i, j, k] = wave_field.psi_new_am[i, j, k]


# ================================================================
# WAVE CENTER INTERACTION  (re-drive WCs on top of the base wave)
# ================================================================
# Called each step AFTER propagate_wave (so psi_am is the new current field). Each
# mode re-imposes the wave centers as localized RADIAL sources on top of the base
# wave; they are selectable so we can A/B which best captures EWT particle behavior.
# A 4th "free" mode = call none (seed once, evolve freely). The drive region is a
# ball of `radius` voxels around each WC; the WC center voxel is a node (r̂ undefined),
# and the outer Dirichlet boundary is left untouched (ψ=0 there).
#
#   dirichlet : hard-pin ψ = A·sin(ωt+off)·r̂           (reflective driven antenna)
#   neumann   : inject velocity (ψ-ψ_prev)/dt = V·r̂     (soft drive, amplitude floats)
#   soft      : additive ψ += A·e^{-r²/2σ²}·cos(ωt+off)·r̂  (continuous radiator, back-reaction)
# ================================================================


@ti.kernel
def interact_wc_dirichlet(
    wave_field: ti.template(),  # type: ignore
    wave_center: ti.template(),  # type: ignore
    elapsed_t_rs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    boost: ti.f32,  # type: ignore
    radius: ti.i32,  # type: ignore
):
    """Hard-pin a driven radial oscillation in a ball around each active WC.

    Overwrites ψ and ψ_prev so the leapfrog sees a consistent driven antenna
    ψ = A·sin(ωt + offset)·r̂ (Dirichlet source; the base wave reflects off the pin).
    """
    omega = 2.0 * ti.math.pi * base_frequency_rHz / wave_field.scale_factor  # rad/rs
    amp = base_amplitude_am * wave_field.scale_factor * boost
    rf = ti.cast(radius, ti.f32)
    for wc in range(wave_center.num_sources):
        if wave_center.active[wc] == 0:
            continue
        c = wave_center.position_grid[wc]
        off = wave_center.offset[wc]
        s_now = ti.sin(omega * elapsed_t_rs + off)
        s_prev = ti.sin(omega * (elapsed_t_rs - dt_rs) + off)
        for di, dj, dk in ti.ndrange(
            (-radius, radius + 1), (-radius, radius + 1), (-radius, radius + 1)
        ):
            i, j, k = c[0] + di, c[1] + dj, c[2] + dk
            inside = (
                (0 < i < wave_field.nx - 1)
                and (0 < j < wave_field.ny - 1)
                and (0 < k < wave_field.nz - 1)
            )
            d = ti.Vector([ti.cast(di, ti.f32), ti.cast(dj, ti.f32), ti.cast(dk, ti.f32)])
            r = d.norm()
            if inside and r <= rf + 1e-3:
                r_hat = d / (r + 1e-6)
                node = ti.select(r < 0.5, 0.0, 1.0)  # node at the WC center voxel
                wave_field.psi_am[i, j, k] = amp * s_now * r_hat * node
                wave_field.psi_prev_am[i, j, k] = amp * s_prev * r_hat * node


@ti.kernel
def interact_wc_neumann(
    wave_field: ti.template(),  # type: ignore
    wave_center: ti.template(),  # type: ignore
    elapsed_t_rs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    boost: ti.f32,  # type: ignore
    radius: ti.i32,  # type: ignore
):
    """Driven OUTGOING (radiating) source in a ball around each active WC.

    Pins a radial traveling wave ψ = A·sin(ωt + offset − k·r)·r̂ (ψ_prev one step
    retarded). Unlike the standing Dirichlet pin, the −k·r spatial phase sends energy
    outward, so this is a flux source (net outward radiation) rather than a reflective
    antenna. It is value-bounded (|ψ| ≤ A) and therefore stable.

    Note: a strict velocity-Neumann pin (fixing (ψ-ψ_prev)/dt) is an unbounded
    integrator at the CFL dt (the drive period is under-resolved) and blows up; this
    bounded radiating pin is the stable flux-source stand-in. A true velocity-Neumann
    would need the sub-CFL dt explored in P5.
    """
    omega = 2.0 * ti.math.pi * base_frequency_rHz / wave_field.scale_factor  # rad/rs
    wavelength_grid = base_wavelength * wave_field.scale_factor / wave_field.dx
    k_grid = 2.0 * ti.math.pi / wavelength_grid  # radians per grid index
    amp = base_amplitude_am * wave_field.scale_factor * boost
    rf = ti.cast(radius, ti.f32)
    for wc in range(wave_center.num_sources):
        if wave_center.active[wc] == 0:
            continue
        c = wave_center.position_grid[wc]
        off = wave_center.offset[wc]
        for di, dj, dk in ti.ndrange(
            (-radius, radius + 1), (-radius, radius + 1), (-radius, radius + 1)
        ):
            i, j, k = c[0] + di, c[1] + dj, c[2] + dk
            inside = (
                (0 < i < wave_field.nx - 1)
                and (0 < j < wave_field.ny - 1)
                and (0 < k < wave_field.nz - 1)
            )
            d = ti.Vector([ti.cast(di, ti.f32), ti.cast(dj, ti.f32), ti.cast(dk, ti.f32)])
            r = d.norm()
            if inside and r <= rf + 1e-3:
                r_hat = d / (r + 1e-6)
                node = ti.select(r < 0.5, 0.0, 1.0)
                phase_now = omega * elapsed_t_rs + off - k_grid * r
                phase_prev = omega * (elapsed_t_rs - dt_rs) + off - k_grid * r
                wave_field.psi_am[i, j, k] = amp * ti.sin(phase_now) * r_hat * node
                wave_field.psi_prev_am[i, j, k] = amp * ti.sin(phase_prev) * r_hat * node


@ti.kernel
def interact_wc_soft(
    wave_field: ti.template(),  # type: ignore
    wave_center: ti.template(),  # type: ignore
    elapsed_t_rs: ti.f32,  # type: ignore
    boost: ti.f32,  # type: ignore
    sigma: ti.f32,  # type: ignore
    radius: ti.i32,  # type: ignore
):
    """Additive Gaussian radial forcing in a ball around each active WC.

    ψ += A·exp(-r²/2σ²)·cos(ωt + offset)·r̂. Does not overwrite ψ, so the field
    back-reacts freely (a continuous radiator). Use a small boost: additive forcing
    pumps energy every step and can drive a resonance if too strong.
    """
    omega = 2.0 * ti.math.pi * base_frequency_rHz / wave_field.scale_factor  # rad/rs
    amp = base_amplitude_am * wave_field.scale_factor * boost
    rf = ti.cast(radius, ti.f32)
    inv2s2 = 1.0 / (2.0 * sigma * sigma)
    for wc in range(wave_center.num_sources):
        if wave_center.active[wc] == 0:
            continue
        c = wave_center.position_grid[wc]
        off = wave_center.offset[wc]
        c_now = ti.cos(omega * elapsed_t_rs + off)
        for di, dj, dk in ti.ndrange(
            (-radius, radius + 1), (-radius, radius + 1), (-radius, radius + 1)
        ):
            i, j, k = c[0] + di, c[1] + dj, c[2] + dk
            inside = (
                (0 < i < wave_field.nx - 1)
                and (0 < j < wave_field.ny - 1)
                and (0 < k < wave_field.nz - 1)
            )
            d = ti.Vector([ti.cast(di, ti.f32), ti.cast(dj, ti.f32), ti.cast(dk, ti.f32)])
            r = d.norm()
            if inside and r <= rf + 1e-3:
                r_hat = d / (r + 1e-6)
                node = ti.select(r < 0.5, 0.0, 1.0)
                env = ti.exp(-(r * r) * inv2s2)
                wave_field.psi_am[i, j, k] += amp * env * c_now * r_hat * node


# ================================================================
# 3-PLANE SAMPLING FOR AVERAGE TRACKERS
# ================================================================
# PERFORMANCE NOTE: Full GPU reduction (atomic_add over all voxels) causes
# severe performance issues due to atomic contention with millions of voxels.
# The 3-plane sampling is a deliberate compromise:
# - Samples ~3N² voxels instead of N³ (e.g., 3% for 100³ grid)
# - Assumes isotropic field distribution (valid for most wave scenarios)
# - Acceptable accuracy vs massive performance gain
# ================================================================

# Cached slice buffers (initialized on first call)
_slice_xy_amp = None
_slice_xy_freq = None
_slice_xy_energy = None
_slice_xz_amp = None
_slice_xz_freq = None
_slice_xz_energy = None
_slice_yz_amp = None
_slice_yz_freq = None
_slice_yz_energy = None


@ti.kernel
def _copy_slice_xy(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    slice_energy: ti.template(),  # type: ignore
    mid_z: ti.i32,  # type: ignore
):
    """Copy center XY slice (fixed z) to 2D buffer."""
    for i, j in slice_amp:
        slice_amp[i, j] = trackers.amp_local_emarms_am[i, j, mid_z]
        slice_freq[i, j] = trackers.freq_local_cross_rHz[i, j, mid_z]
        slice_energy[i, j] = trackers.energy_local_aJ[i, j, mid_z]


@ti.kernel
def _copy_slice_xz(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    slice_energy: ti.template(),  # type: ignore
    mid_y: ti.i32,  # type: ignore
):
    """Copy center XZ slice (fixed y) to 2D buffer."""
    for i, k in slice_amp:
        slice_amp[i, k] = trackers.amp_local_emarms_am[i, mid_y, k]
        slice_freq[i, k] = trackers.freq_local_cross_rHz[i, mid_y, k]
        slice_energy[i, k] = trackers.energy_local_aJ[i, mid_y, k]


@ti.kernel
def _copy_slice_yz(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    slice_energy: ti.template(),  # type: ignore
    mid_x: ti.i32,  # type: ignore
):
    """Copy center YZ slice (fixed x) to 2D buffer."""
    for j, k in slice_amp:
        slice_amp[j, k] = trackers.amp_local_emarms_am[mid_x, j, k]
        slice_freq[j, k] = trackers.freq_local_cross_rHz[mid_x, j, k]
        slice_energy[j, k] = trackers.energy_local_aJ[mid_x, j, k]


def sample_avg_trackers(
    wave_field,
    trackers,
):
    """
    Estimate RMS amplitude and average frequency by sampling 3 orthogonal planes.

    Samples XY, XZ, and YZ center slices to avoid full 3D reduction.
    This is a deliberate performance compromise - full GPU reduction with
    atomic operations causes severe contention with millions of voxels.

    For isotropic fields, center-plane sampling provides representative estimates.

    Args:
        wave_field: WaveField instance containing grid dimensions
        trackers: WaveTrackers instance with per-voxel and average fields
    """
    global _slice_xy_amp, _slice_xy_freq, _slice_xy_energy
    global _slice_xz_amp, _slice_xz_freq, _slice_xz_energy
    global _slice_yz_amp, _slice_yz_freq, _slice_yz_energy

    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz

    # Initialize slice buffers once
    if _slice_xy_amp is None:
        _slice_xy_amp = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xy_freq = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xy_energy = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xz_amp = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_xz_freq = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_xz_energy = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_yz_amp = ti.field(dtype=ti.f32, shape=(ny, nz))
        _slice_yz_freq = ti.field(dtype=ti.f32, shape=(ny, nz))
        _slice_yz_energy = ti.field(dtype=ti.f32, shape=(ny, nz))

    # Copy 3 center slices to 2D buffers (parallel kernels)
    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2
    _copy_slice_xy(trackers, _slice_xy_amp, _slice_xy_freq, _slice_xy_energy, mid_z)
    _copy_slice_xz(trackers, _slice_xz_amp, _slice_xz_freq, _slice_xz_energy, mid_y)
    _copy_slice_yz(trackers, _slice_yz_amp, _slice_yz_freq, _slice_yz_energy, mid_x)

    # Transfer 2D slices to CPU for numpy operations
    # Exclude boundary voxels
    xy_amp = _slice_xy_amp.to_numpy()[1:-1, 1:-1]
    xy_freq = _slice_xy_freq.to_numpy()[1:-1, 1:-1]
    xy_energy = _slice_xy_energy.to_numpy()[1:-1, 1:-1]
    xz_amp = _slice_xz_amp.to_numpy()[1:-1, 1:-1]
    xz_freq = _slice_xz_freq.to_numpy()[1:-1, 1:-1]
    xz_energy = _slice_xz_energy.to_numpy()[1:-1, 1:-1]
    yz_amp = _slice_yz_amp.to_numpy()[1:-1, 1:-1]
    yz_freq = _slice_yz_freq.to_numpy()[1:-1, 1:-1]
    yz_energy = _slice_yz_energy.to_numpy()[1:-1, 1:-1]

    # Compute RMS amplitude: √(⟨A²⟩) for correct energy weighting
    # amp_local_emarms_am contains per-voxel RMS values, square them for energy
    total_amp_squared = (xy_amp**2).sum() + (xz_amp**2).sum() + (yz_amp**2).sum()
    total_freq = xy_freq.sum() + xz_freq.sum() + yz_freq.sum()
    total_energy = xy_energy.sum() + xz_energy.sum() + yz_energy.sum()
    n_samples = xy_amp.size + xz_amp.size + yz_amp.size

    trackers.amp_global_emarms_am[None] = float(np.sqrt(total_amp_squared / n_samples))
    trackers.freq_global_avg_rHz[None] = float(total_freq / n_samples)
    trackers.energy_global_avg_aJ[None] = float(total_energy / n_samples)


# ================================================================
# FLUX MESH VALUES UPDATING
# ================================================================


@ti.kernel
def update_flux_mesh_values(
    wave_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    wave_menu: ti.i32,  # type: ignore
    warp_mesh: ti.i32,  # type: ignore
):
    """
    Update flux mesh colors and vertices by sampling wave properties from voxel grid.

    Samples wave displacement at each plane vertex position and maps it to a color.
    Should be called every frame after wave propagation to update visualization.

    Args:
        wave_field: WaveField instance containing flux mesh fields and displacement data
        trackers: WaveTrackers instance with amplitude/frequency data for color scaling
        wave_menu: Selected Wave displayed with color palette
    """

    # ================================================================
    # XY Plane: Sample at z = fm_plane_z_idx
    # ================================================================
    # Always update all planes (conditionals cause GPU branch divergence)
    for i, j in ti.ndrange(wave_field.nx, wave_field.ny):
        # Sample longitudinal displacement at this voxel
        disp_value = wave_field.psi_am[i, j, wave_field.fm_plane_z_idx].norm()
        amp_value = trackers.amp_local_emarms_am[i, j, wave_field.fm_plane_z_idx]
        freq_value = trackers.freq_local_cross_rHz[i, j, wave_field.fm_plane_z_idx]
        energy_value = trackers.energy_local_aJ[i, j, wave_field.fm_plane_z_idx]
        univ_edge_z = wave_field.universe_size_am[2]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 1:  # orange
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_orange_color(
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                disp_value / univ_edge_z * warp_mesh
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
            )
        elif wave_menu == 2:  # ironbow
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                amp_value, 0, trackers.amp_global_emarms_am[None] * 2
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                amp_value / univ_edge_z * warp_mesh
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
            )
        elif wave_menu == 3:  # blueprint
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_blueprint_color(
                freq_value, 0.0, trackers.freq_global_avg_rHz[None] * 2
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = freq_value / trackers.freq_global_avg_rHz[
                None
            ] / 3000 * warp_mesh + wave_field.flux_mesh_planes[2] * (
                wave_field.nz / wave_field.max_grid_size
            )
        elif wave_menu == 4:  # ironbow
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                energy_value, 0.0, trackers.energy_global_avg_aJ[None] * 2
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                energy_value / univ_edge_z * warp_mesh / 10
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
            )

    # ================================================================
    # XZ Plane: Sample at y = fm_plane_y_idx
    # ================================================================
    for i, k in ti.ndrange(wave_field.nx, wave_field.nz):
        # Sample longitudinal displacement at this voxel
        disp_value = wave_field.psi_am[i, wave_field.fm_plane_y_idx, k].norm()
        amp_value = trackers.amp_local_emarms_am[i, wave_field.fm_plane_y_idx, k]
        freq_value = trackers.freq_local_cross_rHz[i, wave_field.fm_plane_y_idx, k]
        energy_value = trackers.energy_local_aJ[i, wave_field.fm_plane_y_idx, k]
        univ_edge_y = wave_field.universe_size_am[1]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 1:  # orange
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_orange_color(
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                disp_value / univ_edge_y * warp_mesh
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
            )
        elif wave_menu == 2:  # ironbow
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                amp_value, 0, trackers.amp_global_emarms_am[None] * 2
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                amp_value / univ_edge_y * warp_mesh
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
            )
        elif wave_menu == 3:  # blueprint
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_blueprint_color(
                freq_value, 0.0, trackers.freq_global_avg_rHz[None] * 2
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = freq_value / trackers.freq_global_avg_rHz[
                None
            ] / 3000 * warp_mesh + wave_field.flux_mesh_planes[1] * (
                wave_field.ny / wave_field.max_grid_size
            )
        elif wave_menu == 4:  # ironbow
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                energy_value, 0.0, trackers.energy_global_avg_aJ[None] * 2
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                energy_value / univ_edge_y * warp_mesh / 10
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
            )

    # ================================================================
    # YZ Plane: Sample at x = fm_plane_x_idx
    # ================================================================
    for j, k in ti.ndrange(wave_field.ny, wave_field.nz):
        # Sample longitudinal displacement at this voxel
        disp_value = wave_field.psi_am[wave_field.fm_plane_x_idx, j, k].norm()
        amp_value = trackers.amp_local_emarms_am[wave_field.fm_plane_x_idx, j, k]
        freq_value = trackers.freq_local_cross_rHz[wave_field.fm_plane_x_idx, j, k]
        energy_value = trackers.energy_local_aJ[wave_field.fm_plane_x_idx, j, k]
        univ_edge_x = wave_field.universe_size_am[0]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 1:  # orange
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_orange_color(
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                disp_value / univ_edge_x * warp_mesh
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
            )
        elif wave_menu == 2:  # ironbow
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                amp_value, 0, trackers.amp_global_emarms_am[None] * 2
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                amp_value / univ_edge_x * warp_mesh
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
            )
        elif wave_menu == 3:  # blueprint
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_blueprint_color(
                freq_value, 0.0, trackers.freq_global_avg_rHz[None] * 2
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = freq_value / trackers.freq_global_avg_rHz[
                None
            ] / 3000 * warp_mesh + wave_field.flux_mesh_planes[0] * (
                wave_field.nx / wave_field.max_grid_size
            )
        elif wave_menu == 4:  # ironbow
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                energy_value, 0.0, trackers.energy_global_avg_aJ[None] * 2
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                energy_value / univ_edge_x * warp_mesh / 10
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
            )


# ================================================================
# POSITION RENDER
# ================================================================


@ti.kernel
def sample_position_to_render(
    wave_field: ti.template(),  # type: ignore
    amp_boost: ti.f32,  # type: ignore
    stride: ti.i32,  # type: ignore
    num_render: ti.i32,  # type: ignore
):
    """Sample granule positions from z flux_mesh plane with stride, writing to 1D position_render.

    Samples every `stride`-th voxel from the XY plane at the z flux mesh index,
    capping output at `num_render` particles for performance.
    """
    k = int(wave_field.flux_mesh_planes[2] * wave_field.grid_size[2])
    max_dim = ti.cast(wave_field.max_grid_size, ti.f32)
    sampled_ny = (wave_field.ny + stride - 1) // stride  # cols per sampled row

    for render_idx in range(num_render):
        si = render_idx // sampled_ny  # sampled row index
        sj = render_idx % sampled_ny  # sampled col index
        i = si * stride
        j = sj * stride
        displaced = amp_boost * wave_field.psi_am[i, j, k] / wave_field.dx_am + ti.Vector(
            [ti.cast(i, ti.f32), ti.cast(j, ti.f32), ti.cast(k, ti.f32)]
        )
        wave_field.position_render[render_idx] = displaced / max_dim
