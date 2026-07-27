"""
FORCE & MOTION MODULE

Implements force calculation from energy gradients and particle motion integration.

Physics Foundation:
- Energy per voxel: E = ρ · V · (f · A)²  (EWT energy equation)
- Force: F = -∇E  (negative gradient of energy density)
- Motion: Euler integration of F = m · a

Units:
- Energy: aJ (attojoules, 1 aJ = 1e-18 J)
- Force: Newtons (1 aJ/am = 1 N, no conversion needed)
- Mass: qg (quectograms, 1 qg = 1e-33 kg, for f32 precision on GPU)
- Velocity: am/rs (OpenWave scaled units for f32 precision)
- Position: grid indices (float)
- Time: rontoseconds (rs)

Conversion factors:
- 1 aJ / 1 am = 1e-18 J / 1e-18 m = 1 N  (energy gradient → force)
- a_amrs2 = (F_N / m_qg) * 1e-3            (N/qg → am/rs²)
- c = 0.3 am/rs                             (speed of light)

See research/02_force_motion.md for detailed documentation.
"""

import taichi as ti

import numpy as np

from openwave.common import constants

# ================================================================
# Physical Constants (cached for kernel access)
# ================================================================
ATTOMETER = constants.ATTOMETER  # m/am = 1e-18
RONTOSECOND = constants.RONTOSECOND  # s/rs = 1e-27
QUECTOGRAM = constants.QUECTOGRAM  # kg/qg = 1e-33

# Coulomb force constants (for reference comparisons)
COULOMB_CONSTANT = constants.COULOMB_CONSTANT  # N·m²/C², k = 8.99e9
ELEMENTARY_CHARGE = constants.ELEMENTARY_CHARGE  # C, e = 1.60e-19

# EWT particle constants (for EWT force reference)
MEDIUM_DENSITY = constants.MEDIUM_DENSITY  # kg/m³
EWAVE_AMPLITUDE = constants.EWAVE_AMPLITUDE  # m
WAVE_SPEED = constants.WAVE_SPEED  # m/s
EWAVE_LENGTH = constants.EWAVE_LENGTH  # m
ELECTRON_K = constants.ELECTRON_K
ELECTRON_OUTER_SHELL = constants.ELECTRON_OUTER_SHELL
ELECTRON_ORBITAL_G = constants.ELECTRON_ORBITAL_G

# ================================================================
# DEFAULT VELOCITY DAMPING 
# ================================================================
VELOCITY_DAMPING_DEFAULT = 0.999


# ================================================================
# PRESSURE FORCE FROM VACUUM (V_MODE >= 4)
# ================================================================

@ti.func
def compute_density_gradient(
    wave_field: ti.template(),
    pos: ti.types.vector(3, ti.f32),
    r_soliton: ti.f32,
    sigma: ti.f32,
    deficit_depth: ti.f32,
    r_wall: ti.f32,
    wall_height: ti.f32,
    v_mode: ti.i32,
) -> ti.types.vector(3, ti.f32):
    """
    Compute the density gradient ∇ρ(r) at position pos.
    Uses analytic derivative of the density profile.
    
    For V_MODE=4/5: exponential profile
    For V_MODE=6/7: sigmoid profile
    For V_MODE=9/10: Gaussian profile
    For V_MODE=0/1: returns zero (uniform density)
    """
    cx = wave_field.nx * 0.5
    cy = wave_field.ny * 0.5
    cz = wave_field.nz * 0.5
    dr = ti.Vector([pos[0] - cx, pos[1] - cy, pos[2] - cz])
    r = dr.norm() + 1e-6
    r_hat = dr / r
    
    drho_dr = 0.0
    
    if v_mode == 4 or v_mode == 5:
        # Exponential profile: rho = 1 - deficit * exp(-r / (r_wall * 0.3))
        scale = r_wall * 0.3
        drho_dr = deficit_depth * (1.0 / scale) * ti.exp(-r / scale)
        
    elif v_mode == 6 or v_mode == 7:
        # Sigmoid profile: rho = 1 - deficit * 0.5 * (1 - tanh((r - R)/sigma))
        tanh_arg = (r - r_soliton) / sigma
        tanh_val = ti.tanh(tanh_arg)
        sech2 = 1.0 - tanh_val * tanh_val
        drho_dr = deficit_depth * 0.5 * sech2 / sigma
        
        # For V_MODE=7, add wall contribution (Gaussian bump)
        if v_mode == 7:
            wall_sigma = r_wall * 0.1
            d_wall_dr = -(wall_height - 1.0) * (r - r_wall) / (wall_sigma ** 2) * ti.exp(-((r - r_wall) ** 2) / (2.0 * wall_sigma ** 2))
            drho_dr += d_wall_dr
    
    elif v_mode == 9 or v_mode == 10:
        # Gaussian profile: rho = 1 - deficit * exp(-(r/R)^2)
        # derivative: drho_dr = 2 * deficit * r / R^2 * exp(-(r/R)^2)
        drho_dr = 2.0 * deficit_depth * r / (r_soliton ** 2) * ti.exp(-(r / r_soliton) ** 2)
    
    # Gradient: ∇ρ = (dρ/dr) * r_hat
    grad_rho = drho_dr * r_hat
    return grad_rho


@ti.kernel
def add_pressure_force(
    wave_center: ti.template(),
    wave_field: ti.template(),
    pressure_strength: ti.f32,
    v_mode: ti.i32,
    r_soliton: ti.f32,
    sigma: ti.f32,
    deficit_depth: ti.f32,
    r_wall: ti.f32,
    wall_height: ti.f32,
):
    """
    Add vacuum pressure force to each active WC.
    F_pressure = -pressure_strength * ∇ρ(r)
    
    Since ∇ρ points outward, -∇ρ points inward → binding force.
    """
    for wc_idx in range(wave_center.num_sources):
        if wave_center.active[wc_idx] == 0:
            continue
        
        # Skip if density is uniform (V_MODE=0 or 1)
        if v_mode == 0 or v_mode == 1:
            continue
        
        pos = wave_center.position_float[wc_idx]
        grad_rho = compute_density_gradient(
            wave_field, pos, r_soliton, sigma, deficit_depth, r_wall, wall_height, v_mode
        )
        
        # Pressure force: F = -strength * grad_rho
        # Points inward (toward center)
        F_pressure = -pressure_strength * grad_rho
        
        # Add to existing force (from energy gradient)
        wave_center.force[wc_idx][0] += F_pressure[0]
        wave_center.force[wc_idx][1] += F_pressure[1]
        wave_center.force[wc_idx][2] += F_pressure[2]


def compute_ewt_electric_force(
    r: float, K: int = 1, Oe: float = 1.0, glambda: float = 1.0
) -> float:
    """
    Compute EWT electric force between two particles (reference/validation).

    F_e = (4πρ K^7 A^6 c² Oe / 3λ²) × gλ × (Q1×Q2 / r²)

    Args:
        r: Distance between particles in meters
        K: Wave center count (1 for neutrino, 10 for electron)
        Oe: Outer shell multiplier (1.0 for K=1, ~2.14 for electron)
        glambda: Orbital g-factor (1.0 for K=1, ~0.99 for electron)

    Returns:
        Force in Newtons
    """
    coefficient = (
        4.0
        * np.pi
        * MEDIUM_DENSITY
        * (K**7)
        * (EWAVE_AMPLITUDE**6)
        * (WAVE_SPEED**2)
        * Oe
        * glambda
    ) / (3.0 * (EWAVE_LENGTH**2))

    return coefficient / (r**2)


# ================================================================
# Force from Energy Gradient: F = -∇E
# ================================================================

# Gradient sampling: number of voxel shells and weight exponent
# GRADIENT_SAMPLE_RADIUS = 1: standard central difference (single shell)
# GRADIENT_SAMPLE_RADIUS > 1: weighted multi-shell gradient (better well resolution)
# GRADIENT_WEIGHT_FALLOFF = 2: weights as 1/d² (particle energy density ∝ A² ∝ 1/r²)
GRADIENT_SAMPLE_RADIUS = 3  # voxels (increased from 1 for better lock-in well resolution)
GRADIENT_WEIGHT_FALLOFF = 2  # exponent for 1/d^n weighting


@ti.kernel
def compute_force_vector(
    wave_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    wave_center: ti.template(),  # type: ignore
):
    """
    Compute force on each wave center from energy gradient.

    F = -∇E where E is the local energy field (energy_local_aJ).
    Uses weighted central differences: multiple sample shells within
    GRADIENT_SAMPLE_RADIUS, weighted by 1/d² (particle energy density
    falloff ∝ A² ∝ 1/r²). Closer voxels contribute more, matching how
    a particle's own wave structure "feels" the surrounding energy landscape.
    For R=1: identical to standard central difference (single shell).
    Units: aJ/am = N (no conversion needed).

    ┌───────┬────────┬────────────┐
    │ Shell │ Weight │ Percentage │
    ├───────┼────────┼────────────┤
    │ d=1   │ 1.0    │ 73.5%      │
    ├───────┼────────┼────────────┤
    │ d=2   │ 0.25   │ 18.4%      │
    ├───────┼────────┼────────────┤
    │ d=3   │ 0.111  │ 8.2%       │
    └───────┴────────┴────────────┘

    Scale correction: amplitude boost from scale_factor makes energy S² too
    large, so F_physical = F_computed / S².

    Args:
        wave_field: WaveField instance containing grid info
        trackers: Trackers instance with energy_local_aJ field
        wave_center: WaveCenter instance to store computed forces
    """
    dx_am = wave_field.dx_am

    # Scale factor correction: energy scales as S² from amplitude boost
    S = wave_field.scale_factor
    S2 = S * S

    # Precompute weights: 1/d^GRADIENT_WEIGHT_FALLOFF for each shell d = 1..R
    # Physical basis: particle energy density ∝ |ψ|² ∝ 1/r²
    w_sum = ti.cast(0.0, ti.f32)
    for d in ti.static(range(1, GRADIENT_SAMPLE_RADIUS + 1)):
        w_sum += 1.0 / ti.cast(d**GRADIENT_WEIGHT_FALLOFF, ti.f32)

    for wc_idx in range(wave_center.num_sources):
        # Skip inactive (annihilated) WCs
        if wave_center.active[wc_idx] == 0:
            continue

        # Get wave center grid position
        i = wave_center.position_grid[wc_idx][0]
        j = wave_center.position_grid[wc_idx][1]
        k = wave_center.position_grid[wc_idx][2]

        # Initialize force to zero
        F_x = ti.cast(0.0, ti.f32)
        F_y = ti.cast(0.0, ti.f32)
        F_z = ti.cast(0.0, ti.f32)

        # Grid dimensions
        nx = wave_field.nx
        ny = wave_field.ny
        nz = wave_field.nz

        # Boundary check
        if (
            i > GRADIENT_SAMPLE_RADIUS
            and i < nx - GRADIENT_SAMPLE_RADIUS
            and j > GRADIENT_SAMPLE_RADIUS
            and j < ny - GRADIENT_SAMPLE_RADIUS
            and k > GRADIENT_SAMPLE_RADIUS
            and k < nz - GRADIENT_SAMPLE_RADIUS
        ):
            # Weighted gradient: Σ w(d) · (E[+d] - E[-d]) / (2·d·dx) / Σ w(d)
            # w(d) = 1/d^GRADIENT_WEIGHT_FALLOFF (energy density weighting)
            grad_x = ti.cast(0.0, ti.f32)
            grad_y = ti.cast(0.0, ti.f32)
            grad_z = ti.cast(0.0, ti.f32)

            for d in ti.static(range(1, GRADIENT_SAMPLE_RADIUS + 1)):
                w = 1.0 / ti.cast(d**GRADIENT_WEIGHT_FALLOFF, ti.f32)
                dist = 2.0 * d * dx_am

                grad_x += (
                    w
                    * (
                        trackers.energy_local_aJ[i + d, j, k]
                        - trackers.energy_local_aJ[i - d, j, k]
                    )
                    / dist
                )
                grad_y += (
                    w
                    * (
                        trackers.energy_local_aJ[i, j + d, k]
                        - trackers.energy_local_aJ[i, j - d, k]
                    )
                    / dist
                )
                grad_z += (
                    w
                    * (
                        trackers.energy_local_aJ[i, j, k + d]
                        - trackers.energy_local_aJ[i, j, k - d]
                    )
                    / dist
                )

            # F = -∇E / S² (aJ/am = N, with scale correction)
            F_x = -grad_x / (w_sum * S2)
            F_y = -grad_y / (w_sum * S2)
            F_z = -grad_z / (w_sum * S2)

        wave_center.force[wc_idx][0] = F_x
        wave_center.force[wc_idx][1] = F_y
        wave_center.force[wc_idx][2] = F_z


# ================================================================
# Motion Integration (Euler)
# ================================================================


@ti.kernel
def integrate_motion_euler(
    wave_field: ti.template(),  # type: ignore
    wave_center: ti.template(),  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    damping: ti.f32,  # type: ignore  [A2] now passed as parameter
):
    """
    Integrate particle motion using Euler method.

    v_new = v_old + a * dt  (velocity in am/rs)
    x_new = x_old + v_new * dt  (position in grid indices)

    Args:
        wave_field: WaveField instance (for dx voxel size)
        wave_center: WaveCenter instance with force/velocity/position fields
        dt_rs: Timestep in rontoseconds
        damping: per-experiment velocity damping factor (1.0 = no damping)
    """
    # Conversion factor: (N / qg) to am/rs²
    accel_conv_qg = ti.cast(1e-3, ti.f32)  # (F_N / m_qg) * 1e-3 -> am/rs²

    # Voxel size in attometers for position conversion
    dx_am = wave_field.dx / ti.cast(ATTOMETER, ti.f32)

    for wc_idx in range(wave_center.num_sources):
    # Skip inactive (annihilated) WCs
        if wave_center.active[wc_idx] == 0:
            continue

             # Get force (Newtons) and mass (qg - quectograms for GPU precision)
        F_x = wave_center.force[wc_idx][0]
        F_y = wave_center.force[wc_idx][1]
        F_z = wave_center.force[wc_idx][2]
        m_qg = wave_center.mass_qg[wc_idx]

        # Acceleration: a = F/m, then convert (N/qg) to am/rs²
        a_x_amrs = (F_x / m_qg) * accel_conv_qg
        a_y_amrs = (F_y / m_qg) * accel_conv_qg
        a_z_amrs = (F_z / m_qg) * accel_conv_qg

        # Update velocity: v_new = v_old + a * dt (in am/rs)
        wave_center.velocity_amrs[wc_idx][0] += a_x_amrs * dt_rs
        wave_center.velocity_amrs[wc_idx][1] += a_y_amrs * dt_rs
        wave_center.velocity_amrs[wc_idx][2] += a_z_amrs * dt_rs

        # Apply velocity damping (models radiation energy loss)
        wave_center.velocity_amrs[wc_idx][0] *= damping
        wave_center.velocity_amrs[wc_idx][1] *= damping
        wave_center.velocity_amrs[wc_idx][2] *= damping

        # Clamp velocity to speed of light (c = 0.3 am/rs)
        c_amrs = ti.cast(0.3, ti.f32)
        v_mag = ti.sqrt(
            wave_center.velocity_amrs[wc_idx][0] ** 2
            + wave_center.velocity_amrs[wc_idx][1] ** 2
            + wave_center.velocity_amrs[wc_idx][2] ** 2
        )
        if v_mag > c_amrs:
            scale = c_amrs / v_mag
            wave_center.velocity_amrs[wc_idx][0] *= scale
            wave_center.velocity_amrs[wc_idx][1] *= scale
            wave_center.velocity_amrs[wc_idx][2] *= scale

        # Position change in attometers
        dx_am_step = wave_center.velocity_amrs[wc_idx][0] * dt_rs
        dy_am_step = wave_center.velocity_amrs[wc_idx][1] * dt_rs
        dz_am_step = wave_center.velocity_amrs[wc_idx][2] * dt_rs

        # Convert to grid index change
        di = dx_am_step / dx_am
        dj = dy_am_step / dx_am
        dk = dz_am_step / dx_am

        wave_center.position_float[wc_idx][0] += di
        wave_center.position_float[wc_idx][1] += dj
        wave_center.position_float[wc_idx][2] += dk

        # Sync integer position (round instead of truncate)
        wave_center.position_grid[wc_idx][0] = ti.cast(
            ti.round(wave_center.position_float[wc_idx][0]), ti.i32
        )
        wave_center.position_grid[wc_idx][1] = ti.cast(
            ti.round(wave_center.position_float[wc_idx][1]), ti.i32
        )
        wave_center.position_grid[wc_idx][2] = ti.cast(
            ti.round(wave_center.position_float[wc_idx][2]), ti.i32
        )


# ================================================================
# Motion Integration (Leapfrog / Velocity Verlet)
# ================================================================


@ti.kernel
def integrate_motion_leapfrog(
    wave_field: ti.template(),  # type: ignore
    wave_center: ti.template(),  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    damping: ti.f32,  # type: ignore  
):
    """
    Integrate particle motion using Velocity Verlet (leapfrog) method.

    Unlike Euler, leapfrog is symplectic — it conserves energy in oscillatory
    systems, preventing numerical drift that causes particles to escape
    lock-in wells. The method uses half-step velocity updates:

    1. v(t + dt/2) = v(t) + a(t) * dt/2       (half-step kick)
    2. x(t + dt)   = x(t) + v(t + dt/2) * dt   (full-step drift)
    3. compute a(t + dt) from new positions      (done externally)
    4. v(t + dt)   = v(t + dt/2) + a(t + dt) * dt/2  (half-step kick)

    Steps 1+2 are done here. Step 3 is the force computation (external).
    Step 4 is done on the NEXT call (the first half-kick uses the NEW force).

    In practice, we store v at half-steps and do:
    v += a * dt   (full kick, combining two half-kicks across calls)
    x += v * dt   (drift with updated velocity)

    This is equivalent to standard leapfrog and is symplectic.

    Args:
        wave_field: WaveField instance (for dx voxel size)
        wave_center: WaveCenter instance with force/velocity/position fields
        dt_rs: Timestep in rontoseconds
        damping: per-experiment velocity damping factor (1.0 = no damping)
    """
    accel_conv_qg = ti.cast(1e-3, ti.f32)  # (F_N / m_qg) * 1e-3 -> am/rs²
    dx_am = wave_field.dx / ti.cast(ATTOMETER, ti.f32)
    damp = ti.cast(damping, ti.f32)  

    for wc_idx in range(wave_center.num_sources):
        if wave_center.active[wc_idx] == 0:
            continue

        F_x = wave_center.force[wc_idx][0]
        F_y = wave_center.force[wc_idx][1]
        F_z = wave_center.force[wc_idx][2]
        m_qg = wave_center.mass_qg[wc_idx]

        # Acceleration from current force
        a_x = (F_x / m_qg) * accel_conv_qg
        a_y = (F_y / m_qg) * accel_conv_qg
        a_z = (F_z / m_qg) * accel_conv_qg

        # Full velocity kick: v += a * dt (leapfrog: combines two half-kicks)
        wave_center.velocity_amrs[wc_idx][0] += a_x * dt_rs
        wave_center.velocity_amrs[wc_idx][1] += a_y * dt_rs
        wave_center.velocity_amrs[wc_idx][2] += a_z * dt_rs

        # Apply damping (models radiation energy loss)
        wave_center.velocity_amrs[wc_idx][0] *= damp
        wave_center.velocity_amrs[wc_idx][1] *= damp
        wave_center.velocity_amrs[wc_idx][2] *= damp

        # Clamp velocity to speed of light (c = 0.3 am/rs)
        c_amrs = ti.cast(0.3, ti.f32)
        v_mag = ti.sqrt(
            wave_center.velocity_amrs[wc_idx][0] ** 2
            + wave_center.velocity_amrs[wc_idx][1] ** 2
            + wave_center.velocity_amrs[wc_idx][2] ** 2
        )
        if v_mag > c_amrs:
            scale = c_amrs / v_mag
            wave_center.velocity_amrs[wc_idx][0] *= scale
            wave_center.velocity_amrs[wc_idx][1] *= scale
            wave_center.velocity_amrs[wc_idx][2] *= scale

        # Drift: x += v * dt (position update with kicked velocity)
        dx_am_step = wave_center.velocity_amrs[wc_idx][0] * dt_rs
        dy_am_step = wave_center.velocity_amrs[wc_idx][1] * dt_rs
        dz_am_step = wave_center.velocity_amrs[wc_idx][2] * dt_rs

        di = dx_am_step / dx_am
        dj = dy_am_step / dx_am
        dk = dz_am_step / dx_am

        wave_center.position_float[wc_idx][0] += di
        wave_center.position_float[wc_idx][1] += dj
        wave_center.position_float[wc_idx][2] += dk

        # Sync integer position
        wave_center.position_grid[wc_idx][0] = ti.cast(
            ti.round(wave_center.position_float[wc_idx][0]), ti.i32
        )
        wave_center.position_grid[wc_idx][1] = ti.cast(
            ti.round(wave_center.position_float[wc_idx][1]), ti.i32
        )
        wave_center.position_grid[wc_idx][2] = ti.cast(
            ti.round(wave_center.position_float[wc_idx][2]), ti.i32
        )


# ================================================================
# Annihilation Detection
# ================================================================


@ti.kernel
def detect_annihilation(
    wave_center: ti.template(),  # type: ignore
    annihilation_threshold: ti.f32,  # type: ignore
):
    """
    Annihilation naturally occurs from wave physics, but needs numerical precision check.
    Detect and handle particle annihilation when WCs converge to same position.

    When two wave centers with opposite phase (180°) attract and meet:
    1. Their waves cancel perfectly (handled by wave precision rounding)
    2. Snap both WCs to exact same position
    3. Zero velocities and mark inactive

    This ensures annihilation is permanent - no wave reappearance from micro-motion.
    Numerical precision limits may cause slight separation otherwise.

    Args:
        wave_center: WaveCenter instance with position/velocity fields
        annihilation_threshold: Distance in grid units to trigger annihilation
    """
    phase_tolerance = ti.cast(0.17, ti.f32)  # ~10° tolerance
    pi = ti.cast(3.14159265359, ti.f32)

    # Check all pairs (i, j) where i < j
    for i in range(wave_center.num_sources):
        for j in range(i + 1, wave_center.num_sources):
            # Skip if either WC is already inactive
            if wave_center.active[i] == 0 or wave_center.active[j] == 0:
                continue

            # Check if phases are opposite (differ by ~π)
            phase_diff = ti.abs(wave_center.offset[i] - wave_center.offset[j])
            phase_diff_normalized = ti.abs(phase_diff - pi)
            if phase_diff_normalized > phase_tolerance:
                continue

            # Calculate distance between WCs (grid units)
            dx = wave_center.position_float[i][0] - wave_center.position_float[j][0]
            dy = wave_center.position_float[i][1] - wave_center.position_float[j][1]
            dz = wave_center.position_float[i][2] - wave_center.position_float[j][2]
            distance = ti.sqrt(dx * dx + dy * dy + dz * dz)

            if distance < annihilation_threshold:
                # Snap both WCs to midpoint
                mid_x = ti.round(
                    (wave_center.position_float[i][0] + wave_center.position_float[j][0]) / 2.0
                )
                mid_y = ti.round(
                    (wave_center.position_float[i][1] + wave_center.position_float[j][1]) / 2.0
                )
                mid_z = ti.round(
                    (wave_center.position_float[i][2] + wave_center.position_float[j][2]) / 2.0
                )

                for idx in ti.static(range(2)):
                    wc = i if idx == 0 else j
                    wave_center.position_float[wc][0] = mid_x
                    wave_center.position_float[wc][1] = mid_y
                    wave_center.position_float[wc][2] = mid_z
                    wave_center.velocity_amrs[wc][0] = 0.0
                    wave_center.velocity_amrs[wc][1] = 0.0
                    wave_center.velocity_amrs[wc][2] = 0.0
                    wave_center.position_grid[wc][0] = ti.cast(mid_x, ti.i32)
                    wave_center.position_grid[wc][1] = ti.cast(mid_y, ti.i32)
                    wave_center.position_grid[wc][2] = ti.cast(mid_z, ti.i32)
                    wave_center.active[wc] = 0