"""
M5 ENGINE — FIELD OBSERVABLES (engine3_observables)

Derived measurements computed from the field each frame:
  - update_trackers_M          per-voxel ‖M−D‖_F amplitude + ‖Ṁ‖_F clock (Trackers)
  - compute_energyH_density_M  Hamiltonian energy density (FieldObservables)
  - compute_energyF_density    Frank elastic energy density
  - compute_winding_number     topological charge Q on a sphere (CPU numpy)
  - sample_avg_* + helpers     3-plane global aggregates

Depends on engine2_pde (V_M, used by compute_energyH_density_M).
"""

import taichi as ti
import numpy as np

from . import engine2_pde as pde

# ================================================================
# MATRIX-SUBSTRATE TRACKERS (M5.4) — amplitude ‖M−D‖_F + clock ‖Ṁ‖_F
# ================================================================
# The matrix-world replacement for update_trackers' ψ observables (4b
# tracker-redefinitions table). Reuses the SAME Trackers fields, EMA cadence, and
# 3-plane aggregation (sample_avg_trackers) — only the sampled quantity changes:
#   - amplitude = ‖M − D_vac‖_F  (Frobenius deviation from the ẑ-vacuum order
#       parameter D_vac = diag(δ, δ, 1)) → amplitude A (the joint (A, ω) excess)
#   - frequency = ‖Ṁ‖_F = ‖M − M_prev‖_F / dt  (order-parameter rotation rate)
#       → de Broglie clock ω (M5.8 headline). Zero for a static
#         configuration; becomes the genuine clock once M is leapfrogged (M5.6+).
# The ψ_z zero-crossing it replaces was always a convention hack (the update_trackers
# docstring admits it) — the matrix rotation rate is the real thing. This is the
# "measurement infra lands EARLY" deliverable from the M5.3 analysis.


@ti.kernel
def update_trackers_M(
    tensor_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    delta: ti.f32,  # type: ignore
):
    """Matrix-substrate amplitude/frequency trackers (see section header).

    Reads:  tensor_field.M_am, tensor_field.M_prev_am
    Writes: trackers.amp_local_emarms_am  (EMA of ‖M − D_vac‖_F — amplitude A)
            trackers.freq_local_cross_rHz (EMA of ‖Ṁ‖_F        — clock ω)

    Args:
        tensor_field: TensorField (reads M_am, M_prev_am)
        trackers: Trackers (writes the per-voxel amp/freq EMA fields)
        dt_rs: timestep (rs) — divides ‖M − M_prev‖_F to give the rate ‖Ṁ‖_F
        delta: uniaxial minor-axis eigenvalue (tensor_field.lc_delta) → D_vac diagonal
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    decay = ti.cast(0.999, ti.f32)
    alpha_rms = ti.cast(0.005, ti.f32)
    alpha_freq = ti.cast(0.05, ti.f32)
    inv_dt = 1.0 / dt_rs
    # 4×4 ẑ-vacuum D_vac = diag(g, δ, δ, 1) (Duda index-0): time axis (index 0) = g, spatial
    # uniaxial diag(δ, δ, 1) at indices 1,2,3 (director = ẑ, eigenvalue 1 at index 0). The
    # Frobenius deviation spans the time block (g−g=0 in vacuum) so the amplitude tracker stays ≈0.
    g = tensor_field.lc_g
    d_vac = ti.Matrix([[g, 0.0, 0.0, 0.0],
                       [0.0, delta, 0.0, 0.0],
                       [0.0, 0.0, delta, 0.0],
                       [0.0, 0.0, 0.0, 1.0]])

    for i, j, k in ti.ndrange(nx, ny, nz):
        m = tensor_field.M_am[i, j, k]

        # Amplitude: EMA of ‖M − D_vac‖_F² (Frobenius), then √ — same RMS form as
        # update_trackers but on the order-parameter deviation instead of |ψ|².
        dev2 = (m - d_vac).norm_sqr()
        rms2_old = trackers.amp_local_emarms_am[i, j, k] ** 2
        rms2_new = alpha_rms * dev2 + (1.0 - alpha_rms) * rms2_old
        trackers.amp_local_emarms_am[i, j, k] = ti.sqrt(rms2_new) * decay

        # Frequency: EMA of ‖Ṁ‖_F = ‖M − M_prev‖_F / dt (frame rotation rate).
        mdot_mag = (m - tensor_field.M_prev_am[i, j, k]).norm() * inv_dt
        old_freq = trackers.freq_local_cross_rHz[i, j, k]
        trackers.freq_local_cross_rHz[i, j, k] = (
            alpha_freq * mdot_mag + (1.0 - alpha_freq) * old_freq
        ) * decay


@ti.kernel
def compute_energyH_density_M(
    tensor_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    a: ti.f32,  # type: ignore
    b: ti.f32,  # type: ignore
    c: ti.f32,  # type: ignore
    v0: ti.f32,  # type: ignore
    e_scale: ti.f32,  # type: ignore
):
    """Matrix Hamiltonian density (M5.5.4) — the Eq.18 ℋ for the matrix substrate.

        H = ½‖Ṁ‖²_F  +  c²·4 Σ_{μ<ν}‖[M_μ,M_ν]‖²_F  +  (V_M(M) − v0)
            kinetic        curvature (= ¼Σ‖F_μν‖², the page-18 energy)   potential

    Replaces the dormant-ψ placeholder (uniform ¼λ) for matrix xperiments — resolves
    the M5.4 WAVE_MENU=4 carry-over. Matches the evolve_M force so energy is conserved
    (simple ½‖Ṁ‖² kinetic per the 2026-05-26 decision). `e_scale` = physical-energy
    factor (ρ_medium × voxel_volume_am³ × INTERNAL_ENERGY_TO_AJ); pass 1.0 for bare units.

    `v0` = vacuum potential V_M(D_vacuum), subtracted so the displayed energy is measured
    from the vacuum (M5.6.5c): the LdG well bottom is NEGATIVE (V_min = −c·s₂*² for the b=0
    well), so with V on the constant floor (~−1.8e-6) swamps the tiny curvature (~1e-11) and
    the field renders as a uniform "zero". Subtracting v0 zeroes the vacuum and reveals the
    defect structure (curvature rods + the V-deviation core). A constant shift does NOT touch
    the force (−dV_M, unchanged) so dynamics/conservation are identical; pass 0.0 for V off.

    Reads:  tensor_field.M_am, tensor_field.M_prev_am
    Writes: observables.energyH_density_aJ
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_dt = 1.0 / dt_rs
    c2 = c_amrs * c_amrs
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)

    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        m = tensor_field.M_am[i, j, k]
        m_dot = (m - tensor_field.M_prev_am[i, j, k]) * inv_dt
        kinetic = 0.5 * m_dot.norm_sqr()                      # ½‖Ṁ‖²_F

        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        cxy = pde.commutator(mx, my)
        cxz = pde.commutator(mx, mz)
        cyz = pde.commutator(my, mz)
        curvature = 4.0 * (cxy.norm_sqr() + cxz.norm_sqr() + cyz.norm_sqr())

        potential = pde.V_M(m, a, b, c) - v0
        observables.energyH_density_aJ[i, j, k] = e_scale * (kinetic + c2 * curvature + potential)


@ti.kernel
def compute_energyH_density_eta(
    tensor_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    dt_eff: ti.f32,  # type: ignore
    dx_eta: ti.f32,  # type: ignore
    sg: ti.f32,  # type: ignore
    delta: ti.f32,  # type: ignore
    w1: ti.f32,  # type: ignore
    e_scale: ti.f32,  # type: ignore
):
    """M5.24 — the verified-L era energy density for the CANONICAL integrator path:

        H = ½‖Ṁ‖²_F  +  ½ Σ_branch 4 Σ_{i<j} <F_ij, F_ij>_eta  +  V4(M)

    (canonical-kinetic + η-signed curvature on the symmetrized stencil + the
    universal spectral potential). Ṁ = (M − M_prev)/dt_eff, τ-units like the
    canonical evolve pair. NO v0 subtraction: V4 and the curvature are both
    EXACTLY zero on the covariant vacuum diag(−sg, 1, δ, 0) — the display
    floor is the true zero (retires the plain-path v0 hack on this path).

    Reads:  tensor_field.M_am, M_prev_am;  Writes: observables.energyH_density_aJ
    dx_eta = the canonical-unit grid spacing (must match the evolve kernels).
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_dt = 1.0 / dt_eff
    inv_dx = 1.0 / dx_eta
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        m0 = tensor_field.M_am[i, j, k]
        m_dot = (m0 - tensor_field.M_prev_am[i, j, k]) * inv_dt
        kinetic = 0.5 * m_dot.norm_sqr()
        u_eta = 0.0
        for br in ti.static(range(2)):
            ax_ = ti.Matrix.zero(ti.f32, 4, 4)
            ay_ = ti.Matrix.zero(ti.f32, 4, 4)
            az_ = ti.Matrix.zero(ti.f32, 4, 4)
            if ti.static(br == 0):  # fwd (interior: neighbors always exist)
                ax_ = (tensor_field.M_am[i + 1, j, k] - m0) * inv_dx
                ay_ = (tensor_field.M_am[i, j + 1, k] - m0) * inv_dx
                az_ = (tensor_field.M_am[i, j, k + 1] - m0) * inv_dx
            else:  # bwd
                ax_ = (m0 - tensor_field.M_am[i - 1, j, k]) * inv_dx
                ay_ = (m0 - tensor_field.M_am[i, j - 1, k]) * inv_dx
                az_ = (m0 - tensor_field.M_am[i, j, k - 1]) * inv_dx
            fxy = pde.comm_eta44(ax_, ay_)
            fxz = pde.comm_eta44(ax_, az_)
            fyz = pde.comm_eta44(ay_, az_)
            u_eta += 0.5 * 4.0 * (
                pde.inner_eta44(fxy, fxy)
                + pde.inner_eta44(fxz, fxz)
                + pde.inner_eta44(fyz, fyz)
            )
        potential = pde.v4_of(m0, sg, delta, w1)
        observables.energyH_density_aJ[i, j, k] = e_scale * (kinetic + u_eta + potential)


# Note: the global energy aggregate (mean per voxel + grid total) is computed
# by the 3-plane sample_avg_trackers below alongside amp / freq — single pass
# over the slice planes, single GPU↔CPU sync per dashboard cadence. The launcher
# multiplies the returned mean by voxel_count to get the grid-total scalar.


# ================================================================
# FRANK ELASTIC ENERGY DENSITY (M5.1)
# ================================================================
# H_F = (K/2) · |∇n̂|²   per voxel scalar
# |∇n̂|² = Σᵢ Σ_α (∂ᵢ n̂_α)²   (9 central-difference terms: 3 axes × 3 components)
#
# This is the same gradient term that's already implicit in the wave equation's
# c²·∇²ψ Laplacian. The kernel doesn't add new physics — it provides the SCALAR
# E we need for the M5.1 downstream tests:
#   - task 6 (gradient descent):     monitor E falling monotonically
#   - task 7 (Coulomb 1/d fit):       E(d) data points across separations
#   - WAVE_MENU=5 visualization:      color the flux mesh by elastic stress
#
# Naming convention: parallels energyH_density_aJ — quantity is energy density,
# _F suffix tags the Frank formula (vs _H Hamiltonian, future _L Lagrangian).
# K_frank coupling currently hard-coded to 1.0 (Exp 2 baseline); M5.6 plumbs in
# the physical elastic constants alongside the LdG axis hierarchy.


K_FRANK = 1.0  # Frank elastic coupling — Exp 2 baseline; physical scaling in M5.6


@ti.kernel
def compute_energyF_density(
    tensor_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    K_frank: ti.f32,  # type: ignore
):
    """
    Compute per-voxel Frank elastic energy density H_F = (K/2)·|∇n̂|² into
    observables.energyF_density_aJ.

    Reads:  tensor_field.director_nhat  (the matrix-substrate director n̂ = principal
            eigenvector of M; |n̂|=1 by construction. M5.4: repointed from the
            retiring ψ — Frank energy is now the elastic energy of the eigenvector.)
    Writes: observables.energyF_density_aJ

    Same central-difference gradient stencil as compute_energyH_density's
    gradient term, sans the c² factor and sans the kinetic/potential terms.

    Boundary handling: 1-cell halo (gradient stencil); boundary voxels left
    at zero. This is consistent with energyH_density_aJ — both kernels skip
    the outermost voxels.

    PHYSICAL-SCALING TODO — M5.6: the field is named `_aJ` aspirationally but
    actually stores `(K_frank · (am/am)²)` per voxel (dimensionless until the
    LdG elastic constants land in physical units). Sufficient for M5.1 because:
      (a) gradient descent on F depends only on |∇n̂|² gradient (relative)
      (b) Coulomb 1/d fit cares about the FUNCTIONAL form, not absolute scale
    Once M5.6 introduces the biaxial LdG potential with physical elastic
    constants (k1, k2, k3 — splay/twist/bend), apply the appropriate factor
    here so F adds to the H total in the same physical units.
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    half_K = 0.5 * K_frank
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)

    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        # Gradient: 9 central-difference terms — ∂_x/y/z applied to each director component
        d_dx = (tensor_field.director_nhat[i + 1, j, k] - tensor_field.director_nhat[i - 1, j, k]) * inv_2dx
        d_dy = (tensor_field.director_nhat[i, j + 1, k] - tensor_field.director_nhat[i, j - 1, k]) * inv_2dx
        d_dz = (tensor_field.director_nhat[i, j, k + 1] - tensor_field.director_nhat[i, j, k - 1]) * inv_2dx
        grad_n_sqr = d_dx.norm_sqr() + d_dy.norm_sqr() + d_dz.norm_sqr()
        observables.energyF_density_aJ[i, j, k] = half_K * grad_n_sqr


# ================================================================
# EM-FROM-TILTS OBSERVABLES (M5.6.5b — "see EM")
# ================================================================
# Distortion modes of the principal director n̂ that map onto the EM sector
# verified in M5.6.4 (5a §5d): SPLAY ∇·n̂ is Coulomb-charge-like (peaks ±at
# defect cores; for a hedgehog n̂=r̂, ∇·n̂=2/r), TWIST+BEND ∇×n̂ is the B-like
# circulation. Same central-difference stencil as compute_energyF_density.


@ti.kernel
def compute_director_em(
    tensor_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
):
    """Per-voxel ∇·n̂ (signed splay) + ‖∇×n̂‖ (twist+bend magnitude) of director_nhat.

    Reads:  tensor_field.director_nhat
    Writes: observables.director_div_field (signed), observables.director_curl_mag_field (≥0)
    1-cell halo; boundary left at 0 (consistent with energyF/energyH).
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        d_dx = (tensor_field.director_nhat[i + 1, j, k] - tensor_field.director_nhat[i - 1, j, k]) * inv_2dx
        d_dy = (tensor_field.director_nhat[i, j + 1, k] - tensor_field.director_nhat[i, j - 1, k]) * inv_2dx
        d_dz = (tensor_field.director_nhat[i, j, k + 1] - tensor_field.director_nhat[i, j, k - 1]) * inv_2dx
        # div = ∂_x n_x + ∂_y n_y + ∂_z n_z
        observables.director_div_field[i, j, k] = d_dx[0] + d_dy[1] + d_dz[2]
        # curl = (∂_y n_z − ∂_z n_y, ∂_z n_x − ∂_x n_z, ∂_x n_y − ∂_y n_x)
        curl = ti.Vector([d_dy[2] - d_dz[1], d_dz[0] - d_dx[2], d_dx[1] - d_dy[0]])
        observables.director_curl_mag_field[i, j, k] = curl.norm()
        observables.director_curl_field[i, j, k] = curl     # B-direction vector (for glyphs)


@ti.kernel
def compute_director_em_scale(
    tensor_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
):
    """Color-scale maxes from the 3 center planes only (light atomic_max over ~plane voxels,
    NOT a full-grid reduction — avoids the Metal atomic-contention stall, feedback_taichi_metal_atomics).
    Writes observables.director_div_absmax (max|∇·n̂|) + director_curl_max (max‖∇×n̂‖)."""
    observables.director_div_absmax[None] = 1e-12
    observables.director_curl_max[None] = 1e-12
    mid_x, mid_y, mid_z = tensor_field.nx // 2, tensor_field.ny // 2, tensor_field.nz // 2
    for i, j in ti.ndrange((1, tensor_field.nx - 1), (1, tensor_field.ny - 1)):  # XY plane
        ti.atomic_max(observables.director_div_absmax[None], ti.abs(observables.director_div_field[i, j, mid_z]))
        ti.atomic_max(observables.director_curl_max[None], observables.director_curl_mag_field[i, j, mid_z])
    for i, k in ti.ndrange((1, tensor_field.nx - 1), (1, tensor_field.nz - 1)):  # XZ plane
        ti.atomic_max(observables.director_div_absmax[None], ti.abs(observables.director_div_field[i, mid_y, k]))
        ti.atomic_max(observables.director_curl_max[None], observables.director_curl_mag_field[i, mid_y, k])
    for j, k in ti.ndrange((1, tensor_field.ny - 1), (1, tensor_field.nz - 1)):  # YZ plane
        ti.atomic_max(observables.director_div_absmax[None], ti.abs(observables.director_div_field[mid_x, j, k]))
        ti.atomic_max(observables.director_curl_max[None], observables.director_curl_mag_field[mid_x, j, k])


# ================================================================
# VIZ.4 — MAGNETIC-DIPOLE PLACEHOLDER SAMPLE (M5.6.5f stage 1)
# ================================================================
# A *static* hedgehog is a pure ELECTRIC charge (∇×n̂ ≈ 0) — no circulating B,
# no poles — so the bluered N/S coloring (VIZ.2 WM7) and the B glyphs (VIZ.3
# state 3) have nothing to show until a twisting/spinning defect generates a
# real circulating B (the Zitterbewegung clock, M5.8). To validate the RENDER
# path now, this kernel overwrites the curl observables with an ideal analytic
# dipole field. It is a render unit-test source, NOT physics (4b §4.5 / §5.3) —
# the production source swaps in at M5.8 with no render change. Delete or keep
# behind the DIPOLE_SAMPLE flag once the real B exists.


@ti.kernel
def fill_dipole_sample_B(
    tensor_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    m_axis: ti.types.vector(3, ti.f32),  # type: ignore
    cx: ti.f32,  # type: ignore
    cy: ti.f32,  # type: ignore
    cz: ti.f32,  # type: ignore
    r0_vox: ti.f32,  # type: ignore
    amp: ti.f32,  # type: ignore
):
    """Write the ideal magnetic-dipole field B(r) = amp·[3(m̂·r̂)r̂ − m̂]/max(r,r0)³
    (center cx,cy,cz in voxel coords, moment direction m_axis) into the curl
    observables, OVERWRITING compute_director_em's curl output.

    Reads:  nothing (analytic).
    Writes: observables.director_curl_field (B vector, for glyphs + mesh-warp),
            observables.director_curl_mag_field (‖B‖, for orange-magnitude + scale).
    `r0_vox` regularizes the 1/r³ core singularity (B saturates inside r0). `amp`
    folds μ₀/4π·|m| — irrelevant to the self-normalizing color/size scales, kept
    for numerical range. The signed (∇×n̂)·axis coloring (CURL_AXIS = m_axis) then
    shows the dipole signature: RED axial lobes (both poles) + BLUE equatorial band."""
    mhat = m_axis.normalized()
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        r = ti.Vector(
            [ti.cast(i, ti.f32) - cx, ti.cast(j, ti.f32) - cy, ti.cast(k, ti.f32) - cz]
        )
        dist = r.norm()
        rhat = r / (dist + 1e-12)
        reff = ti.max(dist, r0_vox)
        b = amp * (3.0 * mhat.dot(rhat) * rhat - mhat) / (reff * reff * reff)
        observables.director_curl_field[i, j, k] = b
        observables.director_curl_mag_field[i, j, k] = b.norm()


# ================================================================
# WINDING-NUMBER DIAGNOSTIC (M5.1 task 8)
# ================================================================
# Compute the topological charge Q of a director field on a sphere around
# a given center. For a +1 hedgehog Q ≈ +1, anti-hedgehog Q ≈ −1, vacuum
# Q ≈ 0. Pure CPU numpy; called occasionally (once per dashboard frame) as
# a diagnostic, not per-step. Port of Exp 3's `winding_number()`.
#
# Math:
#     Q = (1/4π) ∮_S n̂ · (∂_θ n̂ × ∂_φ n̂) dθ dφ
# Where S is a sphere of given radius around the defect center.
# Method: trilinear-sample ψ on a regular (θ, φ) grid; finite-difference
# the angular derivatives; integrate.


def compute_winding_number(
    psi_np,
    center_vox,
    radius_vox,
    n_theta=32,
    n_phi=64,
):
    """Compute topological winding number Q on a sphere of radius `radius_vox`
    around `center_vox` in voxel coordinates.

    Args:
        psi_np: (nx, ny, nz, 3) numpy array of director-field values
                (typically `tensor_field.director_nhat.to_numpy()`).
        center_vox: (cx, cy, cz) tuple of voxel-coord defect center.
        radius_vox: sphere radius in voxel units. Should be small enough that
                    the sphere fits inside the grid but large enough to be
                    away from the defect core (e.g., 4-10 voxels).
        n_theta: number of polar samples (default 32).
        n_phi: number of azimuthal samples (default 64).

    Returns:
        Q (float): the winding number. ≈ +1 for hedgehog, −1 for anti-hedgehog,
                   0 for vacuum-only enclosure. Discrete grids give Q within
                   ~5% of the integer for a typical relaxed defect.
    """
    nx, ny, nz, _ = psi_np.shape
    cx, cy, cz = center_vox

    # Spherical sample grid
    theta = np.linspace(1e-3, np.pi - 1e-3, n_theta)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    dth = float(theta[1] - theta[0])
    dph = float(phi[1] - phi[0])
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    sx = cx + radius_vox * np.sin(TH) * np.cos(PH)
    sy = cy + radius_vox * np.sin(TH) * np.sin(PH)
    sz = cz + radius_vox * np.cos(TH)

    # Trilinear sample of ψ at the sphere points
    i0 = np.clip(np.floor(sx).astype(int), 0, nx - 2)
    j0 = np.clip(np.floor(sy).astype(int), 0, ny - 2)
    k0 = np.clip(np.floor(sz).astype(int), 0, nz - 2)
    fx = np.clip(sx - i0, 0.0, 1.0)
    fy = np.clip(sy - j0, 0.0, 1.0)
    fz = np.clip(sz - k0, 0.0, 1.0)

    n_s = np.zeros(sx.shape + (3,), dtype=np.float64)
    for dxi in (0, 1):
        for dyi in (0, 1):
            for dzi in (0, 1):
                w = (fx if dxi else 1 - fx) * (fy if dyi else 1 - fy) * (fz if dzi else 1 - fz)
                n_s += w[..., None] * psi_np[i0 + dxi, j0 + dyi, k0 + dzi, :]

    # Re-normalize the sampled directors (interpolation slightly shrinks |n̂|)
    norm = np.linalg.norm(n_s, axis=-1, keepdims=True)
    n_s = n_s / np.maximum(norm, 1e-12)

    # Angular derivatives (central differences; φ wraps periodic)
    dn_dtheta = np.zeros_like(n_s)
    dn_dphi = np.zeros_like(n_s)
    dn_dtheta[1:-1] = (n_s[2:] - n_s[:-2]) / (2 * dth)
    dn_dtheta[0] = (n_s[1] - n_s[0]) / dth
    dn_dtheta[-1] = (n_s[-1] - n_s[-2]) / dth
    dn_dphi[:, :] = (np.roll(n_s, -1, axis=1) - np.roll(n_s, 1, axis=1)) / (2 * dph)

    # Integrand: n̂ · (∂_θ n̂ × ∂_φ n̂). The sinθ factor is implicit in the
    # cross product for a unit-vector field on a parameterized sphere.
    cross = np.cross(dn_dtheta, dn_dphi, axis=-1)
    integrand = (n_s * cross).sum(axis=-1)
    return float(integrand.sum() * dth * dph / (4.0 * np.pi))


def compute_winding_number_guarded(
    psi_np,
    center_vox,
    radii_vox,
    eigenvalues_np=None,
    gap_frac=0.1,
):
    """M5.24 (G9) — the GUARDED winding read per the canonical instrument
    (m5_theory_canonical.md § 5.4): Q at MULTIPLE radii + the eigen-gap
    reliability flag. A single-radius read on a churned state flips branch
    (±1/0 flickers are read artifacts, the § 6 anti-recipe); near-degenerate
    λ₁ ≈ λ₂ regions reshuffle band identity and must DECLINE the read.

    Args:
        psi_np: (nx, ny, nz, 3) director field (director_nhat.to_numpy()).
        center_vox: (cx, cy, cz) defect center, voxel coords.
        radii_vox: iterable of sphere radii (voxels) — use ≥ 3 spanning
            core-adjacent to far-field (e.g. (4, 6, 8, 10)).
        eigenvalues_np: optional (nx, ny, nz, 3) sorted spectrum
            (tensor_field.eigenvalues.to_numpy()); enables the gap guard.
        gap_frac: decline threshold as a fraction of the sphere's median
            λ₁ − λ₂ gap.

    Returns dict: q_mean, q_spread (max − min over radii), per_radius,
    min_gap, declined (True → do NOT quote a charge from this state; the
    caller reports "read declined", never a flickering integer).
    """
    per_radius = []
    min_gap = float("inf")
    declined = False
    nx, ny, nz = psi_np.shape[:3]
    cx, cy, cz = center_vox
    for r in radii_vox:
        per_radius.append(float(compute_winding_number(psi_np, center_vox, r)))
        if eigenvalues_np is not None:
            # nearest-voxel gap sample on the sphere (coarse but sufficient:
            # the guard needs "is any band-crossing region on the sphere",
            # not a precise minimum)
            th = np.linspace(1e-3, np.pi - 1e-3, 16)
            ph = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
            TH, PH = np.meshgrid(th, ph, indexing="ij")
            si = np.clip((cx + r * np.sin(TH) * np.cos(PH)).round().astype(int), 0, nx - 1)
            sj = np.clip((cy + r * np.sin(TH) * np.sin(PH)).round().astype(int), 0, ny - 1)
            sk = np.clip((cz + r * np.cos(TH)).round().astype(int), 0, nz - 1)
            gaps = eigenvalues_np[si, sj, sk, 0] - eigenvalues_np[si, sj, sk, 1]
            g_min, g_med = float(gaps.min()), float(np.median(gaps))
            min_gap = min(min_gap, g_min)
            if g_min < gap_frac * max(g_med, 1e-12):
                declined = True
    q_arr = np.array(per_radius)
    spread = float(q_arr.max() - q_arr.min())
    if spread > 0.2:  # radius-inconsistent read = churned state
        declined = True
    return {
        "q_mean": float(q_arr.mean()),
        "q_spread": spread,
        "per_radius": per_radius,
        "min_gap": (None if eigenvalues_np is None else min_gap),
        "declined": declined,
    }


# ================================================================
# 3-PLANE SAMPLING FOR AVERAGE TRACKERS
# ================================================================
# PERFORMANCE NOTE: Full GPU reduction (atomic_add over all voxels) causes
# severe performance issues due to atomic contention with millions of voxels.
# The 3-plane sampling is a deliberate compromise:
# - Samples ~3N² voxels instead of N³ (e.g., 3% for 100³ grid)
# - Assumes isotropic field distribution (valid for most wave scenarios)
# - Acceptable accuracy vs massive performance gain
#
# AUTO-REDUCE CAVEAT (M5.0h, 2026-05-08): Taichi's "implicit sum reduction"
# pattern (`s = 0; for ... in ti.ndrange: s += ...`) lowers to atomic_add on
# the Metal backend and hits the same contention wall as explicit atomic_add.
# At 63³ × ~5 reductions/step the test pinned the GPU at 100% and stalled
# at ~25 steps/min. Workaround used in m5_0h_dispersion: replace the global
# inner product with sparse point sampling (one voxel read per mode of
# interest, FFT recovers ω from the mixed time series). Don't assume "auto-
# reduction" sidesteps this on Metal — it doesn't.
# ================================================================

# Cached slice buffers — separate caches per domain (Trackers vs FieldObservables)
# so each domain owns its own 3-plane sampling pass. Per `feedback_visual_rendering_priority`,
# each class is self-contained: Trackers samples amp/freq, FieldObservables samples H/F.
_slice_xy_amp = None
_slice_xy_freq = None
_slice_xz_amp = None
_slice_xz_freq = None
_slice_yz_amp = None
_slice_yz_freq = None
_slice_xy_energyH = None
_slice_xy_energyF = None
_slice_xz_energyH = None
_slice_xz_energyF = None
_slice_yz_energyH = None
_slice_yz_energyF = None


# --- TRACKERS slice-copy kernels (amp + freq) -----------------------------


@ti.kernel
def _copy_slice_xy_trackers(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_z: ti.i32,  # type: ignore
):
    """Copy amp / freq at the XY center plane (z=mid_z) to 2D buffers."""
    for i, j in slice_amp:
        slice_amp[i, j] = trackers.amp_local_emarms_am[i, j, mid_z]
        slice_freq[i, j] = trackers.freq_local_cross_rHz[i, j, mid_z]


@ti.kernel
def _copy_slice_xz_trackers(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_y: ti.i32,  # type: ignore
):
    """Copy amp / freq at the XZ center plane (y=mid_y) to 2D buffers."""
    for i, k in slice_amp:
        slice_amp[i, k] = trackers.amp_local_emarms_am[i, mid_y, k]
        slice_freq[i, k] = trackers.freq_local_cross_rHz[i, mid_y, k]


@ti.kernel
def _copy_slice_yz_trackers(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_x: ti.i32,  # type: ignore
):
    """Copy amp / freq at the YZ center plane (x=mid_x) to 2D buffers."""
    for j, k in slice_amp:
        slice_amp[j, k] = trackers.amp_local_emarms_am[mid_x, j, k]
        slice_freq[j, k] = trackers.freq_local_cross_rHz[mid_x, j, k]


# --- FIELD-OBSERVABLES slice-copy kernels (energyH + energyF) -------------


@ti.kernel
def _copy_slice_xy_observables(
    observables: ti.template(),  # type: ignore
    slice_energyH: ti.template(),  # type: ignore
    slice_energyF: ti.template(),  # type: ignore
    mid_z: ti.i32,  # type: ignore
):
    """Copy energyH / energyF at the XY center plane (z=mid_z) to 2D buffers."""
    for i, j in slice_energyH:
        slice_energyH[i, j] = observables.energyH_density_aJ[i, j, mid_z]
        slice_energyF[i, j] = observables.energyF_density_aJ[i, j, mid_z]


@ti.kernel
def _copy_slice_xz_observables(
    observables: ti.template(),  # type: ignore
    slice_energyH: ti.template(),  # type: ignore
    slice_energyF: ti.template(),  # type: ignore
    mid_y: ti.i32,  # type: ignore
):
    """Copy energyH / energyF at the XZ center plane (y=mid_y) to 2D buffers."""
    for i, k in slice_energyH:
        slice_energyH[i, k] = observables.energyH_density_aJ[i, mid_y, k]
        slice_energyF[i, k] = observables.energyF_density_aJ[i, mid_y, k]


@ti.kernel
def _copy_slice_yz_observables(
    observables: ti.template(),  # type: ignore
    slice_energyH: ti.template(),  # type: ignore
    slice_energyF: ti.template(),  # type: ignore
    mid_x: ti.i32,  # type: ignore
):
    """Copy energyH / energyF at the YZ center plane (x=mid_x) to 2D buffers."""
    for j, k in slice_energyH:
        slice_energyH[j, k] = observables.energyH_density_aJ[mid_x, j, k]
        slice_energyF[j, k] = observables.energyF_density_aJ[mid_x, j, k]


# --- Sampler functions ----------------------------------------------------


def sample_avg_trackers(
    tensor_field,
    trackers,
):
    """
    Estimate global Trackers aggregates by sampling 3 orthogonal center planes.

    Computes:
        trackers.amp_global_emarms_am   ← √⟨amp²⟩  (RMS over the 3 planes)
        trackers.freq_global_avg_rHz    ← ⟨freq⟩   (mean over the 3 planes)

    Why 3-plane sampling: full 3D reductions with GPU atomic_add cause severe
    atomic contention at million-voxel scale and stalled the GUI for tens of
    seconds in early M5.0d.2 testing. 3-plane sampling reads ~3·N² voxels
    (3% of a 100³ grid), is exact for spatially-uniform fields, and is a sound
    estimator for smooth wave packets.

    Companion sampler: `sample_avg_observables` (handles energyH / energyF
    via the FieldObservables class). Each domain owns its own pass for clean
    separation of concerns.

    Args:
        tensor_field: TensorField (grid dimensions)
        trackers: Trackers (reads amp/freq per-voxel fields; writes the two
            global aggregates)
    """
    global _slice_xy_amp, _slice_xy_freq
    global _slice_xz_amp, _slice_xz_freq
    global _slice_yz_amp, _slice_yz_freq

    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz

    # Initialize slice buffers once
    if _slice_xy_amp is None:
        _slice_xy_amp = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xy_freq = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xz_amp = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_xz_freq = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_yz_amp = ti.field(dtype=ti.f32, shape=(ny, nz))
        _slice_yz_freq = ti.field(dtype=ti.f32, shape=(ny, nz))

    # Copy 3 center slices to 2D buffers (parallel kernels)
    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2
    _copy_slice_xy_trackers(trackers, _slice_xy_amp, _slice_xy_freq, mid_z)
    _copy_slice_xz_trackers(trackers, _slice_xz_amp, _slice_xz_freq, mid_y)
    _copy_slice_yz_trackers(trackers, _slice_yz_amp, _slice_yz_freq, mid_x)

    # Transfer 2D slices to CPU for numpy operations (exclude boundary voxels)
    xy_amp = _slice_xy_amp.to_numpy()[1:-1, 1:-1]
    xy_freq = _slice_xy_freq.to_numpy()[1:-1, 1:-1]
    xz_amp = _slice_xz_amp.to_numpy()[1:-1, 1:-1]
    xz_freq = _slice_xz_freq.to_numpy()[1:-1, 1:-1]
    yz_amp = _slice_yz_amp.to_numpy()[1:-1, 1:-1]
    yz_freq = _slice_yz_freq.to_numpy()[1:-1, 1:-1]

    n_samples = xy_amp.size + xz_amp.size + yz_amp.size

    # RMS amplitude: √(⟨A²⟩) — energy-weighting convention
    total_amp_squared = (xy_amp**2).sum() + (xz_amp**2).sum() + (yz_amp**2).sum()
    trackers.amp_global_emarms_am[None] = float(np.sqrt(total_amp_squared / n_samples))

    # Plain mean for frequency
    total_freq = xy_freq.sum() + xz_freq.sum() + yz_freq.sum()
    trackers.freq_global_avg_rHz[None] = float(total_freq / n_samples)


def sample_avg_observables(
    tensor_field,
    observables,
):
    """
    Estimate global FieldObservables aggregates by sampling 3 orthogonal planes.

    Computes:
        observables.energyH_global_avg_aJ ← ⟨E_H⟩  (mean Hamiltonian density)
        observables.energyF_global_avg_aJ ← ⟨E_F⟩  (mean Frank elastic density)

    Caller (launcher) derives the grid total trivially:
        energyH_total_aJ = energyH_global_avg_aJ × voxel_count

    Same 3-plane sampling rationale as `sample_avg_trackers`. Separate pass so
    FieldObservables is self-contained.

    Args:
        tensor_field: TensorField (grid dimensions)
        observables: FieldObservables (reads energyH/energyF per-voxel fields;
            writes the two global aggregates)
    """
    global _slice_xy_energyH, _slice_xy_energyF
    global _slice_xz_energyH, _slice_xz_energyF
    global _slice_yz_energyH, _slice_yz_energyF

    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz

    # Initialize slice buffers once
    if _slice_xy_energyH is None:
        _slice_xy_energyH = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xy_energyF = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xz_energyH = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_xz_energyF = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_yz_energyH = ti.field(dtype=ti.f32, shape=(ny, nz))
        _slice_yz_energyF = ti.field(dtype=ti.f32, shape=(ny, nz))

    # Copy 3 center slices to 2D buffers (parallel kernels)
    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2
    _copy_slice_xy_observables(observables, _slice_xy_energyH, _slice_xy_energyF, mid_z)
    _copy_slice_xz_observables(observables, _slice_xz_energyH, _slice_xz_energyF, mid_y)
    _copy_slice_yz_observables(observables, _slice_yz_energyH, _slice_yz_energyF, mid_x)

    # Transfer 2D slices to CPU for numpy operations (exclude boundary voxels)
    xy_energyH = _slice_xy_energyH.to_numpy()[1:-1, 1:-1]
    xy_energyF = _slice_xy_energyF.to_numpy()[1:-1, 1:-1]
    xz_energyH = _slice_xz_energyH.to_numpy()[1:-1, 1:-1]
    xz_energyF = _slice_xz_energyF.to_numpy()[1:-1, 1:-1]
    yz_energyH = _slice_yz_energyH.to_numpy()[1:-1, 1:-1]
    yz_energyF = _slice_yz_energyF.to_numpy()[1:-1, 1:-1]

    n_samples = xy_energyH.size + xz_energyH.size + yz_energyH.size

    # Plain mean for energy density H (per voxel); launcher multiplies by
    # voxel_count to get the grid total.
    total_energyH = xy_energyH.sum() + xz_energyH.sum() + yz_energyH.sum()
    observables.energyH_global_avg_aJ[None] = float(total_energyH / n_samples)

    # Plain mean for Frank elastic energy density F (per voxel).
    total_energyF = xy_energyF.sum() + xz_energyF.sum() + yz_energyF.sum()
    observables.energyF_global_avg_aJ[None] = float(total_energyF / n_samples)
