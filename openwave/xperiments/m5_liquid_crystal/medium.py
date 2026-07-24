"""
LIQUID-CRYSTAL Model Medium Data-Grid

Object Classes @spacetime module.

LIQUID-CRYSTAL evolves the Landau–de Gennes order parameter M(x) via the Eq.18
matrix action (M5.5): the leapfrog ∂²_t M = c²·∇·G(M) − ∂V_M/∂M, where G is the
curvature flux. M is the 4×4 (M5.8) real-symmetric tensor
    M(x) = O(x)·D·O^T(x),   D = diag(g, 1, δ, 0)   (index-0: time/g axis = array index 0)
with O(x) the per-voxel dynamical SO(1,3) frame and D the frozen spectrum (the M5.8
time/boost axis g at index 0 + the spatial 3×3 block diag(1,δ,0) at indices 1,2,3;
Duda convention, see research/_convention_refactor/CONVENTION.md). Time evolution uses a
matrix triple-buffer leapfrog (M_prev_am, M_am, M_new_am).

────────────────────────────────────────────────────────────────────────────
SUBSTRATE HISTORY — Vector(3) ψ → matrix M (M5.4 migration, M5.8 ψ-retire)
────────────────────────────────────────────────────────────────────────────
M5.2 closed as an informative negative: the Vector(3) director ψ cannot host the
paper's physics, so M5.4 migrated the substrate to the LdG order parameter M
(feasibility: research/sandbox_v3/m5_3_matrix_feasibility.py). The M5.8 ψ-retire
then DELETED the legacy Vector(3) ψ triple-buffer and its wave engine (evolve_psi,
the ∇²/∇·/∇× operators, V_psi) entirely. The live substrate is:
  - M_am / M_prev_am / M_new_am  — the 4×4 matrix triple buffer (the substrate)
  - director_nhat / eigenvalues  — DERIVED per-frame via eigen_decompose (the
    principal eigenvector of M's spatial block is the director the render + M5.1
    Coulomb machinery consume). See m5_4b_rendering_features.md.
"""

import math

import taichi as ti

from openwave.common import colormap, constants, utils

# Uniaxial minor-axis eigenvalue for M5.4 (placeholder). The order parameter is
# M = δ·I + (1−δ)·n̂⊗n̂ → eigenvalues (1, δ, δ), principal eigenvector = n̂.
# The physical biaxial δ ~ ℏ (and the gravity boost g) land in M5.6 / M5.8; M5.4
# only needs a uniaxial spectrum to reproduce the M5.1 Coulomb topology on M.
LC_DELTA = 0.5

# M5.8.1 — 4D substrate promotion. The matrix field is 4×4: the time/boost axis at
# INDEX 0 (eigenvalue g, the paper's D=diag(g,1,δ,0); Duda convention, flipped from the
# legacy index-3 ordering 2026-06-21 — research/_convention_refactor/CONVENTION.md) PLUS
# the spatial 3×3 block (indices 1,2,3 = x,y,z — the M5.4/M5.6 order parameter). Callers
# extract the spatial block [1:4,1:4] for the 3×3-Cardano eigensolver (the director is its
# principal eigenvector). In M5.8.1 g is a CONSTANT, boost-DECOUPLED background
# (M = block-diag(g, M_spatial)) so the spatial dynamics are identical to 3D (validated:
# sandbox_v8/m5_8_1_4x4_promotion.py); the Minkowski coupling driving the clock lands in M5.8.2.
MDIM = 4          # matrix substrate dimension (was 3)
LC_G = 8.0        # time-axis (boost/gravity) eigenvalue g ≫ 1, constant in M5.8.1

# VIZ.5 (M5.23) — the "ellipsoid" visualization: one eigenframe glyph per 3D
# angle on an S² shell around each defect center — a FULL-3D view (unlike the
# 3-plane flux-mesh cross-sections); design + decisions:
# research/tasks/m5_23_task_details.md. Buffer ceilings only — the live glyph
# count is the GUI density slider (capped at ELLIPSOID_MAX_DIRS, honoring the
# not-too-dense readability spec), and directions are computed IN-KERNEL as a
# Fibonacci-sphere lattice (taichi-native, no host-side table). Centers:
# dipole configs use 2; headroom for small clusters.
ELLIPSOID_MAX_DIRS = 642
ELLIPSOID_MAX_CENTERS = 4
# Stage B mesh template: a UV sphere generated in-kernel (closed form per
# index — poles + LAT rings × LON meridians; faces = 2·LAT·LON triangles).
# It is a per-glyph SURFACE template (pole clustering is irrelevant here,
# unlike the Fibonacci SAMPLING lattice above); each shell sample maps it
# through M·u to the local eigen-ellipsoid (the m5_6_5b geometry).
ELLIPSOID_TEMPLATE_LAT = 8
ELLIPSOID_TEMPLATE_LON = 12
ELLIPSOID_TVERTS = ELLIPSOID_TEMPLATE_LAT * ELLIPSOID_TEMPLATE_LON + 2
ELLIPSOID_TFACES = 2 * ELLIPSOID_TEMPLATE_LAT * ELLIPSOID_TEMPLATE_LON

# VIZ.5 Stage D — the disclination-rod render: ellipsoid samples ALONG the rod
# axis (the ±z column through each center: the escaped-core disclination pair
# every biaxial seed constructs) plus per-2D-angle RINGS around the cord at
# fixed heights. Fixed v1 layout (course-correctable): ROD_N axis samples
# spanning ± the GUI Radius, RING_COUNT rings (heights ±0.35 / ±0.7 of that
# half-length) × RING_AZ azimuths each.
ELLIPSOID_ROD_N = 41
ELLIPSOID_RING_COUNT = 8  # 4 rows per pole (the reference figure shows 4), on the outer rod sections
ELLIPSOID_RING_AZ = 64  # CEILING azimuths per ring; the live count follows the GUI Count slider (~Count/12)
ELLIPSOID_ROD_SLOTS = ELLIPSOID_ROD_N + ELLIPSOID_RING_COUNT * ELLIPSOID_RING_AZ

# M5.23.2 arm (4) — the energy-density ISOSURFACE: triangle budget for the
# marching-tetrahedra mesh (tri soup; unused slots collapse to the origin).
# A closed surface in a 63³ arena is typically 20-40k triangles (MT emits
# ~2× marching cubes); overflow beyond the budget is COUNTED and the
# launcher WARNS (no silent truncation — the M5.23.2 plan's honesty gate).
ISO_MAX_TRIS = 131072


@ti.func
def _uv_face_locals(f: ti.i32):  # type: ignore
    """Template-local vertex ids (i0, i1, i2) of UV-sphere face `f`: top fan,
    then the quad bands split in two, then the bottom fan. Shared by the shell
    and rod index kernels (two_sided render ⇒ winding order not load-bearing)."""
    lat, lon = ELLIPSOID_TEMPLATE_LAT, ELLIPSOID_TEMPLATE_LON
    i0, i1, i2 = 0, 0, 0
    if f < lon:  # top fan: north pole + ring 0
        b = f
        i0 = 0
        i1 = 1 + b
        i2 = 1 + ((b + 1) % lon)
    elif f < lon + 2 * (lat - 1) * lon:  # quad bands between rings
        q = (f - lon) // 2
        t = (f - lon) - q * 2
        a = q // lon
        b = q - a * lon
        p00 = 1 + a * lon + b
        p01 = 1 + a * lon + ((b + 1) % lon)
        p10 = 1 + (a + 1) * lon + b
        p11 = 1 + (a + 1) * lon + ((b + 1) % lon)
        if t == 0:
            i0, i1, i2 = p00, p10, p01
        else:
            i0, i1, i2 = p01, p10, p11
    else:  # bottom fan: south pole + last ring
        b = f - lon - 2 * (lat - 1) * lon
        ring = 1 + (lat - 1) * lon
        i0 = ELLIPSOID_TVERTS - 1
        i1 = ring + ((b + 1) % lon)
        i2 = ring + b
    return ti.Vector([i0, i1, i2])


@ti.data_oriented
class TensorField:
    """
    Lagrangian tensor field on a cell-centered grid with attometer scaling.

    Matrix triple-buffer leapfrog convention:
        M_prev_am — M at t−dt (history buffer; read by stencil)
        M_am      — M at t    (current buffer; read by stencil)
        M_new_am  — M at t+dt (output buffer; written by evolve_M)

    After each leapfrog step, swap_matrix_buffers() cycles: prev ← curr, curr ← new.
    M_new_am is overwritten by the next step.

    AMR-readiness convention:
        Kernels MUST read grid dimensions via tensor_field.nx / .ny / .nz attributes
        (or the .grid_size tuple). Do NOT bake fixed (nx, ny, nz) constants into
        @ti.kernel signatures — that would prevent the M5.6 / M5.8 AMR retrofit
        from swapping in an octree-based field-storage layer without rewriting
        kernels. Field access via tensor_field.M_am[i, j, k] is the canonical pattern.

    This class:
    - Cell-centered cubic grid
    - Attometer scaling for numerical precision (f32 fields)
    - Computed positions from indices (memory efficient)
    - Wave properties stored at each voxel
    - Asymmetric universe support (nx ≠ ny ≠ nz allowed)

    Initialization Strategy:
    1. User specifies init_universe_size [x, y, z] in meters (can be asymmetric)
    2. Compute universe volume and target voxel count
    3. Calculate cubic voxel size: dx = (volume / target_voxels)^(1/3)
    4. Compute grid dimensions: nx = int(x_size / dx), ny = int(y_size / dx), nz = int(z_size / dx)
    5. Recalculate actual universe size to fit integer voxel counts
    6. Initialize scalar and vector fields with attometer scaling for f32 precision
    """

    def __init__(
        self,
        init_universe_size,
        target_voxels,
        flux_mesh_planes=[0.5, 0.5, 0.5],
        viz_stride=4,
    ):
        """
        Initialize TensorField from universe size with automatic voxel sizing.

        Args:
            init_universe_size: Simulation domain size [x, y, z] in meters.
                Can be asymmetric. Will be adjusted to fit integer voxel counts.
            target_voxels: Desired total voxel count (impacts memory and performance).

        Note:
            Voxel size (dx) is cubic (same for all axes) to preserve wave physics.
            Grid counts (nx, ny, nz) can differ for asymmetric domain shapes.
        """
        # Compute initial grid properties (before rounding and grid symmetry)
        init_universe_volume = (
            init_universe_size[0] * init_universe_size[1] * init_universe_size[2]
        )

        # Calculate cubic voxel size from target voxel count
        # CRITICAL: voxels must remain cubic (same edge length on all axes)
        # This preserves wave physics isotropy. Only the NUMBER of voxels varies per axis.
        self.voxel_volume = init_universe_volume / target_voxels  # cubic voxels
        self.voxel_edge = self.voxel_volume ** (1 / 3)  # same as dx, dx³ = voxel volume
        self.voxel_edge_am = self.voxel_edge / constants.ATTOMETER  # in attometers
        self.dx = self.voxel_edge  # additional alias for simplicity
        self.dx_am = self.voxel_edge_am  # additional alias for simplicity

        # Calculate grid dimensions (number of complete voxels per dimension) - asymmetric
        # Uses nearest odd integer to ensure grid symmetry with unique central voxel:
        # 1. User-specified universe size is arbitrary (any float value)
        # 2. voxel_edge comes from cube root, rarely divides evenly into universe size
        # 3. Ensures integer count needed for array indexing and loop bounds
        # 4. Rounds to nearest odd integer for symmetric grid with central voxel
        # 5. Actual universe size recalculated below to fit integer voxel count
        self.grid_size = [
            utils.round_to_nearest_odd(init_universe_size[0] / self.dx),
            utils.round_to_nearest_odd(init_universe_size[1] / self.dx),
            utils.round_to_nearest_odd(init_universe_size[2] / self.dx),
        ]  # same as (nx, ny, nz)
        self.nx = self.grid_size[0]  # additional alias for simplicity
        self.ny = self.grid_size[1]  # additional alias for simplicity
        self.nz = self.grid_size[2]  # additional alias for simplicity
        self.max_grid_size = max(self.nx, self.ny, self.nz)
        self.min_grid_size = min(self.nx, self.ny, self.nz)

        # Compute total voxels (asymmetric grid)
        self.voxel_count = self.nx * self.ny * self.nz

        # Recompute actual universe dimensions to fit integer number of cubic voxels
        self.universe_size = [self.nx * self.dx, self.ny * self.dx, self.nz * self.dx]
        self.universe_size_am = [size / constants.ATTOMETER for size in self.universe_size]
        self.max_universe_edge = max(self.nx * self.dx, self.ny * self.dx, self.nz * self.dx)
        self.max_universe_edge_am = self.max_universe_edge / constants.ATTOMETER
        self.universe_volume = self.voxel_count * self.voxel_volume

        # ================================================================
        # DATA STRUCTURE & INITIALIZATION
        # ================================================================
        # The propagated field is the matrix order parameter M (declared below);
        # the M5.8 ψ-retire removed the legacy Vector(3) ψ triple-buffer (psi_am /
        # psi_prev_am / psi_new_am) and its swap_buffers — the live substrate is the
        # 4×4 M triple-buffer (M_am / M_prev_am / M_new_am) with swap_matrix_buffers.
        self.position_render = ti.Vector.field(3, dtype=ti.f32, shape=(self.nx * self.ny))  # flat

        # ================================================================
        # MATRIX-FIELD SUBSTRATE (M5.4 order parameter, M5.8.1 4×4 promotion)
        # ================================================================
        # M(x) = O(x)·D·O^T(x), real-symmetric MDIM×MDIM = 4×4 since the M5.8.1
        # promotion: the time/boost axis g at index 0 PLUS the spatial 3×3 block (the M5.4
        # LdG order parameter diag(1,δ,0)) at indices 1,2,3, i.e. D=diag(g,1,δ,0) (Duda
        # index-0, flipped 2026-06-21 — research/_convention_refactor/CONVENTION.md). Stored
        # as the full 16-component ti.Matrix.field(MDIM,MDIM) (M5.4 decision 2026-05-26: the
        # dense matrix matches the proven feasibility spike, ti.sym_eig-ready, no
        # reassembly; the symmetric-packing optimization can swap in behind this
        # accessor later if memory binds). Triple buffer mirrors the psi leapfrog
        # convention — M_prev/M/M_new — for the matrix evolution (M5.5 onward);
        # eigen_decompose reads the spatial [1:4,1:4] block, the time axis (index 0) stays g.
        self.M_am = ti.Matrix.field(MDIM, MDIM, dtype=ti.f32, shape=self.grid_size)  # M at t (4×4: time idx 0 + spatial 1:4)
        self.M_prev_am = ti.Matrix.field(MDIM, MDIM, dtype=ti.f32, shape=self.grid_size)  # M at t−dt
        self.M_new_am = ti.Matrix.field(MDIM, MDIM, dtype=ti.f32, shape=self.grid_size)  # M at t+dt

        # DERIVED per-frame from M_am by engine2_pde.eigen_decompose (the lynchpin
        # kernel — see m5_4b_rendering_features.md). director_nhat is the principal
        # eigenvector (the apolar nematic director the M5.1 Frank/Coulomb/glyph
        # machinery consumes); eigenvalues are the (λ₁≥λ₂≥λ₃) spectrum, used by the
        # ‖M−D‖_F amplitude tracker and the biaxial-ellipsoid render option.
        # These are stateless caches: recomputed every frame, valid even when paused.
        self.director_nhat = ti.Vector.field(3, dtype=ti.f32, shape=self.grid_size)  # eigvec
        # VIZ.3: MIDDLE (δ) eigenvector — the "clock-hand" axis that sweeps around the
        # principal director under the Zitterbewegung twist (the δ-clock-hand glyph).
        # Derived per-frame in eigen_decompose alongside director_nhat; apolar (n̂≡−n̂).
        self.director_mid = ti.Vector.field(3, dtype=ti.f32, shape=self.grid_size)  # mid eigvec
        self.eigenvalues = ti.Vector.field(3, dtype=ti.f32, shape=self.grid_size)  # λ₁≥λ₂≥λ₃/voxel
        # Working buffer for the director-equivalent Coulomb relaxation (M5.4 gate):
        # relax_director_step writes the next director here, then M is rebuilt from it.
        self.director_nhat_new = ti.Vector.field(3, dtype=ti.f32, shape=self.grid_size)
        self.lc_delta = LC_DELTA  # uniaxial minor-axis eigenvalue (M = δI + (1−δ)n̂⊗n̂)
        self.lc_g = LC_G  # time-axis (index 0) eigenvalue g, D=diag(g,1,δ,0); read by kernels via tensor_field.lc_g

        # M5.5.4 — Eq.18 matrix-action leapfrog (simple ½‖Ṁ‖² kinetic + faithful
        # potential U = 4Σ‖[M_μ,M_ν]‖² + V(M)). The curvature flux
        #     G_α = ∂U_curv/∂M_α = 8 Σ_ν [[M_α,M_ν], M_ν]   (symmetric)
        # is computed everywhere by compute_curvature_flux from M_am; evolve_M then
        # takes its divergence Σ_α ∂_α G_α as the force. Two-pass (G field, then
        # divergence) keeps each kernel a 1-cell halo instead of one 2-cell mega-kernel.
        self.curv_flux_x = ti.Matrix.field(MDIM, MDIM, dtype=ti.f32, shape=self.grid_size)
        self.curv_flux_y = ti.Matrix.field(MDIM, MDIM, dtype=ti.f32, shape=self.grid_size)
        self.curv_flux_z = ti.Matrix.field(MDIM, MDIM, dtype=ti.f32, shape=self.grid_size)

        # M5.24 round 3 — FIRE-relaxer per-z-slice partial sums (F·V, F·F, V·V);
        # a two-stage Metal-safe reduction (63 slice threads, no full-grid
        # atomics — the M5.0h contention lesson). Host sums the (3, nz) array.
        self.fire_partials = ti.field(dtype=ti.f32, shape=(3, self.nz))

        # M5.8.2c — 4D Minkowski evolution support (flag-gated; OFF unless the
        # xperiment seeds a DRESSED hedgehog). stable_mask (1.0 = Minkowski-signed
        # flux, 0.0 = Euclidean fallback) is the per-voxel ghost guard: the signed
        # kernel's gradient stiffness Q(∇ψ) loses positive-definiteness on the fuel
        # shell (M5.8.2b-2), where a linearized signed evolution is ill-posed —
        # those voxels evolve with the (always-stable) Euclidean flux instead.
        # Computed ONCE post-seed by engine2_pde.compute_stable_mask (seed-time
        # cost only). M_psi_am is the seed clock tangent ∂M/∂Θ (the (δ,0)-plane
        # rotation pattern) — consumed by the mask kernel + the clock trackers.
        self.stable_mask = ti.field(dtype=ti.f32, shape=self.grid_size)
        self.M_psi_am = ti.Matrix.field(MDIM, MDIM, dtype=ti.f32, shape=self.grid_size)
        # 3-plane-sampled coherent (α,3) velocity sums (engine2_pde.sample_v03_drift)
        self.v03_sums = ti.field(dtype=ti.f32, shape=3)
        # per-voxel V(M) amplitude target t*(x) = Tr(M_sp²) of the SEED (the 2c-0
        # design input "pin the well to the DRESSED amplitude"): the dressed seed
        # then sits exactly AT the V minimum — zero static V-force, confinement
        # kept. Seeded once by engine2_pde.compute_tstar (the 2c-2 V-pump fix).
        self.ldg_tstar = ti.field(dtype=ti.f32, shape=self.grid_size)

        # M5.8.2c OPTION B — the CONSTRAINED spectral-projection integrator state
        # (INTEGRATOR_4D: "constrained"): the Minkowski-SIGNED dynamics evolves
        # the pair (M, P) — P_am is the canonical momentum density, Md_am the
        # solved velocity Ṁ = A⁺P (read by the kinetic flux term + diagnostics),
        # act4d the Dirichlet evolve mask (interior ∧ off-core ∧ off-axis — the
        # 2c-1 validated boundary treatment), sym_basis the 10 orthonormal
        # symmetric 4×4 basis matrices (filled once by the launcher). Validated
        # against the f64 numpy reference in sandbox_v8/m5_8_2cb (B-1 gates).
        self.P_am = ti.Matrix.field(MDIM, MDIM, dtype=ti.f32, shape=self.grid_size)
        self.Md_am = ti.Matrix.field(MDIM, MDIM, dtype=ti.f32, shape=self.grid_size)
        self.act4d = ti.field(dtype=ti.f32, shape=self.grid_size)
        self.sym_basis = ti.Matrix.field(MDIM, MDIM, dtype=ti.f32, shape=10)

        # TODO: check need for velocity field = pressure / density
        # Wave velocity vector field (v = dψ/dt)
        # self.velocity_am = ti.Vector.field(3, dtype=ti.f32, shape=self.grid_size)  # am/s

        # TODO: Implement DERIVED SCALAR FIELDS
        # WAVELENGTH, PERIOD, phase, momentum

        # TODO: Implement DERIVED VECTOR FIELDS (directions normalized to unit vectors)
        # energy_flux, wave_direction, displacement_direction, wave_mode, wave_type

        # ================================================================
        # Grid Lines: data structures & initialization
        # ================================================================
        # Grid: optimized grid lines for rendering
        # Each line spans the entire grid dimension (e.g., from x=0 to x=1 in normalized coords)

        # Line count per direction:
        # - X-direction (parallel to X): (ny+1) × (nz+1) lines
        # - Y-direction (parallel to Y): (nx+1) × (nz+1) lines
        # - Z-direction (parallel to Z): (nx+1) × (ny+1) lines
        # Total vertices = 2 × (sum of lines)
        self.line_count = (
            (self.ny + 1) * (self.nz + 1)  # X-parallel lines
            + (self.nx + 1) * (self.nz + 1)  # Y-parallel lines
            + (self.nx + 1) * (self.ny + 1)  # Z-parallel lines
        )
        self.grid_lines = ti.Vector.field(3, dtype=ti.f32, shape=self.line_count * 2)
        self.edge_lines = ti.Vector.field(3, dtype=ti.f32, shape=24)  # 12 edges, 2 vertices each

        # Populate the grid and edge lines with normalized positions (ready for rendering)
        self.populate_grid_lines()  # initialize grid lines (already normalized)
        self.populate_edge_lines()  # initialize edge lines (already normalized)

        # ================================================================
        # Flux Mesh: data structures & initialization
        # ================================================================
        # Three orthogonal meshes intersecting at universe center
        # Each flux mesh has resolution matching simulation voxel grid
        # Vertices positioned at voxel centers for direct property sampling

        # XY Flux Mesh (at z = center): spans x and y dimensions
        # - Vertices: nx × ny grid points
        # - Indices: (nx-1) × (ny-1) quads × 6 indices (2 triangles each)
        self.fluxmesh_xy_vertices = ti.Vector.field(n=3, dtype=ti.f32, shape=(self.nx, self.ny))
        self.fluxmesh_xy_colors = ti.Vector.field(n=3, dtype=ti.f32, shape=(self.nx, self.ny))
        self.fluxmesh_xy_indices = ti.field(dtype=ti.i32, shape=(self.nx - 1, self.ny - 1, 6))

        # XZ Flux Mesh (at y = center): spans x and z dimensions
        # - Vertices: nx × nz grid points
        # - Indices: (nx-1) × (nz-1) quads × 6 indices (2 triangles each)
        self.fluxmesh_xz_vertices = ti.Vector.field(n=3, dtype=ti.f32, shape=(self.nx, self.nz))
        self.fluxmesh_xz_colors = ti.Vector.field(n=3, dtype=ti.f32, shape=(self.nx, self.nz))
        self.fluxmesh_xz_indices = ti.field(dtype=ti.i32, shape=(self.nx - 1, self.nz - 1, 6))

        # YZ Flux Mesh (at x = center): spans y and z dimensions
        # - Vertices: ny × nz grid points
        # - Indices: (ny-1) × (nz-1) quads × 6 indices (2 triangles each)
        self.fluxmesh_yz_vertices = ti.Vector.field(n=3, dtype=ti.f32, shape=(self.ny, self.nz))
        self.fluxmesh_yz_colors = ti.Vector.field(n=3, dtype=ti.f32, shape=(self.ny, self.nz))
        self.fluxmesh_yz_indices = ti.field(dtype=ti.i32, shape=(self.ny - 1, self.nz - 1, 6))

        # Initialize flux mesh (vertices, indices, colors)
        self.flux_mesh_planes = flux_mesh_planes  # positions in normalized coords
        self.fm_plane_x_idx = int(flux_mesh_planes[0] * self.nx)  # x index of flux mesh plane
        self.fm_plane_y_idx = int(flux_mesh_planes[1] * self.ny)  # y index of flux mesh plane
        self.fm_plane_z_idx = int(flux_mesh_planes[2] * self.nz)  # z index of flux mesh plane
        # Continuous normalized plane positions — the SAME formula create_flux_mesh()
        # uses for the mesh perpendicular coord. Glyphs read these (not idx/max_dim)
        # so glyph origins sit exactly on the mesh sheet for any grid / plane fraction.
        _md = float(self.max_grid_size)
        self.fm_plane_x_pos = flux_mesh_planes[0] * (self.nx / _md)  # x perpendicular, normalized
        self.fm_plane_y_pos = flux_mesh_planes[1] * (self.ny / _md)  # y perpendicular, normalized
        self.fm_plane_z_pos = flux_mesh_planes[2] * (self.nz / _md)  # z perpendicular, normalized

        self.create_flux_mesh()

        # ================================================================
        # Director Glyphs (M5.1): line-segment visualization of n̂ = ψ/|ψ|
        # ================================================================
        # The flux mesh renders one scalar per voxel (magnitude or signed
        # component); for a director field, the *direction* is the structure
        # we need to see. Director glyphs draw a line segment from each
        # sampled voxel to voxel + L·n̂, with color = signed-component RGB
        # so opposite directions look opposite (red↔cyan, green↔magenta,
        # blue↔yellow). Design doc: research/2b_director_glyph_rendering.md.
        #
        # Same 3-plane sampling pattern as flux_mesh, sampled every VIZ_STRIDE
        # voxels. VIZ_STRIDE is the consolidated xparameter that ALSO drives
        # the SHOW_GRANULES sampling stride — same density for both overlays
        # so granules and glyphs visually align (granules show medium motion
        # at each sample point; glyphs show director orientation at the same
        # point). Round-up `(n + stride - 1) // stride` mirrors the granule
        # indexing in _launcher.py so last-row coverage matches.
        self.GLYPH_STRIDE = viz_stride
        nx_s = max(1, (self.nx + viz_stride - 1) // viz_stride)
        ny_s = max(1, (self.ny + viz_stride - 1) // viz_stride)
        nz_s = max(1, (self.nz + viz_stride - 1) // viz_stride)
        self.glyph_nx_s = nx_s
        self.glyph_ny_s = ny_s
        self.glyph_nz_s = nz_s
        self.n_glyphs = nx_s * ny_s + nx_s * nz_s + ny_s * nz_s
        # Two vertices per line segment (start + end); GGUI renders consecutive
        # vertex pairs as a line.
        self.director_glyph_vertices = ti.Vector.field(3, ti.f32, (2 * self.n_glyphs))
        # Per-vertex color matches the line endpoints. Both endpoints get the
        # same color so each glyph reads as a uniform-colored line.
        self.director_glyph_colors = ti.Vector.field(3, ti.f32, (2 * self.n_glyphs))
        # Flat-array offsets into the (2·n_glyphs) buffer per plane:
        self.glyph_offset_xy = 0
        self.glyph_offset_xz = nx_s * ny_s
        self.glyph_offset_yz = nx_s * ny_s + nx_s * nz_s

        # Half-arrowhead barb segments — one extra line per glyph from the
        # shaft tip to a perpendicular offset, so direction (head vs tail) is
        # visually unambiguous. Same indexing scheme + plane offsets as the
        # main shaft buffer; rendered separately by the launcher (gated on
        # engine4_render.SHOW_DIRECTOR_ARROWHEAD) so it can be toggled
        # without re-running the kernel.
        self.director_glyph_arrow_vertices = ti.Vector.field(3, ti.f32, (2 * self.n_glyphs))
        self.director_glyph_arrow_colors = ti.Vector.field(3, ti.f32, (2 * self.n_glyphs))

        # VIZ.4 — single magnetic-MOMENT vector glyph at the defect center (the
        # placeholder dipole sample). 4 vertices: shaft (base→tip) + barb (tip→back).
        # Static (set once by engine4_render.update_moment_glyph when DIPOLE_SAMPLE).
        self.moment_glyph_vertices = ti.Vector.field(3, ti.f32, 4)
        self.moment_glyph_colors = ti.Vector.field(3, ti.f32, 4)

        # VIZ.5 (M5.23) — the "ellipsoid" viz: one eigen-ellipsoid per 3D angle
        # on an S² shell of GUI radius R around each defect center (a FULL-3D
        # view, not a plane cross-section). Centers are filled at seed time by
        # the launcher (voxel coords; count in ellipsoid_n_centers). The line-
        # glyph Stage A variant was retired 2026-07-19 (maintainer verdict:
        # mesh is the better viz); the mesh buffers below are the feature.
        self.ellipsoid_max_dirs = ELLIPSOID_MAX_DIRS
        self.ellipsoid_max_centers = ELLIPSOID_MAX_CENTERS
        self.ellipsoid_centers = ti.Vector.field(3, ti.f32, ELLIPSOID_MAX_CENTERS)  # voxel coords
        self.ellipsoid_n_centers = ti.field(ti.i32, shape=())

        # VIZ.5 Stage B — the M·u ellipsoid MESH: per shell sample the unit-sphere
        # template maps through the local M (spatial block + a small isotropic
        # visibility floor) to the eigen-ellipsoid surface — semi-axes = the
        # eigenvalues along the eigenvectors, NO eigendecomposition needed
        # (m5_6_5b key simplification). One flat vertex/color/index pool spans
        # all (center, direction) slots so the launcher renders the whole shell
        # set with ONE scene.mesh call; unused slots collapse to the origin
        # (zero-area triangles, invisible). Indices are static (filled once);
        # vertices + colors are rewritten per-frame by
        # engine4_render.update_ellipsoid_mesh.
        self.ellipsoid_tverts = ELLIPSOID_TVERTS
        self.ellipsoid_tfaces = ELLIPSOID_TFACES
        self.ellipsoid_template = ti.Vector.field(3, ti.f32, ELLIPSOID_TVERTS)
        n_slots = ELLIPSOID_MAX_CENTERS * ELLIPSOID_MAX_DIRS
        self.ellipsoid_mesh_vertices = ti.Vector.field(3, ti.f32, n_slots * ELLIPSOID_TVERTS)
        self.ellipsoid_mesh_colors = ti.Vector.field(3, ti.f32, n_slots * ELLIPSOID_TVERTS)
        self.ellipsoid_mesh_indices = ti.field(ti.i32, n_slots * ELLIPSOID_TFACES * 3)
        self.populate_ellipsoid_template()
        self.populate_ellipsoid_mesh_indices()

        # VIZ.5 Stage D — the ROD pool: rod-line + rod-ring ellipsoid slots per
        # center (the disclination-rod render), same template + M·u map as the
        # shell pool, rendered as ONE additional scene.mesh call. Indices are
        # static (filled once); vertices + colors are rewritten per-frame by
        # engine4_render.update_rod_ellipsoids (inactive slots collapse).
        self.ellipsoid_rod_n = ELLIPSOID_ROD_N
        self.ellipsoid_ring_count = ELLIPSOID_RING_COUNT
        self.ellipsoid_ring_az = ELLIPSOID_RING_AZ
        self.ellipsoid_rod_slots = ELLIPSOID_ROD_SLOTS
        n_rod = ELLIPSOID_MAX_CENTERS * ELLIPSOID_ROD_SLOTS
        self.ellipsoid_rod_vertices = ti.Vector.field(3, ti.f32, n_rod * ELLIPSOID_TVERTS)
        self.ellipsoid_rod_colors = ti.Vector.field(3, ti.f32, n_rod * ELLIPSOID_TVERTS)
        self.ellipsoid_rod_indices = ti.field(ti.i32, n_rod * ELLIPSOID_TFACES * 3)
        self.populate_ellipsoid_rod_indices()

        # M5.23.2 arm (4) — the energy-density ISOSURFACE buffers: a budgeted
        # triangle-soup mesh (marching tetrahedra over the live energyH
        # density, engine4_render.marching_tetrahedra). iso_tri_count is the
        # RAW emission counter (atomic): a value beyond iso_max_tris means
        # overflow — the launcher warns, writes are clamped to the budget.
        # iso_level_max holds the per-frame interior density max (the GUI
        # Level slider is a fraction of it).
        self.iso_max_tris = ISO_MAX_TRIS
        self.iso_vertices = ti.Vector.field(3, ti.f32, ISO_MAX_TRIS * 3)
        self.iso_colors = ti.Vector.field(3, ti.f32, ISO_MAX_TRIS * 3)
        self.iso_indices = ti.field(ti.i32, ISO_MAX_TRIS * 3)
        self.iso_tri_count = ti.field(ti.i32, shape=())
        self.iso_level_max = ti.field(ti.f32, shape=())
        self.populate_iso_indices()

    @ti.kernel
    def populate_iso_indices(self):
        """M5.23.2 — identity index buffer for the isosurface tri soup (run
        once at init; unused slots point at collapsed origin vertices)."""
        for i in self.iso_indices:
            self.iso_indices[i] = i

    def swap_matrix_buffers(self):
        """Cycle the matrix triple buffer after evolve_M (M5.5.4): M_prev ← M, M ← M_new.

        Called once per evolve_M step. ti.template caches attribute lookups at compile
        time, so use copy_from — NOT pointer rotation (silently breaks the leapfrog).
        """
        self.M_prev_am.copy_from(self.M_am)
        self.M_am.copy_from(self.M_new_am)

    @ti.kernel
    def populate_ellipsoid_template(self):
        """VIZ.5 Stage B — fill the unit-sphere UV template (closed form per index,
        taichi-native): vertex 0 = north pole, then LAT rings × LON meridians,
        last vertex = south pole. Run once at init."""
        lat, lon = ELLIPSOID_TEMPLATE_LAT, ELLIPSOID_TEMPLATE_LON
        for v in range(ELLIPSOID_TVERTS):
            if v == 0:
                self.ellipsoid_template[v] = ti.Vector([0.0, 0.0, 1.0])
            elif v == ELLIPSOID_TVERTS - 1:
                self.ellipsoid_template[v] = ti.Vector([0.0, 0.0, -1.0])
            else:
                a = (v - 1) // lon  # ring index (0..LAT-1)
                b = (v - 1) - a * lon  # meridian index (0..LON-1)
                theta = math.pi * ti.cast(a + 1, ti.f32) / ti.cast(lat + 1, ti.f32)
                phi = 2.0 * math.pi * ti.cast(b, ti.f32) / ti.cast(lon, ti.f32)
                st = ti.sin(theta)
                self.ellipsoid_template[v] = ti.Vector(
                    [st * ti.cos(phi), st * ti.sin(phi), ti.cos(theta)]
                )

    @ti.kernel
    def populate_ellipsoid_mesh_indices(self):
        """VIZ.5 Stage B — static triangle indices for every (center, direction)
        slot of the SHELL pool: the shared `_uv_face_locals` layout, face-local
        vertex ids offset by the slot's vertex block so the whole pool renders
        as one mesh. Run once at init."""
        n_slots = ELLIPSOID_MAX_CENTERS * ELLIPSOID_MAX_DIRS
        for s, f in ti.ndrange(n_slots, ELLIPSOID_TFACES):
            loc = _uv_face_locals(f)
            base_v = s * ELLIPSOID_TVERTS
            base_i = (s * ELLIPSOID_TFACES + f) * 3
            self.ellipsoid_mesh_indices[base_i + 0] = base_v + loc[0]
            self.ellipsoid_mesh_indices[base_i + 1] = base_v + loc[1]
            self.ellipsoid_mesh_indices[base_i + 2] = base_v + loc[2]

    @ti.kernel
    def populate_ellipsoid_rod_indices(self):
        """VIZ.5 Stage D — same static index layout for the ROD pool (rod-line +
        rod-ring slots per center). Run once at init."""
        n_slots = ELLIPSOID_MAX_CENTERS * ELLIPSOID_ROD_SLOTS
        for s, f in ti.ndrange(n_slots, ELLIPSOID_TFACES):
            loc = _uv_face_locals(f)
            base_v = s * ELLIPSOID_TVERTS
            base_i = (s * ELLIPSOID_TFACES + f) * 3
            self.ellipsoid_rod_indices[base_i + 0] = base_v + loc[0]
            self.ellipsoid_rod_indices[base_i + 1] = base_v + loc[1]
            self.ellipsoid_rod_indices[base_i + 2] = base_v + loc[2]

    @ti.kernel
    def populate_grid_lines(self):
        """
        Create optimized grid lines for GGUI rendering and visualization.
        Draws continuous lines across the entire grid face-to-face.

        For a 2×2×2 voxel grid (3×3×3 grid points):
        - 9 lines parallel to X-axis (one per (j,k) pair on YZ plane)
        - 9 lines parallel to Y-axis (one per (i,k) pair on XZ plane)
        - 9 lines parallel to Z-axis (one per (i,j) pair on XY plane)
        - Total: 27 lines × 2 vertices = 54 vertices (vs 108 in non-optimized version)

        Line calculation:
        - X-direction: (ny+1) × (nz+1) lines
        - Y-direction: (nx+1) × (nz+1) lines
        - Z-direction: (nx+1) × (ny+1) lines

        Positions stored directly in normalized coordinates (0-1 range) ready for rendering.
        Uses max_grid_size for uniform normalization across asymmetric grids.
        """
        # Grid dimensions and normalization factor
        max_dim = ti.cast(self.max_grid_size, ti.f32)

        # Calculate line counts per direction
        x_lines = (self.ny + 1) * (self.nz + 1)
        y_lines = (self.nx + 1) * (self.nz + 1)
        # z_lines = (self.nx + 1) * (self.ny + 1)  # implicit, noted as reminder

        # Parallelize over all lines using single outermost loop
        for line_idx in range(self.line_count):
            vertex_idx = line_idx * 2  # Each line has 2 vertices

            if line_idx < x_lines:
                # X-parallel lines: decode (j, k) position
                temp = line_idx
                j = temp // (self.nz + 1)
                k = temp % (self.nz + 1)

                # Line from x=0 to x=nx (normalized to 0 to nx/max)
                y_norm = ti.cast(j, ti.f32) / max_dim
                z_norm = ti.cast(k, ti.f32) / max_dim

                self.grid_lines[vertex_idx] = ti.Vector([0.0, y_norm, z_norm])
                self.grid_lines[vertex_idx + 1] = ti.Vector(
                    [ti.cast(self.nx, ti.f32) / max_dim, y_norm, z_norm]
                )

            elif line_idx < x_lines + y_lines:
                # Y-parallel lines: decode (i, k) position
                temp = line_idx - x_lines
                i = temp // (self.nz + 1)
                k = temp % (self.nz + 1)

                # Line from y=0 to y=ny (normalized to 0 to ny/max)
                x_norm = ti.cast(i, ti.f32) / max_dim
                z_norm = ti.cast(k, ti.f32) / max_dim

                self.grid_lines[vertex_idx] = ti.Vector([x_norm, 0.0, z_norm])
                self.grid_lines[vertex_idx + 1] = ti.Vector(
                    [x_norm, ti.cast(self.ny, ti.f32) / max_dim, z_norm]
                )

            else:
                # Z-parallel lines: decode (i, j) position
                temp = line_idx - x_lines - y_lines
                i = temp // (self.ny + 1)
                j = temp % (self.ny + 1)

                # Line from z=0 to z=nz (normalized to 0 to nz/max)
                x_norm = ti.cast(i, ti.f32) / max_dim
                y_norm = ti.cast(j, ti.f32) / max_dim

                self.grid_lines[vertex_idx] = ti.Vector([x_norm, y_norm, 0.0])
                self.grid_lines[vertex_idx + 1] = ti.Vector(
                    [x_norm, y_norm, ti.cast(self.nz, ti.f32) / max_dim]
                )

    def populate_edge_lines(self):
        """
        Populate bounding box edge lines for GGUI wireframe rendering.

        Creates the 12 edges of the grid bounding box, each defined by 2 vertices.
        Total: 12 edges × 2 vertices = 24 vertices stored in self.edge_lines.

        Edge layout:
        - Edges 0-3: X-parallel (4 edges at YZ face corners)
        - Edges 4-7: Y-parallel (4 edges at XZ face corners)
        - Edges 8-11: Z-parallel (4 edges at XY face corners)

        Coordinates normalized by max_grid_size for uniform rendering of
        asymmetric grids (nx ≠ ny ≠ nz supported).

        Note: CPU assignment used instead of GPU kernel - for 24 fixed values,
        kernel launch overhead exceeds computation time.
        """
        # Normalized extents for asymmetric grid
        max_dim = float(self.max_grid_size)
        x_max = self.nx / max_dim
        y_max = self.ny / max_dim
        z_max = self.nz / max_dim

        # Define 12 edges as (start, end) vertex pairs
        edges = [
            # X-parallel edges (at 4 YZ corners)
            ([0, 0, 0], [x_max, 0, 0]),
            ([0, y_max, 0], [x_max, y_max, 0]),
            ([0, 0, z_max], [x_max, 0, z_max]),
            ([0, y_max, z_max], [x_max, y_max, z_max]),
            # Y-parallel edges (at 4 XZ corners)
            ([0, 0, 0], [0, y_max, 0]),
            ([x_max, 0, 0], [x_max, y_max, 0]),
            ([0, 0, z_max], [0, y_max, z_max]),
            ([x_max, 0, z_max], [x_max, y_max, z_max]),
            # Z-parallel edges (at 4 XY corners)
            ([0, 0, 0], [0, 0, z_max]),
            ([x_max, 0, 0], [x_max, 0, z_max]),
            ([0, y_max, 0], [0, y_max, z_max]),
            ([x_max, y_max, 0], [x_max, y_max, z_max]),
        ]

        # Populate field directly (no kernel needed for 24 values)
        for i, (start, end) in enumerate(edges):
            self.edge_lines[i * 2] = start
            self.edge_lines[i * 2 + 1] = end

    @ti.kernel
    def create_flux_mesh(self):
        """
        Initialize normalized flux mesh for all three orthogonal planes.

        Creates vertex positions and triangle indices for XY, XZ, and YZ flux mesh
        positioned at the xperiment parameter. Each plane is a 2D mesh that samples wave
        properties from the voxel grid.

        Coordinate system:
        - Positions stored in normalized coordinates [0, 1] for rendering
        - Each vertex corresponds to a voxel center for direct property sampling
        - Meshes intersect at xperiment parameter

        Mesh structure:
        - Vertices: Grid of 3D positions matching voxel resolution
        - Indices: Triangle pairs forming quads (2 triangles × 3 vertices = 6 indices)
        - Colors: Initialized to colormap.COLOR_FLUXMESH, updated by update_flux_mesh_values()
        """

        # Planes position in normalized coordinates
        max_dim = ti.cast(self.max_grid_size, ti.f32)
        fm_plane_x = self.flux_mesh_planes[0] * (self.nx / max_dim)  # x position normalized
        fm_plane_y = self.flux_mesh_planes[1] * (self.ny / max_dim)  # y position normalized
        fm_plane_z = self.flux_mesh_planes[2] * (self.nz / max_dim)  # z position normalized

        # ================================================================
        # XY Plane (at z = center): spans (0→1, 0→1, 0.5)
        # ================================================================
        for i, j in ti.ndrange(self.nx, self.ny):
            # Normalized coordinates for rendering
            x_norm = (ti.cast(i, ti.f32) + 0.5) / max_dim
            y_norm = (ti.cast(j, ti.f32) + 0.5) / max_dim

            # Vertex position
            self.fluxmesh_xy_vertices[i, j] = ti.Vector([x_norm, y_norm, fm_plane_z])

            # Initialize color (will be changed by update_flux_mesh_values)
            self.fluxmesh_xy_colors[i, j] = ti.Vector(colormap.COLOR_FLUXMESH[1])

        # Triangle indices for XY plane
        for i, j in ti.ndrange(self.nx - 1, self.ny - 1):
            # Each quad = 2 triangles
            # Triangle 1: (i,j) → (i+1,j) → (i,j+1)
            self.fluxmesh_xy_indices[i, j, 0] = i * self.ny + j
            self.fluxmesh_xy_indices[i, j, 1] = (i + 1) * self.ny + j
            self.fluxmesh_xy_indices[i, j, 2] = i * self.ny + (j + 1)

            # Triangle 2: (i+1,j) → (i+1,j+1) → (i,j+1)
            self.fluxmesh_xy_indices[i, j, 3] = (i + 1) * self.ny + j
            self.fluxmesh_xy_indices[i, j, 4] = (i + 1) * self.ny + (j + 1)
            self.fluxmesh_xy_indices[i, j, 5] = i * self.ny + (j + 1)

        # ================================================================
        # XZ Plane (at y = center): spans (0→1, 0.5, 0→1)
        # ================================================================
        for i, k in ti.ndrange(self.nx, self.nz):
            # Normalized coordinates for rendering
            x_norm = (ti.cast(i, ti.f32) + 0.5) / max_dim
            z_norm = (ti.cast(k, ti.f32) + 0.5) / max_dim

            # Vertex position
            self.fluxmesh_xz_vertices[i, k] = ti.Vector([x_norm, fm_plane_y, z_norm])

            # Initialize color (will be changed by update_flux_mesh_values)
            self.fluxmesh_xz_colors[i, k] = ti.Vector(colormap.COLOR_FLUXMESH[1])

        # Triangle indices for XZ plane
        for i, k in ti.ndrange(self.nx - 1, self.nz - 1):
            # Each quad = 2 triangles
            # Triangle 1: (i,k) → (i+1,k) → (i,k+1)
            self.fluxmesh_xz_indices[i, k, 0] = i * self.nz + k
            self.fluxmesh_xz_indices[i, k, 1] = (i + 1) * self.nz + k
            self.fluxmesh_xz_indices[i, k, 2] = i * self.nz + (k + 1)

            # Triangle 2: (i+1,k) → (i+1,k+1) → (i,k+1)
            self.fluxmesh_xz_indices[i, k, 3] = (i + 1) * self.nz + k
            self.fluxmesh_xz_indices[i, k, 4] = (i + 1) * self.nz + (k + 1)
            self.fluxmesh_xz_indices[i, k, 5] = i * self.nz + (k + 1)

        # ================================================================
        # YZ Plane (at x = center): spans (0.5, 0→1, 0→1)
        # ================================================================
        for j, k in ti.ndrange(self.ny, self.nz):
            # Normalized coordinates for rendering
            y_norm = (ti.cast(j, ti.f32) + 0.5) / max_dim
            z_norm = (ti.cast(k, ti.f32) + 0.5) / max_dim

            # Vertex position
            self.fluxmesh_yz_vertices[j, k] = ti.Vector([fm_plane_x, y_norm, z_norm])

            # Initialize color (will be changed by update_flux_mesh_values)
            self.fluxmesh_yz_colors[j, k] = ti.Vector(colormap.COLOR_FLUXMESH[1])

        # Triangle indices for YZ plane
        for j, k in ti.ndrange(self.ny - 1, self.nz - 1):
            # Each quad = 2 triangles
            # Triangle 1: (j,k) → (j+1,k) → (j,k+1)
            self.fluxmesh_yz_indices[j, k, 0] = j * self.nz + k
            self.fluxmesh_yz_indices[j, k, 1] = (j + 1) * self.nz + k
            self.fluxmesh_yz_indices[j, k, 2] = j * self.nz + (k + 1)

            # Triangle 2: (j+1,k) → (j+1,k+1) → (j,k+1)
            self.fluxmesh_yz_indices[j, k, 3] = (j + 1) * self.nz + k
            self.fluxmesh_yz_indices[j, k, 4] = (j + 1) * self.nz + (k + 1)
            self.fluxmesh_yz_indices[j, k, 5] = j * self.nz + (k + 1)


@ti.data_oriented
class Trackers:
    """
    Wave-property trackers for each voxel — STATEFUL, time-integrated quantities.

    Tracks amplitude envelope (EMA RMS of |ψ|) and frequency (zero-crossing
    detection) at each grid point. These observables have temporal state
    (EMA accumulator, last-crossing timestamps) so they only update meaningfully
    while dynamics are running.

    Separation of concerns (2026-05-11 refactor): per-voxel ENERGY DENSITY
    fields (`energyH_density_aJ`, `energyF_density_aJ`) and their global
    aggregates moved to `FieldObservables` below — those are *instantaneous*
    derived scalars from the current ψ (no temporal state), so they
    conceptually belong in a different class.
    """

    def __init__(self, grid_size):
        """
        Initialize tracker fields for wave-property monitoring.

        Args:
            grid_size: Grid dimensions [nx, ny, nz] matching TensorField.
        """
        # LOCAL FIELDS per voxel (stateful)
        # Amplitude tracks A via EMA of the ‖M−D‖_F deviation and RMS calculation
        # Frequency tracks local oscillation rate via zero-crossing detection
        self.amp_local_emarms_am = ti.field(dtype=ti.f32, shape=grid_size)  # am, rms amp
        self.last_crossing = ti.field(dtype=ti.f32, shape=grid_size)  # rs, last zero crossing
        self.freq_local_cross_rHz = ti.field(dtype=ti.f32, shape=grid_size)  # rHz, local frequency

        # GLOBAL AVERAGES for visualization scaling — initialized to zero; the
        # update_trackers EMA + sample_avg_trackers fill these in during sim.
        # M5 has no universal reference scale, so we let the simulation discover them.
        self.amp_global_emarms_am = ti.field(dtype=ti.f32, shape=())  # RMS all voxels
        self.freq_global_avg_rHz = ti.field(dtype=ti.f32, shape=())  # avg frequency all voxels


@ti.data_oriented
class FieldObservables:
    """
    Per-voxel derived scalar fields — STATELESS, instantaneous from current ψ.

    Each field here is a pure function of the current ψ configuration. No EMA,
    no temporal accumulator, no zero-crossing history — values are recomputed
    fresh every frame by their respective kernels. This means they are valid
    *even when dynamics are paused*: pause the sim, read the observable, get
    the current ψ's value of that observable.

    Companion to `Trackers`; the split was introduced 2026-05-11 to separate
    M3/M4-style wave statistics (stateful) from M5 field-theoretic observables
    (stateless). Each class owns its own 3-plane sampling pass for global
    aggregates (`sample_avg_trackers` vs `sample_avg_observables`).

    Future hooks (post-M5.1) — each lands as: new ti.field + dedicated kernel
    + extend `sample_avg_observables` for the global mean:
      - winding_number_local (M5.1 task 8)
      - divergence_psi (M5.2 ∇·s = 0 diagnostic)
      - curl_psi_magnitude / spin_density_s (M5.2 Close Eq. 23)
      - energy_flux (M5.3+ Poynting analog)
    """

    def __init__(self, grid_size):
        """
        Initialize derived-scalar fields.

        Args:
            grid_size: Grid dimensions [nx, ny, nz] matching TensorField.
        """
        # Per-voxel ENERGY density (HAMILTONIAN formula):
        #     H = ½|ψ̇|² + ½c²|∇ψ|² + V(ψ)
        # Naming: quantity is *energy*; the "_H" prefix denotes the formula
        # (Hamiltonian). Future formulas would parallel: `_L` for Lagrangian
        # density, `_K` for kinetic-only.
        # Populated each frame by engine3_observables.compute_energyH_density.
        # Consumed by:
        #   - sample_avg_observables  →  3-plane-sampled global mean
        #   - launcher WAVE_MENU=4 flux-mesh visualization
        # The potential term is the LdG V(M) (Eq.13) via the V_M hook in engine2_pde;
        # off (a=b=c=0) unless the xperiment sets LDG_STIFFNESS_K > 0.
        self.energyH_density_aJ = ti.field(dtype=ti.f32, shape=grid_size)  # aJ-per-voxel

        # Per-voxel FRANK ELASTIC energy density (M5.1):
        #     H_F = (K/2) · |∇n̂|²
        # Computed by engine3_observables.compute_energyF_density. Same _aJ
        # convention as energyH_density_aJ (the "_F" denotes the Frank
        # formula). The K_frank coupling is hard-coded to 1.0 (Exp 2 baseline);
        # physical value lands with M5.6 LdG elastic constants.
        # Consumed by:
        #   - launcher WAVE_MENU=5 flux-mesh visualization
        #   - M5.1 task 6 gradient-descent diagnostic (monotone decrease)
        #   - M5.1 task 7 Coulomb 1/d fit (sum over volume → E(d))
        self.energyF_density_aJ = ti.field(dtype=ti.f32, shape=grid_size)  # aJ-per-voxel

        # GLOBAL AVERAGES for visualization scaling — filled by
        # sample_avg_observables (M5.1) using its own 3-plane sampling pass.
        self.energyH_global_avg_aJ = ti.field(dtype=ti.f32, shape=())  # mean energy density (H)
        self.energyF_global_avg_aJ = ti.field(dtype=ti.f32, shape=())  # mean energy density (F)

        # EM-from-tilts observables (M5.6.5b, "see EM"): distortion modes of the
        # principal director n̂ that map onto the EM sector verified in M5.6.4 (5a §5d).
        #   director_div   = ∇·n̂  (signed SPLAY → Coulomb-charge-like; +/− at ± defects)
        #   director_curl_mag = ‖∇×n̂‖  (TWIST+BEND → B-like circulation magnitude)
        # Stateless, recomputed each frame by engine3_observables.compute_director_em.
        # Consumed by flux-mesh WAVE_MENU 6 (∇·n̂, bluered diverging) / 7 (‖∇×n̂‖, ironbow).
        self.director_div_field = ti.field(dtype=ti.f32, shape=grid_size)        # 1/am, signed
        self.director_curl_mag_field = ti.field(dtype=ti.f32, shape=grid_size)   # 1/am, ≥0
        self.director_curl_field = ti.Vector.field(3, dtype=ti.f32, shape=grid_size)  # ∇×n̂ vector (B dir)
        self.director_div_absmax = ti.field(dtype=ti.f32, shape=())   # max|∇·n̂| for symmetric color scale
        self.director_curl_max = ti.field(dtype=ti.f32, shape=())     # max‖∇×n̂‖ for color scale


if __name__ == "__main__":
    print("\n================================================================")
    print("SMOKE TEST: DATA-GRID MODULE")
    print("================================================================")

    ti.init(arch=ti.gpu)

    # ================================================================
    # Parameters & Subatomic Objects Instantiation
    # ================================================================

    UNIVERSE_SIZE = [
        2e-15,
        2e-15,
        2e-15,
    ]  # m, simulation domain [x, y, z] dimensions (can be asymmetric)

    tensor_field = TensorField(
        UNIVERSE_SIZE, target_voxels=3.5e8
    )  # 350M voxels (~14GB), 1B voxels (~40GB)

    print(f"\nGrid Statistics:")
    print(
        f"  Requested universe: [{UNIVERSE_SIZE[0]:.1e}, {UNIVERSE_SIZE[1]:.1e}, {UNIVERSE_SIZE[2]:.1e}] m"
    )
    print(
        f"  Actual universe: [{tensor_field.universe_size[0]:.1e}, {tensor_field.universe_size[1]:.1e}, {tensor_field.universe_size[2]:.1e}] m"
    )
    print(f"  Grid size: {tensor_field.nx} x {tensor_field.ny} x {tensor_field.nz} voxels")
    print(f"  Voxel edge: {tensor_field.dx:.2e} m (cubic - same for all axes)")
    print(f"  Voxel count: {tensor_field.voxel_count:,}")
    print(f"  Voxel edge (am): {tensor_field.dx_am:.2f} am")
    print(f"  Universe volume: {tensor_field.universe_volume:.2e} m³")
    print(f"  Note: voxels-per-wavelength resolution is now xperiment-driven")
    print(f"        (declared via TEST_SEED or defect Compton wavelength)")

    print("\n================================================================")
    print("END SMOKE TEST: DATA-GRID MODULE")
    print("================================================================")

    # Properly exit
    import sys

    sys.exit(0)
