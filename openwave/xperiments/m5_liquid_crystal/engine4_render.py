"""
M5 ENGINE — RENDERING / VISUALIZATION (engine4_render)

The viz stack — OpenWave's biggest differential as a simulator:
  - update_ellipsoid_mesh      VIZ.5: M·u eigen-ellipsoid shell (1 per 3D angle)
  - update_rod_ellipsoids      VIZ.5 Stage D: disclination-rod line + ring samples
  - update_director_glyphs     director-orientation line glyphs
  - update_flux_mesh_values    3-plane scalar→color/warp mesh
  - sample_position_to_render  granule point-cloud

Reads pre-computed fields off tensor_field; no engine-func calls yet (M5.4 adds
engine2_pde.eigen_decompose). Full repurposing plan in m5_4b_rendering_features.md.
"""

import math

import taichi as ti

from openwave.common import colormap

# ================================================================
# ELLIPSOID RENDERING
# ================================================================

# VIZ.5 — Fibonacci-lattice azimuth step (the golden angle). Successive shell
# directions advance by this angle while z sweeps (−1,1) uniformly → uniform S²
# coverage at ANY glyph count (the in-kernel, taichi-native direction set).
_ELLIPSOID_GOLDEN_ANGLE = math.pi * (3.0 - math.sqrt(5.0))


# VIZ.5 Stage B — isotropic visibility floor added to M before the M·u map:
# the D = diag(1, δ, 0) family has λ₃ = 0, so the raw ellipsoid degenerates to
# a flat disk; the floor lifts every semi-axis by this amount so the surface
# stays 3D-visible (the m5_6_5b λ_min floor, point 4 of that de-risk).
_ELLIPSOID_MESH_FLOOR = 0.08


@ti.func
def _ellipsoid_vertex(
    tensor_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
    v: ti.i32,  # type: ignore
    pos: ti.types.vector(3, ti.f32),  # type: ignore
    size: ti.f32,  # type: ignore
):
    """One eigen-ellipsoid surface vertex: template vertex `v` mapped through
    the spatial block of M at voxel (i,j,k), the m5_6_5b M·u key simplification
    (semi-axes = eigenvalues along eigenvectors, NO eigendecomposition), plus
    the isotropic visibility floor. Shared by the shell and rod kernels."""
    ut = tensor_field.ellipsoid_template[v]
    m = tensor_field.M_am[i, j, k]
    w = (
        ti.Vector(
            [
                m[1, 1] * ut[0] + m[1, 2] * ut[1] + m[1, 3] * ut[2],
                m[2, 1] * ut[0] + m[2, 2] * ut[1] + m[2, 3] * ut[2],
                m[3, 1] * ut[0] + m[3, 2] * ut[1] + m[3, 3] * ut[2],
            ]
        )
        + _ELLIPSOID_MESH_FLOOR * ut
    )
    return pos + 0.5 * size * w


@ti.kernel
def update_ellipsoid_mesh(
    tensor_field: ti.template(),  # type: ignore
    radius_vox: ti.f32,  # type: ignore
    size: ti.f32,  # type: ignore
    n_active: ti.i32,  # type: ignore
):
    """VIZ.5 (M5.23) — the "ellipsoid" viz (taichi MESH): one sample per 3D
    angle on the S² Fibonacci shell of radius `radius_vox` around each defect
    center; at each sample the unit-sphere template maps through the local M to
    the shaded eigen-ellipsoid surface — a FULL-3D view (not a plane
    cross-section), nearest-voxel sampling.

    The m5_6_5b key simplification: for symmetric M the image of the unit
    sphere IS the ellipsoid with semi-axes = eigenvalues along eigenvectors,
        vertex = p + (size/2) · ((M_sp + floor·I) @ u_template)
    so NO per-glyph eigendecomposition. M_sp is the spatial [1:4,1:4] block of
    the 4×4 substrate (the time axis stays out of the render); `size` matches
    the GUI Size slider (ellipsoid major diameter at vacuum λ₁ = 1).
    Colors: flat director color for now (lighting conveys the shape; a
    clock-phase palette is a later option). Slots beyond (n_centers, n_active)
    collapse to the origin (zero-area triangles, invisible).
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    max_dim = ti.cast(tensor_field.max_grid_size, ti.f32)
    n_centers = tensor_field.ellipsoid_n_centers[None]
    zero_v = ti.Vector([0.0, 0.0, 0.0])
    col = ti.Vector([_GLYPH_DIRECTOR_COLOR[0], _GLYPH_DIRECTOR_COLOR[1], _GLYPH_DIRECTOR_COLOR[2]])
    for c, d, v in ti.ndrange(
        tensor_field.ellipsoid_max_centers,
        tensor_field.ellipsoid_max_dirs,
        tensor_field.ellipsoid_tverts,
    ):
        base = (c * tensor_field.ellipsoid_max_dirs + d) * tensor_field.ellipsoid_tverts + v
        if c < n_centers and d < n_active:
            # Fibonacci-sphere shell point, direction d of n_active
            zf = 1.0 - 2.0 * (ti.cast(d, ti.f32) + 0.5) / ti.cast(n_active, ti.f32)
            rho = ti.sqrt(ti.max(1.0 - zf * zf, 0.0))
            phi = ti.cast(d, ti.f32) * _ELLIPSOID_GOLDEN_ANGLE
            u = ti.Vector([rho * ti.cos(phi), rho * ti.sin(phi), zf])
            ctr = tensor_field.ellipsoid_centers[c]
            s = ctr + radius_vox * u
            i = ti.min(ti.max(ti.cast(ti.round(s[0]), ti.i32), 0), nx - 1)
            j = ti.min(ti.max(ti.cast(ti.round(s[1]), ti.i32), 0), ny - 1)
            k = ti.min(ti.max(ti.cast(ti.round(s[2]), ti.i32), 0), nz - 1)
            pos = ti.Vector(
                [(s[0] + 0.5) / max_dim, (s[1] + 0.5) / max_dim, (s[2] + 0.5) / max_dim]
            )
            tensor_field.ellipsoid_mesh_vertices[base] = _ellipsoid_vertex(
                tensor_field, i, j, k, v, pos, size
            )
            tensor_field.ellipsoid_mesh_colors[base] = col
        else:
            tensor_field.ellipsoid_mesh_vertices[base] = zero_v
            tensor_field.ellipsoid_mesh_colors[base] = zero_v


# VIZ.5 Stage D — rod half-length as a multiple of the shell radius: the rods
# PROTRUDE beyond the shell and the rings sit on the protruding outer sections
# (matching the reference electron-clock composition; also keeps the rings
# clear of shell-ellipsoid occlusion). Consumed by the launcher dispatch.
# 2.4 so the rod tip clears the outermost ring row (2.24R, below).
ELLIPSOID_ROD_SPAN = 2.4
# Ring-row heights in SHELL-RADIUS units: the innermost row sits ONE full row
# gap off the shell (maintainer refinements 2026-07-19: gaps doubled to 0.32R,
# then the stack offset off the shell) → rows at 1.28R / 1.60R / 1.92R / 2.24R
# per pole.
_ROD_RING_H0 = 1.28
_ROD_RING_STEP = 0.32


@ti.kernel
def update_rod_ellipsoids(
    tensor_field: ti.template(),  # type: ignore
    half_len_vox: ti.f32,  # type: ignore
    ring_r_vox: ti.f32,  # type: ignore
    size: ti.f32,  # type: ignore
    n_ring_az: ti.i32,  # type: ignore
    show_rods: ti.i32,  # type: ignore
    show_rings: ti.i32,  # type: ignore
):
    """VIZ.5 (M5.23) Stage D — the disclination-ROD render. The rods are the
    line-defect pair along the spin axis (ẑ in every M5 seed): the hedgehog's
    radial director cannot stay smooth once the clock/spin structure rides it
    (the hairy-ball constraint), so the biaxial seeds construct an escaped,
    eigenvalue-melted core column, and that axis IS the angular-momentum /
    magnetic-dipole axis of the electron picture. Two arms, per (center, slot):

      - slots [0, rod_n): ROD-LINE samples along the axis, heights uniform in
        [−half_len_vox, +half_len_vox] (gated by `show_rods`; delta-cyan so the
        melted cores read apart from the shell); the degenerate / shrunken
        ellipsoid shapes ARE the rod-core melt made visible
      - remaining slots: ROD RINGS, `n_ring_az` ellipsoids per 2D angle (the
        GUI Count slider drives the azimuth density, capped by the
        ellipsoid_ring_az buffer ceiling) on circles of radius `ring_r_vox`
        AROUND the cord: FOUR rows per pole (the reference figure's count) at
        heights `_ROD_RING_H0 + _ROD_RING_STEP·row` in shell-radius units
        (1.28R / 1.60R / 1.92R / 2.24R): the whole stack offset one row gap
        OFF the shell, on the protruding outer rod sections (gated by
        `show_rings`; director light-blue) — the one-value-per-2D-angle vortex
        view placed on the actual defect line

    Same nearest-voxel sampling + `_ellipsoid_vertex` M·u map as the shell;
    inactive slots collapse to the origin (zero-area triangles, invisible).
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    max_dim = ti.cast(tensor_field.max_grid_size, ti.f32)
    n_centers = tensor_field.ellipsoid_n_centers[None]
    rod_n = tensor_field.ellipsoid_rod_n
    zero_v = ti.Vector([0.0, 0.0, 0.0])
    rod_col = ti.Vector([_GLYPH_ROD_COLOR[0], _GLYPH_ROD_COLOR[1], _GLYPH_ROD_COLOR[2]])
    ring_col = ti.Vector(
        [_GLYPH_DIRECTOR_COLOR[0], _GLYPH_DIRECTOR_COLOR[1], _GLYPH_DIRECTOR_COLOR[2]]
    )
    for c, r, v in ti.ndrange(
        tensor_field.ellipsoid_max_centers,
        tensor_field.ellipsoid_rod_slots,
        tensor_field.ellipsoid_tverts,
    ):
        base = (c * tensor_field.ellipsoid_rod_slots + r) * tensor_field.ellipsoid_tverts + v
        active = 0
        px, py, pz = 0.0, 0.0, 0.0
        col = rod_col
        if c < n_centers:
            ctr = tensor_field.ellipsoid_centers[c]
            if r < rod_n:  # rod-line sample
                if show_rods == 1:
                    active = 1
                    h = -half_len_vox + 2.0 * half_len_vox * ti.cast(r, ti.f32) / ti.cast(
                        rod_n - 1, ti.f32
                    )
                    px, py, pz = ctr[0], ctr[1], ctr[2] + h
            else:  # rod-ring sample (slot layout follows the LIVE azimuth count)
                q = r - rod_n
                ring_j = q // n_ring_az
                az = q - ring_j * n_ring_az
                if show_rings == 1 and ring_j < tensor_field.ellipsoid_ring_count:
                    active = 1
                    # 4 rows per pole at (H0 + STEP·row)·R, sign from j//4:
                    # the stack starts one row gap off the shell (R derived
                    # from the passed half-length via the span constant)
                    r_shell = half_len_vox / ELLIPSOID_ROD_SPAN
                    h = (2.0 * ti.cast(ring_j // 4, ti.f32) - 1.0) * (
                        _ROD_RING_H0 + _ROD_RING_STEP * ti.cast(ring_j % 4, ti.f32)
                    )
                    ang = 2.0 * math.pi * ti.cast(az, ti.f32) / ti.cast(n_ring_az, ti.f32)
                    px = ctr[0] + ring_r_vox * ti.cos(ang)
                    py = ctr[1] + ring_r_vox * ti.sin(ang)
                    pz = ctr[2] + h * r_shell
                    col = ring_col
        if active == 1:
            i = ti.min(ti.max(ti.cast(ti.round(px), ti.i32), 0), nx - 1)
            j = ti.min(ti.max(ti.cast(ti.round(py), ti.i32), 0), ny - 1)
            k = ti.min(ti.max(ti.cast(ti.round(pz), ti.i32), 0), nz - 1)
            pos = ti.Vector([(px + 0.5) / max_dim, (py + 0.5) / max_dim, (pz + 0.5) / max_dim])
            tensor_field.ellipsoid_rod_vertices[base] = _ellipsoid_vertex(
                tensor_field, i, j, k, v, pos, size
            )
            tensor_field.ellipsoid_rod_colors[base] = col
        else:
            tensor_field.ellipsoid_rod_vertices[base] = zero_v
            tensor_field.ellipsoid_rod_colors[base] = zero_v


# ================================================================
# GLYPH RENDERING
# ================================================================
# The flux mesh maps every voxel to ONE scalar (color + Z-warp). For a
# director field where |n̂|=1 by construction, scalar magnitude is uninformative
# (always 1 except at singular cores). To see TOPOLOGY, we need to render
# the *direction* — which is what director glyphs do.
#
# Each glyph is a line segment from voxel position to voxel + L · n̂, with
# both endpoints sharing a color computed from the signed components:
#     color = (n̂ + 1) / 2     ∈ [0, 1]³
#
# This signed-component RGB has a useful property: opposite directions get
# RGB-complementary colors (red↔cyan, green↔magenta, blue↔yellow). A +1
# hedgehog shows red on the +x face and cyan on the −x face; a −1 hedgehog
# inverts the pattern. Combined with the line offset (glyphs reach further
# in the n̂ direction), polarity is unambiguous without arrowheads.
#
# Sampling pattern: same 3-plane convention as flux_mesh — XY at fm_plane_z,
# XZ at fm_plane_y, YZ at fm_plane_x. At GLYPH_STRIDE=4 with a 64³ grid,
# total glyph count is ~768 segments; renders in well under 1 ms.
#
# Design doc: research/2b_director_glyph_rendering.md.

# ----------------------------------------------------------------
# GLYPH PARAMETERS — colors + geometry (single place; tweak freely)
# ----------------------------------------------------------------
# Per-glyph-TYPE colors. Each of the four glyph kinds gets its OWN color so they
# never read as the same thing on screen. To re-theme, swap the colormap.* source
# (or drop in a raw (r,g,b) tuple). Resolved to tuples at kernel-compile time.
#
#   - DIRECTOR / DELTA are ORIENTATION glyphs → ALWAYS use their color (they ignore
#     the Glyph-Color toggle; size/color carry no field meaning for an axis).
#   - E / B are FIELD glyphs with two color modes: "single" = the flat color below
#     (keeps weak/far-field glyphs visible — gradient palettes fade to black), and
#     "gradient" = a value-mapped palette (greenyellow for E charge, orange for B
#     magnitude) handled inline in the kernels and NOT governed by these constants.
#     The single-mode colors are kept in the SAME family as each field's gradient
#     (E green↔greenyellow, B orange↔orange) so toggling modes stays coherent, and
#     distinct from the director's light-blue so "single" no longer collides with it.
_GLYPH_DIRECTOR_COLOR = colormap.COLOR_MEDIUM[1]  # light blue — n̂ principal axis (orientation)
_GLYPH_DELTA_COLOR = colormap.COLOR_FIELD[1]  # cyan — δ cross-bar (ellipsoid minor axis)
_GLYPH_ROD_COLOR = colormap.ORANGE[1]  # orange — rod ellipsoid (line samples along the spin axis)
_GLYPH_E_COLOR = colormap.GREEN[1]  # green — E-field single-color
_GLYPH_B_COLOR = colormap.ORANGE[1]  # orange — B-field single-color

# Half-arrowhead toggle (hardcoded for testing — flip to False to compare):
# When True, each glyph gets ONE extra barb at the tip end so head vs tail is
# unambiguous without arrowhead V-shapes. Cost: doubles the line count rendered
# when on (e.g. 768 → 1536 lines at 64³ / stride=4 — still trivial on GPU).
# Buffer is always allocated; the launcher gates the second scene.lines call.
SHOW_DIRECTOR_ARROWHEAD = True
ARROWHEAD_LENGTH_FRAC = 0.35  # barb length relative to shaft length
# Barb angle measured from −n̂ (the backward shaft direction), rotated toward
# the perpendicular axis. NOTE: smaller values give THINNER arrow-head profiles
# (the variable runs opposite to a "sweep from perpendicular" intuition):
#   90° → perpendicular barb (widest, no backward sweep)
#   36° → moderately swept back (thin arrow profile — current default)
#   <30° → very thin (barb collapses toward −n̂)
ARROWHEAD_ANGLE_DEG = 36.0
_arrow_rad = math.radians(ARROWHEAD_ANGLE_DEG)
ARROWHEAD_BACK_COMP = -math.cos(_arrow_rad)  # component along n̂ (negative = backward)
ARROWHEAD_PERP_COMP = math.sin(_arrow_rad)  # component along perpendicular axis


@ti.func
def _write_glyph(
    tensor_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    idx: ti.i32,  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
    pos: ti.types.vector(3, ti.f32),  # type: ignore
    length: ti.f32,  # type: ignore
    div_scale: ti.f32,  # type: ignore
    mode: ti.i32,  # type: ignore
    show_delta: ti.i32,  # type: ignore
    size_mode: ti.i32,  # type: ignore
    color_mode: ti.i32,  # type: ignore
):
    """Write ONE centered glyph (main segment in director_glyph_*; the second
    segment — E barb OR delta cross-bar — in director_glyph_arrow_*).
    mode: 0=Director (principal axis n̂; arrow buffer = delta cross-bar if show_delta
    else blank), 1=E-field (n̂ + +→− barb). See update_director_glyphs."""
    zero_v = ti.Vector([0.0, 0.0, 0.0])
    n_hat = tensor_field.director_nhat[i, j, k]

    if mode == 0:
        # ---- Director: ORIENTATION, agnostic to size/color (always unit length +
        # the fixed _GLYPH_DIRECTOR_COLOR). The main segment is the principal axis n̂.
        # The arrow buffer carries the delta (middle-eigenvector) cross-bar WHEN
        # show_delta=1 → an ellipsoid-wireframe "+" showing the biaxial frame;
        # show_delta=0 → director axis only (arrow buffer blanked). Delta bar is
        # SHORTER (∝ λ₂/λ₁, minor:major ratio) and _GLYPH_DELTA_COLOR to set it apart
        # from the n̂ axis. (Both colors are configurable in the GLYPH PARAMETERS block.)
        dcol = ti.Vector(
            [_GLYPH_DIRECTOR_COLOR[0], _GLYPH_DIRECTOR_COLOR[1], _GLYPH_DIRECTOR_COLOR[2]]
        )
        base = pos - 0.5 * length * n_hat
        tip = pos + 0.5 * length * n_hat
        tensor_field.director_glyph_vertices[idx + 0] = base
        tensor_field.director_glyph_vertices[idx + 1] = tip
        tensor_field.director_glyph_colors[idx + 0] = dcol
        tensor_field.director_glyph_colors[idx + 1] = dcol
        if show_delta == 1:
            # delta cross-bar: middle eigenvector, length scaled by the eigenvalue ratio
            n_mid = tensor_field.director_mid[i, j, k]
            evals = tensor_field.eigenvalues[i, j, k]  # (λ₁≥λ₂≥λ₃)
            ratio = evals[1] / (ti.abs(evals[0]) + 1e-12)  # λ₂/λ₁ ∈ (0,1] → shorter bar
            dlen = length * ti.max(ti.min(ratio, 1.0), 0.12)  # clamp so it stays visible
            dlt = ti.Vector([_GLYPH_DELTA_COLOR[0], _GLYPH_DELTA_COLOR[1], _GLYPH_DELTA_COLOR[2]])
            tensor_field.director_glyph_arrow_vertices[idx + 0] = pos - 0.5 * dlen * n_mid
            tensor_field.director_glyph_arrow_vertices[idx + 1] = pos + 0.5 * dlen * n_mid
            tensor_field.director_glyph_arrow_colors[idx + 0] = dlt
            tensor_field.director_glyph_arrow_colors[idx + 1] = dlt
        else:
            tensor_field.director_glyph_arrow_vertices[idx + 0] = zero_v
            tensor_field.director_glyph_arrow_vertices[idx + 1] = zero_v
            tensor_field.director_glyph_arrow_colors[idx + 0] = zero_v
            tensor_field.director_glyph_arrow_colors[idx + 1] = zero_v
    else:
        # ---- E-field: POLAR vector (E ∝ n̂). Honors size (|∇·n̂| charge density) +
        # color (greenyellow charge gradient); centered shaft + a +→− half-barb.
        div_v = observables.director_div_field[i, j, k]
        shaft = length
        if size_mode == 1:
            shaft = length * ti.min(ti.abs(div_v) / (div_scale + 1e-12), 1.0)
        # single = flat _GLYPH_E_COLOR (green); gradient = greenyellow charge palette
        color = ti.Vector([_GLYPH_E_COLOR[0], _GLYPH_E_COLOR[1], _GLYPH_E_COLOR[2]])
        if color_mode == 1:
            color = colormap.get_greenyellow_color(div_v, -div_scale, div_scale)
        base = pos - 0.5 * shaft * n_hat
        tip = pos + 0.5 * shaft * n_hat
        tensor_field.director_glyph_vertices[idx + 0] = base
        tensor_field.director_glyph_vertices[idx + 1] = tip
        tensor_field.director_glyph_colors[idx + 0] = color
        tensor_field.director_glyph_colors[idx + 1] = color
        # +→− half-barb at the +n̂ tip (local sign gauge-arbitrary; M5.8 fixes orientation)
        ref_e = ti.Vector([0.0, 0.0, 1.0])
        if ti.abs(n_hat[2]) > 0.9:
            ref_e = ti.Vector([1.0, 0.0, 0.0])
        perp_e = n_hat.cross(ref_e).normalized()
        barb_dir = ARROWHEAD_BACK_COMP * n_hat + ARROWHEAD_PERP_COMP * perp_e
        tensor_field.director_glyph_arrow_vertices[idx + 0] = tip
        tensor_field.director_glyph_arrow_vertices[idx + 1] = (
            tip + (ARROWHEAD_LENGTH_FRAC * shaft) * barb_dir
        )
        tensor_field.director_glyph_arrow_colors[idx + 0] = color
        tensor_field.director_glyph_arrow_colors[idx + 1] = color


@ti.func
def _blank_glyph(tensor_field: ti.template(), idx: ti.i32):  # type: ignore
    """Zero both segments of glyph `idx` (off-plane → invisible)."""
    zero_v = ti.Vector([0.0, 0.0, 0.0])
    tensor_field.director_glyph_vertices[idx + 0] = zero_v
    tensor_field.director_glyph_vertices[idx + 1] = zero_v
    tensor_field.director_glyph_colors[idx + 0] = zero_v
    tensor_field.director_glyph_colors[idx + 1] = zero_v
    tensor_field.director_glyph_arrow_vertices[idx + 0] = zero_v
    tensor_field.director_glyph_arrow_vertices[idx + 1] = zero_v
    tensor_field.director_glyph_arrow_colors[idx + 0] = zero_v
    tensor_field.director_glyph_arrow_colors[idx + 1] = zero_v


@ti.kernel
def update_director_glyphs(
    tensor_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    length: ti.f32,  # type: ignore
    show_level: ti.i32,  # type: ignore
    div_scale: ti.f32,  # type: ignore
    size_mode: ti.i32,  # type: ignore
    color_mode: ti.i32,  # type: ignore
    mode: ti.i32,  # type: ignore
    show_delta: ti.i32,  # type: ignore
):
    """Centered director / E-field glyphs on the 3 flux-mesh planes (XY/XZ/YZ).

    Two of the glyph-select states (the B-curl state lives in
    `update_em_vector_glyphs`):
      - **mode=0 Director** — the principal axis `director_nhat` (main segment, full
        `length`, `_GLYPH_DIRECTOR_COLOR`). When `show_delta=1` the middle (delta)
        eigenvector `director_mid` is added as a SHORTER cross-bar (∝ λ₂/λ₁,
        `_GLYPH_DELTA_COLOR`) in the arrow buffer → the biaxial-frame "ellipsoid
        wireframe cross"; when `show_delta=0` only the n̂ axis is drawn (arrow buffer
        blanked). **Agnostic to size_mode/color_mode** — it is orientation, not a
        field, so it is always unit-length + its fixed colors. Both axes apolar
        (centered, no head/tail) ⇒ gauge-stable. The delta bar is the clock-hand axis
        (the would-be Zitterbewegung spin; in free 3D it only tilts/disperses —
        coherent spin needs M5.8 or a drive).
      - **mode=1 E-field** — `director_nhat` as a POLAR field line: honors size
        (shaft ∝ |∇·n̂| charge density) + color (single = `_GLYPH_E_COLOR`, gradient =
        greenyellow charge), with a +→− half-barb. The +→− *orientation* is
        gauge-arbitrary until the M5.8 winding-density fix (consistent with WM6
        honest-but-flipping).

    `show_level` mirrors SHOW_FLUX_MESH: 0 off, 1 XY, 2 +XZ, 3 all. Off-plane glyphs
    are zeroed (invisible 0-length segments).
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    stride = tensor_field.GLYPH_STRIDE
    max_dim = ti.cast(tensor_field.max_grid_size, ti.f32)
    nx_s, ny_s, nz_s = tensor_field.glyph_nx_s, tensor_field.glyph_ny_s, tensor_field.glyph_nz_s
    z_idx, y_idx, x_idx = (
        tensor_field.fm_plane_z_idx,
        tensor_field.fm_plane_y_idx,
        tensor_field.fm_plane_x_idx,
    )

    # ----- XY plane at z = fm_plane_z_idx -----
    for si, sj in ti.ndrange(nx_s, ny_s):
        idx = (tensor_field.glyph_offset_xy + si * ny_s + sj) * 2
        if show_level >= 1:
            i = ti.min(si * stride, nx - 1)
            j = ti.min(sj * stride, ny - 1)
            pos = ti.Vector(
                [
                    (ti.cast(i, ti.f32) + 0.5) / max_dim,
                    (ti.cast(j, ti.f32) + 0.5) / max_dim,
                    tensor_field.fm_plane_z_pos,
                ]
            )
            _write_glyph(
                tensor_field,
                observables,
                idx,
                i,
                j,
                z_idx,
                pos,
                length,
                div_scale,
                mode,
                show_delta,
                size_mode,
                color_mode,
            )
        else:
            _blank_glyph(tensor_field, idx)

    # ----- XZ plane at y = fm_plane_y_idx -----
    for si, sk in ti.ndrange(nx_s, nz_s):
        idx = (tensor_field.glyph_offset_xz + si * nz_s + sk) * 2
        if show_level >= 2:
            i = ti.min(si * stride, nx - 1)
            k = ti.min(sk * stride, nz - 1)
            pos = ti.Vector(
                [
                    (ti.cast(i, ti.f32) + 0.5) / max_dim,
                    tensor_field.fm_plane_y_pos,
                    (ti.cast(k, ti.f32) + 0.5) / max_dim,
                ]
            )
            _write_glyph(
                tensor_field,
                observables,
                idx,
                i,
                y_idx,
                k,
                pos,
                length,
                div_scale,
                mode,
                show_delta,
                size_mode,
                color_mode,
            )
        else:
            _blank_glyph(tensor_field, idx)

    # ----- YZ plane at x = fm_plane_x_idx -----
    for sj, sk in ti.ndrange(ny_s, nz_s):
        idx = (tensor_field.glyph_offset_yz + sj * nz_s + sk) * 2
        if show_level >= 3:
            j = ti.min(sj * stride, ny - 1)
            k = ti.min(sk * stride, nz - 1)
            pos = ti.Vector(
                [
                    tensor_field.fm_plane_x_pos,
                    (ti.cast(j, ti.f32) + 0.5) / max_dim,
                    (ti.cast(k, ti.f32) + 0.5) / max_dim,
                ]
            )
            _write_glyph(
                tensor_field,
                observables,
                idx,
                x_idx,
                j,
                k,
                pos,
                length,
                div_scale,
                mode,
                show_delta,
                size_mode,
                color_mode,
            )
        else:
            _blank_glyph(tensor_field, idx)


@ti.kernel
def update_em_vector_glyphs(
    tensor_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    length: ti.f32,  # type: ignore
    scale: ti.f32,  # type: ignore
    show_level: ti.i32,  # type: ignore
    size_mode: ti.i32,  # type: ignore
    color_mode: ti.i32,  # type: ignore
    curl_axis: ti.types.vector(3, ti.f32),  # type: ignore
    curl_radial: ti.i32,  # type: ignore
    curl_center: ti.types.vector(3, ti.f32),  # type: ignore
):
    """M5.6.5b — B-direction glyphs: half-arrow segments along ∇×n̂ (the curl/circulation vector).

    Reuses the director-glyph BUFFERS (shaft + arrowhead) + 3-plane sampling; the segment points
    along the curl vector observables.director_curl_field. Color modes (color_mode):
      - **single (0)** — flat `_GLYPH_B_COLOR` (orange), keeps weak/far-field glyphs visible.
      - **gradient (1)** — the SIGNED BLUERED N/S projection (`_curl_signed_proj` → `get_bluered`),
        the SAME radial/axial scalar the WM7 flux-mesh N/S coloring uses (radial `(∇×n̂)·r̂` about
        `curl_center` when `curl_radial=1` → red=N/blue=S poles; axial `(∇×n̂)·curl_axis` otherwise).
        So a gradient-colored B glyph reads its *pole* (sign), not just magnitude — matching the
        mesh under it (cf. the E glyph's greenyellow ±charge gradient).
    Shaft length ∝ min(‖curl‖/scale, 1) so the view declutters where there is no circulation
    (static charge ⇒ ~zero-length, invisible; real twist ⇒ visible arrows with a half-barb tip
    showing the B direction). The barb scales with the shaft. `scale` = the shared distortion
    magnitude max(|∇·n̂|, ‖∇×n̂‖) (matches WAVE_MENU 7). show_level mirrors SHOW_DIRECTORS
    (0 off, 1 XY, 2 +XZ, 3 all). Writes director_glyph_vertices/colors + the arrow buffers.
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    stride = tensor_field.GLYPH_STRIDE
    max_dim = ti.cast(tensor_field.max_grid_size, ti.f32)
    nx_s, ny_s, nz_s = tensor_field.glyph_nx_s, tensor_field.glyph_ny_s, tensor_field.glyph_nz_s
    z_idx, y_idx, x_idx = (
        tensor_field.fm_plane_z_idx,
        tensor_field.fm_plane_y_idx,
        tensor_field.fm_plane_x_idx,
    )
    inv_s = 1.0 / (scale + 1e-12)
    zero_v = ti.Vector([0.0, 0.0, 0.0])

    for plane in range(3):
        n_si = nx_s if plane < 2 else ny_s  # plane 0=XY, 1=XZ, 2=YZ
        n_sj = ny_s if plane == 0 else nz_s
        base = tensor_field.glyph_offset_xy
        if plane == 1:
            base = tensor_field.glyph_offset_xz
        if plane == 2:
            base = tensor_field.glyph_offset_yz
        for sa, sb in ti.ndrange(n_si, n_sj):
            idx = (base + sa * n_sj + sb) * 2
            on = (
                (plane == 0 and show_level >= 1)
                or (plane == 1 and show_level >= 2)
                or (plane == 2 and show_level >= 3)
            )
            if on:
                # map (sa, sb) → voxel (i, j, k) on the active plane
                i = ti.min(sa * stride, nx - 1) if plane < 2 else x_idx
                j = y_idx
                k = z_idx
                if plane == 0:
                    j = ti.min(sb * stride, ny - 1)
                if plane == 1:
                    k = ti.min(sb * stride, nz - 1)
                if plane == 2:
                    j = ti.min(sa * stride, ny - 1)
                    k = ti.min(sb * stride, nz - 1)
                curl = observables.director_curl_field[i, j, k]
                mag = curl.norm()
                dirv = curl / (mag + 1e-12)
                # +0.5 on the two in-plane axes (flux-mesh cell center); the
                # perpendicular axis is then snapped to the mesh's continuous plane
                # coord (fm_plane_*_pos) so glyphs sit exactly on the sheet.
                half = ti.Vector([0.5, 0.5, 0.5])
                if plane == 0:  # XY slice → k is perpendicular
                    half[2] = 0.0
                elif plane == 1:  # XZ slice → j is perpendicular
                    half[1] = 0.0
                else:  # YZ slice → i is perpendicular
                    half[0] = 0.0
                pos = (
                    ti.Vector([ti.cast(i, ti.f32), ti.cast(j, ti.f32), ti.cast(k, ti.f32)]) + half
                ) / max_dim
                if plane == 0:
                    pos[2] = tensor_field.fm_plane_z_pos
                elif plane == 1:
                    pos[1] = tensor_field.fm_plane_y_pos
                else:
                    pos[0] = tensor_field.fm_plane_x_pos
                # shaft: unit (dirv≈unit wherever curl detectable) or field magnitude (declutters)
                shaft = length
                if size_mode == 1:
                    shaft = length * ti.min(mag * inv_s, 1.0)
                # CENTERED on the voxel (base = pos − ½·shaft·dirv → tip = pos + ½·shaft·dirv),
                # the same field-line-through-the-point convention as the director + E glyphs
                # (VIZ.1, 4b §4.2). The +→ barb at the tip still marks the circulation direction.
                half_seg = 0.5 * shaft * dirv
                tip = pos + half_seg
                color = ti.Vector(
                    [_GLYPH_B_COLOR[0], _GLYPH_B_COLOR[1], _GLYPH_B_COLOR[2]]
                )  # single = flat _GLYPH_B_COLOR (orange), keeps far-field glyphs visible
                if color_mode == 1:
                    # gradient = SIGNED bluered N/S — same radial/axial projection AND the same
                    # γ-compression-against-own-max as the WM7 mesh, so glyph and mesh match
                    # (red=N/blue=S, visible far from the core). `scale` (shared) drives the SIZE
                    # declutter above; the COLOR uses B's own director_curl_max.
                    proj = _curl_signed_proj(curl, i, j, k, curl_axis, curl_radial, curl_center)
                    gsig = _gamma_signed(proj, observables.director_curl_max[None])
                    color = colormap.get_bluered_color(gsig, -1.0, 1.0)
                tensor_field.director_glyph_vertices[idx + 0] = pos - half_seg
                tensor_field.director_glyph_vertices[idx + 1] = tip
                tensor_field.director_glyph_colors[idx + 0] = color
                tensor_field.director_glyph_colors[idx + 1] = color
                # half-arrow barb at the tip (stable perp; /(norm+eps) avoids NaN at curl≈0)
                ref = ti.Vector([0.0, 0.0, 1.0])
                if ti.abs(dirv[2]) > 0.9:
                    ref = ti.Vector([1.0, 0.0, 0.0])
                perp = dirv.cross(ref)
                perp = perp / (perp.norm() + 1e-12)
                barb_dir = ARROWHEAD_BACK_COMP * dirv + ARROWHEAD_PERP_COMP * perp
                barb_end = tip + (ARROWHEAD_LENGTH_FRAC * shaft) * barb_dir
                tensor_field.director_glyph_arrow_vertices[idx + 0] = tip
                tensor_field.director_glyph_arrow_vertices[idx + 1] = barb_end
                tensor_field.director_glyph_arrow_colors[idx + 0] = color
                tensor_field.director_glyph_arrow_colors[idx + 1] = color
            else:
                tensor_field.director_glyph_vertices[idx + 0] = zero_v
                tensor_field.director_glyph_vertices[idx + 1] = zero_v
                tensor_field.director_glyph_colors[idx + 0] = zero_v
                tensor_field.director_glyph_colors[idx + 1] = zero_v
                tensor_field.director_glyph_arrow_vertices[idx + 0] = zero_v
                tensor_field.director_glyph_arrow_vertices[idx + 1] = zero_v
                tensor_field.director_glyph_arrow_colors[idx + 0] = zero_v
                tensor_field.director_glyph_arrow_colors[idx + 1] = zero_v


@ti.kernel
def update_moment_glyph(
    tensor_field: ti.template(),  # type: ignore
    m_axis: ti.types.vector(3, ti.f32),  # type: ignore
    cx: ti.f32,  # type: ignore
    cy: ti.f32,  # type: ignore
    cz: ti.f32,  # type: ignore
    length: ti.f32,  # type: ignore
    color: ti.types.vector(3, ti.f32),  # type: ignore
):
    """VIZ.4 — a single magnetic-MOMENT vector glyph `μ` at the defect center,
    pointing along m̂ (POLAR → centered shaft + a half-barb at the +m̂ tip).

    A static marker for the dipole-sample placeholder: it labels the moment axis
    the bluered N/S field is organized around. Writes 4 vertices into
    moment_glyph_{vertices,colors} (shaft base→tip, barb tip→back); the launcher
    renders them with scene.lines when DIPOLE_SAMPLE is active. Center (cx,cy,cz)
    in voxel coords; `length` in normalized [0,1] (larger than voxel glyphs so it
    reads as the principal axis)."""
    max_dim = ti.cast(tensor_field.max_grid_size, ti.f32)
    mhat = m_axis.normalized()
    center = (ti.Vector([cx, cy, cz]) + 0.5) / max_dim
    base = center - 0.5 * length * mhat
    tip = center + 0.5 * length * mhat
    tensor_field.moment_glyph_vertices[0] = base
    tensor_field.moment_glyph_vertices[1] = tip
    # half-barb at the +m̂ tip (stable perpendicular reference)
    ref = ti.Vector([0.0, 0.0, 1.0])
    if ti.abs(mhat[2]) > 0.9:
        ref = ti.Vector([1.0, 0.0, 0.0])
    perp = mhat.cross(ref).normalized()
    barb_dir = ARROWHEAD_BACK_COMP * mhat + ARROWHEAD_PERP_COMP * perp
    tensor_field.moment_glyph_vertices[2] = tip
    tensor_field.moment_glyph_vertices[3] = tip + (ARROWHEAD_LENGTH_FRAC * length) * barb_dir
    for n in range(4):
        tensor_field.moment_glyph_colors[n] = color


# ================================================================
# FLUX MESH VALUES UPDATING
# ================================================================
# The flux mesh is a VISUALIZATION layer — it converts simulation-side
# scalars/vectors to per-vertex colors and Z-axis warps so the user can "see"
# what the field is doing. It is not physics; nothing here feeds back into
# evolve_M or any tracker. Treat it as a display driver.
#
# VECTOR-FIELD RENDERING NUANCE (worth remembering when M5.0g+ adds new
# wave_menu options for energy density, curl, divergence, energy flux,
# etc.):
#
# The M substrate yields vector observables (the director n̂, ∇·n̂, ∇×n̂, force
# vectors), but the flux mesh maps every voxel to ONE scalar (one color, one
# Z-warp height). A magnitude projection |v| = v.norm() rectifies a signed field
# — e.g. for a y-polarized v_y = A·sin(k·x):
#
#     |v| = A·|sin(k·x)|         period λ/2 — looks like 2× as many bumps
#     v_y = A· sin(k·x)          period λ — the "true" sinusoid (signed)
#
# Magnitude is the right scalar for energy-density-style observables (always
# positive, ∝ |v|²) but rectifies signed waves and visually doubles their
# spatial frequency. When designing future wave_menu modes, decide per
# observable:
#
#   - magnitude: |v|              → energy/intensity views (always ≥ 0)
#   - signed component: v · ê     → propagation/polarization views (signed)
#   - radial component: v · r̂    → longitudinal vs transverse decomposition
#   - dominant axis: max(|v_i|)·sign  → quick-and-dirty signed view
#
# For the flux mesh specifically (one scalar per vertex), each new wave_menu
# entry should pick one explicitly and not silently fall back to .norm().
# This becomes important once we have curl(n̂), force vectors, energy flux,
# etc. — all of which are vector quantities that need a deliberate scalar
# projection to be rendered.


@ti.func
def _curl_signed_proj(
    curl_vec: ti.types.vector(3, ti.f32),  # type: ignore
    vi: ti.i32,  # type: ignore
    vj: ti.i32,  # type: ignore
    vk: ti.i32,  # type: ignore
    curl_axis: ti.types.vector(3, ti.f32),  # type: ignore
    curl_radial: ti.i32,  # type: ignore
    curl_center: ti.types.vector(3, ti.f32),  # type: ignore
):
    """Signed scalar for the bluered N/S coloring of ∇×n̂ (B). Two projections:

    - **radial** (`curl_radial=1`): `(∇×n̂)·r̂` with `r̂` from `curl_center` (voxel
      coords) → the TRUE magnetic N/S poles — red where B flows OUT (N hemisphere,
      cosθ>0), blue where it flows IN (S hemisphere) → `∝ cosθ`, matching a bar
      magnet (Duda's N-red-top / S-blue-bottom). Needs a defined center.
    - **axial** (`curl_radial=0`): `(∇×n̂)·curl_axis` against a FIXED axis → the
      field's axial COMPONENT. For an ideal dipole this is red at BOTH poles
      (B ∥ m̂ along the whole axis) + blue at the equator — physically real, but
      reads as two red lobes, not a bar magnet. The default for general WM7 runs
      where no single center exists (M5.8 wires radial to the real defect center)."""
    proj = curl_vec.dot(curl_axis)
    if curl_radial == 1:
        r = ti.Vector(
            [
                ti.cast(vi, ti.f32) - curl_center[0],
                ti.cast(vj, ti.f32) - curl_center[1],
                ti.cast(vk, ti.f32) - curl_center[2],
            ]
        )
        rhat = r / (r.norm() + 1e-12)
        proj = curl_vec.dot(rhat)
    return proj


# Bluered signed-color magnitude compression (γ < 1). The B-field (∇×n̂) of a dipole falls
# off as 1/r³ — ~9× steeper than the E charge's 1/r — so a LINEAR map dumps everything beyond
# the core into black (only a thin shell of glyphs/mesh near the poles stays visible). γ<1
# stretches small/mid magnitudes toward the palette extremes so a wide area reads in visible
# color; the peak still clips to the extreme, 0 stays black-center. 1.0 = linear (original).
# 0.4 spreads the dipole from the core out to ~r=24 voxels (probe 2026-05-31). Tune freely.
_BLUERED_GAMMA = 0.4


@ti.func
def _gamma_signed(value: ti.f32, scale: ti.f32) -> ti.f32:  # type: ignore
    """SIGNED, γ-compressed, normalized magnitude ∈ [−1, 1] for bluered N/S coloring.

    sign → red(+) / blue(−); `|value/scale|^γ` lifts mid/far values toward ±1 (the palette
    extremes) while the peak clips to ±1 and 0 stays 0 (black-center). `scale` should be the
    field's OWN max (e.g. `director_curl_max` for B), NOT a shared max — otherwise the peak
    never reaches the extreme. Returns a SCALAR (feed it to `get_bluered_color(_, -1, 1)`),
    shared by the B glyph and the WM7 mesh so their N/S coloring matches."""
    t = ti.math.clamp(value / (scale + 1e-12), -1.0, 1.0)
    sgn = 1.0 - 2.0 * ti.cast(t < 0.0, ti.f32)  # branchless sign: +1 if t≥0 else −1
    return sgn * ti.pow(ti.abs(t), _BLUERED_GAMMA)


@ti.kernel
def update_flux_mesh_values(
    tensor_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    wave_menu: ti.i32,  # type: ignore
    warp_mesh: ti.i32,  # type: ignore
    curl_color: ti.i32,  # type: ignore
    curl_axis: ti.types.vector(3, ti.f32),  # type: ignore
    curl_radial: ti.i32,  # type: ignore
    curl_center: ti.types.vector(3, ti.f32),  # type: ignore
):
    """
    Update flux mesh colors and vertices by sampling wave properties from voxel grid.

    Samples wave displacement at each plane vertex position and maps it to a color.
    Should be called every frame after wave propagation to update visualization.

    Args:
        tensor_field: TensorField instance with flux mesh fields and M-substrate data
        trackers: Trackers — reads amp/freq (WAVE_MENU 2 / 3 + colormap range)
        observables: FieldObservables — reads energyH/energyF (WAVE_MENU 4 / 5)
        wave_menu: Selected Wave displayed with color palette
    """

    # ================================================================
    # XY Plane: Sample at z = fm_plane_z_idx
    # ================================================================
    # Always update all planes (conditionals cause GPU branch divergence)
    # wave_menu == 4 renders the energy density (Hamiltonian) field
    # `observables.energyH_density_aJ`, populated each step by M5.0g.
    # wave_menu == 5 renders the Frank elastic density `observables.energyF_density_aJ`,
    # populated by compute_energyF_density (M5.1 task 5).
    for i, j in ti.ndrange(tensor_field.nx, tensor_field.ny):
        # Sample director (M5.4: orientation, not displacement) + matrix observables
        dir_value = tensor_field.director_nhat[i, j, tensor_field.fm_plane_z_idx]
        amp_value = trackers.amp_local_emarms_am[i, j, tensor_field.fm_plane_z_idx]
        freq_value = trackers.freq_local_cross_rHz[i, j, tensor_field.fm_plane_z_idx]
        energyH_value = observables.energyH_density_aJ[i, j, tensor_field.fm_plane_z_idx]
        energyF_value = observables.energyF_density_aJ[i, j, tensor_field.fm_plane_z_idx]
        univ_edge_z = tensor_field.universe_size_am[2]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 1:  # Orientation deviation ‖n̂−ẑ‖ on orange (wave_menu == 1)
            dev = (dir_value - ti.Vector([0.0, 0.0, 1.0])).norm()  # 0 at vacuum ẑ, →2 at −ẑ
            tensor_field.fluxmesh_xy_colors[i, j] = colormap.get_orange_color(dev, 0.0, 2.0)
            tensor_field.fluxmesh_xy_vertices[i, j][2] = (
                dev / 2.0 * 0.3 * warp_mesh / 300.0
                + tensor_field.flux_mesh_planes[2] * (tensor_field.nz / tensor_field.max_grid_size)
            )
        elif wave_menu == 2:  # ironbow
            tensor_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                amp_value, 0, trackers.amp_global_emarms_am[None] * 2
            )
            tensor_field.fluxmesh_xy_vertices[i, j][2] = (
                amp_value / univ_edge_z * warp_mesh
                + tensor_field.flux_mesh_planes[2] * (tensor_field.nz / tensor_field.max_grid_size)
            )
        elif wave_menu == 3:  # blueprint
            tensor_field.fluxmesh_xy_colors[i, j] = colormap.get_blueprint_color(
                freq_value, 0.0, trackers.freq_global_avg_rHz[None] * 2
            )
            tensor_field.fluxmesh_xy_vertices[i, j][2] = freq_value / trackers.freq_global_avg_rHz[
                None
            ] / 3000 * warp_mesh + tensor_field.flux_mesh_planes[2] * (
                tensor_field.nz / tensor_field.max_grid_size
            )
        elif wave_menu == 4:  # Energy density (Hamiltonian) on ironbow
            H_max = observables.energyH_global_avg_aJ[None] * 4.0 + 1e-10
            tensor_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                energyH_value, 0.0, H_max
            )
            tensor_field.fluxmesh_xy_vertices[i, j][2] = (
                energyH_value / H_max * 0.3 * warp_mesh / 300.0
                + (
                    tensor_field.flux_mesh_planes[2]
                    * (tensor_field.nz / tensor_field.max_grid_size)
                )
            )
        elif wave_menu == 5:  # Frank elastic density on ironbow (defect-focused palette)
            F_max = observables.energyF_global_avg_aJ[None] * 4.0 + 1e-10
            tensor_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                energyF_value, 0.0, F_max
            )
            tensor_field.fluxmesh_xy_vertices[i, j][2] = (
                energyF_value / F_max * 0.3 * warp_mesh / 300.0
                + (
                    tensor_field.flux_mesh_planes[2]
                    * (tensor_field.nz / tensor_field.max_grid_size)
                )
            )
        elif (
            wave_menu == 6
        ):  # EM divergence ∇·n̂ (splay = Coulomb charge) on greenyellow diverging
            div_s = observables.director_div_absmax[None] + 1e-12
            div_v = observables.director_div_field[i, j, tensor_field.fm_plane_z_idx]
            tensor_field.fluxmesh_xy_colors[i, j] = colormap.get_greenyellow_color(
                div_v, -div_s, div_s
            )
            tensor_field.fluxmesh_xy_vertices[i, j][2] = (
                div_v / div_s * 0.3 * warp_mesh / 300.0
                + tensor_field.flux_mesh_planes[2] * (tensor_field.nz / tensor_field.max_grid_size)
            )
        elif wave_menu == 7:  # EM curl ‖∇×n̂‖ (twist+bend = circulation/B) on orange
            # Shared director-distortion scale (M5.6.5b): scale curl against max(|∇·n̂|, ‖∇×n̂‖),
            # NOT its own max — so a static hedgehog's tiny discretization-noise curl reads
            # near-black (no circulation = no B for a static charge); real twist/dynamics lights it.
            curl_s = (
                ti.max(observables.director_div_absmax[None], observables.director_curl_max[None])
                + 1e-12
            )
            curl_v = observables.director_curl_mag_field[i, j, tensor_field.fm_plane_z_idx]
            curl_vec = observables.director_curl_field[i, j, tensor_field.fm_plane_z_idx]
            # VIZ.2/VIZ.4 color toggle: 0=orange magnitude (honest static default),
            # 1=bluered signed — radial (∇×n̂)·r̂ = true N/S poles when curl_radial,
            # else axial (∇×n̂)·curl_axis. See _curl_signed_proj.
            if curl_color == 1:
                # bluered N/S: γ-compressed against B's OWN max (director_curl_max), so the
                # steep 1/r³ dipole stays visible far from the core (not curl_s, which is the
                # shared E-inclusive max used for the warp/size declutter below).
                tensor_field.fluxmesh_xy_colors[i, j] = colormap.get_bluered_color(
                    _gamma_signed(
                        _curl_signed_proj(
                            curl_vec,
                            i,
                            j,
                            tensor_field.fm_plane_z_idx,
                            curl_axis,
                            curl_radial,
                            curl_center,
                        ),
                        observables.director_curl_max[None],
                    ),
                    -1.0,
                    1.0,
                )
            else:
                tensor_field.fluxmesh_xy_colors[i, j] = colormap.get_orange_color(
                    curl_v, 0.0, curl_s
                )
            # VIZ.2 vector-warp: displace the vertex by the RAW ∇×n̂ vector (all 3
            # components, scaled), so the mesh deforms as a *twist in fabric* showing
            # the B-field rotation + handedness — not just a perpendicular magnitude lift.
            warp_amt = 0.3 * warp_mesh / 300.0 / curl_s
            tensor_field.fluxmesh_xy_vertices[i, j] = ti.Vector(
                [
                    (ti.cast(i, ti.f32) + 0.5) / tensor_field.max_grid_size
                    + curl_vec[0] * warp_amt,
                    (ti.cast(j, ti.f32) + 0.5) / tensor_field.max_grid_size
                    + curl_vec[1] * warp_amt,
                    tensor_field.flux_mesh_planes[2]
                    * (tensor_field.nz / tensor_field.max_grid_size)
                    + curl_vec[2] * warp_amt,
                ]
            )

    # ================================================================
    # XZ Plane: Sample at y = fm_plane_y_idx
    # ================================================================
    for i, k in ti.ndrange(tensor_field.nx, tensor_field.nz):
        # Sample director (M5.4: orientation, not displacement) + matrix observables
        dir_value = tensor_field.director_nhat[i, tensor_field.fm_plane_y_idx, k]
        amp_value = trackers.amp_local_emarms_am[i, tensor_field.fm_plane_y_idx, k]
        freq_value = trackers.freq_local_cross_rHz[i, tensor_field.fm_plane_y_idx, k]
        energyH_value = observables.energyH_density_aJ[i, tensor_field.fm_plane_y_idx, k]
        energyF_value = observables.energyF_density_aJ[i, tensor_field.fm_plane_y_idx, k]
        univ_edge_y = tensor_field.universe_size_am[1]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 1:  # Orientation deviation ‖n̂−ẑ‖ on orange (wave_menu == 1)
            dev = (dir_value - ti.Vector([0.0, 0.0, 1.0])).norm()  # 0 at vacuum ẑ, →2 at −ẑ
            tensor_field.fluxmesh_xz_colors[i, k] = colormap.get_orange_color(dev, 0.0, 2.0)
            tensor_field.fluxmesh_xz_vertices[i, k][1] = (
                dev / 2.0 * 0.3 * warp_mesh / 300.0
                + tensor_field.flux_mesh_planes[1] * (tensor_field.ny / tensor_field.max_grid_size)
            )
        elif wave_menu == 2:  # ironbow
            tensor_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                amp_value, 0, trackers.amp_global_emarms_am[None] * 2
            )
            tensor_field.fluxmesh_xz_vertices[i, k][1] = (
                amp_value / univ_edge_y * warp_mesh
                + tensor_field.flux_mesh_planes[1] * (tensor_field.ny / tensor_field.max_grid_size)
            )
        elif wave_menu == 3:  # blueprint
            tensor_field.fluxmesh_xz_colors[i, k] = colormap.get_blueprint_color(
                freq_value, 0.0, trackers.freq_global_avg_rHz[None] * 2
            )
            tensor_field.fluxmesh_xz_vertices[i, k][1] = freq_value / trackers.freq_global_avg_rHz[
                None
            ] / 3000 * warp_mesh + tensor_field.flux_mesh_planes[1] * (
                tensor_field.ny / tensor_field.max_grid_size
            )
        elif wave_menu == 4:  # Energy density (Hamiltonian) on ironbow
            H_max = observables.energyH_global_avg_aJ[None] * 4.0 + 1e-10
            tensor_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                energyH_value, 0.0, H_max
            )
            tensor_field.fluxmesh_xz_vertices[i, k][1] = (
                energyH_value / H_max * 0.3 * warp_mesh / 300.0
                + (
                    tensor_field.flux_mesh_planes[1]
                    * (tensor_field.ny / tensor_field.max_grid_size)
                )
            )
        elif wave_menu == 5:  # Frank elastic density on ironbow
            F_max = observables.energyF_global_avg_aJ[None] * 4.0 + 1e-10
            tensor_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                energyF_value, 0.0, F_max
            )
            tensor_field.fluxmesh_xz_vertices[i, k][1] = (
                energyF_value / F_max * 0.3 * warp_mesh / 300.0
                + (
                    tensor_field.flux_mesh_planes[1]
                    * (tensor_field.ny / tensor_field.max_grid_size)
                )
            )
        elif (
            wave_menu == 6
        ):  # EM divergence ∇·n̂ (splay = Coulomb charge) on greenyellow diverging
            div_s = observables.director_div_absmax[None] + 1e-12
            div_v = observables.director_div_field[i, tensor_field.fm_plane_y_idx, k]
            tensor_field.fluxmesh_xz_colors[i, k] = colormap.get_greenyellow_color(
                div_v, -div_s, div_s
            )
            tensor_field.fluxmesh_xz_vertices[i, k][1] = (
                div_v / div_s * 0.3 * warp_mesh / 300.0
                + tensor_field.flux_mesh_planes[1] * (tensor_field.ny / tensor_field.max_grid_size)
            )
        elif wave_menu == 7:  # EM curl ‖∇×n̂‖ (twist+bend = circulation/B) on orange
            # Shared director-distortion scale (M5.6.5b): scale curl against max(|∇·n̂|, ‖∇×n̂‖),
            # NOT its own max — so a static hedgehog's tiny discretization-noise curl reads
            # near-black (no circulation = no B for a static charge); real twist/dynamics lights it.
            curl_s = (
                ti.max(observables.director_div_absmax[None], observables.director_curl_max[None])
                + 1e-12
            )
            curl_v = observables.director_curl_mag_field[i, tensor_field.fm_plane_y_idx, k]
            curl_vec = observables.director_curl_field[i, tensor_field.fm_plane_y_idx, k]
            if curl_color == 1:  # bluered signed — radial r̂ (N/S poles) or axial
                # γ-compressed against B's own max (see XY plane note)
                tensor_field.fluxmesh_xz_colors[i, k] = colormap.get_bluered_color(
                    _gamma_signed(
                        _curl_signed_proj(
                            curl_vec,
                            i,
                            tensor_field.fm_plane_y_idx,
                            k,
                            curl_axis,
                            curl_radial,
                            curl_center,
                        ),
                        observables.director_curl_max[None],
                    ),
                    -1.0,
                    1.0,
                )
            else:
                tensor_field.fluxmesh_xz_colors[i, k] = colormap.get_orange_color(
                    curl_v, 0.0, curl_s
                )
            # VIZ.2 vector-warp by raw ∇×n̂ (fabric-twist)
            warp_amt = 0.3 * warp_mesh / 300.0 / curl_s
            tensor_field.fluxmesh_xz_vertices[i, k] = ti.Vector(
                [
                    (ti.cast(i, ti.f32) + 0.5) / tensor_field.max_grid_size
                    + curl_vec[0] * warp_amt,
                    tensor_field.flux_mesh_planes[1]
                    * (tensor_field.ny / tensor_field.max_grid_size)
                    + curl_vec[1] * warp_amt,
                    (ti.cast(k, ti.f32) + 0.5) / tensor_field.max_grid_size
                    + curl_vec[2] * warp_amt,
                ]
            )

    # ================================================================
    # YZ Plane: Sample at x = fm_plane_x_idx
    # ================================================================
    for j, k in ti.ndrange(tensor_field.ny, tensor_field.nz):
        # Sample director (M5.4: orientation, not displacement) + matrix observables
        dir_value = tensor_field.director_nhat[tensor_field.fm_plane_x_idx, j, k]
        amp_value = trackers.amp_local_emarms_am[tensor_field.fm_plane_x_idx, j, k]
        freq_value = trackers.freq_local_cross_rHz[tensor_field.fm_plane_x_idx, j, k]
        energyH_value = observables.energyH_density_aJ[tensor_field.fm_plane_x_idx, j, k]
        energyF_value = observables.energyF_density_aJ[tensor_field.fm_plane_x_idx, j, k]
        univ_edge_x = tensor_field.universe_size_am[0]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 1:  # Orientation deviation ‖n̂−ẑ‖ on orange (wave_menu == 1)
            dev = (dir_value - ti.Vector([0.0, 0.0, 1.0])).norm()  # 0 at vacuum ẑ, →2 at −ẑ
            tensor_field.fluxmesh_yz_colors[j, k] = colormap.get_orange_color(dev, 0.0, 2.0)
            tensor_field.fluxmesh_yz_vertices[j, k][0] = (
                dev / 2.0 * 0.3 * warp_mesh / 300.0
                + tensor_field.flux_mesh_planes[0] * (tensor_field.nx / tensor_field.max_grid_size)
            )
        elif wave_menu == 2:  # ironbow
            tensor_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                amp_value, 0, trackers.amp_global_emarms_am[None] * 2
            )
            tensor_field.fluxmesh_yz_vertices[j, k][0] = (
                amp_value / univ_edge_x * warp_mesh
                + tensor_field.flux_mesh_planes[0] * (tensor_field.nx / tensor_field.max_grid_size)
            )
        elif wave_menu == 3:  # blueprint
            tensor_field.fluxmesh_yz_colors[j, k] = colormap.get_blueprint_color(
                freq_value, 0.0, trackers.freq_global_avg_rHz[None] * 2
            )
            tensor_field.fluxmesh_yz_vertices[j, k][0] = freq_value / trackers.freq_global_avg_rHz[
                None
            ] / 3000 * warp_mesh + tensor_field.flux_mesh_planes[0] * (
                tensor_field.nx / tensor_field.max_grid_size
            )
        elif wave_menu == 4:  # Energy density (Hamiltonian) on ironbow
            H_max = observables.energyH_global_avg_aJ[None] * 4.0 + 1e-10
            tensor_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                energyH_value, 0.0, H_max
            )
            tensor_field.fluxmesh_yz_vertices[j, k][0] = (
                energyH_value / H_max * 0.3 * warp_mesh / 300.0
                + (
                    tensor_field.flux_mesh_planes[0]
                    * (tensor_field.nx / tensor_field.max_grid_size)
                )
            )
        elif wave_menu == 5:  # Frank elastic density on ironbow
            F_max = observables.energyF_global_avg_aJ[None] * 4.0 + 1e-10
            tensor_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                energyF_value, 0.0, F_max
            )
            tensor_field.fluxmesh_yz_vertices[j, k][0] = (
                energyF_value / F_max * 0.3 * warp_mesh / 300.0
                + (
                    tensor_field.flux_mesh_planes[0]
                    * (tensor_field.nx / tensor_field.max_grid_size)
                )
            )
        elif (
            wave_menu == 6
        ):  # EM divergence ∇·n̂ (splay = Coulomb charge) on greenyellow diverging
            div_s = observables.director_div_absmax[None] + 1e-12
            div_v = observables.director_div_field[tensor_field.fm_plane_x_idx, j, k]
            tensor_field.fluxmesh_yz_colors[j, k] = colormap.get_greenyellow_color(
                div_v, -div_s, div_s
            )
            tensor_field.fluxmesh_yz_vertices[j, k][0] = (
                div_v / div_s * 0.3 * warp_mesh / 300.0
                + tensor_field.flux_mesh_planes[0] * (tensor_field.nx / tensor_field.max_grid_size)
            )
        elif wave_menu == 7:  # EM curl ‖∇×n̂‖ (twist+bend = circulation/B) on orange
            # Shared director-distortion scale (M5.6.5b): scale curl against max(|∇·n̂|, ‖∇×n̂‖),
            # NOT its own max — so a static hedgehog's tiny discretization-noise curl reads
            # near-black (no circulation = no B for a static charge); real twist/dynamics lights it.
            curl_s = (
                ti.max(observables.director_div_absmax[None], observables.director_curl_max[None])
                + 1e-12
            )
            curl_v = observables.director_curl_mag_field[tensor_field.fm_plane_x_idx, j, k]
            curl_vec = observables.director_curl_field[tensor_field.fm_plane_x_idx, j, k]
            if curl_color == 1:  # bluered signed — radial r̂ (N/S poles) or axial
                # γ-compressed against B's own max (see XY plane note)
                tensor_field.fluxmesh_yz_colors[j, k] = colormap.get_bluered_color(
                    _gamma_signed(
                        _curl_signed_proj(
                            curl_vec,
                            tensor_field.fm_plane_x_idx,
                            j,
                            k,
                            curl_axis,
                            curl_radial,
                            curl_center,
                        ),
                        observables.director_curl_max[None],
                    ),
                    -1.0,
                    1.0,
                )
            else:
                tensor_field.fluxmesh_yz_colors[j, k] = colormap.get_orange_color(
                    curl_v, 0.0, curl_s
                )
            # VIZ.2 vector-warp by raw ∇×n̂ (fabric-twist)
            warp_amt = 0.3 * warp_mesh / 300.0 / curl_s
            tensor_field.fluxmesh_yz_vertices[j, k] = ti.Vector(
                [
                    tensor_field.flux_mesh_planes[0]
                    * (tensor_field.nx / tensor_field.max_grid_size)
                    + curl_vec[0] * warp_amt,
                    (ti.cast(j, ti.f32) + 0.5) / tensor_field.max_grid_size
                    + curl_vec[1] * warp_amt,
                    (ti.cast(k, ti.f32) + 0.5) / tensor_field.max_grid_size
                    + curl_vec[2] * warp_amt,
                ]
            )


# ================================================================
# GRANULE POSITION RENDER
# ================================================================


@ti.kernel
def sample_position_to_render(
    tensor_field: ti.template(),  # type: ignore
    amp_boost: ti.f32,  # type: ignore
    stride: ti.i32,  # type: ignore
    num_render: ti.i32,  # type: ignore
):
    """Sample granule positions from z flux_mesh plane with stride, writing to 1D position_render.

    Samples every `stride`-th voxel from the XY plane at the z flux mesh index,
    capping output at `num_render` particles for performance.

    M5.4 — DIRECTOR POINT-CLOUD (interim): the granule sphere sits at voxel + amp_boost·n̂,
    where n̂ = director_nhat (principal eigenvector of M). The old "granule pushed by the
    displacement wave" has no meaning on an orientation field. This is the cheap 4b interim;
    the full biaxial-ellipsoid showcase is deferred to M5.6 (in uniaxial M5.4 the two minor
    eigen-axes are degenerate, so only a surface-of-revolution reads correctly — premature
    until δ≠g). amp_boost = WARP_MESH (dimensionless director offset; no /dx_am — n̂ is unit).
    """
    k = int(tensor_field.flux_mesh_planes[2] * tensor_field.grid_size[2])
    max_dim = ti.cast(tensor_field.max_grid_size, ti.f32)
    sampled_ny = (tensor_field.ny + stride - 1) // stride  # cols per sampled row

    for render_idx in range(num_render):
        si = render_idx // sampled_ny  # sampled row index
        sj = render_idx % sampled_ny  # sampled col index
        i = si * stride
        j = sj * stride
        displaced = amp_boost * tensor_field.director_nhat[i, j, k] + ti.Vector(
            [ti.cast(i, ti.f32), ti.cast(j, ti.f32), ti.cast(k, ti.f32)]
        )
        tensor_field.position_render[render_idx] = displaced / max_dim


# ================================================================
# M5.23.2 arm (4) — ENERGY-DENSITY ISOSURFACE (marching tetrahedra)
# ================================================================
# Extracts the level surface of the live Hamiltonian density
# (observables.energyH_density_aJ, computed per-frame by the launcher's
# compute_field_observables — on the canonical path V4 + curvature are
# EXACTLY zero on the covariant vacuum, so level sets are meaningful with
# no background subtraction, the M5.24 true-zero property).
#
# METHOD: marching TETRAHEDRA on the 6-tet decomposition of each grid cell
# around the c0-c6 diagonal — chosen over classic marching cubes for a
# TABLE-FREE kernel (the 256-case tri table is replaced by 16 enumerable
# per-tet cases, auditable by reading; the cost is ~2x triangles for the
# same surface, covered by the ISO_MAX_TRIS budget + overflow warning).
# Triangles are emitted as independent soup via an atomic counter; the
# render draws two_sided so winding order is not load-bearing (same rule
# as the rod pool). Vertices land in the shared render space
# (voxel + 0.5)/max_grid_size.

_ISO_COLOR = (1.0, 0.84, 0.30)  # gold — energy surface
_ISO_EPS = 1e-12


@ti.kernel
def iso_density_max(
    tensor_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
):
    """Density max over the interior EXCLUDING a 3-voxel boundary margin →
    tensor_field.iso_level_max (the GUI Level slider is a fraction of this
    per-frame scale). The margin matters on LOADED research states: the
    fixed-J endpoint's density max sits at the disclination-rod / pin-shell
    boundary junction (measured 4.2x the core peak at the m5_23_2 build) —
    normalizing on the grid max would compress the slider's useful range to
    the artifact. The surface itself still marches the full valid region
    (boundary blobs render, clipped at the region edge like any marcher)."""
    tensor_field.iso_level_max[None] = 0.0
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for i, j, k in ti.ndrange((3, nx - 3), (3, ny - 3), (3, nz - 3)):
        ti.atomic_max(tensor_field.iso_level_max[None], observables.energyH_density_aJ[i, j, k])


@ti.func
def _iso_lerp(pa, pb, va: ti.f32, vb: ti.f32, lv: ti.f32):  # type: ignore
    """Edge intersection of the level lv on segment (pa, va)-(pb, vb)."""
    t = 0.5
    if ti.abs(vb - va) > _ISO_EPS:
        t = (lv - va) / (vb - va)
    t = ti.min(ti.max(t, 0.0), 1.0)
    return pa + t * (pb - pa)


@ti.func
def _iso_emit(tensor_field: ti.template(), pa, pb, pc):  # type: ignore
    """Append one triangle to the soup (budget-clamped; the raw counter
    keeps counting so the launcher can report overflow honestly)."""
    idx = ti.atomic_add(tensor_field.iso_tri_count[None], 1)
    if idx < tensor_field.iso_max_tris:
        col = ti.Vector([_ISO_COLOR[0], _ISO_COLOR[1], _ISO_COLOR[2]])
        tensor_field.iso_vertices[3 * idx + 0] = pa
        tensor_field.iso_vertices[3 * idx + 1] = pb
        tensor_field.iso_vertices[3 * idx + 2] = pc
        tensor_field.iso_colors[3 * idx + 0] = col
        tensor_field.iso_colors[3 * idx + 1] = col
        tensor_field.iso_colors[3 * idx + 2] = col


@ti.func
def _iso_tet(tensor_field: ti.template(), p0, p1, p2, p3, v0: ti.f32, v1: ti.f32, v2: ti.f32, v3: ti.f32, lv: ti.f32):  # type: ignore
    """Emit the level-surface triangles of one tetrahedron (16 cases:
    0/4 inside → nothing, 1/3 inside → one triangle, 2 inside → two)."""
    b0, b1, b2, b3 = v0 > lv, v1 > lv, v2 > lv, v3 > lv
    n_in = int(b0) + int(b1) + int(b2) + int(b3)
    if n_in == 1 or n_in == 3:
        # the LONE vertex (inside if n_in == 1, outside if n_in == 3)
        lone_is_in = n_in == 1
        la, va = p0, v0
        ob, vb = p1, v1
        oc, vc = p2, v2
        od, vd = p3, v3
        if b1 == lone_is_in and b0 != lone_is_in:
            la, va = p1, v1
            ob, vb = p0, v0
        elif b2 == lone_is_in and b0 != lone_is_in and b1 != lone_is_in:
            la, va = p2, v2
            oc, vc = p0, v0
        elif b3 == lone_is_in and b0 != lone_is_in and b1 != lone_is_in and b2 != lone_is_in:
            la, va = p3, v3
            od, vd = p0, v0
        qa = _iso_lerp(la, ob, va, vb, lv)
        qb = _iso_lerp(la, oc, va, vc, lv)
        qc = _iso_lerp(la, od, va, vd, lv)
        _iso_emit(tensor_field, qa, qb, qc)
    elif n_in == 2:
        # the inside PAIR (a, b) vs the outside pair (c, d): 6 combinations
        pa, va = p0, v0
        pb, vb = p1, v1
        pc, vc = p2, v2
        pd, vd = p3, v3
        if b0 and b1:
            pass  # (0,1) in / (2,3) out — the default assignment
        elif b0 and b2:
            pb, vb = p2, v2
            pc, vc = p1, v1
        elif b0 and b3:
            pb, vb = p3, v3
            pd, vd = p1, v1
        elif b1 and b2:
            pa, va = p1, v1
            pb, vb = p2, v2
            pc, vc = p0, v0
        elif b1 and b3:
            pa, va = p1, v1
            pb, vb = p3, v3
            pd, vd = p0, v0
        else:  # b2 and b3
            pa, va = p2, v2
            pb, vb = p3, v3
            pc, vc = p0, v0
            pd, vd = p1, v1
        qa = _iso_lerp(pa, pc, va, vc, lv)
        qb = _iso_lerp(pa, pd, va, vd, lv)
        qc = _iso_lerp(pb, pd, vb, vd, lv)
        qd = _iso_lerp(pb, pc, vb, vc, lv)
        _iso_emit(tensor_field, qa, qb, qc)
        _iso_emit(tensor_field, qa, qc, qd)


@ti.kernel
def marching_tetrahedra(
    tensor_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    level: ti.f32,  # type: ignore
):
    """Fill the iso tri-soup buffers with the level surface of the energy
    density. Cells span the density's valid interior [1, n-2] (the boundary
    ring is zero by kernel convention). Collapses unused slots."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    max_dim = ti.cast(tensor_field.max_grid_size, ti.f32)
    tensor_field.iso_tri_count[None] = 0
    for i, j, k in ti.ndrange((1, nx - 2), (1, ny - 2), (1, nz - 2)):
        # cube corners (standard order): c0..c7
        v000 = observables.energyH_density_aJ[i, j, k]
        v100 = observables.energyH_density_aJ[i + 1, j, k]
        v110 = observables.energyH_density_aJ[i + 1, j + 1, k]
        v010 = observables.energyH_density_aJ[i, j + 1, k]
        v001 = observables.energyH_density_aJ[i, j, k + 1]
        v101 = observables.energyH_density_aJ[i + 1, j, k + 1]
        v111 = observables.energyH_density_aJ[i + 1, j + 1, k + 1]
        v011 = observables.energyH_density_aJ[i, j + 1, k + 1]
        # early out: cell entirely on one side
        vmin = ti.min(ti.min(ti.min(v000, v100), ti.min(v110, v010)), ti.min(ti.min(v001, v101), ti.min(v111, v011)))
        vmax = ti.max(ti.max(ti.max(v000, v100), ti.max(v110, v010)), ti.max(ti.max(v001, v101), ti.max(v111, v011)))
        if vmin <= level and level < vmax:
            x0 = (ti.cast(i, ti.f32) + 0.5) / max_dim
            y0 = (ti.cast(j, ti.f32) + 0.5) / max_dim
            z0 = (ti.cast(k, ti.f32) + 0.5) / max_dim
            d = 1.0 / max_dim
            p000 = ti.Vector([x0, y0, z0])
            p100 = ti.Vector([x0 + d, y0, z0])
            p110 = ti.Vector([x0 + d, y0 + d, z0])
            p010 = ti.Vector([x0, y0 + d, z0])
            p001 = ti.Vector([x0, y0, z0 + d])
            p101 = ti.Vector([x0 + d, y0, z0 + d])
            p111 = ti.Vector([x0 + d, y0 + d, z0 + d])
            p011 = ti.Vector([x0, y0 + d, z0 + d])
            # 6-tet decomposition around the c000-c111 diagonal
            _iso_tet(tensor_field, p000, p101, p100, p111, v000, v101, v100, v111, level)
            _iso_tet(tensor_field, p000, p100, p110, p111, v000, v100, v110, v111, level)
            _iso_tet(tensor_field, p000, p110, p010, p111, v000, v110, v010, v111, level)
            _iso_tet(tensor_field, p000, p010, p011, p111, v000, v010, v011, v111, level)
            _iso_tet(tensor_field, p000, p011, p001, p111, v000, v011, v001, v111, level)
            _iso_tet(tensor_field, p000, p001, p101, p111, v000, v001, v101, v111, level)


@ti.kernel
def collapse_iso_tail(tensor_field: ti.template()):  # type: ignore
    """Zero the vertex slots beyond the emitted count (degenerate, invisible).
    Run after marching_tetrahedra; count is budget-clamped here."""
    n_used = ti.min(tensor_field.iso_tri_count[None], tensor_field.iso_max_tris)
    zero_v = ti.Vector([0.0, 0.0, 0.0])
    for t in range(tensor_field.iso_max_tris):
        if t >= n_used:
            tensor_field.iso_vertices[3 * t + 0] = zero_v
            tensor_field.iso_vertices[3 * t + 1] = zero_v
            tensor_field.iso_vertices[3 * t + 2] = zero_v


@ti.kernel
def update_iso_ellipsoids(
    tensor_field: ti.template(),  # type: ignore
    size: ti.f32,  # type: ignore
    n_cover: ti.i32,  # type: ignore
):
    """The author's "surface uniformly covered with ellipsoids" variant:
    the VIZ.5 M·u ellipsoid samples placed at (strided) isosurface triangle
    centroids instead of the S² shell — writes the SHELL pool buffers
    (exclusive with the shell view; the launcher gates them). n_cover =
    live sample budget (the GUI Count slider × centers, slot-clamped)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    max_dim = ti.cast(tensor_field.max_grid_size, ti.f32)
    n_tri = ti.min(tensor_field.iso_tri_count[None], tensor_field.iso_max_tris)
    n_slots = tensor_field.ellipsoid_max_centers * tensor_field.ellipsoid_max_dirs
    n_act = ti.min(ti.min(n_cover, n_slots), n_tri)
    stride = 1
    if n_act > 0:
        stride = ti.max(n_tri // n_act, 1)
    col = ti.Vector([_ISO_COLOR[0], _ISO_COLOR[1], _ISO_COLOR[2]])
    zero_v = ti.Vector([0.0, 0.0, 0.0])
    for s, v in ti.ndrange(n_slots, tensor_field.ellipsoid_tverts):
        base = s * tensor_field.ellipsoid_tverts + v
        t = s * stride
        if s < n_act and t < n_tri:
            ctr = (
                tensor_field.iso_vertices[3 * t + 0]
                + tensor_field.iso_vertices[3 * t + 1]
                + tensor_field.iso_vertices[3 * t + 2]
            ) / 3.0
            i = ti.min(ti.max(ti.cast(ti.round(ctr[0] * max_dim - 0.5), ti.i32), 0), nx - 1)
            j = ti.min(ti.max(ti.cast(ti.round(ctr[1] * max_dim - 0.5), ti.i32), 0), ny - 1)
            k = ti.min(ti.max(ti.cast(ti.round(ctr[2] * max_dim - 0.5), ti.i32), 0), nz - 1)
            tensor_field.ellipsoid_mesh_vertices[base] = _ellipsoid_vertex(
                tensor_field, i, j, k, v, ctr, size
            )
            tensor_field.ellipsoid_mesh_colors[base] = col
        else:
            tensor_field.ellipsoid_mesh_vertices[base] = zero_v
            tensor_field.ellipsoid_mesh_colors[base] = zero_v
