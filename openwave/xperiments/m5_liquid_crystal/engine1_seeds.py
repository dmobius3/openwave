"""
M5 ENGINE — INITIAL-CONDITION SEEDING (engine1_seeds)

Seed kernels that prime the matrix triple-buffer (M_prev_am, M_am, M_new_am)
with a known order-parameter configuration. Split from the former monolithic
lagrangian_engine.py (2026-05-25 SoC refactor). Concern order:
  engine1_seeds → engine2_pde →
  engine3_observables → engine4_render
"""

import taichi as ti

# ================================================================
# MATRIX-FIELD SEEDERS (M5.4 — order parameter M = O·D·O^T)
# ================================================================
# Matrix-substrate analogs of seed_vacuum / seed_hedgehog. The director n̂ is
# computed exactly as in the ψ versions (same blend, same soft-core), then
# embedded in the uniaxial order parameter
#     M = δ·I + (1−δ)·n̂⊗n̂        (eigenvalues 1, δ, δ; principal eigvec = n̂)
# Each seeder writes BOTH the matrix triple buffer (M_prev/M/M_new — triple-buffer
# BC rule, so the matrix leapfrog landing in M5.5 starts BC-consistent) AND the
# derived director_nhat with the correct PHYSICAL sign. Seeding director_nhat
# directly keeps eigen_decompose's sign-continuity branch smooth from t=0 and lets
# the M5.4 director-equivalent Coulomb relaxation operate on a clean director.


@ti.func
def embed4(msp, g):  # type: ignore
    """M5.8.1 — embed a 3×3 spatial order parameter into the 4×4 substrate as
    block-diag(g, M_spatial): time axis = index 0 (eigenvalue g; Duda D=diag(g,1,δ,0)),
    boost-DECOUPLED (the time row/col off-diagonals are 0). The spatial block [1:4,1:4]
    is the M5.4/M5.6 order parameter, extracted by callers for the 3×3-Cardano solver."""
    m4 = ti.Matrix.zero(ti.f32, 4, 4)
    for a in ti.static(range(3)):
        for bb in ti.static(range(3)):
            m4[a + 1, bb + 1] = msp[a, bb]
    m4[0, 0] = g
    return m4


@ti.func
def uniaxial_M(n, delta, g):  # type: ignore
    """Uniaxial order parameter M = δ·I + (1−δ)·n̂⊗n̂ from a unit director n̂, embedded
    in the 4×4 substrate (M5.8.1: block-diag(spatial, g)).

    Spatial eigenvalues (1, δ, δ), principal eigenvector = n̂; the 4th (time) axis = g.
    Verified by the eigen_decompose round-trip (spectrum recovers (1, δ, δ); director
    recovers n̂ at 0.9995).
    """
    eye = ti.Matrix.identity(ti.f32, 3)
    msp = delta * eye + (1.0 - delta) * n.outer_product(n)
    return embed4(msp, g)


@ti.kernel
def seed_vacuum_M(
    tensor_field: ti.template(),  # type: ignore
    delta: ti.f32,  # type: ignore
):
    """Fill the matrix field with the topological-vacuum ground state (n̂ = ẑ).

    M_vac = δ·I + (1−δ)·ẑ⊗ẑ everywhere (uniform orientation = no defects). Writes
    all three matrix buffers + director_nhat (= ẑ). The matrix vacuum seeder.

    Args:
        tensor_field: TensorField (writes M_am, M_prev_am, M_new_am, director_nhat)
        delta: uniaxial minor-axis eigenvalue (tensor_field.lc_delta)
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    z_hat = ti.Vector([0.0, 0.0, 1.0])
    m_vac = uniaxial_M(z_hat, delta, tensor_field.lc_g)
    for i, j, k in ti.ndrange(nx, ny, nz):
        tensor_field.M_am[i, j, k] = m_vac
        tensor_field.M_prev_am[i, j, k] = m_vac
        tensor_field.M_new_am[i, j, k] = m_vac
        tensor_field.director_nhat[i, j, k] = z_hat
        tensor_field.director_nhat_new[i, j, k] = z_hat


@ti.kernel
def seed_hedgehog_M(
    tensor_field: ti.template(),  # type: ignore
    centers_voxel: ti.template(),  # type: ignore
    signs: ti.template(),  # type: ignore
    domain_quarter_voxels: ti.f32,  # type: ignore
    n_defects: ti.i32,  # type: ignore
    delta: ti.f32,  # type: ignore
):
    """Seed N hedgehog defects on the matrix substrate.

    Director construction: weighted radial superposition + vacuum blend +
    renormalization, then embeds the unit director via uniaxial_M and writes the matrix
    triple buffer + director_nhat. The multi-defect API drives the M5.4 Coulomb-pair
    gate directly (centers ±d/2 along x, signs ±1).

    Args:
        tensor_field: TensorField (writes M_am, M_prev_am, M_new_am, director_nhat)
        centers_voxel: ti.field (n_defects, 3) i32 — centers in voxel coords
        signs: ti.field (n_defects,) i32 — +1 outward / −1 inward per defect
        domain_quarter_voxels: f32 — D/4 in voxel units (w_vac falloff radius)
        n_defects: i32 — number of active defects
        delta: uniaxial minor-axis eigenvalue (tensor_field.lc_delta)
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    z_hat = ti.Vector([0.0, 0.0, 1.0])
    soft_core_voxels = ti.cast(0.2, ti.f32)
    proximity_floor = ti.cast(0.5, ti.f32)

    for i, j, k in ti.ndrange(nx, ny, nz):
        x = ti.cast(i, ti.f32)
        y = ti.cast(j, ti.f32)
        z = ti.cast(k, ti.f32)

        n_combined = ti.Vector([0.0, 0.0, 0.0])
        r_nearest = ti.cast(1e10, ti.f32)

        for d in range(n_defects):
            cx = ti.cast(centers_voxel[d, 0], ti.f32)
            cy = ti.cast(centers_voxel[d, 1], ti.f32)
            cz = ti.cast(centers_voxel[d, 2], ti.f32)
            sgn = ti.cast(signs[d], ti.f32)

            ddx = x - cx
            ddy = y - cy
            ddz = z - cz
            r = ti.sqrt(ddx * ddx + ddy * ddy + ddz * ddz)
            r_safe = ti.max(r, soft_core_voxels)

            n_radial = ti.Vector(
                [sgn * ddx / r_safe, sgn * ddy / r_safe, sgn * ddz / r_safe]
            )
            w_defect = 1.0 / (r + proximity_floor)
            n_combined += w_defect * n_radial
            r_nearest = ti.min(r_nearest, r)

        ratio = r_nearest / domain_quarter_voxels
        w_vac = 1.0 / (1.0 + ratio * ratio * ratio * ratio)
        n_blended = w_vac * n_combined + (1.0 - w_vac) * z_hat
        n_unit = n_blended / (n_blended.norm() + 1e-12)

        m = uniaxial_M(n_unit, delta, tensor_field.lc_g)
        tensor_field.M_am[i, j, k] = m
        tensor_field.M_prev_am[i, j, k] = m
        tensor_field.M_new_am[i, j, k] = m
        tensor_field.director_nhat[i, j, k] = n_unit
        tensor_field.director_nhat_new[i, j, k] = n_unit


@ti.kernel
def seed_biaxial_hedgehog_M(
    tensor_field: ti.template(),  # type: ignore
    cx: ti.i32,  # type: ignore
    cy: ti.i32,  # type: ignore
    cz: ti.i32,  # type: ignore
    r0_vox: ti.f32,  # type: ignore
    rhoc_vox: ti.f32,  # type: ignore
    delta: ti.f32,  # type: ignore
):
    """Seed a single BIAXIAL hedgehog: M = O·D(s(r))·Oᵀ, D = diag(1, δ, 0) (M5.6.5a).

    The production port of the sandbox biaxial hedgehog (M5.6.2a frame + M5.6.3b melt,
    `5a §5b/§5c`). Unlike the uniaxial `seed_hedgehog_M` (eigenvalues (1,δ,δ)), this gives
    THREE distinct eigenvalues — the genuine biaxial substrate of the eigenvalue map
    (1 = EM tilt axis, δ ~ ℏ = QM twist axis, 0 = null axis):

      frame O = [r̂ | e_Θ | e_Φ]:  r̂ = principal (EM/1) axis = unit radial;
        e_Φ = azimuthal (the z-axis DISCLINATION), melted by a clamped smoothstep over ρ_c
        (secondary axes shrink inside the disclination core — biaxiality melts like a nematic);
        e_Θ = e_Φ × r̂.
      radial eigenvalue melt s(r) = r/√(r²+r0²):  D → isotropic at the core (singularity
        regularized, Faber `‖q⃗‖→0`), → diag(1,δ,0) far out.   M = Σ_a d_a · col_a ⊗ col_a.

    Writes all three M buffers + director_nhat (= r̂) per the triple-buffer BC rule.

    Args:
        cx, cy, cz: defect center in voxel coords
        r0_vox: radial eigenvalue-melt scale (core size, voxels)
        rhoc_vox: z-axis disclination-melt scale (voxels)
        delta: middle eigenvalue δ (D = diag(1, δ, 0))
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    d_iso = (1.0 + delta) / 3.0                          # isotropic value, trace-preserving
    eps = ti.cast(1e-6, ti.f32)
    for i, j, k in ti.ndrange(nx, ny, nz):
        px = ti.cast(i - cx, ti.f32)
        py = ti.cast(j - cy, ti.f32)
        pz = ti.cast(k - cz, ti.f32)
        r2 = px * px + py * py + pz * pz
        r = ti.sqrt(r2)
        rho = ti.sqrt(px * px + py * py)

        rinv = 1.0 / ti.max(r, eps)
        rhat = ti.Vector([px * rinv, py * rinv, pz * rinv])           # unit radial = EM/1 axis
        ainv = 1.0 / ti.sqrt(rho * rho + eps)
        azim = ti.Vector([-py * ainv, px * ainv, 0.0])                # azimuthal (z-disclination)
        sdisc = ti.min(rho / rhoc_vox, 1.0)
        shrink = sdisc * sdisc * (3.0 - 2.0 * sdisc)                  # clamped smoothstep melt
        ephi = azim * shrink
        ephi = ephi - ephi.dot(rhat) * rhat                          # ⟂ r̂ (no renorm → melts)
        etheta = ephi.cross(rhat)

        srad = r / ti.sqrt(r2 + r0_vox * r0_vox)                     # radial melt (Faber sin α)
        d0 = d_iso + srad * (1.0 - d_iso)
        d1 = d_iso + srad * (delta - d_iso)
        d2 = d_iso + srad * (0.0 - d_iso)
        m_sp = (d0 * rhat.outer_product(rhat)
                + d1 * etheta.outer_product(etheta)
                + d2 * ephi.outer_product(ephi))
        m = embed4(m_sp, tensor_field.lc_g)  # M5.8.1: 4×4 block-diag(spatial, g)
        tensor_field.M_am[i, j, k] = m
        tensor_field.M_prev_am[i, j, k] = m
        tensor_field.M_new_am[i, j, k] = m
        tensor_field.director_nhat[i, j, k] = rhat
        tensor_field.director_nhat_new[i, j, k] = rhat


@ti.kernel
def seed_charged_ring_M(
    tensor_field: ti.template(),  # type: ignore
    cx: ti.i32,  # type: ignore
    cy: ti.i32,  # type: ignore
    cz: ti.i32,  # type: ignore
    a_vox: ti.f32,  # type: ignore
    cord_vox: ti.f32,  # type: ignore
    rhoc_vox: ti.f32,  # type: ignore
    delta: ti.f32,  # type: ignore
):
    """Seed the CHARGED DISCLINATION RING — the hedgehog's same-charge partner
    (M5.21.2, third defect type; Alexander, Chen, Matsumoto & Kamien,
    Rev. Mod. Phys. 84, 497 (2012) § IV.B: disclination loops carry odd
    hedgehog charge — the far sphere reads q = 1 exactly like the point
    hedgehog, but the core is a half-disclination cord of radius a in the
    z = 0 plane and the interior is smoothly ESCAPED along +z).

    Meridional director angle (single-valued in the tensor through the
    half winding; branch cuts rotated so the DIRECTOR seam lies on the
    minimal spanning disk, see below):
        ψ(ρ, z) = ½·[atan2(−z, ρ − a) + atan2(−z, ρ + a)] + π/2
        n̂ = sin ψ · ρ̂ + cos ψ · ẑ     (→ r̂ far away = the hedgehog far
                                        field; → ±ẑ inside the ring)

    ⚠️ DIRECTOR SEAM (representation, not physics): a half-integer
    disclination ring is NON-ORIENTABLE around the cord, so any vector
    director field must flip sign across SOME surface bounded by the
    ring (the Dirac sheet). The atan2 cuts above place that seam on the
    MINIMAL SPANNING DISK inside the ring (ρ < a, z = 0): the whole
    exterior (far field, below-plane region) is seam-free. The tensor
    M is continuous everywhere; only director-DERIVED diagnostics (EM
    div/curl in engine3, div-colored glyphs) show the internal disk.
    (First cut of this kernel used atan2(ρ−a, z), whose seam is a
    half-infinite cylinder at ρ = a, z < 0 — the 2026-07-17 GGUI
    below-plane artifact. Fixed same day.)
    Frame (production biaxial convention, δ on e_Θ, 0 on the azimuth —
    NOTE the M5.21.2 sandbox census used the rotated family, δ on the
    azimuth; both are seeds of the same winding sector):
        e_Φ = azimuth (melted over ρ_c near the z-axis, ⟂ n̂),
        e_Θ = e_Φ × n̂,  D = (1, δ, 0) on (n̂, e_Θ, e_Φ).
    Eigenvalue melt: isotropic at the ring CORD (distance
    d = √((ρ−a)² + z²), smoothstep over cord_vox) — the cord is the only
    singular locus; the origin/axis are regular (escaped interior).

    Static viewing seed: spatial 3×3 embedded block-diag with g (embed4),
    boost-decoupled, no clock tangent. Triple-buffer BC rule.
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    d_iso = (1.0 + delta) / 3.0
    eps = ti.cast(1e-6, ti.f32)
    for i, j, k in ti.ndrange(nx, ny, nz):
        px = ti.cast(i - cx, ti.f32)
        py = ti.cast(j - cy, ti.f32)
        pz = ti.cast(k - cz, ti.f32)
        rho = ti.sqrt(px * px + py * py)
        ainv = 1.0 / ti.sqrt(rho * rho + eps)
        rhohat = ti.Vector([px * ainv, py * ainv, 0.0])
        # seam-on-disk branch choice (see docstring): cuts at z=0, rho<a
        psi = (0.5 * (ti.atan2(-pz, rho - a_vox) + ti.atan2(-pz, rho + a_vox))
               + 1.5707963267948966)
        sp = ti.sin(psi)
        cp = ti.cos(psi)
        nvec = ti.Vector([sp * rhohat[0], sp * rhohat[1], cp])
        azim = ti.Vector([-py * ainv, px * ainv, 0.0])
        sdisc = ti.min(rho / rhoc_vox, 1.0)
        shrink = sdisc * sdisc * (3.0 - 2.0 * sdisc)      # axis melt (house style)
        ephi = azim * shrink
        ephi = ephi - ephi.dot(nvec) * nvec               # ⟂ n̂ (no renorm → melts)
        etheta = ephi.cross(nvec)
        dc = ti.sqrt((rho - a_vox) * (rho - a_vox) + pz * pz)
        sc = ti.min(dc / cord_vox, 1.0)
        smelt = sc * sc * (3.0 - 2.0 * sc)                # cord melt smoothstep
        d0 = d_iso + smelt * (1.0 - d_iso)
        d1 = d_iso + smelt * (delta - d_iso)
        d2 = d_iso + smelt * (0.0 - d_iso)
        m_sp = (d0 * nvec.outer_product(nvec)
                + d1 * etheta.outer_product(etheta)
                + d2 * ephi.outer_product(ephi))
        m = embed4(m_sp, tensor_field.lc_g)
        tensor_field.M_am[i, j, k] = m
        tensor_field.M_prev_am[i, j, k] = m
        tensor_field.M_new_am[i, j, k] = m
        tensor_field.director_nhat[i, j, k] = nvec
        tensor_field.director_nhat_new[i, j, k] = nvec


@ti.kernel
def seed_dressed_hedgehog_M(
    tensor_field: ti.template(),  # type: ignore
    cx: ti.i32,  # type: ignore
    cy: ti.i32,  # type: ignore
    cz: ti.i32,  # type: ignore
    r0_vox: ti.f32,  # type: ignore
    rhoc_vox: ti.f32,  # type: ignore
    delta: ti.f32,  # type: ignore
    b_star: ti.f32,  # type: ignore
    rw_vox: ti.f32,  # type: ignore
    kick_theta: ti.f32,  # type: ignore
):
    """Seed the M5.8.2 BOOST-DRESSED biaxial hedgehog: M = W·D₄·Wᵀ, W = O₄·B(θ).

    The 2b-1 ground state ported to production: the biaxial hedgehog frame + melt
    (identical to seed_biaxial_hedgehog_M) DRESSED with a core-localized boost
    B(θ) mixing the e_Θ eigen-axis (index 2) with the time axis (index 0),
    θ(r) = b*·exp(−(r/r_w)²). At the ground dressing (b* ≈ 0.13) this LOWERS the
    static energy below the bare defect (the GEM dip, 5a §10e winning recipe).
    ⚠️ r̂ stays the spatial principal axis only for modest b* (the boost feeds
    sinh²·g into the e_Θ spatial eigenvalue: 0.3 → ~0.44 at b*=0.13, still < 1).

    Also writes M_psi_am — the clock tangent ∂M/∂Θ for the (δ,0)-plane rotation
    (eigen-axes 1↔2): [G, D] = (d₁−d₂)·(E₁₂+E₂₁) conjugated by W. Consumed by
    compute_stable_mask + the clock trackers.

    kick_theta ≠ 0 starts the clock: M_prev is seeded at clock phase −kick_theta
    (the leapfrog then opens with Θ̇ = kick_theta/dt — caller picks dt scaling).

    Triple-buffer BC rule: writes M_prev (back-rotated if kicked), M_am, M_new.
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    d_iso = (1.0 + delta) / 3.0
    eps = ti.cast(1e-6, ti.f32)
    g_t = tensor_field.lc_g
    ck = ti.cos(kick_theta)
    sk = ti.sin(kick_theta)
    for i, j, k in ti.ndrange(nx, ny, nz):
        px = ti.cast(i - cx, ti.f32)
        py = ti.cast(j - cy, ti.f32)
        pz = ti.cast(k - cz, ti.f32)
        r2 = px * px + py * py + pz * pz
        r = ti.sqrt(r2)
        rho = ti.sqrt(px * px + py * py)

        rinv = 1.0 / ti.max(r, eps)
        rhat = ti.Vector([px * rinv, py * rinv, pz * rinv])
        ainv = 1.0 / ti.sqrt(rho * rho + eps)
        azim = ti.Vector([-py * ainv, px * ainv, 0.0])
        sdisc = ti.min(rho / rhoc_vox, 1.0)
        shrink = sdisc * sdisc * (3.0 - 2.0 * sdisc)
        ephi = azim * shrink
        ephi = ephi - ephi.dot(rhat) * rhat
        etheta = ephi.cross(rhat)

        srad = r / ti.sqrt(r2 + r0_vox * r0_vox)
        d0 = d_iso + srad * (1.0 - d_iso)
        d1 = d_iso + srad * (delta - d_iso)
        d2 = d_iso + srad * (0.0 - d_iso)

        # W eigenvector cols [r̂|e_Θ|e_Φ|t̂] (eig d0,d1,d2,g); boost mixes e_Θ(w1)↔t̂(w3).
        # Index-0 convention: time axis = array index 0 (Duda D=diag(g,1,δ,0)); spatial at 1,2,3.
        th = b_star * ti.exp(-(r2 / (rw_vox * rw_vox)))
        eth = ti.exp(th)                                  # no ti.cosh/sinh — exp identity
        chb = 0.5 * (eth + 1.0 / eth)
        shb = 0.5 * (eth - 1.0 / eth)
        w0 = ti.Vector([0.0, rhat[0], rhat[1], rhat[2]])              # eig d0 (r̂ at idx 1-3)
        w1 = ti.Vector([shb, etheta[0] * chb, etheta[1] * chb, etheta[2] * chb])
        w2 = ti.Vector([0.0, ephi[0], ephi[1], ephi[2]])             # eig d2 (e_Φ at idx 1-3)
        w3 = ti.Vector([chb, etheta[0] * shb, etheta[1] * shb, etheta[2] * shb])
        m = (d0 * w0.outer_product(w0) + d1 * w1.outer_product(w1)
             + d2 * w2.outer_product(w2) + g_t * w3.outer_product(w3))
        # clock tangent ∂M/∂Θ: [G_(1,2), D] = (d1−d2)(E₁₂+E₂₁) conjugated by W
        m_psi = (d1 - d2) * (w1.outer_product(w2) + w2.outer_product(w1))
        # kicked previous step: clock phase −kick_theta rotates eigen-cols 1↔2
        w1k = w1 * ck - w2 * sk
        w2k = w2 * ck + w1 * sk
        m_prev = (d0 * w0.outer_product(w0) + d1 * w1k.outer_product(w1k)
                  + d2 * w2k.outer_product(w2k) + g_t * w3.outer_product(w3))
        tensor_field.M_am[i, j, k] = m
        tensor_field.M_prev_am[i, j, k] = m_prev
        tensor_field.M_new_am[i, j, k] = m
        tensor_field.M_psi_am[i, j, k] = m_psi
        tensor_field.director_nhat[i, j, k] = rhat
        tensor_field.director_nhat_new[i, j, k] = rhat




# ================================================================
# M5.23.2 — RESEARCH-ENDPOINT LOADER (npz → launcher seed)
# ================================================================
# Loads a converged research state (fixed-J isorotation endpoints, the μ/τ
# census states, any saved M array) into a live TensorField as an INITIAL
# CONDITION (not a launcher equilibrium: drift under free evolution of a
# pinned-bc source is physics, not a loader bug).
#
# CONVENTIONS (measured at the loader build, m5_23_2 task_details):
#   - research npz storage keeps M[0,0] = +g (plain storage convention); the
#     production covariant vacuum is diag(−g, 1, δ, 0). V4 measured on the
#     fixed-J endpoint far field: 759.4/cell raw (+g) vs 2e-11 flipped (−g).
#     4×4 data is therefore flipped to −g AT LOAD (info["covariant"] = True
#     tells the canonical activation to skip its own flip); 3×3 data is
#     embedded block-diag at +g (the seeder convention) and takes the normal
#     flip_time_axis at canonical activation.
#   - the TensorField grid is ODD; research grids are even (32³/64³).
#     GRID FIT: launcher grid == n_src − 1 → CROP (drop the last plane per
#     axis, sponge-zone content); launcher grid ≥ n_src → EMBED centered
#     with edge-extend padding (the source boundary IS the far field).
#     Anything else → info["ok"] = False (caller falls back to vacuum).
#     The data array itself is never resampled.


def load_npz_M(tensor_field, npz_path, delta_cfg=None):
    """Load a research M-state npz into all three M buffers (Ṁ = 0).

    Returns an info dict: ok, n_src, blk, fit, covariant, h, delta, bc, warn
    (warn = list of warning strings; caller prints them).
    """
    import numpy as np

    info = {
        "ok": False,
        "n_src": 0,
        "blk": 0,
        "fit": "",
        "covariant": False,
        "h": 0.0,
        "delta": float("nan"),
        "bc": "?",
        "warn": [],
    }
    z = np.load(npz_path)
    if "M" not in z.files:
        info["warn"].append(f"no 'M' array in {npz_path}")
        return info
    m_src = z["M"].astype(np.float32)
    if m_src.ndim != 5 or m_src.shape[3] not in (3, 4) or m_src.shape[3] != m_src.shape[4]:
        info["warn"].append(f"unsupported M shape {m_src.shape}")
        return info
    n_src, blk = int(m_src.shape[0]), int(m_src.shape[3])
    info["n_src"], info["blk"] = n_src, blk
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    if (nx, ny, nz) == (n_src - 1, n_src - 1, n_src - 1):
        m_fit = m_src[:nx, :ny, :nz]
        info["fit"] = "crop"
    elif min(nx, ny, nz) >= n_src:
        ox, oy, oz = (nx - n_src) // 2, (ny - n_src) // 2, (nz - n_src) // 2
        m_fit = np.pad(
            m_src,
            (
                (ox, nx - n_src - ox),
                (oy, ny - n_src - oy),
                (oz, nz - n_src - oz),
                (0, 0),
                (0, 0),
            ),
            mode="edge",
        )
        info["fit"] = f"embed(+{ox},+{oy},+{oz})"
    else:
        info["warn"].append(
            f"grid mismatch: launcher {nx}x{ny}x{nz} vs source {n_src}^3 "
            "(need n_src-1 to crop or >= n_src to embed)"
        )
        return info
    m4 = np.zeros((nx, ny, nz, 4, 4), np.float32)
    if blk == 4:
        m4[:] = m_fit
        m4[..., 0, 0] *= -1.0  # storage +g → covariant −g (measured, header note)
        info["covariant"] = True
    else:  # 3×3 spatial block → the seeder's block-diag embed at +g
        m4[..., 0, 0] = tensor_field.lc_g
        m4[..., 1:, 1:] = m_fit
    tensor_field.M_am.from_numpy(m4)
    tensor_field.M_prev_am.from_numpy(m4)  # Ṁ = 0 (SET-J re-kicks if configured)
    tensor_field.M_new_am.from_numpy(m4)
    if "h" in z.files:
        info["h"] = float(z["h"])
    if "delta" in z.files:
        info["delta"] = float(z["delta"])
        if delta_cfg is not None and abs(info["delta"] - float(delta_cfg)) > 1e-6:
            info["warn"].append(
                f"source delta {info['delta']} != config delta {float(delta_cfg)} "
                "(V4 well mismatch)"
            )
    if "bc" in z.files:
        info["bc"] = str(z["bc"])
    info["ok"] = True
    return info
