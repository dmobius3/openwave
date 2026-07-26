"""
M5 ENGINE — PDE (engine2_pde)

The complete PDE layer: the reusable differential operators AND the equations
integrated with them — merged from the former operators + PDE modules
(2026-05-25 SoC consolidation). Two internal sections below:
  MATRIX OPERATORS       — @ti.func commutator, matrix ∇², antisym Â_μ, eigensolve
  PDE DYNAMICS           — evolve_M (matrix leapfrog), V_M / dV_M (LdG potential),
                           compute_curvature_flux, eigen_decompose,
                           relax_director_step (dissipative gradient descent)

Concern order: engine1_seeds → engine2_pde → engine3_observables → engine4_render.
M5.4 watch-item: if eigen_decompose ends up shared by render + observables + the
PDE, re-extract a small operators module at that point.
"""

import taichi as ti

# ================================================================
# MATRIX OPERATORS (M5.4 — Landau–de Gennes substrate M = O·D·O^T)
# ================================================================
# The Vector(3) toolkit above operates on the retiring ψ field. These operate on
# the real-symmetric 3×3 order parameter M stored in tensor_field.M_am. All three
# primitives were de-risked in research/sandbox_v3/m5_3_matrix_feasibility.py
# (storage round-trip, commutator vs analytic, ti.sym_eig director recovery 0.9995)
# before this production port.
#
#   commutator         — [A, B] = A·B − B·A           pure matrix @ti.func
#   principal_director — M → (principal eigenvector, eigenvalues) via ti.sym_eig
#   matrix_laplacian   — ∇²M component-wise            → 3×3; 1-cell halo
#   antisym_A_mu       — A_μ = [M, ∂_μ M] (paper Eq.19) → (A_x, A_y, A_z); 1-cell halo
#   eigen_decompose    — full-grid kernel: M_am → director_nhat + eigenvalues
#
# INDEX-GENERIC NOTE (M5.8): these are written for 3×3 so the eventual 4×4 SO(1,3)
# promotion (D = diag(g,1,δ,0)) is a type change, not a rewrite. matrix_laplacian
# and antisym_A_mu are dimension-agnostic; only principal_director / eigen_decompose
# assume a 3-vector director (the 4D clock readout gets its own kernel in M5.8).
#
# WHAT M5.4 USES vs DEFERS: the M5.4 Coulomb gate is director-equivalent (relax the
# eigenvector with the existing M5.1 relax, rebuild M), so eigen_decompose is on the
# critical path and commutator / matrix_laplacian / antisym_A_mu are built + unit-
# tested here but NOT yet driving dynamics. The curvature F_μν = ∂_μA_ν − ∂_νA_μ
# (Eq.20) and the full Eq.18 action that consume A_μ land in M5.5.


@ti.func
def commutator(a, b):  # type: ignore
    """Matrix commutator [A, B] = A·B − B·A — building block of Eq.19 A_μ + Eq.20 F_μν."""
    return a @ b - b @ a


@ti.func
def principal_director(m):  # type: ignore
    """M (3×3 symmetric) → (principal eigenvector n̂, eigenvalues vec3 [λ₁≥λ₂≥λ₃]).

    ANALYTIC symmetric-3×3 eigensolver: Cardano closed-form eigenvalues + the
    principal eigenvector as the largest row cross-product of (M − λ₁I). Replaces
    `ti.sym_eig` (M5.6.5a-fix), which is accurate for UNIAXIAL (degenerate) M but
    CATASTROPHICALLY WRONG for BIAXIAL M on Metal/f32 — eigenvalue error ~0.48 vs
    numpy for distinct eigenvalues (and f64 is unavailable: SPIRV codegen fails).
    This analytic solver matches numpy to f32 precision on biaxial M (eigenvalue
    err ~6e-6, director err ~2e-7 over 20k random symmetric matrices); the M5.6
    biaxial substrate (3 distinct eigenvalues) needs it. Returns eigenvalues already
    sorted descending (λ₁≥λ₂≥λ₃).
    """
    a00, a11, a22 = m[0, 0], m[1, 1], m[2, 2]
    a01, a02, a12 = m[0, 1], m[0, 2], m[1, 2]
    q = (a00 + a11 + a22) / 3.0
    p1 = a01 * a01 + a02 * a02 + a12 * a12
    p2 = (a00 - q) ** 2 + (a11 - q) ** 2 + (a22 - q) ** 2 + 2.0 * p1
    p = ti.sqrt(ti.max(p2 / 6.0, 1e-30))  # guard isotropic (p2→0)
    # B = (M − qI)/p ;  r = det(B)/2, clamped to [−1, 1] (numerical safety for acos)
    b00, b11, b22 = (a00 - q) / p, (a11 - q) / p, (a22 - q) / p
    b01, b02, b12 = a01 / p, a02 / p, a12 / p
    detB = (
        b00 * (b11 * b22 - b12 * b12)
        - b01 * (b01 * b22 - b12 * b02)
        + b02 * (b01 * b12 - b11 * b02)
    )
    r = ti.max(-1.0, ti.min(1.0, detB / 2.0))
    phi = ti.acos(r) / 3.0
    two_pi_3 = 2.0943951023931953  # 2π/3
    eig1 = q + 2.0 * p * ti.cos(phi)  # largest
    eig3 = q + 2.0 * p * ti.cos(phi + two_pi_3)  # smallest
    eig2 = 3.0 * q - eig1 - eig3  # middle
    # principal eigenvector = null space of (M − eig1·I): largest of the 3 row cross-products
    row0 = ti.Vector([a00 - eig1, a01, a02])
    row1 = ti.Vector([a01, a11 - eig1, a12])
    row2 = ti.Vector([a02, a12, a22 - eig1])
    c01, c02, c12 = row0.cross(row1), row0.cross(row2), row1.cross(row2)
    n01, n02, n12 = c01.norm_sqr(), c02.norm_sqr(), c12.norm_sqr()
    nvec = c01
    best = n01
    if n02 > best:
        nvec, best = c02, n02
    if n12 > best:
        nvec, best = c12, n12
    n = nvec / (ti.sqrt(best) + 1e-20)
    return n, ti.Vector([eig1, eig2, eig3])


@ti.func
def eigvec_for(m, lam):  # type: ignore
    """Unit eigenvector of symmetric 3×3 `m` for a KNOWN eigenvalue `lam`.

    Null space of (m − lam·I), taken as the largest of the 3 row cross-products
    (same construction as the principal eigenvector in `principal_director`, but
    for an arbitrary supplied eigenvalue — used for the MIDDLE eigenvector, the
    δ "clock-hand" axis, VIZ.3). Robust when `lam` is well-separated; degenerate
    (repeated-eigenvalue) cases give an arbitrary-but-unit vector in the
    eigenspace, which is the correct apolar behaviour for a uniaxial M.
    """
    a00, a11, a22 = m[0, 0], m[1, 1], m[2, 2]
    a01, a02, a12 = m[0, 1], m[0, 2], m[1, 2]
    row0 = ti.Vector([a00 - lam, a01, a02])
    row1 = ti.Vector([a01, a11 - lam, a12])
    row2 = ti.Vector([a02, a12, a22 - lam])
    c01, c02, c12 = row0.cross(row1), row0.cross(row2), row1.cross(row2)
    n01, n02, n12 = c01.norm_sqr(), c02.norm_sqr(), c12.norm_sqr()
    vec = c01
    best = n01
    if n02 > best:
        vec, best = c02, n02
    if n12 > best:
        vec, best = c12, n12
    return vec / (ti.sqrt(best) + 1e-20)


@ti.func
def matrix_laplacian(
    tensor_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """∇²M at voxel [i, j, k] — the 6-point stencil applied entry-wise to M.

    Taichi matrix arithmetic is element-wise for + and scalar ·, so the same
    stencil that gives the Vector(3) Laplacian gives the matrix Laplacian with no
    per-component loop. Reads tensor_field.M_am; returns a 3×3 matrix in [1/am].
    1-cell halo (caller skips boundary). Used by the matrix evolution (M5.5) and
    by operator-verification tests in M5.4.
    """
    face_sum = (
        tensor_field.M_am[i + 1, j, k]
        + tensor_field.M_am[i - 1, j, k]
        + tensor_field.M_am[i, j + 1, k]
        + tensor_field.M_am[i, j - 1, k]
        + tensor_field.M_am[i, j, k + 1]
        + tensor_field.M_am[i, j, k - 1]
    )
    center = tensor_field.M_am[i, j, k]
    return (face_sum - 6.0 * center) / (tensor_field.dx_am**2)


@ti.func
def antisym_A_mu(
    tensor_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """Paper Eq.19 antisymmetric connection A_μ = [M, ∂_μ M] for μ ∈ {x, y, z}.

    Central-differences M along each axis, then commutes with M at the voxel.
    Returns the three antisymmetric 3×3 matrices (A_x, A_y, A_z). The M5.3 spike
    verified the spatial-derivative commutator [∂_xM, ∂_yM] over a field; this is
    the same machinery in the [M, ∂_μM] arrangement. 1-cell halo. Consumed by the
    Eq.20 curvature F_μν and the Eq.18 action — both assembled in M5.5.
    """
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    m = tensor_field.M_am[i, j, k]
    dM_x = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
    dM_y = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
    dM_z = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
    return commutator(m, dM_x), commutator(m, dM_y), commutator(m, dM_z)


@ti.func
def compute_laplacian_director(
    tensor_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """∇²n̂ at voxel [i, j, k] — 6-point stencil on the DERIVED director field.

    Identical form to compute_laplacian, but reads tensor_field.director_nhat (the
    matrix-substrate director, eigenvector of M) rather than the retiring ψ. The
    M5.4 director-equivalent relaxation gradient-descends the Frank energy on this
    field. 1-cell halo (caller skips boundary). (compute_laplacian on ψ stays for
    the retiring wave leapfrog until the M5.4 cleanup deletes it.)
    """
    face_sum = (
        tensor_field.director_nhat[i + 1, j, k]
        + tensor_field.director_nhat[i - 1, j, k]
        + tensor_field.director_nhat[i, j + 1, k]
        + tensor_field.director_nhat[i, j - 1, k]
        + tensor_field.director_nhat[i, j, k + 1]
        + tensor_field.director_nhat[i, j, k - 1]
    )
    center = tensor_field.director_nhat[i, j, k]
    return (face_sum - 6.0 * center) / (tensor_field.dx_am**2)


@ti.kernel
def rebuild_M_from_director(tensor_field: ti.template(), delta: ti.f32):  # type: ignore
    """Rebuild the uniaxial order parameter M = δ·I + (1−δ)·n̂⊗n̂ from director_nhat.

    Called after the director-equivalent relaxation so M_am, the `‖M−D‖_F`/`‖Ṁ‖_F`
    trackers, and the flux-mesh orientation-deviation mode stay consistent with the
    relaxed director. Writes all three M buffers (M_prev = M = M_new) — BC-consistent
    for the M5.5 matrix leapfrog, and `Ṁ = 0` for a static relaxation (clock ω = 0).
    The uniaxial inverse of eigen_decompose.
    """
    eye = ti.Matrix.identity(ti.f32, 3)
    g = tensor_field.lc_g  # time-axis (index 0) eigenvalue (Duda D=diag(g,1,δ,0))
    for i, j, k in tensor_field.director_nhat:
        n = tensor_field.director_nhat[i, j, k]
        msp = delta * eye + (1.0 - delta) * n.outer_product(n)  # 3×3 spatial
        m = ti.Matrix.zero(ti.f32, 4, 4)  # embed block-diag(g, spatial)
        for a in ti.static(range(3)):
            for bb in ti.static(range(3)):
                m[a + 1, bb + 1] = msp[a, bb]
        m[0, 0] = g
        tensor_field.M_am[i, j, k] = m
        tensor_field.M_prev_am[i, j, k] = m
        tensor_field.M_new_am[i, j, k] = m


@ti.kernel
def eigen_decompose(tensor_field: ti.template()):  # type: ignore
    """THE LYNCHPIN (m5_4b_rendering_features.md): refresh director_nhat + eigenvalues from M_am.

    Per voxel: eigen-decompose M_am, store the sorted spectrum (λ₁≥λ₂≥λ₃) in
    `eigenvalues` and the principal eigenvector in `director_nhat`. The whole
    rendering stack (glyphs, vector_warp, granule director) and the redefined
    amplitude tracker (‖M−D‖_F) read these derived fields, so this runs once/frame
    after M changes — never in an inner loop (M5.3 cost: ~1.1× a Vector(3)
    Laplacian at 256³).

    SIGN CONTINUITY: a nematic director is apolar (n ≡ −n), and `ti.sym_eig`
    returns an arbitrary eigenvector sign per voxel. A spurious sign flip between
    neighbours fakes a huge ∇n̂ (wrong Frank energy). We resolve it by aligning the
    fresh eigenvector to the EXISTING director_nhat (flip if n̂·n̂_prev < 0). The
    matrix seeders write director_nhat with the correct physical sign first, so the
    branch stays smooth from seed onward; on a cold field (n̂_prev = 0) the dot is 0
    and no flip occurs.
    """
    for i, j, k in tensor_field.M_am:
        m = tensor_field.M_am[i, j, k]
        # spatial 3×3 block (indices 1,2,3; the time/g axis at index 0 is excluded) — index-0 convention
        msp = ti.Matrix(
            [[m[1, 1], m[1, 2], m[1, 3]], [m[2, 1], m[2, 2], m[2, 3]], [m[3, 1], m[3, 2], m[3, 3]]]
        )
        n, evals = principal_director(msp)
        # sort eigenvalues descending λ₁≥λ₂≥λ₃ (3-element selection sort)
        e0, e1, e2 = evals[0], evals[1], evals[2]
        lo = ti.min(e0, ti.min(e1, e2))
        hi = ti.max(e0, ti.max(e1, e2))
        mid = e0 + e1 + e2 - lo - hi
        tensor_field.eigenvalues[i, j, k] = ti.Vector([hi, mid, lo])
        # sign-continuous principal director
        if n.dot(tensor_field.director_nhat[i, j, k]) < 0.0:
            n = -n
        tensor_field.director_nhat[i, j, k] = n
        # VIZ.3: middle (δ) eigenvector = the "clock-hand" axis that sweeps around
        # the director under the Zitterbewegung twist. Same apolar sign-continuity.
        nm = eigvec_for(msp, mid)
        if nm.dot(tensor_field.director_mid[i, j, k]) < 0.0:
            nm = -nm
        tensor_field.director_mid[i, j, k] = nm


# ================================================================
# MATRIX EVOLUTION (M5.5) — Eq.18 action leapfrog
# ================================================================
# The first LIVE matrix-field dynamics: evolve M under Duda's Eq.18 action with a
# simple ½‖Ṁ‖² kinetic + the FAITHFUL Eq.18 potential (decision 2026-05-26):
#     U(M) = 4 Σ_{μ<ν} ‖[M_μ, M_ν]‖²_F  +  V_M(M)          (M_μ = ∂_μ M)
#     EOM:  ∂²_t M = c²·force_curv − dV_M/dM ,   force_curv = Σ_α ∂_α G_α
#     G_α = ∂U_curv/∂M_α = 8 Σ_ν [[M_α,M_ν], M_ν]          (symmetric)
# Validated in sandbox_v5 (m5_5_2/3); see 5a §9. The faithful curvature kinetic
# (F_μ0² = 4‖[M_μ,Ṁ]‖², degenerate metric) is the M5.6 refinement.
#
# Two-pass per step: compute_curvature_flux (G_α everywhere, 1-cell halo) → evolve_M
# (divergence of G + V-force, leapfrog, 2-cell halo) → tensor_field.swap_matrix_buffers().


@ti.func
def V_M(m, a: ti.f32, b: ti.f32, c: ti.f32):  # type: ignore
    """Eq.13 Landau–de Gennes Higgs-like potential V = a·Tr(M²) − b·Tr(M³) + c·(Tr(M²))².

    Exact (a,b,c) giving the Λ=(1,δ,0) vacuum is Q7 (Duda open); a=b=c=0 → free-curvature
    dynamics (V off). Rotation-invariant ⇒ acts only on the eigenvalue/regularization
    sector (m5_5_3 finding).
    """
    # M5.8.1 — act on the SPATIAL 3×3 block ONLY. The time axis (index 0, eigenvalue g)
    # is boost-decoupled and must NOT feed Tr(M²)/Tr(M³): including g²(=64)/g³(=512)
    # would inflate the LdG potential and its force by orders of magnitude → blow-up.
    msp = ti.Matrix(
        [[m[1, 1], m[1, 2], m[1, 3]], [m[2, 1], m[2, 2], m[2, 3]], [m[3, 1], m[3, 2], m[3, 3]]]
    )
    m2 = msp @ msp
    tr2 = m2.trace()
    tr3 = (m2 @ msp).trace()
    return a * tr2 - b * tr3 + c * tr2 * tr2


@ti.func
def dV_M(m, a: ti.f32, b: ti.f32, c: ti.f32):  # type: ignore
    """∂V_LG/∂M = 2a·M − 3b·M² + 4c·Tr(M²)·M, on the SPATIAL 3×3 block; time row/col
    force = 0 (g decoupled — M5.8.1). Mirrors V_M's spatial-only restriction."""
    msp = ti.Matrix(
        [[m[1, 1], m[1, 2], m[1, 3]], [m[2, 1], m[2, 2], m[2, 3]], [m[3, 1], m[3, 2], m[3, 3]]]
    )
    m2 = msp @ msp
    tr2 = m2.trace()
    dsp = 2.0 * a * msp - 3.0 * b * m2 + 4.0 * c * tr2 * msp  # 3×3
    d4 = ti.Matrix.zero(ti.f32, 4, 4)
    for i in ti.static(range(3)):
        for j in ti.static(range(3)):
            d4[i + 1, j + 1] = dsp[i, j]
    return d4


@ti.kernel
def compute_curvature_flux(tensor_field: ti.template()):  # type: ignore
    """Curvature flux G_α = 8 Σ_ν [[M_α,M_ν],M_ν] = ∂U_curv/∂M_α → curv_flux_{x,y,z}.

    M_α = ∂_α M (central diff of M_am). [[A,B],B] is symmetric (antisym⊗sym), so each
    G_α is symmetric. 1-cell halo (caller skips boundary). Run before evolve_M.
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        # G_α = 8 Σ_{ν≠α} [[M_α,M_ν],M_ν]  (ν=α term is zero)
        tensor_field.curv_flux_x[i, j, k] = 8.0 * (
            commutator(commutator(mx, my), my) + commutator(commutator(mx, mz), mz)
        )
        tensor_field.curv_flux_y[i, j, k] = 8.0 * (
            commutator(commutator(my, mx), mx) + commutator(commutator(my, mz), mz)
        )
        tensor_field.curv_flux_z[i, j, k] = 8.0 * (
            commutator(commutator(mz, mx), mx) + commutator(commutator(mz, my), my)
        )


@ti.kernel
def evolve_M(
    tensor_field: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    a: ti.f32,  # type: ignore
    b: ti.f32,  # type: ignore
    c: ti.f32,  # type: ignore
):
    """Leapfrog M one step under the Eq.18 action: ∂²_t M = c²·(Σ_α ∂_α G_α) − dV_M/dM.

    M_new = 2M − M_prev + (dt)²·[c²·div(G) − dV_M(M)].  Reads M_am, M_prev_am,
    curv_flux_* (from compute_curvature_flux); writes M_new_am. Caller then calls
    tensor_field.swap_matrix_buffers(). 2-cell halo (G is a 1-cell-halo quantity, its
    divergence needs one more). Dirichlet BC: boundary M_new left untouched (seeders
    write all three M buffers, per the triple-buffer BC rule).
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    c2 = c_amrs * c_amrs
    dt2 = dt_rs * dt_rs
    for i, j, k in ti.ndrange((2, nx - 2), (2, ny - 2), (2, nz - 2)):
        div_G = (
            (tensor_field.curv_flux_x[i + 1, j, k] - tensor_field.curv_flux_x[i - 1, j, k])
            + (tensor_field.curv_flux_y[i, j + 1, k] - tensor_field.curv_flux_y[i, j - 1, k])
            + (tensor_field.curv_flux_z[i, j, k + 1] - tensor_field.curv_flux_z[i, j, k - 1])
        ) * inv_2dx
        force = c2 * div_G - dV_M(tensor_field.M_am[i, j, k], a, b, c)
        m_new = 2.0 * tensor_field.M_am[i, j, k] - tensor_field.M_prev_am[i, j, k] + dt2 * force
        # M5.8.1 — freeze the time axis (index 0): a constant, boost-decoupled g
        # background. The curvature force already preserves the block; this also pins
        # it under V-on (where dV_M would otherwise nudge M[0,0]). M5.8.2 replaces this
        # with the real Minkowski time dynamics.
        g_here = tensor_field.M_am[i, j, k][0, 0]
        m_new[1, 0] = 0.0
        m_new[0, 1] = 0.0
        m_new[2, 0] = 0.0
        m_new[0, 2] = 0.0
        m_new[3, 0] = 0.0
        m_new[0, 3] = 0.0
        m_new[0, 0] = g_here
        tensor_field.M_new_am[i, j, k] = m_new


# ================================================================
# M5.8.2c — 4D MINKOWSKI EVOLUTION (flag-gated; the time axis goes LIVE)
# ================================================================
# The signed curvature flux is the 3D production flux with F → ηFη
# (η = diag(−1,1,1,1), time = matrix index 0) — the 5a §10d ℋ sign rule:
# spatial matrix pairs positive, (α,0) pairs negative. Sandbox anchors:
# m5_8_2a (fuel) → 2b-1 (CC) → 2b-2 (field clock) → 2c-1 (full nonlinear,
# f64 numpy — THE cross-validation reference for these f32 kernels).
# Ghost guard: the per-voxel stable_mask (K>0 ∧ Q(∇ψ) positive-definite,
# computed ONCE post-seed) selects signed flux on the stable region and the
# always-stable Euclidean flux on the fuel shell (the 2b-2 λ≈15.6/t linear
# runaway sector — its propulsion physics needs the constrained faithful
# kernel, future work; v1 keeps it Euclidean-stable). Global guard: the
# coherent (α,0) velocity drift is sampled on 3 mid-planes (Metal-safe —
# NO full-grid reductions, the atomics lesson) and subtracted in evolve.


@ti.func
def eta_twist_masked(f, stable: ti.f32):  # type: ignore
    """F → ηFη blended by the stable mask: (α,0)/(0,α) comps × (1 − 2·stable).

    stable=1 → the Minkowski-signed twist; stable=0 → identity (Euclidean
    fallback on the fuel shell)."""
    s = 1.0 - 2.0 * stable
    out = f
    for a_ in ti.static(range(1, 4)):
        out[a_, 0] = f[a_, 0] * s
        out[0, a_] = f[0, a_] * s
    return out


@ti.func
def signed_dot4(a4, b4):  # type: ignore
    """⟨A,B⟩_s = Σ A∘(ηBη): spatial-pair comps +, (α,0) comps − (full-matrix sum)."""
    acc = 0.0
    for p_ in ti.static(range(4)):
        for q_ in ti.static(range(4)):
            sgn = 1.0
            if ti.static((p_ == 0) != (q_ == 0)):
                sgn = -1.0
            acc += a4[p_, q_] * b4[p_, q_] * sgn
    return acc


@ti.kernel
def compute_stable_mask(tensor_field: ti.template()):  # type: ignore
    """Per-voxel ghost guard (run ONCE post-seed): stable_mask = 1 where the
    signed clock kernel is well-posed — K(x) > 0 AND the 3×3 gradient-stiffness
    form Q(∇ψ) is positive-definite (the 2b-2 criterion, exact):

        P_i = [M_i, M_ψ] ,  p_ij = ⟨P_i, P_j⟩_s ,  K = 4 Σ_i p_ii
        Q_ii = Σ_(j≠i) p_jj ,  Q_ij = −p_ij  (i≠j) ;  stable ⇔ K>0 ∧ minEig(Q)>0

    minEig via the analytic Cardano solver (NOT ti.sym_eig — the Metal/f32
    lesson). Boundary voxels default to 0 (Euclidean). Seed-time cost only.
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange(nx, ny, nz):
        tensor_field.stable_mask[i, j, k] = 0.0
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        m_psi = tensor_field.M_psi_am[i, j, k]
        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        p0 = commutator(mx, m_psi)
        p1 = commutator(my, m_psi)
        p2 = commutator(mz, m_psi)
        p00 = signed_dot4(p0, p0)
        p11 = signed_dot4(p1, p1)
        p22 = signed_dot4(p2, p2)
        p01 = signed_dot4(p0, p1)
        p02 = signed_dot4(p0, p2)
        p12 = signed_dot4(p1, p2)
        kin = 4.0 * (p00 + p11 + p22)
        q = ti.Matrix([[p11 + p22, -p01, -p02], [-p01, p00 + p22, -p12], [-p02, -p12, p00 + p11]])
        _, eigs = principal_director(q)  # sorted λ₁≥λ₂≥λ₃
        qscale = ti.abs(q[0, 0]) + ti.abs(q[1, 1]) + ti.abs(q[2, 2]) + 1e-30
        if kin > 0.0 and eigs[2] > 1e-6 * qscale:
            tensor_field.stable_mask[i, j, k] = 1.0


@ti.kernel
def compute_tstar(tensor_field: ti.template()):  # type: ignore
    """Seed the per-voxel V(M) amplitude target t*(x) = Tr(M_sp²) of the CURRENT
    field (call once, right after seeding): the dressed seed becomes the exact
    V-equilibrium — zero static V-force (the 2c-0 'pin the well to the DRESSED
    amplitude' design input; fixes the 2c-2 slow V-pump on the off-minimum
    dressed shell). The well then confines amplitude EXCURSIONS, not the seed."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for i, j, k in ti.ndrange(nx, ny, nz):
        m = tensor_field.M_am[i, j, k]
        msp = ti.Matrix(
            [[m[1, 1], m[1, 2], m[1, 3]], [m[2, 1], m[2, 2], m[2, 3]], [m[3, 1], m[3, 2], m[3, 3]]]
        )
        tensor_field.ldg_tstar[i, j, k] = (msp @ msp).trace()


@ti.func
def dV_M_dressed(m, cc: ti.f32, tstar: ti.f32):  # type: ignore
    """∂V/∂M for the DRESSED well V = cc·(Tr(M_sp²) − t*(x))²: 4·cc·(t−t*)·M_sp,
    spatial block only (time row/col force = 0 — the M5.8.1 rule)."""
    msp = ti.Matrix(
        [[m[1, 1], m[1, 2], m[1, 3]], [m[2, 1], m[2, 2], m[2, 3]], [m[3, 1], m[3, 2], m[3, 3]]]
    )
    t = (msp @ msp).trace()
    dsp = 4.0 * cc * (t - tstar) * msp
    d4 = ti.Matrix.zero(ti.f32, 4, 4)
    for i in ti.static(range(3)):
        for j in ti.static(range(3)):
            d4[i + 1, j + 1] = dsp[i, j]
    return d4


@ti.kernel
def compute_curvature_flux_4d(tensor_field: ti.template()):  # type: ignore
    """The 4D signed curvature flux G_α = 8 Σ_ν [tw(F_αν), M_ν], tw = ηFη on the
    stable region / identity on the fuel shell (stable_mask blend). At b=0 the
    seed has no (α,0) components ⇒ tw is a no-op ⇒ EXACTLY compute_curvature_flux
    (the headless identity check). Same two-pass contract as the 3D kernel."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        st = tensor_field.stable_mask[i, j, k]
        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        fxy = eta_twist_masked(commutator(mx, my), st)
        fxz = eta_twist_masked(commutator(mx, mz), st)
        fyz = eta_twist_masked(commutator(my, mz), st)
        tensor_field.curv_flux_x[i, j, k] = 8.0 * (commutator(fxy, my) + commutator(fxz, mz))
        tensor_field.curv_flux_y[i, j, k] = 8.0 * (
            commutator(-1.0 * fxy, mx) + commutator(fyz, mz)
        )
        tensor_field.curv_flux_z[i, j, k] = 8.0 * (
            commutator(-1.0 * fxz, mx) + commutator(-1.0 * fyz, my)
        )


@ti.kernel
def sample_v03_drift(tensor_field: ti.template(), dt_rs: ti.f32):  # type: ignore
    """Sample the coherent (α,0) velocity drift on 3 mid-planes (Metal-safe:
    plane-scale atomics only — the full-grid-reduction lesson). Results land in
    tensor_field.v03_sums[0..2]; the caller divides by the plane-voxel count."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for a_ in ti.static(range(3)):
        tensor_field.v03_sums[a_] = 0.0
    inv_dt = 1.0 / dt_rs
    im, jm, km = nx // 2, ny // 2, nz // 2
    for j, k in ti.ndrange(ny, nz):
        d = (tensor_field.M_am[im, j, k] - tensor_field.M_prev_am[im, j, k]) * inv_dt
        for a_ in ti.static(range(3)):
            tensor_field.v03_sums[a_] += d[a_ + 1, 0]
    for i, k in ti.ndrange(nx, nz):
        d = (tensor_field.M_am[i, jm, k] - tensor_field.M_prev_am[i, jm, k]) * inv_dt
        for a_ in ti.static(range(3)):
            tensor_field.v03_sums[a_] += d[a_ + 1, 0]
    for i, j in ti.ndrange(nx, ny):
        d = (tensor_field.M_am[i, j, km] - tensor_field.M_prev_am[i, j, km]) * inv_dt
        for a_ in ti.static(range(3)):
            tensor_field.v03_sums[a_] += d[a_ + 1, 0]


@ti.kernel
def evolve_M_4d(
    tensor_field: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    ldg_cc: ti.f32,  # type: ignore
    vm0: ti.f32,  # type: ignore
    vm1: ti.f32,  # type: ignore
    vm2: ti.f32,  # type: ignore
    km: ti.f32,  # type: ignore
):
    """The M5.8.2 leapfrog: the M5.8.1 time-freeze clamp is REPLACED by the soft
    global guard — the time axis is LIVE. Identical to evolve_M except:

    (i) no (α,0) zeroing and no M[0,0] pin (boost dressing + clock evolve freely);
    (ii) the sampled coherent (α,0) drift (vm0..2, from sample_v03_drift) is
        subtracted — the 2b-1 ghost channel (a free GLOBAL dressing mode);
    (iii) DIAGONAL FAITHFUL-LITE INERTIA m(x) = 1 + km·dx²·Σ_i‖M_i‖²_F dividing
        the acceleration — the scalar shadow of the faithful kinetic operator
        A(Ṁ) = 4Σ[η[Ṁ,M_i]η,M_i] that stabilized the 2c-1 spike: inertia grows
        with local structure (steepening modes get HEAVY and decelerate — the
        self-regulation the simple kinetic lacks, which let the signed
        stiffness pump the core, the 2c-2 repro). m → 1 in weak-gradient
        regions ⇒ the 3D/vacuum limit is untouched; km=0 disables.

    V(M) uses the DRESSED well dV_M_dressed with the per-voxel target t*(x)
    (compute_tstar, run once post-seed): the seed is the exact V-equilibrium —
    zero static V-force (the 2c-2 fix for the slow V-pump on the off-minimum
    dressed shell; 2c-0 energetics unchanged: confinement of excursions kept)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    c2 = c_amrs * c_amrs
    dt2 = dt_rs * dt_rs
    dx2 = tensor_field.dx_am * tensor_field.dx_am
    for i, j, k in ti.ndrange((2, nx - 2), (2, ny - 2), (2, nz - 2)):
        div_G = (
            (tensor_field.curv_flux_x[i + 1, j, k] - tensor_field.curv_flux_x[i - 1, j, k])
            + (tensor_field.curv_flux_y[i, j + 1, k] - tensor_field.curv_flux_y[i, j - 1, k])
            + (tensor_field.curv_flux_z[i, j, k + 1] - tensor_field.curv_flux_z[i, j, k - 1])
        ) * inv_2dx
        force = c2 * div_G - dV_M_dressed(
            tensor_field.M_am[i, j, k], ldg_cc, tensor_field.ldg_tstar[i, j, k]
        )
        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        mloc = 1.0 + km * dx2 * (mx.norm_sqr() + my.norm_sqr() + mz.norm_sqr())
        m_new = (
            2.0 * tensor_field.M_am[i, j, k]
            - tensor_field.M_prev_am[i, j, k]
            + dt2 * force * (1.0 / mloc)
        )
        # guard (a): remove the coherent global (α,0) drift accumulated this step
        m_new[1, 0] -= vm0 * dt_rs
        m_new[0, 1] -= vm0 * dt_rs
        m_new[2, 0] -= vm1 * dt_rs
        m_new[0, 2] -= vm1 * dt_rs
        m_new[3, 0] -= vm2 * dt_rs
        m_new[0, 3] -= vm2 * dt_rs
        tensor_field.M_new_am[i, j, k] = m_new


# ================================================================
# M5.8.2c OPTION B — THE CONSTRAINED SPECTRAL-PROJECTION INTEGRATOR
# ================================================================
# RESTORED 2026-06-05 late evening — this section was erased from the working
# tree by a stale-IDE-buffer save (the file was open with a pre-Option-B
# buffer; the GUI "even worse" run executed WITHOUT these kernels). Content
# identical to the B-2 wiring + the velocity-kick fix.
#
# The Minkowski-SIGNED dynamics, stable (INTEGRATOR_4D: "constrained"). The
# structural finding (2c-2): the signed flux under ANY cheap positive inertia
# has slow growing modes — only the faithful kinetic works. Per voxel, the
# inertia operator A(Ṁ) = 4 Σ_i [η[Ṁ,M_i]η, M_i] is eigendecomposed in the
# 10-dim orthonormal symmetric-matrix basis (OWN cyclic Jacobi on local
# matrices — NOT ti.sym_eig, the Metal/f32 lesson); only positive-inertia
# directions evolve (λ > EPS·max|λ|), and P is PROJECTED onto the kept
# subspace every step (frozen directions must not accumulate momentum — the
# 2c-1 v2 constraint-switching pump). State is (M, P): symplectic-Euler
# P += dt·force → clamp → solve Ṁ = A⁺P → M += dt·Ṁ. Validated machine-exact
# against the f64 numpy reference (sandbox_v8/m5_8_2cb_taichi_constrained.py:
# field diff 5.5e-15 @890 steps; f32 1e-5 benign; Metal clean; 8.3 ms/step
# solve at 64³ — B-1 gates [J1][B1][B2][B2m][B3][B4] all PASS, 2026-06-05).
# ⚠️ KNOWN OPEN (2cb-2 long-horizon): the f32 kernel develops a negative-H
# ghost runaway past ~1100 steps even in the validated 24³ config (H plateau
# → collapse; align 1.0 → 0.55 by 6000) — the f64-vs-f32 6000-step
# scheme-vs-precision discrimination is the active investigation.
#
# STATUS (M5.24 audit, 2026-07-19): the canonical registry names this stack
# the HISTORICAL PRECEDENT of the structural disease later proven at M5.20.3
# (the η-null constraint + indefinite H; research m5_theory_canonical.md § 2:
# "every cheap positive-inertia kinetic has slow growing modes — only the
# faithful kinetic works"). Its refined research descendant (the true-L EL
# solve) is DIAGNOSTICS-ONLY by verdict (§ 2 row 1: the free-EL IVP is
# ill-posed). This mode is KEPT as the era 4D path for the dressed-seed
# xperiments; the verified-L era production dynamics is the M5.24 CANONICAL
# section below.
#
# UNITS + NORMALIZATION: these kernels run the 2c-1 natural convention
# (c = 1, τ = c·t) — the caller passes dt_eff = c_amrs·dt_rs (am). Fluxes and
# A both use the 2c-1 4× normalization (the production 8× leapfrog flux is 2×
# this; a COMMON factor cancels exactly in Ṁ = A⁺P, so the dynamics here is
# the validated one verbatim). V enters the τ-units force as dV/(2c²) — the
# launcher passes ldg_cc_4d = LDG_STIFFNESS_K·0.5/dx⁴ (the c² cancels). The
# constrained path is FULLY signed (no stable-mask blend — the spectral
# projection IS the exact pointwise ghost treatment) and ignores
# KM_INERTIA_4D (the faithful A replaces the diagonal shadow).

EPS_EIG_4D = 0.05  # positive-inertia keep threshold (× max|λ| per voxel)
VCAP_4D = 5.0  # ‖Ṁ‖_F backstop (never engaged in the validated 900-step runs)
JACOBI_SWEEPS = 20  # max cyclic sweeps (converges in ~6–8)
JTOL2_4D = (2e-7) ** 2  # relative off-diagonal² stop (f32 eps scale)


@ti.func
def eta_twist(f):  # type: ignore
    """F → ηFη, unmasked — the constrained path is fully Minkowski-signed."""
    out = f
    for a_ in ti.static(range(1, 4)):
        out[a_, 0] = -f[a_, 0]
        out[0, a_] = -f[0, a_]
    return out


@ti.kernel
def init_P_4d(tensor_field: ti.template(), kick: ti.f32):  # type: ignore
    """P₀ = A(Ṁ₀) with Ṁ₀ = kick·M_ψ — the EXACT 2c-1 kick semantics (a
    VELOCITY in τ-units, dt-INDEPENDENT; kick = 0 ⇒ P₀ = 0 exactly). NOT the
    buffer encoding (M − M_prev)/dt: that made the kick velocity scale as
    1/dt — 16× hotter than validated at the constrained dt (the step-31
    v-cap engagement in the 2cb-2 repro). Run ONCE post-seed. Md_am gets Ṁ₀
    on the act region (read by the first flux call — the 2c-1 step-0
    semantics)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange(nx, ny, nz):
        tensor_field.P_am[i, j, k] = ti.Matrix.zero(ti.f32, 4, 4)
        tensor_field.Md_am[i, j, k] = ti.Matrix.zero(ti.f32, 4, 4)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        md = kick * tensor_field.M_psi_am[i, j, k]
        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        a4 = 4.0 * (
            commutator(eta_twist(commutator(md, mx)), mx)
            + commutator(eta_twist(commutator(md, my)), my)
            + commutator(eta_twist(commutator(md, mz)), mz)
        )
        tensor_field.P_am[i, j, k] = a4
        tensor_field.Md_am[i, j, k] = md * tensor_field.act4d[i, j, k]


@ti.kernel
def flux_4d_constrained(tensor_field: ti.template()):  # type: ignore
    """The 2c-1 flux (C1-pinned signs): G_i = ∂U/∂M_i − ∂T/∂M_i with
    dU_i = 4 Σ_{j≠i} [η F̃_ij η, M_j] and dT_i = −4 [η F_0i η, Ṁ] — fully
    signed, NO mask. Writes curv_flux_x/y/z (border stays zero)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        md = tensor_field.Md_am[i, j, k]
        mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        fxy = commutator(mx, my)
        fxz = commutator(mx, mz)
        fyz = commutator(my, mz)
        tensor_field.curv_flux_x[i, j, k] = 4.0 * (
            commutator(eta_twist(fxy), my)
            + commutator(eta_twist(fxz), mz)
            + commutator(eta_twist(commutator(md, mx)), md)
        )
        tensor_field.curv_flux_y[i, j, k] = 4.0 * (
            commutator(eta_twist(-1.0 * fxy), mx)
            + commutator(eta_twist(fyz), mz)
            + commutator(eta_twist(commutator(md, my)), md)
        )
        tensor_field.curv_flux_z[i, j, k] = 4.0 * (
            commutator(eta_twist(-1.0 * fxz), mx)
            + commutator(eta_twist(-1.0 * fyz), my)
            + commutator(eta_twist(commutator(md, mz)), md)
        )


@ti.kernel
def update_P_4d(tensor_field: ti.template(), dt_eff: ti.f32, ldg_cc4d: ti.f32):  # type: ignore
    """P += dt_eff·(div G − dV_dressed) — the τ-units force step (ldg_cc4d
    pre-scaled by the launcher to the 4× convention: K·0.5/dx⁴)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        div_g = (
            (tensor_field.curv_flux_x[i + 1, j, k] - tensor_field.curv_flux_x[i - 1, j, k])
            + (tensor_field.curv_flux_y[i, j + 1, k] - tensor_field.curv_flux_y[i, j - 1, k])
            + (tensor_field.curv_flux_z[i, j, k + 1] - tensor_field.curv_flux_z[i, j, k - 1])
        ) * inv_2dx
        force = div_g - dV_M_dressed(
            tensor_field.M_am[i, j, k], ldg_cc4d, tensor_field.ldg_tstar[i, j, k]
        )
        tensor_field.P_am[i, j, k] += dt_eff * force


@ti.kernel
def sample_p03_drift(tensor_field: ti.template()):  # type: ignore
    """3-mid-plane act-region sums of the (α,0) MOMENTUM components → v03_sums
    (the 2c-1 global clamp, plane-sampled — Metal-safe, NO full-grid atomics).
    The caller divides by the act-plane voxel count (launcher-precomputed)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for a_ in ti.static(range(3)):
        tensor_field.v03_sums[a_] = 0.0
    im, jm, km = nx // 2, ny // 2, nz // 2
    for j, k in ti.ndrange(ny, nz):
        if tensor_field.act4d[im, j, k] > 0.5:
            p = tensor_field.P_am[im, j, k]
            for a_ in ti.static(range(3)):
                tensor_field.v03_sums[a_] += p[a_ + 1, 0]
    for i, k in ti.ndrange(nx, nz):
        if tensor_field.act4d[i, jm, k] > 0.5:
            p = tensor_field.P_am[i, jm, k]
            for a_ in ti.static(range(3)):
                tensor_field.v03_sums[a_] += p[a_ + 1, 0]
    for i, j in ti.ndrange(nx, ny):
        if tensor_field.act4d[i, j, km] > 0.5:
            p = tensor_field.P_am[i, j, km]
            for a_ in ti.static(range(3)):
                tensor_field.v03_sums[a_] += p[a_ + 1, 0]


@ti.kernel
def apply_p03_clamp(tensor_field: ti.template(), m0: ti.f32, m1: ti.f32, m2: ti.f32):  # type: ignore
    """Guard (a): subtract the act-mean (α,0) momentum (act region only), then
    restore P symmetry everywhere — the exact 2c-1 clamp order."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for i, j, k in ti.ndrange(nx, ny, nz):
        if tensor_field.act4d[i, j, k] > 0.5:
            tensor_field.P_am[i, j, k][1, 0] -= m0
            tensor_field.P_am[i, j, k][2, 0] -= m1
            tensor_field.P_am[i, j, k][3, 0] -= m2
        tensor_field.P_am[i, j, k][0, 1] = tensor_field.P_am[i, j, k][1, 0]
        tensor_field.P_am[i, j, k][0, 2] = tensor_field.P_am[i, j, k][2, 0]
        tensor_field.P_am[i, j, k][0, 3] = tensor_field.P_am[i, j, k][3, 0]


@ti.kernel
def solve_constrained_4d(tensor_field: ti.template()):  # type: ignore
    """THE constrained solve (the B-1-validated kernel), fused per voxel:
    build the 10×10 inertia operator A in the sym basis → cyclic Jacobi on
    LOCAL matrices → positive-inertia keep → project P AND solve Ṁ = A⁺P on
    the kept subspace → act mask + ‖Ṁ‖ backstop. Zero-gradient voxels (grid
    border) give A = 0 ⇒ all directions frozen ⇒ P, Ṁ → 0 (the 2c-1 limit)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_2dx = 1.0 / (2.0 * tensor_field.dx_am)
    for i, j, k in ti.ndrange(nx, ny, nz):
        mx = ti.Matrix.zero(ti.f32, 4, 4)
        my = ti.Matrix.zero(ti.f32, 4, 4)
        mz = ti.Matrix.zero(ti.f32, 4, 4)
        if 0 < i < nx - 1:
            mx = (tensor_field.M_am[i + 1, j, k] - tensor_field.M_am[i - 1, j, k]) * inv_2dx
        if 0 < j < ny - 1:
            my = (tensor_field.M_am[i, j + 1, k] - tensor_field.M_am[i, j - 1, k]) * inv_2dx
        if 0 < k < nz - 1:
            mz = (tensor_field.M_am[i, j, k + 1] - tensor_field.M_am[i, j, k - 1]) * inv_2dx
        # --- build A[kk,li] = ⟨A_op(E_li), E_kk⟩, then symmetrize ------------
        aloc = ti.Matrix.zero(ti.f32, 10, 10)
        for li in range(10):
            el = tensor_field.sym_basis[li]
            ae = 4.0 * (
                commutator(eta_twist(commutator(el, mx)), mx)
                + commutator(eta_twist(commutator(el, my)), my)
                + commutator(eta_twist(commutator(el, mz)), mz)
            )
            for kk in range(10):
                aloc[kk, li] = (ae * tensor_field.sym_basis[kk]).sum()
        for p_ in range(10):
            for q_ in range(p_ + 1, 10):
                v = 0.5 * (aloc[p_, q_] + aloc[q_, p_])
                aloc[p_, q_] = v
                aloc[q_, p_] = v
        # --- cyclic Jacobi: A ← JᵀAJ, Q ← QJ ---------------------------------
        qloc = ti.Matrix.identity(ti.f32, 10)
        normf2 = (aloc * aloc).sum()
        off2 = normf2
        sweep = 0
        while sweep < JACOBI_SWEEPS and off2 > JTOL2_4D * normf2:
            for p_ in range(9):
                for q_ in range(p_ + 1, 10):
                    apq = aloc[p_, q_]
                    if ti.abs(apq) > 0:
                        tau = (aloc[q_, q_] - aloc[p_, p_]) / (2.0 * apq)
                        t = 1.0 / (ti.abs(tau) + ti.sqrt(1.0 + tau * tau))
                        if tau < 0:
                            t = -t
                        c = 1.0 / ti.sqrt(1.0 + t * t)
                        s = t * c
                        for r_ in range(10):  # A ← A·J (cols p,q)
                            arp = aloc[r_, p_]
                            arq = aloc[r_, q_]
                            aloc[r_, p_] = c * arp - s * arq
                            aloc[r_, q_] = s * arp + c * arq
                        for r_ in range(10):  # A ← Jᵀ·A (rows p,q)
                            apr = aloc[p_, r_]
                            aqr = aloc[q_, r_]
                            aloc[p_, r_] = c * apr - s * aqr
                            aloc[q_, r_] = s * apr + c * aqr
                        for r_ in range(10):  # Q ← Q·J
                            qrp = qloc[r_, p_]
                            qrq = qloc[r_, q_]
                            qloc[r_, p_] = c * qrp - s * qrq
                            qloc[r_, q_] = s * qrp + c * qrq
            off2 = 0.0
            for p_ in range(9):
                for q_ in range(p_ + 1, 10):
                    off2 += 2.0 * aloc[p_, q_] * aloc[p_, q_]
            sweep += 1
        # --- keep mask + projection + A⁺ solve --------------------------------
        lmax = 0.0
        for a_ in range(10):
            lmax = ti.max(lmax, ti.abs(aloc[a_, a_]))
        pm = tensor_field.P_am[i, j, k]
        pc = ti.Vector.zero(ti.f32, 10)
        for a_ in range(10):
            pc[a_] = (pm * tensor_field.sym_basis[a_]).sum()
        pproj = ti.Vector.zero(ti.f32, 10)
        cdot = ti.Vector.zero(ti.f32, 10)
        for kk in range(10):
            lamk = aloc[kk, kk]
            qp = 0.0
            for a_ in range(10):
                qp += qloc[a_, kk] * pc[a_]
            if lamk > EPS_EIG_4D * (lmax + 1e-30):
                for a_ in range(10):
                    pproj[a_] += qloc[a_, kk] * qp
                    cdot[a_] += qloc[a_, kk] * qp / lamk
        pnew = ti.Matrix.zero(ti.f32, 4, 4)
        mdnew = ti.Matrix.zero(ti.f32, 4, 4)
        for a_ in range(10):
            pnew += pproj[a_] * tensor_field.sym_basis[a_]
            mdnew += cdot[a_] * tensor_field.sym_basis[a_]
        mdnew = mdnew * tensor_field.act4d[i, j, k]
        vn = ti.sqrt((mdnew * mdnew).sum())
        if vn > VCAP_4D:
            mdnew = mdnew * (VCAP_4D / (vn + 1e-30))
        tensor_field.P_am[i, j, k] = pnew
        tensor_field.Md_am[i, j, k] = mdnew


@ti.kernel
def update_M_4d_constrained(tensor_field: ti.template(), dt_eff: ti.f32):  # type: ignore
    """M_new = M + dt_eff·Ṁ (all voxels; Ṁ = 0 outside act ⇒ border/core copy).
    The caller swaps the triple buffer as usual — render paths read M_am."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for i, j, k in ti.ndrange(nx, ny, nz):
        tensor_field.M_new_am[i, j, k] = (
            tensor_field.M_am[i, j, k] + dt_eff * tensor_field.Md_am[i, j, k]
        )


# ================================================================
# M5.24 — THE CANONICAL REGULARIZATION STACK (verified-L era port)
# ================================================================
# The production port of the dynamics formulation of record (research
# m5_theory_canonical.md § 2 row 6): canonical kinetic ½‖Ṁ‖²_F + the
# η-signed curvature energy + the universal spectral potential V4,
# integrated by the position-Verlet triple-buffer leapfrog. Ported from
# the audited M5.21.3 instrument (research/scripts/m5_21_3_a_4d.py:
# comm_eta / inner_eta / e_parts / exact grad; complex-step gated 5e-9,
# SO(1,3) invariance 1e-9, vacuum ≡ 0, 3D block-diag regression 1e-12).
#
#   eta = diag(-1, 1, 1, 1)          [A, B]_eta = A·eta·B − B·eta·A
#   <F, G>_eta = tr(eta·F·eta·Gᵀ)
#   E_u = 4 Σ_{i<j} <F_ij, F_ij>_eta ,  F_ij = [A_i, A_j]_eta ,  A_i = ∂_i M
#   V4  = w Σ_{p=1..4} (tr((M·eta)^p) − C_p)² ,  C_p = sg^p + 1 + δ^p
#   vacuum: M_vac = diag(−sg, 1, δ, 0)  (time-time NEGATIVE; V4(M_vac) = 0 exact)
#
# STENCIL: the M5.21.2b symmetrized pair (½(fwd + bwd) with exact adjoints,
# the well-posedness certificate) — NOT the 2h central diff of the older
# paths (the deep-statics checkerboard anti-recipe). Two flux passes per
# step (one per branch) reusing curv_flux_x/y/z; M_new_am doubles as the
# force accumulator so no extra field is allocated.
#
# LABEL (canonical § 2 row 6): this stack is the RUNNABLE REGULARIZATION,
# "fine for statics-adjacent dynamics questions and films" — exactly the
# launcher's use class. It is NOT the true-L evolution (the free EL IVP is
# ill-posed, § 2 row 1); never present its dynamics as the theory's own.
#
# UNITS: τ-units (c = 1) like the constrained path — the caller passes
# dt_eff = c_amrs·dt_rs. dt cap: the V4 stiff g-mode (ω ≈ 78 at g = 8,
# canonical § 4 gap ladders) bounds dt_eff < 2/78 ≈ 0.026; the certified
# full-3D working step is dt_eff = 0.005 (the 3D twin, canonical § 3; the
# axisym 0.02 is margin-critical in 3D f32 — the M5.24 selftest measured
# NaN there), enforced by the launcher's compute_timestep cap, NOT here.

W1_SPECTRAL = 7.24023879e-4  # the locked potential weight w (WSCALE, canonical § 4)


@ti.func
def eta_left(m):  # type: ignore
    """eta·M for eta = diag(−1,1,1,1): row 0 negated."""
    out = m
    for c_ in ti.static(range(4)):
        out[0, c_] = -m[0, c_]
    return out


@ti.func
def eta_right(m):  # type: ignore
    """M·eta for eta = diag(−1,1,1,1): column 0 negated."""
    out = m
    for r_ in ti.static(range(4)):
        out[r_, 0] = -m[r_, 0]
    return out


@ti.func
def comm_eta44(a, b):  # type: ignore
    """[A, B]_eta = A·eta·B − B·eta·A (the M5.18 verified bracket)."""
    return eta_right(a) @ b - eta_right(b) @ a


@ti.func
def inner_eta44(f, g_mat):  # type: ignore
    """<F, G>_eta = tr(eta·F·eta·Gᵀ) = Σ_ab (eta·F·eta)_ab · G_ab."""
    return (eta_left(eta_right(f)) * g_mat).sum()


@ti.func
def v4_of(m, sg: ti.f32, delta: ti.f32, w1: ti.f32):  # type: ignore
    """The universal spectral potential V4 = w·Σ_{p=1..4} (tr((M·eta)^p) − C_p)².

    C_p = sg^p + 1 + δ^p; exact zero on the covariant vacuum diag(−sg, 1, δ, 0)."""
    me = eta_right(m)
    p1 = me
    p2 = me @ me
    p3 = p2 @ me
    p4 = p3 @ me
    t1 = p1.trace()
    t2 = p2.trace()
    t3 = p3.trace()
    t4 = p4.trace()
    c1 = sg + 1.0 + delta
    c2 = sg * sg + 1.0 + delta * delta
    c3 = sg * sg * sg + 1.0 + delta * delta * delta
    c4 = sg * sg * sg * sg + 1.0 + delta * delta * delta * delta
    return w1 * (
        (t1 - c1) ** 2 + (t2 - c2) ** 2 + (t3 - c3) ** 2 + (t4 - c4) ** 2
    )


@ti.func
def dv4_of(m, sg: ti.f32, delta: ti.f32, w1: ti.f32):  # type: ignore
    """∂V4/∂M (symmetrized): Σ_p 2w(t_p − C_p)·p·(eta·(M·eta)^(p−1))ᵀ — the exact
    reference gradient (m5_21_3_a_4d.grad, complex-step gated)."""
    me = eta_right(m)
    p0 = ti.Matrix.identity(ti.f32, 4)
    p1 = me
    p2 = me @ me
    p3 = p2 @ me
    t1 = p1.trace()
    t2 = p2.trace()
    t3 = (p2 @ me).trace()
    t4 = (p3 @ me).trace()
    c1 = sg + 1.0 + delta
    c2 = sg * sg + 1.0 + delta * delta
    c3 = sg * sg * sg + 1.0 + delta * delta * delta
    c4 = sg * sg * sg * sg + 1.0 + delta * delta * delta * delta
    gv = (
        (2.0 * w1 * (t1 - c1) * 1.0) * eta_left(p0).transpose()
        + (2.0 * w1 * (t2 - c2) * 2.0) * eta_left(p1).transpose()
        + (2.0 * w1 * (t3 - c3) * 3.0) * eta_left(p2).transpose()
        + (2.0 * w1 * (t4 - c4) * 4.0) * eta_left(p3).transpose()
    )
    return 0.5 * (gv + gv.transpose())


@ti.kernel
def compute_eta_flux(tensor_field: ti.template(), branch: ti.i32, dx_eta: ti.f32):  # type: ignore
    """Per-cell η-curvature flux dA_i = ∂E_u/∂A_i for one stencil branch
    (0 = fwd, 1 = bwd) → curv_flux_x/y/z. From the reference grad():

        F_ij = [A_i, A_j]_eta ,  WF = 8·eta·F·eta
        dA_i += WF·(A_j·eta) − (eta·A_j)·WF     (A symmetric)
        dA_j += (eta·A_i)·WF − WF·(A_i·eta)

    A_i = one-sided ∂_i M (fwd/bwd), zero where the branch neighbor is
    outside the grid (the reference d1 edge convention). Valid at every
    cell; the evolve kernels apply the exact adjoint gather.

    dx_eta (M5.24 round 2): the grid spacing in the CANONICAL unit system.
    Pass the research grid unit (config ETA_DX, e.g. 1.5 = the m5_21_2b/3 h)
    to run the launcher as a research-twin: the physical dx_am (~15.6 at 64³
    over 1 fm) weakens the curvature term ~100× against the dimensionless V4
    and distorts the dynamics balance. Default wiring falls back to dx_am."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_dx = 1.0 / dx_eta
    for i, j, k in ti.ndrange(nx, ny, nz):
        m0 = tensor_field.M_am[i, j, k]
        ax_ = ti.Matrix.zero(ti.f32, 4, 4)
        ay_ = ti.Matrix.zero(ti.f32, 4, 4)
        az_ = ti.Matrix.zero(ti.f32, 4, 4)
        if branch == 0:  # fwd
            if i <= nx - 2:
                ax_ = (tensor_field.M_am[i + 1, j, k] - m0) * inv_dx
            if j <= ny - 2:
                ay_ = (tensor_field.M_am[i, j + 1, k] - m0) * inv_dx
            if k <= nz - 2:
                az_ = (tensor_field.M_am[i, j, k + 1] - m0) * inv_dx
        else:  # bwd
            if i >= 1:
                ax_ = (m0 - tensor_field.M_am[i - 1, j, k]) * inv_dx
            if j >= 1:
                ay_ = (m0 - tensor_field.M_am[i, j - 1, k]) * inv_dx
            if k >= 1:
                az_ = (m0 - tensor_field.M_am[i, j, k - 1]) * inv_dx
        fxy = comm_eta44(ax_, ay_)
        fxz = comm_eta44(ax_, az_)
        fyz = comm_eta44(ay_, az_)
        wxy = 8.0 * eta_left(eta_right(fxy))
        wxz = 8.0 * eta_left(eta_right(fxz))
        wyz = 8.0 * eta_left(eta_right(fyz))
        tensor_field.curv_flux_x[i, j, k] = (
            wxy @ eta_right(ay_) - eta_left(ay_) @ wxy
            + wxz @ eta_right(az_) - eta_left(az_) @ wxz
        )
        tensor_field.curv_flux_y[i, j, k] = (
            eta_left(ax_) @ wxy - wxy @ eta_right(ax_)
            + wyz @ eta_right(az_) - eta_left(az_) @ wyz
        )
        tensor_field.curv_flux_z[i, j, k] = (
            eta_left(ax_) @ wxz - wxz @ eta_right(ax_)
            + eta_left(ay_) @ wyz - wyz @ eta_right(ay_)
        )


@ti.kernel
def evolve_M_eta_start(tensor_field: ti.template(), dt_eff: ti.f32, dx_eta: ti.f32):  # type: ignore
    """Leapfrog step 1/2 (fwd branch): M_new = 2M − M_prev − dt²·½·adjᶠʷᵈ(flux).

    adjᶠʷᵈ(g)[x] = (g[x−e_i] − g[x])/h per axis (exact adjoint of the fwd
    difference; flux is zero on the invalid edge planes so the formula holds
    on the whole interior). Boundary shell untouched (Dirichlet, the
    triple-buffer BC rule — same pinning as the 3D twin, canonical § 3).
    dx_eta = the canonical-unit grid spacing (must match compute_eta_flux)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_dx = 1.0 / dx_eta
    dt2 = dt_eff * dt_eff
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        g_fwd = (
            (tensor_field.curv_flux_x[i - 1, j, k] - tensor_field.curv_flux_x[i, j, k])
            + (tensor_field.curv_flux_y[i, j - 1, k] - tensor_field.curv_flux_y[i, j, k])
            + (tensor_field.curv_flux_z[i, j, k - 1] - tensor_field.curv_flux_z[i, j, k])
        ) * inv_dx
        tensor_field.M_new_am[i, j, k] = (
            2.0 * tensor_field.M_am[i, j, k]
            - tensor_field.M_prev_am[i, j, k]
            - dt2 * 0.5 * g_fwd
        )


@ti.kernel
def evolve_M_eta_finish(
    tensor_field: ti.template(),  # type: ignore
    dt_eff: ti.f32,  # type: ignore
    dx_eta: ti.f32,  # type: ignore
    sg: ti.f32,  # type: ignore
    delta: ti.f32,  # type: ignore
    w1: ti.f32,  # type: ignore
):
    """Leapfrog step 2/2 (bwd branch + potential):
    M_new −= dt²·(½·adjᵇʷᵈ(flux) + ∂V4/∂M).  adjᵇʷᵈ(g)[x] = (g[x] − g[x+e_i])/h.
    Caller then swap_matrix_buffers() as usual. dx_eta = canonical-unit spacing."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_dx = 1.0 / dx_eta
    dt2 = dt_eff * dt_eff
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        g_bwd = (
            (tensor_field.curv_flux_x[i, j, k] - tensor_field.curv_flux_x[i + 1, j, k])
            + (tensor_field.curv_flux_y[i, j, k] - tensor_field.curv_flux_y[i, j + 1, k])
            + (tensor_field.curv_flux_z[i, j, k] - tensor_field.curv_flux_z[i, j, k + 1])
        ) * inv_dx
        force = 0.5 * g_bwd + dv4_of(tensor_field.M_am[i, j, k], sg, delta, w1)
        tensor_field.M_new_am[i, j, k] = tensor_field.M_new_am[i, j, k] - dt2 * force


# ---- M5.24 round 3: the FIRE statics relaxer (G8) + the sponge ----------
# FIRE on the SAME pinned-interior η + V4 gradient the dynamics uses (the
# two-pass kernels above), with the TIME ROW + COLUMN FROZEN (incl. [0,0]):
# the canonical § 5.2 frozen-time-row recipe, the bounded statics sector
# (the 3×3 instrument of record is exactly this sector on block-diag data).
# Parameters per the recipe: dt0 = 0.005, dt_max = 0.05, 500-iter chunks.
# The 1/w precondition of the loop recipe is NOT ported (convergence
# accelerator only; noted in the task_details).

FIRE_DT0 = 0.005  # canonical § 5.2 (frozen-time-row FIRE)
FIRE_DT_MAX = 0.05
FIRE_F_INC = 1.1
FIRE_F_DEC = 0.5
FIRE_ALPHA0 = 0.1
FIRE_F_ALPHA = 0.99
FIRE_NDELAY = 5


@ti.func
def _fire_force_at(tensor_field: ti.template(), i: ti.i32, j: ti.i32, k: ti.i32):  # type: ignore
    """FIRE force F = −∂E/∂M = M_new − M (the dt = 1, Ṁ = 0 two-pass extract),
    with the time row/col zeroed (frozen-time-row sector)."""
    f = tensor_field.M_new_am[i, j, k] - tensor_field.M_am[i, j, k]
    for a_ in ti.static(range(4)):
        f[0, a_] = 0.0
        f[a_, 0] = 0.0
    return f


@ti.kernel
def fire_partials_k(tensor_field: ti.template()):  # type: ignore
    """Per-z-slice partial sums over the interior: (Σ F·V, Σ F·F, Σ V·V) →
    fire_partials[0..2, k]. Reads the force from the M_new/M pair (run the
    two-pass with dt = 1 and M_prev = M first) and V from Md_am."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for k in ti.ndrange((1, nz - 1)):
        s_fv = 0.0
        s_ff = 0.0
        s_vv = 0.0
        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                f = _fire_force_at(tensor_field, i, j, k)
                v = tensor_field.Md_am[i, j, k]
                s_fv += (f * v).sum()
                s_ff += (f * f).sum()
                s_vv += (v * v).sum()
        tensor_field.fire_partials[0, k] = s_fv
        tensor_field.fire_partials[1, k] = s_ff
        tensor_field.fire_partials[2, k] = s_vv


@ti.kernel
def fire_update_k(
    tensor_field: ti.template(),  # type: ignore
    dtf: ti.f32,  # type: ignore
    alpha: ti.f32,  # type: ignore
    vmix: ti.f32,  # type: ignore
    vreset: ti.i32,  # type: ignore
):
    """One FIRE step on the interior (boundary pinned): velocity mixing
    V ← (1−α)V + vmix·F (vmix = α·|V|/|F| host-computed), Euler kick
    V += dtf·F, drift M += dtf·V. V lives in Md_am (free outside the
    constrained path). Time row/col frozen throughout."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        f = _fire_force_at(tensor_field, i, j, k)
        v = tensor_field.Md_am[i, j, k]
        if vreset == 1:
            v = ti.Matrix.zero(ti.f32, 4, 4)
        else:
            v = (1.0 - alpha) * v + vmix * f
        v = v + dtf * f
        for a_ in ti.static(range(4)):
            v[0, a_] = 0.0
            v[a_, 0] = 0.0
        tensor_field.Md_am[i, j, k] = v
        tensor_field.M_am[i, j, k] = tensor_field.M_am[i, j, k] + dtf * v


def fire_relax_canonical(tensor_field, dx_eta, sg, delta, w1, n_iters):
    """Host FIRE loop (single source: launcher + selftest both call this).
    Relaxes M toward a stationary state of the pinned-interior η + V4 energy
    on the frozen-time-row sector. Returns dict(f_rms0, f_rms1, iters,
    dt_final). Numpy is used ONLY for the (3, nz) partial-array sums + the
    FIRE scalars (taichi-first: all field work is in the kernels).

    Anti-recipe guards carried: dt capped at FIRE_DT_MAX = 0.05 (the frozen
    sector's cap; free-4D would need 0.02); the sym stencil is already the
    functional (no 2h deep-statics funnel, canonical § 6)."""
    import numpy as np

    tensor_field.Md_am.fill(0.0)
    dtf, alpha, npos = FIRE_DT0, FIRE_ALPHA0, 0
    f_rms0 = f_rms1 = 0.0
    n_int = max((tensor_field.nx - 2) * (tensor_field.ny - 2) * (tensor_field.nz - 2), 1)
    for it in range(n_iters):
        tensor_field.M_prev_am.copy_from(tensor_field.M_am)
        compute_eta_flux(tensor_field, 0, dx_eta)
        evolve_M_eta_start(tensor_field, 1.0, dx_eta)
        compute_eta_flux(tensor_field, 1, dx_eta)
        evolve_M_eta_finish(tensor_field, 1.0, dx_eta, sg, delta, w1)
        fire_partials_k(tensor_field)
        p = tensor_field.fire_partials.to_numpy()
        fv, ff, vv = float(p[0].sum()), float(p[1].sum()), float(p[2].sum())
        f_rms1 = (ff / (16.0 * n_int)) ** 0.5
        if it == 0:
            f_rms0 = f_rms1
        vreset = 0
        if fv > 0.0:
            npos += 1
            if npos > FIRE_NDELAY:
                dtf = min(dtf * FIRE_F_INC, FIRE_DT_MAX)
                alpha *= FIRE_F_ALPHA
        else:
            npos = 0
            dtf *= FIRE_F_DEC
            alpha = FIRE_ALPHA0
            vreset = 1
        vmix = alpha * (vv / ff) ** 0.5 if (ff > 0.0 and vreset == 0) else 0.0
        fire_update_k(tensor_field, dtf, alpha, vmix, vreset)
    # static finalize: all three buffers equal (Ṁ = 0 for a subsequent evolve)
    tensor_field.M_prev_am.copy_from(tensor_field.M_am)
    tensor_field.M_new_am.copy_from(tensor_field.M_am)
    return {"f_rms0": f_rms0, "f_rms1": f_rms1, "iters": n_iters, "dt_final": dtf}


@ti.kernel
def apply_eta_sponge(
    tensor_field: ti.template(),  # type: ignore
    dt_eff: ti.f32,  # type: ignore
    gamma_max: ti.f32,  # type: ignore
    width_cells: ti.f32,  # type: ignore
):
    """Boundary sponge for the canonical path (canonical § 5.3 pattern:
    quadratic ramp, γmax): damped-Verlet correction of the already-forced
    M_new:  M_new ← (M_new + a·M_prev)/(1 + a),  a = ½·γ(x)·dt,
    γ(x) = γmax·q², q = 1 − d/width (d = cell distance to the nearest face).
    γmax = 0 is an exact no-op. Kills the pinned-box reflection re-agitation
    (the t ≈ 12 caveat, canonical § 3) for long live runs; the research
    exact-dissipation LEDGER is not ported (demo device, noted honestly)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        d = ti.cast(i, ti.f32)
        d = ti.min(d, ti.cast(nx - 1 - i, ti.f32))
        d = ti.min(d, ti.cast(j, ti.f32))
        d = ti.min(d, ti.cast(ny - 1 - j, ti.f32))
        d = ti.min(d, ti.cast(k, ti.f32))
        d = ti.min(d, ti.cast(nz - 1 - k, ti.f32))
        if d < width_cells:
            q = 1.0 - d / width_cells
            a = 0.5 * gamma_max * q * q * dt_eff
            tensor_field.M_new_am[i, j, k] = (
                tensor_field.M_new_am[i, j, k] + a * tensor_field.M_prev_am[i, j, k]
            ) * (1.0 / (1.0 + a))


@ti.kernel
def flip_time_axis(tensor_field: ti.template()):  # type: ignore
    """Negate the time-time entry M[0,0] on all three M buffers: converts the
    block-diag seeder embed (+g, the plain-commutator era convention) to the
    covariant vacuum diag(−g, 1, δ, 0) of the verified L (canonical § 1:
    diag(+g, …) is NOT a vacuum — V4 ≈ 759 per cell at g = 8). Run once
    post-seed on the canonical integrator path; block-diag seeds only (the
    dressed seed's boost-mixed time row is a different, era-specific object)."""
    for i, j, k in tensor_field.M_am:
        tensor_field.M_am[i, j, k][0, 0] = -tensor_field.M_am[i, j, k][0, 0]
        tensor_field.M_prev_am[i, j, k][0, 0] = -tensor_field.M_prev_am[i, j, k][0, 0]
        tensor_field.M_new_am[i, j, k][0, 0] = -tensor_field.M_new_am[i, j, k][0, 0]


# ================================================================
# M5.23.1 — FIXED-J ISOROTATION (the M5.21.9 machinery, production port)
# ================================================================
# The constraint-carried ZBW clock of the two-stack consensus (research
# m5_21_9_note.md §§ 1, 5; canonical § 2): free minimization gives ω* = 0 or
# −∞, never a finite clock — the electron carries its rotation as a CONSERVED
# charge J on the CONJUGATION-tangent clock flow (the standard isorotating-
# soliton construction, Coleman → Radu-Volkov):
#
#   a0  = w(r)·[W, M] / ‖·‖_F(global) ,  W = rotation generator about the
#         local leading SPATIAL eigenvector (the director);
#         w(r) = exp(−(r/renv)⁴)  (the research clock envelope)
#   kin(M; a0) = h³ Σ_branches wt Σ_i 4·⟨[a0, A_i]_η, [a0, A_i]_η⟩_η
#         (the clock inertia; E_kin(ω) = ω²·kin — reference kin_of)
#   J   = 2·kin·ω ,  ω* = J/(2·kin) ,  E(J) = E_stat + J²/(4·kin)
#         (the Legendre family; dE/dJ = ω* closes at ~1%, m5_21_9_note § 7)
#   SET-J: Ṁ(0) = ω*·a0  →  M_prev ← M − dt_eff·ω*·a0
#         (position-Verlet implicit-velocity init; O(dt) offset from the
#         reference kick-drift-kick leap(), absorbed by observable-level gates)
#
# CONVENTION PINS (load-bearing, never change silently):
#   - kin convention = the CONJUGATION tangent [W, M] (kin = 0.1206 on the
#     certified state, M5.21.3-adopted; canonical § 6 anti-recipe: the
#     antisymmetric gen_catalog probe flow, kin = 0.297, is NOT physical
#     inertia — it exits the symmetric configuration space).
#   - any ABSOLUTE J / ħ/2 / g statement must use the PHYSICAL-RATE
#     convention (m5_21_5_note § 5): the flow-parameter J carried here is
#     NOT the physical angular momentum (ratio ~771 at the ω = 0.2 rung).
#   - a0 SIGN GAUGE: the leading eigenvector is apolar (±v give opposite
#     local rotations). The research flow inherits numpy eigh's per-voxel
#     signs; production pins a DETERMINISTIC RADIAL gauge (v·r̂ ≥ 0, fixed
#     axis tie-breaks at r ≈ 0). All quadratic reads (kin, E, J·J) are
#     per-voxel-sign invariant; the selftest quantifies the gauge overlap.
#   - the J readout projects Ṁ on the SYMMETRIC global rotation flows
#     [G_k, M] (normalized) — not gen_catalog's antisymmetric rot flows,
#     whose projection tracked a sym4-projected-out spectator (the disclosed
#     M5.21.9 § 5 Mt(0) issue).
#
# Reuses Md_am (a0 scratch — free outside the constrained path and FIRE; a
# FIRE call clobbers it, so SET-J must follow any relax, as the launcher
# flow does) and fire_partials (per-slice reduction rows, the M5.0h
# contention pattern) — no new fields. Selftest:
# research/scripts/m5_23_1_fixedj_engine_selftest.py (S1-S5).


@ti.kernel
def clock_flow_k(tensor_field: ti.template(), dx_eta: ti.f32, renv: ti.f32):  # type: ignore
    """The UNNORMALIZED conjugation-tangent clock flow a0_raw = w(r)·[W, M]
    per voxel → Md_am (full grid; the host normalizes to unit global
    Frobenius norm via md_norm_partials_k + scale_md_k).

    W = spatial rotation generator about the local leading eigenvector of
    the SPATIAL 3×3 block (production Cardano solver), radial sign gauge:
    v·r̂ ≥ 0 (tie-breaks z > y > x at the center voxel). [W, M] with W
    antisymmetric and M symmetric is symmetric: the flow stays in the
    configuration space (the conjugation convention of record)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    cx = 0.5 * ti.cast(nx - 1, ti.f32)
    cy = 0.5 * ti.cast(ny - 1, ti.f32)
    cz = 0.5 * ti.cast(nz - 1, ti.f32)
    for i, j, k in ti.ndrange(nx, ny, nz):
        m4 = tensor_field.M_am[i, j, k]
        m3 = ti.Matrix(
            [
                [m4[1, 1], m4[1, 2], m4[1, 3]],
                [m4[2, 1], m4[2, 2], m4[2, 3]],
                [m4[3, 1], m4[3, 2], m4[3, 3]],
            ]
        )
        v, _lam = principal_director(m3)
        # radial sign gauge (deterministic; apolar ±v ambiguity pinned)
        rx = ti.cast(i, ti.f32) - cx
        ry = ti.cast(j, ti.f32) - cy
        rz = ti.cast(k, ti.f32) - cz
        s = v[0] * rx + v[1] * ry + v[2] * rz
        if s == 0.0:
            s = v[2]
        if s == 0.0:
            s = v[1]
        if s == 0.0:
            s = v[0]
        if s < 0.0:
            v = -v
        w4 = ti.Matrix.zero(ti.f32, 4, 4)
        w4[1, 2] = -v[2]
        w4[1, 3] = v[1]
        w4[2, 1] = v[2]
        w4[2, 3] = -v[0]
        w4[3, 1] = -v[1]
        w4[3, 2] = v[0]
        r_env = ti.sqrt(rx * rx + ry * ry + rz * rz) * dx_eta
        q = (r_env / renv) ** 4
        w_env = ti.exp(-q)
        tensor_field.Md_am[i, j, k] = w_env * (w4 @ m4 - m4 @ w4)


@ti.kernel
def md_norm_partials_k(tensor_field: ti.template()):  # type: ignore
    """Per-z-slice Σ Md·Md → fire_partials[0, k] (full grid: the reference
    normalizes a0 over ALL cells). Host: ‖a0‖_F = sqrt(Σ_k partials)."""
    nz = tensor_field.nz
    for k in range(nz):
        s = 0.0
        for i in range(tensor_field.nx):
            for j in range(tensor_field.ny):
                a = tensor_field.Md_am[i, j, k]
                s += (a * a).sum()
        tensor_field.fire_partials[0, k] = s


@ti.kernel
def scale_md_k(tensor_field: ti.template(), s: ti.f32):  # type: ignore
    for i, j, k in tensor_field.Md_am:
        tensor_field.Md_am[i, j, k] = s * tensor_field.Md_am[i, j, k]


@ti.kernel
def kin_partials_k(tensor_field: ti.template(), dx_eta: ti.f32):  # type: ignore
    """Per-z-slice partial sums of the clock-inertia integrand →
    fire_partials[1, k]:  Σ_branches ½ Σ_i 4·⟨[a0, A_i]_η, [a0, A_i]_η⟩_η
    with a0 = Md_am (normalized), A_i the one-sided ∂_i M per branch (the
    reference d1 edge convention: zero where the neighbor is outside).
    Host: kin = dx_eta³ · Σ_k partials (reference kin_of, sym stencil)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_dx = 1.0 / dx_eta
    for k in range(nz):
        s = 0.0
        for i in range(nx):
            for j in range(ny):
                m0 = tensor_field.M_am[i, j, k]
                a0 = tensor_field.Md_am[i, j, k]
                for br in ti.static(range(2)):
                    ax_ = ti.Matrix.zero(ti.f32, 4, 4)
                    ay_ = ti.Matrix.zero(ti.f32, 4, 4)
                    az_ = ti.Matrix.zero(ti.f32, 4, 4)
                    if ti.static(br == 0):  # fwd
                        if i <= nx - 2:
                            ax_ = (tensor_field.M_am[i + 1, j, k] - m0) * inv_dx
                        if j <= ny - 2:
                            ay_ = (tensor_field.M_am[i, j + 1, k] - m0) * inv_dx
                        if k <= nz - 2:
                            az_ = (tensor_field.M_am[i, j, k + 1] - m0) * inv_dx
                    else:  # bwd
                        if i >= 1:
                            ax_ = (m0 - tensor_field.M_am[i - 1, j, k]) * inv_dx
                        if j >= 1:
                            ay_ = (m0 - tensor_field.M_am[i, j - 1, k]) * inv_dx
                        if k >= 1:
                            az_ = (m0 - tensor_field.M_am[i, j, k - 1]) * inv_dx
                    fx = comm_eta44(a0, ax_)
                    fy = comm_eta44(a0, ay_)
                    fz = comm_eta44(a0, az_)
                    s += (
                        0.5
                        * 4.0
                        * (inner_eta44(fx, fx) + inner_eta44(fy, fy) + inner_eta44(fz, fz))
                    )
        tensor_field.fire_partials[1, k] = s


@ti.kernel
def isorotation_kick_k(
    tensor_field: ti.template(),  # type: ignore
    dt_eff: ti.f32,  # type: ignore
    om: ti.f32,  # type: ignore
    shell: ti.i32,  # type: ignore
):
    """The SET-J velocity init: M_prev ← M − dt_eff·ω*·a0 (a0 = Md_am,
    normalized) on cells at least `shell` cells from every face (production
    Dirichlet shell = 1; the research free mask = 2 at h = 1.5). Encodes
    Ṁ(0) = ω*·a0 in the triple buffer — the next evolve step must run at the
    SAME dt_eff (position-Verlet: the implied velocity scales with dt)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    for i, j, k in ti.ndrange(nx, ny, nz):
        d = ti.min(
            ti.min(i, nx - 1 - i),
            ti.min(ti.min(j, ny - 1 - j), ti.min(k, nz - 1 - k)),
        )
        if d >= shell:
            tensor_field.M_prev_am[i, j, k] = (
                tensor_field.M_am[i, j, k] - dt_eff * om * tensor_field.Md_am[i, j, k]
            )


@ti.kernel
def j_partials_k(tensor_field: ti.template(), axis: ti.i32, dt_eff: ti.f32):  # type: ignore
    """Per-z-slice partials of the J projection on ONE symmetric global
    rotation flow a_rot = [G_axis, M] (axis 0/1/2 = x/y/z; G antisym, flow
    symmetric — the physical-read convention, NOT the probe flow):
    fire_partials[0, k] = Σ Ṁ·a_rot,  fire_partials[1, k] = Σ a_rot·a_rot,
    Ṁ = (M − M_prev)/dt_eff.  Host: J_axis = Σ₀ / sqrt(Σ₁)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_dt = 1.0 / dt_eff
    a = 1
    b = 2
    if axis == 0:  # rotation about x: (y, z) plane
        a = 2
        b = 3
    elif axis == 1:  # about y: (z, x) plane
        a = 3
        b = 1
    else:  # about z: (x, y) plane
        a = 1
        b = 2
    for k in range(nz):
        s_dot = 0.0
        s_nrm = 0.0
        for i in range(nx):
            for j in range(ny):
                m0 = tensor_field.M_am[i, j, k]
                g4 = ti.Matrix.zero(ti.f32, 4, 4)
                g4[a, b] = -1.0
                g4[b, a] = 1.0
                arot = g4 @ m0 - m0 @ g4
                mt = (m0 - tensor_field.M_prev_am[i, j, k]) * inv_dt
                s_dot += (mt * arot).sum()
                s_nrm += (arot * arot).sum()
        tensor_field.fire_partials[0, k] = s_dot
        tensor_field.fire_partials[1, k] = s_nrm


def compute_clock_flow(tensor_field, dx_eta, renv):
    """Host: build the normalized conjugation-tangent clock flow a0 into
    Md_am. Returns the pre-normalization global Frobenius norm (diagnostic;
    ~field scale). Numpy only for the (3, nz) partial sums (taichi-first)."""
    import numpy as np

    clock_flow_k(tensor_field, dx_eta, renv)
    md_norm_partials_k(tensor_field)
    p = tensor_field.fire_partials.to_numpy()
    nrm = float(np.sqrt(max(p[0].sum(), 0.0)))
    scale_md_k(tensor_field, 1.0 / max(nrm, 1e-30))
    return nrm


def kin_canonical(tensor_field, dx_eta):
    """Host: the clock inertia kin(M; a0) with a0 = the current Md_am
    (normalized). Reference: m5_21_3_a_4d.kin_of (sym stencil, h³-weighted)."""
    import numpy as np

    kin_partials_k(tensor_field, dx_eta)
    p = tensor_field.fire_partials.to_numpy()
    return float(dx_eta**3 * p[1].sum())


def set_fixed_j(tensor_field, dx_eta, renv, om_target, dt_eff, shell=1):
    """The launcher SET-J step (RELAX → SET-J → EVOLVE): build a0, measure
    kin, carry J = 2·kin·ω_target, kick Ṁ(0) = ω*·a0 into the triple buffer.
    Returns dict(kin, J, om_star, a0_raw_norm). Run AFTER any FIRE relax
    (FIRE clobbers Md_am); evolve must continue at the same dt_eff."""
    nrm = compute_clock_flow(tensor_field, dx_eta, renv)
    kin = kin_canonical(tensor_field, dx_eta)
    if kin <= 0.0:
        # the indefinite-channel guard: a nonpositive clock inertia means the
        # state left the positive-kin sector — no honest isorotation to set
        return {"kin": kin, "J": 0.0, "om_star": 0.0, "a0_raw_norm": nrm, "set": False}
    j_charge = 2.0 * kin * om_target
    om_star = j_charge / (2.0 * kin)
    isorotation_kick_k(tensor_field, dt_eff, om_star, shell)
    return {"kin": kin, "J": j_charge, "om_star": om_star, "a0_raw_norm": nrm, "set": True}


@ti.kernel
def j_self_partials_k(tensor_field: ti.template(), dt_eff: ti.f32):  # type: ignore
    """Per-z-slice partials of the CARRIED isorotation charge projection:
    fire_partials[2, k] = Σ Ṁ·a0 with a0 = Md_am (normalized, rebuilt on the
    CURRENT state by the host), Ṁ = (M − M_prev)/dt_eff. At SET-J time this
    projection equals ω* exactly (Ṁ = ω*·a0, unit norm) — the meaningful
    hold observable: the GLOBAL rotation flows [G_k, M] near-cancel on
    hedgehog-family states (measured ~1e-4·ω*, selftest S4 info)."""
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz
    inv_dt = 1.0 / dt_eff
    for k in range(nz):
        s = 0.0
        for i in range(nx):
            for j in range(ny):
                mt = (tensor_field.M_am[i, j, k] - tensor_field.M_prev_am[i, j, k]) * inv_dt
                s += (mt * tensor_field.Md_am[i, j, k]).sum()
        tensor_field.fire_partials[2, k] = s


def read_carried_j(tensor_field, dx_eta, renv, dt_eff):
    """Host: the carried isorotation charge J_self = ⟨Ṁ, a0(M)⟩ with a0
    REBUILT on the current state (clobbers the Md_am scratch — same rule as
    set_fixed_j). Equals ω* at SET-J time by construction; its drift is the
    honest hold/decay read. Diagnostics rate."""
    import numpy as np

    compute_clock_flow(tensor_field, dx_eta, renv)
    j_self_partials_k(tensor_field, dt_eff)
    p = tensor_field.fire_partials.to_numpy()
    return float(p[2].sum())


def read_isorotation(tensor_field, dt_eff):
    """Host: the J-vector readout — Ṁ projected on the three normalized
    SYMMETRIC global rotation flows [G_k, M]. Returns dict(Jx, Jy, Jz, Jmag,
    phi_xy). Diagnostics-rate only (3 kernel passes + host sums)."""
    import math

    import numpy as np

    out = []
    for axis in range(3):
        j_partials_k(tensor_field, axis, dt_eff)
        p = tensor_field.fire_partials.to_numpy()
        dot, nrm2 = float(p[0].sum()), float(p[1].sum())
        out.append(dot / max(math.sqrt(max(nrm2, 0.0)), 1e-30))
    jx, jy, jz = out
    return {
        "Jx": jx,
        "Jy": jy,
        "Jz": jz,
        "Jmag": float(np.sqrt(jx * jx + jy * jy + jz * jz)),
        "phi_xy": float(math.atan2(jy, jx)),
    }


# ================================================================
# GRADIENT-DESCENT RELAXATION (M5.1 task 6)
# ================================================================
# Direct port of Exp 2's relax() inner loop to Taichi:
#     dn = ∇²n − (n·∇²n) n         (tangent projection onto unit sphere)
#     n_new = n + τ · dn           (gradient-descent step)
#     n_new ← n_new / |n_new|       (unit-length constraint)
#     n_new[defect_core] = ±ẑ      (soft Dirichlet pin)
#
# Stability: τ < dx²/(2·dim·K) by the standard heat-equation CFL condition.
# For dim=3 in 3D the bound is τ < dx²/(6·K). Use ~50% of the bound for
# headroom; the launcher computes this from tensor_field.dx_am at runtime.
#
# Output buffer convention: writes director_nhat_new. The caller swaps it into
# director_nhat, then rebuild_M_from_director writes all three M buffers equal so
# subsequent Evolve PDE sees Ṁ = 0 (no spurious time-derivative from the relaxation
# update — relaxation is a STATIC operation, not a dynamics step). This is different
# from swap_matrix_buffers, which rotates the M buffers for time evolution.
#
# Boundary handling: director_nhat_new boundary voxels are copied from director_nhat
# (no update) — preserves the fixed-value Dirichlet BC (e.g. ẑ at the universe edge).
#
# Core pin: soft Dirichlet at the SINGLE closest voxel to each defect center.
# Without this, gradient descent can numerically dissolve the defect on a
# discrete grid (topology is not strictly preserved by discretization).
#
# FORWARD-LINK: this kernel is mathematically the γ → ∞ limit of a damped
# wave equation `∂²_t ψ = c²∇²ψ − γ·∂_t ψ`. Although the relaxation use is
# purely numerical (no physics interpretation in M5.1), the same primitives
# — Laplacian stencil, tangent projection, soft-core pin — could be reused
# in future damped-dynamics kernels where γ becomes a PHYSICAL damping
# coefficient (radiation loss, phonon coupling).


@ti.kernel
def relax_director_step(
    tensor_field: ti.template(),  # type: ignore
    tau: ti.f32,  # type: ignore
    pin_centers: ti.template(),  # type: ignore
    pin_signs: ti.template(),  # type: ignore
    n_defects: ti.i32,  # type: ignore
):
    """
    One gradient-descent step on the Frank elastic energy `H_F = (K/2)·|∇n̂|²`.

    Reads:  tensor_field.director_nhat
    Writes: tensor_field.director_nhat_new  (caller copies → director_nhat)

    M5.4: repointed from the retiring ψ buffers onto the matrix-substrate director
    (director_nhat = principal eigenvector of M). The math is byte-identical to the
    M5.1 relaxation — tangent-projected gradient descent on (K/2)|∇n̂|² + soft-core
    pin — so the Coulomb result carries over. Caller rebuilds M from the relaxed
    director via rebuild_M_from_director.

    Args:
        tensor_field: TensorField with the director field in director_nhat
        tau: step size (must satisfy τ < dx²/(6·K) for stability)
        pin_centers: ti.field shape (n_defects, 3) i32 — defect centers in voxel coords
        pin_signs: ti.field shape (n_defects,) i32 — ±1 per defect (pin direction)
        n_defects: number of active defects in pin_centers/signs

    Note: K_frank is implicit in `tau`'s CFL derivation — the kernel itself
    doesn't multiply by K because gradient descent on `(K/2)·|∇n̂|²` gives
    `∂n/∂τ = K·∇²n` and the K can be absorbed into the step size (τ_eff = K·τ).
    For K=1 (M5.1 default) this is τ_eff = τ.
    """
    nx, ny, nz = tensor_field.nx, tensor_field.ny, tensor_field.nz

    # ---- Interior: tangent-projected gradient step + renormalization ----
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        n = tensor_field.director_nhat[i, j, k]
        lap = compute_laplacian_director(tensor_field, i, j, k)
        n_dot_lap = n.dot(lap)
        dn = lap - n_dot_lap * n  # tangent projection (remove component along n)
        n_new = n + tau * dn
        norm = n_new.norm() + 1e-12
        tensor_field.director_nhat_new[i, j, k] = n_new / norm

    # ---- Boundary: preserve Dirichlet BC by copying from current director ----
    for j, k in ti.ndrange(ny, nz):
        tensor_field.director_nhat_new[0, j, k] = tensor_field.director_nhat[0, j, k]
        tensor_field.director_nhat_new[nx - 1, j, k] = tensor_field.director_nhat[nx - 1, j, k]
    for i, k in ti.ndrange(nx, nz):
        tensor_field.director_nhat_new[i, 0, k] = tensor_field.director_nhat[i, 0, k]
        tensor_field.director_nhat_new[i, ny - 1, k] = tensor_field.director_nhat[i, ny - 1, k]
    for i, j in ti.ndrange(nx, ny):
        tensor_field.director_nhat_new[i, j, 0] = tensor_field.director_nhat[i, j, 0]
        tensor_field.director_nhat_new[i, j, nz - 1] = tensor_field.director_nhat[i, j, nz - 1]

    # ---- Pin defect cores (soft Dirichlet at single closest voxel) ----
    for d in range(n_defects):
        ci = pin_centers[d, 0]
        cj = pin_centers[d, 1]
        ck = pin_centers[d, 2]
        sgn = ti.cast(pin_signs[d], ti.f32)
        tensor_field.director_nhat_new[ci, cj, ck] = ti.Vector([0.0, 0.0, sgn])
