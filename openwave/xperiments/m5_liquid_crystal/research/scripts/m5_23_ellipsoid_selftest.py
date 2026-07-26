"""
M5.23 — headless selftest of the VIZ.5 ellipsoid visualization (mesh-only).

The feature renders one M.u eigen-ellipsoid surface per 3D angle on an S² shell
around each defect center (full-3D view; mesh-only since 2026-07-19, the line
and wireframe variants were retired after the live A/B). This gate builds a
small TensorField, seeds a biaxial hedgehog, fills the shell center, runs
update_ellipsoid_mesh and machine-checks (design + decisions:
../tasks/m5_23_task_details.md):

  1. TEMPLATE   — unit-sphere UV template (poles + rings), all vertices unit
  2. INDICES    — every face indexes inside its own slot block, 3 distinct verts
  3. GEOMETRY   — per-direction ellipsoid CENTROID sits exactly on the R-shell
                  and the direction set matches the Fibonacci closed form
                  (the template centroid is 0, so the linear map preserves it)
  4. MESH MAP   — vertex = p + (size/2)*((M_sp + floor*I) @ u) reproduced
                  against a numpy reference at sampled directions
  5. PHYSICS    — on a hedgehog shell the director is radial: |u . n| ~ 1 for
                  the overwhelming majority of directions
  6. HYGIENE    — slots beyond (n_centers, n_active) collapse to the origin,
                  no NaNs, the density knob re-runs cleanly at the 642 ceiling
  7. MULTI-CENTER — the dipole case: two shells, each centered on its defect
  8. RODS (Stage D) — rod-line samples on the axis at uniform heights, the M.u
                  map reproduced, and the PHYSICS read: the rod core ESCAPES TO
                  UNIAXIAL (lam2 - lam3 collapses vs the shell's ~delta split)
  9. ROD RINGS  — ring centroids on circles around the cord at the fixed
                  heights; toggle isolation (rods-only blanks rings and back)

Also renders a Lambert-shaded matplotlib preview of the ellipsoid shell
(research/plots/m5_23_mesh_selftest.png).

USAGE (repo root):
    python -m openwave.xperiments.m5_liquid_crystal.research.scripts.m5_23_ellipsoid_selftest
"""

import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import taichi as ti
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ti.init(arch=ti.cpu, log_level=ti.WARN)

import openwave.xperiments.m5_liquid_crystal.engine1_seeds as seeds  # noqa: E402
import openwave.xperiments.m5_liquid_crystal.engine2_pde as pde  # noqa: E402
import openwave.xperiments.m5_liquid_crystal.engine4_render as viz  # noqa: E402
from openwave.xperiments.m5_liquid_crystal.medium import TensorField  # noqa: E402

# ----------------------------------------------------------------
# Setup: small grid + biaxial hedgehog at the exact center
# ----------------------------------------------------------------
UNIVERSE = [1.0e-15, 1.0e-15, 1.0e-15]  # arbitrary (geometry only)
TARGET_VOXELS = 40**3

tf = TensorField(UNIVERSE, TARGET_VOXELS)
nx, ny, nz = tf.nx, tf.ny, tf.nz
cx, cy, cz = (nx - 1) // 2, (ny - 1) // 2, (nz - 1) // 2
r0_vox = 0.10 * tf.max_grid_size
rhoc_vox = 3.0
BIAX_DELTA = 0.3
seeds.seed_biaxial_hedgehog_M(tf, cx, cy, cz, r0_vox, rhoc_vox, BIAX_DELTA)
pde.eigen_decompose(tf)

tf.ellipsoid_centers[0] = [float(cx), float(cy), float(cz)]
tf.ellipsoid_n_centers[None] = 1

R_FRAC = 0.18  # small enough that the 2.4R rod + 2.24R rings stay in-grid
R_VOX = R_FRAC * tf.max_grid_size
SIZE = 0.036
N_ACTIVE = 162

max_dim = float(tf.max_grid_size)
center_n = (np.array([cx, cy, cz], np.float64) + 0.5) / max_dim  # normalized center
fails = []
total = [0]


def check(name, ok, detail=""):
    total[0] += 1
    tag = "PASS" if ok else "FAIL"
    print(f"[{tag}] {name}" + (f" ({detail})" if detail else ""))
    if not ok:
        fails.append(name)


# Fibonacci reference directions (must match the kernel closed form)
d_idx = np.arange(N_ACTIVE, dtype=np.float64)
zf = 1.0 - 2.0 * (d_idx + 0.5) / N_ACTIVE
rho = np.sqrt(np.maximum(1.0 - zf * zf, 0.0))
phi = d_idx * np.pi * (3.0 - np.sqrt(5.0))
u_ref = np.stack([rho * np.cos(phi), rho * np.sin(phi), zf], axis=1)

# ---- 1. TEMPLATE ----
tmpl = tf.ellipsoid_template.to_numpy().astype(np.float64)
check(
    "template: unit sphere + poles, centroid at origin",
    np.allclose(np.linalg.norm(tmpl, axis=1), 1.0, atol=1e-6)
    and np.allclose(tmpl[0], [0, 0, 1])
    and np.allclose(tmpl[-1], [0, 0, -1])
    and np.allclose(tmpl.mean(axis=0), 0.0, atol=1e-6),
)

# ---- 2. INDICES ----
tverts, tfaces = tf.ellipsoid_tverts, tf.ellipsoid_tfaces
idxs = tf.ellipsoid_mesh_indices.to_numpy().reshape(-1, 3)
n_slots = tf.ellipsoid_max_centers * tf.ellipsoid_max_dirs
slot_of_face = np.repeat(np.arange(n_slots), tfaces)
in_block = (idxs // tverts == slot_of_face[:, None]).all()
distinct = (
    (idxs[:, 0] != idxs[:, 1]) & (idxs[:, 1] != idxs[:, 2]) & (idxs[:, 0] != idxs[:, 2])
).all()
check(
    "indices: in-slot blocks, 3 distinct verts per face",
    bool(in_block and distinct and idxs.min() >= 0 and idxs.max() < n_slots * tverts),
)

# ---- 3. GEOMETRY (centroids on the shell, Fibonacci directions) ----
viz.update_ellipsoid_mesh(tf, R_VOX, SIZE, N_ACTIVE)
mesh_v = tf.ellipsoid_mesh_vertices.to_numpy().astype(np.float64)
cents = mesh_v[: N_ACTIVE * tverts].reshape(N_ACTIVE, tverts, 3).mean(axis=1)
r_norm = np.linalg.norm(cents - center_n, axis=1) * max_dim
check(
    "geometry: ellipsoid centroids on the R-shell",
    np.allclose(r_norm, R_VOX, atol=1e-3),
    f"centroid radius {r_norm.min():.4f}..{r_norm.max():.4f} vox vs R={R_VOX:.4f}",
)
u_meas = (cents - center_n) / np.linalg.norm(cents - center_n, axis=1, keepdims=True)
check(
    "geometry: Fibonacci direction set",
    np.allclose(u_meas, u_ref, atol=1e-4),
    f"max direction mismatch {np.abs(u_meas - u_ref).max():.2e}",
)

# ---- 4. MESH MAP (numpy reproduction) ----
snap = np.clip(
    np.round(np.array([cx, cy, cz], np.float64) + R_VOX * u_ref).astype(int),
    0,
    [nx - 1, ny - 1, nz - 1],
)
ii, jj, kk = snap[:, 0], snap[:, 1], snap[:, 2]
M_np = tf.M_am.to_numpy().astype(np.float64)
floor = viz._ELLIPSOID_MESH_FLOOR
ok_mesh = True
for d in range(0, N_ACTIVE, 17):  # sample every 17th direction
    p = (np.array([cx, cy, cz], np.float64) + R_VOX * u_ref[d] + 0.5) / max_dim
    m_sp = M_np[ii[d], jj[d], kk[d]][1:4, 1:4]
    expect = p + 0.5 * SIZE * (tmpl @ (m_sp + floor * np.eye(3)).T)
    got = mesh_v[d * tverts : (d + 1) * tverts]
    if not np.allclose(got, expect, atol=1e-5):
        ok_mesh = False
        break
check("mesh map: vertex = p + s/2*(M_sp + floor*I)@u for sampled dirs", ok_mesh)

# ---- 5. PHYSICS ----
nhat = tf.director_nhat.to_numpy().astype(np.float64)
n_vox = nhat[ii, jj, kk]
radial = np.abs(np.sum(u_ref * n_vox, axis=1))  # |u . n|, 1 = radial director
frac_radial = float(np.mean(radial > 0.9))
check(
    "physics: hedgehog shell director is radial",
    frac_radial >= 0.90,
    f"{100 * frac_radial:.1f}% of directions with |u.n| > 0.9 "
    f"(disclination-axis neighborhoods exempt)",
)

# ---- 6. HYGIENE ----
check(
    "hygiene: slots beyond n_active collapsed + finite",
    np.all(mesh_v[N_ACTIVE * tverts :] == 0.0) and np.all(np.isfinite(mesh_v)),
)
N2 = 642  # density re-run at the ceiling
viz.update_ellipsoid_mesh(tf, R_VOX, SIZE, N2)
v2 = tf.ellipsoid_mesh_vertices.to_numpy().astype(np.float64)
c2 = v2[: N2 * tverts].reshape(N2, tverts, 3).mean(axis=1)
r2 = np.linalg.norm(c2 - center_n, axis=1) * max_dim
check(
    "hygiene: density knob re-run (n=642) stays on shell",
    np.allclose(r2, R_VOX, atol=1e-3) and np.all(v2[N2 * tverts :] == 0.0),
)

# ---- 7. MULTI-CENTER (the dipole case: one shell per defect) ----
c2v = np.array([cx - 8, cy, cz], np.float64)
tf.ellipsoid_centers[1] = [float(c2v[0]), float(c2v[1]), float(c2v[2])]
tf.ellipsoid_n_centers[None] = 2
viz.update_ellipsoid_mesh(tf, R_VOX, SIZE, N_ACTIVE)
v3 = tf.ellipsoid_mesh_vertices.to_numpy().astype(np.float64)
blk = tf.ellipsoid_max_dirs * tverts  # vertex block per center
ca = v3[: N_ACTIVE * tverts].reshape(N_ACTIVE, tverts, 3).mean(axis=1)
cb = v3[blk : blk + N_ACTIVE * tverts].reshape(N_ACTIVE, tverts, 3).mean(axis=1)
c2n = (c2v + 0.5) / max_dim
ra = np.linalg.norm(ca - center_n, axis=1) * max_dim
rb = np.linalg.norm(cb - c2n, axis=1) * max_dim
check(
    "multi-center: two shells, each on its own center",
    np.allclose(ra, R_VOX, atol=1e-3) and np.allclose(rb, R_VOX, atol=1e-3),
)
tf.ellipsoid_n_centers[None] = 1  # restore for the rod section + preview plot

# ---- 8. RODS (Stage D: the disclination-rod line samples) ----
rod_n = tf.ellipsoid_rod_n
ring_az = min(max(N_ACTIVE // 12, 8), tf.ellipsoid_ring_az)  # mirrors the launcher's Count-driven density
rod_slots = tf.ellipsoid_rod_slots
RING_R = max(3.0, 0.04 * max_dim)  # mirrors the launcher's fixed v1 ring radius
HL = viz.ELLIPSOID_ROD_SPAN * R_VOX  # rod half-length: rods protrude beyond the shell
viz.update_rod_ellipsoids(tf, HL, RING_R, SIZE, ring_az, 1, 1)
rv = tf.ellipsoid_rod_vertices.to_numpy().astype(np.float64)
rod_cents = rv[: rod_n * tverts].reshape(rod_n, tverts, 3).mean(axis=1)
h_ref = np.linspace(-HL, HL, rod_n)
exp_rod = (
    np.stack([np.full(rod_n, cx), np.full(rod_n, cy), cz + h_ref], axis=1) + 0.5
) / max_dim
check(
    "rods: line samples on the axis at uniform heights",
    np.allclose(rod_cents, exp_rod, atol=1e-4),
    f"max centroid err {np.abs(rod_cents - exp_rod).max():.2e}",
)
# M.u reproduction at a mid-rod slot
d_rod = rod_n // 3
pr = exp_rod[d_rod]
snap_rod = np.clip(
    np.round([cx, cy, cz + h_ref[d_rod]]).astype(int), 0, [nx - 1, ny - 1, nz - 1]
)
m_sp = M_np[snap_rod[0], snap_rod[1], snap_rod[2]][1:4, 1:4]
expect = pr + 0.5 * SIZE * (tmpl @ (m_sp + floor * np.eye(3)).T)
check(
    "rods: M.u map reproduced at a rod slot",
    np.allclose(rv[d_rod * tverts : (d_rod + 1) * tverts], expect, atol=1e-5),
)
# PHYSICS: the rod core ESCAPES TO UNIAXIAL — the biaxial split of the two
# minor eigenvalues (lam2 - lam3 ~ delta in the bulk) collapses on the axis
# where the seed's disclination smoothstep switches the biaxial part off.
# (First try tested lam1 - lam2 "melt" and FAILED: the axis is MORE uniaxial,
# not lower-amplitude — the correct core signature, and the natural detection
# criterion for the future disclination-line tracer.)
evals_np = tf.eigenvalues.to_numpy().astype(np.float64)
rod_idx = np.clip(np.round(cz + h_ref).astype(int), 0, nz - 1)
uni_rod = (evals_np[cx, cy, rod_idx, 1] - evals_np[cx, cy, rod_idx, 2]).mean()
uni_shell = (evals_np[ii, jj, kk, 1] - evals_np[ii, jj, kk, 2]).mean()
check(
    "rods physics: rod core escapes to uniaxial (lam2-lam3 collapses)",
    uni_rod < 0.5 * uni_shell,
    f"rod mean lam2-lam3 = {uni_rod:.3f} vs shell {uni_shell:.3f}",
)

# ---- 9. ROD RINGS (per 2D angle around the cord) ----
n_ring_slots = tf.ellipsoid_ring_count * ring_az
ring_block = rv[rod_n * tverts : (rod_n + n_ring_slots) * tverts].reshape(-1, tverts, 3)
ring_cents = ring_block.mean(axis=1)
h_rows_ref = np.array(
    [-1.28, -1.60, -1.92, -2.24, 1.28, 1.60, 1.92, 2.24]
)  # kernel row order (sign by j//4), SHELL-radius units: one row gap off the shell
ok_rings = True
for jr in range(tf.ellipsoid_ring_count):
    seg = ring_cents[jr * ring_az : (jr + 1) * ring_az]
    ang = 2.0 * np.pi * np.arange(ring_az) / ring_az
    expc = (
        np.stack(
            [
                cx + RING_R * np.cos(ang),
                cy + RING_R * np.sin(ang),
                np.full(ring_az, cz + h_rows_ref[jr] * R_VOX),
            ],
            axis=1,
        )
        + 0.5
    ) / max_dim
    if not np.allclose(seg, expc, atol=1e-4):
        ok_rings = False
        break
check("rings: centroids on circles around the cord at the fixed heights", ok_rings)
viz.update_rod_ellipsoids(tf, HL, RING_R, SIZE, ring_az, 1, 0)  # rods only
rv2 = tf.ellipsoid_rod_vertices.to_numpy()
viz.update_rod_ellipsoids(tf, HL, RING_R, SIZE, ring_az, 0, 1)  # rings only
rv3 = tf.ellipsoid_rod_vertices.to_numpy()
check(
    "rings: toggle isolation (rods-only blanks rings and back)",
    np.all(rv2[rod_n * tverts : rod_slots * tverts] == 0.0)
    and bool(np.any(rv2[: rod_n * tverts] != 0.0))  # rods present in rods-only
    and np.all(rv3[: rod_n * tverts] == 0.0)  # rods blank in rings-only
    and bool(np.any(rv3[rod_n * tverts : rod_slots * tverts] != 0.0))
    and np.all(np.isfinite(rv))
    and np.all(rv[rod_slots * tverts :] == 0.0),  # inactive centers collapsed
)


# ----------------------------------------------------------------
# Preview plot (Agg) — Lambert-shaded eigen-ellipsoid shell
# ----------------------------------------------------------------
viz.update_ellipsoid_mesh(tf, R_VOX, SIZE, N_ACTIVE)
mesh_v = tf.ellipsoid_mesh_vertices.to_numpy().astype(np.float64)
tris = mesh_v[idxs[: N_ACTIVE * tfaces]]  # (n_faces, 3, 3) active block
normals = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
normals /= np.linalg.norm(normals, axis=1, keepdims=True) + 1e-15
light = np.array([0.4, 0.3, 0.85])
lam = 0.35 + 0.65 * np.abs(normals @ (light / np.linalg.norm(light)))
face_rgb = np.clip(lam[:, None] * np.array([0.29, 0.64, 0.87]), 0, 1)
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection="3d")
ax.add_collection3d(Poly3DCollection(tris, facecolors=face_rgb, edgecolors="none"))
ax.scatter(*center_n, color="k", s=18)
lim = (center_n[0] - 1.3 * R_FRAC, center_n[0] + 1.3 * R_FRAC)
ax.set_xlim(lim), ax.set_ylim(lim), ax.set_zlim(lim)
ax.set_box_aspect((1, 1, 1))
ax.set_title(
    f"M5.23 VIZ.5 ellipsoid viz (M.u eigen-ellipsoids): biaxial hedgehog, "
    f"{N_ACTIVE} ellipsoids, R={R_FRAC:.2f}"
)
plot_dir = os.path.join(os.path.dirname(__file__), "..", "plots")
os.makedirs(plot_dir, exist_ok=True)
out = os.path.join(plot_dir, "m5_23_mesh_selftest.png")
fig.savefig(out, dpi=110, bbox_inches="tight")
print(f"[plot] {out}")


def _shaded(tris_arr, base_rgb):
    nrm = np.cross(tris_arr[:, 1] - tris_arr[:, 0], tris_arr[:, 2] - tris_arr[:, 0])
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-15
    lt = np.array([0.4, 0.3, 0.85])
    lm = 0.35 + 0.65 * np.abs(nrm @ (lt / np.linalg.norm(lt)))
    return np.clip(lm[:, None] * np.array(base_rgb), 0, 1)


# Combined Stage D preview: shell + disclination rods + rod rings (the
# electron-clock composition: hedgehog shell + the bipolar rod pair)
viz.update_ellipsoid_mesh(tf, R_VOX, SIZE, N_ACTIVE)
sh_v = tf.ellipsoid_mesh_vertices.to_numpy().astype(np.float64)
viz.update_rod_ellipsoids(tf, HL, RING_R, SIZE, ring_az, 1, 1)
rd_v = tf.ellipsoid_rod_vertices.to_numpy().astype(np.float64)
rod_idxs = tf.ellipsoid_rod_indices.to_numpy().reshape(-1, 3)
sh_tris = sh_v[idxs[: N_ACTIVE * tfaces]]
rod_tris_a = rd_v[rod_idxs[: rod_n * tfaces]]  # rod-line block (cyan)
rod_tris_b = rd_v[rod_idxs[rod_n * tfaces : rod_slots * tfaces]]  # ring block
fig3 = plt.figure(figsize=(9, 9))
ax3 = fig3.add_subplot(111, projection="3d")
ax3.add_collection3d(
    Poly3DCollection(sh_tris, facecolors=_shaded(sh_tris, [0.29, 0.64, 0.87]),
                     edgecolors="none", alpha=0.35)
)
ax3.add_collection3d(
    Poly3DCollection(rod_tris_a, facecolors=_shaded(rod_tris_a, [0.0, 0.9, 0.9]),
                     edgecolors="none")
)
ax3.add_collection3d(
    Poly3DCollection(rod_tris_b, facecolors=_shaded(rod_tris_b, [0.29, 0.64, 0.87]),
                     edgecolors="none")
)
ax3.scatter(*center_n, color="k", s=18)
lim3 = (center_n[0] - 1.15 * viz.ELLIPSOID_ROD_SPAN * R_FRAC, center_n[0] + 1.15 * viz.ELLIPSOID_ROD_SPAN * R_FRAC)
ax3.set_xlim(lim3), ax3.set_ylim(lim3), ax3.set_zlim(lim3)
ax3.set_box_aspect((1, 1, 1))
ax3.set_title(
    f"M5.23 VIZ.5 Stage D: shell (translucent) + disclination rods (cyan) + rod rings\n"
    f"biaxial hedgehog, R={R_FRAC:.2f}, rod span +-2.4R, ring rows at 1.28R..2.24R (0.32R gaps)"
)
out3 = os.path.join(plot_dir, "m5_23_rods_selftest.png")
fig3.savefig(out3, dpi=110, bbox_inches="tight")
print(f"[plot] {out3}")

n_pass = total[0] - len(fails)
print(
    f"\nSELFTEST: {n_pass}/{total[0]} checks passed"
    + (f"; FAILED: {fails}" if fails else ": ALL GREEN")
)
sys.exit(1 if fails else 0)
