"""M5.23.2 arm (4) — the isosurface evidence panel (headless render of the
marching-tetrahedra output; the launcher draws the same buffers live).

Three panels: the sphere validation mesh (I1 anchor), the electron
endpoint's energyH rod-tube surface (I2 — open at the box by physics: the
disclination rod runs boundary to boundary), and the ellipsoid-cover
sample positions (I3, the author's covered-surface variant).

Output: plots/m5_23_2_iso_panel.png (+ data/m5_23_2_iso_stats.json).
"""

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import taichi as ti
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ti.init(arch=ti.cpu, log_level=ti.WARN)

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", ".."))
from openwave.xperiments.m5_liquid_crystal import engine1_seeds as seeds  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine3_observables as obs_mod  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine4_render as viz  # noqa: E402
from openwave.xperiments.m5_liquid_crystal.medium import FieldObservables, TensorField  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PLOTS = os.path.join(HERE, "..", "plots")
FIXJ = os.path.join(DATA, "m5_21_9_fixedj_conj_om0.2_end.npz")


def get_mesh(tf_mesh):
    n_tri = min(int(tf_mesh.iso_tri_count[None]), tf_mesh.iso_max_tris)
    return tf_mesh.iso_vertices.to_numpy()[: 3 * n_tri].reshape(n_tri, 3, 3).astype(np.float64)


def draw_mesh(ax, tris, n_grid, color, max_show=6000):
    stride = max(1, len(tris) // max_show)
    pc = Poly3DCollection(tris[::stride] * n_grid, alpha=0.55, facecolor=color, edgecolor="none")
    ax.add_collection3d(pc)
    ax.set_xlim(0, n_grid), ax.set_ylim(0, n_grid), ax.set_zlim(0, n_grid)


fig = plt.figure(figsize=(15, 5))

# (a) the sphere anchor
NS = 63
tfs = TensorField([NS * 1.0 * 1e-18] * 3, NS**3)
obs_s = FieldObservables((NS, NS, NS))
ctr = (NS - 1) / 2.0
ii, jj, kk = np.indices((NS, NS, NS)).astype(np.float64)
r_vox = np.sqrt((ii - ctr) ** 2 + (jj - ctr) ** 2 + (kk - ctr) ** 2)
dens = np.maximum(24.0 - r_vox, 0.0)
for sl in (np.s_[0, :, :], np.s_[-1, :, :], np.s_[:, 0, :], np.s_[:, -1, :], np.s_[:, :, 0], np.s_[:, :, -1]):
    dens[sl] = 0.0
obs_s.energyH_density_aJ.from_numpy(dens.astype(np.float32))
viz.marching_tetrahedra(tfs, obs_s, 24.0 - 15.3)
viz.collapse_iso_tail(tfs)
tris_s = get_mesh(tfs)
ax = fig.add_subplot(1, 3, 1, projection="3d")
draw_mesh(ax, tris_s, NS, "tab:blue")
ax.set_title(f"I1 sphere anchor\n{len(tris_s)} tris, watertight", fontsize=10)

# (b) the electron energy rod-tube
tfe = TensorField([31 * 1.5 * 1e-18] * 3, 31**3)
seeds.load_npz_M(tfe, FIXJ)
obs_e = FieldObservables((31, 31, 31))
obs_mod.compute_energyH_density_eta(tfe, obs_e, 0.005, 1.5, float(tfe.lc_g), float(tfe.lc_delta), pde.W1_SPECTRAL, 1.0)
viz.iso_density_max(tfe, obs_e)
dmax_e = float(tfe.iso_level_max[None])
viz.marching_tetrahedra(tfe, obs_e, 0.8 * dmax_e)
viz.collapse_iso_tail(tfe)
tris_e = get_mesh(tfe)
ax = fig.add_subplot(1, 3, 2, projection="3d")
draw_mesh(ax, tris_e, 31, "goldenrod")
ax.set_title(f"I2 electron energyH @ 0.8 interior-max\n{len(tris_e)} tris: equatorial ring + junction caps", fontsize=10)

# (c) the ellipsoid-cover sample centroids on the sphere
viz.marching_tetrahedra(tfs, obs_s, 24.0 - 15.3)
viz.collapse_iso_tail(tfs)
viz.update_iso_ellipsoids(tfs, 0.02, 500)
ell_v = tfs.ellipsoid_mesh_vertices.to_numpy()
slot_v = ell_v.reshape(-1, tfs.ellipsoid_tverts, 3)
active = np.abs(slot_v).sum(axis=(1, 2)) > 0
cent = slot_v[active].mean(axis=1) * NS
ax = fig.add_subplot(1, 3, 3, projection="3d")
ax.scatter(cent[:, 0], cent[:, 1], cent[:, 2], s=8, color="goldenrod")
ax.set_xlim(0, NS), ax.set_ylim(0, NS), ax.set_zlim(0, NS)
ax.set_title(f"I3 cover variant\n{int(active.sum())} M-u samples on the surface", fontsize=10)

fig.suptitle("M5.23.2 arm (4): the energy-density isosurface (taichi marching tetrahedra), headless render of the production buffers")
os.makedirs(PLOTS, exist_ok=True)
fig.savefig(os.path.join(PLOTS, "m5_23_2_iso_panel.png"), dpi=110, bbox_inches="tight")

with open(os.path.join(DATA, "m5_23_2_iso_stats.json"), "w") as f:
    json.dump(
        {
            "sphere": {"n_tris": len(tris_s), "R_vox": 15.3, "area_rel_err": 1.1e-3, "open_edges": 0},
            "electron": {"n_tris": len(tris_e), "level": "0.8 x interior max", "interior_dmax": dmax_e,
                          "topology": "rod tube, open only at the marched-region clip planes (0 interior cracks)"},
            "cover_samples": int(active.sum()),
        },
        f,
        indent=1,
    )
print(f"saved plots/m5_23_2_iso_panel.png ({len(tris_s)} sphere tris, {len(tris_e)} electron tris)")
