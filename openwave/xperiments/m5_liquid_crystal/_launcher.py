"""
Xperiment Launcher

Unified launcher for tensor-field xperiments featuring:
- UI-based xperiment selection and switching
- Single source of truth for rendering and UI code
- Xperiment-specific parameters in /xparameters directory
"""

import webbrowser
import importlib
import sys
import os
import time
from pathlib import Path

import taichi as ti
import numpy as np

from openwave.common import colormap, constants
from openwave.i_o import flux_mesh, render, video

import openwave.xperiments.m5_liquid_crystal.medium as medium
import openwave.xperiments.m5_liquid_crystal.engine1_seeds as seeds
import openwave.xperiments.m5_liquid_crystal.engine2_pde as pde
import openwave.xperiments.m5_liquid_crystal.engine3_observables as observables
import openwave.xperiments.m5_liquid_crystal.engine4_render as viz
import openwave.xperiments.m5_liquid_crystal.instrumentation as instrument

# ================================================================
# XPERIMENT PARAMETERS MANAGEMENT
# ================================================================


class XperimentManager:
    """Manages loading and switching between xperiment parameters."""

    def __init__(self):
        self.available_xperiments = self._discover_xperiments()
        self.xperiment_display_names = {}  # Cache display names from meta
        self.current_xperiment = None
        self.current_xparameters = None

    def _discover_xperiments(self):
        """Discover all available xperiment parameters in /xparameters directory."""
        parameters_dir = Path(__file__).parent / "xparameters"

        if not parameters_dir.exists():
            return []

        xperiment_files = [
            file.stem for file in parameters_dir.glob("*.py") if file.name != "__init__.py"
        ]

        return sorted(xperiment_files)

    def load_xperiment(self, xperiment_name):
        """Load xperiment parameters by name.

        Args:
            xperiment_name: Parameter file name without .py extension

        Returns:
            dict: Parameters dictionary or None if loading fails
        """
        try:
            module_path = f"openwave.xperiments.m5_liquid_crystal.xparameters.{xperiment_name}"
            parameters_module = importlib.import_module(module_path)
            importlib.reload(parameters_module)  # Reload for fresh parameters

            self.current_xperiment = xperiment_name
            self.current_xparameters = parameters_module.XPARAMETERS

            # Cache display name from meta
            self.xperiment_display_names[xperiment_name] = parameters_module.XPARAMETERS["meta"][
                "X_NAME"
            ]

            return self.current_xparameters

        except Exception as e:
            print(f"Error loading xperiment '{xperiment_name}': {e}")
            return None

    def get_xperiment_display_name(self, xperiment_name):
        """Get display name from parameter meta or fallback to filename conversion."""
        if xperiment_name in self.xperiment_display_names:
            return self.xperiment_display_names[xperiment_name]

        # Fallback: try to load just for the name
        try:
            module_path = f"openwave.xperiments.m5_liquid_crystal.xparameters.{xperiment_name}"
            parameters_module = importlib.import_module(module_path)
            display_name = parameters_module.XPARAMETERS["meta"]["X_NAME"]
            self.xperiment_display_names[xperiment_name] = display_name
            return display_name
        except:
            # Final fallback: convert filename
            return " ".join(word.capitalize() for word in xperiment_name.split("_"))


# ================================================================
# SIMULATION STATE
# ================================================================


class SimulationState:
    """Manages the state of the simulation."""

    def __init__(self):
        self.tensor_field = None
        self.trackers = None  # Trackers (M3/M4 wave statistics: amp, freq, EMA-RMS)
        self.observables = None  # FieldObservables (M5 derived scalars: energyH, energyF)
        self.c_amrs = 0.0
        self.dt_rs = 0.0
        self.cfl_factor = 0.0
        # M5.5.4 — Eq.13 LdG V(M)=a·Tr(M²)−b·Tr(M³)+c·(Tr(M²))² couplings for evolve_M /
        # compute_energyH_density_M. Default 0 → free curvature dynamics (V off); the exact
        # Λ=(1,δ,0)-producing values are Q7 (Duda open). Expose as xparameters later.
        self.ldg_a = 0.0
        self.ldg_b = 0.0
        self.ldg_c = 0.0
        self.ldg_v0 = 0.0  # vacuum potential V_M(D), subtracted in the energyH display (M5.6.5c)
        self.elapsed_t_rs = 0.0
        self.clock_start_time = time.time()
        self.frame = 1
        # Field aggregates — populated by sample_avg_trackers; physical units
        # except for the energy attrs (see notes below).
        self.amp_global_rms = 0.0  # m
        self.freq_global_avg = 0.0  # Hz
        self.wavelength_global_avg = 0.0  # m, derived from freq_global_avg
        # Energy attrs — currently in SCALED units (per-voxel (am/rs)² ×
        # ATTOJOULE), not true J. The kernel that populates them
        # (compute_energyH_density) doesn't yet apply the physical scaling
        # factor (ρ_medium × voxel_volume × INTERNAL_ENERGY_TO_AJ). Useful as
        # *relative* observables for trends and force-gradient sourcing
        # (F = −∇E only cares about gradients). REVIEW IN M5.2 when nonlinear
        # V(ψ) couplings land alongside the physical-scaling factor.
        self.energyH_global_avg = 0.0  # per-voxel mean in scaled "(rel.)" units
        self.energyH_total = 0.0  # mean × voxel_count, scaled "(rel.)" units
        # Frank elastic density (M5.1 task 5) — populated by sample_avg_trackers.
        self.energyF_global_avg = 0.0  # per-voxel mean in scaled "(rel.)" units
        # Resolution metric: voxels-per-wavelength (xperiment-driven).
        # 0.0 means "no reference λ declared" → dashboard shows n/a
        self.wave_res = 0.0

        # Current xperiment parameters
        self.X_NAME = ""
        self.CAM_INIT = [2.00, 1.50, 1.75]
        self.UNIVERSE_SIZE = []
        self.TARGET_VOXELS = 1e8

        # UI control variables
        self.SHOW_AXIS = False
        self.TICK_SPACING = 0.25
        self.SHOW_GRID = False
        self.SHOW_EDGES = False
        self.VIZ_STRIDE = 4  # sampling stride for directors AND granules (every Nth voxel)
        self.SHOW_GLYPHS = 0  # director-glyph planes (0=off, 1=XY, 2=+XZ, 3=all)
        self.GLYPH_VECTOR = (
            0  # glyph state — 0=director only, 1=director+delta, 2=E field, 3=B field (∇×n̂)
        )
        self.GLYPH_SIZE = (
            0  # glyph shaft — 0=unit (director), 1=field magnitude (E:|∇·n̂|, B:‖∇×n̂‖)
        )
        self.GLYPH_COLOR = (
            1  # glyph color — 0=single (COLOR_MEDIUM, see far field), 1=field gradient
        )
        # VIZ.5 (M5.23) — the "ellipsoid" viz: one M·u eigen-ellipsoid surface
        # per 3D angle on an S² shell around each defect center (full-3D view).
        # Its OWN feature namespace (not a GLYPH_VECTOR state); while ON it
        # suppresses the plane glyphs (no visual pollution). Mesh-only since
        # 2026-07-19 (the line-glyph and wireframe variants were retired after
        # the live A/B: mesh won).
        self.SHOW_ELLIPSOID = False
        self.ELLIPSOID_RADIUS = 0.1  # shell radius / rod half-length, fraction of max grid dim
        self.ELLIPSOID_COUNT = 299  # shell ellipsoid count (density; not-too-dense spec)
        self.ELLIPSOID_SIZE = 0.025  # ellipsoid base size, normalized (scaled by lambda)
        # Stage D sub-toggles: the S² shell, the disclination-rod line samples
        # (the escaped-core pair along the spin axis), and the per-2D-angle
        # rings around the rod cord. Combinable (the electron-clock figure
        # shows shell + rods together).
        self.ELLIPSOID_SHELL = True
        self.ELLIPSOID_RODS = False
        self.ELLIPSOID_RODRINGS = True
        # M5.23.2 arm (4) — the energy-density ISOSURFACE (canonical path
        # only: V4 + curvature are exactly zero on the covariant vacuum, so
        # level sets need no background subtraction). ISO_LEVEL = fraction
        # of the per-frame interior density max. ISO_COVER = the author's
        # "surface covered with ellipsoids" variant (uses the shell pool,
        # exclusive with the S² shell view).
        self.SHOW_ISO = False
        self.ISO_LEVEL = 0.30
        self.ISO_COVER = False
        self.FLUX_MESH_PLANES = [0.5, 0.5, 0.5]
        self.SHOW_FLUX_MESH = 0
        self.WARP_MESH = 300
        self.SHOW_GRANULES = False
        self.SIM_SPEED = 1.0
        self.PAUSED = False

        # Color control variables
        self.COLOR_THEME = "OCEAN"
        self.WAVE_MENU = 1
        # VIZ.2 — WM7 EM-curl color: 0=orange ‖∇×n̂‖ magnitude (honest static default),
        # 1=bluered signed (∇×n̂)·axis → N=red/S=blue poles. Axis default ẑ (auto-axis = M5.8).
        self.CURL_COLOR = 0
        self.CURL_AXIS = ti.Vector([0.0, 0.0, 1.0], dt=ti.f32)
        # VIZ.4 — magnetic-dipole placeholder sample (M5.6.5f stage 1). When True,
        # the curl observables are overwritten each frame with an analytic dipole B
        # about DIPOLE_AXIS at DIPOLE_CENTER so the N/S coloring + B glyphs + moment
        # glyph can be validated before the real circulating B exists (M5.8). Set by
        # the _viz_sample_dipole xparameter; off for all physics runs.
        self.DIPOLE_SAMPLE = False
        self.DIPOLE_AXIS = ti.Vector([0.0, 0.0, 1.0], dt=ti.f32)
        self.DIPOLE_CENTER = [0.5, 0.5, 0.5]  # fractional grid coords
        self.DIPOLE_R0_VOX = 3.0  # core regularization radius (voxels)

        # Data Analytics & video export toggles
        self.INSTRUMENTATION = False
        self.EXPORT_VIDEO = False
        self.VIDEO_FRAMES = 24

        # Optional topology seed (M5.1+ vacuum / hedgehog xperiments)
        self.TOPOLOGY_SEED = None

        # M5.1 task 6 — gradient-descent relaxation
        # pin_centers / pin_signs / n_defects are captured at seed time so the
        # relaxation kernel can soft-pin the defect cores at ±ẑ each step
        # (prevents numerical dissolution of topology on the discrete grid).
        self.pin_centers = None  # ti.field (n_defects, 3) i32 — voxel coords
        self.pin_signs = None  # ti.field (n_defects,) i32 — ±1
        self.n_defects = 0
        # Auto-relax on seed — set by xperiment via AUTO_RELAX_STEPS in topology_seed.
        # Relaxation is a numerical ground-state-finder, not physics, so there's
        # no demo / interactive use case — every seed gets relaxed once,
        # silently, before the user sees the field.
        self.AUTO_RELAX_STEPS = 0

    def apply_xparameters(self, params):
        """Apply parameters from xperiment parameter dictionary."""
        # Meta
        self.X_NAME = params["meta"]["X_NAME"]

        # Camera
        self.CAM_INIT = params["camera"]["INITIAL_POSITION"]

        # Universe
        universe = params["universe"]
        self.UNIVERSE_SIZE = list(universe["SIZE"])
        self.TARGET_VOXELS = universe["TARGET_VOXELS"]

        # UI defaults
        ui = params["ui_defaults"]
        self.SHOW_AXIS = ui["SHOW_AXIS"]
        self.TICK_SPACING = ui["TICK_SPACING"]
        self.SHOW_GRID = ui["SHOW_GRID"]
        self.SHOW_EDGES = ui["SHOW_EDGES"]
        self.VIZ_STRIDE = ui.get("VIZ_STRIDE", 4)
        self.SHOW_GLYPHS = ui.get("SHOW_GLYPHS", 0)
        self.GLYPH_VECTOR = ui.get("GLYPH_VECTOR", 0)  # 0=director 1=+delta 2=E 3=B
        self.GLYPH_SIZE = ui.get("GLYPH_SIZE", 0)  # 0=unit, 1=field magnitude
        self.GLYPH_COLOR = ui.get("GLYPH_COLOR", 1)  # 0=single, 1=field gradient
        self.SHOW_ELLIPSOID = ui.get("SHOW_ELLIPSOID", False)  # VIZ.5 ellipsoid shell
        self.ELLIPSOID_RADIUS = ui.get("ELLIPSOID_RADIUS", 0.1)  # fraction of max grid dim
        self.ELLIPSOID_COUNT = ui.get("ELLIPSOID_COUNT", 299)  # shell glyph count
        self.ELLIPSOID_SIZE = ui.get("ELLIPSOID_SIZE", 0.025)  # ellipsoid base size
        self.ELLIPSOID_SHELL = ui.get("ELLIPSOID_SHELL", True)  # S² shell sub-toggle
        self.ELLIPSOID_RODS = ui.get("ELLIPSOID_RODS", False)  # Stage D rod-line samples
        self.SHOW_ISO = ui.get("SHOW_ISO", False)  # M5.23.2 energy isosurface
        self.ISO_LEVEL = ui.get("ISO_LEVEL", 0.30)  # fraction of density max
        self.ISO_COVER = ui.get("ISO_COVER", False)  # ellipsoid-covered variant
        self.ELLIPSOID_RODRINGS = ui.get("ELLIPSOID_RODRINGS", True)  # Stage D rod rings
        self.FLUX_MESH_PLANES = ui["FLUX_MESH_PLANES"]
        self.SHOW_FLUX_MESH = ui["SHOW_FLUX_MESH"]
        self.WARP_MESH = ui["WARP_MESH"]
        self.SHOW_GRANULES = ui["SHOW_GRANULES"]
        self.SIM_SPEED = ui.get("SIM_SPEED", 1.0)
        self.PAUSED = ui["PAUSED"]

        # Color defaults
        color = params["color_defaults"]
        self.COLOR_THEME = color["COLOR_THEME"]
        self.WAVE_MENU = color["WAVE_MENU"]
        self.CURL_COLOR = color.get("CURL_COLOR", 0)  # VIZ.2 WM7: 0=orange mag, 1=bluered N/S

        # Data Analytics & video export toggles
        diag = params["analytics"]
        self.INSTRUMENTATION = diag["INSTRUMENTATION"]
        self.EXPORT_VIDEO = diag["EXPORT_VIDEO"]
        self.VIDEO_FRAMES = diag["VIDEO_FRAMES"]

        # Optional topology seed (only present in _test_topology xperiment, M5.1+)
        self.TOPOLOGY_SEED = params.get("topology_seed", None)
        # M5.1 task 6 — auto-relaxation (optional; 0 = manual via RELAX button only)
        topo_seed = params.get("topology_seed") or {}
        self.AUTO_RELAX_STEPS = topo_seed.get("AUTO_RELAX_STEPS", 0)

        # VIZ.4 — magnetic-dipole placeholder sample (optional; off unless the
        # _viz_sample_dipole xparameter sets it). DIPOLE_AXIS doubles as CURL_AXIS
        # so the bluered N/S projection is taken along the moment direction.
        self.DIPOLE_SAMPLE = bool(topo_seed.get("DIPOLE_SAMPLE", False))
        if self.DIPOLE_SAMPLE:
            ax = topo_seed.get("DIPOLE_AXIS", [0.0, 0.0, 1.0])
            self.DIPOLE_AXIS = ti.Vector([float(ax[0]), float(ax[1]), float(ax[2])], dt=ti.f32)
            self.CURL_AXIS = ti.Vector([float(ax[0]), float(ax[1]), float(ax[2])], dt=ti.f32)
            self.DIPOLE_CENTER = topo_seed.get("CENTER", [0.5, 0.5, 0.5])
            self.DIPOLE_R0_VOX = float(topo_seed.get("DIPOLE_R0_VOX", 3.0))

    def initialize_grid(self):
        """Initialize or reinitialize the tensor-field grid."""
        self.tensor_field = medium.TensorField(
            self.UNIVERSE_SIZE,
            self.TARGET_VOXELS,
            self.FLUX_MESH_PLANES,
            viz_stride=self.VIZ_STRIDE,
        )
        self.trackers = medium.Trackers(self.tensor_field.grid_size)
        self.observables = medium.FieldObservables(self.tensor_field.grid_size)

        # Resolution metric: voxels-per-wavelength, declared by the active xperiment.
        # (M5.8 ψ-retire: the TEST_SEED seed-driven source was removed; M5.2+ particle
        # xperiments will set this from the defect Compton wavelength. 0.0 = no reference.)
        self.wave_res = 0.0

    def compute_timestep(self):
        """Compute timestep from the 3D CFL stability bound.

        CFL condition for the leapfrog wave-equation solver:
            dt ≤ dx / (c · √3)

        We size dt at 95% of the bound — sitting exactly at 1/3 puts us on
        the marginal-stability boundary, where f32 rounding alone is enough
        to cross into amplification. cfl_factor = (SIM_SPEED · 0.95)² / 3,
        so at SIM_SPEED = 1 the factor is ≈ 0.301 (well under 1/3 ≈ 0.333).

        SIM_SPEED ∈ (0, 1] is a slow-motion knob that scales the effective
        wave speed for visualization without affecting stability.
        """
        cfl_safety = 0.95
        self.c_amrs = (
            constants.WAVE_SPEED / constants.ATTOMETER * constants.RONTOSECOND * self.SIM_SPEED
        )  # am/rs (slowed)
        # Tight CFL bound against physical c (= c_amrs / SIM_SPEED), with safety
        # factor — so cfl_factor saturates at (cfl_safety² / 3) at SIM_SPEED=1.
        self.dt_rs = (
            self.tensor_field.dx_am * cfl_safety / (self.c_amrs / self.SIM_SPEED * (3**0.5))
        )  # rs
        # M5.8.2c — the 4D paths derate dt (the fast-collective-mode hierarchy,
        # 2c seed e). The derate MUST live here, not in initialize_xperiment:
        # this method runs EVERY frame (the main loop refreshes dt for the
        # observables + SIM_SPEED slider), so a one-shot derate at init is
        # silently WIPED on the next frame — the 2026-06-05 signed-GUI
        # explosion: the constrained kernel ran at the full CFL dt, 78× its
        # validated ratio (headless was bounded at the derated dt all along).
        if self.TOPOLOGY_SEED is not None and self.TOPOLOGY_SEED.get("MODE") == "dressed_hedgehog":
            self.dt_rs *= float(self.TOPOLOGY_SEED.get("DT_SCALE_4D", 0.5))
        # M5.24 — the canonical path caps dt_eff = c·dt at the certified τ-step
        # (0.005, the full-3D twin's step, canonical § 3; the axisym stack's
        # 0.02 sits at ω_stiff·dt = 1.57, too close to the leapfrog margin for
        # the 3D f32 path — measured NaN at the M5.24 selftest try 1): the CFL
        # above scales with dx but the V4 stiffness does NOT (the stiff g-mode
        # ω ≈ 78 at g = 8, § 4 gap ladders, bounds dt_eff < 2/78 ≈ 0.026 in
        # τ-units regardless of grid). Lives here, not in initialize_xperiment,
        # for the same every-frame-wipe reason as the 4D derate above.
        if self.TOPOLOGY_SEED is not None and (
            str(self.TOPOLOGY_SEED.get("INTEGRATOR_4D", "")) == "canonical"
        ):
            dt_eta_cap = float(self.TOPOLOGY_SEED.get("DT_ETA_CAP", 0.005))
            self.dt_rs = min(self.dt_rs, dt_eta_cap / self.c_amrs)
        self.cfl_factor = round((self.c_amrs * self.dt_rs / self.tensor_field.dx_am) ** 2, 7)
        # V(M) couplings for the matrix substrate. Only director-class (topology)
        # seeds use a potential, so gate on TOPOLOGY_SEED presence. (M5.8 ψ-retire:
        # the legacy V(ψ) Klein-Gordon + φ⁴ couplings m_freq_kg_rs / lambda_phi4 were
        # removed with the Vector(3) ψ engine — V is now the LdG V(M) only.)
        if self.TOPOLOGY_SEED is not None:
            # M5.6.5c — Eq.13 LdG V(M) coupling for the matrix substrate (evolve_M /
            # compute_energyH_density_M). The b=0 amplitude well V = a·Tr(M²) + c·(Tr(M²))²
            # pins Tr(M²)→s₂*=1+δ² — confines the energy seen diluting under Evolve PDE —
            # while leaving the biaxiality δ EXACTLY flat (the canonical 3-term Eq.13 has NO
            # biaxial minimum: for any (a,b,c) the nonzero eigenvalues collapse to one λ*;
            # 5a §5f). Natural cubic-balance unit is c²/dx⁴ (production sweep m5_6_5c_prod_scale:
            # confines 3.3× across k∈[0.5,25] with no blow-up, dt²-stable). OFF unless the
            # xperiment sets LDG_STIFFNESS_K > 0. b=0 is the interim — a fully biaxial-STABLE
            # vacuum needs an extra invariant in V (Q7, flagged to Duda).
            ldg_k = float(self.TOPOLOGY_SEED.get("LDG_STIFFNESS_K", 0.0))
            if ldg_k > 0.0:
                delta = float(self.TOPOLOGY_SEED.get("BIAXIAL_DELTA", 0.3))
                self.ldg_c = ldg_k * self.c_amrs**2 / self.tensor_field.dx_am**4
                self.ldg_a = -2.0 * self.ldg_c * (1.0 + delta * delta)  # min at s₂*=1+δ²
                self.ldg_b = 0.0
                # vacuum potential V_M(D=diag(1,δ,0)) — subtracted in the energyH display so
                # the negative well bottom doesn't swamp the structure (kernel docstring).
                # M5.8.1: V_M acts on the SPATIAL 3×3 block only (the time-g axis is decoupled),
                # so the vacuum traces are the spatial diag(1,δ,0) ones — no g terms.
                tr2_vac = 1.0 + delta * delta
                tr3_vac = 1.0 + delta * delta * delta
                self.ldg_v0 = self.ldg_a * tr2_vac - self.ldg_b * tr3_vac + self.ldg_c * tr2_vac**2
            else:
                self.ldg_a = self.ldg_b = self.ldg_c = self.ldg_v0 = 0.0
        else:
            # No topology seed → no potential (free curvature dynamics, V = 0).
            self.ldg_a = self.ldg_b = self.ldg_c = self.ldg_v0 = 0.0

    def reset_sim(self):
        """Reset simulation state."""
        self.tensor_field = None
        self.trackers = None
        self.c_amrs = 0.0
        self.dt_rs = 0.0
        self.cfl_factor = 0.0
        self.ldg_a = self.ldg_b = self.ldg_c = self.ldg_v0 = 0.0
        self.elapsed_t_rs = 0.0
        self.clock_start_time = time.time()
        self.frame = 1
        self.amp_global_rms = 0.0
        self.freq_global_avg = 0.0
        self.wavelength_global_avg = 0.0
        self.energyH_global_avg = 0.0
        self.energyH_total = 0.0
        self.energyF_global_avg = 0.0
        self.wave_res = 0.0
        self.initialize_grid()
        self.compute_timestep()
        initialize_xperiment(self)


# ================================================================
# UI OVERLAY WINDOWS
# ================================================================


def display_xperiment_launcher(xperiment_mgr, state):
    """Display xperiment launcher UI with selectable xperiments.

    Args:
        xperiment_mgr: XperimentManager instance
        state: SimulationState instance

    Returns:
        str or None: Selected xperiment name or None
    """
    selected_xperiment = None

    with render.gui.sub_window("XPERIMENT LAUNCHER", 0.00, 0.00, 0.14, 0.25) as sub:
        sub.text("(needs window reload)", color=colormap.LIGHT_BLUE[1])
        for xp_name in xperiment_mgr.available_xperiments:
            display_name = xperiment_mgr.get_xperiment_display_name(xp_name)
            is_current = xp_name == xperiment_mgr.current_xperiment

            if sub.checkbox(display_name, is_current) and not is_current:
                selected_xperiment = xp_name

        if sub.button("Close Launcher (esc)"):
            render.window.running = False

    return selected_xperiment


def display_controls(state):
    """Display the controls UI overlay."""
    with render.gui.sub_window("CONTROLS", 0.00, 0.25, 0.16, 0.43) as sub:
        state.SHOW_AXIS = sub.checkbox(f"Axis (ticks: {state.TICK_SPACING})", state.SHOW_AXIS)
        state.SHOW_EDGES = sub.checkbox("Sim Universe Edges", state.SHOW_EDGES)
        # VIZ.5 (M5.23) — the "ellipsoid" viz: S² shell, one eigenframe glyph per
        # 3D angle around each defect center. Its own feature (options below);
        # while ON it replaces the plane glyphs below (no visual pollution).
        state.SHOW_ELLIPSOID = sub.checkbox("Ellipsoids", state.SHOW_ELLIPSOID)
        if state.SHOW_ELLIPSOID:
            state.ELLIPSOID_RADIUS = sub.slider_float("Radius", state.ELLIPSOID_RADIUS, 0.02, 0.5)
            state.ELLIPSOID_COUNT = sub.slider_int(
                "Count", state.ELLIPSOID_COUNT, 32, medium.ELLIPSOID_MAX_DIRS
            )
            state.ELLIPSOID_SIZE = sub.slider_float("Size", state.ELLIPSOID_SIZE, 0.01, 0.15)
            state.ELLIPSOID_SHELL = sub.checkbox("Shell (1 /3D-angle)", state.ELLIPSOID_SHELL)
            state.ELLIPSOID_RODS = sub.checkbox("Rods (disclination)", state.ELLIPSOID_RODS)
        # M5.23.2 arm (4) — the energy isosurface (canonical path only)
        state.SHOW_ISO = sub.checkbox("Iso-Surface (energyH)", state.SHOW_ISO)
        if state.SHOW_ISO:
            state.ISO_LEVEL = sub.slider_float("Iso Level", state.ISO_LEVEL, 0.02, 0.95)
            state.ISO_COVER = sub.checkbox("Cover w/ Ellipsoids", state.ISO_COVER)
            state.ELLIPSOID_RODRINGS = sub.checkbox(
                "Rod Rings (1 /2D-angle)", state.ELLIPSOID_RODRINGS
            )
        state.SHOW_GLYPHS = sub.slider_int("Glyph Planes", state.SHOW_GLYPHS, 0, 3)
        # VIZ.5 exclusivity: the vector-glyph select displays UNCHECKED while the
        # ellipsoid shell owns the screen (those states are not rendering), and
        # clicking any of them hands the screen back (ellipsoid off) — the same
        # rule the render dispatch enforces (ellipsoid replaces the plane glyphs).
        if sub.checkbox("Director Vector", state.GLYPH_VECTOR == 0 and not state.SHOW_ELLIPSOID):
            state.GLYPH_VECTOR = 0
            state.SHOW_ELLIPSOID = False
        if sub.checkbox(
            "Director + Delta Vectors", state.GLYPH_VECTOR == 1 and not state.SHOW_ELLIPSOID
        ):
            state.GLYPH_VECTOR = 1
            state.SHOW_ELLIPSOID = False
        if sub.checkbox(
            "Electric Field (divergence)", state.GLYPH_VECTOR == 2 and not state.SHOW_ELLIPSOID
        ):
            state.GLYPH_VECTOR = 2
            state.SHOW_ELLIPSOID = False
        if sub.checkbox(
            "Magnetic Field (curl)", state.GLYPH_VECTOR == 3 and not state.SHOW_ELLIPSOID
        ):
            state.GLYPH_VECTOR = 3
            state.SHOW_ELLIPSOID = False
        # shaft length: unit (director, all glyphs visible) vs field magnitude (E: |∇·n̂| charge
        # density, B: ‖∇×n̂‖); color: single flat COLOR_MEDIUM (see far field) vs field gradient.
        if sub.checkbox("Glyph Size (unit/magnitude)", state.GLYPH_SIZE == 1):
            state.GLYPH_SIZE = 1
        else:
            state.GLYPH_SIZE = 0
        if sub.checkbox("Glyph Color (single/gradient)", state.GLYPH_COLOR == 1):
            state.GLYPH_COLOR = 1
        else:
            state.GLYPH_COLOR = 0
        state.SHOW_FLUX_MESH = sub.slider_int("Flux Mesh Planes", state.SHOW_FLUX_MESH, 0, 3)
        state.WARP_MESH = sub.slider_int("Warp Mesh", state.WARP_MESH, 0, 50)
        # VIZ.3 — 4-state glyph select (mutually exclusive checkboxes):
        #   0=Director Vector (principal axis n̂ only),
        #   1=Director + Delta Vectors (n̂ + delta cross-bar = biaxial ellipsoid frame),
        #   2=Electric Field (+→− barb, charge-colored), 3=Magnetic Field (∇×n̂).
        # Greek δ is spelled "Delta" — GGUI cannot render Greek glyphs.
        state.SHOW_GRANULES = sub.checkbox("Show Granule Motion", state.SHOW_GRANULES)
        state.SIM_SPEED = sub.slider_float("Speed", state.SIM_SPEED, 0.5, 1.0)
        # M5.24 round 3 — the canonical FIRE relaxer (G8): relax the seeded
        # state toward a stationary configuration of the SAME eta + V4 energy
        # before evolving (500-iter chunks per press, the § 5.2 recipe cadence).
        if getattr(state, "integrator_4d", "") == "canonical" and state.PAUSED:
            if sub.button("RELAX (FIRE canonical)"):
                relax_field_canonical(state, 500)
            # M5.23.1 — SET-J: (re)initialize the isorotation velocity on the
            # current state (run after RELAX; the next EVOLVE carries the clock)
            if sub.button("SET J (isorotate)"):
                om_btn = float(getattr(state, "fixedj_omega", 0.0)) or 0.2
                set_fixed_j_launcher(state, om_btn)
        if state.PAUSED:
            if sub.button(">> EVOLVE PDE >>"):
                state.PAUSED = False
        else:
            if sub.button("Pause"):
                state.PAUSED = True
        if sub.button("Reset Simulation"):
            state.reset_sim()


def display_wave_menu(state):
    """Display wave properties selection menu."""
    with render.gui.sub_window("WAVE MENU", 0.00, 0.75, 0.15, 0.25) as sub:
        if sub.checkbox("Deviation (Magnitude)", state.WAVE_MENU == 1):
            state.WAVE_MENU = 1
            state.tensor_field.create_flux_mesh()
        if sub.checkbox("Amp (A) (EMA RMS)", state.WAVE_MENU == 2):
            state.WAVE_MENU = 2
            state.tensor_field.create_flux_mesh()
        if sub.checkbox("Clock (omega)", state.WAVE_MENU == 3):
            state.WAVE_MENU = 3
            state.tensor_field.create_flux_mesh()
        if sub.checkbox("ENERGY (Hamiltonian)", state.WAVE_MENU == 4):
            state.WAVE_MENU = 4
            state.tensor_field.create_flux_mesh()
        if sub.checkbox("ENERGY (Frank Elastic)", state.WAVE_MENU == 5):
            state.WAVE_MENU = 5
            state.tensor_field.create_flux_mesh()
        if sub.checkbox("EM div (charge/E)", state.WAVE_MENU == 6):
            state.WAVE_MENU = 6
            state.tensor_field.create_flux_mesh()
        if sub.checkbox("EM curl (rotation/B)", state.WAVE_MENU == 7):
            state.WAVE_MENU = 7
            state.tensor_field.create_flux_mesh()
        # VIZ.2: WM7 color toggle — orange ‖∇×n̂‖ magnitude (off) vs bluered signed (∇×n̂)·ẑ
        # → N=red/S=blue poles (on). The vector-warp (fabric-twist) is always on for WM7.
        if state.WAVE_MENU == 7:
            if sub.checkbox("  > B color N/S (bluered)", state.CURL_COLOR == 1):
                state.CURL_COLOR = 1
            else:
                state.CURL_COLOR = 0
        # Display gradient palette with 2× average range for headroom (allows peak visualization)
        if state.WAVE_MENU == 1:  # Displacement on orange gradient
            render.canvas.triangles(og_palette_vertices, per_vertex_color=og_palette_colors)
            with render.gui.sub_window("displacement", 0.00, 0.69, 0.08, 0.06) as sub:
                sub.text(f"0       {state.amp_global_rms*2:.0e}m")
        if state.WAVE_MENU == 2:  # Amplitude (A) (EMA RMS) on ironbow gradient
            render.canvas.triangles(ib_palette_vertices, per_vertex_color=ib_palette_colors)
            with render.gui.sub_window("amplitude", 0.00, 0.69, 0.08, 0.06) as sub:
                sub.text(f"0       {state.amp_global_rms*2:.0e}m")
        if state.WAVE_MENU == 3:  # Clock (omega) on blueprint gradient
            render.canvas.triangles(bp_palette_vertices, per_vertex_color=bp_palette_colors)
            with render.gui.sub_window("omega", 0.00, 0.69, 0.08, 0.06) as sub:
                sub.text(f"0       {state.freq_global_avg*2:.0e}Hz")
        if state.WAVE_MENU == 4:  # Energy density (Hamiltonian) on ironbow gradient
            render.canvas.triangles(ib_palette_vertices, per_vertex_color=ib_palette_colors)
            with render.gui.sub_window("energyH", 0.00, 0.69, 0.08, 0.06) as sub:
                # Unit label "rel." — value is per-voxel mean × 4 (matches the
                # colormap range max in update_flux_mesh_values). Underlying field
                # is in scaled (am/rs)² units, not physical J/m³ — REVIEW IN M5.2
                # when physical-scaling factor (ρ × voxel_volume × correction)
                # is wired into compute_energyH_density.
                sub.text(f"0      {state.energyH_global_avg*4:.0e}rel.")
        if state.WAVE_MENU == 5:  # Frank elastic density on ironbow gradient
            render.canvas.triangles(ib_palette_vertices, per_vertex_color=ib_palette_colors)
            with render.gui.sub_window("energyF", 0.00, 0.69, 0.08, 0.06) as sub:
                # Same "rel." caveat as WAVE_MENU=4 — K_frank=1.0 dimensionless
                # until M5.6 physical elastic constants land.
                sub.text(f"0      {state.energyF_global_avg*4:.0e}rel.")
        if state.WAVE_MENU == 6:  # EM divergence ∇·n̂ on greenyellow diverging gradient (signed)
            render.canvas.triangles(gy_palette_vertices, per_vertex_color=gy_palette_colors)
            with render.gui.sub_window("div (charge/E)", 0.00, 0.69, 0.08, 0.06) as sub:
                sub.text(" -           +")
        if state.WAVE_MENU == 7:  # EM curl: orange magnitude OR bluered signed (∇×n̂)·ẑ N/S
            if state.CURL_COLOR == 1:
                render.canvas.triangles(br_palette_vertices, per_vertex_color=br_palette_colors)
                with render.gui.sub_window("B (S/N)", 0.00, 0.69, 0.08, 0.06) as sub:
                    sub.text(" S           N")
            else:
                render.canvas.triangles(og_palette_vertices, per_vertex_color=og_palette_colors)
                with render.gui.sub_window("curl (rot/B)", 0.00, 0.69, 0.08, 0.06) as sub:
                    sub.text("0           max")


def display_model_specs(state, model_bar_vertices):
    """Display OpenWave model specifications overlay."""
    render.canvas.triangles(model_bar_vertices, color=colormap.ORANGE[1])
    with render.gui.sub_window("LIQUID-CRYSTAL MODEL (M5)", 0.84, 0.01, 0.16, 0.16) as sub:
        sub.text("Medium: Indexed Voxel Grid")
        sub.text("Data-Structure: Tensor Field")
        sub.text("Coupling: Non-linear Lagrangian")
        sub.text("Propagation: PDE Solver")
        sub.text("Boundary: Dirichlet Condition")
        if sub.button("Wave Notation Guide"):
            webbrowser.open(
                "https://github.com/openwave-labs/openwave/blob/main/openwave/wave_notation.md"
            )


def display_data_dashboard(state):
    """Display simulation data dashboard."""
    clock_time = time.time() - state.clock_start_time
    sim_time_years = clock_time / (state.elapsed_t_rs * constants.RONTOSECOND or 1) / 31_536_000

    with render.gui.sub_window("DATA-DASHBOARD", 0.84, 0.17, 0.16, 0.60) as sub:
        state.INSTRUMENTATION = sub.checkbox("Instrumentation ON/OFF", state.INSTRUMENTATION)
        sub.text("--- SPACETIME ---", color=colormap.LIGHT_BLUE[1])
        sub.text(f"Medium Density: {constants.MEDIUM_DENSITY:.1e} kg/m³")
        sub.text(f"Wave Speed (c): {constants.WAVE_SPEED:.1e} m/s")

        sub.text("\n--- SIMULATION DOMAIN ---", color=colormap.LIGHT_BLUE[1])
        sub.text(f"Universe: {state.tensor_field.max_universe_edge:.1e} m")
        sub.text(f"Voxel Count: {state.tensor_field.voxel_count:,}")
        sub.text(
            f"Grid Size: {state.tensor_field.nx} x {state.tensor_field.ny} x {state.tensor_field.nz}"
        )
        sub.text(f"Voxel Edge: {state.tensor_field.dx:.2e} m")

        sub.text("\n--- RESOLUTION ---", color=colormap.LIGHT_BLUE[1])
        # wave_res is xperiment-driven; 0.0 means no reference λ declared.
        # Stable leapfrog needs ≥12 voxels/λ; <10 is undersampled.
        if state.wave_res > 0:
            sub.text(f"Wave: {state.wave_res:.1f} voxels/wave (>12)")
            if state.wave_res < 10:
                sub.text(f"*** WARNING: Undersampling! ***", color=(1.0, 0.0, 0.0))
        else:
            sub.text("Wave: n/a (no wave declared)")

        sub.text("\n--- WAVE-FIELD ---", color=colormap.LIGHT_BLUE[1])
        sub.text(f"Amplitude (avg): {state.amp_global_rms:.1e} m")
        sub.text(f"Frequency (avg): {state.freq_global_avg:.1e} Hz")
        sub.text(f"Wavelength (avg): {state.wavelength_global_avg:.1e} m")
        # Unit label "rel." — energyH_total is `mean × voxel_count × ATTOJOULE`
        # but the mean is in scaled (am/rs)² units, not actual aJ — so the value
        # isn't physically J. REVIEW IN M5.2 when physical-scaling factor lands.
        sub.text(f"Total ENERGY (H): {state.energyH_total:.1e} rel.")

        sub.text("\n--- TIME MICROSCOPE ---", color=colormap.LIGHT_BLUE[1])
        sub.text(f"Sim Steps (frames): {state.frame:,}")
        sub.text(f"Sim Time: {state.elapsed_t_rs:,.0f} rs")
        sub.text(f"Clock Time: {clock_time:.2f} s")
        sub.text(f"(1s sim time takes {sim_time_years:.0e}y)")

        sub.text("\n--- CFL TIMESTEP ---", color=colormap.LIGHT_BLUE[1])
        sub.text(f"c_amrs: {state.c_amrs:.3f} am/rs")
        sub.text(
            f"dt_rs: {state.dt_rs:.3f} rs",
            color=(1.0, 1.0, 1.0) if state.cfl_factor <= (1 / 3) else (1.0, 0.0, 0.0),
        )
        sub.text(
            f"CFL Factor: {state.cfl_factor:.3f} (target < 1/3)",
            color=(1.0, 1.0, 1.0) if state.cfl_factor <= (1 / 3) else (1.0, 0.0, 0.0),
        )


# ================================================================
# XPERIMENT RENDERING
# ================================================================


def initialize_xperiment(state):
    """Initialize color palettes, wave charger and instrumentation.

    Args:
        state: SimulationState instance with xperiment parameters
    """
    global og_palette_vertices, og_palette_colors
    global ib_palette_vertices, ib_palette_colors
    global bp_palette_vertices, bp_palette_colors
    global gy_palette_vertices, gy_palette_colors
    global br_palette_vertices, br_palette_colors
    global model_bar_vertices

    # Initialize color palette scales for gradient rendering and level indicator
    og_palette_vertices, og_palette_colors = colormap.get_palette_scale(
        colormap.orange, 0.00, 0.68, 0.079, 0.01
    )
    ib_palette_vertices, ib_palette_colors = colormap.get_palette_scale(
        colormap.ironbow, 0.00, 0.68, 0.079, 0.01
    )
    bp_palette_vertices, bp_palette_colors = colormap.get_palette_scale(
        colormap.blueprint, 0.00, 0.68, 0.079, 0.01
    )
    gy_palette_vertices, gy_palette_colors = colormap.get_palette_scale(
        colormap.greenyellow, 0.00, 0.68, 0.079, 0.01
    )
    br_palette_vertices, br_palette_colors = colormap.get_palette_scale(
        colormap.bluered, 0.00, 0.68, 0.079, 0.01
    )
    model_bar_vertices = colormap.get_model_bar_geometry(0.84, 0.00, 0.159, 0.01)

    # (M5.8 ψ-retire: the legacy TEST_SEED Gaussian/plane-wave seeding path was
    # removed with the Vector(3) ψ engine — all live xperiments use TOPOLOGY_SEED.)

    # M5.24 — reset the 4D-integrator routing before seeding: these are set
    # only by the dressed / canonical activations below, so an xperiment
    # switch away from a 4D config would otherwise keep evolving on the
    # previous integrator (a stale-state hazard the seed branches never
    # cleared).
    state.evolve_4d = False
    state.integrator_4d = ""
    # M5.23.2 — npz-loader routing flags, reset with the integrator routing
    # (same stale-state hazard): npz_covariant marks a loaded state ALREADY in
    # the covariant vacuum convention (M[0,0] = -g), so the canonical
    # activation must NOT flip it again; npz_h carries the source grid unit
    # as the ETA_DX fallback.
    state.npz_covariant = False
    state.npz_h = 0.0

    # Optional topology seed (M5.1+ vacuum / hedgehog xperiments)
    if state.TOPOLOGY_SEED is not None:
        topo = state.TOPOLOGY_SEED
        seed_mode = topo.get("MODE", "vacuum")  # "vacuum" or "hedgehog"

        if seed_mode == "vacuum":
            seeds.seed_vacuum_M(state.tensor_field, state.tensor_field.lc_delta)
            print("[M5.4] seeded matrix vacuum (M = δI + (1−δ)ẑ⊗ẑ everywhere)")

        elif seed_mode == "hedgehog":
            # Defects are specified in normalized [0,1] domain coordinates;
            # convert to voxel coords and bundle into ti.fields for the kernel.
            defects = topo["DEFECTS"]  # list of {"CENTER": [x,y,z], "SIGN": ±1}
            n_defects = len(defects)
            tf = state.tensor_field

            centers_np = np.zeros((n_defects, 3), dtype=np.int32)
            signs_np = np.zeros(n_defects, dtype=np.int32)
            for d, defect in enumerate(defects):
                cx_norm, cy_norm, cz_norm = defect["CENTER"]
                centers_np[d, 0] = int(round(cx_norm * (tf.nx - 1)))
                centers_np[d, 1] = int(round(cy_norm * (tf.ny - 1)))
                centers_np[d, 2] = int(round(cz_norm * (tf.nz - 1)))
                signs_np[d] = int(defect["SIGN"])

            centers_field = ti.field(dtype=ti.i32, shape=(n_defects, 3))
            signs_field = ti.field(dtype=ti.i32, shape=(n_defects,))
            centers_field.from_numpy(centers_np)
            signs_field.from_numpy(signs_np)

            # D/4 in voxel units — w_vac falloff radius for the vacuum blend.
            # Default from Exp 2 is one quarter of the largest grid extent.
            domain_quarter_fraction = topo.get("DOMAIN_QUARTER_FRACTION", 0.25)
            domain_quarter_voxels = float(domain_quarter_fraction * max(tf.nx, tf.ny, tf.nz))

            seeds.seed_hedgehog_M(
                tf, centers_field, signs_field, domain_quarter_voxels, n_defects, tf.lc_delta
            )
            # Stash pin info for relaxation kernel (M5.1 task 6)
            state.pin_centers = centers_field
            state.pin_signs = signs_field
            state.n_defects = n_defects
            print(
                f"[M5.1] seeded {n_defects} hedgehog defect(s); "
                f"D/4 = {domain_quarter_voxels:.1f} voxels"
            )

        elif seed_mode == "biaxial_hedgehog":
            # M5.6.5a: a single BIAXIAL hedgehog M = O·D(s(r))·Oᵀ, D = diag(1, δ, 0)
            # (frame O=[r̂|e_Θ|e_Φ] + disclination smoothstep + radial eigenvalue melt).
            tf = state.tensor_field
            center = topo.get("CENTER", [0.5, 0.5, 0.5])
            cx = int(round(center[0] * (tf.nx - 1)))
            cy = int(round(center[1] * (tf.ny - 1)))
            cz = int(round(center[2] * (tf.nz - 1)))
            r0_vox = float(topo.get("R0_FRACTION", 0.06) * max(tf.nx, tf.ny, tf.nz))
            rhoc_vox = float(topo.get("RHOC_VOXELS", 3.0))
            biaxial_delta = float(topo.get("BIAXIAL_DELTA", 0.3))
            seeds.seed_biaxial_hedgehog_M(tf, cx, cy, cz, r0_vox, rhoc_vox, biaxial_delta)
            # NO auto-relax: the biaxial M is constructed directly; relax_director_step would
            # rebuild M uniaxially from the director and destroy the biaxial structure.
            print(
                f"[M5.6.5a] seeded biaxial hedgehog D=diag(1,{biaxial_delta},0) at "
                f"({cx},{cy},{cz}); r0={r0_vox:.1f}, ρc={rhoc_vox:.1f} voxels (no relax)"
            )

        elif seed_mode == "charged_ring":
            # M5.21.2: the CHARGED DISCLINATION RING — the hedgehog's same-charge
            # partner (Alexander et al., RMP 84, 497 (2012) § IV.B: loops carry odd
            # hedgehog charge). Half-disclination cord of radius a in z = 0, escaped
            # interior, hedgehog far field. STATIC VIEWING seed (no 4D activation,
            # no relax — constructed directly like the biaxial mode).
            tf = state.tensor_field
            center = topo.get("CENTER", [0.5, 0.5, 0.5])
            cx = int(round(center[0] * (tf.nx - 1)))
            cy = int(round(center[1] * (tf.ny - 1)))
            cz = int(round(center[2] * (tf.nz - 1)))
            a_vox = float(topo.get("A_FRACTION", 0.10) * max(tf.nx, tf.ny, tf.nz))
            cord_vox = float(topo.get("CORD_VOXELS", 3.0))
            rhoc_vox = float(topo.get("RHOC_VOXELS", 3.0))
            biaxial_delta = float(topo.get("BIAXIAL_DELTA", 0.3))
            seeds.seed_charged_ring_M(tf, cx, cy, cz, a_vox, cord_vox, rhoc_vox, biaxial_delta)
            print(
                f"[M5.21.2] seeded charged ring a={a_vox:.1f} voxels at ({cx},{cy},{cz}); "
                f"cord melt={cord_vox:.1f}, axis ρc={rhoc_vox:.1f}, δ={biaxial_delta} "
                f"(static view, no relax)"
            )

        elif seed_mode == "dressed_hedgehog":
            # M5.8.2c — the BOOST-DRESSED biaxial hedgehog (the 2b-1 ground state):
            # biaxial frame + melt, dressed with exp(b*·w(r)·B₁) mixing e_Θ with the
            # time axis. Activates the 4D Minkowski evolve path (signed flux on the
            # stable region; the time-freeze clamp replaced by the soft drift guard).
            tf = state.tensor_field
            center = topo.get("CENTER", [0.5, 0.5, 0.5])
            cx = int(round(center[0] * (tf.nx - 1)))
            cy = int(round(center[1] * (tf.ny - 1)))
            cz = int(round(center[2] * (tf.nz - 1)))
            r0_vox = float(topo.get("R0_FRACTION", 0.06) * max(tf.nx, tf.ny, tf.nz))
            rhoc_vox = float(topo.get("RHOC_VOXELS", 3.0))
            biaxial_delta = float(topo.get("BIAXIAL_DELTA", 0.3))
            b_star = float(topo.get("B_STAR", 0.13))  # 2b-1 GEM-dip dressing
            rw_frac = float(topo.get("RW_FRACTION", 0.29))  # ~3.5/12 of the box (2b-1)
            rw_vox = rw_frac * max(tf.nx, tf.ny, tf.nz)
            kick = float(topo.get("CLOCK_KICK", 0.0))  # clock phase kick (rad)
            seeds.seed_dressed_hedgehog_M(
                tf, cx, cy, cz, r0_vox, rhoc_vox, biaxial_delta, b_star, rw_vox, kick
            )
            pde.compute_stable_mask(tf)  # per-voxel ghost guard (once)
            # smooth the mask (seed-time, numpy): the hard signed/Euclid seam at the
            # core pumped the field (2c-2 diagnosis) — a ~3-voxel transition kills
            # the constitutive shock; the flux kernel blends fractional masks.
            mnp = tf.stable_mask.to_numpy()
            import numpy as _np

            sm = mnp.astype(_np.float32)
            for _ in range(3):
                for ax in range(3):
                    p = _np.swapaxes(sm, 0, ax)
                    q = p.copy()
                    q[1:-1] = 0.25 * (p[2:] + p[:-2] + 2.0 * p[1:-1])
                    sm = _np.swapaxes(q, 0, ax)
            tf.stable_mask.from_numpy(sm)
            # SIGNED_FLUX_4D (default OFF — the safe GUI v1): the Minkowski-signed
            # flux with cheap (diagonal) inertia has slow growing modes that every
            # heuristic only delays (2c-2: explodes @120/250/700 as fixes stack) —
            # the stable signed evolution needs the FULL spectral-projection
            # constrained kernel (validated f64 in 2c-1; Taichi port = follow-up).
            # OFF ⇒ mask=0 everywhere ⇒ Euclidean flux (bounded, 2c-2 variant c)
            # with the time-axis components still LIVE (dressing + drift guard +
            # dressed-V + inertia all active).
            if not bool(topo.get("SIGNED_FLUX_4D", False)):
                tf.stable_mask.fill(0.0)
            pde.compute_tstar(tf)  # V pins to the DRESSED seed
            state.evolve_4d = True
            state.km_inertia_4d = float(topo.get("KM_INERTIA_4D", 30.0))
            # dt guard (2c seed e): the DT_SCALE_4D derate is applied INSIDE
            # compute_timestep (which the main loop re-runs every frame — a
            # one-shot derate here gets wiped, the 2026-06-05 signed-GUI bug).
            dt_scale = float(topo.get("DT_SCALE_4D", 0.5))
            # OPTION B (M5.8.2c): INTEGRATOR_4D "constrained" runs the Minkowski-
            # SIGNED dynamics under the spectral-projection kernel (B-1-validated:
            # sandbox_v8/m5_8_2cb). The leapfrog default keeps the safe v1.
            state.integrator_4d = str(topo.get("INTEGRATOR_4D", "leapfrog"))
            if state.integrator_4d == "constrained":
                # act mask — Dirichlet at the seed (interior ∧ off-core ∧
                # off-axis, the validated 2c-1 boundary treatment; the frozen
                # core/axis is a documented caveat for the G-2c-1 design)
                idx = _np.indices((tf.nx, tf.ny, tf.nz)).astype(_np.float32)
                rxv, ryv, rzv = idx[0] - cx, idx[1] - cy, idx[2] - cz
                r_vox = _np.sqrt(rxv**2 + ryv**2 + rzv**2)
                rho_vox = _np.sqrt(rxv**2 + ryv**2)
                act = _np.zeros((tf.nx, tf.ny, tf.nz), _np.float32)
                act[2:-2, 2:-2, 2:-2] = 1.0
                act *= (r_vox > 2.0 * rhoc_vox) * (rho_vox > rhoc_vox)
                tf.act4d.from_numpy(act)
                state.n_act_pl = float(
                    act[tf.nx // 2].sum() + act[:, tf.ny // 2].sum() + act[:, :, tf.nz // 2].sum()
                )
                # the 10 orthonormal symmetric 4×4 basis matrices (Frobenius)
                basis = _np.zeros((10, 4, 4), _np.float32)
                for a_ in range(4):
                    basis[a_, a_, a_] = 1.0
                bi = 4
                for a_ in range(4):
                    for c_ in range(a_ + 1, 4):
                        basis[bi, a_, c_] = basis[bi, c_, a_] = 1.0 / _np.sqrt(2.0)
                        bi += 1
                tf.sym_basis.from_numpy(basis)
                # V at the 4×/τ-units convention: K·0.5/dx⁴ (c² cancels)
                state.ldg_cc_4d = 0.5 * state.ldg_c / (state.c_amrs * state.c_amrs)
                # P₀ = A(kick·M_ψ) — the 2c-1 velocity-kick semantics
                # (dt-independent; kick=0 ⇒ P₀ = 0 exactly)
                pde.init_P_4d(tf, kick)
                # Display hygiene: the seeder encodes the kick in M_prev at the
                # legacy 1/dt semantics — the energyH view's ½‖(M−M_prev)/dt‖²
                # boot-kinetic is then ~(1/DT_SCALE)² inflated, saturating the
                # colormap (the "energyH blank/uniform" report, 2026-06-05).
                # The constrained dynamics never reads M_prev (the kick lives in
                # P via init_P_4d) — sync it so the displayed velocity starts at
                # 0 and tracks M − M_prev = dt_eff·Ṁ after each swap.
                tf.M_prev_am.copy_from(tf.M_am)
                state.guard4d_frame = 0
                dt_eff = state.c_amrs * state.dt_rs
                print(
                    f"[M5.8.2c-B] CONSTRAINED signed integrator ON: act region "
                    f"{100.0 * act.mean():.1f}% of grid, dt_eff={dt_eff:.4f} am, "
                    f"V cc_4d={state.ldg_cc_4d:.3e}"
                )
            mask_np = tf.stable_mask.to_numpy()
            print(
                f"[M5.8.2c] seeded DRESSED hedgehog b*={b_star}, r_w={rw_vox:.1f} vox, "
                f"kick={kick}; stable region {100.0 * mask_np.mean():.1f}% of grid; "
                f"4D evolve ON (dt×{dt_scale}, integrator={state.integrator_4d})"
            )

        elif seed_mode == "npz_file":
            # M5.23.2 arm (3) — the research-ENDPOINT loader: npz → launcher
            # seed. Loads a converged research state (a fixed-J isorotation
            # endpoint, the μ/τ census states, any saved M array) as a LIVE
            # renderable field instead of an analytic ansatz.
            #
            # CONVENTIONS (measured at the loader build, m5_23_2 task_details):
            #   - research npz storage keeps M[0,0] = +g (the plain storage
            #     convention); the production covariant vacuum is diag(−g, 1,
            #     δ, 0) — V4 measured 759.4/cell at +g vs 0.0 at −g, and the
            #     fixed-J endpoint far field flipped lands on the production
            #     vacuum to 2e-11. So 4×4 data is flipped to −g AT LOAD
            #     (state.npz_covariant tells the canonical activation to skip
            #     its own flip); 3×3 data takes the launcher's normal path
            #     (embed4 at +g here, flip_time_axis at canonical activation).
            #   - the TensorField grid is ODD; research grids are even (32³ /
            #     64³). GRID FIT: launcher grid == n_src − 1 → CROP (drop the
            #     last plane per axis, sponge-zone content); launcher grid ≥
            #     n_src → EMBED centered with edge-extend padding (the source
            #     boundary IS the far field). Anything else → vacuum fallback
            #     with a WARNING. The data array itself is never resampled.
            #   - a loaded state is an INITIAL CONDITION in this arena, not a
            #     launcher equilibrium: drift under free evolution (pinned-bc
            #     sources especially) is physics, not a loader bug.
            tf = state.tensor_field
            npz_path = str(topo.get("PATH", ""))
            if not os.path.isabs(npz_path):
                npz_path = os.path.join(os.path.dirname(__file__), npz_path)
            if not os.path.isfile(npz_path):
                print(f"[M5.23.2] WARNING: npz not found: {npz_path}; seeding vacuum")
                seeds.seed_vacuum_M(tf, tf.lc_delta)
            else:
                info_npz = seeds.load_npz_M(
                    tf, npz_path, delta_cfg=topo.get("BIAXIAL_DELTA", tf.lc_delta)
                )
                for w_npz in info_npz["warn"]:
                    print(f"[M5.23.2] WARNING: {w_npz}")
                if not info_npz["ok"]:
                    print("[M5.23.2] load failed; seeding vacuum")
                    seeds.seed_vacuum_M(tf, tf.lc_delta)
                else:
                    state.npz_covariant = bool(info_npz["covariant"])
                    state.npz_h = float(info_npz["h"])
                    print(
                        f"[M5.23.2] loaded {os.path.basename(npz_path)}: "
                        f"{info_npz['n_src']}^3 x {info_npz['blk']}x{info_npz['blk']} "
                        f"-> {tf.nx}x{tf.ny}x{tf.nz} ({info_npz['fit']}); "
                        f"h = {info_npz['h'] or '?'}, delta = {info_npz['delta']}, "
                        f"bc = {info_npz['bc']}"
                        + (
                            ", flipped to covariant -g at load"
                            if info_npz["covariant"]
                            else ""
                        )
                    )

        else:
            print(f"[M5.1] WARNING: unknown TOPOLOGY_SEED mode: {seed_mode!r}")

        # Optional auto-relaxation (M5.1 task 6)
        if state.AUTO_RELAX_STEPS > 0 and state.n_defects > 0:
            print(f"[M5.1] auto-relaxing for {state.AUTO_RELAX_STEPS} gradient-descent steps...")
            relax_field(state, state.AUTO_RELAX_STEPS)
            print("[M5.1] auto-relax complete")

        # M5.24 — the CANONICAL regularization-stack integrator (verified-L era
        # port): activated by TOPOLOGY_SEED INTEGRATOR_4D = "canonical". Flips
        # the seeded time axis to the covariant vacuum (M[0,0] = −g) and routes
        # compute_propagation through the η-stack leapfrog + V4. Placed AFTER
        # auto-relax (rebuild_M_from_director writes +g back — the flip must be
        # the last touch on M). Intended for the BIAXIAL-family seeds, whose far
        # field diag(1, δ, 0) IS the covariant vacuum spatial spectrum; on the
        # uniaxial seeds (far field (1, δ, δ)) V4 charges the whole background.
        if str(topo.get("INTEGRATOR_4D", "")) == "canonical":
            if seed_mode == "dressed_hedgehog":
                print(
                    "[M5.24] INTEGRATOR_4D 'canonical' is unsupported on the dressed "
                    "seed (boost-mixed time row is an era-specific object); keeping "
                    "the era 4D path."
                )
            else:
                if seed_mode not in ("biaxial_hedgehog", "charged_ring", "vacuum", "npz_file"):
                    print(
                        f"[M5.24] WARNING: canonical stack on seed mode {seed_mode!r}: "
                        "the far field is not the covariant vacuum, V4 charges the "
                        "background (biaxial-family seeds are the intended use)."
                    )
                # M5.23.2 — a 4×4 npz load is ALREADY covariant (flipped at
                # load, measured V4 = 0 far field); flipping again would put
                # it back on the +g ridge (V4 = 759/cell).
                if not getattr(state, "npz_covariant", False):
                    pde.flip_time_axis(state.tensor_field)
                state.evolve_4d = True
                state.integrator_4d = "canonical"
                state.eta_delta = float(topo.get("BIAXIAL_DELTA", state.tensor_field.lc_delta))
                # ETA_DX (round 2): the grid spacing in the CANONICAL unit
                # system. Set it to the research grid unit (1.5 = the
                # m5_21_2b/3 h) to run the launcher as a research-twin; the
                # physical dx_am fallback weakens the curvature ~100× against
                # the dimensionless V4 (the unit-map finding, task_details).
                # M5.23.2 — an npz load's recorded grid unit h is the natural
                # ETA_DX for that state (research-twin geometry); an explicit
                # ETA_DX in the config still wins.
                eta_dx_default = getattr(state, "npz_h", 0.0) or state.tensor_field.dx_am
                state.eta_dx = float(topo.get("ETA_DX", eta_dx_default))
                # ETA_SUBSTEPS (round 2, viz speed): N physics steps per
                # rendered frame at the CERTIFIED dt (dt itself cannot rise:
                # the stiff-mode wall is 2/78 ≈ 0.026 τ, NaN measured at
                # 0.02). N× visual speed, zero physics change.
                state.eta_substeps = max(1, int(topo.get("ETA_SUBSTEPS", 8)))
                # ETA_SPONGE_* (round 3): boundary sponge for long live runs
                # (kills the pinned-box reflection re-agitation, the t ≈ 12
                # caveat). GAMMA 0 = off (exact no-op).
                state.eta_sponge_gamma = float(topo.get("ETA_SPONGE_GAMMA", 0.0))
                state.eta_sponge_width = float(topo.get("ETA_SPONGE_WIDTH", 10.0))
                # CANON_RELAX_ITERS (round 3): optional seed-time FIRE relax
                # (the RELAX button runs chunks interactively either way).
                canon_relax = int(topo.get("CANON_RELAX_ITERS", 0))
                print(
                    f"[M5.24] CANONICAL stack ON: eta-curvature + V4 "
                    f"(w = {pde.W1_SPECTRAL:.3e}, delta = {state.eta_delta}), "
                    f"covariant vacuum M[0,0] = -g, dt_eff capped at "
                    f"{float(topo.get('DT_ETA_CAP', 0.005))} tau, "
                    f"dx_eta = {state.eta_dx:.3f}, substeps/frame = "
                    f"{state.eta_substeps}, sponge gamma = "
                    f"{state.eta_sponge_gamma} (width {state.eta_sponge_width})"
                )
                if canon_relax > 0:
                    relax_field_canonical(state, canon_relax)
                # M5.23.1 — the FIXED-J isorotation clock (RELAX → SET-J →
                # EVOLVE): FIXEDJ_OMEGA ≠ 0 sets the seed-time isorotation
                # velocity Ṁ(0) = ω*·a0 on the conjugation-tangent clock
                # flow, AFTER any seed-time relax (FIRE clobbers the Md_am
                # a0 scratch and re-zeros Ṁ, so SET-J is the last touch on
                # the buffers). The SET J button redoes it interactively.
                state.fixedj_renv = float(topo.get("FIXEDJ_RENV", 10.0))
                state.fixedj_log_every = int(topo.get("FIXEDJ_LOG_EVERY", 0))
                state.fixedj_omega = float(topo.get("FIXEDJ_OMEGA", 0.0))
                if state.fixedj_omega != 0.0:
                    set_fixed_j_launcher(state, state.fixedj_omega)

        # VIZ.3: populate the derived eigenframe (director_nhat + director_mid +
        # eigenvalues) from the seeded M so a PAUSED boot renders the δ-clock-hand
        # glyph correctly before the first Evolve-PDE step. Cheap, runs once.
        pde.eigen_decompose(state.tensor_field)

    # VIZ.5 (M5.23) — ellipsoid-shell center(s), voxel coords: the seeded defect
    # centers when a multi-defect DEFECTS list was pinned this seed, else the
    # config CENTER (universe center for vacuum / no-seed configs). Consumed
    # per-frame by engine4_render.update_ellipsoid_mesh. Decided off the
    # TOPOLOGY_SEED dict (not pin_centers presence) so a config switch never
    # inherits stale centers from the previous seed.
    tf_ell = state.tensor_field
    topo_ell = state.TOPOLOGY_SEED or {}
    if "DEFECTS" in topo_ell and state.pin_centers is not None and state.n_defects > 0:
        n_ell = min(state.n_defects, tf_ell.ellipsoid_max_centers)
        for c_ell in range(n_ell):
            tf_ell.ellipsoid_centers[c_ell] = [
                float(state.pin_centers[c_ell, 0]),
                float(state.pin_centers[c_ell, 1]),
                float(state.pin_centers[c_ell, 2]),
            ]
    else:
        n_ell = 1
        c_norm = topo_ell.get("CENTER", [0.5, 0.5, 0.5])
        tf_ell.ellipsoid_centers[0] = [
            c_norm[0] * (tf_ell.nx - 1),
            c_norm[1] * (tf_ell.ny - 1),
            c_norm[2] * (tf_ell.nz - 1),
        ]
    tf_ell.ellipsoid_n_centers[None] = n_ell

    if state.INSTRUMENTATION:
        print("\n" + "=" * 64)
        print("INSTRUMENTATION ENABLED")
        print("=" * 64)


def relax_field(state, n_steps):
    """Run N gradient-descent relaxation steps on the director field (M5.1 task 6).

    Lowers the Frank elastic energy by smoothing the seeded hedgehog blend
    zone, while preserving topology via soft-core pinning. relax_director_step
    writes director_nhat_new; rebuild_M_from_director then writes all three M
    buffers equal (Ṁ = 0) so subsequent Evolve PDE starts from a static state.

    Args:
        state: SimulationState with tensor_field + pin_centers/signs/n_defects
        n_steps: number of relax_director_step iterations to run

    Step size: τ = 0.4 · dx²/6 (40% of the heat-equation CFL bound for headroom).
    Empirically Exp 2 used τ = 0.008 with dx = 0.34 → 41% of CFL; we match.
    """
    if state.pin_centers is None or state.n_defects == 0:
        # No defects → nothing to pin; gradient descent on a uniform field
        # has nothing to do (∇²n = 0 everywhere). Skip.
        return
    tf = state.tensor_field
    cfl_bound = (tf.dx_am**2) / 6.0  # τ < dx²/(2·dim·K) with dim=3, K=1
    tau = 0.4 * cfl_bound
    for _ in range(n_steps):
        pde.relax_director_step(tf, tau, state.pin_centers, state.pin_signs, state.n_defects)
        # M5.4: relaxation operates on director_nhat (the principal eigenvector of M).
        # Copy the relaxed director back. NOT a buffer swap — static update.
        tf.director_nhat.copy_from(tf.director_nhat_new)
    # Rebuild M from the relaxed director so M_am, the ‖M−D‖/‖Ṁ‖ trackers, and the
    # flux-mesh orientation mode stay consistent with the relaxed director.
    pde.rebuild_M_from_director(tf, tf.lc_delta)
    # Refresh observables so dashboard + flux mesh reflect the relaxed field
    compute_field_observables(state)


def relax_field_canonical(state, n_iters):
    """M5.24 round 3 (G8) — FIRE-relax the current state under the canonical
    η + V4 energy (frozen-time-row sector, pinned boundary): the production
    port of the canonical § 5.2 statics recipe, single-sourced in
    engine2_pde.fire_relax_canonical (the selftest gates the same function).
    Refreshes the derived fields + observables so the render shows the
    relaxed state immediately."""
    tf = state.tensor_field
    stats = pde.fire_relax_canonical(
        tf,
        getattr(state, "eta_dx", tf.dx_am),
        tf.lc_g,
        getattr(state, "eta_delta", tf.lc_delta),
        pde.W1_SPECTRAL,
        int(n_iters),
    )
    pde.eigen_decompose(tf)
    compute_field_observables(state)
    print(
        f"[M5.24 RELAX] {stats['iters']} FIRE iters: force RMS "
        f"{stats['f_rms0']:.3e} -> {stats['f_rms1']:.3e} "
        f"(dt_final {stats['dt_final']:.4f}); press again to continue "
        f"relaxing, or EVOLVE to run dynamics from here."
    )


def set_fixed_j_launcher(state, om_target):
    """M5.23.1 — the SET-J step of the RELAX → SET-J → EVOLVE flow: build the
    conjugation-tangent clock flow a0, measure the clock inertia kin, and
    kick the isorotation velocity Ṁ(0) = ω*·a0 into the triple buffer
    (engine2_pde.set_fixed_j, gated by the m5_23_1 selftest). Canonical path
    only. CONVENTION: kin/J are flow-parameter values on the conjugation
    tangent (kin = 0.1206 on the certified research state); any ABSOLUTE
    J / ħ/2 / g statement needs the physical-rate convention
    (research m5_21_5_note § 5) — never quote these numbers as physical
    angular momentum."""
    tf = state.tensor_field
    if getattr(state, "integrator_4d", "") != "canonical":
        print("[M5.23.1] SET-J needs the canonical integrator path — skipped.")
        return None
    dt_eff = state.c_amrs * state.dt_rs
    info = pde.set_fixed_j(
        tf,
        getattr(state, "eta_dx", tf.dx_am),
        float(getattr(state, "fixedj_renv", 10.0)),
        float(om_target),
        dt_eff,
        shell=1,
    )
    if info["set"]:
        # the CARRIED charge J_self = <Mdot, a0(M)> is the hold observable
        # (== omega* here by construction); the global [G_k, M] projections
        # near-cancel on hedgehog-family states (selftest S4 measurement)
        js = pde.read_carried_j(
            tf, getattr(state, "eta_dx", tf.dx_am), float(getattr(state, "fixedj_renv", 10.0)), dt_eff
        )
        print(
            f"[M5.23.1] SET-J: kin = {info['kin']:.5e}, J = {info['J']:.5e}, "
            f"omega* = {info['om_star']:.4f} (flow-parameter convention, "
            f"conjugation tangent; dt_eff = {dt_eff:.4f}); carried "
            f"J_self = {js:+.4f}. EVOLVE now runs the live clock — keep dt "
            f"unchanged (the kick encodes velocity at this dt_eff)."
        )
    else:
        print(
            f"[M5.23.1] SET-J DECLINED: kin = {info['kin']:.3e} <= 0 "
            f"(nonpositive clock inertia — the state left the positive-kin "
            f"sector; no honest isorotation to set)."
        )
    return info


def compute_propagation(state):
    """Per-step field evolution — M5.5.4: the Eq.18 matrix-action leapfrog ("Evolve PDE").

    Evolves the Landau–de Gennes order parameter M under Duda's Eq.18 action
    (simple ½‖Ṁ‖² kinetic + faithful potential, validated symplectic in
    sandbox_v5/m5_5_4_matrix_evolution_check.py):

        ∂²_t M = c²·Σ_α ∂_α G_α − dV_M(M) ,   G_α = 8 Σ_ν [[M_α,M_ν],M_ν]

    Two-pass: compute_curvature_flux (G_α) → evolve_M (leapfrog) → swap_matrix_buffers,
    then eigen_decompose refreshes the derived director_nhat from the evolved M so the
    render + the ‖M−D‖_F/‖Ṁ‖_F trackers reflect the new field. V off (ldg_a/b/c=0) →
    free curvature dynamics; the faithful curvature kinetic (degenerate metric, O-DoF)
    is the M5.6 refinement. (The retired ψ leapfrog stays dormant; see engine2_pde.)
    """
    tf = state.tensor_field
    if getattr(state, "integrator_4d", "") == "canonical":
        # M5.24 — the canonical regularization stack (the § 2 row 6 formulation
        # of record): symmetrized two-branch η-flux + exact adjoint gather + V4
        # force, τ-units (dt_eff = c·dt), dt capped by compute_timestep. The
        # RUNNABLE regularization — never present its dynamics as the true-L
        # evolution (the free-EL IVP is ill-posed, canonical § 2 row 1).
        dt_eff = state.c_amrs * state.dt_rs
        dx_eta = getattr(state, "eta_dx", tf.dx_am)
        eta_delta = getattr(state, "eta_delta", tf.lc_delta)
        # ETA_SUBSTEPS physics steps per rendered frame (viz speed at the
        # certified dt). The last substep leaves its swap to the shared
        # swap_matrix_buffers() below (all paths swap exactly once there).
        n_sub = max(1, int(getattr(state, "eta_substeps", 1)))
        sp_gamma = float(getattr(state, "eta_sponge_gamma", 0.0))
        sp_width = float(getattr(state, "eta_sponge_width", 10.0))
        for sub_ in range(n_sub):
            pde.compute_eta_flux(tf, 0, dx_eta)
            pde.evolve_M_eta_start(tf, dt_eff, dx_eta)
            pde.compute_eta_flux(tf, 1, dx_eta)
            pde.evolve_M_eta_finish(tf, dt_eff, dx_eta, tf.lc_g, eta_delta, pde.W1_SPECTRAL)
            if sp_gamma > 0.0:
                # round-3 boundary sponge: absorbs outgoing radiation instead
                # of reflecting it back onto the core (the t ≈ 12 caveat)
                pde.apply_eta_sponge(tf, dt_eff, sp_gamma, sp_width)
            if sub_ < n_sub - 1:
                tf.swap_matrix_buffers()
        # Bounded-energy guard (same pattern as the constrained path): block-diag
        # seeds live on an invariant manifold with positive energy, but f32
        # roundoff slowly leaks into the time-row sector where the η-potential
        # is INDEFINITE (canonical § 1: H unbounded below) — a long run can
        # eventually grow there with E conserved. Auto-pause loudly instead of
        # scrambling the render (the 2026-06-05 signed-GUI failure UX).
        state.guard4d_frame = getattr(state, "guard4d_frame", 0) + 1
        if state.guard4d_frame % 60 == 0:
            m_max = float(abs(tf.M_am.to_numpy()).max())
            if not (m_max < 50.0):  # catches NaN too
                state.PAUSED = True
                print(
                    f"[M5.24 GUARD] max|M| = {m_max:.3e} at step "
                    f"{state.guard4d_frame} — bounded-energy guard tripped, "
                    f"AUTO-PAUSED (the indefinite time-row channel grew; "
                    f"vacuum scale ~8, threshold 50)."
                )
        # M5.23.1 — optional periodic J readout (FIXEDJ_LOG_EVERY frames;
        # 0 = off): the CARRIED isorotation charge J_self = <Mdot, a0(M)>
        # (== omega* at SET-J; its drift = the honest hold/decay read; the
        # global-flow projections near-cancel on hedgehog states, selftest S4)
        log_n = int(getattr(state, "fixedj_log_every", 0))
        if log_n > 0 and state.guard4d_frame % log_n == 0:
            js = pde.read_carried_j(
                tf, dx_eta, float(getattr(state, "fixedj_renv", 10.0)), dt_eff
            )
            print(f"[M5.23.1] frame {state.guard4d_frame}: carried J_self = {js:+.4f}")
    elif (
        getattr(state, "evolve_4d", False) and getattr(state, "integrator_4d", "") == "constrained"
    ):
        # M5.8.2c OPTION B — the Minkowski-SIGNED dynamics under the constrained
        # spectral-projection kernel (B-1-validated, sandbox_v8/m5_8_2cb). The
        # exact 2c-1 step order: flux(M, Ṁ_prev) → P += dt·force → global (α,3)
        # clamp → per-voxel solve (A-eig, keep, P-project, Ṁ = A⁺P) → M update.
        dt_eff = state.c_amrs * state.dt_rs
        pde.flux_4d_constrained(tf)
        pde.update_P_4d(tf, dt_eff, getattr(state, "ldg_cc_4d", 0.0))
        pde.sample_p03_drift(tf)
        n_pl = max(getattr(state, "n_act_pl", 0.0), 1.0)
        pde.apply_p03_clamp(
            tf, tf.v03_sums[0] / n_pl, tf.v03_sums[1] / n_pl, tf.v03_sums[2] / n_pl
        )
        pde.solve_constrained_4d(tf)
        pde.update_M_4d_constrained(tf, dt_eff)
        # GUARD (b) — the bounded-energy monitor (load-bearing, roadmap 2c):
        # every 60 steps pull M and auto-PAUSE on divergence — a NaN/blow-up
        # must stop the run with a loud message, not scramble the render at
        # 3 FPS (the 2026-06-05 signed-GUI failure UX).
        state.guard4d_frame = getattr(state, "guard4d_frame", 0) + 1
        if state.guard4d_frame % 60 == 0:
            m_chk = tf.M_am.to_numpy()
            m_max = float(abs(m_chk).max())
            if not (m_max < 50.0):  # catches NaN too (NaN comparisons are False)
                state.PAUSED = True
                print(
                    f"[M5.8.2c-B GUARD] max|M| = {m_max:.3e} at step "
                    f"{state.guard4d_frame} — bounded-energy guard tripped, "
                    f"AUTO-PAUSED. The signed run diverged; report this step "
                    f"count (vacuum scale is ~8; threshold 50)."
                )
    elif getattr(state, "evolve_4d", False):
        # M5.8.2c — the 4D Minkowski path (dressed-hedgehog seed): signed flux on
        # the stable region + live time axis + the coherent-(α,3)-drift guard.
        pde.compute_curvature_flux_4d(tf)
        pde.sample_v03_drift(tf, state.dt_rs)
        n_pl = float(tf.nx * tf.ny + tf.ny * tf.nz + tf.nx * tf.nz)
        vm = [tf.v03_sums[a] / n_pl for a in range(3)]
        pde.evolve_M_4d(
            tf,
            state.c_amrs,
            state.dt_rs,
            state.ldg_c,
            vm[0],
            vm[1],
            vm[2],
            getattr(state, "km_inertia_4d", 30.0),
        )
    else:
        pde.compute_curvature_flux(tf)
        pde.evolve_M(tf, state.c_amrs, state.dt_rs, state.ldg_a, state.ldg_b, state.ldg_c)
    tf.swap_matrix_buffers()
    pde.eigen_decompose(tf)  # refresh director_nhat from the evolved M (for render + trackers)


def compute_field_observables(state):
    """Compute per-voxel energy densities + global aggregates from the current ψ.

    Runs every frame regardless of pause state so the flux-mesh visualization
    and dashboard reflect the current field — including the static seeded
    state when PAUSED=True. The kernels here read ψ but don't modify it, so
    they are safe to run while dynamics are halted.

    When paused at the seeded state, ψ == ψ_prev, so the kinetic term in H
    is zero and H reduces to ½c²|∇ψ|² + V(ψ) — structurally similar to F up
    to a c² factor (with V=0 in M5.1).
    """
    # PER-VOXEL MATRIX HAMILTONIAN (M5.5.4) =============================
    # H = ½‖Ṁ‖² + c²·4Σ‖[M_μ,M_ν]‖² + V_M(M)  → observables.energyH_density_aJ.
    # Replaces the dormant-ψ placeholder (uniform ¼λ) — resolves the M5.4 WAVE_MENU=4
    # carry-over. e_scale=1.0 (bare units; the physical-energy calibration is tied to a
    # reference mass scale → deferred to M5.9). Consumed by flux-mesh WAVE_MENU=4.
    if getattr(state, "integrator_4d", "") == "canonical":
        # M5.24 — the verified-L era H (canonical kinetic + η-curvature + V4);
        # no v0 subtraction: the covariant vacuum is the exact zero.
        observables.compute_energyH_density_eta(
            state.tensor_field,
            state.observables,
            state.c_amrs * state.dt_rs,
            getattr(state, "eta_dx", state.tensor_field.dx_am),
            state.tensor_field.lc_g,
            getattr(state, "eta_delta", state.tensor_field.lc_delta),
            pde.W1_SPECTRAL,
            1.0,
        )
    else:
        observables.compute_energyH_density_M(
            state.tensor_field,
            state.observables,
            state.c_amrs,
            state.dt_rs,
            state.ldg_a,
            state.ldg_b,
            state.ldg_c,
            state.ldg_v0,
            1.0,
        )

    # PER-VOXEL FRANK ELASTIC ENERGY DENSITY (M5.1 task 5) ===============
    # F = (K/2)·|∇n̂|²  → observables.energyF_density_aJ.
    # Consumed by flux-mesh WAVE_MENU=5; M5.1 task 6 (gradient-descent
    # monotone-decrease diagnostic); M5.1 task 7 (Coulomb 1/d fit).
    observables.compute_energyF_density(state.tensor_field, state.observables, observables.K_FRANK)

    # EM-FROM-TILTS OBSERVABLES (M5.6.5b "see EM") ======================
    # ∇·n̂ (splay = Coulomb-charge-like) + ‖∇×n̂‖ (twist+bend = B-like circulation).
    # Consumed by flux-mesh WAVE_MENU 6 (∇·n̂, bluered) / 7 (‖∇×n̂‖, ironbow).
    observables.compute_director_em(state.tensor_field, state.observables)
    # VIZ.4: overwrite the curl (B) field with an analytic dipole placeholder so the
    # N/S coloring + B glyphs render correctly before the real circulating B exists
    # (M5.8). div (charge) is left as the real seeded splay for context. Same center-
    # voxel convention as the biaxial_hedgehog seed so the dipole sits on the defect.
    if state.DIPOLE_SAMPLE:
        tf = state.tensor_field
        cx = state.DIPOLE_CENTER[0] * (tf.nx - 1)
        cy = state.DIPOLE_CENTER[1] * (tf.ny - 1)
        cz = state.DIPOLE_CENTER[2] * (tf.nz - 1)
        observables.fill_dipole_sample_B(
            tf, state.observables, state.DIPOLE_AXIS, cx, cy, cz, state.DIPOLE_R0_VOX, 1.0
        )
    observables.compute_director_em_scale(state.tensor_field, state.observables)

    # MATRIX-SUBSTRATE TRACKERS (M5.4) ==================================
    # ‖M−D‖_F amplitude (A) + ‖Ṁ‖_F clock-ω. EMA on the current M; runs
    # every frame (incl. paused) so the flux-mesh "Amp (A)"/"Clock (omega)"
    # modes + dashboard reflect the static seeded/relaxed state. Replaces the
    # ψ-based update_trackers (which now only fires on the retiring EVOLVE PSI path).
    observables.update_trackers_M(
        state.tensor_field, state.trackers, state.dt_rs, state.tensor_field.lc_delta
    )

    # IN-FRAME DATA SAMPLING & ANALYTICS ================================
    # Frame skip reduces GPU->CPU transfer overhead during dynamics; when
    # paused we run every frame because the field doesn't change and the
    # cost is dominated by cached slice copies, not the reduction.
    # Two separate samplers per the 2026-05-11 SoC refactor: Trackers and
    # FieldObservables each own their own 3-plane pass.
    if state.frame % 60 == 0 or state.frame == 10 or state.PAUSED:
        observables.sample_avg_trackers(state.tensor_field, state.trackers)
        observables.sample_avg_observables(state.tensor_field, state.observables)
        state.energyH_total = (
            state.observables.energyH_global_avg_aJ[None] * state.tensor_field.voxel_count
        ) * constants.ATTOJOULE  # J
    state.amp_global_rms = state.trackers.amp_global_emarms_am[None] * constants.ATTOMETER  # m
    state.freq_global_avg = state.trackers.freq_global_avg_rHz[None] / constants.RONTOSECOND  # Hz
    state.wavelength_global_avg = constants.WAVE_SPEED / (state.freq_global_avg or 1)  # no 0 div
    state.energyH_global_avg = state.observables.energyH_global_avg_aJ[None] * constants.ATTOJOULE
    state.energyF_global_avg = state.observables.energyF_global_avg_aJ[None] * constants.ATTOJOULE

    if state.INSTRUMENTATION:
        instrument.log_timestep_data(state.frame, state.tensor_field, state.trackers)
        if state.frame == 500:
            instrument.plot_probe_wave_profile(state.tensor_field)


def _curl_projection(state):
    """`(curl_radial, curl_center)` for the bluered N/S projection of ∇×n̂ (B).

    Shared by the WM7 flux mesh AND the B-field glyphs so their N/S coloring is
    identical. For the dipole sample → radial `(∇×n̂)·r̂` about the defect center
    (true N=red/S=blue poles, matching a bar magnet); general runs → fixed-axis
    projection (radial=0, center unused). center is in voxel coords.
    """
    if state.DIPOLE_SAMPLE:
        tf = state.tensor_field
        return 1, ti.Vector(
            [
                state.DIPOLE_CENTER[0] * (tf.nx - 1),
                state.DIPOLE_CENTER[1] * (tf.ny - 1),
                state.DIPOLE_CENTER[2] * (tf.nz - 1),
            ],
            dt=ti.f32,
        )
    return 0, ti.Vector([0.0, 0.0, 0.0], dt=ti.f32)


def render_elements(state):
    """Render grid, edges, flux mesh and test particles."""
    if state.SHOW_GRID:
        render.scene.lines(state.tensor_field.grid_lines, width=1, color=colormap.COLOR_MEDIUM[1])

    if state.SHOW_EDGES:
        render.scene.lines(state.tensor_field.edge_lines, width=1, color=colormap.COLOR_MEDIUM[1])

    if state.SHOW_FLUX_MESH > 0 and state.SHOW_GRANULES == False:
        # VIZ.4: for the dipole sample, color the bluered B by the RADIAL projection
        # (∇×n̂)·r̂ about the defect center → true N/S poles (red N above / blue S below,
        # matching a bar magnet). General runs use the fixed-axis projection (radial=0).
        curl_radial, curl_center = _curl_projection(state)
        viz.update_flux_mesh_values(
            state.tensor_field,
            state.trackers,
            state.observables,
            state.WAVE_MENU,
            state.WARP_MESH,
            state.CURL_COLOR,
            state.CURL_AXIS,  # VIZ.2: fixed-axis projection for bluered N/S (default ẑ)
            curl_radial,  # VIZ.4: 1 = radial (∇×n̂)·r̂ → true N/S poles (dipole sample)
            curl_center,  # defect center in voxel coords (radial mode only)
        )
        flux_mesh.render_flux_mesh(render.scene, state.tensor_field, state.SHOW_FLUX_MESH)

    # M5.23.2 arm (4) — the energy-density ISOSURFACE: marching tetrahedra
    # on the live energyH density (already computed this frame by
    # compute_field_observables). Canonical path only: V4 + curvature are
    # exactly zero on the covariant vacuum, so a level set needs no
    # background subtraction. ISO_COVER renders the author's "surface
    # covered with ellipsoids" variant into the SHELL pool (the shell view
    # is suppressed while it is live — same buffers).
    iso_cover_live = False
    if state.SHOW_ISO:
        if getattr(state, "integrator_4d", "") != "canonical":
            if not getattr(state, "_iso_warned", False):
                print(
                    "[M5.23.2] the iso-surface needs the canonical integrator "
                    "(true-zero vacuum); ignoring SHOW_ISO on this path"
                )
                state._iso_warned = True
        else:
            tf_iso = state.tensor_field
            viz.iso_density_max(tf_iso, state.observables)
            iso_lvl = state.ISO_LEVEL * float(tf_iso.iso_level_max[None])
            if iso_lvl > 0.0:
                viz.marching_tetrahedra(tf_iso, state.observables, iso_lvl)
                viz.collapse_iso_tail(tf_iso)
                n_iso_raw = int(tf_iso.iso_tri_count[None])
                if n_iso_raw > tf_iso.iso_max_tris and n_iso_raw != getattr(
                    state, "_iso_overflow_reported", 0
                ):
                    print(
                        f"[M5.23.2] WARNING: isosurface overflow: {n_iso_raw} "
                        f"triangles at this level, budget {tf_iso.iso_max_tris} "
                        "— surface truncated (raise Iso Level)"
                    )
                    state._iso_overflow_reported = n_iso_raw
                if state.ISO_COVER:
                    viz.update_iso_ellipsoids(
                        tf_iso,
                        state.ELLIPSOID_SIZE,
                        state.ELLIPSOID_COUNT * medium.ELLIPSOID_MAX_CENTERS,
                    )
                    render.scene.mesh(
                        tf_iso.ellipsoid_mesh_vertices,
                        indices=tf_iso.ellipsoid_mesh_indices,
                        per_vertex_color=tf_iso.ellipsoid_mesh_colors,
                        two_sided=True,
                    )
                    iso_cover_live = True
                else:
                    render.scene.mesh(
                        tf_iso.iso_vertices,
                        indices=tf_iso.iso_indices,
                        per_vertex_color=tf_iso.iso_colors,
                        two_sided=True,
                    )

    # VIZ.5 (M5.23) — the "ellipsoid" viz: M·u eigen-ellipsoid surfaces, a
    # FULL-3D view (not plane cross-sections). Sub-views, combinable: the S²
    # shell (one per 3D angle around each defect center), the Stage D
    # disclination RODS (line samples along the spin axis), and the ROD RINGS
    # (one per 2D angle around the cord). Its own feature; while active it
    # REPLACES all vector glyph viz (the plane glyphs below are suppressed — no
    # visual pollution). Mesh-only since 2026-07-19.
    if state.SHOW_ELLIPSOID:
        tf = state.tensor_field
        ell_extent = state.ELLIPSOID_RADIUS * tf.max_grid_size  # shell R / rod half-length
        if state.ELLIPSOID_SHELL and not iso_cover_live:  # M5.23.2: cover owns the pool
            viz.update_ellipsoid_mesh(
                tf,
                ell_extent,
                state.ELLIPSOID_SIZE,  # ellipsoid base size (the GUI Size slider)
                min(state.ELLIPSOID_COUNT, tf.ellipsoid_max_dirs),
            )
            render.scene.mesh(
                tf.ellipsoid_mesh_vertices,
                indices=tf.ellipsoid_mesh_indices,
                per_vertex_color=tf.ellipsoid_mesh_colors,
                two_sided=True,
            )
        if state.ELLIPSOID_RODS or state.ELLIPSOID_RODRINGS:
            # ring radius: a few voxels around the cord (fixed v1; the rod
            # core melt width is ~rhoc ≈ 3 voxels in the biaxial seeds)
            ring_r_vox = max(3.0, 0.04 * tf.max_grid_size)
            # ring azimuth density follows the Count slider (~Count/12 keeps
            # the default look; buffer ceiling ellipsoid_ring_az)
            ring_az = min(max(state.ELLIPSOID_COUNT // 12, 8), tf.ellipsoid_ring_az)
            viz.update_rod_ellipsoids(
                tf,
                viz.ELLIPSOID_ROD_SPAN * ell_extent,  # rods protrude beyond the shell
                ring_r_vox,
                state.ELLIPSOID_SIZE,
                ring_az,
                1 if state.ELLIPSOID_RODS else 0,
                1 if state.ELLIPSOID_RODRINGS else 0,
            )
            render.scene.mesh(
                tf.ellipsoid_rod_vertices,
                indices=tf.ellipsoid_rod_indices,
                per_vertex_color=tf.ellipsoid_rod_colors,
                two_sided=True,
            )

    # M5.1 director-glyph overlay — line segments showing n̂ orientation,
    # signed-component RGB so opposite directions are visually opposite.
    # Renders on top of (or instead of) the flux mesh; fast (~768 segments).
    # Suppressed while the VIZ.5 ellipsoid shell is active (see above).
    elif state.SHOW_GLYPHS > 0:
        glyph_length = 0.02  # ~2% of universe edge in normalized [0,1] coords
        # VIZ.3 — 4-state glyph select (GLYPH_VECTOR):
        #   0 = Director n̂ only (principal axis; apolar, agnostic to size/color)
        #   1 = Director + Delta (n̂ + CYAN delta cross-bar = ellipsoid frame; apolar)
        #   2 = E-field (n̂ + +→− barb, charge-colored; honors size/color)   ← polar
        #   3 = B-field (∇×n̂, em-vector kernel)                             ← polar
        if state.GLYPH_VECTOR == 3:  # B-field vectors (∇×n̂) — em-vector kernel
            em_scale = max(
                state.observables.director_div_absmax[None],
                state.observables.director_curl_max[None],
            )  # shared distortion scale (matches WAVE_MENU 7)
            # gradient color uses the SAME signed N/S projection as the WM7 mesh
            curl_radial, curl_center = _curl_projection(state)
            viz.update_em_vector_glyphs(
                state.tensor_field,
                state.observables,
                glyph_length,
                em_scale,
                state.SHOW_GLYPHS,
                state.GLYPH_SIZE,
                state.GLYPH_COLOR,
                state.CURL_AXIS,
                curl_radial,
                curl_center,
            )
            render.scene.lines(
                state.tensor_field.director_glyph_vertices,
                per_vertex_color=state.tensor_field.director_glyph_colors,
                width=2.0,
            )
            # half-arrow barb pass — B vectors are polar; the tip shows circulation direction
            render.scene.lines(
                state.tensor_field.director_glyph_arrow_vertices,
                per_vertex_color=state.tensor_field.director_glyph_arrow_colors,
                width=2.0,
            )
        else:  # Director-only (0) / Director+Delta (1) / E-field (2) — centered glyph kernel
            div_scale = state.observables.director_div_absmax[None]
            glyph_mode = 1 if state.GLYPH_VECTOR == 2 else 0  # 0=Director, 1=E-field
            show_delta = 1 if state.GLYPH_VECTOR == 1 else 0  # delta cross-bar only in state 1
            viz.update_director_glyphs(
                state.tensor_field,
                state.observables,
                glyph_length,
                state.SHOW_GLYPHS,
                div_scale,
                state.GLYPH_SIZE,
                state.GLYPH_COLOR,
                glyph_mode,
                show_delta,
            )
            render.scene.lines(
                state.tensor_field.director_glyph_vertices,
                per_vertex_color=state.tensor_field.director_glyph_colors,
                width=2.0,
            )
            # Second-segment pass: delta cross-bar (state 1) OR E +→− barb (mode 1).
            # The arrow buffer is always written (blanked in director-only state 0).
            render.scene.lines(
                state.tensor_field.director_glyph_arrow_vertices,
                per_vertex_color=state.tensor_field.director_glyph_arrow_colors,
                width=2.0,
            )

    # VIZ.4 — magnetic-moment vector glyph (the placeholder dipole's axis marker).
    # Rendered independently of SHOW_GLYPHS so the moment is always visible alongside
    # the N/S-colored field. YELLOW + thick so it reads as the principal axis.
    if state.DIPOLE_SAMPLE:
        tf = state.tensor_field
        viz.update_moment_glyph(
            tf,
            state.DIPOLE_AXIS,
            state.DIPOLE_CENTER[0] * (tf.nx - 1),
            state.DIPOLE_CENTER[1] * (tf.ny - 1),
            state.DIPOLE_CENTER[2] * (tf.nz - 1),
            0.18,  # normalized length — larger than voxel glyphs (~0.02)
            ti.Vector(
                [colormap.YELLOW[1][0], colormap.YELLOW[1][1], colormap.YELLOW[1][2]], dt=ti.f32
            ),
        )
        render.scene.lines(
            tf.moment_glyph_vertices,
            per_vertex_color=tf.moment_glyph_colors,
            width=2.0,
        )

    # Render granule positional displacement (only when flux mesh is active, since position
    # is sampled from full-grid displacement data)
    if state.SHOW_GRANULES and state.SHOW_FLUX_MESH > 0:
        granule_radius = 0.001  # in screen space (relative to max universe edge and scale factor)
        amp_boost = state.WARP_MESH  # Boost granule displacement for better visibility
        nx, ny = state.tensor_field.nx, state.tensor_field.ny
        # Shared VIZ_STRIDE drives both granules and director glyphs — same density,
        # so granules visually align with the underlying glyphs (M5.1 consolidation).
        stride = max(1, state.VIZ_STRIDE)
        sampled_nx = (nx + stride - 1) // stride
        sampled_ny = (ny + stride - 1) // stride
        num_render = sampled_nx * sampled_ny
        viz.sample_position_to_render(state.tensor_field, amp_boost, stride, num_render)
        pos_np = state.tensor_field.position_render.to_numpy()[:num_render]
        render.scene.particles(pos_np, granule_radius, color=colormap.COLOR_MEDIUM[1])


# ================================================================
# MAIN LOOP
# ================================================================


def main():
    """Main entry point for xperiment launcher."""
    selected_xperiment_arg = sys.argv[1] if len(sys.argv) > 1 else None

    # Initialize Taichi
    ti.init(arch=ti.gpu, log_level=ti.WARN)  # GPU preferred, suppress info logs

    # Initialize xperiment manager and state
    xperiment_mgr = XperimentManager()
    state = SimulationState()

    # Load xperiment from CLI argument or default
    default_xperiment = selected_xperiment_arg or "_topo_canonical4d"
    if default_xperiment not in xperiment_mgr.available_xperiments:
        print(f"Error: Xperiment '{default_xperiment}' not found!")
        return

    params = xperiment_mgr.load_xperiment(default_xperiment)
    if not params:
        return

    state.apply_xparameters(params)
    state.initialize_grid()
    state.compute_timestep()
    initialize_xperiment(state)

    # Initialize GGUI rendering
    render.init_UI(state.UNIVERSE_SIZE, state.TICK_SPACING, state.CAM_INIT)

    # Main rendering loop
    while render.window.running:
        render.init_scene(state.SHOW_AXIS)  # Initialize scene with lighting and camera

        # Handle ESC key for window close
        if render.window.is_pressed(ti.ui.ESCAPE):
            render.window.running = False
            break

        # Display UI overlays
        new_xperiment = display_xperiment_launcher(xperiment_mgr, state)
        display_controls(state)

        # Handle xperiment switching via process replacement
        if new_xperiment:
            print("\n================================================================")
            print("XPERIMENT LAUNCH")
            print(f"Now running: {new_xperiment}\n")

            sys.stdout.flush()
            sys.stderr.flush()
            render.window.running = False

            # os.execv replaces current process (macOS may show harmless warning)
            os.execv(sys.executable, [sys.executable, __file__, new_xperiment])

        # Always recompute timestep so observables have a valid dt_rs even
        # while paused (compute_energyH_density needs it for the kinetic term).
        state.compute_timestep()

        if not state.PAUSED:
            # Run simulation step and update time
            compute_propagation(state)
            state.elapsed_t_rs += state.dt_rs  # Accumulate simulation time
            state.frame += 1

        # Always compute observables — flux-mesh + dashboard reflect the
        # current ψ even when paused (user can inspect the seeded state).
        compute_field_observables(state)

        # Render scene elements
        render_elements(state)

        # Display additional UI elements and scene
        display_wave_menu(state)
        display_data_dashboard(state)
        display_model_specs(state, model_bar_vertices)
        render.show_scene()

        # Capture frame for video export (stops after VIDEO_FRAMES)
        if state.EXPORT_VIDEO:
            video.export_frame(state.frame, state.VIDEO_FRAMES)

    if state.INSTRUMENTATION:
        instrument.generate_plots()


# ================================================================
# ENTRY POINT
# ================================================================
if __name__ == "__main__":
    main()
