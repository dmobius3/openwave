"""
XPERIMENT INSTRUMENTATION (data collection)

This provides zero-overhead data collection that can be toggled on/off
per xperiment.  Timestep data is stored as JSON via the common
openwave.common.json_logger module.
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
from pathlib import Path

from openwave.common import colormap, constants
from openwave.common import json_logger

# ================================================================
# Module-level directories
# ================================================================
_MODULE_DIR = Path(__file__).parent
PLOT_DIR = _MODULE_DIR / "plots"


# ================================================================
# Initialisation (called from launcher)
# ================================================================

def init_instrumentation(state, xperiment_name="unknown", data_dir=None):
    """
    Build the metadata dict and start a json_logger session.
    Must be called once before the main simulation loop.
    """
    if data_dir is None:
        data_dir = _MODULE_DIR / "data"

    meta = {
        "model": "M4",
        "xperiment": xperiment_name,
        "engine": {
            "SEED_MODE": state.SEED_MODE,
            "SEED_BOOST": state.SEED_BOOST,
            "V_MODE": state.V_MODE,
            "V_C1": state.V_C1,
            "V_C2": state.V_C2,
            "WC_INTERACT_MODE": state.WC_INTERACT_MODE,
            "WC_BOOST": state.WC_BOOST,
            "WC_RADIUS": state.WC_RADIUS,
            "WC_SIGMA": state.WC_SIGMA,
        },
        "wave_centers": {
            "K": state.NUM_SOURCES,
            "POSITIONS": state.SOURCES_POSITION,
            "PHASE_OFFSETS_DEG": state.SOURCES_OFFSET_DEG,
        },
        "universe": {
            "EDGE_X": state.UNIVERSE_SIZE[0],
            "EDGE_Y": state.UNIVERSE_SIZE[1],
            "EDGE_Z": state.UNIVERSE_SIZE[2],
            "TARGET_VOXELS": state.TARGET_VOXELS,
        },
        "simulation": {
            "SIM_SPEED": state.SIM_SPEED,
            "dt_rs": state.dt_rs,
            "cfl_factor": state.cfl_factor,
            "PAUSED": state.PAUSED,
        },
        "emc_profile": {
            "R_WALL": state.R_WALL,
            "WALL_HEIGHT": state.WALL_HEIGHT,
            "DEFICIT_DEPTH": state.DEFICIT_DEPTH,
            "R_SOLITON": state.R_SOLITON,           
            "SIGMA": state.SIGMA,                  
            "PRESSURE_STRENGTH": state.PRESSURE_STRENGTH, 
        },
    }

    # Add the K value at top level for filename generation
    meta["K"] = state.NUM_SOURCES

    json_logger.init_session(meta, data_dir=Path(data_dir))


# ================================================================
# Timestep logging
# ================================================================

def log_timestep_data(timestep: int, wave_field, trackers, wave_center=None) -> None:
    """
    Log timestep data including field values and wave-center positions.

    Args:
        timestep: Current simulation frame number
        wave_field: WaveField instance
        trackers: Trackers instance
        wave_center: WaveCenter instance (optional, for WC position logging)
    """
    px, py, pz = wave_field.nx // 2, wave_field.ny // 2, wave_field.nz // 2
    disp = wave_field.psi_am[px, py, pz]
    amp = trackers.amp_local_emarms_am[px, py, pz]
    freq = trackers.freq_local_cross_rHz[px, py, pz]

    def to_float(val):
        try:
            if hasattr(val, '__len__'):
                return float(val[0])
            return float(val)
        except:
            return 0.0

    displacement_am = to_float(disp) / wave_field.scale_factor
    amp_local_emarms_am = to_float(amp) / wave_field.scale_factor
    freq_local_cross_rHz = to_float(freq) * wave_field.scale_factor

    # Build log entry
    log_entry = {
        "timestep": timestep,
        "displacement_am": displacement_am,
        "amp_local_emarms_am": amp_local_emarms_am,
        "freq_local_cross_rHz": freq_local_cross_rHz,
    }

    # Log wave-center positions if wave_center is provided
    if wave_center is not None:
        positions = []
        for i in range(wave_center.num_sources):
            if wave_center.active[i] == 1:
                pos = wave_center.position_float[i]
                positions.append([float(pos[0]), float(pos[1]), float(pos[2])])
            else:
                # Log inactive WCs as None
                positions.append(None)
        log_entry["wc_positions"] = positions

    json_logger.log_timestep(log_entry)


# ================================================================
# Plotting
# ================================================================

def plot_probe_wave_profile(wave_field):
    """
    Plot the longitudinal (x) displacement profile along the x-axis
    through the probe position.
    """
    px, py, pz = wave_field.nx // 2, wave_field.ny // 2, wave_field.nz // 2
    x_indices = np.arange(wave_field.nx)
    displacements = np.zeros(wave_field.nx)

    for i in range(wave_field.nx):
        # psi_am is a 3-vector; extract the x-component (longitudinal)
        displacements[i] = wave_field.psi_am[i, py, pz][0]

    distances = x_indices - px

    plt.style.use("dark_background")
    fig = plt.figure(figsize=(12, 6), facecolor=colormap.DARK_GRAY[1])
    fig.suptitle("OPENWAVE Analytics", fontsize=20, family="Monospace")

    plt.subplot(1, 1, 1)
    plt.plot(
        distances,
        displacements,
        color=colormap.viridis_palette[2][1],
        linewidth=4,
        label="LONGITUDINAL (X)",
    )
    plt.axhline(y=0, color="k", linestyle="--", alpha=0.3)
    plt.axvline(x=0, color="r", linestyle="--", alpha=0.3)
    plt.ylim(-1.0, 6.0)
    plt.xlabel("Distance from Wave-Center (grid indices)", family="Monospace")
    plt.ylabel("Displacement (attometers)", family="Monospace")
    plt.title("WAVE PROFILE", family="Monospace")
    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    save_path = PLOT_DIR / "wave_profile.png"
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print("\nPlot wave_profile saved to:\n", save_path, "\n")


def _read_timestep_data():
    """
    Read timestep data from JSON log for plotting purposes.
    Returns a dict suitable for the existing plotting functions.
    """
    import json
    log_path = json_logger._data_dir / json_logger._filename
    if not log_path.exists():
        print("\nTimestep data log file does not exist.\n")
        return None

    with open(log_path, "r") as f:
        doc = json.load(f)

    data = {
        "timesteps": [],
        "displacements": [],
        "amplitudes": [],
        "frequencies": [],
    }
    for rec in doc["data"]:
        data["timesteps"].append(rec["timestep"])
        data["displacements"].append(rec["displacement_am"])
        data["amplitudes"].append(rec["amp_local_emarms_am"])
        data["frequencies"].append(rec["freq_local_cross_rHz"])
    return data


def plot_probe_values():
    """Plot the logged displacement, amplitude, and frequency over time."""
    data = _read_timestep_data()
    if data is None:
        return

    # Create the plot with 3 subplots (stacked vertically)
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(9, 9), facecolor=colormap.DARK_GRAY[1])
    fig.suptitle("OPENWAVE Analytics", fontsize=20, family="Monospace")

    # Plot 1: Longitudinal Displacement and Amplitude (top)
    plt.subplot(3, 1, 1)
    plt.plot(
        data["timesteps"],
        data["displacements"],
        color=colormap.viridis_palette[2][1],
        linewidth=2,
        label="DISPLACEMENT (am)",
    )
    plt.plot(
        data["timesteps"],
        data["amplitudes"],
        color=colormap.viridis_palette[3][1],
        linewidth=2,
        label="RMS AMPLITUDE (am)",
    )
    plt.axhline(
        y=constants.EWAVE_AMPLITUDE / constants.ATTOMETER,
        color=colormap.viridis_palette[4][1],
        linestyle="--",
        alpha=0.5,
        label="eWAVE AMPLITUDE (am)",
    )
    plt.axhline(y=0, color="w", linestyle="--", alpha=0.3)
    plt.xlabel("Timestep", family="Monospace")
    plt.ylabel("Displacement / Amplitude (am)", family="Monospace")
    plt.title("(LONGITUDINAL) DISPLACEMENT & AMPLITUDE OVER TIME", family="Monospace")
    plt.grid(True, alpha=0.3)
    plt.legend()

    # Plot 2: Transverse Displacement and Amplitude (middle)
    plt.subplot(3, 1, 2)
    plt.plot(
        data["timesteps"],
        data["displacements_T"],
        color=colormap.ironbow_palette[2][1],
        linewidth=2,
        label="DISPLACEMENT (am)",
    )
    plt.plot(
        data["timesteps"],
        data["amplitudes"],
        color=colormap.ironbow_palette[3][1],
        linewidth=2,
        label="RMS AMPLITUDE (am)",
    )
    plt.axhline(
        y=constants.EWAVE_AMPLITUDE / constants.ATTOMETER,
        color=colormap.ironbow_palette[4][1],
        linestyle="--",
        alpha=0.5,
        label="eWAVE AMPLITUDE (am)",
    )
    plt.axhline(y=0, color="w", linestyle="--", alpha=0.3)
    plt.xlabel("Timestep", family="Monospace")
    plt.ylabel("Displacement / Amplitude (am)", family="Monospace")
    plt.title("(TRANSVERSE) DISPLACEMENT & AMPLITUDE OVER TIME", family="Monospace")
    plt.grid(True, alpha=0.3)
    plt.legend()

    # Plot 3: Frequency (bottom)
    plt.subplot(3, 1, 3)
    plt.plot(
        data["timesteps"],
        data["frequencies"],
        color=colormap.blueprint_palette[2][1],
        linewidth=2,
        label="FREQUENCY (rHz)",
    )
    plt.axhline(
        y=constants.EWAVE_FREQUENCY * constants.RONTOSECOND,
        color=colormap.blueprint_palette[1][1],
        linestyle="--",
        alpha=0.5,
        label="eWAVE FREQUENCY (rHz)",
    )
    plt.axhline(y=0, color="w", linestyle="--", alpha=0.3)
    plt.xlabel("Timestep", family="Monospace")
    plt.ylabel("Frequency (rHz)", family="Monospace")
    plt.title("FREQUENCY OVER TIME", family="Monospace")
    plt.grid(True, alpha=0.3)
    plt.legend()

    plt.tight_layout()

    # Save to directory
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    save_path = PLOT_DIR / "probe_values.png"
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print("\nPlot probe values saved to:\n", save_path, "\n")


def generate_plots():
    """Generate all instrumentation plots."""
    json_logger.finalize()
    plot_probe_values()
    plt.show()


def finalize_instrumentation():
    """Flush remaining records. Call at end of simulation."""
    json_logger.finalize()


if __name__ == "__main__":
    plot_probe_values()
    plt.show()