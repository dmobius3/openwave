"""
XPERIMENT INSTRUMENTATION (data collection)

This provides zero-overhead data collection that can be toggled on/off per xperiment.
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
from pathlib import Path

from openwave.common import colormap


# ================================================================
# Module-level Constants (computed once at import)
# ================================================================

_MODULE_DIR = Path(__file__).parent
DATA_DIR = _MODULE_DIR / "data"
PLOT_DIR = _MODULE_DIR / "plots"

# Module-level state
_timestep_buffer = []
_timestep_log_initialized = False
_BUFFER_FLUSH_INTERVAL = 100  # Flush every N timesteps


# ================================================================
# Instrumentation Functions (Zero-Overhead)
# ================================================================


def plot_probe_wave_profile(tensor_field):
    """
    Plot the displacement profile along the x-axis through the probe position.

    Args:
        tensor_field: TensorField instance containing displacement data
    """

    # Define probe position
    px, py, pz = tensor_field.nx // 2, tensor_field.ny // 2, tensor_field.nz // 2

    # Extract the M-substrate "displacement" along x-axis at center (y, z).
    # M5.8 ψ-retire: the Vector(3) ψ field is gone; the director-substrate analog of
    # displacement is the in-plane tilt of the director n̂ away from the ẑ vacuum,
    # √(n_x²+n_y²) — zero in the uniaxial vacuum, rises inside a defect.
    x_indices = np.arange(tensor_field.nx)
    dn = tensor_field.director_nhat.to_numpy()
    displacements = np.hypot(dn[:, py, pz, 0], dn[:, py, pz, 1])

    distances = x_indices - px

    plt.style.use("dark_background")
    fig = plt.figure(figsize=(8, 5), facecolor=colormap.DARK_GRAY[1])
    fig.suptitle("OPENWAVE Analytics", fontsize=20, family="Monospace")

    plt.plot(
        distances,
        displacements,
        color=colormap.viridis_palette[2][1],
        linewidth=2,
        label="director tilt √(nx²+ny²)",
    )
    plt.axhline(y=0, color="w", linestyle="--", alpha=0.3)
    plt.axvline(x=0, color="r", linestyle="--", alpha=0.3)
    plt.xlabel("Distance from probe center (grid indices)", family="Monospace")
    plt.ylabel("director tilt √(nx²+ny²)", family="Monospace")
    plt.title("WAVE PROFILE", family="Monospace")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    # Save to directory
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    save_path = PLOT_DIR / "wave_profile.png"
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print("\nPlot wave_profile saved to:\n", save_path, "\n")


def log_timestep_data(timestep: int, tensor_field, trackers) -> None:
    """Record all timestep data to a buffer, flush periodically to reduce I/O overhead.

    Args:
        timestep: Current simulation timestep
        tensor_field: TensorField instance
        trackers: Trackers instance
    """
    global _timestep_buffer, _timestep_log_initialized

    # Define probe position
    px, py, pz = tensor_field.nx // 2, tensor_field.ny // 2, tensor_field.nz // 2

    # Capture probe values
    # M5.8 ψ-retire: ψ is gone; record the director in-plane tilt √(nx²+ny²) as the
    # M-substrate displacement analog (the ‖M−D‖_F amplitude is logged below too).
    dn = tensor_field.director_nhat[px, py, pz]
    displacement_am = float(np.hypot(dn[0], dn[1]))
    amp_local_emarms_am = float(trackers.amp_local_emarms_am[px, py, pz])
    freq_local_cross_rHz = float(trackers.freq_local_cross_rHz[px, py, pz])

    # Add to buffer
    _timestep_buffer.append(
        [
            timestep,
            displacement_am,
            amp_local_emarms_am,
            freq_local_cross_rHz,
        ]
    )

    # Flush buffer periodically
    if len(_timestep_buffer) >= _BUFFER_FLUSH_INTERVAL:
        _flush_timestep_buffer()


def _flush_timestep_buffer() -> None:
    """Write buffered timestep data to disk."""
    global _timestep_buffer, _timestep_log_initialized

    if not _timestep_buffer:
        return

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    log_path = DATA_DIR / "timestep_data.csv"

    # Write header on first flush
    if not _timestep_log_initialized:
        with open(log_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "timestep",
                    "displacement_am",
                    "amp_local_emarms_am",
                    "freq_local_cross_rHz",
                ]
            )
        _timestep_log_initialized = True

    # Append all buffered rows at once
    with open(log_path, "a", newline="") as f:
        writer = csv.writer(f)
        for row in _timestep_buffer:
            writer.writerow(
                [
                    row[0],
                    f"{row[1]:.6f}",
                    f"{row[2]:.6f}",
                    f"{row[3]:.6f}",
                ]
            )

    _timestep_buffer = []


def _read_timestep_data():
    """Read timestep data from consolidated CSV file.

    Returns:
        dict: Dictionary with lists for each column, or None if file doesn't exist
    """
    log_path = DATA_DIR / "timestep_data.csv"
    if not log_path.exists():
        print("\nTimestep data log file does not exist.\n")
        return None

    data = {
        "timesteps": [],
        "displacements": [],
        "amplitudes": [],
        "frequencies": [],
    }

    with open(log_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data["timesteps"].append(int(row["timestep"]))
            data["displacements"].append(float(row["displacement_am"]))
            data["amplitudes"].append(float(row["amp_local_emarms_am"]))
            data["frequencies"].append(float(row["freq_local_cross_rHz"]))

    return data


def plot_probe_values():
    """Plot the logged displacement, amplitude, and frequency over time."""
    data = _read_timestep_data()
    if data is None:
        return

    # M5.0d.3: dropped EWT-scaled reference axhlines (eWAVE_AMPLITUDE /
    # eWAVE_FREQUENCY) — they don't apply to M5's variable-λ regime, where the
    # reference scale is xperiment-driven (seed wavelength or defect λ_C).
    # Also dropped the broken transverse subplot (data["displacements_T"] was
    # never populated; M5 doesn't decompose L/T at the storage layer).
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(9, 7), facecolor=colormap.DARK_GRAY[1])
    fig.suptitle("OPENWAVE Analytics", fontsize=20, family="Monospace")

    # Plot 1: director-tilt displacement + RMS amplitude over time (top)
    plt.subplot(2, 1, 1)
    plt.plot(
        data["timesteps"],
        data["displacements"],
        color=colormap.viridis_palette[2][1],
        linewidth=2,
        label="director tilt",
    )
    plt.plot(
        data["timesteps"],
        data["amplitudes"],
        color=colormap.viridis_palette[3][1],
        linewidth=2,
        label="RMS AMPLITUDE (am)",
    )
    plt.axhline(y=0, color="w", linestyle="--", alpha=0.3)
    plt.xlabel("Timestep", family="Monospace")
    plt.ylabel("Displacement / Amplitude (am)", family="Monospace")
    plt.title("DISPLACEMENT & AMPLITUDE OVER TIME", family="Monospace")
    plt.grid(True, alpha=0.3)
    plt.legend()

    # Plot 2: Frequency over time (bottom)
    plt.subplot(2, 1, 2)
    plt.plot(
        data["timesteps"],
        data["frequencies"],
        color=colormap.blueprint_palette[2][1],
        linewidth=2,
        label="FREQUENCY (rHz)",
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
    # Flush any remaining buffered data before plotting
    _flush_timestep_buffer()
    plot_probe_values()
    plt.show()


if __name__ == "__main__":
    plot_probe_values()
    plt.show()
