"""
Post-simulation plotting functions.

These are called once after the simulation loop ends.
"""

import matplotlib.pyplot as plt
from pathlib import Path

from openwave.common import colormap, constants
from openwave.common import json_logger
from openwave.xperiments.m4_ewt.utils import sampling

# The model root, one level above utils/, so plots stay where they always were.
PLOT_DIR = Path(__file__).resolve().parent.parent / "plots"


def _read_timestep_data():
    """Read timestep data from JSON log for plotting purposes."""
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
        if "displacement_am" not in rec:
            continue
        data["timesteps"].append(rec["timestep"])
        data["displacements"].append(rec["displacement_am"])
        data["amplitudes"].append(rec["amp_local_emarms_am"])
        data["frequencies"].append(rec["freq_local_cross_rHz"])
    return data


def plot_probe_values():
    """Plot the logged displacement, amplitude, and frequency over time (from JSON)."""
    data = _read_timestep_data()
    if data is None:
        return

    has_transverse = "displacements_T" in data and data["displacements_T"] is not None

    n_plots = 3 if has_transverse else 2
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(9, 3 * n_plots), facecolor=colormap.DARK_GRAY[1])
    fig.suptitle("OPENWAVE Analytics (JSON log)", fontsize=20, family="Monospace")

    # Plot 1: Longitudinal Displacement and Amplitude
    plt.subplot(n_plots, 1, 1)
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
    plt.ylim(bottom=0)

    if has_transverse:
        plt.subplot(n_plots, 1, 2)
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
        plt.ylim(bottom=0)

    plt.subplot(n_plots, 1, n_plots)
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
    plt.ylim(bottom=0)

    plt.tight_layout()
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    save_path = PLOT_DIR / "probe_values.png"
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print("\nPlot probe values saved to:\n", save_path, "\n")


def plot_live_values():
    """Plot the densely sampled live data (in-memory, no JSON involved)."""
    data = sampling.get_plot_data()
    if len(data["timesteps"]) < 2:
        print("\n[plot_live] Not enough live samples to plot.\n")
        return

    plt.style.use("dark_background")
    fig = plt.figure(figsize=(9, 7), facecolor=colormap.DARK_GRAY[1])
    fig.suptitle("OPENWAVE Live Monitor (final)", fontsize=20, family="Monospace")

    # Subplot 1: Displacement and Amplitude
    plt.subplot(2, 1, 1)
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
    plt.title("LIVE DISPLACEMENT & AMPLITUDE", family="Monospace")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.ylim(bottom=0)

    # Subplot 2: Frequency
    plt.subplot(2, 1, 2)
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
    plt.title("LIVE FREQUENCY", family="Monospace")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.ylim(bottom=0)

    plt.tight_layout()
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    save_path = PLOT_DIR / "live_monitor.png"
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print("\nLive monitor plot saved to:\n", save_path, "\n")


def plot_wc_drift():
    """Plot mean pairwise distance drift over time."""
    import json

    log_path = json_logger._data_dir / json_logger._filename
    if not log_path.exists():
        return
    with open(log_path) as f:
        doc = json.load(f)

    drift_data = [d for d in doc["data"] if d.get("mean_drift") is not None]
    if not drift_data:
        print("[plot_wc_drift] No drift data in log.")
        return

    ts = [d["timestep"] for d in drift_data]
    drift = [d["mean_drift"] for d in drift_data]

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(10, 4), facecolor=colormap.DARK_GRAY[1])
    ax.plot(ts, drift, color=colormap.viridis_palette[2][1], linewidth=2)
    ax.set_xlabel("Timestep")
    ax.set_ylabel("Mean Pairwise Drift [vox]")
    ax.set_title("Pairwise Distance Drift")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plt.savefig(PLOT_DIR / "wc_drift.png", dpi=150, bbox_inches="tight")
    print(f"[plot_wc_drift] Saved to {PLOT_DIR / 'wc_drift.png'}")


def plot_wc_active():
    """Plot number of active wave centers over time."""
    import json

    log_path = json_logger._data_dir / json_logger._filename
    if not log_path.exists():
        return
    with open(log_path) as f:
        doc = json.load(f)

    active_data = [d for d in doc["data"] if "active_wc" in d]
    if not active_data:
        print("[plot_wc_active] No active WC data in log.")
        return

    ts = [d["timestep"] for d in active_data]
    active = [d["active_wc"] for d in active_data]

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(10, 4), facecolor=colormap.DARK_GRAY[1])
    ax.plot(ts, active, color=colormap.ironbow_palette[2][1], linewidth=2)
    ax.set_xlabel("Timestep")
    ax.set_ylabel("Active WC count")
    ax.set_title("Active Wave Centers")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    plt.tight_layout()
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plt.savefig(PLOT_DIR / "wc_active.png", dpi=150, bbox_inches="tight")
    print(f"[plot_wc_active] Saved to {PLOT_DIR / 'wc_active.png'}")


def generate_plots():
    """Generate all instrumentation plots. Called after simulation ends."""
    json_logger.finalize()
    plot_probe_values()
    plot_live_values()
    plot_wc_drift()
    plot_wc_active()
    plt.show()
