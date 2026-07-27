#!/usr/bin/env python3
"""
External live monitor viewer. Run as a separate process.
Reads shared JSON file and updates a matplotlib plot.
"""
import sys
import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
from openwave.common import colormap, constants


def main(data_path_str):
    data_path = Path(data_path_str)

    # Wait for the data file to appear (up to 30 seconds)
    waited = 0
    while not data_path.exists() and waited < 30:
        print(f"[monitor] Waiting for data file... ({waited}s)")
        print(f"[monitor] Data file expected at: {data_path}")
        time.sleep(1)
        waited += 1

    if not data_path.exists():
        print(f"[monitor] Data file did not appear after 30s: {data_path}")
        return

    plt.ion()
    plt.style.use("dark_background")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), facecolor=colormap.DARK_GRAY[1])
    fig.suptitle("OPENWAVE Live Monitor (external)", fontsize=20, family="Monospace")

    (line_disp,) = ax1.plot(
        [], [], color=colormap.viridis_palette[2][1], linewidth=2, label="DISPLACEMENT (am)"
    )
    (line_amp,) = ax1.plot(
        [], [], color=colormap.viridis_palette[3][1], linewidth=2, label="RMS AMPLITUDE (am)"
    )
    ax1.axhline(
        y=constants.EWAVE_AMPLITUDE / constants.ATTOMETER,
        color=colormap.viridis_palette[4][1],
        linestyle="--",
        alpha=0.5,
        label="eWAVE AMPLITUDE (am)",
    )
    ax1.axhline(y=0, color="w", linestyle="--", alpha=0.3)
    ax1.set_xlabel("Timestep", family="Monospace")
    ax1.set_ylabel("Displacement / Amplitude (am)", family="Monospace")
    ax1.set_title("(LONGITUDINAL) DISPLACEMENT & AMPLITUDE", family="Monospace")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper right")
    ax1.set_ylim(auto=True)

    (line_freq,) = ax2.plot(
        [], [], color=colormap.blueprint_palette[2][1], linewidth=2, label="FREQUENCY (rHz)"
    )
    ax2.axhline(
        y=constants.EWAVE_FREQUENCY * constants.RONTOSECOND,
        color=colormap.blueprint_palette[1][1],
        linestyle="--",
        alpha=0.5,
        label="eWAVE FREQUENCY (rHz)",
    )
    ax2.axhline(y=0, color="w", linestyle="--", alpha=0.3)
    ax2.set_xlabel("Timestep", family="Monospace")
    ax2.set_ylabel("Frequency (rHz)", family="Monospace")
    ax2.set_title("FREQUENCY", family="Monospace")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper right")
    ax2.set_ylim(auto=True)

    plt.tight_layout()
    fig.show()

    while data_path.exists():
        try:
            with open(data_path, "r") as f:
                data = json.load(f)
            ts = data.get("timesteps", [])
            if ts:
                line_disp.set_data(ts, data.get("displacements", []))
                line_amp.set_data(ts, data.get("amplitudes", []))
                line_freq.set_data(ts, data.get("frequencies", []))
                ax1.set_xlim(ts[0], max(ts[-1], ts[0] + 1))
                ax2.set_xlim(ts[0], max(ts[-1], ts[0] + 1))
                ax1.relim()
                ax1.autoscale_view(scaley=True)
                ax2.relim()
                ax2.autoscale_view(scaley=True)
                fig.canvas.draw()
        except:
            pass
        plt.pause(0.5)

    print("[monitor] Data file removed - exiting.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python live_monitor_viewer.py <path_to_data.json>")
        sys.exit(1)
    main(sys.argv[1])
