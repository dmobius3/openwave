"""
Interactive live monitor window for OpenWave simulations.

Displays displacement, amplitude, and frequency with automatic refresh
every N updates. Uses non-interactive mode to avoid GIL conflicts.
"""

import matplotlib.pyplot as plt
from openwave.common import colormap, constants


class LiveMonitor:
    """Auto-refreshing plot window for monitoring simulation variables."""

    def __init__(self, max_history=300, refresh_every=5):
        self.max_history = max_history
        self.refresh_every = refresh_every
        self.update_counter = 0
        self.timesteps = []
        self.displacements = []
        self.amplitudes = []
        self.frequencies = []

        plt.ioff()
        plt.style.use("dark_background")
        self.fig, (self.ax1, self.ax2) = plt.subplots(
            2, 1, figsize=(10, 7), facecolor=colormap.DARK_GRAY[1]
        )
        self.fig.suptitle("OPENWAVE Live Monitor", fontsize=20, family="Monospace")

        # Subplot 1: Displacement + Amplitude
        (self.line_disp,) = self.ax1.plot(
            [], [], color=colormap.viridis_palette[2][1], linewidth=2, label="DISPLACEMENT (am)"
        )
        (self.line_amp,) = self.ax1.plot(
            [], [], color=colormap.viridis_palette[3][1], linewidth=2, label="RMS AMPLITUDE (am)"
        )
        self.ax1.axhline(
            y=constants.EWAVE_AMPLITUDE / constants.ATTOMETER,
            color=colormap.viridis_palette[4][1],
            linestyle="--",
            alpha=0.5,
            label="eWAVE AMPLITUDE (am)",
        )
        self.ax1.axhline(y=0, color="w", linestyle="--", alpha=0.3)
        self.ax1.set_xlabel("Timestep", family="Monospace")
        self.ax1.set_ylabel("Displacement / Amplitude (am)", family="Monospace")
        self.ax1.set_title("(LONGITUDINAL) DISPLACEMENT & AMPLITUDE", family="Monospace")
        self.ax1.grid(True, alpha=0.3)
        self.ax1.legend(loc="upper right")
        self.ax1.set_ylim(auto=True)

        # Subplot 2: Frequency
        (self.line_freq,) = self.ax2.plot(
            [], [], color=colormap.blueprint_palette[2][1], linewidth=2, label="FREQUENCY (rHz)"
        )
        self.ax2.axhline(
            y=constants.EWAVE_FREQUENCY * constants.RONTOSECOND,
            color=colormap.blueprint_palette[1][1],
            linestyle="--",
            alpha=0.5,
            label="eWAVE FREQUENCY (rHz)",
        )
        self.ax2.axhline(y=0, color="w", linestyle="--", alpha=0.3)
        self.ax2.set_xlabel("Timestep", family="Monospace")
        self.ax2.set_ylabel("Frequency (rHz)", family="Monospace")
        self.ax2.set_title("FREQUENCY", family="Monospace")
        self.ax2.grid(True, alpha=0.3)
        self.ax2.legend(loc="upper right")
        self.ax2.set_ylim(auto=True)

        plt.tight_layout()
        self.fig.show()
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    def update(self, timestep, displacement, amplitude, frequency):
        """Push one sample and refresh the plot window periodically."""
        self.timesteps.append(timestep)
        self.displacements.append(displacement)
        self.amplitudes.append(amplitude)
        self.frequencies.append(frequency)

        # Trim to sliding window
        if len(self.timesteps) > self.max_history:
            self.timesteps = self.timesteps[-self.max_history :]
            self.displacements = self.displacements[-self.max_history :]
            self.amplitudes = self.amplitudes[-self.max_history :]
            self.frequencies = self.frequencies[-self.max_history :]

        self.update_counter += 1
        if self.update_counter % self.refresh_every == 0:
            self._redraw()

    def _redraw(self):
        """Update line data and refresh the canvas."""
        if not self.timesteps:
            return

        self.line_disp.set_data(self.timesteps, self.displacements)
        self.line_amp.set_data(self.timesteps, self.amplitudes)
        self.line_freq.set_data(self.timesteps, self.frequencies)

        x_min = self.timesteps[0]
        x_max = self.timesteps[-1]
        if x_max <= x_min:
            x_max = x_min + 1
        self.ax1.set_xlim(x_min, x_max)
        self.ax2.set_xlim(x_min, x_max)

        self.ax1.relim()
        self.ax1.autoscale_view(scaley=True)
        self.ax2.relim()
        self.ax2.autoscale_view(scaley=True)

        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()

    def close(self):
        """Close the live window."""
        plt.close(self.fig)


_live_monitor = None


def init_live_monitor(max_history=300, refresh_every=5):
    global _live_monitor
    _live_monitor = LiveMonitor(max_history=max_history, refresh_every=refresh_every)


def update_live_monitor(timestep, displacement, amplitude, frequency):
    if _live_monitor is not None:
        _live_monitor.update(timestep, displacement, amplitude, frequency)


def close_live_monitor():
    global _live_monitor
    if _live_monitor is not None:
        _live_monitor.close()
        _live_monitor = None
