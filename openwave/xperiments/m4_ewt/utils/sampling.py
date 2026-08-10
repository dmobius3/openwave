import json
from pathlib import Path

_plot_timesteps = []
_plot_displacements = []
_plot_amplitudes = []
_plot_frequencies = []
_stability_timesteps = []
_stability_drifts = []
_stability_active = []

# File shared with the external live-monitor process
# The model root, one level above utils/. This must match the path the launcher
# hands to the viewer process, or the monitor reads a file nobody writes.
MONITOR_DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "_live_monitor_data.json"
SAVE_EVERY = 10  # write the shared file every 10 samples


def sample_stability_metrics(timestep, mean_drift, active_wc):
    _stability_timesteps.append(timestep)
    _stability_drifts.append(mean_drift if mean_drift is not None else 0.0)
    _stability_active.append(active_wc)

    if len(_stability_timesteps) % SAVE_EVERY == 0:
        save_monitor_data()


def sample_for_plots(timestep, wave_field, trackers):
    px = wave_field.nx // 2 + 1
    py = wave_field.ny // 2
    pz = wave_field.nz // 2

    disp = wave_field.psi_am[px, py, pz]
    amp = trackers.amp_local_emarms_am[px, py, pz]
    freq = trackers.freq_local_cross_rHz[px, py, pz]

    def to_float(val):
        try:
            if hasattr(val, "__len__"):
                return float(val[0])
            return float(val)
        except Exception:
            return 0.0

    d = to_float(disp) / wave_field.scale_factor
    a = to_float(amp) / wave_field.scale_factor
    f = to_float(freq) * wave_field.scale_factor

    _plot_timesteps.append(timestep)
    _plot_displacements.append(d)
    _plot_amplitudes.append(a)
    _plot_frequencies.append(f)

    if len(_plot_timesteps) % SAVE_EVERY == 0:
        save_monitor_data()


def save_monitor_data():
    data = {
        "timesteps": _plot_timesteps,
        "displacements": _plot_displacements,
        "amplitudes": _plot_amplitudes,
        "frequencies": _plot_frequencies,
        "stability_timesteps": _stability_timesteps,
        "mean_drifts": _stability_drifts,
        "active_wcs": _stability_active,
    }
    MONITOR_DATA_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(MONITOR_DATA_PATH, "w") as f:
        json.dump(data, f)


def get_plot_data():
    return {
        "timesteps": _plot_timesteps,
        "displacements": _plot_displacements,
        "amplitudes": _plot_amplitudes,
        "frequencies": _plot_frequencies,
    }
