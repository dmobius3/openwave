"""
XPERIMENT INSTRUMENTATION (data collection)

This provides zero-overhead data collection that can be toggled on/off
per xperiment.  Timestep data is stored as JSON via the common
openwave.common.json_logger module.
"""

from pathlib import Path

from openwave.common import json_logger

# ================================================================
# Module-level directories
# ================================================================
# The model root, one level above utils/, so data and plots stay where they always were.
_MODULE_DIR = Path(__file__).resolve().parent.parent


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
    # Offset probe by +1 voxel along X to avoid the central Dirichlet node
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
                positions.append(None)
        log_entry["wc_positions"] = positions

    json_logger.log_timestep(log_entry)

    import numpy as np


_pairwise_ref = None
_pairwise_ref_set = False


def log_stability_metrics(timestep: int, wave_center) -> tuple:
    """
    Log WC stability metrics: mean pairwise distance drift and active WC count.
    Returns (mean_drift, n_active) so that the caller can feed the live monitor.
    """
    import numpy as np

    global _pairwise_ref, _pairwise_ref_set

    positions = []
    for i in range(wave_center.num_sources):
        if wave_center.active[i] == 0:
            continue
        pos = wave_center.position_float[i]
        x = float(pos[0])
        y = float(pos[1])
        z = float(pos[2])
        if any(np.isnan([x, y, z])):
            continue
        positions.append([x, y, z])

    n_active = len(positions)

    if n_active < 2:
        json_logger.log_timestep(
            {
                "timestep": timestep,
                "mean_drift": None,
                "active_wc": n_active,
            }
        )
        return None, n_active

    positions_np = np.array(positions, dtype=np.float64)
    diff = positions_np[:, np.newaxis, :] - positions_np[np.newaxis, :, :]
    dist = np.sqrt(np.sum(diff**2, axis=-1))

    if not _pairwise_ref_set:
        _pairwise_ref = dist.copy()
        _pairwise_ref_set = True
        mean_drift = 0.0
    else:
        drift_matrix = np.abs(dist - _pairwise_ref)
        i_upper = np.triu_indices_from(drift_matrix, k=1)
        mean_drift = float(np.mean(drift_matrix[i_upper]))

    json_logger.log_timestep(
        {
            "timestep": timestep,
            "mean_drift": mean_drift,
            "active_wc": n_active,
        }
    )
    return mean_drift, n_active
