"""Physics processors. Domain layer on top of the pipeline engine."""

from .wc_types import WC, WCState
from .wc_factory import build_wc_state
from .features import (
    EMCDensityField,
    EMCFluxField,
    PsiBaseField,
    PsiLongField,
    PsiTransField,
    PsiTripleBuffer,
    TrackerFields,
    WaveGrid,
    WaveStats,
)
from .units import (
    NaturalUnitSystem,
    OpenWaveUnitSystem,
    SIUnitSystem,
    UnitSystem,
    make_unit_system,
)
from .allocator import AllocateTrackers, AllocateWaveSpeed, AllocateWaveField
from .emc import UpdateEMCDensityProcessor, UpdateWaveSpeedProcessor
from .seed import SeedPulse, SeedMultiCenter
from .evolution import LaplacianProcessor, LeapfrogProcessor
from .nonlinearity import NonlinearCubic
from .boundary import DirichletBoundaryProcessor
from .measure import AmplitudeTracker
from .visualize import TaichiWindowProcessor

__all__ = [
    # wave centers
    "WC",
    "WCState",
    "build_wc_state",
    # fields
    "PsiTripleBuffer",
    "PsiBaseField",
    "PsiLongField",
    "PsiTransField",
    "EMCDensityField",
    "EMCFluxField",
    "TrackerFields",
    "WaveSpeedField",
    "WaveGrid",
    "WaveStats",
    # units
    "UnitSystem",
    "NaturalUnitSystem",
    "OpenWaveUnitSystem",
    "SIUnitSystem",
    "make_unit_system",
    # processors
    "AllocateWaveField",
    "AllocateTrackers",
    "AllocateWaveSpeed",
    "SeedPulse",
    "SeedMultiCenter",
    "ClearAccelerationProcessor",
    "LaplacianProcessor",
    "LaplacianVariableCoeffProcessor",
    "LeapfrogProcessor",
    "NonlinearCubic",
    "DirichletBoundaryProcessor",
    "AmplitudeTracker",
    "TaichiWindowProcessor",
    "UpdateEMCDensityProcessor",
    "UpdateWaveSpeedProcessor",
]
