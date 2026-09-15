"""Physics processors. Domain layer on top of the pipeline engine."""
from .wc_types import WC, WCState
from .wc_factory import build_wc_state
from .features import PsiField, WaveGrid, WaveStats
from .allocator import AllocateWaveField
from .seed import SeedPulse, SeedMultiCenter
from .evolution import LaplacianProcessor, LeapfrogProcessor
from .nonlinearity import NonlinearCubic
from .boundary import DirichletBoundaryProcessor
from .measure import AmplitudeTracker
from .visualize import TaichiWindowProcessor

__all__ = [
    "WC", "WCState", "build_wc_state",
    "PsiField", "WaveGrid", "WaveStats",
    "AllocateWaveField",
    "SeedPulse", "SeedMultiCenter",
    "LaplacianProcessor", "LeapfrogProcessor",
    "NonlinearCubic",
    "DirichletBoundaryProcessor",
    "AmplitudeTracker",
    "TaichiWindowProcessor",
]