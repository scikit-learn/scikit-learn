# TODO: Remove in 1.2
import warnings

from ..metrics import DistanceMetric
from .. import neighbors

warnings.warn(
    "sklearn.neighbors.DistanceMetric has been moved "
    "to sklearn.metrics.DistanceMetric in 1.0. "
    "This import path will be removed in 1.2",
    category=FutureWarning,
)

# Monkey-patching neighbors to alias sklearn.metrics.DistanceMetric
setattr(neighbors, "DistanceMetric", DistanceMetric)
neighbors.__all__ += ["DistanceMetric"]
