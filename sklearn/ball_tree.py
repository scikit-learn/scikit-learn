import warnings
from .neighbors import BallTree

warnings.warn("BallTree has been moved to sklearn.neighbors.BallTree in v0.9 "
              "sklearn.ball_tree will be removed in v0.11",
              category=DeprecationWarning)
