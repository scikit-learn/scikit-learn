"""
The :mod:`sklearn.covariance` module includes methods and algorithms to
robustly estimate the covariance of features given a set of points. The
precision matrix defined as the inverse of the covariance is also estimated.
Covariance estimation is closely related to the theory of Gaussian Graphical
Models.
"""
from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={
        "_elliptic_envelope": ["EllipticEnvelope"],
        "_empirical_covariance": [
            "log_likelihood",
            "EmpiricalCovariance",
            "empirical_covariance",
        ],
        "_graph_lasso": ["GraphicalLassoCV", "GraphicalLasso", "graphical_lasso"],
        "_robust_covariance": ["MinCovDet", "fast_mcd"],
        "_shrunk_covariance": [
            "ShrunkCovariance",
            "ledoit_wolf_shrinkage",
            "OAS",
            "ledoit_wolf",
            "LedoitWolf",
            "oas",
            "shrunk_covariance",
        ],
    },
)
