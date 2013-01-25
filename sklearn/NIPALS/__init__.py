"""
The :mod:`sklearn.NIPALS` module includes several different projection based
latent variable methods that all are computed using the NIPALS algorithm.
"""

from .NIPALS import PCA

__all__ = ['PCA']
