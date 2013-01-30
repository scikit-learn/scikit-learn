"""
The :mod:`sklearn.NIPALS` module includes several different projection based
latent variable methods that all are computed using the NIPALS algorithm.
"""

from .NIPALS import PCA
from .NIPALS import SVD
from .NIPALS import PLSR
from .NIPALS import PLSC
from .NIPALS import center
from .NIPALS import scale
from .NIPALS import direct

__all__ = ['PCA', 'SVD', 'PLSR', 'PLSC', 'center', 'scale', 'direct']
