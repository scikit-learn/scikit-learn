import warnings

from .cross_decomposition import CCA, PLSSVD, PLSRegression, PLSCanonical


warnings.warn("This module has been moved to cross_decomposition and will be "
              "removed in 0.16", DeprecationWarning)
