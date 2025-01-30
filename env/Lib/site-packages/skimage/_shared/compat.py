"""Compatibility helpers for dependencies."""

from packaging.version import parse

import numpy as np
import scipy as sp


__all__ = [
    "NP_COPY_IF_NEEDED",
    "SCIPY_CG_TOL_PARAM_NAME",
]


NUMPY_LT_2_0_0 = parse(np.__version__) < parse('2.0.0.dev0')

# With NumPy 2.0.0, `copy=False` now raises a ValueError if the copy cannot be
# made. The previous behavior to only copy if needed is provided with `copy=None`.
# During the transition period, use this symbol instead.
# Remove once NumPy 2.0.0 is the minimal required version.
# https://numpy.org/devdocs/release/2.0.0-notes.html#new-copy-keyword-meaning-for-array-and-asarray-constructors
# https://github.com/numpy/numpy/pull/25168
NP_COPY_IF_NEEDED = False if NUMPY_LT_2_0_0 else None


SCIPY_LT_1_12 = parse(sp.__version__) < parse('1.12')

# Starting in SciPy v1.12, 'scipy.sparse.linalg.cg' keyword argument `tol` is
# deprecated in favor of `rtol`.
SCIPY_CG_TOL_PARAM_NAME = "tol" if SCIPY_LT_1_12 else "rtol"
