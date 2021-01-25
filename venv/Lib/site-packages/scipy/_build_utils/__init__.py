import os

import numpy as np
from ._fortran import *
from .system_info import combine_dict


# Don't use the deprecated NumPy C API. Define this to a fixed version instead of
# NPY_API_VERSION in order not to break compilation for released SciPy versions
# when NumPy introduces a new deprecation. Use in setup.py::
#
#   config.add_extension('_name', sources=['source_fname'], **numpy_nodepr_api)
#
numpy_nodepr_api = dict(define_macros=[("NPY_NO_DEPRECATED_API",
                                        "NPY_1_9_API_VERSION")])


def uses_blas64():
    return (os.environ.get("NPY_USE_BLAS_ILP64", "0") != "0")


from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
