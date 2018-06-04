import numpy as np
from ._fortran import *
from scipy._lib._version import NumpyVersion


# Don't use deprecated Numpy C API.  Define this to a fixed version instead of
# NPY_API_VERSION in order not to break compilation for released Scipy versions
# when Numpy introduces a new deprecation.  Use in setup.py::
#
#   config.add_extension('_name', sources=['source_fname'], **numpy_nodepr_api)
#
if NumpyVersion(np.__version__) >= '1.10.0.dev':
    numpy_nodepr_api = dict(define_macros=[("NPY_NO_DEPRECATED_API",
                                            "NPY_1_9_API_VERSION")])
else:
    numpy_nodepr_api = dict()


from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
