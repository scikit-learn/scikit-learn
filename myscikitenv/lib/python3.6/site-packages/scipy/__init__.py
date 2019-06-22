"""
SciPy: A scientific computing package for Python
================================================

Documentation is available in the docstrings and
online at https://docs.scipy.org.

Contents
--------
SciPy imports all the functions from the NumPy namespace, and in
addition provides:

Subpackages
-----------
Using any of these subpackages requires an explicit import.  For example,
``import scipy.cluster``.

::

 cluster                      --- Vector Quantization / Kmeans
 fftpack                      --- Discrete Fourier Transform algorithms
 integrate                    --- Integration routines
 interpolate                  --- Interpolation Tools
 io                           --- Data input and output
 linalg                       --- Linear algebra routines
 linalg.blas                  --- Wrappers to BLAS library
 linalg.lapack                --- Wrappers to LAPACK library
 misc                         --- Various utilities that don't have
                                  another home.
 ndimage                      --- n-dimensional image package
 odr                          --- Orthogonal Distance Regression
 optimize                     --- Optimization Tools
 signal                       --- Signal Processing Tools
 signal.windows               --- Window functions
 sparse                       --- Sparse Matrices
 sparse.linalg                --- Sparse Linear Algebra
 sparse.linalg.dsolve         --- Linear Solvers
 sparse.linalg.dsolve.umfpack --- :Interface to the UMFPACK library:
                                  Conjugate Gradient Method (LOBPCG)
 sparse.linalg.eigen          --- Sparse Eigenvalue Solvers
 sparse.linalg.eigen.lobpcg   --- Locally Optimal Block Preconditioned
                                  Conjugate Gradient Method (LOBPCG)
 spatial                      --- Spatial data structures and algorithms
 special                      --- Special functions
 stats                        --- Statistical Functions

Utility tools
-------------
::

 test              --- Run scipy unittests
 show_config       --- Show scipy build configuration
 show_numpy_config --- Show numpy build configuration
 __version__       --- SciPy version string
 __numpy_version__ --- Numpy version string

"""
from __future__ import division, print_function, absolute_import

__all__ = ['test']

from numpy import show_config as show_numpy_config
if show_numpy_config is None:
    raise ImportError(
        "Cannot import scipy when running from numpy source directory.")
from numpy import __version__ as __numpy_version__

# Import numpy symbols to scipy name space
import numpy as _num
linalg = None
from numpy import *
from numpy.random import rand, randn
from numpy.fft import fft, ifft
from numpy.lib.scimath import *

# Allow distributors to run custom init code
from . import _distributor_init

__all__ += _num.__all__
__all__ += ['randn', 'rand', 'fft', 'ifft']

del _num
# Remove the linalg imported from numpy so that the scipy.linalg package can be
# imported.
del linalg
__all__.remove('linalg')

# We first need to detect if we're being called as part of the scipy
# setup procedure itself in a reliable manner.
try:
    __SCIPY_SETUP__
except NameError:
    __SCIPY_SETUP__ = False


if __SCIPY_SETUP__:
    import sys as _sys
    _sys.stderr.write('Running from scipy source directory.\n')
    del _sys
else:
    try:
        from scipy.__config__ import show as show_config
    except ImportError:
        msg = """Error importing scipy: you cannot import scipy while
        being in scipy source directory; please exit the scipy source
        tree first, and relaunch your python interpreter."""
        raise ImportError(msg)

    from scipy.version import version as __version__
    from scipy._lib._version import NumpyVersion as _NumpyVersion
    if _NumpyVersion(__numpy_version__) < '1.13.3':
        import warnings
        warnings.warn("Numpy 1.13.3 or above is required for this version of "
                      "scipy (detected version %s)" % __numpy_version__,
                      UserWarning)

    del _NumpyVersion

    from scipy._lib._ccallback import LowLevelCallable

    from scipy._lib._testutils import PytestTester
    test = PytestTester(__name__)
    del PytestTester
