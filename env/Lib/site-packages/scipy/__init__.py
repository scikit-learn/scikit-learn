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
 fft                          --- Discrete Fourier transforms
 fftpack                      --- Legacy discrete Fourier transforms
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

# Import numpy symbols to scipy name space (DEPRECATED)
from ._lib.deprecation import _deprecated
import numpy as _num
linalg = None
_msg = ('scipy.{0} is deprecated and will be removed in SciPy 2.0.0, '
        'use numpy.{0} instead')
# deprecate callable objects, skipping classes
for _key in _num.__all__:
    _fun = getattr(_num, _key)
    if callable(_fun) and not isinstance(_fun, type):
        _fun = _deprecated(_msg.format(_key))(_fun)
    globals()[_key] = _fun
from numpy.random import rand, randn
_msg = ('scipy.{0} is deprecated and will be removed in SciPy 2.0.0, '
        'use numpy.random.{0} instead')
rand = _deprecated(_msg.format('rand'))(rand)
randn = _deprecated(_msg.format('randn'))(randn)
from numpy.fft import fft, ifft
# fft is especially problematic, so we deprecate it with a shorter window
fft_msg = ('Using scipy.fft as a function is deprecated and will be '
           'removed in SciPy 1.5.0, use scipy.fft.fft instead.')
# for wrapping in scipy.fft.__call__, so the stacklevel is one off from the
# usual (2)
_dep_fft = _deprecated(fft_msg, stacklevel=3)(fft)
fft = _deprecated(fft_msg)(fft)
ifft = _deprecated('scipy.ifft is deprecated and will be removed in SciPy '
                   '2.0.0, use scipy.fft.ifft instead')(ifft)
import numpy.lib.scimath as _sci
_msg = ('scipy.{0} is deprecated and will be removed in SciPy 2.0.0, '
        'use numpy.lib.scimath.{0} instead')
for _key in _sci.__all__:
    _fun = getattr(_sci, _key)
    if callable(_fun):
        _fun = _deprecated(_msg.format(_key))(_fun)
    globals()[_key] = _fun

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

    # This makes "from scipy import fft" return scipy.fft, not np.fft
    del fft
    from . import fft
