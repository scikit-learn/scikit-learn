from __future__ import division, print_function, absolute_import

import sys
import subprocess


MODULES = [
    "scipy.cluster",
    "scipy.cluster.vq",
    "scipy.cluster.hierarchy",
    "scipy.constants",
    "scipy.fft",
    "scipy.fftpack",
    "scipy.integrate",
    "scipy.interpolate",
    "scipy.io",
    "scipy.io.arff",
    "scipy.io.harwell_boeing",
    "scipy.io.idl",
    "scipy.io.matlab",
    "scipy.io.netcdf",
    "scipy.io.wavfile",
    "scipy.linalg",
    "scipy.linalg.blas",
    "scipy.linalg.cython_blas",
    "scipy.linalg.lapack",
    "scipy.linalg.cython_lapack",
    "scipy.linalg.interpolative",
    "scipy.misc",
    "scipy.ndimage",
    "scipy.odr",
    "scipy.optimize",
    "scipy.signal",
    "scipy.signal.windows",
    "scipy.sparse",
    "scipy.sparse.linalg",
    "scipy.sparse.csgraph",
    "scipy.spatial",
    "scipy.spatial.distance",
    "scipy.special",
    "scipy.stats",
    "scipy.stats.distributions",
    "scipy.stats.mstats",
]


def test_modules_importable():
    # Check that all modules are importable in a new Python
    # process. This is not necessarily true (esp on Python 2) if there
    # are import cycles present.
    for module in MODULES:
        cmd = 'import {}'.format(module)
        subprocess.check_call([sys.executable, '-c', cmd])
