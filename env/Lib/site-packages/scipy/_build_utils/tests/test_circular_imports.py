from __future__ import division, absolute_import, print_function

import sys
import subprocess


PUBLIC_SUBMODULES = [
    'cluster',
    'cluster.hierarchy',
    'cluster.vq',
    'constants',
    'fft',
    'fftpack',
    'fftpack.convolve',
    'integrate',
    'interpolate',
    'io',
    'io.arff',
    'io.wavfile',
    'linalg',
    'linalg.blas',
    'linalg.lapack',
    'linalg.interpolative',
    'misc',
    'ndimage',
    'odr',
    'optimize',
    'signal',
    'sparse',
    'sparse.csgraph',
    'sparse.linalg',
    'spatial',
    'spatial.distance',
    'special',
    'stats',
    'stats.mstats',
]


def test_importing_submodules():
    # Regression test for gh-6793.
    for name in PUBLIC_SUBMODULES:
        try:
            cmd = [sys.executable, '-c', 'import scipy.{0}'.format(name)]
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError:
            raise AssertionError('Importing scipy.{0} failed'.format(name))
