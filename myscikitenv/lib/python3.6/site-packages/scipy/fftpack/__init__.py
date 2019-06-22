"""
==================================================
Discrete Fourier transforms (:mod:`scipy.fftpack`)
==================================================

Fast Fourier Transforms (FFTs)
==============================

.. autosummary::
   :toctree: generated/

   fft - Fast (discrete) Fourier Transform (FFT)
   ifft - Inverse FFT
   fft2 - Two dimensional FFT
   ifft2 - Two dimensional inverse FFT
   fftn - n-dimensional FFT
   ifftn - n-dimensional inverse FFT
   rfft - FFT of strictly real-valued sequence
   irfft - Inverse of rfft
   dct - Discrete cosine transform
   idct - Inverse discrete cosine transform
   dctn - n-dimensional Discrete cosine transform
   idctn - n-dimensional Inverse discrete cosine transform
   dst - Discrete sine transform
   idst - Inverse discrete sine transform
   dstn - n-dimensional Discrete sine transform
   idstn - n-dimensional Inverse discrete sine transform

Differential and pseudo-differential operators
==============================================

.. autosummary::
   :toctree: generated/

   diff - Differentiation and integration of periodic sequences
   tilbert - Tilbert transform:         cs_diff(x,h,h)
   itilbert - Inverse Tilbert transform: sc_diff(x,h,h)
   hilbert - Hilbert transform:         cs_diff(x,inf,inf)
   ihilbert - Inverse Hilbert transform: sc_diff(x,inf,inf)
   cs_diff - cosh/sinh pseudo-derivative of periodic sequences
   sc_diff - sinh/cosh pseudo-derivative of periodic sequences
   ss_diff - sinh/sinh pseudo-derivative of periodic sequences
   cc_diff - cosh/cosh pseudo-derivative of periodic sequences
   shift - Shift periodic sequences

Helper functions
================

.. autosummary::
   :toctree: generated/

   fftshift - Shift the zero-frequency component to the center of the spectrum
   ifftshift - The inverse of `fftshift`
   fftfreq - Return the Discrete Fourier Transform sample frequencies
   rfftfreq - DFT sample frequencies (for usage with rfft, irfft)
   next_fast_len - Find the optimal length to zero-pad an FFT for speed

Note that ``fftshift``, ``ifftshift`` and ``fftfreq`` are numpy functions
exposed by ``fftpack``; importing them from ``numpy`` should be preferred.

Convolutions (:mod:`scipy.fftpack.convolve`)
============================================

.. module:: scipy.fftpack.convolve

.. autosummary::
   :toctree: generated/

   convolve
   convolve_z
   init_convolution_kernel
   destroy_convolve_cache

"""

# List of possibly useful functions in scipy.fftpack._fftpack:
#   drfft
#   zfft
#   zrfft
#   zfftnd
#   destroy_drfft_cache
#   destroy_zfft_cache
#   destroy_zfftnd_cache

from __future__ import division, print_function, absolute_import


__all__ = ['fft','ifft','fftn','ifftn','rfft','irfft',
           'fft2','ifft2',
           'diff',
           'tilbert','itilbert','hilbert','ihilbert',
           'sc_diff','cs_diff','cc_diff','ss_diff',
           'shift',
           'fftfreq', 'rfftfreq',
           'fftshift', 'ifftshift',
           'next_fast_len',
           ]

from .basic import *
from .pseudo_diffs import *
from .helper import *

from numpy.dual import register_func
for k in ['fft', 'ifft', 'fftn', 'ifftn', 'fft2', 'ifft2']:
    register_func(k, eval(k))
del k, register_func

from .realtransforms import *
__all__.extend(['dct', 'idct', 'dst', 'idst', 'dctn', 'idctn', 'dstn',
                'idstn'])

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
