"""
==========================================
Miscellaneous routines (:mod:`scipy.misc`)
==========================================

.. currentmodule:: scipy.misc

Various utilities that don't have another home.

Note that Pillow (https://python-pillow.org/) is not a dependency
of SciPy, but the image manipulation functions indicated in the list
below are not available without it.

.. autosummary::
   :toctree: generated/

   ascent - Get example image for processing
   central_diff_weights - Weights for an n-point central m-th derivative
   derivative - Find the n-th derivative of a function at a point
   face - Get example image for processing
   electrocardiogram - Load an example of a one-dimensional signal.

Deprecated functions:

.. autosummary::
   :toctree: generated/

   bytescale - Byte scales an array (image) [requires Pillow]
   fromimage - Return a copy of a PIL image as a numpy array [requires Pillow]
   imfilter - Simple filtering of an image [requires Pillow]
   imread - Read an image file from a filename [requires Pillow]
   imresize - Resize an image [requires Pillow]
   imrotate - Rotate an image counter-clockwise [requires Pillow]
   imsave - Save an array to an image file [requires Pillow]
   imshow - Simple showing of an image through an external viewer [requires Pillow]
   toimage - Takes a numpy array and returns a PIL image [requires Pillow]


Deprecated aliases:

.. autosummary::
   :toctree: generated/

   comb - Combinations of N things taken k at a time, "N choose k" (imported from `scipy.special`)
   factorial  - The factorial function, ``n! = special.gamma(n+1)``
                (imported from `scipy.special`)
   factorial2 - Double factorial, ``(n!)!`` (imported from `scipy.special`)
   factorialk - ``(...((n!)!)!...)!`` where there are k '!' (imported from `scipy.special`)
   logsumexp - Compute the log of the sum of exponentials of input elements
               (imported from `scipy.special`)
   pade - Pade approximation to function as the ratio of two polynomials.
          (imported from `scipy.interpolate`)
   info - Get help information for a function, class, or module. (imported from `numpy`)
   source - Print function source code. (imported from `numpy`)
   who - Print the Numpy arrays in the given dictionary. (imported from `numpy`)

"""

from __future__ import division, print_function, absolute_import

__all__ = ['who', 'source', 'info', 'doccer', 'pade',
           'comb', 'factorial', 'factorial2', 'factorialk', 'logsumexp']

from . import doccer
from .common import *
from numpy import who as _who, source as _source, info as _info
import numpy as np
from scipy.interpolate._pade import pade as _pade
from scipy.special import (comb as _comb, logsumexp as _lsm,
        factorial as _fact, factorial2 as _fact2, factorialk as _factk)

import sys

_msg = ("Importing `%(name)s` from scipy.misc is deprecated in scipy 1.0.0. Use "
        "`scipy.special.%(name)s` instead.")
comb = np.deprecate(_comb, message=_msg % {"name": _comb.__name__})
logsumexp = np.deprecate(_lsm, message=_msg % {"name": _lsm.__name__})
factorial = np.deprecate(_fact, message=_msg % {"name": _fact.__name__})
factorial2 = np.deprecate(_fact2, message=_msg % {"name": _fact2.__name__})
factorialk = np.deprecate(_factk, message=_msg % {"name": _factk.__name__})

_msg = ("Importing `pade` from scipy.misc is deprecated in scipy 1.0.0. Use "
        "`scipy.interpolate.pade` instead.")
pade = np.deprecate(_pade, message=_msg)

_msg = ("Importing `%(name)s` from scipy.misc is deprecated in scipy 1.0.0. Use "
        "`numpy.%(name)s` instead.")
who = np.deprecate(_who, message=_msg % {"name": "who"})
source = np.deprecate(_source, message=_msg % {"name": "source"})

@np.deprecate(message=_msg % {"name": "info.(..., toplevel='scipy')"})
def info(object=None,maxwidth=76,output=sys.stdout,toplevel='scipy'):
    return _info(object, maxwidth, output, toplevel)


info.__doc__ = _info.__doc__
del sys

try:
    from .pilutil import *
    from . import pilutil
    __all__ += pilutil.__all__
    del pilutil
except ImportError:
    pass

from . import common
__all__ += common.__all__
del common

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
