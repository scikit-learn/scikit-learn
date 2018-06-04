from __future__ import division, print_function, absolute_import

import numpy as np


_have_pil = True
try:
    from scipy.misc.pilutil import imread as _imread
except ImportError:
    _have_pil = False


__all__ = ['imread']


# Use the implementation of `imread` in `scipy.misc.pilutil.imread`.
# If it weren't for the different names of the first arguments of
# ndimage.io.imread and misc.pilutil.imread, we could simplify this file
# by writing
#     from scipy.misc.pilutil import imread
# Unfortunately, because the argument names are different, that
# introduces a backwards incompatibility.

@np.deprecate(message="`imread` is deprecated in SciPy 1.0.0.\n"
                      "Use ``matplotlib.pyplot.imread`` instead.")
def imread(fname, flatten=False, mode=None):
    if _have_pil:
        return _imread(fname, flatten, mode)
    raise ImportError("Could not import the Python Imaging Library (PIL)"
                      " required to load image files.  Please refer to"
                      " http://pillow.readthedocs.org/en/latest/installation.html"
                      " for installation instructions.")


if _have_pil and _imread.__doc__ is not None:
    imread.__doc__ = _imread.__doc__.replace('name : str', 'fname : str')
