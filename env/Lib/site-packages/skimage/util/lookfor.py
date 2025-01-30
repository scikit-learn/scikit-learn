import sys

from .._vendored.numpy_lookfor import lookfor as _lookfor


def lookfor(what):
    """Do a keyword search on scikit-image docstrings and print results.

    .. warning::

        This function may also print results that are not part of
        scikit-image's public API.

    Parameters
    ----------
    what : str
        Words to look for.

    Examples
    --------
    >>> import skimage as ski
    >>> ski.util.lookfor('regular_grid')
    Search results for 'regular_grid'
    ---------------------------------
    skimage.util.regular_grid
        Find `n_points` regularly spaced along `ar_shape`.
    skimage.util.lookfor
        Do a keyword search on scikit-image docstrings and print results.
    """
    return _lookfor(what, sys.modules[__name__.split('.')[0]])
