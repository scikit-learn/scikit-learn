import collections as coll
import numpy as np
from scipy import ndimage as ndi

from ..util import img_as_float
from ..color import guess_spatial_dimensions
from .._shared.utils import warn, convert_to_float


__all__ = ['gaussian']


def _convert_input(image, preserve_range):
    if preserve_range:
        image = image.astype(np.double)
    else:
        image = img_as_float(image)
    return image


def gaussian(image, sigma=1, output=None, mode='nearest', cval=0,
             multichannel=None, preserve_range=False, truncate=4.0):
    """Multi-dimensional Gaussian filter.

    Parameters
    ----------
    image : array-like
        Input image (grayscale or color) to filter.
    sigma : scalar or sequence of scalars, optional
        Standard deviation for Gaussian kernel. The standard
        deviations of the Gaussian filter are given for each axis as a
        sequence, or as a single number, in which case it is equal for
        all axes.
    output : array, optional
        The ``output`` parameter passes an array in which to store the
        filter output.
    mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional
        The `mode` parameter determines how the array borders are
        handled, where `cval` is the value when mode is equal to
        'constant'. Default is 'nearest'.
    cval : scalar, optional
        Value to fill past edges of input if `mode` is 'constant'. Default
        is 0.0
    multichannel : bool, optional (default: None)
        Whether the last axis of the image is to be interpreted as multiple
        channels. If True, each channel is filtered separately (channels are
        not mixed together). Only 3 channels are supported. If `None`,
        the function will attempt to guess this, and raise a warning if
        ambiguous, when the array has shape (M, N, 3).
    preserve_range : bool, optional
        Whether to keep the original range of values. Otherwise, the input
        image is converted according to the conventions of `img_as_float`.
    truncate : float, optional
        Truncate the filter at this many standard deviations.

    Returns
    -------
    filtered_image : ndarray
        the filtered array

    Notes
    -----
    This function is a wrapper around :func:`scipy.ndi.gaussian_filter`.

    Integer arrays are converted to float.

    The multi-dimensional filter is implemented as a sequence of
    one-dimensional convolution filters. The intermediate arrays are
    stored in the same data type as the output. Therefore, for output
    types with a limited precision, the results may be imprecise
    because intermediate results may be stored with insufficient
    precision.

    Examples
    --------

    >>> a = np.zeros((3, 3))
    >>> a[1, 1] = 1
    >>> a
    array([[ 0.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  0.]])
    >>> gaussian(a, sigma=0.4)  # mild smoothing
    array([[ 0.00163116,  0.03712502,  0.00163116],
           [ 0.03712502,  0.84496158,  0.03712502],
           [ 0.00163116,  0.03712502,  0.00163116]])
    >>> gaussian(a, sigma=1)  # more smoothing
    array([[ 0.05855018,  0.09653293,  0.05855018],
           [ 0.09653293,  0.15915589,  0.09653293],
           [ 0.05855018,  0.09653293,  0.05855018]])
    >>> # Several modes are possible for handling boundaries
    >>> gaussian(a, sigma=1, mode='reflect')
    array([[ 0.08767308,  0.12075024,  0.08767308],
           [ 0.12075024,  0.16630671,  0.12075024],
           [ 0.08767308,  0.12075024,  0.08767308]])
    >>> # For RGB images, each is filtered separately
    >>> from skimage.data import astronaut
    >>> image = astronaut()
    >>> filtered_img = gaussian(image, sigma=1, multichannel=True)

    """

    spatial_dims = guess_spatial_dimensions(image)
    if spatial_dims is None and multichannel is None:
        msg = ("Images with dimensions (M, N, 3) are interpreted as 2D+RGB "
               "by default. Use `multichannel=False` to interpret as "
               "3D image with last dimension of length 3.")
        warn(RuntimeWarning(msg))
        multichannel = True
    if np.any(np.asarray(sigma) < 0.0):
        raise ValueError("Sigma values less than zero are not valid")
    if multichannel:
        # do not filter across channels
        if not isinstance(sigma, coll.Iterable):
            sigma = [sigma] * (image.ndim - 1)
        if len(sigma) != image.ndim:
            sigma = np.concatenate((np.asarray(sigma), [0]))
    image = convert_to_float(image, preserve_range)
    return ndi.gaussian_filter(image, sigma, mode=mode, cval=cval,
                               truncate=truncate)
