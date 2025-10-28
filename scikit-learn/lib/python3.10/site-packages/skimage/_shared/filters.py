"""Filters used across multiple skimage submodules.

These are defined here to avoid circular imports.

The unit tests remain under skimage/filters/tests/
"""

from collections.abc import Iterable

import numpy as np
from scipy import ndimage as ndi

from .._shared.utils import (
    _supported_float_type,
    convert_to_float,
)


def gaussian(
    image,
    sigma=1.0,
    *,
    mode='nearest',
    cval=0,
    preserve_range=False,
    truncate=4.0,
    channel_axis=None,
    out=None,
):
    """Multi-dimensional Gaussian filter.

    Parameters
    ----------
    image : ndarray
        Input image (grayscale or color) to filter.
    sigma : scalar or sequence of scalars, optional
        Standard deviation for Gaussian kernel. The standard
        deviations of the Gaussian filter are given for each axis as a
        sequence, or as a single number, in which case it is equal for
        all axes.
    mode : {'reflect', 'constant', 'nearest', 'mirror', 'wrap'}, optional
        The ``mode`` parameter determines how the array borders are
        handled, where ``cval`` is the value when mode is equal to
        'constant'. Default is 'nearest'.
    cval : scalar, optional
        Value to fill past edges of input if ``mode`` is 'constant'. Default
        is 0.0
    preserve_range : bool, optional
        If True, keep the original range of values. Otherwise, the input
        ``image`` is converted according to the conventions of ``img_as_float``
        (Normalized first to values [-1.0 ; 1.0] or [0 ; 1.0] depending on
        dtype of input)

        For more information, see:
        https://scikit-image.org/docs/dev/user_guide/data_types.html
    truncate : float, optional
        Truncate the filter at this many standard deviations.
    channel_axis : int or None, optional
        If None, the image is assumed to be a grayscale (single channel) image.
        Otherwise, this parameter indicates which axis of the array corresponds
        to channels.

        .. versionadded:: 0.19
           `channel_axis` was added in 0.19.
    out : ndarray, optional
        If given, the filtered image will be stored in this array.

        .. versionadded:: 0.23
            `out` was added in 0.23.

    Returns
    -------
    filtered_image : ndarray
        the filtered array

    Notes
    -----
    This function is a wrapper around :func:`scipy.ndimage.gaussian_filter`.

    Integer arrays are converted to float.

    `out` should be of floating-point data type since `gaussian` converts the
    input `image` to float. If `out` is not provided, another array
    will be allocated and returned as the result.

    The multi-dimensional filter is implemented as a sequence of
    one-dimensional convolution filters. The intermediate arrays are
    stored in the same data type as the output. Therefore, for output
    types with a limited precision, the results may be imprecise
    because intermediate results may be stored with insufficient
    precision.

    Examples
    --------
    >>> import skimage as ski
    >>> a = np.zeros((3, 3))
    >>> a[1, 1] = 1
    >>> a
    array([[0., 0., 0.],
           [0., 1., 0.],
           [0., 0., 0.]])
    >>> ski.filters.gaussian(a, sigma=0.4)  # mild smoothing
    array([[0.00163116, 0.03712502, 0.00163116],
           [0.03712502, 0.84496158, 0.03712502],
           [0.00163116, 0.03712502, 0.00163116]])
    >>> ski.filters.gaussian(a, sigma=1)  # more smoothing
    array([[0.05855018, 0.09653293, 0.05855018],
           [0.09653293, 0.15915589, 0.09653293],
           [0.05855018, 0.09653293, 0.05855018]])
    >>> # Several modes are possible for handling boundaries
    >>> ski.filters.gaussian(a, sigma=1, mode='reflect')
    array([[0.08767308, 0.12075024, 0.08767308],
           [0.12075024, 0.16630671, 0.12075024],
           [0.08767308, 0.12075024, 0.08767308]])
    >>> # For RGB images, each is filtered separately
    >>> image = ski.data.astronaut()
    >>> filtered_img = ski.filters.gaussian(image, sigma=1, channel_axis=-1)

    """
    if np.any(np.asarray(sigma) < 0.0):
        raise ValueError("Sigma values less than zero are not valid")
    if channel_axis is not None:
        # do not filter across channels
        if not isinstance(sigma, Iterable):
            sigma = [sigma] * (image.ndim - 1)
        if len(sigma) == image.ndim - 1:
            sigma = list(sigma)
            sigma.insert(channel_axis % image.ndim, 0)
    image = convert_to_float(image, preserve_range)
    float_dtype = _supported_float_type(image.dtype)
    image = image.astype(float_dtype, copy=False)
    if (out is not None) and (not np.issubdtype(out.dtype, np.floating)):
        raise ValueError(f"dtype of `out` must be float; got {out.dtype!r}.")
    return ndi.gaussian_filter(
        image, sigma, output=out, mode=mode, cval=cval, truncate=truncate
    )
