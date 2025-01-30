"""
Grayscale morphological operations
"""

import warnings

import numpy as np
from scipy import ndimage as ndi

from .footprints import _footprint_is_sequence, mirror_footprint, pad_footprint
from .misc import default_footprint
from .._shared.utils import DEPRECATED


__all__ = ['erosion', 'dilation', 'opening', 'closing', 'white_tophat', 'black_tophat']


def _iterate_gray_func(gray_func, image, footprints, out, mode, cval):
    """Helper to call `gray_func` for each footprint in a sequence.

    `gray_func` is a morphology function that accepts `footprint`, `output`,
    `mode` and `cval` keyword arguments (e.g. `scipy.ndimage.grey_erosion`).
    """
    fp, num_iter = footprints[0]
    gray_func(image, footprint=fp, output=out, mode=mode, cval=cval)
    for _ in range(1, num_iter):
        gray_func(out.copy(), footprint=fp, output=out, mode=mode, cval=cval)
    for fp, num_iter in footprints[1:]:
        # Note: out.copy() because the computation cannot be in-place!
        for _ in range(num_iter):
            gray_func(out.copy(), footprint=fp, output=out, mode=mode, cval=cval)
    return out


def _shift_footprint(footprint, shift_x, shift_y):
    """Shift the binary image `footprint` in the left and/or up.

    This only affects 2D footprints with even number of rows
    or columns.

    Parameters
    ----------
    footprint : 2D array, shape (M, N)
        The input footprint.
    shift_x, shift_y : bool or None
        Whether to move `footprint` along each axis. If ``None``, the
        array is not modified along that dimension.

    Returns
    -------
    out : 2D array, shape (M + int(shift_x), N + int(shift_y))
        The shifted footprint.
    """
    footprint = np.asarray(footprint)
    if footprint.ndim != 2:
        # do nothing for 1D or 3D or higher footprints
        return footprint
    m, n = footprint.shape
    if m % 2 == 0:
        extra_row = np.zeros((1, n), footprint.dtype)
        if shift_x:
            footprint = np.vstack((footprint, extra_row))
        else:
            footprint = np.vstack((extra_row, footprint))
        m += 1
    if n % 2 == 0:
        extra_col = np.zeros((m, 1), footprint.dtype)
        if shift_y:
            footprint = np.hstack((footprint, extra_col))
        else:
            footprint = np.hstack((extra_col, footprint))
    return footprint


def _shift_footprints(footprint, shift_x, shift_y):
    """Shifts the footprints, whether it's a single array or a sequence.

    See `_shift_footprint`, which is called for each array in the sequence.
    """
    if shift_x is DEPRECATED and shift_y is DEPRECATED:
        return footprint

    warning_msg = (
        "The parameters `shift_x` and `shift_y` are deprecated since v0.23 and "
        "will be removed in v0.26. Use `pad_footprint` or modify the footprint"
        "manually instead."
    )
    warnings.warn(warning_msg, FutureWarning, stacklevel=4)

    if _footprint_is_sequence(footprint):
        return tuple((_shift_footprint(fp, shift_x, shift_y), n) for fp, n in footprint)
    return _shift_footprint(footprint, shift_x, shift_y)


def _min_max_to_constant_mode(dtype, mode, cval):
    """Replace 'max' and 'min' with appropriate 'cval' and 'constant' mode."""
    if mode == "max":
        mode = "constant"
        if np.issubdtype(dtype, bool):
            cval = True
        elif np.issubdtype(dtype, np.integer):
            cval = np.iinfo(dtype).max
        else:
            cval = np.inf
    elif mode == "min":
        mode = "constant"
        if np.issubdtype(dtype, bool):
            cval = False
        elif np.issubdtype(dtype, np.integer):
            cval = np.iinfo(dtype).min
        else:
            cval = -np.inf
    return mode, cval


_SUPPORTED_MODES = {
    "reflect",
    "constant",
    "nearest",
    "mirror",
    "wrap",
    "max",
    "min",
    "ignore",
}


@default_footprint
def erosion(
    image,
    footprint=None,
    out=None,
    shift_x=DEPRECATED,
    shift_y=DEPRECATED,
    *,
    mode="reflect",
    cval=0.0,
):
    """Return grayscale morphological erosion of an image.

    Morphological erosion sets a pixel at (i,j) to the minimum over all pixels
    in the neighborhood centered at (i,j). Erosion shrinks bright regions and
    enlarges dark regions.

    Parameters
    ----------
    image : ndarray
        Image array.
    footprint : ndarray or tuple, optional
        The neighborhood expressed as a 2-D array of 1's and 0's.
        If None, use a cross-shaped footprint (connectivity=1). The footprint
        can also be provided as a sequence of smaller footprints as described
        in the notes below.
    out : ndarrays, optional
        The array to store the result of the morphology. If None is
        passed, a new array will be allocated.
    mode : str, optional
        The `mode` parameter determines how the array borders are handled.
        Valid modes are: 'reflect', 'constant', 'nearest', 'mirror', 'wrap',
        'max', 'min', or 'ignore'.
        If 'max' or 'ignore', pixels outside the image domain are assumed
        to be the maximum for the image's dtype, which causes them to not
        influence the result. Default is 'reflect'.
    cval : scalar, optional
        Value to fill past edges of input if `mode` is 'constant'. Default
        is 0.0.

        .. versionadded:: 0.23
            `mode` and `cval` were added in 0.23.

    Returns
    -------
    eroded : array, same shape as `image`
        The result of the morphological erosion.

    Other Parameters
    ----------------
    shift_x, shift_y : DEPRECATED

        .. deprecated:: 0.23

    Notes
    -----
    For ``uint8`` (and ``uint16`` up to a certain bit-depth) data, the
    lower algorithm complexity makes the :func:`skimage.filters.rank.minimum`
    function more efficient for larger images and footprints.

    The footprint can also be a provided as a sequence of 2-tuples where the
    first element of each 2-tuple is a footprint ndarray and the second element
    is an integer describing the number of times it should be iterated. For
    example ``footprint=[(np.ones((9, 1)), 1), (np.ones((1, 9)), 1)]``
    would apply a 9x1 footprint followed by a 1x9 footprint resulting in a net
    effect that is the same as ``footprint=np.ones((9, 9))``, but with lower
    computational cost. Most of the builtin footprints such as
    :func:`skimage.morphology.disk` provide an option to automatically generate
    a footprint sequence of this type.

    For even-sized footprints, :func:`skimage.morphology.binary_erosion` and
    this function produce an output that differs: one is shifted by one pixel
    compared to the other.

    Examples
    --------
    >>> # Erosion shrinks bright regions
    >>> import numpy as np
    >>> from skimage.morphology import footprint_rectangle
    >>> bright_square = np.array([[0, 0, 0, 0, 0],
    ...                           [0, 1, 1, 1, 0],
    ...                           [0, 1, 1, 1, 0],
    ...                           [0, 1, 1, 1, 0],
    ...                           [0, 0, 0, 0, 0]], dtype=np.uint8)
    >>> erosion(bright_square, footprint_rectangle((3, 3)))
    array([[0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0],
           [0, 0, 1, 0, 0],
           [0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0]], dtype=uint8)

    """
    if out is None:
        out = np.empty_like(image)

    if mode not in _SUPPORTED_MODES:
        raise ValueError(f"unsupported mode, got {mode!r}")
    if mode == "ignore":
        mode = "max"
    mode, cval = _min_max_to_constant_mode(image.dtype, mode, cval)

    footprint = _shift_footprints(footprint, shift_x, shift_y)
    footprint = pad_footprint(footprint, pad_end=False)
    if not _footprint_is_sequence(footprint):
        footprint = [(footprint, 1)]

    out = _iterate_gray_func(
        gray_func=ndi.grey_erosion,
        image=image,
        footprints=footprint,
        out=out,
        mode=mode,
        cval=cval,
    )
    return out


@default_footprint
def dilation(
    image,
    footprint=None,
    out=None,
    shift_x=DEPRECATED,
    shift_y=DEPRECATED,
    *,
    mode="reflect",
    cval=0.0,
):
    """Return grayscale morphological dilation of an image.

    Morphological dilation sets the value of a pixel to the maximum over all
    pixel values within a local neighborhood centered about it. The values
    where the footprint is 1 define this neighborhood.
    Dilation enlarges bright regions and shrinks dark regions.

    Parameters
    ----------
    image : ndarray
        Image array.
    footprint : ndarray or tuple, optional
        The neighborhood expressed as a 2-D array of 1's and 0's.
        If None, use a cross-shaped footprint (connectivity=1). The footprint
        can also be provided as a sequence of smaller footprints as described
        in the notes below.
    out : ndarray, optional
        The array to store the result of the morphology. If None is
        passed, a new array will be allocated.
    mode : str, optional
        The `mode` parameter determines how the array borders are handled.
        Valid modes are: 'reflect', 'constant', 'nearest', 'mirror', 'wrap',
        'max', 'min', or 'ignore'.
        If 'min' or 'ignore', pixels outside the image domain are assumed
        to be the maximum for the image's dtype, which causes them to not
        influence the result. Default is 'reflect'.
    cval : scalar, optional
        Value to fill past edges of input if `mode` is 'constant'. Default
        is 0.0.

        .. versionadded:: 0.23
            `mode` and `cval` were added in 0.23.

    Returns
    -------
    dilated : uint8 array, same shape and type as `image`
        The result of the morphological dilation.

    Other Parameters
    ----------------
    shift_x, shift_y : DEPRECATED

        .. deprecated:: 0.23

    Notes
    -----
    For ``uint8`` (and ``uint16`` up to a certain bit-depth) data, the lower
    algorithm complexity makes the :func:`skimage.filters.rank.maximum`
    function more efficient for larger images and footprints.

    The footprint can also be a provided as a sequence of 2-tuples where the
    first element of each 2-tuple is a footprint ndarray and the second element
    is an integer describing the number of times it should be iterated. For
    example ``footprint=[(np.ones((9, 1)), 1), (np.ones((1, 9)), 1)]``
    would apply a 9x1 footprint followed by a 1x9 footprint resulting in a net
    effect that is the same as ``footprint=np.ones((9, 9))``, but with lower
    computational cost. Most of the builtin footprints such as
    :func:`skimage.morphology.disk` provide an option to automatically generate
    a footprint sequence of this type.

    For non-symmetric footprints, :func:`skimage.morphology.binary_dilation`
    and :func:`skimage.morphology.dilation` produce an output that differs:
    `binary_dilation` mirrors the footprint, whereas `dilation` does not.

    Examples
    --------
    >>> # Dilation enlarges bright regions
    >>> import numpy as np
    >>> from skimage.morphology import footprint_rectangle
    >>> bright_pixel = np.array([[0, 0, 0, 0, 0],
    ...                          [0, 0, 0, 0, 0],
    ...                          [0, 0, 1, 0, 0],
    ...                          [0, 0, 0, 0, 0],
    ...                          [0, 0, 0, 0, 0]], dtype=np.uint8)
    >>> dilation(bright_pixel, footprint_rectangle((3, 3)))
    array([[0, 0, 0, 0, 0],
           [0, 1, 1, 1, 0],
           [0, 1, 1, 1, 0],
           [0, 1, 1, 1, 0],
           [0, 0, 0, 0, 0]], dtype=uint8)

    """
    if out is None:
        out = np.empty_like(image)

    if mode not in _SUPPORTED_MODES:
        raise ValueError(f"unsupported mode, got {mode!r}")
    if mode == "ignore":
        mode = "min"
    mode, cval = _min_max_to_constant_mode(image.dtype, mode, cval)

    footprint = _shift_footprints(footprint, shift_x, shift_y)
    footprint = pad_footprint(footprint, pad_end=False)
    # Note that `ndi.grey_dilation` mirrors the footprint and this
    # additional inversion should be removed in skimage2, see gh-6676.
    footprint = mirror_footprint(footprint)
    if not _footprint_is_sequence(footprint):
        footprint = [(footprint, 1)]

    out = _iterate_gray_func(
        gray_func=ndi.grey_dilation,
        image=image,
        footprints=footprint,
        out=out,
        mode=mode,
        cval=cval,
    )
    return out


@default_footprint
def opening(image, footprint=None, out=None, *, mode="reflect", cval=0.0):
    """Return grayscale morphological opening of an image.

    The morphological opening of an image is defined as an erosion followed by
    a dilation. Opening can remove small bright spots (i.e. "salt") and connect
    small dark cracks. This tends to "open" up (dark) gaps between (bright)
    features.

    Parameters
    ----------
    image : ndarray
        Image array.
    footprint : ndarray or tuple, optional
        The neighborhood expressed as a 2-D array of 1's and 0's.
        If None, use a cross-shaped footprint (connectivity=1). The footprint
        can also be provided as a sequence of smaller footprints as described
        in the notes below.
    out : ndarray, optional
        The array to store the result of the morphology. If None
        is passed, a new array will be allocated.
    mode : str, optional
        The `mode` parameter determines how the array borders are handled.
        Valid modes are: 'reflect', 'constant', 'nearest', 'mirror', 'wrap',
        'max', 'min', or 'ignore'.
        If 'ignore', pixels outside the image domain are assumed
        to be the maximum for the image's dtype in the erosion, and minimum
        in the dilation, which causes them to not influence the result.
        Default is 'reflect'.
    cval : scalar, optional
        Value to fill past edges of input if `mode` is 'constant'. Default
        is 0.0.

        .. versionadded:: 0.23
            `mode` and `cval` were added in 0.23.

    Returns
    -------
    opening : array, same shape and type as `image`
        The result of the morphological opening.

    Notes
    -----
    The footprint can also be a provided as a sequence of 2-tuples where the
    first element of each 2-tuple is a footprint ndarray and the second element
    is an integer describing the number of times it should be iterated. For
    example ``footprint=[(np.ones((9, 1)), 1), (np.ones((1, 9)), 1)]``
    would apply a 9x1 footprint followed by a 1x9 footprint resulting in a net
    effect that is the same as ``footprint=np.ones((9, 9))``, but with lower
    computational cost. Most of the builtin footprints such as
    :func:`skimage.morphology.disk` provide an option to automatically generate
    a footprint sequence of this type.

    Examples
    --------
    >>> # Open up gap between two bright regions (but also shrink regions)
    >>> import numpy as np
    >>> from skimage.morphology import footprint_rectangle
    >>> bad_connection = np.array([[1, 0, 0, 0, 1],
    ...                            [1, 1, 0, 1, 1],
    ...                            [1, 1, 1, 1, 1],
    ...                            [1, 1, 0, 1, 1],
    ...                            [1, 0, 0, 0, 1]], dtype=np.uint8)
    >>> opening(bad_connection, footprint_rectangle((3, 3)))
    array([[0, 0, 0, 0, 0],
           [1, 1, 0, 1, 1],
           [1, 1, 0, 1, 1],
           [1, 1, 0, 1, 1],
           [0, 0, 0, 0, 0]], dtype=uint8)

    """
    footprint = pad_footprint(footprint, pad_end=False)
    eroded = erosion(image, footprint, mode=mode, cval=cval)
    out = dilation(eroded, mirror_footprint(footprint), out=out, mode=mode, cval=cval)
    return out


@default_footprint
def closing(image, footprint=None, out=None, *, mode="reflect", cval=0.0):
    """Return grayscale morphological closing of an image.

    The morphological closing of an image is defined as a dilation followed by
    an erosion. Closing can remove small dark spots (i.e. "pepper") and connect
    small bright cracks. This tends to "close" up (dark) gaps between (bright)
    features.

    Parameters
    ----------
    image : ndarray
        Image array.
    footprint : ndarray or tuple, optional
        The neighborhood expressed as a 2-D array of 1's and 0's.
        If None, use a cross-shaped footprint (connectivity=1). The footprint
        can also be provided as a sequence of smaller footprints as described
        in the notes below.
    out : ndarray, optional
        The array to store the result of the morphology. If None,
        a new array will be allocated.
    mode : str, optional
        The `mode` parameter determines how the array borders are handled.
        Valid modes are: 'reflect', 'constant', 'nearest', 'mirror', 'wrap',
        'max', 'min', or 'ignore'.
        If 'ignore', pixels outside the image domain are assumed
        to be the maximum for the image's dtype in the erosion, and minimum
        in the dilation, which causes them to not influence the result.
        Default is 'reflect'.
    cval : scalar, optional
        Value to fill past edges of input if `mode` is 'constant'. Default
        is 0.0.

        .. versionadded:: 0.23
            `mode` and `cval` were added in 0.23.

    Returns
    -------
    closing : array, same shape and type as `image`
        The result of the morphological closing.

    Notes
    -----
    The footprint can also be a provided as a sequence of 2-tuples where the
    first element of each 2-tuple is a footprint ndarray and the second element
    is an integer describing the number of times it should be iterated. For
    example ``footprint=[(np.ones((9, 1)), 1), (np.ones((1, 9)), 1)]``
    would apply a 9x1 footprint followed by a 1x9 footprint resulting in a net
    effect that is the same as ``footprint=np.ones((9, 9))``, but with lower
    computational cost. Most of the builtin footprints such as
    :func:`skimage.morphology.disk` provide an option to automatically generate
    a footprint sequence of this type.

    Examples
    --------
    >>> # Close a gap between two bright lines
    >>> import numpy as np
    >>> from skimage.morphology import footprint_rectangle
    >>> broken_line = np.array([[0, 0, 0, 0, 0],
    ...                         [0, 0, 0, 0, 0],
    ...                         [1, 1, 0, 1, 1],
    ...                         [0, 0, 0, 0, 0],
    ...                         [0, 0, 0, 0, 0]], dtype=np.uint8)
    >>> closing(broken_line, footprint_rectangle((3, 3)))
    array([[0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0],
           [1, 1, 1, 1, 1],
           [0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0]], dtype=uint8)

    """
    footprint = pad_footprint(footprint, pad_end=False)
    dilated = dilation(image, footprint, mode=mode, cval=cval)
    out = erosion(dilated, mirror_footprint(footprint), out=out, mode=mode, cval=cval)
    return out


@default_footprint
def white_tophat(image, footprint=None, out=None, *, mode="reflect", cval=0.0):
    """Return white top hat of an image.

    The white top hat of an image is defined as the image minus its
    morphological opening. This operation returns the bright spots of the image
    that are smaller than the footprint.

    Parameters
    ----------
    image : ndarray
        Image array.
    footprint : ndarray or tuple, optional
        The neighborhood expressed as a 2-D array of 1's and 0's.
        If None, use a cross-shaped footprint (connectivity=1). The footprint
        can also be provided as a sequence of smaller footprints as described
        in the notes below.
    out : ndarray, optional
        The array to store the result of the morphology. If None
        is passed, a new array will be allocated.
    mode : str, optional
        The `mode` parameter determines how the array borders are handled.
        Valid modes are: 'reflect', 'constant', 'nearest', 'mirror', 'wrap',
        'max', 'min', or 'ignore'. See :func:`skimage.morphology.opening`.
        Default is 'reflect'.
    cval : scalar, optional
        Value to fill past edges of input if `mode` is 'constant'. Default
        is 0.0.

        .. versionadded:: 0.23
            `mode` and `cval` were added in 0.23.

    Returns
    -------
    out : array, same shape and type as `image`
        The result of the morphological white top hat.

    Notes
    -----
    The footprint can also be a provided as a sequence of 2-tuples where the
    first element of each 2-tuple is a footprint ndarray and the second element
    is an integer describing the number of times it should be iterated. For
    example ``footprint=[(np.ones((9, 1)), 1), (np.ones((1, 9)), 1)]``
    would apply a 9x1 footprint followed by a 1x9 footprint resulting in a net
    effect that is the same as ``footprint=np.ones((9, 9))``, but with lower
    computational cost. Most of the builtin footprints such as
    :func:`skimage.morphology.disk` provide an option to automatically generate
    a footprint sequence of this type.

    See Also
    --------
    black_tophat

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Top-hat_transform

    Examples
    --------
    >>> # Subtract gray background from bright peak
    >>> import numpy as np
    >>> from skimage.morphology import footprint_rectangle
    >>> bright_on_gray = np.array([[2, 3, 3, 3, 2],
    ...                            [3, 4, 5, 4, 3],
    ...                            [3, 5, 9, 5, 3],
    ...                            [3, 4, 5, 4, 3],
    ...                            [2, 3, 3, 3, 2]], dtype=np.uint8)
    >>> white_tophat(bright_on_gray, footprint_rectangle((3, 3)))
    array([[0, 0, 0, 0, 0],
           [0, 0, 1, 0, 0],
           [0, 1, 5, 1, 0],
           [0, 0, 1, 0, 0],
           [0, 0, 0, 0, 0]], dtype=uint8)

    """
    if out is image:
        # We need a temporary image
        opened = opening(image, footprint, mode=mode, cval=cval)
        if np.issubdtype(opened.dtype, bool):
            np.logical_xor(out, opened, out=out)
        else:
            out -= opened
        return out

    # Else write intermediate result into output image
    out = opening(image, footprint, out=out, mode=mode, cval=cval)
    if np.issubdtype(out.dtype, bool):
        np.logical_xor(image, out, out=out)
    else:
        np.subtract(image, out, out=out)
    return out


@default_footprint
def black_tophat(image, footprint=None, out=None, *, mode="reflect", cval=0.0):
    """Return black top hat of an image.

    The black top hat of an image is defined as its morphological closing minus
    the original image. This operation returns the dark spots of the image that
    are smaller than the footprint. Note that dark spots in the
    original image are bright spots after the black top hat.

    Parameters
    ----------
    image : ndarray
        Image array.
    footprint : ndarray or tuple, optional
        The neighborhood expressed as a 2-D array of 1's and 0's.
        If None, use a cross-shaped footprint (connectivity=1). The footprint
        can also be provided as a sequence of smaller footprints as described
        in the notes below.
    out : ndarray, optional
        The array to store the result of the morphology. If None
        is passed, a new array will be allocated.
    mode : str, optional
        The `mode` parameter determines how the array borders are handled.
        Valid modes are: 'reflect', 'constant', 'nearest', 'mirror', 'wrap',
        'max', 'min', or 'ignore'. See :func:`skimage.morphology.closing`.
        Default is 'reflect'.
    cval : scalar, optional
        Value to fill past edges of input if `mode` is 'constant'. Default
        is 0.0.

        .. versionadded:: 0.23
            `mode` and `cval` were added in 0.23.

    Returns
    -------
    out : array, same shape and type as `image`
        The result of the morphological black top hat.

    Notes
    -----
    The footprint can also be a provided as a sequence of 2-tuples where the
    first element of each 2-tuple is a footprint ndarray and the second element
    is an integer describing the number of times it should be iterated. For
    example ``footprint=[(np.ones((9, 1)), 1), (np.ones((1, 9)), 1)]``
    would apply a 9x1 footprint followed by a 1x9 footprint resulting in a net
    effect that is the same as ``footprint=np.ones((9, 9))``, but with lower
    computational cost. Most of the builtin footprints such as
    :func:`skimage.morphology.disk` provide an option to automatically generate
    a footprint sequence of this type.

    See Also
    --------
    white_tophat

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Top-hat_transform

    Examples
    --------
    >>> # Change dark peak to bright peak and subtract background
    >>> import numpy as np
    >>> from skimage.morphology import footprint_rectangle
    >>> dark_on_gray = np.array([[7, 6, 6, 6, 7],
    ...                          [6, 5, 4, 5, 6],
    ...                          [6, 4, 0, 4, 6],
    ...                          [6, 5, 4, 5, 6],
    ...                          [7, 6, 6, 6, 7]], dtype=np.uint8)
    >>> black_tophat(dark_on_gray, footprint_rectangle((3, 3)))
    array([[0, 0, 0, 0, 0],
           [0, 0, 1, 0, 0],
           [0, 1, 5, 1, 0],
           [0, 0, 1, 0, 0],
           [0, 0, 0, 0, 0]], dtype=uint8)

    """
    if out is image:
        # We need a temporary image
        closed = closing(image, footprint, mode=mode, cval=cval)
        if np.issubdtype(closed.dtype, bool):
            np.logical_xor(closed, out, out=out)
        else:
            np.subtract(closed, out, out=out)
        return out

    out = closing(image, footprint, out=out, mode=mode, cval=cval)
    if np.issubdtype(out.dtype, np.bool_):
        np.logical_xor(out, image, out=out)
    else:
        out -= image
    return out
