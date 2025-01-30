"""
Binary morphological operations
"""

import numpy as np
from scipy import ndimage as ndi

from .footprints import _footprint_is_sequence, pad_footprint
from .misc import default_footprint


def _iterate_binary_func(binary_func, image, footprint, out, border_value):
    """Helper to call `binary_func` for each footprint in a sequence.

    binary_func is a binary morphology function that accepts "structure",
    "output" and "iterations" keyword arguments
    (e.g. `scipy.ndimage.binary_erosion`).
    """
    fp, num_iter = footprint[0]
    binary_func(
        image, structure=fp, output=out, iterations=num_iter, border_value=border_value
    )
    for fp, num_iter in footprint[1:]:
        # Note: out.copy() because the computation cannot be in-place!
        #       SciPy <= 1.7 did not automatically make a copy if needed.
        binary_func(
            out.copy(),
            structure=fp,
            output=out,
            iterations=num_iter,
            border_value=border_value,
        )
    return out


# The default_footprint decorator provides a diamond footprint as
# default with the same dimension as the input image and size 3 along each
# axis.
@default_footprint
def binary_erosion(image, footprint=None, out=None, *, mode='ignore'):
    """Return fast binary morphological erosion of an image.

    This function returns the same result as grayscale erosion but performs
    faster for binary images.

    Morphological erosion sets a pixel at ``(i,j)`` to the minimum over all
    pixels in the neighborhood centered at ``(i,j)``. Erosion shrinks bright
    regions and enlarges dark regions.

    Parameters
    ----------
    image : ndarray
        Binary input image.
    footprint : ndarray or tuple, optional
        The neighborhood expressed as a 2-D array of 1's and 0's.
        If None, use a cross-shaped footprint (connectivity=1). The footprint
        can also be provided as a sequence of smaller footprints as described
        in the notes below.
    out : ndarray of bool, optional
        The array to store the result of the morphology. If None is
        passed, a new array will be allocated.
    mode : str, optional
        The `mode` parameter determines how the array borders are handled.
        Valid modes are: 'max', 'min', 'ignore'.
        If 'max' or 'ignore', pixels outside the image domain are assumed
        to be `True`, which causes them to not influence the result.
        Default is 'ignore'.

        .. versionadded:: 0.23
            `mode` was added in 0.23.

    Returns
    -------
    eroded : ndarray of bool or uint
        The result of the morphological erosion taking values in
        ``[False, True]``.

    Notes
    -----
    The footprint can also be a provided as a sequence of 2-tuples where the
    first element of each 2-tuple is a footprint ndarray and the second element
    is an integer describing the number of times it should be iterated. For
    example ``footprint=[(np.ones((9, 1)), 1), (np.ones((1, 9)), 1)]``
    would apply a 9x1 footprint followed by a 1x9 footprint resulting in a net
    effect that is the same as ``footprint=np.ones((9, 9))``, but with lower
    computational cost. Most of the builtin footprints such as
    :func:`skimage.morphology.disk` provide an option to automatically generate a
    footprint sequence of this type.

    For even-sized footprints, :func:`skimage.morphology.erosion` and
    this function produce an output that differs: one is shifted by one pixel
    compared to the other.

    See also
    --------
    skimage.morphology.isotropic_erosion

    """
    if out is None:
        out = np.empty(image.shape, dtype=bool)

    if mode not in {"max", "min", "ignore"}:
        raise ValueError(f"unsupported mode, got {mode!r}")
    border_value = False if mode == 'min' else True

    footprint = pad_footprint(footprint, pad_end=True)
    if not _footprint_is_sequence(footprint):
        footprint = [(footprint, 1)]

    out = _iterate_binary_func(
        binary_func=ndi.binary_erosion,
        image=image,
        footprint=footprint,
        out=out,
        border_value=border_value,
    )
    return out


@default_footprint
def binary_dilation(image, footprint=None, out=None, *, mode='ignore'):
    """Return fast binary morphological dilation of an image.

    This function returns the same result as grayscale dilation but performs
    faster for binary images.

    Morphological dilation sets a pixel at ``(i,j)`` to the maximum over all
    pixels in the neighborhood centered at ``(i,j)``. Dilation enlarges bright
    regions and shrinks dark regions.

    Parameters
    ----------
    image : ndarray
        Binary input image.
    footprint : ndarray or tuple, optional
        The neighborhood expressed as a 2-D array of 1's and 0's.
        If None, use a cross-shaped footprint (connectivity=1). The footprint
        can also be provided as a sequence of smaller footprints as described
        in the notes below.
    out : ndarray of bool, optional
        The array to store the result of the morphology. If None is
        passed, a new array will be allocated.
    mode : str, optional
        The `mode` parameter determines how the array borders are handled.
        Valid modes are: 'max', 'min', 'ignore'.
        If 'min' or 'ignore', pixels outside the image domain are assumed
        to be `False`, which causes them to not influence the result.
        Default is 'ignore'.

        .. versionadded:: 0.23
            `mode` was added in 0.23.

    Returns
    -------
    dilated : ndarray of bool or uint
        The result of the morphological dilation with values in
        ``[False, True]``.

    Notes
    -----
    The footprint can also be a provided as a sequence of 2-tuples where the
    first element of each 2-tuple is a footprint ndarray and the second element
    is an integer describing the number of times it should be iterated. For
    example ``footprint=[(np.ones((9, 1)), 1), (np.ones((1, 9)), 1)]``
    would apply a 9x1 footprint followed by a 1x9 footprint resulting in a net
    effect that is the same as ``footprint=np.ones((9, 9))``, but with lower
    computational cost. Most of the builtin footprints such as
    :func:`skimage.morphology.disk` provide an option to automatically generate a
    footprint sequence of this type.

    For non-symmetric footprints, :func:`skimage.morphology.binary_dilation`
    and :func:`skimage.morphology.dilation` produce an output that differs:
    `binary_dilation` mirrors the footprint, whereas `dilation` does not.

    See also
    --------
    skimage.morphology.isotropic_dilation

    """
    if out is None:
        out = np.empty(image.shape, dtype=bool)

    if mode not in {"max", "min", "ignore"}:
        raise ValueError(f"unsupported mode, got {mode!r}")
    border_value = True if mode == 'max' else False

    footprint = pad_footprint(footprint, pad_end=True)
    if not _footprint_is_sequence(footprint):
        footprint = [(footprint, 1)]

    out = _iterate_binary_func(
        binary_func=ndi.binary_dilation,
        image=image,
        footprint=footprint,
        out=out,
        border_value=border_value,
    )
    return out


@default_footprint
def binary_opening(image, footprint=None, out=None, *, mode='ignore'):
    """Return fast binary morphological opening of an image.

    This function returns the same result as grayscale opening but performs
    faster for binary images.

    The morphological opening on an image is defined as an erosion followed by
    a dilation. Opening can remove small bright spots (i.e. "salt") and connect
    small dark cracks. This tends to "open" up (dark) gaps between (bright)
    features.

    Parameters
    ----------
    image : ndarray
        Binary input image.
    footprint : ndarray or tuple, optional
        The neighborhood expressed as a 2-D array of 1's and 0's.
        If None, use a cross-shaped footprint (connectivity=1). The footprint
        can also be provided as a sequence of smaller footprints as described
        in the notes below.
    out : ndarray of bool, optional
        The array to store the result of the morphology. If None
        is passed, a new array will be allocated.
    mode : str, optional
        The `mode` parameter determines how the array borders are handled.
        Valid modes are: 'max', 'min', 'ignore'.
        If 'ignore', pixels outside the image domain are assumed to be `True`
        for the erosion and `False` for the dilation, which causes them to not
        influence the result. Default is 'ignore'.

        .. versionadded:: 0.23
            `mode` was added in 0.23.

    Returns
    -------
    opening : ndarray of bool
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
    :func:`skimage.morphology.disk` provide an option to automatically generate a
    footprint sequence of this type.

    See also
    --------
    skimage.morphology.isotropic_opening

    """
    tmp = binary_erosion(image, footprint, mode=mode)
    out = binary_dilation(tmp, footprint, out=out, mode=mode)
    return out


@default_footprint
def binary_closing(image, footprint=None, out=None, *, mode='ignore'):
    """Return fast binary morphological closing of an image.

    This function returns the same result as grayscale closing but performs
    faster for binary images.

    The morphological closing on an image is defined as a dilation followed by
    an erosion. Closing can remove small dark spots (i.e. "pepper") and connect
    small bright cracks. This tends to "close" up (dark) gaps between (bright)
    features.

    Parameters
    ----------
    image : ndarray
        Binary input image.
    footprint : ndarray or tuple, optional
        The neighborhood expressed as a 2-D array of 1's and 0's.
        If None, use a cross-shaped footprint (connectivity=1). The footprint
        can also be provided as a sequence of smaller footprints as described
        in the notes below.
    out : ndarray of bool, optional
        The array to store the result of the morphology. If None,
        is passed, a new array will be allocated.
    mode : str, optional
        The `mode` parameter determines how the array borders are handled.
        Valid modes are: 'max', 'min', 'ignore'.
        If 'ignore', pixels outside the image domain are assumed to be `True`
        for the erosion and `False` for the dilation, which causes them to not
        influence the result. Default is 'ignore'.

        .. versionadded:: 0.23
            `mode` was added in 0.23.

    Returns
    -------
    closing : ndarray of bool
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
    :func:`skimage.morphology.disk` provide an option to automatically generate a
    footprint sequence of this type.

    See also
    --------
    skimage.morphology.isotropic_closing

    """
    tmp = binary_dilation(image, footprint, mode=mode)
    out = binary_erosion(tmp, footprint, out=out, mode=mode)
    return out
