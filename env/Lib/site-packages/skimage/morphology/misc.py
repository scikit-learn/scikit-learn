"""Miscellaneous morphology functions."""

import numpy as np
import functools
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

from .._shared.utils import warn
from ._misc_cy import _remove_objects_by_distance


# Our function names don't exactly correspond to ndimages.
# This dictionary translates from our names to scipy's.
funcs = ('erosion', 'dilation', 'opening', 'closing')
skimage2ndimage = {x: 'grey_' + x for x in funcs}

# These function names are the same in ndimage.
funcs = (
    'binary_erosion',
    'binary_dilation',
    'binary_opening',
    'binary_closing',
    'black_tophat',
    'white_tophat',
)
skimage2ndimage.update({x: x for x in funcs})


def default_footprint(func):
    """Decorator to add a default footprint to morphology functions.

    Parameters
    ----------
    func : function
        A morphology function such as erosion, dilation, opening, closing,
        white_tophat, or black_tophat.

    Returns
    -------
    func_out : function
        The function, using a default footprint of same dimension
        as the input image with connectivity 1.

    """

    @functools.wraps(func)
    def func_out(image, footprint=None, *args, **kwargs):
        if footprint is None:
            footprint = ndi.generate_binary_structure(image.ndim, 1)
        return func(image, footprint=footprint, *args, **kwargs)

    return func_out


def _check_dtype_supported(ar):
    # Should use `issubdtype` for bool below, but there's a bug in numpy 1.7
    if not (ar.dtype == bool or np.issubdtype(ar.dtype, np.integer)):
        raise TypeError(
            "Only bool or integer image types are supported. " f"Got {ar.dtype}."
        )


def remove_small_objects(ar, min_size=64, connectivity=1, *, out=None):
    """Remove objects smaller than the specified size.

    Expects ar to be an array with labeled objects, and removes objects
    smaller than min_size. If `ar` is bool, the image is first labeled.
    This leads to potentially different behavior for bool and 0-and-1
    arrays.

    Parameters
    ----------
    ar : ndarray (arbitrary shape, int or bool type)
        The array containing the objects of interest. If the array type is
        int, the ints must be non-negative.
    min_size : int, optional (default: 64)
        The smallest allowable object size.
    connectivity : int, {1, 2, ..., ar.ndim}, optional (default: 1)
        The connectivity defining the neighborhood of a pixel. Used during
        labelling if `ar` is bool.
    out : ndarray
        Array of the same shape as `ar`, into which the output is
        placed. By default, a new array is created.

    Raises
    ------
    TypeError
        If the input array is of an invalid type, such as float or string.
    ValueError
        If the input array contains negative values.

    Returns
    -------
    out : ndarray, same shape and type as input `ar`
        The input array with small connected components removed.

    See Also
    --------
    skimage.morphology.remove_objects_by_distance

    Examples
    --------
    >>> from skimage import morphology
    >>> a = np.array([[0, 0, 0, 1, 0],
    ...               [1, 1, 1, 0, 0],
    ...               [1, 1, 1, 0, 1]], bool)
    >>> b = morphology.remove_small_objects(a, 6)
    >>> b
    array([[False, False, False, False, False],
           [ True,  True,  True, False, False],
           [ True,  True,  True, False, False]])
    >>> c = morphology.remove_small_objects(a, 7, connectivity=2)
    >>> c
    array([[False, False, False,  True, False],
           [ True,  True,  True, False, False],
           [ True,  True,  True, False, False]])
    >>> d = morphology.remove_small_objects(a, 6, out=a)
    >>> d is a
    True

    """
    # Raising type error if not int or bool
    _check_dtype_supported(ar)

    if out is None:
        out = ar.copy()
    else:
        out[:] = ar

    if min_size == 0:  # shortcut for efficiency
        return out

    if out.dtype == bool:
        footprint = ndi.generate_binary_structure(ar.ndim, connectivity)
        ccs = np.zeros_like(ar, dtype=np.int32)
        ndi.label(ar, footprint, output=ccs)
    else:
        ccs = out

    try:
        component_sizes = np.bincount(ccs.ravel())
    except ValueError:
        raise ValueError(
            "Negative value labels are not supported. Try "
            "relabeling the input with `scipy.ndimage.label` or "
            "`skimage.morphology.label`."
        )

    if len(component_sizes) == 2 and out.dtype != bool:
        warn(
            "Only one label was provided to `remove_small_objects`. "
            "Did you mean to use a boolean array?"
        )

    too_small = component_sizes < min_size
    too_small_mask = too_small[ccs]
    out[too_small_mask] = 0

    return out


def remove_small_holes(ar, area_threshold=64, connectivity=1, *, out=None):
    """Remove contiguous holes smaller than the specified size.

    Parameters
    ----------
    ar : ndarray (arbitrary shape, int or bool type)
        The array containing the connected components of interest.
    area_threshold : int, optional (default: 64)
        The maximum area, in pixels, of a contiguous hole that will be filled.
        Replaces `min_size`.
    connectivity : int, {1, 2, ..., ar.ndim}, optional (default: 1)
        The connectivity defining the neighborhood of a pixel.
    out : ndarray
        Array of the same shape as `ar` and bool dtype, into which the
        output is placed. By default, a new array is created.

    Raises
    ------
    TypeError
        If the input array is of an invalid type, such as float or string.
    ValueError
        If the input array contains negative values.

    Returns
    -------
    out : ndarray, same shape and type as input `ar`
        The input array with small holes within connected components removed.

    Examples
    --------
    >>> from skimage import morphology
    >>> a = np.array([[1, 1, 1, 1, 1, 0],
    ...               [1, 1, 1, 0, 1, 0],
    ...               [1, 0, 0, 1, 1, 0],
    ...               [1, 1, 1, 1, 1, 0]], bool)
    >>> b = morphology.remove_small_holes(a, 2)
    >>> b
    array([[ True,  True,  True,  True,  True, False],
           [ True,  True,  True,  True,  True, False],
           [ True, False, False,  True,  True, False],
           [ True,  True,  True,  True,  True, False]])
    >>> c = morphology.remove_small_holes(a, 2, connectivity=2)
    >>> c
    array([[ True,  True,  True,  True,  True, False],
           [ True,  True,  True, False,  True, False],
           [ True, False, False,  True,  True, False],
           [ True,  True,  True,  True,  True, False]])
    >>> d = morphology.remove_small_holes(a, 2, out=a)
    >>> d is a
    True

    Notes
    -----
    If the array type is int, it is assumed that it contains already-labeled
    objects. The labels are not kept in the output image (this function always
    outputs a bool image). It is suggested that labeling is completed after
    using this function.

    """
    _check_dtype_supported(ar)

    # Creates warning if image is an integer image
    if ar.dtype != bool:
        warn(
            "Any labeled images will be returned as a boolean array. "
            "Did you mean to use a boolean array?",
            UserWarning,
        )

    if out is not None:
        if out.dtype != bool:
            raise TypeError("out dtype must be bool")
    else:
        out = ar.astype(bool, copy=True)

    # Creating the inverse of ar
    np.logical_not(ar, out=out)

    # removing small objects from the inverse of ar
    out = remove_small_objects(out, area_threshold, connectivity, out=out)

    np.logical_not(out, out=out)

    return out


def remove_objects_by_distance(
    label_image,
    min_distance,
    *,
    priority=None,
    p_norm=2,
    spacing=None,
    out=None,
):
    """Remove objects, in specified order, until remaining are a minimum distance apart.

    Remove labeled objects from an image until the remaining ones are spaced
    more than a given distance from one another. By default, smaller objects
    are removed first.

    Parameters
    ----------
    label_image : ndarray of integers
        An n-dimensional array containing object labels, e.g. as returned by
        :func:`~.label`. A value of zero is considered background, all other
        object IDs must be positive integers.
    min_distance : int or float
        Remove objects whose distance to other objects is not greater than this
        positive value. Objects with a lower `priority` are removed first.
    priority : ndarray, optional
        Defines the priority with which objects are removed. Expects a
        1-dimensional array of length
        :func:`np.amax(label_image) + 1 <numpy.amax>` that contains the priority
        for each object's label at the respective index. Objects with a lower value
        are removed first until all remaining objects fulfill the distance
        requirement. If not given, priority is given to objects with a higher
        number of samples and their label value second.
    p_norm : int or float, optional
        The Minkowski distance of order p, used to calculate the distance
        between objects. The default ``2`` corresponds to the Euclidean
        distance, ``1`` to the "Manhattan" distance, and ``np.inf`` to the
        Chebyshev distance.
    spacing : sequence of float, optional
        The pixel spacing along each axis of `label_image`. If not specified,
        a grid spacing of unity (1) is implied.
    out : ndarray, optional
        Array of the same shape and dtype as `image`, into which the output is
        placed. By default, a new array is created.

    Returns
    -------
    out : ndarray
        Array of the same shape as `label_image`, for which objects that violate
        the `min_distance` condition were removed.

    See Also
    --------
    skimage.morphology.remove_small_objects
        Remove objects smaller than the specified size.

    Notes
    -----
    The basic steps of this algorithm work as follows:

    1. Find the indices for of all given objects and separate them depending on
       if they point to an object's border or not.
    2. Sort indices by their label value, ensuring that indices which point to
       the same object are next to each other. This optimization allows finding
       all parts of an object, simply by stepping to the neighboring indices.
    3. Sort boundary indices by `priority`. Use a stable-sort to preserve the
       ordering from the previous sorting step. If `priority` is not given,
       use :func:`numpy.bincount` as a fallback.
    4. Construct a :class:`scipy.spatial.cKDTree` from the boundary indices.
    5. Iterate across boundary indices in priority-sorted order, and query the
       kd-tree for objects that are too close. Remove ones that are and don't
       take them into account when evaluating other objects later on.

    The performance of this algorithm depends on the number of samples in
    `label_image` that belong to an object's border.

    Examples
    --------
    >>> import skimage as ski
    >>> ski.morphology.remove_objects_by_distance(np.array([2, 0, 1, 1]), 2)
    array([0, 0, 1, 1])
    >>> ski.morphology.remove_objects_by_distance(
    ...     np.array([2, 0, 1, 1]), 2, priority=np.array([0, 1, 9])
    ... )
    array([2, 0, 0, 0])
    >>> label_image = np.array(
    ...     [[8, 0, 0, 0, 0, 0, 0, 0, 0, 9, 9],
    ...      [8, 8, 8, 0, 0, 0, 0, 0, 0, 9, 9],
    ...      [0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0],
    ...      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    ...      [0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 0],
    ...      [2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    ...      [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    ...      [0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7]]
    ... )
    >>> ski.morphology.remove_objects_by_distance(
    ...     label_image, min_distance=3
    ... )
    array([[8, 0, 0, 0, 0, 0, 0, 0, 0, 9, 9],
           [8, 8, 8, 0, 0, 0, 0, 0, 0, 9, 9],
           [0, 0, 0, 0, 0, 0, 0, 0, 9, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7]])
    """
    if min_distance < 0:
        raise ValueError(f"min_distance must be >= 0, was {min_distance}")
    if not np.issubdtype(label_image.dtype, np.integer):
        raise ValueError(
            f"`label_image` must be of integer dtype, got {label_image.dtype}"
        )
    if out is None:
        out = label_image.copy(order="C")
    elif out is not label_image:
        out[:] = label_image
    # May create a copy if order is not C, account for that later
    out_raveled = out.ravel(order="C")

    if spacing is not None:
        spacing = np.array(spacing)
        if spacing.shape != (out.ndim,) or spacing.min() <= 0:
            raise ValueError(
                "`spacing` must contain exactly one positive factor "
                "for each dimension of `label_image`"
            )

    indices = np.flatnonzero(out_raveled)
    # Optimization: Split indices into those on the object boundaries and inner
    # ones. The KDTree is built only from the boundary indices, which reduces
    # the size of the critical loop significantly! Remaining indices are only
    # used to remove the inner parts of objects as well.
    if (spacing is None or np.all(spacing[0] == spacing)) and p_norm <= 2:
        # For unity spacing we can make the borders more sparse by using a
        # lower connectivity
        footprint = ndi.generate_binary_structure(out.ndim, 1)
    else:
        footprint = ndi.generate_binary_structure(out.ndim, out.ndim)
    border = (
        ndi.maximum_filter(out, footprint=footprint)
        != ndi.minimum_filter(out, footprint=footprint)
    ).ravel()[indices]
    border_indices = indices[border]
    inner_indices = indices[~border]

    if border_indices.size == 0:
        # Image without any or only one object, return early
        return out

    # Sort by label ID first, so that IDs of the same object are contiguous
    # in the sorted index. This allows fast discovery of the whole object by
    # simple iteration up or down the index!
    border_indices = border_indices[np.argsort(out_raveled[border_indices])]
    inner_indices = inner_indices[np.argsort(out_raveled[inner_indices])]

    if priority is None:
        if not np.can_cast(out.dtype, np.intp, casting="safe"):
            # bincount expects intp (32-bit) on WASM or i386, so down-cast to that
            priority = np.bincount(out_raveled.astype(np.intp, copy=False))
        else:
            priority = np.bincount(out_raveled)
    # `priority` can only be indexed by positive object IDs,
    # `border_indices` contains all unique sorted IDs so check the lowest / first
    smallest_id = out_raveled[border_indices[0]]
    if smallest_id < 0:
        raise ValueError(f"found object with negative ID {smallest_id!r}")

    try:
        # Sort by priority second using a stable sort to preserve the contiguous
        # sorting of objects. Because each pixel in an object has the same
        # priority we don't need to worry about separating objects.
        border_indices = border_indices[
            np.argsort(priority[out_raveled[border_indices]], kind="stable")[::-1]
        ]
    except IndexError as error:
        # Use np.amax only for the exception path to provide a nicer error message
        expected_shape = (np.amax(out_raveled) + 1,)
        if priority.shape != expected_shape:
            raise ValueError(
                "shape of `priority` must be (np.amax(label_image) + 1,), "
                f"expected {expected_shape}, got {priority.shape} instead"
            ) from error
        else:
            raise

    # Construct kd-tree from unraveled border indices (optionally scale by `spacing`)
    unraveled_indices = np.unravel_index(border_indices, out.shape)
    if spacing is not None:
        unraveled_indices = tuple(
            unraveled_indices[dim] * spacing[dim] for dim in range(out.ndim)
        )
    kdtree = cKDTree(data=np.asarray(unraveled_indices, dtype=np.float64).T)

    _remove_objects_by_distance(
        out=out_raveled,
        border_indices=border_indices,
        inner_indices=inner_indices,
        kdtree=kdtree,
        min_distance=min_distance,
        p_norm=p_norm,
        shape=label_image.shape,
    )

    if out_raveled.base is not out:
        # `out_raveled` is a copy, re-assign
        out[:] = out_raveled.reshape(out.shape)
    return out
